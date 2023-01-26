## Run phylogenetic scan for mutations enriched in different hosts
## 
## This jupyter notebook contains commands to load the necessary dependencies, and perform the gene-wide scan for mutations enriched in different host species. The way I've structured this code is by coding the actual computational components of the method into other scripts/documents, and then simply calling them here. My intention is that this should avoid making edits to the underlying code, and instead altering parameters like the number of simulated iterations, the input tree, required counts, etc... The real heft of the method is coded in 2 other notebooks: `calculate-enrichment-scores` and `simulate-mutation-gain-loss-markov-chain`. This notebook will read those notebooks in, and run them. 
## 
## The basic principle of this method is to identify mutations that are enriched in different host species. The premise is that if mutations are adaptive/beneficial in a given host species, then that mutation should be found more commonly in one host species than the other. In this project, we are trying to find mutations that are responsible for human adaptation of H5N1 viruses, so our host comparison is between humans and birds. 
## 
## Comparing mutations across host species can be accomplished by simply looking at an alignment, but using a tree-based approach has many benefits. For one, in an alignment, there is no information about how sequences/samples are related to each other. Using a tree-based approach allows you to control for associations based on related samples and to specifically look for independent patterns of recurrent mutations. It allows you to run controls for population structure/sampling, which are important. 
## 
## In this approach, you first enumerate every single mutation that occurs on the phylogenetic tree, and count how many times that mutation is found in the hosts you want to query. Here, for each mutation on the tree, we record how many times it occurs on the tree and how many human and bird tips it is associated with. Because we are not reconstructing ancestral states onto internal nodes on the phylogeny, mutations that arise on nodes are not assigned a host state. Instead, we trace the path from that node to all its descendants, and record tips that have retained that mutation. For each mutation, we then construct a 2x2 contingency table: 
## 
## |host|presence|absence|
## |:------|:-------|:------|
## |host 1|A|B| 
## |host 2|C|D|
## 
## where A, B, C, and D are counts of the mutation's presence and absence in host 1 and host 2. Based on the methodology [here](https://github.com/sheppardlab/pGWAS/blob/master/assomap_given_phylo.py), published in [this paper](https://www.nature.com/articles/s41467-018-07368-7##Sec10) and detailed in lines 245-273 in the code, you can then calculate an odds ratio, describing the degree of association of a particular mutation with a particular host as: `OR = (A * D)/(B * C)`. You can take this one step further and use this 2x2 table to calculate a Fisher's exact test p-value. When you run a fisher's exact test in python, it will output 2 values: an odds ratio and a p-value. The odds ratio output is exactly equivalent to the odds ratio calculated here. For testing purposes, this code currently calcules an odds ratio and p-value by fishers exact test for each mutation. 
## 
## At the end of the first part of the procedure, you then have a list of every mutation on the tree, it's counts in human and avian tips, an odds ratio, and a p-value. 
## 
## A huge, recurring issue for phylogenetic-based methods, and for GWAS methods, is sampling and population structure. There are lots of instances, for example, where a single mutation will be present in a large cluster of human tips, resulting in a large odds ratio. However, we can't be sure that this mutation was not in the avian population from which these humans were infected, because sampling is incomplete. It could absolutely be true that this mutation resulted in a virus more able to infect people, but it could also be in this human cluster simply because in this outbreak, only humans were sampled. 
## 
## One way to control for this issue is to evaluate how skewed counts among these hosts are assuming no selection at all. For example, if we placed mutations onto the tree randomly, with no underlying selection, how frequently would we see ones that are significantly enriched in one host or the other by the odds ratio? By definition, these associations would be due to chance alone. We can then use this "null distribution" to find examples in our true tree that are more skewed than we see in this random, null. To do this, we take the tree and strip it of its mutations. We then simulate a single mutation that toggles on and off on the tree and count how many times it is found in human and avian tips as above. For each simulated tree, we calculate the same odds ratio and p-value. We do this many times to generate a distribution of odds ratios and p-values in these simulated datasets. We can then use these distributions as null distributions, and take the most extreme 5% as our cutoff for statistical significance. 


## Import modules we will need
import pandas as pd 
import time
import multiprocessing as mp
from functools import partial


## Import additional python modules
## All functions, excluding get_iteration_list, are defined in the first three modules. Importing the modules will allow all
## of the functions coded in those notebooks to be available in this one. The first notebook does part 1, and the second
## does part 2. Functions for loading and manipulating trees are contained in the third. Additionally, all config data
## which must be specificed by the user is defined in the config file.
import calculate_enrichment_scores_across_tree_JSON as calenr
import simulate_mutation_gain_loss_markov_chain as simmut
import tree_manager as tm
import config as cfg
import write_files


## One thing I really like doing is appending the current date to all my code and output files. This helps me 
## keep track of changes when I rerun analyses on different days, and also helps prevent me from accidentally 
## overwriting files by accident. This uses the python datetime module to store the current date in YYYY-MM-DD 
## format as a variable (string) that can be appended to output files
from datetime import date
current_date = str(date.today())



def get_iteration_list(iterations, cores):
    """Function to create a list with total iterations evenly split among number of cores for multiprocessing"""
    div, rem = divmod(iterations, cores)
    return [div+1]*rem + [div]*(cores-rem)



if __name__ == "__main__":
    ## Part 1: Infer mutations on tree, and calculate enrichment scores and p-values
    ## 
    ## In this first part of the notebook, we will be reading in a tree, enumerating every mutation on the tree and returning each mutation with an enrichment score (odds ratio) and p-value as assessed by a Fisher's exact test. There are a few required inputs here, which the user should specify which stem from me writing this to be flexible. 
    ## 
    ## 1. **minimum required count:** For every amino acid on the tree, the enrichment score will not be returned if it is present in less than `minimum required count` tips. I played around with this a little bit, and ended up setting it to 0. The reason is that regardless of counts, I wanted to calculate the score for every mutation, and then post-filter afterwards. 
    ## 
    ## 2. **method:** there are 2 possible methods: `proportions` and `counts`. In `counts`, we calculate an odds ratio as: `(A*D)/(B*C)`, where `A`,`B`,`C`, and `D` are counts of tips. In `proportions`, we calculate `(A+D)-(B+C)` where the cells are the proportion of total tips in each category. The `counts` method output is the exact odds ratio calculated with a Fisher's exact test, and is more appropriate. 
    ## 
    ## There are a few outputs: the `times_detected_dict` outputs the number of times that the mutation arose on the tree. The counts in `scores_dict` represent the number of tips with each mutation. 
    ## 
    ## To run this on amino acids, put in the gene name under `gene`. To run on nucleotide mutations, replace gene with `nuc`. 


    ## Load tree, no_muts_tree, and pickled_tree
    if cfg.reload_trees == True:
        print("Reloading trees...")
        json_tree, tree, no_muts_tree, pickled_tree = tm.init_pickled_trees(cfg.output_folder_path)
        print("Finished loading trees")
    else:
        json_tree, tree, no_muts_tree, pickled_tree = tm.init_trees()
    

    ## Gather all mutations. The outputs, aa_muts and nt_muts, are lists of amino acid mutations
    ## and nucleotide mutations on the tree. Every mutation on the tree is included.
    aa_muts, nt_muts = calenr.gather_all_mut_on_tree(tree, cfg.gene)

    ## Determine all host tips. The output, total_host_tips_on_tree, is a dictionary with counts
    ## of the number of tips corresponding to each host on the tree
    total_host_tips_on_tree = calenr.return_all_host_tips(tree, cfg.host1, cfg.host2, cfg.host_annotation)

    ## Calculate the total branch length of the tree, in terms of mutations
    total_tree_branch_length, tree_branch_lengths = calenr.return_total_tree_branch_length(tree)

    ## Calculate enrichment scores for all mutations along the tree. must set method to be counts or proportions; 
    ## the host_counts variable in calculate_enrichmenet_scores is total_host_tips_on_tree
    scores_dict, times_detected_dict, branch_lengths_dict, host_counts_dict2 = calenr.calculate_enrichment_scores(tree, aa_muts, nt_muts, cfg.host1, cfg.host2, cfg.host_annotation, cfg.minimum_required_count, total_host_tips_on_tree, cfg.gene)



    ## Convert tree data to dataframes
    ## Create dataframes from each dictionary
    df1 = pd.DataFrame.from_dict(scores_dict, orient="index")
    df2 = pd.DataFrame.from_dict(times_detected_dict, orient="index", columns=["total_times_detected_on_tree"])
    df3 = pd.DataFrame.from_dict(branch_lengths_dict, orient="index", columns=["branch_length_with_mutation"])
    df4 = pd.DataFrame.from_dict(host_counts_dict2, orient="index", columns=[cfg.host1, cfg.host2, "other"])

    ## Merge dataframes together; pandas join is a merge on the index
    df5 = df1.join(df2.join(df3.join(df4)))





    ## Part 2: simulate mutation gain and loss across the tree to generate a null 
    ## 
    ## The output for `sims_times_detected` will be the number of times in each iteration that the simulated mutation arose. This includes occurrences on internal nodes and on terminal nodes. 
    ## 
    ## Scores of 0 occur when the mutation was never detected/present in host 1. 
    ## 
    ## For this analysis, you will need to specify a few things: 
    ## 1. **iterations:** iterations specifies how many times to simulate mutation gain or loss across the tree.


    ## Get number of cores and create a list to split iterations evenly among them
    cores = mp.cpu_count()
    iter_list = get_iteration_list(cfg.iterations, cores)

    ## Create partial function for simmut.perform_simulations with all arguments excluding iterations
    sim_part = partial(simmut.perform_simulations, pickled_tree, cfg.gene, total_tree_branch_length, cfg.host1, cfg.host2, cfg.host_annotation, cfg.minimum_required_count, total_host_tips_on_tree)

    ## Start timer
    start_time = time.time()

    ## Start multiprocessing pool, run simulations, then close the pool once all cores have finished
    mp.set_start_method('fork')
    pool = mp.Pool()
    sim_data = pool.map(sim_part, iter_list) # Run simmut.perform_simulations using arguments specified in sim_part, with iterations split among cores as specified in iter_list
    pool.close()
    pool.join()

    ## End timer
    total_time_seconds = time.time() - start_time
    total_time_minutes = total_time_seconds/60
    total_time_hours = total_time_minutes/60
    print("This took", total_time_seconds, "seconds (", total_time_minutes," minutes,", total_time_hours," hours) to generate", cfg.iterations, "simulated trees")



    ## Convert simulation data to dataframes
    ## sim_data[core][field][iteration]
        ## core = [0, ..., mp.cpu_count() - 1]
        ## field = [0, 1, 2, 3]; 0 = sim_scores, 1 = sim_times_detected, 2 = branches_that_mutated, 3 = all_branches_dict
        ## iteration = iter_list[core]
    ## These dictionaries are nested, and each core uses same indexing, so need to manually create idx

    ## Create dataframe with sim_scores
    df6list = []
    idx = 0
    for core in range(len(sim_data)):
        for iteration in sim_data[core][0]:
            x = pd.DataFrame.from_dict(sim_data[core][0][iteration], orient="index")
            x['simulation_iteration'] = idx
            idx += 1
            x.reset_index(inplace=True)
            df6list.append(x)
    df6 = pd.concat(df6list)

    ## Create dataframe with sim_times_detected
    df7list = []
    idx = 0
    for core in range(len(sim_data)):
        for iteration in sim_data[core][1]:
            y = pd.DataFrame.from_dict(sim_data[core][1][iteration], orient="index", columns=["times_detected_on_tree"])
            y["simulation_iteration"] = idx
            idx += 1
            y.reset_index(inplace=True)
            df7list.append(y)
    df7 = pd.concat(df7list)

    ## Merge them together; pandas join is a merge on the index
    df8 = df6.merge(df7, on=["simulation_iteration","index"])

    ## If in testing mode, create additional dataframes for simulation validation
    if cfg.testing_mode == True:
        ## Create dataframe with branch lengths and the number of times mutated in all simulations (excludes non-mutated branches)
        sim_branch_length_dict = {}
        sim_branch_mut_times_dict = {}
        for core in range(len(sim_data)):
            for x in sim_data[core][2]:
                sim_branch_length_dict[x] = sim_data[core][2][x]['branch_length']
                sim_branch_mut_times_dict[x] = sim_branch_mut_times_dict.get(x, 0) + sim_data[core][2][x]['times_mutated']
        df9 = pd.DataFrame({'Length':pd.Series(sim_branch_length_dict),'Times':pd.Series(sim_branch_mut_times_dict)})

        ## Create dataframe with all branches, lengths and the number of times mutated in all simulations
        sim_branch_name_dict = {}
        sim_branch_length_dict = {}
        sim_branch_mut_times_dict = {}
        for core in range(len(sim_data)):
            for x in sim_data[core][3]:
                sim_branch_name_dict[x] = x
                sim_branch_length_dict[x] = sim_data[core][3][x]['branch_length']
                sim_branch_mut_times_dict[x] = sim_branch_mut_times_dict.get(x, 0) + sim_data[core][3][x]['times_mutated']
        df10 = pd.DataFrame({'Name':pd.Series(sim_branch_name_dict),'Length':pd.Series(sim_branch_length_dict),'Times':pd.Series(sim_branch_mut_times_dict)})



    ## Write output files
    folder_name = write_files.make_next_folder()
    write_files.write_config(folder_name)
    write_files.write_baltic_tree(folder_name, tree)
    write_files.write_json_tree(folder_name, json_tree)
    if cfg.testing_mode == True:
        write_files.write_dfs(folder_name, df5, df8, df9, df10)
    else:
        write_files.write_dfs(folder_name, df5, df8)