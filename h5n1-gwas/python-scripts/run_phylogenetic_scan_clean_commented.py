#!/usr/bin/env python
# coding: utf-8

# # Run phylogenetic scan for mutations enriched in different hosts
# 
# This jupyter notebook contains commands to load the necessary dependencies, and perform the gene-wide scan for mutations enriched in different host species. The way I've structured this code is by coding the actual computational components of the method into other scripts/documents, and then simply calling them here. My intention is that this should avoid making edits to the underlying code, and instead altering parameters like the number of simulated iterations, the input tree, required counts, etc... The real heft of the method is coded in 2 other notebooks: `calculate-enrichment-scores` and `simulate-mutation-gain-loss-markov-chain`. This notebook will read those notebooks in, and run them. 
# 
# The basic principle of this method is to identify mutations that are enriched in different host species. The premise is that if mutations are adaptive/beneficial in a given host species, then that mutation should be found more commonly in one host species than the other. In this project, we are trying to find mutations that are responsible for human adaptation of H5N1 viruses, so our host comparison is between humans and birds. 
# 
# Comparing mutations across host species can be accomplished by simply looking at an alignment, but using a tree-based approach has many benefits. For one, in an alignment, there is no information about how sequences/samples are related to each other. Using a tree-based approach allows you to control for associations based on related samples and to specifically look for independent patterns of recurrent mutations. It allows you to run controls for population structure/sampling, which are important. 
# 
# In this approach, you first enumerate every single mutation that occurs on the phylogenetic tree, and count how many times that mutation is found in the hosts you want to query. Here, for each mutation on the tree, we record how many times it occurs on the tree and how many human and bird tips it is associated with. Because we are not reconstructing ancestral states onto internal nodes on the phylogeny, mutations that arise on nodes are not assigned a host state. Instead, we trace the path from that node to all its descendants, and record tips that have retained that mutation. For each mutation, we then construct a 2x2 contingency table: 
# 
# |host|presence|absence|
# |:------|:-------|:------|
# |host 1|A|B| 
# |host 2|C|D|
# 
# where A, B, C, and D are counts of the mutation's presence and absence in host 1 and host 2. Based on the methodology [here](https://github.com/sheppardlab/pGWAS/blob/master/assomap_given_phylo.py), published in [this paper](https://www.nature.com/articles/s41467-018-07368-7#Sec10) and detailed in lines 245-273 in the code, you can then calculate an odds ratio, describing the degree of association of a particular mutation with a particular host as: `OR = (A * D)/(B * C)`. You can take this one step further and use this 2x2 table to calculate a Fisher's exact test p-value. When you run a fisher's exact test in python, it will output 2 values: an odds ratio and a p-value. The odds ratio output is exactly equivalent to the odds ratio calculated here. For testing purposes, this code currently calcules an odds ratio and p-value by fishers exact test for each mutation. 
# 
# At the end of the first part of the procedure, you then have a list of every mutation on the tree, it's counts in human and avian tips, an odds ratio, and a p-value. 
# 
# A huge, recurring issue for phylogenetic-based methods, and for GWAS methods, is sampling and population structure. There are lots of instances, for example, where a single mutation will be present in a large cluster of human tips, resulting in a large odds ratio. However, we can't be sure that this mutation was not in the avian population from which these humans were infected, because sampling is incomplete. It could absolutely be true that this mutation resulted in a virus more able to infect people, but it could also be in this human cluster simply because in this outbreak, only humans were sampled. 
# 
# One way to control for this issue is to evaluate how skewed counts among these hosts are assuming no selection at all. For example, if we placed mutations onto the tree randomly, with no underlying selection, how frequently would we see ones that are significantly enriched in one host or the other by the odds ratio? By definition, these associations would be due to chance alone. We can then use this "null distribution" to find examples in our true tree that are more skewed than we see in this random, null. To do this, we take the tree and strip it of its mutations. We then simulate a single mutation that toggles on and off on the tree and count how many times it is found in human and avian tips as above. For each simulated tree, we calculate the same odds ratio and p-value. We do this many times to generate a distribution of odds ratios and p-values in these simulated datasets. We can then use these distributions as null distributions, and take the most extreme 5% as our cutoff for statistical significance. 

"""import modules we will need. Pandas and numpy are for manipulating dataframes, time is for timestamps, glob 
is for manuevering through shell directries, and json is for processing json files"""

import glob
import json
import pandas as pd 
import numpy as np
import time
import multiprocessing as mp
from functools import partial

"""This is a module for running R within a juyter notebook. This is not a recommended way to do plotting, but
one that I have maintained because I'm pretty dependent on ggplot at this point"""
import rpy2



"""Import additional python modules.
All functions, excluding get_iteration_list, are defined in the first two modules. Importing the modules will allow all
of the functions coded in those notebooks to be available in this one. The first notebook does part 1, and the second
does part 2. Additionally, all config data which must be specificed by the user is defined in the config file."""
import calculate_enrichment_scores_across_tree_JSON as calenr
import simulate_mutation_gain_loss_markov_chain as simmut
import config as cfg



"""There are a couple different tools that can be used to parse trees, but my preferred is baltic. Baltic is a tool 
written in python by Gytis Dudas, and available here: https://github.com/evogytis/baltic. If installing via pip (recommended), 
import with `import baltic as bt`. Otherwise, import with imp from a local source file specified in the config file"""
    #pip
#import baltic as bt

    #local
import imp
bt = imp.load_source('baltic', cfg.baltic_path)



"""One thing I really like doing is appending the current date to all my code and output files. This helps me 
keep track of changes when I rerun analyses on different days, and also helps prevent me from accidentally 
overwriting files by accident. This uses the python datetime module to store the current date in YYYY-MM-DD 
format as a variable (string) that can be appended to output files"""
from datetime import date
current_date = str(date.today())



"""Function to create a list with total iterations evenly split among number of cores for multiprocessing"""
def get_iteration_list(iterations, cores):
    d,r = divmod(iterations, cores)
    return([d+1]*r + [d]*(cores-r))



if __name__ == "__main__":
    # # Part 1: Infer mutations on tree, and calculate enrichment scores and p-values
    # 
    # In this first part of the notebook, we will be reading in a tree, enumerating every mutation on the tree and returning each mutation with an enrichment score (odds ratio) and p-value as assessed by a Fisher's exact test. There are a few required inputs here, which the user should specify which stem from me writing this to be flexible. 
    # 
    # 1. **minimum required count:** For every amino acid on the tree, the enrichment score will not be returned if it is present in less than `minimum required count` tips. I played around with this a little bit, and ended up setting it to 0. The reason is that regardless of counts, I wanted to calculate the score for every mutation, and then post-filter afterwards. 
    # 
    # 2. **method:** there are 2 possible methods: `proportions` and `counts`. In `counts`, we calculate an odds ratio as: `(A*D)/(B*C)`, where `A`,`B`,`C`, and `D` are counts of tips. In `proportions`, we calculate `(A+D)-(B+C)` where the cells are the proportion of total tips in each category. The `counts` method output is the exact odds ratio calculated with a Fisher's exact test, and is more appropriate. 
    # 
    # There are a few outputs: the `times_detected_dict` outputs the number of times that the mutation arose on the tree. The counts in `scores_dict` represent the number of tips with each mutation. 
    # 
    # To run this on amino acids, put in the gene name under `gene`. To run on nucleotide mutations, replace gene with `nuc`. 


    """Load tree"""
    tree = calenr.read_in_tree_json(cfg.tree_path)

    """Gather all mutations. The outputs, aa_muts and nt_muts, are lists of amino acid mutations
    and nucleotide mutations on the tree. Every mutation on the tree is included."""
    aa_muts, nt_muts = calenr.gather_all_mut_on_tree(tree, cfg.gene)

    """Determine all host tips. The output, total_host_tips_on_tree, is a dictionary with counts
    of the number of tips corresponding to each host on the tree"""
    total_host_tips_on_tree = calenr.return_all_host_tips(tree, cfg.host1, cfg.host2, cfg.host_annotation)

    """Calculate the total branch length of the tree, in terms of mutations."""
    total_tree_branch_length, tree_branch_lengths = calenr.return_total_tree_branch_length(tree)

    """Calculate enrichment scores for all mutations along the tree. must set method to be counts or proportions; 
    the host_counts variable in calculate_enrichmenet_scores is total_host_tips_on_tree"""
    scores_dict, times_detected_dict, branch_lengths_dict, host_counts_dict2 = calenr.calculate_enrichment_scores(tree, aa_muts, nt_muts, cfg.host1, cfg.host2, cfg.host_annotation, cfg.minimum_required_count, cfg.method, total_host_tips_on_tree, cfg.gene)



    """Convert tree data to dataframes"""
    """Create dataframes from each dictionary"""
    df1 = pd.DataFrame.from_dict(scores_dict, orient="index")
    df2 = pd.DataFrame.from_dict(times_detected_dict, orient="index", columns=["total_times_detected_on_tree"])
    df3 = pd.DataFrame.from_dict(branch_lengths_dict, orient="index", columns=["branch_length_with_mutation"])
    df4 = pd.DataFrame.from_dict(host_counts_dict2, orient="index", columns=[cfg.host1, cfg.host2, "other"])

    """Merge dataframes together; pandas join is a merge on the index"""
    df5 = df1.join(df2.join(df3.join(df4)))

    """Write to csv"""
    output_filename = cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_data_" + current_date + ".tsv"
    df5.to_csv(output_filename, sep="\t", header=True, index_label="mutation")










    # # Part 2: simulate mutation gain and loss across the tree to generate a null 
    # 
    # The output for `sims_times_detected` will be the number of times in each iteration that the simulated mutation arose. This includes occurrences on internal nodes and on terminal nodes. 
    # 
    # Scores of 0 occur when the mutation was never detected/present in host 1. 
    # 
    # For this analysis, you will need to specify a few things: 
    # 1. **iterations:** iterations specifies how many times to simulate mutation gain or loss across the tree.


    """Get number of cores and create a list to split iterations evenly among them"""
    cores = mp.cpu_count()
    iter_list = get_iteration_list(cfg.iterations, cores)

    """Initialize the no_muts_tree"""
    simmut.init_tree(tree, cfg.gene)

    """Create partial function for simmut.perform_simulations with all arguments excluding iterations"""
    sim_part = partial(simmut.perform_simulations, cfg.gene, total_tree_branch_length, cfg.host1, cfg.host2, cfg.host_annotation, cfg.minimum_required_count, cfg.method, total_host_tips_on_tree)

    """Start timer"""
    start_time = time.time()

    """Start multiprocessing pool, run simulations, then close the pool once all cores have finished"""
    mp.set_start_method('fork')
    pool = mp.Pool()
    sim_data = pool.map(sim_part, iter_list) # Run simmut.perform_simulations using arguments specified in sim_part, with iterations split among cores as specified in iter_list
    pool.close()
    pool.join()

    """End timer"""
    total_time_seconds = time.time() - start_time
    total_time_minutes = total_time_seconds/60
    total_time_hours = total_time_minutes/60
    print("This took", total_time_seconds, "seconds (", total_time_minutes," minutes,", total_time_hours," hours) to generate", cfg.iterations, "simulated trees")



    """Convert simulation data to dataframes"""
    """sim_data[core][field][iteration]
        core = [0, ..., mp.cpu_count() - 1]
        field = [0, 1, 2]; 0 = sim_scores, 1 = sim_times_detected, 2 = branches_that_mutated
        iterations = iter_list[core]"""
    """These dictionaries are nested, and each core uses same indexing, so need to manually create idx"""

    """Create dataframe with sim_scores"""
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

    """Create dataframe with sim_times_detected"""
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

    ## I guess for the simulations, I don't record the distribution of host counts for each mutation
    # df8 = pd.DataFrame.from_dict(branch_lengths_dict, orient="index", columns=["branch_length_with_mutation"])
    # df9 = pd.DataFrame.from_dict(host_counts_dict2, orient="index", columns=[host1, host2, "other"])

    """Merge them together; pandas join is a merge on the index"""
    df8 = df6.merge(df7, on=["simulation_iteration","index"])

    """Write to csv"""
    output_filename = cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_simulated_" + current_date + ".tsv"
    df8.to_csv(output_filename, sep="\t", header=True, index=False)