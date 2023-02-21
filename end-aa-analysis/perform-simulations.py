from math import exp
import random
import json
import pandas as pd
from scipy.stats import fisher_exact
import multiprocessing as mp
import time



json_dir = '/Users/jort/coding/h5n1-mutations-rotation/build3-end-aa-analysis/build3-rerun1/HA/' # directory containing tree; results saved here
tree_file = 'flu_avian_h5n1_ha.json' # json tree you want to simulate mutations on
output_file = 'simulation_data.tsv' # output for dataframe with odds ratios and pvalues (.tsv)
host1 = 'Human'
host2 = 'Avian'
iterations = 10000



def read_in_tree_json(input_tree):
    '''read in a tree in json format'''
    with open(input_tree) as json_file:
        json_tree = json.load(json_file)
    return json_tree

def get_branch_lengths(parent, branch_length_dict = None, parenttotaldiv = None):
    '''return dictionary containing branch lengths of a parent branch and all descendents'''
    branch_length_dict = branch_length_dict or {}
    parenttotaldiv = parenttotaldiv or 0
    totaldiv = parent['node_attrs']['div']
    div = totaldiv - parenttotaldiv
    branch_length_dict[parent['name']] = div
    if 'children' in parent:
        for child in parent['children']:
            get_branch_lengths(child, branch_length_dict, totaldiv)
    return branch_length_dict

def simulate_gain_loss(branch_length, total_tree_branch_length):
    '''simulate whether or not to randomly mutate based on an individual branch length and the tree's total branch length;
    return bool'''
    rate = 1/total_tree_branch_length
    probability_stay_same = (1.0+exp(-2.0*branch_length*rate))/2.0
    # pick a random number between 0 and 1.0
    value = random.random()
    if value < probability_stay_same:
        mutation = False   # don't mutate
    else:
        mutation = True
    return mutation

def mutagenize_tree(parent, currentstate = None, branch_length_dict = None, total_branch_length = None):
    '''run simulate_gain_loss on all branches starting at the parent, recursively calling on all children;
    then return list with host and mutation state of all leaves'''
    ## get all branch lengths if this is the initial call
    branch_length_dict = branch_length_dict or get_branch_lengths(parent)
    total_branch_length = total_branch_length or sum(branch_length_dict.values())

    ## set currentstate to 0 if this is the initial call
    ## value of 0 = unmutated, 1 = mutated
    currentstate = currentstate or 0

    ## get length of the current branch
    branch_length = branch_length_dict[parent['name']]

    ## if mutation was simulated, flip currentstate
    if simulate_gain_loss(branch_length, total_branch_length):
        currentstate = 1 - currentstate

    ## if this is a leaf, append its name, host, and currentstate to global sim_results list
    if not 'children' in parent:
        sim_results.append((parent['name'], parent['node_attrs']['host']['value'], currentstate))
    
    ## if not, continue mutagenizing all downstream branches
    else:
        for child in parent['children']:
            mutagenize_tree(child, currentstate, branch_length_dict, total_branch_length)

def run_sims(iterations):
    '''perform n simulations, where n = number of iterations defined, and perform a Fisher's exact test for each iteration;
    then return a dict containing results for all simulations'''
    ## create dictionary of lists to append data to from each simulation iteration
    all_sim_data = {'oddsr': [], 'pvalue': [], 'host1count': [], 'host2count': []}

    for iter in range(iterations):
        ## create blank global list for each iteration to append data after mutagenizing the tree
        global sim_results 
        sim_results = []

        ## mutagenize from root, appending data for that iteration to sim_results
        mutagenize_tree(root)

        ## get 2x2 counts for Fisher's exact test
        p1 = sum([1 for _ in sim_results if _[1] == host1 and _[2] == 1])
        a1 = sum([1 for _ in sim_results if _[1] == host1 and _[2] == 0])
        p2 = sum([1 for _ in sim_results if _[1] == host2 and _[2] == 1])
        a2 = sum([1 for _ in sim_results if _[1] == host2 and _[2] == 0])

        ## add a pseudocount for denominator values that are equal to 0
        if p2 == 0:
            p2 = 1
        if a1 == 0:
            a1 = 1

        ## get odds ratio and pvalue from Fisher's exact test
        oddsr, p = fisher_exact([[p1, a1],[p2, a2]], alternative='two-sided')

        ## append this iteration to all_sim_data
        all_sim_data['oddsr'].append(oddsr)
        all_sim_data['pvalue'].append(p)
        all_sim_data['host1count'].append(p1)
        all_sim_data['host2count'].append(p2)
    
    return all_sim_data

def get_iteration_list(iterations, cores):
    '''return a list with total iterations evenly split among number of cores for multiprocessing'''
    div, rem = divmod(iterations, cores)
    return [div+1]*rem + [div]*(cores-rem)

if __name__ == '__main__':
    ## start timer
    start_time = time.time()

    ## get root of tree
    tree_path = json_dir + tree_file
    json_tree = read_in_tree_json(tree_path)
    root = json_tree['tree']

    ## split iterations among cores
    cores = mp.cpu_count()
    iter_list = get_iteration_list(iterations, cores)

    ## start multiprocessing pool
    mp.set_start_method('fork')
    pool = mp.Pool()

    ## and run the simulations
    pool_sim_data = pool.map(run_sims, iter_list)
    pool.close()
    pool.join()

    ## print timer statement
    print("It took", round(time.time() - start_time, 2), "seconds to run", iterations, "simulations")

    ## make dictionary with multiprocessing data as individual lists instead of lists of lists
    combined_sim_data = {'oddsr': [y for x in pool_sim_data for y in x['oddsr']],
                        'pvalue': [y for x in pool_sim_data for y in x['pvalue']],
                        'host1count': [y for x in pool_sim_data for y in x['host1count']],
                        'host2count': [y for x in pool_sim_data for y in x['host2count']]}

    ## generate dataframe from combined data
    output_df = pd.DataFrame({'iteration': range(1,iterations+1),
                            'oddsratio': combined_sim_data['oddsr'],
                            'pvalue': combined_sim_data['pvalue'],
                            host1.lower()+'count': combined_sim_data['host1count'],
                            host2.lower()+'count': combined_sim_data['host2count']})
    
    ## and save it
    output_path = json_dir + output_file
    output_df.to_csv(output_path, sep='\t', header=True, index=False)