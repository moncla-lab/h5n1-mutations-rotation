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

def get_div(parent, parenttotaldiv = None):
    parenttotaldiv = parenttotaldiv or 0
    totaldiv = parent['node_attrs']['div']
    div = totaldiv - parenttotaldiv
    div_dict[parent['name']] = div
    if 'children' in parent:
        for child in parent['children']:
            get_div(child, totaldiv)

def simulate_gain_loss(branch_length, total_tree_branch_length):
    rate = 1/total_tree_branch_length
    probability_stay_same = (1.0+exp(-2.0*branch_length*rate))/2.0
    # pick a random number between 0 and 1.0
    value = random.random()
    if value < probability_stay_same:
        mutation = False   # don't mutate
    else:
        mutation = True
    return mutation

def mutagenize_tree(parent, currentstate = None):
    currentstate = currentstate or 0
    branch_length = div_dict[parent['name']]
    if simulate_gain_loss(branch_length, total_branch_length):
        currentstate = 1 - currentstate
    if not 'children' in parent:
        sim_results.append((parent['name'], parent['node_attrs']['host']['value'], currentstate))
    else:
        for child in parent['children']:
            mutagenize_tree(child, currentstate)

def run_sims(iterations):
    all_sim_data = {'oddsr': [], 'pvalue': [], 'host1count': [], 'host2count': []}
    for x in range(iterations):
        global sim_results
        sim_results = []
        mutagenize_tree(parent)
        p1 = sum([1 for _ in sim_results if _[1] == host1 and _[2] == 1])
        a1 = sum([1 for _ in sim_results if _[1] == host1 and _[2] == 0])
        p2 = sum([1 for _ in sim_results if _[1] == host2 and _[2] == 1])
        a2 = sum([1 for _ in sim_results if _[1] == host2 and _[2] == 0])
        if p2 == 0:
            p2 = 1
        if a1 == 0:
            a1 = 1
        oddsr, p = fisher_exact([[p1, a1],[p2, a2]], alternative='two-sided')
        all_sim_data['oddsr'].append(oddsr)
        all_sim_data['pvalue'].append(p)
        all_sim_data['host1count'].append(p1)
        all_sim_data['host2count'].append(p2)
    return all_sim_data

def get_iteration_list(iterations, cores):
    """Function to create a list with total iterations evenly split among number of cores for multiprocessing"""
    div, rem = divmod(iterations, cores)
    return [div+1]*rem + [div]*(cores-rem)

if __name__ == '__main__':
    start_time = time.time()

    tree_path = json_dir + tree_file
    json_tree = read_in_tree_json(tree_path)
    parent = json_tree['tree']

    div_dict = {}
    get_div(parent)
    total_branch_length = sum(div_dict.values())

    cores = mp.cpu_count()
    iter_list = get_iteration_list(iterations, cores)
    mp.set_start_method('fork')
    pool = mp.Pool()
    sim_data = pool.map(run_sims, iter_list) # Run simmut.perform_simulations using arguments specified in sim_part, with iterations split among cores as specified in iter_list
    pool.close()
    pool.join()

    print("It took", round(time.time() - start_time, 2), "seconds to run", iterations, "simulations")

    combined_sim_data = {'oddsr': [y for x in sim_data for y in x['oddsr']],
                        'pvalue': [y for x in sim_data for y in x['pvalue']],
                        'host1count': [y for x in sim_data for y in x['host1count']],
                        'host2count': [y for x in sim_data for y in x['host2count']]}

    output_path = json_dir + output_file
    output_df = pd.DataFrame({'iteration': range(1,iterations+1),
                            'oddsratio': combined_sim_data['oddsr'],
                            'pvalue': combined_sim_data['pvalue'],
                            host1.lower()+'count': combined_sim_data['host1count'],
                            host2.lower()+'count': combined_sim_data['host2count']})
    output_df.to_csv(output_path, sep='\t', header=True, index=False)