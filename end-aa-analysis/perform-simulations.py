from math import exp
import random
import json
import pandas as pd
from scipy.stats import fisher_exact


json_dir = '/Users/jort/Desktop/test/build3_pb2_end_aa/' # directory containing tree; results saved here
tree_file = 'flu_avian_h5n1_pb2.json' # json tree you want to simulate mutations on
output_file = 'simulation_data.tsv' # output for dataframe with odds ratios and pvalues (.tsv)
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



tree_path = json_dir + tree_file
json_tree = read_in_tree_json(tree_path)
parent = json_tree['tree']

div_dict = {}
get_div(parent)
total_branch_length = sum(div_dict.values())

all_sim_data = {'oddsr': [], 'pvalue': [], 'humancount': [], 'aviancount': []}
for x in range(iterations):
    sim_results = []
    mutagenize_tree(parent)
    p1 = sum([1 for _ in sim_results if _[1] == 'Human' and _[2] == 1])
    a1 = sum([1 for _ in sim_results if _[1] == 'Human' and _[2] == 0])
    p2 = sum([1 for _ in sim_results if _[1] == 'Avian' and _[2] == 1])
    a2 = sum([1 for _ in sim_results if _[1] == 'Avian' and _[2] == 0])
    if p2 == 0:
        p2 = 1
    if a1 == 0:
        a1 = 1
    oddsr, p = fisher_exact([[p1, a1],[p2, a2]], alternative='two-sided')
    all_sim_data['oddsr'].append(oddsr)
    all_sim_data['pvalue'].append(p)
    all_sim_data['humancount'].append(p1)
    all_sim_data['aviancount'].append(p2)

output_path = json_dir + output_file
output_df = pd.DataFrame({'iteration': range(1,iterations+1), 'oddsratio': all_sim_data['oddsr'], 'pvalue': all_sim_data['pvalue'], 'humancount': all_sim_data['humancount'], 'aviancount': all_sim_data['aviancount']})
output_df.to_csv(output_path, sep='\t', header=True, index=False)