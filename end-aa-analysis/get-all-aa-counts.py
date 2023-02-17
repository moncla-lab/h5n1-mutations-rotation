import json
import imp
import pandas as pd

## specify directory, files, hosts, and gene
baltic_path = '/Users/jort/coding/baltic/baltic/baltic.py' # path to baltic.py file (https://github.com/evogytis/baltic)
gene = 'HA' # gene to analyze
json_dir = f'/Users/jort/coding/h5n1-mutations-rotation/build3-end-aa-analysis/build3-rerun3/{gene.lower()}/' # directory containing json file; output files to be saved here
tree_file = f'flu_avian_h5n1_{gene.lower()}_mutprop.json' # json tree with mutations propagated
host1 = 'Human' # host to screen for enrichment in
host2 = 'Avian' # background host
host1_output_file = 'human_all_aa_counts.tsv' # output for host1 dataframe (.tsv)
host2_output_file = 'avian_all_aa_counts.tsv' # output for host2 dataframe (.tsv)


##### user input above #####


## import baltic
bt = imp.load_source('baltic', baltic_path)

def read_in_tree_json(input_tree):
    '''read in a tree in json format'''
    with open(input_tree) as json_file:
        json_tree = json.load(json_file)
    json_translation = {'absoluteTime':lambda k: k.traits['node_attrs']['num_date']['value'],'name':'name'} ## allows baltic to find correct attributes in JSON, height and name are required at a minimum
    bt_tree, meta = bt.loadJSON(json_tree, json_translation)
    return json_tree, bt_tree






## read in json tree
json_file = json_dir + tree_file
json_tree, tree = read_in_tree_json(json_file)

## get full AA sequence for each leaf
all_seqs = []
for k in tree.Objects:
  seq_dict = {}
  if k.branchType == 'leaf':
    host = k.traits['node_attrs']['host']['value']
    for mut in k.traits['branch_attrs']['mutations'][gene]:
      seq_dict[int(mut[1:len(mut)-1])] = mut[-1]
    all_seqs.append([host, seq_dict])


## generate host count dicts for each amino acid at each position
if gene == 'HA':
  aalength = 568 # 568 AAs in HA
elif gene == 'PB2':
  aalength = 759 # 759 AAs in PB2
all_positions = range(1, aalength+1)
all_aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'other']
host1_dict = dict((pos, dict((aa, 0) for aa in all_aas)) for pos in all_positions)
host2_dict = dict((pos, dict((aa, 0) for aa in all_aas)) for pos in all_positions)

## iterate through each sequence and record counts of AA occurrences in each host
for [host, seq_dict] in all_seqs:
  if host == host1:
    for k, v in seq_dict.items():
      if v in all_aas and k in list(host1_dict):
        host1_dict[k][v] += 1
      elif v not in all_aas and k in list(host1_dict):
        host1_dict[k]['other'] += 1
  elif host == host2:
    for k, v in seq_dict.items():
      if v in all_aas and k in list(host2_dict):
        host2_dict[k][v] += 1
      elif v not in all_aas and k in list(host2_dict):
        host2_dict[k]['other'] += 1

## save host count dataframes as tsv files
host1_output_path = json_dir + host1_output_file
host2_output_path = json_dir + host2_output_file
pd.DataFrame.from_dict(host1_dict).to_csv(host1_output_path, sep='\t', header=True, index=True)
pd.DataFrame.from_dict(host2_dict).to_csv(host2_output_path, sep='\t', header=True, index=True)