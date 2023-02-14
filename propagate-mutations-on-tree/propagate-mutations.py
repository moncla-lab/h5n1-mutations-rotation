import json
import os

## specify directory, files, and gene
json_dir = '/Users/jort/Desktop/test/' # directory containing json files
tree_file = 'flu_avian_h5n1_pb2.json' # json tree you want to propagate mutations on
root_file = 'flu_avian_h5n1_pb2_root-sequence.json' # json file with inferred root sequence
output_file = 'flu_avian_h5n1_pb2_mutprop.json' # name for new json tree; to be saved in json_dir
gene = 'PB2' # gene to propagate muts from


##### user input above #####


## change to directory with files
os.chdir(json_dir)

## functions to read json files
def read_in_tree_json(input_tree):
    '''read in a tree in json format'''
    with open(input_tree) as json_file:
        json_tree = json.load(json_file)
    return json_tree

def get_root_seq(root_file, gene):
    '''get root aa sequence for gene of interest'''
    with open(root_file) as json_file:
        json_root_tree = json.load(json_file)
    return json_root_tree[gene]


## functions to prop muts and save new json tree
def get_position(mut):
   return(mut[1:len(mut)-1])

def assign_root_muts(tree, gene, root_seq):
    '''assign each aa in root sequence as a 'mutation' for the root node (formatted as '_1M')'''
    root_muts = []
    i = 1
    for aa in root_seq:
        root_muts.append(str("_"+str(i)+aa))
        i += 1
    tree['tree']['branch_attrs']['mutations'] = {gene: root_muts}

def propagate_muts(parent):
    '''propagate mutations from parent nodes onto children, excluding those at positions which are mutated in the child'''
    parent_muts = parent['branch_attrs']['mutations'][gene]
    children = parent['children']
    for child in children:
        if gene in child['branch_attrs']['mutations']:
            child_muts = child['branch_attrs']['mutations'][gene]
            muts_to_prop = []
            for pmut in parent_muts:
                if get_position(pmut) not in [get_position(x) for x in child_muts]:
                    muts_to_prop.append(pmut)
            child['branch_attrs']['mutations'][gene] = child['branch_attrs']['mutations'][gene] + muts_to_prop
        else:
            child['branch_attrs']['mutations'] = {gene: parent_muts}
        if 'NODE_' in child['name']:
            propagate_muts(child)


## read in tree, get root sequence, assign root muts, and propagate muts through tree
json_tree = read_in_tree_json(tree_file) #(tree_file, gene)
root_seq = get_root_seq(root_file, gene)
assign_root_muts(json_tree, gene, root_seq)
parent = json_tree['tree']
propagate_muts(parent)


## save json tree with propagated muts
with open(output_file, 'w') as outfile:
    json.dump(json_tree, outfile)