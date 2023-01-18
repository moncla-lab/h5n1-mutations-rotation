import copy
import pickle
import json
import imp

import config as cfg

if cfg.baltic_path == None or cfg.baltic_path == "pip":
    import baltic as bt
else:
    bt = imp.load_source('baltic', cfg.baltic_path)



def init_trees(tree_path = None, gene = None):
    tree_path = tree_path or cfg.tree_path
    gene = gene or cfg.gene
    json_tree, orig_tree = read_in_tree_json(tree_path)
    no_muts_tree = get_no_muts_tree(orig_tree, gene)
    pickled_tree = pickle.dumps(no_muts_tree, -1)
    return json_tree, orig_tree, no_muts_tree, pickled_tree



def read_in_tree_json(input_tree):
    """read in a tree in json format"""

    with open(input_tree) as json_file:
        json_tree = json.load(json_file)

    # Nextstrain tree jsons are broken into 2 components: metadata and tree. We just need the tree
    json_translation = {'absoluteTime':lambda k: k.traits['node_attrs']['num_date']['value'],'name':'name'} ## allows baltic to find correct attributes in JSON, height and name are required at a minimum
    tree, meta = bt.loadJSON(json_tree, json_translation)

    return json_tree, tree



def get_no_muts_tree(orig_tree, gene):
    """this function deletes all amino acid mutations from the tree"""
    
    # we need to make a copy, otherwise this will alter the no muts tree
    tree = copy.deepcopy(orig_tree)
    
    for k in tree.Objects:
        if 'mutations' not in k.traits['branch_attrs']:
            k.traits['branch_attrs']['mutations'] = {}

        if gene not in k.traits['branch_attrs']['mutations']:
            k.traits['branch_attrs']['mutations'][gene] = []
        else:
            k.traits['branch_attrs']['mutations'][gene] = []    
            
        if 'branch_attrs' not in k.parent.traits:
            k.parent.traits = {'branch_attrs':{}}
        if 'mutations' not in k.parent.traits['branch_attrs']:
            k.parent.traits['branch_attrs']['mutations'] = {} 
        if gene not in k.parent.traits['branch_attrs']['mutations']:
            k.parent.traits['branch_attrs']['mutations'][gene] = []
        else:
            k.parent.traits['branch_attrs']['mutations'][gene] = []   
    return tree



def get_clean_tree_copy(pickled_tree):
    return pickle.loads(pickled_tree)