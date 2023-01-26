# Simulate mutation gain and loss as a markov chain
# 
# Based off methodology described in [Pascoe et al](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.13051) and [Sheppard et al](https://www.pnas.org/doi/10.1073/pnas.1305559110)
# 
# In this notebook, we will construct a null distribution for the enrichment scores calculate in `calculate-enrichment-scores-across-tree`. To do so, we want to simulate the gain and loss of mutations as a random process that occurs across the tree. We will do this following the methodology outlined in Pascoe et al. Currently, the way I am doing this is as follows: 
# 
# 1. Starting with the tree in JSON format, make a copy of the tree object. 
# 2. On that copy, delete all existing mutations. Retain the tree's structuree (branch lengths, topology, tips) 
# 3. On that empty tree, traverse from root to tip. At each branch, evaluate whether to randomly add a mutation or not. This is evaluated probablistically as: `probability no mutation = (1.0+exp(-2.0*branch_length*rate))/2.0`, where `rate = 1/total_tree_branch_length`. This is exactly the same as what was described by Sheppard et al and is in his code. 
# 4. Perform a random draw of a number between 0 and 1. If that number is less than `probability no mutation`, do not mutate. If greater than, do mutate. 
# 5. If mutate, determine whether this is a gain of mutation or loss of mutation. To do this, we need to find the current state. To do so, traverse back up the tree to find the most recent ancestral node at which there is a mutation. If there is none, then our current state is unmutated, and this event is a gain eventt. If the current state is mutated, then this event will be a loss event. 
# 4. Add the proper mutation. Either a wild type to mutant (W1M) or mutant to wild type (M1W). Add this as an amino acid mutation on the tree. 
# 5. Return the simulated tree. Calculate enrichment score as: 
# 
# |host|presence|absence|
# |:------|:-------|:------|
# |host 1|A|B| 
# |host 2|C|D|
# 
# where A, B, C, and D are counts of the mutation's presence and absence in host 1 and host 2. The odds ratio is then calculated as: `OR = (A * D)/(B * C)`

from math import exp
import random

import calculate_enrichment_scores_across_tree_JSON as calenr
import tree_manager as tm



def simulate_gain_loss(branch_length, total_tree_branch_length):

    rate = 1/total_tree_branch_length
    
    probability_stay_same = (1.0+exp(-2.0*branch_length*rate))/2.0
    
    # pick a random number between 0 and 1.0
    value = random.random()
    
    if value < probability_stay_same:
        mutation = 0   # don't mutate
    else:
        mutation = 1
        
    return mutation



def return_most_recent_mutated_node(node, gene):
    """given an internal node, traverse back up the tree to find a parental node that has a mutation annotation.
    if you get to the root without finding a mutation, return root. This is necessary for determining the proper 
    starting state for the mutation you are adding"""
    
    if node.traits['branch_attrs']['mutations'][gene] == []:
        if node.parent !=None:
            parent_node = return_most_recent_mutated_node(node.parent, gene)
        else:
            parent_node = node
    else: 
        parent_node = node
    
    return parent_node



def simulate_gain_loss_as_markov_chain(tree, gene, total_tree_branch_length, branches_that_mutated, all_branches):

    # my fake mutation is going to be W1M for wild-type 1 mutant
    for k in tree.Objects:
        
        # collect branch length
        divergence = k.traits['node_attrs']['div']

        # if this happens at the root, set parent divergence to 0
        if k.parent.traits == {}:
            parent_div = 0
        elif 'node_attrs' not in k.parent.traits:
            parent_div = 0
        else:
            parent_div = k.parent.traits['node_attrs']['div']

        branch_length = divergence - parent_div
        
        # find the most recent mutated parent (this could be the parent node, grandparent, etc...)
        most_recent_mutated_parent = return_most_recent_mutated_node(k.parent, gene)
        parent_mut_state = most_recent_mutated_parent.traits['branch_attrs']['mutations'][gene]
        
        # given the length of the current branch and the total tree branch length, perform a random draw to 
        # decide whether to mutate. A result of 1 means mutate, 0 means do not mutate
        mutation = simulate_gain_loss(branch_length, total_tree_branch_length)

        if mutation == 1:  # if we've mutated

            if k.name in all_branches:
                all_branches[k.name]["times_mutated"] += 1
            else:
                all_branches[k.name] = {"branch_length":branch_length, "times_mutated":1}

            # add branch to dictionary for plotting later 
            if k.name in branches_that_mutated:
                branches_that_mutated[k.name]["times_mutated"] += 1

            else:
                branches_that_mutated[k.name] = {"branch_length":branch_length, "times_mutated":1}

            
            #print("we are mutating branch ", k)
            if parent_mut_state == [] or parent_mut_state == ['M1W']:
                k.traits['branch_attrs']['mutations'][gene] = ['W1M']
            elif parent_mut_state == ['W1M']:
                k.traits['branch_attrs']['mutations'][gene] = ['M1W']

        elif mutation == 0:
            if not k.name in all_branches:
                all_branches[k.name] = {"branch_length":branch_length, "times_mutated":0}
                
    return tree, branches_that_mutated, all_branches



def perform_simulations(pickled_tree, gene, total_tree_branch_length, host1, host2,host_annotation, min_required_count, host_counts, iterations):
    times_detected_all = {}
    branches_that_mutated = {}
    scores_dict_all = {}
    all_branches = {}
    
    for i in range(iterations):
        # We need to get a fresh copy of no_muts_tree from pickled_tree.
        # If not, then the new mutations will be appended to the original tree or to previous simulated trees.
        no_muts_tree = tm.get_clean_tree_copy(pickled_tree)
        sim_tree, branches_that_mutated, all_branches = simulate_gain_loss_as_markov_chain(no_muts_tree, gene, total_tree_branch_length, branches_that_mutated, all_branches)
        scores_dict, times_detected_dict, branch_lengths_dict, host_counts_dict2 = calenr.calculate_enrichment_scores(sim_tree, ['W1M'],['W1M'], host1, host2, host_annotation, min_required_count, host_counts, gene)
        times_detected_all[i] = times_detected_dict
        scores_dict_all[i] = scores_dict

    return scores_dict_all, times_detected_all, branches_that_mutated, all_branches