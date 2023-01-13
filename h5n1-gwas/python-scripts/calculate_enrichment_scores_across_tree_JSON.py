#!/usr/bin/env python
# coding: utf-8

# # Calculate enrichment scores across tree
# 
# This is an initial attempt to break this code into reasonable chunks that can be coded, tested, and executed separately. In this notebook, we will be reading in a tree, enumerating all the mutations on that tree, and calculating enrichment scores for each of them. These enrichment scores are based on [this code](https://github.com/sheppardlab/pGWAS/blob/master/assomap_given_phylo.py), written for [this paper](https://www.nature.com/articles/s41467-018-07368-7#Sec10) detailed in lines 245-273 and calculated from the following contingency table as: 
# 
# |host|presence|absence|
# |:------|:-------|:------|
# |host 1|A|B| 
# |host 2|C|D|
# 
# where A, B, C, and D are counts of the mutation's presence and absence in host 1 and host 2. The odds ratio is then calculated as: `OR = (A * D)/(B * C)`
# 
# In this notebook, this code is written for parsing a tree json format, output from Nextstrain. In subsequent notebooks, I will alter this for running on beast trees. 

# ## Infer each mutation that occurs across the tree


import json

import imp
bt = imp.load_source('baltic', '/Users/jort/coding/baltic/baltic/baltic.py')

"""read in a tree in json format"""
def read_in_tree_json(tree_path):
    with open(tree_path) as json_file:
        tree_json = json.load(json_file)

    # Nextstrain tree jsons are broken into 2 components: metadata and tree. We just need the tree
    tree_object=tree_json['tree']
    #meta=tree_json['meta']
    json_translation={'absoluteTime':lambda k: k.traits['node_attrs']['num_date']['value'],'name':'name'} ## allows baltic to find correct attributes in JSON, height and name are required at a minimum

    #tree=bt.loadJSON(tree_object,json_translation)
    tree,meta=bt.loadJSON(tree_json,json_translation)
    
    return(tree)



"""count the number of tips on the tree corresponding to each host category"""
def return_all_host_tips(tree, host1, host2, host_annotation):
    host_counts = {host1:0, host2:0, "other":0}
    
    for k in tree.Objects: 
        if k.branchType == "leaf":
            host = k.traits['node_attrs'][host_annotation]['value']
            if host in [host1, host2]:
                host_counts[host] += 1
            else:
                host_counts['other'] += 1
    return(host_counts)




"""this function traverses the tree, from root to tip, and records 2 quantites:
1. for each branch, it records the branch name and its branch length in a 
dictionary. When parsing a tree in baltic, each branch gets assigned a random, but unique name/identifier
2. it adds up the total branch length on the tree as an integer"""

def return_total_tree_branch_length(tree):
    total_branch_length = 0
    branch_lengths = {}
    
    for k in tree.Objects:
        divergence = k.traits['node_attrs']['div']

        # if this happens at the root, set parent divergence to 0
        if k.parent.traits == {}:
            parent_div = 0
        else:
            parent_div = k.parent.traits['node_attrs']['div']
        
        total_branch_length += divergence - parent_div
        branch_lengths[k] = divergence-parent_div
    
    return(total_branch_length, branch_lengths)




"""Traverse the tree from root to tip. On each branch, tips and nodes, gather every nucleotide and amino acid 
mutation"""

def gather_all_mut_on_tree(tree, gene):
    
    all_aa_muts = []
    all_nt_muts = []
    
    for k in tree.Objects:
        if 'mutations' in k.traits['branch_attrs']:
            if 'nuc' in k.traits['branch_attrs']['mutations']:
                nt_muts = k.traits['branch_attrs']['mutations']['nuc']
                all_nt_muts.extend(nt_muts)
                
            if gene in k.traits['branch_attrs']['mutations']:                
                aa_muts = k.traits['branch_attrs']['mutations'][gene]
                all_aa_muts.extend(aa_muts)
    
    # subset to only unique mutations 
    all_aa_muts = list(set(all_aa_muts))
    all_nt_muts = list(set(all_nt_muts))
    
    return(all_aa_muts, all_nt_muts)




"""return the total number of times that the mutation arises on the phylogeny. This includes instances of mutation 
on internal nodes and on tips and counts each with the same weight"""

def return_number_times_on_tree(tree, mut, gene):
    times_on_tree = 0
    
    for k in tree.Objects:
        if 'mutations' in k.traits['branch_attrs']:
        
            if gene in k.traits['branch_attrs']['mutations']:
                aa_muts = k.traits['branch_attrs']['mutations'][gene]
                if mut in aa_muts: 
                    times_on_tree += 1
                                                
    return(times_on_tree)




"""return the total branch length for all branches on which the mutation arises, added together. This includes 
instances of mutation on internal nodes and on tips"""

def return_branch_length_mut_on_tree(tree, mut, gene):
    branch_length = 0
            
    for k in tree.Objects:
        if 'mutations' in k.traits['branch_attrs']:
            if gene in k.traits['branch_attrs']['mutations']:
                aa_muts = k.traits['branch_attrs']['mutations'][gene]
                
                # if the mutation arises on this branch
                if mut in aa_muts: 
                    divergence = k.traits['node_attrs']['div']
                    
                    # if this happens at the root, set parent divergence to 0
                    if k.parent.traits == {}:
                        parent_div = 0
                    elif 'node_attrs' not in k.parent.traits:
                        parent_div = 0
                    else:
                        parent_div = k.parent.traits['node_attrs']['div']
                
                    local_branch_length = divergence - parent_div
                    branch_length += local_branch_length
                    
    return(branch_length)




"""given a branch and gene, return the mutations present on that branch"""

def return_muts_on_branch(branch, gene):
    muts = []
    if 'branch_attrs' in branch.traits:
        if 'mutations' in branch.traits['branch_attrs']:
            if gene in branch.traits['branch_attrs']['mutations']:
                muts = branch.traits['branch_attrs']['mutations'][gene]
                            
    return(muts)



"""Given a starting internal node, and a tip you would like to end at, traverse the full path from that node to
tip. Along the way, gather mutations that occur along that path. Once you have reached the ending 
tip, return the list of mutations that fell along that path"""

def return_all_muts_on_path_to_tip(starting_node, ending_tip, gene, muts, host_annotation):
    
    # set an empty list of mutations and enumerate the children of the starting node; children can be tips or nodes
    children = starting_node.children
    
    for child in children:
        local_muts = []
        
        """if the child is a leaf: if leaf is the target end tip, add the mutations that occur on that branch to 
        the list and return the list; if leaf is not the target end tip, move on"""
        """if the child is an internal node: first, test whether that child node contains the target tips in its 
        children. child.leaves will output a list of the names of all tips descending from that node. If not, pass. 
        if the node does contain the target end tip in its leaves, keep traversing down that node recursively, 
        collecting mutations as you go"""

        if child.branchType == "leaf":
            if child.name != ending_tip:
                pass
            elif child.name == ending_tip:
                host = child.traits['node_attrs'][host_annotation]['value']
                local_muts = return_muts_on_branch(child, gene)
                muts.extend(local_muts)
                #print(child.name, host, muts)
                return(host, muts)        
        
        elif child.branchType == "node":
            if ending_tip not in child.leaves:
                pass
            else:
                local_muts = return_muts_on_branch(child, gene)
                muts.extend(local_muts)
                host, muts = return_all_muts_on_path_to_tip(child, ending_tip, gene, muts,host_annotation)
    
    return(host, muts)




"""at times, will need to check whether the revertant mutation occcurs downstream. Return the revertant mutation"""

def return_opposite_mutation(mut):
    
    site = mut[1:-1]
    ref = mut[0]
    alt = mut[-1]
    opposite = alt+site+ref
    
    return(opposite)




"""given a tree, mutation, and gene, return the number of times that mutation is present in each host"""

def return_host_distribution_mutation(tree, mut, gene, host1, host2, host_annotation):
    
    host_counts_dict = {host1:0, host2:0, "other":0}
    back_mutation = return_opposite_mutation(mut)
    
    # iterate through tree
    for k in tree.Objects:
        if 'mutations' in k.traits['branch_attrs']:
            if gene in k.traits['branch_attrs']['mutations']:
                aa_muts = k.traits['branch_attrs']['mutations'][gene]  
                
                # if we have reached a node or tip in the tree with the target mutation, enumerate descendants
                if mut in aa_muts: 
                    
                    # if the mutation occurs on a leaf, record the host and move on 
                    if k.branchType == 'leaf':
                        host = k.traits['node_attrs'][host_annotation]['value']
                        if host in [host1,host2]:
                            host_counts_dict[host] += 1
                        else:
                            host_counts_dict["other"]+= 1
                    
                    # else, if the mutation occurs on a node, traverse the children and return host
                    elif k.branchType == "node":
                        all_leaves = k.leaves
                        for leaf in all_leaves: 
                            muts = []
                            host, muts = return_all_muts_on_path_to_tip(k, leaf, gene, muts, host_annotation)
                            if back_mutation in muts: 
                                pass
                            elif back_mutation and mut in muts:  # if both the mutation and backmutation occur, print
                                print(leaf, back_mutation, mut)
                            else:
                                if host in [host1, host2]:
                                    host_counts_dict[host] += 1
                                else: 
                                    host_counts_dict["other"] += 1
                                
    return(host_counts_dict)






# ## Calculate the enrichment scores



"""calculate an enrichment score for an individual mutation, based on the counts across hosts"""

def calculate_enrichment_score_counts(mut_counts_dict, host1, host2, host_counts):
    total_host1_tree = host_counts[host1]
    total_host2_tree = host_counts[host2]
    
    mut_host1 = mut_counts_dict[host1]
    mut_host2 = mut_counts_dict[host2]
    
    # this is calculating this table as counts
    presence_host1 = mut_host1
    absence_host1 = total_host1_tree - mut_host1
    presence_host2 = mut_host2
    absence_host2 = total_host2_tree - mut_host2
    
    if presence_host2 == 0:
        presence_host2 = 1
    if absence_host1 == 0:
        absence_host1 = 1

    # this score is calculated in terms of its enrichment in host 1
    score = (presence_host1 * absence_host2)/(presence_host2 * absence_host1)
    from scipy.stats import fisher_exact
    oddsr, p = fisher_exact([[presence_host1, absence_host1],[presence_host2, absence_host2]], alternative='two-sided')
    
    return(oddsr, p)




"""calculate an enrichment score for an individual mutation, based on the counts across hosts"""

def calculate_enrichment_score_proportions(mut_counts_dict, host1, host2, host_counts):
    total_host1_tree = host_counts[host1]
    total_host2_tree = host_counts[host2]
    
    mut_host1 = mut_counts_dict[host1]
    mut_host2 = mut_counts_dict[host2]
    
    total_tips_in_tree = total_host1_tree + total_host2_tree
    
    # this is calculating this table as proportions
    presence_host1 = (mut_host1)/total_tips_in_tree
    absence_host1 = (total_host1_tree - mut_host1)/total_tips_in_tree
    presence_host2 = (mut_host2)/total_tips_in_tree
    absence_host2 = (total_host2_tree - mut_host2)/total_tips_in_tree
    
#     if presence_host2 == 0:
#         presence_host2 = 1
#     if absence_host1 == 0:
#         absence_host1 = 1

    # this score is calculated in terms of its enrichment in host 1
    #score = (presence_host1 * absence_host2)/(presence_host2 * absence_host1)
    score = (presence_host1 + absence_host2) - (presence_host2 + absence_host1)
    return(score)




"""for a tree and all amino acid mutations, calculate the enrichment scores across the tree"""
def calculate_enrichment_scores(tree, aa_muts, nt_muts, host1, host2, host_annotation, min_required_count, method, host_counts, gene):
    scores_dict = {}
    times_detected_dict = {}
    branch_lengths_dict = {}
    host_counts_dict2 = {}
    
    if method == "counts":
        enrichment_calculation_function = calculate_enrichment_score_counts
    elif method == "proportions":
        enrichment_calculation_function = calculate_enrichment_score_proportions
    
    for a in aa_muts:
        times_detected = return_number_times_on_tree(tree, a, gene)
        times_detected_dict[a] = times_detected

        branch_length_mut = return_branch_length_mut_on_tree(tree, a, gene)
        branch_lengths_dict[a] = branch_length_mut

        host_counts_dict = return_host_distribution_mutation(tree, a, gene, host1, host2, host_annotation)
        host_counts_dict2[a] = host_counts_dict
        total_tips_with_mut = host_counts_dict[host1] + host_counts_dict[host2] + host_counts_dict["other"]

        if total_tips_with_mut >= min_required_count:
            enrichment_score,p_value = enrichment_calculation_function(host_counts_dict, host1, host2, host_counts)
            #print(enrichment_score, total_tips_with_mut, host_counts_dict)
            scores_dict[a] = {"enrichment_score": enrichment_score, "pvalue":p_value}
            
            
    return(scores_dict, times_detected_dict, branch_lengths_dict, host_counts_dict2)











