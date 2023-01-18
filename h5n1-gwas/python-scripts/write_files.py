import pickle
import json
import os
import config as cfg
import tree_manager as tm
from datetime import date

current_date = str(date.today())

def make_next_folder():
    folder_name = cfg.folder_path + cfg.naming_scheme
    i = 0
    while os.path.exists(folder_name+"_"+str(i)):
        i += 1
    folder_name = folder_name+"_"+str(i)
    os.mkdir(folder_name)
    data_folder_name = folder_name + "/data"
    os.mkdir(data_folder_name)
    return folder_name

def write_config(folder_name):
    config_dict = {
        'tree_path': cfg.tree_path,
        'baltic_path': cfg.baltic_path,
        'host_annotation': cfg.host_annotation,
        'gene': cfg.gene,
        'host1': cfg.host1,
        'host2': cfg.host2,
        'method': cfg.method,
        'minimum_required_count': cfg.minimum_required_count,
        'iterations':cfg.iterations
        }
    config_path = folder_name + "/config.txt"
    config_file = open(config_path, "w")
    for key in config_dict.keys():
        config_file.write(F"{key}: {config_dict[key]}")
        if not key == 'iterations':
            config_file.write("\n")
    config_file.close()

def write_baltic_tree(folder_name, baltic_tree):
    pickled_baltic_tree_path = folder_name + "/pickled_baltic_tree.obj"
    pickled_baltic_tree_file = open(pickled_baltic_tree_path, 'wb') 
    pickle.dump(baltic_tree, pickled_baltic_tree_file)
    pickled_baltic_tree_file.close()

def write_json_tree(folder_name, json_tree):
    pickled_json_tree_path = folder_name + "/pickled_json_tree.obj"
    pickled_json_tree_file = open(pickled_json_tree_path, 'wb')
    pickle.dump(json_tree, pickled_json_tree_file)
    pickled_json_tree_file.close()

def write_dfs(folder_name, df5, df8, df9 = None, df10 = None):
    output_filename = folder_name + "/data/" + cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_data_" + current_date + ".tsv"
    df5.to_csv(output_filename, sep="\t", header=True, index_label="mutation")

    output_filename = folder_name + "/data/" + cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_simulated_" + current_date + ".tsv"
    df8.to_csv(output_filename, sep="\t", header=True, index=False)

    if cfg.mode == 'testing':
        output_filename = folder_name + "/data/" + cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_simulated_lengthVStimes_" + current_date + ".tsv"
        df9.to_csv(output_filename, sep="\t", header=True, index=False)

        output_filename = folder_name + "/data/" + cfg.gene + "_" + cfg.host1 + "_vs_" + cfg.host2 + "_simulated_all_branches_" + current_date + ".tsv"
        df10.to_csv(output_filename, sep="\t", header=True, index=False)