import os
import pandas as pd
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.similarity import ModifiedCosine
import sys
from tqdm import tqdm
from MolNotator.others.global_functions import *

def molnet(params : dict):

    # Load parameters
    mz_field = params['mz_field']
    rt_field = params['rt_field']

    in_path= params['mix_out_6_1']
    out_path_full= params['mix_out_7_1']
    out_path_samples= params['mix_out_7_2']
    mass_error= params['mn_mass_error']
    cosine_threshold= params['mn_cosine_threshold']
    matched_peaks= params['mn_matched_peaks']
    
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
    node_table = pd.read_csv(in_path + 'node_table.csv', index_col = "Index")
    edge_table = pd.read_csv(in_path + 'edge_table.csv', index_col = "Index")

    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
    if not os.path.isdir(out_path_samples) :
        os.mkdir(out_path_samples)

    #First, export simplified version of the network
    kept_edges_1 = edge_table[edge_table['status_universal'] == 'add_edge']
    kept_edges_2 = edge_table[edge_table['status_universal'] == 'cos_edge']
    kept_nodes = list(set(list(kept_edges_1['node_1']) + list(kept_edges_1['node_2'])))
    kept_nodes.sort()
    node_table = node_table.loc[kept_nodes]
    node_table[mz_field] = node_table[mz_field].round(4)
    edge_table = kept_edges_1.append(kept_edges_2, ignore_index = True)
    
    node_table.to_csv(out_path_full + "node_table_simple.csv", index_label = "Index")
    edge_table.to_csv(out_path_full + "edge_table_simple.csv", index_label = "Index")



    if params['process_mode'] == "POS":
        mgf_path_pos= params['pos_out_0']
        mgf_file_pos= params['pos_mgf']
        mgf = list(load_from_mgf(mgf_path_pos + mgf_file_pos))
        mgf = [Spectrum_processing(s) for s in mgf]
        molnet_single(node_table, edge_table, mgf, params['process_mode'], params)
        return
    elif params['process_mode'] == "NEG" :
        mgf_path_neg= params['neg_out_0']
        mgf_file_neg= params['neg_mgf']
        mgf = list(load_from_mgf(mgf_path_neg + mgf_file_neg))
        mgf = [Spectrum_processing(s) for s in mgf]
        molnet_single(node_table, edge_table, mgf, params['process_mode'], params)
        return


    mgf_path_neg= params['neg_out_0']
    mgf_path_pos= params['pos_out_0']
    mgf_file_neg= params['neg_mgf']
    mgf_file_pos= params['pos_mgf']

    mgf_pos = list(load_from_mgf(mgf_path_pos + mgf_file_pos))
    mgf_neg = list(load_from_mgf(mgf_path_neg + mgf_file_neg))    
    mgf_pos = [Spectrum_processing(s) for s in mgf_pos]
    mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
       
    # Do molecular families clusters:
    neutral_nodes = list(node_table.index[node_table['status'].str.contains('neutral')])
    neutral_nodes_2 = neutral_nodes.copy()
    neutral_edges_list = []
    total_neutrals = len(neutral_nodes)
    while len(neutral_nodes) > 1 :
        perc = round((1-(len(neutral_nodes)/total_neutrals))*100,1)
        sys.stdout.write("\rLinking neutrals : {0}%".format(perc))
        sys.stdout.flush()
        neutral_1 = neutral_nodes[0]
        adducts_1 = pd.DataFrame(index = list(edge_table['node_2'][edge_table['node_1'] == neutral_1]))
        adducts_1['spec_id'] = list(node_table.loc[adducts_1.index, "spec_id"].astype(int))
        adducts_1['ion_mode'] = list(node_table.loc[adducts_1.index, "ion_mode"])
        adducts_1['adduct'] = list(node_table.loc[adducts_1.index, "Adnotation"])
        neutral_nodes.remove(neutral_1)
        for neutral_2 in neutral_nodes:
            #neutral_2 = neutral_nodes[1]
            tmp_scores = list()
            adducts_2 = pd.DataFrame(index = list(edge_table['node_2'][edge_table['node_1'] == neutral_2]))
            adducts_2['spec_id'] = list(node_table.loc[adducts_2.index, "spec_id"].astype(int))
            adducts_2['ion_mode'] = list(node_table.loc[adducts_2.index, "ion_mode"])
            adducts_2['adduct'] = list(node_table.loc[adducts_2.index, "Adnotation"])
            for mode in adducts_1['ion_mode'].unique():
                #mode = adducts_1['ion_mode'].unique()[0]
                tmp_table_1 = adducts_1[adducts_1['ion_mode'] == mode]
                tmp_table_2 = adducts_2[adducts_2['ion_mode'] == mode]
                if mode == "POS" : mgf_list = mgf_pos
                else: mgf_list = mgf_neg
                for ion_1 in tmp_table_1.index:
                    spectrum_1 = mgf_list[tmp_table_1.loc[ion_1, "spec_id"]]
                    adduct_1 = tmp_table_1.loc[ion_1, 'adduct']
                    for ion_2 in tmp_table_2.index:
                        adduct_2 = tmp_table_2.loc[ion_2, 'adduct']
                        if adduct_1 != adduct_2 : continue
                        spectrum_2 = mgf_list[tmp_table_2.loc[ion_2, "spec_id"]]
                        score, n_matches = modified_cosine.pair(spectrum_1, spectrum_2)
                        if n_matches < matched_peaks : score = 0.0
                        tmp_scores.append((score, n_matches))
            tmp_scores = pd.DataFrame(tmp_scores, columns = ['cos', 'matches'])
            if len(tmp_scores) == 0 : continue
            if tmp_scores['cos'].max() >= cosine_threshold : 
                tmp_scores.sort_values(by = ['cos', 'matches'], ascending = False, inplace = True)
                neutral_edges_list.append((neutral_1, neutral_2, round(tmp_scores['cos'].iloc[0],2), tmp_scores['matches'].iloc[0]))
    
    # Add neutral_edges:
    print("Adding neutral edges...")
    for i in tqdm(range(len(neutral_edges_list))):
        new_idx = max(edge_table.index) + 1
        node_1 = neutral_edges_list[i][0]
        node_2 = neutral_edges_list[i][1]
        cos = round(neutral_edges_list[i][2], 2)
        matches = neutral_edges_list[i][3]
        node_1_mode = node_table.loc[node_1, 'ion_mode']
        node_2_mode = node_table.loc[node_2, 'ion_mode']
        if node_1_mode != node_2_mode or node_1_mode == "MIX":
            edge_mode = "MIX"
        else:
            edge_mode = node_1_mode
        rt_gap = abs(round(node_table.loc[node_1, rt_field] - node_table.loc[node_2, rt_field], 3))
        mz_gap = abs(round(node_table.loc[node_1, mz_field] - node_table.loc[node_2, mz_field], 4))
        edge_table.loc[new_idx] = [None]*len(edge_table.columns)
        edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks", "total_peaks",
                                 "matching_score", "rt_gap", "mz_gap", "status",
                                 "status_universal", "All_annotations", "ion_mode",
                                 "cosine_score"]] = [node_1, node_2, matches, 0, 0.0, rt_gap, mz_gap,
                      edge_mode + "_cosine_neutral_edge", "cosine_neutral_edge", cos, edge_mode, cos]
    
    # Delete undesired edges
    edge_table = edge_table[edge_table['status'].str.contains('cosine_neutral_edge')]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    
    # Add singletons:
    singletons = list(set(neutral_nodes_2) - set(kept_nodes))
    for neutral in tqdm(singletons):
        new_idx = max(edge_table.index) + 1 
        ion_mode = node_table.loc[neutral, "ion_mode"]
        edge_table.loc[new_idx] = [None]*len(edge_table.columns)
        edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks", "total_peaks",
                                 "matching_score", "rt_gap", "mz_gap", "status",
                                 "status_universal", "All_annotations", "ion_mode",
                                 "cosine_score"]] = [neutral, neutral, 0, 0, 0.0, 0.0, 0.0,
                      ion_mode + "_self_edge", "self_edge", 0.0, ion_mode, 0.0]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    node_table = node_table.loc[kept_nodes]
    
    node_table[mz_field] = node_table[mz_field].round(4)
    
    node_table.to_csv(out_path_full + "node_table_families.csv", index_label = "Index")
    edge_table.to_csv(out_path_full + "edge_table_families.csv", index_label = "Index")
    return
