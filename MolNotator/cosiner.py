import os
import pandas as pd
import sys
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.similarity import ModifiedCosine
from tqdm import tqdm
from MolNotator.others.global_functions import *

def cosiner(params : dict):
    
    # Load parameters
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    idx_column = params['index_col']
    in_path= params['mix_out_5_1']
    out_path_full= params['mix_out_6_1']
    out_path_samples= params['mix_out_6_2']
    purge_empty_spectra = params['c_purge_empty_spectra']
    
    mass_error= params['c_mass_error']
    cosine_threshold= params['c_lowcos_threshold']
    cosiner_threshold= params['c_hardcos_threshold']
    matched_peaks= params['c_matched_peaks']
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
    # mkdir
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
    if not os.path.isdir(out_path_samples) :
        os.mkdir(out_path_samples)


    
    # if single mode NEG
    if params['process_mode'] == "NEG":
        
        mzmine_path_neg= params['neg_out_0']
        neg_csv_file= params['neg_csv']
        neg_mgf_file= params['neg_mgf']
        
        # Load node and edge tables
        node_table = pd.read_csv(in_path + "node_table.csv", index_col = "Index")
        edge_table  = pd.read_csv(in_path + "edge_table.csv", index_col = "Index")
        
        mgf = list(load_from_mgf(mzmine_path_neg + neg_mgf_file))
        mgf = [Spectrum_processing(s) for s in mgf]
        
        # Make a Series with MGF indexes as data and feature IDs as indexes
        mgf_data = pd.Series(dtype = int)
        for i in range(len(mgf)):
            mgf_data.loc[int(mgf[i].get(idx_column))] = i
        
        cosiner_single(node_table, edge_table, mgf, mgf_data, params['process_mode'], params)
        return
        
    # If single mode POS
    elif params['process_mode'] == "POS":
        
        mzmine_path_pos= params['pos_out_0']
        pos_csv_file= params['pos_csv']
        pos_mgf_file= params['pos_mgf']
        
        # Load node and edge tables
        node_table = pd.read_csv(in_path + "node_table.csv", index_col = "Index")
        edge_table  = pd.read_csv(in_path + "edge_table.csv", index_col = "Index")
        
        
        mgf = list(load_from_mgf(mzmine_path_pos + pos_mgf_file))
        mgf = [Spectrum_processing(s) for s in mgf]

        # Make a Series with MGF indexes as data and feature IDs as indexes
        mgf_data = pd.Series(dtype = int)
        for i in range(len(mgf)):
            mgf_data.loc[int(mgf[i].get(idx_column))] = i

        cosiner_single(node_table, edge_table, mgf, mgf_data, params['process_mode'], params)
        return


    # If double mode ("BOTH")
    else:
        mzmine_path_neg= params['neg_out_0']
        mzmine_path_pos= params['pos_out_0']
        neg_csv_file= params['neg_csv']
        pos_csv_file= params['pos_csv']
        neg_mgf_file= params['neg_mgf']
        pos_mgf_file= params['pos_mgf']
        neg_mgf = list(load_from_mgf(mzmine_path_neg + neg_mgf_file))
        pos_mgf = list(load_from_mgf(mzmine_path_pos + pos_mgf_file))
        neg_mgf = [Spectrum_processing(s) for s in neg_mgf]
        pos_mgf = [Spectrum_processing(s) for s in pos_mgf]
    
        # Make a Series with MGF indexes as data and feature IDs as indexes
        neg_mgf_data = pd.Series(dtype = int)
        for i in range(len(neg_mgf)):
            neg_mgf_data.loc[int(neg_mgf[i].get(idx_column))] = i
        
        pos_mgf_data = pd.Series(dtype = int)
        for i in range(len(pos_mgf)):
            pos_mgf_data.loc[int(pos_mgf[i].get(idx_column))] = i
    
    
    # Load node and edge tables
    
    node_table = pd.read_csv(in_path + "node_table.csv", index_col = "Index")
    edge_table  = pd.read_csv(in_path + "edge_table.csv", index_col = "Index")
    
    # List the molecular clusters (clusters with at least one neutral node)
    cluster_list = []
    cluster_list_neg = []
    cluster_list_pos = []
    print('Finding molecular clusters...')
    for i in tqdm(node_table['cluster_id'].unique()):
        if sum(node_table['status_universal'][node_table['cluster_id'] == i] == "neutral") > 0:
            cluster_list.append(i)
            if "NEG" in node_table['ion_mode'][node_table['cluster_id'] == i].unique():
                cluster_list_neg.append(i)
            if "POS" in node_table['ion_mode'][node_table['cluster_id'] == i].unique():
                cluster_list_pos.append(i)
    cluster_list.sort()
    cluster_list_neg.sort()
    cluster_list_pos.sort()
    
    # Cluster singletons and precursor ions in non-molecular clusters to ions in molecular clusters
    unclustered_ions = [i for i in node_table.index if node_table.loc[i, "cluster_id"] not in cluster_list]
    remains_ions_neg = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_singleton'])
    remains_ions_neg +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_precursor'])
    remains_ions_pos = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_singleton'])
    remains_ions_pos +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_precursor'])
    
    # Start processing NEG data
    neg_node_table = node_table[node_table['ion_mode'] == "NEG"]
    cluster_ion_list = pd.Series(index = cluster_list_neg, dtype = str)
    for i in cluster_ion_list.index:
        tmp_rows = neg_node_table.index[neg_node_table['cluster_id'] == i]
        cluster_ion_list[i] = '|'.join(neg_node_table.loc[tmp_rows, 'spec_id'].dropna().astype(int).astype(str))
    cluster_ion_list.sort_index(inplace = True)
    full_node_1 = []
    full_node_2 = []
    full_cluster_ids = []
    full_cosine = []
    full_matches = []
    for i in tqdm(remains_ions_neg):
        i_id = int(node_table.loc[i, idx_column])
        i_spectrum = neg_mgf[neg_mgf_data[i_id]]
        tmp_cluster_list = []
        tmp_id_list = []
        tmp_cosine_list = []
        tmp_prod_list = []
        tmp_match_list = []
        for j in cluster_ion_list.index:
            ion_list = list(map(int, cluster_ion_list[j].split('|')))
            cos_list = list()
            match_list = list()
            id_list = list()
            prod_list = list()
            for k in ion_list:
                score, n_matches = modified_cosine.pair(i_spectrum, neg_mgf[k])
                id_list.append(int(neg_mgf[k].get(idx_column)))
                cos_list.append(score)
                match_list.append(n_matches)
                prod_list.append(score * n_matches)
            tmp_prod_list.append(max(prod_list))
            tmp_cosine_list.append(cos_list[prod_list.index(max(prod_list))])
            tmp_match_list.append(match_list[prod_list.index(max(prod_list))])
            tmp_id_list.append(id_list[prod_list.index(max(prod_list))])
            tmp_cluster_list.append(j)
        tmp_table = pd.DataFrame(list(zip(tmp_id_list, tmp_cluster_list, tmp_cosine_list, tmp_prod_list, tmp_match_list)),
                                 columns = [f'ion_{idx_column}', 'cluster_id', 'cosine', 'prod', 'matches'])
        tmp_table = tmp_table.loc[tmp_table['prod'].idxmax()]
        
        counter_feature_id = tmp_table[f'ion_{idx_column}']
        counter_idx = neg_node_table.index[neg_node_table[idx_column] == counter_feature_id][0]
        full_node_1.append(i)
        full_node_2.append(counter_idx)
        full_cluster_ids.append(tmp_table['cluster_id'])
        full_cosine.append(round(tmp_table['cosine'], 2))
        full_matches.append(tmp_table['matches'])
    
    # Start processing POS data
    pos_node_table = node_table[node_table['ion_mode'] == "POS"]
    cluster_ion_list = pd.Series(index = cluster_list_pos, dtype = str)
    for i in cluster_ion_list.index:
        tmp_rows = pos_node_table.index[pos_node_table['cluster_id'] == i]
        cluster_ion_list[i] = '|'.join(pos_node_table.loc[tmp_rows, 'spec_id'].dropna().astype(int).astype(str))
    cluster_ion_list.sort_index(inplace = True)
    for i in tqdm(remains_ions_pos):
        i_id = int(node_table.loc[i, idx_column])
        i_spectrum = pos_mgf[pos_mgf_data[i_id]]
        tmp_cluster_list = []
        tmp_id_list = []
        tmp_cosine_list = []
        tmp_prod_list = []
        tmp_match_list = []
        for j in cluster_ion_list.index:
            ion_list = list(map(int, cluster_ion_list[j].split('|')))
            cos_list = list()
            match_list = list()
            id_list = list()
            prod_list = list()
            for k in ion_list:
                score, n_matches = modified_cosine.pair(i_spectrum, pos_mgf[k])
                id_list.append(int(pos_mgf[k].get(idx_column)))
                cos_list.append(score)
                match_list.append(n_matches)
                prod_list.append(score * n_matches)
            tmp_prod_list.append(max(prod_list))
            tmp_cosine_list.append(cos_list[prod_list.index(max(prod_list))])
            tmp_match_list.append(match_list[prod_list.index(max(prod_list))])
            tmp_id_list.append(id_list[prod_list.index(max(prod_list))])
            tmp_cluster_list.append(j)
        tmp_table = pd.DataFrame(list(zip(tmp_id_list, tmp_cluster_list, tmp_cosine_list, tmp_prod_list, tmp_match_list)),
                                 columns = [f'ion_{idx_column}', 'cluster_id', 'cosine', 'prod', 'matches'])
        tmp_table = tmp_table.loc[tmp_table['prod'].idxmax()]
        
        counter_feature_id = tmp_table[f'ion_{idx_column}']
        counter_idx = pos_node_table.index[pos_node_table[idx_column] == counter_feature_id][0]
        full_node_1.append(i)
        full_node_2.append(counter_idx)
        full_cluster_ids.append(tmp_table['cluster_id'])
        full_cosine.append(round(tmp_table['cosine'], 2))
        full_matches.append(tmp_table['matches'])
        
    cosine_table = pd.DataFrame()
    cosine_table['node_1'] = full_node_1
    cosine_table['node_2'] = full_node_2
    cosine_table['cluster_id'] = full_cluster_ids
    cosine_table['cosine'] = full_cosine
    cosine_table['matches'] = full_matches
    cosine_table = cosine_table[cosine_table['cosine'] >= cosiner_threshold]
    cosine_table = cosine_table[cosine_table['matches'] >= matched_peaks]
    
    edge_table['cosine_score'] = [0.0]*len(edge_table)
    for i in tqdm(cosine_table.index):
        cosined_ion = cosine_table.loc[i, "node_1"]
        linked_ion = cosine_table.loc[i, "node_2"]
        old_edge = edge_table.index[edge_table['node_1'] == cosined_ion][0]
        if (edge_table.loc[old_edge, "status"] == "neg_self_edge") or (edge_table.loc[old_edge, "status"] == "pos_self_edge"):
            edge_table.drop(old_edge, inplace = True)
        node_table.loc[cosined_ion, "status"] = node_table.loc[cosined_ion, "status"][:4] + "cossingleton"
        node_table.loc[cosined_ion, "status_universal"] = "cossingleton"
        node_table.loc[cosined_ion, "cluster_id"] = cosine_table.loc[i, "cluster_id"]
        rt_gap = round(abs(node_table.loc[cosined_ion, rt_field] - node_table.loc[linked_ion, rt_field]),3)
        mz_gap = round(abs(node_table.loc[cosined_ion, mz_field] - node_table.loc[linked_ion, mz_field]),4)
        cosine_score = cosine_table.loc[i, "cosine"]
        ion_mode = node_table.loc[cosined_ion, "ion_mode"]
        n_matches = cosine_table.loc[i, "matches"]
        new_edge = max(edge_table.index) + 1
        edge_table.loc[new_edge] = [None]*len(edge_table.columns)
        edge_table.loc[new_edge, ["node_1", "node_2", "matched_peaks", "total_peaks",
                                  "matching_score", "rt_gap", "mz_gap", "status",
                                  "status_universal", "All_annotations", "ion_mode",
                                  "cosine_score"]] = [linked_ion, cosined_ion,
                    n_matches, 0, 0, rt_gap, mz_gap, ion_mode.lower() + "_singcos_edge",
                    "singcos_edge", cosine_score, ion_mode, cosine_score]
    edge_table.reset_index(drop = True, inplace = True)
    
    neutral_idx = list(node_table.index[node_table['status'].str.contains('neutral')])
    for neutral in tqdm(neutral_idx):
        tmp_edge_table = edge_table[edge_table['node_1'] == neutral]
        pos_ions = list(tmp_edge_table['node_2'][tmp_edge_table['ion_mode'] == "POS"])
        neg_ions = list(tmp_edge_table['node_2'][tmp_edge_table['ion_mode'] == "NEG"])
        tmp_edges = []
        while len(pos_ions) > 1 :
            ion_1 = pos_ions[0]
            ion_1_mgf_idx = int(node_table.loc[ion_1, "spec_id"])
            ion_1_spectrum = pos_mgf[ion_1_mgf_idx]
            ion_1_rt = node_table.loc[ion_1, rt_field]
            ion_1_mz = node_table.loc[ion_1, mz_field]
            pos_ions.remove(ion_1)
            for ion_2 in pos_ions:
                ion_2_mgf_idx = int(node_table.loc[ion_2, "spec_id"])
                ion_2_spectrum = pos_mgf[ion_2_mgf_idx]
                ion_2_rt = node_table.loc[ion_2, rt_field]
                ion_2_mz = node_table.loc[ion_2, mz_field]
                rt_gap = round(ion_1_rt - ion_2_rt, 3)
                mz_gap = round(ion_1_mz - ion_2_mz, 4)
                score, n_matches = modified_cosine.pair(ion_1_spectrum, ion_2_spectrum)
                tmp_edges.append((ion_1, ion_2, round(score, 2), rt_gap, mz_gap, "pos"))
        while len(neg_ions) > 1 :
            ion_1 = neg_ions[0]
            ion_1_mgf_idx = int(node_table.loc[ion_1, "spec_id"])
            ion_1_spectrum = neg_mgf[ion_1_mgf_idx]
            ion_1_rt = node_table.loc[ion_1, rt_field]
            ion_1_mz = node_table.loc[ion_1, mz_field]
            neg_ions.remove(ion_1)
            for ion_2 in neg_ions:
                ion_2_mgf_idx = int(node_table.loc[ion_2, "spec_id"])
                ion_2_spectrum = neg_mgf[ion_2_mgf_idx]
                ion_2_rt = node_table.loc[ion_2, rt_field]
                ion_2_mz = node_table.loc[ion_2, mz_field]
                rt_gap = round(ion_1_rt - ion_2_rt, 3)
                mz_gap = round(ion_1_mz - ion_2_mz, 4)
                score, n_matches = modified_cosine.pair(ion_1_spectrum, ion_2_spectrum)
                tmp_edges.append((ion_1, ion_2, round(score, 2), rt_gap, mz_gap, "neg"))
        for edge in tmp_edges:
            node_1 = edge[0]
            node_2 = edge[1]
            cos = edge[2]
            if cos < cosine_threshold : continue
            rt_gap = edge[3]
            mz_gap = edge[4]
            tmp_ion_mode = edge[5]
            new_idx = max(edge_table.index) + 1
            edge_table.loc[new_idx] = [None]*len(edge_table.columns)
            edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks",
                                     "total_peaks", "matching_score", "rt_gap",
                                     "mz_gap", "status", "status_universal",
                                     "All_annotations", "ion_mode",
                                     "cosine_score"]] = [node_1, node_2, 0, 0, 0,
                                     rt_gap, mz_gap, tmp_ion_mode + "_cos_edge",
                                     "cos_edge", cos, tmp_ion_mode.upper(), cos]
    
    # Produce cosine clusters between singletons and non-molecular clustered precursors/fragments
    non_molecular_clusters = list(set(node_table['cluster_id'].unique()) - set(cluster_list))
    non_molecular_clusters.sort()
    unclustered_ions = [i for i in node_table.index if node_table.loc[i, "cluster_id"] in non_molecular_clusters]
    remains_ions_neg = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_singleton'])
    remains_ions_neg +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_precursor'])
    remains_ions_neg_mgf = [node_table.loc[i, "spec_id"] for i in remains_ions_neg] 
    remains_ions_neg_mgf = list(map(int, remains_ions_neg_mgf))

    remains_ions_pos = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_singleton'])
    remains_ions_pos +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_precursor'])
    remains_ions_pos_mgf = [node_table.loc[i, "spec_id"] for i in remains_ions_pos] 
    remains_ions_pos_mgf = list(map(int, remains_ions_pos_mgf))
    
    # Process NEG data:
    singleton_clusters = list()
    total_nodes = len(remains_ions_neg)
    while len(remains_ions_neg) > 0 :
        perc = round((1-(len(remains_ions_neg)/total_nodes))*100,1)
        sys.stdout.write("\rClustering remaining singletons (NEG) : {0}%".format(perc))
        sys.stdout.flush()
        ion_i = remains_ions_neg[0]
        ion_i_mgf = remains_ions_neg_mgf[0]
        remains_ions_neg.remove(ion_i)
        remains_ions_neg_mgf.remove(ion_i_mgf)
        for ion_j in remains_ions_neg:
            ion_j_mgf = remains_ions_neg_mgf[remains_ions_neg.index(ion_j)]
            score, n_matches = modified_cosine.pair(neg_mgf[ion_i_mgf], neg_mgf[ion_j_mgf])
            singleton_clusters.append((ion_i, ion_j, score, n_matches, "NEG"))

    # Process POS data:
    total_nodes = len(remains_ions_pos)
    while len(remains_ions_pos) > 0 :
        perc = round((1-(len(remains_ions_pos)/total_nodes))*100,1)
        sys.stdout.write("\rClustering remaining singletons (POS) : {0}%".format(perc))
        sys.stdout.flush()
        ion_i = remains_ions_pos[0]
        ion_i_mgf = remains_ions_pos_mgf[0]
        remains_ions_pos.remove(ion_i)
        remains_ions_pos_mgf.remove(ion_i_mgf)
        for ion_j in remains_ions_pos:
            ion_j_mgf = remains_ions_pos_mgf[remains_ions_pos.index(ion_j)]
            score, n_matches = modified_cosine.pair(pos_mgf[ion_i_mgf], pos_mgf[ion_j_mgf])
            singleton_clusters.append((ion_i, ion_j, score, n_matches, "POS"))    
    singleton_clusters = pd.DataFrame(singleton_clusters, columns = ['node_1', 'node_2', 'cos', 'matches', 'ion_mode'])
    singleton_clusters = singleton_clusters[singleton_clusters['cos'] >= cosiner_threshold]
    singleton_clusters = singleton_clusters[singleton_clusters['matches'] >= matched_peaks]
    
    # Define the new clusters
    node_pool = singleton_clusters['node_1'].tolist() + singleton_clusters['node_2'].tolist()
    node_pool = list(set(node_pool))
    node_pool.sort()
    cluster_list = []
    cluster_size_list = []
    total_nodes = len(node_pool)
    while len(node_pool) > 0:
        new_cluster = [node_pool[0]]
        cluster_size = 0
        perc = round((1-(len(node_pool)/total_nodes))*100,1)
        sys.stdout.write("\rDefining cosine singleton clusters : {0}%".format(perc))
        sys.stdout.flush()
        while cluster_size != len(new_cluster):
            cluster_size = len(new_cluster)
            tmp_idx = []
            for i in new_cluster:
                tmp_idx += list(singleton_clusters.index[singleton_clusters['node_1'] == i])
                tmp_idx += list(singleton_clusters.index[singleton_clusters['node_2'] == i])
            new_cluster += list(singleton_clusters.loc[tmp_idx, 'node_1'])
            new_cluster += list(singleton_clusters.loc[tmp_idx, 'node_2'])
            new_cluster = list(set(new_cluster))
        new_cluster.sort()
        node_pool = list(set(node_pool) - set(new_cluster))
        cluster_size_list.append(len(new_cluster))
        cluster_list.append('|'.join(list(map(str, new_cluster))))

    cluster_table= pd.DataFrame()
    cluster_table['cluster'] = cluster_list
    cluster_table['cluster_size'] = cluster_size_list
    cluster_table.sort_values('cluster_size', ascending = False, inplace = True)
    cluster_table.reset_index(drop = True, inplace = True)
    cluster_table.set_index(cluster_table.index + (int(node_table['cluster_id'].max()) + 1), drop = True, inplace = True)
    
    # Report the new data to the edge table:
    print('Reporting cosined singletons to edge table...')
    for i in tqdm(singleton_clusters.index):
        new_edge = edge_table.index.max() + 1
        node_1 = singleton_clusters.loc[i, "node_1"]
        node_2 = singleton_clusters.loc[i, "node_2"]
        del_edge_1 = edge_table[edge_table['node_1'] == node_1]
        del_edge_1 = del_edge_1.index[del_edge_1['node_2'] == node_1]
        del_edge_2 = edge_table[edge_table['node_1'] == node_2]
        del_edge_2 = del_edge_2.index[del_edge_2['node_2'] == node_2]
        if len(del_edge_1) > 0 :
            edge_table.drop(del_edge_1[0], inplace = True)
        if len(del_edge_2) > 0 :
            edge_table.drop(del_edge_2[0], inplace = True)
        tmp_cos = round(singleton_clusters.loc[i, "cos"],2)
        tmp_matches = singleton_clusters.loc[i, "matches"]
        tmp_mode = singleton_clusters.loc[i, "ion_mode"]
        rt_gap = abs(node_table.loc[node_1, rt_field] - node_table.loc[node_2, rt_field])
        mz_gap = abs(node_table.loc[node_1, mz_field] - node_table.loc[node_2, mz_field])
        edge_table.loc[new_edge] = [None]*len(edge_table.columns)
        edge_table.loc[new_edge, ["node_1", "node_2", "matched_peaks",
                                 "total_peaks", "matching_score", "rt_gap",
                                 "mz_gap", "status", "status_universal",
                                 "All_annotations", "ion_mode",
                                 "cosine_score"]] = [node_1, node_2, tmp_matches,
                                 0, 0.0, rt_gap, mz_gap, tmp_mode.lower() + '_singcos_edge',
                                 'cosine_edge', tmp_cos, tmp_mode, tmp_cos]

    # Report the new data to the node table:
    print('Reporting cosined singletons to node table...')
    for i in tqdm(cluster_table.index):
        node_list = list(map(int, cluster_table.loc[i, "cluster"].split('|')))
        tmp_mode = node_table.loc[node_list[0], 'ion_mode']
        for j in node_list:
            node_table.loc[j, "status"] = tmp_mode.lower() + "_cossingleton"
            node_table.loc[j, "status_universal"] = "cossingleton"
            node_table.loc[j, "cluster_id"] = i
    
    if purge_empty_spectra :
        node_table, edge_table = Spectrum_purge(ion_mode = "POS", node_table = node_table, edge_table = edge_table, mgf_file = pos_mgf)
        node_table, edge_table = Spectrum_purge(ion_mode = "NEG", node_table = node_table, edge_table = edge_table, mgf_file = neg_mgf)

        
    node_table[mz_field] = node_table[mz_field].round(4)
    node_table[rt_field] =  node_table[rt_field].round(3)
    edge_table['rt_gap'] = edge_table['rt_gap'].round(3)
    edge_table['mz_gap'] = edge_table['mz_gap'].round(4)
    edge_table['cosine_score'] = edge_table['cosine_score'].round(2)
    
    # Export data
    node_table.to_csv(out_path_full + 'node_table.csv', index_label = "Index")
    edge_table.to_csv(out_path_full + 'edge_table.csv', index_label = "Index")
    
    if params['c_export_samples'] : 
        dr_samplewise_export(neg_csv_file = mzmine_path_neg + neg_csv_file,
                          pos_csv_file = mzmine_path_pos + pos_csv_file,
                          out_path = out_path_samples,
                          edge_table = edge_table,
                          node_table = node_table,
                          params = params)
    return
