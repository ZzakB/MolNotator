def Cosiner(params : dict):

    def Column_correction(table):
        drop_col = [i for i in table.columns if "Unnamed" in i]
        table.drop(drop_col, axis = 1, inplace = True)
        return table
    def Samplewise_export(neg_csv_file, pos_csv_file, out_path, merged_edge_table, merged_node_table) : 
        print("Exporting sample-wise tables...")
        neg_csv = pd.read_csv(neg_csv_file, index_col ="row ID")
        pos_csv = pd.read_csv(pos_csv_file, index_col ="row ID")
        neg_csv = Column_correction(neg_csv)
        pos_csv = Column_correction(pos_csv)
        neg_csv.columns = neg_csv.columns.str.replace(".mzXML Peak area", "").str.replace('NEG_', '')
        pos_csv.columns = pos_csv.columns.str.replace(".mzXML Peak area", "").str.replace('POS_', '')
        neg_csv.drop(["row m/z", "row retention time"], axis = 1, inplace = True)
        pos_csv.drop(["row m/z", "row retention time"], axis = 1, inplace = True)
        samples = list(set(list(neg_csv.columns) + list(pos_csv.columns)))
        samples.sort()
        
        nodes_neg = node_table[node_table['ion_mode'] == "NEG"]
        nodes_neg = nodes_neg['feature_id'][nodes_neg['status_universal'] != "neutral"].astype(int).tolist()
        neg_csv = neg_csv.loc[nodes_neg]
        nodes_pos = node_table[node_table['ion_mode'] == "POS"]
        nodes_pos = nodes_pos['feature_id'][nodes_pos['status_universal'] != "neutral"].astype(int).tolist()
        pos_csv = pos_csv.loc[nodes_pos]
        
        for sample in tqdm(samples):
            #sample = samples[0]
            ion_ids_neg = neg_csv.index[neg_csv[sample] > 0.0]
            ion_ids_pos = pos_csv.index[pos_csv[sample] > 0.0]
            
            #convert feature_ids to the new indexes
            tmp_table = merged_node_table[merged_node_table['status'] != "neg_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "pos_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "mix_neutral"]
            tmp_table_pos = tmp_table[tmp_table['ion_mode'] == "POS"]
            tmp_table_neg = tmp_table[tmp_table['ion_mode'] == "NEG"]
            ion_idx_neg = pd.Series(tmp_table_neg.index, index = tmp_table_neg['feature_id'])
            ion_idx_neg = list(ion_idx_neg[ion_ids_neg])
            ion_idx_pos = pd.Series(tmp_table_pos.index, index = tmp_table_pos['feature_id'])
            ion_idx_pos = list(ion_idx_pos[ion_ids_pos])
            ion_idx_mix = ion_idx_neg + ion_idx_pos
        
            # Get sample neutrals
            neutral_edges = merged_edge_table.loc[merged_edge_table["Adnotation"].dropna().index]
            kept_edges = [i for i in neutral_edges.index if neutral_edges.loc[i, "node_2"] in ion_idx_mix]
    
        
            # Get ion edges
            ion_edges = merged_edge_table[merged_edge_table['status'] != "neg_add_edge"]
            ion_edges = ion_edges[ion_edges['status'] != "pos_add_edge"]
            for i in ion_edges.index:
                if ion_edges.loc[i, "node_1"] in ion_idx_mix:
                    if ion_edges.loc[i, "node_2"] in ion_idx_mix:
                        kept_edges.append(i)
            kept_edges.sort()
            sample_edges = merged_edge_table.loc[kept_edges]
            sample_edges.sort_values('node_1', inplace = True)
            sample_edges.reset_index(inplace = True, drop = True)
            
            kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
            kept_nodes.sort()
            sample_nodes = merged_node_table.loc[kept_nodes].copy()
            sample_nodes.drop(pd.Series(samples) + ".mzXML Peak area", axis = 1, inplace = True)
            sample_nodes[sample] = merged_node_table[sample + ".mzXML Peak area"]
            
            sample_nodes.to_csv(out_path + "MIX_" + sample + "_nodes.csv", index_label = "Index")
            sample_edges.to_csv(out_path + "MIX_" + sample + "_edges.csv", index_label = "Index")
        return
    
    def Spectrum_purge(ion_mode, mgf_file) :
        singletons = node_table.index[node_table['status'] == ion_mode.lower() + "_singleton"].tolist()
        empty_nodes = list()
        print(f"Purging empty {ion_mode} singletons...")
        for i in tqdm(singletons):
            mgf_idx = int(node_table.loc[i, "mgf_index"])
            if len(mgf_file[mgf_idx].peaks.mz) == 0 :
                node_table.drop
                empty_nodes.append(i)
        del_edges = [edge_table.index[edge_table['node_1'] == n][0] for n in empty_nodes]
        node_table.drop(empty_nodes, inplace = True)
        edge_table.drop(del_edges, inplace = True)
        return    
    
    import os
    import pandas as pd
    import sys
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.similarity import ModifiedCosine
    from tqdm import tqdm
    def Spectrum_processing(s):
        s = default_filters(s)
        return s
    
    # Load parameters
    mzmine_path_neg= params['neg_out_0']
    mzmine_path_pos= params['pos_out_0']
    neg_csv_file= params['neg_csv']
    pos_csv_file= params['pos_csv']
    neg_mgf_file= params['neg_mgf']
    pos_mgf_file= params['pos_mgf']
    in_path= params['mix_out_5_1']
    out_path_full= params['mix_out_6_1']
    out_path_samples= params['mix_out_6_2']
    purge_empty_spectra = params['c_purge_empty_spectra']
    
    mass_error= params['c_mass_error']
    cosine_threshold= params['c_lowcos_threshold']
    cosiner_threshold= params['c_hardcos_threshold']
    matched_peaks= params['c_matched_peaks']
    modified_cosine = ModifiedCosine(tolerance=mass_error)

    # Load files
    node_table = pd.read_csv(in_path + "MIX_nodes.csv", index_col = "Index")
    edge_table  = pd.read_csv(in_path + "MIX_edges.csv", index_col = "Index")
    neg_mgf = list(load_from_mgf(mzmine_path_neg + neg_mgf_file))
    pos_mgf = list(load_from_mgf(mzmine_path_pos + pos_mgf_file))
    neg_mgf = [Spectrum_processing(s) for s in neg_mgf]
    pos_mgf = [Spectrum_processing(s) for s in pos_mgf]
    
    
    # Make a Series with MGF indexes as data and feature IDs as indexes
    neg_mgf_data = pd.Series()
    for i in range(len(neg_mgf)):
        neg_mgf_data.loc[int(neg_mgf[i].get("feature_id"))] = i
    
    pos_mgf_data = pd.Series()
    for i in range(len(pos_mgf)):
        pos_mgf_data.loc[int(pos_mgf[i].get("feature_id"))] = i
    
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
    
    neg_node_table = node_table[node_table['ion_mode'] == "NEG"]
    cluster_ion_list = pd.Series(index = cluster_list_neg)
    for i in cluster_ion_list.index:
        tmp_rows = neg_node_table.index[neg_node_table['cluster_id'] == i]
        cluster_ion_list[i] = '|'.join(neg_node_table.loc[tmp_rows, 'mgf_index'].dropna().astype(int).astype(str))
    cluster_ion_list.sort_index(inplace = True)
    full_node_1 = []
    full_node_2 = []
    full_cluster_ids = []
    full_cosine = []
    full_matches = []
    for i in tqdm(remains_ions_neg):
        i_id = int(node_table.loc[i, "feature_id"])
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
                id_list.append(int(neg_mgf[k].get("feature_id")))
                cos_list.append(score)
                match_list.append(n_matches)
                prod_list.append(score * n_matches)
            tmp_prod_list.append(max(prod_list))
            tmp_cosine_list.append(cos_list[prod_list.index(max(prod_list))])
            tmp_match_list.append(match_list[prod_list.index(max(prod_list))])
            tmp_id_list.append(id_list[prod_list.index(max(prod_list))])
            tmp_cluster_list.append(j)
        tmp_table = pd.DataFrame(list(zip(tmp_id_list, tmp_cluster_list, tmp_cosine_list, tmp_prod_list, tmp_match_list)),
                                 columns = ['ion_feature_id', 'cluster_id', 'cosine', 'prod', 'matches'])
        tmp_table = tmp_table.loc[tmp_table['prod'].idxmax()]
        
        counter_feature_id = tmp_table['ion_feature_id']
        counter_idx = neg_node_table.index[neg_node_table['feature_id'] == counter_feature_id][0]
        full_node_1.append(i)
        full_node_2.append(counter_idx)
        full_cluster_ids.append(tmp_table['cluster_id'])
        full_cosine.append(round(tmp_table['cosine'], 2))
        full_matches.append(tmp_table['matches'])
        
    pos_node_table = node_table[node_table['ion_mode'] == "POS"]
    cluster_ion_list = pd.Series(index = cluster_list_pos)
    for i in cluster_ion_list.index:
        tmp_rows = pos_node_table.index[pos_node_table['cluster_id'] == i]
        cluster_ion_list[i] = '|'.join(pos_node_table.loc[tmp_rows, 'mgf_index'].dropna().astype(int).astype(str))
    cluster_ion_list.sort_index(inplace = True)
    for i in tqdm(remains_ions_pos):
        i_id = int(node_table.loc[i, "feature_id"])
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
                id_list.append(int(pos_mgf[k].get("feature_id")))
                cos_list.append(score)
                match_list.append(n_matches)
                prod_list.append(score * n_matches)
            tmp_prod_list.append(max(prod_list))
            tmp_cosine_list.append(cos_list[prod_list.index(max(prod_list))])
            tmp_match_list.append(match_list[prod_list.index(max(prod_list))])
            tmp_id_list.append(id_list[prod_list.index(max(prod_list))])
            tmp_cluster_list.append(j)
        tmp_table = pd.DataFrame(list(zip(tmp_id_list, tmp_cluster_list, tmp_cosine_list, tmp_prod_list, tmp_match_list)),
                                 columns = ['ion_feature_id', 'cluster_id', 'cosine', 'prod', 'matches'])
        tmp_table = tmp_table.loc[tmp_table['prod'].idxmax()]
        
        counter_feature_id = tmp_table['ion_feature_id']
        counter_idx = pos_node_table.index[pos_node_table['feature_id'] == counter_feature_id][0]
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
        rt_gap = round(abs(node_table.loc[cosined_ion, "rt"] - node_table.loc[linked_ion, "rt"]),3)
        mz_gap = round(abs(node_table.loc[cosined_ion, "mz"] - node_table.loc[linked_ion, "mz"]),4)
        cosine_score = cosine_table.loc[i, "cosine"]
        ion_mode = node_table.loc[cosined_ion, "ion_mode"]
        n_matches = cosine_table.loc[i, "matches"]
        new_edge = max(edge_table.index) + 1
        edge_table.loc[new_edge] = [linked_ion, cosined_ion, n_matches, 0, 0, rt_gap, mz_gap,
                      ion_mode.lower() + "_singcos_edge", "singcos_edge", None, None, cosine_score, ion_mode, cosine_score]
    edge_table.reset_index(drop = True, inplace = True)
    
    neutral_idx = list(node_table.index[node_table['status'].str.contains('neutral')])
    for neutral in tqdm(neutral_idx):
        tmp_edge_table = edge_table[edge_table['node_1'] == neutral]
        pos_ions = list(tmp_edge_table['node_2'][tmp_edge_table['ion_mode'] == "POS"])
        neg_ions = list(tmp_edge_table['node_2'][tmp_edge_table['ion_mode'] == "NEG"])
        tmp_edges = []
        while len(pos_ions) > 1 :
            ion_1 = pos_ions[0]
            ion_1_mgf_idx = int(node_table.loc[ion_1, "mgf_index"])
            ion_1_spectrum = pos_mgf[ion_1_mgf_idx]
            ion_1_rt = node_table.loc[ion_1, "rt"]
            ion_1_mz = node_table.loc[ion_1, "mz"]
            pos_ions.remove(ion_1)
            for ion_2 in pos_ions:
                ion_2_mgf_idx = int(node_table.loc[ion_2, "mgf_index"])
                ion_2_spectrum = pos_mgf[ion_2_mgf_idx]
                ion_2_rt = node_table.loc[ion_2, "rt"]
                ion_2_mz = node_table.loc[ion_2, "mz"]
                rt_gap = round(ion_1_rt - ion_2_rt, 3)
                mz_gap = round(ion_1_mz - ion_2_mz, 4)
                score, n_matches = modified_cosine.pair(ion_1_spectrum, ion_2_spectrum)
                tmp_edges.append((ion_1, ion_2, round(score, 2), rt_gap, mz_gap, "pos"))
        while len(neg_ions) > 1 :
            ion_1 = neg_ions[0]
            ion_1_mgf_idx = int(node_table.loc[ion_1, "mgf_index"])
            ion_1_spectrum = neg_mgf[ion_1_mgf_idx]
            ion_1_rt = node_table.loc[ion_1, "rt"]
            ion_1_mz = node_table.loc[ion_1, "mz"]
            neg_ions.remove(ion_1)
            for ion_2 in neg_ions:
                ion_2_mgf_idx = int(node_table.loc[ion_2, "mgf_index"])
                ion_2_spectrum = neg_mgf[ion_2_mgf_idx]
                ion_2_rt = node_table.loc[ion_2, "rt"]
                ion_2_mz = node_table.loc[ion_2, "mz"]
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
            edge_table.loc[new_idx] = [node_1, node_2, 0, 0, 0, rt_gap, mz_gap,
                          tmp_ion_mode + "_cos_edge", "cos_edge", None, None, cos, tmp_ion_mode.upper(),
                          cos]
    
    # Produce cosine clusters between singletons and non-molecular clustered precursors/fragments
    non_molecular_clusters = list(set(node_table['cluster_id'].unique()) - set(cluster_list))
    non_molecular_clusters.sort()
    unclustered_ions = [i for i in node_table.index if node_table.loc[i, "cluster_id"] in non_molecular_clusters]
    remains_ions_neg = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_singleton'])
    remains_ions_neg +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'neg_precursor'])
    remains_ions_neg_mgf = [node_table.loc[i, "mgf_index"] for i in remains_ions_neg] 
    remains_ions_neg_mgf = list(map(int, remains_ions_neg_mgf))

    remains_ions_pos = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_singleton'])
    remains_ions_pos +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status'] == 'pos_precursor'])
    remains_ions_pos_mgf = [node_table.loc[i, "mgf_index"] for i in remains_ions_pos] 
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
        rt_gap = abs(node_table.loc[node_1, 'rt'] - node_table.loc[node_2, 'rt'])
        mz_gap = abs(node_table.loc[node_1, 'mz'] - node_table.loc[node_2, 'mz'])
        edge_table.loc[new_edge] = [node_1, node_2, tmp_matches, 0, 0.0, rt_gap,
                      mz_gap, tmp_mode.lower() + '_singcos_edge', 'cosine_edge',
                      None, None, tmp_cos, tmp_mode, tmp_cos]

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
        Spectrum_purge(ion_mode = "POS", mgf_file = pos_mgf)
        Spectrum_purge(ion_mode = "NEG", mgf_file = neg_mgf)
        
    node_table['mz'] = node_table['mz'].round(4)
    node_table['rt'] =  node_table['rt'].round(3)
    edge_table['rt_gap'] = edge_table['rt_gap'].round(3)
    edge_table['mz_gap'] = edge_table['mz_gap'].round(4)
    edge_table['cosine_score'] = edge_table['cosine_score'].round(2)
    
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
    if not os.path.isdir(out_path_samples) :
        os.mkdir(out_path_samples)
    
    node_table.to_csv(out_path_full + 'MIX_nodes.csv', index_label = "Index")
    edge_table.to_csv(out_path_full + 'MIX_edges.csv', index_label = "Index")
    
    if params['c_export_samples'] : 
        Samplewise_export(neg_csv_file = mzmine_path_neg + neg_csv_file,
                          pos_csv_file = mzmine_path_pos + pos_csv_file,
                          out_path = out_path_samples,
                          merged_edge_table = edge_table,
                          merged_node_table = node_table)
    return