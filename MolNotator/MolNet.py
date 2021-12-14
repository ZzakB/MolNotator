def MolNet(params : dict):

    def Column_correction(table):
        drop_col = [i for i in table.columns if "Unnamed" in i]
        table.drop(drop_col, axis = 1, inplace = True)
        return table

    def Samplewise_export(neg_csv_file, pos_csv_file, node_table_unaltered, edge_table_unaltered, out_path, merged_edge_table, merged_node_table) : 
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
        node_table_unaltered = pd.read_csv(node_table_unaltered, index_col = "Index")
        edge_table_unaltered = pd.read_csv(edge_table_unaltered, index_col = "Index")
        
        nodes_neg = node_table_unaltered[node_table_unaltered['ion_mode'] == "NEG"]
        nodes_neg = nodes_neg['feature_id'][nodes_neg['status_universal'] != "neutral"].astype(int).tolist()
        neg_csv = neg_csv.loc[nodes_neg]
        nodes_pos = node_table_unaltered[node_table_unaltered['ion_mode'] == "POS"]
        nodes_pos = nodes_pos['feature_id'][nodes_pos['status_universal'] != "neutral"].astype(int).tolist()
        pos_csv = pos_csv.loc[nodes_pos]
        
        edge_subtable = edge_table_unaltered.copy()
        edge_subtable = edge_subtable[edge_subtable['status_universal'] == 'add_edge']
        edge_subtable['Index'] = edge_subtable.index
        edge_subtable.set_index('node_2', drop = False, inplace = True)
        edge_subtable.drop('node_2', axis = 1, inplace = True)
        
        
        for sample in tqdm(samples):
            #sample = samples[0]
            ion_ids_neg = neg_csv.index[neg_csv[sample] > 0.0]
            ion_ids_pos = pos_csv.index[pos_csv[sample] > 0.0]
            
            #convert feature_ids to the new indexes
            tmp_node_table = node_table_unaltered.copy()
            tmp_node_table['Index'] = tmp_node_table.index
            tmp_node_table = tmp_node_table[tmp_node_table['status_universal'] != "neutral"]
            tmp_node_table = tmp_node_table.set_index(tmp_node_table['feature_id'].astype(int), drop = False)
            tmp_node_table_neg = tmp_node_table[tmp_node_table['ion_mode'] == "NEG"]
            tmp_node_table_pos = tmp_node_table[tmp_node_table['ion_mode'] == "POS"]   
            
            tmp_node_table_neg = tmp_node_table_neg.loc[ion_ids_neg]
            tmp_node_table_neg.set_index(tmp_node_table_neg['Index'], drop = False, inplace = True)
            tmp_node_table_neg.drop('Index', axis = 1, inplace = True)
            tmp_node_table_pos = tmp_node_table_pos.loc[ion_ids_pos]
            tmp_node_table_pos.set_index(tmp_node_table_pos['Index'], drop = False, inplace = True)
            tmp_node_table_pos.drop('Index', axis = 1, inplace = True)
            
            neg_adducts = set(tmp_node_table_neg.index)
            neg_adducts = list(neg_adducts.intersection(edge_subtable.index))
            tmp_neutrals_neg = edge_subtable.loc[neg_adducts]["node_1"].tolist()

            pos_adducts = set(tmp_node_table_pos.index)
            pos_adducts = list(pos_adducts.intersection(edge_subtable.index))
            tmp_neutrals_pos = edge_subtable.loc[pos_adducts]["node_1"].tolist()
            
            kept_neutrals = list(set(tmp_neutrals_neg + tmp_neutrals_pos))
            kept_edges = list()
            for i in merged_edge_table.index:
                if (merged_edge_table.loc[i, "node_1"] in kept_neutrals) and (merged_edge_table.loc[i, "node_2"] in kept_neutrals):
                    kept_edges.append(i)
            
            kept_edges.sort()
            sample_node_table = merged_node_table.loc[kept_neutrals].copy()
            sample_edge_table = merged_edge_table.loc[kept_edges].copy()

            sample_edge_table.sort_values('node_1', inplace = True)
            sample_edge_table.reset_index(inplace = True, drop = True)
            
            sample_node_table.drop(pd.Series(samples) + ".mzXML Peak area", axis = 1, inplace = True)
            sample_node_table[sample] = merged_node_table[sample + ".mzXML Peak area"]
            
            sample_node_table.to_csv(out_path + "MIX_" + sample + "_nodes.csv", index_label = "Index")
            sample_edge_table.to_csv(out_path + "MIX_" + sample + "_edges.csv", index_label = "Index")
        
        return
    
    def Spectrum_processing(s):
        s = default_filters(s)
        return s

    import os
    import pandas as pd
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.similarity import ModifiedCosine
    import sys
    from tqdm import tqdm

    # Load parameters
    neg_csv_file= params['neg_csv']
    pos_csv_file= params['pos_csv']
    mgf_path_neg= params['neg_out_0']
    mgf_path_pos= params['pos_out_0']
    mgf_file_neg= params['neg_mgf']
    mgf_file_pos= params['pos_mgf']
    in_path= params['mix_out_6_1']
    out_path_full= params['mix_out_7_1']
    out_path_samples= params['mix_out_7_2']
    mass_error= params['mn_mass_error']
    cosine_threshold= params['mn_cosine_threshold']
    matched_peaks= params['mn_matched_peaks']
    
    node_table = pd.read_csv(in_path + 'MIX_nodes.csv', index_col = "Index")
    edge_table = pd.read_csv(in_path + 'MIX_edges.csv', index_col = "Index")
    mgf_pos = list(load_from_mgf(mgf_path_pos + mgf_file_pos))
    mgf_neg = list(load_from_mgf(mgf_path_neg + mgf_file_neg))
    
    mgf_pos = [Spectrum_processing(s) for s in mgf_pos]
    mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
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
    node_table['mz'] = node_table['mz'].round(4)
    edge_table = kept_edges_1.append(kept_edges_2, ignore_index = True)
    
    node_table.to_csv(out_path_full + "MIX_nodes_simple.csv", index_label = "Index")
    edge_table.to_csv(out_path_full + "MIX_edges_simple.csv", index_label = "Index")
    
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
        adducts_1['mgf_index'] = list(node_table.loc[adducts_1.index, "mgf_index"].astype(int))
        adducts_1['ion_mode'] = list(node_table.loc[adducts_1.index, "ion_mode"])
        adducts_1['adduct'] = list(node_table.loc[adducts_1.index, "Adnotation"])
        neutral_nodes.remove(neutral_1)
        for neutral_2 in neutral_nodes:
            #neutral_2 = neutral_nodes[1]
            tmp_scores = list()
            adducts_2 = pd.DataFrame(index = list(edge_table['node_2'][edge_table['node_1'] == neutral_2]))
            adducts_2['mgf_index'] = list(node_table.loc[adducts_2.index, "mgf_index"].astype(int))
            adducts_2['ion_mode'] = list(node_table.loc[adducts_2.index, "ion_mode"])
            adducts_2['adduct'] = list(node_table.loc[adducts_2.index, "Adnotation"])
            for mode in adducts_1['ion_mode'].unique():
                #mode = adducts_1['ion_mode'].unique()[0]
                tmp_table_1 = adducts_1[adducts_1['ion_mode'] == mode]
                tmp_table_2 = adducts_2[adducts_2['ion_mode'] == mode]
                if mode == "POS" : mgf_list = mgf_pos
                else: mgf_list = mgf_neg
                for ion_1 in tmp_table_1.index:
                    spectrum_1 = mgf_list[tmp_table_1.loc[ion_1, "mgf_index"]]
                    adduct_1 = tmp_table_1.loc[ion_1, 'adduct']
                    for ion_2 in tmp_table_2.index:
                        adduct_2 = tmp_table_2.loc[ion_2, 'adduct']
                        if adduct_1 != adduct_2 : continue
                        spectrum_2 = mgf_list[tmp_table_2.loc[ion_2, "mgf_index"]]
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
        rt_gap = abs(round(node_table.loc[node_1, 'rt'] - node_table.loc[node_2, 'rt'], 3))
        mz_gap = abs(round(node_table.loc[node_1, 'mz'] - node_table.loc[node_2, 'mz'], 4))
        edge_table.loc[new_idx] = [node_1, node_2, matches, 0, 0.0, rt_gap, mz_gap,
                      edge_mode + "_cosine_neutral_edge", "cosine_neutral_edge", None, None, cos, edge_mode, cos]
    
    # Delete undesired edges
    edge_table = edge_table[edge_table['status'].str.contains('cosine_neutral_edge')]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    
    # Add singletons:
    singletons = list(set(neutral_nodes_2) - set(kept_nodes))
    for neutral in tqdm(singletons):
        new_idx = max(edge_table.index) + 1 
        ion_mode = node_table.loc[neutral, "ion_mode"]
        edge_table.loc[new_idx] = [neutral, neutral, 0, 0, 0.0, 0.0, 0.0,
                      ion_mode + "_self_edge", "self_edge", None, None, 0.0, ion_mode, 0.0]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    node_table = node_table.loc[kept_nodes]
    
    node_table['mz'] = node_table['mz'].round(4)
    
    node_table.to_csv(out_path_full + "MIX_nodes_families.csv", index_label = "Index")
    edge_table.to_csv(out_path_full + "MIX_edges_families.csv", index_label = "Index")
    
    if params['c_export_samples'] : 
        Samplewise_export(neg_csv_file = mgf_path_neg + neg_csv_file,
                          pos_csv_file = mgf_path_pos + pos_csv_file,
                          node_table_unaltered = in_path + 'MIX_nodes.csv',
                          edge_table_unaltered = in_path + 'MIX_edges.csv',
                          out_path = out_path_samples,
                          merged_edge_table = edge_table,
                          merged_node_table = node_table)

    return
