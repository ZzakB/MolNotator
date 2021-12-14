def MGF_updater(params: dict):
    
    def Spectrum_processing(s):
        s = default_filters(s)
        return s
    import os 
    import pandas as pd
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf
    from matchms.filtering import default_filters
    from tqdm import tqdm
    
    if params['mu_skip']: return
    mgf_file_neg = params['neg_mgf']
    mgf_file_pos = params['pos_mgf']
    in_folder_mgf_neg = params['neg_out_0']
    in_folder_mgf_pos = params['pos_out_0']
    in_folder_csv = params['mix_out_4_1']
    sirius_folder_neg = params['mu_sirius_folder_neg']
    sirius_folder_pos = params['mu_sirius_folder_pos']
    out_mgf_folder = params['mix_mgf_out_4']

    # Load node and edge tables
    print('Loading node and edge tables')
    node_table = pd.read_csv(in_folder_csv + "MIX_nodes.csv", index_col = "Index")
    edge_table = pd.read_csv(in_folder_csv + "MIX_edges.csv", index_col = "Index")

    # Load MGF files
    print('Loading NEG MGF file...')
    mgf_neg = list(load_from_mgf(in_folder_mgf_neg + mgf_file_neg))
    mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
    print('Loading POS MGF file...')
    mgf_pos = list(load_from_mgf(in_folder_mgf_pos + mgf_file_pos))
    mgf_pos = [Spectrum_processing(s) for s in mgf_pos]

    # Get the list of SIRIUS processed data:
    sirius_data_neg = os.listdir(sirius_folder_neg)
    sirius_data_neg.remove('.format')
    sirius_data_neg = pd.DataFrame(sirius_data_neg, columns = ['folder'])
    sirius_data_neg['feature_id'] = sirius_data_neg['folder'].str.split('_').str[-1].astype(int)

    sirius_data_pos = os.listdir(sirius_folder_pos)
    sirius_data_pos.remove('.format')
    sirius_data_pos = pd.DataFrame(sirius_data_pos, columns = ['folder'])
    sirius_data_pos['feature_id'] = sirius_data_pos['folder'].str.split('_').str[-1].astype(int)

    # Acquire the SIRIUS data to the node table, first for neutrals and adducts:
    neutral_idx = node_table.index[node_table['status_universal'] == "neutral"].tolist()
    node_table.insert(node_table.columns.get_loc("TIC"), "formula", [None]*len(node_table))
    print('Establishing consensus chemical formulas for neutral and adducts...')
    for i in tqdm(neutral_idx):
        tmp_edges = edge_table[edge_table['node_1'] == i]
        neg_ions = tmp_edges["node_2"][tmp_edges['ion_mode'] == "NEG"].tolist()
        pos_ions = tmp_edges["node_2"][tmp_edges['ion_mode'] == "POS"].tolist()   
        formula_table = list()
        
        # Get SIRIUS data for the NEG ions
        for j in neg_ions:
            j_adduct = node_table.loc[j, "Adnotation"]
            j_f_id = int(node_table.loc[j, "feature_id"])
            folder = sirius_data_neg['folder'][sirius_data_neg['feature_id'] == j_f_id]
            if len(folder) == 0 : continue
            folder = folder.iloc[0]
            tmp_path = sirius_folder_neg + folder + "/scores/"
            if not os.path.isdir(tmp_path): continue
            tmp_files = os.listdir(tmp_path)
            for f in tmp_files:
                tmp_form = f.split('_')[0]
                tmp_adduct = f.split('_')[1].replace('.info', '')
                if tmp_adduct != j_adduct : continue
                tmp_info = pd.read_csv(tmp_path + f, sep = '\t', names = ['name', 'score'], index_col = "name")
                tmp_score = tmp_info.loc['sirius.scores.SiriusScore', "score"]
                formula_table.append((tmp_form, tmp_score, tmp_adduct))

        # Get SIRIUS data for the POS ions
        for j in pos_ions:
            j_adduct = node_table.loc[j, "Adnotation"]
            j_f_id = int(node_table.loc[j, "feature_id"])
            folder = sirius_data_pos['folder'][sirius_data_pos['feature_id'] == j_f_id]
            if len(folder) == 0 : continue
            folder = folder.iloc[0]
            tmp_path = sirius_folder_pos + folder + "/scores/"
            if not os.path.isdir(tmp_path): continue
            tmp_files = os.listdir(tmp_path)
            for f in tmp_files:
                tmp_form = f.split('_')[0]
                tmp_adduct = f.split('_')[1].replace('.info', '')
                if tmp_adduct != j_adduct : continue
                tmp_info = pd.read_csv(tmp_path + f, sep = '\t', names = ['name', 'score'], index_col = "name")
                tmp_score = tmp_info.loc['sirius.scores.SiriusScore', "score"]
                formula_table.append((tmp_form, tmp_score, tmp_adduct))
        
        # Skip if no data:
        if len(formula_table) == 0 : continue
        
        # Select the highest scoring formula:
        formula_table = pd.DataFrame(formula_table, columns = ['formula', 'score', 'source'])
        unique_formula = formula_table['formula'].unique().tolist()
        tmp_scores = list()
        for formula in unique_formula:
            tmp_scores.append(formula_table['score'][formula_table['formula'] == formula].sum())
        selected_formula = unique_formula[tmp_scores.index(max(tmp_scores))]

        # Propagate the formula:
        node_table.loc[i, "formula"] = selected_formula
        for j in neg_ions : node_table.loc[j, "formula"] = selected_formula
        for j in pos_ions : node_table.loc[j, "formula"] = selected_formula
    
    # Acquire SIRIUS data for non-molecular ions:
    ions = node_table.index[node_table['status_universal'] == "singleton"].tolist()
    ions += node_table.index[node_table['status_universal'] == "fragment"].tolist()
    ions += node_table.index[node_table['status_universal'] == "precursor"].tolist()
    ions_neg = node_table.loc[ions].index[node_table.loc[ions, "ion_mode"] == "NEG"].tolist()
    ions_pos = node_table.loc[ions].index[node_table.loc[ions, "ion_mode"] == "POS"].tolist()
    # Process neg ions:
    print('Establishing formulas for neg singletons...')
    for i in tqdm(ions_neg):
        formula_list = list()
        score_list = list()
        adduct_list = list()
        f_id = int(node_table.loc[i, "feature_id"])
        folder = sirius_data_neg['folder'][sirius_data_neg['feature_id'] == f_id]
        if len(folder) == 0 : continue
        folder = folder.iloc[0]
        tmp_path = sirius_folder_neg + folder + "/scores/"
        if not os.path.isdir(tmp_path): continue
        tmp_files = os.listdir(tmp_path)
        for f in tmp_files:
            tmp_form = f.split('_')[0]
            tmp_adduct = f.split('_')[1].replace('.info', '')
            tmp_info = pd.read_csv(tmp_path + f, sep = '\t', names = ['name', 'score'], index_col = "name")
            tmp_score = tmp_info.loc['sirius.scores.SiriusScore', "score"]
            formula_list.append(tmp_form)
            score_list.append(tmp_score)
            adduct_list.append(tmp_adduct)
        best_score = score_list.index(max(score_list))
        selected_adduct = adduct_list[best_score]
        selected_formula = formula_list[best_score]
        node_table.loc[i, "formula"] = selected_formula
        node_table.loc[i, "Adnotation"] = selected_adduct

    # Process pos ions:
    print('Establishing formulas for pos singletons...')
    for i in tqdm(ions_pos):
        formula_list = list()
        score_list = list()
        adduct_list = list()
        f_id = int(node_table.loc[i, "feature_id"])
        folder = sirius_data_pos['folder'][sirius_data_pos['feature_id'] == f_id]
        if len(folder) == 0 : continue
        folder = folder.iloc[0]
        tmp_path = sirius_folder_pos + folder + "/scores/"
        if not os.path.isdir(tmp_path): continue
        tmp_files = os.listdir(tmp_path)
        for f in tmp_files:
            tmp_form = f.split('_')[0]
            tmp_adduct = f.split('_')[1].replace('.info', '')
            tmp_info = pd.read_csv(tmp_path + f, sep = '\t', names = ['name', 'score'], index_col = "name")
            tmp_score = tmp_info.loc['sirius.scores.SiriusScore', "score"]
            formula_list.append(tmp_form)
            score_list.append(tmp_score)
            adduct_list.append(tmp_adduct)
        best_score = score_list.index(max(score_list))
        selected_adduct = adduct_list[best_score]
        selected_formula = formula_list[best_score]
        node_table.loc[i, "formula"] = selected_formula
        node_table.loc[i, "Adnotation"] = selected_adduct

    # Update the MGF files with MolNotator data:
    node_table_neg = node_table[node_table['ion_mode'] == "NEG"]
    new_neg_mgf = list()
    print('Updating the NEG MGF file...')
    for s in tqdm(mgf_neg):
        tmp_s = s.clone()
        tmp_dict = tmp_s.metadata
        f_id = int(tmp_s.get("feature_id"))
        idx = node_table_neg.index[node_table_neg['feature_id'] == f_id][0]
        adduct = node_table_neg.loc[idx, 'Adnotation']
        if type(adduct) == float : adduct = None
        formula = node_table_neg.loc[idx, 'formula']
        status = node_table_neg.loc[idx, 'status_universal']
        tmp_dict['adduct'] = adduct
        tmp_dict['formula'] = formula
        tmp_dict['status'] = status
        tmp_dict['ionmode'] = "negative"
        tmp_s.metadata = tmp_dict
        new_neg_mgf.append(tmp_s)

    node_table_pos = node_table[node_table['ion_mode'] == "POS"]
    new_pos_mgf = list()
    print('Updating the POS MGF file...')
    for s in tqdm(mgf_pos):
        tmp_s = s.clone()
        tmp_dict = tmp_s.metadata
        f_id = int(tmp_s.get("feature_id"))
        idx = node_table_pos.index[node_table_pos['feature_id'] == f_id][0]
        adduct = node_table_pos.loc[idx, 'Adnotation']
        if type(adduct) == float : adduct = None
        formula = node_table_pos.loc[idx, 'formula']
        status = node_table_pos.loc[idx, 'status_universal']
        tmp_dict['adduct'] = adduct
        tmp_dict['formula'] = formula
        tmp_dict['status'] = status
        tmp_dict['ionmode'] = "positive"
        tmp_s.metadata = tmp_dict
        new_pos_mgf.append(tmp_s)        

    # Export the updated MGFs:
    if not os.path.isdir(out_mgf_folder):
        os.mkdir(out_mgf_folder)
    print("Exporting the updated NEG MGF...")
    save_as_mgf(new_neg_mgf, out_mgf_folder + mgf_file_neg)
    print("Exporting the updated POS MGF...")
    save_as_mgf(new_pos_mgf, out_mgf_folder + mgf_file_pos)
    
    # Export to the mode merger folder the updated node table:
    print("Updating the node table in the mode merger folder...")
    node_table.to_csv(in_folder_csv + "MIX_nodes.csv", index_label = "Index")
    return
