def Dereplicator(params : dict, db_params : dict):

    def Filter_choices_outer(db_params):
        def Three_filters(table):
            return ['|'.join(table[db_params['db_rt_field']]), '|'.join(table[db_params['db_adduct_field']]), '|'.join(table[db_params['db_formula_field_database']])]
        def RT_adducts(table):
            return ['|'.join(table[db_params['db_rt_field']]), '|'.join(table[db_params['db_adduct_field']])]
        def RT_formula(table):
            return ['|'.join(table[db_params['db_rt_field']]), '|'.join(table[db_params['db_formula_field_database']])]
        def Adducts_formula(table):
            return ['|'.join(table[db_params['db_adduct_field']]), '|'.join(table[db_params['db_formula_field_database']])]
        def RT_only(table):
            return ['|'.join(table[db_params['db_rt_field']])]
        def Adduct_only(table):
            return ['|'.join(table[db_params['db_adduct_field']])]
        def Formula_only(table):
            return ['|'.join(table[db_params['db_formula_field_database']])]
        def No_fields(table):
            return []
        if (db_params['db_rt_field'] != None) and (db_params['db_adduct_field'] != None) and (db_params['db_formula_field_database'] != None):
            return Three_filters, ['rt', 'adduct', 'formula']
        elif db_params['db_rt_field'] != None and (db_params['db_adduct_field'] != None):
            return RT_adducts, ['rt', 'adduct']
        elif db_params['db_rt_field'] != None and (db_params['db_formula_field_database'] != None):
            return RT_formula, ['rt', 'formula']
        elif db_params['db_adduct_field'] != None and (db_params['db_formula_field_database'] != None):
            return Adducts_formula, ['adduct', 'formula']
        elif db_params['db_rt_field'] != None:
            return RT_only, ['rt']
        elif db_params['db_adduct_field'] != None:
            return Adduct_only, ['adduct']
        elif db_params['db_formula_field_database'] != None:
            return Adduct_only, ['formula']
        else:
            return No_fields, []
        
    
    def Spectrum_processing(s):
        s = default_filters(s)
        return s

    def Float_prec_mz(s):
        s = s.set(db_params['db_mass_field'], float(s.get(db_params['db_mass_field'])))
        return s

    def Column_correction(table):
        drop_col = [i for i in table.columns if "Unnamed" in i]
        table.drop(drop_col, axis = 1, inplace = True)
        return table    

    def Database_table_mgf(database_mgf, db_params : dict) :
        """Extract the metadata from the MGF file
        """
        def Mz_field_choices():
            def Get_mz_list(s):
                return database_mgf[s].get(db_params['db_mass_field'])[0]
            def Get_mz(s):
                return database_mgf[s].get(db_params['db_mass_field'])
            if isinstance(database_mgf[0].get(db_params['db_mass_field']), list) : 
                return Get_mz_list
            else:
                return Get_mz

        def Filter_choices(db_params):
            def Three_filters(s):
                return [database_mgf[s].get(db_params['db_rt_field']), database_mgf[s].get(db_params['db_adduct_field']), database_mgf[s].get(db_params['db_formula_field_database'])]
            def RT_adduct(s):
                return [database_mgf[s].get(db_params['db_rt_field']), database_mgf[s].get(db_params['db_adduct_field'])]
            def RT_formula(s):
                return [database_mgf[s].get(db_params['db_rt_field']), database_mgf[s].get(db_params['db_formula_field_database'])]
            def Adduct_formula(s):
                return [database_mgf[s].get(db_params['db_adduct_field']), database_mgf[s].get(db_params['db_formula_field_database'])]
            def RT_only(s):
                return [database_mgf[s].get(db_params['db_rt_field'])]
            def Adduct_only(s):
                return [database_mgf[s].get(db_params['db_adduct_field'])]
            def Formula_only(s):
                return [database_mgf[s].get(db_params['db_formula_field_database'])]            
            def No_fields(s):
                return []
            if (db_params['db_rt_field'] != None) and (db_params['db_adduct_field'] != None) and (db_params['db_formula_field_database'] != None):
                return Three_filters, ['rt', 'adduct', 'formula']
            elif (db_params['db_rt_field'] != None) and (db_params['db_adduct_field'] != None):
                return RT_adduct, ['rt', 'adduct']
            elif (db_params['db_rt_field'] != None) and (db_params['db_formula_field_database'] != None):
                return RT_formula, ['rt', 'formula']
            elif (db_params['db_adduct_field'] != None) and (db_params['db_formula_field_database'] != None):
                return Adduct_formula, ['adduct', 'formula']
            elif db_params['db_rt_field'] != None:
                return RT_only, ['rt']
            elif db_params['db_adduct_field'] != None:
                return Adduct_only, ['adduct']
            elif db_params['db_formula_field_database'] != None:
                return Formula_only, ['formula']
            else:
                return No_fields, []
        database_table = list()
        Filter_fields, cols = Filter_choices(db_params)
        Mz_extractor = Mz_field_choices()
        col_names = ['name', 'mz', 'unique_field'] + cols + db_params['db_export_fields'] + ['ion_mode']
        print('Extracting database metadata...')
        for i in tqdm(range(len(database_mgf))):
            name = database_mgf[i].get(db_params['db_name_field'])
            mz = Mz_extractor(i)
            unique_field = database_mgf[i].get(db_unique_field)
            if unique_field == "" : unique_field = None
            filter_fields = Filter_fields(i)
            other_fields = [database_mgf[i].get(field) for field in db_params['db_export_fields']]
            other_fields = [None if f == "" else f for f in other_fields]
            ion_mode = database_mgf[i].get(db_mode_field)
            new_row = [name, mz, unique_field] + filter_fields + other_fields + [ion_mode]
            database_table.append((new_row))
        database_table = pd.DataFrame(database_table, columns = col_names)
        return database_table

    def Database_table_csv(database_csv, db_params : dict) :
        """Extract the metadata from the MGF file
        """
        def Filter_choices(db_params):
            def Three_filters(s):
                return [database_csv.loc[s, db_params['db_rt_field']],
                        database_csv.loc[s, db_params['db_adduct_field']],
                        database_csv.loc[s, db_params['db_formula_field_database']]]
            def RT_adduct(s):
                return [database_csv.loc[s, db_params['db_rt_field']],
                        database_csv.loc[s, db_params['db_adduct_field']]]
            def RT_formula(s):
                return [database_csv.loc[s, db_params['db_rt_field']],
                        database_csv.loc[s, db_params['db_formula_field_database']]]        
            def Adduct_formula(s):
                return [database_csv.loc[s, db_params['db_adduct_field']],
                        database_csv.loc[s, db_params['db_formula_field_database']]]     
            def RT_only(s):
                return [database_csv.loc[s, db_params['db_rt_field']]]
            def Adduct_only(s):
                return [database_csv.loc[s, db_params['db_adduct_field']]]
            def Formula_only(s):
                return [database_csv.loc[s, db_params['db_formula_field_database']]]
            def No_fields(s):
                return []
            if (db_params['db_rt_field'] != None) and (db_params['db_adduct_field'] != None) and (db_params['db_formula_field_database'] != None):
                return Three_filters, ['rt', 'adduct', 'formula']
            elif (db_params['db_rt_field'] != None) and (db_params['db_adduct_field'] != None):
                return RT_adduct, ['rt', 'adduct']
            elif (db_params['db_rt_field'] != None) and (db_params['db_formula_field_database'] != None):
                return RT_formula, ['rt', 'formula']
            elif (db_params['db_adduct_field'] != None) and (db_params['db_formula_field_database'] != None):
                return Adduct_formula, ['adduct', 'formula']
            elif db_params['db_rt_field'] != None:
                return RT_only, ['rt']
            elif db_params['db_adduct_field'] != None:
                return Adduct_only, ['adduct']
            elif db_params['db_formula_field_database'] != None:
                return Formula_only, ['formula']
            else:
                return No_fields, []
        database_table = list()
        Filter_fields, cols = Filter_choices(db_params)
        col_names = ['name', 'mass'] + cols + db_params['db_export_fields']
        print('Extracting database metadata...')
        for i in tqdm(database_csv.index):
            name = database_csv.loc[i, db_params['db_name_field']]
            mass = database_csv.loc[i, db_params['db_mass_field']]
            filter_fields = Filter_fields(i)
            other_fields = [database_csv.loc[i, field] for field in db_params['db_export_fields']]
            new_row = [name, mass] + filter_fields + other_fields
            database_table.append((new_row))
        database_table = pd.DataFrame(database_table, columns = col_names)
        return database_table

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
    
    import os
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.similarity import ModifiedCosine

    # Load parameters
    mgf_file_neg= params['neg_mgf']
    mgf_file_pos= params['pos_mgf']
    neg_csv_file= params['neg_csv']
    pos_csv_file= params['pos_csv']
    input_mgf_neg_path= params['neg_out_0']
    input_mgf_pos_path= params['pos_out_0']
    node_table_path= params['mix_out_4_1']
    out_path_full= params['mix_out_5_1']
    out_path_samples= params['mix_out_5_2']
    
    db_prefix= db_params['db_prefix']
    db_type= db_params['db_type']
    db_folder= db_params['db_folder']
    db_file= db_params['db_file']
    db_cosine= db_params['db_cosine']
    db_matched_peaks= db_params['db_matched_peaks']
    db_prec_error= db_params['db_prec_error']
    db_mass_error= db_params['db_mass_error']
    db_rt_error= db_params['db_rt_error']
    db_hits= db_params['db_hits']
    db_adduct_filter= db_params['db_adduct_filter']
    db_adduct_field= db_params['db_adduct_field']
    db_rt_filter= db_params['db_rt_filter']
    db_mode_field= db_params['db_mode_field']
    db_unique_field= db_params['db_unique_field']
    db_export_fields= db_params['db_export_fields']
    db_prefix= db_params['db_prefix']

    #Create folders and load tables
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
        edge_table = pd.read_csv(node_table_path + 'MIX_edges.csv', index_col = "Index")
        node_table = pd.read_csv(node_table_path + 'MIX_nodes.csv', index_col = "Index",
                                 dtype={'mz': 'float',
                                        'rt': 'float',
                                        'TIC' : 'float',
                                        'charge' : 'int',
                                        'mgf_index' : 'float',
                                        'status' : 'str',
                                        'Adnotation' : 'str',
                                        'ion_mode' : 'str',
                                        'feature_id' : 'float',
                                        'cluster_id' : 'int'})
    else:
        edge_table = pd.read_csv(out_path_full + 'MIX_edges.csv', index_col = "Index")
        node_table = pd.read_csv(out_path_full + 'MIX_nodes.csv', index_col = "Index",
                                 dtype={'mz': 'float',
                                        'rt': 'float',
                                        'TIC' : 'float',
                                        'charge' : 'int',
                                        'mgf_index' : 'float',
                                        'status' : 'str',
                                        'Adnotation' : 'str',
                                        'ion_mode' : 'str',
                                        'feature_id' : 'float',
                                        'cluster_id' : 'int'})
    if not os.path.isdir(out_path_samples) :
        os.mkdir(out_path_samples) 

    # Process tables and remove samples to be stored in a samples DF
    node_table['Adnotation'] = node_table['Adnotation'].replace({np.nan: None})
    samples = list(node_table.columns[node_table.columns.str.contains('Peak area')]) # safe samples in a separate DF
    samples_df = node_table[samples]
    node_table.drop(samples, axis = 1, inplace = True) # delete samples from the node table
    
    ###########################################################################
    if db_type == 'ion': # If database is MGF to dereplicate ion (MS/MS) data
    ###########################################################################
        modified_cosine = ModifiedCosine(tolerance=db_mass_error)
        Filter_fields_outer, filter_cols = Filter_choices_outer(db_params)
        
        # Load MGF files
        print('Loading NEG MGF file...')
        mgf_neg = list(load_from_mgf(input_mgf_neg_path + mgf_file_neg))
        mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
        print('Loading POS MGF file...')
        mgf_pos = list(load_from_mgf(input_mgf_pos_path + mgf_file_pos))
        mgf_pos = [Spectrum_processing(s) for s in mgf_pos]

        # Load the database file
        print('Loading database file and extracting data...')
        database_mgf = list(load_from_mgf(db_folder + db_file))
        database_mgf = [Float_prec_mz(s) for s in database_mgf]
        database_table = Database_table_mgf(database_mgf, db_params)
        
        # Start dereplication (cosine similarity)
        derep_table = list()
        print('Starting ion dereplication (cosine similarity)...')
        for i in tqdm(node_table.index):
            status = node_table.loc[i, "status_universal"]
            if status == "neutral" : 
                new_row = [i, None, None] + [None]*len(filter_cols) + [None]*len(db_export_fields) + [None, None, None]
                derep_table.append(new_row)                    
                continue
            ion_mz = node_table.loc[i, "mz"]
            ion_mgf_idx = int(node_table.loc[i, "mgf_index"])
            ion_mode = node_table.loc[i, "ion_mode"]
            hits = database_table[database_table['mz'].between(ion_mz - db_prec_error, ion_mz + db_prec_error, inclusive = True)].copy()
            
            # Ion mode filter
            if ion_mode == "NEG":
                hits = hits[hits['ion_mode'] == "negative"]
                exp_mgf = mgf_neg
            else:
                hits = hits[hits['ion_mode'] == "positive"]
                exp_mgf = mgf_pos
            
            # Optional filters : Adduct and RT
            if db_adduct_filter and (status == 'adduct'):
                adduct = node_table.loc[i, 'Adnotation']
                if adduct != None:
                    hits = hits[hits['adduct'] == adduct]
            if db_rt_filter:
                rt = node_table.loc[i, 'rt']
                hits = hits[hits['rt'].between(rt - db_rt_error, rt + db_rt_error, inclusive = True)]
            
            # Calculate cosine similarity if hit table is not empty
            similarity_list = list()
            for j in hits.index:
                score, n_matches = modified_cosine.pair(exp_mgf[ion_mgf_idx], database_mgf[j])
                mass_error = abs(ion_mz - hits.loc[j, "mz"])*1000
                prod = score * n_matches
                similarity_list.append((j, score, n_matches, prod, mass_error))
            similarity_list = pd.DataFrame(similarity_list, columns = ["index", "cos", "matches", "prod", "error_mDa"])
            similarity_list.set_index('index', inplace = True)
            hits = hits.join(similarity_list, hits.index)
            
            # Filter using cosine and matched peaks thresholds
            hits = hits[hits['cos'] >= db_cosine]
            hits = hits[hits['matches'] >= db_matched_peaks]

            # If no hits at this point, export None and continue
            if len(hits) == 0 : 
                new_row = [i, None, None] + [None]*len(filter_cols) + [None]*len(db_export_fields) + [None, None, None]
                derep_table.append(new_row)    
                continue
            
            # Filter using unique field:
            unique_mols = hits['unique_field'].dropna().unique()
            mol_idx = hits.index[hits['unique_field'].isnull()].tolist()
            for mol in unique_mols:
                tmp_table_1 = hits[hits['unique_field'] == mol].copy()
                tmp_table_1.sort_values('prod', ascending = False, inplace = True)
                mol_idx.append(tmp_table_1.index[0])
            hits = hits.loc[mol_idx]
            
            # Keep the top db_hits hit in the hits table
            hits.sort_values('prod', ascending = False, inplace = True)
            hits = hits.iloc[:db_hits]
            hits = hits.fillna('')
            
            # Report the results in the derep_table list
            tmp_name = '|'.join(hits['name'])
            tmp_mz = '|'.join(hits['mz'].round(4).astype(str))
            tmp_filter_fields = ['|'.join(hits[col].astype(str)) for col in filter_cols]
            tmp_export_fields = ['|'.join(hits[f].astype(str)) for f in db_export_fields]
            tmp_cos = '|'.join(hits['cos'].round(2).astype(str))
            tmp_matches = '|'.join(hits['matches'].astype(str))
            tmp_error = '|'.join(hits['error_mDa'].round(1).astype(str))
            new_row = [i, tmp_name, tmp_mz] + tmp_filter_fields + tmp_export_fields + [tmp_cos, tmp_matches, tmp_error]
            derep_table.append(new_row)
        
        cols = ['index', 'name', 'mass'] + filter_cols + db_export_fields + ['cos', 'matches', 'error_mDa']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        cols.remove(db_prefix + '_index')
        if len(set(cols).intersection(node_table.columns)) > 0 :
            node_table.drop(cols, axis = 1, inplace = True)
        node_table = node_table.join(derep_table, node_table.index)
        
        # Propagate ion dereplications to neutrals
        neutrals_idx = node_table.index[node_table['status_universal'] == "neutral"]
        print('Propagating ion annotations to neutral nodes...')
        for i in tqdm(neutrals_idx):
            tmp_ions = edge_table['node_2'][edge_table['node_1'] == i].tolist()
            tmp_nodes = node_table.loc[tmp_ions]
            tmp_nodes = tmp_nodes[~tmp_nodes[db_prefix + '_mass'].isnull()]
            if len(tmp_nodes) == 0 : continue
            derep_table = list()
            for j in tmp_nodes.index:
                for k in range(len(tmp_nodes.loc[j, db_prefix + '_mass'].split('|'))):
                    new_row = list()
                    for col in cols:
                        new_row.append(tmp_nodes.loc[j, col].split('|')[k])
                    derep_table.append(new_row)
            derep_table = pd.DataFrame(derep_table, columns = cols)

            derep_table['prod'] = derep_table[db_prefix + '_cos'].astype(float) * derep_table[db_prefix + '_matches'].astype(int)
            unique_mols = derep_table[db_prefix + '_' + db_unique_field].dropna().unique()
            mol_idx = derep_table.index[derep_table[db_prefix + '_' + db_unique_field].isnull()].tolist()
            for mol in unique_mols:
                tmp_table_1 = derep_table[derep_table[db_prefix + '_' + db_unique_field] == mol].copy()
                tmp_table_1.sort_values('prod', ascending = False, inplace = True)
                mol_idx.append(tmp_table_1.index[0])
            derep_table = derep_table.loc[mol_idx]
            derep_table.sort_values("prod", ascending = False, inplace = True)
            derep_table = derep_table.iloc[:db_hits]
            derep_table = derep_table.fillna('')
            for col in cols:
                node_table.loc[i, col] = '|'.join(derep_table[col])
            if db_adduct_field != None:
                node_table.loc[i, db_adduct_field] = None
            
    #############################################################################
    elif db_type == "neutral": # If database is CSV to dereplicate molecule nodes
    #############################################################################
        # Load the database file
        print('Loading database file and extracting data...')
        if db_file[-3:] == 'tsv':
            sep = '\t'
        elif db_file[-3:] == 'csv':
            sep = ','
        if db_params['db_index_field'] != None :
            database_csv = pd.read_csv(db_folder + db_file, sep = sep, index_col = db_params['db_index_field'])
        else:
            database_csv = pd.read_csv(db_folder + db_file, sep = sep)
        database_table = Database_table_csv(database_csv, db_params)
        
        # Start dereplication (on neutrals only)
        derep_table = list()
        print('Starting neutral dereplication...')
        for i in tqdm(node_table.index):
            if node_table.loc[i, "status_universal"] == "neutral":
                mass = node_table.loc[i, "mz"]
                hits = database_table[database_table['mass'].between(mass - db_prec_error, mass + db_prec_error, inclusive = True)].copy()
                if db_adduct_filter :
                    rt = node_table.loc[i, 'rt']
                    hits = hits[hits['rt'].between(rt - db_rt_error, rt + db_rt_error, inclusive = True)]
                if len(hits) == 0 : 
                    new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                    derep_table.append(new_row)    
                    continue
                hits['error_mDa'] = abs(mass - hits['mass'])*1000
                hits.sort_values('error_mDa', ascending = True, inplace = True)
                hits = hits.iloc[:db_hits]
                tmp_names = '|'.join(hits['name'])
                tmp_masses = '|'.join(hits['mass'].round(4).astype(str))
                tmp_error = '|'.join(hits['error_mDa'].round(1).astype(str))
                other_fields = ['|'.join(hits[field].astype(str)) for field in db_params['db_export_fields']]
                new_row = [i, tmp_names, tmp_masses, tmp_error]  + other_fields
                derep_table.append(new_row)
            else:
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
        cols = ['index', 'name', 'mass', 'error_mDa'] + db_params['db_export_fields']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        node_table = node_table.join(derep_table, node_table.index)
        
    #############################################################################
    elif db_type == "formula": # CSV database dereplication based on chem formulas for neutrals and some ions
    #############################################################################

        # Load the database file
        print('Loading database file and extracting data...')
        if db_file[-3:] == 'tsv':
            sep = '\t'
        elif db_file[-3:] == 'csv':
            sep = ','
        if db_params['db_index_field'] != None :
            database_csv = pd.read_csv(db_folder + db_file, sep = sep, index_col = db_params['db_index_field'])
        else:
            database_csv = pd.read_csv(db_folder + db_file, sep = sep)
        database_table = Database_table_csv(database_csv, db_params)
        
        # Start dereplication (on neutrals only)
        derep_table = list()
        print('Starting formula dereplication...')
        for i in tqdm(node_table.index):
            formula = node_table.loc[i, "formula"]
            if type(formula) == float:
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
                continue
            hits = database_table[database_table['formula'] == formula].copy()
            hits = hits.iloc[:db_hits]
            if len(hits) == 0 : 
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
                continue                
            tmp_names = '|'.join(hits['name'])
            tmp_masses = str(hits['mass'].iloc[0].round(4))
            tmp_formula = hits['formula'].iloc[0]
            other_fields = ['|'.join(hits[field].astype(str)) for field in db_params['db_export_fields']]
            new_row = [i, tmp_names, tmp_masses, tmp_formula]  + other_fields
            derep_table.append(new_row)

        cols = ['index', 'name', 'mass', 'formula'] + db_params['db_export_fields']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        node_table = node_table.join(derep_table, node_table.index)
    
    ###########################################################################
    # Prepare data to be exported
    ###########################################################################
    
    # Round values where necessary and return the samples to the dataframe
    node_table = node_table.merge(samples_df, left_index = True, right_index = True)   
    node_table['mz'] = node_table['mz'].round(4)
    node_table['rt'] = node_table['rt'].round(3)
    
    # Export data
    node_table.to_csv(out_path_full + 'MIX_nodes.csv', index_label = "Index")
    edge_table.to_csv(out_path_full + 'MIX_edges.csv', index_label = "Index")
    
    if db_params['db_export_samples'] : 
        Samplewise_export(neg_csv_file = input_mgf_neg_path + neg_csv_file,
                          pos_csv_file = input_mgf_pos_path + pos_csv_file,
                          out_path = out_path_samples,
                          merged_edge_table = edge_table,
                          merged_node_table = node_table)
    return
