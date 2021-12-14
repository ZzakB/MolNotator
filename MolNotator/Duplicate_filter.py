def Duplicate_filter(params : dict, ion_mode : str):

    def Spectrum_processing(s):
        s = default_filters(s)
        return s
    def Mgf_data_extractor(mgf_file) :
        """Takes an MGF file to extract, in the form of a dataframe, each ion's
        feature ID, m/z value, retention time, charge and TIC from the MS2 spectrum.
        
        """
        node_table = pd.DataFrame(index = range(len(mgf_file)), columns = ['feature_id',
                                  'mz', 'rt', 'TIC', 'charge'])
        node_table = node_table.fillna(0.0)
        print('Extracting MGF metadata...')
        for i in tqdm(range(len(mgf_file))) :
            node_table.loc[i, 'feature_id'] = int(mgf_file[i].get('feature_id'))
            node_table.loc[i, 'mz'] = mgf_file[i].get('pepmass')[0]
            node_table.loc[i, 'rt'] = float(mgf_file[i].get('rtinseconds'))
            node_table.loc[i, 'TIC'] = mgf_file[i].peaks.intensities.sum()
            node_table.loc[i, 'charge'] = mgf_file[i].get('charge')[0]
        node_table['feature_id'] = node_table['feature_id'].astype(int)
        return node_table
    
    def Duplicate_finder(node_table, mgf_file):
        modified_cosine = ModifiedCosine(tolerance=mass_error)
        ions_idx = list(node_table.index)
        total_ions = len(ions_idx)
        duplicate_table = pd.DataFrame(columns = ['kept', 'dropped'])
        while len(ions_idx) > 1:
            tmp_duplicate_table = list()
            perc = round((1-(len(ions_idx)/total_ions))*100,1)
            sys.stdout.write("\rFinding duplicates : {0}% ions processed".format(perc))
            sys.stdout.flush()
            ion_pool = [ions_idx[0]]
            pool_size = 0
            pool_table = node_table.loc[ion_pool]
            while len(ion_pool) != pool_size:
                pool_size = len(ion_pool)
                min_rt = pool_table['rt'].min()
                max_rt = pool_table['rt'].max()
                min_mz = pool_table['mz'].min()
                max_mz = pool_table['mz'].max()
                pool_table = node_table[node_table['rt'].between(min_rt - rt_error, max_rt + rt_error)]
                pool_table = pool_table[pool_table['mz'].between(min_mz - mass_error, max_mz + mass_error)]
                ion_pool = list(pool_table.index)
            ions_idx = list(set(ions_idx) - set(ion_pool))
            cos_table = list()
            while len(ion_pool) > 1 :
                ion_1 = ion_pool[0]
                rt_1 = pool_table.loc[ion_1, 'rt']
                mz_1 = pool_table.loc[ion_1, 'mz']
                ion_pool.remove(ion_1)
                for ion_2 in ion_pool:
                    rt_2 = pool_table.loc[ion_2, 'rt']
                    mz_2 = pool_table.loc[ion_2, 'mz']
                    score, n_matches = modified_cosine.pair(mgf_file[ion_1],
                                                            mgf_file[ion_2])
                    cos_table.append((ion_1, ion_2, score, n_matches, abs(mz_1 - mz_2), abs(rt_1 - rt_2)))
            cos_table = pd.DataFrame(cos_table, columns = ['ion_1', 'ion_2', 'cos', 'matched_peaks', 'd_mz', 'd_rt'])
            if len(cos_table) == 0 : continue
            ion_pool = list(pool_table.index)
            while len(ion_pool) > 0 :
                ion_seed = pool_table.loc[ion_pool, 'TIC'].idxmax()
                seed_rt = pool_table.loc[ion_seed, "rt"]
                seed_mz = pool_table.loc[ion_seed, "mz"]

                candidates = cos_table.index[cos_table['ion_1'] == ion_seed].tolist() + cos_table.index[cos_table['ion_2'] == ion_seed].tolist()
                candidates = cos_table.loc[candidates].index[cos_table.loc[candidates, 'cos'] >= 0.8]
                candidates = cos_table.loc[candidates, 'ion_1'].tolist() + cos_table.loc[candidates, 'ion_2'].tolist()
                candidates = list(set(candidates))
                if len(candidates) == 0 : 
                    tmp_duplicate_table.append((ion_seed, []))
                    ion_pool.remove(ion_seed)
                    continue
                candidates.remove(ion_seed)
                candidates.sort()
                candidates = pool_table.loc[candidates]
                
                candidates = candidates[candidates['rt'].between(seed_rt - rt_error, seed_rt + rt_error, inclusive = True)]
                candidates = candidates[candidates['mz'].between(seed_mz - mass_error, seed_mz + mass_error, inclusive = True)]
                candidates = candidates.index.tolist()
                ion_pool.remove(ion_seed)
                ion_pool = list(set(ion_pool) - set(candidates))
                tmp_duplicate_table.append((ion_seed, candidates))
            tmp_duplicate_table = pd.DataFrame(tmp_duplicate_table, columns = ['kept', 'dropped'])
            # Ion affinity
            contested_ions = list(set(flatten(tmp_duplicate_table['dropped'])))
            contested_ions.sort()
            conflict_table = list()
            for i in contested_ions :
                conflicts = list()
                for j in tmp_duplicate_table.index:
                    if i in tmp_duplicate_table.loc[j, "dropped"]:
                        conflicts.append(j)
                if len(conflicts) > 1 :
                    tmp_scores = list()
                    for j in conflicts:
                        score, n_matches = modified_cosine.pair(mgf_file[i],
                                                                mgf_file[tmp_duplicate_table.loc[j, "kept"]])
                        tmp_scores.append(score*n_matches)
                    selected_ions = tmp_scores.index(max(tmp_scores))
                    selected_ions = conflicts[selected_ions]
                    conflict_table.append((i, conflicts, selected_ions))
            conflict_table = pd.DataFrame(conflict_table, columns = ['ion', 'conflicts', 'selected'])
            
            for i in conflict_table.index:
                for j in conflict_table.loc[i, "conflicts"]:
                    if j == conflict_table.loc[i, "selected"] : continue
                    tmp_duplicate_table.loc[j, "dropped"].remove(conflict_table.loc[i, "ion"])
            duplicate_table = duplicate_table.append(tmp_duplicate_table, ignore_index = True)
        
        dup_size = list()
        for i in duplicate_table.index:
            dup_size.append(len(duplicate_table.loc[i, "dropped"]))
        duplicate_table['dup_size'] = dup_size
        duplicate_table = duplicate_table[duplicate_table['dup_size'] > 0]
        duplicate_table.drop('dup_size', axis = 1, inplace = True)
        
        return duplicate_table

    import os
    import sys
    import pandas as pd
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.similarity import ModifiedCosine
    from matchms.exporting import save_as_mgf
    from pandas.core.common import flatten
    from tqdm import tqdm

    # Load parameters
    in_path= params['mzmine_out']
    mass_error= params['df_mass_error']
    rt_error= params['df_rt_error']
    mzmine_suffix= params['mzmine_suffix']
    skip= params['df_skip']
    
    if ion_mode == "NEG":
        csv_file = params['neg_csv']
        mgf_file= params['neg_mgf']
        out_path= params['neg_out_0']
    elif ion_mode == "POS" :
        csv_file = params['pos_csv']
        mgf_file= params['pos_mgf']
        out_path= params['pos_out_0']

    # Save mgf and csv file names for input
    mgf_name = mgf_file
    csv_name = csv_file
    
    # create output dir:
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    
    # Load MZmine mgf and csv files
    print("Loading MGF and CSV files...")
    mgf_file = list(load_from_mgf(f'{in_path}{mgf_file}'))
    csv_file = pd.read_csv(f'{in_path}{csv_file}', index_col = "row ID")
    
    # Format columns with ion modes:
    new_cols = csv_file.columns.tolist()
    new_cols = [ion_mode + "_" + col if (mzmine_suffix in col and col[:4] != f"{ion_mode}_") else col for col in new_cols]
    csv_file.columns = new_cols

    # If this step must be skipped :
    if skip :
        csv_file.to_csv(f'{out_path}{csv_name}')
        save_as_mgf(mgf_file, f'{out_path}{mgf_name}')
        return
    
    # Extract data from the MGF file
    node_table = Mgf_data_extractor(mgf_file)
    
    # Filtering the MGF file
    print('Filtering MGF file...')
    mgf_file = [Spectrum_processing(s) for s in mgf_file]
    
    # Get duplicates
    duplicate_table = Duplicate_finder(node_table, mgf_file)
    
    # Transfer data from removed ions to kept ions
    samples = csv_file.columns
    samples = list(samples.drop(['row m/z', 'row retention time']))
    print("\nTransfering data from removed ions...")
    for i in tqdm(duplicate_table.index):
        idx_i = duplicate_table.loc[i, "kept"]
        id_i = node_table.loc[idx_i, "feature_id"]
        duplicates = duplicate_table.loc[i, "dropped"]
        samples_i = csv_file.loc[id_i, samples] > 0
        samples_i = csv_file.loc[id_i, samples][samples_i]
        for idx_j in duplicates:
            id_j = node_table.loc[idx_j, "feature_id"]
            samples_j = csv_file.loc[id_j, samples] > 0
            samples_j = csv_file.loc[id_j, samples][samples_j]
            for sample in samples_j.index:
                if sample in samples_i.index:
                    csv_file.loc[idx_i, sample] = max(samples_i[sample], samples_j[sample])
                else:
                    csv_file.loc[idx_i, sample] = samples_j[sample]
    
    # Filter kept ions to eliminate some ions that were missed

    dropped_ions = list(flatten(duplicate_table['dropped']))
    kept_ions = list(set(node_table.index) - set(dropped_ions))
    dropped_fids = [node_table.loc[i, 'feature_id'] for i in dropped_ions]
    
    # Export the data
    mgf_file_new = [mgf_file[i] for i in kept_ions]
    csv_file_new = csv_file.drop(dropped_fids)
    print('Exporting MGF and CSV files...')
    save_as_mgf(mgf_file_new, f'{out_path}{mgf_name}')
    csv_file_new.to_csv(f'{out_path}{csv_name}')
    perc = round(100*(len(dropped_ions)/len(mgf_file)),1)
    print('Export finished.')
    print(f'{len(dropped_ions)} ions removed out of {len(mgf_file)} ({perc}%)')
    return
