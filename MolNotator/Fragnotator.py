def Fragnotator(params : dict, ion_mode : str):
    """Takes as input an MGF file from a single sample to connect in-source 
    fragment ions to their in-source precursor.
    
    """
    def Rt_slicer(rt, rt_error, ion_id, input_table) :
        """Given a dataframe with ion data including retention time (rt), will 
        return a subset of that dataframe including only ions co-detected with
        a target ion (using its retention time and a retention time window)
        
        """
        rt_low = rt - rt_error
        rt_high = rt + rt_error
        sliced_table = input_table[input_table['rt'].between(rt_low,
                                  rt_high, inclusive = True)].copy()
        return sliced_table.drop(ion_id)
    
    
    def Mgf_data_extractor(mgf_file) :
        """Takes an MGF file to extract, in the form of a dataframe, each ion's
        feature ID, m/z value, retention time, charge and TIC from the MS2 spectrum.
        
        """
        node_table = pd.DataFrame(index = range(len(mgf_file)), columns = ['feature_id',
                                  'mz', 'rt', 'TIC', 'charge'])
        node_table = node_table.fillna(0.0)
        for i in range(len(mgf_file)) :
            node_table.loc[i, 'feature_id'] = int(mgf_file[i].get('feature_id'))
            node_table.loc[i, 'mz'] = mgf_file[i].get('pepmass')[0]
            node_table.loc[i, 'rt'] = float(mgf_file[i].get('rtinseconds'))
            node_table.loc[i, 'TIC'] = mgf_file[i].peaks.intensities.sum()
            node_table.loc[i, 'charge'] = mgf_file[i].get('charge')[0]
        #node_table.set_index("ionIndex", drop = False, inplace = True)
        return(node_table)
    
    def Fragnotator_edge_table(mass_error, rt_error, min_shared_peaks, score_threshold):
        """Returns an edge table which connects all in-source precursors (node_1)
        to their in-sourc fragments (node_2). For an ion to be considered an
        in-source fragment, it must be co-detected with its ion-source precursor
        within an rt_error window, it must share with it an absolute minimum 
        number of peaks (min_shared_peaks) and relative minimum number of peaks
        (score_threshold). 
        """
        node_1 = list()
        node_2 = list()
        shared_peaks_list = list()
        total_peaks_list = list()
        matching_score_list = list()
        rt_gap_list = list()
        mz_gap_list = list()
        # Each ion is considered to be a parent and their fragment ions are searched
        for i in tqdm(node_table.index) :
            # Get attributes for the current parent ion
            ion1_rt = node_table['rt'][i]
            ion1_mz = node_table['mz'][i]
            
            # Get the m/z values of the MS/MS
            ion1_msms = pd.Series(mgf_file[i].peaks.mz)
            
            # Select the coeluting ions using the retention time window (rt_error)
            # and delete all ions with m/z above the current ion1_mz (fragments
            # should have lower m/z than the parent ion)
            lowermass_table = Rt_slicer(ion1_rt, rt_error, i, node_table)
            lowermass_table = lowermass_table[lowermass_table['mz'] < ion1_mz]
            if len(lowermass_table.index) == 0 : continue
            
            # In the lowermass_table, select the fragment candidates, i.e. with m/z
            # values present in the current parent ion MS/MS (ion1MSMS)
            frag_list = list()
            
            # Iterature through all the lowermass_table to confirm the candidates
            for j in lowermass_table.index:
                # Load the fragment candidate
                ion2_mz = node_table['mz'][j]
                ion2_mz_low = ion2_mz - mass_error
                ion2_mz_high = ion2_mz + mass_error
                match = ion1_msms.between(ion2_mz_low, ion2_mz_high, inclusive = True)
                if match.sum() > 0 : # if the frag candidate m/z is found in MSMS:
                    ion2_msms = mgf_file[j].peaks.mz # extract frag's MSMS to run tests
                    matched_peaks = 0
                    total_peaks = list(ion1_msms)
                    for frag in ion2_msms : # find the number of matched peaks
                        frag_low = frag - mass_error
                        frag_high = frag + mass_error
                        frag_found = ion1_msms.between(frag_low, frag_high, inclusive = True).sum()
                        if frag_found > 0 :
                            matched_peaks += 1
                        else :
                            total_peaks.append(frag)
                    if matched_peaks >= min_shared_peaks : # if number of matches above threshold
                        total_peaks = pd.Series(total_peaks)[total_peaks <= ion2_mz_high]
                        matching_score = round(matched_peaks / len(total_peaks),2)
                        if matching_score >= score_threshold:
                            frag_list.append(True)
                            lowermass_table.loc[j, 'matched_peaks'] = matched_peaks
                            lowermass_table.loc[j, 'total_peaks'] = len(total_peaks)
                            lowermass_table.loc[j, 'matching_score'] = matching_score
                            lowermass_table.loc[j, 'rt_gap'] = ion1_rt - node_table['rt'][j]
                            lowermass_table.loc[j, 'mz_gap'] = ion1_mz - node_table['mz'][j]
                        else : frag_list.append(False)
                    else : frag_list.append(False)
                else : frag_list.append(False)
                
            lowermass_table = lowermass_table[frag_list] # Keep only confirmed frags
            if len(lowermass_table.index) == 0 : continue
        
            # Report the selected fragment candidates to the results lists
            for j in lowermass_table.index :
                node_1.append(i)
                node_2.append(j)
                shared_peaks_list.append(lowermass_table['matched_peaks'][j])
                total_peaks_list.append(lowermass_table['total_peaks'][j])
                matching_score_list.append(lowermass_table['matching_score'][j])
                rt_gap_list.append(lowermass_table['rt_gap'][j])    
                mz_gap_list.append(lowermass_table['mz_gap'][j])
            
        edge_table = list(zip(node_1, node_2, shared_peaks_list, total_peaks_list,
                              matching_score_list, rt_gap_list, mz_gap_list))
        edge_table = pd.DataFrame(edge_table, columns = ['node_1', 'node_2',
                                                         'matched_peaks', 'total_peaks',
                                                         'matching_score', 'rt_gap',
                                                         'mz_gap'])
        return(edge_table)
    
    def Singleton_edges(edge_table):
        """Adds singleton nodes (neither in-source precursors nor fragments) 
        to the edge table produced by the Fragnotator_edge_table 
        function. Their are represented by the same ion in the "node_1" and 
        "node_2" columns.
        """
        paired_nodes = set(list(edge_table['node_1']) + list(edge_table['node_2']))
        singleton_nodes = list(set(node_table.index) - paired_nodes)
        singleton_edge_table = pd.DataFrame(index = range(len(singleton_nodes)), 
                                            columns = edge_table.columns)
        singleton_edge_table['node_1'] = singleton_nodes
        singleton_edge_table['node_2'] = singleton_nodes
        singleton_edge_table['rt_gap'] = [0.0]*len(singleton_edge_table)
        singleton_edge_table['mz_gap'] = [0.0]*len(singleton_edge_table)
        edge_table['status'] = ["frag_edge"]*len(edge_table)
        singleton_edge_table['status'] = ["singleton"]*len(singleton_edge_table)
        edge_table = edge_table.append(singleton_edge_table, ignore_index = False, sort = False)
        return edge_table.reset_index(drop = True)
    
    def Fragnotator_status(edge_table):
        """Updates the node table with the status of each ion, as "precursor",
        "fragment" or "singleton".
        """
        node_table['status'] = ['singleton']*len(node_table.index)
        tmp_table = edge_table[edge_table['status'] == "frag_edge"]
        for i in node_table.index:
            if (tmp_table['node_2'] == i).sum() > 0: 
                node_table.loc[i, 'status'] = "fragment"
            elif (tmp_table['node_1'] == i).sum():
                node_table.loc[i, 'status'] = "precursor"
        return
    
    def Fragnotation():
        """Annotates edges in the edge_table by using a tsv input file (by default,
        "fragnotator_table.tsv" in the params.yaml file). The tsv file contains
        common neutral losses and their mass values. More complex losses can
        be added by the user if necessary.
        """
        import pandas as pd
        frag_table = pd.read_csv("./params/" + fragnotator_table, sep = '\t')
        edge_table['Fragnotation'] = [None]*len(edge_table.index)
        for i in frag_table.index:
            low_mass = frag_table.loc[i, 'mass'] - mass_error
            high_mass = frag_table.loc[i, 'mass'] + mass_error
            temp_edge_table = edge_table[edge_table['mz_gap'].between(low_mass, high_mass, inclusive = True)]
            for j in temp_edge_table.index : 
                edge_table.loc[j, 'Fragnotation'] = frag_table.loc[i, 'loss']
        return(edge_table)
    
    import os
    import pandas as pd
    from matchms.importing import load_from_mgf
    from tqdm import tqdm
    
    # Load parameters
    min_shared_peaks= params['fn_matched_peaks']
    score_threshold = params['fn_score_threshold']
    mass_error= params['fn_mass_error']
    rt_error= params['fn_rt_error']
    fragnotator_table= params['fn_fragtable']    
    if ion_mode == "NEG" :
        in_path= params['neg_out_1']
        out_path= params['neg_out_2']
    elif ion_mode == "POS":
        in_path= params['pos_out_1']
        out_path= params['pos_out_2']
    else:
        print('Ion mode must be either "NEG" or "POS"')
        return

    
    if not os.path.exists(f'{out_path}') :
        os.mkdir(f'{out_path}')

    files = pd.Series(os.listdir(f'{in_path}'))

    iteration = 1
    for file_name in files :
        print("========= File " + str(iteration) + " out of " + str(len(files)) + " =========")
        # Get the node and edge tables :
        out_name_edge = f'{out_path}/{file_name.replace(".mgf" , "_edges.csv")}'
        out_name_node  = f'{out_path}/{file_name.replace(".mgf" , "_nodes.csv")}'
        
        # Load the MGF file
        mgf_file = list(load_from_mgf(f'{in_path}/{file_name}'))
        #mgfFile = list(libmetgem.mgf.read(file_fullname, ignore_unknown = False))
        
        # Create the node_table table : 
        node_table = Mgf_data_extractor(mgf_file)
        
        # Get the parent-fragment pairs:
        edge_table = Fragnotator_edge_table(mass_error, rt_error, min_shared_peaks, score_threshold)
        
        # Add singleton nodes to edge_table
        edge_table = Singleton_edges(edge_table)
            
        # Update the node table with the Parent / fragment / unpaired status
        Fragnotator_status(edge_table)
        
        # Add some annotation to the edge table
        edge_table = Fragnotation()
        
        all_annotations = list()
        for i in edge_table.index:
            if edge_table.loc[i, "Fragnotation"] == None:
                all_annotations.append(round(edge_table.loc[i, 'mz_gap'], 4))
            else:
                all_annotations.append(edge_table.loc[i, 'Fragnotation'])
        edge_table['All_annotations'] = all_annotations
        
        # Export the edge and node tables :
        edge_table.to_csv(out_name_edge, index_label = "Index")
        node_table.to_csv(out_name_node, index_label = "Index")
        iteration += 1
