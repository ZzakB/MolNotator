"""fragnotator.py - fragnotator module for MolNotator"""
import os
import pandas as pd
from matchms.importing import load_from_mgf
from MolNotator.others.spectrum_extractor import spectrum_extractor
from MolNotator.others.fragnotator_edge_table import fragnotator_edge_table
from MolNotator.others.singleton_edges import singleton_edges
from MolNotator.others.reindexer import reindexer

def fragnotator(params : dict, ion_mode : str):
    """
    Takes as input a spectrum file from a single sample to connect in-source 
    fragment ions to their in-source precursor.

    Parameters
    ----------
    params : dict
        Dictionary containing the global parameters for the process..
    ion_mode : str
        Either "POS" or "NEG", ion mode for the data..

    Returns
    -------
    Writes node and edge tables for each sample 

    """

    # Load parameters
    fragnotator_table = params['fn_fragtable']
    mass_error = params['fn_mass_error']
    if ion_mode == "NEG" :
        in_path= params['neg_out_1']
        out_path= params['neg_out_2']
    elif ion_mode == "POS":
        in_path= params['pos_out_1']
        out_path= params['pos_out_2']
    else:
        print('Ion mode must be either "NEG" or "POS"')
        return

    # Create out dir
    if not os.path.exists(f'{out_path}') :
        os.mkdir(f'{out_path}')

    # Get file list and start processing for each sample
    files = pd.Series(os.listdir(f'{in_path}'))
    iteration = 1
    for file_name in files :
        print("========= File " + str(iteration) + " out of " + str(len(files)) + " =========")
        # Get the node and edge tables :
        out_name_edge = f'{out_path}/{file_name.replace(".mgf" , "_edges.csv")}'
        out_name_node  = f'{out_path}/{file_name.replace(".mgf" , "_nodes.csv")}'
        
        # Load the MGF file
        spectrum_list = list(load_from_mgf(f'{in_path}/{file_name}'))
        
        # Create the node_table table : 
        print('Extracting spectrum data...')
        node_table = spectrum_extractor(spectrum_list, True)
        
        # Reindex the node table
        node_table = reindexer(node_table, params)

        # Coerce RT dtype in node_table
        node_table[params["rt_field"]] = node_table[params["rt_field"]].astype(float)
        
        # Get the parent-fragment pairs:
        print('Pairing in-source precursors and fragments...')
        edge_table = fragnotator_edge_table(node_table, spectrum_list, params)
        
        # Add singleton nodes to edge_table
        edge_table = singleton_edges(node_table, edge_table)
            
        # Update the node table with the Parent / fragment / unpaired status
        node_table['status'] = ['singleton']*len(node_table.index)
        tmp_table = edge_table[edge_table['status'] == "frag_edge"]
        prec_nodes = set(tmp_table['node_1'])
        frag_nodes = set(tmp_table['node_2'])
        prec_nodes = prec_nodes - frag_nodes
        prec_nodes = list(prec_nodes)
        frag_nodes = list(frag_nodes)
        node_table.loc[prec_nodes, "status"] = "precursor"
        node_table.loc[frag_nodes, "status"] = "fragment"
        
        # Add neutral losses annotations to the edge table
        frag_table = pd.read_csv("./params/" + fragnotator_table, sep = '\t')
        edge_table['Fragnotation'] = [None]*len(edge_table.index)
        for i in frag_table.index:
            low_mass = frag_table.loc[i, 'mass'] - mass_error
            high_mass = frag_table.loc[i, 'mass'] + mass_error
            temp_edge_table = edge_table[edge_table['mz_gap'].between(low_mass, high_mass, inclusive = "both")]
            for j in temp_edge_table.index : 
                edge_table.loc[j, 'Fragnotation'] = frag_table.loc[i, 'loss']
        
        # Add an all_annotations columns, which will sum up all possible annotations
        all_annotations = list()
        for i in edge_table.index:
            if edge_table.loc[i, "Fragnotation"] == None:
                all_annotations.append(round(edge_table.loc[i, 'mz_gap'], 4))
            else:
                all_annotations.append(edge_table.loc[i, 'Fragnotation'])
        edge_table['All_annotations'] = all_annotations
        
        # Round values
        node_table[params["mz_field"]] = node_table[params["mz_field"]].astype(float).round(4)
        node_table[params["rt_field"]] = node_table[params["rt_field"]].astype(float).round(2)
        edge_table['matching_score'] = edge_table['matching_score'].round(2)
        edge_table['rt_gap'] = edge_table['rt_gap'].round(2)
        edge_table['mz_gap'] = edge_table['mz_gap'].round(4)
        
        
        # Export the edge and node tables :
        edge_table.to_csv(out_name_edge, index_label = "Index")
        node_table.to_csv(out_name_node, index_label = params['index_col'])
        iteration += 1
    return
