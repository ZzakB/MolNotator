"""fragnotator_edge_table.py - fragnotator_edge_table module for fragnotator"""
import pandas as pd
from tqdm import tqdm
from MolNotator.others.rt_slicer import rt_slicer

def fragnotator_edge_table(node_table, mgf_file, params):
    """
    Finds precursor-fragment ion pairs for in-source fragmentation, using the
    metadata in a node table and the spectra in the spectrum file.
    Parameters
    ----------
    node_table : pandas.DataFrame
        Dataframe containing metadata from the fragnotator module.
    mgf_file : list
        List of matchms.Spectrum objects from the fragnotator module.
    params : dict
        Dictionary containing the global parameters for the process.
    Returns
    -------
    edge_table : pandas.DataFrame
        Dataframe containing the precursor-fragment ion pairs.
    """

    # Get parameters
    score_threshold = params['fn_score_threshold']
    min_shared_peaks = params['fn_matched_peaks']
    rt_error = params['fn_rt_error']
    mass_error = params['fn_mass_error']
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    

    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        rt_error = rt_error/60

    # For each ion, search fragment candidates
    edge_table = list()
    for i in tqdm(node_table.index) :
        
        # Get ion 1 data (precursor)
        ion1_spec_id = node_table.loc[i, "spec_id"]
        ion1_rt = node_table.loc[i, rt_field]
        ion1_mz = node_table.loc[i, mz_field]
        ion1_msms = pd.Series(mgf_file[ion1_spec_id].peaks.mz)
        
        # Find fragment candidate ions (below mz, similar RT)
        candidate_table = rt_slicer(ion1_rt, rt_error, i, node_table, rt_field)
        candidate_table = candidate_table[candidate_table[mz_field] < ion1_mz]
        
        # If no candidates are found, skip
        if len(candidate_table.index) == 0 : continue
        
        # Candidates must share their precursor ion with the precursor in MSMS
        for j in candidate_table.index:
            
            # Get ion 2 data (fragment)
            ion2_spec_id = node_table.loc[j, "spec_id"]
            ion2_mz = node_table.loc[j, mz_field]
            ion2_mz_low = ion2_mz - mass_error
            ion2_mz_high = ion2_mz + mass_error
            match = ion1_msms.between(ion2_mz_low, ion2_mz_high, inclusive = "both")
            if match.sum() > 0 : # if the frag candidate m/z is found in MSMS:
                ion2_msms = mgf_file[ion2_spec_id].peaks.mz
                matched_peaks = 0
                total_peaks = list(ion1_msms)
                for frag in ion2_msms : # find the number of matched peaks
                    frag_low = frag - mass_error
                    frag_high = frag + mass_error
                    frag_found = ion1_msms.between(frag_low, frag_high, inclusive = "both").sum()
                    if frag_found > 0 :
                        matched_peaks += 1
                    else :
                        total_peaks.append(frag)
                
                # Check number of matched peaks & matching score to validate frag
                if matched_peaks >= min_shared_peaks : # if number of matches above threshold
                    total_peaks = pd.Series(total_peaks)[total_peaks <= ion2_mz_high]
                    matching_score = round(matched_peaks / len(total_peaks),2)
                    if matching_score >= score_threshold:
                        edge_table.append((i, j, matched_peaks, len(total_peaks),
                                           matching_score,
                                           ion1_rt - node_table[rt_field][j],
                                           ion1_mz - node_table[mz_field][j]))
        
    # Results are stored in the edge table
    edge_table = pd.DataFrame(edge_table, columns = ['node_1', 'node_2',
                                                     'matched_peaks', 'total_peaks',
                                                     'matching_score', 'rt_gap',
                                                     'mz_gap'])

    return edge_table
