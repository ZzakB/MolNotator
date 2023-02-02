"""duplicate_finder.py - duplicate_finder module for the duplicate_filter"""
import sys
import pandas as pd
from matchms.similarity import ModifiedCosine
from pandas.core.common import flatten

def duplicate_finder(node_table, spectrum_list, params, ion_mode):
    """
    Takes a node_table and an spectrum_list from the duplicate_filter module to 
    find duplicate ions (close RT, m/z and cosine values). Returns a list of 
    duplicate ions to be deleted along with the representative ion that is to
    be kept.
    Parameters
    ----------
    node_table : pandas.DataFrame
        Dataframe containing each ion's metadata, including RT, m/z and the 
        relative position in the spectrum file.
    spectrum_list : list
        list of matchms.Spectrum objects, to be used to compare each ions with
        cosine similarity.
    params : dict
        Dictionary containing the global parameters for the process.
    Returns
    -------
    duplicate_table : pandas.DataFrame
        Dataframe containing the representative ions to be kept and the duplicate
        ions to be deleted.
    """
    
    # Get parameters & load cosine function
    mass_error = params['df_mass_error']
    rt_error = params['df_rt_error']
    cos_threshold = params['df_cos_threshold']
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        rt_error = rt_error/60

    # Start finding duplciates using the node_table index to identify ions.
    ions_idx = list(node_table.index)
    total_ions = len(ions_idx)
    duplicate_table = pd.DataFrame(columns = ['kept', 'dropped'])
    while len(ions_idx) > 1:
        perc = round((1-(len(ions_idx)/total_ions))*100,1)
        sys.stdout.write("\rFinding duplicates : {0}% ions processed".format(perc))
        sys.stdout.flush()
        tmp_duplicate_table = list()
        
        # Choose a seed ion to start finding duplicates.
        ion_pool = [ions_idx[0]]
        pool_size = 0
        pool_table = node_table.loc[ion_pool]
        
        # Find all ions with RT and MZ close to the seed ion
        while len(ion_pool) != pool_size:
            pool_size = len(ion_pool)
            min_rt = pool_table[f'{rt_field}'].min()
            max_rt = pool_table[f'{rt_field}'].max()
            min_mz = pool_table[f'{mz_field}'].min()
            max_mz = pool_table[f'{mz_field}'].max()
            pool_table = node_table[node_table[f'{rt_field}'].between(min_rt - rt_error,max_rt + rt_error)]
            pool_table = pool_table[pool_table[f'{mz_field}'].between(min_mz - mass_error, max_mz + mass_error)]
            ion_pool = list(pool_table.index)
            
        # Refresh ions_idx with the duplicates found (ion_pool)
        ions_idx = list(set(ions_idx) - set(ion_pool))
        
        # Start doing a cosine comparison between duplicates
        cos_table = list()
        while len(ion_pool) > 1 :
            ion_1 = ion_pool[0]
            specid_1 = node_table.loc[ion_1, 'spec_id']
            rt_1 = pool_table.loc[ion_1, f'{rt_field}']
            mz_1 = pool_table.loc[ion_1, f'{mz_field}']
            ion_pool.remove(ion_1)
            for ion_2 in ion_pool:
                specid_2 = node_table.loc[ion_2, 'spec_id']
                rt_2 = pool_table.loc[ion_2, f'{rt_field}']
                mz_2 = pool_table.loc[ion_2, f'{mz_field}']
                score, n_matches = modified_cosine.pair(spectrum_list[specid_1],
                                                        spectrum_list[specid_2])
                cos_table.append((ion_1, ion_2, score, n_matches, abs(mz_1 - mz_2),
                                  abs(rt_1 - rt_2)))
        cos_table = pd.DataFrame(cos_table, columns = ['ion_1', 'ion_2', 'cos',
                                                       'matched_peaks', 'd_mz', 'd_rt'])
        
        # Find duplicates based on the cos_threshold
        if len(cos_table) == 0 : 
            continue
        ion_pool = list(pool_table.index)
        while len(ion_pool) > 0 :
            ion_seed = pool_table.loc[ion_pool, 'TIC'].idxmax()
            seed_rt = pool_table.loc[ion_seed, f"{rt_field}"]
            seed_mz = pool_table.loc[ion_seed, f"{mz_field}"]

            candidates = cos_table.index[cos_table['ion_1'] == ion_seed].tolist() + cos_table.index[cos_table['ion_2'] == ion_seed].tolist()
            candidates = cos_table.loc[candidates].index[cos_table.loc[candidates, 'cos'] >= cos_threshold] 
            candidates = cos_table.loc[candidates, 'ion_1'].tolist() + cos_table.loc[candidates, 'ion_2'].tolist()
            candidates = list(set(candidates))
            if len(candidates) == 0 : 
                tmp_duplicate_table.append((ion_seed, []))
                ion_pool.remove(ion_seed)
                continue
            candidates.remove(ion_seed)
            candidates.sort()
            candidates = pool_table.loc[candidates]
            
            # Filter again based on RT and MZ threshold around the seed ion
            candidates = candidates[candidates[f'{rt_field}'].between(seed_rt - rt_error, seed_rt + rt_error, inclusive = "both")]
            candidates = candidates[candidates[f'{mz_field}'].between(seed_mz - mass_error, seed_mz + mass_error, inclusive = "both")]
            candidates = candidates.index.tolist()
            ion_pool.remove(ion_seed)
            ion_pool = list(set(ion_pool) - set(candidates))
            
            # Add duplicates to the duplicate table
            tmp_duplicate_table.append((ion_seed, candidates))
        tmp_duplicate_table = pd.DataFrame(tmp_duplicate_table, columns = ['kept', 'dropped'])

        # Ion affinity : find to which group certain ions belong 
        # based on cosine similarity (affinity)
        contested_ions = list(set(flatten(tmp_duplicate_table['dropped'])))
        contested_ions.sort()
        conflict_table = list()
        for i in contested_ions :
            conflicts = list()
            for j in tmp_duplicate_table.index:
                if i in tmp_duplicate_table.loc[j, "dropped"]:
                    conflicts.append(j)
            if len(conflicts) > 1 :
                spec_id_i = node_table.loc[i, "spec_id"]
                tmp_scores = list()
                for j in conflicts:
                    spec_id_j = node_table.loc[tmp_duplicate_table.loc[j, "kept"], "spec_id"]
                    score, n_matches = modified_cosine.pair(spectrum_list[spec_id_i],
                                                            spectrum_list[spec_id_j])
                    tmp_scores.append(score*n_matches)
                selected_ions = tmp_scores.index(max(tmp_scores))
                selected_ions = conflicts[selected_ions]
                conflict_table.append((i, conflicts, selected_ions))
        conflict_table = pd.DataFrame(conflict_table, columns = ['ion', 'conflicts', 'selected'])
        
        # Use the conflict table to reasign ions to their groupds
        for i in conflict_table.index:
            for j in conflict_table.loc[i, "conflicts"]:
                if j == conflict_table.loc[i, "selected"] :
                    continue
                tmp_duplicate_table.loc[j, "dropped"].remove(conflict_table.loc[i, "ion"])
        duplicate_table = duplicate_table.append(tmp_duplicate_table, ignore_index = True)
    
    dup_size = list()
    for i in duplicate_table.index:
        dup_size.append(len(duplicate_table.loc[i, "dropped"]))
    duplicate_table['dup_size'] = dup_size
    duplicate_table = duplicate_table[duplicate_table['dup_size'] > 0]
    duplicate_table.drop('dup_size', axis = 1, inplace = True)
    
    return duplicate_table

if __name__ == '__main__':
    duplicate_finder()
