"""adnotator_functions.py - module containing functions for the adnotator module"""
import os
import sys
import pandas as pd
from pandas.core.common import flatten
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosine
from matchms.filtering import default_filters
from MolNotator.others.species_validators import *
from tqdm import tqdm
import numpy as np

# Global functions

def Spectrum_processing(s):
    s = default_filters(s)
    return s

# Refresh edge table nodes with new node IDs
def refresh_edge_table(node_table, edge_table, old_id):
    print('Resetting edge table IDs...')
    for i in tqdm(edge_table.index):
        ion_mode = edge_table.loc[i, "ion_mode"]
        if ion_mode == "NEG" : opposed_mode = "POS"
        else : opposed_mode = "NEG"
        
        # Get node table of same mode as edge table
        tmp_node_table = node_table[node_table['ion_mode'] != opposed_mode]

        # Get old IDs, replace by new
        node_1 = edge_table.loc[i, "node_1"]
        node_2 = edge_table.loc[i, "node_2"]
        new_node_1 = tmp_node_table.index[tmp_node_table[old_id] == node_1][0]
        new_node_2 = tmp_node_table.index[tmp_node_table[old_id] == node_2][0]
        edge_table.loc[i, "node_1"] = new_node_1
        edge_table.loc[i, "node_2"] = new_node_2
    return edge_table

# Cluster IDs
def get_cluster_ids(node_table, edge_table):
    node_pool = list(node_table.index)
    singletons = list(edge_table["node_1"][edge_table['status'] == "self_edge"])
    node_pool = list(set(node_pool) - set(singletons))
    cluster_list = []
    cluster_size_list = []
    total_nodes = len(node_pool)
    while len(node_pool) > 0:
        new_cluster = [node_pool[0]]
        cluster_size = 0
        perc = round((1-(len(node_pool)/total_nodes))*100,1)
        sys.stdout.write("\rDefining new clusters : {0}%".format(perc))
        sys.stdout.flush()
        while cluster_size != len(new_cluster):
            cluster_size = len(new_cluster)
            tmp_idx = []
            for i in new_cluster:
                tmp_idx += list(edge_table.index[edge_table['node_1'] == i])
                tmp_idx += list(edge_table.index[edge_table['node_2'] == i])
            new_cluster += list(edge_table.loc[tmp_idx, 'node_1'].astype(int))
            new_cluster += list(edge_table.loc[tmp_idx, 'node_2'].astype(int))
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
    
    # Identify molecular clusters
    cluster_molecular = list()
    for i in cluster_table.index:
        node_list = cluster_table.loc[i, "cluster"].split('|')
        node_list = list(map(int, node_list))
        tmp_table_1 = node_table.loc[node_list]
        if sum(tmp_table_1['status_universal'] == "neutral") > 0 :
            cluster_molecular.append(True)
        else:
            cluster_molecular.append(False)
    cluster_table["molecular_cluster"] = cluster_molecular
    
    # Adding new cluster indexes to the merged node table
    node_table['cluster_id'] = [-1]*len(node_table)
    print('Assigning new cluster indexes...')
    for i in tqdm(cluster_table.index):
        node_list = list(map(int, cluster_table.loc[i, 'cluster'].split('|')))
        for j in node_list :
            node_table.loc[j, 'cluster_id'] = i
    return node_table

def get_samples(node_table, edge_table, csv_table, ion_mode):
    
    # Filter ion_mode and status
    ion_ids = node_table[node_table['ion_mode'] == ion_mode]
    ion_ids = ion_ids.index[ion_ids['status_universal'] != "neutral"]
    for col in csv_table.columns : node_table[col] = 0.0
    
    # Add samples to ions
    print(f"Adding samples to {ion_mode} ions...")
    for ion in tqdm(ion_ids):
        node_table.loc[ion, csv_table.columns] = csv_table.loc[ion, csv_table.columns]

    # Add samples to the samples column (neutrals)
    neutral_ids = node_table.index[node_table['status_universal'] == "neutral"].tolist()
    print(f"Adding samples to {ion_mode} neutrals...")
    for neutral in tqdm(neutral_ids):
        tmp_ions = list(edge_table['node_2'][edge_table['node_1'] == neutral])
        for col in csv_table.columns:
            node_table.loc[neutral, col] += node_table.loc[tmp_ions, col].sum()
    
    return node_table
 

################################################################# Ion functions

def ion_mass_calculator(mol_mass : float, adduct_mass : float, mol_count : int,
                        ion_charge : int) :
    """
    Calculates the m/z value of an ion species, given its supposed molecular
    mass (neutral) and other parameters from the adducts table.

    Parameters
    ----------
    mol_mass : float
        Molecule exact mass.
    adduct_mass : float
        Adduct m/z value as displayed in the adduct table.
    mol_count : int
        M count (M, 2M, 3M...) as in the ion formula, as displayed in the adduct
        table.
    ion_charge : int
        Ion charge for the considered adduct, from the adduct table.

    Returns
    -------
    adduct_mass
        Adduct mass for the considered molecule given a set adduct.
    """
    
    return round((mol_count*mol_mass + adduct_mass)/abs(ion_charge), 4)

def neutral_mass_calculator(ion_mz : float, adduct_mass : float,
                            mol_count : int, ion_charge : int) :
    """
    Calculates the molecular mass of a molecule given an ion species, its 
    m/z value and other parameters from the adducts table.

    Parameters
    ----------
    ion_mz : float
        Ion m/z value.
    adduct_mass : float
        Adduct m/z value as displayed in the adduct table.
    mol_count : int
        M count (M, 2M, 3M...) as in the ion formula, as displayed in the adduct
        table.
    ion_charge : int
        Ion charge for the considered adduct, from the adduct table.

    Returns
    -------
    neutral_mass
        Neutral mass for the considered ion given a set adduct.

    """
    return round(((ion_mz*abs(ion_charge)) - adduct_mass)/mol_count, 4)


def get_ion_ids(spectrum_file : str, params : dict):
    """
    Produce a list containing the unique ion identifiers in the spectrum
    file

    Parameters
    ----------
    spectrum_file : str
        Full file name for the spectrum file.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    ion_ids : list
        List containing the unique identifiers for each ion.

    """
    # Load the spectrum file
    spectrum_file = list(load_from_mgf(spectrum_file))

    # Retrieve the ion IDs
    ion_ids = []
    for spectrum in spectrum_file:
        ion_ids.append(int(spectrum.get(params['index_col'])))

    return ion_ids

def cross_sample_tables(spectrum_file : str, params : dict):
    """
    Create dataframes to report results using the complete original MGF file.
    Each dataframe will contain data for ions (indexes) for each given sample 
    (columns).

    Parameters
    ----------
    spectrum_file : str
        Full file name for the spectrum file.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    cross_annotations : pandas.DataFrame
        Dataframe to be filled with the annotations for each ion in each sample.
    cross_points : pandas.DataFrame
        Dataframe to be filled with the points for each annotation for each ion
        in each sample.
    cross_courts : pandas.DataFrame
        Dataframe to be filled with the court number for each ion in each sample.
    cross_houses : pandas.DataFrame
        Dataframe to be filled with the house number for each ion in each sample.
    cross_rules : pandas.DataFrame
        Dataframe to be filled with the rule points for each ion in each sample.
    cross_neutrals : pandas.DataFrame
        Dataframe to be filled with the neutral mass for each ion in each sample.
    ion_ids : list
        Series containing the unique identifiers for each ion

    """
    # Get ion IDs 
    ion_ids = get_ion_ids(spectrum_file, params)

    # Initialise the result tables
    cross_annotations = pd.DataFrame(index = ion_ids, dtype = str)
    cross_points = pd.DataFrame(index = ion_ids, dtype = float)
    cross_courts = pd.DataFrame(index = ion_ids, dtype = int)
    cross_houses = pd.DataFrame(index = ion_ids, dtype = int)
    cross_rules = pd.DataFrame(index = ion_ids, dtype = int)
    cross_neutrals = pd.DataFrame(index = ion_ids, dtype = int)

    return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals, ion_ids


def neutral_tabler(ion_mz : float, adduct_table_primary):
    """
    Computes all possible molecular masses (hypothetical neutrals) for an ion
    (ion 1) given the different ion species available in the adducts table.

    Parameters
    ----------
    ion_mz : float
        m/z value for the considered ion
    adduct_table_primary : pandas.DataFrame
        Dataframe imported from the primary adduct table.

    Returns
    -------
    neutral_table : pandas.DataFrame
        Dataframe with the neutral masses added for the input m/z value.

    """

    neutral_table = adduct_table_primary.copy()
    neutral_table['neutral_mass'] = [0.0]*len(neutral_table)
    for i in neutral_table.index :
        neutral_table.loc[i, 'neutral_mass'] = neutral_mass_calculator(ion_mz, 
                         neutral_table['Adduct_mass'][i], 
                         neutral_table['Mol_multiplier'][i], 
                         neutral_table['Charge'][i])
    return neutral_table


def ion_hypotheses_tabler(neutral_table):
    """
    Computes all possible ion m/z values (ion 2) given all possible molecular masses
    computed before with the neutral_tabler function (hypothetical neutrals), 
    for a single ion (ion 1).

    Parameters
    ----------
    neutral_table : pandas.DataFrame
        Neutral table produced by the neutral_tabler function

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Dataframe containing the adduct m/z values for all molecules in the neutral
        table.

    """
    
    # Initialise ion hypothesis table
    ion_hypotheses_table = list()

    # Search for each adduct combination in the neutral table
    for i in neutral_table.index :
        for j in neutral_table.index :
            if i == j: continue
            tmp_mz = ion_mass_calculator(neutral_table.loc[i, 'neutral_mass'],
                                         neutral_table.loc[j, 'Adduct_mass'],
                                         neutral_table.loc[j, 'Mol_multiplier'],
                                         neutral_table.loc[j, 'Charge'])
            ion_hypotheses_table.append((i, j, neutral_table.loc[j, 'Complexity'],
                                         tmp_mz))
    return pd.DataFrame(ion_hypotheses_table, columns = ['Ion1_adduct', 'Ion2_adduct',
                                                         'Ion2_complexity', 'Ion2_mz'])

def point_counter(ion_hypotheses_table, coelution_table, params : dict):
    """
    Checks among the coeluted ions (in the coelution table) if any of them
    matches an "ion 2" from the ion_hypotheses_table. Eliminates all ion 2 
    hypotheses that could not be matched to any coeluted ion.

    Parameters
    ----------
    ion_hypotheses_table : pandas.DataFrame
        Dataframe containing the adduct m/z values for all molecules in the
        neutral table.
    coelution_table : pandas.DataFrame
        Dataframe containing coeluted ions.
    mass_error : float
        Mass error in the parameter file.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        ion_hypotheses_table with scored surviving hypotheses.

    """
    
    hits = list()
    ion_indexes = list()
    for i in ion_hypotheses_table.index :
        ion2_mz = ion_hypotheses_table.loc[i, 'Ion2_mz']
        ion2_mz_low = ion2_mz - params['an_mass_error']
        ion2_mz_high = ion2_mz + params['an_mass_error']
        hit_count = coelution_table[params['mz_field']].between(ion2_mz_low,
                                        ion2_mz_high, inclusive = "both").sum()
        if hit_count > 0 :
            hit_ids = list(coelution_table.index[coelution_table[params['mz_field']].between(ion2_mz_low, 
                                        ion2_mz_high, inclusive = "both")])
        else:
            hit_ids = []
        hits.append(hit_count)
        ion_indexes.append(hit_ids)
    ion_hypotheses_table['hit_count'] = hits
    ion_hypotheses_table['hit_indexes'] = ion_indexes
    ion_hypotheses_table = ion_hypotheses_table[ion_hypotheses_table['hit_count'] > 0]
    ion_hypotheses_table.reset_index(drop = True, inplace = True)
    return ion_hypotheses_table

def cosine_validation(ion_1_idx : int, node_table, spectrum_list : list,
                      ion_hypotheses_table, adduct_table_primary, params : dict):
    """
    For each ionisation hypothesis (ion 1 and 2 pairs), reports the
    cosine score from the cosine_table. In case of multiple hits in the
    "hit_indexes" columns, the one with the best cosine score is kept (chosen
    at random if both have the same score : either both completely unrelated
    or the same spectrum). The duplicate hits are transferred to a duplicate_df
    for annotation in the final network, they will take no part in calculation.
    Points are then awared according to the cosine score:
    1+cosine score if the score is above threshold, 0 if it is below.

    Parameters
    ----------
    ion_1_idx : int
        Current ion ID / index.
    node_table : pandas.DataFrame
        Node table with the ion metadata.
    spectrum_list : list
        List of matchms.Spectrum objects.
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table obtained from other functions.
    adduct_table_primary : pandas.DataFrame
        Primary adduct table from the parameter folder.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table with points awarded.
    duplicate_df : pandas.DataFrame
        Duplicates found in the ion hypotheses table.

    """
    
    # Load parameters
    modified_cosine = ModifiedCosine(tolerance=params['an_mass_error'])
    cosine_threshold = params['an_cos_threshold']
    spec_id_1 = node_table.loc[ion_1_idx, 'spec_id']
    
    # Initialise variables
    score_table = list()

    # duplicate_df will contain adducts with roles considered redundant
    duplicate_df = pd.DataFrame(columns = ["ion_1_idx", "ion_1_adduct",
                                           "ion_2_idx", "ion_2_adduct",
                                           "selected_ion", "cosine_score",
                                           "matched_peaks", "different_group"])
    for i in ion_hypotheses_table.index:
        
        # Check adduct groups for ions 1 and 2
        group_1 = adduct_table_primary.loc[ion_hypotheses_table.loc[i, "Ion1_adduct"], "Group"]
        group_2 = adduct_table_primary.loc[ion_hypotheses_table.loc[i, "Ion2_adduct"], "Group"]
        if group_1 == group_2: different_group = False
        else: different_group = True
        
        # If only one ion 2, calculate scores
        if ion_hypotheses_table.loc[i, "hit_count"] == 1 :
            ion_2_idx = ion_hypotheses_table.loc[i,"hit_indexes"][0]
            spec_id_2 = node_table.loc[ion_2_idx, "spec_id"]
            score, n_matches = modified_cosine.pair(spectrum_list[spec_id_1],
                                                    spectrum_list[spec_id_2])
            prod = score * n_matches
        # If several ion 2, select the best one based on scores
        else:
            selected_hit = pd.DataFrame(columns = ['cos', 'peaks', 'prod'])
            for hit in ion_hypotheses_table.loc[i,"hit_indexes"]:
                spec_id_2 = node_table.loc[hit, "spec_id"]
                score, n_matches = modified_cosine.pair(spectrum_list[spec_id_1],
                                                        spectrum_list[spec_id_2])
                selected_hit.loc[hit] = [score, n_matches, score * n_matches]

            selected_hit.sort_values('prod', ascending = False, inplace =True)
            new_hit = selected_hit.index[0]
            score = selected_hit['cos'].iloc[0]
            n_matches = int(selected_hit['peaks'].iloc[0])
            prod = selected_hit['prod'].iloc[0]
            selected_hit.drop(new_hit, inplace = True)
            ion_hypotheses_table.loc[i, "hit_indexes"] = [[new_hit]]
            for j in selected_hit.index:
                tmp_idx = len(duplicate_df)
                duplicate_df.loc[tmp_idx] = [ion_1_idx,
                                    ion_hypotheses_table.loc[i, "Ion1_adduct"],
                                 j, ion_hypotheses_table.loc[i, "Ion2_adduct"],
                                 new_hit, selected_hit.loc[j, 'cos'], selected_hit.loc[j, 'peaks'],
                                 different_group]
        # If score is below threshold, nullify all.
        if score < cosine_threshold: 
            score = 0.0
            prod = 0.0
        score_table.append((different_group, score, n_matches, prod))
    
    # Convert scores to dataframe and merge with ion_hypothesis thable
    score_table = pd.DataFrame(score_table, columns = ['different_group',
                                                       'cosine_score',
                                                       'matched_peaks',
                                                       'product'])
    ion_hypotheses_table = pd.concat([ion_hypotheses_table, score_table], axis = 1)
    
    # Convert the hit_indexes col (which are now guaranteed to be only one), to int
    ion_hypotheses_table['hit_indexes'] = ion_hypotheses_table['hit_indexes'].str[0]
    
    # Filter by score and same group criterium
    tmp_bool = ion_hypotheses_table["different_group"] + ion_hypotheses_table["cosine_score"] >= cosine_threshold
    duplicate_bool = duplicate_df["different_group"] + duplicate_df["cosine_score"] >= cosine_threshold
    return ion_hypotheses_table[tmp_bool], duplicate_df[duplicate_bool]

def species_rules(ion_1_spec_id : int, ion_hypotheses_table, adduct_table_primary, node_table,
                  spectrum_list : list, params : dict, ion_mode : str):
    """
    Awards points to hypotheses for which the annotations for ions 1 or 2 
    could be validated by the ion species rules (checks if the annotation can 
    be confirmed).

    Parameters
    ----------
    ion_1_spec_id : int
        ion 1 relative position in the sample's spectrum file.
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table obtained from other functions.
    adduct_table_primary : pandas.DataFrame
        Primary adduct table from the parameter folder.
    node_table : pandas.DataFrame
        Node table with the ion metadata.
    spectrum_list : list
        List of matchms.Spectrum objects.
    params : dict
        Dictionary containing the global parameters for the process.
    ion_mode : str
        POS or NEG, to be chosen for the species rules.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table with species rule points.

    """

    # Load parameters and prepare the rule_points column
    prec_mass_error = params['an_prec_mass_error']
    ion_hypotheses_table['rule_points'] = [0.0]*len(ion_hypotheses_table)

    # Check each hypothesis
    for hypothesis in ion_hypotheses_table.index:
        
        # Get the ion 1 adduct
        ion_1_adduct = ion_hypotheses_table.loc[hypothesis, "Ion1_adduct"]
        ion_1_adduct = adduct_table_primary.loc[ion_1_adduct, "Adduct_code"]

        # Select the species rule based on the adduct and the ion mode
        Species_rule = Validator_choice(ion_1_adduct, ion_mode)

        # Increment the score to the ion hypotheses table
        ion_hypotheses_table.loc[hypothesis, "rule_points"] += Species_rule(prec_mass_error,
                                                                            ion_1_spec_id,
                                                                            spectrum_list)
        
        # Do the same for ion 2 
        ion_2_adduct = ion_hypotheses_table.loc[hypothesis, "Ion2_adduct"]
        ion_2_adduct = adduct_table_primary.loc[ion_2_adduct, "Adduct_code"]
        ion_2_spec_id = ion_hypotheses_table.loc[hypothesis, "hit_indexes"]
        ion_2_spec_id = node_table.loc[ion_2_spec_id, "spec_id"]
        Species_rule = Validator_choice(ion_2_adduct, ion_mode)
        ion_hypotheses_table.loc[hypothesis, "rule_points"] += Species_rule(prec_mass_error,
                                                                            ion_2_spec_id,
                                                                            spectrum_list)
    return ion_hypotheses_table

def complex_points(neutral_table, ion_hypotheses_table):
    """For a hypothesis bearing an ion 2 with a complex form (neutral complexes,
    i.e. the third componend of the adduct code : Mol|ion|complex), points are
    awared if the uncomplexed form is found among the other hypotheses.

    Parameters
    ----------
    neutral_table : pandas.DataFrame
        Neutral table provided by other MolNotator functions.
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table provided by other MolNotator functions.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypothesis table updated with solvent complex points. 

    """

    ion_hypotheses_table['Complex_points'] = [0]*len(ion_hypotheses_table)
    for i in ion_hypotheses_table.index :
        
        # Get ion 2 adduct index
        ion1_idx = ion_hypotheses_table.loc[i, 'Ion1_adduct']
        
        # Get ion 2 adduct code
        ion2_code = neutral_table.loc[ion_hypotheses_table.loc[i,
                                            'Ion2_adduct'], 'Adduct_code']
        # Check if solvent complex
        ion2_code_split = ion2_code.split('|')
        if ion2_code_split[2] != "" :
            
            # Get the index of the uncomplexed form in the adduct table
            ion2_code_split[2] = ""
            ion2_uncomplexed = '|'.join(ion2_code_split)
            ion2_uncomplexed_idx = neutral_table.index[neutral_table['Adduct_code']==ion2_uncomplexed][0]
            
            # Check if the uncomplexed form was detected using the level_2 table filtered
            bool_1 = ion_hypotheses_table['Ion1_adduct'] == ion1_idx
            bool_2 = ion_hypotheses_table['Ion2_adduct'] == ion2_uncomplexed_idx
            
            # Reward with points level 2 ion hypotheses for which uncomplexed 
            # ions could be detected
            if (bool_1 & bool_2).sum() > 0 : 
                ion_hypotheses_table.loc[i, 'Complex_points'] += 1
    return ion_hypotheses_table

def fragnotator_points(ion1_idx : int, ion_hypotheses_table, edge_table,
                       status : str):
    """
    Awards points for each hypothesis in which the ions 1 and 2 are already
    linked by a fragmentation bond by fragnotator.

    Parameters
    ----------
    ion1_idx : int
        Current ion ID.
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table provided by other MolNotator functions.
    edge_table : pandas.DataFrame
        Edge table provided by MolNotator
    status : str
        Status of the current ion, either fragment or precursor.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table updated with fragnotator points.
    """
    
    # Initialise results lit
    ion_hypotheses_table['fragnotator_points'] = [0]*len(ion_hypotheses_table)
    
    # Set columns to be checked
    if status != "fragment":
        ion_col_1 = 'node_1'
        ion_col_2 = 'node_2'
    else:
        ion_col_1 = 'node_2'
        ion_col_2 = 'node_1'
    
    # Award points
    temp_edges = edge_table[edge_table[ion_col_1] == ion1_idx]
    if len(temp_edges) == 0 : return ion_hypotheses_table
    for j in ion_hypotheses_table.index:
        ion2 = ion_hypotheses_table.loc[j, 'hit_indexes']
        if ion2 in list(temp_edges[ion_col_2]):
            ion_hypotheses_table.loc[j, 'fragnotator_points'] = 1

    return ion_hypotheses_table

def get_adduct_code(ion_hypotheses_table, neutral_table):
    """
    Adds the adduct code for each hypothesis' ion 1.

    Parameters
    ----------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table provided by other MolNotator functions.
    neutral_table : pandas.DataFrame
        Neutral table provided by other MolNotator functions.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table updated with adduct codes.

    """
    adduct_col = list()
    for i in ion_hypotheses_table.index :
        j = ion_hypotheses_table.loc[i, 'Ion1_adduct']
        adduct_col.append(neutral_table.loc[j, 'Adduct_code'])
    ion_hypotheses_table['adduct_code'] = adduct_col
    return ion_hypotheses_table

def weighted_points(ion_hypotheses_table, adduct_table_primary):
    """
    Calculated the points accumulated by each hypothesis, considering cosine 
    similarity, rule points, complex points, fragnotator points and the ion 2
    complexity. 

    Parameters
    ----------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table provided by other MolNotator functions.
    adduct_table_primary : pandas.DataFrame
        Primary adduct table from the parameter folder.

    Returns
    -------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table updated with final scores.

    """
    
    # Add column and sum all points to it
    ion_hypotheses_table['weighted_points'] = [0.0]*len(ion_hypotheses_table)
    
    # List columns with points to be summed
    summed_cols = ['hit_count', 'cosine_score', 'rule_points', 'Complex_points',
              'fragnotator_points']
    
    for i in ion_hypotheses_table.index :
        
        # Get complexities
        adduct = ion_hypotheses_table.loc[i, 'Ion1_adduct']
        c1 = adduct_table_primary.loc[adduct, "Complexity"]
        c2 = ion_hypotheses_table.loc[i, 'Ion2_complexity']
        
        # Calculate score
        score = sum(ion_hypotheses_table.loc[i, summed_cols])
        score = score / (c1 + c2)
        
        # Report points
        ion_hypotheses_table.loc[i, 'weighted_points'] = score
        
    return ion_hypotheses_table

def merged_adducts_table(ion_hypotheses_table, ion_1_idx):
    """
    Merges the ion_hypotheses_table for the current ion (i) being processed
    to all previously processed ions from the sample.

    Parameters
    ----------
    ion_hypotheses_table : pandas.DataFrame
        Ion hypotheses table provided by other MolNotator functions.
    ion_1_idx : int
        Current ion ID / index.

    Returns
    -------
    merged_table_local : pandas.DataFrame
        Summary table for all hypotheses for the considered ion (ion_1_idx).

    """

    # if no hypotheses, skip
    if len(ion_hypotheses_table) == 0 : return ion_hypotheses_table

    # List all unique adducts (unique hypotheses for the current ion)
    unique_adducts = list(ion_hypotheses_table['adduct_code'].unique())

    # Initialise the merged table
    merged_table_local = pd.DataFrame(index = range(len(unique_adducts)))
    merged_table_local['ion_id'] = [ion_1_idx]*len(unique_adducts)
    merged_table_local['adduct'] = unique_adducts
    merged_table_local['relation_codes'] = ['']*len(merged_table_local)
    merged_table_local['hit_count'] = [0]*len(merged_table_local)
    merged_table_local['hit_indexes'] = ['']*len(merged_table_local)
    merged_table_local['weighted_points'] = [0.0]*len(merged_table_local)


    # Fill the table with results
    for i in merged_table_local.index :
        adduct = merged_table_local.loc[i, 'adduct']
        tmp_table_1 = ion_hypotheses_table[ion_hypotheses_table['adduct_code'] == adduct]
        ion1_adduct = tmp_table_1['Ion1_adduct'].astype(str)
        ion2_adduct = tmp_table_1['Ion2_adduct'].astype(str)
        merged_table_local.loc[i, 'relation_codes'] = '|'.join(ion1_adduct+":"+ion2_adduct)
        merged_table_local.loc[i, 'hit_count'] = tmp_table_1['hit_count'].sum()
        merged_table_local.loc[i, 'hit_indexes'] = '|'.join(tmp_table_1['hit_indexes'].astype(str))
        merged_table_local.loc[i, 'weighted_points'] = tmp_table_1['weighted_points'].sum()
    return merged_table_local

########################################################### Reporting functions

def cohort_tabler(neutral_table, merged_table, node_table):
    """
    Produce a cohort_table, with the ion_indexes as rows (ion IDs) and the
    ionisation hypotheses indexes as columns (taken from the merged_table).
    Each ion (row) will have 0, 1 or several ionisation hypotheses (columns). 
    and each of these columns will affect different ions based on the number 
    of ion_2 ionisation hypotheses that were found. Related ionisation hypotheses
    (ion ID + adduct formula) that are not mutually exclusive (i.e. different 
    ionisations for a single ion) are named "Cohorts".

    Parameters
    ----------
    neutral_table : pandas.DataFrame
        Neutral table provided by other MolNotator functions.
    merged_table : pandas.DataFrame
        Summary table for all hypotheses for all ions in the considered sample.
    node_table : pandas.DataFrame
        Node table with the ion metadata.

    Returns
    -------
    cohort_table : pandas.DataFrame
        Dataframe containing all adduct annotations for all ions in the current
        sample. Index = ion IDs, columns = cohorts (several hypothetical adduct
        annotations).

    """
    
    # Initialise the cohort table
    ion_indexes = list(node_table.index)
    cohort_table = pd.DataFrame(data = None, index = ion_indexes,
                                columns = merged_table.index)
    
    print('Producing cohort table...')
    for i in tqdm(merged_table.index) :
        current_ion = int(merged_table.loc[i, 'ion_id'])
        
        # Get relation codes and ions 
        split_relations_code = merged_table.loc[i,'relation_codes'].split('|')
        split_hit_ids = merged_table.loc[i,'hit_indexes'].split('|')
        
        # Get ion 1 adduct
        ion1_adduct = int(split_relations_code[0].split(':')[0])
        ion1_adduct = neutral_table.loc[ion1_adduct, 'Adduct_code']
        
        # Report the adduct in the cohort table
        cohort_table.loc[current_ion, i] = ion1_adduct
        
        # Replace codes by str representation of the adducts
        for j in range(len(split_relations_code)) :
            split_relations_code[j] = split_relations_code[j].split(':')[1]
            split_relations_code[j] = neutral_table.loc[int(split_relations_code[j]),
                                'Adduct_code']
            
        # Get the hit ion ID and report in the cohort table the adduct annotation
        for j in range(len(split_hit_ids)) :
            split_ids = split_hit_ids[j]
            cohort_table.loc[int(split_ids), i] = split_relations_code[j]
            
    return cohort_table

def supercohorts_tabler(cohort_table):
    """
    Merges compatible cohorts into a single supercohort to reduce combinations
    and speed up calculation times. Returns the supercohorts_table and the
    transitions_table which keeps track of the cohorts forming each supercohort.

    Parameters
    ----------
    cohort_table : pandas.DataFrame
        Dataframe containing all adduct annotations for all ions in the current
        sample. Index = ion IDs, columns = cohorts (several hypothetical adduct
        annotations).

    Returns
    -------
    supercohorts_table : pandas.DataFrame
        Dataframe displaying the supercohorts, i.e. mergeable cohorts that could
        be simplified into a single one.
    transition_table : pandas.DataFrame
        Dataframe connecting supercohorts and cohorts by their IDs.

    """

    # Create a transitions table to connect cohorts to supercohorts
    transition_table = pd.DataFrame(columns = ["supercohort", "cohorts"])

    # Initialisethe supercohorts_table
    supercohorts_table = pd.DataFrame(index = cohort_table.index)

    # cohort_pool contains all cohorts to be picked
    cohort_pool = pd.Series(cohort_table.columns, index = cohort_table.columns)
    while len(cohort_pool) > 0 :
        
        # Pick one cohort and start merging it with others to create the supercohort
        i = cohort_pool.iloc[0]
        supercohort = cohort_table[i].dropna()
        mergeable_cohorts = [i]
        all_ions = list(supercohort.index)
        added_ions = all_ions
        len_0 = 0
        while len(all_ions) != len_0:
            len_0 = len(all_ions)
            for j in added_ions:
                annotation = supercohort[j]
                other_cohorts = cohort_table.loc[j].dropna()
                other_cohorts = other_cohorts[list(set(other_cohorts.index) - set(mergeable_cohorts))]
                other_cohorts = other_cohorts.index[other_cohorts == annotation]
                mergeable_cohorts += list(other_cohorts)
                added_ions = []
                for k in other_cohorts :
                    new_ions = list(set(cohort_table[k].dropna().index) - set(all_ions))
                    added_ions += new_ions
                    for l in new_ions:
                        supercohort[l] = cohort_table.loc[l, k]
                all_ions += added_ions
        
        # If no mergeable cohorts, report directly and end
        if len(mergeable_cohorts) == 0 :
            supercohorts_table[len(supercohorts_table.columns)] = cohort_table[i]
            transition_table.loc[len(transition_table)] = [len(transition_table), [i]]
            cohort_pool.drop(i, inplace = True)
            continue
        
        # Otherwise, clean the mergeable cohorts list
        mergeable_cohorts = list(set(mergeable_cohorts))
        mergeable_cohorts.sort()
        
        # Add the new supercohort created as a column
        supercohorts_table[len(supercohorts_table.columns)] = [None]*len(supercohorts_table)
        
        # Fill with the annotations available from the mergeable cohorts
        for j in supercohort.index:
            supercohorts_table.loc[j, len(supercohorts_table.columns) -1] = supercohort[j]
            
        # Report the connexions on the transition table
        transition_table.loc[len(transition_table)] = [len(transition_table), mergeable_cohorts]
        
        # Remove the merged cohorts from the cohort pool
        cohort_pool.drop(mergeable_cohorts, inplace = True)
        
    return supercohorts_table, transition_table

def get_cohort_ions(cohorts : {list, int}, cohort_table) :
    """
    Returns every ion ID associated to the submitted cohorts list.

    Parameters
    ----------
    cohorts : {list, int}
        List of cohorts to be investigated.
    cohort_table : pandas.DataFrame
        Table containing the cohort-ion pairs.

    Returns
    -------
    ion_list : list
        List of all ions involved with the supplied cohort(s).

    """
    
    # If cohort is int, convert to list
    if type(cohorts) != list:
        cohorts = [cohorts]

    # Get all ions involved in each cohort
    ion_list = list()
    for i in cohorts :
        ion_list += list(cohort_table[i].dropna().index)

    # Clean up and export
    ion_list = list(set(ion_list))
    ion_list.sort()
    return ion_list

def get_ion_cohorts(ion_list : {list, int}, cohort_table):
    """
    Returns every cohorts associated with the submitted ions list

    Parameters
    ----------
    ion_list : {list, int}
        List of ions to be investigated.
    cohort_table : pandas.DataFrame
        Table containing the cohort-ion pairs.

    Returns
    -------
    cohort_list : list
        List of all cohorts involved with the supplied ion(s).

    """
    # If ion is int, convert to list
    if type(ion_list) != list:
        ion_list = [ion_list]
        
    # Get all cohorts involved with each ion
    cohort_list = list()
    for i in ion_list :
        cohort_list += list(cohort_table.loc[i].dropna().index)
    
    # Clean up and export
    cohort_list = list(set(cohort_list))
    cohort_list.sort()
    return cohort_list

def get_court(cohort_list : list, cohort_table):
    """
    Produces a court from one or several cohorts. A Court is a collection
    of related ion IDs and cohorts. Two cohorts, each with their own ions,
    sharing at least one ion, will be considered related.

    Parameters
    ----------
    cohort_list : list
        List of cohorts to be investigated.
    cohort_table : pandas.DataFrame
        Dataframe displaying the cohorts or supercohorts.

    Returns
    -------
    ion_list : list
        All ions from the same court.
    cohort_list : list
        All cohorts from the same court.

    """
    
    # Start by saving default ion and cohort count for a new court
    nb_ions = 0
    nb_cohorts = 1
    
    # Boolean - while there are ions and cohorts being added to court : continue
    court_expanding = True
    while court_expanding :
        
        # Get all ions that can be connected to the cohort list
        ion_list = get_cohort_ions(cohort_list, cohort_table)
        
        # Get all cohorts that can be connected to the updated ion list
        cohort_list = get_ion_cohorts(ion_list, cohort_table)

        # If no ions or cohorts are added, stop.
        if nb_ions + nb_cohorts == len(ion_list) + len(cohort_list) :
            court_expanding = False

        # else, add the new elements to the global list and restart
        else :
            nb_ions = len(ion_list)
            nb_cohorts = len(cohort_list)
            
    return ion_list, cohort_list
    


def get_court_table(cohort_table) :
    """
    Produces a court_table from a cohort_table, containing all Courts produced
    from the different cohorts and the associated ions. A Court is a collection
    of related ion IDs and cohorts. Two cohorts, each with their own ions, 
    sharing at least one ion, will be considered related.

    Parameters
    ----------
    cohort_table : pandas.DataFrame
        Dataframe displaying the cohorts or supercohorts.

    Returns
    -------
    court_table : pandas.DataFrame
        Dataframe containing all courts produced by the cohorts.

    """

    # Create a cohort table copy to be emptied
    precourt_table = cohort_table.copy()
    
    # Start first court ID as 0
    court_id = 0 

    # Initialise the output court table
    court_table = pd.DataFrame(columns = ['Cohort_List', 'Ion_List'])

    # Empty precourt table until no indexes or columns remain
    while (len(precourt_table.index) > 0 and len(precourt_table.columns) > 0):

        # Get the first cohort to process
        i = precourt_table.columns[0]
        
        # Get all cohort related ions and store them in the court_table
        ion_list, cohort_list = get_court(i, precourt_table)
        court_table.loc[court_id] = [cohort_list, ion_list]
        
        # remove all involved ions and cohorts from the cohort_to_court_table
        precourt_table.drop(ion_list, axis = 0, inplace = True)
        precourt_table.drop(cohort_list, axis = 1, inplace = True)
        court_id += 1
    
    # Count number of ions and courts in the court table
    if len(court_table.index) > 0:
        court_table['Cohort_count'] = [0]*len(court_table)
        court_table['Ion_count'] = [0]*len(court_table)
        for i in court_table.index :
            court_table.loc[i, 'Cohort_count'] = len(court_table.loc[i, 'Cohort_List'])
            court_table.loc[i, 'Ion_count'] = len(court_table.loc[i, 'Ion_List'])
    
    # Sort based on the ion count and cohort count
    
    court_table.sort_values(["Cohort_count", "Ion_count"], ascending = False, inplace = True)
    
    # Reset index to clean the court IDs and export
    court_table.reset_index(drop = True, inplace = True)
    return court_table

def get_incompatible_cohorts(cohort : int, supercohorts_table):
    """
    Gets supercohorts incompatible with the given supercohort.

    Parameters
    ----------
    cohort : int
        Current cohort ID.
    supercohorts_table : pandas.DataFrame
        Dataframe containing ion IDs x Supercohorts IDs.

    Returns
    -------
    incompatibles : list
        List of cohorts incompatible with cohort.

    """
    # Get the cohort's ions
    cohort_ions = list(supercohorts_table[cohort].dropna().index)
    
    # Store in incompatibles the cohorts with annotations conflicting with the current one
    incompatibles = []
    for ion in cohort_ions:
        annotation = supercohorts_table.loc[ion, cohort]
        incompatible = supercohorts_table.loc[ion].dropna() != annotation
        incompatibles += list(incompatible.index[incompatible])
    return list(set(incompatibles))


def initialise_house_table(seed_cohort : int, court_cohorts : list, supercohorts_table):
    """
    Starts house building by combining compatible cohorts, using an initial
    cohort as a seed.

    Parameters
    ----------
    seed_cohort : int
        Current cohort ID.
    court_cohorts : list
        List of cohorts considered.
    supercohorts_table : pandas.DataFrame
        Dataframe containing ion IDs x Supercohorts IDs.

    Returns
    -------
    next_list : list
        List of compatible cohorts that can be added to each house.
    house_list : list
        List of growing houses.
    compatibles_list : list
        Remaining compatible cohorts for each house.
    left_list : list
        Count of remaining compatible cohorts for each house.

    """
    
    # Get cohorts incompatible with the seed cohort
    incompatibles = get_incompatible_cohorts(seed_cohort, supercohorts_table)

    # Deduce the compatible cohorts
    compatibles = list(set(court_cohorts) - set(incompatibles))

    # Initialise variables for the court expansion process
    next_list = compatibles # List of compatible cohorts that can be added to each house
    house_list = [[seed_cohort]]*len(next_list) # List of growing houses
    compatibles_list = [None]*len(next_list) # slots for the cohorts compatible with the future hourses
    left_list = [0]*len(next_list) # Count of remaining compatible cohorts for each house
    added = [] # List of compatible cohorts already added to houses

    # Initiate the house expansion process with one step
    for i in range(len(next_list)):
        
        # Add each compatible cohort to houses
        added.append(next_list[i])
        
        # Refresh the new incompatibles given that a new cohort is added to the house
        tmp_incompatibles = get_incompatible_cohorts(next_list[i], supercohorts_table)
        tmp_incompatibles += added 

        # Refresh the compatibles list givne the incompatible
        compatibles_list[i] = list(set(compatibles) - set(tmp_incompatibles))

        # Count the remaining compatible cohorts
        left_list[i] = len(compatibles_list[i])
    return next_list, house_list, compatibles_list, left_list

def clean_house_table(next_list : list, house_list : list, compatibles_list : list,
                      left_list : list, houses : list):
    """
    Removes houses that finished expanding (no compatible cohorts left) from
    the house table, and exports them to the houses list.

    Parameters
    ----------
    next_list : list
        List of compatible cohorts that can be added to each house.
    house_list : list
        List of growing houses.
    compatibles_list : list
        Remaining compatible cohorts for each house.
    left_list : list
        Count of remaining compatible cohorts for each house.
    houses : list
        Future house table.

    Returns
    -------
    next_list : list
        List of compatible cohorts that can be added to each house, cleaned.
    house_list : list
        List of growing houses, cleaned.
    compatibles_list : list
        Remaining compatible cohorts for each house, cleaned.
    left_list : list
        Count of remaining compatible cohorts for each house, cleaned.
    houses : list
        Future house table, updated by one step.

    """
    # export idx is a list of houses that finished expanding : no cohorts left 
    export_idx = [i for i in range(len(left_list)) if left_list[i] == 0]
    export_idx.sort(reverse = True)
    
    # Complete the house and add it to houses
    for i in export_idx :
        new_house = house_list[i] + [next_list[i]]
        new_house = list(set(new_house))
        new_house.sort()
        houses.append(new_house) 
        next_list.pop(i)
        house_list.pop(i)
        compatibles_list.pop(i)
        left_list.pop(i)
    return next_list, house_list, compatibles_list, left_list, houses

def expand_house_table(next_list : list, house_list : list, compatibles_list : list,
                       left_list : list, supercohorts_table):
    """
    Expands the houses in the house table by appending to them one of each
    possible compatible cohorts.

    Parameters
    ----------
    next_list : list
        List of compatible cohorts that can be added to each house.
    house_list : list
        List of growing houses.
    compatibles_list : list
        Remaining compatible cohorts for each house.
    left_list : list
        Count of remaining compatible cohorts for each house.
    supercohorts_table : pandas.DataFrame
        Dataframe containing ion IDs x Supercohorts IDs.

    Returns
    -------
    next_list : list
        List of compatible cohorts that can be added to each house, updated by one step.
    house_list : list
        List of growing houses, updated by one step.
    compatibles_list : list
        Remaining compatible cohorts for each house, updated by one step.
    left_list : list
        Count of remaining compatible cohorts for each house, updated by one step.

    """

    drop_list = []
    print("Expanding house table...")
    for i in tqdm(range(len(next_list))):
        
        # For remaining compatible cohorts, add them to a new house, 
        drop_list.append(i) # store old houses idx to remove later
        new_house = house_list[i] + [next_list[i]] # produce the new house
        added = []
        compatibles = compatibles_list[i]
        
        # Expand
        for j in compatibles:
            added.append(j)
            incompatibles = get_incompatible_cohorts(j, supercohorts_table)
            incompatibles += added
            next_list.append(j)
            house_list.append(new_house)
            compatibles_list.append(list(set(compatibles) - set(incompatibles)))
            left_list.append(len(list(set(compatibles) - set(incompatibles))))

    # Drop old houses
    drop_list.sort(reverse = True)
    for i in drop_list:
        next_list.pop(i)
        house_list.pop(i)
        compatibles_list.pop(i)
        left_list.pop(i)
    return next_list, house_list, compatibles_list, left_list

def house_selection(court_table, supercohorts_table, node_table, transition_table,
                    merged_table, spectrum_list : list, params : dict):
    """
    Finds the houses for each given court, i.e. supercohort combinations
    in a court that are compatible with each other. Once all combinations are 
    explored, the most likely house is each cohort is selected based on the points
    cumulated by their cohorts. If there is a draw between two or more houses, they
    are all kept. Final selection will happen when cross examining samples.

    Parameters
    ----------
    court_table : pandas.DataFrame
        Dataframe containing all courts produced by the cohorts.
    supercohorts_table : pandas.DataFrame
        Dataframe containing ion IDs x Supercohorts IDs.
    node_table : pandas.DataFrame
        Node table with the ion metadata.
    transition_table : pandas.DataFrame
        Dataframe connecting supercohorts and cohorts by their IDs.
    merged_table : pandas.DataFrame
        Summary table for all hypotheses for all ions in the considered sample.
    spectrum_list : list
        List of matchms.Spectrum objects.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    court_table : pandas.DataFrame
        Court table updated with houses.

    """
        
    # Declare variables to output results
    modified_cosine = ModifiedCosine(tolerance=params['an_mass_error'])
    selected_houses = []
    house_count = []
    house_points = []
    
    # Unresolved courts : courts to process
    unresolved_courts = list(court_table.index[court_table['Cohort_count'] > 1])

    # Resolved courts : courts with only one cohort, nothing to select
    resolved_courts = list(court_table.index[court_table['Cohort_count'] <= 1])

    # Search all unresolved courts for houses
    for court_id in unresolved_courts:
        houses = []
        
        # Get cohorts and ions in the current court
        court_cohorts = court_table.loc[court_id, "Cohort_List"].copy()
        court_ions = court_table.loc[court_id, "Ion_List"].copy()

        # Create a local cohort table for the current court
        cohorts_table = supercohorts_table.loc[court_ions, court_cohorts].copy()

        # Get contested ions and conflict table / ions with conflicting annotations
        conflict_table = list()
        # i = ion 
        for i in cohorts_table.index:
            # Conflicted = cohort numbers
            conflicted = cohorts_table.loc[i].dropna().index.tolist()
            if len(conflicted) > 1 :
                conflict_table.append((i, conflicted))
        conflict_table = pd.DataFrame(conflict_table, columns = ['conflicted_ion', 'conflicted_cohorts'])

        # Solved the conflicts
        migration = list() 
        is_mirror = list()
        
        # Process each conflicted ion
        for i in conflict_table.index:
            conflict_ion = conflict_table.loc[i, "conflicted_ion"]
            conflict_spec_id = node_table.loc[conflict_ion, 'spec_id']
            cohort_selection = list()
            
            # Process each cohort
            for j in conflict_table.loc[i, "conflicted_cohorts"]:
                cohort_ions = cohorts_table[j].dropna().index.tolist()
                cohort_ions.remove(conflict_ion)
                
                score_list = list()
                
                # Calculate cosine score with ions from each cohort, keep best product
                for k in cohort_ions:
                    k_spec_id = node_table.loc[k, 'spec_id']
                    score, n_matches = modified_cosine.pair(spectrum_list[conflict_spec_id],
                                                            spectrum_list[k_spec_id])
                    if score <= params['an_cos_threshold'] : score = 0
                    score_list.append((score, n_matches, score*n_matches))
                score_list = pd.DataFrame(score_list, columns = ['score', 'matches', 'product'])
                score_list.sort_values('product', ascending = False, inplace = True)
                cohort_selection.append((j,
                                         score_list['score'].iloc[0],
                                         score_list['matches'].iloc[0],
                                         score_list['product'].iloc[0]))
            
            # Check which cohort has the best affinity (highest product score)
            cohort_selection = pd.DataFrame(cohort_selection, columns = ['cohort', 'cos', 'matches', 'product'])
            cohort_selection.sort_values('product', ascending = False, inplace = True)
            
            # List of cohorts blocks to which ions should be migrated
            migration_blocks = list()
            for prod in cohort_selection['product'].unique():
                migration_blocks.append(list(cohort_selection['cohort'][cohort_selection['product'] == prod]))
            migration.append(migration_blocks)
            
            # Mirror : both cohort options are equivalent
            if len(migration_blocks) == 1 :
                is_mirror.append(True)
            else:
                is_mirror.append(False)
        
        # Add the data to the conflict table
        conflict_table['migration'] = migration
        conflict_table['is_mirror'] = is_mirror
        
        
        # Migrate ions to the most probable cohort until migrating = False
        migrating = True
        while migrating :
            
            # Create the migration table
            migration_table = pd.DataFrame(columns = ['ion_list', 'ion_count'])
            for i in cohorts_table.columns:
                migration_table.loc[i] = [cohorts_table[i].dropna().index.tolist(), len(cohorts_table[i].dropna())]

            # For each cohort, remove the migrated ion from the loosing cohort
            for i in conflict_table.index:
                if conflict_table.loc[i, "is_mirror"] : continue
                conflict_ion = conflict_table.loc[i, "conflicted_ion"]
                looser_cohorts = list(flatten(conflict_table.loc[i, "migration"][1:]))
                for j in looser_cohorts:
                    migration_table.loc[j, "ion_list"].remove(conflict_ion)
                    migration_table.loc[j, "ion_count"] = len(migration_table.loc[j, "ion_list"])
            
            # Remove dead cohorts (ion count <= 1)
            dead_cohorts = migration_table.index[migration_table['ion_count'] <= 1]
            cohorts_table.drop(dead_cohorts, axis = 1, inplace = True)
            
            # Continue migrating until there are no dead cohorts
            if len(dead_cohorts) == 0 : 
                migrating = False
                continue
            
            # Set changes to the conflict table
            for i in conflict_table.index:
                
                # Remove dead cohorts from the table
                conflict_table.at[i, "conflicted_cohorts"] = list(set(conflict_table.loc[i, "conflicted_cohorts"]) - set(dead_cohorts))

                # Start migration / removing dead cohorts from the table
                tmp_migration = list()
                for j in range(len(conflict_table.loc[i, "migration"])):
                    tmp_cohorts = conflict_table.loc[i, "migration"][j]
                    tmp_cohorts = list(set(tmp_cohorts) - set(dead_cohorts))
                    if len(tmp_cohorts) > 0 :
                        tmp_migration.append(tmp_cohorts)
                conflict_table.at[i, "migration"] = tmp_migration
                
                
        # Check mirrors or no-choice migrations again:
        for i in conflict_table.index:
            if len(conflict_table.loc[i, "migration"]) == 1 :
                conflict_table.loc[i, "is_mirror"] = True
        
        # Remove loosing cohorts from the conflict table, except for mirrors
        for i in conflict_table.index:
            if conflict_table.loc[i, "is_mirror"] : continue
            loosing_cohorts = list(flatten(conflict_table.loc[i, "migration"][1:]))
            conflict_ion = conflict_table.loc[i, "conflicted_ion"]
            for j in loosing_cohorts:
                cohorts_table.loc[j,conflict_ion] = None
        
        # Check for dead cohorts due to duplicated annotations:
        dead_cohorts = list()
        for i in cohorts_table.columns:
            if len(cohorts_table[i].dropna().unique()) <= 1 : dead_cohorts.append(i)
        cohorts_table.drop(dead_cohorts, axis = 1, inplace = True)
        
        # For surviving cohorts, deleted past conflicting annotations in supercohorts table
        for i in cohorts_table.columns:
            cohort_ions = cohorts_table[i].dropna().index.tolist()
            other_ions = list(set(supercohorts_table[i].dropna().index.tolist()) - set(cohort_ions))
            for j in other_ions:
                supercohorts_table.loc[j, i] = None
        
        # Inert cohorts : cohorts compatible with all others.
        inert_cohorts = list()
        for i in cohorts_table.columns:
            
            # For each cohort, get cohort ions
            cohort_ions = cohorts_table[i].dropna().index.tolist()
            fully_compatible = True
            for ion in cohort_ions:
                
                # Check if ions have multiple annotations (multiple cohorts)
                if len(cohorts_table.loc[ion].dropna().unique()) > 1 : 
                    fully_compatible = False
                    break
            if fully_compatible:
                
                # If ions are annotated by a single cohort, cohort is inert
                inert_cohorts.append(i)
                
        # remove inert cohorts from the cohorts table, won't need house selection
        cohorts_table.drop(inert_cohorts, axis = 1, inplace = True)

        # Build houses
        court_cohorts = cohorts_table.columns.tolist()
        while len(court_cohorts) > 0 :
            seed_cohort = court_cohorts[0]
            court_cohorts.remove(seed_cohort)
            # Initialise house table:
            next_list, house_list, compatibles_list, left_list = initialise_house_table(seed_cohort, court_cohorts, supercohorts_table)
            if len(house_list) == 0 : houses.append([seed_cohort])
            while len(house_list) > 0 :
                next_list, house_list, compatibles_list, left_list, houses = clean_house_table(next_list, house_list, compatibles_list, left_list, houses)
                next_list, house_list, compatibles_list, left_list = expand_house_table(next_list, house_list, compatibles_list, left_list, supercohorts_table)
        
        # Get points for the houses
        house_table = pd.DataFrame(index = range(len(houses)))
        house_table['house'] = houses
        tmp_house_points = []
        tmp_house_ions = []
        print('Assigning points to each house...')
        for i in tqdm(house_table.index):
            tmp_points = 0.0
            tmp_ion_count = []
            
            # Sum all points of each supercohort in each house
            for supercohort in house_table.loc[i, "house"]:
                tmp_ion_count += list(supercohorts_table[supercohort].dropna().index)
                cohorts = transition_table.loc[supercohort, "cohorts"]
                tmp_points_list = []
                for cohort in cohorts :
                    tmp_points_list.append(merged_table.loc[cohort, "weighted_points"])
                tmp_points += max(tmp_points_list)
            tmp_ion_count = list(set(tmp_ion_count))
            tmp_house_ions.append(len(tmp_ion_count))
            tmp_house_points.append(tmp_points)
        house_table['points'] = tmp_house_points
        house_table['ion_count'] = tmp_house_ions
        house_table.sort_values(by = ['points', 'ion_count'], ascending = False, inplace = True)
        
        # Add inert cohorts to each house, also add their points
        if len(house_table) > 0 :
            for i in inert_cohorts:
                tmp_points = []
                tmp_ion_count = len(supercohorts_table[i].dropna())
                cohort_list = transition_table.loc[i, "cohorts"]
                for cohort in cohort_list:
                    tmp_points.append(merged_table.loc[cohort, "weighted_points"])
                tmp_points = max(tmp_points)
                house_table['points'] += tmp_points
                house_table['ion_count'] += tmp_ion_count
                for j in house_table.index :
                    house_table.loc[j, 'house'].append(i)
            for i in house_table.index:
                house_table.at[i, "house"] = [house_table.loc[i, "house"]]
        else :
            tmp_points = []
            tmp_ion_count = 0
            for i in inert_cohorts:
                cohort_list = transition_table.loc[i, "cohorts"]
                tmp_ion_count += len(cohort_list)
                for cohort in cohort_list :
                    tmp_points.append(merged_table.loc[cohort, "weighted_points"])
            tmp_points = max(tmp_points)
            house_table = pd.DataFrame(columns = ['house', 'points', 'ion_count'])
            house_table.loc[0] = [[inert_cohorts], tmp_points, tmp_ion_count]
            
        # Select the best house
        selected_houses.append(house_table['house'].iloc[0])
        house_count.append(len(house_table['house'].iloc[0]))
        house_points.append(house_table['points'].max())
    
    # Once unresolved courts are dealt with, add resovled courts
    for court_id in resolved_courts:
        selected_houses.append([court_table.loc[court_id, "Cohort_List"]])
        house_count.append(1)
        supercohort = court_table.loc[court_id, "Cohort_List"][0]
        cohorts = transition_table.loc[supercohort, "cohorts"]
        tmp_points = []
        for cohort in cohorts:
            tmp_points.append(merged_table.loc[cohort, "weighted_points"])
        house_points.append(max(tmp_points))
        
    # Add results to the court table and export
    court_table["selected_houses"] = selected_houses
    court_table["house_count"] = house_count
    court_table["house_points"] = house_points
    return court_table

def cross_sample_report(court_table, cross_annotations, cross_points, cross_courts,
                        cross_houses, cross_rules, cross_neutrals,
                        sample_base_name : str, supercohorts_table,
                        adduct_table_primary, adduct_table_merged, node_table,
                        ion_mode, duplicate_df, spectrum_list, params):
    """Reports the data for the sample being processed into the cross tables 
    before exporting them as csv files.
    """
    
    # Add sample to tables
    cross_annotations[sample_base_name] = [None]*len(cross_annotations)
    cross_points[sample_base_name] = [None]*len(cross_points)
    cross_courts[sample_base_name] = [None]*len(cross_courts)
    cross_houses[sample_base_name] = [None]*len(cross_courts)
    cross_rules[sample_base_name] = [None]*len(cross_courts)
    cross_neutrals[sample_base_name] = [None]*len(cross_courts)
 
    print("Reporting results...")
    for i in tqdm(court_table.index):
        
        # For each court (i), process each house it contains
        for house in range(len(court_table.loc[i, "selected_houses"])):
            
            # Get each cohort in the house
            cohorts = court_table.loc[i, "selected_houses"][house]
            for cohort in cohorts:
                
                # Get each ion annotation from the cohort
                annotations = supercohorts_table[cohort].dropna()
                for ion in annotations.index:
                    
                    # For each ion, get spec_id
                    spec_id = node_table.loc[ion, 'spec_id']
                    
                    # Get its neutral mass
                    neutral_mass = adduct_table_merged["Adduct_code"] == annotations[ion]
                    neutral_mass = neutral_mass.index[neutral_mass][0]
                    neutral_mass = neutral_mass_calculator(node_table.loc[ion, params['mz_field']],
                                                           adduct_table_merged.loc[neutral_mass, "Adduct_mass"],
                                                           adduct_table_merged.loc[neutral_mass, "Mol_multiplier"],
                                                           adduct_table_merged.loc[neutral_mass, "Charge"])
                    Species_rule = Validator_choice(annotations[ion], ion_mode)
                    rule_validation = Species_rule(params['an_prec_mass_error'], spec_id, spectrum_list)
                    #ion_id = int(node_table.loc[ion, params['index_col']])
                    if cross_annotations.loc[ion, sample_base_name] == None :
                        cross_annotations.loc[ion, sample_base_name] = annotations[ion]
                        cross_houses.loc[ion, sample_base_name] = str(house)
                        cross_rules.loc[ion, sample_base_name] = str(rule_validation)
                        cross_neutrals.loc[ion, sample_base_name] = str(neutral_mass)
                    else:
                        cross_annotations.loc[ion, sample_base_name] += "&" + annotations[ion]
                        cross_houses.loc[ion, sample_base_name] += "&" + str(house)
                        cross_rules.loc[ion, sample_base_name] += "&" + str(rule_validation)
                        cross_neutrals.loc[ion, sample_base_name] += "&" + str(neutral_mass)
                    cross_points.loc[ion, sample_base_name] = court_table.loc[i, "house_points"]
                    cross_courts.loc[ion, sample_base_name] = i

    # Add duplicate annotations :
    for i in duplicate_df.index:
        
        # Get ion 2 (duplicate) ID
        dup_id = duplicate_df.loc[i, "ion_2_idx"]

        # If ion 2 is already annotated, skip
        if type(cross_annotations.loc[dup_id, sample_base_name]) == str : continue

        # Get ion 2 adduct as text
        ion_2_adduct = duplicate_df.loc[i, "ion_2_adduct"]
        ion_2_adduct = adduct_table_primary.loc[ion_2_adduct, "Adduct_code"]

        # Get the selected ion ID
        selected_id = duplicate_df.loc[i, "selected_ion"]

        # Get annotation for the selected ion
        selected_adduct = cross_annotations.loc[selected_id, sample_base_name]
        
        # If the annotations are the same, report the duplicate annotation in the tables
        if ion_2_adduct == selected_adduct:
            cross_annotations.loc[dup_id, sample_base_name] = ion_2_adduct
            cross_houses.loc[dup_id, sample_base_name] = cross_houses.loc[selected_id, sample_base_name]
            cross_points.loc[dup_id, sample_base_name] = cross_points.loc[selected_id, sample_base_name]
            cross_courts.loc[dup_id, sample_base_name] = cross_courts.loc[selected_id, sample_base_name]
            cross_neutrals.loc[dup_id, sample_base_name] = cross_neutrals.loc[selected_id, sample_base_name]
            cross_rules.loc[dup_id, sample_base_name] = cross_rules.loc[selected_id, sample_base_name]
    return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals

def get_secondary_adducts(cross_annotations, cross_points, cross_courts,
                          cross_houses, cross_rules, cross_neutrals, sample_base_name,
                          node_table, spectrum_list, adduct_table_primary, 
                          adduct_table_secondary, ion_mode, params):
    """After processing a sample and annotating ions with the adducts in the 
    primary adducts table, the remaining ions are annotated with adducts from
    the secondary adducts table. This is done on the basis of the already
    computed molecules and consumes less resources. 
    """
    
    rt_error = params['an_rt_error']
    if params['rt_unit'] == "m":
        rt_error = rt_error/60
    
    
    modified_cosine = ModifiedCosine(tolerance=params['an_mass_error'])
    
    # Get annotated ions
    annotated = list(cross_neutrals[sample_base_name].dropna().index)

    # Get unnannotated ions
    unnannotated = list(set(node_table.index.astype(int)) - set(annotated))
    unnannotated_table = node_table.loc[unnannotated]
    
    for ion_id in annotated:
        
        # Get ion attributes
        spec_id = node_table.loc[ion_id, "spec_id"]
        neutrals = cross_neutrals.loc[ion_id, sample_base_name].split('&')
        annotations = cross_annotations.loc[ion_id, sample_base_name].split('&')
        neutral_rt = node_table.loc[ion_id, params['rt_field']]
        coelution_table = unnannotated_table[unnannotated_table[params['rt_field']].between(neutral_rt - rt_error,  neutral_rt + rt_error, inclusive = "both")]
        houses = cross_houses.loc[ion_id, sample_base_name].split('&')

        # Go through each house in case multiple annotations for a single ion
        for k in range(len(houses)): 
            house = houses[k]
            neutral = float(neutrals[int(k)])
            main_adduct = annotations[int(k)]
            main_group = adduct_table_primary["Group"][adduct_table_primary["Adduct_code"] == main_adduct].iloc[0]
            
            # Search for secondary adducts
            for i in adduct_table_secondary.index:
                tmp_adduct = adduct_table_secondary.loc[i, "Adduct_code"]
                tmp_mz= ion_mass_calculator(neutral, adduct_table_secondary.loc[i, "Adduct_mass"],
                                          adduct_table_secondary.loc[i, "Mol_multiplier"],
                                          adduct_table_secondary.loc[i, "Charge"])
                tmp_group = adduct_table_secondary.loc[i, "Group"]
                tmp_table = coelution_table[params['mz_field']].between(tmp_mz - params['an_mass_error'], tmp_mz + params['an_mass_error'], inclusive = "both")

                if sum(tmp_table) != 0:
                    tmp_table = coelution_table[tmp_table]
                    for j in tmp_table.index:
                        valid = False
                        Species_rule = Validator_choice(tmp_adduct, ion_mode)
                        j_spec_id = node_table.loc[j, "spec_id"]
                        if main_group == tmp_group:
                            tmp_cos, n_matches = modified_cosine.pair(spectrum_list[spec_id],
                                                                    spectrum_list[j_spec_id])
                            if tmp_cos >= params['an_cos_threshold'] : valid = True
                        elif Species_rule(params['an_prec_mass_error'], j_spec_id, spectrum_list) != 0:
                            valid = True
                        if valid :
                            cross_courts.loc[j, sample_base_name] = cross_courts.loc[ion_id, sample_base_name]
                            cross_points.loc[j, sample_base_name] = cross_points.loc[ion_id, sample_base_name]
                            
                            if cross_annotations.loc[j, sample_base_name] == None:
                                cross_annotations.loc[j, sample_base_name] = tmp_adduct
                                cross_houses.loc[j, sample_base_name] = house
                                cross_neutrals.loc[j, sample_base_name] = str(neutral)
                                cross_rules.loc[j, sample_base_name] = str(Species_rule(params['an_prec_mass_error'], j_spec_id, spectrum_list))
                            elif tmp_adduct not in cross_annotations.loc[j, sample_base_name].split('&'):
                                cross_annotations.loc[j, sample_base_name] += "&" + tmp_adduct
                                cross_houses.loc[j, sample_base_name] += "&" + house
                                cross_neutrals.loc[j, sample_base_name] += "&" + str(neutral)
                                cross_rules.loc[j, sample_base_name] += "&" + str(Species_rule(params['an_prec_mass_error'], j_spec_id, spectrum_list))
    return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals


######################################################## Cross sample functions

def refresh_status(node_table, edge_table):
    """
    Refresh node status in node table.

    Parameters
    ----------
    node_table : pandas.DataFrame
        Table containing node / ion data.
    edge_table : pandas.DataFrame
        Table containing edge / ion pairs data.

    Returns
    -------
    node_table : pandas.DataFrame
        Node table with updated status.
    edge_table : pandas.DataFrame
        Edge table with updated status.

    """
    node_table['status'] = ['singleton']*len(node_table)
    precursors = list(set(edge_table['node_1']) - set(edge_table['node_2']))
    fragments = list(set(edge_table['node_2']) - set(precursors))
    print('Updating precursor node status:')
    for i in tqdm(precursors):
        node_table.loc[i, "status"] = "precursor"
    print('Updating fragment node status:')
    for i in tqdm(fragments):
        node_table.loc[i, "status"] = "fragment"
    edge_table = edge_table.replace({np.nan: None})
    return node_table, edge_table

def get_merged_tables(input_files, ion_mode : str, params : dict):
    """
    Runs through all samples in input_files (all user-submitted samples) and
    merges the fragnotator outputted tables into a single node and edge table,
    with the "status" updated for each ion in the node table.

    Parameters
    ----------
    input_files : pandas.DataFrame
        Dataframe containing the file names.

    ion_mode : str
        Ion mode of the data.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    merged_node_table : pandas.DataFrame
        Node table produced by merging data from all samples.
    merged_edge_table : pandas.DataFrame
        edge table produced by merging data from all samples.

    """

    # Load parameters
    if ion_mode == "POS":
        in_path_csv= params['pos_out_2']
    else:
        in_path_csv= params['neg_out_2']
    idx_column = params['index_col']
    
    # Initialise merged tables
    merged_node_table = pd.DataFrame()
    merged_edge_table = pd.DataFrame(columns = ['node_1', 'node_2', 'matched_peaks',
                                                'total_peaks', 'matching_score',
                                                'rt_gap', 'mz_gap', 'status',
                                                'Fragnotation', 'All_annotations'])
    
    # Add data from all samples to the merged tables
    print('Merging node and edge tables from all samples:')
    for x in tqdm(input_files.index) :
        
        # Load tables
        file_basename = input_files.at[x, 'base_name']    
        edge_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_edges.csv", 
                                 index_col = 'Index')
        node_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_nodes.csv", 
                                index_col = idx_column)


        # Drop the status col, which will be updated later
        node_table.drop('status', axis = 1, inplace = True)

        # Add only new ions to the merged node table
        new_ions = list(set(node_table.index) - set(merged_node_table.index))
        merged_node_table = merged_node_table.append(node_table.loc[new_ions], ignore_index = False)

        # Add only non-singleton edges
        edge_table = edge_table[edge_table['status'] != "singleton"]
        
        # Create an ID combining the node 1 and node 2 pair
        edge_ids = edge_table['node_1'].astype(str) + "|" + edge_table['node_2'].astype(str)
        edge_table.set_index(edge_ids, inplace = True)
        
        # Only add these unique new IDs to the edge table
        new_edges = list(set(edge_table.index) - set(merged_edge_table.index))
        merged_edge_table = merged_edge_table.append(edge_table.loc[new_edges], ignore_index = False)
    
    # Sort tables 
    merged_node_table.sort_index(inplace = True)
    merged_edge_table.sort_values("node_1", inplace = True)
    merged_edge_table.reset_index(drop = True, inplace = True)
    
    # Refresh edge status
    merged_node_table, merged_edge_table = refresh_status(merged_node_table, merged_edge_table)

    return merged_node_table, merged_edge_table


def get_sample_houses(cross_houses, sample : str):
    """
    Get hourse IDs from a sample using a cross_houses table.

    Parameters
    ----------
    cross_houses : pandas.DataFrame
        Dataframe containing the house IDs in each sample and each ion.
    sample : str
        Sample name (column in cross_houses.

    Returns
    -------
    sample_houses : list
        List of houses for the sample.

    """
    sample_houses = cross_houses[sample].dropna().str.split('&')
    for i in sample_houses.index:
        sample_houses[i] = list(map(int, sample_houses[i]))
    return sample_houses

def resolver_level_2(ion_id : int, sample : str, cross_annotations, cross_points,
                     cross_courts, cross_houses, cross_rules, cross_neutrals):
    """
    Resolve annotations if level 1 resolver was not enough.

    Parameters
    ----------
    ion_id : int
        Ion ID.
    sample : str
        Sample name.
    cross_annotations : pandas.DataFrame
        Table containing annotations for each ion in each sample.
    cross_points : pandas.DataFrame
        Table containing points for each ion annotation in each sample.
    cross_courts : pandas.DataFrame
        Table containing court IDs for each ion in each sample.
    cross_houses : pandas.DataFrame
        Table containing house IDs for each ion in each sample.
    cross_rules : pandas.DataFrame
        Table containing rule points for each ion in each sample.
    cross_neutrals : pandas.DataFrame
        Table containing neutral masses for each ion in each sample.

    Returns
    -------
    Updated cross samples dataframes.

    """
    
    # Get the sample houses and courts
    sample_houses = get_sample_houses(cross_houses, sample)
    sample_courts = cross_courts[sample].dropna()

    # Create a conflict table to store
    conflict_table = pd.DataFrame()
    conflict_table['annotation'] = cross_annotations.loc[ion_id, sample].split('&')
    conflict_table['house_IDs'] = sample_houses[ion_id]
    conflict_table['neutrals'] = cross_neutrals.loc[ion_id, sample].split('&')
    conflict_table['rule_points'] = cross_rules.loc[ion_id, sample].split('&')
    conflict_table['conflicting_ions'] = [[]]*len(conflict_table)
    conflict_table['points'] = [0.0]*len(conflict_table)
    court = sample_courts[ion_id]
    court_ions =set(sample_courts.index[sample_courts == court])
    
    # Get conflicting ions
    conflicting_ions = []
    for i in conflict_table.index:
        house = conflict_table.loc[i, "house_IDs"]
        house_ions = set([k for k in sample_houses.index if house in sample_houses[k]])
        conflicting_ions.append(list(court_ions.intersection(house_ions)))
        
    # ???
    for i in conflict_table.index:
        house_ions = set(conflicting_ions[i])
        other_ions = conflicting_ions.copy()
        other_ions.remove(other_ions[i])
        other_ions = set(flatten(other_ions))
        conflict_table.at[i, "conflicting_ions"] = list(house_ions - other_ions)
    
    # Award points
    for i in conflict_table.index:
        for j in conflict_table.loc[i, "conflicting_ions"]:
            j_annotation = cross_annotations.loc[j, sample]
            if "&" in j_annotation : 
                sys.exit(f"ERROR: conflict inside conflict for ion pair {ion_id}-{j}, sample {sample}")
            conflict_samples = cross_annotations.loc[j].dropna() == j_annotation
            conflict_samples = conflict_samples.index[conflict_samples]
            conflict_table.loc[i, "points"] += cross_points.loc[j, conflict_samples].max()
    
    # Select by points
    conflict_table.sort_values("points", ascending = False, inplace = True)
    corrected_samples = cross_annotations.loc[ion_id] == cross_annotations.loc[ion_id, sample]
    corrected_samples = corrected_samples.index[corrected_samples]
    
    # Apply selection to cross tables
    for i in corrected_samples:
        cross_annotations.loc[ion_id, i] = conflict_table["annotation"].iloc[0]
        cross_houses.loc[ion_id, i] = str(conflict_table["house_IDs"].iloc[0])
        cross_neutrals.loc[ion_id, i] = conflict_table["neutrals"].iloc[0]
        cross_rules.loc[ion_id, i] = conflict_table["rule_points"].iloc[0]
    conflict_table.drop(conflict_table.index[0], inplace = True)
    
    # For dead annotations resulting of the choice, cleanse
    for i in conflict_table.index:
        for j in conflict_table.loc[i, "conflicting_ions"]:
            j_annotation = cross_annotations.loc[j, sample]
            corrected_samples = cross_annotations.loc[j] == j_annotation
            corrected_samples = corrected_samples.index[corrected_samples]
            for k in corrected_samples:
                cross_annotations.loc[j, k] = None
                cross_houses.loc[j, k] = None
                cross_neutrals.loc[j, k] = None
                cross_rules.loc[j, k] = None
                cross_courts.loc[j, k] = None
                cross_points.loc[j, k] = None
    return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals

def resolver_level_1(ion_id : int, sample : str, cross_annotations, cross_points, cross_courts,
                     cross_houses, cross_rules, cross_neutrals):
    """
    

    Parameters
    ----------
    ion_id : int
        Ion ID.
    sample : str
        Sample name.
    cross_annotations : pandas.DataFrame
        Table containing annotations for each ion in each sample.
    cross_points : pandas.DataFrame
        Table containing points for each ion annotation in each sample.
    cross_courts : pandas.DataFrame
        Table containing court IDs for each ion in each sample.
    cross_houses : pandas.DataFrame
        Table containing house IDs for each ion in each sample.
    cross_rules : pandas.DataFrame
        Table containing rule points for each ion in each sample.
    cross_neutrals : pandas.DataFrame
        Table containing neutral masses for each ion in each sample.

    Returns
    -------
    Updated cross samples dataframes.

    """
    
    # Create a conflict table
    conflict_table = pd.DataFrame()
    conflict_table['annotation'] = cross_annotations.loc[ion_id, sample].split('&')
    conflict_table['house'] = cross_houses.loc[ion_id, sample].split('&')
    conflict_table['neutral'] = cross_neutrals.loc[ion_id, sample].split('&')
    conflict_table['rule_points'] = cross_rules.loc[ion_id, sample].split('&')
    conflict_table['points'] = [0.0]*len(conflict_table)
    
    # Get points from all samples to reward the most likely annotation
    for i in conflict_table.index:        
        tmp_samples = cross_annotations.loc[ion_id] == conflict_table.loc[i, "annotation"] 
        tmp_samples = tmp_samples.index[tmp_samples]
        conflict_table.loc[i, "points"] = cross_points.loc[ion_id, tmp_samples].max()
    selection_table = conflict_table[conflict_table['points'] == conflict_table['points'].max()]

    # If a choice yet could not be base, pass onto level 2 resolver
    if len(selection_table) != 1 :
        resolver_level_2(ion_id, sample)
        return
    else:
        # If a choice was made, report it to the tables
        conflict_table.drop(selection_table.index, inplace = True)
        correct_samples = cross_annotations.loc[ion_id] == cross_annotations.loc[ion_id, sample]
        correct_samples = correct_samples.index[correct_samples]
        for tmp_sample in correct_samples:
            cross_annotations.loc[ion_id, tmp_sample] = selection_table['annotation'].iloc[0]
            cross_houses.loc[ion_id, tmp_sample] = selection_table['house'].iloc[0]
            cross_neutrals.loc[ion_id, tmp_sample] = selection_table['neutral'].iloc[0]
            cross_rules.loc[ion_id, tmp_sample] = selection_table['rule_points'].iloc[0]
        
        # For dead annotations resulting of the choice, cleanse
        for i in conflict_table.index:
            correct_samples = cross_annotations.loc[ion_id] == conflict_table.loc[i, "annotation"]
            correct_samples = correct_samples.index[correct_samples]
            for tmp_sample in correct_samples :
                cross_annotations.loc[ion_id, tmp_sample] = None
                cross_houses.loc[ion_id, tmp_sample] = None
                cross_neutrals.loc[ion_id, tmp_sample] = None
                cross_rules.loc[ion_id, tmp_sample] = None
                cross_courts.loc[ion_id, tmp_sample] = None
                cross_points.loc[ion_id, tmp_sample] = None
        return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals

def annotation_resolver(cross_annotations, cross_points, cross_courts, cross_houses,
                        cross_rules, cross_neutrals):
    """
    Resolves the yet sample-specific unresolved ions by cross examining 
    between all samples

    Parameters
    ----------
    cross_annotations : pandas.DataFrame
        Table containing annotations for each ion in each sample.
    cross_points : pandas.DataFrame
        Table containing points for each ion annotation in each sample.
    cross_courts : pandas.DataFrame
        Table containing court IDs for each ion in each sample.
    cross_houses : pandas.DataFrame
        Table containing house IDs for each ion in each sample.
    cross_rules : pandas.DataFrame
        Table containing rule points for each ion in each sample.
    cross_neutrals : pandas.DataFrame
        Table containing neutral masses for each ion in each sample.

    Returns
    -------
    Updated cross samples dataframes.

    """

    print("Cross-examining unresolved annotations...")
    for sample in tqdm(cross_annotations.columns):
        
        # In each sample, check if there are multiple annotations (contain "&")
        tmp_annotations = cross_annotations[sample].dropna()
        tmp_annotations = tmp_annotations[tmp_annotations.str.contains('&')]
        if len(tmp_annotations) == 0 : continue

        # Attempt to chose a single annotation using a resolver function
        for ion_id in tmp_annotations.index:
            resolver_level_1(ion_id, sample)
    
    # Flush out annotations supported only by a single ion (after the above filtering)
    print("Flushing dead annotations...")
    for sample in tqdm(cross_annotations.columns):
        courts = cross_courts[sample].dropna()
        for court in courts.unique():
            court_samples = courts.index[courts==court]
            houses = cross_houses.loc[court_samples, sample]
            for house in houses.unique():
                if (houses == house).sum() <= 1:
                    purge_index = houses[houses == house].index
                    for i in purge_index:
                        cross_annotations.loc[i, sample] = None
                        cross_courts.loc[i, sample] = None
                        cross_houses.loc[i, sample] = None
                        cross_neutrals.loc[i, sample] = None
                        cross_points.loc[i, sample] = None
                        cross_rules.loc[i, sample] = None
    return cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals

def get_cross_courts(cross_annotations, cross_courts, cross_houses):
    """
    Produce a cross-sample court table to select the most likely houses
    valid across samples.

    Parameters
    ----------
    cross_annotations : pandas.DataFrame
        Table containing annotations for each ion in each sample.
    cross_courts : pandas.DataFrame
        Table containing court IDs for each ion in each sample.
    cross_houses : pandas.DataFrame
        Table containing house IDs for each ion in each sample.

    Returns
    -------
    court_table : pandas.DataFrame
        Court table containing courts deduced across samples.
    singletons : list
        List of unannotated ions.

    """

    # Isolate singleton ions / ions without any annotation
    print('Court processing - isolating singleton ions :')
    ion_pool = list(cross_annotations.index)
    court_table = pd.DataFrame(columns = ['ion_list', 'court_size'])
    singletons = list()
    for ion in tqdm(ion_pool):
        if len(cross_annotations.loc[ion].dropna()) == 0:
            singletons.append(ion)

    # Remove singleton ions from the ion pool
    ion_pool = list(set(ion_pool) - set(singletons))
    total_ions = len(ion_pool)

    # Start processing cross courts
    while len(ion_pool) > 0 :
        perc = round((1-(len(ion_pool)/total_ions))*100,1)
        sys.stdout.write("\rCourt processing - creating courts : {0}%".format(perc))
        sys.stdout.flush()
        
        # Initiate court
        court = [ion_pool[0]]
        court_size = 0
        while court_size != len(court):
            
            # Get samples involved with court
            court_size = len(court)
            samples = []
            for ion in court:
                samples += list(cross_annotations.loc[ion].dropna().index)
            if len(samples) == 0: sys.exit("Error : ion associated to 0 samples")
            samples = list(set(samples))
            
            # Add new ions associated to court
            new_ions= court.copy()
            for sample in samples:
                for ion in court:
                    sample_courts = cross_courts[sample].dropna()
                    sample_houses = get_sample_houses(cross_houses, sample)
                    tmp_court = cross_courts.loc[ion, sample]
                    if tmp_court == None : continue
                    house = int(cross_houses.loc[ion, sample])
                    court_ions = set(sample_courts.index[sample_courts == tmp_court])
                    house_ions = set([k for k in sample_houses.index if house in sample_houses[k]])
                    new_ions += list(court_ions.intersection(house_ions))
            court = list(set(new_ions))
        court.sort()
        tmp_idx = len(court_table)
        court_table.loc[tmp_idx, 'ion_list'] = court
        court_table.loc[tmp_idx, 'court_size'] = len(court)
        ion_pool = list(set(ion_pool) - set(court))
    court_table.sort_values('court_size', ascending = False, inplace = True)
    court_table.reset_index(drop = True, inplace = True)
    return court_table, singletons



######################################################### End functions / Hell

def neutral_compatibility(merged_neutral_table, adduct_table_merged):
    """For each neutral in the neutral table, adds a list of compatible neutrals
    found in the same table.
    """
    merged_neutral_table['compatibles'] = [None]*len(merged_neutral_table)
    for i in merged_neutral_table.index:
        merged_neutral_table.at[i, "compatibles"] = list()
    for neutral_1 in merged_neutral_table.index:
        annotations_1 = merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()
        annotations_1 = set(flatten(annotations_1))      
        for neutral_2 in merged_neutral_table.index:
            if neutral_1 == neutral_2 : continue
            annotations_2 = merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()
            annotations_2 = set(flatten(annotations_2))
            if len(annotations_1.intersection(annotations_2)) == 0 :
                merged_neutral_table.loc[neutral_2, 'compatibles'].append(neutral_1)
    for i in merged_neutral_table.index:
        merged_neutral_table.loc[i, "compatibles"] = list(set(merged_neutral_table.loc[i, "compatibles"]))
        merged_neutral_table.loc[i, "compatibles"].sort()        
    return merged_neutral_table

def get_crossneutral_table(samples_table, merged_node_table, cross_annotations,
                           cross_neutrals, params: dict):
    """
    Produces a neutral_table, containing among other things, the computed
    neutral mass for each ion in the samples table.

    Parameters
    ----------
    samples_table : pandas.DataFrame
        Sample table (ion x samples) with annotations
    merged_node_table : pandas.DataFrame
        Node table produced by merging data from all samples.
    cross_annotations : pandas.DataFrame
        Table containing annotations for each ion in each sample.
    cross_neutrals : pandas.DataFrame
        Table containing neutral masses for each ion in each sample.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    neutral_table : pandas.DataFrame
        Table containing data for all neutrals produced.

    """

    # Initiate neutral table
    neutral_table = pd.DataFrame(columns = ['ion_ID', 'mz', 'rt', 'annotation', 'neutral'])

    # Get annotation and neutral information for each ion
    for ion in samples_table.index:
        annotations = samples_table.loc[ion].dropna().unique()
        mz = merged_node_table.loc[ion, params['mz_field']]
        rt = merged_node_table.loc[ion, params['rt_field']]
        
        # Report data to the neutral table
        for annotation in annotations :
            tmp_sample = cross_annotations.columns[cross_annotations.loc[ion] == annotation][0]
            neutral = float(cross_neutrals.loc[ion, tmp_sample])
            neutral_table.loc[len(neutral_table)] = [ion, mz, rt, annotation, neutral]
    return neutral_table

def neutral_merger(neutral_table, adduct_table_merged, params : dict):
    """
    Merges neutrals from the neutral table based on their retention time
    and molecular masses. Returns a merged neutral table, with the neutral's
    retention time, points accumulated, the ions involved and their annotations.

    Parameters
    ----------
    neutral_table : pandas.DataFrame
        Table containing data for all neutrals from a set of ions
    adduct_table_merged : pandas.DataFrame
        Merged primary and secondary adduct tables loaded from the parameters folder.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    neutral_table_2 : pandas.DataFrame
        Merged neutral table

    """

    
    # Load parameters
    rt_error = params['an_rt_error']
    mass_error = params['an_mass_error']
    
    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        rt_error = rt_error/60
    
    # Create the ouput neutral table
    neutral_table_2 = pd.DataFrame(columns = ["mass", "rt", "points"] + list(adduct_table_merged['Adduct_code']))

    # 
    neutral_idx = neutral_table.index
    while len(neutral_idx) > 1:
        
        # Select the first neutral to start merging
        neutral_pool = [neutral_idx[0]]
        pool_size = 0
        pool_table = neutral_table.loc[neutral_pool]
        
        # Start merging starting with the first neutral
        while len(neutral_pool) != pool_size:
            pool_size = len(neutral_pool)
            tmp_rt = (pool_table['rt'].min() + pool_table['rt'].max())/2
            tmp_neutral = (pool_table['neutral'].min() + pool_table['neutral'].max())/2
            pool_table = neutral_table[neutral_table['rt'].between(tmp_rt - rt_error, tmp_rt + rt_error)]
            pool_table = pool_table[pool_table['neutral'].between(tmp_neutral - mass_error, tmp_neutral + mass_error)]
            neutral_pool = list(pool_table.index)
        
        # Add the merged data to the output neutral table
        tmp_idx = len(neutral_table_2)
        neutral_table_2.loc[tmp_idx] = [None]*len(neutral_table_2.columns)
        neutral_table_2.loc[tmp_idx, "mass"] = round(pool_table['neutral'].mean(), 4)
        neutral_table_2.loc[tmp_idx, "rt"] = round(pool_table['rt'].mean(), 2)
        
        # Add annotations to the output table
        for annotation in pool_table['annotation'].unique():
            neutral_table_2.loc[tmp_idx, annotation] = list(pool_table['ion_ID'][pool_table["annotation"] == annotation])

        # Calculate points for the neutral based on adduct complexity
        adducts = list(neutral_table_2.loc[tmp_idx, adduct_table_merged['Adduct_code']].dropna().index)
        tmp_points = 0
        for adduct in adducts:
            tmp_points += 1/(adduct_table_merged["Complexity"][adduct_table_merged["Adduct_code"] == adduct].iloc[0])
        neutral_table_2.loc[tmp_idx, "points"] = tmp_points
        
        # Remove merged neutral from the neutral_idx (left to analyse)
        neutral_idx = list(set(neutral_idx) - set(neutral_pool))
    return neutral_table_2   

def initialise_house_table_neutrals(neutral_seed, merged_neutral_table):
    """Starts house building by combining compatible cohorts, using an initial
    cohort as a neutral_seed.
    """
    next_list = merged_neutral_table.loc[neutral_seed, "compatibles"].copy()
    house_list = [[neutral_seed]]*len(next_list)
    compatibles_list = [None]*len(next_list)
    left_list = [0]*len(next_list)
    added = []
    for i in range(len(next_list)):
        added.append(next_list[i])
        tmp_compatibles = list(set(next_list).intersection(merged_neutral_table.loc[next_list[i], "compatibles"]))
        tmp_compatibles = list(set(tmp_compatibles) - set(added))
        compatibles_list[i] = tmp_compatibles
        left_list[i] = len(compatibles_list[i])
    return next_list, house_list, compatibles_list, left_list

def clean_house_table_neutrals(next_list, house_list, compatibles_list, left_list, houses):
    """Removes houses that finished expanding (no compatible cohorts left) from
    the house table, and exports them to the houses list.
    """
    export_idx = [i for i in range(len(left_list)) if left_list[i] == 0]
    export_idx.sort(reverse = True)
    for i in export_idx:
        new_house = house_list[i] + [next_list[i]]
        new_house = list(set(new_house))
        new_house.sort()
        houses.append(new_house) 
        next_list.pop(i)
        house_list.pop(i)
        compatibles_list.pop(i)
        left_list.pop(i)
    return next_list, house_list, compatibles_list, left_list, houses

def expand_house_table_neutrals(next_list, house_list, compatibles_list, left_list, merged_neutral_table):
    """Expands the houses in the house table by appending to them one of each
    possible compatible cohorts.
    """
    drop_list = []
    for i in range(len(next_list)):
        drop_list.append(i)
        new_house = house_list[i] + [next_list[i]]
        added = []
        compatibles = compatibles_list[i]
        for j in compatibles:
            #j = compatibles[0]
            added.append(j)
            next_list.append(j)
            house_list.append(new_house)
            tmp_compatibles = list(set(compatibles).intersection(merged_neutral_table.loc[j, "compatibles"]))
            tmp_compatibles = list(set(tmp_compatibles) - set(added))
            compatibles_list.append(tmp_compatibles)
            left_list.append(len(set(compatibles).intersection(merged_neutral_table.loc[j, "compatibles"])))
    drop_list.sort(reverse = True)
    for i in drop_list:
        next_list.pop(i)
        house_list.pop(i)
        compatibles_list.pop(i)
        left_list.pop(i)
    return next_list, house_list, compatibles_list, left_list  

def cross_neutral_selection(spectrum_list, cross_court_table, cross_annotations,
                            cross_neutrals, merged_node_table, merged_edge_table,
                            adduct_table_merged, ion_mode, params):
    """Select the most likeley neutrals from the cross_court_table and report
    the data on the merged node and edge tables.
    """
    
    # Load parameters
    mz_field = params['mz_field']
    rt_field = params['rt_field']
    charge_field = params['charge_field']
    prec_mass_error = params['an_prec_mass_error']
    mass_error = params['an_mass_error']
    run_bnr = params['an_run_bnr']
    cosine_threshold = params['an_cos_threshold']
    cosine_hardthreshold = params['an_hardcos_threshold']
    if ion_mode == "POS":
        bnr_list = set(params['an_bnr_pos'])
    else:
        bnr_list = set(params['an_bnr_neg'])
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    
    # For each court, find houses:
    print('')
    print("Court processing - selecting houses and reporting results")
    for court in tqdm(cross_court_table.index):
        
        # Get ions from court
        ion_list = cross_court_table.loc[court, "ion_list"]

        # Get the involved samples
        samples = []
        for ion in ion_list:
            samples += list(cross_annotations.loc[ion].dropna().index)
        samples = list(set(samples))
        samples.sort()
        
        # Create a sample table (ion x samples) and add annotations
        samples_table = pd.DataFrame(index = ion_list, columns = samples)
        for ion in ion_list:
            samples = list(cross_annotations.loc[ion].dropna().index)
            for sample in samples:
                samples_table.loc[ion, sample] = cross_annotations.loc[ion, sample]
    
        # Produce neutral table
        neutral_table = get_crossneutral_table(samples_table, merged_node_table, cross_annotations, cross_neutrals, params)
        
        # Do neutral merging:
        merged_neutral_table = neutral_merger(neutral_table, adduct_table_merged, params)

        # Eliminate stillborn neutrals (1 or less annotations):
        drop_list = list()
        for i in merged_neutral_table.index:
            if len(merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna()) <= 1 :
                drop_list.append(i)
        merged_neutral_table.drop(drop_list, inplace = True)
        
        ##################################################################
        ################################################## BNR BEGINS HERE
        if run_bnr:
            drop_list = list()
            for i in merged_neutral_table.index:
                tmp_cos_list = list()
                adducts = merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna().index.tolist()
                if len(bnr_list.intersection(adducts)) == 0:
                    for adduct_1 in adducts:
                        adduct_ions_1 = merged_neutral_table.loc[i, adduct_1]
                        adducts.remove(adduct_1)
                        for ion_1 in adduct_ions_1:
                            ion_1_spec_id = merged_node_table.loc[ion_1, "spec_id"]
                            for adduct_2 in adducts:
                                adduct_ions_2 = merged_neutral_table.loc[i, adduct_2]
                                for ion_2 in adduct_ions_2:
                                    ion_2_spec_id = merged_node_table.loc[ion_2, "spec_id"]
                                    score, n_matches = modified_cosine.pair(spectrum_list[ion_1_spec_id],
                                                                            spectrum_list[ion_2_spec_id])
                                    if n_matches <= 2 : score = 0.0
                                    tmp_cos_list.append(score)
                    if max(tmp_cos_list) < cosine_hardthreshold:
                        drop_list.append(i)
            merged_neutral_table.drop(drop_list, inplace = True)
                                
        ################################################## BNR ENDS HERE
        ##################################################################
        
        # If onle one neutral or less remain, report
        if len(merged_neutral_table) <= 1 :
            selected_house = list(merged_neutral_table.index)
            merged_neutral_table_selected = merged_neutral_table
        
        # Else : continue selection
        else : 

            # Create the selected neutrals table
            merged_neutral_table_selected = pd.DataFrame()
            
            end_selection = False
            while end_selection == False :
                # Do ion affinity here:
###############################################################################
###################################################### ION AFFINITY BEGINS HERE
                # Ion affinity on neutral remains:
                # Find conflicted ions
                contested_ions = pd.DataFrame(columns = ["conflicted_neutrals"])
                neutral_pool = list(merged_neutral_table.index)
                while len(neutral_pool) > 1 :
                    neutral_1 = neutral_pool[0]
                    ions_1 = set(flatten(merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()))
                    neutral_pool.remove(neutral_1)
                    for neutral_2 in neutral_pool:
                        ions_2 = set(flatten(merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()))
                        tmp_contested_ions = list(ions_1.intersection(ions_2))
                        for ion in tmp_contested_ions :
                            if ion not in contested_ions.index:
                                contested_ions.loc[ion, "conflicted_neutrals"] = [neutral_1, neutral_2]
                                contested_ions.loc[ion, "conflicted_neutrals"].sort()
                            else:
                                contested_ions.loc[ion, "conflicted_neutrals"] += [neutral_1, neutral_2]
                                contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]))
                                contested_ions.loc[ion, "conflicted_neutrals"].sort()
                                
                # Get unresolved neutrals (neutrals with only contested ions)
                conflicted_neutrals = list()
                for neutral in merged_neutral_table.index:
                    resolved_ions = set(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                    resolved_ions = resolved_ions - set(contested_ions.index)
                    if len(resolved_ions) == 0 : conflicted_neutrals.append(neutral)
                    
                # ???
                drop_list = list()
                for ion in contested_ions.index:
                    contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]) - set(conflicted_neutrals))
                    if len(contested_ions.loc[ion, "conflicted_neutrals"]) <= 1 :
                        drop_list.append(ion)
                contested_ions.drop(drop_list, inplace = True)
                
                # For each conflicted ion, do ion affinity:
                same_group_list = list()
                migration_list = list()
                resolved_list = list()
                for ion_1 in contested_ions.index:
                    ion_1_spec_id = merged_node_table.loc[ion_1, "spec_id"]
                    affinity_table = list()
                    resolved = True
                    for neutral in contested_ions.loc[ion_1, "conflicted_neutrals"]:

                        for i in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                            if ion_1 in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']][i] : 
                                ion_1_group = i
                        ion_1_group = adduct_table_merged['Group'][adduct_table_merged['Adduct_code'] == ion_1_group].iloc[0]
                        conflicted_score = False
                        neutral_ions = set(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                        neutral_ions = list(neutral_ions - set(contested_ions.index))
                        if len(neutral_ions) == 0 :
                            neutral_ions = list(flatten(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()))
                            neutral_ions.remove(ion_1)
                            conflicted_score = True
                        score_table = list()
                        # Get best cosine score with associated n_matches
                        for ion_2 in neutral_ions:
                            for i in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                                if ion_2 in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']][i] : 
                                    ion_2_group = i
                            ion_2_group = adduct_table_merged['Group'][adduct_table_merged['Adduct_code'] == ion_2_group].iloc[0]
              
                            # cosine, shared peaks, rule_points, dmz, drt, cohort_size
                            ion_2_spec_id = merged_node_table.loc[ion_2, "spec_id"]
                            score, n_matches = modified_cosine.pair(spectrum_list[ion_1_spec_id],
                                                                    spectrum_list[ion_2_spec_id])
                            if score < cosine_threshold : score = 0.0
                            prod = score * n_matches
                            if ion_1_group == ion_2_group :
                                same_group = True
                            else:
                                same_group = False
                            score_table.append((score, n_matches, prod, same_group))
                        score_table = pd.DataFrame(score_table, columns = ['cos', 'matches', 'prod', 'same_group'])
                        score_table.sort_values("prod", ascending = False, inplace = True)
                        
                        # Get d_mz and d_rt from neutral
                        annotation = merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()
                        ion_count = len(annotation)
                        annotation = [i for i in annotation.index if ion_1 in annotation[i]][0]
                        add_idx = adduct_table_merged.index[adduct_table_merged['Adduct_code'] == annotation][0]
                        ion_neutral = neutral_mass_calculator(merged_node_table.loc[ion_1, mz_field],
                                                adduct_table_merged.loc[add_idx, "Adduct_mass"],
                                                adduct_table_merged.loc[add_idx, "Mol_multiplier"],
                                                adduct_table_merged.loc[add_idx, "Charge"])
                        d_mz = abs(ion_neutral - merged_neutral_table.loc[neutral, "mass"])
                        d_rt = abs(merged_node_table.loc[ion_1, rt_field] - merged_neutral_table.loc[neutral, "rt"])
                        Species_rule = Validator_choice(annotation, ion_mode)
                        rule_point = Species_rule(prec_mass_error, ion_1_spec_id, spectrum_list)
                        if len(score_table) != score_table['same_group'].sum():
                            same_group = False
                        else:
                            same_group = True
                        affinity_table.append((neutral,
                                               score_table['cos'].iloc[0],
                                               score_table['matches'].iloc[0],
                                               rule_point,
                                               (score_table['cos'].iloc[0] * score_table['matches'].iloc[0]) + rule_point,
                                               d_mz,
                                               d_rt,
                                               d_mz * d_rt,
                                               ion_count,
                                               conflicted_score,
                                               same_group))
                    affinity_table = pd.DataFrame(affinity_table, columns = ['neutral', 'cos',
                                                                             'matches', 'rule_points', 'score_1',
                                                                             'd_mz', 'd_rt', 'score_2', 'ion_count',
                                                                             'conflicted_score', 'same_group'])
                    affinity_table.sort_values(by = ['score_1', 'score_2', 'ion_count'], ascending = [False, True, False], inplace = True)
                    max_score = affinity_table['score_1'].max()
                    if sum(affinity_table['score_1'].between(max_score - 0.05, max_score + 0.05, inclusive = "both")) > 1 :
                        resolved = False
                    migration_list.append(list(affinity_table['neutral']))
                    resolved_list.append(resolved)
                    if len(affinity_table) != affinity_table['same_group'].sum():
                        same_group_list.append(False)
                    else:
                        same_group_list.append(True)
                contested_ions['selected_neutrals'] = migration_list
                contested_ions['resolved'] = resolved_list
                contested_ions['same_group'] = same_group_list
                contested_ions = contested_ions[contested_ions['same_group']]
                neutral_ion_count = pd.DataFrame(columns = ['ion_count'])
                for i in merged_neutral_table.index:
                    neutral_ion_count.loc[i, 'ion_count'] = len(list(flatten(merged_neutral_table.loc[i, adduct_table_merged['Adduct_code']].dropna())))
                neutral_ion_count["ion_count"] = neutral_ion_count["ion_count"].astype(int)
                
                eliminating_neutrals = True
                while eliminating_neutrals:                
                    neutral_ion_count["new_ion_count"] = neutral_ion_count["ion_count"].copy()
                    for i in contested_ions.index:
                        if contested_ions.loc[i, "resolved"]:
                            loosing_neutrals = contested_ions.loc[i, "selected_neutrals"][1:]
                            for j in loosing_neutrals:
                                neutral_ion_count.loc[j, "new_ion_count"] = neutral_ion_count.loc[j, "ion_count"] - 1
                    dead_neutrals = neutral_ion_count.index[neutral_ion_count['new_ion_count'] <= 1]
                    if len(dead_neutrals) > 0 :
                        for i in dead_neutrals:
                            neutral_ion_count.drop(i, inplace = True)
                            merged_neutral_table.drop(i, inplace = True)
                            for j in contested_ions.index:
                                if i in contested_ions.loc[j, "conflicted_neutrals"]:
                                    contested_ions.loc[j, "conflicted_neutrals"].remove(i)
                                if i in contested_ions.loc[j, "selected_neutrals"] :
                                    contested_ions.loc[j, "selected_neutrals"].remove(i)
                        drop_list = []
                        for i in contested_ions.index:
                            if len(contested_ions.loc[i, "selected_neutrals"]) == 0 :
                                drop_list.append(i)
                            if len(contested_ions.loc[i, "selected_neutrals"]) == 1 :
                                contested_ions.loc[i, "resolved"] = True
                        contested_ions.drop(drop_list, inplace = True)
                    else:
                        eliminating_neutrals = False
                #Once this is done, do the migration:
                for i in contested_ions.index:
                    for j in contested_ions.loc[i, "selected_neutrals"][1:]:
                        for k in adduct_table_merged['Adduct_code']:
                            if merged_neutral_table.loc[j, k] != None and i in merged_neutral_table.loc[j, k]:
                                merged_neutral_table.loc[j, k].remove(i)
                                if len(merged_neutral_table.loc[j, k]) == 0:
                                    merged_neutral_table.loc[j, k] = None
######################################################## ION AFFINITY ENDS HERE
###############################################################################                
                
                # Check neutral compatibility:
                merged_neutral_table = neutral_compatibility(merged_neutral_table, adduct_table_merged)
                
                # Eject inert neutrals:
                inert_neutrals = list()
                neutral_count = len(merged_neutral_table.index)
                for i in merged_neutral_table.index:
                    inert_neutrals.append((i, len(merged_neutral_table.loc[i, "compatibles"]) + 1))
                inert_neutrals = pd.DataFrame(inert_neutrals, columns = ['neutral', 'compatible_count'])
                inert_neutrals = inert_neutrals[inert_neutrals['compatible_count'] == neutral_count]
                if len(inert_neutrals) > 0 : 
                    merged_neutral_table_selected = merged_neutral_table_selected.append(merged_neutral_table.loc[inert_neutrals['neutral']])
                    merged_neutral_table.drop(inert_neutrals['neutral'], inplace = True)

                # Check neutral compatibility:
                merged_neutral_table = neutral_compatibility(merged_neutral_table, adduct_table_merged)

                if len(merged_neutral_table) == 0 :
                    selected_house = list(merged_neutral_table_selected.index)
                    end_selection = True
                else:
                    
                    # Detect uncompatibles (Neutrals compatible with nothing else):
                    uncompatibles = []
                    for i in merged_neutral_table.index:
                        if len(merged_neutral_table.loc[i, "compatibles"]) == 0:
                            uncompatibles.append([i])   
            
                    # Get houses (combinations of compatible neutrals from the same court)
                    houses = []
                    neutral_pool = list(merged_neutral_table.index)
                    tmp_merged_neutral_table = merged_neutral_table.copy()
                    while len(neutral_pool) > 0:
                        neutral_seed = neutral_pool[0]
                        #for neutral_seed in merged_neutral_table.index:
                        next_list, house_list, compatibles_list, left_list = initialise_house_table_neutrals(neutral_seed, tmp_merged_neutral_table)       
                        if len(house_list) == 0 : houses.append([neutral_seed])
                        while len(house_list) > 0 :
                            next_list, house_list, compatibles_list, left_list, houses = clean_house_table_neutrals(next_list, house_list, compatibles_list, left_list, houses)
                            next_list, house_list, compatibles_list, left_list = expand_house_table_neutrals(next_list, house_list, compatibles_list, left_list, tmp_merged_neutral_table)
                        tmp_merged_neutral_table.drop(neutral_seed, inplace = True)
                        for i in tmp_merged_neutral_table.index:
                            if neutral_seed in tmp_merged_neutral_table.loc[i, "compatibles"]:
                                tmp_merged_neutral_table.loc[i, "compatibles"].remove(neutral_seed)
                        neutral_pool = list(tmp_merged_neutral_table.index)
                    houses += uncompatibles
                    
                    # Create house table for the neutral combinations and assign scores to houses
                    tmp_points_list = list()
                    for i in range(len(houses)):
                        tmp_points = 0.0
                        house = houses[i]
                        for j in house:
                            tmp_points += merged_neutral_table.loc[j, "points"]
                        tmp_points_list.append(tmp_points)
                    house_table = pd.DataFrame()
                    house_table['house'] = houses
                    house_table['points'] = tmp_points_list
                    
                    # Select the best scoring house and report annotation results:
                    selected_house = house_table['points'].idxmax()
                    selected_house = house_table.loc[selected_house, "house"]
                    if len(selected_house) == len(merged_neutral_table.index) :
                        end_selection = True
                    else: 
                        # Ion affinity on neutral remains:
                        # Find conflicted ions
                        contested_ions = pd.DataFrame(columns = ["conflicted_neutrals"])
                        neutral_pool = list(merged_neutral_table.index)
                        while len(neutral_pool) > 1 :
                            neutral_1 = neutral_pool[0]
                            ions_1 = set(flatten(merged_neutral_table.loc[neutral_1, adduct_table_merged['Adduct_code']].dropna()))
                            neutral_pool.remove(neutral_1)
                            for neutral_2 in neutral_pool:
                                ions_2 = set(flatten(merged_neutral_table.loc[neutral_2, adduct_table_merged['Adduct_code']].dropna()))
                                tmp_contested_ions = list(ions_1.intersection(ions_2))
                                for ion in tmp_contested_ions :
                                    if ion not in contested_ions.index:
                                        contested_ions.loc[ion, "conflicted_neutrals"] = [neutral_1, neutral_2]
                                        contested_ions.loc[ion, "conflicted_neutrals"].sort()
                                    else:
                                        contested_ions.loc[ion, "conflicted_neutrals"] += [neutral_1, neutral_2]
                                        contested_ions.loc[ion, "conflicted_neutrals"] = list(set(contested_ions.loc[ion, "conflicted_neutrals"]))
                                        contested_ions.loc[ion, "conflicted_neutrals"].sort()
                        
                        winning_neutral = list()
                        loosing_neutrals = list()
                        for i in contested_ions.index:
                            tmp_win = list(set(selected_house).intersection(contested_ions.loc[i, "conflicted_neutrals"]))
                            tmp_loss = list(set(contested_ions.loc[i, "conflicted_neutrals"]) - set(tmp_win))
                            winning_neutral.append(tmp_win)
                            loosing_neutrals.append(tmp_loss)
                        contested_ions['winning_neutral'] = winning_neutral
                        contested_ions['loosing_neutrals'] = loosing_neutrals
                        for i in contested_ions.index:
                            for neutral in contested_ions.loc[i, "loosing_neutrals"]:
                                for adduct in merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna().index:
                                    merged_neutral_table.loc[neutral, adduct] = list(set(merged_neutral_table.loc[neutral, adduct]) - set([i]))
                                    if len(merged_neutral_table.loc[neutral, adduct]) == 0 :
                                        merged_neutral_table.loc[neutral, adduct] = None
                        dead_neutrals = list()
                        for neutral in merged_neutral_table.index:
                            if len(merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()) <= 1 :
                                dead_neutrals.append(neutral)
                        merged_neutral_table.drop(dead_neutrals, inplace = True)

                                   
        # Report results
        merged_neutral_table = merged_neutral_table_selected
        for neutral in selected_house:
            neutral_mass = merged_neutral_table.loc[neutral, "mass"]
            neutral_rt = merged_neutral_table.loc[neutral, "rt"]
            neutral_idx = merged_node_table.index.max() + 1
            merged_node_table.loc[neutral_idx] = [None]*len(merged_node_table.columns)
            merged_node_table.loc[neutral_idx, mz_field] = neutral_mass
            merged_node_table.loc[neutral_idx, rt_field] = neutral_rt
            merged_node_table.loc[neutral_idx, 'status'] = "neutral"
            merged_node_table.loc[neutral_idx, charge_field] = 0
            merged_node_table.loc[neutral_idx, 'TIC'] = 0.0
            annotations = merged_neutral_table.loc[neutral, adduct_table_merged['Adduct_code']].dropna()
            for annotation in annotations.index:
                ion_list = annotations[annotation]
                for ion in ion_list:
                    mz_gap = round(merged_node_table.loc[ion, mz_field] - neutral_mass, 4)
                    new_edge = merged_edge_table.index.max() + 1
                    merged_edge_table.loc[new_edge] = [neutral_idx, ion, 0, 0, 0.0, 0.0, mz_gap, "add_edge", None, annotation, annotation]
                    merged_node_table.loc[neutral_idx, 'TIC'] += merged_node_table.loc[ion, 'TIC']
                    if merged_node_table.loc[ion, 'Adnotation'] == None :
                        merged_node_table.loc[ion, 'Adnotation'] = annotation
                    else :
                        sys.exit("ERROR HERE FOR SOME REASON : already adnotated ion?")
    return merged_node_table, merged_edge_table

def update_node_table(merged_node_table, merged_edge_table, cross_annotations,
                      cross_rules, adduct_table_merged, params):
    """Updates the node table with the "adduct" status on adnotated nodes.
    Also adds an "adduct_count" columns with the number of ions directly
    connected to the neutral node.
    """
    # Load parameters 
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    
    # Add rule points to the node table
    print("Adding rule points...")
    rule_points = list()
    for i in tqdm(merged_node_table.index):
        adnotation = merged_node_table.loc[i, 'Adnotation']
        if adnotation == None : 
            rule_points.append(0)
        else:
            sample = (cross_annotations.loc[i].dropna() == adnotation).index[0]
            rule_points.append(int(cross_rules.loc[i, sample]))
    merged_node_table.insert(merged_node_table.columns.get_loc("spec_id") + 1, 'rule_points', rule_points)

    # Updated the node table status
    adnoted_ions = merged_node_table['Adnotation'].dropna().index
    print("Updating node table status")
    for i in tqdm(adnoted_ions):
        merged_node_table.loc[i, "status"] = "adduct"
    
    # Add adduct count for neutrals
    merged_node_table['adduct_count'] = [0]*len(merged_node_table)
    for i in merged_node_table.index[merged_node_table['status'] == "neutral"]:
        count = sum(merged_edge_table['node_1'] == i)
        merged_node_table.loc[i, "adduct_count"] = count
    
    # Change adduct codes to adduct strings
    for adduct in tqdm(merged_node_table['Adnotation'].dropna().unique()):
        clean_adduct = adduct_table_merged['Adduct'][adduct_table_merged['Adduct_code'] == adduct].iloc[0]
        change_idx = merged_node_table.index[merged_node_table['Adnotation'] == adduct]
        merged_node_table.loc[change_idx, "Adnotation"] = clean_adduct
            
    # Round values
    merged_node_table[mz_field] = merged_node_table[mz_field].round(4)
    merged_node_table[rt_field] = merged_node_table[rt_field].round(3)
    
    return merged_node_table

def update_edge_table(merged_node_table, merged_edge_table, adduct_table_merged):
    """Updates edge table with singleton nodes
    """
    print("Updating edge table...")
    singleton_ions = merged_node_table.index[merged_node_table["status"] == "singleton"]
    singleton_df = pd.DataFrame(columns = merged_edge_table.columns)
    singleton_df['node_1'] = list(singleton_ions)
    singleton_df['node_2'] = list(singleton_ions)
    singleton_df['matched_peaks'] = [0]*len(singleton_ions)
    singleton_df['total_peaks'] = [0]*len(singleton_ions)
    singleton_df['matching_score'] = [0.0]*len(singleton_ions)
    singleton_df['rt_gap'] = [0.0]*len(singleton_ions)
    singleton_df['mz_gap'] = [0.0]*len(singleton_ions)
    singleton_df['status'] = ["self_edge"]*len(singleton_ions)
    singleton_df['Fragnotation'] = [None]*len(singleton_ions)
    singleton_df['All_annotations'] = [None]*len(singleton_ions)
    singleton_df['Adnotation'] = [None]*len(singleton_ions)
    merged_edge_table = merged_edge_table.append(singleton_df, ignore_index = True)

    print("Replacing adduct codes by adduct formulas...")
    for adduct in tqdm(merged_edge_table['Adnotation'].dropna().unique()):
        clean_adduct = adduct_table_merged['Adduct'][adduct_table_merged['Adduct_code'] == adduct].iloc[0]
        change_idx = merged_edge_table.index[merged_edge_table['Adnotation'] == adduct]
        merged_edge_table.loc[change_idx, "Adnotation"] = clean_adduct
        merged_edge_table.loc[change_idx, "All_annotations"] = clean_adduct

    # Round and correct value types
    merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].astype(float)
    merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].round(4)
    merged_edge_table['rt_gap'] = merged_edge_table['rt_gap'].astype(float)
    merged_edge_table['rt_gap'] = merged_edge_table['rt_gap'].round(3)

    # Change all annotations position
    all_annotations = merged_edge_table["All_annotations"].copy()
    merged_edge_table.drop('All_annotations', axis = 1, inplace = True)
    merged_edge_table['All_annotations'] = all_annotations
    return merged_edge_table






######### MODE MERGER FUNCTIONS

def get_ion_table(neutral_mass : float, adduct_table):
    """
    Computes all possible ion m/z values for an given molecular mass
    given the different ion species available in the adducts table.

    Parameters
    ----------
    neutral_mass : float
        Mass of the neutral node.
    adduct_table : pandas.DataFrame
        Adduct table loaded from the parameters folder

    Returns
    -------
    ion_table : pandas.DataFrame
        Ions that could be generated by the neutral node.

    """
    
    ion_table = adduct_table.copy()
    mz_list = []
    for i in ion_table.index:
        mz_list.append(ion_mass_calculator(neutral_mass,
                                         ion_table.loc[i, "Adduct_mass"],
                                         ion_table.loc[i, "Mol_multiplier"],
                                         ion_table.loc[i, "Charge"]))
    ion_table['ion_mz'] = mz_list
    return ion_table




# Dereplicator functions

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

def Filter_choices(db_params, database_mgf):
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

def Database_table_mgf(database_mgf, db_params : dict) :
    """Extract the metadata from the MGF file
    """
    def Mz_field_choices(s, database_mgf, db_params):
        def Get_mz_list(s, database_mgf, db_params):
            return database_mgf[s].get(db_params['db_mass_field'])[0]
        def Get_mz(s, database_mgf, db_params):
            return database_mgf[s].get(db_params['db_mass_field'])
        if isinstance(database_mgf[0].get(db_params['db_mass_field']), list) : 
            return Get_mz_list
        else:
            return Get_mz

    db_unique_field = db_params['db_unique_field']
    db_mode_field= db_params['db_mode_field']
    database_table = list()
    Filter_fields, cols = Filter_choices(db_params, database_mgf)
    Mz_extractor = Mz_field_choices(0, database_mgf, db_params)
    col_names = ['name', 'mz', 'unique_field'] + cols + db_params['db_export_fields'] + ['ion_mode']
    print('Extracting database metadata...')
    for i in tqdm(range(len(database_mgf))):
        name = database_mgf[i].get(db_params['db_name_field'])
        mz = Mz_extractor(i, database_mgf, db_params)
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



def Float_prec_mz(s, db_params):
    s = s.set(db_params['db_mass_field'], float(s.get(db_params['db_mass_field'])))
    return s

#################################################### Cosiner functions
def cosiner_single(node_table, edge_table, mgf, mgf_data, ion_mode, params):

    idx_column = params['index_col']
    modified_cosine = ModifiedCosine(tolerance=params['c_mass_error'])
    cosiner_threshold= params['c_hardcos_threshold']
    matched_peaks = params['c_matched_peaks']
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    cosine_threshold= params['c_lowcos_threshold']
    out_path_full = params['mix_out_6_1']

    # List the molecular clusters (clusters with at least one neutral node)
    cluster_list = []
    print('Finding molecular clusters...')
    for i in tqdm(node_table['cluster_id'].unique()):
        if sum(node_table['status_universal'][node_table['cluster_id'] == i] == "neutral") > 0:
            cluster_list.append(i)
    cluster_list.sort()

    # Cluster singletons and precursor ions in non-molecular clusters to ions in molecular clusters
    unclustered_ions = [i for i in node_table.index if node_table.loc[i, "cluster_id"] not in cluster_list]
    remains_ions = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status_universal'] == 'singleton'])
    remains_ions +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status_universal'] == 'precursor'])

    # Get spec IDs for nodes
    cluster_ion_list = pd.Series(index = cluster_list, dtype = str)
    for i in cluster_ion_list.index:
        tmp_rows = node_table.index[node_table['cluster_id'] == i]
        cluster_ion_list[i] = '|'.join(node_table.loc[tmp_rows, 'spec_id'].dropna().astype(int).astype(str))
    cluster_ion_list.sort_index(inplace = True)

    if len(cluster_list) > 0 :
        full_node_1 = []
        full_node_2 = []
        full_cluster_ids = []
        full_cosine = []
        full_matches = []
        for i in tqdm(remains_ions):
            i_id = int(node_table.loc[i, idx_column])
            i_spectrum = mgf[mgf_data[i_id]]
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
                    score, n_matches = modified_cosine.pair(i_spectrum, mgf[k])
                    id_list.append(int(mgf[k].get(idx_column)))
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
            counter_idx = node_table.index[node_table[idx_column] == counter_feature_id][0]
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
            if (edge_table.loc[old_edge, "status_universal"] == "self_edge") :
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
            node_2_list = list(tmp_edge_table['node_2'])
            tmp_edges = []
            while len(node_2_list) > 1 :
                ion_1 = node_2_list[0]
                ion_1_mgf_idx = int(node_table.loc[ion_1, "spec_id"])
                ion_1_spectrum = mgf[ion_1_mgf_idx]
                ion_1_rt = node_table.loc[ion_1, rt_field]
                ion_1_mz = node_table.loc[ion_1, mz_field]
                node_2_list.remove(ion_1)
                for ion_2 in node_2_list:
                    ion_2_mgf_idx = int(node_table.loc[ion_2, "spec_id"])
                    ion_2_spectrum = mgf[ion_2_mgf_idx]
                    ion_2_rt = node_table.loc[ion_2, rt_field]
                    ion_2_mz = node_table.loc[ion_2, mz_field]
                    rt_gap = round(ion_1_rt - ion_2_rt, 3)
                    mz_gap = round(ion_1_mz - ion_2_mz, 4)
                    score, n_matches = modified_cosine.pair(ion_1_spectrum, ion_2_spectrum)
                    tmp_edges.append((ion_1, ion_2, round(score, 2), rt_gap, mz_gap, ion_mode.lower()))
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
    remains_ions = list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status_universal'] == 'singleton'])
    remains_ions +=  list(node_table.loc[unclustered_ions].index[node_table.loc[unclustered_ions]['status_universal'] == 'precursor'])
    remains_ions_mgf = [node_table.loc[i, "spec_id"] for i in remains_ions] 
    remains_ions_mgf = list(map(int, remains_ions_mgf))

    # Process data:
    singleton_clusters = list()
    total_nodes = len(remains_ions)
    while len(remains_ions) > 0 :
        perc = round((1-(len(remains_ions)/total_nodes))*100,1)
        sys.stdout.write("\rClustering remaining singletons : {0}%".format(perc))
        sys.stdout.flush()
        ion_i = remains_ions[0]
        ion_i_mgf = remains_ions_mgf[0]
        remains_ions.remove(ion_i)
        remains_ions_mgf.remove(ion_i_mgf)
        for ion_j in remains_ions:
            ion_j_mgf = remains_ions_mgf[remains_ions.index(ion_j)]
            score, n_matches = modified_cosine.pair(mgf[ion_i_mgf], mgf[ion_j_mgf])
            singleton_clusters.append((ion_i, ion_j, score, n_matches, ion_mode))

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

    if params['c_purge_empty_spectra'] :
        node_table, edge_table = Spectrum_purge(ion_mode = ion_mode, node_table = node_table, edge_table = edge_table, mgf_file = mgf)

    node_table[mz_field] = node_table[mz_field].round(4)
    node_table[rt_field] =  node_table[rt_field].round(3)
    edge_table['rt_gap'] = edge_table['rt_gap'].round(3)
    edge_table['mz_gap'] = edge_table['mz_gap'].round(4)
    edge_table['cosine_score'] = edge_table['cosine_score'].round(2)
                                                     
    # Export data
    node_table.to_csv(out_path_full + 'node_table.csv', index_label = "Index")
    edge_table.to_csv(out_path_full + 'edge_table.csv', index_label = "Index")

    if params['c_export_samples'] : 
        samplewise_export(merged_edge_table = edge_table,
                          merged_node_table = node_table,
                          step = "cosiner",
                          ion_mode = ion_mode,
                          params = params)
    return

def Spectrum_purge(ion_mode, node_table, edge_table, mgf_file) :
    singletons = node_table.index[node_table['status'] == ion_mode.lower() + "_singleton"].tolist()
    empty_nodes = list()
    print(f"Purging empty {ion_mode} singletons...")
    for i in tqdm(singletons):
        mgf_idx = int(node_table.loc[i, "spec_id"])
        if len(mgf_file[mgf_idx].peaks.mz) == 0 :
            node_table.drop
            empty_nodes.append(i)
    del_edges = [edge_table.index[edge_table['node_1'] == n][0] for n in empty_nodes]
    node_table.drop(empty_nodes, inplace = True)
    edge_table.drop(del_edges, inplace = True)
    return node_table, edge_table

############################################################## molnet functions
def molnet_single(node_table, edge_table, mgf, ion_mode, params):
    
    mass_error= params['mn_mass_error']
    mz_field = params['mz_field']
    rt_field = params['rt_field']
    matched_peaks= params['mn_matched_peaks']
    modified_cosine = ModifiedCosine(tolerance=mass_error)
    cosine_threshold = params['mn_cosine_threshold']
    out_path_full = params['mix_out_7_1']
    
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
        adducts_1['spec_id'] = list(node_table.loc[adducts_1.index, "spec_id"].astype(int))
        adducts_1['ion_mode'] = list(node_table.loc[adducts_1.index, "ion_mode"])
        adducts_1['adduct'] = list(node_table.loc[adducts_1.index, "Adnotation"])
        neutral_nodes.remove(neutral_1)
        for neutral_2 in neutral_nodes:
            #neutral_2 = neutral_nodes[1]
            tmp_scores = list()
            adducts_2 = pd.DataFrame(index = list(edge_table['node_2'][edge_table['node_1'] == neutral_2]))
            adducts_2['spec_id'] = list(node_table.loc[adducts_2.index, "spec_id"].astype(int))
            adducts_2['ion_mode'] = list(node_table.loc[adducts_2.index, "ion_mode"])
            adducts_2['adduct'] = list(node_table.loc[adducts_2.index, "Adnotation"])
            for mode in adducts_1['ion_mode'].unique():
                #mode = adducts_1['ion_mode'].unique()[0]
                tmp_table_1 = adducts_1[adducts_1['ion_mode'] == mode]
                tmp_table_2 = adducts_2[adducts_2['ion_mode'] == mode]
                for ion_1 in tmp_table_1.index:
                    spectrum_1 = mgf[tmp_table_1.loc[ion_1, "spec_id"]]
                    adduct_1 = tmp_table_1.loc[ion_1, 'adduct']
                    for ion_2 in tmp_table_2.index:
                        adduct_2 = tmp_table_2.loc[ion_2, 'adduct']
                        if adduct_1 != adduct_2 : continue
                        spectrum_2 = mgf[tmp_table_2.loc[ion_2, "spec_id"]]
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
        rt_gap = abs(round(node_table.loc[node_1, rt_field] - node_table.loc[node_2, rt_field], 3))
        mz_gap = abs(round(node_table.loc[node_1, mz_field] - node_table.loc[node_2, mz_field], 4))
        edge_table.loc[new_idx] = [None]*len(edge_table.columns)
        edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks", "total_peaks",
                                 "matching_score", "rt_gap", "mz_gap", "status",
                                 "status_universal", "All_annotations", "ion_mode",
                                 "cosine_score"]] = [node_1, node_2, matches, 0, 0.0, rt_gap, mz_gap,
                      ion_mode + "_cosine_neutral_edge", "cosine_neutral_edge", cos, ion_mode, cos]

    # Delete undesired edges
    edge_table = edge_table[edge_table['status'].str.contains('cosine_neutral_edge')]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    
    # Add singletons:
    singletons = list(set(neutral_nodes_2) - set(kept_nodes))
    for neutral in tqdm(singletons):
        new_idx = max(edge_table.index) + 1 
        ion_mode = node_table.loc[neutral, "ion_mode"]
        edge_table.loc[new_idx] = [None]*len(edge_table.columns)
        edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks", "total_peaks",
                                 "matching_score", "rt_gap", "mz_gap", "status",
                                 "status_universal", "All_annotations", "ion_mode",
                                 "cosine_score"]] = [neutral, neutral, 0, 0, 0.0, 0.0, 0.0,
                      ion_mode + "_self_edge", "self_edge", 0.0, ion_mode, 0.0]
    
    kept_nodes = list(set(list(edge_table['node_1']) + list(edge_table['node_2'])))
    kept_nodes.sort()
    node_table = node_table.loc[kept_nodes]
    
    node_table[mz_field] = node_table[mz_field].round(4)
    
    node_table.to_csv(out_path_full + "node_table_families.csv", index_label = "Index")
    edge_table.to_csv(out_path_full + "edge_table_families.csv", index_label = "Index")
    return
                                                     
                                                     
#################################################### Sample export functions

def samplewise_export(merged_edge_table, merged_node_table, step, ion_mode, params) : 

    col_suffix = params['col_suffix']
    if step == "adnotator":
        idx_column = params['index_col']
        status_col = "status"
        out_idx = idx_column
        if ion_mode == "POS" :
            in_path = params['pos_out_0']
            out_path_samples = params['pos_out_3_2']
            csv_file = params['pos_csv']
        else:
            in_path = params['neg_out_0']
            out_path_samples = params['neg_out_3_2']
            csv_file = params['neg_csv']
    elif step == "dereplicator" :
        idx_column = params['index_col']
        status_col = "status_universal"
        out_idx = "Index"        
        if ion_mode == "POS":
            in_path = params['pos_out_0']
            csv_file = params['pos_csv']
            out_path_samples = params['mix_out_5_2']
        elif ion_mode == "NEG":
            in_path = params['pos_out_0']
            csv_file = params['pos_csv']
            out_path_samples = params['mix_out_5_2']      
    elif step == "cosiner" :
        idx_column = params['index_col']
        status_col = "status_universal"
        out_idx = "Index"
        if ion_mode == "POS":
            in_path = params['pos_out_0']
            csv_file = params['pos_csv']
            out_path_samples = params['mix_out_6_2']
        elif ion_mode == "NEG":
            in_path = params['pos_out_0']
            csv_file = params['pos_csv']
            out_path_samples = params['mix_out_6_2']
            

    print("Exporting sample-wise tables...")
    full_csv = pd.read_csv(f"{in_path}{csv_file}", index_col = idx_column)
    samples = full_csv.columns
    samples = samples[samples.str.contains(col_suffix)]
    samples = list(samples.str.replace(col_suffix, "", regex = False))
    full_csv.columns = full_csv.columns.str.replace(col_suffix, "", regex = False)

    if not os.path.isdir(out_path_samples):
        os.mkdir(out_path_samples)
    
    for sample in tqdm(samples):
        # get sample ions
        ion_ids = full_csv.index[full_csv[sample] > 0.0]
    
        # Get sample neutrals
        kept_edges = list()
        neutral_edges = merged_edge_table[merged_edge_table[status_col] == "add_edge"]
        for i in neutral_edges.index:
            if neutral_edges.loc[i, "node_2"] in ion_ids :
                kept_edges.append(i)
    
        # Get ion edges
        ion_edges = merged_edge_table[merged_edge_table[status_col] != "add_edge"]
        for i in ion_edges.index:
            if ion_edges.loc[i, "node_1"] in ion_ids:
                if ion_edges.loc[i, "node_2"] in ion_ids:
                    kept_edges.append(i)
        kept_edges.sort()
        sample_edges = merged_edge_table.loc[kept_edges]
        sample_edges.sort_values('node_1', inplace = True)
        sample_edges.reset_index(inplace = True, drop = True)
        
        kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
        kept_nodes = list(set(kept_nodes + list(ion_ids)))
        kept_nodes.sort()
        sample_nodes = merged_node_table.loc[kept_nodes]
        
        # Add sample peak intensities / areas to ions
        sample_nodes[sample] = [0.0]*len(sample_nodes)
        sample_nodes.loc[ion_ids, sample] = full_csv.loc[ion_ids, sample]
        
        # Add sample peak intensitites / areas to neutrals
        sample_neutrals = list(set(sample_nodes.index) - set(ion_ids))
        for neutral in sample_neutrals:
            tmp_int = list(sample_edges['node_2'][sample_edges['node_1'] == neutral])
            sample_nodes.loc[neutral, sample] = sample_nodes.loc[tmp_int, sample].sum()
        
        sample_nodes.to_csv(out_path_samples + sample+"_node_table.csv", index_label = out_idx)
        sample_edges.to_csv(out_path_samples + sample+"_edge_table.csv", index_label = "Index")
    return

def mm_samplewise_export(neg_csv_file, pos_csv_file, out_path, edge_table, node_table, params) : 
    idx_column = params['index_col']
    col_suffix = params['col_suffix']

    print("Exporting sample-wise tables...")
    
    # Load full attribute tables
    neg_csv = pd.read_csv(neg_csv_file, index_col = idx_column)
    pos_csv = pd.read_csv(pos_csv_file, index_col = idx_column)
    
    # Coerce attribute tables indexes to int
    neg_csv.index = neg_csv.index.astype(int)
    pos_csv.index = pos_csv.index.astype(int)

    # Filter full attribute tables
    neg_csv = neg_csv.loc[:, neg_csv.columns.str.contains(col_suffix)]
    pos_csv = pos_csv.loc[:, pos_csv.columns.str.contains(col_suffix)]    
    neg_csv.columns = neg_csv.columns.str.replace(col_suffix, "", regex = False).str.replace('NEG_', '', regex = False)
    pos_csv.columns = pos_csv.columns.str.replace(col_suffix, "", regex = False).str.replace('POS_', '', regex = False)
    samples = list(set(list(neg_csv.columns) + list(pos_csv.columns)))
    samples.sort()
    
    for sample in tqdm(samples):
        #sample = samples[0]
        ion_ids_neg = neg_csv.index[neg_csv[sample] > 0.0]
        ion_ids_pos = pos_csv.index[pos_csv[sample] > 0.0]
        
        #convert CSV IDs to the new indexes
        tmp_table = node_table[node_table['status'] != "neg_neutral"]
        tmp_table = tmp_table[tmp_table['status'] != "pos_neutral"]
        tmp_table = tmp_table[tmp_table['status'] != "mix_neutral"]
        tmp_table_pos = tmp_table[tmp_table['ion_mode'] == "POS"]
        tmp_table_neg = tmp_table[tmp_table['ion_mode'] == "NEG"]
        ion_idx_neg = pd.Series(tmp_table_neg.index, index = tmp_table_neg[idx_column])
        ion_idx_neg = list(ion_idx_neg[ion_ids_neg])
        ion_idx_pos = pd.Series(tmp_table_pos.index, index = tmp_table_pos[idx_column])
        ion_idx_pos = list(ion_idx_pos[ion_ids_pos])
        ion_idx_mix = ion_idx_neg + ion_idx_pos
    
        # Get sample neutrals
        neutral_edges = edge_table.loc[edge_table["Adnotation"].dropna().index]
        kept_edges = [i for i in neutral_edges.index if neutral_edges.loc[i, "node_2"] in ion_idx_mix]

    
        # Get ion edges
        ion_edges = edge_table[edge_table['status'] != "neg_add_edge"]
        ion_edges = ion_edges[ion_edges['status'] != "pos_add_edge"]
        for i in ion_edges.index:
            if ion_edges.loc[i, "node_1"] in ion_idx_mix:
                if ion_edges.loc[i, "node_2"] in ion_idx_mix:
                    kept_edges.append(i)
        kept_edges.sort()
        sample_edges = edge_table.loc[kept_edges]
        sample_edges.sort_values('node_1', inplace = True)
        sample_edges.reset_index(inplace = True, drop = True)
        
        kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
        kept_nodes.sort()
        sample_nodes = node_table.loc[kept_nodes].copy()
        sample_nodes.drop(pd.Series(samples), axis = 1, inplace = True)
        sample_nodes[sample] = node_table[sample]
        
        sample_nodes.to_csv(out_path + sample + "_node_table.csv", index_label = "Index")
        sample_edges.to_csv(out_path + sample + "_edge_table.csv", index_label = "Index")
    return


def dr_samplewise_export(neg_csv_file, pos_csv_file, out_path, edge_table, node_table, params) : 
    print("Exporting sample-wise tables...")
    idx_column = params['index_col']
    col_suffix = params['col_suffix']
    
    neg_csv = pd.read_csv(neg_csv_file, index_col = idx_column)
    pos_csv = pd.read_csv(pos_csv_file, index_col = idx_column)
    
    # Coerce attribute tables indexes to int
    neg_csv.index = neg_csv.index.astype(int)
    pos_csv.index = pos_csv.index.astype(int)
    
    neg_csv = neg_csv.loc[:, neg_csv.columns.str.contains(col_suffix, regex = False)]
    pos_csv = pos_csv.loc[:, pos_csv.columns.str.contains(col_suffix, regex = False)]
    
    neg_csv.columns = neg_csv.columns.str.replace(col_suffix, "", regex = False).str.replace('NEG_', '', regex = False)
    pos_csv.columns = pos_csv.columns.str.replace(col_suffix, "", regex = False).str.replace('POS_', '', regex = False)

    samples = list(set(list(neg_csv.columns) + list(pos_csv.columns)))
    samples.sort()
    
    for sample in tqdm(samples):
        #sample = samples[0]
        ion_ids_neg = neg_csv.index[neg_csv[sample] > 0.0]
        ion_ids_pos = pos_csv.index[pos_csv[sample] > 0.0]
        
        #convert feature_ids to the new indexes
        if params['process_mode'] == "BOTH" :
            tmp_table = node_table[node_table['status'] != "neg_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "pos_neutral"]
            tmp_table = tmp_table[tmp_table['status'] != "mix_neutral"]
            tmp_table_pos = tmp_table[tmp_table['ion_mode'] == "POS"]
            tmp_table_neg = tmp_table[tmp_table['ion_mode'] == "NEG"]
            ion_idx_neg = pd.Series(tmp_table_neg.index, index = tmp_table_neg[idx_column])
            ion_idx_neg = list(ion_idx_neg[ion_ids_neg])
            ion_idx_pos = pd.Series(tmp_table_pos.index, index = tmp_table_pos[idx_column])
            ion_idx_pos = list(ion_idx_pos[ion_ids_pos])
            ion_idx_mix = ion_idx_neg + ion_idx_pos
        elif params['process_mode'] == "POS":
            tmp_table_pos = node_table[node_table['status'] != "pos_neutral"]
            ion_idx_pos = pd.Series(tmp_table_pos.index, index = tmp_table_pos[idx_column])
            ion_idx_pos = list(ion_idx_pos[ion_ids_pos])
            ion_idx_mix = ion_idx_pos            
        elif params['process_mode'] == "NEG":
            tmp_table_neg = node_table[node_table['status'] != "neg_neutral"]
            ion_idx_neg = pd.Series(tmp_table_neg.index, index = tmp_table_neg[idx_column])
            ion_idx_neg = list(ion_idx_neg[ion_ids_neg])
            ion_idx_mix = ion_idx_neg
    
        # Get sample neutrals
        neutral_edges = edge_table.loc[edge_table["Adnotation"].dropna().index]
        kept_edges = [i for i in neutral_edges.index if neutral_edges.loc[i, "node_2"] in ion_idx_mix]

    
        # Get ion edges
        ion_edges = edge_table[edge_table['status'] != "neg_add_edge"]
        ion_edges = ion_edges[ion_edges['status'] != "pos_add_edge"]
        for i in ion_edges.index:
            if ion_edges.loc[i, "node_1"] in ion_idx_mix:
                if ion_edges.loc[i, "node_2"] in ion_idx_mix:
                    kept_edges.append(i)
        kept_edges.sort()
        sample_edges = edge_table.loc[kept_edges]
        sample_edges.sort_values('node_1', inplace = True)
        sample_edges.reset_index(inplace = True, drop = True)
        
        kept_nodes = list(set(list(sample_edges['node_1']) + list(sample_edges['node_2'])))
        kept_nodes.sort()
        sample_nodes = node_table.loc[kept_nodes].copy()
        sample_nodes.drop(pd.Series(samples), axis = 1, inplace = True)
        sample_nodes[sample] = node_table[sample]
        
        sample_nodes.to_csv(out_path + sample + "_node_table.csv", index_label = "Index")
        sample_edges.to_csv(out_path + sample + "_edge_table.csv", index_label = "Index")
    return
