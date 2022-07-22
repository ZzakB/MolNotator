"""reindexer.py - reindexer function for MolNotator"""
def reindexer(node_table, params):
    """
    Reindexes tables using the supplied "index_col" argument in the params dict
    and adds a "spec_id" column for the relative position of ions within the
    spectrum file.

    Parameters
    ----------
    node_table : pandas.DataFrame
        Dataframe containing the metadata for each ions.
    params : dict
        Dictionary containing the global parameters for the process.

    Returns
    -------
    node_table : pandas.DataFrame
        Same dataframe but reindexed.

    """
    node_table[params['index_col']] = node_table[params['index_col']].astype(int)
    node_table.set_index(params['index_col'], inplace = True)
    node_table.insert(0, 'spec_id', range(len(node_table)))
    return node_table
