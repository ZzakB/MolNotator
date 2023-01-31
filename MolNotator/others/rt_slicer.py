"""rt_slicer.py - rt_slicer function used by MolNotator"""
def rt_slicer(rt : float, rt_error : float, ion_id, input_table, rt_field : str) :
    """
    Returns a slice of a node table, given an RT, an RT error and a selected
    ion around which coeluted ions are to be search

    Parameters
    ----------
    rt : float
        Retention type of the selected ion.
    rt_error : float
        Retention time window to be used around RT.
    ion_id
        ID of the selected ion.
    input_table : TYPE
        Node table containing the retention time values.
    rt_field : str
        field containing the RT values.

    Returns
    -------
    slice_table : pandas.DataFrame
        Dataframe of ions coeluted with the selected ion.

    """
    rt_low = float(rt) - rt_error
    rt_high = float(rt) + rt_error
    sliced_table = input_table[input_table[rt_field].astype(float).between(rt_low,
                              rt_high, inclusive = "both")].copy()
    return sliced_table.drop(ion_id)
