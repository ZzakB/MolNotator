"""singleton_edges.py - singleton_edges function for the fragnotator module"""
import pandas as pd

def singleton_edges(node_table, edge_table):
    """
    Adds singleton edges to an edge table already containing precursor-fragment
    ion pairs.

    Parameters
    ----------
    node_table : pandas.DataFrame
        Dataframe containing all ions, provided by the fragnotator module
    edge_table : pandas.DataFrame
        Dataframe containing all already paired ions.

    Returns
    -------
    edge_table : pandas.DataFrame
        Dataframe/Edge table update with the singleton edges.

    """

    paired_nodes = set(list(edge_table['node_1']) + list(edge_table['node_2']))
    singleton_nodes = list(set(node_table.index) - paired_nodes)
    singleton_edge_table = pd.DataFrame(0.0, index = range(len(singleton_nodes)), 
                                        columns = edge_table.columns)
    singleton_edge_table['node_1'] = singleton_nodes
    singleton_edge_table['node_2'] = singleton_nodes
    edge_table['status'] = ["frag_edge"]*len(edge_table)
    singleton_edge_table['status'] = ["singleton"]*len(singleton_edge_table)
    edge_table = edge_table.append(singleton_edge_table, ignore_index = False, sort = False)
    edge_table.reset_index(drop = True, inplace = True)
    return edge_table
