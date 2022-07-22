"""spectrum_extractor.py - Spectrum extractor module for MolNotator"""
from tqdm import tqdm
import pandas as pd

def spectrum_extractor(spectrum_list : list, add_tic : bool = False) :
    
    """
    Takes a list of spectra (matchms) to extract the metadata and outputs it as
    a pandas dataframe.
    Parameters
    ----------
    spectrum_list : list
        List of matchms.Spectrum objects.
    add_tic : bool, optional
        Whether or not to add a TIC (Total Ion Current) field for each ion /
        spectrum. The default is False.
    Returns
    -------
    node_table : pandas.DataFrame
        Dataframe containing all the metadata from the list of Spectrum objects.
    """
    
    # Get all columns from the spectrum list
    cols = list(spectrum_list[0].metadata.keys())
    
    # Initialise node table
    node_table = pd.DataFrame(index = range(len(spectrum_list)), columns = cols)
    node_table = node_table.fillna(0.0)
    
    # If add_tic is True
    if add_tic: 
        for i in tqdm(range(len(spectrum_list))) :
            for col in cols :
                if col == "pepmass" : # pepmass is a tuple with matchms, exception.
                    node_table.loc[i, col] = spectrum_list[i].get(col)[0]
                else:
                    node_table.loc[i, col] = spectrum_list[i].get(col)
            node_table.loc[i, "TIC"] = sum(spectrum_list[i].peaks.intensities)
    
    # If no TIC column required
    else:
        for i in tqdm(range(len(spectrum_list))) :
            for col in cols :
                if col == "pepmass" : # pepmass is a tuple with matchms, exception.
                    node_table.loc[i, col] = spectrum_list[i].get(col)[0]
                else:
                    node_table.loc[i, col] = spectrum_list[i].get(col)
    return node_table
