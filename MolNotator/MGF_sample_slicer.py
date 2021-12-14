def mgf_slicer(params : dict, ion_mode : str):
    """Split MGF output from MZmine by samples.
    No spectrum from the MGF file must be empty
    The CSV file from MZmine is expected to have only one columns per sample, 
    i.e. "mzXML Peak area" or "mzXML Peak height".
    When treating positive and negative mode data, it is important to
    differenciate both modes for each sample ("POS" and "NEG"), and to have
    otherwise the same exact same for each sample, i.e. "POS_sample_1" and 
    "NEG_sample_1".
    """
    import os
    from tqdm import tqdm
    import pandas as pd
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf

    # Load parameters
    mzmine_suffix= params['mzmine_suffix']
    if ion_mode == "NEG":
        csv_file= params['neg_csv']
        mgf_file= params['neg_mgf']
        in_path= params['neg_out_0']
        out_path= params['neg_out_1']
    elif ion_mode == "POS" :
        csv_file= params['pos_csv']
        mgf_file= params['pos_mgf']
        in_path= params['pos_out_0']
        out_path= params['pos_out_1']
    else:
        print('Ion mode must be either "NEG" or "POS"')
        return

    # Create out folder
    if not os.path.exists(f'{out_path}') :
        os.mkdir(f'{out_path}')
    
    # read the csv metadata file from mzmine
    csv_file = pd.read_csv(f'{in_path}{csv_file}')
    
    # Get the sample list
    samples = pd.Series(csv_file.columns)
    samples = list(samples.str.replace(mzmine_suffix, '.mgf'))
    csv_file.columns = samples
    samples.remove('row ID')
    samples.remove('row m/z')
    samples.remove('row retention time')
    if 'Unnamed: ' in samples[-1] :
        samples.remove(samples[-1])
    
    # MZmine mgf file
    mgf_file = list(load_from_mgf(f'{in_path}/{mgf_file}'))
    
    # write the mgf files for each individual sample:
    print('Slicing the MGF file into sample MGFs:')
    for sample in tqdm(samples) :
        #sample = samples[0]
        new_mgf = list()
        for i in csv_file.index :
            if csv_file.loc[i, sample] > 0 :
                new_mgf.append(mgf_file[i])
        save_as_mgf(new_mgf, f'{out_path}{sample}')
    
