# This file contains APIs for loading and preprocessing data files
import os
import pandas as pd

def load_preprocess(dataset_filename, as_df=False, drop_sample_id=True):
    '''
    Loads data csv using pandas and returns X and y for training
    Inputs:
    dataset_filename = filename of dataset to be used. Options are as follows:
    1. 'aak81Dataset_cpmdat.csv'  - 81 out of top 100 genes from feature selection using stability
    2. 'dge_05.csv' - differential gene expression, genes with p < 0.05
    3. 'dge_all.csv' - differential gene expression, all genes
    4. 'overlap_81_dge.csv' - overlap of dge with 81 genes from feature selection
    5. 'cncvae_latent_embedding.csv' - cncVAE latent embedding features z0-z15 (Default option)
    There is no check at the moment for incorrect option!!!
    Filenames are resolved relative to the 'data' folder unless an absolute path is given.
    DEPRECATED: n_genes = # top genes considered for feature stability when doing feature selection

    Outputs:
    X = feature vectors
    y = labels
    '''
    if os.path.isabs(dataset_filename):
        data_fr = pd.read_csv(dataset_filename)
    else:
        data_fr = pd.read_csv(os.path.join('data', dataset_filename))

    # Dropping weird first index column, if present
    if 'Unnamed: 0' in data_fr.columns:
        data_fr.drop(['Unnamed: 0'], axis=1, inplace=True)

    # Make a new column with class = 0 for MGS1 and class = 1 for MGS4
    data_fr['class'] = 0
    data_fr.loc[data_fr['mgs_level'] == 'MGS4', 'class'] = 1
    
    if as_df:
        if drop_sample_id:
            data_fr.drop(['sample_id'], axis=1, inplace=True)
        data_fr.drop(['mgs_level'], axis=1, inplace=True)
        return data_fr
    
    # Split data into features and classes
    y = data_fr['class'].to_numpy()
    data_fr.drop(['class', 'mgs_level', 'sample_id'], axis=1, inplace=True)
    X = data_fr.to_numpy()
    
    return X, y

def separate_df_two_groups(dataset_filename):
    '''
    Loads data and returns two dataframes for the two groups of patients based on the MGS column

    Inputs:
    dataset_filename = filename of dataset to be used, resolved relative to the 'data' folder
    unless an absolute path is given

    Outputs:
    control_df = Data frame of controls
    disease_df = Data frame of disease
    '''

    if os.path.isabs(dataset_filename):
        data_fr = pd.read_csv(dataset_filename)
    else:
        data_fr = pd.read_csv(os.path.join('data', dataset_filename))

    # Dropping weird first index column, if present
    if 'Unnamed: 0' in data_fr.columns:
        data_fr.drop(['Unnamed: 0'], axis=1, inplace=True)

    control_df = data_fr[data_fr['mgs_level'] == 'MGS1'].copy(deep=True)
    control_df.drop(['sample_id', 'mgs_level'], axis=1, inplace=True)

    disease_df = data_fr[data_fr['mgs_level'] == 'MGS4'].copy(deep=True)
    disease_df.drop(['sample_id', 'mgs_level'], axis=1, inplace=True)

    return control_df, disease_df

def read_age_data(n_genes):
    metadata_fr = pd.read_csv(os.path.join('data', 'meta_age.csv'))
    
    # Dropping weird first index column
    metadata_fr.drop('Unnamed: 0', 1, inplace=True)
    
    # Drop all other columns not being used
    metadata_fr.drop(['r_id', 'donor', 'race', 'mgs_level'], 1, inplace=True)
    
    return metadata_fr