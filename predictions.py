# -*- coding: utf-8 -*-
"""
Created on Sun May 21 15:33:00 2023

@author: russe
"""

import pandas as pd
import os
import numpy as np
from scipy.interpolate import splrep, splev
import pickle
import glob
import matplotlib.pyplot as plt


# Step 1: Import averages from previous steps & filter genes and glycans 
###############################################################################

# Import data from transcriptomics (df1) and glycomics (df2)
df_transcriptomics = pd.read_csv('gene_averages.csv')
df_glycomics = pd.read_csv('n-glycan_averages.csv')

# Read CSV files of glycans and enzymes
column_HM_filter_df = pd.read_csv('HM_glycan_filter.csv')
row_HM_filter_df = pd.read_csv('HM_gene_filter.csv')

column_U_filter_df = pd.read_csv('U_glycan_filter.csv')
row_U_filter_df = pd.read_csv('U_gene_filter.csv')

column_F_filter_df = pd.read_csv('F_glycan_filter.csv')
row_F_filter_df = pd.read_csv('F_gene_filter.csv')

column_S_filter_df = pd.read_csv('S_glycan_filter.csv')
row_S_filter_df = pd.read_csv('S_gene_filter.csv')

column_SF_filter_df = pd.read_csv('SF_glycan_filter.csv')
row_SF_filter_df = pd.read_csv('SF_gene_filter.csv')

# Extract filter strings from the CSV files
column_HM_filter = column_HM_filter_df['Compounds'].tolist()
row_HM_filter = row_HM_filter_df['Genes'].tolist()

column_U_filter = column_U_filter_df['Compounds'].tolist()
row_U_filter = row_U_filter_df['Genes'].tolist()

column_F_filter = column_F_filter_df['Compounds'].tolist()
row_F_filter = row_F_filter_df['Genes'].tolist()

column_S_filter = column_S_filter_df['Compounds'].tolist()
row_S_filter = row_S_filter_df['Genes'].tolist()

column_SF_filter = column_SF_filter_df['Compounds'].tolist()
row_SF_filter = row_SF_filter_df['Genes'].tolist()

# Filter transcriptomics based on genes
filtered_HM_genes = df_transcriptomics[df_transcriptomics['Genes'].isin(row_HM_filter)]
filtered_U_genes = df_transcriptomics[df_transcriptomics['Genes'].isin(row_U_filter)]
filtered_F_genes = df_transcriptomics[df_transcriptomics['Genes'].isin(row_F_filter)]
filtered_S_genes = df_transcriptomics[df_transcriptomics['Genes'].isin(row_S_filter)]
filtered_SF_genes = df_transcriptomics[df_transcriptomics['Genes'].isin(row_SF_filter)]

filtered_HM_genes.to_csv('HM_genes.csv', index=False)
filtered_U_genes.to_csv('U_genes.csv', index=False)
filtered_F_genes.to_csv('F_genes.csv', index=False)
filtered_S_genes.to_csv('S_genes.csv', index=False)
filtered_SF_genes.to_csv('SF_genes.csv', index=False)

# Filter glycomics based on compounds
filtered_HM_glycans = df_glycomics[df_glycomics['Compounds'].isin(column_HM_filter)]
filtered_U_glycans = df_glycomics[df_glycomics['Compounds'].isin(column_U_filter)]
filtered_F_glycans = df_glycomics[df_glycomics['Compounds'].isin(column_F_filter)]
filtered_S_glycans = df_glycomics[df_glycomics['Compounds'].isin(column_S_filter)]
filtered_SF_glycans = df_glycomics[df_glycomics['Compounds'].isin(column_SF_filter)]

filtered_HM_glycans.to_csv('HM_glycans.csv', index=False)
filtered_U_glycans.to_csv('U_glycans.csv', index=False)
filtered_F_glycans.to_csv('F_glycans.csv', index=False)
filtered_S_glycans.to_csv('S_glycans.csv', index=False)
filtered_SF_glycans.to_csv('SF_glycans.csv', index=False)

###############################################################################


# Step 2: Pair-up genes and glycans per glycan types
###############################################################################

# Import filtered data
filtered_HM_genes_X = pd.read_csv('HM_genes.csv', index_col=0)
filtered_HM_glycans_Y = pd.read_csv('HM_glycans.csv', index_col=0)
transposed_HM_genes_X = filtered_HM_genes_X.T
transposed_HM_glycans_Y = filtered_HM_glycans_Y.T

filtered_U_genes_X = pd.read_csv('U_genes.csv', index_col=0)
filtered_U_glycans_Y = pd.read_csv('U_glycans.csv', index_col=0)
transposed_U_genes_X = filtered_U_genes_X.T
transposed_U_glycans_Y = filtered_U_glycans_Y.T

filtered_F_genes_X = pd.read_csv('F_genes.csv', index_col=0)
filtered_F_glycans_Y = pd.read_csv('F_glycans.csv', index_col=0)
transposed_F_genes_X = filtered_F_genes_X.T
transposed_F_glycans_Y = filtered_F_glycans_Y.T

filtered_S_genes_X = pd.read_csv('S_genes.csv', index_col=0)
filtered_S_glycans_Y = pd.read_csv('S_glycans.csv', index_col=0)
transposed_S_genes_X = filtered_S_genes_X.T
transposed_S_glycans_Y = filtered_S_glycans_Y.T

filtered_SF_genes_X = pd.read_csv('SF_genes.csv', index_col=0)
filtered_SF_glycans_Y = pd.read_csv('SF_glycans.csv', index_col=0)
transposed_SF_genes_X = filtered_SF_genes_X.T
transposed_SF_glycans_Y = filtered_SF_glycans_Y.T

# Create x,y pairs for each Gene-Glycan
paired_HM = pd.DataFrame()
paired_U = pd.DataFrame()
paired_F = pd.DataFrame()
paired_S = pd.DataFrame()
paired_SF = pd.DataFrame()

for col1 in transposed_HM_genes_X.columns:
    for col2 in transposed_HM_glycans_Y.columns:
        paired_HM = pd.DataFrame({f'{col1}': transposed_HM_genes_X[col1], f'{col2}': transposed_HM_glycans_Y[col2]})
        filename_HM = f'{col1}_{col2}.csv'
        export_folder_HM = 'gene-glycan/HM'
        filepath = os.path.join(export_folder_HM, filename_HM)
        paired_HM.to_csv(filepath, index=False)

for col1 in transposed_U_genes_X.columns:
    for col2 in transposed_U_glycans_Y.columns:
        paired_U = pd.DataFrame({f'{col1}': transposed_U_genes_X[col1], f'{col2}': transposed_U_glycans_Y[col2]})
        filename_U = f'{col1}_{col2}.csv'
        export_folder_U = 'gene-glycan/U'
        filepath = os.path.join(export_folder_U, filename_U)
        paired_U.to_csv(filepath, index=False)

for col1 in transposed_F_genes_X.columns:
    for col2 in transposed_F_glycans_Y.columns:
        paired_F = pd.DataFrame({f'{col1}': transposed_F_genes_X[col1], f'{col2}': transposed_F_glycans_Y[col2]})
        filename_F = f'{col1}_{col2}.csv'
        export_folder_F = 'gene-glycan/F'
        filepath = os.path.join(export_folder_F, filename_F)
        paired_F.to_csv(filepath, index=False)

for col1 in transposed_S_genes_X.columns:
    for col2 in transposed_S_glycans_Y.columns:
        paired_S = pd.DataFrame({f'{col1}': transposed_S_genes_X[col1], f'{col2}': transposed_S_glycans_Y[col2]})
        filename_S = f'{col1}_{col2}.csv'
        export_folder_S = 'gene-glycan/S'
        filepath = os.path.join(export_folder_S, filename_S)
        paired_S.to_csv(filepath, index=False)
        
for col1 in transposed_SF_genes_X.columns:
    for col2 in transposed_SF_glycans_Y.columns:
        paired_SF = pd.DataFrame({f'{col1}': transposed_SF_genes_X[col1], f'{col2}': transposed_SF_glycans_Y[col2]})
        filename_SF = f'{col1}_{col2}.csv'
        export_folder_SF = 'gene-glycan/SF'
        filepath = os.path.join(export_folder_SF, filename_SF)
        paired_SF.to_csv(filepath, index=False)
        
###############################################################################


# Step 3: Fit the spline models to each gene-glycan pair and predict unknown values from file
###############################################################################

import pandas as pd
import os
import numpy as np
from scipy.interpolate import splrep, splev
import pickle
import glob
import matplotlib.pyplot as plt

unknown_transcriptomics = pd.read_csv('ccd_genes.csv')

# Filter transcriptomics based on genes
unknown_HM_genes = unknown_transcriptomics[unknown_transcriptomics['Genes'].isin(row_HM_filter)].T
unknown_U_genes = unknown_transcriptomics[unknown_transcriptomics['Genes'].isin(row_U_filter)].T
unknown_F_genes = unknown_transcriptomics[unknown_transcriptomics['Genes'].isin(row_F_filter)].T
unknown_S_genes = unknown_transcriptomics[unknown_transcriptomics['Genes'].isin(row_S_filter)].T
unknown_SF_genes = unknown_transcriptomics[unknown_transcriptomics['Genes'].isin(row_SF_filter)].T

unknown_HM_genes.to_csv('unknown_HM_genes.csv', header=False)
unknown_U_genes.to_csv('unknown_U_genes.csv', header=False)
unknown_F_genes.to_csv('unknown_F_genes.csv', header=False)
unknown_S_genes.to_csv('unknown_S_genes.csv', header=False)
unknown_SF_genes.to_csv('unknown_SF_genes.csv', header=False)

# Load the unknown data from the CSV file
unknown_data = pd.read_csv('unknown_HM_genes.csv', usecols=lambda x:x!=0)
unknown_column = 'MA2A2_HUMAN' # Specify gene to use prediction model on

# Define the folder containing the gene-glycan pairs
GG_folder = 'gene-glycan/HM' # Specify N-glycan type to make predictions on

# Get a list of CSV files in the specified folder
csv_files = glob.glob(os.path.join(GG_folder, '*.csv'))

# Create a list to store spline models
spline_models = []

# Process each CSV file
for file_path in csv_files:
    file_name = os.path.basename(file_path)
    GGpair = pd.read_csv(file_path)

    x_col = GGpair.columns[0]
    y_col = GGpair.columns[1]

    x = GGpair[x_col].values
    y = GGpair[y_col].values

    sorted_indices = np.argsort(x)
    x = x[sorted_indices]
    y = y[sorted_indices]

    spl = splrep(x, y, k=2, task=0, s=32) # Optimize combinations of k (degree) and s (smoothing factor)
    spline_models.append(spl)

# Define the folder to store the spline models
spline_folder = 'gene-glycan/spline'
predictions_folder = f'predictions/{unknown_column}'
os.makedirs(spline_folder, exist_ok=True)
os.makedirs(predictions_folder, exist_ok=True)

# Save each spline model as a separate pickle file and plot the spline curve
for i, spline in enumerate(spline_models):
    file_path = csv_files[i]
    file_name = os.path.basename(file_path)
    spline_file_path = os.path.join(spline_folder, f"spline_{file_name}.pickle")
    image_file_path = os.path.join(spline_folder, f"spline_{file_name}.png")

    with open(spline_file_path, 'wb') as f:
        pickle.dump(spline, f)

    results_df = pd.DataFrame()
    
    x = spline[0]
    y = spline[1]

    x_eval = np.linspace(np.min(x), np.max(x))
    y_eval = splev(x_eval, spline)
    plt.plot(x_eval, y_eval, label=file_name)

    plt.xlabel('Glycogene expression')
    plt.ylabel('N-glycan abundance')
    plt.legend()

    plt.savefig(image_file_path)
    
    x_unknown = unknown_data[unknown_column].values.astype(float)
    y_predicted = splev(x_unknown, spline)
    
    result = pd.DataFrame({f'{unknown_column}': x_unknown, 'y_predicted': y_predicted})
    results_df[f'Spline Model {file_name}'] = result['y_predicted']
    
    # Print the predicted values
    print("Predicting values for spline model", i+1)

    results_df.to_csv(f'predictions/{unknown_column}/{unknown_column}_Spline Model {file_name}.csv', index=False)
    

###############################################################################


# Step 4: Clean-up (optional to remove prediction of unknown gene expression using other models)
###############################################################################
csv_files = glob.glob(os.path.join(predictions_folder, '*.csv'))
for file_path in csv_files:
    file_name = os.path.basename(file_path)
    count = file_name.count(unknown_column)
    if count !=2:
        os.remove(file_path)
        
   