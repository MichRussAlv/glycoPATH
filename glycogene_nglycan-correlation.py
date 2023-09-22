# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:46:44 2023

@author: russe
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Import averages from previous steps & calculate correlation values
###############################################################################

# Import data from transcriptomics (df1) and glycomics (df2)
df1 = pd.read_csv('gene_averages.csv', index_col=0)
df2 = pd.read_csv('n-glycan_averages.csv', index_col=0)

# Calculate the row-wise correlation between the two DataFrames
#df2_transposed = df2.T
#correlation_matrix = df1.apply(lambda row: df2_transposed.corrwith(row), axis=1)

# Export full correlation matrix
#correlation_matrix.to_csv('correlations.csv')
correlation_matrix = pd.read_csv('correlations.csv', index_col=0)

###############################################################################


# Step 2: Filter correlation matrix by N-glycan categories
###############################################################################

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

# Filter DataFrame based on column and row values
filtered_HM = correlation_matrix[column_HM_filter].loc[row_HM_filter]
filtered_U = correlation_matrix[column_U_filter].loc[row_U_filter]
filtered_F = correlation_matrix[column_F_filter].loc[row_F_filter]
filtered_S = correlation_matrix[column_S_filter].loc[row_S_filter]
filtered_SF = correlation_matrix[column_SF_filter].loc[row_SF_filter]

# Export heatmap data
filtered_HM.to_csv('HM_correlations.csv')
filtered_U.to_csv('U_correlations.csv')
filtered_F.to_csv('F_correlations.csv')
filtered_S.to_csv('S_correlations.csv')
filtered_SF.to_csv('SF_correlations.csv')

###############################################################################


# Step 3: Plot filtered & categorized heatmaps with corresponding figures
###############################################################################

cell_height = 0.5

# Locate N-glycan images and convert them to numpy arrays
HM_files = ['9200.png','8200.png','7200.png','6200.png','5200.png','4200.png']
HM_images = [plt.imread(file) for file in HM_files]
U_files = ['4400.png','5400.png','6500.png','7600.png']
U_images = [plt.imread(file) for file in U_files]
F_files = ['5410.png','5420.png','5430.png','6510.png','6520.png','6530.png','6540.png','7610.png','7620.png','7630.png','7640.png','7650.png']
F_images = [plt.imread(file) for file in F_files]
S_files = ['5401.png','5402.png','6501.png','6502.png','6503.png','7601.png','7602.png','7603.png']
S_images = [plt.imread(file) for file in S_files]
SF_files = ['5411.png','5412.png','5421.png','5422.png','5431.png','5432.png','6511.png','6512.png','6513.png','6521.png','6522.png','6523.png','6531.png','6532.png','6541.png','6542.png','6543.png','7611.png','7612.png','7613.png','7614.png','7621.png','7622.png','7623.png','7631.png','7632.png','7633.png','7634.png','7641.png','7642.png','7643.png','7644.png','7651.png','7652.png','7653.png']
SF_images = [plt.imread(file) for file in SF_files]

# Create the figures and axes
HM_figure_height = cell_height * filtered_HM.shape[0]
HM_figure_width = filtered_HM.shape[1]
HM_fig, HM_ax = plt.subplots(figsize=(HM_figure_width, HM_figure_height))

U_figure_height = cell_height * filtered_U.shape[0]
U_figure_width = filtered_U.shape[1]
U_fig, U_ax = plt.subplots(figsize=(U_figure_width, U_figure_height))

F_figure_height = cell_height * filtered_F.shape[0]
F_figure_width = filtered_F.shape[1]
F_fig, F_ax = plt.subplots(figsize=(F_figure_width, F_figure_height))

S_figure_height = cell_height * filtered_S.shape[0]
S_figure_width = filtered_S.shape[1]
S_fig, S_ax = plt.subplots(figsize=(S_figure_width, S_figure_height))

SF_figure_height = cell_height * filtered_SF.shape[0]
SF_figure_width = filtered_SF.shape[1]
SF_fig, SF_ax = plt.subplots(figsize=(SF_figure_width, SF_figure_height))

# Plot the heatmaps
HM_heatmap = HM_ax.imshow(filtered_HM, cmap='coolwarm', aspect='auto')
HM_ax.xaxis.tick_top()
HM_ax.set_xticks(np.arange(filtered_HM.shape[1]))
HM_ax.set_xticklabels(column_HM_filter, rotation=45)
for i, image in enumerate(HM_images):
    HM_ax.text(i, -0.5, '', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='none'))
    HM_ax.imshow(image, extent=[i-0.5,i+0.5,-1,0], aspect='auto')
HM_ax.set_xlim(-0.5, len(HM_images)-0.5)
HM_ax.set_yticks(np.arange(filtered_HM.shape[0]))
HM_ax.set_yticklabels(row_HM_filter)
HM_ax.tick_params(axis='both', labelsize=8)
HM_cbar = HM_fig.colorbar(HM_heatmap)

U_heatmap = U_ax.imshow(filtered_U, cmap='coolwarm', aspect='auto')
U_ax.xaxis.tick_top()
U_ax.set_xticks(np.arange(filtered_U.shape[1]))
U_ax.set_xticklabels(column_U_filter, rotation=45)
for i, image in enumerate(U_images):
    U_ax.text(i, -0.5, '', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='none'))
    U_ax.imshow(image, extent=[i-0.5,i+0.5,-1,0], aspect='auto')
U_ax.set_xlim(-0.5, len(U_images)-0.5)
U_ax.set_yticks(np.arange(filtered_U.shape[0]))
U_ax.set_yticklabels(row_U_filter)
U_ax.tick_params(axis='both', labelsize=8)
U_cbar = U_fig.colorbar(U_heatmap)

F_heatmap = F_ax.imshow(filtered_F, cmap='coolwarm', aspect='auto')
F_ax.xaxis.tick_top()
F_ax.set_xticks(np.arange(filtered_F.shape[1]))
F_ax.set_xticklabels(column_F_filter, rotation=45)
for i, image in enumerate(F_images):
    F_ax.text(i, -0.5, '', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='none'))
    F_ax.imshow(image, extent=[i-0.5,i+0.5,-1,0], aspect='auto')
F_ax.set_xlim(-0.5, len(F_images)-0.5)
F_ax.set_yticks(np.arange(filtered_F.shape[0]))
F_ax.set_yticklabels(row_F_filter)
F_ax.tick_params(axis='both', labelsize=8)
F_cbar = F_fig.colorbar(F_heatmap)

S_heatmap = S_ax.imshow(filtered_S, cmap='coolwarm', aspect='auto')
S_ax.xaxis.tick_top()
S_ax.set_xticks(np.arange(filtered_S.shape[1]))
S_ax.set_xticklabels(column_S_filter, rotation=45)
for i, image in enumerate(S_images):
    S_ax.text(i, -0.5, '', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='none'))
    S_ax.imshow(image, extent=[i-0.5,i+0.5,-1,0], aspect='auto')
S_ax.set_xlim(-0.5, len(S_images)-0.5)
S_ax.set_yticks(np.arange(filtered_S.shape[0]))
S_ax.set_yticklabels(row_S_filter)
S_ax.tick_params(axis='both', labelsize=8)
S_cbar = S_fig.colorbar(S_heatmap)

SF_heatmap = SF_ax.imshow(filtered_SF, cmap='coolwarm', aspect='auto')
SF_ax.xaxis.tick_top()
SF_ax.set_xticks(np.arange(filtered_SF.shape[1]))
SF_ax.set_xticklabels(column_SF_filter, rotation=45)
for i, image in enumerate(SF_images):
    SF_ax.text(i, -0.5, '', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='none'))
    SF_ax.imshow(image, extent=[i-0.5,i+0.5,-1,0], aspect='auto')
SF_ax.set_xlim(-0.5, len(SF_images)-0.5)
SF_ax.set_yticks(np.arange(filtered_SF.shape[0]))
SF_ax.set_yticklabels(row_SF_filter)
SF_ax.tick_params(axis='both', labelsize=8)
SF_cbar = SF_fig.colorbar(SF_heatmap)

# Show & export the heatmaps
HM_fig.savefig('HM_heatmap.png', bbox_inches='tight', format='png', dpi=600)
U_fig.savefig('U_heatmap.png', bbox_inches='tight', format='png', dpi=600)
F_fig.savefig('F_heatmap.png', bbox_inches='tight', format='png', dpi=600)
S_fig.savefig('S_heatmap.png', bbox_inches='tight', format='png', dpi=600)
SF_fig.savefig('SF_heatmap.png', bbox_inches='tight', format='png', dpi=600)