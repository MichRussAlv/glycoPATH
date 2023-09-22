# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:40:15 2023

@author: russe
"""
import pandas as pd
import numpy as np

# Read csv into data frame
df = pd.read_csv('nglycans.csv', index_col=0)
data = df.to_numpy()

# reshape the array into subarrays with 3 columns
subarrays = data.reshape(data.shape[0], -1, 3)

# take the mean of each subarray
averages = np.mean(subarrays, axis=2)
ave = pd.DataFrame(averages)

# add index column using Compounds list
cpd_column = df.index.to_series()
cpd_column = cpd_column.reset_index(drop=True)
cpds = pd.DataFrame(cpd_column)

# Read group labels
#groups = pd.read_csv('groups.csv')
#groups = groups.str.replace(r'\(|\)','')
#groups.columns = groups.columns.str.replace(r'\(|\)', '')
groups = pd.read_csv('groups.csv')
df_groups = groups.columns.tolist()
#print(df_groups)

# concatenate compounds list and averages, then export
ave_df = pd.concat([cpds, ave], axis=1)
ave_df.columns = df_groups
ave_df.to_csv('n-glycan_averages.csv', index=False)

