#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 01:40:02 2023

@author: michaelrussellealvarez
"""

import csv
import graphviz
import pandas as pd


# Step 1: Import correlations from previous step and filter genes based on known enzyme reactions & highest correlations
###############################################################################

# Import data from HM_correlations and perform filtering
filename_HM = 'HM_correlations_mannosidases.csv'
df_HM = pd.read_csv(filename_HM, index_col=0)
lowest_values_HM = df_HM.min().round(3)
corresponding_headers_HM = df_HM.idxmin()

HM_starts = []  
HM_ends = []

with open('HM_reactions.csv', 'r') as file:
    reader_HM = csv.reader(file)
    next(reader_HM)  # Skip header row
    for row in reader_HM:
        HM_starts.append(row[0])  # Append values to the lists
        HM_ends.append(row[1])

result_HM = pd.DataFrame({'start_node': HM_starts, 'end_node': HM_ends, 'Gene': corresponding_headers_HM, 'Correlation': lowest_values_HM})
result_HM.to_csv('HM_edges.csv', index=True)


# Import data from U_correlations_GlcNActransferases and perform filtering
filename_Uglc = 'U_correlations_glcnactransferases.csv'
df_Uglc = pd.read_csv(filename_Uglc, index_col=0)
highest_values_Uglc = df_Uglc.max().round(3)
corresponding_headers_Uglc = df_Uglc.idxmax()

Uglc_starts = []  
Uglc_ends = []

with open('U_reactions_glcnactransferases.csv', 'r') as file:
    reader_Uglc = csv.reader(file)
    next(reader_Uglc)  # Skip header row
    for row in reader_Uglc:
        Uglc_starts.append(row[0])  # Append values to the lists
        Uglc_ends.append(row[1])

result_Uglc = pd.DataFrame({'start_node': Uglc_starts, 'end_node': Uglc_ends, 'Gene': corresponding_headers_Uglc, 'Correlation': highest_values_Uglc})
result_Uglc.to_csv('U_edges_glcnactransferases.csv', index=True)

# Import data from U_correlations_Galactosyltransferases and perform filtering
filename_Ugal = 'U_correlations_galactosyltransferases.csv'
df_Ugal = pd.read_csv(filename_Ugal, index_col=0)
highest_values_Ugal = df_Ugal.max().round(3)
corresponding_headers_Ugal = df_Ugal.idxmax()

Ugal_starts = []  
Ugal_ends = []

with open('U_reactions_galactosyltransferases.csv', 'r') as file:
    reader_Ugal = csv.reader(file)
    next(reader_Ugal)  # Skip header row
    for row in reader_Ugal:
        Ugal_starts.append(row[0])  # Append values to the lists
        Ugal_ends.append(row[1])

result_Ugal = pd.DataFrame({'start_node': Ugal_starts, 'end_node': Ugal_ends, 'Gene': corresponding_headers_Ugal, 'Correlation': highest_values_Ugal})
result_Ugal.to_csv('U_edges_galactosyltransferases.csv', index=True)


# Import data from F_correlations and perform filtering
filename_F = 'F_correlations_fucosyltransferases.csv'
df_F = pd.read_csv(filename_F, index_col=0)
highest_values_F = df_F.max().round(3)
corresponding_headers_F = df_F.idxmax()

F_starts = []  
F_ends = []

with open('F_reactions.csv', 'r') as file:
    reader_F = csv.reader(file)
    next(reader_F)  # Skip header row
    for row in reader_F:
        F_starts.append(row[0])  # Append values to the lists
        F_ends.append(row[1])

result_F = pd.DataFrame({'start_node': F_starts, 'end_node': F_ends, 'Gene': corresponding_headers_F, 'Correlation': highest_values_F})
result_F.to_csv('F_edges.csv', index=True)


# Import data from S_correlations and perform filtering
filename_S = 'S_correlations_sialyltransferases.csv'
df_S = pd.read_csv(filename_S, index_col=0)
highest_values_S = df_S.max().round(3)
corresponding_headers_S = df_S.idxmax()

S_starts = []  
S_ends = []

with open('S_reactions.csv', 'r') as file:
    reader_S = csv.reader(file)
    next(reader_S)  # Skip header row
    for row in reader_S:
        S_starts.append(row[0])  # Append values to the lists
        S_ends.append(row[1])

result_S = pd.DataFrame({'start_node': S_starts, 'end_node': S_ends, 'Gene': corresponding_headers_S, 'Correlation': highest_values_S})
result_S.to_csv('S_edges.csv', index=True)


# Import data from SF_correlations_fucosyltransferases and perform filtering
filename_SFfut = 'SF_correlations_fucosyltransferases.csv'
df_SFfut = pd.read_csv(filename_SFfut, index_col=0)
highest_values_SFfut = df_SFfut.max().round(3)
corresponding_headers_SFfut = df_SFfut.idxmax()

SFfut_starts = []  
SFfut_ends = []

with open('SF_reactions_fucosyltransferases.csv', 'r') as file:
    reader_SFfut = csv.reader(file)
    next(reader_SFfut)  # Skip header row
    for row in reader_SFfut:
        SFfut_starts.append(row[0])  # Append values to the lists
        SFfut_ends.append(row[1])

result_SFfut = pd.DataFrame({'start_node': SFfut_starts, 'end_node': SFfut_ends, 'Gene': corresponding_headers_SFfut, 'Correlation': highest_values_SFfut})
result_SFfut.to_csv('SF_edges_fucosyltransferases.csv', index=True)

# Import data from SF_correlations_sialyltransferases and perform filtering
filename_SFsia = 'SF_correlations_sialyltransferases.csv'
df_SFsia = pd.read_csv(filename_SFsia, index_col=0)
highest_values_SFsia = df_SFsia.max().round(3)
corresponding_headers_SFsia = df_SFsia.idxmax()

SFsia_starts = []  
SFsia_ends = []

with open('SF_reactions_sialyltransferases.csv', 'r') as file:
    reader_SFsia = csv.reader(file)
    next(reader_SFsia)  # Skip header row
    for row in reader_SFsia:
        SFsia_starts.append(row[0])  # Append values to the lists
        SFsia_ends.append(row[1])

result_SFsia = pd.DataFrame({'start_node': SFsia_starts, 'end_node': SFsia_ends, 'Gene': corresponding_headers_SFsia, 'Correlation': highest_values_SFsia})
result_SFsia.to_csv('SF_edges_sialyltransferases.csv', index=True)


###############################################################################


# Step 2: Create network visualizing the correlating glyco-enzymes
###############################################################################

# Create network for HM correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('HM_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('HM_edges.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('HM', view=True, cleanup=True)


# Create network for Uglc correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('U_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('U_edges_glcnactransferases.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('U_glc', view=True, cleanup=True)


# Create network for Ugal correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('U_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('U_edges_galactosyltransferases.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('U_gal', view=True, cleanup=True)

# Create network for F correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('F_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('F_edges.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('F', view=True, cleanup=True)


# Create network for S correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('S_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('S_edges.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('S', view=True, cleanup=True)

# Create network for SFfut correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('SFfut_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('SF_edges_fucosyltransferases.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('SFfut', view=True, cleanup=True)

# Create network for SFsia correlations
graph = graphviz.Digraph(format='png', graph_attr={'rankdir': 'LR'})

with open('SFsia_nodes.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        node_id = row[0]
        image_path = row[1]
        label_text = row[2]
        graph.node(node_id, shape='none', image=image_path, label=label_text, labelloc='t', labeljust='c', fontname='Arial', fontsize='15')

with open('SF_edges_sialyltransferases.csv', 'r') as file:
    reader = csv.reader(file)

    next(reader)  
    for row in reader:
        start_node = row[1]
        end_node = row[2]
        edge_label_above = row[3]
        edge_label_below = row[4]
        
        edge_label = f'{edge_label_above}\\nR= {edge_label_below}'
        graph.edge(start_node, end_node, label=edge_label, fontname='Arial', fontsize='20')

graph.render('SFsia', view=True, cleanup=True)
