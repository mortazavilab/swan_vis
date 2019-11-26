import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import sqlite3

# creates the duplicate index
def create_dupe_index(df, ind_name):
	df[ind_name+'_back'] = df[ind_name]
	return df

# renames old index dupe column in df and resets the index
def reset_dupe_index(df, ind_name):
	df.rename({ind_name: ind_name+'_back'}, inplace=True, axis=1)
	df.reset_index(inplace=True)
	return(df)

# set index, rename dupe index in df
def set_dupe_index(df, ind_name):
	df.set_index(ind_name, inplace=True)
	df.rename({ind_name+'_back': ind_name}, inplace=True, axis=1)
	return(df)

# partner function to label_edges
def set_edge_attrs(x, G, f_df, f_e):
	attr = {(x.v1, x.v2): {f_e: x[f_df]}}
	nx.set_edge_attributes(G, attr)
	return G

# label edges in G based on fields of edge_df
def label_edges(G, edge_df, f_df, f_e):
	edge_df.apply(lambda x: set_edge_attrs(x, G, f_df, f_e), axis=1)
	return G

# parter function to label_nodes
def set_node_attrs(x, G, f_df, f_n):
	attr = {x.vertex_id: {f_n: x[f_df]}}
	nx.set_node_attributes(G, attr)
	return G

# label nodes in G based on fields of loc_df
def label_nodes(G, loc_df, f_df, f_n):
	loc_df.apply(lambda x: set_node_attrs(x, G, f_df, f_n), axis=1)
	return G

# get value associated with keyword in the 9th column of gtf
def get_field_value(key, fields):
    if key not in fields:
        return None
    else:
        return fields.split(key+' "')[1].split()[0].replace('";','')
        






