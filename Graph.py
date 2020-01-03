import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import copy
from collections import defaultdict
import sqlite3
from utils import *

# super class that both SpliceGraph and PlottedGraph inherit from.
# all functions that both subclasses will use are in here.
class Graph:
	def __init__(self):

		self.datasets = []
		self.counts = []
		self.tpm = []

		self.pg = None
		
		self.loc_df = pd.DataFrame(columns=['chrom', 'coord',
									   'strand','vertex_id',
									   'TSS', 'alt_TSS',
									   'TES', 'alt_TES',
									   'internal'])
		self.edge_df = pd.DataFrame(columns=['edge_id', 'edge_type',
									    'strand', 'v1', 'v2'])
		self.t_df = pd.DataFrame(columns=['tid', 'gid',
									 'gname', 'path'])

	##########################################################################
	####################### Related to creating Graph ########################
	##########################################################################

	# update ids according to coordinates in loc_df, edge_df, and t_df
	def update_ids(self):
		id_map = self.get_ordered_id_map()
		self.update_loc_df_ids(id_map)
		self.update_edge_df_ids(id_map)
		self.update_t_df_paths(id_map)

	# get a dictionary mapping vertex id to ordered new vertex id
	def get_ordered_id_map(self):

		# get the strandedness of entry
		strand = self.loc_df.loc[self.loc_df.index[0], 'strand']

		# sort based on coord depending on strandedness
		if strand == '+':
			self.loc_df.sort_values(by='coord', 
									ascending=True,
									inplace=True)
		elif strand == '-':
			self.loc_df.sort_values(by='coord',
									ascending=False,
									inplace=True)

		# dictionary mapping vertex_id to new_id
		self.loc_df['new_id'] = [i for i in range(len(self.loc_df.index))]

		# account for combined nodes, will still be labelled as c#,
		# but can reorder the '#' according to starting coord
		if 'combined' in self.loc_df.columns:
			combined_ids = self.loc_df.loc[self.loc_df.combined, 'vertex_id']
			combined_new_ids = ['c{}'.format(i) for i in range(len(combined_ids))]
			self.loc_df.loc[combined_ids, 'new_id'] = combined_new_ids

		id_map = self.loc_df['new_id'].to_dict()
		self.loc_df.drop('new_id', axis=1, inplace=True)

		return id_map

	# update vertex ids in loc_df
	def update_loc_df_ids(self, id_map):

		# add new id to df and use it to update vertex id columns
		self.loc_df['new_id'] = self.loc_df.apply(
			lambda x: id_map[x.vertex_id], axis=1)
		self.loc_df.vertex_id = self.loc_df.new_id
		self.loc_df.drop('new_id', axis=1, inplace=True)
		self.loc_df.reset_index(drop=True, inplace=True)
		self.loc_df = create_dupe_index(self.loc_df, 'vertex_id')
		self.loc_df = set_dupe_index(self.loc_df, 'vertex_id')

	# update vertex ids in edge_df
	def update_edge_df_ids(self, id_map):

		# remove old edge_id index and replace edge_id, v1, v2 
		# with values from id_map
		self.edge_df.reset_index(drop=True, inplace=True)
		self.edge_df.v1 = self.edge_df.apply(
			lambda x: id_map[x.v1], axis=1)
		self.edge_df.v2 = self.edge_df.apply(
			lambda x: id_map[x.v2], axis=1)
		self.edge_df.edge_id = self.edge_df.apply(
			lambda x: (x.v1, x.v2), axis=1)

		self.edge_df = create_dupe_index(self.edge_df, 'edge_id')
		self.edge_df = set_dupe_index(self.edge_df, 'edge_id')

	# update vertex ids in t_df
	def update_t_df_paths(self, id_map):
		# update vertex ids in the path
		self.t_df.path = self.t_df.apply(
			lambda x: [id_map[n] for n in x.path], axis=1)

	# create the graph object from the dataframes
	def create_graph_from_dfs(self):

		G = nx.DiGraph()

		# add nodes to graph from transcript paths
		paths = self.t_df.path.tolist()
		for path in paths:
			nx.add_path(G, path)

		self.G = G

	# order edge df based on source id
	def order_edge_df(self):
		self.edge_df.sort_values(by='v1', inplace=True)

	##########################################################################
	############################# Other utilities ############################
	##########################################################################

	# check if anything has been added to the graph yet
	def is_empty(self):
		if len(self.datasets) == 0: 
			return True
		else: 
			return False

	# gets the names of the dataset columns in the graph
	# returns None if no datasets have been added
	def get_dataset_cols(self):
		if len(self.datasets) == 0:
			return None
		return self.datasets

	# gets the names of the counts columns in the graph
	# returns None if no counts have been added
	def get_count_cols(self):
		if len(self.counts) == 0:
			return None
		return self.counts

	# gets the names of tpm columns in the graph
	# returns None if no counts have been added
	def get_tpm_cols(self):
		if len(self.tpm) == 0:
			return None
		return self.tpm

	# gets strandedness of transcript from transcript id
	def get_strand_from_tid(self, tid):
		return self.loc_df.loc[self.t_df.loc[tid, 'path'][0], 'strand']

	# get the path from the transcript id
	def get_path_from_tid(self, tid):
		return self.t_df.loc[tid].path

	# get the gene id from the transcript id
	def get_gid_from_tid(self, tid):
		return self.t_df.loc[tid, 'gid']

	# returns the min and max coordinates of an input gene
	def get_gene_min_max(self, gid):

		# get all starts and stops from gene and their coordinates
		paths = self.t_df.loc[self.t_df.gid == gid].path.tolist()
		starts = np.unique([path[0] for path in paths]).tolist()
		stops = np.unique([path[-1] for path in paths]).tolist()
		coords = self.loc_df.loc[starts+stops, 'coord']

		return int(min(coords)), int(max(coords))

	# returns the min and max coordinates of an input transcript
	def get_transcript_min_max(self, tid):

		path = self.t_df.loc[tid, 'path']
		ends = [path[0], path[-1]]
		end_coords = self.loc_df.loc[ends, 'coord']

		return int(min(end_coords)), int(max(end_coords))


	###########################################################################
	################################## Extras #################################
	##########################################################################

