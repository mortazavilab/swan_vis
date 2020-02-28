import networkx as nx
import numpy as np
import pandas as pd
import copy
from swan_vis.utils import *

# super class that both SwanGraph and PlottedGraph inherit from.
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
	################# Related to checking contents of Graph ##################
	##########################################################################

	# check that input datasets are in the Graph:
	def check_datasets(self, datasets):

		# make sure we have an iterable
		if type(datasets) != list:
			datasets = [datasets]

		g_datasets = self.get_dataset_cols()
		for d in datasets:
			if d not in g_datasets:
				raise Exception('Dataset {} not present in graph. '
								'Datasets in graph are {}'.format(d, g_datasets))


	# check that input datasets have abundance data in the Graph:
	def check_abundances(self, datasets):

		# make sure we have an iterable
		if type(datasets) != list:
			datasets = [datasets]

		ab_cols = self.get_count_cols()
		for d in datasets:
			if '{}_counts'.format(d) not in ab_cols:
				raise Exception('Abundance for dataset {} not present in graph. '
								'Datasets with abundance information '
								'in graph are {}'.format(d, ab_cols))

	# check if there are even any datasets
	def check_if_any_datasets(self, task):
		if self.datasets == None:
			raise Exception('No datasets found in graph. '
				'Cannot perform {}'.format(task))

	# # check if the annotation is in the graph
	# def check_annotation(self):
	# 	datasets = self.get_dataset_cols()

	# check if gid is in SwanGraph
	def check_gene(self, gid):
		if gid not in self.t_df.gid.tolist():
			raise Exception('Gene {} not found in SwanGraph.'.format(gid))

	# check if tid is in SwanGraph
	def check_transcript(self, tid):
		if tid not in self.t_df.tid.tolist():
			raise Exception('Transcript {} not found in SwanGraph.'.format(tid))


	##########################################################################
	####################### Related to creating Graph ########################
	##########################################################################

	# update ids according to coordinates in loc_df, edge_df, and t_df
	def update_ids(self, id_map=None):

		# if the user didn't already input an id map, just order
		# according to genomic position
		if not id_map:
			id_map = self.get_ordered_id_map()

		# convert dfs into dicts for the next steps
		self.dfs_to_dicts()

		# update the ids according to id_map
		self.update_loc_df_ids(id_map)
		self.update_edge_df_ids(id_map)
		self.update_t_df_paths(id_map)

		# convert back to dfs
		self.dicts_to_dfs()	

	# get a dictionary mapping vertex id to ordered new vertex id
	def get_ordered_id_map(self):

		# split loc_df into + and - strand parts
		plus_loc_df = self.loc_df.loc[self.loc_df.strand == '+'].copy(deep=True)
		minus_loc_df = self.loc_df.loc[self.loc_df.strand == '-'].copy(deep=True)

		# sort each of the dfs by chrom, coord either ascending
		# or descending based on strand
		plus_loc_df.sort_values(['chrom', 'coord'],
								 ascending=[True, True],
								 inplace=True)
		minus_loc_df.sort_values(['chrom', 'coord'],
								  ascending=[True, False],
								  inplace=True)

		# concatenate the two dfs
		self.loc_df = pd.concat([plus_loc_df, minus_loc_df])

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

		loc_df = {}
		for old_id, new_id in id_map.items():
			loc_df[new_id] = self.loc_df[old_id]
			loc_df[new_id]['vertex_id'] = new_id
		self.loc_df = loc_df

	# update vertex ids in edge_df
	def update_edge_df_ids(self, id_map):

		edge_df = {}
		for edge_id, item in self.edge_df.items():
			new_id = (id_map[edge_id[0]],id_map[edge_id[1]])
			edge_df[new_id] = item
			edge_df[new_id]['v1'] = new_id[0]
			edge_df[new_id]['v2'] = new_id[1]
			edge_df[new_id]['edge_id'] = new_id
		self.edge_df = edge_df

	# update vertex ids in t_df
	def update_t_df_paths(self, id_map):

		t_df = {}
		for tid, item in self.t_df.items():
			path = item['path']
			new_path = []
			for n in path:
				new_path.append(id_map[n])
			t_df[tid] = item
			t_df[tid]['path'] = new_path
		self.t_df = t_df

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
		self.edge_df.sort_values(by=['v1', 'v2'], inplace=True)

	##########################################################################
	####### Functions to switch back and forth between dfs and dicts #########
	##########################################################################

	# convert loc_df, edge_df, and t_df to dictionaries
	def dfs_to_dicts(self):

		self.loc_df = self.loc_df.to_dict('index')
		self.edge_df = self.edge_df.to_dict('index')
		self.t_df = self.t_df.to_dict('index')

	# convert dictionary versions of loc_df, edge_df, and t_df to dfs
	def dicts_to_dfs(self):

		# loc_df
		self.loc_df = pd.DataFrame.from_dict(self.loc_df, orient='index')
		self.loc_df.index.names = ['vertex_id']

		# pandas interprets the tuple as a multiindex so we need to fix it
		self.edge_df = pd.DataFrame.from_dict(self.edge_df, orient='index')
		self.edge_df.reset_index(drop=True, inplace=True)
		self.edge_df = create_dupe_index(self.edge_df, 'edge_id')
		self.edge_df = set_dupe_index(self.edge_df, 'edge_id')

		# t_df
		self.t_df = pd.DataFrame.from_dict(self.t_df, orient='index')
		self.t_df.index.names = ['tid']

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
	def get_dataset_cols(self, include_annotation=True):
		if include_annotation:
			return self.datasets
		elif not include_annotation:
			if 'annotation' not in self.datasets:
				return self.datasets
			else: 
				datasets = copy.deepcopy(self.datasets)
				datasets.remove('annotation')
				return datasets

	# gets the names of the counts columns in the graph
	# returns None if no counts have been added
	# if datasets option given, returns the counts 
	# columns associated with the input datasets
	def get_count_cols(self, datasets=None):

		if datasets:
			if type(datasets) != list:
				datasets = [datasets]
			self.check_abundances(datasets)
			counts_cols = []
			for d in datasets:
				counts_cols.append('{}_counts'.format(d))
			return counts_cols

		return self.counts

	# gets the names of tpm columns in the graph
	# returns None if no counts have been added
	# if datasets option given, returns the counts 
	# columns associated with the input datasets
	def get_tpm_cols(self, datasets=None):

		if datasets:
			if type(datasets) != list:
				datasets = [datasets]
			self.check_abundances(datasets)
			tpm_cols = []
			for d in datasets:
				tpm_cols.append('{}_tpm'.format(d))
			return tpm_cols

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

	# changes loc_df, edge_df, and t_df into only those seen in gene "gid",
	# recreates the graph using this
	def subset_on_gene(self, gid):

		# make sure this gid is even in the SwanGraph
		if gid not in self.t_df.gid.tolist():
			raise Exception('Gene id {} not found in SwanGraph.'.format(gid))

		# subset t_df first, it's the easiest
		self.t_df = self.t_df.loc[self.t_df.gid == gid]

		# subset loc_df based on all the locs that are in the paths from 
		# the already-subset t_df
		paths = self.t_df['path'].tolist()
		locs = [node for path in paths for node in path]
		locs = np.unique(locs)
		self.loc_df = self.loc_df.loc[locs]

		# subset edge_df based on all the edges that are in the paths from 
		# the alread-subset t_df
		edges = [(v1,v2) for path in paths for v1,v2 in zip(path[:-1],path[1:])]
		edges = list(set(edges))
		self.edge_df = self.edge_df.loc[edges]

		# renumber locs
		self.update_ids()

		# finally create the graph
		self.create_graph_from_dfs()


	###########################################################################
	################################## Extras #################################
	##########################################################################

