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

class SpliceGraph:

	def __init__(self):

		self.datasets = []
		
		loc_df = pd.DataFrame(columns=['chrom', 'coord',
									   'strand','vertex_id',
									   'TSS', 'alt_TSS',
									   'TES', 'alt_TES',
									   'internal'])
		edge_df = pd.DataFrame(columns=['edge_id', 'edge_type',
									    'strand', 'v1', 'v2'])
		t_df = pd.DataFrame(columns=['tid', 'gid',
									 'gname', 'path'])

	# check if anything has been added to the graph yet
	def is_empty(self):
		if len(self.datasets) == 0: 
			return True
		else: 
			return False

	###########################################################################
	############## Related to adding datasets and merging #####################
	###########################################################################

	# add annotation to graph 
	def add_annotation(self, gtf=None, db=None):

		# neither gtf nor db given
		if not gtf and not db: 
			raise Exception('Provide a GTF or TALON db')

		# column name for annotation 
		col = 'annotation'

		# use the add_dataset function to add stuff to graph
		if gtf:
			self.add_dataset(col, gtf=gtf)
		elif db:
			self.add_dataset(col, db=db)

	# add dataset into graph from gtf
	def add_dataset(self, col, gtf=None, db=None):

		# make sure that input dataset name is not
		# already in any of the df col spaces
		if col in self.datsets:
			raise Exception('Dataset {} is already in the graph. '
				'Use update_dataset (coming soon) or provide a different name.')
		if col in self.loc_df.columns:
			raise Exception('Dataset name {} conflicts with preexisting '
				'column in loc_df. Choose a different name.'.format(col))
		if col in self.edge_df.columns:
			raise Exception('Dataset name conflicts with preexisting '
				'column in edge_df. Choose a different name.'.format(col))
		if col in self.t_df.columns:
			raise Exception('Dataset name conflicts with preexisting '
				'column in t_df. Choose a different name.'.format(col))

		# first entry is easy 
		if self.is_empty():

			# get loc_df, edge_df, t_df
			if gtf:
				self.create_dfs_gtf(gtf)
			elif db:
				self.create_dfs_db(db)

			# add column to each df to indicate where data came from
			self.loc_df[col] = True
			self.edge_df[col] = True
			self.t_df[col] = True

			self.update_ids()
			self.get_loc_types()
			self.create_graph_from_dfs()

		# adding a new dataset to the graph requires us to merge
		# SpliceGraph objects
		else:
			temp = SpliceGraph()
			if gtf:
				temp.create_dfs_gtf(gtf)
			elif db:
				temp.create_dfs_db(db)

			self.merge_dfs(temp, col)

		# update graph metadata
		self.datasets.append(col)

	# merge dfs from two SpliceGraph objects
	def merge_dfs(self, b, b_col):

		# loc_df merging/updating 
		self.merge_loc_dfs(b, b_col)
		id_map = self.get_merged_id_map()
		self.update_loc_df_ids_merge(id_map)

		# edge_df merging/updating
		b.update_edge_df_ids_merge(id_map)
		self.merge_edge_dfs(b, b_col)

		# t_df merging/updating
		b.update_t_df_paths_merge(id_map)
		self.merge_t_dfs(b, b_col)

	# update t_df with vertex ids in path for the merged dfs
	def update_t_df_paths_merge(self, id_map):
		self.t_df.path = self.t_df.apply(
			lambda x: [id_map[n] for n in x.path], axis=1)

	# merge t_dfs on tid, gid, gname, path
	def merge_t_dfs(self, b, b_col):

		# some df reformatting
		self.t_df.reset_index(drop=True, inplace=True)
		b.t_df.reset_index(drop=True, inplace=True)
		b.t_df[b_col] = True

		# convert paths to tuples so we can merge on them
		self.t_df.path = self.t_df.apply(
			lambda x: tuple(x.path), axis=1)
		b.t_df.path = b.t_df.apply(
			lambda x: tuple(x.path), axis=1)

		# merge on transcript information
		t_df = self.t_df.merge(b.t_df,
			   how='outer',
			   on=['tid', 'gid', 'gname', 'path'],
			   suffixes=['_a', '_b'])

		# convert path back to list
		t_df.path = t_df.apply(
			lambda x: list(x.path), axis=1)

		# assign False to entries that are not in the new dataset, 
		# and to new entries that were not in the prior datasets
		d_cols = self.datasets+[b_col]
		t_df[d_cols] = t_df[d_cols].fillna(value=False, axis=1)

		self.t_df = t_df

	# update edge_df ids and vertex ids from id_map for the merged dfs
	def update_edge_df_ids_merge(self, id_map):
		self.edge_df.v1 = self.edge_df.apply(
			lambda x: id_map[x.v1], axis=1)
		self.edge_df.v2 = self.edge_df.apply(
			lambda x: id_map[x.v2], axis=1)
		self.edge_df.edge_id = self.edge_df.apply(
			lambda x: (x.v1, x.v2), axis=1)

	# merge edge_dfs on edge_id, v1, v2, strand, edge_type
	def merge_edge_dfs(self, b, b_col):

		# some df reformatting
		self.edge_df.reset_index(drop=True, inplace=True)
		b.edge_df.reset_index(drop=True, inplace=True)
		b.edge_df[b_col] = True

		# merge on edge info
		edge_df = self.edge_df.merge(b.edge_df,
				  how='outer',
				  on=['edge_id', 'v1', 'v2', 'edge_type', 'strand'],
				  suffixes=['_a', '_b'])

		# assign False to entries that are not in the new dataset, 
		# and to new entries that were not in the prior datasets
		d_cols = self.datasets+[b_col]
		edge_df[d_cols] = edge_df[d_cols].fillna(value=False, axis=1)

		self.edge_df = edge_df

	# update loc_df vertex ids from id_map for the merged dfs
	def update_loc_df_ids_merge(self, id_map):
		self.loc_df['vertex_id'] = self.loc_df.apply(
			lambda x: int(x.vertex_id_a) if x.vertex_id_b not in id_map.keys()
									else id_map[x.vertex_id_b], axis=1)
		self.loc_df.drop(['vertex_id_a', 'vertex_id_b'], axis=1, inplace=True)

	# merge loc_dfs on coord, chrom, strand
	def merge_loc_dfs(self, b, b_col):

		# some df reformatting
		node_types = ['TSS', 'alt_TSS', 'TES', 'alt_TES', 'internal']

		self.loc_df.drop(node_types, axis=1, inplace=True)
		self.loc_df.reset_index(drop=True, inplace=True)

		b.loc_df.drop(node_types, axis=1, inplace=True)
		b.loc_df.reset_index(drop=True, inplace=True)
		b.loc_df[b_col] = True

		# merge on location info
		loc_df = self.loc_df.merge(b.loc_df,
				 how='outer',
				 on=['chrom', 'coord', 'strand'],
				 suffixes=['_a','_b'])

		# assign False to entries that are not in the new dataset, 
		# and to new entries that were not in prior datasets
		d_cols = self.datasets+[b_col]
		loc_df[d_cols] = loc_df[d_cols].fillna(value=False, axis=1)

		self.loc_df = loc_df

	# returns a dictionary mapping vertex b: vertex a for each
	# vertex in dataset b
	def get_merged_id_map(self):

		id_map = self.loc_df.apply(
			lambda x: [x.vertex_id_b, x.vertex_id_a], axis=1)

		# loop through id_map and assign new ids for 
		# those present in b but not a
		b_ind = int(self.loc_df.vertex_id_a.max() + 1)
		i = 0
		while i < len(id_map):
			if math.isnan(id_map[i][1]):
				id_map[i][1] = b_ind
				b_ind += 1
			# set up entries where there isn't a b id (entries only found in a)
			# to be removed
			elif math.isnan(id_map[i][0]):
				id_map[i] = []
			i += 1

		# remove entries that are only in a but not in b
		# make sure everything is ints
		id_map = [i for i in id_map if len(i) == 2]
		id_map = dict([(int(a), int(b)) for a,b in id_map])

		return id_map

	##########################################################################
	############# Related to creating dfs from gtf or TALON DB ###############
	##########################################################################

	# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
	def create_dfs_gtf(self, gtf):

		# make sure file exists
		if not os.path.exists(gtf):
			raise Exception('GTF file not found. Check path.')

		# get dfs by parsing through gtf
		loc_df = pd.DataFrame(columns=['chrom', 'coord',
									   'strand','vertex_id',
									   'TSS', 'alt_TSS',
									   'TES', 'alt_TES',
									   'internal'])
		loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)

		edge_df = pd.DataFrame(columns=['edge_id', 'edge_type',
									    'strand', 'v1', 'v2'])
		t_df = pd.DataFrame(columns=['tid', 'gid',
									 'gname', 'path'])

		# loop initialization
		vertex_id = 0
		transcript_paths = []
		transcript_path = []

		with open(gtf, 'r') as infile:
			for line in infile:

				# skip header lines
				if '##' in line: continue

				line = line.strip().split('\t')

				# gene entry 
				if line[2] == 'gene':
					curr_gid = get_field_value('gene_id', line[-1])
					curr_gname = get_field_value('gene_name', line[-1])

				# transcript entry
				elif line[2] == 'transcript':
					curr_tid = get_field_value('transcript_id', line[-1])
					
					# start a new transcript path
					if transcript_path != []:

						# add to list of transcript paths and transcript df 
						transcript_paths.append(transcript_path)
						t_df = t_df.append({'tid': prev_tid,
									 'gid': prev_gid,
									 'gname': prev_gname,
									 'path': transcript_path},
									 ignore_index=True)

					transcript_path = []

					# reset some stuff
					terminal_loc = True
					exon = 0
					intron = 1

				# exon entry
				elif line[2] == 'exon':

					# get exon info 
					curr_chr = line[0]
					curr_start = line[3]
					curr_stop = line[4]
					curr_strand = line[6]
					
					if curr_strand == '+': coords = [curr_start, curr_stop]
					else: coords = [curr_stop, curr_start]
					
					for c in coords:

						ind = (curr_chr, int(c), curr_strand)

						# loc not in loc_df already
						if ind not in loc_df.index.tolist():

							# label as not a TSS/TES until further notice
							attr = {'vertex_id': vertex_id,	   
									'TSS': False, 'TES': False,
									'alt_TSS': False,
									'alt_TES': False, 
									'internal': False, 'coord': int(c),
									'strand': curr_strand, 'chrom': curr_chr}

							# update loc_df and increment vertex_id
							loc_df.reset_index(inplace=True)
							loc_df = loc_df.append(attr, ignore_index=True)
							loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)

							curr_loc = int(vertex_id)
							vertex_id += 1

						# loc was already added to graph
						else: curr_loc = int(loc_df.loc[ind].vertex_id)	
		
						# add an edge to previous loc if not terminal loc 
						# and if the edge doesn't already exist
						if not terminal_loc:
							curr_edge = (prev_loc, curr_loc)
							
							if curr_edge not in edge_df.edge_id.to_list():
								attrs = {'edge_id': (curr_edge[0], curr_edge[1]),
									     'v1': curr_edge[0],
										 'v2': curr_edge[1], 
										 'strand': curr_strand}
								if exon: attrs.update({'edge_type': 'exon'})
								elif intron: attrs.update({'edge_type': 'intron'})

								edge_df = edge_df.append(attrs, ignore_index=True)

						# update transcript path with each loc 
						transcript_path.append(curr_loc)
						prev_loc = curr_loc
						prev_tid = curr_tid
						prev_gid = curr_gid
						prev_gname = curr_gname
						terminal_loc = False
						
						# exon or intron
						exon = abs(exon-1)
						intron = abs(intron-1)
						
		# append last transcript info
		transcript_paths.append(transcript_path)
		t_df = t_df.append({'tid': curr_tid,
						    'gid': curr_gid,
						    'gname': curr_gname,
							'path': transcript_path},
							ignore_index=True)

		# label node/edge types and finish formatting dfs correctly
		loc_df.reset_index(inplace=True)
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')

		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')

		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

	# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
	def create_dfs_db(self, db):

		# make sure file exists
		if not os.path.exists(db):
			raise Exception('TALON db file not found. Check path.')

		# open db connection
		conn = sqlite3.connect(db)
		c = conn.cursor()

		# loc_df
		q = 'SELECT loc.* FROM location loc'

		c.execute(q)
		locs = c.fetchall()

		loc_df = pd.DataFrame(locs,
			columns=['location_ID', 'genome_build',
					 'chrom', 'position'])

		# do some df reformatting, add strand
		loc_df.drop('genome_build', axis=1, inplace=True)
		loc_df.rename({'location_ID': 'vertex_id',
					   'position': 'coord'},
					   inplace=True, axis=1)
		loc_df.vertex_id = loc_df.vertex_id.map(int)

		# edge_df
		q = """SELECT e.* 
				FROM edge e 
				JOIN vertex V ON e.v1=v.vertex_ID 
				JOIN gene_annotations ga ON v.gene_ID=ga.ID 
				WHERE ga.attribute='gene_name'
			""" 

		c.execute(q)
		edges = c.fetchall()

		edge_df = pd.DataFrame(edges, 
			columns=['edge_id', 'v1', 'v2',
					 'edge_type', 'strand'])
		edge_df.v1 = edge_df.v1.map(int)
		edge_df.v2 = edge_df.v2.map(int)
		edge_df['talon_edge_id'] = edge_df.edge_id
		edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

		# t_df
		t_df = pd.DataFrame()

		# get tid, gid, gname, and paths
		q = """SELECT ga.value, ta.value,
					  t.start_exon, t.jn_path, t.end_exon,
					  t.start_vertex, t.end_vertex
				FROM gene_annotations ga 
				JOIN transcripts t ON ga.ID=t.gene_ID
				JOIN transcript_annotations ta ON t.transcript_ID=ta.ID
				WHERE ta.attribute='transcript_id'
				AND (ga.attribute='gene_name' 
				OR ga.attribute='gene_id')
			"""

		c.execute(q)
		data = c.fetchall()

		# get fields from each transcript and add to dataframe
		gids, tids, paths = zip(*[(i[0], i[1], i[2:]) for i in data[::2]])
		gnames = [i[0] for i in data[1::2]]
		paths = self.get_db_edge_paths(paths)

		t_df['tid'] = np.asarray(tids)
		t_df['gid'] = np.asarray(gids)
		t_df['gname'] = np.asarray(gnames)
		t_df['path'] = np.asarray(paths)

		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		# furnish the last bit of info in each df
		t_df['path'] = [[int(n) for n in path]
						 for path in self.get_db_vertex_paths(paths, edge_df)]
		loc_df['strand'] = loc_df.apply(lambda x:
				 self.get_db_strand(x, edge_df), axis=1)
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')
		loc_df['internal'] = False
		loc_df['TSS'] = False
		loc_df['alt_TSS'] = False
		loc_df['TES'] = False
		loc_df['alt_TES'] = False

		edge_df.drop('talon_edge_id', axis=1, inplace=True)
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')

		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

	# add node types (internal, TSS, alt TSS, TES, alt_TES) to loc_df
	def get_loc_types(self):

		loc_df = self.loc_df
		t_df = self.t_df

		# label each location as internal off the bat, and not as TSS/TES
		loc_df['internal'] = False
		loc_df['TSS'] = False
		loc_df['TES'] = False
		loc_df['alt_TSS'] = False
		loc_df['alt_TES'] = False

		# label each TSS and TES
		paths = t_df.path.tolist()
		tss = np.unique([path[0] for path in paths])
		loc_df.loc[tss, 'TSS'] = True
		tes = np.unique([path[-1] for path in paths])
		loc_df.loc[tes, 'TES'] = True
		internal = np.unique([n for path in paths for n in path[1:-1]])
		loc_df.loc[internal, 'internal'] = True

		# label each alt TSS and alt TES for each gene
		for g in t_df.gid.unique().tolist():
			gene_entries = t_df.loc[t_df.gid == g]

			# genes that have more than one transcript are alt TSS/TES candidates
			if len(gene_entries.index) != 1: 

				paths = gene_entries.path.tolist()
				tss = [path[0] for path in paths]
				tes = [path[-1] for path in paths]

				# alt TSS/TES
				if len(set(tss)) > 1: 
					loc_df.loc[tss, 'alt_TSS'] = True
				if len(set(tes)) > 1: 
					loc_df.loc[tes, 'alt_TES'] = True

		self.loc_df = loc_df

	# convert talon query into edge path
	def get_db_edge_paths(self, paths):
		edge_paths = []
		for p in paths:
			if p[1] == None:
				edge_paths.append([p[0]])
			else:
				edge_paths.append(
					[p[0], *[int(i) for i in p[1].split(',')], p[2]])
		return edge_paths

	# convert edge path to vertex path
	def get_db_vertex_paths(self, paths, edge_df):
		vertex_paths = []
		for p in paths: 
			path = []
			for i, e in enumerate(p): 
				entry = edge_df.loc[edge_df.talon_edge_id == e]
				if i == 0:
					path.extend([entry.v1.values[0], entry.v2.values[0]])
				else: path.append(entry.v2.values[0])
			vertex_paths.append(path)
		return vertex_paths

	# get the strand of each vertex
	def get_db_strand(self, x, edge_df):
		# use v1 or v2 depending on where vertex is in edge
		try: 
			strand = edge_df.loc[edge_df.v1 == x.vertex_id, 'strand'].values[0]
		except:
			strand = edge_df.loc[edge_df.v2 == x.vertex_id, 'strand'].values[0]
		return strand

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

		# add node attributes from dfs
		G = label_nodes(G, self.loc_df, 'internal', 'internal') 
		G = label_nodes(G, self.loc_df, 'TSS', 'TSS') 
		G = label_nodes(G, self.loc_df, 'alt_TSS', 'alt_TSS') 
		G = label_nodes(G, self.loc_df, 'TES', 'TES')
		G = label_nodes(G, self.loc_df, 'alt_TES', 'alt_TES')
		G = label_nodes(G, self.loc_df, 'coord', 'coord')
		G = label_edges(G, self.edge_df, 'strand', 'strand')
		G = label_edges(G, self.edge_df, 'edge_type', 'edge_type')

		self.G = G

	# # removes the dataset columns from a preexisting SpliceGraph's
	# # dfs ONLY TODO do I want to add support for removing labels from nodes
	# # in the future?
	# def remove_dataset(self, name):

	# 	d_col = 'dataset_'+name
	# 	if d_col not in self.loc_df.columns or d_col not in self.edge_df.columns or d_col not in self.t_df.columns:
	# 		raise Exception('Dataset {} not in the graph'.format(name))

	# 	# remove this column from all of the dfs
	# 	self.loc_df.drop(d_col, axis=1, inplace=True)
	# 	self.edge_df.drop(d_col, axis=1, inplace=True)
	# 	self.t_df.drop(d_col, axis=1, inplace=True)

	# # adds abundance information to the t_df from an input abundance file 
	# # where rows are transcript ids and columns are datasets
	# def add_abundance_dataset(self, file, count_cols, dataset_name):

	# 	# get the counts from the input abundance file
	# 	counts = process_abundance_file(file, count_cols)
	# 	counts.rename({'counts': 'counts_{}'.format(dataset_name)},
	# 				   axis=1, inplace=True)

	# 	# merge on tid and format t_df as necessary
	# 	self.t_df.reset_index(drop=True, inplace=True)
	# 	self.t_df = self.t_df.merge(counts, on='tid', how='left')
	# 	self.t_df.fillna(value=0, inplace=True)
	# 	self.t_df = create_dupe_index(self.t_df, 'tid')
	# 	self.t_df = set_dupe_index(self.t_df, 'tid')

	# # order the transcripts by expression of transcript, transcript id, 
	# # or start/end nodes
	# def order_transcripts(self, order='tid'):

	# 	# order by transcript id
	# 	if order == 'tid':
	# 		ordered_tids = sorted(self.t_df.tid.tolist())
	# 		self.t_df = self.t_df.loc[ordered_tids]

	# 	# order by expression
	# 	elif order == 'expression':
	# 		count_fields = get_count_fields(self.t_df)

	# 		# make sure there are counts in the graph at all
	# 		if count_fields:
	# 			self.t_df['counts_sum'] = self.t_df.apply(lambda x:
	# 				sum(x[count_fields]), axis=1)
	# 			self.t_df.sort_values(by='counts_sum', 
	# 								  ascending=False, 
	# 								  inplace=True)
	# 			self.t_df.drop('counts_sum', axis=1, inplace=True)

	# 	# order by coordinate of tss
	# 	# TODO might be able to roll this in with getting ordered_nodes
	# 	# in PlottedGraph
	# 	elif order == 'tss':
	# 		self.t_df['start_coord'] = self.t_df.apply(lambda x: 
	# 			self.loc_df.loc[x.path[0], 'coord'], axis=1)

	# 		# watch out for strandedness
	# 		if self.loc_df.loc[self.loc_df.index[0], 'strand'] == '-':
	# 			ascending = False
	# 		else: 
	# 			ascending = True
	# 		self.t_df.sort_values(by='start_coord',
	# 							  ascending=ascending,
	# 							  inplace=True)
	# 		self.t_df.drop('start_coord', axis=1, inplace=True)
			
	# 	# order by coordinate of tes
	# 	elif order == 'tes':
	# 		self.t_df['end_coord'] = self.t_df.apply(lambda x: 
	# 			self.loc_df.loc[x.path[-1], 'coord'], axis=1)

	# 		# watch out for strandedness
	# 		if self.loc_df.loc[self.loc_df.index[0], 'strand'] == '-':
	# 			ascending = True
	# 		else: 
	# 			ascending = False
	# 		self.t_df.sort_values(by='end_coord',
	# 							  ascending=ascending,
	# 							  inplace=True)
	# 		self.t_df.drop('end_coord', axis=1, inplace=True)

# # adds graph b with dataset name bname to graph a, which does not 
# # require a completely new column
# def add_graph(a, b, bname):
# 	sg = merge_graphs(a,b,None,bname,store_a=False)
# 	return sg

# # merges two splice graph objects and creates a third
# # store_a is to support just adding b to graph a via 
# # add_graph
# def merge_graphs(a, b, aname, bname, store_a=True):

# 	# remove indices for each df
# 	a.loc_df.reset_index(drop=True, inplace=True)
# 	a.edge_df.reset_index(drop=True, inplace=True)
# 	a.t_df.reset_index(drop=True, inplace=True)
# 	b.loc_df.reset_index(drop=True, inplace=True)
# 	b.edge_df.reset_index(drop=True, inplace=True)
# 	b.t_df.reset_index(drop=True, inplace=True)

# 	# dataset columns
# 	if not store_a:
# 		a_col = 'dataset_TEMP'
# 	else:
# 		a_col = 'dataset_'+aname
# 	b_col = 'dataset_'+bname

# 	# merge loc_dfs based on chrom, coord, strand and update vertex ids
# 	loc_df = merge_loc_dfs(a.loc_df, b.loc_df, a_col, b_col)
# 	id_map = get_vertex_id_map(loc_df, a_col, b_col)
# 	loc_df = assign_new_vertex_ids(loc_df, id_map)
# 	loc_df['vertex_id'] = loc_df['vertex_id'].astype(int)

# 	# merge edge_df based on new vertex ids 
# 	b.edge_df = assign_new_edge_ids(b.edge_df, id_map)
# 	edge_df = merge_edge_dfs(a.edge_df, b.edge_df, a_col, b_col)
	
# 	# merge t_df based on new vertex ids
# 	b.t_df= assign_new_paths(b.t_df, id_map)
# 	t_df = merge_t_dfs(a.t_df, b.t_df, a_col, b_col)

# 	# final df formatting
# 	loc_df.drop(['vertex_id_a', 'vertex_id_b'], inplace=True, axis=1)
# 	loc_df = create_dupe_index(loc_df, 'vertex_id')
# 	loc_df = set_dupe_index(loc_df, 'vertex_id')
# 	loc_df = get_loc_types(loc_df, t_df)

# 	edge_df = create_dupe_index(edge_df, 'edge_id')
# 	edge_df = set_dupe_index(edge_df, 'edge_id')

# 	t_df = create_dupe_index(t_df, 'tid')
# 	t_df = set_dupe_index(t_df, 'tid')

# 	# finally, let's make the merged graph
# 	sg = SpliceGraph(loc_df=loc_df, edge_df=edge_df, t_df=t_df)

# 	# if we didn't want to make a new dataset column for a, 
# 	# remove it here. 
# 	# otherwise, we're free to label nodes/edges with new 
# 	# dataset information
# 	if not store_a:
# 		sg.remove_dataset('TEMP')
# 	# else:
# 	# 	sg.G = label_nodes(sg.G, loc_df, a_col, a_col)
# 	# 	sg.G = label_edges(sg.G, edge_df, a_col, a_col)

# 	# assign labels to nodes and edges based on what dataset they came from
# 	d_fields = get_dataset_fields(df=loc_df)
# 	for d_field in d_fields:
# 		sg.G = label_nodes(sg.G, loc_df, d_field, d_field)
# 		sg.G = label_edges(sg.G, edge_df, d_field, d_field)

# 	# # 
# 	# sg.merged = True

# 	return sg

# # merge transcript dfs on tid, gid, gname, and path
# def merge_t_dfs(a,b,a_col,b_col):

# 	# track which datasets each transcript is in 
# 	a[a_col] = True
# 	b[b_col] = True

# 	# first convert paths to tuples so we can merge on them
# 	a.path = a.apply(lambda x: tuple(x.path), axis=1)
# 	b.path = b.apply(lambda x: tuple(x.path), axis=1)

# 	# merge on path ids as well as transcript-associated names and ids
# 	t_df = a.merge(b, 
# 			how='outer',
# 			on=['tid', 'gid', 'gname', 'path'], 
# 			suffixes=['_a', '_b'])

# 	# convert back to lists for path
# 	t_df.path = t_df.apply(lambda x: list(x.path), axis=1)

# 	# assign False to entries that are not in one dataset or another 
# 	d_fields = get_dataset_fields(df=t_df)
# 	t_df[d_fields] = t_df[d_fields].fillna(value=False, axis=1)

# 	return t_df

# # 
# def assign_new_paths(b, id_map):
# 	b.path = b.apply(lambda x: [id_map[n] for n in x.path], axis=1)
# 	return b

# # merge the edge dfs
# def merge_edge_dfs(a, b, a_col, b_col):

# 	# add column that we can track which dataset this comes from
# 	a[a_col] = True
# 	b[b_col] = True

# 	# merge on edge_id as these have already been updated in merge_graphs
# 	edge_df = a.merge(b,
# 			how='outer',
# 			on=['edge_id','v1','v2','edge_type','strand'],
# 			suffixes=['_a', '_b'])

# 	# assign False to entries that are not in one dataset or another 
# 	d_fields = get_dataset_fields(df=edge_df)
# 	edge_df[d_fields] = edge_df[d_fields].fillna(value=False, axis=1)

# 	return edge_df

# # replace vertex ids according to new values in loc_df 
# def assign_new_edge_ids(b, id_map):
# 	b.v1 = b.apply(lambda x: id_map[x.v1], axis=1)
# 	b.v2 = b.apply(lambda x: id_map[x.v2], axis=1)
# 	b.edge_id = b.apply(lambda x: (x.v1, x.v2), axis=1)
# 	return b

# # merge loc_dfs
# def merge_loc_dfs(a, b, a_col, b_col):

# 	# add column that we can track which dataset this comes from
# 	a[a_col] = True
# 	b[b_col] = True

# 	# remove all node type columns as these will be recomputed
# 	node_types = ['TSS', 'alt_TSS', 'TES', 'alt_TES', 'internal']
# 	a.drop(node_types, axis=1, inplace=True)
# 	b.drop(node_types, axis=1, inplace=True)

# 	# merge on location info
# 	loc_df = a.merge(b,
# 			how='outer',
# 			on=['chrom', 'coord', 'strand'],
# 			suffixes=['_a','_b'])

# 	# assign False to entries that are not in one dataset or another 
# 	d_fields = get_dataset_fields(df=loc_df)
# 	loc_df[d_fields] = loc_df[d_fields].fillna(value=False, axis=1)

# 	return loc_df

# def assign_new_vertex_ids(df, id_map):
# 	# id_map = get_vertex_id_map(df)
# 	df['vertex_id'] = df.apply(lambda x: x.vertex_id_a
# 						 if x.vertex_id_b not in id_map.keys()
# 						 else id_map[x.vertex_id_b], axis=1)
# 	return df

# # determines which dfs this entry was present in before the merge
# def present_in(x):
# 	# # a merge has been done before
# 	# if 'present_in' in x.columns:
# 	# 	datasets = x.present_in
# 	# else:
# 	# 	datasets = []
# 	datasets = []
# 	if x.present_in_a == True:
# 		datasets.append('a')
# 	if x.present_in_b == True:
# 		datasets.append('b')
# 	return datasets

# # returns the mapping of vertices from b to their new ids
# def get_vertex_id_map(df, a_col, b_col):

# 	# vertices in both graph a and b 
# 	ab_ids = df.apply(lambda x: (int(x.vertex_id_b), int(x.vertex_id_a))
# 						if x[a_col] and x[b_col]
# 						else np.nan, axis=1)
# 	ab_ids = [ab_id for ab_id in ab_ids if type(ab_id) == tuple]
# 	vertex_id_map = dict(ab_ids)

# 	# vertices only in graph b
# 	b_ids = df.apply(lambda x: x.vertex_id_b
# 						if not x[a_col] and x[b_col]
# 						else np.nan, axis=1)
# 	b_ids = [int(b_id) for b_id in b_ids if not math.isnan(b_id)]

# 	# new ids for these guys
# 	start_b_id = int(df.vertex_id_a.max()+1)
# 	new_b_ids = [int(i) for i in range(start_b_id, len(b_ids)+start_b_id)]
# 	vertex_id_map.update(dict(zip(b_ids, new_b_ids)))

# 	return vertex_id_map

# # add node types (internal, TSS, alt TSS, TES, alt_TES) to loc_df
# def get_loc_types(loc_df, t_df):

# 	# label each location as internal off the bat, and not as TSS/TES
# 	loc_df['internal'] = False
# 	loc_df['TSS'] = False
# 	loc_df['TES'] = False
# 	loc_df['alt_TSS'] = False
# 	loc_df['alt_TES'] = False

# 	# label each TSS and TES
# 	paths = t_df.path.tolist()
# 	tss = np.unique([path[0] for path in paths])
# 	loc_df.loc[tss, 'TSS'] = True
# 	tes = np.unique([path[-1] for path in paths])
# 	loc_df.loc[tes, 'TES'] = True
# 	internal = np.unique([n for path in paths for n in path[1:-1]])
# 	loc_df.loc[internal, 'internal'] = True

# 	# label each alt TSS and alt TES for each gene
# 	for g in t_df.gid.unique().tolist():
# 		gene_entries = t_df.loc[t_df.gid == g]

# 		# genes that have more than one transcript are alt TSS/TES candidates
# 		if len(gene_entries.index) != 1: 

# 			paths = gene_entries.path.tolist()
# 			tss = [path[0] for path in paths]
# 			tes = [path[-1] for path in paths]

# 			# alt TSS/TES
# 			if len(set(tss)) > 1: 
# 				loc_df.loc[tss, 'alt_TSS'] = True
# 			if len(set(tes)) > 1: 
# 				loc_df.loc[tes, 'alt_TES'] = True

# 	return loc_df

# # returns the fields in a graph that specify which dataset a node or
# # edge belongs to in a merged graph
# def get_dataset_fields(graph=None, df=None):
# 	# TODO does graph have datasets? If not throw an error
# 	if graph is not None:
# 		data = graph.nodes(data=True)[0]
# 		d_fields = [k for k in data.keys() if 'dataset_' in k]
# 	if df is not None:
# 		d_fields = [col for col in df.columns if 'dataset_' in col]
# 	return d_fields

# # returns the fields in the t_df that hold counts of datasets
# def get_count_fields(t_df):
# 	c_fields = [col for col in t_df.columns if 'counts_' in col]
# 	return c_fields

# # returns the (min, max) coordinates of an input gene
# def get_gene_min_max(loc_df, t_df, gid):

# 	# all transcripts from this gene
# 	paths = t_df.loc[t_df.gid == gid].path.tolist()
# 	starts = np.unique([path[0] for path in paths]).tolist()
# 	stops = np.unique([path[-1] for path in paths]).tolist()
# 	coords = loc_df.loc[starts + stops, 'coord']

# 	return int(min(coords)), int(max(coords))



