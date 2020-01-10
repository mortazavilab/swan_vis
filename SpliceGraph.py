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
from PlottedGraph import PlottedGraph 
from Graph import Graph
from Report import *
import time

class SpliceGraph(Graph):

	def __init__(self):
		super().__init__()

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
	def add_dataset(self, col, gtf=None, db=None,
					counts_file=None, count_cols=None):

		# make sure that input dataset name is not
		# already in any of the df col spaces
		if col in self.datasets:
			raise Exception('Dataset {} is already in the graph. '
				'Use update_dataset (coming soon) or provide a different name.'.format(col))
		if col in self.loc_df.columns:
			raise Exception('Dataset name {} conflicts with preexisting '
				'column in loc_df. Choose a different name.'.format(col))
		if col in self.edge_df.columns:
			raise Exception('Dataset name {} conflicts with preexisting '
				'column in edge_df. Choose a different name.'.format(col))
		if col in self.t_df.columns:
			raise Exception('Dataset name {} conflicts with preexisting '
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

		# adding a new dataset to the graph requires us to merge
		# SpliceGraph objects
		else:
			temp = SpliceGraph()
			if gtf:
				temp.create_dfs_gtf(gtf)
			elif db:
				temp.create_dfs_db(db)
			self.merge_dfs(temp, col)

		# order node ids by genomic position, add node types,
		# and create graph
		self.update_ids()
		self.order_edge_df()
		self.get_loc_types()
		self.create_graph_from_dfs()

		# update graph metadata
		self.datasets.append(col)

		# if we're also adding abundances
		if counts_file and count_cols:
			self.add_abundance(counts_file, count_cols, col)

	# adds counts columns to t_df based on columns counts_cols found in 
	# tsv counts_file. Relies on column annot_transcript_id to have the
	# transcript id (tid) that is used to index t_df. 
	# TODO make this more flexible in the future
	def add_abundance(self, counts_file, count_cols, dataset_name):

		# if the dataset we're trying to add counts too doesn't exist
		if dataset_name not in self.datasets:
			raise Exception('Trying to add expression data to a dataset '
							'that is not in the graph. Add dataset to graph first.')

		# get counts from input abundance file 
		abundance_df = process_abundance_file(counts_file, count_cols)
		abundance_df.rename({'tpm': '{}_tpm'.format(dataset_name),
							 'counts': '{}_counts'.format(dataset_name)},
							 axis=1, inplace=True)

		# merge on transcript id (tid) with t_df and make sure it's 
		# formatted correctly
		self.t_df.reset_index(drop=True, inplace=True)
		self.t_df = self.t_df.merge(abundance_df, on='tid', how='left')
		self.t_df.fillna(value=0, inplace=True)
		self.t_df = create_dupe_index(self.t_df, 'tid')
		self.t_df = set_dupe_index(self.t_df, 'tid')

		# finally update object's metadata
		self.counts.append('{}_counts'.format(dataset_name))
		self.tpm.append('{}_tpm'.format(dataset_name))

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

		# set up index again
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

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
	############# Related to creating dfs from GTF or TALON DB ###############
	##########################################################################

	# create loc_df (nodes), edge_df (edges), and t_df (transcripts) from gtf
	# adapted from Dana Wyman and TALON
	def create_dfs_gtf(self, gtf_file):

		# make sure file exists
		if not os.path.exists(gtf_file):
			raise Exception('GTF file not found. Check path.')

		# depending on the strand, determine the stard and stop
		# coords of an intron or exon
		def find_edge_start_stop(v1, v2, strand):
			if strand == '-':
				start = max([v1, v2])
				stop = min([v1, v2])
			elif strand == '+':
				start = min([v1, v2])
				stop = max([v1, v2])
			return start, stop

		# dictionaries to hold unique edges and transcripts
		transcripts = {}
		edges = {}

		start_time = time.time()
		line_num = 0 
		with open(gtf_file) as gtf:
			for line in gtf:

				line_num += 1
				if line_num == 200000:
					break

				# ignore header lines
				if '##' in line:
					continue

				# split each entry
				line = line.strip().split('\t')

				# get some fields from gtf that we care about
				chrom = line[0]
				entry_type = line[2]
				start = line[3]
				stop = line[4]
				strand = line[6]
				fields = line[-1]

				# transcript entry 
				if entry_type == "transcript":
					attributes = get_fields(fields)
					tid = attributes['transcript_id']
					gid = attributes['gene_id']
					gname = attributes['gene_name']

					# add transcript to dictionary 
					transcript = {tid: {'gid': gid,
										'gname': gname,
										'tid': tid,
										'strand': strand,
										'exons': []}}
					transcripts.update(transcript)
					
				# exon entry
				elif entry_type == "exon":
					attributes = get_fields(fields)
					start, stop = find_edge_start_stop(start, stop, strand)
					eid = '{}_{}_{}_{}_exon'.format(chrom, start, stop, strand)
					tid = attributes['transcript_id']	

					# add novel exon to dictionary 
					if eid not in edges:
						edge = {eid: {'eid': eid,
									  'chrom': chrom,
									  'v1': int(start),
									  'v2': int(stop),
									  'strand': strand,
									  'edge_type': 'exon'}}
						edges.update(edge)
			   
			   		# add this exon to the transcript's list of exons
					if tid in transcripts:
						transcripts[tid]['exons'].append(eid)

		# gets the vertex id associated with chrom, coord, and strand info
		def get_vertex_id_from_loc(chrom, coord, strand, loc_df):
			v = loc_df.loc[(chrom, coord, strand), 'vertex_id']
			return v

		print('time to create dictionaries')
		print("--- %s seconds ---" % (time.time() - start_time))
		start_time = time.time()

		# once we have all transcripts, make loc_df
		locs = []
		for edge_id, edge in edges.items():
			chrom = edge['chrom']
			strand = edge['strand']

			v1 = edge['v1']
			v2 = edge['v2']

			# exon start
			locs.append({'chrom': chrom,
						 'coord': v1,
						 'strand': strand})
			# exon end
			locs.append({'chrom': chrom,
						 'coord': v2,
						 'strand': strand})

		loc_df = pd.DataFrame(locs)
		loc_df.drop_duplicates(['chrom', 'coord', 'strand'], inplace=True)

		# sort loc_df based on genomic position
		strand = loc_df.loc[loc_df.index[0], 'strand']
		if strand == '+':
			loc_df.sort_values(by='coord', 
							   ascending=True,
							   inplace=True)
		elif strand == '-':
			loc_df.sort_values(by='coord', 
							   ascending=False, 
							   inplace=True)

		# set index to make querying table easier for populating 
		# t_df and edge_df
		loc_df['vertex_id'] = [i for i in range(len(loc_df.index))]
		loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)

		# testing
		loc_df.to_csv('loc_df_test.csv')

		# add loc_df-indexed path to transcripts
		for _, t in transcripts.items():
			t['path'] = []
			for e in t['exons']:
				exon = edges[e]
				chrom = exon['chrom']
				v1 = exon['v1']
				v2 = exon['v2']
				strand = exon['strand']

				v1 = get_vertex_id_from_loc(chrom, v1, strand, loc_df)
				v2 = get_vertex_id_from_loc(chrom, v2, strand, loc_df)

				t['path'].append(v1)
				t['path'].append(v2)

		# add introns by iterating through transcripts
		for tid, t in transcripts.items():
			strand = t['strand']
			if strand == '-':
				t['exons'].reverse()

			# each exon-exon junction constitutes an intron
			for i in range(len(t['exons'])-1):				

				intron = [t['exons'][i], t['exons'][i+1]]
				v1 = edges[intron[0]]['v2']
				v2 = edges[intron[1]]['v1']
				eid = '{}_{}_{}_{}_intron'.format(chrom, v1, v2, strand)

				# if we haven't seen this intron before, add it
				if eid not in edges:
					edge = {eid: {'eid': eid,
							'chrom': chrom,
							'v1': v1,
							'v2': v2,
							'strand': strand,
							'edge_type': 'intron'}}
					edges.update(edge)

		# finish populating/formatting edge_df
		edges = [edge for eid, edge in edges.items()]
		edge_df = pd.DataFrame(edges)
		edge_df.drop_duplicates(['chrom', 'v1', 'v2', 'strand', 'edge_type'], inplace=True)
		# edge_df['v1'] = edge_df['v1'].apply(pd.to_numeric)
		# edge_df['v2'] = edge_df['v2'].apply(pd.to_numeric)
		# edge_df['v1'] = pd.to_numeric(edge_df['v1'])
		# edge_df['v2'] = pd.to_numeric(edge_df['v2'])

		print(edge_df.dtypes)

		# replace v1 and v2 with the actual vertex ids in edge_df
		# and create the edge ids
		edge_df['new_v1'] = edge_df.apply(lambda x:
			get_vertex_id_from_loc(x.chrom, x.v1, x.strand, loc_df), axis=1)
		edge_df['new_v2'] = edge_df.apply(lambda x: 
			get_vertex_id_from_loc(x.chrom, x.v2, x.strand, loc_df), axis=1)
		edge_df.drop(['v1','v2', 'eid', 'chrom'], axis=1, inplace=True)
		edge_df.rename({'new_v1': 'v1', 'new_v2': 'v2'}, axis=1, inplace=True)
		edge_df['edge_id'] = edge_df.apply(lambda x: (x.v1, x.v2), axis=1)
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')

		# finish populating/formatting t_df
		transcripts = [transcript for tid, transcript in transcripts.items()]
		t_df = pd.DataFrame(transcripts)
		t_df.drop('exons', axis=1, inplace=True)
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		# finish formatting loc_df
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')


		print(loc_df.head())
		print(edge_df.head())
		print(t_df.head())

		# for t in transcripts:
		# 	print(transcripts[t]['path'])

		print("--- %s seconds ---" % (time.time() - start_time))

		return transcripts, edges

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

	##########################################################################
	######################## Other SpliceGraph utilities #####################
	##########################################################################

	# order the transcripts by expression of transcript, transcript id, 
	# or start/end nodes
	def order_transcripts(self, order='tid'):

		# order by transcript id
		if order == 'tid':
			ordered_tids = sorted(self.t_df.tid.tolist())
			self.t_df = self.t_df.loc[ordered_tids]

		# order by expression
		elif order == 'expression':
			tpm_cols = self.get_tpm_cols()

			# make sure there are counts in the graph at all
			if tpm_cols:
				self.t_df['tpm_sum'] = self.t_df.apply(lambda x:
					sum(x[tpm_cols]), axis=1)
				self.t_df.sort_values(by='tpm_sum', 
									  ascending=False, 
									  inplace=True)
				self.t_df.drop('tpm_sum', axis=1, inplace=True)
			else: 
				raise Exception('Cannot order by expression because '
								'there is no expression data.')

		# order by coordinate of tss
		# TODO might be able to roll this in with getting ordered_nodes
		# in PlottedGraph
		elif order == 'tss':
			self.t_df['start_coord'] = self.t_df.apply(lambda x: 
				self.loc_df.loc[x.path[0], 'coord'], axis=1)

			# watch out for strandedness
			if self.loc_df.loc[self.loc_df.index[0], 'strand'] == '-':
				ascending = False
			else: 
				ascending = True
			self.t_df.sort_values(by='start_coord',
								  ascending=ascending,
								  inplace=True)
			self.t_df.drop('start_coord', axis=1, inplace=True)
			
		# order by coordinate of tes
		elif order == 'tes':
			self.t_df['end_coord'] = self.t_df.apply(lambda x: 
				self.loc_df.loc[x.path[-1], 'coord'], axis=1)

			# watch out for strandedness
			if self.loc_df.loc[self.loc_df.index[0], 'strand'] == '-':
				ascending = False
			else: 
				ascending = True
			self.t_df.sort_values(by='end_coord',
								  ascending=ascending,
								  inplace=True)
			self.t_df.drop('end_coord', axis=1, inplace=True)

	##########################################################################
	############################ Plotting utilities ##########################
	##########################################################################


	# plot the SpliceGraph object according to the user's input
	def plot_graph(self, combine=False,
				   indicate_dataset=False,
				   indicate_novel=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel)

		# TODO check to see if we already have a pg object and 
		# reinitialize based on compatibility of object 

		# create PlottedGraph object and plot summary graph
		self.pg = PlottedGraph(self, combine, indicate_dataset, indicate_novel)
		self.pg.plot_graph()


	# plot an input transcript's path through the summary graph 
	def plot_transcript_path(self, tid, combine=False,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)

		# create PlottedGraph object
		self.pg = PlottedGraph(self,
							   combine,
							   indicate_dataset,
							   indicate_novel,
							   tid=tid,
							   browser=browser)

		self.pg.plot_graph()

	# plots each transcript's path through the summary graph, and automatically saves them!
	def plot_each_transcript(self, gid, prefix, combine=False,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)

		# make sure this gid is even in the SpliceGraph
		if gid not in self.t_df.gid.tolist():
			raise Exception('Gene id {} not found in SpliceGraph.'.format(gid))

		# loop through each transcript in the SpliceGraph object
		for tid in self.t_df.loc[self.t_df.gid == gid, 'tid'].tolist():
			self.pg = PlottedGraph(self,
								   combine,
								   indicate_dataset,
								   indicate_novel,
								   tid=tid,
								   gid=gid,
								   browser=browser)
			fname = create_fname(prefix,
								 combine,
								 indicate_dataset,
								 indicate_novel,
								 browser,
								 tid=tid)
			self.pg.plot_graph()
			self.save_fig(fname)

	# saves current figure named oname. clears the figure space so additional
	# plotting can be done
	def save_fig(self, oname):
		plt.axis('off')
		plt.tight_layout()
		plt.savefig(oname, format='png', dpi=200)
		plt.clf()
		plt.close()

	##########################################################################
	############################### Report stuff #############################
	##########################################################################

	# creates a report for each transcript model for a gene according to user input
	def gen_report(self,
				   gid,
				   prefix,
				   datasets='all',
				   tpm=False,
				   include_unexpressed=False,
				   combine=False,
				   indicate_dataset=False, 
				   indicate_novel=False,
				   browser=False,
				   order='expression'):

		# make sure all input datasets are present in graph
		if datasets == 'all':
			datasets = copy.deepcopy(self.datasets)
			datasets.remove('annotation')
		elif not datasets:
			datasets = []
		self.check_datasets(datasets)

		# now check to make sure abundance data is there for the
		# query columns, if user is asking
		if tpm:
			self.check_abundances(datasets)
			report_cols = self.get_tpm_cols(datasets)
		elif datasets:
			report_cols = datasets

		# order transcripts by user's preferences
		self.order_transcripts(order)

		# if user doesn't care about datasets, just show all transcripts
		if not datasets:
			include_unexpressed = True
		# user only wants transcript isoforms that appear in their data
		if not include_unexpressed:

			## TODO could also use the all(nonzero) strategy
			counts_cols = self.get_count_cols(datasets)
			self.t_df['counts_sum'] = self.t_df.apply(lambda x: sum(x[counts_cols]), axis=1)
			report_tids = self.t_df.loc[(self.t_df.gid==gid)&(self.t_df.counts_sum>0),
						  	'tid'].tolist()
			self.t_df.drop('counts_sum', inplace=True, axis=1)
		else:
			report_tids = self.t_df.loc[self.t_df.gid == gid, 'tid'].tolist()

		# make sure our swan/browser graph args work together
		# and plot each transcript with these settings
		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)
		self.plot_each_transcript(gid, prefix, combine,
								  indicate_dataset,
								  indicate_novel,
								  browser=browser)

		# if we're plotting tracks, we need a scale as well
		if not browser:
			report_type = 'swan'
		else:
			self.pg.plot_browser_scale()
			self.save_fig(prefix+'_browser_scale.png')
			report_type = 'browser'

		# create report
		pdf_name = create_fname(prefix, 
					 combine,
					 indicate_dataset,
					 indicate_novel,
					 browser,
					 ftype='report',
					 gid=gid)
		report = Report(prefix, report_type, report_cols)
		report.add_page()

		# loop through each transcript and add it to the report
		for tid in report_tids:
			entry = self.t_df.loc[tid]
			## TODO would be faster if I didn't have to compute these names twice....
			## ie once in plot_each_transcript and once here
			fname = create_fname(prefix,
								 combine,
								 indicate_dataset,
								 indicate_novel, 
								 browser,
								 tid=entry.tid)
			report.add_transcript(entry, fname)
		report.write_pdf(pdf_name)

	##########################################################################
	############################# Error handling #############################
	##########################################################################

	# make sure that the set of arguments work with each other 
	# before we start plotting
	def check_plotting_args(self, combine, indicate_dataset, indicate_novel, browser=False):

		# can only do one or another
		if indicate_dataset and indicate_novel:
			raise Exception('Please choose either indicate_dataset '
							'or indicate_novel, not both.')

		# if indicate_dataset or indicate_novel are chosen, make sure
		# the dataset or annotation data exists in the SpliceGraph
		if indicate_novel and 'annotation' not in self.get_dataset_cols():
			raise Exception('Annotation data not present in graph. Use  '
							'add_annotation before using indicate_novel')
		if indicate_dataset and indicate_dataset not in self.get_dataset_cols():
			raise Exception('Dataset {} not present in the graph. '
							''.format(indicate_dataset))

		# if browser, can't do indicate_novel, combine, or indicate_dataset
		if browser:
			if indicate_novel or indicate_dataset:
				raise Exception('Currently highlighting splice junctions in '
								'browser form is unsupported. Use browser option '
								'without inidicate_novel or inidicate_dataset.')
			if combine:
				raise Exception('Combining non-branching splice junctions '
								'not possible for browser track plotting.')


##########################################################################
################################## Extras ################################
##########################################################################





