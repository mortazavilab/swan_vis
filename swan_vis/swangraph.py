
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import math
import copy
from collections import defaultdict
import sqlite3
import pickle
from swan_vis.utils import *
from swan_vis.graph import Graph
from swan_vis.plottedgraph import PlottedGraph
from swan_vis.report import Report

class SwanGraph(Graph):

	def __init__(self):
		super().__init__()

	###########################################################################
	############## Related to adding datasets and merging #####################
	###########################################################################

	# add annotation to graph 
	def add_annotation(self, fname):

		# column name for annotation 
		col = 'annotation'

		# use the add_dataset function to add stuff to graph
		self.add_dataset(col, fname)

	# add dataset into graph from gtf
	def add_dataset(self, col, fname, dname=None,
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

		# are we dealing with a gtf or a db?
		ftype = gtf_or_db(fname)

		# first entry is easy 
		if self.is_empty():

			# get loc_df, edge_df, t_df
			if ftype == 'gtf':
				self.create_dfs_gtf(fname)
			elif ftype == 'db':
				self.create_dfs_db(fname, dname)

			# add column to each df to indicate where data came from
			self.loc_df[col] = True
			self.edge_df[col] = True
			self.t_df[col] = True

		# adding a new dataset to the graph requires us to merge
		# SwanGraph objects
		else:
			temp = SwanGraph()
			if ftype == 'gtf':
				temp.create_dfs_gtf(fname)
			elif ftype == 'db':
				temp.create_dfs_db(fname, dname)
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

	# merge dfs from two SwanGraph objects
	def merge_dfs(self, b, b_col):

		# merge loc dfs
		# what locations correspond between the datasets?
		self.merge_loc_dfs(b, b_col)
		id_map = self.get_merged_id_map()

		self.loc_df.drop(['vertex_id_a','vertex_id_b'], axis=1, inplace=True)
		self.loc_df['vertex_id'] = self.loc_df.index
		self.loc_df = create_dupe_index(self.loc_df, 'vertex_id')
		self.loc_df = set_dupe_index(self.loc_df, 'vertex_id')
		b.loc_df = create_dupe_index(b.loc_df, 'vertex_id')
		b.loc_df = set_dupe_index(b.loc_df, 'vertex_id')

		# update the ids in b to make edge_df, t_df merging easier
		b.update_ids(id_map=id_map)

		# merge edge_df and t_df
		self.merge_edge_dfs(b, b_col)
		self.merge_t_dfs(b, b_col)

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

		# remake index
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')
		
		self.edge_df = edge_df

	# merge loc_dfs on coord, chrom, strand
	def merge_loc_dfs(self, b, b_col):

		# some df reformatting
		node_types = ['TSS', 'alt_TSS', 'TES', 'alt_TES', 'internal']

		self.loc_df.drop(node_types, axis=1, inplace=True)
		self.loc_df.reset_index(drop=True, inplace=True)

		# b.loc_df.drop(node_types, axis=1, inplace=True)
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
		exons = {}

		with open(gtf_file) as gtf:
			for line in gtf:

				# ignore header lines
				if '##' in line:
					continue

				# split each entry
				line = line.strip().split('\t')

				# get some fields from gtf that we care about
				chrom = line[0]
				entry_type = line[2]
				start = int(line[3])
				stop = int(line[4])
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
					if eid not in exons:
						edge = {eid: {'eid': eid,
									  'chrom': chrom,
									  'v1': start,
									  'v2': stop,
									  'strand': strand}}
						exons.update(edge)
			   
			   		# add this exon to the transcript's list of exons
					if tid in transcripts:
						transcripts[tid]['exons'].append(eid)

		# once we have all transcripts, make loc_df
		locs = {}
		vertex_id = 0
		for edge_id, edge in exons.items():
			chrom = edge['chrom']
			strand = edge['strand']

			v1 = edge['v1']
			v2 = edge['v2']

			# exon start
			key = (chrom, v1, strand)
			if key not in locs:
				locs[key] = vertex_id
				vertex_id += 1
			# exon end
			key = (chrom, v2, strand)
			if key not in locs:
				locs[key] = vertex_id
				vertex_id += 1

		# add locs-indexed path to transcripts, and populate edges
		edges = {}
		for _,t in transcripts.items():
			t['path'] = []
			strand = t['strand']
			t_exons = t['exons']

			for i, exon_id in enumerate(t_exons):

				# pull some information from exon dict
				exon = exons[exon_id]
				chrom = exon['chrom']
				v1 = exon['v1']
				v2 = exon['v2']
				strand = exon['strand']

				# add current exon and subsequent intron 
				# (if not the last exon) for each exon to edges
				key = (chrom, v1, v2, strand)
				v1_key = (chrom, v1, strand)
				v2_key = (chrom, v2, strand)
				edge_id = (locs[v1_key], locs[v2_key])
				if key not in edges:
					edges[key] = {'edge_id': edge_id, 'edge_type': 'exon'}

				# add exon locs to path
				t['path'] += list(edge_id)

				# if this isn't the last exon, we also needa add an intron
				# this consists of v2 of the prev exon and v1 of the next exon
				if i < len(t_exons)-1:
					next_exon = exons[t_exons[i+1]]
					v1 = next_exon['v1']
					key = (chrom, v2, v1, strand)
					v1_key = (chrom, v1, strand)
					edge_id = (locs[v2_key], locs[v1_key])
					if key not in edges:
						edges[key] = {'edge_id': edge_id, 'edge_type': 'intron'}

		# turn transcripts, edges, and locs into dataframes
		locs = [{'chrom': key[0],
				 'coord': key[1],
				 'strand': key[2],
				 'vertex_id': vertex_id} for key, vertex_id in locs.items()]
		loc_df = pd.DataFrame(locs)

		edges = [{'v1': item['edge_id'][0],
				  'v2': item['edge_id'][1], 
				  'strand': key[3],
				  'edge_id': item['edge_id'],
				  'edge_type': item['edge_type']} for key, item in edges.items()]
		edge_df = pd.DataFrame(edges)

		transcripts = [{'tid': key,
						'gid': item['gid'],
						'gname': item['gname'],
						'path': item['path']} for key, item in transcripts.items()]
		t_df = pd.DataFrame(transcripts)

		# final df formatting
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

	# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
	def create_dfs_db(self, db, dname):

		# are we pulling a specific dataset's observed transcripts
		# from the database?
		if dname == None:
			print('No dataset name given.'
				' Will pull everything, including unobserved transcripts,'
				' from {}'.format(db))
		else:
			print('Getting transcripts for {} from {}'.format(db, dname))

		# make sure file exists
		if not os.path.exists(db):
			raise Exception('TALON db file {} not found. Check path.'.format(db))

		# open db connection
		conn = sqlite3.connect(db)
		c = conn.cursor()

		# t_df
		t_df = pd.DataFrame()

		# get tid, gid, gname, and paths
		q = """SELECT ga.value, ta.value,
					  t.start_exon, t.jn_path, t.end_exon,
					  t.start_vertex, t.end_vertex
				FROM gene_annotations ga 
				JOIN transcripts t ON ga.ID=t.gene_ID
				JOIN transcript_annotations ta ON t.transcript_ID=ta.ID
			"""
		if dname: 
			q += """ JOIN observed o ON o.transcript_ID=t.transcript_ID"""
		q += """ WHERE ta.attribute='transcript_id'
				AND (ga.attribute='gene_name' 
				OR ga.attribute='gene_id')
			"""
		if dname: 
			q += """ AND o.dataset='{}'""".format(dname)

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

		# edge_df

		# get the list of edge ids we need to pull from db
		edge_ids = list(set([str(n) for path in paths for n in path]))
		edge_str = '({})'.format(','.join(edge_ids))

		q = """SELECT DISTINCT e.* 
				FROM edge e 
				JOIN vertex V ON e.v1=v.vertex_ID 
			"""
		if dname:
			q += """ WHERE e.edge_ID IN {}""".format(edge_str)

		c.execute(q)
		edges = c.fetchall()

		edge_df = pd.DataFrame(edges, 
			columns=['edge_id', 'v1', 'v2',
					 'edge_type', 'strand'])
		edge_df.v1 = edge_df.v1.map(int)
		edge_df.v2 = edge_df.v2.map(int)
		edge_df['talon_edge_id'] = edge_df.edge_id
		edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

		# loc_df

		# get the list of vertex ids we need to pull from db
		edge_ids = edge_df.edge_id.tolist()
		loc_ids = list(set([str(n) for edge_id in edge_ids for n in edge_id]))
		loc_string = '({})'.format(','.join(loc_ids))


		q = """SELECT loc.* FROM location loc"""
		if dname:
			q += """ WHERE loc.location_ID IN {}""".format(loc_string)

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

		# furnish the last bit of info in each df
		t_df['path'] = [[int(n) for n in path]
						 for path in self.get_db_vertex_paths(paths, edge_df)]

		# get strandedness of each sj
		extra_entries = pd.DataFrame(columns=['vertex_id', 'chrom',
			'coord', 'strand', 'plus_vertex_id'])
		loc_df['strand'] = np.nan
		for index, entry in loc_df.iterrows():

			# use v1 or v2 depending on where vertex is in edge
			loc_id = entry.vertex_id
			strands = edge_df.loc[(edge_df.v1==loc_id) | (edge_df.v2==loc_id), 'strand']
			strands = list(set(strands.tolist()))

			# if this sj is found on both the + and - strand,
			# create a minus strand entry to add after the apply
			# function is done. return the plus strand entry
			if len(strands) > 1:
				minus_entry = entry.copy(deep=True)
				minus_entry.strand = '-'
				minus_entry.vertex_id = len(loc_df.index)+len(extra_entries.index)+1
				minus_entry['plus_vertex_id'] = entry.vertex_id
				extra_entries = extra_entries.append(minus_entry)

				loc_df.loc[index, 'strand'] = '+'
			else:
				loc_df.loc[index, 'strand'] = strands[0]
		loc_df = loc_df.append(extra_entries, ignore_index=True, sort=True)

		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')

		# replace sjs in t_df paths that had to be introduced in extra entries
		double_sjs = extra_entries.plus_vertex_id.tolist()
		for tid, entry in t_df.iterrows():
			path = entry.path
			if any(x in path for x in double_sjs):
				coords = loc_df.loc[path[:2], 'coord'].tolist()

				# this transcript is on the minus strand
				if coords[1] < coords[0]:
					new_path = [extra_entries.loc[extra_entries.plus_vertex_id == n,
						'vertex_id'].tolist()[0] if n in double_sjs else n for n in path]
					t_df.loc[tid, 'path'] = new_path

		# also replace sjs in edge_df that had to be introduced in extra_entries
		edge_df.drop('talon_edge_id', axis=1, inplace=True)
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')
		edge_df.to_csv('input_files/edge_df_db.csv')
		for eid, entry in edge_df.iterrows():
			if any(x in eid for x in double_sjs):
				coords = loc_df.loc[eid, 'coord'].tolist()	

				# this edge is on the minus strand
				if coords[1] < coords[0]:
					new_eid = tuple([extra_entries.loc[extra_entries.plus_vertex_id == n, 'vertex_id'].tolist()[0] if n in double_sjs else n for n in eid])
					edge_df.loc[[eid], 'v1'] = new_eid[0]
					edge_df.loc[[eid], 'v2'] = new_eid[1]

		# reset edge ids to use the ones updated from extra_entries
		edge_df.reset_index(drop=True, inplace=True)
		edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')

		# finally, drop unused column from loc_df
		loc_df.drop('plus_vertex_id', axis=1, inplace=True)

		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

	# add node types (internal, TSS, alt TSS, TES, alt_TES) to loc_df
	def get_loc_types(self):

		def label_node_types(locs, vertex_ids, node_type):
			for vertex_id in vertex_ids:
				locs[vertex_id][node_type] = True
			return locs 

		loc_df = self.loc_df
		t_df = self.t_df

		# create a dictionary to hold loc info to speed stuff up
		locs = loc_df.loc[:, ['vertex_id', 'coord']].copy(deep=True)
		locs['internal'] = False
		locs['TSS'] = False
		locs['TES'] = False
		locs['alt_TSS'] = False
		locs['alt_TES'] = False
		locs.drop(['coord', 'vertex_id'], axis=1, inplace=True)
		locs = locs.to_dict('index')

		# label each TSS and TES
		paths = t_df.path.tolist()
		tss = np.unique([path[0] for path in paths])
		locs = label_node_types(locs, tss, 'TSS')
		tes = np.unique([path[-1] for path in paths])
		locs = label_node_types(locs, tes, 'TES')
		internal = np.unique([n for path in paths for n in path[1:-1]])
		locs = label_node_types(locs, internal, 'internal')

		# also create a dictionary of gid: [path1, ... pathn] to speed up
		path_dict = defaultdict(list)
		for tid, entry in t_df.iterrows():
			path_dict[entry.gid].append(entry.path)

		# label alt TES/TSS
		for gid in path_dict.keys():
			paths = path_dict[gid]
			if len(paths) > 1:
				tss = np.unique([path[0] for path in paths])
				tes = np.unique([path[-1] for path in paths])
				if len(tss) > 1:  
					locs = label_node_types(locs, tss, 'alt_TSS')
				if len(tes) > 1:
					locs = label_node_types(locs, tes, 'alt_TES')

		# create df from locs dict and append old loc_df with node types
		locs = pd.DataFrame.from_dict(locs, orient='index')
		loc_df = pd.concat([loc_df, locs], axis=1)
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')

		self.loc_df = loc_df
		self.t_df = t_df

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

	##########################################################################
	######################## Other SwanGraph utilities #####################
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
				self.t_df['tpm_sum'] = self.t_df[tpm_cols].sum(axis=1)
				self.t_df.sort_values(by='tpm_sum', 
									  ascending=False, 
									  inplace=True)
				self.t_df.drop('tpm_sum', axis=1, inplace=True)
			else: 
				raise Exception('Cannot order by expression because '
								'there is no expression data.')

		# order by coordinate of tss in PlottedGraph
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
	######################## Finding "interesting" genes #####################
	##########################################################################

	# returns a list of genes that are "interesting"
	def find_genes_with_novel_isoforms(self):

		# get all the datasets, make sure we're not counting transcripts 
		# that are only in the annotation
		if 'annotation' not in self.datasets:
			raise Exception('No annotation data in graph. Cannot ',
				'determine isoform novelty.')
		datasets = self.get_dataset_cols(include_annotation=False)
		t_df = self.t_df.copy(deep=True)
		t_df = t_df.loc[t_df[datasets].any(axis=1)]

		# how many known and novel isoforms does each gene have
		t_df['known'] = t_df.annotation
		t_df['novel'] = [not i for i in t_df.annotation.tolist()]
		keep_cols = ['annotation', 'known', 'novel', 'gid']
		g_df = t_df[keep_cols].groupby(['gid']).sum()

		# create 'interestingness' column ranking how many novel 
		# compared to known isoforms there are, also ranked by 
		# number of total isoforms
		g_df.known = g_df.known.astype('int32')
		g_df.novel = g_df.novel.astype('int32')
		g_df['interestingness'] = ((g_df.novel+1)/(g_df.known+1))*(g_df.known+g_df.novel)
		g_df.sort_values(by='interestingness', ascending=False, inplace=True)

		# top 10 in case the user doesn't care about whole df
		genes = g_df.index.tolist()[:10]

		return genes, g_df

	# find genes with higher expression in novel than known isoforms
	def find_genes_with_high_novel_expression(self):
		
		# get all the datasets, make sure we're not counting transcripts 
		# that are only in the annotation
		if 'annotation' not in self.datasets:
			raise Exception('No annotation data in graph. Cannot ',
				'determine isoform novelty.')
		datasets = self.get_dataset_cols(include_annotation=False)
		t_df = self.t_df.copy(deep=True)
		t_df = t_df.loc[t_df[datasets].any(axis=1)]

		# how much expression do known and novel isoforms have?
		t_df['known'] = t_df.annotation
		tpm_cols = self.get_tpm_cols()
		keep_cols = tpm_cols+['known', 'gid']
		g_df = t_df[keep_cols].groupby(['gid', 'known']).sum()
		g_df.reset_index(inplace=True)
		g_df['total_known_exp'] = 0
		g_df['total_novel_exp'] = 0
		g_df.loc[g_df.known == True, 'total_known_exp'] = g_df.loc[g_df.known == True, tpm_cols].sum(axis=1) 
		g_df.loc[g_df.known == False, 'total_novel_exp'] = g_df.loc[g_df.known == False, tpm_cols].sum(axis=1) 
		keep_cols = tpm_cols+['total_known_exp', 'total_novel_exp', 'gid']
		g_df = g_df[keep_cols].groupby('gid').sum()

		# create 'interestingness' column ranking how much expression
		# of the gene is attributable to novel isoforms versus known isoforms
		g_df['interestingness'] = ((g_df.total_novel_exp+1)/(g_df.total_known_exp+1))*np.log2(g_df.total_known_exp+1+g_df.total_novel_exp+1)
		g_df.sort_values(by='interestingness', ascending=False, inplace=True)

		# top 10 in case the user doesn't care about whole df
		genes = g_df.index.tolist()[:10]

		return genes, g_df

	# find differentially expressed genes
	# d1 = list or string containing first set of datasets to analyze
	# d2 = list or string containing second set of datasets to analyze
	# the above two will be compared for differential experssion
	def find_differentially_expressed_genes(self, d1, d2):

		# first check to see if the queried expression data is in the graph
		if type(d1) != list: d1 = [d1]
		if type(d2) != list: d2 = [d2]
		datasets = d1+d2
		self.check_abundances(datasets)

		# group t_df into gene df and sum up abundances
		# both across genes and across datasets
		t_df = self.t_df.copy(deep=True)
		d1_cols = self.get_tpm_cols(d1)
		d2_cols = self.get_tpm_cols(d2)
		keep_cols = d1_cols+d2_cols+['gid']
		g_df = t_df[keep_cols].groupby('gid').sum()

		# make sure we're taking the mean expression
		g_df['d1_mean_tpm'] = g_df[d1_cols].sum(axis=1)/len(d1_cols)
		g_df['d2_mean_tpm'] = g_df[d2_cols].sum(axis=1)/len(d2_cols)

		# add pseudocounts for each gene
		g_df.d1_mean_tpm = g_df.d1_mean_tpm + 1
		g_df.d2_mean_tpm = g_df.d2_mean_tpm + 1

		# calculate log2 fold change and order g_df based on magnitude
		# of fold change
		g_df['log_fold_change'] = np.log2(g_df.d1_mean_tpm/g_df.d2_mean_tpm)
		g_df['abs_log_fold_change'] = np.absolute(g_df.log_fold_change)
		g_df.sort_values(by='abs_log_fold_change', ascending=False, inplace=True)

		# top 10 in case the user doesn't care about whole df
		genes = g_df.index.tolist()[:10]

		return genes, g_df

	# find differentially expressed transcripts
	def find_differentially_expressed_transcripts(self, d1, d2):
		
		# first check to see if the queried expression data is in the graph
		if type(d1) != list: d1 = [d1]
		if type(d2) != list: d2 = [d2]
		datasets = d1+d2
		self.check_abundances(datasets)

		# group t_df into gene df and sum up abundances
		# both across genes and across datasets
		t_df = self.t_df.copy(deep=True)
		d1_cols = self.get_tpm_cols(d1)
		d2_cols = self.get_tpm_cols(d2)
		keep_cols = d1_cols+d2_cols+['gid']
		g_df = t_df[keep_cols].groupby('gid').sum()

		# make sure we're taking the mean expression
		g_df['d1_mean_gene_tpm'] = g_df[d1_cols].sum(axis=1)/len(d1_cols)
		g_df['d2_mean_gene_tpm'] = g_df[d2_cols].sum(axis=1)/len(d2_cols)

		# add pseudocounts for each gene
		g_df.d1_mean_gene_tpm = g_df.d1_mean_gene_tpm + 1
		g_df.d2_mean_gene_tpm = g_df.d2_mean_gene_tpm + 1

		# also calculate transcript isoform-level expression
		# and add pseudocounts
		t_df['d1_mean_tpm'] = t_df[d1_cols].sum(axis=1)/len(d1_cols)+1
		t_df['d2_mean_tpm'] = t_df[d2_cols].sum(axis=1)/len(d2_cols)+1

		# merge g_df and t_df to get gene expression and transcript
		# expression in the same dataframe
		keep_cols = ['gid', 'd1_mean_gene_tpm', 'd2_mean_gene_tpm']
		g_df.reset_index(inplace=True)
		t_df = t_df.merge(g_df[keep_cols], on='gid')


		# calculate the isoform fraction as defined by Vitting-Seerup
		# (2017)
		t_df['d1_if'] = t_df['d1_mean_tpm']/t_df['d1_mean_gene_tpm']
		t_df['d2_if'] = t_df['d2_mean_tpm']/t_df['d2_mean_gene_tpm']

		# calculate the log2 fold change of isoform fraction per transcript
		t_df['log_fold_change'] = np.log2(t_df.d1_if/t_df.d2_if)
		t_df['abs_log_fold_change'] = np.absolute(t_df.log_fold_change)

		# groupby to create g_df (again I guess) and sum up over log2
		# if changes
		keep_cols = ['gid', 'abs_log_fold_change']
		g_df = t_df[keep_cols].groupby('gid').sum()
		g_df.reset_index(inplace=True)

		# sort by largest change
		t_df.sort_values('abs_log_fold_change', ascending=False, inplace=True)
		g_df.sort_values('abs_log_fold_change', ascending=False, inplace=True)
		genes = g_df.head(10).gid.tolist()

		return genes, g_df, t_df

	##########################################################################
	######################## Loading/saving SwanGraphs #####################
	##########################################################################

	# saves a splice graph object in pickle format
	def save_graph(self, prefix):
		picklefile = open(prefix+'.p', 'wb')
		pickle.dump(self, picklefile)
		picklefile.close()

	# loads a splice graph object from pickle form
	def load_graph(self, fname):

		picklefile = open(fname, 'rb')
		graph = pickle.load(picklefile)

		# assign SwanGraph fields from file to self
		self.loc_df = graph.loc_df
		self.edge_df = graph.edge_df
		self.t_df = graph.t_df
		self.datasets = graph.datasets
		self.counts = graph.counts
		self.tpm = graph.tpm
		self.pg = graph.pg
		self.G = graph.G

		picklefile.close()

		print('Graph from {} loaded'.format(fname))

	##########################################################################
	############################ Plotting utilities ##########################
	##########################################################################


	# plot the SwanGraph object according to the user's input
	def plot_graph(self, gid, combine=False,
				   indicate_dataset=False,
				   indicate_novel=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel)

		# TODO check to see if we already have a pg object and 
		# reinitialize based on compatibility of object 

		# create PlottedGraph object and plot summary graph
		self.pg = PlottedGraph(self, combine, indicate_dataset, indicate_novel, gid=gid)
		self.pg.plot_graph()

	# plot an input transcript's path through the summary graph 
	def plot_transcript_path(self, tid, combine=False,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)
		self.check_transcript(tid)

		# create PlottedGraph object
		self.pg = PlottedGraph(self,
							   combine,
							   indicate_dataset,
							   indicate_novel,
							   tid=tid,
							   browser=browser)

		self.pg.plot_graph()

	# plots each input transcript path through its gene's summary graph,
	# and automatically saves them
	def plot_each_transcript(self, tids, prefix, combine=False,
						indicate_dataset=False,
						indicate_novel=False,
						browser=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)

		# loop through each transcript in the SwanGraph object
		for tid in tids:
			self.check_transcript(tid)

		print('Plotting transcript for {}'.format(tid))
		print()
		for tid in tids:
			gid = self.get_gid_from_tid(tid)
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
			print('Saving plot for {} as {}'.format(tid, fname))
			self.save_fig(fname)

	# plots each transcript's path through the summary graph, and automatically saves them!
	def plot_each_transcript_in_gene(self, gid, prefix, combine=False,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):

		self.check_plotting_args(combine, indicate_dataset, indicate_novel, browser)

		# loop through each transcript in the SwanGraph object
		tids = self.t_df.loc[self.t_df.gid == gid, 'tid'].tolist()
		print('Plotting {} transcripts for {}'.format(len(tids), gid))
		print()
		for tid in tids:
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
			print('Saving plot for {} as {}'.format(tid, fname))
			self.save_fig(fname)

	# saves current figure named oname. clears the figure space so additional
	# plotting can be done
	def save_fig(self, oname):
		check_file_loc(oname)
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
				   gids,
				   prefix,
				   datasets='all',
				   dataset_groups=False,
				   dataset_group_names=False,
				   heatmap=False,
				   tpm=False,
				   include_unexpressed=False,
				   combine=False,
				   indicate_dataset=False, 
				   indicate_novel=False,
				   browser=False,
				   order='expression'):

		# check to see if input genes are in the graph
		if type(gids) != list:
			gids = [gids]
		for gid in gids:
			self.check_gene(gid)

		# check to see if these plotting settings will play together
		self.check_plotting_args(combine, indicate_dataset,
			indicate_novel, browser)

		# make sure all input datasets are present in graph
		if datasets == 'all':
			datasets = self.get_dataset_cols(include_annotation=False)
		elif not datasets:
			datasets = []
			order = 'tid'
		else:
			self.check_datasets(datasets)

		# if we have dataset groupings make sure that they are a subset
		# of the datasets already requested
		if dataset_groups:
			if not datasets: 
				raise Exception('Cannot group datasets as none were requested.')
			else:
				all_dgs = [j for i in dataset_groups for j in i]
				self.check_datasets(all_dgs)

				subsumed_datasets = [True if i in datasets else False for i in all_dgs]
				if False in subsumed_datasets:
					bad_dataset = all_dgs[subsumed_datasets.index(False)]
					raise Exception("Grouping dataset {} not present in " 
						"datasets {}.".format(bad_dataset, datasets))


		# now check to make sure abundance data is there for the
		# query columns, if user is asking
		if tpm or heatmap:
			self.check_abundances(datasets)
			tpm_cols = self.get_tpm_cols(datasets)
			report_cols = copy.deepcopy(tpm_cols)
		elif datasets:
			report_cols = datasets

		# order transcripts by user's preferences 
		self.order_transcripts(order)

		# loop through each gid and create the report
		for gid in gids:

			# subset based on gid
			t_df = self.t_df.loc[self.t_df.gid == gid].copy(deep=True)

			# if user doesn't care about datasets, just show all transcripts
			if not datasets:
				include_unexpressed = True
			# user only wants transcript isoforms that appear in their data
			if not include_unexpressed:
				counts_cols = self.get_count_cols(datasets)
				report_tids = t_df.loc[t_df[counts_cols].sum(axis=1)>0, 'tid']
			else:
				report_tids = t_df.loc[t_df.gid == gid, 'tid'].tolist()

			# plot each transcript with these settings
			print('Plotting transcripts for {}'.format(gid))
			self.plot_each_transcript(report_tids, prefix, combine,
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

			# if we're grouping things switch up the report_cols 
			# and how t_df is formatted
			if dataset_groups:

				# use 
				if not dataset_group_names:
					print('No group names given. Will just use Group_#.')
					dataset_group_names = ['Group_{}'.format(i) for i in range(len(dataset_groups))]

				# check if we have the right number of group names
				if len(dataset_groups) != len(dataset_group_names):
					print('Not enough group names given. Will just use Group_#.')

				for i in range(len(dataset_groups)):
					group = dataset_groups[i]
					group_name = dataset_group_names[i]

					if not heatmap and not tpm:
						t_df[group_name] = t_df[group].any(axis=1)
					else:
						tpm_group_cols = self.get_tpm_cols(group)
						t_df[group_name] = t_df[tpm_group_cols].mean(axis=1)
				datasets = dataset_group_names
				report_cols = datasets

			# if we're making a heatmap, need to first 
			# add pseudocount, log and normalize tpm values
			if heatmap:

				log_cols = ['{}_log_tpm'.format(d) for d in datasets]
				norm_log_cols = ['{}_norm_log_tpm'.format(d) for d in datasets]
				t_df[log_cols] = np.log2(t_df[tpm_cols]+1)
				max_val = max(t_df[log_cols].max().tolist())
				min_val = min(t_df[log_cols].min().tolist())
				t_df[norm_log_cols] = (t_df[log_cols]-min_val)/(max_val-min_val)
				report_cols = norm_log_cols

				print(t_df[report_cols])

				# create a colorbar 
				plt.figure(1, figsize=(14,1), frameon=False)
				ax = plt.gca()

				cmap = plt.get_cmap('Spectral_r')
				norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
				cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
				cb.set_label('TPM')
				plt.savefig(prefix+'_colorbar_scale.png', format='png', dpi=200)
				plt.clf()
				plt.close()


			# create report
			print('Generating report for {}'.format(gid))
			pdf_name = create_fname(prefix, 
						 combine,
						 indicate_dataset,
						 indicate_novel,
						 browser,
						 ftype='report',
						 gid=gid)
			report = Report(prefix, report_type, report_cols, datasets)
			report.add_page()

			# loop through each transcript and add it to the report
			for tid in report_tids:
				entry = t_df.loc[tid]
				## TODO would be faster if I didn't have to compute these names twice....
				## ie once in plot_each_transcript and once here
				fname = create_fname(prefix,
									 combine,
									 indicate_dataset,
									 indicate_novel, 
									 browser,
									 tid=entry.tid)
				report.add_transcript(entry, fname, heatmap=heatmap)
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
		# the dataset or annotation data exists in the SwanGraph
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





