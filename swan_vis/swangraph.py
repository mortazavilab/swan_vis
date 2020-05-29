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
import anndata
import diffxpy.api as de
from multiprocessing import Pool
from itertools import repeat
from swan_vis.utils import *
from swan_vis.graph import *
from swan_vis.plottedgraph import PlottedGraph
from swan_vis.report import Report

class SwanGraph(Graph):

	def __init__(self, file=None):

		if not file:
			super().__init__()

			# only a SwanGraph should have a plotted graph
			self.pg = PlottedGraph()

			# only a SwanGraph should have DEG and DET data
			self.deg_test = pd.DataFrame()
			self.deg_test_groups = ''
			self.det_test = pd.DataFrame()
			self.det_test_groups = ''

		else:
			check_file_loc(file, 'SwanGraph')
			self.load_graph(file)

	###########################################################################
	############## Related to adding datasets and merging #####################
	###########################################################################

	# add annotation to graph 
	def add_annotation(self, fname):

		# column name for annotation 
		col = 'annotation'

		# use the add_dataset function to add stuff to graph
		self.add_dataset(col, fname, include_isms=True)

		# call all transcripts from the annotation "Known"
		self.t_df.loc[self.t_df.annotation == True, 'novelty'] = 'Known'
		self.t_df.novelty.fillna('Undefined', inplace=True)

	# add dataset into graph from gtf
	def add_dataset(self, col, fname, dname=None,
					counts_file=None, count_cols=None, 
					include_isms=False):

		# make sure that input dataset name is not
		# already in any of the df col spaces
		if col in self.datasets:
			raise Exception('Dataset {} is already in the graph. '
				'Provide a different name.'.format(col))
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

		print('Adding dataset {} to the SwanGraph.'.format(col))

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

		# remove isms if we have access to that information
		if 'novelty' in self.t_df.columns and not include_isms:
			self.t_df = self.t_df.loc[self.t_df.novelty != 'ISM']

		# order node ids by genomic position, add node types,
		# and create graph
		self.update_ids()
		self.order_edge_df()
		self.order_transcripts()
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

		# print note to user about merging with novelty
		existing_cols = self.t_df.columns
		add_cols = b.t_df.columns
		if 'novelty' not in existing_cols and 'novelty' in add_cols:
			print('Novelty info not found for '
			      'existing data. Transcripts '
			      'without novelty information will be '
			      'labelled "Undefined".')
		elif 'novelty' not in add_cols and 'novelty' in existing_cols:
			print('Novelty info not found for '
			     '{} data. Transcripts '
			     'without novelty information will be '
			     'labelled "Undefined".'.format(b_col))

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
		t_df.path = list(t_df.path)

		# assign False to entries that are not in the new dataset, 
		# and to new entries that were not in the prior datasets
		d_cols = self.datasets+[b_col]
		t_df[d_cols] = t_df[d_cols].fillna(value=False, axis=1)

		# deal with novelties
		t_df_cols = t_df.columns.tolist()
		if 'novelty' in t_df_cols or 'novelty_a' in t_df_cols:
			t_df = self.merge_t_df_novelties(t_df)

		# set up index again
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		self.t_df = t_df

	def merge_t_df_novelties(self, t_df):

		# merged dfs with and without novelty
		if 'novelty' in t_df.columns.tolist():
			t_df.fillna(value={'novelty': 'Undefined'},
				inplace=True)

		# merged dfs where both have novelty types
		elif 'novelty_a' in t_df.columns.tolist():

			# if we already have any undefined entries, fill with nan
			t_df.replace({'novelty_a': {'Undefined': np.nan},
						  'novelty_b': {'Undefined': np.nan}}, 
						  inplace=True)

			# first take values that are only present in one dataset
			t_df['novelty_a'].fillna(t_df['novelty_b'], inplace=True)
			t_df['novelty_b'].fillna(t_df['novelty_a'], inplace=True)
			a = t_df[['tid', 'novelty_a']].copy(deep=True)
			a.rename({'novelty_a': 'novelty'}, axis=1, inplace=True)
			a.reset_index(drop=True, inplace=True)
			b = t_df[['tid', 'novelty_b']].copy(deep=True)
			b.rename({'novelty_b': 'novelty'}, axis=1, inplace=True)
			b.reset_index(drop=True, inplace=True)

			# merge novelties on tid and novelty, then extract
			# transcript ids that are duplicated, which represent
			# those that have conflicting novelty assignments
			nov = a.merge(b, on=['tid', 'novelty'], how='outer')
			amb_tids = nov[nov.tid.duplicated()].tid.tolist()

			# label conflicting transcripts as Ambiguous
			if amb_tids:
				print('Novelty types between datasets conflict. Strongly '
					  'consider using input from the same data source to '
					  'reconcile these. Conflicting isoforms will be '
					  'labelled "Ambiguous".')
				nov.set_index('tid', inplace=True)
				nov.loc[amb_tids, 'novelty'] = 'Ambiguous'
				nov.reset_index(inplace=True)
				nov.drop_duplicates(inplace=True)

			# finally, merge new novelty types into t_df
			t_df.drop(['novelty_a', 'novelty_b'], axis=1, inplace=True)
			t_df = t_df.merge(nov, on='tid')

		return t_df


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
		node_types = ['TSS', 'TES', 'internal']

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
		check_file_loc(gtf_file, 'GTF')

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

		# get novelty types associated with each transcript
		def get_transcript_novelties_gtf(fields):
			if fields['transcript_status'] == 'KNOWN':
				return 'Known'
			elif 'ISM_transcript' in fields:
				return 'ISM'
			elif 'NIC_transcript' in fields:
				return 'NIC'
			elif 'NNC_transcript' in fields:
				return 'NNC'
			elif 'antisense_transcript' in fields:
				return 'Antisense'
			elif 'intergenic_transcript' in fields:
				return 'Intergenic'
			elif 'genomic_transcript' in fields:
				return 'Genomic'

		# dictionaries to hold unique edges and transcripts
		transcripts = {}
		exons = {}

		with open(gtf_file) as gtf:
			for line in gtf:

				# ignore header lines
				if line.startswith('#'):
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

					# check if this gtf has transcript novelty vals
					# for the first transcript entry
					if not transcripts:
						if 'talon_transcript' in attributes:
							from_talon = True
						else:
							from_talon = False

					tid = attributes['transcript_id']
					gid = attributes['gene_id']
					gname = attributes['gene_name']

					# add transcript to dictionary 
					entry = {'gid': gid,
							 'gname': gname,
							 'tid': tid,
							 'strand': strand,
							 'exons': []}

					# if we're using a talon gtf, add a novelty field
					if from_talon:
						novelty = get_transcript_novelties_gtf(attributes)
						entry['novelty'] = novelty

					transcript = {tid: entry}
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

		if from_talon:
			transcripts = [{'tid': key,
						'gid': item['gid'],
						'gname': item['gname'],
						'path': item['path'],
						'novelty': item['novelty']} for key, item in transcripts.items()]
		else:
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
				' Will fetch all observed transcripts'
				' from {}'.format(db))
		else:
			print('Getting transcripts for {} from {}'.format(db, dname))

		# make sure file exists
		check_file_loc(db, 'TALON DB')

		# open db connection
		conn = sqlite3.connect(db)
		c = conn.cursor()

		# t_df

		# convert talon novelty values into human-readable values
		def get_transcript_novelties_db(x):
			if x.attribute == 'transcript_status':
				return 'Known'
			elif x.attribute == 'ISM_transcript':
				return 'ISM'
			elif x.attribute == 'NIC_transcript':
				return 'NIC'
			elif x.attribute == 'NNC_transcript':
				return 'NNC'
			elif x.attribute == 'antisense_transcript':
				return 'Antisense'
			elif x.attribute == 'intergenic_transcript':
				return 'Intergenic'
			elif x.attribute == 'genomic_transcript':
				return 'Genomic'

		t_df = pd.DataFrame()

		# get talon tids from observed transcripts in input dataset
		q = """SELECT transcript_ID, gene_ID, start_exon, jn_path, end_exon,
			   start_vertex, end_vertex
			   FROM transcripts t
			   WHERE transcript_ID in
			   (SELECT DISTINCT transcript_ID from observed)"""
		if dname:
			q += """ WHERE dataset='{}'""".format(dname)

		c.execute(q)
		data = c.fetchall()
		talon_tids, talon_gids, paths = zip(*[(i[0], i[1], i[2:]) for i in data])
		paths = self.get_db_edge_paths(paths)

		t_df['talon_tid'] = np.asarray(talon_tids)
		t_df['talon_gid'] = np.asarray(talon_gids)
		t_df['path'] = np.asarray(paths)

		# get transcript ids and merge 
		talon_tids = format_for_in(talon_tids)
		q = """SELECT ID, value 
			   FROM transcript_annotations
			   WHERE attribute='transcript_id'
			   AND ID in {}""".format(talon_tids)
		tids = pd.read_sql_query(q, conn)
		tids.rename({'value': 'tid', 'ID': 'talon_tid'},
			axis=1, inplace=True)
		t_df = t_df.merge(tids, on='talon_tid')

		# get gene ids and merge
		talon_gids = format_for_in(list(set(talon_gids)))
		q = """SELECT ID, value 
			   FROM gene_annotations
			   WHERE attribute='gene_id'
			   AND ID in {}""".format(talon_gids)
		gids = pd.read_sql_query(q, conn)
		gids.rename({'ID': 'talon_gid', 'value': 'gid'},
			axis=1, inplace=True)
		t_df = t_df.merge(gids, on='talon_gid')

		# get gene names and merge
		q = """SELECT ID, value 
			   FROM gene_annotations
			   WHERE attribute='gene_name'
			   AND ID in {}""".format(talon_gids)
		gnames = pd.read_sql_query(q, conn)
		gnames.rename({'ID': 'talon_gid', 'value': 'gname'},
			axis=1, inplace=True)
		t_df = t_df.merge(gnames, on='talon_gid')
		t_df.drop('talon_gid', axis=1, inplace=True)

		# reorder columns
		t_df = t_df[['tid', 'gid', 'gname', 'path', 'talon_tid']]

		# get novelty types and merge
		novelties = ['transcript_status', 'ISM_transcript', 'NIC_transcript', 
					   'NNC_transcript', 'antisense_transcript', 
					   'genomic_transcript', 'intergenic_transcript']
		novelties = format_for_in(novelties)
		q = """SELECT ta.ID, ta.attribute, ta.value 
			   FROM transcript_annotations ta
			   WHERE ta.attribute IN {}
			   AND ta.ID IN {}
			""".format(novelties, talon_tids)
		novelties = pd.read_sql_query(q, conn)
		novelties = novelties[novelties.value != 'NOVEL']
		novelties['novelty'] = novelties.apply(lambda x:
			get_transcript_novelties_db(x), axis=1)
		novelties.drop(['attribute', 'value'], axis=1, inplace=True)

		t_df = t_df.merge(novelties, how='left',
			left_on='talon_tid', right_on='ID')
		t_df.drop(['talon_tid', 'ID'], axis=1, inplace=True)

		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		# edge_df

		# get the list of edge ids we need to pull from db
		edge_ids = list(set([str(n) for path in paths for n in path]))
		# edge_str = '({})'.format(','.join(edge_ids))
		edge_str = format_for_in(edge_ids)

		q = """SELECT DISTINCT e.edge_ID, e.v1,
				e.v2, e.strand, e.edge_type
				FROM edge e 
				JOIN vertex V ON e.v1=v.vertex_ID 
			WHERE e.edge_ID IN {}""".format(edge_str)

		c.execute(q)
		edges = c.fetchall()

		edge_df = pd.DataFrame(edges, 
			columns=['edge_id', 'v1', 'v2',
					 'strand', 'edge_type'])
		edge_df.v1 = edge_df.v1.map(int)
		edge_df.v2 = edge_df.v2.map(int)
		edge_df['talon_edge_id'] = edge_df.edge_id
		edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

		# loc_df

		# get the list of vertex ids we need to pull from db
		edge_ids = edge_df.edge_id.tolist()
		loc_ids = list(set([str(n) for edge_id in edge_ids for n in edge_id]))
		loc_string = format_for_in(loc_ids)

		q = """SELECT loc.* FROM location loc
			   WHERE loc.location_ID IN {}""".format(loc_string)

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
		t_df['path'] = t_df.apply(lambda x:
			self.get_db_vertex_path(x, edge_df), axis=1)

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
					t_df.at[tid, 'path'] = new_path

		# also replace sjs in edge_df that had to be introduced in extra_entries
		edge_df.drop('talon_edge_id', axis=1, inplace=True)
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')
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

	# add node types (internal, TSS, TES) to loc_df
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
	def get_db_vertex_path(self, x, edge_df):
		path = []
		for i, e in enumerate(x.path): 
			entry = edge_df.loc[edge_df.talon_edge_id == e]
			if i == 0:
				path.extend([entry.v1.values[0], entry.v2.values[0]])
			else: path.append(entry.v2.values[0])
		return path

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

	# find genes with isoforms that contain novel intron retention events
	def find_ir_genes(self):

		# get only novel edges
		if 'annotation' not in self.edge_df.columns:
			raise Exception('Cannot find novel IR events without '
				'annotation in SwanGraph.')

		edge_ids = self.edge_df.loc[ \
			(self.edge_df.annotation == False)& \
			(self.edge_df.edge_type == 'exon'), 'edge_id']
		print('Analyzing {} exonic edges for IR'.format(len(edge_ids)))

		# get subset of transcripts that are novel to look for ir edges in
		nt_df = self.t_df.loc[self.t_df.annotation == False]
		
		# for each edge, see if the subgraph between the edge vertices 
		# contains an exonic edge  
		ir_genes = []
		for i, eid in enumerate(edge_ids):
			sub_nodes = [i for i in range(eid[0]+1,eid[1])]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			sub_edges = self.edge_df.loc[sub_edges]
			sub_edges = sub_edges.loc[sub_edges.edge_type == 'intron']

			if len(sub_edges.index) > 0:

				# transcripts that contain the exon-skipping edge
				cand_t_df = nt_df[[eid in vertex_to_edge_path(x) \
					for x in nt_df.path.values.tolist()]]

				# circumvent the ISM bug
				if len(cand_t_df) == 0:
					continue

				# does at least one of the retained introns belong
				# to the same gene as the retaining edge?
				else:
					# genes that contain the intron-retaining edge edge
					cand_genes = cand_t_df.gid.values.tolist()
					cand_g_df = self.t_df.loc[self.t_df.gid.isin(cand_genes)]

					# check if the retained edges are in one of the 
					# intron-retaining genes (wow this is confusing)
					for gid in cand_genes:
						if gid in ir_genes: continue
						for cand_eid in sub_edges.index:
							temp_df = cand_g_df[[cand_eid in vertex_to_edge_path(x) \
									for x in cand_g_df.path.values.tolist()]]
							if len(temp_df.index) > 0:
								ir_genes.append(gid)

		print('Found {} novel ir events from {} genes.'.format(len(ir_genes), 
			len(list(set(ir_genes)))))
		ir_genes = list(set(ir_genes))
		return ir_genes

	# find genes with isoforms that contain novel exon skipping events
	def find_es_genes(self):

		# get only novel edges
		if 'annotation' not in self.edge_df.columns:
			raise Exception('Cannot find novel IR events without '
				'annotation in SwanGraph.')

		edge_ids = self.edge_df.loc[ \
			(self.edge_df.annotation == False)& \
			(self.edge_df.edge_type == 'intron'), 'edge_id']
		print('Analyzing {} intronic edges for ES'.format(len(edge_ids)))

		# get subset of transcripts that are novel to look for ir edges in
		nt_df = self.t_df.loc[self.t_df.annotation == False]

		# for each edge, see if the subgraph between the edge vertices 
		# contains an exonic edge
		es_genes = []
		for eid in edge_ids:
			sub_nodes = [i for i in range(eid[0]+1,eid[1])]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			sub_edges = self.edge_df.loc[sub_edges]
			sub_edges = sub_edges.loc[sub_edges.edge_type == 'exon']

			if len(sub_edges.index) > 0:

				# transcripts that contain the exon-skipping edge
				skip_t_df = nt_df[[eid in vertex_to_edge_path(x) \
					for x in nt_df.path.values.tolist()]]

				# circumvent the ISM bug
				if len(skip_t_df) == 0:
					continue

				# does at least one of the skipped exons belong
				# to the same gene as the skipping edge?
				else:
					# genes that contain the exon-skipping edge
					skip_genes = skip_t_df.gid.values.tolist()
					skip_g_df = self.t_df.loc[self.t_df.gid.isin(skip_genes)]

					# check if the skipped edges are in one of the 
					# exon-skipping genes (wow this is confusing)
					for gid in skip_genes:
						if gid in es_genes: continue
						for skip_eid in sub_edges.index:
							temp_df = skip_g_df[[skip_eid in vertex_to_edge_path(x) \
									for x in skip_g_df.path.values.tolist()]]
							if len(temp_df.index) > 0:
								es_genes.append(gid)

		print('Found {} novel es events from {} genes.'.format(len(es_genes),
			len(list(set(es_genes)))))
		es_genes = list(set(es_genes))
		return es_genes

	# run differential expression test for genes
	def de_gene_test(self, dataset_groups):
		# format expression data to be used by diffxpy
		ann = self.create_gene_anndata(dataset_groups)

		# test
		test = de.test.wald(
		    data=ann,
		    formula_loc="~ 1 + condition",
		    factor_loc_totest="condition")
		test = test.summary()
		test.rename({'gene': 'gid'}, axis=1, inplace=True)


		# add gene name column
		gnames = self.t_df[['gid', 'gname']].copy(deep=True)
		gnames.reset_index(drop=True, inplace=True)
		test = test.merge(gnames, how='left', on='gid')
		test.drop_duplicates(inplace=True)

		# sort on log2fc
		test = test.reindex(test.log2fc.abs().sort_values(ascending=False).index)

		# assign the summary table to the parent object
		self.deg_test = test
		self.deg_test_groups = dataset_groups

		return test

	# subset de genes on user's qval cutoff, and return the top n_genes
	# if n_genes not given, will return all genes called de by qval cutoff
	def get_de_genes(self, q=0.05, n_genes=None):

		# make sure we have the result of a deg test first!
		if self.deg_test.empty:
			raise Exception('Cannot find DE genes without test results. '
				'Run de_gene_test first.')

		# subset on q value 
		test = self.deg_test.loc[self.deg_test.qval <= q].copy(deep=True)

		# list and the df of the top de genes according qval threshold
		if not n_genes:
			genes = test.gname.tolist()
		else:
			if n_genes < len(test.index):
				n_genes = len(test.index)
				test = test.head(n_genes)
				genes = test.gname.tolist()
		return genes, test

	# run differential expression test for transcripts
	def de_transcript_test(self, dataset_groups):
		# format expression data to be used by diffxpy
		ann = self.create_transcript_anndata(dataset_groups)

		# test
		test = de.test.wald(
		    data=ann,
		    formula_loc="~ 1 + condition",
		    factor_loc_totest="condition")
		test = test.summary()
		test.rename({'gene': 'tid'}, axis=1, inplace=True)

		# add gene name column
		gnames = self.t_df[['tid', 'gid', 'gname']].copy(deep=True)
		gnames.reset_index(drop=True, inplace=True)
		test = test.merge(gnames, how='left', on='tid')

		# sort on log2fc
		test = test.reindex(test.log2fc.abs().sort_values(ascending=False).index)

		# assign the summary table to the parent object
		self.det_test = test
		self.det_test_groups = dataset_groups

		return test

	# subset de transcripts on user's qval cutoff, and return the top n_transcripts
	# if n_transcripts not given, will return all transcripts called de by qval cutoff
	def get_de_transcripts(self, q=0.05, n_transcripts=None):

		# make sure we have the result of a deg test first!
		if self.det_test.empty:
			raise Exception('Cannot find DE transcripts without test results. '
				'Run de_transcript_test first.')

		# subset on q value 
		test = self.det_test.loc[self.det_test.qval <= q].copy(deep=True)

		# list and the df of the top de genes according qval threshold
		if not n_transcripts:
			tids = test.gname.tolist()
		else:
			if n_transcripts < len(test.index):
				n_transcripts = len(test.index)
			n_transcripts = test.head(n_transcripts)
			tids = test.transcript.tolist()
		return tids, test

	# find transcript isoforms that are DE where the 
	# parent gene for the isoform is not DE. should be semi-indicative
	# of isoform switching
	def find_isoform_switching_genes(self, q=0.05, n_genes=None):

		# make sure both deg and det tests have been run
		if self.det_test.empty or self.deg_test.empty:
			raise Exception('Cannot find isoform switches without test results. '
				'Run de_gene_test and de_transcript_test first.')

		# subset for genes that aren't DE
		not_degs = self.deg_test.loc[self.deg_test.qval > q]
		not_degs = not_degs.gid

		# subset for dets
		dets = self.det_test.loc[self.det_test.qval <= q]

		# merge on gene id 
		switches = dets.merge(not_degs, how='inner', on='gid')

		# list and the df of the top de genes according qval threshold
		unique_genes = switches.gid.unique().tolist()
		if not n_genes:
			genes = unique_genes
		else:
			if n_genes < len(unique_genes):
				n_genes = len(unique_genes)
			switches = switches.loc[switches.gid.isin(unique_genes[:n_genes])]
			genes = unique_genes[:n_genes]
		return genes, switches

	def get_de_and_not_de_transcripts(self, dataset_groups):
		ann = self.create_transcript_anndata(dataset_groups)
		results = de.test.wald(data=ann,
			formula_loc="~ 1 + condition",
			factor_loc_totest="condition")
		test = results.summary()
		test.rename({'gene': 'transcript'}, axis=1, inplace=True)

		gnames = self.t_df[['tid', 'gname']].copy(deep=True)
		gnames.reset_index(drop=True, inplace=True)
		test = test.merge(gnames, how='left', left_on='transcript', right_on='tid')
		test = test.reindex(test.log2fc.abs().sort_values(ascending=False).index)

		det = test.loc[test.qval < 0.05]
		not_det = test.loc[test.qval >= 0.05]
		genes_w_det = det.gname.tolist()
		not_det = not_det.loc[not_det.gname.isin(genes_w_det)]
		df = pd.concat([det, not_det])
		df = df.loc[df.gname.duplicated(keep=False)]

		return df

	# return an anndata object that can be used to perform different 
	# differential gene expression tests using the diffxpy module
	def create_gene_anndata(self, dataset_groups):

		# group t_df into gene df and sum up abundances
		# both across genes and across datasets
		t_df = self.t_df.copy(deep=True)
		dataset_cols = []
		all_dataset_cols = []
		for group in dataset_groups:
			tpm_cols = self.get_tpm_cols(group)
			dataset_cols.append(tpm_cols)
			all_dataset_cols.extend(tpm_cols)

		keep_cols = all_dataset_cols+['gid']
		g_df = t_df[keep_cols].groupby('gid').sum()

		# add pseudocounts for each gene
		g_df[all_dataset_cols] = g_df[all_dataset_cols] + 1

		# create obs, var, and x entries for the anndata object
		ann_x = g_df.to_numpy().T 
		ann_var = pd.DataFrame(index=g_df.index)
		ann_obs = pd.DataFrame(columns=['batch'],
							   data=all_dataset_cols)
		ann_obs['condition'] = np.nan
		for i, group in enumerate(dataset_cols):
			ann_obs.loc[ann_obs.batch.isin(group),  'condition'] = i
		ann = anndata.AnnData(X=ann_x, var=ann_var, obs=ann_obs)

		return ann

	# returns an anndata object that can be used to perform different 
	# differential transcript expression tests using diffxpy
	def create_transcript_anndata(self, dataset_groups):

		# group t_df into gene df and sum up abundances
		# both across genes and across datasets
		t_df = self.t_df.copy(deep=True)
		dataset_cols = []
		all_dataset_cols = []
		for group in dataset_groups:
			tpm_cols = self.get_tpm_cols(group)
			dataset_cols.append(tpm_cols)
			all_dataset_cols.extend(tpm_cols)

		# add pseudocounts for each transcript
		t_df[all_dataset_cols] = t_df[all_dataset_cols] + 1

		# create obs, var, and x entries for the anndata object
		ann_x = t_df[all_dataset_cols].to_numpy().T 
		ann_var = pd.DataFrame(index=t_df.index)
		ann_obs = pd.DataFrame(columns=['batch'],
							   data=all_dataset_cols)
		ann_obs['condition'] = np.nan
		for i, group in enumerate(dataset_cols):
			ann_obs.loc[ann_obs.batch.isin(group),  'condition'] = i
		ann = anndata.AnnData(X=ann_x, var=ann_var, obs=ann_obs)

		return ann

	##########################################################################
	######################## Loading/saving SwanGraphs #####################
	##########################################################################

	# saves a splice graph object in pickle format
	def save_graph(self, prefix):
		print('Saving graph as '+prefix+'.p')
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

		self.deg_test = graph.deg_test
		self.deg_test_groups = graph.deg_test_groups
		self.det_test = graph.det_test
		self.det_test_groups = graph.det_test_groups

		picklefile.close()

		print('Graph from {} loaded'.format(fname))

	##########################################################################
	############################ Plotting utilities ##########################
	##########################################################################

	# plot a summary SwanGraph for a gene
	def plot_graph(self, gid,
				   indicate_dataset=False,
				   indicate_novel=False,
				   prefix=None):

		if gid not in self.t_df.gid.tolist():
			gid = self.get_gid_from_gname(gid)

		self.check_plotting_args(indicate_dataset, indicate_novel)
		self.check_gene(gid)

		# reinit PlottedGraph object and plot
		self.pg.init_plot_settings(self, gid=gid,
			indicate_dataset=indicate_dataset, 
			indicate_novel=indicate_novel)
		self.pg.plot_graph()

		# if the user has provided a place to save
		if prefix:
			browser = False # can't plot browser for entire gene
			fname = create_fname(prefix,
								indicate_dataset,
								indicate_novel,
								browser,
								ftype='summary',
								gid=gid)
			self.pg.plot_graph()
			print('Saving summary graph for {} as {}'.format(gid, fname))
			save_fig(fname)

	# plot an input transcript's path through the summary graph 
	def plot_transcript_path(self, tid,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False,
							 prefix=None):

		self.check_plotting_args(indicate_dataset, indicate_novel, browser)
		self.check_transcript(tid)

		# reinit PlottedGraph object and plot
		self.pg.init_plot_settings(self, tid=tid, 
			indicate_dataset=indicate_dataset,
			indicate_novel=indicate_novel,
			browser=browser)
		self.pg.plot_graph()

		# if the user has provided a place to save
		if prefix:
			fname = create_fname(prefix,
								indicate_dataset,
								indicate_novel,
								browser,
								ftype='path',
								tid=tid)
			self.pg.plot_graph()
			print('Saving transcript path graph for {} as {}'.format(tid, fname))
			save_fig(fname)


	# plots each input transcript path through its gene's summary graph,
	# and automatically saves them
	def plot_each_transcript(self, tids, prefix,
						indicate_dataset=False,
						indicate_novel=False,
						browser=False):

		self.check_plotting_args(indicate_dataset, indicate_novel, browser)

		# loop through each transcript in the SwanGraph object
		for tid in tids:
			self.check_transcript(tid)

		for tid in tids:
			self.pg.init_plot_settings(self, tid=tid,
				indicate_dataset=indicate_dataset,
				indicate_novel=indicate_novel,
				browser=browser)
			fname = create_fname(prefix,
								 indicate_dataset,
								 indicate_novel,
								 browser,
								 ftype='path',
								 tid=tid)
			self.pg.plot_graph()
			print('Saving transcript path graph for {} as {}'.format(tid, fname))
			save_fig(fname)

	# plots each transcript's path through the summary graph, and automatically saves them!
	def plot_each_transcript_in_gene(self, gid, prefix,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):

		if gid not in self.t_df.gid.tolist():
			gid = self.get_gid_from_gname(gid)
		self.check_gene(gid)

		self.check_plotting_args(indicate_dataset, indicate_novel, browser)

		# loop through each transcript in the SwanGraph object
		tids = self.t_df.loc[self.t_df.gid == gid, 'tid'].tolist()
		print()
		print('Plotting {} transcripts for {}'.format(len(tids), gid))
		for tid in tids:
			self.pg.init_plot_settings(self, tid=tid,
				indicate_dataset=indicate_dataset,
				indicate_novel=indicate_novel,
				browser=browser)
			fname = create_fname(prefix,
								 indicate_dataset,
								 indicate_novel,
								 browser,
								 ftype='path',
								 tid=tid)
			self.pg.plot_graph()
			print('Saving transcript path graph for {} as {}'.format(tid, fname))
			save_fig(fname)

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
				   novelty=False,
				   heatmap=False,
				   tpm=False,
				   include_qvals=False,
				   q=0.05,
				   include_unexpressed=False,
				   indicate_dataset=False, 
				   indicate_novel=False,
				   browser=False,
				   order='expression'):

		# check to see if input genes are in the graph
		if type(gids) != list:
			gids = [gids]
		for i, gid in enumerate(gids):
			if gid not in self.t_df.gid.tolist():
				gid = self.get_gid_from_gname(gid)
				gids[i] = gid
			self.check_gene(gid)

		# check to see if these plotting settings will play together
		self.check_plotting_args(indicate_dataset,
			indicate_novel, browser)

		# make sure all input datasets are present in graph
		if datasets == 'all':
			datasets = self.get_dataset_cols(include_annotation=False)
		elif not datasets:
			datasets = []
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

		# if we've asked for novelty first check to make sure it's there
		if novelty:
			if not self.has_novelty():
				raise Exception('No novelty information present in the graph. '
					'Add it or do not use the "novelty" report option.')

		# check to make sure abundance data is there for the
		# query columns, if user is asking
		if tpm or heatmap:
			self.check_abundances(datasets)

		# order transcripts by user's preferences 
		if order == 'expression' and not self.get_count_cols():
			order = 'tid'
		self.order_transcripts(order)

		# subset t_df based on relevant tids and expression requirements
		t_df = self.t_df[self.t_df.gid.isin(gids)].copy(deep=True)

		# make sure de has been run if needed
		if include_qvals:
			self.check_de('transcript')
			de_df = self.det_test.copy(deep=True)
			t_df = reset_dupe_index(t_df, 'tid')
			t_df['significant'] = False
			t_df = t_df.merge(de_df[['tid', 'qval']], how='left', on='tid')
			t_df['significant'] = t_df.qval <= q
			t_df = set_dupe_index(t_df, 'tid')

		# if user doesn't care about datasets, just show all transcripts
		if not datasets:
			include_unexpressed = True

		# user only wants transcript isoforms that appear in their data
		if not include_unexpressed:
			counts_cols = self.get_count_cols(datasets)
			t_df = t_df[t_df[counts_cols].sum(axis=1)>0]

		# if we're grouping things switch up the datasets 
		# and how t_df is formatted
		if dataset_groups:

			# no grouped dataset names were given - generate names
			if not dataset_group_names:
				print('No group names given. Will just use Group_#.')
				dataset_group_names = ['Group_{}'.format(i) for i in range(len(dataset_groups))]

			# check if we have the right number of group names
			if len(dataset_groups) != len(dataset_group_names):
				print('Not enough group names given. Will just use Group_#.')
				dataset_group_names = ['Group_{}'.format(i) for i in range(len(dataset_groups))]

			for i in range(len(dataset_groups)):
				group = dataset_groups[i]
				group_name = dataset_group_names[i]

				# true or false
				if not heatmap and not tpm:
					t_df[group_name] = t_df[group].any(axis=1)
				# tpm values
				else:
					tpm_group_cols = self.get_tpm_cols(group)
					t_df[group_name] = t_df[tpm_group_cols].mean(axis=1)
			datasets = dataset_group_names
			# report_cols = dataset_group_names

		# determine report type
		if heatmap:
			data_type = 'heatmap'
		elif tpm:
			data_type = 'tpm'
		else:
			data_type = None

		# determine report type
		if not browser:
			report_type = 'swan'
		else:
			report_type = 'browser'

		# parallel
		# launch report jobs on different threads
		with Pool() as pool:
			pool.starmap(create_gene_report, zip(gids, repeat(self), repeat(t_df),
				repeat(datasets), repeat(data_type), repeat(prefix), repeat(indicate_dataset),
				repeat(indicate_novel), repeat(browser), repeat(report_type),
				repeat(novelty), repeat(heatmap), repeat(include_qvals)))

		# # not parallel
		# # loop through each gid and create the report
		# for gid in gids:

		# 	report_tids = t_df.loc[t_df.gid == gid, 'tid'].tolist()

		# 	# plot each transcript with these settings
		# 	print('Plotting transcripts for {}'.format(gid))
		# 	self.plot_each_transcript(report_tids, prefix,
		# 							  indicate_dataset,
		# 							  indicate_novel,
		# 							  browser=browser)

		# 	# if we're plotting tracks, we need a scale as well
		# 	if not browser:
		# 		report_type = 'swan'
		# 	else:
		# 		self.pg.plot_browser_scale()
		# 		self.save_fig(prefix+'_browser_scale.png')
		# 		report_type = 'browser'

		# 	# subset on gene
		# 	gid_t_df = t_df.loc[t_df.gid == gid].copy(deep=True)

		# 	if heatmap:
		# 		# take log2(tpm) and gene-normalize 
		# 		count_cols = ['{}_counts'.format(d) for d in datasets]
		# 		log_cols = ['{}_log_tpm'.format(d) for d in datasets]
		# 		norm_log_cols = ['{}_norm_log_tpm'.format(d) for d in datasets]
		# 		gid_t_df[log_cols] = np.log2(gid_t_df[count_cols]+1)
		# 		max_val = max(gid_t_df[log_cols].max().tolist())
		# 		min_val = min(gid_t_df[log_cols].min().tolist())
		# 		gid_t_df[norm_log_cols] = (gid_t_df[log_cols]-min_val)/(max_val-min_val)

		# 		# create a colorbar 
		# 		plt.rcParams.update({'font.size': 20})
		# 		fig, ax = plt.subplots(figsize=(14, 1.5))
		# 		fig.subplots_adjust(bottom=0.5)
		# 		fig.patch.set_visible(False)
		# 		ax.patch.set_visible(False)

		# 		cmap = plt.get_cmap('Spectral_r')
		# 		norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

		# 		cb = mpl.colorbar.ColorbarBase(ax,
		# 										cmap=cmap,
		# 		                                norm=norm,
		# 		                                orientation='horizontal')
		# 		cb.set_label('log2(TPM)')
		# 		plt.savefig(prefix+'_colorbar_scale.png', format='png', dpi=200)
		# 		plt.clf()
		# 		plt.close()

		# 	# create report
		# 	print('Generating report for {}'.format(gid))
		# 	pdf_name = create_fname(prefix, 
		# 				 indicate_dataset,
		# 				 indicate_novel,
		# 				 browser,
		# 				 ftype='report',
		# 				 gid=gid)
		# 	report = Report(prefix,
		# 					report_type,
		# 					datasets,
		# 					data_type,
		# 					novelty=novelty,
		# 					heatmap=heatmap)
		# 	report.add_page()

		# 	# loop through each transcript and add it to the report
		# 	for tid in report_tids:
		# 		entry = gid_t_df.loc[tid]
		# 		## TODO would be faster if I didn't have to compute these names twice....
		# 		## ie once in plot_each_transcript and once here
		# 		fname = create_fname(prefix,
		# 							 indicate_dataset,
		# 							 indicate_novel, 
		# 							 browser,
		# 							 tid=entry.tid)
		# 		report.add_transcript(entry, fname)
		# 	report.write_pdf(pdf_name)

	##########################################################################
	############################# Error handling #############################
	##########################################################################

	# make sure that the set of arguments work with each other 
	# before we start plotting
	def check_plotting_args(self,
							indicate_dataset,
							indicate_novel,
							browser=False):

		# can only do one or another
		if indicate_dataset and indicate_novel:
			raise Exception('Please choose either indicate_dataset '
							'or indicate_novel, not both.')

		# if indicate_dataset or indicate_novel are chosen, make sure
		# the dataset or annotation data exists in the SwanGraph
		if indicate_novel and 'annotation' not in self.get_dataset_cols():
			raise Exception('Annotation data not present in graph. Use '
							'add_annotation before using indicate_novel')
		if indicate_dataset and indicate_dataset not in self.get_dataset_cols():
			raise Exception('Dataset {} not present in the graph. '
							''.format(indicate_dataset))

		# if browser, can't do indicate_novel, or indicate_dataset
		if browser:
			if indicate_novel or indicate_dataset:
				raise Exception('Cannot indicate_novel or indicate_dataset '
								'with browser option.')

##########################################################################
################################## Extras ################################
##########################################################################

# generate a report for one gene; used for parallelization
def create_gene_report(gid, sg, t_df, 
	datasets, data_type,
	prefix,
	indicate_dataset, indicate_novel,
	browser, 
	report_type, novelty, heatmap, 
	include_qvals):

	report_tids = t_df.loc[t_df.gid == gid, 'tid'].tolist()

	# plot each transcript with these settings
	print()
	print('Plotting transcripts for {}'.format(gid))
	sg.plot_each_transcript(report_tids, prefix,
							  indicate_dataset,
							  indicate_novel,
							  browser=browser)

	# get a different prefix for saving colorbars and scales
	gid_prefix = prefix+'_{}'.format(gid)

	# if we're plotting tracks, we need a scale as well
	if browser:
		sg.pg.plot_browser_scale()
		save_fig(gid_prefix+'_browser_scale.png')

	# subset on gene
	gid_t_df = t_df.loc[t_df.gid == gid].copy(deep=True)

	if heatmap:
		# take log2(tpm) and gene-normalize 
		count_cols = ['{}_counts'.format(d) for d in datasets]
		log_cols = ['{}_log_tpm'.format(d) for d in datasets]
		norm_log_cols = ['{}_norm_log_tpm'.format(d) for d in datasets]
		gid_t_df[log_cols] = np.log2(gid_t_df[count_cols]+1)
		max_val = max(gid_t_df[log_cols].max().tolist())
		min_val = min(gid_t_df[log_cols].min().tolist())
		gid_t_df[norm_log_cols] = (gid_t_df[log_cols]-min_val)/(max_val-min_val)

		# create a colorbar 
		plt.rcParams.update({'font.size': 20})
		fig, ax = plt.subplots(figsize=(14, 1.5))
		fig.subplots_adjust(bottom=0.5)
		fig.patch.set_visible(False)
		ax.patch.set_visible(False)

		cmap = plt.get_cmap('Spectral_r')
		norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

		cb = mpl.colorbar.ColorbarBase(ax,
										cmap=cmap,
		                                norm=norm,
		                                orientation='horizontal')
		cb.set_label('log2(TPM)')
		plt.savefig(gid_prefix+'_colorbar_scale.png', format='png', dpi=200)
		plt.clf()
		plt.close()

	# create report
	print('Generating report for {}'.format(gid))
	pdf_name = create_fname(prefix, 
				 indicate_dataset,
				 indicate_novel,
				 browser,
				 ftype='report',
				 gid=gid)
	report = Report(gid_prefix,
					report_type,
					datasets,
					data_type,
					novelty=novelty,
					heatmap=heatmap,
					include_qvals=include_qvals)
	report.add_page()

	# loop through each transcript and add it to the report
	for tid in report_tids:
		entry = gid_t_df.loc[tid]
		## TODO would be faster if I didn't have to compute these names twice....
		## ie once in plot_each_transcript and once here
		fname = create_fname(prefix,
							 indicate_dataset,
							 indicate_novel, 
							 browser,
							 ftype='path',
							 tid=entry.tid)
		report.add_transcript(entry, fname)
	report.write_pdf(pdf_name)

