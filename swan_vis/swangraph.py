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
import multiprocessing
from itertools import repeat
from swan_vis.utils import *
from swan_vis.talon_utils import *
from swan_vis.graph import *
from swan_vis.plottedgraph import PlottedGraph
from swan_vis.report import Report

class SwanGraph(Graph):
	""" 
	A graph class to represent a transcriptome and perform
	plotting and analysis from it

		Attributes:

			datasets (list of str):
				Names of datasets in the Graph
			counts (list of str):
				Names of columns holding counts in the Graph
			tpm (list of str):
				Names of columns holding tpm values in the Graph
			loc_df (pandas DataFrame): 
				DataFrame of all unique observed genomic 
				coordinates in the transcriptome
			edge_df (pandas DataFrame):
				DataFrame of all unique observed exonic or intronic
				combinations of splice sites in the transcriptome
			t_df (pandas DataFrame): 
				DataFrame of all unique transcripts found 
				in the transcriptome
			pg (swan PlottedGraph):
				The PlottedGraph holds the information from the most 
				recently made plot
			deg_test (pandas DataFrame): 
				A summary table of the results of a differential gene
				expression test
			deg_test_groups (list of str, len 2):
				The configuration of groupings used to run the differential
				gene expression test
				det_test (pandas DataFrame): 
				A summary table of the results of a differential transcript
				expression test
			det_test_groups (list of str, len 2):
				The configuration of groupings used to run the differential
				transcript expression test
	"""

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

	def add_annotation(self, fname):
		"""
		Adds an annotation from input fname to the SwanGraph.

			Parameters:
				fname (str): Path to annotation GTF
		"""

		# column name for annotation 
		col = 'annotation'

		# use the add_dataset function to add stuff to graph
		self.add_dataset(col, fname, include_isms=True)

		# call all transcripts from the annotation "Known"
		self.t_df.loc[self.t_df.annotation == True, 'novelty'] = 'Known'
		self.t_df.novelty.fillna('Undefined', inplace=True)

	def add_dataset(self, col, fname,
					dataset_name=None,
					whitelist=None,
					annot=None,
					counts_file=None, count_cols=None, 
					tid_col='annot_transcript_id',
					include_isms=False):
		"""
		Add transcripts from a dataset from either a GTF or a TALON database.

			Parameters:

				col (str): Name of column to add data to in the SwanGraph
				fname (str): Path to GTF or TALON db

				# Only for loading from TALON
				dataset_name (str): Dataset name in TALON db to add transcripts from
					Default=None
				whitelist (str): TALON whitelist of transcripts to add.
					Default: None
				annot (str): TALON annotation name in database to 
					add transcripts from
					Default: None

				# Only if also adding abundance
				counts_file (str): Path to tsv counts matrix
					Default=None
				count_cols (str or list of str): Column names in counts_file to use
					Default=None
				tid_col (str): Column name in counts_file containing transcript id
					Default='annot_transcript_id'

				include_isms (bool): Include ISMs from input dataset
					Default=False
		"""

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
				self.create_dfs_db(fname, annot, whitelist, dataset_name)

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
				temp.create_dfs_db(fname, annot, whitelist, dataset_name)
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
			self.add_abundance(counts_file, count_cols, col, tid_col)
 
	def add_abundance(self, counts_file, count_cols,
					  dataset_name, tid_col='annot_transcript_id'):
		"""
		Adds abundance information to an existing dataset in the SwanGraph.

			Parameters:

				counts_file (str): Path to tsv counts matrix
				count_cols (str or list of str): Column names in counts_file to use
				dataset_name (str): Name of SwanGraph dataset to associate abundance with
				tid_col (str): Column name in counts_file containing transcript id
					Default='annot_transcript_id'
		"""

		# if the dataset we're trying to add counts too doesn't exist
		if dataset_name not in self.datasets:
			raise Exception('Trying to add expression data to a dataset '
							'that is not in the graph. Add dataset to graph first.')

		# get counts from input abundance file 
		abundance_df = process_abundance_file(counts_file, count_cols, tid_col)
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
		self.t_df.path = self.t_df.path.map(tuple)
		b.t_df.path = b.t_df.path.map(tuple)

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

		id_map = list(zip(self.loc_df.vertex_id_b,
						  self.loc_df.vertex_id_a))
		id_map = [list(i) for i in id_map]

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
						novelty = get_transcript_novelties(attributes)
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

		# create inverse loc dict to sort paths by
		locs_inv = {v: k for k, v in locs.items()}

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

			# sort the path based on chromosomal coordinates and strand
			# in case there's some weird ordering in the gtf 
			t['path'] = reorder_locs(t['path'], strand, locs_inv)

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

	# create SwanGraph dataframes from a TALON db. Code very ripped from 
	# TALON's create_GTF utility
	def create_dfs_db(self, database, annot, whitelist, dataset):

		# make sure file exists
		check_file_loc(database, 'TALON DB')

		annot = check_annot_validity(annot, database)

		whitelist = handle_filtering(database, 
											annot, 
											True, 
											whitelist, 
											dataset)
		# create separate gene and transcript whitelists
		gene_whitelist = []
		transcript_whitelist = []
		for key,group in itertools.groupby(whitelist,operator.itemgetter(0)):
			gene_whitelist.append(key)
			for id_tuple in list(group):
				transcript_whitelist.append(id_tuple[1])

		# get gene, transcript, and exon annotations
		gene_annotations = get_annotations(database, "gene", annot, 
										   whitelist = gene_whitelist)  
		transcript_annotations = get_annotations(database, "transcript", annot,
												 whitelist = transcript_whitelist) 
		exon_annotations = get_annotations(database, "exon", annot)

		# get transcript data from the database
		gene_2_transcripts = get_gene_2_transcripts(database, 
							 transcript_whitelist)

		# get exon location info from database
		exon_ID_2_location = fetch_exon_locations(database)

		transcripts = {}
		exons = {}

		# loop through genes, transcripts, and exons
		for gene_ID, transcript_tuples in gene_2_transcripts.items():
			curr_annot = gene_annotations[gene_ID]
			gene_annotation_dict = {}
			for annot in curr_annot:
				attribute = annot[3]
				value = annot[4]
				gene_annotation_dict[attribute] = value
		
			# get transcript entries
			for transcript_entry in transcript_tuples:
				transcript_ID = transcript_entry["transcript_ID"]
				curr_transcript_annot = transcript_annotations[transcript_ID]

				transcript_annotation_dict = {}
				for annot in curr_transcript_annot:
					attribute = annot[3]
					value = annot[4]
					transcript_annotation_dict[attribute] = value
		  
				tid = transcript_annotation_dict['transcript_id']
				gid = gene_annotation_dict['gene_id']  
				gname = gene_annotation_dict['gene_name']
				strand = transcript_entry['strand'] 
				novelty = get_transcript_novelties(transcript_annotation_dict)  

				# add transcript to dictionary 
				entry = {'gid': gid,
						 'gname': gname,
						 'tid': tid,
						 'strand': strand,
						 'novelty': novelty,
						 'exons': []}
				transcript = {tid: entry}
				transcripts.update(transcript)
						 
				if transcript_entry["n_exons"] != 1:
					transcript_edges = [str(transcript_entry["start_exon"])] + \
									   str(transcript_entry["jn_path"]).split(",")+ \
									   [str(transcript_entry["end_exon"])]
				else:
					transcript_edges = [transcript_entry["start_exon"]]

				# get exon entries
				for exon_ID in transcript_edges[::2]:
					exon_ID = int(exon_ID)
					curr_exon_annot = exon_annotations[exon_ID]

					exon_annotation_dict = {}
					for annot in curr_exon_annot:
						attribute = annot[3]
						value = annot[4]
						exon_annotation_dict[attribute] = value

					e_tuple = exon_ID_2_location[exon_ID]
					chrom = e_tuple[0]
					start = e_tuple[1]
					stop = e_tuple[2]
					strand = e_tuple[3]
					start, stop = find_edge_start_stop(start, stop, strand)
					eid = '{}_{}_{}_{}_exon'.format(chrom, start, stop, strand)

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
		# print(dict(list(transcripts.items())[:3]))
		for _,t in transcripts.items():
			t['path'] = []
			strand = t['strand']
			t_exons = t['exons']

			for i, exon_id in enumerate(t_exons):
				# print('shouldnt u be in here')
				# exit()

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
					'path': item['path'],
					'novelty': item['novelty']} for key, item in transcripts.items()]

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

	# add node types (internal, TSS, TES) to loc_df
	def get_loc_types(self):

		self.loc_df['internal'] = False
		self.loc_df['TSS'] = False
		self.loc_df['TES'] = False

		# get lists of locations that are used as TSS, TES
		paths = self.t_df.path.tolist()
		internal = list(set([n for path in paths for n in path[1:-1]]))
		tss = [path[0] for path in paths]
		tes = [path[-1] for path in paths]

		# set node types in t_df
		self.loc_df.loc[internal, 'internal'] = True
		self.loc_df.loc[tss, 'TSS'] = True
		self.loc_df.loc[tes, 'TES'] = True

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

	def find_ir_genes(self):
		"""
		Finds all unique genes containing novel intron retention events. 
		Requires that an annotation has been added to the SwanGraph.

			Returns:

				ir_genes (list of str): A list of gene ids from the SwanGraph with 
					at least one novel intron retention event
				ir_transcripts (list of str): A list of transcript ids from the 
					SwanGraph with at least one novel intron retention event
		"""

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
		ir_transcripts = []
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
							tids = temp_df.tid.tolist()
							if len(temp_df.index) > 0:
								ir_genes.append(gid)
								ir_transcripts.extend(tids)

		print('Found {} novel ir events from {} genes.'.format(len(ir_genes), 
			len(list(set(ir_genes)))))
		ir_genes = list(set(ir_genes))
		ir_transcripts = list(set(ir_transcripts))

		return ir_genes, ir_transcripts

	def find_es_genes(self):
		"""
		Finds all unique genes containing novel exon skipping events. 
		Requires that an annotation has been added to the SwanGraph.

			Returns:

				es_genes (list of str): A list of gene ids from the SwanGraph with 
					at least one novel exon skipping event
				es_transcripts (list of str): A list of transcript ids from the 
					SwanGraph with at least one novel exon skipping event
		"""

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
		es_transcripts = []
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
							tids = temp_df.tid.tolist()
							if len(temp_df.index) > 0:
								es_genes.append(gid)
								es_transcripts.extend(tids)

		print('Found {} novel es events from {} genes.'.format(len(es_genes),
			len(list(set(es_genes)))))
		es_genes = list(set(es_genes))
		es_transcripts = list(set(es_transcripts))

		return es_genes, es_transcripts

	def de_gene_test(self, dataset_groups):
		""" 
		Runs a differential expression test on the gene level.

			Parameters:

				dataset_groups (list of list of str, len 2): Grouping of datasets 
					from the SwanGraph to be used in the differential
					expression test
					Example: [['data1','data2'],['data3','data4']]

			Returns: 

				test (pandas DataFrame): A summary table of the differential
					expression test, including p and q-values, as well 
					as log fold change.
		"""

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

	def get_de_genes(self, q=0.05, n_genes=None):
		"""
		Subsets the differential gene expression test summary table based
		on a q-value cutoff. Requires that de_gene_test has already been
		run.

			Parameters:

				q (float): q-value threshold to declare a gene as significant
					Default: 0.05
				n_genes (int): Number of results to return. 
					Default: None (returns all found significant)

			Returns:

				genes (list of str): List of gene names that pass the 
					significance threshold
				test (pandas DataFrame): Summary table of genes that pass the
					significance threshold
		"""

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

	def de_transcript_test(self, dataset_groups):
		""" 
		Runs a differential expression test on the transcript level.

			Parameters:

				dataset_groups (list of list of str, len 2): Grouping of datasets 
					from the SwanGraph to be used in the differential
					expression test
					Example: [['data1','data2'],['data3','data4']]

			Returns: 

				test (pandas DataFrame): A summary table of the differential
					expression test, including p and q-values, as well 
					as log fold change.
		"""

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

	def get_de_transcripts(self, q=0.05, n_transcripts=None):
		"""
		Subsets the differential transcript expression test summary table based
		on a q-value cutoff. Requires that de_transcript_test has already been
		run.

			Parameters:

				q (float): q-value threshold to declare a transcript as significant
					Default: 0.05
				n_transcripts (int): Number of results to return. 
					Default: None (returns all found significant)

			Returns:

				tids (list of str): List of transcript ids that pass the 
					significance threshold
				test (pandas DataFrame): Summary table of transcripts that pass
					the significance threshold
		"""

		# make sure we have the result of a deg test first!
		if self.det_test.empty:
			raise Exception('Cannot find DE transcripts without test results. '
				'Run de_transcript_test first.')

		# subset on q value 
		test = self.det_test.loc[self.det_test.qval <= q].copy(deep=True)

		# list and the df of the top de genes according qval threshold
		if not n_transcripts:
			tids = test.tid.tolist()
		else:
			if n_transcripts < len(test.index):
				n_transcripts = len(test.index)
			n_transcripts = test.head(n_transcripts)
			tids = test.transcript.tolist()
		return tids, test

	def find_isoform_switching_genes(self, q=0.05, n_genes=None):
		""" Finds isoform switching genes; genes that are not differentially
			expressed but contain at least one transcript that is. Requires
			that de_gene_test and de_transcript_test have been run.

			Parameters:

				q (float): q-value threshold to declare a gene/transcript 
					as significant
					Default: 0.05
				n_genes (int): Number of results to return. 
					Default: None (returns all found significant)
			
			Returns:

				genes (list of str): List of gene names that are categorized as 
					isoform switching
				switches (pandas DataFrame): Summary table of genes that are 
					categorized as isoform switching
		"""

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

	def save_graph(self, prefix):
		"""
		Saves the current SwanGraph in pickle format with the .p extension

			Parameters: 

				prefix (str): Path and filename prefix. Resulting file will 
					be saved as prefix.p
		"""
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

	def plot_graph(self, gid,
				   indicate_dataset=False,
				   indicate_novel=False,
				   prefix=None):
		"""
		Plots a gene summary SwanGraph for an input gene.
		Does not automatically save the figure by default!

			Parameters:

				gid (str): Gene ID to plot for (can also be gene name but 
					we've seen non-unique gene names so use at your own risk!)
				indicate_dataset (str): Dataset name from SwanGraph to
					highlight with outlined nodes and dashed edges
					Incompatible with indicate_novel
					Default: False (no highlighting)
				indicate_novel (bool): Highlight novel nodes and edges by 
					outlining them and dashing them respectively
					Incompatible with indicate_dataset
					Default: False
				prefix (str): Path and file prefix to automatically save
					the plotted figure
					Default: None, won't automatically save 
		"""

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

	def plot_transcript_path(self, tid,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False,
							 prefix=None):
		"""
		Plots a path of a single transcript isoform through a gene summary 
		SwanGraph.

			Parameters:

				tid (str): Transcript id of transcript to plot
				indicate_dataset (str): Dataset name from SwanGraph to
					highlight with outlined nodes and dashed edges
					Incompatible with indicate_novel
					Default: False (no highlighting)
				indicate_novel (bool): Highlight novel nodes and edges by 
					outlining them and dashing them respectively
					Incompatible with indicate_dataset
					Default: False
				browser (bool): Plot transcript models in genome browser-
					style format. Incompatible with indicate_dataset and
					indicate_novel
				prefix (str): Path and file prefix to automatically save
					the plotted figure
					Default: None, won't automatically save
		"""

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

	def plot_each_transcript(self, tids, prefix,
						indicate_dataset=False,
						indicate_novel=False,
						browser=False):
		"""
		Plot each input transcript and automatically save figures

			Parameters:

				tids (list of str): List of transcript ids to plot
				prefix (str): Path and file prefix to automatically save
					the plotted figures
				indicate_dataset (str): Dataset name from SwanGraph to
					highlight with outlined nodes and dashed edges
					Incompatible with indicate_novel
					Default: False (no highlighting)
				indicate_novel (bool): Highlight novel nodes and edges by 
					outlining them and dashing them respectively
					Incompatible with indicate_dataset
					Default: False
				browser (bool): Plot transcript models in genome browser-
					style format. Incompatible with indicate_dataset and
					indicate_novel
		"""

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

	def plot_each_transcript_in_gene(self, gid, prefix,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False):
		"""
		Plot each transcript in a given gene and automatically save figures

			Parameters:

				gid (str): Gene id or gene name to plot transcripts from
				prefix (str): Path and file prefix to automatically save
					the plotted figures
				indicate_dataset (str): Dataset name from SwanGraph to
					highlight with outlined nodes and dashed edges
					Incompatible with indicate_novel
					Default: False (no highlighting)
				indicate_novel (bool): Highlight novel nodes and edges by 
					outlining them and dashing them respectively
					Incompatible with indicate_dataset
					Default: False
				browser (bool): Plot transcript models in genome browser-
					style format. Incompatible with indicate_dataset and
					indicate_novel
		"""

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
				   order='expression',
				   threads=1):
		"""
		Generates a PDF report for a given gene or list of genes according
		to the user's input.

			Parameters: 

				gids (str or list of str): Gene ids or names to generate
					reports for
				prefix (str): Path and/or filename prefix to save PDF and
					images used to generate the PDF

				datasets (list of str): Datasets to include in the report
					Default: Include columns for all datasets
				dataset_groups (list of list of str): Datasets to average
					together in the report and display as one column
					Example: [['group1_1','group1_2'],['group2_1','group2_2']]
				dataset_group_names (list of str): Names to give each group 
					given by dataset_groups. Must be the same length as 
					dataset_groups
					Example: ['group1', 'group2']
					Default: Will assign numbers 1 through length(dataset_group)

				novelty (bool): Include a column to dipslay novelty type of
					each transcript. Requires that a TALON GTF or DB has 
					been used to load data in
					Default: False

				heatmap (bool): Display expression values in a heatmap
					format. Requires that abundance information has been 
					added to the SwanGraph
					Default: False
				tpm (bool): Display TPM value of each transcript/dataset 
					combination, instead of presence/absence of each 
					transcript. Requires that abundance information has
					been added to the SwanGraph
					Default:False

				include_qvals (bool): Display q-val of differential expression 
					for each transcript and bold entries found to be
					differentially expressed. Requires that de_transcript_test
					has been run, and that abundance information has been
					added to the SwanGraph
					Default: False
				q (float): Q-value significance threshold to use when 
					bolding transcripts if include_qvals = True.
					Default: 0.05

				include_unexpressed (bool): Add transcript entries to report
					that are not expressed in any input dataset.
					Default: False

				indicate_dataset (str): Dataset name from SwanGraph to
					highlight with outlined nodes and dashed edges
					Incompatible with indicate_novel
					Default: False (no highlighting)
				indicate_novel (bool): Highlight novel nodes and edges by 
					outlining them and dashing them respectively
					Incompatible with indicate_dataset
					Default: False
				browser (bool): Plot transcript models in genome browser-
					style format. Incompatible with indicate_dataset and
					indicate_novel

				order (str): Order to display transcripts in the report.
					Options are 
						'tid': alphabetically by transcript ID
						'expression': cumulative expression from high to low
							Requires that abundance information has been 
							added to the SwanGraph
						'tss': genomic coordinate of transcription start site
						'tes': genomic coordinate of transcription end site
					Default: 'expression' if abundance information is present,
							 'tid' if not

				threads (int): Number of threads to use. Multithreading is 
					recommended when making a large number of gene reports.
					Default: 1. 
		"""

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
			if not self.check_de('transcript'):
				raise Exception('Differential transcript expression test needed '
					'to use include_qvals. Run de_transcript_test.')
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
					group_name += '_counts'
					count_group_cols = self.get_count_cols(group)
					t_df[group_name] = t_df[count_group_cols].mean(axis=1)
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

		# make sure number of threads is compatible with the system
		max_cores = multiprocessing.cpu_count()
		if threads > max_cores:
			threads = max_cores

		# if there are fewer than 10 genes, the overhead for multiprocessing 
		# takes longer 
		if len(gids) < 10:
			for gid in gids:
				_create_gene_report(gid, self, t_df, datasets, data_type,
					prefix, indicate_dataset, indicate_novel, browser,
					report_type, novelty, heatmap, include_qvals) 
		else:
			# launch report jobs on different threads
			with multiprocessing.Pool(threads) as pool:
				pool.starmap(_create_gene_report, zip(gids, repeat(self), repeat(t_df),
					repeat(datasets), repeat(data_type), repeat(prefix), repeat(indicate_dataset),
					repeat(indicate_novel), repeat(browser), repeat(report_type),
					repeat(novelty), repeat(heatmap), repeat(include_qvals)))

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
def _create_gene_report(gid, sg, t_df, 
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

