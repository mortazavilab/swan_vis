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
from statsmodels.stats.multitest import multipletests
import multiprocessing
from itertools import repeat
# from tqdm import tqdm
from swan_vis.utils import *
from swan_vis.talon_utils import *
from swan_vis.graph import *
from swan_vis.plottedgraph import PlottedGraph
from swan_vis.report import Report

pd.options.mode.chained_assignment = None

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
		adata (anndata AnnData):
			Annotated data holding transcript structure and
			expression information
	"""

	def __init__(self, file=None):

		if not file:
			super().__init__()

			# only a SwanGraph should have a plotted graph
			self.pg = PlottedGraph()

		else:
			check_file_loc(file, 'SwanGraph')
			self.load_graph(file)

	###########################################################################
	############## Related to adding datasets and merging #####################
	###########################################################################

	# def add_datasets(self, config, include_isms=False, verbose=False):
	# 	"""
	# 	Add transcripts from multiple datasets from a config TSV file
	#
	# 	Parameters:
	# 		config (str): Path to TSV config file with the following
	# 			columns (indicated by the header):
	#
	# 			Required:
	# 				config (strs): Path to config file
	#
	# 			Optional:
	# 				dataset_name: Dataset name in TALON db to add transcripts from
	# 					Default=None
	# 				pass_list: TALON pass_list of transcripts to add.
	# 					Default: None
	# 				counts_file: Path to tsv counts matrix
	# 					Default=Non
	# 				count_cols: Column names in counts_file to use
	# 					Default=None
	# 				tid_col: Column name in counts_file containing transcript id
	# 					Default='annot_transcript_id'
	#
	# 		include_isms (bool): Include ISMs from input datasets
	# 			Default=False
	#
	# 		verbose (bool): Display progress
	# 			Default: False
	# 	"""
	#
	# 	# make sure the config file exits
	# 	check_file_loc(config, 'config')
	#
	# 	# read in the config file
	# 	df = pd.read_csv(config, sep='\t')
	#
	# 	# check for required columns
	# 	if 'fname' not in df.columns:
	# 		raise Exception('Please provide the "fname" column in '
	# 			'config file for batch SwanGraph initialization.')
	# 	if 'col' not in df.columns:
	# 		raise Exception('Please provide the "col" column in '
	# 			'config file for batch SwanGraph initialization.')
	#
	# 	# are there any unexpected columns?
	# 	expected_cols = ['fname', 'col',
	# 					 'pass_list', 'dataset_name',
	# 					 'counts_file', 'count_cols',
	# 					 'tid_col', 'include_isms',
	# 					 'verbose']
	# 	for c in df.columns.tolist():
	# 		if c not in expected_cols:
	# 			print('Encountered unexpected column name "{}"'.format(c))
	#
	# 	# loop through each entry in config file
	# 	for ind, entry in df.iterrows():
	#
	# 		# get the values for the rest of the arguments
	# 		file = entry['fname']
	# 		col = entry['col']
	#
	# 		kwargs = {}
	# 		for c in df.columns.tolist():
	# 			if c != 'fname' and c != 'col':
	# 				if not pd.isnull(entry[c]):
	# 					kwargs[c] = entry[c]
	#
	# 		# call add_annotation if we got that as a column
	# 		if col == 'annotation':
	# 			self.add_annotation(file, verbose=verbose)
	#
	# 		# otherwise add_dataset
	# 		else:
	# 			self.add_dataset(col, file,
	# 				include_isms=include_isms,
	# 				verbose=verbose, **kwargs)

	def add_annotation(self, fname, verbose=False):
		"""
		Adds an annotation from input fname to the SwanGraph.

		Parameters:
			fname (str): Path to annotation GTF
			verbose (bool): Display progress
				Default: False
		"""

		# is there already an annotation?
		if self.annotation:
			raise ValueError('Annotation already added to SwanGraph')

		# column name for annotation
		col = 'annotation'

		# use the add_dataset function to add stuff to graph
		self.add_dataset(fname, include_isms=True, \
			annotation=True, verbose=verbose)

		# call all transcripts from the annotation "Known"
		self.t_df.loc[self.t_df.annotation == True, 'novelty'] = 'Known'
		self.t_df.novelty.fillna('Undefined', inplace=True)

		# set flag
		self.annotation = True

	def add_transcriptome(self, fname, pass_list=None, include_isms=False, verbose=False):
		"""
		Adds a whole transcriptome from a set of samples. No abundance is
		included here!

		Parameters:
			fname (str): Path to GTF or TALON db
		"""

		# use the add_dataset function to add transcripts to graph
		self.add_dataset(fname, pass_list=pass_list,
			include_isms=include_isms,
			verbose=verbose)

		# fill NaN annotation transcripts with false
		if 'annotation' in self.t_df.columns:
			self.t_df.annotation.fillna(False, inplace=True)

	def add_dataset(self, fname,
					pass_list=None,
					include_isms=False,
					annotation=False,
					verbose=False):
		"""
		Add transcripts from a dataset from either a GTF or a TALON database.

		Parameters:
			fname (str): Path to GTF or TALON db
			include_isms (bool): Include ISMs from input dataset
				Default: False
			annotation (bool): Whether transcripts being added are from
				an annotation. Set automatically from add_annotation.
				Default: False
			verbose (bool): Display progress
				Default: False
		"""

		# are we dealing with a gtf or a db?
		ftype = gtf_or_db(fname)

		if annotation:
			data = 'annotation'
		else:
			data = 'transcriptome'
		print()
		print('Adding {} to the SwanGraph'.format(data))

		# get loc_df, edge_df, t_df
		if ftype == 'gtf':
			check_file_loc(fname, 'GTF')
			self.create_dfs_gtf(fname, verbose)
		elif ftype == 'db':
			check_file_loc(fname, 'TALON DB')
			self.create_dfs_db(fname, pass_list, verbose)

		# TODO - this won't work if we add the annotation after
		# a transcriptome!
		if annotation:
			self.t_df[data] = True

		# remove isms if we have access to that information
		if 'novelty' in self.t_df.columns and not include_isms:
			self.t_df = self.t_df.loc[self.t_df.novelty != 'ISM']

		# order node ids by genomic position, add node types,
		# and create graph
		if verbose:
			print('Reindexing and sorting entries on genomic location...')
		self.update_ids(verbose=verbose)
		self.order_edge_df()
		self.order_transcripts()
		self.get_loc_types()
		self.create_graph_from_dfs()

		if verbose:
			if annotation:
				data = 'annotation'
			else:
				data = 'transcriptome'
			print()
			print('{} to the SwanGraph'.format(data.capitalize()))


	def add_abundance(self, counts_file):
		"""
		Adds abundance from a counts matrix to the SwanGraph. Transcripts in the
		SwanGraph but not in the counts matrix will be assigned 0 counts.
		Transcripts in the abundance matrix but not in the SwanGraph will not
		have expression added.

		Parameters:
			counts_file (str): Path to TSV expression file where first column is
				the transcript ID and following columns name the added datasets and
				their counts in each dataset, OR to a TALON abundance matrix.
		"""

		# read in abundance file
		check_file_loc(counts_file, 'abundance matrix')
		try:
			df = pd.read_csv(counts_file, sep='\t')
		except:
			raise Error('Problem reading expression matrix {}'.format(counts_file))

		# check if abundance matrix is a talon abundance matrix
		cols = ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_transcript_id',
		 	'annot_gene_name', 'annot_transcript_name', 'n_exons', 'length',
			'gene_novelty', 'transcript_novelty', 'ISM_subtype']
		if df.columns.tolist()[:11] == cols:
			df = reformat_talon_abundance(counts_file)

		# rename transcript ID column
		col = df.columns[0]
		df.rename({col: 'tid'}, axis=1, inplace=True)

		# limit to just the transcripts already in the graph
		sg_tids = self.t_df.tid.tolist()
		ab_tids = df.tid.tolist()

		# right merge to keep everything in the t_df already
		df = df.merge(self.t_df['tid'].to_frame(), how='right', left_on='tid', right_index=True)
		df.drop(['tid_x', 'tid_y'], axis=1, inplace=True)

		# fill NaNs with 0 counts
		df.fillna(0, inplace=True)
		df.set_index('tid', inplace=True)

 		# transpose to get adata format
		df = df.T

		# get adata components - obs, var, and X
		var = df.columns.to_frame()
		var.columns = ['tid']
		obs = df.index.to_frame()
		obs.columns = ['dataset']
		X = df.to_numpy()

		# add each dataset to list of "datasets", check if any are already there!
		datasets = obs.dataset.tolist()
		for d in datasets:
			if d in self.datasets:
				raise ValueError('Dataset {} already present in the SwanGraph.'.format(d))
		self.datasets.extend(datasets)

		print()
		if len(datasets) <= 5:
			print('Adding abundance for datasets {} to SwanGraph.'.format(', '.join(datasets)))
		else:
			mini_datasets = datasets[:5]
			n = len(datasets) - len(mini_datasets)
			print('Adding abundance for datasets {}... (and {} more) to SwanGraph'.format(', '.join(mini_datasets), n))

		# if there is preexisting abundance data in the SwanGraph, concatenate
		# otherwise, adata is the new transcript level adata
		if self.is_empty():

			# create transcript-level adata object
			self.adata = anndata.AnnData(var=var, obs=obs, X=X)

			# add counts as layers
			self.adata.layers['counts'] = self.adata.X
			self.adata.layers['tpm'] = calc_tpm(self.adata, self.t_df).to_numpy()
			self.adata.layers['pi'] = calc_pi(self.adata, self.t_df)[0].to_numpy()
		else:

			# concatenate the raw, tpm, and pi data
			X = np.concatenate((self.adata.layers['counts'], X), axis=0)

			# concatenate obs
			obs = pd.concat([self.adata.obs, obs], axis=0)

			# construct a new adata from the concatenated objects
			adata = anndata.AnnData(var=var, obs=obs, X=X)
			adata.layers['counts'] = X

			# some cleanup for unstructured data
			adata.uns = self.adata.uns

			# set anndata to new guy
			self.adata = adata

			# recalculate pi and tpm
			self.adata.layers['tpm'] = calc_tpm(self.adata, self.t_df).to_numpy()
			self.adata.layers['pi'] = calc_pi(self.adata, self.t_df)[0].to_numpy()


	def add_metadata(self, fname, overwrite=False):
		"""
		Adds metadata to the SwanGraph from a tsv.

		Parameters:
			fname (str): Path / filename of tab-separated input file to add as
				metadata for the datasets in the SwanGraph. Must contain column
				'dataset' which contains dataset names that match those already
				in the SwanGraph.
			overwrite (bool): Whether or not to overwrite duplicated columns
				already present in the SwanGraph.
				Default: False
		"""

		# read in abundance file
		check_file_loc(fname, 'metadata file')
		try:
			df = pd.read_csv(fname, sep='\t')
		except:
			raise Error('Problem reading metadata file {}'.format(fname))

		# which columns are duplicate?
		sg_cols = list(set(self.adata.obs.columns.tolist())-set(['dataset']))
		meta_cols = list(set(df.columns.tolist())-set(['dataset']))
		dupe_cols = [c for c in meta_cols if c in sg_cols]

		# drop columns from one table depending on overwrite settings
		if dupe_cols:
			if overwrite:
				self.adata.obs.drop(dupe_cols, axis=1, inplace=True)
			else:
				df.drop(dupe_cols, axis=1, inplace=True)

		# merge df with adata obs table
		self.adata.obs = self.adata.obs.merge(df, how='left', on='dataset')

		# make sure all dtypes in obs table are non intergers
		for ind, entry in self.adata.obs.dtypes.to_frame().iterrows():
			if entry[0] == 'int64':
				self.adata.obs[ind] = self.adata.obs[ind].astype('str')

		# and set index to dataset
		self.adata.obs['index'] = self.adata.obs.dataset
		self.adata.obs.set_index('index', inplace=True)

	##########################################################################
	############# Related to creating dfs from GTF or TALON DB ###############
	##########################################################################

	def get_current_locs(self):
		"""
		Get a dictionary of unique locations already in the SwanGraph, along
		with the number of unique locations. Used to create the dataframes
		in create_dfs_gtf and create_dfs_db.

		Returns:
			locs (dict): Dictionary of unique genomic locations. Keys are
				(chromosome, coordinate). Items are the vertex IDs from the
				SwanGraph.
			n (int): Number of unique genomic locations. -1 if SwanGraph
				is empty.
		"""
		if len(self.loc_df.index) != 0:
			n = self.loc_df.vertex_id.max()

			ids = self.loc_df.vertex_id.tolist()
			chroms = self.loc_df.chrom.tolist()
			coords = self.loc_df.coord.tolist()
			locs = dict([((ch, co), vid) for ch, co, vid in zip(chroms, coords, ids)])
		else:
			locs = {}
			n = -1

		return locs, n

	def get_current_edges(self):
		"""
		Get a dictionary of unique edges already in the SwanGraph, along with
		the number of unique edges. Used to create the dataframes in
		create_dfs_gtf and create_dfs_db.

		Returns:
			edges (dict): Dictionary of unique edges. Keys are (chromosome,
				start (v1) coord, end (v2) coord, strand, edge_type). Items are
				the edge IDs, v1, v2, and edge_type from the SwanGraph.
			n (int): Number of unique edges. -1 if SwanGraph is emtpy.
		"""
		if len(self.edge_df.index) != 0:
			n = self.edge_df.edge_id.max()

			# populate edges dict with those that already exist
			temp = self.loc_df[['chrom', 'coord', 'vertex_id']]
			temp.reset_index(inplace=True, drop=True)

			self.edge_df = self.edge_df.merge(temp, how='left',
				left_on='v1', right_on='vertex_id')
			self.edge_df.rename({'coord': 'v1_coord'}, axis=1, inplace=True)

			temp = self.loc_df[['coord', 'vertex_id']]
			temp.reset_index(inplace=True, drop=True)

			self.edge_df = self.edge_df.merge(temp, how='left',
				left_on='v2', right_on='vertex_id')
			self.edge_df.rename({'coord': 'v2_coord'}, axis=1, inplace=True)
			drop_cols = ['chrom', 'v1_coord', 'vertex_id_x', 'v2_coord', 'vertex_id_y']

			edges = {}
			for ind, entry in self.edge_df.iterrows():
				ch = entry.chrom
				v1_coord = entry.v1_coord
				v2_coord = entry.v2_coord
				strand = entry.strand
				edge_type = entry.edge_type
				v1 = entry.v1
				v2 = entry.v2
				key = (ch,v1_coord,v2_coord,strand,edge_type)
				item = {'edge_id': entry['edge_id'],
						'edge_type': edge_type,
						'v1': v1,
						'v2': v2}
				edges[key] = item
			self.edge_df.drop(drop_cols, axis=1, inplace=True)
		else:
			edges = {}
			n = -1

		return edges, n

	def create_dfs(self, transcripts, exons,
				   from_talon=False):
		"""
		Create loc, edge, and transcript dataframes. Input parameters are from
		either parse_gtf or parse_db. By default retains structure and metadata
		from preexisting transcripts in the SwanGraph.

		Parameters:
			transcripts (dict of dict): Dictionary of transcripts from file.
				Keys (str): transcript id
				Items (dict): gene id, gene name, transcript id (same as key),
				transcript name, strand, and exons belonging to the transcript.
			exons (dict of dict): Dictionary of exons in GTF.
				Keys (str): exon ids (chromosome_v1_v2_strand_exon)
				Items (dict): edge id (same as key), chromosome, v1, v2, strand,
				and edge type (all exon in this case)
			from_talon (bool): Whether or not the GTF was determined to be
				from TALON.
				Default: False

		Returns:
			loc_df (pandas DataFrame): DataFrame of unique genomic locations
				from input TALON DB or GTF as well as those already in the
				SwanGraph.
			edge_df (pandas DataFrame): DataFrame of unique introns and exons
				from input TALON DB or GTF as well as those already in the
				SwanGraph.
			t_df (pandas DataFrame): DataFrame of unique transcripts from input
				TALON DB or GTF as well as those already in the SwanGraph.
		"""

		# turn each dataframe back into a dict
		if type(transcripts) != dict:
			transcripts = transcripts.to_dict(orient='index')
		if type(exons) != dict:
			exons = exons.to_dict(orient='index')

		locs = self.create_loc_dict(exons)
		transcripts, edges = self.create_transcript_edge_dicts(transcripts, exons, locs)

		# turn transcripts, edges, and locs into dataframes
		locs = [{'chrom': key[0],
				 'coord': key[1],
				 'vertex_id': vertex_id} for key, vertex_id in locs.items()]
		loc_df = pd.DataFrame(locs)

		edges = [{'v1': item['v1'],
				  'v2': item['v2'],
				  'strand': key[3],
				  'edge_id': item['edge_id'],
				  'edge_type': item['edge_type']} for key, item in edges.items()]
		edge_df = pd.DataFrame(edges)

		if from_talon:
			transcripts = [{'tid': key,
						'tname': item['tname'],
						'gid': item['gid'],
						'gname': item['gname'],
						'path': item['path'],
						'novelty': item['novelty']} for key, item in transcripts.items()]
		else:
			transcripts = [{'tid': key,
						'tname': item['tname'],
						'gid': item['gid'],
						'gname': item['gname'],
						'path': item['path']} for key, item in transcripts.items()]
		t_df = pd.DataFrame(transcripts)

		# remove entries that are already in the SwanGraph so we just get
		# dfs of new entries
		tids = self.t_df.tid.tolist()
		eids = self.edge_df.edge_id.tolist()
		vids = self.loc_df.vertex_id.tolist()
		t_df = t_df.loc[~t_df.tid.isin(tids)]
		edge_df = edge_df.loc[~edge_df.edge_id.isin(eids)]
		loc_df = loc_df.loc[~loc_df.vertex_id.isin(vids)]

		# format indices of dfs
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')
		edge_df = create_dupe_index(edge_df, 'edge_id')
		edge_df = set_dupe_index(edge_df, 'edge_id')
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		# and concatenate with those that exist
		loc_df = pd.concat([self.loc_df, loc_df])
		edge_df = pd.concat([self.edge_df, edge_df])
		t_df = pd.concat([self.t_df, t_df])

		# account for mixed novelty content
		t_df['novelty'] = t_df['novelty'].fillna('Undefined')

		return loc_df, edge_df, t_df

	def create_loc_dict(self, exons):
		"""
		Create location dictionary using the exons found from a GTF or TALON
		DB.

		Parameters:
			exons (dict): Dictionary of exons found in the GTF or TALON DB.
				Keys (str): chr_v1_c2_strand_exon.
				Items (dict): eid, chrom, v1, v2, strand, edge_type

		Returns:
			locs (dict): Dictionary of locations and their vertex IDs
				in the SwanGraph, including existing and novel locations.
				Keys: (chrom, coord). Items: vertex_id.
		"""

		locs, vertex_id = self.get_current_locs()
		vertex_id += 1
		for eid, edge in exons.items():
			chrom = edge['chrom']
			v1 = edge['v1']
			v2 = edge['v2']

			# exon start
			key = (chrom, v1)
			if key not in locs:
				locs[key] = vertex_id
				vertex_id += 1

			# exon end
			key = (chrom, v2)
			if key not in locs:
				locs[key] = vertex_id
				vertex_id += 1

		return locs

	def create_transcript_edge_dicts(self, transcripts, exons, locs):
		"""
		Create edge dictionary using the exons found from a GTF or TALON DB, and
		add the edge path to the transcript dictionary (only dictionary
		attribute added by this function).

		Parameters:
			transcripts (dict): Dictionary of transcripts found in the GTF or
				TALON DB.
				Keys (str): transcript ID
				Items (dict): gid, gname, tid, tname, strand, novelty (if talon),
					exons (list of 'chr_coord1_coord2_strand_exon').
			exons (dict): Dictionary of exons found in the GTF or TALON DB.
				Keys (str): chr_v1_c2_strand_exon.
				Items (dict): eid, chrom, v1, v2, strand, edge_type

		Returns:
			transcripts (dict): Dictionary of only NEW transcripts.
				Keys (str): transcript ID
			edges (dict): Dictionary of edges and their edge IDs
				in the SwanGraph, including existing and novel edges.
				Keys: TODO
				Items: TODO



		"""

		edges, edge_id = self.get_current_edges()
		edge_id += 1
		for key, t in transcripts.items():
			t['path'] = []
			strand = t['strand']
			t_exons = t['exons']

			# reorder exons that are in weird orders from the GTF
			t_exons = reorder_exons(t_exons)

			for i, exon_id in enumerate(t_exons):

				# pull some information from exon dict
				exon = exons[exon_id]
				chrom = exon['chrom']
				v1 = exon['v1']
				v2 = exon['v2']
				strand = exon['strand']

				# add current exon and subsequent intron
				# (if not the last exon) for each exon to edge
				# also add this exon to the current transcript's path
				key = (chrom, v1, v2, strand, 'exon')
				v1_key = (chrom, v1)
				v2_key = (chrom, v2)
				if key not in edges:
					edges[key] = {'edge_id': edge_id, 'edge_type': 'exon',
								  'v1': locs[v1_key], 'v2': locs[v2_key]}
					t['path'] += [edge_id]
					edge_id += 1
				else:
					t['path'] += [edges[key]['edge_id']]

				# if this isn't the last exon, we also need to add an intron
				# this consists of v2 of the prev exon and v1 of the next exon
				# also add this intron to the current transcript's path
				if i < len(t_exons)-1:
					next_exon = exons[t_exons[i+1]]
					v1 = next_exon['v1']
					key = (chrom, v2, v1, strand, 'intron')
					v1_key = (chrom, v1)
					if key not in edges:
						edges[key] = {'edge_id': edge_id, 'edge_type': 'intron',
									  'v1': locs[v2_key], 'v2': locs[v1_key]}
						t['path'] += [edge_id]
						edge_id += 1
					else:
						t['path'] += [edges[key]['edge_id']]
		return transcripts, edges

	# create loc_df (nodes), edge_df (edges), and t_df (transcripts) from gtf
	# adapted from Dana Wyman and TALON
	def create_dfs_gtf(self, gtf_file, verbose, annotation):
		"""
		Create pandas DataFrames for unique genomic locations, exons, introns,
		and transcripts from an input GTF file. Merge the resultant DataFrames
		with those already in the SwanGraph. Called from add_dataset.

		Parameters:
			gtf_file (str): Path to GTF file
			verbose (bool): Display progress
			annotation (bool): GTF contains annotation
		"""

		t_df, exon_df, from_talon = parse_gtf(gtf_file, verbose)

		# keep track of transcripts from GTF if we're adding annotation
		if annotation:
			annot_tids = t_df.tid.tolist()

		loc_df, edge_df, t_df = self.create_dfs(t_df, exon_df, from_talon)

		# add annotation column
		if annotation:
			t_df['annotation'] = t_df.loc[t_df.tid.isin(annot_tids)]

		# set SwanGraph dfs
		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

		# # concatenate dfs
		# self.loc_df = pd.concat([self.loc_df, loc_df])
		# self.edge_df = pd.concat([self.edge_df, edge_df])
		# self.t_df = pd.concat([self.t_df, t_df])

	# create SwanGraph dataframes from a TALON db. Code very ripped from
	# TALON's create_GTF utility
	def create_dfs_db(self, database, pass_list, verbose):
		"""
		Create pandas DataFrames for unique genomic locations, exons, introns,
		and transcripts from an input TALON database. Merge the resultant DataFrames
		with those already in the SwanGraph. Called from add_dataset.

		Parameters:
			database (str): Path to database file
			pass_list (str): Path to TALON pass list files
			vebose (bool): Display progress
		"""

		# make sure files exist
		if pass_list:
			check_file_loc(pass_list, 'pass list')

		# annot = check_annot_validity(annot, database)

		pass_list = handle_filtering(database, True, pass_list)

		# create separate gene and transcript pass_lists
		gene_pass_list = []
		transcript_pass_list = []
		for key,group in itertools.groupby(pass_list,operator.itemgetter(0)):
			gene_pass_list.append(key)
			for id_tuple in list(group):
				transcript_pass_list.append(id_tuple[1])

		# get gene, transcript, and exon annotations
		gene_annotations = get_annotations(database, "gene",
										   pass_list = gene_pass_list)
		transcript_annotations = get_annotations(database, "transcript",
												 pass_list = transcript_pass_list)
		exon_annotations = get_annotations(database, "exon")

		# get transcript data from the database
		gene_2_transcripts = get_gene_2_transcripts(database,
							 transcript_pass_list)

		# get exon location info from database
		exon_ID_2_location = fetch_exon_locations(database)

		transcripts = {}
		exons = {}

		if verbose:
			n_transcripts = len(transcript_pass_list)
			pbar = tqdm(total=n_transcripts)
			pbar.set_description('Processing transcripts')

		# loop through genes, transcripts, and exons
		for gene_ID, transcript_tuples in gene_2_transcripts.items():
			curr_annot = gene_annotations[gene_ID]
			gene_annotation_dict = {}
			for annot in curr_annot:
				attribute = annot[3]
				value = annot[4]
				gene_annotation_dict[attribute] = value

			# check if there's a gene name field and add one if not
			if 'gene_name' not in gene_annotation_dict:
				gene_annotation_dict['gene_name'] = gene_annotation_dict['gene_id']

			# get transcript entries
			for transcript_entry in transcript_tuples:
				transcript_ID = transcript_entry["transcript_ID"]

				# if this transcript is already in the SwanGraph, skip
				# to the next one
				if transcript_ID in self.t_df.tid.tolist():
					continue

				curr_transcript_annot = transcript_annotations[transcript_ID]

				transcript_annotation_dict = {}
				for annot in curr_transcript_annot:
					attribute = annot[3]
					value = annot[4]
					transcript_annotation_dict[attribute] = value


				if 'transcript_name' not in transcript_annotation_dict:
					transcript_annotation_dict['transcript_name'] = transcript_annotation_dict['transcript_id']
				tid = transcript_annotation_dict['transcript_id']
				tname = transcript_annotation_dict['transcript_name']
				gid = gene_annotation_dict['gene_id']
				gname = gene_annotation_dict['gene_name']
				strand = transcript_entry['strand']
				novelty = get_transcript_novelties(transcript_annotation_dict)

				# add transcript to dictionary
				entry = {'gid': gid,
						 'gname': gname,
						 'tid': tid,
						 'tname': tname,
						 'strand': strand,
						 'novelty': novelty,
						 'exons': []}
				transcript = {tid: entry}
				transcripts.update(transcript)

				if verbose:
					pbar.update(1)

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

		# TODO - pull this shit out and stick it in its own function cause it's
		# the same as what's in create_dfs_gtf
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
					'tname': item['tname'],
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

	def get_loc_types(self):
		"""
		Determine what role (TSS, TES, internal) each unique genomic location
		plays in the transcripts it is used in. Assigns each location in loc_df
		a boolean label for the TSS, TES, and internal columns based on its use.
		"""

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


	def order_transcripts(self, order='tid'):
		"""
		Order the transcripts in the SwanGraph t_df based on an input heuristic.
		Can order alphabetically by transcript ID, by expression of each
		transcript, or by the genomic location of the transcripts' TSSs or TESs.

		Parameters:
			order (str): Method to order transcripts by. Choose from ['tid',
				'expression', 'tss', 'tes']
		"""

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

		# order by coordinate of tss
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
	######################## Analysis tools  #################################
	##########################################################################

	def find_ir_genes(self):
		"""
		Finds all unique genes containing novel intron retention events.
		Requires that an annotation has been added to the SwanGraph.

		Returns:
			ir_genes (list of str): A list of gene ids from the SwanGraph with
				at least one novel intron retention event
			ir_transcripts (list of str): A list of transcript ids from the
				SwanGraph with at least one novel intron retention event
			ir_edges (list of tuples): A list of exonic edges in the
				SwanGraph that retain at least one intronic edge
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
		ir_edges = []
		ir_genes = []
		ir_transcripts = []
		for i, eid in enumerate(edge_ids):
			sub_nodes = [i for i in range(eid[0]+1,eid[1])]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			sub_edges = self.edge_df.loc[sub_edges]
			# find edges that are intronic; if there are none, this is not
			# an intron-retaining edge
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
							tids = cand_t_df.tid.tolist()
							if len(temp_df.index) > 0:
								ir_edges.append(eid)
								ir_genes.append(gid)
								ir_transcripts.extend(tids)

		ir_genes = list(set(ir_genes))
		ir_transcripts = list(set(ir_transcripts))
		ir_edges = list(set(ir_edges))

		print('Found {} novel ir events from {} transcripts.'.format(len(ir_genes),
			len(ir_transcripts)))

		return ir_genes, ir_transcripts, ir_edges

	def find_es_genes(self):
		"""
		Finds all unique genes containing novel exon skipping events.
		Requires that an annotation has been added to the SwanGraph.

		Returns:
			es_genes (list of str): A list of gene ids from the SwanGraph with
				at least one novel exon skipping event
			es_transcripts (list of str): A list of transcript ids from the
				SwanGraph with at least one novel exon skipping event
			es_edges (list of tuples): A list of intronic edges in the
				SwanGraph that skip at least one exonic edge
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
		es_edges = []
		es_genes = []
		es_transcripts = []
		for eid in edge_ids:
			# subgraph consisting of all nodes between the candidate
			# exon-skipping edge coords in order and its edges
			sub_nodes = [i for i in range(eid[0]+1,eid[1])]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			sub_edges = self.edge_df.loc[sub_edges]
			# find edges that are exonic; if there are none, this is not
			# an exon-skipping edge
			sub_edges = sub_edges.loc[sub_edges.edge_type == 'exon']

			if len(sub_edges.index) > 0:

				# transcripts that contain the candidate exon-skipping edge
				skip_t_df = nt_df[[eid in vertex_to_edge_path(x) \
					for x in nt_df.path.values.tolist()]]

				# circumvent the ISM bug
				if len(skip_t_df) == 0:
					continue

				# does at least one of the skipped exons belong
				# to the same gene as the skipping edge?
				else:
					# genes that contain the candidate exon-skipping edge
					skip_genes = skip_t_df.gid.values.tolist()
					skip_g_df = self.t_df.loc[self.t_df.gid.isin(skip_genes)]

					# check if the skipped edges are in one of the
					# exon-skipping genes (wow this is confusing)
					for gid in skip_genes:
						if gid in es_genes: continue
						for skip_eid in sub_edges.index:
							# transcripts with the exons that are skipped
							temp_df = skip_g_df[[skip_eid in vertex_to_edge_path(x) \
									for x in skip_g_df.path.values.tolist()]]
							tids = skip_t_df.tid.tolist()
							if len(temp_df.index) > 0:
								es_edges.append(eid)
								es_genes.append(gid)
								es_transcripts.extend(tids)

		es_genes = list(set(es_genes))
		es_transcripts = list(set(es_transcripts))
		es_edges = list(set(es_edges))

		print('Found {} novel es events in {} transcripts.'.format(len(es_edges),
			len(es_transcripts)))


		return es_genes, es_transcripts, es_edges

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

	# filter differential isoform expression test results based on
	# both an adjusted p value and dpi cutoff
	# TODO make this more object-oriented when migration to adata happens
	# ie store test results in adata.uns? Indexed by obs_col, obs_conditions
	def filter_die_results(self, df, p=0.05, dpi=10):
		"""
		Filters differential isoform expression test results based on adj.
		p-value and change in percent isoform usage (dpi).

		Parameters:
			df (pandas DataFrame): DIE test results, output from
				get_die_genes
			p (float): Adj. p-value threshold to declare a gene as isoform
				switching / having DIE.
				Default: 0.05
			dpi (float): DPI (in percent) value to threshold genes with DIE
				Default: 10

		Returns:
			gids (list of str): List of gene ids that pass the
				significance thresholds
			test (pandas DataFrame): Summary table of genes that pass
				the significance threshold
		"""
		df = df.loc[(df.adj_p_val<=p)&(df.dpi>=dpi)]
		gids = df['index'].tolist()
		return gids, df

	def get_die_genes(self, obs_col='dataset',
					  obs_conditions=None, rc_thresh=10):
		"""
		Finds genes with differential isoform expression between two conditions
		that are in the obs table. If there are more than 2 unique values in
		`obs_col`,

		Parameters:
			obs_col (str): Column name from self.adata.obs table to group on.
				Default: 'dataset'
			obs_conditions (list of str, len 2): Which conditions from obs_col
				to compare? Required if obs_col has more than 2 unqiue values.
			rc_thresh (int): Number of reads required for each conditions
				in order to test the gene.
				Default: 10

		Returns:
			test (pandas DataFrame): A summary table of the differential
				isoform expression test, including p-values and adjusted
				p-values, as well as change in percent isoform usage (dpi) for
				all tested genes.
		"""

		# check if there are more than 2 unique values in obs_col
		if len(self.adata.obs[obs_col].unique().tolist())!=2 and not obs_conditions:
			raise ValueError('Must provide obs_conditions argument when obs_col has'\
				' >2 unique values')
		elif not obs_conditions:
			obs_conditions = self.adata.obs[obs_col].unique().tolist()
		elif len(obs_conditions) != 2:
			raise ValueError('obs_conditions must have exactly 2 values')

		# use calculated values already in the SwanGraph
		if obs_col == 'dataset':
			df = pd.DataFrame(data=self.adata.layers['pi'],
								 index=self.adata.obs.index,
								 columns=self.adata.var.index)
			sums = calc_total_counts(self.adata, obs_col)

		# recalculate pi and aggregate counts
		else:
			df, sums = calc_pi(self.adata, self.t_df, obs_col=obs_col)

		df = df.transpose()
		sums = sums.transpose()

		# limit to obs_conditions
		if obs_conditions:
			if len(obs_conditions) != 2:
				raise Error('obs_conditions must have exactly two values.')
			df = df[obs_conditions]
			sums = sums[obs_conditions]
		else:
			obs_conditions = adata.obs[obs_col].unique().tolist()

		# add gene id
		df = df.merge(self.t_df['gid'], how='left', left_index=True, right_index=True)

		# add counts and cumulative counts across the samples
		count_cols = sums.columns
		sums['total_counts'] = sums[count_cols].sum(axis=1)
		df = df.merge(sums, how='left', left_index=True, right_index=True, suffixes=(None, '_counts'))

		# construct tables for each gene and test!
		gids = df.gid.unique().tolist()
		gene_de_df = pd.DataFrame(index=gids, columns=['p_val', 'dpi'], data=[[np.nan for i in range(2)] for j in range(len(gids))])
		gene_de_df.index.name = 'gid'
		for gene in gids:
			gene_df = df.loc[df.gid==gene]
			gene_df = get_die_gene_table(gene_df, obs_conditions, rc_thresh)
			# if the gene is valid for testing, do so
			if isinstance(gene_df, pd.DataFrame):
				p, dpi = test_gene(gene_df, obs_conditions)
			else:
				p = dpi = np.nan
			gene_de_df.loc[gene, 'p_val'] = p
			gene_de_df.loc[gene, 'dpi'] = dpi

		# remove untestable genes and perform p value
		# Benjamini-Hochberg correction
		gene_de_df.dropna(axis=0, inplace=True)
		p_vals = gene_de_df.p_val.tolist()
		_, adj_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')
		gene_de_df['adj_p_val'] = adj_p_vals
		gene_de_df.reset_index(inplace=True)

		return gene_de_df

# 	def create_gene_anndata(self, dataset_groups):
# 		"""
# 		Creates a gene-level AnnData object containing TPM that's
# 		compatible with diffxpy. Assigns different condition labels
# 		to the given dataset groups.
#
# 		Parameters:
# 			dataset_groups (list of list of str, len 2): Grouping of datasets
# 				from the SwanGraph to be used in the differential
# 				expression test
# 				Example: [['data1','data2'],['data3','data4']]
#
# 		Returns:
# 			ann (AnnData): AnnData object containing gene-level TPM
# 				with different conditions labelled for DE testing
# 		"""
#
# 		# group t_df into gene df and sum up abundances
# 		# both across genes and across datasets
# 		t_df = self.t_df.copy(deep=True)
# 		dataset_cols = []
# 		all_dataset_cols = []
# 		for group in dataset_groups:
# 			tpm_cols = self.get_tpm_cols(group)
# 			dataset_cols.append(tpm_cols)
# 			all_dataset_cols.extend(tpm_cols)
#
# 		keep_cols = all_dataset_cols+['gid']
# 		g_df = t_df[keep_cols].groupby('gid').sum()
#
# 		# add pseudocounts for each gene
# 		g_df[all_dataset_cols] = g_df[all_dataset_cols] + 1
#
# 		# create obs, var, and x entries for the anndata object
# 		ann_x = g_df.to_numpy().T
# 		ann_var = pd.DataFrame(index=g_df.index)
# 		ann_obs = pd.DataFrame(columns=['dataset'],
# 							   data=all_dataset_cols)
# 		ann_obs['condition'] = np.nan
# 		for i, group in enumerate(dataset_cols):
# 			ann_obs.loc[ann_obs.dataset.isin(group),  'condition'] = i
# 		ann = anndata.AnnData(X=ann_x, var=ann_var, obs=ann_obs)
#
# 		return ann
#
# 	def create_transcript_anndata(self, obs_col='dataset',
# 								  obs_conditions=None, how='counts'):
# 		"""
# 		Creates a transcript-level AnnData object containing TPM that's
# 		compatible with diffxpy. Assigns different condition labels
# 		to the given dataset groups.
#
# Parameters:
# 	obs_col (str): Column name from self.adata.obs table to group on.
# 		Default: 'dataset'
# 	obs_conditions (list of str, len 2): Which conditions from obs_col
# 		to compare? Required if obs_col has more than 2 unqiue values.
# 	rc_thresh (int): Number of reads required for each conditions
# 		in order to test the gene.
# 		Default: 10
#
# 		Returns:
# 			ann (AnnData): AnnData object containing transcript-level TPM
# 				with different conditions labelled for DE testing
# 		"""
#
# 		# group t_df
# 		t_df = self.t_df.copy(deep=True)
# 		dataset_cols = []
# 		all_dataset_cols = []
# 		for group in dataset_groups:
# 			if how == 'tpm':
# 				cols = self.get_tpm_cols(group)
# 			elif how == 'counts':
# 				cols = self.get_count_cols(group)
# 			dataset_cols.append(cols)
# 			all_dataset_cols.extend(cols)
#
# 		if how == 'tpm':
# 			# add pseudocounts for each transcript
# 			t_df[all_dataset_cols] = t_df[all_dataset_cols] + 1
#
# 		# create obs, var, and x entries for the anndata object
# 		ann_x = t_df[all_dataset_cols].to_numpy().T
# 		ann_var = t_df[['gid', 'gname']]
# 		# ann_var['tid'] = ann_var.index
# 		# ann_var = pd.DataFrame(index=t_df.index)
# 		ann_obs = pd.DataFrame(columns=['dataset'],
# 							   data=all_dataset_cols)
# 		ann_obs['condition'] = np.nan
# 		for i, group in enumerate(dataset_cols):
# 			ann_obs.loc[ann_obs.dataset.isin(group),  'condition'] = i
# 		ann = anndata.AnnData(X=ann_x, var=ann_var, obs=ann_obs)
#
# 		return ann

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
				   prefix=None,
				   display=False):
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
			display (bool): Display the plot during runtime
				Default: False
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

		# display the plot if option is given
		if display:
			plt.show()

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
							 prefix=None,
							 display=False):
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
			display (bool): Display the plot during runtime
				Default: False
		"""

		self.check_plotting_args(indicate_dataset, indicate_novel, browser)
		self.check_transcript(tid)

		# reinit PlottedGraph object and plot
		self.pg.init_plot_settings(self, tid=tid,
			indicate_dataset=indicate_dataset,
			indicate_novel=indicate_novel,
			browser=browser)
		self.pg.plot_graph()

		# display the plot if option is given
		if display:
			plt.show()

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
				   dpi=False, # TODO this option is named stupidly. change it
				   cmap='Spectral_r',
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
			dpi (bool): Plot proportion isoform usage per condition
				as opposed to log2(tpm)
			cmap (str): Matplotlib color map to display heatmap values
				in.
				Default: 'Spectral_r'
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
		if dpi or tpm or heatmap:
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
					report_type, novelty, heatmap, dpi, cmap, include_qvals)
		else:
			# launch report jobs on different threads
			with multiprocessing.Pool(threads) as pool:
				pool.starmap(_create_gene_report, zip(gids, repeat(self), repeat(t_df),
					repeat(datasets), repeat(data_type), repeat(prefix), repeat(indicate_dataset),
					repeat(indicate_novel), repeat(browser), repeat(report_type),
					repeat(novelty), repeat(heatmap), repeat(dpi),
					repeat(cmap), repeat(include_qvals)))

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
	dpi, cmap,
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

		if not dpi:
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

			try:
				cmap = plt.get_cmap(cmap)
			except:
				raise ValueError('Colormap {} not found'.format(cmap))

			norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

			cb = mpl.colorbar.ColorbarBase(ax,
											cmap=cmap,
											norm=norm,
											orientation='horizontal')
			cb.set_label('log2(TPM)')
			plt.savefig(gid_prefix+'_colorbar_scale.png', format='png', dpi=200)
			plt.clf()
			plt.close()

		# calculate % isoform for each dataset group
		elif dpi:
			count_cols = ['{}_counts'.format(d) for d in datasets]
			dpi_cols = ['{}_dpi'.format(d) for d in datasets]
			gid_t_df[dpi_cols] = gid_t_df[count_cols].div(gid_t_df[count_cols].sum(axis=0), axis=1)

			# create a colorbar
			plt.rcParams.update({'font.size': 20})
			fig, ax = plt.subplots(figsize=(14, 1.5))
			fig.subplots_adjust(bottom=0.5)
			fig.patch.set_visible(False)
			ax.patch.set_visible(False)

			try:
				cmap = plt.get_cmap(cmap)
			except:
				raise ValueError('Colormap {} not found'.format(cmap))

			norm = mpl.colors.Normalize(vmin=0, vmax=1)

			cb = mpl.colorbar.ColorbarBase(ax,
											cmap=cmap,
											norm=norm,
											orientation='horizontal')
			cb.set_label('Proportion of isoform use')
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
					dpi=dpi,
					cmap=cmap,
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
