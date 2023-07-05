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
from tqdm import tqdm
import matplotlib.colors as mc
import colorsys
from scipy import sparse
import scanpy as sc
import pdb

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

	Attributes
	----------
	datasets (list of str):
		Names of datasets in the Graph
	gene_datasets (list of str):
		Names of datasets w/ gene expression in the Graph
	annotation (bool):
		Whether an annotation transcriptome has been added.
	abundance (bool):
		Whether abundance data has been added to the SwanGraph
	gene_abundance (bool):
		Whether gene abundance data has been added to the SwanGraph
	loc_df (pandas DataFrame):
		DataFrame of all unique observed genomic
		coordinates in the transcriptome
	edge_df (pandas DataFrame):
		DataFrame of all unique observed exonic or intronic
		combinations of splice sites in the transcriptome
	t_df (pandas DataFrame):
		DataFrame of all unique transcripts found
		in the transcriptome
	adata (anndata AnnData):
		Annotated data object to hold transcript expression values
		and metadata
	gene_adata (anndata AnnData):
		Annotated data object to hold gene expression values / metadata
	ic_adata (anndata AnnData):
		Annotated data object to hold intron chain expression values / metadata
	tss_adata (anndata AnnData):
		Annotated data object to hold TSS expression values and metadata
	tes_adata (anndata AnnData):
		Annotated data object to hold TES expression values and metadata
	"""

	def __init__(self,
				 sc=False,
				 edge_adata=True,
				 end_adata=True,
				 ic_adata=True):
		"""
		Parameters:
			sc (bool): Whether this is coming from single cell data
			edge_adata (bool): Whether to create the edge_adata table
			end_adata (bool): Whether to create the tss/tes_adata tables
			ic_adata (bool): Whether to create the ic_adata table
		"""

		super().__init__()

		# only a SwanGraph should have a plotted graph
		self.pg = PlottedGraph()

		if sc:
			self.sc = True
		else:
			self.sc = False

		if edge_adata:
			self.make_edge_adata = True
		else:
			self.make_edge_adata = False

		if end_adata:
			self.make_end_adata = True
		else:
			self.make_end_adata = False

		if ic_adata:
			self.make_ic_adata = True
		else:
			self.make_ic_adata = False

	###########################################################################
	############## Related to adding datasets and merging #####################
	###########################################################################

	def add_annotation(self,
					   fname,
					   verbose=False):
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

		# use the add_dataset function to add stuff to graph
		self.add_dataset(fname, include_isms=True, \
			annotation=True, verbose=verbose)

		# call all transcripts from the annotation "Known"
		self.t_df.loc[self.t_df.annotation == True, 'novelty'] = 'Known'
		self.t_df.novelty.fillna('Undefined', inplace=True)

		# set flag
		self.annotation = True

	def add_transcriptome(self,
						  fname,
						  pass_list=None,
						  include_isms=False,
						  verbose=False):

		"""
		Adds a whole transcriptome from a set of samples.

		Parameters:
			fname (str): Path to GTF or TALON db
			pass_list (str): Path to pass list file (if passing a TALON DB)
			include_isms (bool): Include ISMs from input dataset
				Default: False
			verbose (bool): Display progress
				Default: False
		"""

		# use the add_dataset function to add transcripts to graph
		tids = self.add_dataset(fname, pass_list=pass_list,
			include_isms=include_isms,
			verbose=verbose)

		# fill NaN annotation transcripts with false
		if 'annotation' in self.t_df.columns:
			self.t_df.annotation.fillna(False, inplace=True)

	def add_dataset(self,
					fname,
					pass_list=None,
					include_isms=False,
					annotation=False,
					verbose=False):
		"""
		Add transcripts from a dataset from either a GTF or a TALON database.

		Parameters:
			fname (str): Path to GTF or TALON db
			pass_list (str): Path to pass list file of transcript IDs to include
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
			t_df, exon_df, from_talon, from_cerberus = parse_gtf(fname, include_isms, verbose)
		elif ftype == 'db':
			check_file_loc(fname, 'TALON DB')
			observed = True
			t_df, exon_df = parse_db(fname, pass_list, observed,
									 include_isms, verbose)
			from_talon = True
			from_cerberus = False

		# keep track of transcripts from GTF if we're adding annotation
		if annotation:
			annot_tids = t_df.tid.tolist()

		# create the dfs with new and preexisting data and assign them to the
		# df fields of the SwanGraph
		loc_df, edge_df, t_df = self.create_dfs(t_df, exon_df, from_talon, from_cerberus)

		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df

		# add location path
		self.get_loc_path()

		# remove ISM transcripts and locations and edges exclusively from ISM
		# transcripts
		if not include_isms:
			self.remove_isms()

		# label elements from the annotation
		if annotation:
			self.label_annotated(annot_tids)

		# # remove isms if we have access to that information
		# if 'novelty' in self.t_df.columns and not include_isms:
		# 	self.t_df = self.t_df.loc[self.t_df.novelty != 'ISM']

		# if there's already an annotation, label the transcripts added
		# that are not in the annotation
		if self.annotation and not annotation:
			self.t_df.annotation = self.t_df.annotation.fillna(False)
			self.edge_df.annotation = self.edge_df.annotation.fillna(False)
			self.loc_df.annotation = self.loc_df.annotation.fillna(False)

		# order node ids by genomic position, add node types, add loc_path,
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
			print('{} added to the SwanGraph'.format(data.capitalize()))

	def remove_isms(self):
		"""
		Remove ISM transcripts as well as nodes and edges that exclusively
		come from ISM transcripts.
		"""
		if 'novelty' in self.t_df.columns:

			# transcripts
			self.t_df = self.t_df.loc[self.t_df.novelty != 'ISM']

			# edges
			edge_df = pivot_path_list(self.t_df, 'path')
			eids = edge_df.edge_id.unique().tolist()
			self.edge_df = self.edge_df.loc[eids]

			# locs
			loc_df = pivot_path_list(self.t_df, 'loc_path')
			lids = loc_df.vertex_id.unique().tolist()
			self.loc_df = self.loc_df.loc[lids]

	def abundance_to_adata(self,
						   counts_file,
						   how='iso'):
		"""
		Convert expression matrix to adata

		Parameters:
			counts_file (str): Path to TALON or Swan-formatted counts matrix
			how (str): {'iso', 'gene'}

		Returns:
			adata (anndata AnnData): AnnData representation of gene or
				transcript-level expression
		"""

		# read in abundance file
		check_file_loc(counts_file, 'abundance matrix')
		try:
			df = pd.read_csv(counts_file, sep='\t')
		except:
			raise ValueError('Problem reading expression matrix {}'.format(counts_file))

		# check if abundance matrix is a talon transcript or gene abundance matrix
		t_cols = ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_transcript_id',
			'annot_gene_name', 'annot_transcript_name', 'n_exons', 'length',
			'gene_novelty', 'transcript_novelty', 'ISM_subtype']
		g_cols = ['gene_ID', 'annot_gene_id', 'annot_gene_name', 'gene_novelty']
		if df.columns.tolist()[:11] == t_cols or df.columns.tolist()[:4] == g_cols:
			df = reformat_talon_abundance(df, how=how)

		# rename id ID column
		col = df.columns[0]
		if how == 'gene':
			id_col = 'gid'
		elif how == 'iso':
			id_col = 'tid'

		df.rename({col: id_col}, axis=1, inplace=True)

		# sum for gene level
		if how == 'gene':
			df = df.groupby(id_col).sum().reset_index()

		# limit to just the transcripts already in the graph
		if how == 'iso':
			sg_tids = self.t_df.tid.tolist()
			ab_tids = df.tid.tolist()
			tids = list(set(sg_tids)&set(ab_tids))
			df = df.loc[df.tid.isin(tids)]

		# transpose to get adata format
		df.set_index(id_col, inplace=True)
		df = df.T

		# get adata components - obs, var, and X
		var = df.columns.to_frame()
		var.columns = [id_col]
		obs = df.index.to_frame()
		obs.columns = ['dataset']
		X = sparse.csr_matrix(df.to_numpy())

		# create transcript-level adata object and filter out unexpressed transcripts
		adata = anndata.AnnData(var=var, obs=obs, X=X)
		genes, _  = sc.pp.filter_genes(adata, min_counts=1, inplace=False)
		adata = adata[:, genes]
		adata.layers['counts'] = adata.X

		return adata

	def merge_adata_abundance(self,
							  adata,
							  how='iso'):
		"""
		Merge abundance information with preexisting abundance info in the
			SwanGraph

		Parameters:
			adata (anndata AnnData): Output from `abundance_to_adata`
			how (str): {'iso', 'gene'}
		"""

		if how == 'gene':
			dataset_list = self.gene_datasets
			ab_bool = self.has_gene_abundance()
			sg_adata = self.gene_adata
		elif how == 'iso':
			dataset_list = self.datasets
			ab_bool = self.has_abundance()
			sg_adata = self.adata

		# add each dataset to list of "datasets", check if any are already there!
		datasets = adata.obs.dataset.tolist()
		for d in datasets:
			if d in dataset_list:
				raise ValueError('Dataset {} already present in the SwanGraph.'.format(d))
		dataset_list.extend(datasets)

		print()
		if len(datasets) <= 5:
			print('Adding abundance for datasets {} to SwanGraph.'.format(', '.join(datasets)))
		else:
			mini_datasets = datasets[:5]
			n = len(datasets) - len(mini_datasets)
			print('Adding abundance for datasets {}... (and {} more) to SwanGraph'.format(', '.join(mini_datasets), n))

		# if there is preexisting abundance data in the SwanGraph, concatenate
		# otherwise, adata is the new transcript level adata
		if not ab_bool:

			# create transcript-level adata object
			sg_adata = adata

			# add counts as layers
			sg_adata.layers['counts'] = sg_adata.X
			print('Calculating TPM...')
			sg_adata.layers['tpm'] = sparse.csr_matrix(calc_tpm(sg_adata, recalc=True).to_numpy())

			if not self.sc and how == 'iso':
				print('Calculating PI...')
				sg_adata.layers['pi'] = sparse.csr_matrix(calc_pi(sg_adata, self.t_df)[0].to_numpy())
		else:

			# first set current layer to be counts
			sg_adata.X = sg_adata.layers['counts']

			# concatenate existing adata with new one
			# outer join to add all new transcripts (that are from added
			# annotation or transcriptome) to the abundance
			uns = sg_adata.uns
			sg_adata = sg_adata.concatenate(adata, join='outer', index_unique=None)
			sg_adata.uns = uns

			# recalculate pi and tpm
			print('Calculating TPM...')
			sg_adata.layers['tpm'] = sparse.csr_matrix(calc_tpm(sg_adata, recalc=True).to_numpy())

			if not self.sc and how == 'iso':
				print('Calculating PI...')
				sg_adata.layers['pi'] = sparse.csr_matrix(calc_pi(sg_adata, self.t_df)[0].to_numpy())

		if how == 'iso':
			self.adata = sg_adata
		elif how == 'gene':
			self.gene_adata = sg_adata

		# add abundance for edges, TSS per gene, and TES per gene
		if how == 'iso':
			if self.make_edge_adata:
				print('Calculating edge usage...')
				self.create_edge_adata()
			if self.make_end_adata:
				print('Calculating TSS usage...')
				self.create_feat_adata(kind='tss')
				print('Calculating TES usage...')
				self.create_feat_adata(kind='tes')

			if 'ic_id' in self.t_df.columns:
				if self.make_ic_adata:
					print('Calculating IC usage...')
					self.create_feat_adata(kind='ic')

		# set abundance flag to true
		if how == 'iso':
			self.abundance = True
			self.adata = sg_adata
		elif how == 'gene':
			self.gene_abundance = True
			self.gene_adata = sg_adata


	def format_adata(self,
	                 adata_file,
					 how='iso'):
		"""
		Format AnnData file for use with Swan

		Parameters:
			adata_file (str): Path to AnnData file where var index is the
				transcript ID and obs index is the dataset name.
			how (str): {'iso', 'gene'}

		Returns:
			adata (anndata AnnData): AnnData representation of gene or
				transcript-level expression
		"""

		# read in abundance file
		check_file_loc(adata_file, 'AnnData')
		try:
			adata = sc.read(adata_file)
		except:
			raise ValueError('Problem reading AnnData {}'.format(counts_file))

		# format tid for var table
		adata.var['tid'] = adata.var.index
		adata.var.index.name = 'tid'

		# format dataset for obs table
		adata.obs['dataset'] = adata.obs.index

		# sum for gene level
		if how == 'gene':
			# df = df.groupby(id_col).sum().reset_index()
			raise Error('Not implemented yet')
			# will have to think about summing across genes
			# probably use tid<->gid conversions from t_df

		# create transcript-level adata object and filter out unexpressed transcripts
		genes, _  = sc.pp.filter_genes(adata, min_counts=1, inplace=False)
		adata = adata[:, genes].copy()

		# limit to just the transcripts already in the graph
		sg_tids = self.t_df.tid.tolist()
		ab_tids = adata.var.index.tolist()
		tids = list(set(sg_tids)&set(ab_tids))
		adata = adata[:, tids].copy()

		# make sure we're storing as sparse
		if isinstance(adata.X, np.ndarray):
			X = sparse.csr_matrix(adata.X)
			adata.X = X

		# name this layer
		adata.layers['counts'] = adata.X

		return adata

	def add_adata(self,
	              adata_file,
				  how='iso'):
		"""
		Adds abundance / metadata information from an AnnData object into the
		SwanGraph. Transcripts in the SwanGraph but not in the AnnData will be
		assigned 0 counts. Transcripts in the abundance matrix but not in the
		SwanGraph will not have expression added.

		Parameters:
			adata_file (str): Path to AnnData file where var index is the
				transcript ID, obs index is the dataset name, and X is the
				raw counts (ie not already log transformed).
			how (str): {'iso', 'gene'} (Gene-level currently not implemented)
		"""
		adata = self.format_adata(adata_file, how=how)
		self.merge_adata_abundance(adata, how=how)

	def add_abundance(self,
					  counts_file,
					  how='iso'):
		"""
		Adds abundance from a counts matrix to the SwanGraph. Transcripts in the
		SwanGraph but not in the counts matrix will be assigned 0 counts.
		Transcripts in the abundance matrix but not in the SwanGraph will not
		have expression added.

		Parameters:
			counts_file (str): Path to TSV expression file where first column is
				the transcript ID and following columns name the added datasets and
				their counts in each dataset, OR to a TALON abundance matrix.
			how (str): {'iso', 'gene'}
		"""
		adata = self.abundance_to_adata(counts_file, how=how)
		self.merge_adata_abundance(adata, how=how)

	def add_metadata(self,
	                 fname,
					 overwrite=False):
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
			raise ValueError('Problem reading metadata file {}'.format(fname))

		# has abundance been added yet?
		if not self.abundance:
			raise ValueError('Cannot add metadata. No datasets have been added')

		# which columns are duplicate?
		sg_cols = list(set(self.adata.obs.columns.tolist())-set(['dataset']))
		meta_cols = list(set(df.columns.tolist())-set(['dataset']))
		dupe_cols = [c for c in meta_cols if c in sg_cols]

		# drop columns from one table depending on overwrite settings
		# adatas = [self.adata, self.tss_adata, self.tes_adata, self.edge_adata]
		# adatas = [self.adata, self.tss_adata, self.tes_adata,
		# 		  self.gene_adata, self.ic_adata]
		adatas = self.get_adatas()
		if dupe_cols:
			if overwrite:
				for adata in adatas:
					if overwrite:
						adata.obs.drop(dupe_cols, axis=1, inplace=True)
			else:
				df.drop(dupe_cols, axis=1, inplace=True)

		# for each anndata
		for adata in adatas:

			# merge df with adata obs table
			adata.obs = adata.obs.merge(df, how='left', on='dataset')

			# make sure all dtypes in obs table are non intergers
			for ind, entry in adata.obs.dtypes.to_frame().iterrows():
				if entry[0] == 'int64':
					adata.obs[ind] = adata.obs[ind].astype('str')

			# and set index to dataset
			adata.obs['index'] = adata.obs.dataset
			adata.obs.set_index('index', inplace=True)

	##########################################################################
	############# Obtaining abundance of edges, locs, and ends ###############
	##########################################################################
	def create_cerberus_adata(self, mode):
		"""
		Use the ends and intron chains called by Cerberus to compute
		expression values for TSSs, TESs, or intron chains.

		Parameters:
			mode (str): {'ic', 'tss', 'tes'}
		"""
		id_col = '{}_id'.format(mode)
		name_col = '{}_name'.format(mode)

		# merge ic, tss, tes info with transcript counts
		gb_cols = ['gid', 'gname', id_col]
		subset_cols = copy.deepcopy(gb_cols)
		if mode in ['tss', 'tes']:
			end = True
			vert_col = '{}_vertex'.format(mode)
			subset_cols.append('loc_path')
			if mode == 'tss':
				ind = 0
			elif mode == 'tes':
				ind = -1
		else:
			end = False

		df = self.t_df[subset_cols].reset_index()

		# merge with vertex id
		if end:
			df[vert_col] = [path[ind] for path in df.loc_path.values.tolist()]
			df.drop('loc_path', axis=1, inplace=True)
			gb_cols.append(vert_col)

		t_counts = self.get_transcript_abundance(kind='counts')
		t_counts = t_counts.merge(df, how='left', on='tid')

		# gb and sum counts over the different features
		t_counts.drop('tid', axis=1, inplace=True)
		t_counts = t_counts.groupby(gb_cols).sum()

		# get separate components of the table
		X = sparse.csr_matrix(t_counts.transpose().values)
		obs = self.adata.obs

		# add name of thing
		t_counts.reset_index(inplace=True)
		t_counts[name_col] = t_counts.gname+'_'+t_counts[id_col].str.split('_', expand=True)[1]
		gb_cols.append(name_col)
		var = t_counts[gb_cols]
		var.set_index(id_col, inplace=True)

		if end:
			var.rename({vert_col: 'vertex_id'}, axis=1, inplace=True)

		return obs, var, X

	def create_feat_adata(self, kind):
		"""
		Create a tss / tes / ic-level adata object. Enables calculating tss / tes
		usage across samples.

		Parameters:
			kind (str): {'tss', 'ic', 'tes'}
		"""

		# check to see if end information is already present in swangraph
		id_col = '{}_id'.format(kind)
		if id_col in self.t_df.columns:
			print('Using Cerberus IDs to calculate')
			obs, var, X = self.create_cerberus_adata(kind)

		else:
			# limit to only expresed transcripts
			t_df = self.t_df.copy(deep=True)
			t_df = t_df.loc[self.adata.var.index.tolist()]
			df = get_ends(t_df, kind)

			# get a mergeable transcript expression df
			tid = self.adata.var.index.tolist()
			obs = self.adata.obs.index.tolist()
			data = self.adata.layers['counts'].transpose()
			t_exp_df = pd.DataFrame.sparse.from_spmatrix(columns=obs, data=data, index=tid)
			t_exp_df = t_exp_df.merge(t_df, how='left',
				left_index=True, right_index=True)

			# merge counts per transcript with end expression
			df = df.merge(t_exp_df, how='left',
				left_index=True, right_index=True)

			# sort based on vertex id
			df.sort_index(inplace=True, ascending=True)

			# set index to gene ID, gene name, and vertex id
			df.reset_index(drop=True, inplace=True)
			df.set_index(['gid', 'gname', 'vertex_id'], inplace=True)
			df = df[self.datasets]

			# groupby on gene and assign each unique TSS / gene combo an ID
			# use dense representation b/c I've already removed 0 counts and sparse
			# gb operations are known to be slow in pandas
			# https://github.com/pandas-dev/pandas/issues/36123
			# maybe try this? :
			# https://cmdlinetips.com/2019/03/how-to-write-pandas-groupby-function-using-sparse-matrix/
			id_col = '{}_id'.format(kind)
			name_col = '{}_name'.format(kind)
			df = df.copy()
			df.reset_index(inplace=True)
			df[self.datasets] = df[self.datasets].sparse.to_dense()
			df = df.groupby(['gid', 'gname', 'vertex_id']).sum().reset_index()

			df['end_gene_num'] = df.sort_values(['gid', 'vertex_id'],
							ascending=[True, True])\
							.groupby(['gid']) \
							.cumcount() + 1
			df[id_col] = df['gid']+'_'+df['end_gene_num'].astype(str)
			df[name_col] = df['gname']+'_'+df['end_gene_num'].astype(str)
			df.drop('end_gene_num', axis=1, inplace=True)

			# obs, var, and X tables for new data
			var_cols = ['gid', 'gname', 'vertex_id', id_col, name_col]
			var = df[var_cols]
			var.set_index('{}_id'.format(kind), inplace=True)
			df.drop(var_cols, axis=1, inplace=True)
			df = df[self.adata.obs.index.tolist()]
			X = sparse.csr_matrix(df.transpose().values)
			obs = self.adata.obs

		# create anndata
		adata = anndata.AnnData(var=var, obs=obs, X=X)

		# add counts and tpm as layers
		adata.layers['counts'] = adata.X
		adata.layers['tpm'] = sparse.csr_matrix(calc_tpm(adata, recalc=True).to_numpy())
		if not self.sc:
			adata.layers['pi'] = sparse.csr_matrix(calc_pi(adata,
					adata.var)[0].to_numpy())

		# assign adata and clean up unstructured data if needed
		if kind == 'tss':
			if self.has_abundance():
				adata.uns = self.tss_adata.uns
			self.tss_adata = adata

		elif kind == 'tes':
			if self.has_abundance():
				adata.uns = self.tss_adata.uns
			self.tes_adata = adata

		elif kind == 'ic':
			if self.has_abundance():
				adata.uns = self.ic_adata.uns
			self.ic_adata = adata

	def create_edge_adata(self):
		"""
		Create an edge-level adata object. Enables calculating edge usage across
		samples.
		"""

		# get table what edges are in each transcript
		t_df = self.t_df.copy(deep=True)
		t_df = t_df.loc[self.adata.var.index.tolist()]
		edge_exp_df = pivot_path_list(t_df, 'path')

		# get a mergeable transcript expression df
		tid = self.adata.var.index.tolist()
		obs = self.adata.obs.index.tolist()
		data = self.adata.layers['counts'].transpose()
		t_exp_df = pd.DataFrame.sparse.from_spmatrix(columns=obs,
													 data=data,
													 index=tid)

		# merge counts per transcript with edges
		edge_exp_df = edge_exp_df.merge(t_exp_df, how='left',
			left_index=True, right_index=True)

		# sum the counts per transcript / edge / dataset
		edge_exp_df.reset_index(inplace=True, drop=True)
		edge_exp_df.set_index('edge_id', inplace=True)
		edge_exp_df[self.datasets] = edge_exp_df[self.datasets].sparse.to_dense()
		edge_exp_df = edge_exp_df.groupby(by='edge_id', as_index=True).sum()

		# order based on order of edges in self.edge_df
		edge_exp_df = edge_exp_df.merge(self.edge_df[['v1', 'v2']],
			how='left', left_index=True, right_index=True)
		edge_exp_df.sort_values(by=['v1', 'v2'], inplace=True)
		edge_exp_df.drop(['v1', 'v2'], axis=1, inplace=True)

		# drop edges that are unexpressed
		edge_exp_df = edge_exp_df.loc[edge_exp_df.sum(1) > 0]

		# obs, var, and X tables for new data
		var = edge_exp_df.index.to_frame()
		X = sparse.csr_matrix(edge_exp_df.transpose().values)
		obs = self.adata.obs

		# create edge-level adata object
		adata = anndata.AnnData(var=var, obs=obs, X=X)

		# add counts and tpm as layers
		adata.layers['counts'] = adata.X
		adata.layers['tpm'] = sparse.csr_matrix(calc_tpm(adata, recalc=True).to_numpy())

		# can't make pi for edges unless I make a new edge for
		# each gene that the edge is in
		# could just have self.edge_adata var separate from self.edge_df for now tho
	#	 self.edge_adata.layers['pi'] = sparse.csr_matrix(calc_pi(self.adata, self.edge_df)[0].to_numpy())

		# assign adata and clean up unstructured data if needed
		if self.has_abundance():
			adata.uns = self.edge_adata.uns
		self.edge_adata = adata

	##########################################################################
	############# Related to creating dfs from GTF or TALON DB ###############
	##########################################################################

	def get_current_locs(self):
		"""
		Get a dictionary of unique locations already in the SwanGraph, along
		with the number of unique locations. Used to create the dataframes
		in create_dfs.

		Returns:
			locs (dict): Dictionary of unique genomic locations. Keys are
				(chromosome, coordinate). Items are the vertex IDs from the
				SwanGraph.
			n (int): Number of unique genomic locations. -1 if SwanGraph
				is empty.
		"""
		if not self.is_empty():
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

		if not self.is_empty():

			# number of edges currently
			n = self.edge_df.edge_id.max()

			# get a version of edge_df with the coordinates
			edge_df = self.add_edge_coords()

			edges = {}
			for ind, entry in edge_df.iterrows():
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
		else:
			edges = {}
			n = -1

		return edges, n

	def add_edge_coords(self):
		"""
		Add coordinates and chromosome to edge_df from loc_df.

		Returns:
			edge_df (pandas DataFrame): Copy of self.edge_df with genomic
				coordinate information
		"""
		temp = self.loc_df[['chrom', 'coord', 'vertex_id']]
		temp.reset_index(inplace=True, drop=True)

		edge_df = self.edge_df.copy(deep=True)

		edge_df = edge_df.merge(temp, how='left',
			left_on='v1', right_on='vertex_id')
		edge_df.rename({'coord': 'v1_coord'}, axis=1, inplace=True)

		temp = self.loc_df[['coord', 'vertex_id']]
		temp.reset_index(inplace=True, drop=True)

		edge_df = edge_df.merge(temp, how='left',
			left_on='v2', right_on='vertex_id')
		edge_df.rename({'coord': 'v2_coord'}, axis=1, inplace=True)

		return edge_df

	def create_dfs(self,
	               transcripts,
				   exons,
				   from_talon=False,
				   from_cerberus=False):
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
			from_cerberus (bool): Whether or not GTF is from cerberus

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

		cols = ['tid', 'tname', 'gid', 'gname', 'path']
		if from_talon:
			cols.append('novelty')
		if from_cerberus:
			cols += ['tss_id', 'ic_id', 'tes_id']
		t_df = pd.DataFrame(transcripts)
		t_df = t_df.transpose()
		t_df = t_df.reset_index(drop=True)
		t_df = t_df[cols]

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
		if 'novelty' in t_df.columns:
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
			loc (dict): TODO

		Returns:
			transcripts (dict): Dictionary of only NEW transcripts.
				Keys (str): transcript ID
				Items (str): TODO
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

	def label_annotated(self, tids):
		"""
		From a list of transcript IDs that were added as part of the annotation,
		add a boolean column to each DataFrame in the SwanGraph indicating if
		that location/edge/transcript is or is not in the SwanGraph.

		Parameters:
			tids (list of str): List of transcript IDs from t_df that were parts
				of the added annotation
		"""

		# obtain all unique annotated edges
		edges = self.t_df.loc[tids, 'path'].tolist()
		edges = list(set([e for path in edges for e in path]))

		# obtain all unique annotated locations
		v1_locs = self.edge_df.loc[edges, 'v1'].tolist()
		v2_locs = self.edge_df.loc[edges, 'v2'].tolist()
		locs = list(set(v1_locs+v2_locs))

		# set the annotated column to True for these tids, edges, and locs
		self.t_df['annotation'] = False
		self.t_df.loc[tids, 'annotation'] = True
		self.edge_df['annotation'] = False
		self.edge_df.loc[edges, 'annotation'] = True
		self.loc_df['annotation'] = False
		self.loc_df.loc[locs, 'annotation'] = True

	##########################################################################
	######################## Other SwanGraph utilities #####################
	##########################################################################
	def subset_on_gene_sg(self, gid=None, datasets=None):
		"""
		Subset the swan Graph on a given gene and return the subset graph.

		Parameters:
			gid (str): Gene ID to subset on
			datasets (list of str): List of datasets to keep in the subset

		returns:
			subset_sg (swan Graph): Swan Graph subset on the input gene.
		"""

		# didn't ask for either
		if not gid and not datasets:
			return self

		# subset on gene
		if gid:
			# make sure this gid is even in the Graph
			self.check_gene(gid)

			# get the strand
			strand = self.get_strand_from_gid(gid)

			# subset t_df first, it's the easiest
			tids = self.t_df.loc[self.t_df.gid == gid].index.tolist()
			t_df = self.t_df.loc[tids].copy(deep=True)
			t_df['path'] = self.t_df.loc[tids].apply(
					lambda x: copy.deepcopy(x.path), axis=1)
			t_df['loc_path'] = self.t_df.loc[tids].apply(
					lambda x: copy.deepcopy(x.loc_path), axis=1)

			# since we don't keep all transcripts in adata, make
			# sure to pare that down
			tids = list(set(tids)&set(self.adata.var.index.tolist()))

			# subset loc_df based on all the locs that are in the paths from
			# the already-subset t_df
			paths = t_df['loc_path'].tolist()
			locs = [node for path in paths for node in path]
			locs = np.unique(locs)
			loc_df = self.loc_df.loc[locs].copy(deep=True)

			# subset edge_df based on all the edges that are in the paths from
			# the alread-subset t_df
			paths = t_df['path'].tolist()
			edges = [node for path in paths for node in path]
			edges = np.unique(edges)
			edge_df = self.edge_df.loc[edges].copy(deep=True)
		if not gid:
			t_df = self.t_df.copy(deep=True)
			edge_df = self.edge_df.copy(deep=True)
			loc_df = self.loc_df.copy(deep=True)

		# also subset anndata
		# if obs_col and obs_cats:
		# 	# adatas = [self.adata, self.edge_adata,
		# 	# 		  self.tss_adata, self.tes_adata]
		# 	# for adata in adatas:
		# 	# print(obs_col)
		# 	# print(obs_cats)
		# 	obs_vars = self.adata.obs.loc[self.adata.obs[obs_col].isin(obs_cats)]
		# 	obs_vars = obs_vars.index.tolist()
		# 	# print(obs_vars)
		# 	# print(tids)
		# 	adata = self.adata[obs_vars, tids]
		new_adatas = dict()
		# adatas = {'iso': self.adata, 'edge': self.edge_adata,
		# 		  'tss': self.tss_adata, 'tes': self.tes_adata}
		adatas = {'iso': self.adata}
		for key, adata in adatas.items():
			if datasets and gid:
				new_adatas[key] = adata[datasets, tids]
			elif gid:
				new_adatas[key] = adata[:, tids]
			elif datasets:
				new_adatas[key] = adata[datasets, :]
			else:
				new_adatas[key] = adata

		# create a new graph that's been subset
		subset_sg = SwanGraph()
		subset_sg.loc_df = loc_df
		subset_sg.edge_df = edge_df
		subset_sg.t_df = t_df
		subset_sg.adata = new_adatas['iso']
		# subset_sg.edge_adata = new_adatas['edge']
		# subset_sg.tss_adata = new_adatas['tss']
		# subset_sg.tes_adata = new_adatas['tes']
		subset_sg.datasets = subset_sg.adata.obs.index.tolist()
		subset_sg.abundance = self.abundance
		subset_sg.sc = self.sc
		subset_sg.pg = self.pg
		subset_sg.annotation = self.annotation

		# renumber locs if using a gene
		if gid:
			if strand == '-':
				id_map = subset_sg.get_ordered_id_map(rev_strand=True)
				subset_sg.update_ids(id_map)
			else:
				subset_sg.update_ids()

			subset_sg.get_loc_types()

		# finally create the graph
		subset_sg.create_graph_from_dfs()

		return subset_sg

	def order_transcripts_subset(self, t_df, order='tid'):
		"""
		Order the transcripts in an input t_df based on an input heuristic.
		Can order alphabetically by transcript ID, by expression of each
		transcript, or by the genomic location of the transcripts' TSSs or TESs.

		Parameters:
			t_df (pandas DataFrame):
			order (str): Method to order transcripts by. Choose from ['tid',
				'expression', 'tss', 'tes', 'log2tpm']

		Returns:
			t_df (pandas DataFrame): Transcripts ordered by the input heuristic
			tids (list of str): List of transcript IDs ordered by the input
				heuristic
		"""

		# order by transcript id
		# TODO
		if order == 'tid':
			ordered_tids = sorted(t_df.index.tolist())
			t_df = t_df.loc[ordered_tids]

		# order by expression
		elif order == 'expression':
			# make sure there are counts in the graph at all
			if not self.abundance:
				raise Exception('Cannot order by expression because '
								'there is no expression data.')

			t_df['sum'] = t_df.sum(axis=1)
			t_df.sort_values(by='sum', ascending=False, inplace=True)
			t_df.drop('sum', axis=1, inplace=True)

		# order by log2tpm expression
		elif order == 'log2tpm':
			# make sure there are counts in the graph at all
			if not self.abundance:
				raise Exception('Cannot order by expression because '
								'there is no expression data.')

			t_df['sum'] = np.log2(t_df+1).sum(axis=1)
			t_df.sort_values(by='sum', ascending=False, inplace=True)
			t_df.drop('sum', axis=1, inplace=True)

		# order by cumulative pi value
		elif order == 'pi':
			t_df['sum'] = t_df.sum(axis=1)
			t_df.sort_values(by='sum', ascending=False, inplace=True)
			t_df.drop('sum', axis=1, inplace=True)

		# order by tss location
		elif order == 'tss':
			t_df = t_df.merge(self.t_df[['path', 'loc_path']],
				how='left', left_index=True, right_index=True)
			t_df['start_coord'] = t_df.apply(lambda x:
				self.loc_df.loc[x.loc_path[0], 'coord'], axis=1)
			t_df['strand'] = t_df.apply(lambda x:
				self.edge_df.loc[x.path[0], 'strand'], axis=1)
			fwd = t_df.loc[t_df.strand == '+']
			rev = t_df.loc[t_df.strand == '-']
			fwd.sort_values(by='start_coord', ascending=True, inplace=True)
			# print(fwd[['strand', 'start_coord']])
			rev.sort_values(by='start_coord', ascending=False, inplace=True)
			# print(rev[['strand', 'start_coord']])
			t_df = pd.concat([fwd, rev])
			# print(t_df[['strand', 'start_coord']])
			t_df.drop(['start_coord', 'strand', 'path', 'loc_path'], axis=1, inplace=True)

		# order by tes location
		elif order == 'tes':
			t_df = t_df.merge(self.t_df[['path', 'loc_path']],
				how='left', left_index=True, right_index=True)
			t_df['end_coord'] = t_df.apply(lambda x:
				self.loc_df.loc[x.loc_path[-1], 'coord'], axis=1)
			t_df['strand'] = t_df.apply(lambda x:
				self.edge_df.loc[x.path[-1], 'strand'], axis=1)
			fwd = t_df.loc[t_df.strand == '+']
			rev = t_df.loc[t_df.strand == '-']
			fwd.sort_values(by='end_coord', ascending=True, inplace=True)
			rev.sort_values(by='end_coord', ascending=False, inplace=True)
			t_df = pd.concat([fwd, rev])
			t_df.drop(['end_coord', 'strand', 'path', 'loc_path'], axis=1, inplace=True)

		tids = t_df.index.tolist()
		return t_df, tids

	def order_transcripts(self, order='tid'):
		"""
		Order the transcripts in the SwanGraph t_df based on an input heuristic.
		Can order alphabetically by transcript ID, by expression of each
		transcript, or by the genomic location of the transcripts' TSSs or TESs.

		Parameters:
			order (str): Method to order transcripts by. {'tid',
				'expression', 'tss', 'tes'}
				Default: 'tid'
		"""

		# order by transcript id
		if order == 'tid':
			ordered_tids = sorted(self.t_df.tid.tolist())
			self.t_df = self.t_df.loc[ordered_tids]

		# order by expression
		elif order == 'expression':

			# make sure there are counts in the graph at all
			if not self.abundance:
				raise Exception('Cannot order by expression because '
								'there is no expression data.')
			else:

				# find max expressed transcripts
				data = self.adata.layers['tpm'].sum(axis=0)
				cols = self.adata.var.index.tolist()
				temp = pd.DataFrame(index=['tpm'], data=[data], columns=cols).transpose()
				temp.sort_values(by='tpm', ascending=False, inplace=True)
				ordered_tids = temp.index.tolist()
				# print(temp)

				# order t_df and adata
				self.t_df = self.t_df.loc[ordered_tids]
				# print(self.t_df)
				self.adata = self.adata[:, ordered_tids]

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

	def find_ir_genes(self, verbose=False):
		"""
		Finds all unique genes containing novel intron retention events.
		Requires that an annotation has been added to the SwanGraph.

		Parameters:
			verbose (bool): Display output

		Returns:
			ir_df (pandas DataFrame): DataFrame detailing discovered novel
				intron retention edges and the transcripts and genes they
				come from
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

		# verbose output
		if verbose:
			n_progress = len(edge_ids)
			pbar = tqdm(total=n_progress)
			pbar.set_description('Testing each novel edge for intron retention')

		# for each edge, see if the subgraph between the edge vertices
		# contains an exonic edge
		ir_edges = []
		ir_genes = []
		ir_transcripts = []
		# df for gene+transcript+edge combos
		ir_df = pd.DataFrame()
		for i, eid in enumerate(edge_ids):
			# subgraph consisting of all nodes between the candidate
			# intron-retaining edge coords in order and its edges
			entry = self.edge_df.loc[eid]
			v1 = entry.v1
			v2 = entry.v2
			sub_nodes = [i for i in range(v1+1,v2)]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			self.edge_df['tuple_edge_id'] = self.edge_df[['v1', 'v2']].apply(tuple, axis=1)
			sub_edges = self.edge_df.loc[self.edge_df.tuple_edge_id.isin(sub_edges)]
			self.edge_df.drop('tuple_edge_id', axis=1, inplace=True)

			# find edges that are intronic; if there are none, this is not
			# an intron-retaining edge
			sub_edges = sub_edges.loc[sub_edges.edge_type == 'intron']

			if len(sub_edges.index) > 0:

				# transcripts that contain the exon-skipping edge
				cand_t_df = nt_df[[eid in path for path in nt_df.path.values.tolist()]]

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
							temp_df = cand_g_df[[cand_eid in path for path in cand_g_df.path.values.tolist()]]
							tids = cand_t_df.tid.tolist()
							for tid in tids:
								temp = pd.DataFrame(data=[[gid, tid, eid]],
									columns=['gid', 'tid', 'edge_id'])
								ir_df = pd.concat([ir_df, temp])
								ir_edges.append(eid)
								ir_genes.append(gid)
								ir_transcripts.append(tid)
			if verbose:
				pbar.update(1)

		# ir_genes = list(set(ir_genes))
		# ir_transcripts = list(set(ir_transcripts))
		# ir_edges = list(set(ir_edges))
		ir_transcripts = ir_df.tid.unique().tolist()
		print('Found {} novel ir events in {} transcripts.'.format(len(ir_df.index),
			len(ir_transcripts)))

		# store in 'ir' field of self.adata.uns
		self.adata.uns['ir'] = ir_df

		return ir_df

	def find_es_genes(self, verbose=False):
		"""
		Finds all unique genes containing novel exon skipping events.
		Requires that an annotation has been added to the SwanGraph.

		Parameters:
			verbose (bool): Display output

		Returns:
			es_df (pandas DataFrame): DataFrame detailing discovered novel
				exon-skipping edges and the transcripts and genes they
				come from
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

		# verbose output
		if verbose:
			n_progress = len(edge_ids)
			pbar = tqdm(total=n_progress)
			pbar.set_description('Testing each novel edge for exon skipping')

		# for each edge, see if the subgraph between the edge vertices
		# contains an exonic edge
		es_edges = []
		es_genes = []
		es_transcripts = []
		# df for gene+transcript+edge combos
		es_df = pd.DataFrame()
		for eid in edge_ids:
			# subgraph consisting of all nodes between the candidate
			# exon-skipping edge coords in order and its edges
			entry = self.edge_df.loc[eid]
			v1 = entry.v1
			v2 = entry.v2
			sub_nodes = [i for i in range(v1+1,v2)]
			sub_G = self.G.subgraph(sub_nodes)
			sub_edges = list(sub_G.edges())
			self.edge_df['tuple_edge_id'] = self.edge_df[['v1', 'v2']].apply(tuple, axis=1)
			sub_edges = self.edge_df.loc[self.edge_df.tuple_edge_id.isin(sub_edges)]
			self.edge_df.drop('tuple_edge_id', axis=1, inplace=True)

			# find edges that are exonic; if there are none, this is not
			# an exon-skipping edge
			sub_edges = sub_edges.loc[sub_edges.edge_type == 'exon']

			if len(sub_edges.index) > 0:

				# transcripts that contain the candidate exon-skipping edge
				skip_t_df = nt_df[[eid in path for path in nt_df.path.values.tolist()]]

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
							temp_df = skip_g_df[[skip_eid in path for path in skip_g_df.path.values.tolist()]]
							tids = skip_t_df.tid.tolist()
							if len(temp_df.index) > 0:
								for tid in tids:
									temp = pd.DataFrame(data=[[gid, tid, eid]],
										columns=['gid', 'tid', 'edge_id'])
									es_df = pd.concat([es_df, temp])
									es_edges.append(eid)
									es_genes.append(gid)
									es_transcripts.append(tid)
			if verbose:
				pbar.update(1)

		# es_genes = list(set(es_genes))
		# es_transcripts = list(set(es_transcripts))
		# es_edges = list(set(es_edges))
		es_transcripts = es_df.tid.unique().tolist()
		print('Found {} novel es events in {} transcripts.'.format(len(es_df.index),
			len(es_transcripts)))

		# store in 'es' field of self.adata.uns
		self.adata.uns['es'] = es_df

		return es_df

	# def de_gene_test(self, obs_col, obs_conditions=None):
	# 	"""
	# 	Runs a differential expression test on the gene level.
	#
	# 	Parameters:
	# 		obs_col (str): Metadata column from self.adata.obs to perform
	# 			the test on
	# 		obs_conditions (list of str, len 2): Categories from 'obs_col' in
	# 			self.adata.obs to perform the test on.
	#
	# 	Returns:
	# 		test (pandas DataFrame): A summary table of the differential
	# 			expression test, including p and q-values, as well
	# 			as log fold change.
	# 	"""

		# # check if obs_col is even there
		# if obs_col not in self.adata.obs.columns.tolist():
		# 	raise ValueError('Metadata column {} not found'.format(obs_col))
		#
		# # check if there are more than 2 unique values in obs_col
		# if len(self.adata.obs[obs_col].unique().tolist())!=2 and not obs_conditions:
		# 	raise ValueError('Must provide obs_conditions argument when obs_col has'\
		# 		' >2 unique values')
		# elif not obs_conditions:
		# 	obs_conditions = self.adata.obs[obs_col].unique().tolist()
		# elif len(obs_conditions) != 2:
		# 	raise ValueError('obs_conditions must have exactly 2 values')
		#
		# # check if these values of obs_col exist
		# if obs_col and obs_conditions:
		# 	conds = self.adata.obs[obs_col].unique().tolist()
		# 	for cond in obs_conditions:
		# 		if cond not in conds:
		# 			raise ValueError('Value {} not found in metadata column {}.'.format(cond, obs_col))
		#
		# # get an AnnData that's a subset based on obs_conditions
		# if obs_conditions:
		# 	inds = self.adata.obs.loc[self.adata.obs[obs_col].isin(obs_conditions)]
		# 	inds = inds.index.tolist()
		# 	adata = self.adata[inds, :]
		# else:
		# 	adata = self.adata
		#
		# # groupby gene and sum counts
		# data = adata.layers['counts'].toarray()
		# var = adata.var.index.tolist()
		# obs = adata.obs.index.tolist()
		# df = pd.DataFrame(data=data, columns=var, index=obs)
		# df = df.transpose()
		# df = df.merge(self.t_df['gid'], how='left',
		# 	left_index=True, right_index=True)
		# df.reset_index(drop=True, inplace=True)
		# df = df.groupby('gid').sum()
		#
		# # format as adata and add metadata
		# df = df.transpose()
		# X = sparse.csr_matrix(df.values)
		# var = df.columns.to_frame()
		# obs = df.index.to_frame()
		# obs.rename({0: 'dataset'}, axis=1, inplace=True)
		# adata = anndata.AnnData(X=X, obs=obs, var=var)
		# adata.obs = adata.obs.merge(self.adata.obs, how='left', on='dataset')
		# adata.obs.index = adata.obs.dataset
		# adata.layers['counts'] = adata.X
		#
		# # recalculate TPM
		# adata.layers['tpm'] = sparse.csr_matrix(calc_tpm(adata, recalc=True).to_numpy())
		#
		# # set layer to TPM
		# adata.X = adata.layers['tpm']
		#
		# formula_loc = "~ 1 + {}".format(obs_col)
		# factor_loc_totest = obs_col
		#
		# # test
		# test = de.test.wald(
		# 	data=adata,
		# 	formula_loc=formula_loc,
		# 	factor_loc_totest=factor_loc_totest)
		# test = test.summary()
		# test.rename({'gene': 'gid'}, axis=1, inplace=True)
		#
		# # add gene name column
		# gnames = self.t_df[['gid', 'gname']].copy(deep=True)
		# gnames.drop_duplicates(inplace=True)
		# gnames.reset_index(drop=True, inplace=True)
		# test = test.merge(gnames, how='left', on='gid')
		#
		# # sort on log2fc
		# test = test.reindex(test.log2fc.abs().sort_values(ascending=False).index)
		#
		# # add summary table to anndata
		# uns_name = make_uns_key(kind='deg',
		# 			obs_col=obs_col, obs_conditions=obs_conditions)
		#
		# self.adata.uns[uns_name] = test
		#
		# return test

	# def get_de_genes(self, obs_col, obs_conditions=None, q=0.05, log2fc=1):
	# 	"""
	# 	Subsets the differential gene expression test summary table based
	# 	on a q-value and log2fc cutoff. Requires that de_gene_test has
	# 	already been run.
	#
	# 	Parameters:
	# 		obs_col (str): Column name from self.adata.obs table to group on.
	# 			Default: 'dataset'
	# 		obs_conditions (list of str, len 2): Which conditions from obs_col
	# 			to compare? Required if obs_col has more than 2 unique values.
	# 		q (float): q-value threshold to declare a gene as significant
	# 			Default: 0.05
	# 		log2fc (float): log2fc threshold to declare a gene significant
	# 			Default: 1
	#
	# 	Returns:
	# 		test (pandas DataFrame): Summary table of genes that pass the
	# 			significance threshold
	# 	"""
	#
	# 	# check to see if this was even tested
	# 	uns_name = make_uns_key(kind='deg',
	# 				obs_col=obs_col, obs_conditions=obs_conditions)
	# 	try:
	# 		test = self.adata.uns[uns_name]
	# 	except:
	# 		raise Exception('Problem accessing DEG test results for {}, {}'.format(obs_col, obs_conditions))
	#
	# 	# subset on q value
	# 	test = test.loc[test.qval <= q].copy(deep=True)
	# 	test = test.loc[test.log2fc >= log2fc]
	#
	# 	return test

	# def de_transcript_test(self, obs_col, obs_conditions=None):
	# 	"""
	# 	Runs a differential expression test on the transcript level.
	#
	# 	Parameters:
	# 		obs_col (str): Metadata column from self.adata.obs to perform
	# 			the test on
	# 		obs_conditions (list of str, len 2): Categories from 'obs_col' in
	# 			self.adata.obs to perform the test on.
	#
	# 	Returns:
	# 		test (pandas DataFrame): A summary table of the differential
	# 			expression test, including p and q-values, as well
	# 			as log fold change.
	# 	"""
	#
	# 	# check if obs_col is even there
	# 	if obs_col not in self.adata.obs.columns.tolist():
	# 		raise ValueError('Metadata column {} not found'.format(obs_col))
	#
	# 	# check if there are more than 2 unique values in obs_col
	# 	if len(self.adata.obs[obs_col].unique().tolist())!=2 and not obs_conditions:
	# 		raise ValueError('Must provide obs_conditions argument when obs_col has'\
	# 			' >2 unique values')
	# 	elif not obs_conditions:
	# 		obs_conditions = self.adata.obs[obs_col].unique().tolist()
	# 	elif len(obs_conditions) != 2:
	# 		raise ValueError('obs_conditions must have exactly 2 values')
	#
	# 	# check if these values of obs_col exist
	# 	if obs_col and obs_conditions:
	# 		conds = self.adata.obs[obs_col].unique().tolist()
	# 		for cond in obs_conditions:
	# 			if cond not in conds:
	# 				raise ValueError('Value {} not found in metadata column {}.'.format(cond, obs_col))
	#
	# 	# get an AnnData that's a subset based on obs_conditions
	# 	if obs_conditions:
	# 		inds = self.adata.obs.loc[self.adata.obs[obs_col].isin(obs_conditions)]
	# 		inds = inds.index.tolist()
	# 		adata = self.adata[inds, :]
	# 	else:
	# 		adata = self.adata
	#
	# 	# set layer to TPM
	# 	adata.X = adata.layers['tpm']
	#
	# 	formula_loc = "~ 1 + {}".format(obs_col)
	# 	factor_loc_totest = obs_col
	#
	# 	# test
	# 	test = de.test.wald(
	# 		data=adata,
	# 		formula_loc=formula_loc,
	# 		factor_loc_totest=factor_loc_totest)
	# 	test = test.summary()
	# 	test.rename({'gene': 'tid'}, axis=1, inplace=True)
	#
	# 	# add gene name column
	# 	gnames = self.t_df[['tid', 'gid', 'gname']].copy(deep=True)
	# 	gnames.reset_index(drop=True, inplace=True)
	# 	test = test.merge(gnames, how='left', on='tid')
	#
	# 	# sort on log2fc
	# 	test = test.reindex(test.log2fc.abs().sort_values(ascending=False).index)
	#
	# 	# add summary table to anndata
	# 	uns_name = make_uns_key(kind='det',
	# 				obs_col=obs_col, obs_conditions=obs_conditions)
	#
	# 	self.adata.uns[uns_name] = test
	#
	# 	return test

	# def get_de_transcripts(self, obs_col, obs_conditions=None, q=0.05, log2fc=1):
	# 	"""
	# 	Subsets the differential gene expression test summary table based
	# 	on a q-value and log2fc cutoff. Requires that de_gene_test has
	# 	already been run.
	#
	# 	Parameters:
	# 		obs_col (str): Column name from self.adata.obs table to group on.
	# 			Default: 'dataset'
	# 		obs_conditions (list of str, len 2): Which conditions from obs_col
	# 			to compare? Required if obs_col has more than 2 unique values.
	# 		q (float): q-value threshold to declare a gene as significant
	# 			Default: 0.05
	# 		log2fc (float): log2fc threshold to declare a gene significant
	# 			Default: 1
	#
	# 	Returns:
	# 		test (pandas DataFrame): Summary table of genes that pass the
	# 			significance threshold
	# 	"""
	#
	# 	# check to see if this was even tested
	# 	uns_name = make_uns_key(kind='det',
	# 				obs_col=obs_col, obs_conditions=obs_conditions)
	# 	try:
	# 		test = self.adata.uns[uns_name]
	# 	except:
	# 		raise Exception('Problem accessing DET test results for {}, {}'.format(obs_col, obs_conditions))
	#
	# 	# subset on q value
	# 	test = test.loc[test.qval <= q].copy(deep=True)
	# 	test = test.loc[test.log2fc >= log2fc]
	#
	# 	return test

	def die_gene_test(self,
					  kind='iso',
					  obs_col='dataset',
					  obs_conditions=None,
					  rc_thresh=10,
					  verbose=False):
		"""
		Finds genes with differential isoform expression between two conditions
		that are in the obs table. If there are more than 2 unique values in
		`obs_col`, the specific categories must be specified in `obs_conditions`

		Parameters:
			kind (str): What level you would like to run the test on.
				{'iso', 'tss', 'tes', 'ic'}.
				Default: 'iso'
			obs_col (str): Column name from self.adata.obs table to group on.
				Default: 'dataset'
			obs_conditions (list of str, len 2): Which conditions from obs_col
				to compare? Required if obs_col has more than 2 unqiue values.
			rc_thresh (int): Number of reads required for each conditions
				in order to test the gene.
				Default: 10
			verbose (bool): Display progress

		Returns:
			test (pandas DataFrame): A summary table of the differential
				isoform expression test, including p-values and adjusted
				p-values, as well as change in percent isoform usage (dpi) for
				all tested genes.
			test_results (pandas DataFrame): A table indicating the test outcome
				for each gene regardless of whether it was actually tested.
		"""
		# get correct adata
		if kind == 'iso':
			adata = self.adata
			ref_df = self.t_df
		elif kind == 'tss':
			adata = self.tss_adata
			ref_df = self.tss_adata.var
		elif kind == 'tes':
			adata = self.tes_adata
			ref_df = self.tes_adata.var
		elif kind == 'ic':
			adata = self.ic_adata
			ref_df = self.ic_adata.var

		# check if obs_col is even there
		if obs_col not in adata.obs.columns.tolist():
			raise ValueError('Metadata column {} not found'.format(obs_col))

		# check if there are more than 2 unique values in obs_col
		if len(adata.obs[obs_col].unique().tolist())!=2 and not obs_conditions:
			raise ValueError('Must provide obs_conditions argument when obs_col has'\
				' >2 unique values')
		elif not obs_conditions:
			obs_conditions = adata.obs[obs_col].unique().tolist()
		elif len(obs_conditions) != 2:
			raise ValueError('obs_conditions must have exactly 2 values')

		# check if these values of obs_col exist
		if obs_col and obs_conditions:
			conds = adata.obs[obs_col].unique().tolist()
			for cond in obs_conditions:
				if cond not in conds:
					raise ValueError('Value {} not found in metadata column {}.'.format(cond, obs_col))

		# use calculated values already in the SwanGraph
		if obs_col == 'dataset' and not self.sc:
			df = pd.DataFrame(data=adata.layers['pi'].todense(),
								 index=adata.obs.index,
								 columns=adata.var.index)
			sums = calc_total_counts(adata, obs_col)

		# recalculate pi and aggregate counts
		else:
			df, sums = calc_pi(adata, ref_df, obs_col=obs_col)

		df = df.transpose()
		sums = sums.transpose()

		# limit to obs_conditions
		if obs_conditions:
			if len(obs_conditions) != 2:
				raise ValueError('obs_conditions must have exactly two values.')
			df = df[obs_conditions]
			sums = sums[obs_conditions]
		else:
			obs_conditions = adata.obs[obs_col].unique().tolist()

		# add gene id
		df = df.merge(ref_df['gid'], how='left',
			left_index=True, right_index=True)

		# add counts and cumulative counts across the samples
		count_cols = sums.columns
		sums['total_counts'] = sums[count_cols].sum(axis=1)
		df = df.merge(sums, how='left', left_index=True, right_index=True, suffixes=(None, '_counts'))

		# construct tables for each gene and test!
		gids = df.gid.unique().tolist()
		gene_de_df = pd.DataFrame(index=gids, columns=['p_val', 'dpi',
											   'pos_iso_1', 'pos_iso_2',
											   'pos_iso_1_dpi', 'pos_iso_2_dpi',
											   'neg_iso_1', 'neg_iso_2',
											   'neg_iso_1_dpi', 'neg_iso_2_dpi'],
						  data=[[np.nan for i in range(10)] \
								for j in range(len(gids))])
		gene_de_df.index.name = 'gid'

		# TODO - should parallelize this
		if verbose:
			n_genes = len(gids)
			pbar = tqdm(total=n_genes)
			pbar.set_description('Testing for DIE for each gene')

		test_results = pd.DataFrame(columns=['gene', 'test_result'])
		for gene in gids:
			gene_df = df.loc[df.gid==gene]
			gene_df, test_result = get_die_gene_table(gene_df, obs_conditions, rc_thresh)
			data = [[gene, test_result]]
			test_result = pd.DataFrame(data=data, columns=['gene', 'test_result'])
			test_results = pd.concat((test_results, test_result))

			# if the gene is valid for testing, do so
			if isinstance(gene_df, pd.DataFrame):
				entry = test_gene(gene_df, obs_conditions)
			else:
				entry = [np.nan for i in range(len(gene_de_df.columns))]
			try:
				gene_de_df.loc[gene] = entry
			except:
				print(entry)
				print(len(entry))
				print(len(gene_de_df.columns))
				assert 1==0
			if verbose:
				pbar.update(1)

		# remove untestable genes and perform p value
		# Benjamini-Hochberg correction
		gene_de_df.dropna(subset=['p_val'], axis=0, inplace=True)
		p_vals = gene_de_df.p_val.tolist()
		if len(p_vals) > 0:
			_, adj_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')
			gene_de_df['adj_p_val'] = adj_p_vals
		else:
			gene_de_df['adj_p_val'] = np.nan
		gene_de_df.reset_index(inplace=True)

		# add gene names
		gnames = self.t_df[['gname', 'gid']]
		gnames.reset_index(drop=True, inplace=True)
		gnames.drop_duplicates(inplace=True)
		gene_de_df = gene_de_df.merge(gnames, how='left', on='gid')

		# add summary table to anndata
		kind = 'die_{}'.format(kind)
		uns_name = make_uns_key(kind=kind,
								obs_col=obs_col,
								obs_conditions=obs_conditions)

		self.adata.uns[uns_name] = gene_de_df

		return gene_de_df, test_results

	# filter differential isoform expression test results based on
	# both an adjusted p value and dpi cutoff
	def get_die_genes(self, kind='iso',
							obs_col='dataset',
							obs_conditions=None,
							p=0.05,
							dpi=10):
		"""
		Filters differential isoform expression test results based on adj.
		p-value and change in percent isoform usage (dpi).

		Parameters:
			kind (str): {'iso', 'tss', 'tes', 'ic'}
				Default: 'iso'
			obs_col (str): Column name from self.adata.obs table to group on.
				Default: 'dataset'
			obs_conditions (list of str, len 2): Which conditions from obs_col
				to compare? Required if obs_col has more than 2 unique values.
			p (float): Adj. p-value threshold to declare a gene as isoform
				switching / having DIE.
				Default: 0.05
			dpi (float): DPI (in percent) value to threshold genes with DIE
				Default: 10

		Returns:
			test (pandas DataFrame): Summary table of genes that pass
				the significance threshold
		"""

		# check to see if this was even tested
		uns_name = make_uns_key(kind='die',
					obs_col=obs_col,
					obs_conditions=obs_conditions,
					die_kind=kind)
		try:
			test = self.adata.uns[uns_name]
		except:
			raise Exception('Problem accessing DIE test results for {}, {}'.format(obs_col, obs_conditions))

		test = test.loc[(test.adj_p_val<=p)&(test.dpi>=dpi)]

		return test

	def add_multi_groupby(self, groupby):
		"""
		Adds a groupby column that is comprised of multiple other columns. For
		instance, if 'sex' and 'age' are already in the obs table, add an
		additional column that's comprised of sex and age.

		Parameters:
			groupby (list of str): List of column names to turn into a multi
				groupby column

		"""

		# determine what to name column
		col_name = '_'.join(groupby)

		if col_name in self.adata.obs.columns:
			col_name += '_1'
			i = 1
		while col_name in self.adata.obs.columns:
			i += 1
			col_name = col_name[:-1]
			col_name += str(i)

		# adatas = [self.adata, self.tss_adata, \
		# 		  self.tes_adata, self.edge_adata]
		# adatas = [self.adata, self.tss_adata, \
		# 		  self.tes_adata]
		adatas = self.get_adatas()
		for adata in adatas:
			adata.obs[col_name] = ''

		for i, group in enumerate(groupby):
			for adata in adatas:
				if i == 0:
					adata.obs[col_name] = adata.obs[groupby].astype(str)
				else:
					adata.obs[col_name] = adata.obs[col_name] + '_' + adata.obs[group].astype(str)
		return col_name

	def rm_multi_groupby(self, col_name):
		"""
		Removes given col_name from the AnnDatas in the SwanGraph.

		Parameters:
			col_name (str): Column name
		"""

		# after we're done, drop this column
		# adatas = [self.adata, self.tss_adata, \
		# 		  self.tes_adata, self.edge_adata]
		# adatas = [self.adata, self.tss_adata, \
		# 		  self.tes_adata]
		adatas = self.get_adatas()
		for adata in adatas:
			adata.obs.drop(col_name, axis=1, inplace=True)

	##########################################################################
	############################ Saving SwanGraphs ###########################
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

	def get_transcript_abundance(self, prefix=None, kind='counts'):
		"""
		Gets dense transcript expression matrix from the current SwanGraph
		in a DataFrame.

		Parameters:
			prefix (str): Path and filename prefix. Resulting file will
				be saved as prefix_transcript_abundance.tsv
				Default: None (will not save)
			kind (str): Choose "tpm", "counts", or "pi"

		Returns:
			df (pandas DataFrame): Abundance and metadata information about
				each transcript.
		"""

		# get collapsed abundance table from edge_adata
		columns = self.adata.var.index.tolist()
		rows = self.adata.obs.index.tolist()
		if kind == 'counts':
			data = self.adata.layers['counts'].toarray()
		elif kind == 'tpm':
			data = self.adata.layers['tpm'].toarray()
		elif kind == 'pi':
			data = self.adata.layers['pi'].toarray()

		df = pd.DataFrame(index=rows, columns=columns, data=data)
		df = df.transpose()
		df.reset_index(inplace=True)
		df.rename({'index': 'tid'}, axis=1, inplace=True)

		# save file
		if prefix:
			fname = '{}_transcript_abundance.tsv'.format(prefix)
			df.to_csv(fname, sep='\t', index=False)

		return df

	def get_edge_abundance(self, prefix=None, kind='counts'):
		"""
		Gets edge expression from the current SwanGraph in a DataFrame
		complete information about where edge is.

		Parameters:
			prefix (str): Path and filename prefix. Resulting file will
				be saved as prefix_edge_abundance.tsv
				Default: None (will not save)
			kind (str): Choose "tpm" or "counts"

		Returns:
			df (pandas DataFrame): Abundance and metadata information about
				each edge.
		"""

		# add location information to edge_df
		temp = self.edge_df.merge(self.loc_df[['chrom', 'coord']],
					how='left', left_on='v1', right_on='vertex_id')
		temp.rename({'coord': 'start'}, axis=1, inplace=True)
		temp = temp.merge(self.loc_df[['coord']],
					how='left', left_on='v2', right_on='vertex_id')
		temp.rename({'coord': 'stop'}, axis=1, inplace=True)
		temp.drop(['v1', 'v2'], axis=1, inplace=True)

		# get collapsed abundance table from edge_adata
		columns = self.edge_adata.var.index.tolist()
		rows = self.edge_adata.obs.index.tolist()
		if kind == 'counts':
			data = self.edge_adata.layers['counts'].toarray()
		elif kind == 'tpm':
			data = self.edge_adata.layers['tpm'].toarray()

		df = pd.DataFrame(index=rows, columns=columns, data=data)
		df = df.transpose()
		df.reset_index(inplace=True)
		df['index'] = df['index'].astype('int')
		df.set_index('index', inplace=True)

		# merge the info together with the abundance
		df = df.merge(temp, how='left', left_index=True, right_on='edge_id')

		# save file
		if prefix:
			fname = '{}_edge_abundance.tsv'.format(prefix)
			df.to_csv(fname, sep='\t', index=False)

		return df

	def get_tss_abundance(self, prefix=None, kind='counts'):
		"""
		Gets TSS expression from the current SwanGraph in a DataFrame with
		complete information about where TSS is.

		Parameters:
			prefix (str): Path and filename prefix. Resulting file will
				be saved as prefix_tss_abundance.tsv
			kind (str): Choose "tpm" or "counts"

		Returns:
			df (pandas DataFrame): Abundance and metadata information about
				each TSS.
		"""
		df = self.get_end_abundance(prefix, kind, how='tss')
		return df

	def get_tes_abundance(self, prefix=None, kind='counts'):
		"""
		Gets TES expression from the current SwanGraph in a DataFrame
		complete information about where TES is.

		Parameters:
			prefix (str): Path and filename prefix. Resulting file will
				be saved as prefix_tes_abundance.tsv
				Default: None (will not save)
			kind (str): Choose "tpm" or "counts"

		Returns:
			df (pandas DataFrame): Abundance and metadata information about
				each TSS.
		"""
		df = self.get_end_abundance(prefix, kind, how='tes')
		return df

	def get_end_abundance(self, prefix, kind, how='tss'):
		"""
		Gets end expression from the current SwanGraph in TSV format with
		complete information about where end is. Called from get_end_abundance

		Parameters:
			prefix (str): Path and filename prefix. Resulting file will
				be saved as prefix_tes_abundance.tsv
			kind (str): Choose "tpm" or "counts" or "pi"
			how (str): Choose "tss" or "tes"

		Returns:
			df (pandas DataFrame): Abundance and metadata information about
				each end.
		"""

		id_col = f'{how}_id'

		# add location information to end_adata.var
		if how == 'tss':
			adata = self.tss_adata
			temp = self.tss_adata.var.copy(deep=True)
		elif how == 'tes':
			adata = self.tes_adata
			temp = self.tes_adata.var.copy(deep=True)
		temp.reset_index(inplace=True)

		# add location information to end_adata.var
		temp = temp.merge(self.loc_df[['chrom', 'coord']],
					how='left', on='vertex_id')

		# get abundance table from end_adata
		columns = adata.var.index.tolist()
		rows = adata.obs.index.tolist()
		if kind == 'counts':
			data = adata.layers['counts'].toarray()
		elif kind == 'tpm':
			data = adata.layers['tpm'].toarray()
		elif kind == 'pi':
			data = adata.layers['pi'].toarray()

		df = pd.DataFrame(index=rows, columns=columns, data=data)
		df = df.transpose()
		df.reset_index(inplace=True)

		# merge the info together with the abundance
		# df = temp.merge(df, how='right', left_index=True, right_index=True)
		df = temp.merge(df, how='right', left_on=id_col, right_on='index')

		# drop index
		df.drop('index', axis=1, inplace=True)

		# save file
		if prefix:
			fname = '{}_{}_abundance.tsv'.format(prefix, how)
			df.to_csv(fname, sep='\t', index=False)
		return df

	##########################################################################
	############################ Plotting utilities ##########################
	##########################################################################

	def set_plotting_colors(self, cmap=None, default=False):
		"""
		Set plotting colors for SwanGraph and browser models.

		Parameters:
			cmap (dict): Dictionary of metadata value : color (hex code with #
				character or named matplotlib color). Keys should be a subset
				or all of ['tss', 'tes', 'internal', 'exon', 'intron', 'browser']
				Default: None
			default (bool): Whether to revert to default colors
		"""

		def adjust_lightness(color, amount=0.5):
			import colorsys
			import matplotlib.colors as mc
			try:
				c = mc.cnames[color]
			except:
				c = color
			c = colorsys.rgb_to_hls(*mc.to_rgb(c))
			c = colorsys.hls_to_rgb(c[0], .8, c[2])
			c = mc.to_hex(c)
			return c

		# user-defined colors
		if cmap:
			grays = ['exon', 'intron', 'tss', 'tes']
			# if colors are named, get hex code
			for key, item in cmap.items():
				if '#' not in item:
					cmap[key] = mc.cnames[item]
				self.pg.color_dict[key] = item
				if key in grays:
					gray_key = '{}_gray'.format(key)
					gray_color = adjust_lightness(item, 2)
					self.pg.color_dict[gray_key] = gray_color

		# just get the default colors
		if default:
			self.pg.color_dict = self.pg.get_color_dict()

	def set_metadata_colors(self, obs_col, cmap):
		"""
		Set plotting colors for datasets based on a column in the metadata
		table.

		Parameters:
			obs_col (str): Name of metadata column to set colors for
			cmap (dict): Dictionary of metadata value : color (hex code with #
				character or named matplotlib color)
		"""

		cmap = copy.deepcopy(cmap)

		# check if obs_col is even there
		if obs_col not in self.adata.obs.columns.tolist():
			raise Exception('Metadata column {} not found'.format(obs_col))

		# TODO check if all values in cmap are in the obs_col
		# also maybe not my problem

		# map values in order specific to
		self.adata.obs[obs_col] = self.adata.obs[obs_col].astype('category')
		obs_order = list(self.adata.obs_names)
		sample_order = self.adata.obs[obs_col].cat.categories.tolist()
		sample_colors = [cmap[s] for s in sample_order]
		self.adata.uns['{}_colors'.format(obs_col)] = sample_colors

		# if colors are named, get hex code
		for key, item in cmap.items():
			if '#' not in item:
				cmap[key] = mc.cnames[item]

		# also store rgb values in dict for use with gen_report
		for key, item in cmap.items():
			item = item[1:]
			r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
			cmap[key] = (r,g,b)
		self.adata.uns['{}_dict'.format(obs_col)] = cmap

		# also add these to the other adatas
		self.edge_adata.uns['{}_colors'.format(obs_col)] = sample_colors
		self.edge_adata.uns['{}_dict'] = cmap
		self.tss_adata.uns['{}_colors'.format(obs_col)] = sample_colors
		self.tss_adata.uns['{}_dict'] = cmap
		self.tes_adata.uns['{}_colors'.format(obs_col)] = sample_colors
		self.tes_adata.uns['{}_dict'] = cmap
		self.gene_adata.uns['{}_colors'.format(obs_col)] = sample_colors
		self.gene_adata.uns['{}_dict'] = cmap

	def plot_graph(self, gid,
				   indicate_dataset=False,
				   indicate_novel=False,
				   prefix=None):
				   # display=True):
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
				Default: True
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

		# # display the plot if option is given
		# if display:
		# 	plt.show()

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

	def plot_browser(self, tid, **kwargs):
		"""
		Plot browser representation for a given transcript

		Parameters:
			tid (str): Transcript ID of transcript to plot

		Returns:
			ax (matplotlib axes): Axes with transcript plotted
		"""

		# init plot settings
		self.check_transcript(tid)
		self.pg.init_plot_settings(self, tid=tid,
								   indicate_dataset=False,
								   indicate_novel=False,
								   browser=True)

		# plot
		ax = self.pg.plot_browser(**kwargs)

		return ax

	def plot_transcript_path(self, tid,
							 indicate_dataset=False,
							 indicate_novel=False,
							 browser=False,
							 prefix=None):
							 # display=True):
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
				Default: True
		"""

		self.check_plotting_args(indicate_dataset, indicate_novel, browser)
		self.check_transcript(tid)

		# reinit PlottedGraph object and plot
		self.pg.init_plot_settings(self, tid=tid,
			indicate_dataset=indicate_dataset,
			indicate_novel=indicate_novel,
			browser=browser)
		self.pg.plot_graph()

		# # display the plot if option is given
		# if display:
		# 	plt.show()

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
	def gen_report(self,
				   gid,
				   prefix,
				   datasets=None,
				   groupby=None,
				   metadata_cols=None,
				   novelty=False,
   				   layer='tpm', # choose from tpm, pi
				   cmap='Spectral_r',
				   include_qvals=False,
				   q=0.05,
				   log2fc=1,
				   qval_obs_col=None,
				   qval_obs_conditions=None,
				   include_unexpressed=False,
				   indicate_novel=False,
				   display_numbers=False,
				   transcript_col='tid',
				   browser=False,
				   order='expression'):
		"""
		Generates a PDF report for a given gene or list of genes according
		to the user's input.

		Parameters:
			gid (str): Gene id or name to generate
				reports for
			prefix (str): Path and/or filename prefix to save PDF and
				images used to generate the PDF
			datasets (dict of lists): Dictionary of {'metadata_col':
				['metadata_category_1', 'metadata_category_2'...]} to represent
				datasets and their order to include in the report.
				Default: Include columns for all datasets / groupby category
			groupby (str): Column in self.adata.obs to group expression
				values by
				Default: None
			metadata_cols (list of str): Columns from metadata tables to include
				as colored bars. Requires that colors have been set using
				set_metadata_colors
			novelty (bool): Include a column to dipslay novelty type of
				each transcript. Requires that a TALON GTF or DB has
				been used to load data in
				Default: False
			layer (str): Layer to plot expression from. Choose 'tpm' or 'pi'
			cmap (str): Matplotlib color map to display heatmap values
				in.
				Default: 'Spectral_r'
			include_qvals (bool): Display q-val of differential expression
				for each transcript and bold entries found to be
				differentially expressed. Requires that de_transcript_test
				has been run, and that abundance information has been
				added to the SwanGraph
				Default: False
			q (float): Q-value significance threshold to use when
				bolding transcripts if include_qvals = True.
				Default: 0.05
			log2fc (float): Log2fc significance threshold to use when
				bolding transcripts if include_qvals = True
			qval_obs_col (str): Metadata column from self.adata
			include_unexpressed (bool): Add transcript entries to report
				that are not expressed in any input dataset.
				Default: False
			indicate_novel (bool): Emphasize novel nodes and edges by
				outlining them and dashing them respectively
				Incompatible with indicate_dataset
				Default: False
			transcript_col (str): Name of column in sg.t_df to use as
				the transcript's display name
				Default: tid
			browser (bool): Plot transcript models in genome browser-
				style format. Incompatible with indicate_dataset and
				indicate_novel
			display_numbers (bool): Display TPM or pi values atop each cell
				Default: False
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
		"""
		# check if groupby column is present
		multi_groupby = False
		indicate_dataset = False
		if groupby:
			# grouping by more than one column
			if type(groupby) == list and len(groupby) > 1:
				for g in groupby:
					if g not in self.adata.obs.columns.tolist():
						raise Exception('Groupby column {} not found'.format(g))
				groupby = self.add_multi_groupby(groupby)
				multi_groupby = True
			elif groupby not in self.adata.obs.columns.tolist():
				raise Exception('Groupby column {} not found'.format(groupby))


		# check if metadata columns are present
		if metadata_cols:
			for c in metadata_cols:
				if c not in self.adata.obs.columns.tolist():
					raise Exception('Metadata column {} not found'.format(c))

				# if we're grouping by a certain variable, make sure
				# the other metadata cols we plan on plotting have unique
				# mappings to the other columns. if just grouping by dataset,
				# since each dataset is unique, that's ok
				if groupby and groupby != 'dataset':
					if groupby == c:
						continue

					temp = self.adata.obs[[groupby, c, 'dataset']].copy(deep=True)
					temp = temp.groupby([groupby, c]).count().reset_index()
					temp = temp.loc[temp.dataset!=0]
					temp = temp.loc[~temp.dataset.isnull()]

					# if there are duplicates from the metadata column, throw exception
					if temp[groupby].duplicated().any():
							raise Exception('Metadata column {} '.format(c)+\
								'not compatible with groupby column {}. '.format(groupby)+\
								'Groupby column has more than 1 unique possible '+\
								'value from metadata column.')

		# check to see if input gene is in the graph
		if gid not in self.t_df.gid.tolist():
			gid = self.get_gid_from_gname(gid)
		self.check_gene(gid)

		# check to see if these plotting settings will play together
		self.check_plotting_args(indicate_dataset,
			indicate_novel, browser)

		# check to see if transcript column is present in sg.t_df
		if transcript_col not in self.t_df.columns:
			raise Exception('Transcript identifier column {} '.format(transcript_col)+\
							'not present in SwanGraph.')

		# get the list of columns to include from the input datasets dict
		if datasets:
			# get a df that is subset of metadata
			# also sort the datasets based on the order they appear in "datasets"
			i = 0
			sorters = []
			for meta_col, meta_cats in datasets.items():
				if meta_col not in self.adata.obs.columns.tolist():
					raise Exception('Metadata column {} not found'.format(meta_col))
				if type(meta_cats) == str:
					meta_cats = [meta_cats]
				if i == 0:
					temp = self.adata.obs.loc[self.adata.obs[meta_col].isin(meta_cats)]
				else:
					temp = temp.loc[temp[meta_col].isin(meta_cats)]
				sort_ind = dict(zip(meta_cats, range(len(meta_cats))))
				sort_col = '{}_sort'.format(meta_col)
				temp[sort_col] = temp[meta_col].map(sort_ind).astype(int)
				sorters.append(sort_col)
				i += 1

			# sort the df based on the order that different categories appear in "datasets"
			temp.sort_values(by=sorters, inplace=True, ascending=True)
			temp.drop(sorters, axis=1, inplace=True)
			columns = temp.dataset.tolist()
			del temp
		else:
			columns = None

		# if we've asked for novelty first check to make sure it's there
		if novelty:
			if not self.has_novelty():
				raise Exception('No novelty information present in the graph. '
					'Add it or do not use the "novelty" report option.')

		# abundance info to calculate TPM on - subset on datasets that will
		# be included
		if columns or datasets:
			subset_adata = self.subset_on_gene_sg(datasets=columns).adata
		else:
			subset_adata = self.adata

		# small SwanGraph with only this gene's data
		sg = self.subset_on_gene_sg(gid=gid, datasets=columns)

		# if we're grouping data, calculate those new numbers
		# additionally order transcripts
		if groupby:
			if layer == 'tpm':
				# use whole adata to calc tpm
				t_df = tpm_df = calc_tpm(subset_adata, obs_col=groupby).transpose()
			elif layer == 'pi':
				# calc tpm just so we can order based on exp
				tpm_df = calc_tpm(subset_adata, obs_col=groupby).transpose()
				t_df, _ = calc_pi(sg.adata, sg.t_df, obs_col=groupby)
				t_df = t_df.transpose()

		else:
			if layer == 'tpm':
				# use whole adata to calc tpm
				t_df = tpm_df = self.get_tpm().transpose()
				t_df = t_df[subset_adata.obs.dataset.tolist()]
			elif layer == 'pi':
				# calc tpm just so we can order based on exp
				t_df = tpm_df = self.get_tpm().transpose()
				t_df = t_df[subset_adata.obs.dataset.tolist()]
				t_df, _ = calc_pi(sg.adata, sg.t_df)
				t_df = t_df.transpose()

		# order transcripts by user's preferences
		if order == 'pi' and not layer == 'pi':
			order = 'expression'

		if order == 'expression' and self.abundance == False:
			order = 'tid'
		elif order == 'expression':
			order = 'log2tpm'
		tids = self.t_df.loc[self.t_df.gid == gid].index.tolist()
		tids = list(set(tids)&set(tpm_df.index.tolist()))
		tpm_df = tpm_df.loc[tids]
		t_df = t_df.loc[tids]
		if order == 'log2tpm':
			order_df = tpm_df
		elif order == 'pi':
			order_df = t_df
		else:
			order_df = t_df
		_, tids = sg.order_transcripts_subset(order_df, order=order)
		tpm_df = tpm_df.loc[tids]
		t_df = t_df.loc[tids]
	#	 print('finished ordering transcripts')
		del tpm_df

		# remove unexpressed transcripts if desired
		if not include_unexpressed:
			t_df = t_df.loc[t_df.any(axis=1)]

		# make sure de has been run if needed
		if include_qvals:
			uns_key = make_uns_key(kind='det',
								   obs_col=qval_obs_col,
								   obs_conditions=qval_obs_conditions)
			qval_df = self.adata.uns[uns_key].copy(deep=True)
			qval_df['significant'] = (qval_df.qval <= q)&(qval_df.log2fc >= log2fc)
		else:
			qval_df = None

		# get tids in this report
		report_tids = t_df.index.tolist()

		# plot each transcript with these settings
		print()
		print('Plotting transcripts for {}'.format(gid))
		sg.plot_each_transcript(report_tids, prefix,
								  indicate_dataset,
								  indicate_novel,
								  browser=browser)
	#	 print('finished plotting each transcript')

		# get a different prefix for saving colorbars and scales
		gid_prefix = prefix+'_{}'.format(gid)

		# if we're plotting tracks, we need a scale as well
		# also set what type of report this will be, 'swan' or 'browser'
		if browser:
			sg.pg.plot_browser_scale()
			save_fig(gid_prefix+'_browser_scale.png')
			report_type = 'browser'
		else:
			report_type = 'swan'

		# plot colorbar for either tpm or pi
		if layer == 'tpm':

			# take log2(tpm) (add pseudocounts)
			t_df = np.log2(t_df+1)

			# min and max tpm vals
			g_max = t_df.max().max()
			g_min = t_df.min().min()

			# create a colorbar
			plt.rcParams.update({'font.size': 30})
			fig, ax = plt.subplots(figsize=(14, 1.5))
			fig.subplots_adjust(bottom=0.5)
			fig.patch.set_visible(False)
			ax.patch.set_visible(False)

			try:
				cmap = plt.get_cmap(cmap)
			except:
				raise ValueError('Colormap {} not found'.format(cmap))

			norm = mpl.colors.Normalize(vmin=g_min, vmax=g_max)

			cb = mpl.colorbar.ColorbarBase(ax,
								cmap=cmap,
								norm=norm,
								orientation='horizontal')
			cb.set_label('log2(TPM)')
			plt.savefig(gid_prefix+'_colorbar_scale.png', format='png',
				bbox_inches='tight', dpi=200)
			plt.clf()
			plt.close()

		elif layer == 'pi':

			# min and max pi vals
			g_max = 100
			g_min = 0

			# create a colorbar between 0 and 1
			plt.rcParams.update({'font.size': 30})
			fig, ax = plt.subplots(figsize=(14, 1.5))
			fig.subplots_adjust(bottom=0.5)
			fig.patch.set_visible(False)
			ax.patch.set_visible(False)

			try:
				cmap = plt.get_cmap(cmap)
			except:
				raise ValueError('Colormap {} not found'.format(cmap))

			norm = mpl.colors.Normalize(vmin=0, vmax=100)

			cb = mpl.colorbar.ColorbarBase(ax,
								cmap=cmap,
								norm=norm,
								orientation='horizontal')
			cb.set_label('Percent of isoform use (' +'$\pi$'+')')
			plt.savefig(gid_prefix+'_colorbar_scale.png', format='png',
				bbox_inches='tight', dpi=200)
			plt.clf()
			plt.close()

		# merge with sg.t_df to get additional columns
		datasets = t_df.columns
		cols = ['novelty', transcript_col]
		t_df = t_df.merge(sg.t_df[cols], how='left', left_index=True, right_index=True)

		# create report
		print('Generating report for {}'.format(gid))
		pdf_name = create_fname(prefix,
					 indicate_dataset,
					 indicate_novel,
					 browser,
					 ftype='report',
					 gid=gid)
		if transcript_col == 'tid':
			t_disp = 'Transcript ID'
		else:
			t_disp = 'Transcript Name'

		report = Report(gid_prefix,
						report_type,
						sg.adata.obs,
						sg.adata.uns,
						datasets=datasets,
						groupby=groupby,
						metadata_cols=metadata_cols,
						novelty=novelty,
						layer=layer,
						cmap=cmap,
						g_min=g_min,
						g_max=g_max,
						include_qvals=include_qvals,
						qval_df=qval_df,
						display_numbers=display_numbers,
						t_disp=t_disp)
		report.add_page()

		# loop through each transcript and add it to the report

		for ind, entry in t_df.iterrows():
			tid = ind
			# display name for transcript
			if transcript_col == 'tid':
				t_disp = tid
			else:
				t_disp = entry[transcript_col]
				# t_disp = entry['tname']
			# else:
			# 	t_disp = tid
			fname = create_fname(prefix,
								 indicate_dataset,
								 indicate_novel,
								 browser,
								 ftype='path',
								 tid=tid)
			report.add_transcript(entry, fname, t_disp)
		report.write_pdf(pdf_name)

		# remove multi groupby column if necessary
		if multi_groupby:
			self.rm_multi_groupby(groupby)

##########################################################################
############################# Data retrieval #############################
##########################################################################
	def get_tpm(self, kind='iso'):
		"""
		Retrieve TPM per dataset.

		Parameters:
			kind (str): {'iso', 'edge', 'tss', 'tes', 'ic'}
				Default: 'iso'

		Returns:
			df (pandas DataFrame): Pandas dataframe where rows are the different
				conditions from `dataset` and the columns are ids in the
				SwanGraph, and values represent the TPM value per
				isoform/edge/tss/tes/ic per dataset.
		"""
		if kind == 'iso':
			adata = self.adata
		elif kind == 'edge':
			adata = self.edge_adata
		elif kind == 'tss':
			adata = self.tss_adata
		elif kind == 'tes':
			adata = self.tes_adata
		elif kind == 'ic':
			adata = self.ic_adata

		adata.X = adata.layers['tpm']
		df = pd.DataFrame.sparse.from_spmatrix(data=adata.X,
			index=adata.obs['dataset'].tolist(),
			columns=adata.var.index.tolist())
		return df

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
		if indicate_novel and not self.annotation:
			raise Exception('Annotation data not present in graph. Use '
							'add_annotation before using indicate_novel')
		if indicate_dataset and indicate_dataset not in self.datasets:
			raise Exception('Dataset {} not present in the graph. '
							''.format(indicate_dataset))

		# if browser, can't do indicate_novel, or indicate_dataset
		if browser:
			if indicate_novel or indicate_dataset:
				raise Exception('Cannot indicate_novel or indicate_dataset '
								'with browser option.')
