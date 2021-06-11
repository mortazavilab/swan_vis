import networkx as nx
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
from tqdm import tqdm


pd.options.mode.chained_assignment = None

def create_dupe_index(df, ind_name):
	"""
	Creates a duplicate column in the input DataFrame from the input column name

	Parameters:
		df (pandas DataFrame): DataFrame to create duplicate column
		ind_name (str): Name of column to duplicate
	"""
	df[ind_name+'_back'] = df[ind_name]
	return df

def reset_dupe_index(df, ind_name):
	"""
	Reset index and rename duplicate index column

	Parameters:
		df (pandas DataFrame): DataFrame to reset index
		ind_name (str): Name of column to reset
	"""
	df.rename({ind_name: ind_name+'_back'}, inplace=True, axis=1)
	df.reset_index(inplace=True)
	return(df)

# set index, rename dupe index in df
def set_dupe_index(df, ind_name):
	"""
	Set duplicated column from create_dupe_index as the index and rename the
	duplicated column

	Parameters:
		df (pandas DataFrame): DataFrame to set index
		ind_name (str): Name of column to set as index
	"""
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

# creates a dictionary of the last field of a gtf
# adapted from Dana Wyman
def get_fields(fields):

	attributes = {}

	description = fields.strip()
	description = [x.strip() for x in description.split(";")]
	for pair in description:
		if pair == "": continue

		pair = pair.replace('"', '')
		key, val = pair.split()
		attributes[key] = val

	# put in placeholders for important attributes (such as gene_id) if they
	# are absent
	if 'gene_id' not in attributes:
		attributes['gene_id'] = 'NULL'

	return attributes

# check to see if a file save location is valid
def check_dir_loc(loc):
	if '/' in loc:
		d = '/'.join(loc.split('/')[:-1])
		if not os.path.isdir(d):
			raise Exception('Directory {} is not found. '
				'Try a different save location'.format(d))

# check to see if a file exists
def check_file_loc(loc, ftype):
	if not os.path.isfile(loc):
		raise Exception('{} file not found at {}. '
			'Check path.'.format(ftype, loc))

# creates a file name based on input plotting arguments
def create_fname(prefix, indicate_dataset,
				 indicate_novel, browser,
				 ftype='summary', tid=None, gid=None):
	fname = prefix
	if indicate_dataset:
		fname += '_{}'.format(indicate_dataset)
	if indicate_novel:
		fname += '_novel'
	if browser:
		fname += '_browser'
	if tid:
		fname += '_{}'.format(tid)
	if gid:
		fname += '_{}'.format(gid)
	if ftype == 'summary':
		fname += '_summary.png'
	elif ftype == 'path':
		fname += '_path.png'
	elif ftype == 'report':
		fname += '_report.pdf'
	return fname

# checks if a file is a gtf or a db
def gtf_or_db(fname):
	ext = fname.split('.')[-1]
	if ext == 'gtf': return 'gtf'
	elif ext == 'db': return 'db'
	else:
		raise Exception('File type must be gtf or db. '
			'Type received is {}'.format(ext))

# validate that a gtf has the correct fields and info in it
def validate_gtf(fname):
	df = pd.read_csv(fname, sep='\t', usecols=[2,8], comment='#',
		names=['entry_type', 'fields'])
	if df.empty:
		raise Exception("GTF can't load correctly. Could be a header problem.")

	# first make sure that there are transcript and exon entries
	missing_entry_type = False
	missing_entry_types = []
	entry_types = df.entry_type.unique().tolist()
	if 'transcript' not in entry_types:
		missing_entry_type = True
		missing_entry_types.append('transcript')
	if 'exon' not in entry_types:
		missing_entry_type = True
		missing_entry_types.append('exon')
	if missing_entry_type:
		raise Exception('GTF is missing entry types {}'.format(missing_entry_types))

	# next check if gene_id, transcript_id, and gene_name fields exist
	# in each line of the thing
	df = df.loc[df.entry_type == 'transcript']
	df['tid'] = df.fields.str.extract(r'transcript_id "([A-z]+[0-9.]+)"',
		expand=False)
	df['gid'] = df.fields.str.extract(r'gene_id "([A-z]+[0-9.]+)"',
		expand=False)
	df.drop('fields', axis=1, inplace=True)

	missing_fields = []
	if df.tid.isnull().any():
		missing_fields.append('transcript_id')
	if df.gid.isnull().any():
		missing_fields.append('gene_id')
	if missing_fields:
		raise Exception('Last column of GTF is missing entry types {}'.format(missing_fields))

# depending on the strand, determine the start and stop
# coords of an intron or exon
def find_edge_start_stop(v1, v2, strand):
	if strand == '-':
		start = max([v1, v2])
		stop = min([v1, v2])
	elif strand == '+':
		start = min([v1, v2])
		stop = max([v1, v2])
	return start, stop

# reorder exon ids from create_dfs_gtf
def reorder_exons(exon_ids):
	strand = exon_ids[0].split('_')[-2]
	coords = [int(i.split('_')[-4]) for i in exon_ids]
	exons = sorted(zip(exon_ids, coords), key=lambda x: x[1])
	exons = [i[0] for i in exons]
	if strand == '-':
		exons.reverse()
	return exons

# # reorder the locations in a transcript's path based on
# # chromosomal coordinate
# def reorder_locs(path, strand, locs):
# 	coords = [locs[i] for i in path]
# 	path_coords = sorted(zip(path, coords), key=lambda x: x[1])
# 	path = [i[0] for i in path_coords]
# 	coords = [i[1][1] for i in path_coords]
# 	if strand == '-':
# 		path.reverse()
# 	return path

##########################################################################
############### Related to calculating abundance values ##################
##########################################################################

def calc_total_counts(adata, obs_col='dataset'):
	"""
	Calculate cumulative expression per adata entry based on condition given
	by `obs_col`. Default column to use is `adata.obs` index column, `dataset`.

	Parameters:
		adata (anndata AnnData): Annotated data object from the SwanGraph
		obs_col (str): Column name from adata.obs table to group on.
			Default: 'dataset'

	Returns:
		df (pandas DataFrame): Pandas datafrom where rows are the different
			conditions from `obs_col` and the columns are transcript ids in the
			SwanGraph, and values represent the cumulative counts per isoform
			per condition.

	"""
	adata.X = adata.layers['counts']
	df = pd.DataFrame(data=adata.X, index=adata.obs[obs_col].tolist(), \
		columns=adata.var.index.tolist())

	# add up values on condition (row)
	df = df.groupby(level=0).sum()
	# df = df.transpose()

	return df

def calc_pi(adata, t_df, obs_col='dataset'):
	"""
	Calculate the percent isoform per gene per condition given by `obs_col`.
	Default column to use is `adata.obs` index column, `dataset`.

	Parameters:
		adata (anndata AnnData): Annotated data object from the SwanGraph
		t_df (pandas DataFrame): Pandas Dataframe that has index to
			gene id mapping
		obs_col (str): Column name from adata.obs table to group on.
			Default: 'dataset'

	Returns:
		df (pandas DataFrame): Pandas DataFrame where rows are the different
			conditions from `obs_col` and the columns are transcript ids in the
			SwanGraph, and values represent the percent isoform usage per gene
			per condition.
		sums (pandas DataFrame): Pandas DataFrame where rows are the different
			conditions from `obs_col` and the columns are transcript ids in the
			SwanGraph, and values represent the cumulative counts per isoform
			per condition.
	"""

	# calculate cumulative counts across obs_col
	id_col = adata.var.index.name
	conditions = adata.obs[obs_col].unique().tolist()
	df = calc_total_counts(adata, obs_col=obs_col)
	df = df.transpose()
	sums = df.copy(deep=True)
	sums = sums[conditions]
	sums = sums.transpose()

	# add gid
	df = df.merge(t_df['gid'], how='left', left_index=True, right_index=True)

	# calculate total number of reads per gene per condition
	temp = df.copy(deep=True)
	temp.reset_index(drop=True, inplace=True)
	totals = temp.groupby('gid').sum().reset_index()

	# merge back in
	df.reset_index(inplace=True)
	df.rename({'index':id_col}, axis=1, inplace=True)
	df = df.merge(totals, on='gid', suffixes=(None, '_total'))
	del totals

	# calculate percent iso exp for each gene / transcript / condition
	pi_cols = []
	for c in conditions:
		cond_col = '{}_pi'.format(c)
		total_col = '{}_total'.format(c)
		df[cond_col] = (df[c]/df[total_col])*100
		pi_cols.append(cond_col)

	# cleanup: fill nans with 0, set indices, rename cols
	df.fillna(0, inplace=True)

	# formatting
	df.set_index(id_col, inplace=True)
	df = df[pi_cols]
	for col in pi_cols:
		new_col = col[:-3]
		df.rename({col: new_col}, axis=1, inplace=True)

	# reorder columns like adata.obs
	df = df[adata.obs[obs_col].unique().tolist()]
	df = df.transpose()
	# df.index.name = obs_col # maybe put this back in ?

	# reorder in adata.var / t_df order
	df = df[t_df[id_col].tolist()]

	return df, sums

def calc_tpm(adata, t_df, obs_col='dataset'):
	"""
	Calculate the TPM per condition given by `obs_col`.
	Default column to use is `self.adata.obs` index column, `dataset`.

	Parameters:
		adata (anndata AnnData): Annotated data object from the SwanGraph
		t_df (pandas DataFrame): Pandas Dataframe that has index to
			gene id mapping
		obs_col (str): Column name from adata.obs table to group on.
			Default: 'dataset'

	Returns:
		df (pandas DataFrame): Pandas datafrom where rows are the different
			conditions from `obs_col` and the columns are transcript ids in the
			SwanGraph, and values represent the TPM value per isoform per
			condition.
	"""

	# calculate cumulative counts across obs_col
	id_col = adata.var.index.name
	conditions = adata.obs[obs_col].unique().tolist()
	df = calc_total_counts(adata, obs_col=obs_col)
	df = df.transpose()

	# calculate tpm per isoform per condition
	tpm_cols = []
	for c in conditions:
		cond_col = '{}_tpm'.format(c)
		total_col = '{}_total'.format(c)
		df[total_col] = df[c].sum()
		print()
		print(c)
		print(df[total_col])
		df[cond_col] = (df[c]*1000000)/df[total_col]
		tpm_cols.append(cond_col)

	# formatting
	df.index.name = 'tid'
	df = df[tpm_cols]
	for col in tpm_cols:
		new_col = col[:-4]
		df.rename({col: new_col}, axis=1, inplace=True)

	# reorder columns like adata.obs
	df = df[adata.obs[obs_col].unique().tolist()]
	df = df.transpose()

	# reorder in adata.var / t_df order
	df = df[t_df[id_col].tolist()]

	return df

##########################################################################
####################### Related to file parsing ##########################
##########################################################################

def parse_gtf(gtf_file, verbose):
	"""
	Get the unique transcripts and exons that are present in a GTF
	transcriptome.

	Parameters:
		gtf_file (str): File path of GTF
		verbose (bool): Display progress

	Returns:
		transcripts (dict of dict): Dictionary of transcripts in GTF. Keys are
			transcript ids. Items are a dictionary of gene id, gene name,
			transcript id (same as key), transcript name, strand, and exons
			belonging to the transcript.
		exons (dict of dict): Dictionary of exons in GTF. Keys are exon ids
			which consist of chromosome_v1_v2_strand_exon. Items are a
			dictionary of edge id (same as key), chromosome, v1, v2, strand,
			and edge type (all exon in this case) of each exon.
		from_talon (bool): Whether or not the GTF was determined to be
			from TALON
	"""

	# counts the number of transcripts in a given GTF
	# so that we can track progress
	def count_transcripts(gtf_file):
		df = pd.read_csv(gtf_file, sep='\t', usecols=[2],
			names=['entry_type'], comment='#')
		df = df.loc[df.entry_type == 'transcript']
		n = len(df.index)
		return n

	# get the number of transcripts in the file
	n_transcripts = count_transcripts(gtf_file)

	# dictionaries to hold unique edges and transcripts
	transcripts = {}
	exons = {}
	from_talon = False

	# display progess
	if verbose:
		pbar = tqdm(total=n_transcripts)
		counter = 0

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

				# update progress bar
				if verbose:
					counter+=1
					if counter % 100 == 0:
						pbar.update(100)
						pbar.set_description('Processing transcripts')

				attributes = get_fields(fields)

				# check if this gtf has transcript novelty vals
				# for the first transcript entry
				if not transcripts:
					if 'talon_transcript' in attributes:
						from_talon = True

				tid = attributes['transcript_id']
				gid = attributes['gene_id']

				# check if there's a gene/transcript name field and
				# add one if not
				if 'gene_name' not in attributes:
					attributes['gene_name'] = attributes['gene_id']
				if 'transcript_name' not in attributes:
					attributes['transcript_name'] = attributes['transcript_id']

				gname = attributes['gene_name']
				tname = attributes['transcript_name']

				# add transcript to dictionary
				entry = {'gid': gid,
						 'gname': gname,
						 'tid': tid,
						 'tname': tname,
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
								  'strand': strand,
								  'edge_type': 'exon'}}
					exons.update(edge)

				# add this exon to the transcript's list of exons
				if tid in transcripts:
					transcripts[tid]['exons'].append(eid)

	t_df = pd.DataFrame(transcripts).transpose()
	exon_df = pd.DataFrame(exons).transpose()
	return t_df, exon_df, from_talon

##########################################################################
######################## Related to DIE testing ##########################
##########################################################################

def get_die_gene_table(gene_df, conditions, rc):
	"""
	Creates a length n (max 11) table ready for DIE testing for a given gene.
	Returns None for genes deemed untestable if there aren't enough reads per
	condition or if the gene only has one isoform. Removes isoforms that are
	unexpressed in both conditions. Aggregates the counts for the
	lowest-expressed isoforms (11+). Calculates dpi (change in percent isoform
	usage) between conditions.

	Parameters:
		gene_df (pandas DataFrame): DataFrame of transcript counts and percent
			isoform expression per isoform of one gene.
		conditions (list of str, len 2): Names of 2 conditions being tested
		rc (int): Number of reads needed per gene per condition for testing

	Returns:
		gene_df (pandas DataFrame): :ength n table of counts per isoform per
			condition, percent isoform per gene per condition, and change in
			percent isoform across conditions IF the gene is testable
		gene_df (None): Returns None if the gene was deemed untestable.
	"""

	gene_df.sort_values(by='total_counts', ascending=False, inplace=True)

	# limit to just isoforms with > 0 expression in at least one condition
	cond1 = conditions[0]
	cond2 = conditions[1]
	gene_df = gene_df.loc[(gene_df[cond1]>0)|(gene_df[cond2]>0)]

	# limit to genes with more than 1 isoform expressed
	if len(gene_df.index) <= 1:
		return None

	# if there are more than 11 isoforms, agg. the n - 11 least expressed
	# isoforms into one
	if len(gene_df.index) > 11:
		gene_df.reset_index(inplace=True)

		temp = gene_df.iloc[10:].sum()
		temp[id_col] = 'all_other'
		temp.index.name = None
		temp = pd.DataFrame(temp).transpose()

		gene_df = gene_df.iloc[:10]
		gene_df = pd.concat([gene_df, temp])

	# does this gene reach the desired read count threshold?
	for cond in conditions:
		counts_col = cond+'_counts'
		if gene_df[counts_col].sum() < rc:
			return None

	# only do the rest if there's something left
	if gene_df.empty:
		return None

	# compute isoform-level and gene-level delta pis
	gene_df['dpi'] = gene_df[cond1] - gene_df[cond2]

	return gene_df

def test_gene(gene_df, conditions):
	"""
	Performs a chi-squared test between two conditions on their read counts.
	Also calculates the gene's DPI, or change in percent isoform as the sum of
	either the top two positive changes or top two negative changes (whichever
	is greater in magnitude).

	Parameters:
		gene_df (pandas DataFrame): Output from get_die_gene_table.
		conditions (list of str, len 2): Name of condition columns

	Returns:
		p (float): P-value result of chi-squared test on gene
		dpi (float): Overall change in isoform expression
	"""

	counts_cols = [c+'_counts' for c in conditions]

	# get highest 2 positive dpis
	temp = gene_df.sort_values(by='dpi', ascending=False)
	temp = temp.loc[temp.dpi > 0]

	# if there are fewer than 2 isoforms
	if len(temp.index) >= 2:
		pos_dpi = temp.iloc[:2].dpi.sum(axis=0)
	else:
		pos_dpi = temp.dpi.sum(axis=0)

	# get highest 2 negative dpis
	temp = gene_df.sort_values(by='dpi', ascending=True)
	temp = temp.loc[temp.dpi < 0]

	# if there are fewer than 2 isoforms
	if len(temp.index) >= 2:
		neg_dpi = abs(temp.iloc[:2].dpi.sum(axis=0))
	else:
		neg_dpi = abs(temp.dpi.sum(axis=0))

	gene_dpi = max(pos_dpi, neg_dpi)

	# chi squared test
	chi_table = gene_df[counts_cols].to_numpy()
	chi2, p, dof, exp = st.chi2_contingency(chi_table)

	return p, gene_dpi

# turn a list of dataset groups and names for those groups into a
# dictionary
def make_cond_map(groups, group_names):
	cond_map = dict()
	for group, group_name in zip(groups, group_names):
		if type(group) == list:
			for group_item in group:
				cond_map[group_item] = group_name
		else:
			cond_map[group] = group_name
	return cond_map

# get novelty types associated with each transcript
def get_transcript_novelties(fields):
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

# reformat talon abundance file for the generic format expected by swan
def reformat_talon_abundance(fname, ofile=None):
	check_file_loc(fname, 'TALON abundance')

	df = pd.read_csv(fname, sep='\t')
	drop_cols = ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_gene_name',
		'annot_transcript_name', 'n_exons', 'length', 'gene_novelty',
		'transcript_novelty', 'ISM_subtype']
	df.drop(drop_cols, axis=1, inplace=True)

	# if not writing output file just return df
	if not ofile:
		return df

	# otherwise dump to output file
	df.to_csv(ofile, sep='\t', index=False)

# saves current figure named oname. clears the figure space so additional
# plotting can be done
def save_fig(oname):
	check_dir_loc(oname)
	plt.axis('off')
	plt.tight_layout()
	plt.savefig(oname, format='png', dpi=200)
	plt.clf()
	plt.close()
