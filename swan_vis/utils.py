import networkx as nx
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict

# creates the duplicate index
def create_dupe_index(df, ind_name):
	df[ind_name+'_back'] = df[ind_name]
	return df

# renames old index dupe column in df and resets the index
def reset_dupe_index(df, ind_name):
	df.rename({ind_name: ind_name+'_back'}, inplace=True, axis=1)
	df.reset_index(inplace=True)
	return(df)

# set index, rename dupe index in df
def set_dupe_index(df, ind_name):
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

# return a table indexed by transcript id with the appropriate
# abundance
def process_abundance_file(file, cols, tid_col):

	if type(cols) != list: cols = [cols]

	df = pd.read_csv(file, sep='\t')

	# make sure that tid_col is even in the table
	if tid_col not in df.columns:
		raise Exception('Column {} not found in abundance file.'.format(tid_col))

	# make sure that cols are also in the table
	for col in cols:
		if col not in df.columns:
			raise Exception('Dataset column {} not found in abundance file.'.format(col))

	keep_cols = [tid_col]+cols
	df = df[keep_cols]

	# get the counts
	df['counts'] = df[cols].sum(axis=1)

	# get tpms
	for col in cols:
		total_counts = df[col].sum()
		df['{}_tpm'.format(col)] = (df[col]*1000000)/total_counts
	tpm_cols = ['{}_tpm'.format(col) for col in cols]
	df['tpm'] = df[tpm_cols].mean(axis=1)

	# set up for merging
	cols += tpm_cols
	df.drop(cols, axis=1, inplace=True)
	df.rename({tid_col: 'tid'}, inplace=True, axis=1)

	return df

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

	# fields = df.loc[df.entry_type=='exon', 'fields'].tolist()[0]
	# missing_field = False
	# missing_fields = []
	# if not get_field_value('gene_id', fields):
	# 	missing_field = True
	# 	missing_fields.append('gene_id')
	# if not get_field_value('transcript_id', fields):
	# 	missing_field = True
	# 	missing_fields.append('transcript_id')
	# if missing_field:
	# 	raise Exception('Last column of GTF is missing entry types {}'.format(missing_fields))

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

################################################################################
########################### Analysis-related ###################################
################################################################################

# adata: adata with TSS or iso expression
# conditions: len 2 list of strings of conditions to compare
# col: string, which column the condition labels are in
# how: 'tss' or 'iso'
def get_die(adata, conditions, how='tss', rc=15):

    if how == 'tss':
        id_col = 'tss_id'
    elif how == 'iso':
        id_col = 'tid'

    # make df that we can groupby
    col = 'condition'
    colnames = adata.var[id_col].tolist()
    rownames = adata.obs.dataset.tolist()
    raw = adata.X
    df = pd.DataFrame(data=raw, index=rownames, columns=colnames)
    df.reset_index(inplace=True)
    df.rename({'index':'dataset'}, axis=1, inplace=True)
    samp = adata.obs[['dataset', col]]
    df = df.merge(samp, how='left', on='dataset')

    # limit to only the samples that we want in this condition
#     df[col] = df[col].astype('str')
    df = df.loc[df[col].isin(conditions)]

    # groupby sample type and sum over gen
    df.drop('dataset', axis=1, inplace=True)
    df = df.groupby(col).sum().reset_index()

    # melty df
    var_cols = df.columns.tolist()[1:]
    df = df.melt(id_vars=col, value_vars=var_cols)

    # rename some cols
    df.rename({'variable':id_col,'value':'counts'}, axis=1, inplace=True)

    # merge with gene names
    df = df.merge(adata.var, how='left', on=id_col)

#     # get total number of tss or iso / gene
#     bop = df[['gid', id_col]].groupby('gid').count().reset_index()

    # construct tables for each gene and test!
    gids = df.gid.unique().tolist()
    gene_de_df = pd.DataFrame(index=gids, columns=['p_val', 'dpi'], data=[[np.nan for i in range(2)] for j in range(len(gids))])
    for gene in gids:
        gene_df = df.loc[df.gid==gene]
        p, dpi = test_gene(gene_df, conditions, col, id_col, rc=rc)
        gene_de_df.loc[gene, 'p_val'] = p
        gene_de_df.loc[gene, 'dpi'] = dpi

    # correct p values
    gene_de_df.dropna(axis=0, inplace=True)
    p_vals = gene_de_df.p_val.tolist()
    _, adj_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')
    gene_de_df['adj_p_val'] = adj_p_vals

    gene_de_df.reset_index(inplace=True)

    return gene_de_df

# gene_df: pandas dataframe with expression values in each condition for
# each TSS or isoform in a gene
# conditions: list of str of condition names
# rc: threshold of read count per gene in each condition necessary to test
def test_gene(gene_df, conditions, col, id_col, rc=10):

	gene_df = gene_df.pivot(index=col, columns=id_col, values='counts')
	gene_df = gene_df.transpose()

	groups = gene_df.columns.tolist()
	gene_df['total_counts'] = gene_df[groups].sum(axis=1)
	gene_df.sort_values(by='total_counts', ascending=False, inplace=True)

	# limit to just isoforms with > 0 expression in at least one condition
	cond1 = conditions[0]
	cond2 = conditions[1]
	gene_df = gene_df.loc[(gene_df[cond1]>0)|(gene_df[cond2]>0)]

	# limit to genes with more than 1 isoform expressed
	if len(gene_df.index) <= 1:
		return np.nan, np.nan

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
		if gene_df[cond].sum() < rc:
			return np.nan, np.nan

	# only do the rest if there's nothing left
	if gene_df.empty:
		return np.nan, np.nan

	# calculate the percent of each sample each TSS accounts for
	# TODO: replace with new calc_pi function in swangraph.py
	cond_pis = []
	for cond in conditions:
		total_col = '{}_total'.format(cond)
		pi_col = '{}_pi'.format(cond)
		total_count = gene_df[cond].sum()

		cond_pis.append(pi_col)

		gene_df[total_col] = total_count
		gene_df[pi_col] = (gene_df[cond]/gene_df[total_col])*100

	# compute isoform-level and gene-level delta pis
	gene_df['dpi'] = gene_df[cond_pis[0]] - gene_df[cond_pis[1]]
	gene_df['abs_dpi'] = gene_df.dpi.abs()
	gene_dpi = gene_df.iloc[:2].abs_dpi.sum()

	# chi squared test
	chi_table = gene_df[conditions].to_numpy()
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
