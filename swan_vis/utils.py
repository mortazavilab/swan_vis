import networkx as nx
import numpy as np
import pandas as pd
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
	fields = df.loc[df.entry_type=='exon', 'fields'].tolist()[0]
	missing_field = False
	missing_fields = []
	if not get_field_value('gene_id', fields):
		missing_field = True
		missing_fields.append('gene_id')
	if not get_field_value('gene_name', fields):
		missing_field = True		
		missing_fields.append('gene_name')
	if not get_field_value('transcript_id', fields):
		missing_field = True
		missing_fields.append('transcript_id')
	if missing_field:
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

# reorder the locations in a transcript's path based on
# chromosomal coordinate
# TODO
def reorder_locs(path, strand, locs):
	coords = [locs[i] for i in path]
	path_coords = sorted(zip(path, coords), key=lambda x: x[1])
	path = [i[0] for i in path_coords]
	coords = [i[1][1] for i in path_coords]
	if strand == '-':
		path.reverse()
	return path 

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
		
# saves current figure named oname. clears the figure space so additional
# plotting can be done
def save_fig(oname):
	check_dir_loc(oname)
	plt.axis('off')
	plt.tight_layout()
	plt.savefig(oname, format='png', dpi=200)
	plt.clf()
	plt.close()