import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import sqlite3
from utils import *
from plotting_tools import *
import SpliceGraph as sg

#
class PlottedGraph:
	def __init__(self, sg, args):

		G = sg.G.copy()
		loc_df = sg.loc_df.copy(deep=True)
		edge_df = sg.edge_df.copy(deep=True)
		t_df = sg.t_df.copy(deep=True)

		if args['combine']:
			nbps = find_nb_paths(G)
			G, loc_df, edge_df, t_df = agg_nb_nodes(G, loc_df, edge_df, t_df, nbps, args)

		# get positions/sizes of nodes, edges, and labels
		pos, ordered_nodes, node_size, label_size, rad, edge_width = calc_pos_sizes(G)

		# determine edge plotting settings
		edge_style = get_edge_plt_settings(G, ordered_nodes, rad, edge_width, args)
		node_style, sub_node_style = get_node_plt_settings(G, pos, node_size, args)

		# set all fields
		self.G = G
		self.loc_df = loc_df
		self.edge_df = edge_df
		self.t_df = t_df
		self.pos = pos
		self.ordered_nodes = ordered_nodes
		self.node_size = node_size
		self.label_size = label_size
		self.rad = rad
		self.edge_width = edge_width
		self.edge_style = edge_style
		self.node_style = node_style
		self.sub_node_style = sub_node_style 

	# add a new edge plotting style field or update an existing one
	def update_edge_style(self, e, field, value):
		self.edge_style[e].update({field: value})

	# add a new node plotting style field or update an existing one
	def update_node_style(self, n, field, value):
		self.node_style[n].update({field: value})

	# add a new sub node plotting style field or update an existing one
	def update_sub_node_style(self, n, field, value):
		self.sub_node_style[n].update({field: value})

	# def save_fig(self, oname):
	# 	plt.tight_layout()
	# 	plt.savefig(oname, format='png', dpi=200)
	# 	plt.clf()

	# plots a subgraph from a PlottedGraph object
	def plot_subgraph(self, nodelist, args):

		self.G = self.G.subgraph(nodelist)
		plot_graph(self, args)

# starting nodes:
# 1. have only one outgoing edge AND
# 2. have more than one incoming edge | is a TSS | parent is TES AND
# 3. node is not already in a nbp AND
# 4. node is not a TES

def is_nbp_start(G,n,nbps):
	nbp_nodes = [n for p in nbps for n in p]
	if G.out_degree(n) == 1: # 1
		tes_parent = False
		parents = list(G.predecessors(n))
		if parents: 
			for parent in parents:
				if G.nodes[parent]['TES']:
					tes_parent = True
		if G.in_degree(n) >= 1 or G.nodes[n]['TSS'] or tes_parent: # 2
			if n not in nbp_nodes and not G.nodes[n]['TES']: # 3 & 4
				return True
	return False

# end nodes:
# 1. have more than one outgoing edge |
# 2. have more than one incoming edge |
# 3. is a TSS |
# 4. parent node is a TES
def is_nbp_end(G,n):
	if G.out_degree(n) > 1: # 1
		return True
	if G.in_degree(n) > 1: # 2
		return True
	if G.nodes[n]['TSS']: # 3 
		return True
	if G.in_degree(n) == 1: # 4
		parent = list(G.predecessors(n))[0]
		if G.nodes[parent]['TES']:
			return True
	return False

# http://rosalind.info/problems/ba3m/ ty pavel 
# modified to not group up TES/TSSs 
# see is_nbp_start/end to see complete set of conditions
def find_nb_paths(G):

	nbps = []
	for v in G.nodes:
		if is_nbp_start(G,v,nbps):
			if G.out_degree(v) > 0:
				for w in G.successors(v):
					nbp = [v]
					while not is_nbp_end(G,w):
						nbp.append(w)
						succ = list(G.successors(w))
						if succ: w = succ[0]
						else: break
					if len(nbp) > 1: nbps.append(nbp)
	return nbps

# aggregate nonbranching nodes and add to graph. remove old nodes
# update loc_df, edge_df, and t_df to reflect these changes
def agg_nb_nodes(G, loc_df, edge_df, t_df, nbps, args):

	# loop through each transcript path and each non-branching path,
	# replace nodes that will be aggregated 
	paths = t_df.path.tolist() # apparently this gives me a reference to the paths
							   # so when updating I'm directly modifying the t_df
	for path in paths: 
		for combined_index, nbp in enumerate(nbps):

			# get the nodes that need to be deleted from the path
			# in the order they need to be deleted 
			del_nodes = sorted([int(i) for i in list(set(path) & set(nbp))],
				key=lambda x: G.nodes[x]['coord'])

			# remove each node that is in a non-branching path
			for i, n in enumerate(del_nodes):
				if i == 0: 
					insertion_index = path.index(n)
				path.remove(n)
			if not del_nodes: insertion_index = -1

			# add aggregate node tuple instead
			if insertion_index != -1:
				path.insert(insertion_index, 'c{}'.format(combined_index))

	mod_G = nx.DiGraph(G)

	combined_index = 0
	loc_df['agg_path'] = np.nan
	loc_df['combined'] = False
	loc_df['combined_types'] = np.nan
	edge_df = reset_dupe_index(edge_df, 'edge_id')

	# loop through each nonbranching path
	for path in nbps: 

		path = tuple([int(n) for n in path])
		start = path[0]
		stop = path[-1]
		combined_node = 'c{}'.format(combined_index)

		# get the colors for each aggregate node
		if args['color_alt_nodes']:
			combined_types = [j[0] for j in sorted([i for i in 
			   [('alt_TSS',loc_df.loc[path,'alt_TSS'].tolist().count(True)),
			   ('alt_TES',loc_df.loc[path,'alt_TES'].tolist().count(True)),
			   ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
			   if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]
		else: 
			combined_types = [j[0] for j in sorted([i for i in 
			   [('TSS',loc_df.loc[path,'TSS'].tolist().count(True)),
			   ('TES',loc_df.loc[path,'TES'].tolist().count(True)),
			   ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
			   if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]

		# get all incoming edges to first node and
		# all outgoing edges to last node
		data_edges = mod_G.edges(data=True)
		if mod_G.in_degree(start) > 0: 
			in_nodes = list(mod_G.predecessors(start))
			in_edges = [(v1, combined_node) for v1 in in_nodes]

			index_map = {v1: i for i,v1 in enumerate(in_nodes)}
			in_edge_attrs = [md[1] for md in sorted([(v1,d)
							 for v1,v2,d in data_edges if v2==start],
							 key=lambda x: index_map[x[0]])]
		else: 
			in_edges = []
			in_edge_attrs = []

		if G.out_degree(stop) > 0: 
			out_nodes = list(mod_G.successors(stop))
			out_edges = [(combined_node, v2) for v2 in out_nodes]

			index_map = {v2: i for i,v2 in enumerate(out_nodes)}
			out_edge_attrs = [md[1] for md in sorted([(v2,d)
							  for v1,v2,d in data_edges if v1==stop],
							  key=lambda x: index_map[x[0]])]
		else: 
			out_edges = []
			out_edge_attrs = []

		edges = in_edges+out_edges
		edge_attrs = in_edge_attrs+out_edge_attrs

		# get the node makeup of the aggregate node, 
		# and remove the old nodes 
		coord = loc_df.loc[loc_df.vertex_id == start, 'coord'].tolist()[0] # use the first coordinate
		chrom = loc_df.loc[loc_df.vertex_id == start, 'chrom'].tolist()[0]
		strand = loc_df.loc[loc_df.vertex_id == start, 'strand'].tolist()[0]
		tss = alt_tss = tes = alt_tes = internal = False
		# annotated = True
		for n in path:
			if mod_G.nodes[n]['TSS']: tss = True
			if mod_G.nodes[n]['TES']: tes = True
			if mod_G.nodes[n]['alt_TSS']: alt_tss = True
			if mod_G.nodes[n]['alt_TES']: alt_tes = True
			if mod_G.nodes[n]['internal']: internal = True
			# if not mod_G.nodes[n]['annotated']: annotated = False

		# add the aggregate node and all associated edges to modified graph
		node_attrs = {'TSS': tss, 'TES': tes,
					  'alt_TSS': alt_tss, 'alt_TES': alt_tes,
					  'internal': internal, 'coord': coord}
					  # 'annotated': annotated}

		# if we're indicating dataset, which datasets are these nodes present in?
		if args['indicate_dataset']:
			node_attrs = assign_combined_datasets(mod_G, path, node_attrs)

		mod_G.add_node(combined_node)
		nx.set_node_attributes(mod_G, {combined_node: node_attrs})

		for edge, edge_attr in zip(edges, edge_attrs):
			mod_G.add_edge(edge[0], edge[1])
			nx.set_edge_attributes(mod_G, {edge: edge_attr})

		# also add aggregate node and all associated edges to loc_df, edge_df
		node_attrs.update({'chrom': chrom, 'strand': strand,
						   'vertex_id': combined_node,
						   'combined_types': combined_types,
						   'vertex_id_back': combined_node,
						   'agg_path': list(path), 'combined': True})
		for edge, edge_attr in zip(edges, edge_attrs):
			edge_attr.update({'edge_ID': edge, 'strand': strand,
							  'v1': edge[0], 'v2': edge[1]})
			edge_df = edge_df.append(edge_attr, ignore_index=True)
		loc_df = loc_df.append(node_attrs, ignore_index=True)

		# remove all old nodes from the graph
		nb_nodes = [n for n in path]
		for n in nb_nodes:
			mod_G.remove_node(n)

		# increment combined node index
		combined_index += 1

	edge_df = set_dupe_index(edge_df, 'edge_id')

	# finally, label all nodes that are now in the graph with their "combined" 
	# status
	mod_G = label_nodes(mod_G, loc_df, 'combined', 'combined')
	mod_G = label_nodes(mod_G, loc_df, 'combined_types', 'combined_types')
	mod_G = label_nodes(mod_G, loc_df, 'agg_path', 'agg_path')

	return mod_G, loc_df, edge_df, t_df

# calculates the positions and sizes of edges/nodes based on the 
# number of nodes in the graph
def calc_pos_sizes(G):

	pos = defaultdict()

	# order nodes based on genomic coords
	ordered_nodes = [i[0] for i in sorted([[j, int(n['coord'])]
			for j,n in G.nodes(data=True)], key=lambda x: x[1])]

	# check if forward or reverse strand 
	# is the first node in the ordered list a TSS or TES? 
	if G.nodes[ordered_nodes[0]]['TES']:
		ordered_nodes.reverse()
	
	y_coord = 0
	x_coords = np.linspace(-1, 1, len(ordered_nodes))
	for x_coord, n in zip(x_coords, ordered_nodes): 
		pos[n] = [x_coord, y_coord]

	x = len(ordered_nodes)

	# area-related sizes need to be non-linearly scaled
	# calculated by fitting power series curve to handpicked sizes
	node_size = 19248*(x**-1.14) 
	label_size = 43.9*(x**-0.484)

	# linearly-related sizes
	rad = ((11/540)*x)+(73/540)
	edge_width = -x/18 + (121/18)

	return pos, ordered_nodes, node_size, label_size, rad, edge_width

# returns a dictionary indexed by edge ids of plotting styles 
def get_edge_plt_settings(G, ordered_nodes, rad, edge_width, args):

	# plotting styles
	pos_style = 'arc3,rad={}'.format(rad)
	neg_style = 'arc3,rad=-{}'.format(rad)
	straight = None

	# get fields if we're highlighting a certain dataset
	if args['indicate_dataset']:
		d_field = 'dataset_'+args['indicate_dataset']
		d_fields = sg.get_dataset_fields(G)

	edges = list(G.edges)
	ordered_edges = [(n,m) for n,m in zip(ordered_nodes[:-1],ordered_nodes[1:])]
	style_dict = defaultdict()

	# straight edges
	if len(ordered_nodes) > 20:
		straight_edges = []
	else:
		straight_edges = list(set(edges)&set(ordered_edges))
	
	# create a plotting settings dictionary for each edge
	pos = 1
	neg = -1
	for v1,v2,data in G.edges(data=True): 

		e = (v1,v2)

		e_style_dict = {}
		e_style_dict.update({'width': edge_width})

		# straight edge
		if e in straight_edges:
			e_style_dict.update({'connectionstyle': straight})
		else:
			if pos > 0:
				e_style_dict.update({'connectionstyle': pos_style})
				pos *= -1
				neg *= -1
			elif neg > 0: 
				e_style_dict.update({'connectionstyle': neg_style})
				pos *= -1 
				neg *= -1

		# dashed edges if indicate_dataset
		if args['indicate_dataset']:
			if is_unique_to_dataset(data, d_field, d_fields):
				e_style_dict['linestyle'] = 'dashed'
			else: e_style_dict['linestyle'] = None
		else: 
			e_style_dict['linestyle'] = None

		style_dict.update({e: e_style_dict})

	return style_dict

# returns a dictionary indexed by node ids of plotting styles
def get_node_plt_settings(G, pos, node_size, args):

	# get fields if we're highlighting a certain dataset
	if args['indicate_dataset']:
		d_field = 'dataset_'+args['indicate_dataset']	
		d_fields = sg.get_dataset_fields(G)

	# create a plotting settings dictionary for each node
	node_style = defaultdict()
	sub_node_style = defaultdict()
	for n, data in G.nodes(data=True):

		curr_style = defaultdict()
		sub_curr_style = defaultdict()
		curr_style.update({'size': node_size, 'shape': None})

		# combined nodes
		if args['combine'] and data['combined']:
			curr_style.update({'shape': 'H'})

			# sub node
			if len(data['combined_types']) > 1: 
				sub_curr_style = {'size': node_size/2, 'shape': 'H'}

		# node only in the query dataset
		if args['indicate_dataset']:			 
			if is_unique_to_dataset(data, d_field, d_fields):
				curr_style.update({'shape': 'D', 'size': node_size-30})
				if args['combine'] and data['combined']:
					curr_style.update({'shape': 'h', 'size': node_size})

					if len(data['combined_types']) > 1: 
						sub_curr_style.update({'shape': 'h'})

		node_style.update({n: curr_style})
		if sub_curr_style:
			sub_node_style.update({n: sub_curr_style})

	return node_style, sub_node_style

# returns true if this edge or node is unique to the input dataset
def is_unique_to_dataset(data, d_field, d_fields):
	if sum(data[i] for i in d_fields) == 1 and data[d_field]:
		return True
	else:
		return False 

def assign_combined_datasets(G, path, node_attrs):
	d_fields = sg.get_dataset_fields(G)
	for field in d_fields:
		data = [d for n,d in G.nodes(data=True) if n in path]
		if all(d[field] == True for d in data):
			node_attrs.update({field: True})
		else:
			node_attrs.update({field: False})
	return node_attrs