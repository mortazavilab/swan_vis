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
from Graph import Graph

class PlottedGraph(Graph):

	def __init__(self, sg, combine, indicate_dataset, indicate_annotated):

		# save settings so we know how much we need to recompute
		self.combine = combine
		self.indicate_dataset = indicate_dataset 
		self.indicate_annotated = indicate_annotated

		# assign dfs for plotted graph
		self.G = sg.G.copy()
		self.loc_df = sg.loc_df.copy(deep=True)
		self.edge_df = sg.edge_df.copy(deep=True)
		self.t_df = sg.t_df.copy(deep=True)

		# have to deepcopy lists from t_df separately because of how pandas
		# handles lists
		self.t_df['path'] = self.t_df.apply(
			lambda x: copy.deepcopy(x.path), axis=1)

		# used for getting plotting settings 
		self.loc_df['combined'] = False

		# combine non-branching paths
		if combine:
			self.loc_df['sub_color'] = np.nan
			nbps = self.find_nb_paths()
			self.agg_nb_nodes(nbps)
			self.update_ids()

			print()
			print('loc df with reordered ids')
			print(self.loc_df)

		# get positions/sizes of nodes, edges, and labels
		self.calc_pos_sizes()

		# get color/shape plotting settings for nodes, edges
		self.get_plt_settings()

		print(self.loc_df[['color', 'node_shape']])
		print(self.edge_df[['color', 'curve', 'line']])

	# get color/shape plotting settings
	def get_plt_settings(self):

		# plotting styles
		gray = '#DCDCDC'
		yellow = '#F0E442'
		blue = '#0072B2'
		light_blue = '#56B4E9'
		red = '#D55E00'
		orange = '#E69F00'
		pink = '#CC79A7'
		green = '#009E73'

		color_dict = {'intron': pink, 
					  'exon': green,
					  'TSS': blue,
					  'alt_TSS': light_blue,
					  'TES': red,
					  'alt_TES': orange,
					  'internal': yellow}

		# node plotting settings: color
		self.loc_df['color'] = self.loc_df.apply(
			lambda x: get_node_color(x, color_dict), axis=1)

		# node plotting settings: shape
		self.loc_df['node_shape'] = self.loc_df.apply(
			lambda x: get_node_shape(x), axis=1)

		# edge plotting settings: color
		self.edge_df['color'] = self.edge_df.apply(
			lambda x: color_dict[x.edge_type], axis=1)


		edge_dict = {'pos': 'arc3,rad={}'.format(self.rad),
					 'neg': 'arc3,rad=-{}'.format(self.rad),
					 'straight': None}
		print(edge_dict)
		ordered_nodes = self.loc_df.vertex_id.tolist()
		ordered_edges = [(n,m) for n,m in zip(ordered_nodes[:-1],ordered_nodes[1:])]

		# edge plotting settings: curve type
		self.edge_df['raw_index'] = [i for i in range(len(self.edge_df.index))]
		self.edge_df['curve'] = self.edge_df.apply(
			lambda x: get_edge_style(x, ordered_edges, edge_dict), axis=1)
		self.edge_df.drop('raw_index', axis=1, inplace=True)

		# edge plotting settings: line type
		self.edge_df['line'] = self.edge_df.apply(
			lambda x: get_edge_line(x), axis=1)
		
	# calculates the positions and sizes of edges/nodes based on the 
	# number of nodes in the graph
	def calc_pos_sizes(self):

		pos = defaultdict()

		# order nodes based on genomic coords
		ordered_nodes = self.loc_df['vertex_id'].tolist()

		# check if forward or reverse strand 
		# is the first node in the ordered list a TSS or TES? 
		if self.loc_df.loc[ordered_nodes[0], 'TES']:
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
		rad = ((11/540)*x)#+(73/540)
		edge_width = -x/18 + (121/18)

		# assign fields to plotted graph object 
		self.pos = pos
		self.node_size = node_size
		self.label_size = label_size
		self.rad = rad
		self.edge_width = edge_width

	# plots input plotted graph
	def plot_graph(self):

		# plotting stuff
		plt.figure(1, figsize=(14,2.8), frameon=False)
		plt.xlim(-1.05, 1.05)
		plt.ylim(-1.05, 1.05) 

		# plot edges, nodes, and labels
		self.plot_edges()

		# relabel nodes with strings for the purpose of labelling them
		relabel_map = {k:(k if type(k) == str else int(k)) for k in self.G.nodes}
		self.G = nx.relabel_nodes(self.G, relabel_map)
		nx.draw_networkx_labels(self.G,
			self.pos,
			font_size=self.label_size)

		self.plot_nodes()

	# plots edges from edge_df
	def plot_edges(self):
		for _, entry in self.edge_df.iterrows():
			edge = entry.edge_id
			nx.draw_networkx_edges(self.G, self.pos,
				edgelist=[edge],
				width=self.edge_width,
				edge_color=entry.color,
				connectionstyle=entry.curve,
				linestyle=entry.line)

	# plots nodes from loc_df
	def plot_nodes(self):
		for _, entry in self.loc_df.iterrows():
			node = entry.vertex_id
			nx.draw_networkx_nodes(self.G, self.pos,
				nodelist=[node],
				node_color=entry.color,
				node_size=self.node_size,
				node_shape=entry.node_shape)

		# TODO sub node sizes

	###############################################################################
	####################### Combining non-branching paths #########################
	###############################################################################

	# starting nodes:
	# 1. have only one outgoing edge AND
	# 2. have more than one incoming edge | is a TSS | parent is TES AND
	# 3. node is not already in a nbp AND
	# 4. node is not a TES

	def is_nbp_start(self,n,nbps):
		nbp_nodes = [n for p in nbps for n in p]
		if self.G.out_degree(n) == 1: # 1
			tes_parent = False
			parents = list(self.G.predecessors(n))
			if parents: 
				for parent in parents:
					if self.loc_df.loc[parent, 'TES']:
						tes_parent = True
			if self.G.in_degree(n) >= 1 or self.loc_df.loc[n, 'TSS'] or tes_parent: # 2
				if n not in nbp_nodes and not self.loc_df.loc[n, 'TES']: # 3 & 4
					return True
		return False

	# end nodes:
	# 1. have more than one outgoing edge |
	# 2. have more than one incoming edge |
	# 3. is a TSS |
	# 4. parent node is a TES
	def is_nbp_end(self,n):
		if self.G.out_degree(n) > 1: # 1
			return True
		if self.G.in_degree(n) > 1: # 2
			return True
		if self.loc_df.loc[n, 'TSS']: # 3 
			return True
		if self.G.in_degree(n) == 1: # 4
			parent = list(self.G.predecessors(n))[0]
			if self.loc_df.loc[parent, 'TES']:
				return True
		return False

	# http://rosalind.info/problems/ba3m/ ty pavel 
	# modified to not group up TES/TSSs 
	# see is_nbp_start/end to see complete set of conditions
	def find_nb_paths(self):
		nbps = []
		for v in self.G.nodes:
			if self.is_nbp_start(v,nbps):
				if self.G.out_degree(v) > 0:
					for w in self.G.successors(v):
						nbp = [v]
						while not self.is_nbp_end(w):
							nbp.append(w)
							succ = list(self.G.successors(w))
							if succ: w = succ[0]
							else: break
						if len(nbp) > 1: nbps.append(nbp)
		return nbps

	# aggregate nonbranching nodes and add to graph. remove old nodes
	# update loc_df, edge_df, and t_df to reflect these changes
	def agg_nb_nodes(self, nbps):

		G = self.G
		loc_df = self.loc_df
		edge_df = self.edge_df
		t_df = self.t_df

		# loop through each transcript path and each non-branching path,
		# replace nodes that will be aggregated 
		# need to use this because pandas copy DOES NOT deepcopy lists!
		t_df['new_path'] = t_df.apply(lambda x: copy.deepcopy(x.path), axis=1)
		paths = t_df.new_path.tolist() # this gives me a reference to the paths
								       # so when updating I'm directly modifying t_df
		for path in paths: 
			for combined_index, nbp in enumerate(nbps):

				# get the nodes that need to be deleted from the path
				# in the order they need to be deleted 
				del_nodes = sorted([int(i) for i in list(set(path) & set(nbp))],
					key=lambda x: loc_df.loc[x, 'coord'])

				# remove each node that is in a non-branching path
				for i, n in enumerate(del_nodes):
					if i == 0: 
						insertion_index = path.index(n)
					path.remove(n)
				if not del_nodes: insertion_index = -1

				# add aggregate node tuple instead
				if insertion_index != -1:
					path.insert(insertion_index, 'c{}'.format(combined_index))

		# shuffle around path columns in dataframe
		t_df.drop('path', axis=1, inplace=True)
		t_df.rename({'new_path': 'path'}, axis=1, inplace=True)

		mod_G = nx.DiGraph(G)

		combined_index = 0
		loc_df['agg_path'] = np.nan
		loc_df['combined'] = False
		loc_df['combined_types'] = np.nan

		# loop through each nonbranching path
		for path in nbps: 

			print()
			print(path)

			# some beginning of new loop things we need
			path = tuple([int(n) for n in path])
			start = path[0]
			stop = path[-1]
			combined_node = 'c{}'.format(combined_index)

			# get the colors for each aggregate node
			combined_types = [j[0] for j in sorted([i for i in 
			   [('TSS',loc_df.loc[path,'TSS'].tolist().count(True)),
			   ('TES',loc_df.loc[path,'TES'].tolist().count(True)),
			   ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
			   if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]

			# compose the entries for edge_df for new edges that will be added
			# that include the aggregate node
			edge_attrs = {}
			if mod_G.in_degree(start) > 0: 
				in_nodes = list(mod_G.predecessors(start))
				old_edges = [(v1, start) for v1 in in_nodes]
				edge_attrs.update({(edge_id[0], combined_node): item for edge_id,item in edge_df.loc[old_edges,
					['strand', 'edge_type']].to_dict(orient='index').items()})
			if G.out_degree(stop) > 0: 
				out_nodes = list(mod_G.successors(stop))
				old_edges = [(stop, v2) for v2 in out_nodes]
				edge_attrs.update({(combined_node, edge_id[1]): item for edge_id,item in edge_df.loc[old_edges,
					['strand', 'edge_type']].to_dict(orient='index').items()})

			# compose the entry for loc_df for the aggregate node
			coord = loc_df.loc[start, 'coord']
			chrom = loc_df.loc[start, 'chrom']
			strand = loc_df.loc[start, 'strand']
			tss = alt_tss = tes = alt_tes = internal = False
			for n in path:
				if loc_df.loc[n, 'TSS']: tss = True
				if loc_df.loc[n, 'TES']: tes = True
				if loc_df.loc[n, 'alt_TSS']: alt_tss = True
				if loc_df.loc[n, 'alt_TES']: alt_tes = True
				if loc_df.loc[n, 'internal']: internal = True
			loc_attrs = {'chrom': chrom, 'strand': strand,
				   'vertex_id': combined_node,
				   'combined_types': combined_types,
				   'vertex_id_back': combined_node,
				   'agg_path': list(path), 'combined': True,
				   'TSS': tss, 'TES': tes,
				   'alt_TSS': alt_tss, 'alt_TES': alt_tes,
				   'internal': internal, 'coord': coord}

			# # if we're indicating dataset, which datasets are these nodes present in?
			# if args['indicate_dataset']:
			# 	node_attrs = assign_combined_datasets(mod_G, path, node_attrs)

			# add the new node and edges to graph
			mod_G.add_node(combined_node)
			for edge,item in edge_attrs.items():
				mod_G.add_edge(edge[0], edge[1])

			# also add aggregate node and all associated edges to loc_df, edge_df
			loc_df.reset_index(drop=True, inplace=True)
			edge_df.reset_index(drop=True, inplace=True)
			for edge, edge_attr in edge_attrs.items():
				edge_attr.update({'edge_id': edge, 'strand': strand,
								  'v1': edge[0], 'v2': edge[1]})
				edge_df = edge_df.append(edge_attr, ignore_index=True)
			loc_df = loc_df.append(loc_attrs, ignore_index=True)
			loc_df = create_dupe_index(loc_df, 'vertex_id')
			edge_df = create_dupe_index(edge_df, 'edge_id')
			loc_df = set_dupe_index(loc_df, 'vertex_id')
			edge_df = set_dupe_index(edge_df, 'edge_id')

			# remove all old nodes from the graph
			for n in path:
				mod_G.remove_node(n)

			# increment combined node index
			combined_index += 1

		# remove all old nodes/edges from loc_df and edge_df
		del_nodes = [n for path in nbps for n in path]
		del_edges = []
		del_edges += edge_df.loc[edge_df.v1.isin(del_nodes), 'edge_id'].tolist()
		del_edges += edge_df.loc[edge_df.v2.isin(del_nodes), 'edge_id'].tolist()
		del_edges = list(set(del_edges))

		loc_df = loc_df.loc[~loc_df.index.isin(del_nodes)] 
		edge_df = edge_df.loc[~edge_df.index.isin(del_edges)]

		self.G = mod_G
		self.loc_df = loc_df 
		self.edge_df = edge_df 
		self.t_df = t_df

# get the node color #TODO more settings if a path is given perhaps
def get_node_color(x, color_dict):

	# first check if node is combined
	# if x.combined:
	# 	x.sub_color = 

 	# TES 
	if x.alt_TES: return color_dict['alt_TES']
	elif x.TES: return color_dict['TES']
	# TSS
	if x.alt_TSS: return color_dict['alt_TSS']
	elif x.TSS: return color_dict['TSS']
	# internal
	return color_dict['internal']

# get the shape of the node #TODO more settings obvi
def get_node_shape(x):
	# ie hexagonal 
	# diamond
	# etc

	# hexagonal nodes for combined nodes 
	if x.combined:
		return 'H'
	return 'o'

# get the curve/style of the edge
def get_edge_style(x, ordered_edges, edge_dict):

	# over 20 nodes, all should be curved
	if len(ordered_edges) < 20:
		if x.edge_id in ordered_edges:
			return edge_dict['straight']

	if x.raw_index % 2 == 0:
		return edge_dict['pos']
	else:
		return edge_dict['neg']

# get the style of the line for an edge
# TODO add support for indicate_dataset etc
def get_edge_line(x):
	return None

# 	# add a new edge plotting style field or update an existing one
# 	def update_edge_style(self, e, field, value):
# 		self.edge_style[e].update({field: value})

# # add a new node plotting style field or update an existing one
# def update_node_style(self, n, field, value):
# self.node_style[n].update({field: value})

# # add a new sub node plotting style field or update an existing one
# def update_sub_node_style(self, n, field, value):
# self.sub_node_style[n].update({field: value})

# # def save_fig(self, oname):
# # 	plt.tight_layout()
# # 	plt.savefig(oname, format='png', dpi=200)
# # 	plt.clf()

# # plots a subgraph from a PlottedGraph object
# def plot_subgraph(self, nodelist, args):

# self.G = self.G.subgraph(nodelist)
# plot_graph(self, args)

# # starting nodes:
# # 1. have only one outgoing edge AND
# # 2. have more than one incoming edge | is a TSS | parent is TES AND
# # 3. node is not already in a nbp AND
# # 4. node is not a TES

# def is_nbp_start(G,n,nbps):
# nbp_nodes = [n for p in nbps for n in p]
# if G.out_degree(n) == 1: # 1
# tes_parent = False
# parents = list(G.predecessors(n))
# if parents: 
# for parent in parents:
# 	if G.nodes[parent]['TES']:
# 		tes_parent = True
# if G.in_degree(n) >= 1 or G.nodes[n]['TSS'] or tes_parent: # 2
# if n not in nbp_nodes and not G.nodes[n]['TES']: # 3 & 4
# 	return True
# return False

# # end nodes:
# # 1. have more than one outgoing edge |
# # 2. have more than one incoming edge |
# # 3. is a TSS |
# # 4. parent node is a TES
# def is_nbp_end(G,n):
# if G.out_degree(n) > 1: # 1
# return True
# if G.in_degree(n) > 1: # 2
# return True
# if G.nodes[n]['TSS']: # 3 
# return True
# if G.in_degree(n) == 1: # 4
# parent = list(G.predecessors(n))[0]
# if G.nodes[parent]['TES']:
# return True
# return False

# # http://rosalind.info/problems/ba3m/ ty pavel 
# # modified to not group up TES/TSSs 
# # see is_nbp_start/end to see complete set of conditions
# def find_nb_paths(G):

# nbps = []
# for v in G.nodes:
# if is_nbp_start(G,v,nbps):
# if G.out_degree(v) > 0:
# 	for w in G.successors(v):
# 		nbp = [v]
# 		while not is_nbp_end(G,w):
# 			nbp.append(w)
# 			succ = list(G.successors(w))
# 			if succ: w = succ[0]
# 			else: break
# 		if len(nbp) > 1: nbps.append(nbp)
# return nbps

# # aggregate nonbranching nodes and add to graph. remove old nodes
# # update loc_df, edge_df, and t_df to reflect these changes
# def agg_nb_nodes(G, loc_df, edge_df, t_df, nbps, args):

# # loop through each transcript path and each non-branching path,
# # replace nodes that will be aggregated 
# # need to use this because pandas copy DOES NOT deepcopy lists!
# t_df['new_path'] = t_df.apply(lambda x: copy.deepcopy(x.path), axis=1)
# paths = t_df.new_path.tolist() # apparently this gives me a reference to the paths
# 				       # so when updating I'm directly modifying the t_df
# for path in paths: 
# for combined_index, nbp in enumerate(nbps):

# # get the nodes that need to be deleted from the path
# # in the order they need to be deleted 
# del_nodes = sorted([int(i) for i in list(set(path) & set(nbp))],
# 	key=lambda x: G.nodes[x]['coord'])

# # remove each node that is in a non-branching path
# for i, n in enumerate(del_nodes):
# 	if i == 0: 
# 		insertion_index = path.index(n)
# 	path.remove(n)
# if not del_nodes: insertion_index = -1

# # add aggregate node tuple instead
# if insertion_index != -1:
# 	path.insert(insertion_index, 'c{}'.format(combined_index))

# # shuffle around path columns in dataframe
# t_df.drop('path', axis=1, inplace=True)
# t_df.rename({'new_path': 'path'}, axis=1, inplace=True)

# mod_G = nx.DiGraph(G)

# combined_index = 0
# loc_df['agg_path'] = np.nan
# loc_df['combined'] = False
# loc_df['combined_types'] = np.nan
# edge_df = reset_dupe_index(edge_df, 'edge_id')

# # loop through each nonbranching path
# for path in nbps: 

# path = tuple([int(n) for n in path])
# start = path[0]
# stop = path[-1]
# combined_node = 'c{}'.format(combined_index)

# # get the colors for each aggregate node
# if args['color_alt_nodes']:
# combined_types = [j[0] for j in sorted([i for i in 
#    [('alt_TSS',loc_df.loc[path,'alt_TSS'].tolist().count(True)),
#    ('alt_TES',loc_df.loc[path,'alt_TES'].tolist().count(True)),
#    ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
#    if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]
# else: 
# combined_types = [j[0] for j in sorted([i for i in 
#    [('TSS',loc_df.loc[path,'TSS'].tolist().count(True)),
#    ('TES',loc_df.loc[path,'TES'].tolist().count(True)),
#    ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
#    if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]

# # get all incoming edges to first node and
# # all outgoing edges to last node
# data_edges = mod_G.edges(data=True)
# if mod_G.in_degree(start) > 0: 
# in_nodes = list(mod_G.predecessors(start))
# in_edges = [(v1, combined_node) for v1 in in_nodes]

# index_map = {v1: i for i,v1 in enumerate(in_nodes)}
# in_edge_attrs = [md[1] for md in sorted([(v1,d)
# 				 for v1,v2,d in data_edges if v2==start],
# 				 key=lambda x: index_map[x[0]])]
# else: 
# in_edges = []
# in_edge_attrs = []

# if G.out_degree(stop) > 0: 
# out_nodes = list(mod_G.successors(stop))
# out_edges = [(combined_node, v2) for v2 in out_nodes]

# index_map = {v2: i for i,v2 in enumerate(out_nodes)}
# out_edge_attrs = [md[1] for md in sorted([(v2,d)
# 				  for v1,v2,d in data_edges if v1==stop],
# 				  key=lambda x: index_map[x[0]])]
# else: 
# out_edges = []
# out_edge_attrs = []

# edges = in_edges+out_edges
# edge_attrs = in_edge_attrs+out_edge_attrs

# # get the node makeup of the aggregate node, 
# # and remove the old nodes 
# coord = loc_df.loc[loc_df.vertex_id == start, 'coord'].tolist()[0] # use the first coordinate
# chrom = loc_df.loc[loc_df.vertex_id == start, 'chrom'].tolist()[0]
# strand = loc_df.loc[loc_df.vertex_id == start, 'strand'].tolist()[0]
# tss = alt_tss = tes = alt_tes = internal = False
# # annotated = True
# for n in path:
# if mod_G.nodes[n]['TSS']: tss = True
# if mod_G.nodes[n]['TES']: tes = True
# if mod_G.nodes[n]['alt_TSS']: alt_tss = True
# if mod_G.nodes[n]['alt_TES']: alt_tes = True
# if mod_G.nodes[n]['internal']: internal = True
# # if not mod_G.nodes[n]['annotated']: annotated = False

# # add the aggregate node and all associated edges to modified graph
# node_attrs = {'TSS': tss, 'TES': tes,
# 		  'alt_TSS': alt_tss, 'alt_TES': alt_tes,
# 		  'internal': internal, 'coord': coord}
# 		  # 'annotated': annotated}

# # if we're indicating dataset, which datasets are these nodes present in?
# if args['indicate_dataset']:
# node_attrs = assign_combined_datasets(mod_G, path, node_attrs)

# # add the new node and edges to graph
# mod_G.add_node(combined_node)
# nx.set_node_attributes(mod_G, {combined_node: node_attrs})

# for edge, edge_attr in zip(edges, edge_attrs):
# mod_G.add_edge(edge[0], edge[1])
# nx.set_edge_attributes(mod_G, {edge: edge_attr})

# # also add aggregate node and all associated edges to loc_df, edge_df
# node_attrs.update({'chrom': chrom, 'strand': strand,
# 			   'vertex_id': combined_node,
# 			   'combined_types': combined_types,
# 			   'vertex_id_back': combined_node,
# 			   'agg_path': list(path), 'combined': True})
# for edge, edge_attr in zip(edges, edge_attrs):
# edge_attr.update({'edge_ID': edge, 'strand': strand,
# 				  'v1': edge[0], 'v2': edge[1]})
# edge_df = edge_df.append(edge_attr, ignore_index=True)
# loc_df = loc_df.append(node_attrs, ignore_index=True)

# # remove all old nodes from the graph
# nb_nodes = [n for n in path]
# for n in nb_nodes:
# mod_G.remove_node(n)

# # increment combined node index
# combined_index += 1

# # make sure we put dfs back where we left them
# edge_df = set_dupe_index(edge_df, 'edge_id')

# # finally, label all nodes that are now in the graph with their "combined" 
# # status
# mod_G = label_nodes(mod_G, loc_df, 'combined', 'combined')
# mod_G = label_nodes(mod_G, loc_df, 'combined_types', 'combined_types')
# mod_G = label_nodes(mod_G, loc_df, 'agg_path', 'agg_path')

# return mod_G, loc_df, edge_df, t_df

# # returns a dictionary indexed by edge ids of plotting styles 
# def get_edge_plt_settings(G, ordered_nodes, rad, edge_width, args):

# # plotting styles
# pos_style = 'arc3,rad={}'.format(rad)
# neg_style = 'arc3,rad=-{}'.format(rad)
# straight = None

# # get fields if we're highlighting a certain dataset
# if args['indicate_dataset']:
# d_field = 'dataset_'+args['indicate_dataset']
# d_fields = sg.get_dataset_fields(graph=G)

# edges = list(G.edges)
# ordered_edges = [(n,m) for n,m in zip(ordered_nodes[:-1],ordered_nodes[1:])]
# style_dict = defaultdict()

# # straight edges
# if len(ordered_nodes) > 20:
# straight_edges = []
# else:
# straight_edges = list(set(edges)&set(ordered_edges))

# # create a plotting settings dictionary for each edge
# pos = 1
# neg = -1
# for v1,v2,data in G.edges(data=True): 

# e = (v1,v2)

# e_style_dict = {}
# e_style_dict.update({'width': edge_width})

# # straight edge
# if e in straight_edges:
# e_style_dict.update({'connectionstyle': straight})
# # curved edges, alternate between positive and negative
# else:
# if pos > 0:
# 	e_style_dict.update({'connectionstyle': pos_style})
# 	pos *= -1
# 	neg *= -1
# elif neg > 0: 
# 	e_style_dict.update({'connectionstyle': neg_style})
# 	pos *= -1 
# 	neg *= -1

# # dashed edges if indicate_dataset
# if args['indicate_dataset']:
# if is_unique_to_dataset(data, d_field, d_fields):
# 	e_style_dict['linestyle'] = 'dashed'
# else: e_style_dict['linestyle'] = None
# elif args['indicate_novel']:
# if data['dataset_annotation'] == False:
# 	e_style_dict['linestyle'] = 'dashed'
# else:
# 	e_style_dict['linestyle'] = None

# else: 
# e_style_dict['linestyle'] = None

# style_dict.update({e: e_style_dict})

# return style_dict

# # returns a dictionary indexed by node ids of plotting styles
# def get_node_plt_settings(G, pos, node_size, args):

# # get fields if we're highlighting a certain dataset
# if args['indicate_dataset']:
# d_field = 'dataset_'+args['indicate_dataset']	
# d_fields = sg.get_dataset_fields(graph=G)

# # create a plotting settings dictionary for each node
# node_style = defaultdict()
# sub_node_style = defaultdict()
# for n, data in G.nodes(data=True):

# curr_style = defaultdict()
# sub_curr_style = defaultdict()
# curr_style.update({'size': node_size, 'shape': None})

# # combined nodes
# if args['combine'] and data['combined']:
# curr_style.update({'shape': 'H'})

# # sub node
# if len(data['combined_types']) > 1: 
# 	sub_curr_style = {'size': node_size/2, 'shape': 'H'}

# # node only in the query dataset
# if args['indicate_dataset']:			 
# if is_unique_to_dataset(data, d_field, d_fields):
# 	curr_style.update({'shape': 'D', 'size': node_size-30})
# 	if args['combine'] and data['combined']:
# 		curr_style.update({'shape': 'h', 'size': node_size})

# 		if len(data['combined_types']) > 1: 
# 			sub_curr_style.update({'shape': 'h'})
# # novel node
# if args['indicate_novel']:
# if data['dataset_annotation'] == False:
# 	curr_style.update({'shape': 'D', 'size': node_size-30})
# 	if args['combine'] and data['combined']:
# 		curr_style.update({'shape': 'h', 'size': node_size})

# 		if len(data['combined_types']) > 1: 
# 			sub_curr_style.update({'shape': 'h'})

# node_style.update({n: curr_style})
# if sub_curr_style:
# sub_node_style.update({n: sub_curr_style})

# return node_style, sub_node_style

# # returns true if this edge or node is unique to the input dataset
# def is_unique_to_dataset(data, d_field, d_fields):
# if sum(data[i] for i in d_fields) == 1 and data[d_field]:
# return True
# else:
# return False 

# # what datasets does this combined node belong to?
# def assign_combined_datasets(G, path, node_attrs):
# d_fields = sg.get_dataset_fields(graph=G)
# for field in d_fields:
# data = [d for n,d in G.nodes(data=True) if n in path]
# if all(d[field] == True for d in data):
# node_attrs.update({field: True})
# else:
# node_attrs.update({field: False})
# return node_attrs
