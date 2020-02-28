import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import matplotlib.patches as pch
from swan_vis.utils import *
from swan_vis.graph import Graph

class PlottedGraph(Graph):

	def __init__(self,
				 sg,
				 combine,
				 indicate_dataset,
				 indicate_novel,
				 tid=None,
				 browser=False,
				 gid=None):

		# save settings so we know how much we need to recompute
		self.combine = combine
		self.indicate_dataset = indicate_dataset 
		self.indicate_novel = indicate_novel
		# self.path = path
		self.gid = gid # if we're plotting a set of transcripts belonging to a specific gene
		               # or a summary graph of a specific gene
		self.tid = tid
		self.browser = browser

		# other fields to copy over
		self.datasets = sg.datasets

		# assign dfs for plotted graph
		self.G = sg.G.copy()
		self.loc_df = sg.loc_df.copy(deep=True)
		self.edge_df = sg.edge_df.copy(deep=True)
		self.t_df = sg.t_df.copy(deep=True)

		# have to deepcopy lists from t_df separately because of how pandas
		# handles lists
		self.t_df['path'] = self.t_df.apply(
			lambda x: copy.deepcopy(x.path), axis=1)

		if gid: 
			self.subset_on_gene(self.gid)
			self.g_min, self.g_max = self.get_gene_min_max(gid)
		if self.tid: 
			gid = self.get_gid_from_tid(self.tid)
			self.gid = gid
			self.subset_on_gene(gid)
			self.path = self.get_path_from_tid(self.tid)
		else:
			self.path = None

		## TODO is this the best way to handle browser vs. not browser graphs?
		# if we're making a swan graph
		if not browser:

			# used for getting plotting settings 
			self.loc_df['combined'] = False

			# combine non-branching paths
			if combine:
				nbps = self.find_nb_paths()
				self.agg_nb_nodes(nbps)

			# get positions/sizes of nodes, edges, and labels
			self.calc_pos_sizes()

			# get color/shape plotting settings for nodes, edges
			self.get_plt_settings()

		# if we're just making browser track things
		# elif browser:
		# 	1

	###############################################################################
	################### Getting plotting settings for nodes/edges #################
	###############################################################################

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

		gray_light_blue = '#c1ddec'
		gray_yellow = '#f8f6db'
		gray_blue = '#c4e0f0'
		gray_red = '#f1d9c6'
		gray_orange = '#ebdec3'
		gray_pink = '#efdae5'
		gray_green = '#cdf5ea'

		color_dict = {'intron': {'normal': pink, 'gray': gray_pink},
					  'exon': {'normal': green, 'gray': gray_green},
					  'TSS': {'normal': blue, 'gray': gray_blue},
					  'alt_TSS': {'normal': light_blue, 'gray': gray_light_blue},
					  'TES': {'normal': red, 'gray': gray_red},
					  'alt_TES': {'normal': orange, 'gray': gray_orange},
					  'internal': {'normal': yellow, 'gray': gray_yellow}}

		# get the list of datasets we should be looking at to determine
		# node or edge uniqueness
		dataset_cols = self.get_dataset_cols()

		# TODO use this when doing unique_to_dataset
		# if self.indicate_dataset:
		# 	dataset_cols.remove(self.indicate_dataset)

		###############
		#### NODES ####
		###############

		# node plotting settings: color
		self.loc_df = self.loc_df.apply(
			lambda x: self.get_node_color(x, color_dict), axis=1)

		# node plotting settings: shape
		self.loc_df['node_shape'] = self.loc_df.apply(
			lambda x: self.get_node_shape(x, dataset_cols), axis=1)

		###############
		#### EDGES ####
		###############

		# edge plotting settings: color
		self.edge_df['color'] = self.edge_df.apply(
			lambda x: self.get_edge_color(x, color_dict), axis=1)

		edge_dict = {'pos': 'arc3,rad=',
					 'neg': 'arc3,rad=-',
					 'straight': None}
		ordered_nodes = self.get_ordered_nodes()
		ordered_edges = [(n,m) for n,m in zip(ordered_nodes[:-1],ordered_nodes[1:])]

		# edge plotting settings: curve type 
		self.edge_df['raw_index'] = [i for i in range(len(self.edge_df.index))]
		self.edge_df['curve'] = self.edge_df.apply(
			lambda x: self.get_edge_style(x, ordered_edges, edge_dict), axis=1)
		self.edge_df.drop('raw_index', axis=1, inplace=True)

		# edge plotting settings: line type (dashed or not)
		self.edge_df['line'] = self.edge_df.apply(
			lambda x: self.get_edge_line(x, dataset_cols), axis=1)
		
	# calculates the positions and sizes of edges/nodes based on the 
	# number of nodes in the graph
	def calc_pos_sizes(self):

		pos = defaultdict()

		ordered_nodes = self.get_ordered_nodes()

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
		sub_node_size = node_size/2
		label_size = 43.9*(x**-0.484)

		# linearly-related sizes
		edge_width = -x/18 + (121/18)

		# assign fields to plotted graph object 
		self.pos = pos
		self.node_size = node_size
		self.sub_node_size = sub_node_size
		self.label_size = label_size
		self.rad_scale = 0.35
		self.edge_width = edge_width

	# get the shape of the node 
	def get_node_shape(self, x, dataset_cols):

		# default shape
		shape = 'o'

		# hexagonal nodes for combined nodes 
		if x.combined:
			shape = 'H'
			# unique to dataset
			# if is_novel(self.indicate_novel, x) or unique_to_dataset(self.indicate_dataset, x, dataset_cols): 
			# 	shape = 'h'
			# in dataset  TODO make these separate options in the future
			if is_novel(self.indicate_novel, x) or in_dataset(self.indicate_dataset, x):
				shape = 'h'
		# diamond nodes for novel or inidicate_dataset nodes
		# if is_novel(self.indicate_novel, x) or unique_to_dataset(self.indicate_dataset, x, dataset_cols):
		if is_novel(self.indicate_novel, x) or in_dataset(self.indicate_dataset, x):
			shape = 'D'

		return shape

	# get the node color 
	def get_node_color(self, x, color_dict):

		# if we're dealing with a path, color only according
		# to the node's role in the path
		if self.tid:
			self.path = self.get_path_from_tid(self.tid)
			if x.vertex_id == self.path[0]:
				x['color'] = color_dict['alt_TSS']['normal']
				x['sub_color'] = color_dict['alt_TSS']['normal']
			elif x.vertex_id == self.path[-1]:
				x['color'] = color_dict['alt_TES']['normal']
				x['sub_color'] = color_dict['alt_TES']['normal']
			elif x.vertex_id in self.path:
				x['color'] = color_dict['internal']['normal']
				x['sub_color'] = color_dict['internal']['normal']
			
			# if the vertex is not in the path
			else:
				if x.internal: color = color_dict['internal']['gray']
				if x.TSS: color = color_dict['TSS']['gray']
				if x.alt_TSS: color = color_dict['alt_TSS']['gray']
				if x.TES: color = color_dict['TES']['gray']
				if x.alt_TES: color = color_dict['alt_TES']['gray']
				x['color'] = color

		# combined nodes
		elif x.combined:
			types = x.combined_types

			# did one or two types of node go into this node?
			if len(types) == 2:

				x['sub_color'] = color_dict[types[0]]['normal']
				x['color'] = color_dict[types[1]]['normal']
			else: 
				x['sub_color'] = np.nan
				x['color'] = color_dict[types[0]]['normal']

		# non combined nodes
		else:

			# we don't need a sub color
			x['sub_color'] = np.nan

			# colors should label nodes by their 
			# MOST UNIQUE type (yes I know this is subjective) 
			# to me, this means that the label priority for a node is 
			# internal < TSS/alt_TSS < TES/alt_TES
			if x.internal: color = color_dict['internal']['normal']
			if x.TSS: color = color_dict['TSS']['normal']
			if x.alt_TSS: color = color_dict['alt_TSS']['normal']
			if x.TES: color = color_dict['TES']['normal']
			if x.alt_TES: color = color_dict['alt_TES']['normal']
			x['color'] = color

		return x

	# get the color of the edge
	def get_edge_color(self, x, color_dict):

		# firstly, if we're given a path, 
		# only color the edges that are in the path
		if self.path:
			path_edges = [(self.path[i],self.path[i+1])
						   for i in range(len(self.path)-1)]
			if x.edge_id in path_edges:
				color = color_dict[x.edge_type]['normal']
			else:
				color = color_dict[x.edge_type]['gray']

		else:
			color = color_dict[x.edge_type]['normal']

		return color

	# get the curve/style of the edge
	def get_edge_style(self, x, ordered_edges, edge_dict):

		# over 20 nodes, all should be curved
		if len(ordered_edges) < 20:
			if x.edge_id in ordered_edges:
				return edge_dict['straight']

		# make the arcs pretty
		## TODO maybe should have the arc height vary with number of nodes, ie 
		## shorter arc with less nodes? idk, still looks weird
		dist = self.pos[x.v2][0] - self.pos[x.v1][0]
		rad = self.rad_scale/dist

		if x.raw_index % 2 == 0:
			return edge_dict['pos']+str(rad)
		else:
			return edge_dict['neg']+str(rad)

	# get the style of the line for an edge
	def get_edge_line(self, x, dataset_cols):
		style = None
		# dashed if novel or unique to dataset
		# if is_novel(self.indicate_novel, x) or unique_to_dataset(self.indicate_dataset, x, dataset_cols):
		if is_novel(self.indicate_novel, x) or in_dataset(self.indicate_dataset, x):
			style = 'dashed'
		return style

	###############################################################################
	########################## Actually plotting stuff ############################
	###############################################################################

	# plots input plotted graph according to user's choices
	def plot_graph(self):

		# swan graph
		if not self.browser:
			self.plot_swan_graph()
		# browser track
		else:
			self.plot_browser_path()

	###############################################################################
	############################ Swan graph plotting ##############################
	###############################################################################

	# plots swan graph of current plotted graph
	def plot_swan_graph(self):

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
				style=entry.line)

	# plots nodes from loc_df
	def plot_nodes(self):
		for _, entry in self.loc_df.iterrows():
			node = entry.vertex_id
			nx.draw_networkx_nodes(self.G, self.pos,
				nodelist=[node],
				node_color=entry.color,
				node_size=self.node_size,
				node_shape=entry.node_shape)
			if entry.combined:
				if len(entry.combined_types) == 2:
					nx.draw_networkx_nodes(self.G, self.pos,
						nodelist=[node],
						node_color=entry.sub_color,
						node_size=self.sub_node_size,
						node_shape=entry.node_shape)


	###############################################################################
	######################## Browser track style plotting #########################
	###############################################################################

	# plots the browser track representation of the transcript from self.path
	def plot_browser_path(self):

		# which gene does this transcript come from?
		# what are the min/max coords in the gene?
		gid = self.get_gid_from_tid(self.tid)
		g_min, g_max = self.get_gene_min_max(gid)
		t_min, t_max = self.get_transcript_min_max(self.tid)
		g_len = g_max - g_min

		x_min = int(g_min-(6/(g_max-g_min)))
		x_max = int(g_max+(6/(g_max-g_min)))

		strand = self.get_strand_from_tid(self.tid)
		path = self.get_path_from_tid(self.tid)

		# plotting init
		plt.figure(1, figsize=(14,2.8), frameon=False)
		plt.xlim(x_min, x_max)
		plt.ylim(-1.05, 1.05)
		ax = plt.gca()
		teal = '#014753'

		# plot each exon as a rectangle
		y_coord = -0.1
		height = 0.2
		exons = [(v1,v2) for v1,v2 in zip(path[:-1],path[1:])][::2]
		introns = [(v1,v2) for v1,v2 in zip(path[:-1],path[1:])][1::2]
		for v1,v2 in exons:
			x_coord = self.loc_df.loc[v1, 'coord']
			width = self.loc_df.loc[v2, 'coord'] - x_coord
			rect = pch.Rectangle((x_coord,y_coord), width, height, color=teal)
			ax.add_patch(rect)

		# plot each intron as a hashed line
		dist = 0.001*g_len
		plt.plot([t_min+dist, t_max-dist], [0,0], color=teal)

		# reverse plot if we're on the minus strand
		if strand == '-':
			plt.gca().invert_xaxis()

		def get_tick_coords(loc_df, gene_len, exons, g_min, g_max):

			tick_coords = list(np.linspace(g_min, g_max, 40))

			# remove ticks that are before the start of the first exon
			start = loc_df.loc[exons[0][1], 'coord']
			tick_coords = [t for t in tick_coords if t > start]

			# remove ticks that are after the end of the last exon
			stop = loc_df.loc[exons[-1][1], 'coord']
			tick_coords = [t for t in tick_coords if t < stop]

			# remove ticks in and around the area of plotted exons
			dist = 0.001*gene_len
			for v1,v2 in exons:
				start = loc_df.loc[v1, 'coord']
				stop = loc_df.loc[v2, 'coord']
				tick_coords = [t for t in tick_coords
							   if t < start-dist or t > stop+dist]

			# if we only have one intron, and nothing made the cut, just stick
			# a tick in the middle of the intron
			if len(exons) == 2 and not tick_coords:
				tick_coords = [(loc_df.loc[exons[0][1],'coord']+loc_df.loc[exons[1][0],'coord'])/2]

			return tick_coords 

		# get coordinates for evenly-spaced ticks indicating strandedness
		# ala genome browser
		tick_coords = get_tick_coords(self.loc_df, g_len, exons, g_min, g_max)
		plt.plot(tick_coords, [0 for i in range(len(tick_coords))],
			color=teal, marker='4')

	# plots the scale of a browser style graph
	def plot_browser_scale(self):

		# get plot limits
		g_min = self.g_min
		g_max = self.g_max 

		plt.figure(1, figsize=(14,1), frameon=False)
		plt.xlim(g_min, g_max)
		plt.ylim(-1.05, 1.05)
		ax = plt.gca()

		y_coord = -0.5
		height = 1
		gene_len = g_max - g_min
		width = int((10/10000)*gene_len)

		# try different scale lengths until we've reached something optimal
		scale_len = 500
		while 1:
			if scale_len >= gene_len:
				break
			scale_len *= 10

		# dial it back so we show a nice chunk instead of the scale of the entire gene
		scale_len /= 10

		# are we going to indicate this in bp or kb?
		if scale_len < 1000: 
			scale_unit = 'bp'
			text_scale = scale_len
		elif scale_len < 10**6:
			scale_unit = 'kb'
			text_scale =  scale_len / 1000
		else:
			scale_unit = 'mb'
			text_scale = scale_len / 10**6
		text_scale = int(text_scale)

		# center the edges of the scale bar
		gene_mid_pt = (g_max+g_min)/2
		scale_start = gene_mid_pt - (scale_len/2)
		scale_end = gene_mid_pt + (scale_len/2)

		# start bar
		rect = pch.Rectangle((scale_start, y_coord),
			width=width,
			height=height,
			color='k')
		ax.add_patch(rect)

		# end bar
		rect = pch.Rectangle((scale_end, y_coord),
			width=width,
			height=height,
			color='k')
		ax.add_patch(rect)

		# plot the line in between the bars
		plt.plot([scale_start+width, scale_end-width], [0,0],
			color='k',
			linewidth=2)

		# add scale text
		x_coord = scale_start - (gene_len/9)
		y_coord = -0.2
		txt = '{} {}'.format(text_scale, scale_unit)
		ax.text(x_coord, y_coord, txt, fontsize=30)

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
		# for v in self.G.nodes:
		for v in self.loc_df.vertex_id.tolist():
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

			# some beginning of new loop things we need
			path = tuple([int(n) for n in path])
			start = path[0]
			stop = path[-1]
			combined_node = 'c{}'.format(combined_index)

			# get the colors for each aggregate node
			# TODO need to fix how these are colored with alt TSS/TES vs singleton TSS/TES
			combined_types = [j[0] for j in sorted([i for i in 
			   [('alt_TSS',loc_df.loc[path,'alt_TSS'].tolist().count(True)),
			   ('alt_TES',loc_df.loc[path,'alt_TES'].tolist().count(True)),
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

	# gets nodes ordered by genomic position 
	def get_ordered_nodes(self):
		loc_ids = self.loc_df.vertex_id.tolist()
		coords = self.loc_df.coord.tolist()
		ordered_nodes = [i[0] for i in sorted(zip(loc_ids, coords),
			key=lambda x: x[1])]
		return ordered_nodes

# is this entry unique from all other datasets?
# and in the provided dataset?
def unique_to_dataset(dataset_name, x, dataset_cols):
	if not dataset_name: return False
	if x[dataset_name] and not any(x[dataset_cols].tolist()):
		return True
	return False

# is this entry in the provided dataset?
def in_dataset(dataset_name, x):
	if not dataset_name: return False
	return x[dataset_name]

# is this entry annotated?
def is_novel(indicate_novel, x):
	if not indicate_novel: return False
	if x.annotation: return False
	return True



