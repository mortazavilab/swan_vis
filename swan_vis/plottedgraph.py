import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import matplotlib.patches as pch
from swan_vis.utils import *
from swan_vis.graph import *
import time

class PlottedGraph(Graph):

	def __init__(self):
		self.indicate_dataset = False
		self.indicate_novel = False
		self.browser = False
		self.gid = None
		self.tid = None
		self.datasets = []

		# colors
		internal = '#F0E442'
		tss = '#56B4E9'
		tes = '#E69F00'
		intron = '#CC79A7'
		exon = '#009E73'
		browser = '#014753'

		tss_gray = '#c1ddec'
		internal_gray = '#f8f6db'
		tes_gray = '#ebdec3'
		intron_gray = '#efdae5'
		exon_gray = '#cdf5ea'

		self.color_dict = self.get_color_dict()

	def get_color_dict(self):
		"""
		Get the default color settings.

		Returns:
			color_dict (dict of str): Graphed object type to color
		"""

		# colors
		internal = '#F0E442'
		tss = '#56B4E9'
		tes = '#E69F00'
		intron = '#CC79A7'
		exon = '#009E73'
		browser = '#014753'

		tss_gray = '#c1ddec'
		internal_gray = '#f8f6db'
		tes_gray = '#ebdec3'
		intron_gray = '#efdae5'
		exon_gray = '#cdf5ea'

		color_dict = {'internal': internal, 'internal_gray': internal_gray,
						   'tss': tss, 'tss_gray': tss_gray,
						   'tes': tes, 'tes_gray': tes_gray,
						   'exon': exon, 'exon_gray': exon_gray,
						   'intron': intron, 'intron_gray': intron_gray,
						   'node_outline': 'k',
						   'node_outline_gray': '#999999',
						   'browser': browser,
						   None: None}
		return color_dict

	# initialize what needs to be from a previous pg
	# or otherwise
	def init_plot_settings(self, sg,
						   gid=None,
						   tid=None,
						   indicate_dataset=False,
						   indicate_novel=False,
						   browser=False):
		"""
		Initialize plotting settings from preexisting PlottedGraph. Try not to
		resubsample the input SwanGraph when unnecessary.

		Parameters:
			sg (swan SwanGraph): Parent SwanGraph to sample from
			gid (str): Gene that's being plotted
				Default: None
			tid (str): Transcript that's being plotted
				Default: None
			indicate_dataset (bool or str): Name of dataset to stylize nodes and
				edges from as outlined and dashed. If False, then don't do this
				for any nodes/edges
				Default: False
			indicate_novel (bool): Whether or not to stylize novel nodes and
				edges as outlined and dashed.
				Default: False
			browser (bool): Whether to plot the transcript as a browser model
				Default: False
		"""

		# previous plotting parameters
		old_gid = self.gid
		old_datasets = self.datasets

		# new plotting parameters
		self.indicate_dataset = indicate_dataset
		self.indicate_novel = indicate_novel
		self.gid = gid
		self.tid = tid
		self.browser = browser
		self.datasets = sg.datasets

		# determine graph_type
		if browser:
			self.graph_type = 'browser'

		if tid:
			self.graph_type = 'transcript_path'
			self.new_transcript(sg, self.tid, old_gid, old_datasets)
		else:
			self.graph_type = 'summary'

			# we need to update the gene
			if old_gid != self.gid or old_datasets != self.datasets:
				self.new_gene(sg, self.gid)

		if not browser:
			self.calc_node_edge_styles()

	def new_swangraph(self):
		"""
		Determine how the Swan plot for this gene is laid out.
		"""
		self.calc_pos_sizes()
		self.calc_edge_curves()

	def new_transcript(self, sg, tid, old_gid, old_datasets):
		"""
		Get the transcript edge and location path.
		"""
		self.gid = sg.get_gid_from_tid(tid)

		# if we need to update the gene
		if old_gid != self.gid or old_datasets != self.datasets:
			self.new_gene(sg, self.gid)

		self.edge_path = self.get_path_from_tid(self.tid)
		self.loc_path = self.get_loc_path_from_tid(self.tid)

	# update pg to contain information for plotting a new gene
	# includes subsetting the parent swangraph by gene and calculating
	# the node and edge positions and sizes
	def new_gene(self, sg, gid):
		"""
		Update PlottedGraph with information from a new gene.
		"""
		sg = sg.subset_on_gene(gid)
		self.loc_df = sg.loc_df
		self.edge_df = sg.edge_df
		self.t_df = sg.t_df
		self.G = sg.G

		self.g_min, self.g_max = self.get_gene_min_max(gid)
		self.strand = self.get_strand_from_gid(gid)

		# get positions and sizes of nodes/edges in this gene
		self.new_swangraph()

	# settings that can be initialized if we're plotting a new transcript
	# (or a new gene)
	def calc_node_edge_styles(self):
		self.calc_node_colors()
		self.calc_edge_colors()
		if self.graph_type == 'transcript_path':
			self.get_ordered_edges()
		self.calc_edge_linestyles()

	###############################################################################
	################### Getting plotting settings for nodes/edges #################
	###############################################################################

	# calculates the positions and sizes of edges/nodes based on the
	# number of nodes in the graph
	def calc_pos_sizes(self):
		"""
		Determine node position, size, label size, edgewidth (node).
		Determine edge width. Based on the number of nodes in the gene.
		"""

		pos = defaultdict()
		pos_list = []
		ordered_nodes = self.loc_df.vertex_id.tolist()

		y_coord = 0
		x_coords = np.linspace(-1, 1, len(ordered_nodes))
		for x_coord, n in zip(x_coords, ordered_nodes):
			pos[n] = [x_coord, y_coord]
			pos_list.append((x_coord, y_coord))

		x = len(ordered_nodes)

		# area-related sizes need to be non-linearly scaled
		# calculated by fitting power series curve to handpicked sizes
		node_size = 19248*(x**-1.14)
		label_size = 43.9*(x**-0.484)
		line_width = 2

		# linearly-related sizes
		edge_width = -x/18 + (121/18)
		if edge_width < 1:
			edge_width = 1

		# assign fields to plotted graph object
		self.pos = pos
		self.node_size = node_size
		self.label_size = label_size
		self.rad_scale = 0.32
		self.edge_width = edge_width
		self.line_width = line_width

	# calculate the color for each node
	def calc_node_colors(self):
		self.loc_df = self.loc_df.apply(
			lambda x: self.get_node_color(x), axis=1)

	# does the node satisfy indicate_novel or indicate_dataset
	# plotting settings?
	def is_novel_or_in_dataset(self, x):
		if self.indicate_novel and is_novel(x):
			if self.graph_type == 'transcript_path' \
			and x.vertex_id in self.loc_path:
				return 'node_outline', self.line_width
			elif self.graph_type == 'summary':
				return 'node_outline', self.line_width
			else:
				return 'node_outline_gray', self.line_width
		elif self.indicate_dataset and in_dataset(self.indicate_dataset, x):
			if self.graph_type == 'transcript_path' \
			and x.vertex_id in self.loc_path:
				return 'node_outline', self.line_width
			elif self.graph_type == 'summary':
				return 'node_outline', self.line_width
			else:
				return 'node_outline_gray', self.line_width
		return None, None

	# get the node color
	def get_node_color(self, x):

		# transcript path through summary graph
		if self.graph_type == 'transcript_path':

			# vertices in the path should be colored by their roles in
			# the plotted transcript
			if x.vertex_id == self.loc_path[0]:
				x['color'] = 'tss'
			elif x.vertex_id == self.loc_path[-1]:
				x['color'] = 'tes'
			elif x.vertex_id in self.loc_path:
				x['color'] = 'internal'

			# vertices not in the path should be less colored
			# and colored according to their "most unique"
			# role in the current gene ie
			# internal < TSS < TES
			else:
				if x.internal: color = 'internal_gray'
				if x.TSS: color = 'tss_gray'
				if x.TES: color = 'tes_gray'
				x['color'] = color

		# gene summary graph
		else:

			# vertices in the summary graph should be colored
			# according to their "most unique"
			# role in the gene ie
			# internal < TSS < TES
			if x.internal: color = 'internal'
			if x.TSS: color = 'tss'
			if x.TES: color = 'tes'
			x['color'] = color

		# if the node needs to be outlined due to indicate_dataset
		# or indicate_novel
		if self.indicate_novel or self.indicate_dataset:
			ecolor, lwidth = self.is_novel_or_in_dataset(x)
			x['edgecolor'] = ecolor
			x['linewidth'] = lwidth
		else:
			x['edgecolor'] = None
			x['linewidth'] = None

		return x

	def calc_edge_curves(self):
		"""
		Add edge curve style to edge_df in the curve column. Sequential edges
		in small genes (< 20 nodes) will be straight.
		"""

		edge_dict = {'pos': 'arc3,rad=',
					 'neg': 'arc3,rad=-',
					 'straight': None}

		self.edge_df['raw_index'] = [i for i in range(len(self.edge_df.index))]
		self.edge_df['curve'] = self.edge_df.apply(
			lambda x: self.get_edge_curve(x, edge_dict), axis=1)
		self.edge_df.drop('raw_index', axis=1, inplace=True)

	# get the curve/style of the edge
	def get_edge_curve(self, x, edge_dict):

		ordered_nodes = self.loc_df.vertex_id.tolist()

		# over 20 nodes, all should be curved
		if len(ordered_nodes) < 20:
			if abs(x.v1 - x.v2) == 1:
				return edge_dict['straight']

		# make the arcs pretty
		dist = self.pos[x.v2][0] - self.pos[x.v1][0]
		rad = self.rad_scale/dist

		if x.raw_index % 2 == 0:
			return edge_dict['pos']+str(rad)
		else:
			return edge_dict['neg']+str(rad)

	# get the colors of all the edges
	def calc_edge_colors(self):
		self.edge_df['color'] = self.edge_df.apply(
			lambda x: self.get_edge_color(x), axis=1)

	# get the color of the edge
	def get_edge_color(self, x):
		# firstly, if we're given a path,
		# only color the edges that are in the path
		# color_dict = self.color_dict
		if self.tid:
			# path_edges = [(self.path[i],self.path[i+1])
			#				for i in range(len(self.path)-1)]
			if x.edge_id in self.edge_path:
				color = x.edge_type
			else:
				color = x.edge_type+'_gray'
		# otherwise just color them all
		else:
			color = x.edge_type
		return color

	# get the styles of the edges (dashed or not)
	def calc_edge_linestyles(self):
		self.edge_df['line'] = self.edge_df.apply(
			lambda x: self.get_edge_linestyle(x),
			axis=1)

	# get the style of the line for an edge
	def get_edge_linestyle(self, x):
		style = None
		# dashed if novel or unique to dataset
		if self.indicate_novel:
			if is_novel(x):
				style = 'dashed'
		elif self.indicate_dataset:
			if in_dataset(self.indicate_dataset, x):
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

		# remove axis
		plt.axis('off')

	# plots edges from edge_df
	def plot_edges(self):
		for _, entry in self.edge_df.iterrows():
			edge = (entry.v1, entry.v2)
			nx.draw_networkx_edges(self.G, self.pos,
				edgelist=[edge],
				width=self.edge_width,
				edge_color=self.color_dict[entry.color],
				connectionstyle=entry.curve,
				style=entry.line)

	# plots nodes from loc_df
	def plot_nodes(self):
		for _, entry in self.loc_df.iterrows():
			node = entry.vertex_id
			nx.draw_networkx_nodes(self.G, self.pos,
				nodelist=[node],
				node_color=self.color_dict[entry.color],
				node_size=self.node_size,
				edgecolors=self.color_dict[entry.edgecolor],
				linewidths=entry.linewidth)

	###############################################################################
	######################## Browser track style plotting #########################
	###############################################################################

	# plots the browser track representation of the transcript from self.path
	def plot_browser_path(self):

		# which gene does this transcript come from?
		# what are the min/max coords in the gene?
		gid = self.gid
		g_min = self.g_min
		g_max = self.g_max
		t_min, t_max = self.get_transcript_min_max(self.tid)
		g_len = g_max - g_min

		x_min = int(g_min-(6/(g_max-g_min)))
		x_max = int(g_max+(6/(g_max-g_min)))

		strand = self.strand
		loc_path = self.loc_path
		edge_path = self.edge_path

		# plotting init
		plt.figure(1, figsize=(14,2.8), frameon=False)
		# plt.xlim(x_min, x_max)
		plt.xlim(g_min, g_max)
		plt.ylim(-1.05, 1.05)
		ax = plt.gca()
		color = self.color_dict['browser']

		# plot each exon as a rectangle
		y_coord = -0.1
		height = 0.2
		exons = [(v1,v2) for v1,v2 in zip(loc_path[:-1],loc_path[1:])][::2]
		introns = [(v1,v2) for v1,v2 in zip(loc_path[:-1],loc_path[1:])][1::2]
		for v1,v2 in exons:
			x_coord = self.loc_df.loc[v1, 'coord']
			width = self.loc_df.loc[v2, 'coord'] - x_coord
			rect = pch.Rectangle((x_coord,y_coord), width, height, color=color)
			ax.add_patch(rect)

		# plot each intron as a line
		dist = 0.001*g_len
		plt.plot([t_min+dist, t_max-dist], [0,0], color=color)

		# reverse plot if we're on the minus strand
		if strand == '-':
			plt.gca().invert_xaxis()

		# remove axis
		plt.axis('off')

		def get_tick_coords(loc_df, gene_len, exons, g_min, g_max, strand):

			tick_coords = list(np.linspace(g_min, g_max, 40))

			# remove ticks that are before the start of the first exon
			start = loc_df.loc[exons[0][1], 'coord']
			if strand == '+':
				tick_coords = [t for t in tick_coords if t > start]
			else:
				tick_coords = [t for t in tick_coords if t < start]

			# remove ticks that are after the end of the last exon
			stop = loc_df.loc[exons[-1][1], 'coord']
			if strand == '+':
				tick_coords = [t for t in tick_coords if t < stop]
			else:
				tick_coords = [t for t in tick_coords if t > stop]

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
		tick_coords = get_tick_coords(self.loc_df, g_len, exons,
			g_min, g_max, self.strand)
		plt.plot(tick_coords, [0 for i in range(len(tick_coords))],
			color=color, marker='4')

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

		# three digit numbers don't look so good
		if len(str(text_scale)) > 2:
			x_coord = scale_start - (gene_len/7)
		ax.text(x_coord, y_coord, txt, fontsize=30)

		# # check if forward or reverse strand
		# if self.strand == '-':
		#	 ordered_nodes.reverse()
		#
		# return ordered_nodes

	# orders edges by those present in the transcript and those not present in the transcript
	def get_ordered_edges(self):
		self.edge_df['in_transcript'] = self.edge_df.apply(lambda x:
			1 if x.edge_id in self.edge_path else 0, axis=1)
		self.edge_df.sort_values(by='in_transcript', inplace=True)
		self.edge_df.drop('in_transcript', axis=1, inplace=True)

# is this entry unique from all other datasets?
# and in the provided dataset?
def unique_to_dataset(dataset_name, x, dataset_cols):
	if not dataset_name: return False
	if x[dataset_name] and not any(x[dataset_cols].tolist()):
		return True
	return False

# is this entry in the provided dataset?
def in_dataset(dataset_name, x):
	return x[dataset_name]

# is this entry annotated?
def is_novel(x):
	if x.annotation: return False
	return True
