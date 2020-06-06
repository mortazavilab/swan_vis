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
		self.dataset_cols = []

		# plotting colors
		gray = '#DCDCDC'
		yellow = '#F0E442'
		light_blue = '#56B4E9'
		orange = '#E69F00'
		pink = '#CC79A7'
		green = '#009E73'

		gray_light_blue = '#c1ddec'
		gray_yellow = '#f8f6db'
		gray_orange = '#ebdec3'
		gray_pink = '#efdae5'
		gray_green = '#cdf5ea'

		self.color_dict = {'intron': {'normal': pink, 'gray': gray_pink},
					  'exon': {'normal': green, 'gray': gray_green},
					  'TSS': {'normal': light_blue, 'gray': gray_light_blue},
					  'TES': {'normal': orange, 'gray': gray_orange},
					  'internal': {'normal': yellow, 'gray': gray_yellow}}

	# initialize what needs to be from a previous pg 
	# or otherwise
	def init_plot_settings(self, sg,
						   gid=None, 
						   tid=None,
						   indicate_dataset=False,
						   indicate_novel=False,
						   browser=False):

		# save some old stuff
		old_indicate_dataset = self.indicate_dataset
		old_indicate_novel = self.indicate_novel
		old_gid = self.gid
		old_tid = self.tid
		old_browser = self.browser
		old_dataset_cols = self.dataset_cols

		# init some stuff
		self.indicate_dataset = indicate_dataset 
		self.indicate_novel = indicate_novel
		self.gid = gid 
		self.tid = tid
		self.browser = browser
		self.dataset_cols = sg.get_dataset_cols()

		# more human-readable graph types
		if self.tid:
			self.graph_type = 'transcript_path'
		else:
			self.graph_type = 'summary'


		# if we have new datasets
		if len(list(set(old_dataset_cols)-set(self.dataset_cols))) != 0:
			self.new_gene(sg)
			if not self.browser:
				self.calc_node_edge_styles()
				self.get_ordered_edges()

		# plotting the same transcript
		elif old_tid == self.tid and self.tid != None:
			if not browser: 
				if old_indicate_dataset != self.indicate_dataset \
				or old_indicate_novel != self.indicate_novel:
					self.calc_node_edge_styles()

		# plotting a different transcript
		elif old_tid != self.tid and self.tid != None:

			self.gid = sg.get_gid_from_tid(self.tid)
			if old_gid != self.gid:
				self.new_gene(sg)
			else:
				self.path = self.get_path_from_tid(self.tid)
			if not self.browser:
				self.calc_node_edge_styles()
				self.get_ordered_edges()

		# plotting a different gene
		elif old_gid != self.gid:
			self.new_gene(sg)
			if not self.browser:
				self.calc_node_edge_styles()

		elif old_indicate_dataset != self.indicate_dataset \
		or old_indicate_novel != self.indicate_novel:
			self.calc_node_edge_styles()

	# update pg to contain information for plotting a new gene
	def new_gene(self, sg):
		sg = subset_on_gene(sg, self.gid)
		self.loc_df = sg.loc_df
		self.edge_df = sg.edge_df
		self.t_df = sg.t_df
		self.G = sg.G

		self.g_min, self.g_max = self.get_gene_min_max(self.gid)
		self.strand = self.get_strand_from_gid(self.gid)

		if self.tid:
			self.path = self.get_path_from_tid(self.tid)

		if not self.browser:
			self.ordered_nodes = self.get_ordered_nodes()
			self.calc_pos_sizes()
			self.calc_edge_curves()

	# settings that can be initialized if we're plotting a new transcript
	# (or a new gene)
	def calc_node_edge_styles(self):
		self.calc_node_colors()
		self.calc_edge_colors()
		self.calc_edge_linestyles()

	###############################################################################
	################### Getting plotting settings for nodes/edges #################
	###############################################################################
		
	# calculates the positions and sizes of edges/nodes based on the 
	# number of nodes in the graph
	def calc_pos_sizes(self):

		pos = defaultdict()
		ordered_nodes = self.ordered_nodes

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
		edge_width = -x/18 + (121/18)
		if edge_width < 1:
			edge_width = 1

		# assign fields to plotted graph object 
		self.pos = pos
		self.node_size = node_size
		self.label_size = label_size
		self.rad_scale = 0.32
		self.edge_width = edge_width

	# calculate the color for each node
	def calc_node_colors(self):
		self.loc_df = self.loc_df.apply(
			lambda x: self.get_node_color(x), axis=1)

	# does the node satisfy indicate_novel or indicate_dataset
	# plotting settings?
	def is_novel_or_in_dataset(self, x):
		if self.indicate_novel and is_novel(x):
			if self.graph_type == 'transcript_path' \
			and x.vertex_id in self.path:
				return 'k'
			elif self.graph_type == 'summary':
				return 'k'
			else:
				return "#999999"
		elif self.indicate_dataset and in_dataset(self.indicate_dataset, x):
			if self.graph_type == 'transcript_path' \
			and x.vertex_id in self.path:
				return 'k'
			elif self.graph_type == 'summary':
				return 'k'
			else:
				return "#999999"
		return None

	# get the node color 
	def get_node_color(self, x):

		# transcript path through summary graph
		if self.graph_type == 'transcript_path':

			# vertices in the path should be colored by their roles in 
			# the plotted transcript
			if x.vertex_id == self.path[0]:
				x['color'] = self.color_dict['TSS']['normal']
			elif x.vertex_id == self.path[-1]:
				x['color'] = self.color_dict['TES']['normal']
			elif x.vertex_id in self.path:
				x['color'] = self.color_dict['internal']['normal']
			
			# vertices not in the path should be less colored
			# and colored according to their "most unique"
			# role in the current gene ie
			# internal < TSS < TES			
			else:
				if x.internal: color = self.color_dict['internal']['gray']
				if x.TSS: color = self.color_dict['TSS']['gray']
				if x.TES: color = self.color_dict['TES']['gray']
				x['color'] = color

		# gene summary graph
		else:

			# vertices in the summary graph should be colored 
			# according to their "most unique" 
			# role in the gene ie
			# internal < TSS < TES	
			if x.internal: color = self.color_dict['internal']['normal']
			if x.TSS: color = self.color_dict['TSS']['normal']
			if x.TES: color = self.color_dict['TES']['normal']
			x['color'] = color

		# if the node needs to be outlined due to indicate_dataset
		# or indicate_novel
		if self.indicate_novel or self.indicate_dataset:
			x['edgecolor'] = self.is_novel_or_in_dataset(x)
		else:
			x['edgecolor'] = None

		return x

	# get edge curve settings
	def calc_edge_curves(self):
		edge_dict = {'pos': 'arc3,rad=',
					 'neg': 'arc3,rad=-',
					 'straight': None}
		ordered_edges = [(n,m) for n,m in zip(self.ordered_nodes[:-1],
			self.ordered_nodes[1:])]
		self.edge_df['raw_index'] = [i for i in range(len(self.edge_df.index))]
		self.edge_df['curve'] = self.edge_df.apply(
			lambda x: self.get_edge_curve(x, ordered_edges, edge_dict), axis=1)
		self.edge_df.drop('raw_index', axis=1, inplace=True)

	# get the curve/style of the edge
	def get_edge_curve(self, x, ordered_edges, edge_dict):

		# over 20 nodes, all should be curved
		if len(ordered_edges) < 20:
			if x.edge_id in ordered_edges:
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
		color_dict = self.color_dict
		if self.tid:
			path_edges = [(self.path[i],self.path[i+1])
						   for i in range(len(self.path)-1)]
			if x.edge_id in path_edges:
				color = color_dict[x.edge_type]['normal']
			else:
				color = color_dict[x.edge_type]['gray']
		# otherwise just color them all
		else:
			color = color_dict[x.edge_type]['normal']
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
		# if is_novel(self.indicate_novel, x) or unique_to_dataset(self.indicate_dataset, x, dataset_cols):
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
				edgecolors=entry.edgecolor)

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
		path = self.path

		# plotting init
		plt.figure(1, figsize=(14,2.8), frameon=False)
		# plt.xlim(x_min, x_max)
		plt.xlim(g_min, g_max)
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

		# plot each intron as a line
		dist = 0.001*g_len
		plt.plot([t_min+dist, t_max-dist], [0,0], color=teal)

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

	# gets nodes ordered by genomic position 
	def get_ordered_nodes(self):
		loc_ids = self.loc_df.vertex_id.tolist()
		coords = self.loc_df.coord.tolist()
		ordered_nodes = [i[0] for i in sorted(zip(loc_ids, coords),
			key=lambda x: x[1])]

		# check if forward or reverse strand 
		if self.strand == '-':
			ordered_nodes.reverse()
		return ordered_nodes

	# orders edges by those present in the transcript and those not present in the transcript
	def get_ordered_edges(self):
		path_edges = [(self.path[i],self.path[i+1])
			   for i in range(len(self.path)-1)]
		self.edge_df['in_transcript'] = self.edge_df.apply(lambda x:
			1 if x.edge_id in path_edges else 0, axis=1)
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



