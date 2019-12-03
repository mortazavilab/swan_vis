import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import os
import copy
from collections import defaultdict
import sqlite3
import SpliceGraph as sg
import PlottedGraph as pg
from report_tools import *
from utils import *

# updates the style dict with colors during actual plotting calls
def get_node_colors(G, node_style, sub_node_style, args):

	# plotting styles
	gray = '#DCDCDC'
	yellow = '#F0E442'
	blue = '#0072B2'
	light_blue = '#56B4E9'
	red = '#D55E00'
	orange = '#E69F00'
	pink = '#CC79A7'
	green = '#009E73'
	color_dict = {'internal': yellow, 'TSS': blue,
				  'alt_TSS': light_blue, 'TES': red,
				  'alt_TES': orange}

	# assign nodes colors based on plotting settings
	for n,data in G.nodes(data=True):
		node_colors = defaultdict()
		sub_node_colors = defaultdict()
		node = G.nodes[n]

		# gray nodes
		if args['color_nodes'] == False:
			node_colors.update({'color': gray})

		# combined nodes
		elif args['combine'] and data['combined']:

			# all makeup of combined node is same type
			if len(data['combined_types']) == 1:
				color = color_dict[data['combined_types'][0]]
				node_colors.update({'color': color})

			# we need to plot more than one thing for combined node
			else:
				# color
				sub_color = color_dict[data['combined_types'][0]]
				color = color_dict[data['combined_types'][1]]
				node_colors.update({'color': color})
				sub_node_colors.update({'color': sub_color})

				# add combined sub-node color
				sub_node_style[n].update(sub_node_colors)

		# non-combined node
		else:
			node_colors.update({'color': color_dict['internal']})
			if data['TSS']: 
				node_colors.update({'color': color_dict['TSS']})
			if data['alt_TSS'] and args['color_alt_nodes']:
				node_colors.update({'color': color_dict['alt_TSS']})
			if data['TES']:
				node_colors.update({'color': color_dict['TES']})
			if data['alt_TES'] and args['color_alt_nodes']:
				node_colors.update({'color': color_dict['alt_TES']})

		# add color to normal node style
		node_style[n].update(node_colors)

	return node_style, sub_node_style

# returns a dictionary indexed by edge ids of plotting styles 
def get_edge_colors(G, edge_style, args):

	# colors
	pink = '#CC79A7'
	green = '#009E73'
	gray = '#DCDCDC'
	gray = '#DCDCDC'
	yellow = '#F0E442'
	blue = '#0072B2'
	light_blue = '#56B4E9'
	red = '#D55E00'
	orange = '#E69F00'
	color_dict = {'intron': pink, 'exon': green}

	# assign a color for each edge
	edges = list(G.edges)
	for e in edges:
		edge_colors = {}
		if args['color_edges']:
			if G.edges[e]['edge_type'] == 'intron': 
				edge_colors.update({'color': color_dict['intron']})
			elif G.edges[e]['edge_type'] == 'exon':
				edge_colors.update({'color': color_dict['exon']})
		else: 
			edge_colors.update({'color': gray})

		edge_style[e].update(edge_colors)
	return edge_style

# plot edges according to settings in style_dict
def plot_edges(G, pos, edge_style):
	for e in G.edges():
		nx.draw_networkx_edges(G, pos,
			edgelist=[e],
			width=edge_style[e]['width'],
			edge_color=edge_style[e]['color'],
			connectionstyle=edge_style[e]['connectionstyle'],
			linestyle=edge_style[e]['linestyle'])

# plot nodes according to style dicts
def plot_nodes(G, pos, node_style, sub_node_style):

	for n in G.nodes():
		nx.draw_networkx_nodes(G, pos,
			nodelist=[n],
			node_color=node_style[n]['color'],
			node_size=node_style[n]['size'],
			node_shape=node_style[n]['shape'])
		if n in sub_node_style.keys():
			nx.draw_networkx_nodes(G, pos,
				nodelist=[n],
				node_color=sub_node_style[n]['color'],
				node_size=sub_node_style[n]['size'],
				node_shape=sub_node_style[n]['shape'])

# plots input plotted graph
def plot_graph(pg, args):

	# plotting stuff
	plt.figure(1, figsize=(14,2.8), frameon=False)
	plt.xlim(-1.05, 1.05)
	plt.ylim(-1.05, 1.05) 

	# get node/edge colors based on arguments
	node_style, sub_node_style = get_node_colors(pg.G, pg.node_style, pg.sub_node_style, args)
	edge_style = get_edge_colors(pg.G, pg.edge_style, args)

	# plot edges, nodes, and labels
	plot_edges(pg.G, pg.pos, pg.edge_style)

	# relabel nodes with strings for the purpose of labelling them
	relabel_map = {k:(k if type(k) == str else int(k)) for k in pg.G.nodes}
	pg.G = nx.relabel_nodes(pg.G, relabel_map)
	nx.draw_networkx_labels(pg.G, pg.pos, font_size=pg.label_size)

	plot_nodes(pg.G, pg.pos, node_style, sub_node_style)

# plots a transcript path through a preexisiting PlottedGraph graph
def plot_overlaid_path(pg, path, args):

	# first plot the preexisting graph in gray
	args['color_edges'] = args['color_nodes'] = args['color_alt_nodes'] = False
	plot_graph(pg, args)
	
	# get nodes and edges from the input path
	path_nodes = path
	path_edges = [(path[i], path[i+1]) for i in range(len(path)-1)]

	# create a subgraph based on the nodes and edges in this path
	path_pg = copy.deepcopy(pg)
	nodes = [(n,d) for n,d in pg.G.nodes(data=True)
    		 if n in path_nodes]
	edges = [(v1,v2,d) for v1,v2,d in pg.G.edges(data=True)
			 if (v1,v2) in path_edges]
    
	G = nx.DiGraph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)
	path_pg.G = G

	# plot the subgraph of this transcript path on top of 
	# preexisting graph
	args['color_edges'] = args['color_nodes'] = args['color_alt_nodes'] = True
	plot_graph(path_pg, args)
	# save_fig(oname)

# plot a transcript path as we'd see on the genome browser
def plot_path_browser(splice_graph, tid):

	# plotting init stuff
	plt.figure(1, figsize=(14,2.8), frameon=False)
	plt.xlim(-1.05, 1.05)
	plt.ylim(-1.05, 1.05)
	ax = plt.gca() 
	teal = '#014753'

	# fields we'll need
	loc_df = splice_graph.loc_df
	t_df = splice_graph.t_df

	# which gene is this from?
	gid = t_df.loc[tid, 'gid']

	# what are the strand and chrom of the gene?
	chrom = loc_df.loc[t_df.loc[tid, 'path'][0], 'chrom']
	strand = loc_df.loc[t_df.loc[tid, 'path'][0], 'strand']

	# what are the min and max of the gene? use these to define the coord space
	g_min, g_max = sg.get_gene_min_max(loc_df, t_df, gid)
	coord_map = get_coord_map(g_min, g_max, loc_df, chrom, strand)

	# plot each exon as a rectangle
	y_coord = -0.1
	height = 0.2
	path = t_df.loc[tid, 'path']
	exons = [(v1,v2) for v1,v2 in zip(path[:-1],path[1:])][::2]
	introns = [(v1,v2) for v1,v2 in zip(path[:-1],path[1:])][1::2]
	for v1,v2 in exons:
		x_coord = coord_map.loc[v1, 'plot_coord']
		width = coord_map.loc[v2, 'plot_coord'] - x_coord
		rect = pch.Rectangle((x_coord,y_coord),
							  width, height,
							  color=teal)
		ax.add_patch(rect)

	# plot each intron as a hashed line
	for v1,v2 in introns: 
		x_coords = coord_map.loc[[v1,v2], 'plot_coord']
		plt.plot(x_coords, [0,0], color=teal)
	tick_coords = get_tick_coords(exons, coord_map)
	# for t in tick_coords: 
	plt.plot(tick_coords,[0 for i in range(len(tick_coords))],
			 color=teal, marker='4')

# maps vertex id to a matplotlib-digestible coordinate
def get_coord_map(g_min, g_max, loc_df, chrom, strand):

	# reformat loc_df for indexing by coord, chrom, strand
	loc_df = reset_dupe_index(loc_df, 'vertex_id')
	loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)
	inds = loc_df.index.tolist()

	# plot coordinates for each genomic coordinate
	# use .9 to allow for some room on the sides
	plot_coords = list(np.linspace(-0.95, 0.95, g_max-g_min+1)) 

	# reverse if on minus strand 
	if strand == '-':
		plot_coords.reverse()

	# map the locations vertex id to the plotting coord
	g_coords = range(g_min, g_max+1)
	coord_map = pd.DataFrame({'plot_coord': plot_coords, 'g_coord': g_coords})
	coord_map['vertex_id'] = coord_map.apply(lambda x:
									 	     loc_df.loc[(chrom,x.g_coord,strand), 'vertex_id']
									 		 if (chrom,x.g_coord,strand) in inds
									 		 else np.nan, axis=1)

	# reformat loc_df for future use
	loc_df.reset_index(inplace=True)
	loc_df = set_dupe_index(loc_df, 'vertex_id')

	# remove nan rows and set index
	coord_map.dropna(inplace=True)
	coord_map.set_index('vertex_id', inplace=True)

	return coord_map

# get locations of ticks to indicate strandedness
def get_tick_coords(exons, coord_map):
	tick_coords = list(np.linspace(-0.95, 0.95, 40))

	# remove ticks that are before the start of the first exon
	start = coord_map.loc[exons[0][1], 'plot_coord']
	tick_coords = [t for t in tick_coords if t > start]

	# remove ticks that are after the end of the last exon
	stop = coord_map.loc[exons[-1][1], 'plot_coord']
	tick_coords = [t for t in tick_coords if t < stop]

	# remove ticks in and around the area of plotted exons
	for v1,v2 in exons:
		start = coord_map.loc[v1, 'plot_coord']
		stop = coord_map.loc[v2, 'plot_coord']
		tick_coords = [t for t in tick_coords 
					   if t < start-0.002 or t > stop+0.002]
	return tick_coords

# plot each transcript in the splice graph overlaid on the full graph
def plot_each_transcript(splice_graph, args, oprefix, browser=False):
	if not browser:
		plotted_graph = pg.PlottedGraph(splice_graph, args)
		for tid in plotted_graph.t_df.tid.tolist():
			oname = '{}_{}.png'.format(oprefix, tid)
			entry = plotted_graph.t_df.loc[tid]
			path = entry['path']
			plot_overlaid_path(plotted_graph, path, args)
			save_fig(oname)
	else:
		for tid in splice_graph.t_df.tid.tolist():
			oname = '{}_{}_browser.png'.format(oprefix, tid)
			plot_path_browser(splice_graph, tid)
			save_fig(oname)

# saves current figure named oname. clears the figure space so additional
# plotting can be done
def save_fig(oname):
	plt.axis('off')
	plt.tight_layout()
	plt.savefig(oname, format='png', dpi=200)
	plt.clf()
	plt.close()
