import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import sqlite3
import SpliceGraph as sg
import PlottedGraph as pg
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
	for n in G.nodes():
		node_colors = defaultdict()
		sub_node_colors = defaultdict()
		node = G.nodes[n]

		# gray nodes
		if args['color_nodes'] == False:
			node_colors.update({'color': gray})

		# combined nodes
		elif args['combine'] and node['combined']:

			# all makeup of combined node is same type
			if len(node['combined_types']) == 1:
				color = color_dict[node['combined_types'][0]]
				node_colors.update({'color': color})

			# we need to plot more than one thing for combined node
			else:
				# color
				sub_color = color_dict[node['combined_types'][0]]
				color = color_dict[node['combined_types'][1]]
				node_colors.update({'color': color})
				sub_node_colors.update({'color': sub_color})

				# add combined sub-node color
				sub_node_style[n].update(sub_node_colors)

		# non-combined node
		else:
			node_colors.update({'color': color_dict['internal']})
			if node['TSS']: 
				node_colors.update({'color': color_dict['TSS']})
			if node['alt_TSS'] and args['color_alt_nodes']:
				node_colors.update({'color': color_dict['alt_TSS']})
			if node['TES']:
				node_colors.update({'color': color_dict['TES']})
			if node['alt_TES'] and args['color_alt_nodes']:
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
			connectionstyle=edge_style[e]['connectionstyle'])

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
	plt.figure(1, figsize=(14,10), frameon=False)
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

# plots a transcript path through a preexisiting 
# PlottedGraph graph
def plot_overlaid_path(pg, path, args, oname):

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
	path_pg.save_fig(oname)
