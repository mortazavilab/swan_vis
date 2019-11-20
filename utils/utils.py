import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
from collections import defaultdict
import sqlite3

# class SpliceGraph:
# 	def __init__(self, gtf=None, talon_db=None):

# 		if not gtf and not talon_db:
# 			raise Exception('No input to SpliceGraph given.')

# 		# GTF
# 		if gtf:
# 			if not os.path.exists(gtf):
# 				raise Exception('GTF file not found. Check path.')
# 			loc_df, edge_df, t_df = self.gtf_create_dfs(gtf)
# 		if talon_db: 
# 			if not os.path.exists(talon_db):
# 				raise Exception('TALON db file not found. Check path.')
# 			loc_df, edge_df, t_df = self.db_create_dfs(talon_db)

# 		self.G = self.create_graph_from_dfs(loc_df, edge_df, t_df)
# 		self.loc_df = loc_df
# 		self.edge_df = edge_df
# 		self.t_df = t_df

# 	# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
# 	def gtf_create_dfs(self, gtffile):

# 		# get dfs by parsing through gtf
# 		loc_df = pd.DataFrame(columns=['chrom', 'coord',
# 									   'strand','vertex_id',
# 									   'TSS', 'alt_TSS',
# 									   'TES', 'alt_TES',
# 									   'internal'])
# 		loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)

# 		edge_df = pd.DataFrame(columns=['edge_id', 'edge_type',
# 									    'strand', 'v1', 'v2'])
# 		t_df = pd.DataFrame(columns=['tid', 'gid',
# 									 'gname', 'path'])

# 		# loop initialization
# 		vertex_id = 1
# 		transcript_paths = []
# 		transcript_path = []

# 		with open(gtffile, 'r') as gtf:
# 			for line in gtf:

# 				# skip header lines
# 				if '##' in line: continue

# 				line = line.strip().split('\t')

# 				# gene entry 
# 				if line[2] == 'gene':
# 					curr_gid = get_field_value('gene_id', line[-1])
# 					curr_gname = get_field_value('gene_name', line[-1])

# 				# transcript entry
# 				elif line[2] == 'transcript':
# 					curr_tid = get_field_value('transcript_id', line[-1])
					
# 					# start a new transcript path
# 					if transcript_path != []:

# 						# add to list of transcript paths and transcript df 
# 						transcript_paths.append(transcript_path)
# 						t_df = t_df.append({'tid': prev_tid,
# 									 'gid': prev_gid,
# 									 'gname': prev_gname,
# 									 'path': transcript_path},
# 									 ignore_index=True)

# 					transcript_path = []

# 					# reset some stuff
# 					terminal_loc = True
# 					exon = 0
# 					intron = 1

# 				# exon entry
# 				elif line[2] == 'exon':

# 					# get exon info 
# 					curr_chr = line[0]
# 					curr_start = line[3]
# 					curr_stop = line[4]
# 					curr_strand = line[6]
					
# 					if curr_strand == '+': coords = [curr_start, curr_stop]
# 					else: coords = [curr_stop, curr_start]
					
# 					for c in coords:

# 						ind = (curr_chr, int(c), curr_strand)

# 						# loc not in loc_df already
# 						if ind not in loc_df.index.tolist():

# 							# label as not a TSS/TES until further notice
# 							attr = {'vertex_id': vertex_id,	   
# 									'TSS': False, 'TES': False,
# 									'alt_TSS': False,
# 									'alt_TES': False, 
# 									'internal': False, 'coord': int(c),
# 									'strand': curr_strand, 'chrom': curr_chr}

# 							# update loc_df and increment vertex_id
# 							loc_df.reset_index(inplace=True)
# 							loc_df = loc_df.append(attr, ignore_index=True)
# 							loc_df.set_index(['chrom', 'coord', 'strand'], inplace=True)

# 							curr_loc = int(vertex_id)
# 							vertex_id += 1

# 						# loc was already added to graph
# 						else: curr_loc = int(loc_df.loc[ind].vertex_id)	
		
# 						# add an edge to previous loc if not terminal loc 
# 						# and if the edge doesn't already exist
# 						if not terminal_loc:
# 							curr_edge = (prev_loc, curr_loc)
							
# 							if curr_edge not in edge_df.edge_id.to_list():
# 								attrs = {'edge_id': (curr_edge[0], curr_edge[1]),
# 									     'v1': curr_edge[0],
# 										 'v2': curr_edge[1], 
# 										 'strand': curr_strand}
# 								if exon: attrs.update({'edge_type': 'exon'})
# 								elif intron: attrs.update({'edge_type': 'intron'})

# 								edge_df = edge_df.append(attrs, ignore_index=True)

# 						# update transcript path with each loc 
# 						transcript_path.append(curr_loc)
# 						prev_loc = curr_loc
# 						prev_tid = curr_tid
# 						prev_gid = curr_gid
# 						prev_gname = curr_gname
# 						terminal_loc = False
						
# 						# exon or intron
# 						exon = abs(exon-1)
# 						intron = abs(intron-1)
						
# 		# append last transcript info
# 		transcript_paths.append(transcript_path)
# 		t_df = t_df.append({'tid': curr_tid,
# 						    'gid': curr_gid,
# 						    'gname': curr_gname,
# 							'path': transcript_path},
# 							ignore_index=True)

# 		# label node/edge types and finish formatting dfs correctly
# 		loc_df.reset_index(inplace=True)
# 		loc_df = create_dupe_index(loc_df, 'vertex_id')
# 		loc_df = set_dupe_index(loc_df, 'vertex_id')
# 		loc_df = self.get_loc_types(loc_df, t_df)

# 		edge_df['annotated'] = True # we can assume that since we're working from a gtf, it's an annotation?? (maybe)
# 		loc_df['annotated'] = True

# 		t_df = create_dupe_index(t_df, 'tid')
# 		t_df = set_dupe_index(t_df, 'tid')
# 		edge_df.set_index('edge_id', inplace=True)

# 		return loc_df, edge_df, t_df

# 	# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
# 	def db_create_dfs(self, db):

# 		# open db connection
# 		conn = sqlite3.connect(db)
# 		c = conn.cursor()

# 		# loc_df
# 		q = 'SELECT loc.* FROM location loc'

# 		c.execute(q)
# 		locs = c.fetchall()

# 		loc_df = pd.DataFrame(locs,
# 			columns=['location_ID', 'genome_build',
# 					 'chrom', 'position'])

# 		# do some df reformatting, add strand
# 		loc_df.drop('genome_build', axis=1, inplace=True)
# 		loc_df.rename({'location_ID': 'vertex_id',
# 					   'position': 'coord'},
# 					   inplace=True, axis=1)
# 		loc_df.vertex_id = loc_df.vertex_id.map(int)

# 		# edge_df
# 		q = """SELECT e.* 
# 				FROM edge e 
# 				JOIN vertex V ON e.v1=v.vertex_ID 
# 				JOIN gene_annotations ga ON v.gene_ID=ga.ID 
# 				WHERE ga.attribute='gene_name'
# 			""" 

# 		c.execute(q)
# 		edges = c.fetchall()

# 		edge_df = pd.DataFrame(edges, 
# 			columns=['edge_id', 'v1', 'v2',
# 					 'edge_type', 'strand'])
# 		edge_df.v1 = edge_df.v1.map(int)
# 		edge_df.v2 = edge_df.v2.map(int)
# 		edge_df['talon_edge_id'] = edge_df.edge_id
# 		edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

# 		# t_df
# 		t_df = pd.DataFrame()

# 		# get tid, gid, gname, and paths
# 		q = """SELECT ga.value, ta.value,
# 					  t.start_exon, t.jn_path, t.end_exon,
# 					  t.start_vertex, t.end_vertex
# 				FROM gene_annotations ga 
# 				JOIN transcripts t ON ga.ID=t.gene_ID
# 				JOIN transcript_annotations ta ON t.transcript_ID=ta.ID
# 				WHERE ta.attribute='transcript_id'
# 				AND (ga.attribute='gene_name' 
# 				OR ga.attribute='gene_id')
# 			"""

# 		c.execute(q)
# 		data = c.fetchall()

# 		# get fields from each transcript and add to dataframe
# 		gids, tids, paths = zip(*[(i[0], i[1], i[2:]) for i in data[::2]])
# 		gnames = [i[0] for i in data[1::2]]
# 		paths = self.get_edge_paths(paths)

# 		t_df['tid'] = np.asarray(tids)
# 		t_df['gid'] = np.asarray(gids)
# 		t_df['gname'] = np.asarray(gnames)
# 		t_df['path'] = np.asarray(paths)

# 		t_df = create_dupe_index(t_df, 'tid')
# 		t_df = set_dupe_index(t_df, 'tid')

# 		# furnish the last bit of info in each df
# 		t_df['path'] = [[int(n) for n in path]
# 						 for path in self.get_vertex_paths(paths, edge_df)]
# 		loc_df['strand'] = loc_df.apply(lambda x:
# 				 self.get_strand(x, edge_df), axis=1)
# 		loc_df = create_dupe_index(loc_df, 'vertex_id')
# 		loc_df = set_dupe_index(loc_df, 'vertex_id')
# 		loc_df['internal'] = False
# 		loc_df['TSS'] = False
# 		loc_df['alt_TSS'] = False
# 		loc_df['TES'] = False
# 		loc_df['alt_TES'] = False
# 		loc_df['annotated'] = True
# 		loc_df = self.get_loc_types(loc_df, t_df)

# 		edge_df['annotated'] = True
# 		edge_df.drop('talon_edge_id', axis=1, inplace=True)
# 		edge_df = create_dupe_index(edge_df, 'edge_id')
# 		edge_df = set_dupe_index(edge_df, 'edge_id')

# 		return loc_df, edge_df, t_df

# 	# convert talon query into edge path
# 	def get_edge_paths(self, paths):
# 		edge_paths = []
# 		for p in paths:
# 			if p[1] == None:
# 				edge_paths.append([p[0]])
# 			else:
# 				edge_paths.append(
# 					[p[0], *[int(i) for i in p[1].split(',')], p[2]])
# 		return edge_paths

# 	# convert edge path to vertex path
# 	def get_vertex_paths(self, paths, edge_df):
# 		vertex_paths = []
# 		for p in paths: 
# 			path = []
# 			for i, e in enumerate(p): 
# 				entry = edge_df.loc[edge_df.talon_edge_id == e]
# 				if i == 0:
# 					path.extend([entry.v1.values[0], entry.v2.values[0]])
# 				else: path.append(entry.v2.values[0])
# 			vertex_paths.append(path)
# 		return vertex_paths

# 	# get the strand of each vertex
# 	def get_strand(self, x, edge_df):
# 		# use v1 or v2 depending on where vertex is in edge
# 		try: 
# 			strand = edge_df.loc[edge_df.v1 == x.vertex_id, 'strand'].values[0]
# 		except:
# 			strand = edge_df.loc[edge_df.v2 == x.vertex_id, 'strand'].values[0]
# 		return strand

# 	# add node types (internal, TSS, alt TSS, TES, alt_TES) to loc_df
# 	def get_loc_types(self, loc_df, t_df):

# 		# label each location as internal off the bat, and not as TSS/TES
# 		loc_df['internal'] = False
# 		loc_df['TSS'] = False
# 		loc_df['TES'] = False
# 		loc_df['alt_TSS'] = False
# 		loc_df['alt_TES'] = False

# 		# label each TSS and TES
# 		paths = t_df.path.tolist()
# 		tss = np.unique([path[0] for path in paths])
# 		loc_df.loc[tss, 'TSS'] = True
# 		tes = np.unique([path[-1] for path in paths])
# 		loc_df.loc[tes, 'TES'] = True
# 		internal = np.unique([n for path in paths for n in path[1:-1]])
# 		loc_df.loc[internal, 'internal'] = True

# 		# label each alt TSS and alt TES for each gene
# 		for g in t_df.gid.unique().tolist():
# 			gene_entries = t_df.loc[t_df.gid == g]

# 			# genes that have more than one transcript are alt TSS/TES candidates
# 			if len(gene_entries.index) != 1: 

# 				paths = gene_entries.path.tolist()
# 				tss = [path[0] for path in paths]
# 				tes = [path[-1] for path in paths]

# 				# alt TSS/TES
# 				if len(set(tss)) > 1: 
# 					loc_df.loc[tss, 'alt_TSS'] = True
# 				if len(set(tes)) > 1: 
# 					loc_df.loc[tes, 'alt_TES'] = True

# 		return loc_df

# 	# create the graph object from the dataframes
# 	def create_graph_from_dfs(self, loc_df, edge_df, t_df):

# 		# graph initialization
# 		G = nx.DiGraph()

# 		# add nodes to graph from transcript paths
# 		paths = t_df.path.tolist()
# 		for path in paths:
# 			nx.add_path(G, path)

# 		# add node attributes from dfs
# 		G = label_nodes(G, loc_df, 'internal', 'internal') 
# 		G = label_nodes(G, loc_df, 'TSS', 'TSS') 
# 		G = label_nodes(G, loc_df, 'alt_TSS', 'alt_TSS') 
# 		G = label_nodes(G, loc_df, 'TES', 'TES')
# 		G = label_nodes(G, loc_df, 'alt_TES', 'alt_TES')
# 		G = label_nodes(G, loc_df, 'coord', 'coord')
# 		G = label_nodes(G, loc_df, 'annotated', 'annotated')
# 		G = label_edges(G, edge_df, 'annotated', 'annotated')
# 		G = label_edges(G, edge_df, 'strand', 'strand')
# 		G = label_edges(G, edge_df, 'edge_type', 'edge_type')

# 		return G

# #
# class PlottedGraph:
# 	def __init__(self, sg, args):

# 		self.G = sg.G.copy()
# 		self.loc_df = sg.loc_df.copy(deep=True)
# 		self.edge_df = sg.edge_df.copy(deep=True)
# 		self.t_df = sg.t_df.copy(deep=True)

# 		if args['combine']:
# 			nbps = self.find_nb_paths()
# 			self.agg_nb_nodes(nbps, args)

# 		# get positions/sizes of nodes, edges, and labels
# 		self.calc_pos_sizes()

# 		# determine edge plotting settings
# 		self.get_edge_plt_settings(args)
# 		self.get_node_plt_settings(args)

# 		# accessible fields:
# 		# G
# 		# loc_df
# 		# edge_df
# 		# t_df
# 		# edge_style
# 		# node_style
# 		# sub_node_style
# 		# pos
# 		# ordered_nodes
# 		# node_size
# 		# rad
# 		# label_size
# 		# edge_width

# 	# starting nodes:
# 	# 1. have only one outgoing edge AND
# 	# 2. have more than one incoming edge | is a TSS | parent is TES AND
# 	# 3. node is not already in a nbp
# 	def is_nbp_start(self,G,n,nbps):
# 		nbp_nodes = [n for p in nbps for n in p]
# 		if G.out_degree(n) == 1: # 1
# 			tes_parent = False
# 			parents = list(G.predecessors(n))
# 			if parents: 
# 				for parent in parents:
# 					if G.nodes[parent]['TES']:
# 						tes_parent = True
# 			if G.in_degree(n) >= 1 or G.nodes[n]['TSS'] or tes_parent: # 2
# 				if n not in nbp_nodes: # 3
# 					return True
# 		return False

# 	# end nodes:
# 	# 1. have more than one outgoing edge |
# 	# 2. have more than one incoming edge |
# 	# 3. is a TSS |
# 	# 4. parent node is a TES
# 	def is_nbp_end(self,G,n):
# 		if G.out_degree(n) > 1: # 1
# 			return True
# 		if G.in_degree(n) > 1: # 2
# 			return True
# 		if G.nodes[n]['TSS']: # 3 
# 			return True
# 		if G.in_degree(n) == 1: # 4
# 			parent = list(G.predecessors(n))[0]
# 			if G.nodes[parent]['TES']:
# 				return True
# 		return False

# 	# http://rosalind.info/problems/ba3m/ ty pavel 
# 	# modified to not group up TES/TSSs 
# 	# see is_nbp_start/end to see complete set of conditions
# 	def find_nb_paths(self):

# 		G = self.G
# 		t_df = self.t_df

# 		nbps = []
# 		for v in G.nodes:
# 			if self.is_nbp_start(G,v,nbps):
# 				if G.out_degree(v) > 0:
# 					for w in G.successors(v):
# 						nbp = [v]
# 						while not self.is_nbp_end(G,w):
# 							nbp.append(w)
# 							succ = list(G.successors(w))
# 							if succ: w = succ[0]
# 							else: break
# 						if len(nbp) > 1: nbps.append(nbp)

# 		# loop through each transcript path and each non-branching path,
# 		# replace nodes that will be aggregated 
# 		paths = t_df.path.tolist() # apparently this gives me a reference to the paths
# 								   # so when updating I'm directly modifying the t_df
# 		for path in paths: 
# 			for combined_index, nbp in enumerate(nbps):

# 				# get the nodes that need to be deleted from the path
# 				# in the order they need to be deleted 
# 				del_nodes = sorted([int(i) for i in list(set(path) & set(nbp))],
# 					key=lambda x: G.nodes[x]['coord'])

# 				# remove each node that is in a non-branching path
# 				for i, n in enumerate(del_nodes):
# 					if i == 0: 
# 						insertion_index = path.index(n)
# 					path.remove(n)
# 				if not del_nodes: insertion_index = -1

# 				# add aggregate node tuple instead
# 				if insertion_index != -1:
# 					path.insert(insertion_index, 'c{}'.format(combined_index))

# 		self.t_df = t_df
# 		return nbps

# 	# aggregate nonbranching nodes and add to graph. remove old nodes
# 	# update loc_df and edge_df to reflect these changes
# 	def agg_nb_nodes(self, paths, args):

# 		G = self.G
# 		mod_G = nx.DiGraph(self.G)
# 		loc_df = self.loc_df
# 		edge_df = self.edge_df

# 		combined_index = 0
# 		loc_df['agg_path'] = np.nan
# 		loc_df['combined'] = False
# 		loc_df['combined_types'] = np.nan
# 		edge_df.reset_index(inplace=True)

# 		# loop through each nonbranching path
# 		for path in paths: 

# 			path = tuple([int(n) for n in path])
# 			start = path[0]
# 			stop = path[-1]
# 			combined_node = 'c{}'.format(combined_index)

# 			# get the colors for each aggregate node
# 			if args['color_alt_nodes']:
# 				combined_types = [j[0] for j in sorted([i for i in 
# 				   [('alt_TSS',loc_df.loc[path,'alt_TSS'].tolist().count(True)),
# 				   ('alt_TES',loc_df.loc[path,'alt_TES'].tolist().count(True)),
# 				   ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
# 				   if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]
# 			else: 
# 				combined_types = [j[0] for j in sorted([i for i in 
# 				   [('TSS',loc_df.loc[path,'TSS'].tolist().count(True)),
# 				   ('TES',loc_df.loc[path,'TES'].tolist().count(True)),
# 				   ('internal',loc_df.loc[path,'internal'].tolist().count(True))]
# 				   if i[1] != 0], key=lambda x: x[1], reverse=True)][:2]

# 			# get all incoming edges to first node and
# 			# all outgoing edges to last node
# 			data_edges = mod_G.edges(data=True)
# 			if mod_G.in_degree(start) > 0: 
# 				in_nodes = list(mod_G.predecessors(start))
# 				in_edges = [(v1, combined_node) for v1 in in_nodes]

# 				index_map = {v1: i for i,v1 in enumerate(in_nodes)}
# 				in_edge_attrs = [md[1] for md in sorted([(v1,d)
# 								 for v1,v2,d in data_edges if v2==start],
# 								 key=lambda x: index_map[x[0]])]
# 			else: 
# 				in_edges = []
# 				in_edge_attrs = []

# 			if G.out_degree(stop) > 0: 
# 				out_nodes = list(mod_G.successors(stop))
# 				out_edges = [(combined_node, v2) for v2 in out_nodes]

# 				index_map = {v2: i for i,v2 in enumerate(out_nodes)}
# 				out_edge_attrs = [md[1] for md in sorted([(v2,d)
# 								  for v1,v2,d in data_edges if v1==stop],
# 								  key=lambda x: index_map[x[0]])]
# 			else: 
# 				out_edges = []
# 				out_edge_attrs = []

# 			edges = in_edges+out_edges
# 			edge_attrs = in_edge_attrs+out_edge_attrs

# 			# get the node makeup of the aggregate node, 
# 			# and remove the old nodes 
# 			coord = loc_df.loc[loc_df.vertex_id == start, 'coord'].tolist()[0] # use the first coordinate
# 			chrom = loc_df.loc[loc_df.vertex_id == start, 'chrom'].tolist()[0]
# 			strand = loc_df.loc[loc_df.vertex_id == start, 'strand'].tolist()[0]
# 			tss = alt_tss = tes = alt_tes = internal = False
# 			annotated = True
# 			for n in path:
# 				if mod_G.nodes[n]['TSS']: tss = True
# 				if mod_G.nodes[n]['TES']: tes = True
# 				if mod_G.nodes[n]['alt_TSS']: alt_tss = True
# 				if mod_G.nodes[n]['alt_TES']: alt_tes = True
# 				if mod_G.nodes[n]['internal']: internal = True
# 				if not mod_G.nodes[n]['annotated']: annotated = False

# 			# add the aggregate node and all associated edges to modified graph
# 			node_attrs = {'TSS': tss, 'TES': tes,
# 						  'alt_TSS': alt_tss, 'alt_TES': alt_tes,
# 						  'internal': internal, 'coord': coord,
# 						  'annotated': annotated}
# 			mod_G.add_node(combined_node)
# 			nx.set_node_attributes(mod_G, {combined_node: node_attrs})

# 			for edge, edge_attr in zip(edges, edge_attrs):
# 				mod_G.add_edge(edge[0], edge[1])
# 				nx.set_edge_attributes(mod_G, {edge: edge_attr})

# 			# also add aggregate node and all associated edges to loc_df, edge_df
# 			node_attrs.update({'chrom': chrom, 'strand': strand,
# 							   'vertex_id': combined_node,
# 							   'combined_types': combined_types,
# 							   'vertex_id_back': combined_node,
# 							   'agg_path': list(path), 'combined': True})
# 			for edge, edge_attr in zip(edges, edge_attrs):
# 				edge_attr.update({'edge_ID': edge, 'strand': strand,
# 								  'v1': edge[0], 'v2': edge[1]})
# 				edge_df = edge_df.append(edge_attr, ignore_index=True)
# 			loc_df = loc_df.append(node_attrs, ignore_index=True)

# 			# remove all old nodes from the graph
# 			nb_nodes = [n for n in path]
# 			for n in nb_nodes:
# 				mod_G.remove_node(n)

# 			# increment combined node index
# 			combined_index += 1

# 		edge_df.set_index('edge_ID', inplace=True)

# 		# finally, label all nodes that are now in the graph with their "combined" 
# 		# status
# 		mod_G = label_nodes(mod_G, loc_df, 'combined', 'combined')
# 		mod_G = label_nodes(mod_G, loc_df, 'combined_types', 'combined_types')
# 		mod_G = label_nodes(mod_G, loc_df, 'agg_path', 'agg_path')

# 		self.G = mod_G
# 		self.loc_df = loc_df
# 		self.edge_df = edge_df
		
# 	# calculates the positions and sizes of edges/nodes based on the 
# 	# number of nodes in the graph
# 	def calc_pos_sizes(self):

# 		G = self.G

# 		pos = defaultdict()

# 		# order nodes based on genomic coords
# 		ordered_nodes = [i[0] for i in sorted([[j, int(n['coord'])]
# 				for j,n in G.nodes(data=True)], key=lambda x: x[1])]

# 		# check if forward or reverse strand 
# 		# is the first node in the ordered list a TSS or TES? 
# 		if G.nodes[ordered_nodes[0]]['TES']:
# 			ordered_nodes.reverse()
		
# 		y_coord = 0
# 		x_coords = np.linspace(-1, 1, len(ordered_nodes))
# 		for x_coord, n in zip(x_coords, ordered_nodes): 
# 			pos[n] = [x_coord, y_coord]

# 		x = len(ordered_nodes)

# 		# area-related sizes need to be non-linearly scaled
# 		# calculated by fitting power series curve to handpicked sizes
# 		node_size = 19248*(x**-1.14) 
# 		label_size = 43.9*(x**-0.484)

# 		# linearly-related sizes
# 		rad = ((11/540)*x)+(73/540)
# 		edge_width = -x/18 + (121/18)

# 		self.pos = pos
# 		self.ordered_nodes = ordered_nodes
# 		self.node_size = node_size 
# 		self.rad = rad
# 		self.label_size = label_size
# 		self.edge_width = edge_width

# 	def save_fig(self, oname):
# 		plt.tight_layout()
# 		plt.savefig(oname, format='png', dpi=200)
# 		plt.clf()

# 	# plots the whole graph
# 	def plot_graph(self, args):

# 		# plotting stuff
# 		plt.figure(1, figsize=(14,10), frameon=False)
# 		plt.xlim(-1.05, 1.05)
# 		plt.ylim(-1.05, 1.05) 

# 		# get node/edge colors based on arguments
# 		self.get_node_colors(args)
# 		self.get_edge_colors(args)

# 		# plot edges, nodes, and labels
# 		self.plot_edges()

# 		# relabel nodes with strings for the purpose of labelling them
# 		relabel_map = {k:(k if type(k) == str else int(k)) for k in self.G.nodes}
# 		self.G = nx.relabel_nodes(self.G, relabel_map)
# 		nx.draw_networkx_labels(self.G, self.pos, font_size=self.label_size)

# 		self.plot_nodes()

# 	# plots a transcript path through a preexisiting 
# 	# PlottedGraph graph
# 	def plot_overlaid_path(self, path, args, oname):

# 		# first plot the preexisting graph in gray
# 		args['color_edges'] = args['color_nodes'] = args['color_alt_nodes'] = False
# 		self.plot_graph(args)

# 		# get fields from object
# 		G = self.G
		
# 		# get nodes and edges from the input path
# 		path_nodes = path
# 		path_edges = [(path[i], path[i+1]) for i in range(len(path)-1)]

# 		# create a subgraph based on the nodes and edges in this path
# 		path_pg = copy.deepcopy(self)
# 		nodes = [(n,d) for n,d in G.nodes(data=True)
# 	    		 if n in path_nodes]
# 		edges = [(v1,v2,d) for v1,v2,d in G.edges(data=True)
# 				 if (v1,v2) in path_edges]
	    
# 		G = nx.DiGraph()
# 		G.add_nodes_from(nodes)
# 		G.add_edges_from(edges)
# 		path_pg.G = G

# 		# plot the subgraph of this transcript path on top of 
# 		# preexisting graph
# 		args['color_edges'] = args['color_nodes'] = args['color_alt_nodes'] = True
# 		path_pg.plot_graph(args)
# 		path_pg.save_fig(oname)

# 	# plots a subgraph from a PlottedGraph object
# 	def plot_subgraph(self, nodelist, args):

# 		self.G = self.G.subgraph(nodelist)
# 		self.plot_graph()

# 	# plot nodes according to style dicts
# 	def plot_nodes(self):

# 		style_dict = self.node_style
# 		sub_style_dict = self.sub_node_style

# 		for n in self.G.nodes():
# 			nx.draw_networkx_nodes(self.G, self.pos,
# 				nodelist=[n],
# 				node_color=style_dict[n]['color'],
# 				node_size=style_dict[n]['size'],
# 				node_shape=style_dict[n]['shape'])
# 			if n in sub_style_dict.keys():
# 				nx.draw_networkx_nodes(self.G, self.pos,
# 					nodelist=[n],
# 					node_color=sub_style_dict[n]['color'],
# 					node_size=sub_style_dict[n]['size'],
# 					node_shape=sub_style_dict[n]['shape'])

# 	# updates the style dict with colors during actual plotting calls
# 	def get_node_colors(self, args):

# 		# get fields from obj that we care about
# 		G = self.G

# 		# plotting styles
# 		gray = '#999999'
# 		yellow = '#F0E442'
# 		blue = '#0072B2'
# 		light_blue = '#56B4E9'
# 		red = '#D55E00'
# 		orange = '#E69F00'
# 		color_dict = {'internal': yellow, 'TSS': blue,
# 					  'alt_TSS': light_blue, 'TES': red,
# 					  'alt_TES': orange}

# 		# assign nodes colors based on plotting settings
# 		for n in self.G.nodes():
# 			node_colors = defaultdict()
# 			sub_node_colors = defaultdict()
# 			node = G.nodes[n]

# 			# gray nodes
# 			if args['color_nodes'] == False:
# 				node_colors.update({'color': gray})

# 			# combined nodes
# 			elif args['combine'] and node['combined']:

# 				# all makeup of combined node is same type
# 				if len(node['combined_types']) == 1:
# 					color = color_dict[node['combined_types'][0]]
# 					node_colors.update({'color': color})

# 				# we need to plot more than one thing for combined node
# 				else:
# 					# color
# 					sub_color = color_dict[node['combined_types'][0]]
# 					color = color_dict[node['combined_types'][1]]
# 					node_colors.update({'color': color})
# 					sub_node_colors.update({'color': sub_color})

# 					# add combined sub-node color
# 					self.sub_node_style[n].update(sub_node_colors)

# 			# non-combined node
# 			else:
# 				node_colors.update({'color': color_dict['internal']})
# 				if node['TSS']: 
# 					node_colors.update({'color': color_dict['TSS']})
# 				if node['alt_TSS'] and args['color_alt_nodes']:
# 					node_colors.update({'color': color_dict['alt_TSS']})
# 				if node['TES']:
# 					node_colors.update({'color': color_dict['TES']})
# 				if node['alt_TES'] and args['color_alt_nodes']:
# 					node_colors.update({'color': color_dict['alt_TES']})

# 			# add color to normal node style
# 			self.node_style[n].update(node_colors)

# 	# returns a dictionary indexed by node ids of plotting styles
# 	def get_node_plt_settings(self, args):

# 		# get fields from object that we care about
# 		G = self.G
# 		pos = self.pos
# 		node_size = self.node_size

# 		# create a plotting settings dictionary for each node
# 		style_dict = defaultdict()
# 		sub_style_dict = defaultdict()
# 		for n in G.nodes():

# 			node = G.nodes[n]

# 			n_style_dict = defaultdict()
# 			n_style_dict.update({'size': node_size})

# 			# combined nodes 
# 			if args['combine'] and node['combined']:
# 				n_style_dict.update({'shape': 'H'})

# 				if len(node['combined_types']) > 1:
# 					# size
# 					sub_n_style_dict = defaultdict()
# 					sub_n_style_dict.update({'size': node_size/2})
# 					# shape
# 					sub_n_style_dict.update({'shape': 'H'})

# 					# add combined sub-node
# 					sub_style_dict.update({n: sub_n_style_dict})

# 			# non-combined node
# 			else:
# 				n_style_dict.update({'shape': None})

# 			# add normal node 
# 			style_dict.update({n: n_style_dict})


# 		self.node_style = style_dict
# 		self.sub_node_style = sub_style_dict

# 	# plot edges according to settings in style_dict
# 	def plot_edges(self):

# 		style_dict = self.edge_style

# 		for e in self.G.edges():
# 			nx.draw_networkx_edges(self.G, self.pos,
# 				edgelist=[e],
# 				width=style_dict[e]['width'],
# 				edge_color=style_dict[e]['color'],
# 				connectionstyle=style_dict[e]['connectionstyle'])

# 	# returns a dictionary indexed by edge ids of plotting styles 
# 	def get_edge_colors(self, args):

# 		G = self.G

# 		# colors
# 		pink = '#CC79A7'
# 		green = '#009E73'
# 		gray = '#999999'
# 		color_dict = {'intron': pink, 'exon': green}

# 		# assign a color for each edge
# 		edges = list(G.edges)
# 		for e in edges:
# 			e_style_dict = {}
# 			if args['color_edges']:
# 				if G.edges[e]['edge_type'] == 'intron': 
# 					e_style_dict.update({'color': color_dict['intron']})
# 				elif G.edges[e]['edge_type'] == 'exon':
# 					e_style_dict.update({'color': color_dict['exon']})
# 			else: 
# 				e_style_dict.update({'color': gray})

# 			self.edge_style[e].update(e_style_dict)

# 	# returns a dictionary indexed by edge ids of plotting styles 
# 	def get_edge_plt_settings(self, args):

# 		G = self.G
# 		ordered_nodes = self.ordered_nodes
# 		rad = self.rad
# 		edge_width = self.edge_width

# 		# plotting styles
# 		pos_style = 'arc3,rad={}'.format(rad)
# 		neg_style = 'arc3,rad=-{}'.format(rad)
# 		straight = None

# 		edges = list(G.edges)
# 		ordered_edges = [(n,m) for n,m in zip(ordered_nodes[:-1],ordered_nodes[1:])]
# 		style_dict = defaultdict()

# 		# straight edges
# 		if len(ordered_nodes) > 20:
# 			straight_edges = []
# 		else:
# 			straight_edges = list(set(edges)&set(ordered_edges))
		
# 		# create a plotting settings dictionary for each edge
# 		pos = 1
# 		neg = -1
# 		for e in edges: 
# 			e_style_dict = {}
# 			e_style_dict.update({'width': edge_width})

# 			# straight edge
# 			if e in straight_edges:
# 				e_style_dict.update({'connectionstyle': straight})
# 			else:
# 				if pos > 0:
# 					e_style_dict.update({'connectionstyle': pos_style})
# 					pos *= -1
# 					neg *= -1
# 				elif neg > 0: 
# 					e_style_dict.update({'connectionstyle': neg_style})
# 					pos *= -1 
# 					neg *= -1
# 			style_dict.update({e: e_style_dict})

# 		self.edge_style = style_dict

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






