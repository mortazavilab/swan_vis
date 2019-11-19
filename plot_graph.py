import networkx as nx
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from collections import defaultdict
from utils import *

def get_args():

	desc = 'Loads a GTF into a graph'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-g', dest='graph_file',
		help='Pickled graph file (.p)')
	parser.add_argument('--o', dest='oname', default=None,
		help='Output directory/prefix to store figure')
	parser.add_argument('--color_introns', dest='color_introns',
		default=False, action='store_true',
		help='Use this flag for distinguishing exon/intron colors')
	parser.add_argument('--color_TSS', dest='color_TSS',
		default=False, action='store_true',
		help='Use this flag for distinguishing TSS colors')
	parser.add_argument('--color_TES', dest='color_TES',
		default=False, action='store_true',
		help='Use this flag for distinguishing TES colors')
	parser.add_argument('--color_alt_TSS', dest='color_alt_TSS',
		default=False, action='store_true',
		help='Use this flag for distinguishing alt. TSS colors')
	parser.add_argument('--color_alt_TES', dest='color_alt_TES',
		default=False, action='store_true',
		help='Use this flag for distinguishing alt. TES colors')
	parser.add_argument('--agg_nodes', dest='agg_nodes', 
		default=False, action='store_true')

	args = parser.parse_args()

	# if we're coloring alternative TSS/TES's, we should
	# probably just color all of them...
	if not args.color_TSS and args.color_alt_TSS:
		args.color_TSS = True
	if not args.color_TES and args.color_alt_TES: 
		args.color_TES = True

	return args

def make_oname(args):

	oname = args.oname

	if not oname:
		oname = '{}/'.format(os.getcwd())
	else: 
		oname = '{}'.format(oname)

	if args.color_introns: oname += '_exons_introns'

	if args.color_TSS and not args.color_alt_TSS: oname += '_TSS'
	if args.color_TES and not args.color_alt_TES: oname += '_TES'

	if args.color_alt_TSS: oname += '_alt_TSS'
	if args.color_alt_TES: oname += '_alt_TES'

	if args.agg_nodes: oname += '_combined'
	
	oname += '.png'
	return oname

def main():

	args = get_args()

	# get output figure name
	oname = make_oname(args)

	# load in pickled data
	with open(args.graph_file, 'rb') as pfile:
		sg = pickle.load(pfile)

	temp = defaultdict()
	temp['color_edges'] = args.color_introns
	temp['color_nodes'] = args.color_TSS
	temp['color_alt_nodes'] = args.color_alt_TSS
	temp['combine'] = args.agg_nodes
	temp['ann'] = False

	pg = utils.PlottedGraph(sg, temp)
	pg.plot_graph(oname)

if __name__ == '__main__': main()