import networkx as nx
import argparse
import pandas as pd
import pickle
import os
import sys
import SpliceGraph as sg
from utils import *

def get_args():

	desc = 'Loads a GTF into a graph'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-gtf', '-g', dest='gtf',
		help='GTF annotation to load into graph')
	parser.add_argument('--o', dest='oname', default=None,
		help='Output directory/prefix to store graph')

	args = parser.parse_args()
	return args

def make_oname(oname):

	if not oname: oname = os.getcwd()+'/splicing_graph.p'
	else: oname = oname+'_splicing_graph.p'

	return oname

def main():

	args = get_args()
	gtffile = args.gtf
	oname = make_oname(args.oname)

	graph = sg.SpliceGraph(gtf=gtffile)
	
	# save to output file
	with open(oname, 'wb') as ofile:
		pickle.dump(graph, ofile)

if __name__ == '__main__': main()