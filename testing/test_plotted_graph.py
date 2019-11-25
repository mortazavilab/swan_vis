import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 
import networkx as nx

class TestPlottedGraph(object):
	def test_is_unique_to_dataset(self):
		G = nx.Graph()
		G.add_node(0)
		attrs = {0: {'coord': 3, 'internal': True, 'TES': False,
				  'dataset_1': False, 'dataset_lizard': True,
				  'dataset_-1': False}}
		nx.set_node_attributes(G, attrs)

		fields = sg.get_dataset_fields(G)
		data = G.nodes(data=True)[0]
		print(data)

		assert pg.is_unique_to_dataset(data,'dataset_1',fields) == False
		assert pg.is_unique_to_dataset(data,'dataset_lizard',fields) == True
		assert pg.is_unique_to_dataset(data,'dataset_-1',fields) == False

		attrs[0].update({'dataset_-1': True})
		nx.set_node_attributes(G, attrs)
		assert pg.is_unique_to_dataset(data,'dataset_1',fields) == False
		assert pg.is_unique_to_dataset(data,'dataset_lizard',fields) == False
		assert pg.is_unique_to_dataset(data,'dataset_-1',fields) == False


