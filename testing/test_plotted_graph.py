import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
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

	def test_assign_combined_dataset(self):
		G = nx.Graph()
		for i in range(6):
			G.add_node(i)
		paths = [[0,1], [2,3], [4,5]]
		for p in paths:
			nx.add_path(G, p)

		# which datasets are each set of nodes from
		attrs = {0: {'dataset_1': True, 'dataset_2': False},
				 1: {'dataset_1': True, 'dataset_2': False},
				 2: {'dataset_1': False, 'dataset_2': True},
				 3: {'dataset_1': False, 'dataset_2': True},
				 4: {'dataset_1': True, 'dataset_2': True},
				 5: {'dataset_1': True, 'dataset_2': True},}
		nx.set_node_attributes(G, attrs)

		# testing
		control = {(0,1): {'dataset_1': True, 'dataset_2': False},
				   (2,3): {'dataset_1': False, 'dataset_2': True},
				   (4,5): {'dataset_1': True, 'dataset_2': True}}
		test = {}
		for p in paths:
			node_attrs = {}
			node_attrs = pg.assign_combined_datasets(G, p, node_attrs)
			test.update({tuple(p): node_attrs})

		print('control')
		print(control)

		print('test')
		print(test)
		
		for k in test.keys():
			assert test[k] == control[k]



