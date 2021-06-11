import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd


# # TODO:  How to handle adding transcriptome after adding abundance? Need to fill
# nans with 0s and sort anndata object / recompute DPI / TPM? Except for also I just limit
# to transcripts that are already in the SwanGraph, so it should be OK (swangraph.add_abundance())

# TODO
# port all the tests from my test_dev ipynb here

# test different combinations of add_abundance
class TestAbundance(object):

    # add abundance - vanilla
    def test_add_abundance_1(self):
        sg = swan.SwanGraph()
        sg.add_annotation('files/test_full_annotation.gtf')
        sg.add_transcriptome('files/test_full.gtf')

        sg.add_abundance('files/test_ab_1.tsv')

        assert set(sg.datasets) == set(['dataset1', 'dataset2'])

        print(sg.t_df.index)
        print(sg.adata.var.index)
        assert sg.t_df.index.tolist() == sg.adata.var.index.tolist()


        # assert 1 == 0
        ctrl_X = np.transpose(np.array([[5,5],[10,0],[0,10],[10,10],[5,5],[0,0]])).astype(np.float)
        print('test')
        print(sg.adata.X)
        df = pd.DataFrame(index=sg.adata.obs.index, columns=sg.adata.var.index, data=sg.adata.layers['counts'])
        print(df.head())
        print('control')
        print(ctrl_X)
        assert np.array_equal(sg.adata.X, ctrl_X)

        ctrl_tpm = [[166666.6667, 166666.6667],
                    [333333.3333, 0],
                    [0, 333333.3333],
                    [333333.3333, 333333.3333],
                    [166666.6667, 166666.6667],
                    [0, 0]]
        ctrl_tpm = np.transpose(np.array(ctrl_tpm)).astype(np.float)
        print('tpm')
        print('test')
        print(sg.adata.layers['tpm'])
        print('control')
        print(ctrl_tpm)

        df = pd.DataFrame(index=sg.adata.obs.index, columns=sg.adata.var.index, data=sg.adata.layers['tpm'])
        print(df.head())
        ctrl_tpm = np.around(ctrl_tpm)
        test_tpm = np.around(sg.adata.layers['tpm'])
        print('tpm')
        print('test')
        print(test_tpm)
        print('control')
        print(ctrl_tpm)
        print(np.array_equal(ctrl_tpm, test_tpm))
        assert np.array_equal(ctrl_tpm, test_tpm)

# test lower-level functions

# todo - calc_tpm, calc_pi
class TestLowLevel(object):
    pass
