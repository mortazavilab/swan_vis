import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd
import anndata

###########################################################################
###################### Utilities in utils.py ##############################
###########################################################################
class TestUtils(object):

    # test reset_dupe_index
    def test_reset_dupe_index(self):
        data = [[0,1],[2,1]]
        cols = ['col1', 'col2']
        df = pd.DataFrame(data=data, columns=cols)
        df = swan.create_dupe_index(df, 'col2')
        df = swan.set_dupe_index(df, 'col2')
        df = swan.reset_dupe_index(df, 'col2')

        assert 'col2_back' in df.columns
        assert [1,1] == df.col2_back.tolist()

    # test set_dupe_index
    def test_set_dupe_index(self):
        data = [[0,1],[2,1]]
        cols = ['col1', 'col2']
        df = pd.DataFrame(data=data, columns=cols)
        df = swan.create_dupe_index(df, 'col2')
        df = swan.set_dupe_index(df, 'col2')

        assert df.index.name == 'col2'
        assert [1,1] == df.index.tolist()
        assert [1,1] == df.col2.tolist()



    # test create_dupe_index
    def test_create_dupe_index(self):
        data = [[0,1],[2,1]]
        cols = ['col1', 'col2']
        df = pd.DataFrame(data=data, columns=cols)
        df = swan.create_dupe_index(df, 'col2')

        assert 'col2_back' in df.columns
        assert [1,1] == df.col2_back.tolist()
