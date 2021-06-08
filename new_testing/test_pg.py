import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd
import anndata

###########################################################################
###################### Utilities in graph.py ##############################
###########################################################################
class TestGeneSummary(object):

    # test minus strand gene summary
    def test_gene_summary_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_graph('test2_gid')
        sg.pg.edge_df.drop(['curve', 'line'], axis=1, inplace=True)

        ctrl_edge_df = pd.read_csv('files/plot_minus_gene_summary_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v1_old', 'v2_old', 'line'], axis=1, inplace=True)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        ctrl_loc_df = pd.read_csv('files/plot_minus_gene_summary_loc_df.tsv', sep='\t')
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df['edge_color'] = None
        ctrl_loc_df['linewidth'] = None
        ctrl_loc_df.drop(['edgecolor', 'linewidth'], axis=1, inplace=True)
        sg.pg.loc_df.drop(['edgecolor', 'linewidth'], axis=1, inplace=True)

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)






def check_dfs(loc_df, ctrl_loc_df,
              edge_df, ctrl_edge_df):

    # first order to make them comparable
    # sort all values by their IDs
    loc_df.sort_index(inplace=True)
    edge_df.sort_index(inplace=True)
    ctrl_loc_df.sort_index(inplace=True)
    ctrl_edge_df.sort_index(inplace=True)

    # and order columns the same way
    ctrl_loc_df = ctrl_loc_df[loc_df.columns]
    ctrl_edge_df = ctrl_edge_df[edge_df.columns]

    print('test')
    print(loc_df)
    print('control')
    print(ctrl_loc_df)
    print(ctrl_loc_df == loc_df)
    assert (loc_df == ctrl_loc_df).all(axis=0).all()

    print('test')
    print(edge_df)
    print('control')
    print(ctrl_edge_df)
    assert (edge_df == ctrl_edge_df).all(axis=0).all()
