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

    # done: test_new_gene, calc_pos_sizes

    # todo: calc_edge_curves

    # test calc_edge_curves - >= 20 node_size
    def test_calc_edge_curves_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        gid = 'test2_gid'
        sg.pg.new_gene(sg, gid)

        sg.pg.calc_pos_sizes()

        data = [0 for i in range(25)]
        cols = ['vertex_id']
        sg.pg.loc_df = pd.DataFrame(data=data, columns=cols)

        sg.pg.calc_edge_curves()

        # make sure there are no straight edges
        for ind, entry in sg.pg.edge_df.iterrows():
            print(ind)
            print(entry)
            assert 'arc3' in entry.curve

    # test calc_edge_curves - < 20 node_size
    def test_calc_edge_curves_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        gid = 'test2_gid'
        sg.pg.new_gene(sg, gid)
        sg.pg.calc_pos_sizes()
        sg.pg.calc_edge_curves()

        straight_edges = [9,8,12,13,6,5]
        temp = sg.pg.edge_df.loc[straight_edges]
        for ind, entry in temp.iterrows():
            print(ind)
            print(entry)
            assert not entry.curve

        curve_edges = set(sg.pg.edge_df.edge_id.tolist())-set(straight_edges)
        temp = sg.pg.edge_df.loc[curve_edges]
        for ind, entry in temp.iterrows():
            print(ind)
            print(entry)
            assert 'arc3' in entry.curve

    # test calc_pos_sizes - doesn't check any of the individual calculations
    # cause that would be very tedious
    def test_calc_pos_sizes(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        gid = 'test2_gid'
        sg.pg.new_gene(sg, gid)
        sg.pg.calc_pos_sizes()

        # make sure these are going in order
        print(sg.pg.pos)
        curr_x = -100
        for key, item in sg.pg.pos.items():
            x = item[0]
            assert x > curr_x
            curr_x = x

        assert 'pos' in sg.pg.loc_df.columns
        assert sg.pg.node_size
        assert sg.pg.label_size
        assert sg.pg.rad_scale == 0.32
        assert sg.pg.edge_width
        assert sg.pg.line_width

    # test new_gene - minus strand
    def test_new_gene_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        gid = 'test2_gid'
        sg.pg.new_gene(sg, gid)

        # loc_df
        data = [[0, 'chr2', 100, False, True, False],
                [1, 'chr2', 80, True, False, False],
                [2, 'chr2', 75, True, False, False],
                [3, 'chr2', 65, True, False, False],
                [4, 'chr2', 60, True, False, False],
                [5, 'chr2', 50, True, False, True],
                [6, 'chr2', 45, False, False, True]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        # edge_df
        data = [[9, '-', 'exon', 5, 6],
                [8, '-', 'intron', 4, 5],
                [12, '-', 'exon', 4, 5],
                [14, '-', 'intron', 3, 5],
                [7, '-', 'exon', 2, 4],
                [13, '-', 'exon', 2, 3],
                [11, '-', 'intron', 1, 4],
                [6, '-', 'intron', 1, 2],
                [5, '-', 'exon', 0, 1]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

        assert sg.pg.g_min == 45
        assert sg.pg.g_max == 100
        assert sg.pg.strand == '-'

    # # test minus strand gene summary
    # def test_gene_summary_1(self):
    #     sg = swan.SwanGraph()
    #     sg.add_transcriptome('files/test_full.gtf')
    #     sg.plot_graph('test2_gid')
    #     sg.pg.edge_df.drop(['curve', 'line'], axis=1, inplace=True)
    #
    #     ctrl_edge_df = pd.read_csv('files/plot_minus_gene_summary_edge_df.tsv', sep='\t')
    #     ctrl_edge_df.drop(['v1_old', 'v2_old', 'line'], axis=1, inplace=True)
    #     ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
    #     ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')
    #
    #     ctrl_loc_df = pd.read_csv('files/plot_minus_gene_summary_loc_df.tsv', sep='\t')
    #     ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
    #     ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')
    #     ctrl_loc_df['edge_color'] = None
    #     ctrl_loc_df['linewidth'] = None
    #     ctrl_loc_df.drop(['edgecolor', 'linewidth'], axis=1, inplace=True)
    #     sg.pg.loc_df.drop(['edgecolor', 'linewidth'], axis=1, inplace=True)
    #
    #     check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

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
