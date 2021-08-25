import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd

###########################################################################
###################### Related to data analysis ############################
###########################################################################
class TestSGAnalysis(object):

# tests find_ir_genes, find_es_genes, get_die_genes, get_die_gene_table, test_gene

# TODO - add update_ids() call before running both find_ir/es_genes to make sure
# the subgraph shit works -
# TODO - check to make sure this works with + AND - strands because the
# loc_df locations no longer have strand associated with them from get_ordered_id_map!!!!

    # test get_die_gene_table - gene that meets rc
    def test_get_die_table_6(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts']
        data = [[1, 1, .5, .5, 20, 10, 10],
                [2, 1, .5, .5, 20, 10, 10]]
        df = pd.DataFrame(data=data, columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=15)
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts', 'dpi']
        data = [[1, 1, .5, .5, 20, 10, 10, 0],
                [2, 1, .5, .5, 20, 10, 10, 0]]
        ctrl = pd.DataFrame(data=data, columns=columns)
        ctrl.set_index('tid', inplace=True)
        print(df)
        print(ctrl)
        print(ctrl == df)
        assert (ctrl == df).all(axis=0).all()

    # test get_die_gene_table - gene with 11 > n isoforms > 1
    def test_get_die_table_5(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts']
        data = [[1, 1, .5, .5, 20, 10, 10],
                [2, 1, .5, .5, 20, 10, 10]]
        df = pd.DataFrame(data=data, columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts', 'dpi']
        data = [[1, 1, .5, .5, 20, 10, 10, 0],
                [2, 1, .5, .5, 20, 10, 10, 0]]
        ctrl = pd.DataFrame(data=data, columns=columns)
        ctrl.set_index('tid', inplace=True)
        print(df)
        print(ctrl)
        print(ctrl == df)
        assert (ctrl == df).all(axis=0).all()

    # test get_die_gene_table - gene with no expressed isoforms
    def test_get_die_table_4(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts']
        data = [[1, 1, 0, 0, 0, 0, 0],
                [2, 1, 0, 0, 0, 0, 0]]
        df = pd.DataFrame(data=data, columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        assert df == None

    # test get_die_gene_table - gene doesn't have enough reads (rc thresh)
    def test_get_die_table_3(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts']
        data = [[1, 1, .5, .5, 20, 10, 10],
                [2, 1, .5, .5, 20, 10, 10]]
        df = pd.DataFrame(data=data, columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df

        df = swan.get_die_gene_table(df, conditions, rc=30)
        assert df == None

    # test get_die_gene_table - limit to genes with >1 iso
    def test_get_die_table_2(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts', 'dpi']
        data = [1, 1, 1, 1, 30, 15, 15, 0]
        df = pd.DataFrame(data=[data], columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        assert df == None

    # test get_die_gene_table_1 - limit to isoforms with > 0 exp in at least
    # one condition, aggregate last n-11 isos
    def test_get_die_table_1(self):
        conditions = ['cond1', 'cond2']
        df = pd.read_csv('files/test_pi_1.tsv', sep='\t')
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        df.tid = df.tid.astype('str')
        df.set_index('tid', inplace=True)
        ctrl = pd.read_csv('files/test_pi_1_reference.tsv', sep='\t')
        ctrl.tid = ctrl.tid.astype('str')
        ctrl.set_index('tid', inplace=True)
        print(ctrl)
        print(df)
        ctrl = ctrl[df.columns]
        ctrl.dpi = ctrl.dpi.astype('float').round()
        df.dpi = df.dpi.astype('float').round()
        print(ctrl == df)
        assert (ctrl == df).all(axis=0).all()

        # still need better way to test this rip
    # test test_gene - 2 up 2 down
    def test_test_gene_3(self):
        conditions = ['cond1', 'cond2']
        df = pd.read_csv('files/test_pi_1.tsv', sep='\t')
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        ctrl_gene_dpi = 0.41

        pval, gene_dpi = swan.test_gene(df, conditions)
        gene_dpi = np.round(gene_dpi, decimals=2)
        ctrl_gene_dpi = np.round(ctrl_gene_dpi, decimals=2)
        print(gene_dpi)
        print(ctrl_gene_dpi)
        assert gene_dpi == ctrl_gene_dpi

    # test test_gene - <2 up, <2 down
    def test_test_gene_2(self):
        conditions = ['cond1', 'cond2']
        df = pd.read_csv('files/test_pi_1.tsv', sep='\t')
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        ctrl_gene_dpi = 0.22
        tids = [3, 1]
        df = df.loc[df.tid.isin(tids)]
        pval, gene_dpi = swan.test_gene(df, conditions)
        assert gene_dpi == ctrl_gene_dpi

    # test test_gene - all 0
    def test_test_gene_1(self):
        columns = ['tid', 'gid', 'cond1', 'cond2', 'total_counts', 'cond1_counts', 'cond2_counts']
        data = [[1, 1, .5, .5, 20, 10, 10],
                [2, 1, .5, .5, 20, 10, 10]]
        df = pd.DataFrame(data=data, columns=columns)
        conditions = ['cond1', 'cond2']
        df.set_index('tid', inplace=True)
        df = swan.get_die_gene_table(df, conditions, rc=0)
        ctrl_gene_dpi = 0
        pval, gene_dpi = swan.test_gene(df, conditions)
        assert gene_dpi == ctrl_gene_dpi

    # test get_die_genes - obs col doesn't exist
    def test_get_die_genes_5(self):
        sg = get_die_test_sg()
        print(sg.adata.obs)
        obs_col = 'condition'
        obs_conditions = ['PB65_B017', 'PB65_B018']
        id_col = 'tid'
        with pytest.raises(Exception) as e:
            test = sg.get_die_genes(obs_col=obs_col, obs_conditions=obs_conditions, rc_thresh=1)
        assert 'Metadata column' in str(e.value)

    # tests get_die_genes - b/w  dataset conditions
    def test_get_die_genes_4(self):
        sg = get_die_test_sg()
        obs_col = 'dataset'
        obs_conditions = ['PB65_B017', 'PB65_B018']
        id_col = 'tid'
        test = sg.get_die_genes(obs_col=obs_col, obs_conditions=obs_conditions, rc_thresh=1)
        # don't test p vals cause that's tough
        test.drop(['p_val', 'adj_p_val'], axis=1, inplace=True)
        ctrl = pd.read_csv('files/chr11_and_Tcf3_PB65_B017_B018_dpi.tsv', sep='\t')
        test.dpi = test.dpi.round().astype('float')
        ctrl.dpi = ctrl.dpi.round().astype('float')
        print(test)
        print(ctrl)
        print(test == ctrl)
        assert test.equals(ctrl)

    # tests get_die_genes - b/w non dataset conditions
    def test_get_die_genes_3(self):
        sg = get_die_test_sg()
        obs_col = 'cluster'
        id_col = 'tid'
        test = sg.get_die_genes(obs_col=obs_col, rc_thresh=1)

        # don't test p vals cause that's tough
        test.drop(['p_val', 'adj_p_val'], axis=1, inplace=True)
        ctrl = pd.read_csv('files/chr11_and_Tcf3_cluster_dpi.tsv', sep='\t')
        test.dpi = test.dpi.round()
        ctrl.dpi = ctrl.dpi.round()
        print(ctrl)
        print(test)
        print(ctrl == test)
        assert ctrl.equals(test)

    # tests get_die_genes - obs_col has more than 2 values and obs_conditions
    # not provided
    def test_get_die_genes_2(self):
        sg = get_die_test_sg()
        obs_col = 'dataset'
        id_col = 'tid'
        with pytest.raises(Exception) as e:
            test = sg.get_die_genes(obs_col=obs_col, rc_thresh=1)
        assert 'Must provide obs_conditions' in str(e.value)

    # tests get_die_genes - obs_col has more than 2 values and obs_conditions
    # does not have 2 values
    def test_get_die_genes_1(self):
        sg = get_die_test_sg()
        obs_col = 'dataset'
        id_col = 'tid'
        obs_conditions = ['D12']
        with pytest.raises(Exception) as e:
            test = sg.get_die_genes(obs_col=obs_col, obs_conditions=obs_conditions, rc_thresh=1)
        assert 'exactly 2 values' in str(e.value)

    # tests find_es_genes - requires edges to be in order
    def test_find_es_genes(self):
        sg = swan.SwanGraph()
        sg.annotation = True

        # t_df
        data = [[0, [0,1,2,3,4], True, 'g1', 't1'],
                [1, [0,5,4], False, 'g1', 't2']]
        cols = ['vertex_id', 'path', 'annotation', 'gid', 'tid']
        sg.t_df = pd.DataFrame(data=data, columns=cols)

        # edge
        data = [[0, 'exon', True, 0, 1],
                [1, 'intron', True, 1, 2],
                [2, 'exon', True, 2, 3],
                [3, 'intron', True, 3, 4],
                [4, 'exon', True, 4, 5],
                [5, 'intron', False, 1, 4]]
        cols = ['edge_id', 'edge_type', 'annotation', 'v1', 'v2']
        sg.edge_df = pd.DataFrame(data=data, columns=cols)

        # loc
        data = [0,1,2,3,4,5]
        cols = ['vertex_id']
        sg.loc_df = pd.DataFrame(data=data, columns=cols)

        ctrl_gids = ['g1']
        ctrl_tids = ['t2']
        ctrl_edges = [5]

        sg.get_loc_path()
        sg.create_graph_from_dfs()
        gids, tids, edges = sg.find_es_genes()

        print('edges')
        print('control')
        print(ctrl_edges)
        print('test')
        print(edges)
        assert set(ctrl_edges) == set(edges)

        print('transcripts')
        print('control')
        print(ctrl_tids)
        print('test')
        print(tids)
        assert set(ctrl_tids) == set(tids)

        print('genes')
        print('control')
        print(ctrl_gids)
        print('test')
        print(gids)
        assert set(ctrl_gids) == set(gids)

    # tests find_ir_genes - requires edges to be in order...
    # also requires the swangraph
    def test_find_ir_genes(self):
        sg = swan.SwanGraph()
        sg.annotation = True

        # t_df
        data = [[0, [0,1,2], True, 'g1', 't1'],
                [1, [3], False, 'g1', 't2']]
        cols = ['vertex_id', 'path', 'annotation', 'gid', 'tid']
        sg.t_df = pd.DataFrame(data=data, columns=cols)

        # edge
        data = [[0, 'exon', True, 0, 1],
                [1, 'intron', True, 1, 2],
                [2, 'exon', True, 2, 3],
                [3, 'exon', False, 0, 3]]
        cols = ['edge_id', 'edge_type', 'annotation', 'v1', 'v2']
        sg.edge_df = pd.DataFrame(data=data, columns=cols)

        # loc
        data = [0,1,2,3]
        cols = ['vertex_id']
        sg.loc_df = pd.DataFrame(data=data, columns=cols)

        ctrl_edges = [3]
        ctrl_tids = ['t2']
        ctrl_gids = ['g1']

        sg.get_loc_path()
        sg.create_graph_from_dfs()
        gids, tids, edges = sg.find_ir_genes()

        print('edges')
        print('control')
        print(ctrl_edges)
        print('test')
        print(edges)
        assert set(ctrl_edges) == set(edges)

        print('transcripts')
        print('control')
        print(ctrl_tids)
        print('test')
        print(tids)
        assert set(ctrl_tids) == set(tids)

        print('genes')
        print('control')
        print(ctrl_gids)
        print('test')
        print(gids)
        assert set(ctrl_gids) == set(gids)

def get_die_test_sg():
    sg = swan.SwanGraph()
    db = 'files/chr11_and_Tcf3_no_gname.db'
    sg.add_transcriptome(db, include_isms=True)
    ab = 'files/chr11_and_Tcf3_talon_abundance_PB65.tsv'
    sg.add_abundance(ab)
    ab = 'files/chr11_and_Tcf3_talon_abundance_D12.tsv'
    sg.add_abundance(ab)
    meta = 'files/chr11_and_Tcf3_metadata.tsv'
    sg.add_metadata(meta)
    return sg
