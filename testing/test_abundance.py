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

    # add abundance - talon format
    def test_add_abundance_3(self):
        sg = swan.SwanGraph()
        sg.add_annotation("files/test_full_annotation.gtf")
        sg.add_transcriptome("files/test_full.gtf")

        sg.add_abundance("files/test_ab_talon_1.tsv")

        assert set(sg.datasets) == set(["dataset1", "dataset2"])

        ctrl_tids = ["test1", "test2", "test3", "test4", "test5"]
        print(ctrl_tids)
        print(sg.adata.var.index)
        assert set(ctrl_tids) == set(sg.adata.var.index.tolist())

        ctrl_X = np.transpose(
            np.array([[5, 5], [10, 0], [0, 10], [10, 10], [5, 5]])
        ).astype(np.float)
        print("test")
        print(sg.adata.X)
        df = pd.DataFrame(
            index=sg.adata.obs.index,
            columns=sg.adata.var.index,
            data=sg.adata.layers["counts"].toarray(),
        )
        print(df.head())
        print("control")
        print(ctrl_X)
        assert np.array_equal(sg.adata.X.toarray(), ctrl_X)

        ctrl_tpm = [
            [166666.6667, 166666.6667],
            [333333.3333, 0],
            [0, 333333.3333],
            [333333.3333, 333333.3333],
            [166666.6667, 166666.6667],
        ]

        ctrl_tpm = np.transpose(np.array(ctrl_tpm)).astype(np.float)
        ctrl_tpm = np.around(ctrl_tpm)
        test_tpm = np.around(sg.adata.layers["tpm"].toarray())
        print("tpm")
        print("test")
        print(test_tpm)
        print("control")
        print(ctrl_tpm)
        print(np.array_equal(ctrl_tpm, test_tpm))
        assert np.array_equal(ctrl_tpm, test_tpm)

        ctrl_pi = [
            [100, 100],
            [66.6666667, 0],
            [0, 66.6666667],
            [100, 100],
            [33.3333333, 33.3333333],
        ]
        ctrl_pi = np.transpose(np.array(ctrl_pi)).astype(np.float)
        ctrl_pi = np.around(ctrl_pi)
        test_pi = np.around(sg.adata.layers["pi"].toarray())
        print("pi")
        print("test")
        print(test_pi)
        print("ctrl")
        print(ctrl_pi)
        assert np.array_equal(ctrl_pi, test_pi)

        # test tss stuff
        # entries
        ctrl = ["test1_gid_1", "test2_gid_1", "test4_gid_1"]
        test = sg.tss_adata.var.index.tolist()
        print("tss index")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert set(ctrl) == set(test)

        # counts
        ctrl = [[5, 5], [15, 15], [10, 10]]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        test = sg.tss_adata.layers["counts"].toarray()
        print("tss counts")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # tpm
        ctrl = [
            [166666.6667, 166666.6667],
            [500000, 500000],
            [333333.3333, 333333.3333],
        ]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        ctrl = np.around(ctrl)
        test = np.around(sg.tss_adata.layers["tpm"].toarray())
        print("tss tpm")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # pi
        ctrl = [[100, 100], [100, 100], [100, 100]]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        ctrl = np.around(ctrl)
        test = np.around(sg.tss_adata.layers["pi"].toarray())
        print("tss pi")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # test tes stuff
        # entries
        ctrl = ["test1_gid_1", "test2_gid_1", "test2_gid_2", "test4_gid_1"]
        test = sg.tes_adata.var.index.tolist()
        print("tes index")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert set(ctrl) == set(test)

        # counts
        ctrl = [[5, 5], [10, 10], [5, 5], [10, 10]]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        test = sg.tes_adata.layers["counts"].toarray()
        print("tes counts")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # tpm
        ctrl = [
            [166666.6667, 166666.6667],
            [333333.3333, 333333.3333],
            [166666.6667, 166666.6667],
            [333333.3333, 333333.3333],
        ]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        ctrl = np.around(ctrl)
        test = np.around(sg.tes_adata.layers["tpm"].toarray())
        print("tes tpm")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # pi
        ctrl = [
            [100, 100],
            [66.66666667, 66.66666667],
            [33.33333333, 33.33333333],
            [100, 100],
        ]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        ctrl = np.around(ctrl)
        test = np.around(sg.tes_adata.layers["pi"].toarray())
        print("tes pi")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # test edge stuff
        # entries
        ctrl = [0, 1, 2, 3, 4, 10, 9, 8, 12, 15, 7, 14, 11, 6, 5]
        ctrl = [str(c) for c in ctrl]
        test = sg.edge_adata.var.index.tolist()
        print("edge index")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert set(ctrl) == set(test)

        # counts
        ctrl = [
            [5, 5],
            [5, 5],
            [5, 5],
            [5, 5],
            [5, 5],
            [10, 10],
            [10, 10],
            [10, 0],
            [5, 5],
            [0, 10],
            [10, 0],
            [0, 10],
            [5, 5],
            [10, 10],
            [15, 15],
        ]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        test = sg.edge_adata.layers["counts"].toarray()
        print("edge counts")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # tpm
        ctrl = [
            [50000, 50000],
            [50000, 50000],
            [50000, 50000],
            [50000, 50000],
            [50000, 50000],
            [100000, 100000],
            [100000, 100000],
            [100000, 0],
            [50000, 50000],
            [0, 100000],
            [100000, 0],
            [0, 100000],
            [50000, 50000],
            [100000, 100000],
            [150000, 150000],
        ]
        ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        ctrl = np.around(ctrl)
        test = np.around(sg.edge_adata.layers["tpm"].toarray())
        print("edge tpm")
        print("ctrl")
        print(ctrl)
        print("test")
        print(test)
        assert np.array_equal(ctrl, test)

        # # pi
        # ctrl = [[5,5],[5,5],[5,5],[5,5],
        #         [5,5],[1,1],[1,1],[1,0],
        #         [5,5],[0,1],[1,0],[0,1],
        #         [5,5],[1,1],[15,15]]
        # ctrl = np.transpose(np.array(ctrl)).astype(np.float)
        # ctrl = np.around(ctrl)
        # test = np.around(sg.edge_adata.layers['pi'].toarray())
        # print('tes pi')
        # print('ctrl')
        # print(ctrl)
        # print('test')
        # print(test)
        # assert np.array_equal(ctrl, test)

    # add abundance - one after another
    def test_add_abundance_2(self):
        sg = swan.SwanGraph()
        sg.add_annotation("files/test_full_annotation.gtf")
        sg.add_transcriptome("files/test_full.gtf")

        sg.add_abundance("files/test_ab_dataset1.tsv")
        sg.add_abundance("files/test_ab_dataset2.tsv")

        assert set(sg.datasets) == set(["dataset1", "dataset2"])

        ctrl_tids = ["test1", "test2", "test3", "test4", "test5"]
        print(ctrl_tids)
        print(sg.adata.var.index)
        assert set(ctrl_tids) == set(sg.adata.var.index.tolist())

        # need to remove 0s
        ctrl_X = np.transpose(
            np.array([[5, 5], [10, 0], [0, 10], [10, 10], [5, 5]])
        ).astype(np.float)
        print("test")
        print(sg.adata.X)
        df = pd.DataFrame(
            index=sg.adata.obs.index,
            columns=sg.adata.var.index,
            data=sg.adata.layers["counts"].toarray(),
        )
        print(df.head())
        print("control")
        print(ctrl_X)
        assert np.array_equal(sg.adata.X.toarray(), ctrl_X)

        ctrl_tpm = [
            [166666.6667, 166666.6667],
            [333333.3333, 0],
            [0, 333333.3333],
            [333333.3333, 333333.3333],
            [166666.6667, 166666.6667],
        ]

        ctrl_tpm = np.transpose(np.array(ctrl_tpm)).astype(np.float)
        ctrl_tpm = np.around(ctrl_tpm)
        test_tpm = np.around(sg.adata.layers["tpm"].toarray())
        print("tpm")
        print("test")
        print(test_tpm)
        print("control")
        print(ctrl_tpm)
        print(np.array_equal(ctrl_tpm, test_tpm))
        assert np.array_equal(ctrl_tpm, test_tpm)

        ctrl_pi = [
            [100, 100],
            [66.6666667, 0],
            [0, 66.6666667],
            [100, 100],
            [33.3333333, 33.3333333],
        ]
        ctrl_pi = np.transpose(np.array(ctrl_pi)).astype(np.float)
        ctrl_pi = np.around(ctrl_pi)
        test_pi = np.around(sg.adata.layers["pi"].toarray())
        print("pi")
        print("test")
        print(test_pi)
        print("ctrl")
        print(ctrl_pi)
        assert np.array_equal(ctrl_pi, test_pi)

    # add abundance - vanilla
    def test_add_abundance_1(self):
        sg = swan.SwanGraph()
        sg.add_annotation("files/test_full_annotation.gtf")
        sg.add_transcriptome("files/test_full.gtf")

        sg.add_abundance("files/test_ab_1.tsv")

        assert set(sg.datasets) == set(["dataset1", "dataset2"])

        print(sg.t_df.index)
        print(sg.adata.var.index)

        ctrl_tids = ["test1", "test2", "test3", "test4", "test5"]
        print(ctrl_tids)
        print(sg.adata.var.index)
        assert set(ctrl_tids) == set(sg.adata.var.index.tolist())

        ctrl_X = np.transpose(
            np.array([[5, 5], [10, 0], [0, 10], [10, 10], [5, 5]])
        ).astype(np.float)
        print("test")
        print(sg.adata.X)
        df = pd.DataFrame(
            index=sg.adata.obs.index,
            columns=sg.adata.var.index,
            data=sg.adata.layers["counts"].toarray(),
        )
        print(df.head())
        print("control")
        print(ctrl_X)
        assert np.array_equal(sg.adata.X.toarray(), ctrl_X)

        ctrl_tpm = [
            [166666.6667, 166666.6667],
            [333333.3333, 0],
            [0, 333333.3333],
            [333333.3333, 333333.3333],
            [166666.6667, 166666.6667],
        ]

        ctrl_tpm = np.transpose(np.array(ctrl_tpm)).astype(np.float)
        ctrl_tpm = np.around(ctrl_tpm)
        test_tpm = np.around(sg.adata.layers["tpm"].toarray())
        print("tpm")
        print("test")
        print(test_tpm)
        print("control")
        print(ctrl_tpm)
        print(np.array_equal(ctrl_tpm, test_tpm))
        assert np.array_equal(ctrl_tpm, test_tpm)

        ctrl_pi = [
            [100, 100],
            [66.6666667, 0],
            [0, 66.6666667],
            [100, 100],
            [33.3333333, 33.3333333],
        ]
        ctrl_pi = np.transpose(np.array(ctrl_pi)).astype(np.float)
        ctrl_pi = np.around(ctrl_pi)
        test_pi = np.around(sg.adata.layers["pi"].toarray())
        print("pi")
        print("test")
        print(test_pi)
        print("ctrl")
        print(ctrl_pi)
        assert np.array_equal(ctrl_pi, test_pi)


# test lower-level functions

# todo - calc_tpm, calc_pi
class TestLowLevel(object):

    # test calc_pi - use a different obs_col
    def test_calc_pi_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        sg.adata.obs["cluster"] = "c1"
        test_df, test_sums = swan.calc_pi(sg.adata, sg.t_df, obs_col="cluster")
        test_df = test_df.sparse.to_dense()
        test_df = test_df.transpose()
        test_df[["c1"]] = test_df[["c1"]].round()

        data = [
            ["test1", 100],
            ["test2", 33.333],
            ["test3", 33.333],
            ["test4", 100],
            ["test5", 33.33],
        ]
        cols = ["tid", "c1"]
        df = pd.DataFrame(data=data, columns=cols)
        df.set_index("tid", inplace=True)
        df.index.name = None
        df[["c1"]] = df[["c1"]].round()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()

    # test calc_pi - use dataset
    def test_calc_pi_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        test_df, test_sums = swan.calc_pi(sg.adata, sg.t_df, obs_col="dataset")
        test_df = test_df.sparse.to_dense()
        test_df = test_df.transpose()
        test_df[["dataset1", "dataset2"]] = test_df[["dataset1", "dataset2"]].round()

        data = [
            ["test1", 100, 100],
            ["test2", 66.6666667, 0],
            ["test3", 0, 66.6666667],
            ["test4", 100, 100],
            ["test5", 33.3333333, 33.3333333],
        ]
        cols = ["tid", "dataset1", "dataset2"]
        df = pd.DataFrame(data=data, columns=cols)
        df.set_index("tid", inplace=True)
        df.index.name = None
        df[["dataset1", "dataset2"]] = df[["dataset1", "dataset2"]].round()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()

        # don't test test_sums cause that's the direct output from
        # calc_total_counts

    # test calc_tpm - use different obs_col
    def test_calc_tpm_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        sg.adata.obs["cluster"] = ["c1", "c1"]
        test_df = swan.calc_tpm(sg.adata, obs_col="cluster")
        test_df = test_df.sparse.to_dense()
        test_df = test_df.transpose()
        test_df[["c1"]] = test_df[["c1"]].round()

        data = [
            ["test1", 166666.6667],
            ["test2", 166666.6667],
            ["test3", 166666.6667],
            ["test4", 333333.3333],
            ["test5", 166666.6667],
        ]
        cols = ["tid", "c1"]
        df = pd.DataFrame(data=data, columns=cols)
        df[["c1"]] = df[["c1"]].round()
        df.set_index("tid", inplace=True)
        df.index.name = None
        # df = df.transpose()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()

    # test calc_tpm - use dataset col
    def test_calc_tpm_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        test_df = swan.calc_tpm(sg.adata, obs_col="dataset")
        test_df = test_df.sparse.to_dense()
        test_df = test_df.transpose()
        test_df[["dataset1", "dataset2"]] = test_df[["dataset1", "dataset2"]].round()

        data = [
            ["test1", 166666.6667, 166666.6667],
            ["test2", 333333.3333, 0],
            ["test3", 0, 333333.3333],
            ["test4", 333333.3333, 333333.3333],
            ["test5", 166666.6667, 166666.6667],
        ]
        cols = ["tid", "dataset1", "dataset2"]
        df = pd.DataFrame(data=data, columns=cols)
        df[["dataset1", "dataset2"]] = df[["dataset1", "dataset2"]].round()
        df.set_index("tid", inplace=True)
        df.index.name = None
        # df = df.transpose()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()

    # test calc_total_counts - use additional metadata col
    def test_calc_total_counts_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        sg.adata.obs["cluster"] = ["c1", "c1"]

        test_df = swan.calc_total_counts(sg.adata, obs_col="cluster")

        data = [
            ["test1", 10],
            ["test2", 10],
            ["test3", 10],
            ["test4", 20],
            ["test5", 10],
        ]
        cols = ["tid", "c1"]
        df = pd.DataFrame(data=data, columns=cols)
        df.set_index("tid", inplace=True)
        df.index.name = None
        df = df.transpose()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()

    # test calc_total_counts - use dataset col
    def test_calc_total_counts_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome("files/test_full.gtf")
        sg.add_abundance("files/test_ab_1.tsv")
        test_df = swan.calc_total_counts(sg.adata)

        data = [
            ["test1", 5, 5],
            ["test2", 10, 0],
            ["test3", 0, 10],
            ["test4", 10, 10],
            ["test5", 5, 5],
        ]
        cols = ["tid", "dataset1", "dataset2"]
        df = pd.DataFrame(data=data, columns=cols)
        df.set_index("tid", inplace=True)
        df.index.name = None
        df = df.transpose()

        print("test")
        print(test_df)
        print("control")
        print(df)
        print(test_df == df)
        assert (test_df == df).all(axis=0).all()
