import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd
import anndata

###########################################################################
################# Related to input/error handling #########################
###########################################################################


###########################################################################
############################## Colors  ###################################
###########################################################################
class TestColors(object):

    # tests set_metadata_colors

    # test set_metadata_colors - vanilla
    def test_set_metadata_colors_1(self):
        sg = get_die_test_sg()
        cmap = {'GM12878': 'red', 'K562': 'blue'}
        test = sg.set_metadata_colors('sample', cmap)
        assert sg.adata.uns.sample_colors == ['red', 'blue']

    # test set_metadata_colors - obs_col does not exist
    def test_set_metadata_colors_1(self):
        sg = get_die_test_sg()
        cmap = {1: 'red', 2: 'blue'}
        with pytest.raises(Exception) as e:
            test = sg.set_metadata_colors('stage', cmap)
        assert 'Metadata column' in str(e.value)


###########################################################################
################# Related to plotting Swan Plots ##########################
###########################################################################
class TestPlotting(object):

    # done: test_new_gene, calc_pos_sizes, calc_edge_curves, plot_graph,
    # plot_transcript_path

    # init_plot_settings test do not check for indicate_novel / indicate settigns
    # init_plot_settings tests do not check for new dataset addition

    # test init_plot_settings - https://github.com/mortazavilab/swan_vis/issues/8
    # gene summary -> transcript path (same gene) -> gene summary (same gene)
    def test_init_9(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_graph('test2_gid', display=False)
        sg.plot_transcript_path('test5', display=False)
        sg.plot_graph('test2_gid', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting gene summary to
    # gene summary (same gene), also tests working from gene name
    def test_init_8(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_graph('test2_gid', display=False)
        sg.plot_graph('test2_gname', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting gene summary to
    # gene summary (different gene)
    def test_init_7(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_graph('test4_gid', display=False)
        sg.plot_graph('test2_gid', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting transcript path
    # to transcript path (same gene)
    def test_init_6(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_transcript_path('test3', display=False)
        sg.plot_transcript_path('test2', display=False)

        # edge_df
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', None],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        # sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'internal', None, None],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting transcript path
    # to transcript path (different gene)
    def test_init_5(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_transcript_path('test5', display=False)
        sg.plot_transcript_path('test2', display=False)

        # edge_df
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', None],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        # sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'internal', None, None],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting transcript path
    # to gene summary (same gene)
    def test_init_4(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_transcript_path('test2', display=False)
        sg.plot_graph('test2_gid', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting transcript path
    # to gene summary (different gene)
    def test_init_3(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_transcript_path('test1', display=False)
        sg.plot_graph('test2_gid', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting gene summary to
    # transcript path (same gene)
    def test_init_1(self):

        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.plot_graph('test2_gid', display=False)

        sg.plot_transcript_path('test2', display=False)

        # edge_df
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', None],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        # sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'internal', None, None],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test init_plot_settings - going from plotting gene summary to
    # transcript path (different gene)
    def test_init_2(self):

        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.plot_graph('test4_gid', display=False)

        sg.plot_transcript_path('test2', display=False)

        # edge_df
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', None],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        # sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'internal', None, None],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test  - indicate_dataset
    def test_summary_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.datasets = ['dataset_1', 'dataset_2']
        sg.t_df['dataset_1'] = [False, False, False, False, True]
        sg.t_df['dataset_2'] = [True, True, True, True, True]

        sg.edge_df['dataset_1'] = False
        sg.edge_df['dataset_2'] = True
        edges = [5,11,12]
        sg.edge_df.loc[edges, 'dataset_1'] = True

        sg.loc_df['dataset_1'] = False
        sg.loc_df['dataset_2'] = True
        locs = [12,11,8,7]
        sg.loc_df.loc[locs, 'dataset_1'] = True

        sg.plot_graph('test2_gid', display=False, indicate_dataset='dataset_1')

        # edge_df
        sg.pg.edge_df.drop(['dataset_1', 'dataset_2', 'curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', 'dashed'],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', 'dashed'],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', 'dashed']]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        sg.pg.loc_df.drop(['dataset_1', 'dataset_2'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', 'node_outline', 2],
                [1, 'chr2', 80, True, False, False, 'internal', 'node_outline', 2],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', 'node_outline', 2],
                [5, 'chr2', 50, True, False, True, 'TES', 'node_outline', 2],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test calc_node_edge_styles - indicate_novel
    def def_summary_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.annotation = True
        sg.t_df['annotation'] = [True, True, True, True, False]

        sg.edge_df['annotation'] = True
        edges = [5,11,12]
        sg.edge_df.loc[edges, 'annotation'] = False

        sg.loc_df['annotation'] = True
        locs = [12,11,8,7]
        sg.loc_df.loc[locs, 'annotation'] = False

        sg.plot_graph('test2_gid', display=False, indicate_novel=True)

        # edge_df
        sg.pg.edge_df.drop(['annotation', 'curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', 'dashed'],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', 'dashed'],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', 'dashed']]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', 'node_outline', 2],
                [1, 'chr2', 80, True, False, False, 'internal', 'node_outline', 2],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', 'node_outline', 2],
                [5, 'chr2', 50, True, False, True, 'TES', 'node_outline', 2],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test calc_node_edge_styles - vanilla transcript path
    #     transcript_path where role in transcript
    #     is different ie transcripts [1,2,3], [0,1,2,3] 1 needs to be colored int
    #     sth with this guy? [5, 'chr2', 50, True, False, True, 'TES'],
    def test_path_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_transcript_path('test2', display=False)

        # edge_df
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', None],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        # sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'internal', None, None],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test calc_node_edge_styles - indicate_dataset transcript_path
    def test_path_2(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.datasets = ['dataset_1', 'dataset_2']
        sg.t_df['dataset_1'] = [False, False, False, False, True]
        sg.t_df['dataset_2'] = [True, True, True, True, True]

        sg.edge_df['dataset_1'] = False
        sg.edge_df['dataset_2'] = True
        edges = [5,11,12]
        sg.edge_df.loc[edges, 'dataset_1'] = True

        sg.loc_df['dataset_1'] = False
        sg.loc_df['dataset_2'] = True
        locs = [12,11,8,7]
        sg.loc_df.loc[locs, 'dataset_1'] = True

        sg.plot_transcript_path('test2', display=False, indicate_dataset='dataset_1')

        # edge_df
        sg.pg.edge_df.drop(['dataset_1', 'dataset_2', 'curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', 'dashed'],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', 'dashed'],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', 'dashed']]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        sg.pg.loc_df.drop(['dataset_1', 'dataset_2'], axis=1, inplace=True) # not checking this
        data = [[0, 'chr2', 100, False, True, False, 'TSS', 'node_outline', 2],
                [1, 'chr2', 80, True, False, False, 'internal', 'node_outline', 2],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', 'node_outline', 2],
                [5, 'chr2', 50, True, False, True, 'internal', 'node_outline', 2],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

    # test calc_node_edge_styles - indicate_novel transcript_path
    def test_path_3(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')

        sg.annotation = True
        sg.t_df['annotation'] = [True, True, True, True, False]

        sg.edge_df['annotation'] = True
        edges = [5,11,12]
        sg.edge_df.loc[edges, 'annotation'] = False

        sg.loc_df['annotation'] = True
        locs = [12,11,8,7]
        sg.loc_df.loc[locs, 'annotation'] = False

        sg.plot_transcript_path('test2', display=False, indicate_novel=True)

        # edge_df
        sg.pg.edge_df.drop(['annotation', 'curve'], axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon_gray', 'dashed'],
                [14, '-', 'intron', 3, 5, 'intron_gray', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon_gray', None],
                [11, '-', 'intron', 1, 4, 'intron_gray', 'dashed'],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', 'dashed']]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        sg.pg.loc_df.drop(['annotation'], axis=1, inplace=True)
        data = [[0, 'chr2', 100, False, True, False, 'TSS', 'node_outline', 2],
                [1, 'chr2', 80, True, False, False, 'internal', 'node_outline', 2],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal_gray', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', 'node_outline', 2],
                [5, 'chr2', 50, True, False, True, 'internal', 'node_outline', 2],
                [6, 'chr2', 45, False, False, True, 'TES', None, None]]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)


    # test calc_node_edge_styles - vanilla
    def test_summary_1(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_full.gtf')
        sg.plot_graph('test2_gid', display=False)

        # edge_df
        sg.pg.edge_df.drop('curve', axis=1, inplace=True) # not checking this
        data = [[9, '-', 'exon', 5, 6, 'exon', None],
                [8, '-', 'intron', 4, 5, 'intron', None],
                [12, '-', 'exon', 4, 5, 'exon', None],
                [14, '-', 'intron', 3, 5, 'intron', None],
                [7, '-', 'exon', 2, 4, 'exon', None],
                [13, '-', 'exon', 2, 3, 'exon', None],
                [11, '-', 'intron', 1, 4, 'intron', None],
                [6, '-', 'intron', 1, 2, 'intron', None],
                [5, '-', 'exon', 0, 1, 'exon', None]]
        cols = ['edge_id', 'strand', 'edge_type', 'v1', 'v2', 'color', 'line']
        ctrl_edge_df = pd.DataFrame(data=data, columns=cols)
        ctrl_edge_df = swan.create_dupe_index(ctrl_edge_df, 'edge_id')
        ctrl_edge_df = swan.set_dupe_index(ctrl_edge_df, 'edge_id')

        # loc_df
        data = [[0, 'chr2', 100, False, True, False, 'TSS', None, None],
                [1, 'chr2', 80, True, False, False, 'internal', None, None],
                [2, 'chr2', 75, True, False, False, 'internal', None, None],
                [3, 'chr2', 65, True, False, False, 'internal', None, None],
                [4, 'chr2', 60, True, False, False, 'internal', None, None],
                [5, 'chr2', 50, True, False, True, 'TES'],
                [6, 'chr2', 45, False, False, True, 'TES']]
        cols = ['vertex_id', 'chrom', 'coord', 'internal', 'TSS', 'TES', \
                'color', 'edgecolor', 'linewidth']
        ctrl_loc_df = pd.DataFrame(data=data, columns=cols)
        ctrl_loc_df = swan.create_dupe_index(ctrl_loc_df, 'vertex_id')
        ctrl_loc_df = swan.set_dupe_index(ctrl_loc_df, 'vertex_id')

        check_dfs(sg.pg.loc_df, ctrl_loc_df, sg.pg.edge_df, ctrl_edge_df)

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
        sg.pg.edge_df.drop(['curve'], axis=1, inplace=True)
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
    # assert (loc_df == ctrl_loc_df).all(axis=0).all()
    assert loc_df.equals(ctrl_loc_df)

    print('test')
    print(edge_df)
    print('control')
    print(ctrl_edge_df)
    assert edge_df.equals(ctrl_edge_df)

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
