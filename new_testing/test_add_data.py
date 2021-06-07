import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd

###########################################################################
###################### Related to file parsing ############################
###########################################################################
class TestFiles(object):

    # tests GTF parsing
    def test_parse_gtf(self):
        gtf_file = 'files/Canx.gtf'
        t_df, exon_df, from_talon = swan.parse_gtf(gtf_file, False)

        t_df.index.name = 'tid_index'
        t_df = t_df.sort_values(by='tid_index')

        ctrl_t_df = pd.read_csv('files/Canx_transcript.tsv',sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df = ctrl_t_df.sort_values(by='tid_index')

        ctrl_exons = ctrl_t_df.exons.tolist()
        ctrl_exons = [exons.split(',') for exons in ctrl_exons]
        ctrl_t_df['exons'] = ctrl_exons

        print(t_df == ctrl_t_df)
        assert (t_df == ctrl_t_df).all(axis=0).all()

    # # tests TALON DB parsing
    # def test_parse_db(self):
    #     pass
    #     # TODO


###########################################################################
####################### Related to DF creation ############################
###########################################################################
class TestCreateDFs(object):

    # add_edge_coords, get_current_locs, get_current_edges,
    # create_loc_dict, create_transcript_edge_dict create_dfs,

    # tests add_edge_coords
    def test_add_edge_coords(self):
        sg = swan.SwanGraph()
        sg = add_transcriptome_no_reorder_gtf(sg, 'files/test_full.gtf')
        # sg.add_transcriptome('files/test_full.gtf')
        cols = ['edge_id', 'v1', 'v2', 'strand', 'edge_type',
                'v1_coord', 'v2_coord']

        # print(sg.edge_df.head())
        edge_df = sg.add_edge_coords()
        print(edge_df.head())
        edge_df = edge_df[cols]

        ctrl_edge_df = pd.read_csv('files/test_add_edge_coords_result.tsv', sep='\t')
        ctrl_edge_df = ctrl_edge_df[cols]

        # first order to make them comparable
        # sort all values by their IDs
        edge_df.sort_values(by='edge_id', inplace=True)
        ctrl_edge_df.sort_values(by='edge_id', inplace=True)

        # and order columns the same way
        ctrl_edge_df = ctrl_edge_df[edge_df.columns]

        print('test')
        print(edge_df)
        print('control')
        print(ctrl_edge_df)
        assert (edge_df == ctrl_edge_df).all(axis=0).all()

    # tests get_current_locs with an empty swangraph
    def test_get_current_locs_empty_sg(self):
        sg = swan.SwanGraph()
        locs, n = sg.get_current_locs()

        assert locs == {}
        assert n == -1

    # tests get_current_locs with a swangraph with data
    def test_get_current_locs_sg_data(self):
        sg = swan.SwanGraph()
        cols = ['vertex_id', 'chrom', 'coord']
        data = [[0, 1, 2], [1, 1, 3], [2, 3, 50]]
        sg.loc_df = pd.DataFrame(data=data, columns=cols)
        cols = ['tid']
        data = [0]
        sg.t_df = pd.DataFrame(data=data, columns=cols)
        locs, n = sg.get_current_locs()

        ctrl_locs = {(1,2):0, (1,3):1, (3,50):2}
        assert locs == ctrl_locs
        assert n == 2

    # tests get_current_edges with an empty swangraph
    def test_get_current_edges_empty_sg(self):
        sg = swan.SwanGraph()
        edges, n = sg.get_current_edges()

        assert(edges == {})
        assert(n == -1)

    # tests get_current_edges in a sg with data
    def test_get_current_edges_sg_data(self):
        sg = swan.SwanGraph()

        cols = ['vertex_id', 'chrom', 'coord']
        data = [[0, 1, 2], [1, 1, 3], [2, 1, 50]]
        sg.loc_df = pd.DataFrame(data=data, columns=cols)

        cols = ['edge_id', 'v1', 'v2', 'strand', 'edge_type']
        data = [[0, 0, 1, '+', 'exon'],
                [1, 1, 2, '+', 'intron']]
        sg.edge_df = pd.DataFrame(data=data, columns=cols)

        cols = ['tid']
        data = [0]
        sg.t_df = pd.DataFrame(data=data, columns=cols)

        edges, n = sg.get_current_edges()
        ctrl = {(1,2,3,'+','exon'): {'edge_id': 0,
                                     'edge_type': 'exon',
                                     'v1': 0 ,
                                     'v2': 1},
                (1,3,50,'+','intron'): {'edge_id': 1,
                                        'edge_type': 'intron',
                                        'v1': 1,
                                        'v2': 2}}
        assert(edges == ctrl)
        assert(n == 1)

    # test create_loc_dict on an empty sg
    # also checks to make sure exons that use the same loc
    # don't result in dupe entries in loc_df
    def test_create_loc_dict_empty_sg(self):
        _, exons = get_test_transcript_exon_dicts()

        sg = swan.SwanGraph()
        locs = sg.create_loc_dict(exons)
        ctrl_locs = {('chr1',1): 0,
            ('chr1', 20): 1,
            ('chr1', 25): 2,
            ('chr1', 30): 3,
            ('chr1', 35): 4,
            ('chr1', 40): 5,
            ('chr2', 100): 6,
            ('chr2', 80): 7,
            ('chr2', 75): 8,
            ('chr2', 60): 9,
            ('chr2', 50): 10,
            ('chr2', 45): 11,
            ('chr2', 65): 12
             }
        assert(ctrl_locs == locs)

    # tests create_loc_dict when locs already exist in sg
    def test_create_loc_dict_sg_data(self):
        _, exons = get_test_transcript_exon_dicts()

        # dummy preexisting data
        sg = swan.SwanGraph()
        data = [[0, 'chr1', 1], [1, 'chr2', 80]]
        columns = ['vertex_id', 'chrom', 'coord']
        sg.loc_df = pd.DataFrame(data=data, columns=columns)

        cols = ['tid']
        data = [0]
        sg.t_df = pd.DataFrame(data=data, columns=cols)

        locs = sg.create_loc_dict(exons)

        ctrl_locs = {('chr1', 1):0,
            ('chr2', 80): 1,
            ('chr1', 20): 2,
            ('chr1', 25): 3,
            ('chr1', 30): 4,
            ('chr1', 35): 5,
            ('chr1', 40): 6,
            ('chr2', 100): 7,
            ('chr2', 75): 8,
            ('chr2', 60): 9,
            ('chr2', 50): 10,
            ('chr2', 45): 11,
            ('chr2', 65): 12
             }

        print('test')
        print(locs)
        print('control')
        print(ctrl_locs)
        assert(ctrl_locs == locs)

    # tests create_transcript_edge_dict empty swangraph
    def test_create_transcript_edge_dict_emtpy_sg(self):
        transcripts, exons = get_test_transcript_exon_dicts()

        sg = swan.SwanGraph()
        locs = sg.create_loc_dict(exons)

        transcripts, edges = sg.create_transcript_edge_dicts(transcripts, exons, locs)

        # just compare the paths for the transcripts, which is the only
        # part modified by this function
        transcripts = dict([(key, item['path']) for key, item in transcripts.items()])
        ctrl_transcript_paths = {
            'test1': [0,1,2,3,4],
            'test2': [5,6,7,8,9],
            'test3': [5,6,10,11,9],
            'test4': [12],
            'test5': [5,13,14]
        }

        assert(transcripts == ctrl_transcript_paths)

        ctrl_edges = {
            ('chr1', 1, 20, '+', 'exon'): {
                'edge_id': 0,
                'edge_type': 'exon',
                'v1': 0,
                'v2': 1
            },
            ('chr1', 20, 25, '+', 'intron'): {
                'edge_id': 1,
                'edge_type': 'intron',
                'v1': 1,
                'v2': 2
            },
            ('chr1', 25, 30, '+', 'exon'): {
                'edge_id': 2,
                'edge_type': 'exon',
                'v1': 2,
                'v2': 3
            },
            ('chr1', 30, 35, '+', 'intron'): {
                'edge_id': 3,
                'edge_type': 'intron',
                'v1': 3,
                'v2': 4
            },
            ('chr1', 35, 40, '+', 'exon'): {
                'edge_id': 4,
                'edge_type': 'exon',
                'v1': 4,
                'v2': 5
            },
            ('chr2', 100, 80, '-', 'exon'): {
                'edge_id': 5,
                'edge_type': 'exon',
                'v1': 6,
                'v2': 7
            },
            ('chr2', 80, 75, '-', 'intron'): {
                'edge_id': 6,
                'edge_type': 'intron',
                'v1': 7,
                'v2': 8
            },
            ('chr2', 75, 60, '-', 'exon'): {
                'edge_id': 7,
                'edge_type': 'exon' ,
                'v1': 8,
                'v2': 9
            },
            ('chr2', 60, 50, '-', 'intron'): {
                'edge_id': 8,
                'edge_type': 'intron',
                'v1': 9,
                'v2': 10
            },
            ('chr2', 50, 45, '-', 'exon'): {
                'edge_id': 9,
                'edge_type': 'exon',
                'v1': 10,
                'v2': 11
            },
            ('chr2', 75, 65, '-', 'exon'): {
                'edge_id': 10,
                'edge_type': 'exon',
                'v1': 8,
                'v2': 12
            },
            ('chr2', 65, 50, '-', 'intron'): {
                'edge_id': 11,
                'edge_type': 'intron',
                'v1': 12,
                'v2': 10
            },
            ('chr2', 45, 50, '+', 'exon'): {
                'edge_id': 12,
                'edge_type': 'exon',
                'v1': 11,
                'v2': 10
            },
            ('chr2', 80, 60, '-', 'intron'): {
                'edge_id': 13,
                'edge_type': 'intron',
                'v1': 7,
                'v2': 9
            },
            ('chr2', 60, 50, '-', 'exon'): {
                'edge_id': 14,
                'edge_type': 'exon',
                'v1': 9,
                'v2': 10
            }
        }

        assert(edges == ctrl_edges)

    # tests create_transcript_edge_dict with edges already in swangraph
    def test_create_transcript_edge_dict_edge_sg(self):
        transcripts, exons = get_test_transcript_exon_dicts()

        # add some dummy data
        sg = swan.SwanGraph()

        data = [[0, 'chr1', 1],
                [1, 'chr2', 20],
                [2, 'chr2', 100],
                [3, 'chr2', 80]]
        columns = ['vertex_id', 'chrom', 'coord']
        sg.loc_df = pd.DataFrame(data=data, columns=columns)
        cols = ['tid']
        data = [0]
        sg.t_df = pd.DataFrame(data=data, columns=cols)
        locs = sg.create_loc_dict(exons)

        data = [[0, 0, 1, '+', 'exon'],
                [1, 2, 3, '-', 'exon']]
        columns = ['edge_id', 'v1', 'v2', 'strand', 'edge_type']
        sg.edge_df = pd.DataFrame(data=data, columns=columns)

        transcripts, edges = sg.create_transcript_edge_dicts(transcripts, exons, locs)

        # just compare the paths for the transcripts, which is the only
        # part modified by this function
        transcripts = dict([(key, item['path']) for key, item in transcripts.items()])
        ctrl_transcript_paths = {
            'test1': [0,2,3,4,5],
            'test2': [1,6,7,8,9],
            'test3': [1,6,10,11,9],
            'test4': [12],
            'test5': [1,13,14]
        }

        assert(transcripts == ctrl_transcript_paths)

        ctrl_edges = {
            ('chr1', 1, 20, '+', 'exon'): {
                'edge_id': 0,
                'edge_type': 'exon',
                'v1': 0,
                'v2': 1
            },
            ('chr1', 20, 25, '+', 'intron'): {
                'edge_id': 2,
                'edge_type': 'intron',
                'v1': 4,
                'v2': 5
            },
            ('chr1', 25, 30, '+', 'exon'): {
                'edge_id': 3,
                'edge_type': 'exon',
                'v1': 5,
                'v2': 6
            },
            ('chr1', 30, 35, '+', 'intron'): {
                'edge_id': 4,
                'edge_type': 'intron',
                'v1': 6,
                'v2': 7
            },
            ('chr1', 35, 40, '+', 'exon'): {
                'edge_id': 5,
                'edge_type': 'exon',
                'v1': 7,
                'v2': 8
            },
            ('chr2', 100, 80, '-', 'exon'): {
                'edge_id': 1,
                'edge_type': 'exon',
                'v1': 2,
                'v2': 3
            },
            ('chr2', 80, 75, '-', 'intron'): {
                'edge_id': 6,
                'edge_type': 'intron',
                'v1': 3,
                'v2': 9
            },
            ('chr2', 75, 60, '-', 'exon'): {
                'edge_id': 7,
                'edge_type': 'exon' ,
                'v1': 9,
                'v2': 10
            },
            ('chr2', 60, 50, '-', 'intron'): {
                'edge_id': 8,
                'edge_type': 'intron',
                'v1': 10,
                'v2': 11
            },
            ('chr2', 50, 45, '-', 'exon'): {
                'edge_id': 9,
                'edge_type': 'exon',
                'v1': 11,
                'v2': 12
            },
            ('chr2', 75, 65, '-', 'exon'): {
                'edge_id': 10,
                'edge_type': 'exon',
                'v1': 9,
                'v2': 13
            },
            ('chr2', 65, 50, '-', 'intron'): {
                'edge_id': 11,
                'edge_type': 'intron',
                'v1': 13,
                'v2': 11
            },
            ('chr2', 45, 50, '+', 'exon'): {
                'edge_id': 12,
                'edge_type': 'exon',
                'v1': 12,
                'v2': 11
            },
            ('chr2', 80, 60, '-', 'intron'): {
                'edge_id': 13,
                'edge_type': 'intron',
                'v1': 3,
                'v2': 10
            },
            ('chr2', 60, 50, '-', 'exon'): {
                'edge_id': 14,
                'edge_type': 'exon',
                'v1': 10,
                'v2': 11
            }
        }

        assert(edges == ctrl_edges)

    # # tests create_transcript_edge_dict where transcripts already
#     # # exist in the swangraph
#     # def test_create_transcript_edge_dict_edge_t_sg(self):
#     #     pass
#     #     # TODO
#
    # tests create_dfs with an empty sg
    # also ensures that empty dict -> df -> dict conversion doesn't screw up
    def test_create_dfs_empty_sg(self):
        transcripts, exons = get_test_transcript_exon_dicts()
        sg = swan.SwanGraph()

        loc_df, edge_df, t_df = sg.create_dfs(transcripts, exons, False)

        ctrl_loc_df = pd.read_csv('files/test_loc_df.tsv', sep='\t')
        ctrl_loc_df.set_index('vertex_id_index', inplace=True)
        ctrl_loc_df.index.name = 'vertex_id'

        # remove the columns that are there just for debugging purposes
        ctrl_edge_df = pd.read_csv('files/test_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        ctrl_edge_df.set_index('edge_id_index', inplace=True)
        ctrl_edge_df.index.name = 'edge_id'

        # again, remove and reformat columns that are there for debugging
        ctrl_t_df = pd.read_csv('files/test_t_df.tsv', sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df.index.name = 'tid'
        ctrl_t_df.drop(['loc_path', 'novelty'], axis=1, inplace=True)
        ctrl_t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        ctrl_t_df['path'] = ctrl_t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)

        check_dfs(loc_df, ctrl_loc_df, edge_df, ctrl_edge_df, t_df, ctrl_t_df)


    # tests create_dfs when from_talon = True
    def test_create_dfs_empty_sg_from_talon(self):
        transcripts, exons = get_test_transcript_exon_dicts()
        sg = swan.SwanGraph()

        loc_df, edge_df, t_df = sg.create_dfs(transcripts, exons, True)

        ctrl_loc_df = pd.read_csv('files/test_loc_df.tsv', sep='\t')
        ctrl_loc_df.set_index('vertex_id_index', inplace=True)
        ctrl_loc_df.index.name = 'vertex_id'

        # remove the columns that are there just for debugging purposes
        ctrl_edge_df = pd.read_csv('files/test_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        ctrl_edge_df.set_index('edge_id_index', inplace=True)
        ctrl_edge_df.index.name = 'edge_id'

        # again, remove and reformat columns that are there for debugging
        ctrl_t_df = pd.read_csv('files/test_t_df.tsv', sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df.index.name = 'tid'
        ctrl_t_df.drop(['loc_path'], axis=1, inplace=True)
        ctrl_t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        ctrl_t_df['path'] = ctrl_t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        ctrl_t_df = ctrl_t_df[t_df.columns]

        check_dfs(loc_df, ctrl_loc_df, edge_df, ctrl_edge_df, t_df, ctrl_t_df)


    # tests create_dfs in a swangraph with data
    def test_create_dfs_data_sg(self):
        transcripts, exons = get_test_transcript_exon_dicts()

        del transcripts['test2']
        sg = swan.SwanGraph()

        # add dummy data

        # loc_df - format
        loc_df = pd.read_csv('files/test_preexisting_loc_df.tsv', sep='\t')
        loc_df.set_index('vertex_id_index', inplace=True)
        loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        edge_df = pd.read_csv('files/test_preexisting_edge_df.tsv', sep='\t')
        edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        edge_df.set_index('edge_id_index', inplace=True)
        edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        t_df = pd.read_csv('files/test_preexisting_t_df.tsv', sep='\t')
        t_df.set_index('tid_index', inplace=True)
        t_df.index.name = 'tid'
        t_df.drop(['loc_path'], axis=1, inplace=True)
        t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        t_df['path'] = t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        t_df = t_df[t_df.columns]

        sg.loc_df = loc_df
        sg.edge_df = edge_df
        sg.t_df = t_df

        loc_df, edge_df, t_df = sg.create_dfs(transcripts, exons, True)


        # control data

        # loc_df - format
        ctrl_loc_df = pd.read_csv('files/test_preexisting_result_loc_df.tsv', sep='\t')
        ctrl_loc_df.set_index('vertex_id_index', inplace=True)
        ctrl_loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        ctrl_edge_df = pd.read_csv('files/test_preexisting_result_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        ctrl_edge_df.set_index('edge_id_index', inplace=True)
        ctrl_edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        ctrl_t_df = pd.read_csv('files/test_preexisting_result_t_df.tsv', sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df.index.name = 'tid'
        ctrl_t_df.drop(['loc_path'], axis=1, inplace=True)
        ctrl_t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        ctrl_t_df['path'] = ctrl_t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        ctrl_t_df = ctrl_t_df[t_df.columns]

        check_dfs(loc_df, ctrl_loc_df, edge_df, ctrl_edge_df, t_df, ctrl_t_df)

    # tests create_dfs in sg with data where existing data has novelty
    # and added dataset does not
    def test_create_dfs_data_sg_nov1(self):
        transcripts, exons = get_test_transcript_exon_dicts()

        # to do - remove transcript that's already there
        sg = swan.SwanGraph()

        # add dummy data

        # loc_df - format
        loc_df = pd.read_csv('files/test_preexisting_loc_df.tsv', sep='\t')
        loc_df.set_index('vertex_id_index', inplace=True)
        loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        edge_df = pd.read_csv('files/test_preexisting_edge_df.tsv', sep='\t')
        edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        edge_df.set_index('edge_id_index', inplace=True)
        edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        t_df = pd.read_csv('files/test_preexisting_t_df.tsv', sep='\t')
        t_df.set_index('tid_index', inplace=True)
        t_df.index.name = 'tid'
        t_df.drop(['loc_path'], axis=1, inplace=True)
        t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        t_df['path'] = t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        t_df = t_df[t_df.columns]

        sg.loc_df = loc_df
        sg.edge_df = edge_df
        sg.t_df = t_df

        loc_df, edge_df, t_df = sg.create_dfs(transcripts, exons, False)


        # control data

        # loc_df - format
        ctrl_loc_df = pd.read_csv('files/test_preexisting_result_loc_df.tsv', sep='\t')
        ctrl_loc_df.set_index('vertex_id_index', inplace=True)
        ctrl_loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        ctrl_edge_df = pd.read_csv('files/test_preexisting_result_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        ctrl_edge_df.set_index('edge_id_index', inplace=True)
        ctrl_edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        ctrl_t_df = pd.read_csv('files/test_preexisting_result_t_df.tsv', sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df.index.name = 'tid'
        ctrl_t_df.drop(['loc_path'], axis=1, inplace=True)
        ctrl_t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        ctrl_t_df['path'] = ctrl_t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        ctrl_t_df = ctrl_t_df[t_df.columns]

        # remove novelty for entries that are new
        new_tids = ['test1', 'test3', 'test4', 'test5']
        ctrl_t_df.loc[ctrl_t_df.tid.isin(new_tids), 'novelty'] = 'Undefined'

        check_dfs(loc_df, ctrl_loc_df, edge_df, ctrl_edge_df, t_df, ctrl_t_df)


    # tests create_dfs with preexisting data and a duplicate transcript
    # being added
    # also tests that old data (novelty in this case) is not overwritten
    def test_create_dfs_data_sg_dupe(self):
        transcripts, exons = get_test_transcript_exon_dicts()
        sg = swan.SwanGraph()

        # add dummy data

        # loc_df - format
        loc_df = pd.read_csv('files/test_preexisting_loc_df.tsv', sep='\t')
        loc_df.set_index('vertex_id_index', inplace=True)
        loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        edge_df = pd.read_csv('files/test_preexisting_edge_df.tsv', sep='\t')
        edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        edge_df.set_index('edge_id_index', inplace=True)
        edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        t_df = pd.read_csv('files/test_preexisting_t_df.tsv', sep='\t')
        t_df.set_index('tid_index', inplace=True)
        t_df.index.name = 'tid'
        t_df.drop(['loc_path'], axis=1, inplace=True)
        t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        t_df['path'] = t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        t_df = t_df[t_df.columns]

        sg.loc_df = loc_df
        sg.edge_df = edge_df
        sg.t_df = t_df

        loc_df, edge_df, t_df = sg.create_dfs(transcripts, exons, False)

        # control data

        # loc_df - format
        ctrl_loc_df = pd.read_csv('files/test_preexisting_result_loc_df.tsv', sep='\t')
        ctrl_loc_df.set_index('vertex_id_index', inplace=True)
        ctrl_loc_df.index.name = 'vertex_id'

        # edge_df - format and remove the columns that are there
        # just for debugging purposes
        ctrl_edge_df = pd.read_csv('files/test_preexisting_result_edge_df.tsv', sep='\t')
        ctrl_edge_df.drop(['v2_coord', 'v1_coord'], axis=1, inplace=True)
        ctrl_edge_df.set_index('edge_id_index', inplace=True)
        ctrl_edge_df.index.name = 'edge_id'

        # t_df - remove and reformat columns that are there for debugging
        ctrl_t_df = pd.read_csv('files/test_preexisting_result_t_df.tsv', sep='\t')
        ctrl_t_df.set_index('tid_index', inplace=True)
        ctrl_t_df.index.name = 'tid'
        ctrl_t_df.drop(['loc_path'], axis=1, inplace=True)
        ctrl_t_df.rename({'edge_path': 'path'}, axis=1, inplace=True)
        ctrl_t_df['path'] = ctrl_t_df.apply(lambda x: [int(n) for n in x.path.split(',')], axis=1)
        ctrl_t_df = ctrl_t_df[t_df.columns]

        # remove novelty for entries that are new
        new_tids = ['test1', 'test3', 'test4', 'test5']
        ctrl_t_df.loc[ctrl_t_df.tid.isin(new_tids), 'novelty'] = 'Undefined'

        check_dfs(loc_df, ctrl_loc_df, edge_df, ctrl_edge_df, t_df, ctrl_t_df)

# ###########################################################################
# ############### Related higher-level dataset addition #####################
# ###########################################################################
class TestDataset(object):

    # tests: get_loc_types, label_annotated

    # TODO
    #  add_dataset, add_transcriptome, add_annotation tests with the
    # update_ids results also tested!!!! ie coords and paths!

    # test get_loc_types
    def test_get_loc_types(self):
        sg = swan.SwanGraph()
        data = [[0, [0,1], [0,1,2]],
                [1, [2,3], [3,4,5]],
                 [2, [4,5], [6,7,8]]]
        sg.t_df = pd.DataFrame(data=data, columns=['tid', 'path', 'loc_path'])

        data = [[0,0,1], [1,1,2], [2,3,4], [3,4,5],
                [4,6,7], [5,7,8]]
        sg.edge_df = pd.DataFrame(data=data, columns=['edge_id', 'v1', 'v2'])

        data = [0,1,2,3,4,5,6,7,8]
        sg.loc_df = pd.DataFrame(data=data, columns=['vertex_id'])

        sg.get_loc_types()

        ctrl_tss = [0,3,6]
        tss = sg.loc_df.loc[sg.loc_df.TSS == True, 'vertex_id'].tolist()
        assert set(ctrl_tss) == set(tss)
        ctrl_int = [1,4,7]
        int = sg.loc_df.loc[sg.loc_df.internal == True, 'vertex_id'].tolist()
        assert set(ctrl_int) == set(int)
        ctrl_tes = [2,5,8]
        tes = sg.loc_df.loc[sg.loc_df.TES == True, 'vertex_id'].tolist()
        assert set(ctrl_tes) == set(ctrl_tes)


    # label annotated transcripts
    def test_label_annotated(self):
        sg = swan.SwanGraph()
        data = [[0, [0,1]],
                [1, [2,3]],
                 [2, [4,5]]]
        sg.t_df = pd.DataFrame(data=data, columns=['tid', 'path'])

        data = [[0,0,1], [1,1,2], [2,2,3], [3,3,4],
                [4,4,5], [5,5,6]]
        sg.edge_df = pd.DataFrame(data=data, columns=['edge_id', 'v1', 'v2'])

        data = [0,1,2,3,4,5,6]
        sg.loc_df = pd.DataFrame(data=data, columns=['vertex_id'])

        tids = [0,1]
        sg.label_annotated(tids)

        ctrl_tids = [0,1]
        tids = sg.t_df.loc[sg.t_df.annotation == True, 'tid'].tolist()
        assert set(ctrl_tids) == set(tids)

        ctrl_edges = [0,1,2,3]
        edges = sg.edge_df.loc[sg.edge_df.annotation == True, 'edge_id'].tolist()
        assert set(ctrl_edges) == set(edges)

        ctrl_locs = [0,1,2,3,4]
        locs = sg.loc_df.loc[sg.loc_df.annotation == True, 'vertex_id'].tolist()
        assert set(ctrl_locs) == set(locs)

    # add to empty sg, don't add isms
    def test_add_transcriptome(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_novel_talon.gtf', include_isms=False)

        print(sg.t_df)
        assert "ISM" not in sg.t_df.novelty.unique()
        # assert 1 == 0

    # tests if correct error is thrown when adding annotation to
    # sg that already has one
    def test_add_annotation_already(self):
        sg = swan.SwanGraph()
        sg.annotation = True
        with pytest.raises(Exception) as e:
            sg.add_annotation('files/Canx.gtf')
        assert 'Annotation already' in str(e.value)

    # add annotation to empty sg
    def test_add_annotation_empty_sg(self):
        sg = swan.SwanGraph()
        sg.add_annotation('files/test_full.gtf')

        # check annotation columns
        assert all(sg.t_df.annotation.tolist())
        assert all(sg.edge_df.annotation.tolist())
        assert all(sg.loc_df.annotation.tolist())

        # check novelty column in t_df
        assert len(sg.t_df.loc[sg.t_df.novelty=='Known']) == len(sg.t_df.index)

        # check annotation flag
        assert sg.annotation == True

    # add annotation to sg with data already in it
    def test_add_annotation_sg_data(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_novel.gtf')
        sg.add_annotation('files/test_known.gtf')

        # check with coord/chr bc of reindexing fuckery not being
        # remimplemented yet

        # t_df
        annot_tids = ['test1', 'test2', 'test4']
        assert all(sg.t_df.loc[annot_tids, 'annotation'])
        ctrl_novel_tids = ['test3', 'test5']
        novel_tids = sg.t_df.loc[sg.t_df.annotation == False, 'tid'].tolist()
        assert len(set(ctrl_novel_tids)-set(novel_tids)) == 0
        assert len(ctrl_novel_tids) == len(novel_tids)

        # make sure the novelty assignment worked
        annot_tids = sg.t_df.loc[sg.t_df.annotation == True, 'tid'].tolist()
        known_tids = sg.t_df.loc[sg.t_df.novelty == 'Known', 'tid'].tolist()
        assert set(annot_tids) == set(known_tids)

        annot_tids = sg.t_df.loc[sg.t_df.annotation == False, 'tid'].tolist()
        known_tids = sg.t_df.loc[sg.t_df.novelty == 'Undefined', 'tid'].tolist()
        assert set(annot_tids) == set(known_tids)

        # loc_df
        ctrl_novel_locs = [('chr2', 65)]
        temp = sg.loc_df[sg.loc_df.annotation == False]
        chroms = temp.chrom.tolist()
        coords = temp.coord.tolist()
        novel_locs = [(chrom, coord) for chrom, coord in zip(chroms, coords)]
        print('control')
        print(ctrl_novel_locs)
        print('test')
        print(novel_locs)
        assert len(set(ctrl_novel_locs)-set(novel_locs)) == 0
        assert len(novel_locs) == len(ctrl_novel_locs)

        # edge_df
        edge_df = sg.add_edge_coords()
        edge_df = edge_df.loc[edge_df.annotation == False]
        ctrl_novel_edges = [('chr2', 75, 65, '-', 'exon'),
                            ('chr2', 65, 50, '-', 'intron'),
                            ('chr2', 80, 60, '-', 'intron'),
                            ('chr2', 60, 50, '-', 'exon')]
        chroms = edge_df.chrom.tolist()
        v1s = edge_df.v1_coord.tolist()
        v2s = edge_df.v2_coord.tolist()
        strands = edge_df.strand.tolist()
        etypes = edge_df.edge_type.tolist()
        novel_edges = [(chrom,v1,v2,strand,etype) for chrom,v1,v2,strand,etype \
            in zip(chroms,v1s,v2s,strands,etypes)]
        print('control')
        print(ctrl_novel_edges)
        print('test')
        print(novel_edges)
        assert len(set(ctrl_novel_edges)-set(novel_edges)) == 0
        assert len(ctrl_novel_edges) == len(novel_edges)

    # add annotation to sg with data where data contains dupe transcript
    def test_add_annotation_sg_data_dupe_tid(self):
        sg = swan.SwanGraph()
        sg.add_transcriptome('files/test_novel_1.gtf')
        sg.add_annotation('files/test_known.gtf')

        # check with coord/chr bc of reindexing fuckery not being
        # remimplemented yet

        # t_df
        annot_tids = ['test1', 'test2', 'test4']
        assert all(sg.t_df.loc[annot_tids, 'annotation'])
        ctrl_novel_tids = ['test3', 'test5']
        novel_tids = sg.t_df.loc[sg.t_df.annotation == False, 'tid'].tolist()
        assert len(set(ctrl_novel_tids)-set(novel_tids)) == 0
        assert len(ctrl_novel_tids) == len(novel_tids)

        # make sure the novelty assignment worked
        annot_tids = sg.t_df.loc[sg.t_df.annotation == True, 'tid'].tolist()
        known_tids = sg.t_df.loc[sg.t_df.novelty == 'Known', 'tid'].tolist()
        assert set(annot_tids) == set(known_tids)

        annot_tids = sg.t_df.loc[sg.t_df.annotation == False, 'tid'].tolist()
        known_tids = sg.t_df.loc[sg.t_df.novelty == 'Undefined', 'tid'].tolist()
        assert set(annot_tids) == set(known_tids)

        # loc_df
        ctrl_novel_locs = [('chr2', 65)]
        temp = sg.loc_df[sg.loc_df.annotation == False]
        chroms = temp.chrom.tolist()
        coords = temp.coord.tolist()
        novel_locs = [(chrom, coord) for chrom, coord in zip(chroms, coords)]
        print('control')
        print(ctrl_novel_locs)
        print('test')
        print(novel_locs)
        assert len(set(ctrl_novel_locs)-set(novel_locs)) == 0
        assert len(novel_locs) == len(ctrl_novel_locs)

        # edge_df
        edge_df = sg.add_edge_coords()
        edge_df = edge_df.loc[edge_df.annotation == False]
        ctrl_novel_edges = [('chr2', 75, 65, '-', 'exon'),
                            ('chr2', 65, 50, '-', 'intron'),
                            ('chr2', 80, 60, '-', 'intron'),
                            ('chr2', 60, 50, '-', 'exon')]
        chroms = edge_df.chrom.tolist()
        v1s = edge_df.v1_coord.tolist()
        v2s = edge_df.v2_coord.tolist()
        strands = edge_df.strand.tolist()
        etypes = edge_df.edge_type.tolist()
        novel_edges = [(chrom,v1,v2,strand,etype) for chrom,v1,v2,strand,etype \
            in zip(chroms,v1s,v2s,strands,etypes)]
        print('control')
        print(ctrl_novel_edges)
        print('test')
        print(novel_edges)
        assert len(set(ctrl_novel_edges)-set(novel_edges)) == 0
        assert len(ctrl_novel_edges) == len(novel_edges)

def check_dfs(loc_df, ctrl_loc_df,
              edge_df, ctrl_edge_df,
              t_df, ctrl_t_df):

    # first order to make them comparable
    # sort all values by their IDs
    loc_df.sort_index(inplace=True)
    edge_df.sort_index(inplace=True)
    t_df.sort_index(inplace=True)
    ctrl_loc_df.sort_index(inplace=True)
    ctrl_edge_df.sort_index(inplace=True)
    ctrl_t_df.sort_index(inplace=True)

    # and order columns the same way
    ctrl_loc_df = ctrl_loc_df[loc_df.columns]
    ctrl_edge_df = ctrl_edge_df[edge_df.columns]
    ctrl_t_df = ctrl_t_df[t_df.columns]

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

    print('test')
    print(t_df)
    print('control')
    print(ctrl_t_df)
    assert (t_df == ctrl_t_df).all(axis=0).all()


def add_transcriptome_no_reorder_gtf(sg, gtf):
    t_df, exon_df, from_talon = swan.parse_gtf(gtf, False)
    sg.loc_df, sg.edge_df, sg.t_df = sg.create_dfs(t_df, exon_df, from_talon)
    return sg


def get_test_transcript_exon_dicts():

    # features:
    # exons that are "out of order"
    # two transcripts from the same gene
    # transcripts from - strand
    # transcripts from + strand
    transcripts = {
        'test1': { # + strand
            'gid': 'test1_gid',
            'gname': 'test1_gname',
            'tid': 'test1',
            'tname': 'test1_tname',
            'strand': '+',
            'novelty': 'Known',
            'exons': ['chr1_1_20_+_exon',
                      'chr1_35_40_+_exon', # out of order exon (+)
                      'chr1_25_30_+_exon',] # out of order exon (+)
        },
        'test2': { # - strand
            'gid': 'test2_gid',
            'gname': 'test2_gname',
            'tid': 'test2',
            'tname': 'test2_tname',
            'strand': '-',
            'novelty': 'Known',
            'exons': ['chr2_100_80_-_exon', # duplicated exon/locations
                      'chr2_75_60_-_exon', # same edge but intron vs. exon
                      'chr2_50_45_-_exon'] # duplicated exon/locations, same exon different strand
        },
        'test3': { # - strand
            'gid': 'test2_gid',
            'gname': 'test2_gname',
            'tid': 'test3',
            'tname': 'test3_tname',
            'strand': '-',
            'novelty': 'NIC',
            'exons': ['chr2_100_80_-_exon', # duplicated exon/locations
                      'chr2_50_45_-_exon', # out of order exon (-), duplicated exon/locations
                      'chr2_75_65_-_exon'] # out of order exon (-)
        },
        'test4': { # + strand
            'gid': 'test4_gid',
            'gname': 'test4_gname',
            'tid': 'test4',
            'tname': 'test4_tname',
            'strand': '+',
            'novelty': 'Known',
            'exons': ['chr2_45_50_+_exon'] # same exon different strand
        },
        'test5': { # - strand
            'gid': 'test2_gid',
            'gname': 'test2_gname',
            'tid': 'test5',
            'tname': 'test5_tname',
            'strand': '-',
            'novelty': 'ISM',
            'exons': ['chr2_100_80_-_exon', # duplicated exon/locations
                      'chr2_60_50_-_exon'] # same edge but intron vs. exon
        },
    }

    # features
    # locations that are shared across exons
    exons = {
        'chr1_1_20_+_exon': {
            'eid': 'chr1_1_20_+_exon',
            'chrom': 'chr1',
            'v1': 1,
            'v2': 20,
            'strand': '+',
            'edge_type': 'exon'
        },
        'chr1_25_30_+_exon': {
            'eid': 'chr1_25_30_+_exon',
            'chrom': 'chr1',
            'v1': 25,
            'v2': 30,
            'strand': '+',
            'edge_type': 'exon'
        },
        'chr1_35_40_+_exon': {
            'eid': 'chr1_35_40_+_exon',
            'chrom': 'chr1',
            'v1': 35,
            'v2': 40,
            'strand': '+',
            'edge_type': 'exon'
        },
        'chr2_100_80_-_exon': {
            'eid': 'chr2_100_80_-_exon',
            'chrom': 'chr2',
            'v1': 100,
            'v2': 80,
            'strand': '-',
            'edge_type': 'exon'
        },
        'chr2_75_60_-_exon': {
            'eid': 'chr2_75_60_-_exon',
            'chrom': 'chr2',
            'v1': 75,
            'v2': 60,
            'strand': '-',
            'edge_type': 'exon'
        },
        'chr2_50_45_-_exon': {
            'eid': 'chr2_50_45_-_exon',
            'chrom': 'chr2',
            'v1': 50,
            'v2': 45,
            'strand': '-',
            'edge_type': 'exon'
        },
        'chr2_75_65_-_exon': {
            'eid': 'chr2_75_60_-_exon',
            'chrom': 'chr2',
            'v1': 75,
            'v2': 65,
            'strand': '-',
            'edge_type': 'exon'
        },
        'chr2_45_50_+_exon': {
            'eid': 'chr2_45_50_+_exon',
            'chrom': 'chr2',
            'v1': 45,
            'v2': 50,
            'strand': '+',
            'edge_type': 'exon'
            },
        'chr2_100_80_-_exon': {
            'eid': 'chr2_100_80_-_exon',
            'chrom': 'chr2',
            'v1': 100,
            'v2': 80,
            'strand': '-',
            'edge_type': 'exon'
        },
        'chr2_60_50_-_exon': {
            'eid': 'chr2_60_50_-_exon',
            'chrom': 'chr2',
            'v1': 60,
            'v2': 50,
            'strand': '-'
        }
    }

    return transcripts, exons
