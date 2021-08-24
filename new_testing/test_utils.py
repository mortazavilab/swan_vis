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

    # test reformat_talon_abundance
    def test_reformat_talon_abundance(self):
        test = swan.reformat_talon_abundance('files/test_ab_talon_1.tsv')
        data = [['test1', 5, 5],
                ['test2', 10, 0],
                ['test3', 0, 10],
                ['test4', 10, 10],
                ['test5', 5, 5],
                ['test8', 4, 5]]
        cols = ['annot_transcript_id', 'dataset1', 'dataset2']
        ctrl = pd.DataFrame(data=data, columns=cols)
        assert test.equals(ctrl)

    # test reorder_exons - rev
    def test_reorder_exons_2(self):
        test_exons = ['chr2_100_80_-_exon', 'chr2_50_45_-_exon',
                     'chr2_75_60_-_exon']
        test_exons = swan.reorder_exons(test_exons)
        ctrl_exons = ['chr2_100_80_-_exon', 'chr2_75_60_-_exon',
                      'chr2_50_45_-_exon']
        assert test_exons == ctrl_exons

    # test reorder_exons - fwd
    def test_reorder_exons_1(self):
        test_exons = ['chr1_1_20_+_exon', 'chr1_35_40_+_exon',
                     'chr1_25_30_+_exon']
        test_exons = swan.reorder_exons(test_exons)
        ctrl_exons = ['chr1_1_20_+_exon', 'chr1_25_30_+_exon',
                      'chr1_35_40_+_exon']


    # test find_edge_start_stop - minus strand
    def test_find_edge_start_stop_2(self):
        start, stop = swan.find_edge_start_stop(1,2, '-')
        assert 1 == stop
        assert 2 == start

    # test find_edge_start_stop - plus strand
    def test_find_edge_start_stop_1(self):
        start, stop = swan.find_edge_start_stop(1,2, '+')
        assert 1 == start
        assert 2 == stop

    # test gtf_or_db - neither
    def test_gtf_or_db_2(self):
        fname = 'test/test.txt'
        with pytest.raises(Exception) as e:
            swan.gtf_or_db(fname)
        assert 'must be' in str(e.value)

    # test gtf_or_db - gtf
    def test_gtf_or_db_2(self):
        fname = 'test/test.gtf'
        assert swan.gtf_or_db(fname) == 'gtf'

    # test gtf_or_db - db
    def test_gtf_or_db_1(self):
        fname = 'test/test.db'
        assert swan.gtf_or_db(fname) == 'db'

    # test check_file_loc - does not exist
    def test_check_file_loc_2(self):
        with pytest.raises(Exception) as e:
            swan.check_file_loc('files/doesnotexist.gtf', 'gtf')
        assert 'not found' in str(e.value)

    # test check_file_loc - exists
    def test_check_file_loc_1(self):
        swan.check_file_loc('files/test_full.gtf', 'gtf')

    # test check_dir_loc - does not exist
    def test_check_dir_loc_2(self):
        with pytest.raises(Exception) as e:
            swan.check_dir_loc('beepboopdonotexist/helloworld/')
        assert 'not found' in str(e.value)

    # test check_dir_loc - exists
    def test_check_dir_loc_1(self):
        swan.check_dir_loc('files/')

    # test get_fields
    def test_get_fields(self):
        fields = 'gene_id "test1_gid"; gene_name "test1_gname"; transcript_id "test1"; transcript_name "test1_tname"; transcript_status "KNOWN";'
        test_fields = swan.get_fields(fields)
        ctrl_fields = {'gene_id': 'test1_gid', 'gene_name': 'test1_gname',
                       'transcript_id': 'test1', 'transcript_name': 'test1_tname',
                       'transcript_status': 'KNOWN'}
        assert test_fields == ctrl_fields

        fields = 'gene_name "test1_gname"; transcript_id "test1"; transcript_name "test1_tname"; transcript_status "KNOWN";'
        test_fields = swan.get_fields(fields)
        ctrl_fields = {'gene_id': 'NULL', 'gene_name': 'test1_gname',
                       'transcript_id': 'test1', 'transcript_name': 'test1_tname',
                       'transcript_status': 'KNOWN'}
        assert test_fields == ctrl_fields

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
