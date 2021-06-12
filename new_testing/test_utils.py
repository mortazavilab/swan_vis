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
