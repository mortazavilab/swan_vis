
import pytest
import sys
sys.path.append('../utils/')
from utils import *
# @pytest.mark.unit

class TestTest(object):
    def test_1(self):
        """ Plus-strand upstream
        """
        # print('hello world')
        sg = SpliceGraph(talon_db='input_files/combine.db') 
        assert len(sg.t_df.index) == 3
        assert len(sg.loc_df.index) == 12
        assert len(sg.edge_df.index) == 11