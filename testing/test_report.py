import pytest
import sys
import numpy as np
import swan_vis as swan

class TestReport(object):

	# tests to see if gen_report throws errors correctly
	def test_novelty_error(self):
		print('testing to see if gen_report throws an error '
			  'when you ask for the novelty column with no novelty')
		sg = gen_toy_sg()
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG03', 'figures/ensg03', novelty=True)
		print(str(excinfo.value))
		assert 'No novelty information' in str(excinfo.value)

	def test_no_gene(self):
		print('testing to see if gen_report throws an error '
			  'wher you ask for a gene not in the graph')
		sg = gen_toy_sg()
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG04', 'figures/ensg04')
		print(str(excinfo.value))
		assert 'Gene ENSG04 not found' in str(excinfo.value)

	def test_incompatible_plot_args(self):

		print('indicate_novel and indicate_dataset used together')
		sg = gen_toy_sg()
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG03', 'figures/ensg03', indicate_dataset=True, indicate_novel=True)
		print(str(excinfo.value))
		assert 'choose either indicate_dataset' in str(excinfo.value)

		print('indicate_novel without annotation')
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG03', 'figures/ensg03', indicate_novel=True)
		print(str(excinfo.value))
		assert 'Annotation data not present in graph' in str(excinfo.value)

		print('indicate_dataset without that dataset name')
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG03', 'figures/ensg03', indicate_dataset='b')
		print(str(excinfo.value))
		assert 'Dataset b not present' in str(excinfo.value)

		print('browser and indicate_dataset')
		with pytest.raises(Exception) as excinfo:
			sg.gen_report('ENSG03', 'figures/ensg03', browser=True, indicate_dataset='a')
		print(str(excinfo.value))
		assert 'or indicate_dataset' in str(excinfo.value)

		# this is what I'm referring to. Need to address this issue.
		# sg.gen_report('ENSG03', 'figures/ensg03')


def gen_toy_sg():
	sg = swan.SwanGraph()
	sg.add_dataset('a', 'input_files/annot.gtf')
	return sg 




