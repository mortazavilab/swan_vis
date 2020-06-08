import sys
import os
from collections import defaultdict
import swan_vis as swan

class TestUtils(object):
# tests random utilities, from the utils.py file
	def test_process_abundance_file(self):
		file = 'input_files/test_abundance.tsv'
		df = swan.process_abundance_file(file, ['count_1a', 'count_2a'],
			tid_col='annot_transcript_id')

		print(df)
		
		# check counts
		test_pairs = df.apply(lambda x: (x.tid, x.counts), axis=1)
		control_pairs = ((0,3),(1,7),(2,11))
		print('test pairs')
		print(test_pairs)
		check_pairs(control_pairs, test_pairs)

		# check TPM
		tpm_sum = df.tpm.sum()
		print(tpm_sum)
		assert (10**6-1) < tpm_sum < (10**6+1)

	# def test_get_tpm_cols(self):
	

def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control
