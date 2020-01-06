from utils import *
from fpdf import FPDF

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, report_type, datasets):
		super().__init__(orientation='L')

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# the dataset columns that we'll include
		self.datasets = datasets

	# header - should differ based on whether it's a browser report or
	# a swan report
	def header(self):
		self.set_font('Arial', 'B', 10)

		# add extra room for the scale if we're doing the browser version
		if self.report_type == 'swan':
			header_height = 10
		elif self.report_type == 'browser':
			header_height = 20

		self.cell(50, header_height, 'Transcript ID', border=True, align='C')
		for dataset in self.datasets:
			self.cell(50, header_height, '{} TPM'.format(dataset),
					  border=True, align='C')
		self.cell(100, header_height, 'Transcript Model', border=True, align='C')
		self.ln()


	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))


