from utils import *
from fpdf import FPDF

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, prefix, report_type, datasets):
		super().__init__(orientation='L')

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# the dataset columns that we'll include
		self.datasets = datasets

		# prefix for files that we'll pull from 
		self.prefix = prefix

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

		# in case we need to add the browser models
		x = self.get_x()
		y = self.get_y()

		self.cell(100, header_height, 'Transcript Model', border=True, align='C')

		# add scale if we're doing browser models
		if self.report_type == 'browser':
			self.image(self.prefix+'_browser_scale.png',
					   x=x, y=y+12, w=100, h=50/7)
		self.ln()

	# add a transcript model to the report
	def add_transcript(self, entry, oname):
		self.set_font('Arial', '', 10)
		self.cell(50, 20, entry['tid'], border=True, align='C')
		for d in self.datasets:
			self.cell(50, 20, str(round(entry['{}_tpm'.format(d)],2)), border=True, align='C')
		x = self.get_x()
		y = self.get_y()
		self.cell(100, 20, '', border=True)
		self.image(oname, x=x, y=y, w=100, h=20)
		self.ln()


	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))


