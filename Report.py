from utils import *
from SpliceGraph import SpliceGraph
from PlottedGraph import PlottedGraph
from fpdf import FPDF

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, report_type):
		super().__init__()

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.type = report_type

	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))


# generate a report for each transcript model of a specific gene
def gen_report():
