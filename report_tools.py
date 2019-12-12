from utils import *
import plotting_tools as pt
import SpliceGraph as sg 
import PlottedGraph as pg
from fpdf import FPDF

# extension of the FPDF class
class PDF(FPDF):
	def add_datasets(self, datasets):
		self.datasets = datasets

	def header(self):

		# # normie
		# self.set_font('Arial', 'B', 10)
		# self.cell(50, 10, 'Transcript ID', border=True, align='C')
		# for d in self.datasets:
		# 	self.cell(25, 10, '{} TPM'.format(d), border=True, align='C')
		# self.cell(100, 10, 'Transcript Model', border=True, align='C')
		
		# add scale for browser
		self.set_font('Arial', 'B', 10)
		self.cell(50, 20, 'Transcript ID', border=True, align='C')
		for d in self.datasets:
			self.cell(25, 20, '{} TPM'.format(d), border=True, align='C')
		x = self.get_x()
		y = self.get_y()
		self.cell(100, 20, 'Transcript Model', border=True, align='C')
		self.image('figures/scale.png', x=x, y=y+12, w=100, h=7)

		self.ln()

	def add_transcript(self, entry, oname):
		self.set_font('Arial', '', 10)
		self.cell(50, 20, entry['tid'], border=True, align='C')
		for d in self.datasets:
			field = 'tpm_{}'.format(d)
			self.cell(25, 20, str(truncate(entry[field], 3)), border=True, align='C')
		
		x = self.get_x()
		y = self.get_y()
		self.cell(100, 20, '', border=True)
		self.image(oname, x=x, y=y, w=100, h=20)
		self.ln()

	def write_pdf(self, file):
		print('hello world 2')
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))

# generates a report for each transcript model in the graph
def gen_report(splice_graph, args, oprefix, browser=False, order='expression'):

	splice_graph.order_transcripts(order)

	count_fields = sg.get_tpm_fields(splice_graph.t_df)
	dataset_names = [c_f.split('tpm_')[-1] for c_f in count_fields]

	pdf = PDF()
	pdf.set_left_margin(5)
	pdf.add_datasets(dataset_names)
	pdf.add_page()

	# pt.plot_gene_scale(splice_graph)

	# plot each transcript
	if not browser:
		pt.plot_each_transcript(splice_graph, args, oprefix, browser=False)
	else:
		pt.plot_each_transcript(splice_graph, args, oprefix, browser=True)

	for tid in splice_graph.t_df.tid.tolist():
		entry = splice_graph.t_df.loc[tid]
		if not browser:
			oname = '{}_{}.png'.format(oprefix, tid) 
		else:
			oname = '{}_{}_browser.png'.format(oprefix, tid)
		pdf.add_transcript(entry, oname)

	if not browser:
		pdf.write_pdf('{}_report.pdf'.format(oprefix))
	else:
		print('hello world')
		pdf.write_pdf('{}_report_browser.pdf'.format(oprefix))

