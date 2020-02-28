from swan_vis.utils import *
from fpdf import FPDF
import matplotlib.pyplot as pyplot

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, prefix, report_type, report_cols, header_cols):
		super().__init__(orientation='L')

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# the columns that we'll include
		self.report_cols = report_cols
		self.header_cols = header_cols
		
		# prefix for files that we'll pull from 
		self.prefix = prefix

		# color map in case we're making a heatmap
		self.cmap = plt.get_cmap('Spectral_r')

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
		for col in self.header_cols:
			self.cell(25, header_height, col,
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
	def add_transcript(self, entry, oname, heatmap=False):
		self.set_font('Arial', '', 10)
		self.cell(50, 20, entry['tid'], border=True, align='C')
		for col in self.report_cols:

			# heat map coloring
			if heatmap:
				color = self.cmap(entry[col])
				r = color[0]*255
				b = color[1]*255
				g = color[2]*255
				self.set_fill_color(r,b,g)
				border = False
				fill = True
				text = ''
			# TPM	
			elif '_tpm' in col:
				text = str(round(entry[col],2))
				border = True
				fill = False
			# presence/absence
			else:
				text = entry[col]
				if text == True:
					text = 'Yes'
				elif text == False:
					text = 'No'
				border = True
				fill = False 
			self.cell(25, 20, text, border=border, align='C', fill=fill)	
		x = self.get_x()
		y = self.get_y()

		# reset color to white
		self.set_fill_color(255,255,255)

		# embed transcript model
		self.cell(100, 20, '', border=True)
		self.image(oname, x=x, y=y, w=100, h=20)
		self.ln()


	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))


