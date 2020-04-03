from swan_vis.utils import *
from fpdf import FPDF
import matplotlib.pyplot as pyplot

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, prefix, report_type, report_cols, header_cols, heatmap=False):
		super().__init__(orientation='L')

		# change margins
		self.set_margins(0.5, 0.5, 0.5)

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# does this report include a heatmap?
		self.heatmap = heatmap

		# the columns that we'll include
		self.report_cols = report_cols
		self.header_cols = header_cols
		self.n_dataset_cols = len(self.report_cols)
		
		# prefix for files that we'll pull from 
		self.prefix = prefix

		# color map in case we're making a heatmap
		self.cmap = plt.get_cmap('Spectral_r')

	# header - should differ based on whether it's a browser report or
	# a swan report
	def header(self):
		self.set_font('Arial', 'B', 10)

		# add extra room for the scale/colorbar if we're doing 
		# the browser/heatmap version
		if self.report_type == 'swan':
			header_height = 10
		if self.report_type == 'browser' or self.heatmap == True:
			header_height = 20

		# size the dataset IDs differently if we need to include
		# colorbar below
		if self.heatmap == False:
			dataset_height = header_height
		elif self.heatmap == True:
			dataset_height = 10
		
		# transcript ID header
		self.cell(50, header_height, 'Transcript ID', border=True, align='C')

		# in case we need to add the colorbar
		colorbar_x = self.get_x() + float(45/2) + 4.5
		colorbar_y = self.get_y()

		# dataset ID headers
		w_dataset = 146/self.n_dataset_cols
		for col in self.header_cols:
			self.cell(w_dataset, dataset_height, col,
					  border=True, align='C')

		# in case we need to add the browser models
		browser_scale_x = self.get_x()
		browser_scale_y = self.get_y()

		# transcript model header
		self.cell(100, header_height, 'Transcript Model', border=True, align='C')

		# add colorbar if we need to 
		if self.heatmap == True: 
			self.set_xy(50.5, 10.5)
			self.cell(146, 10, border=True)
			self.image(self.prefix+'_colorbar_scale.png',
				x=colorbar_x, 
				y=colorbar_y+10,
				w=90, h=135/14)

		# add scale if we're doing browser models
		if self.report_type == 'browser':
			self.image(self.prefix+'_browser_scale.png',
					   x=browser_scale_x, 
					   y=browser_scale_y+12,
					   w=100, h=50/7)
		self.ln()

	# add a transcript model to the report
	def add_transcript(self, entry, oname):

		# dynamically size dataset cols
		w_dataset = 146/self.n_dataset_cols

		self.set_font('Arial', '', 10)
		self.cell(50, 20, entry['tid'], border=True, align='C')
		for col in self.report_cols:

			# heat map coloring
			if self.heatmap:
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
			self.cell(w_dataset, 20, text, border=border, align='C', fill=fill)	
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