from swan_vis.utils import *
from fpdf import FPDF
import matplotlib.pyplot as pyplot
import seaborn as sns

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self,
				 prefix,
				 report_type, 
				 datasets,
				 data_type, 
				 novelty=False,
				 heatmap=False,
				 include_qvals=False):
		super().__init__(orientation='L')

		# change margins
		self.set_margins(0.5, 0.5, 0.5)

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# booleans of what's in the report
		self.heatmap = heatmap
		self.novelty = novelty 
		self.include_qvals = include_qvals

		# the columns that we'll include
		self.datasets = datasets
		self.report_cols = self.get_report_cols(data_type)
		self.n_dataset_cols = len(self.report_cols)

		# colors for each of the datasets
		colors = sns.color_palette('nipy_spectral',
			len(self.datasets))
		self.dataset_colors = []
		for i, d in enumerate(self.datasets):
			color = colors[i]
			r = color[0]*255	
			g = color[1]*255
			b = color[2]*255
			self.dataset_colors.append((r,g,b))
		
		# prefix for files that we'll pull from 
		self.prefix = prefix

		# color map in case we're making a heatmap
		self.cmap = plt.get_cmap('Spectral_r')

		# settings
		self.entry_height = 20

		# add extra room for the scale/colorbar if we're doing 
		# the browser/heatmap version
		if self.report_type == 'swan':
			self.header_height = 10
		if self.report_type == 'browser':
			self.header_height = 20

		# dataset width is contingent on # of datasets
		# as well as if we're including the novelty column
		if self.novelty:
			self.w_dataset = (146-25)/self.n_dataset_cols
		else:
			self.w_dataset = 146/self.n_dataset_cols

	# grab the relevant columns associated with the report type
	def get_report_cols(self, data_type):
		if not data_type:
			return self.datasets
		elif data_type == 'tpm':
			report_cols = [d+'_tpm' for d in self.datasets]
		elif data_type == 'heatmap':
			report_cols = [d+'_norm_log_tpm' for d in self.datasets]
		return report_cols

	# header - should differ based on whether it's a browser report or
	# a swan report
	def header(self):
		self.set_font('Arial', 'B', 10)
		
		# transcript ID header
		self.cell(50, self.header_height, 'Transcript ID', border=False, align='C')

		# novelty header (if needed)
		if self.novelty:
			self.cell(25, self.header_height, 'Novelty', border=False, align='C')

		# dataset ID headers
		for color, col in zip(self.dataset_colors, self.datasets):
			r = color[0]
			g = color[1]
			b = color[2]
			self.set_fill_color(r, g, b)
			self.cell(self.w_dataset, self.header_height, '',
					  fill=True, align='C')

		# in case we need to add the browser models
		browser_scale_x = self.get_x()
		browser_scale_y = self.get_y()

		# transcript model header
		self.cell(100, self.header_height, 'Transcript Model', border=False, align='C')

		# add scale if we're doing browser models
		if self.report_type == 'browser':
			self.image(self.prefix+'_browser_scale.png',
					   x=browser_scale_x, 
					   y=browser_scale_y+12,
					   w=100, h=50/7)
		self.ln()

	# footer - just add the colorbar if we're using the
	# heatmap option
	def footer(self):
		if self.heatmap:
			if not self.novelty:
				self.set_x(77.5)
			else: 
				self.set_x(90.5)
			self.image(self.prefix+'_colorbar_scale.png',
				w=90, h=135/14)

		# dataset color legends
		start_y = self.get_y()-7
		curr_y = start_y
		curr_x = 0
		self.set_font('Arial', '', 8)
		for i, d in enumerate(self.datasets):

			if i % 6 == 0 and i != 0:
				if i == 12:
					start_y = start_y + 12
				curr_y = start_y
				curr_x += 37
 
			self.set_y(curr_y)
			self.set_x(curr_x)

			color = self.dataset_colors[i]
			r = color[0]
			g = color[1]
			b = color[2]
			self.set_fill_color(r, g, b)

			self.cell(5, 5, '', fill=True)
			self.cell(32, 5, d)

			curr_y = curr_y + 6
		self.set_font('Arial', '', 10)


	# add a transcript model to the report
	def add_transcript(self, entry, oname):

		# entries should not be bolded
		if self.include_qvals:
			if entry.significant:
				self.set_font('Arial', 'B', 10)
			else:
				self.set_font('Arial', '', 10)
		else:
			self.set_font('Arial', '', 10)

		# tid
		self.cell(50, self.entry_height, entry['tid'], border=True, align='C')
		tid_x = self.get_x()
		tid_y = self.get_y()

		# add qvals if needed
		if self.include_qvals:
			curr_x = self.get_x()
			curr_y = self.get_y()
			if entry.significant:
				self.set_font('Arial', 'B', 6)
			else:
				self.set_font('Arial', '', 6)
			self.set_y(tid_y+12)
			text = 'qval = {:.2e}'.format(entry.qval)
			self.cell(50, 4, txt=text, border=False, align='C')
			self.set_font('Arial', '', 10)
			self.set_xy(curr_x, curr_y)

		# novelty, if necessary
		if self.novelty: 
			self.cell(25, self.entry_height, entry['novelty'],
				border=True, align='C')

		# dataset columns
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
			self.cell(self.w_dataset, self.entry_height, text,
				border=border, align='C', fill=fill)	
		x = self.get_x()
		y = self.get_y()

		# reset color to white
		self.set_fill_color(255,255,255)

		# embed transcript model
		self.cell(100, self.entry_height, '', border=True)
		self.image(oname, x=x, y=y, w=100, h=20)
		self.ln()

	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))
