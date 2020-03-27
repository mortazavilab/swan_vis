from swan_vis.utils import *
from fpdf import FPDF
import matplotlib.pyplot as pyplot

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self, prefix, gid, report_type, report_cols, header_cols, heatmap=False, num_transcripts=[]):
		super().__init__(orientation='L')

		# change margins
		self.set_margins(0.5, 0.5, 0.5)

		# gene for this report
		self.gid = gid

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else: 
			self.report_type = report_type

		# does this report include a heatmap?
		self.heatmap = heatmap
		self.num_transcripts = num_transcripts

		# the columns that we'll include
		self.report_cols = report_cols
		self.header_cols = header_cols
		self.n_dataset_cols = len(self.report_cols)
		
		# prefix for files that we'll pull from 
		self.prefix = prefix

		# color map in case we're making a heatmap
		self.cmap = plt.get_cmap('Spectral_r')

		# number of pages we've added
		self.num_pages = 0

	# header - should differ based on whether it's a browser report or
	# a swan report
	def header(self):

		self.set_font('Arial', 'B', 10)

		# increment the number of pages we've made
		# get the number of transcripts on this page
		num_transcripts = self.num_transcripts[self.num_pages]
		self.num_pages += 1

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
		colorbar_x = self.get_x() 
		colorbar_y = self.get_y()

		# dataset ID headers
		w_dataset = 146/self.n_dataset_cols
		for col in self.header_cols:
			self.cell(w_dataset, dataset_height, col,
					  border=True, align='C')

		# transcript model header
		self.cell(100, header_height, 'Transcript Model', border=True, align='C')

		# add colorbar and heatmap if we need to 
		if self.heatmap == True: 

			# colorbar
			self.set_xy(50.5, 10.5)
			self.cell(146, 10, border=True)
			self.image(self.prefix+'_colorbar_scale.png',
				x=colorbar_x+float(45/2) + 4.5, 
				y=colorbar_y+10,
				w=90, h=135/14)

			# heatmap
			heatmap_name = '{}_{}_heatmap_{}.png'.format(self.prefix, self.gid, self.num_pages)
			heatmap_height = num_transcripts*20
			print(num_transcripts)
			print(heatmap_height)
			self.image(heatmap_name,
				x=50.5,
				y=20.5,
				w=146,
				h=heatmap_height)

		# in case we need to add the browser scale
		browser_scale_x = self.get_x()
		browser_scale_y = self.get_y()

		# add scale if we're doing browser models
		if self.report_type == 'browser':
			self.image(self.prefix+'_browser_scale.png',
					   x=browser_scale_x, 
					   y=browser_scale_y+12,
					   w=100, h=50/7)
		self.ln()

	# add a transcript model to the report
	def add_transcript(self, entry, oname):

		self.set_font('Arial', '', 10)
		self.cell(50, 20, entry['tid'], border=True, align='C')

		# tpm and presence/absence
		if not self.heatmap:

			# dynamically size dataset cols
			w_dataset = 146/self.n_dataset_cols

			for col in self.report_cols:

				# fill with tpm number
				if '_tpm' in col:
					text = str(round(entry[col],2))

				# presence/absence
				else:
					text = entry[col]
					if text == True:
						text = 'Yes'
					elif text == False:
						text = 'No'
				self.cell(w_dataset, 20, text, border=True, align='C')	

		# heatmap - don't need columns
		else: 	
			self.cell(146, 20)

		x = self.get_x()
		y = self.get_y()

		# embed transcript model
		self.cell(100, 20, '', border=True)
		self.image(oname, x=x, y=y, w=100, h=20)
		self.ln()


	# writes the pdf to file with the correct formatting
	def write_pdf(self, file):
		with open(file, 'wb') as outfile:
			outfile.write(self.output(dest='S').encode('latin-1'))


