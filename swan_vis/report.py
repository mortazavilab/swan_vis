from swan_vis.utils import *
from fpdf import FPDF
import matplotlib.pyplot as pyplot

# report for genes - extension of FPDF class
class Report(FPDF):

	def __init__(self,
				 prefix,
				 report_type,
				 obs,
				 uns,
				 datasets,
				 groupby,
				 metadata_cols,
				 novelty=False,
				 layer='tpm',
				 cmap='Spectral_r',
				 g_min=None,
				 g_max=None,
				 include_qvals=False,
				 qval_df=None,
				 display_numbers=False,
				 t_disp='Transcript ID'):
		super().__init__(orientation='L')

		# change margins
		self.set_margins(0.5, 0.5, 0.5)

		# set report type, 'browser' or 'swan'
		if report_type != 'browser' and report_type != 'swan':
			raise Exception("Report type must be 'browser' or 'swan'.")
		else:
			self.report_type = report_type

		# booleans of what's in the report
		self.novelty = novelty
		self.include_qvals = include_qvals
		self.display_numbers = display_numbers

		# what we're plotting
		self.datasets = datasets
		self.n_dataset_cols = len(self.datasets)
		self.metadata_cols = metadata_cols
		self.layer = layer
		self.obs = obs
		self.uns = uns
		self.t_disp = t_disp
		if groupby == None:
			self.groupby = 'dataset'
		else:
			self.groupby = groupby
		self.qval_df = qval_df

		# prefix for files that we'll pull from
		self.prefix = prefix

		# color map
		self.g_min = g_min
		self.g_max = g_max
		try:
			self.cmap = plt.get_cmap(cmap)
		except:
			raise ValueError('Colormap {} not found'.format(cmap))

		# settings
		self.entry_height = 20

		# add extra room for the scale/colorbar if we're doing
		# the browser/heatmap version
		# if self.report_type == 'swan':
		# 	self.header_height = 10
		# if self.report_type == 'browser':
		# 	self.header_height = 20
		self.header_height = 20
		if self.metadata_cols:
			# self.meta_height = (self.header_height-1)/len(self.metadata_cols)
			self.meta_height = (self.header_height)/len(self.metadata_cols)


		# dataset width is contingent on # of datasets
		# as well as if we're including the novelty column
		if self.novelty:
			self.w_dataset = (146-25)/self.n_dataset_cols
		else:
			self.w_dataset = 146/self.n_dataset_cols

	# for each metadata column, color report header
	def color_header(self, data_col):
		x = self.get_x()
		y = 0.5

		for ind, meta_col in enumerate(self.metadata_cols):

			self.set_y(y+(ind*self.meta_height))
			self.set_x(x)

			# get meta_col category
			meta_cat = self.obs.loc[self.obs[self.groupby] == data_col, meta_col]
			meta_cat = meta_cat.unique().tolist()[0]

			# get the color and convert to rgb
			# source: https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python
			color = self.uns['{}_dict'.format(meta_col)][meta_cat]
			r = color[0]
			g = color[1]
			b = color[2]
			self.set_fill_color(r,g,b)

			self.cell(self.w_dataset, self.meta_height, fill=True, border=False)

	# header - should differ based on whether it's a browser report or
	# a swan report
	def header(self):
		self.set_font('Arial', 'B', 10)

		# transcript ID header
		self.cell(50, self.header_height, self.t_disp, border=True, align='C')

		# novelty header (if needed)
		if self.novelty:
			self.cell(25, self.header_height, 'Novelty', border=True, align='C')

		# dataset ID headers
		for col in self.datasets:

			# we're using colors to represent different metadata vars.
			if self.metadata_cols:
				self.color_header(col)

			# just using the dataset names
			else:
				self.cell(self.w_dataset, self.header_height, col,
						  border=True, align='C')

		# reset current location in case we colored metadata cols
		self.set_y(0.5)
		self.set_x(196.5)

		# reset to white
		self.set_fill_color(255,255,255)

		# in case we need to add the browser models
		browser_scale_x = self.get_x()
		browser_scale_y = self.get_y()

		# transcript model header
		self.cell(100, self.header_height, 'Transcript Model', border=True, align='C')

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
		if not self.novelty:
			self.set_x(77.5)
		else:
			self.set_x(90.5)
		self.image(self.prefix+'_colorbar_scale.png',
			w=90, h=13.57)

		# dataset color legends
		if self.metadata_cols:
			start_y = self.get_y()-12
			curr_y = start_y
			curr_x = 0
			self.set_font('Arial', '', 10)
			i = 0
			for meta_col in self.metadata_cols:
				if i != 0:
					curr_x += 30
					curr_y = start_y
				self.set_y(curr_y)
				self.set_x(curr_x)
				self.cell(32, 5, meta_col, align='L')
				curr_y = curr_y + 6


				for meta_cat in self.obs[meta_col].unique().tolist():

					# if i % 6 == 0 and i != 0:
					# 	if i == 12:
					# 		start_y = start_y + 12
					# 	curr_y = start_y
					# 	curr_x += 37

					self.set_y(curr_y)
					self.set_x(curr_x)

					# # get meta_col category
					# meta_cat = self.obs.loc[self.obs[self.groupby] == data_col, meta_col]
					# meta_cat = meta_cat.unique().tolist()[0]

					# get the color and convert to rgb
					# source: https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python
					color = self.uns['{}_dict'.format(meta_col)][meta_cat]
					r = color[0]
					g = color[1]
					b = color[2]
					self.set_fill_color(r,g,b)

					self.cell(5, 5, '', fill=True)
					self.cell(32, 5, meta_cat)

					curr_y = curr_y + 6
					i+=1
			self.set_font('Arial', '', 10)

	# add a transcript model to the report
	def add_transcript(self, entry, oname, tid):

		# corresponding entry from qval_df
		if self.include_qvals:
			qval_entry = self.qval_df.loc[self.qval_df.tid == entry.name].head(1).squeeze()

		# entries should not be bolded
		if self.include_qvals:
			if qval_entry.significant:
				self.set_font('Arial', 'B', 10)
				tid += '*'
			else:
				self.set_font('Arial', '', 10)
		else:
			self.set_font('Arial', '', 10)

		# tid
		self.cell(50, self.entry_height, tid, border=True, align='C')
		tid_x = self.get_x()
		tid_y = self.get_y()

		# add qvals if needed
		if self.include_qvals:
			curr_x = self.get_x()
			curr_y = self.get_y()
			if qval_entry.significant:
				self.set_font('Arial', 'B', 6)
			else:
				self.set_font('Arial', '', 6)
			self.set_y(tid_y+12)
			text = 'qval = {:.2e}'.format(qval_entry.qval)
			self.cell(50, 4, txt=text, border=False, align='C')
			self.set_font('Arial', '', 10)
			self.set_xy(curr_x, curr_y)

		# novelty, if necessary
		if self.novelty:
			self.cell(25, self.entry_height, entry['novelty'],
				border=True, align='C')

		# dataset columns
		for col in self.datasets:

			# color each heatmap cell
			norm_val = (entry[col]-self.g_min)/(self.g_max-self.g_min)
			color = self.cmap(norm_val)
			r = color[0]*255
			g = color[1]*255
			b = color[2]*255
			self.set_fill_color(r,g,b)
			border = False
			fill = True

			if self.display_numbers:
				text = str(round(entry[col],2))
				if (r*0.299 + g*0.587 + b*0.114) > 100:
					text_r = text_g = text_b = 0
				else:
					text_r = text_g = text_b = 255
				self.set_text_color(text_r, text_g, text_b)
			else:
				text = ''

			self.cell(self.w_dataset, self.entry_height, text,
				border=border, align='C', fill=fill)
			# reset text color
			self.set_text_color(0,0,0)

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
