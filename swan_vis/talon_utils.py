import copy
import itertools
import operator
import sqlite3

# All functions in this file written by Dana Wyman for TALON, and
# adapted to interface with TALON dbs

# Converts input to string that can be used for IN database query
def format_for_in(l):
	if type(l) is tuple:
		l = list(l)
	if type(l) is str:
		l = [l]
	return "(" + ','.join(['"' + str(x) + '"' for x in l]) + ")"

def fetch_all_transcript_gene_pairs(cursor):
	""" Return gene_ID - transcript_ID tuples from database """

	query = """ SELECT gene_ID, transcript_ID FROM transcripts """
	cursor.execute(query)

	pairs = cursor.fetchall()
	return pairs

def fetch_all_datasets(cursor):
	""" Return a list of all datasets in database """
	cursor.execute("SELECT dataset_name FROM dataset")
	datasets = [str(x[0]) for x in cursor.fetchall()]
	return datasets

def parse_pass_list(pass_list_file):
	""" From the pass_list file, obtain a list of acccepted gene and
		transcript IDs tuples"""
	pass_list = set()
	with open(pass_list_file, 'r') as f:
		for line in f:
			line = line.strip()
			fields = line.split(",")
			gene_ID = fields[0]
			transcript_ID = fields[1]
			try:
				pass_list.add((int(gene_ID), int(transcript_ID)))
			except:
				raise ValueError("Gene/Transcript IDs in pass_list must be integer TALON IDs")
	return pass_list

def parse_datasets(dataset, cursor):
	""" From the dataset file, obtain a list of acccepted dataset names"""
	db_datasets = fetch_all_datasets(cursor)
	if dataset not in db_datasets:
		raise ValueError("Dataset name '%s' not found in database" % dataset)
	return dataset

def handle_filtering(database, observed, pass_list_file):
	""" Determines which transcripts to allow in the analysis. This can be done
		in two different ways. If no pass_list is included, then all of the
		transcripts in the database are included (modified by 'observed'
		option). If a pass_list is provided, then transcripts on that list
		will be included (modified by 'observed' option). This can be
		tuned further by providing a dataset file, but this is optional. """

	conn = sqlite3.connect(database)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()

	# Get list of datasets to use in run
	# if dataset != None:
	# 	datasets = parse_datasets(dataset, cursor)
	# elif observed == True:
	if observed: 
		datasets = fetch_all_datasets(cursor)
	else:
		datasets = None
	# datasets = fetch_all_datasets(cursor)

	# Get initial transcript pass_list
	if pass_list_file != None:
		pass_list = parse_pass_list(pass_list_file)
	else:
		pass_list = fetch_all_transcript_gene_pairs(cursor)

	if datasets != None:
		# Limit the pass_list to transcripts detected in the datasets
		transcripts = [ x[1] for x in pass_list ]
		transcript_str = format_for_in(transcripts)
		dataset_str = format_for_in(datasets)

		query = """ SELECT DISTINCT gene_ID, transcript_ID
					FROM observed
					WHERE transcript_ID IN %s
					AND dataset in %s """
		cursor.execute(query % (transcript_str, dataset_str))
		pass_list = cursor.fetchall()

	conn.close()
	return pass_list

def get_gene_transcript_map(db, pass_list):
	""" Creates a dictionary mapping gene IDs to the transcripts that belong to
		them. The columns in each tuple are:
			0: gene ID
			1: transcript ID
			2: chromosome
			3: start position (min of 5' and 3')
			4: end position (max of 5' and 3')
			5: strand
			6: edge path
			7. n_exons
	"""

	conn = sqlite3.connect(db)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()
	pass_list_string = "(" + ','.join([str(x) for x in pass_list]) + ")"
	query = """
			SELECT
				t.gene_ID,
				t.transcript_ID,
				loc1.chromosome,
				MIN(loc1.position,loc2.position) AS min_pos,
				MAX(loc1.position,loc2.position) AS max_pos,
				genes.strand,
				t.jn_path,
				t.start_exon,
				t.end_exon,
				t.n_exons
			FROM transcripts t
			LEFT JOIN location loc1 ON t.start_vertex = loc1.location_ID
			LEFT JOIN location loc2 ON t.end_vertex = loc2.location_ID
			LEFT JOIN genes ON t.gene_ID = genes.gene_ID
			WHERE t.transcript_ID IN """ + pass_list_string
	cursor.execute(query)
	transcript_tuples = cursor.fetchall()

	# Sort based on gene ID
	sorted_transcript_tuples = sorted(transcript_tuples, key=lambda x: x["gene_ID"])

	gene_groups = {}
	for key,group in itertools.groupby(sorted_transcript_tuples,operator.itemgetter(0)):
		# sort by transcript start position
		gene_groups[key] = sorted(list(group), key=lambda x: x["min_pos"])
	conn.close()

	return gene_groups

def get_annotations(database, feat_type, pass_list=None):
	"""
		Extracts annotations from the gene/transcript/exon annotation table of
		the database (depending on choice of feat_type).

		Returns:
			annotation_dict: dictionary data structure in which the keys are
							 gene/transcript/exon TALON IDs (depending on
							 choice of feat_type) and the value is a list of
							 annotation tuples.
	"""
	# fetch the annotations
	conn = sqlite3.connect(database)
	cursor = conn.cursor()

	table_name = feat_type + "_annotations"

	if pass_list == None:
		query = "SELECT * FROM " + table_name
	else:
		pass_list_string = "(" + ','.join([str(x) for x in pass_list]) + ")"
		query = "SELECT * FROM " + table_name + " WHERE ID IN " + pass_list_string

	cursor.execute(query)
	annotation_tuples = cursor.fetchall()

	# sort based on ID
	sorted_annotations = sorted(annotation_tuples, key=lambda x: x[0])

	# group by ID and store in a dictionary
	ID_groups = {}
	for key,group in itertools.groupby(sorted_annotations,operator.itemgetter(0)):
		ID_groups[key] = list(group)

	return ID_groups

def get_gene_2_transcripts(database, pass_list):
	""" Creates a dictionary mapping gene IDs to the transcripts that belong to
		them. The columns in each tuple are:
			0: gene ID
			1: transcript ID
			2: chromosome
			3: start position (min of 5' and 3')
			4: end position (max of 5' and 3')
			5: strand
			6: edge path
			7. n_exons
	"""

	conn = sqlite3.connect(database)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()
	pass_list_string = "(" + ','.join([str(x) for x in pass_list]) + ")"
	query = """
			SELECT
				t.gene_ID,
				t.transcript_ID,
				loc1.chromosome,
				MIN(loc1.position,loc2.position) AS min_pos,
				MAX(loc1.position,loc2.position) AS max_pos,
				genes.strand,
				t.jn_path,
				t.start_exon,
				t.end_exon,
				t.n_exons
			FROM transcripts t
			LEFT JOIN location loc1 ON t.start_vertex = loc1.location_ID
			LEFT JOIN location loc2 ON t.end_vertex = loc2.location_ID
			LEFT JOIN genes ON t.gene_ID = genes.gene_ID
			WHERE t.transcript_ID IN """ + pass_list_string
	cursor.execute(query)
	transcript_tuples = cursor.fetchall()

	# Sort based on gene ID
	sorted_transcript_tuples = sorted(transcript_tuples, key=lambda x: x["gene_ID"])

	gene_groups = {}
	for key,group in itertools.groupby(sorted_transcript_tuples,operator.itemgetter(0)):
		# Sort by transcript start position
		gene_groups[key] = sorted(list(group), key=lambda x: x["min_pos"])
	conn.close()

	return gene_groups

def fetch_exon_locations(database):
	""" Queries the database to create a dictionary mapping exon IDs to
		the chromosome, start, end, and strand of the exon """

	conn = sqlite3.connect(database)
	cursor = conn.cursor()

	query = """
			SELECT
				e.edge_ID,
				loc1.chromosome,
				MIN(loc1.position,loc2.position),
				MAX(loc1.position,loc2.position),
				e.strand
			 FROM edge e
			 LEFT JOIN location loc1 ON e.v1 = loc1.location_ID
			 LEFT JOIN location loc2 ON e.v2 = loc2.location_ID
			 WHERE e.edge_type = 'exon';"""

	cursor.execute(query)
	exon_location_tuples = cursor.fetchall()

	# Create dictionary
	exon_locations = {}
	for loc_tuple in exon_location_tuples:
		exon_ID = loc_tuple[0]
		exon_locations[exon_ID] = loc_tuple[1:]

	conn.close()
	return exon_locations

# def check_annot_validity(annot, database):
#     """ Make sure that the user has entered a correct annotation name """

#     conn = sqlite3.connect(database)
#     cursor = conn.cursor()

#     cursor.execute("SELECT DISTINCT annot_name FROM gene_annotations")
#     annotations = [str(x[0]) for x in cursor.fetchall()]
#     conn.close()

#     if "TALON" in annotations:
#         annotations.remove("TALON")

#     if annot not in annotations:
#         message = "Annotation name '" + annot + \
#                   "' not found in this database. Try one of the following: " + \
#                   ", ".join(annotations)
#         raise Exception(message)

#     return annot
