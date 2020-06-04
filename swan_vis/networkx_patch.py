def main():
	try:
	  import networkx as nx 
	  import shutil
	  import swan_vis as swan
	  nx_loc = nx.__file__
	  nx_loc = '/'.join(nx_loc.split('/')[:-1])
	  nx_loc += '/drawing/nx_pylab.py'

	  swan_loc = swan.__file__
	  swan_loc = '/'.join(swan_loc.split('/')[:-1])
	  swan_loc += '/nx_pylab.py'
	  shutil.copyfile(swan_loc, nx_loc)
	except:
  		raise Exception('Unable to patch networkx to allow for dashed edges.')	

if __name__ == '__main__': main()