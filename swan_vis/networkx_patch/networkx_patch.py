try:
  import networkx as nx 
  import shutil
  nx_loc = nx.__file__
  nx_loc = '/'.join(nx_loc.split('/')[:-1])
  nx_loc += '/drawing/nx_pylab.py'
  nx_file = 'swan_vis/networkx_patch/nx_pylab.py'
  shutil.copyfile(nx_file, nx_loc)
except:
  raise Exception('Unable to patch networkx to allow for dashed edges.')