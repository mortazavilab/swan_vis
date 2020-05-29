
from setuptools import setup
setup(
  name = 'swan_vis',
  packages = ['swan_vis'],
  version = '0.0.4',
  license='MIT',  description = 'swan is a tool for visualizing and analyzing transcript isoforms',
  author = 'Fairlie Reese',
  author_email = 'fairlie.reese@gmail.com',
  url = 'https://github.com/fairliereese/swan_vis/',
  download_url = 'https://github.com/fairliereese/swan_vis/archive/v0.0.4-alpha.tar.gz',
  keywords = ['swan', 'transcription', 'isoform', 'visualization'],
  install_requires=[
          'networkx',
          'numpy',
          'pandas',
          'matplotlib',
          'fpdf',
          'seaborn',
          'diffxpy',
          'tensorflow',
          'tensorflow-probability'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research ',    
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License', 
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
  ],
)

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