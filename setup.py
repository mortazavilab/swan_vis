
from setuptools import setup
setup(
  name = 'swan_vis',
  packages = ['swan_vis'],
  version = '0.0.2',
  license='MIT',  description = 'swan is a tool for visualizing and analyzing transcript isoforms',
  author = 'Fairlie Reese',
  author_email = 'fairlie.reese@gmail.com',
  url = 'https://github.com/fairliereese/swan_vis/',
  download_url = 'https://github.com/fairliereese/swan_vis/archive/v0.0.2-alpha.tar.gz',
  keywords = ['swan', 'transcription', 'isoform'],
  install_requires=[
          'networkx',
          'numpy',
          'pandas',
          'matplotlib',
          'fpdf'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research ',    
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License', 
    'Programming Language :: Python :: 3.7',
  ],
)