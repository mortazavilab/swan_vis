
from setuptools import setup
setup(
  name = 'swan_vis',
  packages = ['swan_vis'],
  version = 'v2.0',
  license='MIT',
  description = 'swan is a tool for visualizing and analyzing transcript isoforms',
  author = 'Fairlie Reese',
  author_email = 'fairlie.reese@gmail.com',
  url = 'https://github.com/mortazavilab/swan_vis/',
  download_url='https://github.com/mortazavilab/swan_vis/archive/refs/tags/v2.0.tar.gz',
  keywords = ['swan', 'transcription', 'isoform', 'visualization'],
  python_requires='>=3.6',
  install_requires=[
    'fpdf',
    'matplotlib',
    'scipy',
    'statsmodels>=0.12.2',
    'numpy',
    'pandas>=1.1.0',
    'pytest',
    'networkx>=2.5',
    'tqdm',
    'anndata>=0.7',
    'diffxpy==0.7.4',
    ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research ',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
  ]
)
