
from setuptools import setup
setup(
  name = 'swan_vis',
  packages = ['swan_vis'],
  version = 'v3.2',
  license='MIT',
  description = 'swan is a tool for visualizing and analyzing transcript isoforms',
  author = 'Fairlie Reese',
  author_email = 'fairlie.reese@gmail.com',
  url = 'https://github.com/mortazavilab/swan_vis/',
  download_url='https://github.com/mortazavilab/swan_vis/archive/refs/tags/v3.2.tar.gz',
  keywords = ['swan', 'transcription', 'isoform', 'visualization'],
  python_requires='>=3.8',
  install_requires=[
    'fpdf',
    'matplotlib',
    'scipy',
    'statsmodels>=0.12.2',
    'numpy',
    'pandas<2.0',
    'pytest',
    'networkx>=2.5',
    'tqdm',
    'anndata>=0.8',
    'scanpy>=1.9'
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
