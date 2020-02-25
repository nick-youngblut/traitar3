import os
from setuptools import find_packages
from setuptools import setup
import re


VERSIONFILE=os.path.join('traitar', '_version.py')
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else: 
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

long_description = open('README.rst', 'r').read()

data_files = [os.path.join('data', 'colors.txt'),
              os.path.join('data', 'models', 'phypat.tar.gz'),           
              os.path.join('data', 'models', 'phypat+PGL.tar.gz'),
              os.path.join('data', 'sample_data', '1457190.3.RefSeq.faa'),
              os.path.join('data', 'sample_data', '1457190.3.RefSeq.gff'),
              os.path.join('data', 'sample_data', '525367.9.RefSeq.faa'),
              os.path.join('data', 'sample_data', '525367.9.RefSeq.gff'),
              os.path.join('data', 'sample_data', 'samples_gene_gff.txt'),
              os.path.join('data', 'sample_data', 'samples.txt')]              

setup(name='traitar',
      version = verstr,
      description='traitar - The microbial trait analyzer',
      long_description = long_description,
      url = 'http://github.com/aweimann/traitar',
      author='Aaron Weimann',
      author_email='weimann@hhu.de',
      license='GNU General Public License, version 3 (GPL-3.0)',
      packages= ['traitar'],
#      include_package_data = True,
      package_data={'traitar': data_files},           
      scripts = ['traitar/traitar',
                 'traitar/merge_preds.py',
                 'traitar/heatmap.py',
                 'traitar/domtblout2gene_generic.py',
                 'traitar/predict.py',
                 'traitar/hmmer2filtered_best.py',
                 'traitar/hmm2gff.py'],
      zip_safe=False,
      install_requires = ["pandas >= 0.13.1",
                          "matplotlib >= 1.3.1",
                          "numpy >= 1.6",
                          "scipy >= 0.13.3"]
)
