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

long_description = """Traitar3 - The microbial trait analyzer for Python3

See https://github.com/nick-youngblut/traitar3 for more information.
"""

data_files = [os.path.join('data', 'colors.txt'),
              os.path.join('data', 'models', 'phypat.tar.gz'),           
              os.path.join('data', 'models', 'phypat+PGL.tar.gz'),
              os.path.join('data', 'sample_data', '1457190.3.RefSeq.faa'),
              os.path.join('data', 'sample_data', '1457190.3.RefSeq.gff'),
              os.path.join('data', 'sample_data', '525367.9.RefSeq.faa'),
              os.path.join('data', 'sample_data', '525367.9.RefSeq.gff'),
              os.path.join('data', 'sample_data', 'samples_gene_gff.txt'),
              os.path.join('data', 'sample_data', 'samples.txt')]              

setup(name='traitar3',
      version = verstr,
      description='Traitar3 - The microbial trait analyzer for Python3',
      long_description = long_description,
      url = 'https://github.com/nick-youngblut/traitar3',
      author='Nick Youngblut',
      author_email='nyoungb2@gmail.com',
      license='GNU General Public License, version 3 (GPL-3.0)',
      packages = find_packages(),
      package_dir={'traitar' : 'traitar'},       
      package_data={'traitar': data_files},
      entry_points={
        'console_scripts': [
            'traitar = traitar.__main__:main'
        ]
      },
      scripts = ['traitar/merge_preds.py',
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
