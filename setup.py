from setuptools import setup

setup(name='traitar',
        version='0.1.0',
        description='traitAR - The microbial trait analyzer',
        url = 'http://github.com/aweimann/traitAR',
        author='Aaron Weimann',
        author_email='weimann@hhu.de',
        license='GNU General Public License, version 3 (GPL-3.0)',
        packages=['traitar'],
        scripts = ['traitar/phenolyzer.py', 'traitar/merge_preds.py', 'traitar/heatmap.py', 'traitar/domtblout2gene_generic.py', 'traitar/predict.py', 'traitar/hmmer2filtered_best.py'],
        package_data = {"traitar": ["data/sorted_accessions.txt", "data/models/phypat+GGL.tar.gz", "data/models/phypat.tar.gz"]},
        zip_safe=False,
        install_requires = ["pandas >= 0.14.1", "matplotlib >= 1.0", "scipy >= 0.9.0"])
