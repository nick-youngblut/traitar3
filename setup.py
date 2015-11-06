from setuptools import setup

setup(name='traitar',
        version='0.1.6',
        description='traitAR - The microbial trait analyzer',
        url = 'http://github.com/aweimann/traitAR',
        author='Aaron Weimann',
        author_email='weimann@hhu.de',
        license='GNU General Public License, version 3 (GPL-3.0)',
        packages=['traitar'],
        scripts = ['traitar/traitar', 'traitar/merge_preds.py', 'traitar/heatmap.py', 'traitar/domtblout2gene_generic.py', 'traitar/predict.py', 'traitar/hmmer2filtered_best.py'],
        zip_safe=False,
        install_requires = ["pandas >= 0.14.1", "matplotlib >= 1.0", "scipy >= 0.9.0", "numpy >= 1.10.1"])
