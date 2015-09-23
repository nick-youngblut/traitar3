
# traitAR - the microbial trait analyzer
traitAR is a software for characterizing microbial samples from nucleotide or protein sequences. It can accurately phenotype 67 diverse traits related to growth, shape, carbon source utilized etc.

# Installation
traitAR is available for Linux via the python packaging index. We strongly encourage to install it with pip rather than manually by cloning the repository. Install locally by

``pip install -i https://testpypi.python.org/pypi  --user traitar``

You need to adjust the PATH variable in the ~/.bashrc so that the bash finds the executables provided by packages installed with pip

``PATH=$PATH:~/.local/bin/``

You can also install globally with

``sudo pip install -i https://testpypi.python.org/pypi``

traitAR further needs Prodigal and HMMer available on the command line. For parallel execution it further requires GNU parallel.
All three are available as preconfigured package for the major Linux installation e.g. for Debian / Ubuntu install by

``sudo apt-get install parallel prodigal HMMer``

``traitar -h`` will provide help regarding the differnt options of traitAR.
traitAR requires the Pfam 27.0 HMM models. These are not distributed with this package but need to be downloaded

``traitar config <path to folder>``

will download and extract the Pfam models to the specified folder. The download can take a while depending on you internet connection. The Pfam models are extracted automatically to the same folder. This can slow down some machines. 
You may also download and extract the Pfam models manually from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz and run 

``traitar config --local <path to Pfam folder>``

to let traitAR know where.
# Basic usage

``traitar phenotype <in dir>  <sample file> from_nucleotides <out_dir> `` 

will trigger the standard workflow of traitAR, which is to predict open reading frames with Prodigal, annotate the coding sequences with Pfam families using HMMer and finally predict phenotypes from the models for the 67 traits. 

``traitar phenotype <in dir>  <sample file> from_genes <out_dir> `` 
 
assumes that gene prediction has been conducted already externally. In this case analysis will start with the Pfam annotation. If the output directory already exists, traitAR will offer to recompute or resume the individual analysis steps. This option is only available if the process is run interactively.
## Parallel usage
traitAR can benefit from parallel execution. The ``-c`` parameter sets the number of processes used e.g. ``-c 2`` for using two processes

``traitar phenotype <in dir>  <sample file> from_nucleotides out_dir -c 2`` 

This requires installing GNU parallel as noted above.
##Run traitAR with packaged sample data
``traitar phenotype <traitar_dir>/traitar/data/sample_data <traitar_dir>/traitar/data/sample_data/samples.txt from_nucelotides <out_dir> -c 2`` will trigger phenotyping of Mycoplasma pneumoniae 309 and Mycoplasma genitalium G37.Computation should be done within 5 minutes. 
## Output
traitAR provides the gene prediction results in ``<out_dir>/gene_prediction``, the Pfam annotation in ``<out_dir>/pfam_annotation`` and the phenotype prediction in``<out_dir>/phenotype prediction``. The phenotype prediction is summarized in heatmaps individually for each of the two algorithms and combined. Flat files and tab separated files provide access to the raw prediction results.

# Web service
traitAR is also offered as a web service at
http://algbio.cs.uni-duesseldorf.de/webapps/wa-webservice/pipe.php?pr=phenolyzer 
