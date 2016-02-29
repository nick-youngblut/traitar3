# traitar - the microbial trait analyzer
traitar is a software for characterizing microbial samples from nucleotide or protein sequences. It can accurately phenotype [67 diverse traits](traits.tsv).

### Table of Contents  
[Installation](#installation)  
[Basic usage](#basic-usage)  
[Results](#results)  
[Web service](#web-service)  


# Installation
traitar is available for Linux via the python packaging index. We tested the software under Ubuntu 14.04 LTS (trusty). We didn't test older versions but these might work as well. 
Prior to installation with pip make sure the following packages are installed by running

``sudo apt-get install python-scipy python-matplotlib python-pip python-pandas``
and optionally
``sudo apt-get install python-virtualenv``

These package might not be available for your Linux distribution. Please let us know, if you have troubles installing this software, so that we can assist you with the installation.

Once finished we strongly encourage you to install traitar with pip rather than by manually cloning the repository. Install locally by

``pip install traitar-<version>.tar.gz --user``

You need to add the following line to your ~/.bashrc to adjust the PATH variable so that the bash finds the executables needed for running traitar. 

``PATH=$PATH:~/.local/bin/``

You may need to run ``source ~/.bashrc`` once in your current session.

You can also install globally with

``sudo pip install traitar-<version>.tar.gz``
### Creating a virtual environment
You may want to use virtualenv to create a clean environment for traitar i.e. run

```
virtualenv <environment_name> --system-site-packages
source <environment_path>/bin/activate

pip install -U traitar-<version>.tar.gz
PATH=$PATH:<environment_path>/bin/
source ~/.bashrc
```

### Additional requirements
traitar further needs prodigal and hmmsearch available on the command line. For parallel execution it further requires GNU parallel.
All three are available as preconfigured package for many Linux installation e.g. for Debian / Ubuntu. For Ubuntu 14.04 it is available via the trusty-backports.

```
#include trusty backports source in the sources.list
"sudo deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse ">> /etc/apt/sources.list

#update sources
apt-get install update

#install packages
sudo apt-get install parallel prodigal hmmer
```

``traitar -h`` will provide help regarding the different options of traitar.
traitar requires the Pfam 27.0 HMM models. These are not distributed with this package but need to be downloaded

``traitar pfam <path to download folder>``

will download and extract the Pfam models to the specified folder. The download can take a while depending on you internet connection. The Pfam models are extracted automatically to the same folder. This can slow down some machines. 
You may also download and extract the Pfam models manually from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz and run 

``traitar pfam --local <path to Pfam folder>``

to let traitar know where.


# Basic usage

``traitar phenotype <in dir>  <sample file> from_nucleotides <out_dir> `` 

will trigger the standard workflow of traitar, which is to predict open reading frames with Prodigal, annotate the coding sequences provided as nucleotide FASTAs in the <in_dir> for all samples in <sample_file> with Pfam families using HMMer and finally predict phenotypes from the models for the 67 traits. 

The sample file has one column for the sample file names and one for the names as specified by the user. You can also specify a grouping of the samples in the third column, which will be shown in the generated plots. The template looks like following - Please also take a look at the sample file for the packaged example data:

sample1_file_name{tab}sample1_name[{tabl}sample_category1]
sample2_file_name{tab}sample2_name[{tabl}sample_category2]

``traitar phenotype <in dir>  <sample file> from_genes <out_dir> `` 
 
assumes that gene prediction has been conducted already externally. In this case analysis will start with the Pfam annotation. If the output directory already exists, traitar will offer to recompute or resume the individual analysis steps. This option is only available if the process is run interactively.

### Parallel usage
traitar can benefit from parallel execution. The ``-c`` parameter sets the number of processes used e.g. ``-c 2`` for using two processes

``traitar phenotype <in dir>  <sample file> from_nucleotides out_dir -c 2`` 

This requires installing GNU parallel as noted above.

### Run traitar with packaged sample data
``traitar phenotype <traitar_dir>/data/sample_data <traitar_dir>/data/sample_data/samples.txt from_genes <out_dir> -c 2`` will trigger phenotyping of *Listeria grayi DSM_20601* and *Listeria ivanovii WSLC3009*. Computation should be done within 5 minutes. You can find out ``<traitar_dir>`` by running

```
python
>>> import traitar
>>> traitar.__path__
```


# Results
traitar provides the gene prediction results in ``<out_dir>/gene_prediction``, the Pfam annotation in ``<out_dir>/pfam_annotation`` and the phenotype prediction in``<out_dir>/phenotype prediction``.

### Heatmaps
The phenotype prediction is summarized in heatmaps individually for the phyletic pattern classifier in ``heatmap_phypat.png``, for the phylogeny-aware classifier in ``heatmap_phypat_ggl.png`` and for both classifiers combined in ```heatmap_comb.png``` and provide hierarchical clustering dendrograms for phenotypes and the samples.

### Phenotype prediction - Tables and flat files
These heatmaps are based on tab separated text files e.g. ``predictions_majority-votes_combined.txt``. A negative prediction is encoded as 0, a prediction made only by the pure phyletic classifier as 1, one made by the phylogeny-aware classifier by 2 and a prediction supported by both algorithms as 3. ``predictions_flat_majority-votes_combined.txt`` provides a flat version of this table with one prediction per row. The expert user might also want to access the individual results for each algorithm in the respective sub folders ``phypat`` and ``phypat+PGL``.

### Feature tracks
If traitar is run from_nucleotides it will generate a link between the Prodigal gene prediction and predicted phenotypes in ``phypat/feat_gffs`` and ``phypat+PGL/feat_gffs`` (no example in the sample data). The user can visualize gene prediction  phenotype-specific Pfam annotations tracks via GFF files.

#### Feature tracks with *from_genes* option (experimental feature)
If the *from_genes* option is set, the user may specify gene GFF files via an additional column called gene_gff in the sample file. As gene ids are not consistent across gene GFFs from different sources e.g. img, RefSeq or Prodigal the user needs to specify the origin of the gene gff file via the -g / --gene_gff_type parameter. Still there is no guarantee that this works currently. Using samples_gene_gff.txt as the sample file in the above example will generate phenotype-specific Pfam tracks for the two genomes. 

``traitar phenotype . samples_gene_gff.txt from_genes traitar_out -g refseq``


# Web service
We also offer traitar as a web service at
http://algbio.cs.uni-duesseldorf.de/webapps/wa-webservice/pipe.php?pr=phenolyzer 
