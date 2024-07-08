![traitar3](https://github.com/nick-youngblut/traitar3/workflows/traitar3/badge.svg)

# Traitar3 &ndash; the microbial trait analyzer (for Python3)
Traitar3 is a Python3 implementation of [Traitar](https://github.com/aweimann/traitar),
which is a software for characterizing microbial samples from nucleotide or protein sequences.
Traitar(3) can accurately phenotype [67 diverse traits](traits.tsv).

> NOTE: This Python3 implementation is undergoing active development and could change at any time. Use it at your own risk!

### Table of Contents  
[Installation](#installation)  
[Basic usage](#basic-usage)  
[Results](#results)  
[Citing Traitar](#citing-traitar)  

<a name="installation"/>
<a name="basic-usage"/>
<a name="results"/>
<a name="docker"/>
<a name="citing-traitar"/>

# Installation
See the [.travis.yml](.travis.yml) file for installation instructions.

# Singularity

Traitar3 is also available as a singularity image which can be downloaded as follows.

```
# Pull singularity image file
singularity pull library://a_gihawi/traitar3/traitar3
# Test traitar3 using sample data
singularity exec traitar3_latest.sif traitar phenotype -c 1 --overwrite /db /traitar3/tests/data /traitar3/tests/data/samples.txt from_genes traitar_test_results
```

If all went well the results should be available in the ``traitar_test_results/``

Note: If your data directory is in a distant directory, they may need binding when launching the singularity. *i.e.* ``singularity exec --bind data/ traitar3_latest.sif traitar phenotype (options)``

# Basic usage

``traitar phenotype <in dir> <sample file> from_nucleotides <out_dir> `` 

will trigger the standard workflow of Traitar, which is to predict open reading frames with Prodigal, annotate the coding sequences provided as nucleotide FASTAs in the <in_dir> for all samples in <sample_file> with Pfam families using HMMer and finally predict phenotypes from the models for the 67 traits. 

![Alt text](/workflow.png?raw=true "Optional Title")

The sample file has one column for the sample file names and one for the names as specified by the user. You can also specify a grouping of the samples in the third column, which will be shown in the generated plots. The template looks like following - The header row is mandatory; please also take a look at the sample file for the packaged example data:

sample_file_name{tab}sample_name{tab}category  
  sample1_file_name{tab}sample1_name[{tabl}sample_category1]
  sample2_file_name{tab}sample2_name[{tabl}sample_category2]

``traitar phenotype <in dir> <sample file> from_genes <out_dir> `` 
 
assumes that gene prediction has been conducted already externally. In this case analysis will start with the Pfam annotation. If the output directory already exists, Traitar will offer to recompute or resume the individual analysis steps. This option is only available if the process is run interactively.

### Parallel usage
Traitar can benefit from parallel execution:

* The ``-c`` parameter sets the number of processes used for `hmmer` and other subprocesses.
* The ``-p`` parameter sets the number of samples (genomes) processed in parallel.

``traitar phenotype <in dir>  <sample file> from_nucleotides out_dir -c 2`` 

This requires installing GNU parallel as noted above.

### Inspect phenotype classification models
Traitar can be used to inspect the protein families in each phenotype model:

``traitar show 'Glucose fermenter'``

will show the majority features i.e. the Pfam families that contribute to the assignment of the trait Glucose fermenter with *phypat* classifier to some genome sequence. Via --predictor the user may specify the classifier (phypat, phypat+PGL). 


### Run Traitar with packaged sample data
``traitar phenotype <traitar_dir>/data/sample_data <traitar_dir>/data/sample_data/samples.txt from_genes <out_dir> -c 2`` will trigger phenotyping of *Listeria grayi DSM_20601* and *Listeria ivanovii WSLC3009*. Computation should be done within 5 minutes. You can find out ``<traitar_dir>`` by running

```
python
>>> import traitar
>>> traitar.__path__
```


# Results
Traitar provides the gene prediction results in ``<out_dir>/gene_prediction``, the Pfam annotation in ``<out_dir>/pfam_annotation`` and the phenotype prediction in``<out_dir>/phenotype prediction``.

### Heatmaps
The phenotype prediction is summarized in heatmaps individually for the phyletic pattern classifier in ``heatmap_phypat.png``, for the phylogeny-aware classifier in ``heatmap_phypat_ggl.png`` and for both classifiers combined in ```heatmap_comb.png``` and provide hierarchical clustering dendrograms for phenotypes and the samples.

![Alt text](tests/data/traitar_out/phenotype_prediction/heatmap_combined.png?raw=true "Phenotype heatmap")

### Phenotype prediction - Tables and flat files
These heatmaps are based on tab separated text files e.g. ``predictions_majority-vote_combined.txt``. A negative prediction is encoded as 0, a prediction made only by the pure phyletic classifier as 1, one made by the phylogeny-aware classifier by 2 and a prediction supported by both algorithms as 3. ``predictions_flat_majority-votes_combined.txt`` provides a flat version of this table with one prediction per row. The expert user might also want to access the individual results for each algorithm in the respective sub folders ``phypat`` and ``phypat+PGL``.

### Phenotype-relevant protein families and feature tracks
Traitar will link the protein families and predicted phenotypes. The results can be found in ``phypat/feat_gffs`` and ``phypat+PGL/feat_gffs`. If the user picked the 'from nucleotides' option, Traitar will also generate GFF files that link the genes called by Prodidgal with the important protein families. The phenotype-specific protein family annotations tracks can be visualized via GFF files in a genome browser of choice.

#### Feature tracks with *from_genes* option (experimental feature)
If the *from_genes* option is set, the user may specify gene GFF files via an additional column called gene_gff in the sample file. As gene ids are not consistent across gene GFFs from different sources e.g. img, RefSeq or Prodigal the user needs to specify the origin of the gene gff file via the -g / --gene_gff_type parameter. Still there is no guarantee that this works currently. Using samples_gene_gff.txt as the sample file in the above example will generate phenotype-specific Pfam tracks for the two genomes. 

``traitar phenotype . samples_gene_gff.txt from_genes traitar_out -g refseq``

# Citing Traitar

If you use Traitar in your research, please cite our paper:

**From genomes to phenotypes: Traitar, the microbial trait analyzer**  
Aaron Weimann, Kyra Mooren, Jeremy Frank, Phillip B Pope, Andreas Bremges, Alice C McHardy  
*mSystem* (2016) doi:[10.1101/043315](http://dx.doi.org/10.1128/mSystems.00101-16)
