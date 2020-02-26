#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import sys
import re
import json
import glob
import shutil
import logging
import tarfile
import argparse
import subprocess
import multiprocessing as mp
from functools import partial
from distutils.spawn import find_executable
from pkg_resources import resource_filename
## 3rd party
import pandas as ps
## application
from traitar._version import __version__
from traitar.Commands import phenotype as phenotype_cmd
from traitar.Commands import annotate as annotate_cmd
from traitar.Commands import pfam as pfam_cmd
from traitar.Commands import show as show_cmd
from traitar.Commands import evaluate as evaluate_cmd
from traitar.Commands import new as new_cmd
from traitar.Commands import remove as remove_cmd

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# function init
def main(args = None):
    #parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version', action = 'version', version = __version__)
    ## parent parsers
    ### general
    parent1_p = argparse.ArgumentParser(add_help=False)
    parent1_p.add_argument('-c', '--cpus', type = int, default = 1,
                          help='CPUs used for running hmmsearch & other executables' + \
                          ' (default: %(default)s)')
    parent1_p.add_argument('-x', '--parallel', type=int, default = 1,
                          help='Number of samples to process in parallel' + \
                          ' (default: %(default)s)')    
    parent1_p.add_argument('-o', '--overwrite', action='store_true', default=False,
                        help='Overwrite output directories (default: %(default)s)')
    ### annotate + phenotype
    parent2_p = argparse.ArgumentParser(add_help=False)
    parent2_p.add_argument('pfam_dir',
                        help='Download directory of pfam subcommand')
    parent2_p.add_argument("input_dir", help='directory with the input data')
    parent2_p.add_argument('sample2file',
                        help='Mapping from samples to fasta files' + \
                        ' (also see gitHub help):\n' + \
                        ' sample1_file_name{tab}sample1_name\n' + \
                        ' sample2_file_name{tab}sample2_name')
    parent2_p.add_argument('mode',
                        help='Either from_genes if gene prediction amino' + \
                        ' acid fasta is available in otherwise' + \
                        ' from_nucleotides in this case Prodigal is used to' + \
                        ' determine the ORFs from the nucleotide fasta files' + \
                        ' in input_dir',
                          choices=['from_genes', 'from_nucleotides', 'from_annotation_summary'])
    parent2_p.add_argument('output_dir',
                        help='Output directory (default: %(default)s)',
                        default='phenolyzer_output')
    parent2_p.add_argument('-g', '--gene_gff_type',
                          help='If the input is amino acids the type of gene' + \
                          ' prediction GFF file can be specified for mapping' + \
                          ' the phenotype predictions to the genes (default: %(default)s)',
                          default = None,
                          choices = ["genbank", "refseq", "img", "prodigal", "metagenemark"])
    # Defining subparsers
    subparsers = parser.add_subparsers()
    ## phenotype 
    phenotype_cmd.ParseArgs(subparsers=subparsers, parents=[parent1_p, parent2_p])
    ## annotate
    annotate_cmd.ParseArgs(subparsers=subparsers, parents=[parent1_p, parent2_p])
    ## pfam
    pfam_cmd.ParseArgs(subparsers=subparsers)
    ## show
    show_cmd.ParseArgs(subparsers=subparsers)
    ## new
    new_cmd.ParseArgs(subparsers=subparsers)
    ## remove
    remove_cmd.ParseArgs(subparsers=subparsers)
    ## evaluate
    evaluate_cmd.ParseArgs(subparsers=subparsers)

    # unit test args or command line?
    if args:
        args = parser.parse_args(args)
    else:
        args = parser.parse_args()

    # running subcommands
    if len(vars(args)) > 0:
        args.func(args)
    else:
        parser.parse_args(['--help'])
            
if __name__ == "__main__":
    main()
