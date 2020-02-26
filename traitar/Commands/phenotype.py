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
import argparse
## application
from traitar.traitar_cls import Traitar

# init funcs
def phenolyze(args):
    """annotate and then run phenotyping"""
    p = Traitar(args.input_dir, args.output_dir, args.sample2file,
                cpu = args.cpus, parallel = args.parallel,
                heatmap_out = args.rearrange_heatmap,
                heatmap_format = args.heatmap_format,
                no_heatmap_phenotype_clustering = args.no_heatmap_phenotype_clustering,
                no_heatmap_sample_clustering = args.no_heatmap_sample_clustering,
                gene_gff_type = args.gene_gff_type,
                primary_models = args.primary_models,
                secondary_models = args.secondary_models,
                pfam_dir = args.pfam_dir,
                overwrite = args.overwrite)
    
    #check if user wants to supply pre-computed annotation summary
    if not args.mode == "from_annotation_summary":
        logging.info('Running annotate as part of predict')
        p.annotate(args.mode)
    p.phenolyze(args.mode)

def ParseArgs(test_args=None, subparsers=None, parents=None):
    desc = 'Phenotype'
    epi = """DESCRIPTION:
    Annotate genomes and then run phenotyping
    """
    if subparsers:
        parser = subparsers.add_parser('phenotype', description=desc,
                                       epilog=epi, parents = parents, 
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## phenotype
    parser.add_argument("-p", "--primary_models",
                        help='Primary phenotype models archive')
    parser.add_argument("-s", "--secondary_models",
                        help='Secondary phenotype models archive')
    parser.add_argument("-r", "--rearrange_heatmap",
                        help='Recompute the phenotype heatmaps based on a' + \
                        ' subset of previously annotated and phenotyped samples' + \
                        '(default: %(default)s)', default = None)
    parser.add_argument("--no_heatmap_sample_clustering", action = 'store_true',
                        help = 'if option is set, don\'t cluster the phenotype' + \
                        ' heatmaps by samples and keep input ordering')
    parser.add_argument("--no_heatmap_phenotype_clustering", action = 'store_true',
                        help = "if option is set, don't cluster the heatmaps by phenotype and keep input ordering")
    parser.add_argument("-f", "--heatmap_format", choices = ["png", "pdf", "svg", "jpg"],
                        default='pdf', help = "choose file format for the heatmap")    
    parser.set_defaults(func = phenolyze)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
