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
from traitar.PhenotypeCollection import PhenotypeCollection


# init funcs
def evaluate(args):
    evaluation.evaluate.evaluate(args.out, args.gold_standard_f,
                                 args.traitar_pred_f, args.min_samples,
                                 args.are_pt_ids, args.phenotype_archive)

def ParseArgs(test_args=None, subparsers=None):
    desc = 'compare Traitar predictions against a given standard of truth'
    epi = """DESCRIPTION:
    compare Traitar predictions against a given standard of truth
    """
    if subparsers:
        parser = subparsers.add_parser('evaluate', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.add_argument("traitar_pred_f",
                        help = "phenotype prediction matrix as return by Traitar")
    parser.add_argument("gold_standard_f",
                        help = "phenotype matrix with standard of truth")
    parser.add_argument("--are_pt_ids", action = 'store_true',
                        help = 'Set if the gold standard phenotype are indexed' + \
                        ' via phenotype IDs rather than accessions')
    parser.add_argument("--phenotype_archive",
                        help = "Needed if gold standard uses an accession index for mapping")
    parser.add_argument("--min_samples", "-m", default = 5,
                        help='minimum number of positive and negative samples to' + \
                        ' consider phenotypes for calculation of the macro accuracy')
    parser.add_argument("out",
                        help = "output directory")
    parser.set_defaults(func = evaluate)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
