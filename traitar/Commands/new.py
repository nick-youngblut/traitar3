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
from traitar import modify

# init funcs
def new(args):
    """create a new phenotype model archive"""
    modify.new(args.models_dir, args.pf_acc2desc, args.pt_id2acc_desc,
               args.hmm_name, args.hmm_model_f, args.archive_name)

def ParseArgs(test_args=None, subparsers=None):
    desc = 'create new phenotype model archive'
    epi = """DESCRIPTION:
    create new phenotype model archive
    """
    if subparsers:
        parser = subparsers.add_parser('new', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.add_argument("models_dir",
                        help='directory with phenotype models to be included')
    parser.add_argument("pf_acc2desc",
                        help='a mapping between Pfam families and phenotype ids to accessions')
    parser.add_argument("pt_id2acc_desc",
                        help='a mapping between phenotype ids and descriptions')
    parser.add_argument("hmm_name", help='hmm database used',
                        choices = ["dbcan", "pfam"])
    parser.add_argument("hmm_model_f",
                        help='hmm database compatible with the phenotype archive')
    parser.add_argument("archive_name",
                        help='name of the model, which is created')
    parser.set_defaults(func = new)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
