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
def annotate(args):
    """annotate input"""
    p = Traitar(args.input_dir, args.output_dir, args.sample2file,
                cpu = args.cpus, parallel = args.parallel,
                gene_gff_type = args.gene_gff_type,
                overwrite = args.overwrite)
    p.annotate(args.mode)

def ParseArgs(test_args=None, subparsers=None, parents=None):
    desc = 'Annotate'
    epi = """DESCRIPTION:
    Annotate genomes
    """
    if subparsers:
        parser = subparsers.add_parser('annotate', description=desc,
                                       epilog=epi, parents=parents,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.set_defaults(func = annotate)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
