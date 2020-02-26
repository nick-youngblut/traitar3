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
def remove(args):
    """remove phenotypes from phenotype tar archive"""
    modify.remove(args.archive_f, args.phenotypes, args.out_f, args.keep)

def ParseArgs(test_args=None, subparsers=None):
    desc = 'remove phenotypes from a given phenotype archive'
    epi = """DESCRIPTION:
    remove phenotypes from a given phenotype archive
    """
    if subparsers:
        parser = subparsers.add_parser('remove', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.add_argument("archive_f",
                        help = 'phenotype model archive file, which shall be modified')
    parser.add_argument("phenotypes",
                        help = 'phenotypes to be removed')
    parser.add_argument("out_f",
                        help = 'out file for the modified phenotype tar archive')
    parser.add_argument("--keep", action = 'store_true',
                        help = 'instead of remove the given phenotypes' + \
                        ' keep them and forget the rest')
    parser.set_defaults(func = remove)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
