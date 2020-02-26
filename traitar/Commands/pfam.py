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
from traitar import get_external_data

# init funcs
def ParseArgs(test_args=None, subparsers=None):
    desc = 'Download pfam files'
    epi = """DESCRIPTION:
    Download and uncompress pfam files.
    The files are required for gene annotation.
    """
    if subparsers:
        parser = subparsers.add_parser('pfam', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.add_argument("download_dir",
                        help = 'Download Pfam HMMs into this directory' + \
                        ' and untar and unzip them')
    parser.add_argument("--local", "-l", action = 'store_true',
                        help = 'The Pfam HMMs are in the above directory with name "Pfam-A.hmm"')
    parser.set_defaults(func = get_external_data.download)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
