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
from pkg_resources import resource_filename
## application
from traitar.PhenotypeCollection import PhenotypeCollection

primary_default_models = resource_filename('traitar',
                                           os.path.join('data', 'models', 'phypat.tar.gz'))
secondary_default_models = resource_filename('traitar',
                                           os.path.join('data', 'models', 'phypat+PGL.tar.gz'))

# init funcs
def show(args):
    """show features for the given phenotype"""
    if args.models_f is not None:
        pc = [PhenotypeCollection(args.models_f)]
    else:
        if args.predictor == "phypat":
            pc = [PhenotypeCollection(primary_default_models)]
        else:
            pc = [PhenotypeCollection(secondary_default_models)]
    for i in pc:
        x = i.get_selected_features(args.phenotype, args.strategy,
                                    args.include_negative)
        x = x.to_csv(sys.stdout, sep = "\t", float_format='%.3f')

def ParseArgs(test_args=None, subparsers=None):
    desc = 'show features important for classification'
    epi = """DESCRIPTION:
    show features important for classification
    """
    if subparsers:
        parser = subparsers.add_parser('show', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)
    ## args
    parser.add_argument("phenotype",
                        help = "phenotype under investigation")
    parser.add_argument("--predictor",
                        help = "pick phypat or phypat+PGL classifier",
                        choices = ["phypat", "phypat+PGL"], default = "phypat")
    parser.add_argument("--strategy", choices = ["non-zero", "majority"],
                        default = "majority")
    parser.add_argument("-i", "--include_negative", action = 'store_true')
    parser.add_argument("-p", "--models_f",
                        help='phenotype models archive; if set' + \
                        ' look for the target in the phenotype in the' + \
                        'given phenotype collection')
    parser.set_defaults(func = show)

    # test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser
