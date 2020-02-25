#!/usr/bin/env python
from __future__ import print_function
# import
## batteris
import os
import sys
import gzip 
import json
import logging
import urllib2
## application
import traitar

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# function init
def download(args):
    """download Pfam HMMs and write download destination into config file"""
    url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
    attempts = 0
    if not args.local: 
        while attempts < 3:
            try:
                if not os.path.exists(args.download_dir):
                    os.makedirs(args.download_dir)
                logging.info('Download dir: {}'.format(args.download_dir))         
                # downloading
                logging.info('Downloading: {}'.format(url))
                response = urllib2.urlopen(url, timeout = 5)
                with open(os.path.join(args.download_dir, "Pfam-A.hmm.gz"), 'w' ) as f:
                    CHUNK = 1000000
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        else:
                            f.write(chunk)
                # uncompress
                infile = os.path.join(args.download_dir, "Pfam-A.hmm.gz")
                outfile = os.path.join(args.download_dir, "Pfam-A.hmm")
                logging.info('Uncompressing {} to {}'.format(infile, outfile))
                with gzip.open(infile, 'rb') as zf:
                    with open(outfile, 'wb') as out_f:
                        for l in zf:
                            out_f.write(l)
                break
            except urllib2.URLError as e:
                attempts += 1
                print(e)

    # checking output
    if not os.path.exists(os.path.join(args.download_dir, 'Pfam-A.hmm')):
        msg = 'ERROR: no Pfam-A.hmms file in {}'
        sys.exit(msg.format(args.download_dir))

