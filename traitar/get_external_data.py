#!/usr/bin/env python
from __future__ import print_function
# import
## batteris
import os
import sys
import gzip 
import json
import logging
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, URLError
import traitar
from traitar.utils import retry

# function init
@retry(URLError, tries=5, delay=5, backoff=2)
def urlopen_with_retry(url, timeout=10):
    return urlopen(url, timeout)

def download(args):
    """download Pfam HMMs and write download destination into config file"""
    url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
    attempts = 0
    timeout = 20
    if not args.local: 
        while attempts < 3:
            try:
                if not os.path.exists(args.download_dir):
                    os.makedirs(args.download_dir)
                logging.info('Download directory: {}'.format(args.download_dir))         
                # downloading
                logging.info('Downloading: {}'.format(url))
                response = urlopen_with_retry(url, timeout = timeout)
                with open(os.path.join(args.download_dir, "Pfam-A.hmm.gz"), 'wb') as outF:
                    CHUNK = 1000000
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        else:
                            outF.write(chunk)
                # uncompress
                infile = os.path.join(args.download_dir, "Pfam-A.hmm.gz")
                outfile = os.path.join(args.download_dir, "Pfam-A.hmm")
                logging.info('Uncompressing {} to {}'.format(infile, outfile))
                with gzip.open(infile, 'rb') as inF:
                    with open(outfile, 'w') as outF:
                        for line in inF:
                            outF.write(line.decode('utf8'))
                break
            except URLError as e:
                attempts += 1
                timeout *= 2
                print(e)

    # checking output
    if not os.path.exists(os.path.join(args.download_dir, 'Pfam-A.hmm')):
        msg = 'ERROR: no Pfam-A.hmms file in {}'
        sys.exit(msg.format(args.download_dir))

