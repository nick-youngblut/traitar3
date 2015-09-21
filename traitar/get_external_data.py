import urllib2
import json
import sys
import os
import gzip 
import requests
import traitar

def download(args):
    """download Pfam HMMs and write download destination into config file"""
    attempts = 0
    if not args.local: 
        while attempts < 3:
            try:
                response = urllib2.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz", timeout = 5)
                with open(os.path.join(args.download_dest, "Pfam-A.hmm.gz"), 'w' ) as f:
                    CHUNK = 1000000
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        else:
                            f.write(chunk)
                with gzip.open(os.path.join(args.download_dest, "Pfam-A.hmm.gz"), 'rb') as zf:
                    with open(os.path.join(args.download_dest, "Pfam-A.hmm"), 'wb') as out_f:
                        out_f.write(zf.read())
                break
            except urllib2.URLError as e:
                attempts += 1
                print e
    with open(os.path.abspath(os.path.dirname(traitar.__file__), "config.json"), 'w') as config:
        config.write(json.dumps({"pfam_hmms": os.path.join(args.download_dest, "Pfam-A.hmm")}))
