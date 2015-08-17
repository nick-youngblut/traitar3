import urllib2
import json
import sys
import os
import gzip 

def download(args):
    """download Pfam HMMs and write download destination into config file"""
    attempts = 0
    if not args.local: 
        while attempts < 3:
            try:
                response = urllib2.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz", timeout = 5)
                content = response.read()
                f = open(os.path.join(args.download_dest, "Pfam-A.hmm.gz"), 'w' )
                f.write( content)
                f.close()
                with gzip.open(os.path.join(args.download_dest, "Pfam-A.hmm.gz"), 'rb') as zf:
                    with open(os.path.join(args.download_dest, "Pfam-A.hmm"), 'wb') as out_f:
                        out_f.write(zf.read())
                break
            except urllib2.URLError as e:
                attempts += 1
                print e
    with open(os.path.join("/".join(sys.argv[0].split("/")[:-2]), "config.json"), 'w') as config:
        config.write(json.dumps({"pfam_hmms": os.path.join(args.download_dest, "Pfam-A.hmm")}))
