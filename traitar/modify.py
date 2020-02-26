#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import sys
import tarfile
from io import StringIO     
## 3rd party
import pandas as pd
## application
import traitar.PhenotypeCollection

# vars
mfs = ["{}_bias.txt", "{}_feats.txt", "{}_non-zero+weights.txt"]

def validate(model_dir, pts):
    """validate that there is a model for each phenotype"""
    for i in pts.index:
        for j in mfs:
            mfs_file = os.path.join(model_dir, j.format(i))
            if not os.path.exists(mfs_file):
                sys.stderr.write('{} does not exist'.format(mfs_file))
                raise Exception

def remove(archive_f, phenotypes, out_f, keep = False):
    """remove given phenotypes from the pt archive and write a new archive"""
    pts = pd.read_csv(phenotypes, index_col = 0, header = None)
    pts.index = pts.index.values.astype('string')
    ptc = PhenotypeCollection.PhenotypeCollection(archive_f)
    pt2acc = ptc.get_pt2acc()
    pt2acc.index = pt2acc.index.values.astype('string')
    if not keep:
        pt2acc.drop(pts.index, inplace = True)
    else:
        pt2acc = pt2acc.loc[pts.index,]
    pt2acc_s = StringIO.StringIO(pt2acc.to_csv(sep = "\t", encoding = 'utf8'))
    pt2acc_tarinfo = tarfile.TarInfo("pt2acc.txt")
    pt2acc_tarinfo.size = len(pt2acc_s.buf)
    t = ptc.tar
    members = [] 
    members.append((t.extractfile("pf2acc_desc.txt"),
                    t.getmember("pf2acc_desc.txt")))
    members.append((t.extractfile("config.txt"),
                    t.getmember("config.txt")))
    for i in pt2acc.index:
        members.append((t.extractfile("{}_non-zero+weights.txt".format(i)),
                        t.getmember("{}_non-zero+weights.txt".format(i))))
        members.append((t.extractfile("{}_bias.txt".format(i)),
                        t.getmember("{}_bias.txt".format(i))))
        members.append((t.extractfile("{}_feats.txt".format(i)),
                        t.getmember("{}_feats.txt".format(i))))
    out_tar = tarfile.open("{}".format(out_f), "w:gz")
    out_tar.addfile(pt2acc_tarinfo, pt2acc_s)
    for m in members:
        out_tar.addfile(m[1], m[0])
    out_tar.close()

def extend(archive_f, model_dir, pf2acc_desc_f, pts, pfam_v): 
    """extend an existing model archive by additional models"""
    #validate new phenotype models 
    validate(model_dir, pt2acc_desc_f)
    ptc = PhenotypeCollection.PhenotypeCollection(archive_f)
    pt2acc = ptc.get_pt2acc()

def new(models_dir, pf2acc_f, pt2desc_f, hmm_name, hmm_model_f,  archive_name):
    """create new archive with phenotype models"""
    #read in pf and pt accessions 
    pts = pd.read_csv(pt2desc_f, sep = "\t", index_col = 0)
    validate(models_dir, pts)
    create_tar(models_dir, pf2acc_f, pt2desc_f, hmm_name, hmm_model_f,  archive_name)

def create_tar(models_dir, pf2acc_f, pt2desc_f, hmm_name, hmm_model_f,  archive_name):
    #create tar archive
    pt2desc = pd.read_csv(pt2desc_f, sep = "\t", index_col = 0) 
    t = tarfile.open("{}.tar.gz".format(archive_name), "w:gz")
    t.add(pf2acc_f, arcname = "pf2acc_desc.txt")
    t.add(pt2desc_f, arcname = "pt2acc.txt")
    config = [archive_name, hmm_name, hmm_model_f]
    config_df = pd.DataFrame(config,
                             index = ["archive_name", "hmm_name", "hmm_f"],
                             columns = ["value"])
    config_s = StringIO.StringIO(config_df.to_csv(sep = "\t"))
    config_tarinfo = tarfile.TarInfo("config.txt")
    config_tarinfo.size = len(config_s.buf)
    t.addfile(config_tarinfo, config_s)
    for i in pt2desc.index:
        for j in mfs:
            t.add(os.path.join(models_dir, j.format(i)),
                  arcname = os.path.basename(os.path.join(models_dir, j.format(i))))
    t.close()
         
    
    
    
    
