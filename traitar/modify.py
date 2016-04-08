import pandas as ps
import tarfile
import os

mfs = ["%s_bias.txt", "%s_feats.txt","%s_majority_features+weights.txt"]

def validate(model_dir, pts):
    """validate that there is a model for each phenotype"""
    for i in pts.index:
        for j in mfs:
            #if not os.path.exists(os.path.join(model_dir, j % str(int(i)-1))):
            if not os.path.exists(os.path.join(model_dir, j % i)):
                raise Exception

def extend(archive_f, model_dir, pf2acc_desc_f, pts, pfam_v): 
    """extend an existing model archive by additional models"""
    #validate new phenotype models 
    validate(model_dir, pt2acc_desc_f)
    pt2acc = tarfile.open(archive_f, mode = "r:gz").extractfile("pt2acc_f.txt")
    #check if phenotype models not already exist
    pt2desc = ps.read_csv(pt2acc, sep = "\t", index_col = 0, header = None)
    #add phenotype models to the archive

def new(models_dir, pf2acc_f, pt2desc_f, hmm_model_f,  archive_name):
    """create new archive with phenotype models"""
    #read in pf and pt accessions 
    pts = ps.read_csv(pt2desc_f, sep = "\t", index_col = 0)
    validate(models_dir, pts)
    #create tar archive
    t = tarfile.open("%s.tar.gz" % archive_name, "w:gz")
    t.add(pf2acc_f, arcname = "pf2acc_desc.txt")
    t.add(pt2desc_f, arcname = "pt2acc.txt")
    #TODO include HMMER info
    for i in pts.index:
        for j in mfs:
            #t.add(os.path.join(models_dir, j % str(int(i)-1)), arcname = os.path.basename(os.path.join(models_dir, j % str(int(i)-1))))
            t.add(os.path.join(models_dir, j % i), arcname = os.path.basename(os.path.join(models_dir, j % i)))
    t.close()
         
    
    
    
    
