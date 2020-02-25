#!/usr/bin/env python
from __future__ import print_function
"""script to create a summary matrix and a gene2hmm mapping from filtered and aggregated hmmer output files"""
# import
import os
import sys
import logging
## 3rd party
import pandas as ps
## application
from traitar.PhenotypeCollection import PhenotypeCollection

# functions
def gene2hmm(domtblout_list, pt_models, gene2hmm_out = None, is_gene2hmm = False):
    """ Function to create a summary matrix and a gene2hmm mapping from
    filtered and aggregated hmmer output files """
    #read accession file 
    accs = pt_models.get_pf2desc()
    sum_df = ps.DataFrame(ps.np.zeros((len(domtblout_list), accs.shape[0])))
    #set index to the actual files by getting rid of the preceding path
    sum_df.index = [i.split("/")[-1].replace("_filtered_best.dat", "") for i in domtblout_list]
    sum_df.columns = accs.index 
    #create a gene2hmm dictionary
    gene2hmm = {}  
    for f in domtblout_list:
        gene_df = ps.read_csv(f, sep = "\t")
        if not f in gene2hmm:
            f_mod = f.split("/")[-1].replace("_filtered_best.dat", "")
            gene2hmm[f_mod] = {}
        for i in gene_df.index:
            #append genes to gene2hmm dict
            gene = gene_df.loc[i, 'target name']
            if pt_models.get_hmm_name() == "pfam":
                query = gene_df.iloc[i, 4].split('.')[0]
            if pt_models.get_hmm_name() == "dbcan":
                query = gene_df.iloc[i, 3].split('.')[0]
            if is_gene2hmm:
                if gene not in gene2hmm[f_mod]:
                    gene2hmm[f_mod][gene] = [query] 
                else:
                    gene2hmm[f_mod][gene].append(query)
            #add annotations to summary file
            try:
                sum_df.loc[f_mod, query] += 1
            except KeyError:
                print(query)
    #write annotation summary to disk
    if not is_gene2hmm:
        sum_df.to_csv(gene2hmm_out, sep = "\t")
    return sum_df, gene2hmm

def main(archive_f, in_filtered_best_fs, outfile):
    """ Main interface 
    Params:
      archive_f : str, archive (dat) file
      in_filtered_best_fs : list, filtered hmm file paths
      outfile : str, output file path
    """    
    outdir = os.path.split(outfile)[0]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    pt_models = PhenotypeCollection(archive_f)
    sum_df, gene2hmms = gene2hmm(in_filtered_best_fs,
                                 pt_models,
                                 gene2hmm_out = outfile)
    logging.info('File written: {}'.format(outfile))

def best_fs_list(infile):
    x = []
    with open(infile) as inF:
        for line in inF:
            x.append(line.rstrip())
    return x
    
if __name__ == "__main__":
    import argparse
    msg = 'Generate a summary matrix from the filtered best hmmer annotation files'
    parser = argparse.ArgumentParser(msg)
    parser.add_argument("outfile",
                        help='summary matrix output file')
    parser.add_argument("in_filtered_best_fs",
                        help='file with filtered.best file name per row')
    parser.add_argument("archive_f",
                        help = "phenotype archive file")
    args = parser.parse_args()

    args.in_filtered_best_fs = best_fs_list(args.in_filtered_best_fs)
    main(args.archive_f, args.in_filtered_best_fs, args.outfile)
