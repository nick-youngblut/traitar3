#!/usr/bin/env python
from __future__ import print_function
"""parse hmmer output file and generated filtered best file"""
# import
## batteries
import os
import sys
import logging
import argparse
from io import StringIO     
## 3rd party
import pandas as ps

# vars 
hmmer_colnames = ['target name','target accession','tlen','query name','accession',
                  'qlen','E-value','score_overall','bias_overall','#','of','c-Evalue',
                  'i-Evalue','score_domain','bias_domain','from_hmm','to_hmm','ali_from',
                  'ali_to','env_from','env_to','acc','description of target']

def filter_dbcan(m):
    return ((m.iloc[:,13] >= 25) & ((m.iloc[:,12] <= 0.001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] <= 80) | (m.iloc[:,12] <= 0.00001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] > 80))) 

def filter_pfam(m):
    #TODO check thresholds
    return (m.iloc[:,12] <= 0.01) & (m.iloc[:, 13] >= 25)

def apply_thresholds(infile_f, hmm_name, out_filt_f, out_excl_f):
    """parse HMMER output file and apply thresholds for e-valu and bit score)"""
    #preparse lines by replacing white space delimitation by tabs
    #skip header
    infile_f.readline()
    infile_f.readline()
    infile_f.readline()
    #replace whitespace by tabs and skip lines which start with a # char
    cleaned = "".join(filter(lambda x: not x.startswith("#"),
                             ["\t".join(i.split(None, 22)) for i in infile_f.readlines()]))
    #read  tab delimited hmmer output file with pandas via stringIO
    #account for cases where hmmer didn't return any hits
    try:
        m = ps.read_csv(StringIO(cleaned), sep = "\t",  header = None)
    except ValueError:
        m = ps.DataFrame(columns = hmmer_colnames)
    m.columns = hmmer_colnames
    if hmm_name == "pfam":
        keep = filter_pfam(m)
    if hmm_name == "dbcan":
        keep = filter_dbcan(m)
    #apply eval threshold
    m_eval = m.loc[keep, :]
    m_eval = m.loc[keep, :]
    if not out_filt_f is None: 
        m_eval.to_csv(out_filt_f, sep = "\t")
    m_eval_excl = m.loc[~keep, :]    
    if not out_excl_f is None:
        m_eval_excl.to_csv(out_excl_f, sep = "\t")
    return m_eval
    
    
def aggregate_domain_hits(filtered_df, out_f):
    #sort by gene identifier and Pfam
    with open(out_f, 'w') as out_fo:
        ps.DataFrame(filtered_df.columns).T.to_csv(out_f, sep = "\t",
                                                   index = False,
                                                   header = False,
                                                   mode = 'a')
        filtered_df.sort_values(by = ["target name", "query name"], inplace = True)
        if filtered_df.shape[0] > 0:
            current_max = filtered_df.iloc[0,] 
        else:
            #nothing todo
            return
        for i in range(1, filtered_df.shape[0]):
            if current_max.loc["query name"] != filtered_df.iloc[i,].loc["query name"] or \
               current_max.loc["target name"] != filtered_df.iloc[i,].loc["target name"]:
                ps.DataFrame(current_max).T.to_csv(out_f, sep = "\t",
                                                   index = False,
                                                   header = False,
                                                   mode = 'a')
                current_max = filtered_df.iloc[i,]
            else: 
                if current_max.iloc[13] < filtered_df.iloc[i,13]:
                    current_max = filtered_df.iloc[i,]
        #write last domain hit to disk
        ps.DataFrame(current_max).T.to_csv(out_f, sep = "\t",
                                           index = False,
                                           header = False,
                                           mode = 'a')

def make_outdir(filename):
    if filename is None:
        return None
    dirname = os.path.split(filename)[0]
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

def _main_par(x, hmm_name):
    infile_f,out_best_f = x
    main(infile_f, hmm_name, out_best_f)
    return out_best_f
        
def main(infile_f, hmm_name, out_best_f, out_excl_f=None, out_filt_f=None):
    """
    Params:
      infile_f : str, dat input file
      hmm_name : str, hmm model type (eg., pfam)
      out_best_f : str, output dat file
    """
    for F in [out_best_f, out_excl_f, out_filt_f]:
        make_outdir(F)    
    with open(infile_f) as inF:
        filtered_df = apply_thresholds(inF, hmm_name, out_filt_f, out_excl_f) 
    aggregate_domain_hits(filtered_df, out_best_f)
    
    logging.info('File written: {}'.format(out_best_f))
    return out_best_f
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser("parse hmmer output file and generated filtered best file")
    parser.add_argument("infile_f", default = sys.stdin, type = file,
                        help='a space-delimited output file from hmmer')
    parser.add_argument("out_best_f",
                        help='hmmer tab delimited file with hits removed according' + \
                        ' to the thresholds and only best domain hit retained')
    parser.add_argument("hmm_name", choices = ["pfam", "dbcan"],
                        help ='name of the HMM database; this will determine the' +
                        ' filtering applied')
    parser.add_argument("--out_excl_f",
                        help='domain hits filtered due to the applied thresholds')
    parser.add_argument("--out_filt_f",
                        help='hmmer tab delimited file with hits removed according' + \
                        ' to the thresholds')
    args = parser.parse_args()
    main(args.infile_f, args.hmm_name,
         out_filt_f = args.out_filt_f,
         out_excl_f = args.out_excl_f,
         out_best_f = args.out_best_f)
