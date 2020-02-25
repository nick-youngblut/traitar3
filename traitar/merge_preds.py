#!/usr/bin/env python
from __future__ import print_function
# import
import os
import sys
import logging
import argparse
from collections import defaultdict
## 3rd party
import pandas as ps

# functions
def flatten_df(df1, df2, name1, name2, out):
    with open(out, 'w') as f:
        f.write("sample\tphenotype\tphenotype model\n")
        for i in set(df1.index).union(set(df2.index)):
            for j in set(df1.columns).union(set(df2.columns)):
                if i in df1.index and j in df1.columns and  not df1.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df1.loc[i,j]), name1]) + "\n")
                if i in df2.index and j in df2.columns and not df2.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df2.loc[i,j]), name2]) + "\n")
    logging.info('File written: {}'.format(out))

def comb_preds(pred_files, primary_name, secondary_name, out_dir, k):
    #phypat only predictions
    m1_maj = ps.read_csv(pred_files['primary']['majority-vote'],
                         index_col = 0, sep = '\t')
    m1 = ps.read_csv(pred_files['primary']['single-votes'],
                     index_col = 0, sep = '\t')
    m2_maj = ps.read_csv(pred_files['secondary']['majority-vote'],
                         index_col = 0, sep = '\t')
    m2 = ps.read_csv(pred_files['secondary']['single-votes'],
                     index_col = 0, sep = '\t')
    #write to disk a single vote version of the predictions
    #phypat plus ml predictions
    #set NAs to zero
    m1[ps.isnull(m1)] = 0
    m2[ps.isnull(m2)] = 0
    #collect phenotypes that are shared / different in the two data sets
    #combine predictions of phypat and phypat+PGL into one matrix 
    m_columns = set(m1.columns).union(set(m2.columns))
    m = ps.DataFrame(ps.np.zeros((m1.shape[0], len(m_columns))))
    m_maj = ps.DataFrame(ps.np.zeros((m1.shape[0], len(m_columns))))
    m_maj.columns = m.columns = m_columns
    m_maj.index = m1.index
    m.index = m1.index
    #combine single voters
    m.loc[:, m1.columns] = m1
    m.loc[:, m2.columns] = m.loc[:, m2.columns] + m2
    #combine majority votes
    m_maj.loc[:, set(m1.columns).difference(set(m2.columns))] = m1_maj.loc[:, set(m1.columns).difference(set(m2.columns))] * 2
    m_maj.loc[:, set(m2.columns).difference(set(m1.columns))] = m2_maj.loc[:, set(m2.columns).difference(set(m1.columns))] 
    m_temp = m_maj.loc[:, set(m2.columns).intersection(set(m1.columns))]
    if not set(m2.columns).intersection(set(m1.columns)) == set():
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1)] = 3
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 0)] = 1
        m_temp[(m1_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 0) & (m2_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] == 1)] = 2
        m_maj.loc[:, set(m2.columns).intersection(set(m1.columns))] = m_temp

    outfile = os.path.join(out_dir, 'predictions_single-votes_combined.txt')
    m.to_csv(outfile, sep = "\t")
    logging.info('File written: {}'.format(outfile))
    outfile = os.path.join(out_dir, 'predictions_majority-vote_combined.txt')
    m_maj.to_csv(outfile, index_label = None, sep = "\t")
    logging.info('File written: {}'.format(outfile))
    outfile = os.path.join(out_dir, 'predictions_flat_majority-votes_combined.txt')
    flatten_df(m1_maj, m2_maj, primary_name, secondary_name, outfile)
    outfile = os.path.join(out_dir, 'predictions_flat_single-votes_combined.txt')
    flatten_df(m1, m2, primary_name, secondary_name, outfile)

def main(pred_files, primary_name,
         secondary_name, out_dir, voters):
    """ Main interface 
    Params:
      pred_files : {primary|secondary : {majority-vote|conservative-vote|single-votes : value}}
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    comb_preds(pred_files, primary_name, secondary_name, out_dir, voters)
    
def pred_file_get(phypat_dir, phypat_PGL_dir):
    pred_files = defaultdict(dict)
    pred_files['primary']['majority-vote'] = os.path.join(phypat_dir,
                                                          'predictions_majority-vote.txt')
    pred_files['primary']['single-votes'] = os.path.join(phypat_dir,
                                                          'predictions_single-votes.txt')
    pred_files['primary']['majority-vote'] = os.path.join(phypat_PGL_dir,
                                                          'predictions_majority-vote.txt')
    pred_files['primary']['single-votes'] = os.path.join(phypat_PGL_dir,
                                                          'predictions_single-votes.txt')
    return pred_files
    
if __name__ == "__main__":
    msg = 'Combine the misclassified samples of different phenotypes into data matrices'
    parser = argparse.ArgumentParser(msg)
    parser.add_argument("out_dir",
                        help='the output directory')
    parser.add_argument("phypat_dir",
                        help='directory with the phypat predictions')
    parser.add_argument("phypat_PGL_dir",
                        help='directory with the phyapt+PGL predictions')
    parser.add_argument("primary_name",
                        help='name of the primary phenotype model collection')
    parser.add_argument("secondary_name",
                        help='name of the secondary phenotype collection')
    parser.add_argument("-k", "--voters", default = 5, type = int,
                        help='number of classifiers used in the voting committee')
    args = parser.parse_args()

    pred_files = pred_file_get(args.phypat_dir, args.phypat_PGL_dir)
    main(pred_files, args.primary_name, args.secondary_name, args.out_dir, args.voters)
