#!/usr/bin/env python
from __future__ import print_function
"""predict new samples"""

# import
## batteris
import os.path
import sys
import logging
import argparse
import multiprocessing as mp
from functools import partial
## 3rd party
import pandas as ps
## application
from traitar.PhenotypeCollection import PhenotypeCollection

# init functions
def filter_pred(scores, is_majority, k):
    """either do majority vote aggregation or conservative all or nothing vote"""
    if is_majority:
        if (scores > 0).sum() >= (k/2 + 1):
            return scores[scores >= 0].mean()
        else:
            return ps.np.NaN
    else:
        if (scores > 0).all():
            return scores.mean()
        else:
            return ps.np.NaN

def aggregate(pred_df, k, out_dir, pt2acc):
    """employ different prediction strategies"""
    out = ["majority-vote", "conservative-vote", "single-votes"]
    maj_pred_dfs = ps.DataFrame(ps.np.zeros(shape = (pred_df.shape[0], pred_df.shape[1] // k)),
                                columns = [str(i) for i in range(pred_df.shape[1]//k)])
    maj_pred_dfs.index = pred_df.index
    maj_pred_dfs_columns = maj_pred_dfs.columns.tolist()
    for i in range(pred_df.shape[1] // k):
        maj_pred_dfs_columns[i] = pred_df.columns.values[i * k].split("_")[0] 
        maj_pred_dfs.iloc[:, i] = pred_df.iloc[:, (i * k) : (i * k + k)].apply(lambda x: (x > 0).sum(), axis = 1).astype('int')
    maj_pred_dfs.columns = pt2acc.loc[maj_pred_dfs_columns, :].iloc[:, 0]
    #majority vote
    out_files = {}
    outfile = os.path.join(out_dir, 'predictions_{}.txt'.format(out[2]))
    maj_pred_dfs.to_csv(outfile, sep = '\t', float_format='%.0f', encoding = 'utf-8')
    logging.info('File written: {}'.format(outfile))
    out_files[out[2]] = outfile

    outfile = os.path.join(out_dir, 'predictions_{}.txt'.format(out[0]))
    (maj_pred_dfs >= k/2 + 1).astype('int').to_csv(outfile, sep = '\t',
                                                   float_format='%.0f',
                                                   encoding = 'utf-8')
    logging.info('File written: {}'.format(outfile))
    out_files[out[0]] = outfile
    
    outfile = os.path.join(out_dir, 'predictions_{}.txt'.format(out[1]))
    (maj_pred_dfs == k).astype('int').to_csv(outfile, sep = '\t',
                                             float_format='%.0f',
                                             encoding = 'utf-8')
    logging.info('File written: {}'.format(outfile))
    out_files[out[1]] = outfile
    
    return out_files
        
def majority_predict(x, test_data, k, bias_weight = 1):
    """predict the class label based on a committee vote of the models in models""" 
    #TODO if the classification model is trained on non binary data,
    ## we would require the same normalization applied to the training data
    pt,bias,predictors = x
    if bias is None or predictors is None:
        return ps.DataFrame()
    #binarize
    test_data_n = (test_data > 0).astype('int')
    #build prediction matrix with the k best models
    preds = ps.np.zeros((test_data.shape[0], k))
    for i in range(k):
        preds[:, i] = bias.iloc[i, 0] * bias_weight + \
                      predictors.iloc[:, i].dot(test_data_n.loc[:, predictors.iloc[:, i].index].T)
        pred_df = ps.DataFrame(preds, index = test_data.index)
    #set column names
    pred_df.columns = ['{}_{}'.format(pt, predictors.columns[i].split("_")[0]) for i in range(k)]
    return pred_df

def get_pt_info(pt, pt_models):
    try: 
        bias = pt_models.get_bias(pt)
        predictors = pt_models.get_predictors(pt)
    except KeyError:
        bias =  None
        predictors = None
    return [pt, bias, predictors]

def annotate_and_predict(pt_models, summary_f, out_dir, k, cpus=1):
    """load annotation previously generated with HMMER and predictphenotypes
    with phypat and phypat+GGL models"""
    #pred_df = ps.DataFrame()
    #read pfam to description file
    pfam_mapping = pt_models.get_pf2desc()
    pt_mapping = pt_models.get_pt2acc()
    #read annotation file
    m = ps.read_csv(summary_f, sep="\t", index_col = 0)
    #restrict to those pfams contained in the model
    # predicting for each model in parallel
    ## getting ps.dataframes from pt_models() object
    pt_info = []
    for pt in pt_mapping.index:
        pt_info.append(get_pt_info(pt, pt_models))
    ## running majority_predict (in parallel)
    msg = 'Running {} predictions with {} threads'
    logging.info(msg.format(len(pt_info), cpus))
    func = partial(majority_predict, test_data = m, k = k, bias_weight = 1)
    if cpus > 1:
        pool = mp.Pool(cpus)
        preds = pool.map(func, pt_info)
    else:
        preds = map(func, pt_info)
    # formatting predictions
    pred_df = ps.concat(preds, axis = 1)
    pred_df.index = m.index
    # writing output
    outfile = os.path.join(out_dir, 'predictions_raw.txt')
    pred_df.to_csv(outfile, sep = "\t", float_format='%.3f')
    logging.info('File written: {}'.format(outfile))
        
    #aggregate predictions
    out_files= aggregate(pred_df, k, out_dir, pt_mapping)
    return out_files

def main(pt_models, annotation_matrix, out_dir, voters, cpus):
    """ Main interface """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    pt_models = PhenotypeCollection(pt_models)
    out_files = annotate_and_predict(pt_models, annotation_matrix,
                                     out_dir, voters, cpus)
    return out_files
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser("predict phenotypes from hmmer Pfam annotation")
    parser.add_argument("pt_models",
                        help='archive with the phenotype predictors')
    parser.add_argument("out_dir",
                        help='directory for the phenotype predictions')
    parser.add_argument("annotation_matrix",
                        help='summary file with pfams')
    parser.add_argument("-k", "--voters", default = 5, type = int,
                        help='use this number of voters for the classification')
    parser.add_argument("-c", "--cpus", default = 1, type = int,
                        help='number of parallel processes')
    args = parser.parse_args()
    main(args.pt_models, args.annotation_matrix,
         args.out_dir, args.voters, args.cpus)
