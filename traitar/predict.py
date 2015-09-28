#!/usr/bin/env python
import os.path
import pandas as ps
import sys
import tarfile

"""predict new samples"""

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

def aggregate(pred_df, k):
    """employ different prediction strategies"""
    maj_pred_dfs = [ps.DataFrame(ps.np.zeros(shape = (pred_df.shape[0], pred_df.shape[1] / k))) for i in range(4)]
    for i in range(0, pred_df.shape[1], k):
        for j in range(4):
            maj_pred_dfs[j].columns.values[i / k] = pred_df.columns.values[i].split("_")[0] 
            maj_pred_dfs[j].columns = maj_pred_dfs[j].columns.values
            maj_pred_dfs[j].index = pred_df.index
        #majority vote
        maj_pred_dfs[0].iloc[:,i / k] = pred_df.iloc[:, i: i + k].apply(filter_pred, axis = 1, is_majority = True, k = k)
        maj_pred_dfs[1].iloc[:,i / k] = (pred_df.iloc[:, i: i + k].apply(filter_pred, axis = 1, is_majority = True, k = k) > 0).astype('int')
        #conservative vote
        maj_pred_dfs[2].iloc[:,i / k] = pred_df.iloc[:, i: i + k].apply(filter_pred, axis = 1, is_majority = False, k = k) 
        maj_pred_dfs[3].iloc[:,i / k] = (pred_df.iloc[:, i: i + k].apply(filter_pred, axis = 1, is_majority = False, k = k) > 0).astype('int')
    #maj_pred_dfs[0][maj_pred_dfs[0] == 0.0] = None 
    #maj_pred_dfs[1][maj_pred_dfs[1] == 0.0] = None 
    #maj_pred_dfs[2][maj_pred_dfs[2] == 0.0] = None 
    #maj_pred_dfs[3][maj_pred_dfs[3] == 0.0] = None 
    return maj_pred_dfs
    

    
def majority_predict(pt, model_tar, test_data, k, bias_weight = 1):
    """predict the class label based on a committee vote of the models in models""" 
    #TODO if the classification model is trained on non binary data we would require the same normalization applied to the training data 
    #binarize
    test_data_n = (test_data > 0).astype('int')
    #check if classification model exists for the requested phenotype
    try: 
        extracted_f = model_tar.extractfile("%s_feats.txt"%(pt))
        bias_f = model_tar.extractfile("%s_bias.txt"%(pt))
        bias = ps.read_csv(bias_f, sep = "\t", index_col = 0, header = None)
        predictors = ps.read_csv(extracted_f, sep = "\t",  index_col = 0 )
    except KeyError:
        return ps.DataFrame()
    #build prediction matrix with the k best models
    preds = ps.np.zeros((test_data.shape[0], k))
    for i in range(k):
        #preds[:, i] = test_data_n.dot(predictors.iloc[:, i].iloc[test_data_n.columns, ].values)
        #print predictors.index
        #print test_data_n.columns
        preds[:, i] = bias.iloc[i, 0] * bias_weight +  test_data_n.dot(predictors.iloc[:, i].iloc[test_data_n.columns, ].values)
        pred_df = ps.DataFrame(preds, index = test_data.index)
    #set column names
    pred_df.columns = ["%s_%s" %(pt, predictors.columns[i].split("_")[0]) for i in range(k)]
    return pred_df

def annotate_and_predict(pt_range, model_tar, summary_f, pfam_pts_mapping_f, out_dir, k):
    """load annotation previously generated with HMMER and predict phenotypes with phypat and phypat+GGL models"""
    pred_df = ps.DataFrame()
    #read pfam to description file
    pfam_pts_mapping = ps.read_csv(model_tar.extractfile(pfam_pts_mapping_f), header=None, index_col = 0, sep = "\t") 
    pfam_pts_mapping.index = pfam_pts_mapping.index.values - 1 
    #read annotation file
    m = ps.read_csv(summary_f, sep="\t", index_col = 0)
    #restrict to those pfams contained in the model
    m_red = m.loc[:, pfam_pts_mapping.iloc[:-93, 0] ].astype('int')
    m_red.columns = pfam_pts_mapping.index[:-93]
    for pt in range(pt_range[0], pt_range[1] + 1):
        #predict an append to prediction df
        preds = majority_predict(pt, model_tar, m_red, k)
        pred_df = ps.concat([pred_df, preds], axis = 1)
    pred_df.index = m_red.index
    pred_df.to_csv("%s/predictions_raw.txt"%out_dir, sep = "\t", float_format='%.3f')
    #aggregate predictions
    out = ["majority-vote_mean-score", "majority-vote", "scores_conservative-vote_mean-score", "conservative-vote"]
    aggr_dfs = aggregate(pred_df, k)
    for i in range(len(out)):
        aggr_dfs[i].columns = pfam_pts_mapping.loc[aggr_dfs[i].columns, ].iloc[:,0]
        aggr_dfs[i].index = m_red.index
        if i % 2 == 0:
            aggr_dfs[i].to_csv("%s/predictions_%s.txt"%(out_dir, out[i]), sep = "\t", float_format='%.3f')
        else:
            aggr_dfs[i].to_csv("%s/predictions_%s.txt"%(out_dir, out[i]), sep = "\t", float_format='%.0f')
    return pred_df


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("predict phenotypes from hmmer Pfam annotation")
    parser.add_argument("model_tar", help='directory with the phenotype predictors')
    parser.add_argument("out_dir", help='directory for the phenotype predictions')
    parser.add_argument("phenotype_range", help='specify which phenotypes shall be predicted e.g. 8476-8568')
    parser.add_argument("annotation_matrix", help='summary file with pfams')
    parser.add_argument("pfam_pts_mapping_f", help='mapping from pfam and phenotype id to universal accessions')
    parser.add_argument("-k", "--voters", default = 5, help='use this number of voters for the classification', type = int)
    args = parser.parse_args()
    pt1, pt2 = [int(i) for i in args.phenotype_range.split("-")]
    annotate_and_predict((pt1, pt2), tarfile.open(args.model_tar, mode = "r:gz"), args.annotation_matrix,args.pfam_pts_mapping_f,  args.out_dir, args.voters) 
