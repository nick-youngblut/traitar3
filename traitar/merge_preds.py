#!/usr/bin/env python
import pandas as ps
def flatten_df(df1, df2, out):
    with open(out, 'w') as f:
        for i in set(df1.index).union(set(df2.index)):
            for j in set(df1.columns).union(set(df2.columns)):
                if i in df1.index and j in df1.columns and  not df1.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df1.loc[i,j]), "phypat"]) + "\n")
                if i in df2.index and j in df2.columns and not df2.loc[i,j] == 0:
                    f.write("\t".join([str(i), str(j), str(df2.loc[i,j]), "phypat+GGL"]) + "\n")

def comb_preds(phypat_dir, phypat_GGL_dir, out_dir, k):
    #phypat only predictions
    m1_raw = ps.read_csv("%s/predictions_raw.txt"%phypat_dir, index_col = 0, sep = "\t")
    m1 = ps.DataFrame(ps.np.zeros((m1_raw.shape[0], m1_raw.shape[1] / 5)))
    m1.index = m1_raw.index
    for i in range(m1_raw.shape[1] / 5):
        m1.iloc[:, i] = m1_raw.iloc[:, i*5 : i*5 + 5].apply(lambda x: (x > 0).sum(), axis = 1 )
    m1_scores = ps.read_csv("%s/predictions_majority-vote_mean-score.txt"%phypat_dir, index_col = 0, sep = "\t")
    #write to disk a single vote version of the predictions
    m1.columns = m1_scores.columns
    m1.index = m1_scores.index
    m1.to_csv("%s/predictions_single-votes.txt"%phypat_dir, sep = "\t")
    #phypat plus ml predictions
    m2_raw = ps.read_csv("%s/predictions_raw.txt"%phypat_GGL_dir, index_col = 0,  sep = "\t")
    m2 = ps.DataFrame(ps.np.zeros((m2_raw.shape[0], m2_raw.shape[1] / 5)))
    m2.index = m1_raw.index
    #aggregate single predictor outcomes
    for i in range(m2_raw.shape[1] / 5):
        m2.iloc[:, i] = m2_raw.iloc[:, i*5 : i*5 + 5].apply(lambda x: (x > 0).sum(), axis = 1 )
    m2_scores = ps.read_csv("%s/predictions_majority-vote_mean-score.txt"%phypat_GGL_dir, index_col = 0, sep = "\t")
    #write to disk a single vote version of the predictions
    m2.index = m2_scores.index
    m2.columns = m2_scores.columns
    m2.to_csv("%s/predictions_single-votes.txt"%phypat_GGL_dir, sep = "\t")
    #set NAs to zero
    m1[ps.isnull(m1)] = 0
    m2[ps.isnull(m2)] = 0
    m1_scores[ps.isnull(m1_scores)] = 0
    m2_scores[ps.isnull(m2_scores)] = 0
    #generate majority vote matrices
    m1_maj = m1.copy()
    m2_maj = m2.copy()
    m1_maj[m1 < (k/2 + 1)] = 0
    m1_maj[m1 >= (k/2 + 1)] = 1
    m2_maj[m2 < (k/2 + 1)] = 0
    m2_maj[m2 >= (k/2 + 1)] = 1
    #collect phenotypes that are shared / different in the two data sets

    
    #combine predictions of phypat and phypat+GGL into one matrix 
    m = m1 + m2 
    m_maj = m1_maj.copy()
    m_maj[m2_maj == 1] = 2
    m_maj[(m1_maj == 1) & (m2_maj == 1)] = 3
    #m_maj = ((m1_maj + m2_maj) > 0).astype('int')
    m_scores = (m1_scores +  m2_scores)/2
    
    #write named majority, single vote and mean score version of combined, phypat and phypat+GGL 
    m.to_csv("%s/predictions_single-votes_combined.txt"%out_dir, sep = "\t")
    m_maj.to_csv("%s/predictions_majority-vote_combined.txt"%out_dir, index_label = None, sep = "\t")
    m_scores.to_csv("%s/predictions_majority-vote_mean-score_combined.txt"%out_dir, index_label = None, sep = "\t")
    #m1.to_csv("%s/predictions_aggr_phypat_majority-vote_named.csv"%out_dir, index_col = 0, sep = "\t")
    #m1_scores.to_csv("%s/predictions_aggr_phypat_majority-vote_mean-score_named.csv"%out_dir, index_col = 0, sep = "\t")
    #m2.to_csv("%s/predictions_aggr_phypat+ml_majority-vote_named.csv"%out_dir, index_col = 0,  sep = "\t")
    #m2_scores.to_csv("%s/predictions_aggr_phypat+ml_majority-vote_mean-score_named.csv"%out_dir, index_col = 0,  sep = "\t")
    flatten_df(m1_maj, m2_maj, "%s/predictions_flat_majority-votes_combined.txt"%out_dir)
    flatten_df(m1, m2, "%s/predictions_flat_single-votes_combined.txt"%out_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("combine the misclassified samples of different phenotypes into data matrices")
    parser.add_argument("out_dir",help='the output directory')
    parser.add_argument("phypat_dir",help='directory with the phypat predictions')
    parser.add_argument("phypat_GGL_dir",help='directory with the phyapt+GGL predictions')
    parser.add_argument("-k", "--voters", default = 5,  help='number of classifiers used in the voting committee', type = int)
    args = parser.parse_args()
    comb_preds(args.phypat_dir, args.phypat_GGL_dir, args.out_dir, args.voters)
