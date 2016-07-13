import pandas as pd
import tarfile
class PhenotypeCollection:
    def __init__(self, archive_f):
        self.archive_f = archive_f
        self.tar = tarfile.open(archive_f, mode = "r:gz")
        info = self.tar.extractfile("config.txt")
        info_df = pd.read_csv(info, sep = "\t", index_col = 0)
        self.name = info_df.loc["archive_name", "value"]
        self.hmm_f = info_df.loc["hmm_f", "value"]
        self.hmm_name = info_df.loc["hmm_name", "value"]
    def get_pt2acc(self):
        pt2acc_f = self.tar.extractfile("pt2acc.txt")
        pt2acc = pd.read_csv(pt2acc_f, sep = "\t", index_col = 0,  encoding='utf-8')
        pt2acc.index = pt2acc.index.astype('U')
        return pt2acc
    def get_pt2id(self):
        pt2id_f = self.tar.extractfile("pt2acc.txt")
        pt2id = pd.read_csv(pt2id_f, sep = "\t", index_col = 1,  encoding='utf-8')
        pt2id.index = pt2id.index.astype('U')
        return pt2id
    def get_acc2pt(self):
        pt2acc_f = self.tar.extractfile("pt2acc.txt")
        pt2acc = pd.read_csv(pt2acc_f, sep = "\t", index_col = 1,  encoding='utf-8')
        pt2acc.index = pt2acc.index.astype('U')
        return pt2acc
    def get_pf2desc(self):
        pfam_mapping = pd.read_csv(self.tar.extractfile("pf2acc_desc.txt"), index_col = 0, sep = "\t")
        return pfam_mapping
    def get_name(self):
        pt2desc_f = self.tar.extractfile("pt2acc.txt")
        return self.name
    def get_archive_f(self):
        return self.archive_f
    def get_hmm_f(self):
        return self.hmm_f
    def get_hmm_name(self):
        return self.hmm_name
    def get_bias(self, pt):
        bias_f = self.tar.extractfile("%s_bias.txt"%(pt))
        bias = pd.read_csv(bias_f, sep = "\t", index_col = 0, header = None)
        return bias
    def get_predictors(self, pt):
        extracted_f = self.tar.extractfile("%s_feats.txt"%(pt))
        predictors = pd.read_csv(extracted_f, sep = "\t",  index_col = 0 )
        return predictors
    def get_majority_features(self,pt):
        extracted_f = self.tar.extractfile("%s_non-zero+weights.txt" % pt)
        feats = pd.read_csv(extracted_f, sep = "\t", index_col = 0)
        return feats

    def get_selected_features(self, pt, strategy, include_negative):
        pt2id = self.get_pt2id()
        pt_id = pt2id.loc[pt,].iloc[0]
        extracted_f = self.tar.extractfile("%s_feats.txt" % pt_id)
        feats = pd.read_csv(extracted_f, sep = "\t", index_col = 0)
        pos = feats.apply(lambda x: (x > 0).sum(), axis = 1)
        neg = feats.apply(lambda x: (x < 0).sum(), axis = 1)
        if strategy == "non-zero":
            if include_negative:
                out_feats = feats.loc[(pos > 0) | (neg > 0),]
            else: 
                out_feats = feats.loc[(pos > 0),]
        else:
            if include_negative:
                out_feats = feats.loc[((pos > 3) | (neg > 3) ), ]
            else:
                out_feats =  feats.loc[pos > 3, ]
        pf2desc = self.get_pf2desc()
        return pd.concat([out_feats.iloc[:, 0:5], pf2desc.loc[out_feats.index,]], axis = 1)
          

