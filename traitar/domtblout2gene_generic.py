#!/usr/bin/env python
"""script to create a summary matrix and a gene2hmm mapping from filtered and aggregated hmmer output files"""
import pandas as ps
def gene2hmm(domtblout_fs, acc_f, gene2hmm_out):
    """function to create a summary matrix and a gene2hmm mapping from filtered and aggregated hmmer output files"""
    #read accession file 
    accs = ps.read_csv(acc_f, sep = "\t", index_col = 0, header = None).index
    f = open(domtblout_fs, 'r')
    domtblout_list = [i.strip() for i in f.readlines()] 
    f.close()
    sum_df = ps.DataFrame(ps.np.zeros((len(domtblout_list), len(accs))))
    #set index to the actual files by getting rid of the preceding path
    sum_df.index = [i.split("/")[-1].replace("_filtered_best.dat", "") for i in domtblout_list]
    sum_df.columns = accs
    #create a gene2hmm dictionary
    gene2hmm = {}  
    for f in domtblout_list:
        gene_df = ps.read_csv(f, sep = "\t")
        if not f in gene2hmm:
            gene2hmm[f] = {}
        for i in gene_df.index:
            #append genes to gene2hmm dict
            gene = gene_df.loc[i, 'target name']
            query = gene_df.iloc[i, 4].split('.')[0]
            if gene not in gene2hmm[f]:
                gene2hmm[f][gene] = [query] 
            else:
                gene2hmm[f][gene].append(query)
            #add annotations to summary file
            sum_df.loc[f.split("/")[-1].replace("_filtered_best.dat", ""), query] += 1
    #write annotation summary to disk
    sum_df.to_csv(gene2hmm_out, sep = "\t")
    #write gene2hmm to disk 
    #for i in gene2hmm:
    #    with open(gene2hmm_out, 'w') as out: 
    #        out.write("%s\t%s\n" (i, ",".join(gene2hmm[i])))
    return sum_df, gene2hmm
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("generate summary matrix from the filtered best hmmer annotation files")
    parser.add_argument("outfile", help='summary matrix output file')
    parser.add_argument("in_filtered_best_fs", help='file with filtered.best file name per row')
    parser.add_argument("pfam_accession_f", help = "file with Pfam accessions as template for the summary matrix columns")
    args = parser.parse_args()
    sum_df, gene2hmm = gene2hmm(args.in_filtered_best_fs, args.pfam_accession_f, gene2hmm_out = args.outfile)
