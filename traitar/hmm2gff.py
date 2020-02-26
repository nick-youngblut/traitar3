#!/usr/bin/env python
from __future__ import print_function
""" produce gff file from a hmmer hit file """

# import
## batteries
import os
import sys
import tarfile
import warnings
import logging
from io import StringIO     
## 3rd party
import pandas as ps
## application
from traitar.PhenotypeCollection import PhenotypeCollection
warnings.filterwarnings("ignore", category=FutureWarning) 

def get_protein_acc(fasta_id):
    """extract from a fasta id like gi|17233466|ref|NP_490586.1|
    the protein accession NP_490586.1"""
    return fasta_id.split("|")[3]
    
def read_rel_feats(pt_models, pts):
    """read in pfams relevant for all predicted phenotypes"""
    #pfam to phenotype mapping
    pf2pt = {}
    #store relevant pfam lists
    rel_feat_list = []
    #extract pfam lists
    for pt in pts:
        x = pt_models.get_selected_features(pt, strategy = "majority",
                                            include_negative = False)
        rel_feat_list.append(x)
    #process pfam lists
    for pfs, pt in zip(rel_feat_list, pts): 
        for pf in pfs.index:
            if not pf in pf2pt:
                if not type(pfs.index.values) == type(ps.Series()):
                    pf2pt[pf] = [(pt, pfs.loc[pf, "cor"])]
                else: 
                    print(pf, pt)
                    print(pfs.loc[pf, "Pfam_acc"])
            else:
                if not pt in pf2pt[pf]:
                    pf2pt[pf].append((pt, pfs.loc[pf, "cor"]))
    return pf2pt

def read_gff(gff_file, mode):
    """read a gene calling gff"""
    f = open(gff_file)
    gene_dict = {}
    i = 0
    fasta = False
    for l in f:
        i += 1
        # print i
        if l.startswith("#") or l.startswith("\n"):
            continue
        else:
            if mode == "metagenemark":
                read_genemark_entry(l, gene_dict)
            elif mode == "prodigal":
                read_prodigal_entry(l, gene_dict)
            elif mode == "ncbi":
                read_ncbi_entry(l, gene_dict)
            elif mode == "refseq":
                read_refseq_entry(l, gene_dict)
            elif mode == "img":
                read_img_entry(l, gene_dict)
            elif mode == "genbank":
                #in genbank gene gffs, the original nucleotide sequence
                ## can be appended to the gff
                if l.startswith(">"):
                    fasta = True 
                    continue
                if fasta and l.startswith(('a','g', 'c', 't')):
                    continue
                else:
                    fasta = False
                read_genbank_entry(l, gene_dict)
    f.close()
    return gene_dict

def read_genemark_entry(l, gene_dict):
    """read and parse one line from a genemark gff"""
    elems = l.strip().split("\t")
    gene_dict[
        elems[8].split(" ")[0].strip(",").replace("=", "_")] = (
        elems[0].split(" ")[0], int(
            elems[3]), int(
            elems[4]), elems[6])

def read_prodigal_entry(l, gene_dict):
    """read and parse one line from a prodigal gff"""
    elems = l.strip().split("\t")
    attrs = dict(
        [(i.split("=")[0], i.split("=")[1])
        for i in elems[8].strip(";").split(";")])
    gene_dict["{}_{}".format(elems[0],
                             attrs["ID"].split("_")[1])] = (
                                 elems[0], int(elems[3]),
                                 int(elems[4]), elems[6])

def read_img_entry(l, gene_dict):
    """read and parse one line from an IMG gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        try:
            attrs = dict(
                [(i.split("=")[0], i.split("=")[1]) for i in elems[8].strip(";").split(";")])
        except IndexError:
            sys.stderr.write("something went wrong in line; skipping \n{}\n".format(l))
            return 
        gene_dict[attrs["ID"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])


def read_ncbi_entry(l, gene_dict):
    """read and parse one line from a ncbi gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["Name"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])

def read_refseq_entry(l, gene_dict):
    """read and parse one line from a refseq gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["ID"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])

def read_genbank_entry(l, gene_dict):
    """read and parse one line from a genbank gff"""
    elems = l.strip().split("\t")
    if elems[2] == "CDS":
        attrs = dict(
            [(i.split("=")[0], i.split("=")[1])
             for i in elems[8].split(";")])
        gene_dict[
            attrs["locus_tag"]] = (
            elems[0], int(
                elems[3]), int(
                elems[4]), elems[6])


def get_coords(gene_start, gene_end, ali_from, ali_to, strand):
    n_start = ali_from * 3 - 2
    n_end = ali_to * 3 - 2
    # if strand == "+":
    return (gene_start + n_start, gene_start + n_end)


def write_hmm_gff(hmmer_file, out_gff_dir, gene_dict, sample, skip_genes,
                  mode, pt_models, predicted_phenotypes, description = False):
    #read in pfam accession to description mapping
    if description:
        acc2desc = pt_models.get_pf2desc()
        id2desc = pt_models.get_pt2acc()
        #change dtype of index to string
        id2desc.index = id2desc.index.values.astype(str)
    out_table = []
    with open(hmmer_file, 'r') as hmmer_fo:
        # skip param line and header
        hmmer_fo.readline()
        hmmer_fo.readline()
        hmmer_fo.readline()
        rel_feats_dict = None
        # read in relevant which shall be highlighted
        pts = predicted_phenotypes.split(",") 
        if pts is not None:
            pf2pt = read_rel_feats(pt_models, pts)
        if not gene_dict is None:
            # open global output table
            # open out gffs for reading
            infile = "{}/{}_{}_important_features.gff"
            out_gffs = dict([(i, open(infile.format(out_gff_dir, sample, id2desc.loc[i,].iloc[0].replace(" ", "_").replace("(", "_").replace("(", "_")), 'w')) for i in pts])
            infile = "{}/{}_complete_annotation.gff"
            out_gffs["Pfams"] = open(infile.format(out_gff_dir, sample), 'w')
            for i in out_gffs:
                out_gffs[i].write("""##gff-version 3\n""")
        for l in hmmer_fo:
            elems = l.strip().split("\t")
            # change pfam accession PF00001.3 into PF00001
            gid, acc, name, ali_from, ali_to, ieval = [
                get_protein_acc(elems[0]) if mode == "ncbi" or mode == "refseq" else elems[0], 
                    elems[4].split(".")[0] if pt_models.get_hmm_name() == "pfam" else elems[3].split(".")[0],
                    elems[3], int(elems[17]), int(elems[18]), elems[11]]
            if gene_dict is not None and gid not in gene_dict: 
                if not skip_genes:
                    logging.info('"{}" not found in input orf gff file',format(gid))
                    if mode == "ncbi":
                        logging.info("entry might be outdated, skipping")
                        continue
                    if mode == "img":
                        msg = 'might be a ";" delimiter within an attribute; sample is {}'
                        logging.info(msg.format(sample))
                        continue
                else:
                    continue
            colour = ""
            # check if this feature is among the highest ranking features and add
            # colour and rank if so 
            description_string = ""
            if gene_dict is not None and description and acc in acc2desc.index:
                contig_name, gene_start, gene_end, strand = gene_dict[gid]
                hmm_coords = get_coords(gene_start, gene_end, ali_from, ali_to, strand)
                description_string = ";Description={}".format(acc2desc.loc[acc,].iloc[0])
                #gff_line = "%s\tHMMER\t%s\t%s\t%s\t%s\t%s\t.\tParent=%s;Name=%s;Note=%s%s%s\n"
                gff_line = '\t'.join([contig_name, 'HMMER', pt_models.get_name(),
                                      hmm_coords[0], hmm_coords[1], ieval, strand,
                                      '.', 'Parent={}'.format(gid),
                                      'Name={}'.format(acc),
                                      'Note={}{}{}'.format(name, description_string, colour)
                ])
                out_gffs["Pfams"].write(gff_line + '\n')
            if pf2pt is None or acc in pf2pt:
                for i in pf2pt[acc]:
                    if not gene_dict is None:
                        out_gffs[i[0]].write(gff_line)
                    x = [id2desc.loc[i[0],"accession"], gid, acc,
                         acc2desc.loc[acc, "description"], i[1]]
                    out_table.append('\t'.join([str(y) for y in x]))
    if not len(out_table) == 0:
        out_table_df = ps.read_csv(StringIO("".join(out_table)),
                                   sep = "\t", header = None)
        out_table_df.columns = ["Phenotype", "gene_id", "acc", "description", "cor"]
        outfile = "{}/{}_important_features.dat".format(out_gff_dir, sample)
        x = out_table_df.sort_values(by = ["Phenotype", "cor"],
                                     ascending = [True, False])
        x.to_csv(outfile, sep = "\t", index = None)
    #close out gffs
    if not gene_dict is None:
        for i in out_gffs:
            out_gffs[i].close()

def run(in_file, out_gff_f, gene_gff_f, sample, gene_gff_mode, pt_models, predicted_pts = None):
    """run feature mapping"""
    # output directory
    outdir = os.path.split(out_gff_f)[0]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # reading in gene prediction track
    gene_dict = None
    if not gene_gff_f is None:
        gene_dict = read_gff(gene_gff_f, gene_gff_mode)
    # writing hmm-gff
    write_hmm_gff(in_file, out_gff_f, gene_dict, sample, False,
                  gene_gff_mode, pt_models, predicted_pts, True)

if __name__ == "__main__":
    import argparse
    msg = 'Map features contributing to the classfication back to the functional annotation and gene prediction'
    parser = argparse.ArgumentParser(msg)
    parser.add_argument("input_hmm",
                        help= 'input HMMer file')
    parser.add_argument("output_gff_dir",
                        help= 'output GFF file')
    parser.add_argument("sample",
                        help= 'sample file')
    parser.add_argument("model_tar",
                        help = "tar.gz file with relevant features etc.")
    parser.add_argument("predicted_pts",
                        help='phenotypes that were predicted for the given samples')
    parser.add_argument("--gene_gff", default = None,
                        help= "gene prediction track")
    parser.add_argument("--gene_gff_type", default = None,
                        choices = ["img", "prodigal", "ncbi", "refseq",
                                   "metagenemark", "genbank"],
                        help= "origin of the gene prediction (Prodigal, NCBI, metagenemark)")
    args = parser.parse_args()
    pt_models =  PhenotypeCollection(args.model_tar)
    run(args.input_hmm, args.output_gff_dir, args.gene_gff,
        args.sample, args.gene_gff_type,  pt_models, args.predicted_pts)
