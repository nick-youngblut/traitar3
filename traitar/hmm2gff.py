#!/usr/bin/env python
import sys
import pandas as ps
import tarfile
""" produce gff file from a hmmer hit file """


def get_protein_acc(fasta_id):
    """extract from a fasta id like gi|17233466|ref|NP_490586.1| the protein accession NP_490586.1"""
    return fasta_id.split("|")[3]


def read_rel_feats(rel_feats_f):
    f = open(rel_feats_f, 'r')
    rel_feats_dict = {}
    i = 1
    for l in f:
        acc, _ = l.split("\t")
        rel_feats_dict[acc] = i
        i += 1
    f.close()
    return rel_feats_dict

# read in gene prediction gff file as a dictionary


def read_gff(gff_file, mode):
    """read a gene calling gff"""
    f = open(gff_file)
    gene_dict = {}
    i = 0
    for l in f:
        i += 1
        # print i
        if l.startswith("#") or l.startswith("\n"):
            continue
        else:
            if mode == "gm":
                read_genemark_entry(l, gene_dict)
            elif mode == "prodigal":
                read_prodigal_entry(l, gene_dict)
            elif mode == "ncbi":
                read_ncbi_entry(l, gene_dict)
    f.close()
    return gene_dict


def read_genemark_entry(l, gene_dict):
    """read and parse one line from a genemark gff"""
    elems = l.strip().split("\t")
    gene_dict[
        elems[8].replace(
            " ", "_")] = (
        elems[0].split(" ")[0], int(
            elems[3]), int(
            elems[4]), elems[6])

def read_prodigal_entry(l, gene_dict):
    """read and parse one line from a genemark gff"""
    elems = l.strip().split("\t")
    #print elems[8].split(";"), len(elems)
    attrs = dict(
        [(i.split("=")[0], i.split("=")[1])
        for i in elems[8].strip(";").split(";")])
    gene_dict["%s_%s"% (elems[0], attrs["ID"].split("_")[1])] = (
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


def get_coords(gene_start, gene_end, ali_from, ali_to, strand):
    # print gene_start,gene_end,ali_from,ali_to, strand
    n_start = ali_from * 3 - 2
    n_end = ali_to * 3 - 2
    # if strand == "+":
    return (gene_start + n_start, gene_start + n_end)
    # else:
    #    return (gene_end - n_end, gene_end  - n_start)


def write_hmm_gff(hmmer_file, out_gff_f, gene_dict, skip_genes, mode, rel_feats_f, description = False):
    #read in pfam accession to description mapping
    if description:
        id2desc = ps.read_csv(model_tar.extractfile(pfam_pts_mapping_f), header=None, index_col = 0, sep = "\t") 
    with open(hmmer_file, 'r') as hmmer_fo:
        # skip param line and header
        hmmer_fo.readline()
        hmmer_fo.readline()
        hmmer_fo.readline()
        rel_feats_dict = None
        # read in relevant which shall be highlighted
        if rel_feats_f is not None:
            rel_feats_dict = read_rel_feats(rel_feats_f)
        #print rel_feats_dict
        # open out gff for reading
        with open(out_gff_f, 'w') as out_gff:
            out_gff.write("""##gff-version 3\n""")
            for l in hmmer_fo:
                elems = l.strip().split("\t")
                # change pfam accession PF00001.3 into PF00001
                gid, acc, name, ali_from, ali_to, ieval = [
                    get_protein_acc(elems[0]) if mode == "ncbi" else elems[0], elems[4].split(".")[0], elems[3], int(elems[18]), int(elems[18]), elems[11]]

                if gid not in gene_dict:
                    if not skip_genes:
                        print >> sys.stderr, gid, "not found in input orf gff file"
                        if mode == "ncbi":
                            print >> sys.stderr, "entry might be outdated, skipping"
                            continue
                        sys.exit(1)
                    else:
                        continue
                contig_name, gene_start, gene_end, strand = gene_dict[gid]
                hmm_coords = get_coords(gene_start, gene_end, ali_from, ali_to, strand)
                colour = ""
                # check if this feature is among the highest ranking features and add
                # colour and rank if so 
                #print acc, rel_feats_dict 
                if rel_feats_dict is None or acc in rel_feats_dict:
                    description_string = ""
                    if description:
                        description_string = ";Description=%s" % (id2desc.loc[acc,].iloc[0])
                    #if not rel_feats_dict is None:
                    #    colour = ";colour=3;rank=%s" % (rel_feats_dict[acc])
                    out_gff.write(
                        "%s\tHMMER\tPfam\t%s\t%s\t%s\t%s\t.\tParent=%s;Name=%s;Note=%s%s%s\n" %
                        (contig_name, hmm_coords[0], hmm_coords[1], ieval, strand, gid, acc, name, description_string, colour))

# extract the information from the hmmer hit file which goes into the
# output gff file and extract the original position in the genome/contig
def run(in_file, out_gff_f, gene_gff_f,  gene_gff_mode, rel_feats_f=None):
    gene_dict = read_gff(gene_gff_f, gene_gff_mode)
    write_hmm_gff(in_file, out_gff_f, gene_dict, False, gene_gff_mode, rel_feats_f)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("predict phenotypes from hmmer Pfam annotation")
    parser.add_argument("input_hmm", help= 'input HMMer file')
    parser.add_argument("gene_gff", help= 'gene prediction gff file')
    parser.add_argument("output_gff", help= 'output GFF file')
    parser.add_argument("gene_gff_mode", choices = ["prodigal", "ncbi", "metagenemark"], help= "origin of the gene prediction (Prodigal, NCBI, metagenemark)")
    parser.add_argument("--relevant_features", "-r", help='file with some relevant annotation features', default = None)
    args = parser.parse_args()
    run(args.input_hmm, args.output_gff, args.gene_gff, args.gene_gff_mode,  args.relevant_features)
