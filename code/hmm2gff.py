import figfam2gff
import sys
import pandas as ps
# produce gff file from a hmmer hit file


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
    #print gene_dict
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


def write_hmm_gff(
        hmmer_file,
        out_gff_f,
        gene_dict,
        skip_genes,
        mode,
        anno_type,
        rel_feats_f,
        description = True):
    #read in pfam accession to description mapping
    if description:
        id2desc = ps.read_csv("/net/metagenomics/ref_data/HMMs/pfam_27.0/profile_IDs_to_description_mapping.txt", index_col = 0, header = None, sep = "\t")
    f = open(hmmer_file, 'r')
    # skip param line and header
    f.readline()
    f.readline()
    rel_feats_dict = None
    # read in relevant which shall be highlighted
    if rel_feats_f is not None:
        rel_feats_dict = read_rel_feats(rel_feats_f)
        #print rel_feats_dict
    # open out gff for reading
    out_gff = open(out_gff_f, 'w')
    out_gff.write("""##gff-version 3\n""")
    for l in f:
        elems = l.strip().split("\t")
        # change pfam accession PF00001.3 into PF00001
        gid, acc, name, ali_from, ali_to, ieval = [
            get_protein_acc(
                elems[0]) if mode == "ncbi" else elems[0], elems[4].split(".")[0], elems[3], int(
                elems[17]), int(
                elems[18]), elems[11]]

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
        # print acc
        # check if this feature is among the highest ranking features and add
        # colour and rank if so (MAMBA project specific
        #print acc, rel_feats_dict 
        if rel_feats_dict is None or acc in rel_feats_dict:
            description_string = ""
            if description:
                description_string = ";Description=%s" % (id2desc.loc[acc,].iloc[0])
            #if not rel_feats_dict is None:
            #    colour = ";colour=3;rank=%s" % (rel_feats_dict[acc])
            out_gff.write(
                "%s\tHMMER\t%s\t%s\t%s\t%s\t%s\t.\tParent=%s;Name=%s;Note=%s%s%s\n" %
                (contig_name,
                anno_type,
                hmm_coords[0],
                 hmm_coords[1],
                ieval,
                strand,
                gid,
                acc,
                name,
                description_string,
                colour))
    f.close()
    out_gff.close()


def write_motif_gff(motif_file,  out_gff_f, gene_dict,  skip_genes=False):
    """process a prosite motif gff file"""
    # TODO account for ncbi gff
    f = open(motif_file, 'r')
    # skip param line and header
    # open out gff for reading
    out_gff = open(out_gff_f, 'w')
    out_gff.write("""##gff-version 2\n""")
    for l in f:
        elems = l.strip().split("\t")
        # change pfam accession PF00001.3 into PF00001
        gid, ps_scan, acc, start, end, attributes = [
            elems[0], elems[1], elems[2], int(
                elems[3]), int(
                elems[4]), elems[8]]
        if gid not in gene_dict:
            if not skip_genes:
                print gid, "not found in input orf gff file"
                sys.exit(1)
            else:
                continue
        contig_name, gene_start, gene_end, strand = gene_dict[gid]
        motif_coords = get_coords(gene_start, gene_end, start, end, strand)
        # check if this feature is among the highest ranking features and add
        # colour and rank if so (MAMBA project specific
        out_gff.write(
            "%s\t%s\tmotif\t%s\t%s\t.\t%s\t.\t%s ; Parent \"%s\"\n" %
            (contig_name,
             ps_scan,
             motif_coords[0],
             motif_coords[1],
             strand,
             attributes,
             gid))
    f.close()
    out_gff.close()


def write_blast_gff(blast_file,  out_gff_f, gene_dict):
    """process a blast tabular output file"""
    # TODO account for ncbi gff
    f = open(blast_file, 'r')
    # open out gff for reading
    out_gff = open(out_gff_f, 'w')
    out_gff.write("""##gff-version 3\n""")
    # dictionary to only retain best hit
    gid_uq_dict = {}
    for l in f:
        elems = l.strip().split("\t")
        # print elems
        gid, sseqid, qstart, qend, evalue, bitscore, gene_product_subject = [
            elems[0], elems[1], int(
                elems[6]), int(
                elems[7]), elems[10], elems[11], elems[12]]
        if not gid in gid_uq_dict:
            gid_uq_dict[gid] = True
            if "[" in gene_product_subject:
                gene_product = gene_product_subject.split("[")[0].strip()
                species = gene_product_subject.split(
                    "[")[1].replace("]", "").strip()
            else:
                gene_product = gene_product_subject
                species = "NA"
            contig_name, gene_start, gene_end, strand = gene_dict[gid]
            blast_coords = get_coords(
                gene_start,
                gene_end,
                qstart,
                qend,
                strand)
            out_gff.write(
                '%s\tBLAST\tBLAST_gp\t%s\t%s\t%s\t%s\t.\tgene_product=%s;species=%s;sseqid=%s;parent=%s\n' %
                (contig_name,
                 blast_coords[0],
                 blast_coords[1],
                 evalue,
                 strand,
                 gene_product,
                 species.replace,
                 sseqid,
                 gid,
                 ))
    f.close()
    out_gff.close()


# extract the information from the hmmer hit file which goes into the
# output gff file and extract the original position in the genome/contig
def run(
        in_file,
        out_gff_f,
        gene_gff_f,
        anno_mode,
        mode,
        anno_type=None,
        rel_feats_f=None):
    gene_dict = read_gff(gene_gff_f, mode)
    if anno_mode == "B":
        write_blast_gff(in_file, out_gff_f, gene_dict)
    if anno_mode == "H":
        write_hmm_gff(
            in_file,
            out_gff_f,
            gene_dict,
            False,
            mode,
            anno_type,
            rel_feats_f)
    if anno_mode == "F":
        # TODO account for ncbi gff
        figfam2gff.write_fr_ss_gff(in_file, out_gff_f, gene_dict, mode)
    if anno_mode == "M":
        # TODO account for ncbi gff
        write_motif_gff(in_file, out_gff_f, gene_dict)

if __name__ == "__main__":
    # testing
    # uncomment and change parameter if clause to run
    # TODO include test case for figfams
    #run("testing/MO_Vul009-contig00001.domtblout.dat.filtered.best", "testing/MO_Vul009-contig00001.pfam.gff","testing/MO_Vul009-contig00001.mg.gff", anno_mode="H",mode="gm", anno_type="Pfam")
    #run("testing/MO_Vul009-contig00001.prosite.gff", "testing/MO_Vul009-contig00001.contig.prosite.gff","testing/MO_Vul009-contig00001.mg.gff", anno_mode="M",mode="gm")
    #run("testing/MO_Vul009-contig00001.blast.txt", "testing/MO_Vul009-contig00001.contig.blast.gff","testing/MO_Vul009-contig00001.mg.gff", anno_mode="B",mode="gm")
    #run("testing/NC_003197_domtblout.dat.filtered.best", "testing/NC_003197.contig.pfam.gff", "testing/NC_003197.gff", anno_mode="H", mode="ncbi", anno_type="Pfam")
    import getopt
    if len(sys.argv) < 8:
    # if len(sys.argv) < 1:
        print """ wrong number of parameters
Usage:
-i <input annotation file, either blast or hmm>
-o <output gff annotation file>
-g <gff gene prediction file>
-m ncbi or gm for NCBI genomes gff file or genemark gff file
-b H for HMM B for blast or F for figfam
-a name of the annotation e.g. Pfam
-r optionally a name of a file containing the relevant features
e.g.: python /net/metagenomics/projects/cellulose_degraders/code/annotation/hmm2gff.py
-i  /net/metagenomics/data/genomes/annotation/NCBI_20140115/pfam_27.0_hmmer_3.0/RefSeqProId+88073/NC_016894_domtblout.dat.filtered.best
-g /net/metagenomics/data/genomes/public/NCBI20140520_GFFs/Acetobacterium_woodii_DSM_1030_uid88073/NC_016894.gff
-m ncbi -o /tmp/NC_016894_pfam.gff
-b H -a Pfam
"""
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "i:o:g:b:a:r:m:")
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(2)
    in_file = None
    out_gff_f = None
    gene_gff_f = None
    anno_type = None
    anno_mode = None
    rel_feats_f = None
    mode = None
    for o, a in optlist:
        if o == "-i":
            in_file = a
        if o == "-o":
            out_gff_f = a
        if o == "-g":
            gene_gff_f = a
        if o == "-b":
            anno_mode = a
        if o == "-a":
            anno_type = a
        if o == "-r":
            rel_feats_f = a
        if o == "-m":
            mode = a
    if anno_mode == "B":
        run(in_file, out_gff_f, gene_gff_f, anno_mode, mode)
    if anno_mode == "H":
        run(in_file,
            out_gff_f,
            gene_gff_f,
            anno_mode,
            mode,
            anno_type=anno_type,
            rel_feats_f=rel_feats_f)
    if anno_mode == "F":
        run(in_file, out_gff_f, gene_gff_f, anno_mode)
    if anno_mode == "M":
        run(in_file, out_gff_f, gene_gff_f, anno_mode)
