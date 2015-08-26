#!/usr/bin/env python
import pandas as ps
import os
import subprocess
import sys
import shutil
import get_external_data
import json
import os
import sys

def phenolyze(args):
    p = Phenolyzer(args.input_dir, args.output_dir, args.sample2file, args.cpus)
    p.run(args.mode)


class Phenolyzer:

    def __init__(self, input_dir, output_dir, sample2file, cpu = 1):
        self.user_message = "output dir %s already exists; press 1 to continue with data from a previous run; press 2 to remove this directory; press 3 to abort"
        self.error_message =  "directory %s already exists; delete directory or run in interactive mode if the annotation is done and you want to continue from there"
        self.s2f = ps.read_csv(sample2file, dtype = 'string', sep = "\t", header = None)
        self.cpu = cpu
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.phenolyzer_dir = "/".join(sys.argv[0].split("/")[:-2])
        with open(os.path.join(self.phenolyzer_dir, "config.json" ), "r") as cf:
            self.config = json.load(cf)
        #create output dir
        self.check_dir(output_dir)


    def check_dir(self, out_dir):
        if os.path.exists(out_dir):
            if os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno()):
                print  self.user_message % out_dir
                user_input = raw_input()
                while user_input not in ["1", "2", "3"]:
                    user_input = raw_input().strip()
                if user_input == "1":
                    return False 
                if user_input == "2":
                    shutil.rmtree(out_dir)
                    os.mkdir(out_dir)
                    return True 
                if user_input == "3":
                    sys.exit(1)

            else:
                sys.stderr.write( self.error_message % out_dir)
        else: 
            os.mkdir(out_dir)
            return True

    def run(self, mode):
        if mode == "from_nucleotides":
            print "running gene prediction with Prodigal"
            self.run_gene_prediction(self.s2f.iloc[:,0], self.s2f.iloc[:,1])
            self.run_hmmer_annotation(self.s2f.iloc[:, 1], self.s2f.iloc[:,1])
        if mode == "from_genes": 
            print "running hmmer prediction"
            self.run_hmmer_annotation(self.s2f.iloc[:,0], self.s2f.iloc[:,1])
        print "running phenotype prediction"
        self.run_phenotype_prediction(self.s2f.iloc[:,1])
         
    
    def execute_commands(self, commands):
        devnull = open('/dev/null', 'w')
        if self.cpu > 1:
            #run with parallel
            ps.DataFrame(commands).to_csv('/tmp/commands.txt', index = False, header = False) 
            subprocess.call("cat /tmp/commands.txt | parallel -j %s" % self.cpu, shell = True, executable = "/bin/bash",  env = env)
        else:
            #run in sequential order
            for i in commands:
                subprocess.call(i,  executable = "/bin/bash", stdout = devnull, shell = True, env = env)

    def run_gene_prediction(self, in_samples, out_samples):
        #create output directory for the gene prediction 
        gp_dir = os.path.join(self.output_dir, "gene_prediction")
        is_recompute = self.check_dir(gp_dir)
        prodigal = "prodigal < %(in_dir)s/%(in_sample)s > %(gp_dir)s/%(out_sample)s.gff  -a %(gp_dir)s/%(out_sample)s.faa  -f gff"
        prodigal_commands = []
        for i in range(len(in_samples)):
            prodigal_commands.append(prodigal % {"in_dir": self.input_dir, "in_sample": in_samples[i], "out_sample":out_samples[i], "gp_dir":gp_dir})
        if is_recompute:    
            self.execute_commands(prodigal_commands) 

    def run_hmmer_annotation(self, in_samples, out_samples):
        if (in_samples != out_samples).all():
            in_dir = self.input_dir
            file_extension = False 
        else:
            file_extension = True
            in_dir = os.path.join(self.output_dir, "gene_prediction")
        #create output directory for the pfam annotation 
        a_dir = os.path.join(self.output_dir, "pfam_annotation")
        #check if output directory already exists and trigger user input if in interactive mode
        is_recompute = self.check_dir(a_dir)
        #run hmmer annotation
        hmmer =  "hmmsearch --cpu 1 --cut_ga  --domtblout %(a_dir)s/%(out_sample)s_domtblout.dat  %(pfam_hmms)s > /dev/null %(in_dir)s/%(in_sample)s%(file_extension)s"      
        hmmer_commands = []
        for i in range(len(in_samples)):
            hmmer_commands.append(hmmer % {"file_extension": ".faa" if file_extension else "", "in_sample":in_samples[i], "out_sample":out_samples[i], "a_dir":a_dir, "phenolyzer":self.phenolyzer_dir, "in_dir" : in_dir, "pfam_hmms" : self.config["pfam_hmms"]})
        #run gff extraction
        #run filtering and best domain hit aggregation
        filter_and_aggregate = "%(phenolyzer)s/code/hmmer2filtered_best.py %(a_dir)s/%(out_sample)s_domtblout.dat   %(a_dir)s/%(out_sample)s_filtered_best.dat 10e-02 25 "
        fae_commands = []
        for i in range(len(in_samples)):
            fae_commands.append(filter_and_aggregate % {"a_dir":a_dir, "in_sample":in_samples[i], "out_sample":out_samples[i], "phenolyzer":self.phenolyzer_dir})
        #run summary matrix computation
        #write temp sample file to disk
        #TODO modify this piece of code so that there is no temp file required anymore
        #best_fs = ps.DataFrame(["%(a_dir)s/%(sample)s_filtered_best.dat"%{"a_dir" : a_dir, "sample":sample} for sample in out_samples])
        #best_fs.to_csv("/tmp/samples_best.txt", index = None, header = None)
        domtblout2gene_generic = "%(phenolyzer)s/code/domtblout2gene_generic.py %(a_dir)s/summary.dat  <(ls %(a_dir)s/*_filtered_best.dat) %(phenolyzer)s/data/sorted_accessions.txt"%{"a_dir": a_dir, "phenolyzer":self.phenolyzer_dir}
        if is_recompute:
            self.execute_commands(hmmer_commands)
            self.execute_commands(fae_commands)
            subprocess.call(domtblout2gene_generic,executable = "/bin/bash", shell = True, env = env)


    def run_phenotype_prediction(self, in_samples):
        #create output directory for the phenotype prediction 
        pred_dir = os.path.join(self.output_dir, "phenotype_prediction")
        is_recompute = self.check_dir(pred_dir)
        phypat_dir = os.path.join(pred_dir, "phypat")
        phypat_ggl_dir = os.path.join(pred_dir, "phypat+GGL")
        if not os.path.exists(phypat_dir):
            os.mkdir(phypat_dir)
        if not os.path.exists(phypat_ggl_dir):
            os.mkdir(phypat_ggl_dir)
        #run phenotype prediction for phypat and phypat+GGL
        predict_phypat = "%(phenolyzer)s/code/predict.py %(phenolyzer)s/data/models/phypat.tar.gz %(pred_dir)s 8476-8568 %(out_dir)s/pfam_annotation/summary.dat -k 5  pfam_pts_names_nl_desc.txt" % {"out_dir" : self.output_dir, "pred_dir" : phypat_dir, "phenolyzer" : self.phenolyzer_dir} 
        predict_phypat_ggl = "%(phenolyzer)s/code/predict.py %(phenolyzer)s/data/models/phypat+GGL.tar.gz %(pred_dir)s 8682-8774 %(out_dir)s/pfam_annotation/summary.dat -k 5  pfam_pts_names_nl_desc.txt" % {"out_dir" : self.output_dir, "pred_dir" : phypat_ggl_dir, "phenolyzer" : self.phenolyzer_dir} 
        if is_recompute:
            subprocess.call(predict_phypat,executable = "/bin/bash", shell = True, env = env)
            subprocess.call(predict_phypat_ggl, executable = "/bin/bash",shell = True, env = env)
        #combine phypat and phypat+GGL predictions
        merge_preds = "%(phenolyzer)s/code/merge_preds.py %(out_dir)s %(phypat_dir)s %(phypat_ggl_dir)s -k 5" %{"out_dir" : os.path.join(self.output_dir, "phenotype_prediction"), "phypat_dir" : phypat_dir, "phypat_ggl_dir" : phypat_ggl_dir, "phenolyzer" : self.phenolyzer_dir } 
        subprocess.call(merge_preds, executable = "/bin/bash",shell = True, env = env)
        #generate a heatmap from the results
        hm_cmd = "%(phenolyzer)s/code/heatmap.py %(pred_dir)s/predictions_aggr_majority-vote_comb.tsv %(pred_dir)s/heatmap_comb.png" %{"phenolyzer" : self.phenolyzer_dir, "pred_dir" : pred_dir}
        hm_cmd_phypat = "%(phenolyzer)s/code/heatmap.py %(phypat_dir)s/predictions_bin_majority-vote.csv %(pred_dir)s/heatmap_phypat.png" %{"phenolyzer" : self.phenolyzer_dir,"phypat_dir" : phypat_dir, "pred_dir" : pred_dir}
        hm_cmd_phypat_ggl = "%(phenolyzer)s/code/heatmap.py %(phypat_ggl_dir)s/predictions_bin_majority-vote.csv %(pred_dir)s/heatmap_phypat_ggl.png" %{"phenolyzer" : self.phenolyzer_dir,"phypat_ggl_dir" : phypat_ggl_dir, "pred_dir" : pred_dir}
        subprocess.call(hm_cmd, executable = "/bin/bash",shell = True, env = env)
        subprocess.call(hm_cmd_phypat, executable = "/bin/bash",shell = True, env = env)
        subprocess.call(hm_cmd_phypat_ggl, executable = "/bin/bash",shell = True, env = env)
    

    def run_feature_track_generation(self):
        #create output directory for the pfam annotation 
        os.mkdir(os.path.join(output_dir, "feature_tracks"))
        pass



if __name__ == "__main__":
    #add phenolyzer dir to the path
    env = os.environ.copy()
    #env["PATH"] = "/net/metagenomics/projects/phenotypes_20130523/code/phenolyzer/code/:" + env["PATH"]
    #env["PATH"] = "/net/programs/Debian-7-x86_64/Prodigal-2.6.1/:" + env["PATH"]
    #env["PATH"] = "/net/programs/Debian-7-all/parallel-20150122/bin/:" + env["PATH"]
    #import package Prodigal 
    #subprocess.call("importpackage parallel 20150122",executable = "/bin/bash", shell = True)
    #import package parallel 
    #subprocess.call("importpackage Prodigal", executable = "/bin/bash", shell = True)
    import argparse
    parser = argparse.ArgumentParser("run phenolyzer program (try: phenolyzer.py -h for help on the sub programs)")
    subparsers = parser.add_subparsers()
    main_p = subparsers.add_parser("phenotype", help = "run annotation and prediction") 
    main_p.add_argument("input_dir", help='directory with the input data')
    main_p.add_argument("sample2file", help='mapping from samples to fasta files')
    main_p.add_argument("mode", help='either from_genes if gene prediction amino acid fasta is available in input_dir otherwise from_nucleotides in this case Prodigal is used to determine the ORFs from the nucleotide fasta files in input_dir', choices = ["from_genes", "from_nucleotides"])
    main_p.add_argument("output_dir", help='directory for the output directory; will be created if it doesn\'t exist yet', default='phenolyzer_output')
    main_p.add_argument("-c", "--cpus", help='number of cpus used for the individual steps; maximum is number of samples; needs parallel', default = 1)
    main_p.set_defaults(func = phenolyze)
    data_p = subparsers.add_parser("config", help = "download neccessary data to run traitAR, mainly Pfam HMMs") 
    data_p.add_argument("download_dest",  help = "download Pfam HMMs into the given download directory and untar and unzip it")
    data_p.add_argument("--local", "-l", action = 'store_true', help = "the Pfam HMMs are in the above directory with name 'Pfam-A.hmm'")
    data_p.set_defaults(func = get_external_data.download)
    args = parser.parse_args()
    args.func(args)


