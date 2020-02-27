#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import sys
import re
import json
import shutil
import logging
import tarfile
import argparse
import subprocess
import multiprocessing as mp
from functools import partial
from distutils.spawn import find_executable
from pkg_resources import resource_filename
## 3rd party
import pandas as ps
## application
from traitar import hmm2gff
from traitar import utils
from traitar.PhenotypeCollection import PhenotypeCollection
import traitar.domtblout2gene_generic as domtblout2gene_generic
import traitar.hmmer2filtered_best as hmmer2filtered_best
import traitar.merge_preds as merge_preds
import traitar.predict as predict

# package data
primary_default_models = resource_filename('traitar',
                                           os.path.join('data', 'models', 'phypat.tar.gz'))
secondary_default_models = resource_filename('traitar',
                                           os.path.join('data', 'models', 'phypat+PGL.tar.gz'))


# class init
class Traitar:
    def __init__(self, input_dir, output_dir, sample2file,
                 cpu = 1, parallel = 1,
                 heatmap_out = None, heatmap_format = "pdf",
                 no_heatmap_phenotype_clustering = False,
                 no_heatmap_sample_clustering = False,
                 gene_gff_type = None, pfam_dir = None,
                 primary_models = None, secondary_models = None, overwrite = False):
        #self.user_message = "output dir %s already exists; press 1 to continue with data from a previous run; press 2 to remove this directory; press 3 to abort followed by [ENTER]"
        #self.error_message =  "directory %s already exists; delete directory or run in interactive mode if this step is done and you want to continue from there"
        self.overwrite = overwrite
        self.heatmap_format = heatmap_format        
        self.no_heatmap_phenotype_clustering = no_heatmap_phenotype_clustering
        self.no_heatmap_sample_clustering = no_heatmap_sample_clustering
        self.sample2file = sample2file
        self.input_dir = input_dir
        self.gene_gff_type = gene_gff_type
        self.s2f = self.parse_sample_f()
        self.cpu = cpu
        self.parallel = parallel
        self.output_dir = output_dir
        #self.phenolyzer_dir = os.path.abspath(os.path.dirname(traitar.__file__)) 
        #determines whether the primary and secondary phenotype models specify the same set of phenotypes and can be combined in one heatmap
        #initially set to false
        if primary_models is not None:
            self.primary_models = PhenotypeCollection(primary_models)
        else:
            self.primary_models = PhenotypeCollection(primary_default_models)
        if secondary_models is not None and primary_models is not None:
            self.secondary_models = PhenotypeCollection(secondary_models)
        elif self.primary_models.get_name() != "phypat":
            self.secondary_models = None
        else:
            self.secondary_models = PhenotypeCollection(secondary_default_models)
        self.is_gnu_parallel_available = self.is_exe("parallel")
        self.heatmap_out = heatmap_out
        if cpu != 1 and not self.is_gnu_parallel_available:
            msg = 'GNU parallel is not available on the command line'
            msg += ';make sure you installed it properly or decrease number of cpus to 1'
            sys.stderr.write(msg + '\n')
            sys.exit(1)
        # 'config' (now just {hmms : path})
        self.config = {'hmms' : os.path.abspath(pfam_dir)}
        #create output dir
        self.check_dir(output_dir)
        #pred dir
        self.pred_dir = os.path.join(self.output_dir, "phenotype_prediction")
        self.phypat_dir = os.path.join(self.pred_dir, self.primary_models.get_name())
        if not self.secondary_models is None:
            self.phypat_pgl_dir = os.path.join(self.pred_dir, self.secondary_models.get_name())

    def _special_match(self, strg, search = re.compile(r'[^A-Za-z0-9.\-_]').search):
        return not bool(search(strg))
   
    def parse_sample_f(self):
        """read sample file with pandas and make sure the content is appropriate"""
        #check if sample files exist
        if not os.path.exists(self.sample2file):
            sys.exit("sample file %s does not exist" % self.sample2file)
        s2f = ps.read_csv(self.sample2file, dtype = 'string', sep = "\t")
        col_check_cols = ["sample_file_name", "sample_name", "category", "gene_gff"]
        col_check = dict((i, False) for i in col_check_cols)
        for i in s2f.columns:
            if i not in ["sample_file_name", "sample_name", "category", "gene_gff"]:
                sys.exit('{} is not a valid column identifier'.format(i))
            col_check[i] = True        
        for x in ['sample_file_name', 'sample_name']:
            if x not in s2f.columns:
                sys.exit('"{}" column in input sample_file missing'.format(x))
        for i in s2f.loc[:, "sample_file_name"]:
            if not os.path.exists(os.path.join(self.input_dir,i)):
                pass
        for i in s2f.loc[:, 'sample_name']:
            if not self._special_match(i):
                msg = 'Knvalid character in sample name %s; only [a-zA-Z0-9.-_] allowed'
                sys.exit(msg.format(i))
            if len(i) > 41:
                msg = 'sample name too long: {}; names may only have 40 characters'
                sys.ext(msg.format(i))
        if col_check['category']:
            uq = s2f.loc[:, 'category'].unique()
            max_cats = 12
            if len(uq) > max_cats:
                msg = 'reduce the number of sample categories to less than {}'
                sys.exit(msg.format(max_cats))
            for i in uq:
                msg = 'Category error: {}'
                msg += ';sample categories may not be longer than 30 characters'
                if len(i) > 30:
                    sys.exit(msg.format(i))
        if col_check['gene_gff']:
            for i in s2f.loc[:, 'gene_gff']:
                if not os.path.exists(os.path.join(self.input_dir,i)):
                    msg = 'WARNING: sample file {} does not exist in the output directory: {}'
                    msg += '; The program will fail if not in "from_summary_annotation" mode\n'
                    sys.stderr.write(msg.format(self.input_dir, i))
            if self.gene_gff_type is None:
                msg = 'gene gff type needs to be specified with -g / --gene_gff_type'
                msg += '<gene_gff_type> if sample file contains gene_gff column'
                sys.exit(msg)
        return s2f

    def is_exe(self, program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return True 
        else:
            for path in os.environ['PATH'].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return True 
                else:
                    False 

    def check_dir(self, out_dir):
        """ Return True if not present or self.overwrite """
        if os.path.exists(out_dir):
            if self.overwrite:
                logging.info('  `--overwrite` used; deleting dir: {}'.format(out_dir))
                shutil.rmtree(out_dir)
                return True
            else:
                msg = 'Using any existing files in the output directory: {}'
                logging.info(msg.format(out_dir))
                logging.info('  Use `--overwrite` to delete the dir & start fresh')
                return False
        else:
            os.makedirs(out_dir)
            return True

    def file_overwrite(self, infile):
        """Return True if infile doesn't exist or should be overwritten, else False"""
        create_file = True
        msg = ''
        if os.path.isfile(infile):
            if self.overwrite is True:
                msg = '  File exists, but re-creating due to --overwrite: {}'
                create_file = True
            else:
                msg = '  Using existing file (use --overwrite to recreate): {}'
                create_file = False
        else:
            msg = '  File doesn\'t exist, so creating: {}'
            create_file = True
        logging.info(msg.format(infile))
        return create_file
        
    def annotate(self, mode):
        pfam_msg = 'Running annotation with hmmer. This step can take a while.'
        pfam_msg += ' A rough estimate for sequential Pfam annotation of genome'
        pfam_msg += ' samples of ~3 Mbs is 10 min per genome.'
        install_msg = '{} not available on command line'
        install_msg += '; please make sure you have properly installed it'
        #check if executables are available 
        if not find_executable('hmmsearch'):
            sys.stderr.write(install_msg.format('hmmsearch') + '\n')
            sys.exit(1)
        if mode == 'from_nucleotides':
            if not find_executable('prodigal'):
                sys.stderr.write(install_msg.format('prodigal') + '\n')
                sys.exit(1)
            logging.info('Running gene prediction with Prodigal')
            sys.stdout.flush()
            self.run_gene_prediction(self.s2f.loc[:,'sample_file_name'],
                                     self.s2f.loc[:,'sample_name'])
            logging.info(pfam_msg)
            sys.stdout.flush()
            self.run_hmmer_annotation(self.s2f.loc[:,"sample_name"],
                                      self.s2f.loc[:,"sample_name"], mode)
        if mode == 'from_genes': 
            logging.info(pfam_msg)
            sys.stdout.flush()
            self.run_hmmer_annotation(self.s2f.loc[:,'sample_file_name'],
                                      self.s2f.loc[:,'sample_name'], mode)

    def phenolyze(self, mode):
        logging.info('Running phenotype prediction')
        sys.stdout.flush()
        #is_recompute = self.check_dir(self.pred_dir)
        pred_files = self.run_phenotype_prediction()
        logging.info('Running feature track generation')
        if not 'gene_gff' in self.s2f.columns:
            self.run_feature_track_generation(pred_files, self.s2f.loc[:,'sample_name'], mode)
        else:
            self.run_feature_track_generation(pred_files, self.s2f.loc[:,'sample_name'],
                                              mode, self.s2f.loc[: , 'gene_gff'],
                                              self.gene_gff_type)
        sys.stdout.flush()
        logging.info('Running heatmap generation')
        self.run_heatmap_generation(self.s2f.loc[:,'sample_name'])
        sys.stdout.flush()
                        
    def run_gene_prediction(self, in_samples, out_samples):
        #create output directory for the gene prediction 
        gp_dir = os.path.join(self.output_dir, 'gene_prediction')
        is_recompute = self.check_dir(gp_dir)
        prodigal = 'prodigal -a {gp_dir}/{out_sample}.faa  -f gff'
        prodigal += ' < {in_dir}/{in_sample} > {gp_dir}/{out_sample}.gff '
        prodigal_commands = []
        for i in range(len(in_samples)):
            prodigal_commands.append(prodigal.format(in_dir = self.input_dir,
                                                     in_sample = in_samples[i],
                                                     out_sample = out_samples[i],
                                                     gp_dir  = gp_dir))
        if is_recompute:    
            utils.execute_commands(prodigal_commands,
                                   parallel = self.parallel,
                                   joblog = os.path.join(gp_dir, 'joblog.txt'))

    def run_hmmer_annotation(self, in_samples, out_samples, mode):
        file_extension = False
        in_dir = self.input_dir
        if mode == 'from_nucleotides':
            file_extension = True 
            in_dir = os.path.join(self.output_dir, 'gene_prediction')
        #create output directory for the pfam annotation 
        a_dir_base = os.path.join(self.output_dir, 'annotation')
        # check if output directory already exists; recompute if not
        is_recompute = self.check_dir(a_dir_base)
        # models
        if self.secondary_models is not None and \
           self.primary_models.get_hmm_name() != self.secondary_models.get_hmm_name():
            models = [self.primary_models, self.secondary_models]
        else:
            models = [self.primary_models]
        # hmm
        logging.info('Computing hmms')
        for pt_models in models:
            a_dir = os.path.join(a_dir_base, pt_models.get_hmm_name())
            if not os.path.isdir(a_dir):
                os.makedirs(a_dir)
            # hmmsearch
            hmmer = 'hmmsearch --notextw --cpu {cpu} --cut_ga --domtblout {outfile}'
            hmmer += ' {hmms} {infile}'
            hmmer_commands = []
            for i in range(len(in_samples)):
                f_ext = '.faa' if file_extension else ""
                hmms_path = os.path.join(self.config['hmms'], pt_models.get_hmm_f())
                outfile = os.path.join(a_dir, out_samples[i] + '_domtblout.dat')
                infile = os.path.join(in_dir, in_samples[i] + f_ext)
                cmd = hmmer.format(infile = infile,
                                   outfile = outfile,
                                   hmms = hmms_path,
                                   cpu = self.cpu)
                if self.file_overwrite(infile) is True:
                    hmmer_commands.append(cmd)

            x = os.path.join(a_dir, 'joblog_hmmer.txt')
            logging.info('  No. of hmmer commands to run: {}'.format(len(hmmer_commands)))
            utils.execute_commands(hmmer_commands,
                                   parallel = self.parallel,
                                   joblog = x)

            # hmmer2filtered_best
            logging.info('Filtering hmm hits')
            params = []
            dat_best_files = []
            for i in range(len(in_samples)):
                infile = os.path.join(a_dir, out_samples[i] + '_domtblout.dat')
                outfile = os.path.join(a_dir, out_samples[i] + '_filtered_best.dat')
                dat_best_files.append(outfile)
                if self.file_overwrite(outfile):
                    params.append([infile, outfile])                    
            fae_func = partial(hmmer2filtered_best._main_par,
                               hmm_name = pt_models.get_hmm_name())
            logging.info('  No. of hmmer filter commands to run: {}'.format(len(params)))
            if self.parallel > 1:
                pool = mp.Pool(self.parallel)
                pool.map(fae_func, params)
                pool.close()
            else:
                for p in params:
                    fae_func(p)
                
            # summary matrix from all best hmms
            logging.info('Creating a summary matrix')
            dat_file = os.path.join(a_dir, 'summary.dat')                
            domtblout2gene_generic.main(outfile = dat_file,
                                        in_filtered_best_fs = dat_best_files,
                                        archive_f = pt_models.get_archive_f())        

    def run_phenotype_prediction(self):
        logging.info('Predicting with primary models')
        pred_files = {}
        dat_file = '{out_dir}/annotation/{hmm_f}/summary.dat'
        dat_file = dat_file.format(out_dir =  self.output_dir,
                                   hmm_f =  self.primary_models.get_hmm_name())
        pred_files['primary'] = predict.main(pt_models = self.primary_models.get_archive_f(),
                                             annotation_matrix = dat_file, 
                                             out_dir = self.output_dir,
                                             voters = 5,
                                             cpus = self.cpu)
        if not self.secondary_models is None:
            logging.info('Predicting with secondary models')
            dat_file = '{out_dir}/annotation/{hmm_f}/summary.dat'
            dat_file = dat_file.format(out_dir =  self.output_dir,
                                       hmm_f =  self.secondary_models.get_hmm_name())
            x = predict.main(pt_models = self.secondary_models.get_archive_f(),
                             annotation_matrix = dat_file,
                             out_dir =  self.output_dir,
                             voters = 5,
                             cpus = self.cpu)
            pred_files['secondary'] = x
            #combine phypat and phypat+PGL predictions
            logging.info('Merging primary & secondary predictions')
            outdir = os.path.join(self.output_dir, 'phenotype_prediction')
            merge_preds.main(pred_files = pred_files,
                             primary_name = self.primary_models.get_name(),
                             secondary_name = self.secondary_models.get_name(),
                             out_dir = outdir,
                             voters = 5)
        return pred_files
    
    def run_feature_track_generation(self, pred_files, in_samples, mode, gene_gffs = None,
                                     gene_gff_type = 'prodigal'):
        """map the phenotype relevant protein families and write mapping to disk"""
        logging.warning('hmm2gff currently not supported. SKIPPING!')
        return 0
        #read in phypat predictions
        phypat_preds = ps.read_csv(pred_files['primary']['majority-vote'],
                                   index_col = 0, sep = '\t', encoding = 'utf-8')
        phypat_preds.index = phypat_preds.index.values.astype(str)
        #read pfam phenotype id mapping file from model tar
        pt2desc_phypat = self.primary_models.get_acc2pt()
        #secondary models
        if not self.secondary_models is None:
            infile = pred_files['secondary']['majority-vote']
            phypat_pgl_preds = ps.read_csv(infile, index_col = 0,
                                           sep = '\t', encoding = 'utf-8') 
            phypat_pgl_preds.index = phypat_preds.index.values.astype(str)
            pt2desc_phypat_pgl = self.secondary_models.get_acc2pt() 
        #collect predictions and compose hmm2gff command for each samples
        hmm2gff = 'hmm2gff.py {outfile} {out_gff_dir} {sample} {model_tar} {predicted_pts}' 
        gene_gff_extd = ' --gene_gff {gene_gff} --gene_gff_type ' + gene_gff_type
        h2gff_commands = []        
        for i in range(len(in_samples)):
            predicted_pts_phypat = []
            predicted_pts_phypat_pgl = []
            for j in phypat_preds.columns:
                if phypat_preds.loc[in_samples[i], j] == 1:
                    predicted_pts_phypat.append(str(pt2desc_phypat.loc[j,][0]))
            if not self.secondary_models is None:
                for j in phypat_pgl_preds.columns:
                    if not self.secondary_models is None and \
                       phypat_pgl_preds.loc[in_samples[i], j] == 1:
                        predicted_pts_phypat_pgl.append(str(pt2desc_phypat_pgl.loc[j,][0]))
            # phypat
            if not len(predicted_pts_phypat) == 0:
                outfile = os.path.join(self.output_dir, 'annotation',
                                       self.primary_models.get_hmm_name(),
                                       '{}_filtered_best.dat'.format(in_samples[i]))
                cmd = hmm2gff.format(out_gff_dir = '{}/feat_gffs/'.format(self.phypat_dir),
                                     outfile = outfile,
                                     sample = in_samples[i],
                                     model_tar = self.primary_models.get_archive_f(),
                                     predicted_pts = ','.join(predicted_pts_phypat))
                if gene_gffs is not None:
                    x = os.path.join(self.input_dir,  gene_gffs[i])                    
                    cmd += gene_gff_extd.format(gene_gff = x)
                if self.file_overwrite(outfile):
                    h2gff_commands.append(cmd)
            # phypat-pgl
            if not len(predicted_pts_phypat_pgl) == 0:
                if not self.secondary_models is None:
                    outfile = os.path.join(self.output_dir, 'annotation',
                                           self.secondary_models.get_hmm_name(),
                                           '{}_filtered_best.dat'.format(in_samples[i]))
                    cmd = hmm2gff.format(out_gff_dir = os.path.join(self.phypat_pgl_dir,
                                                                    'feat_gffs'),
                                         outfile = outfile,
                                         sample = in_samples[i],
                                         model_tar = self.secondary_models.get_archive_f(),
                                         predicted_pts = ','.join(predicted_pts_phypat))
                    if gene_gffs is not None:
                        x = os.path.join(self.input_dir,  gene_gffs[i])                    
                        cmd += gene_gff_extd.format(gene_gff = x)
                    if self.file_overwrite(outfile):
                        h2gff_commands.append(cmd)
        #create output dirs for feature tracks
        logging.info('  No. of hmm2gff commands: {}'.format(len(h2gff_commands)))
        utils.execute_commands(h2gff_commands, parallel=self.parallel) 
        if mode != 'from_nucleotides' and gene_gffs is None:
            ftco = os.path.join(self.pred_dir, 'feature_track_commands.txt')
            msg = 'Tracks with Pfams relevant for the predictions cannot be ad-hoc generated'
            msg += ' because the input is amino acids and no gene prediction GFF files'
            msg += ' have been generated\n commands are saved to "{}" and can be modified'
            msg += ' and run manually'
            logging.info(msg.format(ftco))
            with open(ftco, 'w') as f:
                f.write('\n'.join(h2gff_commands))

    def run_heatmap_generation(self, in_samples, compatible = False):
        """generate a heatmap from the results"""
        logging.warning('Heatmap currently not supported. SKIPPING!')
        return 0
        #TODO introduce nans via the merge prediction routine
        ## mark the nans in the heatmap with a different color
        if self.heatmap_out is None:
            self.heatmap_out = self.pred_dir
        hm_cmd = 'heatmap.py {pred_dir}/{pred_f} {out}/heatmap_{predictor}.{heatmap_format}'
        hm_cmd += ' {secondary_models} --sample_f {sample_file} {model_archive}'
        hm_cmd += ' {phenolyzer}/data/colors.txt {phenotype_clustering} {sample_clustering}'
        hm_dict = {"phenolyzer" : self.phenolyzer_dir,
                   "out": self.heatmap_out,
                   "sample_file" : self.sample2file,
                   "heatmap_format" : self.heatmap_format,
                   "phenotype_clustering" : "--column_method None" \
                   if self.no_heatmap_phenotype_clustering else "",
                   "sample_clustering" : "--row_method None" \
                   if self.no_heatmap_sample_clustering else "" }
        hm_dict_phypat = {"pred_dir" : self.phypat_dir,
                          "pred_f" : "predictions_majority-vote.txt",
                          "predictor": self.primary_models.get_name(),
                          "model_archive" : self.primary_models.get_archive_f(),
                          "secondary_models" : ""}
        cmds = []
        if self.secondary_models is None:
            hm_dict_phypat.update(hm_dict)
            cmds.append(hm_cmd.format(**hm_dict_phypat))
        else:
            hm_dict_phypat_pgl = {"pred_dir" : self.phypat_pgl_dir,
                                  "pred_f" : "predictions_majority-vote.txt",
                                  "predictor": self.secondary_models.get_name(),
                                  "model_archive" : self.secondary_models.get_archive_f(),
                                  "secondary_models" : ""}
            sec_models = "--secondary_model_tar {}".format(self.secondary_models.get_archive_f())
            hm_dict_comb = {"pred_dir" : self.pred_dir,
                            "pred_f" : "predictions_majority-vote_combined.txt",
                            "predictor": "combined",
                            "model_archive" : self.primary_models.get_archive_f(),
                            "secondary_models" : sec_models}
            for i in [hm_dict_phypat, hm_dict_phypat_pgl, hm_dict_comb]:
                i.update(hm_dict)
                cmds.append(hm_cmd.format(**i))
        utils.execute_commands(cmds, parallel=self.parallel)

