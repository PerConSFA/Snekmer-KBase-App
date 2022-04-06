# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import yaml
import shutil
import subprocess
import sys
import uuid
from pprint import pformat
from Bio import SeqIO

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_uploadmethodsClient import kb_uploadmethods
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

#END_HEADER


class Snekmer:
    '''
    Module Name:
    Snekmer

    Module Description:
    A KBase module: Snekmer
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/abbyjerger/Snekmer.git"
    GIT_COMMIT_HASH = "3e1dab70b6bbcfafd4f8ef93f3a83b829a3d03f2"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_Snekmer_model(self, ctx, params):
        """
        run_Snekmer_model accepts some of the model params for now, and returns results in a KBaseReport
        :param params: instance of type "SnekmerModelParams" (Input
           parameters for Snekmer Model. workspace_name - the name of the
           workspace for input/output kmer - kmer length for features
           alphabet - mapping function for reduced amino acid sequences
           min_rep_thresh - min number of sequences to include feature for
           prefiltering processes - for parallelization) -> structure:
           parameter "workspace_name" of String, parameter "kmer" of Long,
           parameter "alphabet" of String, parameter "min_rep_thresh" of
           Double, parameter "processes" of Long
        :returns: instance of type "SnekmerModelOutput" (Output parameters
           for Snekmer Model. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Snekmer_model

        #print('running snekmer model with params=')
        #pprint(params)

        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['kmer']},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }

        #END run_Snekmer_model

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_Snekmer_model return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_Snekmer_search(self, ctx, params):
        """
        run_Snekmer_search accepts some of the search params for now, and returns results in a KBaseReport
        :param params: instance of type "SnekmerSearchParams" (Input
           parameters for Snekmer Search. staging_file_subdir_path - file in
           staging area workspace_name - the name of the workspace for
           input/output object_ref - Genome object with Protein Translation
           sequence in the Feature k - kmer length for features alphabet -
           mapping function for reduced amino acid sequences min_rep_thresh -
           min number of sequences to include feature for prefiltering
           processes - for parallelization) -> structure: parameter
           "staging_file_subdir_path" of String, parameter "workspace_name"
           of String, parameter "object_ref" of String, parameter "k" of
           Long, parameter "alphabet" of Long, parameter "min_rep_thresh" of
           Long, parameter "processes" of Long
        :returns: instance of type "SnekmerSearchOutput" (Output parameters
           for Snekmer Search. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Snekmer_search

        # check inputs
        staging_file = params['staging_file_subdir_path']
        object_ref = params['object_ref']
        workspace_name = params['workspace_name']
        if 'k' not in params:
            raise ValueError('Parameter kmer is not set in input arguments')
        k = params['k']
        if 'processes' not in params:
            raise ValueError('Parameter processes is not set in input arguments')
        processes = params['processes']
        if 'alphabet' not in params:
            raise ValueError('Parameter alphabet is not set in input arguments')
        alphabet = params['alphabet']
        if 'min_rep_thresh' not in params:
            raise ValueError('Parameter min_rep_thresh is not set in input arguments')
        min_rep_thresh = params['min_rep_thresh']

        # add these param inputs from the UI to the config.yaml
        print('Save UI inputs into config.yaml')
        print("=" * 80)
        new_params = {'k': k, 'alphabet': alphabet,
                      'min_rep_thresh': min_rep_thresh,
                      'processes': processes}
        with open('/kb/module/data/config.yaml', 'r') as file:
            my_config = yaml.safe_load(file)
            my_config.update(new_params)

        # save updated config.yaml to self.shared_folder and
        # start creating directory structure Snekmer needs
        with open(f"{self.shared_folder}/config.yaml", 'w') as file:
            yaml.safe_dump(my_config, file)
        os.makedirs(f"{self.shared_folder}/input")

        # save model_outputs from data to /kb/module/work/tmp
        shutil.copytree("/kb/module/data/model_output", f"{self.shared_folder}/model_output")
        print("="*80)

        # Use input Genome to produce a FASTA file with the protein sequences of the CDSs
        print('Downloading Genome input as protein FASTA file.')
        print("=" * 80)
        genomeUtil = GenomeFileUtil(self.callback_url)
        fasta_file = genomeUtil.genome_proteins_to_fasta({'genome_ref': object_ref})
        #, 'include_functions': 0,'include_aliases': 0
        #sys.stdout.flush()

        # save fasta to the input folder
        shutil.copyfile(fasta_file['file_path'], f"{self.shared_folder}/input/inputfromgenome.fasta")
        print("="*80)

        # after self.shared_folder directory is set up, run commandline section
        print('Run subprocess of snekmer search')
        print("=" * 80)
        cmd_string = "snekmer search"
        cmd_process = subprocess.Popen(cmd_string, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, cwd=self.shared_folder,
                                       shell=True)
        cmd_process.wait()
        print('return code: ' + str(cmd_process.returncode))
        print("="*80)
        output, errors = cmd_process.communicate()
        print("output: " + str(output) + '\n')
        print("errors: " + str(errors) + '\n')
        print("="*80)

        # set up output directory for output files
        output_directory = os.path.join(self.shared_folder, "output/search")
        #os.makedirs(output_directory, exist_ok=True)
        #result_file = os.path.join(output_directory, 'search.zip')
        print("output_directory: " + output_directory)
        print("=" * 80)
        #print("result_file: " + result_file)

        #dfu = DataFileUtil(self.callback_url)
        #report_shock_id = dfu.file_to_shock({'file_path': output_directory,
                                             #'pack': 'zip'})

        # f"{self.shared_folder}/output/search"
        # Step - Build a Report and return
        print('Section: build report data.')
        output_files = list()
        output_files.append({
            'path': output_directory,
            'name': os.path.basename(output_directory),
            'label': os.path.basename(output_directory),
            'description': "Files generated by Snekmer Search"})

        report_params = {
            'message': 'Kmer input was ' + str(k),
            'workspace_name': workspace_name,
            'objects_created': [],
            'file_links': output_files
        }

        report_client = KBaseReport(self.callback_url)
        output = report_client.create_extended_report(report_params)

        # STEP 6: construct the output to send back
        report_output = {'report_name': output['name'], 'report_ref': output['ref']}
        #END run_Snekmer_search

        # At some point might do deeper type checking...
        if not isinstance(report_output, dict):
            raise ValueError('Method run_Snekmer_search return value ' +
                             'output is not type dict as required.')
        # return the results
        return [report_output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
