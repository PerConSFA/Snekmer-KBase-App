# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import yaml
import shutil
import subprocess
import zipfile
import sys
import uuid
from pprint import pformat
from pprint import pprint
from Bio import SeqIO
from datetime import datetime
from pathlib import Path
import pandas as pd

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_uploadmethodsClient import kb_uploadmethods
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils

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
    GIT_COMMIT_HASH = "34a4d55cd05537d6b09cb9e6cc01d2a740413aea"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.workspaceURL = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.DOTFU = KBaseDataObjectToFileUtils(self.callback_url)
        self.wsClient = workspaceService(self.workspaceURL)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_Snekmer_model(self, ctx, params):
        """
        run_Snekmer_model accepts some of the model params for now, and returns results in a KBaseReport
        :param params: instance of type "SnekmerModelParams" (Input
           parameters for Snekmer Model. workspace_name - the name of the
           workspace for input/output k - kmer length for features alphabet -
           mapping function for reduced amino acid sequences min_rep_thresh -
           min number of sequences to include feature for prefiltering
           processes - for parallelization) -> structure: parameter
           "workspace_name" of String, parameter "k" of Long, parameter
           "alphabet" of String, parameter "min_rep_thresh" of Double,
           parameter "processes" of Long
        :returns: instance of type "SnekmerModelOutput" (Output parameters
           for Snekmer Model. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Snekmer_model

        # check inputs
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
        :param params: instance of type "SnekmerSearchParams" -> structure:
           parameter "workspace_name" of String, parameter "object_ref" of
           String, parameter "k" of Long, parameter "alphabet" of Long
        :returns: instance of type "SnekmerSearchOutput" (Output parameters
           for Snekmer Search. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Snekmer_search
        # Print statements to stdout/stderr are captured and available as the App log
        logging.info('Starting run_Snekmer_search function. Params=' + pformat(params))
        # Check inputs
        logging.info('Validating parameters.')
        if 'object_ref' not in params:
            raise ValueError('Parameter object_ref is not set in input arguments')
        object_ref = params['object_ref']
        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'k' not in params:
            raise ValueError('Parameter k is not set in input arguments')
        k = params['k']
        if 'alphabet' not in params:
            raise ValueError('Parameter alphabet is not set in input arguments')
        alphabet = params['alphabet']

        # testing diff between dfu.get_objects and wsClient.get_objects2 (which is used in GenomeSetToFasta)
        obj_dfu_get_obj = self.dfu.get_objects({'object_refs': [object_ref]})
        print("Using DataFileUtil.get_objects: ", obj_dfu_get_obj)

        obj_ws_get_obj2 = self.wsClient.get_objects2({'objects': [{'ref': object_ref}]})
        print("Using WorkspaceClient.get_objects2: ", obj_ws_get_obj2)

        # get GenomeSet name to add into protein fasta file names
        data_obj = self.dfu.get_objects({'object_refs': [object_ref]})['data'][0]
        info = data_obj['info']
        obj_name = str(info[1])

        print('data_obj data for the GenomeSet:', data_obj)
        GenomeSetToFASTA_params = {
            'genomeSet_ref': object_ref,
            'file': obj_name,
            'residue_type': 'protein',
            'feature_type': 'CDS',
            'record_id_pattern': '%%feature_id%%',
            'merge_fasta_files': 'FALSE'
        }
        print("GenomeSetToFasta params: ")
        pprint(GenomeSetToFASTA_params)

        print("Calling GenomeSetToFasta: ")
        GenomeSetToFASTA_retVal = self.DOTFU.GenomeSetToFASTA(GenomeSetToFASTA_params)
        fasta_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list']
        print("=" * 80)
        print("Fasta file path: ")
        print(fasta_file_path)

        # need to rename fasta files since GenomeSetToFasta doesn't use the Genome object's sciname
        genome_ref_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
        print("Genome references should be turned into their scientific name: ")
        print(genome_ref_sci_name)

        genomeSet_object = self.wsClient.get_objects2({'objects': [{'ref': object_ref}]})['data'][0]['data']
        #genomeSet_object['elements'].keys()
        print("=" * 80)
        print("genomeSet_object: ", genomeSet_object)
        # Add params from the UI to the config.yaml
        logging.info('Writing UI inputs into the config.yaml')
        new_params = {'k': k, 'alphabet': alphabet}
        with open('/kb/module/data/config.yaml', 'r') as file:
            my_config = yaml.safe_load(file)
            my_config.update(new_params)

        # save updated config.yaml to self.shared_folder and
        # start creating directory structure Snekmer needs
        with open(f"{self.shared_folder}/config.yaml", 'w') as file:
            yaml.safe_dump(my_config, file)
        os.makedirs(f"{self.shared_folder}/input")

        # save model_outputs from data to /kb/module/work/tmp
        #shutil.copytree("/kb/module/data/model_output", f"{self.shared_folder}/model_output")
        # faster testing
        shutil.copytree("/kb/module/data/small_test_model_output", f"{self.shared_folder}/small_test_model_output")
        print("="*80)
        print("Next copy protein fastas from /kb/module/work/tmp to /kb/module/work/tmp/input")

        # save each protein FASTA to the input folder
        for i in range(len(fasta_file_path)):
            shutil.copy(fasta_file_path[i], f"{self.shared_folder}/input")
        print("="*80)
        print("Copied protein fastas to the input folder for the subprocess step")

        # now that the protein files are in /input
        # remove .params, replace with sci name, and then add .faa extension
        mypath = Path(f"{self.shared_folder}/input")
        for file, i in zip(os.listdir(mypath), genome_ref_sci_name.values()):
            print("Filename in beginning of loop: ", file)
            new_file = file.split('.', 1)[0]
            print("Filename when everything after first period is removed: ", new_file)
            print("genome_ref_sci_name: ", i)
            new_file += "."
            new_file += i
            print("Filename after adding sci name: ", new_file)
            new_file += ".faa"
            print("Filename after adding .faa extension", new_file)
            old = os.path.join(mypath, file)
            new = os.path.join(mypath, new_file)
            if not os.path.exists(new):  # check if the file doesn't exist
                os.rename(old, new)
            print("=" * 80)

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
        result_directory = os.path.join(self.shared_folder, "output", "search", "")
        output_files = list()
        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        os.makedirs(output_directory)
        run_date = datetime.now().strftime("%Y.%m.%d-%I:%M:%S%p")
        result_name = "SnekmerSearch" + str(k) + str(alphabet) + str(run_date) + ".zip"
        result_file = os.path.join(output_directory, result_name)

        print("result directory: " + result_directory)
        print("=" * 80)
        print("output_directory: " + output_directory)
        print("=" * 80)
        print("result_file: " + result_file)
        print("=" * 80)

        # zip output files for the KBase report
        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(result_directory):
                for file in files:
                    if file.endswith('.csv'):
                        zip_file.write(os.path.join(root, file),
                                       os.path.join(os.path.basename(root), file))

        output_files.append({
            'path': result_file,
            'name': os.path.basename(result_file),
            'label': os.path.basename(result_file),
            'description': 'Files generated by Snekmer Search'})

        # analyze csv outputs for the KBase report
        # combine the Search csv outputs into one csv

        print("Starting to analyze the csvs: ")
        filelist = []
        for root, dirs, files in os.walk(result_directory):
            for file in files:
                if file.endswith(".csv"):
                    filelist.append(os.path.join(root, file))

        combined_csv = pd.concat([pd.read_csv(f) for f in filelist])
        combined_csv.to_csv("combined_csv.csv", index=False, encoding='utf-8-sig')

        unique_seq = len(pd.unique(combined_csv['sequence_id']))
        total_seq = len(combined_csv.index)
        TF_counts = combined_csv['in_family'].value_counts()
        print()
        print(TF_counts)

        # prep params for report
        print("=" * 80)
        print("Prep params for report: ")

        genomes_run = list(data_obj['data']['elements'].keys())
        report_message = "Kmer input: {0}\n" \
                         "Alphabet: {1}\n" \
                         "Genomes run: {2}\n" \
                         "Unique number of sequences: {3}\n" \
                         "Total number of sequences: {4}\n" \
                         "Sequences in a family: \n{5}".format(str(k), alphabet, genomes_run,
                                                               unique_seq, total_seq, TF_counts)
        print("Report message:")
        print(report_message)

        # Build a Report and return
        print('Section: build report data.')
        report_params = {
            'message': report_message,
            'workspace_name': workspace_name,
            'objects_created': [],
            'file_links': output_files
        }

        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)

        # construct the output to send back
        output = {'report_name': report_info['name'], 'report_ref': report_info['ref']}
        #END run_Snekmer_search

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_Snekmer_search return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
