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
    GIT_COMMIT_HASH = "6130b72f01a06694fd270930b37c2a608dff491a"

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
           String, parameter "k" of Long, parameter "alphabet" of Long,
           parameter "min_rep_thresh" of Long
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
        if 'min_rep_thresh' not in params:
            raise ValueError('Parameter min_rep_thresh is not set in input arguments')
        min_rep_thresh = params['min_rep_thresh']

        # investigate use of GenomeSets
        # can i turn them into protein fasta files?
        print('object_ref from the Impl')
        pprint(object_ref)

        genomeSet_object = self.wsClient.get_objects2({'objects': [{'ref': object_ref}]})['data'][0]['data']
        print("genomeSet_object, which is what GenomeSetToFasta grabs too: ")
        pprint(genomeSet_object)
        print("=" * 80)

        #genome_ids = genomeSet_object['elements'].keys()
        genome_ids = list(genomeSet_object['elements'])
        print("genome_ids now as list, before the for loop: ")
        pprint(genome_ids)
        print("=" * 80)
        print("length of genome_ids: ", len(genome_ids))
        print("range length of genome_ids: ", range(len(genome_ids)))

        for genome_i in range(len(genome_ids)):
            genome_id = genome_ids[genome_i]
            print("genome_id from in the loop: ", genome_id)


        GenomeSetToFASTA_params = {
            'genomeSet_ref': object_ref,

            'residue_type': 'protein',
            'feature_type': 'CDS',

            'merge_fasta_files': 'FALSE'
        }
        print("GenomeSetToFasta params: ")
        pprint(GenomeSetToFASTA_params)


        GenomeSetToFASTA_retVal = self.DOTFU.GenomeSetToFASTA(GenomeSetToFASTA_params)
        fasta_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list']
        print("=" * 80)
        print("Fasta file path: ")
        print(fasta_file_path)

        # need to rename fasta files since GenomeSetToFasta doesn't use the Genome object's sciname
        genome_ref_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
        print("Genome references should be turned into their scientific name: ")
        print(genome_ref_sci_name)

        sys.exit()
        not_string = [x.strip('') for x in fasta_file_path]
        print("These file paths should not be within quotes: ")
        print(not_string)

        new = []
        for i in not_string:
            new[i] = os.path.basename(i)
        print("Fasta file names with extension still: ")
        print(new)

        no_ext = []
        for i in new:
            no_ext[i] = os.path.splitext(i)[0]
        print("Fasta files names without extension: ")
        print(no_ext)

        #sys.exit()
        # Add params from the UI to the config.yaml
        logging.info('Writing UI inputs into the config.yaml')
        new_params = {'k': k, 'alphabet': alphabet,
                      'min_rep_thresh': min_rep_thresh
                      }
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
        print("Next copy protein fastas from /kb/module/work/tmp to /kb/module/work/tmp/input")

        # Use input Genomes to produce FASTA files with the protein sequences of the CDSs
        #print('Downloading Genome inputs as protein FASTA files.')
        #print("=" * 80)
        #genomeUtil = GenomeFileUtil(self.callback_url)
        #fasta_files = []
        #for i in range(len(object_ref)):
            #fasta_files.append(genomeUtil.genome_proteins_to_fasta({'genome_ref': object_ref[i]}))

        # save each protein FASTA to the input folder
        for i in range(len(fasta_file_path)):
            shutil.copy(fasta_file_path[i], f"{self.shared_folder}/input")
        print("="*80)
        print("Copied protein fastas to the input folder for the subprocess step")

        # now that the protein files are in /input, change all the extensions there
        mypath = Path(f"{self.shared_folder}/input")
        for file in os.listdir(mypath):
            print("Filename in loop: ", file)
            src = os.path.join(mypath, file)
            dst = os.path.join(mypath, file + '.faa')
            if not os.path.exists(dst):  # check if the file doesn't exist
                os.rename(src, dst)

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
                    if file.endswith('.csv') or file.endswith('.png'):
                        zip_file.write(os.path.join(root, file),
                                       os.path.join(os.path.basename(root), file))

        output_files.append({
            'path': result_file,
            'name': os.path.basename(result_file),
            'label': os.path.basename(result_file),
            'description': 'Files generated by Snekmer Search'})

        # Build a Report and return
        print('Section: build report data.')
        report_params = {
            'message': 'Kmer input was ' + str(k),
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
