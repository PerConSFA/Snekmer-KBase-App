# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import yaml
from pprint import pformat
from Bio import SeqIO

from installed_clients.KBaseReportClient import KBaseReport
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
    GIT_COMMIT_HASH = "6b4e77303f0ef97ee64bbc929c15f076ed31a7e7"

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
           parameters for Snekmer Search. workspace_name - the name of the
           workspace for input/output object_ref - Genome object with Protein
           Translation sequence in the Feature k - kmer length for features
           alphabet - mapping function for reduced amino acid sequences
           min_rep_thresh - min number of sequences to include feature for
           prefiltering processes - for parallelization) -> structure:
           parameter "workspace_name" of String, parameter "object_ref" of
           String, parameter "k" of Long, parameter "alphabet" of Long,
           parameter "min_rep_thresh" of Long, parameter "processes" of Long
        :returns: instance of type "SnekmerSearchOutput" (Output parameters
           for Snekmer Search. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_Snekmer_search
        print('Input parameters: ' + pformat(params))
        print(" ")

        # check inputs
        object_ref = params['object_ref']
        if 'k' not in params:
            raise ValueError('Parameter kmer is not set in input arguments')
        k = params['k']
        if 'processes' not in params:
            raise ValueError('Parameter processes is not set in input arguments')
        processes = params['processes']
        workspace_name = params['workspace_name']
        if 'alphabet' not in params:
            raise ValueError('Parameter alphabet is not set in input arguments')
        alphabet = params['alphabet']
        if 'min_rep_thresh' not in params:
            raise ValueError('Parameter min_rep_thresh is not set in input arguments')
        min_rep_thresh = params['min_rep_thresh']

        # add these params to the config.yaml
        new_params = {'k': params['k'], 'alphabet': params['alphabet'], 'min_rep_thresh': params['min_rep_thresh'],
                      'processes': params['processes']}
        #print('New params: ' + pformat(new_params))
        print('add params to config: ')
        with open('/kb/module/data/config.yaml', 'r') as file:
            my_config = yaml.safe_load(file)
            my_config.update(new_params)
        with open('/kb/module/data/config.yaml', 'w') as file:
            yaml.safe_dump(my_config, file)

        # testing that the appended config file is what I intended
        print(' ')
        print('Reading the appended config: ')
        f = open('/kb/module/data/config.yaml', 'r')
        content = f.read()
        print(content)
        f.close()
        print(' ')

        # Step - Build a Report and return
        print('Section: build report data.')
        report_data = {
            'objects_created': [],
            'text_message': 'Kmer input was ' + str(k) + ' using the ' + str(alphabet) + ' alphabet with a rep threshold of ' +
                            str(min_rep_thresh) + ' with ' + str(processes) + ' processes.'
        }
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': report_data, 'workspace_name': workspace_name})

        # STEP 6: contruct the output to send back
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref'],
                  'k': k,
                  'alphabet': alphabet,
                  'min_rep_thresh': min_rep_thresh,
                  'processes': processes
                  }
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
