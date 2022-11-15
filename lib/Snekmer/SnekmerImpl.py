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
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI

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
    GIT_COMMIT_HASH = "ca751677be245833f13674fdbc5e41a9c53bdb6e"

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
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
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
           parameter "output_genome_name" of String
        :returns: instance of type "SnekmerSearchOutput" (Output parameters
           for Snekmer Search. report_name - the name of the
           KBaseReport.Report workspace object. report_ref - the workspace
           reference of the report.) -> structure: parameter "report_name" of
           String, parameter "report_ref" of String, parameter
           "output_genome_ref" of type "genome_ref" (Reference to a Genome
           object in the workspace @id ws KBaseGenomes.Genome)
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
        if 'output_genome_name' not in params:
            raise ValueError('Parameter output_genome_name is not set in input arguments')
        output_genome_name = params['output_genome_name']

        # testing diff between dfu.get_objects and wsClient.get_objects2 (which is used in GenomeSetToFasta)
        obj_dfu_get_obj = self.dfu.get_objects({'object_refs': [object_ref]})
        print("Using DataFileUtil.get_objects: ", obj_dfu_get_obj)
        print("")

        # accessing different parts of dfu.get_objects output
        dfu_elements = obj_dfu_get_obj['data'][0]['data']['elements']
        print("dfu_elements: ", dfu_elements)
        print("")
        dfu_keys = list(dfu_elements)
        # prints the scientific name (at least when testing locally)
        print("dfu_keys should be the genome_ids: ", dfu_keys)

        refs = []
        for i in dfu_keys:
            refs.append(dfu_elements[i]['ref'])
        # prints the reference ids (only) for the genomes within the genomeset
        print("list of genome refs in the genome_set: ", refs)
        print(refs[0])

        ref_list = list(dfu_elements.values())
        print("")
        print("ref_list: ", ref_list)
        print("")

        # testing use of genome_api
        genome_data = []
        for i in refs:
            print('made it into the loop for ', i)
            genome_data.append(self.genome_api.get_genome_v1({"genomes": [{"ref": i}],
                                                'downgrade': 0})["genomes"][0])

        print("")
        print("genome_data length: ", len(genome_data))

        #for i in genome_data:
         #   print("look at some functions: ", i["data"]["features"][0]["functions"])
        print("first genome, first five functions original ")
        print("")

        print("data type of genome_data: ", type(genome_data), "\n")
        print("genome_data[0].keys(): ", genome_data[0].keys(), "\n")
        print("genome_data[0]['info']: ", genome_data[0]['info'], "\n")
        print("genome_data[0]['data'].keys(): ", genome_data[0]['data'].keys(), "\n")
        print("genome_data[0]['data']['features'][0]: ", genome_data[0]['data']['features'][0], "\n")

        # functions has a list in it already, append
        if 'functions' in genome_data[0]['data']['features'][0]:
            print("list- genome_data has functions: ", genome_data[0]['data']['features'][0]['functions'], "\n")
            genome_data[0]['data']['features'][0]['functions'].append('Abby added this')
            print("after Abby appended to 'functions': ", genome_data[0]['data']['features'][0]['functions'], "\n")

        # function is a str, not list yet
        # turn into list, then append to it
        if 'function' in genome_data[0]['data']['features'][0]:
            print("str- genome_data has function: ", genome_data[0]['data']['features'][0]['function'], "\n")
            genome_data[0]['data']['features'][0]['function'] = [genome_data[0]['data']['features'][0]['function']]
            genome_data[0]['data']['features'][0]['function'].append('Abby added this')
            print("after Abby appended to 'function': ", genome_data[0]['data']['features'][0]['function'], "\n")



        sys.exit()
        # look at feature functions for one genome, before and after adding to the function list
        for i in range(5):
            print(genome_data[0]['data']["features"][i]["functions"])
            genome_data[0]['data']["features"][i]["functions"].append("Added by Abby")
            print(genome_data[0]['data']["features"][i]["functions"])
            print("")

        print("")
        # did the functions stay appended? yes!
        print("after functions were appended: ")
        for i in range(5):
            print(genome_data[0]['data']["features"][i]["functions"])
            print("")

        # see if i can save the above edited functions into the actual genome object
        # code from prokka
        # the updated genome
        # "provenance": self.ctx.provenance()
        # in prokka, 'name' param is given by ui
        print("genome name: ", dfu_keys[0], "\n")
        output_workspace_name = 'some name i guess'
        print("workspace_name: ", workspace_name, "\n")
        #print("genome_data[0]['data']", genome_data[0]['data'], "\n")

        print("Start gfu.save_one_genome \n")
        saved = self.gfu.save_one_genome({"workspace": workspace_name,
                                          "name": output_genome_name,
                                          "data": genome_data[0]['data']})
        print("saved: ", saved)
        print("")
        info = saved["info"]
        genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
        print("genome ref: ", genome_ref)
        print("")
        #annotated_genome = namedtuple("annotated_genome",
         #                             "genome_ref function_summary_filepath ontology_summary_filepath stats")

        report_message = "The annotation of one genome, for the first 5 features, was successful"
        # temporary report section, to return the (incorrectly) annotated genome
        print('Section: build report data.')
        print("")
        report_params = {
            'message': report_message,
            'workspace_name': workspace_name,
            'objects_created': [{"ref": genome_ref, "description": "Annotated genome by Abby!"}]
        }
        print("report_params: ", report_params, "\n")

        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)
        print("report info: ", report_info, "\n")
        # construct the output to send back
        output = {'output_genome_ref': genome_ref, 'report_name': report_info['name'], 'report_ref': report_info['ref']}
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
