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

        logging.info("Grabbing the Genome data from the input GenomeSet.")
        # accessing different parts of dfu.get_objects output
        obj_dfu_get_obj = self.dfu.get_objects({'object_refs': [object_ref]})

        # format is 'dfu_elements: {'param0': {'ref': '66073/19/1'}, 'param1': {'ref': '66073/2/1'}}' in appdev
        # the keys are the scientific names (genome_ids) when testing locally
        dfu_elements = obj_dfu_get_obj['data'][0]['data']['elements']
        dfu_keys = list(dfu_elements)

        # get list of refs for the genomes within the genomeset
        refs = []
        for i in dfu_keys:
            refs.append(dfu_elements[i]['ref'])

        # grab the current genome_data
        genome_data = []
        for i in refs:
            print('made it into the loop for ', i)
            genome_data.append(self.genome_api.get_genome_v1({"genomes": [{"ref": i}],
                                                              'downgrade': 0})["genomes"][0])

        logging.info("Annotating the Genomes.")
        # for now annotate the first 5 genes with my lovely message to prove I can do it
        # need if statements for 'functions' vs 'function' because of the differences in genome object versions
        for j in genome_data:
            print('in genome_data loop')
            for i in range(5):
                if 'functions' in j['data']['features'][i]:
                    print("i: ", i)
                    print("has functions: ", j['data']["features"][i]["functions"])
                    j['data']["features"][i]["functions"].append("Added by Abby")
                    print("add new: ", j['data']["features"][i]["functions"])
                    print("")

                if 'function' in j['data']['features'][i]:
                    print("i: ", i)
                    print("has function: ", j['data']['features'][i]['function'])
                    annotationString = 'Abby added this'
                    j['data']['features'][i]['function'] = " , ".join(
                        [j['data']['features'][i]['function'], annotationString])

                    print("add new: ", j['data']['features'][i]['function'])
                    print("")

        logging.info("Saving the annotated Genomes as individual Genome objects.")
        # save the annotated genomes as new genome objects, with new refs
        new_refs = []
        new_names = []
        for i in genome_data:
            # format the organism name into what's acceptable for an object name
            # example- gfu.save_one_genome claimed "Desulfovibrio vulgaris str. 'Miyazaki F'" had an illegal character
            stringName = i['data']['scientific_name']
            obj_name = "_".join(stringName.split())
            obj_name2 = obj_name.replace("'", "_")
            # save the annotated genome object and grab its info
            info = self.gfu.save_one_genome({"workspace": workspace_name,
                                             "name": obj_name2,
                                             "data": i['data']})["info"]
            # get ref of new annotated genome
            save_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
            new_refs.append(save_ref)
            new_names.append(stringName)

        logging.info("Saving the new Genomes into a new GenomeSet object.")
        # save annotated genomes into genomeset object, then get that new ref to pass into the report
        new_gs_data = {'description': 'blah for now', 'elements': dict()}

        for i, j in zip(new_names, new_refs):
            new_gs_data['elements'][i] = {'ref': j}

        wsid = self.dfu.ws_name_to_id(workspace_name)
        obj_info = self.dfu.save_objects({'id': wsid,
                                          'objects': [
                                                    {'type': 'KBaseSearch.GenomeSet',
                                                     'data': new_gs_data,
                                                     'name': output_genome_name,
                                                     'meta': {},
                                                     'provenance': [
                                                            {'service': 'Snekmer',
                                                             'method': 'run_Snekmer_search'
                                                             }]
                                                     }]
                                          })[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I,
         META_I] = list(range(11))
        genomeSet_ref = '{}/{}/{}'.format(obj_info[WSID_I], obj_info[OBJID_I], obj_info[VERSION_I])

        logging.info("Building the output KBaseReport.")
        # temporary report section
        # returns the (incorrectly) annotated genomes in genomeset object
        report_message = "The annotation of 2 genomes, for each of their first 5 features, was successful"
        report_params = {
            'message': report_message,
            'workspace_name': workspace_name,
            'objects_created': [{"ref": genomeSet_ref, "description": "Annotated genome by Abby!"}]
        }

        report_client = KBaseReport(self.callback_url)
        report_info = report_client.create_extended_report(report_params)

        # construct the output to send back
        output = {'output_genome_ref': genomeSet_ref,
                  'report_name': report_info['name'],
                  'report_ref': report_info['ref']}

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
