# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
from pprint import pprint
import shutil

from Snekmer.SnekmerImpl import Snekmer
from Snekmer.SnekmerServer import MethodContext
from Snekmer.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.GenomeFileUtilClient import GenomeFileUtil


class SnekmerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('Snekmer'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'Snekmer',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = Snekmer(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        #suffix = int(time.time() * 1000)
        #cls.wsName = "test_Snekmer_" + str(suffix)
        #ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_Snekmer_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # get obj_ref in form D/D/D
    def get_obj_ref_from_obj_info(self, obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        return '/'.join([str(obj_info[WSID_I]), str(obj_info[OBJID_I]), str(obj_info[VERSION_I])])

    # retrieve stored obj info
    def _get_stored_obj_info(self, obj_type, obj_name, item_i=0):
        infoAttr = obj_type + 'Info_list'  # e.g. 'ama' or 'genome'
        nameAttr = obj_type + 'Name_list'
        if hasattr(self.__class__, infoAttr):
            try:
                info_list = getattr(self.__class__, infoAttr)
                name_list = getattr(self.__class__, nameAttr)
                info = info_list[item_i]
                name = name_list[item_i]
                if info != None:
                    if name != obj_name:
                        info_list[item_i] = None
                        name_list[item_i] = None
                        setattr(self.__class__, infoAttr, info_list)
                        setattr(self.__class__, nameAttr, name_list)
                    else:
                        return info
            except:
                pass
        return None

    # save stored obj info
    def _save_stored_obj_info(self, obj_type, obj_info, obj_name, item_i=0):
        infoAttr = obj_type + 'Info_list'  # e.g. 'ama' or 'genome'
        nameAttr = obj_type + 'Name_list'
        if not hasattr(self.__class__, infoAttr):
            setattr(self.__class__, infoAttr, [])
            setattr(self.__class__, nameAttr, [])

        info_list = getattr(self.__class__, infoAttr)
        name_list = getattr(self.__class__, nameAttr)
        for i in range(item_i + 1):
            try:
                assigned = info_list[i]
            except:
                info_list.append(None)
                name_list.append(None)
        info_list[item_i] = obj_info
        name_list[item_i] = obj_name
        setattr(self.__class__, infoAttr, info_list)
        setattr(self.__class__, nameAttr, name_list)
        return

    # call this method to get the WS object info of a Genome
    #   (will upload the example data if this is the first time the method is called during tests)
    def getGenomeInfo(self, genome_basename, item_i=0):
        info = self._get_stored_obj_info('genome', genome_basename, item_i)
        if info != None:
            return info

        # 1) transform genbank to kbase genome object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        genome_data_file = 'data/genomes/' + genome_basename + '.gbff.gz'
        genome_file = os.path.join(shared_dir, os.path.basename(genome_data_file))
        shutil.copy(genome_data_file, genome_file)

        SERVICE_VER = 'release'
        GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                             token=self.getContext()['token'],
                             service_ver=SERVICE_VER
                             )
        print("UPLOADING genome: " + genome_basename + " to WORKSPACE " + self.getWsName() + " ...")
        genome_upload_result = GFU.genbank_to_genome({'file': {'path': genome_file},
                                                      'workspace_name': self.getWsName(),
                                                      'genome_name': genome_basename
                                                      })
        print('Genome upload result')
        pprint(genome_upload_result)
        genome_ref = genome_upload_result['genome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': genome_ref}]})[0]

        # 2) store it
        self._save_stored_obj_info('genome', new_obj_info, genome_basename, item_i)

        return new_obj_info

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # 66073/19/1 is Shewanella_oneidensis_MR-1 KBaseGenomes.Genome object
    # 66073/2/1 is Escherichia_coli_K-12_MG1655 KBaseGenomes.Genome object
    # 66073/25/3 is MyTestGenomeSet, a KBaseSearch.Genome object, containing the above Genomes
    @unittest.skip("skipped test_run_Snekmer_search")
    def test_run_Snekmer_search(self):
        ref = ["66073/25/3"]
        ret = self.serviceImpl.run_Snekmer_search(
            self.ctx,
            {
             'workspace_name': self.wsName,
             'object_ref': ref,
             'k': 4,
             'alphabet': "miqs"
             }
        )
        # print(result)

    # Test GenomeSet objecct input into Snekmer Search
    def test_Snekmer_search_GenomeSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I,
         META_I] = list(range(11))  # object_info tuple
        genomeSet_name = 'testGenomeSet'

        load_genomes = [
            {'file': 'GCF_001566335.1_ASM156633v1_genomic',
             'sciname': 'E. coli K-12 MG1655'
             },
            {'file': 'GCF_000021385.1_ASM2138v1_genomic',
             'sciname': 'D. vulgaris str. Miyazaki F'
             }
        ]
        for genome_i, genome in enumerate(load_genomes):
            load_genomes[genome_i]['ref'] = self.get_obj_ref_from_obj_info(self.getGenomeInfo(genome['file'], genome_i))

        # create GenomeSet
        testGS = {
            'description': 'two genomes',
            'elements': dict()
        }
        for genome_i, genome in enumerate(load_genomes):
            testGS['elements'][genome['sciname']] = {'ref': genome['ref']}

        print("After adding to the testGS")
        pprint(testGS)

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),
                                                    'objects': [
                                                        {
                                                            'type': 'KBaseSearch.GenomeSet',
                                                            'data': testGS,
                                                            'name': genomeSet_name,
                                                            'meta': {},
                                                            'provenance': [
                                                                {
                                                                    'service': 'Snekmer',
                                                                    'method': 'run_Snekmer_search'
                                                                }
                                                            ]
                                                        }]
                                                    })[0]
        print("Obj_info: ")
        pprint(obj_info)
        print("=" * 80)
        target_genomeSet_ref = self.get_obj_ref_from_obj_info(obj_info)
        print("target_genomeSet_ref: ")
        pprint(target_genomeSet_ref)

        parameters = {'workspace_name': self.getWsName(),
                      'object_ref': target_genomeSet_ref,
                      'k': 6,
                      'alphabet': "standard",
                      'output_genome_name': 'TheOutPutGenomeName'
                      }

        ret = self.getImpl().run_Snekmer_search(self.getContext(), parameters)[0]
