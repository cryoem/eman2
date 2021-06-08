from __future__ import print_function
from __future__ import division

# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

import unittest
from mpi import *
from sphire.libpy import sp_global_def
import numpy
import weakref
import EMAN2db
import EMAN2_cppwrap
from os import path, mkdir


from EMAN2_cppwrap import Transform
from EMAN2 import EMAN2Ctf
import EMAN2

import matplotlib
from matplotlib import pyplot
from pickle import *
import sys
import os
absolute_import = (sys.version_info[0] >= 3)
from bsddb3 import db
import copy


"""
Input to this class is a bdb file with images
"""

bdb_path = "../resources_tests/03_PARTICLES_BDB/mpi_proc_007/EMAN2DB"
bdbname = "TcdA1-0187_frames_ptcls"
bdbfile = 'bdb:'+os.path.join(os.path.dirname(__file__), "../resources_tests/03_PARTICLES_BDB/mpi_proc_007/TcdA1-0187_frames_ptcls")


class Test_functions_outside_class(unittest.TestCase):

    def test_db_convert_path(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        self.assertEqual(newpath, str('bdb:' + bdb_path.replace('/EMAN2DB', '') + '#' + bdbname))

    def test_db_parse_path(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        path,dictname,keys = EMAN2db.db_parse_path(newpath)
        self.assertEqual(dictname , bdbname)
        self.assertEqual(keys, None)

    def test_db_open_dict(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        opendict = EMAN2db.db_open_dict(newpath, ro=True)
        dict_keys = ['name', 'parent', 'dbenv', 'lock', 'file', 'rohint', 'lasttime', 'opencount', 'path', 'txn', 'bdb', 'isro', 'key_translation_dict']
        self.assertEqual(list(opendict.__dict__.keys()) , dict_keys )
        self.assertEqual(opendict.__dict__['name'] , bdbname)

    def test_db_close_dict(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        opendict = EMAN2db.db_close_dict(newpath)
        with self.assertRaises(AttributeError) as cm_new:
            opendict.__dict__['name']
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute '__dict__'")

    def test_db_check_dict(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        status = EMAN2db.db_check_dict(newpath)
        self.assertEqual(status,True)

    def test_db_list_dict(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)
        listdic = EMAN2db.db_list_dicts(newpath)
        self.assertEqual(len(listdic), 3)

    def test_db_read_image_and_images(self):
        """
        case 1
        """
        img = EMAN2.EMData()
        EMAN2.db_read_image(img,bdbfile,0 )
        img1 = EMAN2.db_read_images(bdbfile)
        self.assertTrue(numpy.array_equal(img1[0].get_2dview(), img.get_2dview()))

        """
        case 2
        """
        img = EMAN2.EMData()
        EMAN2.db_read_image(img, bdbfile,0 )
        img2 = EMAN2.EMData().read_images(bdbfile)
        self.assertTrue(numpy.array_equal(img2[0].get_2dview(), img.get_2dview()))

        """
        case 3 
        """
        img = EMAN2.EMData()
        img.read_image(bdbfile, 5)
        img1 = EMAN2.EMData.read_images(bdbfile)
        self.assertTrue(numpy.array_equal(img.get_2dview(), img1[5].get_2dview()))


    def test_db_get_image_count(self):
        case1_bdb_count = EMAN2db.db_get_image_count(bdbfile)
        self.assertEqual(case1_bdb_count, 127)

    def test_get_image_info(self):
        fsp = bdbfile
        iminfo  = EMAN2db.db_get_image_info(fsp)
        self.assertEqual(iminfo[0], 127)
        self.assertEqual(iminfo[1], (352,352,1))

    def test_db_get_all_attributes(self):
        fsp = bdbfile
        atr = EMAN2db.db_get_all_attributes(fsp, 'nx')
        attr_list = [352,352,352,352,352]
        self.assertEqual(atr[0:5], attr_list)


class Test_BDB(unittest.TestCase):

    @unittest.skip('creating dummy files')
    def test_open_db(self):
        dic = EMAN2db.EMAN2DB()
        print('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname)
        dicc = dic.open_db(path= str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        dic.close_dict(dicc)
        self.assertEqual(dicc.path, str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))

    @unittest.skip('creating dummy files')
    def test_bdb_init(self):
        dic = EMAN2db.EMAN2DB()
        dic.__init__(path= str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        self.assertEqual(dic.path, str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))

    @unittest.skip('unrelated')
    def test_bdb_getitem_and_open_dict(self):
        ddb=EMAN2db.EMAN2DB.open_db(bdbfile)
        ddb.open_dict(bdbname, ro=True)
        print(ddb.__dict__)
        data = ddb.__getitem__(key=bdbname)
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        self.assertTrue(numpy.array_equal(a.keys(), data.keys()))

class Test_Dbdict(unittest.TestCase):

    def test_get_attrib(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        img =  a[0]
        print(img.get_attr('MRC.gamma'))
        fsp = str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname)
        atr = EMAN2db.db_get_all_attributes(fsp, 'MRC.gamma')
        self.assertEqual(img.get_attr('MRC.gamma') , atr[0])

    def test_init(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        self.assertEqual(bdbname, a.name)
        self.assertEqual(bdb_path, a.path)
        self.assertEqual({}, a.key_translation_dict)

    def test_realopen(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        a.realopen()
        self.assertEqual(a.keys()[-1], loads(list(a.key_translation_dict.keys())[-1]))

    @unittest.skip("Close function doesnot work, i can still access the keys")
    def test_close(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        a.close
        print(a.keys())
        print(a[0])

    def test_len(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        self.assertEqual(a.__len__(), 127)

    def test_setitem(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        beforedel = a.keys()
        key,value = "disco" , "myname"
        a.__setitem__(key, value)
        a.__delitem__("disco")
        afterdel = a.keys()
        self.assertEqual(beforedel, afterdel)

    def test_getitem(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        self.assertTrue(numpy.array_equal(a[3].get_2dview(), a.__getitem__(3).get_2dview()))

    def test_delitem(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        beforedel = a.keys()
        value , key = "disco" , "myname"
        a.set(value, key)
        a.__delitem__("myname")
        afterdel = a.keys()
        self.assertEqual(beforedel, afterdel )

    def test_contains(self):
        a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
        print(a.__contains__('maxrec'))
        self.assertEqual(a.__contains__(5), True)

    def test_item_type(self):
        a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
        self.assertEqual(a.item_type(0) ,  type(a[0]))

    def test_keys_work(self):
        a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
        keys = a.keys()
        real_keys =[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                    27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
                    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
                    77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
                    101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
                    121, 122, 123, 124, 125, 126, 'maxrec']
        self.assertEqual(keys, real_keys)

    def test_items_work(self):
        a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
        items = a.items()
        keys = a.keys()
        keys_in_items = [key for key,value in items]
        values = a.values()
        self.assertEqual(keys, keys_in_items)

    def test_has_key(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        keys = a.keys()
        self.assertEqual(a.has_key(keys[-1]), True)

    @unittest.skip
    def test_data_path(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        values = a.values()
        keys = a.keys()
        self.assertEqual(a.get_data_path(keys[2]), os.path.join(bdb_path.split('EMAN2DB')
                                    [0], bdbname+'.mrcs'))

    def test_get(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        keys = a.keys()
        items = a.items()
        value_in_items = [value for key, value in items]
        record = a.get(keys[-9])
        self.assertTrue(numpy.array_equal(value_in_items[-9].get_2dview() , record.get_2dview() ))

    # @unittest.skip
    # def test_get_header_items_keys_values(self):
    #     a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
    #     header = a.get_header(0)
    #     header_orig_keys = ['MRC.nzstart', 'apix_x', 'MRC.gamma', 'xform.projection', 'changecount', 'maximum', 'data_path', 'MRC.beta', 'data_n', 'ctf_applied', 'minimum', 'ptcl_source_coord_id', 'mean_nonzero', 'MRC.label0', 'MRC.mapr', 'MRC.maps', 'ny', 'MRC.ispg', 'apix_y', 'MRC.xlen', 'HostEndian', 'nx', 'MRC.nsymbt', 'nz', 'MRC.mapc', 'ptcl_source_coord', 'sigma_nonzero', 'MRC.nxstart', 'ptcl_source_image', 'is_complex_ri', 'relion_max_prob_dist', 'MRC.nz', 'MRC.nx', 'MRC.ny', 'MRC.maximum', 'apix_z', 'MRC.nlabels', 'source_n', 'MRC.mean', 'origin_z', 'is_complex', 'is_complex_x', 'relion_norm_correct', 'chunk_id', 'ImageEndian', 'origin_y', 'origin_x', 'MRC.nystart', 'MRC.zlen', 'MRC.rms', 'resample_ratio', 'MRC.machinestamp', 'data_source', 'ptcl_source_apix', 'adnan_n', 'datatype', 'sigma', 'ptcl_source_relion', 'source_path', 'MRC.minimum', 'ctf', 'square_sum', 'MRC.mz', 'MRC.my', 'MRC.mx', 'MRC.alpha', 'MRC.ylen', 'mean']        # print(list(sorted(header.keys())) == list(sorted(header_orig_keys)))
    #     sb = '[(\'HostEndian\', \'little\'), (\'ImageEndian\', \'little\'), (\'MRC.alpha\', 90.0), (\'MRC.beta\', 90.0), (\'MRC.gamma\', 90.0), (\'MRC.ispg\', 0), (\'MRC.label0\', \'Relion    21-Oct-19  16:12:54\'), (\'MRC.machinestamp\', 16708), (\'MRC.mapc\', 1), (\'MRC.mapr\', 2), (\'MRC.maps\', 3), (\'MRC.maximum\', 4.168662071228027), (\'MRC.mean\', 0.0022823389153927565), (\'MRC.minimum\', -6.85170316696167), (\'MRC.mx\', 360), (\'MRC.my\', 360), (\'MRC.mz\', 1), (\'MRC.nlabels\', 1), (\'MRC.nsymbt\', 0), (\'MRC.nx\', 360), (\'MRC.nxstart\', 0), (\'MRC.ny\', 360), (\'MRC.nystart\', 0), (\'MRC.nz\', 1), (\'MRC.nzstart\', 0), (\'MRC.rms\', 1.0004879236221313), (\'MRC.xlen\', 318.6000061035156), (\'MRC.ylen\', 318.6000061035156), (\'MRC.zlen\', 1.0), (\'adnan_n\', 0), (\'apix_x\', 1.0), (\'apix_y\', 1.0), (\'apix_z\', 1.0), (\'changecount\', 0), (\'chunk_id\', 0), (\'ctf\', EMAN2Ctf().from_dict({"defocus":1.07294,"bfactor":0.0,"ampcont":10.00,"apix":0.885,"voltage":200.00,"cs":1.40}) ...), (\'ctf_applied\', 0), (\'data_n\', 0), (\'data_path\', \'/home/adnan/PycharmProjects/DoseWeighting/Newfolder/20170629_00021_frameImage_ptcls.mrcs\'), (\'data_source\', \'bdb:Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#20170629_00021_frameImage_ptcls\'), (\'datatype\', 7), (\'is_complex\', 0), (\'is_complex_ri\', 1), (\'is_complex_x\', 0), (\'maximum\', 3.379175901412964), (\'mean\', 0.0033558777067810297), (\'mean_nonzero\', 0.0033558777067810297), (\'minimum\', -6.456167697906494), (\'nx\', 360), (\'ny\', 360), (\'nz\', 1), (\'origin_x\', 0.0), (\'origin_y\', 0.0), (\'origin_z\', 0.0), (\'ptcl_source_apix\', 0.885), (\'ptcl_source_coord\', [2937, 3437]), (\'ptcl_source_coord_id\', 0), (\'ptcl_source_image\', \'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc\'), (\'ptcl_source_relion\', \'000001@Extract/job051/MOVIES_RELION/20170629_00021_frameImage.mrcs\'), (\'relion_max_prob_dist\', 0.390257), (\'relion_norm_correct\', 0.504526), (\'resample_ratio\', 1.0), (\'sigma\', 1.0007057189941406), (\'sigma_nonzero\', 1.0007057189941406), (\'source_n\', 0), (\'source_path\', \'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#sphire_relion_stack\'), (\'square_sum\', 129783.453125), (\'xform.projection\', Transform({\'az\':259.319,\'alt\':77.866,\'phi\':23.693,\'tx\':-1.89,\'ty\':-0.30,\'tz\':0.00,\'mirror\':0,\'scale\':1.0000,\'type\':\'eman\'}))]'
    #     val = ['-6.45616769791', '-6.85170316696', '/home/adnan/PycharmProjects/DoseWeighting/Newfolder/20170629_00021_frameImage_ptcls.mrcs', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0.0', '0.0', '0.0', '0.00228233891539', '0.00335587770678', '0.00335587770678', '0.390257', '0.504526', '0.885', '000001@Extract/job051/MOVIES_RELION/20170629_00021_frameImage.mrcs', '1', '1', '1', '1', '1', '1', '1.0', '1.0', '1.0', '1.0', '1.0', '1.00048792362', '1.00070571899', '1.00070571899', '129783.453125', '16708', '2', '3', '3.37917590141', '318.600006104', '318.600006104', '360', '360', '360', '360', '360', '360', '4.16866207123', '7', '90.0', '90.0', '90.0', 'EMAN2Ctf().from_dict({"defocus":1.07294,"bfactor":0.0,"ampcont":10.00,"apix":0.885,"voltage":200.00,"cs":1.40}) ...', 'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc', 'Relion    21-Oct-19  16:12:54', "Transform({'az':259.319,'alt':77.866,'phi':23.693,'tx':-1.89,'ty':-0.30,'tz':0.00,'mirror':0,'scale':1.0000,'type':'eman'})", '[2937, 3437]', 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#sphire_relion_stack', 'bdb:Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#20170629_00021_frameImage_ptcls', 'little', 'little']
    #     self.assertEqual(list(sorted(header.keys())), list(sorted(header_orig_keys)))
    #     self.assertTrue(sb, repr(sorted(header.items())))
    #     self.assertTrue(repr(val), repr(list(sorted(map(str, header.values())))))

    def test_set(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        value , key = "disco" , "myname"
        key_before_set = a.keys()
        a.set(value, key)
        a.__delitem__("myname")
        key_after_set = a.keys()
        self.assertTrue(numpy.array_equal(key_before_set, key_after_set))

    @unittest.skip("for the time being")
    def test_set_header(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        a.close()
        key =  0
        data = a[0]
        data1 = a[a.keys()[0]]
        data2 = a[a.keys()[-3]]
        a.set_header(8, data1)
        self.assertTrue(numpy.array_equal(data2.get_2dview(), data.get_2dview()))


class Test_Stacks(unittest.TestCase):
    def test_stack_substack(self):
        e = EMAN2.EMData()
        e.read_image(bdbfile, 5)
        d = EMAN2.EMData()
        d.read_image(bdbfile, 5)
        img = EMAN2.EMData.read_images(bdbfile)
        self.assertTrue(numpy.array_equal(e.get_2dview(), d.get_2dview()))
        self.assertTrue(numpy.array_equal(img[5].get_2dview(), d.get_2dview()))

    @unittest.skip
    def test_stack_dict_reading(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        img = EMAN2.EMData()
        EMAN2.db_read_image(img,bdbfile,0)
        EMAN2.db_close_dict(a)

    @unittest.skip('skipping writing files ')
    def test_stack_writing(self):
        img = EMAN2.EMData()
        img.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        EMAN2db.db_close_dict(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        newimg = EMAN2.EMData()
        newimg = copy.deepcopy(img)
        if os.path.isfile(os.path.abspath(bdb_path)+'test_img.bdb'):
            os.remove(os.path.abspath(bdb_path)+'test_img.bdb')
        else :
            pass
        newimg.write_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + 'test_img'))
        if os.path.isfile(os.path.abspath(bdb_path)+'test_img.bdb'):
            os.remove(os.path.abspath(bdb_path)+'test_img.bdb')
        else :
            pass

    @unittest.skip('display not necessry to check')
    def test_read_new_file(self):
        imgnew = EMAN2.EMData()
        imgnew.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB', '') + '#' + 'test_img'))
        EMAN2.display(imgnew)
        if os.path.isfile(os.path.abspath(bdb_path)+'test_img.bdb'):
            os.remove(os.path.abspath(bdb_path)+'test_img.bdb')
        else :
            pass




