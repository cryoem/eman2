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
from numpy import allclose, array_equal
from mpi import *
import sp_global_def
import numpy
import weakref
import EMAN2db
import EMAN2_cppwrap
from os import path, mkdir
from pickle import loads as pickle_loads
import threading
import time
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

from EMAN2_cppwrap import Transform
from EMAN2 import EMAN2Ctf
import EMAN2

import matplotlib
from matplotlib import pyplot
# sp_global_def.BATCH = True
# sp_global_def.MPI = False
# mpi_init(0, [])

from pickle import *
import sys
import os
absolute_import = (sys.version_info[0] >= 3)
from bsddb3 import db
import copy


"""
Input to this class is a bdb file with images
"""

bdb_path = '/home/adnan/PycharmProjects/DoseWeighting/Newfolder/EMAN2DB'
bdbname ='20170629_00021_frameImage_ptcls'


class Test_functions_outside_class(unittest.TestCase):

    def test_db_convert_path(self):
        path_file = os.path.join(bdb_path, bdbname)
        # print(path_file)
        newpath = EMAN2db.db_convert_path(path_file)
        # print(newpath)
        self.assertEqual(newpath, str('bdb:' + bdb_path.replace('/EMAN2DB', '') + '#' + bdbname))

    def test_db_parse_path(self):
        path_file = os.path.join(bdb_path, bdbname)
        newpath = EMAN2db.db_convert_path(path_file)

        path,dictname,keys = EMAN2db.db_parse_path(newpath)
        print(path)
        print(dictname)
        print(keys)
        self.assertEqual(path, bdb_path.replace('/EMAN2DB', ''))
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
        self.assertEqual(len(listdic), 24)

    def test_db_read_image_and_images(self):
        import EMAN2
        """
        case 1
        """
        img = EMAN2db.EMData()
        EMAN2db.db_read_image(img, "bdb:/home/adnan/DemoResults/03_Particles#stack",0 )
        img1 = EMAN2db.db_read_images("bdb:/home/adnan/DemoResults/03_Particles#stack")
        self.assertTrue(numpy.array_equal(img1[0].get_2dview(), img.get_2dview()))

        """
        case 2
        """
        img = EMAN2db.EMData()
        EMAN2db.db_read_image(img, "bdb:/home/adnan/DemoResults/03_Particles/mpi_proc_000#TcdA1-0010_frames_ptcls",0 )
        img1 = EMAN2db.db_read_images("bdb:/home/adnan/DemoResults/03_Particles/mpi_proc_000#TcdA1-0010_frames_ptcls")
        self.assertTrue(numpy.array_equal(img1[0].get_2dview(), img.get_2dview()))

        """
        case 3 
        """
        img = EMAN2db.EMData()
        img.read_image("bdb:/home/adnan/DemoResults/04_ISAC#stack_ali2d", 5)
        img1 = EMAN2db.EMData.read_images("bdb:/home/adnan/DemoResults/04_ISAC#case4stack")
        self.assertTrue(numpy.array_equal(e.get_2dview(), d.get_2dview()))

        """
        case 4
        """
        img = EMAN2db.EMData()
        img.read_image("bdb:/home/adnan/DemoResults/04_ISAC#case4stack", 5)
        self.assertTrue(numpy.array_equal(e.get_2dview(), d.get_2dview()))




    def test_db_get_image_count(self):
        case1_bdb_count = EMAN2db.db_get_image_count("bdb:/home/adnan/DemoResults/03_Particles#stack")
        # print(case1_bdb_count)
        cas2_bdb_count = EMAN2db.db_get_image_count(
            "bdb:/home/adnan/DemoResults/03_PARTICLES/mpi_proc_000#TcdA1-0010_frames_ptcls")
        # print(cas2_bdb_count)
        self.assertEqual(case1_bdb_count, 6989)
        self.assertEqual(cas2_bdb_count, 64)

    @unittest.skip("For the time being, i have left this part")
    def test_db_write_image(self):
        pass

    def test_get_image_info(self):
        fsp = "bdb:/home/adnan/DemoResults/03_Particles/mpi_proc_000#TcdA1-0010_frames_ptcls"
        iminfo  = EMAN2db.db_get_image_info(fsp)
        self.assertEqual(iminfo[0], 64)
        self.assertEqual(iminfo[1], (352,352,1))

    def test_db_get_all_attributes(self):
        fsp = "bdb:/home/adnan/DemoResults/03_Particles/mpi_proc_000#TcdA1-0010_frames_ptcls"
        atr = EMAN2db.db_get_all_attributes(fsp, ['nx'])
        attr_list = [{'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}]
        self.assertEqual(atr, attr_list)


    def test_db_get_all_attributes(self):
        fsp = "bdb:/home/adnan/DemoResults/06_SUBSTACK_ANO#isac_substack"
        atr = EMAN2db.db_get_all_attributes(fsp, 'xform.projection')
        print(atr)
        print(type(atr))
        print(type(atr[0]))

        # attr_list = [{'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}, {'nx': 352}]
        # self.assertEqual(atr, attr_list)


class Test_BDB(unittest.TestCase):

    def test_open_db(self):
        dic = EMAN2db.EMAN2DB()
        print('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname)
        dicc = dic.open_db(path= str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        self.assertEqual(dicc.path, str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))

    def test_bdb_init(self):
        dic = EMAN2db.EMAN2DB()
        dic.__init__(path= str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        self.assertEqual(dic.path, str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))

    def test_bdb_getitem_and_open_dict(self):
        ddb=EMAN2db.EMAN2DB.open_db('/home/adnan/PycharmProjects/DoseWeighting/Newfolder')
        ddb.open_dict(bdbname, ro=True)
        # print(ddb.__dict__)
        data = ddb.__getitem__(key = '20170629_00021_frameImage_ptcls')
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        self.assertTrue(numpy.array_equal(a.keys(), data.keys()))



class Test_Dbdict(unittest.TestCase):

    # def test_init_ro_True(self):
    #     a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
    #     print(a.keys())
    #     print(a.items())
    #     print(a.get('maxrec'))

    # def test_get_attrib(self):
    #     a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
    #     img =  a[0]
    #     print(img.get_attr('MRC.gamma'))
    #     print(a.get_header(0))
    #     print(a.get_attr(b'0', 'MRC.gamma'))

    def test_init(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)

        self.assertEqual(bdbname, a.name)
        self.assertEqual(bdb_path, a.path)
        self.assertEqual({}, a.key_translation_dict)

    # def test_updateold(self):
    #     bdb_path = '/home/adnan/DemoResults/03_PARTICLES/NEWSTACK/EMAN2DB'
    #     bdbname = 'stack'
    #     a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
    #     oldfile = '00image_counts.old'
    #     a.updateold(oldfile)

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
        self.assertEqual(a.__len__(), 210)


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

        real_keys =[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 'maxrec']

        self.assertEqual(keys, real_keys)

    def test_items_work(self):
        a = EMAN2db.DBDict(name = bdbname, path = bdb_path, ro=True)
        items = a.items()
        keys = a.keys()
        keys_in_items = [key for key,value in items]
        values = a.values()
        # print(values)
        self.assertEqual(keys, keys_in_items)


    def test_has_key(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        keys = a.keys()
        self.assertEqual(a.has_key(keys[-1]), True)

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


    def test_get_header_items_keys_values(self):
        a = EMAN2db.DBDict(name=bdbname, path=bdb_path, ro=True)
        header = a.get_header(0)
        header_orig_keys = ['MRC.nzstart', 'apix_x', 'MRC.gamma', 'xform.projection', 'changecount', 'maximum', 'data_path', 'MRC.beta', 'data_n', 'ctf_applied', 'minimum', 'ptcl_source_coord_id', 'mean_nonzero', 'MRC.label0', 'MRC.mapr', 'MRC.maps', 'ny', 'MRC.ispg', 'apix_y', 'MRC.xlen', 'HostEndian', 'nx', 'MRC.nsymbt', 'nz', 'MRC.mapc', 'ptcl_source_coord', 'sigma_nonzero', 'MRC.nxstart', 'ptcl_source_image', 'is_complex_ri', 'relion_max_prob_dist', 'MRC.nz', 'MRC.nx', 'MRC.ny', 'MRC.maximum', 'apix_z', 'MRC.nlabels', 'source_n', 'MRC.mean', 'origin_z', 'is_complex', 'is_complex_x', 'relion_norm_correct', 'chunk_id', 'ImageEndian', 'origin_y', 'origin_x', 'MRC.nystart', 'MRC.zlen', 'MRC.rms', 'resample_ratio', 'MRC.machinestamp', 'data_source', 'ptcl_source_apix', 'adnan_n', 'datatype', 'sigma', 'ptcl_source_relion', 'source_path', 'MRC.minimum', 'ctf', 'square_sum', 'MRC.mz', 'MRC.my', 'MRC.mx', 'MRC.alpha', 'MRC.ylen', 'mean']        # print(list(sorted(header.keys())) == list(sorted(header_orig_keys)))
        sb = '[(\'HostEndian\', \'little\'), (\'ImageEndian\', \'little\'), (\'MRC.alpha\', 90.0), (\'MRC.beta\', 90.0), (\'MRC.gamma\', 90.0), (\'MRC.ispg\', 0), (\'MRC.label0\', \'Relion    21-Oct-19  16:12:54\'), (\'MRC.machinestamp\', 16708), (\'MRC.mapc\', 1), (\'MRC.mapr\', 2), (\'MRC.maps\', 3), (\'MRC.maximum\', 4.168662071228027), (\'MRC.mean\', 0.0022823389153927565), (\'MRC.minimum\', -6.85170316696167), (\'MRC.mx\', 360), (\'MRC.my\', 360), (\'MRC.mz\', 1), (\'MRC.nlabels\', 1), (\'MRC.nsymbt\', 0), (\'MRC.nx\', 360), (\'MRC.nxstart\', 0), (\'MRC.ny\', 360), (\'MRC.nystart\', 0), (\'MRC.nz\', 1), (\'MRC.nzstart\', 0), (\'MRC.rms\', 1.0004879236221313), (\'MRC.xlen\', 318.6000061035156), (\'MRC.ylen\', 318.6000061035156), (\'MRC.zlen\', 1.0), (\'adnan_n\', 0), (\'apix_x\', 1.0), (\'apix_y\', 1.0), (\'apix_z\', 1.0), (\'changecount\', 0), (\'chunk_id\', 0), (\'ctf\', EMAN2Ctf().from_dict({"defocus":1.07294,"bfactor":0.0,"ampcont":10.00,"apix":0.885,"voltage":200.00,"cs":1.40}) ...), (\'ctf_applied\', 0), (\'data_n\', 0), (\'data_path\', \'/home/adnan/PycharmProjects/DoseWeighting/Newfolder/20170629_00021_frameImage_ptcls.mrcs\'), (\'data_source\', \'bdb:Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#20170629_00021_frameImage_ptcls\'), (\'datatype\', 7), (\'is_complex\', 0), (\'is_complex_ri\', 1), (\'is_complex_x\', 0), (\'maximum\', 3.379175901412964), (\'mean\', 0.0033558777067810297), (\'mean_nonzero\', 0.0033558777067810297), (\'minimum\', -6.456167697906494), (\'nx\', 360), (\'ny\', 360), (\'nz\', 1), (\'origin_x\', 0.0), (\'origin_y\', 0.0), (\'origin_z\', 0.0), (\'ptcl_source_apix\', 0.885), (\'ptcl_source_coord\', [2937, 3437]), (\'ptcl_source_coord_id\', 0), (\'ptcl_source_image\', \'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc\'), (\'ptcl_source_relion\', \'000001@Extract/job051/MOVIES_RELION/20170629_00021_frameImage.mrcs\'), (\'relion_max_prob_dist\', 0.390257), (\'relion_norm_correct\', 0.504526), (\'resample_ratio\', 1.0), (\'sigma\', 1.0007057189941406), (\'sigma_nonzero\', 1.0007057189941406), (\'source_n\', 0), (\'source_path\', \'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#sphire_relion_stack\'), (\'square_sum\', 129783.453125), (\'xform.projection\', Transform({\'az\':259.319,\'alt\':77.866,\'phi\':23.693,\'tx\':-1.89,\'ty\':-0.30,\'tz\':0.00,\'mirror\':0,\'scale\':1.0000,\'type\':\'eman\'}))]'
        val = ['-6.45616769791', '-6.85170316696', '/home/adnan/PycharmProjects/DoseWeighting/Newfolder/20170629_00021_frameImage_ptcls.mrcs', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0.0', '0.0', '0.0', '0.00228233891539', '0.00335587770678', '0.00335587770678', '0.390257', '0.504526', '0.885', '000001@Extract/job051/MOVIES_RELION/20170629_00021_frameImage.mrcs', '1', '1', '1', '1', '1', '1', '1.0', '1.0', '1.0', '1.0', '1.0', '1.00048792362', '1.00070571899', '1.00070571899', '129783.453125', '16708', '2', '3', '3.37917590141', '318.600006104', '318.600006104', '360', '360', '360', '360', '360', '360', '4.16866207123', '7', '90.0', '90.0', '90.0', 'EMAN2Ctf().from_dict({"defocus":1.07294,"bfactor":0.0,"ampcont":10.00,"apix":0.885,"voltage":200.00,"cs":1.40}) ...', 'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc', 'Relion    21-Oct-19  16:12:54', "Transform({'az':259.319,'alt':77.866,'phi':23.693,'tx':-1.89,'ty':-0.30,'tz':0.00,'mirror':0,'scale':1.0000,'type':'eman'})", '[2937, 3437]', 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#sphire_relion_stack', 'bdb:Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles#20170629_00021_frameImage_ptcls', 'little', 'little']
        # print(repr(list(sorted(map(str, header.values())))))
        # print(header.values())
        # print(repr(sorted(
        #     [entry if not entry.replace('.', '', 1).isdigit() else str(round(float(entry), 6)) for entry in
        #      map(str, header.values())])))
        self.assertEqual(list(sorted(header.keys())), list(sorted(header_orig_keys)))
        self.assertTrue(sb, repr(sorted(header.items())))
        self.assertTrue(repr(val), repr(list(sorted(map(str, header.values())))))


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
        print(a.get_header(0).values())
        data1 = a[a.keys()[0]]
        data2 = a[a.keys()[-3]]
        a.set_header(8, data1)
        # print(a[0])
        # print(a[a.keys()[0]])
        # self.assertTrue(numpy.array_equal(data1.get_2dview(),a[a.keys()[0]].get_2dview()))
        self.assertTrue(numpy.array_equal(data2.get_2dview(), data.get_2dview()))



    # def test_pickle2(self):
    #     import pickle
    #     ad = {'ptcl_source_coord_id': 5, 'MRC.nzstart': 7}
    #     print('1')
    #     print(pickle.dumps('MRC.nzstart', -1))
    #     print('2')
    #     print([key for key in ad.keys()])
    #     print('3')
    #     print([pickle.loads(pickle.dumps(x,-1)) for x in list(ad.keys())])
    #     print('4')
    #     print([pickle.dumps(x ,-1) for x in list(ad.keys())])
    #     py2bdbkeys = [b'\x80\x02U\x0bMRC.nzstartq\x00.', b'\x80\x02U\x14ptcl_source_coord_idq\x00.']
    #     py3newbdbkeys = [b'\x80\x04\x95\x18\x00\x00\x00\x00\x00\x00\x00\x8c\x14ptcl_source_coord_id\x94.', b'\x80\x04\x95\x0f\x00\x00\x00\x00\x00\x00\x00\x8c\x0bMRC.nzstart\x94.']
    #     for old_key in py2bdbkeys:
    #         print(old_key)
    #         print(pickle.loads(old_key))
    #     key_translation_dict = {}
    #     for old_key in py2bdbkeys:
    #         print(old_key)
    #         key_translation_dict[pickle.dumps(pickle.loads(old_key), -1)] = old_key
    #     print(key_translation_dict.keys())
    #     print(list(sorted(key_translation_dict.keys())) == list(sorted(py3newbdbkeys)))




class Test_Stacks(unittest.TestCase):
    def test_stack_substack(self):
        # img = EMAN2db.EMData()
        # EMAN2db.db_read_image(img, "bdb:/home/adnan/DemoResults/06_SUBSTACK_ANO#isac_substack",0 )
        # img1 = EMAN2db.db_read_images("bdb:/home/adnan/DemoResults/06_SUBSTACK_ANO#isac_substack")
        # print(img1[0].get_2dview())

        e = EMAN2_cppwrap.EMData()
        e.read_image("bdb:/home/adnan/DemoResults/04_ISAC#stack_ali2d", 5)
        d = EMAN2_cppwrap.EMData()
        d.read_image("bdb:/home/adnan/DemoResults/04_ISAC#case4stack", 5)
        img = EMAN2db.EMData.read_images("bdb:/home/adnan/DemoResults/04_ISAC#case4stack")
        self.assertTrue(numpy.array_equal(e.get_2dview(), d.get_2dview()))
        self.assertTrue(numpy.array_equal(img[5].get_2dview(), d.get_2dview()))


    def test_stack_dict_reading(self):
        a = EMAN2db.DBDict(name='case4stack', path="/home/adnan/DemoResults/04_ISAC/EMAN2DB", ro=True)
        # a.keys()
        # print(a.items())
        print(a.__getitem__(0))
        img = EMAN2db.EMData()
        img = EMAN2db.db_read_image(img,"bdb:/home/adnan/DemoResults/04_ISAC#case4stack",0)
        EMAN2db.db_close_dict(a)


    def test_stack_writing(self):
        img = EMAN2db.EMData()
        img.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))
        EMAN2db.db_close_dict(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + bdbname))

        newimg = EMAN2db.EMData()
        newimg = copy.deepcopy(img)
        # EMAN2.display(newimg)
        # print(newimg.get_2dview())
        newimg.write_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + 'test_img'))

    def test_read_new_file(self):
        imgnew = EMAN2db.EMData()
        imgnew.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB', '') + '#' + 'test_img'))
        EMAN2.display(imgnew)


    def test_read_stack_with_write(self):
        a = EMAN2db.EMData()
        a.read_data('/home/adnan/PycharmProjects/DoseWeighting/Newfolder/EMAN2DB/test_img_360x360x1', 0)

        print(a.get_3dview())

        # img = EMAN2db.EMData()
        # img.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + 'test_img'))
        # display(img)

        # img = EMAN2db.EMData()
        # img.read_image(str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + 'test_img'), 0)
        # # display(img)
        # a = EMAN2db.DBDict(name='test_img',
        #                    path=str('bdb:' + bdb_path.replace('/EMAN2DB','') + '#' + 'test_img'), ro=True)
        # # print(a.__dict__.keys())
        #
        # a.realopen()
        # print(a.keys())
        # print(a.__dict__.keys())
        # print(a)
        # EMAN2db.db_close_dict(a)


    def test_ctypes_work_ornot(self):
        import ctypes
        import numpy

        base_ptr = 139888466854280
        size = 40368464
        nimastack = 6989
        target_nx = 76

        ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
        buffer = numpy.frombuffer(ptr.contents, dtype="f4")

        # buffer = buffer.reshape(nimastack, target_nx, target_nx)
        print(buffer)
        print(buffer.shape)


