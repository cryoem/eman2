from __future__ import print_function
from __future__ import division

# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
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
import numpy
import EMAN2db
import EMAN2
from pickle import *
import sys
import os
absolute_import = (sys.version_info[0] >= 3)
import copy


"""
Input to this class is a bdb file with images
"""

bdb_path = "resources_tests/03_PARTICLES_BDB/mpi_proc_007/EMAN2DB"
bdbname = "TcdA1-0187_frames_ptcls"
bdbfile = 'bdb:'+os.path.join(os.getcwd(), "resources_tests/03_PARTICLES_BDB/mpi_proc_007/TcdA1-0187_frames_ptcls")

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
        self.assertEqual(len(listdic), 2)

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

    @unittest.skip('skipping ')
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




