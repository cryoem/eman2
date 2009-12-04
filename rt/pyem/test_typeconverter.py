#!/usr/bin/env python

#
# Author: Liwei Peng, 01/30/2005 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from testlib import *
import os
import sys
import numpy

import unittest
import testlib
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestTypeConverter(unittest.TestCase):
    """this is test call those routines in TestUtil class"""

    def test_emobject_to_py(self):
        """test emobject to Python type ....................."""
        nx = 10
        ny = 12
        nz = 2
        img1 = EMData()
        img1.set_size(nx, ny, nz)        
        #can not call thi sfunciton here, Error:
        #*** glibc detected *** double free or corruption (!prev): 0x09257190 ***
#        img2 = TestUtil.emobject_to_py(img1)
#        self.assertEqual(img2.get_xsize(), nx)
#        self.assertEqual(img2.get_ysize(), ny)
#        self.assertEqual(img2.get_zsize(), nz)

        attr_dict = img1.get_attr_dict()
        self.assertEqual(type(attr_dict["minimum"]), type(2.2))
        self.assertEqual(type(attr_dict["nx"]), type(nx))
        
        b1 = True
        b2 = TestUtil.emobject_to_py(b1)
        self.assertEqual(b1,b2)
        
        n1 = 100
        n2 = TestUtil.emobject_to_py(n1)
        self.assertEqual(n1, n2)

        f1 = 3.14
        f2 = TestUtil.emobject_to_py(f1)
        self.assertEqual(f1, f2)

        str1 = "helloworld"
        str2 = TestUtil.emobject_to_py(str1)
        self.assertEqual(str1, str2)

        farray = TestUtil.emobject_farray_to_py()
        farray2 = testlib.get_list("float")
        self.assertEqual(farray, farray2)

        strarray = TestUtil.emobject_strarray_to_py()
        strarray2 = testlib.get_list("string")
        self.assertEqual(strarray, strarray2)

        alt = 1.0329837512591338
        az = 3.7260642381912579
        phi = 5.7671541529246966
        t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
        t2 = TestUtil.emobject_to_py(t1)
        self.assertEqual(t1.get_matrix(), t2.get_matrix())
        
        ctf1 = EMAN1Ctf()
        ctf1.from_vector((1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0))
        ctf2 = TestUtil.emobject_to_py(ctf1)
        self.assertEqual(ctf1.to_string(), ctf2.to_string())

        testfile = "test_emobject_to_py_xydata.txt"
        out = open(testfile, "wb")
        for f in farray2:
            out.write(str(f) + " " + str(f) + "\n")
        out.close()
        xyd = XYData()
        xyd.read_file(testfile)
        xyd2 = TestUtil.emobject_to_py(xyd)
        self.assertEqual(xyd2.get_size(), len(farray2))
        for i in range(len(farray2)):
            self.assertEqual(xyd2.get_x(i), farray2[i])
            self.assertEqual(xyd2.get_y(i), farray2[i])
        os.unlink(testfile)        
   
    def test_emobject(self):
        """test emobject ...................................."""
        num = TestUtil.get_debug_int(0)
        TestUtil.to_emobject({"int": num})

        fnum = TestUtil.get_debug_float(0)
        TestUtil.to_emobject({"float": fnum})

        lnum = long(num)
        TestUtil.to_emobject({"long": lnum})

        fl = get_list("float")
        TestUtil.to_emobject({"floatarray": fl})

        str1 = TestUtil.get_debug_string(0)
        TestUtil.to_emobject({"string": str1})

        e = EMData()
        nx = TestUtil.get_debug_int(0)
        ny = TestUtil.get_debug_int(1)
        nz = TestUtil.get_debug_int(2)
        e.set_size(nx, ny, nz)
        TestUtil.to_emobject({"emdata": e})

        xyd = XYData()
        testfile = "xydata.txt"
        out = open(testfile, "wb")
        for f in fl:
                out.write(str(f) + " " + str(f) + "\n")
        out.close()

        xyd.read_file(testfile)
        TestUtil.to_emobject({"xydata" : xyd})
        os.unlink(testfile)

        strlist = get_list("string")
        TestUtil.to_emobject({"stringarray":strlist})

    def test_Dict(self):
        """test Dict class .................................."""
        edict = get_dict("float")
        edict2 = TestUtil.test_dict(edict)        
        self.assertEqual(edict, edict2)

    def test_point_size(self):
        """test vector and point size ......................."""
        nlist = get_list("int")
        flist = get_list("float")

        vec3i = TestUtil.test_Vec3i(nlist)
        self.assertEqual(vec3i, nlist)

        vec3f = TestUtil.test_Vec3f(flist)
        self.assertEqual(vec3f, flist)

        ip1 = TestUtil.test_IntPoint(nlist)
        self.assertEqual(list(ip1), nlist)

        fp1 = TestUtil.test_FloatPoint(flist)
        self.assertEqual(list(fp1), flist)

        is1 = TestUtil.test_IntSize(nlist)
        self.assertEqual(list(is1), nlist)

        fs1 = TestUtil.test_FloatSize(flist)
        self.assertEqual(list(fs1), flist)

    def test_map(self):
        """test map ........................................."""
        imap = get_dict("int")
        imap2 = TestUtil.test_map_int(imap)
        self.assertEqual(imap, imap2)

        lmap = get_dict("long")
        lmap2 = TestUtil.test_map_long(lmap)
        self.assertEqual(lmap, lmap2)

        fmap = get_dict("float")
        fmap2 = TestUtil.test_map_float(fmap)
        self.assertEqual(fmap, fmap2)

        smap = get_dict("string")
        smap2 = TestUtil.test_map_string(smap)
        self.assertEqual(smap, smap2)

    def test_vector(self):
        """test vector ......................................"""
        nlist = get_list("int")
        flist = get_list("float")
        llist = get_list("long")
        slist = get_list("string")

        nlist2 = TestUtil.test_vector_int(nlist)
        self.assertEqual(nlist, nlist2)

        flist2 = TestUtil.test_vector_float(flist)
        self.assertEqual(flist, flist2)

        llist2 = TestUtil.test_vector_long(llist)
        self.assertEqual(llist, llist2)

        slist2 = TestUtil.test_vector_string(slist)
        self.assertEqual(slist, slist2)

        imgfile1 = "test_vector1.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        e1 = EMData()
        e1.read_image(imgfile1)
        e2 = EMData()
        e2.set_size(10, 20, 5)
        e3 = EMData()
        e3.set_size(10, 20, 5)

        elist = [e1, e2, e3]
        testlib.check_emdata_list(elist, sys.argv[0])
        #can not call this routine, Error:
        #*** glibc detected *** double free or corruption (!prev): 0x0916b8e8 ***
#        elist2 = TestUtil.test_vector_emdata(elist)
#        testlib.check_emdata_list(elist2, sys.argv[0])

        os.unlink(imgfile1)
        
        p1 = Pixel(1,2,3, 1.1)
        p2 = Pixel(4,5,6, 4.4)
        p3 = Pixel(7,8,9, 5.5)

        plist = [p1,p2,p3]
        plist2 = TestUtil.test_vector_pixel(plist)
        self.assertEqual(plist,plist2)

    def test_em2numpy(self):
        """test em2numpy ...................................."""
        imgfile1 = "test_em2numpy_1.mrc"
        nx0 = 100
        ny0 = 200
        TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, nx0, ny0)
        
        e = EMData()
        e.read_image(imgfile1)
        nx = e.get_xsize()
        ny = e.get_ysize()

        a = EMNumPy.em2numpy(e)
        n = ny/2

        for i in range(nx):
            self.assertEqual(e.get_value_at(i, n), a[n][i])

        self.assertEqual(a.shape, (ny0, nx0))
        self.assertEqual(a.dtype, "f")

        for x in range(nx):
            for y in range(n):
                a[y][x] = 0

        testlib.check_emdata(e, sys.argv[0])
        
        b = a.copy()
        e2 = EMNumPy.numpy2em(b)
        testlib.check_emdata(e2, sys.argv[0])

        os.unlink(imgfile1)

    def test_numpy2em(self):
        """test numpy2em .................................... """
        n = 100
        l = range(2*n*n)
        a = numpy.reshape(numpy.array(l, numpy.float32), (2*n, n))

        self.assertEqual(a.shape, (2*n, n))
        self.assertEqual(a.dtype, "float32")

        e = EMNumPy.numpy2em(a)
        testlib.check_emdata(e, sys.argv[0])

        for i in range(n):
            self.assertEqual(e.get_value_at(i, 0), i)
         
        #test the float64 data(Python float) convert to float32(EMData float)     
    def test_numpy2em2(self):
        """test numpy2em for float64 ........................"""
        n1 = numpy.random.ranf( (5,3) )
        #e1 = EMData()
        e1 = EMNumPy.numpy2em(n1)
        n2 = EMNumPy.em2numpy(e1)
        diff = numpy.max(numpy.max(n2 - n1))
        self.assertAlmostEqual(diff, 0, 3)
        
    def test_em2numpy2(self):
        """test em2numpy again .............................."""
        imgfile1 = "test_em2numpy2_1.mrc"
        nx0 = 100
        ny0 = 200
        TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, nx0, ny0)
        
        e = EMData()
        e.read_image(imgfile1)
        nx = e.get_xsize()
        ny = e.get_ysize()

        a = EMNumPy.em2numpy(e)
        os.unlink(imgfile1)
        
    def test_Point_and_Size_class(self):
        """test point and Size class ........................"""
        imgfile1 = "test_Point_and_Size_class_1.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, 32,32,32)

        img1 = EMData()
        img1.read_image(imgfile1)

        ptuple1 = (16,16,16)
        plist1 = list(ptuple1)
        
        img2=img1.get_rotated_clip(Transform(), plist1, 1.0)
        img3=img1.get_rotated_clip(Transform(), ptuple1, 1.0)
        
        testlib.check_emdata(img2, sys.argv[0])
        testlib.check_emdata(img3, sys.argv[0])

        os.unlink(imgfile1)
    
    def test_None_emobject(self):
        """test None as EMObject ............................"""
        e = EMData(32,32)
        e.set_attr('Nothing', None)
        self.assertEqual(e.get_attr('Nothing'), None)

def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTypeConverter)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
