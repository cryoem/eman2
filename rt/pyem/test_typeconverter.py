#!/bin/env python

from EMAN2 import *
from testlib import *
import os
import sys
import Numeric

import unittest
from test import test_support
import testlib

class TestTypeConverter(unittest.TestCase):


    def test_emobject(self):
        num = TestUtil.get_debug_int(0)
        TestUtil.to_emobject({"int": num})

        fnum = TestUtil.get_debug_float(0)
        TestUtil.to_emobject({"float": fnum})

        lnum = long(num)
        TestUtil.to_emobject({"long": lnum})

        fl = get_list("float")
        TestUtil.to_emobject({"floatarray": fl})

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
        edict = get_dict("emobject")
        edict2 = TestUtil.test_dict(edict)
        self.assertEqual(edict, edict2)


    def test_point_size(self):
        nlist = get_list("int")
        flist = get_list("float")

        vec3i = TestUtil.test_Vec3i(nlist)
        self.assertEqual(vec3i.as_list(), nlist)

        vec3f = TestUtil.test_Vec3f(flist)
        self.assertEqual(vec3f.as_list(), flist)

        ip1 = TestUtil.test_IntPoint(nlist)
        self.assertEqual(list(ip1), nlist)

        fp1 = TestUtil.test_FloatPoint(flist)
        self.assertEqual(list(fp1), flist)

        is1 = TestUtil.test_IntSize(nlist)
        self.assertEqual(list(is1), nlist)

        fs1 = TestUtil.test_FloatSize(flist)
        self.assertEqual(list(fs1), flist)

    def test_map(self):

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

        emobjectmap = get_dict("emobject")
        emobjectmap2 = TestUtil.test_map_emobject(emobjectmap)
        self.assertEqual(emobjectmap, emobjectmap2)


    def test_vector(self):

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

        e1 = EMData()
        e1.read_image(TestUtil.get_debug_image("samesize1.mrc"))
        e2 = EMData()
        e2.set_size(10, 20, 5)
        e3 = EMData()

        elist = [e1, e2, e3]
        elist2 = TestUtil.test_vector_emdata(elist)
        testlib.check_emdata_list(elist, sys.argv[0])
        testlib.check_emdata_list(elist2, sys.argv[0])

        p1 = Pixel(1,2,3, 1.1)
        p2 = Pixel(4,5,6, 4.4)
        p3 = Pixel(7,8,9, 5.5)

        plist = [p1,p2,p3]
        plist2 = TestUtil.test_vector_pixel(plist)

        self.assertEqual(plist,plist2)


    def test_em2numpy(self):
        e = EMData()
        e.read_image(TestUtil.get_debug_image("groel2d.mrc"))
        nx = e.get_xsize()
        ny = e.get_ysize()

        a = EMNumPy.em2numpy(e)
        n = ny/2

        for i in range(nx):
            self.assertEqual(e.get_value_at(i, n), a[n][i])

        self.assertEqual(a.shape, (100,200))
        self.assertEqual(a.typecode(), "f")

        for x in range(nx):
            for y in range(n):
                a[y][x] = 0

        testlib.check_emdata(e, sys.argv[0])
        
        e2 = EMData()
        EMNumPy.numpy2em(a, e2)
        testlib.check_emdata(e2, sys.argv[0])


    def test_numpy2em(self):
        n = 100
        l = range(2*n*n)
        a = Numeric.reshape(Numeric.array(l, Numeric.Float32), (2*n, n))

        self.assertEqual(a.shape, (200, 100))
        self.assertEqual(a.typecode(), "f")

        e = EMData()
        EMNumPy.numpy2em(a, e)
        testlib.check_emdata(e, sys.argv[0])

        for i in range(n):
            self.assertEqual(e.get_value_at(i, 0), i)


    def test_Point_and_Size_class(self):        
        imagename = TestUtil.get_debug_image("monomer.mrc")
        img1 = EMData()
        img1.read_image(imagename)

        ptuple1 = (16,16,16)
        plist1 = list(ptuple1)
        
        img2=img1.get_rotated_clip(Transform(plist1, EULER_EMAN, 0,0,0), plist1, 1.0)
        img3=img1.get_rotated_clip(Transform(ptuple1, EULER_EMAN, 0,0,0), ptuple1, 1.0)
        
        testlib.check_emdata(img2, sys.argv[0])
        testlib.check_emdata(img3, sys.argv[0])



def test_main():
    test_support.run_unittest(TestTypeConverter)

if __name__ == '__main__':
    test_main()




