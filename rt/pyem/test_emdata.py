#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import testlib
import sys
import math
import os


class TestEMData(unittest.TestCase):
    

    def test_rotate_translate(self):
        img = TestUtil.get_debug_image("monomer.mrc")
        x=EMData()
        x.read_image(img)

        alt = 1.0329837512591338
        az = 3.7260642381912579
        phi = 5.7671541529246966

        x.rotate_translate(alt, az, phi, 12, 12, 12)
        testlib.check_emdata(x, sys.argv[0])

        x.rotate_translate(alt, az, phi, 16,16,16)
        testlib.check_emdata(x, sys.argv[0])

		
    def test_rotate_2d(self):
        img = TestUtil.get_debug_image("lattice.mrc")

        a = EMData()
        a.read_image(img)
        b=a.copy(0)
        b.rotate(0,0,math.pi/4)
        testlib.check_emdata(b, sys.argv[0])

        b=a.copy(0)
        b.rotate(0,0,math.pi/2)
        testlib.check_emdata(b, sys.argv[0])

    def test_rotate_3d(self):
        img = TestUtil.get_debug_image("3d.mrc")

        a = EMData()
        a.read_image(img)
        b=a.copy(0)
        b.rotate(0,0,math.pi/4)
        testlib.check_emdata(b, sys.argv[0])

        b=a.copy(0)
        b.rotate(0,0,math.pi/2)
        testlib.check_emdata(b, sys.argv[0])

    def test_project(self):        
        volume = EMData()
        volume.read_image(TestUtil.get_debug_image("groel3d.mrc"))
        pi = math.pi
        proj = volume.project("Standard", { "alt" : pi/3, "az" : pi/5, "phi" : 1})
        self.assertEqual(proj.get_xsize(), 100)
        self.assertEqual(proj.get_ysize(), 100)
        self.assertEqual(proj.get_zsize(), 1)
        testlib.check_emdata(proj, sys.argv[0])


    def test_calc_highest_locations(self):
        e = EMData()
        e.read_image(TestUtil.get_debug_image("search.dm3"))
        pixels = e.calc_highest_locations(1200)
        self.assertEqual(len(pixels),  612502)
        p = pixels[0]
        self.assertEqual(pixels[0], Pixel(776,677,0, 1201))


    def test_get_attr_dict(self):
        e = EMData()
        imgfile = "tablet.mrc"
        e.read_image(TestUtil.get_debug_image(imgfile))
        d0 = e.get_attr_dict()

        e += 10
        d1 = e.get_attr_dict()
        
        e += 100
        d2 = e.get_attr_dict()

        cur_attrlist = []
        
        for mydict in (d0, d1, d2):
            for mykey in mydict.keys():
                if (not ("MRC" in mykey)):
                    cur_attrlist.append(mykey + "=" + mydict[mykey].to_str()+"\n")


        attrfile = open(TestUtil.get_debug_image(imgfile+".attr"), "rb")
        old_attrlist = attrfile.readlines()        
        attrfile.close()

        self.assertEqual(cur_attrlist, old_attrlist)

    def test_get_rotated_clip(self):
        imagename = TestUtil.get_debug_image("monomer.mrc")
        a=EMData()
        a.read_image(imagename)
        b=a.get_rotated_clip(Transform([24,24,24], EULER_EMAN,0,0,0),[32,32,32],1.0)

        outfile = "test_get_rotated_clip_" + str(os.getpid()) + ".mrc"
        b.write_image(outfile)


def test_main():
    test_support.run_unittest(TestEMData)

if __name__ == '__main__':
    test_main()



