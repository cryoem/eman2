#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import testlib
import sys
import math
import os


class TestEMData(unittest.TestCase):

    def test_multi_array_2d(self):
        nx = 16
        ny = 32
        infile = "test_multi_array_2d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny)

        e = EMData()
        e.read_image(infile)
        array = e.get_2dview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "f")
        self.assertEqual(array.shape, (ny, nx))

        for i in range(ny):
            for j in range(nx):
                self.assertEqual(e.get_value_at(j,i), array[i][j])
                

        os.unlink(infile)



    def test_multi_array_3d(self):
        nx = 8
        ny = 16
        nz = 4
        infile = "test_multi_array_3d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny, nz)

        e = EMData()
        e.read_image(infile)
        array = e.get_3dview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "f")
        self.assertEqual(array.shape, (nz, ny, nx))

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), array[i][j][k])
                

        os.unlink(infile)

        
    def test_multi_array_c2d(self):
        nx = 16
        ny = 16
        infile = "test_multi_array_c2d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny)

        e = EMData()
        e.read_image(infile)
        fft = e.do_fft()
        nx2 = fft.get_xsize()
        array = fft.get_2dcview()
        
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "F")
        self.assertEqual(array.shape, (ny, nx2/2))

        for i in range(ny):
            for j in range(0, nx2, 2):
                c1 = array[i][j/2]
                testlib.assertfloat(self, fft.get_value_at(j,i), c1.real)
                testlib.assertfloat(self, fft.get_value_at(j+1,i), c1.imag)
            
        os.unlink(infile)

    def test_multi_array_c3d(self):
        nx = 8
        ny = 16
        nz = 4
        infile = "test_multi_array_c3d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny, nz)

        e = EMData()
        e.read_image(infile)
        fft = e.do_fft()

        nx2 = fft.get_xsize()
        
        array = fft.get_3dcview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "F")
        self.assertEqual(array.shape, (nz, ny, nx2/2))

        for i in range(nz):
            for j in range(ny):
                for k in range(0,nx2,2):
                    c1 = array[i][j][k/2]
                    testlib.assertfloat(self, fft.get_value_at(k,j,i), c1.real)
                    testlib.assertfloat(self, fft.get_value_at(k+1,j,i), c1.imag)

        os.unlink(infile)


 
    def test_rotate_translate(self):
        infile = "test_rotate_translate.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 16,16,16)
        
        x=EMData()
        x.read_image(infile)

        alt = 1.0329837512591338
        az = 3.7260642381912579
        phi = 5.7671541529246966

        x.rotate_translate(az, alt, phi, 2, 2, 2)
        testlib.check_emdata(x, sys.argv[0])

        x.rotate_translate(az, alt, phi, 6, 6, 6)
        testlib.check_emdata(x, sys.argv[0])

        os.unlink(infile)
        
		
    def test_rotate_2d(self):
        infile = "test_rotate_2d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 24, 32)

        outfile1 = "test_rotate_2d_out_1.mrc"
        outfile2 = "test_rotate_2d_out_2.mrc"
        
        a = EMData()
        a.read_image(infile)
        b=a.copy(0)
        b.rotate(0, 0, math.pi/4)
        b.write_image(outfile1)
        
        # verify b
        
        b=a.copy(0)
        b.rotate(0, 0, math.pi/2)
        b.write_image(outfile2)

        # verify b

        os.unlink(infile)
        os.unlink(outfile1)
        os.unlink(outfile2)
        
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
        n = 20
        infile = "test_project.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, n, n, n)
        volume = EMData()
        volume.read_image(infile)
        pi = math.pi
        proj = volume.project("Standard", { "alt" : pi/3, "az" : pi/5, "phi" : 1})
        self.assertEqual(proj.get_xsize(), n)
        self.assertEqual(proj.get_ysize(), n)
        self.assertEqual(proj.get_zsize(), 1)
        testlib.check_emdata(proj, sys.argv[0])
        os.unlink(infile)

    def test_calc_highest_locations(self):
        infile = "test_calc_highest_locations.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 40, 60)
        
        e = EMData()
        e.read_image(infile)
        pixels = e.calc_highest_locations(400)
        self.assertEqual(len(pixels),  1144)
        self.assertEqual(pixels[0], Pixel(19, 10, 0, 401.0))

        os.unlink(infile)

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
                    cur_attrlist.append(mykey + "=" + str(mydict[mykey])+"\n")

        attrfile = open(TestUtil.get_debug_image(imgfile+".attr"), "rb")
        old_attrlist = attrfile.readlines()        
        attrfile.close()
        
        
    def test_get_rotated_clip(self):
        imagename = TestUtil.get_debug_image("monomer.mrc")
        a=EMData()
        a.read_image(imagename)
        b=a.get_rotated_clip(Transform([24,24,24], EULER_EMAN,0,0,0),[32,32,32],1.0)

        outfile = "test_get_rotated_clip_" + str(os.getpid()) + ".mrc"
        b.write_image(outfile)
        os.unlink(outfile)
        

	def test_get_clip(self):
		nx = 50
		ny = 50
		filebase = "test_get_clip_" + str(os.getpid())
		infile = filebase + ".mrc"
		TestUtil.make_image_file(infile, MRC, nx, ny, 2)
		
		e = EMData()
		e.read_image(infile)

		region1 = Region(nx/4, ny/4, nx/2, ny/2)
		outfile1 = filebase + "_out1.mrc"
		e2 = e.get_clip(region1)
		e2.write_image(outfile1)

		region2 = Region(-nx/4, -ny/4, 2*nx, 2*ny)
		outfile2 = filebase + "_out2.mrc"
		e3 = e.get_clip(region2)
		
	



def test_main():
    test_support.run_unittest(TestEMData)

if __name__ == '__main__':
    test_main()



