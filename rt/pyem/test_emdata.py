#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import testlib
import sys
import math
import os

class TestEMData(unittest.TestCase):

    def test_operator_unary(self):
        file1 = "test_operator_unary_1.mrc"
        e = EMData()
        e.set_size(100,200,1)
        e.to_one()
        e + 0.5
        e - 0.5
        e * 2
        e / 2.4
        12 + e
        12 - e
        4 * e
        4 / e

        e2 = EMData()
        e2.set_size(100,200,1)
        e2.to_one()

        e + e2
        e - e2
        e * e2
        e / e2
        
        e.write_image(file1)


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

        e *= 2
        for i in range(ny):
            for j in range(nx):
                self.assertEqual(e.get_value_at(j,i), array[i][j])

        array.savespace(1)
        array *= 2
        
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

        e *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), array[i][j][k])

        array.savespace(1)
        array *= 2
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


        e *= 2
        for i in range(ny):
            for j in range(0, nx2, 2):
                c1 = array[i][j/2]
                testlib.assertfloat(self, fft.get_value_at(j,i), c1.real)
                testlib.assertfloat(self, fft.get_value_at(j+1,i), c1.imag)

        array.savespace(1)
        array *= 2
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

        e *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(0,nx2,2):
                    c1 = array[i][j][k/2]
                    testlib.assertfloat(self, fft.get_value_at(k,j,i), c1.real)
                    testlib.assertfloat(self, fft.get_value_at(k+1,j,i), c1.imag)


        array.savespace(1)
        array *= 2
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
        
    def test_get_clip(self):
        nx = 32
        ny = 48
        filebase = "test_get_clip_" + str(os.getpid())
        infile = filebase + ".mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny)

        e = EMData()
        e.read_image(infile)

        region1 = Region(nx/4, ny/4, nx/2, ny/2)
        outfile1 = filebase + "_out1.mrc"
        e2 = e.get_clip(region1)
        e2.write_image(outfile1)

        self.assertEqual(e2.get_xsize(), nx/2)
        self.assertEqual(e2.get_ysize(), ny/2)
        
        region2 = Region(-nx/4, -ny/4, 2*nx, 2*ny)
        outfile2 = filebase + "_out2.mrc"
        e3 = e.get_clip(region2)
        e3.write_image(outfile2)
        
        self.assertEqual(e3.get_xsize(), nx*2)
        self.assertEqual(e3.get_ysize(), ny*2)

        os.unlink(infile)
        os.unlink(outfile1)
        os.unlink(outfile2)

        
    def test_get_rotated_clip(self):
        imagename = TestUtil.get_debug_image("monomer.mrc")
        a=EMData()
        a.read_image(imagename)
        b=a.get_rotated_clip(Transform3D([24,24,24], 0,0,0),[32,32,32],1.0)

        outfile = "test_get_rotated_clip_" + str(os.getpid()) + ".mrc"
        b.write_image(outfile)
        os.unlink(outfile)


    def test_complex_image(self):
        nx = 16
        ny = 32
        infile = "test_complex_image_1.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT_COMPLEX, nx, ny)

        e = EMData()
        e.read_image(infile)
        attrd = e.get_attr_dict()

        self.assertEqual(attrd["nx"], nx+2)
        self.assertEqual(attrd["ny"], ny)
        self.assertEqual(attrd["nz"], 1)
        self.assertEqual(attrd["datatype"], EM_FLOAT_COMPLEX)

        os.unlink(infile)

        
    def test_set_value_at(self):
        nx = 10
        ny = 20
        nz = 2
        infile = "test_set_value_at_in.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny, nz)
        
        e = EMData()
        e.read_image(infile)

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    e.set_value_at(k,j,i, 1)

        os.unlink(infile)
        narray = EMNumPy.em2numpy(e)

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), 1)
                    self.assertEqual(narray[i][j][k], 1)
        
    def test_to_one(self):
        nx = 12
        ny = 24
        nz = 2
        file1 = "test_to_one.mrc"
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny, nz)

        e1 = EMData()
        e1.read_image(file1)
        e1.to_one()

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e1.get_value_at(k,j,i), 1)
        
        os.unlink(file1)
        
    def test_to_zero(self):
        nx = 12
        ny = 24
        nz = 2
        file1 = "test_to_zero.mrc"
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny, nz)

        e1 = EMData()
        e1.read_image(file1)
        e1.to_zero()

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e1.get_value_at(k,j,i), 0)
        
        os.unlink(file1)

    def test_calc_radial_dist(self):
        file1 = "test_calc_radial_dist.mrc"
        TestUtil.make_image_file(file1, MRC)
        
        e1 = EMData()
        e1.read_image(file1)
        ny = e1.get_ysize()
        dlist = e1.calc_radial_dist(ny, 0, 0.5)

        dlist2 = [0.0, 0.5732, 1.1464, 2.0, 4.5135, 5.7996, 9.2307, 12.1255, 16.4826, 19.1460, 25.6796, 30.2887, 36.4654, 41.4641, 49.9667, 1538.3833]

        for i in range(len(dlist)):
            testlib.assertfloat(self, dlist[i], dlist2[i])

        os.unlink(file1)


    def test_stat_locations(self):
        nx = 16
        ny = 24
        file1 = "test_stat_locations.mrc"
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny)
        
        e1 = EMData()
        e1.read_image(file1)

        max_index = e1.calc_max_index()
        self.assertEqual(max_index, 0)

        max_location = e1.calc_max_location()
        self.assertEqual(max_location, (0,0,0))

        max_index2 = max_location[2] * nx * ny + max_location[1] * nx + max_location[0]
        self.assertEqual(max_index, max_index2)
        
        attrdict = e1.get_attr_dict()
        max = attrdict["maximum"]
        self.assertEqual(max, 208)

        max2 = e1.get_value_at(max_location[2], max_location[1], max_location[0])
        self.assertEqual(max, max2)
        
        os.unlink(file1)


    def test_image_overwrite(self):
        file1 = "test_image_overwrite.mrc"
        nx = 16
        ny = 18
        nz = 2
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny, nz)
        e = EMData()
        e.read_image(file1)
        
        nx2 = nx*2
        ny2 = ny*2
        nz2 = nz
        e.set_size(nx2, ny2, nz2)
        e.write_image(file1)
        
        e2 = EMData()
        e2.read_image(file1)
        self.assertEqual(nx2, e2.get_xsize())
        self.assertEqual(ny2, e2.get_ysize())
        self.assertEqual(nz2, e2.get_zsize())

        os.unlink(file1)
        


    def test_ctf(self):
        infile = "test_ctf_in.mrc"
        TestUtil.make_image_file(infile, MRC)
        ctf = SimpleCtf()
        d = {"defocus":1, "bfactor":2, "amplitude":3, "ampcont":4, "noise1":5, "noise2":6, "noise3":7, "noise4":8, "voltage":9, "cs":10,"apix":11}
        
        ctf.from_dict(d)

        e = EMData()
        e.read_image(infile)
        e.set_ctf(ctf)

        outfile = "test_ctf_out.mrc"
        e.write_image(outfile)

        e2 = EMData()
        e2.read_image(outfile)
        ctf2 = e2.get_ctf()

        ctfstr = ctf2.to_string()
        self.assertEqual(ctfstr, "1 2 3 4 5 6 7 8 9 10 11")
        os.unlink(infile)
        os.unlink(outfile)
        

        
def test_main():
    test_support.run_unittest(TestEMData)

if __name__ == '__main__':
    test_main()



