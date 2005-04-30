#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import os
import sys
import testlib
import os


class TestSpiderIO(unittest.TestCase):


    def test_make_spider(self):
        file1 = "test_make_spider_1.spi"
        nx1 = 100
        ny1 = 24
        TestUtil.make_image_file(file1, SPIDER, EM_FLOAT, nx1, nx1)
        err = TestUtil.verify_image_file(file1, SPIDER, EM_FLOAT, nx1, nx1)
        self.assertEqual(err, 0)
        os.unlink(file1)

    def test_overwrite_spider(self):
        file1 = "test_overwrite_spider.spi"
        nx1 = 24
        ny1 = 32
        TestUtil.make_image_file(file1, SPIDER, EM_FLOAT, nx1, ny1)
        TestUtil.make_image_file(file1, SPIDER, EM_FLOAT, nx1*2, ny1*2)
        #os.unlink(file1)
        

class TestHdfIO(unittest.TestCase):

    def test_make_image(self):
        imgfile = "test_make_image_1.h5"
        TestUtil.make_image_file(imgfile, HDF, EM_FLOAT)
        err = TestUtil.verify_image_file(imgfile, HDF, EM_FLOAT)
        self.assertEqual(err, 0)

        os.unlink(imgfile)
        
    def test_read_image(self):
        nx = 20
        ny = 30
        nz = 2
        imgfile = "test_read_image.h5"
        TestUtil.make_image_file(imgfile, HDF, EM_FLOAT, nx, ny, nz)

        e = EMData()
        e.read_image(imgfile)
        attrdict = e.get_attr_dict()
        
        self.assertEqual(attrdict["nx"], nx)
        self.assertEqual(attrdict["ny"], ny)
        self.assertEqual(attrdict["nz"], nz)
        self.assertEqual(attrdict["ImageEndian"], "big")
        self.assertEqual(attrdict["datatype"], EM_FLOAT)
        self.assertEqual(attrdict["is_complex"], 0)
        self.assertEqual(attrdict["maximum"], 325.0)
        self.assertEqual(attrdict["minimum"], 0.0)

        os.unlink(imgfile)
        

    def test_hdf_attr(self):

        infile = "test_hdf_attr_in.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT)

        ctf = SimpleCtf()
        d = {"defocus":1, "bfactor":2}
        ctf.from_dict(d)

        e = EMData()
        e.read_image(infile)
        az = 1.5
        alt = 2.5
        phi = 0.5        
        e.set_rotation(az, alt, phi)
        e.set_ctf(ctf)
        
        outfile = "test_hdf_attr_out_1.h5"
        e.write_image(outfile, 0, HDF)
        
        e2 = EMData()
        e2.read_image(outfile)
        attrdict2 = e2.get_attr_dict()

        self.assertEqual(attrdict2["orientation_convention"], "EMAN")
        self.assertEqual(attrdict2["euler_az"], az)
        self.assertEqual(attrdict2["euler_alt"], alt)
        self.assertEqual(attrdict2["euler_phi"], phi)

        theta = 10.0
        phi = 20.0
        omega = 30.0

        e2.set_attr("orientation_convention", "MRC")
        e2.set_attr("euler_theta", theta)
        e2.set_attr("euler_phi", phi)
        e2.set_attr("euler_omega", omega)

        outfile2 = "test_hdf_attr_out_2.h5"
        e2.write_image(outfile2, 0, HDF)

        e3 = EMData()
        e3.read_image(outfile2)
        attrdict3 = e3.get_attr_dict()
        
        self.assertEqual(attrdict3["orientation_convention"], "MRC")
        self.assertEqual(attrdict3["euler_theta"], theta)
        self.assertEqual(attrdict3["euler_omega"], omega)
        self.assertEqual(attrdict3["euler_phi"], phi)
        
        os.unlink(infile)
        os.unlink(outfile)
        os.unlink(outfile2)

    def test_write_image(self):
        imgfile1 = "test_write_image_1.mrc"
        imgfile2 = "test_write_image_2.mrc"
        imgfile3 = "test_write_image_3.mrc"

        TestUtil.make_image_file(imgfile1, MRC)
        TestUtil.make_image_file(imgfile2, MRC)
        TestUtil.make_image_file(imgfile3, MRC)

        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = EMData()
        e2.read_image(imgfile2)

        e3 = EMData()
        e3.read_image(imgfile3)

        outfile = "test_write_image_out.h5"
        e1.write_image(outfile, 0)
        e2.write_image(outfile, 1)
        e3.write_image(outfile, 2)
        
        os.unlink(imgfile1)
        os.unlink(imgfile2)
        os.unlink(imgfile3)

        os.unlink(outfile)
        

class TestMrcIO(unittest.TestCase):
    
    def test_overwrite(self):
        base = "overwrite_" + str(os.getpid())
        imgfile1 = base + "_1.mrc"
        TestUtil.make_image_file(imgfile1, MRC, EM_FLOAT, 10, 20, 1)
        TestUtil.make_image_file(imgfile1, MRC, EM_FLOAT, 30, 40, 1)
        e = EMData()
        e.read_image(imgfile1)
        self.assertEqual(TestUtil.verify_image_file(imgfile1, MRC, EM_FLOAT, 30, 40, 1), 0)
        
        imgfile2 = base + "_2.mrc"
        TestUtil.make_image_file(imgfile2, MRC, EM_FLOAT, 30, 40, 1)
        TestUtil.make_image_file(imgfile2, MRC, EM_FLOAT, 10, 20, 1)
        self.assertEqual(TestUtil.verify_image_file(imgfile2, MRC, EM_FLOAT, 10, 20, 1), 0)
        
        imgfile3 = base + "_3.mrc"
        TestUtil.make_image_file(imgfile3, MRC, EM_FLOAT, 30, 40, 1)
        TestUtil.make_image_file2(imgfile3, MRC, EM_FLOAT, 30, 40, 1)
        self.assertEqual(TestUtil.verify_image_file2(imgfile3, MRC, EM_FLOAT, 30, 40, 1), 0)

        os.unlink(imgfile1)
        os.unlink(imgfile2)
        os.unlink(imgfile3)
        

    def test_make_image_file(self):
        base = "test_make_image_file"
        img1 = base + "_c.mrc"
        img2 = base + "_s.mrc"
        img3 = base + "_f.mrc"
        img4 = base + "_sc.mrc"
        img5 = base + "_fc.mrc"
        
        TestUtil.make_image_file(img1, MRC, EM_UCHAR)
        self.assertEqual(TestUtil.verify_image_file(img1, MRC, EM_UCHAR), 0)
  
        TestUtil.make_image_file(img2, MRC, EM_USHORT)
        self.assertEqual(TestUtil.verify_image_file(img2, MRC, EM_USHORT), 0)
        
        TestUtil.make_image_file(img3, MRC, EM_FLOAT,64,64)
        self.assertEqual(TestUtil.verify_image_file(img3, MRC, EM_FLOAT, 64,64), 0)
        
        TestUtil.make_image_file(img4, MRC, EM_USHORT_COMPLEX)
        self.assertEqual(TestUtil.verify_image_file(img4, MRC, EM_USHORT_COMPLEX), 0)
        
        TestUtil.make_image_file(img5, MRC, EM_FLOAT_COMPLEX)
        self.assertEqual(TestUtil.verify_image_file(img5, MRC, EM_FLOAT_COMPLEX), 0)

        os.unlink(img1)
        os.unlink(img2)
        os.unlink(img3)
        os.unlink(img4)
        os.unlink(img5)


        img6 = base + "_3d.mrc"
        TestUtil.make_image_file(img6, MRC, EM_FLOAT, 16,16,10)
        os.unlink(img6)



    def test_complex_image(self):
        imgfile1 = "test_complex_image.mrc"
        TestUtil.make_image_file(imgfile1, MRC, EM_FLOAT)

        imgfile2 = "test_complex_image_fft.mrc"

        e = EMData()
        e.read_image(imgfile1)
        fft = e.do_fft()
        fft.write_image(imgfile2)


        self.assertEqual(TestUtil.verify_image_file(imgfile2, MRC, EM_FLOAT_COMPLEX,
                                                    e.get_xsize(), e.get_ysize(),
                                                    e.get_zsize()), 0)
        
        os.unlink(imgfile1)
        os.unlink(imgfile2)


   
    def test_mrcio_label(self):
        pid = str(os.getpid())
        infile = "test_mrcio_label_in_" + pid + ".mrc"
        TestUtil.make_image_file(infile, MRC)

        e = EMData()
        e.read_image(infile)
        label = "ByLiweiPeng"
        label_i = 3
        labelname = "MRC.label" + str(label_i)
        
        e.set_attr(labelname, label)

        outfile="test_mrcio_label_out_" + pid + ".mrc"        
        e.write_image(outfile, 0, MRC)
        
        e2 = EMData()
        e2.read_image(outfile)
        d = e2.get_attr_dict()
        nlabels = int(d["MRC.nlabels"])

        os.unlink(outfile)
        os.unlink(infile)
        
        self.assert_(nlabels > label_i)
        self.assertEqual(d[labelname], label)


class TestImagicIO(unittest.TestCase):
    
    def test_no_ext_filename(self):        
        infile = "test_no_ext_filename.mrc"
        TestUtil.make_image_file(infile, MRC)

        outfile = "test_no_ext_filename_out"
        
        img = EMData()
        img.read_image(infile)
        img.write_image(outfile, 0, IMAGIC)

        (hedfile, imgfile) = testlib.get_imagic_filename_pair(outfile)
        os.unlink(infile)
        os.unlink(hedfile)
        os.unlink(imgfile)

    def test_append_image(self):
        infile = "test_append_image.hed"
        TestUtil.make_image_file(infile, IMAGIC, EM_FLOAT, 16, 16, 10)
               
        e = EMData()
        all_imgs = e.read_images(infile)

        outfile1 = "test_append_image_out_" + str(os.getpid()) + ".hed"
        
        for img in all_imgs:
            img.append_image(outfile1)
        

        (infilehed, infileimg) = testlib.get_imagic_filename_pair(infile)
        (outfilehed, outfileimg) = testlib.get_imagic_filename_pair(outfile1)

        os.unlink(infilehed)
        os.unlink(infileimg)
        os.unlink(outfilehed)
        os.unlink(outfileimg)


    def test_append_to_newfile(self):
        infile = "test_append_to_newfile_in.mrc"
        outfile = "test_append_to_newfile_in.img"
        
        TestUtil.make_image_file(infile, MRC)
        e = EMData()
        e.read_image(infile)
        e.append_image(outfile, IMAGIC)

        # check e

        (outhed, outimg) = testlib.get_imagic_filename_pair(outfile)

        os.unlink(infile)
        os.unlink(outhed)
        os.unlink(outimg)

    def test_append_to_existing_file(self):
        img1 = "test_append_to_existing_file_1.hed"
        TestUtil.make_image_file(img1, IMAGIC)
        e = EMData()
        e.read_image(img1, 0, False, None, True)
        e.append_image(img1, IMAGIC)

        # verify here
        
        (hedfile, imgfile) = testlib.get_imagic_filename_pair(img1)
        os.unlink(hedfile)
        os.unlink(imgfile)

        
    def test_insert_to_newfile(self):
        img1 = "test_insert_to_newfile_in.hed"
        TestUtil.make_image_file(img1, IMAGIC)
        e = EMData()
        e.read_image(img1)
        outfile = "test_insert_to_newfile_out.hed"
        nimg = 4
        e.write_image(outfile, nimg-1, IMAGIC)

        nimg2 = EMUtil.get_image_count(outfile)
        self.assertEqual(nimg2, nimg)
        
        (hedfile1, imgfile1) = testlib.get_imagic_filename_pair(img1)
        (hedfile2, imgfile2) = testlib.get_imagic_filename_pair(outfile)
        
        os.unlink(hedfile1)
        os.unlink(imgfile1)
        os.unlink(hedfile2)
        os.unlink(imgfile2)


    def test_insert_beyond_existing_file(self):
        infile = "insert_beyond_existing_in.hed"
        TestUtil.make_image_file(infile, IMAGIC)
        e = EMData()
        e.read_image(infile, 0, False, None, True)

        nimg1 = EMUtil.get_image_count(infile)
        self.assertEqual(nimg1, 1)

        n2 = 9
        e.write_image(infile, n2, IMAGIC)
        nimg2 = EMUtil.get_image_count(infile)
        self.assertEqual(nimg2, n2+1)

        # todo: verify images

        n3 = 14
        e.write_image(infile, n3, IMAGIC)
        nimg3 = EMUtil.get_image_count(infile)
        self.assertEqual(nimg3, n3+1)

        # todo: verify images
    
        (hedfile, imgfile) = testlib.get_imagic_filename_pair(infile)
        os.unlink(hedfile)
        os.unlink(imgfile)

    def test_insert_inside_existing_file(self):
        infile = "test_insert_inside_existing_file_1.img"
        TestUtil.make_image_file(infile, IMAGIC, EM_FLOAT, 20, 30, 20)
        
        insertfile = "test_insert_inside_existing_file_2.mrc"
        TestUtil.make_image_file(insertfile, MRC, EM_FLOAT, 20, 30)
        e = EMData()
        e.read_image(insertfile)
        e.write_image(infile, 2, IMAGIC)

        # verify result
        
        (hedfile1, imgfile1) = testlib.get_imagic_filename_pair(infile)
        os.unlink(hedfile1)
        os.unlink(imgfile1)
        os.unlink(insertfile)


class TestEmIO(unittest.TestCase):        

    def test_emio_write(self):
        infile = "test_emio_write_in.mrc"
        TestUtil.make_image_file(infile, MRC)
        e = EMData()
        e.read_image(infile)

        outfile = "test_emio_write_out.em"
        e.write_image(outfile, 0, EM)
        TestUtil.check_image(outfile)

        os.unlink(infile)
        os.unlink(outfile)
        testlib.unlink_data_header_files(outfile)
        

class TestImageIO(unittest.TestCase):


        
    def test_emdata_overwriting(self):
        imgfile1 = "test_emdata_overwriting_1.mrc"
        imgfile2 = "test_emdata_overwriting_2.mrc"
        TestUtil.make_image_file(imgfile1, MRC, EM_FLOAT, 16,32,2)
        TestUtil.make_image_file(imgfile2, MRC, EM_FLOAT, 32, 24)

        a=EMData()
        b=EMData()
        
        a.read_image(imgfile1, 0)
        b.read_image(imgfile2, 0)
        b.read_image(imgfile1, 0)

        os.unlink(imgfile1)
        os.unlink(imgfile2)


    def test_image_overwriting(self):
        e = EMData()
        e.set_size(20, 20, 1)
        outfile = "test_image_overwriting_" + str(os.getpid()) + ".mrc"
        e.write_image(outfile)
        e.write_image(outfile)
        os.unlink(outfile)
        
    def create_dummy_region(self,e):
        nx = e.get_xsize()
        ny = e.get_ysize()
        nz = e.get_zsize()

        x0 = nx/4
        y0 = ny/4
        z0 = nz/4
        
        xsize = nx/2
        ysize = ny/2
        zsize = nz/2
        
        if zsize == 0:
            zsize = 1
            
        ndims = e.get_ndim()

        if ndims == 2:
            region_2d = Region(x0, y0, xsize, ysize)
            region_3d = None
        elif ndims == 3:
            region_2d = Region(x0, y0, z0, xsize, ysize, 1)
            region_3d = Region(x0, y0, z0, xsize, ysize, zsize)

        return (region_2d, region_3d)
    

    def region_read_test(self, imgtype, imgfile, outtype = None):
         if not outtype:
            outtype = imgtype
            
         is_3d = False
         if imgtype == IMAGIC:
            is_3d = True

         imgbase = Util.remove_filename_ext(imgfile)
         ext = Util.get_filename_ext(imgfile)

         readfile_2d = imgbase + "_read_region_2d." + ext
         readfile_3d = imgbase + "_read_region_3d." + ext

         e = EMData()
         e.read_image(imgfile, 0, False, None, is_3d)
     
         TestUtil.check_image(imgfile, e)
         
         #testlib.unlink_data_header_files(imgfile)    
         
         (region_2d, region_3d) = self.create_dummy_region(e)

         e2 = EMData()
         e2.read_image(imgfile, 0, False, region_2d, is_3d)
         
         TestUtil.check_image(readfile_2d, e2)
         e2.write_image(readfile_2d, 0, outtype)
         
         if region_3d:
             e4 = EMData()
             e4.read_image(imgfile, 0, False, region_3d, is_3d)
             e4.write_image(readfile_3d, 0, outtype)
             TestUtil.check_image(readfile_3d, e4)

             if outtype == IMAGIC:
                 (hed3d, img3d) = testlib.get_imagic_filename_pair(readfile_3d)
                 os.unlink(hed3d)
                 os.unlink(img3d)
             else:
                 os.unlink(readfile_3d)
             
         if outtype == IMAGIC:
             (hed2d, img2d) = testlib.get_imagic_filename_pair(readfile_2d)
             os.unlink(hed2d)
             os.unlink(img2d)
         else:
             os.unlink(readfile_2d)
             


    def region_write_test(self, imgtype, imgfile, outtype = None):
        if not outtype:
            outtype = imgtype
        
        is_3d = False
        if imgtype == IMAGIC:
            is_3d = True

        imgbase = Util.remove_filename_ext(imgfile)
        ext = Util.get_filename_ext(imgfile)
        writefile_2d = imgbase + "_write_region_2d." + ext
        writefile_3d = imgbase + "_write_region_3d." + ext

        e = EMData()
        e.read_image(imgfile, 0, False, None, is_3d)

        TestUtil.check_image(imgfile, e)
        testlib.unlink_data_header_files(imgfile)
         
        ndims = e.get_ndim()
        e.write_image(writefile_2d, 0, outtype)
        
        if ndims == 3:
            e.write_image(writefile_3d, 0, outtype)

        (region_2d, region_3d) = self.create_dummy_region(e)

        zsize = e.get_zsize()/2
        if zsize == 0:
            zsize = 1
            
        e3 = EMData()
        e3.set_size(e.get_xsize()/2, e.get_ysize()/2, zsize)
        e3.to_zero()

        image_index = 0
        if outtype == SPIDER:
            image_index = e.get_zsize()/2
        
        e3.write_image(writefile_2d, image_index, outtype, False, region_2d)
        TestUtil.check_image(writefile_2d)

        if ndims == 3:
            e3.write_image(writefile_3d, image_index, outtype, False, region_3d)
            TestUtil.check_image(writefile_3d)

            if outtype  == IMAGIC:
                (hed3d, img3d) = testlib.get_imagic_filename_pair(writefile_3d)
                os.unlink(hed3d)
                os.unlink(img3d)
            else:
                os.unlink(writefile_3d)
                
        if outtype == IMAGIC:
            (hed2d, img2d) = testlib.get_imagic_filename_pair(writefile_2d)
            os.unlink(hed2d)
            os.unlink(img2d)
        else:
            os.unlink(writefile_2d)


    def region_read_write_test(self, imgtype, imgfile, outtype = None):
        self.region_read_test(imgtype, imgfile, outtype)
        self.region_write_test(imgtype, imgfile, outtype)
        
        
    def test_mrcio_region(self):
        mrc2d = "test_mrcio_region_2d.mrc"
        mrc3d = "test_mrcio_region_3d.mrc"
        TestUtil.make_image_file(mrc2d, MRC, EM_FLOAT, 32,64)
        TestUtil.make_image_file(mrc3d, MRC, EM_FLOAT, 32,32,12)        
        self.region_read_write_test(MRC, mrc2d)
        self.region_read_write_test(MRC, mrc3d)
        os.unlink(mrc2d)
        os.unlink(mrc3d)

    def test_imagicio_region(self):
        infile = "test_imagicio_region_11.hed"
        TestUtil.make_image_file(infile, IMAGIC, EM_FLOAT, 32,32,64)
        self.region_read_write_test(IMAGIC, infile)

        (hedfile1, imgfile1) = testlib.get_imagic_filename_pair(infile)
        os.unlink(hedfile1)
        os.unlink(imgfile1)


    def test_hdfio_region(self):
        file1 = "test_hdfio_region_1.h5"
        nx = 48
        ny = 64
        nz1 = 1
        TestUtil.make_image_file(file1, HDF, EM_FLOAT, nx, ny, nz1)
        self.region_read_write_test(HDF, file1)

        file2 = "test_hdfio_region_2.h5"
        nz2 = 12
        TestUtil.make_image_file(file2, HDF, EM_FLOAT, nx, ny, nz2)
        self.region_read_write_test(HDF, file2)

        os.unlink(file1)
        os.unlink(file2)
              
"""
    def  test_spiderio_region(self):
        file1 = "test_spiderio_region_1.h5"
        
  
    def test_spiderio_region(self):
        self.region_read_write_test(SINGLE_SPIDER, "spider-single.spi")
        self.region_read_write_test(SPIDER, "spider-stack.spi")

    def test_emio_region(self):
        self.region_read_write_test(EM, "20s2d.em")
        self.region_read_write_test(EM, "stack3d.em")

    def test_pifio_region(self):
        self.region_read_write_test(PIF, "sv-3d.pif")

    def test_xplorio_region(self):
        self.region_read_write_test(XPLOR, "2f.xplor")
"""
        
def test_main():
    TestUtil.set_progname("region")
    #test_support.run_unittest(TestSpiderIO)
    test_support.run_unittest(TestSpiderIO, TestImageIO, TestHdfIO,TestMrcIO, TestImagicIO)

if __name__ == '__main__':
    test_main()


