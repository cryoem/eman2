#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import os

class TestImageIO(unittest.TestCase):

    def test_imagicio(self):        
        e = EMData()
        imgfile = "start.hed"
        all_imgs = e.read_images(TestUtil.get_debug_image(imgfile))
        outfilebase = "start" + str(os.getpid())
        outfile1 = outfilebase + ".hed"
        outfile2 = outfilebase + ".img"
        
        for img in all_imgs:
            img.append_image(outfile1)
        
        os.unlink(outfile1)
        os.unlink(outfile2)

    def test_emdata_overwriting(self):
        img1 = TestUtil.get_debug_image("groel2d.mrc")
        img2 = TestUtil.get_debug_image("3d.mrc")
        a=EMData()
        b=EMData()
        a.read_image(img1,0)
        b.read_image(img2,0)
        b.read_image(img1,0)

    def test_image_overwriting(self):
        e = EMData()
        e.set_size(20, 20, 1)
        outfile = "test_image_overwriting_" + str(os.getpid()) + ".mrc"
        e.write_image(outfile)
        e.write_image(outfile)
        os.unlink(outfile)
        
        
def test_main():
    test_support.run_unittest(TestImageIO)

if __name__ == '__main__':
    test_main()


