#!/bin/env python

from EMAN2 import *
import unittest,os,sys
from test import test_support
import testlib

class TestFilter(unittest.TestCase):

    def test_get_filter_list(self):
        filternames = Filters.get_list()
        self.assertEqual(len(filternames), 78)


    def test_BinarizeFilter(self):
        imgfile1 = "test_BinarizeFilter.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e = EMData()
        e.read_image(imgfile1)
        fnum = 1000
        f1 = Filters.get("threshold.binary", {'value': fnum})
        new_params = f1.get_params()
        self.assertEqual(float(new_params["value"]), fnum)
        f1.process(e)
        testlib.check_emdata(e, sys.argv[0])
        os.unlink(imgfile1)
        

def test_main():
    test_support.run_unittest(TestFilter)

if __name__ == '__main__':
    test_main()


