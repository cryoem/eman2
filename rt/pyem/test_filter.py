#!/bin/env python

from EMAN2 import *
import unittest,os,sys
from test import test_support
import testlib
from pyemtbx.exceptions import *

class TestProcessor(unittest.TestCase):
    """Processor test"""
    
    def test_get_processor_list(self):
        """test get processor list .........................."""
        processor_names = Processors.get_list()
        self.assertEqual(len(processor_names), 107)

        try:
            f2 = Processors.get("_nosuchfilter___")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "NotExistingObjectException")

    def test_BinarizeProcessor(self):
        """test binary processor ............................"""
        imgfile1 = "test_BinarizeFilter.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e = EMData()
        e.read_image(imgfile1)
        fnum = 1000
        f1 = Processors.get("eman1.threshold.binary", {'value': fnum})
        new_params = f1.get_params()
        self.assertEqual(float(new_params["value"]), fnum)
        f1.process(e)
        testlib.check_emdata(e, sys.argv[0])
        os.unlink(imgfile1)

    def test_RangeThreshold(self):
        """test range threshhold processor .................."""
        imgfile1 = "test_RangeThreshold.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e = EMData()
        e.read_image(imgfile1)
        low = 10
        high = 20
        f1 = Processors.get("eman1.threshold.binaryrange", {"low":low, "high":high});
        d1 = f1.get_params()
        self.assertEqual(d1["low"], low)
        self.assertEqual(d1["high"], high)

        xc = 12
        rw = 12.5
        f2 = Processors.get("eman1.mask.ringmean", {"xc":xc, "ring_width":rw})
        d2 = f2.get_params()
        self.assertEqual(d2["xc"], xc)
        self.assertEqual(d2["ring_width"], rw)

        outfile1 = "test_RangeThreshold_out.mrc"
        
        e.process("eman1.threshold.binary", {"value": 200})
        e.write_image(outfile1)

        os.unlink(imgfile1)
        os.unlink(outfile1)
        
class TestCmp(unittest.TestCase):
    """cmp test"""
    
    def test_variance(self):
        """test variance ...................................."""
        imgfile1 = "test_variance_1.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = e1.copy()
        score = e2.cmp("variance", e1, {"keepzero": 0})
        self.assertEqual(score, 0)
        os.unlink(imgfile1)
        
    def test_basic_cmp(self):
        """test basic cmp ..................................."""
        imgfile1 = "test_basic_cmp_1.hed"
        TestUtil.make_image_file(imgfile1, IMAGIC, EM_FLOAT, 16,16,4)

        e1 = EMData()
        e1.read_image(imgfile1, 1)

        e2 = EMData()
        e2.read_image(imgfile1, 2)

        #e1.write_image("test_basic_cmp_out_1.mrc")
        #e2.write_image("test_basic_cmp_out_2.mrc")

        dot_score = e2.cmp("dot", e1, {"evenonly":0})
#        self.assertEqual(dot_score, 19944.0)    #todo: dot score not match, anything wrong?

        variance_score = e2.cmp("variance", e1, {"keepzero":1})
#        self.assertEqual(variance_score, 0)    #todo: score not match, anything wrong?
        
        phase_score = e2.cmp("phase", e1, {})
#        testlib.assertfloat(self, phase_score, 1.6488)    #todo: score not match, anything wrong?
        
        frc_score = e2.cmp("frc", e1, {})
#        testlib.assertfloat(self, frc_score, -0.4011)    #todo: score not match, anything wrong?

        (hed1,img1) = testlib.get_imagic_filename_pair(imgfile1)
        os.unlink(hed1)
        os.unlink(img1)

def test_main():
    test_support.run_unittest(TestProcessor, TestCmp)

if __name__ == '__main__':
    test_main()


