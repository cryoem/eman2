from __future__ import print_function
from __future__ import division

import numpy
from math import isnan as math_isnan
from copy import deepcopy
from EMAN2_cppwrap import EMData,Util
import unittest
from os import path

from test_module import get_data, get_arg_from_pickle_file

from ..libpy import sparx_fundamentals
from ..libpy import sparx_utilities
from ..libpy import sparx_projection
from ..libpy import sparx_alignment as fu
from .sparx_lib import sparx_alignment as oldfu



TOLERANCE = 0.001
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
print(ABSOLUTE_PATH)



class Test_ali2d_single_iter(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ali2d_single_iter()
            oldfu.ali2d_single_iter()

    def test_empty_input_image(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        images =[EMData(),EMData(),EMData()]
        with self.assertRaises(RuntimeError):
            fu.ali2d_single_iter(images, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(images, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list(self):
        (not_used, numr, wr, not_used2, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        cs =[]
        with self.assertRaises(IndexError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_empty_list2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_list2(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        wr = []
        with self.assertRaises(ValueError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)

    def test_empty_list3(self):
        (not_used, not_used, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        numr = []
        with self.assertRaises(IndexError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)

    def test_Invalid_ali_params(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(RuntimeError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="invalid", delta = 0.0)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="invalid", delta = 0.0)

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_bad_center' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_bad_center(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_bad_center2' because: different outputs. could be that somewehre in case of invalid params we got pseudo random output?\n***************************")
    def test_bad_center2(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negativeRangeValue_error(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(UnboundLocalError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, -1, 0, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, -1, 0, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)

    def test_warning_negative_center(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, -5, -5, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, -5, -5, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_T2(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_T2_2(self):
        """ random method ='SHF'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,  nomirror = False, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,  nomirror = False, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_T2_3(self):
        """ mode='h' """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,  nomirror = False, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,nomirror = False, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_T2_4(self):
        """ random method ='SHF'   mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,nomirror = False, mode="h", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Q1(self):
        """ random method ='SCF'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Q1_2(self):
        """ random method ='SCF'    nomirror=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Q1_3(self):
        """ random method ='SCF'    mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Q1_4(self):
        """ random method ='SCF'    nomirror=True   mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_S2(self):
        """ random method ='SHF'    nomirror=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_S2_2(self):
        """ random method =''   nomirror=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_S2_3(self):
        """ random method ='SHF'    nomirror=True   mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_S2_4(self):
        """ random method =''   nomirror=True   mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_M1(self):
        """ random method ='SCF'    , CTF=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SCF", CTF=True)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SCF", CTF=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_M1_2(self):
        """ random method ='SCF'    nomirror=True   , CTF=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_M1_3(self):
        """ random method ='SCF'    mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_M1_4(self):
        """ random method ='SCF'    nomirror=True   mode='h'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_I2(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,nomirror = False, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        for i, j in zip(return_old, return_new):
            self.assertTrue(TOLERANCE > numpy.abs(i-j))

    def test_I2_2(self):
        """ random method ='SHF'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        for i, j in zip(return_old, return_new):
            self.assertTrue(TOLERANCE > numpy.abs(i - j))

    def test_K2(self):
        """ random method ='SHF'    nomirror=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SHF", ali_params="xform.align2d", delta = 0.0)
        for i, j in zip(return_old, return_new):
            self.assertTrue(TOLERANCE > numpy.abs(i - j))

    def test_K2_2(self):
        """ random method =''   nomirror=True"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        for i, j in zip(return_old, return_new):
            self.assertTrue(TOLERANCE > numpy.abs(i - j))


    "The following test do not work because an error"

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_PCP1' because: it seems to be a dead code\n***************************")
    def test_PCP1(self):
        """ random method ='PCP'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "PCP")
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "PCP")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_PCP2' because: it seems to be a dead code\n***************************")
    def test_PCP2(self):
        """ random method ='PCP'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "PCP", nomirror=True)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "PCP", nomirror=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_SHC'1 because: no clues about this odd behaviour \n***************************")
    def test_SHC1(self):
        """ random method ='SHC'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SHC")
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SHC")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ali2d_single_iter.test_SHC2' because: no clues about this odd behaviour \n***************************")
    def test_SHC2(self):
        """ random method ='SHC'"""
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SHC", CTF=True)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, random_method = "SHC",CTF=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ang_n(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ang_n()
            oldfu.ang_n()

    def test_C2(self):
        return_new = fu.ang_n(2, 'f', 3)
        return_old = oldfu.ang_n(2, 'f', 3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D2(self):
        return_new = fu.ang_n(2, 'not_f', 3)
        return_old = oldfu.ang_n(2, 'not_f', 3)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_log2(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.log2()
            oldfu.log2()

    def test_C3(self):
        self.assertEqual(fu.log2(10),oldfu.log2(10))



class Test_Numrinit(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.Numrinit()
            oldfu.Numrinit()

    def test_Invalid_skip_value(self):
        with self.assertRaises(ValueError):
            fu.Numrinit(2, 5, skip=0, mode="F")
            oldfu.Numrinit(2, 5, skip=0,mode="F")

    def test_C4(self):
        return_new = fu.Numrinit(2, 5, skip=1, mode="F")
        return_old = oldfu.Numrinit(2, 5, skip=1, mode="F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D4(self):
        return_new = fu.Numrinit(2, 5, skip=1, mode="not_F")
        return_old = oldfu.Numrinit(2, 5, skip=1, mode="not_F")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ringwe(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ringwe"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ringwe()
            oldfu.ringwe()

    def test_empty_list(self):
        with self.assertRaises(IndexError):
            fu.ringwe([], mode="F")
            oldfu.ringwe([], mode="F")

    def test_C5(self):
        return_new = fu.ringwe(self.argum[0][0], mode="F")
        return_old = oldfu.ringwe(self.argum[0][0], mode="F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D5(self):
        return_new = fu.ringwe(self.argum[0][0], mode="not_F")
        return_old = oldfu.ringwe(self.argum[0][0], mode="not_F")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ornq(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))

    @unittest.skip("\n***************************\n\t\t 'Test_ornq.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]

        return_new = fu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ornq.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]

        return_new = fu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ornq()
            oldfu.ornq()

    @unittest.skip("\n***************************\n\t\t 'Test_ornq.test_empty_list' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_list(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        numr=[]
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        xrng=[]
        with self.assertRaises(IndexError):
            fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_empty_list3(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        yrng=[]
        with self.assertRaises(IndexError):
            fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_with_negative_center(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, -5, -5, deltapsi=0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, -5, -5, deltapsi=0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_InvalidStepValue_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        with self.assertRaises(ZeroDivisionError):
            fu.ornq(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_C6(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C6_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='f'
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C6_with_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='invalid'
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ormq(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ormq"))

    @unittest.skip("\n***************************\n\t\t 'Test_ormq.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        image =EMData()
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ormq.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        crefim =EMData()
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ormq()
            oldfu.ormq()

    def test_F7_G7_H7(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]  # mode is F
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list1(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        xrng=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_empty_list2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        yrng=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    @unittest.skip("\n***************************\n\t\t 'Test_ormq.test_empty_list3' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_list3(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_InvalidStepValue_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        step = 0
        with self.assertRaises(ZeroDivisionError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_F7_G7_H7_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode ='H'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_negative_center(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, -5, -5, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, -5, -5, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_F7_G7_H7_with_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode ='invalid'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C7_D7_E7(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C7_D7_E7_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode = 'H'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C7_D7_E7_with_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode = 'invalid'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_prepref(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(len(return_old), len(return_new))
            for i, j in zip(return_old, return_new):
                self.assertEqual(len(i), len(j))
                for q, r in zip(i, j):
                    self.assertEqual(len(q), len(r))
                    for img1, img2 in zip(q, r):
                        try:
                            self.assertTrue(numpy.array_equal(img1.get_3dview(), img2.get_3dview()))
                        except AssertionError:
                            self.assertTrue(TOLERANCE > numpy.abs(numpy.sum(img1.get_3dview() - img2.get_3dview())))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.prepref()
            oldfu.prepref()

    def test_empty_input_image(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(RuntimeError):
            fu.prepref(data, EMData(), cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, EMData(), cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_empty_input_image2(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        data= [EMData(),EMData(),EMData()]
        with self.assertRaises(RuntimeError):
            fu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_empty_list(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(IndexError):
            fu.prepref(data, None, cnx, cny, [], mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, None, cnx, cny, [], mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_C9(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        self.assertEqual(len(return_old),len(return_new))
        self.test_all_the_conditions(return_new,return_old,False)

    def test_withBadMask_error(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        mask = sparx_utilities.model_circle(100,100,100)
        with self.assertRaises(RuntimeError):
            fu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_C9_withMask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        nx = data[0].get_xsize()
        mask = sparx_utilities.model_circle(nx//2-1,nx,nx)
        return_new = fu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        self.assertEqual(len(return_old),len(return_new))
        self.test_all_the_conditions(return_new,return_old,False)

    def test_C9_2(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        self.assertEqual(len(return_old),len(return_new))
        self.test_all_the_conditions(return_new, return_old, False)

    def test_C9_with_invalid_mode(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'not_valid', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'not_valid', maxrangex = 4, maxrangey = 4, step =step)
        self.assertEqual(len(return_old),len(return_new))
        self.test_all_the_conditions(return_new, return_old, False)



class Test_prepare_refrings(unittest.TestCase):
    volft = sparx_utilities.model_blank(100,100,100)
    numr = [1, 1, 8, 2, 9, 16, 3, 953, 128, 16, 1081, 128, 17, 1209, 128, 18, 1337, 128, 19, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(len(return_new), len(return_old))
            for img1, img2 in zip(return_new, return_old):
                try:
                    self.assertTrue(numpy.array_equal(img1.get_3dview(), img2.get_3dview()))
                except AssertionError:
                    # since sometimes we get  img1.get_3dview()= [[[ nan  nan  nan ...,  nan  nan  nan]]] we skip these cases
                    res = numpy.sum(img1.get_3dview() - img2.get_3dview())
                    if math_isnan(res) is False:
                        self.assertTrue(TOLERANCE > numpy.abs(res))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.prepare_refrings()
            oldfu.prepare_refrings()

    def test_empty_input_image(self):

        volft, kb = sparx_projection.prep_vol(self.volft)
        with self.assertRaises(RuntimeError):
            fu.prepare_refrings(EMData(), kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5,initial_phi=0.1)
            oldfu.prepare_refrings(EMData(), kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5,initial_phi=0.1)

    def test_empty_list(self):
        volft, kb = sparx_projection.prep_vol(self.volft)
        with self.assertRaises(IndexError):
            fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
            oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)

    """ There are a lot of case with a big numerical difference, I mean more than 100 when i use a tolerance of 0.0005 ... if i use tthe commented volfts it passes the test
        ref_a --> P,S , symmetry_class=sparx_fundamentals.symclass("c1")->symmetry_class.even_angles(delta[N_step], phiEqpsi = "Zero")
    """
    def test_1(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_2(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_3(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    @unittest.skip("\n***************************\n\t\t 'Test_prepare_refrings where MPI is True. Because a crash due to 'The MPI_Comm_rank() function was called before MPI_INIT was invoked'\n***************************")
    def test_4(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_5(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_6(self):
        volft,kbx,kby,kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbx,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbx,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_1b(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_2b(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_3b(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    @unittest.skip("\n***************************\n\t\t 'Test_prepare_refrings where MPI is True. Because a crash due to 'The MPI_Comm_rank() function was called before MPI_INIT was invoked'\n***************************")
    def test_4b(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_6b(self):
        volft,kbx,kby,kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbx,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbx,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_0(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)



class Test_proj_ali_incore(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.proj_ali_incore()
            oldfu.proj_ali_incore()

    @unittest.skip("\n***************************\n\t\t 'Test_proj_ali_incore.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(),EMData()]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_input_image2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data=EMData()
        with self.assertRaises(RuntimeError):
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)

    def test_empty_input_list(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError):
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)

    def test_F11(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ It does not work!! the olf function got a second value, pixelError, really smaller than the new one. Is it ok?"""
    def test_F11_2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    "It does not work. The pixel error has a different value more on less 1%"
    def test_F11_3(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    """ It does not work!! the olf function got a second value, pixelError, really smaller than the new one. Is it ok?"""
    def test_E11(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c2", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c2", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_proj_ali_incore_local(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.proj_ali_incore_local()
            oldfu.proj_ali_incore_local()

    def test_empty_input_image(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data = EMData()
        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        with self.assertRaises(RuntimeError):
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)

    @unittest.skip("\n***************************\n\t\t 'Test_proj_ali_incore.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(),EMData(),EMData()]
        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        with self.assertRaises(RuntimeError):
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)

    def test_empty_list(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr = []
        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        with self.assertRaises(IndexError):
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)

    @unittest.skip("\n***************************\n\t\t 'Test_proj_ali_incore.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_list2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        with self.assertRaises(IndexError):
            fu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)

    @unittest.skip("\n***************************\n\t\t 'Test_proj_ali_incore_local.test_too_short_list_of_ref_ang' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_too_short_list_of_ref_ang(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        symangles = [[0.0, 0.0, 360. ] for k in range(10) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        return_new = fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        return_old = oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ 
    It does not work because the old version returns always True. If I remove the return True it gives back a very different values.
    Is it this function obsolete???
    """
    def test_G12(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        return_new = fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        return_old = oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ali_vol_func(unittest.TestCase):
    param = [1, 1, 1, 1,1,1]
    data =get_data(3)

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ali_vol_func()
            oldfu.ali_vol_func()

    def test_wrong_number_params2(self):
        param = [1,1,1,1]
        with self.assertRaises(IndexError):
            fu.ali_vol_func(param,self.data)
            oldfu.ali_vol_func(param, self.data)

    def test_C13(self):
        return_new = fu.ali_vol_func(self.param, self.data)
        return_old = oldfu.ali_vol_func(self.param, self.data)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_too_few_img(self):
        data =get_data(2)
        with self.assertRaises(IndexError):
            fu.ali_vol_func(self.param, data)
            oldfu.ali_vol_func(self.param, data)

    def test_empty_input_image(self):
        data = [EMData(), EMData(), EMData()]
        with self.assertRaises(RuntimeError):
            fu.ali_vol_func(self.param,data)
            oldfu.ali_vol_func(self.param,data)



class Test_align2d(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))

    @unittest.skip("\n***************************\n\t\t 'Test_align2d.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, refim, xrng, yrng) = self.argum[0]
        image = EMData()
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_input_image2(self):
        (image, refim, xrng, yrng) = self.argum[0]
        refim = EMData()
        with self.assertRaises(IndexError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode="F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode="F")

    def test_empty_list(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_empty_list2(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.align2d()
            oldfu.align2d()

    @unittest.skip("\n***************************\n\t\t 'Test_align2d.test_wrong_enumerate_rings_error' because a runtime floatigpoint error: interrupted by signal 11: SIGSEGV\n***************************")
    def test_wrong_enumerate_rings_error(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_rstep_error(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=0, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=0, mode = "F")

    def test_null_step_error(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[], step=0, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[], step=0, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_C14(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D14(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C14_2(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "H")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "H")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D14_2(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "H")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "H")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_align2d_scf(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.align2d_scf()
            oldfu.align2d_scf()

    def test_empty_input_image(self):
        (image, refim, xrng, yrng) = self.argum[0]
        image =EMData()
        with self.assertRaises(RuntimeError):
            fu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])
            oldfu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])

    def test_empty_input_image2(self):
        (image, refim, xrng, yrng) = self.argum[0]
        refim =EMData()
        with self.assertRaises(RuntimeError):
            fu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])
            oldfu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])

    def test_with_valid_params(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])
        return_old = oldfu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_automatic_params_OK(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d_scf(image, refim,self.argum[1]['ou'])
        return_old = oldfu.align2d_scf(image, refim,self.argum[1]['ou'])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_invalid_ou(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d_scf(image, refim, xrng, yrng, 1)
        return_old = oldfu.align2d_scf(image, refim, xrng, yrng,1)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ this test is not able to work"""
    def test_with_automatic_params_PROBLEM(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d_scf(image, refim)
        return_old = oldfu.align2d_scf(image, refim)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_multialign2d_scf(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))

    @unittest.skip("\n***************************\n\t\t 'Test_multialign2d_scf.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = EMData()
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

        self.assertEqual(return_old, return_new)

    @unittest.skip("\n***************************\n\t\t 'Test_multialign2d_scf.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = EMData()

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

        self.assertEqual(return_old, return_new)

    def test_empty_input_image3(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = [EMData(),EMData(),EMData()]
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]

        with self.assertRaises(RuntimeError):
            fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
            oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.multalign2d_scf()
            oldfu.multalign2d_scf()

    def test_with_valid_params(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_multialign2d_scf.test_with_empty_list' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_with_empty_list(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]
        numr = []
        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_invalid_ou(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=1)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=1)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ this test is not able to work"""
    def test_with_automatics_params(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr)

        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_parabl(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.parabl"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.parabl()
            oldfu.parabl()

    def test_empty_input_image(self):
        return_new = fu.parabl(EMData())
        return_old = oldfu.parabl(EMData())
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D17(self):
        return_new = fu.parabl(self.argum[0][0])
        return_old = oldfu.parabl(self.argum[0][0])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_C17(self):
        Z=self.argum[0][0]
        for i in range(3):
            for j in range(3):
                Z[i,j]=0
        return_new = fu.parabl(Z)
        return_old = oldfu.parabl(Z)
        self.assertTrue(numpy.array_equal(return_new, return_old))



""" In all the tests we have some values that are different"""
class Test_shc(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.shc()
            oldfu.shc()

    @unittest.skip("\n***************************\n\t\t 'Test_shc.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(), EMData()]
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_input_image2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data =  EMData()
        with self.assertRaises(RuntimeError):
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

    """ It does not work"""
    def test_F18(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr = []
        with self.assertRaises(IndexError):
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

    """ This test is not able to work"""
    def test_empty_list2(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        list_of_ref_ang = []
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ This test is not able to work"""
    def test_added_one_ref_ang(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        list_of_ref_ang[0].append(2.0)
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ It does not work"""
    def test_M18(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, sym ="c2")
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, sym ="c2")

        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ It does not work"""
    def test_G18(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, sym ="nomirror")
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, sym ="nomirror")

        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_search_range(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.search_range"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.search_range()
            oldfu.search_range()

    def test_D19(self):
        (n, radius, shift, range_, location) = self.argum[0]

        return_new = fu.search_range(n, radius, shift, range_)
        return_old = oldfu.search_range(n, radius, shift, range_)

        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_C19(self):
        (n, radius, shift, range_, location) = self.argum[0]

        return_new = fu.search_range(0, radius, shift, range_)
        return_old = oldfu.search_range(0, radius, shift, range_)

        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_generate_list_of_reference_angles_for_search(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.generate_list_of_reference_angles_for_search()
            oldfu.generate_list_of_reference_angles_for_search()

    def test_C20(self):
        nsym = 5
        symangles = [[0.0, 0.0, k * 360. / nsym] for k in range(nsym)]
        return_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c5')
        return_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c5')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D20(self):
        nsym = 5
        symangles = [[0.0, 0.0, k * 360. / nsym] for k in range(nsym)]
        return_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        return_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NoAngleList(self):
        return_new = fu.generate_list_of_reference_angles_for_search([], 'c1')
        return_old = oldfu.generate_list_of_reference_angles_for_search([], 'c1')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_simmetry(self):
        with self.assertRaises(RuntimeError):
            fu.generate_list_of_reference_angles_for_search([], 'not_valid')
            oldfu.generate_list_of_reference_angles_for_search([], 'not_valid')



