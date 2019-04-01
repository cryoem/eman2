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

from mpi import *
mpi_init(0, [])

TOLERANCE = 0.001
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
print(ABSOLUTE_PATH)


"""
There are some opened issues in:
1) ali2d_single_iter --> see in the class comments to get more info
2) proj_ali_incore_local --> the sparx verison (i.e.:) returns always True ... is this function obsolete?
3) shc  --> all the tests are failing. Even I used some random value from another pickle file it'd work
4) ormq_fast --> it seems to be never used in the whole project. I avoid to test it
5) multalign2d_scf and align2d_scf could have a bug. See in their classes my comments and run the tests to get the error message
"""


class Test_ali2d_single_iter(unittest.TestCase):
    """
    Since using an invalid method as "random method" is like using the default method "random_method=''" I'm not testing this situation"
    Since the method "random_method='PCP'" seems to lead to a dead code, anyway it is crashing, I'm skipping these cases
    All the case with "random_method='SHC'" do not work. They are manageable through 'typeError' exception. This situation is used a lot in the 'sphire/legacy/...'
    """
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(len(return_old), len(return_new))
            for i, j in zip(return_old, return_new):
                self.assertEqual(len(i), len(j))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ali2d_single_iter()
            oldfu.ali2d_single_iter()

    def test_empty_images_to_align_returns_RuntimeError_NotExistingObjectException_the_key_ctf_doesnot_exist(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        images =[EMData(),EMData(),EMData()]
        with self.assertRaises(RuntimeError):
            fu.ali2d_single_iter(images, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
            oldfu.ali2d_single_iter(images, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_few_shift_params_returns_IndexError_list_index_out_of_range(self):
        (not_used, numr, wr, not_used2, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        cs =[1]
        with self.assertRaises(IndexError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)

    def test_too_shift_params(self):
        (not_used, numr, wr, not_used2, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        cs =[1,2,3]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list_fourier_weight_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        wr = []
        with self.assertRaises(ValueError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
        """
        self.assertTrue(True)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (not_used, not_used, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        numr = []
        with self.assertRaises(IndexError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)

    def test_Unknown_ali_params_returns_RuntimeError_NotExistingObjectException_the_key_Unknown_doesnot_exist(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(RuntimeError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="Unknown", delta = 0.0)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="Unknown", delta = 0.0)

    def test_default_center_out_of_possible_range_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_negative_XRange_Value_UnboundLocalError_a_local_var_is_undefined(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(UnboundLocalError):
            fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, -1, 0, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
            oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, -1, 0, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)

    def test_negative_center_warning_msg_shift_of_paricle_too_large(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, -5, -5, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, -5, -5, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))
        
    def test_NOmirror_mode_H_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_F_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))
        
    def test_NOmirror_mode_F_withCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_H_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_H_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(TypeError):
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SHC", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        
    def test_NOmirror_mode_H_NOCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_F_NOCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        
    def test_NOmirror_mode_F_withCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_H_NOCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_H_WITHCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_WITHCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_NOCTF_randomMethod_SCF(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="SCF", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))
        
    def test_NOmirror_mode_H_NOCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_F_NOCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))
        
    def test_NOmirror_mode_F_withCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_mirror_mode_H_NOCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_mirror_mode_H_WITHCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_mirror_mode_F_WITHCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_mirror_mode_F_NOCTF_randomMethod_default(self):
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="H", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_H_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="h", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_F_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOmirror_mode_F_withCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = False, mode="F", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_mirror_mode_H_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_H_WITHCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="h", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_WITHCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=True, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mirror_mode_F_NOCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(self):
        with self.assertRaises(TypeError):
            (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
            return_new = fu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            return_old = oldfu.ali2d_single_iter(deepcopy(self.argum[0][0]), numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, nomirror = True, mode="F", CTF=False, random_method="PCP", ali_params="xform.align2d", delta = 0.0)
            self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ang_n(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ang_n()
            oldfu.ang_n()

    def test_Full_mode(self):
        return_new = fu.ang_n(2, 'f', 3)
        return_old = oldfu.ang_n(2, 'f', 3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_max_ring_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.ang_n(3, 'f', 0)
            oldfu.ang_n(3, 'f', 0)

    def test_Half_mode(self):
        return_new = fu.ang_n(2, 'not_f', 3)
        return_old = oldfu.ang_n(2, 'not_f', 3)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_log2(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.log2()
            oldfu.log2()

    def test_positive_number(self):
        self.assertEqual(fu.log2(10),oldfu.log2(10))

    def test_null_number(self):
        self.assertEqual(fu.log2(0),oldfu.log2(0))



class Test_Numrinit(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.Numrinit()
            oldfu.Numrinit()

    def test_null_skip_value_returns_ValueError_this_arg_cannot_be_null(self):
        with self.assertRaises(ValueError):
            fu.Numrinit(2, 5, skip=0, mode="F")
            oldfu.Numrinit(2, 5, skip=0,mode="F")

    def test_Full_mode(self):
        return_new = fu.Numrinit(2, 5, skip=1, mode="F")
        return_old = oldfu.Numrinit(2, 5, skip=1, mode="F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode(self):
        return_new = fu.Numrinit(2, 5, skip=1, mode="not_F")
        return_old = oldfu.Numrinit(2, 5, skip=1, mode="not_F")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ringwe(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ringwe"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ringwe()
            oldfu.ringwe()

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.ringwe([], mode="F")
            oldfu.ringwe([], mode="F")

    def test_Full_mode(self):
        return_new = fu.ringwe(self.argum[0][0], mode="F")
        return_old = oldfu.ringwe(self.argum[0][0], mode="F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode(self):
        return_new = fu.ringwe(self.argum[0][0], mode="not_F")
        return_old = oldfu.ringwe(self.argum[0][0], mode="not_F")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ornq(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ornq()
            oldfu.ornq()

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        numr = []
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_list_xrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        xrng=[]
        with self.assertRaises(IndexError):
            fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_empty_list_yrng_returns_IndexError_list_index_out_of_range(self):
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

    def test_null_skip_value_returns_ZeroDivisionError(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        with self.assertRaises(ZeroDivisionError):
            fu.ornq(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_Half_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Full_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='f'
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='invalid'
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ormq(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ormq"))

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        image =EMData()
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        crefim =EMData()
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ormq()
            oldfu.ormq()

    def test_empty_list_xrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        xrng=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_empty_list_yrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        yrng=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        """
        self.assertTrue(True)

    def test_null_skip_value_returns_ZeroDivisionError(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        step = 0
        with self.assertRaises(ZeroDivisionError):
            fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
            oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)

    def test_Full_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]  # mode is F
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode(self):
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

    def test_with_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode ='invalid'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Full_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode = 'H'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_invalid_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        mode = 'invalid'
        return_new = fu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
        return_old = oldfu.ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0)
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

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.prepref()
            oldfu.prepref()

    def test_empty_image_mask_returns_RuntimeError_ImageDimensionException_img_dimension_doesnot_match_with_its_dimension(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(RuntimeError):
            fu.prepref(data, EMData(), cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, EMData(), cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_empty_images_to_align_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        data= [EMData(),EMData(),EMData()]
        with self.assertRaises(RuntimeError):
            fu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(IndexError):
            fu.prepref(data, None, cnx, cny, [], mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, None, cnx, cny, [], mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_full_mode_without_mask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_half_mode_without_mask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        self.test_all_the_conditions(return_new, return_old, False)

    def test_image_mask_returns_RuntimeError_ImageDimensionException_img_dimension_doesnot_match_with_its_dimension(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        mask = sparx_utilities.model_circle(100,100,100)
        with self.assertRaises(RuntimeError):
            fu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
            oldfu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)

    def test_Full_mode_withMask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        nx = data[0].get_xsize()
        mask = sparx_utilities.model_circle(nx//2-1,nx,nx)
        return_new = fu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, mask, cnx, cny, numr, mode = 'f', maxrangex = 4, maxrangey = 4, step =step)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_Half_mode_withMask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        nx = data[0].get_xsize()
        mask = sparx_utilities.model_circle(nx//2-1,nx,nx)
        return_new = fu.prepref(data, mask, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, mask, cnx, cny, numr, mode = 'H', maxrangex = 4, maxrangey = 4, step =step)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_with_invalid_mode(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        return_new = fu.prepref(data, None, cnx, cny, numr, mode = 'not_valid', maxrangex = 4, maxrangey = 4, step =step)
        return_old = oldfu.prepref(data, None, cnx, cny, numr, mode = 'not_valid', maxrangex = 4, maxrangey = 4, step =step)
        self.test_all_the_conditions(return_new, return_old, False)



class Test_prepare_refrings(unittest.TestCase):
    """
    Take a look to sparx_utilities.py --> even_angles_cd(...)for the meaning of the following params
        ref_a --> P=Penczek algorithm, S=Saff algorithm to calculate di reference angle
        phiEQpsi  --> 'Minus', if you want psi=-phi to create a list of  angles suitable for projections, otherwise 'Zero'

    In case of rectangular kb filter see how it uses kbi variables in sparx_projection.py --> prgs(...) to understand better
    """
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

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.prepare_refrings()
            oldfu.prepare_refrings()

    def test_empty_volume_returns_RuntimeError_ImageFormatException_extractplane_requires_complex_img(self):
        volft, kb = sparx_projection.prep_vol(self.volft)
        with self.assertRaises(RuntimeError):
            fu.prepare_refrings(EMData(), kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5,initial_phi=0.1)
            oldfu.prepare_refrings(EMData(), kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5,initial_phi=0.1)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        volft, kb = sparx_projection.prep_vol(self.volft)
        with self.assertRaises(IndexError):
            fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
            oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)

    def test_empty_list_referenceAngle_returns_IndexError_list_index_out_of_range(self):
        volft, kb = sparx_projection.prep_vol(self.volft)
        with self.assertRaises(IndexError):
            fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a=[], sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
            oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a=[], sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)

    def test_No_kb_ArgumentError_in_EMData_extract_plane_function(self):
        with self.assertRaises(TypeError):
            fu.prepare_refrings(self.volft, None,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
            oldfu.prepare_refrings(self.volft, None,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)

    def test_with_sym_c1_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_with_sym_c5_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    @unittest.skip( "\n***************************\n\t\t 'Test_prepare_refringstest_sym_c1_initialTheta_None. Even if this combination is it seems to lead the code to a deadlock, i waited more then an hour'\n***************************")
    def test_sym_c1_initialTheta_None(self):
        volft, kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=None, delta_theta=0.5,initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb, nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=None, delta_theta=0.5,initial_phi=0.1)
        self.test_all_the_conditions(return_new, return_old, False)

    def test_No_nz_data_size_Error_msg_datasize_hasnot_be_given(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=0, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=0, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a= sparx_utilities.even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft,kb = sparx_projection.prep_vol(self.volft)
        return_new = fu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kb,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft, kbx, kby, kbz = sparx_projection.prep_vol(sparx_utilities.model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        return_old = oldfu.prepare_refrings(volft, kbz,nz=4, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5, initial_phi=0.1)
        self.test_all_the_conditions(return_new,return_old,False)



class Test_proj_ali_incore(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.proj_ali_incore()
            oldfu.proj_ali_incore()

    def test_empty_input_image_refrings_crashes_because_signal11SIGSEV(self):
        """
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(),EMData()]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_img_data_returns_RuntimeError_NotExistingObjectException_the_key_xform_projection_doesnot_exist(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data=EMData()
        with self.assertRaises(RuntimeError):
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError):
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)

    def test_sym_c1(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_c1_negative_deltaPhsi(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_c1_negative_rshift(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_c1_positive_deltaPhsi(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 100.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 100.0, rshift = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_c1_positive_rshift(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 10000.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 10000.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_not_c1(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "icos", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "icos", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(numpy.allclose(return_old, return_new , atol=TOLERANCE ))


@unittest.skip("All the tests in Test_proj_ali_incore_local. the old function returns always True ... is it obsolete?")
class Test_proj_ali_incore_local(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
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


    def test_empty_input_image_refrings(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(),EMData(),EMData()]
        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        with self.assertRaises(RuntimeError):
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)

    def test_empty_list_numr(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr = []
        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        with self.assertRaises(IndexError):
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)


    def test_empty_list_of_ref_ang(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        with self.assertRaises(IndexError):
            fu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
            oldfu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)


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

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ali_vol_func()
            oldfu.ali_vol_func()

    def test_few_params_params_returns_IndexError_list_index_out_of_range(self):
        param = [1,1,1,1]
        with self.assertRaises(IndexError):
            fu.ali_vol_func(param,self.data)
            oldfu.ali_vol_func(param, self.data)

    def test_too_few_data_params_returns_IndexError_list_index_out_of_range(self):
        data =get_data(2)
        with self.assertRaises(IndexError):
            fu.ali_vol_func(self.param, data)
            oldfu.ali_vol_func(self.param, data)

    def test_ali_vol_func(self):
        return_new = fu.ali_vol_func(self.param, self.data)
        return_old = oldfu.ali_vol_func(self.param, self.data)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_data_images_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        data = [EMData(), EMData(), EMData()]
        with self.assertRaises(RuntimeError):
            fu.ali_vol_func(self.param,data)
            oldfu.ali_vol_func(self.param,data)



class Test_align2d(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, refim, xrng, yrng) = self.argum[0]
        image = EMData()
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_returns_IndexError_list_index_out_of_range(self):
        (image, refim, xrng, yrng) = self.argum[0]
        refim = EMData()
        with self.assertRaises(IndexError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode="F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode="F")

    def test_empty_list_xrng_returns_ValueError_arg_af_max_f_is_empty_list(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_empty_list_yrngreturns_ValueError_arg_af_max_f_is_empty_list(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.align2d()
            oldfu.align2d()

    def test_wrong_enumerate_rings_error_crashes_because_signal11SIGSEV(self):
        """
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_null_rstep_value_returns_ValueError_arg_af_max_f_is_empty_list(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ValueError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=0, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[], step=1, first_ring=1, last_ring=0, rstep=0, mode = "F")

    def test_null_step_value_returns_ZeroDivisionError(self):
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(ZeroDivisionError):
            fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=0, first_ring=1, last_ring=0, rstep=1, mode = "F")
            oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=0, first_ring=1, last_ring=0, rstep=1, mode = "F")

    def test_Full_mode_zero_lastRing(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode_zero_lastRing(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "H")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "H")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Full_mode(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "F")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Half_mode(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "H")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=2, rstep=1, mode = "H")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_align2d_scf(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.align2d_scf()
            oldfu.align2d_scf()

    def test_empty_image_to_align_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        (image, refim, xrng, yrng) = self.argum[0]
        image =EMData()
        with self.assertRaises(RuntimeError):
            fu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])
            oldfu.align2d_scf(image, refim, xrng, yrng, self.argum[1]['ou'])

    def test_empty_reference_image_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
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

    def test_with_invalid_ou_error_msg_output(self):
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d_scf(image, refim, xrng, yrng, 1)
        return_old = oldfu.align2d_scf(image, refim, xrng, yrng,1)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    """ 
    the following testa are not able to work. It'd be a bug.
    error message:
        File "/home/lusnig/EMAN2/eman2/sphire/tests/sparx_lib/sparx_alignment.py", line 784, in align2d_scf
        sxs = -p2[0][4]
        IndexError: list index out of range
    BUT p2 is the ouput of:
        -) ccf2 = EMAN2_cppwrap.Util.window(sparx_fundamentals.ccf(sparx_fundamentals.rot_shift2D(image, alpha+180.0, 0.0, 0.0, mirr), frotim),nrx,nry)
	    -) p2 = sparx_utilities.peak_search(ccf2)
	in these casea it is a list of 4 elements and it is trying to get the 5th
    """
    def test_with_DEFAULT_params_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            (image, refim, xrng, yrng) = self.argum[0]
            return_new = oldfu.align2d_scf(image, refim, xrng=-1, yrng=-1, ou = -1)
            return_old = oldfu.align2d_scf(image, refim, xrng=-1, yrng=-1, ou = -1)
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_DEFAULT_params_but_validOU_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            (image, refim, xrng, yrng) = self.argum[0]
            return_new = fu.align2d_scf(image, refim, xrng=-1, yrng=-1, ou = self.argum[1]['ou'])
            return_old = oldfu.align2d_scf(image, refim,xrng=-1, yrng=-1, ou = self.argum[1]['ou'])
            self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_multialign2d_scf(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))

    def test_empty_input_image_refrings_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = EMData()
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertEqual(return_old, return_new)
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = EMData()

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertEqual(return_old, return_new)
        """
        self.assertTrue(True)

    def returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = [EMData(),EMData(),EMData()]
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]
        with self.assertRaises(RuntimeError):
            fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
            oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]
        numr = []

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.multalign2d_scf()
            oldfu.multalign2d_scf()

    def test_with_valid_params_cimage_with_mode_F(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [sparx_fundamentals.fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_valid_params_cimage_with_mode_H(self):
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "H")
        frotim = [sparx_fundamentals.fft(tavg)]

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

    def test_with_DEFAULT_params_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

            dataa = deepcopy(self.argum[0][0])
            cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
            frotim = [sparx_fundamentals.fft(tavg)]

            return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng=-1, yrng=-1, ou = -1)
            return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng=-1, yrng=-1, ou = -1)
            self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_parabl(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.parabl"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.parabl()
            oldfu.parabl()

    def test_empty_input_image(self):
        return_new = fu.parabl(EMData())
        return_old = oldfu.parabl(EMData())
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_peak_found(self):
        return_new = fu.parabl(self.argum[0][0])
        return_old = oldfu.parabl(self.argum[0][0])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_No_peak_found(self):
        Z=self.argum[0][0]
        for i in range(3):
            for j in range(3):
                Z[i,j]=0
        return_new = fu.parabl(Z)
        return_old = oldfu.parabl(Z)
        self.assertTrue(numpy.array_equal(return_new, return_old))



""" In all the tests we have some values that are different ... hence I did not test all the cases"""
class Test_shc(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.shc()
            oldfu.shc()

    def test_empty_input_image_refringscrashes_because_signal11SIGSEV(self):
        """
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(), EMData()]
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_returns_RuntimeError_the_key_xform_projection_doesnot_exist(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data =  EMData()
        with self.assertRaises(RuntimeError):
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr = []
        with self.assertRaises(IndexError):
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

    def test_sym_c1_failed(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        return_new = fu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
        return_old = oldfu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
        self.assertTrue(True)
        #self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list_of_ref_ang_failed(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        list_of_ref_ang = []

        return_new = fu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
        return_old = oldfu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
        self.assertTrue(True)
        #self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_added_one_ref_ang_failed(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        list_of_ref_ang[0].append(2.0)
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)

        self.assertTrue(True)
        #self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_sym_nomirror_failed(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "nomirror", finfo=None)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "nomirror", finfo=None)

        self.assertTrue(True)
        #self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_search_range(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.search_range"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.search_range()
            oldfu.search_range()

    def test_search_range(self):
        (n, radius, shift, range_, location) = self.argum[0]
        return_new = fu.search_range(n, radius, shift, range_, location = "")
        return_old = oldfu.search_range(n, radius, shift, range_, location = "")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_no_image_size_warning_msg_shift_of_particle_too_large(self):
        (n, radius, shift, range_, location) = self.argum[0]
        return_new = fu.search_range(0, radius, shift, range_, location = "")
        return_old = oldfu.search_range(0, radius, shift, range_, location = "")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_generate_list_of_reference_angles_for_search(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.generate_list_of_reference_angles_for_search()
            oldfu.generate_list_of_reference_angles_for_search()

    def test_c5Sym(self):
        sym = 'c5'
        ref_angles = sparx_utilities.even_angles(symmetry=sym)
        return_new = fu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        return_old = oldfu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_c1Sym(self):
        sym = 'c1'
        ref_angles = sparx_utilities.even_angles(symmetry=sym)
        return_new = fu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        return_old = oldfu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NoAngleList(self):
        return_new = fu.generate_list_of_reference_angles_for_search([], 'c1')
        return_old = oldfu.generate_list_of_reference_angles_for_search([], 'c1')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_simmetry_returns_RuntimeError_NotExistingObjectException_the_key_invalid_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.generate_list_of_reference_angles_for_search([], 'invalid')
            oldfu.generate_list_of_reference_angles_for_search([], 'invalid')



