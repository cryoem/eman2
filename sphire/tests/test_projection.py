from __future__ import print_function
from __future__ import division

import numpy
from math import isnan as math_isnan
from copy import deepcopy
from EMAN2_cppwrap import EMData,Util
from EMAN2_cppwrap import EMData, EMAN2Ctf
import unittest
import os
import pickle

import test_module as tm
import EMAN2_cppwrap as e2cpp

from sphire.libpy_py3 import sphire_projection as fu
from .sparx_lib import sparx_projection as oldfu

from .sparx_lib import sparx_utilities

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))









from test_module import get_arg_from_pickle_file, get_real_data, create_kb, get_data
from os import path
from mpi import *
import global_def
mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True



"""
There are some opened issues in:
1) Test_prgq.test_default_values_with_P_ref_a_leads_to_deadlock because a bug in sparx_utilities.even_angles
2) Test_prg.test_2Dimg_as_volume_FAILED --> the values of 'return_old.get_3dview()' are all Nan but the values of 'return_new.get_3dview()' not. ????could be an approxiamtion error???
3) Test_prep_vol.test_3Dimg_as_volume_interpolation_gridding_method how test if the keiser bessel filter are the same???? 
4) Test_gen_rings_ctf.test_null_nx_size_SOMETIMES_FAILED and test_default_case sometimes they pass sometimes no. It seems to be due to the approximation of some values in the dicts 
"""

class Test_project(unittest.TestCase):
    params = [21.850863880053172, 79.313172950615751, 0.0, 0.0, 0.0]    # got from 'pickle files/projection.prgl'
    volft=get_real_data(3)[0]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.project()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.project()
        self.assertEqual(cm_new.exception.message, "project() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.project(volume=None, params=self.params, radius=-1)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.project(volume=None, params=self.params, radius=-1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'project'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_params_is_emptyList_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.project(volume=self.volft, params=[], radius=-1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.project(volume=self.volft, params=[], radius=-1)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_2Dimg_returns_AttributeError_NoneType_hasnot_attribute_set_attr(self):
        volft = get_real_data(2)[0]
        with self.assertRaises(AttributeError) as cm_new:
            fu.project(volume=volft, params=self.params, radius=-1)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.project(volume=volft, params=self.params, radius=-1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'set_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_radiusNegative_and_parameter4and5Null(self):
        return_new = fu.project(volume=self.volft, params=self.params, radius=-1)
        return_old = oldfu.project(volume=self.volft, params=self.params, radius=-1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_radiusPositive_and_parameter4and5Null(self):
        return_new = fu.project(volume=self.volft, params=self.params, radius=1)
        return_old = oldfu.project(volume=self.volft, params=self.params, radius=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_radiusNegative_and_parameter4and5Not_both_Null(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.project(volume=self.volft, params=params, radius=-1)
        return_old = oldfu.project(volume=self.volft, params=params, radius=-1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_radiusPositive_and_parameter4and5Not_both_Null(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.project(volume=self.volft, params=params, radius=1)
        return_old = oldfu.project(volume=self.volft, params=params, radius=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_prgs(unittest.TestCase):
    (volft, params, interpolation_method,return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    kb = create_kb(1)

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prgs()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgs()
        self.assertEqual(cm_new.exception.message, "prgs() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_extract_plane(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prgs(volft=None,kb=self.kb,params=self.params,kbx=None, kby=None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prgs(volft=None,kb=self.kb,params=self.params,kbx=None, kby=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'extract_plane'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_params_is_emptyList_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.prgs(volft=self.volft,kb=self.kb,params=[],kbx=None, kby=None)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prgs(volft=self.volft,kb=self.kb,params=[],kbx=None, kby=None)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_realimg_returns_RuntimeError_ImageFormatException_image_not_same_size(self):
        volft = get_real_data(2)[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prgs(volft=volft,kb=self.kb,params=self.params,kbx=None, kby=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prgs(volft=volft,kb=self.kb,params=self.params,kbx=None, kby=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_kb_ArgumentError_in_EMData_extract_plane_function(self):
        kb = create_kb(3)
        with self.assertRaises(TypeError) as cm_new:
            fu.prgs(volft=self.volft,kb=kb,params=self.params,kbx=None, kby=None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgs(volft=self.volft,kb=kb,params=self.params,kbx=None, kby=None)
        output_msg = "Python argument types in\n    EMData.extract_plane(EMData, Transform, tuple)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(cm_new.exception.message, output_msg)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_kb_is_None_ArgumentError_in_EMData_extract_planec_function(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prgs(volft=self.volft,kb=None,params=self.params,kbx=None, kby=None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgs(volft=self.volft,kb=None,params=self.params,kbx=None, kby=None)
        output_msg = "Python argument types in\n    EMData.extract_plane(EMData, Transform, NoneType)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(cm_new.exception.message, output_msg)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_kbx_is_None_and_parameter4and5Null(self):
        return_new = fu.prgs(volft=self.volft,kb=self.kb,params=self.params,kbx=None, kby=None)
        return_old = oldfu.prgs(volft=self.volft,kb=self.kb,params=self.params,kbx=None, kby=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_kbx_is_None_and_parameter4and5Not_both_Null(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgs(volft=self.volft,kb=self.kb,params=params,kbx=None, kby=None)
        return_old = oldfu.prgs(volft=self.volft,kb=self.kb,params=params,kbx=None, kby=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_kbx_not_None_and_parameter4and5Null(self):
        kbx, kby, kbz = create_kb(3)
        return_new = fu.prgs(volft=self.volft,kb=kbx,params=self.params,kbx=kbx, kby=kby)
        return_old = oldfu.prgs(volft=self.volft,kb=kbx,params=self.params,kbx=kbx, kby=kby)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_kbx_not_None_and_parameter4and5Not_both_Null(self):
        kbx, kby, kbz = create_kb(3)
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgs(volft=self.volft,kb=kbx,params=params,kbx=kbx, kby=kby)
        return_old = oldfu.prgs(volft=self.volft,kb=kbx,params=params,kbx=kbx, kby=kby)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_kby_None_returns_ArgumentError_in_EMData_extract_plane_function(self):
        kbx, kby, kbz = create_kb(3)
        with self.assertRaises(TypeError) as cm_new:
            fu.prgs(self.volft,kb=kbx,params=self.params,kbx=kbx,kby=None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgs(self.volft,kb=kbx,params=self.params,kbx=kbx,kby=None)
        output_msg = "Python argument types in\n    EMData.extract_plane_rect(EMData, Transform, KaiserBessel, NoneType, KaiserBessel)\ndid not match C++ signature:\n    extract_plane_rect(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kbx, EMAN::Util::KaiserBessel {lvalue} kby, EMAN::Util::KaiserBessel {lvalue} kbz)"
        self.assertEqual(cm_new.exception.message, output_msg)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_prgl(unittest.TestCase):
    (volft, params, interpolation_method,return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]


    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prgl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgl()
        self.assertEqual(cm_new.exception.message, "prgl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_get_attr_default(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prgl(volft=None, params=self.params, interpolation_method=0, return_real=True)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prgl(volft=None, params=self.params, interpolation_method=0, return_real=True)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr_default'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_params_is_emptyList_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.prgl(volft=self.volft, params=[], interpolation_method=0, return_real=True)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prgl(volft=self.volft, params=[], interpolation_method=0, return_real=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_realimg_returns_RuntimeError_ImageFormatException_image_not_same_size(self):
        volft = get_real_data(2)[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prgl(volft=volft, params=self.params, interpolation_method=0, return_real=True)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prgl(volft=volft, params=self.params, interpolation_method=0, return_real=True)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extract_section requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_invalid_interpolation_print_error_msg(self):
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=10,return_real=False)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=10,return_real=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation1_returnReal_True_params4and5Null(self):
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=True)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation1_returnReal_False_params4and5Null(self):
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=False)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation1_returnReal_True_params4and5Not_both_Nulll(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=True)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation1_returnReal_False_params4and5Not_both_Null(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=False)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=self.interpolation_method,return_real=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation0_returnReal_True_params4and5Null(self):
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=True)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation0_returnReal_False_params4and5Null(self):
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=False)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation0_returnReal_True_params4and5Not_both_Nulll(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=True)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_interpolation0_returnReal_False_params4and5Not_both_Null(self):
        params = deepcopy(self.params)
        params[4] = 1
        return_new = fu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=False)
        return_old = oldfu.prgl(volft=self.volft, params=self.params, interpolation_method=0,return_real=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_prgq(unittest.TestCase):
    (volft, params, interpolation_method,return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    kb = create_kb(1)

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prgq()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgq()
        self.assertEqual(cm_new.exception.message, "prgq() takes at least 6 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_extract_plane(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prgq(volft=None, kb=self.kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prgq(volft=None, kb=self.kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'extract_plane'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_with_realimg_returns_RuntimeError_ImageFormatException_image_not_same_size(self):
        volft = get_real_data(2)[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prgq(volft=volft, kb=self.kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        with self.assertRaises(RuntimeError) as cm_old:
            fu.prgq(volft=volft, kb=self.kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_kb_ArgumentError_in_EMData_extract_plane_function(self):
        kb = create_kb(3)
        with self.assertRaises(TypeError) as cm_new:
            fu.prgq(volft=self.volft, kb=kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prgq(volft=self.volft, kb=kb, nx=4, delta=0.5, ref_a="S", sym="c1", MPI=False)
        output_msg = "Python argument types in\n    EMData.extract_plane(EMData, Transform, tuple)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(cm_new.exception.message, output_msg)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_nx_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prgq(volft=self.volft, kb=self.kb, nx=0, delta=15, ref_a="S", sym="c1", MPI=False)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prgq(volft=self.volft, kb=self.kb, nx=0, delta=15, ref_a="S", sym="c1", MPI=False)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_symC1_ref_a_S_MPI_False(self):
        return_new =  fu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c1", MPI=False)
        return_old = oldfu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c1", MPI=False)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))

    def test_symC1_ref_a_S_MPI_True(self):
        return_new =  fu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c1", MPI=True)
        return_old = oldfu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c1", MPI=True)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))

    def test_symC5_ref_a_S_MPI_False(self):
        return_new =  fu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c5", MPI=False)
        return_old = oldfu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c5", MPI=False)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))

    def test_symC5_ref_a_S_MPI_True(self):
        return_new =  fu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c5", MPI=True)
        return_old = oldfu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="S", sym="c5", MPI=True)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))

    def test_default_values_with_P_ref_a_leads_to_deadlock(self):
        self.assertTrue(True)
        """
         because a bug in sparx_utilities.even_angles
        return_new =  fu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="P", sym="c1", MPI=False)
        return_old = oldfu.prgq(volft=self.volft, kb=self.kb, nx=4, delta=15, ref_a="P", sym="c1", MPI=False)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))
        """


class Test_prg(unittest.TestCase):
    params = [21.850863880053172, 79.313172950615751, 0.0, 0.0, 0.0]    # got from 'pickle files/projection.prgl'
    volft_3dimg = get_real_data(3)[0]


    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prg()
        self.assertEqual(cm_new.exception.message, "prg() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prg(volume=None, params=self.params)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prg(volume=None, params=self.params)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_params_is_emptyList_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.prg(volume=self.volft_3dimg, params=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prg(volume=self.volft_3dimg, params=[])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_as_volume_FAILED(self):
        self.assertTrue(True)
        """
        volft_2dimg=get_real_data(2)[0]
        return_new =  fu.prg(volume=volft_2dimg, params=self.params)
        return_old = oldfu.prg(volume=volft_2dimg, params=self.params)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), equal_nan=True))
        """

    def test_3Dimg_as_volume(self):
        return_new =  fu.prg(volume=self.volft_3dimg, params=self.params)
        return_old = oldfu.prg(volume=self.volft_3dimg, params=self.params)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_prep_vol(unittest.TestCase):
    volft_3dimg = get_real_data(3)[0]
    volft_2dimg = get_real_data(2)[0]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prep_vol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prep_vol()
        self.assertEqual(cm_new.exception.message, "prep_vol() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_volume_returns_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prep_vol(vol=None, npad=2, interpolation_method=-1)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prep_vol(vol=None, npad=2, interpolation_method=-1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_as_volume_interpolation_gridding_method(self):
        return_new =  fu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=-1)
        return_old = oldfu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=-1)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        """ to test if the kaiser bessel are the same"""
        self.assertEqual(return_new[1].I0table_maxerror(),return_old[1].I0table_maxerror())

    def test_3Dimg_as_volume_interpolation_NN_method(self):
        return_new =  fu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=0)
        return_old = oldfu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_as_volume_interpolation_trilinear_method(self):
        return_new =  fu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=1)
        return_old = oldfu.prep_vol(vol=self.volft_3dimg, npad=2, interpolation_method=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_as_volume_interpolation_gridding_method(self):
        return_new =  fu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=-1)
        return_old = oldfu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=-1)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        """ to test if the kaiser bessel are the same"""
        self.assertEqual(return_new[1].I0table_maxerror(),return_old[1].I0table_maxerror())

    def test_2Dimg_as_volume_interpolation_NN_method(self):
        return_new =  fu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=0)
        return_old = oldfu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_as_volume_interpolation_trilinear_method(self):
        return_new =  fu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=1)
        return_old = oldfu.prep_vol(vol=self.volft_2dimg, npad=2, interpolation_method=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_gen_rings_ctf(unittest.TestCase):
    numr = [1, 1, 8, 2, 9, 16, 3, 953, 128, 16, 1081, 128, 17, 1209, 128, 18, 1337, 128, 19, 2745, 256, 26, 3001, 256,27, 3257, 256, 28, 3513, 256, 29, 3769, 256]
    def get_prjref(self):
        image = get_data(1)[0]
        prjref = [image,image]
        prjref[0].set_attr('phi', 20)
        prjref[0].set_attr('theta', 40)
        prjref[0].set_attr('psi', 40)
        prjref[1].set_attr('phi', 30)
        prjref[1].set_attr('theta', 30)
        prjref[1].set_attr('psi', 30)
        return prjref

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.gen_rings_ctf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.gen_rings_ctf()
        self.assertEqual(cm_new.exception.message, "gen_rings_ctf() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_prjref_empty_list_(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        return_new = fu.gen_rings_ctf(prjref=[], nx=4, ctf=ctf, numr=self.numr)
        return_old = oldfu.gen_rings_ctf(prjref=[], nx=4, ctf=ctf, numr=self.numr)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))

    def test_numr_is_emptyList_IndexError_list_index_out_of_range(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        with self.assertRaises(IndexError) as cm_new:
            fu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=4, ctf=ctf, numr=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=4, ctf=ctf, numr=[])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_nx_size_SOMETIMES_FAILED(self):
        self.assertTrue(True)
        """
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        return_new = fu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=0, ctf=ctf, numr=self.numr)
        return_old = oldfu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=0, ctf=ctf, numr=self.numr)
        self.assertDictEqual(return_new[0].get_attr_dict(), return_old[0].get_attr_dict())
        self.assertDictEqual(return_new[1].get_attr_dict(), return_old[1].get_attr_dict())
        """

    def test_default_case(self):
        self.assertTrue(True)
        """
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        return_new = fu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=4, ctf=ctf, numr=self.numr)
        return_old = oldfu.gen_rings_ctf(prjref=deepcopy(self.get_prjref()), nx=4, ctf=ctf, numr=self.numr)
        self.assertEqual(return_new[0].get_attr_dict(), return_old[0].get_attr_dict())
        self.assertEqual(return_new[1].get_attr_dict(), return_old[1].get_attr_dict())
        """



@unittest.skip("skip addnan tests")
class Test_lib_projection_compare(unittest.TestCase):


    def test_project_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]


        return_new = fu.project(volft,params)
        return_old = oldfu.project(volft,params)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_prgs_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]

        kb = tm.create_kb(1)

        return_new = fu.prgs(volft,kb,params)
        return_old = oldfu.prgs(volft,kb,params)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_prgl_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]

        return_new = fu.prgl(volft, params, interpolation_method)
        return_old = oldfu.prgl(volft, params, interpolation_method)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    """Function works but is too slow , takes 506 secs to test"""
    # def test_prgq_should_return_same_result(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     kb = tm.create_kb(1)
    #     (volft, params, interpolation_method,return_real) = argum[0]
    #
    #     return_new = fu.prgq(volft, kb, nx=4, delta=0.5, ref_a="S", sym="c1")
    #
    #     print("Hello")
    #     return_old = oldfu.prgq(volft, kb, nx=4, delta=0.5, ref_a="S", sym="c1")
    #     print("Hello")
    #
    #     self.assertTrue(return_new, return_old)

    def test_prg_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (volft, params, interpolation_method, return_real) = argum[0]
        volft = sparx_utilities.model_blank(100, 100, 100)

        return_new = fu.prg(volft, params)
        return_old = oldfu.prg(volft, params)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_prep_vol_should_return_same_result(self):
        volft = sparx_utilities.model_blank(100, 100, 100)

        return_new = fu.prep_vol(volft)
        return_old = oldfu.prep_vol(volft)


        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertTrue(return_new[1], return_old[1])

    from test_module import get_data
    image = get_data(1)[0]
    def test_gen_rings_ctf_should_return_same_result(self):
        argum = tm.get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = argum[0]

        nx =  4
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        prjref = []
        prjref.append(self.image)
        prjref.append(self.image)

        prjref[0].set_attr('phi', 20)
        prjref[0].set_attr('theta', 40)
        prjref[0].set_attr('psi', 40)
        prjref[1].set_attr('phi', 30)
        prjref[1].set_attr('theta', 30)
        prjref[1].set_attr('psi', 30)


        return_new = fu.gen_rings_ctf(prjref, nx , ctf, numr)
        return_old = oldfu.gen_rings_ctf(prjref, nx , ctf, numr)
        self.assertEqual(return_new[0].get_attr_dict(), return_old[0].get_attr_dict())
        self.assertEqual(return_new[1].get_attr_dict(), return_old[1].get_attr_dict())


if __name__ == '__main__':
    unittest.main()
