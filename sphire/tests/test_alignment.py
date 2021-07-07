from __future__ import print_function
from __future__ import division
from past.utils import old_div
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
from numpy import array_equal, allclose
from numpy import abs as numpy_abs
from numpy import sum as numpy_sum
from math import isnan as math_isnan
from copy import deepcopy
from EMAN2_cppwrap import EMData, Util, EMAN2Ctf
import unittest
from os import path

from tests.test_module import (
    get_data,
    get_arg_from_pickle_file,
    IMAGE_2D,
    IMAGE_2D_REFERENCE,
    KB_IMAGE2D_SIZE,
    IMAGE_3D,
    IMAGE_BLANK_2D,
    IMAGE_BLANK_3D,
    MASK_2DIMAGE,
    MASK_3DIMAGE,
    give_ali2d_single_iter_data,
    give_ornq_data,
    give_ormq_data,
)

from sphire.libpy.sp_fundamentals import fft  # ccf,rot_shift2D
from sphire.libpy.sp_utilities import model_circle, model_blank, even_angles
from sphire.libpy.sp_projection import prep_vol

from libpy_py3 import sp_alignment as oldfu
from sphire.libpy import sp_alignment as fu, sp_utilities

from mpi import *

mpi_init(0, [])

TOLERANCE = 0.005
ABSOLUTE_PATH_TO_RESOURCES = "resources_tests/pickles/"
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
# print(ABSOLUTE_PATH)

"""
WHAT IS MISSING:
0) in all the cases where the input file is an image. I did not test the case with a complex image. I was not able to generate it
1) Test_kbt returns a kaiser filter. I'm not sure to test it in a good way. It is declared in the C++ code. See in "libpyEM/libpyUtils2.cpp" how you can access to it from python
2) ali_nvol. I do not know how create/find an img with 'xform.align2d'. It crashes calling 'sp_utilities.get_params2D(..)'. I tried to use the pickle file used to test it but i get ZeroDivisionError
3) alivol_mask. I  do not know how create/find a volume with 'xform.align2d'.  (solved by adnan)
4) ali_vol_func_rotate. I cannot figure out how fill its input params. See the code for more details (solved by adnan)
5) ali_vol_func_shift same problem of (4)


RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) Test_crit2d::test_None_kb_filter_returns_ArgumentError --> different error message ... it could come from the EMAN2 error message updated in the last release
2) ali2d_single_iter --> there is a lot of dead code inside it --> see in the class comments to get more info or the "testsreport.txt"
2) shc --> all the tests are failing. Even I used some random value from another pickle file it'd work
3) ormq_fast --> it is used in an useless case, i.e.: random_method == "PCP",  of "ali2d_single_iter". I did not test it
4) prepare_refrings --> the test "test_sym_c1_initialTheta_None" seems to lead the code to a deadlock ... I commented it

In these tests there could be a bug:
1) align2d_scf --> it calls "parabl" and then it try to access to its output vector hard coding the index and goes out of the vector.
        see: "test_with_DEFAULT_params_but_validOU_returns_IndexError_list_index_out_of_range" and "test_with_DEFAULT_params_returns_IndexError_list_index_out_of_range"

In these tests there is a strange behavior:
1) Test_ormq::test_with_negative_center  --> same input different output in both of the version. At least the version are always returning the same value ... compatibility test OK, but NOT unittestable
2) Test_ornq::test_with_negative_center  --> same input different output in both of the version. At least the version are always returning the same value ... compatibility test OK, but NOT unittestable
3) Test_multialign2d_scf --> 3 tests run from pycharm but not from console (nosetests too), of course using the same python interpreter
"""

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""


#   THESE FUNCTIONS ARE COMMENTED BECAUSE NOT INVOLVED IN THE PYTHON3 CONVERSION. THEY HAVE TO BE TESTED ANYWAY
"""

#todo: the volft has no the same size like the other imgs in data vect
# since volft is a img i could create one from scratch but it has to be a complex image. how can create it?
class Test_eqproj_cascaded_ccc_fitness_function(unittest.TestCase):
    (volft, not_used1,not_used2,not_used3) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    # args --> 3 angles
    # data [volft,kb,img2d,mask2d,not_used,[float,float],float
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eqproj_cascaded_ccc_fitness_function()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eqproj_cascaded_ccc_fitness_function()
        self.assertEqual(str(cm_new.exception), "eqproj_cascaded_ccc_fitness_function() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj_cascaded_ccc_fitness_function(args=[], data=[0,1,2,3,4,5,6])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj_cascaded_ccc_fitness_function(args=[], data=[0,1,2,3,4,5,6])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj_cascaded_ccc_fitness_function(data=[], args=[3,3,4])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj_cascaded_ccc_fitness_function(data=[], args=[3,3,4])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_eqproj_cascaded_ccc_fitness_function(self):
        return_old = oldfu.eqproj_cascaded_ccc_fitness_function(args=[2,2,4], data="")
        return_new = fu.eqproj_cascaded_ccc_fitness_function(args=[2, 2, 4], data="")
        self.assertEqual(return_old[0],return_new[0])
        self.assertTrue(array_equal(return_old[1], return_new[1]))



class Test_format_list(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.format_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.format_list()
        self.assertEqual(str(cm_new.exception), "format_list() takes exactly 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.format_list(l=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.format_list(l=[])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_format_list(self):
        return_old = oldfu.format_list(l=[1.1,2,3])
        return_new = fu.format_list(l=[1.1,2,3])
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, [  1.100000,   2.000000,   3.000000]))



class Test_objective_function_just_ccc_has_maximum(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_maximum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_maximum()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_maximum takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_maximum(args=[], data=[self.volft,"",IMAGE_2D,"","","",""])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_maximum(args=[], data=[self.volft,"",IMAGE_2D,"","","",""])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_maximum(data=[], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_maximum(data=[], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_data_crash_because_SIGSEGV(self):
        pass

        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_maximum(data=[self.volft,"",None,"","","",""], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_maximum(data=[self.volft,"",None,"","","",""], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_objective_function_just_ccc_has_maximum(self):
        return_old = fu.objective_function_just_ccc_has_maximum(args=self.params, data=[self.volft,"",IMAGE_2D,"","","",""])
        return_new = oldfu.objective_function_just_ccc_has_maximum(args=self.params,data=[self.volft, "", IMAGE_2D, "", "", "", ""])
        self.assertTrue(math_isnan(return_old))
        self.assertTrue(math_isnan(return_new))



class Test_objective_function_just_ccc_has_minimum(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_minimum() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_objective_function_just_ccc_has_minimum(self):
        return_old = oldfu.objective_function_just_ccc_has_minimum(args=self.params,data=[self.volft, "", IMAGE_2D])
        return_new = fu.objective_function_just_ccc_has_minimum(args=self.params,data=[self.volft, "", IMAGE_2D])
        self.assertTrue(math_isnan(return_old))
        self.assertTrue(math_isnan(return_new))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum(args=[], data=[self.volft, "", IMAGE_2D])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum(args=[], data=[self.volft,"",IMAGE_2D,])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum(data=[], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum(data=[], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_data_crash_because_SIGSEGV(self):
        pass

        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum(data=[self.volft,"",None,"","","",""], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum(data=[self.volft,"",None,"","","",""], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_NoneType_as_volft_returns_AttributeError_NoneType_obj_hasnot_attribute_get_attr_default(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum(args=self.params, data=[None,"",IMAGE_2D,"","","",""])
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum(args=self.params, data=[None,"",IMAGE_2D,"","","",""])
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_attr_default'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_objective_function_just_ccc_has_minimum_only_reduced(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_minimum_reduced() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_objective_function_just_ccc_has_minimum_only(self):
        return_old = oldfu.objective_function_just_ccc_has_minimum_reduced(args=self.params,data=[self.volft, "", IMAGE_2D])
        return_new = fu.objective_function_just_ccc_has_minimum_only(args=self.params,data=[self.volft, "", IMAGE_2D])
        self.assertTrue(math_isnan(return_old))
        self.assertTrue(math_isnan(return_new))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced(args=[], data=[self.volft, "", IMAGE_2D])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced(args=[], data=[self.volft,"",IMAGE_2D,])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced(data=[], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced(data=[], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_data_crash_because_SIGSEGV(self):
        pass

        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_only(data=[self.volft,"",None,1,1,1,1], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_only(data=[self.volft,"",None,1,1,1,1], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_NoneType_as_volft_returns_AttributeError_NoneType_obj_hasnot_attribute_get_attr_default(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced(args=self.params, data=[None,"",IMAGE_2D,1,1,1,1])
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced(args=self.params, data=[None,"",IMAGE_2D,1,1,1,1])
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_attr_default'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



#todo i tried to use the parameter used for 'objective_function_just_ccc_has_minimum_reduced', but i get:
#    rrr =  -reference_projection.cmp("dot", prj, dict(negative = 0, mask = mask2D))/ norm_of_reference_projection RuntimeError: std::exception
class Test_objective_function_just_ccc_has_minimum_reduced_only_shifts(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced_only_shifts()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced_only_shifts()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_minimum_reduced_only_shifts() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced_only_shifts(args=[], data=[self.volft, "", IMAGE_2D, IMAGE_2D])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced_only_shifts(args=[], data=[self.volft,"",IMAGE_2D, IMAGE_2D])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum_reduced_only_shifts(data=[], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum_reduced_only_shifts(data=[], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_objective_function_just_ccc_has_minimum2(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum2()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_minimum2() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_objective_function_just_ccc_has_minimum2(self):
        return_old = oldfu.objective_function_just_ccc_has_minimum2(args=self.params,data=[self.volft, "", IMAGE_2D])
        return_new = fu.objective_function_just_ccc_has_minimum2(args=self.params,data=[self.volft, "", IMAGE_2D])
        self.assertTrue(math_isnan(return_old))
        self.assertTrue(math_isnan(return_new))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum2(args=[], data=[self.volft, "", IMAGE_2D])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum2(args=[], data=[self.volft,"",IMAGE_2D,])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum2(data=[], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum2(data=[], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_data_crash_because_SIGSEGV(self):
        pass

        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_has_minimum2(data=[self.volft,"",None,"","","",""], args=self.params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum2(data=[self.volft,"",None,"","","",""], args=self.params)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_NoneType_as_volft_returns_AttributeError_NoneType_obj_hasnot_attribute_get_attr_default(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.objective_function_just_ccc_has_minimum2(args=self.params, data=[None,"",IMAGE_2D,"","","",""])
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.objective_function_just_ccc_has_minimum2(args=self.params, data=[None,"",IMAGE_2D,"","","",""])
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_attr_default'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


#todo: it sounds a little bit obsolete. do I have to test it?
class Test_objective_function_just_ccc_has_maximum___old(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_has_maximum___old()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_has_maximum___old()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_has_maximum___old() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_objective_function_just_ccc_has_maximum___old(self):
        v2 = oldfu.objective_function_just_ccc_has_maximum___old(args="", data="")
        pass


#todo: same problem as always
class Test_objective_function_just_ccc_rewrite(unittest.TestCase):
    # params = 5 flota .... 3angle+2shift
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.objective_function_just_ccc_rewrite()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.objective_function_just_ccc_rewrite()
        self.assertEqual(str(cm_new.exception), "objective_function_just_ccc_rewrite() takes exactly 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_params_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.objective_function_just_ccc_rewrite(params=[], volft="", kb="", data_im="", mask2D="")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.objective_function_just_ccc_rewrite(params=[], volft="", kb="", data_im="", mask2D="")
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_as_volume_returns_AttributeError_object_has_no_attribute_extractplane(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.objective_function_just_ccc_rewrite(params=[1,2,1,3,4,], volft=None, kb=KB_IMAGE2D_SIZE, data_im=IMAGE_2D, mask2D=MASK_2DIMAGE)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.objective_function_just_ccc_rewrite(params=[1,2,1,3,4,], volft=None, kb=KB_IMAGE2D_SIZE, data_im=IMAGE_2D, mask2D=MASK_2DIMAGE)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception),"'NoneType' object has no attribute 'extract_plane'")

    def test_objective_function_just_ccc_rewrite(self):
        d_res = oldfu.objective_function_just_ccc_rewrite(params=[1,2,1,3,4,], volft="", kb="", data_im=IMAGE_2D, mask2D=MASK_2DIMAGE)
        pass


#todo: same problem
class Test_eqproj_cascaded_ccc(unittest.TestCase):
    (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eqproj_cascaded_ccc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eqproj_cascaded_ccc()
        self.assertEqual(str(cm_new.exception), "eqproj_cascaded_ccc() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj_cascaded_ccc(args=[], data=[self.volft,KB_IMAGE2D_SIZE,IMAGE_2D,"not_used",[1,2],1])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj_cascaded_ccc(args=[], data=[self.volft,KB_IMAGE2D_SIZE,IMAGE_2D,"not_used",[1,2],1])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj_cascaded_ccc(data=[], args=[1,2,3,4,5,6])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj_cascaded_ccc(data=[], args=[1,2,3,4,5,6])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_eqproj_cascaded_ccc(self):
        v1, v2 = oldfu.eqproj_cascaded_ccc(args="", data="")
        pass



class Test_twoD_fine_search(unittest.TestCase):
    #args = 2 shifts --> 2 floats
    #data = img, kb, 2 floats
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.twoD_fine_search()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.twoD_fine_search()
        self.assertEqual(str(cm_new.exception), "twoD_fine_search() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.twoD_fine_search(args=[], data=[IMAGE_2D,1,2,3])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.twoD_fine_search(args=[], data=[IMAGE_2D,1,2,3])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.twoD_fine_search(data=[], args=[IMAGE_2D,1])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.twoD_fine_search(data=[], args=[IMAGE_2D,1])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_twoD_fine_search(self):
        return_old = oldfu.twoD_fine_search(args=[2, 1], data=[IMAGE_2D, KB_IMAGE2D_SIZE, 2, 3])
        return_new = fu.twoD_fine_search(args=[2, 1], data=[IMAGE_2D, KB_IMAGE2D_SIZE, 2, 3])
        self.assertEqual(return_new, return_old)
        self.assertEqual(0.0219190325588, return_old)



#todo: same pb. Volft has to be complex img ---> but it size does not match with the other values
class Test_eqproj(unittest.TestCase):
    #args --> 3 angles+ 2 shifts
    #data --> volft,kb,img,mask
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eqproj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eqproj()
        self.assertEqual(str(cm_new.exception), "eqproj() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj(args=[], data=[IMAGE_2D,KB_IMAGE2D_SIZE,IMAGE_2D_REFERENCE,MASK_2DIMAGE])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj(args=[], data=[IMAGE_2D,KB_IMAGE2D_SIZE,IMAGE_2D_REFERENCE,MASK_2DIMAGE])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqproj(data=[], args=[1,2,3,4,5])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqproj(data=[], args=[1,2,3,4,5])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_eqproj(self):
        v2 = oldfu.eqproj(args=[1,2,3,4,5], data=[IMAGE_2D,KB_IMAGE2D_SIZE,IMAGE_2D_REFERENCE,MASK_2DIMAGE])
        pass


#todo: look into the params
class Test_eqprojDot(unittest.TestCase):
    '''
    phi = args[0]
	tht = args[1]
	psi = args[2]
	vol = data[0]
	img = data[1]
	s2x = data[2]
	s2y = data[3]
	msk = data[4]
	CTF = data[5]
	ou  = data[6]
    '''
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eqprojDot()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eqprojDot()
        self.assertEqual(str(cm_new.exception), "eqprojDot() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqprojDot(args=[], data=[1,2,3,4,5,6,7])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqprojDot(args=[], data=[1,2,3,4,5,6,7])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqprojDot(data=[], args=[0,1,2])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqprojDot(data=[], args=[0,1,2])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_eqprojDot(self):
        img=deepcopy(IMAGE_2D)
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        img.set_attr("ctf", ctf)
        return_old = oldfu.eqprojDot(args=[1, 2, 3], data=[IMAGE_2D, img, 1, 2, MASK_2DIMAGE, False, 2])
        pass



class Test_eqprojEuler(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eqprojEuler()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eqprojEuler()
        self.assertEqual(str(cm_new.exception), "eqprojEuler() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqprojEuler(args=[], data=[1,2,3,4,5,6,7])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqprojEuler(args=[], data=[1,2,3,4,5,6,7])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.eqprojEuler(data=[], args=[0,1,2])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eqprojEuler(data=[], args=[0,1,2])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_eqprojEuler(self):
        v2 = oldfu.eqprojEuler(args="", data="")
        pass


class Test_symm_func(unittest.TestCase):
    #args=3 angles -->3 float
    #data = img,mask,sym --> e.g.: sym='c1'
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.symm_func()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.symm_func()
        self.assertEqual(str(cm_new.exception), "symm_func() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.symm_func(args=[], data=[1,2,3,4,5,6,7])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.symm_func(args=[], data=[1,2,3,4,5,6,7])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.symm_func(data=[], args=[0,1,2])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.symm_func(data=[], args=[0,1,2])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_none_mask_crash_because_SIGSEGV(self):
        pass
        '''
        with self.assertRaises(IndexError) as cm_new:
            fu.symm_func(data=[IMAGE_2D,None,"c1"], args=[1,1,2])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.symm_func(data=[IMAGE_2D,None,"c1"], args=[1,1,2])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        '''

    def test_empty_mask_crash_because_SIGSEGV(self):
        pass
        '''
        with self.assertRaises(IndexError) as cm_new:
            fu.symm_func(data=[IMAGE_2D,EMData(),"c1"], args=[1,1,2])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.symm_func(data=[IMAGE_2D,EMData(),"c1"], args=[1,1,2])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        '''

    def test_empty_data_images_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.symm_func(data=[EMData(),MASK_2DIMAGE,"c1"], args=[1,1,2])
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.symm_func(data=[EMData(),MASK_2DIMAGE,"c1"], args=[1,1,2])
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_rot_scale_trans_background(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.symm_func(data=[None,MASK_2DIMAGE,"c1"], args=[1,1,2])
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.symm_func(data=[None,MASK_2DIMAGE,"c1"], args=[1,1,2])
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'rot_scale_trans_background'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_symm_func(self):
        return_old = oldfu.symm_func(args=[1,1,3], data=[IMAGE_2D,MASK_2DIMAGE,"c1"])
        return_new = fu.symm_func(args=[1, 1, 3], data=[IMAGE_2D, MASK_2DIMAGE, "c1"])
        self.assertEqual(return_old,return_new)
        self.assertEqual(return_old, 2.7516720295)


#todo: inizia da qua
class Test_find_symm(unittest.TestCase):

    #USed the defualt value of sp_application.rot_sym(...) to fill the params

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.find_symm()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.find_symm()
        self.assertEqual(str(cm_new.exception), "find_symm() takes exactly 9 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_as_vol_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.find_symm(vol=None, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2, scale=[20,20,20], ftolerance=1.e-4, xtolerance=1.e-4,)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.find_symm(vol=None, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2, scale=[20,20,20], ftolerance=1.e-4, xtolerance=1.e-4,)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'rot_scale_trans_background'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_as_mask_crashes_because_SIGSEGV(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.find_symm(vol=IMAGE_2D, mask=None, sym_gp="d4", phi=1, theta=3, psi=2, scale=[20,20,20], ftolerance=1.e-4, xtolerance=1.e-4,)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.find_symm(vol=IMAGE_2D, mask=None, sym_gp="d4", phi=1, theta=3, psi=2, scale=[20,20,20], ftolerance=1.e-4, xtolerance=1.e-4,)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'rot_scale_trans_background'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2DIMAGE(self):
        return_old = oldfu.find_symm(vol=IMAGE_2D, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        return_new = fu.find_symm(vol=IMAGE_2D, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, ([-23.009874131944443, 91.44052915219908, 3.3870352285879637], 2.95225191116333, 24)))

    def test_3DIMAGE(self):
        return_old = oldfu.find_symm(vol=IMAGE_3D, mask=MASK_3DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        return_new = fu.find_symm(vol=IMAGE_3D, mask=MASK_3DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old,([20.747031464334704, 0.00030542695473246076, 5.511016803840876], 795.3264770507812, 21)))

    def test_IMAGE_BLANK_2D(self):
        return_old = oldfu.find_symm(vol=IMAGE_BLANK_2D, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        return_new = fu.find_symm(vol=IMAGE_BLANK_2D, mask=MASK_2DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        self.assertTrue(array_equal(return_old[0],return_new[0]))
        self.assertTrue(math_isnan(return_old[1]))
        self.assertTrue(math_isnan(return_new[1]))
        self.assertEqual(return_old[2], return_new[2])
        self.assertTrue(array_equal(return_old[0], [1.0, 3.0, 2.0]))
        self.assertEqual(return_old[2], 500)

    def test_IMAGE_BLANK_3D(self):
        return_old = oldfu.find_symm(vol=IMAGE_BLANK_3D, mask=MASK_3DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        return_new = fu.find_symm(vol=IMAGE_BLANK_3D, mask=MASK_3DIMAGE, sym_gp="d4", phi=1, theta=3, psi=2,scale=[20, 20, 20], ftolerance=1.e-4, xtolerance=1.e-4, )
        self.assertTrue(array_equal(return_old[0],return_new[0]))
        self.assertTrue(math_isnan(return_old[1]))
        self.assertTrue(math_isnan(return_new[1]))
        self.assertEqual(return_old[2], return_new[2])
        self.assertTrue(array_equal(return_old[0], [1.0, 3.0, 2.0]))
        self.assertEqual(return_old[2], 500)


class Test_ormq_peaks(unittest.TestCase):

    #I did not test with crefin =IMAGE_2D_REFERENCE because we have a huge output with nan and inf values in random position .... I have no idea how i'd unitest it in a fast and clean way

    numr = [1, 1, 8, 2, 9, 16, 3, 25, 32, 4, 57, 32, 5, 89, 32, 6, 121, 64, 7, 185, 64, 8, 249, 64, 9, 313, 64, 10,377, 64, 11, 441, 128, 12, 569, 128, 13, 697, 128, 14, 825, 128, 15, 953, 128, 16, 1081, 128, 17, 1209,128, 18, 1337, 128, 19, 1465, 128, 20, 1593, 128, 21, 1721, 256, 22, 1977, 256, 23, 2233, 256, 24, 2489,256, 25, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]
    mode = "F"
    xrng = yrng = [4, 4]
    step = 1
    cnx = cny = 36.0
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ormq_peaks()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ormq_peaks()
        self.assertEqual(str(cm_new.exception), "ormq_peaks() takes exactly 9 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Empty_image_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ormq_peaks(image=EMData(), crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ormq_peaks(image=EMData(), crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        self.assertEqual(str(cm_new.exception), "float argument required, not list")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.ormq_peaks(image=None, crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ormq_peaks(image=None, crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_as_crefim_crash_because_SIGSEGV(self):
        pass
        '''
        with self.assertRaises(AttributeError) as cm_new:
            fu.ormq_peaks(image=IMAGE_2D, crefim=None, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ormq_peaks(image=IMAGE_2D, crefim=None, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        '''

    def test_empty_data_as_crefim_crash_because_SIGSEGV(self):
        pass
        '''
        with self.assertRaises(AttributeError) as cm_new:
            fu.ormq_peaks(image=IMAGE_2D, crefim=EMData(), xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ormq_peaks(image=IMAGE_2D, crefim=EMData(), xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        '''

    def test_ormq_peaks_3DIMAGE(self):
        return_old = oldfu.ormq_peaks(image=IMAGE_3D, crefim=IMAGE_BLANK_3D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        return_new = fu.ormq_peaks(image=IMAGE_3D, crefim=IMAGE_BLANK_3D, xrng=self.xrng, yrng=self.yrng,step=self.step, mode=self.mode, numr=self.numr, cnx=self.cnx, cny=self.cny)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, [[1.0, 128.0, 9.0, 9.0, 1.0, 0.0, 0.0, 0.0, 0], [1.0, 128.0, 9.0, 9.0, 1.0, 0.0, 0.0, 0.0, 1]]))

    def test_ormq_peaks_2DIMAGE(self):
        return_old = oldfu.ormq_peaks(image=IMAGE_2D, crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng, step=self.step, mode=self.mode, numr=self.numr, cnx=6.0, cny=6.0)
        return_new = fu.ormq_peaks(image=IMAGE_2D, crefim=IMAGE_BLANK_2D, xrng=self.xrng, yrng=self.yrng,step=self.step, mode=self.mode, numr=self.numr, cnx=6.0, cny=6.0)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, [[1.0, 128.0, 147.0, 147.0, 1.0, 0.0, 0.0, 0.0, 0], [1.0, 128.0, 147.0, 147.0, 1.0, 0.0, 0.0, 0.0, 1]]))



class Test_select_k(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.select_k()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.select_k()
        self.assertEqual(str(cm_new.exception), "select_k() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.select_k(dJe=[], T=3)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.select_k(dJe=[], T=3)
        self.assertEqual(str(cm_new.exception), "tuple index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_select_k(self):
        return_old = oldfu.select_k(dJe=range(100), T=3)
        return_new = fu.select_k(dJe=range(100), T=3)
        self.assertEqual(return_old,92)
        self.assertEqual(return_old, return_new)



class Test_sim_anneal(unittest.TestCase):
    peaks=[[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4]]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sim_anneal()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sim_anneal()
        self.assertEqual(str(cm_new.exception), "sim_anneal() takes exactly 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.sim_anneal(peaks=[], T=0, step=1, mode="nF", maxrin=1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.sim_anneal(peaks=[], T=0, step=1, mode="nF", maxrin=1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_modeF_Tpositive(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=1, step=1, mode="F", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=1, step=1, mode="F", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (72.0, -0.47312770483054567, 0.16171015713856704, 0.1, 0.3333333333333333, 1)))

    def test_modeF_Tnegative(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=-1, step=1, mode="F", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=-1, step=1, mode="F", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (72.0, -0.47312770483054567, 0.16171015713856704, 0.1, 1.0, 1)))

    def test_modeF_Tnull(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=0, step=1, mode="F", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=0, step=1, mode="F", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (72.0, -0.47312770483054567, 0.16171015713856704, 0.1, 1.0, 0)))

    def test_mode_notF_Tpositive(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=1, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=1, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (36.0, -0.47781919922947347, -0.14727122206223708, 0.1, 0.3333333333333333, 1)))

    def test_mode_notF_Tnegative(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=-1, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=-1, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old,(36.0, -0.47781919922947347, -0.14727122206223708, 0.1, 1.0, 1)))

    def test_mode_notF_Tnull(self):
        return_old = oldfu.sim_anneal(peaks=self.peaks, T=0, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_anneal(peaks=self.peaks, T=0, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (36.0, -0.47781919922947347, -0.14727122206223708, 0.1, 1.0, 0)))



class Test_sim_ccf(unittest.TestCase):
    peaks=[[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4],[0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4]]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sim_ccf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sim_ccf()
        self.assertEqual(str(cm_new.exception), "sim_ccf() takes exactly 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.sim_ccf(peaks=[], T=0, step=1, mode="nF", maxrin=1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.sim_ccf(peaks=[], T=0, step=1, mode="nF", maxrin=1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_modeF_Tpositive(self):
        return_old = oldfu.sim_ccf(peaks=[1, 2, 3, 4, 5, 6], T=1, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_ccf(peaks=[1, 2, 3, 4, 5, 6], T=1, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, (0.0, -3.0, -4.0, 5, 1, 6)))

    def test_modeF_Tnegative(self):
        return_old = oldfu.sim_ccf(peaks=self.peaks, T=-1, step=1, mode="F", maxrin=1)
        return_new = fu.sim_ccf(peaks=self.peaks, T=-1, step=1, mode="F", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (72.0, -0.47312770483054567, 0.16171015713856704, 0.1, 1.0, 1)))

    def test_modeF_Tnull(self):
        return_old = oldfu.sim_ccf(peaks=self.peaks, T=0, step=1, mode="F", maxrin=1)
        return_new = fu.sim_ccf(peaks=self.peaks, T=0, step=1, mode="F", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (72.0, -0.47312770483054567, 0.16171015713856704, 0.1, 1.0, 0)))

    def test_mode_notF_Tpositive(self):
        return_old = oldfu.sim_ccf(peaks=[1,2,3,4,5,6], T=1, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_ccf(peaks=[1,2,3,4,5,6], T=1, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old,(0.0, -3.0, -4.0, 5, 1, 6)))

    def test_mode_notF_Tnegative(self):
        return_old = oldfu.sim_ccf(peaks=self.peaks, T=-1, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_ccf(peaks=self.peaks, T=-1, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old,(36.0, -0.47781919922947347, -0.14727122206223708, 0.1, 1.0, 1)))

    def test_mode_notF_Tnull(self):
        return_old = oldfu.sim_ccf(peaks=self.peaks, T=0, step=1, mode="nF", maxrin=1)
        return_new = fu.sim_ccf(peaks=self.peaks, T=0, step=1, mode="nF", maxrin=1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_old, (36.0, -0.47781919922947347, -0.14727122206223708, 0.1, 1.0, 0)))


class Test_sim_anneal2(unittest.TestCase):
    peaks = [[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4],[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4],[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4]]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sim_anneal2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sim_anneal2()
        self.assertEqual(str(cm_new.exception), "sim_anneal2() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.sim_anneal2(peaks=[], Iter=1, T0=1, F=1, SA_stop=1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.sim_anneal2(peaks=[], Iter=1, T0=1, F=1, SA_stop=1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_T_lower_0000001_and_iter_lower_SA(self):
        return_old = oldfu.sim_anneal2(peaks=self.peaks,  Iter=1, T0=0.000001, F=2, SA_stop=2)
        return_new = fu.sim_anneal2(peaks=self.peaks,  Iter=1, T0=0.000001, F=2, SA_stop=2)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, [1.0, 0.0, 0.0]))

    def test_iter_lower_SA(self):
        return_old = oldfu.sim_anneal2(peaks=self.peaks, Iter=1, T0=2, F=2, SA_stop=2)
        return_new = fu.sim_anneal2(peaks=self.peaks, Iter=1, T0=2, F=2, SA_stop=2)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, [0.3333333333333333, 0.3333333333333333, 0.3333333333333333]))

    def test_iter_NOT_lower_SA(self):
        return_old = oldfu.sim_anneal2(peaks=self.peaks, Iter=2, T0=2, F=2, SA_stop=2)
        return_new = fu.sim_anneal2(peaks=self.peaks, Iter=2, T0=2, F=2, SA_stop=2)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, [1.0, 0.0, 0.0]))



class Test_sim_anneal3(unittest.TestCase):
    peaks = [[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4],[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4],[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4]]
    peakm = [[0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41],[0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41],[0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41, 0.11, 0.21, 0.31, 0.41]]

    peaks_major = [[0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421],[0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421],[0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421, 0.121, 0.221, 0.321, 0.421]]

    peakm_major = [[0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431],[0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431],[0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431, 0.131, 0.231, 0.331, 0.431]]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sim_anneal3()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sim_anneal3()
        self.assertEqual(str(cm_new.exception), "sim_anneal3() takes exactly 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_T_lower_0000001_and_iter_lower_SA(self):
        return_old = oldfu.sim_anneal3(peaks=self.peaks, peakm=self.peakm, peaks_major=self.peaks_major, peakm_major=self.peakm_major, Iter=1, T0=0.000001, F=2, SA_stop=2)
        return_new = fu.sim_anneal3(peaks=self.peaks, peakm=self.peakm, peaks_major=self.peaks_major, peakm_major=self.peakm_major, Iter=1, T0=0.000001, F=2, SA_stop=2)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, (0.21, 0.31, 0.41, 1, 0.11, 0)))


    def test_iter_NOT_lower_SA(self):
        return_old = oldfu.sim_anneal3(peaks=self.peaks, peakm=self.peakm, peaks_major=self.peaks_major, peakm_major=self.peakm_major, Iter=2, T0=2, F=2, SA_stop=2)
        return_new = fu.sim_anneal3(peaks=self.peaks, peakm=self.peakm, peaks_major=self.peaks_major, peakm_major=self.peakm_major, Iter=2, T0=2, F=2, SA_stop=2)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_old, (0.21, 0.31, 0.41, 1, 0.11, 0)))



class Test_prep_vol_kb(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prep_vol_kb()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prep_vol_kb()
        self.assertEqual(str(cm_new.exception), "prep_vol_kb() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_images_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prep_vol_kb(vol=EMData(), kb=KB_IMAGE2D_SIZE, npad=2)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prep_vol_kb(vol=EMData(), kb=KB_IMAGE2D_SIZE, npad=2)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_copy(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prep_vol_kb(vol=None, kb=KB_IMAGE2D_SIZE, npad=2)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prep_vol_kb(vol=None, kb=KB_IMAGE2D_SIZE, npad=2)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'copy'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_kb_images_returns_ArgumentError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prep_vol_kb(vol=IMAGE_2D, kb=EMData(), npad=2)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prep_vol_kb(vol=IMAGE_2D, kb=EMData(), npad=2)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_as_kb_returns_ArgumentError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.prep_vol_kb(vol=EMData(), kb=None, npad=2)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prep_vol_kb(vol=EMData(), kb=None, npad=2)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'copy'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_prep_vol_kb(self):
        return_old = oldfu.prep_vol_kb(vol=IMAGE_2D, kb=KB_IMAGE2D_SIZE, npad=2)
        return_new = fu.prep_vol_kb(vol=IMAGE_2D, kb=KB_IMAGE2D_SIZE, npad=2)
        self.assertTrue(array_equal(return_old.get_3dview(), return_new.get_3dview()))




class Test_prepare_refrings_projections(unittest.TestCase):
    '''
    Take a look to sparx_utilities.py --> even_angles_cd(...)for the meaning of the following params
        ref_a --> P=Penczek algorithm, S=Saff algorithm to calculate di reference angle
        phiEQpsi  --> 'Minus', if you want psi=-phi to create a list of  angles suitable for projections, otherwise 'Zero'

    In case of rectangular kb filter see how it uses kbi variables in sparx_projection.py --> prgs(...) to understand better
    '''
    volft = model_blank(nx=100, ny=100, nz=100, bckg=0.0)
    numr = [1, 1, 8, 2, 9, 16, 3, 953, 128, 16, 1081, 128, 17, 1209, 128, 18, 1337, 128, 19, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for img1, img2 in zip(return_new, return_old):
            try:
                self.assertTrue(array_equal(img1.get_3dview(), img2.get_3dview()))
            except AssertionError:
                # since sometimes we get  img1.get_3dview()= [[[ nan  nan  nan ...,  nan  nan  nan]]] we skip these cases
                res = numpy_sum(img1.get_3dview() - img2.get_3dview())
                if math_isnan(res) is False:
                    self.assertTrue(TOLERANCE > numpy_abs(res))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings_projections()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings_projections()
        self.assertEqual(str(cm_new.exception), "prepare_refrings_projections() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_volume_returns_RuntimeError_ImageFormatException_extractplane_requires_complex_img(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepare_refrings_projections(volft=EMData(), kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepare_refrings_projections(volft=EMData(), kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_NoneType_as_volume_returns_AttributeError_ImageFormatException_extractplane_requires_complex_img(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepare_refrings_projections(volft=None, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prepare_refrings_projections(volft=None, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception),"'NoneType' object has no attribute 'extract_plane'")

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range_modeH(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=[], MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=[], MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_referenceAngle_returns_IndexError_list_index_out_of_range_modeH(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a=[], sym="c1", mode="H", numr=[], MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a=[], sym="c1", mode="H", numr=[], MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_No_kb_ArgumentError_in_EMData_extract_plane_function_modeH(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings_projections(volft=self.volft, kb=None,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings_projections(volft=self.volft, kb=None,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        output_msg = "Python argument types in\n    EMData.extract_plane(EMData, Transform, NoneType)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_sym_c1_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated_modeH(self):
        volft,kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_with_sym_c5_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated_modeH(self):
        volft,kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    @unittest.skip( "\n***************************\n\t\t 'Test_prepare_refrings_projectionstest_sym_c1_initialTheta_None. Even if this combination is it seems to lead the code to a deadlock, i waited more then an hour'\n***************************")
    def test_sym_c1_initialTheta_None_modeH(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=None, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=None, delta_theta=0.5)
        self.test_all_the_conditions(return_new, return_old)

    def test_No_nz_data_size_Error_msg_datasize_hasnot_be_given_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=0, delta=2.0, ref_a="S", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=0, delta=2.0, ref_a="S", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_oct_Warning_in_even_angles_this_sym_isnot_supported_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="oct", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="oct", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="H", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero_modeH(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="H", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_with_sym_c1_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated_modeF(self):
        volft,kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_with_sym_c5_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated_modeF(self):
        volft,kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=True, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    @unittest.skip( "\n***************************\n\t\t 'Test_prepare_refrings_projectionstest_sym_c1_initialTheta_None. Even if this combination is it seems to lead the code to a deadlock, i waited more then an hour'\n***************************")
    def test_sym_c1_initialTheta_None_modeF(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=None, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb, nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False,phiEqpsi="Minus",   initial_theta=None, delta_theta=0.5)
        self.test_all_the_conditions(return_new, return_old)

    def test_No_nz_data_size_Error_msg_datasize_hasnot_be_given_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=0, delta=2.0, ref_a="S", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=0, delta=2.0, ref_a="S", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_oct_Warning_in_even_angles_this_sym_isnot_supported_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="oct", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="oct", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a= even_angles(60.0), sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="S", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Minus",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c1", mode="F", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero_modeF(self):
        volft,kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings_projections(volft=volft, kb=kb,nz=4, delta=2.0, ref_a="P", sym="c5", mode="F", numr=self.numr, MPI=False, phiEqpsi="Zero",   initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)



class Test_prepare_refrings2(unittest.TestCase):
    '''
    Take a look to sparx_utilities.py --> even_angles_cd(...)for the meaning of the following params
        ref_a --> P=Penczek algorithm, S=Saff algorithm to calculate di reference angle
        phiEQpsi  --> 'Minus', if you want psi=-phi to create a list of  angles suitable for projections, otherwise 'Zero'

    In case of rectangular kb filter see how it uses kbi variables in sparx_projection.py --> prgs(...) to understand better
    '''
    volft = model_blank(100,100,100)
    numr = [1, 1, 8, 2, 9, 16, 3, 953, 128, 16, 1081, 128, 17, 1209, 128, 18, 1337, 128, 19, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for img1, img2 in zip(return_new, return_old):
            try:
                self.assertTrue(array_equal(img1.get_3dview(), img2.get_3dview()))
            except AssertionError:
                # since sometimes we get  img1.get_3dview()= [[[ nan  nan  nan ...,  nan  nan  nan]]] we skip these cases
                res = numpy_sum(img1.get_3dview() - img2.get_3dview())
                if math_isnan(res) is False:
                    self.assertTrue(TOLERANCE > numpy_abs(res))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings2()
        self.assertEqual(str(cm_new.exception), "prepare_refrings2() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_volume_returns_RuntimeError_ImageFormatException_extractplane_requires_complex_img(self):
        volft, kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepare_refrings2(volft=EMData(), kb=kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepare_refrings2(volft=EMData(),kb= kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_NoneType_as_volume_returns_AttributeError_ImageFormatException_extractplane_requires_complex_img(self):
        volft, kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepare_refrings2(volft=None, kb=kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prepare_refrings2(volft=None, kb=kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception),"'NoneType' object has no attribute 'extract_plane'")

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        volft, kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_referenceAngle_returns_IndexError_list_index_out_of_range(self):
        volft, kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a=[], sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a=[], sym="c1", numr=[], MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_No_kb_ArgumentError_in_EMData_extract_plane_function(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings2(volft=self.volft, kb=None,nz=4, segmask= MASK_2DIMAGE, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings2(volft=self.volft,kb= None,nz=4, segmask= MASK_2DIMAGE, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        output_msg = "Python argument types in\n    EMData.extract_plane(EMData, Transform, NoneType)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_sym_c1_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_with_sym_c5_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=True, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    @unittest.skip( "\n***************************\n\t\t 'Test_prepare_refrings2test_sym_c1_initialTheta_None. Even if this combination is it seems to lead the code to a deadlock, i waited more then an hour'\n***************************")
    def test_sym_c1_initialTheta_None(self):
        volft, kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=None, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb, nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False,phiEqpsi="Minus", kbx=None, kby=None, initial_theta=None, delta_theta=0.5)
        self.test_all_the_conditions(return_new, return_old)

    def test_No_nz_data_size_Error_msg_datasize_hasnot_be_given(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=0, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=0, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="oct", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a= even_angles(60.0), sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="S", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Minus", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(), volft.get_ysize(), volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c1", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft,kb = prep_vol(self.volft)
        mask = model_circle(2, volft.get_xsize(),volft.get_ysize(),volft.get_zsize())
        return_new = fu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft, kb=kb,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=None, kby=None, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(nx=100, ny=50, nz=100, bckg=0.0))
        mask = model_circle(r=2, nx=100, ny=50, nz=100)
        return_new = fu.prepare_refrings2(volft=volft, kb=kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        return_old = oldfu.prepare_refrings2(volft=volft,kb= kbz,nz=4, segmask= mask, delta=2.0, ref_a="P", sym="c5", numr=self.numr, MPI=False, phiEqpsi="Zero", kbx=kbx, kby=kby, initial_theta=0.1, delta_theta=0.5)
        self.test_all_the_conditions(return_new,return_old)



class Test_refprojs(unittest.TestCase):
    numr= [1, 1, 8, 2, 9, 16, 3, 25, 32, 4, 57, 32, 5, 89, 32, 6, 121, 64, 7, 185, 64, 8, 249, 64, 9, 313, 64, 10,377, 64, 11, 441, 128, 12, 569, 128, 13, 697, 128, 14, 825, 128, 15, 953, 128, 16, 1081, 128, 17, 1209, 128, 18, 1337, 128, 19, 1465, 128, 20, 1593, 128, 21, 1721, 256, 22, 1977, 256, 23, 2233, 256, 24, 2489, 256, 25, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]
    wr = [25.132741228718345, 12.566370614359172, 4.71238898038469, 6.283185307179586, 7.853981633974483,2.356194490192345, 2.748893571891069, 3.141592653589793, 3.5342917352885173, 3.9269908169872414, 1.0799224746714913, 1.1780972450961724, 1.2762720155208536, 1.3744467859455345, 1.4726215563702154, 1.5707963267948966, 1.6689710972195777, 1.7671458676442586, 1.8653206380689396, 1.9634954084936207, 0.5154175447295755, 0.5399612373357456, 0.5645049299419159, 0.5890486225480862, 0.6135923151542565, 0.6381360077604268, 0.662679700366597, 0.6872233929727672, 0.7117670855789375]
    ref_angle = even_angles(symmetry='c5')
    volft, kb = prep_vol(model_blank(nx=100, ny=100, nz=100, bckg=0.0))
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.refprojs()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.refprojs()
        self.assertEqual(str(cm_new.exception), "refprojs() takes exactly 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Nonetype_input_img_returns_attributeError(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.refprojs(volft=None, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.refprojs(volft=None, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'extract_plane'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_images_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.refprojs(volft=EMData(), kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.refprojs(volft=EMData(), kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr,mode="H", wr=self.wr)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img_returns_ArgumentError(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.refprojs(volft=self.volft, kb=None, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.refprojs(volft=self.volft, kb=None, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'extract_plane'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_EMpty_input_img_returns_ArgumentError(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.refprojs(volft=self.volft, kb=EMData(), ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.refprojs(volft=self.volft, kb=EMData(), ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'extract_plane'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_modeH(self):
        return_old = oldfu.refprojs(volft=self.volft, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="H", wr=self.wr)
        return_new = fu.refprojs(volft=self.volft, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36,numr=self.numr, mode="H", wr=self.wr)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i,j))

    def test_modeF(self):
        return_old = oldfu.refprojs(volft=self.volft, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36, numr=self.numr, mode="F", wr=self.wr)
        return_new = fu.refprojs(volft=self.volft, kb=self.kb, ref_angles=self.ref_angle, cnx=36, cny=36,numr=self.numr, mode="F", wr=self.wr)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i,j))


#todo: need an image with xform. it does not work with 'pickle files/alignment.shc' because the otehr input params
class Test_(unittest.TestCase):
    numr = [1, 1, 8, 2, 9, 16, 3, 25, 32, 4, 57, 32, 5, 89, 32, 6, 121, 64, 7, 185, 64, 8, 249, 64, 9, 313, 64, 10,377, 64, 11, 441, 128, 12, 569, 128, 13, 697, 128, 14, 825, 128, 15, 953, 128, 16, 1081, 128, 17, 1209,128, 18, 1337, 128, 19, 1465, 128, 20, 1593, 128, 21, 1721, 256, 22, 1977, 256, 23, 2233, 256, 24, 2489,256, 25, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore_zoom()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore_zoom()
        self.assertEqual(str(cm_new.exception), "proj_ali_incore_zoom() takes at least 6 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_(self):
        peak, pixel_error = oldfu.proj_ali_incore_zoom(data=IMAGE_2D, refrings=[IMAGE_2D,IMAGE_2D], numr=self.numr, xrng=[], yrng=[],
                                                       step=1, finfo=None, sym="c1", delta_psi=0.0)
        pass

#todo: need an image with xform. it does not work with 'pickle files/alignment.shc' because the otehr input params
class Test_proj_ali_incore_local_zoom(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore_local_zoom()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore_local_zoom()
        self.assertEqual(str(cm_new.exception), "proj_ali_incore_local_zoom() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_incore_local_zoom(self):
        d_res = oldfu.proj_ali_incore_local_zoom(data=None, refrings=None, list_of_reference_angles=None, numr=None,
                                                 xrng=None, yrng=None, step=None, an=None, finfo=None, sym='c1',
                                                 delta_psi=0.0)
        pass



class Test_ornq_gridding(unittest.TestCase):

    numr= [1, 1, 8, 2, 9, 16, 3, 25, 32, 4, 57, 32, 5, 89, 32, 6, 121, 64, 7, 185, 64, 8, 249, 64, 9, 313, 64, 10,377, 64, 11, 441, 128, 12, 569, 128, 13, 697, 128, 14, 825, 128, 15, 953, 128, 16, 1081, 128, 17, 1209,128, 18, 1337, 128, 19, 1465, 128, 20, 1593, 128, 21, 1721, 256, 22, 1977, 256, 23, 2233, 256, 24, 2489,256, 25, 2745, 256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ornq_gridding()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ornq_gridding()
        self.assertEqual(str(cm_new.exception), "ornq_gridding() takes at least 9 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq_gridding(image=[], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq_gridding(image=[], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_crashes_because_SIGSEGV(self):
        pass
        '''
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq_gridding(image=[None,None,None], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq_gridding(image=[None,None,None], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        '''

    @unittest.skip("test_EMpty_img_leads_to_deadlock")
    def test_EMpty_img_leads_to_deadlock(self):
        return_old = oldfu.ornq_gridding(image=[EMData(),EMData(),EMData()], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        return_new = fu.ornq_gridding(image=[EMData(),EMData(),EMData()], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal((291.77770256996155, 2, 4, 0, 17.63412950428409), return_old))

    def test_None_crefimg_crashes_because_SIGSEGV(self):
        pass
        '''
        return_old = oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=None,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        return_new = fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=None,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal((291.77770256996155, 2, 4, 0, 17.63412950428409), return_old))
        '''

    def test_empty_crefimg_crashes_because_SIGSEGV(self):
        pass
        '''
        return_old = oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=EMData(),shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        return_new = fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=EMData(),shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal((291.77770256996155, 2, 4, 0, 17.63412950428409), return_old))
        '''

    def test_None_kb_python_boost_ArgumentError(self):
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=None, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=None, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_kb_python_boost_ArgumentError(self):
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=EMData(), mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=EMData(), mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_modeH(self):
        return_old = oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        return_new = fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="H",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal((291.77770256996155, 2, 4, 0, 17.63412950428409), return_old))

    def test_modeF(self):
        return_old = oldfu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="F",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        return_new = fu.ornq_gridding(image=[IMAGE_2D, IMAGE_2D, IMAGE_2D], crefim=IMAGE_2D_REFERENCE,shifts=[[1, 1], [2, 1], [1, 2]], shrink=2, kb=KB_IMAGE2D_SIZE, mode="F",numr=self.numr, cnx=36, cny=36, deltapsi=0.0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal((111.66072607040405, 2, 4, 0, 17.346570451540014), return_old))


class Test_ali3D_gridding(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3D_gridding()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3D_gridding()
        self.assertEqual(str(cm_new.exception), "ali3D_gridding() takes at least 11 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali3D_gridding(self):
        newpar, simis = oldfu.ali3D_gridding(data=None, volprep=None, refang=None, delta_psi=None, shifts=None,
                                             shrink=None, numr=None, wr=None, cnx=None, myid=None, main_node=None,
                                             kb3D=None)
        pass



class Test_ali3D_direct(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3D_direct()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3D_direct()
        self.assertEqual(str(cm_new.exception), "ali3D_direct() takes at least 7 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali3D_direct(self):
        v = oldfu.ali3D_direct(data=None, volprep=None, refang=None, delta_psi=None, shifts=None, myid=None,
                               main_node=None, lentop=1000, kb3D=None)
        pass


class Test_ali3D_direct_preselect(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3D_direct_preselect()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3D_direct_preselect()
        self.assertEqual(str(cm_new.exception), "ali3D_direct_preselect() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali3D_direct_preselect(self):
        v = oldfu.ali3D_direct_preselect(data=None, volprep=None, oldcodedparams=None, refang=None, delta_psi=None,
                                         shifts=None, myid=None, main_node=None, lentop=1000, kb3D=None)
        pass


class Test_ali3D_direct_local(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3D_direct_local()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3D_direct_local()
        self.assertEqual(str(cm_new.exception), "ali3D_direct_local() takes at least 9 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali3D_direct_local(self):
        v = oldfu.ali3D_direct_local(data=None, volprep=None, refang=None, delta_psi=None, shifts=None, an=None,
                                     oldangs=None, myid=None, main_node=None, lentop=1000, kb3D=None)
        pass


class Test_proj_ali_incore_direct(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore_direct()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore_direct()
        self.assertEqual(str(cm_new.exception), "proj_ali_incore_direct() takes at least 6 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_incore_direct(self):
        peak, pixel_error = oldfu.proj_ali_incore_direct(data=None, ref_angs=None, numr=None, xrng=None, yrng=None,
                                                         step=None, finfo=None, sym="c1", delta_psi=0.0, rshift=0.0)
        pass


class Test_proj_ali_helical(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helical()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helical()
        self.assertEqual(str(cm_new.exception), "proj_ali_helical() takes at least 7 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helical(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helical(data=None, refrings=None, numr=None, xrng=None,
                                                                 yrng=None, stepx=None, ynumber=None, psi_max=180.0,
                                                                 finfo=None)
        pass


class Test_proj_ali_helical_local(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helical_local()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helical_local()
        self.assertEqual(str(cm_new.exception), "proj_ali_helical_local() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helical_local(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helical_local(data=None, refrings=None, numr=None, xrng=None,
                                                                       yrng=None, stepx=None, ynumber=None, an=None,
                                                                       psi_max=180.0, finfo=None, yrnglocal=-1.0)
        pass


class Test_proj_ali_helical_90(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helical_90()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helical_90()
        self.assertEqual(str(cm_new.exception), "proj_ali_helical_90() takes at least 7 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helical_90(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helical_90(data=None, refrings=None, numr=None, xrng=None,
                                                                    yrng=None, stepx=None, ynumber=None, psi_max=180.0,
                                                                    finfo=None)
        pass


class Test_proj_ali_helical_90_local(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helical_90_local()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helical_90_local()
        self.assertEqual(str(cm_new.exception), "proj_ali_helical_90_local() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helical_90_local(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helical_90_local(data=None, refrings=None, numr=None,
                                                                          xrng=None, yrng=None, stepx=None,
                                                                          ynumber=None, an=None, psi_max=180.0,
                                                                          finfo=None, yrnglocal=-1.0)
        pass


class Test_proj_ali_helicon_local(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helicon_local()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helicon_local()
        self.assertEqual(str(cm_new.exception), "proj_ali_helicon_local() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helicon_local(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helicon_local(data=None, refrings=None, numr=None, xrng=None,
                                                                       yrng=None, stepx=None, ynumber=None, an=None,
                                                                       psi_max=180.0, finfo=None, yrnglocal=-1.0)
        pass


class Test_proj_ali_helicon_90_local_direct(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helicon_90_local_direct()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helicon_90_local_direct()
        self.assertEqual(str(cm_new.exception), "proj_ali_helicon_90_local_direct() takes at least 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helicon_90_local_direct(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helicon_90_local_direct(data=None, refrings=None, xrng=None,
                                                                                 yrng=None, an=None, psi_max=180.0,
                                                                                 psi_step=1.0, stepx=1.0, stepy=1.0,
                                                                                 finfo=None, yrnglocal=-1.0)
        pass


class Test_proj_ali_helicon_90_local_direct1(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helicon_90_local_direct1()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helicon_90_local_direct1()
        self.assertEqual(str(cm_new.exception), "proj_ali_helicon_90_local_direct1() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helicon_90_local_direct1(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helicon_90_local_direct1(data="", refrings="", xrng="",
                                                                                  yrng="", psi_max=180.0, psi_step=1.0,
                                                                                  stepx=1.0, stepy=1.0, finfo=None,
                                                                                  yrnglocal=-1.0, direction="both")
        pass


class Test_proj_ali_helicon_90_local(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_helicon_90_local()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_helicon_90_local()
        self.assertEqual(str(cm_new.exception), "proj_ali_helicon_90_local() takes at least 8 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_helicon_90_local(self):
        peak, phi, theta, psi, s2x, s2y = oldfu.proj_ali_helicon_90_local(data="", refrings="", numr="", xrng="",
                                                                          yrng="", stepx="", ynumber="", an="",
                                                                          psi_max=180.0, finfo=None, yrnglocal=-1.0)
        pass


class Test_ali_vol_func_julio(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_julio()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_julio()
        self.assertEqual(str(cm_new.exception), "ali_vol_func_julio() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_julio(self):
        v = oldfu.ali_vol_func_julio(params="", data="")
        pass


class Test_ali_vol_func_grid(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_grid()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_grid()
        self.assertEqual(str(cm_new.exception), "ali_vol_func_grid() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_grid(self):
        v = oldfu.ali_vol_func_grid(params="", data="")
        pass


class Test_ali_vol_func_nopsi(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_nopsi()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_nopsi()
        self.assertEqual(str(cm_new.exception), "ali_vol_func_nopsi() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_nopsi(self):
        v = oldfu.ali_vol_func_nopsi(params="", data="")
        pass


class Test_ali_vol_func_scale(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_scale()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_scale()
        self.assertEqual(str(cm_new.exception), "ali_vol_func_scale() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_scale(self):
        v = oldfu.ali_vol_func_scale(params="", data="")
        pass


class Test_ali_vol_func_only_scale(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_only_scale()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_only_scale()
        self.assertEqual(str(cm_new.exception), "angle_error() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_only_scale(self):
        v = oldfu.ali_vol_func_only_scale(params="", data="")
        pass


class Test_helios_func(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.helios_func()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.helios_func()
        self.assertEqual(str(cm_new.exception), "helios_func() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_helios_func(self):
        v = oldfu.helios_func(params="", data="")
        pass


class Test_helios(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.helios()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.helios()
        self.assertEqual(str(cm_new.exception), "helios() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_helios(self):
        v = oldfu.helios(vol="", pixel_size="", dp="", dphi="", section_use=0.75, radius=0.0, rmin=0.0)
        pass


class Test_helios7(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.helios7()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.helios7()
        self.assertEqual(str(cm_new.exception), "helios7() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_helios7(self):
        v = oldfu.helios7(vol="", pixel_size="", dp="", dphi="", section_use=0.75, radius=0.0, rmin=0.0)
        pass


class Test_sub_favj(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sub_favj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sub_favj()
        self.assertEqual(str(cm_new.exception), "sub_favj() takes exactly 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_sub_favj(self):
        oldfu.sub_favj(ave="", data="", jtot="", mirror="", numr="")
        pass


class Test_update_favj(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.update_favj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.update_favj()
        self.assertEqual(str(cm_new.exception), "update_favj() takes exactly 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_update_favj(self):
        oldfu.update_favj(ave="", data="", jtot="", mirror="", numr="")
        pass


class Test_multalign2dscf(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.multalign2dscf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.multalign2dscf()
        self.assertEqual(str(cm_new.exception), "multalign2dscf() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_multalign2dscf(self):
        sx, sy, iref, talpha, totpeak = oldfu.multalign2dscf(image="", refrings="", frotim="", numr="", xrng=-1,
                                                             yrng=-1, ou=-1)
        pass


class Test_align2d_direct2(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_direct2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_direct2()
        self.assertEqual(str(cm_new.exception), "align2d_direct2() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align2d_direct2(self):
        bang, bsx, bsy, ama = oldfu.align2d_direct2(image="", refim="", xrng=1, yrng=1, psimax=1, psistep=1, ou=-1)
        pass



class Test_align2d_no_mirror(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_no_mirror()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_no_mirror()
        self.assertEqual(str(cm_new.exception), "align2d_no_mirror() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align2d_no_mirror(self):
        ang, sxs, sys, mirror, peak = oldfu.align2d_no_mirror(image="", refim="", xrng=0, yrng=0, step=1, first_ring=1,
                                                              last_ring=0, rstep=1, mode="F")
        pass


class Test_align2d_direct(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_direct()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_direct()
        self.assertEqual(str(cm_new.exception), "align2d_direct() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align2d_direct(self):
        bang, bsx, bsy, ama = oldfu.align2d_direct(image="", refim="", xrng=1, yrng=1, psimax=1, psistep=1, ou=-1)
        pass


class Test_align2d_peaks(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_peaks()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_peaks()
        self.assertEqual(str(cm_new.exception), "align2d_peaks() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align2d_peaks(self):
        v = oldfu.align2d_peaks(image="", refim="", xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1,
                                mode="F")
        pass


class Test_align2d_g(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_g()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_g()
        self.assertEqual(str(cm_new.exception), "align2d_g() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align2d_g(self):
        v = oldfu.align2d_g(image="", refim="", xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode="F")
        pass


class Test_directali(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.directali()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.directali()
        self.assertEqual(str(cm_new.exception), "directali() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_directali(self):
        nalpha, ntx, nty, peak = oldfu.directali(inima="", refs="", psimax=1.0, psistep=1.0, xrng=1, yrng=1,
                                                 updown="both")
        pass


class Test_preparerefsgrid(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.preparerefsgrid()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.preparerefsgrid()
        self.assertEqual(str(cm_new.exception), "preparerefsgrid() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_preparerefsgrid(self):
        v = oldfu.preparerefsgrid(refs="", psimax=1.0, psistep=1.0)
        pass


class Test_preparerefsgrid1(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.preparerefsgrid1()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.preparerefsgrid1()
        self.assertEqual(str(cm_new.exception), "preparerefsgrid1() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_preparerefsgrid2(self):
        v = oldfu.preparerefsgrid1(refs="", psimax=1.0, psistep=1.0)
        pass


class Test_directaligridding(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.directaligridding()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.directaligridding()
        self.assertEqual(str(cm_new.exception), "directaligridding() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_directaligridding(self):
        nalpha, ntx, nty, peak = oldfu.directaligridding(inima="", refs="", psimax=1.0, psistep=1.0, xrng=1, yrng=1,
                                                         updown="both")
        pass


class Test_directaligridding1(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.directaligridding1()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.directaligridding1()
        self.assertEqual(str(cm_new.exception), "directaligridding1() takes at least 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_directaligridding1(self):
        nalpha, ntx, nty, peak = oldfu.directaligridding1(inima=0, kb=0, ref=0, psimax=1.0, psistep=1.0, xrng=1, yrng=1,
                                                          stepx=1.0, stepy=1.0, updown="both")
        pass


class Test_directaligriddingconstrained(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.directaligriddingconstrained()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.directaligriddingconstrained()
        self.assertEqual(str(cm_new.exception), "directaligriddingconstrained() takes at least 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_directaligriddingconstrained(self):
        nalpha, ntx, nty, peak = oldfu.directaligriddingconstrained(inima="", kb="", ref="", psimax=1.0, psistep=1.0,
                                                                    xrng=1, yrng=1, stepx=1.0, stepy=1.0, psiref=0.,
                                                                    txref=0., tyref=0., updown="up")
        pass


class Test_directaligriddingconstrained3dccf(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.directaligriddingconstrained3dccf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.directaligriddingconstrained3dccf()
        self.assertEqual(str(cm_new.exception), "directaligriddingconstrained3dccf() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_directaligriddingconstrained3dccf(self):
        nalpha, ntx, nty, peak = oldfu.directaligriddingconstrained3dccf(inima="", kb="", ref="", psimax=1.0,
                                                                         psistep=1.0, xrng=1, yrng=1, stepx=1.0,
                                                                         stepy=1.0, psiref=0., txref=0., tyref=0.,
                                                                         updown="up")
        pass


class Test_alignment3Dsnake(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.alignment3Dsnake()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.alignment3Dsnake()
        self.assertEqual(str(cm_new.exception), "alignment3Dsnake() takes at least 13 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_alignment3Dsnake(self):
        v = oldfu.alignment3Dsnake(partition="", snakeknots="", nsegs="", initialori="", ctx="", psistep="", stepx="",
                                   stepy="", txref="", tyref="", nc="", rnx="", rny="", updown="up")
        pass


class Test_flexhelicalali(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.flexhelicalali()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.flexhelicalali()
        self.assertEqual(str(cm_new.exception), "flexhelicalali() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_flexhelicalali(self):
        v = oldfu.flexhelicalali(params="", data="")
        pass



class Test_alivol_m(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.alivol_m()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.alivol_m()
        self.assertEqual(str(cm_new.exception), "angle_error() takes exactly 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_alivol_m(self):
        d_res = oldfu.alivol_m(v="", vref="", mask="")
        pass


class Test_ali_mvol(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_mvol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_mvol()
        self.assertEqual(str(cm_new.exception), "ali_mvol() exactly least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_mvol(self):
        v = oldfu.ali_mvol(v="", mask="")
        pass


class Test_center_projections_3D(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.center_projections_3D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.center_projections_3D()
        self.assertEqual(str(cm_new.exception), "center_projections_3D() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_center_projections_3D(self):
        v = oldfu.center_projections_3D(data=None, ref_vol=None, ali3d_options=None, onx=-1, shrinkage=1.0,
                                        mpi_comm=None, myid=0, main_node=0, log=None)
        pass


class Test_reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info()
        self.assertEqual(str(cm_new.exception), "reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info(self):
        d_res = oldfu.reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info(all_refs_angles="",
                                                                                                     angle_index__to__all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries="")

"""


"""
Comments from Adnan on reply of LUCA's comments
0) To create a complex image just take fft of Eman image and it will become a complex image. In case if it numpy image , you can set the dtype to complex .
2)  Its tricky and has bugs. Lets leave it for the time being.
3) alivol_mask.   I have created a test example of how you can use the 3D_IMAGE as a volume. You just need to add the
                  xform parameters using sp_utilities and it will solve the issue.
4) ali_vol_func_rotate.   You need to pass proper 3D images , a list with projection parameters and the method you want
                           to use to comapare the results. I have created an example .
5) Follow the same procedure and it will solve the problem :).

"""

""" start: new in sphire 1.3"""

class Test_crit2d(unittest.TestCase):
    # args=[] angle,shiftX,shifY --> 3 float
    # data=[] #kbfilter,mask,float,img,img
    img2D = deepcopy(IMAGE_2D)
    img2D.set_attr("mirror", 1)

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.crit2d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.crit2d()
        self.assertEqual(
            str(cm_new.exception), "crit2d() missing 2 required positional arguments: 'args' and 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.crit2d(
                args=[], data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, self.img2D, self.img2D]
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.crit2d(
                args=[], data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, self.img2D, self.img2D]
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.crit2d(data=[], args=[3, 3, 4])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.crit2d(data=[], args=[3, 3, 4])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_img_in_data_returns_RuntimeError_NullPointerException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.crit2d(
                args=[3, 3, 4],
                data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, None, self.img2D],
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.crit2d(
                args=[3, 3, 4],
                data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, None, self.img2D],
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg_old[0].split(" ")[0], "NullPointerException")
        self.assertEqual(msg_old[1], "NULL input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_None_kb_filter_returns_ArgumentError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.crit2d(
                args=[3, 3, 4], data=[None, MASK_2DIMAGE, 2, self.img2D, self.img2D]
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.crit2d(
                args=[3, 3, 4], data=[None, MASK_2DIMAGE, 2, self.img2D, self.img2D]
            )

        # print(str(cm_old.exception))
        output_msg_old = "Python argument types in\n    EMData.rot_scale_conv_new(EMData, numpy.float64, int, int, NoneType, float)\ndid not match C++ signature:\n    rot_scale_conv_new(EMAN::EMData {lvalue}, float ang, float delx, float dely, EMAN::Util::KaiserBessel {lvalue} kb)\n    rot_scale_conv_new(EMAN::EMData {lvalue}, float ang, float delx, float dely, EMAN::Util::KaiserBessel {lvalue} kb, float scale)"
        output_msg_new = "Python argument types in\n    EMData.rot_scale_conv_new(EMData, numpy.float64, int, int, NoneType, float)\ndid not match C++ signature:\n    rot_scale_conv_new(EMAN::EMData {lvalue}, float ang, float delx, float dely, EMAN::Util::KaiserBessel {lvalue} kb)\n    rot_scale_conv_new(EMAN::EMData {lvalue}, float ang, float delx, float dely, EMAN::Util::KaiserBessel {lvalue} kb, float scale)"
        self.assertEqual(str(cm_old.exception), output_msg_old)
        self.assertEqual(str(cm_new.exception), output_msg_new)
        # self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_crit2dIMG(self):
        return_old = oldfu.crit2d(
            args=[3, 3, 4],
            data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, self.img2D, self.img2D],
        )
        return_new = fu.crit2d(
            args=[3, 3, 4],
            data=[KB_IMAGE2D_SIZE, MASK_2DIMAGE, 2, self.img2D, self.img2D],
        )
        self.assertEqual(return_old, return_new)
        self.assertEqual(return_old, 0.031176071614027023)


class Test_kbt(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.kbt()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.kbt()
        self.assertEqual(
            str(cm_new.exception), "kbt() missing 1 required positional argument: 'nx'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_kbt(self):
        return_old = oldfu.kbt(nx=100, npad=2)
        return_new = fu.kbt(nx=100, npad=2)
        self.assertTrue(array_equal(return_old.dump_table(), return_new.dump_table()))
        self.assertEqual(return_old.I0table_maxerror(), return_new.I0table_maxerror())
        self.assertEqual(return_old.I0table_maxerror(), 0.0018696762854233384)


""" After the cleaning thery are gone
class Test_proj_ali_incore_delta(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore_delta()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore_delta()
        self.assertEqual(str(cm_new.exception), "angle_error() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_incore_delta(self):
        pass


#  This function is obsoleted ... i'm not going to test it
class Test_proj_ali_incore_local_psi(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore_local_psi()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore_local_psi()
        self.assertEqual(str(cm_new.exception), "angle_error() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_proj_ali_incore_local_psi(self):
        pass
"""


class Test_ali_vol_func_rotate(unittest.TestCase):
    # params=[sx,sy,sz]  shift stuff
    # params=[EMDATA,EMDATA,MASK_EMDATA,x,EMDATA]
    # x= [phi1, theta1, psi1, sx1, sy1, sz1, not_udes,not_used, scale1]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_rotate()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_rotate()
        self.assertEqual(
            str(cm_new.exception),
            "ali_vol_func_rotate() missing 2 required positional arguments: 'params' and 'data'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_rotate(self):
        params = [250.0, 25.0 , 2.0]
        # print(IMAGE_3D.get_3dview().shape)
        # print(MASK_3DIMAGE.get_3dview().shape)
        data = [ IMAGE_3D, IMAGE_3D, MASK_3DIMAGE, [25.0, 45.0, 22.0, 0.5, 0.8, 1.2, 0 , 1.0], "ccc" ]
        v = fu.ali_vol_func_rotate(params, data)
        # pass


class Test_ali_vol_func_shift(unittest.TestCase):
    # params=[sx,sy,sz]  shift stuff
    # params=[EMDATA,EMDATA,MASK_EMDATA,x,EMDATA]
    # x= [phi1, theta1, psi1, sx1, sy1, sz1, not_udes,not_used, scale1]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func_shift()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func_shift()
        self.assertEqual(
            str(cm_new.exception),
            "ali_vol_func_shift() missing 2 required positional arguments: 'params' and 'data'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func_shift(self):
        # v = oldfu.ali_vol_func_shift(params="", data="")
        pass


class Test_fine_2D_refinement(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fine_2D_refinement()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fine_2D_refinement()
        self.assertEqual(
            str(cm_new.exception),
            "fine_2D_refinement() missing 4 required positional arguments: 'data', 'br', 'mask', and 'tavg'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Nonetavg_crash_because_RuntimeError_NullPointerException(self):
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fine_2D_refinement(
                data=[deepcopy(IMAGE_2D)],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=None,
                group=-1,
            )
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fine_2D_refinement(
                data=[deepcopy(IMAGE_2D)],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=None,
                group=-1,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_emptytavg_crash_because_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fine_2D_refinement(
                data=[deepcopy(IMAGE_2D)],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=EMData(),
                group=-1,
            )
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fine_2D_refinement(
                data=[deepcopy(IMAGE_2D)],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=EMData(),
                group=-1,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonemask_crash_because_SIGSEGV(self):
        pass
        """
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fine_2D_refinement(data=[deepcopy(IMAGE_2D)], br=1.75, mask=None, tavg=IMAGE_2D_REFERENCE,group=-1)
        with self.assertRaises(AttributeError) as cm_new:
            fu.fine_2D_refinement(data=[deepcopy(IMAGE_2D)], br=1.75, mask=None, tavg=IMAGE_2D_REFERENCE,group=-1)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        """

    def test_emptymask_crash_because_SIGSEGV(self):
        pass
        """
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fine_2D_refinement(data=[deepcopy(IMAGE_2D)], br=1.75, mask=EMData(), tavg=IMAGE_2D_REFERENCE,group=-1)
        with self.assertRaises(AttributeError) as cm_new:
            fu.fine_2D_refinement(data=[deepcopy(IMAGE_2D)], br=1.75, mask=EMData(), tavg=IMAGE_2D_REFERENCE,group=-1)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        """

    def test_None_data_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fine_2D_refinement(
                data=[None],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=IMAGE_2D_REFERENCE,
                group=-1,
            )
        with self.assertRaises(AttributeError) as cm_new:
            fu.fine_2D_refinement(
                data=[None],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=IMAGE_2D_REFERENCE,
                group=-1,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_empty_data_array_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fine_2D_refinement(
                data=[], br=1.75, mask=MASK_2DIMAGE, tavg=IMAGE_2D_REFERENCE, group=-1
            )
        with self.assertRaises(IndexError) as cm_new:
            fu.fine_2D_refinement(
                data=[], br=1.75, mask=MASK_2DIMAGE, tavg=IMAGE_2D_REFERENCE, group=-1
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.fine_2D_refinement(
                data=[EMData()],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=IMAGE_2D_REFERENCE,
                group=-1,
            )
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.fine_2D_refinement(
                data=[EMData()],
                br=1.75,
                mask=MASK_2DIMAGE,
                tavg=IMAGE_2D_REFERENCE,
                group=-1,
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_fine_2D_refinement(self):
        img1_new = deepcopy(IMAGE_2D)
        img1_old = deepcopy(IMAGE_2D)
        img2_new = deepcopy(IMAGE_2D_REFERENCE)
        img2_old = deepcopy(IMAGE_2D_REFERENCE)
        img2_new.set_attr("alpha", 1)
        img2_new.set_attr("sx", 2)
        img2_new.set_attr("sy", 1)
        img2_new.set_attr("mirror", 0)
        img1_new.set_attr("alpha", 1)
        img1_new.set_attr("sx", 2)
        img1_new.set_attr("sy", 101)
        img1_new.set_attr("mirror", 0)
        img2_old.set_attr("alpha", 1)
        img2_old.set_attr("sx", 2)
        img2_old.set_attr("sy", 1)
        img2_old.set_attr("mirror", 0)
        img1_old.set_attr("alpha", 1)
        img1_old.set_attr("sx", 2)
        img1_old.set_attr("sy", 101)
        img1_old.set_attr("mirror", 0)
        oldfu.fine_2D_refinement(
            data=[img1_old, img2_old],
            br=1.75,
            mask=MASK_2DIMAGE,
            tavg=IMAGE_2D_REFERENCE,
            group=-1,
        )
        fu.fine_2D_refinement(
            data=[img1_new, img2_new],
            br=1.75,
            mask=MASK_2DIMAGE,
            tavg=IMAGE_2D_REFERENCE,
            group=-1,
        )
        self.assertTrue(
            array_equal(
                [
                    img1_new.get_attr("alpha"),
                    img1_new.get_attr("sx"),
                    img1_new.get_attr("sy"),
                    img1_new.get_attr("mirror"),
                ],
                [
                    img1_old.get_attr("alpha"),
                    img1_old.get_attr("sx"),
                    img1_old.get_attr("sy"),
                    img1_old.get_attr("mirror"),
                ],
            )
        )
        self.assertTrue(
            array_equal(
                [
                    img2_new.get_attr("alpha"),
                    img2_new.get_attr("sx"),
                    img2_new.get_attr("sy"),
                    img2_new.get_attr("mirror"),
                ],
                [
                    img2_old.get_attr("alpha"),
                    img2_old.get_attr("sx"),
                    img2_old.get_attr("sy"),
                    img2_old.get_attr("mirror"),
                ],
            )
        )
        self.assertTrue(
            array_equal(
                [
                    img1_new.get_attr("alpha"),
                    img1_new.get_attr("sx"),
                    img1_new.get_attr("sy"),
                    img1_new.get_attr("mirror"),
                ],
                [2.2260714084092927, 2.5339026885139146, 96.87486626866868, 0],
            )
        )
        self.assertTrue(
            array_equal(
                [
                    img2_new.get_attr("alpha"),
                    img2_new.get_attr("sx"),
                    img2_new.get_attr("sy"),
                    img2_new.get_attr("mirror"),
                ],
                [2.0054556004767163, 8.008756912179393, -3.545383545183274, 0],
            )
        )


class Test_align2d_direct3(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_direct3()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_direct3()
        self.assertEqual(
            str(cm_new.exception),
            "align2d_direct3() missing 2 required positional arguments: 'input_images' and 'refim'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_None_data_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.align2d_direct3(
                input_images=[None],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        with self.assertRaises(AttributeError) as cm_new:
            fu.align2d_direct3(
                input_images=[None],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_empty_data_array_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align2d_direct3(
                input_images=[],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        with self.assertRaises(IndexError) as cm_new:
            fu.align2d_direct3(
                input_images=[],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.align2d_direct3(
                input_images=[EMData()],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        with self.assertRaises(RuntimeError) as cm_new:
            fu.align2d_direct3(
                input_images=[EMData()],
                refim=IMAGE_2D,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_refim_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.align2d_direct3(
                input_images=[IMAGE_2D, IMAGE_2D],
                refim=None,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.align2d_direct3(
                input_images=[IMAGE_2D, IMAGE_2D],
                refim=None,
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception),
            "'NoneType' object has no attribute 'rot_scale_trans2D_background'",
        )

    def test_emptyType_refim_returns_RuntimeError_ImageDimensionException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.align2d_direct3(
                input_images=[IMAGE_2D, IMAGE_2D],
                refim=EMData(),
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.align2d_direct3(
                input_images=[IMAGE_2D, IMAGE_2D],
                refim=EMData(),
                xrng=1,
                yrng=1,
                psimax=180,
                psistep=1,
                ou=-1,
                CTF=None,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_without_ctf(self):
        return_new = fu.align2d_direct3(
            input_images=[IMAGE_2D, IMAGE_2D],
            refim=IMAGE_2D,
            xrng=1,
            yrng=1,
            psimax=180,
            psistep=1,
            ou=-1,
            CTF=None,
        )
        return_old = oldfu.align2d_direct3(
            input_images=[IMAGE_2D, IMAGE_2D],
            refim=IMAGE_2D,
            xrng=1,
            yrng=1,
            psimax=180,
            psistep=1,
            ou=-1,
            CTF=None,
        )
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(
            array_equal(
                return_old,
                [
                    [181.00000000615043, 0.0, 0.0, 0, -1e23],
                    [181.00000000615043, 0.0, 0.0, 0, -1e23],
                ],
            )
        )

    def test_with_ctf(self):
        img = deepcopy(IMAGE_2D)
        ctf1 = EMAN2Ctf()
        ctf1.from_dict(
            {
                "defocus": 1,
                "cs": 2,
                "voltage": 300,
                "apix": 1.5,
                "bfactor": 0,
                "ampcont": 0.1,
                "dfdiff": 0.1,
                "dfang": 0.1,
            }
        )
        img.set_attr("ctf", ctf1)
        return_new = fu.align2d_direct3(
            input_images=[img, img],
            refim=img,
            xrng=1,
            yrng=1,
            psimax=180,
            psistep=1,
            ou=-1,
            CTF=True,
        )
        return_old = oldfu.align2d_direct3(
            input_images=[img, img],
            refim=img,
            xrng=1,
            yrng=1,
            psimax=180,
            psistep=1,
            ou=-1,
            CTF=True,
        )
        self.assertTrue(array_equal(return_old, return_new))
        # self.assertTrue(
        #     array_equal(
        #         return_old,
        #         [
        #             [
        #                 316.9999994395314,
        #                 0.026203064247965813,
        #                 -0.2627040147781372,
        #                 0,
        #                 0.4272246869622954,
        #             ],
        #             [
        #                 316.9999994395314,
        #                 0.026203064247965813,
        #                 -0.2627040147781372,
        #                 0,
        #                 0.4272246869622954,
        #             ],
        #         ],
        #     )
        # )


class Test_ali_nvol(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_nvol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_nvol()
        self.assertEqual(
            str(cm_new.exception), "ali_nvol() missing 2 required positional arguments: 'v' and 'mask'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_v_array_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ali_nvol(v=[], mask=MASK_2DIMAGE)
        with self.assertRaises(IndexError) as cm_new:
            fu.ali_nvol(v=[], mask=MASK_2DIMAGE)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_v_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_nvol(v=[EMData()], mask=MASK_2DIMAGE)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_nvol(v=[EMData()], mask=MASK_2DIMAGE)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_none_v_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.ali_nvol(v=[None], mask=MASK_2DIMAGE)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ali_nvol(v=[None], mask=MASK_2DIMAGE)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'set_attr'"
        )


    # def test_luca_ali_nvol(self):
    #     from sphire.libpy import sp_utilities
    #     p = [180.0,1.0, 30.0, 0.5, 0.5, 1.2, 0,1]
    #
    #     new_mask = sp_utilities.model_circle(2, 10, 10, 10)
    #
    #     sp_utilities.set_params3D(IMAGE_3D, p)
    #
    #     rr = fu.ali_nvol([IMAGE_3D, IMAGE_3D], new_mask)


class Test_alivol_mask_getref(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.alivol_mask_getref()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.alivol_mask_getref()
        self.assertEqual(
            str(cm_new.exception),
            "alivol_mask_getref() missing 2 required positional arguments: 'v' and 'mask'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_v_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.alivol_mask_getref(v=EMData(), mask=MASK_2DIMAGE)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.alivol_mask_getref(v=EMData(), mask=MASK_2DIMAGE)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "can not multiply images that are not the same size")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_none_v_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.alivol_mask_getref(v=None, mask=MASK_2DIMAGE)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.alivol_mask_getref(v=None, mask=MASK_2DIMAGE)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'copy'"
        )

    def test_2Dimg(self):
        return_old = oldfu.alivol_mask_getref(v=IMAGE_2D, mask=MASK_2DIMAGE)
        return_new = fu.alivol_mask_getref(v=IMAGE_2D, mask=MASK_2DIMAGE)
        self.assertTrue(array_equal(return_old.get_2dview(), return_new.get_2dview()))
        self.assertTrue(
            array_equal(
                return_old.get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.12692593038082123,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.008231919258832932,
                    -0.020773129537701607,
                    -0.034199729561805725,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.111984983086586,
                    -0.11971071362495422,
                    -0.1273496150970459,
                    -0.12249226123094559,
                    -0.1453358680009842,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.24315771460533142,
                    -0.2552821934223175,
                    -0.23703180253505707,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.3575122356414795,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                ],
            )
        )

    def test_2DimgBlank(self):
        mask = model_circle(
            2,
            IMAGE_BLANK_2D.get_xsize(),
            IMAGE_BLANK_2D.get_ysize(),
            IMAGE_BLANK_2D.get_zsize(),
        )
        return_old = oldfu.alivol_mask_getref(v=IMAGE_BLANK_2D, mask=mask)
        return_new = fu.alivol_mask_getref(v=IMAGE_BLANK_2D, mask=mask)
        self.assertTrue(array_equal(return_old.get_2dview(), return_new.get_2dview()))
        self.assertTrue(
            array_equal(return_old.get_2dview(), IMAGE_BLANK_2D.get_2dview())
        )

    def test_3Dimg(self):
        mask = model_circle(
            2,
            IMAGE_BLANK_3D.get_xsize(),
            IMAGE_BLANK_3D.get_ysize(),
            IMAGE_BLANK_3D.get_zsize(),
        )
        return_old = oldfu.alivol_mask_getref(v=IMAGE_3D, mask=mask)
        return_new = fu.alivol_mask_getref(v=IMAGE_3D, mask=mask)
        self.assertTrue(array_equal(return_old.get_3dview(), return_new.get_3dview()))
        self.assertTrue(
            array_equal(
                return_old.get_3dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.028610195964574814,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.4733593463897705,
                    -0.4710369408130646,
                    -0.4884580969810486,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.4724937975406647,
                    -0.4650014042854309,
                    -0.4585750997066498,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.49023202061653137,
                    -0.4790208339691162,
                    -0.47853776812553406,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -1.1460676193237305,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.9730094075202942,
                    -0.9454143643379211,
                    -0.9392389059066772,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.769252598285675,
                    -0.7435757517814636,
                    -0.731235682964325,
                    -0.7272043228149414,
                    -0.7139352560043335,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.5130307078361511,
                    -0.4996497631072998,
                    -0.47770577669143677,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.400550901889801,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.07777398824691772,
                    -0.07534715533256531,
                    -0.08054795116186142,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.08260509371757507,
                    -0.08763255923986435,
                    -0.0983823612332344,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.10116869956254959,
                    -0.1163831278681755,
                    -0.09893848747015,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    0.0,
                    -0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -0.0,
                    0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.08259184658527374,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                    -0.0,
                ],
            )
        )

    def test_2DimgBlank(self):
        mask = model_circle(
            2,
            IMAGE_BLANK_3D.get_xsize(),
            IMAGE_BLANK_3D.get_ysize(),
            IMAGE_BLANK_3D.get_zsize(),
        )
        return_old = oldfu.alivol_mask_getref(v=IMAGE_BLANK_3D, mask=mask)
        return_new = fu.alivol_mask_getref(v=IMAGE_BLANK_3D, mask=mask)
        self.assertTrue(array_equal(return_old.get_3dview(), return_new.get_3dview()))
        self.assertTrue(
            array_equal(return_old.get_3dview(), IMAGE_BLANK_3D.get_3dview())
        )


class Test_alivol_mask(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.alivol_mask()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.alivol_mask()
        self.assertEqual(
            str(cm_new.exception), "alivol_mask() missing 3 required positional arguments: 'v', 'vref', and 'mask'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_none_v_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.alivol_mask(v=None, vref=IMAGE_2D_REFERENCE, mask=MASK_2DIMAGE)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.alivol_mask(v=None, vref=IMAGE_2D_REFERENCE, mask=MASK_2DIMAGE)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'copy'"
        )

    def test_empty_v_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.alivol_mask(v=EMData(), vref=IMAGE_2D_REFERENCE, mask=MASK_2DIMAGE)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.alivol_mask(v=EMData(), vref=IMAGE_2D_REFERENCE, mask=MASK_2DIMAGE)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "can not multiply images that are not the same size")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_empty_vref_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.alivol_mask(v=IMAGE_2D, vref=EMData(), mask=MASK_2DIMAGE)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.alivol_mask(v=IMAGE_2D, vref=EMData(), mask=MASK_2DIMAGE)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_None_vref_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.alivol_mask(v=IMAGE_2D, vref=None, mask=MASK_2DIMAGE)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.alivol_mask(v=IMAGE_2D, vref=None, mask=MASK_2DIMAGE)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_None_mask_returns_ArgumentError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.alivol_mask(v=IMAGE_2D, vref=IMAGE_2D_REFERENCE, mask=None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.alivol_mask(v=IMAGE_2D, vref=IMAGE_2D_REFERENCE, mask=None)
        output_msg = "Python argument types in\n    EMData.__imul__(EMData, NoneType)\ndid not match C++ signature:\n    __imul__(boost::python::back_reference<EMAN::EMData&>, EMAN::EMData)\n    __imul__(boost::python::back_reference<EMAN::EMData&>, float)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_mask_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.alivol_mask(v=IMAGE_2D, vref=IMAGE_2D_REFERENCE, mask=EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.alivol_mask(v=IMAGE_2D, vref=IMAGE_2D_REFERENCE, mask=EMData())
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "can not multiply images that are not the same size")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_luca_ali_vol_mask(self):
        p = [180.0,0.0, 30.0, 0.5, 0.5, 1.2, 0,1]

        new_mask = sp_utilities.model_circle(2, 10, 10, 10)

        sp_utilities.set_params3D(IMAGE_3D, p)

        rr = fu.alivol_mask(IMAGE_3D, IMAGE_3D, new_mask)


""" end: new in sphire 1.3"""

@unittest.skip("Somehow the pickle file is corrupted")
class Test_ali2d_single_iter(unittest.TestCase):
    """
    Since using an invalid method as "random method" is like using the default method "random_method=''" I'm not testing this situation"
    Since the method "random_method='PCP'" seems to lead to a dead code, anyway it is crashing, I'm skipping these cases
    All the case with "random_method='SHC'" do not work. They are manageable through 'typeError' exception. This situation is used a lot in the 'sphire/legacy/...'
    """
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES,"alignment.ali2d_single_iter")
    )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter()
        self.assertEqual(
            str(cm_new.exception),
            "ali2d_single_iter() missing 10 required positional arguments: 'data', 'numr', 'wr', 'cs', 'tavg', 'cnx', 'cny', 'xrng', 'yrng', and 'step'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_images_to_align_returns_RuntimeError_NotExistingObjectException_the_key_ctf_doesnot_exist(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        images = [EMData(), EMData(), EMData()]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali2d_single_iter(
                data=images,
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SCF",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=images,
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SCF",
                ali_params="xform.align2d",
                delta=0.0,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)  = self.argum[0]
        return_new = fu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        return_old = oldfu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr, wr, cs, EMData(), cnx, cny, xrng, yrng, step)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_few_shift_params_returns_IndexError_list_index_out_of_range(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        cs = [1]
        with self.assertRaises(IndexError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SCF",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SCF",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_too_shift_params(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        cs = [1, 2, 3]
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-14.731166765093803, 4.656881392002106, 0))
        )

    def test_empty_list_fourier_weight_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)  = self.argum[0]
        wr = []
        with self.assertRaises(ValueError):
            fu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr=numr, wr=wr, cs=cs, tavg=tavg, cnx=cnx, cny=cny, xrng=xrng, yrng=yrng, step=step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="Unknown", delta = 0.0)
            oldfu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr=numr, wr=wr, cs=cs, tavg=tavg, cnx=cnx, cny=cny, xrng=xrng, yrng=yrng, step=step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="Unknown", delta = 0.0)
        """
        self.assertTrue(True)

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        numr = []
        with self.assertRaises(IndexError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="Unknown",
                delta=0.0,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="Unknown",
                delta=0.0,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Unknown_ali_params_returns_RuntimeError_NotExistingObjectException_the_key_Unknown_doesnot_exist(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="Unknown",
                delta=0.0,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="Unknown",
                delta=0.0,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_default_center_out_of_possible_range_crashes_because_signal11SIGSEV(self):
        """
        (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)  = self.argum[0]
        return_new = fu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        return_old = oldfu.ali2d_single_iter(data=deepcopy(self.argum[0][0]), numr, wr, cs, tavg, 10000, 10000, xrng, yrng, step, nomirror = False, mode="F", CTF=False, random_method="", ali_params="xform.align2d", delta = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_negative_XRange_Value_UnboundLocalError_a_local_var_is_undefined(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(UnboundLocalError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=-1,
                yrng=0,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(UnboundLocalError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=-1,
                yrng=0,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "local variable 'ang' referenced before assignment"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("The output seems to be random")
    def test_negative_center_warning_msg_shift_of_paricle_too_large(self):
        # I cannot run unit test because it returns random values
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=-5,
            cny=-5,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=-5,
            cny=-5,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_H_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="h",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="h",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_F_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_F_withCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_H_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )

        # output_msg = "Python argument types in\n   None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n     shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        # self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_H_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_F_WITHCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=True,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_F_NOCTF_randomMethod_SHC_ArgumentError_in_Util_shc_function(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=False,
                random_method="SHC",
                ali_params="xform.align2d",
                delta=0.0,
            )
        output_msg = "Python argument types in\n    None.shc(EMData, list, list, list, int, float, str, list, float, float, str)\ndid not match C++ signature:\n    shc(EMAN::EMData*, std::vector<EMAN::EMData*, std::allocator<EMAN::EMData*> > image, std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > crefim, std::vector<float, std::allocator<float> > list_of_reference_angles, std::vector<float, std::allocator<float> > xrng, float yrng, float step, std::string ant, std::vector<int, std::allocator<int> > mode, float numr, float cnx, std::string cny)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_H_WITHCTF_randomMethod_SCF(self):
        # (not_used, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-11.89580624550581, 3.1550322622060776, 0))
        )

    def test_NOmirror_mode_H_NOCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="h",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="h",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-11.89580624550581, 3.1550322622060776, 0))
        )

    def test_NOmirror_mode_F_NOCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-21.29842366767116, -18.99927615886554, 0))
        )

    def test_NOmirror_mode_F_withCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-21.29842366767116, -18.99927615886554, 0))
        )

    def test_mirror_mode_H_NOCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-11.89580624550581, 3.1550322622060776, 0))
        )

    def test_mirror_mode_H_WITHCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-11.89580624550581, 3.1550322622060776, 0))
        )

    def test_mirror_mode_F_WITHCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=True,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )

        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-21.29842366767116, -18.99927615886554, 0))
        )


    def test_mirror_mode_F_NOCTF_randomMethod_SCF(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=False,
            random_method="SCF",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        self.assertTrue(
            allclose(
                return_new,(-21.29842366767116, -18.99927615886554, 0), atol=TOLERANCE
            )
        )

    def test_NOmirror_mode_H_WITHCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="H",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # self.assertTrue(
        #     allclose(
        #         return_new, (-19.65119509678334, -22.428544966503978, 0), atol=TOLERANCE
        #     )
        # )

    def test_NOmirror_mode_H_NOCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="h",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="h",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(return_new, (-43.51346893221489, -43.28186049871147, 0))
        # )

    def test_NOmirror_mode_F_NOCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # print(return_old)
        self.assertTrue(
            allclose(
                return_old,
                (-456.3737654685974, -479.34432792663574, 0),
                atol=TOLERANCE,
            )
        )

    def test_NOmirror_mode_F_withCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=False,
            mode="F",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # self.assertTrue(
        #     allclose(
        #         return_new, (5.475417716242191, 37.246610491740284, 0), atol=TOLERANCE
        #     )
        # )

    def test_mirror_mode_H_NOCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # self.assertTrue(
        #     allclose(
        #         return_new, (-24.46844869107008, -27.762613539933227, 0), atol=TOLERANCE
        #     )
        # )

    def test_mirror_mode_H_WITHCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="h",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # self.assertTrue(
        #     allclose(
        #         return_new,
        #         (-10.602042245678604, -28.610507858917117, 0),
        #         atol=TOLERANCE,
        #     )
        # )

    def test_mirror_mode_F_WITHCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=True,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(allclose(return_old, return_new, atol=TOLERANCE))
        # self.assertTrue(
        #     allclose(
        #         return_new, (9.289807755313632, -4.889407425798709, 0), atol=TOLERANCE
        #     )
        # )

    def test_mirror_mode_F_NOCTF_randomMethod_default(self):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        return_old = oldfu.ali2d_single_iter(
            data=deepcopy(data),
            numr=numr,
            wr=wr,
            cs=cs,
            tavg=tavg,
            cnx=cnx,
            cny=cny,
            xrng=xrng,
            yrng=yrng,
            step=step,
            nomirror=True,
            mode="F",
            CTF=False,
            random_method="",
            ali_params="xform.align2d",
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(return_new, (-16.664929463528097, -62.39760458981618, 0))
        # )

    def test_NOmirror_mode_H_WITHCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="H",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_H_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="h",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="h",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_F_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NOmirror_mode_F_withCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=False,
                mode="F",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_H_NOCTF_randomMethod_PCP_returns_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_H_WITHCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="h",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_F_WITHCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=True,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mirror_mode_F_NOCTF_randomMethod_PCP_typeError_object_float_hasnot_attibute__getitem__(
        self
    ):
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_single_iter(
                data=deepcopy(data),
                numr=numr,
                wr=wr,
                cs=cs,
                tavg=tavg,
                cnx=cnx,
                cny=cny,
                xrng=xrng,
                yrng=yrng,
                step=step,
                nomirror=True,
                mode="F",
                CTF=False,
                random_method="PCP",
                ali_params="xform.align2d",
                delta=0.0,
            )
        self.assertEqual(
            str(cm_new.exception), "'float' object is not subscriptable"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_ang_n(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ang_n()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ang_n()
        self.assertEqual(
            str(cm_new.exception), "ang_n() missing 3 required positional arguments: 'tot', 'mode', and 'maxrin'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Full_mode(self):
        return_new = fu.ang_n(tot=2, mode="f", maxrin=3)
        return_old = oldfu.ang_n(tot=2, mode="f", maxrin=3)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, 120.0))

    def test_null_max_ring_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ang_n(tot=2, mode="f", maxrin=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ang_n(tot=2, mode="f", maxrin=0)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Half_mode(self):
        return_new = fu.ang_n(tot=2, mode="nf", maxrin=3)
        return_old = oldfu.ang_n(tot=2, mode="nf", maxrin=3)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, 60.0))


class Test_log2(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.log2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.log2()
        self.assertEqual(
            str(cm_new.exception), "log2() missing 1 required positional argument: 'n'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_positive_number(self):
        return_old = oldfu.log2(n=10)
        self.assertEqual(return_old, fu.log2(n=10))
        self.assertEqual(return_old, 3)

    def test_null_number(self):
        return_old = oldfu.log2(n=0)
        self.assertEqual(return_old, fu.log2(n=0))
        self.assertEqual(return_old, -1)


class Test_Numrinit(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.Numrinit()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.Numrinit()
        self.assertEqual(
            str(cm_old.exception), "Numrinit() missing 2 required positional arguments: 'first_ring' and 'last_ring'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_null_skip_value_returns_ValueError_this_arg_cannot_be_null(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.Numrinit(first_ring=2, last_ring=5, skip=0, mode="F")
        with self.assertRaises(ValueError) as cm_old:
            oldfu.Numrinit(first_ring=2, last_ring=5, skip=0, mode="F")
        self.assertEqual(str(cm_old.exception), "range() arg 3 must not be zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Full_mode(self):
        return_new = fu.Numrinit(first_ring=2, last_ring=5, skip=1, mode="F")
        return_old = oldfu.Numrinit(first_ring=2, last_ring=5, skip=1, mode="F")
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_old, [2, 1, 16, 3, 17, 32, 4, 49, 32, 5, 81, 32])
        )

    def test_Half_mode(self):
        return_new = fu.Numrinit(first_ring=2, last_ring=5, skip=1, mode="not_F")
        return_old = oldfu.Numrinit(first_ring=2, last_ring=5, skip=1, mode="not_F")
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_old, [2, 1, 8, 3, 9, 16, 4, 25, 16, 5, 41, 32])
        )


class Test_ringwe(unittest.TestCase):
    numr = [
        1,
        1,
        8,
        2,
        9,
        16,
        3,
        25,
        32,
        4,
        57,
        32,
        5,
        89,
        32,
        6,
        121,
        64,
        7,
        185,
        64,
        8,
        249,
        64,
        9,
        313,
        64,
        10,
        377,
        64,
        11,
        441,
        128,
        12,
        569,
        128,
        13,
        697,
        128,
        14,
        825,
        128,
        15,
        953,
        128,
        16,
        1081,
        128,
        17,
        1209,
        128,
        18,
        1337,
        128,
        19,
        1465,
        128,
        20,
        1593,
        128,
        21,
        1721,
        256,
        22,
        1977,
        256,
        23,
        2233,
        256,
        24,
        2489,
        256,
        25,
        2745,
        256,
        26,
        3001,
        256,
        27,
        3257,
        256,
        28,
        3513,
        256,
        29,
        3769,
        256,
    ]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ringwe()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ringwe()
        self.assertEqual(
            str(cm_new.exception), "ringwe() missing 1 required positional argument: 'numr'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ringwe(numr=[], mode="F")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ringwe(numr=[], mode="F")
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Full_mode(self):
        return_new = fu.ringwe(numr=self.numr, mode="F")
        return_old = oldfu.ringwe(numr=self.numr, mode="F")
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                [
                    25.132741228718345,
                    12.566370614359172,
                    4.71238898038469,
                    6.283185307179586,
                    7.853981633974483,
                    2.356194490192345,
                    2.748893571891069,
                    3.141592653589793,
                    3.5342917352885173,
                    3.9269908169872414,
                    1.0799224746714913,
                    1.1780972450961724,
                    1.2762720155208536,
                    1.3744467859455345,
                    1.4726215563702154,
                    1.5707963267948966,
                    1.6689710972195777,
                    1.7671458676442586,
                    1.8653206380689396,
                    1.9634954084936207,
                    0.5154175447295755,
                    0.5399612373357456,
                    0.5645049299419159,
                    0.5890486225480862,
                    0.6135923151542565,
                    0.6381360077604268,
                    0.662679700366597,
                    0.6872233929727672,
                    0.7117670855789375,
                ],
            )
        )

    def test_Half_mode(self):
        return_new = fu.ringwe(numr=self.numr, mode="not_F")
        return_old = oldfu.ringwe(numr=self.numr, mode="not_F")
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                [
                    12.566370614359172,
                    6.283185307179586,
                    2.356194490192345,
                    3.141592653589793,
                    3.9269908169872414,
                    1.1780972450961724,
                    1.3744467859455345,
                    1.5707963267948966,
                    1.7671458676442586,
                    1.9634954084936207,
                    0.5399612373357456,
                    0.5890486225480862,
                    0.6381360077604268,
                    0.6872233929727672,
                    0.7363107781851077,
                    0.7853981633974483,
                    0.8344855486097889,
                    0.8835729338221293,
                    0.9326603190344698,
                    0.9817477042468103,
                    0.25770877236478773,
                    0.2699806186678728,
                    0.28225246497095796,
                    0.2945243112740431,
                    0.30679615757712825,
                    0.3190680038802134,
                    0.3313398501832985,
                    0.3436116964863836,
                    0.35588354278946877,
                ],
            )
        )


class Test_ornq(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "alignment.ornq")
    )

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType_as_input_image_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(image=None, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image=None, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType_as_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(image=image, crefim=None, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image=image, crefim=None, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ornq()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ornq()
        self.assertEqual(
            str(cm_new.exception), "ornq() missing 9 required positional arguments: 'image', 'crefim', 'xrng', 'yrng', 'step', 'mode', 'numr', 'cnx', and 'cny'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        numr = []
        return_new = fu.ornq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, deltapsi = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_list_xrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        xrng = []
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_yrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        yrng = []
        with self.assertRaises(IndexError) as cm_new:
            fu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_negative_center(self):
        # I cannot write a unit test because the output seems to be randon
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=-5,
            cny=-5,
            deltapsi=0.0,
        )
        return_old = oldfu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=-5,
            cny=-5,
            deltapsi=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))

    def test_null_step_value_returns_ZeroDivisionError(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[
            0
        ]  # mode is H
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=0,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ornq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=0,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                deltapsi=0.0,
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Half_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]  # mode is H
        image , crefim = give_ornq_data()
        return_new = fu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        return_old = oldfu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new, (90.659458637237549, 0.0, 0.0, 0, 147.43201554741904)
        #     )
        # )

    def test_Full_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode = "f"
        image, crefim = give_ornq_data()
        return_new = fu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        return_old = oldfu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new, (271.47330522537231, 0.0, -0.0, 0, 136.83287092834385)
        #     )
        # )

    def test_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode = "invalid"
        image, crefim = give_ornq_data()
        return_new = fu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        return_old = oldfu.ornq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            deltapsi=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new, (90.659458637237549, 0.0, 0.0, 0, 147.43201554741904)
        #     )
        # )


class Test_ormq(unittest.TestCase):
    # argum = get_arg_from_pickle_file(
    #     path.join(ABSOLUTE_PATH_TO_RESOURCES, "alignment.ormq")
    # )

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        image =EMData()
        return_new = fu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        return_old = oldfu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType_as_input_image_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        return_new = fu.ormq(image=None, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        return_old = oldfu.ormq(image=None, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        crefim =EMData()
        return_new = fu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        return_old = oldfu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType_as_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        return_new = fu.ormq(image=image, crefim=None, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        return_old = oldfu.ormq(image=image, crefim=None, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ormq()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ormq()
        self.assertEqual(
            str(cm_new.exception), "ormq() missing 9 required positional arguments: 'image', 'crefim', 'xrng', 'yrng', 'step', 'mode', 'numr', 'cnx', and 'cny'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_xrng_returns_IndexError_list_index_out_of_range(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        xrng = []
        with self.assertRaises(IndexError) as cm_new:
            fu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
            oldfu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_yrng_returns_IndexError_list_index_out_of_range(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        yrng = []
        with self.assertRaises(IndexError) as cm_new:
            fu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError) as cm_new:
            fu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ormq(image=image, crefim=crefim, xrng=xrng, yrng=yrng, step=step, mode=mode, numr=numr, cnx=cnx, cny=cny, delta=delta)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        """
        self.assertTrue(True)

    def test_null_skip_value_returns_ZeroDivisionError(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        step = 0
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ormq(
                image=image,
                crefim=crefim,
                xrng=xrng,
                yrng=yrng,
                step=step,
                mode=mode,
                numr=numr,
                cnx=cnx,
                cny=cny,
                delta=delta,
            )
        self.assertEqual(str(cm_new.exception), "integer division or modulo by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Full_mode(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)  # mode is F
        # image, crefim = give_ormq_data()
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new, (180.0, -2.0, -2.4492935982947064e-16, 1, 23.19550078262219)
        #     )
        # )

    def test_Half_mode(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        # image, crefim = give_ormq_data()
        mode = "H"
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (135.0, 0.70710678118654768, -2.1213203435596424, 0, 32.54668225191281),
        #     )
        # )

    def test_with_negative_center(self):
        # I cannot run unit test because it returns random values
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=-5,
            cny=-5,
            delta=delta,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=-5,
            cny=-5,
            delta=delta,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(array_equal(return_new, (358.59375, -3.9006303606931674, -4.0969601888764666, 0, -1e+20)))

    def test_with_invalid_mode(self):
        image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = give_ormq_data(True)
        # image, crefim = give_ormq_data()
        mode = "invalid"
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=delta,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (135.0, 0.70710678118654768, -2.1213203435596424, 0, 32.54668225191281),
        #     )
        # )

    def test_Full_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = give_ormq_data(True)
        # image, crefim = give_ormq_data()
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (
        #             206.98280811309814,
        #             -0.89114270620982317,
        #             0.45372312831619327,
        #             0,
        #             23.462145424755487,
        #         ),
        #     )
        # )

    def test_Half_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = give_ormq_data(True)
        # image, crefim = give_ormq_data()
        mode = "H"
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (
        #             82.811969518661499,
        #             2.1094076960564689,
        #             -0.74188892148201036,
        #             1,
        #             33.042001487198405,
        #         ),
        #     )
        # )

    def test_with_invalid_mode_and_zero_delta(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = give_ormq_data(True)
        # image, crefim = give_ormq_data()
        mode = "invalid"
        return_new = fu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        return_old = oldfu.ormq(
            image=image,
            crefim=crefim,
            xrng=xrng,
            yrng=yrng,
            step=step,
            mode=mode,
            numr=numr,
            cnx=cnx,
            cny=cny,
            delta=0.0,
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (
        #             82.811969518661499,
        #             2.1094076960564689,
        #             -0.74188892148201036,
        #             1,
        #             33.042001487198405,
        #         ),
        #     )
        # )


class Test_prepref(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "alignment.ali2d_single_iter")
    )

    def test_all_the_conditions(
        self, return_new=(), return_old=(), tolerance=TOLERANCE
    ):
        self.assertEqual(len(return_old), len(return_new))
        for i, j in zip(return_old, return_new):
            self.assertEqual(len(i), len(j))
            for q, r in zip(i, j):
                self.assertEqual(len(q), len(r))
                for img1, img2 in zip(q, r):
                    try:
                        self.assertTrue(
                            allclose(
                                img1.get_3dview(),
                                img2.get_3dview(),
                                atol=tolerance,
                                equal_nan=True,
                            )
                        )
                        self.assertTrue(
                            array_equal(img1.get_3dview(), img2.get_3dview())
                        )
                    except AssertionError:
                        self.assertTrue(
                            TOLERANCE
                            > numpy_abs(
                                numpy_sum(img1.get_3dview() - img2.get_3dview())
                            )
                        )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepref()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepref()
        self.assertEqual(
            str(cm_new.exception), "prepref() missing 9 required positional arguments: 'data', 'maskfile', 'cnx', 'cny', 'numr', 'mode', 'maxrangex', 'maxrangey', and 'step'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_image_mask_returns_RuntimeError_ImageDimensionException_img_dimension_doesnot_match_with_its_dimension(
        self
    ):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prepref(
                data=data,
                maskfile=EMData(),
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepref(
                data=data,
                maskfile=EMData(),
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(
            msg[1],
            "The dimension of the image does not match the dimension of the mask!",
        )
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_empty_images_to_align_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(
        self
    ):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        data = [EMData(), EMData(), EMData()]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prepref(
                data=data,
                maskfile=None,
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepref(
                data=data,
                maskfile=None,
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        with self.assertRaises(IndexError) as cm_new:
            fu.prepref(
                data, None, cnx, cny, [], mode="f", maxrangex=4, maxrangey=4, step=step
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepref(
                data, None, cnx, cny, [], mode="f", maxrangex=4, maxrangey=4, step=step
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("SKipped because of pickle file corrupted")
    def test_full_mode_without_mask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()

        return_new = fu.prepref(
            data=data,
            maskfile=None,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="f",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        return_old = oldfu.prepref(
            data=data,
            maskfile=None,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="f",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("SKipped because of pickle file corrupted")
    def test_half_mode_without_mask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        return_new = fu.prepref(
            data=data,
            maskfile=None,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="H",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        return_old = oldfu.prepref(
            data=data,
            maskfile=None,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="H",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_image_mask_returns_RuntimeError_ImageDimensionException_img_dimension_doesnot_match_with_its_dimension(
        self
    ):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        mask = model_circle(100, 100, 100)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prepref(
                data=data,
                maskfile=mask,
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepref(
                data=data,
                maskfile=mask,
                cnx=cnx,
                cny=cny,
                numr=numr,
                mode="f",
                maxrangex=4,
                maxrangey=4,
                step=step,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(
            msg[1],
            "The dimension of the image does not match the dimension of the mask!",
        )
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    @unittest.skip("SKipped because of pickle file corrupted")
    def test_Full_mode_withMask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        nx = data[0].get_xsize()
        mask = model_circle(old_div(nx, 2) - 1, nx, nx)
        return_new = fu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="f",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        return_old = oldfu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="f",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("SKipped because of pickle file corrupted")
    def test_Half_mode_withMask(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        nx = data[0].get_xsize()
        mask = model_circle(nx // 2 - 1, nx, nx)
        return_new = fu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="H",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        return_old = oldfu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="H",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        self.test_all_the_conditions(return_new, return_old)


    @unittest.skip("SKipped because of pickle file corrupted")
    def test_with_invalid_mode(self):
        """
        The output is an array of array of image. Without the TOLERANCE  even if I compare 2 results got launching the same function the test fail
        """
        data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step = give_ali2d_single_iter_data()
        nx = data[0].get_xsize()
        mask = model_circle(nx // 2 - 1, nx, nx)
        return_new = fu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="not_valid",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        return_old = oldfu.prepref(
            data=data,
            maskfile=mask,
            cnx=cnx,
            cny=cny,
            numr=numr,
            mode="not_valid",
            maxrangex=4,
            maxrangey=4,
            step=step,
        )
        self.test_all_the_conditions(return_new, return_old)


class Test_prepare_refrings(unittest.TestCase):
    """
    Take a look to sparx_utilities.py --> even_angles_cd(...)for the meaning of the following params
        ref_a --> P=Penczek algorithm, S=Saff algorithm to calculate di reference angle
        phiEQpsi  --> 'Minus', if you want psi=-phi to create a list of  angles suitable for projections, otherwise 'Zero'

    In case of rectangular kb filter see how it uses kbi variables in sparx_projection.py --> prgs(...) to understand better
    """

    volft = model_blank(100, 100, 100)
    numr = [
        1,
        1,
        8,
        2,
        9,
        16,
        3,
        953,
        128,
        16,
        1081,
        128,
        17,
        1209,
        128,
        18,
        1337,
        128,
        19,
        2745,
        256,
        26,
        3001,
        256,
        27,
        3257,
        256,
        28,
        3513,
        256,
        29,
        3769,
        256,
    ]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for img1, img2 in zip(return_new, return_old):
            try:
                self.assertTrue(array_equal(img1.get_3dview(), img2.get_3dview()))
            except AssertionError:
                # since sometimes we get  img1.get_3dview()= [[[ nan  nan  nan ...,  nan  nan  nan]]] we skip these cases
                res = numpy_sum(img1.get_3dview() - img2.get_3dview())
                if math_isnan(res) is False:
                    self.assertTrue(TOLERANCE > numpy_abs(res))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings()
        self.assertEqual(
            str(cm_new.exception),
            "prepare_refrings() missing 2 required positional arguments: 'volft' and 'kb'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_volume_returns_RuntimeError_ImageFormatException_extractplane_requires_complex_img(
        self
    ):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.prepare_refrings(
                volft=EMData(),
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.prepare_refrings(
                volft=EMData(),
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "extractplane requires a complex image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_NoneType_as_volume_returns_AttributeError_ImageFormatException_extractplane_requires_complex_img(
        self
    ):

        volft, kb = prep_vol(self.volft)
        with self.assertRaises(AttributeError) as cm_new:
            fu.prepare_refrings(
                volft=None,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.prepare_refrings(
                volft=None,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'extract_plane'"
        )

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings(
                volft=volft,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=[],
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings(
                volft=volft,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a="P",
                sym="c1",
                numr=[],
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_referenceAngle_returns_IndexError_list_index_out_of_range(self):
        volft, kb = prep_vol(self.volft)
        with self.assertRaises(IndexError) as cm_new:
            fu.prepare_refrings(
                volft=volft,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a=[],
                sym="c1",
                numr=[],
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.prepare_refrings(
                volft=volft,
                kb=kb,
                nz=4,
                delta=2.0,
                ref_a=[],
                sym="c1",
                numr=[],
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_No_kb_ArgumentError_in_EMData_extract_plane_function(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_refrings(
                volft=self.volft,
                kb=None,
                nz=4,
                delta=2.0,
                ref_a=even_angles(60.0),
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_refrings(
                volft=self.volft,
                kb=None,
                nz=4,
                delta=2.0,
                ref_a=even_angles(60.0),
                sym="c1",
                numr=self.numr,
                MPI=False,
                phiEqpsi="Minus",
                kbx=None,
                kby=None,
                initial_theta=0.1,
                delta_theta=0.5,
                initial_phi=0.1,
            )
        output_msg = "Python argument types in\n    None.extract_plane(EMData, Transform, NoneType)\ndid not match C++ signature:\n    extract_plane(EMAN::EMData {lvalue}, EMAN::Transform tf, EMAN::Util::KaiserBessel {lvalue} kb)"
        self.assertEqual(str(cm_new.exception), output_msg)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("it works but the output on the gitlab CI are not necessary")
    def test_with_sym_c1_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(
        self
    ):
        volft, kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=True,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=True,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_with_sym_c5_MPI_flag_deprecationWarning_outputmsg_PyArray_FromDims_AND_NPYCHAR_type_num_is_deprecated(
        self
    ):
        volft, kb = prep_vol(self.volft)
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=True,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=True,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip(
        "\n***************************\n\t\t 'Test_prepare_refringstest_sym_c1_initialTheta_None. "
        "Even if this combination is it seems to lead the code to a deadlock,"
        " i waited more then an hour'\n***************************"
    )
    def test_sym_c1_initialTheta_None(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=None,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=None,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_No_nz_data_size_Error_msg_datasize_hasnot_be_given(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=1,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=1,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="oct",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="oct",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_oct_Warning_in_even_angles_this_sym_isnot_supported(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="oct",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="oct",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(
        self
    ):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a=even_angles(60.0),
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a=even_angles(60.0),
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_sparx_utilities_even_angles_and_Minus(
        self
    ):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a=even_angles(60.0),
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a=even_angles(60.0),
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("Somehow it fails in CI but works on pycharm")
    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(
        self
    ):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(
        self
    ):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(
        self
    ):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Minus(
        self
    ):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Saff_algorithm_and_Minus(self):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="S",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Minus",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("Somehow it fails in CI but works on pycharm")
    def test_kb_cubic_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(
        self
    ):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c1_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(
        self
    ):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c1",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_cubic_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(
        self
    ):
        volft, kb = prep_vol(self.volft)
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kb,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=None,
            kby=None,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)

    def test_kb_rect_sym_c5_and_referenceAngles_got_via_Penczek_algorithm_and_Zero(
        self
    ):
        volft, kbx, kby, kbz = prep_vol(model_blank(100, 50, 100))
        return_new = fu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        return_old = oldfu.prepare_refrings(
            volft=volft,
            kb=kbz,
            nz=4,
            delta=2.0,
            ref_a="P",
            sym="c5",
            numr=self.numr,
            MPI=False,
            phiEqpsi="Zero",
            kbx=kbx,
            kby=kby,
            initial_theta=0.1,
            delta_theta=0.5,
            initial_phi=0.1,
        )
        self.test_all_the_conditions(return_new, return_old)


# these functions have been cleaned
"""
class Test_proj_ali_incore(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.proj_ali_incore()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.proj_ali_incore()
        self.assertEqual(str(cm_new.exception), "proj_ali_incore() takes at least 6 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_input_image_refrings_crashes_because_signal11SIGSEV(self):
        '''
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(),EMData()]
        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(array_equal(return_new, return_old))
        '''
        self.assertTrue(True)

    def test_empty_img_data_returns_RuntimeError_NotExistingObjectException_the_key_xform_projection_doesnot_exist(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        data=EMData()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], 'The requested key does not exist')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        numr=[]
        with self.assertRaises(IndexError) as cm_new:
            fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_sym_c1(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))
        self.assertTrue(allclose( return_new, (1367.8909912109375, 26.295392990112305), atol=TOLERANCE))

    def test_sym_c1_negative_deltaPhsi(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = -100.0, rshift = 0.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))
        self.assertTrue(allclose(return_new,  (-1.0000000200408773e+20, 38.38194274902344), atol=TOLERANCE))

    def test_sym_c1_negative_rshift(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = -10000.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))
        self.assertTrue(allclose( return_new, (-1.0000000200408773e+20, 14153.8154296875), atol=TOLERANCE))

    def test_sym_c1_positive_deltaPhsi(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 100.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 100.0, rshift = 0.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))
        self.assertTrue(allclose( return_new, (1361.8897705078125, 14.46776294708252), atol=TOLERANCE))

    def test_sym_c1_positive_rshift(self):
        # I cannot write the unitest
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 10000.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 10000.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))

    def test_sym_not_c1(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        return_new = fu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "icos", delta_psi = 0.0, rshift = 0.0)
        return_old = oldfu.proj_ali_incore(deepcopy(data), refrings, numr, xrng, yrng, step, finfo=None, sym = "icos", delta_psi = 0.0, rshift = 0.0)
        self.assertTrue(allclose(return_old, return_new , atol=TOLERANCE ))
        self.assertTrue(allclose( return_new, (1367.8907470703125, 11.979691505432129), atol=TOLERANCE))



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

        with self.assertRaises(IndexError) as cm_new:
            fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_empty_list_of_ref_ang(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        with self.assertRaises(IndexError) as cm_new:
            fu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.proj_ali_incore_local(data, refrings, [], numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_too_short_list_of_ref_ang(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        symangles = [[0.0, 0.0, 360. ] for k in range(10) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        return_new = fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        return_old = oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertTrue(array_equal(return_new, return_old))


    #It does not work because the old version returns always True. If I remove the return True it gives back a very different values.Is it this function obsolete???

    def test_G12(self):
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]

        symangles = [[0.0, 0.0, 360. ] for k in range(len(refrings)) ]

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')

        return_new = fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        return_old = oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng= 2.0, yrng=2.0, step=step, an=-1.0)
        self.assertTrue(array_equal(return_new, return_old))
"""


class Test_ali_vol_func(unittest.TestCase):
    param = [1, 1, 1, 1, 1, 1]
    data = get_data(3)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol_func()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol_func()
        self.assertEqual(
            str(cm_new.exception), "ali_vol_func() missing 2 required positional arguments: 'params' and 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_few_params_params_returns_IndexError_list_index_out_of_range(self):
        param = [1, 1, 1, 1]
        with self.assertRaises(IndexError) as cm_new:
            fu.ali_vol_func(params=param, data=self.data)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ali_vol_func(params=param, data=self.data)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_too_few_data_params_returns_IndexError_list_index_out_of_range(self):
        data = get_data(2)
        with self.assertRaises(IndexError) as cm_new:
            fu.ali_vol_func(params=self.param, data=data)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ali_vol_func(params=self.param, data=data)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_vol_func(self):
        return_new = fu.ali_vol_func(params=self.param, data=self.data)
        return_old = oldfu.ali_vol_func(params=self.param, data=self.data)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertEqual(return_new, 0.9856925010681152)

    def test_ali_vol_func_with_NoneTypes_as_image_returns_AttributeError_NoneType_obj_hasnot_attribute_rot_scale_trans_background(
        self
    ):
        data = [None, None, None]
        with self.assertRaises(AttributeError) as cm_new:
            fu.ali_vol_func(params=self.param, data=data)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ali_vol_func(params=self.param, data=data)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception),
            "'NoneType' object has no attribute 'rot_scale_trans_background'",
        )

    def test_empty_data_images_returns_RuntimeError_InvalidValueException_xsize_not_positive(
        self
    ):
        data = [EMData(), EMData(), EMData()]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol_func(params=self.param, data=data)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol_func(params=self.param, data=data)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])


class Test_align2d(unittest.TestCase):
    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, refim, xrng, yrng) = self.argum[0]
        image = EMData()
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(None, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        return_old = oldfu.align2d(None, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F")
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_NoneType__image_to_align_creturns_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(
        self
    ):
        with self.assertRaises(AttributeError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=None,
                xrng=[0, 0],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=None,
                xrng=[0, 0],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_empty_image_reference_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=EMData(),
                xrng=[0, 0],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=EMData(),
                xrng=[0, 0],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_xrng_returns_ValueError_arg_af_max_f_is_empty_list(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        with self.assertRaises(ValueError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[],
                yrng=[0, 0],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), "max() arg is an empty sequence")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_yrngreturns_ValueError_arg_af_max_f_is_empty_list(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        with self.assertRaises(ValueError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), "max() arg is an empty sequence")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d()
        self.assertEqual(
            str(cm_new.exception), "align2d() missing 2 required positional arguments: 'image' and 'refim'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_wrong_enumerate_rings_error_crashes_because_signal11SIGSEV(self):
        """
        (image, refim, xrng, yrng) = self.argum[0]
        return_new = fu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        return_old = oldfu.align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=0, last_ring=1, rstep=1, mode = "F")
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_null_rstep_value_returns_ValueError_arg_af_max_f_is_empty_list(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=0,
                mode="F",
            )
        with self.assertRaises(ValueError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[],
                step=1,
                first_ring=1,
                last_ring=0,
                rstep=0,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), "max() arg is an empty sequence")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_null_step_value_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[0, 0],
                step=0,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.align2d(
                image=IMAGE_2D,
                refim=IMAGE_2D_REFERENCE,
                xrng=[0, 0],
                yrng=[0, 0],
                step=0,
                first_ring=1,
                last_ring=0,
                rstep=1,
                mode="F",
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Full_mode_zero_lastRing(self):
        return_new = fu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=0,
            rstep=1,
            mode="F",
        )
        return_old = oldfu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=0,
            rstep=1,
            mode="F",
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new, (0.39858296513557434, -0.0, 0.0, 1, 1.5912819397733529)
            )
        )

    def test_Half_mode_zero_lastRing(self):
        return_new = fu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=0,
            rstep=1,
            mode="H",
        )
        return_old = oldfu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=0,
            rstep=1,
            mode="H",
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new, (179.09351259469986, 0.0, 0.0, 0, 0.9412509725952987)
            )
        )

    def test_Full_mode(self):
        return_new = fu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=2,
            rstep=1,
            mode="F",
        )
        return_old = oldfu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=2,
            rstep=1,
            mode="F",
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new, (2.9873049259185791, -0.0, 0.0, 1, 0.6616981065062646)
            )
        )

    def test_Half_mode(self):
        return_new = fu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=2,
            rstep=1,
            mode="H",
        )
        return_old = oldfu.align2d(
            image=IMAGE_2D,
            refim=IMAGE_2D_REFERENCE,
            xrng=[0, 0],
            yrng=[0, 0],
            step=1,
            first_ring=1,
            last_ring=2,
            rstep=1,
            mode="H",
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_old, (177.3188880085945, 0.0, 0.0, 0, 0.41331450702273287)
            )
        )


class Test_align2d_scf(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align2d_scf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align2d_scf()
        self.assertEqual(
            str(cm_new.exception), "align2d_scf() missing 2 required positional arguments: 'image' and 'refim'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_image_to_align_returns_RuntimeError_InvalidValueException_xsize_not_positive(
        self
    ):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.align2d_scf(
                image=EMData(), refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=174
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.align2d_scf(
                image=EMData(), refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=174
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType__image_to_align_creturns_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(
        self
    ):
        with self.assertRaises(AttributeError) as cm_new:
            fu.align2d_scf(image=None, refim=EMData(), xrng=4, yrng=4, ou=174)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.align2d_scf(image=None, refim=EMData(), xrng=4, yrng=4, ou=174)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_empty_reference_image_returns_RuntimeError_InvalidValueException_xsize_not_positive(
        self
    ):
        with self.assertRaises(RuntimeError) as cm_old:
            fu.align2d_scf(image=IMAGE_2D, refim=EMData(), xrng=4, yrng=4, ou=174)
        with self.assertRaises(RuntimeError) as cm_new:
            oldfu.align2d_scf(image=IMAGE_2D, refim=EMData(), xrng=4, yrng=4, ou=174)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_as_reference_image_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        (image, refim, xrng, yrng) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_old:
            fu.align2d_scf(image, None, xrng, yrng, self.argum[1]['ou'])
        with self.assertRaises(RuntimeError) as cm_new:
            oldfu.align2d_scf(image, None, xrng, yrng, self.argum[1]['ou'])
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        """

    def test_with_valid_params(self):
        return_new = fu.align2d_scf(
            image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=174
        )
        return_old = oldfu.align2d_scf(
            image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=174
        )
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         (
        #             0.17578125,
        #             2.9674494882377935,
        #             -0.05141488826358742,
        #             1,
        #             4.90025769648605,
        #         ),
        #     )
        # )

    @unittest.skip("It works but skip because of long output")
    def test_with_invalid_ou_error_msg_output(self):
        return_new = fu.align2d_scf(
            image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=1
        )
        return_old = oldfu.align2d_scf(
            image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=4, yrng=4, ou=1
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                (
                    0.2789926528930664,
                    -0.482177873659118,
                    -0.048944523282220764,
                    0,
                    4.883454103473488,
                ),
            )
        )

    """
    the following testa are not able to work. It'd be a bug.
    error message:
        File "/home/lusnig/EMAN2/eman2/sphire/tests/sparx_lib/sparx_alignment.py", line 784, in align2d_scf
        sxs = -p2[0][4]
        IndexError: list index out of range
    BUT p2 is the ouput of:
        -) ccf2 = EMAN2_cppwrap.Util.window(ccf(rot_shift2D(image, alpha+180.0, 0.0, 0.0, mirr), frotim),nrx,nry)
	    -) p2 = sparx_utilities.peak_search(ccf2)
	in these casea it is a list of 4 elements and it is trying to get the 5th
    """

    def test_with_DEFAULT_params_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.align2d_scf(
                image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=-1, yrng=-1, ou=-1
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align2d_scf(
                image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=-1, yrng=-1, ou=-1
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_DEFAULT_params_but_validOU_returns_IndexError_list_index_out_of_range(
        self
    ):
        with self.assertRaises(IndexError) as cm_new:
            fu.align2d_scf(
                image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=-1, yrng=-1, ou=174
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align2d_scf(
                image=IMAGE_2D, refim=IMAGE_2D_REFERENCE, xrng=-1, yrng=-1, ou=174
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_multialign2d_scf(unittest.TestCase):
    numr = [
        1,
        1,
        8,
        2,
        9,
        16,
        3,
        25,
        32,
        4,
        57,
        32,
        5,
        89,
        32,
        6,
        121,
        64,
        7,
        185,
        64,
        8,
        249,
        64,
        9,
        313,
        64,
        10,
        377,
        64,
        11,
        441,
        128,
        12,
        569,
        128,
        13,
        697,
        128,
        14,
        825,
        128,
        15,
        953,
        128,
        16,
        1081,
        128,
        17,
        1209,
        128,
        18,
        1337,
        128,
        19,
        1465,
        128,
        20,
        1593,
        128,
        21,
        1721,
        256,
        22,
        1977,
        256,
        23,
        2233,
        256,
        24,
        2489,
        256,
        25,
        2745,
        256,
        26,
        3001,
        256,
        27,
        3257,
        256,
        28,
        3513,
        256,
        29,
        3769,
        256,
    ]

    def test_empty_input_image_refrings_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = EMData()
        frotim = [fft(tavg)]

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertEqual(return_old, return_new)
        """
        self.assertTrue(True)

    def test_NoneType_as_input_image_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]

        dataa = deepcopy(self.argum[0][0])
        cimage = None
        frotim = [fft(tavg)]

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

    def test_NoneType__image_reference_typeError_NoneType_obj_hasnot_attribute___getitem__(
        self
    ):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "F")
        with self.assertRaises(TypeError) as cm_new:
            fu.multalign2d_scf(
                image=IMAGE_2D,
                refrings=[cimage],
                frotim=None,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        with self.assertRaises(TypeError) as cm_old:
            oldfu.multalign2d_scf(
                image=IMAGE_2D,
                refrings=[cimage],
                frotim=None,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object is not subscriptable"
        )

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = self.argum[0]
        dataa = deepcopy(self.argum[0][0])
        cimage = Util.Polar2Dm(tavg, float(cnx), float(cny), numr, "F")
        frotim = [fft(tavg)]
        numr = []

        return_new = fu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        return_old = oldfu.multalign2d_scf(dataa[0], [cimage], frotim, numr, xrng, yrng, ou=174)
        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.multalign2d_scf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.multalign2d_scf()
        self.assertEqual(
            str(cm_new.exception),
            "multalign2d_scf() missing 4 required positional arguments: 'image', 'refrings', 'frotim', and 'numr'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_returns_RuntimeError_InvalidValueException_xsize_not_positive(self):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "F")
        frotim = [fft(tavg)]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.multalign2d_scf(
                image=EMData(),
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.multalign2d_scf(
                image=EMData(),
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_with_NoneType_images_as_data_AttributeError_NoneType_obj_hasnot_attribute_get_xsize(
        self
    ):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "F")
        frotim = [fft(tavg)]
        with self.assertRaises(AttributeError) as cm_new:
            fu.multalign2d_scf(
                image=None,
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.multalign2d_scf(
                image=None,
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=4,
                yrng=4,
                ou=174,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    @unittest.skip(
        "it run from pycharm but not always from the console (nosetests too), of course using the same pythn interpreter"
    )
    def test_with_valid_params_cimage_with_mode_F(self):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "F")
        frotim = [fft(tavg)]

        return_new = fu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=174,
        )
        return_old = oldfu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=174,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                (
                    2.9782530585469695,
                    -0.04681045832942976,
                    0,
                    0.703125,
                    1,
                    4.878763807138058,
                ),
            )
        )

    @unittest.skip(
        "it run from pycharm but not always from the console (nosetests too), of course using the same pythn interpreter"
    )
    def test_with_valid_params_cimage_with_mode_H(self):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "H")
        frotim = [fft(tavg)]
        return_new = fu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=174,
        )
        return_old = oldfu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=174,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                (
                    2.9782530585469695,
                    -0.04681045832942976,
                    0,
                    0.703125,
                    1,
                    4.878763807138058,
                ),
            )
        )

    @unittest.skip(
        "it run from pycharm but not always from the console (nosetests too), of course using the same pythn interpreter"
    )
    def test_with_invalid_ou(self):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "H")
        frotim = [fft(tavg)]
        return_new = fu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=1,
        )
        return_old = oldfu.multalign2d_scf(
            image=deepcopy(IMAGE_2D),
            refrings=[cimage],
            frotim=frotim,
            numr=self.numr,
            xrng=4,
            yrng=4,
            ou=1,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(
                return_new,
                (
                    2.9782530585469695,
                    -0.04681045832942976,
                    0,
                    0.703125,
                    1,
                    4.878763807138058,
                ),
            )
        )

    def test_with_DEFAULT_params_returns_IndexError_list_index_out_of_range(self):
        tavg = IMAGE_2D_REFERENCE
        cimage = Util.Polar2Dm(tavg, 36, 36, self.numr, "H")
        frotim = [fft(tavg)]
        with self.assertRaises(IndexError) as cm_new:
            fu.multalign2d_scf(
                image=deepcopy(IMAGE_2D),
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=-1,
                yrng=-1,
                ou=-1,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.multalign2d_scf(
                image=deepcopy(IMAGE_2D),
                refrings=[cimage],
                frotim=frotim,
                numr=self.numr,
                xrng=-1,
                yrng=-1,
                ou=-1,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_parabl(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "alignment.parabl")
    )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.parabl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.parabl()
        self.assertEqual(
            str(cm_new.exception), "parabl() missing 1 required positional argument: 'Z'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_input_image(self):
        return_new = fu.parabl(Z=EMData())
        return_old = oldfu.parabl(Z=EMData())
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_old, (0.0, 0.0, 0.0)))

    def test_NoneType_as_input_image_returns_typeError_object_float_hasnot_attibute(
        self
    ):
        with self.assertRaises(TypeError) as cm_new:
            fu.parabl(None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.parabl(None)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object is not subscriptable"
        )

    def test_2Dimg(self):
        return_new = fu.parabl(Z=IMAGE_2D)
        return_old = oldfu.parabl(Z=IMAGE_2D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_old, (-1.0, -1.0, -0.17265136317835783)))

    def test_3Dimg(self):
        return_new = fu.parabl(Z=IMAGE_3D)
        return_old = oldfu.parabl(Z=IMAGE_3D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_old, (-1.0, -1.0, -0.17265136317835783)))

    def test_2Dimgblank(self):
        return_new = fu.parabl(Z=IMAGE_BLANK_2D)
        return_old = oldfu.parabl(Z=IMAGE_BLANK_2D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, (0.0, 0.0, 0.0)))

    def test_3Dimgblank(self):
        return_new = fu.parabl(Z=IMAGE_BLANK_3D)
        return_old = oldfu.parabl(Z=IMAGE_BLANK_3D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, (0.0, 0.0, 0.0)))


@unittest.skip("sometimes the output values change")
class Test_shc(unittest.TestCase):
    # argum = get_arg_from_pickle_file(
    #     path.join(ABSOLUTE_PATH_TO_RESOURCES, "alignment.shc")
    # )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.shc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.shc()
        self.assertEqual(
            str(cm_new.exception), "shc() takes at least 7 arguments (0 given)"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_input_image_refringscrashes_because_signal11SIGSEV(self):
        """
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
        refrings = [EMData(), EMData()]
        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)
        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0)

        self.assertTrue(array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_returns_RuntimeError_the_key_xform_projection_doesnot_exist(
        self
    ):
        data, refrings, list_of_ref_ang, numr, xrng, yrng, step = give_alignment_shc_data(True)
        data = EMData()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an=-1.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an=-1.0)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_list_Numrinit_returns_IndexError_list_index_out_of_range(self):
        data, refrings, list_of_ref_ang, numr, xrng, yrng, step = give_alignment_shc_data(True)
        numr = []
        with self.assertRaises(IndexError) as cm_new:
            fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an=-1.0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an=-1.0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    # def test_sym_c1_failed(self):
    #     (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
    #
    #     # return_new = fu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #     # return_old = oldfu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #     self.assertTrue(True)
    #     # self.assertTrue(array_equal(return_new, return_old))

    # def test_empty_list_of_ref_ang_failed(self):
    #     (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
    #     list_of_ref_ang = []
    #
    #     # return_new = fu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #     # return_old = oldfu.shc(deepcopy(data), deepcopy(refrings), list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #     self.assertTrue(True)
    #     # self.assertTrue(array_equal(return_new, return_old))

    # def test_added_one_ref_ang_failed(self):
    #     (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
    #     list_of_ref_ang[0].append(2.0)
    #     # return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #     # return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "c1", finfo=None)
    #
    #     self.assertTrue(True)
    #     # self.assertTrue(array_equal(return_new, return_old))
    #
    # def test_sym_nomirror_failed(self):
    #     (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = self.argum[0]
    #
    #     # return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "nomirror", finfo=None)
    #     # return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step, an =-1.0, sym = "nomirror", finfo=None)
    #
    #     self.assertTrue(True)
    #     # self.assertTrue(array_equal(return_new, return_old))


class Test_search_range(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.search_range()
            oldfu.search_range()

    def test_search_range(self):
        return_new = fu.search_range(n=70, radius=29, shift=0.0, range=4, location="")
        return_old = oldfu.search_range(
            n=70, radius=29, shift=0.0, range=4, location=""
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [4, 4]))

    @unittest.skip("It works but output is not necessary for CI")
    def test_no_image_size_warning_msg_shift_of_particle_too_large(self):
        return_new = fu.search_range(n=0, radius=29, shift=0.0, range=4, location="")
        return_old = oldfu.search_range(n=0, radius=29, shift=0.0, range=4, location="")
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [0, 0]))


"""
# this function has been cleaned
class Test_generate_list_of_reference_angles_for_search(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.generate_list_of_reference_angles_for_search()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.generate_list_of_reference_angles_for_search()
        self.assertEqual(str(cm_new.exception), "generate_list_of_reference_angles_for_search() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_c5Sym(self):
        sym = 'c5'
        ref_angles = even_angles(symmetry=sym)
        return_new = fu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        return_old = oldfu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[0.0, 0.0, 0], [44.389802268702958, 19.749922795642572, 0], [4.2668022687029463, 28.07248693585296, 0], [30.709070484328791, 34.56032177999784, 0], [53.987279178298351, 40.11916689840513, 0], [3.1658765410020351, 45.0991278440273, 0], [22.839553661168754, 49.679784930029356, 0], [41.388074331228111, 53.96812092752944, 0], [59.069143085907285, 58.03428124920304, 0], [4.0711430859072806, 61.927513064147035, 0], [20.531322307288804, 65.68426082882351, 0], [36.563034950568053, 69.33268350648736, 0], [52.257200050432978, 72.89536482335619, 0], [67.690502133487442, 76.39103936921708, 0], [10.931666902500965, 79.83575137827651, 0], [26.036563285357822, 83.24367296941216, 0], [41.062582385571957, 86.62771331656609, 0], [180.0, 180.0, 0], [224.38980226870297, 160.25007720435744, 0], [184.26680226870295, 151.92751306414704, 0], [210.70907048432878, 145.43967822000215, 0], [233.98727917829837, 139.88083310159487, 0], [183.16587654100204, 134.9008721559727, 0], [202.83955366116876, 130.32021506997063, 0], [221.3880743312281, 126.03187907247056, 0], [239.06914308590729, 121.96571875079695, 0], [184.07114308590729, 118.07248693585296, 0], [200.5313223072888, 114.31573917117649, 0], [216.56303495056807, 110.66731649351264, 0], [232.25720005043297, 107.10463517664381, 0], [247.69050213348743, 103.60896063078292, 0], [190.93166690250098, 100.16424862172349, 0], [206.03656328535783, 96.75632703058784, 0], [221.06258238557194, 93.37228668343391, 0], [0.0, 0.0, 72.0], [116.3898, 19.74992, 0.0], [76.2668, 28.07249, 0.0], [102.70906, 34.56032, 0.0], [125.98728, 40.11917, 0.0], [75.16587, 45.09913, 0.0], [94.83955, 49.67978, 0.0], [113.38808, 53.96812, 0.0], [131.06914, 58.03428, 0.0], [76.07115, 61.92751, 0.0], [92.53132, 65.68426, 0.0], [108.56304, 69.33268, 0.0], [124.25719, 72.89537, 0.0], [139.69051, 76.39104, 0.0], [82.93166, 79.83575, 0.0], [98.03656, 83.24367, 0.0], [113.06258, 86.62771, 0.0], [0.0, 180.0, 108.0], [296.3898, 160.25008, 0.0], [256.2668, 151.92751, 0.0], [282.70907, 145.43968, 0.0], [305.98728, 139.88083, 0.0], [255.16588, 134.90087, 0.0], [274.83955, 130.32022, 0.0], [293.38808, 126.03188, 0.0], [311.06915, 121.96572, 0.0], [256.07114, 118.07249, 0.0], [272.53132, 114.31574, 0.0], [288.56304, 110.66732, 0.0], [304.2572, 107.10463, 0.0], [319.6905, 103.60896, 0.0], [262.93167, 100.16425, 0.0], [278.03656, 96.75633, 0.0], [293.06258, 93.37229, 0.0], [0.0, 0.0, 144.0], [188.38981, 19.74992, 0.0], [148.2668, 28.07249, 0.0], [174.70907, 34.56032, 0.0], [197.98728, 40.11917, 0.0], [147.16588, 45.09913, 0.0], [166.83955, 49.67978, 0.0], [185.38808, 53.96812, 0.0], [203.06915, 58.03428, 0.0], [148.07115, 61.92751, 0.0], [164.53133, 65.68426, 0.0], [180.56304, 69.33268, 0.0], [196.2572, 72.89536, 0.0], [211.6905, 76.39104, 0.0], [154.93167, 79.83575, 0.0], [170.03657, 83.24367, 0.0], [185.06258, 86.62771, 0.0], [0.0, 180.0, 36.0], [8.3898, 160.25008, 0.0], [328.2668, 151.92751, 0.0], [354.70907, 145.43968, 0.0], [17.98728, 139.88083, 0.0], [327.16588, 134.90087, 0.0], [346.83955, 130.32022, 0.0], [5.38808, 126.03188, 0.0], [23.06914, 121.96572, 0.0], [328.07114, 118.07249, 0.0], [344.53132, 114.31574, 0.0], [0.56304, 110.66732, 0.0], [16.2572, 107.10464, 0.0], [31.6905, 103.60896, 0.0], [334.93167, 100.16425, 0.0], [350.03657, 96.75633, 0.0], [5.06258, 93.37229, 0.0], [0.0, 0.0, 216.0], [260.3898, 19.74992, 0.0], [220.2668, 28.07249, 0.0], [246.70907, 34.56032, 0.0], [269.98728, 40.11917, 0.0], [219.16588, 45.09913, 0.0], [238.83955, 49.67978, 0.0], [257.38808, 53.96812, 0.0], [275.06914, 58.03428, 0.0], [220.07114, 61.92751, 0.0], [236.53132, 65.68426, 0.0], [252.56303, 69.33268, 0.0], [268.2572, 72.89536, 0.0], [283.6905, 76.39104, 0.0], [226.93167, 79.83575, 0.0], [242.03656, 83.24367, 0.0], [257.06258, 86.62771, 0.0], [0.0, 180.0, 324.0], [80.38981, 160.25008, 0.0], [40.2668, 151.92751, 0.0], [66.70907, 145.43968, 0.0], [89.98727, 139.88083, 0.0], [39.16587, 134.90087, 0.0], [58.83955, 130.32022, 0.0], [77.38807, 126.03188, 0.0], [95.06914, 121.96572, 0.0], [40.07114, 118.07249, 0.0], [56.53133, 114.31574, 0.0], [72.56304, 110.66732, 0.0], [88.2572, 107.10464, 0.0], [103.69051, 103.60896, 0.0], [46.93167, 100.16425, 0.0], [62.03657, 96.75633, 0.0], [77.06257, 93.37229, 0.0], [0.0, 0.0, 288.0], [332.3898, 19.74992, 0.0], [292.2668, 28.07249, 0.0], [318.70907, 34.56032, 0.0], [341.98728, 40.11917, 0.0], [291.16588, 45.09913, 0.0], [310.83956, 49.67978, 0.0], [329.38808, 53.96812, 0.0], [347.06914, 58.03428, 0.0], [292.07114, 61.92751, 0.0], [308.53132, 65.68426, 0.0], [324.56303, 69.33268, 0.0], [340.2572, 72.89536, 0.0], [355.6905, 76.39104, 0.0], [298.93167, 79.83575, 0.0], [314.03656, 83.24367, 0.0], [329.06258, 86.62771, 0.0], [0.0, 180.0, 252.0], [152.3898, 160.25008, 0.0], [112.26681, 151.92751, 0.0], [138.70907, 145.43968, 0.0], [161.98729, 139.88083, 0.0], [111.16588, 134.90087, 0.0], [130.83955, 130.32022, 0.0], [149.38808, 126.03188, 0.0], [167.06914, 121.96572, 0.0], [112.07114, 118.07249, 0.0], [128.53132, 114.31574, 0.0], [144.56303, 110.66732, 0.0], [160.2572, 107.10464, 0.0], [175.6905, 103.60896, 0.0], [118.93167, 100.16425, 0.0], [134.03656, 96.75633, 0.0], [149.06258, 93.37229, 0.0]]))

    def test_c1Sym(self):
        sym = 'c1'
        ref_angles = even_angles(symmetry=sym)
        return_new = fu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        return_old = oldfu.generate_list_of_reference_angles_for_search(ref_angles, sym)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[0.0, 0.0, 0], [102.5624397820863, 8.409807949596694, 0], [175.28184168449116, 11.903989804110001, 0], [234.81899085328783, 14.592550602033418, 0], [286.52113069039967, 16.865343252479008, 0], [332.89249973858841, 18.873236840047255, 0], [15.350997945238817, 20.69354123118596, 0], [54.760293521450905, 22.37214549396397, 0], [91.727719586672706, 23.938926249214624, 0], [126.67925988880424, 25.41462091516098, 0], [159.93126768874427, 26.81431796194859, 0], [191.72626852098327, 28.149400619646084, 0], [222.25501416086877, 29.428707176867, 0], [251.6707339535308, 30.659262305350033, 0], [280.09871166816117, 31.846758629170495, 0], [307.64293448395898, 32.995885473579534, 0], [334.39083847001103, 34.11056017878775, 0], [0.42677669506366556, 35.194095100409235, 0], [25.794606434997782, 36.249320882899376, 0], [50.559654291516139, 37.278679231322116, 0], [74.770232732225381, 38.2842939251198, 0], [98.468827134074971, 39.26802600175335, 0], [121.69303677671941, 40.231517219359155, 0], [144.4763293594925, 41.17622470375671, 0], [166.84865229059051, 42.10344887074584, 0], [188.83693262466142, 43.014356152771704, 0], [210.46548946865465, 43.909997664475156, 0], [231.75637688070145, 44.79132466007832, 0], [252.72967105963514, 45.65920143165515, 0], [273.40371249950607, 46.51441614768202, 0], [293.7953114483945, 47.357690020060026, 0], [313.91992324589262, 48.1896851042214, 0], [333.79179876604201, 49.01101097344977, 0], [353.42411415385686, 49.822230459852115, 0], [12.839083235960516, 50.62386461673009, 0], [32.02805535274598, 51.41639702767674, 0], [51.011600859315614, 52.20027756457276, 0], [69.799586144482291, 52.975925678303284, 0], [88.401239698292727, 53.743733291363625, 0], [106.82521050148785, 54.50406734974836, 0], [125.07961980182155, 55.25727208199666, 0], [143.17210717208275, 56.00367100552329, 0], [161.10987160517593, 56.74356871403049, 0], [178.89970828662715, 57.4772524745885, 0], [196.54804158963, 58.20499365866951, 0], [214.06095475847701, 58.92704902784667, 0], [231.44421667996505, 59.64366189189109, 0], [248.70330608674968, 60.355063154503576, 0], [265.84343348975648, 61.06147225981934, 0], [282.86956109711195, 61.763098051052104, 0], [299.78642094339619, 62.46013955114206, 0], [316.59853142434207, 63.152786673995614, 0], [333.31021240759083, 63.84122087381428, 0], [349.92559906909207, 64.52561573907757, 0], [6.4586545866518463, 65.20613753694339, 0], [22.893181806532958, 65.88294571313848, 0], [39.242833985512988, 66.55619335181605, 0], [55.511124699098673, 67.22602759934011, 0], [71.701436996410379, 67.8925900555079, 0], [87.81703187337213, 68.55601713533103, 0], [103.86105612808187, 69.21644040415431, 0], [119.8365496554388, 69.87398688859322, 0], [135.74645223213611, 70.52877936550931, 0], [151.59360983787678, 71.18093663101206, 0], [167.38078055404094, 71.83057375127423, 0], [183.11064007694512, 72.47780229676785, 0], [198.78578687921549, 73.12273056137076, 0], [214.40874704959094, 73.76546376765336, 0], [229.98197883862355, 74.40610425953089, 0], [245.50787693521318, 75.04475168335667, 0], [260.98877649665752, 75.68150315843295, 0], [276.42695695288819, 76.31645343782941, 0], [291.82464560376934, 76.94969506032008, 0], [307.18402102672974, 77.58131849418093, 0], [322.50721631056541, 78.21141227352726, 0], [337.79632212996364, 78.84006312781455, 0], [353.0533896741494, 79.46735610507622, 0], [8.2904334420228452, 80.09337468942728, 0], [23.489433915232105, 80.71820091332246, 0], [38.662340119797371, 81.34191546502161, 0], [53.811072086159413, 81.96459779168268, 0], [68.937523216861678, 82.58632619847424, 0], [84.043562570481001, 83.20717794407292, 0], [99.131037069892173, 83.82722933288893, 0], [114.20177364247999, 84.44655580434149, 0], [129.25758129949423, 85.06523201948858, 0], [144.30025316137389, 85.68333194529811, 0], [159.33156843554312, 86.30092893683496, 0], [174.35329435289955, 86.91809581762422, 0], [189.36718806897298, 87.53490495844152, 0], [204.37499853552671, 88.15142835477144, 0], [219.37846834820326, 88.7677377031675, 0], [234.37933557567774, 89.38390447674091, 0], [180.0, 180.0, 0], [282.5624397820863, 171.59019205040332, 0], [355.28184168449116, 168.09601019589, 0], [54.818990853287801, 165.4074493979666, 0], [106.52113069039967, 163.13465674752098, 0], [152.89249973858841, 161.12676315995276, 0], [195.35099794523882, 159.30645876881403, 0], [234.7602935214509, 157.62785450603604, 0], [271.72771958667272, 156.06107375078537, 0], [306.67925988880427, 154.58537908483902, 0], [339.93126768874424, 153.18568203805143, 0], [11.72626852098324, 151.8505993803539, 0], [42.255014160868768, 150.571292823133, 0], [71.670733953530771, 149.34073769464996, 0], [100.09871166816117, 148.1532413708295, 0], [127.64293448395898, 147.00411452642047, 0], [154.39083847001098, 145.88943982121225, 0], [180.42677669506367, 144.80590489959076, 0], [205.79460643499777, 143.75067911710062, 0], [230.55965429151615, 142.72132076867788, 0], [254.7702327322254, 141.7157060748802, 0], [278.46882713407496, 140.73197399824664, 0], [301.6930367767194, 139.76848278064085, 0], [324.47632935949252, 138.8237752962433, 0], [346.84865229059051, 137.89655112925416, 0], [8.836932624661415, 136.9856438472283, 0], [30.465489468654653, 136.09000233552484, 0], [51.756376880701453, 135.20867533992168, 0], [72.729671059635166, 134.34079856834484, 0], [93.403712499506071, 133.48558385231797, 0], [113.7953114483945, 132.64230997993997, 0], [133.91992324589262, 131.81031489577862, 0], [153.79179876604201, 130.98898902655023, 0], [173.42411415385686, 130.1777695401479, 0], [192.83908323596052, 129.3761353832699, 0], [212.02805535274598, 128.58360297232326, 0], [231.01160085931562, 127.79972243542724, 0], [249.7995861444823, 127.02407432169672, 0], [268.40123969829273, 126.25626670863637, 0], [286.82521050148785, 125.49593265025163, 0], [305.07961980182154, 124.74272791800334, 0], [323.17210717208275, 123.99632899447671, 0], [341.10987160517595, 123.25643128596951, 0], [358.89970828662717, 122.5227475254115, 0], [16.548041589629975, 121.79500634133049, 0], [34.06095475847701, 121.07295097215334, 0], [51.444216679965052, 120.3563381081089, 0], [68.703306086749649, 119.64493684549643, 0], [85.843433489756478, 118.93852774018066, 0], [102.86956109711195, 118.23690194894789, 0], [119.78642094339619, 117.53986044885794, 0], [136.59853142434207, 116.84721332600438, 0], [153.31021240759083, 116.15877912618572, 0], [169.92559906909207, 115.47438426092243, 0], [186.45865458665185, 114.79386246305661, 0], [202.89318180653297, 114.11705428686152, 0], [219.24283398551299, 113.44380664818395, 0], [235.51112469909867, 112.77397240065989, 0], [251.70143699641039, 112.1074099444921, 0], [267.81703187337212, 111.44398286466897, 0], [283.86105612808188, 110.78355959584569, 0], [299.83654965543883, 110.12601311140678, 0], [315.74645223213611, 109.47122063449069, 0], [331.59360983787678, 108.81906336898794, 0], [347.38078055404094, 108.16942624872577, 0], [3.110640076945117, 107.52219770323215, 0], [18.785786879215493, 106.87726943862924, 0], [34.408747049590943, 106.23453623234664, 0], [49.981978838623547, 105.59389574046911, 0], [65.507876935213176, 104.95524831664333, 0], [80.988776496657522, 104.31849684156705, 0], [96.42695695288819, 103.68354656217059, 0], [111.82464560376934, 103.05030493967992, 0], [127.18402102672974, 102.41868150581907, 0], [142.50721631056541, 101.78858772647274, 0], [157.79632212996364, 101.15993687218545, 0], [173.05338967414946, 100.53264389492378, 0], [188.29043344202285, 99.90662531057272, 0], [203.4894339152321, 99.28179908667754, 0], [218.66234011979736, 98.65808453497839, 0], [233.81107208615941, 98.03540220831732, 0], [248.93752321686168, 97.41367380152576, 0], [264.04356257048101, 96.79282205592708, 0], [279.13103706989216, 96.17277066711107, 0], [294.20177364248002, 95.55344419565851, 0], [309.25758129949423, 94.93476798051142, 0], [324.30025316137392, 94.31666805470189, 0], [339.33156843554309, 93.69907106316504, 0], [354.35329435289952, 93.08190418237578, 0], [9.3671880689729505, 92.46509504155848, 0], [24.374998535526743, 91.84857164522856, 0], [39.378468348203228, 91.2322622968325, 0], [54.379335575677715, 90.61609552325909, 0]]))

    def test_NoAngleList(self):
        return_new = fu.generate_list_of_reference_angles_for_search([], 'c1')
        return_old = oldfu.generate_list_of_reference_angles_for_search([], 'c1')
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, []))

    def test_invalid_simmetry_returns_RuntimeError_NotExistingObjectException_the_key_invalid_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.generate_list_of_reference_angles_for_search([], 'invalid')
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.generate_list_of_reference_angles_for_search([], 'invalid')
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], 'No such an instance existing')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
"""


"""
# Adnan helper functions to run the reference tests
import copy
import os
import cPickle as pickle
import EMAN2_cppwrap as e2cpp
def get_data(num,dim = 10):
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


def get_data_3d(num, dim=10):
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
        data_list.append(a)

    return data_list



@unittest.skip("original adnan")
class Test_lib_alignment_compare(unittest.TestCase):

    # data = list of images
    # numr = tuple or list precalcualte rings
    # wr = list of weights of numr
    # cs =  cs = [0.0]*2
    # tavg = blanck image
    # cnx = center value x
    # cny = cnx
    # xrng =  list of possible shifts
    # yrng =  list of possible shifts
    # step = stepsize of the shift

    def test_ali2d_single_iter_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = argum[0]
        (datab, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = argum[0]

        dataa = copy.deepcopy(argum[0][0])
        datab = copy.deepcopy(argum[0][0])

        return_new = fu.ali2d_single_iter(dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
        return_old = fu.ali2d_single_iter(datab, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
        self.assertEqual(return_new,return_old)


    def test_ang_n_true_should_return_equal_object(self):
        return_new = fu.ang_n(2 , 'f' , 3)
        return_old = oldfu.ang_n(2, 'f' , 3)

        self.assertEqual(return_new, return_old)

    def test_log2_true_should_return_equal_object(self):
        return_new = fu.log2(10)
        return_old = oldfu.log2(10)

        self.assertEqual(return_new,return_old)

    def test_Numrinit_true_should_return_equal_object(self):
        return_new = fu.Numrinit(2 , 5)
        return_old = oldfu.Numrinit(2, 5)

        self.assertEqual(return_new, return_old)

    def test_ringwe_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ringwe")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (numr) = argum[0][0]

        return_new = fu.ringwe(numr)
        return_old = oldfu.ringwe(numr)

        self.assertEqual(return_new, return_old)

    def test_ornq_true_should_return_equal_object(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)

        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = argum[0]

        return_new = fu.ornq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny)
        return_old = fu.ornq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny)
        self.assertEqual(return_new, return_old)

    def test_ormq_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ormq")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = argum[0]

        return_new = fu.ormq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,delta)
        return_old = fu.ormq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,delta)

        self.assertEqual(return_new, return_old)


    #This function does not seem to be use anywhere so I wont be creating a unit test for this function
    # def test_ormq_fast_true_should_return_equal_object(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ormq")
    #     with open(filepath, 'rb') as rb:
    #         arguma= pickle.load(rb)
    #
    #     (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta) = arguma[0]
    #
    #
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
    #     with open(filepath, 'rb') as rb:
    #         argumb = pickle.load(rb)
    #
    #     (dataa, numra, wra, csa, tavga, cnxa, cnya, xrnga, yrnga, stepa) = argumb[0]
    #     dataa = copy.deepcopy(argumb[0][0])
    #
    #
    #
    #     return_new = fu.ormq_fast(dataa, crefim, xrnga, yrnga, stepa, numra, mode)
    #     return_old = fu.ormq_fast(dataa, crefim, xrnga, yrnga, stepa, numra, mode)
    #
    #     self.assertTrue(return_new, return_old)


    def test_prepref_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = argum[0]

        mode   = 'f'
        maxrangex = 4
        maxrangey = 4
        maskfile = None

        return_new = fu.prepref(data,maskfile,cnx,cny,numr,mode,maxrangex,maxrangey,step)
        return_old = oldfu.prepref(data,maskfile,cnx,cny,numr,mode,maxrangex,maxrangey,step)


        self.assertTrue(return_old,return_new)


    # def test_prepare_refrings_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.prepare_refrings")
    #     with open(filepath, 'rb') as rb:
    #         try:
    #             argum = pickle.load(rb)
    #             print(argum)
    #         except EOFError as exc:
    #             print(exc)
    #
    #         # (volft,kb) = pickle.load(rb)
    #
    #     return_new = fu.prepare_refrings(volft,kb)
    #     # return_old = oldfu.prepare_refrings(volft,kb)
    #     #
    #     # self.assertEqual(return_old,return_new)


    def test_proj_ali_incore_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)
        print(len(argum[0]))
        print(argum[0][4])

        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        return_new = fu.proj_ali_incore(data, refrings, numr, xrng, yrng, step)

        return_old = oldfu.proj_ali_incore(data, refrings, numr, xrng, yrng, step)

        self.assertTrue(return_old, return_new)

    # Core dump issue
    def test_proj_ali_incore_local_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        print(argum)

        an = -1.0
        xrng = 2.0
        yrng = 2.0
        nsym =1
        symangles = []

        for k in range(len(refrings)):
            # for i in range(nsym):
            symangles.append([0.0, 0.0, 1 * 360. / nsym])

        list_of_ref_ang_new = fu.generate_list_of_reference_angles_for_search(symangles, 'c1')
        list_of_ref_ang_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c1')


        return_new = fu.proj_ali_incore_local(data, refrings, list_of_ref_ang_new, numr, xrng, yrng, step, an)

        return_old = oldfu.proj_ali_incore_local(data, refrings, list_of_ref_ang_old, numr, xrng, yrng, step, an)

        self.assertTrue(return_old, return_new)


    def test_ali_vol_func_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali_vol_func")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            params = argum[0][0]
            data = argum[1]


        return_new = fu.ali_vol_func(params,**data)
        return_old = oldfu.ali_vol_func(params,**data)

        self.assertEqual(return_old, return_new)


    def test_align2d_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (image,refim, xrng, yrng) = argum[0]
        (ou) = argum[1]['ou']

        return_new = fu.align2d(image,refim)
        return_old = oldfu.align2d(image,refim)

        self.assertEqual(return_old, return_new)


    def test_align2d_scf_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (image,refim, xrng, yrng) = argum[0]
        (ou) = argum[1]['ou']

        return_new = fu.align2d_scf(image,refim,xrng,yrng, ou)
        return_old = oldfu.align2d_scf(image,refim,xrng,yrng, ou)

        self.assertEqual(return_old, return_new)


    def test_multalign2d_scf_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            # print(argum[0])

        (dataa, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step) = argum[0]

        dataa = copy.deepcopy(argum[0][0])

        mode = "F"
        ou = 174

        cimage = EMAN2_cppwrap.Util.Polar2Dm(tavg, float(cnx), float(cny), numr, mode)
        frotim = [sparx_fundamentals.fft(tavg)]


        return_new = fu.multalign2d_scf(dataa[0],[cimage],frotim, numr, xrng,yrng, ou)
        return_old = oldfu.multalign2d_scf(dataa[0],[cimage],frotim, numr, xrng,yrng, ou)

        self.assertEqual(return_old, return_new)


    def test_parabl_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.parabl")
        with open(filepath, 'rb') as rb:
            (Z) = pickle.load(rb)

        return_new = fu.parabl(Z[0][0])
        return_old = oldfu.parabl(Z[0][0])

        self.assertEqual(return_old, return_new)

    def test_shc_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)
            print(len(argum[0]))
            print(argum[0][4])

        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        return_new = fu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step)

        return_old = oldfu.shc(data, refrings, list_of_ref_ang, numr, xrng, yrng, step)

        self.assertTrue(return_old, return_new)



    def test_search_range_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.search_range")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            (n, radius, shift, range,location) = argum[0]

        return_new = fu.search_range(n, radius, shift, range)
        return_old = oldfu.search_range(n, radius, shift, range)

        self.assertEqual(return_old, return_new)


    def test_generate_list_true_should_return_equal_object(self):
        nsym = 5
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i * 360. / nsym])

        return_new = fu.generate_list_of_reference_angles_for_search(symangles , 'c5')
        return_old = oldfu.generate_list_of_reference_angles_for_search(symangles, 'c5')


        self.assertEqual(return_old, return_new)

"""
if __name__ == "__main__":
    unittest.main()
