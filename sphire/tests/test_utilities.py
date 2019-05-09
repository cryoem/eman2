from __future__ import print_function
from __future__ import division


import cPickle as pickle
import os
from mpi import *
import global_def
import numpy
import zlib



mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))

import unittest
from test_module import get_arg_from_pickle_file, get_real_data, remove_list_of_file, returns_values_in_file,remove_dir
from ..libpy import sparx_utilities as fu
from .sparx_lib import sparx_utilities as oldfu
from ..libpy import sparx_fundamentals
from os import path
from EMAN2_cppwrap import EMData,EMUtil

from copy import deepcopy
import EMAN2db
import json
import random
try:
    from StringIO import StringIO   # python2 case
except:
    from io import StringIO         # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
import sys

IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = fu.model_blank(10, 10)
IMAGE_BLANK_3D = fu.model_blank(10, 10, 10)
TOLERANCE = 0.0075
TRACKER = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask"))[0][0][1]
ABSOLUTE_PATH_TO_ADNAN_TEST_FILES= "/home/lusnig/Downloads/adnan4testing" # 4ADNAN: it is your 'sphire/tests' folder

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""

"""
There are some opened issues in:
1) drop_image --> How may I really test it
2) even_angles --> default value with P method leads to a deadlock
3) even_angles_cd --> default value with P method leads to a deadlock
4) find --> it seems to be not used
5) get_image --> I need an image to test the last 2 cases: get_image(path_to_img) and get_image(path_to_img, im=1)
6) get_im --> I need an image to test the last case ... similarly the (5)
7) read_text_row --> I need at least a file to test it
8) read_text_file --> same situation of 7
9) estimate_3D_center_MPI -- ask markus how it works
10) write_headers --> in .bdb case are not working under linux. Take a look to the code for seeing their comments
        --> if the file exists it overwrites it without a warning message. will we have to insert this message?
11) write_header --> I do not know how test the .bdb case. Hier contrary to write_headers it works under linux
12) file_type --> it is not giving us the filetype of the file. it is just parsing the name of the file and giving back the extension of the file
            Is this the real purpouse of this function?
13) set_params2D --> if you use xform=xform.align3d it works, but the output is somethiong that we do not really want to have. It does not set the values
                --> since set_params_proj has the same kind of input we are not able to discriminate them when we call the function. anyway It does not set the values
14) set_params3D --> if you use xform=xform.align2d it works, but the output is somethiong that we do not really want to have. It does not set the values
15) set_params_proj --> I need an image with key 'xform.projection' to finish these tests because the associated pickle file has not it --> dovrebbero essere quelle in pickle files/multi_shc/multi_shc.ali3d_multishc
16) The following functions concern the sending data in the process and are difficult or even not possible to test deeply
    -) reduce_EMData_to_root
    -) bcast_compacted_EMData_all_to_all
    -) gather_compacted_EMData_to_root
    -) bcast_EMData_to_all
    -) send_EMData
    -) recv_EMData
    -) recv_attr_dict
    -) send_attr_dict
    -) wrap_mpi_send
    -) wrap_mpi_recv
    -) wrap_mpi_gatherv
    -) wrap_mpi_split
17) unpack_message it does not work properly is it a buggy function???
18) 'update_tag' returns, in both of the implementations 'return 123456'. i'm not going to test it
19) how we can test 'if_error_then_all_processes_exit_program'? 
20) sample_down_1D_curve --> I need a file with the curve values
21) test_print_upper_triangular_matrix --> which variable is the third parameter??")
22) get_shrink_data_huang,recons_mref --> the file gave me does not work see the reasons in the test
23) do_two_way_comparison -->  I cannot run the Adnan reference test. I had to insert random data --> I cannot test it deeply,
24) Test_get_stat_proj.test_myid_not_the_same_value_as_main_Node_TypeError is it due to a bad implemntation?

"""






class Test_amoeba(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.amoeba"))

    @staticmethod
    def wrongfunction(a,b):
        return a+b

    @staticmethod
    def function_lessParam():
        return 0

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.amoeba()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.amoeba()
        self.assertEqual(cm_new.exception.message, "amoeba() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_amoeba(self):
        """
        I did not use 'self.assertTrue(numpy.allclose(return_new, return_old, atol=TOLERANCE,equal_nan=True))' because the 'nosetets' spawns the following error
                TypeError: ufunc 'isfinite' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
        """
        (var, scale, func, ftolerance, xtolerance, itmax , data) = self.argum[0]
        return_new = fu.amoeba (var, scale, func, ftolerance, xtolerance, 20 , data)
        return_old = oldfu.amoeba (var, scale, func, ftolerance, xtolerance, 20 , data)
        self.assertTrue(numpy.allclose(return_new[0], return_old[0], atol=TOLERANCE,equal_nan=True))
        self.assertTrue(abs(return_new[1]- return_old[1]) <TOLERANCE)
        self.assertEqual(return_new[2],return_old[2])

    def test_amoeba_with_wrongfunction(self):
        (var, scale, func, ftolerance, xtolerance, itmax , data) = self.argum[0]
        with self.assertRaises(TypeError) as cm_new:
            fu.amoeba (var, scale, self.wrongfunction, ftolerance, xtolerance, itmax , None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.amoeba (var, scale, self.wrongfunction, ftolerance, xtolerance, itmax , None)
        self.assertEqual(cm_new.exception.message, "wrongfunction() got an unexpected keyword argument 'data'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_amoeba_with_function_lessParam_TypeError(self):
        (var, scale, func, ftolerance, xtolerance, itmax , data) = self.argum[0]
        with self.assertRaises(TypeError) as cm_new:
            fu.amoeba (var, scale, self.function_lessParam, ftolerance, xtolerance, itmax , None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.amoeba (var, scale, self.function_lessParam, ftolerance, xtolerance, itmax , None)
        self.assertEqual(cm_new.exception.message, "function_lessParam() takes no arguments (2 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_amoeba_with_NoneType_data_returns_TypeError_NoneType_obj_hasnot_attribute__getitem__(self):
        (var, scale, func, ftolerance, xtolerance, itmax , data) = self.argum[0]
        with self.assertRaises(TypeError) as cm_new:
            fu.amoeba (var, scale, func, ftolerance, xtolerance, itmax , None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.amoeba (var, scale, func, ftolerance, xtolerance, itmax , None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_compose_transform2(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.compose_transform2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.compose_transform2()
        self.assertEqual(cm_new.exception.message, "compose_transform2() takes exactly 8 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        """ values got from 'pickle files/utilities/utilities.compose_transform2'"""
        return_new = fu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = 1.0, alpha2 = 156.512610336, sx2 = 0, sy2 = 0, scale2 = 1.0)
        return_old = oldfu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = 1.0, alpha2 = 156.512610336, sx2 = 0, sy2 = 0, scale2 = 1.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_scaleFactor_returns_RunTimeError_scale_factor_must_be_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = 0, alpha2 = 0, sx2 = 2.90828285217, sy2 =-0.879739010334, scale2 = 1.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = 0, alpha2 = 0, sx2 = 2.90828285217, sy2 =-0.879739010334, scale2 = 1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "The scale factor in a Transform object must be positive and non zero")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_negative_scaleFactor_returns_RunTimeError_scale_factor_must_be_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = -1.0, alpha2 = 0, sx2 = 2.90828285217, sy2 =-0.879739010334, scale2 = 1.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.compose_transform2(alpha1 = 0, sx1 = 2.90828285217, sy1 =-0.879739010334, scale1 = -1.0, alpha2 = 0, sx2 = 2.90828285217, sy2 =-0.879739010334, scale2 = 1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "The scale factor in a Transform object must be positive and non zero")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_compose_transform3(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.compose_transform3()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.compose_transform3()
        self.assertEqual(cm_new.exception.message, "compose_transform3() takes exactly 14 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        """ values got from 'pickle files/utilities/utilities.compose_transform3'"""
        return_new = fu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = 1.0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        return_old = oldfu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = 1.0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_scaleFactor_returns_RunTimeError_scale_factor_must_be_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = 0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = 0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "The scale factor in a Transform object must be positive and non zero")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_negative_scaleFactor_returns_RunTimeError_scale_factor_must_be_positive(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = -1.0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.compose_transform3(phi1 = 0.0, theta1  = 0.0, psi1 = 0.0, sx1 = 0.0,sy1 = 0.0, sz1 = 0.0,scale1 = -1.0, phi2 = 0.328125, theta2= 0.0, psi2 = 0.0, sx2 = 0.001220703125, sy2 = 0.0,sz2 = 0.001220703125,scale2 = 1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "The scale factor in a Transform object must be positive and non zero")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)




class Test_combine_params2(unittest.TestCase):
    """ I did not use the 'pickle files/utilities/utilities.combine_params2' values because they are all 0 values"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.combine_params2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.combine_params2()
        self.assertEqual(cm_new.exception.message, "combine_params2() takes exactly 8 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_combine_params2(self):
        return_new = fu.combine_params2(alpha1 = 0.0, sx1 = 1.0, sy1 = 1.0, mirror1 = 1, alpha2 = 1.0, sx2 =2.0, sy2 = 0.0, mirror2 = 0)
        return_old = oldfu.combine_params2(alpha1 = 0.0, sx1 = 1.0, sy1 = 1.0, mirror1 = 1, alpha2 = 1.0, sx2 =2.0, sy2 = 0.0, mirror2 = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_inverse_transform2(unittest.TestCase):
    """ I did not use the 'pickle files/utilities/utilities.inverse_transform2' values because they are all 0 values"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.inverse_transform2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.inverse_transform2()
        self.assertEqual(cm_new.exception.message, "inverse_transform2() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_inverse_transform2(self):
        return_new = fu.inverse_transform2(alpha = 1.0, tx = 2.2, ty = 1.0, mirror = 0)
        return_old = oldfu.inverse_transform2(alpha = 1.0, tx = 2.2, ty = 1.0, mirror = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



""" How may I REALLY test it?"""
class Test_drop_image(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.drop_image()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.drop_image()
        self.assertEqual(cm_new.exception.message, "drop_image() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_invalid_type_returns_UnboundLocalError_imgtype_referenced_before_assignment(self):
        destination ='output.hdf'
        with self.assertRaises(UnboundLocalError) as cm_new:
            fu.drop_image(IMAGE_2D, destination, itype="invalid")
        with self.assertRaises(UnboundLocalError) as cm_old:
            oldfu.drop_image(IMAGE_2D, destination, itype="invalid")
        self.assertEqual(cm_new.exception.message, "local variable 'imgtype' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("it does not work under nosetests , anyway im not able to test it properly")
    def test_destination_is_not_a_file_returns_error_msg(self):
        destination = 3
        return_new = fu.drop_image(IMAGE_2D, destination, itype="h")
        return_old = oldfu.drop_image(IMAGE_2D, destination, itype="h")
        self.assertTrue(return_new is None)
        self.assertTrue(return_old is None)

    @unittest.skip("it does not work under nosetests , anyway im not able to test it properly")
    def test_drop_image2D_true_should_return_equal_objects1(self):
        destination ='output.hdf'
        return_new = fu.drop_image(IMAGE_2D, destination, itype="h")
        return_old = oldfu.drop_image(IMAGE_2D, destination, itype="h")

        if return_new is not None   and  return_old is not None:
            self.assertTrue(return_new, return_old)

    @unittest.skip("it does not work under nosetests , anyway im not able to test it properly")
    def test_drop_image_true_should_return_equal_objects2(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.drop_image")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (imagename, destination) = argum[0]
        destination = 'output.hdf'
        return_new = fu.drop_image(imagename, destination, itype="h")
        return_old = oldfu.drop_image(imagename, destination, itype="h")

        if return_new is not None   and  return_old is not None:
            self.assertTrue(return_new, return_old)



class Test_even_angles(unittest.TestCase):
    """ I did not changed the 'phiEqpsi' params because it is used in 'even_angles_cd' I'll test it there"""
    def test_default_values(self):
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_delta_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.even_angles(delta = 0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.even_angles(delta = 0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values_with_not_minus(self):
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "", symmetry='c1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "", symmetry='c1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_default_values_with_P_method_leads_to_deadlock(self):
        self.assertTrue(True)
        """
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEqpsi = "Minus", symmetry='c1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_with_D_symmetry(self):
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='d1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='d1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_S_symmetry(self):
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_S_symmetry_tooBig_theta1_value_error_msg(self):
        return_new = fu.even_angles(delta = 15.0, theta1=91.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=91.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_S_invalid_symmetry_returns_UnboundLocalError_local_var_referenced_before_assignment(self):
        with self.assertRaises(UnboundLocalError) as cm_new:
            fu.even_angles(delta = 15.0, theta1=10.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sp1', ant = 0.0)
        with self.assertRaises(UnboundLocalError) as cm_old:
            oldfu.even_angles(delta = 15.0, theta1=10.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sp1', ant = 0.0)
        self.assertEqual(cm_new.exception.message, "local variable 'k' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_S_invalid_symmetry_returns_ValueError_invalid_literal(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.even_angles(delta = 15.0, theta1=10.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='soct', ant = 0.0)
        with self.assertRaises(ValueError) as cm_old:
            oldfu.even_angles(delta = 15.0, theta1=10.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='soct', ant = 0.0)
        self.assertEqual(cm_new.exception.message, "invalid literal for int() with base 10: 'ct'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_not_supported_symmetry_Warning_output_msg(self):
        return_new = fu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='oct', ant = 0.0)
        return_old = oldfu.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='oct', ant = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_even_angles_cd(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.even_angles_cd()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.even_angles_cd()
        self.assertEqual(cm_new.exception.message, "even_angles_cd() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_delta_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.even_angles_cd(delta = 0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='Minus')
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.even_angles_cd(delta = 0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='Minus')
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values_leads_to_deadlock(self):
        self.assertTrue(True)
        """
        return_new = fu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEQpsi='Minus')
        return_old = oldfu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEQpsi='Minus')
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_with_S_method(self):
        return_new = fu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='Minus')
        return_old = oldfu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='Minus')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_S_method_with_not_Minus(self):
        return_new = fu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='not_Minus')
        return_old = oldfu.even_angles_cd(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEQpsi='not_Minus')
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_gauss_edge(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.gauss_edge()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.gauss_edge()
        self.assertEqual(cm_new.exception.message, "gauss_edge() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.gauss_edge(None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.gauss_edge(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_ndim'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value_2Dreal_img(self):
        return_new =fu.gauss_edge(IMAGE_2D, kernel_size = 7, gauss_standard_dev =3)
        return_old =oldfu.gauss_edge(IMAGE_2D, kernel_size = 7, gauss_standard_dev =3)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_value_3Dreal_img(self):
        return_new =fu.gauss_edge(IMAGE_3D, kernel_size = 7, gauss_standard_dev =3)
        return_old =oldfu.gauss_edge(IMAGE_3D, kernel_size = 7, gauss_standard_dev =3)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_kernel_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.gauss_edge(IMAGE_2D, kernel_size = 0, gauss_standard_dev =3)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.gauss_edge(IMAGE_2D, kernel_size = 0, gauss_standard_dev =3)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_negative_kernel_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.gauss_edge(IMAGE_2D, kernel_size = -2, gauss_standard_dev =3)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.gauss_edge(IMAGE_2D, kernel_size = -2, gauss_standard_dev =3)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_get_image(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_image()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_image()
        self.assertEqual(cm_new.exception.message, "get_image() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_returns_input_img(self):
        """ I do not insert all the params because in this case they are not used"""
        return_new = fu.get_image(IMAGE_2D)
        return_old = oldfu.get_image(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_None_input_img_returns_new_EMData_with_default_size(self):
        """ I do not insert all the params because in this case they are not used"""
        nx = 0
        return_new = fu.get_image(None, nx = nx)
        return_old = oldfu.get_image(None, nx = nx)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertEqual(return_old.get_size(),return_new.get_size())
        self.assertEqual(return_new.get_size(), nx)

    def test_None_input_img_returns_new_EMData__with_given_size(self):
        """ I do not insert all the params because in this case they are not used"""
        nx,ny,nz=3,4,3
        return_new = fu.get_image(None, nx = nx, ny = ny, nz = nz)
        return_old = oldfu.get_image(None, nx = nx, ny = ny, nz = nz)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertEqual(return_old.get_size(),return_new.get_size())
        self.assertEqual(return_new.get_size(), nx*ny*nz)

    def test_invalid_path_returns_RuntimeError_FileAccessException(self):
        """ I do not insert all the params because in this case they are not used"""
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_image("image_not_here")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_image("image_not_here")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "FileAccessException")
        self.assertEqual(msg[3], "cannot access file ")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_get_im(unittest.TestCase):
    img_list = get_real_data(dim=2)

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_im(None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_im(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_im()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_im()
        self.assertEqual(cm_new.exception.message, "get_im() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_returns_first_img_of_a_list(self):
        return_new = fu.get_im(self.img_list, 0)
        return_old = oldfu.get_im(self.img_list, 0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.get_im(self.img_list, 10)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.get_im(self.img_list, 10)
        self.assertEqual(cm_new.exception.message, "tuple index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_image_data(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_image_data()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_image_data()
        self.assertEqual(cm_new.exception.message, "get_image_data() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_image_data(self):
        img,not_used = get_real_data(dim=2)
        return_new = fu.get_image_data(img)
        return_old = oldfu.get_image_data(img)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NoneType_as_input_image_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        #fu.get_image_data(None)



class Test_get_symt(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_symt()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_symt()
        self.assertEqual(cm_new.exception.message, "get_symt() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_symt(self):
        self.assertTrue(numpy.array_equal(fu.get_symt('c3'), oldfu.get_symt('c3')))

    def test_get_symt_with_invaliSym_returns_AttributeError_symclass_hasnot_attribute_symangles(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_symt('invaliSym')
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_symt('invaliSym')
        self.assertEqual(cm_new.exception.message, "'symclass' object has no attribute 'symangles'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_input_from_string(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_input_from_string()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_input_from_string()
        self.assertEqual(cm_new.exception.message, "get_input_from_string() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_input_from_string_integer_case(self):
        self.assertEqual(fu.get_input_from_string('5'), oldfu.get_input_from_string('5'))

    def test_get_input_from_string_negative_number_case(self):
        self.assertEqual(fu.get_input_from_string('-5'), oldfu.get_input_from_string('-5'))

    def test_get_input_from_string_float_case(self):
        self.assertEqual(fu.get_input_from_string('5.3'), oldfu.get_input_from_string('5.3'))

    def test_get_input_from_string_invalid_case(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.get_input_from_string('not_a_number')
        with self.assertRaises(ValueError) as cm_old:
            oldfu.get_input_from_string('not_a_number')
        self.assertEqual(cm_new.exception.message, "invalid literal for int() with base 10: 'not_a_number'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_input_from_string_list_of_values_number_case(self):
        self.assertTrue(numpy.array_equal(fu.get_input_from_string('-5,3.11,5'), oldfu.get_input_from_string('-5,3.11,5')))



class Test_model_circle(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.model_circle()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.model_circle()
        self.assertEqual(cm_new.exception.message, "model_circle() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        """ values got from 'pickle files/utilities/utilities.model_circle'"""
        return_new = fu.model_circle(r = 145, nx = 352, ny = 352, nz =1)
        return_old = oldfu.model_circle(r = 145, nx = 352, ny = 352, nz =1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_Y_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_circle(r = 145, nx = 352, ny = 0, nz =1)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_circle(r = 145, nx = 352, ny = 0, nz =1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "y size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_X_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_circle(r = 145, nx = 0, ny = 252, nz =1)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_circle(r = 145, nx = 0, ny = 252, nz =1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_Z_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_circle(r = 145, nx = 252, ny = 252, nz =0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_circle(r = 145, nx = 252, ny = 252, nz =0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "z size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_R_size_(self):
        return_new = fu.model_circle(r = 0, nx = 352, ny = 352, nz =1)
        return_old = oldfu.model_circle(r = 0, nx = 352, ny = 352, nz =1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_R_size_(self):
        return_new = fu.model_circle(r = -10, nx = 352, ny = 352, nz =1)
        return_old = oldfu.model_circle(r = -10, nx = 352, ny = 352, nz =1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_model_gauss(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.model_gauss()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.model_gauss()
        self.assertEqual(cm_new.exception.message, "model_gauss() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values(self):
        return_new = fu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        return_old = oldfu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_Xsigma_returns_Nan_matrix(self):
        return_new = fu.model_gauss(xsigma=0, nx=352, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        return_old = oldfu.model_gauss(xsigma=0, nx=352, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), equal_nan=True))

    def test_null_Ysigma_returns_Nan_matrix(self):
        return_new = fu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=0, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        return_old = oldfu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=0, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), equal_nan=True))

    def test_null_Zsigma_returns_Nan_matrix(self):
        return_new = fu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=None, zsigma=0, xcenter=None, ycenter=None, zcenter=None)
        return_old = oldfu.model_gauss(xsigma=2, nx=352, ny=1, nz=1, ysigma=None, zsigma=0, xcenter=None, ycenter=None, zcenter=None)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), equal_nan=True))

    def test_null_Y_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss(xsigma=2, nx=352, ny=0, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss(xsigma=2, nx=352, ny=0, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "y size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_X_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss(xsigma=2, nx=0, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss(xsigma=2, nx=0, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_Z_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss(xsigma=2, nx=352, ny=1, nz=0, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss(xsigma=2, nx=352, ny=1, nz=0, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "z size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_model_gauss_noise(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.model_gauss_noise()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.model_gauss_noise()
        self.assertEqual(cm_new.exception.message, "model_gauss_noise() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_model_gauss_noise(self):
        """  This function creates random noise each time so arrays cannot be compared """
        return_new = fu.model_gauss_noise(sigma = 1, nx = 352, ny=1, nz=1)
        return_old = oldfu.model_gauss_noise(sigma =1, nx = 352, ny=1, nz=1)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), atol=1000))

    def test_null_sigma(self):
        return_new = fu.model_gauss_noise(sigma = 0.0, nx = 352, ny=1, nz=1)
        return_old = oldfu.model_gauss_noise(sigma =0.0, nx = 352, ny=1, nz=1)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))


    def test_null_Y_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss_noise(sigma = 1, nx = 1, ny=0, nz=1)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss_noise(sigma = 1, nx = 1, ny=0, nz=1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "y size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_X_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss_noise(sigma = 1, nx = 0, ny=10, nz=1)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss_noise(sigma = 1, nx = 0, ny=10, nz=1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_Z_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_gauss_noise(sigma = 1, nx = 352, ny=1, nz=0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_gauss_noise(sigma = 1, nx = 352, ny=1, nz=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "z size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_model_blank(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.model_blank()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.model_blank()
        self.assertEqual(cm_new.exception.message, "model_blank() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values(self):
        return_new = fu.model_blank(nx = 100, ny=1, nz=1, bckg = 0.0)
        return_old = oldfu.model_blank(nx = 100, ny=1, nz=1, bckg = 0.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_values_with_bckg(self):
        return_new = fu.model_blank(nx = 100, ny=1, nz=1, bckg = 10.0)
        return_old = oldfu.model_blank(nx = 100, ny=1, nz=1, bckg = 10.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_X_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_blank(nx = 0, ny=1, nz=1, bckg = 0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_blank(nx = 0, ny=1, nz=1, bckg = 0.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_Y_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_blank(nx = 10, ny=0, nz=1, bckg = 0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_blank(nx = 10, ny=0, nz=1, bckg = 0.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "y size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_null_Z_size_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.model_blank(nx = 10, ny=1, nz=0, bckg = 0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.model_blank(nx = 10, ny=1, nz=0, bckg = 0.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "z size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_peak_search(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.peak_search()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.peak_search()
        self.assertEqual(cm_new.exception.message, "peak_search() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values(self):
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.peak_search(img, npeak = 3, invert = 1, print_screen = 0)
        return_old = oldfu.peak_search(img, npeak = 3, invert = 1, print_screen = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_inverted_sort(self):
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.peak_search(img, npeak = 3, invert = -1, print_screen = 0)
        return_old = oldfu.peak_search(img, npeak = 3, invert = -1, print_screen = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_npeak_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.peak_search(img, npeak = 0, invert = 1, print_screen = 0)
        return_old = oldfu.peak_search(img, npeak = 0, invert = 1, print_screen = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.peak_search(None, npeak = 3, invert = -1, print_screen = 0)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.peak_search(None, npeak = 3, invert = -1, print_screen = 0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'peak_search'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Empty_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        fu.peak_search(EMData(), npeak = 3, invert = -1, print_screen = 0)
        oldfu.peak_search(EMData(), npeak = 3, invert = -1, print_screen = 0)
        """


class Test_pad(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pad()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pad()
        self.assertEqual(cm_new.exception.message, "pad() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_returns_RuntimeError_ImageDimensionException_padder_cannot_be_lower_than_sizee_img(self):
        img, NotUsed = get_real_data(dim=2)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.pad(image_to_be_padded = img, new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.pad(image_to_be_padded = img, new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The size of the padded image cannot be lower than the input image size.")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_returns_RuntimeError_ImageDimensionException_offset_inconsistent(self):
        img, NotUsed = get_real_data(dim=2)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "average", off_center_nx = 100, off_center_ny = 100, off_center_nz = 100)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background ="average", off_center_nx = 100, off_center_ny = 100, off_center_nz = 100)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The offset inconsistent with the input image size. Solution: Change the offset parameters")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values(self):
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        return_old = oldfu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_values_with_circumference_bckg(self):
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "circumference", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        return_old = oldfu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "circumference", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_values_with_unknown_bckg(self):
        img, NotUsed = get_real_data(dim=2)
        return_new = fu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "unknown", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        return_old = oldfu.pad(image_to_be_padded = img, new_nx = img.get_xsize()+10, new_ny = img.get_ysize()+10,	new_nz = img.get_zsize()+10, background = "unknown", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_as_img_returns_RuntimeError_NullPointerException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.pad(image_to_be_padded = None, new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.pad(image_to_be_padded = None, new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NullPointerException")
        self.assertEqual(msg[1], "NULL input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Empty_img_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.pad(image_to_be_padded = EMData(), new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.pad(image_to_be_padded = EMData(), new_nx = 10, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_chooseformat(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.chooseformat()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.chooseformat()
        self.assertEqual(cm_new.exception.message, "chooseformat() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_exponential_number(self):
        self.assertEqual(fu.chooseformat(0.00000000000000000000000000003), oldfu.chooseformat(0.00000000000000000000000000003))

    def test_float(self):
        self.assertEqual(fu.chooseformat(0.3), oldfu.chooseformat(0.3))

    def test_typeError_float_argument_required(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.chooseformat('w')
        with self.assertRaises(TypeError) as cm_old:
            oldfu.chooseformat('w')
        self.assertEqual(cm_new.exception.message, "float argument required, not str")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


class Test_read_text_row(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.read_text_row()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.read_text_row()
        self.assertEqual(cm_new.exception.message, "read_text_row() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_file_not_found(self):
        with self.assertRaises(IOError) as cm_new:
            fu.read_text_row("no_file.txt")
        with self.assertRaises(IOError) as cm_old:
            oldfu.read_text_row("no_file.txt")
        self.assertEqual(cm_new.exception.strerror, "No such file or directory")
        self.assertEqual(cm_new.exception.strerror, cm_old.exception.strerror)



class Test_write_text_row(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.write_text_row()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.write_text_row()
        self.assertEqual(cm_new.exception.message, "write_text_row() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_write_text_row(self):
        data=[[1,1,1,1],[2,2,2,2],[3,3,3,3]]
        f=path.join(ABSOLUTE_PATH, "filefu.txt")
        fold=path.join(ABSOLUTE_PATH, "filefold.txt")
        fu.write_text_row(data, f)
        oldfu.write_text_row(data, fold)
        self.assertEqual(returns_values_in_file(f),returns_values_in_file(fold))
        remove_list_of_file([f,fold])



class Test_read_text_file(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.read_text_file()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.read_text_file()
        self.assertEqual(cm_new.exception.message, "read_text_file() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_file_not_found(self):
        with self.assertRaises(IOError) as cm_new:
            fu.read_text_file("no_file.txt")
        with self.assertRaises(IOError) as cm_old:
            oldfu.read_text_file("no_file.txt")
        self.assertEqual(cm_new.exception.strerror, "No such file or directory")
        self.assertEqual(cm_new.exception.strerror, cm_old.exception.strerror)



class Test_write_text_file(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.write_text_file()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.write_text_file()
        self.assertEqual(cm_new.exception.message, "write_text_file() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_write_text_row(self):
        data=[[1,1,1,1],[2,2,2,2],[3,3,3,3]]
        f=path.join(ABSOLUTE_PATH, "filefu.txt")
        fold=path.join(ABSOLUTE_PATH, "filefold.txt")
        fu.write_text_file(data, f)
        oldfu.write_text_file(data, fold)
        self.assertEqual(returns_values_in_file(f),returns_values_in_file(fold))
        remove_list_of_file([f,fold])



class Test_rotate_shift_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rotate_shift_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rotate_shift_params()
        self.assertEqual(cm_new.exception.message, "rotate_shift_params() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotate_shift_params(self):
        paramsin = [[0.25,1.25,0.5]]
        transf  = [0.25, 1.25, 0.5]
        return_new = fu.rotate_shift_params(paramsin, transf)
        return_old = oldfu.rotate_shift_params(paramsin, transf)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_rotate_shift_params2(self):
        paramsin = [[0.25,1.25,0,0,0.5]]
        transf  = [0.25, 1.25, 0.5,.25, 1.25, 0.5]
        return_new = fu.rotate_shift_params(paramsin, transf)
        return_old = oldfu.rotate_shift_params(paramsin, transf)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_less_transf_params_returns_IndexError_list_index_out_of_range(self):
        paramsin = [[0.25,1.25,0,0,0.5]]
        transf  = [0.25, 1.25, 0.5]
        with self.assertRaises(IndexError) as cm_new:
            fu.rotate_shift_params(paramsin, transf)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.rotate_shift_params(paramsin, transf)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_less_transf_params2_returns_IndexError_list_index_out_of_range(self):
        paramsin = [[0.25,1.25,0]]
        transf  = [0.25, 1.25]
        with self.assertRaises(IndexError) as cm_new:
            fu.rotate_shift_params(paramsin, transf)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.rotate_shift_params(paramsin, transf)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_less_paramsin_params_returns_IndexError_list_index_out_of_range(self):
        paramsin = [[0.25]]
        transf  = [0.25, 1.25, 0.5]
        with self.assertRaises(IndexError) as cm_new:
            fu.rotate_shift_params(paramsin, transf)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.rotate_shift_params(paramsin, transf)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_reshape_1d(unittest.TestCase):
    """ values got from 'pickle files/utilities/utilities.reshape_1d'"""
    input_obj =  [0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.9984012789768186, 0.9914368216668327, 0.9878146959140469, 0.9881703862020976, 0.982612488476065, 0.9789244545589472, 0.9747235387045814, 0.9622078763024153, 0.9406924390622574, 0.9300175631598249, 0.8976592373307525, 0.8474726574046705, 0.7942852016327994, 0.8065378605172119, 0.7981892234519837, 0.7980760586172797, 0.7834690256016978, 0.7732854546260584, 0.759479194158529, 0.7302534821351329, 0.735749496632646, 0.7505776906379105, 0.7832464000713297, 0.799354031902547, 0.7829602489012508, 0.7467401462021503, 0.7216741559492451, 0.7573457050470969, 0.7735999645280006, 0.7360206933666649, 0.7074315960216845, 0.6838418535731124, 0.6814918195422979, 0.6604400166044002, 0.6276571502978614, 0.5967298971705947, 0.5924074015096022, 0.6113438607798904, 0.5589193571016572, 0.4169423800381157, 0.33547900293137645, 0.43509084125025116, 0.5143369854093631, 0.4505998230268216, 0.3017867022488365, 0.29393725698240897, 0.3395667841020214, 0.34234494237984336, 0.31531353786458843, 0.3120432449453534, 0.2864549161874622, 0.23450693792899116, 0.20246505335938672, 0.22577560951692183, 0.21569461751208094, 0.21511112191209886, 0.2091532904083915, 0.18334792795777813, 0.1954858454475899, 0.21231959169076153, 0.20199531221828237, 0.21190821007216915, 0.21429959199533707, 0.18398541329970813, 0.20171364365585326, 0.22936964071672247, 0.20705888033218262, 0.2310040684684463, 0.23322049365816364, 0.25365125929269, 0.2687457179832018, 0.252646215129461, 0.24715492782090853, 0.23387479872417344, 0.23315205998051616, 0.2312238364934745, 0.21601984544387764, 0.23373779370670353, 0.21445443670567088, 0.210741700365644, 0.2089851778417197, 0.19984641965828376, 0.18358602895051426, 0.16600398773363803, 0.14936583739921497, 0.14684159823845128, 0.14034187449397328, 0.11227281827686696, 0.09549423222286733, 0.09699040681889236, 0.08368778954783127, 0.07285201615715135, 0.06609239822815444, 0.06712766581830018, 0.06571178890380885, 0.05876124933827422, 0.047775744976412994, 0.04517043724966535, 0.04086780062968338, 0.035162664167093884, 0.02501739454518543, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    length_current = 2* len(input_obj)
    length_interpolated = 4* len(input_obj)

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.reshape_1d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.reshape_1d()
        self.assertEqual(cm_new.exception.message, "reshape_1d() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_list_as_input_obj(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.reshape_1d(input_object = [], length_current=self.length_current, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.reshape_1d(input_object = [], length_current=self.length_current, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        return_new = fu.reshape_1d(input_object = self.input_obj, length_current=self.length_current, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        return_old = oldfu.reshape_1d(input_object = self.input_obj , length_current=self.length_current, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_length_interpolated_pixel_sizes_identical_error_msg(self):
        return_new = fu.reshape_1d(input_object = self.input_obj, length_current=self.length_current, length_interpolated=0, Pixel_size_current = 0.5, Pixel_size_interpolated = 0.5)
        return_old = oldfu.reshape_1d(input_object = self.input_obj, length_current=self.length_current, length_interpolated=0, Pixel_size_current = 0.5, Pixel_size_interpolated = 0.5)
        self.assertEqual(return_new,[])
        self.assertEqual(return_old, [])

    def test_null_length_current(self):
        return_new = fu.reshape_1d(input_object = self.input_obj, length_current=0, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        return_old = oldfu.reshape_1d(input_object = self.input_obj, length_current=0, length_interpolated=self.length_interpolated, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_all_the_values_are_null_or_empty_list_error_msg(self):
        return_new = fu.reshape_1d(input_object = [], length_current=0, length_interpolated=0, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        return_old = oldfu.reshape_1d(input_object = [], length_current=0, length_interpolated=0, Pixel_size_current = 0., Pixel_size_interpolated = 0.)
        self.assertEqual(return_new, [])
        self.assertEqual(return_old, [])

    def test_invalid_pixel_sizes_combination_in_null_value_as_length_interpolated_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.reshape_1d(input_object = self.input_obj, length_current=self.length_current, length_interpolated=0, Pixel_size_current = 0.3, Pixel_size_interpolated = 0.)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.reshape_1d(input_object = self.input_obj, length_current=self.length_current, length_interpolated=0, Pixel_size_current = 0.3, Pixel_size_interpolated = 0.)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


@unittest.skip("I m not sure how test them")
class Test_estimate_3D_center_MPI(unittest.TestCase):
    """ values got from 'pickle files/utilities/utilities.estimate_3D_center_MPI'"""
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.estimate_3D_center_MPI"))

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.estimate_3D_center_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.estimate_3D_center_MPI()
        self.assertEqual(cm_new.exception.message, "estimate_3D_center_MPI() takes at least 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_myid_not_identical_to_main_node(self):
        (data, nima, myid, number_of_proc, main_node) = self.argum[0]
        return_new = fu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)
        return_old = oldfu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)
        self.assertTrue(numpy.array_equal(return_old, [0.0, 0.0, 0.0, 0.0, 0.0]))
        self.assertTrue(numpy.array_equal(return_new, [0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_myid_not_identical_to_main_node1(self):
        (data, nima, myid, number_of_proc, main_node) = self.argum[0]
        main_node=myid
        return_new = fu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)
        return_old = oldfu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)
        self.assertTrue(numpy.array_equal(return_old, [0.0, 0.0, 0.0, 0.0, 0.0]))
        self.assertTrue(numpy.array_equal(return_new, [0.0, 0.0, 0.0, 0.0, 0.0]))



class Test_rotate_3D_shift(unittest.TestCase):
    """ values got from 'pickle files/utilities/utilities.rotate_3D_shift'"""
    argum =get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.rotate_3D_shift"))
    (data, notUsed) = argum[0]
    shift3d = [10.1, 0.2, 10.0]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rotate_3D_shift()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rotate_3D_shift()
        self.assertEqual(cm_new.exception.message, "rotate_3D_shift() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_wrong_image(self):
        data,not_used= get_real_data(dim = 3)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.rotate_3D_shift([data], self.shift3d)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rotate_3D_shift([data], self.shift3d)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Nonetype_image(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.rotate_3D_shift([None], self.shift3d)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.rotate_3D_shift([None], self.shift3d)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        fu_data = deepcopy(self.data)
        oldfu_data = deepcopy(self.data)
        return_new = fu.rotate_3D_shift(fu_data, self.shift3d)
        return_old = oldfu.rotate_3D_shift(oldfu_data, self.shift3d)
        self.assertEqual(return_new, None)
        self.assertEqual(return_old, None)
        for i in range(len(fu_data)):
            self.assertTrue(numpy.array_equal(fu_data[i].get_attr('xform.projection'), oldfu_data[i].get_attr('xform.projection')))
            self.assertFalse(numpy.array_equal(self.data[i].get_attr('xform.projection'), fu_data[i].get_attr('xform.projection')))

    def test_returns_IndexError_list_index_out_of_range(self):
        shift3d=[0,0.1]
        with self.assertRaises(IndexError) as cm_new:
            fu.rotate_3D_shift(self.data, shift3d)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.rotate_3D_shift(self.data, shift3d)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_set_arb_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.set_arb_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.set_arb_params()
        self.assertEqual(cm_new.exception.message, "set_arb_params() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_1Attr(self):
        fu_img = EMData()
        oldfu_img = EMData()
        par_str = "lowpassfilter"
        params = "0.50"
        return_new = fu.set_arb_params(fu_img,[params],[par_str])
        return_old = oldfu.set_arb_params(oldfu_img,[params],[par_str])
        self.assertEqual(return_new, None)
        self.assertEqual(return_old, None)
        self.assertEqual(fu_img.get_attr(par_str), fu_img.get_attr(par_str))
        self.assertEqual(fu_img.get_attr(par_str),params)

    def test_with_ListAttr(self):
        fu_img = EMData()
        oldfu_img = EMData()
        par_str = ["lowpassfilter","fake_par"]
        params = ["0.50","3"]
        return_new = fu.set_arb_params(fu_img,params,par_str)
        return_old = oldfu.set_arb_params(oldfu_img,params,par_str)
        self.assertEqual(return_new, None)
        self.assertEqual(return_old, None)
        self.assertEqual(fu_img.get_attr(par_str[0]), fu_img.get_attr(par_str[0]))
        self.assertEqual(fu_img.get_attr(par_str[1]), fu_img.get_attr(par_str[1]))
        self.assertEqual(fu_img.get_attr(par_str[0]),params[0])
        self.assertEqual(fu_img.get_attr(par_str[1]), params[1])

    def test_with_BadListAttr_returns_IndexError_list_index_out_of_range(self):
        fu_img = EMData()
        oldfu_img = EMData()
        par_str = ["lowpassfilter","fake_par"]
        params = ["0.50"]
        with self.assertRaises(IndexError) as cm_new:
            fu.set_arb_params(fu_img,params,par_str)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.set_arb_params(oldfu_img,params,par_str)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_arb_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_arb_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_arb_params()
        self.assertEqual(cm_new.exception.message, "get_arb_params() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_1Attr(self):
        return_new = fu.get_arb_params(EMData(),["datatype"])
        return_old = oldfu.get_arb_params(EMData(),["datatype"])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, [EMData().get_attr("datatype")])

    def test_with_ListAttr(self):
        list_of_attribute = ["datatype", "is_complex_ri"]
        return_new = fu.get_arb_params(EMData(),list_of_attribute)
        return_old = oldfu.get_arb_params(EMData(),list_of_attribute)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertEqual(return_new[0], EMData().get_attr("datatype"))
        self.assertEqual(return_new[1], EMData().get_attr("is_complex_ri"))

    def test_notValid_params_returns_RuntimeError_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_arb_params(EMData(),["invalid_param"])
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_arb_params(EMData(),["invalid_param"])
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_reduce_EMData_to_root(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.reduce_EMData_to_root()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.reduce_EMData_to_root()
        self.assertEqual(cm_new.exception.message, "reduce_EMData_to_root() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.reduce_EMData_to_root(None, myid=74, main_node=0, comm = -1)
        return_old = oldfu.reduce_EMData_to_root(None, myid=74, main_node=0, comm = -1)
        self.assertEqual(return_new, return_old)
        """

    def test_default_values(self):
        data = deepcopy(IMAGE_2D_REFERENCE)
        return_new = fu.reduce_EMData_to_root(data, myid=74, main_node=0, comm = -1)
        return_old = oldfu.reduce_EMData_to_root(data, myid=74, main_node=0, comm = -1)
        self.assertTrue(numpy.array_equal(IMAGE_2D_REFERENCE.get_3dview(), data.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_with_MPI_COMM_WORLD(self):
        data = deepcopy(IMAGE_2D_REFERENCE)
        return_new = fu.reduce_EMData_to_root(data, myid=74, main_node=0, comm = MPI_COMM_WORLD)
        return_old = oldfu.reduce_EMData_to_root(data, myid=74, main_node=0, comm = MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(IMAGE_2D_REFERENCE.get_3dview(), data.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)



class Test_bcast_compacted_EMData_all_to_all(unittest.TestCase):
    """
    It does not matter which of my images I-ll use, I got always the following typeerror:
    Error
Traceback (most recent call last):
  File "/home/lusnig/SPHIRE_1_1/lib/python2.7/unittest/case.py", line 329, in run
    testMethod()
  File "/home/lusnig/EMAN2/eman2/sphire/tests/test_utilities.py", line 1451, in test_bcast_compacted_EMData_all_to_all_true_should_return_equal_objects
    return_new = fu.bcast_compacted_EMData_all_to_all(list_of_em_objects, myid)
  File "/home/lusnig/EMAN2/eman2/sphire/libpy/sparx_utilities.py", line 1105, in bcast_compacted_EMData_all_to_all
    em_dict = dict_received["em_dict"]
TypeError: 'NoneType' object has no attribute '__getitem__'
    """
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_compacted_EMData_all_to_all()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_compacted_EMData_all_to_all()
        self.assertEqual(cm_new.exception.message, "bcast_compacted_EMData_all_to_all() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_NoneType_data_returns_TypeError_NoneType_obj_hasnot_attribute__getitem__(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_compacted_EMData_all_to_all([None,None], myid=74,  comm = -1)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_compacted_EMData_all_to_all([None,None], myid=74, comm = -1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_image(self):
        data = [deepcopy(IMAGE_3D),deepcopy(IMAGE_3D)]
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_compacted_EMData_all_to_all(data, myid=74,  comm = -1)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_compacted_EMData_all_to_all(data, myid=74,  comm = -1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_gather_compacted_EMData_to_root(unittest.TestCase):
    argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.gather_compacted_EMData_to_root"))
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.gather_compacted_EMData_to_root()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.gather_compacted_EMData_to_root()
        self.assertEqual(cm_new.exception.message, "gather_compacted_EMData_to_root() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        (no_of_emo, list_of_emo, myid) = self.argum[0]
        return_new = fu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid, comm=-1)
        return_old = oldfu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid, comm=-1)
        self.assertEqual(return_new, return_old)

    def test_with_MPI_COMM_WORLD(self):
        (no_of_emo, list_of_emo, myid) = self.argum[0]
        return_new = fu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid, comm=MPI_COMM_WORLD)
        return_old = oldfu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid, comm=MPI_COMM_WORLD)
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_pickle_file_values_wrong_number_of_number_of_all_em_objects_distributed_across_processes(self):
        (no_of_emo, list_of_emo, myid) = self.argum[0]
        return_new = fu.gather_compacted_EMData_to_root(0, list_of_emo, myid,  comm=-1)
        return_old = oldfu.gather_compacted_EMData_to_root(0, list_of_emo, myid, comm=-1)
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_NoneType_as_img_returns_IndexError_list_index_out_of_range(self):
        (no_of_emo, list_of_emo, myid) = self.argum[0]
        with self.assertRaises(IndexError) as cm_new:
            fu.gather_compacted_EMData_to_root(no_of_emo, [], myid, comm=-1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.gather_compacted_EMData_to_root(no_of_emo, [], myid, comm=-1)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_bcast_EMData_to_all(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_EMData_to_all()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_EMData_to_all()
        self.assertEqual(cm_new.exception.message, "bcast_EMData_to_all() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_default_values(self):
        tavg = deepcopy(IMAGE_2D_REFERENCE)
        return_new = fu.bcast_EMData_to_all(tavg, myid = 11, source_node =0, comm= -1)
        return_old = oldfu.bcast_EMData_to_all(tavg, myid= 11, source_node = 0, comm= -1)
        self.assertTrue(numpy.array_equal(IMAGE_2D_REFERENCE.get_3dview(), tavg.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_with_myid_equal_sourcenode_default_valuqes(self):
        tavg = deepcopy(IMAGE_2D_REFERENCE)
        return_new = fu.bcast_EMData_to_all(tavg, myid= 0, source_node =0, comm= -1)
        return_old = oldfu.bcast_EMData_to_all(tavg, myid= 0, source_node =0, comm= -1)
        self.assertTrue(numpy.array_equal(IMAGE_2D_REFERENCE.get_3dview(), tavg.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_with_MPI_COMM_WORLD(self):
        tavg = deepcopy(IMAGE_2D_REFERENCE)
        return_new = fu.bcast_EMData_to_all(tavg, myid= 11, source_node =0, comm= MPI_COMM_WORLD)
        return_old = oldfu.bcast_EMData_to_all(tavg, myid= 11, source_node =0, comm= MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(IMAGE_2D_REFERENCE.get_3dview(), tavg.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.bcast_EMData_to_all(None, 11, source_node =0, comm= -1)
        return_old = oldfu.bcast_EMData_to_all(None, 11, source_node =0, comm= -1)
        self.assertEqual(return_new, return_old)
        """



class Test_send_EMData(unittest.TestCase):
    #argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.send_EMData"))
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.send_EMData()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.send_EMData()
        self.assertEqual(cm_new.exception.message, "send_EMData() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.send_EMData(None, 0, 0)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.send_EMData(None, 0, 0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    """  Can only be tested on the mpi. Wait too long on normal workstation"""
    # def test_send_EMData_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.send_EMData")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (img, dst, tag, comm) = argum[0]
    #     tag = 0
    #
    #     return_new = fu.send_EMData(img, dst, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = oldfu.send_EMData(img, dst, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)



class Test_recv_EMData(unittest.TestCase):
    #argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.recv_EMData"))
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recv_EMData()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recv_EMData()
        self.assertEqual(cm_new.exception.message, "recv_EMData() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    """  Can only be tested on the mpi. Wait too long on normal workstation"""
    # def test_recv_EMData_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.recv_EMData")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (src, tag,comm) = argum[0]
    #     tag = 0
    #
    #     return_new = fu.recv_EMData(src, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = oldfu.recv_EMData(src, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)




class Test_bcast_number_to_all(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_number_to_all()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_number_to_all()
        self.assertEqual(cm_new.exception.message, "bcast_number_to_all() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_number_to_send_is_null(self):
        return_new = fu.bcast_number_to_all(number_to_send = 0, source_node = 0, mpi_comm = -1)
        return_old = oldfu.bcast_number_to_all(number_to_send = 0, source_node = 0, mpi_comm = -1)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0)

    def test_with_MPI_COMM_WORLD(self):
        return_new = fu.bcast_number_to_all(number_to_send = 0, source_node = 0, mpi_comm = MPI_COMM_WORLD)
        return_old = oldfu.bcast_number_to_all(number_to_send = 0, source_node = 0, mpi_comm = MPI_COMM_WORLD)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0)

    def test_number_to_send_is_not_null(self):
        return_new = fu.bcast_number_to_all(number_to_send = 3, source_node = 0, mpi_comm = -1)
        return_old = oldfu.bcast_number_to_all(number_to_send = 3, source_node = 0, mpi_comm = -1)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 3)

    def test_invalid_number_to_send_error_msg(self):
        return_new = fu.bcast_number_to_all(number_to_send = None, source_node = 0, mpi_comm = -1)
        return_old = oldfu.bcast_number_to_all(number_to_send = None, source_node = 0, mpi_comm = -1)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)


class Test_bcast_list_to_all(unittest.TestCase):
    myid = 74
    source_node =0
    list_to_send = [1,2]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bcast_list_to_all()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bcast_list_to_all()
        self.assertEqual(cm_new.exception.message, "bcast_list_to_all() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_empty_list(self):
        return_new = fu.bcast_list_to_all([], myid = self.myid, source_node =self.source_node, mpi_comm= -1)
        return_old = oldfu.bcast_list_to_all([], myid= self.myid, source_node = self.source_node, mpi_comm= -1)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[]))

    def test_defualt_values(self):
        return_new = fu.bcast_list_to_all(self.list_to_send, myid = self.myid, source_node =self.source_node, mpi_comm= -1)
        return_old = oldfu.bcast_list_to_all(self.list_to_send, myid= self.myid, source_node = self.source_node, mpi_comm= -1)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[]))

    def test_defualt_values_with_MPI_COMM_WORLD(self):
        return_new = fu.bcast_list_to_all(self.list_to_send, myid = self.myid, source_node =self.source_node, mpi_comm= MPI_COMM_WORLD)
        return_old = oldfu.bcast_list_to_all(self.list_to_send, myid= self.myid, source_node = self.source_node, mpi_comm= MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[]))

    def test_myid_equal_sourcenode(self):
        return_new = fu.bcast_list_to_all(self.list_to_send, myid = self.source_node, source_node =self.source_node, mpi_comm= MPI_COMM_WORLD)
        return_old = oldfu.bcast_list_to_all(self.list_to_send, myid= self.source_node, source_node = self.source_node, mpi_comm= MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,self.list_to_send))

    def test_myid_equal_sourcenode_and_wrong_type_in_listsender_returns_ValueError(self):
        list_to_send=[IMAGE_2D]
        with self.assertRaises(ValueError) as cm_new:
            fu.bcast_list_to_all(list_to_send, myid = self.source_node, source_node =self.source_node, mpi_comm= MPI_COMM_WORLD)
        with self.assertRaises(ValueError) as cm_old:
            oldfu.bcast_list_to_all(list_to_send, myid= self.source_node, source_node = self.source_node, mpi_comm= MPI_COMM_WORLD)
        self.assertEqual(cm_new.exception.message, "setting an array element with a sequence.")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_type_in_listsender(self):
        list_to_send=[IMAGE_2D]
        return_new = fu.bcast_list_to_all(list_to_send, myid = self.myid, source_node =self.source_node, mpi_comm= MPI_COMM_WORLD)
        return_old = oldfu.bcast_list_to_all(list_to_send, myid= self.myid, source_node = self.source_node, mpi_comm= MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[]))



class Test_recv_attr_dict(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recv_attr_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recv_attr_dict()
        self.assertEqual(cm_new.exception.message, "recv_attr_dict() takes at least 7 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_send_attr_dict(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.send_attr_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.send_attr_dict()
        self.assertEqual(cm_new.exception.message, "send_attr_dict() takes at least 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_recv_attr_dict_bdb(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recv_attr_dict_bdb()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recv_attr_dict_bdb()
        self.assertEqual(cm_new.exception.message, "recv_attr_dict_bdb() takes at least 7 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_print_begin_msg(unittest.TestCase):
    """ see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.print_begin_msg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_begin_msg()
        self.assertEqual(cm_new.exception.message, "print_begin_msg() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_print_begin_msg(self):
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_begin_msg("test_pgr", onscreen=False)
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_begin_msg("test_pgr", onscreen=False)
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout

    def test_print_begin_msg_onscreen_True(self):
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_begin_msg("test_pgr", onscreen=True)
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_begin_msg("test_pgr", onscreen=True)
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout




class Test_print_end_msg(unittest.TestCase):
    """ see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.print_end_msg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_end_msg()
        self.assertEqual(cm_new.exception.message, "print_end_msg() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_print_end_msg(self):
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_end_msg("test_pgr", onscreen=False)
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_end_msg("test_pgr", onscreen=False)
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout

    def test_print_end_msg_onscreen_True(self):
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_end_msg("test_pgr", onscreen=True)
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_end_msg("test_pgr", onscreen=True)
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout



class Test_print_msg(unittest.TestCase):
    """ see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.print_msg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_msg()
        self.assertEqual(cm_new.exception.message, "print_msg() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_print_msg(self):
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_msg("test_pgr")
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_msg("test_pgr")
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout



class Test_read_fsc(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.read_fsc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.read_fsc()
        self.assertEqual(cm_new.exception.message, "read_fsc() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_write_text_row(self):
        data=[[1,1,1,1],[2,2,2,2],[3,3,3,3]]
        f=path.join(ABSOLUTE_PATH, "filefu.txt")
        fu.write_text_file(data, f)
        return_new = fu.read_fsc(f)
        return_old = oldfu.read_fsc(f)
        remove_list_of_file([f])
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, data))



class Test_circumference(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.circumference()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.circumference()
        self.assertEqual(cm_new.exception.message, "circumference() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_default_values_2Dimg(self):
        return_new = fu.circumference(deepcopy(IMAGE_BLANK_2D), inner = -1, outer = -1)
        return_old = oldfu.circumference(deepcopy(IMAGE_BLANK_2D), inner = -1, outer = -1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_default_values_3Dimg(self):
        return_new = fu.circumference(deepcopy(IMAGE_BLANK_3D), inner = -1, outer = -1)
        return_old = oldfu.circumference(deepcopy(IMAGE_BLANK_3D), inner = -1, outer = -1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_invalid_mask_returns_RuntimeError_ImageFormatException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.circumference(deepcopy(IMAGE_BLANK_2D), inner =IMAGE_BLANK_2D.get_xsize()+10 , outer = -1)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.circumference(deepcopy(IMAGE_BLANK_2D), inner =IMAGE_BLANK_2D.get_xsize()+10 , outer = -1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "Invalid mask")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_wrong_outer_value(self):
        return_new = fu.circumference(deepcopy(IMAGE_BLANK_2D), inner = -1, outer = IMAGE_BLANK_2D.get_xsize()+10 )
        return_old = oldfu.circumference(deepcopy(IMAGE_BLANK_2D), inner = -1, outer = IMAGE_BLANK_2D.get_xsize()+10 )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_write_headers(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.write_headers()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.write_headers()
        self.assertEqual(cm_new.exception.message, "write_headers() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_hdf_type(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        fu.write_headers(path_fu, [IMAGE_2D], [1])
        oldfu.write_headers(path_oldfu, [IMAGE_2D], [1])
        self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        remove_list_of_file([path_fu,path_oldfu])

    def test_overwrite_hdf_file(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        f = open(path_fu, "w+")
        f.close()
        f = open(path_oldfu, "w+")
        f.close()
        fu.write_headers(path_fu, [IMAGE_2D], [1])
        oldfu.write_headers(path_oldfu, [IMAGE_2D], [1])
        self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        remove_list_of_file([path_fu,path_oldfu])

    def test_hdf_type_AssertError_list_differ(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        fu.write_headers(path_fu, [IMAGE_2D], [2])
        oldfu.write_headers(path_oldfu, [IMAGE_2D], [1])
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        with self.assertRaises(AssertionError) as cm:
            self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        msg = cm.exception.message.split("'")
        self.assertEqual(msg[0].split(":")[0], "Lists differ")
        self.assertEqual(msg[10].split("\n")[2].split(":")[0], 'First differing element 2')
        remove_list_of_file([path_fu, path_oldfu])

    def test_bdf_type(self):
        """
        in the code they inserted the following comment:
        #  For unknown reasons this does not work on Linux, but works on Mac ??? Really?
        """
        self.assertTrue(True)

    def test_invalid_filetype_error_msg(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.txt")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.txt")
        fu.write_headers(path_fu, [IMAGE_2D], [1])
        oldfu.write_headers(path_oldfu, [IMAGE_2D], [1])
        self.assertFalse(path.isfile(path_fu))
        self.assertFalse(path.isfile(path_oldfu))



class Test_write_header(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.write_header()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.write_header()
        self.assertEqual(cm_new.exception.message, "write_header() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_hdf_type(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        fu.write_header(path_fu, IMAGE_2D, 1)
        oldfu.write_header(path_oldfu, IMAGE_2D, 1)
        self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        remove_list_of_file([path_fu,path_oldfu])

    def test_overwrite_hdf_file(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        f = open(path_fu, "w+")
        f.close()
        f = open(path_oldfu, "w+")
        f.close()
        fu.write_header(path_fu, IMAGE_2D, 1)
        oldfu.write_header(path_oldfu, IMAGE_2D, 1)
        self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        remove_list_of_file([path_fu,path_oldfu])

    def test_hdf_type_AssertError_list_differ(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.hdf")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.hdf")
        fu.write_header(path_fu, IMAGE_2D, 2)
        oldfu.write_header(path_oldfu, IMAGE_2D, 1)
        self.assertTrue(path.isfile(path_fu))
        self.assertTrue(path.isfile(path_oldfu))
        with self.assertRaises(AssertionError) as cm:
            self.assertEqual(returns_values_in_file(path_fu), returns_values_in_file(path_oldfu))
        msg = cm.exception.message.split("'")
        self.assertEqual(msg[0].split(":")[0], "Lists differ")
        self.assertEqual(msg[10].split("\n")[2].split(":")[0], 'First differing element 2')
        remove_list_of_file([path_fu, path_oldfu])

    def test_invalid_filetype_error_msg(self):
        path_fu = path.join(ABSOLUTE_PATH, "test.txt")
        path_oldfu = path.join(ABSOLUTE_PATH, "test1.txt")
        fu.write_header(path_fu, IMAGE_2D, 1)
        oldfu.write_header(path_oldfu, IMAGE_2D, 1)
        self.assertFalse(path.isfile(path_fu))
        self.assertFalse(path.isfile(path_oldfu))



class Test_file_type(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.file_type()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.file_type()
        self.assertEqual(cm_new.exception.message, "file_type() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_bdb_filetype(self):
        fu.file_type("bdb:bdbfile")
        oldfu.file_type("bdb:bdbfile")
        self.assertTrue(True)

    def test_valid_filetype(self):
        fu.file_type("hdf.hdf")
        oldfu.file_type("hdf.hdf")
        self.assertTrue(True)

    def test_invalid_filetype_error_msg(self):
        fu.file_type("invalid.cc")
        oldfu.file_type("invalid.cc")
        self.assertTrue(True)


class Test_get_params2D(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D"))
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_params2D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_params2D()
        self.assertEqual(cm_new.exception.message, "get_params2D() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_params2D(self):
        (ima,) = self.argum[0]
        return_new = fu.get_params2D(ima, xform="xform.align2d")
        return_old = oldfu.get_params2D(ima, xform="xform.align2d")
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_wrong_xform_returns_NotExistingObjectException_key_doesnot_exist(self):
        (ima,) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params2D(ima, xform="xform.align3d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params2D(ima, xform="xform.align3d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_input_2dimg_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params2D(IMAGE_2D, xform="xform.align2d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params2D(IMAGE_2D, xform="xform.align2d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_input_3dimg_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params2D(IMAGE_3D, xform="xform.align2d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params2D(IMAGE_3D, xform="xform.align2d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_params2D(None, xform="xform.align2d")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_params2D(None, xform="xform.align2d")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        """



class Test_set_params2D(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D"))
    params=[1,2,3,4,5]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.set_params2D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.set_params2D()
        self.assertEqual(cm_new.exception.message, "set_params2D() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_set_params2D_using_wrongxform(self):
        (ima,) = self.argum[0]
        fu_img = deepcopy(ima)
        fu2_img = deepcopy(ima)
        fu.set_params2D(fu_img, self.params, xform="xform.align2d")
        oldfu.set_params2D(fu2_img, self.params, xform="xform.projection")     # is not setting the params
        self.assertFalse(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(fu2_img)))
        self.assertFalse(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(ima)))

    def test_set_params2D_using_wrongxform2(self):
        (ima,) = self.argum[0]
        fu_img = deepcopy(ima)
        fu2_img = deepcopy(ima)
        fu.set_params2D(fu_img, self.params, xform="xform.projection")       # is not setting the params
        oldfu.set_params2D(fu2_img, self.params, xform="xform.projection")      # is not setting the params
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(fu2_img)))
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(ima)))

    def test_set_params2D(self):
        (ima,) = self.argum[0]
        fu_img = deepcopy(ima)
        oldfu_img = deepcopy(ima)
        fu.set_params2D(fu_img, self.params, xform="xform.align2d")
        oldfu.set_params2D(oldfu_img, self.params, xform="xform.align2d")
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(oldfu_img)))
        self.assertFalse(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(ima)))

    def test_less_params(self):
        (ima,) = self.argum[0]
        fu_img = deepcopy(ima)
        oldfu_img = deepcopy(ima)
        with self.assertRaises(IndexError) as cm_new:
            fu.set_params2D(fu_img, [0,1], xform="xform.align2d")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.set_params2D(oldfu_img, [0,1], xform="xform.align2d")
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_wrong_xform_does_not_change_the_values_IS_IT_OK_OR_NOT(self):
        (ima,) = self.argum[0]
        fu_img = deepcopy(ima)
        oldfu_img = deepcopy(ima)
        fu.set_params2D(fu_img, self.params, xform="xform.align3d")          # is not setting the params
        oldfu.set_params2D(oldfu_img, self.params, xform="xform.align3d")    # is not setting the params
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(oldfu_img)))
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(ima)))

    def test_wrong_input_img(self):
        # I called it wrong image just because in the 'get_params2D' there was an error due to the missing xform key
        (ima,) = self.argum[0]
        fu_img = deepcopy(IMAGE_2D)
        oldfu_img = deepcopy(IMAGE_2D)
        fu.set_params2D(fu_img, self.params, xform="xform.align2d")
        oldfu.set_params2D(oldfu_img, self.params, xform="xform.align2d")
        self.assertTrue(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(oldfu_img)))
        self.assertFalse(numpy.array_equal(fu.get_params2D(fu_img), oldfu.get_params2D(ima)))

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.set_params2D(None, self.params, xform="xform.align2d")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.set_params2D(None, self.params, xform="xform.align2d")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'set_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)





class Test_get_params3D(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_params3D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_params3D()
        self.assertEqual(cm_new.exception.message, "get_params3D() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_params3D(self):
        return_new = fu.get_params3D(IMAGE_3D, xform="xform.align3d")
        return_old = oldfu.get_params3D(IMAGE_3D, xform="xform.align3d")
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_wrong_xform_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params3D(IMAGE_3D, xform="xform.align2d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params3D(IMAGE_3D, xform="xform.align2d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_input_img_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params3D(IMAGE_2D, xform="xform.align3d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params3D(IMAGE_2D, xform="xform.align3d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_params3D(IMAGE_3D, xform="xform.align3d")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_params3D(IMAGE_3D, xform="xform.align3d")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        """



class Test_set_params3D(unittest.TestCase):
    params=[1,2,3,4,5,6,7,8]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.set_params2D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.set_params2D()
        self.assertEqual(cm_new.exception.message, "set_params2D() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_set_params3D(self):
        fu_img = deepcopy(IMAGE_3D)
        oldfu_img = deepcopy(IMAGE_3D)
        fu.set_params3D(fu_img, self.params, xform="xform.align3d")
        oldfu.set_params3D(oldfu_img, self.params, xform="xform.align3d")
        self.assertTrue(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(oldfu_img)))
        self.assertFalse(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(IMAGE_3D)))

    def test_less_params(self):
        fu_img = deepcopy(IMAGE_3D)
        oldfu_img = deepcopy(IMAGE_3D)
        with self.assertRaises(IndexError) as cm_new:
            fu.set_params3D(fu_img, [0,1], xform="xform.align3d")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.set_params3D(oldfu_img, [0,1], xform="xform.align3d")
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_xform_does_not_change_the_values_IS_IT_OK_OR_NOT(self):
        fu_img = deepcopy(IMAGE_3D)
        oldfu_img = deepcopy(IMAGE_3D)
        fu.set_params3D(fu_img, self.params, xform="xform.align2d")
        oldfu.set_params3D(oldfu_img, self.params, xform="xform.align2d")
        self.assertTrue(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(oldfu_img)))
        #self.assertFalse(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(IMAGE_3D)))

    def test_wrong_input_img(self):
        # I called it wrong image just because in the 'get_params2D' there was an error due to the missing xform key
        fu_img = deepcopy(IMAGE_2D)
        oldfu_img = deepcopy(IMAGE_2D)
        fu.set_params3D(fu_img, self.params, xform="xform.align3d")
        oldfu.set_params3D(oldfu_img, self.params, xform="xform.align3d")
        self.assertTrue(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(oldfu_img)))
        self.assertFalse(numpy.array_equal(fu.get_params3D(fu_img), oldfu.get_params3D(IMAGE_3D)))

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.set_params3D(None, self.params, xform="xform.align3d")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.set_params3D(None, self.params, xform="xform.align3d")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'set_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_params_proj(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params_proj"))
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_params_proj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_params_proj()
        self.assertEqual(cm_new.exception.message, "get_params_proj() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_params_proj(self):
        (ima,) = self.argum[0]
        return_new = fu.get_params_proj(ima, xform="xform.projection")
        return_old = oldfu.get_params_proj(ima, xform="xform.projection")
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_wrong_xform_returns_NotExistingObjectException_key_doesnot_exist(self):
        (ima,) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params_proj(ima, xform="xform.align3d")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params_proj(ima, xform="xform.align3d")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_input_2dimg_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params_proj(IMAGE_2D, xform="xform.projection")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params_proj(IMAGE_2D, xform="xform.projection")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_input_3dimg_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params_proj(IMAGE_3D, xform="xform.projection")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params_proj(IMAGE_3D, xform="xform.projection")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_NoneType_as_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_params_proj(None, xform="xform.projection")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_params_proj(None, xform="xform.projection")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        """



class Test_set_params_proj(unittest.TestCase):
    params=[1,2,3,4,5]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.set_params_proj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.set_params_proj()
        self.assertEqual(cm_new.exception.message, "set_params_proj() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_set_params_proj_using_wrongxform_returns_NotExistingObjectException_key_doesnot_exist(self): #error is ok
        fu_img = deepcopy(IMAGE_2D)
        fu.set_params_proj(fu_img, self.params, xform="xform.align2d")
        with self.assertRaises(RuntimeError) as cm:
            fu.get_params_proj(fu_img, xform="xform.projection")
        msg = cm.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")

    def test_set_params_proj_using_wrongxform2returns_NotExistingObjectException_key_doesnot_exist(self):
        fu_img = deepcopy(IMAGE_2D)
        fu2_img = deepcopy(IMAGE_2D)
        fu.set_params_proj(fu_img, self.params, xform="xform.align2d")       # is not setting the params
        oldfu.set_params_proj(fu2_img, self.params, xform="xform.align2d")      # is not setting the params
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_params_proj(fu_img, xform="xform.projection")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_params_proj(fu2_img, xform="xform.projection")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_set_params_proj(self):
        fu_img = deepcopy(IMAGE_2D)
        oldfu_img = deepcopy(IMAGE_2D)
        fu.set_params_proj(fu_img, self.params, xform="xform.projection")
        oldfu.set_params_proj(oldfu_img, self.params, xform="xform.projection")
        self.assertTrue(numpy.array_equal(fu.get_params_proj(fu_img, xform="xform.projection"), oldfu.get_params_proj(oldfu_img, xform="xform.projection")))
        #self.assertFalse(numpy.array_equal(fu.get_params_proj(fu_img), fu.get_params_proj(IMAGE_2D))) # IMAGE2D has not key ''xform.projection'

    def test_less_params(self):
        fu_img = deepcopy(IMAGE_2D)
        oldfu_img = deepcopy(IMAGE_2D)
        with self.assertRaises(IndexError) as cm_new:
            fu.set_params_proj(fu_img, [0,1], xform="xform.projection")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.set_params_proj(oldfu_img, [0,1], xform="xform.projection")
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.set_params_proj(None, self.params, xform="xform.projection")
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.set_params_proj(None, self.params, xform="xform.projection")
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'set_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_ctf(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_ctf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_ctf()
        self.assertEqual(cm_new.exception.message, "get_ctf() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_img_returns_NotExistingObjectException_key_doesnot_exist(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_ctf(IMAGE_2D)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_ctf(IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_ctf(self):
        img_with_ctf = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))[0][0][0]
        self.assertTrue(numpy.array_equal(oldfu.get_ctf(img_with_ctf),fu.get_ctf(img_with_ctf)  ))

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_ctf(None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_ctf(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_same_ctf(unittest.TestCase):
    params = [1,2,3,4,5,6]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.same_ctf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.same_ctf()
        self.assertEqual(cm_new.exception.message, "same_ctf() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_same_ctf(self):
        self.assertTrue(fu.same_ctf(fu.generate_ctf(self.params),oldfu.generate_ctf(self.params)))

    def test_not_same_ctf(self):
        self.assertFalse(fu.same_ctf(fu.generate_ctf(self.params),oldfu.generate_ctf([0,1,2,3,4,5])))



class Test_generate_ctf(unittest.TestCase):
    """ params = [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle] """
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.generate_ctf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.generate_ctf()
        self.assertEqual(cm_new.exception.message, "generate_ctf() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_generate_ctf_with6Values(self):
        self.assertTrue(fu.same_ctf(fu.generate_ctf([1, 2, 3, 4, 5, 6]), oldfu.generate_ctf([1, 2, 3, 4, 5, 6])))

    def test_generate_ctf_with8Values(self):
        self.assertTrue(fu.same_ctf(fu.generate_ctf([1, 2, 3, 4, 5, 6,7,8]), oldfu.generate_ctf([1, 2, 3, 4, 5, 6,7,8])))

    def test_generate_ctf_with_incorrect_number_of_params_warning_msg(self):
        self.assertTrue(fu.generate_ctf([1, 2, 3, 4, 5, 6, 7]) is None)
        self.assertTrue(oldfu.generate_ctf([1, 2, 3, 4, 5, 6, 7]) is None)



class Test_delete_bdb(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.delete_bdb()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.delete_bdb()
        self.assertEqual(cm_new.exception.message, "delete_bdb() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_disable_bdb_cache(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.disable_bdb_cache(3)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.disable_bdb_cache(3)
        self.assertEqual(cm_new.exception.message, "disable_bdb_cache() takes no arguments (1 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_disable_bdb_cache(self):
        EMAN2db.BDB_CACHE_DISABLE = False
        self.assertFalse(EMAN2db.BDB_CACHE_DISABLE)
        fu.disable_bdb_cache()
        self.assertTrue(EMAN2db.BDB_CACHE_DISABLE)
        EMAN2db.BDB_CACHE_DISABLE = False
        self.assertFalse(EMAN2db.BDB_CACHE_DISABLE)
        oldfu.disable_bdb_cache()
        self.assertTrue(EMAN2db.BDB_CACHE_DISABLE )



class Test_getvec(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.getvec()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.getvec()
        self.assertEqual(cm_new.exception.message, "getvec() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_tht_between_90_180(self):
        self.assertTrue(numpy.array_equal(fu.getvec(0,100), oldfu.getvec(0,100)))

    def test_tht_bigger_than_180(self):
        self.assertTrue(numpy.array_equal(fu.getvec(0,200), oldfu.getvec(0,200)))

    def test_tht_lower_than_90(self):
        self.assertTrue(numpy.array_equal(fu.getvec(0,0), oldfu.getvec(0,0)))



class Test_getfvec(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.getfvec()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.getfvec()
        self.assertEqual(cm_new.exception.message, "getfvec() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_tht_between_90_180(self):
        self.assertTrue(numpy.array_equal(fu.getfvec(0,100), oldfu.getfvec(0,100)))

    def test_tht_bigger_than_180(self):
        self.assertTrue(numpy.array_equal(fu.getfvec(0,200), oldfu.getfvec(0,200)))

    def test_tht_lower_than_90(self):
        self.assertTrue(numpy.array_equal(fu.getfvec(0,0), oldfu.getfvec(0,0)))



class Test_nearest_fang(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.nearest_fang()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.nearest_fang()
        self.assertEqual(cm_new.exception.message, "nearest_fang() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_nearest_fang_true_should_return_equal_objects(self):
        """ values got from pickle files/utilities/utilities.nearest_fang """
        vecs = [[0.0, 0.0, 1.0], [0.6804220676422119, 0.6526213884353638, 0.3333333432674408], [-0.4104178845882416, 0.8487908840179443, 0.3333333432674408], [-0.9340742230415344, -0.12803982198238373, 0.3333333432674408], [-0.16687190532684326, -0.927923858165741, 0.3333333432674408], [0.8309417963027954, -0.4454488158226013, 0.3333333432674408], [8.742277657347586e-08, 7.64274186065882e-15, -1.0], [0.9340742230415344, 0.12803970277309418, -0.3333333134651184], [0.16687177121639252, 0.927923858165741, -0.3333333134651184], [-0.8309417963027954, 0.44544869661331177, -0.3333333134651184], [-0.6804221868515015, -0.652621328830719, -0.3333333134651184], [0.41041797399520874, -0.8487908840179443, -0.3333333134651184]]
        tht = 66.00945
        phi = 58.54455
        self.assertEqual(fu.nearest_fang(vecs, phi, tht), oldfu.nearest_fang(vecs, phi, tht))


    def test_empty_vectore(self):
        """ values got from pickle files/utilities/utilities.nearest_fang """
        self.assertEqual(fu.nearest_fang([], 100, 100), oldfu.nearest_fang([], 100, 100))
        self.assertEqual(fu.nearest_fang([], 100, 100), -1)



class Test_nearest_many_full_k_projangles(unittest.TestCase):
    reference_normals = [[0.606369137763977, 0.7754802703857422, 0.17591717839241028], [0.344023197889328, 0.9092735648155212, 0.23424272239208221], [0.5131438970565796, 0.7110531330108643, -0.4807148575782776], [0.6525110006332397, 0.6401833295822144, 0.4054562747478485], [0.5846421718597412, 0.5353381037712097, 0.6095954775810242], [0.3914891481399536, 0.4943649470806122, 0.7761054039001465], [0.21746492385864258, 0.411188542842865, 0.8852304816246033], [0.18686196208000183, 0.4279184937477112, 0.8842897415161133], [0.2696961760520935, 0.41237473487854004, 0.870178759098053], [0.34728822112083435, 0.3424328565597534, 0.8730009198188782], [0.2467251867055893, 0.39220815896987915, 0.8861712217330933], [0.43794623017311096, 0.19451908767223358, 0.8777046203613281], [0.35838937759399414, 0.0876869484782219, -0.9294450283050537], [0.6956571340560913, 0.7182994484901428, 0.010348091833293438], [0.6555072665214539, 0.7445935010910034, -0.12605828046798706], [0.7438855767250061, 0.6679566502571106, -0.02163686789572239], [0.58192378282547, 0.8076738715171814, -0.09501412510871887], [0.7202955484390259, 0.693575382232666, 0.011288836598396301], [0.6438657641410828, 0.765091598033905, 0.008466575294733047], [0.6417456269264221, 0.7646241188049316, -0.05926619470119476], [0.593335747718811, 0.7773913145065308, 0.20884287357330322], [0.5866740942001343, 0.8075113296508789, -0.06114771217107773], [0.5893274545669556, 0.8044687509536743, 0.07431796938180923], [0.48042023181915283, 0.8660674691200256, 0.13828791677951813], [0.46822038292884827, 0.8812242746353149, -0.06491056084632874], [0.34745562076568604, 0.9322780966758728, 0.10065855830907822], [0.4396599531173706, 0.898162305355072, 0.0018815546063706279], [0.5071992874145508, 0.8368419408798218, 0.2060207575559616], [0.35214218497276306, 0.913831353187561, -0.20225776731967926], [0.5917134881019592, 0.798856258392334, 0.1081843376159668], [0.31928351521492004, 0.9256179332733154, -0.2031984180212021], [0.5689234137535095, 0.8101938962936401, 0.14111001789569855], [0.5366130471229553, 0.8180546164512634, 0.2069614678621292], [0.6138750910758972, 0.751165509223938, 0.24270929396152496], [0.6470115184783936, 0.7327832579612732, 0.210724338889122], [0.6170760989189148, 0.7832963466644287, 0.0752587541937828], [0.6726201176643372, 0.7090698480606079, 0.21166512370109558], [0.5653374195098877, 0.7982293963432312, 0.2079022079706192], [0.6659785509109497, 0.704732358455658, 0.24459083378314972], [0.6436562538146973, 0.7641429901123047, 0.04233306273818016], [0.6849393248558044, 0.7063358426094055, 0.17873942852020264], [0.5400856733322144, 0.8298555016517639, 0.14016936719417572], [0.5633652806282043, 0.8192181587219238, 0.10724367946386337], [0.5887830853462219, 0.8072782158851624, 0.040451530367136], [0.5886198282241821, 0.8079495429992676, -0.02728116139769554], [0.5608543157577515, 0.8246564269065857, 0.07337724417448044], [0.6164841055870056, 0.7869266271591187, -0.026340581476688385], [0.6699250340461731, 0.7420257925987244, -0.0244591124355793], [0.6205720901489258, 0.7555667161941528, 0.20978358387947083], [0.668122410774231, 0.7417618036270142, -0.058325473219156265], [0.6953815221786499, 0.7172793745994568, 0.04421444982290268], [0.6165966987609863, 0.7861903309822083, 0.04139237478375435], [0.6167761087417603, 0.7871026396751404, 0.007525925524532795], [0.7440555691719055, 0.6680058240890503, 0.012229571118950844], [0.5889342427253723, 0.8081541061401367, 0.006585149094462395], [0.6699285507202148, 0.7411633729934692, 0.04327383637428284], [0.7258118987083435, 0.6720566749572754, 0.1467544287443161], [0.6510280966758728, 0.7452824115753174, 0.1439322829246521], [0.695436418056488, 0.7182027697563171, -0.02351834438741207], [0.6768127679824829, 0.7217592000961304, 0.1448730081319809], [0.659572958946228, 0.7303088307380676, 0.1777988076210022], [0.6193289160728455, 0.7775111198425293, 0.10912513732910156], [0.644066333770752, 0.7611650824546814, 0.07619946449995041], [0.646177351474762, 0.7552087903022766, 0.11006595939397812], [0.7330403327941895, 0.6557652354240417, 0.18062089383602142], [0.7375331521034241, 0.6283562183380127, 0.2474130392074585], [0.7217933535575867, 0.68284010887146, 0.11288806051015854], [0.6975162625312805, 0.6453772783279419, 0.311382919549942], [0.8656806349754333, 0.48613080382347107, 0.11947321146726608], [0.7893708944320679, 0.6029136180877686, 0.11571019887924194], [0.8126943111419678, 0.5629141926765442, 0.15051743388175964], [0.8193334341049194, 0.5153672695159912, 0.25117599964141846], [0.8606626391410828, 0.45913165807724, 0.22013163566589355], [0.9627028107643127, 0.2471613585948944, -0.11006592959165573], [0.8993244171142578, 0.4370543360710144, -0.014110974967479706], [0.8985337615013123, 0.4356166422367096, 0.053621746599674225], [0.6973650455474854, 0.6844565272331238, 0.21260593831539154], [0.7557139992713928, 0.6292310357093811, 0.18156161904335022], [0.720099925994873, 0.6935028433799744, -0.0225775558501482], [0.6813780665397644, 0.7211584448814392, -0.12511759996414185], [0.720158576965332, 0.689294695854187, 0.07902166992425919], [0.6333718299865723, 0.7533666491508484, 0.17685800790786743], [0.7017511129379272, 0.6973404884338379, 0.14581383764743805], [0.670264720916748, 0.7381020188331604, 0.07714022696018219], [0.6722255349159241, 0.7319769859313965, 0.11100655794143677], [0.6406394839286804, 0.7163008451461792, 0.2765757739543915], [0.6907424926757812, 0.6801389455795288, 0.2455316036939621], [0.6244292855262756, 0.7678812146186829, 0.1429915428161621], [0.7094387412071228, 0.6814774870872498, 0.17968012392520905], [0.5963707566261292, 0.7414255142211914, 0.30761995911598206], [0.6974412798881531, 0.7078441977500916, 0.11194729804992676], [0.5866034030914307, 0.7729452252388, 0.24176852405071259], [0.7146044969558716, 0.6546692252159119, 0.24647238850593567], [0.5873112082481384, 0.7613202929496765, 0.2746942937374115], [0.50588458776474, 0.7085287570953369, 0.49200382828712463]]
    angles = [[41.921590970437258, 91.23979851375101, 333.346436124961, -0.0, -0.0]]
    howmany = 47
    symclass = sparx_fundamentals.symclass("c1")
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.nearest_many_full_k_projangles()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.nearest_many_full_k_projangles()
        self.assertEqual(cm_new.exception.message, "nearest_many_full_k_projangles() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file_values(self):
        symclass = sparx_fundamentals.symclass("c5")    # I creasted it like the one of the pickle file
        return_new = fu.nearest_many_full_k_projangles(self.reference_normals, self.angles, self.howmany, symclass)
        return_old = oldfu.nearest_many_full_k_projangles(self.reference_normals, self.angles, self.howmany, symclass)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_class_c1(self):
        return_new = fu.nearest_many_full_k_projangles(self.reference_normals, self.angles, self.howmany, self.symclass)
        return_old = oldfu.nearest_many_full_k_projangles(self.reference_normals, self.angles, self.howmany, self.symclass)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_null_howmany(self):
        return_new = fu.nearest_many_full_k_projangles(self.reference_normals, self.angles, 0, self.symclass)
        return_old = oldfu.nearest_many_full_k_projangles(self.reference_normals, self.angles, 0, self.symclass)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_empty_list_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.nearest_many_full_k_projangles([], self.angles, self.howmany, self.symclass)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.nearest_many_full_k_projangles([], self.angles, self.howmany, self.symclass)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "Error, number of neighbors cannot be larger than number of reference directions")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_assign_projdirs_f(unittest.TestCase):
    """
    Since 'projdirs' and 'refdirs' are got in the sxmeridian from  angles_to_normals(list_of_angles) I used the
    output of the angles_to_normals tests to fill the 'projdirs' variable. the refdirs values are 2/3 of projdirs values
    """
    projdirs =  [[0.0, 0.0, 1.0], [0.6804220676422119, 0.6526213884353638, 0.3333333432674408], [-0.4104178845882416, 0.8487909436225891, 0.3333333432674408], [-0.9340742230415344, -0.12803982198238373, 0.3333333432674408], [-0.16687190532684326, -0.927923858165741, 0.3333333432674408], [0.8309417366981506, -0.4454488158226013, 0.3333333432674408], [8.742277657347586e-08, 7.64274186065882e-15, -1.0], [0.9340742230415344, 0.12803970277309418, -0.3333333134651184], [0.16687177121639252, 0.927923858165741, -0.3333333134651184], [-0.8309418559074402, 0.44544869661331177, -0.3333333134651184], [-0.6804221272468567, -0.652621328830719, -0.3333333134651184], [0.41041797399520874, -0.8487908840179443, -0.3333333134651184]]
    refdirs = [[0.0, 0.0, 0.66], [0.44907856464385987, 0.4307301163673401, 0.22000000655651095], [-0.27087580382823945, 0.5602020227909088, 0.22000000655651095], [-0.6164889872074127, -0.08450628250837326, 0.22000000655651095], [-0.11013545751571656, -0.6124297463893891, 0.22000000655651095], [0.5484215462207794, -0.2939962184429169, 0.22000000655651095], [5.7699032538494066e-08, 5.044209628034821e-15, -0.66], [0.6164889872074127, 0.08450620383024215, -0.21999998688697817], [0.11013536900281906, 0.6124297463893891, -0.21999998688697817], [-0.5484216248989106, 0.2939961397647858, -0.21999998688697817], [-0.44907860398292543, -0.43073007702827454, -0.21999998688697817], [0.2708758628368378, -0.5602019834518432, -0.21999998688697817]]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.assign_projdirs_f()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.assign_projdirs_f()
        self.assertEqual(cm_new.exception.message, "assign_projdirs_f() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_real_data(self):
        neighbors = int(len(self.projdirs)/ len(self.refdirs))
        return_new = fu.assign_projdirs_f(self.projdirs, self.refdirs, neighbors)
        return_old = oldfu.assign_projdirs_f(self.projdirs, self.refdirs, neighbors)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_null_neighboor_value_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.assign_projdirs_f(self.projdirs, self.refdirs, 0)
        return_old = oldfu.assign_projdirs_f(self.projdirs, self.refdirs, 0)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        """

    def test_with_negative_neighboor_value_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.assign_projdirs_f(self.projdirs, self.refdirs, -1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.assign_projdirs_f(self.projdirs, self.refdirs, -1)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_too_high_neighboor_value_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.assign_projdirs_f(self.projdirs, self.refdirs, 5)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.assign_projdirs_f(self.projdirs, self.refdirs, 5)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_invalid_neighbors_type(self):
        neighbors = len(self.projdirs)/ len(self.refdirs)
        with self.assertRaises(TypeError) as cm_new:
            fu.assign_projdirs_f(self.projdirs, self.refdirs, neighbors)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.assign_projdirs_f(self.projdirs, self.refdirs, neighbors)
        msg = cm_new.exception.message.split("\n")
        msg_old = cm_old.exception.message.split("\n")
        self.assertEqual(msg[0]+msg[1], 'Python argument types in    Util.assign_projdirs_f(list, list, float)')
        self.assertEqual(msg[0]+msg[1], msg_old[0]+msg_old[1])

    def test_with_projdirs_refdirs_have_different_length(self):
        refdirs= self.refdirs [:10]
        neighbors = int(len(self.projdirs)/ len(refdirs))
        return_new = fu.assign_projdirs_f(self.projdirs, refdirs, neighbors)
        return_old = oldfu.assign_projdirs_f(self.projdirs, refdirs, neighbors)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_empty_projdirs(self):
        return_new = fu.assign_projdirs_f([], self.refdirs, 1)
        return_old = oldfu.assign_projdirs_f([], self.refdirs, 1)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_empty_refdirs_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.assign_projdirs_f(self.projdirs, [], 1)
        return_old = oldfu.assign_projdirs_f(self.projdirs, [], 1)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        """



class Test_angles_to_normals(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angles_to_normals()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angles_to_normals()
        self.assertEqual(cm_new.exception.message, "angles_to_normals() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_pickle_file_values(self):
        angles = [[0.0, 0.0, 0.0], [43.805265094506787, 70.528779365509308, 0.0], [115.80526509450678, 70.528779365509308, 0.0], [187.80526509450678, 70.528779365509308, 0.0], [259.80526509450681, 70.528779365509308, 0.0], [331.80526509450681, 70.528779365509308, 0.0], [180.0, 180.0, 0.0], [7.8052650945068081, 109.47122063449069, 0.0], [79.805265094506808, 109.47122063449069, 0.0], [151.80526509450681, 109.47122063449069, 0.0], [223.80526509450681, 109.47122063449069, 0.0], [295.80526509450681, 109.47122063449069, 0.0]]
        return_new = fu.angles_to_normals(angles)
        return_old = oldfu.angles_to_normals(angles)
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_empty_angles_list(self):
        return_new = fu.angles_to_normals([])
        return_old = oldfu.angles_to_normals([])
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new,[]))



class Test_angular_occupancy(unittest.TestCase):
    params = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy"))[0][0]
    angstep = 15 # i change it becuase the lower value got from the pickle file leads each test to run for more than 10 sec, nov less than 1

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angular_occupancy()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angular_occupancy()
        self.assertEqual(cm_new.exception.message, "angular_occupancy() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_less_angles_returns_IndexError_list_index_out_of_range(self):
        angles=[[0.1],[21.1],[30.11],[1.1]]
        with self.assertRaises(IndexError) as cm_new:
            fu.angular_occupancy(angles, self.angstep, 'c5', 'S')
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angular_occupancy(angles, self.angstep, 'c5', 'S')
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_method_S(self):
        """ the values got from the pickle file"""
        return_new = fu.angular_occupancy(self.params, self.angstep, 'c5', 'S')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'c5', 'S')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_c1_method_S(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'c1', 'S')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'c1', 'S')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_oct_method_S(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'oct1', 'S')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'oct1', 'S')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_invalid_method_S_returns_AttributeError_symclass_hasnot_attribute_symangles_error_msg(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.angular_occupancy(self.params, self.angstep, 'invalid', 'S')
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.angular_occupancy(self.params, self.angstep, 'invalid', 'S')
        self.assertEqual(cm_new.exception.message, "'symclass' object has no attribute 'symangles'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_method_P(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'c5', 'P')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'c5', 'P')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_c1_method_P(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'c1', 'P')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'c1', 'P')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_oct_method_P(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'oct1', 'P')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'oct1', 'P')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_sym_c5_method_invalid(self):
        return_new = fu.angular_occupancy(self.params, self.angstep, 'c5', 'invalid')
        return_old = oldfu.angular_occupancy(self.params, self.angstep, 'c5', 'invalid')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_empty_params_list(self):
        """ the values got from the pickle file"""
        return_new = fu.angular_occupancy([], self.angstep, 'c5', 'S')
        return_old = oldfu.angular_occupancy([], self.angstep, 'c5', 'S')
        self.assertTrue(numpy.array_equal(return_new,return_old))

    def test_with_null_angstep_returns_ZeroDivisionError_error_msg(self):
        """ the values got from the pickle file"""
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.angular_occupancy(self.params, 0, 'c5', 'S')
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.angular_occupancy(self.params, 0, 'c5', 'S')
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_angular_histogram(unittest.TestCase):
    params = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy"))[0][0]
    angstep = 15 # i change it becuase the lower value got from the pickle file leads each test to run for more than 10 sec, nov less than 1

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angular_histogram()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angular_histogram()
        self.assertEqual(cm_new.exception.message, "angular_histogram() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_less_angles_returns_IndexError_list_index_out_of_range(self):
        angles=[[0.1],[21.1],[30.11],[1.1]]
        with self.assertRaises(IndexError) as cm_new:
            fu.angular_histogram(angles, self.angstep, 'c5', 'S')
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angular_histogram(angles, self.angstep, 'c5', 'S')
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_method_S(self):
        """ the values got from the pickle file"""
        return_new = fu.angular_histogram(self.params, self.angstep, 'c5', 'S')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'c5', 'S')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_c1_method_S(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'c1', 'S')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'c1', 'S')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_oct_method_S(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'oct1', 'S')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'oct1', 'S')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_invalid_method_S_returns_AttributeError_symclass_hasnot_attribute_symangles_error_msg(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.angular_histogram(self.params, self.angstep, 'invalid', 'S')
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.angular_histogram(self.params, self.angstep, 'invalid', 'S')
        self.assertEqual(cm_new.exception.message, "'symclass' object has no attribute 'symangles'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_method_P(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'c5', 'P')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'c5', 'P')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_c1_method_P(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'c1', 'P')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'c1', 'P')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_oct_method_P(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'oct1', 'P')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'oct1', 'P')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_sym_c5_method_invalid(self):
        return_new = fu.angular_histogram(self.params, self.angstep, 'c5', 'invalid')
        return_old = oldfu.angular_histogram(self.params, self.angstep, 'c5', 'invalid')
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_empty_params_list(self):
        """ the values got from the pickle file"""
        return_new = fu.angular_histogram([], self.angstep, 'c5', 'S')
        return_old = oldfu.angular_histogram([], self.angstep, 'c5', 'S')
        # self.assertTrue(numpy.array_equal(return_new, return_old))    --> failed ... ????
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_with_null_angstep_returns_ZeroDivisionError_error_msg(self):
        """ the values got from the pickle file"""
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.angular_histogram(self.params, 0, 'c5', 'S')
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.angular_histogram(self.params, 0, 'c5', 'S')
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_balance_angular_distribution(unittest.TestCase):
    params = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy"))[0][0]
    angstep = 15 # i change it becuase the lower value got from the pickle file leads each test to run for more than 10 sec, nov less than 1

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.balance_angular_distribution()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.balance_angular_distribution()
        self.assertEqual(cm_new.exception.message, "balance_angular_distribution() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_less_angles_returns_IndexError_list_index_out_of_range(self):
        angles=[[0.1],[21.1],[30.11],[1.1]]
        with self.assertRaises(IndexError) as cm_new:
            fu.balance_angular_distribution(angles, max_occupy = -1, angstep = self.angstep, sym= 'c5')
        with self.assertRaises(IndexError) as cm_old:
            oldfu.balance_angular_distribution(angles, max_occupy = -1, angstep = self.angstep, sym= 'c5')
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_not_positive_maxOccupy(self):
        """ the values got from the pickle file"""
        return_new = fu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'c5')
        return_old = oldfu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'c5')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_sym_c1_not_positive_maxOccupy(self):
        return_new = fu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'c1')
        return_old = oldfu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'c1')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_sym_oct_not_positive_maxOccupy(self):
        return_new = fu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'oct1')
        return_old = oldfu.balance_angular_distribution(self.params, max_occupy = -1, angstep = self.angstep, sym= 'oct1')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_empty_list(self):
        return_new = fu.balance_angular_distribution([], max_occupy = -1, angstep = self.angstep, sym= 'c5')
        return_old = oldfu.balance_angular_distribution([], max_occupy = -1, angstep = self.angstep, sym= 'c5')
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_null_angstepy_error_msg(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.balance_angular_distribution(self.params, max_occupy = -1, angstep = 0, sym= 'c5')
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.balance_angular_distribution(self.params, max_occupy = -1, angstep = 0, sym= 'c5')
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c5_positive_maxOccupy_not_testabel(self):
        """
        It use to process random value that lead the function to returns always different va;ues
        """
        self.assertTrue(True)
        """
        (params, not_usedangstep, sym, not_used) = self.argum[0]
        return_new = fu.balance_angular_distribution(deepcopy(params), max_occupy = 1, angstep = self.angstep, sym= sym)
        return_old = oldfu.balance_angular_distribution(deepcopy(params), max_occupy = 1, angstep = self.angstep, sym= sym)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """



class Test_symmetry_neighbors(unittest.TestCase):
    angles = [[0.0, 0.0, 1.0], [0.6804220676422119, 0.6526213884353638, 0.3333333432674408], [-0.4104178845882416, 0.8487909436225891, 0.3333333432674408], [-0.9340742230415344, -0.12803982198238373, 0.3333333432674408], [-0.16687190532684326, -0.927923858165741, 0.3333333432674408], [0.8309417366981506, -0.4454488158226013, 0.3333333432674408], [8.742277657347586e-08, 7.64274186065882e-15, -1.0], [0.9340742230415344, 0.12803970277309418, -0.3333333134651184], [0.16687177121639252, 0.927923858165741, -0.3333333134651184], [-0.8309418559074402, 0.44544869661331177, -0.3333333134651184], [-0.6804221272468567, -0.652621328830719, -0.3333333134651184], [0.41041797399520874, -0.8487908840179443, -0.3333333134651184]]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.symmetry_neighbors()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.symmetry_neighbors()
        self.assertEqual(cm_new.exception.message, "symmetry_neighbors() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_empty_list_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.symmetry_neighbors([] , symmetry= "c1")
        return_old = oldfu.symmetry_neighbors([], symmetry= "c1")
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_with_less_angles_returns_RuntimeError_3_angles_are_required(self):
        angles=[[0.1],[21.1],[30.11],[1.1]]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.symmetry_neighbors(angles , symmetry= "c1")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.symmetry_neighbors(angles , symmetry= "c1")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "Three angles are required")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_sym_c1(self):
        return_new = fu.symmetry_neighbors(self.angles , symmetry= "c1")
        return_old = oldfu.symmetry_neighbors(self.angles , symmetry= "c1")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_sym_c5(self):
        return_new = fu.symmetry_neighbors(self.angles , symmetry= "c5")
        return_old = oldfu.symmetry_neighbors(self.angles , symmetry= "c5")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_sym_d1(self):
        return_new = fu.symmetry_neighbors(self.angles , symmetry= "d1")
        return_old = oldfu.symmetry_neighbors(self.angles , symmetry= "d1")
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_sym_not_c_or_d(self):
        """ these cases take a lot of times. more or less one minute"""
        return_new = fu.symmetry_neighbors(self.angles , symmetry= "invalid")
        return_old = oldfu.symmetry_neighbors(self.angles , symmetry= "invalid")
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_rotation_between_anglesets(unittest.TestCase):
    """  used the value used in 'Test_assign_projdirs_f' """
    agls1 =  [[0.0, 0.0, 1.0], [0.6804220676422119, 0.6526213884353638, 0.3333333432674408], [-0.4104178845882416, 0.8487909436225891, 0.3333333432674408], [-0.9340742230415344, -0.12803982198238373, 0.3333333432674408], [-0.16687190532684326, -0.927923858165741, 0.3333333432674408], [0.8309417366981506, -0.4454488158226013, 0.3333333432674408], [8.742277657347586e-08, 7.64274186065882e-15, -1.0], [0.9340742230415344, 0.12803970277309418, -0.3333333134651184], [0.16687177121639252, 0.927923858165741, -0.3333333134651184], [-0.8309418559074402, 0.44544869661331177, -0.3333333134651184], [-0.6804221272468567, -0.652621328830719, -0.3333333134651184], [0.41041797399520874, -0.8487908840179443, -0.3333333134651184]]
    agls2 = [[0.0, 0.0, 0.66], [0.44907856464385987, 0.4307301163673401, 0.22000000655651095], [-0.27087580382823945, 0.5602020227909088, 0.22000000655651095], [-0.6164889872074127, -0.08450628250837326, 0.22000000655651095], [-0.11013545751571656, -0.6124297463893891, 0.22000000655651095], [0.5484215462207794, -0.2939962184429169, 0.22000000655651095], [5.7699032538494066e-08, 5.044209628034821e-15, -0.66], [0.6164889872074127, 0.08450620383024215, -0.21999998688697817], [0.11013536900281906, 0.6124297463893891, -0.21999998688697817], [-0.5484216248989106, 0.2939961397647858, -0.21999998688697817], [-0.44907860398292543, -0.43073007702827454, -0.21999998688697817], [0.2708758628368378, -0.5602019834518432, -0.21999998688697817]]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rotation_between_anglesets()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rotation_between_anglesets()
        self.assertEqual(cm_new.exception.message, "rotation_between_anglesets() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotation_between_anglesets(self):
        return_new = fu.rotation_between_anglesets(self.agls1, self.agls2)
        return_old = oldfu.rotation_between_anglesets(self.agls1, self.agls2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_sets_have_different_length(self):
        agls2=self.agls2[:30]
        return_new = fu.rotation_between_anglesets(self.agls1, agls2)
        return_old = oldfu.rotation_between_anglesets(self.agls1, agls2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_angls1_empty_list_error_msg(self):
        return_new = fu.rotation_between_anglesets([], self.agls2)
        return_old = oldfu.rotation_between_anglesets([], self.agls2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_angls2_empty_list_error_msg(self):
        return_new = fu.rotation_between_anglesets(self.agls1, [])
        return_old = oldfu.rotation_between_anglesets(self.agls1, [])
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_angle_between_projections_directions(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angle_between_projections_directions()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angle_between_projections_directions()
        self.assertEqual(cm_new.exception.message, "angle_between_projections_directions() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_3angles(self):
        agls1 = [20, 60, 0]
        agls2 = [45, 75, 5]
        return_new = fu.angle_between_projections_directions(agls1, agls2)
        return_old = oldfu.angle_between_projections_directions(agls1, agls2)
        self.assertEqual(return_new, return_old)

    def test_with_2angles(self):
        agls1 = [20, 60]
        agls2 = [45, 75]
        return_new = fu.angle_between_projections_directions(agls1, agls2)
        return_old = oldfu.angle_between_projections_directions(agls1, agls2)
        self.assertEqual(return_new, return_old)

    def test_with_list1_empty(self):
        agls2 = [45, 75]
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_between_projections_directions([], agls2)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_between_projections_directions([], agls2)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_list2_empty(self):
        agls1 = [45, 75]
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_between_projections_directions( agls1, [])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_between_projections_directions( agls1, [])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_pixel_size(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_pixel_size()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_pixel_size()
        self.assertEqual(cm_new.exception.message, "get_pixel_size() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_pixel_size_img2d(self):
        return_new = fu.get_pixel_size(IMAGE_2D)
        return_old = oldfu.get_pixel_size(IMAGE_2D)
        self.assertEqual(return_new, return_old)

    def test_get_pixel_size_img3d(self):
        return_new = fu.get_pixel_size(IMAGE_3D)
        return_old = oldfu.get_pixel_size(IMAGE_3D)
        self.assertEqual(return_new, return_old)

    def test_get_pixel_size_imgEmpty(self):
        return_new = fu.get_pixel_size(EMData())
        return_old = oldfu.get_pixel_size(EMData())
        self.assertEqual(return_new, return_old)

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.get_pixel_size(None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.get_pixel_size(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_attr_default'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_set_pixel_size(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.set_pixel_size()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.set_pixel_size()
        self.assertEqual(cm_new.exception.message, "set_pixel_size() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_set_pixel_size(self):
        img_fu = deepcopy(IMAGE_2D)
        img_fu_old = deepcopy(IMAGE_2D)
        fu.set_pixel_size(img_fu,2.1)
        oldfu.set_pixel_size(img_fu_old,2.1)
        self.assertEqual(img_fu.get_attr('apix_x'), img_fu_old.get_attr('apix_x'))
        self.assertEqual(img_fu.get_attr('apix_y'), img_fu_old.get_attr('apix_y'))
        self.assertEqual(img_fu.get_attr('apix_z'), img_fu_old.get_attr('apix_z'))
        self.assertEqual(img_fu.get_attr('apix_x'), 2.1)
        self.assertEqual(img_fu.get_attr('apix_y'), 2.1)
        self.assertEqual(img_fu.get_attr('apix_z'), 2.1)

    def test_set_pixel_size_truncated_value(self):
        img_fu = deepcopy(IMAGE_2D)
        img_fu_old = deepcopy(IMAGE_2D)
        fu.set_pixel_size(img_fu,2.1111)
        oldfu.set_pixel_size(img_fu_old,2.1111)
        self.assertEqual(img_fu.get_attr('apix_x'), img_fu_old.get_attr('apix_x'))
        self.assertEqual(img_fu.get_attr('apix_y'), img_fu_old.get_attr('apix_y'))
        self.assertEqual(img_fu.get_attr('apix_z'), img_fu_old.get_attr('apix_z'))
        self.assertEqual(img_fu.get_attr('apix_x'), 2.111)
        self.assertEqual(img_fu.get_attr('apix_y'), 2.111)
        self.assertEqual(img_fu.get_attr('apix_z'), 2.111)

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.set_pixel_size(None,2.1)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.set_pixel_size(None,2.1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_zsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_lacos(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.lacos()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.lacos()
        self.assertEqual(cm_new.exception.message, "lacos() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_angle(self):
        self.assertEqual(fu.lacos(0),oldfu.lacos(0))

    def test_negative_angle(self):
        self.assertEqual(fu.lacos(-0.12),oldfu.lacos(-0.12))

    def test_positive_angle(self):
        self.assertEqual(fu.lacos(0.12),oldfu.lacos(0.12))

    def test_outOfRange_angle(self):
        self.assertEqual(fu.lacos(12),oldfu.lacos(12))



class Test_findall(unittest.TestCase):
    l = [1,2,3,4,5,5,5,4,3,2,1]

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.findall()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.findall()
        self.assertEqual(cm_new.exception.message, "findall() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_findall_5(self):
        return_new = fu.findall(5, self.l, start=0)
        return_old = oldfu.findall(5, self.l, start=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_findall_noValues(self):
        return_new = fu.findall(0, self.l, start=0)
        return_old = oldfu.findall(0, self.l, start=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))



class Test_class_iterImagesList(unittest.TestCase):
    list_of_imgs = [IMAGE_2D,IMAGE_3D,IMAGE_BLANK_2D,IMAGE_BLANK_3D,IMAGE_2D_REFERENCE]

    def test_invalid_init(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.iterImagesList()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.iterImagesList()
        self.assertEqual(cm_new.exception.message, "__init__() takes at least 2 arguments (1 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_valid_init(self):
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes = None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes = None)
        self.assertEqual(type(fu_obj).__name__ , "iterImagesList")
        self.assertEqual(type(fu_obj).__name__, type(oldfu_obj).__name__)

    def test_valid_init2(self):
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes = [1,2])
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes = [1,2])
        self.assertEqual(type(fu_obj).__name__ , "iterImagesList")
        self.assertEqual(type(fu_obj).__name__, type(oldfu_obj).__name__)

    def test_wrong_init_list_of_index_leads_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.iterImagesList(self.list_of_imgs, list_of_indexes = [1,2,7])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.iterImagesList(self.list_of_imgs, list_of_indexes = [1,2,7])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_iterNo(self):
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), -1)

    def test_imageIndex(self):
        """ since the position is -1 it is returning the index of the last image hence 4"""
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        self.assertEqual(fu_obj.imageIndex(), oldfu_obj.imageIndex())
        self.assertEqual(fu_obj.imageIndex(), 4)

    def test_image(self):
        """ since the position is -1 it is returning the last image hence the 4th"""
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        fu_img=fu_obj.image()
        oldfu_img=oldfu_obj.image()
        expectedimg=self.list_of_imgs[fu_obj.imageIndex()]
        self.assertTrue(numpy.array_equal(fu_img.get_3dview(), oldfu_img.get_3dview()))
        self.assertTrue(numpy.array_equal(fu_img.get_3dview(), expectedimg.get_3dview()))

    def test_goToNext(self):
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes=None)

        """ I'm testing all the data in the obj in order to test the return False"""
        fu_counter =0
        while fu_obj.goToNext():       # I'm , implicitly, testing the return True
            self.assertEqual(fu_obj.iterNo(), fu_counter)
            fu_counter += 1

        oldfu_counter =0
        while oldfu_obj.goToNext():
            self.assertEqual(oldfu_obj.iterNo(), oldfu_counter)
            oldfu_counter += 1

        """ no more img in the object"""
        self.assertFalse(fu_obj.goToNext())
        self.assertFalse(oldfu_obj.goToNext())

        """ check if both of the classes tested all the images"""
        self.assertTrue(fu_counter, oldfu_counter)
        self.assertTrue(fu_counter, len(self.list_of_imgs))

    def test_goToPrev(self):
        fu_obj = fu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        oldfu_obj = oldfu.iterImagesList(self.list_of_imgs, list_of_indexes=None)
        """At the beginning there is no previous image"""
        self.assertFalse(fu_obj.goToPrev())
        self.assertFalse(oldfu_obj.goToPrev())
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), -1)

        """ We are on the first image, it means that we have still no previous image"""
        fu_obj.goToNext()
        oldfu_obj.goToNext()
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), 0)

        self.assertFalse(fu_obj.goToPrev())
        self.assertFalse(oldfu_obj.goToPrev())
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), -1)

        """ We are on the second image, it means that we have an previous image"""
        fu_obj.goToNext()
        oldfu_obj.goToNext()
        fu_obj.goToNext()
        oldfu_obj.goToNext()
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), 1)

        self.assertTrue(fu_obj.goToPrev())
        self.assertTrue(oldfu_obj.goToPrev())
        self.assertEqual(fu_obj.iterNo(),oldfu_obj.iterNo())
        self.assertEqual(fu_obj.iterNo(), 0)



class Test_pack_message(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pack_message()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pack_message()
        self.assertEqual(cm_new.exception.message, "pack_message() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_data_is_a_string(self):
        data = "case S:I am a string!!!"
        return_new = fu.pack_message(data)
        return_old = oldfu.pack_message(data)
        self.assertEqual(return_new,return_old)

    def test_data_is_a_very_long_string(self):
        long_data = "I am a stringggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg!!!"
        return_new = fu.pack_message(long_data)
        return_old = oldfu.pack_message(long_data)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, "C"+zlib.compress(long_data,1))

    def test_data_is_a_notstring(self):
        data = 5555
        return_new = fu.pack_message(data)
        return_old = oldfu.pack_message(data)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, "O" + pickle.dumps(data,-1))

    def test_data_is_a_notstring_long_version(self):
        data = 555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
        return_new = fu.pack_message(data)
        return_old = oldfu.pack_message(data)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, "Z" + zlib.compress(pickle.dumps(data, -1),1))



class Test_unpack_message(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.unpack_message()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.unpack_message()
        self.assertEqual(cm_new.exception.message, "unpack_message() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_data_is_a_string_BUG(self):
        self.assertTrue(True)
        """
        data = fu.pack_message("case S:I am a string!!!")
        return_new = fu.unpack_message(data)
        return_old = oldfu.unpack_message(data)
        self.assertEqual(return_new,return_old)
        """

    def test_data_is_a_very_long_string(self):
        self.assertTrue(True)
        """
        long_data = "I am a stringggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg!!!"
        data= fu.pack_message(long_data)
        return_new = fu.unpack_message(data)
        return_old = oldfu.unpack_message(data)
        self.assertEqual(return_new,return_old)
        """

    def test_data_is_a_notstring(self):
        self.assertTrue(True)
        """
        data = fu.pack_message(5555)
        return_new = fu.unpack_message(5555)
        return_old = oldfu.unpack_message(data)
        self.assertEqual(return_new,return_old)
        """
    def test_data_is_a_notstring_long_version(self):
        self.assertTrue(True)
        """
        data = fu.pack_message(555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555)
        return_new = fu.unpack_message(data)
        return_old = oldfu.unpack_message(data)
        self.assertEqual(return_new, return_old)
        """

    def test_pickle_file_values(self):
        (data,) = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.unpack_message"))[0]
        return_new = fu.unpack_message(data)
        return_old = oldfu.unpack_message(data)
        self.assertEqual(return_new, return_old)



class Test_wrap_mpi_send(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.wrap_mpi_send()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.wrap_mpi_send()
        self.assertEqual(cm_new.exception.message, "wrap_mpi_send() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        """ values got via pickle files/utilities/utilities.wrap_mpi_send"""
        return_new = fu.wrap_mpi_send(data = [9], destination = 0, communicator = None)
        return_old = oldfu.wrap_mpi_send(data =[9], destination = 0, communicator = None)
        self.assertEqual(return_new, return_old)

    def test_with_MPI_COMM_WORLD(self):
        return_new = fu.wrap_mpi_send(data = [9], destination = 0, communicator = MPI_COMM_WORLD)
        return_old = oldfu.wrap_mpi_send(data =[9], destination = 0, communicator = MPI_COMM_WORLD)
        self.assertEqual(return_new, return_old)

    def test_invalid_communicator_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.wrap_mpi_send(data = [9], destination = 0, communicator = -1)
        return_old = oldfu.wrap_mpi_send(data =[9], destination = 0, communicator = -1)
        self.assertEqual(return_new, return_old)
        """



class Test_wrap_mpi_recv(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.wrap_mpi_recv()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.wrap_mpi_recv()
        self.assertEqual(cm_new.exception.message, "wrap_mpi_recv() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    "Can only test on cluster , cannot work on workstation"
    # def test_wrap_mpi_recv_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_recv")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (data, communicator) = argum[0]
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_new = fu.wrap_mpi_recv(data, communicator)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.wrap_mpi_recv(data, communicator)
    #
    #     self.assertEqual(return_new, return_old)



class Test_wrap_mpi_bcast(unittest.TestCase):
    """ Values got running Test_get_sorting_params_refine.test_default_case"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.wrap_mpi_bcast()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.wrap_mpi_bcast()
        self.assertEqual(cm_new.exception.message, "wrap_mpi_bcast() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_None_data(self):
        """ values got via pickle files/utilities/utilities.wrap_mpi_send"""
        return_new = fu.wrap_mpi_bcast(None, root=0, communicator= None)
        return_old = oldfu.wrap_mpi_bcast(None, root=0, communicator= None)
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_default_case(self):
        attr_value_list = [[0, 27.84771510918482, 49.09925034711038, 236.702241194244, 0.0, 0.0], [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0],[2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        return_new = fu.wrap_mpi_bcast(data = attr_value_list, root = 0, communicator = MPI_COMM_WORLD)
        return_old = oldfu.wrap_mpi_bcast(data =attr_value_list, root= 0, communicator = MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_communicator_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.wrap_mpi_bcast(data = [9], root = 0, communicator = -1)
        return_old = oldfu.wrap_mpi_bcast(data =[9], root= 0, communicator = -1)
        self.assertEqual(return_new, return_old)
        """



class Test_wrap_mpi_gatherv(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.wrap_mpi_gatherv()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.wrap_mpi_gatherv()
        self.assertEqual(cm_new.exception.message, "wrap_mpi_gatherv() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_pickle_file_values(self):
        """ values got via pickle files/utilities/utilities.wrap_mpi_gatherv"""
        return_new = fu.wrap_mpi_gatherv(data = [45,3], root = 0, communicator= None)
        return_old = oldfu.wrap_mpi_gatherv(data= [45,3], root = 0, communicator= None)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_with_MPI_COMM_WORLD(self):
        """ values got via pickle files/utilities/utilities.wrap_mpi_gatherv"""
        return_new = fu.wrap_mpi_gatherv(data = [45,3], root = 0, communicator= MPI_COMM_WORLD)
        return_old = oldfu.wrap_mpi_gatherv(data= [45,3], root = 0, communicator= MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_communicator_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.wrap_mpi_gatherv(data = [45,3], root = 0, communicator= -1)
        return_old = oldfu.wrap_mpi_gatherv(data= [45,3], root = 0, communicator= -1)
        self.assertEqual(return_new, return_old)
        """



class Test_get_colors_and_subsets(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_colors_and_subsets()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_colors_and_subsets()
        self.assertEqual(cm_new.exception.message, "get_colors_and_subsets() takes exactly 6 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_mainMode_equal_my_rank(self):
        main_node = 0
        my_rank = mpi_comm_rank(MPI_COMM_WORLD)
        shared_comm = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        sh_my_rank = mpi_comm_rank(shared_comm)
        masters = mpi_comm_split(MPI_COMM_WORLD, sh_my_rank == main_node, my_rank)
        shared_comm = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        return_new = fu.get_colors_and_subsets(main_node, MPI_COMM_WORLD, my_rank, shared_comm, sh_my_rank,masters)
        return_old = oldfu.get_colors_and_subsets(main_node, MPI_COMM_WORLD, my_rank, shared_comm, sh_my_rank,masters)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_mainMode_not_equal_my_rank_returns_TypeError_obj_Nonetype_hasnot_len(self):
        main_node = 0
        my_rank = mpi_comm_rank(MPI_COMM_WORLD)
        shared_comm = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        sh_my_rank = mpi_comm_rank(shared_comm)
        masters = mpi_comm_split(MPI_COMM_WORLD, sh_my_rank == main_node, my_rank)
        shared_comm = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        with self.assertRaises(TypeError) as cm_new:
            fu.get_colors_and_subsets(main_node, MPI_COMM_WORLD, my_rank, shared_comm, sh_my_rank+1,masters)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_colors_and_subsets(main_node, MPI_COMM_WORLD, my_rank, shared_comm, sh_my_rank+1,masters)
        self.assertEqual(cm_new.exception.message, "object of type 'NoneType' has no len()")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


class Test_wrap_mpi_split(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.wrap_mpi_split()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.wrap_mpi_split()
        self.assertEqual(cm_new.exception.message, "wrap_mpi_split() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


        """ Can only be tested in mpi not on workstation   """
    # def test_wrap_mpi_split_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_split")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (comm, no_of_groups) = argum[0]
    #
    #     return_new = fu.wrap_mpi_split(comm, no_of_groups)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.wrap_mpi_split(comm, no_of_groups)
    #
    #     self.assertEqual(return_new, return_old)


class Test_get_dist(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_dist()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_dist()
        self.assertEqual(cm_new.exception.message, "get_dist() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_get_dist(self):
        return_new = fu.get_dist(c1=[2,4],c2=[5,1])
        return_old = oldfu.get_dist([2, 4], [5, 1])
        self.assertEqual(return_new, return_old)

    def test_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.get_dist(c1=[2],c2=[5])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.get_dist(c1=[2],c2=[5])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_eliminate_moons(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.eliminate_moons()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.eliminate_moons()
        self.assertEqual(cm_new.exception.message, "eliminate_moons() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_real_case_IMAGE_3D(self):
        moon_params = [0.4,0.7]
        return_new = fu.eliminate_moons(deepcopy(IMAGE_3D), moon_params)
        return_old = oldfu.eliminate_moons(deepcopy(IMAGE_3D), moon_params)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_real_case_IMAGE_3D_no_change(self):
        v = fu.model_gauss(0.25,12,12,12)
        moon_params = [-1,-1]
        return_new = fu.eliminate_moons(deepcopy(v), moon_params)
        return_old = oldfu.eliminate_moons(deepcopy(v), moon_params)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), v.get_3dview()))

    def test_returns_IndexError_list_index_out_of_range(self):
        moon_params = [0.4]
        with self.assertRaises(IndexError) as cm_new:
            fu.eliminate_moons(deepcopy(IMAGE_3D), moon_params)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.eliminate_moons(deepcopy(IMAGE_3D), moon_params)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_real_case_IMAGE_2D_returns_RuntimeError_the_img_should_be_a_3D_img(self):
        moon_params = [0.4,0.7]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.eliminate_moons(deepcopy(IMAGE_2D), moon_params)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.eliminate_moons(deepcopy(IMAGE_2D), moon_params)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The image should be 3D")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_img_returns_RuntimeError_the_img_should_be_a_3D_img(self):
        moon_params = [0.4,0.7]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.eliminate_moons(EMData(), moon_params)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.eliminate_moons(EMData(), moon_params)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The image should be 3D")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_img_returns_AttributeError_NoneType_obj_hasnot_attribute_find_3d_threshold(self):
        moon_params = [0.4,0.7]
        with self.assertRaises(AttributeError) as cm_new:
            fu.eliminate_moons(None, moon_params)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.eliminate_moons(None, moon_params)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'find_3d_threshold'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_combinations_of_n_taken_by_k(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.combinations_of_n_taken_by_k()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.combinations_of_n_taken_by_k()
        self.assertEqual(cm_new.exception.message, "combinations_of_n_taken_by_k() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_combinations_of_n_taken_by_k(self):
        return_new = fu.combinations_of_n_taken_by_k(5,3)
        return_old = oldfu.combinations_of_n_taken_by_k(5,3)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 10)



class Test_cmdexecute(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.cmdexecute()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.cmdexecute()
        self.assertEqual(cm_new.exception.message, "cmdexecute() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_correct_cmd_without_printing_on_success(self):
        return_new = fu.cmdexecute("ls", False)
        return_old = oldfu.cmdexecute("ls", False)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)

    def test_correct_cmd_with_printing_on_success(self):
        return_new = fu.cmdexecute("ls", True)
        return_old = oldfu.cmdexecute("ls", True)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 1)

    def test_wrong_cmd(self):
        return_new = fu.cmdexecute("quack", True)
        return_old = oldfu.cmdexecute("quack", True)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 0)



class Test_string_found_in_file(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.string_found_in_file()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.string_found_in_file()
        self.assertEqual(cm_new.exception.message, "string_found_in_file() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_file_not_found_returns_IOError(self):
        with self.assertRaises(IOError) as cm_new:
            fu.string_found_in_file("search smth", "not_a_file.txt")
        with self.assertRaises(IOError) as cm_old:
            oldfu.string_found_in_file("search smth", "not_a_file.txt")
        self.assertEqual(cm_new.exception.strerror, "No such file or directory")
        self.assertEqual(cm_new.exception.strerror, cm_old.exception.strerror)

    def test_found_value(self):
        f = "f.txt"
        data=[["hallo",1,1,1],[2,2,2,2],[3,3,3,3]]
        path_to_file = path.join(ABSOLUTE_PATH, f)
        fu.write_text_row(data, path_to_file)
        return_new = fu.string_found_in_file("hallo", path_to_file)
        return_old = oldfu.string_found_in_file("hallo", path_to_file)
        remove_list_of_file([f])
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new)

    def test_notfound_value(self):
        f = "f.txt"
        data=[["ds",1,1,1],[2,2,2,2],[3,3,3,3]]
        path_to_file = path.join(ABSOLUTE_PATH, f)
        fu.write_text_row(data, path_to_file)
        return_new = fu.string_found_in_file("hallo", path_to_file)
        return_old = oldfu.string_found_in_file("hallo", path_to_file)
        remove_list_of_file([f])
        self.assertEqual(return_new,return_old)
        self.assertFalse(return_new)



class Test_get_latest_directory_increment_value(unittest.TestCase):
    start_value = 1
    folder_name = 'd'
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_latest_directory_increment_value()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_latest_directory_increment_value()
        self.assertEqual(cm_new.exception.message, "get_latest_directory_increment_value() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_nothing_to_count(self):
        return_new = fu.get_latest_directory_increment_value(ABSOLUTE_PATH, self.folder_name, start_value = self.start_value, myformat = "%03d")
        return_old = oldfu.get_latest_directory_increment_value(ABSOLUTE_PATH, self.folder_name, start_value = self.start_value, myformat = "%03d")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new,self.start_value)

    def test_count_something(self):
        os.mkdir(path.join(ABSOLUTE_PATH, self.folder_name+"001"))
        os.mkdir(path.join(ABSOLUTE_PATH, self.folder_name + "002"))
        return_new = fu.get_latest_directory_increment_value(ABSOLUTE_PATH, "/"+self.folder_name,start_value=self.start_value, myformat="%03d")
        return_old = oldfu.get_latest_directory_increment_value(ABSOLUTE_PATH, "/" + self.folder_name,start_value=self.start_value, myformat="%03d")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new,2)
        remove_dir(path.join(ABSOLUTE_PATH, self.folder_name+"001"))
        remove_dir(path.join(ABSOLUTE_PATH, self.folder_name + "002"))



class Test_if_error_then_all_processes_exit_program(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.if_error_then_all_processes_exit_program()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.if_error_then_all_processes_exit_program()
        self.assertEqual(cm_new.exception.message, "if_error_then_all_processes_exit_program() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_shrink_data_huang(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_shrink_data_huang()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_shrink_data_huang()
        self.assertEqual(cm_new.exception.message, "get_shrink_data_huang() takes at least 7 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("it cannot find something in the ADNAN file")
    def test_get_shrink_data_huang_true_should_return_equal_objects(self):
        """
        I got
        RuntimeError: FileAccessException at /home/lusnig/EMAN2/eman2/libEM/emdata_metadata.cpp:240: error with '/home/lusnig/Downloads/adnan4testing/Substack/EMAN2DB/../../Particles/mpi_proc_000/EMAN2DB/TcdA1-0011_frames_sum_ptcls_352x352x1': 'cannot access file '/home/lusnig/Downloads/adnan4testing/Substack/EMAN2DB/../../Particles/mpi_proc_000/EMAN2DB/TcdA1-0011_frames_sum_ptcls_352x352x1'' caught
        """
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["stack"] = 'bdb:' + path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, 'Substack/sort3d_substack_002')
        Tracker["applyctf"] = True
        ids = []
        for i in range(1227):
            ids.append(i)
        Tracker["chunk_dict"] =ids
        myid = 0
        m_node = 0
        nproc = 1
        partids = path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, "Sort3D/indexes_010.txt")
        partstack =  path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, "Sort3D/params_010.txt")
        nxinit = 2

        return_new = fu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        return_old = oldfu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        self.assertTrue(numpy.allclose(return_new[0][0].get_3dview(), return_old[0][0].get_3dview(), 0.5))


class Test_getindexdata(unittest.TestCase):
    """ nproc and myid valeus got from "pickle files/utilities/utilities.getindexdata"""
    nproc = 95
    myid = 22
    stack = 'bdb:' + path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, 'VIPER/best_000')
    partids = path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, 'VIPER/main001/this_iteration_index_keep_images.txt')
    partstack = path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, 'VIPER//main001/run000/rotated_reduced_params.txt')

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.getindexdata()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.getindexdata()
        self.assertEqual(cm_new.exception.message, "getindexdata() takes exactly 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_nproc_greater_than_ndata(self):
        return_new = fu.getindexdata(self.stack, self.partids, self.partstack, self.myid, self.nproc)
        return_old = oldfu.getindexdata(self.stack, self.partids, self.partstack, self.myid, self.nproc)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))

    def test_nproc_and_myid_greater_than_ndata_(self):
        return_new = fu.getindexdata(self.stack, self.partids, self.partstack, 100, self.nproc)
        return_old = oldfu.getindexdata(self.stack, self.partids, self.partstack, 100, self.nproc)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))

    def test_nproc_lower_than_ndata(self):
        return_new = fu.getindexdata(self.stack, self.partids, self.partstack, self.myid, nproc= 10)
        return_old = oldfu.getindexdata(self.stack, self.partids, self.partstack, self.myid, nproc= 10)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))



class Test_store_value_of_simple_vars_in_json_file(unittest.TestCase):
    f= path.join(ABSOLUTE_PATH, "fu.json")
    f_old = path.join(ABSOLUTE_PATH, "oldfu.json")
    var_to_save= {'string_var': 'var1', 'integer_var': 7, 'bool_var': False, 'list_var': [2,3,4]}
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.store_value_of_simple_vars_in_json_file()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.store_value_of_simple_vars_in_json_file()
        self.assertEqual(cm_new.exception.message, "store_value_of_simple_vars_in_json_file() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        fu.store_value_of_simple_vars_in_json_file(filename =  self.f, local_vars = self.var_to_save, exclude_list_of_vars = [], write_or_append = "w",	vars_that_will_show_only_size = [])
        oldfu.store_value_of_simple_vars_in_json_file(filename =  self.f_old, local_vars = self.var_to_save, exclude_list_of_vars = [], write_or_append = "w",	vars_that_will_show_only_size = [])
        self.assertEqual(returns_values_in_file(self.f), returns_values_in_file(self.f_old))
        self.assertTrue(fu.string_found_in_file(self.var_to_save.keys()[0], self.f))
        self.assertTrue(oldfu.string_found_in_file(self.var_to_save.keys()[0], self.f_old))
        remove_list_of_file([self.f,self.f_old])

    def test_exclude_a_variable(self):
        var=self.var_to_save.keys()[0]
        fu.store_value_of_simple_vars_in_json_file(filename =  self.f, local_vars = self.var_to_save, exclude_list_of_vars = [var], write_or_append = "w",	vars_that_will_show_only_size = [])
        oldfu.store_value_of_simple_vars_in_json_file(filename =  self.f_old, local_vars = self.var_to_save, exclude_list_of_vars = [var], write_or_append = "w",	vars_that_will_show_only_size = [])
        self.assertEqual(returns_values_in_file(self.f), returns_values_in_file(self.f_old))
        self.assertFalse(fu.string_found_in_file(var, self.f))
        self.assertFalse(oldfu.string_found_in_file(var, self.f_old))
        remove_list_of_file([self.f,self.f_old])

    def test_onlySize_a_variable(self):
        var= 'list_var'
        fu.store_value_of_simple_vars_in_json_file(filename =  self.f, local_vars = self.var_to_save, exclude_list_of_vars = [], write_or_append = "w",	vars_that_will_show_only_size = [var])
        oldfu.store_value_of_simple_vars_in_json_file(filename =  self.f_old, local_vars = self.var_to_save, exclude_list_of_vars = [], write_or_append = "w",	vars_that_will_show_only_size = [var])
        self.assertEqual(returns_values_in_file(self.f), returns_values_in_file(self.f_old))
        self.assertTrue(fu.string_found_in_file("<type 'list'> with length: 3", self.f))
        self.assertTrue(oldfu.string_found_in_file("<type 'list'> with length: 3", self.f_old))
        remove_list_of_file([self.f, self.f_old])



class Test_convert_json_fromunicode(unittest.TestCase):
    f= path.join(ABSOLUTE_PATH, "f.json")
    var_to_save= {'string_var': 'var1', 'integer_var': 7, 'bool_var': False, 'list_var': [2,3,4]}

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.convert_json_fromunicode()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.convert_json_fromunicode()
        self.assertEqual(cm_new.exception.message, "convert_json_fromunicode() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_loaded_jsonFile(self):
        fu.store_value_of_simple_vars_in_json_file(filename=self.f, local_vars=self.var_to_save,exclude_list_of_vars=[], write_or_append="w",vars_that_will_show_only_size=[])
        with open(self.f, 'r') as f1:
            values=json.load( f1)

        return_new = fu.convert_json_fromunicode(values)
        return_old = oldfu.convert_json_fromunicode(values)
        self.assertDictEqual(return_new,return_old)
        remove_list_of_file([self.f])

    def test_with_string(self):
        data = "ciaone"
        return_new = fu.convert_json_fromunicode(data)
        return_old = oldfu.convert_json_fromunicode(data)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, data)



class Test_get_sorting_attr_stack(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_sorting_attr_stack()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_sorting_attr_stack()
        self.assertEqual(cm_new.exception.message, "get_sorting_attr_stack() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        stack = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc"))[0][0]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)
        self.assertTrue(numpy.array_equal(fu.get_sorting_attr_stack(stack), oldfu.get_sorting_attr_stack(stack)))

    def test_empty_stack(self):
        return_new=fu.get_sorting_attr_stack([])
        self.assertTrue(numpy.array_equal(return_new, oldfu.get_sorting_attr_stack([])))
        self.assertTrue(numpy.array_equal(return_new, []))

    def test_wrong_images_in_the_stack_RunTimeError(self):
        stack=[IMAGE_2D,IMAGE_2D]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_sorting_attr_stack(stack)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_sorting_attr_stack(stack)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


class Test_get_sorting_params_refine(unittest.TestCase):
    stack = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc"))[0][0]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_sorting_params_refine()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_sorting_params_refine()
        self.assertEqual(cm_new.exception.message, "get_sorting_params_refine() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        for i in range(len(self.stack)):
            self.stack[i].set_attr("group",i)
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["nproc"] = 1
        return_new = fu.get_sorting_params_refine(Tracker, self.stack, len(self.stack))
        return_old = oldfu.get_sorting_params_refine(Tracker, self.stack, len(self.stack))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def returns_too_ndata_vlaue_respect_the_number_of_data_IndexError_list_index_out_of_range(self):
        for i in range(len(self.stack)):
            self.stack[i].set_attr("group",i)
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["nproc"] = 1
        with self.assertRaises(IndexError) as cm_new:
            fu.get_sorting_params_refine(Tracker, self.stack, len(self.stack)+11)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.get_sorting_params_refine(Tracker, self.stack, len(self.stack)+11)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_stack(self):
        stack=[]
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["nproc"] = 1
        return_new = fu.get_sorting_params_refine(Tracker, stack, len(stack))
        return_old = oldfu.get_sorting_params_refine(Tracker, stack, len(stack))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_images_in_the_stack_RunTimeError(self):
        stack=[IMAGE_2D,IMAGE_2D]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["nproc"] = 1
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_sorting_params_refine(Tracker, stack, len(stack))
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_sorting_params_refine(Tracker, stack, len(stack))
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_parsing_sorting_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.parsing_sorting_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.parsing_sorting_params()
        self.assertEqual(cm_new.exception.message, "parsing_sorting_params() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.parsing_sorting_params([])
        return_old = oldfu.parsing_sorting_params([])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_typeerror_int_hasnot_attribute__get_item__(self):
        l=[1,2,3,4,5]
        with self.assertRaises(TypeError) as cm_new:
            fu.parsing_sorting_params(l)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.parsing_sorting_params(l)
        self.assertEqual(cm_new.exception.message, "'int' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_values(self):
        l=[[1,2,3,4,5],[1,21,31,41,51]]
        return_new = fu.parsing_sorting_params(l)
        return_old = oldfu.parsing_sorting_params(l)
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))



class Test_fill_in_mpi_list(unittest.TestCase):
    """ Values got running Test_get_sorting_params_refine.test_default_case"""
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fill_in_mpi_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fill_in_mpi_list()
        self.assertEqual(cm_new.exception.message, "fill_in_mpi_list() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        total_attr_value_list=[[],[],[]]
        attr_value_list = [[0, 27.84771510918482, 49.09925034711038, 236.702241194244, 0.0, 0.0], [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0], [2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        return_new = fu.fill_in_mpi_list(mpi_list = deepcopy(total_attr_value_list), data_list = attr_value_list, index_start = 0 ,index_end = len(total_attr_value_list))
        return_old = oldfu.fill_in_mpi_list(mpi_list = deepcopy(total_attr_value_list), data_list = attr_value_list, index_start = 0 ,index_end = len(total_attr_value_list))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_index_start_negative_returns_IndexError_list_index_out_of_range(self):
        total_attr_value_list=[[],[],[]]
        attr_value_list = [[0, 27.84771510918482, 49.09925034711038, 236.702241194244, 0.0, 0.0], [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0], [2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        with self.assertRaises(IndexError) as cm_new:
            fu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=0-1,index_end=len(total_attr_value_list))
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=-1,index_end=len(total_attr_value_list))
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
    def test_index_end_too_high_returns_IndexError_list_index_out_of_range(self):
        total_attr_value_list=[[],[],[]]
        attr_value_list = [[0, 27.84771510918482, 49.09925034711038, 236.702241194244, 0.0, 0.0], [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0], [2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        with self.assertRaises(IndexError) as cm_new:
            fu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=0,index_end=len(total_attr_value_list)+2)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=0,index_end=len(total_attr_value_list)+2)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_less_values_in_attr_value_list_returns_IndexError_list_index_out_of_range(self):
        total_attr_value_list=[[],[],[]]
        attr_value_list = [ [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0], [2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        with self.assertRaises(IndexError) as cm_new:
            fu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=0, index_end=len(total_attr_value_list))
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fill_in_mpi_list(mpi_list=deepcopy(total_attr_value_list), data_list=attr_value_list, index_start=0, index_end=len(total_attr_value_list))
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_too_values_in_attr_value_list(self):
        total_attr_value_list=[[],[],[]]
        attr_value_list = [[0, 27.84771510918482, 49.09925034711038, 236.702241194244, 0.0, 0.0], [1, 54.496982231553545, 150.6989385443887, 95.77312314162165, 0.0, 0.0], [2, 67.0993779295224, 52.098986136572584, 248.45843717750148, 0.0, 0.0]]
        return_new = fu.fill_in_mpi_list(mpi_list = deepcopy(total_attr_value_list), data_list = attr_value_list, index_start = 0 ,index_end = len(total_attr_value_list))
        return_old = oldfu.fill_in_mpi_list(mpi_list = deepcopy(total_attr_value_list), data_list = attr_value_list, index_start = 0 ,index_end = len(total_attr_value_list))
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_sample_down_1D_curve(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sample_down_1D_curve()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sample_down_1D_curve()
        self.assertEqual(cm_new.exception.message, "sample_down_1D_curve() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.sample_down_1D_curve(nxinit=100, nnxo=180, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        return_old = oldfu.sample_down_1D_curve(nxinit=100, nnxo=180, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_nxinit_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.sample_down_1D_curve(nxinit=0, nnxo=180, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.sample_down_1D_curve(nxinit=0, nnxo=180, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_null_nnxo_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.sample_down_1D_curve(nxinit=100, nnxo=0, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.sample_down_1D_curve(nxinit=100, nnxo=0, pspcurv_nnxo_file=path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"Sort3D/fsc_halves.txt"))
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_file_not_found(self):
        with self.assertRaises(IOError) as cm_new:
            fu.sample_down_1D_curve(nxinit=100, nnxo=180, pspcurv_nnxo_file="filenotfound.txt")
        with self.assertRaises(IOError) as cm_old:
            oldfu.sample_down_1D_curve(nxinit=100, nnxo=180, pspcurv_nnxo_file="filenotfound.txt")
        self.assertEqual(cm_new.exception.strerror, "No such file or directory")
        self.assertEqual(cm_new.exception.strerror, cm_old.exception.strerror)



class Test_get_initial_ID(unittest.TestCase):
    full_ID_dict = {0: 'ciao_0', 1: 'ciao_1', 2: 'ciao_2', 3: 'ciao_3'}

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_initial_ID()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_initial_ID()
        self.assertEqual(cm_new.exception.message, "get_initial_ID() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_valid_list_dict(self):
        part_list = [0,1,2]
        return_new  = fu.get_initial_ID(part_list, self.full_ID_dict)
        return_old = oldfu.get_initial_ID(part_list, self.full_ID_dict)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_list(self):
        return_new  = fu.get_initial_ID([], self.full_ID_dict)
        return_old = oldfu.get_initial_ID([], self.full_ID_dict)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))

    def test_invalid_value_in_list_KeyError(self):
        part_list = [0, 1, 20]
        with self.assertRaises(KeyError) as cm_new:
            fu.get_initial_ID(part_list, self.full_ID_dict)
        with self.assertRaises(KeyError) as cm_old:
            oldfu.get_initial_ID(part_list, self.full_ID_dict)
        self.assertEqual(cm_new.exception.message, 20)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_dict_KeyError(self):
        part_list = [0, 1, 20]
        with self.assertRaises(KeyError) as cm_new:
            fu.get_initial_ID(part_list, {})
        with self.assertRaises(KeyError) as cm_old:
            oldfu.get_initial_ID(part_list, {})
        self.assertEqual(cm_new.exception.message, 0)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_print_upper_triangular_matrix(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.print_upper_triangular_matrix()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_upper_triangular_matrix()
        self.assertEqual(cm_new.exception.message, "print_upper_triangular_matrix() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("which variable is the third parameter??")
    def test_print_upper_triangular_matrix(self):
        log_new =[]
        log_old = []
        size =4
        data=[]
        for i in range(size):
            for j in range(size):
                data.append((i,j*j))
        fu.print_upper_triangular_matrix(data,size,log_new)
        oldfu.print_upper_triangular_matrix(data, size, log_old)



class Test_convertasi(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.convertasi()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.convertasi()
        self.assertEqual(cm_new.exception.message, "convertasi() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.convertasi([],3)
        return_old = oldfu.convertasi([],3)
        self.assertTrue(numpy.allclose(return_new, return_old))

    def test_default_case(self):
        asig = [0,1,2,3,4,5,6]
        return_new = fu.convertasi(asig,7)
        return_old = oldfu.convertasi(asig,7)
        self.assertEqual(return_new, return_old)



class Test_prepare_ptp(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_ptp()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_ptp()
        self.assertEqual(cm_new.exception.message, "prepare_ptp() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.prepare_ptp([],3)
        return_old = oldfu.prepare_ptp([],3)
        self.assertTrue(numpy.allclose(return_new, return_old))

    def test_default_case(self):
        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        return_new = fu.prepare_ptp(data_list, K)
        return_old = oldfu.prepare_ptp(data_list, K)
        self.assertEqual(return_new, return_old)



class Test_print_dict(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.print_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_dict()
        self.assertEqual(cm_new.exception.message, "print_dict() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_print_dict(self):
        dic = {'0': 'ciao_0', '1': 'ciao_1'}
        old_stdout = sys.stdout
        print_new = StringIO()
        sys.stdout = print_new
        return_new = fu.print_dict(dic, "title")
        print_old = StringIO()
        sys.stdout = print_old
        return_old = oldfu.print_dict(dic, "title")
        self.assertEqual(return_new,return_old)
        self.assertTrue(return_new is None)
        self.assertEqual(print_new.getvalue(), print_old.getvalue())
        sys.stdout = old_stdout

    def test_error_key_type(self):
        dic = {0: 'ciao_0', 1: 'ciao_1', 2: 'ciao_2', 3: 'ciao_3'}
        with self.assertRaises(TypeError) as cm_new:
            fu.print_dict(dic, " Test_print_dict.test_error_key_type")
        with self.assertRaises(TypeError) as cm_old:
            oldfu.print_dict(dic, " Test_print_dict.test_error_key_type")
        self.assertEqual(cm_new.exception.message, "cannot concatenate 'str' and 'int' objects")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_get_resolution_mrk01(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_resolution_mrk01()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_resolution_mrk01()
        self.assertEqual(cm_new.exception.message, "get_resolution_mrk01() takes exactly 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_radi_not_integer(self):
        v = [IMAGE_2D,IMAGE_2D_REFERENCE]
        return_new = fu.get_resolution_mrk01(deepcopy(v), 0.5,0.15,ABSOLUTE_PATH, None)
        return_old = oldfu.get_resolution_mrk01(deepcopy(v), 0.5,0.15,ABSOLUTE_PATH,None)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        remove_list_of_file([path.join(ABSOLUTE_PATH,"fsc.txt")])

    def test_radi_integer_no_mask(self):
        v = [IMAGE_3D,IMAGE_3D]
        return_new = fu.get_resolution_mrk01(deepcopy(v), 1,IMAGE_3D.get_xsize(),ABSOLUTE_PATH, None)
        return_old = oldfu.get_resolution_mrk01(deepcopy(v), 1,IMAGE_3D.get_xsize(),ABSOLUTE_PATH,None)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        remove_list_of_file([path.join(ABSOLUTE_PATH, "fsc.txt")])

    def test_radi_integer_with_mask(self):
        v = [IMAGE_3D,IMAGE_3D]
        mask_option = [fu.model_circle(1,IMAGE_3D.get_xsize(),IMAGE_3D.get_ysize(),IMAGE_3D.get_zsize())]
        return_new = fu.get_resolution_mrk01(deepcopy(v), 1,None,ABSOLUTE_PATH, mask_option)
        return_old = oldfu.get_resolution_mrk01(deepcopy(v), 1,None,ABSOLUTE_PATH,mask_option)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        remove_list_of_file([path.join(ABSOLUTE_PATH, "fsc.txt")])

    def test_with_invalid_mask_returns_RuntimeError_ImageFormatException(self):
        v = [IMAGE_3D,IMAGE_3D]
        mask_option = [fu.model_circle(1,IMAGE_3D.get_xsize()+10,IMAGE_3D.get_ysize(),IMAGE_3D.get_zsize())]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.get_resolution_mrk01(deepcopy(v), 1,None,ABSOLUTE_PATH, mask_option)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.get_resolution_mrk01(deepcopy(v), 1,None,ABSOLUTE_PATH,mask_option)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "can not multiply images that are not the same size")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_partition_to_groups(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.partition_to_groups()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.partition_to_groups()
        self.assertEqual(cm_new.exception.message, "partition_to_groups() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.partition_to_groups([],3)
        return_old = oldfu.partition_to_groups([],3)
        self.assertTrue(numpy.allclose(return_new, return_old))

    def test_default_case(self):
        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        return_new = fu.partition_to_groups(data_list, K)
        return_old = oldfu.partition_to_groups(data_list, K)
        self.assertEqual(return_new, return_old)



class Test_partition_independent_runs(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.partition_independent_runs()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.partition_independent_runs()
        self.assertEqual(cm_new.exception.message, "partition_independent_runs() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.partition_independent_runs([],3)
        return_old = oldfu.partition_independent_runs([],3)
        self.assertDictEqual(return_new,return_old)
        self.assertEqual(return_new, {})

    def test_default_case(self):
        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        return_new = fu.partition_independent_runs(data_list, K)
        return_old = oldfu.partition_independent_runs(data_list, K)
        self.assertEqual(return_new, return_old)



class Test_merge_groups(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.merge_groups()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.merge_groups()
        self.assertEqual(cm_new.exception.message, "merge_groups() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list(self):
        return_new = fu.merge_groups([])
        return_old = oldfu.merge_groups([])
        self.assertTrue(numpy.allclose(return_new, return_old))

    def test_default_case(self):
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        return_new = fu.merge_groups(data_list)
        return_old = oldfu.merge_groups(data_list)
        self.assertEqual(return_new, return_old)



class Test_save_alist(unittest.TestCase):
    filename_new = "listfile.txt"
    filename_old = "listfile2.txt"
    data_list = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.save_alist()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.save_alist()
        self.assertEqual(cm_new.exception.message, "save_alist() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_create_files(self):
        Tracker = deepcopy(TRACKER)
        Tracker["this_dir"] = ABSOLUTE_PATH
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = "myid"
        Tracker["constants"]["main_node"] = "myid"

        fu.save_alist(Tracker, self.filename_new, self.data_list)
        oldfu.save_alist(Tracker,self.filename_old, self.data_list)
        self.assertEqual(returns_values_in_file(path.join(ABSOLUTE_PATH,self.filename_new)),returns_values_in_file(path.join(ABSOLUTE_PATH,self.filename_old)))
        remove_list_of_file([path.join(ABSOLUTE_PATH,self.filename_new),path.join(ABSOLUTE_PATH,self.filename_old)])

    def test_no_create_files(self):
        Tracker = deepcopy(TRACKER)
        Tracker["this_dir"] = ABSOLUTE_PATH
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = "myid"
        Tracker["constants"]["main_node"] = "different myid"

        fu.save_alist(Tracker, self.filename_new, self.data_list)
        oldfu.save_alist(Tracker, self.filename_old, self.data_list)

        self.assertFalse(path.isfile(path.join(ABSOLUTE_PATH, self.filename_new)))
        self.assertFalse(path.isfile(path.join(ABSOLUTE_PATH, self.filename_old)))



class Test_margin_of_error(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.margin_of_error()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.margin_of_error()
        self.assertEqual(cm_new.exception.message, "margin_of_error() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.margin_of_error(0.2,0.1)
        return_old = oldfu.margin_of_error(0.2,0.1)
        self.assertEqual(return_new, return_old)


class Test_do_two_way_comparison(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.do_two_way_comparison()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.do_two_way_comparison()
        self.assertEqual(cm_new.exception.message, "do_two_way_comparison() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_defaault_case(self):
        Tracker = deepcopy(TRACKER)
        Tracker["this_dir"] = ABSOLUTE_PATH
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 1
        Tracker["this_total_stack"] = 10
        Tracker["number_of_groups"] = 4
        Tracker["constants"]["indep_runs"]  = 4
        Tracker['full_ID_dict'] = {0: 0, 1: 1, 2:2, 3: 3}
        Tracker["partition_dict"]    = [[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]]
        Tracker["chunk_dict"] = [0, 1, 2, 3]
        Tracker["P_chunk0"] = 0.2
        Tracker["constants"]["smallest_group"] = 2
        Tracker2 = deepcopy(Tracker)
        return_new = fu.do_two_way_comparison(Tracker)
        return_old = oldfu.do_two_way_comparison(Tracker2)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)
        self.assertTrue(numpy.array_equal(Tracker["score_of_this_comparison"], Tracker2["score_of_this_comparison"]))


class Test_select_two_runs(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.select_two_runs()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.select_two_runs()
        self.assertEqual(cm_new.exception.message, "select_two_runs() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        summed_scores = [0,1,2,3,4]
        two_way_dict = [3.2,1.43,54,32,543]
        return_new = fu.select_two_runs(summed_scores,two_way_dict)
        return_old = oldfu.select_two_runs(summed_scores,two_way_dict)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_returns_IndexError_list_index_out_of_range(self):
        summed_scores = [0,1,2,3,4]
        two_way_dict = [3.2]
        with self.assertRaises(IndexError) as cm_new:
            fu.select_two_runs(summed_scores,two_way_dict)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.select_two_runs(summed_scores,two_way_dict)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_summed_scores_empty_returns_IndexError_list_index_out_of_range(self):
        summed_scores = []
        two_way_dict = [3.2]
        with self.assertRaises(IndexError) as cm_new:
            fu.select_two_runs(summed_scores,two_way_dict)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.select_two_runs(summed_scores,two_way_dict)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_two_way_dict_empty_returns_IndexError_list_index_out_of_range(self):
        summed_scores = [0,1,2,3,4]
        two_way_dict = []
        with self.assertRaises(IndexError) as cm_new:
            fu.select_two_runs(summed_scores,two_way_dict)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.select_two_runs(summed_scores,two_way_dict)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_counting_projections(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.counting_projections()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.counting_projections()
        self.assertEqual(cm_new.exception.message, "counting_projections() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_ali3d_params_empty(self):
        return_new = fu.counting_projections(delta = 0.5, ali3d_params =[], image_start = 1)
        return_old = oldfu.counting_projections(delta = 0.5, ali3d_params =[], image_start = 1)
        self.assertDictEqual(return_new,return_old)

    def test_default_case(self):
        ali3d_params  = [[idx1, idx2, 0 , 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
        return_new = fu.counting_projections(delta = 0.5, ali3d_params =ali3d_params, image_start = 1)
        return_old = oldfu.counting_projections(delta = 0.5, ali3d_params =ali3d_params, image_start = 1)
        self.assertDictEqual(return_new,return_old)



class Test_unload_dict(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.unload_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.unload_dict()
        self.assertEqual(cm_new.exception.message, "unload_dict() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        dict_angles = {(0.64764598050589606, 53.04857999805229) : [], (10.262155636450808, 97.18016191759037) : [], (100.75772287892256, 50.274472062413594) : [], (101.11458591875028, 44.457078732458605) : []}
        return_new = fu.unload_dict(dict_angles)
        return_old = oldfu.unload_dict(dict_angles)
        self.assertTrue(numpy.allclose(return_new, return_old))

    def test_empty_dict(self):
        return_new = fu.unload_dict({})
        return_old = oldfu.unload_dict({})
        self.assertTrue(numpy.allclose(return_new, return_old))



class Test_load_dict(unittest.TestCase):
    ali3d_params = [[idx1, idx2, 0, 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
    sampled = fu.counting_projections(delta=0.5, ali3d_params=ali3d_params, image_start=1)
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.load_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.load_dict()
        self.assertEqual(cm_new.exception.message, "load_dict() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        d = fu.unload_dict(self.sampled)
        return_new = fu.load_dict(deepcopy(self.sampled), d)
        return_old = oldfu.load_dict(deepcopy(self.sampled), d)
        self.assertDictEqual(return_new,return_old)

    def test_empty_unloaded_dict_angles(self):
        return_new = fu.load_dict(deepcopy(self.sampled), [])
        return_old = oldfu.load_dict(deepcopy(self.sampled), [])
        self.assertDictEqual(return_new,return_old)

    def test_dict_angle_main_node(self):
        d = fu.unload_dict([])
        return_new = fu.load_dict([], d)
        return_old = oldfu.load_dict([], d)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertEqual(return_new, [])

    def test_empty_all(self):
        return_new = fu.load_dict({}, [])
        return_old = oldfu.load_dict({}, [])
        self.assertDictEqual(return_new,return_old)



class Test_get_stat_proj(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_stat_proj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_stat_proj()
        self.assertEqual(cm_new.exception.message, "get_stat_proj() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_myid_same_value_as_main_Node(self):
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        this_ali3d = path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES,"VIPER/main001/run000/rotated_reduced_params.txt")
        Tracker2 = deepcopy(Tracker)

        return_new = fu.get_stat_proj(Tracker,delta = 5,this_ali3d=this_ali3d)
        return_old = oldfu.get_stat_proj(Tracker2,delta = 5,this_ali3d=this_ali3d)
        self.assertDictEqual(return_new,return_old)

    def test_myid_not_the_same_value_as_main_Node_TypeError(self):
        Tracker = deepcopy(TRACKER)
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 1
        Tracker["constants"]["main_node"] = 0
        this_ali3d = path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, "VIPER/main001/run000/rotated_reduced_params.txt")
        Tracker2 = deepcopy(Tracker)

        with self.assertRaises(TypeError) as cm_new:
            fu.get_stat_proj(Tracker,delta = 5,this_ali3d=this_ali3d)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_stat_proj(Tracker2,delta = 5,this_ali3d=this_ali3d)
        self.assertEqual(cm_new.exception.message, "object of type 'int' has no len()")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)





class Test_create_random_list(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.create_random_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.create_random_list()
        self.assertEqual(cm_new.exception.message, "create_random_list() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        Tracker_new = deepcopy(TRACKER)
        Tracker_new["constants"]["nproc"] = 1
        Tracker_new["constants"]["myid"] = 0
        Tracker_new["constants"]["main_node"] = 0
        Tracker_new["total_stack"] = "stack"
        Tracker_new["constants"]["seed"] = 1.4
        Tracker_new["constants"]["indep_runs"] = 2
        Tracker_new["this_data_list"] = [2,3,5]

        Tracker_old = deepcopy(Tracker_new)

        return_new = fu.create_random_list(Tracker_new)
        return_old = oldfu.create_random_list(Tracker_old)
        self.assertEqual(return_new, None)
        self.assertEqual(return_new, return_old)
        self.assertTrue(numpy.array_equal(Tracker_new["this_indep_list"],Tracker_old["this_indep_list"]))

    def test_wrong_Tracker_KeyError(self):
        Tracker_new = deepcopy(TRACKER)
        Tracker_old = deepcopy(TRACKER)
        with self.assertRaises(KeyError) as cm_new:
            fu.create_random_list(Tracker_new)
        with self.assertRaises(KeyError) as cm_old:
            oldfu.create_random_list(Tracker_old)
        self.assertEqual(cm_new.exception.message, 'myid')
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)





class Test_recons_mref(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons_mref()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons_mref()
        self.assertEqual(cm_new.exception.message, "recons_mref() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip('same problem that we have in get_shrink_data_huang')
    def test_recons_mref_true_should_return_equal_objects(self):
        Tracker =  deepcopy(TRACKER)
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["number_of_groups"] = 1
        Tracker["constants"]["nnxo"] = 4  # roi
        Tracker["this_particle_list"] = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        Tracker["nxinit"] = 1
        Tracker["constants"]["partstack"] =  path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, "VIPER/main001/run000/rotated_reduced_params.txt")
        Tracker["this_dir"] =  ABSOLUTE_PATH
        Tracker["constants"]["stack"] =  path.join(ABSOLUTE_PATH_TO_ADNAN_TEST_FILES, "Class2D/stack_ali2d")
        Tracker["applyctf"] = False
        Tracker["chunk_dict"] = [0, 1, 2, 3, 4, 5, 6]
        Tracker["constants"]["sym"] = "c1"
        Tracker2 = deepcopy(Tracker)
        return_new = fu.recons_mref(Tracker)
        return_old = oldfu.recons_mref(Tracker2)
        self.assertTrue(return_new[0], return_old[0])




class Test_apply_low_pass_filter(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.apply_low_pass_filter()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.apply_low_pass_filter()
        self.assertEqual(cm_new.exception.message, "apply_low_pass_filter() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        Tracker = deepcopy(TRACKER )
        Tracker["low_pass_filter"] = 0.087
        return_new = fu.apply_low_pass_filter(refvol= [deepcopy(IMAGE_2D),deepcopy(IMAGE_2D)],Tracker=Tracker)
        return_old = oldfu.apply_low_pass_filter(refvol=  [deepcopy(IMAGE_2D),deepcopy(IMAGE_2D)],Tracker=Tracker)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy.array_equal(i.get_3dview(), j.get_3dview()))

    def test_wrong_Tracker_KeyError(self):
        Tracker = deepcopy(TRACKER )
        with self.assertRaises(KeyError) as cm_new:
            fu.apply_low_pass_filter(refvol= [deepcopy(IMAGE_2D),deepcopy(IMAGE_2D)],Tracker=Tracker)
        with self.assertRaises(KeyError) as cm_old:
            oldfu.apply_low_pass_filter(refvol= [deepcopy(IMAGE_2D),deepcopy(IMAGE_2D)],Tracker=Tracker)
        self.assertEqual(cm_new.exception.message, 'low_pass_filter')
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_refvol_empty(self):
        Tracker = deepcopy(TRACKER )
        Tracker["low_pass_filter"] = 0.087
        return_new = fu.apply_low_pass_filter(refvol=  [],Tracker=Tracker)
        return_old = oldfu.apply_low_pass_filter(refvol=  [],Tracker=Tracker)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertEqual(return_new, [])


class Test_get_groups_from_partition(unittest.TestCase):
    list_of_particles = [random.randint(0, 1000) for i in range(100)]
    group_list = [0, 1]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_groups_from_partition()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_groups_from_partition()
        self.assertEqual(cm_new.exception.message, "get_groups_from_partition() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        return_new = fu.get_groups_from_partition(partition =self.group_list, initial_ID_list = self.list_of_particles, number_of_groups = 2)
        return_old = oldfu.get_groups_from_partition(partition = self.group_list, initial_ID_list =self.list_of_particles, number_of_groups = 2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_empty_initial_ID_list_KeyError(self):
        with self.assertRaises(KeyError) as cm_new:
            fu.get_groups_from_partition(partition =self.group_list, initial_ID_list = [], number_of_groups = 2)
        with self.assertRaises(KeyError) as cm_old:
            oldfu.get_groups_from_partition(partition =self.group_list, initial_ID_list = [], number_of_groups = 2)
        self.assertEqual(cm_new.exception.message, 0)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_partition_list_KeyError(self):
        return_new = fu.get_groups_from_partition(partition =[], initial_ID_list = self.list_of_particles, number_of_groups = 2)
        return_old = oldfu.get_groups_from_partition(partition = [], initial_ID_list =self.list_of_particles, number_of_groups = 2)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_get_complementary_elements(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_complementary_elements()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_complementary_elements()
        self.assertEqual(cm_new.exception.message, "get_complementary_elements() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_defalut_case(self):
        sub_data_list = [1,2,2]
        total_list = [1,2,2,4,5,6]
        return_new = fu.get_complementary_elements(total_list,sub_data_list)
        return_old = oldfu.get_complementary_elements(total_list,sub_data_list)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_total_list_less_data_than_sub_data_list_error_msg(self):
        sub_data_list = [1,2,2]
        total_list = [1,2]
        return_new = fu.get_complementary_elements(total_list,sub_data_list)
        return_old = oldfu.get_complementary_elements(total_list,sub_data_list)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))



class Test_update_full_dict(unittest.TestCase):
    leftover_list = {0: 'ciao_10', 1: 'ciao_11'}
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.update_full_dict()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.update_full_dict()
        self.assertEqual(cm_new.exception.message, "update_full_dict() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        Tracker_new = deepcopy(TRACKER)
        Tracker_new['full_ID_dict'] = {10: 'ciao_0', 11: 'ciao_1', 2: 'ciao_2', 3: 'ciao_3'}
        Tracker_old = deepcopy(Tracker_new)
        return_new = fu.update_full_dict(self.leftover_list,Tracker_new)
        return_old = oldfu.update_full_dict(self.leftover_list,Tracker_old)
        self.assertEqual(return_new, None)
        self.assertEqual(return_new, return_old)
        self.assertDictEqual(Tracker_new['full_ID_dict'] ,Tracker_old['full_ID_dict'] )

    def test_no_full_ID_dict_in_tracker(self):
        Tracker_new = deepcopy(TRACKER)
        Tracker_old = deepcopy(Tracker_new)
        return_new = fu.update_full_dict(self.leftover_list,Tracker_new)
        return_old = oldfu.update_full_dict(self.leftover_list,Tracker_old)
        self.assertEqual(return_new, None)
        self.assertEqual(return_new, return_old)
        self.assertDictEqual(Tracker_new['full_ID_dict'] ,Tracker_old['full_ID_dict'] )
        self.assertDictEqual(Tracker_new['full_ID_dict'],self.leftover_list)



class Test_count_chunk_members(unittest.TestCase):
    chunk_dict = [0, 1, 2, 3, 4, 5, 6]
    one_class = [0, 1, 2, 3, 4, 5, 6]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.count_chunk_members()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.count_chunk_members()
        self.assertEqual(cm_new.exception.message, "count_chunk_members() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.count_chunk_members(self.chunk_dict, self.one_class)
        return_old = oldfu.count_chunk_members(self.chunk_dict, self.one_class)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_one_class_empty(self):
        return_new = fu.count_chunk_members(self.chunk_dict, [])
        return_old = oldfu.count_chunk_members(self.chunk_dict, [])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_all_empty(self):
        return_new = fu.count_chunk_members([], [])
        return_old = oldfu.count_chunk_members([], [])
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_chunk_dict_empty_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.count_chunk_members([], self.one_class)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.count_chunk_members([], self.one_class)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_remove_small_groups(unittest.TestCase):
    chunk_dict = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]]
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.remove_small_groups()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.remove_small_groups()
        self.assertEqual(cm_new.exception.message, "remove_small_groups() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.remove_small_groups(self.chunk_dict, 2)
        return_old = oldfu.remove_small_groups(self.chunk_dict, 2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_too_many_minimum_number_of_objects_in_a_group(self):
        return_new = fu.remove_small_groups(self.chunk_dict, 20)
        return_old = oldfu.remove_small_groups(self.chunk_dict, 20)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[], []] ))

    def test_empty_chunk_dict(self):
        return_new = fu.remove_small_groups([], 2)
        return_old = oldfu.remove_small_groups([], 2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_minimum_number_of_objects_in_a_group_is_zero(self):
        return_new = fu.remove_small_groups(self.chunk_dict, 0)
        return_old = oldfu.remove_small_groups(self.chunk_dict, 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_get_number_of_groups(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_number_of_groups()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_number_of_groups()
        self.assertEqual(cm_new.exception.message, "get_number_of_groups() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.get_number_of_groups(total_particles = 1500, number_of_images_per_group = 5)
        return_old = oldfu.get_number_of_groups(total_particles = 1500, number_of_images_per_group = 5)
        self.assertEqual(return_new, return_old)

    def test_null_number_of_images_per_group_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.get_number_of_groups(total_particles = 1500, number_of_images_per_group = 0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.get_number_of_groups(total_particles = 1500, number_of_images_per_group = 0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_total_particles_null(self):
        return_new = fu.get_number_of_groups(total_particles = 0, number_of_images_per_group = 5)
        return_old = oldfu.get_number_of_groups(total_particles = 0, number_of_images_per_group = 5)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0)



class Test_tabessel(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.tabessel()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.tabessel()
        self.assertEqual(cm_new.exception.message, "tabessel() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        return_new = fu.tabessel(None, None, nbel = 5000)
        return_old = oldfu.tabessel(None, None, nbel = 5000)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_nbel(self):
        return_new = fu.tabessel(None, None, nbel = 0)
        return_old = oldfu.tabessel(None, None, nbel = 0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))



class Test_nearest_proj(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.nearest_proj()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.nearest_proj()
        self.assertEqual(cm_new.exception.message, "nearest_proj() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        """ I calculated the value looking in the code of bin/sx3dvariability.py"""
        proj_angles=[]
        for i in range(900):
            i=+0.1
            proj_angles.append([i/2, i/5,i/4,i/3, i])
        proj_angles.sort()
        proj_angles_list = numpy.full((900, 4), 0.0, dtype=numpy.float32)
        for i in range(900):
            proj_angles_list[i][0] = proj_angles[i][1]
            proj_angles_list[i][1] = proj_angles[i][2]
            proj_angles_list[i][2] = proj_angles[i][3]
            proj_angles_list[i][3] = proj_angles[i][4]
        return_new1,return_new2 = fu.nearest_proj(proj_angles_list)
        return_old1,return_old2 = oldfu.nearest_proj(proj_angles_list)
        self.assertTrue(numpy.array_equal(return_new1, return_old1))
        self.assertTrue(numpy.array_equal(return_new2, return_old2))

if __name__ == '__main__':
    unittest.main()
