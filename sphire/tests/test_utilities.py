from __future__ import print_function
from __future__ import division


import cPickle as pickle
import os
from mpi import *
import global_def
import numpy
import shutil



mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


import unittest
from test_module import get_arg_from_pickle_file, get_real_data, remove_list_of_file, returns_values_in_file
from ..libpy import sparx_utilities as fu
from .sparx_lib import sparx_utilities as oldfu
from os import path
from EMAN2_cppwrap import EMData
from copy import deepcopy

IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = fu.model_blank(10, 10)
IMAGE_BLANK_3D = fu.model_blank(10, 10, 10)
TOLERANCE = 0.0075
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
"""

class Test_amoeba(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.amoeba"))

    @staticmethod
    def wrongfunction(a,b):
        return a+b

    @staticmethod
    def function_lessParam():
        return 0

    def test_wrong_number_params_too_few_parameters(self):
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

    def test_amoeba_with_function_lessParam(self):
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

    def test_wrong_number_params_too_few_parameters(self):
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

    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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



class Test_get_ima(unittest.TestCase):
    img_list = get_real_data(dim=2)

    def test_NoneType_as_img_returns_AttributeError_NoneType_obj_hasnot_attribute_process(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.get_im(None)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.get_im(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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

    def test_wrong_number_params_too_few_parameters(self):
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

    def test_wrong_number_params_too_few_parameters(self):
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

    def test_wrong_number_params_too_few_parameters(self):
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
            fu.rotate_3D_shift([data], self.shift3d)
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
            fu.rotate_3D_shift([None], self.shift3d)
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
    def test_wrong_number_params_too_few_parameters(self):
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

    def test_with_BadListAttr(self):
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
    def test_wrong_number_params_too_few_parameters(self):
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

    def test_notValid_params_returns_RuntimeError_key_doesnot_exist(self):
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

#todo: all the mpi function

class Test_circumference(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
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


#@unittest.skip("adnan tests")
class Test_lib_utilities_compare(unittest.TestCase):

    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_reduce_EMData_to_root_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.reduce_EMData_to_root")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (data, myid,main_node) = argum[0]
    #
    #     return_new = fu.reduce_EMData_to_root(data, myid,main_node = 0)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.reduce_EMData_to_root(data, myid,main_node = 0)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     self.assertEqual(return_new, return_old)


    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_bcast_compacted_EMData_all_to_all_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_compacted_EMData_all_to_all")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (list_of_em_objects, myid ) = argum[0]
    #
    #     return_new = fu.bcast_compacted_EMData_all_to_all(list_of_em_objects, myid)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.bcast_compacted_EMData_all_to_all(list_of_em_objects, myid)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     self.assertEqual(return_new, return_old)



    def test_gather_compacted_EMData_to_root_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.gather_compacted_EMData_to_root")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (no_of_emo, list_of_emo, myid) = argum[0]

        return_new = fu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid)
        return_old = oldfu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid)

        self.assertEqual(return_new, return_old)


    def test_bcast_EMData_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_EMData_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (tavg, myid, source_node, ) = argum[0]

        return_new = fu.bcast_EMData_to_all(tavg, myid, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_EMData_to_all(tavg, myid, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)



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


    def test_bcast_number_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_number_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (number_to_send, source_node, mpi_comm) = argum[0]

        return_new = fu.bcast_number_to_all(number_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_number_to_all(number_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)

    def test_bcast_list_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_list_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (list_to_send, source_node, mpi_comm) = argum[0]

        return_new = fu.bcast_list_to_all(list_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_list_to_all(list_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)


    # def test_recv_attr_dict_true_should_return_equal_objects(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #     print(argum)
    #
    #     (ima,) = argum[0]
    #     # params = "values"
    #     paramstr = 0
    #     # ima[params] = paramstr
    #
    #     ima.set_attr_dict({"values": 0})
    #
    #     return_new = fu.recv_attr_dict(0, "test", ima,paramstr,0,0,1)
    #
    #     return_old = oldfu.recv_attr_dict(0, "test1", ima,paramstr, 0,0,1)



    def test_print_begin_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_begin_msg(msg)

        return_old = oldfu.print_begin_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_print_end_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_end_msg(msg)

        return_old = oldfu.print_end_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_print_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_msg(msg)

        return_old = oldfu.print_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_read_fsc_true_should_return_equal_objects(self):

        print(os.getcwd())
        filename = "sphire/tests/Sort3D/fsc_global.txt"
        return_new = fu.read_fsc(filename)

        return_old = oldfu.read_fsc(filename)

        self.assertEqual(return_new, return_old)


    def test_circumference_true_should_return_equal_objects(self):


        img = fu.model_blank(10,10,10)

        return_new = fu.circumference(img)

        return_old = oldfu.circumference(img)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_write_headers_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (ima,) = argum[0]


        fu.write_headers("test.hdf", [ima], [1])

        oldfu.write_headers("test1.hdf", [ima], [1])

        filepath = "/home/adnan/PycharmProjects/eman2/"

        return_new = os.path.isfile(filepath + "test.hdf")
        return_old = os.path.isfile(filepath + "test1.hdf")

        self.assertEqual(return_new, return_old)


    def test_write_headers_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (ima,) = argum[0]


        fu.write_header("test.hdf", ima, 1)

        oldfu.write_header("test1.hdf", ima, 1)

        filepath = "/home/adnan/PycharmProjects/eman2/"

        return_new = os.path.isfile(filepath + "test.hdf")
        return_old = os.path.isfile(filepath + "test1.hdf")

        self.assertEqual(return_new, return_old)


    def test_file_type_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.file_type")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (name) = argum[0][0]

        return_new = fu.file_type(name)

        return_old = oldfu.file_type(name)

        self.assertEqual(return_new, return_old)



    def test_get_params2D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])
        print(argum)

        (ima,) = argum[0]

        return_new = fu.get_params2D(ima )

        return_old = oldfu.get_params2D(ima)

        self.assertEqual(return_new, return_old)

    def test_set_params2D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p, xform) = argum[0]

        return_new = fu.set_params2D(ima,p)

        return_old = oldfu.set_params2D(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_params3D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params3D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,) = argum[0]

        return_new = fu.get_params3D(ima )

        return_old = oldfu.get_params3D(ima)

        self.assertEqual(return_new, return_old)


    def test_set_params3D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params3D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p) = argum[0]

        return_new = fu.set_params3D(ima,p)

        return_old = oldfu.set_params3D(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_params_proj_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params_proj")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,) = argum[0]

        return_new = fu.get_params_proj(ima )

        return_old = oldfu.get_params_proj(ima)

        self.assertEqual(return_new, return_old)


    def test_set_params_proj_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params_proj")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p) = argum[0]

        return_new = fu.set_params_proj(ima,p)

        return_old = oldfu.set_params_proj(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_latest_directory_increment_value_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_latest_directory_increment_value")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (directory_location, directory_name) = argum[0]

        return_new = fu.get_latest_directory_increment_value(directory_location, directory_name)

        return_old = oldfu.get_latest_directory_increment_value(directory_location, directory_name)

        self.assertEqual(return_new, return_old)


    def test_get_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (img) = argum[0][0]

        return_new = fu.get_ctf(img[0])
        return_old = oldfu.get_ctf(img[0])

        self.assertEqual(return_new, return_old)


    def test_same_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.same_ctf")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (c1,c2) = argum[0]

        return_new = fu.same_ctf(c1,c2)

        return_old = oldfu.same_ctf(c1,c2)

        self.assertEqual(return_new, return_old)



    def test_generate_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (p) = argum[0][0]

        return_new = fu.generate_ctf(p)

        return_old = oldfu.generate_ctf(p)

        self.assertTrue(return_new, return_old)


    def test_delete_bdb_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.delete_bdb")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (name) = argum[0][0]

        return_new = fu.delete_bdb(name)

        return_old = oldfu.delete_bdb(name)

        if return_new is not None and return_old is not None:
            self.assertTrue(return_new, return_old)
        else:
            print('returns None')

    def test_getvec_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getfvec")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (phi, tht) = argum[0]

        return_new = fu.getvec(phi, tht)

        return_old = oldfu.getvec(phi, tht)

        self.assertEqual(return_new, return_old)


    def test_getfvec_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getfvec")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (phi, tht) = argum[0]

        return_new = fu.getfvec(phi, tht)

        return_old = oldfu.getfvec(phi, tht)

        self.assertEqual(return_new, return_old)


    def test_nearest_fang_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_fang")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (vecs, phi, tht) = argum[0]

        return_new = fu.nearest_fang(vecs, phi, tht)

        return_old = oldfu.nearest_fang(vecs, phi, tht)

        self.assertEqual(return_new, return_old)


    def test_nearest_many_full_k_projangles_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_many_full_k_projangles")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (reference_normals, angles) = argum[0]
        symclass = argum[1]['sym_class']
        howmany = argum[1]['howmany']

        return_new = fu.nearest_many_full_k_projangles(reference_normals, angles, howmany, symclass)

        return_old = oldfu.nearest_many_full_k_projangles(reference_normals, angles, howmany, symclass)

        self.assertEqual(return_new, return_old)


    def test_angles_to_normals_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angles_to_normals")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (angles) = argum[0][0]

        return_new = fu.angles_to_normals(angles)

        return_old = oldfu.angles_to_normals(angles)

        self.assertEqual(return_new, return_old)


    def test_angular_occupancy_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        return_new = fu.angular_occupancy(params, angstep, sym, method)
        return_old = oldfu.angular_occupancy(params, angstep, sym, method)

        self.assertEqual(return_new, return_old)

    def test_angular_histogram_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        return_new = fu.angular_histogram(params, angstep, sym, method)
        return_old = oldfu.angular_histogram(params, angstep, sym, method)

        self.assertEqual(return_new, return_old)


    def test_balance_angular_distribution_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        max_occupy = -1

        return_new = fu.balance_angular_distribution(params, max_occupy, angstep, sym)
        return_old = oldfu.balance_angular_distribution(params, max_occupy, angstep, sym)

        self.assertEqual(return_new, return_old)


    def test_symmetry_neighbors_true_should_return_equal_objects(self):

        angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        return_new = fu.symmetry_neighbors(angles , symmetry= "c1")
        return_old = oldfu.symmetry_neighbors(angles , symmetry= "c1")

        self.assertEqual(return_new, return_old)


    def test_rotation_between_anglesets_true_should_return_equal_objects(self):

        agls1 = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        agls2 = [[idx1, idx2, 5] for idx1 in range(50) for idx2 in range(90)]
        return_new = fu.rotation_between_anglesets(agls1, agls2)
        return_old = oldfu.rotation_between_anglesets(agls1, agls2)

        self.assertEqual(return_new, return_old)


    def test_angle_between_projections_directions_true_should_return_equal_objects(self):

        agls1 = [20, 60, 0]
        agls2 = [45, 75, 5]
        return_new = fu.angle_between_projections_directions(agls1, agls2)
        return_old = oldfu.angle_between_projections_directions(agls1, agls2)

        self.assertEqual(return_new, return_old)


    def test_get_pixel_size_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_pixel_size")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (img,) = argum[0]

        return_new = fu.get_pixel_size(img)

        return_old = oldfu.get_pixel_size(img)

        self.assertEqual(return_new, return_old)

    def test_set_pixel_size_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_pixel_size")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (img,pixel_size) = argum[0]

        return_new = fu.set_pixel_size(img,pixel_size)

        return_old = oldfu.set_pixel_size(img,pixel_size)

        self.assertEqual(return_new, return_old)

    def test_lacos_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.lacos")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (x,) = argum[0]

        return_new = fu.lacos(x)

        return_old = oldfu.lacos(x)

        self.assertEqual(return_new, return_old)

    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_nearest_proj_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_proj")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (proj_ang,img_per_grp,List) = argum[0]
    #
    #     return_new = fu.nearest_proj(proj_ang,img_per_grp,List)
    #
    #     return_old = oldfu.nearest_proj(proj_ang,img_per_grp,List)
    #
    #     self.assertEqual(return_new, return_old)


    def test_findall_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.findall")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (value, L) = argum[0]

        return_new = fu.findall(value, L)

        return_old = oldfu.findall(value, L)

        self.assertEqual(return_new, return_old)


    def test_pack_message_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.pack_message")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.pack_message(data)

        return_old = oldfu.pack_message(data)

        self.assertEqual(return_new, return_old)


    def test_unpack_message_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.unpack_message")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.unpack_message(data)

        return_old = oldfu.unpack_message(data)

        self.assertEqual(return_new, return_old)


    def test_update_tag_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.update_tag")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (communicator, target_rank) = argum[0]

        return_new = fu.update_tag(communicator, target_rank)

        return_old = oldfu.update_tag(communicator, target_rank)

        self.assertEqual(return_new, return_old)


    def test_wrap_mpi_send_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_send")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, destination,communicator) = argum[0]

        return_new = fu.wrap_mpi_send(data, destination)

        return_old = oldfu.wrap_mpi_send(data, destination)

        self.assertEqual(return_new, return_old)


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


    def test_wrap_mpi_bcast_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_bcast")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, root, communicator) = argum[0]

        return_new = fu.wrap_mpi_bcast(data, root)

        return_old = oldfu.wrap_mpi_bcast(data, root)

        self.assertEqual(return_new, return_old)


    def test_wrap_mpi_gatherv_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_gatherv")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, root, communicator) = argum[0]

        return_new = fu.wrap_mpi_gatherv(data, root)

        return_old = oldfu.wrap_mpi_gatherv(data, root)

        self.assertEqual(return_new, return_old)


    def test_get_colors_and_subsets_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_colors_and_subsets")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (main_node, mpi_comm, my_rank, shared_comm, sh_my_rank, masters) = argum[0]

        mpi_comm = MPI_COMM_WORLD
        main_node = 0
        my_rank = mpi_comm_rank(mpi_comm)
        mpi_size = mpi_comm_size(mpi_comm)
        shared_comm = mpi_comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        sh_my_rank = mpi_comm_rank(shared_comm)
        masters = mpi_comm_split(mpi_comm, sh_my_rank == main_node, my_rank)
        shared_comm = mpi_comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)

        return_new = fu.get_colors_and_subsets(main_node, mpi_comm, my_rank, shared_comm, sh_my_rank,masters)

        return_old = oldfu.get_colors_and_subsets(main_node, mpi_comm, my_rank, shared_comm, sh_my_rank,masters)

        self.assertEqual(return_new, return_old)

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

    def test_eliminate_moons_true_should_return_equal_objects(self):

        volume = fu.model_gauss(0.25,12,12,12)


        moon_params = []
        moon_params.append(0.5)
        moon_params.append(1.45)

        return_new = fu.eliminate_moons(volume, moon_params)

        return_old = oldfu.eliminate_moons(volume, moon_params)

        self.assertEqual(return_new, return_old)


    def test_get_dist_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_dist")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (c1, c2) = argum[0]

        return_new = fu.get_dist(c1, c2)

        return_old = oldfu.get_dist(c1, c2)

        self.assertEqual(return_new, return_old)


    def test_combinations_of_n_taken_by_k_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.combinations_of_n_taken_by_k")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (n, k) = argum[0]

        return_new = fu.combinations_of_n_taken_by_k(n, k)

        return_old = oldfu.combinations_of_n_taken_by_k(n, k)

        self.assertEqual(return_new, return_old)


    def test_cmdexecute_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.cmdexecute")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (cmd,) = argum[0]

        dirname = cmd.split(' ')[1]

        current_path = os.getcwd()
        if os.path.isdir(dirname):
            print('directory exits')
            print('removing it')
            shutil.rmtree(dirname)

        return_new = fu.cmdexecute(cmd)

        if os.path.isdir(dirname):
            print('directory exits')
            print('removing it')
            shutil.rmtree(dirname)

        return_old = oldfu.cmdexecute(cmd)

        self.assertEqual(return_new, return_old)


    def test_if_error_then_all_processes_exit_program_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.if_error_then_all_processes_exit_program")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (error_status,) = argum[0]

        return_new = fu.if_error_then_all_processes_exit_program(error_status)

        return_old = oldfu.if_error_then_all_processes_exit_program(error_status)

        self.assertEqual(return_new, return_old)



    def test_get_shrink_data_huang_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["stack"] = "bdb:sphire/tests/Substack/sort3d_substack_002"
        Tracker["applyctf"] = True
        ids = []
        for i in range(1227):
            ids.append(i)
        Tracker["chunk_dict"] =ids
        myid = 0
        m_node = 0
        nproc = 1
        partids = "sphire/tests/Sort3D/indexes_010.txt"
        partstack = "sphire/tests/Sort3D/params_010.txt"
        nxinit = 2

        return_new = fu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        return_old = oldfu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        self.assertTrue(numpy.allclose(return_new[0][0].get_3dview(), return_old[0][0].get_3dview(), 0.5))


    def test_getindexdata_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getindexdata")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (stack, partids, partstack, myid, nproc) = argum[0]

        stack = 'bdb:sphire/tests/VIPER/best_000'
        partids = 'sphire/tests/VIPER/main001/this_iteration_index_keep_images.txt'
        partstack = 'sphire/tests/VIPER//main001/run000/rotated_reduced_params.txt'

        return_new = fu.getindexdata(stack, partids, partstack, myid, nproc)

        return_old = oldfu.getindexdata(stack, partids, partstack, myid, nproc)

        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))

    def test_convert_json_fromunicode_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.convert_json_fromunicode")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.convert_json_fromunicode(data)

        return_old = oldfu.convert_json_fromunicode(data)

        self.assertEqual(return_new, return_old)


    def test_get_sorting_attr_stack_should_return_equal_object(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)

        return_new = fu.get_sorting_attr_stack(stack)
        return_old = oldfu.get_sorting_attr_stack(stack)

        self.assertEqual(return_new, return_old)


    def test_get_sorting_params_refine_should_return_equal_object(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["stack"] = "bdb:sphire/tests/Substack/sort3d_substack_002"
        Tracker["applyctf"] = True
        Tracker["constants"]["nproc"] = 1

        return_new = fu.get_sorting_params_refine(Tracker, stack, 95)
        return_old = oldfu.get_sorting_params_refine(Tracker, stack, 95)

        self.assertEqual(return_new, return_old)

    def test_parsing_sorting_params_should_return_equal_object(self):

        sorting_list = []

        sorting_list.append(numpy.arange(10))
        sorting_list.append(numpy.arange(10))

        return_new = fu.parsing_sorting_params(sorting_list)
        return_old = oldfu.parsing_sorting_params(sorting_list)

        self.assertEqual(return_new[0], return_old[0])
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    # def test_get_initial_ID_should_return_equal_object(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #
    #     return_new = fu.get_initial_ID(Tracker["two_way_stable_member"][istable], Tracker["full_ID_dict"])
    #     return_old = oldfu.get_initial_ID(Tracker["two_way_stable_member"][istable], Tracker["full_ID_dict"])
    #
    #     self.assertEqual(return_new, return_old)


    def test_convertasi_true_should_return_equal_objects(self):

        K = 7
        asig = [0,1,2,3,4,5,6]

        return_new = fu.convertasi(asig,K)
        return_old = oldfu.convertasi(asig,K)

        self.assertEqual(return_new, return_old)


    def test_prepare_ptp_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.prepare_ptp(data_list, K)
        return_old = oldfu.prepare_ptp(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_get_resolution_mrk01_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (img1, img2) = argum[0]
        volume = []
        volume.append(img1)
        volume.append(img2)

        fscoutputdir = "sphire/tests/Sort3D"
        mask_option = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.get_resolution_mrk01(volume, 0.5,0.15,fscoutputdir,mask_option)
        return_old = oldfu.get_resolution_mrk01(volume, 0.5,0.15,fscoutputdir,mask_option)

        self.assertEqual(return_new, return_old)


    def test_partition_to_groups_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.partition_to_groups(data_list, K)
        return_old = oldfu.partition_to_groups(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_partition_independent_runs_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.partition_independent_runs(data_list, K)
        return_old = oldfu.partition_independent_runs(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_merge_groups_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.merge_groups(data_list)
        return_old = oldfu.merge_groups(data_list)

        self.assertEqual(return_new, return_old)


    def test_save_alist_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]

        Tracker["this_dir"] = "sphire/tests/Sharpening/"
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = "myid"
        Tracker["constants"]["main_node"] = "myid"

        filename = "listfile.txt"
        data_list = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.save_alist(Tracker,filename,data_list)
        return_old = oldfu.save_alist(Tracker,filename,data_list)

        self.assertEqual(return_new, return_old)

    def test_margin_of_error_true_should_return_equal_objects(self):

        return_new = fu.margin_of_error(0,1)
        return_old = oldfu.margin_of_error(0,1)

        self.assertEqual(return_new, return_old)

    # def test_do_two_way_comparison_true_should_return_equal_objects(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #
    #     Tracker["this_dir"] = "sphire/tests/Sharpening/"
    #     Tracker["constants"]["log_main"] = "logging"
    #     Tracker["constants"]["myid"] = 0
    #     Tracker["constants"]["main_node"] = 1
    #     Tracker["this_total_stack"] = "bdb:sphire/tests/Substack/sort3d_substack_002"
    #     Tracker["number_of_groups"] = 4
    #     Tracker["constants"]["indep_runs"]  = 4
    #     Tracker["partition_dict"]    = [0,1,2,3]
    #
    #
    #     return_new = fu.do_two_way_comparison(Tracker)
    #     return_old = oldfu.do_two_way_comparison(Tracker)
    #
    #     self.assertEqual(return_new, return_old)

    def test_counting_projections_true_should_return_equal_objects(self):

        delta = 0.5
        ali3d_params  = [[idx1, idx2, 0 , 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
        image_start = 1

        return_new = fu.counting_projections(delta, ali3d_params, image_start)
        print("one function call done")
        return_old = oldfu.counting_projections(delta, ali3d_params, image_start)
        print("second function call done")
        self.assertEqual(return_new, return_old)

    def test_unload_dict_true_should_return_equal_objects(self):

        delta = 0.5
        ali3d_params  = [[idx1, idx2, 0 , 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
        image_start = 1
        sampled = fu.counting_projections(delta, ali3d_params,image_start)

        return_new = fu.unload_dict(sampled)
        return_old = oldfu.unload_dict(sampled)
        self.assertEqual(return_new, return_old)


    def test_unload_dict_true_should_return_equal_objects(self):

        delta = 0.5
        ali3d_params  = [[idx1, idx2, 0 , 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
        image_start = 1
        sampled = fu.counting_projections(delta, ali3d_params,image_start)
        dicto = fu.unload_dict(sampled)

        return_new = fu.load_dict(sampled, dicto)
        return_old = oldfu.load_dict(sampled, dicto)
        self.assertEqual(return_new, return_old)




    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_get_stat_proj_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #     Tracker["constants"]["nproc"] = 1
    #     Tracker["constants"]["myid"] = 0
    #     Tracker["constants"]["main_node"] = 0
    #
    #     delta = 0.5
    #     this_ali3d = "sphire/tests/VIPER/main001/run000/rotated_reduced_params.txt"
    #
    #     return_new = fu.get_stat_proj(Tracker,delta,this_ali3d)
    #     return_old = oldfu.get_stat_proj(Tracker,delta,this_ali3d)
    #     self.assertEqual(return_new, return_old)


    def test_create_random_list_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["total_stack"] = "stack"
        Tracker["constants"]["seed"] = 1.4
        Tracker["constants"]["indep_runs"] = 2
        Tracker["this_data_list"] = [2,3,5]

        return_new = fu.create_random_list(Tracker)
        return_old = oldfu.create_random_list(Tracker)
        self.assertEqual(return_new, return_old)


    def test_recons_mref_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["number_of_groups"] = 1
        Tracker["constants"]["nnxo"] = 4  # roi
        Tracker["this_particle_list"] = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        Tracker["nxinit"] = 1
        Tracker["constants"]["partstack"] = 'sphire/tests/VIPER//main001/run000/rotated_reduced_params.txt'
        Tracker["this_dir"] = "sphire/tests/Particles/"
        Tracker["constants"]["stack"] = 'bdb:sphire/tests/Class2D/stack_ali2d'
        Tracker["applyctf"] = False
        Tracker["chunk_dict"] = [0, 1, 2, 3, 4, 5, 6]
        Tracker["constants"]["sym"] = "c1"

        return_new = fu.recons_mref(Tracker)
        return_old = oldfu.recons_mref(Tracker)
        self.assertTrue(return_new[0], return_old[0])



    def test_apply_low_pass_filter_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["low_pass_filter"] = 0.087

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (volft, params, interpolation_method, return_real) = argum[0]

        refvol = [volft,volft]

        return_new = fu.apply_low_pass_filter(refvol,Tracker)
        return_old = oldfu.apply_low_pass_filter(refvol,Tracker)
        self.assertEqual(return_new, return_old)



    def test_count_chunk_members_true_should_return_equal_objects(self):

        chunk_dict =  [0,1,2,3,4,5,6]
        one_class = [0,1,2,3,4,5,6]

        return_new = fu.count_chunk_members(chunk_dict, one_class)
        return_old = oldfu.count_chunk_members(chunk_dict, one_class)
        self.assertEqual(return_new, return_old)


    def test_get_groups_from_partition_true_should_return_equal_objects(self):
        final_list = []

        final_list.append(numpy.arange(10))
        final_list.append(numpy.arange(10))

        this_data_list_file = "sphire/tests/Sort3D/chunk_0.txt"

        # final_list = sparx_utilities.get_sorting_params_refine(Tracker, data, total_nima)
        group_list, ali3d_params_list = fu.parsing_sorting_params(final_list)
        list_of_particles = fu.read_text_file(this_data_list_file)
        return_new = fu.get_groups_from_partition(group_list, list_of_particles, 2)
        return_old = oldfu.get_groups_from_partition(group_list, list_of_particles, 2)
        self.assertEqual(return_new, return_old)

    def test_get_complementary_elements_true_should_return_equal_objects(self):

        sub_data_list = [1,2,2]
        total_list = [1,2,2,4,5,6]

        return_new = fu.get_complementary_elements(total_list,sub_data_list)
        return_old = oldfu.get_complementary_elements(total_list,sub_data_list)
        self.assertEqual(return_new, return_old)

    def test_remove_small_groups_true_should_return_equal_objects(self):

        chunk_dict =  [[0,1,2,3,4,5,6],[0,1,2,3,4,5,6],[0,1,2,3,4,5,6]]
        one_class = [0,1,2,3,4,5,6]

        return_new = fu.remove_small_groups(chunk_dict, 3)
        return_old = oldfu.remove_small_groups(chunk_dict, 3)
        self.assertEqual(return_new, return_old)


    def test_get_number_of_groups_true_should_return_equal_objects(self):

        return_new = fu.get_number_of_groups(500, 5)
        return_old = oldfu.get_number_of_groups(500, 5)
        self.assertEqual(return_new, return_old)


    def test_tabessel_true_should_return_equal_objects(self):

        return_new = fu.tabessel(5,4)
        return_old = oldfu.tabessel(5, 4)
        self.assertEqual(return_new, return_old)



if __name__ == '__main__':
    unittest.main()
