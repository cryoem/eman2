from __future__ import print_function
from __future__ import division

import unittest

import cPickle as pickle
import os
import sys
from mpi import *
import global_def

mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from sphire.libpy_py3 import sphire_pixel_error as fu

from sphire.tests.sparx_lib import sparx_pixel_error as oldfu


class Test_lib_pixel_error_compare(unittest.TestCase):

    def test_pixel_error_2D_true_should_return_equal_objects(self):
        print("Hello testing")
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.pixel_error_2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params1, ali_params2, r) = argum[0]

        return_new = fu.pixel_error_2D(ali_params1, ali_params2, r)
        return_old = oldfu.pixel_error_2D(ali_params1, ali_params2, r)


        self.assertEqual(return_new, return_old)


    def test_max_3D_pixel_error_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.max_3D_pixel_error")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (t1, t2, r) = argum[0]

        return_new = fu.max_3D_pixel_error(t1, t2, r)
        return_old = oldfu.max_3D_pixel_error(t1, t2, r)

        self.assertEqual(return_new, return_old)


    def test_angle_ave_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_ave")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0][0])

        (angle1) = argum[0][0]

        return_new = fu.angle_ave(angle1)
        return_old = oldfu.angle_ave(angle1)

        self.assertEqual(return_new, return_old)


    def test_angle_diff_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_diff")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (angle1,angle2) = argum[0]

        return_new = fu.angle_diff(angle1,angle2)
        return_old = oldfu.angle_diff(angle1,angle2)

        self.assertEqual(return_new, return_old)


    def test_angle_diff_sym_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_diff_sym")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (angle1,angle2, simi) = argum[0]

        return_new = fu.angle_diff_sym(angle1,angle2,simi)
        return_old = oldfu.angle_diff_sym(angle1,angle2,simi)

        self.assertEqual(return_new, return_old)


    def test_align_diff_params_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.align_diff_params")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params1, ali_params2) = argum[0]

        return_new = fu.align_diff_params(ali_params1, ali_params2)
        return_old = oldfu.align_diff_params(ali_params1, ali_params2)

        self.assertEqual(return_new, return_old)



    def test_multi_align_stability_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.multi_align_stability")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params, mir_stab_thld, grp_err_thld, err_thld,print_individual, d) = argum[0]

        return_new = fu.multi_align_stability(ali_params, mir_stab_thld, grp_err_thld, err_thld, print_individual, d)
        return_old = oldfu.multi_align_stability(ali_params, mir_stab_thld, grp_err_thld, err_thld, print_individual, d)

        self.assertEqual(return_new, return_old)






if __name__ == '__main__':
    unittest.main()
