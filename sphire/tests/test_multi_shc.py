from __future__ import print_function
from __future__ import division
import unittest


import numpy
import copy
import math
# import EMAN2_cppwrap as e2cpp
import sys



from ..libpy import sparx_multi_shc as fu

from .sparx_lib import sparx_multi_shc as oldfu

import sparx_fundamentals as sf

from ..libpy import sparx_utilities as ut


import cPickle as pickle
import os
import sys
from mpi import *
import global_def


mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


class Test_lib_multi_shc_compare(unittest.TestCase):

    def test_orient_params_true_should_return_equal_objects(self):
        print(ABSOLUTE_PATH)
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc.orient_params")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, refparams, indexes) = argum[0]
        symmetry_class = argum[1]['symmetry_class']


        return_new = fu.orient_params(params, refparams, indexes,symmetry_class)

        return_old = oldfu.orient_params(params, refparams, indexes,symmetry_class)

        self.assertEqual(return_new, return_old)


    def test_find_common_subset_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc.find_common_subset")
        import sparx_fundamentals
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (projs,target_threshold,minimal_subset_size,symmetry_class) = argum[0]

        return_new = fu.find_common_subset(projs,target_threshold,minimal_subset_size,symmetry_class)

        return_old = oldfu.find_common_subset(projs, target_threshold, minimal_subset_size,symmetry_class)

        self.assertEqual(return_new, return_old)



    def test_ali3d_multishc_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc.ali3d_multishc")
        import sparx_fundamentals
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]

        (dd) = argum[1]

        print('stack values are', stack)
        print('refvol values are', ref_vol)
        print('ali3d_option are ', ali3d_options)
        print('symmetry_class are', symmetry_class)
        print('argument 1 are' , dd)

        mpi_barrier(MPI_COMM_WORLD)

        return_new = fu.ali3d_multishc(stack, ref_vol, ali3d_options, symmetry_class)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali3d_multishc(stack, ref_vol, ali3d_options, symmetry_class)

        mpi_barrier(MPI_COMM_WORLD)

        mpi_finalize()

        # self.assertTrue(return_new, return_old)

if __name__ == '__main__':
    unittest.main()
