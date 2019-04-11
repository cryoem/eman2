from __future__ import print_function
from __future__ import division
import unittest


import numpy
import copy
import math
# import EMAN2_cppwrap as e2cpp
import sys


from sphire.libpy import sparx_multi_shc as fu

from sphire.libpy import sparx_multi_shc as oldfu

# from sphire.libpy import sparx_fundamentals
# sys.modules['sparx_fundamentals'] = sparx_fundamentals

# from ..libpy import sparx_utilities as ut

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
        from sphire.libpy import sparx_fundamentals
        print(ABSOLUTE_PATH)
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.orient_params")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, refparams, indexes) = argum[0]
        symmetry_class = argum[1]['symmetry_class']


        return_new = fu.orient_params(params, refparams,indexes,symmetry_class)

        return_old = oldfu.orient_params(params, refparams,indexes,symmetry_class)

        self.assertEqual(return_new, return_old)


    def test_find_common_subset_true_should_return_equal_objects(self):
        print(ABSOLUTE_PATH)
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.find_common_subset")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (projs,target_threshold,minimal_subset_size,symmetry_class) = argum[0]

        print(symmetry_class)

        print(symmetry_class.sym)

        return_new = fu.find_common_subset(projs,target_threshold,minimal_subset_size,symmetry_class)

        return_old = oldfu.find_common_subset(projs, target_threshold, minimal_subset_size,symmetry_class)

        self.assertEqual(return_new, return_old)


    """
    Cannot Work without proper value of mpi_comm . have to ask markus for this
    """
    def test_ali3d_multishc_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc")
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

        return_new = fu.ali3d_multishc(stack, ref_vol, ali3d_options, symmetry_class, number_of_runs=1)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali3d_multishc(stack, ref_vol, ali3d_options, symmetry_class, number_of_runs=1)

        mpi_barrier(MPI_COMM_WORLD)

        # mpi_finalize()

        if (return_old is not None  and return_new is not None) :
            self.assertTrue(return_new, return_old)



    # def test_ali3d_multishc_2_true_should_return_equal_objects(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc_2")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum)
    #
    #     (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]
    #
    #     (dd) = argum[1]
    #
    #     print('stack values are', stack)
    #     print('refvol values are', ref_vol)
    #     print('ali3d_option are ', ali3d_options)
    #     print('symmetry_class are', symmetry_class)
    #     print('argument 1 are' , dd)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_new = fu.ali3d_multishc_2(stack, ref_vol, ali3d_options, symmetry_class)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = oldfu.ali3d_multishc_2(stack, ref_vol, ali3d_options, symmetry_class)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #

    """
    Cannot Work without proper value of mpi_comm . have to ask markus for this
    """
    # def test_multi_shc_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc.multi_shc")
    #     import sparx_fundamentals
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #
    #     (all_projs, subset, runs_count, ali3d_options) = argum[0]
    #
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc.find_common_subset")
    #     import sparx_fundamentals
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (projs, target_threshold, minimal_subset_size, symmetry_class) = argum[0]
    #
    #     all_projs = projs
    #     n = len(projs[0])
    #     subset = list(range(n))
    #
    #     print(type(subset))
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_new = fu.multi_shc(all_projs, subset, runs_count, ali3d_options, mpi_comm=MPI_COMM_WORLD)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = fu.multi_shc(all_projs, subset, runs_count, ali3d_options, mpi_comm=MPI_COMM_WORLD)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     if (return_old is not None  and return_new is not None) :
    #         self.assertTrue(return_new, return_old)


    def test_mirror_and_reduce_dsym_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.find_common_subset")

        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (projs, target_threshold, minimal_subset_size, symmetry_class) = argum[0]
        n = len(projs[0])
        subset = list(range(n))

        return_new = fu.mirror_and_reduce_dsym(projs, subset, symmetry_class)
        return_old = oldfu.mirror_and_reduce_dsym(projs, subset, symmetry_class)

        if (return_old is not None and return_new is not None):
            self.assertTrue(return_new, return_old)


    """
    Cannot Work without proper value of mpi_comm . have to ask markus for this
    """

    def test_do_volume_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume")

        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,options,iter) = argum[0]

        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.do_volume(data,options,iter, mpi_comm = MPI_COMM_WORLD)

        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.do_volume(data, options, iter, mpi_comm = MPI_COMM_WORLD)

        if (return_old is not None  and return_new is not None) :
            self.assertTrue(return_new, return_old)


if __name__ == '__main__':
    unittest.main()
    mpi_finalize()
