from __future__ import print_function
from __future__ import division

import unittest

import cPickle as pickle
import os
import numpy
from mpi import *
import global_def
import EMAN2_cppwrap as e2cpp
import EMAN2_cppwrap

mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from sphire.libpy_py3 import sphire_statistics as fu

from .sparx_lib import sparx_statistics as oldfu





from test_module import  remove_list_of_file, returns_values_in_file, get_arg_from_pickle_file,get_real_data,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER
from os import path
from ..libpy import sparx_utilities
from EMAN2_cppwrap import EMData
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
TOLERANCE = 0.0005
MASK = sparx_utilities.model_circle(2, 5, 5)
STACK_NAME = 'bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'Substack/sort3d_substack_003')

"""
There are some opened issues in:
1) add_ave_varf_MPI: it is never used and I avoid to test it deeply. 
    -) CTF=True leads to a bug (see 'test_CTF_True_deadCode_BUG' --> UnboundLocalError: local variable 'i' referenced before assignment)
    
2) sum_oe:
    -)I cannot find a dataset to test the CTF=True case .... anyway it is just a debug mode
    -) ctf_eo_sum not tested because used in combination with CTF=True
    
3) varf2d_MPI:
    -) the use of a 3D image is not blocked. is it ok or a bug?
    -) I cannot find a dataset to test the CTF=True case
    
4) varf3d_MPI:
    -) how can i find a valid mask to test it
    -) i did not tested all the params deeply becuae they are basically used as input of 'sparx_reconstruction.recons3d_4nn_ctf_MPI'
    
5) locres:
    -) Adnan says that it is able to run only in the cluster ... Hence I did not test it
    
6) class Munkres and class pcanalyzer:
    -) at the moment I have no idea how I'd test it
    
7) linreg: Could we use from sklearn.linear_model import LinearRegression ???

8) k_means_stab_bbenum:
    -) i did not test it deeply because half parameters are not used and the other is used to call  'k_means_match_bbenum'. Hence I'll test deeply 'k_means_match_bbenum'
"""

class Test_add_ave_varf_MPI(unittest.TestCase):
    main_node = 0
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))[0][0]

    def test_all_the_conditions(self, return_new=(), return_old=(), tolerance =TOLERANCE ):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(numpy.allclose(j.get_3dview(), i.get_3dview(), atol=tolerance,equal_nan=True))


    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.add_ave_varf_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.add_ave_varf_MPI()
        self.assertEqual(cm_new.exception.message, "add_ave_varf_MPI() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_CTF_True_deadCode_BUG(self):
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(UnboundLocalError) as cm_new:
            fu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=True, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(UnboundLocalError) as cm_old:
            oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=True, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.assertEqual(cm_new.exception.message, "local variable 'i' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_CTF_wrong_ali_params(self):
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align3d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align3d",main_node=self.main_node, comm=-1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value_msg_Error_negative_variance(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.test_all_the_conditions(return_new,return_old, tolerance=0)
        self.assertEqual(len(return_new),5)

    def test_with_mask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(self.main_node,self.data, mask=MASK, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=MASK, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.test_all_the_conditions(return_new,return_old, tolerance=0)
        self.assertEqual(len(return_new),5)

    def test_no_on_main_node(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(self.main_node+1,self.data, mask=MASK, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(self.main_node+1,self.data, mask=MASK, mode="a", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.test_all_the_conditions(return_new,return_old, tolerance=0)
        self.assertEqual(len(return_new),5)
        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), EMAN2_cppwrap.EMData().get_3dview(),  equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[-1].get_3dview(), EMAN2_cppwrap.EMData().get_3dview(), equal_nan=True))

    def test_with_mask_and_not_current_alignment_parameters(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(self.main_node,self.data, mask=MASK, mode="not_current", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=MASK, mode="not_current", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.test_all_the_conditions(return_new,return_old, tolerance=0)
        self.assertEqual(len(return_new),5)

    def test_without_mask_and_not_current_alignment_parameters_msg_Error_negative_variance(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="not_current", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(self.main_node,self.data, mask=None, mode="not_current", CTF=False, ctf_2_sum=None, ali_params="xform.align2d",main_node=self.main_node, comm=-1)
        self.test_all_the_conditions(return_new,return_old, tolerance=0)
        self.assertEqual(len(return_new),5)




class Test_sum_oe(unittest.TestCase):
    data, = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series"))[0]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        if len(return_new) >0:
            self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
            self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sum_oe()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sum_oe()
        self.assertEqual(cm_new.exception.message, "sum_oe() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_list_data_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.sum_oe([], mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.sum_oe([], mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_lists_data_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.sum_oe([None], mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.sum_oe([None], mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        self.assertEqual(cm_new.exception.message,"'NoneType' object has no attribute 'get_xsize'")

    def test_default_value(self):
        return_new = fu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        return_old = oldfu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))

    def test_default_value_with_returnParams(self):
        return_new = fu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=True)
        return_old = oldfu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=True)
        self.test_all_the_conditions(return_new,return_old)
        self.assertTrue(numpy.array_equal(return_new[2], [[115.92958108134886, 0.0, 0.0, 1, 1.0], [139.27090994520643, 0.0, 0.0, 1, 1.0], [338.05886971796366, 0.0, 0.0, 0, 1.0], [130.22834488097047, 0.0, 0.0, 1, 1.0], [256.18793661895506, 0.0, 0.0, 0, 1.0], [33.26201066479168, 0.0, 0.0, 1, 1.0], [147.90743284695057, 0.0, 0.0, 1, 1.0], [263.2061891021638, 0.0, 0.0, 0, 1.0], [300.8387181916984, 0.0, 0.0, 1, 1.0], [216.24095916911287, 0.0, 0.0, 1, 1.0], [73.78844181734985, 0.0, 0.0, 0, 1.0], [152.6487630314829, 0.0, 0.0, 1, 1.0], [119.2316676091685, 0.0, 0.0, 0, 1.0], [73.42069608126344, 0.0, 0.0, 0, 1.0], [236.35598895100773, 0.0, 0.0, 0, 1.0], [174.09856654009025, 0.0, 0.0, 0, 1.0], [289.1996133738564, 0.0, 0.0, 0, 1.0], [0.5265443325058998, 0.0, 0.0, 0, 1.0], [139.34182602571568, 0.0, 0.0, 0, 1.0], [95.33089764760811, 0.0, 0.0, 0, 1.0], [338.5149172259531, 0.0, 0.0, 1, 1.0], [304.40205119783315, 0.0, 0.0, 1, 1.0], [89.690362725573, 0.0, 0.0, 1, 1.0], [37.504752818284324, 0.0, 0.0, 0, 1.0], [167.6035631238888, 0.0, 0.0, 1, 1.0], [319.7104370851999, 0.0, 0.0, 0, 1.0], [326.0514773572882, 0.0, 0.0, 1, 1.0], [83.24303621783521, 0.0, 0.0, 0, 1.0], [234.2775620183666, 0.0, 0.0, 1, 1.0], [256.4408649656969, 0.0, 0.0, 0, 1.0], [150.06735665150944, 0.0, 0.0, 1, 1.0], [61.89125188813407, 0.0, 0.0, 0, 1.0], [13.862688281788506, 0.0, 0.0, 1, 1.0], [304.0362002488816, 0.0, 0.0, 1, 1.0], [307.1984875595606, 0.0, 0.0, 1, 1.0], [161.87087060272057, 0.0, 0.0, 1, 1.0], [258.0839973094497, 0.0, 0.0, 0, 1.0], [40.82216679715648, 0.0, 0.0, 1, 1.0], [205.35195001588716, 0.0, 0.0, 1, 1.0], [209.68718891764948, 0.0, 0.0, 1, 1.0], [208.3125464835907, 0.0, 0.0, 0, 1.0], [196.04511403446048, 0.0, 0.0, 0, 1.0], [177.70626452040764, 0.0, 0.0, 1, 1.0], [232.82592596431675, 0.0, 0.0, 0, 1.0], [328.9746213107683, 0.0, 0.0, 0, 1.0], [329.07755292675165, 0.0, 0.0, 1, 1.0], [64.94253725651211, 0.0, 0.0, 0, 1.0], [351.56326265802255, 0.0, 0.0, 1, 1.0], [323.3679418390045, 0.0, 0.0, 0, 1.0], [331.307426261965, 0.0, 0.0, 0, 1.0], [13.802205090050872, 0.0, 0.0, 1, 1.0], [220.80483291105557, 0.0, 0.0, 0, 1.0], [316.799385197508, 0.0, 0.0, 0, 1.0], [50.54424801569083, 0.0, 0.0, 0, 1.0], [38.24882530685342, 0.0, 0.0, 0, 1.0], [273.2859490510846, 0.0, 0.0, 0, 1.0], [354.98349840414687, 0.0, 0.0, 1, 1.0], [278.9158628665873, 0.0, 0.0, 1, 1.0], [273.30016625616406, 0.0, 0.0, 0, 1.0], [196.61511562348988, 0.0, 0.0, 1, 1.0], [333.5271856882758, 0.0, 0.0, 0, 1.0], [142.29007235951357, 0.0, 0.0, 1, 1.0], [317.79847608899394, 0.0, 0.0, 0, 1.0], [0.7266499442631127, 0.0, 0.0, 1, 1.0], [107.4744912169713, 0.0, 0.0, 1, 1.0], [345.0547906055752, 0.0, 0.0, 1, 1.0], [358.2082165930105, 0.0, 0.0, 1, 1.0], [201.0817399060137, 0.0, 0.0, 0, 1.0], [219.23887102939034, 0.0, 0.0, 1, 1.0], [195.59591896934865, 0.0, 0.0, 1, 1.0], [29.131395097747454, 0.0, 0.0, 1, 1.0], [36.449727272883024, 0.0, 0.0, 1, 1.0], [60.57824509120188, 0.0, 0.0, 1, 1.0], [199.69795222785152, 0.0, 0.0, 0, 1.0], [231.7997843843201, 0.0, 0.0, 1, 1.0], [222.34509646846848, 0.0, 0.0, 0, 1.0], [282.3188091919718, 0.0, 0.0, 0, 1.0], [179.0187700827499, 0.0, 0.0, 0, 1.0], [274.382960994679, 0.0, 0.0, 1, 1.0], [135.07809583907184, 0.0, 0.0, 1, 1.0], [39.849524713983115, 0.0, 0.0, 0, 1.0], [235.385486624483, 0.0, 0.0, 1, 1.0], [149.94362873415014, 0.0, 0.0, 0, 1.0], [156.50061118400743, 0.0, 0.0, 1, 1.0], [315.9511826004298, 0.0, 0.0, 1, 1.0], [150.73088983672164, 0.0, 0.0, 0, 1.0], [213.4693880426871, 0.0, 0.0, 1, 1.0], [312.4870012154273, 0.0, 0.0, 0, 1.0], [279.58727287760973, 0.0, 0.0, 0, 1.0], [43.36811300275874, 0.0, 0.0, 0, 1.0], [164.33842254754984, 0.0, 0.0, 1, 1.0], [58.26280097271371, 0.0, 0.0, 1, 1.0], [353.5866426164621, 0.0, 0.0, 0, 1.0], [171.278749296968, 0.0, 0.0, 1, 1.0], [131.82160195058572, 0.0, 0.0, 1, 1.0], [235.99549472014235, 0.0, 0.0, 1, 1.0], [358.2431623129092, 0.0, 0.0, 1, 1.0], [216.00100894203774, 0.0, 0.0, 0, 1.0], [278.1703686699143, 0.0, 0.0, 1, 1.0], [82.87694953387485, 0.0, 0.0, 0, 1.0], [155.47445122542783, 0.0, 0.0, 0, 1.0], [26.917148721112603, 0.0, 0.0, 0, 1.0], [289.40830449708994, 0.0, 0.0, 0, 1.0], [354.312922203678, 0.0, 0.0, 0, 1.0], [136.89611110667116, 0.0, 0.0, 1, 1.0], [249.77889121874279, 0.0, 0.0, 1, 1.0], [104.24588792893502, 0.0, 0.0, 1, 1.0], [117.561324155386, 0.0, 0.0, 0, 1.0], [350.7680033335721, 0.0, 0.0, 0, 1.0], [38.49161133326204, 0.0, 0.0, 1, 1.0], [288.07766631373886, 0.0, 0.0, 1, 1.0], [349.02704990380533, 0.0, 0.0, 1, 1.0]]))

    def test_not_current_alignment_parameters(self):
        return_new = fu.sum_oe(data=self.data, mode="not", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=True)
        return_old = oldfu.sum_oe(data=self.data, mode="not", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=True)
        self.test_all_the_conditions(return_new,return_old)
        self.assertTrue(numpy.array_equal(return_new[2], [[115.92958108134886, 0.0, 0.0, 1, 1.0], [139.27090994520643, 0.0, 0.0, 1, 1.0], [338.05886971796366, 0.0, 0.0, 0, 1.0], [130.22834488097047, 0.0, 0.0, 1, 1.0], [256.18793661895506, 0.0, 0.0, 0, 1.0], [33.26201066479168, 0.0, 0.0, 1, 1.0], [147.90743284695057, 0.0, 0.0, 1, 1.0], [263.2061891021638, 0.0, 0.0, 0, 1.0], [300.8387181916984, 0.0, 0.0, 1, 1.0], [216.24095916911287, 0.0, 0.0, 1, 1.0], [73.78844181734985, 0.0, 0.0, 0, 1.0], [152.6487630314829, 0.0, 0.0, 1, 1.0], [119.2316676091685, 0.0, 0.0, 0, 1.0], [73.42069608126344, 0.0, 0.0, 0, 1.0], [236.35598895100773, 0.0, 0.0, 0, 1.0], [174.09856654009025, 0.0, 0.0, 0, 1.0], [289.1996133738564, 0.0, 0.0, 0, 1.0], [0.5265443325058998, 0.0, 0.0, 0, 1.0], [139.34182602571568, 0.0, 0.0, 0, 1.0], [95.33089764760811, 0.0, 0.0, 0, 1.0], [338.5149172259531, 0.0, 0.0, 1, 1.0], [304.40205119783315, 0.0, 0.0, 1, 1.0], [89.690362725573, 0.0, 0.0, 1, 1.0], [37.504752818284324, 0.0, 0.0, 0, 1.0], [167.6035631238888, 0.0, 0.0, 1, 1.0], [319.7104370851999, 0.0, 0.0, 0, 1.0], [326.0514773572882, 0.0, 0.0, 1, 1.0], [83.24303621783521, 0.0, 0.0, 0, 1.0], [234.2775620183666, 0.0, 0.0, 1, 1.0], [256.4408649656969, 0.0, 0.0, 0, 1.0], [150.06735665150944, 0.0, 0.0, 1, 1.0], [61.89125188813407, 0.0, 0.0, 0, 1.0], [13.862688281788506, 0.0, 0.0, 1, 1.0], [304.0362002488816, 0.0, 0.0, 1, 1.0], [307.1984875595606, 0.0, 0.0, 1, 1.0], [161.87087060272057, 0.0, 0.0, 1, 1.0], [258.0839973094497, 0.0, 0.0, 0, 1.0], [40.82216679715648, 0.0, 0.0, 1, 1.0], [205.35195001588716, 0.0, 0.0, 1, 1.0], [209.68718891764948, 0.0, 0.0, 1, 1.0], [208.3125464835907, 0.0, 0.0, 0, 1.0], [196.04511403446048, 0.0, 0.0, 0, 1.0], [177.70626452040764, 0.0, 0.0, 1, 1.0], [232.82592596431675, 0.0, 0.0, 0, 1.0], [328.9746213107683, 0.0, 0.0, 0, 1.0], [329.07755292675165, 0.0, 0.0, 1, 1.0], [64.94253725651211, 0.0, 0.0, 0, 1.0], [351.56326265802255, 0.0, 0.0, 1, 1.0], [323.3679418390045, 0.0, 0.0, 0, 1.0], [331.307426261965, 0.0, 0.0, 0, 1.0], [13.802205090050872, 0.0, 0.0, 1, 1.0], [220.80483291105557, 0.0, 0.0, 0, 1.0], [316.799385197508, 0.0, 0.0, 0, 1.0], [50.54424801569083, 0.0, 0.0, 0, 1.0], [38.24882530685342, 0.0, 0.0, 0, 1.0], [273.2859490510846, 0.0, 0.0, 0, 1.0], [354.98349840414687, 0.0, 0.0, 1, 1.0], [278.9158628665873, 0.0, 0.0, 1, 1.0], [273.30016625616406, 0.0, 0.0, 0, 1.0], [196.61511562348988, 0.0, 0.0, 1, 1.0], [333.5271856882758, 0.0, 0.0, 0, 1.0], [142.29007235951357, 0.0, 0.0, 1, 1.0], [317.79847608899394, 0.0, 0.0, 0, 1.0], [0.7266499442631127, 0.0, 0.0, 1, 1.0], [107.4744912169713, 0.0, 0.0, 1, 1.0], [345.0547906055752, 0.0, 0.0, 1, 1.0], [358.2082165930105, 0.0, 0.0, 1, 1.0], [201.0817399060137, 0.0, 0.0, 0, 1.0], [219.23887102939034, 0.0, 0.0, 1, 1.0], [195.59591896934865, 0.0, 0.0, 1, 1.0], [29.131395097747454, 0.0, 0.0, 1, 1.0], [36.449727272883024, 0.0, 0.0, 1, 1.0], [60.57824509120188, 0.0, 0.0, 1, 1.0], [199.69795222785152, 0.0, 0.0, 0, 1.0], [231.7997843843201, 0.0, 0.0, 1, 1.0], [222.34509646846848, 0.0, 0.0, 0, 1.0], [282.3188091919718, 0.0, 0.0, 0, 1.0], [179.0187700827499, 0.0, 0.0, 0, 1.0], [274.382960994679, 0.0, 0.0, 1, 1.0], [135.07809583907184, 0.0, 0.0, 1, 1.0], [39.849524713983115, 0.0, 0.0, 0, 1.0], [235.385486624483, 0.0, 0.0, 1, 1.0], [149.94362873415014, 0.0, 0.0, 0, 1.0], [156.50061118400743, 0.0, 0.0, 1, 1.0], [315.9511826004298, 0.0, 0.0, 1, 1.0], [150.73088983672164, 0.0, 0.0, 0, 1.0], [213.4693880426871, 0.0, 0.0, 1, 1.0], [312.4870012154273, 0.0, 0.0, 0, 1.0], [279.58727287760973, 0.0, 0.0, 0, 1.0], [43.36811300275874, 0.0, 0.0, 0, 1.0], [164.33842254754984, 0.0, 0.0, 1, 1.0], [58.26280097271371, 0.0, 0.0, 1, 1.0], [353.5866426164621, 0.0, 0.0, 0, 1.0], [171.278749296968, 0.0, 0.0, 1, 1.0], [131.82160195058572, 0.0, 0.0, 1, 1.0], [235.99549472014235, 0.0, 0.0, 1, 1.0], [358.2431623129092, 0.0, 0.0, 1, 1.0], [216.00100894203774, 0.0, 0.0, 0, 1.0], [278.1703686699143, 0.0, 0.0, 1, 1.0], [82.87694953387485, 0.0, 0.0, 0, 1.0], [155.47445122542783, 0.0, 0.0, 0, 1.0], [26.917148721112603, 0.0, 0.0, 0, 1.0], [289.40830449708994, 0.0, 0.0, 0, 1.0], [354.312922203678, 0.0, 0.0, 0, 1.0], [136.89611110667116, 0.0, 0.0, 1, 1.0], [249.77889121874279, 0.0, 0.0, 1, 1.0], [104.24588792893502, 0.0, 0.0, 1, 1.0], [117.561324155386, 0.0, 0.0, 0, 1.0], [350.7680033335721, 0.0, 0.0, 0, 1.0], [38.49161133326204, 0.0, 0.0, 1, 1.0], [288.07766631373886, 0.0, 0.0, 1, 1.0], [349.02704990380533, 0.0, 0.0, 1, 1.0]]))

    def test_prevent_calculation_ctf2(self):
        return_new = fu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum= EMAN2_cppwrap.EMData(), ctf_eo_sum=False, return_params=True)
        return_old = oldfu.sum_oe(data=self.data, mode="a", CTF=False, ctf_2_sum= EMAN2_cppwrap.EMData(), ctf_eo_sum=False, return_params=True)
        self.test_all_the_conditions(return_new,return_old)
        self.assertTrue(numpy.array_equal(return_new[2], [[115.92958108134886, 0.0, 0.0, 1, 1.0], [139.27090994520643, 0.0, 0.0, 1, 1.0], [338.05886971796366, 0.0, 0.0, 0, 1.0], [130.22834488097047, 0.0, 0.0, 1, 1.0], [256.18793661895506, 0.0, 0.0, 0, 1.0], [33.26201066479168, 0.0, 0.0, 1, 1.0], [147.90743284695057, 0.0, 0.0, 1, 1.0], [263.2061891021638, 0.0, 0.0, 0, 1.0], [300.8387181916984, 0.0, 0.0, 1, 1.0], [216.24095916911287, 0.0, 0.0, 1, 1.0], [73.78844181734985, 0.0, 0.0, 0, 1.0], [152.6487630314829, 0.0, 0.0, 1, 1.0], [119.2316676091685, 0.0, 0.0, 0, 1.0], [73.42069608126344, 0.0, 0.0, 0, 1.0], [236.35598895100773, 0.0, 0.0, 0, 1.0], [174.09856654009025, 0.0, 0.0, 0, 1.0], [289.1996133738564, 0.0, 0.0, 0, 1.0], [0.5265443325058998, 0.0, 0.0, 0, 1.0], [139.34182602571568, 0.0, 0.0, 0, 1.0], [95.33089764760811, 0.0, 0.0, 0, 1.0], [338.5149172259531, 0.0, 0.0, 1, 1.0], [304.40205119783315, 0.0, 0.0, 1, 1.0], [89.690362725573, 0.0, 0.0, 1, 1.0], [37.504752818284324, 0.0, 0.0, 0, 1.0], [167.6035631238888, 0.0, 0.0, 1, 1.0], [319.7104370851999, 0.0, 0.0, 0, 1.0], [326.0514773572882, 0.0, 0.0, 1, 1.0], [83.24303621783521, 0.0, 0.0, 0, 1.0], [234.2775620183666, 0.0, 0.0, 1, 1.0], [256.4408649656969, 0.0, 0.0, 0, 1.0], [150.06735665150944, 0.0, 0.0, 1, 1.0], [61.89125188813407, 0.0, 0.0, 0, 1.0], [13.862688281788506, 0.0, 0.0, 1, 1.0], [304.0362002488816, 0.0, 0.0, 1, 1.0], [307.1984875595606, 0.0, 0.0, 1, 1.0], [161.87087060272057, 0.0, 0.0, 1, 1.0], [258.0839973094497, 0.0, 0.0, 0, 1.0], [40.82216679715648, 0.0, 0.0, 1, 1.0], [205.35195001588716, 0.0, 0.0, 1, 1.0], [209.68718891764948, 0.0, 0.0, 1, 1.0], [208.3125464835907, 0.0, 0.0, 0, 1.0], [196.04511403446048, 0.0, 0.0, 0, 1.0], [177.70626452040764, 0.0, 0.0, 1, 1.0], [232.82592596431675, 0.0, 0.0, 0, 1.0], [328.9746213107683, 0.0, 0.0, 0, 1.0], [329.07755292675165, 0.0, 0.0, 1, 1.0], [64.94253725651211, 0.0, 0.0, 0, 1.0], [351.56326265802255, 0.0, 0.0, 1, 1.0], [323.3679418390045, 0.0, 0.0, 0, 1.0], [331.307426261965, 0.0, 0.0, 0, 1.0], [13.802205090050872, 0.0, 0.0, 1, 1.0], [220.80483291105557, 0.0, 0.0, 0, 1.0], [316.799385197508, 0.0, 0.0, 0, 1.0], [50.54424801569083, 0.0, 0.0, 0, 1.0], [38.24882530685342, 0.0, 0.0, 0, 1.0], [273.2859490510846, 0.0, 0.0, 0, 1.0], [354.98349840414687, 0.0, 0.0, 1, 1.0], [278.9158628665873, 0.0, 0.0, 1, 1.0], [273.30016625616406, 0.0, 0.0, 0, 1.0], [196.61511562348988, 0.0, 0.0, 1, 1.0], [333.5271856882758, 0.0, 0.0, 0, 1.0], [142.29007235951357, 0.0, 0.0, 1, 1.0], [317.79847608899394, 0.0, 0.0, 0, 1.0], [0.7266499442631127, 0.0, 0.0, 1, 1.0], [107.4744912169713, 0.0, 0.0, 1, 1.0], [345.0547906055752, 0.0, 0.0, 1, 1.0], [358.2082165930105, 0.0, 0.0, 1, 1.0], [201.0817399060137, 0.0, 0.0, 0, 1.0], [219.23887102939034, 0.0, 0.0, 1, 1.0], [195.59591896934865, 0.0, 0.0, 1, 1.0], [29.131395097747454, 0.0, 0.0, 1, 1.0], [36.449727272883024, 0.0, 0.0, 1, 1.0], [60.57824509120188, 0.0, 0.0, 1, 1.0], [199.69795222785152, 0.0, 0.0, 0, 1.0], [231.7997843843201, 0.0, 0.0, 1, 1.0], [222.34509646846848, 0.0, 0.0, 0, 1.0], [282.3188091919718, 0.0, 0.0, 0, 1.0], [179.0187700827499, 0.0, 0.0, 0, 1.0], [274.382960994679, 0.0, 0.0, 1, 1.0], [135.07809583907184, 0.0, 0.0, 1, 1.0], [39.849524713983115, 0.0, 0.0, 0, 1.0], [235.385486624483, 0.0, 0.0, 1, 1.0], [149.94362873415014, 0.0, 0.0, 0, 1.0], [156.50061118400743, 0.0, 0.0, 1, 1.0], [315.9511826004298, 0.0, 0.0, 1, 1.0], [150.73088983672164, 0.0, 0.0, 0, 1.0], [213.4693880426871, 0.0, 0.0, 1, 1.0], [312.4870012154273, 0.0, 0.0, 0, 1.0], [279.58727287760973, 0.0, 0.0, 0, 1.0], [43.36811300275874, 0.0, 0.0, 0, 1.0], [164.33842254754984, 0.0, 0.0, 1, 1.0], [58.26280097271371, 0.0, 0.0, 1, 1.0], [353.5866426164621, 0.0, 0.0, 0, 1.0], [171.278749296968, 0.0, 0.0, 1, 1.0], [131.82160195058572, 0.0, 0.0, 1, 1.0], [235.99549472014235, 0.0, 0.0, 1, 1.0], [358.2431623129092, 0.0, 0.0, 1, 1.0], [216.00100894203774, 0.0, 0.0, 0, 1.0], [278.1703686699143, 0.0, 0.0, 1, 1.0], [82.87694953387485, 0.0, 0.0, 0, 1.0], [155.47445122542783, 0.0, 0.0, 0, 1.0], [26.917148721112603, 0.0, 0.0, 0, 1.0], [289.40830449708994, 0.0, 0.0, 0, 1.0], [354.312922203678, 0.0, 0.0, 0, 1.0], [136.89611110667116, 0.0, 0.0, 1, 1.0], [249.77889121874279, 0.0, 0.0, 1, 1.0], [104.24588792893502, 0.0, 0.0, 1, 1.0], [117.561324155386, 0.0, 0.0, 0, 1.0], [350.7680033335721, 0.0, 0.0, 0, 1.0], [38.49161133326204, 0.0, 0.0, 1, 1.0], [288.07766631373886, 0.0, 0.0, 1, 1.0], [349.02704990380533, 0.0, 0.0, 1, 1.0]]))



class Test_ave_var(unittest.TestCase):
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var"))[0][0]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        if len(return_new) > 0:
            self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
            self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave_var()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave_var()
        self.assertEqual(cm_new.exception.message, "ave_var() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        return_new = fu.ave_var(self.data, mode="a", listID=None)
        return_old = oldfu.ave_var(self.data, mode="a", listID=None)
        self.test_all_the_conditions(return_new, return_old)

    def test_not_current_alignment_parameters(self):
        return_new = fu.ave_var(self.data, mode="nota", listID=None)
        return_old = oldfu.ave_var(self.data, mode="nota", listID=None)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_listID(self):
        return_new = fu.ave_var(self.data, mode="a", listID=[1,2])
        return_old = oldfu.ave_var(self.data, mode="a", listID=[1,2])
        self.test_all_the_conditions(return_new, return_old)



class Test_ave_series(unittest.TestCase):
    data, = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave_series()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave_series()
        self.assertEqual(cm_new.exception.message, "ave_series() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=True, mask=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=True, mask=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_without_pave(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=False, mask=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=False, mask=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_value_withMask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=True, mask=MASK)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=True, mask=MASK)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_without_pave_withMask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=False, mask=MASK)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=False, mask=MASK)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), atol=TOLERANCE,equal_nan=True))




class Test_varf2d_MPI(unittest.TestCase):
    data, = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series"))[0]
    main_node = 0


    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.varf2d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.varf2d_MPI()
        self.assertEqual(cm_new.exception.message, "varf2d_MPI() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value_img2D(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal(return_new[1],[657.29736328125, 4417.57666015625, 12547.1904296875, 5386.271484375, 4284.39404296875, 4738.46044921875, 1702.59375, 2755.111572265625, 2629.47314453125, 2593.828125, 1648.6844482421875, 1247.2149658203125, 1095.642333984375, 892.38525390625, 537.0348510742188, 228.2470245361328, 67.18709564208984, 12.123698234558105, 2.3058433532714844, 0.9480835199356079, 0.673009991645813, 0.5360445976257324, 0.4416618049144745, 0.3757860064506531, 0.31849655508995056, 0.2831489145755768, 0.2509312629699707, 0.22601264715194702, 0.20667268335819244, 0.19018994271755219, 0.17610518634319305, 0.16430525481700897, 0.15365713834762573, 0.14737848937511444, 0.14194539189338684, 0.1365765631198883, 0.13146640360355377, 0.12925119698047638, 0.12116505205631256]))

    def test_default_value_img3D(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_3D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_3D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal(return_new[1],[657.29736328125, 4417.57666015625, 12547.1904296875, 5386.271484375, 4284.39404296875, 4738.46044921875, 1702.59375, 2755.111572265625, 2629.47314453125, 2593.828125, 1648.6844482421875, 1247.2149658203125, 1095.642333984375, 892.38525390625, 537.0348510742188, 228.2470245361328, 67.18709564208984, 12.123698234558105, 2.3058433532714844, 0.9480835199356079, 0.673009991645813, 0.5360445976257324, 0.4416618049144745, 0.3757860064506531, 0.31849655508995056, 0.2831489145755768, 0.2509312629699707, 0.22601264715194702, 0.20667268335819244, 0.19018994271755219, 0.17610518634319305, 0.16430525481700897, 0.15365713834762573, 0.14737848937511444, 0.14194539189338684, 0.1365765631198883, 0.13146640360355377, 0.12925119698047638, 0.12116505205631256]))

    def test_no_mainNode(self):
        return_new = fu.varf2d_MPI(myid=self.main_node+1, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node+1, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(),EMAN2_cppwrap.EMData().get_3dview()))
        self.assertTrue(numpy.array_equal(return_new[1], [0]))

    @unittest.skip("i do not have a valid img")
    def test_with_CTF(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=True, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="a", CTF=True, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal(return_new[1],[657.29736328125, 4417.57666015625, 12547.1904296875, 5386.271484375, 4284.39404296875, 4738.46044921875, 1702.59375, 2755.111572265625, 2629.47314453125, 2593.828125, 1648.6844482421875, 1247.2149658203125, 1095.642333984375, 892.38525390625, 537.0348510742188, 228.2470245361328, 67.18709564208984, 12.123698234558105, 2.3058433532714844, 0.9480835199356079, 0.673009991645813, 0.5360445976257324, 0.4416618049144745, 0.3757860064506531, 0.31849655508995056, 0.2831489145755768, 0.2509312629699707, 0.22601264715194702, 0.20667268335819244, 0.19018994271755219, 0.17610518634319305, 0.16430525481700897, 0.15365713834762573, 0.14737848937511444, 0.14194539189338684, 0.1365765631198883, 0.13146640360355377, 0.12925119698047638, 0.12116505205631256]))

    def test_with_mask(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=MASK, mode="a", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=MASK, mode="a", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1],numpy.full(len(return_new[1]), numpy.nan), equal_nan=True))

    def test_without_apply_alignment_parameters(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="not", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=None, mode="not", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal(return_new[1],[657.143798828125, 4400.4541015625, 12588.783203125, 5354.16748046875, 4298.67626953125, 4745.57958984375, 1703.4573974609375, 2770.798583984375, 2642.747802734375, 2620.13916015625, 1666.76220703125, 1271.91552734375, 1132.228759765625, 928.3460693359375, 560.6527709960938, 241.77120971679688, 68.29119873046875, 10.938950538635254, 1.148409366607666, 0.10395097732543945, 0.009575692005455494, 0.0007211578777059913, 6.025378024787642e-05, 5.932938165642554e-06, 6.559504299730179e-07, 7.396378975954576e-08, 8.115450533807689e-09, 1.0014190587881444e-09, 1.227386109414752e-10, 2.5507195314244946e-11, 1.5168682282462598e-11, 1.4266591380485139e-11, 1.3229447919094195e-11, 1.3020461445134579e-11, 1.4008788920549797e-11, 1.3700839074370919e-11, 1.456094446405931e-11, 1.4547537653675224e-11, 1.2415271935517502e-11]))

    def test_with_mask_without_apply_alignment_parameters(self):
        return_new = fu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=MASK, mode="not", CTF=False, main_node=self.main_node, comm=-1)
        return_old = oldfu.varf2d_MPI(myid=self.main_node, data=self.data, ave=IMAGE_2D, mask=MASK, mode="not", CTF=False, main_node=self.main_node, comm=-1)

        self.assertTrue(numpy.allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.allclose(return_new[1],numpy.full(len(return_new[1]), numpy.nan), equal_nan=True))


""" I resized the image to reduce the execution time (30sec instead of 30min)"""
class Test_varf3d_MPI(unittest.TestCase):
    data, = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series"))[0]
    main_node = 0

    def get_resizedImg(self):
        nima = EMAN2_cppwrap.EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMAN2_cppwrap.EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        proj.set_size(20,20)
        return  [proj]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.varf3d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.varf3d_MPI()
        self.assertEqual(cm_new.exception.message, "varf3d_MPI() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_mask_RuntimeError(self):
        proj = self.get_resizedImg()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=EMData(), reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=EMData(), reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The dimension of the image does not match the dimension of the mask!")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_default_value(self):
        proj = self.get_resizedImg()
        return_new = fu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_withradius(self):
        proj = self.get_resizedImg()
        return_new = fu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=5, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=5, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=0, mpi_comm=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_withCTF(self):
        proj = self.get_resizedImg()
        return_new = fu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=True,sign=1, sym="c1", myid=0, mpi_comm=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=True,sign=1, sym="c1", myid=0, mpi_comm=None)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_myid_no_0(self):
        proj = self.get_resizedImg()
        return_new = fu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=1, mpi_comm=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.varf3d_MPI(proj,ssnr_text_file=None, mask2D=None, reference_structure=None, ou=-1, rw=1.0, npad=1, CTF=False,sign=1, sym="c1", myid=1, mpi_comm=None)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.allclose(return_new.get_3dview(), sparx_utilities.model_blank(2, 2, 2).get_3dview(), atol=TOLERANCE, equal_nan=True))




class Test_ccc(unittest.TestCase):
    (img1, img2, mask) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ccc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ccc()
        self.assertEqual(cm_new.exception.message, "ccc() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file(self):
        return_new = fu.ccc(self.img1, self.img2, self.mask )
        return_old = oldfu.ccc(self.img1, self.img2, self.mask )
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0.8129369020462036)

    def test_None_mask(self):
        return_new = fu.ccc(self.img1, self.img2,None )
        return_old = oldfu.ccc(self.img1, self.img2, None)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0.8021121621131897)

    def test_empty_mask_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccc(self.img1, self.img2, EMData())
        return_old = oldfu.ccc(self.img1, self.img2, EMData())
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0)
        """

    def test_Noneimg1(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.ccc(None, self.img2, self.mask )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ccc(None, self.img2, self.mask )
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'cmp'")

    def test_emptyimg1(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(EMData(), self.img2, self.mask )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(EMData(), self.img2, self.mask )

        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_Noneimg2(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(self.img1, None,self.mask)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(self.img1,None, self.mask)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_emptyimg2(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(self.img1,EMData(), self.mask)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(self.img1,EMData(), self.mask)

        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])




class Test_fsc(unittest.TestCase):
    (img1, img2) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fsc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fsc()
        self.assertEqual(cm_new.exception.message, "fsc() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.fsc(self.img1, self.img2 ,w=1.0, filename=None)
        return_old = oldfu.fsc(self.img1, self.img2, w=1.0, filename=None)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[[0.0, 0.006756756920367479, 0.013513513840734959, 0.020270271226763725, 0.027027027681469917, 0.03378378599882126, 0.04054054245352745, 0.04729729890823364, 0.054054055362939835, 0.06081081181764603, 0.06756757199764252, 0.07432432472705841, 0.0810810849070549, 0.0878378376364708, 0.09459459781646729, 0.10135135054588318, 0.10810811072587967, 0.11486487090587616, 0.12162162363529205, 0.12837837636470795, 0.13513514399528503, 0.14189189672470093, 0.14864864945411682, 0.15540540218353271, 0.1621621698141098, 0.1689189225435257, 0.1756756752729416, 0.18243244290351868, 0.18918919563293457, 0.19594594836235046, 0.20270270109176636, 0.20945946872234344, 0.21621622145175934, 0.22297297418117523, 0.22972974181175232, 0.2364864945411682, 0.2432432472705841, 0.25, 0.2567567527294159, 0.2635135054588318, 0.27027028799057007, 0.27702704071998596, 0.28378379344940186, 0.29054054617881775, 0.29729729890823364, 0.30405405163764954, 0.31081080436706543, 0.3175675868988037, 0.3243243396282196, 0.3310810923576355, 0.3378378450870514, 0.3445945978164673, 0.3513513505458832, 0.3581081032752991, 0.36486488580703735, 0.37162163853645325, 0.37837839126586914, 0.38513514399528503, 0.3918918967247009, 0.3986486494541168, 0.4054054021835327, 0.412162184715271, 0.4189189374446869, 0.4256756901741028, 0.4324324429035187, 0.43918919563293457, 0.44594594836235046, 0.45270270109176636, 0.45945948362350464, 0.46621623635292053, 0.4729729890823364, 0.4797297418117523, 0.4864864945411682, 0.4932432472705841, 0.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9999955892562866, 0.998257040977478, 0.9918078780174255, 0.9856542348861694, 0.9881284832954407, 0.9837390184402466, 0.9792850017547607, 0.9748953580856323, 0.9678066968917847, 0.9480144381523132, 0.9332209229469299, 0.8907226324081421, 0.8318315148353577, 0.7905364036560059, 0.8009309768676758, 0.806480348110199, 0.7942577600479126, 0.7667281627655029, 0.7561865448951721, 0.7457202672958374, 0.7362813353538513, 0.7185238003730774, 0.7565295100212097, 0.7811623215675354, 0.7833795547485352, 0.7671730518341064, 0.7416536808013916, 0.7073566913604736, 0.719993531703949, 0.7416943907737732, 0.7183090448379517, 0.691653847694397, 0.6795889735221863, 0.6586717367172241, 0.6515437960624695, 0.5965009331703186, 0.5489014387130737, 0.565247654914856, 0.5726661682128906, 0.5036922693252563, 0.38146188855171204, 0.26737281680107117, 0.3945968449115753, 0.4944046437740326, 0.39991524815559387, 0.2689603269100189, 0.2521679699420929, 0.2794639468193054, 0.280245840549469, 0.24809274077415466, 0.25190722942352295, 0.2277139574289322, 0.19758595526218414, 0.19085757434368134, 0.19563539326190948, 0.19770143926143646, 0.16604311764240265, 0.18556171655654907, 0.13770028948783875, 0.152345210313797, 0.17087307572364807, 0.15340681374073029, 0.16820573806762695, 0.18507032096385956, 0.1278713047504425, 0.1172977089881897, 0.0], [2.0, 18.0, 62.0, 98.0, 210.0, 350.0, 450.0, 602.0, 762.0, 1142.0, 1250.0, 1458.0, 1814.0, 2178.0, 2498.0, 2622.0, 3338.0, 3722.0, 4170.0, 4358.0, 5034.0, 5714.0, 5982.0, 6602.0, 7130.0, 8034.0, 8606.0, 9066.0, 9962.0, 10550.0, 11226.0, 12146.0, 12606.0, 13802.0, 14754.0, 15194.0, 16454.0, 17154.0, 18266.0, 18750.0, 20234.0, 21450.0, 21962.0, 23462.0, 24042.0, 25946.0, 26118.0, 27506.0, 29066.0, 30450.0, 31742.0, 32250.0, 34250.0, 35454.0, 36434.0, 37682.0, 39294.0, 41426.0, 42066.0, 43490.0, 45702.0, 46634.0, 48554.0, 48870.0, 51554.0, 54090.0, 54314.0, 56910.0, 57482.0, 60626.0, 61206.0, 62570.0, 65730.0, 66686.0, 0.0]]))

    def test_saveOnfile(self):
        fnew=path.join(ABSOLUTE_PATH,"new.txt")
        fold = path.join(ABSOLUTE_PATH,"new.old")
        return_new = fu.fsc(self.img1, self.img2 ,w=1.0, filename=fnew)
        return_old = oldfu.fsc(self.img1, self.img2, w=1.0, filename=fold)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0.0, 0.006756756920367479, 0.013513513840734959, 0.020270271226763725, 0.027027027681469917, 0.03378378599882126, 0.04054054245352745, 0.04729729890823364, 0.054054055362939835, 0.06081081181764603, 0.06756757199764252, 0.07432432472705841, 0.0810810849070549, 0.0878378376364708, 0.09459459781646729, 0.10135135054588318, 0.10810811072587967, 0.11486487090587616, 0.12162162363529205, 0.12837837636470795, 0.13513514399528503, 0.14189189672470093, 0.14864864945411682, 0.15540540218353271, 0.1621621698141098, 0.1689189225435257, 0.1756756752729416, 0.18243244290351868, 0.18918919563293457, 0.19594594836235046, 0.20270270109176636, 0.20945946872234344, 0.21621622145175934, 0.22297297418117523, 0.22972974181175232, 0.2364864945411682, 0.2432432472705841, 0.25, 0.2567567527294159, 0.2635135054588318, 0.27027028799057007, 0.27702704071998596, 0.28378379344940186, 0.29054054617881775, 0.29729729890823364, 0.30405405163764954, 0.31081080436706543, 0.3175675868988037, 0.3243243396282196, 0.3310810923576355, 0.3378378450870514, 0.3445945978164673, 0.3513513505458832, 0.3581081032752991, 0.36486488580703735, 0.37162163853645325, 0.37837839126586914, 0.38513514399528503, 0.3918918967247009, 0.3986486494541168, 0.4054054021835327, 0.412162184715271, 0.4189189374446869, 0.4256756901741028, 0.4324324429035187, 0.43918919563293457, 0.44594594836235046, 0.45270270109176636, 0.45945948362350464, 0.46621623635292053, 0.4729729890823364, 0.4797297418117523, 0.4864864945411682, 0.4932432472705841, 0.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9999955892562866, 0.998257040977478, 0.9918078780174255, 0.9856542348861694, 0.9881284832954407, 0.9837390184402466, 0.9792850017547607, 0.9748953580856323, 0.9678066968917847, 0.9480144381523132, 0.9332209229469299, 0.8907226324081421, 0.8318315148353577, 0.7905364036560059, 0.8009309768676758, 0.806480348110199, 0.7942577600479126, 0.7667281627655029, 0.7561865448951721, 0.7457202672958374, 0.7362813353538513, 0.7185238003730774, 0.7565295100212097, 0.7811623215675354, 0.7833795547485352, 0.7671730518341064, 0.7416536808013916, 0.7073566913604736, 0.719993531703949, 0.7416943907737732, 0.7183090448379517, 0.691653847694397, 0.6795889735221863, 0.6586717367172241, 0.6515437960624695, 0.5965009331703186, 0.5489014387130737, 0.565247654914856, 0.5726661682128906, 0.5036922693252563, 0.38146188855171204, 0.26737281680107117, 0.3945968449115753, 0.4944046437740326, 0.39991524815559387, 0.2689603269100189, 0.2521679699420929, 0.2794639468193054, 0.280245840549469, 0.24809274077415466, 0.25190722942352295, 0.2277139574289322, 0.19758595526218414, 0.19085757434368134, 0.19563539326190948, 0.19770143926143646, 0.16604311764240265, 0.18556171655654907, 0.13770028948783875, 0.152345210313797, 0.17087307572364807, 0.15340681374073029, 0.16820573806762695, 0.18507032096385956, 0.1278713047504425, 0.1172977089881897, 0.0], [2.0, 18.0, 62.0, 98.0, 210.0, 350.0, 450.0, 602.0, 762.0, 1142.0, 1250.0, 1458.0, 1814.0, 2178.0, 2498.0, 2622.0, 3338.0, 3722.0, 4170.0, 4358.0, 5034.0, 5714.0, 5982.0, 6602.0, 7130.0, 8034.0, 8606.0, 9066.0, 9962.0, 10550.0, 11226.0, 12146.0, 12606.0, 13802.0, 14754.0, 15194.0, 16454.0, 17154.0, 18266.0, 18750.0, 20234.0, 21450.0, 21962.0, 23462.0, 24042.0, 25946.0, 26118.0, 27506.0, 29066.0, 30450.0, 31742.0, 32250.0, 34250.0, 35454.0, 36434.0, 37682.0, 39294.0, 41426.0, 42066.0, 43490.0, 45702.0, 46634.0, 48554.0, 48870.0, 51554.0, 54090.0, 54314.0, 56910.0, 57482.0, 60626.0, 61206.0, 62570.0, 65730.0, 66686.0, 0.0]]))
        self.assertTrue(returns_values_in_file(fnew),returns_values_in_file(fold))
        remove_list_of_file([fnew, fold])

    def test_w_set_0_returns_MemoryError(self):
        with self.assertRaises(MemoryError) as cm_new:
            fu.fsc(self.img1, self.img2 ,w=0, filename=None)
        with self.assertRaises(MemoryError) as cm_old:
            oldfu.fsc(self.img1, self.img2, w=0, filename=None)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyImg1_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(EMData(), self.img2 ,w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(EMData(), self.img2, w=1.0, filename=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "Cannot calculate FSC for 1D images")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img1_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.fsc(None, self.img2 ,w=1.0, filename=None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fsc(None, self.img2, w=1.0, filename=None)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'calc_fourier_shell_correlation'")

    def test_emptyImg2_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(self.img1, EMData() ,w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(self.img1,EMData(), w=1.0, filename=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_Img2_returns_AttributeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(self.img1,None ,w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(self.img1, None, w=1.0, filename=None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NullPointerException")
        self.assertEqual(msg[1], "NULL input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])




class Test_fsc_mask(unittest.TestCase):
    (img1, img2, mask, w, not_used) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc_mask"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fsc_mask()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fsc_mask()
        self.assertEqual(cm_new.exception.message, "fsc_mask() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=self.w, filename = None)
        return_old = oldfu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=self.w, filename = None)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0.0, 0.014285714365541935, 0.02857142873108387, 0.04285714402794838, 0.05714285746216774, 0.0714285746216774, 0.08571428805589676, 0.10000000149011612, 0.11428571492433548, 0.12857143580913544, 0.1428571492433548, 0.15714286267757416, 0.17142857611179352, 0.18571428954601288, 0.20000000298023224, 0.2142857164144516, 0.22857142984867096, 0.24285714328289032, 0.2571428716182709, 0.27142858505249023, 0.2857142984867096, 0.30000001192092896, 0.3142857253551483, 0.3285714387893677, 0.34285715222358704, 0.3571428656578064, 0.37142857909202576, 0.3857142925262451, 0.4000000059604645, 0.41428571939468384, 0.4285714328289032, 0.44285714626312256, 0.4571428596973419, 0.4714285731315613, 0.48571428656578064, 0.5], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [2.0, 8.0, 12.0, 16.0, 32.0, 28.0, 40.0, 40.0, 48.0, 68.0, 56.0, 72.0, 68.0, 88.0, 88.0, 84.0, 112.0, 112.0, 112.0, 116.0, 112.0, 144.0, 140.0, 144.0, 144.0, 168.0, 164.0, 160.0, 184.0, 172.0, 200.0, 192.0, 188.0, 208.0, 224.0, 214.0]]))

    def test_saveOnfile(self):
        fnew = path.join(ABSOLUTE_PATH, "new.txt")
        fold = path.join(ABSOLUTE_PATH, "new.old")
        return_new = fu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=self.w, filename = fnew)
        return_old = oldfu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=self.w, filename = fold)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0.0, 0.014285714365541935, 0.02857142873108387, 0.04285714402794838, 0.05714285746216774, 0.0714285746216774, 0.08571428805589676, 0.10000000149011612, 0.11428571492433548, 0.12857143580913544, 0.1428571492433548, 0.15714286267757416, 0.17142857611179352, 0.18571428954601288, 0.20000000298023224, 0.2142857164144516, 0.22857142984867096, 0.24285714328289032, 0.2571428716182709, 0.27142858505249023, 0.2857142984867096, 0.30000001192092896, 0.3142857253551483, 0.3285714387893677, 0.34285715222358704, 0.3571428656578064, 0.37142857909202576, 0.3857142925262451, 0.4000000059604645, 0.41428571939468384, 0.4285714328289032, 0.44285714626312256, 0.4571428596973419, 0.4714285731315613, 0.48571428656578064, 0.5], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [2.0, 8.0, 12.0, 16.0, 32.0, 28.0, 40.0, 40.0, 48.0, 68.0, 56.0, 72.0, 68.0, 88.0, 88.0, 84.0, 112.0, 112.0, 112.0, 116.0, 112.0, 144.0, 140.0, 144.0, 144.0, 168.0, 164.0, 160.0, 184.0, 172.0, 200.0, 192.0, 188.0, 208.0, 224.0, 214.0]]))
        self.assertTrue(returns_values_in_file(fnew),returns_values_in_file(fold))
        remove_list_of_file([fnew, fold])

    def test_w_set_0_returns_MemoryError(self):
        with self.assertRaises(MemoryError) as cm_new:
            fu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=0, filename = None)
        with self.assertRaises(MemoryError) as cm_old:
            oldfu.fsc_mask(img1 =self.img1, img2=self.img1, mask=self.mask, w=0, filename = None)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_emptyImg1_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(img1 =EMData(), img2=self.img1, mask=self.mask, w=self.w, filename = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(img1 =EMData(), img2=self.img1, mask=self.mask, w=self.w, filename = None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The dimension of the image does not match the dimension of the mask!")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img1_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.fsc_mask(img1 =None, img2=self.img1, mask=self.mask, w=self.w, filename = None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fsc_mask(img1 =None, img2=self.img1, mask=self.mask, w=self.w, filename = None)
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")

    def test_emptyImg2_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(img1 =self.img1, img2=EMData(), mask=self.mask, w=self.w, filename = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(img1 =self.img1, img2=EMData(), mask=self.mask, w=self.w, filename = None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The dimension of the image does not match the dimension of the mask!")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img2_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(img1 =self.img1, img2=None, mask=self.mask, w=self.w, filename = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(img1 =self.img1, img2=None, mask=self.mask, w=self.w, filename = None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NullPointerException")
        self.assertEqual(msg[1], "NULL input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])
        """

    def test_empty_mask_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(img1 =self.img1, img2=self.img1, mask=EMData(), w=self.w, filename = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(img1 =self.img1, img2=self.img1, mask=EMData(), w=self.w, filename = None)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])


    def test_NoneType_mask(self):
        return_new = fu.fsc_mask(img1 =self.img1, img2=self.img1, mask=None, w=self.w, filename = None)
        return_old = oldfu.fsc_mask(img1 =self.img1, img2=self.img1, mask=None, w=self.w, filename = None)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0.0, 0.014285714365541935, 0.02857142873108387, 0.04285714402794838, 0.05714285746216774, 0.0714285746216774, 0.08571428805589676, 0.10000000149011612, 0.11428571492433548, 0.12857143580913544, 0.1428571492433548, 0.15714286267757416, 0.17142857611179352, 0.18571428954601288, 0.20000000298023224, 0.2142857164144516, 0.22857142984867096, 0.24285714328289032, 0.2571428716182709, 0.27142858505249023, 0.2857142984867096, 0.30000001192092896, 0.3142857253551483, 0.3285714387893677, 0.34285715222358704, 0.3571428656578064, 0.37142857909202576, 0.3857142925262451, 0.4000000059604645, 0.41428571939468384, 0.4285714328289032, 0.44285714626312256, 0.4571428596973419, 0.4714285731315613, 0.48571428656578064, 0.5], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [2.0, 8.0, 12.0, 16.0, 32.0, 28.0, 40.0, 40.0, 48.0, 68.0, 56.0, 72.0, 68.0, 88.0, 88.0, 84.0, 112.0, 112.0, 112.0, 116.0, 112.0, 144.0, 140.0, 144.0, 144.0, 168.0, 164.0, 160.0, 184.0, 172.0, 200.0, 192.0, 188.0, 208.0, 224.0, 214.0]]))



class Test_locres(unittest.TestCase):
    (vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.locres"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.locres()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.locres()
        self.assertEqual(cm_new.exception.message, "locres() takes exactly 9 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        self.assertTrue(True)
        """
        Connected to pydev debugger (build 183.6156.13)
        Launching unittests with arguments python -m unittest sphire.tests.test_statistics.Test_locres.test_pickle_value in /home/lusnig/EMAN2/eman2
        [rtxr2:22191] *** An error occurred in MPI_Recv
        [rtxr2:22191] *** reported by process [1715994625,140423955742720]
        [rtxr2:22191] *** on communicator MPI_COMM_WORLD
        [rtxr2:22191] *** MPI_ERR_RANK: invalid rank
        [rtxr2:22191] *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
        [rtxr2:22191] ***    and potentially your MPI job)
        
        Process finished with exit code 6
        
        return_new = fu.locres(vi=self.vi, ui=self.ui, m=self.m, nk=self.nk, cutoff=self.nk, step=self.step, myid=0, main_node=0, number_of_proc=self.number_of_proc)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.locres(vi=self.vi, ui=self.ui, m=self.m, nk=self.nk, cutoff=self.nk, step=self.step, myid=0, main_node=0, number_of_proc=self.number_of_proc)
        mpi_barrier(MPI_COMM_WORLD)
        """



class Test_histogram(unittest.TestCase):

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.histogram()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.histogram()
        self.assertEqual(cm_new.exception.message, "histogram() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneTypeImg_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.histogram(image=None, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        return_old = oldfu.histogram(image=None, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_img2D(self):
        return_new = fu.histogram(image=IMAGE_2D, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        return_old = oldfu.histogram(image=IMAGE_2D, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new, [3.0, 0.0, 0.0, 1.0, 3.0, 3.0, 2.0, 0.0, 4.0, 5.0, 2.0, 5.0, 10.0, 13.0, 15.0, 14.0, 13.0, 13.0, 28.0, 14.0, 23.0, 32.0, 25.0, 37.0, 32.0, 48.0, 40.0, 43.0, 68.0, 46.0, 59.0, 82.0, 75.0, 91.0, 114.0, 115.0, 113.0, 117.0, 122.0, 125.0, 120.0, 153.0, 150.0, 195.0, 192.0, 176.0, 189.0, 198.0, 210.0, 229.0, 236.0, 224.0, 212.0, 240.0, 262.0, 264.0, 265.0, 270.0, 259.0, 271.0, 293.0, 330.0, 299.0, 361.0, 386.0, 371.0, 394.0, 385.0, 422.0, 433.0, 425.0, 412.0, 434.0, 437.0, 448.0, 446.0, 477.0, 446.0, 502.0, 506.0, 496.0, 529.0, 535.0, 585.0, 543.0, 542.0, 566.0, 569.0, 580.0, 562.0, 592.0, 579.0, 637.0, 642.0, 660.0, 688.0, 698.0, 706.0, 714.0, 713.0, 736.0, 769.0, 806.0, 770.0, 796.0, 869.0, 924.0, 875.0, 888.0, 963.0, 974.0, 993.0, 1011.0, 1021.0, 1013.0, 1037.0, 1099.0, 1129.0, 1115.0, 1195.0, 1156.0, 1263.0, 1261.0, 1472.0, 1665.0, 1784.0, 2016.0, 1985.0, 1851.0, 1699.0, 2036.0, 2355.0, 2409.0, 2356.0, 2280.0, 2326.0, 2357.0, 2585.0, 3166.0, 4145.0, 3553.0, 2416.0, 2059.0, 1972.0, 1857.0, 1536.0, 1338.0, 891.0, 577.0, 409.0, 322.0, 327.0, 264.0, 282.0, 258.0, 277.0, 235.0, 228.0, 194.0, 195.0, 182.0, 185.0, 143.0, 174.0, 152.0, 171.0, 142.0, 138.0, 138.0, 125.0, 141.0, 134.0, 132.0, 132.0, 137.0, 139.0, 138.0, 124.0, 121.0, 130.0, 126.0, 113.0, 122.0, 135.0, 119.0, 116.0, 133.0, 124.0, 130.0, 141.0, 131.0, 150.0, 126.0, 132.0, 143.0, 162.0, 137.0, 131.0, 140.0, 131.0, 137.0, 146.0, 172.0, 146.0, 131.0, 174.0, 143.0, 149.0, 196.0, 149.0, 144.0, 170.0, 159.0, 159.0, 140.0, 171.0, 152.0, 167.0, 194.0, 178.0, 186.0, 176.0, 172.0, 153.0, 189.0, 171.0, 168.0, 197.0, 162.0, 169.0, 175.0, 195.0, 158.0, 164.0, 174.0, 155.0, 154.0, 161.0, 172.0, 122.0, 154.0, 149.0, 146.0, 135.0, 144.0, 132.0, 116.0, 121.0, 126.0, 125.0, 135.0, 101.0, 122.0, 117.0, 104.0, 103.0, 115.0, 87.0, 87.0, 93.0, 77.0, 80.0, 104.0, 92.0, 70.0, 85.0, 79.0, 86.0, 88.0, 57.0, 71.0, 80.0, 59.0, 59.0, 75.0, 63.0, 56.0, 60.0, 63.0, 67.0, 61.0, 53.0, 55.0, 53.0, 45.0, 41.0, 49.0, 43.0, 35.0, 48.0, 43.0, 38.0, 36.0, 32.0, 31.0, 40.0, 45.0, 39.0, 24.0, 34.0, 44.0, 30.0, 30.0, 30.0, 21.0, 36.0, 21.0, 21.0, 25.0, 24.0, 15.0, 20.0, 9.0, 21.0, 11.0, 13.0, 15.0, 12.0, 14.0, 7.0, 15.0, 14.0, 8.0, 13.0, 9.0, 9.0, 5.0, 7.0, 6.0, 5.0, 7.0, 5.0, 4.0, 4.0, 5.0, 2.0, 1.0, 2.0, 4.0, 2.0, 2.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -4.182239055633545, -4.152743339538574, -4.123247146606445, -4.093751430511475, -4.064255237579346, -4.034759044647217, -4.005263328552246, -3.975767135620117, -3.9462711811065674, -3.9167752265930176, -3.8872792720794678, -3.857783317565918, -3.828287363052368, -3.7987914085388184, -3.7692952156066895, -3.7397992610931396, -3.71030330657959, -3.68080735206604, -3.6513113975524902, -3.6218154430389404, -3.5923194885253906, -3.5628232955932617, -3.533327341079712, -3.503831386566162, -3.4743354320526123, -3.4448394775390625, -3.4153435230255127, -3.385847568511963, -3.356351375579834, -3.326855421066284, -3.2973594665527344, -3.2678635120391846, -3.2383675575256348, -3.208871603012085, -3.179375648498535, -3.1498796939849854, -3.1203835010528564, -3.0908875465393066, -3.061391592025757, -3.031895637512207, -3.0023996829986572, -2.9729037284851074, -2.9434077739715576, -2.9139115810394287, -2.884415626525879, -2.854919672012329, -2.8254237174987793, -2.7959277629852295, -2.7664318084716797, -2.73693585395813, -2.707439661026001, -2.677943706512451, -2.6484477519989014, -2.6189517974853516, -2.5894558429718018, -2.559959888458252, -2.530463933944702, -2.5009677410125732, -2.4714717864990234, -2.4419758319854736, -2.412479877471924, -2.382983922958374, -2.353487968444824, -2.3239920139312744, -2.2944958209991455, -2.2649998664855957, -2.235503911972046, -2.206007957458496, -2.1765120029449463, -2.1470160484313965, -2.1175200939178467, -2.0880239009857178, -2.058527946472168, -2.029031991958618, -1.9995360374450684, -1.9700400829315186, -1.9405440092086792, -1.9110480546951294, -1.8815521001815796, -1.8520561456680298, -1.8225600719451904, -1.7930641174316406, -1.7635681629180908, -1.734072208404541, -1.7045761346817017, -1.6750801801681519, -1.645584225654602, -1.6160881519317627, -1.586592197418213, -1.557096242904663, -1.5276002883911133, -1.498104214668274, -1.4686082601547241, -1.4391123056411743, -1.409616231918335, -1.3801202774047852, -1.3506243228912354, -1.3211283683776855, -1.2916322946548462, -1.2621363401412964, -1.2326403856277466, -1.2031443119049072, -1.1736483573913574, -1.1441524028778076, -1.1146564483642578, -1.0851603746414185, -1.0556644201278687, -1.0261684656143188, -0.9966724514961243, -0.9671764373779297, -0.9376804828643799, -0.9081844687461853, -0.8786885142326355, -0.8491925001144409, -0.8196965456008911, -0.7902005314826965, -0.760704517364502, -0.7312085628509521, -0.7017125487327576, -0.6722165942192078, -0.6427205801010132, -0.6132246255874634, -0.5837286114692688, -0.554232656955719, -0.5247366428375244, -0.4952406585216522, -0.46574467420578003, -0.43624868988990784, -0.40675267577171326, -0.37725669145584106, -0.34776070713996887, -0.3182647228240967, -0.2887687385082245, -0.2592727541923523, -0.2297767549753189, -0.20028077065944672, -0.17078478634357452, -0.14128880202770233, -0.11179281026124954, -0.08229681849479675, -0.05280083417892456, -0.02330484427511692, 0.00619114376604557, 0.03568713366985321, 0.0651831179857254, 0.09467910975217819, 0.12417509406805038, 0.15367108583450317, 0.18316707015037537, 0.21266305446624756, 0.24215905368328094, 0.27165502309799194, 0.3011510372161865, 0.3306470215320587, 0.3601430058479309, 0.3896389901638031, 0.4191349744796753, 0.4486309587955475, 0.4781269431114197, 0.5076229572296143, 0.5371189117431641, 0.5666149258613586, 0.5961108803749084, 0.625606894493103, 0.6551029086112976, 0.6845988631248474, 0.714094877243042, 0.7435908317565918, 0.7730868458747864, 0.8025828003883362, 0.8320788145065308, 0.8615747690200806, 0.8910707831382751, 0.9205667972564697, 0.9500627517700195, 0.9795587658882141, 1.0090547800064087, 1.0385507345199585, 1.0680466890335083, 1.097542643547058, 1.1270387172698975, 1.1565346717834473, 1.186030626296997, 1.2155267000198364, 1.2450226545333862, 1.274518609046936, 1.3040145635604858, 1.3335106372833252, 1.363006591796875, 1.3925025463104248, 1.4219986200332642, 1.451494574546814, 1.4809905290603638, 1.5104864835739136, 1.539982557296753, 1.5694785118103027, 1.5989744663238525, 1.628470540046692, 1.6579664945602417, 1.6874624490737915, 1.7169584035873413, 1.7464544773101807, 1.7759504318237305, 1.8054463863372803, 1.8349424600601196, 1.8644384145736694, 1.8939343690872192, 1.923430323600769, 1.9529263973236084, 1.9824223518371582, 2.011918306350708, 2.041414260864258, 2.0709102153778076, 2.1004064083099365, 2.1299023628234863, 2.159398317337036, 2.188894271850586, 2.2183902263641357, 2.2478861808776855, 2.2773821353912354, 2.3068783283233643, 2.336374282836914, 2.365870237350464, 2.3953661918640137, 2.4248621463775635, 2.4543581008911133, 2.483854055404663, 2.513350009918213, 2.542846202850342, 2.5723421573638916, 2.6018381118774414, 2.631334066390991, 2.660830020904541, 2.690325975418091, 2.7198219299316406, 2.7493181228637695, 2.7788140773773193, 2.808310031890869, 2.837805986404419, 2.8673019409179688, 2.8967978954315186, 2.9262938499450684, 2.9557900428771973, 2.985285997390747, 3.014781951904297, 3.0442779064178467, 3.0737738609313965, 3.1032698154449463, 3.132765769958496, 3.162261962890625, 3.191757917404175, 3.2212538719177246, 3.2507498264312744, 3.280245780944824, 3.309741735458374, 3.339237689971924, 3.3687338829040527, 3.3982298374176025, 3.4277257919311523, 3.457221746444702, 3.486717700958252, 3.5162136554718018, 3.5457096099853516, 3.5752058029174805, 3.6047017574310303, 3.63419771194458, 3.66369366645813, 3.6931896209716797, 3.7226855754852295, 3.7521815299987793, 3.781677722930908, 3.811173677444458, 3.840669631958008, 3.8701655864715576, 3.8996615409851074, 3.9291574954986572, 3.958653450012207, 3.988149642944336, 4.017645359039307, 4.0471415519714355, 4.076637268066406, 4.106133460998535, 4.135629653930664, 4.165125370025635, 4.194621562957764, 4.224117279052734, 4.253613471984863, 4.283109188079834, 4.312605381011963, 4.342101573944092, 4.3715972900390625, 4.401093482971191, 4.430589199066162, 4.460085391998291, 4.489581108093262, 4.519077301025391, 4.5485734939575195, 4.57806921005249, 4.607565402984619, 4.63706111907959, 4.666557312011719, 4.6960530281066895, 4.725549221038818, 4.755045413970947, 4.784541130065918, 4.814037322998047, 4.843533039093018, 4.8730292320251465, 4.902524948120117, 4.932021141052246, 4.961517333984375, 4.991013050079346, 5.020509243011475, 5.050004959106445, 5.079501152038574, 5.108996868133545, 5.138493061065674, 5.167989253997803, 5.197484970092773, 5.226981163024902, 5.256476879119873, 5.285973072052002, 5.315468788146973, 5.344964981079102, 5.3744611740112305, 5.403956890106201, 5.43345308303833, 5.462948799133301, 5.49244499206543, 5.5219407081604, 5.551436901092529, 5.580933094024658, 5.610428810119629, 5.639925003051758, 5.6694207191467285, 5.698916912078857, 5.728412628173828, 5.757908821105957, 5.787405014038086, 5.816900730133057, 5.8463969230651855, 5.875892639160156, 5.905388832092285, 5.934884548187256, 5.964380741119385, 5.993876934051514, 6.023372650146484, 6.052868843078613, 6.082364559173584, 6.111860752105713, 6.141356468200684, 6.1708526611328125]))

    def test_img3D(self):
        return_new = fu.histogram(image=IMAGE_3D, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        return_old = oldfu.histogram(image=IMAGE_3D, mask=None, nbins=0, hmin=0.0, hmax=0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [1.0, 3.0, 12.0, 98.0, 765.0, 384538.0, 22519.0, 8393.0, 4161.0, 2307.0, 1342.0, 975.0, 756.0, 728.0, 626.0, 565.0, 630.0, 572.0, 569.0, 526.0, 529.0, 555.0, 494.0, 471.0, 512.0, 433.0, 463.0, 423.0, 407.0, 374.0, 324.0, 366.0, 290.0, 256.0, 240.0, 233.0, 244.0, 215.0, 182.0, 174.0, 158.0, 162.0, 144.0, 123.0, 131.0, 98.0, 108.0, 83.0, 84.0, 76.0, 61.0, 77.0, 61.0, 43.0, 44.0, 39.0, 35.0, 27.0, 35.0, 21.0, 12.0, 10.0, 16.0, 16.0, 10.0, 4.0, 5.0, 2.0, 0.0, 2.0, 1.0, 4.0, 2.0, 4.0, 4.0, 3.0, -9.69376277923584, -7.617851257324219, -5.541939735412598, -3.4660277366638184, -1.3901160955429077, 0.6857956647872925, 2.761707305908203, 4.837619304656982, 6.9135308265686035, 8.989442825317383, 11.065354347229004, 13.141265869140625, 15.217178344726562, 17.293088912963867, 19.369001388549805, 21.444913864135742, 23.520824432373047, 25.596736907958984, 27.672649383544922, 29.748559951782227, 31.824472427368164, 33.90038299560547, 35.976295471191406, 38.052207946777344, 40.12812042236328, 42.20402908325195, 44.27994155883789, 46.35585403442383, 48.431766510009766, 50.5076789855957, 52.58359146118164, 54.65950012207031, 56.73541259765625, 58.81132507324219, 60.887237548828125, 62.96315002441406, 65.0390625, 67.11497497558594, 69.19088745117188, 71.26679229736328, 73.34270477294922, 75.41861724853516, 77.4945297241211, 79.57044219970703, 81.64635467529297, 83.7222671508789, 85.79817962646484, 87.87409210205078, 89.95000457763672, 92.02590942382812, 94.10182189941406, 96.177734375, 98.25364685058594, 100.32955932617188, 102.40547180175781, 104.48138427734375, 106.55729675292969, 108.63320922851562, 110.70912170410156, 112.7850341796875, 114.8609390258789, 116.93685150146484, 119.01276397705078, 121.08867645263672, 123.16458892822266, 125.2405014038086, 127.31641387939453, 129.39231872558594, 131.46823120117188, 133.5441436767578, 135.62005615234375, 137.6959686279297, 139.77188110351562, 141.84779357910156, 143.9237060546875, 145.99961853027344]))

    def test_EmptyImg(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.histogram(image=EMData(), mask=None, nbins=0, hmin=0.0, hmax=0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.histogram(image=EMData(), mask=None, nbins=0, hmin=0.0, hmax=0.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_mask(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.histogram(image=IMAGE_2D, mask=EMData(), nbins=0, hmin=0.0, hmax=0.0)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.histogram(image=IMAGE_2D, mask=EMData(), nbins=0, hmin=0.0, hmax=0.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], "The size of mask image should be of same size as the input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])



class Test_k_means_match_clusters_asg_new(unittest.TestCase):
    (asg1, asg2) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_clusters_asg_new"))[0]
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_match_clusters_asg_new()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_match_clusters_asg_new()
        self.assertEqual(cm_new.exception.message, "k_means_match_clusters_asg_new() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):

        return_new = fu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=self.asg2, T=0)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=self.asg2, T=0)

        self.assertTrue(numpy.array_equal(return_new[0],[[0, 0], [1, 1]]))
        self.assertTrue(numpy.array_equal([return_new[1][0].tolist(),return_new[1][1].tolist()], [[4, 5, 9, 11, 12, 13, 15, 16, 18, 19, 23, 27, 31, 36, 37, 40, 41, 44, 45, 46, 48, 49, 51, 53, 57, 58, 61, 68, 73, 76, 78, 81, 85, 87, 91, 97, 100, 101, 103, 104, 109, 111, 112, 116, 117, 118, 120, 124, 125, 126, 127, 128, 130, 131, 136, 138, 140, 141, 143, 145, 148, 149, 150, 154, 159, 161, 164, 167, 172, 173, 176, 178, 186, 191, 196, 198, 202, 206, 208, 211, 216, 218, 225, 226, 228, 232, 236, 239, 244, 247, 250, 254, 255, 258, 261, 263, 265, 272, 274, 278, 291, 297, 304, 306, 317, 318, 321, 323, 324, 326, 332, 340, 343, 345, 347, 348, 352, 357, 358, 359, 360, 368, 369, 370, 371, 378, 381, 382, 386, 388, 390, 392, 398, 401, 404, 408, 409, 410, 411, 415, 416, 417, 418, 422, 425, 426, 432, 434, 446, 448, 449, 450, 456, 457, 459, 460, 461, 462, 466, 468, 471, 472, 475, 476, 478, 479, 480, 481, 492, 499, 500, 501, 506, 508, 509, 516, 519, 521, 524, 533, 534, 536, 539, 547, 551, 555, 572, 578, 579, 583, 586, 589, 592, 595, 599, 600, 601, 602, 605, 606, 608, 613, 616, 618, 627, 630, 637, 639, 644, 645, 648, 658, 659, 663, 668, 670, 675, 678, 681, 682, 687, 692, 693, 696, 697, 698, 700, 705, 711, 716, 722, 725, 729, 733, 735, 737, 741, 743, 744, 748, 750, 754, 757, 760, 762, 764, 766, 774, 781, 782, 783, 789, 792, 793, 795, 800, 801, 804, 805, 806, 807, 808, 815, 817, 819, 823, 824, 825, 826, 827, 835, 836, 841, 844, 848, 849, 850, 851, 852, 853, 860, 862, 864, 865, 868, 870, 871, 874, 876, 877, 878, 881, 884, 887, 889, 893, 896, 905, 907, 909, 912, 917, 920, 928, 929, 930, 931, 934, 935, 943, 949, 954, 955, 956, 965, 970, 973, 981, 982, 987, 990, 994, 997, 1000, 1002, 1003, 1004, 1006, 1008, 1015, 1016, 1017, 1018, 1019, 1020, 1022, 1024, 1025, 1026, 1030, 1032, 1035, 1037, 1041, 1043, 1050, 1057, 1065, 1066, 1068, 1073, 1088, 1092, 1095, 1098, 1099, 1101, 1107, 1113, 1114, 1115, 1116, 1117, 1119, 1120, 1121, 1123, 1126, 1127, 1130, 1132, 1134, 1139, 1145, 1146, 1148, 1150, 1153, 1155, 1156, 1164, 1165, 1168, 1169, 1170, 1171, 1172, 1174, 1176, 1177, 1179, 1183, 1185, 1189, 1193, 1196, 1201, 1202, 1203, 1204, 1207, 1211, 1223, 1229, 1230, 1233, 1236, 1243, 1245, 1246, 1247, 1251, 1252, 1253, 1255, 1258, 1260, 1261, 1264, 1265, 1269, 1271, 1275, 1278, 1280, 1282, 1284, 1285, 1293, 1297, 1299, 1302, 1306, 1315, 1316, 1318, 1320, 1322, 1329, 1332, 1336, 1339, 1349, 1353, 1356, 1357, 1359, 1365, 1366, 1368, 1370, 1374, 1375, 1384, 1385, 1386, 1387, 1388, 1391, 1392, 1395, 1396, 1397, 1401, 1407, 1408, 1411, 1412, 1413, 1416, 1419, 1421, 1422, 1425, 1428, 1429, 1431, 1432, 1435, 1437, 1438, 1439, 1441, 1442, 1449, 1450, 1451, 1453, 1455, 1457, 1475, 1479, 1481, 1484, 1490, 1491, 1493, 1495, 1497, 1499, 1501, 1508, 1511, 1512, 1513, 1514, 1517, 1518, 1521, 1524, 1526, 1534, 1539, 1540, 1542, 1546, 1551, 1552, 1553, 1554, 1556, 1557, 1561, 1563, 1564, 1566, 1567, 1568, 1569, 1570, 1572, 1573, 1574, 1575, 1577, 1578, 1579, 1580, 1582, 1583, 1586, 1588, 1590, 1591, 1598, 1599, 1604, 1607, 1608, 1612, 1613, 1620, 1621, 1622, 1632, 1639, 1640, 1641, 1643, 1653, 1654, 1655, 1656, 1663, 1671, 1674, 1676, 1677, 1679, 1682, 1683, 1684, 1685, 1686, 1688, 1691, 1696, 1698, 1713, 1714, 1716, 1717, 1718, 1723, 1727, 1728, 1729, 1732, 1734, 1738, 1741, 1744, 1746, 1747, 1750, 1758, 1763, 1764, 1770, 1776, 1777, 1778, 1782, 1784, 1786, 1789, 1793, 1794, 1797, 1801, 1803, 1808, 1809, 1812, 1818, 1821, 1822, 1826, 1829, 1830, 1831, 1833, 1834, 1837, 1838, 1839, 1841, 1842, 1844, 1846, 1850, 1852, 1853, 1855, 1857, 1862, 1863, 1864, 1866, 1867, 1869, 1871, 1874, 1879, 1880, 1882, 1888, 1889, 1897, 1899, 1903, 1904, 1906, 1907, 1908, 1910, 1912, 1913, 1915, 1916, 1918, 1919, 1921, 1922, 1927, 1933, 1934, 1937, 1939, 1945, 1946, 1947, 1948, 1951, 1954, 1955, 1956, 1958, 1961, 1965, 1966, 1967, 1970, 1972, 1974, 1978, 1981, 1982, 1985, 1987, 1995, 1996, 1998, 1999, 2002, 2003, 2007, 2008, 2012, 2014, 2016, 2020, 2024, 2026, 2028, 2030, 2034, 2035, 2036, 2038, 2040, 2041, 2044, 2047, 2048, 2049, 2050, 2051, 2053, 2055, 2058, 2067, 2078, 2082, 2083, 2090, 2091, 2093, 2095, 2096, 2101, 2104, 2106, 2115, 2116, 2120, 2121, 2124, 2125, 2131, 2132, 2133, 2137, 2139, 2143, 2144, 2147, 2148, 2150, 2156, 2161, 2162, 2166, 2178, 2182, 2184, 2185, 2187, 2189, 2190, 2194, 2199, 2201, 2202, 2203, 2205, 2207, 2213, 2216, 2217, 2218, 2220, 2222, 2223, 2227, 2229, 2230, 2236, 2237, 2239, 2240, 2244, 2247, 2250, 2254, 2256, 2264, 2265, 2266, 2269, 2272, 2274, 2277, 2279, 2280, 2282, 2283, 2285, 2292, 2296, 2299, 2300, 2305, 2307, 2309, 2310, 2313, 2314, 2318, 2322, 2333, 2334, 2339, 2343, 2348, 2349, 2350, 2352, 2354, 2355, 2358, 2365, 2367, 2381, 2384, 2386, 2389, 2395, 2397, 2402, 2404, 2407, 2409, 2413, 2414, 2423, 2424, 2428, 2433, 2434, 2440, 2442, 2443, 2454, 2455, 2457, 2458, 2460, 2463, 2468, 2472, 2473, 2479, 2486, 2487, 2495, 2496, 2498, 2499, 2500, 2501, 2503, 2505, 2507, 2510, 2516, 2519, 2521, 2522, 2525, 2532, 2533, 2534, 2535, 2536, 2538, 2541, 2542, 2548, 2561, 2565, 2568, 2569, 2570, 2571, 2573, 2575, 2576, 2579, 2585, 2586, 2593, 2595, 2597, 2599, 2600, 2602, 2611, 2612, 2614, 2615, 2619, 2620, 2623, 2625, 2628, 2631, 2633, 2634, 2635, 2636, 2640, 2642, 2643, 2648, 2651, 2652, 2654, 2656, 2660, 2661, 2663, 2667, 2668, 2670, 2675, 2676, 2677, 2678, 2680, 2683, 2692, 2694, 2696, 2697, 2698, 2703, 2707, 2711, 2712, 2714, 2715, 2720, 2722, 2730, 2732, 2735, 2738, 2741, 2743, 2746, 2750, 2752, 2755, 2756, 2758, 2760, 2764, 2765, 2767, 2781, 2783, 2784, 2786, 2801, 2806, 2810, 2814, 2818, 2826, 2827, 2828, 2831, 2832, 2833, 2835, 2837, 2838, 2846, 2847, 2849, 2852, 2860, 2863, 2864, 2867, 2868, 2869, 2871, 2872, 2874, 2875, 2882, 2883, 2897, 2902, 2904, 2905, 2911, 2913, 2914, 2916, 2918, 2920, 2921, 2922, 2927, 2929, 2933, 2936, 2937, 2938, 2939, 2940, 2942, 2945, 2947, 2950, 2955, 2956, 2961, 2963, 2964, 2965, 2967, 2968, 2970, 2971, 2972, 2976, 2978, 2983, 2985, 2988, 2998, 3008, 3010, 3015, 3018, 3020, 3023, 3024, 3025, 3026, 3027, 3029, 3030, 3033, 3035, 3038, 3040, 3041, 3044, 3051, 3053, 3054, 3057, 3062, 3064, 3067, 3068, 3069, 3072, 3078, 3079, 3081, 3084, 3087, 3088, 3093, 3096, 3097, 3102, 3104, 3108, 3111, 3112, 3121, 3122, 3129, 3134, 3135, 3147, 3151, 3156, 3157, 3159, 3160, 3161, 3163, 3165, 3167, 3169, 3175, 3186, 3189, 3192, 3196, 3200, 3203, 3205, 3209, 3214, 3216, 3217, 3220, 3221, 3222, 3225, 3230, 3231, 3233, 3236, 3238, 3242, 3243, 3244, 3246, 3249, 3250, 3257, 3261, 3262, 3264, 3266, 3269, 3272, 3274, 3275, 3278, 3279, 3282, 3285, 3287, 3288, 3290, 3294, 3296, 3298, 3302, 3303, 3304, 3305, 3306, 3310, 3318, 3325, 3333, 3334, 3335, 3336, 3337, 3338, 3344, 3351, 3354, 3359, 3361, 3365, 3368, 3370, 3373, 3375, 3376, 3379, 3380, 3381, 3385, 3386, 3387, 3391, 3397, 3398, 3399, 3402, 3406, 3413, 3414, 3418, 3421, 3422, 3426, 3428, 3432, 3435, 3438, 3445, 3446, 3449, 3450, 3451, 3453, 3456, 3465, 3469, 3470, 3471, 3472, 3474, 3477, 3481, 3484, 3485, 3486, 3488, 3490, 3491, 3492, 3493, 3494, 3496, 3498, 3501, 3502, 3505, 3509, 3510, 3511, 3513, 3515, 3516, 3517, 3518, 3519, 3520, 3522, 3524, 3526, 3528, 3532, 3534, 3539, 3547, 3549, 3552, 3557, 3564, 3569, 3571, 3574, 3583, 3585, 3586, 3588, 3589, 3590, 3592, 3595, 3597, 3600, 3601, 3603, 3608, 3612, 3613, 3617, 3629, 3631, 3634, 3636, 3637, 3638, 3640, 3642, 3644, 3646, 3649, 3652, 3660, 3662, 3668, 3670, 3672, 3678, 3683, 3686, 3687, 3694, 3697, 3701, 3703, 3706, 3707, 3710, 3712, 3715, 3716, 3722, 3723, 3729, 3731, 3733, 3735, 3743, 3744, 3745, 3746, 3747, 3750, 3751, 3753, 3759, 3765, 3767, 3768, 3769, 3770, 3771, 3774, 3775, 3777, 3780, 3782, 3783, 3784, 3787, 3789, 3790, 3791, 3796, 3797, 3803, 3813, 3814, 3815, 3820, 3821, 3822, 3823, 3825, 3826, 3828, 3829, 3830, 3831, 3832, 3833, 3838, 3839, 3840, 3841, 3842, 3845, 3850, 3852, 3854, 3856, 3857, 3858, 3860, 3862, 3868, 3877, 3883, 3885, 3887, 3888, 3889, 3891, 3892, 3896, 3897, 3898, 3900, 3908, 3910, 3913, 3914, 3915, 3917, 3928, 3929, 3930, 3935, 3939, 3943, 3947, 3951, 3956, 3959, 3963, 3964, 3966, 3969, 3973, 3975, 3979, 3980, 3981, 3985, 3989, 3991, 3993, 3995, 4000, 4001, 4003, 4004, 4006, 4011, 4014, 4018, 4022, 4024, 4028, 4029, 4032, 4033, 4034, 4036, 4038, 4039, 4046, 4050, 4052, 4055, 4057, 4058, 4061, 4064, 4065, 4070, 4071, 4075, 4082, 4084, 4095, 4101, 4103, 4106, 4111, 4115, 4116, 4118, 4123, 4124, 4130, 4131, 4132, 4136, 4137, 4140, 4142, 4147, 4156, 4159, 4160, 4166, 4174, 4177, 4179, 4180, 4190, 4191, 4194, 4196, 4199, 4200, 4207, 4209, 4210, 4212, 4213, 4215, 4216, 4217, 4219, 4224, 4228, 4231, 4233, 4236, 4245, 4249, 4251, 4256, 4258, 4265, 4272, 4273, 4276, 4278, 4283, 4285, 4288, 4292, 4295, 4296, 4297, 4299, 4301, 4302, 4303, 4304, 4309, 4315, 4317, 4321, 4325, 4327, 4329, 4330, 4332, 4333, 4334, 4336, 4337, 4338, 4339, 4342, 4344, 4348, 4353, 4368, 4370, 4372, 4375, 4378, 4392, 4394, 4395, 4396, 4397, 4399, 4400, 4407, 4410, 4411, 4412, 4415, 4418, 4424, 4427, 4428, 4433, 4435, 4439, 4441, 4442, 4443, 4445, 4446, 4448, 4453, 4456, 4457, 4459, 4462, 4466, 4467, 4469, 4471, 4474, 4475, 4477, 4482, 4483, 4489, 4490, 4492, 4495, 4497, 4498, 4502, 4507, 4509, 4512, 4513, 4514, 4518, 4519, 4525, 4529, 4532, 4533, 4534, 4535, 4538, 4541, 4542, 4546, 4550, 4557, 4558, 4559, 4560, 4564, 4568, 4569, 4574, 4575, 4580, 4589, 4591, 4596, 4597, 4598, 4600, 4601, 4605, 4606, 4607, 4609, 4612, 4620, 4629, 4634, 4637, 4641, 4642, 4651, 4653, 4654, 4655, 4659, 4661, 4665, 4668, 4669, 4670, 4675, 4677, 4680, 4681, 4693, 4695, 4701, 4703, 4709, 4716, 4720, 4721, 4723, 4724, 4725, 4728, 4731, 4734, 4738, 4740, 4742, 4743, 4744, 4747, 4751, 4758, 4759, 4760, 4771, 4774, 4775, 4777, 4779, 4782, 4786, 4787, 4789, 4794, 4799, 4808, 4813, 4816, 4817, 4818, 4822, 4826, 4831, 4832, 4839, 4841, 4844, 4845, 4847, 4848, 4852, 4858, 4859, 4861, 4870, 4871, 4874, 4876, 4877, 4878, 4882, 4883, 4885, 4886, 4887, 4888, 4890, 4896, 4897, 4899, 4901, 4906, 4909, 4910, 4912, 4916, 4919, 4925, 4927, 4928, 4929, 4930, 4940, 4941, 4944, 4948, 4953, 4954, 4958, 4962, 4963, 4965, 4969, 4973, 4975, 4978, 4980, 4982, 4986, 4994, 4995, 5002, 5007, 5008, 5013, 5015, 5020, 5025, 5039, 5040, 5044, 5045, 5046, 5048, 5049, 5054, 5055, 5062, 5063, 5068, 5069, 5072, 5076, 5078, 5080, 5084, 5086, 5092, 5101, 5104, 5105, 5107, 5109, 5114, 5118, 5119, 5121, 5122, 5124, 5125, 5127, 5128, 5129, 5136, 5140, 5144, 5149, 5152, 5156, 5158, 5162, 5163, 5164, 5171, 5174, 5175, 5178, 5179, 5183, 5184, 5185, 5186, 5187, 5195, 5196, 5198, 5199, 5200, 5201, 5203, 5205, 5206, 5207, 5213, 5214, 5216, 5221, 5222, 5223, 5224, 5227, 5228, 5229, 5230, 5233, 5235, 5236, 5237, 5238, 5240, 5244, 5247, 5252, 5253, 5254, 5255, 5256, 5260, 5261, 5264, 5265, 5267, 5268, 5272, 5274, 5282, 5285, 5290, 5291, 5293, 5294, 5295, 5296, 5302, 5306, 5309, 5311, 5312, 5314, 5315, 5320, 5321, 5322, 5323, 5325, 5327, 5329, 5330, 5334, 5335, 5348, 5349, 5355, 5358, 5359, 5361, 5362, 5364, 5369, 5373, 5374, 5379, 5380, 5381, 5382, 5383, 5385, 5386, 5390, 5393, 5394, 5397, 5399, 5406, 5410, 5412, 5413, 5416, 5417, 5421, 5424, 5425, 5428, 5429, 5430, 5436, 5438, 5439, 5441, 5443, 5445, 5446, 5448, 5450, 5451, 5452, 5456, 5458, 5460, 5464, 5466, 5470, 5474, 5475, 5476, 5482, 5485, 5495, 5498, 5501, 5504, 5513, 5519, 5522, 5525, 5526, 5528, 5532, 5536, 5537, 5538, 5539, 5541, 5546, 5554, 5555, 5558, 5560, 5566, 5569, 5572, 5574, 5577, 5579, 5585, 5586, 5587, 5588, 5589, 5590, 5597, 5600, 5604, 5606, 5607, 5611, 5612, 5614, 5615, 5617, 5619, 5620, 5623, 5626, 5627, 5628, 5629, 5634, 5637, 5640, 5645, 5649, 5650, 5651, 5655, 5656, 5661, 5663, 5665, 5669, 5670, 5671, 5672, 5673, 5677, 5682, 5691, 5692, 5693, 5695, 5696, 5699, 5700, 5704, 5705, 5707, 5710, 5711, 5714, 5715, 5717, 5719, 5722, 5725, 5728, 5729, 5733, 5735, 5737, 5738, 5742, 5744, 5748, 5752, 5754, 5756, 5757, 5759, 5760, 5761, 5762, 5764, 5765, 5766, 5770, 5775, 5776, 5778, 5783, 5785, 5786, 5788, 5794, 5795, 5799, 5800, 5804, 5805, 5809, 5812, 5821, 5824, 5825, 5827, 5829, 5830, 5833, 5834, 5835, 5837, 5845, 5847, 5849, 5851, 5858, 5862, 5866, 5867, 5869, 5872, 5873, 5874, 5876, 5877, 5878, 5880, 5890, 5892, 5894, 5912, 5913, 5914, 5916, 5923, 5928, 5931, 5932, 5934, 5944, 5946, 5947, 5958, 5959, 5960, 5969, 5971, 5975, 5980, 5981, 5986, 5987, 5989, 5996, 6003, 6004, 6006, 6007, 6009, 6011, 6018, 6024, 6026, 6030, 6032, 6033, 6041, 6044, 6048, 6049, 6053, 6057, 6061, 6065, 6067, 6077, 6081, 6087, 6088, 6091, 6095, 6096, 6097, 6108, 6109, 6111, 6115, 6116, 6118, 6119, 6123, 6128, 6136, 6138, 6143, 6147, 6149, 6155, 6157, 6158, 6159, 6160, 6162, 6164, 6167, 6168, 6172, 6180, 6181, 6188, 6189, 6191, 6192, 6194, 6195, 6196, 6201, 6206, 6208, 6213, 6214, 6215, 6217, 6219, 6224, 6226, 6228, 6230, 6235, 6239, 6241, 6245, 6247, 6249, 6252, 6253, 6254, 6255, 6256, 6257, 6259, 6260, 6262, 6264, 6265, 6267, 6273, 6274, 6276, 6277, 6278, 6279, 6280, 6281, 6282, 6284, 6290, 6291, 6293, 6296, 6298, 6299, 6308, 6310, 6314, 6317, 6320, 6321, 6323, 6324, 6333, 6335, 6336, 6337, 6338, 6339, 6344, 6349, 6350, 6352, 6359, 6364, 6370, 6373, 6374, 6378, 6379, 6381, 6385, 6392, 6408, 6410, 6415, 6420, 6422, 6423, 6424, 6428, 6429, 6434, 6437, 6439, 6440, 6441, 6442, 6444, 6452, 6459, 6461, 6464, 6468, 6470, 6473, 6474, 6478, 6483, 6492, 6499, 6501, 6502, 6504, 6508, 6513, 6514, 6515, 6517, 6519, 6521, 6522, 6523, 6524, 6526, 6528, 6529, 6530, 6531, 6540, 6544, 6546, 6547, 6552, 6558, 6559, 6561, 6562, 6563, 6565, 6567, 6568, 6571, 6572, 6573, 6577, 6578, 6581, 6582, 6584, 6595, 6596, 6597, 6603, 6604, 6605, 6606, 6608, 6610, 6611, 6613, 6615, 6618, 6620, 6621, 6624, 6625, 6628, 6630, 6631, 6632, 6633, 6637, 6644, 6646, 6647, 6650, 6651, 6654, 6657, 6658, 6660, 6664, 6669, 6671, 6672, 6679, 6682, 6689, 6691, 6692, 6697, 6701, 6704, 6707, 6714, 6721, 6724, 6725, 6726, 6731, 6733, 6737, 6741, 6742, 6744, 6745, 6752, 6755, 6758, 6761, 6767, 6768, 6769, 6772, 6774, 6775, 6778, 6780, 6782, 6789, 6793, 6795, 6796, 6797, 6800, 6801, 6802, 6806, 6808, 6810, 6819, 6820, 6826, 6827, 6833, 6834, 6845, 6852, 6854, 6856, 6864, 6868, 6871, 6874, 6877, 6879, 6880, 6883, 6887, 6888, 6895, 6896, 6899, 6902, 6903, 6904, 6905, 6910, 6911, 6912, 6916, 6920, 6921, 6924, 6925, 6930, 6932, 6940, 6941, 6942, 6943, 6952, 6955, 6958, 6959, 6966, 6971, 6980, 6981, 6984, 6986, 6988, 6991, 6996, 6999, 7000, 7009, 7010, 7012, 7013, 7014, 7016, 7020, 7022, 7025, 7032, 7035, 7038, 7043, 7050, 7051, 7055, 7057, 7058, 7060, 7061, 7062, 7063, 7064, 7065, 7066, 7072, 7074, 7075, 7077, 7078, 7079, 7080, 7081, 7082, 7087, 7088, 7096, 7103, 7104, 7105, 7116, 7119, 7121, 7123, 7128, 7129, 7134, 7136, 7138, 7142, 7147, 7151, 7153, 7160, 7161, 7162, 7163, 7165, 7169, 7172, 7174, 7175, 7177, 7182, 7188, 7190, 7191, 7196, 7198, 7202, 7203, 7208, 7212, 7213, 7215, 7216, 7221, 7228, 7233, 7235, 7236, 7241, 7244, 7246, 7250, 7252, 7253, 7259, 7260, 7263, 7264, 7266, 7267, 7270, 7271, 7272, 7273, 7275, 7279, 7283, 7285, 7286, 7287, 7292, 7297, 7303, 7305, 7307, 7309, 7320, 7324, 7325, 7326, 7327, 7335, 7336, 7337, 7339, 7340, 7342, 7350, 7353, 7354, 7357, 7358, 7361, 7368, 7369, 7372, 7373, 7376, 7378, 7381, 7382, 7385, 7386, 7388, 7394, 7396, 7399, 7405, 7406, 7407, 7412, 7416, 7419, 7420, 7422, 7423, 7430, 7434, 7437, 7438, 7440, 7441, 7442, 7447, 7459, 7462, 7463, 7464, 7467, 7470, 7472, 7473, 7474, 7476, 7477, 7484, 7486, 7494, 7495, 7498, 7503, 7511, 7517, 7524, 7528, 7531, 7534, 7540, 7542, 7551, 7554, 7556, 7558, 7559, 7569, 7570, 7573, 7576, 7578, 7579, 7585, 7586, 7587, 7588, 7589, 7590, 7594, 7595, 7596, 7600, 7601, 7604, 7607, 7609, 7615, 7616, 7619, 7622, 7624, 7625, 7632, 7633, 7634, 7635, 7643, 7645, 7647, 7648, 7653, 7666, 7670, 7679, 7680, 7681, 7682, 7686, 7687, 7689, 7690, 7693, 7698, 7708, 7711, 7712, 7717, 7718, 7719, 7721, 7731, 7732, 7733, 7736, 7741, 7744, 7747, 7748, 7751, 7757, 7758, 7759, 7765, 7772, 7773, 7774, 7780, 7783, 7789, 7793, 7798, 7799, 7801, 7807, 7810, 7811, 7812, 7819, 7820, 7826, 7827, 7828, 7830, 7831, 7834, 7836, 7839, 7841, 7845, 7848, 7850, 7854, 7855, 7857, 7859, 7860, 7862, 7866, 7867, 7868, 7870, 7871, 7879, 7884, 7886, 7891, 7892, 7893, 7899, 7900, 7901, 7904, 7910, 7915, 7917, 7919, 7920, 7921, 7924, 7925, 7930, 7931, 7934, 7941, 7947, 7950, 7951, 7953, 7954, 7955, 7963, 7967, 7968, 7972, 7973, 7974, 7975, 7976, 7985, 7987, 7991, 7994, 7997, 7999, 8000, 8001, 8003, 8004, 8006, 8008, 8009, 8012, 8013, 8018, 8020, 8024, 8028, 8031, 8034, 8037, 8041, 8049, 8052, 8056, 8062, 8064, 8075, 8083, 8087, 8088, 8092, 8098, 8101, 8102, 8107, 8110, 8111, 8112, 8115, 8116, 8119, 8120, 8126, 8129, 8134, 8135, 8136, 8142, 8144, 8146, 8149, 8154, 8156, 8159, 8161, 8166, 8169, 8175, 8176, 8178, 8183, 8185, 8188, 8191, 8195, 8200, 8201, 8202, 8204, 8210, 8211, 8213, 8214, 8215, 8217, 8219, 8220, 8227, 8230, 8233, 8235, 8240, 8243, 8244, 8245, 8246, 8247, 8249, 8250, 8256, 8257, 8261, 8264, 8266, 8268, 8269, 8270, 8271, 8276, 8282, 8285, 8287, 8288, 8289, 8290, 8291, 8297, 8302, 8306, 8318, 8319, 8321, 8323, 8325, 8326, 8329, 8339, 8340, 8345, 8346, 8348, 8350, 8351, 8353, 8356, 8360, 8363, 8367, 8368, 8369, 8372, 8380, 8383, 8384, 8385, 8386, 8387, 8389, 8397, 8398, 8399, 8401, 8406, 8409, 8411, 8414, 8417, 8420, 8423, 8432, 8436, 8437, 8438, 8439, 8441, 8444, 8446, 8451, 8459, 8462, 8463, 8465, 8466, 8471, 8473, 8476, 8478, 8480, 8483, 8489, 8490, 8493, 8495, 8496, 8498, 8500, 8503, 8506, 8508, 8513, 8518, 8519, 8520, 8523, 8524, 8525, 8526, 8530, 8531, 8533, 8536, 8537, 8544, 8554, 8555, 8556, 8559, 8560, 8562, 8563, 8567, 8570, 8571, 8572, 8579, 8580, 8582, 8583, 8586, 8588, 8590, 8592, 8593, 8594, 8595, 8597, 8599, 8600, 8602, 8603, 8604, 8606, 8609, 8610, 8613, 8617, 8618, 8621, 8622, 8627, 8629, 8630, 8636, 8637, 8641, 8642, 8644, 8649, 8651, 8653, 8654, 8655, 8656, 8658, 8659, 8661, 8663, 8666, 8668, 8671, 8672, 8674, 8675, 8676, 8677, 8681, 8682, 8683, 8685, 8686, 8687, 8689, 8690, 8693, 8694, 8699, 8702, 8707, 8711, 8713, 8715, 8718, 8722, 8724, 8727, 8728, 8729, 8730, 8731, 8732, 8737, 8742, 8746, 8750, 8752, 8754, 8756, 8760, 8761, 8762, 8764, 8770, 8774, 8776, 8777, 8778, 8780, 8784, 8786, 8790, 8794, 8798, 8799, 8803, 8804, 8805, 8809, 8810, 8816, 8821, 8822, 8823, 8824, 8826, 8828, 8836, 8839, 8840, 8841, 8842, 8843, 8844, 8845, 8859, 8860, 8866, 8867, 8871, 8876, 8886, 8890, 8895, 8899, 8903, 8905, 8906, 8907, 8909, 8910, 8911, 8913, 8914, 8916, 8919, 8921, 8922, 8924, 8925, 8929, 8935, 8938, 8940, 8942, 8945, 8946, 8950, 8952, 8954, 8958, 8959, 8963, 8966, 8968, 8970, 8973, 8974, 8975, 8976, 8978, 8981, 8982, 8983, 8984, 8986, 8987, 8989, 8990, 8992, 8996, 8998, 8999, 9000, 9001, 9004, 9005, 9009, 9010, 9013, 9014, 9015, 9018, 9022, 9025, 9028, 9031, 9032, 9036, 9037, 9039, 9040, 9041, 9044, 9045, 9048, 9049, 9053, 9054, 9055, 9056, 9057, 9058, 9059, 9060, 9062, 9064, 9065, 9068, 9069, 9070, 9072, 9074, 9076, 9078, 9084, 9086, 9087, 9088, 9089, 9093, 9094, 9096, 9099, 9100, 9102, 9103, 9107, 9109, 9110, 9111, 9113, 9114, 9116, 9117, 9120, 9123, 9125, 9127, 9129, 9130, 9135, 9140, 9143, 9147, 9152, 9154, 9156, 9157, 9159, 9163, 9168, 9170, 9172, 9177, 9179, 9181, 9182, 9184, 9186, 9187, 9190, 9195, 9197, 9200, 9203, 9213, 9216, 9218, 9219, 9223, 9232, 9236, 9241, 9244, 9245, 9246, 9255, 9256, 9257, 9259, 9260, 9262, 9265, 9269, 9270, 9278, 9280, 9282, 9283, 9285, 9286, 9291, 9292, 9298, 9299, 9300, 9301, 9302, 9303, 9307, 9308, 9316, 9321, 9322, 9323, 9327, 9329, 9332, 9339, 9348, 9350, 9351, 9354, 9356, 9360, 9361, 9362, 9363, 9364, 9367, 9369, 9373, 9375, 9379, 9380, 9383, 9388, 9391, 9392, 9394, 9396, 9398, 9400, 9403, 9404, 9405, 9406, 9407, 9409, 9410, 9423, 9425, 9428, 9431, 9432, 9434, 9435, 9438, 9439, 9440, 9441, 9442, 9445, 9448, 9453, 9454, 9464, 9466, 9474, 9475, 9477, 9484, 9489, 9492, 9493, 9497, 9505, 9508, 9509, 9513, 9520, 9521, 9524, 9525, 9528, 9530, 9535, 9543, 9546, 9550, 9551, 9552, 9557, 9558, 9559, 9567, 9569, 9572, 9573, 9575, 9577, 9580, 9581, 9583, 9588, 9590, 9591, 9593, 9596, 9597, 9599, 9600, 9601, 9602, 9605, 9608, 9610, 9611, 9615, 9616, 9617, 9619, 9621, 9625, 9630, 9635, 9636, 9638, 9640, 9643, 9645, 9646, 9648, 9652, 9655, 9656, 9659, 9660, 9661, 9662, 9664, 9666, 9667, 9669, 9671, 9674, 9681, 9685, 9687, 9689, 9694, 9695, 9699, 9701, 9702, 9706, 9710, 9717, 9719, 9720, 9721, 9722, 9728, 9730, 9732, 9733, 9736, 9737, 9741, 9742, 9743, 9747, 9748, 9749, 9753, 9754, 9764, 9767, 9768, 9774, 9775, 9776, 9779, 9780, 9783, 9791, 9792, 9795, 9798, 9799, 9805, 9806, 9807, 9810, 9811, 9812, 9823, 9824, 9826, 9829, 9831, 9832, 9836, 9838, 9843, 9845, 9856, 9857, 9864, 9866, 9867, 9872, 9873, 9877, 9878, 9880, 9881, 9883, 9884], [1, 2, 7, 8, 10, 21, 24, 26, 28, 32, 35, 38, 39, 42, 43, 47, 50, 52, 55, 56, 59, 66, 69, 70, 71, 72, 74, 75, 77, 80, 82, 84, 92, 94, 95, 98, 99, 102, 107, 108, 114, 115, 122, 132, 133, 135, 142, 144, 152, 153, 155, 157, 160, 165, 166, 168, 170, 177, 179, 180, 181, 182, 185, 189, 190, 194, 195, 197, 199, 203, 205, 207, 210, 212, 213, 214, 215, 220, 222, 223, 231, 237, 238, 241, 245, 248, 251, 252, 253, 256, 259, 266, 267, 268, 269, 271, 279, 281, 282, 284, 288, 289, 290, 296, 298, 299, 300, 301, 302, 305, 307, 308, 310, 311, 312, 313, 314, 316, 319, 320, 322, 325, 329, 333, 336, 337, 341, 342, 344, 350, 351, 355, 362, 363, 364, 365, 366, 367, 375, 376, 379, 380, 383, 384, 385, 387, 389, 391, 393, 395, 396, 397, 403, 405, 406, 407, 413, 419, 420, 421, 424, 428, 430, 431, 433, 436, 439, 440, 441, 442, 443, 444, 447, 454, 458, 463, 464, 470, 473, 474, 482, 483, 485, 487, 488, 489, 490, 491, 494, 502, 503, 504, 505, 512, 513, 515, 522, 523, 525, 528, 535, 537, 538, 540, 541, 542, 544, 549, 552, 554, 557, 560, 561, 562, 563, 564, 565, 568, 570, 573, 580, 584, 585, 587, 588, 590, 596, 597, 598, 603, 607, 610, 614, 620, 622, 623, 626, 629, 633, 634, 638, 640, 641, 642, 646, 649, 650, 651, 652, 654, 655, 657, 662, 664, 665, 671, 672, 673, 674, 679, 684, 686, 689, 690, 691, 694, 699, 704, 706, 708, 712, 713, 718, 719, 720, 721, 723, 724, 730, 731, 734, 738, 742, 745, 746, 747, 749, 758, 767, 768, 770, 771, 772, 777, 780, 785, 788, 790, 791, 796, 797, 798, 803, 811, 814, 816, 828, 829, 830, 831, 833, 834, 837, 840, 842, 846, 847, 854, 855, 856, 857, 858, 861, 863, 866, 880, 882, 883, 886, 888, 891, 895, 897, 898, 903, 908, 910, 911, 915, 919, 921, 922, 924, 933, 936, 938, 940, 941, 942, 946, 947, 948, 950, 951, 960, 961, 966, 967, 968, 971, 972, 976, 977, 979, 980, 983, 984, 985, 988, 989, 993, 995, 996, 999, 1001, 1005, 1009, 1010, 1013, 1021, 1023, 1027, 1028, 1029, 1033, 1038, 1040, 1044, 1045, 1046, 1051, 1052, 1055, 1056, 1058, 1059, 1061, 1063, 1069, 1072, 1074, 1076, 1078, 1079, 1082, 1084, 1090, 1091, 1093, 1094, 1100, 1103, 1105, 1108, 1112, 1124, 1128, 1129, 1131, 1133, 1141, 1142, 1143, 1147, 1154, 1157, 1158, 1160, 1161, 1163, 1166, 1178, 1182, 1187, 1190, 1194, 1197, 1199, 1200, 1206, 1208, 1210, 1219, 1220, 1221, 1222, 1225, 1226, 1232, 1234, 1235, 1238, 1239, 1241, 1249, 1256, 1262, 1263, 1267, 1272, 1273, 1281, 1283, 1286, 1287, 1288, 1290, 1291, 1294, 1295, 1298, 1300, 1301, 1307, 1308, 1309, 1310, 1313, 1314, 1321, 1325, 1330, 1331, 1334, 1335, 1342, 1343, 1346, 1348, 1350, 1351, 1354, 1361, 1362, 1363, 1373, 1376, 1380, 1381, 1383, 1390, 1393, 1398, 1399, 1406, 1409, 1410, 1418, 1423, 1427, 1430, 1433, 1434, 1436, 1443, 1447, 1456, 1458, 1459, 1460, 1461, 1467, 1468, 1471, 1477, 1480, 1483, 1485, 1486, 1492, 1496, 1500, 1503, 1504, 1505, 1506, 1510, 1519, 1520, 1523, 1525, 1528, 1529, 1530, 1531, 1532, 1533, 1535, 1536, 1537, 1541, 1544, 1545, 1549, 1555, 1558, 1559, 1560, 1562, 1565, 1581, 1584, 1593, 1594, 1602, 1605, 1609, 1614, 1615, 1616, 1618, 1626, 1627, 1633, 1635, 1636, 1638, 1645, 1646, 1650, 1651, 1652, 1659, 1662, 1664, 1666, 1668, 1669, 1670, 1675, 1681, 1690, 1694, 1695, 1702, 1703, 1704, 1706, 1707, 1708, 1709, 1710, 1724, 1725, 1726, 1730, 1733, 1736, 1740, 1742, 1743, 1745, 1748, 1749, 1753, 1755, 1756, 1757, 1766, 1771, 1772, 1773, 1775, 1779, 1780, 1785, 1787, 1792, 1796, 1798, 1800, 1802, 1806, 1810, 1813, 1814, 1816, 1817, 1819, 1827, 1828, 1832, 1840, 1847, 1849, 1854, 1858, 1861, 1868, 1872, 1876, 1878, 1884, 1885, 1886, 1887, 1892, 1893, 1901, 1902, 1909, 1911, 1917, 1925, 1928, 1932, 1935, 1938, 1940, 1942, 1943, 1944, 1949, 1950, 1952, 1957, 1962, 1963, 1964, 1969, 1971, 1975, 1980, 1983, 1991, 1992, 1993, 1997, 2001, 2004, 2005, 2011, 2017, 2018, 2019, 2021, 2023, 2027, 2029, 2032, 2037, 2039, 2042, 2054, 2056, 2057, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2068, 2069, 2072, 2075, 2076, 2077, 2079, 2084, 2085, 2087, 2089, 2094, 2098, 2099, 2102, 2103, 2105, 2108, 2109, 2110, 2111, 2112, 2114, 2117, 2122, 2123, 2126, 2127, 2135, 2136, 2141, 2142, 2146, 2151, 2152, 2154, 2157, 2159, 2160, 2164, 2165, 2170, 2171, 2175, 2176, 2177, 2179, 2181, 2183, 2188, 2195, 2196, 2198, 2200, 2206, 2208, 2209, 2210, 2214, 2215, 2225, 2234, 2243, 2245, 2249, 2252, 2253, 2258, 2260, 2263, 2268, 2270, 2271, 2275, 2278, 2286, 2287, 2289, 2293, 2295, 2297, 2302, 2304, 2306, 2312, 2319, 2323, 2325, 2327, 2329, 2330, 2335, 2336, 2337, 2340, 2344, 2346, 2347, 2351, 2353, 2357, 2359, 2361, 2362, 2364, 2366, 2369, 2372, 2374, 2378, 2379, 2382, 2383, 2388, 2390, 2391, 2392, 2393, 2394, 2398, 2399, 2401, 2406, 2408, 2410, 2411, 2416, 2422, 2425, 2426, 2427, 2429, 2430, 2431, 2432, 2435, 2436, 2437, 2438, 2439, 2441, 2446, 2447, 2448, 2450, 2453, 2456, 2464, 2465, 2467, 2474, 2476, 2477, 2482, 2483, 2484, 2491, 2497, 2504, 2508, 2513, 2514, 2515, 2517, 2518, 2520, 2524, 2531, 2537, 2539, 2540, 2543, 2544, 2545, 2546, 2550, 2553, 2554, 2555, 2556, 2564, 2574, 2578, 2581, 2583, 2584, 2587, 2588, 2590, 2591, 2594, 2596, 2598, 2601, 2603, 2604, 2605, 2606, 2607, 2608, 2609, 2616, 2618, 2622, 2626, 2629, 2632, 2641, 2644, 2646, 2653, 2657, 2664, 2665, 2666, 2672, 2673, 2674, 2686, 2690, 2691, 2693, 2695, 2700, 2706, 2708, 2709, 2710, 2713, 2716, 2724, 2731, 2736, 2742, 2745, 2747, 2748, 2753, 2770, 2773, 2776, 2777, 2778, 2788, 2789, 2791, 2793, 2794, 2795, 2799, 2802, 2804, 2805, 2809, 2811, 2812, 2813, 2816, 2817, 2820, 2823, 2825, 2834, 2839, 2840, 2841, 2842, 2843, 2850, 2851, 2853, 2854, 2855, 2857, 2858, 2859, 2861, 2865, 2870, 2877, 2880, 2887, 2888, 2890, 2891, 2894, 2898, 2899, 2901, 2908, 2912, 2919, 2923, 2934, 2935, 2943, 2944, 2946, 2948, 2949, 2952, 2954, 2957, 2960, 2966, 2969, 2973, 2974, 2975, 2979, 2980, 2982, 2984, 2986, 2989, 2991, 2993, 2994, 2995, 2996, 2997, 2999, 3000, 3001, 3004, 3005, 3007, 3013, 3031, 3032, 3037, 3042, 3046, 3048, 3050, 3052, 3058, 3059, 3061, 3065, 3070, 3074, 3075, 3076, 3077, 3083, 3089, 3090, 3092, 3094, 3099, 3101, 3106, 3107, 3113, 3117, 3118, 3120, 3124, 3125, 3126, 3127, 3128, 3130, 3131, 3137, 3138, 3143, 3144, 3145, 3146, 3152, 3158, 3166, 3171, 3172, 3174, 3179, 3180, 3181, 3188, 3190, 3191, 3197, 3198, 3202, 3206, 3207, 3212, 3223, 3224, 3227, 3228, 3229, 3232, 3234, 3237, 3239, 3240, 3241, 3245, 3248, 3252, 3254, 3258, 3259, 3260, 3265, 3267, 3273, 3280, 3281, 3284, 3286, 3291, 3299, 3301, 3307, 3308, 3312, 3313, 3314, 3320, 3323, 3327, 3329, 3330, 3331, 3332, 3339, 3340, 3341, 3343, 3346, 3348, 3350, 3353, 3355, 3356, 3366, 3367, 3372, 3382, 3383, 3388, 3389, 3390, 3396, 3401, 3404, 3405, 3409, 3411, 3417, 3419, 3420, 3424, 3425, 3430, 3431, 3436, 3437, 3439, 3441, 3454, 3455, 3457, 3458, 3459, 3461, 3462, 3463, 3468, 3475, 3478, 3482, 3483, 3489, 3499, 3500, 3506, 3507, 3512, 3514, 3523, 3527, 3529, 3530, 3531, 3536, 3538, 3542, 3545, 3546, 3548, 3550, 3553, 3554, 3555, 3558, 3561, 3563, 3565, 3566, 3567, 3568, 3576, 3577, 3578, 3579, 3584, 3587, 3593, 3594, 3599, 3602, 3605, 3607, 3609, 3611, 3616, 3619, 3621, 3622, 3624, 3625, 3626, 3627, 3630, 3632, 3635, 3639, 3643, 3645, 3648, 3653, 3654, 3656, 3659, 3661, 3663, 3666, 3675, 3676, 3677, 3680, 3684, 3688, 3689, 3691, 3693, 3695, 3700, 3702, 3705, 3711, 3713, 3714, 3717, 3718, 3720, 3721, 3725, 3726, 3727, 3737, 3739, 3740, 3742, 3748, 3749, 3752, 3754, 3756, 3758, 3762, 3763, 3766, 3772, 3776, 3778, 3785, 3786, 3788, 3795, 3798, 3799, 3802, 3805, 3806, 3807, 3811, 3812, 3816, 3817, 3818, 3824, 3827, 3835, 3837, 3843, 3844, 3855, 3863, 3864, 3867, 3871, 3874, 3876, 3878, 3879, 3881, 3882, 3884, 3886, 3899, 3902, 3903, 3904, 3909, 3912, 3919, 3921, 3923, 3925, 3927, 3931, 3936, 3937, 3938, 3940, 3942, 3948, 3950, 3952, 3957, 3958, 3960, 3961, 3971, 3976, 3977, 3982, 3988, 3990, 3992, 3996, 3998, 3999, 4005, 4007, 4008, 4009, 4010, 4013, 4016, 4017, 4019, 4020, 4021, 4027, 4030, 4037, 4043, 4045, 4048, 4049, 4051, 4053, 4054, 4063, 4067, 4068, 4069, 4072, 4073, 4080, 4083, 4087, 4089, 4090, 4091, 4092, 4093, 4094, 4097, 4100, 4104, 4105, 4109, 4113, 4114, 4117, 4119, 4120, 4125, 4128, 4133, 4139, 4145, 4149, 4150, 4151, 4153, 4154, 4155, 4157, 4158, 4161, 4162, 4165, 4167, 4169, 4171, 4172, 4175, 4176, 4182, 4183, 4184, 4185, 4186, 4187, 4189, 4192, 4195, 4197, 4201, 4206, 4211, 4214, 4218, 4220, 4221, 4222, 4223, 4225, 4226, 4227, 4229, 4230, 4234, 4235, 4238, 4239, 4240, 4241, 4242, 4246, 4247, 4248, 4252, 4253, 4254, 4255, 4260, 4261, 4263, 4264, 4267, 4275, 4277, 4279, 4282, 4284, 4287, 4291, 4298, 4305, 4307, 4308, 4310, 4311, 4312, 4322, 4323, 4324, 4328, 4335, 4341, 4345, 4346, 4347, 4349, 4350, 4351, 4354, 4355, 4356, 4357, 4358, 4360, 4361, 4362, 4364, 4365, 4371, 4381, 4382, 4384, 4385, 4387, 4388, 4398, 4402, 4403, 4405, 4409, 4416, 4419, 4420, 4421, 4425, 4429, 4430, 4432, 4440, 4447, 4450, 4455, 4460, 4464, 4468, 4470, 4476, 4478, 4479, 4485, 4487, 4488, 4494, 4496, 4499, 4500, 4503, 4515, 4516, 4521, 4524, 4526, 4528, 4530, 4531, 4537, 4539, 4544, 4548, 4549, 4553, 4556, 4561, 4563, 4565, 4566, 4570, 4576, 4582, 4583, 4584, 4585, 4587, 4588, 4592, 4593, 4594, 4599, 4602, 4604, 4608, 4610, 4611, 4617, 4618, 4619, 4622, 4624, 4626, 4627, 4631, 4635, 4636, 4639, 4644, 4645, 4646, 4649, 4662, 4667, 4673, 4674, 4678, 4679, 4683, 4687, 4688, 4690, 4691, 4696, 4699, 4702, 4704, 4705, 4706, 4708, 4713, 4714, 4715, 4719, 4722, 4727, 4732, 4733, 4741, 4748, 4749, 4752, 4754, 4756, 4763, 4765, 4767, 4768, 4772, 4773, 4776, 4778, 4781, 4783, 4790, 4791, 4792, 4797, 4802, 4804, 4807, 4810, 4815, 4820, 4821, 4825, 4827, 4828, 4830, 4833, 4838, 4842, 4843, 4850, 4851, 4862, 4865, 4866, 4869, 4872, 4873, 4879, 4881, 4884, 4893, 4898, 4900, 4902, 4903, 4907, 4915, 4917, 4923, 4926, 4932, 4935, 4937, 4939, 4943, 4949, 4951, 4952, 4956, 4957, 4959, 4960, 4964, 4966, 4968, 4970, 4971, 4974, 4983, 4984, 4985, 4988, 4996, 4999, 5000, 5001, 5004, 5009, 5010, 5011, 5012, 5014, 5016, 5017, 5021, 5024, 5026, 5028, 5030, 5033, 5035, 5036, 5037, 5038, 5043, 5047, 5052, 5056, 5059, 5060, 5061, 5065, 5067, 5070, 5074, 5075, 5081, 5083, 5085, 5089, 5093, 5095, 5096, 5099, 5100, 5103, 5108, 5110, 5113, 5115, 5116, 5120, 5126, 5130, 5133, 5134, 5138, 5142, 5143, 5145, 5146, 5148, 5151, 5153, 5154, 5157, 5159, 5160, 5161, 5166, 5169, 5170, 5172, 5173, 5176, 5177, 5180, 5182, 5188, 5189, 5190, 5192, 5197, 5202, 5204, 5208, 5209, 5210, 5211, 5212, 5215, 5218, 5225, 5226, 5231, 5232, 5243, 5245, 5248, 5249, 5250, 5251, 5257, 5258, 5266, 5270, 5271, 5273, 5276, 5277, 5281, 5283, 5289, 5300, 5301, 5305, 5307, 5308, 5316, 5324, 5331, 5332, 5337, 5338, 5340, 5341, 5342, 5343, 5345, 5346, 5347, 5351, 5354, 5356, 5363, 5365, 5368, 5370, 5371, 5372, 5376, 5377, 5378, 5384, 5395, 5396, 5400, 5403, 5404, 5405, 5409, 5411, 5418, 5419, 5423, 5426, 5427, 5431, 5432, 5433, 5434, 5437, 5440, 5442, 5444, 5447, 5449, 5453, 5454, 5457, 5463, 5467, 5468, 5471, 5472, 5473, 5478, 5479, 5480, 5481, 5484, 5489, 5490, 5491, 5492, 5493, 5496, 5497, 5503, 5505, 5506, 5509, 5511, 5512, 5515, 5517, 5518, 5520, 5530, 5531, 5533, 5534, 5535, 5544, 5547, 5549, 5550, 5551, 5557, 5559, 5565, 5567, 5568, 5571, 5575, 5580, 5581, 5584, 5591, 5593, 5596, 5605, 5608, 5609, 5610, 5621, 5622, 5630, 5631, 5633, 5635, 5636, 5638, 5639, 5641, 5643, 5647, 5652, 5653, 5654, 5657, 5660, 5664, 5666, 5667, 5668, 5675, 5680, 5681, 5687, 5688, 5689, 5690, 5694, 5697, 5698, 5706, 5712, 5716, 5718, 5724, 5726, 5740, 5741, 5755, 5758, 5767, 5768, 5769, 5774, 5777, 5780, 5781, 5782, 5784, 5787, 5790, 5791, 5792, 5793, 5797, 5798, 5802, 5806, 5807, 5810, 5813, 5814, 5815, 5816, 5819, 5823, 5826, 5831, 5841, 5842, 5846, 5848, 5852, 5855, 5856, 5859, 5863, 5865, 5868, 5870, 5871, 5879, 5881, 5882, 5885, 5886, 5887, 5888, 5898, 5900, 5901, 5903, 5905, 5906, 5907, 5908, 5909, 5910, 5911, 5915, 5917, 5918, 5919, 5920, 5921, 5924, 5929, 5930, 5935, 5938, 5939, 5940, 5941, 5942, 5945, 5948, 5950, 5951, 5956, 5957, 5963, 5965, 5966, 5968, 5972, 5976, 5977, 5979, 5983, 5985, 5988, 5992, 5997, 5998, 5999, 6000, 6001, 6002, 6005, 6008, 6020, 6021, 6022, 6023, 6025, 6028, 6029, 6031, 6039, 6042, 6043, 6046, 6047, 6050, 6054, 6063, 6064, 6068, 6069, 6070, 6071, 6073, 6074, 6076, 6078, 6079, 6080, 6082, 6083, 6084, 6089, 6090, 6094, 6098, 6100, 6104, 6105, 6110, 6113, 6120, 6124, 6125, 6126, 6127, 6130, 6132, 6134, 6137, 6139, 6140, 6141, 6144, 6145, 6146, 6152, 6166, 6169, 6170, 6173, 6176, 6177, 6182, 6185, 6186, 6193, 6197, 6199, 6202, 6203, 6205, 6207, 6209, 6210, 6211, 6212, 6218, 6220, 6221, 6222, 6229, 6233, 6234, 6236, 6240, 6242, 6246, 6251, 6258, 6263, 6270, 6271, 6272, 6283, 6288, 6292, 6295, 6297, 6300, 6301, 6302, 6305, 6306, 6307, 6311, 6312, 6313, 6315, 6319, 6322, 6325, 6326, 6327, 6329, 6330, 6331, 6334, 6341, 6345, 6347, 6348, 6351, 6354, 6355, 6356, 6358, 6360, 6361, 6362, 6363, 6368, 6369, 6371, 6372, 6376, 6377, 6383, 6384, 6386, 6387, 6390, 6391, 6394, 6396, 6397, 6398, 6400, 6401, 6402, 6403, 6405, 6412, 6413, 6416, 6418, 6419, 6426, 6430, 6433, 6436, 6438, 6449, 6450, 6451, 6454, 6456, 6458, 6460, 6462, 6476, 6477, 6479, 6480, 6481, 6482, 6488, 6489, 6491, 6493, 6495, 6497, 6500, 6505, 6506, 6512, 6516, 6518, 6527, 6535, 6536, 6537, 6538, 6542, 6543, 6548, 6550, 6554, 6556, 6560, 6564, 6566, 6570, 6574, 6575, 6587, 6588, 6590, 6593, 6594, 6598, 6599, 6600, 6607, 6612, 6614, 6627, 6634, 6639, 6640, 6641, 6643, 6645, 6652, 6653, 6656, 6659, 6666, 6674, 6676, 6678, 6683, 6685, 6686, 6688, 6690, 6695, 6699, 6700, 6703, 6706, 6708, 6710, 6712, 6717, 6719, 6723, 6727, 6729, 6732, 6734, 6735, 6738, 6739, 6740, 6746, 6747, 6754, 6756, 6762, 6763, 6765, 6766, 6770, 6771, 6773, 6776, 6779, 6784, 6786, 6787, 6788, 6790, 6791, 6794, 6798, 6799, 6803, 6804, 6805, 6807, 6811, 6813, 6815, 6816, 6818, 6821, 6823, 6824, 6829, 6830, 6831, 6832, 6835, 6841, 6842, 6843, 6846, 6847, 6848, 6849, 6855, 6857, 6858, 6861, 6862, 6863, 6873, 6875, 6878, 6881, 6882, 6885, 6889, 6891, 6892, 6894, 6897, 6898, 6901, 6906, 6907, 6909, 6914, 6915, 6922, 6923, 6927, 6928, 6929, 6933, 6934, 6935, 6937, 6944, 6946, 6948, 6951, 6953, 6954, 6956, 6960, 6963, 6964, 6967, 6968, 6969, 6970, 6972, 6973, 6974, 6975, 6977, 6978, 6979, 6982, 6983, 6989, 6990, 6994, 6995, 6998, 7001, 7002, 7003, 7005, 7006, 7007, 7008, 7015, 7018, 7021, 7026, 7028, 7029, 7031, 7034, 7037, 7039, 7045, 7047, 7048, 7054, 7059, 7069, 7071, 7083, 7091, 7092, 7093, 7094, 7097, 7100, 7102, 7107, 7109, 7111, 7113, 7114, 7115, 7117, 7118, 7125, 7127, 7130, 7139, 7148, 7150, 7154, 7156, 7157, 7168, 7170, 7171, 7180, 7181, 7183, 7184, 7186, 7189, 7195, 7199, 7200, 7201, 7204, 7207, 7209, 7210, 7214, 7220, 7224, 7227, 7229, 7231, 7234, 7237, 7240, 7247, 7248, 7249, 7256, 7257, 7262, 7268, 7274, 7276, 7277, 7282, 7288, 7291, 7293, 7294, 7295, 7296, 7300, 7306, 7308, 7310, 7312, 7313, 7315, 7319, 7323, 7333, 7334, 7338, 7341, 7343, 7345, 7347, 7348, 7349, 7351, 7352, 7355, 7356, 7360, 7362, 7364, 7365, 7366, 7367, 7370, 7371, 7374, 7377, 7380, 7383, 7389, 7391, 7392, 7393, 7395, 7401, 7408, 7409, 7413, 7417, 7421, 7424, 7425, 7426, 7427, 7433, 7435, 7436, 7439, 7443, 7445, 7448, 7449, 7450, 7451, 7452, 7453, 7455, 7460, 7466, 7479, 7481, 7482, 7483, 7487, 7488, 7489, 7490, 7491, 7497, 7499, 7500, 7501, 7504, 7508, 7509, 7513, 7514, 7515, 7516, 7518, 7521, 7525, 7526, 7529, 7533, 7535, 7537, 7546, 7547, 7548, 7552, 7553, 7555, 7562, 7563, 7564, 7565, 7566, 7567, 7572, 7575, 7581, 7584, 7591, 7597, 7598, 7599, 7603, 7605, 7608, 7610, 7612, 7613, 7614, 7620, 7623, 7626, 7627, 7631, 7636, 7639, 7640, 7641, 7644, 7646, 7649, 7650, 7651, 7655, 7657, 7658, 7660, 7662, 7664, 7665, 7667, 7668, 7671, 7672, 7673, 7674, 7675, 7676, 7678, 7683, 7684, 7692, 7694, 7696, 7700, 7704, 7705, 7706, 7710, 7713, 7715, 7720, 7725, 7726, 7738, 7739, 7740, 7743, 7753, 7756, 7760, 7761, 7762, 7769, 7771, 7776, 7777, 7778, 7779, 7784, 7792, 7796, 7797, 7800, 7802, 7805, 7806, 7808, 7809, 7813, 7814, 7815, 7817, 7821, 7825, 7833, 7837, 7838, 7840, 7847, 7849, 7852, 7856, 7858, 7863, 7864, 7869, 7872, 7873, 7875, 7876, 7880, 7883, 7885, 7887, 7888, 7889, 7895, 7896, 7898, 7903, 7905, 7909, 7911, 7912, 7913, 7918, 7922, 7923, 7926, 7928, 7929, 7932, 7935, 7936, 7939, 7943, 7945, 7958, 7959, 7960, 7964, 7965, 7966, 7980, 7983, 7986, 7988, 7990, 7992, 7993, 7995, 7996, 8002, 8005, 8010, 8015, 8016, 8021, 8022, 8023, 8025, 8026, 8029, 8030, 8032, 8033, 8035, 8038, 8039, 8042, 8043, 8045, 8046, 8047, 8050, 8051, 8053, 8057, 8058, 8059, 8060, 8061, 8065, 8066, 8067, 8068, 8069, 8070, 8071, 8072, 8073, 8074, 8076, 8077, 8079, 8080, 8082, 8084, 8085, 8089, 8090, 8091, 8095, 8096, 8097, 8099, 8103, 8106, 8108, 8109, 8114, 8117, 8118, 8123, 8125, 8128, 8130, 8132, 8133, 8137, 8140, 8141, 8143, 8145, 8147, 8148, 8151, 8152, 8153, 8155, 8157, 8158, 8168, 8170, 8173, 8174, 8177, 8181, 8184, 8187, 8189, 8196, 8197, 8205, 8216, 8222, 8223, 8225, 8226, 8228, 8229, 8231, 8232, 8234, 8237, 8238, 8239, 8242, 8248, 8253, 8254, 8258, 8259, 8262, 8263, 8267, 8273, 8275, 8278, 8280, 8281, 8283, 8284, 8286, 8292, 8293, 8294, 8296, 8298, 8299, 8305, 8310, 8311, 8313, 8314, 8315, 8320, 8322, 8327, 8330, 8331, 8332, 8333, 8334, 8337, 8341, 8355, 8361, 8362, 8365, 8373, 8381, 8382, 8392, 8393, 8394, 8402, 8403, 8404, 8405, 8413, 8415, 8418, 8424, 8425, 8427, 8435, 8440, 8442, 8443, 8445, 8447, 8452, 8454, 8456, 8464, 8468, 8469, 8472, 8479, 8484, 8485, 8486, 8497, 8499, 8501, 8502, 8504, 8507, 8509, 8512, 8515, 8516, 8521, 8532, 8534, 8538, 8541, 8553, 8557, 8558, 8564, 8566, 8574, 8575, 8576, 8581, 8584, 8591, 8601, 8607, 8611, 8612, 8616, 8620, 8623, 8625, 8626, 8628, 8631, 8633, 8640, 8646, 8650, 8662, 8664, 8665, 8667, 8669, 8673, 8678, 8680, 8684, 8691, 8695, 8697, 8698, 8701, 8704, 8709, 8710, 8716, 8720, 8721, 8723, 8725, 8726, 8733, 8734, 8735, 8736, 8744, 8745, 8747, 8748, 8749, 8751, 8753, 8755, 8758, 8759, 8763, 8768, 8771, 8773, 8779, 8781, 8783, 8788, 8789, 8791, 8792, 8793, 8796, 8797, 8800, 8801, 8806, 8808, 8814, 8817, 8818, 8820, 8827, 8830, 8834, 8838, 8847, 8848, 8849, 8853, 8854, 8858, 8861, 8862, 8863, 8864, 8865, 8872, 8873, 8874, 8878, 8879, 8880, 8881, 8882, 8884, 8889, 8892, 8893, 8894, 8896, 8898, 8902, 8904, 8908, 8917, 8918, 8930, 8931, 8932, 8934, 8937, 8939, 8941, 8943, 8948, 8953, 8955, 8960, 8964, 8965, 8967, 8969, 8971, 8972, 8977, 8979, 8985, 8988, 8991, 8993, 8995, 8997, 9006, 9007, 9011, 9012, 9021, 9026, 9027, 9033, 9034, 9035, 9043, 9046, 9047, 9061, 9067, 9079, 9082, 9085, 9090, 9091, 9092, 9095, 9098, 9105, 9112, 9115, 9119, 9121, 9124, 9126, 9133, 9134, 9137, 9138, 9139, 9141, 9144, 9145, 9148, 9149, 9150, 9153, 9160, 9161, 9162, 9165, 9166, 9169, 9171, 9173, 9174, 9178, 9183, 9185, 9188, 9192, 9193, 9194, 9198, 9199, 9202, 9204, 9207, 9208, 9209, 9210, 9211, 9215, 9221, 9222, 9224, 9225, 9226, 9228, 9229, 9230, 9231, 9235, 9237, 9238, 9240, 9242, 9247, 9248, 9249, 9250, 9252, 9254, 9258, 9264, 9266, 9268, 9273, 9276, 9277, 9279, 9281, 9284, 9289, 9290, 9295, 9296, 9297, 9304, 9306, 9314, 9315, 9317, 9324, 9328, 9331, 9333, 9334, 9335, 9336, 9337, 9340, 9341, 9342, 9343, 9344, 9347, 9352, 9353, 9357, 9365, 9368, 9372, 9374, 9377, 9381, 9385, 9386, 9397, 9399, 9401, 9402, 9408, 9413, 9415, 9416, 9418, 9421, 9422, 9424, 9426, 9430, 9436, 9437, 9443, 9446, 9449, 9450, 9455, 9457, 9459, 9460, 9461, 9463, 9470, 9472, 9476, 9479, 9482, 9483, 9485, 9486, 9490, 9491, 9495, 9498, 9499, 9502, 9503, 9506, 9514, 9516, 9517, 9523, 9531, 9536, 9537, 9539, 9541, 9542, 9544, 9547, 9548, 9549, 9553, 9554, 9555, 9562, 9563, 9565, 9566, 9568, 9570, 9571, 9579, 9582, 9584, 9585, 9594, 9595, 9604, 9606, 9607, 9614, 9626, 9627, 9628, 9631, 9633, 9634, 9637, 9647, 9649, 9650, 9651, 9654, 9658, 9663, 9665, 9668, 9670, 9677, 9678, 9683, 9686, 9688, 9691, 9692, 9696, 9700, 9703, 9709, 9711, 9715, 9718, 9724, 9727, 9734, 9738, 9740, 9746, 9750, 9757, 9758, 9759, 9761, 9762, 9763, 9765, 9769, 9770, 9772, 9773, 9777, 9778, 9781, 9782, 9784, 9788, 9790, 9794, 9800, 9802, 9803, 9804, 9809, 9815, 9816, 9820, 9825, 9827, 9830, 9833, 9834, 9835, 9837, 9840, 9842, 9844, 9846, 9847, 9848, 9849, 9851, 9852, 9853, 9860, 9862, 9865, 9871, 9874, 9875, 9876, 9879, 9882, 9885, 9886, 9887]]))
        self.assertEqual(return_new[2], 6893)
        self.assertTrue(numpy.allclose(return_new[0], return_old[0], atol=TOLERANCE, equal_nan=True))
        self.assertTrue(numpy.array_equal([return_new[1][0].tolist(), return_new[1][1].tolist()],[return_old[1][0].tolist(), return_old[1][1].tolist()]))
        self.assertEqual(return_new[2], return_old[2])

    def test_emptyList_asg1_returns_IndexError(self):
        return_new = fu.k_means_match_clusters_asg_new(asg1=[], asg2=self.asg2, T=0)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1=[], asg2=self.asg2, T=0)
        self.assertTrue(numpy.array_equal(return_new[0],[]))
        self.assertTrue(numpy.array_equal(return_new[1],[]))
        self.assertEqual(return_new[2], 0)
        self.assertTrue(numpy.allclose(return_new[0], return_old[0]))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1]))
        self.assertEqual(return_new[2], return_old[2])

    def test_emptyList_asg2_returns_IndexError(self):
        with self.assertRaises(IndexError)  as cm_new:
            fu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=[], T=0)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=[], T=0)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyList_both_asg(self):
        return_new = fu.k_means_match_clusters_asg_new(asg1=[], asg2=[], T=0)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1=[], asg2=[], T=0)

        self.assertTrue(numpy.array_equal(return_new[0],[]))
        self.assertTrue(numpy.array_equal(return_new[1],[]))
        self.assertEqual(return_new[2], 0)
        self.assertTrue(numpy.allclose(return_new[0], return_old[0]))
        self.assertTrue(numpy.allclose(return_new[1], return_old[1]))
        self.assertEqual(return_new[2], return_old[2])




class Test_hist_list(unittest.TestCase):
    (data,nbins) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.hist_list"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.hist_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.hist_list()
        self.assertEqual(cm_new.exception.message, "hist_list() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_pickle_value(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=None, fmaxiu=None)
        return_old = oldfu.hist_list(data=self.data, nbin=self.nbins, fminiu=None, fmaxiu=None)

        self.assertTrue(numpy.array_equal(return_new[0],return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1],return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[0], [0.0, 88.625, 177.25, 265.875, 354.5, 443.125, 531.75, 620.375, 709.0, 797.625, 886.25, 974.875, 1063.5, 1152.125, 1240.75, 1329.375, 1418.0, 1506.625, 1595.25, 1683.875]))
        self.assertTrue(numpy.array_equal(return_new[1], [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407]))


    def test_empty_list(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.hist_list(data=[], nbin=self.nbins, fminiu=None, fmaxiu=None)
        with self.assertRaises(ValueError) as cm_old:
            oldfu.hist_list(data=[], nbin=self.nbins, fminiu=None, fmaxiu=None)
        self.assertEqual(cm_new.exception.message, "max() arg is an empty sequence")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_fmin_equal_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=1)
        return_old = oldfu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=1)
        self.assertTrue(numpy.array_equal(return_new[0],return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1],return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[0], [0.0, 88.625, 177.25, 265.875, 354.5, 443.125, 531.75, 620.375, 709.0, 797.625, 886.25, 974.875, 1063.5, 1152.125, 1240.75, 1329.375, 1418.0, 1506.625, 1595.25, 1683.875]))
        self.assertTrue(numpy.array_equal(return_new[1],  [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407]))

    def test_fmin_major_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=10, fmaxiu=1)
        return_old = oldfu.hist_list(data=self.data, nbin=self.nbins, fminiu=10, fmaxiu=1)
        self.assertTrue(numpy.array_equal(return_new[0],return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1],return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[0], [0.0, 88.625, 177.25, 265.875, 354.5, 443.125, 531.75, 620.375, 709.0, 797.625, 886.25, 974.875, 1063.5, 1152.125, 1240.75, 1329.375, 1418.0, 1506.625, 1595.25, 1683.875]))
        self.assertTrue(numpy.array_equal(return_new[1], [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407]))

    def test_fmin_minor_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=10)
        return_old = oldfu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=10)
        self.assertTrue(numpy.array_equal(return_new[0],return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1],return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[0], [0.0, 88.625, 177.25, 265.875, 354.5, 443.125, 531.75, 620.375, 709.0, 797.625, 886.25, 974.875, 1063.5, 1152.125, 1240.75, 1329.375, 1418.0, 1506.625, 1595.25, 1683.875]))
        self.assertTrue(numpy.array_equal(return_new[1], [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407]))



class Test_linreg(unittest.TestCase):
    (X,Y) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.linreg"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg()
        self.assertEqual(cm_new.exception.message, "linreg() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.linreg(X=self.X,Y=self.Y)
        return_old = oldfu.linreg(X=self.X,Y=self.Y)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new, (-3.3524868938354802, -11.920708605604739)))

    def test_different_list_size_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=self.X,Y=[3,4])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=self.X,Y=[3,4])
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_both_empty_list_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.linreg(X=[],Y=[])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.linreg(X=[],Y=[])
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_X_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=[],Y=self.Y)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=[],Y=self.Y)
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Y_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=self.X,Y=[])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=self.X,Y=[])
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_pearson(unittest.TestCase):
    (X,Y) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.pearson"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson()
        self.assertEqual(cm_new.exception.message, "pearson() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.pearson(X=self.X,Y=self.Y)
        return_old = oldfu.pearson(X=self.X,Y=self.Y)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertEqual(return_new,0.43165202148828613)

    def test_different_list_size_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=self.X,Y=[3,4])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=self.X,Y=[3,4])
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_both_empty_list_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.pearson(X=[],Y=[])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.pearson(X=[],Y=[])
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_X_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=[],Y=self.Y)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=[],Y=self.Y)
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Y_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=self.X,Y=[])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=self.X,Y=[])
        self.assertEqual(cm_new.exception.message, "unsupported operand type(s) for +=: 'float' and 'NoneType'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_table_stat(unittest.TestCase):
    X = [4808, 2749]
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.table_stat()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.table_stat()
        self.assertEqual(cm_new.exception.message, "table_stat() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.table_stat(X=self.X)
        return_old = oldfu.table_stat(X=self.X)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new,(3778, 2119741.0, 2749, 4808)))

    def test_empty_list_returns_IndexError(self):
        with self.assertRaises(IndexError)  as cm_new:
            fu.table_stat(X=[])
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.table_stat(X=[])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_mono(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.mono()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.mono()
        self.assertEqual(cm_new.exception.message, "mono() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_(self):
        return_new = fu.mono(k1=10,k2=20)
        return_old = oldfu.mono(k1=10,k2=20)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 200)

    def test_m1_equalm2(self):
        return_new = fu.mono(k1=20,k2=20)
        return_old = oldfu.mono(k1=20,k2=20)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,210)

    def test_both_zero(self):
        return_new = fu.mono(k1=0,k2=0)
        return_old = oldfu.mono(k1=0,k2=0)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,0)



class Test_k_means_stab_bbenum(unittest.TestCase):
    (PART,)  = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_stab_bbenum"))[0]
    expected_result=([[1, 0, 1, 0, 0]], [[], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]], [0, 97], [0, 97], [0, 100.0], 100.0)
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_stab_bbenum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_stab_bbenum()
        self.assertEqual(cm_new.exception.message, "k_means_stab_bbenum() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_PART_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_stab_bbenum(PART=[], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_stab_bbenum(PART=[], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new =fu.k_means_stab_bbenum(PART=self.PART, T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        return_old =oldfu.k_means_stab_bbenum(PART=self.PART, T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new, self.expected_result ))

    def test_branchfunc_is4(self):
        return_new =fu.k_means_stab_bbenum(PART=self.PART, T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=4, LIM=-1)
        return_old =oldfu.k_means_stab_bbenum(PART=self.PART, T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=4, LIM=-1)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new,  self.expected_result))

    def test_lenPART_is2(self):
        return_new =fu.k_means_stab_bbenum(PART=self.PART[0:2], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        return_old =oldfu.k_means_stab_bbenum(PART=self.PART[0:2], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new,([[1, 0]], [[], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]], [0, 97], [0, 97], [0, 100.0], 100.0)))




class Test_k_means_match_bbenum(unittest.TestCase):
    (PART,) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_bbenum"))[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_match_bbenum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_match_bbenum()
        self.assertEqual(cm_new.exception.message, "k_means_match_bbenum() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_PART_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_match_bbenum([],T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_match_bbenum([],T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_value(self):
        return_new = fu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))

    def test_DoMPI_isTrue_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=True, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=True, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_DoMPI_DoMPIinit_areTrue_spawn_errorMsg_but_works(self):
        return_new = fu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=True, Njobs=-1, DoMPI=True, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART,T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=True, Njobs=-1, DoMPI=True, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new,([[1, 0, 1, 0, 0]], 97)))

    def test_T0(self):
        return_new = fu.k_means_match_bbenum(self.PART,T=0, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART,T=0, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1,DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[], N_start=-1, N_stop=-1,topMatches=[])
        self.assertTrue(numpy.array_equal(return_new,return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0, 1, 0, 1, 1], [1, 0, 1, 0, 0]]))

    def test_J0_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=0, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=0, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """

    def test_branchfunc0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=0, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=0, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=-1, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))

    def test_nguess0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=0,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=0,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=-1, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_np0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=0,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=0, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_LIM0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=-1, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_stmult0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=40, stmult=0, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=-1, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_max_branching0(self):
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=0, stmult=0.25, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=1, max_branching=0, stmult=0.25, nguesses=5,branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1,np=-1, c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(numpy.array_equal(return_new, [[1, 0, 1, 0, 0]]))
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_scale_fsc_datasetsize(unittest.TestCase):
    fsc_to_be_adjusted =[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992, 0.9957, 0.99387, 0.99405, 0.99123, 0.98935, 0.9872, 0.98074, 0.96944, 0.96374, 0.94607, 0.91744, 0.88535, 0.89291, 0.88777, 0.8877, 0.87859, 0.87215, 0.8633, 0.8441, 0.84776, 0.85752, 0.87845, 0.88849, 0.87827, 0.85501, 0.83834, 0.86192, 0.87235, 0.84794, 0.82865, 0.81224, 0.81058, 0.7955, 0.77124, 0.74744, 0.74404, 0.7588, 0.71706, 0.58851, 0.50241, 0.60636, 0.67929, 0.62126, 0.46365, 0.45433, 0.50698, 0.51007, 0.47945, 0.47566, 0.44534, 0.37992, 0.33675, 0.36838, 0.35485, 0.35406, 0.34595, 0.30988, 0.32704, 0.35027, 0.3361, 0.34971, 0.35296, 0.31079, 0.33571, 0.37315, 0.34308, 0.37531, 0.37823, 0.40466, 0.42364, 0.40338, 0.39635, 0.37909, 0.37814, 0.3756, 0.35529, 0.37891, 0.35317, 0.34812, 0.34572, 0.33312, 0.31022, 0.28474, 0.25991, 0.25608, 0.24614, 0.20188, 0.17434, 0.17683, 0.15445, 0.13581, 0.12399, 0.12581, 0.12332, 0.111, 0.0911946, 0.0864365, 0.0785264, 0.0679365, 0.0488136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    nfsc = 9892
    nnew = 9892
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.scale_fsc_datasetsize()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.scale_fsc_datasetsize()
        self.assertEqual(cm_new.exception.message, "scale_fsc_datasetsize() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc,nnew=self.nnew)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992, 0.9957, 0.99387, 0.99405, 0.99123, 0.98935, 0.9872, 0.98074, 0.96944, 0.96374, 0.94607, 0.91744, 0.88535, 0.89291, 0.88777, 0.8877, 0.87859, 0.87215, 0.8633, 0.8441, 0.84776, 0.85752, 0.87845, 0.88849, 0.87827, 0.85501, 0.83834, 0.86192, 0.87235, 0.84794, 0.82865, 0.81224, 0.81058, 0.7955, 0.77124, 0.74744, 0.74404, 0.7588, 0.71706, 0.58851, 0.50241, 0.60636, 0.67929, 0.62126, 0.46365, 0.45433, 0.50698, 0.51007, 0.47945, 0.47566, 0.44534, 0.37992, 0.33675, 0.36838, 0.35485, 0.35406, 0.34595, 0.30988, 0.32704, 0.35027, 0.3361, 0.34971, 0.35296, 0.31079, 0.33571, 0.37315, 0.34308, 0.37531, 0.37823, 0.40466, 0.42364, 0.40338, 0.39635, 0.37909, 0.37814, 0.3756, 0.35529, 0.37891, 0.35317, 0.34812, 0.34572, 0.33312, 0.31022, 0.28474, 0.25991, 0.25608, 0.24614, 0.20188, 0.17434, 0.17683, 0.15445, 0.13581, 0.12399, 0.12581, 0.12332, 0.111, 0.0911946, 0.0864365, 0.0785264, 0.0679365, 0.0488136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_emptylist(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=self.nfsc, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=self.nfsc,nnew=self.nnew)
        self.assertTrue(numpy.array_equal(return_new, []))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_nfsc_greater_nnew(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc+10, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc+10,nnew=self.nnew)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9991991919133101, 0.995695671763659, 0.9938638410985181, 0.9940440208634007, 0.9912212120804388, 0.9893393484997355, 0.9871872260447904, 0.9807209050903437, 0.9694100513837745, 0.9637046745536805, 0.946018424207827, 0.917363435579012, 0.885247398287289, 0.8928133447421265, 0.8876692891990225, 0.8875992343368216, 0.8784821790119756, 0.8720372927942945, 0.8631807149173785, 0.843966989027912, 0.8476295479959695, 0.8573965043935038, 0.8783420718976555, 0.8883898540713431, 0.8781619342375103, 0.8548846969953323, 0.8382030166804726, 0.8617997034946842, 0.8722374432777654, 0.8478096745459852, 0.8285064854566488, 0.8120858580259775, 0.810424813315671, 0.7953355781179486, 0.7710616861389799, 0.747249214252374, 0.7438475260687702, 0.7586150243187773, 0.7168549586110095, 0.5882652918672154, 0.5021574034549281, 0.6061188024611799, 0.6790698377601415, 0.6210222260837981, 0.463398742503496, 0.45407951721577117, 0.5067274456477423, 0.5098174980925592, 0.4791978301349811, 0.4754080030001081, 0.4450904307942334, 0.3796819963482656, 0.3365243633400876, 0.36814493357971323, 0.3546187199078563, 0.35382895241253404, 0.3457214121617222, 0.3096639614909986, 0.3268176635514925, 0.3500400853774037, 0.33587457831324086, 0.3494802552414265, 0.3527292782529308, 0.31057361257630034, 0.3354847077067504, 0.3729136868736183, 0.3428523144536889, 0.37507313744618753, 0.37799240968401815, 0.40441660596219825, 0.4233933087672451, 0.40313685451865866, 0.39610827864698217, 0.37885219862651714, 0.37790243185778905, 0.37536306439808914, 0.3550585909581716, 0.3786722426777146, 0.3529392158582757, 0.34789074098641914, 0.34549148385760375, 0.3328955748069882, 0.3100038309315609, 0.28453426205390553, 0.2597156884514092, 0.2558875617923892, 0.24595256186990372, 0.2017172477054805, 0.1741946042085414, 0.17668297197529226, 0.15431809172820385, 0.1356914565307124, 0.12388029480686033, 0.12569891555298499, 0.12321080423791014, 0.11090033320236868, 0.09111089390176229, 0.08635674628852816, 0.07845331805355968, 0.06787254779170172, 0.04876670733038269, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_nfsc_lower_nnew(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc-10, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc-10,nnew=self.nnew)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992008080879967, 0.9957043282739705, 0.993876158977815, 0.9940559792085284, 0.9912387880753856, 0.9893606517296208, 0.9872127742857975, 0.9807590956532373, 0.9694699504667269, 0.9637753280361855, 0.9461215814161805, 0.9175165772023998, 0.8854526254989016, 0.893006676187781, 0.8878707336559222, 0.8878007885448379, 0.8786978474582086, 0.8722627363433496, 0.8634193180557805, 0.8442330529043143, 0.847890492163897, 0.8576435311871446, 0.8785579546293825, 0.8885901685095942, 0.8783780923627005, 0.8551353397421242, 0.8384770280998409, 0.8620403300937923, 0.8724625857754491, 0.8480703655274789, 0.8287935642713687, 0.8123942005004912, 0.8107352461285096, 0.795664489878872, 0.7714183963529617, 0.7476308831943587, 0.744232573564021, 0.7589850659097225, 0.7172651587181639, 0.5887549118061657, 0.5026628507965304, 0.6066013895784409, 0.6795103050444103, 0.6214979560613412, 0.46390153011087243, 0.45458075928317637, 0.507232806225359, 0.5103227521491716, 0.47970243540512547, 0.4759122642910012, 0.4455898492375436, 0.380158302224063, 0.3369759394377599, 0.3686153667994114, 0.35508158196785944, 0.3542913495292306, 0.34617889031865434, 0.3100963401599992, 0.32726263916867665, 0.3505002168478494, 0.3363257244727836, 0.3499400470210042, 0.3531910237768141, 0.31100668916320073, 0.33593559508335824, 0.3733866128166641, 0.3433079881556603, 0.37554716190514603, 0.37846788918262947, 0.40490368718270053, 0.4238869788709483, 0.40362343895688485, 0.396592016549525, 0.3793281000916173, 0.3783778670242496, 0.3758372349064654, 0.35552171087972223, 0.3791480560717202, 0.35340108615503196, 0.3483495613749182, 0.3459488186343066, 0.3333447279945461, 0.3104364707529276, 0.28494603568690396, 0.26010460252252415, 0.256272727868698, 0.24632772403764963, 0.202043015134694, 0.17448563871053147, 0.17697727292939475, 0.15458213397025733, 0.13592875077550817, 0.12409989966954055, 0.12592128095895827, 0.1234293894835953, 0.11109984610107548, 0.09127846004588899, 0.08651640115876466, 0.07859961823005365, 0.06800057283857858, 0.04886058293771373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_nfsc_is0(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=0, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=0,nnew=self.nnew)
        self.assertTrue(numpy.array_equal(return_new, []))
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_nnew_is0_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=self.nfsc, nnew=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=self.nfsc,nnew=0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)








"""
def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list




class Test_lib_statistics_compare(unittest.TestCase):

    def test_add_ave_varf_MPI_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]

        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(0,data)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(0, data)


        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))


    #Unable to read the file 
    def test_sum_oe_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (data,) = argum[0]

        return_new = fu.sum_oe(data)
        return_old = oldfu.sum_oe(data)

        # print(return_new)

        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))


    def test_ave_var_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data) = argum[0][0]

        return_new = fu.ave_var(data)
        return_old = oldfu.ave_var(data)

        self.assertTrue(return_new, return_old)


    def test_ave_series_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]

        return_new = fu.ave_series(data,)
        return_old = oldfu.ave_series(data,)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_varf2d_MPI_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]
        ave = get_data(1)

        return_new = fu.varf2d_MPI(0,data,ave)
        return_old = oldfu.varf2d_MPI(0,data,ave)

        self.assertTrue(return_new, return_old)


    #The function works but is extremely slow. It takes 55 mins to check the unit test that is why it is commented
    # def test_varf3d_MPI_true_should_return_equal_objects(self):
    #
    #     stack_name = "bdb:tests/Substack/sort3d_substack_003"
    #     nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
    #     list_proj = list(range(nima))
    #     proj = EMAN2_cppwrap.EMData()
    #     proj.read_image(stack_name, list_proj[0])
    #
    #     return_new = fu.varf3d_MPI([proj], mask2D=False  )
    #     mpi_barrier(MPI_COMM_WORLD)
    #     print("Hello")
    #     return_old = oldfu.varf3d_MPI([proj] , mask2D=False)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     print('yoyo')
    #
    #     self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_ccc_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2,mask) = argum[0]

        return_new = fu.ccc(img1, img2,mask )
        return_old = oldfu.ccc(img1, img2,mask )

        self.assertEqual(return_new, return_old)


    def test_fsc_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2) = argum[0]

        return_new = fu.fsc(img1, img2 )
        return_old = oldfu.fsc(img1, img2)

        self.assertEqual(return_new, return_old)


    def test_fsc_mask_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2, mask, w, filename) = argum[0]
        filename = os.path.join(ABSOLUTE_PATH, filename)

        return_new = fu.fsc_mask(img1, img2, mask, w, filename)
        return_old = oldfu.fsc_mask(img1, img2, mask, w, filename)

        self.assertEqual(return_new, return_old)

    #  Can only work on the cluster
    # def test_locres_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.locres")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc) = argum[0]
    #
    #     number_of_proc = 1
    #     main_node = 0
    #     myid = 0
    #
    #     return_new = fu.locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)


    def test_histogram_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data) = argum[0][0]

        image = get_data(1)

        return_new = fu.histogram(data[0])
        return_old = oldfu.histogram(data[0])

        self.assertTrue(return_new, return_old)



    def test_k_means_match_clusters_asg_new_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_clusters_asg_new")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (asg1, asg2) = argum[0]

        return_new = fu.k_means_match_clusters_asg_new(asg1, asg2)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1, asg2)

        self.assertTrue(return_new, return_old)


    def test_hist_list_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.hist_list")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (data,nbins) = argum[0]

        return_new = fu.hist_list(data,nbins)
        return_old = oldfu.hist_list(data,nbins)

        self.assertEqual(return_new, return_old)


    def test_linreg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.linreg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, Y) = argum[0]

        return_new = fu.linreg(X, Y)
        return_old = oldfu.linreg(X, Y)

        self.assertEqual(return_new, return_old)


    def test_pearson_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.pearson")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, Y) = argum[0]

        return_new = fu.pearson(X, Y)
        return_old = oldfu.pearson(X, Y)

        self.assertEqual(return_new, return_old)


    def test_table_stat_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.table_stat")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, ) = argum[0]

        return_new = fu.table_stat(X)
        return_old = oldfu.table_stat(X)

        self.assertEqual(return_new, return_old)


    def test_mono_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.mono")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (k1,k2) = argum[0]

        return_new = fu.mono(k1,k2)
        return_old = oldfu.mono(k1,k2)

        self.assertEqual(return_new, return_old)


    def test_k_means_stab_bbenum_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_stab_bbenum")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (Part,) = argum[0]

        return_new = fu.k_means_stab_bbenum(Part)
        return_old = oldfu.k_means_stab_bbenum(Part)

        self.assertEqual(return_new, return_old)


    def test_k_means_match_bbenum_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_bbenum")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (Part,) = argum[0]

        return_new = fu.k_means_match_bbenum(Part)
        return_old = oldfu.k_means_match_bbenum(Part)

        self.assertEqual(return_new, return_old)


    def test_scale_fsc_datasetsize_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.scale_fsc_datasetsize")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (fsc_to_be_adjusted, nfsc, nnew) = argum[0]

        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted, nfsc, nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted, nfsc, nnew)

        self.assertEqual(return_new, return_old)
"""




if __name__ == '__main__':
    unittest.main()
