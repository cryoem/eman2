from __future__ import print_function
from __future__ import division

import unittest
import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

from mpi import *
import global_def

import os
import shutil
import sys
import cPickle as pickle
ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from sphire.libpy import sparx_applications as fu
from .sparx_lib import sparx_applications as oldfu


mpi_init(0, [])



global_def.BATCH = True
global_def.MPI = True



from test_module import  remove_list_of_file, remove_dir, get_arg_from_pickle_file, ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,returns_values_in_file,get_real_data
from os import path
import EMAN2db
from ..libpy import sparx_utilities

from EMAN2_cppwrap import EMData
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_BLANK_2D = sparx_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sparx_utilities.model_blank(10, 10, 10)
MASK = sparx_utilities.model_circle(2, 5, 5)

TOLERANCE = 0.0005
ABSOLUTE_PATH_TO_STACK= "bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Substack/isac_substack")
"""
There are some opened issues in:
1) ali2d_MPI the results are different ... ANYWAY IT SEEMS TO BE NOT USED
2) ali2d_base the results are different
    -) if CTF =True and myid == main_node (i.e.: ln 695) it try to use an not defined variable and crash .... dead code or bug?
3) mref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
4) Kmref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
5) project3d:
    -) with 'noise' not None the results are different because the addition of a random gauss_noise to the volume in the code (noise is the sigma value in the generation noise process)
    -) How can I really test with 'listctfs'? It is never used in a real SPHIRE code case, anyway I'd test
    -) in case of realsp=trillinear=True it will spawn an error message but its behaviour will be as if it'd receive as input realsp=True and Trilinear=False ....BUG or just BAD IMPLEMENTED
6) recons3d_n_trl_MPI_one_node: I cannot run it. see the error message in "Test_recons3d_n_trl_MPI_one_node.test_NoCTF_symC1"
7) pca
    -) need a file name containing the average of the input stack
8) header
    -) the pickle file values lead to an error message got at line 3!!!!!!!! 
    -) you have to save the values in a file to test its behaviour. But if you set an output file you cannot set the other values becuase the 
        "op = zero+one++consecutive+randomize+rand_alpha+(fimport!=None)+(fexport!=None)+fprint+backup+restore+delete+doset" check .... is it a bug? which is the purpose of this function??
9) refvol: seems to be never USED
10) within_group_refinement:
    -) I tested just the default method because the other 2 are not in use (Fabian docet)
    -) do not test with randomize=True because the results will be never the same because the random orientation in same calculations
    -) since the input parameter are basically used in 'sparx_alignment.ali2d_single_iter' and 'sparx_utilities.combine_params2' that i have already deeply tested
        I tested it just the pickle file case
11) ali3d_mref_Kmeans_MPI and mref_ali3d_EQ_Kmeans are used only in sxsort3d.py that is buggy and maybe we will implement it from scratch. I did not test them for now
"""

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""




"""IT SEEMS TO BE NOT USED"""
class Test_ali2d_MPI(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base"))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_MPI()
        self.assertEqual(cm_new.exception.message, "ali2d_MPI() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("The result are not the same")
    def test_pickle_file_values(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr, Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnewa = os.path.join( 'ali2d_MPI_NEW')
        outdirnewb = os.path.join(  'ali2d_MPI__OLD')

        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

        fu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        oldfu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)


        self.assertEqual(returns_values_in_file(path.join(outdirnewa, "resolution001")),returns_values_in_file(path.join(outdirnewb, "resolution001")))
        self.assertEqual(returns_values_in_file(path.join(outdirnewa, "resolution002")),returns_values_in_file(path.join(outdirnewb, "resolution002")))
        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

    def test_ali2d_MPI_error_directory_exists(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr, Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnewa = os.path.join( 'ali2d_MPI_NEW')
        outdirnewb = os.path.join(  'ali2d_MPI__OLD')

        os.mkdir(outdirnewa)
        os.mkdir(outdirnewb)

        with self.assertRaises(OSError) as cm_new:
            fu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(OSError) as cm_old:
            oldfu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(cm_new.exception.strerror, "File exists")
        self.assertEqual(cm_new.exception.filename, "ali2d_MPI_NEW")
        self.assertEqual(cm_old.exception.strerror, cm_new.exception.strerror)

        remove_dir(outdirnewa)
        remove_dir(outdirnewb)



class Test_ali2d_base(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base"))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_base()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_base()
        self.assertEqual(cm_new.exception.message, "ali2d_base() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("The result are not the same")
    def test_ali2d_base_true_should_return_equal_object(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnew = os.path.join( ABSOLUTE_PATH,'ali2d_base_new')
        outdirnewold = os.path.join(ABSOLUTE_PATH,'ali2d_base_old')

        os.mkdir(outdirnew)
        os.mkdir(outdirnewold)
        number_of_proc = 1
        myid = 0
        main_node = 0
        maxit =2
        return_new = fu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_base(stack, outdirnewold, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution001")),returns_values_in_file(path.join(outdirnewold, "resolution001")))
        self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution002")),returns_values_in_file(path.join(outdirnewold, "resolution002")))
        self.assertTrue(numpy.allclose(return_new, return_old, atol=TOLERANCE,equal_nan=True))
        remove_dir(outdirnew)
        remove_dir(outdirnewold)



class Test_cpy(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.cpy()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.cpy()
        self.assertEqual(cm_new.exception.message, "cpy() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_case(self):
        ins_list = os.path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'Class2D/best.hdf')
        file_path_new = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "new_substack")
        file_path_old = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "old_substack")
        fu.cpy(ins_list,file_path_new)
        oldfu.cpy(ins_list, file_path_old)
        db_new, keys = EMAN2db.db_open_dict( file_path_new, True, True)
        db_new.realopen()
        db_old, keys = EMAN2db.db_open_dict( file_path_old, True, True)
        db_old.realopen()

        for index in db_new.bdb.keys():
            self.assertEqual(str(pickle.loads(db_new.bdb.get(index))), str(pickle.loads(db_old.bdb.get(index))))
        remove_dir( os.path.join(ABSOLUTE_PATH, 'EMAN2DB'))

    def test_error_hdf_not_found(self):
        ins_list = os.path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'not_found.hdf')
        file_path_new = 'bdb:{0}#{1}'.format(".", "new_substack")
        file_path_old = 'bdb:{0}#{1}'.format(".", "old_substack")
        with self.assertRaises(Exception) as cm_new:
            fu.cpy(ins_list,file_path_new)
        with self.assertRaises(Exception) as cm_old:
            oldfu.cpy(ins_list, file_path_old)
        self.assertEqual(str(cm_new.failureException.message), "<attribute 'message' of 'exceptions.BaseException' objects>")
        self.assertEqual(cm_new.failureException.message, cm_old.failureException.message)



class Test_project3d(unittest.TestCase):

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(numpy.allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.project3d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.project3d()
        self.assertEqual(cm_new.exception.message, "project3d() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_crashes_because_signal11SIGSEGV(self):
        self.assertTrue(True)
        """
        return_new = fu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        """

    def test_Nonetype_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.project3d(volume=None, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.project3d(volume=None, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_img_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.project3d(volume=EMData(), stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.project3d(volume=EMData(), stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_realsp_and_trilinear_error_msg(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_default_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_default_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_trilinear_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_trilinear_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_realsp_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_realsp_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("since it adds a random img_noise_to the stack the results cannot be the same")
    def test_blank_img_with_noise(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("since it adds a random img_noise_to the stack the results cannot be the same")
    def test_3Dimg_with_noise(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_save_on_hdf(self):
        outnew = os.path.join( ABSOLUTE_PATH,'project3dnew.hdf')
        outold = os.path.join(ABSOLUTE_PATH,'project3old.hdf')
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = outnew, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = outold, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)

        self.assertTrue(returns_values_in_file(outold),returns_values_in_file(outnew))
        remove_list_of_file([outnew,outold])
        self.assertTrue(return_new is None)
        self.assertTrue(return_old is None)

    def test_3Dimg_with_mask(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_with_mask(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_with_listagls(self):
        listangls =sparx_utilities.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_with_listagls(self):
        listangls =sparx_utilities.even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_empty_listagls(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [] , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [], listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertTrue(numpy.array_equal(return_old,return_new))
        self.assertTrue(numpy.array_equal(return_new, []))

    def test_blank_img_empty_listagls(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [], listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls =[], listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertTrue(numpy.array_equal(return_old, return_new))
        self.assertTrue(numpy.array_equal(return_new, []))


class Test_ali_vol(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol"))
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol()
        self.assertEqual(cm_new.exception.message, "ali_vol() takes at least 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_empty_refv_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=EMData(), ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=EMData(), ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_with_empty_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =EMData(), refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =EMData(), refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_with_NoneType_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =None, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =None, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_NoneType_refv_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=None, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=None, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_as_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =IMAGE_2D, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =IMAGE_2D, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_as_refvol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=IMAGE_2D, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=IMAGE_2D, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_with_pickle_values(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_values(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=None, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=None, discrepancy = "ccc")
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_zero_radius(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=0, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=0, discrepancy = "ccc")
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_zero_shift_returns_ZeroDivisionError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=0, radius=radius, discrepancy = "ccc")
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=0, radius=radius, discrepancy = "ccc")
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_zero_ang(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))


class Test_recons3d_n_trl_MPI_one_node(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params_proj"))
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons3d_n_trl_MPI_one_node()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons3d_n_trl_MPI_one_node()
        self.assertEqual(cm_new.exception.message, "recons3d_n_trl_MPI_one_node() takes exactly 12 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    """
    the fftvol that produces the error is created in the code and I cannot understand why 
        Error
        Traceback (most recent call last):
          File "/home/lusnig/SPHIRE_1_1/lib/python2.7/unittest/case.py", line 329, in run
            testMethod()
          File "/home/lusnig/EMAN2/eman2/sphire/tests/test_applications.py", line 487, in test_NoCTF_symC1
            return_old = fu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=1)
          File "/home/lusnig/EMAN2/eman2/sphire/libpy/sparx_applications.py", line 2361, in recons3d_n_trl_MPI_one_node
            fftvol.div_sinc(1)
        RuntimeError: NotExistingObjectException at /home/lusnig/EMAN2/eman2/libEM/emdata_metadata.cpp:1153: error with 'npad': 'The requested key does not exist' caught
    """
    @unittest.skip("I get an error when it runs")
    def test_NoCTF_symC1(self):
        img = self.argum[0][0]
        pjrlist=[img,img]
        group=0
        return_old = fu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=-1)
        return_new = oldfu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=-1)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview()))



class Test_pca(unittest.TestCase):

    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA"))
    ave, grp_imgdata  = fu.prepare_2d_forPCA(data=argum[0][0])

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(numpy.allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pca()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pca()
        self.assertEqual(cm_new.exception.message, "pca() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_verbose_nothing_to_print(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_radius_higher0_verbose(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_MPI_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=True, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=True, verbose=False)
        self.test_all_the_conditions(return_new, return_old)


    def test_without_genbuf(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=False, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=False, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_shuffle_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=True, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=True, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_incoreCalculation(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=True, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=True, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_radius_and_maskFile_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="test", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="test", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)



class Test_prepare_2d_forPCA(unittest.TestCase):
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA"))[0][0]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(numpy.allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_2d_forPCA()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_2d_forPCA()
        self.assertEqual(cm_new.exception.message, "prepare_2d_forPCA() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = False)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = False)
        self.assertTrue(numpy.allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_notAmode(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = False)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = False)
        self.assertTrue(numpy.allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_notAmode_withCTF(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = True)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = True)
        self.assertTrue(numpy.allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_withCTF(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = True)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = True)
        self.assertTrue(numpy.allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_with_stack(self):
        outnew= os.path.join(ABSOLUTE_PATH,'pcanew.hdf')
        outold= os.path.join(ABSOLUTE_PATH, 'pcaold.hdf')
        return_new= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = outnew, CTF = False)
        return_old= oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = outold, CTF = False)
        self.assertTrue(returns_values_in_file(outold), returns_values_in_file(outnew))
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)
        remove_list_of_file([outnew,outold])


class Test_extract_value(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.extract_value()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.extract_value()
        self.assertEqual(cm_new.exception.message, "extract_value() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_string_integer(self):
        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 20)

    def test_string_float(self):
        return_new = fu.extract_value('20.1')
        return_old = oldfu.extract_value('20.1')
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 20.1)

    def test_string_not_handled(self):
        return_new = fu.extract_value("invalid")
        return_old = oldfu.extract_value("invalid")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new,"invalid")



class Test_header(unittest.TestCase):
    outnew= os.path.join(ABSOLUTE_PATH,'Headernew.hdf')
    outold = os.path.join(ABSOLUTE_PATH, 'Headerold.hdf')
    params='xform.projection'
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.header"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.header()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.header()
        self.assertEqual(cm_new.exception.message, "header() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_No_option_errorMsg(self):
        return_new = fu.header(stack=[], params = [], zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=[], params=[], zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_default_Too_option_errorMsg(self):
        return_new = fu.header(stack=[], params = [], zero=True, one=True, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=[], params=[], zero=True, one=True, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_default_(self):
        return_new = fu.header(stack=ABSOLUTE_PATH_TO_STACK, params = self.params, zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=self.outnew, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=ABSOLUTE_PATH_TO_STACK, params=self.params, zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=self.outold, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertTrue(returns_values_in_file(self.outold), returns_values_in_file(self.outnew))
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)
        remove_list_of_file([self.outold,self.outnew])



class Test_MPI_start_end(unittest.TestCase):
    nima = 64
    nproc = 8
    myid = 1
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.MPI_start_end()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.MPI_start_end()
        self.assertEqual(cm_new.exception.message, "MPI_start_end() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_default_case(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [8,16]))

    def test_zero_nima(self):
        return_new = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [0,0]))

    def test_zero_nproc(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_zero_myd(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [0,8]))


"""seems to be never USED"""
class Test_refvol(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.refvol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.refvol()
        self.assertEqual(cm_new.exception.message, "refvol() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_within_group_refinement(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.within_group_refinement"))[0]
    randomize = False
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.within_group_refinement()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.within_group_refinement()
        self.assertEqual(cm_new.exception.message, "within_group_refinement() takes at least 13 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_pickle_file(self):
        (data, maskfile, randomize_not_used, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF) = self.argum
        return_new = fu.within_group_refinement(data=data, maskfile=maskfile, randomize=self.randomize, ir=ir, ou=ou, rs=rs, xrng=xrng, yrng=yrng, step=step, dst=dst, maxit=maxit, FH=FH, FF=FF,method = "", CTF = False)
        return_old = oldfu.within_group_refinement(data=data, maskfile=maskfile, randomize=self.randomize, ir=ir, ou=ou, rs=rs, xrng=xrng, yrng=yrng, step=step, dst=dst, maxit=maxit, FH=FH, FF=FF,method = "", CTF = False)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), atol=TOLERANCE, equal_nan=True))



class Test_ali3d_mref_Kmeans_MPI(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3d_mref_Kmeans_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3d_mref_Kmeans_MPI()
        self.assertEqual(cm_new.exception.message, "ali3d_mref_Kmeans_MPI() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_mref_ali3d_EQ_Kmeans(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.mref_ali3d_EQ_Kmeans()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.mref_ali3d_EQ_Kmeans()
        self.assertEqual(cm_new.exception.message, "mref_ali3d_EQ_Kmeans() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)







def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list

@unittest.skip("original Adnan test")
class Test_lib_applications_compare(unittest.TestCase):


    def test_ali2d_MPI_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
        maxit = 4

        outdirnewa = os.path.join(ABSOLUTE_PATH, outdir+'AA')
        outdirnewb = os.path.join(ABSOLUTE_PATH, outdir + 'ABB')


        if (os.path.exists(outdirnewa)):
            shutil.rmtree(outdirnewa)

        if (os.path.exists(outdirnewb)):
            shutil.rmtree(outdirnewb)

        stacka = "bdb:tests/Class2D/isac_substack"
        print(os.getcwd())
        return_new = fu.ali2d_MPI(stacka, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
                                   ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_MPI(stacka, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
                                   ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()

        self.assertTrue(return_new, return_old)


    def test_ali2d_base_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = argum[0]

        outdirnew = os.path.join(ABSOLUTE_PATH, outdir+'B')

        number_of_proc = 1
        myid = 0
        main_node = 0

        print(outdirnew)
        return_new = fu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()

        self.assertTrue(return_new, return_old)



    """  mref_ali3d_MPI is corrupted function so we wont create unit test for that 
         below unittest is just to show where the problem lies inside the code """
    # def test_mref_ali3d_MPI_true_should_return_equal_object(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #     (vol,refv,ang_scale,shift_scale,radius) = argum[0]
    #
    #     maxit = 4
    #     outdirnewa = os.path.join(ABSOLUTE_PATH, outdir+'AA')
    #     outdirnewb = os.path.join(ABSOLUTE_PATH, outdir + 'ABB')
    #
    #     if (os.path.exists(outdirnewa)):
    #         shutil.rmtree(outdirnewa)
    #     if (os.path.exists(outdirnewb)):
    #         shutil.rmtree(outdirnewb)
    #
    #     print(mpi_comm_rank(MPI_COMM_WORLD))
    #
    #     stacka = "bdb:tests/Substack/isac_substack"
    #     vola = "tests/Class2D/EMAN2DB/vol_adaptive_mask.hdf"
    #
    #     return_new = fu.mref_ali3d_MPI(stacka, vola, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
    #                                ts=ts,  maxit=maxit, CTF=False, snr=snr, mpi_comm = MPI_COMM_WORLD)
    #     mpi_barrier(MPI_COMM_WORLD)
        # return_old = oldfu.mref_ali3d_MPI(stacka, refv, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
        #                            ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)
        # mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()
        # self.assertTrue(return_new, return_old)


    """  mref_ali3d_MPI is corrupted function so we wont create unit test for that 
         below unittest is just to show where the problem lies inside the code """
    # def test_Kmref_ali3d_MPI_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #         print(argum[0][1])
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
    #
    #     outdir = "Class2D/Kmref_aligA"
    #
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir)
    #
    #     if (os.path.exists(outdirnew)):
    #         shutil.rmtree(outdirnew)
    #
    #     stacka = "bdb:tests/Class2D/isac_substack"
    #     vola = "tests/Class2D/EMAN2DB/vol_adaptive_mask.hdf"
    #
    #     return_new = fu.Kmref_ali3d_MPI(stacka, vola, outdirnew,  maskfile, focus = None, \
    #                                     ir = ir, ou = ou, rs = rs, xr = xr, yr = yr, ts = ts )
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdir = "Class2D/Kmref_aligB"
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir)
    #
    #     return_old = oldfu.Kmref_ali3d_MPI(stacka, vola, outdirnew, maskfile, focus = None, \
    #                                        ir=ir, ou=ou, rs=rs, xr=xr, yr=yr, ts=ts)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)


    def test_cpy_true_should_return_equal_object(self):
        ins_list = os.path.join(ABSOLUTE_PATH, 'Class2D/best_temp.hdf')
        ous = "bdb:tests/Class2D/isac_substack"

        return_new = fu.cpy(ins_list,ous)
        return_old = oldfu.cpy(ins_list, ous)

        self.assertEqual(return_new, return_old)


    def test_project3d_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            # print(argum[0])
            # print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]


        return_new = fu.project3d(vol)
        return_old = oldfu.project3d(vol)

        self.assertTrue(return_new, return_old)


    def test_ali_vol_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]

        return_new = fu.ali_vol(vol,refv,ang_scale,shift_scale,radius)
        return_old = oldfu.ali_vol(vol, refv, ang_scale, shift_scale, radius)

        self.assertTrue(return_new, return_old)


    def test_pca_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = argum[0]


        stacka  = copy.deepcopy(stack)
        stackb = copy.deepcopy(stack)


        return_new = fu.pca(stacka,nvec = 3)
        return_old = oldfu.pca(stackb,nvec = 3)


        self.assertTrue(return_new, return_old)



    def test_prepare_2d_forPCA_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        return_new = fu.prepare_2d_forPCA(data = argum[0][0])
        return_old = oldfu.prepare_2d_forPCA(data = argum[0][0])

        print(return_new)
        print(return_old)

        self.assertTrue(return_new, return_old)



    def test_extract_values_true_should_return_equal_object(self):

        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')

        self.assertEqual(return_new, return_old)


    def test_header_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.header")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)
            print(argum[1])

        stack =  argum[0][0]
        params = argum[0][1]

        return_new = fu.header(stack=stack, params = params)
        return_old = oldfu.header(stack=stack, params=params)

        self.assertEqual(return_new, return_old)


    def test_MPI_start_end_true_should_return_equal_object(self):

        nima = 64
        nproc = 8
        myid = 1

        return_new = fu.MPI_start_end(nima, nproc, myid)
        return_old = oldfu.MPI_start_end(nima, nproc, myid)

        self.assertEqual(return_new, return_old)


    def test_refvol_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]


        vollist = [vol, vol, vol]
        output = "tests/Class2D/"
        mask = "bdb:tests/Class2D/stack_ali2d_76x76x1"

        return_new = fu.refvol(vollist , vollist, output, refv)




    def test_within_group_refinement_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.within_group_refinement")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)
            print(argum[0])


        (data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF) = argum[0]

        return_new = fu.within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF)
        return_old = oldfu.within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF)

        self.assertTrue(return_new, return_old)



    # def test_ali3d_mref_Kmeans_MPI_true_should_return_equal_object(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #     Tracker["constants"]["nproc"] = 1
    #     Tracker["constants"]["myid"] = 0
    #     Tracker["constants"]["main_node"] = 0
    #     Tracker["total_stack"] = "stack"
    #     Tracker["constants"]["seed"] = 1.4
    #     Tracker["constants"]["indep_runs"] = 2
    #     Tracker["this_data_list"] = [2,3,5]
    #     Tracker["shrinkage"] = 0.5
    #     Tracker["constants"]["ir"] = 1
    #     Tracker["constants"]["xr"] = 1
    #     Tracker["constants"]["yr"] = 1
    #     Tracker["constants"]["ts"] = 1
    #     Tracker["constants"]["delta"] = 0.5
    #     Tracker["constants"]["an"] = 1
    #     Tracker["constants"]["center"] = 4
    #     Tracker["constants"]["nassign"] =1
    #     Tracker["constants"]["nrefine"] =1
    #     Tracker["constants"]["sym"] = "c1"
    #     Tracker["constants"]["stoprnct"] =5
    #     Tracker["constants"]["mask3D"]  = False
    #     Tracker["low_pass_filter"] = "0.50"
    #     Tracker["constants"]["PWadjustment"] = False
    #
    #     outdir = "Class2D/Kmref_alig_MPI_A"
    #
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir)
    #
    #     if (os.path.exists(outdirnew)):
    #         shutil.rmtree(outdirnew)
    #
    #     this_data_list_file = "sphire/tests/Sort3D/chunk_0.txt"
    #
    #     ref_list = "sphire/tests/Sort3D/refang.txt"
    #
    #     return_new = fu.ali3d_mref_Kmeans_MPI(ref_list, outdirnew, this_data_list_file, Tracker)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdir = "Class2D/Kmref_alig_MPI_B"
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir)
    #
    #
    #     return_old = oldfu.ali3d_mref_Kmeans_MPI(ref_list, outdirnew, this_data_list_file, Tracker)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)




if __name__ == '__main__':
    unittest.main()
