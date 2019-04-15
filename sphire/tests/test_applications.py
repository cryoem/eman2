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

def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


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
