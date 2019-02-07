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
import sys
import cPickle as pickle
ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from ..libpy import sparx_applications as fu
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


    # def test_ali2d_MPI_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #         print(argum[1])
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
    #
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir+'AA')
    #     print(outdirnew)
    #     print(outdir)
    #     # if os.path.isdir(os.path.join(ABSOLUTE_PATH , 'Class2D')):
    #     #     os.mkdir(ABSOLUTE_PATH + '/Class2D'+ '/2ddalignment')
    #
    #     number_of_proc = 1
    #     myid = 82
    #     main_node = 0
    #
    #     return_new = fu.ali2d_MPI(stack, outdirnew, maskfile=maskfile, ir=ir, ou=ou, rs=rs, xr=xr, yr=yr, \
    #                                ts=ts, nomirror=nomirror, dst=dst, center=center, maxit=maxit, CTF=True, snr=snr, \
    #                                Fourvar=Fourvar, user_func_name=user_func_name, random_method=random_method)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir + 'AB')
    #     return_old = oldfu.ali2d_MPI(stack, outdirnew, maskfile=maskfile, ir=ir, ou=ou, rs=rs, xr=xr, yr=yr, \
    #                                   ts=ts, nomirror=nomirror, dst=dst, center=center, maxit=maxit, CTF=True, snr=snr, \
    #                                   Fourvar=Fourvar, user_func_name=user_func_name, random_method=random_method)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)


    def test_ali2d_base_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[1])

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = argum[0]

        outdirnew = os.path.join(ABSOLUTE_PATH, outdir+'B')
        print(outdirnew)
        print(outdir)
        # if os.path.isdir(os.path.join(ABSOLUTE_PATH , 'Class2D')):
        #     os.mkdir(ABSOLUTE_PATH + '/Class2D'+ '/2ddalignment')

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
        mpi_finalize()

        self.assertTrue(return_new, return_old)

    # def test_Kmref_ali3d_MPI_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         # print(argum[0])
    #         # print(argum[0][1])
    #
    #     (vol,refv,ang_scale,shift_scale,radius) = argum[0]
    #
    #
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
    #     return_new = fu.Kmref_ali3d_MPI(stack, vol, outdirnew,  maskfile, focus = None, \
    #                                     ir = ir, ou = ou, rs = rs, xr = xr, yr = yr, ts = ts )
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdir = "Class2D/Kmref_aligB"
    #     outdirnew = os.path.join(ABSOLUTE_PATH, outdir)
    #
    #     return_old = oldfu.Kmref_ali3d_MPI(stack, vol, outdirnew, maskfile, focus = None, \
    #                                        ir=ir, ou=ou, rs=rs, xr=xr, yr=yr, ts=ts)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)


    # def test_cpy_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.cpy")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0][1])
    #
    #     # ins_list = argum[0]
    #
    #     ins_list = os.path.join(ABSOLUTE_PATH,  'Class2D/best.hdf')
    #     ous = os.path.join(ABSOLUTE_PATH, 'Class2D/stack_ali2d')
    #
    #
    #     ins_list = 'best.hdf'
    #     ous = 'best_000'
    #
    #     return_new = fu.cpy(ins_list,ous)
    #     return_old = oldfu.cpy(ins_list, ous)
    #
    #     self.assertEqual(return_new, return_old)


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


    # def test_pca_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.header")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum)
    #         print(argum[1])
    #
    #     stack =  argum[0][0]
    #     params = argum[0][1]
    #
    #     return_new = fu.pca(stack)
    #     return_old = oldfu.pca(stack)
    #
    #     self.assertEqual(return_new, return_old)



    def test_prepare_2d_forPCA_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        return_new = fu.prepare_2d_forPCA(data = argum[0][0])
        return_old = oldfu.prepare_2d_forPCA(data = argum[0][0])

        self.assertTrue(return_new, return_old)



    def test_extract_values_true_should_return_equal_object(self):

        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')

        self.assertTrue(return_new, return_old)


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


    # def test_refvol_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #         print(argum[0][1])
    #
    #     (vol,refv,ang_scale,shift_scale,radius) = argum[0]
    #
    #
    #     vollist = [vol, vol, vol]
    #
    #     return_new = fu.refvol(vollist , vollist)




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

if __name__ == '__main__':
    unittest.main()
