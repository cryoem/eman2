from __future__ import print_function
from __future__ import division

import unittest
import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

import os
import cPickle as pickle
ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from ..libpy import sparx_applications as fu
from .sparx_lib import sparx_applications as oldfu


def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


class Test_lb_compare(unittest.TestCase):


    # def test_ali2d_MPI_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #         print(argum[1])
    #         stack = argum[0][0]
    #         outdir = argum[0][1]
    #
    #         (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #          Fourvar, user_func_name, random_method, log ,num_of_proc, myid, main_mode, mpicom) = argum[0]
    #
    #     # outdir = 'files/Class2D/2dalignment'
    #
    #     return_new = fu.ali2d_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #          Fourvar)
    #     return_old = oldfu.ali2d_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #          Fourvar)
    #
    #     self.assertTrue(return_new, return_old)
    #
    # def test_ali2d_base_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0][1])
    #         print(argum[0][5])
    #         stack = argum[0][0]
    #         outdir = argum[0][1]
    #         outdir = 'files/Class2D/2dalignment'
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, num_of_proc, myid, main_mode, mpicom) = argum[0]
    #
    #     outdir = 'files/Class2D/2dalignment'
    #
    #     return_new = fu.ali2d_base(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #                               Fourvar, user_func_name, random_method, log, num_of_proc, myid, main_mode, mpicom)
    #     return_old = oldfu.ali2d_base(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomir, dst, center, maxit, CTF, snr,
    #                                  Fourvar, user_func_name, random_method, log, num_of_proc, myid, main_mode, mpicom)
    #
    #     self.assertEqual(return_new, return_old)

    # def test_cpy_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.cpy")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0][0])
    #         print(argum[0][1])
    #
    #     # (ins_list, ous) = argum[0]
    #
    #     ins_list = 'files/temp'
    #     ous = 'files/best_000'
    #
    #     return_new = fu.cpy(ins_list,ous)
    #     # return_old = oldfu.cpy(ins_list, ous)
    #     #
    #     # self.assertTrue(return_new, return_old)


    def test_ali_vol_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]

        return_new = fu.ali_vol(vol,refv,ang_scale,shift_scale,radius)
        return_old = oldfu.ali_vol(vol, refv, ang_scale, shift_scale, radius)

        self.assertTrue(return_new, return_old)


    def test_extract_values_true_should_return_equal_object(self):

        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')

        self.assertTrue(return_new, return_old)


    # def test_header_true_should_return_equal_object(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.header")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum)
    #         print(argum[1])
    #
    #     return_new = fu.header(stack='files/', fprint = True, params = 'original_image_index')


    def test_within_group_refinement_true_should_return_equal_object(self):
        filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/applications.within_group_refinement")
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
