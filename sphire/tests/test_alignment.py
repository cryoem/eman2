from __future__ import print_function
from __future__ import division

import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

import unittest


from ..libpy import sparx_alignment as fu
from .sparx_lib import sparx_alignment as oldfu

def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


class Test_lib_compare(unittest.TestCase):
    """"
    # data = list of images
    # numr = tuple or list precalcualte rings
    # wr = list of weights of numr
    # cs =  cs = [0.0]*2
    # tavg = blanck image
    # cnx = center value x
    # cny = cnx
    # xrng =  list of possible shifts
    # yrng =  list of possible shifts
    # step = stepsize of the shift
    """
    def test_ali2d_single_iter_true_should_return_equal_objects(self):
        a,b,c = get_data(3)
        data = [a]
        numr = [int(entry) for entry in numpy.arange(0, 10).tolist()]
        wr   = [0.0] * int(len(numr)/3)    #list(numpy.arange(0, 4))
        xrng = [int(entry) for entry in numpy.arange(0, 1).tolist()]
        yrng = [int(entry) for entry in numpy.arange(0, 1).tolist()]
        cs   = [0.0] * 2
        cnx  = 6
        cny  = 6
        step = 2
        tavg = a
        return_new = fu.ali2d_single_iter(data,numr,wr,cs,tavg,cnx,cny,xrng,yrng,step)
        return_old = oldfu.ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
        self.assertEqual(return_new,return_old)


    def test_ang_n_true_should_return_equal_object(self):
        return_new = fu.ang_n(2 , 'f' , 3)
        return_old = oldfu.ang_n(2, 'f' , 3)

        self.assertEqual(return_new, return_old)

    def test_log2_true_should_return_equal_object(self):
        return_new = fu.log2(10)
        return_old = oldfu.log2(10)

        self.assertEqual(return_new,return_old)

    def test_Numrinit_true_should_return_equal_object(self):
        return_new = fu.Numrinit(2 , 5)
        return_old = oldfu.Numrinit(2, 5)

        self.assertEqual(return_new, return_old)

    def test_ringwe_true_should_return_equal_object(self):
        numr = numpy.linspace(0, 9, 10)
        numr = list(numr)
        return_new = fu.ringwe(numr)
        return_old = oldfu.ringwe(numr)

        self.assertEqual(return_new, return_old)

    # def test_ornq(self):
    #     # will do something later
    #     a,b,c = get_data(3)
    #     image = [a,b,c]
    #     crefim = a
    #     xrng = list(numpy.arange(0, 4))
    #     yrng = list(numpy.arange(0, 4))
    #     step = 2
    #     mode = 'f'
    #     numr = numpy.arange(0, 4)
    #     numr = list(numr)




    # def test_ormq(self):
    #     # will do something later


    # def test_ormq_fast(self):
    #     # will do something later



