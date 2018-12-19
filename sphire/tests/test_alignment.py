from __future__ import print_function
from __future__ import division

import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

import unittest

import os
import pickle
ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))

from ..libpy import sparx_alignment as fu
from .sparx_lib import sparx_alignment as oldfu
from ..libpy import sparx_utilities as ut


def get_data(num,dim = 10):
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


def get_data_3d(num, dim=10):
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
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
        filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/alignment.ornq")
        with open(filepath, 'rb') as rb:
            (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi) = pickle.load(rb)

        tavg, = get_data(1,352)
        cs = 2
        wr = 512
        # cs = [0.0] * 2
        # wr = [0.0] * int(len(numr) / 3)

        print(numpy.shape(image.get_3dview()))
        print(numpy.shape(tavg.get_3dview()))
        print(len(numr))

        return_new = fu.ali2d_single_iter(image, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step,)
        # return_old = oldfu.ali2d_single_iter(image, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, mode = mode)
        # self.assertEqual(return_new,return_old)


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

    def test_ornq(self):

        filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/alignment.ornq")
        with open(filepath, 'rb') as rb:
            (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi) = pickle.load(rb)

        return_new = fu.ornq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi)
        return_old = fu.ornq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi)
        self.assertEqual(return_new, return_old)

    def test_ormq(self):
        filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/alignment.ornq")
        with open(filepath, 'rb') as rb:
            (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi) = pickle.load(rb)

        return_new = fu.ormq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi)
        return_old = fu.ormq(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi)

        self.assertEqual(return_new, return_old)

    # def test_ormq_fast(self):

    #     filepath = os.path.join(ABSOLUTE_PATH, "files/picklefiles/alignment.ornq")
    #     with open(filepath, 'rb') as rb:
    #         (image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi) = pickle.load(rb)
    #
    #     return_new = fu.ormq_fast(image, crefim, xrng, yrng, step, numr, mode)
        # return_old = fu.ormq_fast(image, crefim, xrng, yrng, step, numr, mode)
        #
        # self.assertEqual(return_new, return_old)


    # def test_prepref(self):
    # dont know what should be maskfile
    #     a,b,c = get_data(3)
    #     data = [a,b,c]
    #     maskfile = b
    #     mode   = 'f'
    #     numr = [int(entry) for entry in numpy.arange(0, 20).tolist()]
    #     cnx  = 2
    #     cny  = 2
    #     maxrangex = 5
    #     maxrangey = 5
    #     step = 3
    #
    #
    #     return_new = fu.prepref(data,maskfile,cnx,cny,numr,mode,maxrangex,maxrangey,step)
    #     return_old = fu.prepref(data,maskfile,cnx,cny,numr,mode,maxrangex,maxrangey,step)
    #
    #     self.assertEqual(return_new, return_old)




