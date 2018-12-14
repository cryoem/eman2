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
    # def test_ali2d_single_iter_true_should_return_equal_objects(self):
    #     a,b,c = get_data(3)
    #     data = [a]
    #     numr = [int(entry) for entry in numpy.arange(0, 6).tolist()]
    #     wr   = [0.0] * int(len(numr)/3)    #list(numpy.arange(0, 4))
    #     xrng = [int(entry) for entry in numpy.arange(0, 1).tolist()]
    #     yrng = [int(entry) for entry in numpy.arange(0, 1).tolist()]
    #     cs   = [0.0] * 2
    #     cnx  = 2
    #     cny  = 2
    #     step = 16
    #     tavg = e2cpp.EMData(10, 10)
    #     tavg = numpy.arange(10 * 10, dtype=numpy.float32).reshape(10, 10)
    #     return_new = fu.ali2d_single_iter(data,numr,wr,cs,tavg,cnx,cny,xrng,yrng,step)
    #     print("hello ")


        # return_old = oldfu.ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step)
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

    # def test_ornq(self):
    #     # will do something later
    #     a,b = get_data(2)
    #     image  = a
    #     crefim = b
    #     xrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #     yrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #     step   = 2
    #     mode   = 'f'
    #     numr = [int(entry) for entry in numpy.arange(0, 6).tolist()]
    #     cnx  = 2
    #     cny  = 2
    #     return_new = fu.ornq(a,b,xrng,yrng,step,mode,numr,cnx,cny)
    #     return_old = fu.ornq(a, b, xrng, yrng, step, mode, numr, cnx, cny)
    #
    #     self.assertEqual(return_new, return_old)
    #
    # def test_ormq(self):
    #     a,b = get_data(2)
    #     image  = a
    #     crefim = b
    #     xrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #     yrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #     step   = 2
    #     mode   = 'f'
    #     numr = [int(entry) for entry in numpy.arange(0, 6).tolist()]
    #     cnx  = 2
    #     cny  = 2
    #     return_new = fu.ormq(a,b,xrng,yrng,step,mode,numr,cnx,cny)
    #     return_old = fu.ormq(a, b, xrng, yrng, step, mode, numr, cnx, cny)
    #
    #     self.assertEqual(return_new, return_old)

    # def test_ormq_fast(self):
    #     a,b,c,d,e,f,g,h,i = get_data(9)
    #     image = []
    #     image.append(a)
    #     image.append(b)
    #     image.append(c)
    #     image.append(d)
    #     image.append(e)
    #     image.append(f)
    #     image.append(g)
    #     image.append(h)
    #     crefim = i
    #     # xrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #     # yrng = [int(entry) for entry in numpy.arange(0, 4).tolist()]
    #
    #     xrng = 4
    #     yrng = 4
    #     step = 1
    #     numr = [int(entry) for entry in numpy.arange(0, 6).tolist()]
    #     mode = 'f'
    #
    #     return_new = fu.ormq_fast(image, crefim, xrng, yrng, step, numr, mode)
    #     return_old = fu.ormq_fast(image, crefim, xrng, yrng, step, numr, mode)
    #
    #     self.assertEqual(return_new, return_old)


    # def test_prepref(self):
    #     a,b,c = get_data(3)
    #     data = [a,b,c]
    #     maskfile = ????
    #     mode   = 'f'
    #     numr = [int(entry) for entry in numpy.arange(0, 6).tolist()]
    #     cnx  = 2
    #     cny  = 2
    #     maxrangex = 5
    #     maxrangey = 5
    #     step = 3
    #
    #
    #     return_new = fu.prepref(data,maskfile,cnx,cny,numr,mode,maxrangex,maxrangey,step)
    #     return_old = fu.prepref(data, maskfile, cnx, cny, numr, mode, maxrangex, maxrangey, step)





