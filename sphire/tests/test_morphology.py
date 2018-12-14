from __future__ import print_function
from __future__ import division

import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

from ..libpy import sparx_morphology as fu
from .sparx_lib import sparx_morphology as oldfu
from ..libpy import sparx_utilities as ut


import unittest

def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


def get_data_3d(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
        data_list.append(a)

    return data_list

def get_data_gauss_noise():
    dim = 10
    return ut.model_gauss_noise(0.25 , dim,dim,dim)



class MyTestCase(unittest.TestCase):

    def test_binarize_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.binarize(image)
        return_old = oldfu.binarize(image)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_collapse_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.collapse(image)
        return_old = oldfu.collapse(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_dilation_true_should_return_equal_object(self):
        image = ut.model_blank(10,10)

        return_new = fu.dilation(image)
        return_old = oldfu.dilation(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_erosion_true_should_return_equal_object(self):
        image = ut.model_blank(10,10)

        return_new = fu.erosion(image)
        return_old = oldfu.erosion(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_power_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.power(image)
        return_old = oldfu.power(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_square_root_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.square_root(image)
        return_old = oldfu.square_root(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_square_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.square(image)
        return_old = oldfu.square(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))





if __name__ == '__main__':
    unittest.main()
