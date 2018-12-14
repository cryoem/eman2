from __future__ import print_function
from __future__ import division

import unittest
import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

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


class MyTestCase(unittest.TestCase):
    # def test_ali2d_MPI_true_should_return_equal_object(self):
    #     return_new = fu.ali2d_MPI(2 , 'f' , 3)
    #     return_old = oldfu.ali2d_MPI(2, 'f' , 3)




if __name__ == '__main__':
    unittest.main()
