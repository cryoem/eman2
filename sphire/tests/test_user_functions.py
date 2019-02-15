from __future__ import print_function
from __future__ import division

import unittest

import cPickle as pickle
import os
import sys
from mpi import *
import global_def

mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from ..libpy import sparx_user_functions as fu

from .sparx_lib import sparx_user_functions as oldfu


class Test_lib_user_functions_compare(unittest.TestCase):
    def test_amoeba_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error.ref_ali2d")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params1, ali_params2, r) = argum[0]

        return_new = fu.ref_ali2d(ali_params1, ali_params2, r)
        return_old = oldfu.ref_ali2d(ali_params1, ali_params2, r)

        self.assertEqual(return_new, return_old)



if __name__ == '__main__':
    unittest.main()
