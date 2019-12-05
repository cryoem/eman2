from __future__ import division
from past.utils import old_div
import numpy
import copy
import global_def

import unittest
import os
import shutil
from sparx.libpy import utilities as fu

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))
print(ABSOLUTE_PATH)


class MyTestCase(unittest.TestCase):
   def test_angular_distribution_returns_same_results(self):

       params_file = ABSOLUTE_PATH + "/final_params_032.txt"
       output_folder_new = "Angular_distribution_New"
       output_folder_old = "Angular_distribution_Old"
       prefix = "angdis"
       method = "P"
       pixel_size = 1.14
       delta = 3.75
       symmetry = "icos"
       box_size = 320
       particle_radius = 140
       dpi  = 72

       if os.path.isdir(output_folder_new):
           shutil.rmtree(output_folder_new)

       if os.path.isdir(output_folder_old):
           shutil.rmtree(output_folder_old)

       import time
       start = time.time()
       return_new = fu.angular_distribution(params_file, output_folder_new, prefix, method, \
                                    pixel_size, delta, symmetry, box_size,particle_radius, \
                                            dpi, do_print=True)
       print(time.time()-start)
       start = time.time()
       return_old = fu.angular_distribution_old(params_file, output_folder_old, prefix, method, \
                                    pixel_size, delta, symmetry, box_size, particle_radius,\
                                            dpi, do_print=True)
       print(time.time() - start)
       self.assertEqual(return_new, return_old)




if __name__ == '__main__':
    unittest.main()
