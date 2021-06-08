# from __future__ import print_function
# from __future__ import division
#
#
#
# from numpy import array_equal
#
# from sphire.bin_py3 import sp_process as oldfu
# from ..sphire.bin import sp_process as fu
#
# from os import path
# from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH
# import unittest
# from sp_utilities import get_im
#
# try:
#     # python 3.4+ should use builtin unittest.mock not mock package
#     from unittest.mock import patch
# except ImportError:
#     from mock import patch
#
# try:
#     from StringIO import StringIO  # python2 case
# except ImportError:
#     # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
#     from io import StringIO
# import sys
#
#
#
# """
# WHAT IS MISSING:
# All the helper function tests for the following reasons:
# -) pca, tsp because are not used ( NB: tsp will crash because a lot of function not declared)
# -) Distance, TotalDistance,reverse,transpt because are used in tsp
#
# RESULT AND KNOWN ISSUES
#
#
# """
#
# class Test_run(unittest.TestCase):
#     old_output_folder="combineMapsOld"
#     new_output_folder = "combineMapsNew"
#     fname = "vol_combined.hdf"
#
#
#     # @classmethod
#     # def tearDownClass(cls):
#     #     remove_dir(cls.new_output_folder)
#     #     remove_dir(cls.old_output_folder)
#
#     # it is the run of the tutorial (pg.65)
#     # At pg.93 there is another run. I did not test it because it performs the same operation done in this test
#     # At pg.98 there is another run. I did not test it because it performs the same operation done in this test
#     def test_combinemaps_Halfset(self):
#         testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_process.py"),
#                         "--combinemaps",
#                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_0_unfil_028.hdf"),
#                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_1_unfil_028.hdf"),
#                         "--output_dir="+ self.new_output_folder,
#                         "--pixel_size=1.14",
#                         "--do_adaptive_mask",
#                         "--threshold=0.02",
#                         "--mtf="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"FalconIImtf.txt")]
#         testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_process.py"),
#                         "--combinemaps",
#                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_0_unfil_028.hdf"),
#                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_1_unfil_028.hdf"),
#                         "--output_dir="+ self.old_output_folder,
#                         "--pixel_size=1.14",
#                         "--do_adaptive_mask",
#                         "--threshold=0.02",
#                         "--mtf="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"FalconIImtf.txt")]
#
#         with patch.object(sys, 'argv', testargs_new):
#             fu.main()
#         with patch.object(sys, 'argv', testargs_old):
#             oldfu.main()
#
#         return_new = get_im(path.join(self.new_output_folder, self.fname))
#         return_old = get_im(path.join(self.old_output_folder, self.fname))
#         self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
#         # self.assertTrue(allclose(return_new.get_3dview().flatten().tolist()[3641969:3642050],[8.926989539759234e-05, -9.56452640821226e-05, -0.00022921698109712452, -2.0858768039033748e-05, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0]))
#
#
