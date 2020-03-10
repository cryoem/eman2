from __future__ import print_function
from __future__ import division



from numpy import allclose,array_equal
from sphire.bin_py2 import sp_filterlocal as oldfu
from sphire.bin import sp_filterlocal as fu

from os import path
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_list_of_file
import unittest
from sp_utilities import get_im

try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO
import sys



class Test_run(unittest.TestCase):
    def test_(self):
        old_final="old_final.hdf"
        new_final = "new_final.hdf"
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_filterlocal.py"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_combined.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf"),
                        old_final,
                        "--radius=145 "]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_filterlocal.py"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_combined.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES", "localres.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf"),
                        new_final,
                        "--radius=145 "]

        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()


        return_old = get_im(old_final)
        return_new = get_im(new_final)
        self.assertTrue(array_equal(return_old.get_3dview(),return_new.get_3dview()))
        self.assertTrue(allclose(return_new.get_3dview().flatten().tolist()[4034427:4034527],[3.201385334250517e-05, -0.0033583296462893486, -0.0033261957578361034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],atol=0.01))
        remove_list_of_file([old_final,new_final])