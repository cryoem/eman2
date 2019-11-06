from __future__ import print_function
from __future__ import division



from sphire.utils.SPHIRE.bin import sp_cter as fu
from sphire.bin import sp_cter as oldfu
from os import path,listdir
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH,remove_list_of_file
import unittest




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


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        remove_list_of_file([f for f in listdir(".") if "sp_cter_logfile" in f])

    def test_lowest_resolution_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc'"), "'/home/lusnig/Downloads/luca5nov'", '--apix=1.0', '--f_start=0.3' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'"), "'/home/lusnig/Downloads/luca5nov'", '--apix=1.0','--f_start=0.3']
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => f_start should be in Angstrom')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])

    def test_highest_resolution_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc'"), "'/home/lusnig/Downloads/luca5nov'", '--apix=1.0', '--f_stop=0.3']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'"), "'/home/lusnig/Downloads/luca5nov'", '--apix=1.0', '--f_stop=0.3']
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => f_stop should be in Angstrom')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])

    def test_too_few_params_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc'")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'")]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit) as cnew:
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old) :
            with self.assertRaises(SystemExit) as cold:
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout

        self.assertEqual(str(cnew.exception),str(cold.exception))
        self.assertEqual(str(cnew.exception),"None")

