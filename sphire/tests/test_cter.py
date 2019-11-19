from __future__ import print_function
from __future__ import division



from sphire.utils.SPHIRE.bin import sp_cter as fu
from sphire.bin import sp_cter as oldfu
from os import path,listdir
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH,remove_list_of_file,remove_dir
import unittest

"""
WHAT IS MISSING:
0) we have to run it using mpirun
1) the results, the txt and the hdf have always differences (sometimes small and sometimes big) even when I run the same
    script version ... How I can perform compability and unittest properly?


RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) Test_run see what_missing_1


In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:
1) Test_Error_cases::test_negative_radius_error 
"""



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


class Test_run(unittest.TestCase):
    old_output_folder="CterOld"
    new_output_folder = "CterNew"

    def remove_folders(self):
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)


    def test_cter_mrk(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-001*_frames_sum.mrc'"),self.new_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'"),self.old_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()
        self.remove_folders()


    def test_cter_vpp(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-001*_frames_sum.mrc'"),self.new_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--vpp", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'"),self.old_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--vpp", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()
        self.remove_folders()

    def test_cter_pap(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-001*_frames_sum.mrc'"),self.new_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--pap", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "CorrectedSums", "corrsum","TcdA1-*_frames_sum.mrc'"),self.old_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "DriftAsses","Tutorial_selected.txt'"), "--apix=1.14", "--Cs=0", "--pap", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()
        self.remove_folders()
