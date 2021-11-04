from __future__ import print_function
from __future__ import division




from bin_py3 import sp_recons3d_n as oldfu
from sphire.bin import sp_recons3d_n as fu

from os import path
### Comment from Adnan......this file has to be modified because now we want to run everything on CI gitlab
ABSOLUTE_PATH_TO_RESOURCES = "resources_tests/pickles/"
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

#Absolute path to the bin folders for compatibility tests purpose.  since in the bin folder we have to set the parameter via cmd I have to run them.
ABSOLUTE_SPHIRE_PATH= ABSOLUTE_PATH.split("/tests")[0]
ABSOLUTE_OLDBIN_PATH = path.join(ABSOLUTE_SPHIRE_PATH, "bin_py3")
ABSOLUTE_BIN_PATH = path.join(ABSOLUTE_SPHIRE_PATH,"bin")
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

#todo: why the new version does not work????
"""
WHAT IS MISSING:
I just test the error case because the script collects the input values from the gui and then call "recons3d_n" or "recons3d_trl_MPI" form sp_applications .


WHAT IS WRONG AFTER MOVING TO PYTHON3.
The skipped tests are failing because an MPI error.
The errors happen in the "sphire.bin" folder not in  "bin_py3"

    error message:
    test_recons3d_n.py::Test_Error_cases::test_error_group_and_list_together
        *** The MPI_Finalize() function was called after MPI_FINALIZE was invoked.
        *** This is disallowed by the MPI standard.
        *** Your MPI job will now abort.
        [rtxr2:124665] Local abort after MPI_FINALIZE started completed successfully, but am not able to aggregate error messages, and not able to guarantee that all other processes were killed!
        No data to report.
"""

class Test_Error_cases(unittest.TestCase):
    @unittest.skip("MPI error message ONLY in the sphire.bin version")
    def test_too_few_input_values(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_recons3d_n.py")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_recons3d_n.py")]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => Incomplete list of arguments')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])

    @unittest.skip("MPI error message ONLY in the sphire.bin version")
    def test_error_group_and_list_together(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_recons3d_n.py"), "input_stack","out_volume", "--list='list'", '1','11','1', "--group=1"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_recons3d_n.py"), "input_stack","out_volume", "--list='list'", '1','11','1', "--group=1"]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => options group and list cannot be used together')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])

    def test_error_interpolation_method(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_recons3d_n.py"), "input_stack","out_volume", "--list='list'", '1','11','1', "--interpolation_method='invalid_method'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_recons3d_n.py"), "input_stack","out_volume", "--list='list'", '1','11','1', "--interpolation_method='invalid_method'"]
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
        self.assertEqual(print_new.getvalue().split('\n')[4].split("ERROR")[1],' => Wrong interpolation method. The current options are 4nn, and tril. 4nn is the default one.')
        self.assertEqual(print_old.getvalue().split('\n')[1].split("ERROR")[1],' => Wrong interpolation method. The current options are 4nn, and tril. 4nn is the default one.')
        self.assertEqual(print_new.getvalue().split('\n')[4].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])
