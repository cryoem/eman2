from __future__ import print_function
from __future__ import division




from sphire.bin_py3 import sp_recons3d_n as oldfu
from sphire.bin import sp_recons3d_n as fu

from os import path,listdir
from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH,remove_list_of_file
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

"""
WHAT IS MISSING:
I just test the error case because the script collects the input values from the gui and then call "recons3d_n" or "recons3d_trl_MPI" form sp_applications .

"""


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => Wrong interpolation method. The current options are 4nn, and tril. 4nn is the default one.')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])
