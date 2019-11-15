from __future__ import print_function
from __future__ import division




from sphire.bin import sp_viper as oldfu
from sphire.utils.SPHIRE.bin import sp_viper as fu
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

"""
WHAT IS MISSING:
I just test the error case because the script collects the input values from the gui and then call "multi_shc" from sp_multi_shc.

NB:
If for some reason you want to run the test you have to use mpirun with --nrun option activated

"""


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):
    def test_too_few_input_values(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_viper.py")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_viper.py")]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[7].split("ERROR")[1],' => Invalid number of parameters used. Please see usage information above.')
        self.assertEqual(print_new.getvalue().split('\n')[7].split("ERROR")[1],print_old.getvalue().split('\n')[7].split("ERROR")[1])

    def test_invalid_nruns(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_viper.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_viper.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov"]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])