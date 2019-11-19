from __future__ import print_function
from __future__ import division




from sphire.bin import sp_mask as oldfu
from sphire.utils.SPHIRE.bin import sp_mask as fu
from os import path
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH
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
1)  I just test the error case because the script collects the input values from the gui and then call
        "sp_morphology.adaptive_mask_scipy", sparx_filter.filt_tanland some basic functions from "sp_utilities"
2) cannot catch the error message in "test_too_few_input_values". no clue why it is happening

NB:
The error message are addressed on the standard error channel  

"""


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):

    def test_too_few_input_values(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_mask.py")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_mask.py")]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stderr
                print_new = StringIO()
                sys.stderr = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stderr = print_old
                oldfu.main()
        sys.stderr = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[17],'sp_mask.py: error: too few arguments')
        self.assertEqual(print_new.getvalue().split('\n')[17],print_old.getvalue().split('\n')[17])

    def test_mol_mass_invalid(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_mask.py"),"DS","DS","--mol_mass=-1.0"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_mask.py"),"DS","DS","--mol_mass=-1.0"]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stderr
                print_new = StringIO()
                sys.stderr = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stderr = print_old
                oldfu.main()
        sys.stderr = old_stdout
        a=print_new.getvalue().split('\n')[17]

        self.assertEqual(print_new.getvalue().split('\n')[17],'sp_mask.py: error: Minimum value for --mol_mass is 0')
        self.assertEqual(print_new.getvalue().split('\n')[17],print_old.getvalue().split('\n')[17])


    def test_invalid_shape(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_mask.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "localews.hdf") ,"output_folderMASK","--mol_mass=1.0",
                        "--second_mask='/home/lusnig/Downloads/SphireDemoResults/Sharpening-after-Meridien-Substack-Local_000/vol_adaptive_mask.hdf","--second_mask_shape='sphere'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_mask.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "localews.hdf") ,"output_folderMASK","--mol_mass=1.0", "--second_mask_shape='sphere'"]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stderr
                print_new = StringIO()
                sys.stderr = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stderr = print_old
                oldfu.main()
        sys.stderr = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[17].replace('"',''),"sp_mask.py: error: argument --second_mask_shape/--sms: invalid choice: 'sphere' (choose from 'cylinder', 'sphere', 'cube')")
        self.assertEqual(print_new.getvalue().split('\n')[17],print_old.getvalue().split('\n')[17])