from __future__ import print_function
from __future__ import division



from numpy import array_equal, allclose
from sphire.bin import sp_mask as oldfu
from sphire.utils.SPHIRE.bin import sp_mask as fu
from os import path
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir
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

        self.assertEqual(print_new.getvalue().split('\n')[17],'sp_mask.py: error: Minimum value for --mol_mass is 0')
        self.assertEqual(print_new.getvalue().split('\n')[17],print_old.getvalue().split('\n')[17])


    def test_invalid_shape(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_mask.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres.hdf") ,"output_folderMASK","--mol_mass=1.0",
                        "--second_mask="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres_ang.hdf"),"--second_mask_shape='sphere'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_mask.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres.hdf") ,"output_folderMASK","--mol_mass=1.0",
                        "--second_mask="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres_ang.hdf"),"--second_mask_shape='sphere'"]

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



class Test_run(unittest.TestCase):
    old_output_folder="maskOld"
    new_output_folder = "maskNew"
    filename = "sp_mask_mask.hdf"
    hdf_ordered_class_averages = "ordered_class_averages.hdf"


    @classmethod
    def tearDownClass(cls):
        remove_dir(cls.new_output_folder)
        remove_dir(cls.old_output_folder)


    def test_(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_mask.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"09_ADJUSTMENT","vol3d_ref_moon_eliminated.hdf"),
                         "--pixel_size=1.14",
                         self.new_output_folder,
                         "--mol_mass=1400",
                         "--edge_width=10"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_mask.py"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "09_ADJUSTMENT",
                                  "vol3d_ref_moon_eliminated.hdf"),
                        "--pixel_size=1.14",
                        self.old_output_folder,
                        "--mol_mass=1400",
                        "--edge_width=10"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()

        old_value =get_im(path.join(self.old_output_folder,self.filename))
        new_value = get_im(path.join(self.new_output_folder,self.filename))
        self.assertTrue(array_equal(old_value.get_3dview(),new_value.get_3dview()))
        self.assertTrue(allclose( new_value.get_3dview().flatten().tolist()[3022799:3022899], [0.00024921807926148176, 0.0005635463166981936, 0.0005635463166981936, 0.00024921807926148176, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
