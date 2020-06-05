from __future__ import print_function
from __future__ import division
import unittest
import sys
from os import path, listdir
from time import sleep

from sphire.bin import sp_window as fu
from sphire.bin_py3 import sp_window as oldfu

from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_list_of_file,remove_dir
from numpy import array_equal
from mrcfile import open as mrcfile_open
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

import mpi
import sp_global_def

mpi.mpi_init(0, [])
sp_global_def.BATCH = True
sp_global_def.MPI = False

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO

"""
WHAT IS MISSING:
0) get_time_stamp_suffix --> returns a timestamp ... not testable
1) cannot test without valid input data. The tutorial 1.3 pg 26 are not ok??

RESULT AND KNOWN ISSUES
1) 

In these tests there is a bug --> syntax error:
1)

In these tests there is a strange behavior:
1) 
"""

class Test_helperFunctions(unittest.TestCase):
    def test_is_float_True(self):
        self.assertTrue(fu.is_float(3))
        self.assertTrue(oldfu.is_float(3))

    def test_is_float_False(self):
        self.assertFalse(fu.is_float("d"))
        self.assertFalse(oldfu.is_float("d"))

    def test_get_cmd_line(self):
        cmdLine=["this", "is", "a", "test"]
        with patch.object(sys, 'argv', cmdLine):
            return_new = fu.get_cmd_line()
            return_old = oldfu.get_cmd_line()
            self.assertEqual(return_new,return_old)
            self.assertEqual(return_new, 'Shell line command: this  is  a  test  ')

    def test_estimate_angle(self):
        return_new = fu.estimate_angle(coords_a=[15,32], coords_b=[18,38])
        return_old = oldfu.estimate_angle(coords_a=[15,32], coords_b=[18,38])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, -63.43494882292201)

    def test_estimate_angle_same_coords(self):
        return_new = fu.estimate_angle(coords_a=[15,32], coords_b=[15,32])
        return_old = oldfu.estimate_angle(coords_a=[15,32], coords_b=[15,32])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, -0.0)

    def test_estimate_angle_index_error_coords_a(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_angle(coords_a=[], coords_b=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_angle(coords_a=[], coords_b=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_estimate_angle_index_error_coords_b(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_angle(coords_b=[], coords_a=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_angle(coords_b=[], coords_a=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_estimate_distance(self):
        return_new = fu.estimate_distance(coords_a=[15,32], coords_b=[18,38])
        return_old = oldfu.estimate_distance(coords_a=[15,32], coords_b=[18,38])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 6.7082039324993694)

    def test_estimate_distance_same_coords(self):
        return_new = fu.estimate_distance(coords_a=[15,32], coords_b=[15,32])
        return_old = oldfu.estimate_distance(coords_a=[15,32], coords_b=[15,32])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, -0.0)

    def test_estimate_distance_index_error_coords_a(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_distance(coords_a=[], coords_b=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_distance(coords_a=[], coords_b=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_estimate_distance_index_error_coords_b(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_distance(coords_b=[], coords_a=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_distance(coords_b=[], coords_a=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        sleep(0.1)      #to give the script the time to close the last logfile.
        remove_list_of_file([f for f in listdir(".") if "sp_window_logfile" in f])

    def test_wrong_CTF_param_source_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","CBOX","TcdA1-*_frames.cbox"),  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","CBOX","TcdA1-*_frames.cbox"),  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
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

        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Specified CTER partres file is not found. Please check input_ctf_params_source argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[2],print_old.getvalue().split('\n')[2])

    def test_wrong_input_coordinates_path_pattern_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Input coordinates file name pattern must contain wild card (*). Please check input_coordinates_pattern argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[1].split('\n')[0].split(' ')[9] , print_old.getvalue().split('\n')[1].split('\n')[0].split(' ')[9] )

    def test_wrong_input_micrograph_path_pattern_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),"nofile", "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),"nofile", "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9]    ,print_old.getvalue().split('\n')[1].split(' ')[9])

    def test_existing_output_dir_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,'--box_size=352' ]
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

        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[8],
                         print_old.getvalue().split('\n')[1].split(' ')[8])
        # # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Output directory exists. Please change the name and restart the program.')
        # self.assertEqual(print_new.getvalue().split('\n')[3].split(' ')[7],print_old.getvalue().split('\n')[3].split(' ')[7])

    def test_wrong_selection_list_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--selection_list=invalid1', '--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--selection_list=invalid1','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[6],
                         print_old.getvalue().split('\n')[1].split(' ')[6])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: File specified by selection_list option does not exists. Please check selection_list option. Run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_coordinates_format_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--coordinates_format=invalid_coordinates_format', '--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--coordinates_format=invalid_coordinates_format','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --coordinates_format=invalid_coordinates_format. Please run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_box_sizes_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=0' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=0' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --box_size=0. The box size must be an interger larger than zero. Please run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])


    def test_resample_ratio_higher_1_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=110' ,'--resample_ratio=2' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=110','--resample_ratio=2' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --resample_ratio=2.0. Please run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_resample_not_higher_0_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=110' ,'--resample_ratio=0' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.cbox"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=110','--resample_ratio=0' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --resample_ratio=0.0. Please run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_no_micrograph_files_are_found(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_not_a_file.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_not_a_file.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[6].split('(')[0],'** Error: No micrograph files are found in the directory specified by micrograph path pattern ')
        # self.assertEqual(print_new.getvalue().split('\n')[6].split("(")[0],print_old.getvalue().split('\n')[6].split("(")[0])

    def test_coordinates_file_notFound_format_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","CBOX","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","CBOX","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[11].split("(")[0],'** Error: No coordinates files are found in the directory specified by coordinates file path pattern ')
        # self.assertEqual(print_new.getvalue().split('\n')[11].split("(")[0],print_old.getvalue().split('\n')[11].split("(")[0])

    def test_partres_file_not_found(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select_not_found.txt"), "test_windows_new",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select_not_found.txt"), "test_windows_old",'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Specified CTER partres file is not found. Please check input_ctf_params_source argument. Run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_error_invalid_option_value(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_new","--coordinates_format='eman1'",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_old","--coordinates_format='eman1'",'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],"** Error: Invalid option value: --coordinates_format='eman1'. Please run sp_window.py -h for help.")
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_error_invalid_resample_ratio(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_new","--resample_ratio=3",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_old","--resample_ratio=3",'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split(' ')[9],
                         print_old.getvalue().split('\n')[1].split(' ')[9])
        # self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --resample_ratio=3.0. Please run sp_window.py -h for help.')
        # self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_error_invalid_box_size(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_new","--box_size=3.5",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), "test_windows_old","--box_size=3.5",'--box_size=352' ]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stderr
                print_new = StringIO()
                sys.stderr = print_new
                fu.main()
            with patch.object(sys, 'argv', testargs_old):
                with self.assertRaises(SystemExit):
                    print_old = StringIO()
                    sys.stderr= print_old
                    oldfu.main()
        sys.stderr = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[42],"sp_window.py: error: option --box_size: invalid integer value: '3.5'")
        self.assertEqual(print_new.getvalue().split('\n')[42],print_old.getvalue().split('\n')[42])

    def test_error_files_do_not_match(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "02_CRYOLO", "EMAN","TcdA1-001*_frames.box"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_partres_select.txt"), "outfolder_new", '--box_size=352']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "02_CRYOLO", "EMAN","TcdA1-001*_frames.box"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_partres_select.txt"), "outfolder_old", '--box_size=352']
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
            with patch.object(sys, 'argv', testargs_old):
                with self.assertRaises(SystemExit):
                    print_old = StringIO()
                    sys.stdout= print_old
                    oldfu.main()
        sys.stdout = old_stdout

        out_parser=print_new.getvalue().split('\n')[16]
        err_new=out_parser[0] + out_parser[1]
        out_parser2=print_old.getvalue().split('\n')[16]
        err_old=out_parser2[0]+out_parser2[1]
        # self.assertEqual(err_new,'** Error: A micrograph name TcdA1-0100_frames.mrc) in the CTER partres file  does not match with input micrograph basename pattern TcdA1-001*_frames.mrc) ')
        self.assertEqual(err_new,err_old)



class Test_run(unittest.TestCase):
    old_output_folder="WindowOld"
    new_output_folder = "WindowNew"
    filename = "TcdA1-0010_frames_ptcls.mrcs"

    @classmethod
    def tearDownClass(cls):
        sleep(0.1)      #to give the script the time to close the last logfile.
        remove_list_of_file([f for f in listdir(".") if "sp_window_logfile" in f])

    def test_(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), self.new_output_folder,'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-*_frames.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"02_CRYOLO","EMAN","TcdA1-*_frames.box"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"01_CTER","Tutorial_partres_select.txt"), self.old_output_folder,'--box_size=352' ]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()

        old_data=mrcfile_open(path.join(self.old_output_folder,self.filename)).data[0]
        new_data = mrcfile_open(path.join(self.new_output_folder, self.filename)).data[0]
        self.assertTrue(array_equal(old_data,new_data))
        self.assertTrue(array_equal(new_data.flatten().tolist()[0:100],[2.033139944076538, 2.118870258331299, 2.0796854496002197, 1.8459678888320923, 1.7696151733398438, 1.4505984783172607, 1.2712483406066895, 1.221462607383728, 0.6416369080543518, 0.5872614979743958, -0.6979288458824158, -1.0111565589904785, -0.7067759037017822, -0.8540730476379395, -1.3486729860305786, -1.7476952075958252, -1.3114615678787231, -0.8387683629989624, -0.787508487701416, -0.9037533402442932, -0.663044810295105, -0.294295996427536, -0.09013121575117111, 0.04511815309524536, 0.19695419073104858, -0.033459026366472244, -0.3311372697353363, -0.8275147676467896, -1.0807580947875977, -1.0741568803787231, -0.533863365650177, -0.19527646899223328, -0.36622801423072815, -0.4379308223724365, -0.4848000705242157, -0.49066489934921265, -0.4228202700614929, 0.41493555903434753, 0.506386935710907, 0.2449536770582199, 0.3998281955718994, 0.2207028865814209, 0.40251466631889343, 0.3158993124961853, 0.4267796576023102, 1.6939674615859985, 1.8226191997528076, 1.889028549194336, 1.858549952507019, 1.552672028541565, 1.4727593660354614, 0.594339907169342, 0.2293643206357956, 0.16122445464134216, -0.050342343747615814, -0.2755216360092163, 0.3465840220451355, 0.2246960550546646, 0.16776084899902344, 0.22618092596530914, -0.4091913104057312, -0.3639237880706787, -0.6263998746871948, -0.37294498085975647, -0.7555119395256042, -0.7194178104400635, -0.8766475319862366, -0.7931090593338013, -0.9084152579307556, -1.1437935829162598, -1.3992165327072144, -1.2049322128295898, -1.1964422464370728, -0.5358992218971252, -0.19423657655715942, -0.054526664316654205, -0.26979440450668335, -1.20499849319458, -1.6301522254943848, -1.9685299396514893, -1.4616751670837402, -0.6779497861862183, -0.34084251523017883, 0.30665281414985657, 0.7742570638656616, 0.8332902789115906, 1.022173523902893, 1.3131401538848877, 1.2987929582595825, 1.7227187156677246, 1.8645423650741577, 1.0463528633117676, -0.2233831286430359, -0.5092124938964844, -0.24659070372581482, 0.08108771592378616, -0.021659117192029953, -0.19913935661315918, -0.7791047692298889, -0.6084468960762024]))
        remove_dir(self.old_output_folder)
        remove_dir(self.new_output_folder)