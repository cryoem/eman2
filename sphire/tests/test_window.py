from __future__ import print_function
from __future__ import division
import unittest
import sys
from os import path, listdir
from time import sleep
from sphire.utils.SPHIRE.bin import sp_window as fu
from sphire.bin import sp_window as oldfu
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH,remove_list_of_file

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
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Specified CTER partres file is not found. Please check input_ctf_params_source argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_input_coordinates_path_pattern_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Input coordinates file name pattern must contain wild card (*). Please check input_coordinates_pattern argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_input_micrograph_path_pattern_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),"nofile", "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),"nofile", "nofile",  "nofileCCTER",'lucaprovawindow','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_existing_output_dir_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,'--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Output directory exists. Please change the name and restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_selection_list_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--selection_list=invalid1', '--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--selection_list=invalid1','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: File specified by selection_list option does not exists. Please check selection_list option. Run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_coordinates_format_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--coordinates_format=invalid_coordinates_format', '--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--coordinates_format=invalid_coordinates_format','--box_size=352' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --coordinates_format=invalid_coordinates_format. Please run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_wrong_box_sizes_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=0' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=0' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --box_size=0. The box size must be an interger larger than zero. Please run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])


    def test_resample_ratio_higher_1_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=110' ,'--resample_ratio=2' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=110','--resample_ratio=2' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --resample_ratio=2.0. Please run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])

    def test_resample_not_higher_0_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=110' ,'--resample_ratio=0' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=110','--resample_ratio=0' ]
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
        self.assertEqual(print_new.getvalue().split('\n')[3],'** Error: Invalid option value: --resample_ratio=0.0. Please run sp_window.py -h for help.')
        self.assertEqual(print_new.getvalue().split('\n')[3],print_old.getvalue().split('\n')[3])



    def test_coordinates_file_notFound_format_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
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
        a=print_new.getvalue().split('\n')[11].split("(")[0]
        self.assertEqual(print_new.getvalue().split('\n')[11].split("(")[0],'** Error: No coordinates files are found in the directory specified by coordinates file path pattern ')
        self.assertEqual(print_new.getvalue().split('\n')[11].split("(")[0],print_old.getvalue().split('\n')[11].split("(")[0])


#testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_window.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CorrectedSums","corrsum","TcdA1-*_frames_sum.mrc"), path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Coordinates","corrsum_dose_filtered-TcdA1-*_frames_sum.mrc"),  path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"CTFEST","Tutorial_partres_select.txt"), "test_windows",'--box_size=352' ]
#--selection_list='invalid1' --coordinates_format='invalidformat'