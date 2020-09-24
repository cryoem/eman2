from __future__ import print_function
from __future__ import division

from sphire.bin_py3 import sp_proj_compare as oldfu
from ..sphire.bin import sp_proj_compare as fu
from os import path
from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir
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
from numpy import allclose

"""
WHAT IS MISSING:
0) run with 'meridien' and 'projmatch projs cases and all the helper functions called by these cases
1) 'prepare_outdir_log' is just creating the output directory for log stuff. I 'm not going to test it
2) cannot run from pycharm because it cannot find 'e2proc2d.py'

RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!


In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:

"""

'''Since it call the e2proc2d.py via comandline if you run it via pycharm it will be not able to find the file and it 
will crash.
Call it via console:
-) pytest test_proj_compare.py::Test_run
-) nostests test_proj_compare.py:Test_run
'''



class Test_run(unittest.TestCase):
    old_output_folder="compar2Dold"
    new_output_folder = "compar2Dnew"


    @classmethod
    def tearDownClass(cls):
        remove_dir(cls.new_output_folder)
        remove_dir(cls.old_output_folder)


    def test_(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         self.new_output_folder,
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        testargs_old =  [path.join(ABSOLUTE_OLDBIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         self.old_output_folder,
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        with patch.object(sys, 'argv', testargs_new):
            options = fu.parse_command_line()
            outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
            if options.mode == "viper":
                selectdoc = options.classselect
            elif options.mode == "projmatch":
                selectdoc = None
            elif options.mode == "meridien":
                selectdoc = options.partselect

            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main_proj_compare(
                options.classavgs,
                options.vol3d,
                outdir,
                options,
                mode=options.mode,
                prjmethod=options.prjmethod,
                classangles=options.classangles,
                partangles=options.partangles,
                selectdoc=selectdoc,
                displayYN=options.display,
                debug=False
            )
        with patch.object(sys, 'argv', testargs_old):
            options = fu.parse_command_line()
            outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
            if options.mode == "viper":
                selectdoc = options.classselect
            elif options.mode == "projmatch":
                selectdoc = None
            elif options.mode == "meridien":
                selectdoc = options.partselect
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main_proj_compare(
                options.classavgs,
                options.vol3d,
                outdir,
                options,
                mode=options.mode,
                prjmethod=options.prjmethod,
                classangles=options.classangles,
                partangles=options.partangles,
                selectdoc=selectdoc,
                displayYN=options.display,
                debug=False
            )
        sys.stdout = old_stdout

        averageCCC_old=float(print_old.getvalue().split('\n')[17].split(" ")[-1])
        averageCCC_new = float(print_new.getvalue().split('\n')[17].split(" ")[-1])
        self.assertTrue(allclose([averageCCC_new],[averageCCC_old],atol=0.000001))
        self.assertTrue(allclose([averageCCC_new],[0.819339885874],atol=0.000001))


class Test_helperFunctions(unittest.TestCase):
    def test_check(self):
        with self.assertRaises(SystemExit):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.check(file="file_not_found.txt", verbose=False)
        with self.assertRaises(SystemExit):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.check(file="file_not_found.txt", verbose=False)
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[0].split("ERROR")[1],"!! file_not_found.txt doesn't exist!")
        self.assertEqual(print_new.getvalue().split('\n')[0].split("ERROR")[1],print_old.getvalue().split('\n')[0].split("ERROR")[1])



class Test_Error_cases(unittest.TestCase):
    def test_error_no_input_alignment_params(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         "compar2Dluca",
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        testargs_old =  [path.join(ABSOLUTE_OLDBIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         "compar2Dlucaold",
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                options = fu.parse_command_line()
                outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main_proj_compare(
                    options.classavgs,
                    options.vol3d,
                    outdir,
                    options,
                    mode="viper",
                    prjmethod=options.prjmethod,
                    classangles=None,
                    partangles=options.partangles,
                    selectdoc=options.classselect,
                    displayYN=options.display,
                )
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                options = fu.parse_command_line()
                outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main_proj_compare(
                    options.classavgs,
                    options.vol3d,
                    outdir,
                    options,
                    mode="viper",
                    prjmethod=options.prjmethod,
                    classangles=None,
                    partangles=options.partangles,
                    selectdoc=options.classselect,
                    displayYN=options.display,
                )
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[2], 'ERROR!! Input alignment parameters not specified.')
        self.assertEqual(print_new.getvalue().split('\n')[2], print_old.getvalue().split('\n')[2])

    def test_error_volume_and_stack_have_different_dimension(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"14_MASK_FOR_CTF","sp_mask_mask.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         "compar2Dluca",
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        testargs_old =  [path.join(ABSOLUTE_OLDBIN_PATH, "sp_proj_compare.py"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"14_MASK_FOR_CTF","sp_mask_mask.hdf"),
                         path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_volume.hdf"),
                         "compar2Dlucaold",
                         "--classangles="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "run001","rotated_reduced_params.txt"),
                         "--classselect="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "07_RVIPER","main001", "this_iteration_index_keep_images.txt")
                         ]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                options = fu.parse_command_line()
                outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
                if options.mode == "viper":
                    selectdoc = options.classselect
                elif options.mode == "projmatch":
                    selectdoc = None
                elif options.mode == "meridien":
                    selectdoc = options.partselect

                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main_proj_compare(
                    options.classavgs,
                    options.vol3d,
                    outdir,
                    options,
                    mode=options.mode,
                    prjmethod=options.prjmethod,
                    classangles=options.classangles,
                    partangles=options.partangles,
                    selectdoc=selectdoc,
                    displayYN=options.display,
                )
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                options = fu.parse_command_line()
                outdir = path.dirname(path.realpath(options.classavgs)) if not options.outdir else options.outdir
                if options.mode == "viper":
                    selectdoc = options.classselect
                elif options.mode == "projmatch":
                    selectdoc = None
                elif options.mode == "meridien":
                    selectdoc = options.partselect

                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main_proj_compare(
                    options.classavgs,
                    options.vol3d,
                    outdir,
                    options,
                    mode=options.mode,
                    prjmethod=options.prjmethod,
                    classangles=options.classangles,
                    partangles=options.partangles,
                    selectdoc=selectdoc,
                    displayYN=options.display,
                )
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[2], "ERROR!! Dimension of input volume doesn't match that of image stack: 76 vs. 352")
        self.assertEqual(print_new.getvalue().split('\n')[2], print_old.getvalue().split('\n')[2])


