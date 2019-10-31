from __future__ import print_function
from __future__ import division

from ..libpy import sp_morphology as fu
from ..libpy import sp_morphology as oldfu

from ..libpy import sp_utilities
from ..libpy import sp_utilities as oldsparx_utilities


import unittest

from test_module import get_data, remove_dir, get_real_data
import numpy
from os import path, mkdir

from mpi import *

mpi_init(0, [])

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

"""
change it when you run the tests with your path.In this folder I copied 'TcdA1-0010_frames.mrc' got from the sphire tutorial i.e.: 'SphireDemoResults/CorrectedSums/corrsum':
"""
ABSOLUTE_PATH_TO_MRC_FILES = (
    "/home/adnan/Downloads/sphire_1_0_precalculated_results/SphireDemoResults/"
)
TOLERANCE = 0.0075

IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = sp_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sp_utilities.model_blank(10, 10, 10)
MASK = sp_utilities.model_circle(2, 5, 5)


class Test_MORPHOLOGY(unittest.TestCase):
    """
    1) Since the process finishes with an not-specified exit code, we cannot test it uniquely
    2) the nosetests are not able to run the SystemExit raise. It seems to be a known bug https://code.google.com/archive/p/python-nose/issues?page=5
    """

    """ default params got from sxcter.py and Test_defocusgett"""
    defocus = 1
    cs = 2
    voltage = 300
    pixel_size = 1.0
    bfactor = 0
    amp_contrast = 0.1  # it is the 'ac' input user param
    wn = 512
    i_start = 0
    i_stop = 10
    vpp_options = [0.3, 9.0, 0.1, 5.0, 175.0, 5.0]
    image1 = get_data(1, 256)[0]
    selection_list = "image.mrc"
    input_image_path = path.join(ABSOLUTE_PATH_TO_MRC_FILES, "TcdA1-*_frames_sum.mrc")
    output_directory = path.join(ABSOLUTE_PATH_TO_MRC_FILES, "cter_mrk_results")

    """ CTER_VPP cases """

    def test_cter_vpp_invalid_pixelSize_SYSEXIT(self):
        """
        output message: ERROR!!! Pixel size (0.000000) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            oldfu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_vpp_invalid_windowSize_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            oldfu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_vpp_no_image_input_in_path_SYSEXIT(self):
        """
        output message: ERROR!!! Input image file path (.) for All Micrographs Mode must be a  path pattern containing wild card (*). Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_vpp(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_vpp(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_vpp_output_dir_already_exist_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        mkdir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                vpp_options=self.vpp_options,
            )
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_vpp_stackModeTrue_with_invalid_input_image_path(self):
        """
        output message: ERROR!!! Stack file path specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode should not contain wild card (*). Please check input_image_path argument.
                        ERROR!!! Stack file specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode does not exist. Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                vpp_options=self.vpp_options,
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_vpp(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                vpp_options=self.vpp_options,
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    """ CTER_PAP cases """

    def test_cter_pap_invalid_pixelSize_SYSEXIT(self):
        """
        output message: ERROR!!! Pixel size (0.000000) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_pap_invalid_windowSize_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_pap_no_image_input_in_path_SYSEXIT(self):
        """
        output message: ERROR!!! Input image file path (.) for All Micrographs Mode must be a  path pattern containing wild card (*). Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_pap(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_pap(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_pap_output_dir_already_exist_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        mkdir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        self.assertEqual(cm_new.exception.code, None)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_pap_stackModeTrue_with_invalid_input_image_path(self):
        """
        output message: ERROR!!! Stack file path specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode should not contain wild card (*). Please check input_image_path argument.
                        ERROR!!! Stack file specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode does not exist. Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        self.assertEqual(cm_new.exception.code, None)
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_pap(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    """ CTER_MRK cases """

    def test_cter_mrk_invalid_pixelSize_SYSEXIT(self):
        """
        output message: ERROR!!! Pixel size (0.000000) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=0,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_mrk_invalid_windowSize_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=0,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_mrk_no_image_input_in_path_SYSEXIT(self):
        """
        output message: ERROR!!! Input image file path (.) for All Micrographs Mode must be a  path pattern containing wild card (*). Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_mrk(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_mrk(
                ".",
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_mrk_output_dir_already_exist_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        mkdir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        with self.assertRaises(SystemExit) as cm_old:
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
            )
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)

    def test_cter_mrk_stackModeTrue_with_invalid_input_image_path(self):
        """
        output message: ERROR!!! Stack file path specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode should not contain wild card (*). Please check input_image_path argument.
                        ERROR!!! Stack file specified by input_image_path (/home/lusnig/Downloads/mrc_files_for_unit_test/TcdA1-*_frames_sum.mrc) for Stack Mode does not exist. Please check input_image_path argument.
        """
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_new:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        remove_dir(self.output_directory)
        with self.assertRaises(SystemExit) as cm_old:
            mpi_barrier(MPI_COMM_WORLD)
            fu.cter_mrk(
                self.input_image_path,
                self.output_directory,
                selection_list=None,
                wn=self.wn,
                pixel_size=self.pixel_size,
                Cs=self.cs,
                voltage=self.voltage,
                f_start=self.i_start,
                f_stop=self.i_stop,
                kboot=3,
                overlap_x=50,
                overlap_y=50,
                edge_x=0,
                edge_y=0,
                check_consistency=False,
                stack_mode=True,
                debug_mode=False,
                program_name="cter_vpp() in morphology.py",
                RUNNING_UNDER_MPI=True,
                main_mpi_proc=0,
                my_mpi_proc_id=0,
                n_mpi_procs=1,
            )
        remove_dir(self.output_directory)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)


class Test_UTILITIES(unittest.TestCase):
    def test_error_img_per_group_larger_than_the_number_of_particles_exit(self):
        proj_angles = []
        for i in range(9):
            i = +0.1
            proj_angles.append([i / 2, i / 5, i / 4, i / 3, i])
        proj_angles.sort()
        proj_angles_list = numpy.full((9, 4), 0.0, dtype=numpy.float32)
        for i in range(9):
            proj_angles_list[i][0] = proj_angles[i][1]
            proj_angles_list[i][1] = proj_angles[i][2]
            proj_angles_list[i][2] = proj_angles[i][3]
            proj_angles_list[i][3] = proj_angles[i][4]

        with self.assertRaises(SystemExit) as cm_new:
            sparx_utilities.nearest_proj(proj_angles_list)
        with self.assertRaises(SystemExit) as cm_old:
            oldsparx_utilities.nearest_proj(proj_angles_list)
        self.assertEqual(cm_new.exception.code, None)
        self.assertEqual(cm_new.exception.code, cm_old.exception.code)
