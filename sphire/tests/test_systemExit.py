from __future__ import print_function
from __future__ import division

from ..libpy import sparx_morphology as fu

from ..libpy import sparx_utilities


import unittest

from test_module import get_data, remove_dir,get_real_data

from os import path, mkdir



ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

"""
change it when you run the tests with your path.In this folder I copied 'TcdA1-0010_frames.mrc' got from the sphire tutorial i.e.: 'SphireDemoResults/CorrectedSums/corrsum':
"""
ABSOLUTE_PATH_TO_MRC_FILES="/home/lusnig/Downloads/mrc_files_for_unit_test"
TOLERANCE = 0.0075

IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = sparx_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sparx_utilities.model_blank(10, 10, 10)
MASK = sparx_utilities.model_circle(2, 5, 5)


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
    amp_contrast = 0.1   # it is the 'ac' input user param
    wn = 512
    i_start = 0
    i_stop = 10
    vpp_options = [0.3, 9.0, 0.1, 5.0, 175.0, 5.0]
    image1 = get_data(1, 256)[0]
    selection_list = 'image.mrc'
    input_image_path = path.join(ABSOLUTE_PATH_TO_MRC_FILES, "TcdA1-*_frames_sum.mrc")
    output_directory = path.join(ABSOLUTE_PATH_TO_MRC_FILES, "cter_mrk_results")


    def test_cter_vpp_invalid_pixelSize_SYSEXIT(self):
        """
        output message: ERROR!!! Pixel size (0.000000) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option.
        """
        with self.assertRaises(SystemExit) as cm:
            fu.cter_vpp(self.input_image_path, self.output_directory, selection_list=None, wn=self.wn,
                        pixel_size=0, Cs=self.cs, voltage=self.voltage, f_start=self.i_start,
                        f_stop=self.i_stop, vpp_options=self.vpp_options)
        self.assertEqual(cm.exception.code, None)

    def test_cter_vpp_invalid_windowSize_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        with self.assertRaises(SystemExit) as cm:
            fu.cter_vpp(self.input_image_path, self.output_directory, selection_list=None, wn=0,
                        pixel_size=self.pixel_size, Cs=self.cs, voltage=self.voltage, f_start=self.i_start,
                        f_stop=self.i_stop, vpp_options=self.vpp_options)
        self.assertEqual(cm.exception.code, None)

    def test_cter_vpp_no_image_input_in_path_SYSEXIT(self):
        """
        output message: ERROR!!! Input image file path (.) for All Micrographs Mode must be a  path pattern containing wild card (*). Please check input_image_path argument.
        """
        with self.assertRaises(SystemExit) as cm:
            fu.cter_vpp(".", self.output_directory, selection_list=None, wn=self.wn,
                        pixel_size=self.pixel_size, Cs=self.cs, voltage=self.voltage, f_start=self.i_start,
                        f_stop=self.i_stop, vpp_options=self.vpp_options)
        self.assertEqual(cm.exception.code,None)

    def test_cter_vpp_output_dir_already_exist_SYSEXIT(self):
        """
        output message: ERROR!!! Output directory (/home/lusnig/Downloads/mrc_files_for_unit_test/cter_mrk_results) exists already. Please check output_directory argument.
        """
        remove_dir(self.output_directory)
        mkdir(self.output_directory)
        with self.assertRaises(SystemExit) as cm:
            fu.cter_vpp(self.input_image_path, self.output_directory, selection_list=None, wn=self.wn,
                        pixel_size=self.pixel_size, Cs=self.cs, voltage=self.voltage, f_start=self.i_start,
                        f_stop=self.i_stop, vpp_options=self.vpp_options)
        self.assertEqual(cm.exception.code, None)
        remove_dir(self.output_directory)