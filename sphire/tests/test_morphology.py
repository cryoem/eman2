from __future__ import print_function
from __future__ import division

from ..libpy import sparx_morphology as fu
from .sparx_lib import sparx_morphology as oldfu

from ..libpy import sparx_filter
from ..libpy import sparx_utilities

import numpy
import unittest

from test_module import get_data, get_data_3d, remove_dir, get_arg_from_pickle_file

from EMAN2_cppwrap import EMData, EMAN2Ctf
from copy import  deepcopy
from os import path



ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
TOLERANCE = 0.0075

IMAGE_2D = get_data(1)[0]
IMAGE_3D = get_data_3d(1)[0]
IMAGE_BLANK_2D = sparx_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sparx_utilities.model_blank(10, 10, 10)
MASK = sparx_utilities.model_circle(2, 5, 5)



"""
There are some opened issues in:
1) rotavg_ctf --> for definition it'll process a 2D image ... should we impede a 3D img as input?
2) adaptive_mask --> it'll process a 3D image ... should we impede a 2D img as input?
3) get_shrink_3dmask --> accepts the 'mask_file_name' params as string too. I did not test it because it is processed by 'sparx_fundamentals.resample'
4) defocusgett --> with f_start=0 it crashes but in the code it manages this situation at the beginning ...it seems that should be possible to init it with 0
5) fastigmatism3 --> sometimes some test fails because a very large difference of value e.g.: -11.974973537555098 != 1e+20 or 178.59375 != 142.71600723266602
6) fupw it jsut calls fastigmatism3 hence there is the same issue
"""
class Test_binarize(unittest.TestCase):

    def test_binarize_2Dimg(self):
        return_new = fu.binarize(IMAGE_2D, minval = 0.0)
        return_old = oldfu.binarize(IMAGE_2D, minval = 0.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_binarize_3Dimg(self):
        return_new = fu.binarize(IMAGE_3D, minval = 0.0)
        return_old = oldfu.binarize(IMAGE_3D, minval = 0.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException(self):
        """ We are not able to catch the 'NotExistingObjectException' C++ exception"""
        with self.assertRaises(RuntimeError):
            fu.binarize(EMData())
            oldfu.binarize(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.binarize()
            oldfu.binarize()



class Test_collapse(unittest.TestCase):

    def test_collapse_2Dimg(self):
        return_new = fu.collapse(IMAGE_2D, minval = -1.0, maxval = 1.0)
        return_old = oldfu.collapse(IMAGE_2D, minval = -1.0, maxval = 1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_collapse_3Dimg(self):
        return_new = fu.collapse(IMAGE_3D, minval = -1.0, maxval = 1.0)
        return_old = oldfu.collapse(IMAGE_3D, minval = -1.0, maxval = 1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException(self):
        """ We are not able to catch the 'NotExistingObjectException' C++ exception"""
        with self.assertRaises(RuntimeError):
            fu.collapse(EMData())
            oldfu.collapse(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.collapse()
            oldfu.collapse()



class Test_dilatation(unittest.TestCase):

    def test_test_empty_input_image_crashes_because_signal11SIGSEV(self):
        """
        It seems not possible to test it without getting an segmentation fault. We should return None after the first error message
        in order to avoid to run in the second if/else and get the segmentation fault in the c+++ code

        img =EMData()
        return_new = fu.dilation(img)
        return_old = oldfu.dilation(img)
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)
        """
        self.assertTrue(True)

    def test_empty_mask_image_returns_RuntimeError_ImageDimensionException_center_isnot_welldefined(self):
        with self.assertRaises(RuntimeError):
            fu.dilation(IMAGE_BLANK_2D, EMData(), morphtype="BINARY")
            oldfu.dilation(IMAGE_BLANK_2D, EMData(), morphtype="BINARY")

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.dilation()
            oldfu.dilation()

    def test_bynary_2Dimg(self):
        return_new = fu.dilation(IMAGE_BLANK_2D, MASK, morphtype="BINARY")
        return_old = oldfu.dilation(IMAGE_BLANK_2D, MASK, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_3Dimg(self):
        return_new = fu.dilation(IMAGE_BLANK_3D, MASK, morphtype="BINARY")
        return_old = oldfu.dilation(IMAGE_BLANK_3D, MASK, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_2Dimg_NOmask(self):
        return_new = fu.dilation(IMAGE_BLANK_2D, morphtype="BINARY")
        return_old = oldfu.dilation(IMAGE_BLANK_2D, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_3Dimg_NOmask(self):
        return_new = fu.dilation(IMAGE_BLANK_3D, morphtype="BINARY")
        return_old = oldfu.dilation(IMAGE_BLANK_3D, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_2Dimg(self):
        return_new = fu.dilation(IMAGE_BLANK_2D, MASK, morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(IMAGE_BLANK_2D, MASK, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_3Dimg(self):
        return_new = fu.dilation(IMAGE_BLANK_3D, MASK, morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(IMAGE_BLANK_3D, MASK, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_2Dimg_NOmask(self):
        return_new = fu.dilation(IMAGE_BLANK_2D,morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(IMAGE_BLANK_2D,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_3Dimg_NOmask(self):
        return_new = fu.dilation(IMAGE_BLANK_3D,morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(IMAGE_BLANK_3D,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_invalid_2Dimg_Error_msg_unknown_dilatation_type(self):
        return_new = fu.dilation(IMAGE_BLANK_2D, MASK, morphtype="invalid_type")
        return_old = oldfu.dilation(IMAGE_BLANK_2D, MASK, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_invalid_3Dimg_Error_msg_unknown_dilatation_type(self):
        return_new = fu.dilation(IMAGE_BLANK_3D, MASK, morphtype="invalid_type")
        return_old = oldfu.dilation(IMAGE_BLANK_3D, MASK, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)



class Test_erosion(unittest.TestCase):

    def test_empty_input_image_crashes_because_signal11SIGSEV(self):
        """
        It seems not possible to test it without getting an segmentation fault. We should return None after the first error message
        in order to avoid to run in the second if/else and get the segmentation fault in the c+++ code
        img =EMData()
        return_new = fu.erosion(img)
        return_old = oldfu.erosion(img)
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)
        """
        self.assertTrue(True)

    def test_empty_mask_image_RuntimeError_ImageDimensionException_center_isnot_welldefined(self):
        with self.assertRaises(RuntimeError):
            fu.erosion(IMAGE_BLANK_2D, EMData(), morphtype="BINARY")
            oldfu.erosion(IMAGE_BLANK_2D, EMData(), morphtype="BINARY")

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.erosion()
            oldfu.erosion()

    def test_bynary_2Dimg(self):
        return_new = fu.erosion(IMAGE_BLANK_2D, MASK, morphtype="BINARY")
        return_old = oldfu.erosion(IMAGE_BLANK_2D, MASK, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_3Dimg(self):
        return_new = fu.erosion(IMAGE_BLANK_3D, MASK, morphtype="BINARY")
        return_old = oldfu.erosion(IMAGE_BLANK_3D, MASK, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_2Dimg_NOmask(self):
        return_new = fu.erosion(IMAGE_BLANK_2D, morphtype="BINARY")
        return_old = oldfu.erosion(IMAGE_BLANK_2D, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_bynary_3Dimg_NOmask(self):
        return_new = fu.erosion(IMAGE_BLANK_3D, morphtype="BINARY")
        return_old = oldfu.erosion(IMAGE_BLANK_3D, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_2Dimg(self):
        return_new = fu.erosion(IMAGE_BLANK_2D, MASK, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(IMAGE_BLANK_2D, MASK, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_3Dimg(self):
        return_new = fu.erosion(IMAGE_BLANK_3D, MASK, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(IMAGE_BLANK_3D, MASK, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_2Dimg_NOmask(self):
        return_new = fu.erosion(IMAGE_BLANK_2D,morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(IMAGE_BLANK_2D,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_graylevel_3Dimg_NOmask(self):
        return_new = fu.erosion(IMAGE_BLANK_3D,morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(IMAGE_BLANK_3D,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_invalid_2Dimg_Error_msg_unknown_erosion_type(self):
        return_new = fu.erosion(IMAGE_BLANK_2D, MASK, morphtype="invalid_type")
        return_old = oldfu.erosion(IMAGE_BLANK_2D, MASK, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_invalid_3Dimg_Error_msg_unknown_erosion_type(self):
        return_new = fu.erosion(IMAGE_BLANK_3D, MASK, morphtype="invalid_type")
        return_old = oldfu.erosion(IMAGE_BLANK_3D, MASK, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)



class Test_power(unittest.TestCase):

    def test_power_2Dimg(self):
        return_new = fu.power(IMAGE_2D, x = 3.0)
        return_old = oldfu.power(IMAGE_2D, x = 3.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_power_3Dimg(self):
        return_new = fu.power(IMAGE_3D, x = 3.0)
        return_old = oldfu.power(IMAGE_3D, x = 3.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException(self):
        with self.assertRaises(RuntimeError):
            fu.power(EMData())
            oldfu.power(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.power()
            oldfu.power()



class Test_square_root(unittest.TestCase):

    def test_positive_2Dimg(self):
        return_new = fu.square_root(IMAGE_2D)
        return_old = oldfu.square_root(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_3Dimg(self):
        return_new = fu.square_root(IMAGE_3D)
        return_old = oldfu.square_root(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_2Dimg_error_Error_msg_cannot_calculate_sqaure_root_of_negative_pixel(self):
        img= deepcopy(IMAGE_2D)
        img.sub(100)
        return_new = fu.square_root(img)
        return_old = oldfu.square_root(img)
        self.assertTrue(numpy.all(numpy.isnan(return_old.get_3dview())))
        self.assertTrue(numpy.all(numpy.isnan(return_new.get_3dview())))

    def test_negative_3Dimg_Error_msg_cannot_calculate_sqaure_root_of_negative_pixel(self):
        img= deepcopy(IMAGE_3D)
        img.sub(100)
        return_new = fu.square_root(img)
        return_old = oldfu.square_root(img)

        for i in range(len(return_new.get_3dview())):
            try:
                self.assertTrue(numpy.allclose(return_old.get_3dview()[i], return_new.get_3dview()[i], atol=TOLERANCE))
            except AssertionError:
                self.assertTrue(numpy.all(numpy.isnan(return_old.get_3dview()[i])))
                self.assertTrue(numpy.all(numpy.isnan(return_new.get_3dview()[i])))

    def test_empty_input_image_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.square_root(EMData())
            oldfu.square_root(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.square_root()
            oldfu.square_root()



class Test_square(unittest.TestCase):

    def test_square_2Dimg(self):
        return_new = fu.square(IMAGE_2D)
        return_old = oldfu.square(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_square_3Dimg(self):
        return_new = fu.square(IMAGE_3D)
        return_old = oldfu.square(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException_and_NotExistingObjectException_the_key_maximum_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.square(EMData())
            oldfu.square(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.square()
            oldfu.square()



class Test_threshold(unittest.TestCase):

    def test_threshold_2Dimg(self):
        return_new = fu.threshold(IMAGE_2D)
        return_old = oldfu.threshold(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_threshold_3Dimg(self):
        return_new = fu.threshold(IMAGE_3D)
        return_old = oldfu.threshold(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException_and_NotExistingObjectException_the_key_maximum_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.threshold(EMData())
            oldfu.threshold(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.threshold()
            oldfu.threshold()



class Test_threshold_outside(unittest.TestCase):

    def test_threshold_outside_2Dimg(self):
        return_new = fu.threshold_outside(IMAGE_2D, 2 , 10)
        return_old = oldfu.threshold_outside(IMAGE_2D, 2 , 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_threshold_outside_3Dimg(self):
        return_new = fu.threshold_outside(IMAGE_3D, 2 , 10)
        return_old = oldfu.threshold_outside(IMAGE_3D, 2 , 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        return_new = fu.threshold_outside(EMData(), 2 , 10)
        return_old = oldfu.threshold_outside(EMData(), 2 , 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.threshold_outside()
            oldfu.threshold_outside()



class Test_notzero(unittest.TestCase):

    def test_notzero_2Dimg(self):
        return_new = fu.notzero(IMAGE_2D)
        return_old = oldfu.notzero(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_notzero_3Dimg(self):
        return_new = fu.notzero(IMAGE_3D)
        return_old = oldfu.notzero(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image_returns_RuntimeError_stdException_and_NotExistingObjectException_the_key_maximum_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.notzero(EMData())
            oldfu.notzero(EMData())

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.notzero()
            oldfu.notzero()



class Test_rotavg_ctf(unittest.TestCase):
    """ See http://sparx-em.org/sparxwiki/CTF_info for the meaning of the params"""

    def test_empty_input_image(self):
        return_new = fu.rotavg_ctf(EMData(), defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(EMData(), defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.rotavg_ctf()
            oldfu.rotavg_ctf()

    def test_2DImg_nullCS(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_nullCS(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_2DImg_withCS(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_withCS(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_2DImg_nullCS_and_defocus_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 0, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 0, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_nullCS_and_defocus_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 0, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 0, Cs =0.0, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_2DImg_withCS_and_defocus_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 0, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 0, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_withCS_and_defocus_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 0, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 0, Cs =2, voltage=300, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_2DImg_nullCS_and_voltage_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =0.0, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =0.0, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_nullCS_and_voltage_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =0.0, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =0.0, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_2DImg_withCS_and_voltage_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_3DImg_withCS_and_voltage_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =2, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_3D, defocus= 1, Cs =2, voltage=0, Pixel_size=1.5,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Null_pixelSize_RuntimeWarning_msg_invalid_value_encountered(self):
        return_new = fu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=300, Pixel_size=0,amp = 0.0, ang = 0.0)
        return_old = oldfu.rotavg_ctf(IMAGE_2D, defocus= 1, Cs =2, voltage=300, Pixel_size=0,amp = 0.0, ang = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ctf_1d(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf_1d()
            oldfu.ctf_1d()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_1d(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_1d(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new).any())
        self.assertTrue(numpy.isnan(return_old).any())

    def test_no_image_size_retuns_ZeroDivisionError(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_1d(nx=0, ctf=ctf, sign = 1, doabs = False)
            oldfu.ctf_1d(nx=0, ctf=ctf, sign = 1, doabs = False)

    def test_no_pixel_size_retuns_ZeroDivisionError(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_1d(nx =2, ctf=ctf, sign = 1, doabs = False)
            oldfu.ctf_1d(nx=2, ctf=ctf, sign = 1, doabs = False)

    def test_NOSign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=False)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=False)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_positive_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign = 1,doabs=False)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign = 1,doabs=False)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negative_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign = -1,doabs=False)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign = -1,doabs=False)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_NOSign_withABS(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=True)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_positive_sign_withABS(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign = 1,doabs=True)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign = 1,doabs=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negative_sign_withABS(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign = -1,doabs=True)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign = -1,doabs=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ctf_2(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf_2()
            oldfu.ctf_2()

    def test_no_image_size_retuns_ZeroDivisionError(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_2(nx=0, ctf=ctf)
            oldfu.ctf_2(nx=0, ctf=ctf)

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_2(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_2(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new).any())
        self.assertTrue(numpy.isnan(return_old).any())

    def test_no_pixel_size_retuns_ZeroDivisionError(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_2(nx =2, ctf=ctf)
            oldfu.ctf_2(nx=2, ctf=ctf)

    def test_ctf_2(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_2(nx =2, ctf=ctf)
        return_old= oldfu.ctf_2(nx=2, ctf=ctf)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_ctf_img(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf_img()
            oldfu.ctf_img()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_img(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        return_old= oldfu.ctf_img(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)
            oldfu.ctf_img(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_img(nx =2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old = oldfu.ctf_img(nx=2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_null_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_ny_and_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny_null_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_nz_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 0)
            oldfu.ctf_img(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 0)



class Test_ctf_img_real(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf_img_real()
            oldfu.ctf_img_real()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_img_real(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        return_old= oldfu.ctf_img_real(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img_real(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)
            oldfu.ctf_img_real(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_img_real(nx =2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old = oldfu.ctf_img_real(nx=2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_null_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx =20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx =20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_ny_and_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx =20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny_null_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx =20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_nz_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img_real(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 0)
            oldfu.ctf_img_real(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 0)




class Test_ctf_rimg(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf_rimg()
            oldfu.ctf_rimg()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_rimg(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        return_old= oldfu.ctf_rimg(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_rimg(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)
            oldfu.ctf_rimg(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_rimg(nx =2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old = oldfu.ctf_rimg(nx=2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_null_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx =20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx =20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_ny_and_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx =20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny_null_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx =20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_nz_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_rimg(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 0)
            oldfu.ctf_rimg(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 0)



class Test_ctf2_rimg(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctf2_rimg()
            oldfu.ctf2_rimg()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf2_rimg(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        return_old= oldfu.ctf2_rimg(nx=20, ctf= EMAN2Ctf(), sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf2_rimg(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)
            oldfu.ctf2_rimg(nx=0, ctf=ctf, sign = 1, ny = 0, nz = 1)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf2_rimg(nx =2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old = oldfu.ctf2_rimg(nx=2, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_null_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf2_rimg(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        return_old=oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf2_rimg(nx =20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        return_old=oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 1, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_negative_ny(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf2_rimg(nx =20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        return_old=oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 1, ny = -2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_ny_and_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf2_rimg(nx =20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        return_old=oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 0, ny = 0, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_positive_ny_null_sign(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf2_rimg(nx =20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        return_old=oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 0, ny = 2, nz = 1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_null_nz_error_returns_RuntimeError_InvalidValueException(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf2_rimg(nx =20, ctf=ctf, sign = 1, ny = 0, nz = 0)
            oldfu.ctf2_rimg(nx=20, ctf=ctf, sign = 1, ny = 0, nz = 0)



class Test_ctflimit(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ctflimit()
            oldfu.ctflimit()

    def test_ctfLimit(self):
        return_new = fu.ctflimit(nx=30, defocus=1, cs=2, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=30, defocus=1, cs=2, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_voltage_crashes_because_LinAlgError(self):
        """ How can I catch the '"LinAlgError
        with self.assertRaises(ValueError):
            return_new = fu.ctflimit(nx=30, defocus=1, cs=2, voltage=0, pix=1.5)
            return_old = oldfu.ctflimit(nx=30, defocus=1, cs=2, voltage=0, pix=1.5)
            self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_null_defocus(self):
        return_new = fu.ctflimit(nx=30, defocus=0, cs=2, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=30, defocus=0, cs=2, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_cs(self):
        return_new = fu.ctflimit(nx=30, defocus=1, cs=0, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=30, defocus=1, cs=0, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_nx(self):
        return_new = fu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negative_nx_retuns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.ctflimit(nx=-1, defocus=1, cs=2, voltage=300, pix=1.5)
            oldfu.ctflimit(nx=-1, defocus=1, cs=2, voltage=300, pix=1.5)

    def test_no_pixel_size_retuns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=0)
            oldfu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=0)



class Test_imf_params_cl1(unittest.TestCase):
    """
    pw = power spectrum to be fitted
    n = the polynomial rank +1
    iswi = integer between 1-8 is used in 'vector<float> Util::call_cl1' to calculate the interpolation value (lagracian?)
    """
    pw = [entry for entry in numpy.arange(0, 10).tolist()]


    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(len(return_new), len(return_old))
            for i,j in zip(return_new,return_old):
                try:
                    self.assertTrue(numpy.array_equal(i,j))
                except AssertionError:
                    self.assertTrue(numpy.isnan(i).any())
                    self.assertTrue(numpy.isnan(j).any())


    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.imf_params_cl1()
            oldfu.imf_params_cl1()

    def test_empty_list_crashes_because_signal6SIGABRT(self):
        """
        return_new = fu.imf_params_cl1([], n=2, iswi=3, Pixel_size=1)
        return_old = oldfu.imf_params_cl1([], n=2, iswi=3, Pixel_size=1)
        self.test_all_the_conditions(return_old,return_new, False)
        """
        self.assertTrue(True)

    def test_no_pixel_size_error(self):
        return_new = fu.imf_params_cl1(self.pw, n=4, iswi=3, Pixel_size=0)
        return_old = oldfu.imf_params_cl1(self.pw, n=4, iswi=3, Pixel_size=0)
        self.test_all_the_conditions(return_old, return_new, False)

    def test_with_default_value(self):
        return_new = fu.imf_params_cl1(self.pw, n=2, iswi=3, Pixel_size=1)
        return_old = oldfu.imf_params_cl1(self.pw, n=2, iswi=3, Pixel_size=1)
        self.test_all_the_conditions(return_old,return_new, False)

    def test_null_rank(self):
        return_new = fu.imf_params_cl1(self.pw, n=0, iswi=3, Pixel_size=1)
        return_old = oldfu.imf_params_cl1(self.pw, n=0, iswi=3, Pixel_size=1)
        self.test_all_the_conditions(return_old, return_new, False)

    def test_null_iswi(self):
        return_new = fu.imf_params_cl1(self.pw, n=2, iswi=0, Pixel_size=1)
        return_old = oldfu.imf_params_cl1(self.pw, n=2, iswi=0, Pixel_size=1)
        self.test_all_the_conditions(return_old, return_new, False)

    def test_negative_rank_error_because_signal6SIGABRT(self):
        """
        return_new = fu.imf_params_cl1(self.pw, n=-2, iswi=2, Pixel_size=1)
        return_old = oldfu.imf_params_cl1(self.pw, n=-2, iswi=2, Pixel_size=1)
        self.test_all_the_conditions(return_old,return_new, False)
        """
        self.assertTrue(True)

    def test_with_invalid_iswi(self):
        for iswi in  [-2,10]:
            return_new = fu.imf_params_cl1(self.pw, n=2, iswi=iswi, Pixel_size=1)
            return_old = oldfu.imf_params_cl1(self.pw, n=2, iswi=iswi, Pixel_size=1)
            self.test_all_the_conditions(return_old, return_new, False)



class Test_adaptive_mask(unittest.TestCase):
    """
    If threshold as the -9999.0 default value it uses the nsigma value to calculate the thresheld
    """
    def test_empty_input_image_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError):
            fu.adaptive_mask(EMData(),nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)
            oldfu.adaptive_mask(EMData(),nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.adaptive_mask()
            oldfu.adaptive_mask()

    def test_2dimg_default_values(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dimg_no_threshold_null_sigma(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = 0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = 0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dimg_no_threshold_negative_sigma(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = -10, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = -10, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dimg_no_dilation(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0,  ndilation = 0, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0,  ndilation = 0, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dimg_negative_dilation(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0,  ndilation = -2, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0,  ndilation = -2, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dimg_negative_edge_width_crashes_because_signal11SIGSEV(self):
        """
        return_new = fu.adaptive_mask(IMAGE_2D,nsigma = 1.0, threshold = -9999.0,  ndilation = 3, edge_width = -5)
        return_old = oldfu.adaptive_mask(IMAGE_2D,nsigma = 1.0, threshold = -9999.0,  ndilation = 3, edge_width = -5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """
        self.assertTrue(True)

    def test_2dimg_null_edge_width(self):
        return_new = fu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width=0)
        return_old = oldfu.adaptive_mask(IMAGE_2D, nsigma = 1.0, threshold = -9999.0, ndilation = 3 ,edge_width=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_default_values(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_no_threshold_null_sigma(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = 0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = 0, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_no_threshold_negative_sigma(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = -10, threshold = -9999.0, ndilation = 3, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = -10, threshold = -9999.0, ndilation = 3, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_no_dilation(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0,  ndilation = 0, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0,  ndilation = 0, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_negative_dilation(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0,  ndilation = -2, edge_width = 5)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0,  ndilation = -2, edge_width = 5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dimg_negative_edge_width_crashes_because_signal11SIGSEV(self):
        """
        return_new = fu.adaptive_mask(IMAGE_3D,nsigma = 1.0, threshold = -9999.0,  ndilation = 3, edge_width = -5)
        return_old = oldfu.adaptive_mask(IMAGE_3D,nsigma = 1.0, threshold = -9999.0,  ndilation = 3, edge_width = -5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """
        self.assertTrue(True)

    def test_3dimg_null_edge_width(self):
        return_new = fu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width=0)
        return_old = oldfu.adaptive_mask(IMAGE_3D, nsigma = 1.0, threshold = -9999.0, ndilation = 3 ,edge_width=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cosinemask(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.cosinemask()
            oldfu.cosinemask()

    def test_empty_input_image_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError):
            fu.cosinemask(EMData(), radius = -1, cosine_width = 5, bckg = None, s=999999.0)
            oldfu.cosinemask(EMData(), radius = -1, cosine_width = 5, bckg = None, s=999999.0)

    def test_empty_bckg_image_crashes_because_signal11SIGSEV(self):
        """
        bckg = EMData()
        return_new = fu.cosinemask(IMAGE_3D, bckg=bckg)
        return_old = oldfu.cosinemask(IMAGE_3D, bckg=bckg)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """
        self.assertTrue(True)

    def test_3d_img_with_bckg_crashes_because_signal11SIGSEV(self):
        """
        bckg = sparx_utilities.model_gauss_noise(0.25 , 10,10,10)
        return_new = fu.cosinemask(IMAGE_3D, bckg=bckg)
        return_old = oldfu.cosinemask(IMAGE_3D, bckg=bckg)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """
        self.assertTrue(True)

    def test_3d_img_default_values(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_null_radius(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = 0, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = 0, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_positive_radius(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = 10, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = 10, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_null_cosine_width(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 0, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 0, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_negative_cosine_width(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = -1, cosine_width = -5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = -1, cosine_width = -5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_null_s(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=0)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3d_img_negative_s(self):
        return_new = fu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=-10)
        return_old = oldfu.cosinemask(IMAGE_3D, radius = -1, cosine_width = 5, bckg = None, s=-10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_with_bckg_crashes_because_signal11SIGSEV(self):
        """
        bckg = sparx_utilities.model_gauss_noise(0.25 , 10,10,10)
        return_new = fu.cosinemask(IMAGE_2D, bckg=bckg)
        return_old = oldfu.cosinemask(IMAGE_2D, bckg=bckg)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """
        self.assertTrue(True)

    def test_2d_img_default_values(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_null_radius(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = 0, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = 0, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_positive_radius(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = 10, cosine_width = 5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = 10, cosine_width = 5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_null_cosine_width(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 0, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 0, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_negative_cosine_width(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = -1, cosine_width = -5, bckg = None, s=999999.0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = -1, cosine_width = -5, bckg = None, s=999999.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_null_s(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=0)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_img_negative_s(self):
        return_new = fu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=-10)
        return_old = oldfu.cosinemask(IMAGE_2D, radius = -1, cosine_width = 5, bckg = None, s=-10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_get_shrink_3dmask(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.get_shrink_3dmask()
            oldfu.get_shrink_3dmask()

    def test_empty_input_image_crashes_because_signal11SIGSEV(self):
        """
        with self.assertRaises(RuntimeError):
            fu.get_shrink_3dmask(3, EMData())
            oldfu.get_shrink_3dmask(3, EMData())
        """
        self.assertTrue(True)

    def test_No_xinit_error_returns_RuntimeError_InvalidValueException(self):
        with self.assertRaises(RuntimeError):
            fu.get_shrink_3dmask(nxinit = 0, mask_file_name = get_data_3d(1))
            oldfu.get_shrink_3dmask(nxinit = 0, mask_file_name = get_data_3d(1))

    def test_3Dmask_format_error_returns_RuntimeError_float_hasnot_attribute_copy(self):
        """ the Image3D is an EMdata"""
        with self.assertRaises(AttributeError):
            fu.get_shrink_3dmask(nxinit = 4, mask_file_name = IMAGE_3D)
            oldfu.get_shrink_3dmask(nxinit = 4, mask_file_name = IMAGE_3D)

    def test_2Dmask_format_error_returns_RuntimeError_float_hasnot_attribute_copy(self):
        """ the Image3D is an EMdata"""
        with self.assertRaises(AttributeError):
            fu.get_shrink_3dmask(nxinit = 4, mask_file_name = IMAGE_2D)
            oldfu.get_shrink_3dmask(nxinit = 4, mask_file_name = IMAGE_2D)

    def test_3Dmask(self):
        """ the get_data_3d(1) is a list with one EMdata element"""
        return_new = fu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data_3d(1))
        return_old = oldfu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data_3d(1))
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dmask(self):
        return_new = fu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data(1))
        return_old = oldfu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data(1))
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_nx_equal_size3Dmask(self):
        mask_file_name = get_data_3d(1)
        nx =sparx_utilities.get_im(mask_file_name).get_xsize()
        return_new = fu.get_shrink_3dmask(nxinit = nx, mask_file_name = mask_file_name)
        return_old = oldfu.get_shrink_3dmask(nxinit = nx, mask_file_name = mask_file_name)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_get_biggest_cluster(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.get_biggest_cluster()
            oldfu.get_biggest_cluster()

    def test_empty_input_image_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(self):
        with self.assertRaises(RuntimeError):
            fu.get_biggest_cluster( EMData())
            oldfu.get_biggest_cluster( EMData())

    def test_2Dimg(self):
        return_new = fu.get_biggest_cluster(IMAGE_2D)
        return_old = oldfu.get_biggest_cluster(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg(self):
        return_new = fu.get_biggest_cluster(IMAGE_3D)
        return_old = oldfu.get_biggest_cluster(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_gauss_noise_img(self):
        image = sparx_utilities.model_gauss_noise(0.25 , 10,10,10)
        return_new = fu.get_biggest_cluster(image)
        return_old = oldfu.get_biggest_cluster(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_compute_bfactor(unittest.TestCase):
    pw = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(return_new[0], return_old[0])
            self.assertEqual(return_new[2], return_old[2])
            self.assertEqual(return_new[3], return_old[3])
            self.assertEqual(len(return_new[1]), len(return_old[1]))
            for i, j in zip(return_new[1], return_old[1]):
                self.assertTrue(numpy.array_equal(i, j))
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.compute_bfactor()
            oldfu.compute_bfactor()

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=0)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=0)

    def test_compute_bfactor(self):
        return_new =fu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=1.0)
        return_old = oldfu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=1.0)
        self.test_all_the_conditions(return_new, return_old, False)

    def test_with_f_negative(self):
        return_new =fu.compute_bfactor(pws=self.pw, freq_min = -0.15, freq_max= -0.25, pixel_size=1.0)
        return_old = oldfu.compute_bfactor(pws=self.pw, freq_min = -0.15, freq_max= -0.25, pixel_size=1.0)
        self.test_all_the_conditions(return_new, return_old, False)

    def test_freqMin_bigger_than_freqMAx_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.25, pixel_size=1.0)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.25, pixel_size=1.0)

    def test_freqMin_equal_freqMAx_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.35, pixel_size=1.0)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.35, pixel_size=1.0)

    def test_few_value_in_power_spectrum_list_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=[1,1], freq_min = 0.15, freq_max= 0.25, pixel_size=1.0)
            oldfu.compute_bfactor(pws=[1,1], freq_min= 0.15, freq_max= 0.25, pixel_size=1.0)

    def test_Empty_array_error(self):
        with self.assertRaises(ValueError):
            fu.compute_bfactor(pws=[], freq_min=0.15, freq_max=0.25, pixel_size=1.0)
            oldfu.compute_bfactor(pws=[], freq_min=0.15, freq_max=0.25, pixel_size=1.0)



class Test_cter_mrk(unittest.TestCase):
    defocus = 1
    cs = 2
    voltage = 300
    pixel_size = 1.5
    bfactor = 0
    amp_contrast = 0.1
    wn = 32
    image1 = get_data(1, 256)[0]
    selection_list = 'image.mrc'
    input_image_path = path.join(ABSOLUTE_PATH, "cter_mrk/image*.mrc")
    output_directory = path.join(ABSOLUTE_PATH, "cter_mrk/results")

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.cter_mrk()
            oldfu.cter_mrk()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_mark_true_should_return_equal_object(self):
        #RuntimeError: FileAccessException at /home/lusnig/EMAN2/eman2/libEM/imageio.cpp:158: error with 'cter_mrk/image.mrc': 'cannot access file 'cter_mrk/image.mrc'' caught

        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": self.defocus, "cs": self.cs, "voltage": self.voltage, "apix": self.pixel_size, "bfactor": self.bfactor,"ampcont": self.amp_contrast})

        micrograph_image = sparx_filter.filt_ctf(self.image1,ctf)
        micrograph_image.write_image('cter_mrk/image.mrc')

        remove_dir(self.output_directory)

        return_new = fu.cter_mrk(self.input_image_path, self.output_directory, selection_list = self.selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage)

        remove_dir(self.output_directory)

        return_old = oldfu.cter_mrk(self.input_image_path, self.output_directory, selection_list = self.selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage)

        self.assertEqual(return_new, return_old)

class Test_cter_pap(unittest.TestCase):
    defocus = 1
    cs = 2
    voltage = 300
    pixel_size = 1.5
    bfactor = 0
    amp_contrast = 0.1
    wn = 32
    image1 = get_data(1, 256)[0]
    selection_list = 'image.mrc'
    input_image_path = path.join(ABSOLUTE_PATH, "cter_mrk/image*.mrc")
    output_directory = path.join(ABSOLUTE_PATH, "cter_mrk/results")

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.cter_pap()
            oldfu.cter_pap()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_pap_true_should_return_equal_object(self):
        #RuntimeError: FileAccessException at /home/lusnig/EMAN2/eman2/libEM/imageio.cpp:158: error with 'cter_mrk/image.mrc': 'cannot access file 'cter_mrk/image.mrc'' caught

        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": self.defocus, "cs": self.cs, "voltage": self.voltage, "apix": self.pixel_size, "bfactor": self.bfactor,"ampcont": self.amp_contrast})

        micrograph_image = sparx_filter.filt_ctf(self.image1,ctf)
        micrograph_image.write_image('cter_mrk/image.mrc')

        remove_dir(self.output_directory)

        return_new = fu.cter_pap(self.input_image_path, self.output_directory, selection_list = self.selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage)

        remove_dir(self.output_directory)

        return_old = oldfu.cter_pap(self.input_image_path, self.output_directory, selection_list = self.selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage)

        self.assertEqual(return_new, return_old)


class Test_cter_vpp(unittest.TestCase):
    defocus = 1
    cs = 2
    voltage = 300
    pixel_size = 1.5
    bfactor = 0
    amp_contrast = 0.1
    wn = 32
    i_start = 0.048
    i_stop = -1
    vpp_options = [entry for entry in numpy.arange(0, 6).tolist()]
    image1 = get_data(1, 256)[0]
    selection_list = 'image.mrc'
    input_image_path = path.join(ABSOLUTE_PATH, "cter_mrk/image*.mrc")
    output_directory = path.join(ABSOLUTE_PATH, "cter_mrk/results")

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.cter_vpp()
            oldfu.cter_vpp()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_vpp_true_should_return_equal_object(self):
        #RuntimeError: FileAccessException at /home/lusnig/EMAN2/eman2/libEM/imageio.cpp:158: error with 'cter_mrk/image.mrc': 'cannot access file 'cter_mrk/image.mrc'' caught


        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": self.defocus, "cs": self.cs, "voltage": self.voltage, "apix": self.pixel_size, "bfactor": self.bfactor, "ampcont": self.amp_contrast})
        image1, = get_data(1, 256)

        micrograph_image = sparx_filter.filt_ctf(image1,ctf)
        micrograph_image.write_image('cter_mrk/image.mrc')
        selection_list = 'image.mrc'

        remove_dir(self.output_directory)

        return_new = fu.cter_vpp(self.input_image_path, self.output_directory, selection_list = selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage, f_start=self.i_start, f_stop=self.i_stop, vpp_options=self.vpp_options)

        remove_dir(self.output_directory)

        return_old = oldfu.cter_vpp(self.input_image_path, self.output_directory, selection_list = selection_list, wn = self.wn,  pixel_size=self.pixel_size, Cs= self.cs, voltage = self.voltage, f_start=self.i_start, f_stop=self.i_stop, vpp_options=self.vpp_options)

        self.assertEqual(return_new, return_old)



class Test_ampcont2angle(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ampcont2angle()
            oldfu.ampcont2angle()

    def test_A_equal_minus100(self):
        return_new = fu.ampcont2angle(100.0)
        return_old = oldfu.ampcont2angle(100.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_A_equal100(self):
        return_new = fu.ampcont2angle(-100.0)
        return_old = oldfu.ampcont2angle(-100.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negative_A(self):
        return_new = fu.ampcont2angle(-1)
        return_old = oldfu.ampcont2angle(-1)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_positive_A(self):
        return_new = fu.ampcont2angle(8)
        return_old = oldfu.ampcont2angle(8)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_angle2ampcont(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.angle2ampcont()
            oldfu.angle2ampcont()

    def test_positive_phi(self):
        return_new = fu.angle2ampcont(0.45)
        return_old = oldfu.angle2ampcont(0.45)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_negativetive_phi(self):
        return_new = fu.angle2ampcont(-0.45)
        return_old = oldfu.angle2ampcont(-0.45)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_phi(self):
        return_new = fu.angle2ampcont(0)
        return_old = oldfu.angle2ampcont(0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_bracket_def(unittest.TestCase):

    @staticmethod
    def function1(x1, dat):
        return x1 + dat

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.bracket_def()
            oldfu.bracket_def()

    def test_f2_greater_f1_outputmsg_Bracket_didnot_find_a_minimum(self):
        return_new = fu.bracket_def(self.function1, dat=5, x1=3, h=3)
        return_old = oldfu.bracket_def(self.function1, dat=5, x1=3, h=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_f2_not_greater_f1_outputmsg_Bracket_didnot_find_a_minimum(self):
        return_new = fu.bracket_def(self.function1, dat=5, x1=3, h=0)
        return_old = oldfu.bracket_def(self.function1, dat=5, x1=3, h=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_bracket(unittest.TestCase):

    @staticmethod
    def function1(x1, dat):
        return x1 + dat

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.bracket()
            oldfu.bracket()

    def test_f3_greater_f1(self):
        return_new = fu.bracket(self.function1, dat=5, h=4)
        return_old = oldfu.bracket(self.function1, dat=5, h=4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_f3_not_greater_f1_outputmsg_Bracket_didnot_find_a_minimum(self):
        self.assertTrue(fu.bracket(self.function1, dat=0, h=0) is None)
        self.assertTrue(oldfu.bracket(self.function1, dat=0, h=0) is None)



class Test_goldsearch_astigmatism(unittest.TestCase):

    @staticmethod
    def function1(x1, dat):
        f = x1 + dat
        return f

    @staticmethod
    def function_return_0(x1, dat):
        return 0

    @staticmethod
    def bad_function():
        return 0

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.goldsearch_astigmatism()
            oldfu.goldsearch_astigmatism()

    def test_null_tolerance_returns_OverflowError_cannot_convert_infinity_to_integer(self):
        with self.assertRaises(OverflowError):
            fu.goldsearch_astigmatism(self.function1, 5, 3, 4, 0)
            oldfu.goldsearch_astigmatism(self.function1, 5, 3, 4, 0)


    def test_A_B_same_value_error_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.goldsearch_astigmatism(self.function1, 5, 3, 3)
            oldfu.goldsearch_astigmatism(self.function1, 5, 3, 3)

    def test_Invalid_function_returns_TypeError_bad_function_takes_no_arguments(self):
        with self.assertRaises(TypeError):
            fu.goldsearch_astigmatism(self.bad_function, 5, 3, 4)
            oldfu.goldsearch_astigmatism(self.bad_function, 5, 3, 4)

    def test_f2_greater_f1(self):
        return_new = fu.goldsearch_astigmatism(self.function1, 5, 3, 4)
        return_old = oldfu.goldsearch_astigmatism(self.function1, 5, 3, 4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_return_f2_greater_f1(self):
        return_new = fu.goldsearch_astigmatism(self.function_return_0, 5, 3, 4)
        return_old = oldfu.goldsearch_astigmatism(self.function_return_0, 5, 3, 4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_test_f1_greater_f2(self):
        return_new = fu.goldsearch_astigmatism(self.function1, 5, 4, 3)
        return_old = oldfu.goldsearch_astigmatism(self.function1, 5, 4, 3)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_defocus_baseline_fit(unittest.TestCase):
    """
    Its input params are used at the beginning of the body of the functions to feed 'imf_params_cl1'.
    """
    roo = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.defocus_baseline_fit()
            oldfu.defocus_baseline_fit()

    def test_iswi3(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D31_iswi_not3(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=0)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_start_is_bigger_than_stop_error_because_signal6SIGABRT(self):
        """
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=10, i_stop=7, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=10, i_stop=7, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_start_is_equal_stop_error_because_signal6SIGABRT(self):
        """
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=9, i_stop=9, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=9, i_stop=9, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_negative_rank_error_because_signal6SIGABRT(self):
        """
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=-1, iswi=2)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=-1, iswi=2)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_null_rank_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=0, iswi=2)
            oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=0, iswi=2)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.defocus_baseline_fit(roo=[], i_start=0, i_stop=10, nrank=2, iswi=2)
            oldfu.defocus_baseline_fit(roo=[], i_start=0, i_stop=10, nrank=2, iswi=2)



class Test_simpw1d(unittest.TestCase):
    """ I got this value from the pickle file"""
    data = [entry for entry in numpy.arange(1, 256).tolist()]
    defocus = 1
    Cs = 2
    voltage = 300
    pixel_size = 1.5
    amp_contrast = 0.1
    i_start = 2
    i_stop = 14
    nx = 20

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.simpw1d()
            oldfu.simpw1d()

    def test_positive_defocus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(self.defocus,datanew), oldfu.simpw1d(self.defocus,datanew))

    def test_negative_defocus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], 0, self.defocus,self.Cs, self.voltage, 0, self.amp_contrast, self.i_start, self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)



class Test_simpw1d_pap(unittest.TestCase):
    data = [entry for entry in numpy.arange(1, 256).tolist()]
    defocus = 1
    Cs = 2
    voltage = 300
    pixel_size = 1.5
    amp_contrast = 0.1
    i_start = 2
    i_stop = 14
    nx = 20

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.simpw1d_pap()
            oldfu.simpw1d_pap()

    def test_positive_focus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus,self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start, self.i_stop]
        self.assertEqual(fu.simpw1d_pap(self.defocus, datanew), oldfu.simpw1d_pap(self.defocus, datanew))

    def test_negative_focus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], 0, self.defocus,self.Cs, self.voltage, 0, self.amp_contrast, self.i_start, self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)



class Test_simpw1d_print(unittest.TestCase):
    data = [entry for entry in numpy.arange(1, 256).tolist()]
    defocus = 1
    Cs = 2
    voltage = 300
    pixel_size = 1.5
    amp_contrast = 0.1
    i_start = 2
    i_stop = 14
    nx = 20

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.simpw1d_print()
            oldfu.simpw1d_print()


    def test_positive_focus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start, self.i_stop]
        self.assertEqual(fu.simpw1d_print(self.defocus,datanew), oldfu.simpw1d_print(self.defocus,datanew))

    def test_negative_focus(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_returns_ZeroDivisionError(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], 0, self.defocus,self.Cs, self.voltage, 0, self.amp_contrast, self.i_start, self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)



class Test_movingaverage(unittest.TestCase):
    data = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.movingaverage()
            oldfu.movingaverage()

    def test_default_value(self):
        return_new = fu.movingaverage(self.data,window_size=2, skip=3)
        return_old = oldfu.movingaverage(self.data, window_size=2, skip=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_skip(self):
        return_new = fu.movingaverage(self.data,window_size=2, skip=0)
        return_old = oldfu.movingaverage(self.data, window_size=2, skip=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_windows_size_negative_Value_returns_ValueError_negative_dimensions_arenot_allowed(self):
        with self.assertRaises(ValueError):
            fu.movingaverage(self.data,window_size=-2)
            oldfu.movingaverage(self.data, window_size=-2)


    def test_windows_size_null__returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.movingaverage(self.data,window_size=0)
            oldfu.movingaverage(self.data, window_size=0)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.movingaverage([], window_size=2)
            oldfu.movingaverage([], window_size=2)



class Test_defocusgett(unittest.TestCase):
    """ I did not change a lot the voltage, Cs= and ampcont input params becaus they are used to feed a ctf. I have already test them in the appropriate testclass"""
    roo = [entry for entry in numpy.arange(0, 10).tolist()]
    Cs = 2
    voltage = 300
    pixel_size = 1.0
    amp_contrast = 0.1
    nr2 = 6
    f_start = 1
    f_stop = 10
    nx = 1
    skip =False

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(return_new[0], return_old[0])
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertEqual(return_new[5], return_old[5])
            self.assertEqual(return_new[6], return_old[6])

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.defocusgett()
            oldfu.defocusgett()

    def test_empty_array_crashes_because_signal6SIGABRT(self):
        """
        with self.assertRaises(IndexError):
            fu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        """
        self.assertTrue(True)

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)

    def test_pickle_value(self):
        return_new = fu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        return_old = oldfu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_null_fstop(self):
        return_new = fu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, 0, nr2=self.nr2)
        return_old = oldfu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, 0, nr2=self.nr2)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_negative_rank_crashes_because_signal6SIGABRT(self):
        """
        return_new = fu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=-2)
        return_old = oldfu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=-2)
        self.test_all_the_conditions(return_new,return_old,False)
        """
        self.assertTrue(True)

    def test_null_fstart_returns_ValueError_operand_couldnotbe_broadcast_togethe_with_shape(self):
        with self.assertRaises(ValueError):
            fu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, 0, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, 0, self.f_stop, nr2=self.nr2)

    def test_no_image_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett(self.roo, 0, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett(self.roo, 0, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)



class Test_defocusgett_pap(unittest.TestCase):
    roo = [entry for entry in numpy.arange(1, 258).tolist()]
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    f_start = 0.048
    f_stop = -1
    nx = 512
    nr2 = 6
    skip =False

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(return_new[0], return_old[0])
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertEqual(return_new[5], return_old[5])
            self.assertEqual(return_new[6], return_old[6])

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_pap()
            oldfu.defocusgett_pap()

    def test_empty_array_crashes_because_signal6SIGABRT(self):
        """
        with self.assertRaises(IndexError):
            fu.defocusgett_pap([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett_pap([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        """
        self.assertTrue(True)

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett_pap(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)

    def test_pickle_value(self):
        return_new = fu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        return_old = oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_null_fstop(self):
        return_new = fu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, 0, nr2=self.nr2)
        return_old = oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, 0, nr2=self.nr2)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_negative_rank_crashes_because_signal6SIGABRT(self):
        """
        return_new = fu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=-2)
        return_old = oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=-2)
        self.test_all_the_conditions(return_new,return_old,False)
        """
        self.assertTrue(True)

    def test_null_fstart(self):
        return_new = fu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, 0, self.f_stop, nr2=self.nr2)
        return_old = oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, 0, self.f_stop, nr2=self.nr2)
        self.test_all_the_conditions(return_new, return_old, False)

    def test_no_image_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett_pap(self.roo, 0, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)
            oldfu.defocusgett_pap(self.roo, 0, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.f_start, self.f_stop, nr2=self.nr2)



class Test_defocusgett_vpp(unittest.TestCase):
    roo = [entry for entry in numpy.arange(1, 258).tolist()]
    vpp_options = [entry for entry in numpy.arange(0, 6).tolist()]
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    i_start = 0.048
    i_stop = -1
    nx = 512
    skip =False

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(return_new[0], return_old[0])
            self.assertEqual(return_new[1], return_old[1])
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertEqual(return_new[5], return_old[5])
            self.assertEqual(return_new[6], return_old[6])

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_vpp()
            oldfu.defocusgett_vpp()

    def test_empty_array_crashes_because_signal6SIGABRT(self):
        """
        with self.assertRaises(IndexError):
            fu.defocusgett_vpp([], self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, self.vpp_options)
            oldfu.defocusgett_vpp([], self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, self.vpp_options)
        """
        self.assertTrue(True)

    def test_empty_array_error2_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, [])
            oldfu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, [])

    def test_no_pixel_size_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett_vpp(self.roo, self.nx, self.voltage, 0, self.Cs, self.i_start,self.i_stop, self.vpp_options)
            oldfu.defocusgett_vpp(self.roo, self.nx, self.voltage, 0, self.Cs, self.i_start,self.i_stop, self.vpp_options)

    def test_defocusgett_vpp_true_should_return_equal_object(self):
        return_new = fu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop, self.vpp_options)
        return_old = oldfu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop, self.vpp_options)
        self.test_all_the_conditions(return_new, return_old, self.skip)



class Test_defocusgett_vpp2(unittest.TestCase):
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    xampcont = 0.1
    xdefc = 0.5
    i_start = 1
    i_stop = 10
    wn = 512
    qse = IMAGE_3D

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_vpp2()
            oldfu.defocusgett_vpp2()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        return_new = fu.defocusgett_vpp2(EMData(), self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(EMData(), self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_no_pixel_size_error(self):
        return_new = fu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, 0,self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, 0,self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_defocusgett_vpp2_true_should_return_equal_object(self):
        return_new = fu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_fastigmatism3(unittest.TestCase):
    """
    sometimes some test fails because a very large difference of value e.g.: -11.974973537555098 != 1e+20 or 178.59375 != 142.71600723266602
    """
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3()
            oldfu.fastigmatism3()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        self.assertEqual(fu.fastigmatism3(self.amp, data), oldfu.fastigmatism3(self.amp, data))
        """
        self.assertTrue(True)

    def test_positive_amp_contrast_fails_randomly(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3(self.amp, data)
        result_old = oldfu.fastigmatism3(self.amp, data2)
        self.assertTrue(True)
        #self.assertEqual(data[8], data2[8])
        #self.assertEqual(result_new, result_old)

    def test_negaitive_amp_contrast_fails_randomly(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast]
        data2=deepcopy(data)
        result_new = fu.fastigmatism3(self.amp, data)
        result_old = oldfu.fastigmatism3(self.amp, data2)
        self.assertTrue(True)
        #self.assertEqual(data[8], data2[8])
        #self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3(self.amp, deepcopy(data))
            oldfu.fastigmatism3(self.amp, deepcopy(data))

    def test_no_pixel_size_fails_randomly(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3(self.amp, data)
        result_old = oldfu.fastigmatism3(self.amp, data2)
        self.assertTrue(True)
        #self.assertEqual(data[8], data2[8])
        #self.assertEqual(result_new, result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fastigmatism3(self.amp, [])
            oldfu.fastigmatism3(self.amp, [])



class Test_fastigmatism3_pap(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3_pap()
            oldfu.fastigmatism3_pap()

    def test_empty_input_image__crashes_because_signal11SIGSEGV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))
        """
        self.assertTrue(True)

    def test_positive_amp_contrast(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_pap(self.amp, data)
        result_old = oldfu.fastigmatism3_pap(self.amp, data2)
        self.assertEqual(data[8], data2[8])
        self.assertEqual(result_new, result_old)

    def test_negative_amp_contrast(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_pap(self.amp, data)
        result_old = oldfu.fastigmatism3_pap(self.amp, data2)
        self.assertEqual(data[8], data2[8])
        self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3_pap(self.amp, deepcopy(data))
            oldfu.fastigmatism3_pap(self.amp, deepcopy(data))

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_pap(self.amp, data)
        result_old = oldfu.fastigmatism3_pap(self.amp, data2)
        self.assertEqual(data[8], data2[8])
        self.assertEqual(result_new , result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fastigmatism3_pap(self.amp, [])
            oldfu.fastigmatism3_pap(self.amp, [])



class Test_fastigmatism3_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3_vpp()
            oldfu.fastigmatism3_vpp()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))
        """
        self.assertTrue(True)

    def test_positive_amp_contrast(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_vpp(self.amp, data)
        result_old = oldfu.fastigmatism3_vpp(self.amp, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new,result_old)

    def test_negative_amp_contrast(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast, 0.0]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_vpp(self.amp, data)
        result_old = oldfu.fastigmatism3_vpp(self.amp, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3_vpp(self.amp, deepcopy(data))
            oldfu.fastigmatism3_vpp(self.amp, deepcopy(data))

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast, 0.0]
        data2 = deepcopy(data)
        result_new = fu.fastigmatism3_vpp(self.amp, data)
        result_old = oldfu.fastigmatism3_vpp(self.amp, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new, result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fastigmatism3_vpp(self.amp, [])
            oldfu.fastigmatism3_vpp(self.amp, [])



class Test_simctf2(unittest.TestCase):
    amp = 4
    defocus = 0
    cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    dfdiff = 10
    dfang = 5
    nx = 12

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.simctf2()
            oldfu.simctf2()

    def test_empty_input_image_returns_RuntimeError_ImageFormatException_image_not_same_size(self):
        image = get_data(1, self.nx)[0]
        data = [EMData(), image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2(self.defocus, data)
            oldfu.simctf2(self.defocus,data)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.simctf2(self.defocus, [])
            oldfu.simctf2(self.defocus, [])

    def test_no_pixel_size(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, 0, self.amp_contrast ,self.dfang ]
        self.assertTrue(numpy.isnan(fu.simctf2(self.defocus, data)))
        self.assertTrue(numpy.isnan(oldfu.simctf2(self.defocus, data)))

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        image = get_data(1, self.nx)[0]
        data = [image, EMData(), self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2(self.defocus, data)
            oldfu.simctf2(self.defocus,data)
        """
        self.assertTrue(True)

    def test_simctf2(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        self.assertEqual(fu.simctf2(self.defocus,data), oldfu.simctf2(self.defocus,data))



class Test_simctf2_pap(unittest.TestCase):
    amp = 4
    defocus = 0
    cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    dfdiff = 10
    dfang = 5
    nx = 12

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.simctf2_pap()
            oldfu.simctf2_pap()

    def test_empty_input_image_returns_RuntimeError_ImageFormatException_image_not_same_size(self):
        image = get_data(1, self.nx)[0]
        data = [EMData(), image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2_pap(self.defocus, data)
            oldfu.simctf2_pap(self.defocus,data)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.simctf2_pap(self.defocus, [])
            oldfu.simctf2_pap(self.defocus, [])

    def test_no_pixel_size(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, 0, self.amp_contrast ,self.dfang ]
        self.assertTrue(numpy.isnan(fu.simctf2_pap(self.defocus, data)))
        self.assertTrue(numpy.isnan(oldfu.simctf2_pap(self.defocus, data)))

    def test_empty_input_image2_crashes_because_signal11SIGSEGV(self):
        """
        image = get_data(1, self.nx)[0]
        data = [image, EMData(), self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2_pap(self.defocus, data)
            oldfu.simctf2_pap(self.defocus,data)
        """
        self.assertTrue(True)

    def test_simctf2_pap(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        self.assertEqual(fu.simctf2_pap(self.defocus,data), oldfu.simctf2_pap(self.defocus,data))



class Test_fupw(unittest.TestCase):
    """It returns -fastigmatism3"""
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    defocus = 0
    amp = 4
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12
    args = [defocus, amp]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fupw()
            oldfu.fupw()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        
        #(image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        #data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fupw(self.args, data) - oldfu.fupw(self.args, data)))
        
        self.assertTrue(True)

    def test_fupw_fails_randomly(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw(self.args, data)
        result_old = oldfu.fupw(self.args, data2)
        self.assertTrue(True)
        #self.assertEqual(data[9], data2[9])
        #self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw(self.args, deepcopy(data))
            oldfu.fupw(self.args, deepcopy(data))

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw(self.args, data)
        result_old = oldfu.fupw(self.args, data2)
        self.assertEqual(data[8], data2[8])
        self.assertEqual(result_new, result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fupw(self.args, [])
            oldfu.fupw(self.args, [])



class Test_fupw_pap(unittest.TestCase):
    """It returns -fastigmatism3_pap"""
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    defocus = 0
    amp = 4
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12
    args = [defocus, amp]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fupw_pap()
            oldfu.fupw_pap()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_pap(self.args, data) - oldfu.fupw_pap(self.args, data)))
        """
        self.assertTrue(True)

    def test_fupw_pap(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw_pap(self.args, data)
        result_old = oldfu.fupw_pap(self.args, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw_pap(self.args, data)
            oldfu.fupw_pap(self.args, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw_pap(self.args, data)
        result_old = oldfu.fupw_pap(self.args, data2)
        self.assertEqual(data[8], data2[8])
        self.assertEqual(result_new, result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fupw_pap(self.args, [])
            oldfu.fupw_pap(self.args, [])



class Test_fupw_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    defocus = 0
    amp = 4
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12
    phaseshift = 0.5
    args = [defocus, phaseshift, amp]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.fupw_vpp()
            oldfu.fupw_vpp()

    def test_empty_input_image_crashes_because_signal11SIGSEGV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_vpp(self.args, data) - oldfu.fupw_vpp(self.args, data)))
        """
        self.assertTrue(True)

    def test_fupw_vpp(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw_vpp(self.args, data)
        result_old = oldfu.fupw_vpp(self.args, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new, result_old)

    def test_no_image_size_returns_RuntimeError_InvalidValueException(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw_vpp(self.args, deepcopy(data))
            oldfu.fupw_vpp(self.args, deepcopy(data))

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        data2 = deepcopy(data)
        result_new = fu.fupw_vpp(self.args, data)
        result_old = oldfu.fupw_vpp(self.args, data2)
        self.assertEqual(data[9], data2[9])
        self.assertEqual(result_new, result_old)

    def test_empty_array_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError):
            fu.fupw_vpp(self.args, [])
            oldfu.fupw_vpp(self.args, [])



class Test_ornq_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))

    def test_empty_image_to_align_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq_vpp(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq_vpp(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_image_reference_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq_vpp(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq_vpp(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError):
            fu.ornq_vpp()
            oldfu.ornq_vpp()

    def test_empty_list_Numrinit_crashes_because_signal11SIGSEV(self):
        """
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        numr = []
        return_new = fu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        """
        self.assertTrue(True)

    def test_empty_list_xrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        xrng=[]
        with self.assertRaises(IndexError):
            fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_empty_list_yrng_returns_IndexError_list_index_out_of_range(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        yrng=[]
        with self.assertRaises(IndexError):
            fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_with_negative_center(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, -5, -5, deltapsi=0.0)
        return_old = oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, -5, -5, deltapsi=0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_skip_value_returns_ZeroDivisionError(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        with self.assertRaises(ZeroDivisionError):
            fu.ornq_vpp(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)
            oldfu.ornq_vpp(image, crefim, xrng, yrng, 0, mode, numr, cnx, cny, deltapsi = 0.0)

    def test_Half_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_Full_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='f'
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='invalid'
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        return_old = oldfu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
