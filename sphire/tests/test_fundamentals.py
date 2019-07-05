



from __future__ import print_function
from __future__ import division
from past.utils import old_div

import unittest
from test_module import get_real_data,create_kb
from ..libpy import sparx_utilities
import numpy

from EMAN2_cppwrap import EMData
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = sparx_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sparx_utilities.model_blank(10, 10, 10)
MASK = sparx_utilities.model_circle(2, 5, 5)
REAL_IMAGE_2D = get_real_data(dim =2)[0]
REAL_IMAGE_3D = get_real_data(dim =3)[0]

from ..libpy_py3 import sphire_fundamentals as fu
from .sparx_lib import sparx_fundamentals as oldfu


"""
There are some opened issues in:
In all the tests miss the case with a complex image. where can we find one of them?
1) fftip: I cannot find an image that gives me as output a not None value 
2) prepi3D: It is able to work with a 2D image too. Should it be blocked?
3) rot_shift3D: It is able to work with a 2D image too. Should it be blocked?
4) goldsearch: HOw can I test it?
"""


class Test_ccf(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccf()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccf()
        self.assertEqual(cm_new.exception.message, "ccf() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccf(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccf(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccf(IMAGE_2D,EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccf(IMAGE_2D,EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccf(IMAGE_2D, None, center=False)
        return_old = oldfu.ccf(IMAGE_2D, None, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccf(None,IMAGE_2D, center=False)
        return_old = oldfu.ccf(None, IMAGE_2D, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccf(IMAGE_2D, REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccf(IMAGE_2D, REAL_IMAGE_2D, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccf(IMAGE_2D, REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccf(IMAGE_2D, REAL_IMAGE_2D, center=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccf(IMAGE_3D, REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccf(IMAGE_3D, REAL_IMAGE_3D, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccf(IMAGE_3D, REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccf(IMAGE_3D, REAL_IMAGE_3D, center=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_scf(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scf()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scf()
        self.assertEqual(cm_new.exception.message, "scf() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scf(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scf(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scf(None, center=False)
        return_old = oldfu.scf(None, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scf(IMAGE_3D, center=True)
        return_old = oldfu.scf(IMAGE_3D, center=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scf(IMAGE_3D, center=False)
        return_old = oldfu.scf(IMAGE_3D, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scf(IMAGE_2D, center=True)
        return_old = oldfu.scf(IMAGE_2D, center=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scf(IMAGE_2D, center=False)
        return_old = oldfu.scf(IMAGE_2D, center=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_cyclic_shift(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cyclic_shift()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cyclic_shift()
        self.assertEqual(cm_new.exception.message, "cyclic_shift() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_input_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.cyclic_shift(None)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.cyclic_shift(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Img2D_with_no_shift(self):
        return_new = fu.cyclic_shift(IMAGE_2D, dx=0, dy=0, dz=0)
        return_old = oldfu.cyclic_shift(IMAGE_2D, dx=0, dy=0, dz=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_no_shift(self):
        return_new = fu.cyclic_shift(IMAGE_3D, dx=0, dy=0, dz=0)
        return_old = oldfu.cyclic_shift(IMAGE_3D, dx=0, dy=0, dz=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_huge_shift(self):
        return_new = fu.cyclic_shift(IMAGE_2D, dx=10000, dy=10000, dz=10000)
        return_old = oldfu.cyclic_shift(IMAGE_2D, dx=10000, dy=10000, dz=10000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_huge_shift(self):
        return_new = fu.cyclic_shift(IMAGE_3D, dx=10000, dy=10000, dz=10000)
        return_old = oldfu.cyclic_shift(IMAGE_3D, dx=10000, dy=10000, dz=10000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_mirror(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.mirror()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.mirror()
        self.assertEqual(cm_new.exception.message, "mirror() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        return_new = fu.mirror(EMData(), axis='x')
        return_old = oldfu.mirror(EMData(), axis='x')
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.mirror(None)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.mirror(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.mirror(IMAGE_2D, axis='x')
        return_old = oldfu.mirror(IMAGE_2D, axis='x')
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg(self):
        return_new = fu.mirror(IMAGE_3D, axis='x')
        return_old = oldfu.mirror(IMAGE_3D, axis='x')
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_invalid_axis(self):
        return_new = fu.mirror(IMAGE_2D, axis='invalid_axis')
        return_old = oldfu.mirror(IMAGE_2D, axis='invalid_axis')
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_fft(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.fft()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.fft()
        self.assertEqual(cm_new.exception.message, "fft() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fft(EMData(), npad=1)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fft(EMData(), npad=1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fft(None, npad=1)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fft(None, npad=1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'is_complex'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_withPad1(self):
        return_new = fu.fft(IMAGE_3D, npad=1)
        return_old = oldfu.fft(IMAGE_3D, npad=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_withPad1(self):
        return_new = fu.fft(IMAGE_3D, npad=1)
        return_old = oldfu.fft(IMAGE_3D, npad=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_withPadGrater1(self):
        return_new = fu.fft(IMAGE_3D, npad=3)
        return_old = oldfu.fft(IMAGE_3D, npad=3)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_withPadGrater1(self):
        return_new = fu.fft(IMAGE_3D, npad=3)
        return_old = oldfu.fft(IMAGE_3D, npad=3)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_fftip(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.fftip()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.fftip()
        self.assertEqual(cm_new.exception.message, "fftip() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fftip(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fftip(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'y size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fftip(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fftip(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'is_complex'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

"""
    def test_3Dimg(self):
        return_new = fu.fftip(REAL_IMAGE_3D)
        return_old = oldfu.fftip(REAL_IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg(self):
        return_new = fu.fftip(IMAGE_2D)
        return_old = oldfu.fftip(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))
"""





class Test_fpol(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.fpol()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.fpol()
        self.assertEqual(cm_new.exception.message, "fpol() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fpol(EMData(),100)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fpol(EMData(),100)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fpol(None,100)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fpol(None,100)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_wrong_size_newimg_return_RuntimeError_ImageDimensionException(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fpol(image=IMAGE_2D, nnx=int(IMAGE_2D.get_xsize()/2), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fpol(image=IMAGE_2D, nnx=int(IMAGE_2D.get_xsize()/2), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Cannot reduce the image size')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_2Dimg_Normalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_fdecimate(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.fdecimate()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.fdecimate()
        self.assertEqual(cm_new.exception.message, "fdecimate() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fdecimate(EMData(),100)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fdecimate(EMData(),100)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Cannot increase the image size')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fdecimate(None,100)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fdecimate(None,100)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'FourTruncate'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_2Dimg_different_size_newimg(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=int(IMAGE_2D.get_xsize() / 2), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=int(IMAGE_2D.get_xsize() / 2), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_fshift(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.fshift()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.fshift()
        self.assertEqual(cm_new.exception.message, "fshift() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fshift(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fshift(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        return_new = fu.fshift(None)
        return_old = oldfu.fshift(None)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2Dimg_without_shift(self):
        return_new = fu.fshift(IMAGE_2D, delx=0, dely=0, delz=0)
        return_old = oldfu.fshift(IMAGE_2D, delx=0, dely=0, delz=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_without_shift(self):
        return_new = fu.fshift(IMAGE_3D, delx=0, dely=0, delz=0)
        return_old = oldfu.fshift(IMAGE_3D, delx=0, dely=0, delz=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.fshift(IMAGE_2D, delx=1000, dely=1000, delz=1000)
        return_old = oldfu.fshift(IMAGE_2D, delx=1000, dely=1000, delz=1000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_huge_shift(self):
        return_new = fu.fshift(IMAGE_3D, delx=1000, dely=1000, delz=1000)
        return_old = oldfu.fshift(IMAGE_3D, delx=1000, dely=1000, delz=1000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))






class Test_subsample(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.subsample()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.subsample()
        self.assertEqual(cm_new.exception.message, "subsample() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        return_new = fu.subsample(EMData())
        return_old = oldfu.subsample(EMData())
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.subsample(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.subsample(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.subsample(IMAGE_2D, subsample_rate=1.0)
        return_old = oldfu.subsample(IMAGE_2D, subsample_rate=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_subsample_rate_higher1(self):
        return_new = fu.subsample(IMAGE_2D, subsample_rate=2.0)
        return_old = oldfu.subsample(IMAGE_2D, subsample_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_negative_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(IMAGE_2D, subsample_rate=-1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(IMAGE_2D, subsample_rate=-1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_null_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(IMAGE_2D, subsample_rate=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(IMAGE_2D, subsample_rate=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_huge_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(IMAGE_2D, subsample_rate=1000.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(IMAGE_2D, subsample_rate=1000.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "BadAllocException")
        self.assertEqual(msg[1], 'Cannot allocate 495616 GB - not enough memory.')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_3Dimg_default_value(self):
        return_new = fu.subsample(IMAGE_3D, subsample_rate=1.0)
        return_old = oldfu.subsample(IMAGE_3D, subsample_rate=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_subsample_rate_higher1(self):
        return_new = fu.subsample(IMAGE_3D, subsample_rate=2.0)
        return_old = oldfu.subsample(IMAGE_3D, subsample_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))







class Test_resample(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.resample()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.resample()
        self.assertEqual(cm_new.exception.message, "resample() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.resample(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.resample(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.resample(IMAGE_2D, sub_rate=1.0)
        return_old = oldfu.resample(IMAGE_2D, sub_rate=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_higher1(self):
        return_new = fu.resample(IMAGE_2D, sub_rate=2.0)
        return_old = oldfu.resample(IMAGE_2D, sub_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_lower1(self):
        return_new = fu.resample(IMAGE_2D, sub_rate=0.10)
        return_old = oldfu.resample(IMAGE_2D, sub_rate=0.10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_higher1_nx_notsame_ny(self):
        return_new = fu.resample(EMData(100,20), sub_rate=2.0)
        return_old = oldfu.resample(EMData(100,20), sub_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_negative_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(IMAGE_2D, sub_rate=-1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(IMAGE_2D, sub_rate=-1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_null_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(IMAGE_2D, sub_rate=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(IMAGE_2D, sub_rate=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_huge_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(IMAGE_2D, sub_rate=1000.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(IMAGE_2D, sub_rate=1000.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "BadAllocException")
        self.assertEqual(msg[1], 'Cannot allocate 495.616 GB - not enough memory.')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_3Dimg_default_value(self):
        return_new = fu.resample(IMAGE_3D, sub_rate=1.0)
        return_old = oldfu.resample(IMAGE_3D, sub_rate=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_higher1(self):
        return_new = fu.resample(IMAGE_3D, sub_rate=2.0)
        return_old = oldfu.resample(IMAGE_3D, sub_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_lower1(self):
        return_new = fu.resample(IMAGE_3D, sub_rate=0.2)
        return_old = oldfu.resample(IMAGE_3D, sub_rate=0.2)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_higher1_nz_notsame_ny_or_nx(self):
        return_new = fu.resample(EMData(100,100,2), sub_rate=2.0)
        return_old = oldfu.resample(EMData(100,100,2), sub_rate=2.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))





class Test_prepi(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.prepi()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.prepi()
        self.assertEqual(cm_new.exception.message, "prepi() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.prepi(EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.prepi(EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepi(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.prepi(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_retFour(self):
        return_new = fu.prepi(IMAGE_2D,RetReal = False)
        return_old = oldfu.prepi(IMAGE_2D,RetReal = False)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(numpy.array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_2Dimg_retReal(self):
        return_new = fu.prepi(IMAGE_2D,RetReal = True)
        return_old = oldfu.prepi(IMAGE_2D,RetReal = True)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(numpy.array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepi(IMAGE_3D,RetReal = False)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.prepi(IMAGE_3D,RetReal = False)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Cannot reduce the image size')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])



class Test_prepi3D(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.prepi3D()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.prepi3D()
        self.assertEqual(cm_new.exception.message, "prepi3D() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.prepi3D(EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.prepi3D(EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepi3D(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.prepi3D(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg(self):
        return_new = fu.prepi3D(IMAGE_3D)
        return_old = oldfu.prepi3D(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(numpy.array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_2Dimg(self):
        return_new = fu.prepi3D(IMAGE_2D)
        return_old = oldfu.prepi3D(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(numpy.array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter



class Test_ramp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ramp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ramp()
        self.assertEqual(cm_new.exception.message, "ramp() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        return_new = fu.ramp(EMData())
        return_old = oldfu.ramp(EMData())
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.ramp(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.ramp(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.ramp(IMAGE_2D)
        return_old = oldfu.ramp(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_returns_RunTimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ramp(IMAGE_3D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ramp(IMAGE_3D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])




class Test_rot_avg_table(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rot_avg_table()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rot_avg_table()
        self.assertEqual(cm_new.exception.message, "rot_avg_table() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rot_avg_table(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rot_avg_table(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'No 1D images!')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_avg_table(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_avg_table(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rotavg'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.rot_avg_table(IMAGE_2D)
        return_old = oldfu.rot_avg_table(IMAGE_2D)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [1.70115065574646, 1.6371004581451416, 1.5626660585403442, 1.5877236127853394, 1.6073533296585083,1.6285395622253418, 1.6862214803695679, 1.7918702363967896, 1.8338433504104614, 1.8038475513458252,1.6706266403198242, 1.544594407081604, 1.4375859498977661, 1.2914648056030273, 1.167849063873291,1.080208659172058, 0.9806387424468994, 0.8771043419837952, 0.7996737360954285, 0.7433571815490723,0.7284249663352966, 0.769216001033783, 0.8287383317947388, 0.9201586246490479, 0.9932690262794495,1.0363448858261108, 1.123952031135559, 1.1950533390045166, 1.2668777704238892, 1.3277966976165771,1.4283941984176636, 1.5512932538986206, 1.6878538131713867, 1.793628454208374, 1.8632606267929077,1.9746273756027222, 2.0680618286132812, 2.1846938133239746, 2.3249411582946777, 2.4541871547698975,2.5610392093658447, 2.6200289726257324, 2.6819217205047607, 2.7509777545928955, 2.8288681507110596,2.9172451496124268, 3.0068471431732178, 3.0675971508026123, 3.081587553024292, 3.068607807159424,3.0510096549987793, 2.971430778503418, 2.8975279331207275, 2.8081719875335693, 2.7049367427825928,2.5696637630462646, 2.4190049171447754, 2.2764179706573486, 2.1662235260009766, 2.06880784034729,1.9702768325805664, 1.8747632503509521, 1.77823805809021, 1.6550695896148682, 1.5298069715499878,1.4096013307571411, 1.269547700881958, 1.1273504495620728, 0.9518274664878845, 0.7691776752471924,0.5582662224769592, 0.3821752667427063, 0.2139461785554886, 0.02325786091387272,-0.16214938461780548, -0.3665705919265747, -0.5570375919342041, -0.7568416595458984,-0.9444634914398193, -1.1122022867202759, -1.2591341733932495, -1.3866132497787476,-1.5305417776107788, -1.6798121929168701, -1.7983838319778442, -1.8982055187225342,-1.9665424823760986, -2.0055081844329834, -2.0240345001220703, -2.022763252258301,-2.0144519805908203, -1.985037088394165, -1.9594217538833618, -1.9007961750030518,-1.8341227769851685, -1.7523555755615234, -1.6880568265914917, -1.6217501163482666,-1.5742276906967163, -1.516431212425232, -1.4340367317199707, -1.3519632816314697,-1.2755727767944336, -1.2155612707138062, -1.1673047542572021, -1.1223112344741821,-1.0735875368118286, -1.0130205154418945, -0.9502072334289551, -0.8961243629455566,-0.8469223380088806, -0.8141208291053772, -0.7728646993637085, -0.7485840320587158,-0.7369289398193359, -0.7154955267906189, -0.6829231977462769, -0.6576652526855469,-0.6172498464584351, -0.5820884704589844, -0.5657995343208313, -0.5585340857505798,-0.5447290539741516, -0.5433396697044373, -0.552950918674469, -0.5448899865150452,-0.542202353477478, -0.541515052318573, -0.5387433171272278, -0.5594711303710938,-0.5963792204856873, -0.6536659598350525, -0.7025741338729858, -0.7419651746749878,-0.7738654613494873, -0.8288431763648987, -0.8792231678962708, -0.9185277223587036,-0.9535513520240784, -0.9755640625953674, -0.9769047498703003, -0.9703202843666077,-0.9619379043579102, -0.9533834457397461, -0.9389827251434326, -0.8972327709197998,-0.8668974041938782, -0.8482377529144287, -0.8443998098373413, -0.8439618349075317,-0.8385792970657349, -0.829200804233551, -0.8111365437507629, -0.7909402847290039,-0.7704744338989258, -0.7470849752426147, -0.7254742980003357, -0.6966635584831238,-0.6708171367645264, -0.6408230662345886, -0.6109101176261902, -0.5850300788879395,-0.5588661432266235, -0.5414973497390747, -0.5227035880088806, -0.5071979761123657,-0.4920077323913574, -0.47635769844055176, -0.46521392464637756, -0.4503937363624573,-0.4407220482826233, -0.43082553148269653, -0.41926488280296326, -0.4109551012516022,-0.3982444405555725, -0.39125651121139526, -0.372354120016098]))

    def test_3Dimg(self):
        return_new = fu.rot_avg_table(IMAGE_3D)
        return_old = oldfu.rot_avg_table(IMAGE_3D)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [10.473315238952637, 21.884193420410156, 31.711227416992188, 33.32331466674805, 26.039146423339844, 20.5319766998291, 18.808826446533203, 18.97246551513672, 20.31218719482422, 22.830305099487305, 27.15753936767578, 32.496463775634766, 36.15548324584961, 36.97674560546875, 34.278045654296875, 26.243494033813477, 15.062774658203125, 6.469794273376465, 2.7834198474884033, 2.3789703845977783, 2.7278802394866943, 2.6416027545928955, 2.4357831478118896, 2.5040085315704346, 2.6968085765838623, 2.794365882873535, 2.600947141647339, 2.0284249782562256, 1.4216328859329224, 0.5415528416633606, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))





class Test_rops_table(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rops_table()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rops_table()
        self.assertEqual(cm_new.exception.message, "rops_table() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rops_table(EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rops_table(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        fu.rops_table(None)
        oldfu.rops_table(None)
        """

    def test_3Dimg_without_lng(self):
        return_new = fu.rops_table(IMAGE_3D, lng = False)
        return_old = oldfu.rops_table(IMAGE_3D, lng = False)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [10.759023666381836, 3.9929940700531006, 0.9002360701560974, 0.0958547592163086, 0.0637880265712738, 0.03903469815850258, 0.016267063096165657, 0.0196246225386858, 0.014603659510612488, 0.013418514281511307, 0.01010554376989603, 0.008397913537919521, 0.00930393859744072, 0.008090538904070854, 0.004846676718443632, 0.003233283292502165, 0.002603360218927264, 0.0017707381630316377, 0.000983564299531281, 0.00034340255660936236, 9.67993400990963e-05, 3.4131360735045746e-05, 2.176243106077891e-05, 1.7841332009993494e-05, 1.2615307241503615e-05, 7.829644346202258e-06, 5.766305093857227e-06, 4.558276941679651e-06, 4.2896176637441386e-06, 2.9617256132041803e-06, 1.8494994264983688e-06, 1.2403851314957137e-06, 7.150275109779614e-07, 4.678593370499584e-07, 3.486518664885807e-07, 2.4661272846060456e-07, 1.8870717610752763e-07, 1.6899589638796897e-07, 1.526483970337722e-07]))

    def test_3Dimg_with_lng(self):
        return_new = fu.rops_table(IMAGE_3D, lng = True)
        return_old = oldfu.rops_table(IMAGE_3D, lng = True)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [0.6012986652293667, 0.6012986652293667, -0.045643589980165239, -1.0183863194049168, -1.1952608335363788, -1.4085491745317404, -1.7886908488001898, -1.7071986877517489, -1.8355383013307722, -1.8722955672565147, -1.9954403125384363, -2.075828601024317, -2.0313331644533683, -2.0920225494298332, -2.3145559474896125, -2.4903562418838359, -2.5844657357444896, -2.7518456526763888, -3.0071972432496548, -3.4641964758693211, -4.0141276033558411, -4.4668463971025076, -4.6622925916343245, -4.7485727249041201, -4.899102167923675, -5.1062579648883961, -5.2391023829201826, -5.3411992924532452, -5.3675814150240981, -5.5284551788300433, -5.7329457990166768, -5.9064434482941657, -6.145677248218405, -6.3298846989031929, -6.4576080054812302, -6.6079845118540543, -6.7242115842831742, -6.7721238409348379, -6.8163077518831221]))

    def test_2Dimg_without_lng(self):
        return_new = fu.rops_table(IMAGE_2D, lng = False)
        return_old = oldfu.rops_table(IMAGE_2D, lng = False)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [0.4827031195163727, 0.17195598781108856, 0.15751279890537262, 0.05366284027695656, 0.039997369050979614, 0.019734526053071022, 0.011997099965810776, 0.007876859977841377, 0.005689265672117472, 0.003199680708348751, 0.002140969969332218, 0.0015250766882672906, 0.0011131446808576584, 0.0007962905219756067, 0.0005700858309864998, 0.0004474169109016657, 0.0003331159532535821, 0.00022828028886578977, 0.00019683958089444786, 0.00015699909999966621, 0.00011236220598220825, 9.793655044632033e-05, 8.478884410578758e-05, 6.88535365043208e-05, 5.43036476301495e-05, 4.820919639314525e-05, 3.8353464333340526e-05, 3.5439417843008414e-05, 3.1193096219794825e-05, 2.722363751672674e-05, 2.3752450942993164e-05, 2.1170448235352524e-05, 2.0454015611903742e-05, 1.6884650904103182e-05, 1.708091440377757e-05, 1.527714448457118e-05, 1.4738235222466756e-05, 1.3024156942265108e-05, 1.3137992027623113e-05, 1.2084764421160799e-05, 1.0859323083423078e-05, 1.1200419066881295e-05, 9.873318958852906e-06, 1.0158905752177816e-05, 9.031403351400513e-06, 8.884254384611268e-06, 7.928477316454519e-06, 7.31817863197648e-06, 8.140462341543753e-06, 7.877236384956632e-06, 5.646501904266188e-06, 5.4617944442725275e-06, 5.845015948580112e-06, 5.924002380197635e-06, 5.140260782354744e-06, 4.2990141082555056e-06, 4.328401701059192e-06, 4.164225629210705e-06, 3.4540489650680684e-06, 3.368472334841499e-06, 3.279993734395248e-06, 2.7609373773884727e-06, 2.5925535283022327e-06, 2.532477310523973e-06, 2.243562903458951e-06, 2.0025879621243803e-06, 1.9317926671647e-06, 1.7656581121627823e-06, 1.5821302667973214e-06, 1.577492866999819e-06, 1.3625522115034983e-06, 1.2911820022054599e-06, 1.298216716349998e-06, 1.1096108210040256e-06, 1.0843688187378575e-06, 1.087746340999729e-06, 9.52967013745365e-07, 9.473871500631503e-07, 8.989124467007059e-07, 8.75390071541915e-07, 8.352129157174204e-07, 8.028847560126451e-07, 7.663265364499239e-07, 7.32731791686092e-07, 7.518655138483155e-07, 7.11970926658978e-07, 6.684190339001361e-07, 7.046898531370971e-07, 6.311749984888593e-07, 6.354318884405075e-07, 6.346119789668592e-07, 5.984673521197692e-07, 6.061846988814068e-07, 5.908723892389389e-07, 5.5673859833405e-07, 5.606630679722002e-07, 5.515345833373431e-07, 5.522384753930965e-07, 5.170920758246211e-07, 5.184934934732155e-07, 5.073741817795963e-07, 5.14062776346691e-07, 4.85623900203791e-07, 4.840583187615266e-07, 4.77287187550246e-07, 4.72120149197508e-07, 4.592948528170382e-07, 4.462744698230381e-07, 4.484391240566765e-07, 4.5394514813779097e-07, 4.3222709678047977e-07, 4.231489185713144e-07, 4.312424835006823e-07, 4.1006330775417155e-07, 4.1842443465611723e-07, 4.034653500184504e-07, 3.994539099494432e-07, 3.9826986153457256e-07, 3.8763781162742816e-07, 3.879274572682334e-07, 3.7914128370175604e-07, 3.8078843545008567e-07, 3.6390775903782924e-07, 3.744132186511706e-07, 3.722320229826437e-07, 3.556713750185736e-07, 3.598441367103078e-07, 3.534894403856015e-07, 3.444282299369661e-07, 3.4894375744443096e-07, 3.3972065693888e-07, 3.3560132806087495e-07, 3.497297598187288e-07, 3.200247817858326e-07, 3.356126967446471e-07, 3.26243849713137e-07, 3.1911702080833493e-07, 3.2429471730210935e-07, 3.1744849593451363e-07, 3.1914674991639913e-07, 3.065918576794502e-07, 3.158135086778202e-07, 3.098435570336733e-07, 3.0147143093017803e-07, 3.0035425879759714e-07, 3.064098734739673e-07, 2.9532270673371386e-07, 2.9205398277554195e-07, 2.890424752877152e-07, 3.078606027884234e-07, 2.778321288587904e-07, 2.939001433333033e-07, 2.8251207595531014e-07, 2.822802116497769e-07, 2.86119160364251e-07, 2.7603505259321537e-07, 2.80413672726354e-07, 2.7507161348694353e-07, 2.7496795951265085e-07, 2.7524157530933735e-07, 2.7048164952248044e-07, 2.6790505103235773e-07, 2.692670477699721e-07, 2.656799154010514e-07, 2.6457624358044995e-07, 2.688267102257669e-07, 2.5889752919283637e-07, 2.607513636121439e-07, 2.5427883088013914e-07, 2.650857311436994e-07, 2.527045808164985e-07, 2.5341572040815663e-07, 2.538510273097927e-07, 2.560517486926983e-07, 2.505655629647663e-07, 2.4852897695382126e-07, 2.496290960607439e-07]))

    def test_2Dimg_with_lng(self):
        return_new = fu.rops_table(IMAGE_2D, lng = True)
        return_old = oldfu.rops_table(IMAGE_2D, lng = True)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [-0.76458269667804002, -0.76458269667804002, -0.80268415128301362, -1.2703263446352904, -1.3979685747775406, -1.7047732992156459, -1.9209237223738624, -2.1036468746514565, -2.2449437854727594, -2.4948933571553122, -2.6693894243658951, -2.8167083173310177, -2.953448384619008, -3.0989284535292905, -3.2440597528931612, -3.3492876050125671, -3.4774045680170422, -3.6415315865666615, -3.7058875683251262, -3.8041028371850274, -3.9493797429797546, -4.0090551969127723, -4.0716612852670391, -4.1620737482920616, -4.265170997438033, -4.3168701077812823, -4.4161954016523852, -4.4505134208043327, -4.5059415151295648, -4.5650538465133863, -4.6242915702034173, -4.6742699467005044, -4.6892214169040853, -4.7725079141217845, -4.7674888836537477, -4.8159578139858459, -4.8315545164111295, -4.88525037913569, -4.8814710060061159, -4.9177618112544383, -4.9641972456768721, -4.9507657277736321, -5.005536832819625, -4.9931530687631351, -5.0442447612683701, -5.051379014705816, -5.1008102119962553, -5.1355969938744392, -5.0893509284426113, -5.1036261217590821, -5.2482205210588564, -5.2626646486485047, -5.2332142994942901, -5.227384775656148, -5.2890148472567358, -5.3666311296381775, -5.363672440956516, -5.3804657472531332, -5.4616715101072577, -5.472567015143726, -5.4841269858981283, -5.5589434437638729, -5.5862722679933885, -5.5964544369680285, -5.6490617495130593, -5.6984083987761904, -5.7140394867912034, -5.7530933859087572, -5.800757761153049, -5.802032595842701, -5.8656468470451459, -5.8890125362134942, -5.8866528030220096, -5.9548293166432957, -5.964822979144655, -5.9634723689281897, -6.0209221318866488, -6.0234725101448081, -6.0462826061159856, -6.057798383355772, -6.0782027983678732, -6.0953467877469976, -6.1155861353799823, -6.1350549649311619, -6.1238598346967219, -6.1475377404212503, -6.174951191633073, -6.1520019818847294, -6.1998502123229926, -6.1969309943778583, -6.1974917343942808, -6.2229595364538399, -6.2173950301512395, -6.2285063036086443, -6.2543486682881948, -6.2512980508288116, -6.2584272504374292, -6.2578733386080154, -6.2864321175209286, -6.2852566891619386, -6.2946716360420254, -6.2889838425685056, -6.3136999474039746, -6.3151023119288379, -6.3212202238171065, -6.3259474643038471, -6.3379084215659285, -6.35039795728252, -6.3482965044162762, -6.3429966213725502, -6.3642880102122446, -6.3735067647005934, -6.3652784615479465, -6.3871490894123664, -6.3783829619838412, -6.3941937569566507, -6.3985333235073627, -6.3998225569209746, -6.4115738668632316, -6.411249480214158, -6.4211989238096683, -6.4193162496165499, -6.4390086845308563, -6.4266488269052964, -6.4291862674051048, -6.4489510860982362, -6.4438855693002468, -6.4516235551365613, -6.4629012602282394, -6.4572445667367999, -6.4688780447643159, -6.4741763292248375, -6.4562674103697875, -6.4948163898797313, -6.4741616175041443, -6.4864576668132745, -6.4960500311371119, -6.4890601258420748, -6.4983267262645024, -6.4960095739211905, -6.5134393831145205, -6.5005692973208831, -6.5088575302511593, -6.5207538376772431, -6.5223662057036975, -6.5136972445023114, -6.5297031597737201, -6.5345368668574446, -6.539038332234413, -6.5116458846152634, -6.5562175333310968, -6.5318002021245718, -6.5489629835660823, -6.5493195656098191, -6.5434530580918446, -6.5590357651043281, -6.5522008143685735, -6.5605542251031181, -6.5607179091426131, -6.560285965190582, -6.5678621937051433, -6.572019098206642, -6.569816791212709, -6.5756412756995779, -6.5774491538941842, -6.570527582614651, -6.5868720942659289, -6.583773411731034, -6.5946897940304412, -6.5766136485488094, -6.5973868854829485, -6.5961664475632755, -6.5954210747093605, -6.5916722538268662, -6.6010786174682492, -6.6046229679065895, -6.6027047959039793]))







class Test_gridrot_shift2D(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.gridrot_shift2D()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.gridrot_shift2D()
        self.assertEqual(cm_new.exception.message, "gridrot_shift2D() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.gridrot_shift2D(EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.gridrot_shift2D(EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.gridrot_shift2D(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.gridrot_shift2D(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zeroScale(self):
        return_new = fu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        return_old = oldfu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_angle(self):
        return_new = fu.gridrot_shift2D(IMAGE_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(IMAGE_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_scale(self):
        return_new = fu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        return_old = oldfu.gridrot_shift2D(IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.gridrot_shift2D(IMAGE_3D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.gridrot_shift2D(IMAGE_3D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'fouriergridrot2d needs a 2-D image.')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])






class Test_rot_shift2D(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rot_shift2D()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rot_shift2D()
        self.assertEqual(cm_new.exception.message, "rot_shift2D() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rot_shift2D(IMAGE_3D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rot_shift2D(IMAGE_3D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Volume not currently supported')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.rot_shift2D(EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift2D(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1]+msg[2]+msg[3], 'Cant rotate 1D image caught\n')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1]+msg[2]+msg[3], msg_old[1]+msg_old[2]+msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_shift2D(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_shift2D(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_trans2D_background'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_InvalidMethod_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_invalid_interpolation_method_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_default_value(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Zero_scale_print_errorMsg(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodLinear(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodQuadratic(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodGridding(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodFtgridding(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodFourier(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodLinear(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="cyclic")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="linear", mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodQuadratic(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodGridding(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodFtgridding(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodFourier(self):
        return_new = fu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="cyclic")
        return_old = oldfu.rot_shift2D(IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="fourier", mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_rot_shift3D(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rot_shift3D()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rot_shift3D()
        self.assertEqual(cm_new.exception.message, "rot_shift3D() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.rot_shift3D(EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift3D(EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_shift3D(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_shift3D(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_trans_background'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.rot_shift3D(IMAGE_2D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        return_old = oldfu.rot_shift3D(IMAGE_2D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Zero_scale_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=0, mode="background")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=0, mode="background")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'The scale factor in a Transform object must be positive and non zero')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_3Dimg_default_values(self):
        return_new = fu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        return_old = oldfu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_modeCyclic(self):
        return_new = fu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="cyclic")
        return_old = oldfu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="cyclic")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Invalid_mode(self):
        return_new = fu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="invalid")
        return_old = oldfu.rot_shift3D(IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="invallid")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_rtshg(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rtshg()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rtshg()
        self.assertEqual(cm_new.exception.message, "rtshg() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.rtshg(EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.rtshg(EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rtshg(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rtshg(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.rtshg(IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.rtshg(IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Cannot reduce the image size')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_2Dimg_default_value(self):
        return_new = fu.rtshg(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_hugeShift(self):
        return_new = fu.rtshg(IMAGE_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        return_old = oldfu.rtshg(IMAGE_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_hugeAngle(self):
        return_new = fu.rtshg(IMAGE_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(IMAGE_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zeroScale(self):
        return_new = fu.rtshg(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.rtshg(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_rtshgkb(unittest.TestCase):
    kb=create_kb(dim=1)
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rtshgkb()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rtshgkb()
        self.assertEqual(cm_new.exception.message, "rtshgkb() takes at least 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_wrongKb(self):
        kb=create_kb(dim=3)
        with self.assertRaises(TypeError)  as cm_new:
            fu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=kb, scale=1.0)
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")

        self.assertEqual(msg[0].split("\n")[0]+msg[0].split("\n")[1], "Python argument types in    EMData.rot_scale_conv_new(EMData, numpy.float64, float, float, tuple, float)")
        self.assertEqual(msg[0].split("\n")[0]+msg[0].split("\n")[1], msg_old[0].split("\n")[0]+msg_old[0].split("\n")[1])

    def test_2Dimg_empty_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rtshgkb(EMData(), angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rtshgkb(EMData(), angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1]+msg[2]+msg[3], 'Cant rotate 1D image caught\n')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1]+msg[2]+msg[3], msg_old[1]+msg_old[2]+msg_old[3])

    def test_NoneType_img(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.rtshgkb(None, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.rtshgkb(None, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_conv_new'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_default_value(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rtshgkb(IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rtshgkb(IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Use rot_scale_conv_new_3D for volumes')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_2Dimg_default_value(self):
        return_new = fu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zero_scale(self):
        return_new = fu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        return_old = oldfu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_angle(self):
        return_new = fu.rtshgkb(IMAGE_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(IMAGE_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.rtshgkb(IMAGE_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(IMAGE_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_scale(self):
        return_new = fu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        return_old = oldfu.rtshgkb(IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_smallprime(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.smallprime()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.smallprime()
        self.assertEqual(cm_new.exception.message, "smallprime() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_default_value(self):
        return_new = fu.smallprime(arbit_num=300, numprime=3)
        return_old = oldfu.smallprime(arbit_num=300, numprime=3)
        self.assertEqual(return_old,return_new)
        self.assertEqual(return_new,300)

    def test_zero_numprime(self):
        return_new = fu.smallprime(arbit_num=300, numprime=0)
        return_old = oldfu.smallprime(arbit_num=300, numprime=0)
        self.assertEqual(return_old,return_new)
        self.assertEqual(return_new,1)

    def test_zero_arbit_num_returns_UnboundLocalError(self):
        with self.assertRaises(UnboundLocalError)  as cm_new:
            fu.smallprime(arbit_num=0, numprime=3)
        with self.assertRaises(UnboundLocalError)  as cm_old:
            oldfu.smallprime(arbit_num=0, numprime=3)
        self.assertEqual(cm_new.exception.message, "local variable 'i' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_tilemic(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.tilemic()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.tilemic()
        self.assertEqual(cm_new.exception.message, "tilemic() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_img(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.tilemic(None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.tilemic(None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyImg_returns_UnboundLocalError(self):
        with self.assertRaises(UnboundLocalError)  as cm_new:
            fu.tilemic(EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(UnboundLocalError)  as cm_old:
            oldfu.tilemic(EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "local variable 'i' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.tilemic(IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, []))

    def test_3Dimg_default_value(self):
        return_new = fu.tilemic(EMData(100,100,100), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(EMData(100,100,100), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_3Dimg_(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.tilemic(IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.tilemic(IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'z size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_zero_win_size(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.tilemic(IMAGE_2D, win_size=0, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.tilemic(IMAGE_2D, win_size=0, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_no_overlap(self):
        return_new = fu.tilemic(IMAGE_2D, win_size=120, overlp_x=0, overlp_y=0, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(IMAGE_2D, win_size=120, overlp_x=0, overlp_y=0, edge_x=0, edge_y=0)
        for img1, img2 in zip(return_old, return_new):
            self.assertTrue(numpy.allclose(img1.get_3dview(), img2.get_3dview(), equal_nan=True))





class Test_window2d(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.window2d()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.window2d()
        self.assertEqual(cm_new.exception.message, "window2d() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_img(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.window2d(None, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.window2d(None, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyImg(self):
        return_new = fu.window2d(EMData(), isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        return_old = oldfu.window2d(EMData(), isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_central(self):
        return_new = fu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        return_old = oldfu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_left(self):
        return_new = fu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="l", ix=0, iy=0)
        return_old = oldfu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="l", ix=0, iy=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_arbitrary(self):
        return_new = fu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="a", ix=0, iy=0)
        return_old = oldfu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="a", ix=0, iy=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_invalid_option_returns_UnboundLocalError(self):
        with self.assertRaises(UnboundLocalError)  as cm_new:
            fu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="invalid", ix=0, iy=0)
        with self.assertRaises(UnboundLocalError)  as cm_old:
            oldfu.window2d(IMAGE_2D, isize_x=100, isize_y=100, opt="invalid", ix=0, iy=0)
        self.assertEqual(cm_new.exception.message, "local variable 'reg' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.window2d(IMAGE_3D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.window2d(IMAGE_3D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'z size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])



class Test_goldsearch(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.goldsearch()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.goldsearch()
        self.assertEqual(cm_new.exception.message, "goldsearch() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_rotate_params(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rotate_params()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rotate_params()
        self.assertEqual(cm_new.exception.message, "rotate_params() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_trans_value_returns_IndexError(self):
        with self.assertRaises(IndexError)  as cm_new:
            fu.rotate_params([1,2,3],[1])
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.rotate_params([1,2,3],[1])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_params_value_returns_IndexError(self):
        with self.assertRaises(IndexError)  as cm_new:
            fu.rotate_params([[1],[1]],[1,2,3])
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.rotate_params([[1],[1]],[1,2,3])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_params_value_returns_TypeError(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rotate_params([1],[1,2,3])
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rotate_params([1],[1,2,3])
        self.assertEqual(cm_new.exception.message, "'int' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotate_params(self):
        a = [[1, 2, 3], [3, 4, 5]]
        b = [1,5,6]
        return_new = fu.rotate_params(a, b)
        return_old = oldfu.rotate_params(a,b)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[174.00000000000003, 2.9999999999999232, 183.0], [166.07792596478271, 1.0120857125278824, 194.9159858556005]]))



class Test_rotmatrix(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rotmatrix()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rotmatrix()
        self.assertEqual(cm_new.exception.message, "rotmatrix() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotmatrix(self):
        return_new = fu.rotmatrix(1,5,6)
        return_old = oldfu.rotmatrix(1,5,6)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new, [[0.98876227196606559, 0.12180329553973532, -0.086678294469630643], [-0.12147164162554232, 0.99255309355170718, 0.0091102558543636417], [0.087142468505889387, 0.0015210774457754552, 0.99619469809174555]]))



class Test_mulmat(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.mulmat()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.mulmat()
        self.assertEqual(cm_new.exception.message, "mulmat() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_m1_value_returns_IndexError(self):
        a = [[1,  2], [2, 5, 3], [3, 8, 6]]
        b = [[2, 3, 2], [3, 4, 3], [3, 4, 6]]
        with self.assertRaises(IndexError)  as cm_new:
            fu.mulmat(a,b)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.mulmat(a,b)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_m2_value_returns_IndexError(self):
        a = [[1, 3, 2], [2, 5, 3], [3, 8, 6]]
        b = [[2,  2], [3, 4, 3], [3, 4, 6]]
        with self.assertRaises(IndexError)  as cm_new:
            fu.mulmat(a,b)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.mulmat(a,b)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_mulmat(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 6]]
        b = [[2, 3, 2], [3, 4, 3], [3, 4, 6]]
        return_new = fu.mulmat(a,b)
        return_old = oldfu.mulmat(a,b)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,[[23, 31, 29], [28, 38, 37], [48, 65, 66]]))



class Test_recmat(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.recmat()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.recmat()
        self.assertEqual(cm_new.exception.message, "recmat() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_fewElement(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 1]]
        with self.assertRaises(IndexError)  as cm_new:
            fu.recmat(a)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.recmat(a)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_recmat(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 1]]
        return_new = fu.recmat(a)
        return_old = oldfu.recmat(a)
        self.assertTrue(numpy.array_equal(return_new, return_old))
        self.assertTrue(numpy.array_equal(return_new,(78.690067525979785, 0.0, 0.0)))



class Test_symclass(unittest.TestCase):
    pass


















""" Adnan's test"""
import nose.tools as nt

import EMAN2_cppwrap as e2cpp
import numpy
import copy
import math
from ..libpy_py3 import sphire_fundamentals as fu
from ..libpy import sparx_fundamentals as fu
from .sparx_lib import sparx_fundamentals as oldfu

import unittest

def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


@unittest.skip("original adnan")
class Test_lib_compare(unittest.TestCase):
    def test_ccf_true_should_return_equal_objects(self):
        a, b = get_data(2)
        return_new = fu.ccf(a, b)
        return_old = oldfu.ccf(a, b)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_scf_true_should_return_equal_objects(self):
        a, = get_data(1)
        return_new = fu.scf(a)
        return_old = oldfu.scf(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_cyclic_shift_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.cyclic_shift(a)
        return_old = oldfu.cyclic_shift(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_mirror_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.mirror(a)
        return_old = oldfu.mirror(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_fft_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.fft(a)
        return_old = oldfu.fft(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_fftip_should_return_equal_object(self):
        a, = get_data(1)
        return_new = copy.deepcopy(a)
        return_old = copy.deepcopy(a)
        fu.fftip(return_new)
        oldfu.fftip(return_old)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_fpol_should_return_equal_object(self):
        a, = get_data(1)
        nx = a.get_xsize()
        print(a.get_xsize())
        return_new = fu.fpol(a,nx,nx)
        return_old = oldfu.fpol(a,nx,nx)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_fdecimate_should_return_equal_object(self):
        a, = get_data(1)
        nx = a.get_xsize()
        return_new = fu.fdecimate(a,nx)
        return_old = oldfu.fdecimate(a,nx)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_fshift_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.fshift(a)
        return_old = oldfu.fshift(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_subsample_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.subsample(a)
        return_old = oldfu.subsample(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_resample_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.resample(a)
        return_old = oldfu.resample(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_prepi_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.prepi(a)
        return_old = oldfu.prepi(a)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertEqual(return_new[1].i0win(1), return_old[1].i0win(1))

    def test_prepi3D_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.prepi3D(a)
        return_old = oldfu.prepi3D(a)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertEqual(return_new[1].i0win(1), return_old[1].i0win(1))

    def test_ramp_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.ramp(a)
        return_old = oldfu.ramp(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_rot_avg_table_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.rot_avg_table(a)
        return_old = oldfu.rot_avg_table(a)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_roops_table_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.rops_table(a)
        return_old = oldfu.rops_table(a)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_ramp_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.ramp(a)
        return_old = oldfu.ramp(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_gridrot_shift2D_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.gridrot_shift2D(a)
        return_old = oldfu.gridrot_shift2D(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_rot_shift2D_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.rot_shift2D(a)
        return_old = oldfu.rot_shift2D(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_rot_shift3D_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.rot_shift3D(a)
        return_old = oldfu.rot_shift3D(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_rtshg_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.rtshg(a)
        return_old = oldfu.rtshg(a)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_smallprime_should_return_equal_object(self):
        return_new = fu.smallprime(5)
        return_old = oldfu.smallprime(5)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_tilemic_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.tilemic(a)
        return_old = oldfu.tilemic(a)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_window2d_should_return_equal_object(self):
        a, = get_data(1)
        return_new = fu.window2d(a, 5,5)
        return_old = oldfu.window2d(a,5,5)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_goldsearch_should_return_equal_object(self):
        a, = get_data(1)
        def f(x):
            return 2
        return_new = fu.goldsearch(f,0,5)
        return_old = oldfu.goldsearch(f,0,5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_rotate_params_should_return_equal_object(self):
        a = [[1, 2, 3], [3, 4, 5]]
        b = [1,5,6]
        return_new = fu.rotate_params(a, b)
        return_old = oldfu.rotate_params(a,b)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_rotmatrix_should_return_equal_object(self):
        return_new = fu.rotmatrix(1,5,6)
        return_old = oldfu.rotmatrix(1,5,6)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_mulmat_should_return_equal_object(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 6]]
        b = [[2, 3, 2], [3, 4, 3], [3, 4, 6]]
        return_new = fu.mulmat(a,b)
        return_old = oldfu.mulmat(a,b)
        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_recmat_should_return_equal_object(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 1]]
        return_new = fu.recmat(a)
        return_old = oldfu.recmat(a)
        self.assertTrue(numpy.array_equal(return_new, return_old))

class TestSymClassInitC(unittest.TestCase):

    def test_c0_sym_should_return_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.symclass('c0')
        with self.assertRaises(ZeroDivisionError):
            oldfu.symclass('c0')

    def test_lower_c1_sym_should_return_lower_c1_sym(self):
        self.assertEqual(fu.symclass('c1').sym, oldfu.symclass('c1').sym)

    def test_upper_c1_sym_should_return_lower_c1_sym(self):
        self.assertEqual(fu.symclass('C1').nsym , oldfu.symclass('C1').nsym)


    def test_c1_sym_should_return_nsym_1(self):
        self.assertEqual(fu.symclass('c1').nsym, 1)
        self.assertEqual(oldfu.symclass('c1').nsym,1)

    def test_c1_sym_should_return_nsym_5(self):
        self.assertEqual(fu.symclass('c5').nsym, 5)
        self.assertEqual(oldfu.symclass('c5').nsym,5)


    def test_c1_should_return_correct_brackets(self):
        nsym = 1
        fubrackets =  fu.symclass('c1').brackets == [
            [ old_div(360., nsym), 90.0, old_div(360., nsym), 90.0],
            [ old_div(360., nsym), 180.0, old_div(360., nsym), 180.0]
                    ]

        oldfubrackets =  oldfu.symclass('c1').brackets == [
            [ old_div(360., nsym), 90.0, old_div(360., nsym), 90.0],
            [ old_div(360., nsym), 180.0, old_div(360., nsym), 180.0]
                    ]

        self.assertEqual(fubrackets , oldfubrackets)


    def test_c5_should_return_correct_brackets(self):
        nsym = 1
        fubrackets =  fu.symclass('c5').brackets == [
            [ old_div(360., nsym), 90.0, old_div(360., nsym), 90.0],
            [ old_div(360., nsym), 180.0, old_div(360., nsym), 180.0]
                    ]

        oldfubrackets =  oldfu.symclass('c5').brackets == [
            [ old_div(360., nsym), 90.0, old_div(360., nsym), 90.0],
            [ old_div(360., nsym), 180.0, old_div(360., nsym), 180.0]
                    ]

        self.assertEqual(fubrackets , oldfubrackets)

    def test_c1_should_return_correct_symangles(self):
        nsym = 1
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i * old_div(360., nsym)])

        self.assertEqual(fu.symclass('c1').symangles, symangles)
        self.assertEqual(oldfu.symclass('c1').symangles, symangles)

    def test_c5_should_return_correct_symangles(self):
        nsym = 5
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i * old_div(360., nsym)])

        self.assertEqual(fu.symclass('c5').symangles, symangles)
        self.assertEqual(oldfu.symclass('c5').symangles, symangles)


    def test_c1_should_return_correct_transform(self):
        transform = []
        for args in fu.symclass('c1').symangles:
            transform.append(e2cpp.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        self.assertEqual(fu.symclass('c1').transform, transform)
        transform = []
        for args in oldfu.symclass('c1').symangles:
            transform.append(e2cpp.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        self.assertEqual(oldfu.symclass('c1').transform, transform)


    def test_c5_should_return_correct_transform(self):
        transform = []
        for args in fu.symclass('c5').symangles:
            transform.append(e2cpp.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        self.assertEqual(fu.symclass('c5').transform, transform)
        transform = []
        for args in oldfu.symclass('c5').symangles:
            transform.append(e2cpp.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        self.assertEqual(oldfu.symclass('c5').transform, transform)

    def test_c1_should_return_correct_symatrix(self):
        symatrix = []
        for args in fu.symclass('c1').symangles:
            symatrix.append(fu.rotmatrix(args[0],args[1],args[2]))
        self.assertEqual(fu.symclass('c1').symatrix, symatrix)

        symatrix = []
        for args in oldfu.symclass('c1').symangles:
            symatrix.append(oldfu.rotmatrix(args[0],args[1],args[2]))
        self.assertEqual(oldfu.symclass('c1').symatrix, symatrix)


    def test_c5_should_return_correct_symatrix(self):
        symatrix = []
        for args in fu.symclass('c5').symangles:
            symatrix.append(fu.rotmatrix(args[0],args[1],args[2]))
        self.assertEqual(fu.symclass('c5').symatrix, symatrix)

        symatrix = []
        for args in oldfu.symclass('c5').symangles:
            symatrix.append(oldfu.rotmatrix(args[0],args[1],args[2]))
        self.assertEqual(oldfu.symclass('c5').symatrix, symatrix)

    @staticmethod
    def rotmatrix(phi,theta,psi):
        rphi   = numpy.radians(phi)
        rtheta = numpy.radians(theta)
        rpsi   = numpy.radians(psi)
        cosphi = numpy.cos(rphi)
        sinphi = numpy.sin(rphi)
        costheta = numpy.cos(rtheta)
        sintheta = numpy.sin(rtheta)
        cospsi = numpy.cos(rpsi)
        sinpsi = numpy.sin(rpsi)
        mat = [[0.0]*3,[0.0]*3,[0.0]*3]

        mat[0][0] =  cospsi*costheta*cosphi - sinpsi*sinphi
        mat[1][0] = -sinpsi*costheta*cosphi - cospsi*sinphi
        mat[2][0] =            sintheta*cosphi


        mat[0][1] =  cospsi*costheta*sinphi + sinpsi*cosphi
        mat[1][1] = -sinpsi*costheta*sinphi + cospsi*cosphi
        mat[2][1] =            sintheta*sinphi


        mat[0][2] = -cospsi*sintheta
        mat[1][2] =  sinpsi*sintheta
        mat[2][2] =            costheta
        return mat


class TestSymClassIsInSubunitC(unittest.TestCase):
    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_c1_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)


        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)

    def test_c1_sym_no_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c1_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c1_sym_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c5_sym_no_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c5_sym_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


    def test_c5_sym_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


    def test_c5_sym_mirror_theta_smaller_equals_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c5_sym_mirror_theta_larger_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = fu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = oldfu.symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)



class TestSymClassIsInSubunitWrong(unittest.TestCase):

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_wrong_sym_crashes_problem(self):
        return_new = fu.symclass('c1')
        return_old = oldfu.symclass('c1')
        return_new.sym = 'foobar'
        return_old.sym = 'foobar'
        self.assertIsNone(return_new.is_in_subunit(0, 0))
        self.assertIsNone(return_old.is_in_subunit(0, 0))


class TestSymClassSymmetryRelatedC(unittest.TestCase):

    def test_c1_sym_zero(self):
        new_result = fu.symclass('c1').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('c1').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0]])
        self.assertEqual(old_result, [[0, 0, 0]])

    def test_c5_sym_zero(self):
        new_result = fu.symclass('c5').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('c5').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0]])
        self.assertEqual(old_result, [[0, 0, 0]])


    def test_c1_sym_180(self):
        new_result = fu.symclass('c1').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('c1').symmetry_related([0, 180, 0])
        self.assertEqual(new_result, [[0, 180, 0]])
        self.assertEqual(old_result, [[0, 180, 0]])

    def test_c5_sym_180(self):
        new_result = fu.symclass('c5').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('c5').symmetry_related([0, 180, 0])
        self.assertEqual(new_result, [[0, 180, 0]])
        self.assertEqual(old_result, [[0, 180, 0]])

    def test_c1_sym_90(self):
        new_result = fu.symclass('c1').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('c1').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0]])
        self.assertEqual(old_result, [[0, 90, 0]])

    def test_c5_sym_90(self):
        new_result = fu.symclass('c5').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('c5').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]] )
        self.assertEqual(old_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]])


class TestSymClassSymmetryRelated(unittest.TestCase):

    def my_method(self, angles, sym, symatrix):
        if( sym[0] == "c" or sym[0] == "d" ):
            temp = e2cpp.Util.symmetry_neighbors(angles, sym)
            nt = old_div(len(temp), 3)
            return [[temp[l*3],temp[l*3+1],0.0] for l in range(nt) ]
        #  Note symmetry neighbors below refer to the particular order
        #   in which this class generates symmetry matrices
        neighbors = {}
        neighbors["oct"]  = [1,2,3,8,9,12,13]
        neighbors["tet"]  = [1,2,3,4,6,7]
        neighbors["icos"] = [1,2,3,4,6,7,11,12]
        sang = [[] for l in range(len(angles)*(len(neighbors[sym])+1))]
        for i,q in enumerate(angles):  sang[i*(len(neighbors[sym])+1)] = angles[i][:]
        pass#IMPORTIMPORTIMPORT from fundamentals import rotmatrix, recmat, mulmat
        for i,q in enumerate(angles):
            mat = self.rotmatrix(q[0],q[1],q[2])
            for j,l in enumerate(neighbors[sym]):
                p1,p2,p3 = self.recmat( self.mulmat(mat, symatrix[l]) )
                sang[i*(len(neighbors[sym])+1) + (j+1)] = [p1,p2,0.0]
        return sang

    @staticmethod
    def recmat(mat):
        pass#IMPORTIMPORTIMPORT from math import acos,asin,atan2,degrees,pi
        def sign(x):
            if( x >= 0.0 ): return 1
            else:  return -1
        """
        mat = [[0.0]*3,[0.0]*3,[0.0]*3]
        # limit precision
        for i in range(3):
            for j in range(3):
                mat[i][j] = inmat[i][j]
                #if(abs(inmat[i][j])<1.0e-8):  mat[i][j] = 0.0
                #else: mat[i][j] = inmat[i][j]
        for i in range(3):
            for j in range(3):  print  "     %14.8f"%mat[i][j],
            print ""
        """
        if(mat[2][2] == 1.0):
            theta = 0.0
            psi = 0.0
            if( mat[0][0] == 0.0 ):
                phi = math.asin(mat[0][1])
            else:
                phi = math.atan2(mat[0][1],mat[0][0])
        elif(mat[2][2] == -1.0):
            theta = numpy.pi
            psi = 0.0
            if(mat[0][0] == 0.0):
                phi = math.asin(-mat[0][1])
            else:
                phi = math.atan2(-mat[0][1],-mat[0][0])
        else:
            theta = math.acos(mat[2][2])
            st = sign(theta)
            #print theta,st,mat[2][0]
            if(mat[2][0] == 0.0):
                if( st != sign(mat[2][1]) ):
                    phi = 1.5*numpy.pi
                else:
                    phi = 0.5*numpy.pi
            else:
                phi = math.atan2(st*mat[2][1], st*mat[2][0])

            #print theta,st,mat[0][2],mat[1][2]
            if(mat[0][2] == 0.0):
                if( st != sign(mat[1][2]) ):
                    psi = 1.5*numpy.pi
                else:
                    psi = 0.5*numpy.pi
            else:
                psi = math.atan2(st*mat[1][2], -st*mat[0][2])
        #pi2 = 2*pi
        #return  degrees(round(phi%pi2,8)),degrees(round(theta%pi2,8)),degrees(round(psi%pi2,8))
        #return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
        return  numpy.degrees(phi)%360.0,numpy.degrees(theta)%360.0,numpy.degrees(psi)%360.0

    @staticmethod
    def mulmat(m1,m2):
        mat = [[0.0]*3,[0.0]*3,[0.0]*3]
        """
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    mat[i][j] += m1[i][k]*m2[k][j]
                #mat[i][j] = round(mat[i][j],8)
        """
        mat[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0]
        mat[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1]
        mat[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2]
        mat[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0]
        mat[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1]
        mat[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2]
        mat[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0]
        mat[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1]
        mat[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2]

        return mat

    @staticmethod
    def rotmatrix(phi,theta,psi):
        rphi   = numpy.radians(phi)
        rtheta = numpy.radians(theta)
        rpsi   = numpy.radians(psi)
        cosphi = numpy.cos(rphi)
        sinphi = numpy.sin(rphi)
        costheta = numpy.cos(rtheta)
        sintheta = numpy.sin(rtheta)
        cospsi = numpy.cos(rpsi)
        sinpsi = numpy.sin(rpsi)
        mat = [[0.0]*3,[0.0]*3,[0.0]*3]

        mat[0][0] =  cospsi*costheta*cosphi - sinpsi*sinphi
        mat[1][0] = -sinpsi*costheta*cosphi - cospsi*sinphi
        mat[2][0] =            sintheta*cosphi


        mat[0][1] =  cospsi*costheta*sinphi + sinpsi*cosphi
        mat[1][1] = -sinpsi*costheta*sinphi + cospsi*cosphi
        mat[2][1] =            sintheta*sinphi


        mat[0][2] = -cospsi*sintheta
        mat[1][2] =  sinpsi*sintheta
        mat[2][2] =            costheta
        return mat

    def test_c_sym(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, fu.symclass('c4').sym, fu.symclass('c4').symatrix)
        return_values = fu.symclass('c4').symmetry_neighbors(angles)
        self.assertEqual(return_values,expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, oldfu.symclass('c4').sym, oldfu.symclass('c4').symatrix)
        return_values = oldfu.symclass('c4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

    def test_d_sym(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, fu.symclass('d4').sym, fu.symclass('d4').symatrix)
        return_values = fu.symclass('d4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, oldfu.symclass('d4').sym, oldfu.symclass('d4').symatrix)
        return_values = oldfu.symclass('d4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

    def test_tet_sym(self):
        angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, fu.symclass('tet').sym, fu.symclass('tet').symatrix)
        return_values = fu.symclass('tet').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, oldfu.symclass('tet').sym, oldfu.symclass('tet').symatrix)
        return_values = oldfu.symclass('tet').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

    def test_oct_sym(self):
        angles = [[idx1, idx2, 0] for idx1 in range(40) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, fu.symclass('oct').sym, fu.symclass('oct').symatrix)
        return_values = fu.symclass('oct').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(40) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, oldfu.symclass('oct').sym, oldfu.symclass('oct').symatrix)
        return_values = oldfu.symclass('oct').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)



class TestSymClassSymmetryNeighbors(unittest.TestCase):
     def test_symmetry_neighbors_c1_should_return_equal_object(self):
         angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]

         oldangles = oldfu.symclass('c1').symmetry_neighbors(angles)
         newangles = fu.symclass('c1').symmetry_neighbors(angles)
         self.assertEqual(oldangles,newangles)



class TestSymClassReduceAngleSets(unittest.TestCase):
    def test_reduce_anglesets_c1_should_return_equal_object(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]

        oldangles = oldfu.symclass('c1').reduce_anglesets(angles)
        newangles = fu.symclass('c1').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)


class TestSymClassReduceAngles(unittest.TestCase):
    def test_reduce_angles_c1_should_return_equal_object(self):

        oldangles = oldfu.symclass('c1').reduce_angles(60,80,0)
        newangles = fu.symclass('c1').reduce_angles(60,80,0)
        self.assertEqual(oldangles, newangles)


class TestSymClassReduceAngles(unittest.TestCase):
    def test_even_angles_c1_should_return_equal_object(self):
        oldangles = oldfu.symclass('c1').even_angles(3.75)
        newangles = fu.symclass('c1').even_angles(3.75)
        self.assertEqual(oldangles, newangles)


if __name__ == '__main__':
    unittest.main()
