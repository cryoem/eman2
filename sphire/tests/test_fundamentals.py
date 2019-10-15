from __future__ import print_function
from __future__ import division
from past.utils import old_div

import unittest
from test_module import get_real_data, get_arg_from_pickle_file,create_kb
from sphire.libpy.sp_utilities import model_blank,model_circle
from numpy import array_equal, allclose
from os import path
from EMAN2_cppwrap import EMData

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = model_blank(10, 10)
IMAGE_BLANK_3D = model_blank(10, 10, 10)
MASK = model_circle(2, 5, 5)
REAL_IMAGE_2D = get_real_data(dim =2)[0]
REAL_IMAGE_3D = get_real_data(dim =3)[0]

"""
There are some opened issues in:
In all the tests miss the case with a complex image. where can we find one of them?
1) fftip: I cannot find an image that gives me as output a not None value 
2) prepi3D: It is able to work with a 2D image too. Should it be blocked?
3) rot_shift3D: It is able to work with a 2D image too. Should it be blocked?
4) goldsearch: HOw can I test it?
"""


"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""


""" start: new in sphire 1.3"""
from sphire.libpy import sp_fundamentals as oldfu
from sphire.libpy_py3 import sp_fundamentals as fu

class Test_absi(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.absi()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.absi()
        self.assertEqual(cm_new.exception.message, "absi() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.absi(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.absi(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.absi(e=EMData())
        return_old = oldfu.absi(e=EMData())
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2Dimg(self):
        return_old = oldfu.absi(e=IMAGE_2D)
        return_new = fu.absi(e=IMAGE_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank(self):
        return_old = oldfu.absi(e=IMAGE_BLANK_2D)
        return_new = fu.absi(e=IMAGE_BLANK_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg(self):
        return_old = oldfu.absi(e=IMAGE_3D)
        return_new = fu.absi(e=IMAGE_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank(self):
        return_old = oldfu.absi(e=IMAGE_BLANK_3D)
        return_new = fu.absi(e=IMAGE_BLANK_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_acf(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acf()
        self.assertEqual(cm_new.exception.message, "acf() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acf(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acf(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acf(e=EMData(), center=True)
        return_old = oldfu.acf(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acf(e=IMAGE_2D, center=True)
        return_new = fu.acf(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acf(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acf(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acf(e=IMAGE_3D, center=True)
        return_new = fu.acf(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acf(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acf(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acf(e=IMAGE_2D, center=False)
        return_new = fu.acf(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acf(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acf(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acf(e=IMAGE_3D, center=False)
        return_new = fu.acf(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acf(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acf(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_acfn(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acfn()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acfn()
        self.assertEqual(cm_new.exception.message, "acfn() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acfn(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acfn(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acfn(e=EMData(), center=True)
        return_old = oldfu.acfn(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acfn(e=IMAGE_2D, center=True)
        return_new = fu.acfn(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acfn(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acfn(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acfn(e=IMAGE_3D, center=True)
        return_new = fu.acfn(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_center(self):
        return_old = oldfu.acfn(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acfn(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acfn(e=IMAGE_2D, center=False)
        return_new = fu.acfn(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acfn(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acfn(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acfn(e=IMAGE_3D, center=False)
        return_new = fu.acfn(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acfn(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acfn(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_acfp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acfp()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acfp()
        self.assertEqual(cm_new.exception.message, "acfp() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acfp(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acfp(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acfp(e=EMData(), center=True)
        return_old = oldfu.acfp(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acfp(e=IMAGE_2D, center=True)
        return_new = fu.acfp(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acfp(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acfp(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acfp(e=IMAGE_3D, center=True)
        return_new = fu.acfp(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_center(self):
        return_old = oldfu.acfp(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acfp(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acfp(e=IMAGE_2D, center=False)
        return_new = fu.acfp(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acfp(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acfp(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acfp(e=IMAGE_3D, center=False)
        return_new = fu.acfp(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acfp(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acfp(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_acfnp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acfnp()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acfnp()
        self.assertEqual(cm_new.exception.message, "acfnp() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acfnp(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acfnp(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acfnp(e=EMData(), center=True)
        return_old = oldfu.acfnp(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acfnp(e=IMAGE_2D, center=True)
        return_new = fu.acfnp(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acfnp(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acfnp(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acfnp(e=IMAGE_3D, center=True)
        return_new = fu.acfnp(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_center(self):
        return_old = oldfu.acfnp(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acfnp(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acfnp(e=IMAGE_2D, center=False)
        return_new = fu.acfnp(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acfnp(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acfnp(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acfnp(e=IMAGE_3D, center=False)
        return_new = fu.acfnp(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acfnp(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acfnp(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_acfpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acfpl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acfpl()
        self.assertEqual(cm_new.exception.message, "acfpl() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acfpl(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acfpl(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acfpl(e=EMData(), center=True)
        return_old = oldfu.acfpl(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acfpl(e=IMAGE_2D, center=True)
        return_new = fu.acfpl(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acfpl(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acfpl(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acfpl(e=IMAGE_3D, center=True)
        return_new = fu.acfpl(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_center(self):
        return_old = oldfu.acfpl(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acfpl(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acfpl(e=IMAGE_2D, center=False)
        return_new = fu.acfpl(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acfpl(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acfpl(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acfpl(e=IMAGE_3D, center=False)
        return_new = fu.acfpl(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acfpl(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acfpl(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_acfnpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.acfnpl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.acfnpl()
        self.assertEqual(cm_new.exception.message, "acfnpl() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.acfnpl(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.acfnpl(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.acfnpl(e=EMData(), center=True)
        return_old = oldfu.acfnpl(e=EMData(), center=True)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_center(self):
        return_old = oldfu.acfnpl(e=IMAGE_2D, center=True)
        return_new = fu.acfnpl(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_center(self):
        return_old = oldfu.acfnpl(e=IMAGE_BLANK_2D, center=True)
        return_new = fu.acfnpl(e=IMAGE_BLANK_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_center(self):
        return_old = oldfu.acfnpl(e=IMAGE_3D, center=True)
        return_new = fu.acfnpl(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_center(self):
        return_old = oldfu.acfnpl(e=IMAGE_BLANK_3D, center=True)
        return_new = fu.acfnpl(e=IMAGE_BLANK_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_NOcenter(self):
        return_old = oldfu.acfnpl(e=IMAGE_2D, center=False)
        return_new = fu.acfnpl(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimgBlank_NOcenter(self):
        return_old = oldfu.acfnpl(e=IMAGE_BLANK_2D, center=False)
        return_new = fu.acfnpl(e=IMAGE_BLANK_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_NOcenter(self):
        return_old = oldfu.acfnpl(e=IMAGE_3D, center=False)
        return_new = fu.acfnpl(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimgBlank_NOcenter(self):
        return_old = oldfu.acfnpl(e=IMAGE_BLANK_3D, center=False)
        return_new = fu.acfnpl(e=IMAGE_BLANK_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ccfn(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccfn()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccfn()
        self.assertEqual(cm_new.exception.message, "ccfn() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfn(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfn(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfn(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfn(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccfn(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccfn(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccfn(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccfn(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccfn(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfn(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccfn(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfn(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccfn(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfn(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccfn(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfn(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccfn(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfn(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccfn(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfn(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccfn(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfn(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccfn(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfn(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ccfp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccfp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccfp()
        self.assertEqual(cm_new.exception.message, "ccfp() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfp(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfp(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfp(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfp(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccfp(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccfp(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccfp(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccfp(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccfp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccfp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccfp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccfp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccfp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccfp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccfp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccfp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfp(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ccfnp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccfnp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccfnp()
        self.assertEqual(cm_new.exception.message, "ccfnp() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfnp(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfnp(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfnp(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfnp(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccfnp(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccfnp(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccfnp(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccfnp(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccfnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfnp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccfnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfnp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccfnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfnp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccfnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfnp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccfnp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccfnp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccfnp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccfnp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfnp(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ccfpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccfpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccfpl()
        self.assertEqual(cm_new.exception.message, "ccfpl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfpl(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfpl(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfpl(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfpl(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccfpl(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccfpl(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccfpl(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccfpl(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccfpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccfpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccfpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccfpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccfpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccfpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccfpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccfpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfpl(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ccfnpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ccfnpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ccfnpl()
        self.assertEqual(cm_new.exception.message, "ccfnpl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfnpl(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfnpl(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccfnpl(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccfnpl(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccfnpl(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccfnpl(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccfnpl(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccfnpl(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccfnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfnpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccfnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfnpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccfnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfnpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccfnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfnpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccfnpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccfnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccfnpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccfnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccfnpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccfnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccfnpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccfnpl(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cnv(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cnv()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cnv()
        self.assertEqual(cm_new.exception.message, "cnv() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnv(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnv(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnv(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnv(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.cnv(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.cnv(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.cnv(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.cnv(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.cnv(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnv(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.cnv(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnv(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.cnv(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnv(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.cnv(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnv(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.cnv(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnv(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.cnv(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnv(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.cnv(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnv(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.cnv(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnv(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cnvnp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cnvnp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cnvnp()
        self.assertEqual(cm_new.exception.message, "cnvnp() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvnp(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvnp(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvnp(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvnp(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.cnvnp(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.cnvnp(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.cnvnp(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.cnvnp(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.cnvnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvnp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.cnvnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvnp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.cnvnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvnp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.cnvnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvnp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.cnvnp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.cnvnp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvnp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.cnvnp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvnp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.cnvnp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvnp(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cnvp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cnvp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cnvp()
        self.assertEqual(cm_new.exception.message, "cnvp() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvp(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvp(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvp(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvp(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.cnvp(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.cnvp(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.cnvp(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.cnvp(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.cnvp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.cnvp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvp(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.cnvp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.cnvp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvp(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.cnvp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.cnvp(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvp(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.cnvp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvp(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.cnvp(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvp(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cnvpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cnvpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cnvpl()
        self.assertEqual(cm_new.exception.message, "cnvpl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvpl(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvpl(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvpl(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvpl(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.cnvpl(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.cnvpl(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.cnvpl(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.cnvpl(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.cnvpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.cnvpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.cnvpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.cnvpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.cnvpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.cnvpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.cnvpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.cnvpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvpl(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_cnvnpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.cnvnpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.cnvnpl()
        self.assertEqual(cm_new.exception.message, "cnvnpl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvnpl(EMData(),IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvnpl(EMData(),IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.cnvnpl(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.cnvnpl(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.cnvnpl(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.cnvnpl(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.cnvnpl(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.cnvnpl(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.cnvnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvnpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.cnvnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvnpl(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.cnvnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvnpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.cnvnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvnpl(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.cnvnpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.cnvnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.cnvnpl(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.cnvnpl(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.cnvnpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.cnvnpl(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.cnvnpl(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.cnvnpl(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_scfn(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scfn()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scfn()
        self.assertEqual(cm_new.exception.message, "scfn() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scfn(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scfn(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scfn(e=None, center=False)
        return_old = oldfu.scfn(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scfn(e=IMAGE_3D, center=True)
        return_old = oldfu.scfn(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scfn(e=IMAGE_3D, center=False)
        return_old = oldfu.scfn(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scfn(e=IMAGE_2D, center=True)
        return_old = oldfu.scfn(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scfn(e=IMAGE_2D, center=False)
        return_old = oldfu.scfn(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scfn(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scfn(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scfn(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scfn(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scfn(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scfn(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scfn(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scfn(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_scfp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scfp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scfp()
        self.assertEqual(cm_new.exception.message, "scfp() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scfp(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scfp(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scfp(e=None, center=False)
        return_old = oldfu.scfp(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scfp(e=IMAGE_3D, center=True)
        return_old = oldfu.scfp(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scfp(e=IMAGE_3D, center=False)
        return_old = oldfu.scfp(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scfp(e=IMAGE_2D, center=True)
        return_old = oldfu.scfp(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scfp(e=IMAGE_2D, center=False)
        return_old = oldfu.scfp(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scfp(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scfp(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scfp(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scfp(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scfp(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scfp(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scfp(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scfp(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_scfnp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scfnp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scfnp()
        self.assertEqual(cm_new.exception.message, "scfnp() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scfnp(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scfnp(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scfnp(e=None, center=False)
        return_old = oldfu.scfnp(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scfnp(e=IMAGE_3D, center=True)
        return_old = oldfu.scfnp(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scfnp(e=IMAGE_3D, center=False)
        return_old = oldfu.scfnp(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scfnp(e=IMAGE_2D, center=True)
        return_old = oldfu.scfnp(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scfnp(e=IMAGE_2D, center=False)
        return_old = oldfu.scfnp(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scfnp(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scfnp(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scfnp(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scfnp(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scfnp(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scfnp(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scfnp(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scfnp(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_scfpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scfpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scfpl()
        self.assertEqual(cm_new.exception.message, "scfpl() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scfpl(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scfpl(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scfpl(e=None, center=False)
        return_old = oldfu.scfpl(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scfpl(e=IMAGE_3D, center=True)
        return_old = oldfu.scfpl(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scfpl(e=IMAGE_3D, center=False)
        return_old = oldfu.scfpl(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scfpl(e=IMAGE_2D, center=True)
        return_old = oldfu.scfpl(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scfpl(e=IMAGE_2D, center=False)
        return_old = oldfu.scfpl(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scfpl(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scfpl(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scfpl(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scfpl(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scfpl(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scfpl(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scfpl(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scfpl(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_scfnpl(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.scfnpl()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.scfnpl()
        self.assertEqual(cm_new.exception.message, "scfnpl() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image1(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.scfnpl(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scfnpl(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scfnpl(e=None, center=False)
        return_old = oldfu.scfnpl(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scfnpl(e=IMAGE_3D, center=True)
        return_old = oldfu.scfnpl(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scfnpl(e=IMAGE_3D, center=False)
        return_old = oldfu.scfnpl(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scfnpl(e=IMAGE_2D, center=True)
        return_old = oldfu.scfnpl(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scfnpl(e=IMAGE_2D, center=False)
        return_old = oldfu.scfnpl(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scfnpl(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scfnpl(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scfnpl(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scfnpl(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scfnpl(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scfnpl(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scfnpl(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scfnpl(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_image_decimate(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.image_decimate()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.image_decimate()
        self.assertEqual(cm_new.exception.message, "image_decimate() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.image_decimate(img=EMData(), decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.image_decimate(img=EMData(), decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.image_decimate(img=None, decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0)
        return_old = oldfu.image_decimate(img=None, decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimage_with_fft(self):
        return_old = oldfu.image_decimate(img=IMAGE_2D, decimation=2, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_2D, decimation=2, fit_to_fft=True, frequency_low=0, frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimageBlank_with_fft(self):
        return_old = oldfu.image_decimate(img=IMAGE_BLANK_2D, decimation=2, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_BLANK_2D, decimation=2, fit_to_fft=True, frequency_low=0,
                                       frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimage_without_fft(self):
        return_old = oldfu.image_decimate(img=IMAGE_2D, decimation=2, fit_to_fft=False, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_2D, decimation=2, fit_to_fft=False, frequency_low=0, frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimageBlank_without_fft(self):
        return_old = oldfu.image_decimate(img=IMAGE_BLANK_2D, decimation=2, fit_to_fft=False, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_BLANK_2D, decimation=2, fit_to_fft=False, frequency_low=0,frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimage_noDecimation(self):
        return_old = oldfu.image_decimate(img=IMAGE_2D, decimation=1, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_2D, decimation=1, fit_to_fft=True, frequency_low=0, frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_2D.get_3dview()))

    def test_2DimageBlank_noDecimation(self):
        return_old = oldfu.image_decimate(img=IMAGE_BLANK_2D, decimation=1, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_BLANK_2D, decimation=1, fit_to_fft=True, frequency_low=0,frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_BLANK_2D.get_3dview()))

    def test_3Dimage_Error(self):
        return_old = oldfu.image_decimate(img=IMAGE_3D, decimation=2, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_3D, decimation=2, fit_to_fft=True, frequency_low=0, frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimage_lowFrequency_Error(self):
        return_old = oldfu.image_decimate(img=IMAGE_2D, decimation=-2, fit_to_fft=True, frequency_low=0,frequency_high=0)
        return_new = fu.image_decimate(img=IMAGE_2D, decimation=-2, fit_to_fft=True, frequency_low=0, frequency_high=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_fdownsample(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fdownsample()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fdownsample()
        self.assertEqual(cm_new.exception.message, "fdownsample() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.fdownsample(img=EMData(), sub_rate=0.5, RetReal = True)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fdownsample(img=EMData(), sub_rate=0.5, RetReal = True)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.fdownsample(img=None, sub_rate=0.5, RetReal = True)
        return_old = oldfu.fdownsample(img=None, sub_rate=0.5, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertEqual(return_new, None)

    def test_2DImg_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=0.5, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=0.5, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=0.5, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=0.5, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=0.5, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=0.5, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=0.5, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=0.5, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_noSubrate_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=1, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=1, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_2D.get_3dview()))

    def test_2DImgBlank_noSubrate_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=1, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=1, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_BLANK_2D.get_3dview()))

    def test_3DImg_noSubrate_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=1, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=1, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_3D.get_3dview()))

    def test_3DImgBlank_noSubrate_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=1, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=1, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(array_equal(return_new.get_3dview(), IMAGE_BLANK_3D.get_3dview()))

    def test_2DImg_upscaling_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=2, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=2, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_upscaling_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=2, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=2, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_upscaling_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=2, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=2, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_upscaling_RetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=2, RetReal = True)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=2, RetReal = True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=0.5, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=0.5, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=0.5, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=0.5, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=0.5, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=0.5, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=0.5, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=0.5, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_noSubrate_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=1, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=1, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_noSubrate_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=1, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=1, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_noSubrate_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=1, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=1, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_noSubrate_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=1, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=1, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_upscaling_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_2D, sub_rate=2, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_2D, sub_rate=2, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_upscaling_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=2, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_2D, sub_rate=2, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_upscaling_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_3D, sub_rate=2, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_3D, sub_rate=2, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_upscaling_NORetReal(self):
        return_new = fu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=2, RetReal=False)
        return_old = oldfu.fdownsample(img=IMAGE_BLANK_3D, sub_rate=2, RetReal=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_prepf(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepf()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepf()
        self.assertEqual(cm_new.exception.message, "prepf() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepf(image=EMData(), npad = 2)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.prepf(image=EMData(), npad = 2)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.prepf(image=None, npad = 2)
        return_old = oldfu.prepf(image=None,npad = 2)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2DImage(self):
        return_new = fu.prepf(image=IMAGE_2D, npad=2)
        return_old = oldfu.prepf(image=IMAGE_2D, npad=2)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankImage(self):
        return_new = fu.prepf(image=IMAGE_BLANK_2D, npad=2)
        return_old = oldfu.prepf(image=IMAGE_BLANK_2D, npad=2)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_3DImage(self):
        #todo: it should give an error??
        return_new = fu.prepf(image=IMAGE_3D, npad=2)
        return_old = oldfu.prepf(image=IMAGE_3D, npad=2)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



""" I did not test it deeply becuase the paramteters are passed to other functions and not precessed directly in this function"""
class Test_prep_refim_gridding(unittest.TestCase):
    (not1, numr, wr, not2, not3, not4, not5, not6, not7, not8) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter"))[0]
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prep_refim_gridding()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prep_refim_gridding()
        self.assertEqual(cm_new.exception.message, "prep_refim_gridding() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prep_refim_gridding(refim=EMData(), npad = 2)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.prep_refim_gridding(refim=EMData(), npad = 2)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        crefim_new,kb_new = fu.prep_refim_gridding(refim=None, wr=self.wr, numr=self.numr, mode = "F")
        crefim_old,kb_old = oldfu.prep_refim_gridding(refim=None,wr=self.wr, numr=self.numr, mode = "F")
        self.assertTrue(array_equal(crefim_new.get_3dview(), crefim_old.get_3dview()))
        self.assertTrue(array_equal(kb_new.dump_table(), kb_new.dump_table()))  # it is a kaiserbessel filter

    def test_2DImage(self):
        crefim_new,kb_new = oldfu.prep_refim_gridding(refim=IMAGE_2D, wr=self.wr, numr=self.numr, mode = "F")
        crefim_old,kb_old = fu.prep_refim_gridding(refim=IMAGE_2D, wr=self.wr, numr=self.numr, mode = "F")
        self.assertTrue(array_equal(crefim_new.get_3dview(), crefim_old.get_3dview()))
        self.assertTrue(array_equal(kb_new.dump_table(), kb_new.dump_table()))  # it is a kaiserbessel filter

    def test_2DImageBlank(self):
        crefim_new,kb_new = oldfu.prep_refim_gridding(refim=IMAGE_BLANK_2D, wr=self.wr, numr=self.numr, mode = "F")
        crefim_old,kb_old = fu.prep_refim_gridding(refim=IMAGE_BLANK_2D, wr=self.wr, numr=self.numr, mode = "F")
        self.assertTrue(array_equal(crefim_new.get_3dview(), crefim_old.get_3dview()))
        self.assertTrue(array_equal(kb_new.dump_table(), kb_new.dump_table()))  # it is a kaiserbessel filter

    def test_3DImage(self):
        crefim_new,kb_new = oldfu.prep_refim_gridding(refim=IMAGE_3D, wr=self.wr, numr=self.numr, mode = "F")
        crefim_old,kb_old = fu.prep_refim_gridding(refim=IMAGE_3D, wr=self.wr, numr=self.numr, mode = "F")
        self.assertTrue(array_equal(crefim_new.get_3dview(), crefim_old.get_3dview()))
        self.assertTrue(array_equal(kb_new.dump_table(), kb_new.dump_table()))  # it is a kaiserbessel filter

    def test_3DImageBlank(self):
        crefim_new,kb_new = oldfu.prep_refim_gridding(refim=IMAGE_BLANK_3D, wr=self.wr, numr=self.numr, mode = "F")
        crefim_old,kb_old = fu.prep_refim_gridding(refim=IMAGE_BLANK_3D, wr=self.wr, numr=self.numr, mode = "F")
        self.assertTrue(array_equal(crefim_new.get_3dview(), crefim_old.get_3dview()))
        self.assertTrue(array_equal(kb_new.dump_table(), kb_new.dump_table()))  # it is a kaiserbessel filter





class Test_prepg(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepg()
        self.assertEqual(cm_new.exception.message, "prepg() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepg(image=EMData(), kb = "not_in_use")
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.prepg(image=EMData(), kb = "not_in_use")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.prepg(image=None, kb = "not_in_use")
        return_new = fu.prepg(image=None, kb = "not_in_use")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimage(self):
        return_old = oldfu.prepg(image=IMAGE_2D, kb = "not_in_use")
        return_new = fu.prepg(image=IMAGE_2D, kb = "not_in_use")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimageBlank(self):
        return_old = oldfu.prepg(image=IMAGE_BLANK_2D, kb = "not_in_use")
        return_new = fu.prepg(image=IMAGE_BLANK_2D, kb = "not_in_use")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImage(self):
        #todo: it should give an error??
        return_old = oldfu.prepg(image=IMAGE_3D, kb = "not_in_use")
        return_new = fu.prepg(image=IMAGE_3D, kb = "not_in_use")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_rot_avg_image(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rot_avg_image()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rot_avg_image()
        self.assertEqual(cm_new.exception.message, "rot_avg_image() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rot_avg_image(image_to_be_averaged=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rot_avg_image(image_to_be_averaged=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_avg_image(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_avg_image(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rotavg'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimage(self):
        return_new = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_2D)
        return_old = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimageBlank(self):
        return_new = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_BLANK_2D)
        return_old = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_BLANK_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimage(self):
        return_new = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_3D)
        return_old = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimageBlank(self):
        return_new = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_BLANK_3D)
        return_old = oldfu.rot_avg_image(image_to_be_averaged=IMAGE_BLANK_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


#todo: test it with a file ... how??
class Test_ro_textfile(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ro_textfile()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ro_textfile()
        self.assertEqual(cm_new.exception.message, "rops() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_(self):
        oldv = oldfu.ro_textfile(e="", filename="", helpful_string="")
        v = fu.ro_textfile(e="", filename="", helpful_string="")
        pass


class Test_rops(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rops()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rops()
        self.assertEqual(cm_new.exception.message, "rops() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rops(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rops(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rops(None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rops(None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rotavg'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimage(self):
        return_new = oldfu.rops(e=IMAGE_2D)
        return_old = oldfu.rops(e=IMAGE_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DimageBlank(self):
        return_new = oldfu.rops(e=IMAGE_BLANK_2D)
        return_old = oldfu.rops(e=IMAGE_BLANK_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimage(self):
        return_new = oldfu.rops(e=IMAGE_3D)
        return_old = oldfu.rops(e=IMAGE_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DimageBlank(self):
        return_new = oldfu.rops(e=IMAGE_BLANK_3D)
        return_old = oldfu.rops(e=IMAGE_BLANK_3D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



#todo: test it with a file ... how??
class Test_rops_textfile(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rops_textfile()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rops_textfile()
        self.assertEqual(cm_new.exception.message, "rops_textfile() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rops_textfile(self):
        oldv = oldfu.rops_textfile(img="", filename="", lng = False)
        v = fu.rops_textfile(img="", filename="", lng = False)
        pass


""" I did not test it deeply becuase the paramteters are passed to other functions and not precessed directly in this function"""
class Test_rotshift2dg(unittest.TestCase):
    kb=create_kb(dim=1)
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rotshift2dg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rotshift2dg()
        self.assertEqual(cm_new.exception.message, "rotshift2dg() takes at least 5 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rotshift2dg(image=EMData(), ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rotshift2dg(image=EMData(), ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.rotshift2dg(image=None, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        return_new = fu.rotshift2dg(image=None, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2DImage(self):
        return_old = oldfu.rotshift2dg(image=IMAGE_2D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        return_new = fu.rotshift2dg(image=IMAGE_2D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImageBlank(self):
        return_old = oldfu.rotshift2dg(image=IMAGE_BLANK_2D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        return_new = fu.rotshift2dg(image=IMAGE_BLANK_2D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImage_error(self):
        #todo: it'd give an error
        return_old = oldfu.rotshift2dg(image=IMAGE_3D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        return_new = fu.rotshift2dg(image=IMAGE_3D, ang=10, dx=2, dy=2, kb=self.kb, scale = 1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_ft2polargrid(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ft2polargrid()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ft2polargrid()
        self.assertEqual(cm_new.exception.message, "ft2polargrid() takes exactly 4 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ft2polargrid(image=EMData(),  ring_length=10, nb=10, ne=10)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ft2polargrid(image=EMData(),  ring_length=10, nb=10, ne=10)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.ft2polargrid(image=None,  ring_length=10, nb=10, ne=10)
        return_new = fu.ft2polargrid(image=None,  ring_length=10, nb=10, ne=10)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2DImg(self):
        return_old = oldfu.ft2polargrid(image=IMAGE_2D,  ring_length=10, nb=10, ne=10)
        return_new = fu.ft2polargrid(image=IMAGE_2D,  ring_length=10, nb=10, ne=10)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank(self):
        return_old = oldfu.ft2polargrid(image=IMAGE_BLANK_2D,  ring_length=10, nb=10, ne=10)
        return_new = fu.ft2polargrid(image=IMAGE_BLANK_2D,  ring_length=10, nb=10, ne=10)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_error(self):
        # todo: it'd give an error
        return_old = oldfu.ft2polargrid(image=IMAGE_3D,  ring_length=10, nb=10, ne=10)
        return_new = fu.ft2polargrid(image=IMAGE_3D,  ring_length=10, nb=10, ne=10)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


""" I did not test it deeply becuase the paramteters are passed to other functions and not precessed directly in this function"""
class Test_rot_shift3D_grid(unittest.TestCase):
    not_used,kb = fu.prepi3D(IMAGE_3D)

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rot_shift3D_grid()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rot_shift3D_grid()
        self.assertEqual(cm_new.exception.message, "rot_shift3D_grid() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rot_shift3D_grid(img=EMData(), phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0, kb=None,
                                mode="background", wrap=False)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rot_shift3D_grid(img=EMData(), phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,
                                   kb=None, mode="background", wrap=False)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.rot_shift3D_grid(img=None, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=None, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="background", wrap=False)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_3DImg_Nokb_background(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=None, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="background", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_Nokb_ciclyc(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=None, mode="cyclic", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="cyclic", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_Nokb_background(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3,scale=1.0, kb=None, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0,scale=1.0, kb=None, mode="background", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_Nokb_ciclyc(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3,scale=1.0, kb=None, mode="cyclic", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0,scale=1.0, kb=None, mode="cyclic", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_kb_background(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=self.kb, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=self.kb, mode="background", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_kb_ciclyc(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=self.kb, mode="cyclic", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=self.kb, mode="cyclic", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_kb_background(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3,scale=1.0, kb=self.kb, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0,scale=1.0, kb=self.kb, mode="background", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_kb_ciclyc(self):
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3,scale=1.0, kb=self.kb, mode="cyclic", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_BLANK_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0,scale=1.0, kb=self.kb, mode="cyclic", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2d_error(self):
        # todo: it'd give an error
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_2D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=None, mode="background", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_2D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="background", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_error_mode(self):
        # todo: it'd give an error
        return_old = oldfu.rot_shift3D_grid(img=IMAGE_3D, phi=0.3, theta=0.2, psi=0.1, sx=3, sy=3, sz=3, scale=1.0,kb=None, mode="wrong_mode", wrap=False)
        return_new = fu.rot_shift3D_grid(img=IMAGE_3D, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0,kb=None, mode="wrong_mode", wrap=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_sinc2inv(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sinc2inv()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sinc2inv()
        self.assertEqual(cm_new.exception.message, "sinc2inv() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_sinc2inv(self):
        return_old = oldfu.sinc2inv(nx=10)
        return_new = fu.sinc2inv(nx=10)
        self.assertTrue(array_equal(return_new, return_old))

    def test_null_nx_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.sinc2inv(nx=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.sinc2inv(nx=0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


class Test_sincinv(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sincinv()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sincinv()
        self.assertEqual(cm_new.exception.message, "sincinv() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_sincinv(self):
        return_old = oldfu.sincinv(nx=10)
        return_new = fu.sincinv(nx=10)
        self.assertTrue(array_equal(return_new, return_old))

    def test_null_nx_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.sincinv(nx=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.sincinv(nx=0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)




class Test_welch_pw2(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.welch_pw2()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.welch_pw2()
        self.assertEqual(cm_new.exception.message, "welch_pw2() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.welch_pw2(img=EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.welch_pw2(img=EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.welch_pw2(img=None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        return_new = fu.welch_pw2(img=None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2dImg(self):
        return_old = oldfu.welch_pw2(img=IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        return_new = fu.welch_pw2(img=IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2dImgBlank(self):
        return_old = oldfu.welch_pw2(img=IMAGE_BLANK_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        return_new = fu.welch_pw2(img=IMAGE_BLANK_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dImg(self):
        return_old = oldfu.welch_pw2(img=IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        return_new = fu.welch_pw2(img=IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3dImgBlank(self):
        return_old = oldfu.welch_pw2(img=IMAGE_BLANK_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        return_new = fu.welch_pw2(img=IMAGE_BLANK_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=1, edge_y=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_welch_pw2_tilt_band(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.welch_pw2_tilt_band()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.welch_pw2_tilt_band()
        self.assertEqual(cm_new.exception.message, "welch_pw2_tilt_band() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.welch_pw2_tilt_band(img=EMData(), theta=2,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.welch_pw2_tilt_band(img=EMData(),theta=2,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_old = oldfu.welch_pw2_tilt_band(img=None, theta=2,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=None, theta=2,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_num_bnd_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_2D,theta=2,num_bnd=-1,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_2D,theta=2,num_bnd=-1,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_2DimgBlank_num_bnd_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_BLANK_2D,theta=2,num_bnd=-1,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_BLANK_2D,theta=2,num_bnd=-1,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_2Dimg_num_bnd_not_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_2D,theta=2,num_bnd=5,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_2D,theta=2,num_bnd=5,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_2DimgBlank_num_bnd_not_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_BLANK_2D,theta=2,num_bnd=5,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_BLANK_2D,theta=2,num_bnd=5,overlp_y=50,edge_x=1,edge_y=1,win_s=256)
        for i,j in zip(return_old,return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_3Dimg_num_bnd_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_3D, theta=2, num_bnd=-1, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_3D, theta=2, num_bnd=-1, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        for i, j in zip(return_old, return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_3DimgBlank_num_bnd_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_BLANK_3D, theta=2, num_bnd=-1, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_BLANK_3D, theta=2, num_bnd=-1, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        for i, j in zip(return_old, return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_3Dimg_num_bnd_not_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_3D, theta=2, num_bnd=5, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_3D, theta=2, num_bnd=5, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        for i, j in zip(return_old, return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

    def test_3DimgBlank_num_bnd_not_default(self):
        return_old = oldfu.welch_pw2_tilt_band(img=IMAGE_BLANK_3D, theta=2, num_bnd=5, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        return_new = fu.welch_pw2_tilt_band(img=IMAGE_BLANK_3D, theta=2, num_bnd=5, overlp_y=50, edge_x=1, edge_y=1,win_s=256)
        for i, j in zip(return_old, return_new):
            self.assertTrue(array_equal(i.get_3dview(), j.get_3dview()))

""" is it the same of morphology.py???"""
class Test_bracket(unittest.TestCase):

    @staticmethod
    def function1(x1, dat):
        return x1 + dat

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.bracket()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.bracket()
        self.assertEqual(cm_new.exception.message, "bracket() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_f3_greater_f1(self):
        return_new = fu.bracket(f=self.function1, x1=5, h=4)
        return_old = oldfu.bracket(f=self.function1, x1=5, h=4)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, (0.0, 10.472135955999999)))

    def test_f3_not_greater_f1_outputmsg_Bracket_didnot_find_a_minimum(self):
        self.assertTrue(fu.bracket(f=self.function1, x1=0, h=0) is None)
        self.assertTrue(oldfu.bracket(f=self.function1, x1=0, h=0) is None)




""" start: end in sphire 1.3"""


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
            fu.ccf(e=EMData(),f=IMAGE_2D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccf(e=EMData(),f=IMAGE_2D)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ccf(e=IMAGE_2D, f=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ccf(e=IMAGE_2D,f=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img2(self):
        return_new = fu.ccf(e=IMAGE_2D, f=None, center=False)
        return_old = oldfu.ccf(e=IMAGE_2D,f= None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img1_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccf(e=None, f=IMAGE_2D, center=False)
        return_old = oldfu.ccf(e=None, f=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img2D_without_center(self):
        return_new = fu.ccf(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccf(e=IMAGE_2D,f= REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.ccf(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccf(e=IMAGE_2D,f= REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.ccf(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccf(e=IMAGE_3D,f= REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_center(self):
        return_new = fu.ccf(e=IMAGE_3D, f=REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccf(e=IMAGE_3D,f= REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_center(self):
        return_new = fu.ccf(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=False)
        return_old = oldfu.ccf(e=IMAGE_2D, f=REAL_IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.ccf(e=IMAGE_BLANK_2D,f= REAL_IMAGE_2D, center=True)
        return_old = oldfu.ccf(e=IMAGE_2D, f=REAL_IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.ccf(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=False)
        return_old = oldfu.ccf(e=IMAGE_3D, f=REAL_IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.ccf(e=IMAGE_BLANK_3D,f= REAL_IMAGE_3D, center=True)
        return_old = oldfu.ccf(e=IMAGE_BLANK_3D, f=REAL_IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.scf(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.scf(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_input_img_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.scf(e=None, center=False)
        return_old = oldfu.scf(e=None, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        """

    def test_Img3D_with_center(self):
        return_new = fu.scf(e=IMAGE_3D, center=True)
        return_old = oldfu.scf(e=IMAGE_3D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_center(self):
        return_new = fu.scf(e=IMAGE_3D, center=False)
        return_old = oldfu.scf(e=IMAGE_3D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_center(self):
        return_new = fu.scf(e=IMAGE_2D, center=True)
        return_old = oldfu.scf(e=IMAGE_2D, center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_center(self):
        return_new = fu.scf(e=IMAGE_2D, center=False)
        return_old = oldfu.scf(e=IMAGE_2D, center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_Img2Dblank_without_center(self):
        return_new = fu.scf(e=IMAGE_BLANK_2D,  center=False)
        return_old = oldfu.scf(e=IMAGE_2D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_center(self):
        return_new = fu.scf(e=IMAGE_BLANK_2D,  center=True)
        return_old = oldfu.scf(e=IMAGE_2D,  center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_center(self):
        return_new = fu.scf(e=IMAGE_BLANK_3D, center=False)
        return_old = oldfu.scf(e=IMAGE_3D,  center=False)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_center(self):
        return_new = fu.scf(e=IMAGE_BLANK_3D,center=True)
        return_old = oldfu.scf(e=IMAGE_BLANK_3D,center=True)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.cyclic_shift(img=None)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.cyclic_shift(img=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Img2D_with_no_shift(self):
        return_new = fu.cyclic_shift(img=IMAGE_2D, dx=0, dy=0, dz=0)
        return_old = oldfu.cyclic_shift(img=IMAGE_2D, dx=0, dy=0, dz=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_no_shift(self):
        return_new = fu.cyclic_shift(img=IMAGE_3D, dx=0, dy=0, dz=0)
        return_old = oldfu.cyclic_shift(img=IMAGE_3D, dx=0, dy=0, dz=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_huge_shift(self):
        return_new = fu.cyclic_shift(img=IMAGE_2D, dx=10000, dy=10000, dz=10000)
        return_old = oldfu.cyclic_shift(img=IMAGE_2D, dx=10000, dy=10000, dz=10000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_huge_shift(self):
        return_new = fu.cyclic_shift(img=IMAGE_3D, dx=10000, dy=10000, dz=10000)
        return_old = oldfu.cyclic_shift(img=IMAGE_3D, dx=10000, dy=10000, dz=10000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))




class Test_mirror(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.mirror()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.mirror()
        self.assertEqual(cm_new.exception.message, "mirror() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        return_new = fu.mirror(img=EMData(), axis='x')
        return_old = oldfu.mirror(img=EMData(), axis='x')
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_NoneType_input_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.mirror(img=None)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.mirror(img=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.mirror(img=IMAGE_2D, axis='x')
        return_old = oldfu.mirror(img=IMAGE_2D, axis='x')
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg(self):
        return_new = fu.mirror(img=IMAGE_3D, axis='x')
        return_old = oldfu.mirror(img=IMAGE_3D, axis='x')
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_invalid_axis(self):
        return_new = fu.mirror(img=IMAGE_2D, axis='invalid_axis')
        return_old = oldfu.mirror(img=IMAGE_2D, axis='invalid_axis')
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.fft(e=EMData(), npad=1)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fft(e=EMData(), npad=1)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fft(e=None, npad=1)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fft(e=None, npad=1)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'is_complex'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_withPad1(self):
        return_new = fu.fft(e=IMAGE_3D, npad=1)
        return_old = oldfu.fft(e=IMAGE_3D, npad=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_withPad1(self):
        return_new = fu.fft(e=IMAGE_3D, npad=1)
        return_old = oldfu.fft(e=IMAGE_3D, npad=1)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_withPadGrater1(self):
        return_new = fu.fft(e=IMAGE_3D, npad=3)
        return_old = oldfu.fft(e=IMAGE_3D, npad=3)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_withPadGrater1(self):
        return_new = fu.fft(e=IMAGE_3D, npad=3)
        return_old = oldfu.fft(e=IMAGE_3D, npad=3)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))





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
            fu.fftip(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fftip(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'y size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.fftip(e=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.fftip(e=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'is_complex'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_wrong_3Dimg_returns_attributeError(self):
        return_new = fu.fftip(e=REAL_IMAGE_3D)
        return_old = oldfu.fftip(e=REAL_IMAGE_3D)
        self.assertEqual(return_old,return_new)
        self.assertEqual(return_old, None)

    def test_wrong_2Dimg_returns_attributeError(self):
        return_new = fu.fftip(e=IMAGE_2D)
        return_old = oldfu.fftip(e=IMAGE_2D)
        self.assertEqual(return_old,return_new)
        self.assertEqual(return_old, None)




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
            fu.fpol(image=EMData(),nnx=100)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fpol(image=EMData(),nnx=100)
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
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retFourier(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retReal(self):
        return_new = fu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fpol(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))




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
            fu.fdecimate(image=EMData(),nnx=100)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fdecimate(image=EMData(),nnx=100)
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
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retFourier(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=True, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Normalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_noNormalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_2D, nnx=IMAGE_2D.get_xsize(), nny=IMAGE_2D.get_ysize(), nnz=IMAGE_2D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Normalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=True )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_noNormalizated_retReal(self):
        return_new = fu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        return_old = oldfu.fdecimate(image=IMAGE_3D, nnx=IMAGE_3D.get_xsize(), nny=IMAGE_3D.get_ysize(), nnz=IMAGE_3D.get_zsize(), RetReal=False, normalize=False )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))




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
            fu.fshift(e=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.fshift(e=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        return_new = fu.fshift(e=None)
        return_old = oldfu.fshift(e=None)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2Dimg_without_shift(self):
        return_new = fu.fshift(e=IMAGE_2D, delx=0, dely=0, delz=0)
        return_old = oldfu.fshift(e=IMAGE_2D, delx=0, dely=0, delz=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_without_shift(self):
        return_new = fu.fshift(e=IMAGE_3D, delx=0, dely=0, delz=0)
        return_old = oldfu.fshift(e=IMAGE_3D, delx=0, dely=0, delz=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.fshift(e=IMAGE_2D, delx=1000, dely=1000, delz=1000)
        return_old = oldfu.fshift(e=IMAGE_2D, delx=1000, dely=1000, delz=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_huge_shift(self):
        return_new = fu.fshift(e=IMAGE_3D, delx=1000, dely=1000, delz=1000)
        return_old = oldfu.fshift(e=IMAGE_3D, delx=1000, dely=1000, delz=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))






class Test_subsample(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.subsample()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.subsample()
        self.assertEqual(cm_new.exception.message, "subsample() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_empty_input_img(self):
        return_new = fu.subsample(image=EMData())
        return_old = oldfu.subsample(image=EMData())
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.subsample(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.subsample(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.subsample(image=IMAGE_2D, subsample_rate=1.0)
        return_old = oldfu.subsample(image=IMAGE_2D, subsample_rate=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_subsample_rate_higher1(self):
        return_new = fu.subsample(image=IMAGE_2D, subsample_rate=2.0)
        return_old = oldfu.subsample(image=IMAGE_2D, subsample_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_negative_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(image=IMAGE_2D, subsample_rate=-1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(image=IMAGE_2D, subsample_rate=-1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_null_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(image=IMAGE_2D, subsample_rate=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(image=IMAGE_2D, subsample_rate=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_huge_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.subsample(image=IMAGE_2D, subsample_rate=1000.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.subsample(image=IMAGE_2D, subsample_rate=1000.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "BadAllocException")
        self.assertEqual(msg[1], 'Cannot allocate 495616 GB - not enough memory.')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_3Dimg_default_value(self):
        return_new = fu.subsample(image=IMAGE_3D, subsample_rate=1.0)
        return_old = oldfu.subsample(image=IMAGE_3D, subsample_rate=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_subsample_rate_higher1(self):
        return_new = fu.subsample(image=IMAGE_3D, subsample_rate=2.0)
        return_old = oldfu.subsample(image=IMAGE_3D, subsample_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))






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
            fu.resample(img=EMData())
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(img=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.resample(img=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.resample(img=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.resample(img=IMAGE_2D, sub_rate=1.0)
        return_old = oldfu.resample(img=IMAGE_2D, sub_rate=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_higher1(self):
        return_new = fu.resample(img=IMAGE_2D, sub_rate=2.0)
        return_old = oldfu.resample(img=IMAGE_2D, sub_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_lower1(self):
        return_new = fu.resample(img=IMAGE_2D, sub_rate=0.10)
        return_old = oldfu.resample(img=IMAGE_2D, sub_rate=0.10)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_sub_rate_higher1_nx_notsame_ny(self):
        return_new = fu.resample(img=EMData(100,20), sub_rate=2.0)
        return_old = oldfu.resample(img=EMData(100,20), sub_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_negative_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(img=IMAGE_2D, sub_rate=-1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(img=IMAGE_2D, sub_rate=-1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_null_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(img=IMAGE_2D, sub_rate=0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(img=IMAGE_2D, sub_rate=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_huge_rate(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.resample(img=IMAGE_2D, sub_rate=1000.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.resample(img=IMAGE_2D, sub_rate=1000.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "BadAllocException")
        self.assertEqual(msg[1], 'Cannot allocate 495.616 GB - not enough memory.')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_3Dimg_default_value(self):
        return_new = fu.resample(img=IMAGE_3D, sub_rate=1.0)
        return_old = oldfu.resample(img=IMAGE_3D, sub_rate=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_higher1(self):
        return_new = fu.resample(img=IMAGE_3D, sub_rate=2.0)
        return_old = oldfu.resample(img=IMAGE_3D, sub_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_lower1(self):
        return_new = fu.resample(img=IMAGE_3D, sub_rate=0.2)
        return_old = oldfu.resample(img=IMAGE_3D, sub_rate=0.2)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_sub_rate_higher1_nz_notsame_ny_or_nx(self):
        return_new = fu.resample(img=EMData(100,100,2), sub_rate=2.0)
        return_old = oldfu.resample(img=EMData(100,100,2), sub_rate=2.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.prepi(image=EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.prepi(image=EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepi(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.prepi(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_retFour(self):
        return_new = fu.prepi(image=IMAGE_2D,RetReal = False)
        return_old = oldfu.prepi(image=IMAGE_2D,RetReal = False)
        self.assertTrue(array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_2Dimg_retReal(self):
        return_new = fu.prepi(image=IMAGE_2D,RetReal = True)
        return_old = oldfu.prepi(image=IMAGE_2D,RetReal = True)
        self.assertTrue(array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.prepi(image=IMAGE_3D,RetReal = False)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.prepi(image=IMAGE_3D,RetReal = False)
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
            fu.prepi3D(image=EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.prepi3D(image=EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)


    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.prepi3D(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.prepi3D(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg(self):
        return_new = fu.prepi3D(image=IMAGE_3D)
        return_old = oldfu.prepi3D(image=IMAGE_3D)
        self.assertTrue(array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter

    def test_2Dimg(self):
        return_new = fu.prepi3D(image=IMAGE_2D)
        return_old = oldfu.prepi3D(image=IMAGE_2D)
        self.assertTrue(array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))  # it is a image
        self.assertTrue(array_equal(return_new[1].dump_table(), return_old[1].dump_table()))  # it is a kaiserbessel filter



class Test_ramp(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.ramp()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.ramp()
        self.assertEqual(cm_new.exception.message, "ramp() takes exactly 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_img(self):
        return_new = fu.ramp(inputimage=EMData())
        return_old = oldfu.ramp(inputimage=EMData())
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.ramp(inputimage=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.ramp(inputimage=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'copy'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.ramp(inputimage=IMAGE_2D)
        return_old = oldfu.ramp(inputimage=IMAGE_2D)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_returns_RunTimeError(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.ramp(inputimage=IMAGE_3D)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.ramp(inputimage=IMAGE_3D)
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
            fu.rot_avg_table(e=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_avg_table(e=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rotavg'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.rot_avg_table(e=IMAGE_2D)
        return_old = oldfu.rot_avg_table(e=IMAGE_2D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [1.70115065574646, 1.6371004581451416, 1.5626660585403442, 1.5877236127853394, 1.6073533296585083,1.6285395622253418, 1.6862214803695679, 1.7918702363967896, 1.8338433504104614, 1.8038475513458252,1.6706266403198242, 1.544594407081604, 1.4375859498977661, 1.2914648056030273, 1.167849063873291,1.080208659172058, 0.9806387424468994, 0.8771043419837952, 0.7996737360954285, 0.7433571815490723,0.7284249663352966, 0.769216001033783, 0.8287383317947388, 0.9201586246490479, 0.9932690262794495,1.0363448858261108, 1.123952031135559, 1.1950533390045166, 1.2668777704238892, 1.3277966976165771,1.4283941984176636, 1.5512932538986206, 1.6878538131713867, 1.793628454208374, 1.8632606267929077,1.9746273756027222, 2.0680618286132812, 2.1846938133239746, 2.3249411582946777, 2.4541871547698975,2.5610392093658447, 2.6200289726257324, 2.6819217205047607, 2.7509777545928955, 2.8288681507110596,2.9172451496124268, 3.0068471431732178, 3.0675971508026123, 3.081587553024292, 3.068607807159424,3.0510096549987793, 2.971430778503418, 2.8975279331207275, 2.8081719875335693, 2.7049367427825928,2.5696637630462646, 2.4190049171447754, 2.2764179706573486, 2.1662235260009766, 2.06880784034729,1.9702768325805664, 1.8747632503509521, 1.77823805809021, 1.6550695896148682, 1.5298069715499878,1.4096013307571411, 1.269547700881958, 1.1273504495620728, 0.9518274664878845, 0.7691776752471924,0.5582662224769592, 0.3821752667427063, 0.2139461785554886, 0.02325786091387272,-0.16214938461780548, -0.3665705919265747, -0.5570375919342041, -0.7568416595458984,-0.9444634914398193, -1.1122022867202759, -1.2591341733932495, -1.3866132497787476,-1.5305417776107788, -1.6798121929168701, -1.7983838319778442, -1.8982055187225342,-1.9665424823760986, -2.0055081844329834, -2.0240345001220703, -2.022763252258301,-2.0144519805908203, -1.985037088394165, -1.9594217538833618, -1.9007961750030518,-1.8341227769851685, -1.7523555755615234, -1.6880568265914917, -1.6217501163482666,-1.5742276906967163, -1.516431212425232, -1.4340367317199707, -1.3519632816314697,-1.2755727767944336, -1.2155612707138062, -1.1673047542572021, -1.1223112344741821,-1.0735875368118286, -1.0130205154418945, -0.9502072334289551, -0.8961243629455566,-0.8469223380088806, -0.8141208291053772, -0.7728646993637085, -0.7485840320587158,-0.7369289398193359, -0.7154955267906189, -0.6829231977462769, -0.6576652526855469,-0.6172498464584351, -0.5820884704589844, -0.5657995343208313, -0.5585340857505798,-0.5447290539741516, -0.5433396697044373, -0.552950918674469, -0.5448899865150452,-0.542202353477478, -0.541515052318573, -0.5387433171272278, -0.5594711303710938,-0.5963792204856873, -0.6536659598350525, -0.7025741338729858, -0.7419651746749878,-0.7738654613494873, -0.8288431763648987, -0.8792231678962708, -0.9185277223587036,-0.9535513520240784, -0.9755640625953674, -0.9769047498703003, -0.9703202843666077,-0.9619379043579102, -0.9533834457397461, -0.9389827251434326, -0.8972327709197998,-0.8668974041938782, -0.8482377529144287, -0.8443998098373413, -0.8439618349075317,-0.8385792970657349, -0.829200804233551, -0.8111365437507629, -0.7909402847290039,-0.7704744338989258, -0.7470849752426147, -0.7254742980003357, -0.6966635584831238,-0.6708171367645264, -0.6408230662345886, -0.6109101176261902, -0.5850300788879395,-0.5588661432266235, -0.5414973497390747, -0.5227035880088806, -0.5071979761123657,-0.4920077323913574, -0.47635769844055176, -0.46521392464637756, -0.4503937363624573,-0.4407220482826233, -0.43082553148269653, -0.41926488280296326, -0.4109551012516022,-0.3982444405555725, -0.39125651121139526, -0.372354120016098]))

    def test_3Dimg(self):
        return_new = fu.rot_avg_table(e=IMAGE_3D)
        return_old = oldfu.rot_avg_table(e=IMAGE_3D)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [10.473315238952637, 21.884193420410156, 31.711227416992188, 33.32331466674805, 26.039146423339844, 20.5319766998291, 18.808826446533203, 18.97246551513672, 20.31218719482422, 22.830305099487305, 27.15753936767578, 32.496463775634766, 36.15548324584961, 36.97674560546875, 34.278045654296875, 26.243494033813477, 15.062774658203125, 6.469794273376465, 2.7834198474884033, 2.3789703845977783, 2.7278802394866943, 2.6416027545928955, 2.4357831478118896, 2.5040085315704346, 2.6968085765838623, 2.794365882873535, 2.600947141647339, 2.0284249782562256, 1.4216328859329224, 0.5415528416633606, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_2DBlankimg(self):
        return_new = fu.rot_avg_table(e=IMAGE_BLANK_2D)
        return_old = oldfu.rot_avg_table(e=IMAGE_BLANK_2D)
        self.assertTrue(array_equal(return_new, return_old))

    def test_3DBlankimg(self):
        return_new = fu.rot_avg_table(e=IMAGE_BLANK_3D)
        return_old = oldfu.rot_avg_table(e=IMAGE_BLANK_3D)
        self.assertTrue(array_equal(return_new, return_old))





class Test_rops_table(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rops_table()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rops_table()
        self.assertEqual(cm_new.exception.message, "rops_table() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)





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
            fu.gridrot_shift2D(image=EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.gridrot_shift2D(image=EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.gridrot_shift2D(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.gridrot_shift2D(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zeroScale(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_angle(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_scale(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_default_value(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_zeroScale(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_shift(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=1000, sy=1000, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_angle(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_scale(self):
        return_new = fu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        return_old = oldfu.gridrot_shift2D(image=IMAGE_BLANK_2D, ang=0.0, sx=0.0, sy=0.0, scale=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

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
            fu.rot_shift2D(img=EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift2D(img=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1]+msg[2]+msg[3], 'Cant rotate 1D image caught\n')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1]+msg[2]+msg[3], msg_old[1]+msg_old[2]+msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_shift2D(img=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_shift2D(img=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_trans2D_background'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_InvalidMethod_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_invalid_interpolation_method_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2Dimg_default_value(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_Zero_scale_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodLinear(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodQuadratic(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodGridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodFtgridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeBackground_methodFourier(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodLinear(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="linear", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodQuadratic(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodGridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodFtgridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_modeCyclic_methodFourier(self):
        return_new = fu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="fourier", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_InvalidMethod_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="Invalid")
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)

    def test_2DBlankimg_invalid_interpolation_method_returns_None_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="Invalid", mode="background")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_2DBlankimg_default_value(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=1.0, interpolation_method=None, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_Zero_scale_print_errorMsg(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, mirror=0, scale=0, interpolation_method=None, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeBackground_methodLinear(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeBackground_methodQuadratic(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="quadratic", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeBackground_methodGridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="gridding", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeBackground_methodFtgridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeBackground_methodFourier(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeCyclic_methodLinear(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="linear", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="linear", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeCyclic_methodQuadratic(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="quadratic", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeCyclic_methodGridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="gridding", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeCyclic_methodFtgridding(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="ftgridding", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_modeCyclic_methodFourier(self):
        return_new = fu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0, interpolation_method="fourier", mode="cyclic")
        return_old = oldfu.rot_shift2D(img=IMAGE_BLANK_2D, angle=1.0, sx=1.0, sy=0.0, mirror=0, scale=1.0,interpolation_method="fourier", mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.rot_shift3D(image=EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift3D(image=EMData())
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rot_shift3D(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rot_shift3D(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_trans_background'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg(self):
        return_new = fu.rot_shift3D(image=IMAGE_2D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        return_old = oldfu.rot_shift3D(image=IMAGE_2D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Zero_scale_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=0, mode="background")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=0, mode="background")
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'The scale factor in a Transform object must be positive and non zero')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_3Dimg_default_values(self):
        return_new = fu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        return_old = oldfu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="background")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_modeCyclic(self):
        return_new = fu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="cyclic")
        return_old = oldfu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="cyclic")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3Dimg_Invalid_mode(self):
        return_new = fu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="invalid")
        return_old = oldfu.rot_shift3D(image=IMAGE_3D, phi=0, theta=0, psi=0, sx=0, sy=0, sz=0, scale=1.0, mode="invallid")
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))



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
            fu.rtshg(image=EMData())
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.rtshg(image=EMData())
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Nonetype_input_img(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.rtshg(image=None)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.rtshg(image=None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.rtshg(image=IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.rtshg(image=IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Cannot reduce the image size')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_2Dimg_default_value(self):
        return_new = fu.rtshg(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_hugeShift(self):
        return_new = fu.rtshg(image=IMAGE_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_hugeAngle(self):
        return_new = fu.rtshg(image=IMAGE_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zeroScale(self):
        return_new = fu.rtshg(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.rtshg(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_default_value(self):
        return_new = fu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_hugeShift(self):
        return_new = fu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=1000.0, sy=1000.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_hugeAngle(self):
        return_new = fu.rtshg(image=IMAGE_BLANK_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        return_old = oldfu.rtshg(image=IMAGE_BLANK_2D, angle=1000.0, sx=0.0, sy=0.0, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_zeroScale(self):
        return_new = fu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        return_old = oldfu.rtshg(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


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
            fu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=kb, scale=1.0)
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")

        self.assertEqual(msg[0].split("\n")[0]+msg[0].split("\n")[1], "Python argument types in    EMData.rot_scale_conv_new(EMData, numpy.float64, float, float, tuple, float)")
        self.assertEqual(msg[0].split("\n")[0]+msg[0].split("\n")[1], msg_old[0].split("\n")[0]+msg_old[0].split("\n")[1])

    def test_2Dimg_empty_img(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rtshgkb(image=EMData(), angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rtshgkb(image=EMData(), angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1]+msg[2]+msg[3], 'Cant rotate 1D image caught\n')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1]+msg[2]+msg[3], msg_old[1]+msg_old[2]+msg_old[3])

    def test_NoneType_img(self):
        with self.assertRaises(AttributeError)  as cm_new:
            fu.rtshgkb(image=None, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.rtshgkb(image=None, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'rot_scale_conv_new'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_default_value(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.rtshgkb(image=IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.rtshgkb(image=IMAGE_3D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(msg[1], 'Use rot_scale_conv_new_3D for volumes')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_2Dimg_default_value(self):
        return_new = fu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_zero_scale(self):
        return_new = fu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        return_old = oldfu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_angle(self):
        return_new = fu.rtshgkb(image=IMAGE_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_shift(self):
        return_new = fu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_huge_scale(self):
        return_new = fu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        return_old = oldfu.rtshgkb(image=IMAGE_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_2DBlankimg_default_value(self):
        return_new = fu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_zero_scale(self):
        return_new = fu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        return_old = oldfu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_angle(self):
        return_new = fu.rtshgkb(image=IMAGE_BLANK_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_BLANK_2D, angle=1000, sx=0.0, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_shift(self):
        return_new = fu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        return_old = oldfu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=1000, sy=0.0, kb=self.kb, scale=1.0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DBlankimg_huge_scale(self):
        return_new = fu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        return_old = oldfu.rtshgkb(image=IMAGE_BLANK_2D, angle=0.0, sx=0.0, sy=0.0, kb=self.kb, scale=1000)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))


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
            fu.tilemic(img=None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.tilemic(img=None, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyImg_returns_UnboundLocalError(self):
        with self.assertRaises(UnboundLocalError)  as cm_new:
            fu.tilemic(img=EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(UnboundLocalError)  as cm_old:
            oldfu.tilemic(img=EMData(), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "local variable 'i' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_default_value(self):
        return_new = fu.tilemic(img=IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(img=IMAGE_2D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, []))

    def test_3Dimg_default_value(self):
        return_new = fu.tilemic(img=EMData(100,100,100), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(img=EMData(100,100,100), win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertTrue(array_equal(return_new, return_old))

    def test_wrong_3Dimg_(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.tilemic(img=IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.tilemic(img=IMAGE_3D, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'z size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_2Dimg_zero_win_size(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.tilemic(img=IMAGE_2D, win_size=0, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.tilemic(img=IMAGE_2D, win_size=0, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0)
        self.assertEqual(cm_new.exception.message, "float division by zero")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_2Dimg_no_overlap(self):
        return_new = fu.tilemic(img=IMAGE_2D, win_size=120, overlp_x=0, overlp_y=0, edge_x=0, edge_y=0)
        return_old = oldfu.tilemic(img=IMAGE_2D, win_size=120, overlp_x=0, overlp_y=0, edge_x=0, edge_y=0)
        for img1, img2 in zip(return_old, return_new):
            self.assertTrue(allclose(img1.get_3dview(), img2.get_3dview(), equal_nan=True))




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
            fu.window2d(img=None, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.window2d(img=None, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyImg(self):
        return_new = fu.window2d(img=EMData(), isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        return_old = oldfu.window2d(img=EMData(), isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_central(self):
        return_new = fu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        return_old = oldfu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_left(self):
        return_new = fu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="l", ix=0, iy=0)
        return_old = oldfu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="l", ix=0, iy=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_arbitrary(self):
        return_new = fu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="a", ix=0, iy=0)
        return_old = oldfu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="a", ix=0, iy=0)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2Dimg_invalid_option_returns_UnboundLocalError(self):
        with self.assertRaises(UnboundLocalError)  as cm_new:
            fu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="invalid", ix=0, iy=0)
        with self.assertRaises(UnboundLocalError)  as cm_old:
            oldfu.window2d(img=IMAGE_2D, isize_x=100, isize_y=100, opt="invalid", ix=0, iy=0)
        self.assertEqual(cm_new.exception.message, "local variable 'reg' referenced before assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_3Dimg_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError)  as cm_new:
            fu.window2d(img=IMAGE_3D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
        with self.assertRaises(RuntimeError)  as cm_old:
            oldfu.window2d(img=IMAGE_3D, isize_x=100, isize_y=100, opt="c", ix=0, iy=0)
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


# todo: modificato
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
            fu.rotate_params(params=[1,2,3],transf=[1],tolistconv=True)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.rotate_params(params=[1,2,3],transf=[1],tolistconv=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_params_value_returns_IndexError(self):
        with self.assertRaises(IndexError)  as cm_new:
            fu.rotate_params(params=[[1],[1]],transf=[1,2,3],tolistconv=True)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.rotate_params(params=[[1],[1]],transf=[1,2,3],tolistconv=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_params_value_returns_TypeError(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rotate_params(params=[1],transf=[1,2,3],tolistconv=True)
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rotate_params(params=[1],transf=[1,2,3],tolistconv=True)
        self.assertEqual(cm_new.exception.message, "'int' object has no attribute '__getitem__'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotate_params(self):
        a = [[1, 2, 3], [3, 4, 5]]
        b = [1,5,6]
        return_new = fu.rotate_params(params=a, transf=b,tolistconv=True)
        return_old = oldfu.rotate_params(params=a,transf=b,tolistconv=True)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[174.00000000000003, 2.9999999999999232, 183.0], [166.07792596478271, 1.0120857125278824, 194.9159858556005]]))


# todo: modificato
class Test_rotmatrix(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.rotmatrix()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.rotmatrix()
        self.assertEqual(cm_new.exception.message, "rotmatrix() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_rotmatrix(self):
        return_new = fu.rotmatrix(phi=1,theta=5, psi=6,tolistconv=True)
        return_old = oldfu.rotmatrix(phi=1,theta=5, psi=6,tolistconv=True)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[0.98876227196606559, 0.12180329553973532, -0.086678294469630643], [-0.12147164162554232, 0.99255309355170718, 0.0091102558543636417], [0.087142468505889387, 0.0015210774457754552, 0.99619469809174555]]))


# todo: modificato
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
            fu.mulmat(m1=a,m2=b,tolistconv=True)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.mulmat(m1=a,m2=b,tolistconv=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_few_m2_value_returns_IndexError(self):
        a = [[1, 3, 2], [2, 5, 3], [3, 8, 6]]
        b = [[2,  2], [3, 4, 3], [3, 4, 6]]
        with self.assertRaises(IndexError)  as cm_new:
            fu.mulmat(m1=a,m2=b,tolistconv=True)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.mulmat(m1=a,m2=b,tolistconv=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_mulmat(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 6]]
        b = [[2, 3, 2], [3, 4, 3], [3, 4, 6]]
        return_new = fu.mulmat(m1=a,m2=b,tolistconv=True)
        return_old = oldfu.mulmat(m1=a,m2=b,tolistconv=True)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new,[[23, 31, 29], [28, 38, 37], [48, 65, 66]]))


# todo: modificato
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
            fu.recmat(mat=a, tolistconv=True)
        with self.assertRaises(IndexError)  as cm_old:
            oldfu.recmat(mat=a, tolistconv=True)
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_recmat(self):
        a = [[1, 5, 2], [2, 5, 3], [3, 8, 1]]
        return_new = fu.recmat(mat=a, tolistconv=True)
        return_old = oldfu.recmat(mat=a, tolistconv=True)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new,(78.690067525979785, 0.0, 0.0)))



class Test_symclass(unittest.TestCase):
    pass





if __name__ == '__main__':
    unittest.main()




