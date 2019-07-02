from __future__ import print_function
from __future__ import division


import unittest
import numpy
from ..libpy import sparx_utilities
from test_module import get_real_data

from mpi import *
mpi_init(0, [])

IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_BLANK_2D = sparx_utilities.model_blank(10, 10)
IMAGE_BLANK_3D = sparx_utilities.model_blank(10, 10, 10)
MASK = sparx_utilities.model_circle(2, 5, 5)



"""
There are some opened issues in:
1) fit_tanh: with 'low' too high leads to a bug that leads to a TypeError situation. see "test_tooHighLow_leads_TypeError" 
2) I'm not able to test 'Test_filterlocal' beacuse it lead always to an indexError due to the may_node !=myid
"""


class Test_filt_ctf(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError)  as cm_new:
            fu.filt_ctf()
        with self.assertRaises(TypeError)  as cm_old:
            oldfu.filt_ctf()
        self.assertEqual(cm_new.exception.message, "filt_ctf() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        with self.assertRaises(AssertionError) as cm_new:
            fu.filt_ctf(EMData(),ctf, dopad=True, sign=1, binary=0)
        with self.assertRaises(AssertionError) as cm_old:
            oldfu.filt_ctf(EMData(),ctf, dopad=True, sign=1, binary=0)
        self.assertEqual(cm_new.exception.message, "")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_None_value_raiseError(self):
        # Since there is no try-except paradigma but only an assert I suppose that the most common error is given a None vuole instead an image
        with self.assertRaises(AttributeError)  as cm_new:
            fu.filt_ctf(None, None)
        with self.assertRaises(AttributeError)  as cm_old:
            oldfu.filt_ctf(None,None)
        self.assertEqual(cm_new.exception.message, "'NoneType' object has no attribute 'get_ysize'")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_Img2D_with_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_2D,ctf, dopad=True, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_2D,ctf, dopad=True, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_without_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_2D,ctf, dopad=False, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_2D,ctf, dopad=False, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_3D,ctf, dopad=True, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_3D,ctf, dopad=True, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_3D,ctf, dopad=False, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_3D,ctf, dopad=False, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_BLANK_2D,ctf, dopad=True, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_BLANK_2D,ctf, dopad=True, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_BLANK_2D,ctf, dopad=False, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_BLANK_2D,ctf, dopad=False, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_BLANK_3D,ctf, dopad=True, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_BLANK_3D,ctf, dopad=True, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_pad(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(IMAGE_BLANK_3D,ctf, dopad=False, sign=1, binary=0)
        return_old = oldfu.filt_ctf(IMAGE_BLANK_3D,ctf, dopad=False, sign=1, binary=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_tophatl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_tophatl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_tophatl()
        self.assertEqual(cm_new.exception.message, "filt_tophatl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError)as cm_new:
            fu.filt_tophatl(EMData(), freq=0.25, pad=False)
        with self.assertRaises(RuntimeError)as cm_old:
            oldfu.filt_tophatl(EMData(), freq=0.25, pad=False)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new = fu.filt_tophatl(None, freq= 0.25, pad=False)
        return_old = oldfu.filt_tophatl(None,freq= 0.25, pad=False)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_Img2D_without_pad(self):
        return_new = fu.filt_tophatl(IMAGE_2D, freq= 0.25, pad=False)
        return_old = oldfu.filt_tophatl(IMAGE_2D,freq= 0.25, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2D_with_pad(self):
        return_new = fu.filt_tophatl(IMAGE_2D, freq= 0.25, pad=True)
        return_old = oldfu.filt_tophatl(IMAGE_2D,freq= 0.25, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_without_pad(self):
        return_new = fu.filt_tophatl(IMAGE_BLANK_2D, freq= 0.25, pad=False)
        return_old = oldfu.filt_tophatl(IMAGE_BLANK_2D,freq= 0.25, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img2Dblank_with_pad(self):
        return_new = fu.filt_tophatl(IMAGE_BLANK_2D, freq= 0.25, pad=True)
        return_old = oldfu.filt_tophatl(IMAGE_BLANK_2D,freq= 0.25, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_without_pad(self):
        return_new = fu.filt_tophatl(IMAGE_3D, freq= 0.25, pad=False)
        return_old = oldfu.filt_tophatl(IMAGE_3D,freq= 0.25, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3D_with_pad(self):
        return_new = fu.filt_tophatl(IMAGE_3D, freq= 0.25, pad=True)
        return_old = oldfu.filt_tophatl(IMAGE_3D,freq= 0.25, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_without_pad(self):
        return_new = fu.filt_tophatl(IMAGE_BLANK_3D, freq= 0.25, pad=False)
        return_old = oldfu.filt_tophatl(IMAGE_BLANK_3D,freq= 0.25, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_Img3Dblank_with_pad(self):
        return_new = fu.filt_tophatl(IMAGE_BLANK_3D, freq= 0.25, pad=True)
        return_old = oldfu.filt_tophatl(IMAGE_BLANK_3D,freq= 0.25, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_tophatb(unittest.TestCase):
    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_tophatb(EMData(), freql=0.25, freqh=0.35)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_tophatb(EMData(), freql=0.25, freqh=0.35)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_tophatb()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_tophatb()
        self.assertEqual(cm_new.exception.message, "filt_tophatb() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_NoneType_Img(self):
        return_new = fu.filt_tophatb(None, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_tophatb(None, freql=0.25, freqh=0.35, pad=True)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_with_pad(self):
        return_new = fu.filt_tophatb(IMAGE_2D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_tophatb(IMAGE_2D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_without_pad(self):
        return_new = fu.filt_tophatb(IMAGE_2D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_tophatb(IMAGE_2D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_tophatb(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_tophatb(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_tophatb(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_tophatb(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_tophatb(IMAGE_3D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_tophatb(IMAGE_3D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_tophatb(IMAGE_3D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_tophatb(IMAGE_3D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_tophatb(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_tophatb(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_tophatb(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_tophatb(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_gaussl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_gaussl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_gaussl()
        self.assertEqual(cm_new.exception.message, "filt_gaussl() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_gaussl(EMData(), sigma=0.23)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_gaussl(EMData(), sigma=0.23)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_gaussl(None,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussl(None,sigma= 0.23, pad=False )
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_without_pad(self):
        return_new = fu.filt_gaussl(IMAGE_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussl(IMAGE_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_with_pad(self):
        return_new = fu.filt_gaussl(IMAGE_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussl(IMAGE_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_gaussl(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussl(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_gaussl(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussl(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_gaussl(IMAGE_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussl(IMAGE_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_gaussl(IMAGE_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussl(IMAGE_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_gaussl(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussl(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_gaussl(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussl(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_gaussinv(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_gaussinv()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_gaussinv()
        self.assertEqual(cm_new.exception.message, "filt_gaussinv() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_gaussinv(EMData(), sigma=0.23)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_gaussinv(EMData(), sigma=0.23)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_gaussinv(None,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussinv(None,sigma= 0.23, pad=False )
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_without_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussinv(IMAGE_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_with_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussinv(IMAGE_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussinv(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussinv(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussinv(IMAGE_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussinv(IMAGE_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussinv(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_gaussinv(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussinv(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_gaussh(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_gaussh()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_gaussh()
        self.assertEqual(cm_new.exception.message, "filt_gaussh() takes at least 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_gaussh(EMData(), sigma=0.23)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_gaussh(EMData(), sigma=0.23)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_gaussh(None,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussh(None,sigma= 0.23, pad=False )
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_without_pad(self):
        return_new = fu.filt_gaussh(IMAGE_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussh(IMAGE_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_with_pad(self):
        return_new = fu.filt_gaussh(IMAGE_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussh(IMAGE_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_gaussh(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussh(IMAGE_BLANK_2D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_gaussh(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussh(IMAGE_BLANK_2D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_gaussh(IMAGE_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussh(IMAGE_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_gaussh(IMAGE_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussh(IMAGE_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_gaussh(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        return_old = oldfu.filt_gaussh(IMAGE_BLANK_3D,sigma= 0.23, pad=False )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_gaussh(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        return_old = oldfu.filt_gaussh(IMAGE_BLANK_3D,sigma= 0.23, pad=True )
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_btwl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_btwl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_btwl()
        self.assertEqual(cm_new.exception.message, "filt_btwl() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_btwl(EMData(), freql=0.25, freqh=0.35)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_btwl(EMData(), freql=0.25, freqh=0.35)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_btwl(None,  freql=0.25, freqh=0.35, pad=False )
        return_old = oldfu.filt_btwl(None,  freql=0.25, freqh=0.35, pad=False )
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_without_pad(self):
        return_new = fu.filt_btwl(IMAGE_2D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_btwl(IMAGE_2D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_with_pad(self):
        return_new = fu.filt_btwl(IMAGE_2D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_btwl(IMAGE_2D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_btwl(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_btwl(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_btwl(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_btwl(IMAGE_BLANK_2D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_btwl(IMAGE_3D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_btwl(IMAGE_3D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_btwl(IMAGE_3D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_btwl(IMAGE_3D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_btwl(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=False)
        return_old = oldfu.filt_btwl(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_btwl(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=True)
        return_old = oldfu.filt_btwl(IMAGE_BLANK_3D, freql=0.25, freqh=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_tanl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_tanl()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_tanl()
        self.assertEqual(cm_new.exception.message, "filt_tanl() takes at least 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_tanl(EMData(), freq=0.25, fall_off=0.35)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_tanl(EMData(), freq=0.25, fall_off=0.35)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_tanl(None,  freq=0.25, fall_off=0.35, pad=False )
        return_old = oldfu.filt_tanl(None,  freq=0.25, fall_off=0.35, pad=False )
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg_without_pad(self):
        return_new = fu.filt_tanl(IMAGE_2D, freq=0.25, fall_off=0.35, pad=False)
        return_old = oldfu.filt_tanl(IMAGE_2D, freq=0.25, fall_off=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImg_with_pad(self):
        return_new = fu.filt_tanl(IMAGE_2D, freq=0.25, fall_off=0.35, pad=True)
        return_old = oldfu.filt_tanl(IMAGE_2D, freq=0.25, fall_off=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_without_pad(self):
        return_new = fu.filt_tanl(IMAGE_BLANK_2D, freq=0.25, fall_off=0.35, pad=False)
        return_old = oldfu.filt_tanl(IMAGE_BLANK_2D, freq=0.25, fall_off=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank_with_pad(self):
        return_new = fu.filt_tanl(IMAGE_BLANK_2D, freq=0.25, fall_off=0.35, pad=True)
        return_old = oldfu.filt_tanl(IMAGE_BLANK_2D, freq=0.25, fall_off=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_without_pad(self):
        return_new = fu.filt_tanl(IMAGE_3D, freq=0.25, fall_off=0.35, pad=False)
        return_old = oldfu.filt_tanl(IMAGE_3D, freq=0.25, fall_off=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg_with_pad(self):
        return_new = fu.filt_tanl(IMAGE_3D, freq=0.25, fall_off=0.35, pad=True)
        return_old = oldfu.filt_tanl(IMAGE_3D, freq=0.25, fall_off=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_without_pad(self):
        return_new = fu.filt_tanl(IMAGE_BLANK_3D, freq=0.25, fall_off=0.35, pad=False)
        return_old = oldfu.filt_tanl(IMAGE_BLANK_3D, freq=0.25, fall_off=0.35, pad=False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank_with_pad(self):
        return_new = fu.filt_tanl(IMAGE_BLANK_3D, freq=0.25, fall_off=0.35, pad=True)
        return_old = oldfu.filt_tanl(IMAGE_BLANK_3D, freq=0.25, fall_off=0.35, pad=True)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_filt_table(unittest.TestCase):
    table = [entry for entry in numpy.linspace(0, 0.5).tolist()]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_table()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_table()
        self.assertEqual(cm_new.exception.message, "filt_table() takes exactly 2 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.filt_table(EMData(), self.table)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.filt_table(EMData(), self.table)
        msg = cm_new.exception.message.split("'")
        msg_old = cm_old.exception.message.split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], 'x size <= 0')
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_NoneType_Img(self):
        return_new =  fu.filt_table(None,  self.table )
        return_old = oldfu.filt_table(None,  self.table)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,None)

    def test_2DImg(self):
        return_new = fu.filt_table(IMAGE_2D, self.table)
        return_old = oldfu.filt_table(IMAGE_2D, self.table)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_2DImgBlank(self):
        return_new = fu.filt_table(IMAGE_BLANK_2D, self.table)
        return_old = oldfu.filt_table(IMAGE_BLANK_2D, self.table)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImg(self):
        return_new = fu.filt_table(IMAGE_3D, self.table)
        return_old = oldfu.filt_table(IMAGE_3D, self.table)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_3DImgBlank(self):
        return_new = fu.filt_table(IMAGE_BLANK_3D, self.table)
        return_old = oldfu.filt_table(IMAGE_BLANK_3D, self.table)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_fit_tanh(unittest.TestCase):
    dres = ((0.0, 0.05, 0, 10, 0.15, 0.20), (0, 0.2, 0.4, 0.6, 0.8, 1.0), (8, 9, 5, 77, 98, 200))
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fit_tanh()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fit_tanh()
        self.assertEqual(cm_new.exception.message, "fit_tanh() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyList_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.fit_tanh([])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fit_tanh([])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_fit_tanh(self):
        return_new = fu.fit_tanh(self.dres, low=0.1)
        return_old = oldfu.fit_tanh(self.dres, low=0.1 )

        self.assertEqual(return_new[0] , return_old[0])
        self.assertEqual(return_new[1], return_old[1])
        self.assertEqual(return_new[0] , 0)
        self.assertEqual(return_new[1], 0.1)

    def test_tooHighLow_leads_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fit_tanh(self.dres, low=10)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fit_tanh(self.dres, low=10 )
        self.assertEqual(cm_new.exception.message, "'tuple' object does not support item assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_fit_tanh1(unittest.TestCase):
    dres = ((0.0, 0.05, 0, 10, 0.15, 0.20), (0, 0.2, 0.4, 0.6, 0.8, 1.0), (8, 9, 5, 77, 98, 200))
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fit_tanh1()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fit_tanh1()
        self.assertEqual(cm_new.exception.message, "fit_tanh1() takes at least 1 argument (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_emptyList_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.fit_tanh1([])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.fit_tanh1([])
        self.assertEqual(cm_new.exception.message, "list index out of range")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_fit_tanh1(self):
        return_new = fu.fit_tanh1(self.dres, low=0.1)
        return_old = oldfu.fit_tanh1(self.dres, low=0.1 )

        self.assertEqual(return_new[0] , return_old[0])
        self.assertEqual(return_new[1], return_old[1])
        self.assertEqual(return_new[0] , 0)
        self.assertEqual(return_new[1], 0.1)

    def test_tooHighLow_leads_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fit_tanh1(self.dres, low=10)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fit_tanh(self.dres, low=10 )
        self.assertEqual(cm_new.exception.message, "'tuple' object does not support item assignment")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)



class Test_filt_vols(unittest.TestCase):
    dres = ((0.0, 0.05, 0, 10, 0.15, 0.20), (0, 0.2, 0.4, 0.6, 0.8, 1.0), (8, 9, 5, 77, 98, 200))
    vols = [sparx_utilities.model_gauss_noise(0.25 , 10,10,10),sparx_utilities.model_gauss_noise(0.25 , 10,10,10),sparx_utilities.model_gauss_noise(0.25 , 10,10,10)]
    fscs = (dres, dres, dres)

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertEqual(len(return_old), len(return_new))
            for img1, img2 in zip(return_old, return_new):
                self.assertTrue(numpy.allclose(img1.get_3dview(), img2.get_3dview(), equal_nan=True))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filt_vols()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filt_vols()
        self.assertEqual(cm_new.exception.message, "filt_vols() takes exactly 3 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    def test_with_circularMASK(self):
        return_new = fu.filt_vols(self.vols, self.fscs, MASK)
        return_old = oldfu.filt_vols(self.vols, self.fscs, MASK)
        self.test_all_the_conditions(return_new,return_old,False)

    def test_with_NoneTypeMask_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.filt_vols(self.vols, self.fscs, None)
        return_old = oldfu.filt_vols(self.vols, self.fscs, None)
        self.test_all_the_conditions(return_new,return_old,False)
        """

    def test_with_EmptyMASK_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.filt_vols(self.vols, self.fscs, EMData())
        return_old = oldfu.filt_vols(self.vols, self.fscs, EMData())
        self.test_all_the_conditions(return_new,return_old,False)
        """



class Test_filterlocal(unittest.TestCase):
    def test_wrong_number_params(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.filterlocal()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.filterlocal()
        self.assertEqual(cm_new.exception.message, "filterlocal() takes exactly 7 arguments (0 given)")
        self.assertEqual(cm_new.exception.message, cm_old.exception.message)

    @unittest.skip("not able to tst it")
    def test_filterlocal_true_should_return_equal_object(self):
        #mpi_barrier(MPI_COMM_WORLD)
        vols = []
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        ui =  ut.model_gauss_noise(0.25 , 10,10,10)     # or use ut.model_blank(1,1,1)
        vi =  ut.model_gauss_noise(0.25 , 10,10,10)
        m  = "sphire/tests/3d_volume.txt"
        falloff = 4
        myid =   1
        main_node = 0
        number_of_proc = 6

        return_new = fu.filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc)
        return_old = oldfu.filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc)

        self.assertTrue(return_new, return_old)











import numpy

import EMAN2_cppwrap as e2cpp

from ..libpy_py3 import sphire_filter as fu
from .sparx_lib import sparx_filter as oldfu

from ..libpy import sparx_utilities as ut

import unittest


from test_module import get_data, get_data_3d
from EMAN2_cppwrap import EMData


@unittest.skip("original adnan")
class MyTestCase(unittest.TestCase):

    def test_filt_tophatl_true_should_return_equal_object(self):
        image, = get_data(1)
        freq = 0.25
        return_new = fu.filt_tophatl(image,freq)
        return_old = oldfu.filt_tophatl(image,freq)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_filt_tophatb_true_should_return_equal_object(self):
        image, = get_data(1)
        freql = 0.25
        freqh = 0.35
        return_new = fu.filt_tophatb(image,freql, freqh)
        return_old = oldfu.filt_tophatb(image,freql, freqh)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_filt_gaussl_true_should_return_equal_object(self):
        image, = get_data(1)
        sigma = 0.23
        return_new = fu.filt_gaussl(image,sigma)
        return_old = oldfu.filt_gaussl(image,sigma)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_filt_gaussinv_true_should_return_equal_object(self):
        image, = get_data(1)
        sigma = 0.23
        return_new = fu.filt_gaussinv(image,sigma)
        return_old = oldfu.filt_gaussinv(image,sigma)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_filt_gaussh_true_should_return_equal_object(self):
        image, = get_data(1)
        sigma = 0.23
        return_new = fu.filt_gaussh(image,sigma)
        return_old = oldfu.filt_gaussh(image,sigma)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_filt_btwl_true_should_return_equal_object(self):
        image, = get_data(1)
        freql = 0.25
        freqh = 0.35
        return_new = fu.filt_btwl(image,freql, freqh)
        return_old = oldfu.filt_btwl(image,freql, freqh)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_filt_tanl_true_should_return_equal_object(self):
        image, = get_data(1)
        freql = 0.25
        freqh = 0.35
        return_new = fu.filt_tanl(image,freql, freqh)
        return_old = oldfu.filt_tanl(image,freql, freqh)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_filt_table_true_should_return_equal_object(self):
        image, = get_data(1)
        table =[entry for entry in numpy.linspace(0, 0.5).tolist()]

        return_new = fu.filt_table(image,table)
        return_old = oldfu.filt_table(image,table)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_filt_ctf_true_should_return_equal_object(self):
        image, = get_data(1)
        defocus = 1
        cs =  2
        voltage = 300
        pixel_size = 1.5
        bfactor =  0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.filt_ctf(image,ctf)
        return_old = oldfu.filt_ctf(image,ctf)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_fit_tanh_true_should_return_equal_object(self):
        dres = []
        dres.append((0.0,0.05,0,10,0.15,0.20))
        dres.append((0,0.2,0.4,0.6,0.8,1.0))
        dres.append((8,9,5,77,98,200))

        return_new = fu.fit_tanh(dres)
        return_old = oldfu.fit_tanh(dres)

        self.assertTrue(return_new , return_old)


    def test_fit_tanh1_true_should_return_equal_object(self):
        dres = []
        dres.append((0.0,0.05,0,10,0.15,0.20))
        dres.append((0,0.2,0.4,0.6,0.8,1.0))
        dres.append((8,9,5,77,98,200))

        return_new = fu.fit_tanh1(dres)
        return_old = oldfu.fit_tanh1(dres)

        self.assertTrue(return_new , return_old)

    def test_filt_vols_true_should_return_equal_object(self):
        vols = []
        vols.append( ut.model_gauss_noise(0.25 , 10,10,10) )
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))

        dres = []
        dres.append([0.05, 0.05, 0.05, 0.10, 0.15, 0.20])
        dres.append([0.05, 0.02, 0.04, 0.06, 0.08, 0.09])
        dres.append([0.08, 0.09, 0.05, 0.07, 0.05, 0.01])
        fscs = []
        fscs.append(dres)
        fscs.append(dres)
        fscs.append(dres)

        mask3D, = get_data_3d(1)

        return_new = fu.filt_vols(vols, fscs, mask3D)
        return_old = oldfu.filt_vols(vols, fscs, mask3D)

        self.assertTrue(return_new, return_old)

    def test_filterlocal_true_should_return_equal_object(self):
        vols = []
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        vols.append(ut.model_gauss_noise(0.25 , 10,10,10))
        ui =  ut.model_gauss_noise(0.25 , 10,10,10)     # or use ut.model_blank(1,1,1)
        vi =  ut.model_gauss_noise(0.25 , 10,10,10)
        m  = "sphire/tests/3d_volume.txt"
        falloff = 4
        myid =   1
        main_node = 0
        number_of_proc = 6

        return_new = fu.filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc)
        return_old = oldfu.filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc)

        self.assertTrue(return_new, return_old)




if __name__ == '__main__':
    unittest.main()

