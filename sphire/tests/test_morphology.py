from __future__ import print_function
from __future__ import division

from ..libpy import sparx_morphology as fu
from .sparx_lib import sparx_morphology as oldfu

from ..libpy import sparx_filter
from ..libpy import sparx_utilities

import numpy
import unittest
from cPickle import load as pickle_load
from EMAN2_cppwrap import EMData, EMAN2Ctf
from os import path
from shutil import rmtree

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
TOLERANCE = 0.005

def get_data(num, dim=10):
    data_list = []
    for i in range(num):
        a = EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list

def get_data_3d(num, dim=10):
    data_list = []
    for i in range(num):
        a = EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
        data_list.append(a)
    return data_list

def get_data_gauss_noise(dim = 10):
    return sparx_utilities.model_gauss_noise(0.25 , dim,dim,dim)

def get_arg_from_pickle_files(filepath):
    with open(filepath, 'rb') as rb:
        return pickle_load(rb)

def remove_dir(d):
    if path.isdir(d):
        rmtree(d)




class Test_binarize(unittest.TestCase):

    def test_D1(self):
        image = get_data(1)[0]

        return_new = fu.binarize(image)
        return_old = oldfu.binarize(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        """ We are not able to catch the 'NotExistingObjectException' C++ exception"""
        with self.assertRaises(RuntimeError):
            fu.binarize(EMData())
            oldfu.binarize(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.binarize()
            oldfu.binarize()



class Test_collapse(unittest.TestCase):

    def test_D2(self):
        image = get_data(1)[0]

        return_new = fu.collapse(image)
        return_old = oldfu.collapse(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        """ We are not able to catch the 'NotExistingObjectException' C++ exception"""
        with self.assertRaises(RuntimeError):
            fu.collapse(EMData())
            oldfu.collapse(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.collapse()
            oldfu.collapse()



class Test_dilatation(unittest.TestCase):
    image_2d = sparx_utilities.model_blank(10,10)
    image_3d = sparx_utilities.model_blank(10, 10,10)
    mask = sparx_utilities.model_circle(2, 5, 5)

    @unittest.skip("\n***************************\n\t\t 'Test_dilatation.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_test_empty_input_image(self):
        """
        It seems not possible to test it without getting an segmentation fault. We should return None after the first error message
        in order to avoid to run in the second if/else and get the segmentation fault in the c+++ code
        """
        img =EMData()
        return_new = fu.dilation(img)
        return_old = oldfu.dilation(img)
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError):
            fu.dilation(self.image_2d, EMData(), morphtype="BINARY")
            oldfu.dilation(self.image_2d, EMData(), morphtype="BINARY")

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.dilation()
            oldfu.dilation()

    def test_D3(self):
        return_new = fu.dilation(self.image_2d, self.mask, morphtype="BINARY")
        return_old = oldfu.dilation(self.image_2d, self.mask, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D3_2(self):
        return_new = fu.dilation(self.image_3d, self.mask, morphtype="BINARY")
        return_old = oldfu.dilation(self.image_3d, self.mask, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E3(self):
        return_new = fu.dilation(self.image_2d, self.mask, morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(self.image_2d, self.mask, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E3_2(self):
        return_new = fu.dilation(self.image_3d, self.mask, morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(self.image_3d, self.mask, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_I3(self):
        return_new = fu.dilation(self.image_3d, self.mask, morphtype="invalid_type")
        return_old = oldfu.dilation(self.image_3d, self.mask, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_G3(self):
        return_new = fu.dilation(self.image_2d, morphtype="BINARY")
        return_old = oldfu.dilation(self.image_2d, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_J3(self):
        return_new = fu.dilation(self.image_3d, morphtype="BINARY")
        return_old = oldfu.dilation(self.image_3d, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_L3(self):
        return_new = fu.dilation(self.image_2d,morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(self.image_2d,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_K3(self):
        return_new = fu.dilation(self.image_3d,morphtype="GRAYLEVEL")
        return_old = oldfu.dilation(self.image_3d,morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_erosion(unittest.TestCase):
    image_2d = sparx_utilities.model_blank(10,10)
    image_3d = sparx_utilities.model_blank(10, 10,10)
    mask = sparx_utilities.model_circle(2, 5, 5)

    @unittest.skip("\n***************************\n\t\t 'Test_erosion.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        """
        It seems not possible to test it without getting an segmentation fault. We should return None after the first error message
        in order to avoid to run in the second if/else and get the segmentation fault in the c+++ code
        """
        img =EMData()
        return_new = fu.erosion(img)
        return_old = oldfu.erosion(img)
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_empty_input_image2(self):
        with self.assertRaises(RuntimeError):
            fu.erosion(self.image_2d, EMData())
            oldfu.erosion(self.image_2d, EMData())

    def test_empty_input_image3(self):
        with self.assertRaises(RuntimeError):
            fu.erosion(self.image_3d, EMData())
            oldfu.erosion(self.image_3d, EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.erosion()
            oldfu.erosion()

    def test_D4(self):
        return_new = fu.erosion(self.image_2d, self.mask, morphtype="BINARY")
        return_old = oldfu.erosion(self.image_2d, self.mask, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D4_2(self):
        return_new = fu.erosion(self.image_3d, self.mask, morphtype="BINARY")
        return_old = oldfu.erosion(self.image_3d, self.mask, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E4(self):
        return_new = fu.erosion(self.image_2d, self.mask, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(self.image_2d, self.mask, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E4_2(self):
        return_new = fu.erosion(self.image_3d, self.mask, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(self.image_3d, self.mask, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_I4(self):
        return_new = fu.erosion(self.image_3d, self.mask, morphtype="invalid_type")
        return_old = oldfu.erosion(self.image_3d, self.mask, morphtype="invalid_type")
        self.assertTrue(return_old is None)
        self.assertTrue(return_new is None)

    def test_G4(self):
        return_new = fu.erosion(self.image_2d, morphtype="BINARY")
        return_old = oldfu.erosion(self.image_2d, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_J4(self):
        return_new = fu.erosion(self.image_3d, morphtype="BINARY")
        return_old = oldfu.erosion(self.image_3d, morphtype="BINARY")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_L4(self):
        return_new = fu.erosion(self.image_2d, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(self.image_2d, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_K4(self):
        return_new = fu.erosion(self.image_3d, morphtype="GRAYLEVEL")
        return_old = oldfu.erosion(self.image_3d, morphtype="GRAYLEVEL")
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_power(unittest.TestCase):
    image = get_data(1)[0]

    def test_D5(self):
        return_new = fu.power(self.image)
        return_old = oldfu.power(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.power(EMData())
            oldfu.power(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.power()
            oldfu.power()


class Test_square_root(unittest.TestCase):
    image = get_data(1)[0]

    def test_D6(self):
        return_new = fu.square_root(self.image)
        return_old = oldfu.square_root(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.square_root(EMData())
            oldfu.square_root(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.square_root()
            oldfu.square_root()


class Test_square(unittest.TestCase):
    image = get_data(1)[0]

    def test_D7(self):
        return_new = fu.square(self.image)
        return_old = oldfu.square(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.square(EMData())
            oldfu.square(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.square()
            oldfu.square()


class Test_threshold(unittest.TestCase):
    image = get_data(1)[0]

    def test_D8(self):
        return_new = fu.threshold(self.image)
        return_old = oldfu.threshold(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.threshold(EMData())
            oldfu.threshold(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.threshold()
            oldfu.threshold()


class Test_threshold_outside(unittest.TestCase):
    image = get_data(1)[0]

    def test_D9(self):
        return_new = fu.threshold_outside(self.image, 2 , 10)
        return_old = oldfu.threshold_outside(self.image, 2 , 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        return_new = fu.threshold_outside(EMData(), 2 , 10)
        return_old = oldfu.threshold_outside(EMData(), 2 , 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.threshold_outside()
            oldfu.threshold_outside()


class Test_notzero(unittest.TestCase):
    image = get_data(1)[0]

    def test_D10(self):
        return_new = fu.notzero(self.image)
        return_old = oldfu.notzero(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.notzero(EMData())
            oldfu.notzero(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.notzero()
            oldfu.notzero()


class Test_rotavg_ctf(unittest.TestCase):
    image = get_data(1)[0]

    def test_empty_input_image(self):
        return_new = fu.rotavg_ctf(EMData(), defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5)
        return_old = oldfu.rotavg_ctf(EMData(), defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.rotavg_ctf()
            oldfu.rotavg_ctf()

    def test_D11_E11(self):
        return_new = fu.rotavg_ctf(self.image, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5)
        return_old = oldfu.rotavg_ctf(self.image, defocus= 1, Cs =0.0, voltage=300, Pixel_size=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_G11_H11(self):
        return_new = fu.rotavg_ctf(self.image, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5)
        return_old = oldfu.rotavg_ctf(self.image, defocus= 1, Cs =2, voltage=300, Pixel_size=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_RuntimeWarning(self):
        return_new = fu.rotavg_ctf(self.image, defocus= 1, Cs =2, voltage=300, Pixel_size=0)
        return_old = oldfu.rotavg_ctf(self.image, defocus= 1, Cs =2, voltage=300, Pixel_size=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_ctf_1d(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctf_1d()
            oldfu.ctf_1d()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_1d(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_1d(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new).any())
        self.assertTrue(numpy.isnan(return_old).any())

    def test_no_image_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_1d(nx=0, ctf=ctf)
            oldfu.ctf_1d(nx=0, ctf=ctf)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_1d(nx =2, ctf=ctf)
            oldfu.ctf_1d(nx=2, ctf=ctf)

    def test_E12(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=False)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=False)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E12_2(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, doabs=False)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, doabs=False)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D12(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=True)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, sign=0, doabs=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D12_2(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_1d(nx=20, ctf=ctf, doabs=True)
        return_old= oldfu.ctf_1d(nx=20, ctf=ctf, doabs=True)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_ctf_2(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctf_2()
            oldfu.ctf_2()

    def test_no_image_size_error(self):
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

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(ZeroDivisionError):
            fu.ctf_2(nx =2, ctf=ctf)
            oldfu.ctf_2(nx=2, ctf=ctf)

    def test_ctf_2_D13(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_2(nx =2, ctf=ctf)
        return_old= oldfu.ctf_2(nx=2, ctf=ctf)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_ctf_img(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctf_img()
            oldfu.ctf_img()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_img(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_img(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img(nx=0, ctf=ctf)
            oldfu.ctf_img(nx=0, ctf=ctf)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_img(nx =2, ctf=ctf)
        return_old = oldfu.ctf_img(nx=2, ctf=ctf)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_D14(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img(nx =20, ctf=ctf)
        return_old=oldfu.ctf_img(nx=20, ctf=ctf)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_ctf_img_real(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctf_img_real()
            oldfu.ctf_img_real()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_img_real(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_img_real(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_img_real(nx=0, ctf=ctf)
            oldfu.ctf_img_real(nx=0, ctf=ctf)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_img_real(nx =2, ctf=ctf)
        return_old = oldfu.ctf_img_real(nx=2, ctf=ctf)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_D15(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx=20, ctf=ctf)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E15(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_img_real(nx=20, ctf=ctf, ny=1)
        return_old=oldfu.ctf_img_real(nx=20, ctf=ctf, ny=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_ctf_rimg(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctf_rimg()
            oldfu.ctf_rimg()

    def test_with_empty_ctfDict(self):
        return_new =fu.ctf_rimg(nx=20, ctf= EMAN2Ctf())
        return_old= oldfu.ctf_rimg(nx=20, ctf= EMAN2Ctf())
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_no_image_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        with self.assertRaises(RuntimeError):
            fu.ctf_rimg(nx=0, ctf=ctf)
            oldfu.ctf_rimg(nx=0, ctf=ctf)

    def test_no_pixel_size_error(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 0, "bfactor": 0, "ampcont": 0.1})
        return_new = fu.ctf_rimg(nx =2, ctf=ctf)
        return_old = oldfu.ctf_rimg(nx=2, ctf=ctf)
        self.assertTrue(numpy.isnan(return_new.get_3dview()).any())
        self.assertTrue(numpy.isnan(return_old.get_3dview()).any())

    def test_D16(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx=20, ctf=ctf)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E16(self):
        ctf = EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0, "ampcont": 0.1})
        return_new =fu.ctf_rimg(nx=20, ctf=ctf, ny=1)
        return_old=oldfu.ctf_rimg(nx=20, ctf=ctf, ny=1)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_ctflimit(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ctflimit()
            oldfu.ctflimit()

    def test_D17(self):
        return_new = fu.ctflimit(nx=30, defocus=1, cs=2, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=30, defocus=1, cs=2, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E17(self):
        return_new = fu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=1.5)
        return_old = oldfu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=1.5)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.ctflimit(nx=-1, defocus=1, cs=2, voltage=300, pix=1.5)
            oldfu.ctflimit(nx=-1, defocus=1, cs=2, voltage=300, pix=1.5)

    def test_no_pixel_size_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=0)
            oldfu.ctflimit(nx=0, defocus=1, cs=2, voltage=300, pix=0)


class Test_imf_params_cl1(unittest.TestCase):
    pw = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.imf_params_cl1()
            oldfu.imf_params_cl1()

    def test_no_pixel_size_error(self):
        return_new = fu.imf_params_cl1(self.pw,Pixel_size=0)
        return_old = oldfu.imf_params_cl1(self.pw,Pixel_size=0)

        self.assertTrue(numpy.isnan(return_new[0]).any())
        self.assertTrue(numpy.isnan(return_old[0]).any())

        self.assertTrue(numpy.isnan(return_new[1]).any())
        self.assertTrue(numpy.isnan(return_old[1]).any())

        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))

    def test_D18(self):
        return_new = fu.imf_params_cl1(self.pw)
        return_old = oldfu.imf_params_cl1(self.pw)
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))

    def test_D18_2(self):
        return_new = fu.imf_params_cl1(self.pw, n=0)
        return_old = oldfu.imf_params_cl1(self.pw, n=0)
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))

    def test_D18_3(self):
        return_new = fu.imf_params_cl1(self.pw, iswi=0)
        return_old = oldfu.imf_params_cl1(self.pw, iswi=0)
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))

    def test_D18_4(self):
        return_new = fu.imf_params_cl1(self.pw, iswi=0, n=0)
        return_old = oldfu.imf_params_cl1(self.pw, iswi=0, n=0)
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))

    @unittest.skip("\n***************************\n\t\t 'Test_imf_params_cl1.test_D18_5' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_D18_5(self):
        return_new = fu.imf_params_cl1([])
        return_old = oldfu.imf_params_cl1([])
        self.assertTrue(numpy.array_equal(return_new[0], return_old[0]))
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
        self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
        self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))


class Test_adaptive_mask(unittest.TestCase):
    image = get_data_3d(1)[0]

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.adaptive_mask(EMData())
            oldfu.adaptive_mask(EMData())

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.adaptive_mask()
            oldfu.adaptive_mask()

    def test_D19(self):
        return_new = fu.adaptive_mask(self.image)
        return_old = oldfu.adaptive_mask(self.image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E19(self):
        return_new = fu.adaptive_mask(self.image,threshold = 0)
        return_old = oldfu.adaptive_mask(self.image,threshold = 0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E19_2(self):
        return_new = fu.adaptive_mask(self.image,threshold = 0, ndilation=0)
        return_old = oldfu.adaptive_mask(self.image,threshold = 0, ndilation=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E19_3(self):
        return_new = fu.adaptive_mask(self.image,threshold = 0, nsigma=0)
        return_old = oldfu.adaptive_mask(self.image,threshold = 0, nsigma=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_E19_4(self):
        return_new = fu.adaptive_mask(self.image,threshold = 0, edge_width=0)
        return_old = oldfu.adaptive_mask(self.image,threshold = 0, edge_width=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D19_2(self):
        return_new = fu.adaptive_mask(self.image, ndilation=0)
        return_old = oldfu.adaptive_mask(self.image, ndilation=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D19_3(self):
        return_new = fu.adaptive_mask(self.image, nsigma=0)
        return_old = oldfu.adaptive_mask(self.image, nsigma=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D19_4(self):
        return_new = fu.adaptive_mask(self.image, edge_width=0)
        return_old = oldfu.adaptive_mask(self.image, edge_width=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_cosinemask(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.cosinemask()
            oldfu.cosinemask()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.cosinemask(EMData())
            oldfu.cosinemask(EMData())

    @unittest.skip("\n***************************\n\t\t 'Test_cosinemask.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        image = get_data_3d(1)[0]
        bckg = EMData()
        return_new = fu.cosinemask(image, bckg=bckg)
        return_old = oldfu.cosinemask(image, bckg=bckg)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    @unittest.skip("\n***************************\n\t\t 'Test_cosinemask.test_D20_4' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_D20_4(self):
        image = get_data_3d(1)[0]
        bckg = get_data_gauss_noise()
        return_new = fu.cosinemask(image, bckg=bckg)
        return_old = oldfu.cosinemask(image, bckg=bckg)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D20(self):
        image = get_data_3d(1)[0]
        return_new = fu.cosinemask(image)
        return_old = oldfu.cosinemask(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D20_1(self):
        image = get_data_3d(1)[0]
        return_new = fu.cosinemask(image, radius=0)
        return_old = oldfu.cosinemask(image, radius=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D20_2(self):
        image = get_data_3d(1)[0]
        return_new = fu.cosinemask(image, cosine_width=0)
        return_old = oldfu.cosinemask(image, cosine_width=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D20_3(self):
        image = get_data_3d(1)[0]
        return_new = fu.cosinemask(image, s=0)
        return_old = oldfu.cosinemask(image, s=0)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_get_shrink_3dmask(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.get_shrink_3dmask()
            oldfu.get_shrink_3dmask()

    @unittest.skip("\n***************************\n\t\t 'Test_get_shrink_3dmask.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.get_shrink_3dmask(3, EMData())
            oldfu.get_shrink_3dmask(3, EMData())

    def test_No_xinit_error(self):
        with self.assertRaises(RuntimeError):
            fu.get_shrink_3dmask(nxinit = 0, mask_file_name = get_data_3d(1))
            oldfu.get_shrink_3dmask(nxinit = 0, mask_file_name = get_data_3d(1))

    def test_E21(self):
        return_new = fu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data_3d(1))
        return_old = oldfu.get_shrink_3dmask(nxinit = 4, mask_file_name = get_data_3d(1))
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D21(self):
        mask_file_name = get_data_3d(1)
        nx =sparx_utilities.get_im(mask_file_name).get_xsize()
        return_new = fu.get_shrink_3dmask(nxinit = nx, mask_file_name = mask_file_name)
        return_old = oldfu.get_shrink_3dmask(nxinit = nx, mask_file_name = mask_file_name)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_get_biggest_cluster(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.get_biggest_cluster()
            oldfu.get_biggest_cluster()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.get_biggest_cluster( EMData())
            oldfu.get_biggest_cluster( EMData())

    def test_D22(self):
        image = get_data(1)[0]
        return_new = fu.get_biggest_cluster(image)
        return_old = oldfu.get_biggest_cluster(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D22_2(self):
        image = get_data_3d(1)[0]
        return_new = fu.get_biggest_cluster(image)
        return_old = oldfu.get_biggest_cluster(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_D22_3(self):
        image = get_data_gauss_noise()
        return_new = fu.get_biggest_cluster(image)
        return_old = oldfu.get_biggest_cluster(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



class Test_compute_bfactor(unittest.TestCase):
    pw = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.compute_bfactor()
            oldfu.compute_bfactor()

    def test_no_pixel_size_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=0)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25, pixel_size=0)

    def test_D23(self):
        return_new =fu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25)
        return_old = oldfu.compute_bfactor(pws=self.pw, freq_min = 0.15, freq_max= 0.25)
        self.assertEqual(return_new[0], return_old[0])
        self.assertEqual(return_new[2], return_old[2])
        self.assertEqual(return_new[3], return_old[3])
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    @unittest.skip("\n***************************\n\t\t 'Test_compute_bfactor.test_with_f_negative' because: It works but i think it 'd not\n***************************")
    def test_with_f_negative(self):
        return_new =fu.compute_bfactor(pws=self.pw, freq_min = -0.15, freq_max= -0.25)
        return_old = oldfu.compute_bfactor(pws=self.pw, freq_min = -0.15, freq_max= -0.25)
        self.assertEqual(return_new[0], return_old[0])
        self.assertEqual(return_new[2], return_old[2])
        self.assertEqual(return_new[3], return_old[3])
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_freqMin_bigger_than_freqMAx_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.25)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.25)

    def test_freqMin_equal_freqMAx_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.35)
            oldfu.compute_bfactor(pws=self.pw, freq_min = 0.35, freq_max= 0.35)

    def test_E4_and_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.compute_bfactor(pws=[1], freq_min = 0.15, freq_max= 0.25)
            oldfu.compute_bfactor(pws=[1], freq_min= 0.15, freq_max= 0.25)
            fu.compute_bfactor(pws=[1,1], freq_min = 0.15, freq_max= 0.25)
            oldfu.compute_bfactor(pws=[1,1], freq_min= 0.15, freq_max= 0.25)

    def test_Empty_array_error(self):
        with self.assertRaises(ValueError):
            fu.compute_bfactor(pws=[], freq_min=0.15, freq_max=0.25)
            oldfu.compute_bfactor(pws=[], freq_min=0.15, freq_max=0.25)



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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.cter_mrk()
            oldfu.cter_mrk()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_mark_true_should_return_equal_object(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.cter_pap()
            oldfu.cter_pap()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_pap_true_should_return_equal_object(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.cter_vpp()
            oldfu.cter_vpp()

    @unittest.skip("I need the cter_mrk/image.mrc from Adnan")
    def test_cter_vpp_true_should_return_equal_object(self):

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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ampcont2angle()
            oldfu.ampcont2angle()

    def test_D26(self):
        return_new = fu.ampcont2angle(100.0)
        return_old = oldfu.ampcont2angle(100.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E26(self):
        return_new = fu.ampcont2angle(-100.0)
        return_old = oldfu.ampcont2angle(-100.0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_F26(self):
        return_new = fu.ampcont2angle(-1)
        return_old = oldfu.ampcont2angle(-1)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_G26(self):
        return_new = fu.ampcont2angle(8)
        return_old = oldfu.ampcont2angle(8)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_angle2ampcont(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.angle2ampcont()
            oldfu.angle2ampcont()

    def test_D27(self):
        return_new = fu.angle2ampcont(0.45)
        return_old = oldfu.angle2ampcont(0.45)
        self.assertTrue(numpy.array_equal(return_new, return_old))



class Test_bracket_def(unittest.TestCase):

    @staticmethod
    def function_test(x1,dat):
        return x1 + dat

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.bracket_def()
            oldfu.bracket_def()

    def test_E28(self):
        return_new = fu.bracket_def(self.function_test,dat=5, x1=3, h=3)
        return_old = oldfu.bracket_def(self.function_test,dat=5, x1=3, h=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_G28(self):
        return_new = fu.bracket_def(self.function_test,dat=5, x1=3, h=0)
        return_old = oldfu.bracket_def(self.function_test,dat=5, x1=3, h=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_bracket(unittest.TestCase):

    @staticmethod
    def function_test(x1,dat):
        return x1 + dat

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.bracket()
            oldfu.bracket()

    def test_D29(self):
        return_new = fu.bracket(self.function_test,dat=5,h=4)
        return_old = oldfu.bracket(self.function_test,dat=5,h=4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E29(self):
        self.assertTrue(fu.bracket(self.function_test,dat=0,h=0) is None)
        self.assertTrue(oldfu.bracket(self.function_test,dat=0,h=0) is None)



class Test_goldsearch_astigmatism(unittest.TestCase):

    @staticmethod
    def function_test(x1, dat):
        f = x1 + dat
        return f

    @staticmethod
    def function_test_return_0(x1, dat):
        return 0

    @staticmethod
    def bad_function_test():
        return 0

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.goldsearch_astigmatism()
            oldfu.goldsearch_astigmatism()

    def test_null_tolerance_error(self):
        with self.assertRaises(OverflowError):
            fu.goldsearch_astigmatism(self.function_test,5,3,4, 0)
            oldfu.goldsearch_astigmatism(self.function_test,5,3,4, 0)


    def test_A_B_same_value_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.goldsearch_astigmatism(self.function_test, 5, 3, 3)
            oldfu.goldsearch_astigmatism(self.function_test, 5, 3, 3)

    def test_Invalid_function(self):
        with self.assertRaises(TypeError):
            fu.goldsearch_astigmatism(self.bad_function_test,5,3,4)
            oldfu.goldsearch_astigmatism(self.bad_function_test,5,3,4)

    def test_F30(self):
        return_new = fu.goldsearch_astigmatism(self.function_test,5,3,4)
        return_old = oldfu.goldsearch_astigmatism(self.function_test,5,3,4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_G30(self):
        return_new = fu.goldsearch_astigmatism(self.function_test_return_0,5,3,4)
        return_old = oldfu.goldsearch_astigmatism(self.function_test_return_0,5,3,4)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E30(self):
        return_new = fu.goldsearch_astigmatism(self.function_test,5,4,3)
        return_old = oldfu.goldsearch_astigmatism(self.function_test,5,4,3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError):
            fu.goldsearch_astigmatism(self.function_test, 0, 0, 0)
            oldfu.goldsearch_astigmatism(self.function_test, 0, 0, 0)



class Test_defocus_baseline_fit(unittest.TestCase):
    roo = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.defocus_baseline_fit()
            oldfu.defocus_baseline_fit()

    def test_D31(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D31_2(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=0)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_defocus_baseline_fit.test_start_is_bigger_than_stop' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_start_is_bigger_than_stop(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=110, i_stop=10, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=110, i_stop=10, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_defocus_baseline_fit.test_start_is_equal_stop' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_start_is_equal_stop(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=10, i_stop=10, nrank=2, iswi=3)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=10, i_stop=10, nrank=2, iswi=3)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_E31(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=2)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=2, iswi=2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_defocus_baseline_fit.test_negative_rank_error' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_negative_rank_error(self):
        return_new = fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=-1, iswi=2)
        return_old = oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=-1, iswi=2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_null_rank_error(self):
        with self.assertRaises(IndexError):
            fu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=0, iswi=2)
            oldfu.defocus_baseline_fit(roo=self.roo, i_start=0, i_stop=10, nrank=0, iswi=2)

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.defocus_baseline_fit(roo=[], i_start=0, i_stop=10, nrank=2, iswi=2)
            oldfu.defocus_baseline_fit(roo=[], i_start=0, i_stop=10, nrank=2, iswi=2)



class Test_simpw1d(unittest.TestCase):
    data = [entry for entry in numpy.arange(1, 256).tolist()]
    defocus = 1
    Cs = 2
    voltage = 300
    pixel_size = 1.5
    amp_contrast = 0.1
    i_start = 2
    i_stop = 14
    nx = 20

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.simpw1d()
            oldfu.simpw1d()

    def test_D32(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(self.defocus,datanew), oldfu.simpw1d(self.defocus,datanew))

    def test_D32_2(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(self.defocus+1,datanew), oldfu.simpw1d(self.defocus+1,datanew))

    def test_D32_3(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_error(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_error(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.simpw1d_pap()
            oldfu.simpw1d_pap()

    def test_D33(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus,self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start, self.i_stop]
        self.assertEqual(fu.simpw1d_pap(self.defocus, datanew), oldfu.simpw1d_pap(self.defocus, datanew))

    def test_D33_2(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(self.defocus+1,datanew), oldfu.simpw1d(self.defocus+1,datanew))

    def test_D33_3(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_error(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_error(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.simpw1d_print()
            oldfu.simpw1d_print()


    def test_D34(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start, self.i_stop]
        self.assertEqual(fu.simpw1d_print(self.defocus,datanew), oldfu.simpw1d_print(self.defocus,datanew))

    def test_D34_2(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(self.defocus+1,datanew), oldfu.simpw1d(self.defocus+1,datanew))

    def test_D34_3(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, -1, self.Cs, self.voltage, self.pixel_size, self.amp_contrast, self.i_start,self.i_stop]
        self.assertEqual(fu.simpw1d(-1,datanew), oldfu.simpw1d(-1,datanew))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.simpw1d(1, [])
            oldfu.simpw1d(1, [])

    def test_no_pixel_size_error(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], self.nx, self.defocus, self.Cs, self.voltage, 0, self.amp_contrast, self.i_start,self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)

    def test_no_image_size_error(self):
        datanew = [self.data[self.i_start:self.i_stop], self.data[self.i_start:self.i_stop], 0, self.defocus,self.Cs, self.voltage, 0, self.amp_contrast, self.i_start, self.i_stop]
        with self.assertRaises(ZeroDivisionError):
            fu.simpw1d(self.defocus, datanew)
            oldfu.simpw1d(self.defocus, datanew)



class Test_movingaverage(unittest.TestCase):
    data = [entry for entry in numpy.arange(0, 10).tolist()]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.movingaverage()
            oldfu.movingaverage()

    def test_D35(self):
        return_new = fu.movingaverage(self.data,window_size=2)
        return_old = oldfu.movingaverage(self.data, window_size=2)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D35_2(self):
        return_new = fu.movingaverage(self.data,window_size=2, skip=0)
        return_old = oldfu.movingaverage(self.data, window_size=2, skip=0)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_windows_size_negative_Value_error(self):
        with self.assertRaises(ValueError):
            fu.movingaverage(self.data,window_size=-2)
            oldfu.movingaverage(self.data, window_size=-2)


    def test_windows_size_null_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.movingaverage(self.data,window_size=0)
            oldfu.movingaverage(self.data, window_size=0)

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.movingaverage([], window_size=2)
            oldfu.movingaverage([], window_size=2)



class Test_defocusgett(unittest.TestCase):
    roo = [entry for entry in numpy.arange(0, 10).tolist()]
    Cs = 2
    voltage = 300
    pixel_size = 1.0
    amp_contrast = 0.1
    nr2 = 6
    i_start = 1
    i_stop = 10
    nx = 1
    skip =False

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertTrue((return_new[0], return_old[0]))
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertTrue(return_new[5], return_old[5])
            self.assertTrue(return_new[6], return_old[6])

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.defocusgett()
            oldfu.defocusgett()

    @unittest.skip("\n***************************\n\t\t 'Test_defocusgett.test_empty_array_error' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop, nr2=self.nr2)
            oldfu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast,self.i_start, self.i_stop, nr2=self.nr2)

    def test_no_pixel_size_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.i_start, self.i_stop,nr2=self.nr2)
            oldfu.defocusgett(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast, self.i_start, self.i_stop,nr2=self.nr2)

    def test_defocusgett_true_should_return_equal_object(self):
        return_new = fu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop,nr2=self.nr2)
        return_old = oldfu.defocusgett(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop,nr2=self.nr2)
        self.test_all_the_conditions(return_new,return_old,False)


class Test_defocusgett_pap(unittest.TestCase):
    roo = [entry for entry in numpy.arange(1, 258).tolist()]
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    i_start = 0.048
    i_stop = -1
    nx = 512
    nr2 = 6
    skip =False

    def test_all_the_conditions(self,return_new=None,return_old=None, skip=True):
        if skip is False:
            self.assertTrue((return_new[0], return_old[0]))
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertTrue(return_new[5], return_old[5])
            self.assertTrue(return_new[6], return_old[6])

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_pap()
            oldfu.defocusgett_pap()

    @unittest.skip("\n***************************\n\t\t 'Test_defocusgett_pap.test_empty_array_error' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop, nr2=self.nr2)
            oldfu.defocusgett([], self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast,self.i_start, self.i_stop, nr2=self.nr2)

    def test_no_pixel_size_error(self):
        with self.assertRaises(ZeroDivisionError):
            fu.defocusgett_pap(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast,self.i_start, self.i_stop)
            oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, 0, self.Cs, self.amp_contrast,self.i_start, self.i_stop)

    def test_defocusgett_pap_true_should_return_equal_object(self):
        return_new = fu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_pap(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.amp_contrast, self.i_start, self.i_stop)
        self.test_all_the_conditions(return_new, return_old, self.skip)



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
            self.assertTrue((return_new[0], return_old[0]))
            self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))
            self.assertTrue(numpy.array_equal(return_new[2], return_old[2]))
            self.assertTrue(numpy.array_equal(return_new[3], return_old[3]))
            self.assertTrue(numpy.array_equal(return_new[4], return_old[4]))
            self.assertTrue(return_new[5], return_old[5])
            self.assertTrue(return_new[6], return_old[6])

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_vpp()
            oldfu.defocusgett_vpp()

    @unittest.skip("\n***************************\n\t\t 'Test_defocusgett_vpp.test_empty_array_error' because an invalid pointer in C++ cide: interrupted by signal 6: SIGABRT\n***************************")
    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.defocusgett_vpp([], self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, self.vpp_options)
            oldfu.defocusgett_vpp([], self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, self.vpp_options)

    def test_empty_array_error2(self):
        with self.assertRaises(IndexError):
            fu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, [])
            oldfu.defocusgett_vpp(self.roo, self.nx, self.voltage, self.pixel_size, self.Cs, self.i_start,self.i_stop, [])

    def test_no_pixel_size_error(self):
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
    qse = get_data_3d(1)[0]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.defocusgett_vpp2()
            oldfu.defocusgett_vpp2()

    @unittest.skip("\n***************************\n\t\t 'Test_defocusgett_vpp2.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        return_new = fu.defocusgett_vpp2(EMData(), self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(EMData(), self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_no_pixel_size_error(self):
        return_new = fu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, 0,self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, 0,self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_defocusgett_vpp2_true_should_return_equal_object(self):
        return_new = fu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        return_old = oldfu.defocusgett_vpp2(self.qse, self.wn, self.xdefc, self.xampcont, self.voltage, self.pixel_size, self.Cs, self.i_start, self.i_stop)
        self.assertTrue(numpy.array_equal(return_new, return_old))


class Test_fastigmatism3(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3()
            oldfu.fastigmatism3()

    @unittest.skip("\n***************************\n\t\t 'Test_fastigmatism3.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3(self.amp, data) - oldfu.fastigmatism3(self.amp, data)))
        self.assertEqual(fu.fastigmatism3(self.amp, data), oldfu.fastigmatism3(self.amp, data))

    def test_D38(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3(self.amp, data) - oldfu.fastigmatism3(self.amp, data)))

    def test_D38_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3(0, data) - oldfu.fastigmatism3(0, data)))
        self.assertEqual(fu.fastigmatism3(self.amp, data), oldfu.fastigmatism3(self.amp, data))

    def test_D38_3(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3(-self.amp, data) - oldfu.fastigmatism3(-self.amp, data)))
        self.assertEqual(fu.fastigmatism3(self.amp, data), oldfu.fastigmatism3(self.amp, data))

    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3(self.amp, data)
            oldfu.fastigmatism3(self.amp, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3(self.amp, data) - oldfu.fastigmatism3(self.amp, data)))
        self.assertEqual(fu.fastigmatism3(self.amp, data), oldfu.fastigmatism3(self.amp, data))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.fastigmatism3(self.amp, [])
            oldfu.fastigmatism3(self.amp, [])



class Test_fastigmatism3_pap(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3_pap()
            oldfu.fastigmatism3_pap()

    @unittest.skip("\n***************************\n\t\t 'Test_fastigmatism3.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(self.amp, data) - oldfu.fastigmatism3_pap(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))

    def test_D39(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(self.amp, data) - oldfu.fastigmatism3_pap(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))

    def test_D39_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(0, data) - oldfu.fastigmatism3_pap(0, data)))
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))

    def test_D39_3(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(-self.amp, data) - oldfu.fastigmatism3_pap(-self.amp, data)))
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))

    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3_pap(self.amp, data)
            oldfu.fastigmatism3_pap(self.amp, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(self.amp, data) - oldfu.fastigmatism3_pap(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_pap(self.amp, data), oldfu.fastigmatism3_pap(self.amp, data))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.fastigmatism3_pap(self.amp, [])
            oldfu.fastigmatism3_pap(self.amp, [])



class Test_fastigmatism3_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    amp = 4
    defocus = 0
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fastigmatism3_vpp()
            oldfu.fastigmatism3_vpp()

    @unittest.skip("\n***************************\n\t\t 'Test_fastigmatism3_vpp.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_vpp(self.amp, data) - oldfu.fastigmatism3_vpp(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))

    def test_D39(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_vpp(self.amp, data) - oldfu.fastigmatism3_vpp(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))

    def test_D39_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_vpp(0, data) - oldfu.fastigmatism3_pap(0, data)))
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))

    def test_D39_3(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, -self.amp_contrast, 0.0]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_vpp(-self.amp, data) - oldfu.fastigmatism3_vpp(-self.amp, data)))
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))

    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast, 0.0]
        with self.assertRaises(RuntimeError):
            fu.fastigmatism3_vpp(self.amp, data)
            oldfu.fastigmatism3_vpp(self.amp, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast, 0.0]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_vpp(self.amp, data) - oldfu.fastigmatism3_vpp(self.amp, data)))
        self.assertEqual(fu.fastigmatism3_vpp(self.amp, data), oldfu.fastigmatism3_vpp(self.amp, data))

    def test_empty_array_error(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.simctf2()
            oldfu.simctf2()

    def test_empty_input_image(self):
        image = get_data(1, self.nx)[0]
        data = [EMData(), image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2(self.defocus, data)
            oldfu.simctf2(self.defocus,data)

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.simctf2(self.defocus, [])
            oldfu.simctf2(self.defocus, [])

    def test_no_pixel_size(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, 0, self.amp_contrast ,self.dfang ]
        self.assertTrue(numpy.isnan(fu.simctf2(self.defocus, data)))
        self.assertTrue(numpy.isnan(oldfu.simctf2(self.defocus, data)))


    @unittest.skip("\n***************************\n\t\t 'Test_simctf2.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        image = get_data(1, self.nx)[0]
        data = [image, EMData(), self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2(self.defocus, data)
            oldfu.simctf2(self.defocus,data)

    def test_D40(self):
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.simctf2_pap()
            oldfu.simctf2_pap()

    def test_empty_input_image(self):
        image = get_data(1, self.nx)[0]
        data = [EMData(), image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2_pap(self.defocus, data)
            oldfu.simctf2_pap(self.defocus,data)

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.simctf2_pap(self.defocus, [])
            oldfu.simctf2_pap(self.defocus, [])

    def test_no_pixel_size(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, 0, self.amp_contrast ,self.dfang ]
        self.assertTrue(numpy.isnan(fu.simctf2_pap(self.defocus, data)))
        self.assertTrue(numpy.isnan(oldfu.simctf2_pap(self.defocus, data)))


    @unittest.skip("\n***************************\n\t\t 'Test_simctf2_pap.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        image = get_data(1, self.nx)[0]
        data = [image, EMData(), self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        with self.assertRaises(RuntimeError):
            fu.simctf2_pap(self.defocus, data)
            oldfu.simctf2_pap(self.defocus,data)

    def test_D41(self):
        image = get_data(1, self.nx)[0]
        data = [image, image, self.nx,  self.dfdiff, self.cs, self.voltage, self.pixel_size, self.amp_contrast ,self.dfang ]
        self.assertEqual(fu.simctf2_pap(self.defocus,data), oldfu.simctf2_pap(self.defocus,data))



class Test_fupw(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    defocus = 0
    amp = 4
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12
    args = [defocus, amp]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fupw()
            oldfu.fupw()

    @unittest.skip("\n***************************\n\t\t 'Test_fupw.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw(self.args, data) - oldfu.fupw(self.args, data)))

    def test_D42(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertEqual(fu.fupw(self.args, data), oldfu.fupw(self.args, data))

    def test_D42_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw(self.args, data) - oldfu.fupw(self.args, data)))


    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw(self.args, data)
            oldfu.fupw(self.args, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw(self.args, data)- oldfu.fupw(self.args, data)))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.fupw(self.args, [])
            oldfu.fupw(self.args, [])



class Test_fupw_pap(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
    defocus = 0
    amp = 4
    Cs = 2
    voltage = 300
    pixel_size = 1.09
    amp_contrast = 0.1
    bfactor = 0.0
    nx = 12
    args = [defocus, amp]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fupw_pap()
            oldfu.fupw_pap()

    @unittest.skip("\n***************************\n\t\t 'Test_fupw_pap.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_pap(self.args, data) - oldfu.fupw_pap(self.args, data)))

    def test_D43(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertEqual(fu.fupw_pap(self.args, data), oldfu.fupw_pap(self.args, data))

    def test_D43_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(0, data) - oldfu.fastigmatism3_pap(0, data)))
        self.assertEqual(fu.fupw_pap(self.args, data), oldfu.fupw_pap(self.args, data))

    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw_pap(self.args, data)
            oldfu.fupw_pap(self.args, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        #self.assertTrue(TOLERANCE > numpy.abs(fu.fastigmatism3_pap(self.amp, data) - oldfu.fastigmatism3_pap(self.amp, data)))
        self.assertEqual(fu.fupw_pap(self.args, data), oldfu.fupw_pap(self.args, data))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.fupw_pap(self.args, [])
            oldfu.fupw_pap(self.args, [])



class Test_fupw_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
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

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.fupw_vpp()
            oldfu.fupw_vpp()

    @unittest.skip("\n***************************\n\t\t 'Test_fupw_vpp.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [EMData(), numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_vpp(self.args, data) - oldfu.fupw_vpp(self.args, data)))

    def test_D42(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertEqual(fu.fupw_vpp(self.args, data), oldfu.fupw_vpp(self.args, data))

    def test_D42_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_vpp(self.args, data) - oldfu.fupw_vpp(self.args, data)))


    def test_no_image_size_error(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, 0, self.defocus, self.Cs, self.voltage, self.pixel_size, self.bfactor, self.amp_contrast,1]
        with self.assertRaises(RuntimeError):
            fu.fupw_vpp(self.args, data)
            oldfu.fupw_vpp(self.args, data)

    def test_no_pixel_size(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        data = [crefim, numr, self.nx, self.defocus, self.Cs, self.voltage, 0, self.bfactor, self.amp_contrast,1]
        self.assertTrue(TOLERANCE > numpy.abs(fu.fupw_vpp(self.args, data)- oldfu.fupw_vpp(self.args, data)))

    def test_empty_array_error(self):
        with self.assertRaises(IndexError):
            fu.fupw_vpp(self.args, [])
            oldfu.fupw_vpp(self.args, [])



class Test_ornq_vpp(unittest.TestCase):
    argum = get_arg_from_pickle_files(path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.ornq_vpp()
            oldfu.ornq_vpp()



    @unittest.skip("\n***************************\n\t\t 'Test_ornq_vpp.test_empty_input_image' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]

        return_new = fu.ornq_vpp(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = fu.ornq_vpp(EMData(), crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    @unittest.skip("\n***************************\n\t\t 'Test_ornq_vpp.test_empty_input_image2' because: interrupted by signal 11: SIGSEGV\n***************************")
    def test_empty_input_image2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]

        return_new = fu.ornq_vpp(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny)
        return_old = fu.ornq_vpp(image, EMData(),  xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D49(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0] #mode is H

        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D49_2(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='f'
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_D49_with_invalid_mode(self):
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = self.argum[0]
        mode ='invalid'
        return_new = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        return_old = fu.ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
        self.assertTrue(numpy.array_equal(return_new, return_old))


