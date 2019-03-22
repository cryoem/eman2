from __future__ import print_function
from __future__ import division


import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

from ..libpy_py3 import sphire_filter as fu
from .sparx_lib import sparx_filter as oldfu

from ..libpy import sparx_utilities as ut

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


def get_data_3d(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
        data_list.append(a)

    return data_list

def get_data_gauss_noise():
    dim = 10
    return ut.model_gauss_noise(0.25 , dim,dim,dim)

from EMAN2_cppwrap import EMData

class Test_filt_ctd(unittest.TestCase):
    image = get_data(1)[0]

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_ctf()
            oldfu.filt_ctf()

    def test_empty_input_image(self):
        image =EMData()
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})
        with self.assertRaises(AssertionError):
            fu.filt_ctf(image,ctf)
            oldfu.filt_ctf(image,ctf)

    def test_E9(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(self.image,ctf)
        return_old = oldfu.filt_ctf(self.image,ctf)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_F9(self):
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        return_new = fu.filt_ctf(self.image,ctf, False)
        return_old = oldfu.filt_ctf(self.image,ctf, False)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_raiseError(self):
        # Since there is no try-except paradigma but only an assert I suppose that the most common error is given a None vuole instead an image
        with self.assertRaises(AttributeError) :
            fu.filt_ctf(None,None)


class Test_filt_tophatl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_tophatl()
            oldfu.filt_tophatl()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_tophatl(EMData(), freq=0.25)
            oldfu.filt_tophatl(EMData(), freq=0.25)

    def test_D1(self):
        image, = get_data(1)
        return_new = fu.filt_tophatl(image, freq= 0.25)
        return_old = oldfu.filt_tophatl(image,freq= 0.25)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_tophatb(unittest.TestCase):

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_tophatb(EMData(), freql=0.25, freqh=0.35)
            oldfu.filt_tophatb(EMData(), freql=0.25, freqh=0.35)

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_tophatb()
            oldfu.filt_tophatb()

    def test_D2(self):
        image, = get_data(1)
        return_new = fu.filt_tophatb(image,freql= 0.25, freqh= 0.35)
        return_old = oldfu.filt_tophatb(image,freql= 0.25, freqh= 0.35)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_gaussl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_gaussl()
            oldfu.filt_gaussl()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_gaussl(EMData(), sigma=0.23)
            oldfu.filt_gaussl(EMData(), sigma=0.23)

    def test_D3(self):
        image, = get_data(1)
        return_new = fu.filt_gaussl(image,sigma= 0.23 )
        return_old = oldfu.filt_gaussl(image,sigma= 0.23)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_gaussinv(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_gaussinv()
            oldfu.filt_gaussinv()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_gaussinv(EMData(), sigma=0.23)
            oldfu.filt_gaussinv(EMData(), sigma=0.23)

    def test_D4(self):
        image, = get_data(1)
        return_new = fu.filt_gaussinv(image,sigma= 0.23)
        return_old = oldfu.filt_gaussinv(image,sigma= 0.23)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_gaussh(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_gaussh()
            oldfu.filt_gaussh()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_gaussh(EMData(), sigma=0.23)
            oldfu.filt_gaussh(EMData(), sigma=0.23)

    def test_filt_gaussh_true_should_return_equal_object(self):
        image, = get_data(1)
        return_new = fu.filt_gaussh(image,sigma= 0.23)
        return_old = oldfu.filt_gaussh(image,sigma= 0.23)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_btwl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_btwl()
            oldfu.filt_btwl()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_btwl(EMData(), freql=0.25, freqh=0.35)
            oldfu.filt_btwl(EMData(), freql=0.25, freqh=0.35)

    def test_D6(self):
        image, = get_data(1)
        return_new = fu.filt_btwl(image,freql= 0.25, freqh= 0.35)
        return_old = oldfu.filt_btwl(image,freql= 0.25, freqh= 0.35)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_tanl(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_tanl()
            oldfu.filt_tanl()

    def test_empty_input_image(self):
        with self.assertRaises(RuntimeError):
            fu.filt_tanl(EMData(), freq=0.25, fall_off=0.35)
            oldfu.filt_tanl(EMData(), freq=0.25, fall_off=0.35)

    def test_D7(self):
        image, = get_data(1)
        return_new = fu.filt_tanl(image,freq= 0.25, fall_off=0.35)
        return_old = oldfu.filt_tanl(image,freq= 0.25, fall_off=0.35)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


class Test_filt_table(unittest.TestCase):

    def test_wrong_number_params(self):
        with self.assertRaises(TypeError):
            fu.filt_table()
            oldfu.filt_table()

    def test_empty_input_image(self):
        table = [entry for entry in numpy.linspace(0, 0.5).tolist()]
        with self.assertRaises(RuntimeError):
            fu.filt_table(EMData(), table)
            oldfu.filt_table(EMData(), table)

    def test_filt_table_true_should_return_equal_object(self):
        image, = get_data(1)
        table =[entry for entry in numpy.linspace(0, 0.5).tolist()]

        return_new = fu.filt_table(image,table)
        return_old = oldfu.filt_table(image,table)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))














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
        vols.append(get_data_gauss_noise())
        vols.append(get_data_gauss_noise())
        vols.append(get_data_gauss_noise())

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
        vols.append(get_data_gauss_noise())
        vols.append(get_data_gauss_noise())
        vols.append(get_data_gauss_noise())
        ui =  get_data_gauss_noise()     # or use ut.model_blank(1,1,1)
        vi =  get_data_gauss_noise()
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

