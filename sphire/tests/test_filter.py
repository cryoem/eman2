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

