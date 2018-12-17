from __future__ import print_function
from __future__ import division

import numpy
import copy
import math
import EMAN2_cppwrap as e2cpp

from ..libpy import sparx_morphology as fu
from .sparx_lib import sparx_morphology as oldfu

from ..libpy import sparx_filter as ft
from ..libpy import sparx_utilities as ut

import pickle
import os
import shutil
ABSOLUTE_PATH =  os.path.dirname(os.path.realpath(__file__))


import unittest

def get_data(num, dim=10):
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


def get_data_3d(num, dim=10):
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

    def test_binarize_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.binarize(image)
        return_old = oldfu.binarize(image)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_collapse_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.collapse(image)
        return_old = oldfu.collapse(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_dilation_true_should_return_equal_object(self):
        image = ut.model_blank(10,10)

        return_new = fu.dilation(image)
        return_old = oldfu.dilation(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_erosion_true_should_return_equal_object(self):
        image = ut.model_blank(10,10)

        return_new = fu.erosion(image)
        return_old = oldfu.erosion(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_power_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.power(image)
        return_old = oldfu.power(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_square_root_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.square_root(image)
        return_old = oldfu.square_root(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_square_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.square(image)
        return_old = oldfu.square(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_threshold_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.threshold(image)
        return_old = oldfu.threshold(image)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_threshold_outside_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.threshold_outside(image, 2, 10)
        return_old = oldfu.threshold_outside(image, 2, 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_notzero_true_should_return_equal_object(self):
        image, = get_data(1)

        return_new = fu.threshold_outside(image, 2, 10)
        return_old = oldfu.threshold_outside(image, 2, 10)
        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_rotavg_ctf_true_should_return_equal_object(self):
        image, = get_data(1)
        defocus = 1
        Cs = 2
        voltage = 300
        Pixel_size = 1.5

        return_new = fu.rotavg_ctf(image, defocus, Cs, voltage, Pixel_size)
        return_old = oldfu.rotavg_ctf(image, defocus, Cs, voltage, Pixel_size)

        self.assertTrue(return_new , return_old)

    def test_ctf_1d_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs =  2
        voltage = 300
        pixel_size = 1.5
        bfactor =  0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf_1d(nx, ctf)
        return_old= fu.ctf_1d(nx, ctf)

        self.assertTrue(return_new, return_old)


    def test_ctf_2_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs =  2
        voltage = 300
        pixel_size = 1.5
        bfactor =  0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf_2(nx, ctf)
        return_old= fu.ctf_2(nx, ctf)

        self.assertTrue(return_new, return_old)

    def test_ctf_img_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 1.5
        bfactor = 0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf_img(nx, ctf)
        return_old = fu.ctf_img(nx, ctf)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_ctf_img_real_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 1.5
        bfactor = 0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf_img_real(nx, ctf)
        return_old = fu.ctf_img_real(nx, ctf)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_ctf_rimg_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 1.5
        bfactor = 0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf_rimg(nx, ctf)
        return_old = fu.ctf_rimg(nx, ctf)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_ctf2_rimg_true_should_return_equal_object(self):
        nx = 20
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 1.5
        bfactor = 0
        amp_contrast = 0.1
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})

        return_new = fu.ctf2_rimg(nx, ctf)
        return_old = fu.ctf2_rimg(nx, ctf)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_ctflimit_true_should_return_equal_object(self):
        nx = 30
        defocus = 1
        Cs = 2
        voltage = 300
        Pixel_size = 1.5

        return_new = fu.ctflimit(nx, defocus, Cs, voltage, Pixel_size)
        return_old = oldfu.ctflimit(nx, defocus, Cs, voltage, Pixel_size)

        self.assertTrue(return_new , return_old)


    def test_imf_params_cl1_true_should_return_equal_object(self):
        pw = [entry for entry in numpy.arange(0,10).tolist()]

        return_new = fu.imf_params_cl1(pw)
        return_old = oldfu.imf_params_cl1(pw)

        self.assertTrue(return_new, return_old)

    def test_adaptive_mask_true_should_return_equal_object(self):
        image, = get_data_3d(1)

        return_new = fu.adaptive_mask(image)
        return_old = oldfu.adaptive_mask(image)

        self.assertTrue(return_new, return_old)


    def test_cosinemask_true_should_return_equal_object(self):
        image, = get_data_3d(1)

        return_new = fu.cosinemask(image)
        return_old = oldfu.cosinemask(image)

        self.assertTrue(return_new, return_old)


    def test_shrink_3dmask_true_should_return_equal_object(self):
        nxinit = 4
        mask_file_name = get_data_3d(1)
        return_new = fu.get_shrink_3dmask(nxinit,mask_file_name)
        return_old = oldfu.get_shrink_3dmask(nxinit,mask_file_name)

        self.assertTrue(return_new, return_old)


    def test_get_biggest_cluster_true_should_return_equal_object(self):
        image = ut.model_blank(10, 10)

        return_new = fu.get_biggest_cluster(image)
        return_old = oldfu.get_biggest_cluster(image)

        self.assertTrue(return_new, return_old)


    def test_computer_bfactor_true_should_return_equal_object(self):
        pw = [entry for entry in numpy.arange(0, 10).tolist()]
        freq_min = 0.15
        freq_max = 0.25

        return_new = fu.compute_bfactor(pw, freq_min, freq_max)
        return_old = oldfu.compute_bfactor(pw, freq_min, freq_max)

        self.assertTrue(return_new, return_old)

    def test_cter_mrk_true_should_return_equal_object(self):
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 0.5
        bfactor = 0
        amp_contrast = 0.1
        wn = 32
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})
        image1, = get_data(1, 256)

        micrograph_image = ft.filt_ctf(image1,ctf)
        micrograph_image.write_image('sphire/tests/files/cter_mrk/image.mrc')
        selection_list = 'image.mrc'

        input_image_path = os.path.join(ABSOLUTE_PATH, "files/cter_mrk/image*.mrc")
        output_directory = os.path.join(ABSOLUTE_PATH, "files/cter_mrk/results")
        if os.path.isdir(output_directory):
            shutil.rmtree(output_directory)

        return_new = fu.cter_mrk(input_image_path, output_directory, selection_list= selection_list,\
                                 wn = wn, pixel_size = pixel_size, Cs= cs, voltage = voltage)

        if os.path.isdir(output_directory):
            shutil.rmtree(output_directory)

        return_old = oldfu.cter_mrk(input_image_path, output_directory, selection_list = selection_list,\
                                    wn = wn,  pixel_size=pixel_size, Cs= cs, voltage = voltage)

        self.assertEqual(return_new, return_old)


    def test_cter_pep_true_should_return_equal_object(self):
        defocus = 1
        cs = 2
        voltage = 300
        pixel_size = 1.5
        bfactor = 0
        amp_contrast = 0.1
        wn = 32
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": defocus, "cs": cs, "voltage": voltage, "apix": pixel_size, "bfactor": bfactor,
                       "ampcont": amp_contrast})
        image1, = get_data(1, 256)

        micrograph_image = ft.filt_ctf(image1,ctf)
        micrograph_image.write_image('sphire/tests/files/cter_mrk/image.mrc')
        selection_list = 'image.mrc'

        input_image_path = os.path.join(ABSOLUTE_PATH, "files/cter_mrk/image*.mrc")
        output_directory = os.path.join(ABSOLUTE_PATH, "files/cter_mrk/results")
        if os.path.isdir(output_directory):
            shutil.rmtree(output_directory)

        return_new = fu.cter_pap(input_image_path, output_directory, selection_list = selection_list,\
                                    wn = wn,  pixel_size=pixel_size, Cs= cs, voltage = voltage)

        if os.path.isdir(output_directory):
            shutil.rmtree(output_directory)

        return_old = oldfu.cter_pap(input_image_path, output_directory, selection_list = selection_list,\
                                    wn = wn,  pixel_size=pixel_size, Cs= cs, voltage = voltage)

        self.assertEqual(return_new, return_old)


    def test_ampcont2angle_true_should_return_equal_object(self):
        return_new = fu.ampcont2angle(8)
        return_old = oldfu.ampcont2angle(8)

        self.assertTrue(return_new, return_old)

    def test_angle2ampcont_true_should_return_equal_object(self):
        return_new = fu.angle2ampcont(0.45)
        return_old = oldfu.angle2ampcont(0.45)

        self.assertTrue(return_new, return_old)

    def test_bracket_def_true_should_return_equal_object(self):
        def function_test(x1,dat):
            f = x1 + dat
            return f
        dat = 5
        x1 = 3
        h  = 3

        return_new = fu.bracket_def(function_test,dat, x1, h)
        return_old = oldfu.bracket_def(function_test,dat, x1, h)

        self.assertTrue(return_new, return_old)

    def test_bracket_true_should_return_equal_object(self):
        def function_test(x1,dat):
            f = x1 + dat
            return f
        dat = 5
        x1 = 3
        h  = 3

        return_new = fu.bracket(function_test,dat,h)
        return_old = oldfu.bracket(function_test,dat,h)

        self.assertTrue(return_new, return_old)

    def test_goldsearch_astigmatism_true_should_return_equal_object(self):
        def function_test(x1,dat):
            f = x1 + dat
            return f
        dat = 5
        x1 = 3
        a  = 3
        b = 4

        return_new = fu.goldsearch_astigmatism(function_test,dat,a,b)
        return_old = oldfu.goldsearch_astigmatism(function_test,dat,a,b)

        self.assertTrue(return_new, return_old)

    def test_defocus_baseline_fit_true_should_return_equal_object(self):
        roo = [entry for entry in numpy.arange(0, 10).tolist()]
        i_start = 0
        i_stop = 10
        nrank = 2
        iswi = 3

        return_new = fu.defocus_baseline_fit(roo, i_start, i_stop, nrank, iswi)
        return_old = oldfu.defocus_baseline_fit(roo, i_start, i_stop, nrank, iswi)

        self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_simpw1d_true_should_return_equal_object(self):
        data = [entry for entry in numpy.arange(1, 256).tolist()]
        defocus = 1
        Cs = 2
        voltage = 300
        pixel_size = 1.5
        amp_contrast = 0.1
        i_start = 2
        i_stop = 14
        nx = 20

        datanew = [data[i_start:i_stop], data[i_start:i_stop], nx, defocus, \
                Cs, voltage, pixel_size, amp_contrast, i_start,i_stop]

        return_new = fu.simpw1d(defocus,datanew)
        return_old = oldfu.simpw1d(defocus,datanew)

        self.assertTrue(return_new, return_old)


    def test_simpw1d_pap_true_should_return_equal_object(self):
        data = [entry for entry in numpy.arange(1, 256).tolist()]
        defocus = 1
        Cs = 2
        voltage = 300
        pixel_size = 1.5
        amp_contrast = 0.1
        i_start = 2
        i_stop = 14
        nx = 20

        datanew = [data[i_start:i_stop], data[i_start:i_stop], nx, defocus, \
                Cs, voltage, pixel_size, amp_contrast, i_start,i_stop]

        return_new = fu.simpw1d_pap(defocus,datanew)
        return_old = oldfu.simpw1d_pap(defocus,datanew)

        self.assertTrue(return_new, return_old)


    def test_simpw1d_print_true_should_return_equal_object(self):

        data = [entry for entry in numpy.arange(1, 256).tolist()]
        defocus = 1
        Cs = 2
        voltage = 300
        pixel_size = 1.5
        amp_contrast = 0.1
        i_start = 2
        i_stop = 14
        nx = 20

        datanew = [data[i_start:i_stop], data[i_start:i_stop], nx, defocus, \
                Cs, voltage, pixel_size, amp_contrast, i_start,i_stop]
        return_new = fu.simpw1d_print(defocus,datanew)
        return_old = oldfu.simpw1d_print(defocus,datanew)

        self.assertTrue(return_new, return_old)


    def test_movingaverage_true_should_return_equal_object(self):
        data = [entry for entry in numpy.arange(0, 10).tolist()]
        window_size = 2

        return_new = fu.movingaverage(data,window_size)
        return_old = oldfu.movingaverage(data, window_size)

        self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_defocusgett_true_should_return_equal_object(self):
        roo = [entry for entry in numpy.arange(0, 10).tolist()]
        Cs = 2
        voltage = 300
        pixel_size = 1.0
        amp_contrast = 0.1
        nr2 = 6
        i_start = 1
        i_stop = 10
        nx = 1

        return_new = fu.defocusgett(roo,nx,voltage,pixel_size,Cs,amp_contrast,i_start,i_stop,nr2=nr2)
        return_old = oldfu.defocusgett(roo,nx,voltage,pixel_size,Cs,amp_contrast,i_start,i_stop,nr2=nr2)

        self.assertTrue(return_new, return_old)

    def test_defocusgett_pap_true_should_return_equal_object(self):
            roo = [entry for entry in numpy.arange(1, 258).tolist()]
            Cs = 2
            voltage = 300
            pixel_size = 1.09
            amp_contrast = 0.1
            i_start = 0.048
            i_stop = -1
            nx = 512
            nr2 = 6

            return_new = fu.defocusgett_pap(roo,nx,voltage,pixel_size,Cs,amp_contrast, \
               i_start,i_stop)
            return_old = oldfu.defocusgett_pap(roo,nx,voltage,pixel_size,Cs,amp_contrast, \
               i_start,i_stop)

            self.assertTrue(return_new, return_old)


    # def test_fastigmatism3_true_should_return_equal_object(self):
    #     data = [entry for entry in numpy.arange(0, 10).tolist()]
    #     amp = 4
    #
    #     return_new = fu.fastigmatism3(amp,data)
    #     return_old = oldfu.fastigmatism3(amp,data)
    #
    #     self.assertTrue(return_new, return_old)






if __name__ == '__main__':
    unittest.main()
