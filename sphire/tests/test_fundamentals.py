from __future__ import print_function
from __future__ import division
from past.utils import old_div

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
