from __future__ import print_function
from __future__ import division
from past.utils import old_div

import nose.tools as nt

import EMAN2_cppwrap as e2cpp
import numpy
import copy
import math
from ..libpy_py3 import sphire_fundamentals as fu
# from ..libpy import sparx_fundamentals as fu
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
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

        oldfubrackets =  oldfu.symclass('c1').brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

        self.assertEqual(fubrackets , oldfubrackets)


    def test_c5_should_return_correct_brackets(self):
        nsym = 1
        fubrackets =  fu.symclass('c5').brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

        oldfubrackets =  oldfu.symclass('c5').brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

        self.assertEqual(fubrackets , oldfubrackets)

    def test_c1_should_return_correct_symangles(self):
        nsym = 1
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i * 360. / nsym])

        self.assertEqual(fu.symclass('c1').symangles, symangles)
        self.assertEqual(oldfu.symclass('c1').symangles, symangles)

    def test_c5_should_return_correct_symangles(self):
        nsym = 5
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i * 360. / nsym])

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
        assert results == expected_results

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
        assert results == expected_results









if __name__ == '__main__':
    unittest.main()
