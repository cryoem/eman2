from past.utils import old_div

import nose.tools as nt

import EMAN2_cppwrap as e2cpp
import EMAN2_cppwrap
import numpy
import copy
import math

from sparx.libpy import fundamentals as fu
import fundamentals as oldfu



import global_def
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


# class symclass_mod(object):
#     def __init__(self, sym):
#         """
#           sym: cn, dn, oct, tet, icos
#         """
#         pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt, pi
#         # from utilities import get_sym, get_symt
#         pass  # IMPORTIMPORTIMPORT from string import lower
#         self.sym = sym.lower()
#         if (self.sym[0] == "c"):
#             self.nsym = int(self.sym[1:])
#             if (self.nsym < 1):  global_def.ERROR("For Cn symmetry, we need n>0", "symclass", 1)
#             self.brackets = [[old_div(360., self.nsym), 90.0, old_div(360., self.nsym), 90.0],
#                              [old_div(360., self.nsym), 180.0, old_div(360., self.nsym), 180.0]]
#             self.symangles = []
#             for i in range(self.nsym):
#                 self.symangles.append([0.0, 0.0, i * old_div(360., self.nsym)])
#
#         elif (self.sym[0] == "d"):
#             self.nsym = 2 * int(self.sym[1:])
#             if (self.nsym < 1):  global_def.ERROR("For Dn symmetry, we need n>0", "symclass", 1)
#             self.brackets = [[old_div(360., self.nsym), 90.0, old_div(360., self.nsym), 90.0],
#                              [old_div(360., self.nsym) * 2, 90.0, old_div(360., self.nsym) * 2, 90.0]]
#             self.symangles = []
#             for i in range(old_div(self.nsym, 2)):
#                 self.symangles.append([0.0, 0.0, i * old_div(360., self.nsym) * 2])
#             for i in reversed(range(old_div(self.nsym, 2))):
#                 self.symangles.append(
#                     [0.0, 180.0, (i * old_div(360., self.nsym) * 2 + 180.0 * (int(self.sym[1:]) % 2)) % 360.0])
#
#         elif (self.sym[:3] == "oct"):
#             self.nsym = 24
#             ncap = 4
#             cap_sig = old_div(360.0, ncap)
#             alpha = numpy.degrees(math.acos(
#                 old_div(1.0, (numpy.sqrt(3.0) * numpy.tan(
#                     2 * old_div(old_div(numpy.pi, ncap), 2.0))))))  # also platonic_params["alt_max"]
#             theta = numpy.degrees(0.5 * math.acos(old_div(numpy.cos(numpy.radians(cap_sig)),
#                                                           (1.0 - numpy.cos(numpy.radians(
#                                                               cap_sig))))))  # also platonic_params["theta_c_on_two"]
#             self.brackets = [[old_div(180., ncap), theta, cap_sig, alpha], [old_div(360., ncap), theta, cap_sig, alpha]]
#             self.symangles = [[0.0, 0.0, float(i)] for i in range(0, 271, 90)]
#             for i in range(0, 271, 90):
#                 for j in range(0, 271, 90):
#                     self.symangles.append([float(j), 90.0, float(i)])
#             for i in range(0, 271, 90):  self.symangles.append([0.0, 180.0, float(i)])
#
#         elif (self.sym[:3] == "tet"):
#             self.nsym = 12
#             ncap = 3
#             cap_sig = old_div(360.0, ncap)
#             alpha = numpy.degrees(math.acos(old_div(1.0, (numpy.sqrt(3.0) * numpy.tan(
#                 2 * old_div(old_div(numpy.pi, ncap), 2.0))))))  # also platonic_params["alt_max"]
#             theta = numpy.degrees(0.5 * math.acos(old_div(numpy.cos(numpy.radians(cap_sig)),
#                                                           (1.0 - numpy.cos(numpy.radians(
#                                                               cap_sig))))))  # also platonic_params["theta_c_on_two"]
#             self.brackets = [[old_div(360.0, ncap), theta, cap_sig, alpha],
#                              [old_div(360.0, ncap), theta, cap_sig, alpha]]
#             lvl1 = numpy.degrees(math.acos(old_div(-1.0, 3.0)))  # There  are 3 faces at this angle
#             self.symangles = [[0., 0., 0.], [0., 0., 120.], [0., 0., 240.]]
#             for l1 in range(0, 241, 120):
#                 for l2 in range(60, 301, 120):
#                     self.symangles.append([float(l1), lvl1, float(l2)])
#
#             """Multiline Comment4"""
#
#         elif (self.sym[:4] == "icos"):
#             self.nsym = 60
#             ncap = 5
#             cap_sig = old_div(360.0, ncap)
#             alpha = numpy.degrees(math.acos(old_div(1.0, (numpy.sqrt(3.0) * numpy.tan(
#                 2 * old_div(old_div(numpy.pi, ncap), 2.0))))))  # also platonic_params["alt_max"]
#             theta = numpy.degrees(0.5 * math.acos(old_div(numpy.cos(numpy.radians(cap_sig)),
#                                                           (1.0 - numpy.cos(numpy.radians(
#                                                               cap_sig))))))  # also platonic_params["theta_c_on_two"]
#             self.brackets = [[36., theta, cap_sig, alpha], [72., theta, cap_sig, alpha]]
#             lvl1 = numpy.degrees(math.atan(2.0))  # there are 5 pentagons with centers at this height (angle)
#             lvl2 = 180.0 - lvl1  # there are 5 pentagons with centers at this height (angle)
#             self.symangles = [[0.0, 0.0, float(i)] for i in range(0, 288 + 1, 72)]
#             for l1 in range(0, 288 + 1, 72):
#                 for l2 in range(36, 324 + 1, 72):
#                     self.symangles.append([float(l1), lvl1, float(l2)])
#             for l1 in range(36, 324 + 1, 72):
#                 for l2 in range(0, 288 + 1, 72):
#                     self.symangles.append([float(l1), lvl2, float(l2)])
#             for i in range(0, 288 + 1, 72):  self.symangles.append([0.0, 180.0, float(i)])
#
#         else:
#             global_def.ERROR("Unknown symmetry", "symclass", 1)
#
#         #
#         self.transform = []
#         for args in self.symangles:
#             self.transform.append(
#                 EMAN2_cppwrap.Transform({"type": "spider", "phi": args[0], "theta": args[1], "psi": args[2]}))
#         self.symatrix = self.rotmatrix(self.symangles)
#         self.round = 12
#
#
#
#     def symmetry_related(self, angles, return_mirror=0, neighbors = None,  tolistconv=True):
#
#         symatrix = numpy.array(self.symatrix, numpy.float64)
#         if neighbors is None:
#             symatrix = symatrix
#         else:
#             symatrix = symatrix[neighbors]
#         n_neighbors = symatrix.shape[0]
#         mirror_list = []
#         angles = numpy.atleast_2d(numpy.array(angles))
#         if return_mirror == 1:
#             nsym = n_neighbors
#             mult = -1
#             mask = numpy.ones(nsym*angles.shape[0], dtype=numpy.bool)
#             mirror_list.append([mult, mask])
#         elif return_mirror == 0:
#             nsym = n_neighbors
#             mult = 1
#             mask = numpy.ones(nsym*angles.shape[0], dtype=numpy.bool)
#             mirror_list.append([mult, mask])
#         else:
#             nsym = 2 * n_neighbors
#             mult = 1
#             mask = numpy.zeros(nsym*angles.shape[0], dtype=numpy.bool)
#             for i in range(n_neighbors):
#                 mask[i::2*n_neighbors] = True
#             mirror_list.append([mult, mask])
#             mirror_list.append([-mult, ~mask])
#
#         if(self.sym[0] == "c"):
#             inside_values = (
#                 (self.brackets[0][0], 0, 360.0, 0),
#                 (self.brackets[0][0], 180, 360.0, 0),
#             )
#         elif(self.sym[0] == "d"):
#             inside_values = (
#                 (self.brackets[0][0], 0, 360.0, 0),
#                 (self.brackets[0][0], 180, 360.0, 0),
#                 (self.brackets[0][0], 0, 360.0, self.brackets[0][0]),
#                 (self.brackets[0][0], 180.0, 360.0,self.brackets[0][0]),
#                 (numpy.nan, 90.0, 180.0, 0),
#             )
#         elif (self.sym == "tet") :
#             inside_values = (
#                 (self.brackets[0][0], self.brackets[0][1], 180, 0),
#                 (self.brackets[0][0], 180 - self.brackets[0][1], 180, 60),
#                 (self.brackets[0][0], 0, self.brackets[0][0], 0),
#                 (self.brackets[0][0], 180 - self.brackets[0][3], self.brackets[0][0], 0),
#                 (self.brackets[0][0], 180, self.brackets[0][0], 0),
#                 (self.brackets[0][0], self.brackets[0][3], self.brackets[0][0], 60),
#             )
#         elif(self.sym == "oct") :
#             inside_values = (
#                 (self.brackets[0][2], 180, self.brackets[0][2], 0),
#                 (self.brackets[0][2], 0, self.brackets[0][2], 0),
#                 (self.brackets[0][2], 2 * self.brackets[0][1], self.brackets[0][2], 0),
#                 (self.brackets[0][2], 2 * self.brackets[0][1], 180, 45),
#                 (self.brackets[0][2], 3 * self.brackets[0][1], 180, 0),
#                 (self.brackets[0][2], self.brackets[0][1], 180, 0),
#                 (self.brackets[0][2], self.brackets[0][3], 120, 45),
#                 (self.brackets[0][2], 180 - self.brackets[0][3], 120, 45),
#             )
#         elif(self.sym == "icos"):
#             inside_values = (
#                 (self.brackets[0][2], 180, self.brackets[0][2], 0),
#                 (self.brackets[0][2], 0, self.brackets[0][2], 0),
#                 (self.brackets[0][2], 2 * self.brackets[0][1], self.brackets[0][2], 0),
#                 (self.brackets[0][2], 180 - 2 * self.brackets[0][1], self.brackets[0][2], self.brackets[0][0]),
#                 (self.brackets[0][2], self.brackets[0][3], 60, self.brackets[0][0]),
#                 (self.brackets[0][2], self.brackets[0][3]+2*self.brackets[0][1], 120, 0),
#                 (self.brackets[0][2], 180 - self.brackets[0][3] - 2 * self.brackets[0][1], 120, self.brackets[0][0]),
#                 (self.brackets[0][2], 180 - self.brackets[0][3], 120, 0),
#                 (self.brackets[0][2], self.brackets[0][1], 180, 0),
#                 (self.brackets[0][2], 90 - self.brackets[0][1], 180, self.brackets[0][0]),
#                 (self.brackets[0][0], 90, 180, self.brackets[0][0]/2.0),
#                 (self.brackets[0][2], 180 - self.brackets[0][1], 180, self.brackets[0][0]),
#                 (self.brackets[0][2], 90 + self.brackets[0][1], 180, 0),
#             )
#         else :
#             raise NameError("Symmetry unknown")
#
#         sang_new_raw = numpy.atleast_2d(numpy.array(angles, numpy.float64)).repeat(nsym, axis=0)
#         final_masks = []
#         for multiplier, sang_mask in mirror_list:
#
#             if return_mirror not in (0, 1) and self.sym[0] == 'd' and multiplier == -1:
#                 theta_0_or_180 = (sang_new_raw[:, 1] == 0) | (sang_new_raw[:, 1] == 180)
#                 sang_mask[theta_0_or_180] = False
#
#             if return_mirror not in (0, 1) and self.sym[0] == 'i' and multiplier == -1:
#                 theta_0 = (sang_new_raw[:, 1] == 0)
#                 sang_mask[theta_0] = False
#
#
#             sang_mod = sang_new_raw[sang_mask]
#
#             matrices = self.rotmatrix(sang_mod, tolistconv=False)
#
#             matrices_mod = numpy.sum(
#                 numpy.transpose(matrices, (0, 2, 1)).reshape(
#                     matrices.shape[0],
#                     matrices.shape[1],
#                     matrices.shape[2],
#                     1
#                 ) *
#                 numpy.tile(
#                     multiplier * symatrix,
#                     (sang_mod.shape[0] // n_neighbors, 1, 1)).reshape(
#                     matrices.shape[0], matrices.shape[1], 1, matrices.shape[2], ), -3)
#
#             sang_new = self.recmat(matrices_mod, tolistconv=False)
#             theta_0_or_180 = (sang_new[:,1] == 0) | (sang_new[:,1] == 180)
#             if return_mirror not in (0, 1) and self.sym[0] != 'c' and multiplier == -1:
#                 sang_new[~theta_0_or_180, 2] += 180
#                 sang_new[~theta_0_or_180, 2] %= 360
#
#             masks_good = []
#             masks_bad = []
#
#             for phi, theta, psi, offset in inside_values:
#
#                 if not numpy.isnan(phi):
#                     phi_0_180 = numpy.round(sang_new[:, 0] - offset, 6) < numpy.round(phi, 6)
#                     phi_not_0_180 = 0 == numpy.round(sang_new[:,0] - offset, self.round) % numpy.round(phi, self.round)
#                     phi_good = numpy.logical_xor(
#                         phi_0_180 & theta_0_or_180,
#                         phi_not_0_180 & ~theta_0_or_180
#                     )
#                 else:
#                     phi_good = numpy.ones(sang_new.shape[0], dtype=numpy.bool)
#                 theta_good = numpy.round(sang_new[:,1], self.round) == numpy.round(theta, self.round)
#                 psi_good = numpy.round(sang_new[:,2], self.round) < numpy.round(psi, self.round)
#                 masks_good.append(phi_good & theta_good & psi_good)
#                 if not numpy.isnan(phi):
#                     phi_bad_0_180 = numpy.round(sang_new[:, 0] - offset, self.round) >= numpy.round(phi, self.round)
#                     phi_bad = numpy.logical_xor(
#                         phi_bad_0_180 & theta_0_or_180,
#                         phi_not_0_180 & ~theta_0_or_180
#                     )
#                 else:
#                     phi_bad = numpy.ones(sang_new.shape[0], dtype=numpy.bool)
#
#                 psi_bad_not_0_180 = numpy.round(sang_new[:,2], self.round) >= numpy.round(psi, self.round)
#                 psi_bad = numpy.logical_xor(
#                     psi_good & theta_0_or_180,
#                     psi_bad_not_0_180 & ~theta_0_or_180
#                 )
#
#                 masks_bad.append(phi_bad & theta_good & psi_bad)
#
#             mask_good = numpy.zeros(sang_new.shape[0], numpy.bool)
#             for entry in masks_good:
#                 mask_good = numpy.logical_or(mask_good, entry)
#
#             mask_bad = numpy.zeros(sang_new.shape[0], numpy.bool)
#             for entry in masks_bad:
#                 mask_bad = numpy.logical_or(mask_bad, entry)
#
#             mask_not_special = ~numpy.logical_or(
#                 numpy.logical_xor(mask_good, mask_bad),
#                 numpy.logical_and(mask_good, mask_bad)
#             )
#             maski = numpy.logical_or(mask_good, mask_not_special)
#             output_mask = numpy.zeros(sang_new_raw.shape[0], dtype=numpy.bool)
#             output_mask[sang_mask] = maski
#
#             sang_new_raw[sang_mask] = sang_new
#             final_masks.append(output_mask)
#
#         final_mask = numpy.zeros(nsym*angles.shape[0], dtype=numpy.bool)
#         for entry in final_masks:
#             final_mask = numpy.logical_or(final_mask, entry)
#
#         if return_mirror not in (0, 1):
#             semi_final_mask = numpy.zeros(sang_new_raw.shape[0], dtype=numpy.bool)
#             mask1 = numpy.zeros(sang_new_raw.shape[0], dtype=numpy.bool)
#             sang_new_cpy = sang_new_raw.copy()
#             sang_new_cpy[:, 2] = 0
#             for i in range(sang_new_raw.shape[0] // nsym):
#                 mask1[i*nsym:(i+1)*nsym] = True
#                 _, idx = numpy.unique(sang_new_cpy[mask1, :], axis=0, return_index=True)
#                 for entry in idx:
#                     semi_final_mask[i*nsym+entry] = True
#                 mask1[...] = False
#         else:
#             semi_final_mask = numpy.ones(sang_new_raw.shape[0], dtype=numpy.bool)
#         # print(semi_final_mask)
#         sang_new = sang_new_raw[final_mask & semi_final_mask]
#         sang_new %= 360
#
#
#
#
#
#
#         if tolistconv:
#             return sang_new.tolist()
#         else:
#             return sang_new
#
#
#     def symmetry_neighbors(self, angles, return_mirror=0, tolistconv = True):
#
#         if self.sym[0] == 'c':
#             if int(self.sym[1:]) > 2:
#                 neighbors = [0, 1, -1]
#             elif int(self.sym[1:]) == 2:
#                 neighbors = [0, 1]
#             elif int(self.sym[1:]) == 1:
#                 neighbors = [0]
#
#         if self.sym[0] == 'd':
#             if int(self.sym[1:]) >= 3 and int(self.sym[1:]) % 2 != 0:
#                 neighbors = [0, 1, self.nsym//2-1, self.nsym//2, -2,-1]
#             elif int(self.sym[1:]) >= 3 and int(self.sym[1:]) % 2 == 0:
#                 offset = (self.nsym//2  - 4 ) // 2
#                 neighbors = [0, 1, self.nsym//2 -1, self.nsym//2 + offset + 1,self.nsym//2 + offset +2  ,self.nsym//2 + offset + 3]
#             elif int(self.sym[1:]) == 2:
#                 neighbors = [0, 1, 2]
#             elif int(self.sym[1:]) == 1:
#                 neighbors = [0, 1]
#
#         elif self.sym == 'oct':
#             neighbors  = [0,1,2,3,8,9,12,13]
#         elif self.sym == 'tet':
#             neighbors  = [0,1,2,3,4,6,7]
#         elif self.sym == 'icos':
#             neighbors = [0,1,2,3,4,6,7,11,12]
#         sang_mod = self.symmetry_related(angles, return_mirror=return_mirror, neighbors=neighbors, tolistconv=tolistconv)
#         return sang_mod
#
#
#
#     @staticmethod
#     def recmat(mat, out=None, tolistconv = True):
#         pass#IMPORTIMPORTIMPORT from math import acos,asin,atan2,degrees,pi
#         def sign(x):
#             return_array = numpy.sign(x)
#             return_array[return_array == 0] = 1
#             return return_array
#
#         mask_2_2_1 = mat[:, 2, 2] == 1.0
#         mask_0_0_0 = mat[:, 0, 0] == 0.0
#         mask_2_2_m1 = mat[:, 2, 2] == -1.0
#         mask_2_0_0 = mat[:, 2, 0] == 0.0
#         mask_0_2_0 = mat[:, 0, 2] == 0.0
#         # mask_2_2_0 = mat[:, 2, 2] == 0.0
#         theta_2_2 = numpy.arccos(mat[:, 2, 2])
#         st = sign(theta_2_2)
#         if out is None:
#             output_array = numpy.empty((mat.shape[0], 3), dtype=numpy.float64)
#         else:
#             output_array = out
#
#         output_array[mask_2_2_1 & mask_0_0_0, 0] = numpy.degrees(numpy.arcsin(mat[mask_2_2_1 & mask_0_0_0, 0, 1]))
#         output_array[mask_2_2_1 & ~mask_0_0_0, 0] = numpy.degrees(numpy.arctan2(
#             mat[mask_2_2_1 & ~mask_0_0_0, 0, 1],
#             mat[mask_2_2_1 & ~mask_0_0_0, 0, 0]
#         ))
#         output_array[mask_2_2_1 & ~mask_2_2_m1, 1] = numpy.degrees(0.0)  # theta
#         output_array[mask_2_2_1 & ~mask_2_2_m1 , 2] = numpy.degrees(0.0)  # psi
#
#         output_array[mask_2_2_m1 &  mask_0_0_0, 0] = numpy.degrees(numpy.arcsin(-mat[mask_2_2_m1 & mask_0_0_0, 0, 1] ))
#         output_array[mask_2_2_m1 & ~mask_0_0_0, 0] = numpy.degrees(numpy.arctan2(
#             -mat[mask_2_2_m1 & ~mask_0_0_0, 0, 1],
#             -mat[mask_2_2_m1 & ~mask_0_0_0, 0, 0]
#         ))
#         output_array[mask_2_2_m1 & ~mask_2_2_1, 1] = numpy.degrees(numpy.pi)
#         output_array[mask_2_2_m1 & ~mask_2_2_1, 2] = numpy.degrees(0.0)
#
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & mask_2_0_0 & (st != sign(mat[:,2,1])), 0] = numpy.degrees(1.5*numpy.pi)
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & mask_2_0_0 & (st == sign(mat[:,2,1])), 0] = numpy.degrees(0.5*numpy.pi)
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_2_0_0, 0] = numpy.degrees(numpy.arctan2(
#             st[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_2_0_0] * mat[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_2_0_0, 2, 1],
#             st[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_2_0_0] * mat[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_2_0_0, 2, 0]
#         ))
#         output_array[~mask_2_2_1 & ~mask_2_2_m1, 1] = numpy.degrees(theta_2_2[~mask_2_2_1 & ~mask_2_2_m1])
#
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & mask_0_2_0 & (st != sign(mat[:,1,2])), 2] = numpy.degrees(1.5*numpy.pi)
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & mask_0_2_0 & (st == sign(mat[:, 1, 2])), 2] = numpy.degrees(0.5 * numpy.pi)
#
#         output_array[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_0_2_0 , 2] = numpy.degrees(numpy.arctan2(
#             st[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_0_2_0]  * mat[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_0_2_0,  1, 2],
#             -st[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_0_2_0] * mat[~mask_2_2_1 & ~mask_2_2_m1 & ~mask_0_2_0, 0, 2]))
#
#         numpy.round(output_array, 12, out=output_array)
#         output_array %= 360.0
#
#         if tolistconv:
#             return output_array.tolist()
#         else:
#             return output_array
#
#
#     @staticmethod
#     def rotmatrix(angles , tolistconv = True):
#
#         newmat = numpy.zeros((len(angles), 3, 3),dtype = numpy.float64)
#         index = numpy.arange(len(angles))
#
#         cosphi = numpy.cos(numpy.radians(numpy.array(angles)[:, 0]))
#         costheta = numpy.cos(numpy.radians(numpy.array(angles)[:, 1]))
#         cospsi = numpy.cos(numpy.radians(numpy.array(angles)[ : ,2] ))
#
#         sinphi = numpy.sin(numpy.radians(numpy.array(angles)[:, 0]))
#         sintheta = numpy.sin(numpy.radians(numpy.array(angles)[:, 1]))
#         sinpsi = numpy.sin(numpy.radians(numpy.array(angles)[: ,2] ))
#
#         newmat[:,0,0] =  cospsi[index]*costheta[index]*cosphi[index] - sinpsi[index]*sinphi[index]
#         newmat[:,1,0] =     -sinpsi[index]*costheta[index]*cosphi[index] - cospsi[index]*sinphi[index]
#         newmat[:,2,0] =           sintheta[index]*cosphi[index]
#         newmat[:,0,1] =  cospsi[index]*costheta[index]*sinphi[index] + sinpsi[index]*cosphi[index]
#         newmat[:,1,1] = -sinpsi[index]*costheta[index]*sinphi[index] + cospsi[index]*cosphi[index]
#         newmat[:,2,1] =            sintheta[index]*sinphi[index]
#         newmat[:,0,2] = -cospsi[index]*sintheta[index]
#         newmat[:,1,2] =  sinpsi[index]*sintheta[index]
#         newmat[:,2,2] =            costheta[index]
#
#         if tolistconv:
#             return newmat.tolist()
#         else:
#             return newmat
#
#     def is_in_subunit(self, phi_orig, theta_orig=None, inc_mirror=1, tolistconv = True):
#
#         if theta_orig is None:
#             angles = phi_orig
#             return_single = True
#         else:
#             angles = [phi_orig, theta_orig, 0]
#             return_single = False
#
#         """
#         Input:  Before it was a projection direction specified by (phi, theta).
#                 Now it is a projection direction specified by phi(i) , theta(i) for all angles
#                 inc_mirror = 1 consider mirror directions as unique
#                 inc_mirror = 0 consider mirror directions as outside of unique range.
#         Output: True if input projection direction is in the first asymmetric subunit,
#                 False otherwise.
#         """
#         pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
#
#         angles = numpy.atleast_2d(numpy.array(angles))
#
#         condstat = numpy.zeros(numpy.shape(angles)[0], dtype = bool  )
#
#         phi = angles[:,0]
#         phi_0 =   phi >= 0.0
#         phi_ld_br_inmirr_0 = phi < self.brackets[inc_mirror][0]
#         phi_ld_br_1_0 = phi < self.brackets[inc_mirror][0]
#         theta = numpy.round(angles[:, 1], self.round)
#         theta_ldeq_br_incmirr_1 = theta  <= self.brackets[inc_mirror][1]
#         theta_ldeq_br_incmirr_3 = theta <= self.brackets[inc_mirror][3]
#         theta_180 =  (numpy.logical_and(theta ==180 , inc_mirror))
#         theta_0 = theta == 0
#
#
#         if self.sym[0] == "c" :
#             condstat[phi_0  & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_1] = True
#             condstat[theta_180] = True
#             condstat[theta_0] = True
#             condstat[~phi_0 & ~phi_ld_br_inmirr_0 & ~theta_ldeq_br_incmirr_1 & ~theta_180  & ~theta_0] = False
#
#         elif self.sym[0] == "d" and (old_div(self.nsym, 2)) % 2 == 0:
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_1] = True
#             condstat[theta_0] = True
#             condstat[~phi_0 & ~phi_ld_br_inmirr_0 & ~theta_ldeq_br_incmirr_1 & ~theta_0] = False
#
#         elif self.sym[0] == "d" and (old_div(self.nsym, 2)) % 2 == 1:
#             phib = old_div(360.0, self.nsym)
#             condstat[numpy.logical_and( (theta_ldeq_br_incmirr_1 & phi_0 & phi_ld_br_1_0), inc_mirror)] = True
#             condstat[ theta_ldeq_br_incmirr_1 & phi_0 & phi_ld_br_1_0 & \
#                 ( numpy.logical_or(numpy.logical_and(  (phi >= old_div(phib, 2)) , (phi < phib)) , numpy.logical_and( (phi >= phib), (phi <= phib + old_div(phib, 2)) ))) ] = True
#             condstat[theta_ldeq_br_incmirr_1 & phi_0 & phi_ld_br_1_0 & theta_0] = True
#             condstat[~(theta_ldeq_br_incmirr_1 & phi_0 & phi_ld_br_1_0) & theta_0] = True
#             condstat[~(theta_ldeq_br_incmirr_1 & phi_0 & phi_ld_br_1_0) &  ~theta_0] = False
#
#         elif ((self.sym[:3] == "oct") or (self.sym[:4] == "icos")):
#             tmphi = numpy.minimum(phi, self.brackets[inc_mirror][2] - phi)
#             baldwin_lower_alt_bound = old_div(
#                 (
#                     old_div(
#                         numpy.sin(numpy.radians(old_div(
#                             self.brackets[inc_mirror][2],
#                             2.0
#                         ) - tmphi)),
#                         numpy.tan(numpy.radians(self.brackets[inc_mirror][1]))
#                     ) +
#                     old_div(
#                         numpy.sin(numpy.radians(tmphi)),
#                         numpy.tan(numpy.radians(self.brackets[inc_mirror][3]))
#                     )
#                 ),
#                 numpy.sin(numpy.radians(old_div(
#                     self.brackets[inc_mirror][2],2.0)
#                 ))
#             )
#             baldwin_lower_alt_bound = numpy.degrees(numpy.arctan(old_div(1.0, baldwin_lower_alt_bound)))
#
#             numpy.round(baldwin_lower_alt_bound, self.round, out=baldwin_lower_alt_bound)
#             # print(baldwin_lower_alt_bound)
#             condstat[phi_0 & phi_ld_br_inmirr_0  & theta_ldeq_br_incmirr_3 & (baldwin_lower_alt_bound >= theta)] = True
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3 & ~(baldwin_lower_alt_bound >= theta)] = False
#             condstat[~phi_0 & ~phi_ld_br_inmirr_0 & ~theta_ldeq_br_incmirr_3 & ~theta_0 ] = False
#             condstat[theta_0] = True
#             condstat[
#                 (theta <= numpy.round(self.brackets[inc_mirror][3], self.round))
#                 & (phi == numpy.round(self.brackets[inc_mirror][0], self.round))
#                 ] = True
#             condstat[
#                 (theta == baldwin_lower_alt_bound)
#                 & (phi > self.brackets[0][0])
#                 ] = False
#             condstat[phi == numpy.round(self.brackets[inc_mirror][2], self.round)] = False
#
#             condstat[  theta_0  & (phi ==  numpy.round(self.brackets[inc_mirror][3] ,6)) ] = False
#
#         elif (self.sym[:3] == "tet"):
#             print('ok')
#             tmphi = numpy.minimum(phi, self.brackets[inc_mirror][2] - phi)
#             baldwin_lower_alt_bound_1 = \
#                 old_div(
#                     (old_div(numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) - tmphi)),
#                              numpy.tan(numpy.radians(self.brackets[inc_mirror][1]))) \
#                      + old_div(numpy.sin(numpy.radians(tmphi)), numpy.tan(numpy.radians(self.brackets[inc_mirror][3])))) \
#                     , numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0))))
#             is_zero = baldwin_lower_alt_bound_1 == 0
#             baldwin_lower_alt_bound_1[~is_zero] = numpy.degrees(numpy.arctan(old_div(1.0, baldwin_lower_alt_bound_1[~is_zero])))
#             baldwin_lower_alt_bound_1[is_zero] = self.brackets[inc_mirror][3]
#             baldwin_upper_alt_bound_2 = \
#                 old_div(
#                     (old_div((numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) - tmphi))),
#                              (numpy.tan(numpy.radians(self.brackets[inc_mirror][1]))))
#                      + old_div((numpy.sin(numpy.radians(tmphi))),
#                                numpy.tan(numpy.radians(old_div(self.brackets[inc_mirror][3], 2.0))))) \
#                     , (numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0)))))
#             baldwin_upper_alt_bound_2 = numpy.degrees(numpy.arctan(old_div(1.0, baldwin_upper_alt_bound_2)))
#
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3  & \
#                       numpy.logical_and((baldwin_lower_alt_bound_1 > theta) , inc_mirror)] = True
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3  & \
#                      ~numpy.logical_and((baldwin_lower_alt_bound_1 > theta) , inc_mirror) & (numpy.round(baldwin_upper_alt_bound_2, self.round) <= numpy.round(theta, self.round))] = False
#             condstat[ phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3  & \
#                      ~numpy.logical_and((baldwin_lower_alt_bound_1 > theta) , inc_mirror) & ~(numpy.round(baldwin_upper_alt_bound_2, self.round) <= numpy.round(theta, self.round))] = True
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3  & ~(baldwin_lower_alt_bound_1 > theta) &  theta_0] = True
#             condstat[phi_0 & phi_ld_br_inmirr_0 & theta_ldeq_br_incmirr_3  & ~(baldwin_lower_alt_bound_1 > theta) &  ~theta_0] = False
#             condstat[~phi_0 & ~phi_ld_br_inmirr_0 & ~theta_ldeq_br_incmirr_3 &  theta_0] = True
#             condstat[~phi_0 & ~phi_ld_br_inmirr_0 & ~theta_ldeq_br_incmirr_3 & ~theta_0] = False
#             condstat[theta_0] = True
#             condstat[
#                 (numpy.round(theta, self.round) == numpy.round(baldwin_lower_alt_bound_1, self.round))
#                 & (phi > self.brackets[0][0] / 2.0)
#                 ] = False
#             condstat[
#                 (numpy.round(theta, self.round) == numpy.round(baldwin_upper_alt_bound_2, self.round))
#                 & (theta < self.brackets[inc_mirror][3])
#                 ] = True
#         else:
#             global_def.ERROR("unknown symmetry", "symclass: is_in_subunit", 1)
#         if tolistconv:
#             condstat = condstat.tolist()
#         else:
#             condstat = condstat
#
#         if return_single:
#             return condstat
#         else:
#             return condstat[0][0]
#
#
#     def reduce_anglesets(self, angles,inc_mirror=1, tolistconv = True):
#         """
#           Input is either list or lists [[phi,thet,psi],[],[]] or a triplet [phi,thet,psi]
#                 inc_mirror = 1 consider mirror directions as unique
#                 inc_mirror = 0 consider mirror directions as outside of unique range.
#           It will map all triplets to the first asymmetric subunit.
#         """
#         if inc_mirror == 1:
#            return_mirror = 0
#         else:
#            return_mirror = 2
#
#         sym_angles = self.symmetry_related(angles, return_mirror = return_mirror, tolistconv=False)
#         sym_angles = sym_angles.tolist()
#         subunits   = self.is_in_subunit(sym_angles, inc_mirror = inc_mirror)
#         reduced_anglesets = numpy.array(sym_angles)[subunits]
#
#         if tolistconv:
#             return reduced_anglesets.tolist()
#         else:
#             return reduced_anglesets



class TestSymClassIsInSubunitC(unittest.TestCase):
    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_c1_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)

    def test_c1_sym_no_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        # [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles)
        self.assertTrue(results, expected_results)

    def test_c1_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_c1_sym_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        # [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles)
        self.assertTrue(results, expected_results)

    def test_c5_sym_no_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
        self.assertEqual(results, expected_results)

    def test_c5_sym_no_mirror_theta_larger_90_andsmaller_180_degrees_phi_smaller_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(180)]
        # [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)

        self.assertTrue(results, expected_results)

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
        self.assertTrue(results, expected_results)

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        # [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c5').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
        # for i, k, j in zip(results, expected_results, angles):
        #     if i != k:
        #         print(self.output_template_angles.format(i, k, j))
        self.assertTrue(results, expected_results)

    def test_c5_sym_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
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
            is_in = symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
        self.assertTrue(results, expected_results)


    def test_c5_sym_mirror_theta_smaller_equals_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
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
            is_in =  symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles)
        self.assertTrue(results, expected_results)


class TestSymClassIsInSubunit(unittest.TestCase):

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_wrong_sym_crashes_problem(self):
        return_new = fu.symclass('c1')
        return_old = oldfu.symclass('c1')
        return_new.sym = 'foobar'
        return_old.sym = 'foobar'
        self.assertIsNone(return_new.is_in_subunit(0, 0))
        self.assertIsNone(return_old.is_in_subunit(0, 0))

    # -------------------------------------------------------[ Test Cases for C symmetry ]
    def test_newmethod_c1_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c1').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)


    def test_newmethod_c1_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c1').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)

    def test_newmethod_c4_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in =  symclass('c4').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('c4').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)


    def test_newmethod_c4_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        # [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c4').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('c4').is_in_subunit(angles, inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results,expected_results)



    def test_newmethod_c1_sym_no_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        # [angles.append([entry, 180, 0]) for entry in range(360)]
        print (fu.symclass('c1').brackets)
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c1').is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles,inc_mirror=0)
        self.assertTrue(results, expected_results)


    def test_newmethod_c1_sym_mirror_theta_larger_90_degrees_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass('c1').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c1').is_in_subunit(angles,inc_mirror=1)
        self.assertTrue(results, expected_results)




    def test_newmethod_c5_sym_mirror_theta_smaller_equals_90_degrees_phi_bigger_equals_72_should_return_False(self):
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('c5').is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = fu.symclass('c5').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)



    # -------------------------------------------------------[ Test Cases for D symmetry ]
    def test_newmethod_d1_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('d1').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('d1').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)



    def test_newmethod_d1_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]

        results = []
        for phi, theta, psi in angles:
            is_in = symclass('d1').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('d1').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


    def test_newmethod_d5_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('d5').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('d5').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)



    def test_newmethod_d5_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('d5').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('d5').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    # -------------------------------------------------------[ Test Cases for tet symmetry ]
    def test_newmethod_tet_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('tet').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('tet').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_newmethod_tet_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('tet').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('tet').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


    # -------------------------------------------------------[ Test Cases for Oct symmetry ]
    def test_newmethod_oct_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('oct').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('oct').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_newmethod_oct_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('oct').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('oct').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    # -------------------------------------------------------[ Test Cases for Icos symmetry ]
    def test_newmethod_icos_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('icos').is_in_subunit(phi, theta, inc_mirror=1)
            results.append(is_in)
        expected_results = fu.symclass('icos').is_in_subunit(angles,inc_mirror=1)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)

    def test_newmethod_icos_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        angles = [[entry, thet, 0.0] for entry in range(120) for thet in range(55)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass('icos').is_in_subunit(phi, theta, inc_mirror=0)
            results.append(is_in)
        expected_results = fu.symclass('icos').is_in_subunit(angles,inc_mirror=0)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        self.assertEqual(results, expected_results)


class TestSymClassSymmetryRelatedC(unittest.TestCase):

    def test_c1_sym_zero(self):
        new_result = fu.symclass('c1').symmetry_related([0, 0, 0])
        old_result = symclass('c1').symmetry_related([0, 0, 0])

        self.assertEqual(new_result,[[0, 0, 0]])
        self.assertEqual(old_result, [[0, 0, 0]])

    def test_c5_sym_zero(self):
        new_result = fu.symclass('c5').symmetry_related([0, 0, 0])
        old_result = symclass('c5').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0]])
        self.assertEqual(old_result, [[0, 0, 0]])


    def test_c6_sym_zero(self):
        new_result = fu.symclass('c6').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('c6').symmetry_related([0, 0, 0])
        self.assertEqual(new_result,[[0, 0, 0]])
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

    def test_c6_sym_180(self):
        new_result = fu.symclass('c6').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('c6').symmetry_related([0, 180, 0])
        self.assertEqual(new_result,[[0, 180, 0]])
        self.assertEqual(old_result, [[0, 180, 0]])


    def test_c1_sym_90(self):
        new_result = fu.symclass('c1').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('c1').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0]])
        self.assertEqual(old_result, [[0, 90, 0]])

    def test_c5_sym_90(self):
        new_result = fu.symclass('c5').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('c5').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]]   )
        self.assertEqual(old_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]]   )

    def test_c6_sym_90(self):
        new_result = fu.symclass('c6').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('c6').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0], [60, 90, 0], [60*2, 90, 0], [60*3, 90, 0], [60*4, 90, 0] , [60*5, 90, 0] ]   )
        self.assertEqual(old_result, [[0, 90, 0], [60, 90, 0], [60*2, 90, 0], [60*3, 90, 0], [60*4, 90, 0] , [60*5, 90, 0] ]   )


    def test_c1_sym_42(self):
        new_result = fu.symclass('c1').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('c1').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0]])
        self.assertEqual(old_result, [[0, 42, 0]])

    def test_c5_sym_42(self):
        new_result = symclass_mod('c5').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('c5').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0], [72, 42, 0], [72*2, 42, 0], [72*3, 42, 0], [72*4, 42, 0]]   )
        self.assertEqual(old_result, [[0, 42, 0], [72, 42, 0], [72*2, 42, 0], [72*3, 42, 0], [72*4, 42, 0]]   )

    def test_c6_sym_42(self):
        new_result = fu.symclass('c6').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('c6').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0], [60, 42, 0], [60*2, 42, 0], [60*3, 42, 0], [60*4, 42, 0] , [60*5, 42, 0] ]   )
        self.assertEqual(old_result, [[0, 42, 0], [60, 42, 0], [60*2, 42, 0], [60*3, 42, 0], [60*4, 42, 0] , [60*5, 42, 0] ]   )




class TestSymClassSymmetryRelatedD(unittest.TestCase):

    def test_d1_sym_zero(self):
        new_result = fu.symclass('d1').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('d1').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0], [0, 180, 180]] )
        # self.assertEqual(old_result, [[0, 0, 0], [0, 180, 0]] )


    def test_d5_sym_zero(self):
        new_result = fu.symclass('d5').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('d5').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0], [0, 180, 180]] )
        self.assertEqual(old_result, [[0, 0, 0], [0, 180, 180]] )


    def test_d6_sym_zero(self):
        new_result = fu.symclass('d6').symmetry_related([0, 0, 0])
        old_result = oldfu.symclass('d6').symmetry_related([0, 0, 0])
        self.assertEqual(new_result, [[0, 0, 0], [0, 180, 180]] )
        # self.assertEqual(old_result, [[0, 0, 0], [0, 180, 180]] )


    def test_d1_sym_90(self):
        new_result = fu.symclass('d1').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('d1').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0]])
        self.assertEqual(old_result, [[0, 90, 0]])

    def test_d5_sym_90(self):
        new_result = fu.symclass('d5').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('d5').symmetry_related([0, 90, 0])

        self.assertEqual(new_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]]   )
        self.assertEqual(old_result, [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]]   )

    def test_d6_sym_90(self):
        new_result = fu.symclass('d6').symmetry_related([0, 90, 0])
        old_result = oldfu.symclass('d6').symmetry_related([0, 90, 0])
        self.assertEqual(new_result, [[0, 90, 0], [60, 90, 0], [60*2, 90, 0], [60*3, 90, 0], [60*4, 90, 0] , [60*5, 90, 0] ]   )
        self.assertEqual(old_result, [[0, 90, 0], [60, 90, 0], [60*2, 90, 0], [60*3, 90, 0], [60*4, 90, 0] , [60*5, 90, 0] ]   )

    """ Old version is wrong, have corrected it  """
    def test_d1_sym_180(self):
        new_result = fu.symclass('d1').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('d1').symmetry_related([0, 180, 0])
        self.assertEqual(new_result, [[0, 180, 0], [0, 0, 180]] )
        # self.assertEqual(old_result, [[0, 180, 0], [0,  180, 180]] )

    """ Old version is wrong, have corrected it  """
    def test_d5_sym_180(self):
        new_result = fu.symclass('d5').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('d5').symmetry_related([0, 180, 0])
        self.assertEqual(new_result, [[0, 180, 0], [0, 0, 180]] )
        # self.assertEqual(old_result, [[0, 180, 0], [0, 180, 180]] )


    def test_d6_sym_180(self):
        new_result = fu.symclass('d6').symmetry_related([0, 180, 0])
        old_result = oldfu.symclass('d6').symmetry_related([0, 180, 0])
        self.assertEqual(new_result, [[0, 180, 0], [0, 0, 180]] )
        # self.assertEqual(old_result, [[0, 0, 0], [0, 180, 0]] )

    def test_d1_sym_42(self):
        new_result = fu.symclass('d1').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('d1').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0], [0.0, 138.0, 180.0]])
        self.assertEqual(old_result, [[0, 42, 0], [0.0, 138.0, 180.0]])

    def test_d5_sym_42(self):
        new_result = fu.symclass('d5').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('d5').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0], [72.0, 42, 0],  [144.0, 42, 0],  [216.0, 42, 0],  [288.0, 42, 0], [0.0, 138.0, 180.0], \
                                      [288.0, 138.0, 180.0],  [216.0, 138.0, 180.0],  [144.0, 138.0, 180.0], [72.0, 138.0, 180.0]] )
        self.assertEqual(old_result, [[0, 42, 0], [72.0, 42, 0],  [144.0, 42, 0],  [216.0, 42, 0],  [288.0, 42, 0], [0.0, 138.0, 180.0], \
                                      [288.0, 138.0, 180.0],  [216.0, 138.0, 180.0],  [144.0, 138.0, 180.0], [72.0, 138.0, 180.0]] )
    def test_d6_sym_42(self):
        new_result = fu.symclass('d6').symmetry_related([0, 42, 0])
        old_result = oldfu.symclass('d6').symmetry_related([0, 42, 0])
        self.assertEqual(new_result, [[0, 42, 0], [60.0, 42, 0],  [120.0, 42, 0],  [180.0, 42, 0],  [240.0, 42, 0],  [300.0, 42, 0], \
                                     [0.0, 138.0, 0.0], [300.0, 138.0, 0.0],  [240.0, 138.0, 0.0],  [180.0, 138.0, 0.0],  [120.0, 138.0, 0.0], [60.0, 138.0, 0.0]]  )
        self.assertEqual(old_result, [[0, 42, 0], [60.0, 42, 0],  [120.0, 42, 0],  [180.0, 42, 0],  [240.0, 42, 0],  [300.0, 42, 0], \
                                     [0.0, 138.0, 0.0], [300.0, 138.0, 0.0],  [240.0, 138.0, 0.0],  [180.0, 138.0, 0.0],  [120.0, 138.0, 0.0], [60.0, 138.0, 0.0]]  )




class TestSymClassSymmetryNeighbors(unittest.TestCase):

    def test_c_sym_c4_neighbors(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('c4').symmetry_neighbors(angles)
        return_values = fu.symclass('c4').symmetry_neighbors(angles)
        self.assertEqual(return_values,expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('c4').symmetry_neighbors(angles)
        return_values = oldfu.symclass('c4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)


    def test_c_sym_c5_neighbors(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('c5').symmetry_neighbors(angles)
        return_values = fu.symclass('c5').symmetry_neighbors(angles)
        self.assertEqual(return_values,expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('c5').symmetry_neighbors(angles)
        return_values = oldfu.symclass('c5').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

    def test_d_sym_d4_neighbors(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d4').symmetry_neighbors(angles)
        return_values = fu.symclass('d4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d4').symmetry_neighbors(angles)
        return_values = oldfu.symclass('d4').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

    def test_d_sym_d5_neighbors(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d5').symmetry_neighbors(angles)
        return_values = fu.symclass('d5').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d5').symmetry_neighbors(angles)
        return_values = oldfu.symclass('d5').symmetry_neighbors(angles)
        self.assertEqual(return_values, expected_return_values)


    def test_tet_sym_neighbors(self):
        # angles = [[idx1, idx2, idx3] for idx1 in range(180) for idx2 in range(180) for idx3 in range(120) ]
        angles = [[0, 0 ,0], [0,180,0], [0,90,0], [90,0,0], [90,90,0], [90,180,0] ]
        expected_return_values = symclass_mod('tet').symmetry_neighbors(angles)
        return_values = fu.symclass('tet').symmetry_neighbors(angles)

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())

        # # angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[0, 0, 0], [0, 180, 0], [0, 90, 0], [90, 0, 0], [90, 90, 0], [90, 180, 0]]
        # expected_return_values = symclass_mod('tet').symmetry_neighbors(angles)
        # return_values = oldfu.symclass('tet').symmetry_neighbors(angles)
        # self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())


    def test_oct_sym_neighbors(self):
        # angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        angles = [[0, 0, 0], [0, 180, 0], [0, 90, 0], [90, 0, 0], [90, 90, 0], [90, 180, 0]]
        expected_return_values = symclass_mod('oct').symmetry_neighbors(angles)
        return_values = fu.symclass('oct').symmetry_neighbors(angles)
        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())



        # angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        # expected_return_values = symclass_mod('oct').symmetry_neighbors(angles)
        # return_values = oldfu.symclass('oct').symmetry_neighbors(angles)
        # self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())
        #

    def test_icos_sym_neighbors(self):
        # angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        angles = [[0, 0, 0], [0, 180, 0], [0, 90, 0], [90, 0, 0], [90, 90, 0], [90, 180, 0]]
        expected_return_values = symclass_mod('icos').symmetry_neighbors(angles)
        return_values = fu.symclass('icos').symmetry_neighbors(angles)
        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())

        # angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        # expected_return_values = symclass_mod('icos').symmetry_neighbors(angles)
        # return_values = oldfu.symclass('icos').symmetry_neighbors(angles)
        # self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol = 1  ).all())





class TestSymClassSymmetryRelated(unittest.TestCase):

    def test_c_sym_related_c4(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        # angles = [[0, 0, 0], [0, 180, 0], [90, 45, 29]]

        expected_return_values = symclass_mod('c4').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('c4').symmetry_related(angle))

        print(angles)
        print(expected_return_values)

        self.assertEqual(return_values,expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]

        expected_return_values = symclass_mod('c4').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('c4').symmetry_related(angle))

        self.assertEqual(return_values, expected_return_values)

    def test_c_sym_related_c5(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        # angles = [[0, 0, 0], [0, 180, 0], [90, 45, 29]]

        expected_return_values = symclass_mod('c5').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('c5').symmetry_related(angle))

        print(angles)
        print(expected_return_values)

        self.assertEqual(return_values, expected_return_values)

        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]

        expected_return_values = symclass_mod('c5').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('c5').symmetry_related(angle))

        self.assertEqual(return_values, expected_return_values)


    def test_d_sym_related_d4(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]
        # angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]

        import time; start=time.time()
        expected_return_values = symclass_mod('d4').symmetry_related(angles)
        print(time.time()-start); start=time.time()
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('d4').symmetry_related(angle))
        print(time.time() - start);
        start = time.time()
        self.assertEqual(sorted(return_values),sorted(expected_return_values))

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d4').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('d4').symmetry_related(angle))

        self.assertEqual(sorted(return_values),sorted(expected_return_values))

    def test_d_sym_related_d5(self):
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]
        # angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]

        import time; start=time.time()
        expected_return_values = symclass_mod('d5').symmetry_related(angles)
        print(time.time()-start); start=time.time()
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('d5').symmetry_related(angle))
        print(time.time() - start);
        start = time.time()
        self.assertEqual(sorted(return_values),sorted(expected_return_values))

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        expected_return_values = symclass_mod('d5').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('d5').symmetry_related(angle))

        self.assertEqual(sorted(return_values),sorted(expected_return_values))


    def test_tet_sym_related(self):
        # angles = [[idx1, idx2,idx3] for idx1 in range(90) for idx2 in range(90) for idx3 in range(5)]
        # angles = [[0, 0, 0], [0, 180, 0], [0, 90, 0], [90, 0, 0], [90, 90, 0], [90, 180, 0]]
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]

        import time; start=time.time()
        expected_return_values = symclass_mod('tet').symmetry_related(angles)
        print(time.time()-start); start=time.time()
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('tet').symmetry_related(angle))
        print(time.time() - start);
        start = time.time()

        # for i in range (len(return_values)):
        #      if ~(numpy.isclose(numpy.array(return_values)[i], numpy.array(expected_return_values)[i] , atol = 1 ).all())   :
        #          print(i)
        #          print(angles[i//symclass_mod('tet').nsym])
        #          print(expected_return_values[i])
        #          print(return_values[i])

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())


        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]
        expected_return_values = symclass_mod('tet').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('tet').symmetry_related(angle))

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())

    def test_oct_sym_related(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]

        import time; start=time.time()
        expected_return_values = symclass_mod('oct').symmetry_related(angles)
        print(time.time()-start); start=time.time()
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('oct').symmetry_related(angle))
        print(time.time() - start);
        start = time.time()

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]
        expected_return_values = symclass_mod('oct').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('oct').symmetry_related(angle))

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())

    def test_icos_sym_related(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]

        import time; start=time.time()
        expected_return_values = symclass_mod('icos').symmetry_related(angles)
        print(time.time()-start); start=time.time()
        return_values = []
        for angle in angles:
            return_values.extend(fu.symclass('icos').symmetry_related(angle))
        print(time.time() - start);
        start = time.time()

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())

        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(180)]
        # angles = [[90, 0, 0], [90, 90, 29], [90, 45, 29]]
        expected_return_values = symclass_mod('icos').symmetry_related(angles)
        return_values = []
        for angle in angles:
            return_values.extend(oldfu.symclass('icos').symmetry_related(angle))

        self.assertTrue(numpy.isclose(numpy.array(return_values), numpy.array(expected_return_values), atol=1).all())


class TestSymClassReduceAngleSets(unittest.TestCase):

    def test_reduce_anglesets_c1_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('c1').reduce_anglesets(angles)
        newangles = fu.symclass('c1').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_c4_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('c4').reduce_anglesets(angles)
        newangles = fu.symclass('c4').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_c5_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('c5').reduce_anglesets(angles)
        newangles = fu.symclass('c5').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_d1_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('d1').reduce_anglesets(angles)
        newangles = fu.symclass('d1').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_d4_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('d4').reduce_anglesets(angles)
        newangles = fu.symclass('d4').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_d5_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('d5').reduce_anglesets(angles)
        newangles = fu.symclass('d5').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_tet_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(180)]

        oldangles = oldfu.symclass('tet').reduce_anglesets(angles, inc_mirror= 0)
        newangles = fu.symclass('tet').reduce_anglesets(angles, inc_mirror= 0)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_icos_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(50) for theta in range(91)]

        oldangles = oldfu.symclass('icos').reduce_anglesets(angles)
        newangles = fu.symclass('icos').reduce_anglesets(angles)
        self.assertEqual(oldangles, newangles)

    def test_reduce_anglesets_oct_should_return_equal_object(self):
        angles = [[phi, theta, 0] for phi in range(360) for theta in range(91)]

        oldangles = oldfu.symclass('oct').reduce_anglesets(angles, inc_mirror= 0)
        newangles = fu.symclass('oct').reduce_anglesets(angles, inc_mirror=0)
        self.assertEqual(oldangles, newangles)


    def test_reduce_anglesets_new_oct_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        # angles = [[entry, thet, psi]for entry in range(120) for thet in range(0,55) for psi in range(20)   ]

        angles = [[idx1, idx2, 0] for idx1 in range(72) for idx2 in range(39)]
        angles = [[60, 35.264389682754654, 0]]
        # angles = [[60,  70.528779365509308 * idx/360., 0] for idx in range(360)]
        sym = 'tet'
        s1 = fu.symclass(sym)
        s2 = symclass_mod(sym)
        results = []
        expected_results = []
        sum = 0
        for ang in angles:
            r= s1.reduce_anglesets(ang, inc_mirror=0)
            results.append(r)
            a = s2.reduce_anglesets([ang], inc_mirror=0)
            sum += len(a)

        expected_results = s2.reduce_anglesets(angles, inc_mirror=1)
        print(numpy.array(angles).shape)
        print(numpy.array(results).shape)
        print(numpy.array(expected_results).shape)

        self.assertTrue(numpy.isclose(numpy.array(results), numpy.array(expected_results), atol=1).all())

        print("knock knock")

        # angles = [[entry, thet, psi] for entry in range(120) for thet in range(1,55)  for psi in range(34)]
        # # [angles.append([entry, 0, 0]) for entry in range(360)]
        # results = []
        # for phi, theta, psi in angles:
        #     is_in = oldfu.symclass('oct').reduce_anglesets(angles, inc_mirror=1)
        #     results.append(is_in)
        # expected_results = self.reduce_anglesets_new(angles, oldfu.symclass('oct').sym, oldfu.symclass('oct').nsym,
        #                                           oldfu.symclass('oct').brackets, oldfu.symclass('oct').symatrix, inc_mirror=1)
        # self.assertTrue(numpy.array_equal(numpy.array(results), numpy.array(expected_results)))


class symclass(object):
    pass  # IMPORTIMPORTIMPORT import numpy as np

    def __init__(self, sym):
        """
          sym: cn, dn, oct, tet, icos
        """
        pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt, pi
        # from utilities import get_sym, get_symt
        pass  # IMPORTIMPORTIMPORT from string import lower
        self.sym = sym.lower()
        if (self.sym[0] == "c"):
            self.nsym = int(self.sym[1:])
            if (self.nsym < 1):  global_def.ERROR("For Cn symmetry, we need n>0", "symclass", 1)
            self.brackets = [[ old_div(360., self.nsym), 90.0, old_div(360., self.nsym), 90.0],
                             [ old_div(360., self.nsym), 180.0, old_div(360., self.nsym), 180.0]]
            self.symangles = []
            for i in range(self.nsym):
                self.symangles.append([0.0, 0.0, i * old_div(360., self.nsym)])

        elif (self.sym[0] == "d"):
            self.nsym = 2 * int(self.sym[1:])
            if (self.nsym < 1):  global_def.ERROR("For Dn symmetry, we need n>0", "symclass", 1)
            self.brackets = [[old_div(360., self.nsym), 90.0, old_div(360., self.nsym), 90.0],
                             [old_div(360., self.nsym) * 2, 90.0, old_div(360., self.nsym) * 2, 90.0]]
            self.symangles = []
            for i in range(old_div(self.nsym, 2)):
                self.symangles.append([0.0, 0.0, i * old_div(360., self.nsym) * 2])
            for i in range(old_div(self.nsym, 2)):
                self.symangles.append(
                    [0.0, 180.0, (i * old_div(360., self.nsym) * 2 + 180.0 * (int(self.sym[1:]) % 2)) % 360.0])

        elif (self.sym[:3] == "oct"):
            self.nsym = 24
            ncap = 4
            cap_sig = old_div(360.0,  ncap)
            alpha = numpy.degrees(math.acos(
               old_div(1.0 , (numpy.sqrt(3.0) * numpy.tan(2 * old_div(old_div(numpy.pi,ncap), 2.0))))) )  # also platonic_params["alt_max"]
            theta = numpy.degrees(0.5 * math.acos(old_div(numpy.cos(numpy.radians(cap_sig)) ,
                                                          (1.0 - numpy.cos(numpy.radians(cap_sig))))))  # also platonic_params["theta_c_on_two"]
            self.brackets = [[ old_div(180., ncap), theta, cap_sig, alpha], [old_div(360., ncap), theta, cap_sig, alpha]]
            self.symangles = [[0.0, 0.0, float(i)] for i in range(0, 271, 90)]
            for i in range(0, 271, 90):
                for j in range(0, 271, 90):
                    self.symangles.append([float(j), 90.0, float(i)])
            for i in range(0, 271, 90):  self.symangles.append([0.0, 180.0, float(i)])

        elif (self.sym[:3] == "tet"):
            self.nsym = 12
            ncap = 3
            cap_sig = old_div(360.0,  ncap)
            alpha = numpy.degrees(math.acos(old_div(1.0,(numpy.sqrt(3.0) * numpy.tan(2 * old_div(old_div(numpy.pi, ncap), 2.0))))))  # also platonic_params["alt_max"]
            theta = numpy.degrees(0.5 * math.acos(old_div(numpy.cos(numpy.radians(cap_sig)),
                                                           (1.0 - numpy.cos(numpy.radians(cap_sig))))))  # also platonic_params["theta_c_on_two"]
            self.brackets = [[old_div(360.0,  ncap), theta, cap_sig, alpha], [old_div(360.0,  ncap), theta, cap_sig, alpha]]
            lvl1 = numpy.degrees(math.acos(old_div(-1.0, 3.0)))  # There  are 3 faces at this angle
            self.symangles = [[0., 0., 0.], [0., 0., 120.], [0., 0., 240.]]
            for l1 in range(0, 241, 120):
                for l2 in range(60, 301, 120):
                    self.symangles.append([float(l1), lvl1, float(l2)])

            """Multiline Comment4"""

        elif (self.sym[:4] == "icos"):
            self.nsym = 60
            ncap = 5
            cap_sig = old_div(360.0,  ncap)
            alpha = numpy.degrees(math.acos( old_div(1.0, (numpy.sqrt(3.0) * numpy.tan(2 * old_div( old_div(numpy.pi, ncap), 2.0))  ))))  # also platonic_params["alt_max"]
            theta = numpy.degrees(0.5 * math.acos(old_div( numpy.cos(numpy.radians(cap_sig)),
                                                    (1.0 - numpy.cos(numpy.radians(cap_sig)))) ))  # also platonic_params["theta_c_on_two"]
            self.brackets = [[36., theta, cap_sig, alpha], [72., theta, cap_sig, alpha]]
            lvl1 = numpy.degrees(math.atan(2.0))  # there are 5 pentagons with centers at this height (angle)
            lvl2 = 180.0 - lvl1  # there are 5 pentagons with centers at this height (angle)
            self.symangles = [[0.0, 0.0, float(i)] for i in range(0, 288 + 1, 72)]
            for l1 in range(0, 288 + 1, 72):
                for l2 in range(36, 324 + 1, 72):
                    self.symangles.append([float(l1), lvl1, float(l2)])
            for l1 in range(36, 324 + 1, 72):
                for l2 in range(0, 288 + 1, 72):
                    self.symangles.append([float(l1), lvl2, float(l2)])
            for i in range(0, 288 + 1, 72):  self.symangles.append([0.0, 180.0, float(i)])

        else:
            global_def.ERROR("Unknown symmetry", "symclass", 1)

        #
        self.transform = []
        self.symatrix = []
        for args in self.symangles:
            self.transform.append(
                EMAN2_cppwrap.Transform({"type": "spider", "phi": args[0], "theta": args[1], "psi": args[2]}))
            self.symatrix.append(self.rotmatrix(args[0], args[1], args[2]))

    def is_in_subunit(self, phi, theta, inc_mirror=1):
        """
        Input: a projection direction specified by (phi, theta).
                inc_mirror = 1 consider mirror directions as unique
                inc_mirror = 0 consider mirror directions as outside of unique range.
        Output: True if input projection direction is in the first asymmetric subunit,
                False otherwise.
        """
        pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
        if self.sym[0] == "c":
            if phi >= 0.0 and phi < self.brackets[inc_mirror][0] and theta <= self.brackets[inc_mirror][1]:
                return True
            elif theta == 180 and inc_mirror:
                return True
            elif theta == 0:
                return True
            else:
                return False

        elif self.sym[0] == "d" and (old_div(self.nsym,2)) % 2 == 0:
            if phi >= 0.0 and phi < self.brackets[inc_mirror][0] and theta <= self.brackets[inc_mirror][1]:
                return True
            elif theta == 0:
                return True
            else:
                return False

        elif self.sym[0] == "d" and (old_div(self.nsym,2)) % 2 == 1:
            if theta <= self.brackets[inc_mirror][1]:
                phib = old_div(360.0, self.nsym)
                if phi >= 0.0 and phi < self.brackets[1][0]:
                    if inc_mirror == 1:
                        return True
                    elif (phi >= old_div(phib, 2) and phi < phib) or (phi >= phib and phi <= phib + old_div(phib, 2)):
                        return True
                    elif theta == 0:
                        return True
            if theta == 0:
                return True
            else:
                return False

        elif ((self.sym[:3] == "oct") or (self.sym[:4] == "icos")):
            if (phi >= 0.0 and phi < self.brackets[inc_mirror][0] and theta <= self.brackets[inc_mirror][3]):
                tmphi = min(phi, self.brackets[inc_mirror][2] - phi)
                baldwin_lower_alt_bound = \
                  old_div(( old_div(numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) - tmphi)), numpy.tan(
                        numpy.radians(self.brackets[inc_mirror][1])) )
                            + old_div(numpy.sin(numpy.radians(tmphi)) , numpy.tan(numpy.radians(self.brackets[inc_mirror][3])))),
                          numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0))) )
                baldwin_lower_alt_bound = numpy.degrees(math.atan(old_div(1.0, baldwin_lower_alt_bound)))
                if (baldwin_lower_alt_bound > theta):
                    return True
                else:
                    return False
            elif theta == 0:
                return True
            else:
                # print "phi",self.brackets
                return False

        elif (self.sym[:3] == "tet"):
            if (phi >= 0.0 and phi < self.brackets[inc_mirror][0] and theta <= self.brackets[inc_mirror][3]):
                tmphi = min(phi, self.brackets[inc_mirror][2] - phi)
                baldwin_lower_alt_bound = \
                   old_div(
                       (old_div(numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) - tmphi)),
                               numpy.tan(numpy.radians(self.brackets[inc_mirror][1]))) \
                     + old_div(numpy.sin(numpy.radians(tmphi)), numpy.tan(numpy.radians(self.brackets[inc_mirror][3]))))  \
                    , numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0))))
                baldwin_lower_alt_bound = numpy.degrees(math.atan(old_div(1.0, baldwin_lower_alt_bound)))
                # print(  "  baldwin_lower_alt_bound ",phi,theta,baldwin_lower_alt_bound,self.brackets[inc_mirror])

                if (baldwin_lower_alt_bound > theta):
                    if (inc_mirror == 1):
                        return True
                    else:
                        baldwin_upper_alt_bound = \
                          old_div (
                              (old_div((numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) - tmphi ))),
                              (numpy.tan(numpy.radians(self.brackets[inc_mirror][1])) ))
                             + old_div((numpy.sin(numpy.radians(tmphi))) ,
                             numpy.tan(numpy.radians(old_div(self.brackets[inc_mirror][3], 2.0)   )))) \
                            , (numpy.sin(numpy.radians(old_div(self.brackets[inc_mirror][2], 2.0) ))))
                        baldwin_upper_alt_bound = numpy.degrees(math.atan(old_div(1.0 , baldwin_upper_alt_bound)))
                        # print(  "  baldwin_upper_alt_bound ",phi,theta,baldwin_upper_alt_bound,self.brackets[inc_mirror])

                        if (baldwin_upper_alt_bound < theta):
                            return False
                        else:
                            return True
                elif theta == 0:
                    return True
                else:
                    return False
            elif theta == 0:
                return True
            else:
                # print "phi",self.brackets
                return False
        else:
            global_def.ERROR("unknown symmetry", "symclass: is_in_subunit", 1)

    def symmetry_related(self, angles):
        """
        Generate all symmetry related instances of input angles
        Input:  is a list of list of n triplets
            [[phi0,theta0,psi0],[phi1,theta1,psi1],...]
        Output: is a list of list triplets whose length is k times larger than that of input,
                where k is symmetry multiplicity.  First k are symmetry related version of the first triplet,
                second k are symmetry related version of the second triplet and so on...
           [[phi0,theta0,0],[phi0,theta0,0]_SYM0,...,[phi1,theta1,],[phi1,theta1,]_SYM1,...]
        """
        redang = [angles[:]]
        if (self.sym[0] == "c"):
            qt = old_div(360.0, self.nsym)
            if angles[1] != 0 and angles[1] != 180:
                for l in range(1, self.nsym):
                    redang.append([(angles[0] + l * qt) % 360.0, angles[1], angles[2]])
        elif (self.sym[0] == "d"):
            nsm = old_div(self.nsym,  2)
            qt = old_div(360.0,  nsm)
            if angles[1] == 0:
                redang.append([0, 180, (angles[2] + 180.0 * (nsm % 2)) % 360.0])
            else:
                for l in range(1, nsm):
                    redang.append([(angles[0] + l * qt) % 360.0, angles[1], angles[2]])
                if angles[1] != 90:
                    for l in range(nsm, self.nsym):
                        redang.append([(360.0 - redang[l - nsm][0]) % 360.0, 180.0 - angles[1],
                                       (angles[2] + 180.0 * (nsm % 2)) % 360.0])

        else:
            pass  # IMPORTIMPORTIMPORT from fundamentals import rotmatrix, recmat, mulmat
            mat = self.rotmatrix(angles[0], angles[1], angles[2])
            for l in range(1, self.nsym):
                p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                redang.append([p1, p2, p3])

        return redang

    def symmetry_neighbors(self, angles):
        """
        Generate symmetry related instances of input angles only for asymmetric regions adjacent to the zero's one
          input is a list of lists  [[phi0,theta0,0],[phi1,theta1,0],...]
        Output: is a list of list triplets whose length is k times larger than that of input,
                where k is the number of adjacent regions for a given point group symmetry.
                First k are symmetry related version of the first triplet,
                second k are symmetry related version of the second triplet and so on...
          output is [[phi0,theta0,0],[phi0,theta0,0]_SYM1,...,[phi1,theta1,],[phi1,theta1,]_SYM1,...]
        """
        if (self.sym[0] == "c" or self.sym[0] == "d"):
            temp = EMAN2_cppwrap.Util.symmetry_neighbors(angles, self.sym)
            nt = old_div(len(temp), 3)
            return [[temp[l * 3], temp[l * 3 + 1], 0.0] for l in range(nt)]
        #  Note symmetry neighbors below refer to the particular order
        #   in which this class generates symmetry matrices
        neighbors = {}
        neighbors["oct"] = [1, 2, 3, 8, 9, 12, 13]
        neighbors["tet"] = [1, 2, 3, 4, 6, 7]
        neighbors["icos"] = [1, 2, 3, 4, 6, 7, 11, 12]
        sang = [[] for l in range(len(angles) * (len(neighbors[self.sym]) + 1))]
        for i, q in enumerate(angles):  sang[i * (len(neighbors[self.sym]) + 1)] = angles[i][:]
        pass  # IMPORTIMPORTIMPORT from fundamentals import rotmatrix, recmat, mulmat
        for i, q in enumerate(angles):
            mat = self.rotmatrix(q[0], q[1], q[2])
            for j, l in enumerate(neighbors[self.sym]):
                p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                sang[i * (len(neighbors[self.sym]) + 1) + (j + 1)] = [p1, p2, 0.0]
        return sang

    def reduce_anglesets(self, angles, inc_mirror=1):
        """
          Input is either list ot lists [[phi,thet,psi],[],[]] or a triplet [phi,thet,psi]
                inc_mirror = 1 consider mirror directions as unique
                inc_mirror = 0 consider mirror directions as outside of unique range.
          It will map all triplets to the first asymmetric subunit.
        """
        pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
        pass  # IMPORTIMPORTIMPORT import types
        is_platonic_sym = self.sym[0] == "o" or self.sym[0] == "i"
        if (self.sym[0] == "c"):
            qs = old_div(360.0, self.nsym)
        elif (self.sym[0] == "d"):
            qs = old_div(720.0, self.nsym)
        if (type(angles[0]) is list):
            toprocess = angles
            lis = True
        else:
            toprocess = [angles]
            lis = False
        redang = []
        for q in toprocess:
            phi = q[0];
            theta = q[1];
            psi = q[2]
            if is_platonic_sym:
                if (not self.is_in_subunit(phi, theta, 1)):
                    mat = self.rotmatrix(phi, theta, psi)
                    for l in range(self.nsym):
                        p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                        # print(p1,p2,p3)
                        if (self.is_in_subunit(p1, p2, 1)):
                            phi = p1;
                            theta = p2;
                            psi = p3
                            # print("  FOUND ")
                            break
                    if (inc_mirror == 0):
                        if (phi >= self.brackets[0][0]):
                            phi = self.brackets[1][0] - phi
                            if (l > 0): psi = (360.0 - psi) % 360.0
                elif (inc_mirror == 0):
                    if (phi >= self.brackets[0][0]):
                        phi = self.brackets[1][0] - phi
                        psi = (360.0 - psi) % 360.0
                    elif (inc_mirror == 0):
                        if (phi >= self.brackets[0][0]):  phi = self.brackets[1][0] - phi
            elif (self.sym[0] == "t"):
                if (not self.is_in_subunit(phi, theta, inc_mirror)):
                    mat = self.rotmatrix(phi, theta, psi)
                    fifi = False
                    for l in range(self.nsym):
                        p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                        if (self.is_in_subunit(p1, p2, inc_mirror)):
                            if (inc_mirror):
                                phi = p1;
                                theta = p2;
                                psi = p3
                                fifi = True
                                break
                            else:
                                if (self.is_in_subunit(p1, p2, 0)):
                                    phi = p1;
                                    theta = p2;
                                    psi = p3
                                    fifi = True
                                    break
                    if (inc_mirror == 1 and not fifi): print("  FAILED no mirror ")
                    if (not fifi):
                        phi = (180.0 + phi) % 360.0;
                        theta = 180.0 - theta;
                        psi = (180.0 - psi) % 360.0
                        mat = self.rotmatrix(phi, theta, psi)
                        for l in range(self.nsym):
                            p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                            if (self.is_in_subunit(p1, p2, 0)):
                                phi = p1;
                                theta = p2;
                                psi = p3
                                fifi = True
                                break

                    if (not fifi):  print("  FAILED mirror ")
            else:
                if (theta > 90.0 and inc_mirror == 0):
                    phi = (180.0 + phi) % 360.0;
                    theta = 180.0 - theta;
                    psi = (180.0 - psi) % 360.0
                phi = phi % qs
                if (self.sym[0] == "d"):
                    if (inc_mirror == 0):
                        if (old_div(self.nsym, 2) % 2 == 0):
                            if (phi >= (old_div(qs, 2))):
                                phi = qs - phi
                                psi = (360.0 - psi) % 360.0
                        else:
                            if ( phi >=  (old_div(old_div(360.0, self.nsym ),2))  and phi < old_div(360.0, self.nsym)):
                                phi = old_div(360.0, self.nsym) - phi
                                psi = 360.0 - psi
                            elif (phi >=  old_div(360.0, self.nsym) + old_div(old_div(360.0, self.nsym ),2)  and phi < old_div(720.0, self.nsym)):
                                phi = old_div(720.0, self.nsym) - phi + old_div(360.0, self.nsym)
                                psi = (360.0 - psi) % 360.0


            redang.append([phi, theta, psi])

        if lis:
            return redang
        else:
            return redang[0]

    def reduce_angles(self, phiin, thetain, psiin, inc_mirror=1):
        """
          It will map three input angles to the first asymmetric subunit.
                inc_mirror = 1 consider mirror directions as unique
                inc_mirror = 0 consider mirror directions as outside of unique range.
        """
        pass  # IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
        is_platonic_sym = self.sym[0] == "o" or self.sym[0] == "i"
        if (self.sym[0] == "c"):
            qs = old_div(360.0, self.nsym)
        elif (self.sym[0] == "d"):
            qs = old_div(720.0, self.nsym)
        if is_platonic_sym:
            if (not self.is_in_subunit(phiin, thetain, 1)):
                mat = self.rotmatrix(phiin, thetain, psiin)
                phi = phiin;
                theta = thetain;
                psi = psiin
                for l in range(self.nsym):
                    p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                    # print(p1,p2,p3)
                    if (self.is_in_subunit(p1, p2, 1)):
                        phi = p1;
                        theta = p2;
                        psi = p3
                        # print("  FOUND ")
                        break
                if (inc_mirror == 0):
                    if (phi >= self.brackets[0][0]):
                        phi = self.brackets[1][0] - phi
                        if (l > 0): psi = (360.0 - psi) % 360.0
            elif (inc_mirror == 0):
                phi = phiin;
                theta = thetain;
                psi = psiin
                if (phi >= self.brackets[0][0]):
                    phi = self.brackets[1][0] - phi
                    psi = (360.0 - psi) % 360.0
            else:
                phi = phiin;
                theta = thetain;
                psi = psiin
        elif (self.sym[0] == "t"):
            phi = phiin;
            theta = thetain;
            psi = psiin
            if (not self.is_in_subunit(phi, theta, inc_mirror)):
                mat = self.rotmatrix(phi, theta, psi)
                for l in range(self.nsym):
                    p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                    # print(" ASYM  ",l,p1,p2,p3 )
                    if (self.is_in_subunit(p1, p2, inc_mirror)):
                        if (inc_mirror):
                            return p1, p2, p3
                        else:
                            if (self.is_in_subunit(p1, p2, 0)):  return p1, p2, p3
                if (inc_mirror == 1): print("  FAILED no mirror ")
                phi = (180.0 + phi) % 360.0;
                theta = 180.0 - theta;
                psi = (180.0 - psi) % 360.0
                mat = self.rotmatrix(phi, theta, psi)
                for l in range(self.nsym):
                    p1, p2, p3 = self.recmat(self.mulmat(mat, self.symatrix[l]))
                    # print(" MIR  ",l,p1,p2,p3 )
                    if (self.is_in_subunit(p1, p2, 0)):  return p1, p2, p3

                print("  FAILED mirror ")
        else:
            if (thetain > 90.0 and inc_mirror == 0):
                phi = (180.0 + phiin) % 360.0;
                theta = 180.0 - thetain;
                psi = (180.0 - psiin) % 360.0
            else:
                phi = phiin;
                theta = thetain;
                psi = psiin
            phi = phi % qs
            if (self.sym[0] == "d"):
                if (inc_mirror == 0):
                    if (old_div(self.nsym, 2) % 2 == 0):
                        if (phi >= (old_div(qs, 2))):
                            phi = qs - phi
                            psi = 360.0 - psi
                    else:
                        if (phi >=  (old_div(old_div(360.0, self.nsym ),2))  and phi < (old_div(360.0, self.nsym)) ):
                            phi = old_div(360.0, self.nsym) - phi
                            psi = 360.0 - psi
                        elif (phi >= old_div(360.0, self.nsym ) + old_div(old_div(360.0, self.nsym ),2) and phi < old_div(720.0, self.nsym)):
                            phi = old_div(720.0, self.nsym) - phi + old_div(360.0, self.nsym)
                            psi = 360.0 - psi

        return phi, theta, psi

    """Multiline Comment5"""

    def even_angles(self, delta=15.0, theta1=-1.0, theta2=-1.0, phi1=-1.0, phi2=-1.0, \
                    method='S', phiEqpsi="Zero", inc_mirror=1):
        """Create a list of Euler angles suitable for projections.
           method is either 'S' - for Saff algorithm
                       or   'P' - for Penczek '94 algorithm
                         'S' assumes phi1<phi2 and phi2-phi1>> delta ;
           symmetry  - if this is set to point-group symmetry (cn or dn) or helical symmetry with point-group symmetry (scn or sdn), \
                         it will yield angles from the asymmetric unit, not the specified range;
        """

        pass  # IMPORTIMPORTIMPORT from math      import pi, sqrt, cos, acos, tan, sin, radians, degrees
        pass  # IMPORTIMPORTIMPORT from utilities import even_angles_cd
        angles = []
        phi2_org = phi2
        if (phi2_org < 0.0):  phi2_org = self.brackets[1][0] - 1.0e-7  # exclude right border of unit
        theta2_org = theta2
        if (theta2_org < 0.0): theta2_org = self.brackets[1][3]
        if (phi2 < phi1 or theta2 < theta1 or delta <= 0.0):  global_def.ERROR("even_angles",
                                                                                     "incorrect parameters (phi1,phi2,theta1,theta2,delta): %f   %f   %f   %f   %f" % (
                                                                                     phi1, phi2, theta1, theta2, delta),
                                                                                     1)
        if (phi1 < 0.0):  phi1 = 0.0
        if (phi2 < 0.0):  phi2 = self.brackets[inc_mirror][0] - 1.0e-7  # exclude right border of unit
        if (theta1 < 0.0): theta1 = 0.0
        if (theta2 < 0.0): theta2 = self.brackets[inc_mirror][3]
        # print " parameters (phi1,phi2,theta,theta2,delta): %f   %f   %f   %f   %f"%(phi1,phi2,theta1,theta2,delta)
        #
        if (self.sym[0] != "s"):
            """Multiline Comment6"""
            angles = []
            zzis_platonic_sym = self.sym[0] == "o" or self.sym[0] == "t" or self.sym[0] == "i"
            if (method == 'P'):
                theta = theta1
                while (theta <= theta2):
                    phi = phi1
                    if (theta == 0.0 or theta == 180.0):
                        detphi = 2 * phi2
                    else:
                        detphi = old_div(delta , numpy.sin(numpy.radians(theta)))
                    while (phi < phi2):
                        if (self.is_in_subunit(phi, theta, inc_mirror)):
                            angles.append([phi, theta, 0.0])
                        else:
                            angles.append([phi, theta, 0.0])
                        phi += detphi
                    theta += delta
            else:
                # I have to use original phi2 and theta2 to compute Deltaz and wedgeFactor as otherwise
                # points for include mirror differ from do not include mirror.
                Deltaz = numpy.cos(numpy.radians(theta2_org)) - numpy.cos(numpy.radians(theta1))
                s = delta * old_div(numpy.pi, 180.0)
                NFactor = old_div(3.6,  s)
                wedgeFactor = abs(Deltaz * old_div((phi2_org - phi1), 720.0))
                NumPoints = int(NFactor * NFactor * wedgeFactor)
                angles.append([phi1, theta1, 0.0])
                # initialize loop
                phistep = phi2_org - phi1
                z1 = numpy.cos(numpy.radians(theta1))
                phi = phi1
                for k in range(1, NumPoints - 1):
                    z = z1 + old_div(Deltaz * k, (NumPoints - 1))
                    r = numpy.sqrt(1.0 - z * z)
                    phi = phi1 + (phi +  old_div(delta, r) - phi1)%phistep
                    theta = numpy.degrees(math.acos(z))
                    if (theta > 180.0):  break
                    if (not self.is_in_subunit(phi, theta, inc_mirror)): continue
                    angles.append([phi, theta, 0.0])
            # angles.append([p2,t2,0])  # This is incorrect, as the last angle is really the border, not the element we need. PAP 01/15/07
            if (phiEqpsi == 'Minus'):
                for k in range(len(angles)): angles[k][2] = (720.0 - angles[k][0]) % 360.0
            if ((self.sym[0] == "c" or self.sym[0] == "d") and (
                    (theta2 == 180.) or (theta2 >= 180. and delta == 180.0))):  angles.append([0.0, 180.0, 0.0])

        """Multiline Comment7"""
        return angles
    @staticmethod
    def rotmatrix(phi, theta, psi):
        pass  # IMPORTIMPORTIMPORT from math import sin,cos,radians
        rphi = numpy.radians(phi)
        rtheta = numpy.radians(theta)
        rpsi = numpy.radians(psi)
        cosphi = numpy.cos(rphi)
        sinphi = numpy.sin(rphi)
        costheta = numpy.cos(rtheta)
        sintheta = numpy.sin(rtheta)
        cospsi = numpy.cos(rpsi)
        sinpsi = numpy.sin(rpsi)
        mat = [[0.0] * 3, [0.0] * 3, [0.0] * 3]

        mat[0][0] = cospsi * costheta * cosphi - sinpsi * sinphi
        mat[1][0] = -sinpsi * costheta * cosphi - cospsi * sinphi
        mat[2][0] = sintheta * cosphi

        mat[0][1] = cospsi * costheta * sinphi + sinpsi * cosphi
        mat[1][1] = -sinpsi * costheta * sinphi + cospsi * cosphi
        mat[2][1] = sintheta * sinphi

        mat[0][2] = -cospsi * sintheta
        mat[1][2] = sinpsi * sintheta
        mat[2][2] = costheta
        return mat

    @staticmethod
    def mulmat(m1, m2):
        mat = [[0.0] * 3, [0.0] * 3, [0.0] * 3]
        """
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    mat[i][j] += m1[i][k]*m2[k][j]
                #mat[i][j] = round(mat[i][j],8)
        """
        mat[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0]
        mat[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1]
        mat[0][2] = m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2]
        mat[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0]
        mat[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1]
        mat[1][2] = m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2]
        mat[2][0] = m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0]
        mat[2][1] = m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1]
        mat[2][2] = m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2]

        return mat

    @staticmethod
    def recmat(mat):
        pass  # IMPORTIMPORTIMPORT from math import acos,asin,atan2,degrees,pi

        def sign(x):
            if (x >= 0.0):
                return 1
            else:
                return -1

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
        if (mat[2][2] == 1.0):
            theta = 0.0
            psi = 0.0
            if (mat[0][0] == 0.0):
                phi = math.asin(mat[0][1])
            else:
                phi = math.atan2(mat[0][1], mat[0][0])
        elif (mat[2][2] == -1.0):
            theta = numpy.pi
            psi = 0.0
            if (mat[0][0] == 0.0):
                phi = math.asin(-mat[0][1])
            else:
                phi = math.atan2(-mat[0][1], -mat[0][0])
        else:
            theta = math.acos(mat[2][2])
            st = sign(theta)
            # print theta,st,mat[2][0]
            if (mat[2][0] == 0.0):
                if (st != sign(mat[2][1])):
                    phi = 1.5 * numpy.pi
                else:
                    phi = 0.5 * numpy.pi
            else:
                phi = math.atan2(st * mat[2][1], st * mat[2][0])

            # print theta,st,mat[0][2],mat[1][2]
            if (mat[0][2] == 0.0):
                if (st != sign(mat[1][2])):
                    psi = 1.5 * numpy.pi
                else:
                    psi = 0.5 * numpy.pi
            else:
                psi = math.atan2(st * mat[1][2], -st * mat[0][2])
        # pi2 = 2*pi
        # return  degrees(round(phi%pi2,8)),degrees(round(theta%pi2,8)),degrees(round(psi%pi2,8))
        # return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
        return numpy.degrees(phi) % 360.0, numpy.degrees(theta) % 360.0, numpy.degrees(psi) % 360.0


if __name__ == '__main__':
    unittest.main()



