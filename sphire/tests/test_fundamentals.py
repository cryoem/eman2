from __future__ import print_function
from __future__ import division

import nose.tools as nt

import EMAN2_cppwrap
import numpy
import math
from ..libpy import sparx_fundamentals as fu


########
########
########


class TestSymClassInitUnknown:

    @nt.raises(AttributeError)
    def test_unknown_sym_should_return_error(self):
        symclass = fu.symclass('foobar')


########
########
########


class TestSymClassInitC:

    @nt.raises(ZeroDivisionError)
    def test_c0_sym_should_return_error(self):
        symclass = fu.symclass('c0')

    def test_lower_c1_sym_should_return_lower_c1_sym(self):
        symclass = fu.symclass('c1')
        print(symclass.sym)
        assert symclass.sym == 'c1'

    def test_upper_c1_sym_should_return_lower_c1_sym(self):
        symclass = fu.symclass('C1')
        print(symclass.sym)
        assert symclass.sym == 'c1'

    def test_c1_sym_should_return_nsym_1(self):
        symclass = fu.symclass('c1')
        print(symclass.nsym)
        assert symclass.nsym == 1

    def test_c1_sym_should_return_nsym_5(self):
        symclass = fu.symclass('c5')
        print(symclass.nsym)
        assert symclass.nsym == 5

    def test_c1_should_return_correct_brackets(self):
        symclass = fu.symclass('c1')
        print(symclass.brackets)
        nsym = 1
        assert symclass.brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

    def test_c5_should_return_correct_brackets(self):
        symclass = fu.symclass('c5')
        print(symclass.brackets)
        nsym = 5
        assert symclass.brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [360./nsym, 180.0, 360./nsym, 180.0]
            ]

    def test_c1_should_return_correct_symangles(self):
        symclass = fu.symclass('c1')
        print(symclass.symangles)
        nsym = 1
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i*360./nsym])
        assert symclass.symangles == symangles

    def test_c5_should_return_correct_symangles(self):
        symclass = fu.symclass('c5')
        print(symclass.symangles)
        nsym = 5
        symangles = []
        for i in range(nsym):
            symangles.append([0.0, 0.0, i*360./nsym])
        assert symclass.symangles == symangles

    def test_c1_should_return_correct_transform(self):
        symclass = fu.symclass('c1')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_c5_should_return_correct_transform(self):
        symclass = fu.symclass('c5')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_c1_should_return_correct_symatrix(self):
        symclass = fu.symclass('c1')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

    def test_c5_should_return_correct_symatrix(self):
        symclass = fu.symclass('c5')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

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


########
########
########


class TestSymClassInitD:

    @nt.raises(ZeroDivisionError)
    def test_d0_sym_should_return_error(self):
        symclass = fu.symclass('d0')

    def test_lower_d1_sym_should_return_lower_d1_sym(self):
        symclass = fu.symclass('d1')
        print(symclass.sym)
        assert symclass.sym == 'd1'

    def test_upper_d1_sym_should_return_lower_d1_sym(self):
        symclass = fu.symclass('D1')
        print(symclass.sym)
        assert symclass.sym == 'd1'

    def test_d1_sym_should_return_nsym_2(self):
        symclass = fu.symclass('d1')
        print(symclass.nsym)
        assert symclass.nsym == 2

    def test_d1_sym_should_return_nsym_10(self):
        symclass = fu.symclass('d5')
        print(symclass.nsym)
        assert symclass.nsym == 10

    def test_d1_should_return_correct_brackets(self):
        symclass = fu.symclass('d1')
        print(symclass.brackets)
        nsym = 2*1
        assert symclass.brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [2*360./nsym, 90.0, 2*360./nsym, 90.0]
            ]

    def test_d5_should_return_correct_brackets(self):
        symclass = fu.symclass('d5')
        print(symclass.brackets)
        nsym = 2*5
        assert symclass.brackets == [
            [360./nsym, 90.0, 360./nsym, 90.0],
            [2*360./nsym, 90.0, 2*360./nsym, 90.0]
            ]

    def test_d1_should_return_correct_symangles(self):
        symclass = fu.symclass('d1')
        print(symclass.symangles)
        nsym = 2*1
        symangles = []
        for i in range(nsym//2):
            symangles.append([0.0, 0.0, i*360./nsym*2])
        for i in range(nsym//2):
            symangles.append([0.0, 180.0, (i*360./nsym*2+180.0*(int(symclass.sym[1:])%2))%360.0])
        assert symclass.symangles == symangles

    def test_d5_should_return_correct_symangles(self):
        symclass = fu.symclass('d5')
        print(symclass.symangles)
        nsym = 2*5
        symangles = []
        for i in range(nsym//2):
            symangles.append([0.0, 0.0, i*360./nsym*2])
        for i in range(nsym//2):
            symangles.append([0.0, 180.0, (i*360./nsym*2+180.0*(int(symclass.sym[1:])%2))%360.0])
        assert symclass.symangles == symangles

    def test_d1_should_return_correct_transform(self):
        symclass = fu.symclass('d1')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_d5_should_return_correct_transform(self):
        symclass = fu.symclass('d5')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_d1_should_return_correct_symatrix(self):
        symclass = fu.symclass('d1')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

    def test_d5_should_return_correct_symatrix(self):
        symclass = fu.symclass('d5')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

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


########
########
########


class TestSymClassInitOct:

    def test_lower_oct_sym_should_return_lower_oct_sym(self):
        symclass = fu.symclass('oct')
        print(symclass.sym)
        assert symclass.sym == 'oct'

    def test_upper_oct_sym_should_return_lower_oct_sym(self):
        symclass = fu.symclass('OCT')
        print(symclass.sym)
        assert symclass.sym == 'oct'

    def test_oct_sym_should_return_nsym_24(self):
        symclass = fu.symclass('OCT')
        print(symclass.nsym)
        assert symclass.nsym == 24

    def test_oct_should_return_correct_brackets(self):
        symclass = fu.symclass('oct')
        print(symclass.brackets)
        nsym = 24
        brackets = []
        ncap = 4
        cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
        alpha = numpy.degrees(math.acos(1.0/(numpy.sqrt(3.0)*numpy.tan(2*numpy.pi/ncap/2.0)))) # also platonic_params["alt_max"]
        theta = numpy.degrees(0.5*math.acos( numpy.cos(numpy.radians(cap_sig))/(1.0-numpy.cos(numpy.radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
        brackets = [[180./ncap,theta,cap_sig,alpha],[360./ncap,theta,cap_sig,alpha]]
        assert symclass.brackets == brackets

    def test_oct_should_return_correct_symangles(self):
        symclass = fu.symclass('oct')
        print(symclass.symangles)
        nsym = 24
        brackets = []
        ncap = 4
        cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
        alpha = numpy.degrees(math.acos(1.0/(numpy.sqrt(3.0)*numpy.tan(2*numpy.pi/ncap/2.0)))) # also platonic_params["alt_max"]
        theta = numpy.degrees(0.5*math.acos( numpy.cos(numpy.radians(cap_sig))/(1.0-numpy.cos(numpy.radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
        symangles = [[0.0,0.0,float(i)] for i in range(0,271,90)]
        for i in range(0,271,90):
            for j in range(0,271,90):
                symangles.append([float(j),90.0,float(i)])
        for i in range(0,271,90):  symangles.append([0.0,180.0,float(i)])
        assert symclass.symangles == symangles

    def test_oct_should_return_correct_transform(self):
        symclass = fu.symclass('oct')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_oct_should_return_correct_symatrix(self):
        symclass = fu.symclass('oct')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

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


########
########
########


class TestSymClassInitTet:

    def test_lower_tet_sym_should_return_lower_tet_sym(self):
        symclass = fu.symclass('tet')
        print(symclass.sym)
        assert symclass.sym == 'tet'

    def test_upper_tet_sym_should_return_lower_tet_sym(self):
        symclass = fu.symclass('tet')
        print(symclass.sym)
        assert symclass.sym == 'tet'

    def test_tet_sym_should_return_nsym_12(self):
        symclass = fu.symclass('tet')
        print(symclass.nsym)
        assert symclass.nsym == 12

    def test_tet_should_return_correct_brackets(self):
        symclass = fu.symclass('tet')
        print(symclass.brackets)
        nsym = 12
        brackets = []
        ncap = 3
        cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
        alpha = numpy.degrees(math.acos(1.0/(numpy.sqrt(3.0)*numpy.tan(2*numpy.pi/ncap/2.0)))) # also platonic_params["alt_max"]
        theta = numpy.degrees(0.5*math.acos( numpy.cos(numpy.radians(cap_sig))/(1.0-numpy.cos(numpy.radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
        brackets = [[360.0/ncap,theta,cap_sig,alpha],[360.0/ncap,theta,cap_sig,alpha]]
        assert symclass.brackets == brackets

    def test_tet_should_return_correct_symangles(self):
        symclass = fu.symclass('tet')
        print(symclass.symangles)
        nsym = 12
        symangles = []
        lvl1 = numpy.degrees(math.acos(-1.0/3.0)) # There  are 3 faces at this angle
        symangles = [ [0.,0.,0.], [0., 0., 120.], [0., 0., 240.]]
        for l1 in range(0,241,120):
            for l2 in range(60,301,120):
                symangles.append([float(l1),lvl1,float(l2)])
        assert symclass.symangles == symangles

    def test_tet_should_return_correct_transform(self):
        symclass = fu.symclass('tet')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_tet_should_return_correct_symatrix(self):
        symclass = fu.symclass('tet')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

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


########
########
########


class TestSymClassInitIcos:

    def test_lower_icos_sym_should_return_lower_icos_sym(self):
        symclass = fu.symclass('icos')
        print(symclass.sym)
        assert symclass.sym == 'icos'

    def test_upper_icos_sym_should_return_lower_icos_sym(self):
        symclass = fu.symclass('icos')
        print(symclass.sym)
        assert symclass.sym == 'icos'

    def test_icos_sym_should_return_nsym_60(self):
        symclass = fu.symclass('icos')
        print(symclass.nsym)
        assert symclass.nsym == 60

    def test_icos_should_return_correct_brackets(self):
        symclass = fu.symclass('icos')
        print(symclass.brackets)
        nsym = 60
        brackets = []
        ncap = 5
        cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
        alpha = numpy.degrees(math.acos(1.0/(numpy.sqrt(3.0)*numpy.tan(2*numpy.pi/ncap/2.0)))) # also platonic_params["alt_max"]
        theta = numpy.degrees(0.5*math.acos( numpy.cos(numpy.radians(cap_sig))/(1.0-numpy.cos(numpy.radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
        brackets = [[36.,theta,cap_sig,alpha],[72.,theta,cap_sig,alpha]]
        assert symclass.brackets == brackets

    def test_icos_should_return_correct_symangles(self):
        symclass = fu.symclass('icos')
        print(symclass.symangles)
        nsym = 60
        symangles = []
        lvl1= numpy.degrees(math.atan(2.0))  #there are 5 pentagons with centers at this height (angle)
        lvl2 = 180.0 - lvl1      #there are 5 pentagons with centers at this height (angle)
        symangles = [[0.0,0.0,float(i)] for i in range(0,288+1,72)]
        for l1 in range(0,288+1,72):
            for l2 in range(36,324+1,72):
                symangles.append([float(l1),lvl1,float(l2)])
        for l1 in range(36,324+1,72):
            for l2 in range(0,288+1,72):
                symangles.append([float(l1),lvl2,float(l2)])
        for i in range(0,288+1,72):  symangles.append([0.0,180.0,float(i)])
        assert symclass.symangles == symangles

    def test_icos_should_return_correct_transform(self):
        symclass = fu.symclass('icos')
        print(symclass.transform)
        transform = []
        for args in symclass.symangles:
            transform.append(EMAN2_cppwrap.Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
        assert symclass.transform == transform

    def test_icos_should_return_correct_symatrix(self):
        symclass = fu.symclass('icos')
        print(symclass.symatrix)
        symatrix = []
        for args in symclass.symangles:
            symatrix.append(self.rotmatrix(args[0],args[1],args[2]))
        assert symclass.symatrix == symatrix

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


########
########
########


class TestSymClassIsInSubunitC:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_c1_sym_no_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        symclass = fu.symclass('c1')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c1_sym_no_mirror_theta_larger_90_degrees_should_return_False(self):
        symclass = fu.symclass('c1')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c1_sym_mirror_theta_smaller_equals_90_degrees_should_return_True(self):
        symclass = fu.symclass('c1')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c1_sym_mirror_theta_larger_90_degrees_should_return_False(self):
        symclass = fu.symclass('c1')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_no_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_no_mirror_theta_larger_90_degrees_phi_bigger_equals_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_mirror_theta_smaller_equals_90_phi_smaller_72_degrees_should_return_True(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_mirror_theta_larger_90_degrees_phi_smaller_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_mirror_theta_smaller_equals_90_degrees_phi_bigger_equals_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_c5_sym_mirror_theta_larger_90_degrees_phi_bigger_equals_72_should_return_False(self):
        symclass = fu.symclass('c5')
        angles = [[entry, thet, 0] for entry in range(72, 360) for thet in range(90, 180)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results


########
########
########


class TestSymClassIsInSubunitDEven:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'
    

    def test_d4_sym_no_mirror_theta_smaller_equals_90_degrees_phi_smaller_45_should_return_True(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(int(360/8)) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d4_sym_no_mirror_theta_smaller_equals_90_degrees_phi_larger_45_smaller_90_should_return_False(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(int(360/8), int(360/4)) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d4_sym_no_mirror_theta_larger_90_degrees_should_return_False(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d4_sym_mirror_theta_smaller_equals_90_degrees_phi_smaller_45_should_return_True(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(int(360/8)) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d4_sym_mirror_theta_smaller_equals_90_degrees_phi_larger_45_smaller_90_should_return_True(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(int(360/8), int(360/4)) for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d4_sym_mirror_theta_larger_90_degrees_should_return_False(self):
        symclass = fu.symclass('d4')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results


########
########
########


class TestSymClassIsInSubunitDOdd:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'
    

    def test_d5_phi_ge_18_le_54_should_return_True(self):
        symclass = fu.symclass('d5')
        nsym = 2 * 5
        hnsym = 2 * nsym
        lower = int(360/nsym) - int(360/hnsym)
        higher = int(360/nsym) + int(360/hnsym) + 1
        angles = [[entry, thet, 0] for entry in range(lower, higher) for thet in range(0, 91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        print(symclass.brackets)
        assert results == expected_results

    def test_d5_phi_le_18_ge_54_l_72_should_return_False(self):
        symclass = fu.symclass('d5')
        nsym = 2 * 5
        hnsym = 2 * nsym
        lower = int(360/nsym) - int(360/hnsym)
        higher = int(360/nsym) + int(360/hnsym) + 1
        exclude = set(range(lower, higher))
        angles = [[entry, thet, 0] for entry in range(int(360/5)) if entry not in exclude for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d5_phi_ge_72_should_return_False(self):
        symclass = fu.symclass('d5')
        nsym = 2 * 5
        hnsym = 2 * nsym
        lower = int(360/nsym) - int(360/hnsym)
        higher = int(360/nsym) + int(360/hnsym) + 1
        exclude = set(range(lower, higher))
        angles = [[entry, thet, 0] for entry in range(int(360/5), 360) if entry not in exclude for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d5_theta_g_90_should_return_False(self):
        symclass = fu.symclass('d5')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 0)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d5_mirror_phi_l_72_should_return_True(self):
        symclass = fu.symclass('d5')
        angles = [[entry, thet, 0] for entry in range(int(360/5)) for thet in range(91)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [True] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d5_mirror_theta_g_90_should_return_False(self):
        symclass = fu.symclass('d5')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(90, 180)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        results = []
        for phi, theta, psi in angles:
            if theta != 180:
                theta = theta + 0.1
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results

    def test_d5_mirror_phi_g_72_should_return_False(self):
        symclass = fu.symclass('d5')
        nsym = 2 * 5
        hnsym = 2 * nsym
        lower = int(360/nsym) - int(360/hnsym)
        higher = int(360/nsym) + int(360/hnsym) + 1
        exclude = set(range(lower, higher))
        angles = [[entry, thet, 0] for entry in range(int(360/5), 360) if entry not in exclude for thet in range(1, 91)]
        results = []
        for phi, theta, psi in angles:
            is_in = symclass.is_in_subunit(phi, theta, 1)
            results.append(is_in)
        expected_results = [False] * len(angles)
        for i, k, j in zip(results, expected_results, angles):
            if i != k:
                print(self.output_template_angles.format(i, k, j))
        assert results == expected_results


########
########
########


class TestSymClassIsInSubunitOct:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    @staticmethod
    def generate_result_oct(phi, theta, inc_mirror, brackets):
        if phi >= 0.0 and phi < brackets[inc_mirror][0] and theta <= brackets[inc_mirror][3]:
            tmphi = min(phi, brackets[inc_mirror][2] - phi)
            baldwin_lower_alt_bound = \
                (numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0-tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][1])) + \
                numpy.sin(numpy.radians(tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][3])))/numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0))
            baldwin_lower_alt_bound = numpy.degrees(math.atan(1.0/baldwin_lower_alt_bound))
            if baldwin_lower_alt_bound > theta:
                return True
            else:
                return False
        elif theta == 0:
            return True
        else:
            return False

    def test_oct_no_mirror(self):
        symclass = fu.symclass('oct')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_oct(phi, theta, 0, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 0))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results

    def test_oct_mirror(self):
        symclass = fu.symclass('oct')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_oct(phi, theta, 1, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 1))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results


########
########
########


class TestSymClassIsInSubunitIcos:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    @staticmethod
    def generate_result_icos(phi, theta, inc_mirror, brackets):
        if phi >= 0.0 and phi < brackets[inc_mirror][0] and theta <= brackets[inc_mirror][3]:
            tmphi = min(phi, brackets[inc_mirror][2] - phi)
            baldwin_lower_alt_bound = \
                (numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0-tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][1])) + \
                numpy.sin(numpy.radians(tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][3])))/numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0))
            baldwin_lower_alt_bound = numpy.degrees(math.atan(1.0/baldwin_lower_alt_bound))
            if baldwin_lower_alt_bound > theta:
                return True
            else:
                return False
        elif theta == 0:
            return True
        else:
            return False

    def test_icos_no_mirror(self):
        symclass = fu.symclass('icos')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_icos(phi, theta, 0, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 0))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results

    def test_icos_mirror(self):
        symclass = fu.symclass('icos')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_icos(phi, theta, 1, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 1))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results


########
########
########


class TestSymClassIsInSubunitTet:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    @staticmethod
    def generate_result_tet(phi, theta, inc_mirror, brackets):
        if phi >= 0.0 and phi < brackets[inc_mirror][0] and theta <= brackets[inc_mirror][3]:
            tmphi = min(phi, brackets[inc_mirror][2]-phi)
            baldwin_lower_alt_bound = \
                (numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0-tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][1])) + \
                numpy.sin(numpy.radians(tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][3])))/numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0))
            baldwin_lower_alt_bound = numpy.degrees(math.atan(1.0/baldwin_lower_alt_bound))
            if baldwin_lower_alt_bound > theta:
                if inc_mirror == 1:
                    return True
                else:
                    baldwin_upper_alt_bound = \
                        (numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0-tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][1])) + \
                        numpy.sin(numpy.radians(tmphi))/numpy.tan(numpy.radians(brackets[inc_mirror][3]/2.0)))/numpy.sin(numpy.radians(brackets[inc_mirror][2]/2.0))
                    baldwin_upper_alt_bound = numpy.degrees(math.atan(1.0/baldwin_upper_alt_bound))
                    if baldwin_upper_alt_bound < theta:
                        return False
                    else:
                        return True
            else:
                return False
        elif theta == 0:
            return True
        else:
            return False

    def test_tet_no_mirror(self):
        symclass = fu.symclass('tet')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_tet(phi, theta, 0, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 0))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results

    def test_tet_mirror(self):
        symclass = fu.symclass('tet')
        angles = [[entry, thet, 0] for entry in range(360) for thet in range(1, 180)]
        [angles.append([entry, 0, 0]) for entry in range(360)]
        [angles.append([entry, 180, 0]) for entry in range(360)]
        expected_results = []
        results = []
        for idx, entry in enumerate(angles):
            phi, theta, psi = entry
            expected_results.append(self.generate_result_tet(phi, theta, 1, symclass.brackets))
            results.append(symclass.is_in_subunit(phi, theta, 1))
            if theta == 0:
                assert expected_results[idx]
                assert results[idx]
        assert expected_results == results



########
########
########


class TestSymClassIsInSubunitWrong:

    output_template_angles = 'Got: {0} ; Expected: {1} ; Angle: {2}'

    def test_wrong_sym_crashes_problem(self):
        symclass = fu.symclass('c1')
        symclass.sym = 'foobar'
        assert symclass.is_in_subunit(0, 0) is None



########
########
########


class TestSymClassSymmetryRelatedC:

    def test_c1_sym_zero(self):
        symclass = fu.symclass('c1')
        result = symclass.symmetry_related([0, 0, 0])
        print(result)
        assert result == [[0, 0, 0]]

    def test_c5_sym_zero(self):
        symclass = fu.symclass('c5')
        result = symclass.symmetry_related([0, 0, 0])
        print(result)
        assert result == [[0, 0, 0]]

    def test_c1_sym_180(self):
        symclass = fu.symclass('c1')
        result = symclass.symmetry_related([0, 180, 0])
        print(result)
        assert result == [[0, 180, 0]]

    def test_c5_sym_180(self):
        symclass = fu.symclass('c5')
        result = symclass.symmetry_related([0, 180, 0])
        print(result)
        assert result == [[0, 180, 0]]

    def test_c1_sym_90(self):
        symclass = fu.symclass('c1')
        result = symclass.symmetry_related([0, 90, 0])
        print(result)
        assert result == [[0, 90, 0]]

    def test_c5_sym_90(self):
        symclass = fu.symclass('c5')
        result = symclass.symmetry_related([0, 90, 0])
        print(result)
        assert result == [[0, 90, 0], [72, 90, 0], [72*2, 90, 0], [72*3, 90, 0], [72*4, 90, 0]]



########
########
########


class TestSymClassSymmetryRelatedD:

    def test_d4_sym_zero(self):
        symclass = fu.symclass('d4')
        result = symclass.symmetry_related([0, 0, 0])
        print(result)
        assert result == [[0, 0, 0], [0, 180, 0]]

    def test_d5_sym_zero(self):
        symclass = fu.symclass('d5')
        result = symclass.symmetry_related([0, 0, 0])
        print(result)
        assert result == [[0, 0, 0], [0, 180, 180]]

    def test_d1_sym_90(self):
        symclass = fu.symclass('d1')
        result = symclass.symmetry_related([0, 90, 0])
        print(result)
        assert result == [[0, 90, 0]]

    def test_d4_sym_90(self):
        symclass = fu.symclass('d4')
        result = symclass.symmetry_related([0, 90, 0])
        print(result)
        expected_results = [[90*i, 90, 0] for i in range(4)]
        assert result == expected_results

    def test_d5_sym_90(self):
        symclass = fu.symclass('d5')
        result = symclass.symmetry_related([0, 90, 0])
        print(result)
        expected_results = [[72*i, 90, 0] for i in range(5)]
        assert result == expected_results

    def test_d1_sym_45(self):
        symclass = fu.symclass('d1')
        result = symclass.symmetry_related([0, 45, 0])
        print(result)
        assert sorted(result) == [[0, 45, 0], [0, 45+90, 180]]

    def test_d4_sym_45(self):
        symclass = fu.symclass('d4')
        result = symclass.symmetry_related([0, 45, 0])
        print(result)
        expected_results = [[90*i, 45, 0] for i in range(4)]
        expected_results.extend([[90*i, 45+90, 0] for i in range(4)])
        assert sorted(result) == sorted(expected_results)

    def test_d5_sym_45(self):
        symclass = fu.symclass('d5')
        result = symclass.symmetry_related([0, 45, 0])
        print(result)
        expected_results = [[72*i, 45, 0] for i in range(5)]
        expected_results.extend([[72*i, 45+90, 180] for i in range(5)])
        assert sorted(result) == sorted(expected_results)


########
########
########


class TestSymClassSymmetryRelatedTetIcosOct:

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

    def test_returns_correct_tet(self):
        symclass = fu.symclass('tet')
        angles = [0, 0, 0]
        return_value = symclass.symmetry_related(angles)

        expected_return = [angles]
        mat = self.rotmatrix(angles[0],angles[1],angles[2])
        for l in range(1, 12):
            p1,p2,p3 = self.recmat( self.mulmat( mat , symclass.symatrix[l]) )
            expected_return.append([p1,p2,p3])
        print(expected_return)
        print(return_value)
        assert expected_return == return_value

    def test_returns_correct_oct(self):
        symclass = fu.symclass('oct')
        angles = [0, 0, 0]
        return_value = symclass.symmetry_related(angles)

        expected_return = [angles]
        mat = self.rotmatrix(angles[0],angles[1],angles[2])
        for l in range(1, 24):
            p1,p2,p3 = self.recmat( self.mulmat( mat , symclass.symatrix[l]) )
            expected_return.append([p1,p2,p3])
        print(expected_return)
        print(return_value)
        assert expected_return == return_value

    def test_returns_correct_oct(self):
        symclass = fu.symclass('icos')
        angles = [0, 0, 0]
        return_value = symclass.symmetry_related(angles)

        expected_return = [angles]
        mat = self.rotmatrix(angles[0],angles[1],angles[2])
        for l in range(1, 60):
            p1,p2,p3 = self.recmat( self.mulmat( mat , symclass.symatrix[l]) )
            expected_return.append([p1,p2,p3])
        print(expected_return)
        print(return_value)
        assert expected_return == return_value


########
########
########


class TestSymClassSymmetryRelated:

    def my_method(self, angles, sym, symatrix):
        if( sym[0] == "c" or sym[0] == "d" ):
            temp = EMAN2_cppwrap.Util.symmetry_neighbors(angles, sym)
            nt = len(temp)//3
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
                p1,p2,p3 = self.recmat( self.mulmat( mat , symatrix[l]) )
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
        symclass = fu.symclass('c4')
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, symclass.sym, symclass.symatrix)
        return_values = symclass.symmetry_neighbors(angles)
        assert return_values == expected_return_values

    def test_d_sym(self):
        symclass = fu.symclass('d4')
        angles = [[idx1, idx2, 0] for idx1 in range(90) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, symclass.sym, symclass.symatrix)
        return_values = symclass.symmetry_neighbors(angles)
        assert return_values == expected_return_values

    def test_tet_sym(self):
        symclass = fu.symclass('tet')
        angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, symclass.sym, symclass.symatrix)
        return_values = symclass.symmetry_neighbors(angles)
        assert return_values == expected_return_values

    def test_oct_sym(self):
        symclass = fu.symclass('oct')
        angles = [[idx1, idx2, 0] for idx1 in range(40) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, symclass.sym, symclass.symatrix)
        return_values = symclass.symmetry_neighbors(angles)
        assert return_values == expected_return_values

    def test_icos_sym(self):
        symclass = fu.symclass('icos')
        angles = [[idx1, idx2, 0] for idx1 in range(30) for idx2 in range(90)]
        expected_return_values = self.my_method(angles, symclass.sym, symclass.symatrix)
        return_values = symclass.symmetry_neighbors(angles)
        assert return_values == expected_return_values


########
########
########


class TestReduceAnglesets:

    def test_c_sym_no_mirror(self):
        angles = [[idx1, idx2, 0] for idx1 in range(360) for idx2 in range(181)]
