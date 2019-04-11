from __future__ import print_function
from __future__ import division

import numpy
from math import isnan as math_isnan
from copy import deepcopy
from EMAN2_cppwrap import EMData,Util
from EMAN2_cppwrap import EMData, EMAN2Ctf
import unittest
import os
import pickle

import test_module as tm
import EMAN2_cppwrap as e2cpp

from sphire.libpy_py3 import sphire_projection as fu
from .sparx_lib import sparx_projection as oldfu

from .sparx_lib import sparx_utilities

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))

class Test_lib_projection_compare(unittest.TestCase):


    def test_project_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]


        return_new = fu.project(volft,params)
        return_old = oldfu.project(volft,params)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_prgs_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]

        kb = tm.create_kb(1)

        return_new = fu.prgs(volft,kb,params)
        return_old = oldfu.prgs(volft,kb,params)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_prgl_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)


        (volft, params, interpolation_method,return_real) = argum[0]

        return_new = fu.prgl(volft, params, interpolation_method)
        return_old = oldfu.prgl(volft, params, interpolation_method)


        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))

    """Function works but is too slow , takes 506 secs to test"""
    # def test_prgq_should_return_same_result(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     kb = tm.create_kb(1)
    #     (volft, params, interpolation_method,return_real) = argum[0]
    #
    #     return_new = fu.prgq(volft, kb, nx=4, delta=0.5, ref_a="S", sym="c1")
    #
    #     print("Hello")
    #     return_old = oldfu.prgq(volft, kb, nx=4, delta=0.5, ref_a="S", sym="c1")
    #     print("Hello")
    #
    #     self.assertTrue(return_new, return_old)

    def test_prg_should_return_same_result(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (volft, params, interpolation_method, return_real) = argum[0]
        volft = sparx_utilities.model_blank(100, 100, 100)

        return_new = fu.prg(volft, params)
        return_old = oldfu.prg(volft, params)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_prep_vol_should_return_same_result(self):
        volft = sparx_utilities.model_blank(100, 100, 100)

        return_new = fu.prep_vol(volft)
        return_old = oldfu.prep_vol(volft)


        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertTrue(return_new[1], return_old[1])

    from test_module import get_data
    image = get_data(1)[0]
    def test_gen_rings_ctf_should_return_same_result(self):
        argum = tm.get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ornq"))
        (image, crefim, xrng, yrng, step, mode, numr, cnx, cny) = argum[0]

        nx =  4
        ctf = e2cpp.EMAN2Ctf()
        ctf.from_dict({"defocus": 1, "cs": 2, "voltage": 300, "apix": 1.5, "bfactor": 0,"ampcont": 0.1})

        prjref = []
        prjref.append(self.image)
        prjref.append(self.image)

        prjref[0].set_attr('phi', 20)
        prjref[0].set_attr('theta', 40)
        prjref[0].set_attr('psi', 40)
        prjref[1].set_attr('phi', 30)
        prjref[1].set_attr('theta', 30)
        prjref[1].set_attr('psi', 30)


        return_new = fu.gen_rings_ctf(prjref, nx , ctf, numr)
        return_old = oldfu.gen_rings_ctf(prjref, nx , ctf, numr)
        self.assertEqual(return_new[0].get_attr_dict(), return_old[0].get_attr_dict())
        self.assertEqual(return_new[1].get_attr_dict(), return_old[1].get_attr_dict())


if __name__ == '__main__':
    unittest.main()
