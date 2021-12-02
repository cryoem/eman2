from __future__ import print_function
from __future__ import division

from os import path
from numpy import allclose,array_equal,full,nan
from copy import deepcopy

from sphire.libpy import sp_user_functions as fu
from libpy_py3 import sp_user_functions as oldfu
from tests.test_module import get_arg_from_pickle_file,ABSOLUTE_PATH_TO_RESOURCES

import unittest


"""
'build_user_function" function is used by the class 'factory_class' which is basically a factory method pattern used to call the other functions implemented in the file.
Since ALL the other functions are not used in the code I have not tested anything
All the functions implemented but not used:
-) ref_ali2d( ref_data )
-) ref_ali2d_c( ref_data )
-) julien(ref_data)
-) ref_ali2d_m( ref_data )
-) ref_ali3dm( refdata )
-) ref_sort3d(refdata)
-) ref_ali3dm_ali_50S( refdata )
-) ref_random( ref_data )
-) ref_ali3d( ref_data )
-) helical( ref_data )
-) helical2( ref_data )
-) reference3( ref_data )
-) reference4( ref_data )
-) ref_aliB_cone( ref_data )
-) ref_7grp( ref_data )
-) spruce_up (ref_data)
-) spruce_up_variance(ref_data)
-) minfilt( fscc )
-) ref_ali3dm_new( refdata )
-) spruce_up_var_m( refdata )
-) steady( refdata )
-) constant( refdata )
-) temp_dovolume( refdata )
-) dovolume( refdata )
-) do_volume_mask( refdata )
-) do_volume_mrk02( refdata )
-) do_volume_mrk03( refdata )
-) do_volume_mrk04( refdata )
-) do_volume_mrk05( refdata )

"""

refdata = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH_TO_RESOURCES, "user_functions.ref_ali2d_new"))[0][0]

class Test_lib_user_functions_compare(unittest.TestCase):

    def test_ref_ali2d_should_return_equal_objects(self):

        return_new = fu.ref_ali2d(deepcopy(refdata))
        return_old = oldfu.ref_ali2d(deepcopy(refdata))
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

        self.assertTrue(allclose(return_new[0].get_3dview(), full((70,70),nan),equal_nan=True))
        self.assertTrue(array_equal(return_new[1], [0.0, 0.0]))


    def test_ref_ali2d_c_true_should_return_equal_objects(self):
        return_new = fu.ref_ali2d_c(refdata)
        return_old = oldfu.ref_ali2d_c(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_julien_true_should_return_equal_objects(self):
        return_new = fu.julien(refdata)
        return_old = oldfu.julien(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_ref_ali2d_m_true_should_return_equal_objects(self):
        return_new = fu.ref_ali2d_m(refdata)
        return_old = oldfu.ref_ali2d_m(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(), atol=1.e-5, equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    # def test_ref_ali3dm_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.ref_ali2d")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (refdata) = argum[0][0]
    #
    #     print('a', refdata[0])
    #     print('b', refdata[1])
    #     print('b', refdata[2])
    #     print('d', refdata[3])
    #
    #     refdata.extend("varf")
    #     refdata.extend("mask")
    #
    #     return_new = fu.ref_ali3dm(refdata)
    #     return_old = oldfu.ref_ali3dm(refdata)

    # self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
    # self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_ref_ali3d_true_should_return_equal_objects(self):
        return_new = fu.ref_ali3d(refdata)
        return_old = oldfu.ref_ali3d(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_helical_true_should_return_equal_objects(self):
        return_new = fu.helical(refdata)
        return_old = oldfu.helical(refdata)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(),atol=1.e-5,equal_nan=True))

    def test_helical2_true_should_return_equal_objects(self):
        return_new = fu.helical2(refdata)
        return_old = oldfu.helical2(refdata)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), atol=1.e-5, equal_nan=True))

    def test_reference3_true_should_return_equal_objects(self):
        return_new = fu.reference3(refdata)
        return_old = oldfu.reference3(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_reference4_true_should_return_equal_objects(self):
        return_new = fu.reference4(refdata)
        return_old = oldfu.reference4(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    # def test_ref_aliB_cone_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.ref_ali2d")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (refdata) = argum[0][0]
    #     refdata[1] = 1
    #
    #     return_new = fu.ref_aliB_cone(refdata)
    #     return_old = oldfu.ref_aliB_cone(refdata)
    #
    #     self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
    #     self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_ref_7grp_true_ZeroDivisionError(self):
        refdata_=deepcopy(refdata)
        refdata_[1] = 1
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ref_7grp(refdata_)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ref_7grp(refdata_)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    # def test_spruce_up_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.ref_ali2d")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (refdata) = argum[0][0]
    #
    #     return_new = fu.spruce_up(refdata)
    #     return_old = oldfu.spruce_up(refdata)
    #
    #     self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
    #     self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    def test_steady_true_should_return_equal_objects(self):
        return_new = fu.steady(refdata)
        return_old = oldfu.steady(refdata)
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_constant_true_should_return_equal_objects(self):
        return_new = fu.constant(deepcopy(refdata))
        return_old = oldfu.constant(deepcopy(refdata))
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_temp_dovolume_true_should_return_equal_objects(self):
        return_new = fu.temp_dovolume(deepcopy(refdata))
        return_old = oldfu.temp_dovolume(deepcopy(refdata))
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    def test_dovolume_true_should_return_equal_objects(self):
        return_new = fu.dovolume(deepcopy(refdata))
        return_old = oldfu.dovolume(deepcopy(refdata))
        self.assertTrue(allclose(return_new[0].get_3dview(), return_old[0].get_3dview(),atol=1.e-5,equal_nan=True))
        self.assertTrue(array_equal(return_new[1], return_old[1]))

    @unittest.skip('dont have mask.hdf')
    def test_do_volume_mask_true_should_return_equal_objects(self):
        filepath = os.path.join(
            ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask"
        )
        with open(filepath, "rb") as rb:
            argum = pickle.load(rb)
        refdata = []
        refdata.append(argum[0][0][0])
        refdata.append(argum[0][0][1])
        refdata.append(argum[0][0][2])

        refdata[1]["constants"][
            "mask3D"
        ] = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.do_volume_mask(deepcopy(refdata))
        return_old = oldfu.do_volume_mask(deepcopy(refdata))

        self.assertTrue(
            numpy.array_equal(return_new.get_3dview(), return_old.get_3dview())
        )

    @unittest.skip('dont have mask.hdf')
    def test_do_volume_mrk03_true_should_return_equal_objects(self):
        with open(filepath_mask, "rb") as rb:
            argum = pickle.load(rb, encoding='latin1')

        argum[0][0][1]["lowpass"] = 0.08
        argum[0][0][1]["falloff"] = 0.1
        refdata = []
        refdata.append(argum[0][0][0])
        refdata.append(argum[0][0][1])
        refdata.append(0)
        refdata.append(12)

        refdata[1]["constants"][
            "mask3D"
        ] = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.do_volume_mrk03(deepcopy(refdata))
        return_old = oldfu.do_volume_mrk03(deepcopy(refdata))

        self.assertTrue(
            numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview())
        )

    @unittest.skip('dont have mask.hdf')
    def test_do_volume_mrk04_true_should_return_equal_objects(self):
        filepath = os.path.join(
            ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask"
        )
        with open(filepath, "rb") as rb:
            argum = pickle.load(rb)

        argum[0][0][1]["lowpass"] = 0.08
        argum[0][0][1]["falloff"] = 0.1
        refdata = []
        refdata.append(argum[0][0][0])
        refdata.append(argum[0][0][1])
        refdata.append(0)
        refdata.append(12)

        refdata[1]["constants"][
            "mask3D"
        ] = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.do_volume_mrk04(deepcopy(refdata))
        return_old = oldfu.do_volume_mrk04(deepcopy(refdata))

        self.assertTrue(
            numpy.array_equal(return_new.get_3dview(), return_old.get_3dview())
        )

    @unittest.skip('dont have mask.hdf')
    def test_do_volume_mrk05_true_should_return_equal_objects(self):
        filepath = os.path.join(
            ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask"
        )
        with open(filepath, "rb") as rb:
            argum = pickle.load(rb)

        refdata = []
        refdata.append(argum[0][0][0])
        refdata.append(argum[0][0][1])

        refdata[1]["constants"][
            "mask3D"
        ] = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.do_volume_mrk05(deepcopy(refdata))
        return_old = oldfu.do_volume_mrk05(deepcopy(refdata))

        self.assertTrue(
            numpy.array_equal(return_new.get_3dview(), return_old.get_3dview())
        )



if __name__ == "__main__":
    unittest.main()
