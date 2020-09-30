from __future__ import print_function
from __future__ import division



from numpy import allclose,array_equal

from sphire.bin_py3 import sp_locres as oldfu
from ..sphire.bin import sp_locres as fu

from os import path
from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir,IMAGE_3D
import unittest
from sp_utilities import get_im

try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO
import sys




"""
WHAT IS MISSING:
1) 'output_volume' because it just creates the hdf files

RESULT AND KNOWN ISSUES


"""
class Test_helperFunction(unittest.TestCase):
    def test_makeAngRes(self):
        return_new= fu.makeAngRes(freqvol=IMAGE_3D, nx=IMAGE_3D.get_xsize(), ny=IMAGE_3D.get_ysize(), nz=IMAGE_3D.get_zsize(), pxSize=1.14, freq_to_real=True)
        return_old = oldfu.makeAngRes(freqvol=IMAGE_3D, nx=IMAGE_3D.get_xsize(), ny=IMAGE_3D.get_ysize(),nz=IMAGE_3D.get_zsize(), pxSize=1.14, freq_to_real=True)
        self.assertTrue(array_equal(return_old.get_3dview(), return_new.get_3dview()))
        self.assertTrue(array_equal(return_old.get_3dview().flatten().tolist()[0:100],[119.94329071044922, 44.04085922241211, 74.16339874267578, 38.44627380371094, 44.49072265625, 47.507747650146484, 48.8928108215332, 31.434528350830078, 26.989473342895508, 21.403867721557617, 16.29381561279297, 21.045682907104492, 22.35525894165039, 32.074127197265625, 31.6229190826416, 16.2719783782959, 20.086626052856445, 16.940921783447266, 12.672057151794434, 14.242024421691895, 15.819912910461426, 15.925355911254883, 13.410531044006348, 14.119293212890625, 12.700398445129395, 11.933272361755371, 11.71253490447998, 9.380097389221191, 11.660004615783691, 18.62543296813965, 18.851110458374023, 13.68883228302002, 14.266586303710938, 15.70208740234375, 10.971940994262695, 8.981616020202637, 12.670086860656738, 19.85761070251465, 20.27568244934082, 20.640796661376953, 81.96916961669922, 159.50694274902344, 75.41804504394531, 45228.44921875, 138.48532104492188, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

class Test_run(unittest.TestCase):
    out_dir_old="localResOld"
    out_dir_new = "localResNew"
    fname_ang = "localres_ang.hdf"
    fname = "localres.hdf"


    @classmethod
    def tearDownClass(cls):
        remove_dir(cls.out_dir_new)
        remove_dir(cls.out_dir_old)

    # it is the run of the tutorial (pg.65)
    # At pg.93 there is another run. I did not test it because it performs the same operation done in this test
    def test_(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "localRes.py"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_0_unfil_028.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_1_unfil_028.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf"),
                        self.out_dir_new,
                        "--out_ang_res",
                        "--apix=1.14"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_process.py"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_0_unfil_028.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","vol_1_unfil_028.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf"),
                        self.out_dir_old,
                        "--out_ang_res",
                        "--apix=1.14"]


        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()

        old_localRes_ang = get_im(path.join(self.out_dir_old, self.fname_ang))
        old_localRes = get_im(path.join(self.out_dir_old, self.fname))
        new_localRes_ang = get_im(path.join(self.out_dir_new, self.fname_ang))
        new_localRes = get_im(path.join(self.out_dir_new, self.fname))
        self.assertTrue(allclose(old_localRes_ang.get_3dview(),new_localRes_ang.get_3dview(),atol=0.05))
        self.assertTrue(allclose(old_localRes.get_3dview(), new_localRes.get_3dview(), atol=0.05))
        self.assertTrue(allclose(new_localRes_ang.get_3dview().flatten().tolist()[4034427:4034527],[6.859487056732178, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],atol=0.05))
        self.assertTrue(allclose(new_localRes.get_3dview().flatten().tolist()[4034427:4034527],[0.1661931872367859, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],atol=0.05))