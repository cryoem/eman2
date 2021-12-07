from __future__ import print_function
from __future__ import division



from numpy import allclose,array_equal
from os import path
from tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_BIN_PATH,remove_list_of_file
import unittest
from sphire.libpy.sp_utilities import get_im
import shutil
import subprocess

ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW = "Adnan removed it"

MPI_PATH = shutil.which("mpi_run")
NUM_PROC = 8
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


@unittest.skip("nov_21 IT FAILS because there is not '12_POSTREFINER'' folder in Adnan files ")
class Test_run(unittest.TestCase):
    def test_calls(self):
        old_final="old_final.hdf"
        new_final = "new_final.hdf"

        testargs_new = [str(MPI_PATH) +
                         " -np " + str(NUM_PROC)+
                         " "+path.join(ABSOLUTE_BIN_PATH,"sp_filterlocal.py")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_combined.hdf")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres.hdf")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf")+
                         " " + new_final+
                         " --radius=145 --MPI"]

        testargs_old = [str(MPI_PATH) +
                         " -np " + str(NUM_PROC)+
                         " "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_filterlocal.py")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_combined.hdf")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "17_LOCAL_RES","localres.hdf")+
                         " " + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "12_POSTREFINER","vol_adaptive_mask.hdf")+
                         " " + old_final+
                         " --radius=145 --MPI"]

        a = subprocess.run(args =testargs_new, shell=True, stderr=subprocess.STDOUT)

        b = subprocess.run(args=testargs_old, shell=True, stderr=subprocess.STDOUT)

        return_old = get_im(old_final)
        return_new = get_im(new_final)
        self.assertTrue(array_equal(return_old.get_3dview(),return_new.get_3dview()))
        self.assertTrue(allclose(return_new.get_3dview().flatten().tolist()[4034427:4034527],[3.201385334250517e-05, -0.0033583296462893486, -0.0033261957578361034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],atol=0.01))
        remove_list_of_file([old_final,new_final])

