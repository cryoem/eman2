from __future__ import print_function
from __future__ import division



from numpy import allclose,array_equal



from os import path
from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_list_of_file
import unittest
from sp_utilities import get_im


import subprocess
MPI_PATH = "/home/adnan/applications/sphire/miniconda3/envs/py3_v5/bin/mpirun"
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


# from sphire.bin_py3 import sp_filterlocal as oldfu
# from sphire.bin import sp_filterlocal as fu

import sp_global_def
# import mpi
# mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
# mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)




class Test_run(unittest.TestCase):
    def test_calls(self):
        old_final="old_final.hdf"
        new_final = "new_final.hdf"

        testargs_new = [str(MPI_PATH) +
                         " -np"
                         " 4"
                         " /home/adnan/PycharmProjects/python3conversion/sphire/bin/sp_filterlocal.py"
                         " /home/adnan/DemoResults/12_POSTREFINER/vol_combined.hdf"
                         " /home/adnan/DemoResults/17_LOCAL_RES/localres.hdf"
                         " /home/adnan/DemoResults/12_POSTREFINER/vol_adaptive_mask.hdf"
                         " new_final.hdf"
                         " --radius=145"
                         " --MPI"]

        testargs_old = [str(MPI_PATH) +
                         " -np"
                         " 4"
                         " /home/adnan/PycharmProjects/python3conversion/sphire/bin_py3/sp_filterlocal.py"
                         " /home/adnan/DemoResults/12_POSTREFINER/vol_combined.hdf"
                         " /home/adnan/DemoResults/17_LOCAL_RES/localres.hdf"
                         " /home/adnan/DemoResults/12_POSTREFINER/vol_adaptive_mask.hdf"
                         " old_final.hdf"
                         " --radius=145"
                         " --MPI"]


        print(testargs_new)

        print(testargs_old)

        a = subprocess.run(args =testargs_new, shell=True, stderr=subprocess.STDOUT)

        b = subprocess.run(args=testargs_old, shell=True, stderr=subprocess.STDOUT)

        return_old = get_im(old_final)
        return_new = get_im(new_final)
        self.assertTrue(array_equal(return_old.get_3dview(),return_new.get_3dview()))
        self.assertTrue(allclose(return_new.get_3dview().flatten().tolist()[4034427:4034527],[3.201385334250517e-05, -0.0033583296462893486, -0.0033261957578361034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],atol=0.01))
        remove_list_of_file([old_final,new_final])

