from __future__ import print_function
from __future__ import division


from numpy import allclose
from sphire.libpy.sp_utilities import get_im

from os import path
from tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_BIN_PATH, remove_dir
import unittest

#todo: need these data
ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW = "removed by Adnan"

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

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))
TOLERANCE=190


from sphire.libpy import sp_global_def
import shutil

sp_global_def.BATCH = True
sp_global_def.MPI = True

MPI_PATH = shutil.which("mpi_run")
NUM_PROC = 6  # has to be a multiple of 3

import subprocess
"""
WHAT IS MISSING:
I just test the error case because the script collects the input values from the gui and then call "multi_shc" from sp_multi_shc.
"""

@unittest.skip("nov_30 IT FAILS because there are not some folders in Adnan files, and the MPI_PATH is unknown")
class Test_Error_cases(unittest.TestCase):
    def test_too_few_input_values(self):
        b=subprocess.run(args=[path.join(ABSOLUTE_BIN_PATH, "sp_viper.py")], shell=True,  capture_output=True)
        a=subprocess.run(args=[path.join(ABSOLUTE_OLDBIN_PATH, "sp_viper.py")], shell=True,  capture_output=True)

        self.assertTrue(' => Invalid number of parameters used. Please see usage information above.' in b.stdout.decode('utf8'))
        self.assertTrue(' => Invalid number of parameters used. Please see usage information above.' in a.stdout.decode('utf8'))


    def test_invalid_nruns(self):
        testargs_new= (
            path.join(ABSOLUTE_BIN_PATH, "sp_viper.py")
            +' bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_Particles#stack")
            +" out_new_dir")

        testargs_old= (
            path.join(ABSOLUTE_BIN_PATH, "sp_viper.py")
            +' bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_Particles#stack")
            +" out_old_dir")

        b=subprocess.run(args=[testargs_old], shell=True,  capture_output=True)
        a=subprocess.run(args=[testargs_new], shell=True,  capture_output=True)

        self.assertTrue(' => Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.' in b.stdout.decode('utf8'))
        self.assertTrue(' => Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.' in a.stdout.decode('utf8'))


#todo: check if it works. the tolerance value is huge
@unittest.skip("nov_30 IT FAILS because there are not some folders in Adnan files, and the MPI_PATH is unknown")
class Test_run(unittest.TestCase):
    def test_(self):
        out_dir_old = "oldviper"
        out_dir_new = "newviper"

        filename_vol = "volf.hdf"
        filename_refvol = "refvol2.hdf"

        testcommand_new = (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_viper.py")
            +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","class_averages.hdf")
            +" "+out_dir_old
            +" --sym=c5")


        testcommand_old= (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_BIN_PATH,"sp_viper.py")
            +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","class_averages.hdf")
            +" "+out_dir_new
            + " --sym=c5")



        subprocess.run(args =[testcommand_new], shell=True, stderr=subprocess.STDOUT)
        subprocess.run(args=[testcommand_old], shell=True, stderr=subprocess.STDOUT)



        return_new_avg = get_im( path.join(ABSOLUTE_PATH,out_dir_new,filename_vol) )
        return_new_var = get_im( path.join(ABSOLUTE_PATH,out_dir_new,filename_refvol))

        return_old_avg = get_im( path.join(ABSOLUTE_PATH,out_dir_old,filename_vol) )
        return_old_var = get_im( path.join(ABSOLUTE_PATH,out_dir_old,filename_refvol))

        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE))
        self.assertTrue(allclose(return_old_var.get_3dview(), return_new_var.get_3dview(), atol=TOLERANCE))

        remove_dir(out_dir_new)
        remove_dir(out_dir_old)




"""
platform linux2 -- Python 2.7.14, pytest-4.6.6, py-1.8.0, pluggy-0.13.0
rootdir: /home/lusnig/src_sphire_1_3/eman2/sphire/tests
plugins: cov-2.8.1
collected 1 item

test_viper.py [Start] : 2019-11-25 11:44:45
2019-11-25 11:44:45 write_command => /home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py /home/lusnig/SphireDemoResults2/04_ISAC/best.hdf oldviper --sym=c5

2019-11-25 11:45:12 logLine =>  Start VIPER1
2019-11-25 11:45:12 logLine =>  ITERATION #  1
2019-11-25 11:46:34 logLine =>  Time to prepare rings: 81.971819
2019-11-25 11:46:34 logLine =>  Time of alignment = 0.320021

2019-11-25 11:46:34 logLine =>  = Pixel error        Number of images in all runs
2019-11-25 11:46:34 logLine =>        4.651                        4
2019-11-25 11:46:34 logLine =>        7.318                        5
2019-11-25 11:46:34 logLine =>        9.985                        9
2019-11-25 11:46:34 logLine =>       12.653                        4
2019-11-25 11:46:34 logLine =>       15.320                       14
2019-11-25 11:46:34 logLine =>       17.988                       17
2019-11-25 11:46:34 logLine =>       20.655                       16
2019-11-25 11:46:34 logLine =>       23.323                       16
2019-11-25 11:46:34 logLine =>       25.990                       22
2019-11-25 11:46:34 logLine =>       28.658                       23
2019-11-25 11:46:34 logLine =>       31.325                       23
2019-11-25 11:46:34 logLine =>       33.992                       20
2019-11-25 11:46:34 logLine =>       36.660                       31
2019-11-25 11:46:34 logLine =>       39.327                       29
2019-11-25 11:46:34 logLine =>       41.995                       41
2019-11-25 11:46:34 logLine =>       44.662                       49
2019-11-25 11:46:34 logLine =>       47.330                       44
2019-11-25 11:46:34 logLine =>       49.997                       45
2019-11-25 11:46:34 logLine =>       52.665                       60
2019-11-25 11:46:34 logLine =>       55.332                      146
2019-11-25 11:46:34 logLine =>  =================================================
2019-11-25 11:46:34 logLine =>  Percent of positions with pixel error below 1.0 =  0 %    Mutations:  False
2019-11-25 11:46:34 logLine =>   Average center x =      0.000        Center y =      0.000        Center z =      0.000

2019-11-25 11:46:34 logLine =>  For symmetry group cn (n>1), we only center the volume in z-direction

Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 211, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 206, in main
    out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 211, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 206, in main
    out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 211, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 206, in main
    out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
-------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code.. Per user-direction, the job has been aborted.
-------------------------------------------------------
Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 211, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/bin/sp_viper.py", line 206, in main
    out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
--------------------------------------------------------------------------
mpirun detected that one or more processes exited with non-zero status, thus causing
the job to be terminated. The first process to do so was:

  Process name: [[16163,1],2]
  Exit code:    1




        (sphire1_3) lusnig@ ~/src_sphire_1_3/eman2/sphire/tests (add_tests_sphire1_3) $pytest test_viper.py::Test_run::test_ --capture=no
=========================================================================================================================== test session starts ============================================================================================================================
platform linux2 -- Python 2.7.14, pytest-4.6.6, py-1.8.0, pluggy-0.13.0
rootdir: /home/lusnig/src_sphire_1_3/eman2/sphire/tests
plugins: cov-2.8.1
collected 1 item

test_viper.py [Start] : 2019-11-25 11:36:29
2019-11-25 11:36:29 write_command => /home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py /home/lusnig/SphireDemoResults2/04_ISAC/best.hdf newviper --sym=c5

2019-11-25 11:36:42 logLine =>  Start VIPER1
2019-11-25 11:36:42 logLine =>  ITERATION #  1
2019-11-25 11:38:04 logLine =>  Time to prepare rings: 81.609615
2019-11-25 11:38:04 logLine =>  Time of alignment = 0.284481

2019-11-25 11:38:04 logLine =>  = Pixel error        Number of images in all runs
2019-11-25 11:38:04 logLine =>        4.408                        2
2019-11-25 11:38:04 logLine =>        7.088                        3
2019-11-25 11:38:04 logLine =>        9.767                        3
2019-11-25 11:38:04 logLine =>       12.447                        7
2019-11-25 11:38:04 logLine =>       15.126                        8
2019-11-25 11:38:04 logLine =>       17.806                       10
2019-11-25 11:38:04 logLine =>       20.486                       19
2019-11-25 11:38:04 logLine =>       23.165                       18
2019-11-25 11:38:04 logLine =>       25.845                       18
2019-11-25 11:38:04 logLine =>       28.524                       24
2019-11-25 11:38:04 logLine =>       31.204                       28
2019-11-25 11:38:04 logLine =>       33.884                       29
2019-11-25 11:38:04 logLine =>       36.563                       26
2019-11-25 11:38:04 logLine =>       39.243                       32
2019-11-25 11:38:04 logLine =>       41.922                       37
2019-11-25 11:38:04 logLine =>       44.602                       50
2019-11-25 11:38:04 logLine =>       47.282                       44
2019-11-25 11:38:04 logLine =>       49.961                       51
2019-11-25 11:38:04 logLine =>       52.641                       66
2019-11-25 11:38:04 logLine =>       55.320                      143
2019-11-25 11:38:04 logLine =>  =================================================
2019-11-25 11:38:04 logLine =>  Percent of positions with pixel error below 1.0 =  0 %    Mutations:  False
2019-11-25 11:38:04 logLine =>   Average center x =      0.000        Center y =      0.000        Center z =      0.000

2019-11-25 11:38:04 logLine =>  For symmetry group cn (n>1), we only center the volume in z-direction

Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 336, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 330, in main
    all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
-------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code.. Per user-direction, the job has been aborted.
-------------------------------------------------------
Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 336, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 330, in main
    all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
Traceback (most recent call last):
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 336, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 330, in main
    all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2366, in do_volume
    if options.filament_width != -1:
AttributeError: Values instance has no attribute 'filament_width'
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 336, in <module>
    main(sys.argv[1:])
  File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/bin/sp_viper.py", line 330, in main
    all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 1449, in multi_shc
    out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 799, in ali3d_multishc
    vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_multi_shc.py", line 2356, in do_volume
    else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, snr=snr, npad=npad, mpi_comm=mpi_comm)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_reconstruction.py", line 248, in recons3d_4nn_MPI
    insert_slices(r, prj)
  File "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/lib/python2.7/site-packages/sp_reconstruction.py", line 84, in insert_slices
    reconstructor.insert_slice( proj, xforms[i], weights[i] )
Shutdown complete, exiting
--------------------------------------------------------------------------
mpirun detected that one or more processes exited with non-zero status, thus causing
the job to be terminated. The first process to do so was:

  Process name: [[16008,1],2]
  Exit code:    1
--------------------------------------------------------------------------
.

"""


