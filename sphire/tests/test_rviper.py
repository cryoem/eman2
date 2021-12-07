from __future__ import print_function
from __future__ import division

from numpy import allclose
from sphire.libpy.sp_utilities import get_im


from tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_BIN_PATH,remove_dir
import unittest
from os import path
from bin_py3 import sp_rviper as oldfu
from sphire.bin import sp_rviper as fu

import shutil

ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW = " Adnan removed this folder"

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

import subprocess
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))



MPI_PATH = shutil.which("mpi_run")
NUM_PROC = 8

"""
WHAT IS MISSING:
1) plot_errors_between_any_number_of_projections is plotting stuff. i'm not going to test it
2) All the helper functions. Since I have to run it using mpirun I have, at the moment, diffilculties to fill
    the input values


RESULT AND KNOWN ISSUES
The compatibility test of the run test failed even if I increase the tollerance value for a x10 factor

"""

'''
class Test_helperFunctions(unittest.TestCase):
    def test_calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(self):
        return_new = fu.calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            masterdir=,
            rviper_iter=,
            criterion_name=,
            symc=,
            runs_ite=)
        return_old = old.calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            masterdir=,
            rviper_iter=,
            criterion_name=,
            symc=,
            runs_ite=)



    def test_identify_outliers(self):
        return_new = fu.identify_outliers(
            myid=,
            main_node=,
            rviper_iter=,
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            masterdir=,
            bdb_stack_location=,
            outlier_percentile=,
            criterion_name=,
            outlier_index_threshold_method=,
            angle_threshold=,
            symc=,
            options=,
            runs_iter=, )
        return_old = oldfu.identify_outliers(
            myid=,
            main_node=,
            rviper_iter=,
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            masterdir=,
            bdb_stack_location=,
            outlier_percentile=,
            criterion_name=,
            outlier_index_threshold_method=,
            angle_threshold=,
            symc=,
            options=,
            runs_iter=, )


    def test_find_index_of_discontinuity_in_derivative(self):
        return_new = fu.find_index_of_discontinuity_in_derivative(
            error_curve_func=,
            list_of_projection_indices=,
            mainoutputdir=,
            outlier_percentile=,
            runs_iter=, )
        return_old = oldfu.find_index_of_discontinuity_in_derivative(
            error_curve_func=,
            list_of_projection_indices=,
            mainoutputdir=,
            outlier_percentile=,
            runs_iter=, )

    # it is called in the 'calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination' helper's function
    def test_measure_for_outlier_criterion(self):
        return_new = fu.measure_for_outlier_criterion(
            criterion_name=, masterdir=, rviper_iter=, list_of_viper_run_indices=, symc=)
        return_old = oldfu.measure_for_outlier_criterion(
            criterion_name=, masterdir=, rviper_iter=, list_of_viper_run_indices=, symc=)

    def test_found_outliers(self):
        return_new = fu.found_outliers(
            list_of_projection_indices=,
            outlier_percentile=,
            rviper_iter=,
            masterdir=,
            bdb_stack_location=,
            outlier_index_threshold_method=,
            angle_threshold=,
            symc=,
            options=,
            runs_iter=, )
        return_old = oldfu.found_outliers(
            list_of_projection_indices=,
            outlier_percentile=,
            rviper_iter=,
            masterdir=,
            bdb_stack_location=,
            outlier_index_threshold_method=,
            angle_threshold=,
            symc=,
            options=,
            runs_iter=, )

    def test_calculate_volumes_after_rotation_and_save_them(self):
        fu.calculate_volumes_after_rotation_and_save_them(
            ali3d_options=,
            rviper_iter=,
            masterdir=,
            bdb_stack_location=,
            mpi_rank=,
            mpi_size=,
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            mpi_comm=-1)
        oldfu.calculate_volumes_after_rotation_and_save_them(
            ali3d_options=,
            rviper_iter=,
            masterdir=,
            bdb_stack_location=,
            mpi_rank=,
            mpi_size=,
            no_of_viper_runs_analyzed_together=,
            no_of_viper_runs_analyzed_together_from_user_options=,
            mpi_comm=-1)

    def test_get_already_processed_viper_runs_False(self):
        return_new = fu.get_already_processed_viper_runs(run_get_already_processed_viper_runs=False)
        return_old = oldfu.get_already_processed_viper_runs(run_get_already_processed_viper_runs=False)
        self.assertIsNone(return_new)
        self.assertIsNone(return_old)

    @unittest.skip("crash because StopIteration error)
    def test_get_already_processed_viper_runs_True_(self):
        #crash because "StopIteration error" in next(os.walk(location_location))
        return_new = fu.get_already_processed_viper_runs(run_get_already_processed_viper_runs=True)
        return_old = oldfu.get_already_processed_viper_runs(run_get_already_processed_viper_runs=True)
'''

class Test_Error_cases(unittest.TestCase):
    @unittest.skip("Path not found,i ll fix it ")
    def test_invalid_number_params(self):
        a = subprocess.run(args=[path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py")], shell=True, capture_output=True)
        b = subprocess.run(args=[path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py")], shell=True,  capture_output=True)
        self.assertTrue(' => Invalid number of parameters used. Please see usage information above.' in b.stdout.decode('utf8'))
        self.assertTrue(' => Invalid number of parameters used. Please see usage information above.' in a.stdout.decode('utf8'))

    @unittest.skip("Adnan removed the ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW folder link, No MPI_path")
    def test_invalid_numbers_of_processor(self):
        testargs_new = (
                MPI_PATH
                + " -np 5"
                + " " + path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py")
                        +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf")
                        +" newrvipe  --radius=29")
        testargs_old =  (
                MPI_PATH
                + " -np 5"
                + " " + path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py")
                        +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf")
                        +" oldrvipe  --radius=29")

        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)
        self.assertTrue(" => Number of processes needs to be a multiple of the number of quasi-independent runs (shc) within each viper run. Total quasi-independent runs by default are 3, you can change it by specifying --n_shc_runs option (in sxviper this option is called --nruns). Also, to improve communication time it is recommended that the number of processes divided by the number of quasi-independent runs is a power of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has)." in b.stdout.decode('utf8'))
        self.assertTrue(" => Number of processes needs to be a multiple of the number of quasi-independent runs (shc) within each viper run. Total quasi-independent runs by default are 3, you can change it by specifying --n_shc_runs option (in sxviper this option is called --nruns). Also, to improve communication time it is recommended that the number of processes divided by the number of quasi-independent runs is a power of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has)." in a.stdout.decode('utf8'))



# the bin_py3 version is buggy, we cannot run the compatiubily test
@unittest.skip("Adnan removed the ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW folder link, No MPI_path")
class Test_run(unittest.TestCase):
    #it is not the same run of the tutorial because I minimized the iteration to speed up the test.
    def test_(self):
        out_dir_old = "oldrviper"
        out_dir_new = "newrviper"
        filename_avg = "average_volume.hdf"
        filename_var = "variance_volume.hdf"


        testargs_new= (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_BIN_PATH,"sp_rviper.py")
            +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","class_averages.hdf")
            +" "+out_dir_new
            + " --sym=c5"
            + " --n_rv_runs=1"
            + " --n_shc_runs=1")

        testargs_old = (
                        MPI_PATH
                        + " -np "
                        +str(NUM_PROC)
                        +" "+path.join(ABSOLUTE_OLDBIN_PATH, "sp_rviper.py")
                        +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","class_averages.hdf")
                        + " " + out_dir_old
                        + " --sym=c5"
                        + " --n_rv_runs=1"
                        + " --n_shc_runs=1")


        subprocess.run(args=[testargs_new], shell=True,  capture_output=True)


        return_new_avg = get_im( path.join(ABSOLUTE_PATH,out_dir_new,"main001",filename_avg) )
        return_new_var = get_im( path.join(ABSOLUTE_PATH,out_dir_new,"main001",filename_var))

        self.assertTrue(allclose(return_new_avg.get_3dview().flatten().tolist()[54381:54481], [0.001278397161513567, 0.00135961570776999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.002701209858059883, 0.0160306915640831, 0.0298649650067091, 0.02618340216577053, 0.013211097568273544, 0.006482064723968506, 0.004847334697842598, 0.014219394885003567, 0.029471751302480698, 0.024112213402986526, 0.005161711946129799, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], atol=0.05))
        self.assertTrue(allclose(return_new_var.get_3dview().flatten().tolist()[49057:49157], [6.550220366108306e-10, 7.117464519978967e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.6259397706572827e-09, 5.740563938161358e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], atol=0.005))
        remove_dir(out_dir_new)
        '''
        subprocess.run(args=[testargs_old], shell=True,  capture_output=True)
        return_old_avg = get_im( path.join(out_dir_old,filename_avg) )
        return_old_var = get_im( path.join(out_dir_old,filename_var) )
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=0.1))
        self.assertTrue(allclose(return_old_var.get_3dview(), return_new_var.get_3dview(), atol=0.1))

        #remove_dir(out_dir_new)
        #remove_dir(out_dir_old)
        '''
