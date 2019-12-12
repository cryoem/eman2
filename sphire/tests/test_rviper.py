from __future__ import print_function
from __future__ import division

from numpy import allclose
from sp_utilities import get_im
from sphire.bin import sp_rviper as oldfu
from sphire.utils.SPHIRE.bin import sp_rviper as fu
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir
import unittest
from os import path,system as os_system

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

MPI_PATH = "/home/lusnig/SPHIRE_1_1/envs/sphire1_3/bin/mpirun" #"/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/"
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

#testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Class2D","best.hdf"), out_dir_old, "--radius=29"]
#testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_rviper.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Class2D","best.hdf"), out_dir_new, "--radius=29"]
#        out_dir_old="oldrvipe"
        #out_dir_new = "newrvipe"
class Test_Error_cases(unittest.TestCase):

    def test_invalid_number_params(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_rviper.py")]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            fu.main()
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[11].split("ERROR")[1],' => Invalid number of parameters used. Please see usage information above.')
        self.assertEqual(print_new.getvalue().split('\n')[11].split("ERROR")[1],print_old.getvalue().split('\n')[11].split("ERROR")[1])

    def test_invalid_numbers_of_processor(self):
        out_dir_old = "oldrvipe"
        out_dir_new = "newrvipe"
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_rviper.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"), out_dir_old, "--radius=29"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_rviper.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf"), out_dir_new, "--radius=29"]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        old=print_old.getvalue().split('\n')[5].split("ERROR")[1]
        nw=print_new.getvalue().split('\n')[5].split("ERROR")[1]
        self.assertEqual(old," => Number of processes needs to be a multiple of the number of quasi-independent runs (shc) within each viper run. Total quasi-independent runs by default are 3, you can change it by specifying --n_shc_runs option (in sxviper this option is called --nruns). Also, to improve communication time it is recommended that the number of processes divided by the number of quasi-independent runs is a power of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has).")
        self.assertEqual(old,nw)


# it fails!!!
class Test_run(unittest.TestCase):
    #it is not the same run of the tutorial because I minimized the iteration to speed up the test. It takes anyway an hour
    def test_(self):
        out_dir_old = "oldrviper"
        out_dir_new = "newrviper"
        filename_avg = "average_volume_001.hdf"
        filename_var = "variance_volume_001.hdf"
        os_system(
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_rviper.py")
            +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf")
            +" "+out_dir_old
            +" --sym=c5"
            + " --n_rv_runs=1"
            +" --n_shc_runs=1")
        os_system(
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_BIN_PATH,"sp_rviper.py")
            +" "+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"04_ISAC","best.hdf")
            +" "+out_dir_new
            + " --sym=c5"
            +" --n_rv_runs=1"
            +" --n_shc_runs=1")

        return_new_avg = get_im( path.join(out_dir_new,filename_avg) )
        return_new_var = get_im( path.join(out_dir_new,filename_var) )
        return_old_avg = get_im( path.join(out_dir_old,filename_avg) )
        return_old_var = get_im( path.join(out_dir_old,filename_var) )
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=0.1))
        self.assertTrue(allclose(return_old_var.get_3dview(), return_new_var.get_3dview(), atol=0.1))

        #remove_dir(out_dir_new)
        #remove_dir(out_dir_old)