from __future__ import print_function
from __future__ import division




from sphire.bin import sp_isac2 as oldfu
from sphire.utils.SPHIRE.bin import sp_isac2 as fu
from os import path,listdir
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,ABSOLUTE_BIN_PATH,remove_list_of_file
import unittest

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


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):

    def test_no_radius_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov", ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov",]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            return_new=fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            return_old=oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 1)
        self.assertEqual(print_new.getvalue().split('\n')[1],' ==radius== mandatory option is missing.')
        self.assertEqual(print_new.getvalue().split('\n')[1],print_old.getvalue().split('\n')[1])

    # todo: there is an error in li = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0] when I go throught it in the second 'main' call. it does not depend on the lib version. it seems to be an mpi syncronization situation ... ask markus
    def test_negative_radius_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov", "--radius=-145"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov2","--radius=-145"]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[20].split("ERROR")[1],' => Particle radius has to be provided!')
        self.assertEqual(print_new.getvalue().split('\n')[20].split("ERROR")[1],print_old.getvalue().split('\n')[20].split("ERROR")[1])

    def test_minimum_group_size_higher_group_size_errowr(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov", "--radius=145", "--minimum_grp_size=3", "--img_per_grp=2"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov","--radius=145", "--minimum_grp_size=3", "--img_per_grp=2"]
        with patch.object(sys, 'argv', testargs_new):
            old_stdout = sys.stdout
            print_new = StringIO()
            sys.stdout = print_new
            return_new=fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            print_old = StringIO()
            sys.stdout = print_old
            return_old=oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 1)
        self.assertEqual(print_new.getvalue().split('\n')[1],'ERROR! Minimum group size (3) is larger than the actual group size (2). Oh dear :(')
        self.assertEqual(print_new.getvalue().split('\n')[1],print_old.getvalue().split('\n')[1])

    #todo: there is an error in li = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0] when I go throught it in the second 'main' call. it does not depend on the lib version. it seems to be an mpi syncronization situation ... ask markus
    def test_CTF_and_VPP_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov", "--radius=145", "--CTF", "--VPP"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "/home/lusnig/Downloads/luca5nov2","--radius=145", "--CTF", "--VPP"]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main(testargs_new[1:])
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main(testargs_old[1:])
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[20].split("ERROR")[1],'=> Options CTF and VPP cannot be used together')
        self.assertEqual(print_new.getvalue().split('\n')[20].split("ERROR")[1],print_old.getvalue().split('\n')[20].split("ERROR")[1])

