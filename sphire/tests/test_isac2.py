from __future__ import print_function
from __future__ import division


import subprocess
MPI_PATH = "/home/adnan/applications/sphire/miniconda3/envs/py3_v5/bin/mpirun"
from numpy import allclose

# from sphire.bin_py3 import sp_isac2 as oldfu
# from sphire.bin import sp_isac2 as fu

from sphire.libpy.sp_utilities import get_im
from os import path
from sphire.tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir
import unittest
import numpy
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
1) 2 tests in the 'Test_Error_cases' class are not working because the call to 'mpi.mpi_bcast'. See the #todo 
2) since I have to run it via console using mpirun I'm not able to debug it ... hence I did not implement tests for
    the helper functions


RESULT AND KNOWN ISSUES
1) we have to run it using mpirun


In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:
1) the clean version cannot run!!!!!!!! 
2) Test_Error_cases::test_negative_radius_error 
"""


MPI_PATH = "/home/adnan/applications/sphire/miniconda3/envs/py3_v5/bin/mpirun" #"/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/"
NUM_PROC = 8
class Test_run(unittest.TestCase):

    old_output_folder="IsacOld"
    new_output_folder = "IsacNew"
    filename = "class_averages.hdf"


    @classmethod
    def tearDownClass(cls):
        remove_dir(cls.new_output_folder)
        remove_dir(cls.old_output_folder)

    def test_run(self):

        testcommand_new = (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_BIN_PATH,"sp_isac2.py")
            +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_Particles#stack'")
            +" '"+self.new_output_folder+"' "
            +" --radius=145"
            +" --CTF")

        testcommand_old= (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_isac2.py")
            +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_Particles#stack'")
            +" '"+self.old_output_folder+"' "
            +" --radius=145"
            +" --CTF")

        a = subprocess.run(args =[testcommand_new], shell=True, stderr=subprocess.STDOUT)
        b = subprocess.run(args=[testcommand_old], shell=True, stderr=subprocess.STDOUT)

        # old_value = get_im(path.join(self.old_output_folder, self.filename))
        # new_value = get_im(path.join(self.new_output_folder, self.filename))
        # self.assertTrue(numpy.array_equal(old_value.get_2dview(),new_value.get_2dview()))
        # self.assertTrue(allclose( old_value.get_2dview().flatten().tolist()[0:100], [0.02826664224267006, 0.027966178953647614, 0.02758362889289856, 0.027607476338744164, 0.028002280741930008, 0.028992081061005592, 0.02917819656431675, 0.025063568726181984, 0.02249879390001297, 0.02685454860329628, 0.03184734284877777, 0.03102857992053032, 0.03166065365076065, 0.0392359234392643, 0.04360540211200714, 0.036976393312215805, 0.01727971062064171, 0.023975739255547523, 0.039678484201431274, 0.03920724242925644, 0.048566970974206924, 0.05904018506407738, 0.06332758814096451, 0.059697940945625305, 0.02878061681985855, 0.020571228116750717, 0.05693618208169937, 0.05775655433535576, 0.04874090477824211, 0.03842538222670555, 0.017679473385214806, 0.007761257700622082, 0.012288063764572144, 0.023420047014951706, 0.018611961975693703, 0.005126635078340769, 0.01053755171597004, 0.04059397801756859, 0.09011068940162659, 0.12636150419712067, 0.07826534658670425, 0.010442654602229595, -0.013606780208647251, -0.014976754784584045, 0.02649846114218235, 0.06628189980983734, 0.06066082417964935, 0.06835232675075531, 0.030094420537352562, -0.018862448632717133, 0.00870553869754076, 0.030758436769247055, 0.029359430074691772, 0.04237156733870506, 0.05532151088118553, 0.07659796625375748, 0.05731084942817688, 0.021425651386380196, 0.014865261502563953, 0.0036516557447612286, -0.005537862423807383, 0.011143092066049576, 0.011356080882251263, -0.008499382995069027, -0.015915805473923683, 0.001295937690883875, 0.03629670664668083, 0.038976095616817474, 0.033428825438022614, 0.032993897795677185, 0.02841532602906227, 0.02451731264591217, 0.022545520216226578, 0.02724434807896614, 0.02838853746652603, 0.027603670954704285, 0.028169872239232063, 0.028828293085098267, 0.02889416739344597, 0.028997056186199188, 0.028535181656479836, 0.026410352438688278, 0.026367953047156334, 0.03121563233435154, 0.03588540852069855, 0.03955190256237984, 0.035150784999132156, 0.03051012195646763, 0.032955560833215714, 0.02942657098174095, 0.018812870606780052, 0.008279817178845406, 0.003518986515700817, 0.029718324542045593, 0.043439824134111404, 0.037052590399980545, 0.06802663952112198, 0.08852826058864594, 0.07487991452217102, 0.051908522844314575], atol=0.5))


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):

    def test_no_radius_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack")
                         +" "+"outnewfolder")



        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py")
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack")
                        +" "+"outoldfolder")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)

        print_new = a.stdout.decode('utf8').split('\n')
        print_old = b.stdout.decode('utf8').split('\n')

        print(print_new)
        # with patch.object(sys, 'argv', testargs_new):
        #     old_stdout = sys.stdout
        #     print_new = StringIO()
        #     sys.stdout = print_new
        #     return_new=fu.main(testargs_new[1:])
        # with patch.object(sys, 'argv', testargs_old):
        #     print_old = StringIO()
        #     sys.stdout = print_old
        #     return_old=oldfu.main(testargs_old[1:])
        # sys.stdout = old_stdout
        # self.assertEqual(return_new,return_old)
        # self.assertEqual(return_new, 1)
        # self.assertEqual(print_new.getvalue().split('\n')[1],' ==radius== mandatory option is missing.')
        # self.assertEqual(print_new.getvalue().split('\n')[1],print_old.getvalue().split('\n')[1])

    # todo: there is an error in li = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0] when I go throught it in the second 'main' call. it does not depend on the lib version. it seems to be an mpi syncronization situation ... ask markus
    @unittest.skip("exitcode 1")
    def test_negative_radius_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"), "outnewfolder", "--radius=-145"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"), "outoldfolder","--radius=-145"]
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
        testargs_new = (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack")
                        +" "+"outnewfolder"
                        +" "+"--radius=145"
                        +" "+"--minimum_grp_size=3"
                        +" "+"--img_per_grp=2")


        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack")
                        +" "+"outoldfolder"
                        +" "+"--radius=145"
                        +" "+"--minimum_grp_size=3"
                        +" "+"--img_per_grp=2")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)

        print_new = a.stdout.decode('utf8').split('\n')
        print_old = b.stdout.decode('utf8').split('\n')


        # with patch.object(sys, 'argv', testargs_new):
        #     old_stdout = sys.stdout
        #     print_new = StringIO()
        #     sys.stdout = print_new
        #     return_new=fu.main(testargs_new[1:])
        # with patch.object(sys, 'argv', testargs_old):
        #     print_old = StringIO()
        #     sys.stdout = print_old
        #     return_old=oldfu.main(testargs_old[1:])
        # sys.stdout = old_stdout
        # self.assertEqual(return_new,return_old)
        # self.assertEqual(return_new, 1)
        # self.assertEqual(print_new.getvalue().split('\n')[1],'ERROR! Minimum group size (3) is larger than the actual group size (2). Oh dear :(')
        # self.assertEqual(print_new.getvalue().split('\n')[1],print_old.getvalue().split('\n')[1])

    #todo: there is an error in li = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0] when I go throught it in the second 'main' call. it does not depend on the lib version. it seems to be an mpi syncronization situation ... ask markus
    @unittest.skip("exitcode 1")
    def test_CTF_and_VPP_error(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"), "outnewfolder", "--radius=145", "--CTF", "--VPP"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py"),'bdb:'+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"), "outoldfolder","--radius=145", "--CTF", "--VPP"]
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
