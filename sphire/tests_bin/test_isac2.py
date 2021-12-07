from __future__ import print_function
from __future__ import division


import subprocess

from numpy import allclose
import shutil
from sphire.libpy.sp_utilities import get_im
from os import path
from tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_RESOURCES,ABSOLUTE_BIN_PATH,remove_dir
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


MPI_PATH = shutil.which("mpirun")
NUM_PROC = 8
TOLERANCE =2.5

@unittest.skip("nov_21 IT FAILS because invalid MPI_PATH ")
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
            +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_Particles#stack'")
            +" '"+self.new_output_folder+"' "
            +" --radius=145"
            +" --CTF")

        testcommand_old= (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_isac2.py")
            +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_Particles#stack'")
            +" '"+self.old_output_folder+"' "
            +" --radius=145"
            +" --CTF")

        subprocess.run(args =[testcommand_new], shell=True, capture_output=True)
        subprocess.run(args=[testcommand_old], shell=True, capture_output=True)

        old_value = get_im(path.join(self.old_output_folder, self.filename))
        new_value = get_im(path.join(self.new_output_folder, self.filename))

        self.assertTrue(allclose(old_value.get_2dview(),new_value.get_2dview(), atol=TOLERANCE))
        self.assertTrue(allclose( new_value.get_2dview().flatten().tolist()[0:100], [0.020844653248786926, 0.022487156093120575, 0.02208145707845688, 0.018550248816609383, 0.018665246665477753, 0.019494876265525818, 0.02418450079858303, 0.027667080983519554, 0.02900763787329197, 0.025944046676158905, 0.031390730291604996, 0.027120821177959442, 0.00910606887191534, -0.0004518263158388436, 0.009620467200875282, 0.025462428107857704, 0.02613203041255474, 0.026837026700377464, 0.019651882350444794, -0.0013422714546322823, -0.004115818999707699, 0.03214756399393082, 0.03178049996495247, 0.011443180963397026, 0.025857364758849144, 0.026638126000761986, 0.03655451536178589, 0.05426323413848877, 0.06268023699522018, 0.03984970971941948, 0.013704168610274792, 0.008193270303308964, 0.005416750907897949, -0.012422064319252968, -0.012422808445990086, 0.014235804788768291, 0.024870628491044044, -0.005099047441035509, -0.009270750917494297, 0.015586614608764648, -0.0036001892294734716, -0.02533501200377941, -0.004470369778573513, 0.02382775768637657, 0.05031612142920494, 0.04690878838300705, 0.02393484115600586, 0.027014095336198807, 0.039149295538663864, 0.03573647513985634, 0.02213861420750618, 0.032027628272771835, 0.03532850742340088, 0.03128813952207565, 0.01826963573694229, 0.028490642085671425, 0.046303264796733856, 0.01782909594476223, -0.00606165686622262, -0.016502195969223976, -0.005856842268258333, 0.019041359424591064, 0.02800695225596428, 0.026184026151895523, 0.023107700049877167, 0.02628708817064762, 0.022019468247890472, 0.019245078787207603, 0.020709989592432976, 0.021585211157798767, 0.02417927235364914, 0.021674703806638718, 0.0199696384370327, 0.021289950236678123, 0.020935336127877235, 0.021501578390598297, 0.022525610402226448, 0.02044103667140007, 0.01824958808720112, 0.023783020675182343, 0.025899477303028107, 0.023173417896032333, 0.0289889145642519, 0.03435958921909332, 0.03157936781644821, 0.03505755960941315, 0.04005855321884155, 0.034827277064323425, 0.022344911471009254, 0.008706190623342991, 0.01249337662011385, 0.019858719781041145, 0.021026115864515305, 0.029216332361102104, 0.010308701545000076, -0.01953706704080105, -0.012493203394114971, 0.022878523916006088, 0.022390954196453094, 0.0009752019541338086], atol=TOLERANCE))


""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
@unittest.skip("nov_21 IT FAILS in local because missing file ... it could be because my config")
class Test_Error_cases(unittest.TestCase):

    def test_no_radius_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                         +" "+"outnewfolder")



        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py")
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                        +" "+"outoldfolder")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)

        self.assertTrue('==radius== mandatory option is missing' in b.stdout.decode('utf8'))
        self.assertTrue('==radius== mandatory option is missing' in a.stdout.decode('utf8'))

    def test_negative_radius_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                         +" "+"outnewfolder"
                         +" " + "--radius=-145")
        testargs_old =  (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                         +" "+"outoldfolder"
                         +" " + "--radius=-145")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)


        self.assertTrue('ERROR => Particle radius has to be provided!' in b.stdout.decode('utf8'))
        self.assertTrue('ERROR => Particle radius has to be provided!' in a.stdout.decode('utf8'))

        self.assertTrue('li = mpi.mpi_bcast(li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD)' in b.stderr.decode('utf8'))
        self.assertTrue('li = mpi.mpi_bcast(li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD)' in a.stderr.decode('utf8'))


    def test_minimum_group_size_higher_group_size_error(self):
        testargs_new = (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                        +" "+"outnewfolder"
                        +" "+"--radius=145"
                        +" "+"--minimum_grp_size=3"
                        +" "+"--img_per_grp=2")


        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py")
                        +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                        +" "+"outoldfolder"
                        +" "+"--radius=145"
                        +" "+"--minimum_grp_size=3"
                        +" "+"--img_per_grp=2")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)
        self.assertTrue('ERROR! Minimum group size (3) is larger than the actual group size (2)' in b.stdout.decode('utf8'))
        self.assertTrue('ERROR! Minimum group size (3) is larger than the actual group size (2)' in a.stdout.decode('utf8'))

    def test_CTF_and_VPP_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                         +" "+"outnewfolder"
                         +" " + "--radius=145"
                         +" " + "--CTF"
                         +" " + "--VPP")
        testargs_old =  (path.join(ABSOLUTE_OLDBIN_PATH, "sp_isac2.py")
                         +" "+'bdb:'+path.join(ABSOLUTE_PATH_TO_RESOURCES,"03_PARTICLES#stack")
                         +" " + "--radius=145"
                         +" " + "--CTF"
                         +" " + "--VPP")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)

        self.assertTrue('Options CTF and VPP cannot be used together' in b.stdout.decode('utf8'))
        self.assertTrue('Options CTF and VPP cannot be used together' in a.stdout.decode('utf8'))