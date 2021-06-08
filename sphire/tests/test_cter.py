from __future__ import print_function
from __future__ import division


from time import sleep

from sphire.bin import sp_cter as fu
from sphire.bin_py3 import sp_cter as oldfu


from sphire.libpy.sp_utilities import get_im
from numpy import array_equal
from os import path,listdir
from .test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_list_of_file,remove_dir
import unittest

from sphire.libpy import sp_global_def

"""
WHAT IS MISSING:
0) we have to run it using mpirun
1) Due to internal randomization of the statistical resampling procedure used to compute tparameters and their 
    standard deviation, results of different runs will differ slightly.
    For this reason, after speaking with Tapu, I decided to check the image in the "micthumb" folder



RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) Test_run see what_missing_1


In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:
1) Test_Error_cases::test_negative_radius_error 
"""



try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch
except :
    pass

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO
import sys



""" see https://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/"""
class Test_Error_cases(unittest.TestCase):
    old_output_folder="CterOld"
    new_output_folder = "CterNew"

    def remove_folders(self):
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)

    @classmethod
    def tearDownClass(cls):
        remove_list_of_file([f for f in listdir(".") if "sp_cter_logfile" in f])

    def test_lowest_resolution_error(self):
        sp_global_def.BATCH = True
        self.remove_folders()
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc"),self.new_output_folder, '--apix=1.0', '--f_start=0.3' ]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc"), self.old_output_folder, '--apix=1.0','--f_start=0.3']
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        sp_global_def.BATCH= True
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => f_start should be in Angstrom')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])
        self.remove_folders()

    def test_highest_resolution_error(self):
        sp_global_def.BATCH = True
        self.remove_folders()
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc"), self.new_output_folder, '--apix=1.0', '--f_stop=0.3']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc"), self.old_output_folder, '--apix=1.0', '--f_stop=0.3']
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit):
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            with self.assertRaises(SystemExit):
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => f_stop should be in Angstrom')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])
        self.remove_folders()

    def test_too_few_params_error(self):
        sp_global_def.BATCH = True
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"CorrectedSums","corrsum_dw","TcdA1-001*_frames.mrc")]
        with patch.object(sys, 'argv', testargs_new):
            with self.assertRaises(SystemExit) as cnew:
                old_stdout = sys.stdout
                print_new = StringIO()
                sys.stdout = print_new
                fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old) :
            with self.assertRaises(SystemExit) as cold:
                print_old = StringIO()
                sys.stdout = print_old
                oldfu.main()
        sys.stdout = old_stdout

        self.assertEqual(str(cnew.exception),str(cold.exception))
        self.assertEqual(str(cnew.exception),"None")



class Test_run(unittest.TestCase):
    old_output_folder="CterOld"
    new_output_folder = "CterNew"
    filename= "micthumb/TcdA1-0010_frames_thumb.hdf"

    @classmethod
    def tearDownClass(cls):
        sleep(0.1)      #to give the script the time to close the last logfile.
        remove_list_of_file([f for f in listdir(".") if "sp_cter_logfile" in f])

    def remove_folders(self):
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)


    def test_cter_mrk(self):
        sp_global_def.BATCH = True
        self.remove_folders()
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.new_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.old_output_folder,"--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()
        old_value=get_im(path.join(self.old_output_folder,self.filename))
        new_value = get_im(path.join(self.new_output_folder, self.filename))
        self.assertTrue(array_equal(old_value.get_2dview(),new_value.get_2dview()))
        self.assertTrue(array_equal( new_value.get_2dview().flatten().tolist()[0:100],[-151.91610717773438, -5.15772008895874, -7.199570178985596, 16.844715118408203, -86.01467895507812, 224.5313720703125, 255.96446228027344, 127.28541564941406, -107.42439270019531, -146.97732543945312, 45.27460861206055, 134.89849853515625, 193.75506591796875, 1.6738426685333252, -291.7831726074219, -261.452880859375, -87.77507019042969, -253.7869873046875, -14.231243133544922, 105.49187469482422, 54.15253448486328, 202.34719848632812, -145.9459686279297, -310.8766784667969, -168.98582458496094, -246.27511596679688, -267.9120788574219, -83.6864013671875, 277.79608154296875, 151.8765106201172, 13.589789390563965, 11.07779312133789, -84.60417175292969, 100.78224182128906, -182.1418914794922, -1.2914609909057617, -105.66407775878906, -8.37000560760498, 96.46575164794922, -81.93226623535156, -29.234130859375, -128.64105224609375, 80.63739776611328, -11.062923431396484, 107.84751892089844, -108.09345245361328, 17.696636199951172, -93.53311157226562, 67.41581726074219, 273.3633728027344, -28.720521926879883, -34.706932067871094, -195.11680603027344, -200.7438507080078, -46.45978546142578, -22.808040618896484, -91.24805450439453, 113.92597961425781, 115.37574005126953, -19.18140411376953, -39.19940948486328, -80.63643646240234, 291.9253845214844, 32.61348342895508, -112.9101333618164, -130.7708282470703, -70.37080383300781, -260.5411071777344, -163.13958740234375, -141.8155975341797, -4.3243408203125, -79.0434341430664, -118.28507232666016, 245.7770233154297, -86.85401153564453, -84.1547622680664, -5.078440189361572, -53.13208770751953, -47.206058502197266, 53.22468948364258, 75.35106658935547, 77.2694320678711, 113.52538299560547, -53.97866439819336, 1.6344045400619507, 57.784603118896484, -78.45945739746094, -69.10759735107422, -83.55438232421875, 63.753501892089844, 75.83992767333984, -36.30785369873047, 97.44556427001953, 100.17158508300781, 57.608177185058594, -89.75988006591797, 141.65589904785156, 14.772521018981934, -143.1587371826172, -54.35596466064453]))
        self.remove_folders()


    def test_cter_vpp(self):
        sp_global_def.BATCH = True
        self.remove_folders()
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.new_output_folder+"vpp","--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--vpp", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.old_output_folder+"vpp","--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--vpp", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()

        old_value = get_im(path.join(self.old_output_folder+"vpp", self.filename))
        new_value = get_im(path.join(self.new_output_folder+"vpp", self.filename))
        self.assertTrue(array_equal(old_value.get_2dview(), new_value.get_2dview()))
        self.assertTrue(array_equal(new_value.get_2dview().flatten().tolist()[0:100],[-151.91610717773438, -5.15772008895874, -7.199570178985596, 16.844715118408203, -86.01467895507812, 224.5313720703125, 255.96446228027344, 127.28541564941406, -107.42439270019531, -146.97732543945312, 45.27460861206055, 134.89849853515625, 193.75506591796875, 1.6738426685333252, -291.7831726074219, -261.452880859375, -87.77507019042969, -253.7869873046875, -14.231243133544922, 105.49187469482422, 54.15253448486328, 202.34719848632812, -145.9459686279297, -310.8766784667969, -168.98582458496094, -246.27511596679688, -267.9120788574219, -83.6864013671875, 277.79608154296875, 151.8765106201172, 13.589789390563965, 11.07779312133789, -84.60417175292969, 100.78224182128906, -182.1418914794922, -1.2914609909057617, -105.66407775878906, -8.37000560760498, 96.46575164794922, -81.93226623535156, -29.234130859375, -128.64105224609375, 80.63739776611328, -11.062923431396484, 107.84751892089844, -108.09345245361328, 17.696636199951172, -93.53311157226562, 67.41581726074219, 273.3633728027344, -28.720521926879883, -34.706932067871094, -195.11680603027344, -200.7438507080078, -46.45978546142578, -22.808040618896484, -91.24805450439453, 113.92597961425781, 115.37574005126953, -19.18140411376953, -39.19940948486328, -80.63643646240234, 291.9253845214844, 32.61348342895508, -112.9101333618164, -130.7708282470703, -70.37080383300781, -260.5411071777344, -163.13958740234375, -141.8155975341797, -4.3243408203125, -79.0434341430664, -118.28507232666016, 245.7770233154297, -86.85401153564453, -84.1547622680664, -5.078440189361572, -53.13208770751953, -47.206058502197266, 53.22468948364258, 75.35106658935547, 77.2694320678711, 113.52538299560547, -53.97866439819336, 1.6344045400619507, 57.784603118896484, -78.45945739746094, -69.10759735107422, -83.55438232421875, 63.753501892089844, 75.83992767333984, -36.30785369873047, 97.44556427001953, 100.17158508300781, 57.608177185058594, -89.75988006591797, 141.65589904785156, 14.772521018981934, -143.1587371826172, -54.35596466064453]))
        self.remove_folders()

    def test_cter_pap(self):
        sp_global_def.BATCH = True
        self.remove_folders()
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.new_output_folder+"pap","--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--pap", "--f_start=40", "--f_stop=34"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_cter.py"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "CorrectedSums", "corrsum_dw","TcdA1-001*_frames.mrc"),self.old_output_folder+"pap","--selection_list="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "01_CTER","Tutorial_micrographs_select.txt"), "--apix=1.14", "--Cs=0", "--pap", "--f_start=40", "--f_stop=34"]
        with patch.object(sys, 'argv', testargs_new):
            fu.main()
        sp_global_def.BATCH = True
        with patch.object(sys, 'argv', testargs_old):
            oldfu.main()

        old_value = get_im(path.join(self.old_output_folder+"pap", self.filename))
        new_value = get_im(path.join(self.new_output_folder+"pap", self.filename))
        self.assertTrue(array_equal(old_value.get_2dview(), new_value.get_2dview()))
        self.assertTrue(array_equal(new_value.get_2dview().flatten().tolist()[0:100], [-151.91610717773438, -5.15772008895874, -7.199570178985596, 16.844715118408203, -86.01467895507812, 224.5313720703125, 255.96446228027344, 127.28541564941406, -107.42439270019531, -146.97732543945312, 45.27460861206055, 134.89849853515625, 193.75506591796875, 1.6738426685333252, -291.7831726074219, -261.452880859375, -87.77507019042969, -253.7869873046875, -14.231243133544922, 105.49187469482422, 54.15253448486328, 202.34719848632812, -145.9459686279297, -310.8766784667969, -168.98582458496094, -246.27511596679688, -267.9120788574219, -83.6864013671875, 277.79608154296875, 151.8765106201172, 13.589789390563965, 11.07779312133789, -84.60417175292969, 100.78224182128906, -182.1418914794922, -1.2914609909057617, -105.66407775878906, -8.37000560760498, 96.46575164794922, -81.93226623535156, -29.234130859375, -128.64105224609375, 80.63739776611328, -11.062923431396484, 107.84751892089844, -108.09345245361328, 17.696636199951172, -93.53311157226562, 67.41581726074219, 273.3633728027344, -28.720521926879883, -34.706932067871094, -195.11680603027344, -200.7438507080078, -46.45978546142578, -22.808040618896484, -91.24805450439453, 113.92597961425781, 115.37574005126953, -19.18140411376953, -39.19940948486328, -80.63643646240234, 291.9253845214844, 32.61348342895508, -112.9101333618164, -130.7708282470703, -70.37080383300781, -260.5411071777344, -163.13958740234375, -141.8155975341797, -4.3243408203125, -79.0434341430664, -118.28507232666016, 245.7770233154297, -86.85401153564453, -84.1547622680664, -5.078440189361572, -53.13208770751953, -47.206058502197266, 53.22468948364258, 75.35106658935547, 77.2694320678711, 113.52538299560547, -53.97866439819336, 1.6344045400619507, 57.784603118896484, -78.45945739746094, -69.10759735107422, -83.55438232421875, 63.753501892089844, 75.83992767333984, -36.30785369873047, 97.44556427001953, 100.17158508300781, 57.608177185058594, -89.75988006591797, 141.65589904785156, 14.772521018981934, -143.1587371826172, -54.35596466064453]))
        self.remove_folders()
