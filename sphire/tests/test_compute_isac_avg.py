from __future__ import print_function
from __future__ import division


from sphire.bin import sp_compute_isac_avg as oldfu
from sphire.utils.SPHIRE.bin import sp_compute_isac_avg as fu
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

    #there is not a test, just the right syntax to run it. The idea is to create the data running isac2 over Fabian's tiny dataset.
    def test_(self):
        testargs_new =  [path.join(ABSOLUTE_BIN_PATH, "sp_compute_isac_avg.py"), "--stack=bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"), "--isac_dir="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Class2D"),"--output_dir=BeautLuca","--pixel_size=-1.0","--local_alignment"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_compute_isac_avg.py"),"--stack=bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Particles#stack"),"--isac_dir="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"Class2D"),"--output_dir=BeautLuca2","--pixel_size=-1.0","--local_alignment"]
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
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],' => f_start should be in Angstrom')
        self.assertEqual(print_new.getvalue().split('\n')[1].split("ERROR")[1],print_old.getvalue().split('\n')[1].split("ERROR")[1])
