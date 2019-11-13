from __future__ import print_function
from __future__ import division




from sphire.bin import sp_pipe as oldfu
from sphire.utils.SPHIRE.bin import sp_pipe as fu
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


"""
WHAT IS MISSING:
0) get_time_stamp_suffix --> returns a timestamp ... not testable
1) Test_isac_substack cannot run it under pycharm

RESULT AND KNOWN ISSUES
1) 

In these tests there is a bug --> syntax error:
1)

In these tests there is a strange behavior:
"""


class Test_helperFunctions(unittest.TestCase):
    def test_get_cmd_line(self):
        cmdLine=["this", "is", "a", "test"]
        with patch.object(sys, 'argv', cmdLine):
            return_new = fu.get_cmd_line()
            return_old = oldfu.get_cmd_line()
            self.assertEqual(return_new,return_old)
            self.assertEqual(return_new, 'Shell line command: this  is  a  test  ')

class Test_Error_cases_isac_substack(unittest.TestCase):
    def test_bdb_not_existing(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack','bdb:' + "Bad_path_toDB",path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut") ,' /home/lusnig/Downloads/LucaTest2substack']#," --isac_class_avgs_path ='/home/lusnig/Downloads/LucaTest_beaut/best.hdf'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + "Bad_path_toDB" ,path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut") , '/home/lusnig/Downloads/LucaTest2substackold']
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
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => Input BDB image stack file does not exist. Please check the file path and restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])

    def test_run_output_dir_not_existing(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),"output_dir_not_existing",' /home/lusnig/Downloads/LucaTest2substack']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),"output_dir_not_existing", '/home/lusnig/Downloads/LucaTest2substackold']
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
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => ISAC or Beautifier run output directory does not exist. Please check the directory path and restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])

    def test_output_dir_existing(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut") ,' /']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut") , '/']
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
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => Output directory exists. Please change the name and restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])

    def test_default_ordered_class_averagesHDF_not_found(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),"/",' /home/lusnig/Downloads/LucaTedsdsst2substack']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),"/",'/home/lusnig/Downloads/LucaTedxdsst2substackold']
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

        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => ISAC or Beautifier run output directory does not contain the default ISAC class average stack file (/ordered_class_averages.hdf). Please check the directory path or specify ISAC class average stack file, then restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])

    def test_specified_ordered_class_averagesHDF_not_found(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),' /home/lusnig/Downloads/LucaTedsdsst2substack',"--isac_class_avgs_path="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"invalidDB.hdf")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),'/home/lusnig/Downloads/LucaTedxdsst2substackold',"--isac_class_avgs_path="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,"invalidDB.hdf")]
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
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => The specifed ISAC class average stack file does not exist. Please check the file path and restart the program.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])



    def test_substack_basename_empty_string(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),' /home/lusnig/Downloads/LucaTedsdsst2substack',"--substack_basename= "]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),'/home/lusnig/Downloads/LucaTedxdsst2substackold',"--substack_basename= "]
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
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],' => Substack basename cannot be empty string or only white spaces.')
        self.assertEqual(print_new.getvalue().split('\n')[5].split("ERROR")[1],print_old.getvalue().split('\n')[5].split("ERROR")[1])

# At line 859 it tries to run this cmd line but it can find, because i'm running it from pycharm the e2dbd.py files
# cmd_line'e2bdb.py bdb:/home/lusnig/Downloads/SphireDemoResults/Particles#stack --makevstack=bdb:/home/lusnig/Downloads/LucaTedsdsst2substack#isac_substack --list=/home/lusnig/Downloads/LucaTedsdsst2substack/isac_substack_particle_id_list.txt'
@unittest.skip("cannot run it from pycharm")
class Test_isac_substack(unittest.TestCase):
    """ see pag 42 of the sphire1.3 tutorial"""
    def test_(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),' /home/lusnig/Downloads/LucaTedsdsst2substack']#," --isac_class_avgs_path ='/home/lusnig/Downloads/LucaTest_beaut/best.hdf'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Particles#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Beaut"),'/home/lusnig/Downloads/LucaTedxdsst2substackold']
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