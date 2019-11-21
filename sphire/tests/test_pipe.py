from __future__ import print_function
from __future__ import division


from numpy import array_equal

from sphire.bin import sp_pipe as oldfu
from sphire.utils.SPHIRE.bin import sp_pipe as fu
from os import path
from test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir
import unittest
from sphire.libpy.sp_utilities import get_im
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack','bdb:' + "Bad_path_toDB",path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC") ,'new_output_folder']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack','bdb:' + "Bad_path_toDB" ,path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC") , 'old_output_folder']
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack', "bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),"output_dir_not_existing",'new_output_folder']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"),'isac_substack', "bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),"output_dir_not_existing",'old_output_folder']
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"),'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC") ,' /']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC") ,' /']
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),"/",'new_output_folder']
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),"/",'new_output_folder']
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),'NEW_OUTPUT',"--isac_class_avgs_path="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"invalidDB.hdf")]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),'old_outpuT',"--isac_class_avgs_path="+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"invalidDB.hdf")]
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
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),'NEW_OUTPUT',"--substack_basename= "]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),'OLD_OUTPUT',"--substack_basename= "]
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


class Test_isac_substack(unittest.TestCase):
    old_output_folder="substackOld"
    new_output_folder = "substackNew"
    def test_run(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),self.new_output_folder]#," --isac_class_avgs_path ='/home/lusnig/Downloads/LucaTest_beaut/best.hdf'"]
        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_pipe.py"), 'isac_substack',"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"03_PARTICLES#stack"),path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "04_ISAC"),self.old_output_folder]
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
        old_value =get_im(stackname="bdb:"+self.old_output_folder+"#isac_substack",im=0)
        new_value = get_im(stackname="bdb:" + self.new_output_folder + "#isac_substack", im=0)
        self.assertTrue(array_equal(old_value.get_2dview(),new_value.get_2dview()))
        self.assertTrue(array_equal( old_value.get_2dview().flatten().tolist()[0:100], [-0.7735422849655151, 0.04221929982304573, 0.9302629828453064, 1.15420401096344, 0.5812007784843445, 0.38344860076904297, 0.2833080291748047, 0.5844124555587769, 0.59711754322052, 0.8901447057723999, 1.249489665031433, 1.6613267660140991, 1.9065828323364258, 1.6155502796173096, 1.8675326108932495, 1.7501505613327026, 1.1296839714050293, 0.5398319959640503, -0.26309671998023987, -0.4878522455692291, -0.9763734340667725, -1.1236287355422974, -1.4459232091903687, -1.1048630475997925, -0.25124940276145935, -0.5049065351486206, -0.4947991967201233, -0.444855272769928, -0.12026157975196838, 0.33534929156303406, -0.03821266070008278, -0.6292971968650818, -0.9841951727867126, -0.5810408592224121, -0.11035419255495071, -0.016610929742455482, 0.6521958708763123, 0.7927666306495667, 1.0315550565719604, 1.305680274963379, 1.018584966659546, 1.4733165502548218, 0.7437952756881714, 0.2895500063896179, 0.5073809623718262, 0.37767210602760315, 0.2228148877620697, 0.4134519398212433, 0.2176506370306015, 0.9907333850860596, 0.9602268934249878, 1.1047149896621704, 1.0669927597045898, 1.2039128541946411, 1.2368600368499756, 1.0187309980392456, 0.058885738253593445, -0.0037675732746720314, -0.19165875017642975, -0.2862400412559509, -0.39332082867622375, -0.37311962246894836, -0.46995630860328674, -0.6182337403297424, -0.18015365302562714, 0.07094576954841614, 0.40471479296684265, 0.21657632291316986, 0.7946628928184509, 1.6166480779647827, 1.65879487991333, 0.9631286263465881, 0.8932328820228577, 0.3658081591129303, -0.19450199604034424, -1.180352807044983, -1.0464118719100952, -0.967138946056366, -0.6414624452590942, -0.18117772042751312, -0.020910970866680145, 0.2834717631340027, 0.7001526355743408, 1.42698073387146, 1.5674123764038086, 1.5099875926971436, 1.6790717840194702, 1.8142207860946655, 2.0849273204803467, 1.8612970113754272, 1.4969502687454224, 0.930846631526947, -0.12830914556980133, -0.18110060691833496, -0.6524917483329773, -0.4672186076641083, -0.3151601254940033, -0.050327446311712265, 0.022716891020536423, 0.6298686861991882]))
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)