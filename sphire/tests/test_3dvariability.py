from __future__ import print_function
from __future__ import division

from numpy import allclose,array_equal

from os import path
# from sphire.tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_BIN_PATH,remove_dir,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW
import unittest
from sphire.libpy.sp_utilities import get_im
import EMAN2db
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


"""
WHAT IS MISSING:
1) I am not able to reproduce the windows negative size error (see ln279)

RESULT AND KNOWN ISSUES:
The "SymStack" folder used at pag 73 of the tutorial1_3 is not present in the new "sphireDemoResults" folder. If you
    copy it from the old one,that you can download from the homepage (25/11/19), some data miss and it crash
    Maybe We could generate these values by 'test_symmetrize'. In case we will do that when we generalize all the tests


In these tests there is a bug --> syntax error:
-) Test_helperFunctions::test_VAR_VAR3D_error if you activate the flag "VAR" you will have a crash because this bug:
    UnboundLocalError: local variable 'index_of_proj' referenced before assignment (ln:572)
-) Test_run::test_ in the tutorial (pg73) example crashes in both the libs at ln (690 or 546) with the following error:
    UnboundLocalError: local variable 'dummy' referenced before assignment
    You have to run it using the mpirun but in the tutorial is not specified and the script does not spawn an error 
    at the beginning of the code
    


In these tests there is a strange behavior:
-) 'Test_run::test_symmetrize' it cannot recognize the db of the precalculated results because in the new tutorial 
    the folder name's are renamed. Replacing "03_PARTICLES" in "03_Particles" it works
"""
MPI_PATH = "/home/lusnig/SPHIRE_1_1/envs/sphire_py3transition/bin/mpirun" #"/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/"
NUM_PROC = 8
class Test_helperFunctions(unittest.TestCase):
    old_output_folder = "3DvariabilitySimmetryOld"
    new_output_folder = "3DvariabilitySimmetryNew"

    def test_c1Sym_error_new(self):
        testargs_new = (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" " +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+ self.new_output_folder
                        +" "+ "--sym=c1")
        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" " +" 'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+ self.old_output_folder
                        +" "+ "--sym=c1")

        a = subprocess.run(args =[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args =[testargs_old], shell=True, capture_output=True)

        self.assertTrue(' => There is no need to symmetrize stack for C1 symmetry' in a.stdout.decode('utf8'))
        self.assertTrue(' => There is no need to symmetrize stack for C1 symmetry' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)

    def test_c1Sym_error_old(self):
        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" "+" 'bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir=" + self.old_output_folder
                        +" "+"--sym=c1")
        testargs_new = (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" "+" 'bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir=" + self.new_output_folder
                        +" "+"--sym=c1")
        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True, capture_output=True)


        self.assertTrue(' => There is no need to symmetrize stack for C1 symmetry' in b.stdout.decode('utf8'))
        self.assertTrue(' => There is no need to symmetrize stack for C1 symmetry' in a.stdout.decode('utf8'))
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)

    def test_inputStack_symmetrize_MPI_error(self):
        testargs_new = (MPI_PATH
                        + " -np "
                        +str(NUM_PROC)
                        +" "+path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" "+" 'bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir=" + self.new_output_folder
                        +" "+"--sym=c5")
        testargs_old = (MPI_PATH
                        + " -np "
                        +str(NUM_PROC)
                        +" "+path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        + " " + "--symmetrize"
                        +" "+" 'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+self.old_output_folder
                        +" "+"--sym=c5")

        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)



        self.assertTrue(' => Cannot use more than one CPU for symmetry preparation' in a.stdout.decode('utf8'))
        self.assertTrue(' => Cannot use more than one CPU for symmetry preparation' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)

    def test_LowPass_param_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                         +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                         +" "+"--output_dir="+self.new_output_folder
                         +" "+"--sym=c5"
                         +" --fl=0.1"
                         +" --aa=0")
        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+self.old_output_folder
                        +" "+"--sym=c5"
                        +" --fl=0.1"
                        +" --aa=0")

        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)

        self.assertTrue(' => Fall off has to be given for the low-pass filter' in a.stdout.decode('utf8'))
        self.assertTrue(' => Fall off has to be given for the low-pass filter' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)


    def test_VAR_VAR2D_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                         +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                         +" "+"--output_dir="+self.new_output_folder
                         +" "+"--var2D='given_file'"
                         +" "+"--VAR")
        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+self.old_output_folder
                        +" "+"--var2D='given_file'"
                        +" "+"--VAR")
        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)

        self.assertTrue(' => When VAR is set, the program cannot output ave2D, ave3D or var2D' in a.stdout.decode('utf8'))
        self.assertTrue(' => When VAR is set, the program cannot output ave2D, ave3D or var2D' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)


    def test_decimate_higher1_error(self):
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH,"sp_3dvariability.py")
                         +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                         +" "+"--output_dir="+self.new_output_folder
                         +" "+"--decimate=1.1")

        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+self.old_output_folder
                        +" "+"--decimate=1.1")

        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)

        self.assertTrue(' => Decimate rate should be a value between 0.0 and 1.0' in a.stdout.decode('utf8'))
        self.assertTrue(' => Decimate rate should be a value between 0.0 and 1.0' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)


    def test_decimate_lower0_error(self):
        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                         +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                         +" "+"--output_dir="+self.new_output_folder
                         +" "+"--decimate=-1.1")

        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+self.old_output_folder
                        +" "+"--decimate=-1.1")

        a = subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        b = subprocess.run(args=[testargs_old], shell=True,  capture_output=True)

        self.assertTrue(' => Decimate rate should be a value between 0.0 and 1.0' in a.stdout.decode('utf8'))
        self.assertTrue(' => Decimate rate should be a value between 0.0 and 1.0' in b.stdout.decode('utf8'))

        remove_dir(self.new_output_folder)
        remove_dir(self.old_output_folder)




class Test_run(unittest.TestCase):

    def test_VAR_VAR3D(self):
        new_output_folder = "3DvariabilitySimmetry_VAR_new"
        old_output_folder = "3DvariabilitySimmetry_VAR_Old"
        testargs_new = (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                        + " " + "'bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,
                                                    "06_SUBSTACK_ANO#isac_substack'")
                        + " " + "--output_dir=" + new_output_folder
                        + " " + "--var3D='given_file.hdf'"
                        + " " + "--VAR")

        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        + " " + "'bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,
                                                    "06_SUBSTACK_ANO#isac_substack'")
                        + " " + "--output_dir=" + old_output_folder
                        + " " + "--var3D='given_file.hdf'"
                        + " " + "--VAR")

        subprocess.run(args=[testargs_new], shell=True, capture_output=True)
        subprocess.run(args=[testargs_old], shell=True, capture_output=True)

        newimg = EMAN2db.EMData()
        oldimg = EMAN2db.EMData()
        newimg.read_image(fsp=path.join(new_output_folder, 'given_file.hdf'))
        oldimg.read_image(fsp=path.join(old_output_folder, 'given_file.hdf'))

        self.assertTrue(array_equal(newimg.get_3dview().flatten().tolist()[0:100],[-3.398010085220449e-05, -2.7903670343221165e-05, -2.64776263065869e-05, -3.099508830928244e-05, -4.180952601018362e-05, -5.889337626285851e-05, -8.216228889068589e-05, -0.0001114098631660454, -0.00014643699978478253, -0.00018704932881519198, -0.0002329741109861061, -0.0002830602170433849, -0.00033638920285739005, -0.0003925539494957775, -0.00045063687139190733, -0.0005096903187222779, -0.0005686617805622518, -0.0006270494195632637, -0.0006828202167525887, -0.0007355636917054653, -0.0007843654602766037, -0.0008282316266559064, -0.0008665817440487444, -0.0008984902524389327, -0.0009236905025318265, -0.0009416076354682446, -0.0009520584717392921, -0.0009527443908154964, -0.0009448318160139024, -0.0009287822176702321, -0.0009044426842592657, -0.0008719016332179308, -0.0008315170998685062, -0.0007839857717044652, -0.0007296598632819951, -0.000669682864099741, -0.0006020283326506615, -0.000529955665115267, -0.00045405683340504766, -0.0003751835902221501, -0.00029369888943620026, -0.00021090128575451672, -0.00012747617438435555, -4.43613862444181e-05, 3.792068673647009e-05, 0.00011861506209243089, 0.0001968519063666463, 0.000272340519586578, 0.00034437805879861116, 0.00041267750202678144, 0.00047643273137509823, 0.0005357097834348679, 0.0005901738186366856, 0.000641072925645858, 0.0006845941534265876, 0.0007241867715492845, 0.0007594855851493776, 0.0007906834944151342, 0.00081777130253613, 0.0008409013389609754, 0.0008603507303632796, 0.0008765358361415565, 0.0008905157446861267, 0.0009002792648971081, 0.0009084320045076311, 0.000915086769964546, 0.0009204641683027148, 0.0009247256675735116, 0.0009280719677917659, 0.0009306868887506425, 0.0009327546576969326, 0.0009340904070995748, 0.0009345625294372439, 0.0009347845334559679, 0.0009346535662189126, 0.000934049254283309, 0.0009328546584583819, 0.000930961687117815, 0.0009282570681534708, 0.0009246706031262875, 0.0009183128131553531, 0.0009107337682507932, 0.0009011240326799452, 0.0008893080521374941, 0.0008750110282562673, 0.0008581893052905798, 0.0008386565605178475, 0.0008164481841959059, 0.0007915245951153338, 0.0007615668582729995, 0.0007291691144928336, 0.0006934280972927809, 0.0006543724448420107, 0.0006121385376900434, 0.0005671223625540733, 0.0005194474360905588, 0.0004693828523159027, 0.0004160360258538276, 0.0003614358720369637, 0.00030481320573017]))
        self.assertTrue(array_equal(newimg.get_3dview(), oldimg.get_3dview()))

        remove_dir(new_output_folder)
        remove_dir(old_output_folder)

    def test_symmetrize(self):
        old_output_folder = "3DvariabilitySimmetry_sym_Old"
        new_output_folder = "3DvariabilitySimmetry_sym_New"
        remove_dir(new_output_folder)
        remove_dir(old_output_folder)
        filename = "sdata"
        testargs_new =  (path.join(ABSOLUTE_BIN_PATH, "sp_3dvariability.py")
                         +" "+"--symmetrize"
                         +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                         +" "+"--output_dir="+new_output_folder
                         +" "+"--sym=c5")
        testargs_old = (path.join(ABSOLUTE_OLDBIN_PATH, "sp_3dvariability.py")
                        +" "+"--symmetrize"
                        +" "+"'bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack'")
                        +" "+"--output_dir="+old_output_folder
                        +" "+"--sym=c5")

        subprocess.run(args=[testargs_new], shell=True,capture_output=True)
        subprocess.run(args=[testargs_old], shell=True,  capture_output=True)


        return_old = get_im('bdb:{0}#{1}'.format(old_output_folder, filename))
        return_new = get_im('bdb:{0}#{1}'.format(new_output_folder, filename))
        self.assertTrue(allclose(return_new.get_2dview(), return_old.get_2dview(), atol=0.1))
        self.assertTrue(allclose(return_new.get_2dview().flatten()[0:100], [-0.38441434502601624, -0.4573197066783905, -0.39960816502571106, -0.10431570559740067, 0.018583765253424644, 0.37694355845451355, 0.4353419840335846, 0.06498070061206818, -0.38467326760292053, -0.7100721001625061, -0.7571740746498108, -0.9328497648239136, -0.8550546169281006, -0.21433517336845398, 0.34148913621902466, 0.30916473269462585, 0.39394184947013855, 0.4447155296802521, -0.032833874225616455, 0.4121330976486206, 0.638403594493866, 0.34945985674858093, 0.2751978039741516, -0.14872148633003235, -0.5307353138923645, -0.6574574112892151, -0.8945914506912231, -1.107301950454712, -1.3676011562347412, -1.2225056886672974, -1.2273787260055542, -0.7121877074241638, 0.06733409315347672, 0.29879996180534363, 0.6046097874641418, 0.496036559343338, 0.7214235663414001, 0.7652679681777954, 0.9319247603416443, 0.11259084939956665, -0.3054805397987366, -0.8183746933937073, -1.4462037086486816, -1.4775571823120117, -1.2116808891296387, -0.911469042301178, -0.9934141039848328, -0.9891195297241211, -1.0073580741882324, -1.234678864479065, -1.4292954206466675, -1.5580313205718994, -1.446936845779419, -1.67134428024292, -2.0199873447418213, -2.0543317794799805, -2.2354137897491455, -1.99867844581604, -1.4405670166015625, -1.7088592052459717, -2.060110569000244, -1.6585056781768799, -1.0386717319488525, -0.4783007502555847, -0.12102964520454407, 0.16874085366725922, -0.12371158599853516, 0.01841016113758087, -0.6442297101020813, -0.6972493529319763, -0.1051267459988594, 0.016293810680508614, -0.10689342767000198, -0.472159206867218, -0.4886413514614105, -0.09828135371208191, -0.7409976720809937, -0.8893253207206726, -1.1477444171905518, -0.35942375659942627, 0.030692437663674355, 0.5380961298942566, 1.06252121925354, 0.6961002349853516, 0.7563707232475281, 1.0197871923446655, 0.9021824598312378, 0.34891727566719055, 0.402149498462677, 0.634938657283783, 0.6044892072677612, 0.24482420086860657, -0.4132026731967926, -0.9945327043533325, -1.260424256324768, -1.6913681030273438, -1.7999876737594604, -1.5891680717468262, -1.3196191787719727, -0.9893324375152588], atol=0.1))
        remove_dir(new_output_folder)
        remove_dir(old_output_folder)

    #pg73 tutorial
    @unittest.skip("need valid data")
    def test_(self):
        old_output_folder = "3DvarOld"
        new_output_folder = "3DvarNew"
        filename3DvarNew = "3DvarNew.hdf"
        filename3DavgNew = "3DAvgNew.hdf"
        filename3DvarOld = "3DvarOld.hdf"
        filename3DavgOld = "3DAvgOld.hdf"

        remove_dir(new_output_folder)
        remove_dir(old_output_folder)

        testargs_new = (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_OLDBIN_PATH,"sp_3dvariability.py")
            +" "+"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack")
            +" "+ "--output_dir="+new_output_folder
            +" --var3D="+filename3DvarNew
            + " --ave3D="+filename3DavgNew
            + " --sym=c5"
            +" --CTF"
            + " --var2D=''"
            + " --ave2D=''"
            + " --window=290"
        )

        testargs_old = (
            MPI_PATH
            + " -np "
            +str(NUM_PROC)
            +" "+path.join(ABSOLUTE_BIN_PATH,"sp_3dvariability.py")
            +" "+"bdb:"+path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_ANO#isac_substack")
            +" "+ "--output_dir="+old_output_folder
            +" --var3D="+filename3DvarOld
            + " --ave3D="+filename3DavgOld
            + " --sym=c5"
            +" --CTF"
            + " --var2D=''"
            + " --ave2D=''"
            + " --window=290"
        )

        a = subprocess.run(args=[testargs_new], shell=True, stderr=subprocess.STDOUT)
        b = subprocess.run(args=[testargs_old], shell=True,  stderr=subprocess.STDOUT)

        print_new = a.stdout.decode('utf8').split('\n')
        print_oldv = b.stdout.decode('utf8').split('\n')


        self.assertEqual(print_new[16].split("ERROR")[1],
                         ' => Input stack is not prepared for symmetry, please follow instructions')
        remove_dir(new_output_folder)
        remove_dir(old_output_folder)