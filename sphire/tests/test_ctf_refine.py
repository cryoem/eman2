from __future__ import print_function
from __future__ import division

from numpy import array_equal

from sphire.bin_py3 import sp_ctf_refine as oldfu
from sphire.bin import sp_ctf_refine as fu

from os import path
from sphire.tests.test_module import ABSOLUTE_OLDBIN_PATH,ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,ABSOLUTE_BIN_PATH,remove_dir,IMAGE_2D,get_arg_from_pickle_file,ABSOLUTE_PATH, give_ali_vol_data, give_ormq_data
import unittest
from sp_utilities import get_im
from EMAN2_cppwrap import EMAN2Ctf
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
import numpy
import sys


"""
WHAT IS MISSING:
0) compute_average I need to create images with 'xform.align2d' key to 
1) this f has been cleaned: "calc_similarity_complex"
2) 

RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) Test_helperFunctions.test_adjust_pw_to_model_2DimgBlank_return_zeroDivisionException. In the cleaned version of the libs the exception is not raised and the function returns an image 
2) Test_Error_cases::test_negative_pixelSize seems to create an image with random values .... but the compatibility test is ok. why?
3) The Test_run case failed !!!!!!! 

In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:
1) 
"""


class Test_helperFunctions(unittest.TestCase):
    def test_create_ctf_list(self):
        expected_defocus = [0.8999999761581421, 0.9200000166893005, 0.9399999976158142, 0.9599999785423279, 0.9800000190734863, 1.0, 1.0199999809265137, 1.0399999618530273, 1.059999942779541, 1.0800000429153442, 1.100000023841858, 1.1200000047683716]
        ctf = EMAN2Ctf()
        ctf.from_dict(
            {
                "defocus": 1,
                "cs": 2,
                "voltage": 300,
                "apix": 1.5,
                "bfactor": 0,
                "ampcont": 0.1,
                "dfdiff": 0.1,
                "dfang": 0.1,
            }
        )

        # dict([(key, value) for key, value in ctf.to_dict().items() if key not in ('background', 'snr')])
        return_new =fu.create_ctf_list(current_ctf=ctf, def_search_range=0.1, def_step_size=0.02)
        return_old = oldfu.create_ctf_list(current_ctf=ctf, def_search_range=0.1, def_step_size=0.02)
        for i,j in zip(return_new,return_old):
            self.assertTrue(i.defocus, j.defocus)
        for i,j in zip(return_new,expected_defocus):
            self.assertTrue(i.defocus, j)

    def test_calc_2d_projection(self):
        (volft, params, interpolation_method, return_real) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/projection.prgl"))[0]
        volft, refv = give_ali_vol_data()
        volft = volft.do_fft()
        return_new=fu.calc_2d_projection(particle_volume=volft, projection_angle_and_shift=params, interpolation=interpolation_method, return_real=return_real)
        return_old = oldfu.calc_2d_projection(particle_volume=volft, projection_angle_and_shift=params,interpolation=interpolation_method, return_real=return_real)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        # self.assertTrue(array_equal(return_new.get_2dview().flatten().tolist()[:100],[101.89811706542969, 91.55075073242188, 72.54855346679688, 60.872314453125, 78.21893310546875, 91.90130615234375, 77.81035614013672, 82.16301727294922, 114.83625793457031, 137.77029418945312, 175.0357208251953, 251.92584228515625, 264.3248291015625, 178.52642822265625, 100.36066436767578, 111.91815948486328, 139.7021484375, 142.91371154785156, 134.08544921875, 116.979736328125, 74.44189453125, 40.43936538696289, 24.7244930267334, 34.545494079589844, 41.635902404785156, 41.44943618774414, 32.24332046508789, 25.620906829833984, 18.462657928466797, 12.844551086425781, 11.637168884277344, 5.591159820556641, 7.756763458251953, 18.637340545654297, 26.93962860107422, 23.060644149780273, 18.233295440673828, 3.2267913818359375, -12.685356140136719, 0.15061187744140625, 34.825313568115234, 43.70667266845703, 28.631973266601562, 33.28607177734375, 56.45341491699219, 70.11874389648438, 91.97492218017578, 125.14620208740234, 131.40933227539062, 123.06419372558594, 154.791015625, 200.84963989257812, 190.6787109375, 124.27879333496094, 83.90874481201172, 75.40718841552734, 88.04712677001953, 103.45071411132812, 119.39303588867188, 105.75433349609375, 103.92453002929688, 104.42353820800781, 96.89095306396484, 106.53619384765625, 100.96305847167969, 87.70775604248047, 72.27996826171875, 62.253379821777344, 79.1793212890625, 81.55354309082031, 71.97140502929688, 106.90995025634766, 154.16908264160156, 163.37646484375, 166.71383666992188, 202.56167602539062, 203.20172119140625, 137.35975646972656, 69.41737365722656, 74.89479064941406, 112.64189147949219, 126.77521514892578, 100.93851470947266, 81.2967529296875, 64.08192443847656, 53.83599853515625, 54.185874938964844, 71.47483825683594, 63.05082702636719, 40.84798812866211, 30.750411987304688, 36.67400360107422, 34.88325500488281, 33.18653869628906, 32.885746002197266, 30.127445220947266, 33.49412536621094, 45.50977325439453, 60.73741149902344, 72.96221160888672]))

    def test_calc_similarity_real(self):
        (img1, img2, mask) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc"))[0]
        img1, img2 = give_ali_vol_data()
        # volft = volft.do_fft()
        return_new=fu.calc_similarity_real(particle_image=img1, projection_ctf_applied_2d=img2, mask=mask)
        return_old = oldfu.calc_similarity_real(particle_image=img1, projection_ctf_applied_2d=img2, mask=mask)
        self.assertEqual(return_new, return_old)
        # self.assertEqual(return_new, 0.8129369020462036)

    def test_apply_ctf(self):
        ctf = EMAN2Ctf()
        ctf.from_dict(
            {
                "defocus": 1,
                "cs": 2,
                "voltage": 300,
                "apix": 1.5,
                "bfactor": 0,
                "ampcont": 0.1,
                "dfdiff": 0.1,
                "dfang": 0.1,
            }
        )
        volft, refv = give_ormq_data()
        return_new = fu.apply_ctf(projection_2d=volft, ctf=ctf)
        return_old = oldfu.apply_ctf(projection_2d=volft, ctf=ctf)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        # self.assertTrue(array_equal(return_new.get_2dview().flatten(),  [0.024286171421408653, 0.058696337044239044, 0.10402560234069824, 0.12992535531520844, 0.09275151789188385, 0.03502361848950386, 0.11739957332611084, 0.06164858862757683, 0.07032382488250732, 0.03957751765847206, 0.021595070138573647, 0.06361360102891922, 0.08738928288221359, 0.1192685067653656, 0.10259803384542465, 0.040069177746772766, 0.10591607540845871, 0.09205082803964615, 0.12429622560739517, 0.07694292068481445, 0.09733756631612778, 0.1372169703245163, 0.14484426379203796, 0.14918211102485657, 0.11907980591058731, 0.10208652168512344, 0.1392429918050766, 0.13015149533748627, 0.11991294473409653, 0.08955641835927963, 0.06140293926000595, 0.10304044187068939, 0.12554489076137543, 0.11262138932943344, 0.06353975832462311, 0.08587218821048737, 0.11244022101163864, 0.09934824705123901, 0.09497562795877457, 0.061048153787851334, 0.07756707817316055, 0.10022807121276855, 0.12350466102361679, 0.1090707778930664, 0.08835660666227341, 0.13008004426956177, 0.1098373606801033, 0.09089691191911697, 0.08573641628026962, 0.05659455060958862, 0.004100685473531485, 0.008962074294686317, 0.032936085015535355, 0.030323926359415054, -0.023060372099280357, 0.013550573028624058, 0.0064983367919921875, -0.01710440404713154, -0.01073240116238594, -0.017573906108736992, -0.03412095457315445, -0.05568891391158104, -0.052081622183322906, -0.04967673122882843, -0.05739298835396767, -0.031187398359179497, -0.040683239698410034, -0.0780254602432251, -0.07197131961584091, -0.048417337238788605, -0.06368181854486465, -0.08691969513893127, -0.07733665406703949, -0.06406906992197037, -0.0645952969789505, -0.06410295516252518, -0.07824750989675522, -0.09381715208292007, -0.10522931814193726, -0.07656799256801605, -0.05408206209540367, -0.10620465874671936, -0.06528617441654205, -0.06386560201644897, -0.02619960717856884, -0.06233495473861694, -0.06505712121725082, -0.051096897572278976, -0.05767035484313965, -0.023533765226602554, -0.033616144210100174, -0.060482509434223175, -0.02472449466586113, -0.06295625865459442, -0.01411000732332468, -0.035392340272665024, -0.06724553555250168, -0.07538473606109619, -0.06100975722074509, -0.02385104075074196]))

    def test_create_mask(self):
        return_new=fu.create_mask(size=IMAGE_2D.get_2dview().shape, index_particle_pixels=numpy.where(IMAGE_2D.get_2dview() > 0.1), fraction_size=0.25)
        return_old = oldfu.create_mask(size=IMAGE_2D.get_2dview().shape,
                                    index_particle_pixels=numpy.where(IMAGE_2D.get_2dview() > 0.1), fraction_size=0.25)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        self.assertTrue(array_equal( return_new.get_2dview().flatten(), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_find_index_best_projection(self):
        (img1, img2, mask) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc"))[0]
        img1, img2 = give_ali_vol_data()
        return_new = fu.find_index_best_projection(particle_image=img1, reprojections=[img2,img1], mask=mask)
        return_old = oldfu.find_index_best_projection(particle_image=img1, reprojections=[img2, img1], mask=mask)
        self.assertEqual(return_old,return_new)
        self.assertEqual(1, return_new)

    @unittest.skip("do not work")
    def test_optimize(self):
        (img1, img2, mask) = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc"))[0]
        ctf = EMAN2Ctf()
        ctf.from_dict(
            {
                "defocus": 1,
                "cs": 2,
                "voltage": 300,
                "apix": 1.5,
                "bfactor": 0,
                "ampcont": 0.1,
                "dfdiff": 0.1,
                "dfang": 0.1,
            }
        )
        ctf_list =fu.create_ctf_list(current_ctf=ctf, def_search_range=0.1, def_step_size=0.02)
        return_new=fu.optimize(particle_image=img1, projection_2d=[img1,img1], ctf_list=ctf_list, masks=[mask,mask])


    def test_refine_defocus_with_error_est(self):
        pass




class Test_run(unittest.TestCase):
    old_output_folder="ctfRefineManualOld"
    new_output_folder = "ctfRefineManualNew"
    bdb_name= "ctf_refined"

    @classmethod
    def tearDownClass(cls):
        remove_dir(cls.new_output_folder)
        remove_dir(cls.old_output_folder)

    def test_(self):
        testargs_new = [path.join(ABSOLUTE_BIN_PATH, "sp_ctf_refine.py"),
                        "manual",
                        "bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_NEW#isac_substack"),
                        self.new_output_folder,
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"12_POSTREFINER","vol_combined.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","final_params_028.txt")]

        print("2nd round")

        testargs_old = [path.join(ABSOLUTE_OLDBIN_PATH, "sp_ctf_refine.py"),
                        "manual",
                        "bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"06_SUBSTACK_NEW#isac_substack"),
                        self.old_output_folder,
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW,"12_POSTREFINER","vol_combined.hdf"),
                        path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW, "11_MERIDIEN","final_params_028.txt")]
        with patch.object(sys, 'argv', testargs_new):
            fu._main_()
        with patch.object(sys, 'argv', testargs_old):
            oldfu._main_()
        fpath_new = 'bdb:{0}#{1}'.format(self.new_output_folder, self.bdb_name)
        fpath_old = 'bdb:{0}#{1}'.format(self.old_output_folder, self.bdb_name)
        return_old = get_im(fpath_old)
        return_new = get_im(fpath_new)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        # self.assertTrue(array_equal(return_new.get_3dview().flatten().tolist()[0:100],[-0.38441434502601624, -0.4573197066783905, -0.39960816502571106, -0.10431570559740067, 0.018583765253424644, 0.37694355845451355, 0.4353419840335846, 0.06498070061206818, -0.38467326760292053, -0.7100721001625061, -0.7571740746498108, -0.9328497648239136, -0.8550546169281006, -0.21433517336845398, 0.34148913621902466, 0.30916473269462585, 0.39394184947013855, 0.4447155296802521, -0.032833874225616455, 0.4121330976486206, 0.638403594493866, 0.34945985674858093, 0.2751978039741516, -0.14872148633003235, -0.5307353138923645, -0.6574574112892151, -0.8945914506912231, -1.107301950454712, -1.3676011562347412, -1.2225056886672974, -1.2273787260055542, -0.7121877074241638, 0.06733409315347672, 0.29879996180534363, 0.6046097874641418, 0.496036559343338, 0.7214235663414001, 0.7652679681777954, 0.9319247603416443, 0.11259084939956665, -0.3054805397987366, -0.8183746933937073, -1.4462037086486816, -1.4775571823120117, -1.2116808891296387, -0.911469042301178, -0.9934141039848328, -0.9891195297241211, -1.0073580741882324, -1.234678864479065, -1.4292954206466675, -1.5580313205718994, -1.446936845779419, -1.67134428024292, -2.0199873447418213, -2.0543317794799805, -2.2354137897491455, -1.99867844581604, -1.4405670166015625, -1.7088592052459717, -2.060110569000244, -1.6585056781768799, -1.0386717319488525, -0.4783007502555847, -0.12102964520454407, 0.16874085366725922, -0.12371158599853516, 0.01841016113758087, -0.6442297101020813, -0.6972493529319763, -0.1051267459988594, 0.016293810680508614, -0.10689342767000198, -0.472159206867218, -0.4886413514614105, -0.09828135371208191, -0.7409976720809937, -0.8893253207206726, -1.1477444171905518, -0.35942375659942627, 0.030692437663674355, 0.5380961298942566, 1.06252121925354, 0.6961002349853516, 0.7563707232475281, 1.0197871923446655, 0.9021824598312378, 0.34891727566719055, 0.402149498462677, 0.634938657283783, 0.6044892072677612, 0.24482420086860657, -0.4132026731967926, -0.9945327043533325, -1.260424256324768, -1.6913681030273438, -1.7999876737594604, -1.5891680717468262, -1.3196191787719727, -0.9893324375152588]))



