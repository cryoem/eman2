from __future__ import print_function
from __future__ import division
import unittest
import sys

from sphire.bin import sp_rewindow as fu
from sphire.bin_py3 import sp_rewindow as oldfu

from .test_module import IMAGE_2D,IMAGE_3D
from copy import deepcopy
from numpy import array_equal
try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

import mpi
from mpi import *
mpi.mpi_init(0, [])

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO

"""
WHAT IS MISSING:
0) get_time_stamp_suffix --> returns a timestamp ... not testable

RESULT AND KNOWN ISSUES
1) 

In these tests there is a bug --> syntax error:
1)

In these tests there is a strange behavior:
1) 
"""

class Test_helperFunctions(unittest.TestCase):
    def test_is_float_True(self):
        self.assertTrue(fu.is_float(3))
        self.assertTrue(oldfu.is_float(3))

    def test_is_float_False(self):
        self.assertFalse(fu.is_float("d"))
        self.assertFalse(oldfu.is_float("d"))

    def test_get_cmd_line(self):
        cmdLine=["this", "is", "a", "test"]
        with patch.object(sys, 'argv', cmdLine):
            return_new = fu.get_cmd_line()
            return_old = oldfu.get_cmd_line()
            self.assertEqual(return_new,return_old)
            self.assertEqual(return_new, 'Shell line command: this  is  a  test  ')

    def test_mrk_resample2d_subrate1(self):
        return_new = fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=1, target_size = None)
        return_old = oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=1, target_size = None)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        self.assertTrue(array_equal(return_new.get_2dview().flatten(),[0.009504491463303566, 0.025885052978992462, 0.015371464192867279, 0.029651766642928123,
                   0.025623319670557976, 0.023996084928512573, 0.023316310718655586, 0.03626585379242897,
                   0.042238689959049225, 0.053261399269104004, 0.06996519863605499, 0.05416787788271904,
                   0.050994712859392166, 0.03554266691207886, 0.03604980185627937, 0.07005909085273743,
                   0.056754179298877716, 0.06729267537593842, 0.0899617150425911, 0.08004479855298996,
                   0.07206107676029205, 0.07158395648002625, 0.08500781655311584, 0.08074058592319489,
                   0.08976095914840698, 0.09553121030330658, 0.09733162075281143, 0.12153391540050507,
                   0.09777011722326279, 0.0612066276371479, 0.060473889112472534, 0.0832795649766922,
                   0.07990699261426926, 0.0726018100976944, 0.10390139371156693, 0.12692593038082123,
                   0.08997570723295212, 0.05740871652960777, 0.05622498691082001, 0.05523042380809784,
                   0.013907668180763721, 0.0071470243856310844, 0.01511574536561966, 2.5205374186043628e-05,
                   0.008231919258832932, -0.020773129537701607, -0.034199729561805725, -0.04089483618736267,
                   -0.042460259050130844, -0.06925757229328156, -0.06893884390592575, -0.08000176399946213,
                   -0.11662115156650543, -0.111984983086586, -0.11971071362495422, -0.1273496150970459,
                   -0.12249226123094559, -0.1453358680009842, -0.14758040010929108, -0.15034900605678558,
                   -0.17081016302108765, -0.2014905959367752, -0.2121349573135376, -0.22736789286136627,
                   -0.24315771460533142, -0.2552821934223175, -0.23703180253505707, -0.2393375188112259,
                   -0.2672199606895447, -0.28808265924453735, -0.3236375153064728, -0.3262620270252228,
                   -0.35172849893569946, -0.3602631986141205, -0.35741564631462097, -0.3575122356414795,
                   -0.38925597071647644, -0.377326101064682, -0.38598355650901794, -0.39209896326065063,
                   -0.3882087767124176, -0.3639817535877228, -0.3711523711681366, -0.37047016620635986,
                   -0.39362388849258423, -0.40711337327957153, -0.3925972580909729, -0.4149233400821686,
                   -0.41900205612182617, -0.4641905426979065, -0.46107935905456543, -0.46086275577545166,
                   -0.4773290157318115, -0.473482221364975, -0.4543262720108032, -0.44096702337265015,
                   -0.4387476146221161, -0.4229215085506439, -0.4376510977745056, -0.4369300603866577]))

    def test_mrk_resample2d_negative_targetsize_assertionError(self):
        with self.assertRaises(AssertionError) as cm_new:
            fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = -10)
        with self.assertRaises(AssertionError) as cm_old:
            oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = -10)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception), "")

    def test_mrk_resample2d_null_targetsize(self):
        with self.assertRaises(AssertionError) as cm_new:
            fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = 0)
        with self.assertRaises(AssertionError) as cm_old:
            oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = 0)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception), "")

    def test_mrk_resample2d_positive_targetsize(self):
        return_new = fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = 5)
        return_old = oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = 5)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        a = return_new.get_2dview().flatten().tolist()
        self.assertTrue(array_equal(return_new.get_2dview().flatten(),[0.0821991041302681, 0.0839809998869896, 0.1071862056851387, 0.10582967847585678, 0.11003056913614273, 0.06691155582666397, 0.06422650068998337, 0.09262319654226303, 0.07385403662919998, 0.03432386741042137, -0.04923679307103157, -0.04758632183074951, -0.061087723821401596, -0.07956168800592422, -0.08740868419408798, -0.16807080805301666, -0.1731511950492859, -0.19433154165744781, -0.180768221616745, -0.182774156332016, -0.31768250465393066, -0.3294888138771057, -0.3275138735771179, -0.34638598561286926, -0.3418603241443634]))

    def test_mrk_resample2d_subrate_less1_returns(self):
        return_new = fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = None)
        return_old = oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0.9, target_size = None)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        self.assertTrue(array_equal(return_new.get_2dview().flatten(),[-0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.047721877694129944, -0.00979264173656702, -0.004350133705884218, -0.011892050504684448, -0.0016624846030026674, 0.0011983346194028854, -0.0002359417558182031, 0.0050877770408988, 0.01729704812169075, 0.03196311742067337, -0.047721877694129944, 0.07785306125879288, 0.06396913528442383, 0.06371123343706131, 0.04520949348807335, 0.0620175376534462, 0.0731816217303276, 0.0753120705485344, 0.09960953146219254, 0.08418302237987518, -0.047721877694129944, 0.06732919067144394, 0.07741532474756241, 0.0821991041302681, 0.0839809998869896, 0.1071862056851387, 0.10582967847585678, 0.11003056913614273, 0.09924037009477615, 0.06709948182106018, -0.047721877694129944, 0.050387803465127945, 0.06995890289545059, 0.06691155582666397, 0.06422650068998337, 0.09262319654226303, 0.07385403662919998, 0.03432386741042137, 0.028595253825187683, 0.020939644426107407, -0.047721877694129944, -0.027883881703019142, -0.03170311823487282, -0.04923679307103157, -0.04758632183074951, -0.061087723821401596, -0.07956168800592422, -0.08740868419408798, -0.08973980695009232, -0.10579903423786163, -0.047721877694129944, -0.12305420637130737, -0.1417393535375595, -0.16807080805301666, -0.1731511950492859, -0.19433154165744781, -0.180768221616745, -0.182774156332016, -0.2067779004573822, -0.22201569378376007, -0.047721877694129944, -0.28302663564682007, -0.2964749038219452, -0.31768250465393066, -0.3294888138771057, -0.3275138735771179, -0.34638598561286926, -0.3418603241443634, -0.34814080595970154, -0.367236465215683, -0.047721877694129944, -0.3802380859851837, -0.35358381271362305, -0.368375688791275, -0.37149152159690857, -0.3932974934577942, -0.3932115137577057, -0.40853428840637207, -0.41270506381988525, -0.4473392069339752, -0.047721877694129944, -0.4593888521194458, -0.4584948420524597, -0.4698401093482971, -0.4616740942001343, -0.4503171443939209, -0.4427659511566162, -0.4284581243991852, -0.43833982944488525, -0.45278245210647583]))

    def test_mrk_resample2d_subrate_null_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0, target_size = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=0, target_size = None)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mrk_resample2d_subrate_negative_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=-0.9, target_size = None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=-0.9, target_size = None)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_mrk_resample2d_subrate_higher1(self):
        return_new = fu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=1.9, target_size=None)
        return_old = oldfu.mrk_resample2d(img=deepcopy(IMAGE_2D), sub_rate=1.9, target_size=None)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))
        # self.assertTrue(array_equal(return_new.get_2dview().flatten(),[0.08320573717355728, 0.07271669059991837, 0.06818031519651413, 0.09391328692436218, 0.12102993577718735, 0.11594352126121521, 0.09705159068107605, 0.09719530493021011, 0.10561344772577286, 0.0967702567577362, 0.08243812620639801, 0.07083393633365631, 0.07310925424098969, 0.10512372851371765, 0.13259288668632507, 0.12432464212179184, 0.09990441799163818, 0.08922552317380905, 0.08367384225130081, 0.06763394176959991, 0.07941532135009766, 0.06784776598215103, 0.06707008928060532, 0.08599218726158142, 0.09862350672483444, 0.08668579906225204, 0.06695523113012314, 0.05619390681385994, 0.044884003698825836, 0.028007399290800095, 0.017739417031407356, 0.01256449893116951, 0.010886979289352894, 0.0130070261657238, 0.007318354211747646, -0.006462703458964825, -0.016188614070415497, -0.02113327383995056, -0.03109278529882431, -0.04217332974076271, -0.07314638793468475, -0.0707002580165863, -0.0658802017569542, -0.07014332711696625, -0.08426815271377563, -0.09429328888654709, -0.09492906928062439, -0.0977659597992897, -0.10696263611316681, -0.11053655296564102, -0.11953753978013992, -0.11697637289762497, -0.1076805517077446, -0.114369235932827, -0.12980273365974426, -0.13085418939590454, -0.12293270975351334, -0.12762132287025452, -0.13854701817035675, -0.13642024993896484, -0.151829794049263, -0.15250827372074127, -0.14399021863937378, -0.15832467377185822, -0.17718826234340668, -0.16823138296604156, -0.15090616047382355, -0.15815560519695282, -0.17099638283252716, -0.16258449852466583, -0.2475888729095459, -0.24331074953079224, -0.23039957880973816, -0.252596914768219, -0.27765825390815735, -0.2613894045352936, -0.23789778351783752, -0.251621276140213, -0.26763278245925903, -0.24817988276481628, -0.3585677444934845, -0.34278038144111633, -0.3169642388820648, -0.3409063518047333, -0.3708537220954895, -0.35105839371681213, -0.3286631107330322, -0.3563586175441742, -0.3792392909526825, -0.34514787793159485, -0.38008394837379456, -0.35952889919281006, -0.324637770652771, -0.3495219647884369, -0.38339370489120483, -0.3623982071876526, -0.343816339969635, -0.38571256399154663, -0.41599899530410767, -0.3707023561000824]))

    def test_mrk_resample2d__with_3Dimage_assertionError(self):
        with self.assertRaises(AssertionError) as cm_new:
            fu.mrk_resample2d(img=deepcopy(IMAGE_3D), sub_rate=1.9, target_size=None)
        with self.assertRaises(AssertionError) as cm_old:
            oldfu.mrk_resample2d(img=deepcopy(IMAGE_3D), sub_rate=1.9, target_size=None)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(str(cm_new.exception), "")