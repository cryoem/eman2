from __future__ import print_function
from __future__ import division
from past.utils import old_div
import unittest

from numpy import array_equal

from sphire.libpy_py3 import sp_isac as oldfu
from sphire.libpy import sp_isac as fu

from .test_module import IMAGE_2D
"""
These functions have been cleaned
class Test_iter_isac(unittest.TestCase):
    def test_iter_isac(self):
        oldv = oldfu.iter_isac(
            stack="",
            ir="",
            ou="",
            rs="",
            xr="",
            yr="",
            ts="",
            maxit="",
            CTF="",
            snr="",
            dst="",
            FL="",
            FH="",
            FF="",
            init_iter="",
            main_iter="",
            iter_reali="",
            match_first="",
            max_round="",
            match_second="",
            stab_ali="",
            thld_err="",
            indep_run="",
            thld_grp="",
            img_per_grp="",
            generation="",
            candidatesexist=False,
            random_seed=None,
            new=False,
        )
        v = fu.iter_isac(
            stack="",
            ir="",
            ou="",
            rs="",
            xr="",
            yr="",
            ts="",
            maxit="",
            CTF="",
            snr="",
            dst="",
            FL="",
            FH="",
            FF="",
            init_iter="",
            main_iter="",
            iter_reali="",
            match_first="",
            max_round="",
            match_second="",
            stab_ali="",
            thld_err="",
            indep_run="",
            thld_grp="",
            img_per_grp="",
            generation="",
            candidatesexist=False,
            random_seed=None,
            new=False,
        )
        pass


class Test_isac_MPI(unittest.TestCase):
    def test_isac_MPI(self):
        oldv = oldfu.isac_MPI(
            stack="",
            refim="",
            maskfile=None,
            outname="avim",
            ir=1,
            ou=-1,
            rs=1,
            xrng=0,
            yrng=0,
            step=1,
            maxit=30,
            isac_iter=10,
            CTF=False,
            snr=1.0,
            rand_seed=-1,
            color=0,
            comm=-1,
            stability=False,
            stab_ali=5,
            iter_reali=1,
            thld_err=1.732,
            FL=0.1,
            FH=0.3,
            FF=0.2,
            dst=90.0,
            method="",
        )
        v = fu.isac_MPI(
            stack="",
            refim="",
            maskfile=None,
            outname="avim",
            ir=1,
            ou=-1,
            rs=1,
            xrng=0,
            yrng=0,
            step=1,
            maxit=30,
            isac_iter=10,
            CTF=False,
            snr=1.0,
            rand_seed=-1,
            color=0,
            comm=-1,
            stability=False,
            stab_ali=5,
            iter_reali=1,
            thld_err=1.732,
            FL=0.1,
            FH=0.3,
            FF=0.2,
            dst=90.0,
            method="",
        )
        pass


class Test_match_independent_runs(unittest.TestCase):
    def test_match_independent_runs(self):
        oldv = oldfu.match_independent_runs(data="", refi="", n_group="", T="")
        v = fu.match_independent_runs(data="", refi="", n_group="", T="")
        pass


class Test_match_2_way(unittest.TestCase):
    def test_match_2_way(self):
        oldv = oldfu.match_2_way(
            data="",
            refi="",
            indep_run="",
            thld_grp="",
            FH="",
            FF="",
            find_unique=True,
            wayness=2,
            suffix="",
        )
        v = fu.match_2_way(
            data="",
            refi="",
            indep_run="",
            thld_grp="",
            FH="",
            FF="",
            find_unique=True,
            wayness=2,
            suffix="",
        )
        pass


class Test_get_unique_averages(unittest.TestCase):
    def test_generate_random_averages(self):
        oldv = oldfu.get_unique_averages(data="", indep_run="", m_th=0.45)
        v = fu.get_unique_averages(data="", indep_run="", m_th=0.45)
        pass
"""


class Test_generate_random_averages(unittest.TestCase):
    def test_null_K(self):
        return_old = oldfu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=0, rand_seed=-1)
        return_new = fu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=0, rand_seed=-1)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal([], return_new))

    def test(self):
        return_old = oldfu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=2, rand_seed=-1)
        return_new = fu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=2, rand_seed=-1)
        for i,j in zip(return_new,return_old):
            self.assertTrue(array_equal(i.get_2dview(),j.get_2dview()))
        self.assertTrue(array_equal(return_new[0].get_2dview().flatten(),[0.009504491463303566, 0.025885052978992462, 0.015371464192867279, 0.029651766642928123, 0.025623319670557976, 0.023996084928512573, 0.023316310718655586, 0.03626585379242897, 0.042238689959049225, 0.053261399269104004, 0.06996519863605499, 0.05416787788271904, 0.050994712859392166, 0.03554266691207886, 0.03604980185627937, 0.07005909085273743, 0.056754179298877716, 0.06729267537593842, 0.0899617150425911, 0.08004479855298996, 0.07206107676029205, 0.07158395648002625, 0.08500781655311584, 0.08074058592319489, 0.08976095914840698, 0.09553121030330658, 0.09733162075281143, 0.12153391540050507, 0.09777011722326279, 0.0612066276371479, 0.060473889112472534, 0.0832795649766922, 0.07990699261426926, 0.0726018100976944, 0.10390139371156693, 0.12692593038082123, 0.08997570723295212, 0.05740871652960777, 0.05622498691082001, 0.05523042380809784, 0.013907668180763721, 0.0071470243856310844, 0.01511574536561966, 2.5205374186043628e-05, 0.008231919258832932, -0.020773129537701607, -0.034199729561805725, -0.04089483618736267, -0.042460259050130844, -0.06925757229328156, -0.06893884390592575, -0.08000176399946213, -0.11662115156650543, -0.111984983086586, -0.11971071362495422, -0.1273496150970459, -0.12249226123094559, -0.1453358680009842, -0.14758040010929108, -0.15034900605678558, -0.17081016302108765, -0.2014905959367752, -0.2121349573135376, -0.22736789286136627, -0.24315771460533142, -0.2552821934223175, -0.23703180253505707, -0.2393375188112259, -0.2672199606895447, -0.28808265924453735, -0.3236375153064728, -0.3262620270252228, -0.35172849893569946, -0.3602631986141205, -0.35741564631462097, -0.3575122356414795, -0.38925597071647644, -0.377326101064682, -0.38598355650901794, -0.39209896326065063, -0.3882087767124176, -0.3639817535877228, -0.3711523711681366, -0.37047016620635986, -0.39362388849258423, -0.40711337327957153, -0.3925972580909729, -0.4149233400821686, -0.41900205612182617, -0.4641905426979065, -0.46107935905456543, -0.46086275577545166, -0.4773290157318115, -0.473482221364975, -0.4543262720108032, -0.44096702337265015, -0.4387476146221161, -0.4229215085506439, -0.4376510977745056, -0.4369300603866577]))
        self.assertTrue(array_equal(return_new[1].get_2dview().flatten(),[0.009504491463303566, 0.025885052978992462, 0.015371464192867279, 0.029651766642928123, 0.025623319670557976, 0.023996084928512573, 0.023316310718655586, 0.03626585379242897, 0.042238689959049225, 0.053261399269104004, 0.06996519863605499, 0.05416787788271904, 0.050994712859392166, 0.03554266691207886, 0.03604980185627937, 0.07005909085273743, 0.056754179298877716, 0.06729267537593842, 0.0899617150425911, 0.08004479855298996, 0.07206107676029205, 0.07158395648002625, 0.08500781655311584, 0.08074058592319489, 0.08976095914840698, 0.09553121030330658, 0.09733162075281143, 0.12153391540050507, 0.09777011722326279, 0.0612066276371479, 0.060473889112472534, 0.0832795649766922, 0.07990699261426926, 0.0726018100976944, 0.10390139371156693, 0.12692593038082123, 0.08997570723295212, 0.05740871652960777, 0.05622498691082001, 0.05523042380809784, 0.013907668180763721, 0.0071470243856310844, 0.01511574536561966, 2.5205374186043628e-05, 0.008231919258832932, -0.020773129537701607, -0.034199729561805725, -0.04089483618736267, -0.042460259050130844, -0.06925757229328156, -0.06893884390592575, -0.08000176399946213, -0.11662115156650543, -0.111984983086586, -0.11971071362495422, -0.1273496150970459, -0.12249226123094559, -0.1453358680009842, -0.14758040010929108, -0.15034900605678558, -0.17081016302108765, -0.2014905959367752, -0.2121349573135376, -0.22736789286136627, -0.24315771460533142, -0.2552821934223175, -0.23703180253505707, -0.2393375188112259, -0.2672199606895447, -0.28808265924453735, -0.3236375153064728, -0.3262620270252228, -0.35172849893569946, -0.3602631986141205, -0.35741564631462097, -0.3575122356414795, -0.38925597071647644, -0.377326101064682, -0.38598355650901794, -0.39209896326065063, -0.3882087767124176, -0.3639817535877228, -0.3711523711681366, -0.37047016620635986, -0.39362388849258423, -0.40711337327957153, -0.3925972580909729, -0.4149233400821686, -0.41900205612182617, -0.4641905426979065, -0.46107935905456543, -0.46086275577545166, -0.4773290157318115, -0.473482221364975, -0.4543262720108032, -0.44096702337265015, -0.4387476146221161, -0.4229215085506439, -0.4376510977745056, -0.4369300603866577]))



    def test_K_higher_lenData_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=5, rand_seed=-1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.generate_random_averages(data=[IMAGE_2D,IMAGE_2D,IMAGE_2D], K=5, rand_seed=-1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_K_higher_lenData_returns_IndexError_list_index_out_of_range2(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.generate_random_averages(data=[], K=5, rand_seed=-1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.generate_random_averages(data=[], K=5, rand_seed=-1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))