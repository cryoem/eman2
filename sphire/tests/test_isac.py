from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy
import unittest
from test_module import get_real_data,create_kb
from sphire.libpy import sp_isac as oldfu
from sphire.libpy_py3 import sp_isac as fu


class Test_iter_isac(unittest.TestCase):
    def test_iter_isac(self):
        oldv = oldfu.iter_isac(stack="", ir="", ou="", rs="", xr="", yr="", ts="", maxit="", CTF="", snr="", dst="", FL="", FH="", FF="", init_iter="", main_iter="", iter_reali="", match_first="", max_round="", match_second="", stab_ali="", thld_err="", indep_run="", thld_grp="", img_per_grp="", generation="", candidatesexist = False, random_seed=None, new = False)
        v = fu.iter_isac(stack="", ir="", ou="", rs="", xr="", yr="", ts="", maxit="", CTF="", snr="", dst="", FL="", FH="", FF="", init_iter="", main_iter="", iter_reali="", match_first="", max_round="", match_second="", stab_ali="", thld_err="", indep_run="", thld_grp="", img_per_grp="", generation="", candidatesexist = False, random_seed=None, new = False)
        pass

class Test_isac_MPI(unittest.TestCase):
    def test_isac_MPI(self):
        oldv = oldfu.isac_MPI(stack="", refim="", maskfile = None, outname = "avim", ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, maxit=30, isac_iter=10, CTF=False, snr=1.0, rand_seed=-1, color=0, comm=-1, stability=False, stab_ali=5, iter_reali=1, thld_err=1.732, FL=0.1, FH=0.3, FF=0.2, dst=90.0, method = "")
        v = fu.isac_MPI(stack="", refim="", maskfile = None, outname = "avim", ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, maxit=30, isac_iter=10, CTF=False, snr=1.0, rand_seed=-1, color=0, comm=-1, stability=False, stab_ali=5, iter_reali=1, thld_err=1.732, FL=0.1, FH=0.3, FF=0.2, dst=90.0, method = "")
        pass

class Test_match_independent_runs(unittest.TestCase):
    def test_match_independent_runs(self):
        oldv = oldfu.match_independent_runs(data="", refi="", n_group="", T="")
        v = fu.match_independent_runs(data="", refi="", n_group="", T="")
        pass


class Test_match_2_way(unittest.TestCase):
    def test_match_2_way(self):
        oldv = oldfu.match_2_way(data="", refi="", indep_run="", thld_grp="", FH="", FF="", find_unique=True, wayness=2,suffix="")
        v = fu.match_2_way(data="", refi="", indep_run="", thld_grp="", FH="", FF="", find_unique=True, wayness=2,suffix="")
        pass


class Test_generate_random_averages(unittest.TestCase):
    def test_generate_random_averages(self):
        oldv = oldfu.generate_random_averages(data="", K="", rand_seed = -1)
        v = fu.generate_random_averages(data="", K="", rand_seed = -1)
        pass


class Test_get_unique_averages(unittest.TestCase):
    def test_generate_random_averages(self):
        oldv = oldfu.get_unique_averages(data="", indep_run="", m_th=0.45)
        v = fu.get_unique_averages(data="", indep_run="", m_th=0.45)
        pass

