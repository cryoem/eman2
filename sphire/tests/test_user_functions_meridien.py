import unittest

from sphire.libpy import sp_user_functions_meridien as oldfu
#from sphire.libpy_py3 import sp_user_functions_meridien as fu



class Test_do_volume_mask(unittest.TestCase):
    def test_do_volume_mask(self):
        vol = oldfu.do_volume_mask(ref_data="")
        pass

class Test_compute_search_params(unittest.TestCase):
    def test_compute_search_params(self):
        new_range, step = oldfu.compute_search_params(acc_trans="", shifter="", old_range="")
        pass

class Test_ai_spa(unittest.TestCase):
    def test_ai_spa(self):
        keepgoing = oldfu.ai_spa( Tracker="", fff="", anger="", shifter="", do_local="", chout = False)
        pass

class Test_ai_filament(unittest.TestCase):
    def test_ai_filament(self):
        keepgoing = oldfu.ai_filament( Tracker="", fff="", anger="", shifter="", do_local="", chout = False)
        pass
