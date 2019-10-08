import unittest

from sphire.libpy import sp_ctf_refine_io as oldfu
#from sphire.libpy_py3 import sp_ctf_refine_io as fu

class Test_read_meridien_data(unittest.TestCase):
    def test_read_meridien_data(self):
        d_res = oldfu.read_meridien_data(meridien_path="") #  "final_params","first_halfmap", "second_halfmap","chunk0", "chunk1"
        pass

class Test_write_virtual_bdb_stack(unittest.TestCase):
    def test_write_virtual_bdb_stack(self):
        oldfu.write_virtual_bdb_stack(output_stack_path="", origin_stack_path="", refined_ctfs="", number_of_particles=None)
        pass

class Test_read_meridien_params(unittest.TestCase):
    def test_read_meridien_params(self):
        values = oldfu.read_meridien_params(path="")
        pass

class Test_read_particle(unittest.TestCase):
    def test_read_particle(self):
        particle_image = oldfu.read_particle(path="", index="", header_only=False)
        pass

class Test_prepare_volume(unittest.TestCase):
    def test_prepare_volume(self):
        vol = oldfu.prepare_volume(volume_path="", mask=None, resolution=None, pixel_size=None)
        pass

class Test_read_volume(unittest.TestCase):
    def test_read_volume(self):
        vol1, vol2, mask_vol = oldfu.read_volume(path_vol_1="", path_vol_2="", path_mask=None, resolution=None, pixel_size=None)
        pass

class Test_write_statistics(unittest.TestCase):
    def test_write_statistics(self):
        oldfu.write_statistics(output_stats_path="", mic_stats="")
        pass
