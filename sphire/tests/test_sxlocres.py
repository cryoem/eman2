import copy
from sphire.libpy import sp_utilities
import numpy
from sphire.bin import sp_locres

import unittest
class TestMakeAngRes(unittest.TestCase):
    def test_input_image_returns_values_freq_to_real(self):
        nx = 10
        ny = 10
        nz = 10
        apix = 2
        test_data = numpy.arange(nx * ny * nz, dtype=numpy.float32).reshape(
            nx, ny, nz
        ) / (2 * nx * ny * nz)
        mask = test_data > 0

        blank_data = sp_utilities.model_blank(nx, ny, nz)
        data_in = blank_data.get_3dview()
        data_in[...] = copy.deepcopy(test_data)
        returned_data = sp_locres.makeAngRes(blank_data, nx, ny, nz, apix)

        expected_data = copy.deepcopy(test_data)
        expected_data[mask] = apix / test_data[mask]
        assert numpy.array_equal(expected_data, returned_data.get_3dview())

    @unittest.skip("Both values are returned differently")
    def test_input_image_returns_values_real_to_freq(self):
        nx = 10
        ny = 10
        nz = 10
        apix = 2
        test_data = numpy.arange(nx * ny * nz, dtype=numpy.float32).reshape(nx, ny, nz)
        mask = test_data >= 2 * apix

        blank_data = sp_utilities.model_blank(nx, ny, nz)
        data_in = blank_data.get_3dview()
        data_in[...] = copy.deepcopy(test_data)
        returned_data = sp_locres.makeAngRes(blank_data, nx, ny, nz, apix, False)

        expected_data = copy.deepcopy(test_data)
        expected_data[mask] = apix / expected_data[mask]
        expected_data[~mask] = 0
        print("EXP", expected_data)
        print("RET", returned_data.get_3dview())
        assert numpy.array_equal(expected_data, returned_data.get_3dview())
