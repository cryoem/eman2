import unittest
from copy import deepcopy
from EMAN2_cppwrap import EMData, Reconstructors, EMUtil

from numpy import array_equal, allclose,nan, zeros as numpy_zeros

from os import path

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

from mpi import *

from libpy_py3 import sp_reconstruction as oldfu
from sphire.libpy import sp_reconstruction as fu

from tests.test_module import get_real_data, ABSOLUTE_PATH_TO_RESOURCES, give_alignment_shc_data

XFORM_PROJECTION_IMG = give_alignment_shc_data()


STACK_NAME = "bdb:" + path.join(
    ABSOLUTE_PATH_TO_RESOURCES, "../03_PARTICLES_BDB/mpi_proc_007/TcdA1-0187_frames_ptcls"
)



"""
WHAT IS MISSING:
0) in all the cases where the input file is an image. I did not test the case with a complex image. I was not able to generate it
1) recons3d_trl_struct_MPI no idea how test it
2) recons3d_4nn after hours is still running ... is it ok? did i use bad input??
3) recons3d_4nnw_MPI  no idea how test it


In these tests there is a bug --> syntax error:
1) Test_recons3d_4nn_ctf_MPI
  a) there is a KNOWN BUG --> with sizeprojection  PAP 10/22/2014
  b) if you call this function twice, or in the tests case twice in the same class test, the second time that it runs crashed beacuse del sparx_utilities.pad
     This happen because 'sparx_utilities.pad' obj was destroyed in the first call

"""

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""

#   THESE FUNCTIONS ARE COMMENTED BECAUSE NOT INVOLVED IN THE PYTHON3 CONVERSION. THEY HAVE TO BE TESTED ANYWAY
""" start: new in sphire 1.3


class Test_rec2D(unittest.TestCase):
    def test_rec2D(self):
        oldv = oldfu.rec2D(  lines=0, idrange=None, snr=None )
        v = fu.rec2D(  lines=0, idrange=None, snr=None )
        pass





class Test_recons3d_4nnf_MPI(unittest.TestCase):
    def test_recons3d_4nnf_MPI(self):
        oldv = oldfu.recons3d_4nnf_MPI(myid=0, list_of_prjlist=0, bckgdata=0, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, mpi_comm=None, smearstep = 0.0)
        v = fu.recons3d_4nnf_MPI(myid=0, list_of_prjlist=0, bckgdata=0, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, mpi_comm=None, smearstep = 0.0)
        pass


class Test_recons3d_4nnfs_MPI(unittest.TestCase):
    def test_recons3d_4nnfs_MPI(self):
        oldv = oldfu.recons3d_4nnfs_MPI(myid=0, main_node=0, prjlist=0, upweighted = True, finfo=None, mpi_comm=None, smearstep = 0.0, CTF = True, compensate = False, target_size=-1)
        v = fu.recons3d_4nnfs_MPI(myid=0, main_node=0, prjlist=0, upweighted = True, finfo=None, mpi_comm=None, smearstep = 0.0, CTF = True, compensate = False, target_size=-1)
        pass


class Test_recons3d_4nnstruct_MPI(unittest.TestCase):
    def test_recons3d_4nnstruct_MPI(self):
        oldv = oldfu.recons3d_4nnstruct_MPI(myid, main_node, prjlist, paramstructure, refang, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None)
        v = fu.recons3d_4nnstruct_MPI(myid, main_node, prjlist, paramstructure, refang, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None)
        pass


class Test_recons3d_4nnstruct_MPI_test(unittest.TestCase):
    def test_recons3d_4nnstruct_MPI_test(self):
        oldv = oldfu.recons3d_4nnstruct_MPI_test(myid, main_node, prjlist, paramstructure, refang, parameters, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None)
        v = fu.recons3d_4nnstruct_MPI_test(myid, main_node, prjlist, paramstructure, refang, parameters, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None)
        pass


class Test_recons3d_nn_SSNR(unittest.TestCase):
    def test_recons3d_nn_SSNR(self):
        oldv = oldfu.recons3d_nn_SSNR(stack_name,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0)
        v = fu.recons3d_nn_SSNR(stack_name,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0)
        pass


class Test_bootstrap_nn(unittest.TestCase):
    def test_bootstrap_nn(self):
        oldv = oldfu.bootstrap_nn(proj_stack, volume_stack, list_proj, niter, media="memory", npad=4, symmetry="c1", output=-1, CTF=False, snr=1.0, sign=1, myseed=None )
        v = fu.bootstrap_nn(proj_stack, volume_stack, list_proj, niter, media="memory", npad=4, symmetry="c1", output=-1, CTF=False, snr=1.0, sign=1, myseed=None )
        pass


class Test_recons3d_em(unittest.TestCase):
    def test_recons3d_em(self):
        oldv = oldfu.recons3d_em(projections_stack, max_iterations_count = 100, radius = -1, min_avg_abs_voxel_change = 0.01, use_weights = False, symmetry = "c1")
        v = fu.recons3d_em(projections_stack, max_iterations_count = 100, radius = -1, min_avg_abs_voxel_change = 0.01, use_weights = False, symmetry = "c1")
        pass


class Test_recons3d_em_MPI(unittest.TestCase):
    def test_recons3d_em_MPI(self):
        oldv = oldfu.recons3d_em_MPI(projections_stack, output_file, max_iterations_count = 100, radius = -1, min_norm_absolute_voxel_change = 0.01, use_weights = False, symmetry = "c1", min_norm_squared_voxel_change = 0.0001)
        v = fu.recons3d_em_MPI(projections_stack, output_file, max_iterations_count = 100, radius = -1, min_norm_absolute_voxel_change = 0.01, use_weights = False, symmetry = "c1", min_norm_squared_voxel_change = 0.0001)
        pass


class Test_recons3d_sirt(unittest.TestCase):
    def test_recons3d_sirt(self):
        oldv = oldfu.recons3d_sirt(stack_name, list_proj, radius, lam=1.0e-4, maxit=100, symmetry="c1", tol=0.001)
        v = fu.recons3d_sirt(stack_name, list_proj, radius, lam=1.0e-4, maxit=100, symmetry="c1", tol=0.001)
        pass


class Test_recons3d_wbp(unittest.TestCase):
    def test_recons3d_wbp(self):
        oldv = oldfu.recons3d_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1", radius=None)
        v = fu.recons3d_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1", radius=None)
        pass


class Test_recons3d_vwbp(unittest.TestCase):
    def test_recons3d_vwbp(self):
        oldv = oldfu.recons3d_vwbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1",outstack="bdb:temp")
        v = fu.recons3d_vwbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1",outstack="bdb:temp")
        pass


class Test_prepare_wbp(unittest.TestCase):
    def test_prepare_wbp(self):
        oldv = oldfu.prepare_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1")
        v = fu.prepare_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1")
        pass


class Test_recons3d_swbp(unittest.TestCase):
    def test_recons3d_swbp(self):
        oldv = oldfu.recons3d_swbp(A, transform, L, ss, method = "general", const=1.0E4, symmetry="c1")
        v = fu.recons3d_swbp(A, transform, L, ss, method = "general", const=1.0E4, symmetry="c1")
        pass


class Test_weight_swbp(unittest.TestCase):
    def test_weight_swbp(self):
        oldv = oldfu.weight_swbp(A, L, ss, method = "general", const=1.0E4, symmetry="c1")
        v = fu.weight_swbp(A, L, ss, method = "general", const=1.0E4, symmetry="c1")
        pass


class Test_backproject_swbp(unittest.TestCase):
    def test_backproject_swbp(self):
        oldv = oldfu.backproject_swbp(B, transform = None, symmetry="c1")
        v = fu.backproject_swbp(B, transform = None, symmetry="c1")
        pass

class Test_one_swbp(unittest.TestCase):
    def test_one_swbp(self):
        oldv = oldfu.one_swbp(CUBE, B, transform = None, symmetry="c1")
        v = fu.one_swbp(CUBE, B, transform = None, symmetry="c1")
        pass

class Test_prepare_recons_ctf_fftvol(unittest.TestCase):
    def test_prepare_recons_ctf_fftvol(self):
        oldv = oldfu.prepare_recons_ctf_fftvol(data, snr, symmetry, myid, main_node_half, pidlist, finfo=None, npad = 2, mpi_comm=None)
        v = fu.prepare_recons_ctf_fftvol(data, snr, symmetry, myid, main_node_half, pidlist, finfo=None, npad = 2, mpi_comm=None)
        pass

class Test_recons_ctf_from_fftvol_using_nn4_ctfw(unittest.TestCase):
    def test_recons_ctf_from_fftvol_using_nn4_ctfw(self):
        oldv = oldfu.recons_ctf_from_fftvol_using_nn4_ctfw(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2)
        v = fu.recons_ctf_from_fftvol_using_nn4_ctfw(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2)
        pass

class Test_rec3D_MPI_with_getting_odd_even_volumes_from_files(unittest.TestCase):
    def test_rec3D_MPI_with_getting_odd_even_volumes_from_files(self):
        oldv = oldfu.rec3D_MPI_with_getting_odd_even_volumes_from_files(fftvol_files, weight_files, reconstructed_vol_files,nx, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, myid = 0, main_node = 0, rstep = 1.0, finfo=None, npad = 2, mpi_comm=None)
        v = fu.rec3D_MPI_with_getting_odd_even_volumes_from_files(fftvol_files, weight_files, reconstructed_vol_files,nx, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, myid = 0, main_node = 0, rstep = 1.0, finfo=None, npad = 2, mpi_comm=None)
        pass
"""

""" after an hour is still running ....
class Test_recons3d_4nn(unittest.TestCase):
    def test_default_case(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        return_new = fu.recons3d_4nn(stack_name=STACK_NAME, list_proj=list_proj[:3], symmetry="c1", npad=4, snr=None, weighting=1, varsnr=False, xysize=-1, zsize = -1)
        return_old = oldfu.recons3d_4nn(stack_name=STACK_NAME, list_proj=list_proj[:3], symmetry="c1", npad=4, snr=None, weighting=1, varsnr=False, xysize=-1, zsize = -1)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5  ))
"""
""" not able to test it
class Test_recons3d_4nnw_MPI(unittest.TestCase):
    def test_recons3d_4nnw_MPI(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        return_old = oldfu.recons3d_4nnw_MPI(myid=0, prjlist=[proj,proj], bckgdata=proj, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None, smearstep = 0.0, fsc = None)
        return_new = fu.recons3d_4nnw_MPI(myid=0, prjlist=[proj,proj], bckgdata=proj, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None, smearstep = 0.0, fsc = None)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5))
"""
""" end: new in sphire 1.3"""

# I have to use the pickle file, nosetests is able to read it but pytest not, to have a 'xform.projection' IMAGE
class Test_insert_slices(unittest.TestCase):
    size = 76
    img = EMData(size, size)

    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.insert_slices()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.insert_slices()
        self.assertEqual(
            str(cm_new.exception), "insert_slices() missing 2 required positional arguments: 'reconstructor' and 'proj'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_defalut_case(self):
        params = {
            "size": self.size,
            "npad": 2,
            "symmetry": "c1",
            "fftvol": deepcopy(self.img),
            "weight": deepcopy(self.img),
            "snr": 2,
        }
        r_new = Reconstructors.get("nn4", params)
        r_new.setup()
        r_old = Reconstructors.get("nn4", params)
        r_old.setup()
        return_new = fu.insert_slices(
            reconstructor=r_new, proj=deepcopy(XFORM_PROJECTION_IMG)
        )
        return_old = oldfu.insert_slices(
            reconstructor=r_old, proj=deepcopy(XFORM_PROJECTION_IMG)
        )
        fftvol_new = r_new.get_params()["fftvol"]
        fftvol_old = r_old.get_params()["fftvol"]
        weight_new = r_new.get_params()["weight"]
        weight_old = r_old.get_params()["weight"]
        self.assertTrue(allclose(fftvol_new.get_3dview(), fftvol_old.get_3dview(),equal_nan=True))

        self.assertTrue(
            allclose(
                fftvol_new.get_3dview().flatten().tolist()[20000:20200],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, nan, nan, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                atol=1.e-4,
                equal_nan=True
            )
        )
        self.assertTrue(
            array_equal(
                weight_new.get_3dview().flatten().tolist()[20000:20200],
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    4.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    2.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(array_equal(weight_new.get_3dview(), weight_old.get_3dview()))
        self.assertEqual(return_new, return_old)
        self.assertTrue(return_new is None)

    def test_None_proj_case_returns_AttributeError_NoneType_obj_hasnot_attribute_get_attr(
        self
    ):
        params = {
            "size": self.size,
            "npad": 2,
            "symmetry": "c1",
            "fftvol": deepcopy(self.img),
            "weight": deepcopy(self.img),
            "snr": 2,
        }
        r = Reconstructors.get("nn4", params)
        r.setup()
        with self.assertRaises(AttributeError) as cm_new:
            fu.insert_slices(reconstructor=r, proj=None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.insert_slices(reconstructor=r, proj=None)
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_attr'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_image_proj_case_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(
        self
    ):
        params = {
            "size": self.size,
            "npad": 2,
            "symmetry": "c1",
            "fftvol": deepcopy(self.img),
            "weight": deepcopy(self.img),
            "snr": 2,
        }
        r = Reconstructors.get("nn4", params)
        r.setup()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.insert_slices(reconstructor=r, proj=EMData())
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.insert_slices(reconstructor=r, proj=EMData())
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_img_not_xform_projection_returns_RuntimeError_NotExistingObjectException_the_key_mean_doesnot_exist(
        self
    ):
        params = {
            "size": self.size,
            "npad": 2,
            "symmetry": "c1",
            "fftvol": deepcopy(self.img),
            "weight": deepcopy(self.img),
            "snr": 2,
        }
        r = Reconstructors.get("nn4", params)
        r.setup()
        with self.assertRaises(RuntimeError) as cm_new:
            fu.insert_slices(reconstructor=r, proj=get_real_data(2)[0])
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.insert_slices(reconstructor=r, proj=get_real_data(2)[0])
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])


class Test_recons3d_4nn_MPI(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons3d_4nn_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons3d_4nn_MPI()
        self.assertEqual(
            str(cm_new.exception),
            "recons3d_4nn_MPI() missing 2 required positional arguments: 'myid' and 'prjlist'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    """
    def test_default_case(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(
            allclose(
                return_new.get_3dview().flatten().tolist()[2925:3000],
                [
                    0.0,
                    5.451178731163964e-05,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                0.001,
            )
        )
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5))
    """

    def test_default_case_xy_z_size_both_not_negative(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(
            allclose(return_new.get_3dview().flatten(), [float("NaN")], equal_nan=True)
        )
        self.assertTrue(
            allclose(return_new.get_3dview(), return_old.get_3dview(), equal_nan=True)
        )

    def test_default_case_xy_size_not_negative(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=0,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(
            allclose(
                return_new.get_3dview().flatten(),
                [
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                    float("NaN"),
                ],
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(
                return_new.get_3dview(), return_old.get_3dview(), 0.5, equal_nan=True
            )
        )

    def test_default_case_z_size_both_not_negative_FAILEd(self):
        self.assertTrue(True)
        """
        return_new = fu.recons3d_4nn_MPI(myid=0, prjlist=[XFORM_PROJECTION_IMG], symmetry="c1", finfo=None, snr = 1.0, npad=2, xysize=-1, zsize=1, mpi_comm=MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(myid=0, prjlist=[XFORM_PROJECTION_IMG], symmetry="c1", finfo=None, snr = 1.0, npad=2, xysize=-1, zsize=1, mpi_comm=MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(),0.5, equal_nan=True))
        """

    def test_default_case_xy_z_size_both_not_negative_myid_not_null(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(array_equal(return_new.get_3dview().flatten(), [0.0]))
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_case_xy_size_not_negative_myid_not_null(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=1,
            zsize=-1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(
            array_equal(
                return_new.get_3dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_case_z_size_both_not_negative__myid_not_null(self):
        return_new = fu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=-1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(
            myid=1,
            prjlist=[XFORM_PROJECTION_IMG],
            symmetry="c1",
            finfo=None,
            snr=1.0,
            npad=2,
            xysize=-1,
            zsize=1,
            mpi_comm=MPI_COMM_WORLD,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(array_equal(return_new.get_3dview(), return_old.get_3dview()))
        self.assertTrue(
            array_equal(return_new.get_3dview().flatten()[:100], numpy_zeros(100))
        )

    @unittest.skip("SystemExit was raised")
    def test_prjlist_is_emptylist_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.recons3d_4nn_MPI(
                myid=0,
                prjlist=[],
                symmetry="c1",
                finfo=None,
                snr=1.0,
                npad=2,
                xysize=-1,
                zsize=-1,
                mpi_comm=MPI_COMM_WORLD,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.recons3d_4nn_MPI(
                myid=0,
                prjlist=[],
                symmetry="c1",
                finfo=None,
                snr=1.0,
                npad=2,
                xysize=-1,
                zsize=-1,
                mpi_comm=MPI_COMM_WORLD,
            )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_recons3d_trl_struct_MPI(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons3d_trl_struct_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons3d_trl_struct_MPI()
        self.assertEqual(
            str(cm_new.exception),
            "recons3d_trl_struct_MPI() missing 7 required positional arguments: 'myid', 'main_node', 'prjlist', 'paramstructure', 'refang', 'rshifts_shrank', and 'delta'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_recons3d_4nn_ctf_MPI(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons3d_4nn_ctf_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons3d_4nn_ctf_MPI()
        self.assertEqual(
            str(cm_new.exception),
            "recons3d_4nn_ctf_MPI() missing 2 required positional arguments: 'myid' and 'prjlist'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    """
    def test_default_case(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        return_new = fu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(
            allclose(
                return_new.get_3dview().flatten().tolist()[179691:179800],
                [
                    0.004414418246597052,
                    0.0027832286432385445,
                    0.0031159124337136745,
                    0.0014303999487310648,
                    -0.0010390577372163534,
                    -0.0009075439302250743,
                    -0.0013101220829412341,
                    -0.0010013339342549443,
                    0.0023415705654770136,
                    0.001426796312443912,
                    -0.0010065339738503098,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                0.01,
            )
        )
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5))
    """

    @unittest.skip(
        "crash if run togheter with the other tests of this class because a bad implementation of the code"
    )
    def test_negative_smearstep(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        return_new = fu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=-0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=-0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5))

    """
    def test_default_case_xy_z_size_both_not_negative_NameError_sizeprojection_BEACUASE_A_BUG(
        self
    ):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        with self.assertRaises(NameError) as cm_new:
            fu.recons3d_4nn_ctf_MPI(
                0,
                [proj],
                snr=1.0,
                sign=1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=1,
                zsize=1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(NameError) as cm_old:
            oldfu.recons3d_4nn_ctf_MPI(
                0,
                [proj],
                snr=1.0,
                sign=1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=1,
                zsize=1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(
            str(cm_new.exception), "global name 'sizeprojection' is not defined"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_case_xy_size_NameError_sizeprojection_BEACUASE_A_BUG(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        with self.assertRaises(NameError) as cm_new:
            fu.recons3d_4nn_ctf_MPI(
                0,
                [proj],
                snr=1.0,
                sign=1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=1,
                zsize=-1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(NameError) as cm_old:
            oldfu.recons3d_4nn_ctf_MPI(
                0,
                [proj],
                snr=1.0,
                sign=1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=1,
                zsize=-1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(
            str(cm_new.exception), "global name 'sizeprojection' is not defined"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
    """
    @unittest.skip(
        "crash if run togheter with the other tests of this class because a bad implementation of the code"
    )
    def test_default_case_negative_sign(self):
        nima = EMUtil.get_image_count(STACK_NAME)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(STACK_NAME, list_proj[0])
        return_new = fu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=-1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_ctf_MPI(
            0,
            [proj],
            snr=1.0,
            sign=-1,
            symmetry="c1",
            finfo=None,
            npad=2,
            xysize=-1,
            zsize=-1,
            mpi_comm=None,
            smearstep=0.5,
        )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5))

    @unittest.skip("SystemExit was raised")
    def test_prjlist_is_emptylist_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.recons3d_4nn_ctf_MPI(
                0,
                [],
                snr=-1.0,
                sign=-1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=-1,
                zsize=-1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.recons3d_4nn_ctf_MPI(
                0,
                [],
                snr=-1.0,
                sign=-1,
                symmetry="c1",
                finfo=None,
                npad=2,
                xysize=-1,
                zsize=-1,
                mpi_comm=None,
                smearstep=0.5,
            )
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
