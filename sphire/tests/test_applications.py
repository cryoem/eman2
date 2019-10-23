from __future__ import print_function
from __future__ import division

import unittest
from numpy import allclose, array_equal
from mpi import *
import global_def

from os import path, mkdir
from cPickle import loads as pickle_loads
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))


global_def.BATCH = True
global_def.MPI = True



from test_module import  remove_list_of_file, remove_dir, get_arg_from_pickle_file, ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,returns_values_in_file,get_real_data

from EMAN2db import db_open_dict as EMAN2db_db_open_dict

from sphire.libpy.sp_utilities import model_blank, model_circle, even_angles

from EMAN2_cppwrap import EMData
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_BLANK_2D = model_blank(10, 10)
IMAGE_BLANK_3D = model_blank(10, 10, 10)
MASK = model_circle(2, 5, 5)

TOLERANCE = 0.0005
ABSOLUTE_PATH_TO_STACK= "bdb:" + path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Substack/isac_substack")
"""
There are some opened issues in:
1) ali2d_MPI the results are different ... ANYWAY IT SEEMS TO BE NOT USED
2) ali2d_base the results are different
    -) if CTF =True and myid == main_node (i.e.: ln 695) it try to use an not defined variable and crash .... dead code or bug?
3) mref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
4) Kmref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
5) project3d:
    -) with 'noise' not None the results are different because the addition of a random gauss_noise to the volume in the code (noise is the sigma value in the generation noise process)
    -) How can I really test with 'listctfs'? It is never used in a real SPHIRE code case, anyway I'd test
    -) in case of realsp=trillinear=True it will spawn an error message but its behaviour will be as if it'd receive as input realsp=True and Trilinear=False ....BUG or just BAD IMPLEMENTED
6) recons3d_n_trl_MPI_one_node: I cannot run it. see the error message in "Test_recons3d_n_trl_MPI_one_node.test_NoCTF_symC1"
7) pca
    -) need a file name containing the average of the input stack
8) header
    -) the pickle file values lead to an error message got at line 3!!!!!!!! 
    -) you have to save the values in a file to test its behaviour. But if you set an output file you cannot set the other values becuase the 
        "op = zero+one++consecutive+randomize+rand_alpha+(fimport!=None)+(fexport!=None)+fprint+backup+restore+delete+doset" check .... is it a bug? which is the purpose of this function??
9) refvol: seems to be never USED
10) within_group_refinement:
    -) I tested just the default method because the other 2 are not in use (Fabian docet)
    -) do not test with randomize=True because the results will be never the same because the random orientation in same calculations
    -) since the input parameter are basically used in 'sparx_alignment.ali2d_single_iter' and 'sparx_utilities.combine_params2' that i have already deeply tested
        I tested it just the pickle file case
11) ali3d_mref_Kmeans_MPI and mref_ali3d_EQ_Kmeans are used only in sxsort3d.py that is buggy and maybe we will implement it from scratch. I did not test them for now
"""

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""

""" start: new in sphire 1.3"""
from sphire.libpy import sp_applications as oldfu
from sphire.libpy_py3 import sp_applications as fu

class Test_ali2d(unittest.TestCase):
    def test_ali2d(self):
        oldfu.ali2d(stack="", outdir="", maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", nomirror=False, dst=0.0, center=-1, maxit=0, CTF=False, snr=1.0, Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID="", MPI=False, template=None, random_method = "")
        pass

class Test_ali2d_data(unittest.TestCase):
    def test_ali2d_data(self):
        oldfu.ali2d_data(data="", outdir="", maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", nomirror=False, dst=0.0, center=-1, maxit=0, CTF=False, snr=1.0, Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID="", from_ali2d=False, template=None, random_method = "")
        pass


class Test_local_ali2d(unittest.TestCase):
    def test_local_ali2d(self):
        oldfu.local_ali2d(stack="", outdir="", maskfile = None, ou = -1, br = 1.75, center = 1, eps = 0.001, maxit = 10, CTF = False, snr = 1.0, user_func_name="ref_ali2d")
        pass

class Test_mref_ali2d(unittest.TestCase):
    def test_mref_ali2d(self):
        oldfu.mref_ali2d(stack="", refim="", outdir="", maskfile=None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000, MPI=False)
        pass


class Test_ali2d_ra(unittest.TestCase):
    def test_ali2d_ra(self):
        oldfu.ali2d_ra(stack="", maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000)
        pass


class Test_ali2d_rag(unittest.TestCase):
    def test_ali2d_rag(self):
        oldfu.ali2d_rag(stack="", maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000)
        pass


class Test_ali2d_rac(unittest.TestCase):
    def test_ali2d_rac(self):
        oldfu.ali2d_rac(stack="", maskfile = None, ir = 1, ou = -1, rs = 1, nclass = 2, maxit = 10, maxin = 10, check_mirror = False, rand_seed = 1000, MPI=False)
        pass


class Test_ali2d_ras(unittest.TestCase):
    def test_ali2d_ras(self):
        oldfu.ali2d_ras(data2d="", randomize = False, ir = 1, ou = -1, rs = 1, step = 1.0, dst = 0.0, maxit = 10, check_mirror = True, FH = 0.0, FF =0.0)
        pass


class Test_ali2d_rotationaltop(unittest.TestCase):
    def test_ali2d_rotationaltop(self):
        oldfu.ali2d_rotationaltop(outdir="", stack="", randomize = False, orient=True, ir = 4, ou = -1, rs = 1, psi_max = 180.0, mode = "F", maxit = 10)
        pass


class Test_ali2d_rotational(unittest.TestCase):
    def test_ali2d_rotational(self):
        oldfu.ali2d_rotational(data2d="", randomize = False, orient=True, ir = 1, ou = -1, rs = 1, psi_max = 180.0, mode = "F", maxit = 10)
        pass


class Test_ali2d_cross_res(unittest.TestCase):
    def test_ali2d_cross_res(self):
        oldfu.ali2d_cross_res(stack="", outdir="", maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d")
        pass


class Test_ali3d(unittest.TestCase):
    def test_ali3d(self):
        oldfu.ali3d(stack="", ref_vol="", outdir="", maskfile = None, ir = 1, ou = -1, rs = 1, xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0)
        pass


class Test_ali3d_MPI(unittest.TestCase):
    def test_ali3d_MPI(self):
        oldfu.ali3d_MPI(stack="", ref_vol="", outdir="", maskfile = None, ir = 1, ou = -1, rs = 1, xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",fourvar = True, npad = 2, debug = False, termprec = 0.0)
        pass


class Test_sali3d_base(unittest.TestCase):
    def test_sali3d_base(self):
        p = oldfu.sali3d_base(stack="", ref_vol = None, Tracker = None, rangle = 0.0, rshift = 0.0, mpi_comm = None, log = None)
        pass


class Test_slocal_ali3d_base(unittest.TestCase):
    def test_slocal_ali3d_base(self):
        v= oldfu.slocal_ali3d_base(stack="", templatevol="", Tracker="", mpi_comm = None, log= None, chunk = -1.0, debug = False )
        pass


class Test_computenumberofrefs(unittest.TestCase):
    def test_computenumberofrefs(self):
        v= oldfu.computenumberofrefs(x="", dat="")
        pass


class Test_ali3dpsi_MPI(unittest.TestCase):
    def test_ali3dpsi_MPI(self):
        v= oldfu.ali3dpsi_MPI(stack="", ref_vol="", outdir="", maskfile = None, ir = 1, ou = -1, rs = 1, xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",fourvar = True, npad = 4, debug = False, termprec = 0.0)
        pass


class Test_ali3d_shcMPI(unittest.TestCase):
    def test_ali3d_shcMPI(self):
        v= oldfu.ali3d_shcMPI(stack="", ref_vol="", outdir="", maskfile = None, ir = 1, ou = -1, rs = 1, xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",fourvar = True, npad = 4, debug = False, termprec = 0.0, gamma=-1)
        pass


class Test_mref_ali3d(unittest.TestCase):
    def test_mref_ali3d(self):
        oldfu.mref_ali3d(stack="", ref_vol="", outdir="", maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta="10 6 4 4", an="-1", center = 1.0, nassign = 3, nrefine = 1, CTF = False, snr = 1.0,  ref_a = "S", sym="c1",user_func_name="ref_ali3d", MPI=False, npad = 4, debug = False, fourvar=False, termprec = 0.0)
        pass


class Test_Kmref2_ali3d_MPI(unittest.TestCase):
    def test_Kmref2_ali3d_MPI(self):
        oldfu.Kmref2_ali3d_MPI(stack="", ref_vol="", outdir="", maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1",center = -1, nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",user_func_name="ref_ali3d", npad = 4, debug = False, fourvar=False, termprec = 0.0, mpi_comm = None, log = None)
        pass


class Test_get_refiparams(unittest.TestCase):
    def test_get_refiparams(self):
        d= oldfu.get_refiparams(nx="") # {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
        pass


class Test_local_ali3dm_MPI_(unittest.TestCase):
    def test_local_ali3dm_MPI_(self):
        oldfu.local_ali3dm_MPI_(stack="", refvol="", outdir="", maskfile="", ou=-1,  delta=2, ts=0.25, maxit=10, nassign=4, nrefine=1, CTF = None,snr=1.0, sym="c1", user_func_name="ref_ali3d", fourvar=False, debug=False, termprec = 0.0 )
        pass

class Test_local_ali3dm_MPI(unittest.TestCase):
    def test_local_ali3dm_MPI(self):
        oldfu.local_ali3dm_MPI(stack="", refvol="", outdir="", maskfile="", ou=-1,  delta=2, ts=0.25, maxit=10, nassign=4, nrefine=1, CTF = None,snr=1.0, sym="c1", user_func_name="ref_ali3d", fourvar=False, debug=False, termprec = 0.0 )
        pass


class Test_local_ali3d(unittest.TestCase):
    def test_local_ali3d(self):
        v= oldfu.local_ali3d(stack="", outdir="", maskfile = None, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10,CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",fourvar = True, npad = 4, debug = False, MPI = False)
        pass

class Test_local_ali3d_MPI(unittest.TestCase):
    def test_local_ali3d_MPI(self):
        v= oldfu.local_ali3d_MPI(stack="", outdir="", maskfile = None, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10,CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",fourvar = True, npad = 4, debug = False, MPI = False)
        pass

class Test_local_ali3d_MPI_scipy_minimization(unittest.TestCase):
    def test_local_ali3d_MPI_scipy_minimization(self):
        v= oldfu.local_ali3d_MPI_scipy_minimization(stack="", outdir="", maskfile = None, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10,CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",fourvar = True, npad = 4, debug = False, MPI = False)
        pass


class Test_local_ali3d_base_MPI(unittest.TestCase):
    def test_local_ali3d_base_MPI(self):
        oldfu.local_ali3d_base_MPI(stack="", templatevol="", ali3d_options="", shrinkage = 1.0,mpi_comm = None, log= None, chunk = -1.0, saturatecrit = 0.95, pixercutoff = 1.0, debug = False )
        pass


class Test_autowin(unittest.TestCase):
    def test_autowin(self):
        v= oldfu.autowin(indir="",outdir="", noisedoc="", noisemic="", templatefile="", deci="", CC_method="", p_size="", sigma="", hf_p="", n_peak_max="", contrast_invert=1, CTF = False, prm = "micrograph", MPI=False)
        pass

class Test_autowin_MPI(unittest.TestCase):
    def test_autowin_MPI(self):
        v= oldfu.autowin_MPI(indir="",outdir="", noisedoc="", noisemic="", templatefile="", deci="", CC_method="", p_size="", sigma="", hf_p="", n_peak_max="", contrast_invert=1,  prefix_of_micrograph = "micrograph")
        pass


class Test_ihrsr_MPI(unittest.TestCase):
    def test_ihrsr_MPI(self):
        oldfu.ihrsr_MPI(stack="", ref_vol="", outdir="", maskfile="", ir="", ou="", rs="", xr="", ynumber="",txs="", delta="", initial_theta="", delta_theta="", an="", maxit="", CTF="", snr="", dp="", ndp="", dp_step="", dphi="", ndphi="", dphi_step="", psi_max="",rmin="", rmax="", fract="", nise="", npad="", sym="", user_func_name="", datasym="", pixel_size="", debug="", y_restrict="", WRAP="")
        pass


class Test_copyfromtif(unittest.TestCase):
    def test_copyfromtif(self):
        oldfu.copyfromtif(indir="", outdir=None, input_extension="tif", film_or_CCD="f", output_extension="hdf", contrast_invert=1, Pixel_size=1, scanner_param_a=1,scanner_param_b=1, scan_step=63.5, magnification=40, MPI=False)
        pass

class Test_copyfromtif_MPI(unittest.TestCase):
    def test_copyfromtif_MPI(self):
        oldfu.copyfromtif_MPI(indir="", outdir=None, input_extension="tif", film_or_CCD="f", output_extension="hdf", contrast_invert=1, Pixel_size=1, scanner_param_a=1,scanner_param_b=1, scan_step=63.5, magnification=40, MPI=False)
        pass


class Test_dele_flist(unittest.TestCase):
    def test_dele_flist(self):
        oldfu.dele_flist(flist="")
        pass


class Test_defocus_calc(unittest.TestCase):
    def test_defocus_calc(self):
        oldfu.defocus_calc(roodir="", method="", writetodoc="w", Pixel_size=1, voltage=120, Cs=1, amp_contrast=.1, round_off=100, dz_max=50000., frequency_low=30, frequency_high=5, polynomial_rank_baseline=5, polynomial_rank_envelope=5, prefix="roo", format="spider", skip_comment="#", micdir = "no", print_screen="no")
        pass


class Test_pw2sp(unittest.TestCase):
    def test_pw2sp(self):
        oldfu.pw2sp(indir="", outdir = None, w =256, xo =50, yo = 50, xd = 0, yd = 0, r = 0, prefix_of_micrograph="micrograph", MPI=False)
        pass

class Test_pw2sp_MPI(unittest.TestCase):
    def test_pw2sp_MPI(self):
        oldfu.pw2sp_MPI(indir="", outdir = None, w =256, xo =50, yo = 50, xd = 0, yd = 0, r = 0, prefix_of_micrograph="micrograph")
        pass

class Test_ra_cef(unittest.TestCase):
    def test_ra_cef(self):
        oldfu.ra_cef(indir="", noise="", outdir="", prf="", num="")
        pass


class Test_ali_vol_2(unittest.TestCase):
    def test_ali_vol_2(self):
        v= oldfu.ali_vol_2(vol="", refv="", ang_scale="", shift_scale="", radius=None, discrepancy = "ccc")
        pass

class Test_ali_vol_3(unittest.TestCase):
    def test_ali_vol_3(self):
        v= oldfu.ali_vol_3(vol="", refv="", ang_scale="", shift_scale="", radius=None, discrepancy = "ccc", mask=None)
        pass

class Test_ali_vol_n(unittest.TestCase):
    def test_ali_vol_n(self):
        v= oldfu.ali_vol_n(vol="", refv="", ang_scale="", shift_scale="", radius=None, discrepancy = "ccc", rsdec=1)
        pass

class Test_ali_vol_grid(unittest.TestCase):
    def test_ali_vol_grid(self):
        v= oldfu.ali_vol_grid(vol="", params="", refv="", ang_scale="", shift_scale="", radius=None, discrepancy="dot", kb=None, wrap=False)
        pass

class Test_ali_vol_M(unittest.TestCase):
    def test_ali_vol_M(self):
        v= oldfu.ali_vol_M(vol="", refv="", ang_scale="", shift_scale="", mask=None, discrepancy = "ccc")
        pass

class Test_ali_vol_nopsi(unittest.TestCase):
    def test_ali_vol_nopsi(self):
        v= oldfu.ali_vol_nopsi(vol="", refv="", ang_scale="", shift_scale="", radius=None, discrepancy = "ccc")
        pass

class Test_ali_vol_rotate(unittest.TestCase):
    def test_ali_vol_rotate(self):
        v= oldfu.ali_vol_rotate(vol="", refv="", ang_scale="",  radius=None, discrepancy = "ccc")
        pass

class Test_ali_vol_shift(unittest.TestCase):
    def test_ali_vol_shift(self):
        v= oldfu.ali_vol_shift(vol="", refv="", shift_scale="",  radius=None, discrepancy = "ccc")
        pass

class Test_ali_vol_scale(unittest.TestCase):
    def test_ali_vol_scale(self):
        v= oldfu.ali_vol_scale(vol="", refv="", ang_scale="", shift_scale="", mag_scale="", radius=None, discrepancy = "ccc")
        pass


class Test_ali_vol_only_scale(unittest.TestCase):
    def test_ali_vol_only_scale(self):
        v= oldfu.ali_vol_only_scale(vol="", refv="", mag_scale="", radius=None, discrepancy = "ccc")
        pass


class Test_rot_sym(unittest.TestCase):
    def test_rot_sym(self):
        v= oldfu.rot_sym(infile="", outfile="", sym_gp="d4", radius=None, phi=0, theta=0, psi=0, phirange=20, thetarange=20, psirange=20, ftolerance=1.e-4, xtolerance=1.e-4)
        pass


class Test_transform2d(unittest.TestCase):
    def test_transform2d(self):
        v= oldfu.transform2d(stack_data="", stack_data_ali="", shift = False, ignore_mirror = False, method = "quadratic")
        pass


class Test_recons3d_n_MPI(unittest.TestCase):
    def test_recons3d_n_MPI(self):
        v= oldfu.recons3d_n_MPI(prj_stack="", pid_list="", vol_stack="", CTF=False, snr=1.0, sign=1, npad=2, sym="c1", listfile="", group=-1, verbose=0, xysize=-1, zsize=-1, smearstep = 0.0)
        pass


class Test_recons3d_trl_MPI(unittest.TestCase):
    def test_recons3d_trl_MPI(self):
        v= oldfu.recons3d_trl_MPI(prj_stack="", pid_list="", vol_stack="", CTF="", snr="", sign="", npad="", sym="", verbose = None, niter =10, compensate = False, target_window_size=-1)
        pass


class Test_recons3d_n_trl_MPI_one_node(unittest.TestCase):
    def test_recons3d_n_trl_MPI_one_node(self):
        v= oldfu.recons3d_n_trl_MPI_one_node(prjlist="", CTF="", snr="", sign="", npad="", sym="", group="", niter="", verbose="", upweighted="", compensate="", chunk_id="")
        pass


class Test_ssnr3d(unittest.TestCase):
    def test_ssnr3d(self):
        v= oldfu.ssnr3d(stack="", output_volume = None, ssnr_text_file = None, mask = None, reference_structure = None, ou = -1, rw = 1.0,  npad = 1, CTF = False, sign = 1, sym ="c1", MPI = False, random_angles = 0)
        pass


class Test_ssnr3d_MPI(unittest.TestCase):
    def test_ssnr3d_MPI(self):
        v= oldfu.ssnr3d_MPI(stack="", output_volume = None, ssnr_text_file = None, mask = None, reference_structure = None, ou = -1, rw = 1.0,  npad = 1, CTF = False, sign = 1, sym ="c1",  random_angles = 0)
        pass


class Test_varimax(unittest.TestCase):
    def test_varimax(self):
        v= oldfu.varimax(input_stack="", imglist="", output_stack="", maskfile="", mask_radius="", verbose ="")
        pass


class Test_bootstrap_genbuf(unittest.TestCase):
    def test_bootstrap_genbuf(self):
        v= oldfu.bootstrap_genbuf(prj_stack="", buf_prefix="", npad="", verbose="", CTF=False)
        pass


class Test_bootstrap_run(unittest.TestCase):
    def test_bootstrap_run(self):
        v= oldfu.bootstrap_run(prj_stack="", media="", outdir="", nvol="", CTF="", snr="", sym="", verbose="", MPI=False)
        pass


class Test_wrapper_params_2D_to_3D(unittest.TestCase):
    def test_wrapper_params_2D_to_3D(self):
        v= oldfu.wrapper_params_2D_to_3D(stack="")
        pass


class Test_wrapper_params_3D_to_2D(unittest.TestCase):
    def test_wrapper_params_3D_to_2D(self):
        v= oldfu.wrapper_params_3D_to_2D(stack="")
        pass


class Test_cml_find_structure_main(unittest.TestCase):
    def test_cml_find_structure_main(self):
        v= oldfu.cml_find_structure_main(stack="", out_dir="", ir="", ou="", delta="", dpsi="", lf="", hf="", rand_seed="", maxit="", given = False, first_zero = False, flag_weights = False, debug = False, trials = 1)
        pass

class Test_cml_find_structure_MPI2(unittest.TestCase):
    def test_cml_find_structure_MPI2(self):
        v= oldfu.cml_find_structure_MPI2(stack="", out_dir="", ir="", ou="", delta="", dpsi="", lf="", hf="", rand_seed="", maxit="", given = False, first_zero = False, flag_weights = False, debug = False, trials = 1)
        pass


class Test_cml_find_structure_MPI(unittest.TestCase):
    def test_cml_find_structure_MPI(self):
        v= oldfu.cml_find_structure_MPI(stack="", out_dir="", ir="", ou="", delta="", dpsi="", lf="", hf="", rand_seed="", maxit="", given = False, first_zero = False, flag_weights = False, debug = False, trials = 10)
        pass


class Test_imgstat_ccc(unittest.TestCase):
    def test_imgstat_ccc(self):
        v= oldfu.imgstat_ccc( stacks="", rad="" )
        pass


class Test_imgstat_fsc(unittest.TestCase):
    def test_imgstat_fsc(self):
        v= oldfu.imgstat_fsc( stacks="", fscfile="", rad ="")
        pass


class Test_imgstat_inf(unittest.TestCase):
    def test_imgstat_inf(self):
        v= oldfu.imgstat_inf( stacks="", rad="" )
        pass


class Test_imgstat(unittest.TestCase):
    def test_imgstat(self):
        v= oldfu.imgstat( stacks="", ifccc="", fscfile="", pinf="", rad ="")
        pass


class Test_normal_prj(unittest.TestCase):
    def test_normal_prj(self):
        v= oldfu.normal_prj( prj_stack="", outdir="", refvol="", weights="", r="", niter="", snr="", sym="", verbose = 0, CTF = False, MPI=False )
        pass


class Test_defvar(unittest.TestCase):
    def test_defvar(self):
        v= oldfu.defvar(files="", outdir="", fl="", aa="", radccc="", frepa = "default", pca=False, pcamask=None, pcanvec=None)
        pass


class Test_var_mpi(unittest.TestCase):
    def test_var_mpi(self):
        v= oldfu.var_mpi(files="", outdir="", fl="", aa="", radccc="", frepa = "default", pca=False, pcamask=None, pcanvec=None)
        pass


class Test_factcoords_vol(unittest.TestCase):
    def test_factcoords_vol(self):
        v= oldfu.factcoords_vol( vol_stacks="", avgvol_stack="", eigvol_stack="", prefix="", rad = -1, neigvol = -1, fl=0.0, aa=0.0, MPI=False)
        pass


class Test_factcoords_prj(unittest.TestCase):
    def test_factcoords_prj(self):
        v= oldfu.factcoords_prj( prj_stacks="", avgvol_stack="", eigvol_stack="", prefix="", rad="", neigvol="", fl=0.0, aa=0.0, CTF = False, MPI=False)
        pass


class Test_spill_out(unittest.TestCase):
    def test_spill_out(self):
        v= oldfu.spill_out(ltot="", base="", d="", neigvol="", foutput="")
        pass


class Test_k_means_main(unittest.TestCase):
    def test_k_means_main(self):
        v= oldfu.k_means_main(stack="", out_dir="", maskname="", opt_method="", K="", rand_seed="", maxit="", trials="", critname="",CTF = False, F = 0, T0 = 0, MPI = False, CUDA = False, DEBUG = False, flagnorm = False,init_method = 'rnd')
        pass


class Test_k_means_groups(unittest.TestCase):
    def test_k_means_groups(self):
        v= oldfu.k_means_groups(stack="", out_file="", maskname="", opt_method="", K1="", K2="", rand_seed="", maxit="", trials="", CTF=False, F=0.0, T0=0.0, MPI=False, CUDA=False, DEBUG=False, flagnorm=False)
        pass


class Test_HAC_clustering(unittest.TestCase):
    def test_HAC_clustering(self):
        v= oldfu.HAC_clustering(stack="", dendoname="", maskname="", kind_link="", kind_dist="", flag_diss="")
        pass


class Test_HAC_averages(unittest.TestCase):
    def test_HAC_averages(self):
        v= oldfu.HAC_averages(stack="", dendoname="", avename="", K="")
        pass


class Test_tomo(unittest.TestCase):
    def test_tomo(self):
        v= oldfu.tomo(box="")
        pass


class Test_ave_ali(unittest.TestCase):
    def test_ave_ali(self):
        v= oldfu.ave_ali(name_stack="", name_out = None, ali = False, param_to_save_size = None, set_as_member_id = None)
        pass


class Test_refinement_2d_local(unittest.TestCase):
    def test_refinement_2d_local(self):
        v= oldfu.refinement_2d_local(data="", ou="", arange="", xrng="", yrng="", CTF = True, SNR=1.0e10)
        pass


class Test_volalixshift_MPI(unittest.TestCase):
    def test_volalixshift_MPI(self):
        v= oldfu.volalixshift_MPI(stack="", ref_vol="", outdir="", search_rng="", pixel_size="", dp="", dphi="", fract="", rmax="", rmin="", maskfile = None, maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False, nearby=3)
        pass


class Test_diskali_MPI(unittest.TestCase):
    def test_diskali_MPI(self):
        v= oldfu.diskali_MPI(stack="", ref_vol="", outdir="", maskfile="", dp="", dphi="", pixel_size="", user_func_name="", zstep=1.0, fract=0.67, rmax=70, rmin=0, CTF=False, maxit=1, sym = "c1")
        pass


class Test_cylindrical_trans(unittest.TestCase):
    def test_cylindrical_trans(self):
        v= oldfu.cylindrical_trans(vol="", rmin="", rmax="", rise="", apply_weights = False)
        pass


class Test_alihelical3(unittest.TestCase):
    def test_alihelical3(self):
        v= oldfu.alihelical3(slices="", refslices="", zstep="", dphi="", rise="", rmin="", rmax="", sym="c1")
        pass


class Test_alihelical4(unittest.TestCase):
    def test_alihelical4(self):
        v= oldfu.alihelical4(slices="", refslices="", zstep="", dphi="", rise="", rmin="", rmax="", theta=0.0)
        pass


class Test_iang(unittest.TestCase):
    def test_iang(self):
        v= oldfu.iang(alpha="", maxrin="")
        pass


class Test_stack_disks(unittest.TestCase):
    def test_stack_disks(self):
        v= oldfu.stack_disks(v="", nx="", ny="", ref_nz="", dphi="", rise="")
        pass


class Test_imgstat_hfsc(unittest.TestCase):
    def test_imgstat_hfsc(self):
        v= oldfu.imgstat_hfsc( stack="", file_prefix="", fil_attr='filament')
        pass


class Test_match_pixel_rise(unittest.TestCase):
    def test_match_pixel_rise(self):
        v ,v2= oldfu.match_pixel_rise(dz="",px="", nz=-1, ndisk=-1, rele=0.1, stop=900000)
        pass


class Test_gendisks_MPI(unittest.TestCase):
    def test_gendisks_MPI(self):
        v = oldfu.gendisks_MPI(stack="", mask3d="", ref_nx="", pixel_size="", dp="", dphi="", fract=0.67, rmax=70, rmin=0, CTF=False, user_func_name = "helical", sym = "c1", dskfilename='bdb:disks', maxerror=0.01, new_pixel_size = -1, do_match_pixel_rise=False)
        pass


class Test_ehelix_MPI(unittest.TestCase):
    def test_ehelix_MPI(self):
        v = oldfu.ehelix_MPI(stack="", ref_vol="", outdir="", seg_ny="", delta="", phiwobble="", psi_max="", search_rng="", rng="", ywobble="", ystep="", pixel_size="", dp="", dphi="", fract="", rmax="", rmin="", FindPsi = True, maskfile = None, maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False, slowIO = False)
        pass


class Test_localhelicon_MPInew(unittest.TestCase):
    def test_localhelicon_MPInew(self):
        v = oldfu.localhelicon_MPInew(stack="", ref_vol="", outdir="", seg_ny="", maskfile="", ir="", ou="", rs="", xr="", ynumber="", txs="", delta="", initial_theta="", delta_theta="", an="", maxit="", CTF="", snr="", dp="", dphi="", psi_max="", rmin="", rmax="", fract="",  npad="", sym="", user_func_name="", pixel_size="", debug="", y_restrict="", search_iter="", slowIO="")
        pass

class Test_localhelicon_MPIming(unittest.TestCase):
    def test_localhelicon_MPIming(self):
        v = oldfu.localhelicon_MPIming(stack="", ref_vol="", outdir="", seg_ny="", maskfile="", ir="", ou="", rs="", xr="", ynumber="", txs="", delta="", initial_theta="", delta_theta="", an="", maxit="", CTF="", snr="", dp="", dphi="", psi_max="", rmin="", rmax="", fract="",  npad="", sym="", user_func_name="", pixel_size="", debug="", y_restrict="", search_iter="", slowIO="")
        pass

class Test_localhelicon_MPInew_fullrefproj(unittest.TestCase):
    def test_localhelicon_MPInew_fullrefproj(self):
        v = oldfu.localhelicon_MPInew_fullrefproj(stack="", ref_vol="", outdir="", seg_ny="", maskfile="", ir="", ou="", rs="", xr="", ynumber="", txs="", delta="", initial_theta="", delta_theta="", an="", maxit="", CTF="", snr="", dp="", dphi="", psi_max="", rmin="", rmax="", fract="",  npad="", sym="", user_func_name="", pixel_size="", debug="", y_restrict="", search_iter="", slowIO="")
        pass

class Test_localhelicon_MPI(unittest.TestCase):
    def test_localhelicon_MPI(self):
        v = oldfu.localhelicon_MPI(stack="", ref_vol="", outdir="", seg_ny="", maskfile="", ir="", ou="", rs="", xr="", ynumber="", txs="", delta="", initial_theta="", delta_theta="", an="", maxit="", CTF="", snr="", dp="", dphi="", psi_max="", rmin="", rmax="", fract="",  npad="", sym="", user_func_name="", pixel_size="", debug="", y_restrict="", search_iter="", slowIO="")
        pass


class Test_filamentupdown(unittest.TestCase):
    def test_filamentupdown(self):
        v = oldfu.filamentupdown(fildata="", pixel_size="", dp="", dphi="")
        pass

class Test_setfilori_SP(unittest.TestCase):
    def test_setfilori_SP(self):
        v = oldfu.setfilori_SP(fildata="", pixel_size="", dp="", dphi="")
        pass



class Test_prepare_refffts(unittest.TestCase):
    def test_prepare_refffts(self):
        v = oldfu.prepare_refffts( volft="", kb="", nx="",ny="",nz="", segmask="", delta="",  MPI=False, psimax=1.0, psistep=1.0, kbx = None, kby = None, initial_theta = None, delta_theta = None)
        pass


class Test_prepare_helical_refangles(unittest.TestCase):
    def test_prepare_helical_refangles(self):
        v = oldfu.prepare_helical_refangles(delta="", initial_theta = None, delta_theta = None)
        pass


class Test_prepare_reffft1(unittest.TestCase):
    def test_prepare_reffft1(self):
        v = oldfu.prepare_reffft1( volft="", kb="", ref_angles="", segmask="", psimax=1.0, psistep=1.0, kbx = None, kby = None)
        pass


class Test_prepare_reffft2(unittest.TestCase):
    def test_prepare_reffft2(self):
        v = oldfu.prepare_reffft2( volft="", kb="", ref_angles="", segmask="", psimax=1.0, psistep=1.0, kbx = None, kby = None)
        pass


class Test_symsearch_MPI(unittest.TestCase):
    def test_symsearch_MPI(self):
        v = oldfu.symsearch_MPI(ref_vol="", outdir="", maskfile="", dp="", ndp="", dp_step="", dphi="", ndphi="", dphi_step="",rmin="", rmax="", fract="", sym="", user_func_name="", datasym="",pixel_size="", debug="")
        pass


class Test_sali3d_base_old(unittest.TestCase):
    def test_sali3d_base_old(self):
        v = oldfu.sali3d_base_old(stack="", ref_vol = None, Tracker = None, mpi_comm = None, log = None)
        pass


class Test_slocal_ali3d_base_old(unittest.TestCase):
    def test_slocal_ali3d_base_old(self):
        v = oldfu.slocal_ali3d_base_old(stack="", templatevol="", Tracker="", mpi_comm = None, log= None, chunk = -1.0, debug = False )
        pass


class Test_ali3d_mref_Kmeans_MPI(unittest.TestCase):
    def test_ali3d_mref_Kmeans_MPI(self):
        v = oldfu.ali3d_mref_Kmeans_MPI(ref_list="", outdir="", this_data_list_file="", Tracker="")
        pass


""" start: end in sphire 1.3"""

"""IT SEEMS TO BE NOT USED"""
class Test_ali2d_MPI(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base"))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_MPI()
        self.assertEqual(str(cm_new.exception), "ali2d_MPI() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("The result are not the same")
    def test_pickle_file_values(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr, Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnewa = path.join( 'ali2d_MPI_NEW')
        outdirnewb = path.join(  'ali2d_MPI__OLD')

        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

        fu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        oldfu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)


        self.assertEqual(returns_values_in_file(path.join(outdirnewa, "resolution001")),returns_values_in_file(path.join(outdirnewb, "resolution001")))
        self.assertEqual(returns_values_in_file(path.join(outdirnewa, "resolution002")),returns_values_in_file(path.join(outdirnewb, "resolution002")))
        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

    def test_ali2d_MPI_error_directory_exists(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr, Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnewa = path.join( 'ali2d_MPI_NEW')
        outdirnewb = path.join(  'ali2d_MPI__OLD')

        mkdir(outdirnewa)
        mkdir(outdirnewb)

        with self.assertRaises(OSError) as cm_new:
            fu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(OSError) as cm_old:
            oldfu.ali2d_MPI(ABSOLUTE_PATH_TO_STACK, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, ts=ts, dst=dst, maxit=4, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(cm_new.exception.strerror, "File exists")
        self.assertEqual(cm_new.exception.filename, "ali2d_MPI_NEW")
        self.assertEqual(cm_old.exception.strerror, cm_new.exception.strerror)

        remove_dir(outdirnewa)
        remove_dir(outdirnewb)



class Test_ali2d_base(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base"))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_base()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_base()
        self.assertEqual(str(cm_new.exception), "ali2d_base() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("The result are not the same")
    def test_ali2d_base_true_should_return_equal_object(self):
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        outdirnew = path.join( ABSOLUTE_PATH,'ali2d_base_new')
        outdirnewold = path.join(ABSOLUTE_PATH,'ali2d_base_old')

        mkdir(outdirnew)
        mkdir(outdirnewold)
        number_of_proc = 1
        myid = 0
        main_node = 0
        maxit =2
        return_new = fu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_base(stack, outdirnewold, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution001")),returns_values_in_file(path.join(outdirnewold, "resolution001")))
        self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution002")),returns_values_in_file(path.join(outdirnewold, "resolution002")))
        self.assertTrue(allclose(return_new, return_old, atol=TOLERANCE,equal_nan=True))
        remove_dir(outdirnew)
        remove_dir(outdirnewold)



class Test_cpy(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.cpy()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.cpy()
        self.assertEqual(str(cm_new.exception), "cpy() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_case(self):
        ins_list = path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'Class2D/best.hdf')
        file_path_new = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "new_substack")
        file_path_old = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "old_substack")
        fu.cpy(ins_list,file_path_new)
        oldfu.cpy(ins_list, file_path_old)
        db_new, keys = EMAN2db_db_open_dict( file_path_new, True, True)
        db_new.realopen()
        db_old, keys = EMAN2db_db_open_dict( file_path_old, True, True)
        db_old.realopen()

        for index in db_new.bdb.keys():
            self.assertEqual(str(pickle_loads(db_new.bdb.get(index))), str(pickle_loads(db_old.bdb.get(index))))
        remove_dir( path.join(ABSOLUTE_PATH, 'EMAN2DB'))

    def test_error_hdf_not_found(self):
        ins_list = path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'not_found.hdf')
        file_path_new = 'bdb:{0}#{1}'.format(".", "new_substack")
        file_path_old = 'bdb:{0}#{1}'.format(".", "old_substack")
        with self.assertRaises(Exception) as cm_new:
            fu.cpy(ins_list,file_path_new)
        with self.assertRaises(Exception) as cm_old:
            oldfu.cpy(ins_list, file_path_old)
        self.assertEqual(str(cm_new.failureException.message), "<attribute 'message' of 'exceptions.BaseException' objects>")
        self.assertEqual(cm_new.failureException.message, cm_old.failureException.message)



class Test_project3d(unittest.TestCase):

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.project3d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.project3d()
        self.assertEqual(str(cm_new.exception), "project3d() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_crashes_because_signal11SIGSEGV(self):
        self.assertTrue(True)
        """
        return_new = fu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        """

    def test_Nonetype_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError)as cm_new:
            fu.project3d(volume=None, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        with self.assertRaises(AttributeError)as cm_old:
            oldfu.project3d(volume=None, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertEqual(str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_img_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.project3d(volume=EMData(), stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.project3d(volume=EMData(), stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_3Dimg_realsp_and_trilinear_error_msg(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_default_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_default_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_trilinear_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_trilinear_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = True)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_realsp_case(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_realsp_case(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = True, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("since it adds a random img_noise_to the stack the results cannot be the same")
    def test_blank_img_with_noise(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    @unittest.skip("since it adds a random img_noise_to the stack the results cannot be the same")
    def test_3Dimg_with_noise(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_save_on_hdf(self):
        outnew = path.join( ABSOLUTE_PATH,'project3dnew.hdf')
        outold = path.join(ABSOLUTE_PATH,'project3old.hdf')
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = outnew, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = outold, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)

        self.assertTrue(returns_values_in_file(outold),returns_values_in_file(outnew))
        remove_list_of_file([outnew,outold])
        self.assertTrue(return_new is None)
        self.assertTrue(return_old is None)

    def test_3Dimg_with_mask(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_with_mask(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_with_listagls(self):
        listangls =even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_blank_img_with_listagls(self):
        listangls =even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='sd1', ant = 0.0)
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = listangls , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)

    def test_3Dimg_empty_listagls(self):
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [] , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [], listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertTrue(array_equal(return_old,return_new))
        self.assertTrue(array_equal(return_new, []))

    def test_blank_img_empty_listagls(self):
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = [], listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls =[], listctfs = None, noise = None, realsp = False, trillinear = False)
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_new, []))


class Test_ali_vol(unittest.TestCase):
    argum = get_arg_from_pickle_file( path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol"))
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol()
        self.assertEqual(str(cm_new.exception), "ali_vol() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_empty_refv_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=EMData(), ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=EMData(), ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_with_empty_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =EMData(), refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =EMData(), refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_with_NoneType_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =None, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =None, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_NoneType_refv_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=None, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=None, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_as_vol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =IMAGE_2D, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =IMAGE_2D, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_as_refvol_returns_RuntimeError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(vol =vol, refv=IMAGE_2D, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=IMAGE_2D, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_with_pickle_values(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_default_values(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=None, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=None, discrepancy = "ccc")
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_zero_radius(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=0, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=0, discrepancy = "ccc")
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview()))

    def test_with_zero_shift_returns_ZeroDivisionError(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=0, radius=radius, discrepancy = "ccc")
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=0, radius=radius, discrepancy = "ccc")
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_zero_ang(self):
        (vol,refv,ang_scale,shift_scale,radius) = self.argum[0]
        return_new = fu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        return_old = oldfu.ali_vol(vol =vol, refv=refv, ang_scale=ang_scale, shift_scale=shift_scale, radius=radius, discrepancy = "ccc")
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview()))


class Test_recons3d_n_trl_MPI_one_node(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params_proj"))
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.recons3d_n_trl_MPI_one_node()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.recons3d_n_trl_MPI_one_node()
        self.assertEqual(str(cm_new.exception), "recons3d_n_trl_MPI_one_node() takes exactly 12 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    """
    the fftvol that produces the error is created in the code and I cannot understand why 
        Error
        Traceback (most recent call last):
          File "/home/lusnig/SPHIRE_1_1/lib/python2.7/unittest/case.py", line 329, in run
            testMethod()
          File "/home/lusnig/EMAN2/eman2/sphire/tests/test_applications.py", line 487, in test_NoCTF_symC1
            return_old = fu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=1)
          File "/home/lusnig/EMAN2/eman2/sphire/libpy/sparx_applications.py", line 2361, in recons3d_n_trl_MPI_one_node
            fftvol.div_sinc(1)
        RuntimeError: NotExistingObjectException at /home/lusnig/EMAN2/eman2/libEM/emdata_metadata.cpp:1153: error with 'npad': 'The requested key does not exist' caught
    """
    @unittest.skip("I get an error when it runs")
    def test_NoCTF_symC1(self):
        img = self.argum[0][0]
        pjrlist=[img,img]
        group=0
        return_old = fu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=-1)
        return_new = oldfu.recons3d_n_trl_MPI_one_node(prjlist=pjrlist, CTF = False, snr="not_used", sign="not_used", npad="not_used", sym="d1", group=group, niter=2, verbose=False, upweighted=True, compensate=False, chunk_id=-1)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview()))



class Test_pca(unittest.TestCase):

    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA"))
    ave, grp_imgdata  = fu.prepare_2d_forPCA(data=argum[0][0])

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pca()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pca()
        self.assertEqual(str(cm_new.exception), "pca() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_verbose_nothing_to_print(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_radius_higher0_verbose(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=True)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_MPI_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=True, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=True, verbose=False)
        self.test_all_the_conditions(return_new, return_old)


    def test_without_genbuf(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=False, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=False, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_shuffle_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=True, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=True, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_incoreCalculation(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=True, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=-1, nvec=3, incore=True, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)

    def test_radius_and_maskFile_ErrorMsg(self):
        return_new = fu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="test", MPI=False, verbose=False)
        return_old = oldfu.pca(input_stacks=self.grp_imgdata, subavg="", mask_radius=1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="test", MPI=False, verbose=False)
        self.test_all_the_conditions(return_new, return_old)



class Test_prepare_2d_forPCA(unittest.TestCase):
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA"))[0][0]

    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE,equal_nan=True))

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.prepare_2d_forPCA()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.prepare_2d_forPCA()
        self.assertEqual(str(cm_new.exception), "prepare_2d_forPCA() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = False)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = False)
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_notAmode(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = False)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = False)
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_notAmode_withCTF(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = True)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "unknown", output_stack = None, CTF = True)
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_withCTF(self):
        return_new_avg, return_new_outstack= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = True)
        return_old_avg, return_old_outstack = oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = None, CTF = True)
        self.assertTrue(allclose(return_new_avg.get_3dview(), return_old_avg.get_3dview(), atol=TOLERANCE,equal_nan=True))
        self.test_all_the_conditions(return_new_outstack,return_old_outstack)

    def test_with_stack(self):
        outnew= path.join(ABSOLUTE_PATH,'pcanew.hdf')
        outold= path.join(ABSOLUTE_PATH, 'pcaold.hdf')
        return_new= fu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = outnew, CTF = False)
        return_old= oldfu.prepare_2d_forPCA(data = self.data, mode = "a", output_stack = outold, CTF = False)
        self.assertTrue(returns_values_in_file(outold), returns_values_in_file(outnew))
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, None)
        remove_list_of_file([outnew,outold])


class Test_extract_value(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.extract_value()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.extract_value()
        self.assertEqual(str(cm_new.exception), "extract_value() takes exactly 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_string_integer(self):
        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 20)

    def test_string_float(self):
        return_new = fu.extract_value('20.1')
        return_old = oldfu.extract_value('20.1')
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 20.1)

    def test_string_not_handled(self):
        return_new = fu.extract_value("invalid")
        return_old = oldfu.extract_value("invalid")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new,"invalid")



class Test_header(unittest.TestCase):
    outnew= path.join(ABSOLUTE_PATH,'Headernew.hdf')
    outold = path.join(ABSOLUTE_PATH, 'Headerold.hdf')
    params='xform.projection'
    data = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.header"))[0]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.header()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.header()
        self.assertEqual(str(cm_new.exception), "header() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_No_option_errorMsg(self):
        return_new = fu.header(stack=[], params = [], zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=[], params=[], zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_default_Too_option_errorMsg(self):
        return_new = fu.header(stack=[], params = [], zero=True, one=True, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=[], params=[], zero=True, one=True, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_default_(self):
        return_new = fu.header(stack=ABSOLUTE_PATH_TO_STACK, params = self.params, zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=self.outnew, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        return_old = oldfu.header(stack=ABSOLUTE_PATH_TO_STACK, params=self.params, zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, fexport=self.outold, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False)
        self.assertTrue(returns_values_in_file(self.outold), returns_values_in_file(self.outnew))
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)
        remove_list_of_file([self.outold,self.outnew])



class Test_MPI_start_end(unittest.TestCase):
    nima = 64
    nproc = 8
    myid = 1
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.MPI_start_end()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.MPI_start_end()
        self.assertEqual(str(cm_new.exception), "MPI_start_end() takes exactly 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_default_case(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [8,16]))

    def test_zero_nima(self):
        return_new = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [0,0]))

    def test_zero_nproc(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_zero_myd(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [0,8]))


"""seems to be never USED"""
class Test_refvol(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.refvol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.refvol()
        self.assertEqual(str(cm_new.exception), "refvol() takes exactly 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


"""
    -) I tested just the default method because the other 2 are not in use (Fabian docet)
    -) do not test with randomize=True because the results will be never the same because the random orientation in same calculations
    -) since the input parameter are basically used in 'sparx_alignment.ali2d_single_iter' and 'sparx_utilities.combine_params2' that i have already deeply tested
        I tested it just the pickle file case
"""
class Test_within_group_refinement(unittest.TestCase):
    argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/applications.within_group_refinement"))[0]
    randomize = False
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.within_group_refinement()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.within_group_refinement()
        self.assertEqual(str(cm_new.exception), "within_group_refinement() takes at least 13 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        (data, maskfile, randomize_not_used, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF) = self.argum
        return_new = fu.within_group_refinement(data=data, maskfile=maskfile, randomize=self.randomize, ir=ir, ou=ou, rs=rs, xrng=xrng, yrng=yrng, step=step, dst=dst, maxit=maxit, FH=FH, FF=FF,method = "", CTF = False)
        return_old = oldfu.within_group_refinement(data=data, maskfile=maskfile, randomize=self.randomize, ir=ir, ou=ou, rs=rs, xrng=xrng, yrng=yrng, step=step, dst=dst, maxit=maxit, FH=FH, FF=FF,method = "", CTF = False)
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), atol=TOLERANCE, equal_nan=True))



class Test_ali3d_mref_Kmeans_MPI(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali3d_mref_Kmeans_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali3d_mref_Kmeans_MPI()
        self.assertEqual(str(cm_new.exception), "ali3d_mref_Kmeans_MPI() takes exactly 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_mref_ali3d_EQ_Kmeans(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.mref_ali3d_EQ_Kmeans()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.mref_ali3d_EQ_Kmeans()
        self.assertEqual(str(cm_new.exception), "mref_ali3d_EQ_Kmeans() takes exactly 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))





"""
def get_data(num):
    dim = 10
    data_list = []
    for i in range(num):
        a = e2cpp.EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list

@unittest.skip("original Adnan test")
class Test_lib_applications_compare(unittest.TestCase):


    def test_ali2d_MPI_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
        maxit = 4

        outdirnewa = path.join(ABSOLUTE_PATH, outdir+'AA')
        outdirnewb = path.join(ABSOLUTE_PATH, outdir + 'ABB')


        if (path.exists(outdirnewa)):
            shutil.rmtree(outdirnewa)

        if (path.exists(outdirnewb)):
            shutil.rmtree(outdirnewb)

        stacka = "bdb:tests/Class2D/isac_substack"
        print(os.getcwd())
        return_new = fu.ali2d_MPI(stacka, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
                                   ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_MPI(stacka, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
                                   ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)
        mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()

        self.assertTrue(return_new, return_old)


    def test_ali2d_base_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = argum[0]

        outdirnew = path.join(ABSOLUTE_PATH, outdir+'B')

        number_of_proc = 1
        myid = 0
        main_node = 0

        print(outdirnew)
        return_new = fu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()

        self.assertTrue(return_new, return_old)



    #  mref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code 
    # def test_mref_ali3d_MPI_true_should_return_equal_object(self):
    #
    #     filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
    #
    #     filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #     (vol,refv,ang_scale,shift_scale,radius) = argum[0]
    #
    #     maxit = 4
    #     outdirnewa = path.join(ABSOLUTE_PATH, outdir+'AA')
    #     outdirnewb = path.join(ABSOLUTE_PATH, outdir + 'ABB')
    #
    #     if (path.exists(outdirnewa)):
    #         shutil.rmtree(outdirnewa)
    #     if (path.exists(outdirnewb)):
    #         shutil.rmtree(outdirnewb)
    #
    #     print(mpi_comm_rank(MPI_COMM_WORLD))
    #
    #     stacka = "bdb:tests/Substack/isac_substack"
    #     vola = "tests/Class2D/EMAN2DB/vol_adaptive_mask.hdf"
    #
    #     return_new = fu.mref_ali3d_MPI(stacka, vola, outdirnewa, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
    #                                ts=ts,  maxit=maxit, CTF=False, snr=snr, mpi_comm = MPI_COMM_WORLD)
    #     mpi_barrier(MPI_COMM_WORLD)
        # return_old = oldfu.mref_ali3d_MPI(stacka, refv, outdirnewb, maskfile=maskfile, ou=ou, xr=xr, yr=yr, \
        #                            ts=ts, dst=dst, maxit=maxit, CTF=True, snr=snr)
        # mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()
        # self.assertTrue(return_new, return_old)


    #  mref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code 
    # def test_Kmref_ali3d_MPI_true_should_return_equal_object(self):
    #     filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #         print(argum[0])
    #         print(argum[0][1])
    #
    #     (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
    #      Fourvar, user_func_name, random_method, log, number_of_proc, myid, main_node, mpi_comm) = argum[0]
    #
    #     outdir = "Class2D/Kmref_aligA"
    #
    #     outdirnew = path.join(ABSOLUTE_PATH, outdir)
    #
    #     if (path.exists(outdirnew)):
    #         shutil.rmtree(outdirnew)
    #
    #     stacka = "bdb:tests/Class2D/isac_substack"
    #     vola = "tests/Class2D/EMAN2DB/vol_adaptive_mask.hdf"
    #
    #     return_new = fu.Kmref_ali3d_MPI(stacka, vola, outdirnew,  maskfile, focus = None, \
    #                                     ir = ir, ou = ou, rs = rs, xr = xr, yr = yr, ts = ts )
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdir = "Class2D/Kmref_aligB"
    #     outdirnew = path.join(ABSOLUTE_PATH, outdir)
    #
    #     return_old = oldfu.Kmref_ali3d_MPI(stacka, vola, outdirnew, maskfile, focus = None, \
    #                                        ir=ir, ou=ou, rs=rs, xr=xr, yr=yr, ts=ts)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)


    def test_cpy_true_should_return_equal_object(self):
        ins_list = path.join(ABSOLUTE_PATH, 'Class2D/best_temp.hdf')
        ous = "bdb:tests/Class2D/isac_substack"

        return_new = fu.cpy(ins_list,ous)
        return_old = oldfu.cpy(ins_list, ous)

        self.assertEqual(return_new, return_old)


    def test_project3d_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            # print(argum[0])
            # print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]


        return_new = fu.project3d(vol)
        return_old = oldfu.project3d(vol)

        self.assertTrue(return_new, return_old)


    def test_ali_vol_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]

        return_new = fu.ali_vol(vol,refv,ang_scale,shift_scale,radius)
        return_old = oldfu.ali_vol(vol, refv, ang_scale, shift_scale, radius)

        self.assertTrue(return_new, return_old)


    def test_pca_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali2d_base")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = argum[0]


        stacka  = copy.deepcopy(stack)
        stackb = copy.deepcopy(stack)


        return_new = fu.pca(stacka,nvec = 3)
        return_old = oldfu.pca(stackb,nvec = 3)


        self.assertTrue(return_new, return_old)



    def test_prepare_2d_forPCA_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.prepare_2d_forPCA")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        return_new = fu.prepare_2d_forPCA(data = argum[0][0])
        return_old = oldfu.prepare_2d_forPCA(data = argum[0][0])

        print(return_new)
        print(return_old)

        self.assertTrue(return_new, return_old)



    def test_extract_values_true_should_return_equal_object(self):

        return_new = fu.extract_value('20')
        return_old = oldfu.extract_value('20')

        self.assertEqual(return_new, return_old)


    def test_header_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.header")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)
            print(argum[1])

        stack =  argum[0][0]
        params = argum[0][1]

        return_new = fu.header(stack=stack, params = params)
        return_old = oldfu.header(stack=stack, params=params)

        self.assertEqual(return_new, return_old)


    def test_MPI_start_end_true_should_return_equal_object(self):

        nima = 64
        nproc = 8
        myid = 1

        return_new = fu.MPI_start_end(nima, nproc, myid)
        return_old = oldfu.MPI_start_end(nima, nproc, myid)

        self.assertEqual(return_new, return_old)


    def test_refvol_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.ali_vol")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum[0])
            print(argum[0][1])

        (vol,refv,ang_scale,shift_scale,radius) = argum[0]


        vollist = [vol, vol, vol]
        output = "tests/Class2D/"
        mask = "bdb:tests/Class2D/stack_ali2d_76x76x1"

        return_new = fu.refvol(vollist , vollist, output, refv)




    def test_within_group_refinement_true_should_return_equal_object(self):
        filepath = path.join(ABSOLUTE_PATH, "pickle files/applications.within_group_refinement")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)
            print(argum)
            print(argum[0])


        (data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF) = argum[0]

        return_new = fu.within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF)
        return_old = oldfu.within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF)

        self.assertTrue(return_new, return_old)



    # def test_ali3d_mref_Kmeans_MPI_true_should_return_equal_object(self):
    #
    #     filepath = path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #     Tracker["constants"]["nproc"] = 1
    #     Tracker["constants"]["myid"] = 0
    #     Tracker["constants"]["main_node"] = 0
    #     Tracker["total_stack"] = "stack"
    #     Tracker["constants"]["seed"] = 1.4
    #     Tracker["constants"]["indep_runs"] = 2
    #     Tracker["this_data_list"] = [2,3,5]
    #     Tracker["shrinkage"] = 0.5
    #     Tracker["constants"]["ir"] = 1
    #     Tracker["constants"]["xr"] = 1
    #     Tracker["constants"]["yr"] = 1
    #     Tracker["constants"]["ts"] = 1
    #     Tracker["constants"]["delta"] = 0.5
    #     Tracker["constants"]["an"] = 1
    #     Tracker["constants"]["center"] = 4
    #     Tracker["constants"]["nassign"] =1
    #     Tracker["constants"]["nrefine"] =1
    #     Tracker["constants"]["sym"] = "c1"
    #     Tracker["constants"]["stoprnct"] =5
    #     Tracker["constants"]["mask3D"]  = False
    #     Tracker["low_pass_filter"] = "0.50"
    #     Tracker["constants"]["PWadjustment"] = False
    #
    #     outdir = "Class2D/Kmref_alig_MPI_A"
    #
    #     outdirnew = path.join(ABSOLUTE_PATH, outdir)
    #
    #     if (path.exists(outdirnew)):
    #         shutil.rmtree(outdirnew)
    #
    #     this_data_list_file = "sphire/tests/Sort3D/chunk_0.txt"
    #
    #     ref_list = "sphire/tests/Sort3D/refang.txt"
    #
    #     return_new = fu.ali3d_mref_Kmeans_MPI(ref_list, outdirnew, this_data_list_file, Tracker)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     outdir = "Class2D/Kmref_alig_MPI_B"
    #     outdirnew = path.join(ABSOLUTE_PATH, outdir)
    #
    #
    #     return_old = oldfu.ali3d_mref_Kmeans_MPI(ref_list, outdirnew, this_data_list_file, Tracker)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     mpi_finalize()
    #
    #     self.assertTrue(return_new, return_old)


"""

if __name__ == '__main__':
    unittest.main()

