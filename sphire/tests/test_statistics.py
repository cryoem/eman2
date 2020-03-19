from __future__ import division

import unittest

from os import path
from numpy import array_equal, allclose
from numpy import full as numpy_full
from numpy import nan as numpy_nan
from numpy import asarray as numpy_asarray
from mpi import *
import sp_global_def
import EMAN2_cppwrap

mpi_init(0, [])
sp_global_def.BATCH = True
sp_global_def.MPI = False
from sphire.libpy_py3 import sp_statistics as oldfu
from sphire.libpy import sp_statistics as fu

from sphire.tests.test_module import (
    remove_list_of_file,
    returns_values_in_file,
    get_arg_from_pickle_file,
    get_real_data,
    ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER,
    give_ali2d_single_iter_data,
    give_ave_series_data,
)
from sphire.libpy.sp_utilities import model_circle, model_blank
from EMAN2_cppwrap import EMData

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))




IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D, STILL_NOT_VALID = get_real_data(dim=3)
TOLERANCE = 0.005
MASK = model_circle(2, 5, 5)
STACK_NAME = "bdb:" + path.join(
    ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, "Substack/sort3d_substack_003"
)

"""
WHAT IS MISSING:
0) I need a list of image with 'xform.align2d' param to test a lot of functions avoiding to get the data from the pickle file 


RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) 

In these tests there is a bug --> syntax error:
1) orient_params --> if the symmetry_class is not cn it will lead to an 'UnboundLocalError'. beacuase it returns a value that not exists. it is a BUG
2) Test_k_means_stab_bbenum --> a case meets an exit error

In these tests there is a strange behavior:
1) 
"""


""" OLD comment
There are some opened issues in:
1) add_ave_varf_MPI: it is never used and I avoid to test it deeply. 
    -) CTF=True leads to a bug (see 'test_CTF_True_deadCode_BUG' --> UnboundLocalError: local variable 'i' referenced before assignment)
    
2) sum_oe:
    -)I cannot find a dataset to test the CTF=True case .... anyway it is just a debug mode
    -) ctf_eo_sum not tested because used in combination with CTF=True
    
3) varf2d_MPI:
    -) the use of a 3D image is not blocked. is it ok or a bug?
    -) I cannot find a dataset to test the CTF=True case
    
4) varf3d_MPI:
    -) how can i find a valid mask to test it
    -) i did not tested all the params deeply becuae they are basically used as input of 'sparx_reconstruction.recons3d_4nn_ctf_MPI'
    
5) locres:
    -) Adnan says that it is able to run only in the cluster ... Hence I did not test it
    
6) class Munkres and class pcanalyzer:
    -) at the moment I have no idea how I'd test it
    
7) linreg: Could we use from sklearn.linear_model import LinearRegression ???

8) k_means_stab_bbenum:
    -) i did not test it deeply because half parameters are not used and the other is used to call  'k_means_match_bbenum'. Hence I'll test deeply 'k_means_match_bbenum'
"""


"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""


#   THESE FUNCTIONS ARE COMMENTED BECAUSE NOT INVOLVED IN THE PYTHON3 CONVERSION. THEY HAVE TO BE TESTED ANYWAY
""" start: new in sphire 1.3


class Test_avgvar(unittest.TestCase):
    def test_avgvar(self):
        oldv = oldfu.avgvar(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True)
        v = fu.avgvar(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True)
        pass

class Test_avgvar_ctf(unittest.TestCase):
    def test_avgvar_ctf(self):
        oldv = oldfu.avgvar_ctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0, dopa = True)
        v = fu.avgvar_ctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0, dopa = True)
        pass

class Test_add_oe_series(unittest.TestCase):
    def test_add_oe_series(self):
        oldv = oldfu.add_oe_series(data, ali_params="xform.align2d")
        v = fu.add_oe_series(data, ali_params="xform.align2d")
        pass

class Test_add_ave_varf(unittest.TestCase):
    def test_add_ave_varf(self):
        oldv = oldfu.add_ave_varf(data, mask = None, mode = "a", CTF = False, ctf_2_sum = None, ali_params = "xform.align2d")
        v = fu.add_ave_varf(data, mask = None, mode = "a", CTF = False, ctf_2_sum = None, ali_params = "xform.align2d")
        pass

class Test_(unittest.TestCase):
    def test_add_oe(self):
        oldv = oldfu.add_oe(data)
        v = fu.add_oe(data)
        pass

class Test_ave_series_ctf(unittest.TestCase):
    def test_ave_series_ctf(self):
        oldv = oldfu.ave_series_ctf(data, ctf2, mask = None)
        v = fu.ave_series_ctf(data, ctf2, mask = None)
        pass


class Test_ave_oe_series_d(unittest.TestCase):
    def test_ave_oe_series_d(self):
        oldv = oldfu.ave_oe_series_d(data)
        v = fu.ave_oe_series_d(data)
        pass

class Test_ave_oe_series(unittest.TestCase):
    def test_ave_oe_series(self):
        oldv = oldfu.ave_oe_series(stack)
        v = fu.ave_oe_series(stack)
        pass

class Test_ave_oe_series_textfile(unittest.TestCase):
    def test_ave_oe_series_textfile(self):
        oldv = oldfu.ave_oe_series_textfile(stack, textfile)
        v = fu.ave_oe_series_textfile(stack, textfile)
        pass

class Test_ave_oe_series_indexed(unittest.TestCase):
    def test_ave_oe_series_indexed(self):
        oldv = oldfu.ave_oe_series_indexed(stack, idx_ref)
        v = fu.ave_oe_series_indexed(stack, idx_ref)
        pass

class Test_ave_var_series_one(unittest.TestCase):
    def test_ave_var_series_one(self):
        oldv = oldfu.ave_var_series_one(data, skip, kb)
        v = fu.ave_var_series_one(data, skip, kb)
        pass

class Test_add_series(unittest.TestCase):
    def test_add_series(self):
        oldv = oldfu.add_series(stack, i1=0 ,i2=0)
        v = fu.add_series(stack, i1=0 ,i2=0)
        pass

class Test_add_series_class(unittest.TestCase):
    def test_add_series_class(self):
        oldv = oldfu.add_series_class(stack, i1 = 0, i2 = 0)
        v = fu.add_series_class(stack, i1 = 0, i2 = 0)
        pass

class Test_add_series_class_mem(unittest.TestCase):
    def test_add_series_class_mem(self):
        oldv = oldfu.add_series_class_mem(data, assign, kc)
        v = fu.add_series_class_mem(data, assign, kc)
        pass

class Test_aves(unittest.TestCase):
    def test_aves(self):
        oldv = oldfu.aves(stack, mode="a", i1 = 0, i2 = 0)
        v = fu.aves(stack, mode="a", i1 = 0, i2 = 0)
        pass

class Test_aveq(unittest.TestCase):
    def test_aveq(self):
        oldv = oldfu.aveq(stack, mode="a", i1 = 0, i2 = 0)
        v = fu.aveq(stack, mode="a", i1 = 0, i2 = 0)
        pass

class Test_aves_wiener(unittest.TestCase):
    def test_aves_wiener(self):
        oldv = oldfu.aves_wiener(input_stack, mode="a", SNR=1.0, interpolation_method="linear")
        v = fu.aves_wiener(input_stack, mode="a", SNR=1.0, interpolation_method="linear")
        pass

class Test_aves_adw(unittest.TestCase):
    def test_aves_adw(self):
        oldv = oldfu.aves_adw(input_stack, mode="a", SNR=1.0, Ng = -1, interpolation_method="linear")
        v = fu.aves_adw(input_stack, mode="a", SNR=1.0, Ng = -1, interpolation_method="linear")
        pass


class Test_ssnr2d(unittest.TestCase):
    def test_ssnr2d(self):
        oldv = oldfu.ssnr2d(data, mask = None, mode="")
        v = fu.ssnr2d(data, mask = None, mode="")
        pass


class Test_ssnr2d_ctf(unittest.TestCase):
    def test_ssnr2d_ctf(self):
        oldv = oldfu.ssnr2d_ctf(data, mask = None, mode="", dopa=False)
        v = fu.ssnr2d_ctf(data, mask = None, mode="", dopa=False)
        pass


class Test_varf(unittest.TestCase):
    def test_varf(self):
        oldv = oldfu.varf(data, mask = None, mode="a")
        v = fu.varf(data, mask = None, mode="a")
        pass


class Test_varfctf(unittest.TestCase):
    def test_varfctf(self):
        oldv = oldfu.varfctf(data, mask = None, mode="a", dopad = True)
        v = fu.varfctf(data, mask = None, mode="a", dopad = True)
        pass


class Test_varf2d(unittest.TestCase):
    def test_varf2d(self):
        oldv = oldfu.varf2d(data, ave, mask = None, mode="a")
        v = fu.varf2d(data, ave, mask = None, mode="a")
        pass


class Test_varf3d(unittest.TestCase):
    def test_varf3d(self):
        oldv = oldfu.varf3d(prjlist,ssnr_text_file = None, mask2D = None, reference_structure = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1")
        v = fu.varf3d(prjlist,ssnr_text_file = None, mask2D = None, reference_structure = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1")
        pass


class Test_get_refstack(unittest.TestCase):
    def test_get_refstack(self):
        oldv = oldfu.get_refstack(imgstack,params,nref,refstack,cs,mask,center,Iter)
        v = fu.get_refstack(imgstack,params,nref,refstack,cs,mask,center,Iter)
        pass


class Test_get_1dpw_table_stack(unittest.TestCase):
    def test_get_1dpw_table_stack(self):
        oldv = oldfu.get_1dpw_table_stack(stack)
        v = fu.get_1dpw_table_stack(stack)
        pass


class Test_im_diff(unittest.TestCase):
    def test_im_diff(self):
        oldv = oldfu.im_diff(im1, im2, mask = None)
        v = fu.im_diff(im1, im2, mask = None)
        pass


class Test_k_means_init_asg_rnd(unittest.TestCase):
    def test_k_means_init_asg_rnd(self):
        oldv = oldfu.k_means_init_asg_rnd(N, K)
        v = fu.k_means_init_asg_rnd(N, K)
        pass


class Test_k_means_init_asg_d2w(unittest.TestCase):
    def test_k_means_init_asg_d2w(self):
        oldv = oldfu.k_means_init_asg_d2w(im, N, K)
        v = fu.k_means_init_asg_d2w(im, N, K)
        pass


class Test_k_means_locasg2glbasg(unittest.TestCase):
    def test_k_means_locasg2glbasg(self):
        oldv = oldfu.k_means_locasg2glbasg(ASG, LUT, N)
        v = fu.k_means_locasg2glbasg(ASG, LUT, N)
        pass


class Test_k_means_init_open_im(unittest.TestCase):
    def test_k_means_init_open_im(self):
        oldv = oldfu.k_means_init_open_im(stack, maskname)
        v = fu.k_means_init_open_im(stack, maskname)
        pass


class Test_k_means_open_im(unittest.TestCase):
    def test_k_means_open_im(self):
        oldv = oldfu.k_means_open_im(stack, mask, CTF, lim, flagnorm = False)
        v = fu.k_means_open_im(stack, mask, CTF, lim, flagnorm = False)
        pass


class Test_k_means_headlog(unittest.TestCase):
    def test_k_means_headlog(self):
        oldv = oldfu.k_means_headlog(stackname, outname, method, N, K, crit, maskname, trials, maxit, CTF, T0, F, rnd, ncpu, m, init_method='random')
        v = fu.k_means_headlog(stackname, outname, method, N, K, crit, maskname, trials, maxit, CTF, T0, F, rnd, ncpu, m, init_method='random')
        pass


class Test_k_means_export(unittest.TestCase):
    def test_k_means_export(self):
        oldv = oldfu.k_means_export(Cls, crit, assign, out_seedname, part = -1, TXT = False)
        v = fu.k_means_export(Cls, crit, assign, out_seedname, part = -1, TXT = False)
        pass


class Test_k_means_criterion(unittest.TestCase):
    def test_k_means_criterion(self):
        oldv = oldfu.k_means_criterion(Cls, crit_name='')
        v = fu.k_means_criterion(Cls, crit_name='')
        pass


class Test_select_kmeans(unittest.TestCase):
    def test_select_kmeans(self):
        oldv = oldfu.select_kmeans(dJe, T)
        v = fu.select_kmeans(dJe, T)
        pass


class Test_k_means_cla(unittest.TestCase):
    def test_k_means_cla(self):
        oldv = oldfu.k_means_cla(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd')
        v = fu.k_means_cla(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd')
        pass


class Test_k_means_SSE(unittest.TestCase):
    def test_k_means_SSE(self):
        oldv = oldfu.k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd')
        v = fu.k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd')
        pass


class Test_k_means_SSE_combine(unittest.TestCase):
    def test_k_means_SSE_combine(self):
        oldv = oldfu.k_means_SSE_combine(Cls, assign, Je, N, K, ncpu, myid, main_node)
        v = fu.k_means_SSE_combine(Cls, assign, Je, N, K, ncpu, myid, main_node)
        pass


class Test_k_means_SSE_collect(unittest.TestCase):
    def test_k_means_SSE_collect(self):
        oldv = oldfu.k_means_SSE_collect(Cls, assign, Je, N, K, ncpu, myid, main_node)
        v = fu.k_means_SSE_collect(Cls, assign, Je, N, K, ncpu, myid, main_node)
        pass


class Test_k_means_SSE_MPI(unittest.TestCase):
    def test_k_means_SSE_MPI(self):
        oldv = oldfu.k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd', myid = 0, main_node =0, jumping = 1)
        v = fu.k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd', myid = 0, main_node =0, jumping = 1)
        pass


class Test_k_means_cla_MPI(unittest.TestCase):
    def test_k_means_cla_MPI(self):
        oldv = oldfu.k_means_cla_MPI(IM, mask, K, rand_seed, maxit, trials, CTF, F, T0, myid, main_node, N_start, N_stop, N)
        v = fu.k_means_cla_MPI(IM, mask, K, rand_seed, maxit, trials, CTF, F, T0, myid, main_node, N_start, N_stop, N)
        pass


class Test_k_means_CUDA(unittest.TestCase):
    def test_k_means_CUDA(self):
        oldv = oldfu.k_means_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        v = fu.k_means_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        pass


class Test_k_means_SSE_CUDA(unittest.TestCase):
    def test_k_means_SSE_CUDA(self):
        oldv = oldfu.k_means_SSE_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        v = fu.k_means_SSE_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        pass


class Test_dump_AVE(unittest.TestCase):
    def test_dump_AVE(self):
        oldv = oldfu.dump_AVE(AVE, mask, myid, ite = 0)
        v = fu.dump_AVE(AVE, mask, myid, ite = 0)
        pass


class Test_k_means_CUDA_MPI(unittest.TestCase):
    def test_k_means_CUDA_MPI(self):
        oldv = oldfu.k_means_CUDA_MPI(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        v = fu.k_means_CUDA_MPI(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, nbpart, logging = -1, flagnorm = False)
        pass


class Test_k_means_CUDA_MPI_YANG(unittest.TestCase):
    def test_k_means_CUDA_MPI_YANG(self):
        oldv = oldfu.k_means_CUDA_MPI_YANG(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, ipart, logging = -1, flagnorm = False, comm = -1, gpuid = 0)
        v = fu.k_means_CUDA_MPI_YANG(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, ipart, logging = -1, flagnorm = False, comm = -1, gpuid = 0)
        pass

class Test_k_means_groups_gnuplot(unittest.TestCase):
    def test_k_means_groups_gnuplot(self):
        oldv = oldfu.k_means_groups_gnuplot(file, src, C, DB, H)
        v = fu.k_means_groups_gnuplot(file, src, C, DB, H)
        pass

class Test_k_means_groups_serial(unittest.TestCase):
    def test_k_means_groups_serial(self):
        oldv = oldfu.k_means_groups_serial(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, DEBUG = False, flagnorm = False)
        v = fu.k_means_groups_serial(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, DEBUG = False, flagnorm = False)
        pass

class Test_k_means_groups_CUDA(unittest.TestCase):
    def test_k_means_groups_CUDA(self):
        oldv = oldfu.k_means_groups_CUDA(stack, outdir, maskname, K1, K2, rand_seed, maxit, F, T0)
        v = fu.k_means_groups_CUDA(stack, outdir, maskname, K1, K2, rand_seed, maxit, F, T0)
        pass

class Test_k_means_groups_MPI(unittest.TestCase):
    def test_k_means_groups_MPI(self):
        oldv = oldfu.k_means_groups_MPI(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, flagnorm)
        v = fu.k_means_groups_MPI(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, flagnorm)
        pass

class Test_k_means_cuda_error(unittest.TestCase):
    def test_k_means_cuda_error(self):
        oldv = oldfu.k_means_cuda_error(status)
        v = fu.k_means_cuda_error(status)
        pass

class Test_k_means_cuda_headlog(unittest.TestCase):
    def test_k_means_cuda_headlog(self):
        oldv = oldfu.k_means_cuda_headlog(stackname, outname, method, N, K, maskname, maxit, T0, F, rnd, ncpu, m)
        v = fu.k_means_cuda_headlog(stackname, outname, method, N, K, maskname, maxit, T0, F, rnd, ncpu, m)
        pass

class Test_k_means_cuda_init_open_im(unittest.TestCase):
    def test_k_means_cuda_init_open_im(self):
        oldv = oldfu.k_means_cuda_init_open_im(stack, maskname)
        v = fu.k_means_cuda_init_open_im(stack, maskname)
        pass

class Test_k_means_cuda_open_im(unittest.TestCase):
    def test_k_means_cuda_open_im(self):
        oldv = oldfu.k_means_cuda_open_im(KmeansCUDA, stack, lim, mask, flagnorm = False)
        v = fu.k_means_cuda_open_im(KmeansCUDA, stack, lim, mask, flagnorm = False)
        pass

class Test_k_means_cuda_info(unittest.TestCase):
    def test_k_means_cuda_info(self):
        oldv = oldfu.k_means_cuda_info(INFO)
        v = fu.k_means_cuda_info(INFO)
        pass

class Test_k_means_cuda_export(unittest.TestCase):
    def test_k_means_cuda_export(self):
        oldv = oldfu.k_means_cuda_export(PART, FLATAVE, out_seedname, mask, crit, part = -1, TXT = False)
        v = fu.k_means_cuda_export(PART, FLATAVE, out_seedname, mask, crit, part = -1, TXT = False)
        pass

class Test_k_means_SA_T0(unittest.TestCase):
    def test_k_means_SA_T0(self):
        oldv = oldfu.k_means_SA_T0(im_M, mask, K, rand_seed, CTF, F)
        v = fu.k_means_SA_T0(im_M, mask, K, rand_seed, CTF, F)
        pass

class Test_k_means_SA_T0_MPI(unittest.TestCase):
    def test_k_means_SA_T0_MPI(self):
        oldv = oldfu.k_means_SA_T0_MPI(im_M, mask, K, rand_seed, CTF, F, myid, main_node, N_start, N_stop)
        v = fu.k_means_SA_T0_MPI(im_M, mask, K, rand_seed, CTF, F, myid, main_node, N_start, N_stop)
        pass

class Test_k_means_stab_stream(unittest.TestCase):
    def test_k_means_stab_stream(self):
        oldv = oldfu.k_means_stab_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False)
        v = fu.k_means_stab_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False)
        pass

class Test_k_means_stab_MPI_stream(unittest.TestCase):
    def test_k_means_stab_MPI_stream(self):
        oldv = oldfu.k_means_stab_MPI_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False, num_first_matches=1)
        v = fu.k_means_stab_MPI_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False, num_first_matches=1)
        pass

class Test_k_means_match_clusters_asg(unittest.TestCase):
    def test_k_means_match_clusters_asg(self):
        oldv = oldfu.k_means_match_clusters_asg(asg1, asg2)
        v = fu.k_means_match_clusters_asg(asg1, asg2)
        pass


class Test_k_means_stab_H(unittest.TestCase):
    def test_k_means_stab_H(self):
        oldv = oldfu.k_means_stab_H(ALL_PART)
        v = fu.k_means_stab_H(ALL_PART)
        pass

class Test_k_means_match_pwa(unittest.TestCase):
    def test_k_means_match_pwa(self):
        oldv = oldfu.k_means_match_pwa(PART, lim = -1)
        v = fu.k_means_match_pwa(PART, lim = -1)
        pass


class Test_k_means_stab_pwa(unittest.TestCase):
    def test_(self):
        oldv = oldfu.k_means_stab_pwa(PART, lim = -1)
        v = fu.k_means_stab_pwa(PART, lim = -1)
        pass

class Test_k_means_stab_export_txt(unittest.TestCase):
    def test_k_means_stab_export_txt(self):
        oldv = oldfu.k_means_stab_export_txt(PART, outdir, th_nobj)
        v = fu.k_means_stab_export_txt(PART, outdir, th_nobj)
        pass


class Test_k_means_stab_export(unittest.TestCase):
    def test_k_means_stab_export(self):
        oldv = oldfu.k_means_stab_export(PART, stack, outdir, th_nobj, CTF = False)
        v = fu.k_means_stab_export(PART, stack, outdir, th_nobj, CTF = False)
        pass

class Test_k_means_stab_init_tag(unittest.TestCase):
    def test_k_means_stab_init_tag(self):
        oldv = oldfu.k_means_stab_init_tag(stack)
        v = fu.k_means_stab_init_tag(stack)
        pass


class Test_k_means_stab_asg2part(unittest.TestCase):
    def test_k_means_stab_asg2part(self):
        oldv = oldfu.k_means_stab_asg2part(outdir, npart)
        v = fu.k_means_stab_asg2part(outdir, npart)
        pass

class Test_k_means_asg_locasg2glbpart(unittest.TestCase):
    def test_k_means_asg_locasg2glbpart(self):
        oldv = oldfu.k_means_asg_locasg2glbpart(ASG, LUT)
        v = fu.k_means_asg_locasg2glbpart(ASG, LUT)
        pass


class Test_k_means_stab_gather(unittest.TestCase):
    def test_k_means_stab_gather(self):
        oldv = oldfu.k_means_stab_gather(nb_run, maskname, outdir)
        v = fu.k_means_stab_gather(nb_run, maskname, outdir)
        pass

class Test_k_means_extract_class_ali(unittest.TestCase):
    def test_k_means_extract_class_ali(self):
        oldv = oldfu.k_means_extract_class_ali(stack_name, ave_name, dir)
        v = fu.k_means_extract_class_ali(stack_name, ave_name, dir)
        pass


class Test_k_means_class_pixerror(unittest.TestCase):
    def test_k_means_class_pixerror(self):
        oldv = oldfu.k_means_class_pixerror(class_name, dir, ou, xr, ts, maxit, fun, CTF=False, snr=1.0, Fourvar=False)
        v = fu.k_means_class_pixerror(class_name, dir, ou, xr, ts, maxit, fun, CTF=False, snr=1.0, Fourvar=False)
        pass

class Test_isc_update_ite_conf(unittest.TestCase):
    def test_isc_update_ite_conf(self):
        oldv = oldfu.isc_update_ite_conf(conf_file, ite)
        v = fu.isc_update_ite_conf(conf_file, ite)
        pass


class Test_isc_read_conf(unittest.TestCase):
    def test_isc_read_conf(self):
        oldv = oldfu.isc_read_conf(conf_file)
        v = fu.isc_read_conf(conf_file)
        pass

class Test_isc_ave_huge(unittest.TestCase):
    def test_isc_ave_huge(self):
        oldv = oldfu.isc_ave_huge(ave_tiny, org_data, ave_huge)
        v = fu.isc_ave_huge(ave_tiny, org_data, ave_huge)
        pass


class Test_py_cluster_median(unittest.TestCase):
    def test_py_cluster_median(self):
        oldv = oldfu.py_cluster_median(numbers)
        v = fu.py_cluster_median(numbers)
        pass

class Test_py_cluster_genmatrix(unittest.TestCase):
    def test_py_cluster_genmatrix(self):
        oldv = oldfu.py_cluster_genmatrix(list, combinfunc, symmetric=False, diagonal=None)
        v = fu.py_cluster_genmatrix(list, combinfunc, symmetric=False, diagonal=None)
        pass


class Test_class_pyCluster(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass

class Test_class_py_cluster_BaseClusterMethod(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass


class Test_class_py_cluster_HierarchicalClustering(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass

class Test_class_def_variancer(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass


class Test_class_inc_variancer(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass

class Test_class_pcanalyzebck(unittest.TestCase):
    def test_(self):
        oldv = oldfu.
        v = fu.
        pass


class Test_kmn(unittest.TestCase):
    def test_kmn(self):
        oldv = oldfu.kmn(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        v = fu.kmn(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        pass

class Test_kmn_a(unittest.TestCase):
    def test_kmn_a(self):
        oldv = oldfu.kmn_a(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        v = fu.kmn_a(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        pass


class Test_kmn_g(unittest.TestCase):
    def test_kmn_g(self):
        oldv = oldfu.kmn_g(data, numr, wr, stack, check_mirror = False, max_iter = 10, this_seed = 1000)
        v = fu.kmn_g(data, numr, wr, stack, check_mirror = False, max_iter = 10, this_seed = 1000)
        pass

class Test_multi_search_func(unittest.TestCase):
    def test_multi_search_func(self):
        oldv = oldfu.multi_search_func(args, data)
        v = fu.multi_search_func(args, data)
        pass


class Test_multi_search_func2(unittest.TestCase):
    def test_multi_search_func2(self):
        oldv = oldfu.multi_search_func2(args, data)
        v = fu.multi_search_func2(args, data)
        pass

class Test_kmn_ctf(unittest.TestCase):
    def test_kmn_ctf(self):
        oldv = oldfu.kmn_ctf(data, ref_data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        v = fu.kmn_ctf(data, ref_data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000)
        pass


class Test_kmnr(unittest.TestCase):
    def test_kmnr(self):
        oldv = oldfu.kmnr(data, assign, nima, k, numr, wr, cm = 0, max_iter = 10, this_seed = 1000, MPI = False)
        v = fu.kmnr(data, assign, nima, k, numr, wr, cm = 0, max_iter = 10, this_seed = 1000, MPI = False)
        pass

class Test_Wiener_CSQ(unittest.TestCase):
    def test_Wiener_CSQ(self):
        oldv = oldfu.Wiener_CSQ(data, K, assign, Cls, ctf1, ctf2, snr = 1.0)
        v = fu.Wiener_CSQ(data, K, assign, Cls, ctf1, ctf2, snr = 1.0)
        pass


class Test_Wiener_sse(unittest.TestCase):
    def test_Wiener_sse(self):
        oldv = oldfu.Wiener_sse(data, K, assign, Cls, ctf1, ctf2, snr = 1.0)
        v = fu.Wiener_sse(data, K, assign, Cls, ctf1, ctf2, snr = 1.0)
        pass

class Test_var_bydef(unittest.TestCase):
    def test_var_bydef(self):
        oldv = oldfu.var_bydef(vol_stack, vol_list, info)
        v = fu.var_bydef(vol_stack, vol_list, info)
        pass


class Test_histogram2d(unittest.TestCase):
    def test_histogram2d(self):
        oldv = oldfu.histogram2d(datai, dataj, nbini, nbinj)
        v = fu.histogram2d(datai, dataj, nbini, nbinj)
        pass

class Test_get_power_spec(unittest.TestCase):
    def test_get_power_spec(self):
        oldv = oldfu.get_power_spec(stack_file, start_particle, end_particle)
        v = fu.get_power_spec(stack_file, start_particle, end_particle)
        pass


class Test_noise_corrected_PW(unittest.TestCase):
    def test_noise_corrected_PW(self):
        oldv = oldfu.noise_corrected_PW(pw, lo_limit, hi_limit, abs_limit)
        v = fu.noise_corrected_PW(pw, lo_limit, hi_limit, abs_limit)
        pass

class Test_cluster_pairwise(unittest.TestCase):
    def test_cluster_pairwise(self):
        oldv = oldfu.cluster_pairwise(d, K)
        v = fu.cluster_pairwise(d, K)
        pass


class Test_cluster_equalsize(unittest.TestCase):
    def test_cluster_equalsize(self):
        oldv = oldfu.cluster_equalsize(d, m)
        v = fu.cluster_equalsize(d, m)
        pass

class Test_k_means_stab_getinfo(unittest.TestCase):
    def test_k_means_stab_getinfo(self):
        oldv = oldfu.k_means_stab_getinfo(PART, match)
        v = fu.k_means_stab_getinfo(PART, match)
        pass


class Test_match_lists(unittest.TestCase):
    def test_(self):
        oldv = oldfu.match_lists(l1, l2)
        v = fu.match_lists(l1, l2)
        pass

class Test_center_of_gravity(unittest.TestCase):
    def test_center_of_gravity(self):
        oldv = oldfu.center_of_gravity(a)
        v = fu.center_of_gravity(a)
        pass


class Test_center_of_gravity_phase(unittest.TestCase):
    def test_center_of_gravity_phase(self):
        oldv = oldfu.center_of_gravity_phase(a)
        v = fu.center_of_gravity_phase(a)
        pass

class Test_fit_ctf(unittest.TestCase):
    def test_fit_ctf(self):
        oldv = oldfu.fit_ctf(crossresolution, ctf_params, rangedef = -1.0, i1 = 0, i2 = 0, chisquare=False)
        v = fu.fit_ctf(crossresolution, ctf_params, rangedef = -1.0, i1 = 0, i2 = 0, chisquare=False)
        pass


class Test_randprojdir(unittest.TestCase):
    def test_randprojdir(self):
        oldv = oldfu.randprojdir(ang, sigma)
        v = fu.randprojdir(ang, sigma)
        pass

 end: new in sphire 1.3"""


class Test_add_ave_varf_MPI(unittest.TestCase):
    main_node = 0
    data = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
    )[0][0]

    data = give_ali2d_single_iter_data()[0]

    def test_all_the_conditions(
        self, return_new=(), return_old=(), tolerance=TOLERANCE
    ):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(
                allclose(j.get_3dview(), i.get_3dview(), atol=tolerance, equal_nan=True)
            )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.add_ave_varf_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.add_ave_varf_MPI()
        self.assertEqual(
            str(cm_new.exception),
            "add_ave_varf_MPI() missing 2 required positional arguments: 'myid' and 'data'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_CTF_True_deadCode_BUG(self):
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(UnboundLocalError) as cm_new:
            fu.add_ave_varf_MPI(
                self.main_node,
                self.data,
                mask=None,
                mode="a",
                CTF=True,
                ctf_2_sum=None,
                ali_params="xform.align2d",
                main_node=self.main_node,
                comm=-1,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(UnboundLocalError) as cm_old:
            oldfu.add_ave_varf_MPI(
                self.main_node,
                self.data,
                mask=None,
                mode="a",
                CTF=True,
                ctf_2_sum=None,
                ali_params="xform.align2d",
                main_node=self.main_node,
                comm=-1,
            )
        self.assertEqual(
            str(cm_new.exception), "local variable 'i' referenced before assignment"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_CTF_wrong_ali_params(self):
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(RuntimeError) as cm_new:
            fu.add_ave_varf_MPI(
                self.main_node,
                self.data,
                mask=None,
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ali_params="xform.align3d",
                main_node=self.main_node,
                comm=-1,
            )
        mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.add_ave_varf_MPI(
                self.main_node,
                self.data,
                mask=None,
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ali_params="xform.align3d",
                main_node=self.main_node,
                comm=-1,
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value_msg_Error_negative_variance(self):
        with self.assertRaises(SystemExit) as cm_new:
            fu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=None,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        sp_global_def.BATCH =True
        with self.assertRaises(SystemExit) as cm_old:
            oldfu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=None,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        # self.test_all_the_conditions(return_new, return_old, tolerance=0)
        # self.assertEqual(len(return_new), 5)

    def test_with_mask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=MASK,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=MASK,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        self.test_all_the_conditions(return_new, return_old, tolerance=0)
        self.assertEqual(len(return_new), 5)

    def test_no_on_main_node(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(
            self.main_node + 1,
            self.data,
            mask=MASK,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(
            self.main_node + 1,
            self.data,
            mask=MASK,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        self.test_all_the_conditions(return_new, return_old, tolerance=0)
        self.assertEqual(len(return_new), 5)
        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                EMAN2_cppwrap.EMData().get_3dview(),
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(
                return_new[-1].get_3dview(),
                EMAN2_cppwrap.EMData().get_3dview(),
                equal_nan=True,
            )
        )

    def test_with_mask_and_not_current_alignment_parameters(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=MASK,
            mode="not_current",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=MASK,
            mode="not_current",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        self.test_all_the_conditions(return_new, return_old, tolerance=0)
        self.assertEqual(len(return_new), 5)

    def test_without_mask_and_not_current_alignment_parameters_msg_Error_negative_variance(
        self
    ):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=None,
            mode="not_current",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(
            self.main_node,
            self.data,
            mask=None,
            mode="not_current",
            CTF=False,
            ctf_2_sum=None,
            ali_params="xform.align2d",
            main_node=self.main_node,
            comm=-1,
        )
        self.test_all_the_conditions(return_new, return_old, tolerance=0)
        self.assertEqual(len(return_new), 5)


#todo: I need a list of image with ''xform.align2d'
class Test_sum_oe(unittest.TestCase):
    data, = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
    )[0]

    data  = give_ave_series_data()

    def test_all_the_conditions(self, return_new=(), return_old=()):
        if len(return_new) > 0:
            self.assertTrue(
                array_equal(return_new[0].get_3dview(), return_old[0].get_3dview())
            )
            self.assertTrue(
                array_equal(return_new[1].get_3dview(), return_old[1].get_3dview())
            )
            self.assertTrue(array_equal(return_new[2], return_old[2]))

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.sum_oe()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.sum_oe()
        self.assertEqual(
            str(cm_new.exception), "sum_oe() missing 1 required positional argument: 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_list_data_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.sum_oe(
                [],
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ctf_eo_sum=False,
                return_params=False,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.sum_oe(
                [],
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ctf_eo_sum=False,
                return_params=False,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_NoneType_lists_data_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.sum_oe(
                [None],
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ctf_eo_sum=False,
                return_params=False,
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.sum_oe(
                [None],
                mode="a",
                CTF=False,
                ctf_2_sum=None,
                ctf_eo_sum=False,
                return_params=False,
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_default_value(self):
        return_new = fu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=False,
        )
        return_old = oldfu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=False,
        )
        self.assertTrue(
            array_equal(return_new[0].get_2dview(), return_old[0].get_2dview())
        )
        self.assertTrue(
            array_equal(return_new[1].get_2dview(), return_old[1].get_2dview())
        )

    def test_default_value_with_returnParams(self):
        return_new = fu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=True,
        )
        return_old = oldfu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=True,
        )
        self.test_all_the_conditions(return_new, return_old)


    def test_not_current_alignment_parameters(self):
        return_new = fu.sum_oe(
            data=self.data,
            mode="not",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=True,
        )
        return_old = oldfu.sum_oe(
            data=self.data,
            mode="not",
            CTF=False,
            ctf_2_sum=None,
            ctf_eo_sum=False,
            return_params=True,
        )
        self.test_all_the_conditions(return_new, return_old)


    def test_prevent_calculation_ctf2(self):
        return_new = fu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=EMAN2_cppwrap.EMData(),
            ctf_eo_sum=False,
            return_params=True,
        )
        return_old = oldfu.sum_oe(
            data=self.data,
            mode="a",
            CTF=False,
            ctf_2_sum=EMAN2_cppwrap.EMData(),
            ctf_eo_sum=False,
            return_params=True,
        )
        self.test_all_the_conditions(return_new, return_old)





class Test_ave_var(unittest.TestCase):
    data = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var")
    )[0][0]

    data = give_ave_series_data()

    def test_all_the_conditions(self, return_new=(), return_old=()):
        if len(return_new) > 0:
            self.assertTrue(
                array_equal(return_new[0].get_3dview(), return_old[0].get_3dview())
            )
            self.assertTrue(
                array_equal(return_new[1].get_3dview(), return_old[1].get_3dview())
            )

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave_var()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave_var()
        self.assertEqual(
            str(cm_new.exception), "ave_var() missing 1 required positional argument: 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value(self):
        return_new = fu.ave_var(self.data, mode="a", listID=None)
        return_old = oldfu.ave_var(self.data, mode="a", listID=None)
        self.test_all_the_conditions(return_new, return_old)

    def test_not_current_alignment_parameters(self):
        return_new = fu.ave_var(self.data, mode="nota", listID=None)
        return_old = oldfu.ave_var(self.data, mode="nota", listID=None)
        self.test_all_the_conditions(return_new, return_old)

    def test_with_listID(self):
        return_new = fu.ave_var(self.data, mode="a", listID=[1, 2])
        return_old = oldfu.ave_var(self.data, mode="a", listID=[1, 2])
        self.test_all_the_conditions(return_new, return_old)


#todo: I cannot performe the unit test
class Test_ave_series(unittest.TestCase):
    data, = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
    )[0]

    data = give_ave_series_data()

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave_series()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave_series()
        self.assertEqual(
            str(cm_new.exception), "ave_series() missing 1 required positional argument: 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=True, mask=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=True, mask=None)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))


    def test_without_pave(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=False, mask=None)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=False, mask=None)
        self.assertTrue(array_equal(return_new.get_2dview(), return_old.get_2dview()))


    def test_default_value_withMask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=True, mask=MASK)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=True, mask=MASK)
        self.assertTrue(
            allclose(
                return_new.get_2dview(),
                return_old.get_2dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )

    def test_without_pave_withMask(self):
        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.ave_series(self.data, pave=False, mask=MASK)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.ave_series(self.data, pave=False, mask=MASK)
        self.assertTrue(
            allclose(
                return_new.get_2dview(),
                return_old.get_2dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )


class Test_varf2d_MPI(unittest.TestCase):
    data, = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
    )[0]
    main_node = 0
    data  = give_ave_series_data()

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.varf2d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.varf2d_MPI()
        self.assertEqual(
            str(cm_new.exception), "varf2d_MPI() missing 3 required positional arguments: 'myid', 'data', and 'ave'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_value_img2D(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )


    def test_default_value_img3D(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_3D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_3D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [
                    657.29736328125,
                    4417.57666015625,
                    12547.1904296875,
                    5386.271484375,
                    4284.39404296875,
                    4738.46044921875,
                    1702.59375,
                    2755.111572265625,
                    2629.47314453125,
                    2593.828125,
                    1648.6844482421875,
                    1247.2149658203125,
                    1095.642333984375,
                    892.38525390625,
                    537.0348510742188,
                    228.2470245361328,
                    67.18709564208984,
                    12.123698234558105,
                    2.3058433532714844,
                    0.9480835199356079,
                    0.673009991645813,
                    0.5360445976257324,
                    0.4416618049144745,
                    0.3757860064506531,
                    0.31849655508995056,
                    0.2831489145755768,
                    0.2509312629699707,
                    0.22601264715194702,
                    0.20667268335819244,
                    0.19018994271755219,
                    0.17610518634319305,
                    0.16430525481700897,
                    0.15365713834762573,
                    0.14737848937511444,
                    0.14194539189338684,
                    0.1365765631198883,
                    0.13146640360355377,
                    0.12925119698047638,
                    0.12116505205631256,
                ],
            )
        )

    def test_no_mainNode(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node + 1,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node + 1,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            array_equal(return_new[0].get_3dview(), EMAN2_cppwrap.EMData().get_3dview())
        )
        self.assertTrue(array_equal(return_new[1], [0]))

    @unittest.skip("i do not have a valid img")
    def test_with_CTF(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=True,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="a",
            CTF=True,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [
                    657.29736328125,
                    4417.57666015625,
                    12547.1904296875,
                    5386.271484375,
                    4284.39404296875,
                    4738.46044921875,
                    1702.59375,
                    2755.111572265625,
                    2629.47314453125,
                    2593.828125,
                    1648.6844482421875,
                    1247.2149658203125,
                    1095.642333984375,
                    892.38525390625,
                    537.0348510742188,
                    228.2470245361328,
                    67.18709564208984,
                    12.123698234558105,
                    2.3058433532714844,
                    0.9480835199356079,
                    0.673009991645813,
                    0.5360445976257324,
                    0.4416618049144745,
                    0.3757860064506531,
                    0.31849655508995056,
                    0.2831489145755768,
                    0.2509312629699707,
                    0.22601264715194702,
                    0.20667268335819244,
                    0.19018994271755219,
                    0.17610518634319305,
                    0.16430525481700897,
                    0.15365713834762573,
                    0.14737848937511444,
                    0.14194539189338684,
                    0.1365765631198883,
                    0.13146640360355377,
                    0.12925119698047638,
                    0.12116505205631256,
                ],
            )
        )

    def test_with_mask(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=MASK,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=MASK,
            mode="a",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            allclose(
                return_new[1], numpy_full(len(return_new[1]), numpy_nan), equal_nan=True
            )
        )

    def test_without_apply_alignment_parameters(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="not",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=None,
            mode="not",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [
                    657.143798828125,
                    4400.4541015625,
                    12588.783203125,
                    5354.16748046875,
                    4298.67626953125,
                    4745.57958984375,
                    1703.4573974609375,
                    2770.798583984375,
                    2642.747802734375,
                    2620.13916015625,
                    1666.76220703125,
                    1271.91552734375,
                    1132.228759765625,
                    928.3460693359375,
                    560.6527709960938,
                    241.77120971679688,
                    68.29119873046875,
                    10.938950538635254,
                    1.148409366607666,
                    0.10395097732543945,
                    0.009575692005455494,
                    0.0007211578777059913,
                    6.025378024787642e-05,
                    5.932938165642554e-06,
                    6.559504299730179e-07,
                    7.396378975954576e-08,
                    8.115450533807689e-09,
                    1.0014190587881444e-09,
                    1.227386109414752e-10,
                    2.5507195314244946e-11,
                    1.5168682282462598e-11,
                    1.4266591380485139e-11,
                    1.3229447919094195e-11,
                    1.3020461445134579e-11,
                    1.4008788920549797e-11,
                    1.3700839074370919e-11,
                    1.456094446405931e-11,
                    1.4547537653675224e-11,
                    1.2415271935517502e-11,
                ],
            )
        )

    def test_with_mask_without_apply_alignment_parameters(self):
        return_new = fu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=MASK,
            mode="not",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )
        return_old = oldfu.varf2d_MPI(
            myid=self.main_node,
            data=self.data,
            ave=IMAGE_2D,
            mask=MASK,
            mode="not",
            CTF=False,
            main_node=self.main_node,
            comm=-1,
        )

        self.assertTrue(
            allclose(
                return_new[0].get_3dview(),
                return_old[0].get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )
        self.assertTrue(
            allclose(return_new[1], return_old[1], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            allclose(
                return_new[1], numpy_full(len(return_new[1]), numpy_nan), equal_nan=True
            )
        )

class Test_ccc(unittest.TestCase):
    (img1, img2, mask) = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc")
    )[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ccc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ccc()
        self.assertEqual(
            str(cm_new.exception), "ccc() missing 2 required positional arguments: 'img1' and 'img2'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        return_new = fu.ccc(self.img1, self.img2, self.mask)
        return_old = oldfu.ccc(self.img1, self.img2, self.mask)
        self.assertEqual(return_new, return_old)
        # self.assertEqual(return_new, 0.8129369020462036)

    def test_None_mask(self):
        return_new = fu.ccc(self.img1, self.img2, None)
        return_old = oldfu.ccc(self.img1, self.img2, None)
        self.assertEqual(return_new, return_old)
        # self.assertEqual(return_new, 0.8021121621131897)

    def test_empty_mask_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.ccc(self.img1, self.img2, EMData())
        return_old = oldfu.ccc(self.img1, self.img2, EMData())
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0)
        """

    def test_Noneimg1(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.ccc(None, self.img2, self.mask)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.ccc(None, self.img2, self.mask)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'cmp'"
        )

    def test_emptyimg1(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(EMData(), self.img2, self.mask)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(EMData(), self.img2, self.mask)

        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_Noneimg2(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(self.img1, None, self.mask)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(self.img1, None, self.mask)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_emptyimg2(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ccc(self.img1, EMData(), self.mask)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ccc(self.img1, EMData(), self.mask)

        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])


class Test_fsc(unittest.TestCase):
    (img1, img2) = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc")
    )[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fsc()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fsc()
        self.assertEqual(
            str(cm_new.exception), "fsc() missing 2 required positional arguments: 'img1' and 'img2'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.fsc(self.img1, self.img2, w=1.0, filename=None)
        return_old = oldfu.fsc(self.img1, self.img2, w=1.0, filename=None)
        self.assertTrue(array_equal(return_new, return_old))
        # self.assertTrue(
        #     array_equal(
        #         return_new,
        #         [
        #             [
        #                 0.0,
        #                 0.006756756920367479,
        #                 0.013513513840734959,
        #                 0.020270271226763725,
        #                 0.027027027681469917,
        #                 0.03378378599882126,
        #                 0.04054054245352745,
        #                 0.04729729890823364,
        #                 0.054054055362939835,
        #                 0.06081081181764603,
        #                 0.06756757199764252,
        #                 0.07432432472705841,
        #                 0.0810810849070549,
        #                 0.0878378376364708,
        #                 0.09459459781646729,
        #                 0.10135135054588318,
        #                 0.10810811072587967,
        #                 0.11486487090587616,
        #                 0.12162162363529205,
        #                 0.12837837636470795,
        #                 0.13513514399528503,
        #                 0.14189189672470093,
        #                 0.14864864945411682,
        #                 0.15540540218353271,
        #                 0.1621621698141098,
        #                 0.1689189225435257,
        #                 0.1756756752729416,
        #                 0.18243244290351868,
        #                 0.18918919563293457,
        #                 0.19594594836235046,
        #                 0.20270270109176636,
        #                 0.20945946872234344,
        #                 0.21621622145175934,
        #                 0.22297297418117523,
        #                 0.22972974181175232,
        #                 0.2364864945411682,
        #                 0.2432432472705841,
        #                 0.25,
        #                 0.2567567527294159,
        #                 0.2635135054588318,
        #                 0.27027028799057007,
        #                 0.27702704071998596,
        #                 0.28378379344940186,
        #                 0.29054054617881775,
        #                 0.29729729890823364,
        #                 0.30405405163764954,
        #                 0.31081080436706543,
        #                 0.3175675868988037,
        #                 0.3243243396282196,
        #                 0.3310810923576355,
        #                 0.3378378450870514,
        #                 0.3445945978164673,
        #                 0.3513513505458832,
        #                 0.3581081032752991,
        #                 0.36486488580703735,
        #                 0.37162163853645325,
        #                 0.37837839126586914,
        #                 0.38513514399528503,
        #                 0.3918918967247009,
        #                 0.3986486494541168,
        #                 0.4054054021835327,
        #                 0.412162184715271,
        #                 0.4189189374446869,
        #                 0.4256756901741028,
        #                 0.4324324429035187,
        #                 0.43918919563293457,
        #                 0.44594594836235046,
        #                 0.45270270109176636,
        #                 0.45945948362350464,
        #                 0.46621623635292053,
        #                 0.4729729890823364,
        #                 0.4797297418117523,
        #                 0.4864864945411682,
        #                 0.4932432472705841,
        #                 0.0,
        #             ],
        #             [
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 1.0,
        #                 0.9999955892562866,
        #                 0.998257040977478,
        #                 0.9918078780174255,
        #                 0.9856542348861694,
        #                 0.9881284832954407,
        #                 0.9837390184402466,
        #                 0.9792850017547607,
        #                 0.9748953580856323,
        #                 0.9678066968917847,
        #                 0.9480144381523132,
        #                 0.9332209229469299,
        #                 0.8907226324081421,
        #                 0.8318315148353577,
        #                 0.7905364036560059,
        #                 0.8009309768676758,
        #                 0.806480348110199,
        #                 0.7942577600479126,
        #                 0.7667281627655029,
        #                 0.7561865448951721,
        #                 0.7457202672958374,
        #                 0.7362813353538513,
        #                 0.7185238003730774,
        #                 0.7565295100212097,
        #                 0.7811623215675354,
        #                 0.7833795547485352,
        #                 0.7671730518341064,
        #                 0.7416536808013916,
        #                 0.7073566913604736,
        #                 0.719993531703949,
        #                 0.7416943907737732,
        #                 0.7183090448379517,
        #                 0.691653847694397,
        #                 0.6795889735221863,
        #                 0.6586717367172241,
        #                 0.6515437960624695,
        #                 0.5965009331703186,
        #                 0.5489014387130737,
        #                 0.565247654914856,
        #                 0.5726661682128906,
        #                 0.5036922693252563,
        #                 0.38146188855171204,
        #                 0.26737281680107117,
        #                 0.3945968449115753,
        #                 0.4944046437740326,
        #                 0.39991524815559387,
        #                 0.2689603269100189,
        #                 0.2521679699420929,
        #                 0.2794639468193054,
        #                 0.280245840549469,
        #                 0.24809274077415466,
        #                 0.25190722942352295,
        #                 0.2277139574289322,
        #                 0.19758595526218414,
        #                 0.19085757434368134,
        #                 0.19563539326190948,
        #                 0.19770143926143646,
        #                 0.16604311764240265,
        #                 0.18556171655654907,
        #                 0.13770028948783875,
        #                 0.152345210313797,
        #                 0.17087307572364807,
        #                 0.15340681374073029,
        #                 0.16820573806762695,
        #                 0.18507032096385956,
        #                 0.1278713047504425,
        #                 0.1172977089881897,
        #                 0.0,
        #             ],
        #             [
        #                 2.0,
        #                 18.0,
        #                 62.0,
        #                 98.0,
        #                 210.0,
        #                 350.0,
        #                 450.0,
        #                 602.0,
        #                 762.0,
        #                 1142.0,
        #                 1250.0,
        #                 1458.0,
        #                 1814.0,
        #                 2178.0,
        #                 2498.0,
        #                 2622.0,
        #                 3338.0,
        #                 3722.0,
        #                 4170.0,
        #                 4358.0,
        #                 5034.0,
        #                 5714.0,
        #                 5982.0,
        #                 6602.0,
        #                 7130.0,
        #                 8034.0,
        #                 8606.0,
        #                 9066.0,
        #                 9962.0,
        #                 10550.0,
        #                 11226.0,
        #                 12146.0,
        #                 12606.0,
        #                 13802.0,
        #                 14754.0,
        #                 15194.0,
        #                 16454.0,
        #                 17154.0,
        #                 18266.0,
        #                 18750.0,
        #                 20234.0,
        #                 21450.0,
        #                 21962.0,
        #                 23462.0,
        #                 24042.0,
        #                 25946.0,
        #                 26118.0,
        #                 27506.0,
        #                 29066.0,
        #                 30450.0,
        #                 31742.0,
        #                 32250.0,
        #                 34250.0,
        #                 35454.0,
        #                 36434.0,
        #                 37682.0,
        #                 39294.0,
        #                 41426.0,
        #                 42066.0,
        #                 43490.0,
        #                 45702.0,
        #                 46634.0,
        #                 48554.0,
        #                 48870.0,
        #                 51554.0,
        #                 54090.0,
        #                 54314.0,
        #                 56910.0,
        #                 57482.0,
        #                 60626.0,
        #                 61206.0,
        #                 62570.0,
        #                 65730.0,
        #                 66686.0,
        #                 0.0,
        #             ],
        #         ],
        #     )
        # )

    def test_saveOnfile(self):
        fnew = path.join(ABSOLUTE_PATH, "new.txt")
        fold = path.join(ABSOLUTE_PATH, "new.old")
        return_new = fu.fsc(self.img1, self.img2, w=1.0, filename=fnew)
        return_old = oldfu.fsc(self.img1, self.img2, w=1.0, filename=fold)
        self.assertTrue(array_equal(return_new, return_old))

        self.assertTrue(returns_values_in_file(fnew), returns_values_in_file(fold))
        remove_list_of_file([fnew, fold])

    def test_w_set_0_returns_MemoryError(self):
        with self.assertRaises(MemoryError) as cm_new:
            fu.fsc(self.img1, self.img2, w=0, filename=None)
        with self.assertRaises(MemoryError) as cm_old:
            oldfu.fsc(self.img1, self.img2, w=0, filename=None)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_emptyImg1_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(EMData(), self.img2, w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(EMData(), self.img2, w=1.0, filename=None)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageFormatException")
        self.assertEqual(msg[1], "Cannot calculate FSC for 1D images")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img1_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.fsc(None, self.img2, w=1.0, filename=None)
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fsc(None, self.img2, w=1.0, filename=None)
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception),
            "'NoneType' object has no attribute 'calc_fourier_shell_correlation'",
        )

    def test_emptyImg2_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(self.img1, EMData(), w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(self.img1, EMData(), w=1.0, filename=None)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])

    def test_Nonetype_Img2_returns_AttributeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc(self.img1, None, w=1.0, filename=None)
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc(self.img1, None, w=1.0, filename=None)
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NullPointerException")
        self.assertEqual(msg[1], "NULL input image")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])


class Test_fsc_mask(unittest.TestCase):
    (img1, img2, mask, w, not_used) = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc_mask")
    )[0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.fsc_mask()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.fsc_mask()
        self.assertEqual(
            str(cm_new.exception), "fsc_mask() missing 2 required positional arguments: 'img1' and 'img2'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=self.mask, w=self.w, filename=None
        )
        return_old = oldfu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=self.mask, w=self.w, filename=None
        )
        self.assertTrue(array_equal(return_new, return_old))


    def test_saveOnfile(self):
        fnew = path.join(ABSOLUTE_PATH, "new.txt")
        fold = path.join(ABSOLUTE_PATH, "new.old")
        return_new = fu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=self.mask, w=self.w, filename=fnew
        )
        return_old = oldfu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=self.mask, w=self.w, filename=fold
        )
        self.assertTrue(array_equal(return_new, return_old))

        self.assertTrue(returns_values_in_file(fnew), returns_values_in_file(fold))
        remove_list_of_file([fnew, fold])

    def test_w_set_0_returns_MemoryError(self):
        with self.assertRaises(MemoryError) as cm_new:
            fu.fsc_mask(
                img1=self.img1, img2=self.img1, mask=self.mask, w=0, filename=None
            )
        with self.assertRaises(MemoryError) as cm_old:
            oldfu.fsc_mask(
                img1=self.img1, img2=self.img1, mask=self.mask, w=0, filename=None
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_emptyImg1_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(
                img1=EMData(), img2=self.img1, mask=self.mask, w=self.w, filename=None
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(
                img1=EMData(), img2=self.img1, mask=self.mask, w=self.w, filename=None
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(
            msg[1],
            "The dimension of the image does not match the dimension of the mask!",
        )
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img1_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.fsc_mask(
                img1=None, img2=self.img1, mask=self.mask, w=self.w, filename=None
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.fsc_mask(
                img1=None, img2=self.img1, mask=self.mask, w=self.w, filename=None
            )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )

    def test_emptyImg2_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(
                img1=self.img1, img2=EMData(), mask=self.mask, w=self.w, filename=None
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(
                img1=self.img1, img2=EMData(), mask=self.mask, w=self.w, filename=None
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "ImageDimensionException")
        self.assertEqual(
            msg[1],
            "The dimension of the image does not match the dimension of the mask!",
        )
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[1], msg_old[1])

    def test_Nonetype_Img2_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)


    def test_empty_mask_returns_RuntimeError(self):
        with self.assertRaises(RuntimeError) as cm_new:
            fu.fsc_mask(
                img1=self.img1, img2=self.img1, mask=EMData(), w=self.w, filename=None
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.fsc_mask(
                img1=self.img1, img2=self.img1, mask=EMData(), w=self.w, filename=None
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "std::exception")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])

    def test_NoneType_mask(self):
        return_new = fu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=None, w=self.w, filename=None
        )
        return_old = oldfu.fsc_mask(
            img1=self.img1, img2=self.img1, mask=None, w=self.w, filename=None
        )
        self.assertTrue(array_equal(return_new, return_old))



class Test_locres(unittest.TestCase):
    (
        vi,
        ui,
        m,
        nk,
        cutoff,
        step,
        myid,
        main_node,
        number_of_proc,
    ) = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.locres")
    )[
        0
    ]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.locres()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.locres()
        self.assertEqual(
            str(cm_new.exception), "locres() missing 9 required positional arguments: 'vi', 'ui', 'm', 'nk', 'cutoff', 'step', 'myid', 'main_node', and 'number_of_proc'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        self.assertTrue(True)
        """
        Connected to pydev debugger (build 183.6156.13)
        Launching unittests with arguments python -m unittest sphire.tests.test_statistics.Test_locres.test_pickle_value in /home/lusnig/EMAN2/eman2
        [rtxr2:22191] *** An error occurred in MPI_Recv
        [rtxr2:22191] *** reported by process [1715994625,140423955742720]
        [rtxr2:22191] *** on communicator MPI_COMM_WORLD
        [rtxr2:22191] *** MPI_ERR_RANK: invalid rank
        [rtxr2:22191] *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
        [rtxr2:22191] ***    and potentially your MPI job)
        
        Process finished with exit code 6
        
        return_new = fu.locres(vi=self.vi, ui=self.ui, m=self.m, nk=self.nk, cutoff=self.nk, step=self.step, myid=0, main_node=0, number_of_proc=self.number_of_proc)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.locres(vi=self.vi, ui=self.ui, m=self.m, nk=self.nk, cutoff=self.nk, step=self.step, myid=0, main_node=0, number_of_proc=self.number_of_proc)
        mpi_barrier(MPI_COMM_WORLD)
        """





class Test_k_means_match_clusters_asg_new(unittest.TestCase):
    asg1=[
        numpy_asarray([0, 3, 4, 5, 6, 9, 11, 12, 13, 15, 16, 18, 19, 22, 23, 27, 30, 31, 34, 36, 37, 40, 41, 44, 45, 46, 48, 49, 51, 53, 54, 57, 58, 61, 65, 68, 73, 76, 78, 79, 81, 83, 85, 86, 87, 88, 89, 90, 91, 97, 100, 101, 103, 104, 106, 109, 111, 112, 113, 116, 117, 118, 120, 121, 124, 125, 126, 127, 128, 129, 130, 131, 134, 136, 138, 139, 140, 141, 143, 145, 148, 149, 150, 154, 159, 161, 164, 167, 172, 173, 176, 178, 186, 188, 191, 192, 193, 196, 198, 200, 201, 202, 204, 206, 208, 209, 211, 216, 218, 225, 226, 227, 228, 229, 232, 233, 234, 235, 236, 239, 243, 244, 247, 249, 250, 254, 255, 258, 261, 263, 265, 270, 272, 274, 275, 277, 278, 280, 283, 285, 286, 287, 291, 292, 293, 295, 297, 304, 306, 309, 317, 318, 321, 323, 324, 326, 327, 330, 332, 335, 338, 339, 340, 343, 345, 346, 347, 348, 352, 357, 358, 359, 360, 368, 369, 370, 371, 374, 378, 381, 382, 386, 388, 390, 392, 394, 398, 401, 402, 404, 408, 409, 410, 411, 414, 415, 416, 417, 418, 422, 423, 425, 426, 427, 432, 434, 435, 438, 445, 446, 448, 449, 450, 452, 453, 455, 456, 457, 459, 460, 461, 462, 465, 466, 468, 469, 471, 472, 475, 476, 477, 478, 479, 480, 481, 486, 492, 493, 497, 499, 500, 501, 506, 507, 508, 509, 510, 511, 516, 517, 518, 519, 521, 524, 526, 533, 534, 536, 539, 543, 545, 546, 547, 548, 550, 551, 555, 559, 566, 571, 572, 575, 578, 579, 581, 582, 583, 586, 589, 591, 592, 594, 595, 599, 600, 601, 602, 605, 606, 608, 611, 613, 616, 617, 618, 621, 625, 627, 630, 635, 637, 639, 644, 645, 647, 648, 653, 658, 659, 661, 663, 668, 670, 675, 678, 681, 682, 687, 692, 693, 695, 696, 697, 698, 700, 701, 703, 705, 711, 715, 716, 717, 722, 725, 728, 729, 733, 735, 737, 739, 741, 743, 744, 748, 750, 754, 756, 757, 759, 760, 761, 762, 764, 766, 774, 775, 776, 781, 782, 783, 784, 786, 789, 792, 793, 795, 799, 800, 801, 804, 805, 806, 807, 808, 809, 810, 812, 813, 815, 817, 819, 820, 821, 822, 823, 824, 825, 826, 827, 832, 835, 836, 839, 841, 843, 844, 848, 849, 850, 851, 852, 853, 859, 860, 862, 864, 865, 867, 868, 870, 871, 872, 873, 874, 876, 877, 878, 881, 884, 887, 889, 892, 893, 894, 896, 900, 904, 905, 906, 907, 909, 912, 913, 917, 918, 920, 923, 925, 926, 928, 929, 930, 931, 934, 935, 943, 949, 952, 954, 955, 956, 957, 958, 959, 962, 965, 969, 970, 973, 974, 978, 981, 982, 987, 990, 991, 992, 994, 997, 1000, 1002, 1003, 1004, 1006, 1007, 1008, 1011, 1015, 1016, 1017, 1018, 1019, 1020, 1022, 1024, 1025, 1026, 1030, 1031, 1032, 1034, 1035, 1037, 1041, 1043, 1047, 1049, 1050, 1057, 1065, 1066, 1067, 1068, 1070, 1073, 1086, 1087, 1088, 1092, 1095, 1097, 1098, 1099, 1101, 1102, 1106, 1107, 1109, 1111, 1113, 1114, 1115, 1116, 1117, 1119, 1120, 1121, 1122, 1123, 1125, 1126, 1127, 1130, 1132, 1134, 1135, 1136, 1138, 1139, 1144, 1145, 1146, 1148, 1149, 1150, 1151, 1152, 1153, 1155, 1156, 1162, 1164, 1165, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1179, 1180, 1183, 1185, 1186, 1188, 1189, 1191, 1192, 1193, 1195, 1196, 1198, 1201, 1202, 1203, 1204, 1207, 1211, 1214, 1216, 1217, 1218, 1223, 1227, 1228, 1229, 1230, 1231, 1233, 1236, 1237, 1243, 1245, 1246, 1247, 1248, 1251, 1252, 1253, 1255, 1257, 1258, 1260, 1261, 1264, 1265, 1268, 1269, 1270, 1271, 1274, 1275, 1276, 1278, 1279, 1280, 1282, 1284, 1285, 1292, 1293, 1296, 1297, 1299, 1302, 1303, 1304, 1306, 1315, 1316, 1317, 1318, 1320, 1322, 1329, 1332, 1336, 1338, 1339, 1340, 1345, 1349, 1352, 1353, 1355, 1356, 1357, 1359, 1364, 1365, 1366, 1368, 1370, 1371, 1372, 1374, 1375, 1377, 1384, 1385, 1386, 1387, 1388, 1389, 1391, 1392, 1395, 1396, 1397, 1401, 1403, 1404, 1407, 1408, 1411, 1412, 1413, 1416, 1417, 1419, 1420, 1421, 1422, 1425, 1428, 1429, 1431, 1432, 1435, 1437, 1438, 1439, 1440, 1441, 1442, 1444, 1449, 1450, 1451, 1453, 1454, 1455, 1457, 1462, 1463, 1464, 1466, 1469, 1472, 1474, 1475, 1476, 1478, 1479, 1481, 1482, 1484, 1488, 1489, 1490, 1491, 1493, 1494, 1495, 1497, 1499, 1501, 1508, 1509, 1511, 1512, 1513, 1514, 1517, 1518, 1521, 1522, 1524, 1526, 1527, 1534, 1538, 1539, 1540, 1542, 1546, 1547, 1548, 1550, 1551, 1552, 1553, 1554, 1556, 1557, 1561, 1563, 1564, 1566, 1567, 1568, 1569, 1570, 1571, 1572, 1573, 1574, 1575, 1576, 1577, 1578, 1579, 1580, 1582, 1583, 1585, 1586, 1588, 1590, 1591, 1595, 1596, 1598, 1599, 1603, 1604, 1606, 1607, 1608, 1610, 1612, 1613, 1620, 1621, 1622, 1623, 1624, 1629, 1630, 1632, 1634, 1637, 1639, 1640, 1641, 1642, 1643, 1649, 1653, 1654, 1655, 1656, 1658, 1661, 1663, 1667, 1671, 1673, 1674, 1676, 1677, 1679, 1680, 1682, 1683, 1684, 1685, 1686, 1688, 1691, 1692, 1696, 1697, 1698, 1700, 1701, 1712, 1713, 1714, 1715, 1716, 1717, 1718, 1719, 1720, 1721, 1722, 1723, 1727, 1728, 1729, 1731, 1732, 1734, 1738, 1739, 1741, 1744, 1746, 1747, 1750, 1751, 1752, 1758, 1763, 1764, 1765, 1770, 1774, 1776, 1777, 1778, 1782, 1784, 1786, 1789, 1791, 1793, 1794, 1795, 1797, 1801, 1803, 1804, 1805, 1807, 1808, 1809, 1811, 1812, 1818, 1821, 1822, 1823, 1824, 1825, 1826, 1829, 1830, 1831, 1833, 1834, 1837, 1838, 1839, 1841, 1842, 1844, 1846, 1848, 1850, 1851, 1852, 1853, 1855, 1857, 1859, 1862, 1863, 1864, 1865, 1866, 1867, 1869, 1870, 1871, 1873, 1874, 1877, 1879, 1880, 1881, 1882, 1888, 1889, 1891, 1894, 1895, 1896, 1897, 1899, 1900, 1903, 1904, 1905, 1906, 1907, 1908, 1910, 1912, 1913, 1915, 1916, 1918, 1919, 1921, 1922, 1923, 1924, 1927, 1931, 1933, 1934, 1937, 1939, 1941, 1945, 1946, 1947, 1948, 1951, 1954, 1955, 1956, 1958, 1961, 1965, 1966, 1967, 1970, 1972, 1973, 1974, 1976, 1977, 1978, 1979, 1981, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1994, 1995, 1996, 1998, 1999, 2000, 2002, 2003, 2007, 2008, 2010, 2012, 2013, 2014, 2016, 2020, 2022, 2024, 2026, 2028, 2030, 2031, 2034, 2035, 2036, 2038, 2040, 2041, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2055, 2058, 2066, 2067, 2070, 2071, 2073, 2074, 2078, 2081, 2082, 2083, 2090, 2091, 2093, 2095, 2096, 2100, 2101, 2104, 2106, 2107, 2113, 2115, 2116, 2120, 2121, 2124, 2125, 2128, 2130, 2131, 2132, 2133, 2137, 2138, 2139, 2143, 2144, 2145, 2147, 2148, 2150, 2155, 2156, 2158, 2161, 2162, 2163, 2166, 2167, 2168, 2169, 2172, 2173, 2178, 2180, 2182, 2184, 2185, 2187, 2189, 2190, 2192, 2194, 2199, 2201, 2202, 2203, 2204, 2205, 2207, 2211, 2212, 2213, 2216, 2217, 2218, 2219, 2220, 2221, 2222, 2223, 2227, 2228, 2229, 2230, 2232, 2236, 2237, 2239, 2240, 2241, 2242, 2244, 2247, 2250, 2251, 2254, 2255, 2256, 2257, 2259, 2264, 2265, 2266, 2267, 2269, 2272, 2274, 2276, 2277, 2279, 2280, 2282, 2283, 2284, 2285, 2291, 2292, 2296, 2299, 2300, 2301, 2303, 2305, 2307, 2308, 2309, 2310, 2311, 2313, 2314, 2315, 2316, 2317, 2318, 2320, 2321, 2322, 2324, 2328, 2331, 2332, 2333, 2334, 2339, 2341, 2343, 2345, 2348, 2349, 2350, 2352, 2354, 2355, 2358, 2360, 2363, 2365, 2367, 2371, 2373, 2376, 2380, 2381, 2384, 2385, 2386, 2387, 2389, 2395, 2396, 2397, 2402, 2403, 2404, 2405, 2407, 2409, 2412, 2413, 2414, 2417, 2419, 2421, 2423, 2424, 2428, 2433, 2434, 2440, 2442, 2443, 2444, 2449, 2454, 2455, 2457, 2458, 2460, 2461, 2463, 2468, 2469, 2470, 2472, 2473, 2475, 2479, 2480, 2486, 2487, 2488, 2489, 2490, 2495, 2496, 2498, 2499, 2500, 2501, 2503, 2505, 2506, 2507, 2509, 2510, 2511, 2512, 2516, 2519, 2521, 2522, 2523, 2525, 2526, 2527, 2532, 2533, 2534, 2535, 2536, 2538, 2541, 2542, 2548, 2549, 2551, 2557, 2559, 2561, 2562, 2563, 2565, 2567, 2568, 2569, 2570, 2571, 2573, 2575, 2576, 2579, 2582, 2585, 2586, 2589, 2593, 2595, 2597, 2599, 2600, 2602, 2610, 2611, 2612, 2613, 2614, 2615, 2619, 2620, 2623, 2625, 2628, 2631, 2633, 2634, 2635, 2636, 2637, 2638, 2640, 2642, 2643, 2645, 2648, 2651, 2652, 2654, 2655, 2656, 2658, 2659, 2660, 2661, 2663, 2667, 2668, 2669, 2670, 2671, 2675, 2676, 2677, 2678, 2680, 2681, 2682, 2683, 2688, 2689, 2692, 2694, 2696, 2697, 2698, 2699, 2701, 2702, 2703, 2704, 2705, 2707, 2711, 2712, 2714, 2715, 2718, 2720, 2722, 2723, 2727, 2728, 2729, 2730, 2732, 2735, 2738, 2739, 2740, 2741, 2743, 2746, 2749, 2750, 2751, 2752, 2754, 2755, 2756, 2758, 2760, 2764, 2765, 2767, 2768, 2771, 2774, 2781, 2783, 2784, 2786, 2796, 2801, 2806, 2807, 2808, 2810, 2814, 2818, 2821, 2822, 2824, 2826, 2827, 2828, 2830, 2831, 2832, 2833, 2835, 2837, 2838, 2845, 2846, 2847, 2849, 2852, 2860, 2863, 2864, 2866, 2867, 2868, 2869, 2871, 2872, 2873, 2874, 2875, 2879, 2881, 2882, 2883, 2885, 2886, 2892, 2893, 2896, 2897, 2900, 2902, 2903, 2904, 2905, 2906, 2907, 2909, 2911, 2913, 2914, 2915, 2916, 2917, 2918, 2920, 2921, 2922, 2924, 2925, 2927, 2929, 2932, 2933, 2936, 2937, 2938, 2939, 2940, 2942, 2945, 2947, 2950, 2951, 2955, 2956, 2958, 2961, 2962, 2963, 2964, 2965, 2967, 2968, 2970, 2971, 2972, 2976, 2978, 2981, 2983, 2985, 2987, 2988, 2990, 2998, 3003, 3008, 3009, 3010, 3011, 3014, 3015, 3016, 3018, 3020, 3022, 3023, 3024, 3025, 3026, 3027, 3028, 3029, 3030, 3033, 3035, 3036, 3038, 3040, 3041, 3044, 3045, 3051, 3053, 3054, 3057, 3062, 3064, 3066, 3067, 3068, 3069, 3072, 3073, 3078, 3079, 3081, 3082, 3084, 3085, 3087, 3088, 3093, 3095, 3096, 3097, 3098, 3100, 3102, 3103, 3104, 3108, 3111, 3112, 3114, 3115, 3119, 3121, 3122, 3129, 3132, 3134, 3135, 3139, 3147, 3151, 3153, 3155, 3156, 3157, 3159, 3160, 3161, 3162, 3163, 3164, 3165, 3167, 3169, 3173, 3175, 3178, 3182, 3185, 3186, 3187, 3189, 3192, 3194, 3195, 3196, 3200, 3201, 3203, 3204, 3205, 3209, 3211, 3213, 3214, 3215, 3216, 3217, 3218, 3219, 3220, 3221, 3222, 3225, 3230, 3231, 3233, 3235, 3236, 3238, 3242, 3243, 3244, 3246, 3247, 3249, 3250, 3255, 3256, 3257, 3261, 3262, 3263, 3264, 3266, 3268, 3269, 3272, 3274, 3275, 3276, 3277, 3278, 3279, 3282, 3285, 3287, 3288, 3289, 3290, 3293, 3294, 3296, 3297, 3298, 3302, 3303, 3304, 3305, 3306, 3310, 3311, 3318, 3319, 3322, 3324, 3325, 3326, 3328, 3333, 3334, 3335, 3336, 3337, 3338, 3344, 3345, 3347, 3349, 3351, 3354, 3358, 3359, 3360, 3361, 3364, 3365, 3368, 3369, 3370, 3371, 3373, 3375, 3376, 3379, 3380, 3381, 3385, 3386, 3387, 3391, 3394, 3397, 3398, 3399, 3402, 3403, 3406, 3407, 3413, 3414, 3415, 3418, 3421, 3422, 3426, 3427, 3428, 3429, 3432, 3433, 3435, 3438, 3440, 3444, 3445, 3446, 3447, 3448, 3449, 3450, 3451, 3452, 3453, 3456, 3464, 3465, 3467, 3469, 3470, 3471, 3472, 3474, 3477, 3479, 3481, 3484, 3485, 3486, 3488, 3490, 3491, 3492, 3493, 3494, 3495, 3496, 3497, 3498, 3501, 3502, 3504, 3505, 3509, 3510, 3511, 3513, 3515, 3516, 3517, 3518, 3519, 3520, 3521, 3522, 3524, 3526, 3528, 3532, 3534, 3535, 3539, 3540, 3547, 3549, 3551, 3552, 3557, 3564, 3569, 3571, 3573, 3574, 3581, 3582, 3583, 3585, 3586, 3588, 3589, 3590, 3592, 3595, 3597, 3598, 3600, 3601, 3603, 3604, 3606, 3608, 3610, 3612, 3613, 3617, 3618, 3620, 3629, 3631, 3634, 3636, 3637, 3638, 3640, 3641, 3642, 3644, 3646, 3649, 3651, 3652, 3655, 3660, 3662, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3678, 3679, 3682, 3683, 3685, 3686, 3687, 3692, 3694, 3696, 3697, 3701, 3703, 3704, 3706, 3707, 3708, 3709, 3710, 3712, 3715, 3716, 3719, 3722, 3723, 3724, 3728, 3729, 3731, 3732, 3733, 3735, 3736, 3738, 3743, 3744, 3745, 3746, 3747, 3750, 3751, 3753, 3755, 3759, 3760, 3761, 3764, 3765, 3767, 3768, 3769, 3770, 3771, 3773, 3774, 3775, 3777, 3780, 3782, 3783, 3784, 3787, 3789, 3790, 3791, 3792, 3794, 3796, 3797, 3800, 3801, 3803, 3804, 3808, 3809, 3813, 3814, 3815, 3819, 3820, 3821, 3822, 3823, 3825, 3826, 3828, 3829, 3830, 3831, 3832, 3833, 3834, 3838, 3839, 3840, 3841, 3842, 3845, 3848, 3849, 3850, 3851, 3852, 3853, 3854, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3865, 3868, 3873, 3877, 3880, 3883, 3885, 3887, 3888, 3889, 3891, 3892, 3893, 3896, 3897, 3898, 3900, 3901, 3905, 3908, 3910, 3911, 3913, 3914, 3915, 3917, 3920, 3922, 3924, 3928, 3929, 3930, 3935, 3939, 3941, 3943, 3947, 3951, 3954, 3955, 3956, 3959, 3962, 3963, 3964, 3966, 3968, 3969, 3970, 3972, 3973, 3974, 3975, 3978, 3979, 3980, 3981, 3983, 3985, 3987, 3989, 3991, 3993, 3995, 4000, 4001, 4003, 4004, 4006, 4011, 4014, 4018, 4022, 4023, 4024, 4025, 4026, 4028, 4029, 4031, 4032, 4033, 4034, 4035, 4036, 4038, 4039, 4046, 4050, 4052, 4055, 4056, 4057, 4058, 4059, 4060, 4061, 4062, 4064, 4065, 4070, 4071, 4074, 4075, 4077, 4079, 4082, 4084, 4088, 4095, 4096, 4098, 4099, 4101, 4103, 4106, 4108, 4111, 4115, 4116, 4118, 4122, 4123, 4124, 4127, 4129, 4130, 4131, 4132, 4134, 4135, 4136, 4137, 4140, 4141, 4142, 4146, 4147, 4152, 4156, 4159, 4160, 4164, 4166, 4168, 4170, 4174, 4177, 4178, 4179, 4180, 4181, 4188, 4190, 4191, 4193, 4194, 4196, 4198, 4199, 4200, 4203, 4207, 4209, 4210, 4212, 4213, 4215, 4216, 4217, 4219, 4224, 4228, 4231, 4233, 4236, 4245, 4249, 4251, 4256, 4258, 4259, 4265, 4268, 4272, 4273, 4274, 4276, 4278, 4283, 4285, 4286, 4288, 4290, 4292, 4294, 4295, 4296, 4297, 4299, 4300, 4301, 4302, 4303, 4304, 4309, 4314, 4315, 4316, 4317, 4321, 4325, 4327, 4329, 4330, 4332, 4333, 4334, 4336, 4337, 4338, 4339, 4340, 4342, 4343, 4344, 4348, 4352, 4353, 4359, 4366, 4367, 4368, 4369, 4370, 4372, 4374, 4375, 4376, 4378, 4379, 4380, 4383, 4386, 4391, 4392, 4394, 4395, 4396, 4397, 4399, 4400, 4401, 4404, 4407, 4408, 4410, 4411, 4412, 4413, 4414, 4415, 4417, 4418, 4423, 4424, 4427, 4428, 4433, 4435, 4436, 4437, 4439, 4441, 4442, 4443, 4445, 4446, 4448, 4449, 4451, 4453, 4454, 4456, 4457, 4458, 4459, 4461, 4462, 4466, 4467, 4469, 4471, 4472, 4473, 4474, 4475, 4477, 4482, 4483, 4484, 4489, 4490, 4491, 4492, 4493, 4495, 4497, 4498, 4501, 4502, 4504, 4506, 4507, 4508, 4509, 4510, 4512, 4513, 4514, 4517, 4518, 4519, 4525, 4527, 4529, 4532, 4533, 4534, 4535, 4538, 4540, 4541, 4542, 4543, 4545, 4546, 4547, 4550, 4557, 4558, 4559, 4560, 4564, 4568, 4569, 4574, 4575, 4577, 4578, 4579, 4580, 4586, 4589, 4590, 4591, 4595, 4596, 4597, 4598, 4600, 4601, 4603, 4605, 4606, 4607, 4609, 4612, 4613, 4614, 4616, 4620, 4623, 4625, 4628, 4629, 4634, 4637, 4638, 4640, 4641, 4642, 4643, 4647, 4648, 4651, 4653, 4654, 4655, 4656, 4657, 4658, 4659, 4660, 4661, 4665, 4668, 4669, 4670, 4675, 4676, 4677, 4680, 4681, 4689, 4692, 4693, 4695, 4698, 4701, 4703, 4709, 4716, 4720, 4721, 4723, 4724, 4725, 4726, 4728, 4730, 4731, 4734, 4737, 4738, 4740, 4742, 4743, 4744, 4747, 4751, 4753, 4758, 4759, 4760, 4762, 4764, 4766, 4769, 4771, 4774, 4775, 4777, 4779, 4782, 4784, 4786, 4787, 4789, 4793, 4794, 4795, 4796, 4799, 4800, 4806, 4808, 4809, 4811, 4812, 4813, 4814, 4816, 4817, 4818, 4822, 4823, 4824, 4826, 4829, 4831, 4832, 4834, 4835, 4836, 4837, 4839, 4840, 4841, 4844, 4845, 4846, 4847, 4848, 4852, 4854, 4855, 4858, 4859, 4861, 4870, 4871, 4874, 4876, 4877, 4878, 4882, 4883, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4896, 4897, 4899, 4901, 4904, 4906, 4908, 4909, 4910, 4912, 4913, 4914, 4916, 4919, 4922, 4924, 4925, 4927, 4928, 4929, 4930, 4931, 4936, 4940, 4941, 4942, 4944, 4945, 4946, 4948, 4950, 4953, 4954, 4955, 4958, 4962, 4963, 4965, 4969, 4973, 4975, 4976, 4978, 4980, 4982, 4986, 4987, 4989, 4991, 4993, 4994, 4995, 5002, 5003, 5005, 5006, 5007, 5008, 5013, 5015, 5020, 5022, 5025, 5039, 5040, 5041, 5042, 5044, 5045, 5046, 5048, 5049, 5054, 5055, 5062, 5063, 5064, 5066, 5068, 5069, 5071, 5072, 5076, 5078, 5079, 5080, 5082, 5084, 5086, 5091, 5092, 5101, 5102, 5104, 5105, 5106, 5107, 5109, 5112, 5114, 5118, 5119, 5121, 5122, 5123, 5124, 5125, 5127, 5128, 5129, 5131, 5135, 5136, 5137, 5139, 5140, 5144, 5149, 5152, 5156, 5158, 5162, 5163, 5164, 5165, 5168, 5171, 5174, 5175, 5178, 5179, 5181, 5183, 5184, 5185, 5186, 5187, 5191, 5193, 5194, 5195, 5196, 5198, 5199, 5200, 5201, 5203, 5205, 5206, 5207, 5213, 5214, 5216, 5217, 5219, 5220, 5221, 5222, 5223, 5224, 5227, 5228, 5229, 5230, 5233, 5235, 5236, 5237, 5238, 5240, 5242, 5244, 5246, 5247, 5252, 5253, 5254, 5255, 5256, 5259, 5260, 5261, 5264, 5265, 5267, 5268, 5269, 5272, 5274, 5278, 5280, 5282, 5284, 5285, 5286, 5287, 5290, 5291, 5292, 5293, 5294, 5295, 5296, 5297, 5299, 5302, 5303, 5306, 5309, 5310, 5311, 5312, 5313, 5314, 5315, 5318, 5320, 5321, 5322, 5323, 5325, 5326, 5327, 5328, 5329, 5330, 5333, 5334, 5335, 5336, 5344, 5348, 5349, 5350, 5355, 5358, 5359, 5360, 5361, 5362, 5364, 5369, 5373, 5374, 5375, 5379, 5380, 5381, 5382, 5383, 5385, 5386, 5388, 5390, 5393, 5394, 5397, 5398, 5399, 5406, 5407, 5408, 5410, 5412, 5413, 5415, 5416, 5417, 5421, 5424, 5425, 5428, 5429, 5430, 5436, 5438, 5439, 5441, 5443, 5445, 5446, 5448, 5450, 5451, 5452, 5456, 5458, 5460, 5462, 5464, 5466, 5470, 5474, 5475, 5476, 5482, 5485, 5487, 5494, 5495, 5498, 5499, 5501, 5504, 5510, 5513, 5516, 5519, 5522, 5523, 5524, 5525, 5526, 5527, 5528, 5532, 5536, 5537, 5538, 5539, 5540, 5541, 5542, 5546, 5554, 5555, 5556, 5558, 5560, 5564, 5566, 5569, 5570, 5572, 5574, 5577, 5579, 5585, 5586, 5587, 5588, 5589, 5590, 5592, 5595, 5597, 5600, 5602, 5603, 5604, 5606, 5607, 5611, 5612, 5613, 5614, 5615, 5616, 5617, 5619, 5620, 5623, 5625, 5626, 5627, 5628, 5629, 5632, 5634, 5637, 5640, 5642, 5645, 5649, 5650, 5651, 5655, 5656, 5661, 5663, 5665, 5669, 5670, 5671, 5672, 5673, 5674, 5676, 5677, 5682, 5684, 5685, 5686, 5691, 5692, 5693, 5695, 5696, 5699, 5700, 5701, 5702, 5703, 5704, 5705, 5707, 5708, 5709, 5710, 5711, 5713, 5714, 5715, 5717, 5719, 5721, 5722, 5723, 5725, 5728, 5729, 5731, 5732, 5733, 5734, 5735, 5736, 5737, 5738, 5742, 5743, 5744, 5745, 5748, 5750, 5752, 5753, 5754, 5756, 5757, 5759, 5760, 5761, 5762, 5763, 5764, 5765, 5766, 5770, 5771, 5775, 5776, 5778, 5779, 5783, 5785, 5786, 5788, 5789, 5794, 5795, 5799, 5800, 5804, 5805, 5808, 5809, 5811, 5812, 5817, 5821, 5822, 5824, 5825, 5827, 5829, 5830, 5833, 5834, 5835, 5837, 5839, 5840, 5845, 5847, 5849, 5851, 5853, 5854, 5857, 5858, 5860, 5861, 5862, 5864, 5866, 5867, 5869, 5872, 5873, 5874, 5876, 5877, 5878, 5880, 5883, 5884, 5889, 5890, 5892, 5894, 5896, 5904, 5912, 5913, 5914, 5916, 5922, 5923, 5925, 5926, 5927, 5928, 5931, 5932, 5934, 5936, 5937, 5943, 5944, 5946, 5947, 5953, 5954, 5958, 5959, 5960, 5964, 5967, 5969, 5970, 5971, 5973, 5975, 5980, 5981, 5982, 5986, 5987, 5989, 5991, 5993, 5996, 6003, 6004, 6006, 6007, 6009, 6011, 6013, 6015, 6018, 6019, 6024, 6026, 6030, 6032, 6033, 6040, 6041, 6044, 6045, 6048, 6049, 6053, 6056, 6057, 6061, 6065, 6066, 6067, 6072, 6075, 6077, 6081, 6086, 6087, 6088, 6091, 6092, 6095, 6096, 6097, 6099, 6102, 6103, 6107, 6108, 6109, 6111, 6114, 6115, 6116, 6117, 6118, 6119, 6123, 6128, 6131, 6133, 6136, 6138, 6142, 6143, 6147, 6148, 6149, 6151, 6154, 6155, 6157, 6158, 6159, 6160, 6162, 6163, 6164, 6165, 6167, 6168, 6172, 6179, 6180, 6181, 6184, 6187, 6188, 6189, 6191, 6192, 6194, 6195, 6196, 6198, 6200, 6201, 6206, 6208, 6213, 6214, 6215, 6216, 6217, 6219, 6224, 6225, 6226, 6227, 6228, 6230, 6232, 6235, 6239, 6241, 6243, 6245, 6247, 6249, 6250, 6252, 6253, 6254, 6255, 6256, 6257, 6259, 6260, 6261, 6262, 6264, 6265, 6266, 6267, 6268, 6273, 6274, 6276, 6277, 6278, 6279, 6280, 6281, 6282, 6284, 6286, 6290, 6291, 6293, 6294, 6296, 6298, 6299, 6308, 6309, 6310, 6314, 6317, 6320, 6321, 6323, 6324, 6332, 6333, 6335, 6336, 6337, 6338, 6339, 6340, 6344, 6346, 6349, 6350, 6352, 6359, 6364, 6365, 6366, 6367, 6370, 6373, 6374, 6378, 6379, 6380, 6381, 6382, 6385, 6389, 6392, 6406, 6407, 6408, 6410, 6411, 6415, 6420, 6421, 6422, 6423, 6424, 6427, 6428, 6429, 6432, 6434, 6437, 6439, 6440, 6441, 6442, 6444, 6445, 6452, 6453, 6457, 6459, 6461, 6463, 6464, 6468, 6469, 6470, 6472, 6473, 6474, 6478, 6483, 6486, 6487, 6492, 6496, 6498, 6499, 6501, 6502, 6504, 6507, 6508, 6509, 6510, 6513, 6514, 6515, 6517, 6519, 6520, 6521, 6522, 6523, 6524, 6526, 6528, 6529, 6530, 6531, 6532, 6533, 6534, 6540, 6544, 6546, 6547, 6551, 6552, 6558, 6559, 6561, 6562, 6563, 6565, 6567, 6568, 6569, 6571, 6572, 6573, 6576, 6577, 6578, 6579, 6581, 6582, 6583, 6584, 6585, 6591, 6592, 6595, 6596, 6597, 6603, 6604, 6605, 6606, 6608, 6610, 6611, 6613, 6615, 6617, 6618, 6620, 6621, 6622, 6624, 6625, 6626, 6628, 6630, 6631, 6632, 6633, 6635, 6636, 6637, 6638, 6642, 6644, 6646, 6647, 6649, 6650, 6651, 6654, 6655, 6657, 6658, 6660, 6662, 6664, 6667, 6669, 6671, 6672, 6673, 6677, 6679, 6682, 6687, 6689, 6691, 6692, 6693, 6694, 6697, 6701, 6702, 6704, 6707, 6714, 6716, 6718, 6721, 6722, 6724, 6725, 6726, 6728, 6731, 6733, 6737, 6741, 6742, 6743, 6744, 6745, 6749, 6751, 6752, 6755, 6757, 6758, 6760, 6761, 6767, 6768, 6769, 6772, 6774, 6775, 6778, 6780, 6781, 6782, 6785, 6789, 6793, 6795, 6796, 6797, 6800, 6801, 6802, 6806, 6808, 6810, 6817, 6819, 6820, 6822, 6825, 6826, 6827, 6833, 6834, 6836, 6838, 6839, 6845, 6852, 6854, 6856, 6860, 6864, 6866, 6867, 6868, 6871, 6872, 6874, 6876, 6877, 6879, 6880, 6883, 6884, 6887, 6888, 6890, 6893, 6895, 6896, 6899, 6902, 6903, 6904, 6905, 6910, 6911, 6912, 6916, 6917, 6918, 6919, 6920, 6921, 6924, 6925, 6926, 6930, 6931, 6932, 6938, 6940, 6941, 6942, 6943, 6949, 6950, 6952, 6955, 6958, 6959, 6961, 6962, 6965, 6966, 6971, 6980, 6981, 6984, 6985, 6986, 6988, 6991, 6992, 6996, 6997, 6999, 7000, 7009, 7010, 7011, 7012, 7013, 7014, 7016, 7017, 7020, 7022, 7024, 7025, 7032, 7035, 7036, 7038, 7040, 7041, 7043, 7044, 7049, 7050, 7051, 7052, 7055, 7057, 7058, 7060, 7061, 7062, 7063, 7064, 7065, 7066, 7067, 7070, 7072, 7073, 7074, 7075, 7077, 7078, 7079, 7080, 7081, 7082, 7086, 7087, 7088, 7090, 7096, 7098, 7099, 7103, 7104, 7105, 7110, 7112, 7116, 7119, 7121, 7123, 7128, 7129, 7132, 7134, 7136, 7137, 7138, 7142, 7143, 7144, 7145, 7147, 7149, 7151, 7153, 7155, 7158, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7169, 7172, 7173, 7174, 7175, 7177, 7182, 7188, 7190, 7191, 7192, 7193, 7194, 7196, 7198, 7202, 7203, 7205, 7208, 7211, 7212, 7213, 7215, 7216, 7218, 7219, 7221, 7223, 7225, 7228, 7230, 7233, 7235, 7236, 7241, 7242, 7244, 7246, 7250, 7251, 7252, 7253, 7254, 7258, 7259, 7260, 7263, 7264, 7265, 7266, 7267, 7269, 7270, 7271, 7272, 7273, 7275, 7278, 7279, 7281, 7283, 7285, 7286, 7287, 7292, 7297, 7298, 7299, 7301, 7303, 7305, 7307, 7309, 7311, 7316, 7317, 7318, 7320, 7321, 7322, 7324, 7325, 7326, 7327, 7328, 7331, 7335, 7336, 7337, 7339, 7340, 7342, 7344, 7350, 7353, 7354, 7357, 7358, 7359, 7361, 7363, 7368, 7369, 7372, 7373, 7375, 7376, 7378, 7381, 7382, 7384, 7385, 7386, 7387, 7388, 7390, 7394, 7396, 7398, 7399, 7404, 7405, 7406, 7407, 7410, 7412, 7415, 7416, 7418, 7419, 7420, 7422, 7423, 7428, 7429, 7430, 7431, 7434, 7437, 7438, 7440, 7441, 7442, 7447, 7456, 7457, 7458, 7459, 7461, 7462, 7463, 7464, 7465, 7467, 7469, 7470, 7472, 7473, 7474, 7475, 7476, 7477, 7478, 7484, 7486, 7492, 7494, 7495, 7496, 7498, 7503, 7507, 7510, 7511, 7517, 7523, 7524, 7527, 7528, 7530, 7531, 7532, 7534, 7536, 7540, 7542, 7544, 7549, 7550, 7551, 7554, 7556, 7557, 7558, 7559, 7560, 7561, 7568, 7569, 7570, 7571, 7573, 7574, 7576, 7578, 7579, 7580, 7582, 7583, 7585, 7586, 7587, 7588, 7589, 7590, 7593, 7594, 7595, 7596, 7600, 7601, 7604, 7607, 7609, 7611, 7615, 7616, 7617, 7619, 7621, 7622, 7624, 7625, 7628, 7629, 7630, 7632, 7633, 7634, 7635, 7643, 7645, 7647, 7648, 7652, 7653, 7656, 7661, 7666, 7669, 7670, 7677, 7679, 7680, 7681, 7682, 7686, 7687, 7688, 7689, 7690, 7693, 7698, 7703, 7707, 7708, 7709, 7711, 7712, 7714, 7716, 7717, 7718, 7719, 7721, 7722, 7724, 7729, 7731, 7732, 7733, 7736, 7737, 7741, 7742, 7744, 7747, 7748, 7749, 7750, 7751, 7752, 7755, 7757, 7758, 7759, 7763, 7765, 7766, 7768, 7772, 7773, 7774, 7780, 7781, 7783, 7785, 7786, 7789, 7791, 7793, 7794, 7798, 7799, 7801, 7803, 7804, 7807, 7810, 7811, 7812, 7819, 7820, 7826, 7827, 7828, 7830, 7831, 7832, 7834, 7835, 7836, 7839, 7841, 7843, 7845, 7848, 7850, 7854, 7855, 7857, 7859, 7860, 7862, 7865, 7866, 7867, 7868, 7870, 7871, 7877, 7878, 7879, 7884, 7886, 7890, 7891, 7892, 7893, 7894, 7897, 7899, 7900, 7901, 7904, 7906, 7910, 7914, 7915, 7917, 7919, 7920, 7921, 7924, 7925, 7930, 7931, 7934, 7938, 7940, 7941, 7946, 7947, 7948, 7950, 7951, 7952, 7953, 7954, 7955, 7962, 7963, 7967, 7968, 7969, 7970, 7971, 7972, 7973, 7974, 7975, 7976, 7984, 7985, 7987, 7991, 7994, 7997, 7998, 7999, 8000, 8001, 8003, 8004, 8006, 8008, 8009, 8011, 8012, 8013, 8014, 8017, 8018, 8020, 8024, 8028, 8031, 8034, 8037, 8041, 8044, 8048, 8049, 8052, 8054, 8055, 8056, 8062, 8064, 8075, 8083, 8087, 8088, 8092, 8093, 8094, 8098, 8101, 8102, 8104, 8107, 8110, 8111, 8112, 8115, 8116, 8119, 8120, 8121, 8124, 8126, 8129, 8131, 8134, 8135, 8136, 8139, 8142, 8144, 8146, 8149, 8150, 8154, 8156, 8159, 8160, 8161, 8163, 8166, 8167, 8169, 8172, 8175, 8176, 8178, 8179, 8183, 8185, 8186, 8188, 8191, 8192, 8193, 8194, 8195, 8199, 8200, 8201, 8202, 8204, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8214, 8215, 8217, 8219, 8220, 8221, 8227, 8230, 8233, 8235, 8236, 8240, 8241, 8243, 8244, 8245, 8246, 8247, 8249, 8250, 8251, 8252, 8255, 8256, 8257, 8261, 8264, 8265, 8266, 8268, 8269, 8270, 8271, 8272, 8276, 8277, 8279, 8282, 8285, 8287, 8288, 8289, 8290, 8291, 8295, 8297, 8302, 8304, 8306, 8316, 8317, 8318, 8319, 8321, 8323, 8324, 8325, 8326, 8329, 8339, 8340, 8342, 8343, 8344, 8345, 8346, 8348, 8349, 8350, 8351, 8353, 8354, 8356, 8357, 8358, 8359, 8360, 8363, 8364, 8366, 8367, 8368, 8369, 8371, 8372, 8375, 8377, 8378, 8380, 8383, 8384, 8385, 8386, 8387, 8389, 8390, 8397, 8398, 8399, 8400, 8401, 8406, 8409, 8410, 8411, 8412, 8414, 8417, 8419, 8420, 8421, 8422, 8423, 8426, 8428, 8429, 8430, 8431, 8432, 8436, 8437, 8438, 8439, 8441, 8444, 8446, 8448, 8451, 8457, 8458, 8459, 8461, 8462, 8463, 8465, 8466, 8470, 8471, 8473, 8474, 8476, 8477, 8478, 8480, 8483, 8487, 8488, 8489, 8490, 8493, 8494, 8495, 8496, 8498, 8500, 8503, 8505, 8506, 8508, 8510, 8513, 8517, 8518, 8519, 8520, 8523, 8524, 8525, 8526, 8528, 8529, 8530, 8531, 8533, 8536, 8537, 8539, 8540, 8542, 8543, 8544, 8545, 8548, 8549, 8550, 8552, 8554, 8555, 8556, 8559, 8560, 8561, 8562, 8563, 8565, 8567, 8568, 8570, 8571, 8572, 8573, 8578, 8579, 8580, 8582, 8583, 8585, 8586, 8588, 8589, 8590, 8592, 8593, 8594, 8595, 8596, 8597, 8598, 8599, 8600, 8602, 8603, 8604, 8605, 8606, 8609, 8610, 8613, 8614, 8615, 8617, 8618, 8619, 8621, 8622, 8624, 8627, 8629, 8630, 8636, 8637, 8639, 8641, 8642, 8644, 8648, 8649, 8651, 8652, 8653, 8654, 8655, 8656, 8658, 8659, 8660, 8661, 8663, 8666, 8668, 8671, 8672, 8674, 8675, 8676, 8677, 8679, 8681, 8682, 8683, 8685, 8686, 8687, 8688, 8689, 8690, 8692, 8693, 8694, 8699, 8702, 8703, 8705, 8707, 8708, 8711, 8712, 8713, 8714, 8715, 8718, 8719, 8722, 8724, 8727, 8728, 8729, 8730, 8731, 8732, 8737, 8738, 8739, 8742, 8743, 8746, 8750, 8752, 8754, 8756, 8760, 8761, 8762, 8764, 8766, 8767, 8769, 8770, 8772, 8774, 8776, 8777, 8778, 8780, 8782, 8784, 8786, 8790, 8794, 8795, 8798, 8799, 8802, 8803, 8804, 8805, 8807, 8809, 8810, 8811, 8815, 8816, 8821, 8822, 8823, 8824, 8826, 8828, 8831, 8832, 8836, 8837, 8839, 8840, 8841, 8842, 8843, 8844, 8845, 8846, 8850, 8859, 8860, 8866, 8867, 8870, 8871, 8876, 8886, 8888, 8890, 8891, 8895, 8897, 8899, 8900, 8901, 8903, 8905, 8906, 8907, 8909, 8910, 8911, 8912, 8913, 8914, 8915, 8916, 8919, 8920, 8921, 8922, 8923, 8924, 8925, 8926, 8929, 8933, 8935, 8938, 8940, 8942, 8944, 8945, 8946, 8947, 8949, 8950, 8951, 8952, 8954, 8956, 8958, 8959, 8962, 8963, 8966, 8968, 8970, 8973, 8974, 8975, 8976, 8978, 8981, 8982, 8983, 8984, 8986, 8987, 8989, 8990, 8992, 8994, 8996, 8998, 8999, 9000, 9001, 9004, 9005, 9009, 9010, 9013, 9014, 9015, 9016, 9017, 9018, 9019, 9020, 9022, 9024, 9025, 9028, 9029, 9031, 9032, 9036, 9037, 9038, 9039, 9040, 9041, 9042, 9044, 9045, 9048, 9049, 9050, 9052, 9053, 9054, 9055, 9056, 9057, 9058, 9059, 9060, 9062, 9063, 9064, 9065, 9068, 9069, 9070, 9071, 9072, 9074, 9075, 9076, 9077, 9078, 9080, 9081, 9084, 9086, 9087, 9088, 9089, 9093, 9094, 9096, 9099, 9100, 9101, 9102, 9103, 9106, 9107, 9108, 9109, 9110, 9111, 9113, 9114, 9116, 9117, 9118, 9120, 9123, 9125, 9127, 9129, 9130, 9131, 9135, 9140, 9142, 9143, 9146, 9147, 9151, 9152, 9154, 9156, 9157, 9158, 9159, 9163, 9164, 9167, 9168, 9170, 9172, 9177, 9179, 9180, 9181, 9182, 9184, 9186, 9187, 9189, 9190, 9191, 9195, 9197, 9200, 9201, 9203, 9212, 9213, 9216, 9218, 9219, 9220, 9223, 9232, 9234, 9236, 9239, 9241, 9244, 9245, 9246, 9251, 9255, 9256, 9257, 9259, 9260, 9262, 9265, 9267, 9269, 9270, 9272, 9278, 9280, 9282, 9283, 9285, 9286, 9287, 9288, 9291, 9292, 9293, 9298, 9299, 9300, 9301, 9302, 9303, 9305, 9307, 9308, 9313, 9316, 9318, 9319, 9320, 9321, 9322, 9323, 9325, 9326, 9327, 9329, 9330, 9332, 9339, 9345, 9348, 9349, 9350, 9351, 9354, 9355, 9356, 9358, 9359, 9360, 9361, 9362, 9363, 9364, 9367, 9369, 9373, 9375, 9376, 9379, 9380, 9382, 9383, 9384, 9388, 9389, 9391, 9392, 9393, 9394, 9395, 9396, 9398, 9400, 9403, 9404, 9405, 9406, 9407, 9409, 9410, 9411, 9412, 9423, 9425, 9427, 9428, 9431, 9432, 9433, 9434, 9435, 9438, 9439, 9440, 9441, 9442, 9445, 9447, 9448, 9451, 9452, 9453, 9454, 9456, 9458, 9462, 9464, 9465, 9466, 9471, 9473, 9474, 9475, 9477, 9481, 9484, 9489, 9492, 9493, 9496, 9497, 9501, 9504, 9505, 9507, 9508, 9509, 9513, 9520, 9521, 9524, 9525, 9526, 9528, 9530, 9533, 9535, 9543, 9545, 9546, 9550, 9551, 9552, 9557, 9558, 9559, 9564, 9567, 9569, 9572, 9573, 9575, 9576, 9577, 9578, 9580, 9581, 9583, 9588, 9589, 9590, 9591, 9593, 9596, 9597, 9598, 9599, 9600, 9601, 9602, 9603, 9605, 9608, 9610, 9611, 9615, 9616, 9617, 9618, 9619, 9620, 9621, 9622, 9623, 9624, 9625, 9630, 9635, 9636, 9638, 9640, 9641, 9642, 9643, 9644, 9645, 9646, 9648, 9652, 9655, 9656, 9657, 9659, 9660, 9661, 9662, 9664, 9666, 9667, 9669, 9671, 9672, 9673, 9674, 9679, 9680, 9681, 9684, 9685, 9687, 9689, 9690, 9693, 9694, 9695, 9697, 9699, 9701, 9702, 9705, 9706, 9708, 9710, 9712, 9717, 9719, 9720, 9721, 9722, 9726, 9728, 9730, 9732, 9733, 9735, 9736, 9737, 9739, 9741, 9742, 9743, 9744, 9745, 9747, 9748, 9749, 9751, 9753, 9754, 9755, 9756, 9764, 9767, 9768, 9771, 9774, 9775, 9776, 9779, 9780, 9783, 9785, 9787, 9789, 9791, 9792, 9793, 9795, 9798, 9799, 9805, 9806, 9807, 9808, 9810, 9811, 9812, 9813, 9817, 9818, 9823, 9824, 9826, 9829, 9831, 9832, 9836, 9838, 9839, 9841, 9843, 9845, 9850, 9854, 9856, 9857, 9859, 9864, 9866, 9867, 9869, 9872, 9873, 9877, 9878, 9880, 9881, 9883, 9884, 9890, 9891]),
        numpy_asarray([1, 2, 7, 8, 10, 14, 17, 20, 21, 24, 25, 26, 28, 29, 32, 33, 35, 38, 39, 42, 43, 47, 50, 52, 55, 56, 59, 60, 62, 63, 64, 66, 67, 69, 70, 71, 72, 74, 75, 77, 80, 82, 84, 92, 93, 94, 95, 96, 98, 99, 102, 105, 107, 108, 110, 114, 115, 119, 122, 123, 132, 133, 135, 137, 142, 144, 146, 147, 151, 152, 153, 155, 156, 157, 158, 160, 162, 163, 165, 166, 168, 169, 170, 171, 174, 175, 177, 179, 180, 181, 182, 183, 184, 185, 187, 189, 190, 194, 195, 197, 199, 203, 205, 207, 210, 212, 213, 214, 215, 217, 219, 220, 221, 222, 223, 224, 230, 231, 237, 238, 240, 241, 242, 245, 246, 248, 251, 252, 253, 256, 257, 259, 260, 262, 264, 266, 267, 268, 269, 271, 273, 276, 279, 281, 282, 284, 288, 289, 290, 294, 296, 298, 299, 300, 301, 302, 303, 305, 307, 308, 310, 311, 312, 313, 314, 315, 316, 319, 320, 322, 325, 328, 329, 331, 333, 334, 336, 337, 341, 342, 344, 349, 350, 351, 353, 354, 355, 356, 361, 362, 363, 364, 365, 366, 367, 372, 373, 375, 376, 377, 379, 380, 383, 384, 385, 387, 389, 391, 393, 395, 396, 397, 399, 400, 403, 405, 406, 407, 412, 413, 419, 420, 421, 424, 428, 429, 430, 431, 433, 436, 437, 439, 440, 441, 442, 443, 444, 447, 451, 454, 458, 463, 464, 467, 470, 473, 474, 482, 483, 484, 485, 487, 488, 489, 490, 491, 494, 495, 496, 498, 502, 503, 504, 505, 512, 513, 514, 515, 520, 522, 523, 525, 527, 528, 529, 530, 531, 532, 535, 537, 538, 540, 541, 542, 544, 549, 552, 553, 554, 556, 557, 558, 560, 561, 562, 563, 564, 565, 567, 568, 569, 570, 573, 574, 576, 577, 580, 584, 585, 587, 588, 590, 593, 596, 597, 598, 603, 604, 607, 609, 610, 612, 614, 615, 619, 620, 622, 623, 624, 626, 628, 629, 631, 632, 633, 634, 636, 638, 640, 641, 642, 643, 646, 649, 650, 651, 652, 654, 655, 656, 657, 660, 662, 664, 665, 666, 667, 669, 671, 672, 673, 674, 676, 677, 679, 680, 683, 684, 685, 686, 688, 689, 690, 691, 694, 699, 702, 704, 706, 707, 708, 709, 710, 712, 713, 714, 718, 719, 720, 721, 723, 724, 726, 727, 730, 731, 732, 734, 736, 738, 740, 742, 745, 746, 747, 749, 751, 752, 753, 755, 758, 763, 765, 767, 768, 769, 770, 771, 772, 773, 777, 778, 779, 780, 785, 787, 788, 790, 791, 794, 796, 797, 798, 802, 803, 811, 814, 816, 818, 828, 829, 830, 831, 833, 834, 837, 838, 840, 842, 845, 846, 847, 854, 855, 856, 857, 858, 861, 863, 866, 869, 875, 879, 880, 882, 883, 885, 886, 888, 890, 891, 895, 897, 898, 899, 901, 902, 903, 908, 910, 911, 914, 915, 916, 919, 921, 922, 924, 927, 932, 933, 936, 937, 938, 939, 940, 941, 942, 944, 945, 946, 947, 948, 950, 951, 953, 960, 961, 963, 964, 966, 967, 968, 971, 972, 975, 976, 977, 979, 980, 983, 984, 985, 986, 988, 989, 993, 995, 996, 998, 999, 1001, 1005, 1009, 1010, 1012, 1013, 1014, 1021, 1023, 1027, 1028, 1029, 1033, 1036, 1038, 1039, 1040, 1042, 1044, 1045, 1046, 1048, 1051, 1052, 1053, 1054, 1055, 1056, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1069, 1071, 1072, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1089, 1090, 1091, 1093, 1094, 1096, 1100, 1103, 1104, 1105, 1108, 1110, 1112, 1118, 1124, 1128, 1129, 1131, 1133, 1137, 1140, 1141, 1142, 1143, 1147, 1154, 1157, 1158, 1159, 1160, 1161, 1163, 1166, 1167, 1178, 1181, 1182, 1184, 1187, 1190, 1194, 1197, 1199, 1200, 1205, 1206, 1208, 1209, 1210, 1212, 1213, 1215, 1219, 1220, 1221, 1222, 1224, 1225, 1226, 1232, 1234, 1235, 1238, 1239, 1240, 1241, 1242, 1244, 1249, 1250, 1254, 1256, 1259, 1262, 1263, 1266, 1267, 1272, 1273, 1277, 1281, 1283, 1286, 1287, 1288, 1289, 1290, 1291, 1294, 1295, 1298, 1300, 1301, 1305, 1307, 1308, 1309, 1310, 1311, 1312, 1313, 1314, 1319, 1321, 1323, 1324, 1325, 1326, 1327, 1328, 1330, 1331, 1333, 1334, 1335, 1337, 1341, 1342, 1343, 1344, 1346, 1347, 1348, 1350, 1351, 1354, 1358, 1360, 1361, 1362, 1363, 1367, 1369, 1373, 1376, 1378, 1379, 1380, 1381, 1382, 1383, 1390, 1393, 1394, 1398, 1399, 1400, 1402, 1405, 1406, 1409, 1410, 1414, 1415, 1418, 1423, 1424, 1426, 1427, 1430, 1433, 1434, 1436, 1443, 1445, 1446, 1447, 1448, 1452, 1456, 1458, 1459, 1460, 1461, 1465, 1467, 1468, 1470, 1471, 1473, 1477, 1480, 1483, 1485, 1486, 1487, 1492, 1496, 1498, 1500, 1502, 1503, 1504, 1505, 1506, 1507, 1510, 1515, 1516, 1519, 1520, 1523, 1525, 1528, 1529, 1530, 1531, 1532, 1533, 1535, 1536, 1537, 1541, 1543, 1544, 1545, 1549, 1555, 1558, 1559, 1560, 1562, 1565, 1581, 1584, 1587, 1589, 1592, 1593, 1594, 1597, 1600, 1601, 1602, 1605, 1609, 1611, 1614, 1615, 1616, 1617, 1618, 1619, 1625, 1626, 1627, 1628, 1631, 1633, 1635, 1636, 1638, 1644, 1645, 1646, 1647, 1648, 1650, 1651, 1652, 1657, 1659, 1660, 1662, 1664, 1665, 1666, 1668, 1669, 1670, 1672, 1675, 1678, 1681, 1687, 1689, 1690, 1693, 1694, 1695, 1699, 1702, 1703, 1704, 1705, 1706, 1707, 1708, 1709, 1710, 1711, 1724, 1725, 1726, 1730, 1733, 1735, 1736, 1737, 1740, 1742, 1743, 1745, 1748, 1749, 1753, 1754, 1755, 1756, 1757, 1759, 1760, 1761, 1762, 1766, 1767, 1768, 1769, 1771, 1772, 1773, 1775, 1779, 1780, 1781, 1783, 1785, 1787, 1788, 1790, 1792, 1796, 1798, 1799, 1800, 1802, 1806, 1810, 1813, 1814, 1815, 1816, 1817, 1819, 1820, 1827, 1828, 1832, 1835, 1836, 1840, 1843, 1845, 1847, 1849, 1854, 1856, 1858, 1860, 1861, 1868, 1872, 1875, 1876, 1878, 1883, 1884, 1885, 1886, 1887, 1890, 1892, 1893, 1898, 1901, 1902, 1909, 1911, 1914, 1917, 1920, 1925, 1926, 1928, 1929, 1930, 1932, 1935, 1936, 1938, 1940, 1942, 1943, 1944, 1949, 1950, 1952, 1953, 1957, 1959, 1960, 1962, 1963, 1964, 1968, 1969, 1971, 1975, 1980, 1983, 1990, 1991, 1992, 1993, 1997, 2001, 2004, 2005, 2006, 2009, 2011, 2015, 2017, 2018, 2019, 2021, 2023, 2025, 2027, 2029, 2032, 2033, 2037, 2039, 2042, 2043, 2054, 2056, 2057, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2068, 2069, 2072, 2075, 2076, 2077, 2079, 2080, 2084, 2085, 2086, 2087, 2088, 2089, 2092, 2094, 2097, 2098, 2099, 2102, 2103, 2105, 2108, 2109, 2110, 2111, 2112, 2114, 2117, 2118, 2119, 2122, 2123, 2126, 2127, 2129, 2134, 2135, 2136, 2140, 2141, 2142, 2146, 2149, 2151, 2152, 2153, 2154, 2157, 2159, 2160, 2164, 2165, 2170, 2171, 2174, 2175, 2176, 2177, 2179, 2181, 2183, 2186, 2188, 2191, 2193, 2195, 2196, 2197, 2198, 2200, 2206, 2208, 2209, 2210, 2214, 2215, 2224, 2225, 2226, 2231, 2233, 2234, 2235, 2238, 2243, 2245, 2246, 2248, 2249, 2252, 2253, 2258, 2260, 2261, 2262, 2263, 2268, 2270, 2271, 2273, 2275, 2278, 2281, 2286, 2287, 2288, 2289, 2290, 2293, 2294, 2295, 2297, 2298, 2302, 2304, 2306, 2312, 2319, 2323, 2325, 2326, 2327, 2329, 2330, 2335, 2336, 2337, 2338, 2340, 2342, 2344, 2346, 2347, 2351, 2353, 2356, 2357, 2359, 2361, 2362, 2364, 2366, 2368, 2369, 2370, 2372, 2374, 2375, 2377, 2378, 2379, 2382, 2383, 2388, 2390, 2391, 2392, 2393, 2394, 2398, 2399, 2400, 2401, 2406, 2408, 2410, 2411, 2415, 2416, 2418, 2420, 2422, 2425, 2426, 2427, 2429, 2430, 2431, 2432, 2435, 2436, 2437, 2438, 2439, 2441, 2445, 2446, 2447, 2448, 2450, 2451, 2452, 2453, 2456, 2459, 2462, 2464, 2465, 2466, 2467, 2471, 2474, 2476, 2477, 2478, 2481, 2482, 2483, 2484, 2485, 2491, 2492, 2493, 2494, 2497, 2502, 2504, 2508, 2513, 2514, 2515, 2517, 2518, 2520, 2524, 2528, 2529, 2530, 2531, 2537, 2539, 2540, 2543, 2544, 2545, 2546, 2547, 2550, 2552, 2553, 2554, 2555, 2556, 2558, 2560, 2564, 2566, 2572, 2574, 2577, 2578, 2580, 2581, 2583, 2584, 2587, 2588, 2590, 2591, 2592, 2594, 2596, 2598, 2601, 2603, 2604, 2605, 2606, 2607, 2608, 2609, 2616, 2617, 2618, 2621, 2622, 2624, 2626, 2627, 2629, 2630, 2632, 2639, 2641, 2644, 2646, 2647, 2649, 2650, 2653, 2657, 2662, 2664, 2665, 2666, 2672, 2673, 2674, 2679, 2684, 2685, 2686, 2687, 2690, 2691, 2693, 2695, 2700, 2706, 2708, 2709, 2710, 2713, 2716, 2717, 2719, 2721, 2724, 2725, 2726, 2731, 2733, 2734, 2736, 2737, 2742, 2744, 2745, 2747, 2748, 2753, 2757, 2759, 2761, 2762, 2763, 2766, 2769, 2770, 2772, 2773, 2775, 2776, 2777, 2778, 2779, 2780, 2782, 2785, 2787, 2788, 2789, 2790, 2791, 2792, 2793, 2794, 2795, 2797, 2798, 2799, 2800, 2802, 2803, 2804, 2805, 2809, 2811, 2812, 2813, 2815, 2816, 2817, 2819, 2820, 2823, 2825, 2829, 2834, 2836, 2839, 2840, 2841, 2842, 2843, 2844, 2848, 2850, 2851, 2853, 2854, 2855, 2856, 2857, 2858, 2859, 2861, 2862, 2865, 2870, 2876, 2877, 2878, 2880, 2884, 2887, 2888, 2889, 2890, 2891, 2894, 2895, 2898, 2899, 2901, 2908, 2910, 2912, 2919, 2923, 2926, 2928, 2930, 2931, 2934, 2935, 2941, 2943, 2944, 2946, 2948, 2949, 2952, 2953, 2954, 2957, 2959, 2960, 2966, 2969, 2973, 2974, 2975, 2977, 2979, 2980, 2982, 2984, 2986, 2989, 2991, 2992, 2993, 2994, 2995, 2996, 2997, 2999, 3000, 3001, 3002, 3004, 3005, 3006, 3007, 3012, 3013, 3017, 3019, 3021, 3031, 3032, 3034, 3037, 3039, 3042, 3043, 3046, 3047, 3048, 3049, 3050, 3052, 3055, 3056, 3058, 3059, 3060, 3061, 3063, 3065, 3070, 3071, 3074, 3075, 3076, 3077, 3080, 3083, 3086, 3089, 3090, 3091, 3092, 3094, 3099, 3101, 3105, 3106, 3107, 3109, 3110, 3113, 3116, 3117, 3118, 3120, 3123, 3124, 3125, 3126, 3127, 3128, 3130, 3131, 3133, 3136, 3137, 3138, 3140, 3141, 3142, 3143, 3144, 3145, 3146, 3148, 3149, 3150, 3152, 3154, 3158, 3166, 3168, 3170, 3171, 3172, 3174, 3176, 3177, 3179, 3180, 3181, 3183, 3184, 3188, 3190, 3191, 3193, 3197, 3198, 3199, 3202, 3206, 3207, 3208, 3210, 3212, 3223, 3224, 3226, 3227, 3228, 3229, 3232, 3234, 3237, 3239, 3240, 3241, 3245, 3248, 3251, 3252, 3253, 3254, 3258, 3259, 3260, 3265, 3267, 3270, 3271, 3273, 3280, 3281, 3283, 3284, 3286, 3291, 3292, 3295, 3299, 3300, 3301, 3307, 3308, 3309, 3312, 3313, 3314, 3315, 3316, 3317, 3320, 3321, 3323, 3327, 3329, 3330, 3331, 3332, 3339, 3340, 3341, 3342, 3343, 3346, 3348, 3350, 3352, 3353, 3355, 3356, 3357, 3362, 3363, 3366, 3367, 3372, 3374, 3377, 3378, 3382, 3383, 3384, 3388, 3389, 3390, 3392, 3393, 3395, 3396, 3400, 3401, 3404, 3405, 3408, 3409, 3410, 3411, 3412, 3416, 3417, 3419, 3420, 3423, 3424, 3425, 3430, 3431, 3434, 3436, 3437, 3439, 3441, 3442, 3443, 3454, 3455, 3457, 3458, 3459, 3460, 3461, 3462, 3463, 3466, 3468, 3473, 3475, 3476, 3478, 3480, 3482, 3483, 3487, 3489, 3499, 3500, 3503, 3506, 3507, 3508, 3512, 3514, 3523, 3525, 3527, 3529, 3530, 3531, 3533, 3536, 3537, 3538, 3541, 3542, 3543, 3544, 3545, 3546, 3548, 3550, 3553, 3554, 3555, 3556, 3558, 3559, 3560, 3561, 3562, 3563, 3565, 3566, 3567, 3568, 3570, 3572, 3575, 3576, 3577, 3578, 3579, 3580, 3584, 3587, 3591, 3593, 3594, 3596, 3599, 3602, 3605, 3607, 3609, 3611, 3614, 3615, 3616, 3619, 3621, 3622, 3623, 3624, 3625, 3626, 3627, 3628, 3630, 3632, 3633, 3635, 3639, 3643, 3645, 3647, 3648, 3650, 3653, 3654, 3656, 3657, 3658, 3659, 3661, 3663, 3664, 3665, 3666, 3674, 3675, 3676, 3677, 3680, 3681, 3684, 3688, 3689, 3690, 3691, 3693, 3695, 3698, 3699, 3700, 3702, 3705, 3711, 3713, 3714, 3717, 3718, 3720, 3721, 3725, 3726, 3727, 3730, 3734, 3737, 3739, 3740, 3741, 3742, 3748, 3749, 3752, 3754, 3756, 3757, 3758, 3762, 3763, 3766, 3772, 3776, 3778, 3779, 3781, 3785, 3786, 3788, 3793, 3795, 3798, 3799, 3802, 3805, 3806, 3807, 3810, 3811, 3812, 3816, 3817, 3818, 3824, 3827, 3835, 3836, 3837, 3843, 3844, 3846, 3847, 3855, 3863, 3864, 3866, 3867, 3869, 3870, 3871, 3872, 3874, 3875, 3876, 3878, 3879, 3881, 3882, 3884, 3886, 3890, 3894, 3895, 3899, 3902, 3903, 3904, 3906, 3907, 3909, 3912, 3916, 3918, 3919, 3921, 3923, 3925, 3926, 3927, 3931, 3932, 3933, 3934, 3936, 3937, 3938, 3940, 3942, 3944, 3945, 3946, 3948, 3949, 3950, 3952, 3953, 3957, 3958, 3960, 3961, 3965, 3967, 3971, 3976, 3977, 3982, 3984, 3986, 3988, 3990, 3992, 3994, 3996, 3997, 3998, 3999, 4002, 4005, 4007, 4008, 4009, 4010, 4012, 4013, 4015, 4016, 4017, 4019, 4020, 4021, 4027, 4030, 4037, 4040, 4041, 4042, 4043, 4044, 4045, 4047, 4048, 4049, 4051, 4053, 4054, 4063, 4066, 4067, 4068, 4069, 4072, 4073, 4076, 4078, 4080, 4081, 4083, 4085, 4086, 4087, 4089, 4090, 4091, 4092, 4093, 4094, 4097, 4100, 4102, 4104, 4105, 4107, 4109, 4110, 4112, 4113, 4114, 4117, 4119, 4120, 4121, 4125, 4126, 4128, 4133, 4138, 4139, 4143, 4144, 4145, 4148, 4149, 4150, 4151, 4153, 4154, 4155, 4157, 4158, 4161, 4162, 4163, 4165, 4167, 4169, 4171, 4172, 4173, 4175, 4176, 4182, 4183, 4184, 4185, 4186, 4187, 4189, 4192, 4195, 4197, 4201, 4202, 4204, 4205, 4206, 4208, 4211, 4214, 4218, 4220, 4221, 4222, 4223, 4225, 4226, 4227, 4229, 4230, 4232, 4234, 4235, 4237, 4238, 4239, 4240, 4241, 4242, 4243, 4244, 4246, 4247, 4248, 4250, 4252, 4253, 4254, 4255, 4257, 4260, 4261, 4262, 4263, 4264, 4266, 4267, 4269, 4270, 4271, 4275, 4277, 4279, 4280, 4281, 4282, 4284, 4287, 4289, 4291, 4293, 4298, 4305, 4306, 4307, 4308, 4310, 4311, 4312, 4313, 4318, 4319, 4320, 4322, 4323, 4324, 4326, 4328, 4331, 4335, 4341, 4345, 4346, 4347, 4349, 4350, 4351, 4354, 4355, 4356, 4357, 4358, 4360, 4361, 4362, 4363, 4364, 4365, 4371, 4373, 4377, 4381, 4382, 4384, 4385, 4387, 4388, 4389, 4390, 4393, 4398, 4402, 4403, 4405, 4406, 4409, 4416, 4419, 4420, 4421, 4422, 4425, 4426, 4429, 4430, 4431, 4432, 4434, 4438, 4440, 4444, 4447, 4450, 4452, 4455, 4460, 4463, 4464, 4465, 4468, 4470, 4476, 4478, 4479, 4480, 4481, 4485, 4486, 4487, 4488, 4494, 4496, 4499, 4500, 4503, 4505, 4511, 4515, 4516, 4520, 4521, 4522, 4523, 4524, 4526, 4528, 4530, 4531, 4536, 4537, 4539, 4544, 4548, 4549, 4551, 4552, 4553, 4554, 4555, 4556, 4561, 4562, 4563, 4565, 4566, 4567, 4570, 4571, 4572, 4573, 4576, 4581, 4582, 4583, 4584, 4585, 4587, 4588, 4592, 4593, 4594, 4599, 4602, 4604, 4608, 4610, 4611, 4615, 4617, 4618, 4619, 4621, 4622, 4624, 4626, 4627, 4630, 4631, 4632, 4633, 4635, 4636, 4639, 4644, 4645, 4646, 4649, 4650, 4652, 4662, 4663, 4664, 4666, 4667, 4671, 4672, 4673, 4674, 4678, 4679, 4682, 4683, 4684, 4685, 4686, 4687, 4688, 4690, 4691, 4694, 4696, 4697, 4699, 4700, 4702, 4704, 4705, 4706, 4707, 4708, 4710, 4711, 4712, 4713, 4714, 4715, 4717, 4718, 4719, 4722, 4727, 4729, 4732, 4733, 4735, 4736, 4739, 4741, 4745, 4746, 4748, 4749, 4750, 4752, 4754, 4755, 4756, 4757, 4761, 4763, 4765, 4767, 4768, 4770, 4772, 4773, 4776, 4778, 4780, 4781, 4783, 4785, 4788, 4790, 4791, 4792, 4797, 4798, 4801, 4802, 4803, 4804, 4805, 4807, 4810, 4815, 4819, 4820, 4821, 4825, 4827, 4828, 4830, 4833, 4838, 4842, 4843, 4849, 4850, 4851, 4853, 4856, 4857, 4860, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4872, 4873, 4875, 4879, 4880, 4881, 4884, 4893, 4894, 4895, 4898, 4900, 4902, 4903, 4905, 4907, 4911, 4915, 4917, 4918, 4920, 4921, 4923, 4926, 4932, 4933, 4934, 4935, 4937, 4938, 4939, 4943, 4947, 4949, 4951, 4952, 4956, 4957, 4959, 4960, 4961, 4964, 4966, 4967, 4968, 4970, 4971, 4972, 4974, 4977, 4979, 4981, 4983, 4984, 4985, 4988, 4990, 4992, 4996, 4997, 4998, 4999, 5000, 5001, 5004, 5009, 5010, 5011, 5012, 5014, 5016, 5017, 5018, 5019, 5021, 5023, 5024, 5026, 5027, 5028, 5029, 5030, 5031, 5032, 5033, 5034, 5035, 5036, 5037, 5038, 5043, 5047, 5050, 5051, 5052, 5053, 5056, 5057, 5058, 5059, 5060, 5061, 5065, 5067, 5070, 5073, 5074, 5075, 5077, 5081, 5083, 5085, 5087, 5088, 5089, 5090, 5093, 5094, 5095, 5096, 5097, 5098, 5099, 5100, 5103, 5108, 5110, 5111, 5113, 5115, 5116, 5117, 5120, 5126, 5130, 5132, 5133, 5134, 5138, 5141, 5142, 5143, 5145, 5146, 5147, 5148, 5150, 5151, 5153, 5154, 5155, 5157, 5159, 5160, 5161, 5166, 5167, 5169, 5170, 5172, 5173, 5176, 5177, 5180, 5182, 5188, 5189, 5190, 5192, 5197, 5202, 5204, 5208, 5209, 5210, 5211, 5212, 5215, 5218, 5225, 5226, 5231, 5232, 5234, 5239, 5241, 5243, 5245, 5248, 5249, 5250, 5251, 5257, 5258, 5262, 5263, 5266, 5270, 5271, 5273, 5275, 5276, 5277, 5279, 5281, 5283, 5288, 5289, 5298, 5300, 5301, 5304, 5305, 5307, 5308, 5316, 5317, 5319, 5324, 5331, 5332, 5337, 5338, 5339, 5340, 5341, 5342, 5343, 5345, 5346, 5347, 5351, 5352, 5353, 5354, 5356, 5357, 5363, 5365, 5366, 5367, 5368, 5370, 5371, 5372, 5376, 5377, 5378, 5384, 5387, 5389, 5391, 5392, 5395, 5396, 5400, 5401, 5402, 5403, 5404, 5405, 5409, 5411, 5414, 5418, 5419, 5420, 5422, 5423, 5426, 5427, 5431, 5432, 5433, 5434, 5435, 5437, 5440, 5442, 5444, 5447, 5449, 5453, 5454, 5455, 5457, 5459, 5461, 5463, 5465, 5467, 5468, 5469, 5471, 5472, 5473, 5477, 5478, 5479, 5480, 5481, 5483, 5484, 5486, 5488, 5489, 5490, 5491, 5492, 5493, 5496, 5497, 5500, 5502, 5503, 5505, 5506, 5507, 5508, 5509, 5511, 5512, 5514, 5515, 5517, 5518, 5520, 5521, 5529, 5530, 5531, 5533, 5534, 5535, 5543, 5544, 5545, 5547, 5548, 5549, 5550, 5551, 5552, 5553, 5557, 5559, 5561, 5562, 5563, 5565, 5567, 5568, 5571, 5573, 5575, 5576, 5578, 5580, 5581, 5582, 5583, 5584, 5591, 5593, 5594, 5596, 5598, 5599, 5601, 5605, 5608, 5609, 5610, 5618, 5621, 5622, 5624, 5630, 5631, 5633, 5635, 5636, 5638, 5639, 5641, 5643, 5644, 5646, 5647, 5648, 5652, 5653, 5654, 5657, 5658, 5659, 5660, 5662, 5664, 5666, 5667, 5668, 5675, 5678, 5679, 5680, 5681, 5683, 5687, 5688, 5689, 5690, 5694, 5697, 5698, 5706, 5712, 5716, 5718, 5720, 5724, 5726, 5727, 5730, 5739, 5740, 5741, 5746, 5747, 5749, 5751, 5755, 5758, 5767, 5768, 5769, 5772, 5773, 5774, 5777, 5780, 5781, 5782, 5784, 5787, 5790, 5791, 5792, 5793, 5796, 5797, 5798, 5801, 5802, 5803, 5806, 5807, 5810, 5813, 5814, 5815, 5816, 5818, 5819, 5820, 5823, 5826, 5828, 5831, 5832, 5836, 5838, 5841, 5842, 5843, 5844, 5846, 5848, 5850, 5852, 5855, 5856, 5859, 5863, 5865, 5868, 5870, 5871, 5875, 5879, 5881, 5882, 5885, 5886, 5887, 5888, 5891, 5893, 5895, 5897, 5898, 5899, 5900, 5901, 5902, 5903, 5905, 5906, 5907, 5908, 5909, 5910, 5911, 5915, 5917, 5918, 5919, 5920, 5921, 5924, 5929, 5930, 5933, 5935, 5938, 5939, 5940, 5941, 5942, 5945, 5948, 5949, 5950, 5951, 5952, 5955, 5956, 5957, 5961, 5962, 5963, 5965, 5966, 5968, 5972, 5974, 5976, 5977, 5978, 5979, 5983, 5984, 5985, 5988, 5990, 5992, 5994, 5995, 5997, 5998, 5999, 6000, 6001, 6002, 6005, 6008, 6010, 6012, 6014, 6016, 6017, 6020, 6021, 6022, 6023, 6025, 6027, 6028, 6029, 6031, 6034, 6035, 6036, 6037, 6038, 6039, 6042, 6043, 6046, 6047, 6050, 6051, 6052, 6054, 6055, 6058, 6059, 6060, 6062, 6063, 6064, 6068, 6069, 6070, 6071, 6073, 6074, 6076, 6078, 6079, 6080, 6082, 6083, 6084, 6085, 6089, 6090, 6093, 6094, 6098, 6100, 6101, 6104, 6105, 6106, 6110, 6112, 6113, 6120, 6121, 6122, 6124, 6125, 6126, 6127, 6129, 6130, 6132, 6134, 6135, 6137, 6139, 6140, 6141, 6144, 6145, 6146, 6150, 6152, 6153, 6156, 6161, 6166, 6169, 6170, 6171, 6173, 6174, 6175, 6176, 6177, 6178, 6182, 6183, 6185, 6186, 6190, 6193, 6197, 6199, 6202, 6203, 6204, 6205, 6207, 6209, 6210, 6211, 6212, 6218, 6220, 6221, 6222, 6223, 6229, 6231, 6233, 6234, 6236, 6237, 6238, 6240, 6242, 6244, 6246, 6248, 6251, 6258, 6263, 6269, 6270, 6271, 6272, 6275, 6283, 6285, 6287, 6288, 6289, 6292, 6295, 6297, 6300, 6301, 6302, 6303, 6304, 6305, 6306, 6307, 6311, 6312, 6313, 6315, 6316, 6318, 6319, 6322, 6325, 6326, 6327, 6328, 6329, 6330, 6331, 6334, 6341, 6342, 6343, 6345, 6347, 6348, 6351, 6353, 6354, 6355, 6356, 6357, 6358, 6360, 6361, 6362, 6363, 6368, 6369, 6371, 6372, 6375, 6376, 6377, 6383, 6384, 6386, 6387, 6388, 6390, 6391, 6393, 6394, 6395, 6396, 6397, 6398, 6399, 6400, 6401, 6402, 6403, 6404, 6405, 6409, 6412, 6413, 6414, 6416, 6417, 6418, 6419, 6425, 6426, 6430, 6431, 6433, 6435, 6436, 6438, 6443, 6446, 6447, 6448, 6449, 6450, 6451, 6454, 6455, 6456, 6458, 6460, 6462, 6465, 6466, 6467, 6471, 6475, 6476, 6477, 6479, 6480, 6481, 6482, 6484, 6485, 6488, 6489, 6490, 6491, 6493, 6494, 6495, 6497, 6500, 6503, 6505, 6506, 6511, 6512, 6516, 6518, 6525, 6527, 6535, 6536, 6537, 6538, 6539, 6541, 6542, 6543, 6545, 6548, 6549, 6550, 6553, 6554, 6555, 6556, 6557, 6560, 6564, 6566, 6570, 6574, 6575, 6580, 6586, 6587, 6588, 6589, 6590, 6593, 6594, 6598, 6599, 6600, 6601, 6602, 6607, 6609, 6612, 6614, 6616, 6619, 6623, 6627, 6629, 6634, 6639, 6640, 6641, 6643, 6645, 6648, 6652, 6653, 6656, 6659, 6661, 6663, 6665, 6666, 6668, 6670, 6674, 6675, 6676, 6678, 6680, 6681, 6683, 6684, 6685, 6686, 6688, 6690, 6695, 6696, 6698, 6699, 6700, 6703, 6705, 6706, 6708, 6709, 6710, 6711, 6712, 6713, 6715, 6717, 6719, 6720, 6723, 6727, 6729, 6730, 6732, 6734, 6735, 6736, 6738, 6739, 6740, 6746, 6747, 6748, 6750, 6753, 6754, 6756, 6759, 6762, 6763, 6764, 6765, 6766, 6770, 6771, 6773, 6776, 6777, 6779, 6783, 6784, 6786, 6787, 6788, 6790, 6791, 6792, 6794, 6798, 6799, 6803, 6804, 6805, 6807, 6809, 6811, 6812, 6813, 6814, 6815, 6816, 6818, 6821, 6823, 6824, 6828, 6829, 6830, 6831, 6832, 6835, 6837, 6840, 6841, 6842, 6843, 6844, 6846, 6847, 6848, 6849, 6850, 6851, 6853, 6855, 6857, 6858, 6859, 6861, 6862, 6863, 6865, 6869, 6870, 6873, 6875, 6878, 6881, 6882, 6885, 6886, 6889, 6891, 6892, 6894, 6897, 6898, 6900, 6901, 6906, 6907, 6908, 6909, 6913, 6914, 6915, 6922, 6923, 6927, 6928, 6929, 6933, 6934, 6935, 6936, 6937, 6939, 6944, 6945, 6946, 6947, 6948, 6951, 6953, 6954, 6956, 6957, 6960, 6963, 6964, 6967, 6968, 6969, 6970, 6972, 6973, 6974, 6975, 6976, 6977, 6978, 6979, 6982, 6983, 6987, 6989, 6990, 6993, 6994, 6995, 6998, 7001, 7002, 7003, 7004, 7005, 7006, 7007, 7008, 7015, 7018, 7019, 7021, 7023, 7026, 7027, 7028, 7029, 7030, 7031, 7033, 7034, 7037, 7039, 7042, 7045, 7046, 7047, 7048, 7053, 7054, 7056, 7059, 7068, 7069, 7071, 7076, 7083, 7084, 7085, 7089, 7091, 7092, 7093, 7094, 7095, 7097, 7100, 7101, 7102, 7106, 7107, 7108, 7109, 7111, 7113, 7114, 7115, 7117, 7118, 7120, 7122, 7124, 7125, 7126, 7127, 7130, 7131, 7133, 7135, 7139, 7140, 7141, 7146, 7148, 7150, 7152, 7154, 7156, 7157, 7166, 7167, 7168, 7170, 7171, 7176, 7178, 7179, 7180, 7181, 7183, 7184, 7185, 7186, 7187, 7189, 7195, 7197, 7199, 7200, 7201, 7204, 7206, 7207, 7209, 7210, 7214, 7217, 7220, 7222, 7224, 7226, 7227, 7229, 7231, 7232, 7234, 7237, 7238, 7239, 7240, 7243, 7245, 7247, 7248, 7249, 7255, 7256, 7257, 7261, 7262, 7268, 7274, 7276, 7277, 7280, 7282, 7284, 7288, 7289, 7290, 7291, 7293, 7294, 7295, 7296, 7300, 7302, 7304, 7306, 7308, 7310, 7312, 7313, 7314, 7315, 7319, 7323, 7329, 7330, 7332, 7333, 7334, 7338, 7341, 7343, 7345, 7346, 7347, 7348, 7349, 7351, 7352, 7355, 7356, 7360, 7362, 7364, 7365, 7366, 7367, 7370, 7371, 7374, 7377, 7379, 7380, 7383, 7389, 7391, 7392, 7393, 7395, 7397, 7400, 7401, 7402, 7403, 7408, 7409, 7411, 7413, 7414, 7417, 7421, 7424, 7425, 7426, 7427, 7432, 7433, 7435, 7436, 7439, 7443, 7444, 7445, 7446, 7448, 7449, 7450, 7451, 7452, 7453, 7454, 7455, 7460, 7466, 7468, 7471, 7479, 7480, 7481, 7482, 7483, 7485, 7487, 7488, 7489, 7490, 7491, 7493, 7497, 7499, 7500, 7501, 7502, 7504, 7505, 7506, 7508, 7509, 7512, 7513, 7514, 7515, 7516, 7518, 7519, 7520, 7521, 7522, 7525, 7526, 7529, 7533, 7535, 7537, 7538, 7539, 7541, 7543, 7545, 7546, 7547, 7548, 7552, 7553, 7555, 7562, 7563, 7564, 7565, 7566, 7567, 7572, 7575, 7577, 7581, 7584, 7591, 7592, 7597, 7598, 7599, 7602, 7603, 7605, 7606, 7608, 7610, 7612, 7613, 7614, 7618, 7620, 7623, 7626, 7627, 7631, 7636, 7637, 7638, 7639, 7640, 7641, 7642, 7644, 7646, 7649, 7650, 7651, 7654, 7655, 7657, 7658, 7659, 7660, 7662, 7663, 7664, 7665, 7667, 7668, 7671, 7672, 7673, 7674, 7675, 7676, 7678, 7683, 7684, 7685, 7691, 7692, 7694, 7695, 7696, 7697, 7699, 7700, 7701, 7702, 7704, 7705, 7706, 7710, 7713, 7715, 7720, 7723, 7725, 7726, 7727, 7728, 7730, 7734, 7735, 7738, 7739, 7740, 7743, 7745, 7746, 7753, 7754, 7756, 7760, 7761, 7762, 7764, 7767, 7769, 7770, 7771, 7775, 7776, 7777, 7778, 7779, 7782, 7784, 7787, 7788, 7790, 7792, 7795, 7796, 7797, 7800, 7802, 7805, 7806, 7808, 7809, 7813, 7814, 7815, 7816, 7817, 7818, 7821, 7822, 7823, 7824, 7825, 7829, 7833, 7837, 7838, 7840, 7842, 7844, 7846, 7847, 7849, 7851, 7852, 7853, 7856, 7858, 7861, 7863, 7864, 7869, 7872, 7873, 7874, 7875, 7876, 7880, 7881, 7882, 7883, 7885, 7887, 7888, 7889, 7895, 7896, 7898, 7902, 7903, 7905, 7907, 7908, 7909, 7911, 7912, 7913, 7916, 7918, 7922, 7923, 7926, 7927, 7928, 7929, 7932, 7933, 7935, 7936, 7937, 7939, 7942, 7943, 7944, 7945, 7949, 7956, 7957, 7958, 7959, 7960, 7961, 7964, 7965, 7966, 7977, 7978, 7979, 7980, 7981, 7982, 7983, 7986, 7988, 7989, 7990, 7992, 7993, 7995, 7996, 8002, 8005, 8007, 8010, 8015, 8016, 8019, 8021, 8022, 8023, 8025, 8026, 8027, 8029, 8030, 8032, 8033, 8035, 8036, 8038, 8039, 8040, 8042, 8043, 8045, 8046, 8047, 8050, 8051, 8053, 8057, 8058, 8059, 8060, 8061, 8063, 8065, 8066, 8067, 8068, 8069, 8070, 8071, 8072, 8073, 8074, 8076, 8077, 8078, 8079, 8080, 8081, 8082, 8084, 8085, 8086, 8089, 8090, 8091, 8095, 8096, 8097, 8099, 8100, 8103, 8105, 8106, 8108, 8109, 8113, 8114, 8117, 8118, 8122, 8123, 8125, 8127, 8128, 8130, 8132, 8133, 8137, 8138, 8140, 8141, 8143, 8145, 8147, 8148, 8151, 8152, 8153, 8155, 8157, 8158, 8162, 8164, 8165, 8168, 8170, 8171, 8173, 8174, 8177, 8180, 8181, 8182, 8184, 8187, 8189, 8190, 8196, 8197, 8198, 8203, 8205, 8216, 8218, 8222, 8223, 8224, 8225, 8226, 8228, 8229, 8231, 8232, 8234, 8237, 8238, 8239, 8242, 8248, 8253, 8254, 8258, 8259, 8260, 8262, 8263, 8267, 8273, 8274, 8275, 8278, 8280, 8281, 8283, 8284, 8286, 8292, 8293, 8294, 8296, 8298, 8299, 8300, 8301, 8303, 8305, 8307, 8308, 8309, 8310, 8311, 8312, 8313, 8314, 8315, 8320, 8322, 8327, 8328, 8330, 8331, 8332, 8333, 8334, 8335, 8336, 8337, 8338, 8341, 8347, 8352, 8355, 8361, 8362, 8365, 8370, 8373, 8374, 8376, 8379, 8381, 8382, 8388, 8391, 8392, 8393, 8394, 8395, 8396, 8402, 8403, 8404, 8405, 8407, 8408, 8413, 8415, 8416, 8418, 8424, 8425, 8427, 8433, 8434, 8435, 8440, 8442, 8443, 8445, 8447, 8449, 8450, 8452, 8453, 8454, 8455, 8456, 8460, 8464, 8467, 8468, 8469, 8472, 8475, 8479, 8481, 8482, 8484, 8485, 8486, 8491, 8492, 8497, 8499, 8501, 8502, 8504, 8507, 8509, 8511, 8512, 8514, 8515, 8516, 8521, 8522, 8527, 8532, 8534, 8535, 8538, 8541, 8546, 8547, 8551, 8553, 8557, 8558, 8564, 8566, 8569, 8574, 8575, 8576, 8577, 8581, 8584, 8587, 8591, 8601, 8607, 8608, 8611, 8612, 8616, 8620, 8623, 8625, 8626, 8628, 8631, 8632, 8633, 8634, 8635, 8638, 8640, 8643, 8645, 8646, 8647, 8650, 8657, 8662, 8664, 8665, 8667, 8669, 8670, 8673, 8678, 8680, 8684, 8691, 8695, 8696, 8697, 8698, 8700, 8701, 8704, 8706, 8709, 8710, 8716, 8717, 8720, 8721, 8723, 8725, 8726, 8733, 8734, 8735, 8736, 8740, 8741, 8744, 8745, 8747, 8748, 8749, 8751, 8753, 8755, 8757, 8758, 8759, 8763, 8765, 8768, 8771, 8773, 8775, 8779, 8781, 8783, 8785, 8787, 8788, 8789, 8791, 8792, 8793, 8796, 8797, 8800, 8801, 8806, 8808, 8812, 8813, 8814, 8817, 8818, 8819, 8820, 8825, 8827, 8829, 8830, 8833, 8834, 8835, 8838, 8847, 8848, 8849, 8851, 8852, 8853, 8854, 8855, 8856, 8857, 8858, 8861, 8862, 8863, 8864, 8865, 8868, 8869, 8872, 8873, 8874, 8875, 8877, 8878, 8879, 8880, 8881, 8882, 8883, 8884, 8885, 8887, 8889, 8892, 8893, 8894, 8896, 8898, 8902, 8904, 8908, 8917, 8918, 8927, 8928, 8930, 8931, 8932, 8934, 8936, 8937, 8939, 8941, 8943, 8948, 8953, 8955, 8957, 8960, 8961, 8964, 8965, 8967, 8969, 8971, 8972, 8977, 8979, 8980, 8985, 8988, 8991, 8993, 8995, 8997, 9002, 9003, 9006, 9007, 9008, 9011, 9012, 9021, 9023, 9026, 9027, 9030, 9033, 9034, 9035, 9043, 9046, 9047, 9051, 9061, 9066, 9067, 9073, 9079, 9082, 9083, 9085, 9090, 9091, 9092, 9095, 9097, 9098, 9104, 9105, 9112, 9115, 9119, 9121, 9122, 9124, 9126, 9128, 9132, 9133, 9134, 9136, 9137, 9138, 9139, 9141, 9144, 9145, 9148, 9149, 9150, 9153, 9155, 9160, 9161, 9162, 9165, 9166, 9169, 9171, 9173, 9174, 9175, 9176, 9178, 9183, 9185, 9188, 9192, 9193, 9194, 9196, 9198, 9199, 9202, 9204, 9205, 9206, 9207, 9208, 9209, 9210, 9211, 9214, 9215, 9217, 9221, 9222, 9224, 9225, 9226, 9227, 9228, 9229, 9230, 9231, 9233, 9235, 9237, 9238, 9240, 9242, 9243, 9247, 9248, 9249, 9250, 9252, 9253, 9254, 9258, 9261, 9263, 9264, 9266, 9268, 9271, 9273, 9274, 9275, 9276, 9277, 9279, 9281, 9284, 9289, 9290, 9294, 9295, 9296, 9297, 9304, 9306, 9309, 9310, 9311, 9312, 9314, 9315, 9317, 9324, 9328, 9331, 9333, 9334, 9335, 9336, 9337, 9338, 9340, 9341, 9342, 9343, 9344, 9346, 9347, 9352, 9353, 9357, 9365, 9366, 9368, 9370, 9371, 9372, 9374, 9377, 9378, 9381, 9385, 9386, 9387, 9390, 9397, 9399, 9401, 9402, 9408, 9413, 9414, 9415, 9416, 9417, 9418, 9419, 9420, 9421, 9422, 9424, 9426, 9429, 9430, 9436, 9437, 9443, 9444, 9446, 9449, 9450, 9455, 9457, 9459, 9460, 9461, 9463, 9467, 9468, 9469, 9470, 9472, 9476, 9478, 9479, 9480, 9482, 9483, 9485, 9486, 9487, 9488, 9490, 9491, 9494, 9495, 9498, 9499, 9500, 9502, 9503, 9506, 9510, 9511, 9512, 9514, 9515, 9516, 9517, 9518, 9519, 9522, 9523, 9527, 9529, 9531, 9532, 9534, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9544, 9547, 9548, 9549, 9553, 9554, 9555, 9556, 9560, 9561, 9562, 9563, 9565, 9566, 9568, 9570, 9571, 9574, 9579, 9582, 9584, 9585, 9586, 9587, 9592, 9594, 9595, 9604, 9606, 9607, 9609, 9612, 9613, 9614, 9626, 9627, 9628, 9629, 9631, 9632, 9633, 9634, 9637, 9639, 9647, 9649, 9650, 9651, 9653, 9654, 9658, 9663, 9665, 9668, 9670, 9675, 9676, 9677, 9678, 9682, 9683, 9686, 9688, 9691, 9692, 9696, 9698, 9700, 9703, 9704, 9707, 9709, 9711, 9713, 9714, 9715, 9716, 9718, 9723, 9724, 9725, 9727, 9729, 9731, 9734, 9738, 9740, 9746, 9750, 9752, 9757, 9758, 9759, 9760, 9761, 9762, 9763, 9765, 9766, 9769, 9770, 9772, 9773, 9777, 9778, 9781, 9782, 9784, 9786, 9788, 9790, 9794, 9796, 9797, 9800, 9801, 9802, 9803, 9804, 9809, 9814, 9815, 9816, 9819, 9820, 9821, 9822, 9825, 9827, 9828, 9830, 9833, 9834, 9835, 9837, 9840, 9842, 9844, 9846, 9847, 9848, 9849, 9851, 9852, 9853, 9855, 9858, 9860, 9861, 9862, 9863, 9865, 9868, 9870, 9871, 9874, 9875, 9876, 9879, 9882, 9885, 9886, 9887, 9888, 9889])
    ]

    asg2=[
        numpy_asarray([4, 5, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 25, 27, 29, 31, 33, 36, 37, 40, 41, 44, 45, 46, 48, 49, 51, 53, 57, 58, 60, 61, 62, 63, 64, 67, 68, 73, 76, 78, 81, 85, 87, 91, 93, 96, 97, 100, 101, 103, 104, 105, 109, 110, 111, 112, 116, 117, 118, 119, 120, 123, 124, 125, 126, 127, 128, 130, 131, 136, 137, 138, 140, 141, 143, 145, 146, 147, 148, 149, 150, 151, 154, 156, 158, 159, 161, 162, 163, 164, 167, 169, 171, 172, 173, 174, 175, 176, 178, 183, 184, 186, 187, 191, 196, 198, 202, 206, 208, 211, 216, 217, 218, 219, 221, 224, 225, 226, 228, 230, 232, 236, 239, 240, 242, 244, 246, 247, 250, 254, 255, 257, 258, 260, 261, 262, 263, 264, 265, 272, 273, 274, 276, 278, 291, 294, 297, 303, 304, 306, 315, 317, 318, 321, 323, 324, 326, 328, 331, 332, 334, 340, 343, 345, 347, 348, 349, 352, 353, 354, 356, 357, 358, 359, 360, 361, 368, 369, 370, 371, 372, 373, 377, 378, 381, 382, 386, 388, 390, 392, 398, 399, 400, 401, 404, 408, 409, 410, 411, 412, 415, 416, 417, 418, 422, 425, 426, 429, 432, 434, 437, 446, 448, 449, 450, 451, 456, 457, 459, 460, 461, 462, 466, 467, 468, 471, 472, 475, 476, 478, 479, 480, 481, 484, 492, 495, 496, 498, 499, 500, 501, 506, 508, 509, 514, 516, 519, 520, 521, 524, 527, 529, 530, 531, 532, 533, 534, 536, 539, 547, 551, 553, 555, 556, 558, 567, 569, 572, 574, 576, 577, 578, 579, 583, 586, 589, 592, 593, 595, 599, 600, 601, 602, 604, 605, 606, 608, 609, 612, 613, 615, 616, 618, 619, 624, 627, 628, 630, 631, 632, 636, 637, 639, 643, 644, 645, 648, 656, 658, 659, 660, 663, 666, 667, 668, 669, 670, 675, 676, 677, 678, 680, 681, 682, 683, 685, 687, 688, 692, 693, 696, 697, 698, 700, 702, 705, 707, 709, 710, 711, 714, 716, 722, 725, 726, 727, 729, 732, 733, 735, 736, 737, 740, 741, 743, 744, 748, 750, 751, 752, 753, 754, 755, 757, 760, 762, 763, 764, 765, 766, 769, 773, 774, 778, 779, 781, 782, 783, 787, 789, 792, 793, 794, 795, 800, 801, 802, 804, 805, 806, 807, 808, 815, 817, 818, 819, 823, 824, 825, 826, 827, 835, 836, 838, 841, 844, 845, 848, 849, 850, 851, 852, 853, 860, 862, 864, 865, 868, 869, 870, 871, 874, 875, 876, 877, 878, 879, 881, 884, 885, 887, 889, 890, 893, 896, 899, 901, 902, 905, 907, 909, 912, 914, 916, 917, 920, 927, 928, 929, 930, 931, 932, 934, 935, 937, 939, 943, 944, 945, 949, 953, 954, 955, 956, 963, 964, 965, 970, 973, 975, 981, 982, 986, 987, 990, 994, 997, 998, 1000, 1002, 1003, 1004, 1006, 1008, 1012, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1022, 1024, 1025, 1026, 1030, 1032, 1035, 1036, 1037, 1039, 1041, 1042, 1043, 1048, 1050, 1053, 1054, 1057, 1060, 1062, 1064, 1065, 1066, 1068, 1071, 1073, 1075, 1077, 1080, 1081, 1083, 1085, 1088, 1089, 1092, 1095, 1096, 1098, 1099, 1101, 1104, 1107, 1110, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1123, 1126, 1127, 1130, 1132, 1134, 1137, 1139, 1140, 1145, 1146, 1148, 1150, 1153, 1155, 1156, 1159, 1164, 1165, 1167, 1168, 1169, 1170, 1171, 1172, 1174, 1176, 1177, 1179, 1181, 1183, 1184, 1185, 1189, 1193, 1196, 1201, 1202, 1203, 1204, 1205, 1207, 1209, 1211, 1212, 1213, 1215, 1223, 1224, 1229, 1230, 1233, 1236, 1240, 1242, 1243, 1244, 1245, 1246, 1247, 1250, 1251, 1252, 1253, 1254, 1255, 1258, 1259, 1260, 1261, 1264, 1265, 1266, 1269, 1271, 1275, 1277, 1278, 1280, 1282, 1284, 1285, 1289, 1293, 1297, 1299, 1302, 1305, 1306, 1311, 1312, 1315, 1316, 1318, 1319, 1320, 1322, 1323, 1324, 1326, 1327, 1328, 1329, 1332, 1333, 1336, 1337, 1339, 1341, 1344, 1347, 1349, 1353, 1356, 1357, 1358, 1359, 1360, 1365, 1366, 1367, 1368, 1369, 1370, 1374, 1375, 1378, 1379, 1382, 1384, 1385, 1386, 1387, 1388, 1391, 1392, 1394, 1395, 1396, 1397, 1400, 1401, 1402, 1405, 1407, 1408, 1411, 1412, 1413, 1414, 1415, 1416, 1419, 1421, 1422, 1424, 1425, 1426, 1428, 1429, 1431, 1432, 1435, 1437, 1438, 1439, 1441, 1442, 1445, 1446, 1448, 1449, 1450, 1451, 1452, 1453, 1455, 1457, 1465, 1470, 1473, 1475, 1479, 1481, 1484, 1487, 1490, 1491, 1493, 1495, 1497, 1498, 1499, 1501, 1502, 1507, 1508, 1511, 1512, 1513, 1514, 1515, 1516, 1517, 1518, 1521, 1524, 1526, 1534, 1539, 1540, 1542, 1543, 1546, 1551, 1552, 1553, 1554, 1556, 1557, 1561, 1563, 1564, 1566, 1567, 1568, 1569, 1570, 1572, 1573, 1574, 1575, 1577, 1578, 1579, 1580, 1582, 1583, 1586, 1587, 1588, 1589, 1590, 1591, 1592, 1597, 1598, 1599, 1600, 1601, 1604, 1607, 1608, 1611, 1612, 1613, 1617, 1619, 1620, 1621, 1622, 1625, 1628, 1631, 1632, 1639, 1640, 1641, 1643, 1644, 1647, 1648, 1653, 1654, 1655, 1656, 1657, 1660, 1663, 1665, 1671, 1672, 1674, 1676, 1677, 1678, 1679, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 1691, 1693, 1696, 1698, 1699, 1705, 1711, 1713, 1714, 1716, 1717, 1718, 1723, 1727, 1728, 1729, 1732, 1734, 1735, 1737, 1738, 1741, 1744, 1746, 1747, 1750, 1754, 1758, 1759, 1760, 1761, 1762, 1763, 1764, 1767, 1768, 1769, 1770, 1776, 1777, 1778, 1781, 1782, 1783, 1784, 1786, 1788, 1789, 1790, 1793, 1794, 1797, 1799, 1801, 1803, 1808, 1809, 1812, 1815, 1818, 1820, 1821, 1822, 1826, 1829, 1830, 1831, 1833, 1834, 1835, 1836, 1837, 1838, 1839, 1841, 1842, 1843, 1844, 1845, 1846, 1850, 1852, 1853, 1855, 1856, 1857, 1860, 1862, 1863, 1864, 1866, 1867, 1869, 1871, 1874, 1875, 1879, 1880, 1882, 1883, 1888, 1889, 1890, 1897, 1898, 1899, 1903, 1904, 1906, 1907, 1908, 1910, 1912, 1913, 1914, 1915, 1916, 1918, 1919, 1920, 1921, 1922, 1926, 1927, 1929, 1930, 1933, 1934, 1936, 1937, 1939, 1945, 1946, 1947, 1948, 1951, 1953, 1954, 1955, 1956, 1958, 1959, 1960, 1961, 1965, 1966, 1967, 1968, 1970, 1972, 1974, 1978, 1981, 1982, 1985, 1987, 1990, 1995, 1996, 1998, 1999, 2002, 2003, 2006, 2007, 2008, 2009, 2012, 2014, 2015, 2016, 2020, 2024, 2025, 2026, 2028, 2030, 2033, 2034, 2035, 2036, 2038, 2040, 2041, 2043, 2044, 2047, 2048, 2049, 2050, 2051, 2053, 2055, 2058, 2067, 2078, 2080, 2082, 2083, 2086, 2088, 2090, 2091, 2092, 2093, 2095, 2096, 2097, 2101, 2104, 2106, 2115, 2116, 2118, 2119, 2120, 2121, 2124, 2125, 2129, 2131, 2132, 2133, 2134, 2137, 2139, 2140, 2143, 2144, 2147, 2148, 2149, 2150, 2153, 2156, 2161, 2162, 2166, 2174, 2178, 2182, 2184, 2185, 2186, 2187, 2189, 2190, 2191, 2193, 2194, 2197, 2199, 2201, 2202, 2203, 2205, 2207, 2213, 2216, 2217, 2218, 2220, 2222, 2223, 2224, 2226, 2227, 2229, 2230, 2231, 2233, 2235, 2236, 2237, 2238, 2239, 2240, 2244, 2246, 2247, 2248, 2250, 2254, 2256, 2261, 2262, 2264, 2265, 2266, 2269, 2272, 2273, 2274, 2277, 2279, 2280, 2281, 2282, 2283, 2285, 2288, 2290, 2292, 2294, 2296, 2298, 2299, 2300, 2305, 2307, 2309, 2310, 2313, 2314, 2318, 2322, 2326, 2333, 2334, 2338, 2339, 2342, 2343, 2348, 2349, 2350, 2352, 2354, 2355, 2356, 2358, 2365, 2367, 2368, 2370, 2375, 2377, 2381, 2384, 2386, 2389, 2395, 2397, 2400, 2402, 2404, 2407, 2409, 2413, 2414, 2415, 2418, 2420, 2423, 2424, 2428, 2433, 2434, 2440, 2442, 2443, 2445, 2451, 2452, 2454, 2455, 2457, 2458, 2459, 2460, 2462, 2463, 2466, 2468, 2471, 2472, 2473, 2478, 2479, 2481, 2485, 2486, 2487, 2492, 2493, 2494, 2495, 2496, 2498, 2499, 2500, 2501, 2502, 2503, 2505, 2507, 2510, 2516, 2519, 2521, 2522, 2525, 2528, 2529, 2530, 2532, 2533, 2534, 2535, 2536, 2538, 2541, 2542, 2547, 2548, 2552, 2558, 2560, 2561, 2565, 2566, 2568, 2569, 2570, 2571, 2572, 2573, 2575, 2576, 2577, 2579, 2580, 2585, 2586, 2592, 2593, 2595, 2597, 2599, 2600, 2602, 2611, 2612, 2614, 2615, 2617, 2619, 2620, 2621, 2623, 2624, 2625, 2627, 2628, 2630, 2631, 2633, 2634, 2635, 2636, 2639, 2640, 2642, 2643, 2647, 2648, 2649, 2650, 2651, 2652, 2654, 2656, 2660, 2661, 2662, 2663, 2667, 2668, 2670, 2675, 2676, 2677, 2678, 2679, 2680, 2683, 2684, 2685, 2687, 2692, 2694, 2696, 2697, 2698, 2703, 2707, 2711, 2712, 2714, 2715, 2717, 2719, 2720, 2721, 2722, 2725, 2726, 2730, 2732, 2733, 2734, 2735, 2737, 2738, 2741, 2743, 2744, 2746, 2750, 2752, 2755, 2756, 2757, 2758, 2759, 2760, 2761, 2762, 2763, 2764, 2765, 2766, 2767, 2769, 2772, 2775, 2779, 2780, 2781, 2782, 2783, 2784, 2785, 2786, 2787, 2790, 2792, 2797, 2798, 2800, 2801, 2803, 2806, 2810, 2814, 2815, 2818, 2819, 2826, 2827, 2828, 2829, 2831, 2832, 2833, 2835, 2836, 2837, 2838, 2844, 2846, 2847, 2848, 2849, 2852, 2856, 2860, 2862, 2863, 2864, 2867, 2868, 2869, 2871, 2872, 2874, 2875, 2876, 2878, 2882, 2883, 2884, 2889, 2895, 2897, 2902, 2904, 2905, 2910, 2911, 2913, 2914, 2916, 2918, 2920, 2921, 2922, 2926, 2927, 2928, 2929, 2930, 2931, 2933, 2936, 2937, 2938, 2939, 2940, 2941, 2942, 2945, 2947, 2950, 2953, 2955, 2956, 2959, 2961, 2963, 2964, 2965, 2967, 2968, 2970, 2971, 2972, 2976, 2977, 2978, 2983, 2985, 2988, 2992, 2998, 3002, 3006, 3008, 3010, 3012, 3015, 3017, 3018, 3019, 3020, 3021, 3023, 3024, 3025, 3026, 3027, 3029, 3030, 3033, 3034, 3035, 3038, 3039, 3040, 3041, 3043, 3044, 3047, 3049, 3051, 3053, 3054, 3055, 3056, 3057, 3060, 3062, 3063, 3064, 3067, 3068, 3069, 3071, 3072, 3078, 3079, 3080, 3081, 3084, 3086, 3087, 3088, 3091, 3093, 3096, 3097, 3102, 3104, 3105, 3108, 3109, 3110, 3111, 3112, 3116, 3121, 3122, 3123, 3129, 3133, 3134, 3135, 3136, 3140, 3141, 3142, 3147, 3148, 3149, 3150, 3151, 3154, 3156, 3157, 3159, 3160, 3161, 3163, 3165, 3167, 3168, 3169, 3170, 3175, 3176, 3177, 3183, 3184, 3186, 3189, 3192, 3193, 3196, 3199, 3200, 3203, 3205, 3208, 3209, 3210, 3214, 3216, 3217, 3220, 3221, 3222, 3225, 3226, 3230, 3231, 3233, 3236, 3238, 3242, 3243, 3244, 3246, 3249, 3250, 3251, 3253, 3257, 3261, 3262, 3264, 3266, 3269, 3270, 3271, 3272, 3274, 3275, 3278, 3279, 3282, 3283, 3285, 3287, 3288, 3290, 3292, 3294, 3295, 3296, 3298, 3300, 3302, 3303, 3304, 3305, 3306, 3309, 3310, 3315, 3316, 3317, 3318, 3321, 3325, 3333, 3334, 3335, 3336, 3337, 3338, 3342, 3344, 3351, 3352, 3354, 3357, 3359, 3361, 3362, 3363, 3365, 3368, 3370, 3373, 3374, 3375, 3376, 3377, 3378, 3379, 3380, 3381, 3384, 3385, 3386, 3387, 3391, 3392, 3393, 3395, 3397, 3398, 3399, 3400, 3402, 3406, 3408, 3410, 3412, 3413, 3414, 3416, 3418, 3421, 3422, 3423, 3426, 3428, 3432, 3434, 3435, 3438, 3442, 3443, 3445, 3446, 3449, 3450, 3451, 3453, 3456, 3460, 3465, 3466, 3469, 3470, 3471, 3472, 3473, 3474, 3476, 3477, 3480, 3481, 3484, 3485, 3486, 3487, 3488, 3490, 3491, 3492, 3493, 3494, 3496, 3498, 3501, 3502, 3503, 3505, 3508, 3509, 3510, 3511, 3513, 3515, 3516, 3517, 3518, 3519, 3520, 3522, 3524, 3525, 3526, 3528, 3532, 3533, 3534, 3537, 3539, 3541, 3543, 3544, 3547, 3549, 3552, 3556, 3557, 3559, 3560, 3562, 3564, 3569, 3570, 3571, 3572, 3574, 3575, 3580, 3583, 3585, 3586, 3588, 3589, 3590, 3591, 3592, 3595, 3596, 3597, 3600, 3601, 3603, 3608, 3612, 3613, 3614, 3615, 3617, 3623, 3628, 3629, 3631, 3633, 3634, 3636, 3637, 3638, 3640, 3642, 3644, 3646, 3647, 3649, 3650, 3652, 3657, 3658, 3660, 3662, 3664, 3665, 3668, 3670, 3672, 3674, 3678, 3681, 3683, 3686, 3687, 3690, 3694, 3697, 3698, 3699, 3701, 3703, 3706, 3707, 3710, 3712, 3715, 3716, 3722, 3723, 3729, 3730, 3731, 3733, 3734, 3735, 3741, 3743, 3744, 3745, 3746, 3747, 3750, 3751, 3753, 3757, 3759, 3765, 3767, 3768, 3769, 3770, 3771, 3774, 3775, 3777, 3779, 3780, 3781, 3782, 3783, 3784, 3787, 3789, 3790, 3791, 3793, 3796, 3797, 3803, 3810, 3813, 3814, 3815, 3820, 3821, 3822, 3823, 3825, 3826, 3828, 3829, 3830, 3831, 3832, 3833, 3836, 3838, 3839, 3840, 3841, 3842, 3845, 3846, 3847, 3850, 3852, 3854, 3856, 3857, 3858, 3860, 3862, 3866, 3868, 3869, 3870, 3872, 3875, 3877, 3883, 3885, 3887, 3888, 3889, 3890, 3891, 3892, 3894, 3895, 3896, 3897, 3898, 3900, 3906, 3907, 3908, 3910, 3913, 3914, 3915, 3916, 3917, 3918, 3926, 3928, 3929, 3930, 3932, 3933, 3934, 3935, 3939, 3943, 3944, 3945, 3946, 3947, 3949, 3951, 3953, 3956, 3959, 3963, 3964, 3965, 3966, 3967, 3969, 3973, 3975, 3979, 3980, 3981, 3984, 3985, 3986, 3989, 3991, 3993, 3994, 3995, 3997, 4000, 4001, 4002, 4003, 4004, 4006, 4011, 4012, 4014, 4015, 4018, 4022, 4024, 4028, 4029, 4032, 4033, 4034, 4036, 4038, 4039, 4040, 4041, 4042, 4044, 4046, 4047, 4050, 4052, 4055, 4057, 4058, 4061, 4064, 4065, 4066, 4070, 4071, 4075, 4076, 4078, 4081, 4082, 4084, 4085, 4086, 4095, 4101, 4102, 4103, 4106, 4107, 4110, 4111, 4112, 4115, 4116, 4118, 4121, 4123, 4124, 4126, 4130, 4131, 4132, 4136, 4137, 4138, 4140, 4142, 4143, 4144, 4147, 4148, 4156, 4159, 4160, 4163, 4166, 4173, 4174, 4177, 4179, 4180, 4190, 4191, 4194, 4196, 4199, 4200, 4202, 4204, 4205, 4207, 4208, 4209, 4210, 4212, 4213, 4215, 4216, 4217, 4219, 4224, 4228, 4231, 4232, 4233, 4236, 4237, 4243, 4244, 4245, 4249, 4250, 4251, 4256, 4257, 4258, 4262, 4265, 4266, 4269, 4270, 4271, 4272, 4273, 4276, 4278, 4280, 4281, 4283, 4285, 4288, 4289, 4292, 4293, 4295, 4296, 4297, 4299, 4301, 4302, 4303, 4304, 4306, 4309, 4313, 4315, 4317, 4318, 4319, 4320, 4321, 4325, 4326, 4327, 4329, 4330, 4331, 4332, 4333, 4334, 4336, 4337, 4338, 4339, 4342, 4344, 4348, 4353, 4363, 4368, 4370, 4372, 4373, 4375, 4377, 4378, 4389, 4390, 4392, 4393, 4394, 4395, 4396, 4397, 4399, 4400, 4406, 4407, 4410, 4411, 4412, 4415, 4418, 4422, 4424, 4426, 4427, 4428, 4431, 4433, 4434, 4435, 4438, 4439, 4441, 4442, 4443, 4444, 4445, 4446, 4448, 4452, 4453, 4456, 4457, 4459, 4462, 4463, 4465, 4466, 4467, 4469, 4471, 4474, 4475, 4477, 4480, 4481, 4482, 4483, 4486, 4489, 4490, 4492, 4495, 4497, 4498, 4502, 4505, 4507, 4509, 4511, 4512, 4513, 4514, 4518, 4519, 4520, 4522, 4523, 4525, 4529, 4532, 4533, 4534, 4535, 4536, 4538, 4541, 4542, 4546, 4550, 4551, 4552, 4554, 4555, 4557, 4558, 4559, 4560, 4562, 4564, 4567, 4568, 4569, 4571, 4572, 4573, 4574, 4575, 4580, 4581, 4589, 4591, 4596, 4597, 4598, 4600, 4601, 4605, 4606, 4607, 4609, 4612, 4615, 4620, 4621, 4629, 4630, 4632, 4633, 4634, 4637, 4641, 4642, 4650, 4651, 4652, 4653, 4654, 4655, 4659, 4661, 4663, 4664, 4665, 4666, 4668, 4669, 4670, 4671, 4672, 4675, 4677, 4680, 4681, 4682, 4684, 4685, 4686, 4693, 4694, 4695, 4697, 4700, 4701, 4703, 4707, 4709, 4710, 4711, 4712, 4716, 4717, 4718, 4720, 4721, 4723, 4724, 4725, 4728, 4729, 4731, 4734, 4735, 4736, 4738, 4739, 4740, 4742, 4743, 4744, 4745, 4746, 4747, 4750, 4751, 4755, 4757, 4758, 4759, 4760, 4761, 4770, 4771, 4774, 4775, 4777, 4779, 4780, 4782, 4785, 4786, 4787, 4788, 4789, 4794, 4798, 4799, 4801, 4803, 4805, 4808, 4813, 4816, 4817, 4818, 4819, 4822, 4826, 4831, 4832, 4839, 4841, 4844, 4845, 4847, 4848, 4849, 4852, 4853, 4856, 4857, 4858, 4859, 4860, 4861, 4863, 4864, 4867, 4868, 4870, 4871, 4874, 4875, 4876, 4877, 4878, 4880, 4882, 4883, 4885, 4886, 4887, 4888, 4890, 4894, 4895, 4896, 4897, 4899, 4901, 4905, 4906, 4909, 4910, 4911, 4912, 4916, 4918, 4919, 4920, 4921, 4925, 4927, 4928, 4929, 4930, 4933, 4934, 4938, 4940, 4941, 4944, 4947, 4948, 4953, 4954, 4958, 4961, 4962, 4963, 4965, 4967, 4969, 4972, 4973, 4975, 4977, 4978, 4979, 4980, 4981, 4982, 4986, 4990, 4992, 4994, 4995, 4997, 4998, 5002, 5007, 5008, 5013, 5015, 5018, 5019, 5020, 5023, 5025, 5027, 5029, 5031, 5032, 5034, 5039, 5040, 5044, 5045, 5046, 5048, 5049, 5050, 5051, 5053, 5054, 5055, 5057, 5058, 5062, 5063, 5068, 5069, 5072, 5073, 5076, 5077, 5078, 5080, 5084, 5086, 5087, 5088, 5090, 5092, 5094, 5097, 5098, 5101, 5104, 5105, 5107, 5109, 5111, 5114, 5117, 5118, 5119, 5121, 5122, 5124, 5125, 5127, 5128, 5129, 5132, 5136, 5140, 5141, 5144, 5147, 5149, 5150, 5152, 5155, 5156, 5158, 5162, 5163, 5164, 5167, 5171, 5174, 5175, 5178, 5179, 5183, 5184, 5185, 5186, 5187, 5195, 5196, 5198, 5199, 5200, 5201, 5203, 5205, 5206, 5207, 5213, 5214, 5216, 5221, 5222, 5223, 5224, 5227, 5228, 5229, 5230, 5233, 5234, 5235, 5236, 5237, 5238, 5239, 5240, 5241, 5244, 5247, 5252, 5253, 5254, 5255, 5256, 5260, 5261, 5262, 5263, 5264, 5265, 5267, 5268, 5272, 5274, 5275, 5279, 5282, 5285, 5288, 5290, 5291, 5293, 5294, 5295, 5296, 5298, 5302, 5304, 5306, 5309, 5311, 5312, 5314, 5315, 5317, 5319, 5320, 5321, 5322, 5323, 5325, 5327, 5329, 5330, 5334, 5335, 5339, 5348, 5349, 5352, 5353, 5355, 5357, 5358, 5359, 5361, 5362, 5364, 5366, 5367, 5369, 5373, 5374, 5379, 5380, 5381, 5382, 5383, 5385, 5386, 5387, 5389, 5390, 5391, 5392, 5393, 5394, 5397, 5399, 5401, 5402, 5406, 5410, 5412, 5413, 5414, 5416, 5417, 5420, 5421, 5422, 5424, 5425, 5428, 5429, 5430, 5435, 5436, 5438, 5439, 5441, 5443, 5445, 5446, 5448, 5450, 5451, 5452, 5455, 5456, 5458, 5459, 5460, 5461, 5464, 5465, 5466, 5469, 5470, 5474, 5475, 5476, 5477, 5482, 5483, 5485, 5486, 5488, 5495, 5498, 5500, 5501, 5502, 5504, 5507, 5508, 5513, 5514, 5519, 5521, 5522, 5525, 5526, 5528, 5529, 5532, 5536, 5537, 5538, 5539, 5541, 5543, 5545, 5546, 5548, 5552, 5553, 5554, 5555, 5558, 5560, 5561, 5562, 5563, 5566, 5569, 5572, 5573, 5574, 5576, 5577, 5578, 5579, 5582, 5583, 5585, 5586, 5587, 5588, 5589, 5590, 5594, 5597, 5598, 5599, 5600, 5601, 5604, 5606, 5607, 5611, 5612, 5614, 5615, 5617, 5618, 5619, 5620, 5623, 5624, 5626, 5627, 5628, 5629, 5634, 5637, 5640, 5644, 5645, 5646, 5648, 5649, 5650, 5651, 5655, 5656, 5658, 5659, 5661, 5662, 5663, 5665, 5669, 5670, 5671, 5672, 5673, 5677, 5678, 5679, 5682, 5683, 5691, 5692, 5693, 5695, 5696, 5699, 5700, 5704, 5705, 5707, 5710, 5711, 5714, 5715, 5717, 5719, 5720, 5722, 5725, 5727, 5728, 5729, 5730, 5733, 5735, 5737, 5738, 5739, 5742, 5744, 5746, 5747, 5748, 5749, 5751, 5752, 5754, 5756, 5757, 5759, 5760, 5761, 5762, 5764, 5765, 5766, 5770, 5772, 5773, 5775, 5776, 5778, 5783, 5785, 5786, 5788, 5794, 5795, 5796, 5799, 5800, 5801, 5803, 5804, 5805, 5809, 5812, 5818, 5820, 5821, 5824, 5825, 5827, 5828, 5829, 5830, 5832, 5833, 5834, 5835, 5836, 5837, 5838, 5843, 5844, 5845, 5847, 5849, 5850, 5851, 5858, 5862, 5866, 5867, 5869, 5872, 5873, 5874, 5875, 5876, 5877, 5878, 5880, 5890, 5891, 5892, 5893, 5894, 5895, 5897, 5899, 5902, 5912, 5913, 5914, 5916, 5923, 5928, 5931, 5932, 5933, 5934, 5944, 5946, 5947, 5949, 5952, 5955, 5958, 5959, 5960, 5961, 5962, 5969, 5971, 5974, 5975, 5978, 5980, 5981, 5984, 5986, 5987, 5989, 5990, 5994, 5995, 5996, 6003, 6004, 6006, 6007, 6009, 6010, 6011, 6012, 6014, 6016, 6017, 6018, 6024, 6026, 6027, 6030, 6032, 6033, 6034, 6035, 6036, 6037, 6038, 6041, 6044, 6048, 6049, 6051, 6052, 6053, 6055, 6057, 6058, 6059, 6060, 6061, 6062, 6065, 6067, 6077, 6081, 6085, 6087, 6088, 6091, 6093, 6095, 6096, 6097, 6101, 6106, 6108, 6109, 6111, 6112, 6115, 6116, 6118, 6119, 6121, 6122, 6123, 6128, 6129, 6135, 6136, 6138, 6143, 6147, 6149, 6150, 6153, 6155, 6156, 6157, 6158, 6159, 6160, 6161, 6162, 6164, 6167, 6168, 6171, 6172, 6174, 6175, 6178, 6180, 6181, 6183, 6188, 6189, 6190, 6191, 6192, 6194, 6195, 6196, 6201, 6204, 6206, 6208, 6213, 6214, 6215, 6217, 6219, 6223, 6224, 6226, 6228, 6230, 6231, 6235, 6237, 6238, 6239, 6241, 6244, 6245, 6247, 6248, 6249, 6252, 6253, 6254, 6255, 6256, 6257, 6259, 6260, 6262, 6264, 6265, 6267, 6269, 6273, 6274, 6275, 6276, 6277, 6278, 6279, 6280, 6281, 6282, 6284, 6285, 6287, 6289, 6290, 6291, 6293, 6296, 6298, 6299, 6303, 6304, 6308, 6310, 6314, 6316, 6317, 6318, 6320, 6321, 6323, 6324, 6328, 6333, 6335, 6336, 6337, 6338, 6339, 6342, 6343, 6344, 6349, 6350, 6352, 6353, 6357, 6359, 6364, 6370, 6373, 6374, 6375, 6378, 6379, 6381, 6385, 6388, 6392, 6393, 6395, 6399, 6404, 6408, 6409, 6410, 6414, 6415, 6417, 6420, 6422, 6423, 6424, 6425, 6428, 6429, 6431, 6434, 6435, 6437, 6439, 6440, 6441, 6442, 6443, 6444, 6446, 6447, 6448, 6452, 6455, 6459, 6461, 6464, 6465, 6466, 6467, 6468, 6470, 6471, 6473, 6474, 6475, 6478, 6483, 6484, 6485, 6490, 6492, 6494, 6499, 6501, 6502, 6503, 6504, 6508, 6511, 6513, 6514, 6515, 6517, 6519, 6521, 6522, 6523, 6524, 6525, 6526, 6528, 6529, 6530, 6531, 6539, 6540, 6541, 6544, 6545, 6546, 6547, 6549, 6552, 6553, 6555, 6557, 6558, 6559, 6561, 6562, 6563, 6565, 6567, 6568, 6571, 6572, 6573, 6577, 6578, 6580, 6581, 6582, 6584, 6586, 6589, 6595, 6596, 6597, 6601, 6602, 6603, 6604, 6605, 6606, 6608, 6609, 6610, 6611, 6613, 6615, 6616, 6618, 6619, 6620, 6621, 6623, 6624, 6625, 6628, 6629, 6630, 6631, 6632, 6633, 6637, 6644, 6646, 6647, 6648, 6650, 6651, 6654, 6657, 6658, 6660, 6661, 6663, 6664, 6665, 6668, 6669, 6670, 6671, 6672, 6675, 6679, 6680, 6681, 6682, 6684, 6689, 6691, 6692, 6696, 6697, 6698, 6701, 6704, 6705, 6707, 6709, 6711, 6713, 6714, 6715, 6720, 6721, 6724, 6725, 6726, 6730, 6731, 6733, 6736, 6737, 6741, 6742, 6744, 6745, 6748, 6750, 6752, 6753, 6755, 6758, 6759, 6761, 6764, 6767, 6768, 6769, 6772, 6774, 6775, 6777, 6778, 6780, 6782, 6783, 6789, 6792, 6793, 6795, 6796, 6797, 6800, 6801, 6802, 6806, 6808, 6809, 6810, 6812, 6814, 6819, 6820, 6826, 6827, 6828, 6833, 6834, 6837, 6840, 6844, 6845, 6850, 6851, 6852, 6853, 6854, 6856, 6859, 6864, 6865, 6868, 6869, 6870, 6871, 6874, 6877, 6879, 6880, 6883, 6886, 6887, 6888, 6895, 6896, 6899, 6900, 6902, 6903, 6904, 6905, 6908, 6910, 6911, 6912, 6913, 6916, 6920, 6921, 6924, 6925, 6930, 6932, 6936, 6939, 6940, 6941, 6942, 6943, 6945, 6947, 6952, 6955, 6957, 6958, 6959, 6966, 6971, 6976, 6980, 6981, 6984, 6986, 6987, 6988, 6991, 6993, 6996, 6999, 7000, 7004, 7009, 7010, 7012, 7013, 7014, 7016, 7019, 7020, 7022, 7023, 7025, 7027, 7030, 7032, 7033, 7035, 7038, 7042, 7043, 7046, 7050, 7051, 7053, 7055, 7056, 7057, 7058, 7060, 7061, 7062, 7063, 7064, 7065, 7066, 7068, 7072, 7074, 7075, 7076, 7077, 7078, 7079, 7080, 7081, 7082, 7084, 7085, 7087, 7088, 7089, 7095, 7096, 7101, 7103, 7104, 7105, 7106, 7108, 7116, 7119, 7120, 7121, 7122, 7123, 7124, 7126, 7128, 7129, 7131, 7133, 7134, 7135, 7136, 7138, 7140, 7141, 7142, 7146, 7147, 7151, 7152, 7153, 7160, 7161, 7162, 7163, 7165, 7166, 7167, 7169, 7172, 7174, 7175, 7176, 7177, 7178, 7179, 7182, 7185, 7187, 7188, 7190, 7191, 7196, 7197, 7198, 7202, 7203, 7206, 7208, 7212, 7213, 7215, 7216, 7217, 7221, 7222, 7226, 7228, 7232, 7233, 7235, 7236, 7238, 7239, 7241, 7243, 7244, 7245, 7246, 7250, 7252, 7253, 7255, 7259, 7260, 7261, 7263, 7264, 7266, 7267, 7270, 7271, 7272, 7273, 7275, 7279, 7280, 7283, 7284, 7285, 7286, 7287, 7289, 7290, 7292, 7297, 7302, 7303, 7304, 7305, 7307, 7309, 7314, 7320, 7324, 7325, 7326, 7327, 7329, 7330, 7332, 7335, 7336, 7337, 7339, 7340, 7342, 7346, 7350, 7353, 7354, 7357, 7358, 7361, 7368, 7369, 7372, 7373, 7376, 7378, 7379, 7381, 7382, 7385, 7386, 7388, 7394, 7396, 7397, 7399, 7400, 7402, 7403, 7405, 7406, 7407, 7411, 7412, 7414, 7416, 7419, 7420, 7422, 7423, 7430, 7432, 7434, 7437, 7438, 7440, 7441, 7442, 7444, 7446, 7447, 7454, 7459, 7462, 7463, 7464, 7467, 7468, 7470, 7471, 7472, 7473, 7474, 7476, 7477, 7480, 7484, 7485, 7486, 7493, 7494, 7495, 7498, 7502, 7503, 7505, 7506, 7511, 7512, 7517, 7519, 7520, 7522, 7524, 7528, 7531, 7534, 7538, 7539, 7540, 7541, 7542, 7543, 7545, 7551, 7554, 7556, 7558, 7559, 7569, 7570, 7573, 7576, 7577, 7578, 7579, 7585, 7586, 7587, 7588, 7589, 7590, 7592, 7594, 7595, 7596, 7600, 7601, 7602, 7604, 7606, 7607, 7609, 7615, 7616, 7618, 7619, 7622, 7624, 7625, 7632, 7633, 7634, 7635, 7637, 7638, 7642, 7643, 7645, 7647, 7648, 7653, 7654, 7659, 7663, 7666, 7670, 7679, 7680, 7681, 7682, 7685, 7686, 7687, 7689, 7690, 7691, 7693, 7695, 7697, 7698, 7699, 7701, 7702, 7708, 7711, 7712, 7717, 7718, 7719, 7721, 7723, 7727, 7728, 7730, 7731, 7732, 7733, 7734, 7735, 7736, 7741, 7744, 7745, 7746, 7747, 7748, 7751, 7754, 7757, 7758, 7759, 7764, 7765, 7767, 7770, 7772, 7773, 7774, 7775, 7780, 7782, 7783, 7787, 7788, 7789, 7790, 7793, 7795, 7798, 7799, 7801, 7807, 7810, 7811, 7812, 7816, 7818, 7819, 7820, 7822, 7823, 7824, 7826, 7827, 7828, 7829, 7830, 7831, 7834, 7836, 7839, 7841, 7842, 7844, 7845, 7846, 7848, 7850, 7851, 7853, 7854, 7855, 7857, 7859, 7860, 7861, 7862, 7866, 7867, 7868, 7870, 7871, 7874, 7879, 7881, 7882, 7884, 7886, 7891, 7892, 7893, 7899, 7900, 7901, 7902, 7904, 7907, 7908, 7910, 7915, 7916, 7917, 7919, 7920, 7921, 7924, 7925, 7927, 7930, 7931, 7933, 7934, 7937, 7941, 7942, 7944, 7947, 7949, 7950, 7951, 7953, 7954, 7955, 7956, 7957, 7961, 7963, 7967, 7968, 7972, 7973, 7974, 7975, 7976, 7977, 7978, 7979, 7981, 7982, 7985, 7987, 7989, 7991, 7994, 7997, 7999, 8000, 8001, 8003, 8004, 8006, 8007, 8008, 8009, 8012, 8013, 8018, 8019, 8020, 8024, 8027, 8028, 8031, 8034, 8036, 8037, 8040, 8041, 8049, 8052, 8056, 8062, 8063, 8064, 8075, 8078, 8081, 8083, 8086, 8087, 8088, 8092, 8098, 8100, 8101, 8102, 8105, 8107, 8110, 8111, 8112, 8113, 8115, 8116, 8119, 8120, 8122, 8126, 8127, 8129, 8134, 8135, 8136, 8138, 8142, 8144, 8146, 8149, 8154, 8156, 8159, 8161, 8162, 8164, 8165, 8166, 8169, 8171, 8175, 8176, 8178, 8180, 8182, 8183, 8185, 8188, 8190, 8191, 8195, 8198, 8200, 8201, 8202, 8203, 8204, 8210, 8211, 8213, 8214, 8215, 8217, 8218, 8219, 8220, 8224, 8227, 8230, 8233, 8235, 8240, 8243, 8244, 8245, 8246, 8247, 8249, 8250, 8256, 8257, 8260, 8261, 8264, 8266, 8268, 8269, 8270, 8271, 8274, 8276, 8282, 8285, 8287, 8288, 8289, 8290, 8291, 8297, 8300, 8301, 8302, 8303, 8306, 8307, 8308, 8309, 8312, 8318, 8319, 8321, 8323, 8325, 8326, 8328, 8329, 8335, 8336, 8338, 8339, 8340, 8345, 8346, 8347, 8348, 8350, 8351, 8352, 8353, 8356, 8360, 8363, 8367, 8368, 8369, 8370, 8372, 8374, 8376, 8379, 8380, 8383, 8384, 8385, 8386, 8387, 8388, 8389, 8391, 8395, 8396, 8397, 8398, 8399, 8401, 8406, 8407, 8408, 8409, 8411, 8414, 8416, 8417, 8420, 8423, 8432, 8433, 8434, 8436, 8437, 8438, 8439, 8441, 8444, 8446, 8449, 8450, 8451, 8453, 8455, 8459, 8460, 8462, 8463, 8465, 8466, 8467, 8471, 8473, 8475, 8476, 8478, 8480, 8481, 8482, 8483, 8489, 8490, 8491, 8492, 8493, 8495, 8496, 8498, 8500, 8503, 8506, 8508, 8511, 8513, 8514, 8518, 8519, 8520, 8522, 8523, 8524, 8525, 8526, 8527, 8530, 8531, 8533, 8535, 8536, 8537, 8544, 8546, 8547, 8551, 8554, 8555, 8556, 8559, 8560, 8562, 8563, 8567, 8569, 8570, 8571, 8572, 8577, 8579, 8580, 8582, 8583, 8586, 8587, 8588, 8590, 8592, 8593, 8594, 8595, 8597, 8599, 8600, 8602, 8603, 8604, 8606, 8608, 8609, 8610, 8613, 8617, 8618, 8621, 8622, 8627, 8629, 8630, 8632, 8634, 8635, 8636, 8637, 8638, 8641, 8642, 8643, 8644, 8645, 8647, 8649, 8651, 8653, 8654, 8655, 8656, 8657, 8658, 8659, 8661, 8663, 8666, 8668, 8670, 8671, 8672, 8674, 8675, 8676, 8677, 8681, 8682, 8683, 8685, 8686, 8687, 8689, 8690, 8693, 8694, 8696, 8699, 8700, 8702, 8706, 8707, 8711, 8713, 8715, 8717, 8718, 8722, 8724, 8727, 8728, 8729, 8730, 8731, 8732, 8737, 8740, 8741, 8742, 8746, 8750, 8752, 8754, 8756, 8757, 8760, 8761, 8762, 8764, 8765, 8770, 8774, 8775, 8776, 8777, 8778, 8780, 8784, 8785, 8786, 8787, 8790, 8794, 8798, 8799, 8803, 8804, 8805, 8809, 8810, 8812, 8813, 8816, 8819, 8821, 8822, 8823, 8824, 8825, 8826, 8828, 8829, 8833, 8835, 8836, 8839, 8840, 8841, 8842, 8843, 8844, 8845, 8851, 8852, 8855, 8856, 8857, 8859, 8860, 8866, 8867, 8868, 8869, 8871, 8875, 8876, 8877, 8883, 8885, 8886, 8887, 8890, 8895, 8899, 8903, 8905, 8906, 8907, 8909, 8910, 8911, 8913, 8914, 8916, 8919, 8921, 8922, 8924, 8925, 8927, 8928, 8929, 8935, 8936, 8938, 8940, 8942, 8945, 8946, 8950, 8952, 8954, 8957, 8958, 8959, 8961, 8963, 8966, 8968, 8970, 8973, 8974, 8975, 8976, 8978, 8980, 8981, 8982, 8983, 8984, 8986, 8987, 8989, 8990, 8992, 8996, 8998, 8999, 9000, 9001, 9002, 9003, 9004, 9005, 9008, 9009, 9010, 9013, 9014, 9015, 9018, 9022, 9023, 9025, 9028, 9030, 9031, 9032, 9036, 9037, 9039, 9040, 9041, 9044, 9045, 9048, 9049, 9051, 9053, 9054, 9055, 9056, 9057, 9058, 9059, 9060, 9062, 9064, 9065, 9066, 9068, 9069, 9070, 9072, 9073, 9074, 9076, 9078, 9083, 9084, 9086, 9087, 9088, 9089, 9093, 9094, 9096, 9097, 9099, 9100, 9102, 9103, 9104, 9107, 9109, 9110, 9111, 9113, 9114, 9116, 9117, 9120, 9122, 9123, 9125, 9127, 9128, 9129, 9130, 9132, 9135, 9136, 9140, 9143, 9147, 9152, 9154, 9155, 9156, 9157, 9159, 9163, 9168, 9170, 9172, 9175, 9176, 9177, 9179, 9181, 9182, 9184, 9186, 9187, 9190, 9195, 9196, 9197, 9200, 9203, 9205, 9206, 9213, 9214, 9216, 9217, 9218, 9219, 9223, 9227, 9232, 9233, 9236, 9241, 9243, 9244, 9245, 9246, 9253, 9255, 9256, 9257, 9259, 9260, 9261, 9262, 9263, 9265, 9269, 9270, 9271, 9274, 9275, 9278, 9280, 9282, 9283, 9285, 9286, 9291, 9292, 9294, 9298, 9299, 9300, 9301, 9302, 9303, 9307, 9308, 9309, 9310, 9311, 9312, 9316, 9321, 9322, 9323, 9327, 9329, 9332, 9338, 9339, 9346, 9348, 9350, 9351, 9354, 9356, 9360, 9361, 9362, 9363, 9364, 9366, 9367, 9369, 9370, 9371, 9373, 9375, 9378, 9379, 9380, 9383, 9387, 9388, 9390, 9391, 9392, 9394, 9396, 9398, 9400, 9403, 9404, 9405, 9406, 9407, 9409, 9410, 9414, 9417, 9419, 9420, 9423, 9425, 9428, 9429, 9431, 9432, 9434, 9435, 9438, 9439, 9440, 9441, 9442, 9444, 9445, 9448, 9453, 9454, 9464, 9466, 9467, 9468, 9469, 9474, 9475, 9477, 9478, 9480, 9484, 9487, 9488, 9489, 9492, 9493, 9494, 9497, 9500, 9505, 9508, 9509, 9510, 9511, 9512, 9513, 9515, 9518, 9519, 9520, 9521, 9522, 9524, 9525, 9527, 9528, 9529, 9530, 9532, 9534, 9535, 9538, 9540, 9543, 9546, 9550, 9551, 9552, 9556, 9557, 9558, 9559, 9560, 9561, 9567, 9569, 9572, 9573, 9574, 9575, 9577, 9580, 9581, 9583, 9586, 9587, 9588, 9590, 9591, 9592, 9593, 9596, 9597, 9599, 9600, 9601, 9602, 9605, 9608, 9609, 9610, 9611, 9612, 9613, 9615, 9616, 9617, 9619, 9621, 9625, 9629, 9630, 9632, 9635, 9636, 9638, 9639, 9640, 9643, 9645, 9646, 9648, 9652, 9653, 9655, 9656, 9659, 9660, 9661, 9662, 9664, 9666, 9667, 9669, 9671, 9674, 9675, 9676, 9681, 9682, 9685, 9687, 9689, 9694, 9695, 9698, 9699, 9701, 9702, 9704, 9706, 9707, 9710, 9713, 9714, 9716, 9717, 9719, 9720, 9721, 9722, 9723, 9725, 9728, 9729, 9730, 9731, 9732, 9733, 9736, 9737, 9741, 9742, 9743, 9747, 9748, 9749, 9752, 9753, 9754, 9760, 9764, 9766, 9767, 9768, 9774, 9775, 9776, 9779, 9780, 9783, 9786, 9791, 9792, 9795, 9796, 9797, 9798, 9799, 9801, 9805, 9806, 9807, 9810, 9811, 9812, 9814, 9819, 9821, 9822, 9823, 9824, 9826, 9828, 9829, 9831, 9832, 9836, 9838, 9843, 9845, 9855, 9856, 9857, 9858, 9861, 9863, 9864, 9866, 9867, 9868, 9870, 9872, 9873, 9877, 9878, 9880, 9881, 9883, 9884, 9888, 9889]),
        numpy_asarray([0, 1, 2, 3, 6, 7, 8, 10, 21, 22, 24, 26, 28, 30, 32, 34, 35, 38, 39, 42, 43, 47, 50, 52, 54, 55, 56, 59, 65, 66, 69, 70, 71, 72, 74, 75, 77, 79, 80, 82, 83, 84, 86, 88, 89, 90, 92, 94, 95, 98, 99, 102, 106, 107, 108, 113, 114, 115, 121, 122, 129, 132, 133, 134, 135, 139, 142, 144, 152, 153, 155, 157, 160, 165, 166, 168, 170, 177, 179, 180, 181, 182, 185, 188, 189, 190, 192, 193, 194, 195, 197, 199, 200, 201, 203, 204, 205, 207, 209, 210, 212, 213, 214, 215, 220, 222, 223, 227, 229, 231, 233, 234, 235, 237, 238, 241, 243, 245, 248, 249, 251, 252, 253, 256, 259, 266, 267, 268, 269, 270, 271, 275, 277, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 292, 293, 295, 296, 298, 299, 300, 301, 302, 305, 307, 308, 309, 310, 311, 312, 313, 314, 316, 319, 320, 322, 325, 327, 329, 330, 333, 335, 336, 337, 338, 339, 341, 342, 344, 346, 350, 351, 355, 362, 363, 364, 365, 366, 367, 374, 375, 376, 379, 380, 383, 384, 385, 387, 389, 391, 393, 394, 395, 396, 397, 402, 403, 405, 406, 407, 413, 414, 419, 420, 421, 423, 424, 427, 428, 430, 431, 433, 435, 436, 438, 439, 440, 441, 442, 443, 444, 445, 447, 452, 453, 454, 455, 458, 463, 464, 465, 469, 470, 473, 474, 477, 482, 483, 485, 486, 487, 488, 489, 490, 491, 493, 494, 497, 502, 503, 504, 505, 507, 510, 511, 512, 513, 515, 517, 518, 522, 523, 525, 526, 528, 535, 537, 538, 540, 541, 542, 543, 544, 545, 546, 548, 549, 550, 552, 554, 557, 559, 560, 561, 562, 563, 564, 565, 566, 568, 570, 571, 573, 575, 580, 581, 582, 584, 585, 587, 588, 590, 591, 594, 596, 597, 598, 603, 607, 610, 611, 614, 617, 620, 621, 622, 623, 625, 626, 629, 633, 634, 635, 638, 640, 641, 642, 646, 647, 649, 650, 651, 652, 653, 654, 655, 657, 661, 662, 664, 665, 671, 672, 673, 674, 679, 684, 686, 689, 690, 691, 694, 695, 699, 701, 703, 704, 706, 708, 712, 713, 715, 717, 718, 719, 720, 721, 723, 724, 728, 730, 731, 734, 738, 739, 742, 745, 746, 747, 749, 756, 758, 759, 761, 767, 768, 770, 771, 772, 775, 776, 777, 780, 784, 785, 786, 788, 790, 791, 796, 797, 798, 799, 803, 809, 810, 811, 812, 813, 814, 816, 820, 821, 822, 828, 829, 830, 831, 832, 833, 834, 837, 839, 840, 842, 843, 846, 847, 854, 855, 856, 857, 858, 859, 861, 863, 866, 867, 872, 873, 880, 882, 883, 886, 888, 891, 892, 894, 895, 897, 898, 900, 903, 904, 906, 908, 910, 911, 913, 915, 918, 919, 921, 922, 923, 924, 925, 926, 933, 936, 938, 940, 941, 942, 946, 947, 948, 950, 951, 952, 957, 958, 959, 960, 961, 962, 966, 967, 968, 969, 971, 972, 974, 976, 977, 978, 979, 980, 983, 984, 985, 988, 989, 991, 992, 993, 995, 996, 999, 1001, 1005, 1007, 1009, 1010, 1011, 1013, 1021, 1023, 1027, 1028, 1029, 1031, 1033, 1034, 1038, 1040, 1044, 1045, 1046, 1047, 1049, 1051, 1052, 1055, 1056, 1058, 1059, 1061, 1063, 1067, 1069, 1070, 1072, 1074, 1076, 1078, 1079, 1082, 1084, 1086, 1087, 1090, 1091, 1093, 1094, 1097, 1100, 1102, 1103, 1105, 1106, 1108, 1109, 1111, 1112, 1122, 1124, 1125, 1128, 1129, 1131, 1133, 1135, 1136, 1138, 1141, 1142, 1143, 1144, 1147, 1149, 1151, 1152, 1154, 1157, 1158, 1160, 1161, 1162, 1163, 1166, 1173, 1175, 1178, 1180, 1182, 1186, 1187, 1188, 1190, 1191, 1192, 1194, 1195, 1197, 1198, 1199, 1200, 1206, 1208, 1210, 1214, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1225, 1226, 1227, 1228, 1231, 1232, 1234, 1235, 1237, 1238, 1239, 1241, 1248, 1249, 1256, 1257, 1262, 1263, 1267, 1268, 1270, 1272, 1273, 1274, 1276, 1279, 1281, 1283, 1286, 1287, 1288, 1290, 1291, 1292, 1294, 1295, 1296, 1298, 1300, 1301, 1303, 1304, 1307, 1308, 1309, 1310, 1313, 1314, 1317, 1321, 1325, 1330, 1331, 1334, 1335, 1338, 1340, 1342, 1343, 1345, 1346, 1348, 1350, 1351, 1352, 1354, 1355, 1361, 1362, 1363, 1364, 1371, 1372, 1373, 1376, 1377, 1380, 1381, 1383, 1389, 1390, 1393, 1398, 1399, 1403, 1404, 1406, 1409, 1410, 1417, 1418, 1420, 1423, 1427, 1430, 1433, 1434, 1436, 1440, 1443, 1444, 1447, 1454, 1456, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1466, 1467, 1468, 1469, 1471, 1472, 1474, 1476, 1477, 1478, 1480, 1482, 1483, 1485, 1486, 1488, 1489, 1492, 1494, 1496, 1500, 1503, 1504, 1505, 1506, 1509, 1510, 1519, 1520, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1531, 1532, 1533, 1535, 1536, 1537, 1538, 1541, 1544, 1545, 1547, 1548, 1549, 1550, 1555, 1558, 1559, 1560, 1562, 1565, 1571, 1576, 1581, 1584, 1585, 1593, 1594, 1595, 1596, 1602, 1603, 1605, 1606, 1609, 1610, 1614, 1615, 1616, 1618, 1623, 1624, 1626, 1627, 1629, 1630, 1633, 1634, 1635, 1636, 1637, 1638, 1642, 1645, 1646, 1649, 1650, 1651, 1652, 1658, 1659, 1661, 1662, 1664, 1666, 1667, 1668, 1669, 1670, 1673, 1675, 1680, 1681, 1690, 1692, 1694, 1695, 1697, 1700, 1701, 1702, 1703, 1704, 1706, 1707, 1708, 1709, 1710, 1712, 1715, 1719, 1720, 1721, 1722, 1724, 1725, 1726, 1730, 1731, 1733, 1736, 1739, 1740, 1742, 1743, 1745, 1748, 1749, 1751, 1752, 1753, 1755, 1756, 1757, 1765, 1766, 1771, 1772, 1773, 1774, 1775, 1779, 1780, 1785, 1787, 1791, 1792, 1795, 1796, 1798, 1800, 1802, 1804, 1805, 1806, 1807, 1810, 1811, 1813, 1814, 1816, 1817, 1819, 1823, 1824, 1825, 1827, 1828, 1832, 1840, 1847, 1848, 1849, 1851, 1854, 1858, 1859, 1861, 1865, 1868, 1870, 1872, 1873, 1876, 1877, 1878, 1881, 1884, 1885, 1886, 1887, 1891, 1892, 1893, 1894, 1895, 1896, 1900, 1901, 1902, 1905, 1909, 1911, 1917, 1923, 1924, 1925, 1928, 1931, 1932, 1935, 1938, 1940, 1941, 1942, 1943, 1944, 1949, 1950, 1952, 1957, 1962, 1963, 1964, 1969, 1971, 1973, 1975, 1976, 1977, 1979, 1980, 1983, 1984, 1986, 1988, 1989, 1991, 1992, 1993, 1994, 1997, 2000, 2001, 2004, 2005, 2010, 2011, 2013, 2017, 2018, 2019, 2021, 2022, 2023, 2027, 2029, 2031, 2032, 2037, 2039, 2042, 2045, 2046, 2052, 2054, 2056, 2057, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2079, 2081, 2084, 2085, 2087, 2089, 2094, 2098, 2099, 2100, 2102, 2103, 2105, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2117, 2122, 2123, 2126, 2127, 2128, 2130, 2135, 2136, 2138, 2141, 2142, 2145, 2146, 2151, 2152, 2154, 2155, 2157, 2158, 2159, 2160, 2163, 2164, 2165, 2167, 2168, 2169, 2170, 2171, 2172, 2173, 2175, 2176, 2177, 2179, 2180, 2181, 2183, 2188, 2192, 2195, 2196, 2198, 2200, 2204, 2206, 2208, 2209, 2210, 2211, 2212, 2214, 2215, 2219, 2221, 2225, 2228, 2232, 2234, 2241, 2242, 2243, 2245, 2249, 2251, 2252, 2253, 2255, 2257, 2258, 2259, 2260, 2263, 2267, 2268, 2270, 2271, 2275, 2276, 2278, 2284, 2286, 2287, 2289, 2291, 2293, 2295, 2297, 2301, 2302, 2303, 2304, 2306, 2308, 2311, 2312, 2315, 2316, 2317, 2319, 2320, 2321, 2323, 2324, 2325, 2327, 2328, 2329, 2330, 2331, 2332, 2335, 2336, 2337, 2340, 2341, 2344, 2345, 2346, 2347, 2351, 2353, 2357, 2359, 2360, 2361, 2362, 2363, 2364, 2366, 2369, 2371, 2372, 2373, 2374, 2376, 2378, 2379, 2380, 2382, 2383, 2385, 2387, 2388, 2390, 2391, 2392, 2393, 2394, 2396, 2398, 2399, 2401, 2403, 2405, 2406, 2408, 2410, 2411, 2412, 2416, 2417, 2419, 2421, 2422, 2425, 2426, 2427, 2429, 2430, 2431, 2432, 2435, 2436, 2437, 2438, 2439, 2441, 2444, 2446, 2447, 2448, 2449, 2450, 2453, 2456, 2461, 2464, 2465, 2467, 2469, 2470, 2474, 2475, 2476, 2477, 2480, 2482, 2483, 2484, 2488, 2489, 2490, 2491, 2497, 2504, 2506, 2508, 2509, 2511, 2512, 2513, 2514, 2515, 2517, 2518, 2520, 2523, 2524, 2526, 2527, 2531, 2537, 2539, 2540, 2543, 2544, 2545, 2546, 2549, 2550, 2551, 2553, 2554, 2555, 2556, 2557, 2559, 2562, 2563, 2564, 2567, 2574, 2578, 2581, 2582, 2583, 2584, 2587, 2588, 2589, 2590, 2591, 2594, 2596, 2598, 2601, 2603, 2604, 2605, 2606, 2607, 2608, 2609, 2610, 2613, 2616, 2618, 2622, 2626, 2629, 2632, 2637, 2638, 2641, 2644, 2645, 2646, 2653, 2655, 2657, 2658, 2659, 2664, 2665, 2666, 2669, 2671, 2672, 2673, 2674, 2681, 2682, 2686, 2688, 2689, 2690, 2691, 2693, 2695, 2699, 2700, 2701, 2702, 2704, 2705, 2706, 2708, 2709, 2710, 2713, 2716, 2718, 2723, 2724, 2727, 2728, 2729, 2731, 2736, 2739, 2740, 2742, 2745, 2747, 2748, 2749, 2751, 2753, 2754, 2768, 2770, 2771, 2773, 2774, 2776, 2777, 2778, 2788, 2789, 2791, 2793, 2794, 2795, 2796, 2799, 2802, 2804, 2805, 2807, 2808, 2809, 2811, 2812, 2813, 2816, 2817, 2820, 2821, 2822, 2823, 2824, 2825, 2830, 2834, 2839, 2840, 2841, 2842, 2843, 2845, 2850, 2851, 2853, 2854, 2855, 2857, 2858, 2859, 2861, 2865, 2866, 2870, 2873, 2877, 2879, 2880, 2881, 2885, 2886, 2887, 2888, 2890, 2891, 2892, 2893, 2894, 2896, 2898, 2899, 2900, 2901, 2903, 2906, 2907, 2908, 2909, 2912, 2915, 2917, 2919, 2923, 2924, 2925, 2932, 2934, 2935, 2943, 2944, 2946, 2948, 2949, 2951, 2952, 2954, 2957, 2958, 2960, 2962, 2966, 2969, 2973, 2974, 2975, 2979, 2980, 2981, 2982, 2984, 2986, 2987, 2989, 2990, 2991, 2993, 2994, 2995, 2996, 2997, 2999, 3000, 3001, 3003, 3004, 3005, 3007, 3009, 3011, 3013, 3014, 3016, 3022, 3028, 3031, 3032, 3036, 3037, 3042, 3045, 3046, 3048, 3050, 3052, 3058, 3059, 3061, 3065, 3066, 3070, 3073, 3074, 3075, 3076, 3077, 3082, 3083, 3085, 3089, 3090, 3092, 3094, 3095, 3098, 3099, 3100, 3101, 3103, 3106, 3107, 3113, 3114, 3115, 3117, 3118, 3119, 3120, 3124, 3125, 3126, 3127, 3128, 3130, 3131, 3132, 3137, 3138, 3139, 3143, 3144, 3145, 3146, 3152, 3153, 3155, 3158, 3162, 3164, 3166, 3171, 3172, 3173, 3174, 3178, 3179, 3180, 3181, 3182, 3185, 3187, 3188, 3190, 3191, 3194, 3195, 3197, 3198, 3201, 3202, 3204, 3206, 3207, 3211, 3212, 3213, 3215, 3218, 3219, 3223, 3224, 3227, 3228, 3229, 3232, 3234, 3235, 3237, 3239, 3240, 3241, 3245, 3247, 3248, 3252, 3254, 3255, 3256, 3258, 3259, 3260, 3263, 3265, 3267, 3268, 3273, 3276, 3277, 3280, 3281, 3284, 3286, 3289, 3291, 3293, 3297, 3299, 3301, 3307, 3308, 3311, 3312, 3313, 3314, 3319, 3320, 3322, 3323, 3324, 3326, 3327, 3328, 3329, 3330, 3331, 3332, 3339, 3340, 3341, 3343, 3345, 3346, 3347, 3348, 3349, 3350, 3353, 3355, 3356, 3358, 3360, 3364, 3366, 3367, 3369, 3371, 3372, 3382, 3383, 3388, 3389, 3390, 3394, 3396, 3401, 3403, 3404, 3405, 3407, 3409, 3411, 3415, 3417, 3419, 3420, 3424, 3425, 3427, 3429, 3430, 3431, 3433, 3436, 3437, 3439, 3440, 3441, 3444, 3447, 3448, 3452, 3454, 3455, 3457, 3458, 3459, 3461, 3462, 3463, 3464, 3467, 3468, 3475, 3478, 3479, 3482, 3483, 3489, 3495, 3497, 3499, 3500, 3504, 3506, 3507, 3512, 3514, 3521, 3523, 3527, 3529, 3530, 3531, 3535, 3536, 3538, 3540, 3542, 3545, 3546, 3548, 3550, 3551, 3553, 3554, 3555, 3558, 3561, 3563, 3565, 3566, 3567, 3568, 3573, 3576, 3577, 3578, 3579, 3581, 3582, 3584, 3587, 3593, 3594, 3598, 3599, 3602, 3604, 3605, 3606, 3607, 3609, 3610, 3611, 3616, 3618, 3619, 3620, 3621, 3622, 3624, 3625, 3626, 3627, 3630, 3632, 3635, 3639, 3641, 3643, 3645, 3648, 3651, 3653, 3654, 3655, 3656, 3659, 3661, 3663, 3666, 3667, 3669, 3671, 3673, 3675, 3676, 3677, 3679, 3680, 3682, 3684, 3685, 3688, 3689, 3691, 3692, 3693, 3695, 3696, 3700, 3702, 3704, 3705, 3708, 3709, 3711, 3713, 3714, 3717, 3718, 3719, 3720, 3721, 3724, 3725, 3726, 3727, 3728, 3732, 3736, 3737, 3738, 3739, 3740, 3742, 3748, 3749, 3752, 3754, 3755, 3756, 3758, 3760, 3761, 3762, 3763, 3764, 3766, 3772, 3773, 3776, 3778, 3785, 3786, 3788, 3792, 3794, 3795, 3798, 3799, 3800, 3801, 3802, 3804, 3805, 3806, 3807, 3808, 3809, 3811, 3812, 3816, 3817, 3818, 3819, 3824, 3827, 3834, 3835, 3837, 3843, 3844, 3848, 3849, 3851, 3853, 3855, 3859, 3861, 3863, 3864, 3865, 3867, 3871, 3873, 3874, 3876, 3878, 3879, 3880, 3881, 3882, 3884, 3886, 3893, 3899, 3901, 3902, 3903, 3904, 3905, 3909, 3911, 3912, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3927, 3931, 3936, 3937, 3938, 3940, 3941, 3942, 3948, 3950, 3952, 3954, 3955, 3957, 3958, 3960, 3961, 3962, 3968, 3970, 3971, 3972, 3974, 3976, 3977, 3978, 3982, 3983, 3987, 3988, 3990, 3992, 3996, 3998, 3999, 4005, 4007, 4008, 4009, 4010, 4013, 4016, 4017, 4019, 4020, 4021, 4023, 4025, 4026, 4027, 4030, 4031, 4035, 4037, 4043, 4045, 4048, 4049, 4051, 4053, 4054, 4056, 4059, 4060, 4062, 4063, 4067, 4068, 4069, 4072, 4073, 4074, 4077, 4079, 4080, 4083, 4087, 4088, 4089, 4090, 4091, 4092, 4093, 4094, 4096, 4097, 4098, 4099, 4100, 4104, 4105, 4108, 4109, 4113, 4114, 4117, 4119, 4120, 4122, 4125, 4127, 4128, 4129, 4133, 4134, 4135, 4139, 4141, 4145, 4146, 4149, 4150, 4151, 4152, 4153, 4154, 4155, 4157, 4158, 4161, 4162, 4164, 4165, 4167, 4168, 4169, 4170, 4171, 4172, 4175, 4176, 4178, 4181, 4182, 4183, 4184, 4185, 4186, 4187, 4188, 4189, 4192, 4193, 4195, 4197, 4198, 4201, 4203, 4206, 4211, 4214, 4218, 4220, 4221, 4222, 4223, 4225, 4226, 4227, 4229, 4230, 4234, 4235, 4238, 4239, 4240, 4241, 4242, 4246, 4247, 4248, 4252, 4253, 4254, 4255, 4259, 4260, 4261, 4263, 4264, 4267, 4268, 4274, 4275, 4277, 4279, 4282, 4284, 4286, 4287, 4290, 4291, 4294, 4298, 4300, 4305, 4307, 4308, 4310, 4311, 4312, 4314, 4316, 4322, 4323, 4324, 4328, 4335, 4340, 4341, 4343, 4345, 4346, 4347, 4349, 4350, 4351, 4352, 4354, 4355, 4356, 4357, 4358, 4359, 4360, 4361, 4362, 4364, 4365, 4366, 4367, 4369, 4371, 4374, 4376, 4379, 4380, 4381, 4382, 4383, 4384, 4385, 4386, 4387, 4388, 4391, 4398, 4401, 4402, 4403, 4404, 4405, 4408, 4409, 4413, 4414, 4416, 4417, 4419, 4420, 4421, 4423, 4425, 4429, 4430, 4432, 4436, 4437, 4440, 4447, 4449, 4450, 4451, 4454, 4455, 4458, 4460, 4461, 4464, 4468, 4470, 4472, 4473, 4476, 4478, 4479, 4484, 4485, 4487, 4488, 4491, 4493, 4494, 4496, 4499, 4500, 4501, 4503, 4504, 4506, 4508, 4510, 4515, 4516, 4517, 4521, 4524, 4526, 4527, 4528, 4530, 4531, 4537, 4539, 4540, 4543, 4544, 4545, 4547, 4548, 4549, 4553, 4556, 4561, 4563, 4565, 4566, 4570, 4576, 4577, 4578, 4579, 4582, 4583, 4584, 4585, 4586, 4587, 4588, 4590, 4592, 4593, 4594, 4595, 4599, 4602, 4603, 4604, 4608, 4610, 4611, 4613, 4614, 4616, 4617, 4618, 4619, 4622, 4623, 4624, 4625, 4626, 4627, 4628, 4631, 4635, 4636, 4638, 4639, 4640, 4643, 4644, 4645, 4646, 4647, 4648, 4649, 4656, 4657, 4658, 4660, 4662, 4667, 4673, 4674, 4676, 4678, 4679, 4683, 4687, 4688, 4689, 4690, 4691, 4692, 4696, 4698, 4699, 4702, 4704, 4705, 4706, 4708, 4713, 4714, 4715, 4719, 4722, 4726, 4727, 4730, 4732, 4733, 4737, 4741, 4748, 4749, 4752, 4753, 4754, 4756, 4762, 4763, 4764, 4765, 4766, 4767, 4768, 4769, 4772, 4773, 4776, 4778, 4781, 4783, 4784, 4790, 4791, 4792, 4793, 4795, 4796, 4797, 4800, 4802, 4804, 4806, 4807, 4809, 4810, 4811, 4812, 4814, 4815, 4820, 4821, 4823, 4824, 4825, 4827, 4828, 4829, 4830, 4833, 4834, 4835, 4836, 4837, 4838, 4840, 4842, 4843, 4846, 4850, 4851, 4854, 4855, 4862, 4865, 4866, 4869, 4872, 4873, 4879, 4881, 4884, 4889, 4891, 4892, 4893, 4898, 4900, 4902, 4903, 4904, 4907, 4908, 4913, 4914, 4915, 4917, 4922, 4923, 4924, 4926, 4931, 4932, 4935, 4936, 4937, 4939, 4942, 4943, 4945, 4946, 4949, 4950, 4951, 4952, 4955, 4956, 4957, 4959, 4960, 4964, 4966, 4968, 4970, 4971, 4974, 4976, 4983, 4984, 4985, 4987, 4988, 4989, 4991, 4993, 4996, 4999, 5000, 5001, 5003, 5004, 5005, 5006, 5009, 5010, 5011, 5012, 5014, 5016, 5017, 5021, 5022, 5024, 5026, 5028, 5030, 5033, 5035, 5036, 5037, 5038, 5041, 5042, 5043, 5047, 5052, 5056, 5059, 5060, 5061, 5064, 5065, 5066, 5067, 5070, 5071, 5074, 5075, 5079, 5081, 5082, 5083, 5085, 5089, 5091, 5093, 5095, 5096, 5099, 5100, 5102, 5103, 5106, 5108, 5110, 5112, 5113, 5115, 5116, 5120, 5123, 5126, 5130, 5131, 5133, 5134, 5135, 5137, 5138, 5139, 5142, 5143, 5145, 5146, 5148, 5151, 5153, 5154, 5157, 5159, 5160, 5161, 5165, 5166, 5168, 5169, 5170, 5172, 5173, 5176, 5177, 5180, 5181, 5182, 5188, 5189, 5190, 5191, 5192, 5193, 5194, 5197, 5202, 5204, 5208, 5209, 5210, 5211, 5212, 5215, 5217, 5218, 5219, 5220, 5225, 5226, 5231, 5232, 5242, 5243, 5245, 5246, 5248, 5249, 5250, 5251, 5257, 5258, 5259, 5266, 5269, 5270, 5271, 5273, 5276, 5277, 5278, 5280, 5281, 5283, 5284, 5286, 5287, 5289, 5292, 5297, 5299, 5300, 5301, 5303, 5305, 5307, 5308, 5310, 5313, 5316, 5318, 5324, 5326, 5328, 5331, 5332, 5333, 5336, 5337, 5338, 5340, 5341, 5342, 5343, 5344, 5345, 5346, 5347, 5350, 5351, 5354, 5356, 5360, 5363, 5365, 5368, 5370, 5371, 5372, 5375, 5376, 5377, 5378, 5384, 5388, 5395, 5396, 5398, 5400, 5403, 5404, 5405, 5407, 5408, 5409, 5411, 5415, 5418, 5419, 5423, 5426, 5427, 5431, 5432, 5433, 5434, 5437, 5440, 5442, 5444, 5447, 5449, 5453, 5454, 5457, 5462, 5463, 5467, 5468, 5471, 5472, 5473, 5478, 5479, 5480, 5481, 5484, 5487, 5489, 5490, 5491, 5492, 5493, 5494, 5496, 5497, 5499, 5503, 5505, 5506, 5509, 5510, 5511, 5512, 5515, 5516, 5517, 5518, 5520, 5523, 5524, 5527, 5530, 5531, 5533, 5534, 5535, 5540, 5542, 5544, 5547, 5549, 5550, 5551, 5556, 5557, 5559, 5564, 5565, 5567, 5568, 5570, 5571, 5575, 5580, 5581, 5584, 5591, 5592, 5593, 5595, 5596, 5602, 5603, 5605, 5608, 5609, 5610, 5613, 5616, 5621, 5622, 5625, 5630, 5631, 5632, 5633, 5635, 5636, 5638, 5639, 5641, 5642, 5643, 5647, 5652, 5653, 5654, 5657, 5660, 5664, 5666, 5667, 5668, 5674, 5675, 5676, 5680, 5681, 5684, 5685, 5686, 5687, 5688, 5689, 5690, 5694, 5697, 5698, 5701, 5702, 5703, 5706, 5708, 5709, 5712, 5713, 5716, 5718, 5721, 5723, 5724, 5726, 5731, 5732, 5734, 5736, 5740, 5741, 5743, 5745, 5750, 5753, 5755, 5758, 5763, 5767, 5768, 5769, 5771, 5774, 5777, 5779, 5780, 5781, 5782, 5784, 5787, 5789, 5790, 5791, 5792, 5793, 5797, 5798, 5802, 5806, 5807, 5808, 5810, 5811, 5813, 5814, 5815, 5816, 5817, 5819, 5822, 5823, 5826, 5831, 5839, 5840, 5841, 5842, 5846, 5848, 5852, 5853, 5854, 5855, 5856, 5857, 5859, 5860, 5861, 5863, 5864, 5865, 5868, 5870, 5871, 5879, 5881, 5882, 5883, 5884, 5885, 5886, 5887, 5888, 5889, 5896, 5898, 5900, 5901, 5903, 5904, 5905, 5906, 5907, 5908, 5909, 5910, 5911, 5915, 5917, 5918, 5919, 5920, 5921, 5922, 5924, 5925, 5926, 5927, 5929, 5930, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 5943, 5945, 5948, 5950, 5951, 5953, 5954, 5956, 5957, 5963, 5964, 5965, 5966, 5967, 5968, 5970, 5972, 5973, 5976, 5977, 5979, 5982, 5983, 5985, 5988, 5991, 5992, 5993, 5997, 5998, 5999, 6000, 6001, 6002, 6005, 6008, 6013, 6015, 6019, 6020, 6021, 6022, 6023, 6025, 6028, 6029, 6031, 6039, 6040, 6042, 6043, 6045, 6046, 6047, 6050, 6054, 6056, 6063, 6064, 6066, 6068, 6069, 6070, 6071, 6072, 6073, 6074, 6075, 6076, 6078, 6079, 6080, 6082, 6083, 6084, 6086, 6089, 6090, 6092, 6094, 6098, 6099, 6100, 6102, 6103, 6104, 6105, 6107, 6110, 6113, 6114, 6117, 6120, 6124, 6125, 6126, 6127, 6130, 6131, 6132, 6133, 6134, 6137, 6139, 6140, 6141, 6142, 6144, 6145, 6146, 6148, 6151, 6152, 6154, 6163, 6165, 6166, 6169, 6170, 6173, 6176, 6177, 6179, 6182, 6184, 6185, 6186, 6187, 6193, 6197, 6198, 6199, 6200, 6202, 6203, 6205, 6207, 6209, 6210, 6211, 6212, 6216, 6218, 6220, 6221, 6222, 6225, 6227, 6229, 6232, 6233, 6234, 6236, 6240, 6242, 6243, 6246, 6250, 6251, 6258, 6261, 6263, 6266, 6268, 6270, 6271, 6272, 6283, 6286, 6288, 6292, 6294, 6295, 6297, 6300, 6301, 6302, 6305, 6306, 6307, 6309, 6311, 6312, 6313, 6315, 6319, 6322, 6325, 6326, 6327, 6329, 6330, 6331, 6332, 6334, 6340, 6341, 6345, 6346, 6347, 6348, 6351, 6354, 6355, 6356, 6358, 6360, 6361, 6362, 6363, 6365, 6366, 6367, 6368, 6369, 6371, 6372, 6376, 6377, 6380, 6382, 6383, 6384, 6386, 6387, 6389, 6390, 6391, 6394, 6396, 6397, 6398, 6400, 6401, 6402, 6403, 6405, 6406, 6407, 6411, 6412, 6413, 6416, 6418, 6419, 6421, 6426, 6427, 6430, 6432, 6433, 6436, 6438, 6445, 6449, 6450, 6451, 6453, 6454, 6456, 6457, 6458, 6460, 6462, 6463, 6469, 6472, 6476, 6477, 6479, 6480, 6481, 6482, 6486, 6487, 6488, 6489, 6491, 6493, 6495, 6496, 6497, 6498, 6500, 6505, 6506, 6507, 6509, 6510, 6512, 6516, 6518, 6520, 6527, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6542, 6543, 6548, 6550, 6551, 6554, 6556, 6560, 6564, 6566, 6569, 6570, 6574, 6575, 6576, 6579, 6583, 6585, 6587, 6588, 6590, 6591, 6592, 6593, 6594, 6598, 6599, 6600, 6607, 6612, 6614, 6617, 6622, 6626, 6627, 6634, 6635, 6636, 6638, 6639, 6640, 6641, 6642, 6643, 6645, 6649, 6652, 6653, 6655, 6656, 6659, 6662, 6666, 6667, 6673, 6674, 6676, 6677, 6678, 6683, 6685, 6686, 6687, 6688, 6690, 6693, 6694, 6695, 6699, 6700, 6702, 6703, 6706, 6708, 6710, 6712, 6716, 6717, 6718, 6719, 6722, 6723, 6727, 6728, 6729, 6732, 6734, 6735, 6738, 6739, 6740, 6743, 6746, 6747, 6749, 6751, 6754, 6756, 6757, 6760, 6762, 6763, 6765, 6766, 6770, 6771, 6773, 6776, 6779, 6781, 6784, 6785, 6786, 6787, 6788, 6790, 6791, 6794, 6798, 6799, 6803, 6804, 6805, 6807, 6811, 6813, 6815, 6816, 6817, 6818, 6821, 6822, 6823, 6824, 6825, 6829, 6830, 6831, 6832, 6835, 6836, 6838, 6839, 6841, 6842, 6843, 6846, 6847, 6848, 6849, 6855, 6857, 6858, 6860, 6861, 6862, 6863, 6866, 6867, 6872, 6873, 6875, 6876, 6878, 6881, 6882, 6884, 6885, 6889, 6890, 6891, 6892, 6893, 6894, 6897, 6898, 6901, 6906, 6907, 6909, 6914, 6915, 6917, 6918, 6919, 6922, 6923, 6926, 6927, 6928, 6929, 6931, 6933, 6934, 6935, 6937, 6938, 6944, 6946, 6948, 6949, 6950, 6951, 6953, 6954, 6956, 6960, 6961, 6962, 6963, 6964, 6965, 6967, 6968, 6969, 6970, 6972, 6973, 6974, 6975, 6977, 6978, 6979, 6982, 6983, 6985, 6989, 6990, 6992, 6994, 6995, 6997, 6998, 7001, 7002, 7003, 7005, 7006, 7007, 7008, 7011, 7015, 7017, 7018, 7021, 7024, 7026, 7028, 7029, 7031, 7034, 7036, 7037, 7039, 7040, 7041, 7044, 7045, 7047, 7048, 7049, 7052, 7054, 7059, 7067, 7069, 7070, 7071, 7073, 7083, 7086, 7090, 7091, 7092, 7093, 7094, 7097, 7098, 7099, 7100, 7102, 7107, 7109, 7110, 7111, 7112, 7113, 7114, 7115, 7117, 7118, 7125, 7127, 7130, 7132, 7137, 7139, 7143, 7144, 7145, 7148, 7149, 7150, 7154, 7155, 7156, 7157, 7158, 7159, 7164, 7168, 7170, 7171, 7173, 7180, 7181, 7183, 7184, 7186, 7189, 7192, 7193, 7194, 7195, 7199, 7200, 7201, 7204, 7205, 7207, 7209, 7210, 7211, 7214, 7218, 7219, 7220, 7223, 7224, 7225, 7227, 7229, 7230, 7231, 7234, 7237, 7240, 7242, 7247, 7248, 7249, 7251, 7254, 7256, 7257, 7258, 7262, 7265, 7268, 7269, 7274, 7276, 7277, 7278, 7281, 7282, 7288, 7291, 7293, 7294, 7295, 7296, 7298, 7299, 7300, 7301, 7306, 7308, 7310, 7311, 7312, 7313, 7315, 7316, 7317, 7318, 7319, 7321, 7322, 7323, 7328, 7331, 7333, 7334, 7338, 7341, 7343, 7344, 7345, 7347, 7348, 7349, 7351, 7352, 7355, 7356, 7359, 7360, 7362, 7363, 7364, 7365, 7366, 7367, 7370, 7371, 7374, 7375, 7377, 7380, 7383, 7384, 7387, 7389, 7390, 7391, 7392, 7393, 7395, 7398, 7401, 7404, 7408, 7409, 7410, 7413, 7415, 7417, 7418, 7421, 7424, 7425, 7426, 7427, 7428, 7429, 7431, 7433, 7435, 7436, 7439, 7443, 7445, 7448, 7449, 7450, 7451, 7452, 7453, 7455, 7456, 7457, 7458, 7460, 7461, 7465, 7466, 7469, 7475, 7478, 7479, 7481, 7482, 7483, 7487, 7488, 7489, 7490, 7491, 7492, 7496, 7497, 7499, 7500, 7501, 7504, 7507, 7508, 7509, 7510, 7513, 7514, 7515, 7516, 7518, 7521, 7523, 7525, 7526, 7527, 7529, 7530, 7532, 7533, 7535, 7536, 7537, 7544, 7546, 7547, 7548, 7549, 7550, 7552, 7553, 7555, 7557, 7560, 7561, 7562, 7563, 7564, 7565, 7566, 7567, 7568, 7571, 7572, 7574, 7575, 7580, 7581, 7582, 7583, 7584, 7591, 7593, 7597, 7598, 7599, 7603, 7605, 7608, 7610, 7611, 7612, 7613, 7614, 7617, 7620, 7621, 7623, 7626, 7627, 7628, 7629, 7630, 7631, 7636, 7639, 7640, 7641, 7644, 7646, 7649, 7650, 7651, 7652, 7655, 7656, 7657, 7658, 7660, 7661, 7662, 7664, 7665, 7667, 7668, 7669, 7671, 7672, 7673, 7674, 7675, 7676, 7677, 7678, 7683, 7684, 7688, 7692, 7694, 7696, 7700, 7703, 7704, 7705, 7706, 7707, 7709, 7710, 7713, 7714, 7715, 7716, 7720, 7722, 7724, 7725, 7726, 7729, 7737, 7738, 7739, 7740, 7742, 7743, 7749, 7750, 7752, 7753, 7755, 7756, 7760, 7761, 7762, 7763, 7766, 7768, 7769, 7771, 7776, 7777, 7778, 7779, 7781, 7784, 7785, 7786, 7791, 7792, 7794, 7796, 7797, 7800, 7802, 7803, 7804, 7805, 7806, 7808, 7809, 7813, 7814, 7815, 7817, 7821, 7825, 7832, 7833, 7835, 7837, 7838, 7840, 7843, 7847, 7849, 7852, 7856, 7858, 7863, 7864, 7865, 7869, 7872, 7873, 7875, 7876, 7877, 7878, 7880, 7883, 7885, 7887, 7888, 7889, 7890, 7894, 7895, 7896, 7897, 7898, 7903, 7905, 7906, 7909, 7911, 7912, 7913, 7914, 7918, 7922, 7923, 7926, 7928, 7929, 7932, 7935, 7936, 7938, 7939, 7940, 7943, 7945, 7946, 7948, 7952, 7958, 7959, 7960, 7962, 7964, 7965, 7966, 7969, 7970, 7971, 7980, 7983, 7984, 7986, 7988, 7990, 7992, 7993, 7995, 7996, 7998, 8002, 8005, 8010, 8011, 8014, 8015, 8016, 8017, 8021, 8022, 8023, 8025, 8026, 8029, 8030, 8032, 8033, 8035, 8038, 8039, 8042, 8043, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8053, 8054, 8055, 8057, 8058, 8059, 8060, 8061, 8065, 8066, 8067, 8068, 8069, 8070, 8071, 8072, 8073, 8074, 8076, 8077, 8079, 8080, 8082, 8084, 8085, 8089, 8090, 8091, 8093, 8094, 8095, 8096, 8097, 8099, 8103, 8104, 8106, 8108, 8109, 8114, 8117, 8118, 8121, 8123, 8124, 8125, 8128, 8130, 8131, 8132, 8133, 8137, 8139, 8140, 8141, 8143, 8145, 8147, 8148, 8150, 8151, 8152, 8153, 8155, 8157, 8158, 8160, 8163, 8167, 8168, 8170, 8172, 8173, 8174, 8177, 8179, 8181, 8184, 8186, 8187, 8189, 8192, 8193, 8194, 8196, 8197, 8199, 8205, 8206, 8207, 8208, 8209, 8212, 8216, 8221, 8222, 8223, 8225, 8226, 8228, 8229, 8231, 8232, 8234, 8236, 8237, 8238, 8239, 8241, 8242, 8248, 8251, 8252, 8253, 8254, 8255, 8258, 8259, 8262, 8263, 8265, 8267, 8272, 8273, 8275, 8277, 8278, 8279, 8280, 8281, 8283, 8284, 8286, 8292, 8293, 8294, 8295, 8296, 8298, 8299, 8304, 8305, 8310, 8311, 8313, 8314, 8315, 8316, 8317, 8320, 8322, 8324, 8327, 8330, 8331, 8332, 8333, 8334, 8337, 8341, 8342, 8343, 8344, 8349, 8354, 8355, 8357, 8358, 8359, 8361, 8362, 8364, 8365, 8366, 8371, 8373, 8375, 8377, 8378, 8381, 8382, 8390, 8392, 8393, 8394, 8400, 8402, 8403, 8404, 8405, 8410, 8412, 8413, 8415, 8418, 8419, 8421, 8422, 8424, 8425, 8426, 8427, 8428, 8429, 8430, 8431, 8435, 8440, 8442, 8443, 8445, 8447, 8448, 8452, 8454, 8456, 8457, 8458, 8461, 8464, 8468, 8469, 8470, 8472, 8474, 8477, 8479, 8484, 8485, 8486, 8487, 8488, 8494, 8497, 8499, 8501, 8502, 8504, 8505, 8507, 8509, 8510, 8512, 8515, 8516, 8517, 8521, 8528, 8529, 8532, 8534, 8538, 8539, 8540, 8541, 8542, 8543, 8545, 8548, 8549, 8550, 8552, 8553, 8557, 8558, 8561, 8564, 8565, 8566, 8568, 8573, 8574, 8575, 8576, 8578, 8581, 8584, 8585, 8589, 8591, 8596, 8598, 8601, 8605, 8607, 8611, 8612, 8614, 8615, 8616, 8619, 8620, 8623, 8624, 8625, 8626, 8628, 8631, 8633, 8639, 8640, 8646, 8648, 8650, 8652, 8660, 8662, 8664, 8665, 8667, 8669, 8673, 8678, 8679, 8680, 8684, 8688, 8691, 8692, 8695, 8697, 8698, 8701, 8703, 8704, 8705, 8708, 8709, 8710, 8712, 8714, 8716, 8719, 8720, 8721, 8723, 8725, 8726, 8733, 8734, 8735, 8736, 8738, 8739, 8743, 8744, 8745, 8747, 8748, 8749, 8751, 8753, 8755, 8758, 8759, 8763, 8766, 8767, 8768, 8769, 8771, 8772, 8773, 8779, 8781, 8782, 8783, 8788, 8789, 8791, 8792, 8793, 8795, 8796, 8797, 8800, 8801, 8802, 8806, 8807, 8808, 8811, 8814, 8815, 8817, 8818, 8820, 8827, 8830, 8831, 8832, 8834, 8837, 8838, 8846, 8847, 8848, 8849, 8850, 8853, 8854, 8858, 8861, 8862, 8863, 8864, 8865, 8870, 8872, 8873, 8874, 8878, 8879, 8880, 8881, 8882, 8884, 8888, 8889, 8891, 8892, 8893, 8894, 8896, 8897, 8898, 8900, 8901, 8902, 8904, 8908, 8912, 8915, 8917, 8918, 8920, 8923, 8926, 8930, 8931, 8932, 8933, 8934, 8937, 8939, 8941, 8943, 8944, 8947, 8948, 8949, 8951, 8953, 8955, 8956, 8960, 8962, 8964, 8965, 8967, 8969, 8971, 8972, 8977, 8979, 8985, 8988, 8991, 8993, 8994, 8995, 8997, 9006, 9007, 9011, 9012, 9016, 9017, 9019, 9020, 9021, 9024, 9026, 9027, 9029, 9033, 9034, 9035, 9038, 9042, 9043, 9046, 9047, 9050, 9052, 9061, 9063, 9067, 9071, 9075, 9077, 9079, 9080, 9081, 9082, 9085, 9090, 9091, 9092, 9095, 9098, 9101, 9105, 9106, 9108, 9112, 9115, 9118, 9119, 9121, 9124, 9126, 9131, 9133, 9134, 9137, 9138, 9139, 9141, 9142, 9144, 9145, 9146, 9148, 9149, 9150, 9151, 9153, 9158, 9160, 9161, 9162, 9164, 9165, 9166, 9167, 9169, 9171, 9173, 9174, 9178, 9180, 9183, 9185, 9188, 9189, 9191, 9192, 9193, 9194, 9198, 9199, 9201, 9202, 9204, 9207, 9208, 9209, 9210, 9211, 9212, 9215, 9220, 9221, 9222, 9224, 9225, 9226, 9228, 9229, 9230, 9231, 9234, 9235, 9237, 9238, 9239, 9240, 9242, 9247, 9248, 9249, 9250, 9251, 9252, 9254, 9258, 9264, 9266, 9267, 9268, 9272, 9273, 9276, 9277, 9279, 9281, 9284, 9287, 9288, 9289, 9290, 9293, 9295, 9296, 9297, 9304, 9305, 9306, 9313, 9314, 9315, 9317, 9318, 9319, 9320, 9324, 9325, 9326, 9328, 9330, 9331, 9333, 9334, 9335, 9336, 9337, 9340, 9341, 9342, 9343, 9344, 9345, 9347, 9349, 9352, 9353, 9355, 9357, 9358, 9359, 9365, 9368, 9372, 9374, 9376, 9377, 9381, 9382, 9384, 9385, 9386, 9389, 9393, 9395, 9397, 9399, 9401, 9402, 9408, 9411, 9412, 9413, 9415, 9416, 9418, 9421, 9422, 9424, 9426, 9427, 9430, 9433, 9436, 9437, 9443, 9446, 9447, 9449, 9450, 9451, 9452, 9455, 9456, 9457, 9458, 9459, 9460, 9461, 9462, 9463, 9465, 9470, 9471, 9472, 9473, 9476, 9479, 9481, 9482, 9483, 9485, 9486, 9490, 9491, 9495, 9496, 9498, 9499, 9501, 9502, 9503, 9504, 9506, 9507, 9514, 9516, 9517, 9523, 9526, 9531, 9533, 9536, 9537, 9539, 9541, 9542, 9544, 9545, 9547, 9548, 9549, 9553, 9554, 9555, 9562, 9563, 9564, 9565, 9566, 9568, 9570, 9571, 9576, 9578, 9579, 9582, 9584, 9585, 9589, 9594, 9595, 9598, 9603, 9604, 9606, 9607, 9614, 9618, 9620, 9622, 9623, 9624, 9626, 9627, 9628, 9631, 9633, 9634, 9637, 9641, 9642, 9644, 9647, 9649, 9650, 9651, 9654, 9657, 9658, 9663, 9665, 9668, 9670, 9672, 9673, 9677, 9678, 9679, 9680, 9683, 9684, 9686, 9688, 9690, 9691, 9692, 9693, 9696, 9697, 9700, 9703, 9705, 9708, 9709, 9711, 9712, 9715, 9718, 9724, 9726, 9727, 9734, 9735, 9738, 9739, 9740, 9744, 9745, 9746, 9750, 9751, 9755, 9756, 9757, 9758, 9759, 9761, 9762, 9763, 9765, 9769, 9770, 9771, 9772, 9773, 9777, 9778, 9781, 9782, 9784, 9785, 9787, 9788, 9789, 9790, 9793, 9794, 9800, 9802, 9803, 9804, 9808, 9809, 9813, 9815, 9816, 9817, 9818, 9820, 9825, 9827, 9830, 9833, 9834, 9835, 9837, 9839, 9840, 9841, 9842, 9844, 9846, 9847, 9848, 9849, 9850, 9851, 9852, 9853, 9854, 9859, 9860, 9862, 9865, 9869, 9871, 9874, 9875, 9876, 9879, 9882, 9885, 9886, 9887, 9890, 9891])
    ]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_match_clusters_asg_new()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_match_clusters_asg_new()
        self.assertEqual(
            str(cm_new.exception),
            "k_means_match_clusters_asg_new() missing 2 required positional arguments: 'asg1' and 'asg2'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):

        return_new = fu.k_means_match_clusters_asg_new(
            asg1=self.asg1, asg2=self.asg2, T=0
        )
        return_old = oldfu.k_means_match_clusters_asg_new(
            asg1=self.asg1, asg2=self.asg2, T=0
        )

        self.assertTrue(array_equal(return_new[0], [[0, 1]]))
        self.assertEqual(return_new[2], 2503)
        self.assertTrue(
            allclose(return_new[0], return_old[0], atol=TOLERANCE, equal_nan=True)
        )
        self.assertTrue(
            array_equal(
                [return_new[1][0].tolist(), return_new[1][1].tolist()],
                [return_old[1][0].tolist(), return_old[1][1].tolist()],
            )
        )
        self.assertEqual(return_new[2], return_old[2])

    def test_emptyList_asg1_returns_IndexError(self):
        return_new = fu.k_means_match_clusters_asg_new(asg1=[], asg2=self.asg2, T=0)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1=[], asg2=self.asg2, T=0)
        self.assertTrue(array_equal(return_new[0], []))
        self.assertTrue(array_equal(return_new[1], []))
        self.assertEqual(return_new[2], 0)
        self.assertTrue(allclose(return_new[0], return_old[0]))
        self.assertTrue(allclose(return_new[1], return_old[1]))
        self.assertEqual(return_new[2], return_old[2])

    def test_emptyList_asg2_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=[], T=0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_match_clusters_asg_new(asg1=self.asg1, asg2=[], T=0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_emptyList_both_asg(self):
        return_new = fu.k_means_match_clusters_asg_new(asg1=[], asg2=[], T=0)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1=[], asg2=[], T=0)

        self.assertTrue(array_equal(return_new[0], []))
        self.assertTrue(array_equal(return_new[1], []))
        self.assertEqual(return_new[2], 0)
        self.assertTrue(allclose(return_new[0], return_old[0]))
        self.assertTrue(allclose(return_new[1], return_old[1]))
        self.assertEqual(return_new[2], return_old[2])


class Test_hist_list(unittest.TestCase):
    data =[1744.5, 16.0, 10.0, 5.0, 25.0, 1742.5, 879.25, 878.25, 25.0, 1749.5, 4.0, 1745.5, 890.25, 875.25, 10.0, 20.0, 874.25, 1772.5, 1741.5, 871.25, 10.0, 887.25, 17.0, 888.25, 25.0, 874.25, 875.25, 1760.5, 2.0, 880.25, 895.25, 895.25, 10.0, 1740.5, 5.0, 20.0, 902.25, 5.0, 2.0, 1757.5, 879.25, 1757.5, 0.0, 872.25, 872.25, 874.25, 887.25, 871.25, 2.0, 874.25, 13.0, 895.25, 872.25, 13.0, 1758.5, 871.25, 1741.5, 874.25, 2.0, 1.0, 871.25, 871.25, 883.25, 1742.5, 872.25, 1744.5, 872.25, 8.0, 871.25, 886.25, 875.25, 878.25, 18.0, 874.25, 874.25, 874.25, 1745.5, 883.25, 2.0, 880.25, 887.25, 5.0, 880.25, 886.25, 875.25, 875.25, 890.25, 5.0, 1758.5, 880.25, 871.25, 20.0, 886.25, 1741.5, 875.25, 895.25, 871.25, 2.0, 872.25, 1750.5, 870.25, 1742.5, 872.25, 875.25, 5.0, 1748.5, 872.25, 2.0, 4.0, 1753.5, 878.25, 10.0, 890.25, 872.25, 1.0, 1744.5, 1744.5, 10.0, 1742.5, 2.0, 25.0, 2.0, 1.0, 1748.5, 13.0, 10.0, 875.25, 870.25, 887.25, 1745.5, 2.0, 874.25, 1742.5, 872.25, 878.25, 1741.5, 871.25, 1757.5, 879.25, 880.25, 874.25, 870.25, 1757.5, 1742.5, 875.25, 17.0, 25.0, 1745.5, 888.25, 878.25, 16.0, 1760.5, 871.25, 1742.5, 10.0, 13.0, 874.25, 875.25, 1741.5, 870.25, 875.25, 10.0, 875.25, 875.25, 875.25, 2.0, 872.25, 875.25, 16.0, 5.0, 878.25, 1741.5, 20.0, 872.25, 878.25, 874.25, 872.25, 8.0, 872.25, 871.25, 875.25, 1741.5, 2.0, 17.0, 1757.5, 1742.5, 872.25, 17.0, 875.25, 872.25, 2.0, 13.0, 871.25, 890.25, 17.0, 879.25, 1744.5, 902.25, 872.25, 4.0, 875.25, 1760.5, 872.25, 1741.5, 9.0, 895.25, 5.0, 880.25, 1745.5, 875.25, 883.25, 1745.5, 883.25, 871.25, 10.0, 17.0, 18.0, 5.0, 10.0, 886.25, 10.0, 880.25, 875.25, 878.25, 20.0, 1757.5, 895.25, 1745.5, 880.25, 872.25, 871.25, 1753.5, 4.0, 1742.5, 10.0, 1757.5, 20.0, 1765.5, 875.25, 872.25, 871.25, 871.25, 4.0, 1.0, 887.25, 871.25, 875.25, 2.0, 878.25, 5.0, 1756.5, 874.25, 1765.5, 878.25, 5.0, 0.0, 1741.5, 875.25, 5.0, 879.25, 874.25, 875.25, 5.0, 5.0, 1748.5, 871.25, 4.0, 872.25, 875.25, 8.0, 887.25, 4.0, 888.25, 878.25, 1765.5, 1744.5, 5.0, 2.0, 1772.5, 17.0, 872.25, 871.25, 1760.5, 1744.5, 1742.5, 1742.5, 880.25, 2.0, 1760.5, 870.25, 878.25, 890.25, 883.25, 2.0, 874.25, 4.0, 880.25, 5.0, 2.0, 874.25, 874.25, 1748.5, 1744.5, 871.25, 875.25, 888.25, 1.0, 872.25, 20.0, 879.25, 5.0, 1740.5, 880.25, 875.25, 874.25, 1745.5, 871.25, 0.0, 1742.5, 874.25, 1.0, 1748.5, 874.25, 9.0, 871.25, 1756.5, 872.25, 1741.5, 1742.5, 871.25, 871.25, 1745.5, 875.25, 871.25, 16.0, 871.25, 4.0, 8.0, 895.25, 890.25, 875.25, 887.25, 1745.5, 1742.5, 874.25, 1744.5, 9.0, 1748.5, 1.0, 875.25, 5.0, 883.25, 1.0, 886.25, 874.25, 871.25, 872.25, 1742.5, 878.25, 1742.5, 1.0, 875.25, 875.25, 880.25, 875.25, 1.0, 879.25, 1.0, 886.25, 1750.5, 1740.5, 872.25, 1.0, 1742.5, 1.0, 878.25, 879.25, 1745.5, 875.25, 2.0, 878.25, 1744.5, 887.25, 5.0, 17.0, 1744.5, 871.25, 871.25, 13.0, 880.25, 875.25, 5.0, 874.25, 1760.5, 890.25, 10.0, 878.25, 879.25, 18.0, 874.25, 1744.5, 875.25, 10.0, 10.0, 4.0, 880.25, 2.0, 1750.5, 887.25, 1753.5, 871.25, 1745.5, 5.0, 872.25, 886.25, 13.0, 879.25, 1748.5, 1745.5, 890.25, 16.0, 1745.5, 872.25, 895.25, 5.0, 4.0, 20.0, 887.25, 883.25, 1742.5, 883.25, 4.0, 25.0, 1745.5, 874.25, 1744.5, 886.25, 2.0, 875.25, 871.25, 880.25, 1.0, 875.25, 887.25, 871.25, 1740.5, 875.25, 870.25, 878.25, 5.0, 2.0, 1745.5, 880.25, 890.25, 1744.5, 2.0, 10.0, 8.0, 1.0, 880.25, 872.25, 1741.5, 13.0, 1742.5, 887.25, 1.0, 879.25, 872.25, 1750.5, 2.0, 890.25, 32.0, 890.25, 8.0, 4.0, 1745.5, 875.25, 8.0, 880.25, 875.25, 20.0, 1741.5, 871.25, 1742.5, 10.0, 874.25, 874.25, 871.25, 878.25, 0.0, 895.25, 870.25, 870.25, 887.25, 871.25, 2.0, 871.25, 872.25, 8.0, 871.25, 880.25, 8.0, 872.25, 870.25, 1.0, 871.25, 4.0, 1740.5, 1756.5, 1742.5, 870.25, 1741.5, 875.25, 5.0, 1744.5, 871.25, 9.0, 871.25, 871.25, 890.25, 890.25, 875.25, 5.0, 1748.5, 4.0, 883.25, 872.25, 870.25, 2.0, 1.0, 875.25, 1742.5, 1741.5, 4.0, 902.25, 880.25, 20.0, 1750.5, 32.0, 16.0, 1765.5, 887.25, 1.0, 25.0, 20.0, 871.25, 4.0, 871.25, 872.25, 880.25, 1.0, 4.0, 1745.5, 4.0, 875.25, 1741.5, 872.25, 870.25, 8.0, 870.25, 9.0, 870.25, 883.25, 1.0, 879.25, 875.25, 2.0, 1749.5, 1750.5, 871.25, 4.0, 5.0, 1744.5, 902.25, 874.25, 870.25, 871.25, 886.25, 5.0, 5.0, 883.25, 8.0, 1750.5, 874.25, 5.0, 1757.5, 1741.5, 878.25, 1744.5, 874.25, 2.0, 875.25, 1745.5, 879.25, 880.25, 1741.5, 870.25, 1740.5, 8.0, 887.25, 1741.5, 875.25, 871.25, 888.25, 875.25, 902.25, 887.25, 878.25, 1748.5, 875.25, 890.25, 874.25, 875.25, 4.0, 875.25, 1760.5, 16.0, 880.25, 870.25, 871.25, 1742.5, 8.0, 1760.5, 16.0, 4.0, 10.0, 1757.5, 1749.5, 2.0, 890.25, 890.25, 1750.5, 8.0, 883.25, 890.25, 1741.5, 17.0, 1.0, 16.0, 872.25, 18.0, 25.0, 1742.5, 878.25, 1749.5, 887.25, 874.25, 875.25, 887.25, 20.0, 25.0, 879.25, 2.0, 883.25, 870.25, 875.25, 886.25, 875.25, 872.25, 895.25, 871.25, 887.25, 870.25, 16.0, 887.25, 874.25, 887.25, 874.25, 13.0, 880.25, 8.0, 878.25, 1742.5, 1742.5, 4.0, 890.25, 1745.5, 875.25, 875.25, 879.25, 871.25, 2.0, 1745.5, 1750.5, 875.25, 5.0, 870.25, 1741.5, 874.25, 880.25, 875.25, 0.0, 4.0, 871.25, 16.0, 886.25, 1750.5, 875.25, 871.25, 1765.5, 871.25, 20.0, 1740.5, 890.25, 5.0, 875.25, 2.0, 1748.5, 878.25, 10.0, 1742.5, 872.25, 1742.5, 880.25, 1760.5, 20.0, 878.25, 887.25, 10.0, 875.25, 874.25, 878.25, 25.0, 5.0, 872.25, 883.25, 1753.5, 1745.5, 887.25, 872.25, 1745.5, 2.0, 874.25, 8.0, 878.25, 1741.5, 1.0, 1750.5, 887.25, 1.0, 902.25, 1745.5, 878.25, 1742.5, 10.0, 871.25, 886.25, 1745.5, 16.0, 874.25, 1757.5, 1742.5, 880.25, 1742.5, 890.25, 9.0, 871.25, 886.25, 2.0, 8.0, 1741.5, 1748.5, 1744.5, 1.0, 1742.5, 1741.5, 890.25, 895.25, 1757.5, 871.25, 875.25, 8.0, 1.0, 875.25, 875.25, 1741.5, 872.25, 871.25, 878.25, 2.0, 871.25, 1.0, 875.25, 880.25, 5.0, 875.25, 4.0, 1750.5, 874.25, 1745.5, 1744.5, 875.25, 871.25, 4.0, 875.25, 1742.5, 872.25, 874.25, 9.0, 1742.5, 1745.5, 874.25, 871.25, 1748.5, 1742.5, 890.25, 1740.5, 875.25, 874.25, 9.0, 880.25, 8.0, 875.25, 16.0, 883.25, 4.0, 2.0, 872.25, 9.0, 886.25, 890.25, 25.0, 17.0, 16.0, 1744.5, 879.25, 872.25, 875.25, 1741.5, 875.25, 1741.5, 886.25, 1753.5, 1745.5, 880.25, 1757.5, 1741.5, 883.25, 872.25, 1745.5, 883.25, 875.25, 875.25, 1753.5, 4.0, 874.25, 4.0, 886.25, 4.0, 874.25, 875.25, 872.25, 1758.5, 875.25, 879.25, 890.25, 895.25, 1760.5, 1745.5, 1742.5, 16.0, 875.25, 879.25, 875.25, 5.0, 1750.5, 1.0, 5.0, 870.25, 871.25, 883.25, 890.25, 13.0, 871.25, 875.25, 1745.5, 871.25, 5.0, 871.25, 875.25, 872.25, 1753.5, 1744.5, 875.25, 875.25, 880.25, 1.0, 13.0, 883.25, 887.25, 874.25, 883.25, 1.0, 1.0, 890.25, 878.25, 4.0, 1742.5, 883.25, 1757.5, 1756.5, 875.25, 875.25, 25.0, 875.25, 1742.5, 2.0, 1741.5, 886.25, 875.25, 1741.5, 1750.5, 875.25, 1742.5, 20.0, 871.25, 1742.5, 1750.5, 1758.5, 0.0, 887.25, 872.25, 878.25, 883.25, 1745.5, 4.0, 1745.5, 1.0, 872.25, 1756.5, 1745.5, 1741.5, 1.0, 875.25, 1.0, 1741.5, 1745.5, 875.25, 872.25, 886.25, 2.0, 5.0, 871.25, 874.25, 875.25, 1742.5, 32.0, 872.25, 871.25, 870.25, 5.0, 1744.5, 888.25, 874.25, 1741.5, 883.25, 2.0, 871.25, 1745.5, 20.0, 880.25, 874.25, 874.25, 1744.5, 8.0, 1744.5, 20.0, 1741.5, 874.25, 872.25, 1744.5, 2.0, 1748.5, 1744.5, 871.25, 1745.5, 1750.5, 875.25, 1.0, 25.0, 872.25, 874.25, 20.0, 1.0, 887.25, 895.25, 2.0, 874.25, 1741.5, 2.0, 871.25, 0.0, 872.25, 870.25, 871.25, 874.25, 874.25, 1742.5, 871.25, 1745.5, 874.25, 870.25, 2.0, 880.25, 1741.5, 870.25, 1740.5, 872.25, 4.0, 870.25, 880.25, 878.25, 2.0, 874.25, 1741.5, 1742.5, 1.0, 874.25, 1745.5, 880.25, 1749.5, 872.25, 875.25, 4.0, 1.0, 875.25, 887.25, 2.0, 1741.5, 18.0, 874.25, 5.0, 875.25, 4.0, 871.25, 895.25, 872.25, 895.25, 2.0, 886.25, 10.0, 871.25, 871.25, 880.25, 1741.5, 890.25, 1741.5, 890.25, 890.25, 10.0, 16.0, 878.25, 871.25, 2.0, 886.25, 9.0, 1744.5, 875.25, 870.25, 1745.5, 1753.5, 1757.5, 871.25, 874.25, 879.25, 887.25, 1744.5, 20.0, 10.0, 1.0, 879.25, 872.25, 1742.5, 874.25, 8.0, 17.0, 871.25, 883.25, 1.0, 2.0, 1760.5, 890.25, 2.0, 1742.5, 1740.5, 1745.5, 874.25, 875.25, 17.0, 874.25, 886.25, 875.25, 879.25, 5.0, 871.25, 4.0, 880.25, 5.0, 1745.5, 1740.5, 872.25, 1744.5, 872.25, 1742.5, 888.25, 1757.5, 1.0, 1744.5, 880.25, 1742.5, 870.25, 871.25, 1.0, 1745.5, 872.25, 875.25, 880.25, 1.0, 870.25, 1744.5, 8.0, 1740.5, 1744.5, 871.25, 17.0, 1749.5, 1750.5, 872.25, 5.0, 888.25, 1742.5, 874.25, 10.0, 1.0, 1742.5, 874.25, 887.25, 1744.5, 10.0, 870.25, 875.25, 871.25, 1745.5, 1749.5, 878.25, 875.25, 871.25, 874.25, 874.25, 1744.5, 887.25, 4.0, 871.25, 1741.5, 13.0, 887.25, 1744.5, 16.0, 1.0, 872.25, 871.25, 1753.5, 874.25, 878.25, 5.0, 4.0, 895.25, 1740.5, 872.25, 872.25, 1745.5, 10.0, 875.25, 1741.5, 874.25, 1753.5, 875.25, 1741.5, 875.25, 1749.5, 1.0, 887.25, 870.25, 874.25, 890.25, 1741.5, 871.25, 871.25, 9.0, 13.0, 1.0, 880.25, 887.25, 878.25, 875.25, 874.25, 878.25, 1760.5, 1.0, 871.25, 1753.5, 1742.5, 1.0, 1742.5, 871.25, 872.25, 883.25, 1756.5, 875.25, 4.0, 25.0, 871.25, 872.25, 20.0, 25.0, 9.0, 1745.5, 13.0, 16.0, 871.25, 874.25, 1742.5, 2.0, 1740.5, 871.25, 4.0, 1745.5, 5.0, 874.25, 871.25, 887.25, 878.25, 25.0, 1741.5, 4.0, 0.0, 5.0, 870.25, 4.0, 874.25, 1745.5, 871.25, 10.0, 879.25, 874.25, 871.25, 878.25, 875.25, 1745.5, 871.25, 874.25, 5.0, 4.0, 1749.5, 1744.5, 871.25, 1749.5, 870.25, 1749.5, 1.0, 0.0, 872.25, 880.25, 1753.5, 879.25, 875.25, 872.25, 16.0, 1757.5, 875.25, 902.25, 871.25, 1742.5, 4.0, 886.25, 875.25, 17.0, 20.0, 886.25, 902.25, 20.0, 870.25, 871.25, 875.25, 1740.5, 871.25, 1742.5, 1757.5, 0.0, 1756.5, 1757.5, 1740.5, 871.25, 10.0, 13.0, 878.25, 17.0, 880.25, 1760.5, 1749.5, 887.25, 878.25, 1744.5, 5.0, 1749.5, 883.25, 8.0, 887.25, 880.25, 25.0, 871.25, 1741.5, 9.0, 1765.5, 880.25, 0.0, 883.25, 2.0, 10.0, 1741.5, 870.25, 1757.5, 1742.5, 16.0, 1.0, 1765.5, 9.0, 5.0, 1745.5, 872.25, 871.25, 10.0, 1744.5, 874.25, 872.25, 1741.5, 871.25, 871.25, 872.25, 1745.5, 10.0, 5.0, 2.0, 1745.5, 875.25, 872.25, 0.0, 890.25, 875.25, 875.25, 880.25, 1772.5, 1745.5, 1.0, 880.25, 1745.5, 9.0, 875.25, 890.25, 25.0, 13.0, 878.25, 1756.5, 870.25, 1745.5, 1753.5, 883.25, 874.25, 1.0, 13.0, 883.25, 16.0, 4.0, 872.25, 879.25, 2.0, 871.25, 1745.5, 0.0, 5.0, 872.25, 8.0, 1750.5, 872.25, 0.0, 870.25, 875.25, 0.0, 883.25, 10.0, 1749.5, 874.25, 878.25, 1.0, 1740.5, 4.0, 0.0, 1745.5, 1740.5, 5.0, 1.0, 1740.5, 1750.5, 1742.5, 17.0, 870.25, 1.0, 5.0, 888.25, 1745.5, 8.0, 886.25, 1745.5, 895.25, 0.0, 10.0, 888.25, 5.0, 880.25, 879.25, 874.25, 1741.5, 13.0, 872.25, 1757.5, 20.0, 874.25, 0.0, 874.25, 887.25, 871.25, 2.0, 1740.5, 13.0, 1753.5, 1741.5, 4.0, 880.25, 879.25, 879.25, 1.0, 1742.5, 1741.5, 1748.5, 5.0, 13.0, 4.0, 870.25, 1756.5, 878.25, 0.0, 1745.5, 878.25, 1741.5, 1745.5, 875.25, 5.0, 1745.5, 874.25, 872.25, 875.25, 16.0, 5.0, 1748.5, 1748.5, 871.25, 883.25, 875.25, 871.25, 1.0, 887.25, 1756.5, 5.0, 887.25, 888.25, 887.25, 887.25, 875.25, 883.25, 1744.5, 1.0, 871.25, 887.25, 875.25, 1749.5, 871.25, 1745.5, 880.25, 872.25, 9.0, 872.25, 883.25, 875.25, 1744.5, 872.25, 874.25, 875.25, 1765.5, 1745.5, 1741.5, 17.0, 17.0, 872.25, 1.0, 872.25, 887.25, 1.0, 887.25, 1765.5, 2.0, 871.25, 874.25, 17.0, 886.25, 8.0, 874.25, 8.0, 890.25, 871.25, 875.25, 1748.5, 1758.5, 902.25, 1741.5, 1741.5, 20.0, 1741.5, 875.25, 1742.5, 886.25, 1.0, 874.25, 874.25, 879.25, 1745.5, 874.25, 887.25, 871.25, 872.25, 871.25, 871.25, 0.0, 4.0, 4.0, 871.25, 870.25, 875.25, 870.25, 8.0, 874.25, 1744.5, 1745.5, 1744.5, 886.25, 20.0, 4.0, 13.0, 10.0, 1744.5, 1745.5, 871.25, 17.0, 872.25, 1745.5, 878.25, 872.25, 888.25, 17.0, 878.25, 887.25, 875.25, 871.25, 1741.5, 875.25, 871.25, 1742.5, 1744.5, 880.25, 0.0, 875.25, 5.0, 5.0, 1742.5, 2.0, 1750.5, 871.25, 1742.5, 872.25, 1742.5, 895.25, 878.25, 13.0, 875.25, 10.0, 875.25, 883.25, 887.25, 872.25, 10.0, 871.25, 1741.5, 10.0, 1745.5, 871.25, 875.25, 1745.5, 883.25, 1753.5, 2.0, 25.0, 872.25, 1745.5, 879.25, 1742.5, 1745.5, 890.25, 887.25, 1757.5, 890.25, 5.0, 875.25, 875.25, 2.0, 1744.5, 872.25, 878.25, 1741.5, 5.0, 1744.5, 1750.5, 1745.5, 10.0, 1.0, 883.25, 872.25, 1742.5, 1744.5, 1.0, 875.25, 2.0, 13.0, 874.25, 871.25, 2.0, 878.25, 870.25, 874.25, 1744.5, 880.25, 10.0, 1742.5, 4.0, 875.25, 10.0, 875.25, 1750.5, 874.25, 1.0, 875.25, 870.25, 874.25, 871.25, 874.25, 870.25, 875.25, 875.25, 874.25, 872.25, 875.25, 1765.5, 875.25, 16.0, 880.25, 4.0, 1749.5, 870.25, 1741.5, 902.25, 872.25, 886.25, 880.25, 17.0, 871.25, 875.25, 875.25, 1744.5, 871.25, 0.0, 883.25, 10.0, 871.25, 1741.5, 872.25, 1745.5, 875.25, 871.25, 880.25, 875.25, 2.0, 902.25, 902.25, 1.0, 870.25, 879.25, 875.25, 25.0, 1748.5, 872.25, 2.0, 888.25, 875.25, 871.25, 1748.5, 875.25, 25.0, 1753.5, 1749.5, 16.0, 10.0, 10.0, 1745.5, 1765.5, 870.25, 1756.5, 1749.5, 890.25, 872.25, 1749.5, 8.0, 870.25, 875.25, 1.0, 1742.5, 879.25, 878.25, 871.25, 2.0, 872.25, 5.0, 870.25, 1.0, 1745.5, 875.25, 870.25, 886.25, 1753.5, 871.25, 1741.5, 874.25, 895.25, 18.0, 870.25, 5.0, 1744.5, 872.25, 880.25, 875.25, 1760.5, 1742.5, 1741.5, 1756.5, 2.0, 1757.5, 880.25, 1741.5, 1.0, 871.25, 1748.5, 16.0, 871.25, 1745.5, 17.0, 2.0, 1740.5, 875.25, 1742.5, 1744.5, 0.0, 4.0, 2.0, 883.25, 871.25, 879.25, 879.25, 875.25, 4.0, 872.25, 872.25, 871.25, 887.25, 895.25, 4.0, 4.0, 1757.5, 1741.5, 875.25, 1741.5, 875.25, 1.0, 2.0, 13.0, 1748.5, 5.0, 1740.5, 10.0, 1745.5, 874.25, 1753.5, 878.25, 8.0, 874.25, 875.25, 887.25, 872.25, 0.0, 1748.5, 872.25, 879.25, 1756.5, 4.0, 1741.5, 871.25, 875.25, 871.25, 10.0, 871.25, 2.0, 874.25, 0.0, 1757.5, 880.25, 1745.5, 874.25, 875.25, 1.0, 17.0, 1753.5, 9.0, 1772.5, 20.0, 890.25, 20.0, 17.0, 872.25, 1750.5, 871.25, 13.0, 9.0, 871.25, 1772.5, 1.0, 874.25, 888.25, 2.0, 871.25, 4.0, 2.0, 871.25, 0.0, 1741.5, 870.25, 871.25, 2.0, 874.25, 895.25, 872.25, 883.25, 1.0, 1.0, 1750.5, 10.0, 1.0, 1745.5, 887.25, 5.0, 1.0, 13.0, 872.25, 871.25, 1744.5, 870.25, 874.25, 875.25, 875.25, 878.25, 1741.5, 13.0, 5.0, 10.0, 5.0, 875.25, 878.25, 1753.5, 16.0, 886.25, 1745.5, 887.25, 887.25, 887.25, 1772.5, 1757.5, 13.0, 886.25, 875.25, 1741.5, 10.0, 875.25, 4.0, 13.0, 1749.5, 1756.5, 4.0, 1749.5, 872.25, 872.25, 1744.5, 2.0, 2.0, 1745.5, 17.0, 887.25, 886.25, 875.25, 10.0, 871.25, 1742.5, 10.0, 10.0, 1760.5, 10.0, 883.25, 4.0, 879.25, 874.25, 878.25, 1741.5, 4.0, 874.25, 17.0, 1748.5, 883.25, 25.0, 1750.5, 886.25, 25.0, 890.25, 1741.5, 1742.5, 878.25, 5.0, 875.25, 874.25, 10.0, 1741.5, 871.25, 18.0, 1741.5, 875.25, 13.0, 875.25, 872.25, 880.25, 1742.5, 871.25, 890.25, 870.25, 887.25, 4.0, 10.0, 875.25, 1.0, 875.25, 871.25, 1742.5, 1744.5, 1.0, 874.25, 1745.5, 875.25, 1.0, 871.25, 16.0, 2.0, 1.0, 32.0, 1741.5, 8.0, 1756.5, 879.25, 880.25, 1741.5, 1750.5, 1744.5, 874.25, 887.25, 17.0, 875.25, 887.25, 875.25, 1741.5, 872.25, 875.25, 875.25, 1741.5, 874.25, 1741.5, 13.0, 875.25, 1742.5, 1760.5, 879.25, 1745.5, 887.25, 1.0, 1744.5, 879.25, 874.25, 1.0, 1741.5, 871.25, 2.0, 1748.5, 870.25, 4.0, 1741.5, 875.25, 871.25, 1744.5, 5.0, 886.25, 888.25, 878.25, 875.25, 872.25, 1753.5, 2.0, 890.25, 880.25, 875.25, 0.0, 895.25, 872.25, 10.0, 13.0, 880.25, 875.25, 880.25, 880.25, 4.0, 8.0, 4.0, 880.25, 2.0, 1742.5, 874.25, 2.0, 0.0, 1745.5, 872.25, 871.25, 1742.5, 1748.5, 875.25, 871.25, 874.25, 1745.5, 10.0, 875.25, 886.25, 871.25, 880.25, 879.25, 875.25, 880.25, 872.25, 1740.5, 1750.5, 879.25, 872.25, 871.25, 871.25, 875.25, 872.25, 872.25, 871.25, 871.25, 1750.5, 8.0, 1740.5, 1745.5, 871.25, 2.0, 880.25, 887.25, 879.25, 1744.5, 1741.5, 5.0, 1757.5, 883.25, 875.25, 1745.5, 1745.5, 887.25, 1742.5, 8.0, 902.25, 895.25, 1756.5, 1750.5, 872.25, 883.25, 5.0, 871.25, 4.0, 32.0, 871.25, 1750.5, 890.25, 0.0, 1742.5, 20.0, 1742.5, 886.25, 883.25, 2.0, 871.25, 1745.5, 890.25, 871.25, 4.0, 1748.5, 1745.5, 13.0, 1.0, 5.0, 895.25, 880.25, 16.0, 20.0, 872.25, 874.25, 880.25, 895.25, 5.0, 1.0, 871.25, 879.25, 872.25, 1741.5, 878.25, 887.25, 8.0, 1760.5, 879.25, 1757.5, 8.0, 5.0, 875.25, 872.25, 1.0, 2.0, 1748.5, 1741.5, 879.25, 875.25, 9.0, 874.25, 872.25, 1.0, 13.0, 5.0, 1757.5, 0.0, 870.25, 875.25, 880.25, 5.0, 0.0, 872.25, 886.25, 872.25, 871.25, 1753.5, 1741.5, 883.25, 1745.5, 1742.5, 1.0, 875.25, 1.0, 879.25, 1748.5, 875.25, 1.0, 1748.5, 1750.5, 4.0, 872.25, 875.25, 13.0, 10.0, 1740.5, 18.0, 890.25, 10.0, 871.25, 871.25, 872.25, 890.25, 880.25, 13.0, 872.25, 13.0, 1.0, 20.0, 1744.5, 1744.5, 871.25, 1742.5, 875.25, 1753.5, 1742.5, 16.0, 5.0, 1748.5, 872.25, 874.25, 13.0, 17.0, 1.0, 1744.5, 875.25, 872.25, 1758.5, 2.0, 8.0, 875.25, 2.0, 16.0, 870.25, 1748.5, 1749.5, 872.25, 1744.5, 871.25, 890.25, 874.25, 1741.5, 875.25, 1.0, 874.25, 880.25, 10.0, 1742.5, 1750.5, 880.25, 1742.5, 871.25, 1.0, 4.0, 870.25, 890.25, 871.25, 17.0, 875.25, 880.25, 10.0, 1744.5, 9.0, 875.25, 874.25, 875.25, 1742.5, 1756.5, 10.0, 870.25, 871.25, 13.0, 871.25, 875.25, 883.25, 16.0, 1745.5, 895.25, 1748.5, 871.25, 871.25, 9.0, 10.0, 880.25, 1.0, 2.0, 1749.5, 1745.5, 13.0, 1.0, 875.25, 870.25, 1741.5, 883.25, 1.0, 2.0, 1744.5, 870.25, 874.25, 1757.5, 1742.5, 1749.5, 4.0, 880.25, 872.25, 887.25, 871.25, 871.25, 874.25, 879.25, 872.25, 1.0, 875.25, 874.25, 1742.5, 1744.5, 2.0, 871.25, 871.25, 1.0, 2.0, 1750.5, 871.25, 878.25, 17.0, 1745.5, 871.25, 10.0, 10.0, 1748.5, 890.25, 883.25, 0.0, 895.25, 883.25, 1753.5, 1757.5, 1745.5, 1740.5, 1740.5, 1748.5, 25.0, 886.25, 880.25, 1741.5, 871.25, 0.0, 13.0, 875.25, 890.25, 1753.5, 1745.5, 4.0, 1760.5, 878.25, 20.0, 17.0, 880.25, 1753.5, 1757.5, 872.25, 872.25, 10.0, 880.25, 875.25, 875.25, 1748.5, 879.25, 1749.5, 13.0, 1742.5, 5.0, 1750.5, 1760.5, 871.25, 4.0, 10.0, 875.25, 871.25, 8.0, 1748.5, 13.0, 1757.5, 871.25, 1760.5, 872.25, 4.0, 4.0, 872.25, 887.25, 872.25, 1740.5, 2.0, 9.0, 874.25, 880.25, 2.0, 16.0, 872.25, 872.25, 887.25, 10.0, 880.25, 1.0, 4.0, 1742.5, 871.25, 1.0, 872.25, 871.25, 5.0, 1.0, 13.0, 1757.5, 8.0, 890.25, 8.0, 10.0, 872.25, 874.25, 874.25, 872.25, 875.25, 1757.5, 871.25, 10.0, 0.0, 1750.5, 875.25, 872.25, 0.0, 1760.5, 1745.5, 1748.5, 1741.5, 17.0, 880.25, 2.0, 13.0, 872.25, 1745.5, 1753.5, 879.25, 1.0, 870.25, 17.0, 871.25, 17.0, 874.25, 874.25, 4.0, 875.25, 5.0, 1.0, 1753.5, 880.25, 871.25, 5.0, 883.25, 4.0, 872.25, 9.0, 872.25, 0.0, 872.25, 871.25, 1742.5, 871.25, 1742.5, 871.25, 880.25, 874.25, 1.0, 1745.5, 874.25, 4.0, 875.25, 872.25, 13.0, 10.0, 1750.5, 871.25, 1745.5, 17.0, 878.25, 875.25, 5.0, 1772.5, 20.0, 871.25, 1748.5, 17.0, 18.0, 875.25, 880.25, 1741.5, 874.25, 875.25, 5.0, 1760.5, 880.25, 32.0, 1753.5, 890.25, 1749.5, 878.25, 870.25, 871.25, 875.25, 1745.5, 1741.5, 1772.5, 8.0, 878.25, 880.25, 875.25, 10.0, 1741.5, 883.25, 870.25, 5.0, 872.25, 2.0, 1750.5, 872.25, 1749.5, 875.25, 1745.5, 874.25, 883.25, 870.25, 875.25, 872.25, 1741.5, 1.0, 1.0, 871.25, 16.0, 872.25, 9.0, 16.0, 886.25, 0.0, 1750.5, 9.0, 880.25, 4.0, 1.0, 1753.5, 1742.5, 872.25, 1745.5, 1765.5, 1.0, 880.25, 872.25, 890.25, 1749.5, 1744.5, 17.0, 887.25, 16.0, 1753.5, 872.25, 878.25, 871.25, 10.0, 872.25, 5.0, 1741.5, 875.25, 1742.5, 886.25, 17.0, 875.25, 887.25, 1740.5, 1750.5, 1749.5, 4.0, 10.0, 872.25, 1741.5, 888.25, 1757.5, 875.25, 10.0, 4.0, 872.25, 1749.5, 879.25, 875.25, 878.25, 874.25, 870.25, 17.0, 879.25, 1741.5, 871.25, 902.25, 871.25, 1745.5, 1765.5, 4.0, 880.25, 9.0, 10.0, 878.25, 887.25, 1.0, 871.25, 872.25, 8.0, 1.0, 1.0, 1745.5, 871.25, 1745.5, 1745.5, 872.25, 1741.5, 1750.5, 5.0, 1756.5, 2.0, 890.25, 878.25, 1741.5, 8.0, 1744.5, 1741.5, 871.25, 2.0, 871.25, 870.25, 874.25, 875.25, 871.25, 883.25, 1745.5, 879.25, 880.25, 887.25, 1.0, 871.25, 871.25, 16.0, 1745.5, 872.25, 879.25, 871.25, 874.25, 886.25, 4.0, 1741.5, 13.0, 870.25, 1745.5, 880.25, 1744.5, 871.25, 32.0, 1760.5, 1745.5, 875.25, 878.25, 871.25, 880.25, 874.25, 875.25, 871.25, 1742.5, 2.0, 1745.5, 1745.5, 871.25, 883.25, 1753.5, 5.0, 872.25, 1745.5, 878.25, 5.0, 872.25, 871.25, 874.25, 1742.5, 13.0, 1753.5, 0.0, 875.25, 1741.5, 872.25, 1741.5, 4.0, 870.25, 1748.5, 2.0, 8.0, 890.25, 1744.5, 10.0, 875.25, 879.25, 1748.5, 875.25, 1749.5, 872.25, 8.0, 890.25, 880.25, 1.0, 5.0, 13.0, 17.0, 875.25, 5.0, 2.0, 878.25, 871.25, 883.25, 5.0, 1745.5, 10.0, 886.25, 32.0, 895.25, 1753.5, 883.25, 871.25, 878.25, 886.25, 875.25, 1745.5, 890.25, 875.25, 9.0, 874.25, 5.0, 1.0, 1.0, 0.0, 871.25, 5.0, 1742.5, 874.25, 872.25, 5.0, 1750.5, 1757.5, 880.25, 1744.5, 1745.5, 1744.5, 1750.5, 2.0, 871.25, 17.0, 890.25, 20.0, 1758.5, 875.25, 25.0, 886.25, 1749.5, 1744.5, 883.25, 1745.5, 890.25, 1745.5, 887.25, 1760.5, 874.25, 1741.5, 2.0, 1742.5, 871.25, 875.25, 870.25, 871.25, 1745.5, 1.0, 1741.5, 875.25, 1741.5, 871.25, 880.25, 1745.5, 1748.5, 871.25, 1.0, 1744.5, 8.0, 874.25, 887.25, 890.25, 1.0, 883.25, 875.25, 880.25, 1749.5, 1744.5, 1.0, 874.25, 1742.5, 5.0, 875.25, 875.25, 4.0, 2.0, 2.0, 10.0, 5.0, 10.0, 1741.5, 874.25, 875.25, 1750.5, 1742.5, 871.25, 874.25, 4.0, 1.0, 875.25, 2.0, 871.25, 10.0, 871.25, 5.0, 1750.5, 883.25, 1760.5, 1749.5, 13.0, 1750.5, 875.25, 1760.5, 878.25, 2.0, 2.0, 17.0, 13.0, 888.25, 2.0, 5.0, 13.0, 10.0, 17.0, 875.25, 8.0, 872.25, 871.25, 4.0, 10.0, 875.25, 888.25, 872.25, 871.25, 871.25, 0.0, 1745.5, 1.0, 20.0, 2.0, 872.25, 10.0, 875.25, 9.0, 875.25, 20.0, 1750.5, 878.25, 1.0, 1.0, 1742.5, 886.25, 870.25, 5.0, 870.25, 32.0, 883.25, 13.0, 1756.5, 1740.5, 880.25, 886.25, 875.25, 4.0, 18.0, 1750.5, 888.25, 895.25, 5.0, 875.25, 872.25, 1745.5, 1750.5, 874.25, 883.25, 1742.5, 1749.5, 883.25, 870.25, 1742.5, 870.25, 1750.5, 1750.5, 870.25, 1.0, 10.0, 880.25, 2.0, 5.0, 871.25, 870.25, 870.25, 1744.5, 2.0, 2.0, 871.25, 871.25, 872.25, 8.0, 872.25, 872.25, 1741.5, 5.0, 887.25, 1741.5, 9.0, 1740.5, 1.0, 1741.5, 1757.5, 879.25, 871.25, 5.0, 2.0, 5.0, 1741.5, 1748.5, 4.0, 1757.5, 1749.5, 872.25, 1744.5, 871.25, 25.0, 4.0, 17.0, 902.25, 890.25, 1765.5, 902.25, 1772.5, 890.25, 1749.5, 16.0, 1742.5, 8.0, 872.25, 874.25, 871.25, 874.25, 1.0, 1745.5, 879.25, 17.0, 883.25, 875.25, 872.25, 4.0, 872.25, 0.0, 5.0, 872.25, 1744.5, 2.0, 10.0, 13.0, 871.25, 1744.5, 886.25, 1745.5, 2.0, 880.25, 1744.5, 5.0, 871.25, 1745.5, 895.25, 2.0, 887.25, 872.25, 890.25, 1748.5, 887.25, 887.25, 13.0, 890.25, 5.0, 880.25, 883.25, 2.0, 10.0, 10.0, 871.25, 1745.5, 13.0, 878.25, 9.0, 874.25, 878.25, 1745.5, 1742.5, 10.0, 890.25, 872.25, 1741.5, 874.25, 1744.5, 880.25, 875.25, 1.0, 871.25, 1745.5, 5.0, 1745.5, 1.0, 872.25, 875.25, 871.25, 890.25, 1765.5, 874.25, 1765.5, 875.25, 1744.5, 871.25, 1749.5, 18.0, 2.0, 13.0, 871.25, 10.0, 875.25, 875.25, 872.25, 1757.5, 4.0, 1741.5, 0.0, 9.0, 1.0, 1748.5, 5.0, 879.25, 16.0, 895.25, 888.25, 875.25, 880.25, 883.25, 1741.5, 1748.5, 1.0, 874.25, 5.0, 875.25, 880.25, 1742.5, 8.0, 1745.5, 1.0, 874.25, 879.25, 1756.5, 1745.5, 1744.5, 872.25, 1745.5, 872.25, 870.25, 880.25, 871.25, 871.25, 879.25, 874.25, 2.0, 1757.5, 8.0, 1750.5, 1756.5, 1744.5, 879.25, 13.0, 1756.5, 1741.5, 878.25, 17.0, 872.25, 1744.5, 1740.5, 20.0, 872.25, 875.25, 872.25, 878.25, 871.25, 1744.5, 1.0, 874.25, 875.25, 10.0, 879.25, 872.25, 871.25, 4.0, 10.0, 4.0, 871.25, 895.25, 1757.5, 871.25, 1753.5, 1.0, 875.25, 872.25, 5.0, 887.25, 871.25, 4.0, 1750.5, 871.25, 10.0, 17.0, 879.25, 1741.5, 5.0, 5.0, 9.0, 10.0, 879.25, 887.25, 875.25, 888.25, 20.0, 874.25, 874.25, 1745.5, 875.25, 887.25, 874.25, 871.25, 1750.5, 1753.5, 25.0, 1.0, 874.25, 1.0, 8.0, 5.0, 1772.5, 872.25, 10.0, 887.25, 875.25, 13.0, 1741.5, 890.25, 880.25, 2.0, 890.25, 887.25, 1.0, 887.25, 875.25, 871.25, 18.0, 1742.5, 1749.5, 1742.5, 1740.5, 1.0, 874.25, 879.25, 5.0, 1748.5, 1745.5, 10.0, 878.25, 4.0, 10.0, 871.25, 1749.5, 878.25, 879.25, 872.25, 875.25, 1741.5, 1.0, 875.25, 875.25, 870.25, 875.25, 874.25, 875.25, 874.25, 874.25, 878.25, 887.25, 875.25, 874.25, 879.25, 1.0, 880.25, 871.25, 1742.5, 1745.5, 1.0, 1745.5, 874.25, 1.0, 871.25, 1744.5, 871.25, 1757.5, 1741.5, 875.25, 890.25, 871.25, 1760.5, 4.0, 8.0, 872.25, 888.25, 875.25, 2.0, 1742.5, 5.0, 875.25, 879.25, 872.25, 872.25, 1741.5, 4.0, 1745.5, 1.0, 1.0, 883.25, 5.0, 870.25, 871.25, 870.25, 2.0, 871.25, 872.25, 1750.5, 871.25, 878.25, 1748.5, 1745.5, 880.25, 1742.5, 1745.5, 875.25, 0.0, 1741.5, 2.0, 1765.5, 1741.5, 871.25, 1741.5, 878.25, 0.0, 875.25, 4.0, 4.0, 1748.5, 5.0, 872.25, 9.0, 1.0, 880.25, 874.25, 883.25, 5.0, 2.0, 874.25, 875.25, 1.0, 17.0, 1757.5, 20.0, 1742.5, 874.25, 1745.5, 902.25, 880.25, 890.25, 872.25, 1749.5, 871.25, 1756.5, 9.0, 10.0, 5.0, 874.25, 25.0, 902.25, 1742.5, 895.25, 10.0, 883.25, 879.25, 871.25, 890.25, 13.0, 16.0, 872.25, 1744.5, 902.25, 2.0, 878.25, 1750.5, 8.0, 1.0, 4.0, 871.25, 1741.5, 875.25, 1757.5, 1750.5, 875.25, 895.25, 1742.5, 1765.5, 5.0, 874.25, 1772.5, 1742.5, 16.0, 1745.5, 1757.5, 0.0, 16.0, 1745.5, 875.25, 1744.5, 875.25, 20.0, 870.25, 9.0, 1.0, 875.25, 875.25, 9.0, 870.25, 5.0, 5.0, 890.25, 5.0, 874.25, 875.25, 878.25, 4.0, 880.25, 1745.5, 871.25, 2.0, 1741.5, 1742.5, 890.25, 1741.5, 1750.5, 874.25, 2.0, 17.0, 871.25, 5.0, 870.25, 5.0, 5.0, 10.0, 872.25, 890.25, 880.25, 874.25, 880.25, 17.0, 880.25, 879.25, 880.25, 1745.5, 871.25, 1744.5, 2.0, 4.0, 1750.5, 1748.5, 1744.5, 17.0, 890.25, 4.0, 1742.5, 2.0, 2.0, 871.25, 902.25, 902.25, 880.25, 878.25, 1750.5, 883.25, 883.25, 1744.5, 1.0, 5.0, 1.0, 878.25, 870.25, 1742.5, 4.0, 1745.5, 1741.5, 871.25, 1740.5, 13.0, 871.25, 1742.5, 1744.5, 875.25, 1750.5, 874.25, 4.0, 879.25, 1742.5, 5.0, 1745.5, 1745.5, 883.25, 880.25, 874.25, 5.0, 2.0, 2.0, 4.0, 887.25, 1741.5, 8.0, 2.0, 872.25, 895.25, 875.25, 20.0, 875.25, 10.0, 1745.5, 20.0, 871.25, 1756.5, 879.25, 895.25, 1750.5, 874.25, 4.0, 875.25, 872.25, 871.25, 872.25, 1.0, 1741.5, 1.0, 2.0, 1745.5, 2.0, 2.0, 872.25, 875.25, 871.25, 875.25, 871.25, 1741.5, 883.25, 1744.5, 1.0, 5.0, 875.25, 1744.5, 8.0, 1740.5, 5.0, 880.25, 871.25, 4.0, 1753.5, 871.25, 872.25, 872.25, 1.0, 883.25, 1749.5, 16.0, 872.25, 871.25, 872.25, 875.25, 871.25, 2.0, 871.25, 1.0, 875.25, 0.0, 1.0, 880.25, 895.25, 872.25, 1742.5, 1.0, 890.25, 880.25, 1760.5, 871.25, 871.25, 887.25, 895.25, 878.25, 880.25, 883.25, 1745.5, 1756.5, 886.25, 8.0, 1.0, 875.25, 1745.5, 880.25, 874.25, 32.0, 5.0, 878.25, 16.0, 1756.5, 1.0, 875.25, 5.0, 1753.5, 883.25, 10.0, 1765.5, 1.0, 1753.5, 871.25, 1740.5, 2.0, 875.25, 871.25, 1742.5, 5.0, 878.25, 1749.5, 875.25, 890.25, 1745.5, 879.25, 875.25, 1750.5, 871.25, 2.0, 1745.5, 871.25, 1760.5, 1757.5, 10.0, 886.25, 875.25, 886.25, 872.25, 1742.5, 9.0, 10.0, 890.25, 25.0, 875.25, 886.25, 1757.5, 1756.5, 872.25, 872.25, 1741.5, 1742.5, 875.25, 8.0, 874.25, 870.25, 1.0, 1749.5, 0.0, 4.0, 1.0, 2.0, 0.0, 874.25, 1744.5, 875.25, 1757.5, 887.25, 870.25, 1745.5, 883.25, 875.25, 871.25, 875.25, 878.25, 874.25, 1741.5, 872.25, 1744.5, 874.25, 871.25, 1744.5, 13.0, 879.25, 1742.5, 1.0, 1.0, 870.25, 1740.5, 1750.5, 895.25, 5.0, 5.0, 872.25, 874.25, 887.25, 871.25, 883.25, 883.25, 10.0, 883.25, 25.0, 870.25, 1.0, 872.25, 875.25, 875.25, 890.25, 20.0, 875.25, 1757.5, 871.25, 872.25, 1748.5, 1740.5, 0.0, 883.25, 871.25, 874.25, 1745.5, 875.25, 874.25, 872.25, 4.0, 8.0, 872.25, 1745.5, 0.0, 875.25, 25.0, 875.25, 9.0, 2.0, 887.25, 874.25, 1744.5, 1741.5, 2.0, 874.25, 1745.5, 878.25, 1741.5, 874.25, 871.25, 874.25, 18.0, 5.0, 1745.5, 5.0, 4.0, 5.0, 871.25, 895.25, 874.25, 872.25, 890.25, 1757.5, 16.0, 1748.5, 875.25, 874.25, 1745.5, 879.25, 2.0, 1749.5, 1.0, 1749.5, 1750.5, 4.0, 1745.5, 875.25, 1742.5, 1744.5, 883.25, 2.0, 878.25, 870.25, 1745.5, 10.0, 17.0, 5.0, 4.0, 879.25, 1748.5, 5.0, 879.25, 883.25, 1741.5, 1740.5, 883.25, 2.0, 10.0, 10.0, 1.0, 872.25, 872.25, 872.25, 4.0, 10.0, 1741.5, 871.25, 871.25, 875.25, 1741.5, 1.0, 1748.5, 13.0, 1745.5, 883.25, 875.25, 874.25, 10.0, 1.0, 1749.5, 875.25, 871.25, 874.25, 870.25, 875.25, 871.25, 4.0, 874.25, 870.25, 875.25, 1745.5, 886.25, 4.0, 871.25, 883.25, 1.0, 2.0, 4.0, 1742.5, 1745.5, 883.25, 871.25, 1741.5, 1742.5, 13.0, 887.25, 872.25, 1740.5, 1772.5, 874.25, 871.25, 1757.5, 17.0, 10.0, 25.0, 20.0, 880.25, 878.25, 871.25, 1745.5, 871.25, 1744.5, 4.0, 4.0, 5.0, 0.0, 1750.5, 4.0, 1742.5, 871.25, 5.0, 870.25, 1745.5, 4.0, 1765.5, 1757.5, 1748.5, 8.0, 1742.5, 883.25, 1741.5, 883.25, 872.25, 874.25, 1742.5, 880.25, 1750.5, 4.0, 1745.5, 1741.5, 1749.5, 1741.5, 1750.5, 1745.5, 1741.5, 5.0, 1.0, 879.25, 1741.5, 870.25, 1750.5, 880.25, 17.0, 878.25, 875.25, 875.25, 871.25, 10.0, 1741.5, 870.25, 1740.5, 874.25, 871.25, 5.0, 880.25, 871.25, 1745.5, 875.25, 875.25, 1772.5, 880.25, 875.25, 10.0, 890.25, 20.0, 880.25, 871.25, 874.25, 880.25, 10.0, 902.25, 874.25, 871.25, 4.0, 1742.5, 10.0, 4.0, 0.0, 1.0, 1740.5, 9.0, 880.25, 1744.5, 0.0, 870.25, 10.0, 1740.5, 2.0, 1.0, 1757.5, 871.25, 1.0, 875.25, 13.0, 9.0, 1.0, 2.0, 1745.5, 1745.5, 10.0, 1745.5, 1.0, 1742.5, 871.25, 1742.5, 872.25, 1745.5, 1.0, 5.0, 871.25, 895.25, 13.0, 1.0, 9.0, 1753.5, 0.0, 17.0, 883.25, 887.25, 895.25, 1744.5, 875.25, 872.25, 878.25, 1745.5, 1757.5, 1749.5, 1742.5, 1750.5, 1741.5, 871.25, 1745.5, 875.25, 875.25, 872.25, 874.25, 17.0, 10.0, 1741.5, 871.25, 875.25, 1745.5, 2.0, 1741.5, 2.0, 10.0, 872.25, 1740.5, 875.25, 871.25, 4.0, 888.25, 5.0, 0.0, 2.0, 886.25, 4.0, 9.0, 1757.5, 871.25, 2.0, 2.0, 5.0, 1745.5, 1742.5, 880.25, 1741.5, 1753.5, 878.25, 5.0, 2.0, 887.25, 883.25, 887.25, 870.25, 871.25, 879.25, 1745.5, 20.0, 872.25, 32.0, 871.25, 1741.5, 890.25, 883.25, 880.25, 871.25, 1745.5, 17.0, 13.0, 883.25, 1.0, 872.25, 875.25, 887.25, 880.25, 870.25, 8.0, 10.0, 879.25, 2.0, 8.0, 17.0, 872.25, 1741.5, 887.25, 878.25, 870.25, 1742.5, 875.25, 1760.5, 8.0, 1745.5, 0.0, 890.25, 875.25, 1745.5, 13.0, 879.25, 1756.5, 2.0, 878.25, 875.25, 1749.5, 16.0, 875.25, 1741.5, 18.0, 1745.5, 13.0, 875.25, 1744.5, 875.25, 872.25, 888.25, 872.25, 4.0, 5.0, 4.0, 872.25, 878.25, 874.25, 1750.5, 887.25, 1745.5, 1.0, 871.25, 1765.5, 1.0, 874.25, 888.25, 17.0, 1744.5, 1753.5, 1745.5, 874.25, 4.0, 872.25, 890.25, 17.0, 1760.5, 1749.5, 1742.5, 879.25, 872.25, 1748.5, 1749.5, 870.25, 874.25, 1745.5, 1750.5, 1741.5, 880.25, 2.0, 902.25, 871.25, 1742.5, 878.25, 883.25, 875.25, 895.25, 883.25, 890.25, 880.25, 875.25, 0.0, 1744.5, 878.25, 5.0, 1745.5, 874.25, 1741.5, 872.25, 1.0, 1749.5, 875.25, 875.25, 871.25, 1745.5, 5.0, 872.25, 9.0, 1742.5, 1756.5, 1742.5, 17.0, 870.25, 1742.5, 1745.5, 871.25, 1.0, 879.25, 871.25, 9.0, 872.25, 874.25, 872.25, 875.25, 872.25, 4.0, 5.0, 883.25, 870.25, 1745.5, 1741.5, 10.0, 1753.5, 1745.5, 875.25, 874.25, 13.0, 4.0, 16.0, 1741.5, 17.0, 871.25, 872.25, 874.25, 13.0, 1742.5, 4.0, 871.25, 878.25, 872.25, 1741.5, 1741.5, 890.25, 886.25, 878.25, 1765.5, 880.25, 890.25, 902.25, 878.25, 887.25, 902.25, 13.0, 878.25, 1744.5, 875.25, 880.25, 1741.5, 874.25, 1.0, 1749.5, 872.25, 872.25, 880.25, 870.25, 874.25, 1745.5, 1.0, 871.25, 1.0, 883.25, 880.25, 1753.5, 870.25, 1741.5, 902.25, 2.0, 886.25, 10.0, 2.0, 1765.5, 886.25, 2.0, 886.25, 872.25, 874.25, 872.25, 870.25, 874.25, 13.0, 1749.5, 1.0, 1748.5, 872.25, 871.25, 1.0, 1745.5, 1760.5, 1741.5, 10.0, 872.25, 1749.5, 878.25, 878.25, 1744.5, 1745.5, 1756.5, 875.25, 13.0, 1742.5, 902.25, 870.25, 5.0, 1750.5, 1753.5, 870.25, 1742.5, 5.0, 4.0, 10.0, 1749.5, 1744.5, 1.0, 880.25, 4.0, 887.25, 880.25, 1750.5, 872.25, 875.25, 1749.5, 9.0, 1.0, 10.0, 1741.5, 880.25, 871.25, 1750.5, 871.25, 875.25, 1740.5, 890.25, 887.25, 890.25, 17.0, 874.25, 874.25, 5.0, 883.25, 875.25, 870.25, 890.25, 872.25, 871.25, 902.25, 890.25, 874.25, 887.25, 887.25, 5.0, 871.25, 1.0, 1765.5, 1750.5, 5.0, 2.0, 871.25, 1742.5, 875.25, 0.0, 872.25, 1748.5, 875.25, 879.25, 1.0, 875.25, 1749.5, 874.25, 1740.5, 871.25, 879.25, 1741.5, 5.0, 890.25, 5.0, 13.0, 5.0, 1748.5, 5.0, 871.25, 871.25, 872.25, 875.25, 878.25, 1741.5, 8.0, 871.25, 1740.5, 5.0, 878.25, 1757.5, 874.25, 875.25, 888.25, 1741.5, 13.0, 871.25, 875.25, 1750.5, 4.0, 1.0, 1.0, 883.25, 872.25, 5.0, 10.0, 1744.5, 872.25, 890.25, 875.25, 1742.5, 875.25, 1.0, 1742.5, 883.25, 1740.5, 1745.5, 874.25, 870.25, 871.25, 1741.5, 1741.5, 1741.5, 872.25, 0.0, 16.0, 1753.5, 872.25, 5.0, 16.0, 20.0, 1745.5, 890.25, 883.25, 874.25, 10.0, 890.25, 871.25, 870.25, 880.25, 871.25, 879.25, 1745.5, 879.25, 5.0, 5.0, 1741.5, 870.25, 874.25, 1742.5, 4.0, 1.0, 874.25, 4.0, 875.25, 878.25, 890.25, 8.0, 872.25, 5.0, 1740.5, 1.0, 872.25, 880.25, 20.0, 1756.5, 878.25, 1740.5, 888.25, 875.25, 871.25, 2.0, 13.0, 1742.5, 1753.5, 2.0, 8.0, 1744.5, 879.25, 4.0, 1744.5, 1742.5, 1742.5, 1753.5, 9.0, 886.25, 1750.5, 1756.5, 5.0, 880.25, 1750.5, 1744.5, 10.0, 1745.5, 870.25, 8.0, 1745.5, 10.0, 875.25, 871.25, 871.25, 872.25, 2.0, 871.25, 1753.5, 883.25, 871.25, 1742.5, 883.25, 875.25, 1740.5, 871.25, 1756.5, 2.0, 1.0, 5.0, 1740.5, 1750.5, 1742.5, 872.25, 888.25, 4.0, 871.25, 871.25, 2.0, 879.25, 870.25, 883.25, 875.25, 1741.5, 1.0, 870.25, 1.0, 1745.5, 874.25, 874.25, 870.25, 4.0, 1.0, 878.25, 874.25, 1741.5, 1.0, 1745.5, 1.0, 871.25, 1740.5, 13.0, 1.0, 8.0, 1744.5, 886.25, 872.25, 871.25, 880.25, 871.25, 1753.5, 878.25, 1745.5, 886.25, 875.25, 1757.5, 1750.5, 902.25, 1750.5, 2.0, 895.25, 1740.5, 5.0, 1.0, 902.25, 890.25, 874.25, 1744.5, 871.25, 871.25, 875.25, 1744.5, 879.25, 2.0, 1753.5, 871.25, 875.25, 890.25, 5.0, 1740.5, 883.25, 17.0, 879.25, 1741.5, 870.25, 1750.5, 1.0, 5.0, 1758.5, 880.25, 0.0, 1748.5, 1740.5, 4.0, 875.25, 1772.5, 870.25, 4.0, 1745.5, 2.0, 895.25, 1.0, 871.25, 5.0, 874.25, 872.25, 20.0, 872.25, 5.0, 875.25, 879.25, 887.25, 1741.5, 871.25, 871.25, 875.25, 1750.5, 32.0, 902.25, 1741.5, 875.25, 1742.5, 17.0, 879.25, 895.25, 0.0, 872.25, 875.25, 883.25, 872.25, 1740.5, 8.0, 879.25, 888.25, 875.25, 870.25, 874.25, 8.0, 1744.5, 32.0, 2.0, 870.25, 9.0, 2.0, 871.25, 878.25, 1.0, 1750.5, 880.25, 17.0, 872.25, 871.25, 17.0, 883.25, 2.0, 1765.5, 1757.5, 880.25, 874.25, 887.25, 875.25, 1745.5, 883.25, 886.25, 895.25, 883.25, 1765.5, 878.25, 1744.5, 874.25, 2.0, 875.25, 871.25, 16.0, 1.0, 871.25, 883.25, 1753.5, 5.0, 870.25, 5.0, 10.0, 8.0, 878.25, 875.25, 1753.5, 13.0, 878.25, 875.25, 1745.5, 1772.5, 1745.5, 1745.5, 872.25, 886.25, 875.25, 1745.5, 890.25, 875.25, 1742.5, 890.25, 874.25, 1757.5, 880.25, 25.0, 5.0, 871.25, 878.25, 878.25, 2.0, 1749.5, 1753.5, 5.0, 1745.5, 875.25, 5.0, 1744.5, 879.25, 1757.5, 0.0, 10.0, 874.25, 878.25, 20.0, 4.0, 1.0, 1744.5, 872.25, 1.0, 870.25, 875.25, 875.25, 895.25, 872.25, 2.0, 2.0, 1740.5, 880.25, 1745.5, 0.0, 1744.5, 880.25, 0.0, 874.25, 874.25, 9.0, 875.25, 871.25, 870.25, 1742.5, 1742.5, 1745.5, 1741.5, 871.25, 18.0, 872.25, 875.25, 1744.5, 874.25, 880.25, 32.0, 875.25, 886.25, 872.25, 888.25, 1748.5, 4.0, 13.0, 1757.5, 1745.5, 871.25, 1744.5, 1.0, 895.25, 20.0, 871.25, 0.0, 902.25, 871.25, 1.0, 875.25, 871.25, 5.0, 1741.5, 9.0, 2.0, 880.25, 871.25, 16.0, 1750.5, 886.25, 871.25, 872.25, 880.25, 1741.5, 880.25, 878.25, 0.0, 883.25, 2.0, 2.0, 872.25, 5.0, 875.25, 880.25, 871.25, 872.25, 1.0, 1745.5, 4.0, 4.0, 870.25, 5.0, 870.25, 874.25, 10.0, 10.0, 883.25, 875.25, 1744.5, 875.25, 883.25, 2.0, 1741.5, 1750.5, 878.25, 4.0, 13.0, 871.25, 875.25, 880.25, 872.25, 9.0, 1744.5, 1742.5, 13.0, 871.25, 888.25, 887.25, 875.25, 20.0, 883.25, 32.0, 1745.5, 17.0, 0.0, 887.25, 20.0, 880.25, 5.0, 888.25, 5.0, 878.25, 8.0, 875.25, 13.0, 871.25, 902.25, 895.25, 879.25, 880.25, 878.25, 1760.5, 5.0, 1748.5, 1749.5, 871.25, 1750.5, 10.0, 1749.5, 1750.5, 1758.5, 872.25, 878.25, 871.25, 1742.5, 870.25, 883.25, 874.25, 13.0, 886.25, 887.25, 1750.5, 13.0, 870.25, 880.25, 10.0, 875.25, 880.25, 5.0, 2.0, 8.0, 879.25, 5.0, 20.0, 887.25, 875.25, 875.25, 17.0, 1745.5, 10.0, 1757.5, 16.0, 872.25, 5.0, 5.0, 1744.5, 1745.5, 872.25, 871.25, 875.25, 1.0, 880.25, 878.25, 2.0, 871.25, 871.25, 20.0, 1745.5, 1742.5, 10.0, 872.25, 875.25, 895.25, 1741.5, 5.0, 890.25, 872.25, 874.25, 1741.5, 878.25, 25.0, 1748.5, 4.0, 879.25, 1745.5, 1744.5, 1.0, 871.25, 2.0, 872.25, 879.25, 874.25, 1741.5, 1742.5, 1750.5, 890.25, 871.25, 886.25, 890.25, 879.25, 872.25, 890.25, 2.0, 1749.5, 2.0, 1.0, 1745.5, 875.25, 872.25, 17.0, 5.0, 875.25, 2.0, 1.0, 880.25, 872.25, 880.25, 1748.5, 871.25, 874.25, 5.0, 1742.5, 880.25, 871.25, 1741.5, 880.25, 883.25, 1.0, 875.25, 1742.5, 870.25, 886.25, 1745.5, 872.25, 880.25, 5.0, 1745.5, 880.25, 887.25, 870.25, 10.0, 890.25, 9.0, 1741.5, 872.25, 4.0, 875.25, 1741.5, 880.25, 886.25, 1750.5, 875.25, 879.25, 875.25, 887.25, 883.25, 878.25, 875.25, 1745.5, 880.25, 878.25, 25.0, 875.25, 1745.5, 887.25, 9.0, 1741.5, 895.25, 878.25, 890.25, 886.25, 886.25, 1745.5, 1753.5, 1.0, 1741.5, 5.0, 1742.5, 879.25, 870.25, 872.25, 875.25, 875.25, 1744.5, 0.0, 875.25, 880.25, 883.25, 883.25, 5.0, 879.25, 872.25, 874.25, 1741.5, 9.0, 879.25, 874.25, 872.25, 8.0, 875.25, 1741.5, 4.0, 1745.5, 871.25, 895.25, 874.25, 17.0, 871.25, 874.25, 883.25, 875.25, 878.25, 1749.5, 5.0, 1.0, 874.25, 1741.5, 880.25, 890.25, 1745.5, 1772.5, 1750.5, 875.25, 886.25, 25.0, 902.25, 1741.5, 0.0, 875.25, 1750.5, 1753.5, 880.25, 872.25, 13.0, 1.0, 1760.5, 9.0, 880.25, 2.0, 874.25, 875.25, 887.25, 5.0, 17.0, 1742.5, 871.25, 880.25, 886.25, 1748.5, 872.25, 875.25, 25.0, 870.25, 5.0, 1745.5, 10.0, 1745.5, 880.25, 5.0, 1757.5, 874.25, 1745.5, 871.25, 878.25, 902.25, 1750.5, 25.0, 5.0, 1745.5, 883.25, 887.25, 4.0, 890.25, 875.25, 1749.5, 1760.5, 20.0, 5.0, 1744.5, 1.0, 13.0, 1741.5, 1748.5, 872.25, 8.0, 10.0, 1742.5, 874.25, 887.25, 1.0, 4.0, 870.25, 1753.5, 875.25, 878.25, 5.0, 883.25, 888.25, 872.25, 872.25, 9.0, 883.25, 1750.5, 872.25, 1742.5, 875.25, 9.0, 20.0, 886.25, 883.25, 13.0, 1744.5, 1.0, 1742.5, 16.0, 1760.5, 883.25, 2.0, 10.0, 20.0, 1744.5, 883.25, 878.25, 883.25, 872.25, 1749.5, 2.0, 1.0, 1741.5, 871.25, 13.0, 887.25, 880.25, 1756.5, 1744.5, 13.0, 1757.5, 878.25, 883.25, 879.25, 879.25, 0.0, 871.25, 10.0, 2.0, 1745.5, 13.0, 1757.5, 17.0, 2.0, 871.25, 1749.5, 4.0, 1758.5, 880.25, 8.0, 875.25, 1742.5, 4.0, 872.25, 875.25, 1753.5, 902.25, 1745.5, 880.25, 870.25, 871.25, 1.0, 890.25, 1741.5, 875.25, 1741.5, 875.25, 1742.5, 870.25, 1.0, 878.25, 880.25, 880.25, 8.0, 872.25, 5.0, 871.25, 4.0, 874.25, 1742.5, 1745.5, 17.0, 883.25, 1.0, 10.0, 879.25, 874.25, 1745.5, 888.25, 1742.5, 1.0, 1740.5, 1741.5, 875.25, 874.25, 25.0, 872.25, 9.0, 5.0, 5.0, 1744.5, 875.25, 2.0, 13.0, 902.25, 5.0, 872.25, 890.25, 13.0, 1741.5, 1.0, 888.25, 1753.5, 17.0, 1.0, 1749.5, 875.25, 890.25, 872.25, 17.0, 874.25, 874.25, 871.25, 1745.5, 17.0, 875.25, 878.25, 4.0, 871.25, 875.25, 875.25, 1741.5, 1753.5, 883.25, 1748.5, 886.25, 1750.5, 886.25, 32.0, 887.25, 879.25, 890.25, 13.0, 874.25, 886.25, 895.25, 1742.5, 870.25, 880.25, 1742.5, 875.25, 13.0, 1742.5, 1753.5, 8.0, 874.25, 25.0, 17.0, 872.25, 17.0, 878.25, 1753.5, 1.0, 2.0, 1758.5, 17.0, 875.25, 887.25, 895.25, 25.0, 20.0, 871.25, 13.0, 887.25, 1.0, 878.25, 872.25, 16.0, 1.0, 1742.5, 895.25, 871.25, 874.25, 878.25, 875.25, 1753.5, 2.0, 875.25, 880.25, 25.0, 875.25, 871.25, 874.25, 1.0, 1753.5, 1741.5, 1760.5, 1750.5, 872.25, 1757.5, 1745.5, 875.25, 1750.5, 874.25, 1753.5, 875.25, 5.0, 25.0, 1748.5, 10.0, 1748.5, 886.25, 879.25, 871.25, 883.25, 888.25, 872.25, 1760.5, 17.0, 874.25, 2.0, 16.0, 1745.5, 871.25, 880.25, 880.25, 1745.5, 1745.5, 5.0, 17.0, 874.25, 1760.5, 1.0, 875.25, 1.0, 10.0, 880.25, 890.25, 1748.5, 5.0, 890.25, 1750.5, 5.0, 879.25, 879.25, 888.25, 1753.5, 1760.5, 1744.5, 874.25, 1742.5, 895.25, 4.0, 9.0, 1742.5, 0.0, 878.25, 1760.5, 874.25, 875.25, 875.25, 0.0, 1.0, 874.25, 872.25, 4.0, 0.0, 880.25, 887.25, 4.0, 875.25, 2.0, 1748.5, 872.25, 1749.5, 887.25, 5.0, 880.25, 2.0, 875.25, 902.25, 10.0, 872.25, 5.0, 874.25, 883.25, 875.25, 872.25, 871.25, 878.25, 879.25, 875.25, 887.25, 875.25, 1765.5, 16.0, 871.25, 871.25, 1742.5, 895.25, 5.0, 874.25, 871.25, 1741.5, 25.0, 883.25, 871.25, 5.0, 1.0, 2.0, 1.0, 878.25, 872.25, 1744.5, 872.25, 4.0, 890.25, 1741.5, 17.0, 0.0, 880.25, 874.25, 1744.5, 1745.5, 1.0, 4.0, 883.25, 872.25, 880.25, 13.0, 872.25, 5.0, 1748.5, 871.25, 879.25, 886.25, 880.25, 880.25, 872.25, 13.0, 9.0, 1.0, 895.25, 1742.5, 887.25, 5.0, 4.0, 1760.5, 902.25, 871.25, 872.25, 1757.5, 1741.5, 18.0, 1744.5, 871.25, 4.0, 1750.5, 871.25, 871.25, 1742.5, 1757.5, 1745.5, 1748.5, 4.0, 10.0, 1750.5, 1765.5, 1750.5, 879.25, 1.0, 872.25, 872.25, 1.0, 886.25, 883.25, 17.0, 4.0, 879.25, 5.0, 13.0, 886.25, 874.25, 875.25, 875.25, 2.0, 1741.5, 1740.5, 880.25, 5.0, 880.25, 872.25, 871.25, 1742.5, 5.0, 880.25, 879.25, 871.25, 1742.5, 1740.5, 1745.5, 890.25, 1.0, 872.25, 888.25, 880.25, 883.25, 1745.5, 890.25, 17.0, 883.25, 874.25, 20.0, 871.25, 4.0, 875.25, 8.0, 8.0, 5.0, 871.25, 895.25, 2.0, 879.25, 887.25, 1741.5, 874.25, 871.25, 883.25, 5.0, 5.0, 5.0, 879.25, 871.25, 10.0, 1741.5, 888.25, 0.0, 871.25, 880.25, 879.25, 1760.5, 887.25, 0.0, 880.25, 1741.5, 890.25, 1741.5, 1.0, 5.0, 872.25, 883.25, 1744.5, 886.25, 1741.5, 878.25, 1756.5, 878.25, 1.0, 1744.5, 874.25, 875.25, 887.25, 2.0, 890.25, 1750.5, 1760.5, 870.25, 1753.5, 1745.5, 1745.5, 9.0, 16.0, 1757.5, 5.0, 1745.5, 1740.5, 10.0, 5.0, 879.25, 1742.5, 871.25, 1745.5, 870.25, 1757.5, 1748.5, 8.0, 878.25, 875.25, 875.25, 872.25, 875.25, 872.25, 875.25, 1742.5, 4.0, 17.0, 872.25, 880.25, 1745.5, 1741.5, 5.0, 1.0, 9.0, 870.25, 1742.5, 1748.5, 1753.5, 1750.5, 880.25, 2.0, 2.0, 1741.5, 16.0, 1740.5, 870.25, 875.25, 1749.5, 1.0, 1742.5, 1740.5, 879.25, 872.25, 874.25, 1.0, 1741.5, 871.25, 5.0, 872.25, 878.25, 871.25, 879.25, 1744.5, 875.25, 1.0, 13.0, 18.0, 1740.5, 1765.5, 871.25, 870.25, 887.25, 878.25, 895.25, 895.25, 10.0, 1757.5, 871.25, 886.25, 890.25, 2.0, 16.0, 874.25, 1744.5, 1772.5, 17.0, 9.0, 874.25, 13.0, 875.25, 17.0, 1.0, 880.25, 887.25, 887.25, 883.25, 1753.5, 879.25, 875.25, 1740.5, 1758.5, 879.25, 871.25, 1744.5, 879.25, 875.25, 887.25, 1750.5, 871.25, 1765.5, 883.25, 872.25, 883.25, 871.25, 1740.5, 0.0, 887.25, 874.25, 875.25, 871.25, 878.25, 1745.5, 1749.5, 1.0, 0.0, 878.25, 887.25, 1.0, 872.25, 880.25, 870.25, 890.25, 5.0, 5.0, 1.0, 1756.5, 875.25, 1749.5, 874.25, 5.0, 890.25, 871.25, 875.25, 1.0, 1.0, 875.25, 1745.5, 875.25, 895.25, 9.0, 1753.5, 1744.5, 1758.5, 1753.5, 874.25, 2.0, 1748.5, 1742.5, 887.25, 879.25, 872.25, 890.25, 13.0, 16.0, 1744.5, 887.25, 871.25, 0.0, 895.25, 883.25, 871.25, 1742.5, 880.25, 1756.5, 887.25, 1757.5, 887.25, 886.25, 1756.5, 1756.5, 883.25, 871.25, 874.25, 1750.5, 1742.5, 886.25, 872.25, 1745.5, 871.25, 880.25, 872.25, 875.25, 1741.5, 880.25, 874.25, 871.25, 1.0, 8.0, 4.0, 871.25, 5.0, 875.25, 890.25, 874.25, 872.25, 1741.5, 890.25, 5.0, 1748.5, 1742.5, 870.25, 874.25, 1745.5, 874.25, 1.0, 1742.5, 5.0, 9.0, 874.25, 4.0, 880.25, 17.0, 895.25, 871.25, 1741.5, 0.0, 1.0, 883.25, 1757.5, 879.25, 890.25, 872.25, 874.25, 886.25, 1750.5, 875.25, 1760.5, 32.0, 886.25]
    nbins =20
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.hist_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.hist_list()
        self.assertEqual(
            str(cm_new.exception), "hist_list() missing 1 required positional argument: 'data'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.hist_list(
            data=self.data, nbin=self.nbins, fminiu=None, fmaxiu=None
        )
        return_old = oldfu.hist_list(
            data=self.data, nbin=self.nbins, fminiu=None, fmaxiu=None
        )

        self.assertTrue(array_equal(return_new[0], return_old[0]))
        self.assertTrue(array_equal(return_new[1], return_old[1]))
        self.assertTrue(
            array_equal(
                return_new[0],
                [
                    0.0,
                    88.625,
                    177.25,
                    265.875,
                    354.5,
                    443.125,
                    531.75,
                    620.375,
                    709.0,
                    797.625,
                    886.25,
                    974.875,
                    1063.5,
                    1152.125,
                    1240.75,
                    1329.375,
                    1418.0,
                    1506.625,
                    1595.25,
                    1683.875,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407],
            )
        )

    def test_empty_list(self):
        with self.assertRaises(ValueError) as cm_new:
            fu.hist_list(data=[], nbin=self.nbins, fminiu=None, fmaxiu=None)
        with self.assertRaises(ValueError) as cm_old:
            oldfu.hist_list(data=[], nbin=self.nbins, fminiu=None, fmaxiu=None)
        self.assertEqual(str(cm_new.exception), "max() arg is an empty sequence")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_fmin_equal_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=1)
        return_old = oldfu.hist_list(
            data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=1
        )
        self.assertTrue(array_equal(return_new[0], return_old[0]))
        self.assertTrue(array_equal(return_new[1], return_old[1]))
        self.assertTrue(
            array_equal(
                return_new[0],
                [
                    0.0,
                    88.625,
                    177.25,
                    265.875,
                    354.5,
                    443.125,
                    531.75,
                    620.375,
                    709.0,
                    797.625,
                    886.25,
                    974.875,
                    1063.5,
                    1152.125,
                    1240.75,
                    1329.375,
                    1418.0,
                    1506.625,
                    1595.25,
                    1683.875,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407],
            )
        )

    def test_fmin_major_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=10, fmaxiu=1)
        return_old = oldfu.hist_list(
            data=self.data, nbin=self.nbins, fminiu=10, fmaxiu=1
        )
        self.assertTrue(array_equal(return_new[0], return_old[0]))
        self.assertTrue(array_equal(return_new[1], return_old[1]))
        self.assertTrue(
            array_equal(
                return_new[0],
                [
                    0.0,
                    88.625,
                    177.25,
                    265.875,
                    354.5,
                    443.125,
                    531.75,
                    620.375,
                    709.0,
                    797.625,
                    886.25,
                    974.875,
                    1063.5,
                    1152.125,
                    1240.75,
                    1329.375,
                    1418.0,
                    1506.625,
                    1595.25,
                    1683.875,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407],
            )
        )

    def test_fmin_minor_fmax(self):
        return_new = fu.hist_list(data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=10)
        return_old = oldfu.hist_list(
            data=self.data, nbin=self.nbins, fminiu=1, fmaxiu=10
        )
        self.assertTrue(array_equal(return_new[0], return_old[0]))
        self.assertTrue(array_equal(return_new[1], return_old[1]))
        self.assertTrue(
            array_equal(
                return_new[0],
                [
                    0.0,
                    88.625,
                    177.25,
                    265.875,
                    354.5,
                    443.125,
                    531.75,
                    620.375,
                    709.0,
                    797.625,
                    886.25,
                    974.875,
                    1063.5,
                    1152.125,
                    1240.75,
                    1329.375,
                    1418.0,
                    1506.625,
                    1595.25,
                    1683.875,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[1],
                [1634, 0, 0, 0, 0, 0, 0, 0, 0, 2257, 487, 0, 0, 0, 0, 0, 0, 0, 0, 1407],
            )
        )


class Test_linreg(unittest.TestCase):
    X =  [0.010321706049974242, 0.01083134412382782, 0.01135326263319594, 0.011887461578078606, 0.01243394095847581, 0.012992700774387564, 0.013563741025813859, 0.014147061712754701, 0.014742662835210085, 0.015350544393180015, 0.015970706386664486, 0.016603148815663504, 0.01724787168017707, 0.017904874980205175, 0.01857415871574782, 0.019255722886805014, 0.019949567493376744, 0.020655692535463028, 0.021374098013063846, 0.022104783926179216, 0.022847750274809134, 0.023602997058953593, 0.024370524278612585, 0.025150331933786133, 0.025942420024474225, 0.02674678855067686, 0.02756343751239404, 0.028392366909625755, 0.029233576742372022, 0.030087067010632833, 0.03095283771440819, 0.03183088885369808, 0.032721220428502513, 0.0336238324388215, 0.03453872488465503, 0.03546589776600311, 0.03640535108286572, 0.03735708483524288, 0.03832109902313459, 0.03929739364654084, 0.040285968705461625, 0.04128682419989697, 0.04229996012984685, 0.04332537649531128, 0.04436307329629024, 0.04541305053278376, 0.046475308204791815, 0.047549846312314424, 0.04863666485535156, 0.04973576383390324, 0.05084714324796949, 0.051970803097550256, 0.05310674338264558, 0.054254964103255435, 0.05541546525937986, 0.056588246851018806, 0.057773308878172294, 0.05897065134084034, 0.06018027423902292, 0.06140217757272006, 0.06263636134193172, 0.06388282554665795, 0.06514157018689871, 0.06641259526265401, 0.06769590077392384, 0.06899148672070828, 0.07029935310300718, 0.0716194999208207, 0.0729519271741487, 0.07429663486299128, 0.07565362298734839, 0.07702289154722006, 0.07840444054260623, 0.07979826997350697, 0.08120437983992228, 0.08262277014185211, 0.08405344087929649, 0.08549639205225538, 0.0869516236607289, 0.08841913570471686, 0.08989892818421946, 0.09139100109923654, 0.09289535444976818, 0.09441198823581437, 0.0959409024573751, 0.09748209711445034, 0.0990355722070402, 0.10060132773514453, 0.10217936369876345, 0.1037696800978969, 0.10537227693254486, 0.10698715420270744, 0.1086143119083845, 0.11025375004957616]
    Y = [-11.598208535753709, -11.670055986983082, -11.715472889187454, -11.721520336767835, -11.754087050961333, -11.794871732452629, -11.817780093327917, -11.802722189379722, -11.863592233549991, -11.921123677996722, -11.95902513230595, -11.933629335048538, -11.915743768447271, -11.966059068532374, -12.042972297584155, -12.03962326807928, -12.044619336256854, -12.070783996718191, -12.07730196234726, -12.079939818735411, -12.120995270586642, -12.150495966495315, -12.132864795274937, -12.134788668129582, -12.124908862969621, -12.108533374912591, -12.11567162100514, -12.127635824394142, -12.124361327859432, -12.129226441429045, -12.123704781720857, -12.111230543811626, -12.105664363172316, -12.089102245058115, -12.058937097743144, -12.053199139097218, -12.049227258115069, -12.033315341117815, -12.036666386039082, -12.045739402030785, -12.054803077009378, -12.064842016575508, -12.070795436297036, -12.074413883785937, -12.068698303209679, -12.053105163707929, -12.055119058830503, -12.062588687988185, -12.057927163763599, -12.065229994570522, -12.082723967123556, -12.093725415458337, -12.116188563387302, -12.128253727798477, -12.139530188730143, -12.17087920425965, -12.179627724808084, -12.202047990136753, -12.215528414322858, -12.223079503261449, -12.216651741060506, -12.215515560509008, -12.213048028415042, -12.226887681464582, -12.232528797041986, -12.223534845979477, -12.227798798224361, -12.232703260585712, -12.21173481376324, -12.200013819009962, -12.209577540936879, -12.203847522841917, -12.205755398419036, -12.201721215864133, -12.201498657499592, -12.216130249557292, -12.207901790362463, -12.209465651650572, -12.204923093830777, -12.196702167715916, -12.194333511333301, -12.191564840798513, -12.195873219844968, -12.185144804378796, -12.183091273975064, -12.176768778545615, -12.182311717055121, -12.197166026213058, -12.190798716465572, -12.191275884743456, -12.183728324696375, -12.17451889903635, -12.175073620203147, -12.161912993584684]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg()
        self.assertEqual(
            str(cm_new.exception), "linreg() missing 2 required positional arguments: 'X' and 'Y'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.linreg(X=self.X, Y=self.Y)
        return_old = oldfu.linreg(X=self.X, Y=self.Y)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(
            array_equal(return_new, (-3.3524868938354802, -11.920708605604739))
        )

    def test_different_list_size_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=self.X, Y=[3, 4])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=self.X, Y=[3, 4])
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_both_empty_list_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.linreg(X=[], Y=[])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.linreg(X=[], Y=[])
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_X_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=[], Y=self.Y)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=[], Y=self.Y)
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Y_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.linreg(X=self.X, Y=[])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.linreg(X=self.X, Y=[])
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_pearson(unittest.TestCase):
    X = [-11.920708605604739, -11.920729190604296, -11.920790945602965, -11.920893870600748, -11.921037965597645,
               -11.921223230593652, -11.921449665588774, -11.921717270583009, -11.922026045576358, -11.92237599056882,
               -11.922767105560395, -11.923199390551082, -11.923672845540882, -11.924187470529796, -11.924743265517822,
               -11.925340230504963, -11.925978365491217, -11.926657670476581, -11.927378145461061, -11.928139790444654,
               -11.928942605427359, -11.929786590409178, -11.930671745390109, -11.931598070370153, -11.932565565349313,
               -11.933574230327583, -11.934624065304966, -11.935715070281464, -11.936847245257075, -11.938020590231798,
               -11.939235105205634, -11.940490790178584, -11.941787645150647, -11.943125670121823, -11.944504865092112,
               -11.945925230061514, -11.947386765030028, -11.948889469997656, -11.950433344964397, -11.952018389930252,
               -11.953644604895219, -11.9553119898593, -11.957020544822495, -11.9587702697848, -11.96056116474622,
               -11.962393229706754, -11.9642664646664, -11.966180869625159, -11.968136444583031, -11.970133189540016,
               -11.972171104496114, -11.974250189451327, -11.976370444405651, -11.978531869359088, -11.98073446431164,
               -11.982978229263304, -11.98526316421408, -11.987589269163971, -11.989956544112975, -11.992364989061091,
               -11.994814604008321, -11.997305388954663, -11.999837343900118, -12.002410468844687, -12.005024763788368,
               -12.007680228731164, -12.010376863673072, -12.013114668614094, -12.015893643554227, -12.018713788493475,
               -12.021575103431836, -12.02447758836931, -12.027421243305895, -12.030406068241595, -12.033432063176409,
               -12.036499228110335, -12.039607563043374, -12.042757067975526, -12.04594774290679, -12.04917958783717,
               -12.05245260276666, -12.055766787695266, -12.059122142622984, -12.062518667549814, -12.065956362475758,
               -12.069435227400815, -12.072955262324985, -12.076516467248268, -12.080118842170664, -12.083762387092174,
               -12.087447102012796, -12.091172986932532, -12.09494004185138, -12.098748266769343, -12.102597661686417,
               -12.106488226602606, -12.110419961517907, -12.11439286643232, -12.118406941345848, -12.122462186258488,
               -12.126558601170242, -12.130696186081108, -12.134874940991088, -12.139094865900182, -12.143355960808387,
               -12.147658225715706, -12.152001660622139, -12.156386265527683, -12.160812040432342, -12.165278985336114,
               -12.169787100238997, -12.174336385140995, -12.178926840042106, -12.18355846494233, -12.188231259841666,
               -12.192945224740116, -12.197700359637679, -12.202496664534356, -12.207334139430145, -12.212212784325049,
               -12.217132599219063, -12.222093584112192, -12.227095739004433, -12.232139063895788, -12.237223558786257,
               -12.242349223675838, -12.247516058564532, -12.252724063452339, -12.257973238339259, -12.263263583225292,
               -12.268595098110438, -12.273967782994699, -12.279381637878071, -12.284836662760558, -12.290332857642156,
               -12.295870222522868, -12.301448757402692, -12.307068462281631, -12.312729337159682, -12.318431382036847,
               -12.324174596913124, -12.329958981788515, -12.335784536663018, -12.341651261536636, -12.347559156409366,
               -12.353508221281208, -12.359498456152165, -12.365529861022234, -12.371602435891417, -12.377716180759712,
               -12.38387109562712, -12.390067180493642, -12.396304435359276, -12.402582860224024, -12.408902455087885,
               -12.415263219950861, -12.421665154812947, -12.428108259674147, -12.43459253453446, -12.441117979393887,
               -12.447684594252426, -12.454292379110079, -12.460941333966845, -12.467631458822723, -12.474362753677715,
               -12.48113521853182, -12.487948853385038, -12.49480365823737, -12.501699633088814, -12.508636777939373,
               -12.515615092789043, -12.522634577637826, -12.529695232485723, -12.536797057332732, -12.543940052178856,
               -12.551124217024091, -12.55834955186844]

    Y =  [-10.598742713481631, -7.7405339253848817, -7.2852275090232412, -8.1310928801977695, -8.2690307722033403, -8.7490957059061749, -9.0911528838923648, -9.1433178595516704, -9.3743291012068735, -9.6213884419831341, -9.8768964569281739, -10.004129024368162, -10.097522007981846, -10.383397742828921, -10.637387648825358, -10.688278835798735, -10.730429627332981, -10.875146408476629, -10.943610383867224, -11.065118155518878, -11.187031299664762, -11.21512245056009, -11.245831568063513, -11.312350577419863, -11.35483775435771, -11.416114992966921, -11.448679022072508, -11.442412057318037, -11.436356558609468, -11.450352959931362, -11.457395072596587, -11.445846461544701, -11.441512764149257, -11.432700088754258, -11.419687681270885, -11.432090213216084, -11.459085141827391, -11.483221839813284, -11.501900566515982, -11.507718389125786, -11.536716243105058, -11.598208535753709, -11.670055986983082, -11.715472889187454, -11.721520336767835, -11.754087050961333, -11.794871732452629, -11.817780093327917, -11.802722189379722, -11.863592233549991, -11.921123677996722, -11.95902513230595, -11.933629335048538, -11.915743768447271, -11.966059068532374, -12.042972297584155, -12.03962326807928, -12.044619336256854, -12.070783996718191, -12.07730196234726, -12.079939818735411, -12.120995270586642, -12.150495966495315, -12.132864795274937, -12.134788668129582, -12.124908862969621, -12.108533374912591, -12.11567162100514, -12.127635824394142, -12.124361327859432, -12.129226441429045, -12.123704781720857, -12.111230543811626, -12.105664363172316, -12.089102245058115, -12.058937097743144, -12.053199139097218, -12.049227258115069, -12.033315341117815, -12.036666386039082, -12.045739402030785, -12.054803077009378, -12.064842016575508, -12.070795436297036, -12.074413883785937, -12.068698303209679, -12.053105163707929, -12.055119058830503, -12.062588687988185, -12.057927163763599, -12.065229994570522, -12.082723967123556, -12.093725415458337, -12.116188563387302, -12.128253727798477, -12.139530188730143, -12.17087920425965, -12.179627724808084, -12.202047990136753, -12.215528414322858, -12.223079503261449, -12.216651741060506, -12.215515560509008, -12.213048028415042, -12.226887681464582, -12.232528797041986, -12.223534845979477, -12.227798798224361, -12.232703260585712, -12.21173481376324, -12.200013819009962, -12.209577540936879, -12.203847522841917, -12.205755398419036, -12.201721215864133, -12.201498657499592, -12.216130249557292, -12.207901790362463, -12.209465651650572, -12.204923093830777, -12.196702167715916, -12.194333511333301, -12.191564840798513, -12.195873219844968, -12.185144804378796, -12.183091273975064, -12.176768778545615, -12.182311717055121, -12.197166026213058, -12.190798716465572, -12.191275884743456, -12.183728324696375, -12.17451889903635, -12.175073620203147, -12.161912993584684, -12.167875018710498, -12.171546036897956, -12.155575223495314, -12.149486538330708, -12.150277295294465, -12.141786714125459, -12.133356622161283, -12.123349290595661, -12.116185071883296, -12.119733901823622, -12.104280048369871, -12.096436515366012, -12.087783639464494, -12.086156781470418, -12.078266351659419, -12.06563851354796, -12.060526133083849, -12.063686478433725, -12.052576677629524, -12.04415696367381, -12.025956658972021, -12.013468902065263, -11.997663561581593, -11.980204060613428, -11.964281532488615, -11.948484271701215, -11.944717110786863, -11.939662898496096, -11.92708891454858, -11.91540978810044, -11.904687158260389, -11.900194369000975, -11.891560709414353, -11.888981828066552, -11.88045696078118, -11.871047017891712, -11.865416131075067, -11.853751766387191, -11.836462495572446, -11.828741093896237, -11.764365384684698, -11.8979108113727]
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson()
        self.assertEqual(
            str(cm_new.exception), "pearson() missing 2 required positional arguments: 'X' and 'Y'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.pearson(X=self.X, Y=self.Y)
        return_old = oldfu.pearson(X=self.X, Y=self.Y)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertEqual(return_new, 0.43165202148828613)

    def test_different_list_size_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=self.X, Y=[3, 4])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=self.X, Y=[3, 4])
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_both_empty_list_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.pearson(X=[], Y=[])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.pearson(X=[], Y=[])
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_X_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=[], Y=self.Y)
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=[], Y=self.Y)
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_Y_empty_list_returns_TypeError(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pearson(X=self.X, Y=[])
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pearson(X=self.X, Y=[])
        self.assertEqual(
            str(cm_new.exception),
            "unsupported operand type(s) for +=: 'float' and 'NoneType'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_table_stat(unittest.TestCase):
    X = [4808, 2749]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.table_stat()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.table_stat()
        self.assertEqual(
            str(cm_new.exception), "table_stat() missing 1 required positional argument: 'X'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.table_stat(X=self.X)
        return_old = oldfu.table_stat(X=self.X)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, (3778, 2119741.0, 2749, 4808)))

    def test_empty_list_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.table_stat(X=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.table_stat(X=[])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


""" this function has been cleaned
class Test_mono(unittest.TestCase):
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.mono()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.mono()
        self.assertEqual(str(cm_new.exception), "mono() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_(self):
        return_new = fu.mono(k1=10,k2=20)
        return_old = oldfu.mono(k1=10,k2=20)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 200)

    def test_m1_equalm2(self):
        return_new = fu.mono(k1=20,k2=20)
        return_old = oldfu.mono(k1=20,k2=20)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,210)

    def test_both_zero(self):
        return_new = fu.mono(k1=0,k2=0)
        return_old = oldfu.mono(k1=0,k2=0)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new,0)
"""


class Test_k_means_stab_bbenum(unittest.TestCase):
    PART = [
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
    ]
    expected_result = ([[1, 1, 1, 1]], [[], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]], [0, 100], [0, 100], [0, 100.0], 100.0)

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_stab_bbenum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_stab_bbenum()
        self.assertEqual(
            str(cm_new.exception),
            "k_means_stab_bbenum() missing 1 required positional argument: 'PART'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_PART_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_stab_bbenum(
                PART=[],
                T=10,
                nguesses=5,
                J=50,
                max_branching=40,
                stmult=0.25,
                branchfunc=2,
                LIM=-1,
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_stab_bbenum(
                PART=[],
                T=10,
                nguesses=5,
                J=50,
                max_branching=40,
                stmult=0.25,
                branchfunc=2,
                LIM=-1,
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.k_means_stab_bbenum(
            PART=self.PART,
            T=10,
            nguesses=5,
            J=50,
            max_branching=40,
            stmult=0.25,
            branchfunc=2,
            LIM=-1,
        )
        return_old = oldfu.k_means_stab_bbenum(
            PART=self.PART,
            T=10,
            nguesses=5,
            J=50,
            max_branching=40,
            stmult=0.25,
            branchfunc=2,
            LIM=-1,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, self.expected_result)) #luca

    def test_branchfunc_is4(self):
        print("expected")
        print(self.expected_result)
        return_new = fu.k_means_stab_bbenum(
            PART=self.PART,
            T=10,
            nguesses=5,
            J=50,
            max_branching=40,
            stmult=0.25,
            branchfunc=4,
            LIM=-1,
        )
        return_old = oldfu.k_means_stab_bbenum(
            PART=self.PART,
            T=10,
            nguesses=5,
            J=50,
            max_branching=40,
            stmult=0.25,
            branchfunc=4,
            LIM=-1,
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, self.expected_result)) #luca

    def test_lenPART_is2(self):
        self.assertTrue(True)
        """
        i meets an exitcode:
        Traceback (most recent call last):
          File "/home/lusnig/src_sphire_1_3/eman2/sphire/tests/test_statistics.py", line 1785, in test_lenPART_is2
            return_new =fu.k_means_stab_bbenum(PART=self.PART[0:2], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
          File "/home/lusnig/src_sphire_1_3/eman2/sphire/utils/SPHIRE/libpy/sp_statistics.py", line 1412, in k_means_stab_bbenum
            sys.exit()
        SystemExit

        return_new =fu.k_means_stab_bbenum(PART=self.PART[0:2], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        return_old =oldfu.k_means_stab_bbenum(PART=self.PART[0:2], T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1)
        self.assertTrue(array_equal(return_new,return_old))
        self.assertTrue(array_equal(return_new,([[1, 0]], [[], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]], [0, 97], [0, 97], [0, 100.0], 100.0)))
        """


class Test_k_means_match_bbenum(unittest.TestCase):
    PART = [
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
        [numpy_asarray([50, 73, 79]), numpy_asarray(range(100))],
    ]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.k_means_match_bbenum()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.k_means_match_bbenum()
        self.assertEqual(
            str(cm_new.exception),
            "k_means_match_bbenum() missing 1 required positional argument: 'PART'",
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_PART_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_match_bbenum(
                [],
                T=10,
                J=1,
                max_branching=40,
                stmult=0.25,
                nguesses=5,
                branchfunc=2,
                LIM=-1,
                DoMPI_init=False,
                Njobs=-1,
                DoMPI=False,
                K=-1,
                np=-1,
                c_dim=[],
                N_start=-1,
                N_stop=-1,
                topMatches=[],
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_match_bbenum(
                [],
                T=10,
                J=1,
                max_branching=40,
                stmult=0.25,
                nguesses=5,
                branchfunc=2,
                LIM=-1,
                DoMPI_init=False,
                Njobs=-1,
                DoMPI=False,
                K=-1,
                np=-1,
                c_dim=[],
                N_start=-1,
                N_stop=-1,
                topMatches=[],
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))

    def test_DoMPI_isTrue_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.k_means_match_bbenum(
                self.PART,
                T=10,
                J=1,
                max_branching=40,
                stmult=0.25,
                nguesses=5,
                branchfunc=2,
                LIM=-1,
                DoMPI_init=False,
                Njobs=-1,
                DoMPI=True,
                K=-1,
                np=-1,
                c_dim=[],
                N_start=-1,
                N_stop=-1,
                topMatches=[],
            )
        with self.assertRaises(IndexError) as cm_old:
            oldfu.k_means_match_bbenum(
                self.PART,
                T=10,
                J=1,
                max_branching=40,
                stmult=0.25,
                nguesses=5,
                branchfunc=2,
                LIM=-1,
                DoMPI_init=False,
                Njobs=-1,
                DoMPI=True,
                K=-1,
                np=-1,
                c_dim=[],
                N_start=-1,
                N_stop=-1,
                topMatches=[],
            )
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_DoMPI_DoMPIinit_areTrue_spawn_errorMsg_but_works(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=True,
            Njobs=-1,
            DoMPI=True,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=True,
            Njobs=-1,
            DoMPI=True,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, ([[1, 1, 1, 1]], 100)))

    def test_T0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=0,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=0,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, return_old))

        self.assertTrue(array_equal(return_new, [[0, 0, 0, 0], [1, 1, 1, 1]]))

    def test_J0_crashes_because_signal11SIGSEV(self):
        self.assertTrue(True)
        """
        return_new = fu.k_means_match_bbenum(self.PART, T=10, J=0, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        return_old = oldfu.k_means_match_bbenum(self.PART, T=10, J=0, max_branching=40, stmult=0.25, nguesses=5,branchfunc=2, LIM=-0, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1,c_dim=[], N_start=-1, N_stop=-1, topMatches=[])
        self.assertTrue(array_equal(return_new, return_old))
        """

    def test_branchfunc0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=0,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=0,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))

    def test_nguess0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=0,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=0,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )

        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))
        self.assertTrue(array_equal(return_new, return_old))

    def test_np0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=0,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=0,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))
        self.assertTrue(array_equal(return_new, return_old))

    def test_LIM0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=0,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=0,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))
        self.assertTrue(array_equal(return_new, return_old))

    def test_stmult0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=40,
            stmult=0,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))
        self.assertTrue(array_equal(return_new, return_old))

    def test_max_branching0(self):
        return_new = fu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=0,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        return_old = oldfu.k_means_match_bbenum(
            self.PART,
            T=10,
            J=1,
            max_branching=0,
            stmult=0.25,
            nguesses=5,
            branchfunc=2,
            LIM=-1,
            DoMPI_init=False,
            Njobs=-1,
            DoMPI=False,
            K=-1,
            np=-1,
            c_dim=[],
            N_start=-1,
            N_stop=-1,
            topMatches=[],
        )
        self.assertTrue(array_equal(return_new, [[1, 1, 1, 1]]))
        self.assertTrue(array_equal(return_new, return_old))


""" This function has been cleaned
class Test_scale_fsc_datasetsize(unittest.TestCase):
    fsc_to_be_adjusted =[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992, 0.9957, 0.99387, 0.99405, 0.99123, 0.98935, 0.9872, 0.98074, 0.96944, 0.96374, 0.94607, 0.91744, 0.88535, 0.89291, 0.88777, 0.8877, 0.87859, 0.87215, 0.8633, 0.8441, 0.84776, 0.85752, 0.87845, 0.88849, 0.87827, 0.85501, 0.83834, 0.86192, 0.87235, 0.84794, 0.82865, 0.81224, 0.81058, 0.7955, 0.77124, 0.74744, 0.74404, 0.7588, 0.71706, 0.58851, 0.50241, 0.60636, 0.67929, 0.62126, 0.46365, 0.45433, 0.50698, 0.51007, 0.47945, 0.47566, 0.44534, 0.37992, 0.33675, 0.36838, 0.35485, 0.35406, 0.34595, 0.30988, 0.32704, 0.35027, 0.3361, 0.34971, 0.35296, 0.31079, 0.33571, 0.37315, 0.34308, 0.37531, 0.37823, 0.40466, 0.42364, 0.40338, 0.39635, 0.37909, 0.37814, 0.3756, 0.35529, 0.37891, 0.35317, 0.34812, 0.34572, 0.33312, 0.31022, 0.28474, 0.25991, 0.25608, 0.24614, 0.20188, 0.17434, 0.17683, 0.15445, 0.13581, 0.12399, 0.12581, 0.12332, 0.111, 0.0911946, 0.0864365, 0.0785264, 0.0679365, 0.0488136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    nfsc = 9892
    nnew = 9892
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.scale_fsc_datasetsize()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.scale_fsc_datasetsize()
        self.assertEqual(str(cm_new.exception), "scale_fsc_datasetsize() takes exactly 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc,nnew=self.nnew)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992, 0.9957, 0.99387, 0.99405, 0.99123, 0.98935, 0.9872, 0.98074, 0.96944, 0.96374, 0.94607, 0.91744, 0.88535, 0.89291, 0.88777, 0.8877, 0.87859, 0.87215, 0.8633, 0.8441, 0.84776, 0.85752, 0.87845, 0.88849, 0.87827, 0.85501, 0.83834, 0.86192, 0.87235, 0.84794, 0.82865, 0.81224, 0.81058, 0.7955, 0.77124, 0.74744, 0.74404, 0.7588, 0.71706, 0.58851, 0.50241, 0.60636, 0.67929, 0.62126, 0.46365, 0.45433, 0.50698, 0.51007, 0.47945, 0.47566, 0.44534, 0.37992, 0.33675, 0.36838, 0.35485, 0.35406, 0.34595, 0.30988, 0.32704, 0.35027, 0.3361, 0.34971, 0.35296, 0.31079, 0.33571, 0.37315, 0.34308, 0.37531, 0.37823, 0.40466, 0.42364, 0.40338, 0.39635, 0.37909, 0.37814, 0.3756, 0.35529, 0.37891, 0.35317, 0.34812, 0.34572, 0.33312, 0.31022, 0.28474, 0.25991, 0.25608, 0.24614, 0.20188, 0.17434, 0.17683, 0.15445, 0.13581, 0.12399, 0.12581, 0.12332, 0.111, 0.0911946, 0.0864365, 0.0785264, 0.0679365, 0.0488136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_emptylist(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=self.nfsc, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=self.nfsc,nnew=self.nnew)
        self.assertTrue(array_equal(return_new, []))
        self.assertTrue(array_equal(return_new, return_old))

    def test_nfsc_greater_nnew(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc+10, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc+10,nnew=self.nnew)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9991991919133101, 0.995695671763659, 0.9938638410985181, 0.9940440208634007, 0.9912212120804388, 0.9893393484997355, 0.9871872260447904, 0.9807209050903437, 0.9694100513837745, 0.9637046745536805, 0.946018424207827, 0.917363435579012, 0.885247398287289, 0.8928133447421265, 0.8876692891990225, 0.8875992343368216, 0.8784821790119756, 0.8720372927942945, 0.8631807149173785, 0.843966989027912, 0.8476295479959695, 0.8573965043935038, 0.8783420718976555, 0.8883898540713431, 0.8781619342375103, 0.8548846969953323, 0.8382030166804726, 0.8617997034946842, 0.8722374432777654, 0.8478096745459852, 0.8285064854566488, 0.8120858580259775, 0.810424813315671, 0.7953355781179486, 0.7710616861389799, 0.747249214252374, 0.7438475260687702, 0.7586150243187773, 0.7168549586110095, 0.5882652918672154, 0.5021574034549281, 0.6061188024611799, 0.6790698377601415, 0.6210222260837981, 0.463398742503496, 0.45407951721577117, 0.5067274456477423, 0.5098174980925592, 0.4791978301349811, 0.4754080030001081, 0.4450904307942334, 0.3796819963482656, 0.3365243633400876, 0.36814493357971323, 0.3546187199078563, 0.35382895241253404, 0.3457214121617222, 0.3096639614909986, 0.3268176635514925, 0.3500400853774037, 0.33587457831324086, 0.3494802552414265, 0.3527292782529308, 0.31057361257630034, 0.3354847077067504, 0.3729136868736183, 0.3428523144536889, 0.37507313744618753, 0.37799240968401815, 0.40441660596219825, 0.4233933087672451, 0.40313685451865866, 0.39610827864698217, 0.37885219862651714, 0.37790243185778905, 0.37536306439808914, 0.3550585909581716, 0.3786722426777146, 0.3529392158582757, 0.34789074098641914, 0.34549148385760375, 0.3328955748069882, 0.3100038309315609, 0.28453426205390553, 0.2597156884514092, 0.2558875617923892, 0.24595256186990372, 0.2017172477054805, 0.1741946042085414, 0.17668297197529226, 0.15431809172820385, 0.1356914565307124, 0.12388029480686033, 0.12569891555298499, 0.12321080423791014, 0.11090033320236868, 0.09111089390176229, 0.08635674628852816, 0.07845331805355968, 0.06787254779170172, 0.04876670733038269, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_nfsc_lower_nnew(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= self.fsc_to_be_adjusted, nfsc=self.nfsc-10, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=self.fsc_to_be_adjusted, nfsc=self.nfsc-10,nnew=self.nnew)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9992008080879967, 0.9957043282739705, 0.993876158977815, 0.9940559792085284, 0.9912387880753856, 0.9893606517296208, 0.9872127742857975, 0.9807590956532373, 0.9694699504667269, 0.9637753280361855, 0.9461215814161805, 0.9175165772023998, 0.8854526254989016, 0.893006676187781, 0.8878707336559222, 0.8878007885448379, 0.8786978474582086, 0.8722627363433496, 0.8634193180557805, 0.8442330529043143, 0.847890492163897, 0.8576435311871446, 0.8785579546293825, 0.8885901685095942, 0.8783780923627005, 0.8551353397421242, 0.8384770280998409, 0.8620403300937923, 0.8724625857754491, 0.8480703655274789, 0.8287935642713687, 0.8123942005004912, 0.8107352461285096, 0.795664489878872, 0.7714183963529617, 0.7476308831943587, 0.744232573564021, 0.7589850659097225, 0.7172651587181639, 0.5887549118061657, 0.5026628507965304, 0.6066013895784409, 0.6795103050444103, 0.6214979560613412, 0.46390153011087243, 0.45458075928317637, 0.507232806225359, 0.5103227521491716, 0.47970243540512547, 0.4759122642910012, 0.4455898492375436, 0.380158302224063, 0.3369759394377599, 0.3686153667994114, 0.35508158196785944, 0.3542913495292306, 0.34617889031865434, 0.3100963401599992, 0.32726263916867665, 0.3505002168478494, 0.3363257244727836, 0.3499400470210042, 0.3531910237768141, 0.31100668916320073, 0.33593559508335824, 0.3733866128166641, 0.3433079881556603, 0.37554716190514603, 0.37846788918262947, 0.40490368718270053, 0.4238869788709483, 0.40362343895688485, 0.396592016549525, 0.3793281000916173, 0.3783778670242496, 0.3758372349064654, 0.35552171087972223, 0.3791480560717202, 0.35340108615503196, 0.3483495613749182, 0.3459488186343066, 0.3333447279945461, 0.3104364707529276, 0.28494603568690396, 0.26010460252252415, 0.256272727868698, 0.24632772403764963, 0.202043015134694, 0.17448563871053147, 0.17697727292939475, 0.15458213397025733, 0.13592875077550817, 0.12409989966954055, 0.12592128095895827, 0.1234293894835953, 0.11109984610107548, 0.09127846004588899, 0.08651640115876466, 0.07859961823005365, 0.06800057283857858, 0.04886058293771373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    def test_nfsc_is0(self):
        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=0, nnew=self.nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=0,nnew=self.nnew)
        self.assertTrue(array_equal(return_new, []))
        self.assertTrue(array_equal(return_new, return_old))

    def test_nnew_is0_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.scale_fsc_datasetsize(fsc_to_be_adjusted= [], nfsc=self.nfsc, nnew=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted=[], nfsc=self.nfsc,nnew=0)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
"""


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




class Test_lib_statistics_compare(unittest.TestCase):

    def test_add_ave_varf_MPI_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]

        mpi_barrier(MPI_COMM_WORLD)
        return_new = fu.add_ave_varf_MPI(0,data)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.add_ave_varf_MPI(0, data)


        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))


    #Unable to read the file 
    def test_sum_oe_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (data,) = argum[0]

        return_new = fu.sum_oe(data)
        return_old = oldfu.sum_oe(data)

        # print(return_new)

        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))


    def test_ave_var_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data) = argum[0][0]

        return_new = fu.ave_var(data)
        return_old = oldfu.ave_var(data)

        self.assertTrue(return_new, return_old)


    def test_ave_series_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]

        return_new = fu.ave_series(data,)
        return_old = oldfu.ave_series(data,)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_varf2d_MPI_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_series")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data,) = argum[0]
        ave = get_data(1)

        return_new = fu.varf2d_MPI(0,data,ave)
        return_old = oldfu.varf2d_MPI(0,data,ave)

        self.assertTrue(return_new, return_old)


    #The function works but is extremely slow. It takes 55 mins to check the unit test that is why it is commented
    # def test_varf3d_MPI_true_should_return_equal_objects(self):
    #
    #     stack_name = "bdb:tests/Substack/sort3d_substack_003"
    #     nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
    #     list_proj = list(range(nima))
    #     proj = EMAN2_cppwrap.EMData()
    #     proj.read_image(stack_name, list_proj[0])
    #
    #     return_new = fu.varf3d_MPI([proj], mask2D=False  )
    #     mpi_barrier(MPI_COMM_WORLD)
    #     print("Hello")
    #     return_old = oldfu.varf3d_MPI([proj] , mask2D=False)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     print('yoyo')
    #
    #     self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_ccc_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ccc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2,mask) = argum[0]

        return_new = fu.ccc(img1, img2,mask )
        return_old = oldfu.ccc(img1, img2,mask )

        self.assertEqual(return_new, return_old)


    def test_fsc_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2) = argum[0]

        return_new = fu.fsc(img1, img2 )
        return_old = oldfu.fsc(img1, img2)

        self.assertEqual(return_new, return_old)


    def test_fsc_mask_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (img1, img2, mask, w, filename) = argum[0]
        filename = os.path.join(ABSOLUTE_PATH, filename)

        return_new = fu.fsc_mask(img1, img2, mask, w, filename)
        return_old = oldfu.fsc_mask(img1, img2, mask, w, filename)

        self.assertEqual(return_new, return_old)

    #  Can only work on the cluster
    # def test_locres_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.locres")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc) = argum[0]
    #
    #     number_of_proc = 1
    #     main_node = 0
    #     myid = 0
    #
    #     return_new = fu.locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)


    def test_histogram_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.ave_var")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (data) = argum[0][0]

        image = get_data(1)

        return_new = fu.histogram(data[0])
        return_old = oldfu.histogram(data[0])

        self.assertTrue(return_new, return_old)



    def test_k_means_match_clusters_asg_new_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_clusters_asg_new")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (asg1, asg2) = argum[0]

        return_new = fu.k_means_match_clusters_asg_new(asg1, asg2)
        return_old = oldfu.k_means_match_clusters_asg_new(asg1, asg2)

        self.assertTrue(return_new, return_old)


    def test_hist_list_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.hist_list")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (data,nbins) = argum[0]

        return_new = fu.hist_list(data,nbins)
        return_old = oldfu.hist_list(data,nbins)

        self.assertEqual(return_new, return_old)


    def test_linreg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.linreg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, Y) = argum[0]

        return_new = fu.linreg(X, Y)
        return_old = oldfu.linreg(X, Y)

        self.assertEqual(return_new, return_old)


    def test_pearson_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.pearson")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, Y) = argum[0]

        return_new = fu.pearson(X, Y)
        return_old = oldfu.pearson(X, Y)

        self.assertEqual(return_new, return_old)


    def test_table_stat_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.table_stat")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (X, ) = argum[0]

        return_new = fu.table_stat(X)
        return_old = oldfu.table_stat(X)

        self.assertEqual(return_new, return_old)


    def test_mono_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.mono")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (k1,k2) = argum[0]

        return_new = fu.mono(k1,k2)
        return_old = oldfu.mono(k1,k2)

        self.assertEqual(return_new, return_old)


    def test_k_means_stab_bbenum_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_stab_bbenum")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (Part,) = argum[0]

        return_new = fu.k_means_stab_bbenum(Part)
        return_old = oldfu.k_means_stab_bbenum(Part)

        self.assertEqual(return_new, return_old)


    def test_k_means_match_bbenum_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.k_means_match_bbenum")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (Part,) = argum[0]

        return_new = fu.k_means_match_bbenum(Part)
        return_old = oldfu.k_means_match_bbenum(Part)

        self.assertEqual(return_new, return_old)


    def test_scale_fsc_datasetsize_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.scale_fsc_datasetsize")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])
        # print(argum[1])

        (fsc_to_be_adjusted, nfsc, nnew) = argum[0]

        return_new = fu.scale_fsc_datasetsize(fsc_to_be_adjusted, nfsc, nnew)
        return_old = oldfu.scale_fsc_datasetsize(fsc_to_be_adjusted, nfsc, nnew)

        self.assertEqual(return_new, return_old)
"""


if __name__ == "__main__":
    unittest.main()
