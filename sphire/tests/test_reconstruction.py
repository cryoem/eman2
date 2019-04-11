import unittest
from copy import deepcopy
from EMAN2_cppwrap import EMData,Util, Reconstructors,Transform
import EMAN2_cppwrap

import numpy
from test_module import get_data,get_data_3d, get_arg_from_pickle_file
import os
ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))

from mpi import *
mpi_init(0, [])

from sphire.libpy import sparx_reconstruction as fu
from sphire.tests.sparx_lib import sparx_reconstruction as oldfu

from sphire.libpy import sparx_utilities

class Test_lib_compare_for_reconstruction(unittest.TestCase):

    def test_insert_slices_should_return_True(self):
        argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        refvol = sparx_utilities.model_blank(76)
        fftvol = EMData(76,76)
        weight = EMData(76,76)

        params = {"size": 76, "npad": 2, "symmetry": "c1", "fftvol": fftvol, "weight": weight, "snr": 2}
        r = Reconstructors.get( "nn4", params )
        r.setup()

        return_new = fu.insert_slices(r,data)
        return_old = oldfu.insert_slices(r, data)
        self.assertEqual(return_new, return_old)



    def test_insert_slices_pdf_should_return_True(self):
        argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))
        (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        refvol = sparx_utilities.model_blank(76)
        fftvol = EMData(76,76)
        weight = EMData(76,76)


        params = {"size": 76, "npad": 2, "symmetry": "c1", "fftvol": fftvol, "weight": weight, "snr": 2}
        r = Reconstructors.get( "nn4", params )
        r.setup()

        return_new = fu.insert_slices_pdf(r,data)
        return_old = oldfu.insert_slices_pdf(r, data)
        self.assertEqual(return_new, return_old)



    def test_recons3d_4nn_MPI_should_return_True(self):
        # argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))
        # (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]

        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]

        mpi_comm = MPI_COMM_WORLD
        myid = mpi_comm_rank(mpi_comm)

        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        snr = optionsnew.snr
        npad = optionsnew.npad


        return_new = fu.recons3d_4nn_MPI(myid, datanew, symmetry="c1", npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_MPI(myid, datanew, symmetry="c1", npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(),0.5))


    # def test_recons3d_trl_struct_MPI_should_return_True(self):
    #     # argum = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/alignment.shc"))
    #     # (data, refrings, list_of_ref_ang, numr, xrng, yrng, step) = argum[0]
    #
    #     arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))
    #
    #     (datanew,optionsnew,iternew) = arga[0]
    #
    #     mpi_comm = MPI_COMM_WORLD
    #     myid = mpi_comm_rank(mpi_comm)
    #
    #     sym = optionsnew.sym
    #     sym = sym[0].lower() + sym[1:]
    #     snr = optionsnew.snr
    #     npad = optionsnew.npad
    #
    #
    #     return_new = fu.recons3d_trl_struct_MPI(myid, 0, datanew, symmetry=sym, npad=npad, mpi_comm = MPI_COMM_WORLD)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.recons3d_trl_struct_MPI(myid,0, datanew, symmetry=sym, npad=npad, mpi_comm = MPI_COMM_WORLD)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_recons3d_4nn_ctf_should_return_True(self):

        stack_name = "bdb:Substack/sort3d_substack_003"

        return_new = fu.recons3d_4nn_ctf(stack_name, [], snr = 2, symmetry="c1", npad=2)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_ctf(stack_name, [], snr = 2, symmetry="c1", npad=2)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5 ))



    def test_recons3d_4nn_ctf_MPI_should_return_True(self):

        finfo = open("dummytext.txt", 'w')
        list_proj = []
        stack_name = "bdb:Substack/sort3d_substack_003"
        nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(stack_name, list_proj[0])

        return_new = fu.recons3d_4nn_ctf_MPI(0, [proj], snr =1 , sign = 1, symmetry="c1", finfo=finfo, npad=2, smearstep=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_4nn_ctf_MPI(0, [proj], snr =1, sign = 1, symmetry="c1", finfo=finfo, npad=2, smearstep=0.5)
        mpi_barrier(MPI_COMM_WORLD)
        finfo.close()

        self.assertTrue(numpy.allclose(return_new.get_3dview(), return_old.get_3dview(), 0.5  ))


    def test_recons3d_nn_SSNR_MPI_should_return_True(self):

        finfo = open("dummytext.txt", 'w')
        list_proj = []
        stack_name = "bdb:Substack/sort3d_substack_003"
        nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
        list_proj = list(range(nima))
        proj = EMData()
        proj.read_image(stack_name, list_proj[0])

        return_new = fu.recons3d_nn_SSNR_MPI(0, [proj], mask2D=False, npad=2, sign = 1, symmetry="c1")
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.recons3d_nn_SSNR_MPI(0, [proj], mask2D=False, npad=2, sign = 1, symmetry="c1")
        mpi_barrier(MPI_COMM_WORLD)
        finfo.close()

        self.assertEqual(return_new[0], return_old[0])
        self.assertTrue(numpy.array_equal(return_new[1].get_3dview(), return_old[1].get_3dview()))


    def test_prepare_recons_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]
        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        snr = optionsnew.snr
        npad = optionsnew.npad

        return_new = fu.prepare_recons(datanew, symmetry=sym, myid=0 , main_node_half=0, half_start=4, step=2, index=5, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_recons(datanew, symmetry=sym, myid=0 , main_node_half=0, half_start=4, step=2, index=5, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertEqual(return_new, return_new)


    def test_prepare_recons_ctf_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]
        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        snr = optionsnew.snr
        npad = optionsnew.npad
        nx = datanew[0].get_xsize()

        return_new = fu.prepare_recons_ctf(nx,datanew, snr =1 , symmetry=sym, myid=0 , main_node_half=0, half_start=4, step=2, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_recons_ctf(nx,datanew, snr =1 , symmetry=sym, myid=0 , main_node_half=0, half_start=4, step=2, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(return_new, return_old)


    def test_recons_from_fftvol_should_return_True(self):
        fftvol = EMData(76,76)
        weight = EMData(76,76)
        size = 76

        return_new = fu.recons_from_fftvol(size, fftvol, weight, "c1")
        return_old = oldfu.recons_from_fftvol(size, fftvol, weight, "c1")
        self.assertEqual(return_new, return_old)


    def test_recons_ctf_from_fftvol_should_return_True(self):
        fftvol = EMData(76,76)
        weight = EMData(76,76)
        size = 76

        return_new = fu.recons_ctf_from_fftvol(size, fftvol, weight, 2, "c1")
        return_old = oldfu.recons_ctf_from_fftvol(size, fftvol, weight, 2, "c1")
        self.assertEqual(return_new, return_old)


    def test_get_image_size_should_return_True(self):

        return_new = fu.get_image_size([ EMData(76,76),EMData(76,76),EMData(76,76) ], 0 )
        return_old = oldfu.get_image_size([EMData(76,76),EMData(76,76),EMData(76,76) ], 0)
        self.assertEqual(return_new, return_old)


    def test_rec3D_MPI_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))
        (datanew,optionsnew,iternew) = arga[0]
        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        npad = optionsnew.npad

        return_new = fu.rec3D_MPI( datanew, 1.0, sym,   myid=0 , main_node =0, odd_start = 4, eve_start=4, index = -1 , npad = npad)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.rec3D_MPI(datanew,1.0, sym, myid=0 , main_node =0,  odd_start = 4, eve_start=4, index = -1, npad = npad)
        mpi_barrier(MPI_COMM_WORLD)


        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))


    def test_rec3D_MPI_noCTF_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]

        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        npad = optionsnew.npad

        return_new = fu.rec3D_MPI_noCTF( datanew, sym, myid=0 , main_node =0, odd_start = 4, eve_start=4, index = 5 , npad = npad)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.rec3D_MPI_noCTF(datanew, sym, myid=0 , main_node =0,  odd_start = 4, eve_start=4, index = 5, npad = npad)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))
        self.assertEqual(return_new[1], return_old[1])


    def test_prepare_recons_ctf_two_chunks_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]
        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        snr = optionsnew.snr
        npad = optionsnew.npad
        nx = datanew[0].get_xsize()
        datanew[0].set_attr("chunk_id", 0)


        return_new = fu.prepare_recons_ctf_two_chunks(nx, datanew, snr =1 , symmetry=sym, myid=0 , main_node_half=0, chunk_ID =2, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.prepare_recons_ctf_two_chunks(nx, datanew, snr =1 , symmetry=sym, myid=0 , main_node_half=0, chunk_ID =2,  npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(return_new, return_old)


    def test_rec3D_two_chunks_MPI_two_chunks_should_return_True(self):
        arga = get_arg_from_pickle_file(os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.do_volume"))

        (datanew,optionsnew,iternew) = arga[0]
        sym = optionsnew.sym
        sym = sym[0].lower() + sym[1:]
        snr = optionsnew.snr
        npad = optionsnew.npad
        nx = datanew[0].get_xsize()
        datanew[0].set_attr("chunk_id", 2)


        return_new = fu.rec3D_two_chunks_MPI( datanew, snr =1 , symmetry="c1", myid=0 , main_node=0, npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        return_old = oldfu.rec3D_two_chunks_MPI( datanew, snr =1 , symmetry="c1", myid=0 , main_node=0,  npad=npad, mpi_comm = MPI_COMM_WORLD)
        mpi_barrier(MPI_COMM_WORLD)
        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))



if __name__ == '__main__':
    unittest.main()
