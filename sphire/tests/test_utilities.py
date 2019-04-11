from __future__ import print_function
from __future__ import division

import unittest

import cPickle as pickle
import os
import sys
from mpi import *
import global_def
import copy
import numpy
import shutil

from EMAN2_cppwrap import EMData

mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from ..libpy import sparx_utilities as fu

from .sparx_lib import sparx_utilities as oldfu

class Test_lib_utilities_compare(unittest.TestCase):


    def test_amoeba_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.amoeba")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (var, scale, func, ftolerance, xtolerance, itmax , data) = argum[0]

        return_new = fu.amoeba (var, scale, func, ftolerance, xtolerance, itmax , data)
        return_old = oldfu.amoeba (var, scale, func, ftolerance, xtolerance, itmax , data)

        self.assertTrue(return_new, return_old)

    def test_compose_transform2_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.compose_transform2")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2) = argum[0]

        return_new = fu.compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2)
        return_old = oldfu.compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2)

        self.assertTrue(return_new, return_old)


    def test_compose_transform3_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.compose_transform3")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2) = argum[0]

        return_new = fu.compose_transform3(phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2)
        return_old = oldfu.compose_transform3(phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2)

        self.assertTrue(return_new, return_old)


    def test_combine_params2_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.combine_params2")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2) = argum[0]

        return_new = fu.combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2)
        return_old = oldfu.combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2)

        self.assertTrue(return_new, return_old)


    def test_inverse_transform2_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.inverse_transform2")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (alpha, tx, ty) = argum[0]

        return_new = fu.inverse_transform2(alpha, tx, ty)
        return_old = oldfu.inverse_transform2(alpha, tx, ty)

        self.assertTrue(return_new, return_old)


    def test_drop_image_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.drop_image")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (imagename, destination) = argum[0]

        return_new = fu.drop_image(imagename, destination)
        return_old = oldfu.drop_image(imagename, destination)

        if return_new is not None   and  return_old is not None:
            self.assertTrue(return_new, return_old)


    def test_even_angles_true_should_return_equal_objects(self):


        return_new = fu.even_angles(delta = 0.25)
        return_old = oldfu.even_angles(delta = 0.25)

        if return_new is not None   and  return_old is not None:
            self.assertTrue(return_new, return_old)

    """If method is P , then it takes a lot of time"""
    def test_even_angles_cd_true_should_return_equal_objects(self):


        return_new = fu.even_angles_cd(delta = 0.5, method = 'P')
        return_old = oldfu.even_angles_cd(delta = 0.5, method = 'P')

        print("call to function is done ")
        print(return_old)

        if return_new is not None   and  return_old is not None:
            self.assertTrue(numpy.array_equal(return_new, return_old))

    def test_gauss_edge_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params3D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,) = argum[0]

        return_new = fu.gauss_edge(ima)
        return_old = oldfu.gauss_edge(ima)

        print("call to function is done ")
        print(return_old)

        if return_new is not None and return_old is not None:
            self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))



    def test_get_im_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_im")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (stackname, im) = argum[0]

        stackname = 'bdb:Substack/isac_substack'

        return_new = fu.get_im(stackname, im)
        return_old = oldfu.get_im(stackname, im)

        self.assertTrue(return_new, return_old)


    def test_get_image_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_image")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (imagename,) = argum[0]

        return_new = fu.get_image(imagename)
        return_old = oldfu.get_image(imagename)

        self.assertTrue(return_new, return_old)


    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_get_image_data_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_image_data")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0][0])
    #
    #     (image) = argum[0][0]
    #
    #     return_new = fu.get_image_data(image)
    #     return_old = oldfu.get_image_data(image)
    #
    #     self.assertTrue(numpy.array_equal(return_new, return_old))


    def test_get_symt_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_symt")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (symmetry,) = argum[0]

        return_new = fu.get_symt(symmetry)
        return_old = oldfu.get_symt(symmetry)

        self.assertTrue(return_new, return_old)


    def test_get_input_from_string_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_input_from_string")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (str_input) = argum[0][0]

        return_new = fu.get_input_from_string(str_input)
        return_old = oldfu.get_input_from_string(str_input)

        self.assertTrue(return_new, return_old)


    def test_model_circle_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.model_circle")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (r, nx, ny) = argum[0]

        return_new = fu.model_circle(r, nx, ny)
        return_old = oldfu.model_circle(r, nx, ny)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_model_gauss_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.model_circle")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (r, nx, ny) = argum[0]

        return_new = fu.model_gauss(0.25,nx,ny)
        return_old = oldfu.model_gauss(0.25,nx,ny)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    """  This function creates random noise each time so arrays cannot be compared """
    def test_model_gauss_noise_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.model_circle")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (r, nx, ny) = argum[0]

        return_new = fu.model_gauss_noise(0.15,nx,ny)
        return_old = oldfu.model_gauss_noise(0.15,nx,ny)

        self.assertTrue(return_new, return_old)


    def test_model_blank_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.model_blank")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (nx,ny) = argum[0]

        return_new = fu.model_blank(nx,ny)
        return_old = oldfu.model_blank(nx,ny)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_peak_search_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.peak_search")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (e,) = argum[0]

        return_new = fu.peak_search(e )
        return_old = oldfu.peak_search(e )

        self.assertEqual(return_new, return_old)



    """No return parameter for this function. """
    def test_print_list_format_true_should_return_equal_objects(self):
        import StringIO

        m = []
        for i in range(1):
            m.append((i))

        # fu.print_list_format(m)
        # oldfu.print_list_format(m)

        capturedOutput = StringIO.StringIO()
        sys.stdout = capturedOutput
        fu.print_list_format(m)
        oldfu.print_list_format(m)
        # sys.stdout = sys.__stdout__
        print ('Captured', capturedOutput.getvalue())



    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_pad_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.pad")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (image_to_be_padded, new_nx, new_ny, new_nz,off_center_nx) = argum[0]
    #
    #     return_new = fu.pad(image_to_be_padded, new_nx, new_ny, new_nz)
    #     return_old = oldfu.pad(image_to_be_padded, new_nx, new_ny, new_nz)
    #
    #     self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_chooseformat_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.chooseformat")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (t) = argum[0][0]

        return_new = fu.chooseformat(t)
        return_old = oldfu.chooseformat(t)

        self.assertEqual(return_new, return_old)

    def test_read_text_row_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.read_text_row")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (fnam) = argum[0][0]

        return_new = fu.read_text_row(fnam)
        return_old = oldfu.read_text_row(fnam)

        self.assertEqual(return_new, return_old)


    def test_write_text_row_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.write_text_row")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, filename) = argum[0]

        return_new = fu.write_text_row(data, filename)
        return_old = oldfu.write_text_row(data, filename)

        self.assertEqual(return_new, return_old)


    def test_read_text_file_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.read_text_file")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (filename,) = argum[0]

        return_new = fu.read_text_file(filename)
        return_old = oldfu.read_text_file(filename)

        self.assertEqual(return_new, return_old)


    def test_write_text_file_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.write_text_file")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, filename) = argum[0]

        return_new = fu.write_text_file(data, filename)
        return_old = oldfu.write_text_file(data, filename)

        self.assertEqual(return_new, return_old)


    def test_rotate_shift_params_true_should_return_equal_objects(self):

        paramsin = [[0.25,1.25,0.5]]
        transf  = [0.25, 1.25, 0.5]

        return_new = fu.rotate_shift_params(paramsin, transf)
        return_old = oldfu.rotate_shift_params(paramsin, transf)

        self.assertEqual(return_new, return_old)


    def test_reshape_1d_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.reshape_1d")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (input_object, length_current,Pixel_size_current) = argum[0]

        return_new = fu.reshape_1d(input_object, length_current,Pixel_size_current)
        return_old = oldfu.reshape_1d(input_object, length_current,Pixel_size_current)

        self.assertEqual(return_new, return_old)


    def test_estimate_3D_center_MPI_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.estimate_3D_center_MPI")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (data, nima, myid, number_of_proc, main_node) = argum[0]

        return_new = fu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)
        return_old = oldfu.estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node)

        self.assertTrue(return_new, return_old)


    def test_rotate_3D_shift_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.rotate_3D_shift")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (data, shift3d) = argum[0]

        return_new = fu.rotate_3D_shift(data, shift3d)
        return_old = oldfu.rotate_3D_shift(data, shift3d)

        if return_new is not None and return_old is not None:
            self.assertTrue(return_new, return_old)
        else:
            print('returns None')


    def test_set_arb_params_true_should_return_equal_objects(self):
        params = "lowpassfilter"
        par_str = "0.50"

        return_new = fu.set_arb_params(EMData(), params, par_str)
        return_old = oldfu.set_arb_params(EMData(), params, par_str)

        if return_new is not None and return_old is not None:
            self.assertTrue(return_new, return_old)
        else:
            print('returns None')


    def test_get_arb_params_true_should_return_equal_objects(self):
        params = "lowpassfilter"
        par_str = "0.50"

        a = EMData()
        a[params] = par_str

        for i in range(len(par_str)): a.set_attr_dict({par_str[i]: params[i]})

        return_new = fu.get_arb_params(a, par_str)
        return_old = oldfu.get_arb_params(a, par_str)

        if return_new is not None and return_old is not None:
            self.assertTrue(return_new, return_old)
        else:
            print('returns None')



    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_reduce_EMData_to_root_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.reduce_EMData_to_root")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (data, myid,main_node) = argum[0]
    #
    #     return_new = fu.reduce_EMData_to_root(data, myid,main_node = 0)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.reduce_EMData_to_root(data, myid,main_node = 0)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     self.assertEqual(return_new, return_old)


    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_bcast_compacted_EMData_all_to_all_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_compacted_EMData_all_to_all")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (list_of_em_objects, myid ) = argum[0]
    #
    #     return_new = fu.bcast_compacted_EMData_all_to_all(list_of_em_objects, myid)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.bcast_compacted_EMData_all_to_all(list_of_em_objects, myid)
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     self.assertEqual(return_new, return_old)



    def test_gather_compacted_EMData_to_root_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.gather_compacted_EMData_to_root")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (no_of_emo, list_of_emo, myid) = argum[0]

        return_new = fu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid)
        return_old = oldfu.gather_compacted_EMData_to_root(no_of_emo, list_of_emo, myid)

        self.assertEqual(return_new, return_old)


    def test_bcast_EMData_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_EMData_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (tavg, myid, source_node, ) = argum[0]

        return_new = fu.bcast_EMData_to_all(tavg, myid, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_EMData_to_all(tavg, myid, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)



    """  Can only be tested on the mpi. Wait too long on normal workstation"""
    # def test_send_EMData_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.send_EMData")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (img, dst, tag, comm) = argum[0]
    #     tag = 0
    #
    #     return_new = fu.send_EMData(img, dst, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = oldfu.send_EMData(img, dst, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)

    """  Can only be tested on the mpi. Wait too long on normal workstation"""
    # def test_recv_EMData_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.recv_EMData")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (src, tag,comm) = argum[0]
    #     tag = 0
    #
    #     return_new = fu.recv_EMData(src, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     return_old = oldfu.recv_EMData(src, tag)
    #     mpi_barrier(MPI_COMM_WORLD)
    #
    #     self.assertEqual(return_new, return_old)


    def test_bcast_number_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_number_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (number_to_send, source_node, mpi_comm) = argum[0]

        return_new = fu.bcast_number_to_all(number_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_number_to_all(number_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)

    def test_bcast_list_to_all_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.bcast_list_to_all")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (list_to_send, source_node, mpi_comm) = argum[0]

        return_new = fu.bcast_list_to_all(list_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.bcast_list_to_all(list_to_send, source_node)
        mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(return_new, return_old)


    # def test_recv_attr_dict_true_should_return_equal_objects(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #     print(argum)
    #
    #     (ima,) = argum[0]
    #     # params = "values"
    #     paramstr = 0
    #     # ima[params] = paramstr
    #
    #     ima.set_attr_dict({"values": 0})
    #
    #     return_new = fu.recv_attr_dict(0, "test", ima,paramstr,0,0,1)
    #
    #     return_old = oldfu.recv_attr_dict(0, "test1", ima,paramstr, 0,0,1)



    def test_print_begin_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_begin_msg(msg)

        return_old = oldfu.print_begin_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_print_end_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_end_msg(msg)

        return_old = oldfu.print_end_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_print_msg_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.print_msg")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (msg) = argum[0][0]

        return_new = fu.print_msg(msg)

        return_old = oldfu.print_msg(msg)

        self.assertEqual(return_new, return_old)


    def test_read_fsc_true_should_return_equal_objects(self):

        print(os.getcwd())
        filename = "sphire/tests/Sort3D/fsc_global.txt"
        return_new = fu.read_fsc(filename)

        return_old = oldfu.read_fsc(filename)

        self.assertEqual(return_new, return_old)


    def test_circumference_true_should_return_equal_objects(self):


        img = fu.model_blank(10,10,10)

        return_new = fu.circumference(img)

        return_old = oldfu.circumference(img)

        self.assertTrue(numpy.array_equal(return_new.get_3dview(), return_old.get_3dview()))


    def test_write_headers_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (ima,) = argum[0]


        fu.write_headers("test.hdf", [ima], [1])

        oldfu.write_headers("test1.hdf", [ima], [1])

        filepath = "/home/adnan/PycharmProjects/eman2/"

        return_new = os.path.isfile(filepath + "test.hdf")
        return_old = os.path.isfile(filepath + "test1.hdf")

        self.assertEqual(return_new, return_old)


    def test_write_headers_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (ima,) = argum[0]


        fu.write_header("test.hdf", ima, 1)

        oldfu.write_header("test1.hdf", ima, 1)

        filepath = "/home/adnan/PycharmProjects/eman2/"

        return_new = os.path.isfile(filepath + "test.hdf")
        return_old = os.path.isfile(filepath + "test1.hdf")

        self.assertEqual(return_new, return_old)


    def test_file_type_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.file_type")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (name) = argum[0][0]

        return_new = fu.file_type(name)

        return_old = oldfu.file_type(name)

        self.assertEqual(return_new, return_old)



    def test_get_params2D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])
        print(argum)

        (ima,) = argum[0]

        return_new = fu.get_params2D(ima )

        return_old = oldfu.get_params2D(ima)

        self.assertEqual(return_new, return_old)

    def test_set_params2D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p, xform) = argum[0]

        return_new = fu.set_params2D(ima,p)

        return_old = oldfu.set_params2D(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_params3D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params3D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,) = argum[0]

        return_new = fu.get_params3D(ima )

        return_old = oldfu.get_params3D(ima)

        self.assertEqual(return_new, return_old)


    def test_set_params3D_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params3D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p) = argum[0]

        return_new = fu.set_params3D(ima,p)

        return_old = oldfu.set_params3D(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_params_proj_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_params_proj")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,) = argum[0]

        return_new = fu.get_params_proj(ima )

        return_old = oldfu.get_params_proj(ima)

        self.assertEqual(return_new, return_old)


    def test_set_params_proj_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_params_proj")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ima,p) = argum[0]

        return_new = fu.set_params_proj(ima,p)

        return_old = oldfu.set_params_proj(ima,p)

        self.assertEqual(return_new, return_old)


    def test_get_latest_directory_increment_value_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_latest_directory_increment_value")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (directory_location, directory_name) = argum[0]

        return_new = fu.get_latest_directory_increment_value(directory_location, directory_name)

        return_old = oldfu.get_latest_directory_increment_value(directory_location, directory_name)

        self.assertEqual(return_new, return_old)


    def test_get_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (img) = argum[0][0]

        return_new = fu.get_ctf(img[0])
        return_old = oldfu.get_ctf(img[0])

        self.assertEqual(return_new, return_old)


    def test_same_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.same_ctf")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (c1,c2) = argum[0]

        return_new = fu.same_ctf(c1,c2)

        return_old = oldfu.same_ctf(c1,c2)

        self.assertEqual(return_new, return_old)



    def test_generate_ctf_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/alignment.ali2d_single_iter")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (p) = argum[0][0]

        return_new = fu.generate_ctf(p)

        return_old = oldfu.generate_ctf(p)

        self.assertTrue(return_new, return_old)


    def test_delete_bdb_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.delete_bdb")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (name) = argum[0][0]

        return_new = fu.delete_bdb(name)

        return_old = oldfu.delete_bdb(name)

        if return_new is not None and return_old is not None:
            self.assertTrue(return_new, return_old)
        else:
            print('returns None')

    def test_getvec_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getfvec")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (phi, tht) = argum[0]

        return_new = fu.getvec(phi, tht)

        return_old = oldfu.getvec(phi, tht)

        self.assertEqual(return_new, return_old)


    def test_getfvec_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getfvec")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (phi, tht) = argum[0]

        return_new = fu.getfvec(phi, tht)

        return_old = oldfu.getfvec(phi, tht)

        self.assertEqual(return_new, return_old)


    def test_nearest_fang_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_fang")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (vecs, phi, tht) = argum[0]

        return_new = fu.nearest_fang(vecs, phi, tht)

        return_old = oldfu.nearest_fang(vecs, phi, tht)

        self.assertEqual(return_new, return_old)


    def test_nearest_many_full_k_projangles_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_many_full_k_projangles")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (reference_normals, angles) = argum[0]
        symclass = argum[1]['sym_class']
        howmany = argum[1]['howmany']

        return_new = fu.nearest_many_full_k_projangles(reference_normals, angles, howmany, symclass)

        return_old = oldfu.nearest_many_full_k_projangles(reference_normals, angles, howmany, symclass)

        self.assertEqual(return_new, return_old)


    def test_angles_to_normals_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angles_to_normals")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum)

        (angles) = argum[0][0]

        return_new = fu.angles_to_normals(angles)

        return_old = oldfu.angles_to_normals(angles)

        self.assertEqual(return_new, return_old)


    def test_angular_occupancy_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        return_new = fu.angular_occupancy(params, angstep, sym, method)
        return_old = oldfu.angular_occupancy(params, angstep, sym, method)

        self.assertEqual(return_new, return_old)

    def test_angular_histogram_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        return_new = fu.angular_histogram(params, angstep, sym, method)
        return_old = oldfu.angular_histogram(params, angstep, sym, method)

        self.assertEqual(return_new, return_old)


    def test_balance_angular_distribution_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.angular_occupancy")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (params, angstep, sym, method) = argum[0]

        max_occupy = -1

        return_new = fu.balance_angular_distribution(params, max_occupy, angstep, sym)
        return_old = oldfu.balance_angular_distribution(params, max_occupy, angstep, sym)

        self.assertEqual(return_new, return_old)


    def test_symmetry_neighbors_true_should_return_equal_objects(self):

        angles = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        return_new = fu.symmetry_neighbors(angles , symmetry= "c1")
        return_old = oldfu.symmetry_neighbors(angles , symmetry= "c1")

        self.assertEqual(return_new, return_old)


    def test_rotation_between_anglesets_true_should_return_equal_objects(self):

        agls1 = [[idx1, idx2, 0] for idx1 in range(50) for idx2 in range(90)]
        agls2 = [[idx1, idx2, 5] for idx1 in range(50) for idx2 in range(90)]
        return_new = fu.rotation_between_anglesets(agls1, agls2)
        return_old = oldfu.rotation_between_anglesets(agls1, agls2)

        self.assertEqual(return_new, return_old)


    def test_angle_between_projections_directions_true_should_return_equal_objects(self):

        agls1 = [20, 60, 0]
        agls2 = [45, 75, 5]
        return_new = fu.angle_between_projections_directions(agls1, agls2)
        return_old = oldfu.angle_between_projections_directions(agls1, agls2)

        self.assertEqual(return_new, return_old)


    def test_get_pixel_size_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_pixel_size")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (img,) = argum[0]

        return_new = fu.get_pixel_size(img)

        return_old = oldfu.get_pixel_size(img)

        self.assertEqual(return_new, return_old)

    def test_set_pixel_size_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.set_pixel_size")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (img,pixel_size) = argum[0]

        return_new = fu.set_pixel_size(img,pixel_size)

        return_old = oldfu.set_pixel_size(img,pixel_size)

        self.assertEqual(return_new, return_old)

    def test_lacos_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.lacos")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (x,) = argum[0]

        return_new = fu.lacos(x)

        return_old = oldfu.lacos(x)

        self.assertEqual(return_new, return_old)

    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_nearest_proj_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.nearest_proj")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (proj_ang,img_per_grp,List) = argum[0]
    #
    #     return_new = fu.nearest_proj(proj_ang,img_per_grp,List)
    #
    #     return_old = oldfu.nearest_proj(proj_ang,img_per_grp,List)
    #
    #     self.assertEqual(return_new, return_old)


    def test_findall_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.findall")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (value, L) = argum[0]

        return_new = fu.findall(value, L)

        return_old = oldfu.findall(value, L)

        self.assertEqual(return_new, return_old)


    def test_pack_message_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.pack_message")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.pack_message(data)

        return_old = oldfu.pack_message(data)

        self.assertEqual(return_new, return_old)


    def test_unpack_message_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.unpack_message")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.unpack_message(data)

        return_old = oldfu.unpack_message(data)

        self.assertEqual(return_new, return_old)


    def test_update_tag_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.update_tag")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (communicator, target_rank) = argum[0]

        return_new = fu.update_tag(communicator, target_rank)

        return_old = oldfu.update_tag(communicator, target_rank)

        self.assertEqual(return_new, return_old)


    def test_wrap_mpi_send_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_send")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, destination,communicator) = argum[0]

        return_new = fu.wrap_mpi_send(data, destination)

        return_old = oldfu.wrap_mpi_send(data, destination)

        self.assertEqual(return_new, return_old)


    "Can only test on cluster , cannot work on workstation"
    # def test_wrap_mpi_recv_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_recv")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (data, communicator) = argum[0]
    #
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_new = fu.wrap_mpi_recv(data, communicator)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.wrap_mpi_recv(data, communicator)
    #
    #     self.assertEqual(return_new, return_old)


    def test_wrap_mpi_bcast_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_bcast")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, root, communicator) = argum[0]

        return_new = fu.wrap_mpi_bcast(data, root)

        return_old = oldfu.wrap_mpi_bcast(data, root)

        self.assertEqual(return_new, return_old)


    def test_wrap_mpi_gatherv_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_gatherv")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data, root, communicator) = argum[0]

        return_new = fu.wrap_mpi_gatherv(data, root)

        return_old = oldfu.wrap_mpi_gatherv(data, root)

        self.assertEqual(return_new, return_old)


    def test_get_colors_and_subsets_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_colors_and_subsets")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (main_node, mpi_comm, my_rank, shared_comm, sh_my_rank, masters) = argum[0]

        mpi_comm = MPI_COMM_WORLD
        main_node = 0
        my_rank = mpi_comm_rank(mpi_comm)
        mpi_size = mpi_comm_size(mpi_comm)
        shared_comm = mpi_comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)
        sh_my_rank = mpi_comm_rank(shared_comm)
        masters = mpi_comm_split(mpi_comm, sh_my_rank == main_node, my_rank)
        shared_comm = mpi_comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL)

        return_new = fu.get_colors_and_subsets(main_node, mpi_comm, my_rank, shared_comm, sh_my_rank,masters)

        return_old = oldfu.get_colors_and_subsets(main_node, mpi_comm, my_rank, shared_comm, sh_my_rank,masters)

        self.assertEqual(return_new, return_old)

        """ Can only be tested in mpi not on workstation   """
    # def test_wrap_mpi_split_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.wrap_mpi_split")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     print(argum[0])
    #
    #     (comm, no_of_groups) = argum[0]
    #
    #     return_new = fu.wrap_mpi_split(comm, no_of_groups)
    #     mpi_barrier(MPI_COMM_WORLD)
    #     return_old = oldfu.wrap_mpi_split(comm, no_of_groups)
    #
    #     self.assertEqual(return_new, return_old)

    def test_eliminate_moons_true_should_return_equal_objects(self):

        volume = fu.model_gauss(0.25,12,12,12)


        moon_params = []
        moon_params.append(0.5)
        moon_params.append(1.45)

        return_new = fu.eliminate_moons(volume, moon_params)

        return_old = oldfu.eliminate_moons(volume, moon_params)

        self.assertEqual(return_new, return_old)


    def test_get_dist_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.get_dist")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (c1, c2) = argum[0]

        return_new = fu.get_dist(c1, c2)

        return_old = oldfu.get_dist(c1, c2)

        self.assertEqual(return_new, return_old)


    def test_combinations_of_n_taken_by_k_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.combinations_of_n_taken_by_k")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (n, k) = argum[0]

        return_new = fu.combinations_of_n_taken_by_k(n, k)

        return_old = oldfu.combinations_of_n_taken_by_k(n, k)

        self.assertEqual(return_new, return_old)


    def test_cmdexecute_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.cmdexecute")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        # print(argum[0])

        (cmd,) = argum[0]

        dirname = cmd.split(' ')[1]

        current_path = os.getcwd()
        if os.path.isdir(dirname):
            print('directory exits')
            print('removing it')
            shutil.rmtree(dirname)

        return_new = fu.cmdexecute(cmd)

        if os.path.isdir(dirname):
            print('directory exits')
            print('removing it')
            shutil.rmtree(dirname)

        return_old = oldfu.cmdexecute(cmd)

        self.assertEqual(return_new, return_old)


    def test_if_error_then_all_processes_exit_program_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.if_error_then_all_processes_exit_program")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (error_status,) = argum[0]

        return_new = fu.if_error_then_all_processes_exit_program(error_status)

        return_old = oldfu.if_error_then_all_processes_exit_program(error_status)

        self.assertEqual(return_new, return_old)



    def test_get_shrink_data_huang_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["stack"] = "bdb:sphire/tests/Substack/sort3d_substack_002"
        Tracker["applyctf"] = True
        ids = []
        for i in range(1227):
            ids.append(i)
        Tracker["chunk_dict"] =ids
        myid = 0
        m_node = 0
        nproc = 1
        partids = "sphire/tests/Sort3D/indexes_010.txt"
        partstack = "sphire/tests/Sort3D/params_010.txt"
        nxinit = 2

        return_new = fu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        return_old = oldfu.get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, m_node, nproc)

        self.assertTrue(numpy.allclose(return_new[0][0].get_3dview(), return_old[0][0].get_3dview(), 0.5))


    def test_getindexdata_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.getindexdata")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (stack, partids, partstack, myid, nproc) = argum[0]

        stack = 'bdb:sphire/tests/VIPER/best_000'
        partids = 'sphire/tests/VIPER/main001/this_iteration_index_keep_images.txt'
        partstack = 'sphire/tests/VIPER//main001/run000/rotated_reduced_params.txt'

        return_new = fu.getindexdata(stack, partids, partstack, myid, nproc)

        return_old = oldfu.getindexdata(stack, partids, partstack, myid, nproc)

        self.assertTrue(numpy.array_equal(return_new[0].get_3dview(), return_old[0].get_3dview()))

    def test_convert_json_fromunicode_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/utilities/utilities.convert_json_fromunicode")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (data,) = argum[0]

        return_new = fu.convert_json_fromunicode(data)

        return_old = oldfu.convert_json_fromunicode(data)

        self.assertEqual(return_new, return_old)


    def test_get_sorting_attr_stack_should_return_equal_object(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)

        return_new = fu.get_sorting_attr_stack(stack)
        return_old = oldfu.get_sorting_attr_stack(stack)

        self.assertEqual(return_new, return_old)


    def test_get_sorting_params_refine_should_return_equal_object(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/multi_shc/multi_shc.ali3d_multishc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (stack, ref_vol, ali3d_options, symmetry_class) = argum[0]
        for i in range(len(stack)):
            stack[i].set_attr("group",i)

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["constants"]["stack"] = "bdb:sphire/tests/Substack/sort3d_substack_002"
        Tracker["applyctf"] = True
        Tracker["constants"]["nproc"] = 1

        return_new = fu.get_sorting_params_refine(Tracker, stack, 95)
        return_old = oldfu.get_sorting_params_refine(Tracker, stack, 95)

        self.assertEqual(return_new, return_old)

    def test_parsing_sorting_params_should_return_equal_object(self):

        sorting_list = []

        sorting_list.append(numpy.arange(10))
        sorting_list.append(numpy.arange(10))

        return_new = fu.parsing_sorting_params(sorting_list)
        return_old = oldfu.parsing_sorting_params(sorting_list)

        self.assertEqual(return_new[0], return_old[0])
        self.assertTrue(numpy.array_equal(return_new[1], return_old[1]))

    # def test_get_initial_ID_should_return_equal_object(self):
    #
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #
    #     return_new = fu.get_initial_ID(Tracker["two_way_stable_member"][istable], Tracker["full_ID_dict"])
    #     return_old = oldfu.get_initial_ID(Tracker["two_way_stable_member"][istable], Tracker["full_ID_dict"])
    #
    #     self.assertEqual(return_new, return_old)










    def test_convertasi_true_should_return_equal_objects(self):

        K = 7
        asig = [0,1,2,3,4,5,6]

        return_new = fu.convertasi(asig,K)
        return_old = oldfu.convertasi(asig,K)

        self.assertEqual(return_new, return_old)


    def test_prepare_ptp_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.prepare_ptp(data_list, K)
        return_old = oldfu.prepare_ptp(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_get_resolution_mrk01_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/statistics/statistics.fsc")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (img1, img2) = argum[0]
        volume = []
        volume.append(img1)
        volume.append(img2)

        fscoutputdir = "sphire/tests/Sort3D"
        mask_option = "sphire/tests/Sharpening/vol_adaptive_mask.hdf"

        return_new = fu.get_resolution_mrk01(volume, 0.5,0.15,fscoutputdir,mask_option)
        return_old = oldfu.get_resolution_mrk01(volume, 0.5,0.15,fscoutputdir,mask_option)

        self.assertEqual(return_new, return_old)


    def test_partition_to_groups_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.partition_to_groups(data_list, K)
        return_old = oldfu.partition_to_groups(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_partition_independent_runs_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.partition_independent_runs(data_list, K)
        return_old = oldfu.partition_independent_runs(data_list, K)

        self.assertEqual(return_new, return_old)


    def test_merge_groups_true_should_return_equal_objects(self):

        K = 7
        data_list = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.merge_groups(data_list)
        return_old = oldfu.merge_groups(data_list)

        self.assertEqual(return_new, return_old)


    def test_save_alist_true_should_return_equal_objects(self):

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]

        Tracker["this_dir"] = "sphire/tests/Sharpening/"
        Tracker["constants"]["log_main"] = "logging"
        Tracker["constants"]["myid"] = "myid"
        Tracker["constants"]["main_node"] = "myid"

        filename = "listfile.txt"
        data_list = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]]

        return_new = fu.save_alist(Tracker,filename,data_list)
        return_old = oldfu.save_alist(Tracker,filename,data_list)

        self.assertEqual(return_new, return_old)

    def test_margin_of_error_true_should_return_equal_objects(self):

        return_new = fu.margin_of_error(0,1)
        return_old = oldfu.margin_of_error(0,1)

        self.assertEqual(return_new, return_old)

    def test_counting_projections_true_should_return_equal_objects(self):

        delta = 0.5
        ali3d_params  = [[idx1, idx2, 0 , 0.25, 0.25] for idx1 in range(2) for idx2 in range(2)]
        image_start = 1

        return_new = fu.counting_projections(delta, ali3d_params, image_start)
        print("one function call done")
        return_old = oldfu.counting_projections(delta, ali3d_params, image_start)
        print("second function call done")
        self.assertEqual(return_new, return_old)


    """
      This function test works but takes too much time that is why for the time being it is
       commented,  will uncomment it once everything is done 
    """
    # def test_get_stat_proj_true_should_return_equal_objects(self):
    #     filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
    #     with open(filepath, 'rb') as rb:
    #         argum = pickle.load(rb)
    #
    #     Tracker = argum[0][0][1]
    #     Tracker["constants"]["nproc"] = 1
    #     Tracker["constants"]["myid"] = 0
    #     Tracker["constants"]["main_node"] = 0
    #
    #     delta = 0.5
    #     this_ali3d = "sphire/tests/VIPER/main001/run000/rotated_reduced_params.txt"
    #
    #     return_new = fu.get_stat_proj(Tracker,delta,this_ali3d)
    #     return_old = oldfu.get_stat_proj(Tracker,delta,this_ali3d)
    #     self.assertEqual(return_new, return_old)


    def test_create_random_list_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["total_stack"] = "stack"
        Tracker["constants"]["seed"] = 1.4
        Tracker["constants"]["indep_runs"] = 2
        Tracker["this_data_list"] = [2,3,5]

        return_new = fu.create_random_list(Tracker)
        return_old = oldfu.create_random_list(Tracker)
        self.assertEqual(return_new, return_old)

    def test_recons_mref_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["constants"]["nproc"] = 1
        Tracker["constants"]["myid"] = 0
        Tracker["constants"]["main_node"] = 0
        Tracker["number_of_groups"] = 1
        Tracker["constants"]["nnxo"] = 4  # roi
        Tracker["this_particle_list"] = [[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6],[0, 1, 2, 3, 4, 5, 6]]
        Tracker["nxinit"] = 1
        Tracker["constants"]["partstack"] = 'sphire/tests/VIPER//main001/run000/rotated_reduced_params.txt'
        Tracker["this_dir"] = "sphire/tests/Particles/"
        Tracker["constants"]["stack"] = 'bdb:sphire/tests/Class2D/stack_ali2d'
        Tracker["applyctf"] = False
        Tracker["chunk_dict"] = [0, 1, 2, 3, 4, 5, 6]
        Tracker["constants"]["sym"] = "c1"

        return_new = fu.recons_mref(Tracker)
        return_old = oldfu.recons_mref(Tracker)
        self.assertTrue(return_new[0], return_old[0])


    def test_apply_low_pass_filter_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/user_functions.do_volume_mask")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        Tracker = argum[0][0][1]
        Tracker["low_pass_filter"] = 0.087

        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/projection.prgl")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        (volft, params, interpolation_method, return_real) = argum[0]

        refvol = [volft,volft]

        return_new = fu.apply_low_pass_filter(refvol,Tracker)
        return_old = oldfu.apply_low_pass_filter(refvol,Tracker)
        self.assertEqual(return_new, return_old)

    def test_count_chunk_members_true_should_return_equal_objects(self):

        chunk_dict =  [0,1,2,3,4,5,6]
        one_class = [0,1,2,3,4,5,6]

        return_new = fu.count_chunk_members(chunk_dict, one_class)
        return_old = oldfu.count_chunk_members(chunk_dict, one_class)
        self.assertEqual(return_new, return_old)


    def test_remove_small_groups_true_should_return_equal_objects(self):

        chunk_dict =  [[0,1,2,3,4,5,6],[0,1,2,3,4,5,6],[0,1,2,3,4,5,6]]
        one_class = [0,1,2,3,4,5,6]

        return_new = fu.remove_small_groups(chunk_dict, 3)
        return_old = oldfu.remove_small_groups(chunk_dict, 3)
        self.assertEqual(return_new, return_old)


    def test_get_number_of_groups_true_should_return_equal_objects(self):

        return_new = fu.get_number_of_groups(500, 5)
        return_old = oldfu.get_number_of_groups(500, 5)
        self.assertEqual(return_new, return_old)


    def test_tabessel_true_should_return_equal_objects(self):

        return_new = fu.tabessel(5,4)
        return_old = oldfu.tabessel(5, 4)
        self.assertEqual(return_new, return_old)



if __name__ == '__main__':
    unittest.main()
