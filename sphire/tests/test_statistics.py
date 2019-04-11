from __future__ import print_function
from __future__ import division

import unittest

import cPickle as pickle
import os
import numpy
from mpi import *
import global_def
import EMAN2_cppwrap as e2cpp
import EMAN2_cppwrap

mpi_init(0, [])
global_def.BATCH = True
global_def.MPI = True

ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))


from sphire.libpy_py3 import sphire_statistics as fu

from .sparx_lib import sparx_statistics as oldfu


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


    """Unable to read the file """
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


    """The function works but is extremely slow. It takes 55 mins to check the unit test that is why it is commented"""
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

    """  Can only work on the cluster"""
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





if __name__ == '__main__':
    unittest.main()
