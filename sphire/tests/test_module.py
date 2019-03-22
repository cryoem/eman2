"""
helper function used in the tests
List:
1) create_kb --> create a Kaiser Bessel filter
2) get_data --> create a list of 2d EMdata image
3) get_data_3d --> create a list of 3d EMdata image
4) get_arg_from_pickle_file --> returns the values saved in a given pickle file
5) remove_dir --> removed a given directory with its files
"""

import numpy

from EMAN2_cppwrap import Util,EMData
from cPickle import load as pickle_load
from os import path
from shutil import rmtree



def get_data(num,dim = 10):
    """
    Create a list of 2D EMData image
    :param num: number of the images to create
    :param dim: dimension of the image
    :return: list of EMData()
    """
    data_list = []
    for i in range(num):
        a = EMData(dim, dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim, dtype=numpy.float32).reshape(dim, dim) + i
        data_list.append(a)
    return data_list


def get_data_3d(num, dim=10):
    """
    Create a list of 3D EMData image
    :param num: number of the images to create
    :param dim: dimension of the image
    :return: list of EMData()
    """
    data_list = []
    for i in range(num):
        a = EMData(dim, dim,dim)
        data_a = a.get_3dview()
        data_a[...] = numpy.arange(dim * dim * dim, dtype=numpy.float32).reshape(dim, dim, dim) + i
        data_list.append(a)
    return data_list


def get_arg_from_pickle_file(filepath):
    """
    Returns the arguments of the given pickle file
    :param filepath: path of the pickle file
    :return: args saved in the file
    """
    with open(filepath, 'rb') as rb:
        return pickle_load(rb)

def remove_dir(d):
    """
    Remove the given directory with all its files
    :param d: path to the directory
    """
    if path.isdir(d):
        rmtree(d)

def create_kb(dim,  sizex=100 ,sizey=80 ,sizez=70):
    """
    Return a kb filter. If you need more info about the code, please take a look of 'prep_vol' in sparx_projection.py
    if dim =3 returns kbx,kby, kbz the interpolants along x, y and z direction (tabulated Kaiser-Bessel function) because the volume is rectangular
    if dim =1 returns kb the  interpolants (tabulated Kaiser-Bessel function) because the volume is cubic.
    :param dim: dimension of the kb filter
    :param sizex: in case of dim =3 is the xsize. All the axis in case of dim =1
    :param sizey: in case of dim =3 is the xsize
    :param sizez: in case of dim =3 is the xsize
    :return: kb or kbx,kby, kbz
    """
    K = 6
    alpha = 1.75
    npad =2
    if dim == 3:
        Nx = sizex * npad
        Ny = sizey * npad
        Nz = sizez * npad
        kbx = Util.KaiserBessel(alpha, K, sizex / 2, K / (2. * Nx), Nx)
        kby = Util.KaiserBessel(alpha, K, sizey / 2, K / (2. * Ny), Ny)
        kbz = Util.KaiserBessel(alpha, K, sizez / 2, K / (2. * Nz), Nz)
        return kbx, kby, kbz
    elif dim == 1:
        N = sizex * npad
        return Util.KaiserBessel(alpha, K, sizex / 2, K / (2. * N), N)
    print ("Error: The value of dim has to be 1 or 3. You inserted " + str(dim))
    exit(-1)


