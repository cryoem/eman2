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
from os import path,remove
from shutil import rmtree


#argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))
ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

def get_real_data(dim =2):
    """
    In order to run our tests with valid value this function returns the value got from pickle file.
    It returns the iamge and the reference image
    If dim=2 returns 2d img if 3 returns 3d image
    :param dim: dimension of the output image
    :return: image,refim #todo: in 3d case the refim is None ...i have to look for a better pickle file, maybe'projection.prgl'
    """
    if dim == 2:
        argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf"))
        image, refim, xrng, yrng = argum[0]
        return image,refim
    if dim == 3:
        argum = get_arg_from_pickle_file(path.join(ABSOLUTE_PATH, "pickle files/alignment.ali_vol_func"))
        return argum[1].get('data')[0], None

    print ("ERROR: the dimension has to be 2 or 3. Given "+str(dim))
    exit(-1)

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

def remove_list_of_file(l):
    """
    Remove the given directory with all its files
    :param l: list of file to remove
    """
    for f in l:
        if path.isfile(f):
            remove(f)

def returns_values_in_file(f):
    """
    read a file and returns all its lines
    :param f: path to file
    :return: contained values
    """
    if path.isfile(f):
        f1 = open(f, 'r')
        values_f1 = f1.readlines()
        f1.close()
        return values_f1
    print ("ERROR> the given file '"+str(f)+"' is not present!")
    exit(-1)

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


