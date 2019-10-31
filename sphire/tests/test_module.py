"""
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
"""


"""
helper function used in the tests
List:
1) create_kb --> create a Kaiser Bessel filter
2) get_data --> create a list of 2d EMdata image
3) get_data_3d --> create a list of 3d EMdata image
4) get_arg_from_pickle_file --> returns the values saved in a given pickle file
5) remove_dir --> removed a given directory with its files
6) remove_list_of_file --> removed a given lists of files
7) returns_values_in_file --> read the file and give back the text
"""

from os import system as os_system
from copy import deepcopy
from numpy import arange, float32 as np_float32
from sphire.libpy.sp_utilities import model_blank, model_circle
from EMAN2_cppwrap import Util, EMData
from cPickle import load as pickle_load
from os import path, remove


""" 
In order to run automatically all the tests download the precalculated results from the tutorial page http://sphire.mpg.de/wiki/doku.php?id=downloads:sphire_1_0.
or directly from http://sphire.mpg.de/wiki/doku.php?id=downloads:sphire_1_0
And set the variable 'ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER' to its path on your HD
"""

# ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER = "/home/adnan/Downloads/sphire_1_0_precalculated_results/SphireDemoResults"
ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER = "/home/lusnig/Downloads/SphireDemoResults"

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))


def get_real_data(dim=2):
    """
    In order to run our tests with valid value this function returns the value got from pickle file.
    It returns the iamge and the reference image
    If dim=2 returns 2d img if 3 returns 3d image
    :param dim: dimension of the output image
    :return: image,refim #todo: in 3d case the refim is None ...i have to look for a better pickle file, maybe'projection.prgl'
    """
    if dim == 2:
        argum = get_arg_from_pickle_file(
            path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf")
        )
        image, refim, xrng, yrng = argum[0]
        return image, refim
    if dim == 3:
        argum = get_arg_from_pickle_file(
            path.join(ABSOLUTE_PATH, "pickle files/alignment.ali_vol_func")
        )
        return argum[1].get("data")[0], None

    print("ERROR: the dimension has to be 2 or 3. Given " + str(dim))
    exit(-1)


def get_data(num, dim=10):
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
        data_a[...] = arange(dim * dim, dtype=np_float32).reshape(dim, dim) + i
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
        a = EMData(dim, dim, dim)
        data_a = a.get_3dview()
        data_a[...] = (
            arange(dim * dim * dim, dtype=np_float32).reshape(dim, dim, dim) + i
        )
        data_list.append(a)
    return data_list


def get_arg_from_pickle_file(filepath):
    """
    Returns the arguments of the given pickle file
    :param filepath: path of the pickle file
    :return: args saved in the file
    """
    with open(filepath, "rb") as rb:
        return pickle_load(rb)


def remove_dir(d):
    """
    I do not know why "shutil.rmtree" does not work
    Remove the given directory with all its files
    :param d: path to the directory
    """
    if path.isdir(d):
        os_system("rm -rf " + str(d))


def remove_list_of_file(l):
    """
    Remove the given directory with all its files
    :param l: list of file to remove
    """
    for f in l:
        if path.isfile(f):
            remove(f)


def returns_values_in_file(f, mode="r"):
    """
    read a file and returns all its lines
    :param f: path to file
    :param mode: how open the file. Default: read file
    :return: contained values
    """
    if path.isfile(f):
        f1 = open(f, mode)
        values_f1 = f1.readlines()
        f1.close()
        return values_f1
    print("ERROR> the given file '" + str(f) + "' is not present!")
    exit(-1)


def create_kb(dim, sizex=100, sizey=80, sizez=70):
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
    npad = 2
    if dim == 3:
        Nx = sizex * npad
        Ny = sizey * npad
        Nz = sizez * npad
        kbx = Util.KaiserBessel(alpha, K, sizex / 2, K / (2.0 * Nx), Nx)
        kby = Util.KaiserBessel(alpha, K, sizey / 2, K / (2.0 * Ny), Ny)
        kbz = Util.KaiserBessel(alpha, K, sizez / 2, K / (2.0 * Nz), Nz)
        return kbx, kby, kbz
    elif dim == 1:
        N = sizex * npad
        return Util.KaiserBessel(alpha, K, sizex / 2, K / (2.0 * N), N)
    print("Error: The value of dim has to be 1 or 3. You inserted " + str(dim))
    exit(-1)


""" In order to unittest the function which output is an EMData() we have to shrink the original images."""
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)
IMAGE_3D = deepcopy(IMAGE_2D)
IMAGE_2D.set_size(10, 10)
IMAGE_2D_REFERENCE.set_size(10, 10)
IMAGE_BLANK_2D = model_blank(nx=20, ny=20, bckg=0.0)

"""The resized 3DIMAGE is blank --> hence I fake the 2DImage to be a 3DIMAGE """

IMAGE_3D.set_size(10, 10, 10)
IMAGE_BLANK_3D = model_blank(nx=10, ny=10, nz=10, bckg=0.0)

MASK = model_circle(r=2, nx=5, ny=5, nz=1)
MASK_2DIMAGE = model_circle(r=2, nx=IMAGE_2D.get_xsize(), ny=IMAGE_2D.get_ysize(), nz=1)
MASK_IMAGE_BLANK_2D = model_circle(
    r=2, nx=IMAGE_BLANK_2D.get_xsize(), ny=IMAGE_BLANK_2D.get_ysize(), nz=1
)
MASK_3DIMAGE = model_circle(
    r=2, nx=IMAGE_3D.get_xsize(), ny=IMAGE_3D.get_ysize(), nz=IMAGE_3D.get_zsize()
)
MASK_IMAGE_BLANK_3D = model_circle(
    r=2,
    nx=IMAGE_BLANK_3D.get_xsize(),
    ny=IMAGE_BLANK_3D.get_ysize(),
    nz=IMAGE_BLANK_3D.get_zsize(),
)

KB_IMAGE2D_SIZE = create_kb(dim=1, sizex=IMAGE_2D.get_xsize())
KB_IMAGE3D_SIZE_CUBIC = create_kb(dim=1, sizex=IMAGE_3D.get_xsize())
