# encoding: utf-8
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
from pickle import load as pickle_load
from os import path, remove

# import io
# with io.open("my_utf8_file.txt", "r", encoding="utf-8") as my_file:
#      my_unicode_string = my_file.read()

""" 
In order to run automatically all the tests download the precalculated results from the tutorial page 
http://sphire.mpg.de/wiki/doku.php?id=downloads:sphire_1_0.
or directly from http://sphire.mpg.de/wiki/doku.php?id=downloads:sphire_1_0
And set the variable 'ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER' to its path on your HD
"""

# ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER = "/home/adnan/Downloads/sphire_1_0_precalculated_results/SphireDemoResults"
ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER = "/home/adnan/DemoResults/"

'''
the sphire1.3 tutorial do not use the values downloadable from the SPHIRE website. You have to download them from billy: 
sphire-devel/SPHIRE_DEMO_RESULTS/SphireDemoResults2 
'''
ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER_NEW = "/home/adnan/DemoResults/"

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

#Absolute path to the bin folders for compatibility tests purpose.  since in the bin folder we have to set the parameter via cmd I have to run them.
ABSOLUTE_SPHIRE_PATH= ABSOLUTE_PATH.split("/tests")[0]
ABSOLUTE_OLDBIN_PATH = path.join(ABSOLUTE_SPHIRE_PATH, "bin_py3")
ABSOLUTE_BIN_PATH = path.join(ABSOLUTE_SPHIRE_PATH,"bin")

#absolute of python
ABSOLUTE_PATH_PYTHON = "/home/adnan/applications/sphire/miniconda3/envs/py3_v5/bin/python"


def get_real_data(dim=2):
    """
    In order to run our tests with valid value this function returns the value got from pickle file.
    It returns the iamge and the reference image
    If dim=2 returns 2d img if 3 returns 3d image
    :param dim: dimension of the output image
    :return: image,refim #todo: in 3d case the refim is None ...i have to look for a better pickle file, maybe'projection.prgl'
    """
    from EMAN2 import EMNumPy
    if dim == 2:
        # argum = get_arg_from_pickle_file(
        #     path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf_new")
        # )
        # image, refim, xrng, yrng = argum[0]
        argum = get_arg_from_pickle_file(
            path.join(ABSOLUTE_PATH, "pickle files/alignment.align2d_scf_new")
        )

        arrimage,dictimage, arrrefim, dictrefim, xrng, yrng = argum


        image = EMNumPy.numpy2em(arrimage)
        image.set_attr_dict(dictimage)

        refim = EMNumPy.numpy2em(arrrefim)
        refim.set_attr_dict(dictrefim)
        return image, refim

    if dim == 3:
        argum = get_arg_from_pickle_file(
            path.join(ABSOLUTE_PATH, "pickle files/alignment.ali_vol_func_new")
        )

        imagearr = EMNumPy.numpy2em(argum[0])
        imagearr.set_attr_dict(argum[1])
        return imagearr, None
    # if dim == 3:
    #     argum = get_arg_from_pickle_file(
    #         path.join(ABSOLUTE_PATH, "pickle files/alignment.ali_vol_func")
    #     )
    #     return argum[1].get("data")[0], None

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
    import io
    import pickle

    with io.open(filepath, "rb") as rb:
        # u = pickle._Unpickler(rb)
        # try:
        #     u.encoding = 'latin1'
        # except:
        #     u.encoding = 'utf-8'
        # return u.load()
        return pickle.load(rb, encoding='latin1')


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


def returns_values_in_file(f, mode="rb"):
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

import copy
""" In order to unittest the function which output is an EMData() we have to shrink the original images."""
IMAGE_2D, IMAGE_2D_REFERENCE = get_real_data(dim=2)

img2d_data = IMAGE_2D.get_2dview()
img2d_dict = IMAGE_2D.get_attr_dict()
from EMAN2 import EMNumPy

IMAGE_3D =EMNumPy.numpy2em(copy.deepcopy(img2d_data))
IMAGE_3D.set_attr_dict(img2d_dict)
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




import io
import pickle
import sp_utilities
import os
import dill


def give_ali2d_single_iter_data():
    numr = [1, 1, 8, 2, 9, 16, 3, 25, 32, 4, 57, 32, 5, 89, 32, 6, 121, 64, 7, 185, 64, 8, 249, 64, 9, 313, 64, 10, 377,
            64, 11, 441, 128, 12, 569, 128, 13, 697, 128, 14, 825, 128, 15, 953, 128, 16, 1081, 128, 17, 1209, 128, 18,
            1337, 128, 19, 1465, 128, 20, 1593, 128, 21, 1721, 256, 22, 1977, 256, 23, 2233, 256, 24, 2489, 256, 25,
            2745,
            256, 26, 3001, 256, 27, 3257, 256, 28, 3513, 256, 29, 3769, 256]
    wr = [25.132741228718345, 12.566370614359172, 4.71238898038469, 6.283185307179586, 7.853981633974483,
          2.356194490192345, 2.748893571891069, 3.141592653589793, 3.5342917352885173, 3.9269908169872414,
          1.0799224746714913, 1.1780972450961724, 1.2762720155208536, 1.3744467859455345, 1.4726215563702154,
          1.5707963267948966, 1.6689710972195777, 1.7671458676442586, 1.8653206380689396, 1.9634954084936207,
          0.5154175447295755, 0.5399612373357456, 0.5645049299419159, 0.5890486225480862, 0.6135923151542565,
          0.6381360077604268, 0.662679700366597, 0.6872233929727672, 0.7117670855789375]
    cs = [0.0, 0.0]
    cnx = 36
    cny = 36
    xrng = 4
    yrng = 4
    step = 1

    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/alignment.ali2d_single_iter"

    data_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_data_values.pkl"
    data_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_data_head.pkl"
    tavg_value_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_tavg_values.pkl"
    tavg_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_tavg_head.pkl"

    dill._dill._reverse_typemap["ObjectType"] = object

    # print(data_values_name)
    # print(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_values_name))

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_values_name), "rb") as f1:
        data_arras = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_head_name), "rb") as f2:
        data_heads = pickle.load(f2, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + tavg_value_name), "rb") as f3:
        tavg_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" +tavg_head_name), "rb") as f4:
        tavg_head = pickle.load(f4, encoding='latin1')

    data = []
    for i in range(len(data_arras)):
        data.append(sp_utilities.numpy2em_python(data_arras[i]))
        for key in data_heads[i]:
            data[i][key] = data_heads[i][key]

    tavg = sp_utilities.numpy2em_python(tavg_data)
    for keys in tavg_head:
        tavg[keys] = tavg_head[keys]

    return data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step


def give_ornq_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/alignment.ornq"

    image_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_image_values.pkl"
    image_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_image_head.pkl"

    crefim_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_crefim_values.pkl"
    crefim_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "_crefim_head.pkl"


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + image_values_name), "rb") as f3:
        image_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + image_head_name), "rb") as f4:
        image_head = pickle.load(f4, encoding='latin1')


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + crefim_values_name), "rb") as f1:
        crefim_data = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + crefim_head_name), "rb") as f2:
        crefim_head = pickle.load(f2, encoding='latin1')



    image = sp_utilities.numpy2em_python(image_data)
    for keys in image_head:
        image[keys] = image_head[keys]

    crefim = sp_utilities.numpy2em_python(crefim_data)
    for keys in crefim_head:
        crefim[keys] = crefim_head[keys]

    return image , crefim


def give_ormq_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/alignment.ormq"

    image_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "ormq_image_values.pkl"
    image_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "ormq_image_head.pkl"

    crefim_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "ormq_crefim_values.pkl"
    crefim_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "ormq_crefim_head.pkl"


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + image_values_name), "rb") as f3:
        image_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + image_head_name), "rb") as f4:
        image_head = pickle.load(f4, encoding='latin1')


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + crefim_values_name), "rb") as f1:
        crefim_data = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + crefim_head_name), "rb") as f2:
        crefim_head = pickle.load(f2, encoding='latin1')



    image = sp_utilities.numpy2em_python(image_data)
    for keys in image_head:
        image[keys] = image_head[keys]

    crefim = sp_utilities.numpy2em_python(crefim_data)
    for keys in crefim_head:
        crefim[keys] = crefim_head[keys]

    return image , crefim



def give_ali_vol_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/applications.ali_vol"

    vol_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "vol_values.pkl"
    vol_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "vol_head.pkl"

    refv_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "refv_values.pkl"
    refv_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "refv_head.pkl"


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + vol_values_name), "rb") as f3:
        vol_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + vol_head_name), "rb") as f4:
        vol_head = pickle.load(f4, encoding='latin1')


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + refv_values_name), "rb") as f1:
        refv_data = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + refv_head_name), "rb") as f2:
        refv_head = pickle.load(f2, encoding='latin1')



    vol = sp_utilities.numpy2em_python(vol_data)
    for keys in vol_head:
        vol[keys] = vol_head[keys]

    refv = sp_utilities.numpy2em_python(refv_data)
    for keys in refv_head:
        refv[keys] = refv_head[keys]

    return vol , refv


def give_ali2d_base_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/applications.ali2d_base"
    stack_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "stack_values.pkl"
    stack_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "stack_head.pkl"
    dill._dill._reverse_typemap["ObjectType"] = object

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + stack_values_name), "rb") as f1:
        data_arras = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + stack_head_name), "rb") as f2:
        data_heads = pickle.load(f2, encoding='latin1')

    stack = []
    for i in range(len(data_arras)):
        stack.append(sp_utilities.numpy2em_python(data_arras[i]))
        for key in data_heads[i]:
            stack[i][key] = data_heads[i][key]

    return stack


def give_rotate_3D_shift_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/utilities/utilities.rotate_3D_shift"

    data_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "data_values.pkl"
    data_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "data_head.pkl"
    dill._dill._reverse_typemap["ObjectType"] = object

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_values_name), "rb") as f1:
        data_arras = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_head_name), "rb") as f2:
        data_heads = pickle.load(f2, encoding='latin1')

    data = []
    for i in range(len(data_arras)):
        data.append(sp_utilities.numpy2em_python(data_arras[i]))
        for key in data_heads[i]:
            data[i][key] = data_heads[i][key]

    return data


def give_projection_prgl_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/projection.prgl"

    volft_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "volft_values.pkl"
    volft_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "volft_head.pkl"


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + volft_values_name), "rb") as f3:
        volft_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + volft_head_name), "rb") as f4:
        volft_head = pickle.load(f4, encoding='latin1')


    volft = sp_utilities.numpy2em_python(volft_data)
    for keys in volft_head:
        volft[keys] = volft_head[keys]


    return volft


def give_alignment_shc_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/alignment.shc"

    xfomimg_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "xfomimg_values.pkl"
    xfomimg_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "xfomimg_head.pkl"


    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + xfomimg_values_name), "rb") as f3:
        xfomimg_data = pickle.load(f3, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + xfomimg_head_name), "rb") as f4:
        xfomimg_head = pickle.load(f4, encoding='latin1')


    xformimg = sp_utilities.numpy2em_python(xfomimg_data)
    for keys in xfomimg_head:
        xformimg[keys] = xfomimg_head[keys]


    return xformimg



def give_ave_series_data():
    old_pkl = "/home/adnan/PycharmProjects/python3conversion/sphire/tests/pickle files/statistics/statistics.ave_series"
    data_ave_values_name = os.path.splitext(os.path.basename(old_pkl))[0] + "data_ave_values.pkl"
    data_ave_head_name = os.path.splitext(os.path.basename(old_pkl))[0] + "data_ave_head.pkl"
    dill._dill._reverse_typemap["ObjectType"] = object

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_ave_values_name), "rb") as f1:
        data_arras = pickle.load(f1, encoding='latin1')

    with io.open(os.path.join(ABSOLUTE_PATH, "pickle files/" + data_ave_head_name), "rb") as f2:
        data_heads = pickle.load(f2, encoding='latin1')

    stack = []
    for i in range(len(data_arras)):
        stack.append(sp_utilities.numpy2em_python(data_arras[i]))
        for key in data_heads[i]:
            stack[i][key] = data_heads[i][key]

    return stack