#! /usr/bin/env python

# Author: Adnan Ali, 17/06/2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Markus Stabrin 17/06/2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 17/06/2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Copyright (c) 2019 MPI Dortmund

# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
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

from __future__ import print_function
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy as sp
import pandas as pd
import copy
import mrcfile


from EMAN2 import *
import EMAN2_cppwrap
from EMAN2 import EMNumPy

import sp_utilities
import sp_projection
import sp_statistics
import sp_filter
import mpi
import sp_applications
import argparse
import sp_fundamentals
# from progressbar import *

import scipy
from scipy.optimize import minimize, Bounds, fmin_l_bfgs_b

import sp_interpolation
from tqdm import tqdm

location =os.getcwd()
RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
main_mpi_proc = 0
if RUNNING_UNDER_MPI:
    pass  # IMPORTIMPORTIMPORT from mpi import mpi_init
    pass  # IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM

    mpi.mpi_init(0, [])
    my_mpi_proc_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    n_mpi_procs = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    shared_comm = mpi.mpi_comm_split_type(mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED, 0, mpi.MPI_INFO_NULL)
    my_node_proc_id = mpi.mpi_comm_rank(shared_comm)
else:
    my_mpi_proc_id = 0
    n_mpi_procs = 1
    my_node_proc_id = 0

size = n_mpi_procs
rank  = my_node_proc_id


no_of_micrographs = 24
N_ITER =25
shift = 1

try:
    ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))
    GLOBAL_PATH = os.path.abspath(os.path.join(__file__ ,"../../../"))
except:
    GLOBAL_PATH = os.getcwd()


def readimg(path):
    with mrcfile.open(path) as mrc:
        dat = mrc.data
    return dat


def numpy2em_python(numpy_array, out=None):
    """
	Create an EMData object based on a numpy array by reference.
	The output EMData object will have the reversed order of dimensions.
	x,y,z -> z,y,x
	Arguments:
	numpy_array - Array to convert
	Return:
	EMData object
	"""
    if out is None:
        shape = numpy_array.shape[::-1]
        if len(shape) == 1:
            shape = (shape[0], 1)
        return_array = EMData(*shape)
    else:
        return_array = out
    return_view = EMNumPy.em2numpy(return_array)
    return_view[...] = numpy_array
    return_array.update()
    if out is None:
        return return_array


def return_movie_names(input_image_path):
    mic_pattern = input_image_path
    mic_basename_pattern = os.path.basename(mic_pattern)
    global_entry_dict = {}
    subkey_input_mic_path = "Input Micrograph Path"

    mic_basename_tokens = mic_basename_pattern.split('*')
    mic_id_substr_head_idx = len(mic_basename_tokens[0])

    input_mic_path_list = glob.glob(mic_pattern)

    for input_mic_path in input_mic_path_list:
        # Find tail index of  id substring and extract the substring from the  name
        input_mic_basename = os.path.basename(input_mic_path)
        mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
        mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
        if not mic_id_substr in global_entry_dict:
            global_entry_dict[mic_id_substr] = {}
        global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path

    print(" ")
    print("\n Summary of dataset consistency check...")
    print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
    print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))

    valid_mic_id_substr_list = []
    for mic_id_substr in global_entry_dict:
        mic_id_entry = global_entry_dict[mic_id_substr]
        valid_mic_id_substr_list.append(mic_id_substr)

    input_file_path_list = []
    for mic_id_substr in sorted(valid_mic_id_substr_list):
        mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
        input_file_path_list.append(mic_path)

    del global_entry_dict
    del valid_mic_id_substr_list

    print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
    return input_file_path_list


"""
Reads the x and y shifts values per frame in a micrograph  
"""
def returns_values_in_file(f, mode = 'r'):
    """
    read a file and returns all its lines
    :param f: path to file
    :param mode: how open the file. Default: read file
    :return: contained values
    """
    if os.path.isfile(f):
        with open(f, mode) as f1:
            values_f1 = f1.readlines()
        return values_f1
    print ("ERROR> the given file '"+str(f)+"' is not present!")
    exit(-1)


def read_meta_shifts(f):
    x = []
    y  = []
    for row in returns_values_in_file(f):
        if "image #" in row:
            v=row.replace('\n','').split('=')[1].replace(' ', '').split(',')
            x.append(float(v[0]))
            y.append(float(v[1]))
        elif "Frame (" in row:
            v = row.split()
            x.append(float(v[-2]))
            y.append(float(v[-1]))
    return x,y

def giveshiftperparticles(trackingfile , frames):
    df = pd.read_csv(trackingfile)
    particle_number = int(df['data_general'].values[0].rsplit()[1])
    shiftingperptcl = np.zeros((particle_number, frames, 2))

    for i in range(particle_number):
        for j in range(frames):
            if i == 0:
                k = 5 + j
            if i > 0:
                k = i * frames + ((4 * i) + j) + 5
            shiftingperptcl[i, j, 0] = float(df['data_general'].values[k].rsplit()[0])
            shiftingperptcl[i, j, 1] = float(df['data_general'].values[k].rsplit()[1])
    return shiftingperptcl



def givemotioncorrshifts(filename):
    shift_file = pd.read_csv(filename)
    # header = shift_file['data_general'][0:16]
    shift_x = []
    shift_y = []
    for i in range(17, len(shift_file['data_general'])):
        shift_x.append(float(shift_file['data_general'][i].split()[1]))
        shift_y.append(float(shift_file['data_general'][i].split()[2]))

    return shift_x, shift_y

def applyshifts(image, shifx , shify):
    newimage = sp_fundamentals.rot_shift2D(image, sx=shifx, sy=shify, interpolation_method='quadratic')
    return newimage

def read_all_attributes_from_stack(stack):
    no_of_imgs_once = EMUtil.get_image_count(stack)  # Counting how many are there in the stack
    # -------Extracting the information from the substack
    ptcl_source_images_once = EMUtil.get_all_attributes(stack, 'ptcl_source_image')
    project_params_all_once = EMUtil.get_all_attributes(stack, "xform.projection")
    particle_coordinates_all_once = EMUtil.get_all_attributes(stack, "ptcl_source_coord")
    ctf_params_all_once = EMUtil.get_all_attributes(stack, "ctf")
    nx_all_once = EMUtil.get_all_attributes(stack, 'nx')
    ny_all_once = EMUtil.get_all_attributes(stack, 'ny')
    nz_all_once = EMUtil.get_all_attributes(stack, 'nz')
    chunk_params_all_once = EMUtil.get_all_attributes(stack, "chunk_id")

    return no_of_imgs_once, ptcl_source_images_once, project_params_all_once, \
           particle_coordinates_all_once, ctf_params_all_once, nx_all_once, ny_all_once, nz_all_once,chunk_params_all_once


def find_particles_info_from_movie(stack, movie_name, no_of_imgs, ptcl_source_images, project_params_all,
                                   particle_coordinates_all, ctf_params_all, nx_all, ny_all, nz_all,
                                    adnan_n_all, chunk_params_all, show_first = False):

    #-----------------------------------------------   CTF related attributes
    """
    defocus = defocus associated with the image, positive value corresponds to underfocus
    cs =  spherical aberration constant [mm].
    voltage = accelerating voltage of the microscope [kV]
    apix = angular pixel size
    bfactor = The parameter in Gaussian like envelope function, which roughly explains
                Fourier factor dumping of the image.
    ampcont = amplitude contrast
    dfdiff  = astigmatism amplitude
    dfang =  astigmatism angle
    """

    # -------------------------------------------  2D orientation / orientation attributes
    """
    phi =  Eulerian angle for 3D reconstruction (azimuthal) 
    theta = Eulerian angle for 3D reconstruction (tilt) 
    psi = Eulerian angle for 3D reconstruction (in-plane rotation of projection) 
    tx =  shift in x direction
    ty = shift in y direction
    """
    print("Number of images in the substack are %d" % len(ptcl_source_images))

    project_params_per_movie = []
    particle_coordinates_per_movie = []
    ctf_params_per_movie = []
    nx_per_movie = []
    ny_per_movie = []
    nz_per_movie = []
    adnan_per_movie = []
    chunk_per_movie = []

    if str(os.path.basename(movie_name)).split('.')[-1] == 'mrcs':
        for i in range(no_of_imgs):
            if (
                    str(os.path.basename(movie_name)) ==
                    str(os.path.basename(ptcl_source_images[i])) or
                    str(os.path.basename(movie_name)) ==
                    str(os.path.basename(ptcl_source_images[i]))+'s'
                    ):
                project_params_per_movie.append(project_params_all[i])
                particle_coordinates_per_movie.append(particle_coordinates_all[i])
                ctf_params_per_movie.append(ctf_params_all[i])
                nx_per_movie.append(nx_all[i])
                ny_per_movie.append(ny_all[i])
                nz_per_movie.append(nz_all[i])
                adnan_per_movie.append(adnan_n_all[i])
                chunk_per_movie.append(chunk_params_all[i])

    elif str(os.path.basename(movie_name)).split('.')[-1] == 'tiff':
        for i in range(no_of_imgs):
            if (
                    str(os.path.basename(movie_name)).split('.tiff')[0] ==
                    str(os.path.basename(ptcl_source_images[i])).split('.mrc')[0] or
                    str(os.path.basename(movie_name)).split('.tiff')[0] ==
                    str(os.path.basename(ptcl_source_images[i])).split('.mrcs')[0]
            ):
                project_params_per_movie.append(project_params_all[i])
                particle_coordinates_per_movie.append(particle_coordinates_all[i])
                ctf_params_per_movie.append(ctf_params_all[i])
                nx_per_movie.append(nx_all[i])
                ny_per_movie.append(ny_all[i])
                nz_per_movie.append(nz_all[i])
                adnan_per_movie.append(adnan_n_all[i])
                chunk_per_movie.append(chunk_params_all[i])
    if np.array(ctf_params_per_movie).size == 0:
        return False

    if str(os.path.basename(movie_name)).split('.')[-1] == 'mrcs':
        print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)),
                                                             len(project_params_per_movie)))
        print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
        print("Projection parameters for 1st particle in the stack are ", project_params_per_movie[0].get_params('spider'))

    elif str(os.path.basename(movie_name)).split('.')[-1] == 'tiff':
        print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)),
                                                             len(project_params_per_movie)))
        print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
        print("Projection parameters for 1st particle in the stack are ", project_params_per_movie[0].get_params('spider'))

    if show_first:
        ima = EMAN2_cppwrap.EMData()
        ima.read_image(stack, 0, False)
        plt.ion()
        plt.figure()
        plt.imshow(ima.get_2dview(), cmap=plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()

    return project_params_per_movie, particle_coordinates_per_movie, ctf_params_per_movie, \
           nx_per_movie, ny_per_movie, nz_per_movie, adnan_per_movie, chunk_per_movie

"""
Reading a reference map
"""
def get_2D_project_for_all_ptcl_from_reference(volume_ft , kb_fu, project_params_in_stack, frames_length , show = False):
    # print("Reference projections begin")
    project_2D_per_frame = []
    for i in tqdm(range(len(project_params_in_stack)), desc= "Reference Projection"):
        params_substack = project_params_in_stack[i].get_params('spider')
        params_for_each_image = [params_substack['phi'], params_substack['theta'], params_substack['psi'],
                                 params_substack['tx'] , params_substack['ty']  ]   #  + shifts[i][j][0]  + shifts[i][j][1]
        # project_2D_per_frame.append(sp_projection.prgs(volume_ft, kb_fu, params_for_each_image))

        project_2D_per_frame.append(sp_projection.prgl(volume_ft, params_for_each_image, interpolation_method=1))
        del params_substack
        del params_for_each_image

    return project_2D_per_frame
    # print("Reference projections end")


def givemeindices(img):
    w = img.shape[0] //2  + 1
    h = img.shape[1]
    x_index = np.arange(h)
    y_index = np.arange(w)
    x_index_expand = np.repeat(x_index, w)
    y_index_expand = np.tile(y_index, h)
    x_index_expand_cpy = x_index_expand.copy()
    x_index_expand_cpy[x_index_expand_cpy >= w] -= h
    indices = np.sqrt(y_index_expand ** 2 + x_index_expand_cpy ** 2).astype(int)
    mask = indices < w
    indices = indices[mask]
    x_index_expand = x_index_expand[mask]
    y_index_expand = y_index_expand[mask]

    del x_index
    del y_index
    del x_index_expand_cpy

    return indices, x_index_expand, y_index_expand

def givemeindices_new(img):
    w_half = img.shape[0] //2  + 1
    h,w = img.shape[1],img.shape[1]
    x_index = np.arange(h)
    y_index = np.arange(w)
    x_index_expand = np.repeat(x_index, w)
    y_index_expand = np.tile(y_index, h)
    x_index_expand_cpy = x_index_expand.copy()
    y_index_expand_cpy = y_index_expand.copy()
    x_index_expand_cpy[x_index_expand_cpy >= w_half] -= h
    y_index_expand_cpy[y_index_expand_cpy >= w_half] -= h
    indices = np.sqrt(y_index_expand_cpy ** 2 + x_index_expand_cpy ** 2).astype(int)
    mask = indices < w_half
    indices = indices[mask]
    x_index_expand = x_index_expand[mask]
    y_index_expand = y_index_expand[mask]

    return indices, x_index_expand, y_index_expand

def noiseNormalize (img, sigma2):
    area = 0.25 * np.pi * float(361 * 360)
    indx, xind, yind = givemeindices(img)
    for num in np.unique(indx):
        mask = num == indx
        if sigma2[num] == 0.0:
            img[xind[mask] , yind[mask]] = 0.0
        else :
            img[xind[mask], yind[mask]] /= np.sqrt(sigma2[num] * area)

    del indx
    del xind
    del yind
    return img

def noiseNormalize_new(img, sigma2):
    area = 0.25 * np.pi * float(361 * 360)
    indx, xind, yind = givemeindices_new(img)
    for num in np.unique(indx):
        mask = num == indx
        if sigma2[num] == 0.0:
            img[xind[mask], yind[mask]] = 0.0
        else:
            img[xind[mask], yind[mask]] /= np.sqrt(sigma2[num] * area)

    del indx
    del xind
    del yind
    return img


def powerspectrum_per_image(norm_img):
    w = norm_img.shape[0] //2  + 1
    wgh = np.zeros(w)
    outpu = np.zeros(w)
    indx , xind, yind = givemeindices(norm_img)
    for num in np.unique(indx):
        mask = num == indx
        outpu[num] += np.sum(norm_img.real[xind[mask], yind[mask]])
        wgh[num] += np.sum(mask)

    del indx
    del xind
    del yind
    return outpu, wgh

def movie_vari_powspect( movie, gainimage, variance ,count, zsize, part_cord,
                        nxx, nyy, cen_xx, cen_yy ) :
    wh = nxx // 2 + 1
    powspec  = np.zeros((len(part_cord), wh ))
    powspec_wgh = np.zeros((len(part_cord), wh))
    current_frame = EMData()
    angout_pix = 0.885
    for i in tqdm(range(zsize), desc="Movie Power Spectra"):
        current_frame.read_image(movie, i)
        current_frame = Util.muln_img(current_frame, gainimage)
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1, np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                w = fft_part_img.shape[0]
                h = fft_part_img.shape[1]
                fft_part_img[0, 0] = 0.0
                np.divide(fft_part_img, float(nxx * nyy), out=fft_part_img)
                scale = np.sqrt(np.divide(w * h * variance[j], count[j] * zsize))
                np.divide(fft_part_img, scale, out=fft_part_img)
                fft_part_img_cp = fft_part_img.copy()
                fft_part_img_cp = fft_part_img_cp.conj()
                weight0_power = np.multiply(fft_part_img, fft_part_img_cp)
                out_p, wgh_p  =  powerspectrum_per_image(weight0_power)
                powspec[j] += out_p
                powspec_wgh[j]  += wgh_p

                del ptcl
                del fft_part_img
                del scale
                del fft_part_img_cp
                del weight0_power
                del out_p
                del wgh_p
            except RuntimeError:
                continue
    del current_frame

    cum_powspec = np.sum(powspec,  axis=0)
    cum_powspec_wgh = np.sum(powspec_wgh, axis=0)
    np.divide(cum_powspec, cum_powspec_wgh, out=cum_powspec)

    del powspec
    del powspec_wgh
    del cum_powspec_wgh

    return cum_powspec


def calc_fsc_per_part(part_img ,ref_img, fcc, wg0, wg1, fr_no, part_no, index, x_index, y_index, sigma22, var_m , count_m, zsize):
    fft_part_img = np.fft.fft2(part_img.get_2dview())
    nx = fft_part_img.shape[0]
    ny = fft_part_img.shape[1]
    fft_part_img[0,0] = 0.0
    np.divide(fft_part_img, float(nx*ny), out = fft_part_img)
    scale = np.sqrt(np.divide(nx * ny * var_m[part_no], count_m[part_no] * zsize))
    np.divide(fft_part_img, scale, out=fft_part_img)
    fft_part_img = noiseNormalize(fft_part_img, sigma22)
    fft_ref_img = np.fft.fft2(ref_img.get_2dview())
    np.divide(fft_ref_img, float(nx*ny), out = fft_ref_img)
    weight1_power = np.abs(fft_ref_img) ** 2
    weight0_power = np.abs(fft_part_img) ** 2
    result_img = np.multiply(fft_part_img, fft_ref_img.conj())

    for num in np.unique(index):
        mask = num == index
        fcc[fr_no, num] += np.sum(result_img.real[y_index[mask], x_index[mask]])
        wg0[fr_no, num] += np.sum(weight0_power.real[y_index[mask], x_index[mask]])
        wg1[fr_no, num] += np.sum(weight1_power.real[y_index[mask], x_index[mask]])

    del fft_part_img
    del fft_ref_img
    del result_img
    del weight0_power
    del weight1_power
    del scale


def vectorize_indiv_var(fourier_img_norm, var_array, cnt_array, p_index):
    w = fourier_img_norm.shape[0]
    h = fourier_img_norm.shape[1]
    w = fourier_img_norm.shape[0] // 2 + 1

    x_index = np.arange(h)
    y_index = np.arange(w)
    x_index_expand = np.repeat(x_index, w)
    y_index_expand = np.tile(y_index, h)

    count_img = np.ones(w*h).reshape((h, w))
    count_img[0,0] = 0
    fourier_img_norm[:, y_index_expand[y_index_expand > 0]]  *= 2.0
    count_img[:, y_index_expand[y_index_expand > 0]] *= 2.0

    var_array[p_index] += np.sum(fourier_img_norm[x_index_expand, y_index_expand])
    cnt_array[p_index] += np.sum(count_img)

    del x_index
    del y_index
    del x_index_expand
    del y_index_expand
    del count_img

    return var_array, cnt_array


def init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref):
    current_frame = EMData()
    var_sum = np.zeros(len(part_cord))
    cnt_sum = np.zeros(len(part_cord))
    angout_pix = 0.885
    for i in tqdm(range(zsize) , desc="Initialize variance array"):
        current_frame.read_image(movie_name,i)
        current_frame = Util.muln_img(current_frame, gainref)
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1, np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                np.divide(fft_part_img, float(nxx * nyy), out=fft_part_img)
                fft_part_img[0, 0] = 0.0
                np.multiply(fft_part_img, fft_part_img.conj(), out=fft_part_img)
                var_sum, cnt_sum = vectorize_indiv_var(fft_part_img.real, var_sum, cnt_sum, j)
                del fft_part_img
                del ptcl
            except RuntimeError:
                continue
    del current_frame
    return var_sum , cnt_sum


def damageWeights(dosevalue, news, angpixval):
    a = 3.39999999999
    b = -1.060000000
    c = -0.540000000
    s = news
    kc = s//2 + 1
    angpix = angpixval
    rho = 1.0 / (2.0 * (kc - 1) * angpix)
    output = np.zeros((kc, s))
    for y in range(s):
        for x in range(kc):
            xf = x
            if y < kc :
                yf = y
            else:
                yf = y-s
            if x ==0 and y ==0:
                output[0, 0] = 1.0
            else:
                k = np.sqrt(xf*xf  +  yf*yf)
                tau = a * np.power(rho * k, b) + c
                output[x,y] = np.exp(np.divide(-dosevalue, tau))
    return output

def create_doses(dose_per_frame, no_of_frames, voltage):
    dose_frames = np.zeros(no_of_frames)
    frames_range = np.arange(no_of_frames)
    for iframe in range(no_of_frames):
        dose_frames[iframe] = dose_per_frame * (frames_range[iframe] + 1)
        if np.abs(voltage - 200) <= 2 :
            dose_frames[iframe] /= 0.8
        elif np.abs(voltage - 100) <= 2:
            dose_frames[iframe] /= 0.64
    return dose_frames


def givemerealbfactors(pixelval, fc, boxsize):
    # bfile = np.loadtxt('/home/adnan/PycharmProjects/Dose_Fresh/sphire/bin/bfactors.star')
    # slope = bfile[:,1]
    # intercept = bfile[:,2]

    bfile = np.loadtxt('/home/adnan/PycharmProjects/DoseWeighting/Polis001/metadata/BFactors.txt')
    slope = bfile[0]
    intercept = bfile[1]

    bfacy = np.exp(intercept)
    angpix_ref = pixelval
    sh_ref = boxsize//2 + 1
    cf = 8.0 * ((angpix_ref)**2)  * ((sh_ref)**2)
    bfacx = np.sqrt(np.divide(-cf, slope))
    del bfile
    return bfacx, bfacy


def computeweights(bfax, bfay, fc,  boxsize, normalize = True):
    kc2 = boxsize
    kc = boxsize//2 +1
    output = np.zeros((fc, kc, kc2))
    yy = 0
    for f in range (fc):
        for y in range (0, kc2):
            for x in range(0, kc):
                if y < kc :
                    yy = y
                else:
                    yy == y - kc2
                r2 = (x*x) + (yy*yy)
                output[f,x,y] =  bfay[f] * np.exp(np.divide((-0.5 * r2), bfax[f] * bfax[f]))
    if normalize :
        for y in range(0, kc2):
            for x in range(0, kc):
                sum = 0.0
                for f in range(fc):
                    sum += output[f,x,y]

                for f in range(fc):
                    output[f, x, y] /= sum
    return output


def angtopix(a , s, angpix):
    return s * angpix / a

def pixtoang(p, s, angpix):
    return s * angpix / p


def fitBkFactors(fcc_fsc, k0 , k1):
    kc = fcc_fsc.shape[1]
    fc = fcc_fsc.shape[0]

    maxSumf = 0.0
    bestF = 0
    steps = 20
    for f in range(fc):
        sumf = 0.0
        for k in range(kc):
            sumf += fcc_fsc[f,k]

        if (sumf > maxSumf):
            maxSumf = sumf
            bestF = f

    scale0 = np.zeros(kc)
    wgh = np.ones(kc)
    for k in range(kc):
        scale0[k] = max(0.0 , fcc_fsc[bestF, k])

    sigb = np.zeros((fc, 2))

    for it in range (5):
        for f in range (fc):
            q0 = 0.1
            depth0 = 4
            sig00 = 1.0
            sig11 = 10.0 * kc
            sigb[f] = findSigmaKRec(fcc_fsc, f , scale0, wgh, k0, k1, sig00, sig11, steps , depth0, q0)

        for k in range (kc):
            num = 0.0
            denom = 0.0
            for f in range (fc):
                p = fcc_fsc[f,k]
                q = sigb[f][1] *  np.exp( np.divide(-0.5 * k * k , (sigb[f][0] * sigb[f][0])))

                num += q * p
                denom += q*q

            eps = 1e-20
            if denom > eps :
                scale0[k] = np.divide(num , denom)
            else:
                scale0[k] = np.divide(num , eps)
            scale0[k] = max(0.0, scale0[k])
    return sigb , scale0


def findSigmaKRec(fcc_fit, f, envelope, weight, k0, k1, sig0, sig1, steps, depth, q):
    minErr = np.finfo(np.double).max
    bestSig = sig0
    bestA = 1.0
    eps = np.double(1e-20)

    for s in range(steps):
        sig = sig0 + np.divide(s * (sig1 - sig0), (steps - 1))
        sig2 = sig * sig
        if sig2 < eps:
            continue

        num = 0.0
        denom = 0.0
        for k in range(k0, k1):
            p = fcc_fit[f, k]
            qn = np.float(envelope[k] * np.exp(np.divide(-0.5 * k * k, sig2)))
            num += (qn * p)
            denom += (qn * qn)
        if denom > eps:
            a = np.divide(num, denom)
        else:
            a = np.divide(num, eps)
        sum = 0.0

        for k in range(k0, k1):
            d = fcc_fit[f, k] - envelope[k] * a * np.exp(np.divide(-0.5 * k * k, sig2))
            sum += (weight[k] * d * d)

        if sum < minErr:
            minErr = sum
            bestSig = sig
            bestA = a

    if depth > 0:
        hrange = 0.5 * (sig1 - sig0)
        snext0 = bestSig - q * hrange
        snext1 = bestSig + q * hrange
        if snext0 < eps:
            snext0 = eps
        return findSigmaKRec(fcc_fit, f, envelope, weight, k0, k1, snext0, snext1, steps, depth - 1, q)
    return bestSig, bestA

def findminmaxresol(fsc_path) :
    fsc_angs = np.loadtxt(fsc_path)
    angrsol = fsc_angs[:,1]
    fscvals = fsc_angs[:,2]
    for i in range(len(fsc_angs)) :
        if fscvals[i] < 0.00:
            print(fscvals[i])
            return angrsol[i]

def normalizeSigVel(sig_vel, doseperframe, angpix):
    try:
        params_scale_by_dose = doseperframe * sig_vel / angpix
    except:
        params_scale_by_dose = sig_vel / angpix

    return params_scale_by_dose

def normalizeSigDiv(sig_div, angpix):

    params_scale_by_dose = sig_div / angpix

    return params_scale_by_dose

def normalizeSigAcc(sig_acc, doseperframe, angpix):
    try:
        params_scale_by_dose = doseperframe * sig_acc / angpix
    except:
        params_scale_by_dose = sig_acc / angpix

    return params_scale_by_dose


def movieCC(movie_name, nx, ny, sigma22, var_m, count_m, gbdose, angout_pix, zsize, gainref, part_cord,
            ctf_params_part, cen_xx, cen_yy, refer_1, refer_2, CC_ptcl, chunk_arr):
    """
    For noise normalization for image and reference image
    """
    current_frame = EMData()
    angout_pix = 0.885

    for i in tqdm(range(zsize), desc="Cross correlation "):
        current_frame.read_image(movie_name,i)
        current_frame = Util.muln_img(current_frame, gainref)

        dose = gbdose * (i)
        damageweights = damageWeights(dose, nx, angout_pix)
        doseline = damageweights[:, 0].tolist()
        for j in range(len(part_cord)):
            try:

                ptcl = Util.window(current_frame, nx, ny, 1, np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                ptcl = sp_filter.filt_table(ptcl, doseline)
                if chunk_arr[j] == 0 :
                    refer = refer_1[j]
                else :
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)

                beam_tiltx = -0.096435999999999994
                beam_tilty = 0.26850000000000002
                wavelength = 0.025079493567048063
                Cs = 1.3999999999999999

                refer = apply_tilt(refer, beam_tiltx, beam_tilty, wavelength, Cs, angout_pix)

                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                # fft_part_img[0, 0] = 0.0
                np.divide(fft_part_img, float(nx * ny), out=fft_part_img)
                scale = np.sqrt(np.divide(nx * ny * var_m[j], count_m[j] * zsize))
                np.divide(fft_part_img, scale, out=fft_part_img)
                fft_part_img = noiseNormalize(fft_part_img, sigma22)

                fft_ref_img = np.fft.fft2(refer.get_2dview())
                np.divide(fft_ref_img, float(nx * ny), out=fft_ref_img)
                fft_ref_img[0, 0] = 0.0
                fft_ref_img = noiseNormalize(fft_ref_img, sigma22)

                np.multiply(fft_part_img, fft_ref_img.conjugate(), out = fft_part_img)

                ptcl_np_FF = np.fft.ifft2(fft_part_img)
                ptcl_np_FF = np.multiply(ptcl_np_FF.real, float(nx*nx*nx*nx))

                CC_ptcl[i,j] = ptcl_np_FF

                del ptcl_np_FF
                del fft_part_img
                del fft_ref_img
                del scale
                del refer
                del ptcl
            except RuntimeError:
                continue

    return CC_ptcl

def gpmotion(sig_vel, sig_acc, sig_div, globaldose, angout_pix, shift_x, shift_y, part_cord, zsize, correlate):
    sig_vel_px = normalizeSigVel(sig_vel, globaldose, angout_pix)
    sig_acc_px = normalizeSigAcc(sig_acc, globaldose, angout_pix)
    sig_div_px = normalizeSigDiv(sig_div, angout_pix)
    sv2 = sig_vel_px * sig_vel_px
    sd2 = sig_div_px * sig_div_px
    sd1 = sig_div_px
    pc = len(part_cord)
    A = np.zeros((pc, pc))

    expKer = True
    for i in range(pc):
        for j in range(i, pc):
            dd = np.sum(np.abs(np.array(part_cord[i]) - np.array(part_cord[j])) ** 2)  # pos = [x,y]
            if expKer == True:
                k = sv2 * np.exp(- np.sqrt(dd / sd2))
            else:
                k = sv2 * np.exp(-0.5 * dd / sd1)
            A[i, j] = k
            A[j, i] = k

    U, S, Vt = np.linalg.svd(A, full_matrices=True)
    eps = 1e-10

    dc = pc
    for d in range(pc):
        if (S[d] < eps):
            dc = d
            break

    # v_as_rel = [0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 40, 43, 47, 49,
    #             53, 55, 57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
    #             73, 74, 75, 76, 78, 79, 82, 83, 84, 87, 88, 91, 92, 95, 97, 99, 100,
    #             102, 104, 105, 107, 108, 110, 112, 115, 117, 119, 121, 123, 125, 127,
    #             128, 130, 132, 133, 134, 135, 137, 139, 141, 143, 145, 146, 150, 151,
    #             152, 153, 155, 158, 159, 161, 162, 166, 167, 171, 172, 173, 175, 177,
    #             179, 180, 183, 184, 185, 186, 188, 190, 192, 193, 194, 195, 197, 198,
    #             199, 202, 203, 208, 209]
    # for d in range(pc):
    #     if d in v_as_rel:
    #         Vt[d, :] = Vt[d, :]
    #     else:
    #         Vt[d, :] = Vt[d, :] * -1

    basis = np.zeros((dc, pc))
    eigenVals = np.zeros(dc)

    for d in range(dc):
        l = np.sqrt(S[d])
        for p in range(pc):
            basis[d, p] = l * Vt[d, p]

    for d in range(dc):
        eigenVals[d] = S[d]

    initialshifts = np.zeros((pc,zsize,2))
    for i in range(pc):
        initialshifts[i, :, :] = np.array([np.array(shift_x)[:]*-1, np.array(shift_y)[:]*-1 ]).swapaxes(0, 1)

    pad = 0
    newposi = np.zeros((pc, zsize + pad, 2))
    partCt = 2 * (pc + pc * (zsize - 1))
    initialCoeff = np.zeros(partCt)
    initialCoeff  = posToParams(initialshifts, initialCoeff, pc, zsize, basis, eigenVals, dc)

    # er = gpmotion_f(initialCoeff, newposi, pc, zsize, basis ,correlate, sig_acc_px, eigenVals)
    # grad_vec = gpmotion_grad(initialCoeff, newposi , pc , zsize, basis, correlate , sig_acc_px, eigenVals)

    # print(grad_vec)
    options_min = {}
    options_min['disp'] =  True
    options_min['maxcor']  = 200
    options_min['ftol'] = 0.0001
    options_min['gtol'] = 0.90
    options_min['eps'] =  1e-05
    options_min['iprint']  = 80
    options_min['maxfun'] =  10000
    options_min['maxiter'] = 10000
    options_min['maxls'] = 40

    lb = 9.99e-21
    ub = 10e19
    bou = (lb, ub)
    bb = np.tile(bou, partCt).reshape(partCt, 2).tolist()

    result = minimize(fun = gpmotion_f, x0= initialCoeff , args = (newposi, pc, zsize, basis ,correlate, sig_acc_px, eigenVals, dc) ,
                      method='L-BFGS-B', jac = gpmotion_grad, tol = 9.99e-17, options = options_min)

    # x, f ,d = fmin_l_bfgs_b(func = gpmotion_f, x0= initialCoeff , fprime = gpmotion_grad, args = (newposi, pc, zsize, basis ,correlate, sig_acc_px, eigenVals) ,
    #                   bounds = bb , m=6,   maxiter = 10000, iprint = 80, maxls = 40 )

    newshifts = np.zeros((pc, zsize, 2))
    newshifts = ParamsTopos(result.x ,newshifts ,pc, zsize, basis, dc )

    return newshifts

def posToParams(pos, x, pc, zsize, basis, eigenVals, dc):

    for q in range(pc):
        x[2 * q] = pos[q][0][0]
        x[2 * q + 1] = pos[q][0][1]

    for f in range(zsize - 1):
        for d in range(dc):
            cx, cy = 0.0, 0.0
            for p in range(pc):
                vx = pos[p][f+1][0]   - pos[p][f][0]
                vy = pos[p][f+1][1]  - pos[p][f][1]

                cx += vx * basis[d, p]
                cy += vy * basis[d, p]

            x[2 * (pc + dc * f + d)] = cx / eigenVals[d]
            x[2 * (pc + dc * f + d) + 1] = cy / eigenVals[d]

    return x

def ParamsTopos(x_vec, pos, pc, zsize, basis, dc ):

    for p in range(pc):
        pp = [x_vec[2*p],  x_vec[2*p+1]]
        for f in range(zsize):
            pos[p][f][0] = pp[0]
            pos[p][f][1] = pp[1]
            if f < zsize -1 :
                vx, vy  = 0.0, 0.0
                for d in range(dc):
                    cx = x_vec[2*(pc + dc*f + d)]
                    cy = x_vec[2*(pc + dc*f + d)+1]
                    vx += cx * basis[d,p]
                    vy += cy * basis[d,p]
                pp[0] += vx
                pp[1] += vy
    return pos

def gpmotion_f(x_vec, pos, pc, zsize, basis, correlation, sig_acc_px, eigenVals, dc):
    pos = ParamsTopos(x_vec, pos, pc, zsize, basis, dc)
    thre = 12
    e_t = np.zeros(thre)
    count =0

    for p in range(pc):
        for f in range (zsize):
            epf = sp_interpolation.cubicXY(correlation[p][f], pos[p][f][0], pos[p][f][1])
            e_t[count] -=  epf
        if count >= thre-1:
            count = 0
        else:
            count += 1
    count = 0

    for f in range(zsize-1):
        for d in range (dc):
            cx = x_vec[2* (pc + dc*f +d) ]
            cy = x_vec[2* (pc + dc*f +d) +1 ]
            e_t[count] += cx*cx  + cy*cy
        if count >= thre - 1:
            count = 0
        else:
            count += 1

    count = 0
    if ( sig_acc_px > 0.0):
        for f in range(zsize -2):
            for d in range(dc):
                cx0 = x_vec[2 * (pc + dc * f + d)]
                cy0 = x_vec[2 * (pc + dc * f + d) + 1]
                cx1 = x_vec[2 * (pc + dc * (f+1) + d)]
                cy1 = x_vec[2 * (pc + dc * (f+1) + d) + 1]
                dcx = cx1 - cx0
                dcy = cy1 - cy0
                e_t[count] += eigenVals[d]*(dcx*dcx + dcy*dcy) / (sig_acc_px*sig_acc_px)

            if count >= thre - 1:
                count = 0
            else:
                count += 1

    e_tot  = 0.0
    e_tot = sum(e_t)

    return e_tot

def gpmotion_grad(x_vec, pos, pc, zsize, basis, correlation,
                  sig_acc_px, eigenVals, dc):

    pad = 0
    partCt = len(x_vec)
    thre = 12
    gradDest  = np.zeros((thre, partCt+ pad) )
    ccg_pf = np.zeros((pc, zsize + pad, 2))

    pos = ParamsTopos(x_vec, pos, pc, zsize, basis, dc)
    for p in range(pc):
        for f in range (zsize):
            vrx , vry = sp_interpolation.cubicXYgrad(correlation[p][f], pos[p][f][0], pos[p][f][1])
            ccg_pf[p][f]  = [vrx , vry]


    count = 0
    for p in range(pc):
        for f in range(zsize):
            gradDest[count][2*p] -=   ccg_pf[p][f][0]
            gradDest[count][2*p+1] -= ccg_pf[p][f][1]
            if count >= thre - 1:
                count = 0
            else:
                count += 1

    count = 0
    for d in range(dc):
        for p in range(pc):
            bpd = basis[d, p]
            gx, gy = 0.0, 0.0
            for f in range(zsize -2 , -1, -1):
                gx += bpd * ccg_pf[p][f+1][0]
                gy += bpd * ccg_pf[p][f+1][1]
                gradDest[count][2 * (pc + dc * f + d) ]  -= gx
                gradDest[count][2 * (pc + dc * f + d) + 1] -= gy
            if count >= thre - 1:
                count = 0
            else:
                count += 1

    count = 0
    for f in range(zsize-1):
        for d in range(pc):
            gradDest[count][2 * (pc + dc * f + d) ] +=  2.0 *  x_vec[2 * (pc + dc * f + d) ]
            gradDest[count][2 * (pc + dc * f + d) + 1] +=  2.0 *  x_vec[2 * (pc + dc * f + d) +1 ]

        if count >= thre - 1:
            count = 0
        else:
            count += 1

    count = 0
    if ( sig_acc_px > 0.0):
        sa2 = sig_acc_px * sig_acc_px

        for f in range(zsize-2):
            for d in range(dc):
                cx0 = x_vec[2 * (pc + dc * f + d)]
                cy0 = x_vec[2 * (pc + dc * f + d) + 1]
                cx1 = x_vec[2 * (pc + dc * (f + 1) + d)]
                cy1 = x_vec[2 * (pc + dc * (f + 1) + d) + 1]

                dcx = cx1 - cx0
                dcy = cy1 - cy0

                gradDest[count][2 * (pc + dc * f + d)] -= 2.0 * eigenVals[d] * dcx / sa2
                gradDest[count][2 * (pc + dc * f + d)+1] -= 2.0 * eigenVals[d] * dcy / sa2
                gradDest[count][2 * (pc + dc * (f+1) + d)] += 2.0 * eigenVals[d] * dcx / sa2
                gradDest[count][2 * (pc + dc * (f+1) + d)+1] += 2.0 * eigenVals[d] * dcy / sa2

            if count >= thre - 1:
                count = 0
            else:
                count += 1

    gradDest_in = np.zeros(partCt)

    for t in range(thre):
        for i in range(partCt):
            gradDest_in[i] += gradDest[t][i]

    return gradDest_in



def get_fcc_all_particles(movie_name, refer_1, refer_2, nxx, nyy, part_cord, zsize,
                          gainrefname, chunk_arr, ctf_params_part, frame_shift_x, frame_shift_y, newshifts):

    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    cen_zz = zsize // 2
    freqs=  nyy//2 + 1

    fsc_values = np.zeros((zsize,freqs  ))
    weight_0 = np.zeros((zsize, freqs))
    weight_1 = np.zeros((zsize, freqs ))

    x_index = np.arange(freqs)
    y_index = np.arange(nyy)
    x_index_expand = np.repeat(x_index, nyy)
    y_index_expand = np.tile(y_index, freqs)
    y_index_expand_cpy = y_index_expand.copy()
    y_index_expand_cpy[y_index_expand_cpy >= freqs] -= nyy
    indices = np.round(np.sqrt(x_index_expand ** 2 + y_index_expand_cpy ** 2)).astype(int)
    mask = indices < freqs
    indices = indices[mask]
    x_index_expand = x_index_expand[mask]
    y_index_expand = y_index_expand[mask]

    del x_index
    del y_index
    del y_index_expand_cpy
    del mask

    gainref = EMData()
    gainref.read_image(gainrefname,0)
    Util.mul_scalar(gainref, -1.0)
    globaldose = 1.277
    angout_pix = 0.885
    sig_vel = 0.9435
    sig_acc = 2.55
    sig_div = 8235


    var_movie, cnt_movie =  init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
    sigma = movie_vari_powspect(movie_name,  gainref, var_movie, cnt_movie, zsize, part_cord,
                       nxx, nyy, cen_xx , cen_yy)

    CCFS = np.zeros((zsize, len(part_cord), nxx, nyy))
    CCFS = movieCC(movie_name, nxx, nyy, sigma, var_movie, cnt_movie, globaldose, angout_pix,
                   zsize, gainref, part_cord,ctf_params_part, cen_xx, cen_yy,
                   refer_1, refer_2, CCFS, chunk_arr)

    CCFS = CCFS.swapaxes(0, 1)

    newshifts  = gpmotion(sig_vel , sig_acc, sig_div, globaldose, angout_pix,frame_shift_x, frame_shift_y,
             part_cord, zsize, CCFS)

    for i in tqdm(range(zsize), desc="FCC Loop"):
        current_frame.read_image(movie_name,i)
        current_frame = Util.muln_img(current_frame, gainref)
        dose = globaldose * (i)
        damageweights = damageWeights(dose, nxx, angout_pix)
        doseline = damageweights[:, 0].tolist()
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1,
                                   np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                ptcl = applyshifts(ptcl, -newshifts[j][i][0], -newshifts[j][i][1])

                if chunk_arr[j] == 1:
                    refer = refer_1[j]
                else:
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)

                beam_tiltx = -0.096435999999999994
                beam_tilty = 0.26850000000000002
                wavelength = 0.025079493567048063
                Cs = 1.3999999999999999
                refer = apply_tilt(refer, beam_tiltx, beam_tilty, wavelength, Cs, angout_pix)

                calc_fsc_per_part(ptcl, refer, fsc_values, weight_0, weight_1, i, j, indices,
                                  x_index_expand, y_index_expand, sigma, var_movie, cnt_movie, zsize)
                del ptcl
                del refer
            except RuntimeError:
                continue

    del current_frame
    del var_movie
    del cnt_movie
    del gainref
    # del newshifts
    del sigma
    # del CCFS

    return fsc_values, weight_0, weight_1, newshifts


def get_all_polish_ptcl_imgs( movie_name,
                              nxx, nyy,
                              part_cord,
                              zsize,
                              gainrefname,
                              partshifts,
                              bfx,
                              bfy):
    """Getting dimension information"""
    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    gainref = EMData()
    gainref.read_image(gainrefname,0)
    Util.mul_scalar(gainref, -1.0)

    # bfx, bfy = givemerealbfactors(angout_pix, zsize, nxx)
    weights = computeweights(bfx, bfy, zsize, nxx, True)

    # var_movie, cnt_movie = init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
    #

    angout_pix = 0.885
    particle_imgs_in_movie = []
    for i in tqdm(range(zsize) , desc= "Applying Weights and get polish particles"):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)
        # current_frame = applyshifts(current_frame, frame_shift_x[i], frame_shift_y[i])
        line = weights[i, :, 0].tolist()
        crop_imgs = []
        new_coord = []
        for j in range(len(part_cord)):
            try:
                box_img = Util.window(current_frame, nxx, nyy, 1,
                                   np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                box_img = applyshifts(box_img, -partshifts[j][i][0], -partshifts[j][i][1])
                box_img = sp_filter.filt_table(box_img, line)
                
                crop_imgs.append(box_img)
                new_coord.append(j)
                del box_img
            except RuntimeError:
                # bad_particles.append(j)
                continue
        if i == zsize - 1:
            new_coord = new_coord
        else:
            del new_coord
        particle_imgs_in_movie.append(crop_imgs)
        del crop_imgs
        del line
    del current_frame
    del gainref
    del weights
    return particle_imgs_in_movie , new_coord


def apply_tilt(IN_img, tilt_x, tilt_y, wavelength, Cs, angpix, verbose=False):
    img = IN_img
    imgarray = img.get_2dview()
    imgfft = np.fft.fft2(imgarray)
    xdim = imgfft.shape[0]
    ydim = imgfft.shape[1]
    boxsize = angpix * xdim
    halfx = xdim // 2 + 1
    # Convert from milliradians to degrees
    factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize)
    for ycoord in range(ydim):
        if ycoord < halfx:  # will this work with non-square images?
            ip = ycoord
        else:
            ip = ycoord - ydim

        # Loop through x
        for xcoord in range(halfx):
            jp = xcoord

            # Calculate phase shift in degrees
            delta_phase = factor * (ip * ip + jp * jp) * (ip * tilt_y + jp * tilt_x)

            realval = imgfft[ycoord, xcoord].real
            imagval = imgfft[ycoord, xcoord].imag
            mag = np.sqrt((realval * realval) + (imagval * imagval))
            old_phase = np.arctan2(imagval, realval)

            # Convert phase shift from degrees to radians
            new_phase = old_phase + np.deg2rad(delta_phase)

            newrealval = mag * np.cos(new_phase)
            newimagval = mag * np.sin(new_phase)
            imgfft[ycoord, xcoord] = newrealval + (1j * newimagval)
            if verbose:
                # print(ycoord, xcoord, jp, ip, realval, imagval, mag, delta_phase, newrealval, newimagval)
                print(jp, ip, mag, np.rad2deg(old_phase), delta_phase, newrealval, newimagval)
            # print(jp, ip, realval, delta_phase)

    real_array = np.fft.ifft2(imgfft)
    eman_img = EMNumPy.numpy2em(real_array.real)
    return eman_img


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    'movies_folder',
    type=str,
    help='Movies folder'
)

parser.add_argument(
    'log_folder',
    type=str,
    help='Log files folder'
)

parser.add_argument(
    'stack_file',
    type=str,
    help='Substack file'
)

parser.add_argument(
    'reference_vol_0',
    type=str,
    help = 'reference volume 0'
)

parser.add_argument(
    'reference_vol_1',
    type=str,
    help = 'reference volume 1'
)

parser.add_argument(
    'ref_mask',
    type=str,
    help='Reference mask of the particle'
)

parser.add_argument(
    'saved_folder',
    type=str,
    help='Location of folder where new mrcs files and bdb will be stored'
)

parser.add_argument(
    'no_of_mic',
    type=int,
    help='Number of micrographs to be processed'
)

parser.add_argument(
    'ang_pix',
    type=float,
    help='Angular pixel value in angstorms'
)

parser.add_argument(
    'min_res',
    type=float,
    help='Minimum resolution value in angstorms'
)

parser.add_argument(
    'max_res',
    type=float,
    help='Maximum resolution value in angstorms'
)

parser.add_argument(
    'half_curve',
    type=str,
    help='FSC values from the Post refiner'
)

args = parser.parse_args()
no_of_micrographs = args.no_of_mic
Gainreffile = os.path.join( os.path.abspath(os.path.join(str(args.movies_folder), os.pardir)), 'gain.mrc')
movie_names = return_movie_names(str(args.movies_folder))

#-------------------------------Reading corrected sums log files for reading the global shifts applied
log_movie_path = args.log_folder
stackfilename = "bdb:" + args.stack_file
stack_absolute_path = args.saved_folder

print("Reading all the parameters from stack")
no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, \
            ctf_params_all, nx_all, ny_all, nz_all, chunks_all = read_all_attributes_from_stack(stackfilename)

print('no of images', no_of_imgs)
adnan_all = np.arange(EMUtil.get_image_count(stackfilename)).tolist()

print("Finding Non existing .mrc files in substack with respect to the movies")
movie_name_x = []
ptcl_source_images_xx =[]

for i in range(no_of_imgs):
    ptcl_source_images_xx.append(str(os.path.basename(ptcl_source_images[i])).split('.mrc')[0])

count = 0
for i in range (no_of_micrographs) :
    temp = np.where(np.array(ptcl_source_images_xx) == str(os.path.basename(movie_names[i])).split('.tiff')[0])
    if len(temp[0]) != 0:
        movie_name_x.append(movie_names[i])
        count+=1

movie_names = movie_name_x
no_of_micrographs = count

del ptcl_source_images_xx

ref_EM_1 = sp_utilities.get_im(args.reference_vol_0)
ref_EM_2 = sp_utilities.get_im(args.reference_vol_1)

mask_volume = sp_utilities.get_im(args.ref_mask)

ref_EM_1 *= mask_volume
ref_EM_2 *= mask_volume

volft_1 = sp_projection.prep_vol(ref_EM_1, npad=2, interpolation_method=1)
volft_2 = sp_projection.prep_vol(ref_EM_2, npad=2, interpolation_method=1)


kb_1 = 0
kb_2 = 0
del mask_volume
del ref_EM_1
del ref_EM_2

global_fcc = []
global_wg0 = []
global_wg1 = []
global_shift = []

min_res = args.min_res

if args.max_res == -1:
    max_res = findminmaxresol(args.half_curve)
else:
    max_res = args.max_res


# min_res = 20
max_res = 2.97757

ima_start, ima_end = sp_applications.MPI_start_end(no_of_micrographs, n_mpi_procs, my_mpi_proc_id)

for micro in enumerate(movie_names[ima_start:ima_end]):

    print("Applying GLOBAL shifts", micro)
    if str(os.path.basename(micro[1])).split('.')[-1] == 'mrcs':
        logfile = os.path.join(os.path.abspath(os.path.join(log_movie_path, os.pardir)),
                               micro[1].split('.')[0].split('/')[-1] + '.log')
        print(logfile)
    elif str(os.path.basename(micro[1])).split('.')[-1] == 'tiff':
        logfile = os.path.join(os.path.abspath(os.path.join(log_movie_path, os.pardir)),
                           micro[1].split('.')[0].split('/')[-1]  + '.star')
        print(logfile)

    zsize = EMUtil.get_image_count(micro[1])

    shift_x, shift_y = givemotioncorrshifts(logfile)

    trackfile = os.path.join(os.path.abspath(os.path.join(log_movie_path, os.pardir)),
                             micro[1].split('.')[0].split('/')[-1] + '_tracks' + '.star')

    try:
        shiftperptcl = giveshiftperparticles(trackfile, zsize)
    except:
        print("Per particles shifts files not available")

    # print(shift_x, shift_y)
    print("Sorted the particles from the stack ")
    project_params, particle_coordinates, ctf_params, \
    nx, ny, nz, adnan_n, chunks = find_particles_info_from_movie(stackfilename,
                                                                       micro[1],
                                                                       no_of_imgs,
                                                                       ptcl_source_images,
                                                                       project_params_all,
                                                                       particle_coordinates_all,
                                                                       ctf_params_all,
                                                                       nx_all,
                                                                       ny_all,
                                                                       nz_all,
                                                                       adnan_all,
                                                                       chunks_all,
                                                                       show_first=False)


    del adnan_n
    # del project_params
    # shift_file_name = os.path.join(os.path.join(str(args.saved_folder), 'metadata'),
    #                        micro[1].split('.')[0].split('/')[-1] + '_partshift.txt')

    # shiftperptcl = np.loadtxt(shift_file_name).reshape(len(particle_coordinates), zsize, 2)

    ref_project_2D_ptcl_vol1 = get_2D_project_for_all_ptcl_from_reference(volft_1, kb_1, project_params, zsize, show=False)
    ref_project_2D_ptcl_vol2 = get_2D_project_for_all_ptcl_from_reference(volft_2, kb_2, project_params, zsize, show=False)

    fsc_val, weight0, weight1, part_shifts = get_fcc_all_particles(micro[1], ref_project_2D_ptcl_vol1, ref_project_2D_ptcl_vol2,
                                                      nx[0], ny[0], particle_coordinates, zsize, Gainreffile,
                                                      chunks, ctf_params, shift_x, shift_y, shiftperptcl)

    np.savetxt(os.path.join(
        os.path.join(str(args.saved_folder), 'metadata' ),
        micro[1].split('.')[0].split('/')[-1] + '_FCC.txt'), fsc_val)

    np.savetxt(os.path.join(
        os.path.join(str(args.saved_folder), 'metadata'),
        micro[1].split('.')[0].split('/')[-1] + '_wg0.txt'), weight0)

    np.savetxt(os.path.join(
        os.path.join(str(args.saved_folder), 'metadata'),
        micro[1].split('.')[0].split('/')[-1] + '_wg1.txt'), weight1)

    count = 0
    f = open(os.path.join(os.path.join(str(args.saved_folder), 'metadata'),
                           micro[1].split('.')[0].split('/')[-1] + '_partshift.txt'), 'w')
    for sli in part_shifts:
        np.savetxt(f, sli, fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_' + str(count))
        count += 1
    f.close()

    # count = 0
    # f = open('test.txt', 'w')
    # for sli in a:
    #     np.savetxt(f, sli, fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_' + str(count))
    #     count+=1
    # f.close()

    global_fcc.append(fsc_val)
    global_wg0.append(weight0)
    global_wg1.append(weight1)
    global_shift.append(part_shifts)

    del fsc_val
    del weight0
    del weight1
    del part_shifts
    del ref_project_2D_ptcl_vol1
    del ref_project_2D_ptcl_vol2
    del particle_coordinates
    del ctf_params
    del shift_x
    del shift_y
    del project_params

fsc_values_per_micrograph = sp_utilities.wrap_mpi_gatherv(global_fcc, 0, mpi.MPI_COMM_WORLD)
weight0_per_micrograph = sp_utilities.wrap_mpi_gatherv(global_wg0, 0, mpi.MPI_COMM_WORLD)
weight1_per_micrograph = sp_utilities.wrap_mpi_gatherv(global_wg1, 0, mpi.MPI_COMM_WORLD)
part_shifts_per_micrograph = sp_utilities.wrap_mpi_gatherv(global_shift, 0, mpi.MPI_COMM_WORLD)



del global_fcc
del global_wg0
del global_wg1
del global_shift

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

if main_mpi_proc == my_mpi_proc_id :
    fsc_values_per_micrograph = np.array(fsc_values_per_micrograph)
    weight0_per_micrograph  = np.array(weight0_per_micrograph)
    weight1_per_micrograph = np.array(weight1_per_micrograph)

    fsc_sum_per_frame = np.sum(np.array(fsc_values_per_micrograph), axis=0)
    weight0_per_micrograph = np.sum(np.array(weight0_per_micrograph), axis=0)
    weight1_per_micrograph = np.sum(np.array(weight1_per_micrograph), axis=0)

    weight = np.sqrt(weight0_per_micrograph + weight1_per_micrograph)
    FCC = np.divide(fsc_sum_per_frame, weight)

    sh_ref = FCC.shape[1]
    s_ref = 2* (sh_ref - 1)

    k0 = int(angtopix(min_res, s_ref, args.ang_pix))
    k1 = int(angtopix(max_res, s_ref, args.ang_pix ))

    bfacsxx , bfacsyy = fitBkFactors(FCC, k0 , k1)

    cf =  8.0 * ((args.ang_pix)**2)  * ((sh_ref)**2)
    b = np.zeros(FCC.shape[0])
    k = np.zeros(FCC.shape[0])

    for f in range(FCC.shape[0]):
        s = bfacsxx[f, 0]
        b[f] = np.divide(-cf, s * s)
        k[f] = np.log(bfacsxx[f, 1])

    bfactors = np.array((b, k))
    np.savetxt(os.path.join(
        os.path.join(str(args.saved_folder), 'metadata'),
        'BFactors.txt'), bfactors)

    del b
    del k
    bfacxx = bfacsxx[:, 0].tolist()
    bfacsyy = bfacsxx[:, 1].tolist()

else:
    bfacxx = []
    bfacsyy = []

bfacx = sp_utilities.bcast_list_to_all(bfacxx, my_mpi_proc_id, main_mpi_proc, mpi.MPI_COMM_WORLD)
bfacy = sp_utilities.bcast_list_to_all(bfacsyy, my_mpi_proc_id, main_mpi_proc, mpi.MPI_COMM_WORLD)




bfacx,bfacy = np.loadtxt('/home/adnan/PycharmProjects/DoseWeighting/Polis001/metadata/BFactors.txt')

for micro in enumerate(movie_names[ima_start:ima_end]):
    # print("Applying GLOBAL shifts", micro[0])
    # if str(os.path.basename(micro[1])).split('.')[-1] == 'mrcs':
    #     logfile = os.path.join(os.path.abspath(os.path.join(log_movie_path, os.pardir)),
    #                            micro[1].split('.')[0].split('/')[-1] + '.log')
    #     print(logfile)
    # elif str(os.path.basename(micro[1])).split('.')[-1] == 'tiff':
    #     logfile = os.path.join(os.path.abspath(os.path.join(log_movie_path, os.pardir)),
    #                        micro[1].split('.')[0].split('/')[-1]  + '.star')
    #     print(logfile)

    # shift_x, shift_y = givemotioncorrshifts(logfile)
    zsize = EMUtil.get_image_count(micro[1])



    print("Sorted the particles from the stack ")
    project_params, particle_coordinates, ctf_params, \
    nx, ny, nz,  adnan_n, chunks = find_particles_info_from_movie(stackfilename,
                                                                               micro[1],
                                                                               no_of_imgs,
                                                                               ptcl_source_images,
                                                                               project_params_all,
                                                                               particle_coordinates_all,
                                                                               ctf_params_all,
                                                                               nx_all,
                                                                               ny_all,
                                                                               nz_all,
                                                                               adnan_all,
                                                                               chunks_all,
                                                                               show_first=False)


    del project_params
    del ctf_params
    del chunks
    print("Applying dose weighting")

    shift_file_name = os.path.join(os.path.join(str(args.saved_folder), 'metadata'),
                           micro[1].split('.')[0].split('/')[-1] + '_partshift.txt')

    shiftperptcl = np.loadtxt(shift_file_name).reshape(len(particle_coordinates), zsize, 2)

    particle_imgs_dosed, old_ind_coord = get_all_polish_ptcl_imgs(micro[1],
                                                                    nx[0],
                                                                    ny[0],
                                                                    particle_coordinates,
                                                                    zsize,
                                                                    Gainreffile,
                                                                    shiftperptcl,
                                                                    bfacx,
                                                                    bfacy)


    particle_imgs_dosed = np.array(particle_imgs_dosed).swapaxes(0, 1)

    print("Dose weighting done, summing starts")
    mask = sp_utilities.model_circle(nx[0] // 2, nx[0], nx[0])
    ave_particle_dosed = []
    for i in range(len(particle_imgs_dosed)):
        ave_particle_dosed.append(sum(particle_imgs_dosed[i]))
    del particle_imgs_dosed
    del particle_coordinates

    print("Writing into mrcs files", len(ave_particle_dosed))

    local_stack_path = "bdb:%s" % os.path.join(str(os.path.join(str(stack_absolute_path), 'micrographs')) , micro[1].split('.')[0].split('/')[-1] + "_ptcls")
    local_mrc_path = os.path.join(str(os.path.join(str(stack_absolute_path), 'micrographs')) , micro[1].split('.')[0].split('/')[-1] + "_ptcls.mrcs")
    local_bdb_stack = db_open_dict(local_stack_path)
    old_stack = db_open_dict(stackfilename, ro=True)
    print('Polish particles saved in ', local_mrc_path)

    for i in range(len(ave_particle_dosed)):
        index_old = adnan_n[old_ind_coord[i]]
        old_dict = old_stack.get(index_old, nodata=True).get_attr_dict()
        old_dict['data_path'] = local_mrc_path
        old_dict['data_n'] = int(i)
        old_dict['source_n'] = int(i)
        old_dict['ptcl_source_coord_id'] = i
        local_bdb_stack[i] = old_dict
        ave_particle_dosed[i] = sp_fundamentals.resample(ave_particle_dosed[i], 1 )  # (0.885/1.24)
        st = Util.infomask(ave_particle_dosed[i], mask, False)
        ave_particle_dosed[i] += 2 * st[0]
        st = Util.infomask(ave_particle_dosed[i], mask, False)
        ave_particle_dosed[i] -= st[0]
        ave_particle_dosed[i] /= st[1]
        ave_particle_dosed[i].append_image(local_mrc_path)

        db_close_dict(local_stack_path)
        db_close_dict(stackfilename)

    del local_bdb_stack
    del old_stack
    del ave_particle_dosed
    del old_ind_coord

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
mpi.mpi_finalize()



"""
I have placed it here so that i dont confuse two things.
def get_fcc_all_particles(movie_name, refer_1, refer_2, nxx, nyy, part_cord, zsize,
                          gainrefname, part_shifts, chunk_arr, ctf_params_part, frame_shift_x, frame_shift_y):

    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    cen_zz = zsize // 2
    freqs=  nyy//2 + 1

    fsc_values = np.zeros((zsize,freqs  ))
    weight_0 = np.zeros((zsize, freqs))
    weight_1 = np.zeros((zsize, freqs ))

    x_index = np.arange(freqs)
    y_index = np.arange(nyy)
    x_index_expand = np.repeat(x_index, nyy)
    y_index_expand = np.tile(y_index, freqs)
    y_index_expand_cpy = y_index_expand.copy()
    y_index_expand_cpy[y_index_expand_cpy >= freqs] -= nyy
    indices = np.round(np.sqrt(x_index_expand ** 2 + y_index_expand_cpy ** 2)).astype(int)
    mask = indices < freqs
    indices = indices[mask]
    x_index_expand = x_index_expand[mask]
    y_index_expand = y_index_expand[mask]

    del x_index
    del y_index
    del y_index_expand_cpy
    del mask

    gainref = EMData()
    gainref.read_image(gainrefname,0)
    Util.mul_scalar(gainref, -1.0)
    globaldose = 1.277
    angout_pix = 0.885

"""
#Calculating the fourier shell correlation of all the particle images with respect to 2-D reference

"""
Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
creating a window around to extract the same particle from each frame
"""
"""
#----------------------- Particle cordinate
def get_all_polish_ptcl_imgs( movie_name,
                              nxx, nyy,
                              part_cord,
                              zsize,
                              gainrefname,
                              partshifts,
                              bfx,
                              bfy,
                              frame_shift_x,
                              frame_shift_y):
"""
#Getting dimension information
"""
    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    gainref = EMData()
    gainref.read_image(gainrefname,0)
    Util.mul_scalar(gainref, -1.0)

    # bfx, bfy = givemerealbfactors(angout_pix, zsize, nxx)
    weights = computeweights(bfx, bfy, zsize, nxx, True)

    particle_imgs_in_movie = []
    for i in tqdm(range(zsize) , desc= "Applying Weights and get polish particles"):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)
        current_frame = applyshifts(current_frame, frame_shift_x[i], frame_shift_y[i])
        line = weights[i, :, 0].tolist()
        crop_imgs = []
        new_coord = []
        for j in range(len(part_cord)):
            try:
                box_img = Util.window(current_frame, nxx, nyy, 1, part_cord[j][0] - cen_xx,
                                             part_cord[j][1] - cen_yy ,0)
                # box_img = applyshifts(box_img, -partshifts[j][i][0], -partshifts[j][i][1])
                box_img = sp_filter.filt_table(box_img, line)
                crop_imgs.append(box_img)
                new_coord.append(j)
                del box_img
            except RuntimeError:
                # bad_particles.append(j)
                continue
        if i == zsize - 1:
            new_coord = new_coord
        else:
            del new_coord
        particle_imgs_in_movie.append(crop_imgs)
        del crop_imgs
        del line
    del current_frame
    del gainref
    del weights
    return particle_imgs_in_movie , new_coord
"""


"""
movies_path = '/home/adnan/PycharmProjects/DoseWeighting/MOVIES_RELION/'
movie = os.path.join(movies_path, '20170629_00021_frameImage.tiff')
logfilepath = '/home/adnan/PycharmProjects/DoseWeighting/Relion_log/'
logfile = os.path.join(logfilepath, '20170629_00021_frameImage.star')
logfile_track = os.path.join(logfilepath, '20170629_00021_frameImage_tracks.star')
Gainreffile = os.path.join(movies_path, 'gain.mrc')


stackfilename = "bdb:" + '/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv7/MotionCorr/job030/MOVIES_RELION/Particles/sphire_relion_stack'

print("Reading all the parameters from stack")
source_n_ind_all = EMUtil.get_all_attributes(stackfilename, "source_n")
adnan_all = EMUtil.get_all_attributes(stackfilename, 'adnan_n')
no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, \
            ctf_params_all, nx_all, ny_all, nz_all, chunks_all = read_all_attributes_from_stack(stackfilename)


ref_vol_path1 = '/home/adnan/PycharmProjects/DoseWeighting/run_half1_class001_unfil.mrc'
ref_vol_path2 = '/home/adnan/PycharmProjects/DoseWeighting/run_half2_class001_unfil.mrc'
mask_path = '/home/adnan/PycharmProjects/DoseWeighting/mask.mrc'

ref_EM_1 = sp_utilities.get_im(ref_vol_path1)
ref_EM_2 = sp_utilities.get_im(ref_vol_path2)

mask_volume = sp_utilities.get_im(mask_path)

ref_EM_1 *= mask_volume
ref_EM_2 *= mask_volume

volft_1 = sp_projection.prep_vol(ref_EM_1, npad=2, interpolation_method= 1)
volft_2  = sp_projection.prep_vol(ref_EM_2, npad=2, interpolation_method= 1)

kb_1 = 0
kb_2 = 0

del mask_volume
del ref_EM_1
del ref_EM_2
"""


####

# import mrcfile
# import matplotlib.pyplot as plt
# def readimg(path):
#     with mrcfile.open(path) as mrc:
#         dat = mrc.data
#     return dat
#
#
# ref_imgs = readimg('/home/adnan/PycharmProjects/newrelion/Polish/job013/reference_00021.mrcs')

#
# plt.figure()
# plt.imshow(ref_imgs[0])