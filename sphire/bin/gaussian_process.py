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

from scipy.fftpack import ifft2 as iff
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
# widgets = ['Progress: ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
#            ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
#
# pbar = ProgressBar(widgets=widgets, maxval=UnknownLength)

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


no_of_micrographs = 17
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
    particle_number = int(df['data_general'].get_values()[0].rsplit()[1])
    shiftingperptcl = np.zeros((particle_number, frames, 2))

    for i in range(particle_number):
        for j in range(frames):
            if i == 0:
                k = 5 + j
            if i > 0:
                k = i * frames + ((4 * i) + j) + 5
            shiftingperptcl[i, j, 0] = float(df['data_general'].get_values()[k].rsplit()[0])
            shiftingperptcl[i, j, 1] = float(df['data_general'].get_values()[k].rsplit()[1])
    return shiftingperptcl

def givemotioncorrshifts(filename):
    shift_file = pd.read_csv(filename)
    header = shift_file['data_general'][0:16]
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
        # print("Dimensions x for all particles are ", nx_per_movie)
        # print("Dimensions y for all particles are ", ny_per_movie)
        # print("Dimensions z for all particles are ", nz_per_movie)
    elif str(os.path.basename(movie_name)).split('.')[-1] == 'tiff':
        # print('Ctf shape',np.array(ctf_params_per_movie).shape)
        print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)),
                                                             len(project_params_per_movie)))
        print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
        print("Projection parameters for 1st particle in the stack are ", project_params_per_movie[0].get_params('spider'))
        # print("Dimensions x for all particles are ", nx_per_movie[0])
        # print("Dimensions y for all particles are ", ny_per_movie[0])
        # print("Dimensions z for all particles are ", nz_per_movie[0])

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
    # project_2D_per_movie = []
    # pbar.start()
    # for j in tqdm(range(frames_length) , desc= "Reference Projection"):
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
        # project_2D_per_movie.append(project_2D_per_frame)
        # pbar.update()
        # del project_2D_per_frame
    # pbar.finish()
    # print("Reference projections end")

    # if show:
    #     plt.ion()
    #     plt.figure()
    #     plt.imshow(project_2D_per_movie[0][0].get_2dview(), cmap = plt.get_cmap('Greys'))
    #     plt.colorbar()
    #     plt.show()
    # return project_2D_per_movie

def givemeindices(img):
    w = img.shape[0] //2  + 1
    h = img.shape[1]
    x_index = np.arange(h)
    y_index = np.arange(w)
    x_index_expand = np.repeat(x_index, w)
    y_index_expand = np.tile(y_index, h)
    x_index_expand_cpy = x_index_expand.copy()
    # y_index_expand_cpy = y_index_expand.copy()
    x_index_expand_cpy[x_index_expand_cpy >= w] -= h
    # y_index_expand_cpy[y_index_expand_cpy >= w] -= h
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


def extract_image(inp_img, pixel_size, boxsize, cord_x, cord_y ):
    inp = inp_img.get_2dview()
    out =  np.zeros((boxsize,boxsize))
    xp0 = np.int(pixel_size * cord_x / pixel_size) - boxsize//2
    yp0 = np.int(pixel_size * cord_y / pixel_size) - boxsize//2

    x0 = np.int(np.round(xp0 *pixel_size/ pixel_size))
    y0 = np.int(np.round(yp0 *pixel_size/ pixel_size))

    for y in range(boxsize):
        for x in range(boxsize):
            xx = x0+x
            yy = y0+y
            if (xx < 0):
                xx = 0
            elif (xx >= inp.shape[0]):
                xx = inp.shape[0] - 1
            else:
                xx =xx

            if (yy < 0):
                yy = 0
            elif (yy >= inp.shape[1]):
                yy = inp.shape[1] - 1
            else:
                yy =yy

            out[y][x] = inp[yy][xx]


    out_em = EMNumPy.numpy2em(out)
    return out_em










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
    sh_ref = boxsize/2 + 1
    cf = 8.0 * ((angpix_ref)**2)  * ((sh_ref)**2)
    bfacx = np.sqrt(np.divide(-cf, slope))
    del bfile
    return bfacx, bfacy

# @jit(nopython=True)
def computeweights(bfax, bfay, fc,  boxsize, normalize = True):
    kc2 = boxsize
    kc = boxsize/2 +1
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

def readimg(path):
    with mrcfile.open(path) as mrc:
        dat = mrc.data
    return dat

# @jit(parallel= True)
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

def estim_one_particle(movie_name, refer_1, refer_2, nxx, nyy, part_cord, zsize,
                          gainrefname, chunk_arr, ctf_params_part, frame_shift_x, frame_shift_y):

    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    cen_zz = zsize // 2
    freqs=  nyy//2 + 1
    # fsc_values = np.zeros((zsize,freqs  ))
    # weight_0 = np.zeros((zsize, freqs))
    # weight_1 = np.zeros((zsize, freqs ))
    # x_index = np.arange(freqs)
    # y_index = np.arange(nyy)
    # x_index_expand = np.repeat(x_index, nyy)
    # y_index_expand = np.tile(y_index, freqs)
    # y_index_expand_cpy = y_index_expand.copy()
    # y_index_expand_cpy[y_index_expand_cpy >= freqs] -= nyy
    # indices = np.round(np.sqrt(x_index_expand ** 2 + y_index_expand_cpy ** 2)).astype(int)
    # mask = indices < freqs
    # indices = indices[mask]
    # x_index_expand = x_index_expand[mask]
    # y_index_expand = y_index_expand[mask]
    # del x_index
    # del y_index
    # del y_index_expand_cpy
    # del mask

    gainref = EMData()
    gainref.read_image(gainrefname,0)
    Util.mul_scalar(gainref, -1.0)
    globaldose = 1.277
    angout_pix = 0.885


    """
    Calculating the fourier shell correlation of all the particle images with respect to 2-D reference 
    projection of 3-D volume
    """

    ref_imgs = readimg('/home/adnan/PycharmProjects/newrelion/Polish/job013/reference_00021.mrcs')

    var_movie, cnt_movie =  init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
    sigma = movie_vari_powspect(movie_name,  gainref, var_movie, cnt_movie, zsize, part_cord,
                       nxx, nyy, cen_xx , cen_yy)
    #
    # CCFS = np.zeros((zsize, len(part_cord), nxx, nyy))
    # CCFS = movieCC(movie_name, nxx, nyy, sigma, var_movie, cnt_movie, globaldose, angout_pix,
    #                zsize, gainref, part_cord,ctf_params_part, cen_xx, cen_yy,
    #                refer_1, refer_2, CCFS, chunk_arr, ref_imgs)
    #
    # CCFS = CCFS.swapaxes(0,1)

    peak_values = np.zeros((len(part_cord), zsize, 2))

    for i in tqdm(range(zsize), desc= "Normalize noise in image and reference image"):
        current_frame.read_image(movie_name,i)
        current_frame = Util.muln_img(current_frame, gainref)
        current_frame = applyshifts(current_frame , frame_shift_x[i] , frame_shift_y[i] )

        dose = globaldose * (i)
        damageweights = damageWeights(dose, nxx, angout_pix)
        doseline = damageweights[:, 0].tolist()

        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1, part_cord[j][0] - cen_xx,
                                             part_cord[j][1] - cen_yy ,0)

                if chunk_arr[j] == 1 :
                    refer = refer_1[j]
                else:
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)

                ptcl = sp_filter.filt_table(ptcl, doseline)

                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                np.divide(fft_part_img, float(nxx * nyy), out=fft_part_img)
                scale = np.sqrt(np.divide(nxx * nyy * var_movie[j], cnt_movie[j] * zsize))
                np.divide(fft_part_img, scale, out=fft_part_img)
                fft_part_img = noiseNormalize_new(fft_part_img, sigma)

                fft_ref_img = np.fft.fft2(refer.get_2dview())
                np.divide(fft_ref_img, float(nxx * nyy), out=fft_ref_img)
                fft_ref_img[0, 0] = 0.0
                fft_ref_img = noiseNormalize_new(fft_ref_img, sigma)

                ptcl = EMNumPy.numpy2em(np.fft.ifft2(fft_part_img).real)
                refer = EMNumPy.numpy2em(np.fft.ifft2(fft_ref_img).real)

                ccmap = sp_fundamentals.ccf(ptcl, refer)
                peaks = sp_utilities.peak_search(ccmap, 3)

                peak_values[j][i][0] = peaks[0][-2]
                peak_values[j][i][1] = peaks[0][-1]

                # ptcl = numpy2em_python(fft_part_img)
                # refer = numpy2em_python(fft_ref_img)
                # ccmap = sp_fundamentals.ccf(ptcl , refer)
                # peak_values = sp_utilities.peak_search(ccmap, 3)
                # ptcl = applyshifts(ptcl , -peak_values[0][-2], -peak_values[0][-1])
                # box_img = sp_filter.filt_table(box_img, doseline)
                # calc_fsc_per_part(ptcl, refer, fsc_values, weight_0, weight_1, i, j, indices,
                #                   x_index_expand, y_index_expand, sigma, var_movie, cnt_movie, zsize)
                del ptcl
                del refer
            except RuntimeError:
                continue
    del current_frame
    # del var_movie
    # del cnt_movie
    del gainref
    # del sigma


    return peak_values

    # return fsc_values, weight_0, weight_1



def movieCC(movie_name, nx, ny, sigma22, var_m, count_m, gbdose, angout_pix, zsize, gainref, part_cord,
            ctf_params_part, cen_xx, cen_yy, refer_1, refer_2, CC_ptcl, chunk_arr, rel_ref):
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
                # ptcl = sp_filter.filt_table(ptcl, doseline)
                # if chunk_arr[j] == 0 :
                #     refer = refer_1[j]
                # else :
                #     refer = refer_2[j]
                # refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=0, binary=False)

                refer = rel_ref[j]


                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                np.divide(fft_part_img, float(nx * ny), out=fft_part_img)
                # fft_part_img[0, 0] = 0.0
                scale = np.sqrt(np.divide(nx * ny * var_m[j], count_m[j] * zsize))
                np.divide(fft_part_img, scale, out=fft_part_img)
                fft_part_img = noiseNormalize(fft_part_img, sigma22)

                fft_ref_img = np.fft.fft2(refer)
                np.divide(fft_ref_img, float(nx * ny), out=fft_ref_img)
                fft_ref_img[0, 0] = 0.0
                fft_ref_img = noiseNormalize_new(fft_ref_img, sigma22)

                np.multiply(fft_part_img, fft_ref_img.conjugate(), out = fft_part_img)

                ptcl_np_FF = np.fft.ifft2(fft_part_img)
                ptcl_np_FF = np.multiply(ptcl_np_FF.real, float(nx*nx*nx*nx) )

                CC_ptcl[i,j] = ptcl_np_FF

                del ptcl_np_FF
                del fft_part_img
                del fft_ref_img
                del scale
                del refer
                del ptcl
            except RuntimeError:
                print("Particle skipped")
                continue

    return CC_ptcl


movies_path = '/home/adnan/PycharmProjects/DoseWeighting/MOVIES_RELION/'
movie = os.path.join(movies_path, '20170629_00021_frameImage.tiff')
logfilepath = '/home/adnan/PycharmProjects/DoseWeighting/Relion_log/'
logfile = os.path.join(logfilepath, '20170629_00021_frameImage.star')
logfile_track = os.path.join(logfilepath, '20170629_00021_frameImage_tracks.star')
Gainreffile = os.path.join(movies_path, 'gain.mrc')


stackfilename = "bdb:" + '/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv9/MotionCorr/job030/MOVIES_RELION/Particles/sphire_relion_stack'

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

min_res = 20
max_res = 2.97757
# if args.max_res == -1:
#     max_res = findminmaxresol(args.half_curve)
# else:
#     max_res = args.max_res

# ref_volume_1 = readimg(ref_vol_path1)
# ref_volume_2 = readimg(ref_vol_path2)
#
# mask_volume = readimg(mask_path)
#
# ref_volume_1 = np.multiply(ref_volume_1, mask_volume)
# ref_volume_2 = np.multiply(ref_volume_2, mask_volume)
#
# ref_EM_1 = EMData()
# ref_EM_2 = EMData()
#
# ref_EM_1 = EMNumPy.numpy2em(ref_volume_1)
# ref_EM_2 = EMNumPy.numpy2em(ref_volume_2)
#
# volft_1, kb_1 = sp_projection.prep_vol(ref_EM_1, npad=2, interpolation_method=-1)
# volft_2, kb_2 = sp_projection.prep_vol(ref_EM_2, npad=2, interpolation_method=-1)

zsize = EMUtil.get_image_count(movie)
# shiftperptcl = giveshiftperparticles(logfile_track, zsize)
# print(shiftperptcl)

shift_file = pd.read_csv(logfile)
header = shift_file['data_general'][0:16]
shift_x = []
shift_y = []
for i in range(17, len(shift_file['data_general'])):
    shift_x.append(float(shift_file['data_general'][i].split()[1]))
    shift_y.append(float(shift_file['data_general'][i].split()[2]))

# print(shift_x, shift_y)

print("Sorted the particles from the stack ")
project_params, particle_coordinates, ctf_params, \
nx, ny, nz, adnan_n, chunks = find_particles_info_from_movie(stackfilename,
                                                                   movie,
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


ref_project_2D_ptcl_vol1 = get_2D_project_for_all_ptcl_from_reference(volft_1, kb_1, project_params, zsize, show=False)
ref_project_2D_ptcl_vol2 = get_2D_project_for_all_ptcl_from_reference(volft_2, kb_2, project_params, zsize, show=False)


peak_values = estim_one_particle( movie,
                    ref_project_2D_ptcl_vol1,
                    ref_project_2D_ptcl_vol2,
                    nx[0], ny[0],
                    particle_coordinates,
                    zsize,
                    Gainreffile,
                    chunks,
                    ctf_params,
                    shift_x, shift_y)

count = 0
f = open('newcordinates.txt', 'w')
for sli in peak_values:
    np.savetxt(f, sli, fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_' + str(count))
    count += 1
f.close()

print('CTF on Reference projections')






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
# projection of 3-D volume
"""
    var_movie, cnt_movie =  init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
    sigma = movie_vari_powspect(movie_name,  gainref, var_movie, cnt_movie, zsize, part_cord,
                       nxx, nyy, cen_xx , cen_yy)

    for i in tqdm(range(zsize), desc= "FCC Loop"):
        # bad_particles = []
        current_frame.read_image(movie_name,i)
        current_frame = Util.muln_img(current_frame, gainref)
        current_frame = applyshifts(current_frame , frame_shift_x[i] , frame_shift_y[i] )
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1, part_cord[j][0] - cen_xx,
                                             part_cord[j][1] - cen_yy ,0)

                # ptcl = applyshifts(ptcl, -part_shifts[j][i][0], -part_shifts[j][i][1])
                if chunk_arr[j] == 1 :
                    refer = refer_1[j]
                else :
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)

                # ccmap = sp_fundamentals.ccf(ptcl , refer)
                #
                # peak_values = sp_utilities.peak_search(ccmap, 3)
                #
                # ptcl = applyshifts(ptcl , -peak_values[0][-2], -peak_values[0][-1])

                calc_fsc_per_part(ptcl, refer, fsc_values, weight_0, weight_1, i, j, indices,
                                  x_index_expand, y_index_expand, sigma, var_movie, cnt_movie, zsize)
                del ptcl
                del refer
            except RuntimeError:
                # bad_particles.append(j)
                continue
        # del bad_particles
    del current_frame
    del var_movie
    del cnt_movie
    del gainref
    # del bad_particles
    del sigma

    return fsc_values, weight_0, weight_1
"""

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