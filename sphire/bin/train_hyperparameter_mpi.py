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
import random

location = os.getcwd()
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
rank = my_node_proc_id

no_of_micrographs = 24
N_ITER = 25
shift = 1

try:
    ABSOLUTE_PATH = os.path.dirname(os.path.realpath(__file__))
    GLOBAL_PATH = os.path.abspath(os.path.join(__file__, "../../../"))
except:
    GLOBAL_PATH = os.getcwd()


def readimg(path):
    with mrcfile.open(path) as mrc:
        dat = mrc.data
    return dat


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


def returns_values_in_file(f, mode='r'):
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
    print("ERROR> the given file '" + str(f) + "' is not present!")
    exit(-1)


def read_meta_shifts(f):
    x = []
    y = []
    for row in returns_values_in_file(f):
        if "image #" in row:
            v = row.replace('\n', '').split('=')[1].replace(' ', '').split(',')
            x.append(float(v[0]))
            y.append(float(v[1]))
        elif "Frame (" in row:
            v = row.split()
            x.append(float(v[-2]))
            y.append(float(v[-1]))
    return x, y


def giveshiftperparticles(trackingfile, frames):
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


def applyshifts(image, shifx, shify):
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
           particle_coordinates_all_once, ctf_params_all_once, nx_all_once, ny_all_once, nz_all_once, chunk_params_all_once


def find_particles_info_from_movie(stack, movie_name, no_of_imgs, ptcl_source_images, project_params_all,
                                   particle_coordinates_all, ctf_params_all, nx_all, ny_all, nz_all,
                                   adnan_n_all, chunk_params_all, show_first=False):
    # -----------------------------------------------   CTF related attributes
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
                    str(os.path.basename(ptcl_source_images[i])) + 's'
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
        print("Projection parameters for 1st particle in the stack are ",
              project_params_per_movie[0].get_params('spider'))

    elif str(os.path.basename(movie_name)).split('.')[-1] == 'tiff':
        print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)),
                                                             len(project_params_per_movie)))
        print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
        print("Projection parameters for 1st particle in the stack are ",
              project_params_per_movie[0].get_params('spider'))

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


def get_2D_project_for_all_ptcl_from_reference(volume_ft, kb_fu, project_params_in_stack, frames_length, show=False):
    # print("Reference projections begin")
    project_2D_per_frame = []
    for i in tqdm(range(len(project_params_in_stack)), desc="Reference Projection"):
        params_substack = project_params_in_stack[i].get_params('spider')
        params_for_each_image = [params_substack['phi'], params_substack['theta'], params_substack['psi'],
                                 params_substack['tx'], params_substack['ty']]  # + shifts[i][j][0]  + shifts[i][j][1]
        # project_2D_per_frame.append(sp_projection.prgs(volume_ft, kb_fu, params_for_each_image))

        project_2D_per_frame.append(sp_projection.prgl(volume_ft, params_for_each_image, interpolation_method=1))
        del params_substack
        del params_for_each_image

    return project_2D_per_frame
    # print("Reference projections end")


def givemeindices(img):
    w = img.shape[0] // 2 + 1
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


def noiseNormalize(img, sigma2):
    area = 0.25 * np.pi * float(361 * 360)
    indx, xind, yind = givemeindices(img)
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
    w = norm_img.shape[0] // 2 + 1
    wgh = np.zeros(w)
    outpu = np.zeros(w)
    indx, xind, yind = givemeindices(norm_img)
    for num in np.unique(indx):
        mask = num == indx
        outpu[num] += np.sum(norm_img.real[xind[mask], yind[mask]])
        wgh[num] += np.sum(mask)

    del indx
    del xind
    del yind
    return outpu, wgh


def movie_vari_powspect(movie, gainimage, variance, count, zsize, part_cord,
                        nxx, nyy, cen_xx, cen_yy):
    wh = nxx // 2 + 1
    powspec = np.zeros((len(part_cord), wh))
    powspec_wgh = np.zeros((len(part_cord), wh))
    current_frame = EMData()
    angout_pix = 0.885
    for i in tqdm(range(zsize), desc="Movie Power Spectra"):
        current_frame.read_image(movie, i)
        current_frame = Util.muln_img(current_frame, gainimage)
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1,
                                   np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
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
                out_p, wgh_p = powerspectrum_per_image(weight0_power)
                powspec[j] += out_p
                powspec_wgh[j] += wgh_p

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

    cum_powspec = np.sum(powspec, axis=0)
    cum_powspec_wgh = np.sum(powspec_wgh, axis=0)
    np.divide(cum_powspec, cum_powspec_wgh, out=cum_powspec)

    del powspec
    del powspec_wgh
    del cum_powspec_wgh

    return cum_powspec


def calc_fsc_per_part(part_img, ref_img, fcc, wg0, wg1, fr_no, part_no, index, x_index, y_index, sigma22, var_m,
                      count_m, zsize):
    fft_part_img = np.fft.fft2(part_img.get_2dview())
    nx = fft_part_img.shape[0]
    ny = fft_part_img.shape[1]
    fft_part_img[0, 0] = 0.0
    np.divide(fft_part_img, float(nx * ny), out=fft_part_img)
    scale = np.sqrt(np.divide(nx * ny * var_m[part_no], count_m[part_no] * zsize))
    np.divide(fft_part_img, scale, out=fft_part_img)
    fft_part_img = noiseNormalize(fft_part_img, sigma22)
    fft_ref_img = np.fft.fft2(ref_img.get_2dview())
    np.divide(fft_ref_img, float(nx * ny), out=fft_ref_img)
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

    count_img = np.ones(w * h).reshape((h, w))
    count_img[0, 0] = 0
    fourier_img_norm[:, y_index_expand[y_index_expand > 0]] *= 2.0
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
    for i in tqdm(range(zsize), desc="Initialize variance array"):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nxx, nyy, 1,
                                   np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
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
    return var_sum, cnt_sum


def damageWeights(dosevalue, news, angpixval):
    a = 3.39999999999
    b = -1.060000000
    c = -0.540000000
    s = news
    kc = s // 2 + 1
    angpix = angpixval
    rho = 1.0 / (2.0 * (kc - 1) * angpix)
    output = np.zeros((kc, s))
    for y in range(s):
        for x in range(kc):
            xf = x
            if y < kc:
                yf = y
            else:
                yf = y - s
            if x == 0 and y == 0:
                output[0, 0] = 1.0
            else:
                k = np.sqrt(xf * xf + yf * yf)
                tau = a * np.power(rho * k, b) + c
                output[x, y] = np.exp(np.divide(-dosevalue, tau))
    return output


def create_doses(dose_per_frame, no_of_frames, voltage):
    dose_frames = np.zeros(no_of_frames)
    frames_range = np.arange(no_of_frames)
    for iframe in range(no_of_frames):
        dose_frames[iframe] = dose_per_frame * (frames_range[iframe] + 1)
        if np.abs(voltage - 200) <= 2:
            dose_frames[iframe] /= 0.8
        elif np.abs(voltage - 100) <= 2:
            dose_frames[iframe] /= 0.64
    return dose_frames


def givemerealbfactors(pixelval, fc, boxsize):
    bfile = np.loadtxt('/home/adnan/PycharmProjects/DoseWeighting/Polis001/metadata/BFactors.txt')
    slope = bfile[0]
    intercept = bfile[1]
    bfacy = np.exp(intercept)
    angpix_ref = pixelval
    sh_ref = boxsize // 2 + 1
    cf = 8.0 * ((angpix_ref) ** 2) * ((sh_ref) ** 2)
    bfacx = np.sqrt(np.divide(-cf, slope))
    del bfile
    return bfacx, bfacy


def computeweights(bfax, bfay, fc, boxsize, normalize=True):
    kc2 = boxsize
    kc = boxsize // 2 + 1
    output = np.zeros((fc, kc, kc2))
    yy = 0
    for f in range(fc):
        for y in range(0, kc2):
            for x in range(0, kc):
                if y < kc:
                    yy = y
                else:
                    yy == y - kc2
                r2 = (x * x) + (yy * yy)
                output[f, x, y] = bfay[f] * np.exp(np.divide((-0.5 * r2), bfax[f] * bfax[f]))
    if normalize:
        for y in range(0, kc2):
            for x in range(0, kc):
                sum = 0.0
                for f in range(fc):
                    sum += output[f, x, y]

                for f in range(fc):
                    output[f, x, y] /= sum
    return output


def angtopix(a, s, angpix):
    return s * angpix / a


def pixtoang(p, s, angpix):
    return s * angpix / p


def fitBkFactors(fcc_fsc, k0, k1):
    kc = fcc_fsc.shape[1]
    fc = fcc_fsc.shape[0]

    maxSumf = 0.0
    bestF = 0
    steps = 20
    for f in range(fc):
        sumf = 0.0
        for k in range(kc):
            sumf += fcc_fsc[f, k]

        if (sumf > maxSumf):
            maxSumf = sumf
            bestF = f

    scale0 = np.zeros(kc)
    wgh = np.ones(kc)
    for k in range(kc):
        scale0[k] = max(0.0, fcc_fsc[bestF, k])

    sigb = np.zeros((fc, 2))

    for it in range(5):
        for f in range(fc):
            q0 = 0.1
            depth0 = 4
            sig00 = 1.0
            sig11 = 10.0 * kc
            sigb[f] = findSigmaKRec(fcc_fsc, f, scale0, wgh, k0, k1, sig00, sig11, steps, depth0, q0)

        for k in range(kc):
            num = 0.0
            denom = 0.0
            for f in range(fc):
                p = fcc_fsc[f, k]
                q = sigb[f][1] * np.exp(np.divide(-0.5 * k * k, (sigb[f][0] * sigb[f][0])))

                num += q * p
                denom += q * q

            eps = 1e-20
            if denom > eps:
                scale0[k] = np.divide(num, denom)
            else:
                scale0[k] = np.divide(num, eps)
            scale0[k] = max(0.0, scale0[k])
    return sigb, scale0


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


def findminmaxresol(fsc_path):
    fsc_angs = np.loadtxt(fsc_path)
    angrsol = fsc_angs[:, 1]
    fscvals = fsc_angs[:, 2]
    for i in range(len(fsc_angs)):
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
            ctf_params_part, cen_xx, cen_yy, refer_1, refer_2, CC_ptcl, chunk_arr, damageweights):
    """
    For noise normalization for image and reference image
    """
    current_frame = EMData()
    angout_pix = 0.885

    for i in tqdm(range(zsize), desc="Cross correlation "):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)

        # dose = gbdose * (i)
        # damageweights = damageWeights(dose, nx, angout_pix)
        doseline = damageweights[i, :, 0].tolist()
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nx, ny, 1, np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                ptcl = sp_filter.filt_table(ptcl, doseline)
                if chunk_arr[j] == 0:
                    refer = refer_1[j]
                else:
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)
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
                np.multiply(fft_part_img, fft_ref_img.conjugate(), out=fft_part_img)
                ptcl_np_FF = np.fft.ifft2(fft_part_img)
                ptcl_np_FF = np.multiply(ptcl_np_FF.real, float(nx * nx * nx * nx))

                CC_ptcl[i, j] = ptcl_np_FF

                del ptcl_np_FF
                del fft_part_img
                del fft_ref_img
                del scale
                del refer
                del ptcl
            except RuntimeError:
                continue

    return CC_ptcl


def movieCC_train(movie_name, nx, ny, sigma22, var_m, count_m, gbdose, angout_pix, zsize, gainref, part_cord,
                  ctf_params_part, cen_xx, cen_yy, refer_1, refer_2, CC_ptcl, movie_total, chunk_arr, damageweights,
                  accPix, accCords, maxRangeP):
    """
    For noise normalization for image and reference image
    """
    current_frame = EMData()
    angout_pix = 0.885
    for i in tqdm(range(zsize), desc="Cross correlation "):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)
        doseline = damageweights[i, :, 0].tolist()
        for j in range(len(part_cord)):
            try:
                ptcl = Util.window(current_frame, nx, ny, 1, np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                   np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)
                ptcl = sp_filter.filt_table(ptcl, doseline)

                ptcl2 = Util.window(current_frame, nx, ny, 1,
                                    np.int(angout_pix * part_cord[j][0] / angout_pix) - cen_xx,
                                    np.int(angout_pix * part_cord[j][1] / angout_pix) - cen_yy, 0)

                if chunk_arr[j] == 0:
                    refer = refer_1[j]
                else:
                    refer = refer_2[j]

                refer = sp_filter.filt_ctf(refer, ctf_params_part[j], sign=1, binary=False)
                fft_part_img = np.fft.fft2(ptcl.get_2dview())
                np.divide(fft_part_img, float(nx * ny), out=fft_part_img)
                scale = np.sqrt(np.divide(nx * ny * var_m[j], count_m[j] * zsize))
                np.divide(fft_part_img, scale, out=fft_part_img)
                fft_part_img = noiseNormalize(fft_part_img, sigma22)

                fft_part_img_2 = np.fft.fft2(ptcl2.get_2dview())
                np.divide(fft_part_img_2, float(nx * ny), out=fft_part_img_2)
                scale = np.sqrt(np.divide(nx * ny * var_m[j], count_m[j] * zsize))
                np.divide(fft_part_img_2, scale, out=fft_part_img_2)
                fft_part_img_2 = noiseNormalize(fft_part_img_2, sigma22)
                movie_total[i, j] = accelerate(fft_part_img_2, accPix, accCords)

                fft_ref_img = np.fft.fft2(refer.get_2dview())
                np.divide(fft_ref_img, float(nx * ny), out=fft_ref_img)
                fft_ref_img[0, 0] = 0.0
                fft_ref_img = noiseNormalize(fft_ref_img, sigma22)

                np.multiply(fft_part_img, fft_ref_img.conjugate(), out=fft_part_img)
                ptcl_np_FF = np.fft.ifft2(fft_part_img)
                ptcl_np_FF = np.multiply(ptcl_np_FF.real, float(nx * nx * nx * nx) // 10)
                CC_ptcl[i, j] = crop_corner(ptcl_np_FF, maxRangeP, maxRangeP)

                del ptcl_np_FF
                del fft_part_img
                del fft_ref_img
                del fft_part_img_2
                del scale
                del refer
                del ptcl
                del ptcl2
            except RuntimeError:
                continue

    return CC_ptcl, movie_total


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

    basis = np.zeros((dc, pc))
    eigenVals = np.zeros(dc)

    for d in range(dc):
        l = np.sqrt(S[d])
        for p in range(pc):
            basis[d, p] = l * Vt[d, p]

    for d in range(dc):
        eigenVals[d] = S[d]

    initialshifts = np.zeros((pc, zsize, 2))
    for i in range(pc):
        initialshifts[i, :, :] = np.array([np.array(shift_x)[:] * -1, np.array(shift_y)[:] * -1]).swapaxes(0, 1)

    pad = 0
    newposi = np.zeros((pc, zsize + pad, 2))
    partCt = 2 * (pc + pc * (zsize - 1))
    initialCoeff = np.zeros(partCt)
    initialCoeff = posToParams(initialshifts, initialCoeff, pc, zsize, basis, eigenVals, dc)

    options_min = {}
    options_min['disp'] = True
    options_min['maxcor'] = 200
    options_min['ftol'] = 0.0001
    options_min['gtol'] = 0.90
    options_min['eps'] = 1e-05
    options_min['iprint'] = 80
    options_min['maxfun'] = 10000
    options_min['maxiter'] = 10000
    options_min['maxls'] = 40

    lb = 9.99e-21
    ub = 10e19
    bou = (lb, ub)
    bb = np.tile(bou, partCt).reshape(partCt, 2).tolist()

    result = minimize(fun=gpmotion_f, x0=initialCoeff,
                      args=(newposi, pc, zsize, basis, correlate, sig_acc_px, eigenVals, dc),
                      method='L-BFGS-B', jac=gpmotion_grad, tol=9.99e-17, options=options_min)

    newshifts = np.zeros((pc, zsize, 2))
    newshifts = ParamsTopos(result.x, newshifts, pc, zsize, basis, dc)
    return newshifts


def posToParams(pos, x, pc, zsize, basis, eigenVals, dc):
    for q in range(pc):
        x[2 * q] = pos[q][0][0]
        x[2 * q + 1] = pos[q][0][1]
    for f in range(zsize - 1):
        for d in range(dc):
            cx, cy = 0.0, 0.0
            for p in range(pc):
                vx = pos[p][f + 1][0] - pos[p][f][0]
                vy = pos[p][f + 1][1] - pos[p][f][1]
                cx += vx * basis[d, p]
                cy += vy * basis[d, p]
            x[2 * (pc + dc * f + d)] = cx / eigenVals[d]
            x[2 * (pc + dc * f + d) + 1] = cy / eigenVals[d]
    return x


def ParamsTopos(x_vec, pos, pc, zsize, basis, dc):
    for p in range(pc):
        pp = [x_vec[2 * p], x_vec[2 * p + 1]]
        for f in range(zsize):
            pos[p][f][0] = pp[0]
            pos[p][f][1] = pp[1]
            if f < zsize - 1:
                vx, vy = 0.0, 0.0
                for d in range(dc):
                    cx = x_vec[2 * (pc + dc * f + d)]
                    cy = x_vec[2 * (pc + dc * f + d) + 1]
                    vx += cx * basis[d, p]
                    vy += cy * basis[d, p]
                pp[0] += vx
                pp[1] += vy
    return pos


def gpmotion_f(x_vec, pos, pc, zsize, basis, correlation, sig_acc_px, eigenVals, dc):
    pos = ParamsTopos(x_vec, pos, pc, zsize, basis, dc)
    thre = 12
    e_t = np.zeros(thre)
    count = 0
    for p in range(pc):
        for f in range(zsize):
            epf = sp_interpolation.cubicXY(correlation[p][f], pos[p][f][0], pos[p][f][1])
            e_t[count] -= epf
        if count >= thre - 1:
            count = 0
        else:
            count += 1
    count = 0

    for f in range(zsize - 1):
        for d in range(dc):
            cx = x_vec[2 * (pc + dc * f + d)]
            cy = x_vec[2 * (pc + dc * f + d) + 1]
            e_t[count] += cx * cx + cy * cy
        if count >= thre - 1:
            count = 0
        else:
            count += 1
    count = 0
    if (sig_acc_px > 0.0):
        for f in range(zsize - 2):
            for d in range(dc):
                cx0 = x_vec[2 * (pc + dc * f + d)]
                cy0 = x_vec[2 * (pc + dc * f + d) + 1]
                cx1 = x_vec[2 * (pc + dc * (f + 1) + d)]
                cy1 = x_vec[2 * (pc + dc * (f + 1) + d) + 1]
                dcx = cx1 - cx0
                dcy = cy1 - cy0
                e_t[count] += eigenVals[d] * (dcx * dcx + dcy * dcy) / (sig_acc_px * sig_acc_px)
            if count >= thre - 1:
                count = 0
            else:
                count += 1
    e_tot = 0.0
    e_tot = sum(e_t)
    return e_tot


def gpmotion_grad(x_vec, pos, pc, zsize, basis, correlation,
                  sig_acc_px, eigenVals, dc):
    pad = 0
    partCt = len(x_vec)
    thre = 12
    gradDest = np.zeros((thre, partCt + pad))
    ccg_pf = np.zeros((pc, zsize + pad, 2))
    pos = ParamsTopos(x_vec, pos, pc, zsize, basis, dc)
    for p in range(pc):
        for f in range(zsize):
            vrx, vry = sp_interpolation.cubicXYgrad(correlation[p][f], pos[p][f][0], pos[p][f][1])
            ccg_pf[p][f] = [vrx, vry]
    count = 0
    for p in range(pc):
        for f in range(zsize):
            gradDest[count][2 * p] -= ccg_pf[p][f][0]
            gradDest[count][2 * p + 1] -= ccg_pf[p][f][1]
            if count >= thre - 1:
                count = 0
            else:
                count += 1
    count = 0
    for d in range(dc):
        for p in range(pc):
            bpd = basis[d, p]
            gx, gy = 0.0, 0.0
            for f in range(zsize - 2, -1, -1):
                gx += bpd * ccg_pf[p][f + 1][0]
                gy += bpd * ccg_pf[p][f + 1][1]
                gradDest[count][2 * (pc + dc * f + d)] -= gx
                gradDest[count][2 * (pc + dc * f + d) + 1] -= gy
            if count >= thre - 1:
                count = 0
            else:
                count += 1
    count = 0
    for f in range(zsize - 1):
        for d in range(pc):
            gradDest[count][2 * (pc + dc * f + d)] += 2.0 * x_vec[2 * (pc + dc * f + d)]
            gradDest[count][2 * (pc + dc * f + d) + 1] += 2.0 * x_vec[2 * (pc + dc * f + d) + 1]
        if count >= thre - 1:
            count = 0
        else:
            count += 1
    count = 0
    if (sig_acc_px > 0.0):
        sa2 = sig_acc_px * sig_acc_px

        for f in range(zsize - 2):
            for d in range(dc):
                cx0 = x_vec[2 * (pc + dc * f + d)]
                cy0 = x_vec[2 * (pc + dc * f + d) + 1]
                cx1 = x_vec[2 * (pc + dc * (f + 1) + d)]
                cy1 = x_vec[2 * (pc + dc * (f + 1) + d) + 1]
                dcx = cx1 - cx0
                dcy = cy1 - cy0
                gradDest[count][2 * (pc + dc * f + d)] -= 2.0 * eigenVals[d] * dcx / sa2
                gradDest[count][2 * (pc + dc * f + d) + 1] -= 2.0 * eigenVals[d] * dcy / sa2
                gradDest[count][2 * (pc + dc * (f + 1) + d)] += 2.0 * eigenVals[d] * dcx / sa2
                gradDest[count][2 * (pc + dc * (f + 1) + d) + 1] += 2.0 * eigenVals[d] * dcy / sa2
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
                          gainrefname, chunk_arr, ctf_params_part, frame_shift_x, frame_shift_y):
    current_frame = EMData()
    current_frame.read_image(movie_name, 0)
    xsize = current_frame.get_xsize()
    ysize = current_frame.get_ysize()
    cen_xx = xsize // 2
    cen_yy = ysize // 2
    cen_zz = zsize // 2
    freqs = nyy // 2 + 1
    fsc_values = np.zeros((zsize, freqs))
    weight_0 = np.zeros((zsize, freqs))
    weight_1 = np.zeros((zsize, freqs))
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
    gainref.read_image(gainrefname, 0)
    Util.mul_scalar(gainref, -1.0)
    globaldose = 1.277
    angout_pix = 0.885
    sig_vel = 0.9435
    sig_acc = 2.55
    sig_div = 8235

    alignDoseWeights = np.zeros((zsize, boxsize // 2 + 1, boxsize))
    for i in range(zsize):
        dose = globaldose * (i)
        damageweights = damageWeights(dose, boxsize, angout_pix)
        alignDoseWeights[i] = damageweights

    var_movie, cnt_movie = init_variance_array(movie_name, nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
    sigma = movie_vari_powspect(movie_name, gainref, var_movie, cnt_movie, zsize, part_cord,
                                nxx, nyy, cen_xx, cen_yy)

    CCFS = np.zeros((zsize, len(part_cord), nxx, nyy))
    CCFS = movieCC(movie_name, nxx, nyy, sigma, var_movie, cnt_movie, globaldose, angout_pix,
                   zsize, gainref, part_cord, ctf_params_part, cen_xx, cen_yy,
                   refer_1, refer_2, CCFS, chunk_arr, alignDoseWeights)
    CCFS = CCFS.swapaxes(0, 1)

    newshifts = gpmotion(sig_vel, sig_acc, sig_div, globaldose, angout_pix, frame_shift_x, frame_shift_y,
                         part_cord, zsize, CCFS)

    for i in tqdm(range(zsize), desc="FCC Loop"):
        current_frame.read_image(movie_name, i)
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
    del sigma
    del CCFS

    return fsc_values, weight_0, weight_1, newshifts


def get_all_polish_ptcl_imgs(movie_name,
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
    gainref.read_image(gainrefname, 0)
    Util.mul_scalar(gainref, -1.0)

    weights = computeweights(bfx, bfy, zsize, nxx, True)
    angout_pix = 0.885
    particle_imgs_in_movie = []
    for i in tqdm(range(zsize), desc="Applying Weights and get polish particles"):
        current_frame.read_image(movie_name, i)
        current_frame = Util.muln_img(current_frame, gainref)
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
    return particle_imgs_in_movie, new_coord


def MotionParamEstima(movie_names, stackfilename, boxsize, zsize, angout_pix, mc,
                      volft_1, volft_2, gainrefname, logfile_loc):
    align_frac = 0.5
    eval_frac = 0.5
    minParticles = 5000
    Sv_init = 0.6
    Sd_init = 10000
    Sa_init = 3
    inistep = 3000
    conv = 30
    maxIters = 100
    maxRange = 50
    seed = 23

    kout = 106

    k_cutoff = np.int(kout * np.sqrt(align_frac) + 0.5)
    k_eval = np.int(kout * np.sqrt(1.0 - eval_frac) + 0.5)

    k_cutoff_angst = angtopix(k_cutoff, boxsize, angout_pix)
    k_eval_angst = angtopix(k_eval, boxsize, angout_pix)

    print("Maximum frequency to consider the alignment {} A ({} px)".format(k_cutoff_angst, k_cutoff))

    print("frequency range to consider for evaluation {} - {} A ({} - {} px)".format(k_eval_angst,
                                                                                     pixtoang(kout, boxsize,
                                                                                              angout_pix), k_eval,
                                                                                     kout))
    mc = 24  # micrgraph count
    random.seed(23)
    randNums = np.zeros(mc)

    for m in range(mc):
        randNums[m] = random.random()

    order = np.argsort(randNums)
    pc = 0

    sel_files = []
    movie_files = []

    for i in range(len(order)):
        m = order[i]
        stackname = stackfilename.split(stackfilename.split('/')[-1])[0] + \
                    movie_names[m].split('/')[-1].split('.tiff')[0] + '_ptcls'
        pcm = EMUtil.get_image_count(stackname)

        pc += pcm
        sel_files.append(stackname)
        movie_files.append(movie_names[m])

        if pc > minParticles:
            print("{} particles found in {} micrographs".format(pc, mc))
            break

    mc = len(order)
    globaldose = 1.277
    dose_data, accPix, accCords = prepAlignment(globaldose, angout_pix, boxsize, zsize,
                                                                           k_cutoff, k_eval + 2, kout,
                                                                           mc, sel_files, gainrefname, movie_files,
                                                                           maxRange, volft_1, volft_2)

    initial = motiontoProblem([Sv_init, Sd_init, Sa_init])

    n = len(initial)
    m = n + 1
    simplex = np.zeros((m, n))
    simplex[0] = initial
    initialStep = 3000

    for j in range(1, m):
        simplex[j] = initial
        simplex[j][j - 1] += initialStep

    options_min = {}
    options_min['disp'] = True
    options_min['maxiter'] = 100
    options_min['maxfev'] = 100
    options_min['initial_simplex'] = simplex
    sc , tt = evaluateParams(initial, globaldose, angout_pix, mc, sel_files, movie_files, accPix, accCords, dose_data, logfile_loc, 2*maxRange )

    print("1st iteration values ", sc)
    print("1st iteration result", tt)


    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_finalize()




    # optimize call
    vd = minimize(fun=evaluateParams, x0=initial, tol= 10e-8, args=(globaldose, angout_pix, mc, sel_files, movie_files,
                                                                accPix, accCords, dose_data,
                                                                logfile_loc, 2*maxRange), method='Nelder-Mead', options=options_min)

    print("Result of parameter is ", vd.x)
    newvd = problemtoMotion(vd.x)
    velscale = 10000.0
    divScale = 1.0
    accScale = 1000.0
    conv = 30
    nrm = [newvd[0] * velscale, newvd[1] * divScale, newvd[2] * accScale]

    rnd = [conv * 0.5 * (int(2.0 * nrm[0] / conv + 0.5)) / velscale,
           conv * 0.5 * (int(2.0 * nrm[1] / conv + 0.5)) / divScale,
           conv * 0.5 * (int(2.0 * nrm[2] / conv + 0.5)) / accScale]

    print("nrm values is ", nrm)
    print("Final parameters is", rnd)


def updateTsc(tracks, pc, zsize, boxsize, accPix, accCords, obs, pred, damage):
    outT = [0.0, 0.0, 0.0]
    pc = tracks.shape[0]
    sh = boxsize // 2 + 1
    for p in range(pc):
        for f in range(zsize):
            shifts = [0.0, 0.0]
            shifts[0] = tracks[p][f][0] / boxsize
            shifts[1] = tracks[p][f][1] / boxsize
            for i in range(accPix):
                x = accCords[i][0]
                if accCords[i][1] < sh:
                    y = accCords[i][1]
                else:
                    y = accCords[i][1] - boxsize

                dotp = 2 * np.pi * (x * shifts[0] + y * shifts[1])
                b = np.sin(dotp)
                a = np.cos(dotp)

                t2 = obs[p][f][i]
                c = t2.real
                d = t2.imag
                ac = a * c
                bd = b * d
                ab_cd = (a + b) * (c + d)

                z_obs = (ac - bd) + (ab_cd - ac - bd) * 1j

                zpred = pred[p][i]
                dmg = damage[f][i]

                outT[0] += dmg * (zpred.real * z_obs.real + zpred.imag * z_obs.imag)
                outT[1] += dmg * (np.abs(z_obs) ** 2)
                outT[2] += dmg * (np.abs(zpred) ** 2)
    return outT


def evaluateParams(x, dose, angpix, gc, sel_files, movie_all,  accPix, accCords, damage,
                   logfile_loc, maxRangeP):
    vda = problemtoMotion(x)

    sig_v_vals_px = normalizeSigVel(vda[0], dose, angpix)
    sig_d_vals_px = normalizeSigDiv(vda[1], angpix)
    sig_a_vals_px = normalizeSigAcc(vda[2], dose, angpix)
    sv2 = sig_v_vals_px * sig_v_vals_px
    sd2 = sig_d_vals_px * sig_d_vals_px
    sd1 = sig_d_vals_px

    pctot = 0
    tscAs = [0.0, 0.0, 0.0]

    ima_start, ima_end = sp_applications.MPI_start_end(gc, n_mpi_procs, my_mpi_proc_id)

    for micro in enumerate(sel_files[ima_start:ima_end]):
        g = sel_files.index(micro[1])
        pc = EMUtil.get_image_count(micro[1])
        part_cord = EMUtil.get_all_attributes(micro[1], "ptcl_source_coord")
        if pc < 2:
            continue

        zsize = EMUtil.get_image_count(movie_all[g])


        saved_folder = '/home/adnan/PycharmProjects/DoseWeighting/Polis001/'
        corr_file_name = os.path.join(os.path.join(str(saved_folder), 'metadata'),
                                       micro[1].split('.')[0].split('/')[-1] + '_train_CC.txt')

        correlate = np.loadtxt(corr_file_name)

        correlate = correlate.reshape(pc, zsize, maxRangeP, maxRangeP)

        movie_file_name = os.path.join(os.path.join(str(saved_folder), 'metadata'),
                                       micro[1].split('.')[0].split('/')[-1] + '_train_movie.txt')

        movie = np.loadtxt(movie_file_name).view(complex)
        movie = movie.reshape(pc, zsize, accPix)


        ref_file_name = os.path.join(os.path.join(str(saved_folder), 'metadata'),
                                       micro[1].split('.')[0].split('/')[-1] + '_train_reference.txt')

        reference = np.loadtxt(ref_file_name).view(complex)
        reference = reference.reshape(pc, accPix)


        logfile = os.path.join(logfile_loc, movie_all[g].split('/')[-1].split('.tiff')[0] + '.star')
        shift_x, shift_y = givemotioncorrshifts(logfile)
        expKer = True
        A = np.zeros((pc, pc))
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

        basis = np.zeros((dc, pc))
        eigenVals = np.zeros(dc)

        for d in range(dc):
            l = np.sqrt(S[d])
            for p in range(pc):
                basis[d, p] = l * Vt[d, p]

        for d in range(dc):
            eigenVals[d] = S[d]

        initialshifts = np.zeros((pc, zsize, 2))
        for i in range(pc):
            initialshifts[i, :, :] = np.array([np.array(shift_x)[:] * -1, np.array(shift_y)[:] * -1]).swapaxes(0, 1)
        pad = 0
        newposi = np.zeros((pc, zsize + pad, 2))
        partCt = 2 * (pc + pc * (zsize - 1))
        initialCoeff = np.zeros(partCt)
        initialCoeff = posToParams(initialshifts, initialCoeff, pc, zsize, basis, eigenVals, dc)
        options_min = {}
        options_min['disp'] = True
        options_min['maxcor'] = 200
        options_min['ftol'] = 0.0001
        options_min['gtol'] = 0.90
        options_min['eps'] = 1e-05
        options_min['iprint'] = 80
        options_min['maxfun'] = 10000
        options_min['maxiter'] = 10000
        options_min['maxls'] = 40
        lb = 9.99e-21
        ub = 10e19
        bou = (lb, ub)
        bb = np.tile(bou, partCt).reshape(partCt, 2).tolist()

        result = minimize(fun=gpmotion_f, x0=initialCoeff,
                          args=(newposi, pc, zsize, basis, correlate, sig_a_vals_px, eigenVals, dc),
                          method='L-BFGS-B', jac=gpmotion_grad, tol=9.99e-17, options=options_min)
        tracks = np.zeros((pc, zsize, 2))
        tracks = ParamsTopos(result.x, tracks, pc, zsize, basis, dc)

        tsc = updateTsc(tracks, pc, zsize, boxsize, accPix, accCords, movie, reference, damage)
        tscAs[0] += tsc[0]
        tscAs[1] += tsc[1]
        tscAs[2] += tsc[2]


        del movie
        del reference
        del correlate


        del part_cord
        del tracks
        del result
        del options_min
        del shift_x
        del shift_y
        del A
        del U
        del S
        del Vt
        del eigenVals
        del basis
        del initialshifts
        del initialCoeff
        del newposi
        del partCt
        del bb

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    wg = tscAs[1] * tscAs[2]
    if wg > 0.0:
        tp = tscAs[0] / np.sqrt(wg)

    return -tp, tscAs


def motiontoProblem(vd):
    s = vd[2]
    accThresh = 500.0
    accEps = 1e-5
    velscale = 10000.0
    divScale = 1.0
    accScale = 1000.0
    if s < accThresh:
        if s > accEps:
            w = 1.0 / s
        else:
            w = 1.0 / accEps
    else:
        w = 0.0

    return [vd[0] * velscale, vd[1] * divScale, vd[2] * accScale]


def problemtoMotion(x):
    velscale = 10000.0
    divScale = 1.0
    accScale = 1000.0
    accThresh = 500.0
    accEps = 1e-5
    w = x[2] / accScale

    if x[2] > accEps:
        newx2 = x[2] / accScale
    else:
        newx2 = accEps / accScale

    return [x[0] / velscale, x[1] / divScale, newx2]


def Butterworthfilter(img, radIn, radOut):
    w = img.shape[0]
    h = img.shape[1]
    out = np.zeros((w, h))

    for y in range(h):
        for x in range(w):
            xx = x
            if y <= h / 2:
                yy = y
            else:
                yy = y - h

            r = np.sqrt(xx * xx + yy * yy)
            if r < radIn:
                out[x, y] = img[x, y]
            elif r < radOut:
                t = (r - radIn) / (radOut - radIn)
                a = 0.5 * (1.0 + np.cos(np.pi * t))
                out[x, y] = a * img[x, y]
            else:
                out[x, y] = 0.0
    return out


def prepAlignment(globaldose, angout_pix, boxsize, zsize, k_cutoff, k0, k1, mc, sel_files,
                  gainrefname, movie_names, maxRange, volft_1, volft_2):
    alignDoseWeights = np.zeros((zsize, boxsize // 2 + 1, boxsize))
    for i in range(zsize):
        dose = globaldose * (i)
        damageweig = damageWeights(dose, boxsize, angout_pix)
        alignDoseWeights[i] = Butterworthfilter(damageweig, k_cutoff - 1, k_cutoff + 1)

    accCords = []
    num = 0

    for y in range(boxsize):
        for x in range(boxsize // 2 + 1):
            xx = x
            if y < boxsize // 2 + 1:
                yy = y
            else:
                yy = y - boxsize

            r = round(np.sqrt(xx * xx + yy * yy))
            if r >= k0 and r < k1:
                accCords.append([x, y])
                num += 1

    accPix = num
    align_newgh = np.zeros((zsize, accPix))

    for f in range(zsize):
        dose = globaldose * (f)
        damagewe = damageWeights(dose, boxsize, angout_pix)
        for i in range(accPix):
            align_newgh[f, i] = damagewe[accCords[i][0], accCords[i][1]]

    pctot = 0
    maxRangeP = 2 * maxRange
    current_frame = EMData()
    gainref = EMData()
    gainref.read_image(gainrefname, 0)
    Util.mul_scalar(gainref, -1.0)

    # CC_crop = []
    # refer_all = []
    # movies_all = []

    kb_1 = 0
    kb_2 = 0

    saved_folder = '/home/adnan/PycharmProjects/DoseWeighting/Polis001/'


    ima_start, ima_end = sp_applications.MPI_start_end(mc, n_mpi_procs, my_mpi_proc_id)

    for micro in enumerate(sel_files[ima_start:ima_end]):
    # for g in tqdm(range(mc), desc="Pre_Alignment"):
        g = sel_files.index(micro[1])

        print("\n Running pre-alignment in parallel", g, micro[1])
        pc = EMUtil.get_image_count(micro[1])
        if pc < 2:
            continue

        pctot += pc

        part_cord = EMUtil.get_all_attributes(micro[1], "ptcl_source_coord")
        ctf_params_part = EMUtil.get_all_attributes(micro[1], "ctf")
        chunk_arr = EMUtil.get_all_attributes(micro[1], "chunk_id")
        project_params = EMUtil.get_all_attributes(micro[1], "xform.projection")
        nx = EMUtil.get_all_attributes(micro[1], 'nx')
        nxx = nx[0]
        ny = EMUtil.get_all_attributes(micro[1], 'ny')
        nyy = ny[0]
        current_frame.read_image(movie_names[g], 0)
        xsize = current_frame.get_xsize()
        ysize = current_frame.get_ysize()
        zsize = EMUtil.get_image_count(movie_names[g])
        cen_xx = xsize // 2
        cen_yy = ysize // 2

        refer_1 = get_2D_project_for_all_ptcl_from_reference(volft_1, kb_1, project_params, zsize,
                                                             show=False)
        refer_2 = get_2D_project_for_all_ptcl_from_reference(volft_2, kb_2, project_params, zsize,
                                                             show=False)

        var_movie, cnt_movie = init_variance_array(movie_names[g], nxx, nyy, cen_xx, cen_yy, part_cord, zsize, gainref)
        sigma = movie_vari_powspect(movie_names[g], gainref, var_movie, cnt_movie, zsize, part_cord,
                                    nxx, nyy, cen_xx, cen_yy)
        CCFS = np.zeros((zsize, len(part_cord), maxRangeP, maxRangeP))
        movie_mic = np.zeros((zsize, len(part_cord), accPix), dtype=complex)
        CCFS, movie_mic = movieCC_train(movie_names[g], nxx, nyy, sigma, var_movie, cnt_movie, globaldose, angout_pix,
                                        zsize, gainref, part_cord, ctf_params_part, cen_xx, cen_yy,
                                        refer_1, refer_2, CCFS, movie_mic, chunk_arr, alignDoseWeights,
                                        accPix, accCords, maxRangeP)

        CCFS = CCFS.swapaxes(0, 1)
        movie_mic = movie_mic.swapaxes(0, 1)

        count_p = 0
        count_f = 0
        f = open(os.path.join(os.path.join(str(saved_folder), 'metadata'),
                               micro[1].split('.')[0].split('/')[-1] + '_train_CC.txt'), 'w')
        for part in CCFS:
            count_p += 1
            for frame in part:
                np.savetxt(f, frame, fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_{}_Frame_{}'.format(count_p, count_f))
                count_f += 1
        f.close()

        count_p = 0
        count_f = 0
        f = open(os.path.join(os.path.join(str(saved_folder), 'metadata'),
                              micro[1].split('.')[0].split('/')[-1] + '_train_movie.txt'), 'w')
        for part in movie_mic:
            count_p += 1
            for frame in part:
                np.savetxt(f, frame.view(float), fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_{}_Frame_{}'.format(count_p, count_f))
                count_f += 1
        f.close()



        del movie_mic
        del CCFS
        del var_movie
        del cnt_movie
        del sigma

        refer_mic = np.zeros((len(part_cord), accPix), dtype=complex)
        for p in range(len(part_cord)):
            if chunk_arr[p] == 1:
                refer = refer_1[p]
            else:
                refer = refer_2[p]
            refer = sp_filter.filt_ctf(refer, ctf_params_part[p], sign=1, binary=False)
            fft_ref_img = np.fft.fft2(refer.get_2dview())
            np.divide(fft_ref_img, float(nxx * nyy), out=fft_ref_img)
            fft_ref_img[0, 0] = 0.0
            refer_mic[p] = accelerate(fft_ref_img, accPix, accCords)

            del refer
            del fft_ref_img

        count_p = 0
        f = open(os.path.join(os.path.join(str(saved_folder), 'metadata'),
                              micro[1].split('.')[0].split('/')[-1] + '_train_reference.txt'), 'w')
        for part in refer_mic:
                np.savetxt(f, part.view(float), fmt='%20.6f', delimiter='  ', newline='\n', header='Particle_{}'.format(count_p))
                count_p += 1
        f.close()

        del refer_mic
        del refer_1
        del refer_2
        del ctf_params_part
        del part_cord
        del project_params
        del chunk_arr
        del nx
        del ny


    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # CC_crop_m = sp_utilities.wrap_mpi_gatherv(CC_crop, 0, mpi.MPI_COMM_WORLD)
    # movies_all_m = sp_utilities.wrap_mpi_gatherv(movies_all, 0, mpi.MPI_COMM_WORLD)
    # refer_all_m = sp_utilities.wrap_mpi_gatherv(refer_all, 0, mpi.MPI_COMM_WORLD)
    #
    # CC_crop = sp_utilities.bcast_list_to_all(CC_crop_m, my_mpi_proc_id, main_mpi_proc, mpi.MPI_COMM_WORLD)
    # movies_all = sp_utilities.bcast_list_to_all(movies_all_m, my_mpi_proc_id, main_mpi_proc, mpi.MPI_COMM_WORLD)
    # refer_all = sp_utilities.bcast_list_to_all(refer_all_m, my_mpi_proc_id, main_mpi_proc, mpi.MPI_COMM_WORLD)

    return align_newgh, accPix, accCords


def accelerate(in_img, accPix, accCords, comp=True):
    if comp == True:
        out_img = np.zeros(accPix, dtype=complex)
    else:
        out_img = np.zeros(accPix)
    for i in range(accPix):
        out_img[i] = in_img[accCords[i][1], accCords[i][0]]
    return out_img


def crop_corner(correl, w, h):
    # Todo: Here i have to write the function to crop fourier image
    Crop_images = np.zeros((w, h))
    w1 = correl.shape[0]
    h1 = correl.shape[1]
    for y in range(h1):
        for x in range(w1):
            if x < w1 / 2:
                x1 = x
            else:
                x1 = x - w1
            if y < h1 / 2:
                y1 = y
            else:
                y1 = y - h1
            if (x1 < w // 2 and y1 < h // 2 and x1 >= -w // 2 and y1 >= -h // 2):
                if x1 < 0:
                    x0 = x1 + w
                else:
                    x0 = x1
                if y1 < 0:
                    y0 = y1 + h
                else:
                    y0 = y1
                Crop_images[x0][y0] = correl[x][y]

    return Crop_images


movie_names = return_movie_names('/home/adnan/PycharmProjects/DoseWeighting/MOVIES_RELION/20170629_*_frameImage.tiff')
stackfilename = 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv6/MotionCorr/job030/MOVIES_RELION/Particles/sphire_relion_stack'
logfile_loc = '/home/adnan/PycharmProjects/DoseWeighting/Relion_log/'
boxsize = 360
zsize = 24
angout_pix = 0.885
mc = 24
globaldose = 1.277

reference_vol_0 = '/home/adnan/PycharmProjects/DoseWeighting/run_half1_class001_unfil.mrc'
reference_vol_1 = '/home/adnan/PycharmProjects/DoseWeighting/run_half2_class001_unfil.mrc'
ref_mask = '/home/adnan/PycharmProjects/DoseWeighting/mask.mrc'
gainrefname = '/home/adnan/PycharmProjects/DoseWeighting/MOVIES_RELION/gain.mrc'

ref_EM_1 = sp_utilities.get_im(reference_vol_0)
ref_EM_2 = sp_utilities.get_im(reference_vol_1)
mask_volume = sp_utilities.get_im(ref_mask)

ref_EM_1 *= mask_volume
ref_EM_2 *= mask_volume

volft_1 = sp_projection.prep_vol(ref_EM_1, npad=2, interpolation_method=1)
volft_2 = sp_projection.prep_vol(ref_EM_2, npad=2, interpolation_method=1)

del ref_EM_1
del ref_EM_2
del mask_volume


# prepAlignment(globaldose, angout_pix, boxsize, zsize)
MotionParamEstima(movie_names, stackfilename, boxsize, zsize, angout_pix,
                  mc, volft_1, volft_2, gainrefname, logfile_loc)

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
mpi.mpi_finalize()
"""


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
    help='reference volume 0'
)

parser.add_argument(
    'reference_vol_1',
    type=str,
    help='reference volume 1'
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
Gainreffile = os.path.join(os.path.abspath(os.path.join(str(args.movies_folder), os.pardir)), 'gain.mrc')
movie_names = return_movie_names(str(args.movies_folder))

# -------------------------------Reading corrected sums log files for reading the global shifts applied
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
ptcl_source_images_xx = []

for i in range(no_of_imgs):
    ptcl_source_images_xx.append(str(os.path.basename(ptcl_source_images[i])).split('.mrc')[0])

count = 0
for i in range(no_of_micrographs):
    temp = np.where(np.array(ptcl_source_images_xx) == str(os.path.basename(movie_names[i])).split('.tiff')[0])
    if len(temp[0]) != 0:
        movie_name_x.append(movie_names[i])
        count += 1

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

"""




