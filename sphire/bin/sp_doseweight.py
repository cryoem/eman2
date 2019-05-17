from __future__ import print_function
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import cPickle as pickle
from EMAN2 import *
import EMAN2_cppwrap


import sp_utilities
import sp_projection
import sp_statistics
import sp_filter
import mpi
import sp_applications

import sp_fundamentals

from scipy import fftpack


location =os.getcwd()


RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ

no_of_micrographs = 112

main_mpi_proc = 0
if RUNNING_UNDER_MPI:
    pass  # IMPORTIMPORTIMPORT from mpi import mpi_init
    pass  # IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM

    mpi.mpi_init(0, [])
    my_mpi_proc_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    n_mpi_procs = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
else:
    my_mpi_proc_id = 0
    n_mpi_procs = 1


ima_start , ima_end = sp_applications.MPI_start_end(no_of_micrographs, n_mpi_procs, my_mpi_proc_id)


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

    print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
    return input_file_path_list


def return_images_from_movie(movie_name, show_first = False):
    mov_imgs = EMAN2_cppwrap.EMData.read_images(movie_name)
    if show_first:
        plt.ion()
        plt.figure()
        plt.imshow(mov_imgs[0].get_3dview()[0], cmap=plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()
    print("\n No of images in the selected movie %s are %d" % (os.path.basename(movie_name) ,len(mov_imgs[0].get_3dview())))
    print("Shape of the Movie i.e. dimension of the frames and no of frames", np.shape(mov_imgs[0].get_3dview()))
    print("No of images in the movie",mov_imgs[0].get_zsize() )
    print("X dimension of image in the movie", mov_imgs[0].get_xsize())
    print("X dimension of image in the movie", mov_imgs[0].get_ysize())
    return mov_imgs[0]



"""
Reading a particle stack to find all the parameters saved in the header file
"""

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

    return no_of_imgs_once, ptcl_source_images_once, project_params_all_once, particle_coordinates_all_once, ctf_params_all_once, nx_all_once, ny_all_once, nz_all_once


def find_particles_info_from_movie(stack, movie_name, no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, ctf_params_all, nx_all, ny_all, nz_all, show_first = False):

    #-----------------------------------------------   CTF related attributes
    """
    defocus = defocus associated with the image, positive value corresponds to underfocus
    cs =  spherical aberration constant [mm].
    voltage = accelerating voltage of the microscope [kV]
    apix =
    bfactor = The parameter in Gaussian like envelope function, which roughly explains Fourier factor dumping of the image.
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

    for i in range(no_of_imgs):
        if (str(os.path.basename(movie_name)) == str(os.path.basename(ptcl_source_images[i]))):
            project_params_per_movie.append(project_params_all[i])
            particle_coordinates_per_movie.append(particle_coordinates_all[i])
            ctf_params_per_movie.append(ctf_params_all[i])
            nx_per_movie.append(nx_all[i])
            ny_per_movie.append(ny_all[i])
            nz_per_movie.append(nz_all[i])

    print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)), len(project_params_per_movie)))
    print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
    print("Projection parameters for 1st particle in the stack are ", project_params_per_movie[0].get_params('spider'))
    print("Dimensions x for all particles are ", nx_per_movie)
    print("Dimensions y for all particles are ", ny_per_movie)
    print("Dimensions z for all particles are ", nz_per_movie)


    if show_first:
        ima = EMAN2_cppwrap.EMData()
        ima.read_image(stack, 0, False)
        plt.ion()
        plt.figure()
        plt.imshow(ima.get_2dview(), cmap=plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()

    return project_params_per_movie, particle_coordinates_per_movie, ctf_params_per_movie, nx_per_movie, ny_per_movie, nz_per_movie

"""
Reading a reference map
"""
def get_2D_project_for_all_ptcl_from_reference(volume_ft , project_params_in_stack, show = False):
    project_2D_per_movie = []
    for i in range(len(project_params_in_stack)):
        params_substack = project_params_in_stack[i].get_params('spider')
        params_for_each_image = [params_substack['phi'], params_substack['theta'], params_substack['psi'],
                                 params_substack['tx'], params_substack['ty']]
        project_2D_per_movie.append(sp_projection.prgl(volume_ft, params_for_each_image, interpolation_method=1))
    if show:
        plt.ion()
        plt.figure()
        plt.imshow(project_2D_per_movie[0].get_2dview(), cmap = plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()
    return project_2D_per_movie


"""
Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
creating a window around to extract the same particle from each frame
"""
#----------------------- Particle cordinate
def get_all_reduce_ptcl_imgs(frames_i, maski, nxx, nyy, part_cord, ctf_para, cen_xx, cen_yy, cen_zz):
    particle_imgs_in_movie = []
    for j in range(len(part_cord)):
        crop_imgs = []
        for i in range(2*cen_zz):
            crop_imgs.append(Util.window(frames_i, nxx, nyy, 1, part_cord[j][0] - cen_xx,
                                         part_cord[j][1] - cen_yy, i - cen_zz))
            # st = Util.infomask(crop_imgs[i], None, False)  # find statistics of the image returns avg, std, min, max
            # Util.mul_scalar(crop_imgs[i], -1.0)
            # crop_imgs[i] += 2 * st[0]  # st[0]-> avergage  , shifting the negative range to positive
            #
            # st = Util.infomask(crop_imgs[i], maski, False)
            # crop_imgs[i] -= st[0]  # st[0]-> avergage
            # crop_imgs[i] /= st[1]  # st[1]-> standard deviation
            #
            # st = Util.infomask(crop_imgs[i], maski, False)
            # crop_imgs[i] -= st[0]  # st[0]-> avergage
            # crop_imgs[i] = sp_filter.filt_ctf(crop_imgs[i], ctf_para[j], binary=False)
            #
            # st = Util.infomask(crop_imgs[i], maski, False)
            # crop_imgs[i] -= st[0]  # st[0]-> avergage
            # crop_imgs[i] = sp_filter.filt_ctf(crop_imgs[i], ctf_para[j], binary=True)
            #
            # st = Util.infomask(crop_imgs[i], maski, False)
            # crop_imgs[i] -= st[0]  # st[0]-> avergage
            # st = Util.infomask(crop_imgs[i], maski, False)
            # crop_imgs[i] /= st[1]  # st[1]-> standard deviation
        particle_imgs_in_movie.append(crop_imgs)
    return particle_imgs_in_movie


def moving_avg_filter(fsc_curve):
    for i in range(2, len(fsc_curve) -3 ) :
        fsc_curve[i] = (fsc_curve[i] + fsc_curve[i-1] + fsc_curve[i-2] + fsc_curve[i+1] + fsc_curve[i+2])/5
    return fsc_curve

def smooth(x,window_len):
    import numpy as np
    # s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    w = np.ones(window_len, 'd')
    y = np.convolve(w / w.sum(), x, mode='same')
    return y

from scipy.optimize import curve_fit
def fitthecurve(data, inital_guess):
    try:
        coeff, var_matrix = curve_fit(fitfunc, data[0], data[1], p0 = inital_guess)
    except RuntimeError as err:
        print('Handling run-time error:', err)
        coeff = [0, 0, 0]
    print(coeff[0], coeff[1], coeff[2])

    expo = fitfunc(data[0], *coeff)
    return [data[0], expo, coeff, var_matrix]


def fitfunc(x, a, b,c,d ):
    return -a * np.exp(c + (4*b * (np.array(x)*np.array(x)))) + d


def fitslopefunc(x, m,c):
    return m*x + c
#%%

"""
Reading of Movie frames in .mrc format and display one frame
"""
ABSOLUTE_PATH_TO_MRC_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/MOVIES/"
input_image_path = os.path.join(ABSOLUTE_PATH_TO_MRC_FOLDER, "TcdA1-*_frames.mrc")

movie_names = return_movie_names(input_image_path)
# print("\n Following movies were detected : \n",movie_names)
for micro in enumerate(movie_names[ima_start:ima_end]):
    print(micro)

fsc_values = []
fsc_avgs = []
frequencies = []
fsc_s_all = []


for micro in enumerate(movie_names[ima_start:ima_end]):
    # micro = (0, '/home/adnan/PycharmProjects/DoseWeighting/MOVIES/TcdA1-0013_frames.mrc')
    frames = return_images_from_movie(micro[1], show_first = False)

    """
    Reading a particle stack to find all the parameters saved in the header file
    """

    stackfilename = "bdb:/home/adnan/PycharmProjects/DoseWeighting/Substack/isac_substack"
    no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, ctf_params_all, nx_all, ny_all, nz_all = read_all_attributes_from_stack(stackfilename)
    project_params, particle_coordinates, ctf_params, nx, ny, nz = find_particles_info_from_movie(stackfilename, micro[1], \
                                                                                      no_of_imgs,ptcl_source_images,project_params_all,particle_coordinates_all,ctf_params_all, nx_all, ny_all, nz_all,show_first=False)

    """
    Reading a reference map
    """

    #-------------- Loading the reference volume
    ref_vol_filename = "/home/adnan/PycharmProjects/DoseWeighting/vol_combined.hdf"
    ref_volume = sp_utilities.get_im(ref_vol_filename)
    #---------------Preparing the volume for calculation of projections
    volft  =  sp_projection.prep_vol(ref_volume, npad = 2, interpolation_method=1) #transforms the volume in fourier space and then expands the volume with npad and centers the volume and returns the volume


    ref_project_2D_ptcl_all = get_2D_project_for_all_ptcl_from_reference(volft, project_params, show = False) #Projection of 3D volume in 2-D for all the particles in all frames in one movie
    mask = sp_utilities.model_circle(nx[0] / 2, nx[0], nx[0])  # nx/2 is for the radius
    for i in range(len(ref_project_2D_ptcl_all)):
        # st = Util.infomask(ref_project_2D_ptcl_all[i], mask, False)
        # ref_project_2D_ptcl_all[i]-= st[0]  # st[0]-> average
        # ref_project_2D_ptcl_all[i] /= st[1]  # st[1]-> standard deviation

        # st = Util.infomask(ref_project_2D_ptcl_all[i], mask, False)
        # ref_project_2D_ptcl_all[i] -= st[0]  # st[0]-> avergage
        # ref_project_2D_ptcl_all[i] = sp_filter.filt_ctf(ref_project_2D_ptcl_all[i], ctf_params[i], binary=False)

        # st = Util.infomask(ref_project_2D_ptcl_all[i], mask, False)
        # ref_project_2D_ptcl_all[i] -= st[0]  # st[0]-> avergage
        ref_project_2D_ptcl_all[i] = sp_filter.filt_ctf(ref_project_2D_ptcl_all[i], ctf_params[i], binary=True)

    """
    Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
    creating a window around to extract the same particle from each frame
    """
    cen_x = frames.get_xsize() // 2
    cen_y = frames.get_ysize() // 2
    cen_z = frames.get_zsize() // 2

    particle_imgs= get_all_reduce_ptcl_imgs(frames, mask, nx[0], ny[0], particle_coordinates, ctf_params, cen_x, cen_y, cen_z)

    zsize = frames.get_zsize()
    del frames
    """
    Calculating the fourier shell correlation of all the particle images with respect to 2-D reference projection of 3-D volume
    """
    fsc_s = []
    for i in range (zsize):
        fsc_frames = []
        for j in range(len(particle_imgs)):
            fsc_frames.append(sp_statistics.fsc(particle_imgs[j][i], ref_project_2D_ptcl_all[j]))
        fsc_s.append(fsc_frames)

    fsc_final = []
    for i in range (zsize):
        fsc_sum = [entry/len(fsc_s[i]) for entry in fsc_s[i][0][1]]
        for fsc in fsc_s[i][1:]:  # one frame ahead for averageing
            for idx, ok in enumerate(fsc[1]):
                fsc_sum[idx] += ok/len(fsc_s[i])
        fsc_final.append(fsc_sum)

    fsc_final_avg =  []
    for idx in range (0,len(fsc_final)-3):
        avv = []
        for p in range(len(fsc_final[idx])):
            avv.append((fsc_final[idx][p] + fsc_final[idx+1][p] + fsc_final[idx+2][p] + fsc_final[idx+3][p] )/4)
            # print(idx, p)
        fsc_final_avg.append(avv)

    for idx in range(len(fsc_final) - 3, len(fsc_final)):
        avv = []
        for p in range(len(fsc_final[idx])):
            avv.append((fsc_final[idx][p] + fsc_final[idx-1][p] + fsc_final[idx-2][p] + fsc_final[idx-3][p] )/4)
            # print(idx, p)
        fsc_final_avg.append(avv)

    fsc_values.append(fsc_final)
    fsc_avgs.append(fsc_final_avg)
    frequencies.append(fsc_s[0][0][0])
    # fsc_s_all.append(fsc_s)





fsc_values_per_micrograph = sp_utilities.wrap_mpi_gatherv(fsc_values, 0, mpi.MPI_COMM_WORLD)
fsc_avgs_per_micrograph = sp_utilities.wrap_mpi_gatherv(fsc_avgs, 0, mpi.MPI_COMM_WORLD)
freq_per_micrograph = sp_utilities.wrap_mpi_gatherv(frequencies, 0, mpi.MPI_COMM_WORLD)
# raw_fsc_values =  sp_utilities.wrap_mpi_gatherv(fsc_s_all, 0 , mpi.MPI_COMM_WORLD)

print(np.array(fsc_values_per_micrograph).shape)
print(np.array(fsc_avgs_per_micrograph).shape)
print(np.array(freq_per_micrograph).shape)
"""
Writing data in pickle files

"""
mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

if my_mpi_proc_id == main_mpi_proc:
    with open(str(no_of_micrographs) + '_Micrograph_plot_values', 'wb') as wb:
        for value in fsc_values_per_micrograph:
            pickle.dump(value, wb)
    wb.close()


    with open(str(no_of_micrographs) + '_Micrograph_avg_plot_values', 'wb') as wb:
        for value in fsc_avgs_per_micrograph:
            pickle.dump(value, wb)
    wb.close()


    with open(str(no_of_micrographs) + '_Micrograph_frequencies_plot_values', 'wb') as wb:
        for value in freq_per_micrograph:
            pickle.dump(value, wb)
    wb.close()


    # with open(str(no_of_micrographs) + '_Micrograph_raw_fsc_values', 'wb') as wb:
    #     for value in raw_fsc_values:
    #         pickle.dump(value, wb)
    # wb.close()


print("I am finish")

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
mpi.mpi_finalize()

#%%

# """

fsc_values_per_micrograph = []
fsc_avgs_per_micrograph = []
freq_per_micrograph = []


fsc_all_data_per_micrograph = []

with open(location + '/sphire/bin/' + str(no_of_micrographs) + '_Micrograph_plot_values', 'rb') as rb:
    while True:
        try:
            fsc_values_per_micrograph.append(pickle.load(rb))
        except EOFError:
            break

rb.close()


with open(location + '/sphire/bin/' + str(no_of_micrographs) + '_Micrograph_avg_plot_values', 'rb') as rb:
    while True:
        try:
            fsc_avgs_per_micrograph.append(pickle.load(rb))
        except EOFError:
            break
rb.close()


with open(location + '/sphire/bin/' + str(no_of_micrographs) + '_Micrograph_frequencies_plot_values', 'rb') as rb:
    while True:
        try:
            freq_per_micrograph.append(pickle.load(rb))
        except EOFError:
            break
rb.close()


# with open(location + '/sphire/bin/' + str(no_of_micrographs) + '_Micrograph_raw_fsc_values', 'rb') as rb:
#     while True:
#         try:
#             fsc_all_data_per_micrograph.append(pickle.load(rb))
#         except EOFError:
#             break
# rb.close()






print(np.array(fsc_values_per_micrograph).shape)
print(np.array(fsc_avgs_per_micrograph).shape)
print(np.array(freq_per_micrograph).shape)
# print(np.array(fsc_all_data_per_micrograph).shape)

fsc_values_per_micrograph = np.array(fsc_values_per_micrograph)
fsc_avgs_per_micrograph = np.array(fsc_avgs_per_micrograph)
freq_per_micrograph = np.array(freq_per_micrograph)
# fsc_all_data_per_micrograph = np.array(fsc_all_data_per_micrograph)



fsc_sum_per_frame = []
for frames_ind in range(fsc_values_per_micrograph.shape[1]):
    fsc_sum =  [entry / len(fsc_values_per_micrograph) for entry in fsc_values_per_micrograph[0][frames_ind]]
    for micrograph in fsc_values_per_micrograph[1:]:
        for ind , values in enumerate(micrograph[frames_ind]):
            fsc_sum[ind] += values/len(fsc_values_per_micrograph)
    fsc_sum_per_frame.append(fsc_sum)



fsc_sum_avg_per_frame = []
for frames_ind in range(fsc_avgs_per_micrograph.shape[1]):
    fsc_sum =  [entry / len(fsc_avgs_per_micrograph) for entry in fsc_avgs_per_micrograph[0][frames_ind]]
    for micrograph in fsc_avgs_per_micrograph[1:]:
        for ind , values in enumerate(micrograph[frames_ind]):
            fsc_sum[ind] += values/len(fsc_avgs_per_micrograph)
    fsc_sum_avg_per_frame.append(fsc_sum)



# fsc_sum_per_frame = []
# for frames_ind in range(fsc_values_per_micrograph.shape[1]):
#     fsc_sum =  [entry / len(fsc_values_per_micrograph) for entry in fsc_values_per_micrograph[0][frames_ind]]
#     for micrograph in fsc_values_per_micrograph[1:]:
#         for ind in range( np.array(micrograph).shape[1]  ):
#             fsc_sum[ind] += micrograph[frames_ind][ind]/len(fsc_values_per_micrograph)
#     fsc_sum_per_frame.append(fsc_sum)




B_values = []
C_values = []
coeff_list = []
offset_start = 5
offset_end  = 110
window_len = 1



fsc_sum_per_frame = np.array(fsc_sum_per_frame)[:] * -1
fsc_sum_avg_per_frame = np.array(fsc_sum_avg_per_frame)[:] * -1


for i in range (len(fsc_sum_per_frame)):
    print(i)
    data_x = freq_per_micrograph[0][offset_start:offset_end]
    data_y = smooth(fsc_sum_per_frame[i],window_len)
    data_y = data_y[offset_start:offset_end]
    # data_y[data_y <= 0] = 0.0000000001
    coeff, var_matrix = curve_fit(fitfunc, data_x, data_y, p0 =[-0.15, -75, 0.25, 0.035]  ) #bounds=[[-0.20,  -300  ,-400 ], [-0.05, 0 , 400]]
    B_values.append(coeff[1])
    C_values.append(coeff[2])
    coeff_list.append(coeff)




i = 14
fitcurv = fitfunc(freq_per_micrograph[0][offset_start:offset_end], *coeff_list[i])
plt.figure()
plt.plot(freq_per_micrograph[0][offset_start:offset_end], smooth(fsc_sum_per_frame[i][offset_start:offset_end], window_len), label= 'Orig_' + str(i))
plt.plot(freq_per_micrograph[0][offset_start:offset_end], fitcurv, label='fit_' + str(i))
plt.legend()
plt.show()



i = 0

my_mpi_proc_id = 0
main_mpi_proc = 0

if my_mpi_proc_id == main_mpi_proc:
    fig , ax = plt.subplots(nrows = 4, ncols=6 )
    for row in ax:
        for col in row:
            if i < len(fsc_sum_avg_per_frame):
                print(i)
                fitcurv = fitfunc(freq_per_micrograph[0][offset_start:offset_end], *coeff_list[i])
                col.plot(freq_per_micrograph[0][offset_start:offset_end], smooth(fsc_sum_per_frame[i][offset_start:offset_end], window_len), label= 'Orig_' + str(i))
                # col.plot(freq_per_micrograph[0][offset_start:offset_end], smooth(fsc_sum_avg_per_frame[i][offset_start:offset_end], window_len), label= 'Avg_' + str(i))
                col.plot(freq_per_micrograph[0][offset_start:offset_end], fitcurv, label='fit_' + str(i))
                i+= 1
                col.legend()
                # col.set_yscale('log')
    # fig.savefig(str(no_of_micrographs) + "_Micrograph_averaged" + ".pdf", dpi=1000)
    plt.show()




    i = 0
    offset_end = 176
    fig , ax = plt.subplots(nrows = 4, ncols=6 )
    for row in ax:
        for col in row:
            if i < len(fsc_sum_avg_per_frame):
                print(i)
                data_x = np.array(freq_per_micrograph[0][offset_start:offset_end])
                data_y = smooth(fsc_sum_per_frame[i][offset_start:offset_end], window_len)
                data_y[data_y <= 0] = data_y[data_y <= 0] * -1
                data_y = np.log(data_y)
                data_x = data_x* data_x
                # coeff, var_matrix = curve_fit(fitfunc(), data_x[offset_start:offset_end], data_y[offset_start:offset_end], p0 = coeff_list[i])
                coeff, var_matrix = curve_fit(fitslopefunc, data_x[offset_start:offset_end], data_y[offset_start:offset_end])
                col.plot(data_x[offset_start:offset_end] ,  data_y[offset_start:offset_end], label= 'Orig_' + str(i))
                col.plot(data_x[offset_start:offset_end] ,  fitslopefunc(data_x[offset_start:offset_end], *coeff), label='fit_' + str(i))
                i+= 1
                col.legend()
    plt.show()



"""

frames_range = 24
freq_range = len(freq_per_micrograph[0][offset_start:offset_end])
iterations = 5


c_list = [np.random.random() for i in range(frames_range)]
b_list = [np.random.random() for i in range(frames_range)]
d_list = [np.random.random() for i in range(freq_range)]


for i in range(iterations):
    for f in range (frames_range):
        data_x = np.array(freq_per_micrograph[0][offset_start:offset_end])
        data_y = smooth(fsc_sum_per_frame[f][offset_start:offset_end], window_len)
        data_y[data_y <= 0] = data_y[data_y <= 0] * -1
        data_y = np.log(data_y)
        data_x = data_x * data_x

        for k in range (freq_range):
            estim_d = lambda x, param_d : np.log(param_d) + c_list[f] + 4 * b_list[f] * x
            coeff , covar = curve_fit(estim_d,  data_x , data_y )
            d_list[k] = coeff[0]

        print("f", f)
        estim_bc = lambda x, param_c, param_b :  np.log(d_list) + param_c + 4 * param_b * x
        coeff, covar = curve_fit(estim_bc,  data_x , data_y)
        c_list[f] = coeff[0]
        b_list[f] = coeff[1]




plt.figure()
for f in range (2):
    data_xx = np.array(freq_per_micrograph[0][offset_start:offset_end])
    data_yy = smooth(fsc_sum_per_frame[f][offset_start:offset_end], window_len)
    data_yy[data_yy <= 0] = data_yy[data_yy <= 0] * -1
    data_yy = np.log(data_yy)
    data_xx = data_xx * data_xx


    estim_d = lambda x: np.log(d_list[0]) + c_list[f] + 4 * b_list[f] * x
    plt.plot(data_xx, estim_d(data_xx), label='1_{0}'.format(f))
    plt.plot(data_xx, data_yy, label='2_{0}'.format(f))

plt.legend()


"""

"""

frames_range = 24
freq_range = len(freq_per_micrograph[0][offset_start:offset_end])
iterations = 5


c_list = [np.random.random() for i in range(frames_range)]
b_list = [np.random.random() for i in range(frames_range)]
d_list = [np.random.random() for i in range(freq_range)]


for i in range(iterations):
    for f in range (frames_range):
        data_x = np.array(freq_per_micrograph[0][offset_start:offset_end])
        data_y = smooth(fsc_sum_per_frame[f][offset_start:offset_end], window_len)
        data_y[data_y <= 0] = data_y[data_y <= 0] * -1
        # data_y = np.log(data_y)
        # data_x = data_x * data_x

        for k in range (freq_range):
            estim_d = lambda x, param_d :  param_d * np.exp(c_list[f] + 4 * b_list[f] * np.array(x)*np.array(x))
            coeff , covar = curve_fit(estim_d,  data_x , data_y )
            d_list[k] = coeff[0]

        print("f", f)
        estim_bc = lambda x, param_c, param_b :  d_list * np.exp(param_c + 4 * param_b * np.array(x)*np.array(x))
        coeff, covar = curve_fit(estim_bc,  data_x , data_y)
        c_list[f] = coeff[0]
        b_list[f] = coeff[1]


plt.figure()
for f in range (2):
    data_xx = np.array(freq_per_micrograph[0][offset_start:offset_end])
    data_yy = smooth(fsc_sum_per_frame[f][offset_start:offset_end], window_len)
    data_yy[data_yy <= 0] = data_yy[data_yy <= 0] * -1
    # data_yy = np.log(data_yy)
    # data_xx = data_xx * data_xx


    estim_d = lambda x: d_list[0] * np.exp(c_list[f] + 4 * b_list[f] * np.array(x) * np.array(x))
    plt.plot(data_xx, estim_d(data_xx), label='1_{0}'.format(f))
    plt.plot(data_xx, data_yy, label='2_{0}'.format(f))

plt.legend()

"""


#%%
offset_start = 0
offset_end = -1
frames_range = 24
freq_range = len(freq_per_micrograph[0][offset_start:offset_end])
N_ITER =5



c_list = np.array([np.random.random() for i in range(frames_range)])
b_list = np.array([np.random.random() for i in range(frames_range)])
d_list = np.array([np.random.random() for i in range(freq_range)])

FCC_FITTED = np.zeros((frames_range,freq_range))

myk = []
for k_abs in range(freq_range):
        k = k_abs*1.0/freq_range
        myk.append(k)


for iteration in range(N_ITER):
    for k_index, k in enumerate(myk):
        fcc_per_k = fsc_sum_per_frame[:,k_index]
        f_d = lambda u,d: d * np.exp(c_list[u] + 4*b_list[u]*k**2)
        popt, pconv = curve_fit(f_d, np.arange(frames_range).tolist(), fcc_per_k)
        d_list[k_index] = popt[0]

    for f in range(frames_range):
        fcc_per_f =   fsc_sum_per_frame[f,:]
        f_c_b = lambda u,c,b: d_list[u] * np.exp(c + 4*b*(np.array(u)*1.0/freq_range)**2)
        fcc_per_f[fcc_per_f <= 0 ] = 0.01
        popt, pconv = curve_fit(f_c_b, range(freq_range), fcc_per_f[offset_start:offset_end], bounds = ((-np.inf,-np.inf), (np.inf,-50) ) )
        c_list[f] = popt[0]
        b_list[f] = popt[1]



for f in range(frames_range):
    for k_abs in range(freq_range):
        k = k_abs*1.0/freq_range
        FCC_FITTED[f,k_abs] = d_list[k_abs] * np.exp(c_list[f] + 4*b_list[f]*k**2)


i= 0
fig , ax = plt.subplots(nrows = 4, ncols=6 )
for row in ax:
    for col in row:
        if i < len(fsc_sum_avg_per_frame):
            col.plot(freq_per_micrograph[0][offset_start:offset_end], smooth(fsc_sum_per_frame[i][offset_start:offset_end], window_len), label= 'Orig_' + str(i))
            col.plot(freq_per_micrograph[0][offset_start:offset_end],  FCC_FITTED[i], label='fit_' + str(i))
            i+= 1
            col.legend()
plt.show()


freq_k = freq_per_micrograph[0].tolist()
freq_k.extend((np.arange(0.50 +  0.0111473913, 0.70 + 0.0111473913, 0.0111473913)))

freq_k = np.array(freq_k)

sum =[]
weight_per_frame = []
sum_k = np.zeros(np.shape(freq_k))
for j in range(frames_range):
    sum_k += np.exp(np.array(c_list)[j] + 4 * np.array(b_list)[j] * freq_k ** 2)

for i in range (frames_range):
    weight_per_frame.append(np.divide(np.exp(np.array(c_list)[i] + 4 * np.array(b_list)[i] * freq_k**2) , np.array(sum_k)))


plt.figure()
plt.plot(freq_k, weight_per_frame[0])
plt.xlabel('frequencies')
plt.ylabel('weights')


plt.figure()
plt.imshow(weight_per_frame[::-1])
plt.colorbar()


plt.figure()
plt.plot(np.arange(frames_range),b_list, 'o-',label = 'B-factors')
plt.legend()

plt.figure()
plt.plot(np.arange(frames_range),c_list, 'o-', label = 'C_values')
plt.legend()

plt.figure()
plt.plot(np.arange(freq_range), d_list, 'o-', label = 'D_values')
plt.legend()






#%%

plt.figure()
plt.plot(np.arange(24), B_values, 'o-')
plt.legend()
plt.xlabel('Frames', fontsize = 24)
plt.ylabel('B-Factor', fontsize = 24)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.show()


# i = 0
# fig , ax = plt.subplots(nrows = 4, ncols=6 )
# for row in ax:
#     for col in row:
#         if i < len(fsc_sum_avg_per_frame):
#             print("Hello",i)
#             fsc_final_i = smooth(fsc_sum_per_frame[i], window_len)
#             fsc_final_avg_i = smooth(fsc_sum_avg_per_frame[i], window_len)
#             rel_damp_sqrt = np.sqrt(
#                 np.divide(fsc_final_i - fsc_final_i * fsc_final_avg_i, fsc_final_avg_i - fsc_final_i * fsc_final_avg_i))
#             reldamp_i = np.log(rel_damp_sqrt)
#
#             col.plot(freq_per_micrograph[0][0:len(reldamp_i)], reldamp_i, label='Frame' + str(i))
#             i+= 1
#             col.legend()
#
# plt.show()
# """



# freq_k = freq_per_micrograph[0].tolist()
# freq_k.extend((np.arange(0.50 +  0.0111473913, 0.70 + 0.0111473913, 0.0111473913)))
#

# freq_k.extend((np.arange(0.50 +  0.00284091, 0.70 +  0.00284091, 0.00284091)))

# freq_k = np.array(freq_k)
# sum =[]
# weight_per_frame = []
# for i in range (len(B_values)):
#     sum_k = np.zeros(np.shape(freq_k))
#     for j in range( len(B_values)):
#         if j==i:
#             pass
#         else:
#             sum_k += np.exp(C_values[j] + 4 * B_values[j] * freq_k * freq_k)
#
#     weight_per_frame.append(np.divide(np.exp(C_values[i] + 4 * B_values[i] * freq_k * freq_k) , sum_k))



# plt.figure()
# plt.plot(freq_k, weight_per_frame[0])
# plt.xlabel('frequencies')
# plt.ylabel('weights')




#%%

def zero_pad(img, size_new):
    pad_extends = []
    dif_shape_y = size_new[0] - img.shape[0]
    dif_shape_x = size_new[1] - img.shape[1]

    pad_extends.append((dif_shape_y // 2, dif_shape_y // 2 + dif_shape_y % 2))
    pad_extends.append((dif_shape_x // 2, dif_shape_x // 2 + dif_shape_x % 2))

    padded = np.pad(img, pad_extends, "symmetric")

    return padded, pad_extends


def next_power_of2(number):
    return int(np.power(2, np.ceil(np.log2(number))))



micro = (0, '/home/adnan/PycharmProjects/DoseWeighting/MOVIES/TcdA1-0013_frames.mrc')
frames = return_images_from_movie(micro[1], show_first = False)

cen_x = frames.get_xsize() // 2
cen_y = frames.get_ysize() // 2
cen_z = frames.get_zsize() // 2


stackfilename = "bdb:/home/adnan/PycharmProjects/DoseWeighting/Substack/isac_substack"
no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, ctf_params_all, nx_all, ny_all, nz_all = read_all_attributes_from_stack(stackfilename)
project_params, particle_coordinates, ctf_params, nx, ny, nz = find_particles_info_from_movie(stackfilename, micro[1], \
                                                                                              no_of_imgs,
                                                                                              ptcl_source_images,
                                                                                              project_params_all,
                                                                                              particle_coordinates_all,
                                                                                              ctf_params_all, nx_all,
                                                                                              ny_all, nz_all,
                                                                                              show_first=False)

mask = sp_utilities.model_circle(nx[0] / 2, nx[0], nx[0])

particle_imgs = get_all_reduce_ptcl_imgs(frames, mask, nx[0], ny[0], particle_coordinates, ctf_params, cen_x, cen_y, cen_z)

particle_imgs = np.array(particle_imgs).swapaxes(0,1)



# four_img = sp_fundamentals.fft(particle_imgs[0][0], npad = 1)
#
#
# plt.figure()
# plt.imshow(four_img.real().get_2dview() ,cmap = plt.get_cmap('Greys'))
# plt.colorbar()

# four_img_new_1 = sp_fundamentals.prepf(particle_imgs[0][0], npad = 1)
#
#
# plt.figure()
# plt.imshow(four_img_new_1.real().get_2dview() ,cmap = plt.get_cmap('Greys'))
# plt.colorbar()


# four_img_new_2 = scipy.fft(particle_imgs[0][0].get_2dview())
#
#
#
# plt.figure()
# plt.imshow(four_img_new_2.real ,cmap = plt.get_cmap('Greys'))
# plt.colorbar()

img = particle_imgs[0][0].get_2dview()

# next_shape_y = next_power_of2(img.shape[0])
# next_shape_x = next_power_of2(img.shape[1])
# new_image_size = (next_shape_y, next_shape_x)
# padded, pad_extends = zero_pad(img, new_image_size)


four_img_new_3 = np.fft.fft2(img)

four_img_new_3  = np.fft.fftshift(four_img_new_3)


plt.figure()
plt.imshow(four_img_new_3.real ,cmap = plt.get_cmap('Greys'))
plt.colorbar()
plt.clim(-4000,4000)


# mask = np.arange(352* 352).reshape((352,352))
#
# indices = np.array(zip(*np.where(mask==mask)))


mask = np.zeros((352,352))

row, col = np.meshgrid(range(352), range(352), indexing = 'ij')
for i in range(np.shape(row)[0]):
    for j in range(np.shape(col)[0]):
        mask[i][j] = np.sqrt((row[i][j] - np.shape(row)[0]/2)**2   + (col[i][j] - np.shape(col)[0]/2)**2 )


# mask_norm = (mask / np.max(mask)) * 0.5


normal_to = mask[ mask.shape[0]/2][mask.shape[1]-1]


mask_norm = (mask / normal_to) * 0.5


plt.figure()
plt.imshow(mask_norm)
plt.colorbar()


index_frame = 0
mask_applied = np.zeros((352,352))
for i in range (mask.shape[0]):
    for j in range(mask.shape[1]):
        near_value = (np.abs(weight_per_frame[index_frame] - mask_norm[i][j])).argmin()
        mask_applied[i][j] = weight_per_frame[index_frame][near_value]


plt.figure()
plt.imshow(mask_applied,cmap = plt.get_cmap('Greys'))
plt.colorbar()

#%%


# i = 0
# window_len = 6
# fig , ax = plt.subplots(nrows = 4, ncols=6 )
# for row in ax:
#     for col in row:
#         if i < len(fsc_s):
#             print(i)
#             # fitcurv = fitfunc(fsc_s[i][0][0][offset_start:offset_end], *coeff_list[i])
#             col.plot(fsc_s[0][0][0],smooth(fsc_final[i], window_len), label='original' + str(i))
#             col.plot(fsc_s[0][0][0],smooth(fsc_final_avg[i], window_len), label='Average' + str(i))
#             i+= 1
#             col.legend()
#
# plt.show()


"""
#%%
plt.ioff()
plt.figure()
for i in range (zsize):
    plt.plot(fsc_s[i][0][0], smooth(fsc_final[i], window_len), label='Frames' + str(i), linewidth = 3.5)
plt.ylabel('FSC', fontsize = 24)
plt.xlabel('frequencies', fontsize = 24)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.legend()
plt.show()
"""

# fsc_final = fsc_final_avg


#%%
#Trying something new which is similar to what relion does

""""
fsc_final = np.array(fsc_final)
fsc_final_avg = np.array(fsc_final_avg)

i = 0
# window_len = 6
fig , ax = plt.subplots(nrows = 4, ncols=6 )
for row in ax:
    for col in row:
        if i < len(fsc_final):
            print(i)
            fsc_final_i = smooth(fsc_final[i][3:], window_len)
            fsc_final_avg_i = smooth(fsc_final_avg[i][3:], window_len)
            fsc_final_i = fsc_final_i[:][(fsc_final_i > 0.143) & (fsc_final_i < 1)]
            fsc_final_avg_i = fsc_final_avg_i[:][(fsc_final_avg_i > 0.143) & (fsc_final_avg_i < 1)]

            fsc_final_i = fsc_final_i[0:len(fsc_final_avg_i)]
            rel_damp_sqrt = np.sqrt(
                np.divide(fsc_final_i - fsc_final_i * fsc_final_avg_i, fsc_final_avg_i - fsc_final_i * fsc_final_avg_i))
            reldamp_i = np.log(rel_damp_sqrt)

            col.plot(np.array(fsc_s)[0][0][0][0:len(reldamp_i)], reldamp_i, label='Frame' + str(i))
            # col.plot(fsc_s[0][0][0],smooth(fsc_final_avg[i], window_len), label='Average' + str(i))
            i+= 1
            col.legend()

plt.show()
"""



#%%


# i= 4
#
# offset = 6
# fit_params = fitthecurve([fsc_s[i][0][0][offset:] ,  smooth(fsc_final[i][offset:] , window_len) ], [-1.24822667, -6.90053422, -2.63049508])
# print(fit_params[2])
#
# plt.figure()
# plt.plot(fsc_s[i][0][0]  ,  smooth(fsc_final[i] , window_len) , label = 'original')
# plt.plot(fsc_s[i][0][0][offset:]  ,  smooth(fsc_final[i][offset:] , window_len) , label = 'compress')
# plt.plot( fit_params[0] , fit_params[1], label = 'fit')
# plt.ylabel('FSC', fontsize = 24)
# plt.xlabel('frequencies', fontsize = 24)
# plt.xticks(fontsize = 24)
# plt.yticks(fontsize = 24)
# plt.legend()
# plt.show()
#
# plt.figure()
# plt.plot(np.array(fsc_s[i][0][0]) * np.array(fsc_s[i][0][0]) ,  smooth(fsc_final[i] , window_len) , label = 'original')
# plt.plot(np.array(fsc_s[i][0][0][offset:]) *np.array(fsc_s[i][0][0][offset:]) ,  smooth(fsc_final[i][offset:] , window_len) , label = 'compress')
# plt.plot(np.array(fit_params[0]) * np.array(fit_params[0]) , fit_params[1], label = 'fit')
# plt.ylabel('FSC', fontsize = 24)
# plt.xlabel('frequencies', fontsize = 24)
# plt.xticks(fontsize = 24)
# plt.yticks(fontsize = 24)
# plt.legend()
# plt.show()


# def fitfunc_lin(x, a, b):
#     return a*x+b
#
# data_x = np.array(fsc_s[i][0][0])
# data_y = smooth(fsc_final[i] , window_len)
# data_y[data_y <=0] = 0.0000000001
# data_y = np.log(data_y)
# coeff, var_matrix = curve_fit(fitfunc_lin, data_x, data_y)
# plt.figure()
# plt.plot(data_x ,  data_y , label = 'original')
# plt.plot(data_x ,  fitfunc_lin(data_x, *coeff) , label = 'fit')
# plt.legend()
# plt.show()

"""
B_values = []
coeff_list = []
offset_start = 0
offset_end  = 20

bump = []
for i in range(len(fsc_final)):
    logvalues = np.log(smooth(fsc_final[i], window_len))
    bump.append(np.where(np.isnan(logvalues))[0][0])



int_B = 0
for i in range (len(fsc_s)):
    print(i)
    data_x = np.array(fsc_s[i][0][0][offset_start:bump[i]])
    data_y = smooth(fsc_final[i][offset_start:bump[i]], window_len)
    data_y[data_y <= 0] = 0.0000000001
    coeff, var_matrix = curve_fit(fitfunc, data_x, data_y)
    B_values.append(coeff[1])
    coeff_list.append(coeff)


plt.figure()
plt.plot(np.arange(24), B_values, 'o-')
plt.legend()
plt.xlabel('Frames', fontsize = 24)
plt.ylabel('B-Factor', fontsize = 24)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.show()

"""


# for i in range (0,len(fsc_s),4):
#     plt.figure()
#     fitcurv = fitfunc(fsc_s[i][0][0][p:110], *coeff_list[i])
#     plt.plot(np.array(fsc_s[i][0][0][p:110]), smooth(fsc_final[i][p:110] , window_len) , label = 'original' + str(i) )
#     plt.plot(np.array(fsc_s[i][0][0][p:110]), fitcurv, label='fit'+ str(i))
#     plt.xlabel('Frames', fontsize=24)
#     plt.ylabel('B-Factor', fontsize=24)
#     plt.xticks(fontsize=24)
#     plt.yticks(fontsize=24)
#     plt.legend()
#     plt.show()

"""
B_values_linear = []
coeff_list = []
for i in range (len(fsc_s)):
    print(i)
    data_x = np.array(fsc_s[i][0][0][offset_start:bump[i]])
    data_y = smooth(fsc_final[i][offset_start:bump[i]], window_len)
    data_y[data_y <= 0] = 0.0000000001
    data_y = np.log(data_y)
    coeff, var_matrix = curve_fit(fitslopefunc, data_x, data_y)
    B_values_linear.append(coeff[0])
    coeff_list.append(coeff)

i = 0
fig , ax = plt.subplots(nrows = 4, ncols=6 )
for row in ax:
    for col in row:
        if i < len(fsc_s):
            fitcurv = fitslopefunc(np.array(fsc_s[i][0][0][offset_start:bump[i]]), *coeff_list[i])
            col.plot(np.array(fsc_s[i][0][0][offset_start:bump[i]]), np.log(smooth(fsc_final[i][offset_start:bump[i]], window_len)),
             'o-',  label='original' + str(i), linewidth = 6)
            col.plot(np.array(fsc_s[i][0][0][offset_start:bump[i]]), fitcurv, 'o-',label='fit' + str(i))
            i+= 1
            col.legend()
            # col.set_yscale('log')
            # col.set_ylim(10e-4, 1)

plt.show()
"""

# data_x = np.array(fsc_s[i][0][0])**2
# data_y = smooth(fsc_final[i] , window_len)
# data_y[data_y <=0] = 0.0000000001
# data_y = np.log(data_y)


# plt.figure()
# plt.plot(1/np.array(fsc_s)[0,0,0,bump])
# plt.show()

"""
bumpifsc = np.array(fsc_s)[0,0,0,bump]
minfsc = np.min(bumpifsc)
fsc_s_norm = bumpifsc -  minfsc
# maxfsc = np.max(np.array(fsc_s_norm)[:])
maxfsc = 1



minbval = np.min(np.array(B_values)[:])
bval_norm = np.array(B_values)[:] - minbval
# maxbval = np.max(np.array(bval_norm)[:])
maxbval = 1




minbval_lin = np.min(np.array(B_values_linear)[:])
bval_lin_norm = np.array(B_values_linear)[:] - minbval_lin
# maxbval_lin = np.max(np.array(bval_lin_norm)[:])
maxbval_lin = 1



coeff, var_matrix = curve_fit(fitslopefunc,  np.arange(24), bval_norm)
fitline = coeff[0] * np.arange(24) + coeff[1]

plt.figure()
plt.plot(np.arange(24), fsc_s_norm[:], 'o-', label = 'bump')
plt.plot(np.arange(24), bval_norm, 'o-', label = 'Bval')
plt.plot(np.arange(24), fitline, 'o-', label = 'Fit')
plt.plot(np.arange(24), bval_lin_norm, 'o-', label = 'bvalFit')
plt.xlabel('Frames',fontsize = 24)
plt.ylabel('1st bump at frequency',fontsize = 24)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.legend(fontsize = 24)
plt.show()





SNR = np.divide(2* np.array(fsc_final)[4,:] , 1- np.array(fsc_final)[4,:])


plt.figure()
plt.plot(np.array(fsc_s)[4,0,0,:], np.log(smooth(np.array(SNR)[:], window_len)) , 'o')
plt.show()

"""

#%%
# fsc_all = []
# plt.figure()
# idx = 0
# for i in range(0, len(particle_imgs[0]), 11):
#     fsc_all.append(sp_statistics.fsc(particle_imgs[0][i], project_ptcl_all[0]))
#     plt.plot(fsc_all[idx][0], moving_avg_filter(fsc_all[idx][1]), label=str(i))
#     idx+=1
#
#
# fsc_avg = sp_statistics.fsc(particle_avg, projection_2D)
#
# plt.plot(fsc_avg[0], fsc_avg[1], label = "Averaged")
# plt.plot(fsc_avg[0], moving_avg_filter(fsc_avg[1]), label = "filter Averaged")
# plt.legend()
# plt.show()

# fsc_all = []
# idx = 0
#
# plt.figure()
# for j in range (3):
#     for i in range(0, len(particle_imgs[0]), 11):
#         fsc_all.append(sp_statistics.fsc(particle_imgs[j][i], project_2D_all[j]))
#         plt.plot(fsc_all[idx][0], moving_avg_filter(fsc_all[idx][1]), label=str(j)+ "," + str(i))
#         idx+=1
# plt.legend()


# fsc_avg_all = []
# plt.figure()
# idx = 0
# for avg in range (10):
#     fsc_all.append(sp_statistics.fsc(particle_avg_list[avg], project_2D_all[avg]))
#     plt.plot(fsc_all[idx][0], moving_avg_filter(fsc_all[idx][1]), label= "Average" + str(avg))
#     idx += 1

# plt.legend()
# plt.show()





# fsc_avg = []
# fsc_per_frame= []
# idx = 0
# window_len = 5
# plt.figure()
# for i in range(0, len(particle_imgs[0]), 11):
#     fsc_per_frame.append(sp_statistics.fsc(particle_imgs[i][0], project_2D_all[i]))
#     plt.plot(fsc_all[idx][0], fsc_all[idx][1], 'o-', label = "wtout smoth" + str(i))
#     # fsc_all[idx][0] = np.r_[fsc_all[idx][0][window_len - 1:0:-1], fsc_all[idx][0], fsc_all[idx][0][-2:-window_len - 1:-1]]
#
#
#     plt.plot(fsc_all[idx][0], smooth(fsc_all[idx][1], window_len),'o-', label="with smooth" + str(i))
#     idx += 1
#
#
# plt.legend()
# plt.show()
