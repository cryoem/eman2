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

from scipy.optimize import curve_fit
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

print(ima_start, ima_end)



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
Reads the x and y values per frame in a micrograph  
"""
def returns_values_in_file(f, mode = 'r'):
    """
    read a file and returns all its lines
    :param f: path to file
    :param mode: how open the file. Default: read file
    :return: contained values
    """
    if os.path.isfile(f):
        f1 = open(f, mode)
        values_f1 = f1.readlines()
        f1.close()
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
    return x,y



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
            # print(i, part_cord[j][0] - int(sx[i]), nxx, part_cord[j][1] - int(sy[i]), nyy, cen_xx, cen_yy, cen_zz, int(sx[i]) , int(sy[i])  )
            crop_imgs.append(Util.window(frames_i, nxx, nyy, 1, part_cord[j][0] - cen_xx,
                                         part_cord[j][1] - cen_yy , i - cen_zz))
        particle_imgs_in_movie.append(crop_imgs)
    return particle_imgs_in_movie


def moving_avg_filter(fsc_curve):
    for i in range(2, len(fsc_curve) -3 ) :
        fsc_curve[i] = (fsc_curve[i] + fsc_curve[i-1] + fsc_curve[i-2] + fsc_curve[i+1] + fsc_curve[i+2])/5
    return fsc_curve

def smooth(x,window_len):
    import numpy as np
    w = np.ones(window_len, 'd')
    y = np.convolve(w / w.sum(), x, mode='same')
    return y


def fitfunc(x, a, b,c,d ):
    return -a * np.exp(c + (4*b * (np.array(x)*np.array(x)))) + d


#%%

"""
Reading of Movie frames in .mrc format and display one frame
"""
ABSOLUTE_PATH_TO_MRC_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/MOVIES/"
input_image_path = os.path.join(ABSOLUTE_PATH_TO_MRC_FOLDER, "TcdA1-*_frames.mrc")

movie_names = return_movie_names(input_image_path)

ABSOLUTE_PATH_TO_LOG_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/corrsum_dw_log/"
log_movie_path = os.path.join(ABSOLUTE_PATH_TO_LOG_FOLDER, "TcdA1-*_frames.log")
shift_movie_files = return_movie_names(log_movie_path)

# -------------- Loading the reference volume
ref_vol_filename = "/home/adnan/PycharmProjects/DoseWeighting/vol_combined.hdf"
ref_volume = sp_utilities.get_im(ref_vol_filename)
# ---------------Preparing the volume for calculation of projections
volft = sp_projection.prep_vol(ref_volume, npad=2,
                               interpolation_method=1)  # transforms the volume in fourier space and then expands the volume with npad and centers the volume and returns the volume

read_meta_shifts
fsc_values = []
fsc_avgs = []
frequencies = []
fsc_raw_all = []

for micro in enumerate(movie_names[ima_start:ima_end]):

    # micro = (0, '/home/adnan/PycharmProjects/DoseWeighting/MOVIES/TcdA1-0100_frames.mrc')

    frames = return_images_from_movie(micro[1], show_first = False)

    logfile = ABSOLUTE_PATH_TO_LOG_FOLDER + micro[1].split('.')[0].split('/')[-1] + '.log'
    shift_x , shift_y = read_meta_shifts(logfile)

    for i in range (frames.get_zsize()):
        reg = EMAN2_cppwrap.Region(0,0,i,frames.get_xsize(),frames.get_ysize(),1)

        shift_img = sp_fundamentals.fshift(frames.get_clip(reg),shift_x[i], shift_y[i] )
        frames.insert_clip(shift_img,(0,0,i) )

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



    ref_project_2D_ptcl_all = get_2D_project_for_all_ptcl_from_reference(volft, project_params, show = False) #Projection of 3D volume in 2-D for all the particles in all frames in one movie
    mask = sp_utilities.model_circle(nx[0] / 2, nx[0], nx[0])  # nx/2 is for the radius
    for i in range(len(ref_project_2D_ptcl_all)):
        ref_project_2D_ptcl_all[i] = sp_filter.filt_ctf(ref_project_2D_ptcl_all[i], ctf_params[i], binary=True)

    """
    Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
    creating a window around to extract the same particle from each frame
    """
    cen_x = frames.get_xsize() // 2
    cen_y = frames.get_ysize() // 2
    cen_z = frames.get_zsize() // 2

    particle_imgs = get_all_reduce_ptcl_imgs( frames, mask, nx[0], ny[0], particle_coordinates, ctf_params, cen_x, cen_y, cen_z)

    zsize = frames.get_zsize()
    del frames


    """
    Calculating the fourier shell correlation of all the particle images with respect to 2-D reference projection of 3-D volume
    """
    fsc_freq = []
    fsc_val = []
    for i in range (zsize):
        fsc_frames_freq = []
        fsc_frames_val = []
        for j in range(len(particle_imgs)):
            fsc_frames_freq.append(sp_statistics.fsc(particle_imgs[j][i], ref_project_2D_ptcl_all[j])[0])
            fsc_frames_val.append(np.array(sp_statistics.fsc(particle_imgs[j][i], ref_project_2D_ptcl_all[j]))[1])
        fsc_freq.append(fsc_frames_freq)
        fsc_val.append(fsc_frames_val)


    fsc_final = []
    for i in range (zsize):
        fsc_sum = [entry/len(fsc_val[i]) for entry in fsc_val[i][0]]
        for fsc in fsc_val[i][1:]:  # one frame ahead for averageing
            for idx, ok in enumerate(fsc):
                fsc_sum[idx] += ok/len(fsc_val[i])
        fsc_final.append(fsc_sum)

    fsc_final_avg =  []
    for idx in range (0,len(fsc_final)-3):
        avv = []
        for p in range(len(fsc_final[idx])):
            avv.append((fsc_final[idx][p] + fsc_final[idx+1][p] + fsc_final[idx+2][p] + fsc_final[idx+3][p] )/4)
        fsc_final_avg.append(avv)

    for idx in range(len(fsc_final) - 3, len(fsc_final)):
        avv = []
        for p in range(len(fsc_final[idx])):
            avv.append((fsc_final[idx][p] + fsc_final[idx-1][p] + fsc_final[idx-2][p] + fsc_final[idx-3][p] )/4)
        fsc_final_avg.append(avv)

    fsc_values.append(fsc_final)
    fsc_avgs.append(fsc_final_avg)
    frequencies.append(fsc_freq[0][0])
    fsc_raw_all.append(np.array(fsc_val))

    del particle_imgs
    del ref_project_2D_ptcl_all

fsc_values_per_micrograph = sp_utilities.wrap_mpi_gatherv(fsc_values, 0, mpi.MPI_COMM_WORLD)
fsc_avgs_per_micrograph = sp_utilities.wrap_mpi_gatherv(fsc_avgs, 0, mpi.MPI_COMM_WORLD)
freq_per_micrograph = sp_utilities.wrap_mpi_gatherv(frequencies, 0, mpi.MPI_COMM_WORLD)
fsc_raw =  sp_utilities.wrap_mpi_gatherv(fsc_raw_all, 0 , mpi.MPI_COMM_WORLD)

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


    with open(str(no_of_micrographs) + '_Micrograph_raw_fsc_values', 'wb') as wb:
        for value in fsc_raw:
            pickle.dump(value, wb)
    wb.close()


print("I am finish")

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
mpi.mpi_finalize()

#%%

# """

fsc_values_per_micrograph = []
fsc_avgs_per_micrograph = []
freq_per_micrograph = []
fsc_raw = []


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


with open(location + '/sphire/bin/' + str(no_of_micrographs) + '_Micrograph_raw_fsc_values', 'rb') as rb:
    while True:
        try:
            fsc_raw.append(pickle.load(rb))
        except EOFError:
            break
rb.close()



print(np.array(fsc_values_per_micrograph).shape)
print(np.array(fsc_avgs_per_micrograph).shape)
print(np.array(freq_per_micrograph).shape)


fsc_values_per_micrograph = np.array(fsc_values_per_micrograph)
fsc_avgs_per_micrograph = np.array(fsc_avgs_per_micrograph)
freq_per_micrograph = np.array(freq_per_micrograph)




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

fsc_sum_per_frame = np.array(fsc_sum_per_frame)[:] * -1
fsc_sum_avg_per_frame = np.array(fsc_sum_avg_per_frame)[:] * -1


#%%
offset_start = 0
offset_end = -1
frames_range = 24
window_len = 1
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

ABSOLUTE_PATH_TO_LOG_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/corrsum_dw_log/"
logfile = ABSOLUTE_PATH_TO_LOG_FOLDER + micro[1].split('.')[0].split('/')[-1] + '.log'

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

shift_x, shift_y = read_meta_shifts(logfile)
mask = sp_utilities.model_circle(nx[0] / 2, nx[0], nx[0])
particle_imgs = get_all_reduce_ptcl_imgs(frames, mask, nx[0], ny[0], particle_coordinates, ctf_params, cen_x, cen_y, cen_z)

particle_imgs = np.array(particle_imgs).swapaxes(0,1)
img = particle_imgs[0][0].get_2dview()
four_img_new_3 = np.fft.fft2(img)
four_img_new_3  = np.fft.fftshift(four_img_new_3)

mask_distance = np.zeros((352,352))
row, col = np.meshgrid(range(352), range(352), indexing = 'ij')
for i in range(np.shape(row)[0]):
    for j in range(np.shape(col)[0]):
        mask_distance[i][j] = np.sqrt((row[i][j] - np.shape(row)[0]/2)**2   + (col[i][j] - np.shape(col)[0]/2)**2 )

normal_to = mask_distance[ mask_distance.shape[0]/2][mask_distance.shape[1]-1]
mask_norm = (mask_distance / normal_to) * 0.5


index_frame = 0
mask_applied = np.zeros((352,352))
for i in range (mask_distance.shape[0]):
    for j in range(mask_distance.shape[1]):
        near_value = (np.abs(weight_per_frame[index_frame] - mask_norm[i][j])).argmin()
        mask_applied[i][j] = weight_per_frame[index_frame][near_value]


plt.figure()
plt.imshow(four_img_new_3.real ,cmap = plt.get_cmap('Greys'))
plt.colorbar()
plt.clim(-4000,4000)

plt.figure()
plt.imshow(mask_norm)
plt.colorbar()

plt.figure()
plt.imshow(mask_applied)
plt.colorbar()

new_img = four_img_new_3.real  * mask_applied


new_img_shift = np.fft.ifftshift(new_img)
new_img_shift = np.fft.ifft2(new_img_shift)

plt.figure()
plt.imshow(new_img_shift.real,cmap = plt.get_cmap('Greys'))
plt.colorbar()


plt.figure()
plt.imshow(img,cmap = plt.get_cmap('Greys'))
plt.colorbar()






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



    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(5):
    #     if i == 0:
    #         shift_x_app = shift_x
    #         shift_y_app = shift_y
    #         name = 'add'
    #     elif i == 1:
    #         shift_x_app = [-entry for entry in shift_x]
    #         shift_y_app = [-entry for entry in shift_y]
    #         name = 'subtract'
    #     elif i == 2:
    #         shift_x_app = [entry for entry in shift_x]
    #         shift_y_app = [-entry for entry in shift_y]
    #         name = 'pxmy'
    #     elif i == 3:
    #         shift_x_app = [-entry for entry in shift_x]
    #         shift_y_app = [entry for entry in shift_y]
    #         name = 'mxpy'
    #     elif i == 4:
    #         shift_x_app = [0 for entry in shift_x]
    #         shift_y_app = [0 for entry in shift_y]
    #         name = 'normal'
    #     particle_imgs = get_all_reduce_ptcl_imgs(
    #         frames, mask, nx[0], ny[0], [particle_coordinates[num]],
    #         ctf_params, cen_x, cen_y, cen_z, shift_x_app , shift_y_app
    #     )
    #     print(particle_imgs)
    #     print(np.array(particle_coordinates).shape  )
    #
    #     avg_part = sum(particle_imgs[0])
    #     fsc = sp_statistics.fsc(avg_part, ref_project_2D_ptcl_all[num])
    #     print(name, fsc[1])
    #     ax.plot(fsc[0], fsc[1], label=str(i))
    #
    #     plt.figure()
    #     plt.imshow(avg_part.get_2dview()[::-1],cmap = plt.get_cmap('Greys'))
    #     plt.colorbar()
    #     plt.title(str(micro[1]) )
    #     plt.savefig(str(micro[1]) + '_{0}.png'.format(name))
    #     plt.clf()
    #
    # ax.legend()
    # fig.savefig(str(micro[1]) + '_fsc.png')


"""