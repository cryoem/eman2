from __future__ import print_function
import os
import glob
from EMAN2 import *
import EMAN2_cppwrap
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sp_utilities
import sp_projection
import sp_statistics
import sp_filter


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


def return_images_from_movie(movie_name_list, index, show_first = False):
    mov_imgs = EMAN2_cppwrap.EMData.read_images(movie_name_list[index])
    if show_first:
        plt.figure()
        plt.imshow(mov_imgs[0].get_3dview()[0], cmap=plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()
    print("\n No of images in the selected movie %s are %d" % (os.path.basename(movie_name_list[index]) ,len(mov_imgs[0].get_3dview() ) ))
    return mov_imgs


"""
Reading of Movie frames in .mrc format and display one frame
"""
ABSOLUTE_PATH_TO_MRC_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/MOVIES/"
input_image_path = os.path.join(ABSOLUTE_PATH_TO_MRC_FOLDER, "TcdA1-*_frames.mrc")

movie_names = return_movie_names(input_image_path)
print("\n Following movies were detected : \n",movie_names)
frames = return_images_from_movie(movie_names,0, show_first = True)


#%%
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

    return no_of_imgs_once, ptcl_source_images_once, project_params_all_once, particle_coordinates_all_once, ctf_params_all_once



def find_particles_info_from_movie(stack, movie_name, no_of_imgs, ptcl_source_images, project_params_all, particle_coordinates_all, ctf_params_all, show_first = False):

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

    for i in range(no_of_imgs):
        if (str(os.path.basename(movie_name)) == str(os.path.basename(ptcl_source_images[i]))):
            project_params_per_movie.append(project_params_all[i])
            particle_coordinates_per_movie.append(particle_coordinates_all[i])
            ctf_params_per_movie.append(ctf_params_all[i])

    print("Number of particles detected in %s are %d" % (str(os.path.basename(movie_name)), len(project_params_per_movie)))
    print("Ctf estimation parameters for 1st particle in the stack are ", ctf_params_per_movie[0].to_dict())
    print("Projection parameters for 1st particle in the stack are ", project_params_per_movie[0].get_params('spider'))


    if show_first:
        ima = EMAN2_cppwrap.EMData()
        ima.read_image(stack, 0, False)
        plt.figure()
        plt.imshow(ima.get_2dview(), cmap=plt.get_cmap('Greys'))
        plt.colorbar()
        plt.show()

    return project_params_per_movie, particle_coordinates_per_movie, ctf_params_per_movie

stackfilename = "bdb:/home/adnan/PycharmProjects/DoseWeighting/Substack/isac_substack"
no_of_imgs,ptcl_source_images,project_params_all,particle_coordinates_all,ctf_params_all = read_all_attributes_from_stack(stackfilename)
project_params, particle_coordinates, ctf_params = find_particles_info_from_movie(stackfilename, movie_names[0], \
                                                                                  no_of_imgs,ptcl_source_images,project_params_all,particle_coordinates_all,ctf_params_all,show_first=True)


#%%
"""
Reading a reference map
"""
#-------------- Loading the reference volume
ref_vol_filename = "/home/adnan/PycharmProjects/DoseWeighting/vol_combined.hdf"
ref_volume = sp_utilities.get_im(ref_vol_filename)

#---------------Preparing the volume for calculation of projections
# What it does is that it transforms the volume in fourier space and then expands the volume with npad and centers the volume and returns the volume
volft  =  sp_projection.prep_vol(ref_volume, npad = 2, interpolation_method=1)

project_2D_all = []
for i in range(len(project_params)):
    params_substack = project_params[i].get_params('spider')
    params_for_each_image =[params_substack['phi'], params_substack['theta'], params_substack['psi'], params_substack['tx'], params_substack['ty']]
    project_2D_all.append(sp_projection.prgl(volft, params_for_each_image, interpolation_method=1))
print("All projections", project_2D_all)

# plt.figure()
# plt.imshow(projection_2D.get_2dview(), cmap = plt.get_cmap('Greys'))
# plt.colorbar()
# plt.show()

"""
Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
creating a window around to extract the same particle from each frame
"""
#----------------------- Particle cordinates
dict_projec = ima.get_attr_dict()
part_cordinates_x , part_cordinates_y = dict_projec["ptcl_source_coord"]

nx = dict_projec["nx"]
ny = dict_projec["ny"]
nz = dict_projec["nz"]

print(part_cordinates_x, part_cordinates_y)
print(nx, ny, nz)

frames = frames[0]
print(np.shape(frames.get_3dview()))
print(frames.get_xsize() ,frames.get_ysize(), frames.get_zsize())

cen_x = frames.get_xsize() // 2
cen_y = frames.get_ysize() // 2
cen_z = frames.get_zsize() // 2

# print("Concentrate")
# print(ctf_params[0].to_dict())
# print(dict_projec["ctf"])


particle_imgs = []
mask = sp_utilities.model_circle(nx/2, nx, nx)  # nx/2 is for the radius


for j in range(len(particle_coordinates)):
    crop_imgs = []
    for i in range(frames.get_zsize()):
        crop_imgs.append(Util.window(frames, nx, ny, 1, particle_coordinates[j][0] - cen_x, particle_coordinates[j][1] - cen_y, i - cen_z))

        # particle_avg_image.append(Util.window(frames, nx, ny, 1, part_cordinates_x - cen_x, part_cordinates_y - cen_y, i - cen_z))
        st = Util.infomask(crop_imgs[i], None, False)    # find statistics of the image returns avg, std, min, max
        Util.mul_scalar(crop_imgs[i], -1.0)
        crop_imgs[i] += 2 * st[0]  # st[0]-> avergage  , shifting the negative range to positive

        st = Util.infomask(crop_imgs[i], mask, False)
        crop_imgs[i] -= st[0]  # st[0]-> avergage
        crop_imgs[i] /= st[1]  # st[1]-> standard deviation

        st = Util.infomask(crop_imgs[i], mask, False)
        crop_imgs[i] -= st[0]   # st[0]-> avergage
        crop_imgs[i] = sp_filter.filt_ctf(crop_imgs[i], ctf_params[j], binary = True)

        st = Util.infomask(crop_imgs[i], mask, False)
        crop_imgs[i] -= st[0]   # st[0]-> avergage
        st = Util.infomask(crop_imgs[i], mask, False)
        crop_imgs[i] /= st[1]   # st[1]-> standard deviation
    particle_imgs.append(crop_imgs)



# print("stats for proj2D", Util.infomask(projection_2D, mask, False))
#
# st = Util.infomask(projection_2D, mask, False)
# projection_2D-= st[0]  # st[0]-> average
# projection_2D /= st[1]  # st[1]-> standard deviation
#
# print("stats for proj2D", Util.infomask(projection_2D, mask, False))


for i in range(len(project_2D_all)):
    st = Util.infomask(project_2D_all[i], mask, False)
    project_2D_all[i]-= st[0]  # st[0]-> average
    project_2D_all[i] /= st[1]  # st[1]-> standard deviation

# ---------------------------- Calculating the average of all images to get a better contract for testing whether the
# ---------------------------- correct Particle was taken or not

particle_avg_list = []
for i in range(len(particle_imgs)):
    particle_avg_list.append(sum(particle_imgs[i]))



plt.figure()
plt.imshow(particle_avg_list[0].get_2dview(), cmap = plt.get_cmap('Greys'))
plt.colorbar()
# plt.show()



#%%
"""
Calculating the fourier shell correlation of all the particle images with respect to 2-D reference projection of 3-D volume
"""

def moving_avg_filter(fsc_curve):
    for i in range(2, len(fsc_curve) -3 ) :
        fsc_curve[i] = (fsc_curve[i] + fsc_curve[i-1] + fsc_curve[i-2] + fsc_curve[i+1] + fsc_curve[i+2])/5
    return fsc_curve


fsc_all = []
plt.figure()
idx = 0
for i in range(0, len(particle_imgs[0]), 11):
    fsc_all.append(sp_statistics.fsc(particle_imgs[0][i], projection_2D))
    plt.plot(fsc_all[idx][0], moving_avg_filter(fsc_all[idx][1]), label=str(i))
    idx+=1


fsc_avg = sp_statistics.fsc(particle_avg, projection_2D)

plt.plot(fsc_avg[0], fsc_avg[1], label = "Averaged")
plt.plot(fsc_avg[0], moving_avg_filter(fsc_avg[1]), label = "filter Averaged")
plt.legend()
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

plt.legend()
# plt.show()

def smooth(x,window_len):
    import numpy as np
    # s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    w = np.ones(window_len, 'd')
    y = np.convolve(w / w.sum(), x, mode='same')
    return y



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


#%%
fsc_s = []
for i in range (frames.get_zsize()):
    fsc_frames = []
    for j in range(len(particle_imgs)):
        fsc_frames.append(sp_statistics.fsc(particle_imgs[j][i], project_2D_all[j]))
    fsc_s.append(fsc_frames)

# print([entry for entry in fsc_s[0][0][1]])


fsc_final = []
for i in range (frames.get_zsize()):
    fsc_sum = [entry/len(fsc_s[i]) for entry in fsc_s[i][0][1]]
    for fsc in fsc_s[i][1:]:
        for idx, ok in enumerate(fsc[1]):
            fsc_sum[idx] += ok/len(fsc_s[i])
    fsc_final.append(fsc_sum)


window_len = 5
plt.figure()
for i in range (frames.get_zsize()):
    # plt.plot(fsc_s[i][0][0], fsc_final[i], label="without smoothing"+ str(i))
    plt.plot(fsc_s[i][0][0], smooth(fsc_final[i], window_len), label=str(i), linewidth = 3.5)
#
plt.legend()
plt.show()