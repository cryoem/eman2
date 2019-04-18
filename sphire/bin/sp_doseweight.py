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


def return_movie_framenames(input_image_path):
    mic_pattern = input_image_path
    mic_basename_pattern = os.path.basename(mic_pattern)
    global_entry_dict = {}
    subkey_input_mic_path = "Input Micrograph Path"

    mic_basename_tokens = mic_basename_pattern.split('*')
    mic_id_substr_head_idx = len(mic_basename_tokens[0])

    input_mic_path_list = glob.glob(mic_pattern)

    print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))

    for input_mic_path in input_mic_path_list:
        # Find tail index of  id substring and extract the substring from the  name
        input_mic_basename = os.path.basename(input_mic_path)
        mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
        mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
        if not mic_id_substr in global_entry_dict:
            # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
            global_entry_dict[mic_id_substr] = {}
        global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path

    print(" ")
    print("Summary of dataset consistency check...")
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


    print("Following mrc files were detected from location", input_image_path)
    for input_mic_path in input_file_path_list:
        print(" ",  os.path.basename(input_mic_path))

    return input_file_path_list


"""
Reading of Movie frames in .mrc format and display one frame
"""

ABSOLUTE_PATH_TO_MRC_FOLDER= "/home/adnan/PycharmProjects/DoseWeighting/MOVIES/"
input_image_path = os.path.join(ABSOLUTE_PATH_TO_MRC_FOLDER, "TcdA1-*_frames.mrc")


frame_names = return_movie_framenames(input_image_path)
frames = []
for ima in range(len(frame_names)):
    frames.append(EMAN2_cppwrap.EMData.read_images(frame_names[ima]))
    break


print(frame_names)
plt.figure()
plt.imshow(frames[0][0].get_3dview()[0], cmap = plt.get_cmap('Greys'))
plt.colorbar()
# plt.show()


"""
Reading a particle stack to find all the parameters saved in the header file
"""
stackfilename = "bdb:/home/adnan/PycharmProjects/DoseWeighting/Substack/isac_substack"

ima = EMAN2_cppwrap.EMData()
ima.read_image(stackfilename, 0, False)

#----------------------------------------------------------------
#  Reading the projection parameters
#-----------------------------------------------------------------

# -------------------------------------------  2D orientation / alignment attributes
"""
alpha = rotation angle
sx = shift in x direction in image plane (2-D)
sy = shift in y direction in image plane (2-D)
mirror =  0 ,1  whether the image is mirrored or not
scale = scale of the image
"""
alpha,sx, sy, mirror , scale  = sp_utilities.get_params2D(ima, xform = "xform.align2d")

print(alpha,sx, sy, mirror , scale )

# -------------------------------------------  2D orientation / orientation attributes
"""
phi =  Eulerian angle for 3D reconstruction (azimuthal) 
theta = Eulerian angle for 3D reconstruction (tilt) 
psi = Eulerian angle for 3D reconstruction (in-plane rotation of projection) 
tx =  shift in x direction
ty = shift in y direction
"""
phi, theta, psi, tx, ty  = sp_utilities.get_params_proj(ima, xform = "xform.projection")

print(phi, theta, psi, tx, ty)

#-----------------------------------------------   CTF related attributes
"""
defocus = defocus associated with the image, positive value corresponds to underfocus
cs =  spherical aberration constant [mm]. 
voltage = accelerating voltage of the microscope [kV] 
apix = 
bfactor = The parameter in Guassian like envelope function, which roughly explains Fourier factor dumping of the image. 
ampcont = amplitude contrast
dfdiff  = astigmatism amplitude
dfang =  astigmatism angle
"""

defocus, cs, voltage, apix, bfactor, ampcont, dfdiff , dfang = sp_utilities.get_ctf(ima)

print(defocus, cs, voltage, apix, bfactor, ampcont, dfdiff , dfang)

dict_projec = ima.get_attr_dict()

# for keys,values in dict_projec.items():
#     print("\n")
#     print(keys,"\n", values)

"""
Reading a reference map
"""

#-------------- Loading the reference volume
ref_vol_filename = "/home/adnan/PycharmProjects/DoseWeighting/vol_combined.hdf"
ref_volume = sp_utilities.get_im(ref_vol_filename)

#---------------Preparing the volume for calculation of projections
# What it does is that it transforms the volume in fourier space and then expands the volume with npad and centers the volume and returns the volume
volft  =  sp_projection.prep_vol(ref_volume, npad = 2, interpolation_method=1)

#---------------- Calculating the 2-D projection of a 3-D volume
params = [phi, theta, psi, -tx, -ty]
projection_2D  = sp_projection.prgl(volft, params, interpolation_method=1)    # Will return a 2-D projection


plt.figure()
plt.imshow(ima.get_2dview(), cmap = plt.get_cmap('Greys'))
plt.colorbar()
plt.clim(-4,4)

plt.figure()
plt.imshow(projection_2D.get_2dview(), cmap = plt.get_cmap('Greys'))
plt.colorbar()
# plt.show()


"""
Extracting particle image from the movie data. First getting the particle cordinates from the dictionary and then 
creating a window around to extract the same particle from each frame
"""

#----------------------- Particle cordinates
part_cordinates_x , part_cordinates_y = dict_projec["ptcl_source_coord"]

nx = dict_projec["nx"]
ny = dict_projec["ny"]
nz = dict_projec["nz"]

print(part_cordinates_x, part_cordinates_y)
print(nx, ny, nz)

frames = frames[0][0]
print(np.shape(frames.get_3dview()))
print(frames.get_xsize() ,frames.get_ysize(), frames.get_zsize())

cen_x = frames.get_xsize() // 2
cen_y = frames.get_ysize() // 2
cen_z = frames.get_zsize() // 2



# Extracting the particle from each image
particle_imgs = []
particle_avgs_list = []
mask = sp_utilities.model_circle(nx/2, nx, nx)
for i in range(frames.get_zsize()):
    particle_imgs.append(Util.window(frames, nx, ny, 1, part_cordinates_x - cen_x, part_cordinates_y - cen_y, i - cen_z))

    st = Util.infomask(particle_imgs[i], None, False)    # find statistics of the image returns avg, std, min, max
    Util.mul_scalar(particle_imgs[i], -1.0)
    particle_imgs[i] += 2 * st[0]  # st[0]-> avergage  , shifting the negative range to positive

    st = Util.infomask(particle_imgs[i], mask, False)
    particle_imgs[i] -= st[0]  # st[0]-> avergage
    particle_imgs[i] /= st[1]  # st[1]-> standard deviation

    st = Util.infomask(particle_imgs[i], mask, False)
    particle_imgs[i] -= st[0]   # st[0]-> avergage
    particle_imgs[i] = sp_filter.filt_ctf(particle_imgs[i], dict_projec["ctf"], binary = True)

    st = Util.infomask(particle_imgs[i], mask, False)
    particle_imgs[i] -= st[0]   # st[0]-> avergage
    st = Util.infomask(particle_imgs[i], mask, False)
    particle_imgs[i] /= st[1]   # st[1]-> standard deviation
    # print("stats for particle", Util.infomask(particle_imgs[i], mask, False))




print("stats for proj2D", Util.infomask(projection_2D, mask, False))

st = Util.infomask(projection_2D, mask, False)
projection_2D-= st[0]  # st[0]-> avergage
projection_2D /= st[1]  # st[1]-> standard deviation


print("stats for proj2D", Util.infomask(projection_2D, mask, False))

# ---------------------------- Calculating the average of all images to get a better contract for testing whether the
# ---------------------------- correct Particle was taken or not
particle_avg = sum(particle_imgs)

plt.figure()
plt.imshow(particle_avg.get_2dview(), cmap = plt.get_cmap('Greys'))
plt.colorbar()
# plt.show()


"""
Calculating the fourier shell correlation of all the particle images with respect to 2-D reference projection of 3-D volume
"""
fsc_all = []
plt.figure()
idx = 0
for i in range(0, len(particle_imgs), 11):
    fsc_all.append(sp_statistics.fsc(particle_imgs[i], projection_2D))
    plt.plot(fsc_all[idx][0], fsc_all[idx][1], label=str(i))
    idx+=1


fsc_avg = sp_statistics.fsc(particle_avg, projection_2D)

plt.plot(fsc_avg[0], fsc_avg[1], label = "Averaged")

plt.legend()
plt.show()



