from __future__ import print_function
import os
import glob
import EMAN2_cppwrap
import matplotlib.pyplot as plt
from matplotlib import cm

import sp_utilities

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
# plt.show()


"""
Reading a particle stack to find all the parameters saved in the header file
"""

stackfilename ="bdb:/home/adnan/PycharmProjects/eman2/sphire/tests/Class2D/stack_ali2d"
    # "bdb:/home/adnan/PycharmProjects/DoseWeighting/Substack_with_params/isac_substack"

# ftp = sp_utilities.file_type(stackfilename)
#
#
# nima = EMAN2_cppwrap.EMUtil.get_image_count(stackfilename)
# list_of_particles = list(range(nima))
# part_stack = EMAN2_cppwrap.EMData.read_images(stackfilename)

