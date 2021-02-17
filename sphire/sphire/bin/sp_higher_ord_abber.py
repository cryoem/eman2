#!/usr/bin/env python
##
import tqdm
import os
import subprocess
import sys
import argparse
import time
import json
from glob import glob
import numpy as np
import pandas as pd
import shutil

try:
    from ..libpy import sp_global_def
    from ..libpy import sp_utilities
except:
    from sphire.libpy import sp_global_def
    from sphire.libpy import sp_utilities

from pyStarDB import sp_pystardb as star

"""
import mpi
os.unsetenv('OMPI_COMM_WORLD_RANK')
RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
if RUNNING_UNDER_MPI :
    mpi.mpi_init(0, [])
    rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    size = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
else :
    rank  =0
    size = 1


env = os.environ
new_env = {k: v for k, v in env.items() if "MPI" not in k}
"""

def parse_parameters(args):
    parser = argparse.ArgumentParser(args, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "post_refiner",
        type=str,
        help="post refiner directory",
    )

    parser.add_argument(
        "corr_mic",
        type=str,
        help="Motion corr corrected micrographs",
    )

    parser.add_argument(
        "Output_folder",
        type=str,
        help="output_directory",
    )

    parser.add_argument(
        "--estimate_magnification",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--estimate_beamtilt",
        type=str,
        help="",
        default = False,
    )

    parser.add_argument(
        "--estimate_trefoil",
        type=str,
        help="",
        default=False

    )
    parser.add_argument(
        "--estimate_order_aberation",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--perform_CTF_params_fit",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_defcous_micrograph",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_defcous_particle",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_astigmatism_micrograph",
        type=str,
        help="",
        default=False
    )

    parser.add_argument(
        "--submission_template",
        type = str,
        help ="",
        default = ""
    )

    parser.add_argument(
        "--submission_command",
        type=str,
        help="",
        default=""
    )

    parser.add_argument(
        "--mpi_procs",
        type= int,
        help ="",
        default = 1
    )

    parser.add_argument(
        "--no_of_threads",
        type=int,
        help="",
        default= 1,
    )

    return parser.parse_args()


def get_final_iter(tracker_file):
    with open(tracker_file) as jsonfile:
        Tracker = json.load(jsonfile)
    return Tracker['constants']['best']


def get_bdb_loc(tracker_file):
    with open(tracker_file) as jsonfile:
        Tracker = json.load(jsonfile)
    return Tracker['constants']['stack']


def sanity_function_check(input_folder):
    # This is important and will be written as soon as the post refiner part is done.
    pass


def get_final_iter_files(input_folder):
    sanity_function_check(input_folder)
    if os.path.isdir(input_folder):
        if os.path.exists(os.path.join(input_folder, 'Tracker_final.json')):
            tracker_file = os.path.join(input_folder, 'Tracker_final.json')
            use_iter = 1
        else:
            main_directory = sorted(glob(os.path.join(input_folder, 'main*')))[-1]
            tracker_file = glob(os.path.join(main_directory, 'Tracker_*.json'))[0]
            use_iter = -1

        if use_iter == 1:
            final_iter = get_final_iter(str(tracker_file))
            bdbfile = get_bdb_loc(tracker_file)
            tracker_file = os.path.join(str(input_folder), 'main{:03d}'.format(final_iter))

        else:
            final_iter = get_final_iter(str(tracker_file))
            bdbfile = get_bdb_loc(tracker_file)
            tracker_file = os.path.join(str(input_folder), 'main{:03d}'.format(final_iter))

    template = '{0}_{1:03d}.{2}'
    required_data = [['params', 'txt'], ['chunk_0', 'txt'], ['chunk_1', 'txt']]
    required_files = []
    for i, j in required_data:
        required_files.append(os.path.join(tracker_file, template.format(i, final_iter, j)))
    return required_files, bdbfile


def parse_postrefiner(args_post):
    parser_post = argparse.ArgumentParser(args_post, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_post.add_argument(
        "command_python",
        type=str,
        help="",
    )

    parser_post.add_argument(
        "--combinemaps",
        type=str,
        help="",
        nargs=2
    )

    parser_post.add_argument(
        "--output_dir",
        type=str,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--pixel_size",
        type=float,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--mask",
        default="",
        type=str,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--do_adaptive_mask",
        help="",
        action='store_true',
    )

    parser_post.add_argument(
        "--threshold",
        type=float,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--edge_width",
        type=int,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--ndilation",
        type=int,
        help="",
    )

    parser_post.add_argument(
        "--mtf",
        type=str,
        help="",
    )

    parser_post.add_argument(
        "--B_enhance",
        type=float,
        help="",
        nargs=1
    )

    parser_post.add_argument(
        "--f",
        type=float,
        help="",
        nargs=1
    )

    return parser_post.parse_args(args_post.split())


def run(args):
    options = parse_parameters(args)

    #################################################
    ############Post process call#################
    #################################################
    ###### Post process part starts here
    if os.path.isdir(options.post_refiner):
        with open(os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'command.txt')), 'r') as commandfile:
            read_command_txt = commandfile.read()

        post_refine_options = parse_postrefiner(read_command_txt)

        if post_refine_options.mask is not "":
            use_mask = os.path.join(os.getcwd(), os.path.join(options.post_refiner, post_refine_options.mask))
        else:
            use_mask = os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'vol_adaptive_mask.hdf'))

        try:
            shutil.rmtree(os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),'PostProcess')))
        except FileNotFoundError:
            os.makedirs(os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),'PostProcess')))
            pass

        pp_star_file = star.StarFile(os.path.join(os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),'PostProcess')),
                                                  'postprocess.star'))

        half_map1 = os.path.join(os.getcwd(), post_refine_options.combinemaps[0])
        half_map2 = os.path.join(os.getcwd(), post_refine_options.combinemaps[1])

        if half_map1.endswith('hdf'):
            half_1_call = (
                    "e2proc3d.py"
                    + " " + str(half_map1)
                    + " " + str(half_map1).replace('hdf', 'mrc')
            )

            half_2_call = (
                    "e2proc3d.py"
                    + " " + str(half_map2)
                    + " " + str(half_map2).replace('hdf', 'mrc')
            )
            subprocess.run(args=[half_1_call], shell=True, text=True)
            subprocess.run(args=[half_2_call], shell=True, text=True)
            half_map1 = half_map1.replace('hdf', 'mrc')
            half_map2 = half_map2.replace('hdf', 'mrc')


        if use_mask.endswith('hdf'):
            mask_call = (
                    "e2proc3d.py"
                    + " " + str(use_mask)
                    + " " + str(use_mask).replace('hdf', 'mrc')
            )
            subprocess.run(args=[mask_call], shell=True, text=True)
            use_mask = use_mask.replace('hdf', 'mrc')

        general_pp = pd.DataFrame([[half_map1, half_map2, use_mask]],
                                  columns=['_rlnUnfilteredMapHalf1', '_rlnUnfilteredMapHalf2', '_rlnMaskName'])

        pp_star_file.update('general', general_pp, False)

        fsc_halves = np.loadtxt(os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'halves.txt')))
        fsc_masked_halves = np.loadtxt(
            os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'masked_halves.txt')))
        fsc_full = np.loadtxt(os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'full.txt')))

        spectral_index = fsc_halves[:, 0]
        angs_resolution = fsc_halves[:, 1]
        resolution = np.divide(1, angs_resolution)
        fsc_corrected = fsc_full[:, 2]
        fsc_unmasked_maps = fsc_halves[:, 2]
        fsc_masked_maps = fsc_masked_halves[:, 2]

        fsc_data_pp = pd.DataFrame(np.array([spectral_index, resolution, angs_resolution, fsc_corrected,
                                             fsc_unmasked_maps, fsc_masked_maps]).swapaxes(0, 1).tolist(),
                                   columns=['_rlnSpectralIndex', '_rlnResolution', '_rlnAngstromResolution',
                                            '_rlnFourierShellCorrelationCorrected',
                                            '_rlnFourierShellCorrelationUnmaskedMaps',
                                            '_rlnFourierShellCorrelationMaskedMaps'
                                            ])

        for col in fsc_data_pp.columns:
            if col == '_rlnSpectralIndex':
                fsc_data_pp[col] = fsc_data_pp[col].map(lambda x: int(x))
            else:
                fsc_data_pp[col] = fsc_data_pp[col].map(lambda x: '{0:0.4f}'.format(x))

        pp_star_file.update('fsc', fsc_data_pp, True)

        guiner = np.loadtxt(os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'guinierlines.txt')),
                            skiprows=1)
        if post_refine_options.mtf is not "":
            resol_sq = guiner[:, 0]
            log_amp_orig = guiner[:, 1]
            log_amp_weight = guiner[:, 2]
            log_amp_sharpened = guiner[:, 3]

            guiner_data_pp = pd.DataFrame(
                np.array([resol_sq, log_amp_orig, log_amp_weight, log_amp_sharpened]).swapaxes(0, 1).tolist(),
                columns=['_rlnResolutionSquared', '_rlnLogAmplitudesOriginal', '_rlnLogAmplitudesWeighted',
                         '_rlnLogAmplitudesSharpened'
                         ])
        else :
            resol_sq = guiner[:, 0]
            log_amp_orig = guiner[:, 1]
            log_amp_weight = guiner[:, 2]

            guiner_data_pp = pd.DataFrame(
                np.array([resol_sq, log_amp_orig, log_amp_weight]).swapaxes(0, 1).tolist(),
                columns=['_rlnResolutionSquared', '_rlnLogAmplitudesOriginal', '_rlnLogAmplitudesSharpened'
                         ])

        pp_star_file.update('guiner', guiner_data_pp, True)
        pp_star_file.write_star_file()

        ##### Post Process part ends here

        #################################################
        ############SPHIRE 2 RELION#################
        #################################################
        ###### SPHIRE 2 RELION Parts starts here
        meridien_folder = post_refine_options.combinemaps[0].split('/vol')[0]
        meridien_folder = os.path.join(os.getcwd(), meridien_folder)
        iter_files, bdbfile = get_final_iter_files(meridien_folder)
        ##Note sphire2relion requires the 3dparams file and the chunk

        if os.path.isdir(meridien_folder):
            sph2rel_call = (
                    "sp_sphire2relion.py"
                    + " " + os.path.join(os.getcwd(), os.path.join(str(options.Output_folder), "BDB2STAR"))
                    + " " + "--particle_stack=" + str(bdbfile)
                    + " " + "--params_3d_file=" + str(iter_files[0])
                    + " " + "--params_3d_chunk_file_0=" + str(iter_files[1])
                    + " " + "--params_3d_chunk_file_1=" + str(iter_files[2])
            )

            subprocess.run(args=[sph2rel_call], shell=True, text=True)


            #### This wont be necessary in cases where the entire pipeline was run with MotionCorr shifts
            bdb_star = star.StarFile(os.path.join(os.path.join(os.getcwd(), os.path.join(str(options.Output_folder), "BDB2STAR")),
                                                  'sphire2relion.star'))


def main():

    try:
        sp_global_def.print_timestamp("Start")
        run(sys.argv[1:])
    finally:
        sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()




####
import os

polishing_call = "relion_polishing_executable" \
                 + " --i " + os.path.join(os.getcwd(), os.path.join(str("Output_folder"),
                                                                    "BDB2STAR/sphire2relion.star")) \
                 + " " + "--f " + os.path.join(os.getcwd(), os.path.join(str("Output_folder"),
                                                                         "PostProcess/postprocess.star")) \
                 + " " + "--corr_mic " + os.path.join("final_motion_path") \
                 + " " + "--first_frame " + str("first_frame") \
                 + " " + "--last_frame " + str("last_frame") \
                 + " " + "--o " + str(os.path.join(os.getcwd(), str("Output_folder"))) \
                 + " " + "--params_file " + str("training_params") \
                 + " " + "--combine_frames" \
                 + " " + "--bfac_minfreq " + str("bfac_minfreq") \
                 + " " + "--bfac_maxfreq " + str("bfac_maxfreq") \
                 + " " + "--angpix_ref " + str("pixel_size[0]") \
                 + " " + "--j " + str("no_of_threads")

rel2sph_call = "\n\n" + "sp_relion2sphire.py" \
               + " " + os.path.join(os.getcwd(), os.path.join(str("Output_folder"), "shiny.star")) \
               + " " + "Polish_Stack" \
               + " " + "--relion_project_dir='.'" \
               + " " + "--box_size=-1"

polishing_rel2sphire_command = []

print(polishing_rel2sphire_command)

with open("test_appending_commands", "w") as w:
    w.write("".join(polishing_call + rel2sph_call))












