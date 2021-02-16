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
        "--training_params",
        type=str,
        help="",
        default=None
    )

    parser.add_argument(
        "--first_frame",
        type=int,
        help="",
        default = 1,
    )

    parser.add_argument(
        "--last_frame",
        type=int,
        help="",
        default=-1

    )
    parser.add_argument(
        "--bfac_minfreq",
        type=int,
        help="",
        default=20
    )

    parser.add_argument(
        "--bfac_maxfreq",
        type=int,
        help="",
        default=-1
    )

    parser.add_argument(
        "--min_no_particles",
        type=int,
        help="",
        default=5000
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

            """
            old_micrograph_name = bdb_star[""]['_rlnMicrographName']
            newloc = os.path.join(os.path.dirname(options.corr_mic), 'Movies/')
            new_micrograph_name = old_micrograph_name.apply(lambda x: os.path.join(
                                                            os.path.join(os.getcwd(),
                                                                                   newloc)
                                                                                  , os.path.basename(x)
                                                                                   ))
            bdb_star[""]['_rlnMicrographName'] = new_micrograph_name

            bdb_star.write_star_file(overwrite=True)

            """

        ####### SPHIRE 2 RELION Parts ends here

        #################################################
        ############Corrected micrograph#################
        #################################################
        #######Corrected micrograph part starts here
        main_location = os.getcwd()
        # This will be later set in case we have motion corr data
        final_motion_path = os.path.join(main_location, options.corr_mic)
        time.sleep(5)



        if options.training_params != None:
            polishing_call = "/mnt/beegfs/software/em/relion/relion-3.1.0_gcc_7.2.0/bin/relion_motion_refine_mpi"\
                             + " --i " + os.path.join(os.getcwd(),os.path.join(str(options.Output_folder),
                                                                               "BDB2STAR/sphire2relion.star"))\
                             + " " + "--f " + os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),
                                                                            "PostProcess/postprocess.star"))\
                             + " " + "--corr_mic " + os.path.join(final_motion_path)\
                             + " " + "--first_frame " + str(options.first_frame)\
                             + " " + "--last_frame " + str(options.last_frame)\
                             + " " + "--o " + str(os.path.join(os.getcwd(), str(options.Output_folder)))\
                             + " " + "--params_file " + str(options.training_params)\
                             + " " + "--combine_frames"\
                             + " " + "--bfac_minfreq " + str(options.bfac_minfreq)\
                             + " " + "--bfac_maxfreq " + str(options.bfac_maxfreq)\
                             + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])\
                             + " " + "--j " + str(options.no_of_threads)
            try:
                with open(options.submission_template) as read:
                    lines = read.readlines()
            except Exception as e:
                sp_global_def.ERROR(str(e) + '\nCannot open mpi_submission template!', action=1)

            sp_global_def.sxprint("Lines are", lines)
            sp_global_def.sxprint("options.submission_template is ", options.submission_template)

            cmd_lines = []
            for idx, entry in enumerate(lines):
                sp_global_def.sxprint("entry values are ", entry)
                if "XXX_SXCMD_LINE_XXX" in entry and "mpirun" in entry:
                    cmd_lines.append(idx)

            if not cmd_lines:
                sp_global_def.sxprint("Could not find a suitable command line for exchange.")
                sp_global_def.sxprint("The line should contain XXX_SXCMD_LINE_XXX.")
                sys.exit(1)

            line = (lines[cmd_lines[-1]].replace("XXX_SXCMD_LINE_XXX", polishing_call))

            mod_sub_script = "".join(lines).replace("XXX_SXMPI_NPROC_XXX", str(options.mpi_procs)
                                                    ).replace("XXX_SXMPI_JOB_NAME_XXX", "sp_polishing"
                                                              ).replace(lines[cmd_lines[-1]], line
                                                                        ).replace("mpirun", "/mnt/beegfs/software/em/openmpi/openmpi-2.0.4_gcc_7.2.0_cuda_10.0/bin/mpirun")




            out_submission = "{0}/polishing_submission_script.sh".format(str(options.Output_folder))
            with open(out_submission, "w") as w:
                w.write("".join(mod_sub_script))

            sp_global_def.sxprint(
                subprocess.check_output(
                    options.submission_command.split() + [out_submission]
                )
            )

        else:
            print("Parameter file not provided, hence training is performed")
            polishing_call = (
                    "relion_motion_refine"
                    + " --i "
                    + os.path.join(os.getcwd(),os.path.join(str(options.Output_folder), "BDB2STAR/sphire2relion.star") )
                    + " " + "--f " + os.path.join(os.getcwd(),os.path.join(str(options.Output_folder),"PostProcess/postprocess.star"))
                    + " " + "--corr_mic " + os.path.join(os.getcwd(), os.path.join(final_motion_path))
                    + " " + "--first_frame " + str(options.first_frame)
                    + " " + "--last_frame " + str(options.last_frame)
                    + " " + "--o " + str(options.Output_folder)
                    + " " + "--min_p " + str(options.min_no_particles)
                    + " " + "--eval_frac " + str(0.5)
                    + " " + "--align_frac " + str(0.5)
                    + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])
                    + " " + "--params3"
                    + " " + "--j " + str(options.no_of_threads)
                    + " " + "--pipeline_control " + str(options.Output_folder)
            )
            subprocess.run(args=[polishing_call], shell=True, text=True)



            #     #######Corrected micrograph part end here
            # ####  For applying polishing on the particles and estimating bfactor
            # if options.training_params != None:
            #     if options.mpi_procs > 1:
            #         polishing_call = (
            #                 "mpirun"
            #                 + " " + "-np"
            #                 + " " + str(options.mpi_procs)
            #                 + " " + "relion_motion_refine_mpi"
            #                 + " --i "
            #                 + os.path.join(os.getcwd(),os.path.join(str(options.Output_folder),"BDB2STAR/sphire2relion.star"))
            #                 + " " + "--f " + os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),"PostProcess/postprocess.star"))
            #                 + " " + "--corr_mic " + os.path.join(final_motion_path)
            #                 + " " + "--first_frame " + str(options.first_frame)
            #                 + " " + "--last_frame " + str(options.last_frame)
            #                 + " " + "--o " + str(os.path.join(os.getcwd(), str(options.Output_folder)))
            #                 + " " + "--params_file " + str(options.training_params)
            #                 + " " + "--combine_frames"
            #                 + " " + "--bfac_minfreq " + str(options.bfac_minfreq)
            #                 + " " + "--bfac_maxfreq " + str(options.bfac_maxfreq)
            #                 + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])
            #                 + " " + "--j " + str(options.no_of_threads)
            #         )
            #     else :
            #         polishing_call = (
            #                 "relion_motion_refine"
            #                 + " --i "
            #                 + os.path.join(os.getcwd(),os.path.join(str(options.Output_folder),"BDB2STAR/sphire2relion.star"))
            #                 + " " + "--f " + os.path.join(os.getcwd(), os.path.join(str(options.Output_folder),"PostProcess/postprocess.star"))
            #                 + " " + "--corr_mic " + os.path.join(final_motion_path)
            #                 + " " + "--first_frame " + str(options.first_frame)
            #                 + " " + "--last_frame " + str(options.last_frame)
            #                 + " " + "--o " + str(os.path.join(os.getcwd(), str(options.Output_folder)))
            #                 + " " + "--params_file " + str(options.training_params)
            #                 + " " + "--combine_frames"
            #                 + " " + "--bfac_minfreq " + str(options.bfac_minfreq)
            #                 + " " + "--bfac_maxfreq " + str(options.bfac_maxfreq)
            #                 + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])
            #                 + " " + "--j " + str(options.no_of_threads)
            #         )
            #     rel2sph_call = (
            #         "sp_relion2sphire.py"
            #         + " " + os.path.join(os.getcwd(), os.path.join(str(options.Output_folder), "shiny.star"))
            #         + " " + "Polish_Stack"
            #         + " " + "--relion_project_dir='.'"
            #         + " " + "--box_size=-1"
            #     )

                # print("Polishing command is called", polishing_call)
                # subprocess.run(args=[polishing_call],  shell=True, text= True, env=new_env)
                # subprocess.run(args=[rel2sph_call], shell=True, text=True)
            ### For running the training part of polishing

##
def main():

    try:
        sp_global_def.print_timestamp("Start")
        # if rank == 0:
        run(sys.argv[1:])

        # mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    finally:
        # mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        # mpi.mpi_finalize()
        sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()

"""
corr_sums_log = os.path.join(main_location, 'CorrectedSums/corrsum_dw_log/')
corr_sums_mrc = os.path.join(main_location, 'CorrectedSums/corrsum_dw/')

corr_sum_star = os.path.join(main_location, options.corr_mic)
logfiles =  glob('{}*.star'.format(os.path.join(os.path.dirname(corr_sum_star), 'Movies/')))
if logfiles == []:
    logfiles = glob('{}*.log'.format(corr_sums_log))
    mrcfiles = glob('{}*.mrc'.format(corr_sums_mrc))
else:
    mrcfiles = glob('{}*.mrc'.format(os.path.join(os.path.dirname(corr_sum_star), 'Movies/')))

if logfiles[0].endswith('.star'):
    final_motion_path = corr_sum_star
else:
    ### there is no guarantee that conversion from unblur data will work so better to use MotionCorr
    final_motion_path = os.path.join(os.path.join(os.getcwd(), 'MotionCorr/corrected_micrographs.star'))
    # try:
    #     shutil.rmtree(os.path.join(os.getcwd(), 'MotionCorr'))
    # except FileNotFoundError:
    #     os.makedirs(os.path.join(os.getcwd(), 'MotionCorr'))
    #     pass

    early_data = []
    late_data = []
    total_data = []

    for test_logfile in logfiles:

        new_log_name = os.path.basename(test_logfile).split('.log')[0]
        pp_star_motion = star.StarFile(os.path.join(os.path.join(os.getcwd(), 'MotionCorr'),
                                                    '{}.star'.format(new_log_name)))

        with open(test_logfile, 'r') as logfile:
            all_lines = logfile.readlines()
        shifts = []
        for line in all_lines:
            if line[:5] == 'image':
                shifts.append([float(item) for item in line.split('=')[1:][0].split(',')])
            elif line[:20] == 'Input stack filename':
                motion_Input_stack = line.split(':')[1].replace(' ', '').replace('\n', '')
                motion_Input_stack = os.path.join(os.getcwd(), motion_Input_stack)
            elif line[:14] == 'Output binning':
                motion_output_bin = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line[:5] == 'Pixel':
                motion_pixel_size = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line[:8] == 'Exposure':
                motion_exposure = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line[:12] == 'Pre-exposure':
                motion_pre_exposure = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line[:12] == 'Acceleration':
                motion_voltage = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line[:5] == 'First':
                motion_first_frame = line.split(':')[1].replace(' ', '').replace('\n', '')

        stack_metadata = EMAN2db.db_get_image_info(motion_Input_stack)
        xsize = stack_metadata[1][0]
        ysize = stack_metadata[1][1]
        zsize = stack_metadata[1][2]
        motion_meta_data_pp = pd.DataFrame([[xsize, ysize, zsize, motion_Input_stack, motion_output_bin,
                                             motion_pixel_size, motion_exposure,
                                             motion_pre_exposure, motion_voltage, motion_first_frame, 0]],
                                           columns=['_rlnImageSizeX', '_rlnImageSizeY', '_rlnImageSizeZ',
                                                    '_rlnMicrographMovieName', '_rlnMicrographBinning',
                                                    '_rlnMicrographOriginalPixelSize', '_rlnMicrographDoseRate',
                                                    '_rlnMicrographPreExposure', '_rlnVoltage',
                                                    '_rlnMicrographStartFrame', '_rlnMotionModelVersion'])

        pp_star_motion.update('general', motion_meta_data_pp, False)

        motion_data_pp = pd.DataFrame()

        motion_data_pp['_rlnMicrographFrameNumber'] = np.arange(1, len(shifts) + 1)
        motion_data_pp['_rlnMicrographShiftX'] = np.array(shifts)[:, 0] - np.array(shifts)[0, 0]
        motion_data_pp['_rlnMicrographShiftY'] = np.array(shifts)[:, 1] - np.array(shifts)[0, 1]

        pp_star_motion.update('global_shift', motion_data_pp, True)
        pp_star_motion.write_star_file()

        ### because dose_rate = 2.5 ,.... number of frames for early = dose_rate * n < 4 where n is the number of frames

        cutoff_frame = (4 - float(motion_pre_exposure)) / float(motion_exposure)
        cutoff = 0
        for i in range(2, len(shifts)):
            if i <= int(cutoff_frame):
                cuttoff = i
            else:
                pass

        total_diff = np.array(shifts)[1:] - np.array(shifts)[:-1]
        total = np.sum(np.sqrt(np.inner(total_diff, total_diff).diagonal()))
        try:
            if cutoff == 0:
                early = 0
            else:
                early = np.sum(np.sqrt(np.inner(total_diff[:cutoff], total_diff[:cutoff]).diagonal()))
        except ValueError:
            if cutoff == 0:
                early = 0
            else:
                early = np.sum(np.sqrt(np.inner(total_diff[:cutoff], total_diff[:cutoff])))

        late = total - early
        early_data.append(early)
        late_data.append(late)
        total_data.append(total)

    newlogfiles = glob('{}*.star'.format(os.path.join(os.getcwd(), 'MotionCorr/')))

    corr_micro_star = star.StarFile(os.path.join(os.path.join(os.getcwd(), 'MotionCorr'),
                                                 'corrected_micrographs.star'))
    micrograph_files = pd.DataFrame(np.array([mrcfiles, newlogfiles, total_data, early_data,
                                              late_data]).swapaxes(0, 1).tolist(),
                                    columns=['_rlnMicrographName', '_rlnMicrographMetadata',
                                             '_rlnAccumMotionTotal', '_rlnAccumMotionEarly',
                                             '_rlnAccumMotionLate'])
    corr_micro_star.update('', micrograph_files, True)
    corr_micro_star.write_star_file()
"""
