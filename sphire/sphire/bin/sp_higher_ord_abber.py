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


def parse_parameters(args):
    parser = argparse.ArgumentParser(args, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "post_refiner",
        type=str,
        help="post refiner directory",
    )

    parser.add_argument(
        "Output_folder",
        type=str,
        help="output_directory",
    )

    parser.add_argument(
        "--estimate_magnification",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--estimate_beamtilt",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--estimate_trefoil",
        action="store_true",
        help="If set to True, then relion_ctf_refine will also estimate the trefoil (3-fold astigmatism) per optics group."
             "This option is only recommended for data sets that extend beyond 3.5 Angstrom resolution",
        default=False
    )

    parser.add_argument(
        "--estimate_order_aberation",
        action="store_true",
        help="4th order aberration",
        default=False
    )

    parser.add_argument(
        "--perform_CTF_params_fit",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_defocus_micrograph",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_defocus_particle",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_astigmatism_micrograph",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_astigmatism_particle",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_bfactor_micrograph",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_bfactor_particle",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_phase_shift_micrograph",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--fit_phase_shift_particle",
        action="store_true",
        help="",
        default=False
    )

    parser.add_argument(
        "--min_res_fit",
        type=float,
        help="",
        default=30.0
    )

    parser.add_argument(
        "--submission_template",
        type = str,
        help ="",
        default = None
    )

    parser.add_argument(
        "--submission_command",
        type=str,
        help="",
        default="sbatch"
    )

    parser.add_argument(
        "--relion_mpirun_executable",
        type=str,
        help="",
        default="mpirun"
    )

    parser.add_argument(
        "--relion_ctfrefine_executable",
        type=str,
        help="",
        default="relion_ctf_refine_mpi"
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

    parser.add_argument(
        "--mrc_reloc_folder",
        type=str,
        help="",
        default= None
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
        default=None,
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
        default = None,
        nargs=1
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
        with open(os.path.join(options.post_refiner, 'command.txt'), 'r') as commandfile:
            read_command_txt = commandfile.read()

        post_refine_options = parse_postrefiner(read_command_txt)

        if post_refine_options.mask != None:
            use_mask = post_refine_options.mask[0]
        else:
            use_mask = os.path.join(options.post_refiner, 'vol_adaptive_mask.hdf')

        try:
            shutil.rmtree(os.path.join(str(options.Output_folder),'PostProcess'))
        except FileNotFoundError:
            os.makedirs( os.path.join(str(options.Output_folder),'PostProcess'))
            pass

        pp_star_file = star.StarFile(os.path.join( os.path.join(str(options.Output_folder),'PostProcess'),
                                                   'postprocess.star'))

        half_map1 = post_refine_options.combinemaps[0]
        half_map2 = post_refine_options.combinemaps[1]

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

        fsc_halves = np.loadtxt( os.path.join(options.post_refiner, 'halves.txt'))
        fsc_masked_halves = np.loadtxt(os.path.join(options.post_refiner, 'masked_halves.txt'))
        fsc_full = np.loadtxt( os.path.join(options.post_refiner, 'full.txt'))

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

        guiner = np.loadtxt(os.path.join(options.post_refiner, 'guinierlines.txt'), skiprows=1)
        if post_refine_options.mtf != None:
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
        iter_files, bdbfile = get_final_iter_files(meridien_folder)
        ##Note sphire2relion requires the 3dparams file and the chunk

        if os.path.isdir(meridien_folder):
            sph2rel_call = (
                    "sp_sphire2relion.py"
                    + " " + os.path.join(str(options.Output_folder), "BDB2STAR")
                    + " " + "--particle_stack=" + str(bdbfile)
                    + " " + "--params_3d_file=" + str(iter_files[0])
                    + " " + "--params_3d_chunk_file_0=" + str(iter_files[1])
                    + " " + "--params_3d_chunk_file_1=" + str(iter_files[2])
            )

            subprocess.run(args=[sph2rel_call], shell=True, text=True)

            if options.mrc_reloc_folder != None:
                bdb_star = star.StarFile(os.path.join( os.path.join(str(options.Output_folder), "BDB2STAR"),
                                 'sphire2relion.star'))
                old_micrograph_name = bdb_star[""]['_rlnMicrographName']
                new_micrograph_name = old_micrograph_name.apply(lambda x: os.path.join(options.mrc_reloc_folder
                                                                                       , os.path.basename(x)
                                                                                       ))
                bdb_star[""]['_rlnMicrographName'] = new_micrograph_name

                bdb_star.write_star_file(overwrite=True)


        """
        Till here it is generic script where we convert the stack and create the postprocess star file.
        """
        flag_defocus ='f'
        flag_astig = 'f'
        flag_bfactor = 'f'
        flag_phase_shift = 'f'
        total_flag = ""

        if options.estimate_magnification:
            total_flag += "--fit_aniso" + " " + "--kmin_mag " + str(options.min_res_fit)
        else:
            if options.perform_CTF_params_fit:
                if options.fit_defocus_micrograph:
                    flag_defocus = 'm'
                    total_flag += "--fit_defocus" + " " + "--kmin_defocus " + str(options.min_res_fit)
                elif options.fit_defocus_particle:
                    flag_defocus = 'p'
                    total_flag += "--fit_defocus" + " " + "--kmin_defocus " + str(options.min_res_fit)
                else:
                    pass

                if options.fit_astigmatism_micrograph:
                    flag_astig = 'm'
                elif options.fit_astigmatism_particle:
                    flag_astig = 'p'
                else :
                    pass

                if options.fit_bfactor_micrograph:
                    flag_bfactor = 'm'
                elif options.fit_bfactor_particle:
                    flag_bfactor = 'p'
                else :
                    pass

                if options.fit_phase_shift_micrograph:
                    flag_phase_shift = 'm'
                elif options.fit_phase_shift_particle:
                    flag_phase_shift = 'p'
                else :
                    pass

                fit_mode = flag_phase_shift + flag_defocus + flag_astig + 'f' + flag_bfactor
                total_flag += " " + "--fit_mode " + str(fit_mode)
            else :
                pass

            if options.estimate_beamtilt and not options.estimate_trefoil:
                total_flag += " " +  "--fit_beamtilt" + " " + "--kmin_tilt " + str(options.min_res_fit)
            elif options.estimate_beamtilt and options.estimate_trefoil :
                total_flag += " " + "--fit_beamtilt" + " " + "--kmin_tilt " + str(options.min_res_fit) + " " + "--odd_aberr_max_n 3"
            else :
                pass

            if options.estimate_order_aberation:
                total_flag += " " + "--fit_aberr"
            else :
                pass

        print("flags for the all commands is", total_flag)
        print("pixel size obtained is" , post_refine_options.pixel_size[0])
        time.sleep(5)

        ### now we need to decide we want to run it on a single PC workstation or on cluster
        if options.submission_template != None:
            ctfrefine_call = options.relion_ctfrefine_executable\
                             + " --i " + os.path.join(options.Output_folder,"BDB2STAR/sphire2relion.star")\
                             + " " + "--f " + os.path.join(str(options.Output_folder), "PostProcess/postprocess.star")\
                             + " " + "--o " + str(options.Output_folder) \
                             + " " + total_flag \
                             + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0]) \
                             + " " + "--j " + str(options.no_of_threads)

            try:
                with open(options.submission_template) as read:
                    lines = read.readlines()
            except Exception as e:
                sp_global_def.ERROR(str(e) + '\nCannot open mpi_submission template!', action=1)

            cmd_lines = []
            for idx, entry in enumerate(lines):
                if "XXX_SXCMD_LINE_XXX" in entry and "mpirun" in entry:
                    cmd_lines.append(idx)

            if not cmd_lines:
                sp_global_def.sxprint("Could not find a suitable command line for exchange.")
                sp_global_def.sxprint("The line should contain XXX_SXCMD_LINE_XXX.")
                sys.exit(1)

            line = (lines[cmd_lines[-1]].replace("XXX_SXCMD_LINE_XXX", ctfrefine_call ))

            mod_sub_script = "".join(lines).replace("XXX_SXMPI_NPROC_XXX", str(options.mpi_procs)
                                                    ).replace("XXX_SXMPI_JOB_NAME_XXX", "sp_higher_ord_abber"
                                                              ).replace(lines[cmd_lines[-1]], line
                                                                        ).replace("mpirun", options.relion_mpirun_executable)

            out_submission = "{0}/ctfrefine_submission_script.sh".format(str(options.Output_folder))
            with open(out_submission, "w") as w:
                w.write("".join(mod_sub_script))

            sp_global_def.sxprint(
                subprocess.check_output(
                    options.submission_command.split() + [out_submission]
                )
            )

        else:
            if options.mpi_procs > 1 :
                # if we want to run it on a workstation
                import mpi
                os.unsetenv('OMPI_COMM_WORLD_RANK')
                RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
                if RUNNING_UNDER_MPI:
                    mpi.mpi_init(0, [])
                    rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
                    size = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
                else:
                    rank = 0
                    size = 1

                env = os.environ
                new_env = {k: v for k, v in env.items() if "MPI" not in k}


                ctfrefine_call = (
                        "mpirun"
                        + " " + "-np"
                        + " " + str(options.mpi_procs)
                        + " " + options.relion_ctfrefine_executable
                        + " --i " + os.path.join(str(options.Output_folder), "BDB2STAR/sphire2relion.star")
                        + " " + "--f " + os.path.join(str(options.Output_folder), "PostProcess/postprocess.star")
                        + " " + "--o " + str(options.Output_folder)
                        + " " + total_flag
                        + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])
                        + " " + "--j " + str(options.no_of_threads)
                )
                # rel2sph_call = (
                #     "sp_relion2sphire.py"
                #     + " " +  os.path.join(str(options.Output_folder), "shiny.star")
                #     + " " + "Polish_Stack"
                #     + " " + "--relion_project_dir='.'"
                #     + " " + "--box_size=-1"
                # )

                print("Ctf refine with mpi command is called", ctfrefine_call)
                subprocess.run(args=[ctfrefine_call],  shell=True, text= True, env=new_env)
                # subprocess.run(args=[rel2sph_call], shell=True, text=True)
            else :
                ctfrefine_call = (
                        "relion_ctf_refine"
                        + " --i " + os.path.join(str(options.Output_folder), "BDB2STAR/sphire2relion.star")
                        + " " + "--f " + os.path.join(str(options.Output_folder), "PostProcess/postprocess.star")
                        + " " + "--o " +  str(options.Output_folder)
                        + " " + total_flag
                        + " " + "--angpix_ref " + str(post_refine_options.pixel_size[0])
                        + " " + "--j " + str(options.no_of_threads)
                )
                print("Ctf refine without mpi command is called", ctfrefine_call)
                subprocess.run(args=[ctfrefine_call], shell=True, text=True)



def main():

    try:
        sp_global_def.print_timestamp("Start")
        run(sys.argv[1:])
    finally:
        sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()

####

# from pyStarDB import sp_pystardb as star
# import os
# import pandas as pd
# given_file = os.path.join('/home/adnan/DemoResults/Relion_CTF_REFINE_Test_v1/BDB2STAR/',
#                                       'sphire2relion.star')
# outputfile = os.path.join('/home/adnan/DemoResults/Relion_CTF_REFINE_Test_v1/BDB2STAR/','test_newfile.star')
# starfile  = star.StarFile(given_file)
#
# a = [[0, 1], [2, 3]]
# dd = pd.DataFrame(a, columns=['_c1', '_c2'])
# starfile["newtag"] = dd
#
# starfile.write_star_file(outputfile, overwrite=True)












