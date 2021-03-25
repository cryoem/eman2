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
        "reference_map",
        type=str,
        help="reference map",
    )

    parser.add_argument(
        "--reference_mask",
        type=str,
        help="reference mask, optional",
        default = None
    )

    parser.add_argument(
        "--abs_greyscale_map",
        action="store_true",
        help="Ref map is on absolute greyscale",
        default = False
    )

    parser.add_argument(
        "--ini_high",
        type=float,
        help="initial low pass filter in angstorms",
        default = 60
    )

    parser.add_argument(
        "--sym",
        type=str,
        help="symmetry use for classification",
        default = "C1"
    )

    parser.add_argument(
        "--do_ctf",
        action="store_true",
        help="Do CTF correction",
        default = True
    )

    parser.add_argument(
        "--ctf_corr_ref",
        action="store_true",
        help="Has reference been CTF corrected",
        default = False,
    )

    parser.add_argument(
        "--ctf_ignore_peak",
        action="store_true",
        help="Ignore CTFs unit first peak",
        default = False,
    )

    parser.add_argument(
        "--no_of_class",
        type=int,
        help="Number of classes",
        default = 1
    )

    parser.add_argument(
        "--tau_val",
        type=float,
        help="Regularisation parameter T",
        default = 4
    )

    parser.add_argument(
        "--no_of_iter",
        type=int,
        help="Number of iterations",
        default = 25
    )

    parser.add_argument(
        "--use_fast_sets",
        action="store_true",
        help="Use fast subsets (for large datasets)",
        default = False,
    )

    parser.add_argument(
        "--diam_mas",
        type=int,
        help="Maske diameter in angstorm",
        default = 200
    )

    parser.add_argument(
        "--zeros_mas",
        action="store_true",
        help="Mask individual particles with zeros",
        default = True,
    )

    parser.add_argument(
        "--limit_resol_estep",
        type=int,
        help="Limit resolution E-step to (A)",
        default = -1
    )

    parser.add_argument(
        "--skip_img_align",
        action="store_true",
        help="Perform image alignment or not",
    )

    parser.add_argument(
        "--heal_pix_order",
        type=float,
        help="Angular sampling interval",
        default = 7.5
    )

    parser.add_argument(
        "--off_range",
        type=int,
        help="Offset search range (pix)",
        default = 5
    )

    parser.add_argument(
        "--off_step",
        type=int,
        help="Offset search step (pix)",
        default = 1
    )

    parser.add_argument(
        "--ang_search",
        action="store_true",
        help="Perform local angular searches",
        default = False,
    )

    parser.add_argument(
        "--ang_search_range",
        type=int,
        help="Local angular search",
        default = 5
    )

    parser.add_argument(
        "--ang_search_relax_sym",
        type=str,
        help="Relax Symmetry",
        default = None
    )

    parser.add_argument(
        "-coarse_sampling",
        action="store_true",
        help="Allow coarser sampling",
        default = False,
    )

    parser.add_argument(
        "-para_io",
        action="store_true",
        help="Use parallel disc I/O",
        default = True,
    )

    parser.add_argument(
        "--no_of_pool_part",
        type=int,
        help="Number of pooled particles",
        default = 3
    )

    parser.add_argument(
        "--skip_pad",
        action="store_true",
        help="Skip padding",
        default=False,
    )

    parser.add_argument(
        "--skip_grid",
        action="store_true",
        help="Skip padding",
        default=True,
    )

    parser.add_argument(
        "--pre_read_img",
        action="store_true",
        help="Pre-read all particles into RAM",
        default=False,
    )

    parser.add_argument(
        "--scratch_dir",
        type=str,
        help="Copy particles to search directory",
        default = None
    )

    parser.add_argument(
        "--combine_iter_disc",
        action="store_true",
        help="Combine iterations through disc",
        default=False,
    )

    parser.add_argument(
        "--use_gpu",
        action="store_true",
        help="Use GPU acceleration",
        default=False,
    )

    parser.add_argument(
        "--which_gpu",
        type=str,
        help="Which GPU to use",
        default = None
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
        "--relion_3dclassification_executable",
        type=str,
        help="",
        default="relion_refine_mpi"
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
        help="mrc relocation folder directory",
        default = None
    )

    ###### Here the helical part starts
    parser.add_argument(
        "--helical_recons",
        action="store_true",
        help="Do helical reconstruction",
        default = False
    )


    parser.add_argument(
        "--inner_diam",
        type=float,
        help="Helical Tube inner diameter",
        default= -1,
    )


    parser.add_argument(
        "--outer_diam",
        type=float,
        help="Helical Tube outer diameter",
        default= -1,
    )

    parser.add_argument(
        "--sig_tilt",
        type=float,
        help="Angular search range for tilt (deg)",
        default= 15,
    )

    parser.add_argument(
        "--sig_psi",
        type=float,
        help="Angular search range for psi (deg)",
        default= 10,
    )

    parser.add_argument(
        "--sig_rot",
        type=float,
        help="Angular search range for rot (deg)",
        default= -1,
    )

    parser.add_argument(
        "--sigma_dist",
        type=float,
        help="Range factor of local average",
        default= 1.5,
    )

    parser.add_argument(
        "--keep_tilt_fix",
        action="store_true",
        help="Keep tilt-prior fixed",
        default= True,
    )

    parser.add_argument(
        "--apply_helical_sym",
        action="store_true",
        help="Apply helical symmetry",
        default= False,
    )

    parser.add_argument(
        "--unique_asym_unit",
        type=int,
        help="Number of unique asymmetrical units",
        default= 1,
    )

    parser.add_argument(
        "--initial_twist",
        type=float,
        help="Initial twise (deg)",
        default= 0,
    )

    parser.add_argument(
        "--initial_rise",
        type=float,
        help="Initial rise (A)",
        default= 0,
    )

    parser.add_argument(
        "--z_percent",
        type=int,
        help="Central Z length (%)",
        default= 30,
    )

    parser.add_argument(
        "--do_local_search",
        action="store_true",
        help="Do local searches of symmetry",
        default= False,
    )

    parser.add_argument(
        "--twist_min",
        type=float,
        help="Twist search min (deg)",
        default= 0,
    )

    parser.add_argument(
        "--twist_max",
        type=float,
        help="Twist search max (deg)",
        default= 0,
    )

    parser.add_argument(
        "--twist_inistep",
        type=float,
        help="Twist search step (deg)",
        default= 0,
    )

    parser.add_argument(
        "--rise_min",
        type=float,
        help="Rise search min (A)",
        default= 0,
    )

    parser.add_argument(
        "--rise_max",
        type=float,
        help="Rise search max (A)",
        default= 0,
    )

    parser.add_argument(
        "--rise_inistep",
        type=float,
        help="Rise search max (A)",
        default= 0,
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
        default = None
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

    if os.path.isdir(options.post_refiner):
        with open(os.path.join(options.post_refiner, 'command.txt'), 'r') as commandfile:
            read_command_txt = commandfile.read()

        post_refine_options = parse_postrefiner(read_command_txt)

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
            bdb_star = star.StarFile(
                os.path.join(os.path.join(str(options.Output_folder), "BDB2STAR"), 'sphire2relion.star'))
            old_micrograph_name = bdb_star[""]['_rlnMicrographName']
            new_micrograph_name = old_micrograph_name.apply(lambda x: os.path.join(options.mrc_reloc_folder
                                                                                   , os.path.basename(x)
                                                                                   ))
            bdb_star[""]['_rlnMicrographName'] = new_micrograph_name
            bdb_star.write_star_file(overwrite=True)


        total_str = ""
        if options.reference_mask != None :
            total_str += " " + "--solvent_mask" + " " + str(options.reference_mask)
        else:
            pass

        ##Now reference part starts
        if options.ini_high:
            total_str += " " + "--ini_high" + " " + options.ini_high
        else:
            pass
        if options.sym :
            total_str += " " + "--sym" + " " + options.sym
        else:
            pass
        if options.abs_greyscale_map :
            pass
        else:
            total_str += " " + "--firstiter_cc"

        #Now sCTF part starts
        if options.do_ctf:
            total_str += " " + "--ctf"
            if options.ctf_corr_ref:
                total_str += " " + "--ctf_corrected_ref"
            if options.ctf_ignore_peak :
                total_str += " " + "--ctf_intact_first_peak"
        else:
            pass

        ### Now Optimization part starts
        if options.no_of_class:
            total_str += " " + "--K" + " " + options.no_of_class
        else:
            pass
        if options.tau_val:
            total_str += " " + "--tau2_fudge" + " " + options.tau_val
        else:
            pass
        if options.no_of_iter:
            total_str += " " + "--iter" + " " + options.no_of_iter
        else:
            pass
        if options.use_fast_sets:
            total_str += " " + "--fast_subsets"
        else:
            pass
        if options.diam_mas:
            total_str += " " + "--particle_diameter" + " " + options.diam_mas
        else:
            pass
        if options.zeros_mas:
            total_str += " " +  "--zero_mask"
        else:
            pass
        if options.limit_resol_estep:
            total_str += " " + "--strict_highres_exp" + " " + options.limit_resol_estep
        else:
            pass

        ### Now Sampling part starts
        if options.skip_img_align :
            total_str += " " + "--skip_align" + " " + options.skip_img_align
        else:
            if options.heal_pix_order :
                sample_value = 2
                if options.heal_pix_order == 0.1 :
                    sample_value = 8
                elif options.heal_pix_order == 0.2 :
                    sample_value = 7
                elif options.heal_pix_order == 0.5 :
                    sample_value = 6
                elif options.heal_pix_order == 0.9 :
                    sample_value = 5
                elif options.heal_pix_order == 1.8 :
                    sample_value = 4
                elif options.heal_pix_order == 3.7 :
                    sample_value = 3
                elif options.heal_pix_order == 7.5 :
                    sample_value = 2
                elif options.heal_pix_order == 15 :
                    sample_value = 1
                elif options.heal_pix_order == 30 :
                    sample_value = 0
                else :
                    print("Please specify the right value from the drop down menu, the value can be either \
                          0.1, 0.2, 0.5, 0.9, 1.8, 3.7, 7.5, 15, 30")
                total_str += " " + "--healpix_order" + " " + sample_value
            else:
                pass
            if options.off_range :
                total_str += " " + "--offset_range" + " " + options.off_range
            else:
                pass
            if options.off_step :
                total_str += " " + "--offset_step" + " " + options.off_step*2
            else:
                pass
            if options.ang_search :
                total_str += " " + "--sigma_ang" + " " + options.ang_search_range*0.3333
                if ang_search_relax_sym :
                    total_str += " " + "--relax_sym" + " " + str(ang_search_relax_sym)
                else :
                    pass
            else:
                pass
        if options.coarse_sampling :
            total_str += " " + "--allow_coarser_sampling"
        else:
            pass

        ### Now Compute part starts
        if options.para_io :
            pass
        else:
            total_str += " " + "--no_parallel_disc_io"
        if options.no_of_pool_part :
            total_str += " " + "--pool" + " " + options.no_of_pool_part
        else:
            pass
        if options.skip_pad :
            total_str += " " + "--pad" + " " + 1
        else:
            total_str += " " + "--pad" + " " + 2
        if options.skip_grid :
            total_str += " " + "--skip_gridding"
        else:
            pass
        if options.pre_read_img :
            total_str += " " + "--preread_images"
        else:
            pass
        if options.scratch_dir :
            total_str += " " + "--scratch_dir" + " " + str(options.scratch_dir)
        else:
            pass
        if options.combine_iter_disc :
            pass
        else:
            total_str += " " + "--dont_combine_weights_via_disc"

        if options.use_gpu :
            if options.which_gpu !=None :
                total_str += " " + "--gpu" + " " + str(options.which_gpu)
        else:
            pass

        ### Adding unknown flags which i can still not debugg
        total_str += " " +"--flatten_solvent --oversampling 1 --norm  --scale"

        ### Now helical parts start
        if options.helical_recons :
            total_str += " " + "--helix"
            total_str += " " + "--helical_outer_diameter" + " " + str(options.outer_diam)
            if options.inner_diam != -1:
                total_str += " " + "--helical_inner_diameter" + " " + str(options.inner_diam)
            else :
                pass
            if options.sig_tilt >= 3 :
                total_str += " " + "--sigma_tilt" + " " + str(options.sig_tilt / 3)
            else :
                total_str += " " + "--sigma_tilt" + " " + str(0)
            if options.sig_rot >= 3 :
                total_str += " " + "--sigma_rot" + " " + str(options.sig_rot / 3)
            else :
                total_str += " " + "--sigma_rot" + " " + str(0)
            if options.sig_psi >= 3 :
                total_str += " " + "--sigma_psi" + " " + str(options.sig_psi / 3)
            else :
                total_str += " " + "--sigma_rot" + " " + str(0)
            if options.sigma_dist >=1 :
                total_str += " " + "--helical_sigma_distance" + " " + str(options.sigma_dist / 3)
            else :
                pass
            if options.keep_tilt_fix :
                total_str += " " + "--helical_keep_tilt_prior_fixed"
            else :
                pass

            if options.apply_helical_sym :

                total_str += " " + "--helical_nr_asu" + " " + str(options.unique_asym_unit)
                total_str += " " + "--helical_twist_initial" + " " + str(options.initial_twist)
                total_str += " " + "--helical_rise_initial" + " " + str(options.initial_rise)
                total_str += " " + "--helical_z_percentage" + " " + str(options.z_percent/100)
            else :
                pass

            if options.do_local_search :
                total_str += " " + "--helical_symmetry_search"
                total_str += " " + "--helical_twist_min" + " " + str(options.twist_min)
                total_str += " " + "--helical_twist_max" + " " + str(options.twist_max)
                total_str += " " + "--helical_twist_inistep" + " " + str(options.twist_inistep)
                total_str += " " + "--helical_rise_min" + " " + str(options.rise_min)
                total_str += " " + "--helical_rise_max" + " " + str(options.rise_max)
                if options.rise_inistep != 0 :
                    total_str += " " + "--helical_rise_inistep" + " " + str(options.rise_inistep)
                else :
                    pass
            else:
                pass
        else:
            pass


        print("flags for the all commands is", total_str)
        print("pixel size obtained is" , post_refine_options.pixel_size[0])
        time.sleep(5)

        ### now we need to decide we want to run it on a single PC workstation or on cluster
        if options.submission_template != None:
            classification_call = options.relion_3dclassification_executable \
                             + " " + "--o " + str(options.Output_folder) \
                             + " --i " + os.path.join(options.Output_folder, "BDB2STAR/sphire2relion.star") \
                             + " " + "--ref " + str(options.reference_map) \
                             + " " + total_str \
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

            line = (lines[cmd_lines[-1]].replace("XXX_SXCMD_LINE_XXX", classification_call))

            mod_sub_script = "".join(lines).replace("XXX_SXMPI_NPROC_XXX", str(options.mpi_procs)
                                                    ).replace("XXX_SXMPI_JOB_NAME_XXX", "sp_relion_3dclassifi"
                                                              ).replace(lines[cmd_lines[-1]], line
                                                                        ).replace("mpirun",
                                                                                  options.relion_mpirun_executable)

            out_submission = "{0}/classification3d_submission_script.sh".format(str(options.Output_folder))
            with open(out_submission, "w") as w:
                w.write("".join(mod_sub_script))

            sp_global_def.sxprint(
                subprocess.check_output(
                    options.submission_command.split() + [out_submission]
                )
            )

        else:
            if options.mpi_procs > 1:
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

                classification_call = (
                        "mpirun"
                        + " " + "-np"
                        + " " + str(options.mpi_procs)
                        + " " + options.relion_3dclassification_executable
                        + " " + "--o " + str(options.Output_folder)
                        + " --i " + os.path.join(str(options.Output_folder), "BDB2STAR/sphire2relion.star")
                        + " " + "--ref " + str(options.reference_map)
                        + " " + total_str
                        + " " + "--j " + str(options.no_of_threads)
                )
                # rel2sph_call = (
                #     "sp_relion2sphire.py"
                #     + " " +  os.path.join(str(options.Output_folder), "shiny.star")
                #     + " " + "Polish_Stack"
                #     + " " + "--relion_project_dir='.'"
                #     + " " + "--box_size=-1"
                # )

                print("3d classification with mpi command is called", classification_call)
                subprocess.run(args=[classification_call], shell=True, text=True, env=new_env)
                # subprocess.run(args=[rel2sph_call], shell=True, text=True)
            else:
                classification_call = (
                        options.relion_3dclassification_executable
                        + " " + "--o " + str(options.Output_folder)
                        + " --i " + os.path.join(str(options.Output_folder), "BDB2STAR/sphire2relion.star")
                        + " " + "--ref " + str(options.reference_map)
                        + " " + total_flag
                        + " " + "--j " + str(options.no_of_threads)
                )
                print("3dclasification without mpi command is called", classification_call)
                subprocess.run(args=[classification_call], shell=True, text=True)


###

def main():

    try:
        sp_global_def.print_timestamp("Start")
        run(sys.argv[1:])
    finally:
        sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()