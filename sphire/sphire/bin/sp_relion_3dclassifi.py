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
        "refernce_map",
        type=str,
        help="reference map",
    )

    parser.add_argument(
        "--refernce_mask",
        type=str,
        help="reference mask, optional",
    )

    parser.add_argument(
        "--ini_high",
        type=float,
        help="initial low pass filter in angstorms",
    )

    parser.add_argument(
        "--sym",
        type=str,
        help="symmetry use for classification",
    )

    parser.add_argument(
        "--do_ctf",
        action="store_true",
        help="Do CTF correction",
        default = False
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
    )

    parser.add_argument(
        "--tau_val",
        type=float,
        help="Regularisation parameter T",
    )

    parser.add_argument(
        "--no_of_iter",
        type=int,
        help="Number of iterations",
    )

    parser.add_argument(
        "--use_fast_sets",
        action="store_true",
        help="Use fast subsets (for large datasets)",
        default = False,
    )

    parser.add_argument(
        "--mask_diam",
        type=int,
        help="Maske diameter in angstorm",
    )

    parser.add_argument(
        "--mask_zeros",
        action="store_true",
        help="Mask individual particles with zeros",
        default = False,
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
        default = ""
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
        default = False,
    )

    parser.add_argument(
        "--no_of_pool_part",
        type=int,
        help="Number of pooled particles",
        default = 5
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
        default=False,
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
        default = ""
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
        default = ""
    )

    parser.add_argument(
        "--mrc_reloc_folder",
        type=str,
        help="mrc relocation folder directory",
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

    if os.path.isdir(options.post_refiner):
        with open(os.path.join(os.getcwd(), os.path.join(options.post_refiner, 'command.txt')), 'r') as commandfile:
            read_command_txt = commandfile.read()

        post_refine_options = parse_postrefiner(read_command_txt)

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

        if options.mrc_reloc_folder is not "":
            bdb_star = star.StarFile(
                os.path.join(os.path.join(os.getcwd(), os.path.join(str(options.Output_folder), "BDB2STAR")),
                             'sphire2relion.star'))
            old_micrograph_name = bdb_star[""]['_rlnMicrographName']
            new_micrograph_name = old_micrograph_name.apply(lambda x: os.path.join(options.mrc_reloc_folder
                                                                                   , os.path.basename(x)
                                                                                   ))
            bdb_star[""]['_rlnMicrographName'] = new_micrograph_name

            bdb_star.write_star_file(overwrite=True)


def main():

    try:
        sp_global_def.print_timestamp("Start")
        run(sys.argv[1:])
    finally:
        sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()