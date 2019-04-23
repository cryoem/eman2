#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: Markus Stabrin 2019/03/22 (markus.stabrin@mpi-dortmund.mpg.de)
# Copyright (c) 2019 MPI Dortmund
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#


import re
import argparse
import collections as co


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=9999))

    group = parser.add_argument_group('Global settings (required)')
    group.add_argument('--apix', type=float, default=1.0, help='Pixel size in A/pixel.')
    group.add_argument('--mol_mass', type=float, default=250.0, help='Molecular mass of the protein in kDa. Used to calculate the masking density threshold.')
    group.add_argument('--particle_radius', type=int, default=80, help='Particle radius in pixels. Used for normalization.')
    group.add_argument('--box_size', type=int, default=200, help='Particle box size in pixels.')
    group.add_argument('--symmetry', type=str, default='c1', help='Symmetry of the particle.')
    group.add_argument('--voltage', type=float, default=300.0, help='Microscope voltage in kV.')
    group.add_argument('--negative_stain', action='store_true', default=False, help='Input is negative stain.')
    group.add_argument('--phase_plate', action='store_true', default=False, help='Input is phase_plate.')

    group = parser.add_argument_group('Unblur settings (required to run movie alignment)')
    group.add_argument('--unblur_path', type=str, default='/Unblur/executable', help='Path pointing to the unblur executable.')
    group.add_argument('--unblur_mic_pattern', type=str, default='Movies/*.mrc', help='Pattern of the micrographs to use for motion correction.')
    group.add_argument('--unblur_exp_per_frame', type=float, default=2.5, help='Exposure per frame. Used for dose adjustment.')
    group.add_argument('--unblur_gain_file', type=str, default='/Gain/file.mrc', help='File containing the information for gain correction. Not required if the movies are already gain corrected.')

    group = parser.add_argument_group('Unblur settings (optional)')
    group.add_argument('--skip_unblur', action='store_true', default=False, help='Do not run motion correction')
    group.add_argument('--unblur_output_dir', type=str, default='00a_UNBLUR', help='Unblur output directory.')
    group.add_argument('--unblur_first_frame', type=int, default=0, help='First frame to use for motion correction.')
    group.add_argument('--unblur_last_frame', type=int, default=-1, help='Last frame to use for motion correction. -1 means last frame.')
    group.add_argument('--unblur_addition', type=str, default='', help='Additional parameters that are not part of the required ones.')

    group = parser.add_argument_group('CTER settings (required to run CTF estimation)')
    group.add_argument('--cter_cs', type=float, default=2.7, help='Spherical aberration of the microscope.')

    group = parser.add_argument_group('CTER settings (optional)')
    group.add_argument('--skip_cter', action='store_true', default=False, help='Do not run CTF estimation.')
    group.add_argument('--cter_output_dir', type=str, default='01a_CTER', help='CTER output directory.')
    group.add_argument('--cter_mic_pattern', type=str, default='Mics/*.mrc', help='Micrograph pattern in case unblur is skipped.')
    group.add_argument('--cter_window_size', type=int, default=1024, help='CTF estimation window size.')
    group.add_argument('--cter_addition', type=str, default='', help='Additional parameters that are not part of the required ones.')

    group = parser.add_argument_group('CRYOLO settings (required to run particle picking)')
    group.add_argument('--cryolo_predict_path', type=str, default='/Path/cryolo_predict.py', help='Path to the cryolo predict executable.')
    group.add_argument('--cryolo_config_path', type=str, default='/Path/config', help='Path to the cryolo config Path')
    group.add_argument('--cryolo_gpu', type=str, default='0', help='Cryolo GPU list.')

    group = parser.add_argument_group('CRYOLO settings (optional)')
    group.add_argument('--skip_cryolo', action='store_true', default=False, help='Do not run particle picking.')
    group.add_argument('--cryolo_output_dir', type=str, default='02a_CRYOLO_PREDICT', help='CRYOLO output directory.')
    group.add_argument('--cryolo_mic_pattern', type=str, default='Mics/*.mrc', help='Micrograph pattern in case unblur is skipped.')
    group.add_argument('--cryolo_addition', type=str, default='', help='Additional parameters that are not part of the required ones.')

    group = parser.add_argument_group('WINDOW settings (optional)')
    group.add_argument('--skip_window', action='store_true', default=False, help='Do not run particle extraction.')
    group.add_argument('--window_box_pattern', type=str, default='Boxes/*.box', help='Window box file pattern.')
    group.add_argument('--window_mic_pattern', type=str, default='Mics/*.mrc', help='Window mrc file pattern.')
    group.add_argument('--window_partres', type=str, default='CTER/partres.txt', help='CTER partres file. In case of negative stain put this value to the pixel size.')
    group.add_argument('--window_output_dir', type=str, default='02b_WINDOW', help='Window output directory.')
    group.add_argument('--window_addition', type=str, default='', help='Additional parameters that are not part of the required ones.')

    return parser.parse_args()


def get_unblur_cmd(status_dict, **kwargs):
    cmd = []
    cmd.append('sp_unblur.py')
    cmd.append('XXX_SP_UNBLUR_PATH_XXX')
    cmd.append('XXX_SP_UNBLUR_MICROGRAPH_PATTERN_XXX')
    cmd.append('XXX_SP_UNBLUR_OUTPUT_DIR_XXX')
    cmd.append('--pixel_size=XXX_SP_PIXEL_SIZE_XXX')
    cmd.append('--voltage=XXX_SP_VOLTAGE_XXX')
    cmd.append('--exposure_per_frame=XXX_SP_UNBLUR_EXP_PER_FRAME_XXX')
    cmd.append('--gain_file=XXX_SP_UNBLUR_GAIN_FILE_XXX')
    cmd.append('--first_frame=XXX_SP_UNBLUR_FIRST_FRAME_XXX')
    cmd.append('--last_frame=XXX_SP_UNBLUR_LAST_FRAME_XXX')
    cmd.append('XXX_SP_UNBLUR_ADDITION_XXX')
    return cmd


def get_cter_cmd(status_dict, phase_plate, **kwargs):
    cmd = []
    cmd.append('sp_cter.py')
    if status_dict['do_unblur']:
        cmd.append('XXX_SP_UNBLUR_OUTPUT_DIR_XXX/corrsum_dw/*.mrc')
    else:
        cmd.append('XXX_SP_CTER_MICROGRAPH_PATTERN_XXX')
    cmd.append('XXX_SP_CTER_OUTPUT_DIR_XXX')
    cmd.append('--wn=XXX_SP_CTER_WINDOW_SIZE')
    cmd.append('--apix=XXX_SP_PIXEL_SIZE_XXX')
    cmd.append('--Cs=XXX_SP_CTER_CS_XXX')
    cmd.append('--voltage=XXX_SP_VOLTAGE_XXX')
    cmd.append('--pws')
    if phase_plate:
        cmd.append('--vpp')
    cmd.append('XXX_SP_CTER_ADDITION_XXX')
    return cmd


def get_cryolo_predict(status_dict, **kwargs):
    cmd = []
    cmd.append('sp_cryolo_predict.py')
    cmd.append('XXX_SP_CRYOLO_CONFIG_PATH_XXX')
    if status_dict['do_unblur']:
        cmd.append('XXX_SP_UNBLUR_OUTPUT_DIR_XXX/corrsum_dw')
    elif status_dict['do_cter']:
        cmd.append('XXX_SP_CTER_MICROGRAPH_PATTERN_XXX')
    else:
        cmd.append('XXX_SP_CRYOLO_MICROGRAPH_PATTERN_XXX')
    cmd.append('XXX_SP_CRYOLO_OUTPUT_DIR_XXX')
    cmd.append('XXX_SP_CRYOLO_PREDICT_PATH_XXX')
    cmd.append('--gpu=XXX_SP_CRYOLO_GPU_XXX')
    cmd.append('XXX_SP_CRYOLO_ADDITION_XXX')
    return cmd


def get_window(status_dict, negative_stain, **kwargs):
    cmd = []
    cmd.append('sp_window.py')
    if status_dict['do_unblur']:
        cmd.append('XXX_SP_UNBLUR_OUTPUT_DIR_XXX/corrsum_dw/*.mrc')
    elif status_dict['do_cter']:
        cmd.append('XXX_SP_CTER_MICROGRAPH_PATTERN_XXX')
    elif status_dict['do_cryolo']:
        cmd.append('XXX_SP_CRYOLO_MICROGRAPH_PATTERN_XXX')
    else:
        cmd.append('XXX_SP_WINDOW_MICROGRAPH_PATTERN_XXX')

    if status_dict['do_cryolo']:
        cmd.append('XXX_SP_CRYOLO_OUTPUT_DIR_XXX/EMAN/*.box')
    else:
        cmd.append('XXX_SP_WINDOW_BOX_PATTERN_XXX')

    if status_dict['do_cter']:
        cmd.append('XXX_SP_CTER_OUTPUT_DIR_XXX/partres.txt')
    else:
        cmd.append('XXX_SP_WINDOW_PARTRES_XXX')

    cmd.append('XXX_SP_WINDOW_OUTPUT_DIR_XXX')

    cmd.append('--box_size=XXX_SP_BOX_SIZE_XXX')
    if negative_stain:
        cmd.append('--skip_invert')
    cmd.append('XXX_SP_WINDOW_ADDITION_XXX')

    cmd.append(';')
    cmd.append('e2bdb.py')
    cmd.append('XXX_SP_WINDOW_OUTPUT_DIR_XXX/mpi_proc_*')
    cmd.append('--makevstack=bdb:XXX_SP_WINDOW_OUTPUT_DIR_XXX/stack')
    return cmd


def get_isac2(status_dict, phase_plate, negative_stain, **kwargs):
    cmd = []
    cmd.append('sp_isac2.py')
    if status_dict['do_window']:
        cmd.append('bdb:XXX_SP_WINDOW_OUTPUT_DIR_XXX/stack')
    else:
        cmd.append('XXX_SP_ISAC_STACK_XXX')
    cmd.append('XXX_SP_ISAC_OUTPUT_DIR_XXX')
    cmd.append('--radius=XXX_SP_PARTICLE_RADIUS_XXX')
    cmd.append('--img_per_grp=XXX_SP_ISAC_IMG_PER_GRP_XXX')
    if phase_plate:
        cmd.append('--VPP')
    elif not negative_stain:
        cmd.append('--CTF')
    cmd.append('XXX_SP_ISAC_ADDITION_XXX')
    cmd.append(';')
    cmd.append('sp_pipe.py')
    cmd.append('isac_substack')
    if status_dict['do_window']:
        cmd.append('bdb:XXX_SP_WINDOW_OUTPUT_DIR_XXX/stack')
    else:
        cmd.append('XXX_SP_ISAC_STACK_XXX')
    cmd.append('XXX_SP_ISAC_OUTPUT_DIR_XXX')
    cmd.append('XXX_SP_SUBSTACK_OUTPUT_DIR_XXX')
    return cmd


def get_rviper(status_dict, **kwargs):
    cmd = []
    cmd.append('sp_rviper.py')
    if status_dict['do_isac2']:
        cmd.append('XXX_SP_ISAC_OUTPUT_DIR_XXX/ordered_class_averages.hdf')
    else:
        cmd.append('XXX_SP_RVIPER_INPUT_STACK_XXX')
    cmd.append('XXX_SP_RVIPER_OUTPUT_DIR_XXX')
    cmd.append('--sym=XXX_SP_SYMMETRY_XXX')
    cmd.append('XXX_SP_RVIPER_ADDITION_XXX')
    return cmd


def get_adjustment(status_dict, **kwargs):
    cmd = []
    if status_dict['do_rviper']:
        cmd.append('sp_pipe.py')
        cmd.append('moon_elimination')
        cmd.append('XXX_SP_RVIPER_OUTPUT_DIR_XXX/main001/run000/refvol2.hdf')
        cmd.append('XXX_SP_ADJUSTMENT_OUTPUT_DIR_XXX')
        cmd.append('--pixel_size=XXX_SP_PIXEL_SIZE_XXX')
        if status_dict['do_isac2']:
            cmd.append('--resample_ratio=XXX_SP_ISAC_OUTPUT_DIR_XXX')
        else:
            cmd.append('--resample_ratio=XXX_SP_ADJUSTMENT_RESAMPLE_RATIO_XXX')
        cmd.append('--mol_mass=XXX_SP_MOL_MASS_XXX')
        cmd.append('XXX_SP_ADJUSTMENT_ADDITION_XXX')
    return cmd

def get_mask_rviper(status_dict, fill_rviper_mask, **kwargs):
    cmd = []
    if status_dict['do_rviper']:
        cmd.append('sp_mask.py')
        cmd.append('XXX_SP_ADJUSTMENT_OUTPUT_DIR_XXX/vol3d_ref_moon_eliminated.hdf')
        cmd.append('XXX_SP_MASK_RVIPER_OUTPUT_DIR_XXX')
        cmd.append('--mol_mass=XXX_SP_MOL_MASS_XXX')
        cmd.append('--ndilation=XXX_SP_MASK_RVIPER_NDILAITON_XXX')
        cmd.append('--edge_width=XXX_SP_MASK_RVIPER_SOFT_EDGE_XXX')
        if fill_rviper_mask:
            cmd.append('--fill_mask')
        cmd.append('XXX_SP_MASK_RVIPER_ADDITION_XXX')
    return cmd


def get_meridien(status_dict, **kwargs):
    cmd = []
    cmd.append('sp_meridien.py')
    if status_dict['do_isac2']:
        cmd.append('bdb:XXX_SP_SUBSTACK_OUTPUT_DIR_XXX/isac_substack')
    else:
        cmd.append('XXX_SP_MERIDIEN_INPUT_STACK_XXX')
    cmd.append('XXX_SP_MERIDIEN_OUTPUT_DIR')
    if status_dict['do_rviper']:
        cmd.append('XXX_SP_ADJUSTMENT_OUTPUT_DIR_XXX/vol3d_ref_moon_eliminated.hdf')
    else:
        cmd.append('XXX_SP_MERIDIEN_INPUT_VOLUME_XXX')
    cmd.append('--radius=XXX_SP_PARTICLE_RADIUS_XXX')
    cmd.append('--symmetry=XXX_SP_SYMMETRY_XXX')
    cmd.append('--mask3D=XXX_SP_MASK_RVIPER_OUTPUT_DIR_XXX/sxmask_mask.hdf')
    if status_dict['do_isac2']:
        cmd.append('--initialshifts')
        cmd.append('--skip_prealignment')
    cmd.append('XXX_SP_MERIDIEN_ADDITION_XXX')
    return cmd


def show_variables():
    with open(__file__) as r:
        lines = r.readlines()
    xxx_comp = re.compile('XXX_SP_[^ \'/]*_XXX')
    matches = []
    for line in lines:
        for match in xxx_comp.findall(line):
            matches.append(match)
    import numpy as np
    for entry in np.sort(np.unique(matches)):
        print('UNIQUE VAR:', entry)


def get_defaults():
    default_dict = {}
    default_dict['XXX_SP_BOX_SIZE_XXX'] = '352'
    default_dict['XXX_SP_MOL_MASS_XXX'] = '1400'
    default_dict['XXX_SP_PARTICLE_RADIUS_XXX'] = '145'
    default_dict['XXX_SP_PIXEL_SIZE_XXX'] = '1.14'
    default_dict['XXX_SP_SYMMETRY_XXX'] = 'c5'
    default_dict['XXX_SP_VOLTAGE_XXX'] = '300'

    default_dict['XXX_SP_ADJUSTMENT_ADDITION_XXX'] = ''
    default_dict['XXX_SP_ADJUSTMENT_OUTPUT_DIR_XXX'] = '04b_RVIPER_ADJUSTMENT'
    default_dict['XXX_SP_ADJUSTMENT_RESAMPLE_RATIO_XXX'] = '1.0'

    default_dict['XXX_SP_CRYOLO_ADDITION_XXX'] = ''
    default_dict['XXX_SP_CRYOLO_CONFIG_PATH_XXX'] = '/Path/to/cryolo'
    default_dict['XXX_SP_CRYOLO_GPU_XXX'] = '-1'
    default_dict['XXX_SP_CRYOLO_MICROGRAPH_PATTERN_XXX'] = 'bla/*.mrc'
    default_dict['XXX_SP_CRYOLO_OUTPUT_DIR_XXX'] = '02a_CRYOLO_PREDICT'
    default_dict['XXX_SP_CRYOLO_PREDICT_PATH_XXX'] = '/Path/to/cryolo_predict'

    default_dict['XXX_SP_CTER_ADDITION_XXX'] = ''
    default_dict['XXX_SP_CTER_CS_XXX'] = '2.7'
    default_dict['XXX_SP_CTER_MICROGRAPH_PATTERN_XXX'] = ''
    default_dict['XXX_SP_CTER_OUTPUT_DIR_XXX'] = '01a_CTER'
    default_dict['XXX_SP_CTER_WINDOW_SIZE'] = '1024'

    default_dict['XXX_SP_ISAC_ADDITION_XXX'] = ''
    default_dict['XXX_SP_ISAC_IMG_PER_GRP_XXX'] = '100'
    default_dict['XXX_SP_ISAC_OUTPUT_DIR_XXX'] = '03a_ISAC'
    default_dict['XXX_SP_ISAC_STACK_XXX'] = 'bdb:stack'

    default_dict['XXX_SP_MASK_RVIPER_ADDITION_XXX'] = ''
    default_dict['XXX_SP_MASK_RVIPER_NDILAITON_XXX'] = '3'
    default_dict['XXX_SP_MASK_RVIPER_OUTPUT_DIR_XXX'] = '04c_RVIPER_MASK'
    default_dict['XXX_SP_MASK_RVIPER_SOFT_EDGE_XXX'] = '10'

    default_dict['XXX_SP_MERIDIEN_ADDITION_XXX'] = ''
    default_dict['XXX_SP_MERIDIEN_INPUT_STACK_XXX'] = 'bdb:stack'
    default_dict['XXX_SP_MERIDIEN_INPUT_VOLUME_XXX'] = 'input_volume'
    default_dict['XXX_SP_MERIDIEN_OUTPUT_DIR'] = '05a_MERIDIEN'

    default_dict['XXX_SP_RVIPER_ADDITION_XXX'] = ''
    default_dict['XXX_SP_RVIPER_INPUT_STACK_XXX'] = 'bdb:classes'
    default_dict['XXX_SP_RVIPER_OUTPUT_DIR_XXX'] = '04a_RVIPER'

    default_dict['XXX_SP_SUBSTACK_OUTPUT_DIR_XXX'] = '03b_SUBSTACK'

    default_dict['XXX_SP_UNBLUR_ADDITION_XXX'] = ''
    default_dict['XXX_SP_UNBLUR_EXP_PER_FRAME_XXX'] = '2.5'
    default_dict['XXX_SP_UNBLUR_FIRST_FRAME_XXX'] = '0'
    default_dict['XXX_SP_UNBLUR_GAIN_FILE_XXX'] = '/Path/to/Gain'
    default_dict['XXX_SP_UNBLUR_LAST_FRAME_XXX'] = '-1'
    default_dict['XXX_SP_UNBLUR_MICROGRAPH_PATTERN_XXX'] = '/PATTERN*.mrc'
    default_dict['XXX_SP_UNBLUR_OUTPUT_DIR_XXX'] = '00a_UNBLUR'
    default_dict['XXX_SP_UNBLUR_PATH_XXX'] = '/Path/to/unblur'

    default_dict['XXX_SP_WINDOW_ADDITION_XXX'] = ''
    default_dict['XXX_SP_WINDOW_BOX_PATTERN_XXX'] = 'box*.box'
    default_dict['XXX_SP_WINDOW_MICROGRAPH_PATTERN_XXX'] = 'mrc*.mrc'
    default_dict['XXX_SP_WINDOW_OUTPUT_DIR_XXX'] = '02b_WINDOW'
    default_dict['XXX_SP_WINDOW_PARTRES_XXX'] = 'Wuhu/partres'
    return default_dict


def main(args_as_dict):
    #show_variables()
    function_dict = co.OrderedDict()
    function_dict['do_unblur'] = get_unblur_cmd
    function_dict['do_cter'] = get_cter_cmd
    function_dict['do_cryolo'] = get_cryolo_predict
    function_dict['do_window'] = get_window
    function_dict['do_isac2'] = get_isac2
    function_dict['do_rviper'] = get_rviper
    function_dict['do_adjustment'] = get_adjustment
    function_dict['do_mask_rviper'] = get_mask_rviper
    function_dict['do_meridien'] = get_meridien

    do_dict = {}
    do_dict['do_unblur'] = True
    do_dict['do_cter'] = True
    do_dict['do_cryolo'] = True
    do_dict['do_window'] = True
    do_dict['do_isac2'] = True
    do_dict['do_rviper'] = True
    do_dict['do_adjustment'] = True
    do_dict['do_mask_rviper'] = True
    do_dict['do_meridien'] = True

    phase_plate = False
    negative_stain = False
    fill_rviper_mask = False

    mpi_procs = 96
    mpi_submission = 'test_submission.sh'
    out_submission = 'out_submission.sh'

    default_replaces = get_defaults()

    cmds = []
    for key in function_dict:
        if do_dict[key]:
            cmds.append(function_dict[key](phase_plate=phase_plate, negative_stain=negative_stain, fill_rviper_mask=fill_rviper_mask, status_dict=do_dict))

    with open(mpi_submission) as read:
        lines = read.readlines()
    found_lines = []
    for idx, entry in enumerate(lines):
        if 'XXX_SXMPI_NPROC_XXX' in entry and 'XXX_SXCMD_LINE_XXX' in entry and 'mpirun' in entry:
            found_lines.append(idx)
    line = lines[found_lines[-1]].replace('{', '{{').replace('}', '}}').replace('XXX_SXMPI_NPROC_XXX', str(mpi_procs)).replace('XXX_SXCMD_LINE_XXX', '{0}')
    final_line = '\n'.join([line.format(' '.join(entry)) for entry in cmds])

    for key, value in default_replaces.items():
        final_line = final_line.replace(key, value)
    lines[found_lines[-1]] = final_line

    with open(out_submission, 'w') as w:
        w.write('\n'.join(lines))


if __name__ == '__main__':
    main(vars(parse_args()))
