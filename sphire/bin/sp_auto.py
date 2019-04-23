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


import os
import sys
import re
import argparse
import subprocess
import collections as co

import sp_global_def


def parse_args():
	parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=9999))
	parser.add_argument('output_directory', help='Output directory to store the outputs')

	group = parser.add_argument_group('MPI settings (required)')
	group.add_argument('--mpi_procs', type=int, default=2, help='Number of processors to use.')
	group.add_argument('--mpi_job_name', type=str, default='auto_sphire', help='Job name of the submitted job.')
	group.add_argument('--mpi_submission_command', type=str, default='sbatch', help='Submission command, e.g. sbatch, qsub, ...')
	group.add_argument('--mpi_submission_template', type=str, default='TEMPLATES/submit.sh', help='Submission template.')

	group = parser.add_argument_group('Global settings (required)')
	group.add_argument('--apix', dest='XXX_SP_PIXEL_SIZE_XXX', type=float, default=1.0, help='Pixel size in A/pixel.')
	group.add_argument('--mol_mass', dest='XXX_SP_MOL_MASS_XXX', type=float, default=250.0, help='Molecular mass of the protein in kDa. Used to calculate the masking density threshold.')
	group.add_argument('--radius', dest='XXX_SP_PARTICLE_RADIUS_XXX', type=int, default=80, help='Particle radius in pixels. Used for normalization.')
	group.add_argument('--box_size', dest='XXX_SP_BOX_SIZE_XXX', type=int, default=200, help='Particle box size in pixels.')
	group.add_argument('--symmetry', dest='XXX_SP_SYMMETRY_XXX', type=str, default='c1', help='Symmetry of the particle.')
	group.add_argument('--voltage', dest='XXX_SP_VOLTAGE_XXX', type=float, default=300.0, help='Microscope voltage in kV.')
	group.add_argument('--negative_stain', action='store_true', default=False, help='Input is negative stain.')
	group.add_argument('--phase_plate', action='store_true', default=False, help='Input is phase_plate.')
	group.add_argument('--fill_rviper_mask', action='store_true', default=False, help='Fill RVIPER mask.')

	group = parser.add_argument_group('Unblur settings (required to run movie alignment)')
	group.add_argument('--unblur_path', dest='XXX_SP_UNBLUR_PATH_XXX', type=str, default=None, help='Path pointing to the unblur executable.')
	group.add_argument('--unblur_mic_pattern', dest='XXX_SP_UNBLUR_MICROGRAPH_PATTERN_XXX', type=str, default=None, help='Pattern of the micrographs to use for motion correction.')
	group.add_argument('--unblur_exp_per_frame', dest='XXX_SP_UNBLUR_EXP_PER_FRAME_XXX', type=float, default=None, help='Exposure per frame. Used for dose adjustment.')
	group.add_argument('--unblur_gain_file', dest='XXX_SP_UNBLUR_GAIN_FILE_XXX', type=str, default=None, help='File containing the information for gain correction. Not required if the movies are already gain corrected.')

	group = parser.add_argument_group('Unblur settings (optional)')
	group.add_argument('--skip_unblur', action='store_true', default=False, help='Do not run motion correction')
	group.add_argument('--unblur_output_dir', dest='XXX_SP_UNBLUR_OUTPUT_DIR_XXX', type=str, default='00a_UNBLUR', help='Unblur output directory.')
	group.add_argument('--unblur_addition', dest='XXX_SP_UNBLUR_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('CTER settings (required to run CTF estimation)')
	group.add_argument('--cter_cs', dest='XXX_SP_CTER_CS_XXX', type=float, default=2.7, help='Spherical aberration of the microscope.')

	group = parser.add_argument_group('CTER settings (optional)')
	group.add_argument('--skip_cter', action='store_true', default=False, help='Do not run CTF estimation.')
	group.add_argument('--cter_output_dir', dest='XXX_SP_CTER_OUTPUT_DIR_XXX', type=str, default='01a_CTER', help='CTER output directory.')
	group.add_argument('--cter_mic_pattern', dest='XXX_SP_CTER_MICROGRAPH_PATTERN_XXX', type=str, default='Mics/*.mrc', help='Micrograph pattern in case unblur is skipped.')
	group.add_argument('--cter_window_size', dest='XXX_SP_CTER_WINDOW_SIZE', type=int, default=1024, help='CTF estimation window size.')
	group.add_argument('--cter_addition', dest='XXX_SP_CTER_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('CRYOLO settings (required to run particle picking)')
	group.add_argument('--cryolo_predict_path', dest='XXX_SP_CRYOLO_PREDICT_PATH_XXX', type=str, default='/Path/cryolo_predict.py', help='Path to the cryolo predict executable.')
	group.add_argument('--cryolo_config_path', dest='XXX_SP_CRYOLO_CONFIG_PATH_XXX', type=str, default=None, help='Path to the cryolo config file')
	group.add_argument('--cryolo_model_path', dest='XXX_SP_CRYOLO_MODEL_PATH_XXX', type=str, default=None, help='Path to the cryolo model file')
	group.add_argument('--cryolo_gpu', dest='XXX_SP_CRYOLO_GPU_XXX', type=str, default='0', help='Cryolo GPU list.')

	group = parser.add_argument_group('CRYOLO settings (optional)')
	group.add_argument('--skip_cryolo', action='store_true', default=False, help='Do not run particle picking.')
	group.add_argument('--cryolo_output_dir', dest='XXX_SP_CRYOLO_OUTPUT_DIR_XXX', type=str, default='02a_CRYOLO_PREDICT', help='CRYOLO output directory.')
	group.add_argument('--cryolo_mic_path', dest='XXX_SP_CRYOLO_MICROGRAPH_PATH_XXX', type=str, default='Mics/*.mrc', help='Micrograph pattern in case unblur is skipped.')
	group.add_argument('--cryolo_addition', dest='XXX_SP_CRYOLO_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('WINDOW settings (optional)')
	group.add_argument('--skip_window', action='store_true', default=False, help='Do not run particle extraction.')
	group.add_argument('--window_box_pattern', dest='XXX_SP_WINDOW_BOX_PATTERN_XXX', type=str, default='Boxes/*.box', help='Window box file pattern.')
	group.add_argument('--window_mic_pattern', dest='XXX_SP_WINDOW_MICROGRAPH_PATTERN_XXX', type=str, default='Mics/*.mrc', help='Window mrc file pattern.')
	group.add_argument('--window_partres', dest='XXX_SP_WINDOW_PARTRES_XXX', type=str, default='CTER/partres.txt', help='CTER partres file. In case of negative stain put this value to the pixel size.')
	group.add_argument('--window_output_dir', dest='XXX_SP_WINDOW_OUTPUT_DIR_XXX', type=str, default='02b_WINDOW', help='Window output directory.')
	group.add_argument('--window_addition', dest='XXX_SP_WINDOW_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('ISAC2 settings (required to run 2d classification)')
	group.add_argument('--isac2_img_per_grp', dest='XXX_SP_ISAC_IMG_PER_GRP_XXX', type=int, default=100, help='Img per group for the ISAC run.')

	group = parser.add_argument_group('ISAC2 settings (optional)')
	group.add_argument('--skip_isac2', action='store_true', default=False, help='Do not run 2d classification.')
	group.add_argument('--isac2_input_stack', dest='XXX_SP_ISAC_STACK_XXX', type=str, default='bdb:path#stack', help='Path to the Input stack for ISAC')
	group.add_argument('--isac2_output_dir', dest='XXX_SP_ISAC_OUTPUT_DIR_XXX', type=str, default='03a_ISAC', help='ISAC2 output directory.')
	group.add_argument('--isac2_addition', dest='XXX_SP_ISAC_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('Substack ISAC2 settings (optional)')
	group.add_argument('--substack_output_dir', dest='XXX_SP_SUBSTACK_OUTPUT_DIR_XXX', type=str, default='03b_SUBSTACK', help='Substack ISAC2 output directory.')

	group = parser.add_argument_group('RVIPER settings (optional)')
	group.add_argument('--skip_rviper', action='store_true', default=False, help='Do not run 3d ab-initio reconstruction.')
	group.add_argument('--rviper_input_stack', dest='XXX_SP_RVIPER_INPUT_STACK_XXX', type=str, default='bdb:path#stack', help='Path to the input stack for RVIPER')
	group.add_argument('--rviper_output_dir', dest='XXX_SP_RVIPER_OUTPUT_DIR_XXX', type=str, default='04a_RVIPER', help='RVIPER output directory.')
	group.add_argument('--rviper_addition', dest='XXX_SP_RVIPER_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('RVIPER volume adjustment settings (optional)')
	group.add_argument('--skip_adjust_rviper', action='store_true', default=False, help='Skip adjusting a volume.')
	group.add_argument('--adjust_rviper_resample', dest='XXX_SP_ADJUSTMENT_RESAMPLE_RATIO_XXX', type=str, default='bdb:path#stack', help='Resample ratio for RVIPER.')
	group.add_argument('--adjust_rviper_output_dir', dest='XXX_SP_ADJUSTMENT_OUTPUT_DIR_XXX', type=str, default='04b_RVIPER_ADJUSTMENT', help='RVIPER volume adjustment output directory.')
	group.add_argument('--adjust_rviper_addition', dest='XXX_SP_ADJUSTMENT_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('RVIPER mask settings (optional)')
	group.add_argument('--skip_mask_rviper', action='store_true', default=False, help='Skip creating a mask.')
	group.add_argument('--mask_rviper_ndilation', dest='XXX_SP_MASK_RVIPER_NDILAITON_XXX', type=int, default=3, help='Number of dilations of the mask. 1 Dilation adds about 2 pixel to the binary volume.')
	group.add_argument('--mask_rviper_soft_edge', dest='XXX_SP_MASK_RVIPER_SOFT_EDGE_XXX', type=int, default=10, help='Number of pixels for the soft edge.')
	group.add_argument('--mask_rviper_output_dir', dest='XXX_SP_MASK_RVIPER_OUTPUT_DIR_XXX', type=str, default='04c_RVIPER_MASK', help='RVIPER mask output directory.')
	group.add_argument('--mask_rviper_addition', dest='XXX_SP_MASK_RVIPER_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	group = parser.add_argument_group('Meridien settings (optional)')
	group.add_argument('--skip_meridien', action='store_true', default=False, help='Do not run 3d refinement.')
	group.add_argument('--meridien_input_volume', dest='XXX_SP_MERIDIEN_INPUT_VOLUME_XXX', type=str, default='ref_vol.hdf', help='Path to the ref_vol.hdf file')
	group.add_argument('--meridien_input_mask', dest='XXX_SP_MERIDIEN_INPUT_MASK_XXX', type=str, default='mask.hdf', help='Path to the mask.hdf file')
	group.add_argument('--meridien_input_stack', dest='XXX_SP_MERIDIEN_INPUT_STACK_XXX', type=str, default='bdb:path#stack', help='Path to the Input stack for Meridien')
	group.add_argument('--meridien_output_dir', dest='XXX_SP_MERIDIEN_OUTPUT_DIR', type=str, default='05a_MERIDIEN', help='Meridien output directory.')
	group.add_argument('--meridien_addition', dest='XXX_SP_MERIDIEN_ADDITION_XXX', type=str, default='', help='Additional parameters that are not part of the required ones.')

	args = parser.parse_args()
	return args


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
		cmd.append('XXX_SP_CRYOLO_MICROGRAPH_PATH_XXX')
	cmd.append('XXX_SP_CRYOLO_MODEL_PATH_XXX')
	cmd.append('XXX_SP_CRYOLO_OUTPUT_DIR_XXX')
	cmd.append('--gpu=XXX_SP_CRYOLO_GPU_XXX')
	cmd.append('--cryolo_predict_path=XXX_SP_CRYOLO_PREDICT_PATH_XXX')
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
		cmd.append('XXX_SP_CRYOLO_MICROGRAPH_PATH_XXX/*.mrc')
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
	return cmd

def get_window_stack(**kwargs):
	cmd = []
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
	return cmd

def get_isac2_substack(status_dict, **kwargs):
	cmd = []
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
	if status_dict['do_rviper'] and status_dict['do_adjust_rviper']:
		cmd.append('sp_pipe.py')
		cmd.append('moon_eliminator')
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
	if status_dict['do_rviper'] and status_dict['do_mask_rviper'] and status_dict['do_adjust_rviper']:
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
	if status_dict['do_mask_rviper']:
		cmd.append('--mask3D=XXX_SP_MASK_RVIPER_OUTPUT_DIR_XXX/sxmask_mask.hdf')
	else:
		cmd.append('--mask3D=XXX_SP_MERIDIEN_INPUT_MASK_XXX')
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
	default_dict['XXX_SP_UNBLUR_GAIN_FILE_XXX'] = '/Path/to/Gain'
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
	if os.path.exists(args_as_dict['output_directory']):
		sp_global_def.sxprint('Output directory already exists! Please choose another one.')
		sys.exit(1)

	#show_variables()
	function_dict = co.OrderedDict()
	# [Function name, MPI support]
	function_dict['do_unblur'] = [get_unblur_cmd, True]
	function_dict['do_cter'] = [get_cter_cmd, True]
	function_dict['do_cryolo'] = [get_cryolo_predict, False]
	function_dict['do_window'] = [get_window, True]
	function_dict['do_window_stack'] = [get_window_stack, False]
	function_dict['do_isac2'] = [get_isac2, True]
	function_dict['do_isac2_substack'] = [get_isac2_substack, False]
	function_dict['do_rviper'] = [get_rviper, True]
	function_dict['do_adjust_rviper'] = [get_adjustment, False]
	function_dict['do_mask_rviper'] = [get_mask_rviper, False]
	function_dict['do_meridien'] = [get_meridien, True]

	do_dict = {}
	for key, value in args_as_dict.items():
		if key.startswith('skip') and isinstance(value, bool):
			do_dict[key.replace('skip', 'do')] = bool(not value)
			if key == 'skip_window':
				do_dict['{0}_stack'.format(key.replace('skip', 'do'))] = bool(not value)
			elif key == 'skip_isac2':
				do_dict['{0}_substack'.format(key.replace('skip', 'do'))] = bool(not value)

	phase_plate = args_as_dict['phase_plate']
	negative_stain = args_as_dict['negative_stain']
	fill_rviper_mask = args_as_dict['fill_rviper_mask']

	mpi_procs = args_as_dict['mpi_procs']
	mpi_submission = args_as_dict['mpi_submission_template']
	out_submission = '{0}/submission_script.sh'.format(args_as_dict['output_directory'])

	prev_line = ' echo $(date) {0}\n'
	check_line = "if [[ ${{?}} != 0 ]]; then echo $(date) '{0}: Failure!'; exit 1; else echo $(date) '{0}: Success!'; fi\n"

	cmds = []
	for key in function_dict:
		if do_dict[key]:
			cmds.append(prev_line.format(key))
			return_value = [entry for entry in function_dict[key][0](phase_plate=phase_plate, negative_stain=negative_stain, fill_rviper_mask=fill_rviper_mask, status_dict=do_dict) if entry.strip()]
			return_value.insert(0, [function_dict[key][1], key])
			cmds.append(return_value)
			cmds.append(check_line.format(key))

	with open(mpi_submission) as read:
		lines = read.readlines()
	found_lines = []
	for idx, entry in enumerate(lines):
		if 'XXX_SXCMD_LINE_XXX' in entry and 'mpirun' in entry:
			found_lines.append(idx)

	if not found_lines:
		sp_global_def.sxprint('Could not find a suitable command line for exchange.')
		sp_global_def.sxprint('The line should contain XXX_SXMPI_NPROC_XXX and XXX_SXCMD_LINE_XXX.')
		sys.exit(1)

	line = lines[found_lines[-1]].replace('{', '{{').replace('}', '}}').replace('XXX_SXMPI_NPROC_XXX', str(mpi_procs)).replace('XXX_SXCMD_LINE_XXX', '\'{0}\' >> {1}_out.txt 2>>{1}_err.txt')
	line_no_mpi = '\'{0}\' >>{1}_out.txt 2>>{1}_err.txt\n'

	final_line = '\n'.join([
		line.format('\' \''.join(entry[1:]), os.path.join(args_as_dict['output_directory'], entry[0][1]))
		if isinstance(entry, list) and entry[0][0]
		else
			line_no_mpi.format('\' \''.join(entry[1:]), os.path.join(args_as_dict['output_directory'], entry[0][1]))
		if isinstance(entry, list)
		else
			entry
		for entry in cmds
		])

	for key, value in args_as_dict.items():
		if key.startswith('XXX'):
			final_line = final_line.replace('\'XXX_SP_WINDOW_OUTPUT_DIR_XXX/mpi_proc_*\'', 'XXX_SP_WINDOW_OUTPUT_DIR_XXX/mpi_proc_*')
			if 'OUTPUT_DIR' in key:
				value = os.path.join(args_as_dict['output_directory'], value)
			final_line = final_line.replace(key, str(value))
	lines[found_lines[-1]] = final_line.replace(" ''", '').replace('\';\'', ';')

	os.mkdir(args_as_dict['output_directory'])
	sp_global_def.write_command(args_as_dict['output_directory'])
	with open(out_submission, 'w') as w:
		w.write(''.join(lines).replace('XXX_SXMPI_NPROC_XXX', str(mpi_procs)).replace('XXX_SXMPI_JOB_NAME_XXX', args_as_dict['mpi_job_name']))
	sp_global_def.sxprint(subprocess.check_output([args_as_dict['mpi_submission_command'], out_submission]))


if __name__ == '__main__':
	sp_global_def.print_timestamp( "Start" )
	main(vars(parse_args()))
	sp_global_def.print_timestamp( "Finish" )