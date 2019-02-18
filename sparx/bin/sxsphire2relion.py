#!/usr/bin/env python

# Author: Markus Stabrin 2018/09/28 (markus.stabrin@mpi-dortmund.mpg.de)
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete SPHIRE and EMAN2 software packages have some GPL dependencies,
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
# ========================================================================================
from __future__ import print_function, division

import os
import argparse

import numpy as np

import EMAN2_cppwrap
import global_def


def parse_args():
	"""
	Parse command line arguments.

	The command line arguments:
	--particle_stack
	--partres_file
	can be used next to each other, while the command line arguments:
	--params_2d_file
	--params_3d_file
	--params_3d_index_file
	--list
	--exlist
	do require the --particle_stack to be set.

	--list and --exlist cannot be used at the same time.

	Arguments:
	None

	Returns:
	Parsed arguments
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('output_dir', type=str, help='Output directory')
	parser.add_argument('--output_name', type=str, default='sphire2relion.star', help='Output star file name')
	parser.add_argument('--force', action='store_true', default=False, help='Overwrite existing star file.')
	parser.add_argument('--partres_file', type=str, help='Partres file')
	parser.add_argument('--particle_stack', type=str, help='Particle stack in bdb or hdf format')
	parser.add_argument('--params_2d_file', type=str, help='2D alignment parameters. Requires --particle_stack. Cannot be used together with --params_3d_file.')
	parser.add_argument('--params_3d_file', type=str, help='3D projection parameters. Requires --particle_stack. Cannot be used together with --params_2d_file.')
	parser.add_argument('--params_3d_index_file', type=str, help='Index file for the 3d params. Used to find the associated particle stack entry in the params file. In the meridien directories, this file is either called chunk or index. Requires --particle_stack. Requires --params_3d_file.')
	parser.add_argument('--params_3d_chunk_files', type=str, nargs=2, help='Chunk files for the 3d params. Used to extract the _rlnRandomSubset information. In the meridien directories, this file is called chunk. Requires --particle_stack. Requires --params_3d_file.')
	parser.add_argument('--list', type=str, help='List of particles to include. Requires --particle_stack. Cannot be used together with --exlist.')
	parser.add_argument('--exlist', type=str, help='List of particles to exclude. Requires --particle_stack. Cannot be used together with --list.')
	return  parser.parse_args()


def main(args):
	"""
	Main function contains the logic behind the program.

	Arguments:
	args - Command line argument options

	Returns:
	None
	"""
	global_def.BATCH = True
	output_file = sanity_checks(args)

	output_dtype = []
	array_list = []
	create_stack = False
	particle_data = None
	partres_data = None
	params_2d_data = None
	params_3d_data = None
	params_3d_subset_data = None

	if args.particle_stack:
		print('Import particle stack')
		particle_data, create_stack = import_particle_stack(args.particle_stack, args.output_dir)
		output_dtype.extend(particle_data.dtype.descr)

	if args.partres_file:
		print('Import partres file')
		partres_data = import_partres_file(args.partres_file)
		output_dtype.extend(partres_data.dtype.descr)

	if args.params_2d_file:
		print('Import params 2d file')
		params_2d_data = import_params(args.params_2d_file, dim='2d')
		output_dtype.extend(params_2d_data.dtype.descr)

	if args.params_3d_file:
		print('Import params 3d file')
		params_3d_data = import_params(args.params_3d_file, dim='3d')
		output_dtype.extend(params_3d_data.dtype.descr)

	if args.params_3d_chunk_files:
		print('Import params 3d chunk files')
		params_3d_subset_data = np.empty(params_3d_data.shape[0], dtype=[('_rlnRandomSubset', '<i8')])
		params_3d_subset_data.fill(np.nan)
		params_import = []
		for idx, file_name in enumerate(args.params_3d_chunk_files):
			chunk_import = np.genfromtxt(file_name, int)
			params_3d_subset_data['_rlnRandomSubset'][chunk_import] = idx
			params_import.extend(chunk_import.tolist())
		output_dtype.extend(params_3d_subset_data.dtype.descr)
		params_3d_subset_data = params_3d_subset_data[params_index_data]
		assert params_3d_subset_data.shape[0] == params_3d_data.shape[0]
		assert np.unique(params_import).shape[0] == params_3d_data.shape[0]

	if args.params_3d_index_file:
		print('Import params 3d index')
		params_index_data = np.genfromtxt(args.params_3d_index_file, dtype=int)
		assert params_3d_data.shape[0] == params_index_data.shape[0]
		assert np.unique(params_index_data).shape[0] == params_3d_data.shape[0]
	elif args.particle_stack:
		params_index_data = np.arange(particle_data.shape[0])
	else:
		params_index_data = np.arange(partres_data.shape[0])

	mask_array = np.ones(particle_data.shape[0], dtype=np.bool)
	if args.list or args.exlist:
		print('Import list/exlist information')
		mask_array = create_particle_data_mask(args.list, args.exlist, particle_data.shape[0])

	output_data = np.empty(params_index_data.shape[0], dtype=sorted(list(set(output_dtype))))
	particle_data_params = particle_data[params_index_data]
	mask_array_params = mask_array[params_index_data]

	array_list.append(particle_data_params)
	array_list.append(params_2d_data)
	array_list.append(params_3d_data)
	array_list.append(params_3d_subset_data)

	print('Adjust header')
	for array in array_list:
		if array is not None:
			for name in array.dtype.names:
				output_data[name] = array[name]

	if partres_data is not None:
		for row in partres_data:
			mask = output_data['_rlnMicrographName'] == row['_rlnMicrographName']
			for name in partres_data.dtype.names:
				output_data[name][mask] = row[name]

	final_output = output_data[mask_array_params]

	print('Write star file')
	header = ['', 'data_', '', 'loop_']
	header.extend(['{0} #{1}'.format(name, idx+1) for idx, name in enumerate(final_output.dtype.names)])
	dtype_dict = final_output.dtype.fields
	fmt = []
	for name in final_output.dtype.names:
		max_length = len(max([str(entry).split('.')[0] for entry in final_output[name]], key=len)) + 2
		if 'float' in str(dtype_dict[name][0]):
			fmt.append('%{0}.6f'.format(max_length+7))
		elif 'int' in str(dtype_dict[name][0]):
			fmt.append('%{0}d'.format(max_length))
		elif '|S' in str(dtype_dict[name][0]):
			fmt.append('%{0}s'.format(max_length))
		else:
			assert False
	np.savetxt(output_file, final_output, fmt=' '.join(fmt), header='\n'.join(header), comments='')

	if create_stack:
		print('Create particle stacks')
		create_particle_stack(args.particle_stack, args.output_dir, particle_data)

	print('Done!')
	global_def.BATCH = False


def create_particle_stack(particle_stack, output_dir, particle_data):
	"""
	Create a mrcs particle stack based on the bdb/hdf input stack.

	Arguments:
	particle_stack - Particle stack name
	output_dir - Output directory
	particle_data - Particle_data array

	Returns:
	None
	"""
	print('|_Get particle ID and particle names')
	ptcl_ids = [int(entry.split('@')[0]) for entry in particle_data['_rlnImageName']]
	ptcl_names = [entry.split('@')[1] for entry in particle_data['_rlnImageName']]

	print('|_Write images')
	for particle_idx in range(particle_data.shape[0]):
		if particle_idx % 10000 == 0:
			print(particle_idx, ' of ', particle_data.shape[0])
		emdata = EMAN2_cppwrap.EMData(particle_stack, particle_idx)

		output_name = os.path.join(output_dir, ptcl_names[particle_idx])
		try:
			os.makedirs(os.path.dirname(output_name))
		except OSError:
			pass
		emdata.write_image(output_name, -1)


def create_particle_data_mask(list_file, exlist_file, mask_size):
	"""
	Create masked array based on the particle selection files.

	Arguments:
	list_file - File containing the indices of particles to keep.
	exlist_file - File containing the indices of particles to remove.
	mask_size - Maximum size of the masked array. Is equals the number of particles.

	Returns:
	Masked boolean array.
	"""
	mask_array = np.zeros(mask_size, dtype=int)
	if list_file:
		mask_array[np.genfromtxt(list_file, int)] = 1
	elif exlist_file:
		mask_array[np.genfromtxt(exlist_file, int)] = -1
		mask_array += 1
	else:
		assert False
	return mask_array.astype(np.bool)


def import_params(params_file, dim):
	"""
	Import 2D  or 3D parameters.
	The 2d params file has the format:
	angle, shift x, shift y, mirror

	The 3d params file has the format:
	angle_phi, angle_theta, angle_psi, tx, ty

	Arguments:
	params_2d_file - File name of the 2d parameter files

	Returns:
	parameter array
	"""
	if dim == '2d':
		dtype_import_list = [
			('angle_psi', float),
			('shift_x', float),
			('shift_y', float),
			('mirror', int),
			]
		print('What happens with mirror?')
	elif dim == '3d':
		dtype_import_list = [
			('angle_rot', float),
			('angle_theta', float),
			('angle_psi', float),
			('shift_x', float),
			('shift_y', float),
			]

	input_data = np.genfromtxt(params_file, dtype=dtype_import_list)

	dtype_output_list = [
		('_rlnOriginX', float),
		('_rlnOriginY', float),
		('_rlnAngleRot', float),
		('_rlnAngleTilt', float),
		('_rlnAnglePsi', float),
		]
	params_array = np.empty(len(input_data), dtype=sorted(dtype_output_list))

	params_array['_rlnOriginX'] = input_data['shift_x']
	params_array['_rlnOriginY'] = input_data['shift_y']
	params_array['_rlnAnglePsi'] = input_data['angle_psi']

	if dim == '2d':
		params_array['_rlnAngleTilt'] = 0
		params_array['_rlnAngleRot'] = 0
	elif dim == '3d':
		params_array['_rlnAngleTilt'] = input_data['angle_theta']
		params_array['_rlnAngleRot'] = input_data['angle_rot']

	return params_array


def import_partres_file(partres_file):
	"""
	Import the information from a SPHIRE partres file.

	Arguments:
	partres_file - Partres file name

	Returns:
	Array containing the ctf information.
	"""
	with open(partres_file, 'r') as r:
		number_of_columns = len(r.readline().split())

	if number_of_columns == 22:
		columns = [0, 1, 2, 3, 6, 7, 17, 19, 20, 21]
		dtype_import_list = [
			('defocus', float),
			('cs', float),
			('voltage', float),
			('pixel_size', float),
			('astig_amp', float),
			('astig_angle', float),
			('max_resolution', float),
			('amplitude_contrast', float),
			('phase_shift', float),
			('micrograph_name', '|S1000'),
			]
		dtype_output_list = [
			('_rlnDefocusU', float),
			('_rlnDefocusV', float),
			('_rlnDefocusAngle', float),
			('_rlnMicrographName', '|S1000'),
			('_rlnDetectorPixelSize', float),
			('_rlnMagnification', float),
			('_rlnCtfMaxResolution', float),
			('_rlnPhaseShift', float),
			('_rlnAmplitudeContrast', float),
			('_rlnSphericalAberration', float),
			('_rlnVoltage', float),
			]
	else:
		global_def.ERROR('Number of columns in partres file not known: {0}'.format(number_of_columns), 'sxsphire2relion', 1)

	assert len(columns) == len(dtype_import_list)
	partres_import_array = np.genfromtxt(partres_file, dtype=dtype_import_list, usecols=columns)
	partres_array = np.empty(partres_import_array.shape[0], sorted(dtype_output_list))

	partres_array['_rlnDefocusU'] = (20000 * partres_import_array['defocus'] - 10000 * partres_import_array['astig_amp']) / 2
	partres_array['_rlnDefocusV'] = 20000 * partres_import_array['defocus'] - partres_array['_rlnDefocusU']
	partres_array['_rlnDefocusAngle'] = 45 - partres_import_array['astig_angle']
	partres_array['_rlnMicrographName'] = partres_import_array['micrograph_name']
	partres_array['_rlnAmplitudeContrast'] = partres_import_array['amplitude_contrast'] / 100
	partres_array['_rlnVoltage'] = partres_import_array['voltage']
	partres_array['_rlnSphericalAberration'] = partres_import_array['cs']
	partres_array['_rlnPhaseShift'] = partres_import_array['phase_shift']
	partres_array['_rlnDetectorPixelSize'] = partres_import_array['pixel_size']
	partres_array['_rlnMagnification'] = 10000
	partres_array['_rlnCtfMaxResolution'] = 1 / partres_import_array['max_resolution']

	return partres_array


def sanity_checks(args):
	"""
	Check if the settings are valid.

	Arguments:
	args - Command line arguments parsed via argparse

	Returns:
	Output file name
	"""
	file_exists_check = [
		args.particle_stack,
		args.partres_file,
		args.params_2d_file,
		args.params_3d_file,
		args.params_3d_index_file,
		args.list,
		args.exlist,
		]
	stack_dependency_check = [
		args.params_2d_file,
		args.params_3d_file,
		args.params_3d_index_file,
		args.list,
		args.exlist,
		]

	if not args.particle_stack and not args.partres_file:
			global_def.ERROR(
				'Particle_stack or partres_file option needs to be present!'.format(option),
				'sxsphire2relion',
				1
				)


	for option in stack_dependency_check:
		if option and not args.particle_stack:
			global_def.ERROR(
				'{0} requires particle stack option!'.format(option),
				'sxsphire2relion',
				1
				)

	for idx, option in enumerate(file_exists_check):
		if option:
			if option.startswith('bdb:'):
				if '#' in option:
					raw_dirnames, basename = option.split('#')
					dirnames = raw_dirnames[4:]
				else:
					dirnames = os.path.dirname(option[4:])
					if not dirnames:
						dirnames = '.'
					basename = os.path.basename(option[4:])
				option = '{0}/EMAN2DB/{1}.bdb'.format(dirnames,basename)
			if not os.path.isfile(option):
				global_def.ERROR(
					'{0} stack must exist!'.format(option),
					'sxsphire2relion',
					1
					)

	if args.list and args.exlist:
		global_def.ERROR(
			'Arguments list and exlist cannot be used at the same time.',
			'sxsphire2relion',
			1
			)

	if args.params_2d_file and args.params_3d_file:
		global_def.ERROR(
			'Arguments params_2d_file and params_3d_file cannot be used at the same time.',
			'sxsphire2relion',
			1
			)

	if args.params_3d_index_file and not args.params_3d_file:
		global_def.ERROR(
			'Arguments params_3d_index_file requires params_3d_file to be set.',
			'sxsphire2relion',
			1
			)

	if args.params_3d_chunk_files and not args.params_3d_file:
		global_def.ERROR(
			'Arguments params_3d_chunk_files requires params_3d_file to be set.',
			'sxsphire2relion',
			1
			)

	try:
		os.mkdir(args.output_dir)
	except OSError:
		pass

	output_path = os.path.join(args.output_dir, args.output_name)
	if os.path.exists(output_path) and args.force:
		pass
	elif os.path.exists(output_path) and not args.force:
		global_def.ERROR(
			'Output file {0} must not exist! Use the --force flag to overwrite existing files'.format(output_path),
			'sxsphire2relion',
			1
			)
	else:
		pass

	return output_path


def create_stack_dtype(particle_dict):
	"""
	Create a dtype list based on the particle_dict.

	Argument:
	particle_dict - Dictionary containing the images header information

	Returns:
	Dtype list, header name list
	"""
	original_name = {}
	for key in particle_dict:
		if key == 'ptcl_source_coord':
			original_name[key] = [
				('_rlnCoordinateX', float),
				('_rlnCoordinateY', float),
				]

		elif key == 'ptcl_source_apix':
			original_name[key] = [
				('_rlnDetectorPixelSize', float),
				('_rlnMagnification', float)
				]

		elif key == 'ptcl_source_image':
			original_name[key] = [
				('_rlnMicrographName', '|S1000'),
				]

		elif key == 'data_path':
			original_name[key] = [
				('_rlnImageName', '|S1000'),
				]

		elif key == 'ptcl_source_coord_id':
			original_name[key] = [
				('ptcl_source_coord_id', int)
				]

		elif key == 'filament':
			original_name[key] = [
				('_rlnHelicalTubeID', int)
				]

		elif key == 'ctf':
			original_name[key] = [
				('_rlnDefocusU', float),
				('_rlnDefocusV', float),
				('_rlnDefocusAngle', float),
				('_rlnAmplitudeContrast', float),
				('_rlnVoltage', float),
				('_rlnPhaseShift', float),
				('_rlnSphericalAberration', float)
				]

		elif key == 'xform.align2d':
			original_name[key] = [
				('_rlnOriginX', float),
				('_rlnOriginY', float)
				]

		elif key == 'xform.projection':
			original_name[key] = [
				('_rlnOriginX', float),
				('_rlnOriginY', float),
				('_rlnAngleRot', float),
				('_rlnAngleTilt', float),
				('_rlnAnglePsi', float)
				]

	dtype_list = []
	for key in original_name:
		for dtype in original_name[key]:
			if dtype not in dtype_list:
				dtype_list.append(dtype)
	return sorted(dtype_list), original_name.keys()


def create_stack_array(dtype_list, header_dict, output_dir):
	"""
	Import the bdb file into an numpy array.

	Arguments:
	dtype_list - List containing the dtypes
	header_dict - Dictionary containing the actual header entries.

	Returns:
	Particle array
	"""
	final_dtype_list = [entry for entry in dtype_list if entry[0].startswith('_rln')]
	particle_array = np.empty(len(header_dict.values()[0]), dtype=final_dtype_list)
	create_stack = False

	for key in header_dict:
		if key == 'ptcl_source_coord':
			coord_x = [entry[0] for entry in header_dict[key]]
			coord_y = [entry[1] for entry in header_dict[key]]
			particle_array['_rlnCoordinateX'] = coord_x
			particle_array['_rlnCoordinateY'] = coord_y

		elif key == 'ptcl_source_apix':
			particle_array['_rlnDetectorPixelSize'] = header_dict[key]
			particle_array['_rlnMagnification'] = 10000

		elif key == 'filament':
			data = [int(entry[-5:])+1 for entry in header_dict[key]]
			particle_array['_rlnHelicalTubeID'] = data

		elif key == 'ctf':
			dict_list = [entry.to_dict() for entry in header_dict[key]]
			defocus = np.array([entry['defocus'] for entry in dict_list])
			astigmatism_amp = np.array([entry['dfdiff'] for entry in dict_list])
			particle_array['_rlnDefocusAngle'] = 45 - np.array([entry['dfang'] for entry in dict_list])
			particle_array['_rlnAmplitudeContrast'] = np.array([entry['ampcont'] for entry in dict_list]) / 100
			particle_array['_rlnVoltage'] = np.array([entry['voltage'] for entry in dict_list])
			particle_array['_rlnSphericalAberration'] = np.array([entry['cs'] for entry in dict_list])
			particle_array['_rlnDefocusU'] = (20000 * defocus - 10000 * astigmatism_amp) / 2
			particle_array['_rlnDefocusV'] = 20000 * defocus - particle_array['_rlnDefocusU']
			particle_array['_rlnPhaseShift'] = 0

		elif key == 'xform.projection' or key == 'xform.align2d':
			dict_list = [entry.get_params('mrc') for entry in header_dict[key]]
			particle_array['_rlnOriginX'] = np.array([entry['tx'] for entry in dict_list])
			particle_array['_rlnOriginY'] = np.array([entry['ty'] for entry in dict_list])
			particle_array['_rlnAngleRot'] = np.array([entry['phi'] for entry in dict_list])
			particle_array['_rlnAngleTilt'] = np.array([entry['theta'] for entry in dict_list])
			particle_array['_rlnAnglePsi'] = np.array([entry['omega'] for entry in dict_list])

		elif key == 'data_path':
			if header_dict['data_path'][0].endswith('.mrcs'):
				data = [
					'{0:05d}@{1}'.format(entry1+1, os.path.relpath(entry2, output_dir))
					for entry1, entry2 in
					zip(
						header_dict['ptcl_source_coord_id'],
						header_dict['data_path']
						)
					]
				particle_array['_rlnImageName'] = data

			else:
				create_stack = True
				for name in np.unique(header_dict['ptcl_source_image']):
					mask = header_dict['ptcl_source_image'] == name
					particle_array['_rlnImageName'][mask] = [
						'{0:05d}@{1}s'.format(
							entry+1,
							os.path.relpath(
								os.path.join(
									output_dir,
									'Particles',
									os.path.basename(name)
									),
								output_dir
								)
							)
						for entry in
						np.arange(np.sum(mask))
						]

		elif key == 'ptcl_source_image':
			particle_array['_rlnMicrographName'] = [
				entry if entry.startswith('/')
				else os.path.relpath(entry, output_dir)
				for entry in header_dict[key]
				]

	return particle_array, create_stack


def import_particle_stack(particle_stack, output_dir):
	"""
	Import the particle stack.

	Arguments:
	particle_stack - Path to the particle stack

	Returns:
	Particle array
	"""
	particle_header = EMAN2_cppwrap.EMData()
	particle_header.read_image(
		particle_stack,
		0,
		True
		)

	dtype_list, name_list = create_stack_dtype(particle_header.get_attr_dict())

	header_dict = {}
	for name in name_list:
		header_dict[name] = np.array(EMAN2_cppwrap.EMUtil.get_all_attributes(particle_stack, name))

	stack_array, create_stack = create_stack_array(dtype_list, header_dict, output_dir)

	return stack_array, create_stack


if __name__ == '__main__':
	main(parse_args())
