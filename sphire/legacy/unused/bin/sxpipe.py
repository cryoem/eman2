'''1
	if args.generate_mask:
		print(" ")
		print_progress("Generating mooon elimnated soft-edged 3D mask from the 3D binary volume corresponding to the specified molecular mass or density threshold with specified option parameters...")
		print_progress("  Soft-edge type : {}".format(args.edge_type ))
		gm_mask3d_moon_eliminated = None
		gm_bin3d_mol_mass_dilated = None
		gm_soft_edging_start_time = time()
		if args.edge_type == "cosine":
			# Use cosine soft-edged which is same as PostRefiner
			gm_binarize_threshold = 0.5
			if args.debug:
				print_progress("MRK_DEBUG: ")
				print_progress("MRK_DEBUG: Util.adaptive_mask()")
				print_progress("MRK_DEBUG:   gm_binarize_threshold  := {}".format(gm_binarize_threshold))
				print_progress("MRK_DEBUG:   args.gm_dilation       := {}".format(args.gm_dilation))
				print_progress("MRK_DEBUG:   args.gm_edge_width     := {}".format(args.gm_edge_width))
			gm_mask3d_moon_eliminated = Util.adaptive_mask(bin3d_mol_mass, gm_binarize_threshold, args.gm_dilation, args.gm_edge_width)
		elif args.edge_type == "gauss":
			gm_gauss_kernel_radius = args.gm_edge_width - args.gm_dilation
			gm_mask3d_moon_eliminated, gm_bin3d_mol_mass_dilated = mrk_sphere_gauss_edge(bin3d_mol_mass, args.gm_dilation, gm_gauss_kernel_radius, args.gm_edge_sigma, args.debug)
		else:
		print_progress("  Totally, {} soft-edging of 3D mask took {:7.2f} sec...".format(args.edge_type.upper(), time() - gm_soft_edging_start_time))
		
		if args.debug:
			if gm_bin3d_mol_mass_dilated is not None:
				gm_bin3d_mol_mass_dilated_file_path = os.path.join(args.output_directory, "mrkdebug{:02d}_gm_bin3d_mol_mass_dilated.hdf".format(debug_output_id))
				gm_bin3d_mol_mass_dilated.write_image(gm_bin3d_mol_mass_dilated_file_path)
				debug_output_id += 1
		
		gm_mask3d_moon_eliminated_file_path = os.path.join(args.output_directory, "{}_mask_moon_eliminated.hdf".format(args.outputs_root))
		print(" ")
		print_progress("Saving moon eliminated 3D mask to {}...".format(gm_mask3d_moon_eliminated_file_path))
		gm_mask3d_moon_eliminated.write_image(gm_mask3d_moon_eliminated_file_path)
	'''	
'''2
	print(" ")
	print_progress("Summary of processing...")
	print_progress("  Provided expected molecular mass [kDa]      : {}".format(args.mol_mass))
	"""
	print_progress("  Applied density threshold                   : {}".format(density_threshold))
	print_progress("  Computed molecular mass [kDa] of density    : {}".format(computed_mol_mass_from_density))
	print_progress("  Percentage of this molecular mass [kDa]     : {}".format(computed_mol_mass_from_density/args.mol_mass * 100))
	print_progress("  Computed volume [Voxels] of density         : {}".format(computed_vol_voxels_from_density))
	print_progress("  Percentage of this volume [Voxels]          : {}".format(computed_vol_voxels_from_density/computed_vol_voxels_from_mass * 100))
	"""
	print_progress("  Saved moon eliminated 3D reference to       : {}".format(ref3d_moon_eliminated_file_path))
	if args.generate_mask:
		print_progress("  Saved mooon elimnated soft-edged 3D mask to : {}".format(gm_mask3d_moon_eliminated_file_path))
	'''
"""3
	lina = numpy.argsort(radius_array)
	sorted_radius = radius_array[lina[::-1]]
	array_x = numpy.arange(sorted_radius.shape[0])
	angles_no_mirror = angles_no_mirror[lina[::-1]]
	nonzero_mask = list(nonzero_mask[0][lina[::-1]])

	"""
"""4
	print(array_x)
	print(sorted_radius)
	print(len(angles_no_mirror))
	print(angles_no_mirror)
	"""
