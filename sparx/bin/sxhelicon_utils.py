#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
#


def main():
	import os
	import sys
	from optparse import OptionParser
	from global_def import SPARXVERSION
	import global_def
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage2 = progname + """ inputfile outputfile [options]
        Examples:
		
		Generate disks from filament based reconstructions:
		
			sxheader.py stk.hdf --params=xform.projection --import=params.txt
			sxheader.py stk.hdf --params=active --one
			mpirun -np 2 sxhelicon_utils.py stk.hdf --gendisk='bdb:disk' --ref_nx=100 --ref_ny=100 --ref_nz=200 --apix=1.84 --dp=27.6 --dphi=166.715 --fract=0.67 --rmin=0 --rmax=64 --function="[.,nofunc,helical3c]" --sym="c1" --MPI

        Helicise input volume and save the result to output volume:
        
            sxhelicon_utils.py input_vol.hdf output_vol.hdf --helicise --dp=27.6 --dphi=166.5 --fract=0.65 --rmax=70 --rmin=1 --apix=1.84         
			
            sxhelicon_utils.py bdb:big_stack --hfsc='flst_' --filament_attr=filament
			
            sxhelicon_utils.py bdb:big_stack --predict_helical=helical_params.txt --dp=27.6 --dphi=166.5 --apix=1.84
            
            sxhelicon_utils.py disk_to_stack.hdf --stackdisk=stacked_disks.hdf --dphi=166.5 --dp=27.6 --ref_nx=160 --ref_ny=160 --ref_nz=225 --apix=1.84
"""
	parser = OptionParser(usage2,version=SPARXVERSION)
	#parser.add_option("--ir",                 type="float", 	     default= -1,                 help="inner radius for rotational correlation > 0 (set to 1) (Angstroms)")
	parser.add_option("--ou",                 type="float", 	     default= -1,                 help="outer radius for rotational 2D correlation < int(nx/2)-1 (set to the radius of the particle) (Angstroms)")
	parser.add_option("--rs",                 type="int",   		 default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",                 type="string",		 default= " 4  2 1  1   1",   help="range for translation search in x direction, search is +/-xr (Angstroms) ")
	parser.add_option("--txs",                type="string",		 default= "1 1 1 0.5 0.25",   help="step size of the translation search in x directions, search is -xr, -xr+ts, 0, xr-ts, xr (Angstroms)")
	parser.add_option("--delta",              type="string",		 default= " 10 6 4  3   2",   help="angular step of reference projections")
	parser.add_option("--an",                 type="string",		 default= "-1",               help="angular neighborhood for local searches")
	parser.add_option("--maxit",              type="int",            default= 30,                 help="maximum number of iterations performed for each angular step (set to 30) ")
	parser.add_option("--CTF",                action="store_true",   default=False,      		  help="CTF correction")
	parser.add_option("--snr",                type="float",          default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--MPI",                action="store_true",   default=False,               help="use MPI version")
	#parser.add_option("--fourvar",           action="store_true",   default=False,               help="compute Fourier variance")
	parser.add_option("--apix",               type="float",			 default= -1.0,               help="pixel size in Angstroms")   
	parser.add_option("--dp",                 type="float",			 default= -1.0,               help="delta z - translation in Angstroms")   
	parser.add_option("--dphi",               type="float",			 default= -1.0,               help="delta phi - rotation in degrees")  
		  
	parser.add_option("--rmin",               type="float", 		 default= 0.0,                help="minimal radius for hsearch (Angstroms)")   
	parser.add_option("--rmax",               type="float", 		 default= 80.0,               help="maximal radius for hsearch (Angstroms)")
	parser.add_option("--fract",              type="float", 		 default= 0.7,                help="fraction of the volume used for helical search")
	parser.add_option("--sym",                type="string",		 default= "c1",               help="symmetry of the structure")
	parser.add_option("--function",           type="string",		 default="helical",  	      help="name of the reference preparation function")
	parser.add_option("--npad",               type="int",   		 default= 2,                  help="padding size for 3D reconstruction")
	parser.add_option("--debug",              action="store_true",   default=False,               help="debug")
	
	parser.add_option("--volalixshift",       action="store_true",   default=False,               help="Use volalixshift refinement")
	parser.add_option("--searchxshift",       type="float",		     default= 0.0,                help="search range for x-shift determination: +/- searchxshift (Angstroms)")
	parser.add_option("--nearby",             type="float",		     default= 6.0,                help="neighborhood within which to search for peaks in 1D ccf for x-shift search (Angstroms)")

	# diskali
	parser.add_option("--diskali",            action="store_true",   default=False,               help="volume alignment")
	parser.add_option("--zstep",              type="float",          default= 1,                  help="Step size for translational search along z (Angstroms)")   

	# helicise
	parser.add_option("--helicise",           action="store_true",	 default=False,               help="helicise input volume and save results to output volume")
	parser.add_option( "--hfsc",              type="string",      	 default="",                  help="generate two list of image indices used to split segment stack into two for helical fsc calculation. The two lists will be stored in two text files named using file_prefix with '_even' and '_odd' suffixes respectively." )
	parser.add_option( "--filament_attr",     type="string",      	 default="filament",          help="attribute under which filament identification is stored" )
	parser.add_option( "--predict_helical",   type="string",      	 default="",                  help="Generate projection parameters consistent with helical symmetry")

	# input options for generating disks
	parser.add_option("--gendisk",            type="string",		 default="",                  help="Name of file under which generated disks will be saved to") 
	parser.add_option("--ref_nx",             type="int",   		 default= 1,                  help="nx=ny volume size" ) 
	parser.add_option("--ref_nz",             type="int",   		 default= 1,                  help="nz volume size - computed disks will be nx x ny x rise/apix" ) 
	parser.add_option("--new_pixel_size",           type="float", 		 default= -1,                 help="desired pixel size of the output disks. The default is -1, in which case there is no resampling (unless --match_pixel_rise flag is True).")
	parser.add_option("--maxerror",           type="float", 		 default= 0.1,                help="proportional to the maximum amount of error to tolerate between (dp/new_pixel_size) and int(dp/new_pixel_size ), where new_pixel_size is the pixel size calculated when the option --match_pixel_rise flag is True.")
	parser.add_option("--match_pixel_rise",   action="store_true",	 default=False,               help="calculate new pixel size such that the rise is approximately integer number of pixels given the new pixel size. This will be the pixel size of the output disks.")
	
	# get consistency
	parser.add_option("--consistency",        type="string",		 default="",                  help="Name of parameters to get consistency statisticsf for") 
	parser.add_option("--phithr",             type="float", 		 default= 2.0,                help="phi threshold for consistency check")  
	parser.add_option("--ythr",               type="float", 		 default= 2.0,                help="y threshold (in Angstroms) for consistency check")  
	parser.add_option("--segthr",             type="int", 		     default= 3,                  help="minimum number of segments/filament for consistency check")  

	# stack disks
	parser.add_option("--stackdisk",          type="string",		 default="",                  help="Name of file under which output volume will be saved to.")
	parser.add_option("--ref_ny",             type="int",   		 default=-1,                  help="ny of output volume size. Default is ref_nx" ) 

	# symmetry search
	parser.add_option("--symsearch",          action="store_true",	 default=False, 	  	      help="Do helical symmetry search." ) 
	parser.add_option("--ndp",                type="int",            default= 12,                 help="In symmetrization search, number of delta z steps equas to 2*ndp+1") 
	parser.add_option("--ndphi",              type="int",            default= 12,                 help="In symmetrization search,number of dphi steps equas to 2*ndphi+1")  
	parser.add_option("--dp_step",            type="float",          default= 0.1,                help="delta z (Angstroms) step  for symmetrization")  
	parser.add_option("--dphi_step",          type="float",          default= 0.1,                help="dphi step for symmetrization")
	parser.add_option("--datasym",            type="string",		 default="datasym.txt",       help="symdoc")
	
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1 or len(args) > 5:
		print "Various helical reconstruction related functionalities: " + usage2
		print "Please run '" + progname + " -h' for detailed options"
	else:
		
		if len(options.hfsc) > 0:
			if len(args) != 1:
				print  "Incorrect number of parameters"
				sys.exit()
			from applications import imgstat_hfsc
			imgstat_hfsc( args[0], options.hfsc, options.filament_attr)
			sys.exit()
		
		# Convert input arguments in the units/format as expected by ihrsr_MPI in applications.
		if options.apix < 0:
			print "Please enter pixel size"
			sys.exit()

		if len(options.stackdisk) > 0:
			if len(args) != 1:
				print  "Incorrect number of parameters"
				sys.exit()
			dpp = (float(options.dp)/options.apix)
			rise = int(dpp)
			if(float(rise) != dpp):
				print "  dpp has to be integer multiplicity of the pixel size"
				sys.exit()
			from utilities import get_im
			v = get_im(args[0])
			from applications import stack_disks
			ref_ny = options.ref_ny
			if ref_ny < 0:
				ref_ny = options.ref_nx
			sv = stack_disks(v, options.ref_nx, ref_ny, options.ref_nz, options.dphi, rise)
			sv.write_image(options.stackdisk)
			sys.exit()
		
		if len(options.consistency) > 0:
			if len(args) != 1:
				print  "Incorrect number of parameters"
				sys.exit()
			from development import consistency_params	
			consistency_params(args[0], options.consistency, options.dphi, options.dp, options.apix,phithr=options.phithr, ythr=options.ythr, THR=options.segthr)
			sys.exit()

		rminp = int((float(options.rmin)/options.apix) + 0.5)
		rmaxp = int((float(options.rmax)/options.apix) + 0.5)
		
		from utilities import get_input_from_string, get_im

		xr = get_input_from_string(options.xr)
		txs = get_input_from_string(options.txs)

		irp = 1
		if options.ou < 0:  oup = -1
		else:               oup = int( (options.ou/options.apix) + 0.5)
		xrp = ''
		txsp = ''
		
		for i in xrange(len(xr)):
			xrp += " "+str(float(xr[i])/options.apix)
		for i in xrange(len(txs)):
			txsp += " "+str(float(txs[i])/options.apix)

		searchxshiftp = int( (options.searchxshift/options.apix) + 0.5)
		nearbyp = int( (options.nearby/options.apix) + 0.5)
		zstepp = int( (options.zstep/options.apix) + 0.5)

		if options.MPI:
			from mpi import mpi_init, mpi_finalize
			sys.argv = mpi_init(len(sys.argv), sys.argv)

		if len(options.predict_helical) > 0:
			if len(args) != 1:
				print  "Incorrect number of parameters"
				sys.exit()
			if options.dp < 0:
				print "Helical symmetry paramter rise --dp should not be negative"
				sys.exit()
			from applications import predict_helical_params
			predict_helical_params(args[0], options.dp, options.dphi, options.apix, options.predict_helical)
			sys.exit()

		if options.helicise:	
			if len(args) != 2:
				print "Incorrect number of parameters"
				sys.exit()
			if options.dp < 0:
				print "Helical symmetry paramter rise --dp should not be negative"
				sys.exit()
			from utilities import get_im, sym_vol
			vol = get_im(args[0])
			vol = sym_vol(vol, options.sym)
			hvol = vol.helicise(options.apix, options.dp, options.dphi, options.fract, rmaxp, rminp)
			hvol = sym_vol(hvol, options.sym)
			hvol.write_image(args[1])
			sys.exit()


		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.volalixshift:
			if options.maxit > 1:
				print "Inner iteration for x-shift determinatin is restricted to 1"
				sys.exit()
			if len(args) < 4:  mask = None
			else:               mask = args[3]
			from applications import volalixshift_MPI
			global_def.BATCH = True
			volalixshift_MPI(args[0], args[1], args[2], searchxshiftp, options.apix, options.dp, options.dphi, options.fract, rmaxp, rminp, mask, options.maxit, options.CTF, options.snr, options.sym,  options.function, options.npad, options.debug, nearbyp)
			global_def.BATCH = False

		if options.diskali:
			#if options.maxit > 1:
			#	print "Inner iteration for disk alignment is restricted to 1"
			#	sys.exit()
			if len(args) < 4:  mask = None
			else:               mask = args[3]
			global_def.BATCH = True
			if(options.sym[:1] == "d" or options.sym[:1] == "D" ):
				from development import diskaliD_MPI
				diskaliD_MPI(args[0], args[1], args[2], mask, options.dp, options.dphi, options.apix, options.function, zstepp, options.fract, rmaxp, rminp, options.CTF, options.maxit, options.sym)
			else:
				from applications import diskali_MPI
				diskali_MPI(args[0], args[1], args[2], mask, options.dp, options.dphi, options.apix, options.function, zstepp, options.fract, rmaxp, rminp, options.CTF, options.maxit, options.sym)
			global_def.BATCH = False
		
		if options.symsearch:
			from development import symsearch_MPI
			if len(args) < 3:	
				mask = None
			else:
				mask= args[2]
			global_def.BATCH = True
			symsearch_MPI(args[0], args[1], mask, options.dp, options.ndp, options.dp_step, options.dphi, options.ndphi, options.dphi_step, rminp, rmaxp, options.fract, options.sym, options.function, options.datasym, options.apix, options.debug)
			global_def.BATCH = False
			
		elif len(options.gendisk)> 0:
			from applications import gendisks_MPI
			global_def.BATCH = True
			if len(args) == 1:  mask3d = None
			else:               mask3d = args[1]
			if options.dp < 0:
				print "Helical symmetry paramter rise --dp must be explictly set!"
				sys.exit()
			gendisks_MPI(args[0], mask3d, options.ref_nx, options.ref_nx, options.ref_nz, options.apix, options.dp, options.dphi, options.fract, rmaxp, rminp, options.CTF, options.function, options.sym, options.gendisk, options.maxerror, options.new_pixel_size, options.match_pixel_rise)
			global_def.BATCH = False
		
		if options.MPI:
			from mpi import mpi_finalize
			mpi_finalize()

if __name__ == "__main__":
	main()
