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
	from global_def import SPARXVERSION, ERROR
	import global_def
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + " stack ref_vol outdir  <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --ynumber=y_numbers  --txs=translational_search_stepx  --delta=angular_step --an=angular_neighborhood --maxit=max_iter --CTF --snr=1.0  --sym=c1 --datasym=symdoc"
	
	parser = OptionParser(usage,version=SPARXVERSION)
	#parser.add_option("--ir",                 type="float", 	     default= -1,                 help="Inner radius for psi angle search > 0 (set to 1) (Angstroms)")
	parser.add_option("--ou",                 type="float", 	     default= -1,                 help="Outer radius for psi angle search < int(nx*pixel_size/2)-1 (Angstroms)")
	parser.add_option("--rs",                 type="int",   		 default= 1,                  help="Step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",                 type="string",		 default= " 4  2 1  1   1",   help="Range for translation search in x direction, search within +/-xr (Angstroms) ")
	parser.add_option("--txs",                type="string",		 default= "1 1 1 0.5 0.25",   help="Step size of the translation search in x directions, search is -xr, -xr+ts, 0, xr-ts, xr (Angstroms)")
	parser.add_option("--y_restrict",         type="string",		 default= "-1 -1 -1 -1 -1",   help="Range for translational search in y-direction, search is +/-y_restrict in Angstroms. This only applies to local search, i.e., when an is not -1. If y_restrict < 0, then the y search range is set such that it is the same ratio to dp as angular search range is to dphi. For regular ihrsr, y search range is the full range when y_restrict< 0. Default is -1.")
	parser.add_option("--ynumber",            type="string",		 default= "4 8 16 32 32",     help="Even number of the steps for the search in y direction, search is (-dpp/2,-dpp/2+dpp/ny,,..,0,..,dpp/2-dpp/ny dpp/2]")
	parser.add_option("--delta",              type="string",		 default= "10 6 4  3  2",     help="Angular step of reference projections")
	parser.add_option("--an",                 type="string",		 default= "-1",               help="Angular neighborhood for local searches")
	parser.add_option("--maxit",              type="int",            default= 30,                 help="Maximum number of iterations performed for each angular step (set to 30) ")
	parser.add_option("--searchit",           type="int",            default= 1,                  help="Number of iterations to predict/search before doing reconstruction and updating of reference volume. Default is 1. If maxit=3 and searchit=2, then for each of the 3 inner iterations, 2 iterations of prediction/search will be performed before generating reconstruction.")
	parser.add_option("--CTF",                action="store_true",   default=False,      		  help="CTF correction")
	parser.add_option("--snr",                type="float",          default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--slowIO",             action="store_true",   default=False,               help="use slowIO version")
	#parser.add_option("--fourvar",           action="store_true",   default=False,               help="compute Fourier variance")
	parser.add_option("--apix",               type="float",			 default= -1.0,               help="Pixel size in Angstroms")   
	parser.add_option("--dp",                 type="float",			 default= -1.0,               help="Helical symmetry axial rise (Angstroms)")   
	parser.add_option("--dphi",               type="float",			 default= -1.0,               help="Helical symmetry azimuthal angle")  
	#parser.add_option("--MA",                 action="store_true",   default=False,      		  help="predict consistent parameters based on moving average")
	
	parser.add_option("--psi_max",            type="float", 		 default= 10.0,               help="Maximum psi - how far rotation in plane can can deviate from 90 or 270 degrees")   
	parser.add_option("--rmin",               type="float", 		 default= 0.0,                help="Min radius for application of helical symmetry (Angstroms)")   
	parser.add_option("--rmax",               type="float", 		 default= 80.0,               help="Max radius for application of helical symmetry (Angstroms)")
	parser.add_option("--fract",              type="float", 		 default= 0.7,                help="Fraction of volume used for application of helical symmetry")
	parser.add_option("--sym",                type="string",		 default= "c1",               help="Point-group symmetry of the filament")
	parser.add_option("--function",           type="string",		 default="helical",  	      help="Name of the reference preparation function (Default: helical)")
	parser.add_option("--npad",               type="int",   		 default= 2,                  help="Padding size for 3D reconstruction (default=2)")
	parser.add_option("--debug",              action="store_true",   default=False,               help="debug")
	parser.add_option("--initial_theta",      type="float",		     default=90.0,                help="Intial theta for out-of-plane tilt search, the range will be (initial theta to 90.0 in steps of delta) (default = 90, no out-of-plane tilt)")
	parser.add_option("--delta_theta",        type="float",		     default=1.0,                 help="Delta theta for out-of-plane tilt search (default = 1)")
	#parser.add_option("--boundaryavg",        action="store_true",   default=False,      		  help="boundaryavg")
	#parser.add_option("--MA_WRAP",            type="int",            default= 0,                  help="do wrapping in MA if MA_WRAP=1, else no wrapping in MA. Default is 0.")
	parser.add_option("--seg_ny",             type="int",            default= 256,                help="y dimension of desired segment size, should be related to fract in that fract ~ seg_ny/ny, where ny is dimension of input projections. (pixels)")
	parser.add_option("--new",                action="store_true",   default=False,               help="use new version")
	parser.add_option("--snake",              action="store_true",   default=False,               help="use snake method")	
	parser.add_option("--snakeknots",         type="int",            default= -1,                 help="maximal number of knots for each filament snake. If take default value -1, it will take nseg//2+1, where nseg is the number of segments in the filament")
	
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 3 or len(args) > 4:
		print "usage: " + usage + "\n"
		print "Please run '" + progname + " -h' for detailed options"
	else:
		global_def.BATCH = True
		# Convert input arguments in the units/format as expected by ihrsr_MPI in applications.
		if options.apix < 0:
			ERROR("Please specify pixel size apix","sxheliconlocal",1)
		if options.dp < 0 or options.dphi < 0:
			ERROR("Please specify helical symmetry parameters dp and dphi","sxheliconlocal",1)
		if options.an <= 0 :
			ERROR("Angular search range (an) has to be given.  Only local searches are permitted.","sxheliconlocal",1)

		print  " This code is under development, some instabilities are possible 12/28/2014"

		rminp = int((float(options.rmin)/options.apix) + 0.5)
		rmaxp = int((float(options.rmax)/options.apix) + 0.5)
		
		from utilities import get_input_from_string, get_im

		xr = get_input_from_string(options.xr)
		txs = get_input_from_string(options.txs)
		y_restrict = get_input_from_string(options.y_restrict)

		irp = 1
		if options.ou < 0:  oup = -1
		else:               oup = int( (options.ou/options.apix) + 0.5)
		xrp = ""
		txsp = ""
		y_restrict2 = ""

		for i in xrange(len(xr)):    xrp += str(float(xr[i])/options.apix)+" "
		xrp = xrp[:-1]
		for i in xrange(len(txs)):  txsp += str(float(txs[i])/options.apix)+" "
		txsp = txsp[:-1]
		# now y_restrict has the same format as x search range .... has to change ihrsr accordingly
		for i in xrange(len(y_restrict)): y_restrict2 +=  str(float(y_restrict[i])/options.apix)+" "
		y_restrict2 = y_restrict2[:-1]

		from mpi import mpi_init, mpi_finalize
		sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications import localhelicon_MPI, localhelicon_MPInew, localhelicon_MPIming
		if len(args) < 4:  mask = None
		else:              mask = args[3]
		if options.new:  localhelicon_MPInew(args[0], args[1], args[2], options.seg_ny, mask, irp, oup, options.rs, xrp, options.ynumber, \
			txsp, options.delta, options.initial_theta, options.delta_theta, options.an, options.maxit, options.CTF, options.snr, \
				options.dp, options.dphi, options.psi_max, \
			rminp, rmaxp, options.fract, options.npad,options.sym, options.function,\
			options.apix, options.debug, y_restrict2, options.searchit, slowIO.options)
		elif options.snake:	localhelicon_MPIming(args[0], args[1], args[2], options.seg_ny, mask, irp, oup, options.rs, xrp, options.ynumber, \
			txsp, options.delta, options.initial_theta, options.delta_theta, options.an, options.maxit, options.CTF, options.snr, \
				options.dp, options.dphi, options.psi_max, \
			rminp, rmaxp, options.fract, options.npad,options.sym, options.function,\
			options.apix, options.debug, y_restrict2, options.searchit, options.snakeknots, slowIO.options)	
		else:  localhelicon_MPI(args[0], args[1], args[2], options.seg_ny, mask, irp, oup, options.rs, xrp, options.ynumber, \
			txsp, options.delta, options.initial_theta, options.delta_theta, options.an, options.maxit, options.CTF, options.snr, \
				options.dp, options.dphi, options.psi_max, \
			rminp, rmaxp, options.fract, options.npad,options.sym, options.function,\
			options.apix, options.debug, y_restrict2, options.searchit, slowIO.options)
		global_def.BATCH = False
	
		from mpi import mpi_finalize
		mpi_finalize()

if __name__ == "__main__":
	main()
