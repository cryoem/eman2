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
	usage = progname + " stack ref_vol outdir  <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --ynumber=y_numbers  --txs=translational_search_stepx  --delta=angular_step --an=angular_neighborhood --center=1 --maxit=max_iter --CTF --snr=1.0  --ref_a=S --sym=c1 --datasym=symdoc --new"
	
	parser = OptionParser(usage,version=SPARXVERSION)
	#parser.add_option("--ir",                 type="float", 	     default= -1,                 help="inner radius for rotational correlation > 0 (set to 1) (Angstroms)")
	parser.add_option("--ou",                 type="float", 	     default= -1,                 help="outer radius for rotational 2D correlation < int(nx/2)-1 (set to the radius of the particle) (Angstroms)")
	parser.add_option("--rs",                 type="int",   		 default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",                 type="string",		 default= " 4  2 1  1   1",   help="range for translation search in x direction, search is +/-xr (Angstroms) ")
	parser.add_option("--txs",                type="string",		 default= "1 1 1 0.5 0.25",   help="step size of the translation search in x directions, search is -xr, -xr+ts, 0, xr-ts, xr (Angstroms)")
	parser.add_option("--y_restrict",         type="string",		 default= "-1 -1 -1 -1 -1",   help="range for translational search in y-direction, search is +/-y_restrict in Angstroms. This only applies to local search, i.e., when an is not -1. If y_restrict < 0, then for ihrsrlocalcons (option --localcons local search with consistency), the y search range is set such that it is the same ratio to dp as angular search range is to dphi. For regular ihrsr, y search range is the full range when y_restrict< 0. Default is -1.")
	parser.add_option("--ynumber",            type="string",		 default= "4 8 16 32 32",     help="even number of the translation search in y direction, search is (-dpp/2,-dpp/2+dpp/ny,,..,0,..,dpp/2-dpp/ny dpp/2]")
	parser.add_option("--delta",              type="string",		 default= " 10 6 4  3   2",   help="angular step of reference projections")
	parser.add_option("--an",                 type="string",		 default= "-1",               help="angular neighborhood for local searches (default -1, meaning do exhaustive search)")
	parser.add_option("--maxit",              type="int",            default= 30,                 help="maximum number of iterations performed for each angular step (default 30) ")
	parser.add_option("--CTF",                action="store_true",   default=False,      		  help="CTF correction")
	parser.add_option("--snr",                type="float",          default= 1.0,                help="Signal-to-Noise Ratio of the data (default 1)")	
	parser.add_option("--MPI",                action="store_true",   default=True,               help="use MPI version")
	#parser.add_option("--fourvar",           action="store_true",   default=False,               help="compute Fourier variance")
	parser.add_option("--apix",               type="float",			 default= -1.0,               help="pixel size in Angstroms")   
	parser.add_option("--dp",                 type="float",			 default= -1.0,               help="delta z - translation in Angstroms")   
	parser.add_option("--dphi",               type="float",			 default= -1.0,               help="delta phi - rotation in degrees")  
		
	parser.add_option("--ndp",                type="int",            default= 12,                 help="In symmetrization search, number of delta z steps equals to 2*ndp+1") 
	parser.add_option("--ndphi",              type="int",            default= 12,                 help="In symmetrization search,number of dphi steps equas to 2*ndphi+1")  
	parser.add_option("--dp_step",            type="float",          default= 0.1,                help="delta z (Angstroms) step  for symmetrization")  
	parser.add_option("--dphi_step",          type="float",          default= 0.1,                help="dphi step for symmetrization")
	   
	parser.add_option("--psi_max",            type="float", 		 default= 10.0,               help="maximum psi - how far rotation in plane can can deviate from 90 or 270 degrees (default 10)")   
	parser.add_option("--rmin",               type="float", 		 default= 0.0,                help="minimal radius for hsearch (Angstroms)")   
	parser.add_option("--rmax",               type="float", 		 default= 80.0,               help="maximal radius for hsearch (Angstroms)")
	parser.add_option("--fract",              type="float", 		 default= 0.7,                help="fraction of the volume used for helical search")
	parser.add_option("--sym",                type="string",		 default= "c1",               help="symmetry of the structure")
	parser.add_option("--function",           type="string",		 default="helical",  	      help="name of the reference preparation function")
	parser.add_option("--datasym",            type="string",		 default="datasym.txt",       help="symdoc")
	parser.add_option("--nise",               type="int",   		 default= 200,                help="start symmetrization after nise steps (default 200)")
	parser.add_option("--npad",               type="int",   		 default= 2,                  help="padding size for 3D reconstruction, (default 2)")
	parser.add_option("--debug",              action="store_true",   default=False,               help="debug")
	parser.add_option("--new",                action="store_true",   default=False,               help="use rectangular recon and projection version")
	parser.add_option("--initial_theta",      type="float",		     default=90.0,                help="intial theta for reference projection (default 90)")
	parser.add_option("--delta_theta",        type="float",		     default=1.0,                 help="delta theta for reference projection (default 1.0)")
	parser.add_option("--WRAP",               type="int",  		     default= 1,                  help="do helical wrapping (default 1, meaning yes)")

	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1 or len(args) > 5:
		print "usage: " + usage + "\n"
		print "Please run '" + progname + " -h' for detailed options"
	else:
		# Convert input arguments in the units/format as expected by ihrsr_MPI in applications.
		if options.apix < 0:
			print "Please enter pixel size"
			sys.exit()

		rminp = int((float(options.rmin)/options.apix) + 0.5)
		rmaxp = int((float(options.rmax)/options.apix) + 0.5)
		
		from utilities import get_input_from_string, get_im

		xr = get_input_from_string(options.xr)
		txs = get_input_from_string(options.txs)
		y_restrict = get_input_from_string(options.y_restrict)

		irp = 1
		if options.ou < 0:  oup = -1
		else:               oup = int( (options.ou/options.apix) + 0.5)
		xrp = ''
		txsp = ''
		y_restrict2 = ''
		
		for i in xrange(len(xr)):
			xrp += " "+str(float(xr[i])/options.apix)
		for i in xrange(len(txs)):
			txsp += " "+str(float(txs[i])/options.apix)
		# now y_restrict has the same format as x search range .... has to change ihrsr accordingly
		for i in xrange(len(y_restrict)):
			y_restrict2 += " "+str(float(y_restrict[i])/options.apix)

		if options.MPI:
			from mpi import mpi_init, mpi_finalize
			sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications import ihrsr
		global_def.BATCH = True
		if len(args) < 4:  mask = None
		else:               mask = args[3]
		ihrsr(args[0], args[1], args[2], mask, irp, oup, options.rs, xrp, options.ynumber, txsp, options.delta, options.initial_theta, options.delta_theta, options.an, options.maxit, options.CTF, options.snr, options.dp, options.ndp, options.dp_step, options.dphi, options.ndphi, options.dphi_step, options.psi_max, rminp, rmaxp, options.fract, options.nise, options.npad,options.sym, options.function, options.datasym, options.apix, options.debug, options.MPI, options.WRAP, y_restrict2) 
		global_def.BATCH = False

		if options.MPI:
			from mpi import mpi_finalize
			mpi_finalize()

if __name__ == "__main__":
	main()
