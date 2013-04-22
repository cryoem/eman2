#!/usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

# clean up the code, make documentation

import os
import global_def
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ Input Output [options]
	
	Helicise the atom coordinates of input pdb file according to input helical symmetry parameters.
	Input: pdb file containing atom coordinates to be helicised and helical symmetry parameters dp and dphi. 
	Output: pdb file containing helicised atom coordinates
	
		sxhelical_demo.py 3MFP_1SU.pdb rnew.pdb --heli --dp=27.6 --dphi=166.5 \n
	
	Generate projections of helical filament from reference volume. 
	Input: Reference Volume, output directory 
	Output: Three micrographs stored in output directory		
				 
		sxhelical_demo.py tmp.hdf  mic --generate_micrograph --CTF --apix=1.84	
	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	# helicise the Atom coordinates
	parser.add_option("--heli",                   action="store_true",      default=False,      		  	 help="Helicise the atom coordinates of input pdb file according to input helical symmetry parameters. \n Input: pdb file containing atom coordinates to be helicised and helical symmetry parameters dp and dphi. \n Output: pdb file containing helicised atom coordinates")
	parser.add_option("--dp",                     type="float",			    default= -1.0,              	 help="delta z - translation in Angstroms")   
	parser.add_option("--dphi",                   type="float",			    default= -1.0,              	 help="delta phi - rotation in degrees")  
	
	# generate micrographs of helical filament
	parser.add_option("--generate_micrograph",    action="store_true",      default=False,      		  	 help="Generate micrograph of helical filament from reference volume. Input: Reference Volume, output directory, Output: Three micrographs containing helical filaments stored in output directory")
	parser.add_option("--CTF",              	  action="store_true",  	default=False,   				 help="Use CTF correction")
	parser.add_option("--apix",               	  type="float",			 	default= -1,               	     help="pixel size in Angstroms")   
	parser.add_option("--rand_seed",              type="int",			    default=14567,              	 help="the seed used for generating random numbers (default 14567) for adding noise to the generated micrographs.")
	parser.add_option("--Cs",               	  type="float",			 	default= 2.0,               	 help="Microscope Cs (spherical aberation)")
	parser.add_option("--voltage",				  type="float",				default=200.0, 					 help="Microscope voltage in KV")
	parser.add_option("--ac",					  type="float",				default=10.0, 					 help="Amplitude contrast (percentage, default=10)")
	
	# generate initial volume
	parser.add_option("--generate_initvol",       action="store_true",      default=False,      		  	 help="Generate initial volume of noisy cylinder.")
	parser.add_option("--boxsize",                type="string",		    default="100,100,200",           help="String containing x , y, z dimensions (separated by comma) in pixels")
	parser.add_option("--rad",                    type="int",			    default=35,              	 	 help="Radius of initial volume in pixels")
	
	(options, args) = parser.parse_args()
	if len(args) > 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		if options.heli:
			if options.dp < 0 or options.dphi < 0:
				print "Please enter helical symmetry parameters dp and dphi."
				sys.exit()
			helicise_pdb(args[0], args[1], options.dp, options.dphi)
		
		if options.generate_micrograph:
			if options.apix <= 0:
				print "Please enter pixel size."
				sys.exit()
			generate_helimic(args[0], args[1], options.apix, options.CTF, options.Cs, options.voltage, options.ac, options.rand_seed)
		
		if options.generate_initvol:
			from utilities import model_cylinder, model_gauss_noise
			outvol = args[0]
			boxdims = options.boxsize.split(',')
			if len(boxdims) < 1 or len(boxdims) > 3:
				print "Enter box size as string containing x , y, z dimensions (separated by comma) in pixels. E.g.: --boxsize='100,100,200'"
				sys.exit()
			nx= int(boxdims[0])
			if len(boxdims) == 1:
				ny = nx
				nz = nx
			else:
				ny = int(boxdims[1])
				if len(boxdims) == 3:
					nz = int(boxdims[2])
					
			(model_cylinder(options.rad,nx, ny, nz)*model_gauss_noise(1.0, nx, ny, nz) ).write_image(outvol)
		
			
def helicise_pdb(inpdb, outpdb, dp, dphi):
	from math import cos, sin, pi
	from copy import deepcopy
	from numpy import zeros, float32, dot

	dp = dp*-1.0
	infile =open(inpdb,"r")
	pall = infile.readlines()
	infile.close()

	p = []

	pos = []
	for i in xrange( len(pall) ):
		
		if( (pall[i])[:4] == 'ATOM'):
			p.append( pall[i] )
			pos.append(i)
	n = len(p)
	nperiod = 50
	X = zeros( (3,len(p) ), dtype=float32 )
	X_new = zeros( (3,len(p) ), dtype=float32 )
	for i in xrange( len(p) ):
	
		element = deepcopy( p[i] )
		X[0,i]=float(element[30:38])
		X[1,i]=float(element[38:46])	
		X[2,i]=float(element[46:54])
	
	pnew = []
	
	for j in xrange(-nperiod, nperiod+1):
		for i in xrange( n ):
			pnew.append( deepcopy(p[i]) )
	
	for j in xrange(-nperiod, nperiod+1):
		if j != 0:
			rd = pi*j*dphi/180
			m = zeros( (3,3 ), dtype=float32 )
			t = zeros( (3,1 ), dtype=float32 )
			m[0][0] = cos(rd)
			m[0][1] = -sin(rd)
			m[1][0] = sin(rd)
			m[1][1] = cos(rd)
			m[2][2] = 1.0
			t[0,0]=0.0
			t[1,0]=0.0
			t[2,0]=j*dp
			X_new = dot(m, X) + t
			for i in xrange( n ):
				pnew[j*n+i] = pnew[j*n+i].replace( p[i][30:38], "%8.3f"%( float(X_new[0,i]) ) )
				pnew[j*n+i] = pnew[j*n+i].replace( p[i][38:46], "%8.3f"%( float(X_new[1,i]) ) )
				pnew[j*n+i] = pnew[j*n+i].replace( p[i][46:54], "%8.3f"%( float(X_new[2,i]) ) )
	
	outfile=open(outpdb,"w")
	outfile.writelines(pall[0:pos[0]-1])
	outfile.writelines(pnew)
	outfile.writelines(pall[n-1:len(pall)])
	outfile.close()

def generate_helimic(refvol, outdir, pixel, CTF=False, Cs=2.0,voltage = 200.0, ampcont = 10.0,rand_seed=14567):
	
	from utilities	 import model_blank, model_gauss, model_gauss_noise, pad
	from random 	 import random
	from projection  import prgs, prep_vol
	from filter	     import filt_gaussl, filt_ctf
	from EMAN2 	     import EMAN2Ctf
	
	if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', "sxhelical_demo", 1)
	os.mkdir(outdir)
	seed(14567)
	angles =[]
	for i in xrange(3):
		angles.append( [0.0+60.0*i, 90.0-i*5, 90.0, 0.0, 0.0] )
	
	nangle = len(angles)
	modelvol = EMData()
	modelvol.read_image(refvol)
	
	nx = modelvol.get_xsize()
	ny = modelvol.get_ysize()
	nz = modelvol.get_zsize()
	
	iprj    = 0
	width = 500
	xstart = 0
	ystart = 0
	
	for idef in xrange(3,6):
	
		mic = model_blank(2048, 2048)
		defocus = idef*0.2
		if CTF :
			ctf = EMAN2Ctf()
			ctf.from_dict( {"defocus":defocus, "cs":Cs, "voltage":voltage, "apix":pixel, "ampcont":ampcont, "bfactor":0.0} )
	
		sigma = 1.5 + random() # 1.5-2.5
		addon = model_gauss(sigma, nx, ny, nz, sigma, sigma, nx/2, ny/2, nz/2 )
		scale = 2500 * (0.5+random())
		model = modelvol + scale*addon
		volfts, kbx, kby, kbz = prep_vol( model )
		i = idef - 3
		for k in xrange(1):
			dphi = 0.0 #8.0*(random()-0.5)
			dtht = 0.0 #6.0*(random()-0.5)
			psi  = 90+10*( i-1 )
	
			phi = angles[idef-3][0]+dphi
			tht = angles[idef-3][1]-dtht
	
			s2x = 2.5*(i-1)
			s2y = 2.5*(i-1)
	
	
			proj = prgs(volfts, kbz, [phi, tht, psi, -s2x, -s2y], kbx, kby)
			proj = Util.window(proj, 320, nz)		
			mic += pad(proj, 2048, 2048, 1, 0.0, int(750*(i-1)), int(20*(i-1)), 0)
	
		mic += model_gauss_noise(30.0,2048,2048)
		if CTF :
			#apply CTF
			mic = filt_ctf(mic, ctf)
	
		mic += filt_gaussl(model_gauss_noise(17.5,2048,2048), 0.3)
		
		mic.write_image("%s/mic%1d.hdf"%(outdir, idef-3),0)
	
if __name__ == "__main__":
	main()
