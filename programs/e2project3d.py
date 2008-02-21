#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Additional Author: David Woolford 2007-2008 (woolford@bcm.edu)
#
# Copyright (c) 2000-2006 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# initial version of project3d; more features/documentation to come
# try "e2project3d --help"

# References
# 1. Baldwin, P.R. and Penczek, P.A. 2007. The Transform Class in SPARX and EMAN2. J. Struct. Biol. 157, 250-261.
# 2. http://blake.bcm.edu/emanwiki/EMAN2/Symmetry

# 
# 1. Asymmetric units are accurately demarcated and covered by the projection algorithm. This can
# be tested using 3D plotting in Matlab.
# 2. By default the entire asymmetric is projected over, but excluding the mirror portion is supported.
# The accuracy of the demarcation of the asymmetric unit can be tested by using the --verifymirror
# argument, which subtracts (mirrored) mirror projections in the local asymmetric unit from the 
# the equivalent original projections - the result should be zero but this is not the case due
# to interpolation differences, however the images should have a mean of about zero and a standard
# deviation that is relatively small. Visual inspection of the results (in default output files) can also help.
# 3. Random orientation generation is supported - all three euler angles are random
# 4. Perturbation of projections generated in the asymmetric unit is supported in regular reconstructions runs.
# 5. The user is able to generate projections that include in-plane or "phi" rotations
# and this is achieved using the phitoo argument
# 6. The user can smear in-plane projections when the "smear" argument is specified in addition to the "phitoo" argument
# what else..????

import sys, math, os, random
from EMAN2 import *
from optparse import OptionParser
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi
DEBUG = False
WEN_JIANG = False
EMAN1_OCT = False
MIRROR_DEBUG = True
NO_MIRROR = False

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog image [options] 
	Projects in real space over the asymmetric unit using the angular separation as specified by prop and the symmetry as specified by sym."""
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos")
	parser.add_option("--orientgen", dest="orientgen", help="The orientation generator to use. See e2help.py orientgen")
	parser.add_option("--out", dest = "outfile", default = "e2proj.img", help = "Output file. Default is 'e2proj.img'")
	# add --perturb
	parser.add_option("--smear", dest = "smear", type = "int", default=0,help="Used in conjunction with --phitoo, this will rotationally smear between phi steps. The user must specify the amount of smearing (typically 2-10)")
	parser.add_option("--projector", dest = "projector", default = "standard",help = "Projector to use")
	#parser.add_option("--verifymirror",action="store_true",help="Used for testing the accuracy of mirror projects",default=False)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_option("--append", "-a",dest="append",default=False, action="store_true",help="Append to the output file")
	parser.add_option("--verbose","-v", dest="verbose", default=False, action="store_true",help="Toggle verbose mode - prints extra infromation to the command line while executing")
	parser.add_option("--check","-c", default=False, action="store_true",help="Checks to see if the command line arguments will work.")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	
	(options, args) = parser.parse_args()
	
	if ( options.check ): options.verbose = True
	
	if len(args) < 1:
		parser.error("Error: No input file given")
	
	options.model = args[0]
	error = check(options,True)
	
	if ( options.verbose ):
		if (error):
			print "e2project3d.py command line arguments test.... FAILED"
		else:
			if (options.verbose):
				print "e2project3.py command line arguments test.... PASSED"
			
	if error : exit(1)
	if options.check: exit(0)
	
	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if ( os.path.exists(options.outfile )):
		if ( options.force ):
			remove_file(options.outfile)

	logger=E2init(sys.argv)
	eulers = []
	
	data = EMData()
	data.read_image(args[0])
	
	sym_object = parsesym(options.sym)
	[og_name,og_args] = parsemodopt(options.orientgen)
	eulers = sym_object.gen_orientations(og_name, og_args)
		
	# generate and save all the projections to disk - that's it, that main job is done
	if ( options.verbose ):
		print "Generating and saving projections..."
	generate_and_save_projections(options, data, eulers, options.smear)
	
	if ( options.verbose ):
		print "%s...done" %progname
	
	E2end(logger)
	exit(1)
	
	
	#TODO - determine whether or not to add this functionality
	#if options.verifymirror:
		#if options.sym:
			#verify_mirror_test(data, eulers, options.sym, options.projector)
		#else:
			#print "Warning: verify mirror only works when a symmetry has been specified. No action taken."

# this actually doesn't work - the distribution is not uniform (FIXME)
def get_random_orientations( totalprojections, nomirror ):
		
	i = 0
	eulers = []
	while ( i < totalprojections ):
		#get a random 3D vector on [-1,1]
		x = 2*random.random() - 1
		y = 2*random.random() - 1
		if nomirror:
			# if no mirroring is specified, then only consider projections in positive z
			print "no mirroring"
			z = random.random()
		else:
			z = 2*random.random() - 1
		
		length = math.sqrt( x**2 + y**2 + z**2)
		
		#if the point is beyond the unit sphere then we should not 
		#consider it, or else the sampling of Euler angles will be
		#biased
		if length > 1:
			continue
		
		#normalize
		x = x/length
		y = y/length
		z = z/length
		
		# because the point is randomly distributed over the unit sphere,
		# its associated eulers will also be random.
		# note that a point on the unit sphere only implies two angles:
		az = math.atan2(y,x)
		alt = math.acos(z)
			
		randomphi = random.random()*360.0
		transform = Transform3D(EULER_EMAN,az*180.0/pi,alt*180.0/pi,randomphi)
		eulers.append(transform)
		
		i += 1
	
	return eulers

def generate_and_save_projections(options, data, eulers, smear=0):
	for i,euler in enumerate(eulers):
		d = {"t3d":euler}
		p=data.project(options.projector,d)
		p.set_rotation(euler)
		
		if smear:
			pass
			#smear_iterator = deg2rad*euler[2] + deg2rad*phiprop/(smear+1)
			#while ( smear_iterator < euler[2] + phiprop*deg2rad ):
				#ptmp=data.project(options.projector,{"alt" : euler[0],"az" : euler[1],"phi" : smear_iterator * rad2deg})
				
				#p.add(ptmp)
				#smear_iterator +=  deg2rad*phiprop/(smear+1)
		try: 
			p.write_image(options.outfile,-1)
		except:
			print "Error: Cannot write to file %s"%options.outfile
			exit(1)
		
		if (options.verbose):
			d = euler.get_rotation()
			print "%d\t%4.2f\t%4.2f\t%4.2f" % (i, d["az"], d["alt"], d["phi"])


#def verify_mirror_test(data, eulers, symmetry, projector):
	
	#sym_object = get_sym_object( symmetry )
	
	#for i,euler in enumerate(eulers) :
		#a = {"alt" : euler[0] * rad2deg,"az" : euler[1] * rad2deg,"phi" : euler[2] * rad2deg}
		#t3d = Transform3D(EULER_EMAN, a)
		#b = {"t3d": t3d }
		#p=data.project(options.projector,b)
		
		#mirrorEuler = sym_object.asym_unit_mirror_orientation(euler)
			
		
		##Get the projection in this orientations
		#p_mirror = data.project(projector,{"alt" : mirrorEuler[0] * rad2deg,"az" : mirrorEuler[1]* rad2deg,"phi" : mirrorEuler[2] * rad2deg})
		
		### Actually do the mirroring
		#if sym_object.is_d_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'y'})
		#elif sym_object.is_c_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'x'})
		##FIXME: The mirror orientation is dependent on the platonic symmetry see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry
		#elif sym_object.is_platonic_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'y'})
			
		
		### Calculate and write the difference to disk
		#p_difference = p_mirror-p
		#p_difference.write_image(symmetry+"mirror_debug_difference.img",-1)
		#p_mirror.write_image(symmetry+"mirror_debug.img",-1)
		
		### Print debug information
		#print "Orientation %d\t%4.2f\t%4.2f\t%4.2f" % (i, euler[0] * rad2deg, euler[1] * rad2deg, euler[2] * rad2deg)
		#print "Mirror %d\t%4.2f\t%4.2f\t%4.2f" % (i, mirrorEuler[0] * rad2deg, mirrorEuler[1] * rad2deg, mirrorEuler[2] * rad2deg)

def check(options, verbose=False):
	
	error = False
	
	if ( not options.sym ):
		if verbose:
			print "Error: you must specify the sym argument"
		error = True
	else:
		try: sym = parsesym(options.sym)
		except Exception, inst:
			if ( verbose ):
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args:
			error = True
	
	if ( not options.orientgen ):
		if verbose:
			print "Error: you must specify the orientgen argument"
		error = True
	elif ( check_eman2_type(options.orientgen,OrientGens,"Orientgen") == False ):
		error = True
	
	if ( check_eman2_type(options.projector,Projectors,"Projector") == False ):
		error = True

	if not os.path.exists(options.model):
		if verbose:
			print "Error: 3D image %s does not exist" %options.model
		error = True
	
	if ( options.force and options.append):
		if verbose:
			print "Error: cannot specify both append and force"
		error = True
		
	if ( options.nofilecheck == False and os.path.exists(options.outfile )):
		if ( not options.force ):
			if verbose:
				print "Error: output file exists, use -f to overwrite or -a to append. No action taken"
			error = True
	
	return error

if __name__=="__main__":
	main()
