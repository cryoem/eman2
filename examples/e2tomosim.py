#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - 2011, Last update: 02/2011
====================

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
'''

from optparse import OptionParser
from EMAN2 import *
from sys import argv
import EMAN2
import heapq
import operator
import random
import numpy

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <output> [options]

	This program produces simulated sub volumes in random orientations from a given PDB or EM file (mrc or hdf).
	"""
			
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--input", type="string", help="""The name of the input volume from which simulated subtomograms will be generated. 
							The output will be in HDF format, since volume stack support is required. The input CAN be PDB, MRC or and HDF stack. 
							If the input file is PDB or MRC, a version of the supplied model will be written out in HDF format.
							If the input file is a stack, simulatd subvolumes will be generated from each model in the stack and written to different output stacks.
							For example, if the input file contains models A and B, two output stacks with simulated subvolumes will be generated.""", default=None)
	
	#parser.add_option("--output", type="string", help="The name of the simulated stack. MUST be since volume stack support is required.", default=None)
			
	parser.add_option("--filter",type="string",help="A processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
	
	parser.add_option("--shrink", type="int",default=1,help="Optionally shrink the input volume before the simulation if you want binned/down-sampled subtomograms.")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n",type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()	
	
	if options.filter:
		options.filter = parsemodopt(options.filter)
	
	if options.shrink:
		options.shrink= parsemodopt(options.shrink)
	
	if '.pdb' in options.input or '.mrc' in options.input:
		wrongmodel = options.input
		hdfmodel = wrongmodel.replace('.pdb','.hdf')
		hdfmodel = wrongmodel.replace('.mrc','.hdf')
		os.system('e2pdb2mrc.py ' + pdbmodel + ' ' + hdfmodel)
		options.input = hdfmodel
	
	nrefs = EMUtil.get_image_count(options.input)
	
	for i in range(nrefs):
	
		model = EMData(options.input,i)
		if model['nx'] != model['ny'] or model['nx'] != model['nz'] or model['ny'] != model['nz']:
			os.system('e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + max(model['nx'], model['ny'], model['nz'])	
	
		if options.filter != None:
			model.process_inplace(options.filter[0],options.filter[1])
		
		if options.shrink != None and options.shrink > 1 :
			model = model.process("math.meanshrink",{"n":options["shrink"]})
		else:
			model = model
			
		tomosim(model)
					
	return()
	
	
	

'''
ANGDIST - Calculates the distance between two transforms. Used to verify the answers provided by the aligner for simmulated data,
and to determine the proximity of the 10 best peaks proposed, for experimental and simmulated data.
'''
def angdist(t1,t2):
	t2inv = t2.inverse()
	product = t2inv * t1
	product_SPIN = product.get_rotation('spin')
	angular_distance = round(float(product_SPIN["Omega"]),2)
	print "The angular distance is", angular_distance
	return(angular_distance)






'''
====================
RANDOMIZER - Takes a model (.hdf, .mrc or .pdb) makes it suitable for SPT (even, multiple of 8 box size) 
and produces a user defined n-number of particles from the model which are randomly rotated and translated within user-defined ranges
====================
'''
def randomizer(parameters,model = ''):
	print "I am inside the RANDOMIZER"
	
	random_stack_name = 'rand_orient_' + parameters['ID'] + '.hdf'
	
	reference= ''
	
	n = parameters['rand_num']
	trans_range = parameters['trans_range']
	rot_range = parameters['rot_range']
	rot_step = parameters['rot_step']
		
	'''
	The model can be supplied directly (an EMData file) if the function is imported. Otherwise, it's fetched from --model=
	If it's PDB or MRC it needs to be converted to HDF and resized if the box isn't cubical and even
	'''
	if model == '':
		print "NO MODEL WAS PASSED, so I have to load it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", parameters['model']
		if '.pdb' or '.mrc' in parameters['model']:
			print "IT contains EITHER pdb or mrc"
			
			mrcname = parameters['model'].replace('.pdb','.mrc')
			hdfname = mrcname.replace('.mrc','.hdf')

			if '.pdb' in parameters['model']:
				print "It contained PDB!!!!!!!!"
				os.system('e2pdb2mrc.py ' + parameters['model'] + ' ' + mrcname + ' && e2proc3d.py ' + mrcname + ' ' + hdfname)
			if '.mrc' in parameters['model']:
				print "It contained MRC!!!!!!!!"
				os.system('e2proc3d.py ' + mrcname + ' ' + hdfname)
			
			
			print "The hdf model to load is", hdfname
					
			model = EMData(hdfname,0)

		else:
			model = EMData(parameters['model'],0)

	if model['nx'] != model['ny'] or model['nx'] != model['nz'] or model['ny'] != model['nz']:
		model = cubebox(model)
	
	if parameters['shrink']:
		factor = 1.0/parameters['shrink']
		clip = model['nx']/parameters['shrink']
		model = model.process('xform.scale',{'scale':factor,'clip':clip})
	
	if model['nx'] % 2 or parameters['pad']:
		model = fixbox(parameters,model)	
		
	if parameters['alignment_type'] == 'refbased':		#Write out the "reference" for refbased alignments, in the model's default position
		reference = model.copy()
		reference.write_image('reference.hdf')	
	
	random_particles = []
	initlines = []			#Array to store the initial conditions (rotations and translations) as "lines" in a specific format to be written to a file
	
	if parameters['verbose']:
		print "You have requested to generate %d particles with random orientations and translations" %(n)
	for i in range(n):
		if parameters['verbose']:
			print "I will generate particle number", i
		
		b = model.copy()
		b['origin_x'] = 0
		b['origin_y'] = 0
		b['origin_z'] = 0
		print "The sizes of the copied model are", b['nx'], b['ny'], b['nz']
		print "And its type is", type(b)
		
		az=alt=phi=x=y=z=0
		
		initial_trans = Vec3f(0,0,0)
		if trans_range != 0:
			x = random.randrange(trans_range)
			y = random.randrange(trans_range)
			z = random.randrange(trans_range)
			b['rand_x'] = x
			b['rand_y'] = y
			b['rand_z'] = z
			print "This is how much b was translated", x, y, z
			initial_trans = Vec3f(x,y,z)
			b.translate(initial_trans)
		
		initial_rot = (0,0,0)
		
		if parameters['rot_range'] and parameters['rot_step']:
			factor = parameters['rot_range']/parameters['rot_step']
			az = random.randrange(factor) * parameters['rot_step']
			alt = random.randrange(factor / 2) * parameters['rot_step'] 
			phi = random.randrange(factor) * parameters['rot_step']
			b['rand_az'] = az
			b['rand_alt'] = alt
			b['rand_phi'] = phi
			b.rotate(az,alt,phi)
			print "This is how much b was rotated", az, alt, phi

		b['spt_initial'] = [az,alt,phi,x,y,z]
	
		b.write_image(random_stack_name,i)
		random_particles.append(b)		
		
		line = "particle #%s az=%s alt=%s phi=%s tx=%s ty=%s tz=%s\n" % ( str(i).zfill(len(str(n+1))), str(az).zfill(4),str(alt).zfill(4),str(phi).zfill(4),str(x).zfill(4),str(y).zfill(4),str(z).zfill(4) )
		initlines.append(line)

		if parameters['verbose']:
			print "These are the initial conditions for this particle", line

	init_file = open('initial_conditions_' + parameters['ID'] + '.txt','w')
	init_file.writelines(initlines)
	init_file.close()
	
	return(random_particles,reference)

'''
NEEDS TESTING

TOMO SIM - Takes a set of particles and turns them into simulated sub-tomograms.
[That is, it generates projections for each particle within a user-defined range in altitude, at a user-defined tilt-step,
adds noise and ctf to each projection, and recounstructs a new 3D volume from them
'''	
def tomo_sim(parameters,particles):
		
	lower_bound = parameters['tilt_range'] * -1
	upper_bound = parameters['tilt_range']
	
	slices = int(round((upper_bound - lower_bound)/int(parameters['tilt_step'])))
	if parameters['verbose']:
		print "There are these many particles", len(particles)
		print "And these many slices", slices
	
	for i in range(len(particles)):
		if parameters['verbose']:
			print "Projecting and adding noise to particle number", i

		apix = particles[i]['apix_x']
		
		px = random.uniform(-1*parameters['hole_diameter']/2, parameters['hole_diameter']/2)		#random distance of the particle center from the tilt axis
			
		alt = lower_bound
		raw_projections = []
		ctfed_projections = []
		for j in range(slices):
			t = Transform({'type':'eman','az':0,'alt':alt,'phi':0})
			
			dz = -1 * px * numpy.sin(alt)			#For particles left of the tilt axis, px is negative. With negative alt [left end of the ice down, right end up], dz should be negative.
			defocus = parameters['defocus'] + dz
						
			prj = particles[i].project("standard",t)
			prj.set_attr('xform.projection',t)
			prj.write_image('prj_' + parameters['ID'] + '.hdf',j)		#Write projections stack
			raw_projections.append(prj)
			
			ctf = EMAN2Ctf()
			ctf.from_dict({'defocus':defocus,'bfactor':100,'ampcont':0.05,'apix':apix,'voltage':parameters['voltage'],'cs':parameters['cs']})
			
			prj_fft = prj.do_fft()
			prj_ctf = prj_fft.copy()
			
			ctf.compute_2d_complex(prj_ctf,Ctf.CtfType.CTF_AMP)
			prj_ctf.mult(-1)			#Reverse the contrast, as in "authentic" cryoEM data		
			
			prj_fft.mult(prj_ctf)
			
			prj_r = prj_fft.do_ift()	#Go back to real space
			
			prj_n = prj.process('math.addnoise',{'noise':10})
			prj_n.write_image('prj_ctfed' + parameters['ID'] + '.hdf',j)	
			
			ctfed_projections.append(prj_n)		
						
			alt += int(parameters['tilt_step'])
		
		box = particles[i].get_xsize()
		r = Reconstructors.get(parameters['reconstructor'],{'size':(box,box,box),'sym':'c1','mode':'gauss_2','verbose':True})
		r.setup()
		
		for p in ctfed_projections:
			p = r.preprocess_slice(p,p['xform.projection'])
			r.insert_slice(p,p['xform.projection'],1.0)
		
		rec = r.finish(True)
		mname = parameters['model'].split('/')[-1].split('.')[0]
		name = 'rec_' + mname + '#' + str(i).zfill(len(str(len(particles)))) + '.hdf'
		rec.write_image(name)
	return()	











if __name__ == '__main__':
	main()
