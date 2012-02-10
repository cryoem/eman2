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
	usage = """%prog [options]

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
	
	
	parser.add_option("--nptcls", type="int",default=10,help="Number of simulated subtomograms tu generate per referece.")
	parser.add_option("--txrange", type="int",default=None,help="""Maximum number of pixels to randomly translate each subtomogram in X. The random translation will be picked between -txrange and +txrange. 
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_option("--tyrange", type="int",default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Y. The random translation will be picked between -tyrange and +tyrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_option("--tzrange", type="int",default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Z. The random translation will be picked between -tzrange and +tzrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_option("--transrange", type="int",default=4,help="""Maximum number of pixels to randomly translate each subtomogram in all X, Y and Z. 
									The random translation will be picked between -transrage and +transrange; --txrange, --tyrange and --tzrange overwrite --transrange for each specified direction.""")
	
	parser.add_option("--tiltstep", type="int",default=10,help="Degrees between each image in the simulated tilt series for each subtomogram.")
	parser.add_option("--tiltrange", type="int",default=10,help="""Maximum angular value at which the highest tilt picture will be simulated. Projections will be simulated from -tiltrange to +titlrange. 
									For example, if simulating a tilt series collected from -60 to 60 degrees, enter a --tiltrange value of 60. 
									Note that this parameter will determine the size of the missing wedge.""")
	parser.add_option("--nptcls", type="int",default=10,help="Number of simulated subtomograms tu generate per referece.")
	
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
	
	tag = ''
	
	for i in range(nrefs):
		
		randptcls = []
		if nrefs > 1:
			tag = str(i).zfill(len(str(nrefs)))
	
		model = EMData(options.input,i)
		if model['nx'] != model['ny'] or model['nx'] != model['nz'] or model['ny'] != model['nz']:
		
			newsize = max(model['nx'], model['ny'], model['nz'])
			if newsize % 2:
				newsize += 1
						
			os.system('e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + str(mewsize) + ' --first=' + str(i) + ' --last=' + str(i))	
	
		model = EMData(options.input,i)

		if options.filter != None:
			model.process_inplace(options.filter[0],options.filter[1])
		
		if options.shrink != None and options.shrink > 1 :
			model = model.process("math.meanshrink",{"n":options["shrink"]})
		else:
			model = model
		
		model.process_inplace('normalize')
			
		randptcls = randomizer(options, model, tag)
		
		subtomosim(options,randptcls)
					
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
RANDOMIZER - Takes a file in .hdf format and generates a set of particles randomly rotated and translated.
====================
'''
def randomizer(options, model, tag):
	
	stackname = options.input.replace('.hdf','_st' + tag + '_n' + str(options.n).zfill(len(str(options.n))) + '.hdf')

	print "I am inside the RANDOMIZER"
	
	if parameters['verbose']:
		print "You have requested to generate %d particles with random orientations and translations" %(options.n)
	
	randptcls = []
	for i in range(options.n):
		if parameters['verbose']:
			print "I will generate particle #", i
		
		b = model.copy()
		b['origin_x'] = 0							#Make sure the origin is set to zero, to avoid display issues with Chimera
		b['origin_y'] = 0
		b['origin_z'] = 0

		rand_orient = OrientGens.get("rand",{"n":1, "phitoo":1})		#Generate a random orientation (randomizes all 3 euler angles)
		c1_sym = Symmetries.get("c1")
		random_transform = rand_orient.gen_orientations(c1_sym)[0]

		randtx = random.randrange(-1 * options.transrange, options.transrange)	#Generate random translations
		randty = random.randrange(-1 * options.transrange, options.transrange)
		randtz = random.randrange(-1 * options.transrange, options.transrange)
		
		if options.txrange:
			randtx = random.randrange(-1 * options.txrange, options.txrange)	
		if options.tyrange:
			randty = random.randrange(-1 * options.tyrange, options.tyrange)
		if options.tzrange:
			randtz = random.randrange(-1 * options.tzrange, options.tzrange)
	
		random_transform.translate(randtx, randty, randtz)

		b.transform(random_transform)		
		
		b['spt_randT'] = random_transform
		b['xform.align3d'] = Transform()					#This parameter should be set to the identity transform since it can be used later to determine whether
											#alignment programs can "undo" the random rotation in spt_randT accurately or not
		#b.write_image(stackname,i)
		
		randptcls.append(b)
		
		if parameters['verbose']:
			print "The random transform applied to it was", rand 
	
	return(randptcls)
	
'''
====================
SUBTOMOSIM takes a set of particles in .hdf format and generates a simulated sub-tomogram for each, using user-defined parameters for tomographic simulation.
The function generates projections for each particle with a user-defined tilt step and missing wedge size (or data collection range),
adds noise and ctf to each projection, and recounstructs a new 3D volume from the simulated tilt series.
====================
'''	
def subtomosim(options,ptcls):
		
	lower_bound = -1 * options.tiltrange
	upper_bound = options.tiltrange
	
	nslices = int(round((upper_bound - lower_bound)/ options.tiltstep))
	
	if parameters['verbose']:
		print "There are these many particles in the set", len(ptcls)
		print "And these many slices to simulate each subtomogram", nslices
	
	for i in range(len(ptcls)):
		if parameters['verbose']:
			print "Projecting and adding noise to particle #", i

		apix = ptcls[i]['apix_x']
		
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
