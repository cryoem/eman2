#!/usr/bin/env python

#
# Author: Steve Ludtke (6/23/2013)  rewrote John Flanagan's program from scratch
# Copyright (c) 2000-2011 Baylor College of Medicine
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

import os
from EMAN2 import *
import math

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	
	WARNING: While functional, this program is not well tested, and likely needs development before being useful

	This program is designed to generate an initial reconstruction via monte carlo methodology.
	Basically, this just reads in class averages and assigns them random Euler angles. Next it computes how
	well the class avergaes agree with each other, similar to common lines"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="montecarloheader", help='Options below this label are specific to e2montecarlorecon', title="### e2montecarlorecon options ###", row=1, col=0, rowspan=1, colspan=3)
	parser.add_argument("--input",type=str,default=None,help="averages/particles images to process",browser="EMBrowserWidget(withmodal=True,multiselect=False)", guitype='filebox', row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--path",type=str,help="Name of path for output file",default='initial_models')
	parser.add_argument("--output",type=str,default="mc.mrc",help="Name of computed reconstruction, default=mc.mrc", guitype='strbox', row=2, col=0, rowspan=1, colspan=3)
	parser.add_argument("--shrink",type=int,default=0,help="Amount to shrink the CAs, default=0, no shrinking", guitype='shrinkbox', row=4, col=2, rowspan=1, colspan=1)
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry. For now it is not recommened to use this because in genmeral ab initio models are not aligned to symmtry axis. To apply symmetry use e2symsearch3d\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	global options
	(options, args) = parser.parse_args()

	#print help
	if options.classavg == None:
		parser.print_help()
		exit(0)

	# Initialize CUDA iof needed
	if options.cuda: initializeCUDAdevice()

	inimodeldir = os.path.join(".",options.path)
	if not os.access(inimodeldir, os.R_OK):
		os.mkdir(options.path)

	logid=E2init(sys.argv,options.ppid)

	calist = EMData.read_images(options.classavg)	# Load the CAs

	# shrink the CAs if needed
	if options.shrink:
		for ca in calist:
			ca.process_inplace('math.meanshrink',{'n':options.shrink})

	trials = len(calist)*options.mccoeff		# How many MC trials?
	Util.set_randnum_seed(Util.get_randnum_seed())	# Initialize random num generator

	# Setup the reconstructor
	if options.cuda: EMData.switchoncuda()
	size = (int(calist[0].get_attr('nx')),int(calist[0].get_attr('nx')),int(calist[0].get_attr('nx')))	# The images must all be the same size or we will get an exception
	reconstructor = ("fourier", {"size":size,"sym":options.sym,"verbose":options.verbose})
	freconstructor=Reconstructors.get(reconstructor[0], reconstructor[1])
	freconstructor.setup()

	print "\n\nTrying "+str(trials)+" Monte Carlo trials"
	besttlist = []
	bestscore = 0.0
	for mctrial in xrange(int(trials)):
		if options.verbose==1 :
			if((mctrial % 50) == 0):
				print "MC trial", mctrial
		elif options.verbose :
			print "MC trial %5d: "%mctrial,
			sys.stdout.flush()

		score = 0.0
		tlist = []
		for canum, ca in enumerate(calist):
			# Make random angles, but make sure that we don't cluster around the poles on the sphere
			# See http://mathworld.wolfram.com/SpherePointPicking.html for more details
			az = Util.get_frand(0,360) 					# theta
			alt  = math.degrees(math.acos(2*Util.get_frand(0,1) - 1))	# phi
			phi = Util.get_frand(0,360)					# kappa
			#phi = 0.0							# for testing only
			t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
			tlist.append([t,0])

			# Now insert into volume, at first we do not do much preprocessing
			freconstructor.insert_slice(ca,t,1)

		# Now compute the CA simmilarity score
		for canum, ca in enumerate(calist):
			freconstructor.determine_slice_agreement(ca,tlist[canum][0],1,1)
			tlist[canum][1] = ca.get_attr("reconstruct_absqual") # save the score for later use in refinement
			score += tlist[canum][1]


		freconstructor.clear()

		if(score > bestscore):
			bestscore = score
			besttlist = tlist

		if options.verbose>1 : print "%f (%f)"%(score,bestscore)

	refiner = SAsca("simulated annealing, single class average steps\n")
	#refiner = SA("simulated annealing\n")
	refiner.refinerecon(calist, besttlist, reconstructor)

	# Now do reconstruction
	bfreconstructor=Reconstructors.get(reconstructor[0], reconstructor[1]) # lets make a new reconstrcutor
	bfreconstructor.setup()
	for canum, ca in enumerate(calist):
		#ca.set_attr("xform.projection", besttlist[canum])
		bfreconstructor.insert_slice(ca, besttlist[canum][0],1)
	recon = bfreconstructor.finish(1)
	if options.cuda: EMData.switchoffcuda()

	#write output
	if options.output[-3:] == "mrc":
		recon.set_attr('UCSF.chimera',1)
	recon.write_image(os.path.join(inimodeldir,options.output))
	E2end(logid)

# Strategy pattern, allows other algoithms to be plugged in
class Refine:
	def __init__(self, name):
		self.name = name
	def refinerecon(self, calist, blist, reconstructor):
		raise NotImplementedError("Subclass must implement abstract method")
# Do per image SA
class SAsca(Refine):
	def refinerecon(self, calist, blist, reconstructor):
		print "Running: "+self.name

		#setup a reconstructor for use in refinement
		rreconstructor=Reconstructors.get(reconstructor[0], reconstructor[1]) 	# lets make a new reconstrcutor
		rreconstructor.setup()
		#for canum, ca in enumerate(calist):
			#rreconstructor.insert_slice(ca, blist[canum][0],1)

		temp = options.initemp
		K = options.numsasteps*options.numtemps
		for tstep in xrange(int(options.numtemps)):
			print "Temperature is", temp
			for i in xrange(int(options.numsasteps)):
				for canum, ca in enumerate(calist):
					rreconstructor.insert_slice(ca, blist[canum][0],1)

				for canum, ca in enumerate(calist):
					#step 1 remove the ca from the recon
					rreconstructor.insert_slice(ca, blist[canum][0],-1)

					#step 2 change the Eulers at random
					az = Util.get_frand(0,360) 					# theta (az)
					alt  = math.degrees(math.acos(2*Util.get_frand(0,1) - 1))	# phi (alt)
					phi = Util.get_frand(0,360)					# kappa (phi)
					t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})

					#step 3 determine its slice agreement
					rreconstructor.determine_slice_agreement(ca,t,1,0)
					e = ca.get_attr("reconstruct_absqual")
					de = e - blist[canum][1]

					#step 4 from the change in objetcive funtion keep the change if we go uphill, if we go downhill keep the changes some of the time
					if self.pseudoboltzmann(-de, temp):					 # we actually minimize so the negative is taken
						#print "going uphill from "+str(blist[canum][1])+" to "+str(e)
						rreconstructor.insert_slice(ca, t,1)
						blist[canum][0] = t
						blist[canum][1] = e
					else:
						rreconstructor.insert_slice(ca, blist[canum][0],1)
				rreconstructor.clear() # we need to do this to avoid intepolation errors building up
			temp=options.initemp*math.pow((1-tstep*options.numsasteps/K),options.cooling) # This annealing schedule comes from Numerical Recipes, second edition, pg 554
	# should we make the change?
	def pseudoboltzmann(self, de, temp):
		return de < 0 or Util.get_frand(0,1) < math.exp(-de/temp)

# This algorithm is a bit Rubbish!!!!
# DO all image SA
class SA(Refine):
	def refinerecon(self, calist, blist, reconstructor):
		print "Running: "+self.name

		#first get the initial energy
		inienergy = 0
		for i, data in enumerate(blist):
			inienergy += blist[i][1]

		#setup a reconstructor for use in refinement
		rreconstructor=Reconstructors.get(reconstructor[0], reconstructor[1]) 	# lets make a new reconstrcutor
		rreconstructor.setup()

		searchfract = 0.5
		temp = options.initemp
		for i in xrange(int(options.numsasteps)):
			energy = 0
			newt = []
			for canum, ca in enumerate(calist):
				#step 1 perturb the Eulers, and fill the recon
				daz = Util.get_gauss_rand(0,360*searchfract) 					# deltatheta
				dalt = Util.get_gauss_rand(0,1*searchfract)					# deltaphi
				dphi = Util.get_gauss_rand(0,360*searchfract)					# deltakappa
				az = blist[canum][0].get_rotation("eman")["az"] + daz				# theta (az)
				currentalt = blist[canum][0].get_rotation("eman")["alt"]			# get the current alt angle
				v = (math.cos(currentalt) + 1)/2 + dalt						# back convert to uniform variate and add the perturbation
				v =  math.acos(math.cos(v*math.pi))/math.pi					# this ensures that v E [-1, 1] and wraps around if v goes outside its range
				alt = math.degrees(math.acos(2*v - 1))						# Finally compute phi (alt)
				phi = blist[canum][0].get_rotation("eman")["phi"] + dphi			# kappa (phi)
				#az = Util.get_frand(0,360) 					# theta
				#alt  = math.degrees(math.acos(2*Util.get_frand(0,1) - 1))	# phi
				#phi = Util.get_frand(0,360)					# kappa
				t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
				newt.append([t,0])

				#step 2 fill the recon
				rreconstructor.insert_slice(ca, t,1)

			for canum, ca in enumerate(calist):
				#step 3 determine slice aggrement and sum
				rreconstructor.determine_slice_agreement(ca,newt[canum][0],1,1)
				newt[canum][1] = ca.get_attr("reconstruct_absqual")
				energy += newt[canum][1]


			rreconstructor.clear()
			de = energy - inienergy
			if self.pseudoboltzmann(-de, temp):
				print "Found a better match"
				blist[:] = newt[:] # deep copy going on here
				inienergy = energy
			#print inienergy, energy
	# should we make the change?
	def pseudoboltzmann(self, de, temp):
		return de < 0 or Util.get_frand(0,1) < math.exp(-de/temp)

if __name__ == "__main__":
	main()
