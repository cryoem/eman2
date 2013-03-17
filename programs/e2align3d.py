#!/usr/bin/env python

#
# Author: John Flanagan, 10/08/2010 (jfflanag@bcm.edu)
# Copyright (c) 2010 Baylor College of Medicine
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

from EMAN2 import *
import sys, os
	
def main():
  
	progname = os.path.basename(sys.argv[0])
	usage = """prog fixed_model moving_model output_model [options] 
	
	This program is designed to rotationally and translationally align two 3D models 
	Usually the two models are shrunk down to speed things up, then a global exhaustive 
	search is down by the refine.3d.sphere and then this rough alignment is refined using
	refine.3d using the full size maps. The refiner is much quicker than the global aligner
	as it uses a simplex algoritm to bring the alignment downhill, but if your inital global
	alignment is too rough then the refiner might get stuck in a local minima."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2align.py
	parser.add_argument("--shrink",type=int,default=1,help="Fractional amount to shrink the maps by, default=1")
	parser.add_argument("--preprocess",metavar="processor_name(param1=value1:param2=value2)",type=str,default=None,action="append",help="preprocess maps before alignment")
	parser.add_argument("--maskrad",type=int,default=-1,help="Mask the recon using a spherical Gaussian mask (-1 = None), default=-1.0")
	parser.add_argument("--maskfoff",type=float,default=0.1,help="Fall offf of the Gaussian mask, default=0.1")
	parser.add_argument("--nsolns",type=int,default=1,help="number of peaks in the global search to refine, default=1.0")
	parser.add_argument("--famps",type=float,default=1,help="fraction of Fourier amps to exclude from recons. 0 means that this option is not used, default=0.0")
	parser.add_argument("--prec",type=float,default=0.01,help="Precison to determine what solutions are the 'same' used only statistics output, default=0.01")
	parser.add_argument("--cuda",action="store_true", help="Use CUDA for the alignment step.",default=False)
	#options form the sphere alinger
	parser.add_argument("--delta",type=float,default=30.0,help="step size for the orrientation generator, default=30.0")
	parser.add_argument("--dphi",type=float,default=30.0,help="step size for the inplane angle phi, default=30.0")
	parser.add_argument("--phi0",type=float,default=0.0,help="lower bound for the inplane angle phi, default=0.0")
	parser.add_argument("--phi1",type=float,default=359.0,help="Upper bound for the inplane angle phi, default=359.0")
	parser.add_argument("--search",type=int,default=10,help="maximum extent of the translational search, default=10")
	parser.add_argument("--sym",type=str,default='c1',help="model symmetry (using sym, if present, speeds thing up a lot), default='c1'")
	parser.add_argument("--cmp",type=str,default='ccc',help="comparitor and params to use for the 3D refiner, default='ccc'")
	parser.add_argument("--dotrans",type=int,default=1,help="Do translational search, default=1")
	#options associated with  the simplex 3D refiner
	parser.add_argument("--ralign",type=str,default='refine_3d:spin_coeff=1',help="aligner to use for refine alignement, default='refine_3d:spin_coeff=1'")
	parser.add_argument("--rcmp",type=str,default='ccc',help="comparitor and params to use for the 3D refiner, default='ccc'")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose","-v",type=int,default=0,help="Level of verboseness, default=0")

	global options
	(options, args) = parser.parse_args()
	if options.cmp : options.cmp=parsemodopt(options.cmp)
	if options.ralign : options.ralign=parsemodopt(options.ralign)
	if options.rcmp : options.rcmp=parsemodopt(options.rcmp)
	if options.cuda: initializeCUDAdevice()
               
	logid=E2init(sys.argv,options.ppid)
	#read in the input maps and check for sanity
	fixed = EMData()
	moving = EMData()
        
	try:
		fixed.read_image(args[0])
	except:
		print "Not able to read file %s" % args[0]
		exit(1)
	    
	try:
		moving.read_image(args[1])
	except:
		print "Not able to read file %s" % args[1]
		exit(1)
	    
	if (fixed.get_attr('nx') != fixed.get_attr('ny') != fixed.get_a+ options.precttr('nz')):
		print "Fixed map must have cubic dimensions!"
		exit(1)
	    
	if (moving.get_attr('nx') != moving.get_attr('ny') != moving.get_attr('nz')):
		print "Fixed map must have cubic dimensions!"
		exit(1)		 
	    
	if (moving.get_attr('nx') != fixed.get_attr('nx')):
		print "Fixed and model maps must have the same dimensions!"
		exit(1)
			
	#preprocess maps
	if options.preprocess != None:
		for p in options.preprocess:
			try:
				(processorname, param_dict) = parsemodopt(p)
				if not param_dict : param_dict={}
				fixed.process_inplace(str(processorname), param_dict)
				moving.process_inplace(str(processorname), param_dict)
			except:
				print "warning - application of the pre processor",p," failed. Continuing anyway"

	#denoise recons
	if options.famps > 0:
		tmp = fixed.do_fft()
		tmp.process_inplace('threshold.binary.fourier',{'value':tmp.get_amplitude_thres(options.famps)})
		fixed = tmp.do_ift()
		tmp = moving.do_fft()
		tmp.process_inplace('threshold.binary.fourier',{'value':tmp.get_amplitude_thres(options.famps)})
		moving = tmp.do_ift()

	#mask out all the junk
	if options.maskrad > 0:
		fixed.process_inplace('mask.gaussian.nonuniform', {'radius_x':options.maskrad,'radius_y':options.maskrad,'radius_z':options.maskrad,'gauss_width':options.maskfoff})
		moving.process_inplace('mask.gaussian.nonuniform', {'radius_x':options.maskrad,'radius_y':options.maskrad,'radius_z':options.maskrad,'gauss_width':options.maskfoff})

	# shrink the maps, if desired
	if options.shrink > 1:
		sfixed = fixed.process('math.medianshrink', {'n':options.shrink})
		smoving = moving.process('math.medianshrink', {'n':options.shrink})
		options.search = options.search/options.shrink		# must adjust the search range
	else:
		sfixed = fixed
		smoving = moving
        
	# A hack to get around a Chimera bug
	fixed.set_attr('UCSF.chimera',1)
	fixed.write_image('filtered_fixed.mrc')

	# IF the dot product is being used then we need to normalize, othrwise the Quaternion aligner will crash
	if options.cmp[0] == "dot":
		cmpdict = options.cmp[1]
		cmpdict['normalize'] = 1
		options.cmp = (options.cmp[0], cmpdict)
	
	# do the global search
	bestscore = 0
	bestmodel = 0
	galignedref = []
	if options.cuda: EMData.switchoncuda()
	nbest = smoving.xform_align_nbest('rotate_translate_3d', sfixed, {'delta':options.delta,'dotrans':options.dotrans,'dphi':options.dphi,'search':options.search, 'phi0':options.phi0, 'phi1':options.phi1, 'sym':options.sym, 'verbose':options.verbose}, options.nsolns, options.cmp[0],options.cmp[1])
	
	# refine each solution found are write out the best one
	for i, n in enumerate(nbest):
		ralingdict = options.ralign[1]
		ralingdict['xform.align3d'] = n['xform.align3d']
		options.ralign = (options.ralign[0], ralingdict)
		print options.ralign[0], options.ralign[1]
		print options.rcmp[0], options.rcmp[1]
		galignedref.append(smoving.align(options.ralign[0], sfixed, options.ralign[1], options.rcmp[0], options.rcmp[1]))
		score = galignedref[i].get_attr('score')
		if score < bestscore:
			#bestscore = score
			bestmodel = i
		if options.verbose > 0: print "Peak Num: ", i, " Transform: ", n["xform.align3d"], " Ini Score: ", n["score"], " Final Score: ", score
	if options.cuda: EMData.switchoffcuda()
	
	# Find out how many peaks are 'the same' and print stats to the screen
	if options.nsolns > 1:
		thesame = 0
		for m in galignedref:
			if (m.get_attr('score') < bestscore + options.prec and m.get_attr('score') > bestscore - options.prec):
				thesame += 1
		print str(thesame)+" solns refined to the 'same' point within a precision of "+str(options.prec)
        
	#apply the transform to the original model
	moving.read_image(args[1])
	ft = galignedref[bestmodel].get_attr("xform.align3d")
	tft = ft.get_trans()
	ft.set_trans([tft[0]*options.shrink,tft[1]*options.shrink,tft[2]*options.shrink]) # Rescale the translation pars, otherwise it will be crap!
	moving.process_inplace("xform",{"transform":ft})
        
	#now write out the aligned model
	outfile=args[2]
	    
	#if output is mrc format or bdb 
	if outfile[-4:].lower() == ".mrc":
		moving.set_attr('UCSF.chimera',1)
		galignedref[bestmodel].set_attr('UCSF.chimera',1)
	if outfile[:4] == "bdb:":
		filtoutfile = outfile+"_filtered"
	else:
		filtoutfile = "filtered"+outfile
  
	galignedref[bestmodel].write_image(filtoutfile, 0)
	moving.write_image(outfile, 0)
	E2end(logid)
	
if __name__ == "__main__":
	main()
