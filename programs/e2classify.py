#!/usr/bin/env python

#
# Author: David Woolford, 9/12/2007 (woolford@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
from optparse import OptionParser
from math import *
import os
import sys



def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <simMatrixIn> <simMatrixOut> [options]
	Takes a similarity matrix that has been created between reprojections-input (col) and particles-input (row) stacks of 2-D images.
	Typically these have been created via e2simmx such as
	 e2simmx.py proj.hed part.hed simMatrix.hed --saveali --align=rotate_translate:maxshift=5
	The similarity matrix is a stack of 4 images: the similarities is the first image and the translational shifts and rotations are the final 3 images.
	The output, simMatrixOut, is 5 images: with the 2nd thru 4th being identical to that of simMatrixIn (the shifts).
	The first image gives a list of weights as rows. Each weight represents the amount to which the given particle belongs. The projection to which it belongs is written in the 5th image.
	"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--sep", type="int", help="The number of classes a particle can contribute towards (default is 5)", default=5)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_option("--verbose", "-v",dest="verbose",default=False, action="store_true",help="Toggle verbose mode - prints extra infromation to the command line while executing")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	
	#      Setion 0.
#      read in the 4 images from <simMatrixIn> which correspond to 
#        0. similarity 
#        1. Transx 
#        2. Transy 
#        3. Rotation	
	
	(options, args) = parser.parse_args()
	
	if options.nofilecheck: options.check = True

	if (len(args)<2 ): parser.error("Input and output files required")
	
	if (options.check): 
		options.verbose = True # turn verbose on if the user is only checking...
		
	if ( options.nofilecheck == False ):
		options.simmxfile = args[0]
		options.outfile = args[1]
		
	error = check(options,True)
	
	if (options.verbose):
		if (error):
			print "e2classify.py command line arguments test.... FAILED"
		else:
			print "e2classify.py command line arguments test.... PASSED"
		
	if error : exit(1)
	if options.check: exit(0)

	
	E2n=E2init(sys.argv)

	if os.path.exists(args[1]):
		if (options.force):
			remove_file(args[1])

	num_sim =  EMUtil.get_image_count(args[0])
	if (num_sim != 5):
		print "Error, the similarity matrix did not contain  5 images - be sure to use the --saveali argument when running e2simmx.py"
		exit(1)

	if ( options.verbose ):
		print "Reading image data..."
	Simimg=EMData();
	Simimg.read_image(args[0],0)
	TransxImg=EMData();
	TransxImg.read_image(args[0],1)
	TransyImg=EMData();
	TransyImg.read_image(args[0],2)
	RotImg=EMData();
	RotImg.read_image(args[0],3)
	MirrorImg=EMData();
	MirrorImg.read_image(args[0],4)
	
	NumProj= Simimg.get_xsize()
	NumPart= Simimg.get_ysize()
	
#	Section 1.
#     get the average alignment for each particle, and the maximum; renormalize
#     the values of the similarities to get simMatB
#           create lists that will eventually be used to write out the images
	if ( options.verbose ):
		print "Getting similarity statistics..."
	msimMat=range(NumPart)
	for iPart in range(NumPart):
		vv=[Simimg.get_value_at(iProj, iPart) for iProj in range(NumProj)];
		msimMat[iPart]=sum(vv)/NumProj;

	
	NumPeaks=min(options.sep,NumPart);
	
	simMatA       = [ range(NumPart) for j in range(NumProj)]
	simMatB       = [ range(NumPart) for j in range(NumProj)]
	ReturnPart    = [ range(NumPart) for j in range(NumPeaks)]
	ReturnVal     = [ range(NumPart) for j in range(NumPeaks)]
	ReturnWt      = [ range(NumPart) for j in range(NumPeaks)]
	ReturnTransx  = [ range(NumPart) for j in range(NumPeaks)]
	ReturnTransy  = [ range(NumPart) for j in range(NumPeaks)]
	ReturnRot     = [ range(NumPart) for j in range(NumPeaks)]
	ReturnMirror  = [ range(NumPart) for j in range(NumPeaks)]
	
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			simMatA[iProj][iPart]= - ( Simimg.get_value_at(iProj,iPart) - msimMat[iPart]);


	maxSim = -10000000
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			if ( simMatA[iProj][iPart] > maxSim ):
				maxSim = simMatA[iProj][iPart]
				
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			simMatB[iProj][iPart]=10* simMatA[iProj][iPart]/maxSim;
	
	
#	Section 2.
#       i) sort the alignment data [vvCp], 
#      ii) find the projections corresponding to the sorted list [RelProj]
#      iii) get corresponding  transx, transy, rot
#      iv) rescale the weights to sum to unity (over the NumPeaks values)
	
	if ( options.verbose ):
		print "Determining particle classes and contribution weights..."
	
	for iPart in range(NumPart):
		vv=[ simMatB[j][iPart] for j in range(NumProj)];
		vvCp=[ vv[j] for j in range(NumProj)];
		vvCp.sort();
		vvCp.reverse();
		vvIndices=[ vv.index(vvCp[j]) for j in range(NumPeaks)];
		for iPeaks in range(NumPeaks):
			RelProj             =     vvIndices[iPeaks];
			ReturnPart[iPeaks][iPart] = RelProj;
			ReturnVal[iPeaks][iPart]  = vvCp[iPeaks];
			ReturnTransx[iPeaks][iPart]= TransxImg.get_value_at(RelProj,iPart);
			ReturnTransy[iPeaks][iPart]= TransyImg.get_value_at(RelProj,iPart);
			ReturnRot[iPeaks][iPart]= RotImg.get_value_at(RelProj,iPart);
			ReturnMirror[iPeaks][iPart]= MirrorImg.get_value_at(RelProj,iPart);
	
	
	
	for iPart in range(NumPart):
		vv=[ ReturnVal[j][iPart] for j in range(NumPeaks)];
		sumvv= sum(vv);
		if ( sumvv != 0 ):
			for iPeaks in range(NumPeaks):
				ReturnWt[iPeaks][iPart]= ReturnVal[iPeaks][iPart]/sumvv;	
#	Section 3.
#       write to outFile	
	if ( options.verbose ):
		print "Writing output data..."
	
	writeListToImage(ReturnPart,   args[1]   ,0,NumPeaks,NumPart)
	writeListToImage(ReturnWt,     args[1]   ,1,NumPeaks,NumPart)
	writeListToImage(ReturnTransx, args[1]   ,2,NumPeaks,NumPart)
	writeListToImage(ReturnTransy, args[1]   ,3,NumPeaks,NumPart)
	writeListToImage(ReturnRot,    args[1]   ,4,NumPeaks,NumPart)
	writeListToImage(ReturnMirror, args[1]   ,5,NumPeaks,NumPart)
	
	E2end(E2n)

	if ( options.verbose ):
		print "e2classify...done"

def writeListToImage(inList,outFileString,location,NumPeaks,NumPart):
	"""Compares one image (target) to a list of many images (reflist). Returns """
	
	OutImg   =EMData();
	OutImg.set_size(NumPeaks,NumPart)   ; 
	OutImg.to_zero();
	for iPart in range(NumPart):
		for iPeaks in range(NumPeaks):
			OutImg.set_value_at(iPeaks,iPart, inList[iPeaks][iPart] );
			
	OutImg.write_image(outFileString,location)
	return 
	

def check(options,verbose):
	error = False
	
	if (options.sep < 1):
		if verbose:
			print "Error: the --sep argument must be greater than zero, currently it is %d" %(options.sep)
		error = True
	

	
	if ( options.nofilecheck == False ):
		if os.path.exists(options.outfile):
			if (not options.force):
				if verbose:
					print "File %s exists, will not write over, exiting" %options.outfile
				error = True
		
		if not os.path.exists(options.simmxfile) and not db_check_dict(options.simmxfile):
			if verbose:
				print "Error: the similarity matrix file (%s) was not found, cannot run e2classify.py" %(options.simmxfile)
			error = True
		else:
			num_sim =  EMUtil.get_image_count(options.simmxfile)
			if (num_sim != 5):
				if verbose:
					print "Error, the similarity matrix did not contain 5 images - be sure to use the --saveali argument when running e2simmx.py"
				error = True
	
	return error

if __name__ == "__main__":
    main()
