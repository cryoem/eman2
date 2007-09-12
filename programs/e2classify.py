#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/03/2007 (sludtke@bcm.edu)
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
	usage = """%prog [options] <simMatrixIn> <simMatrixOut>
	Takes a similarity matrix that has been created between reprojections-input (col) and particles-input (row) stacks of 2-D images.
	Typically these have been created via e2simmx such as
	 e2simmx.py proj.hed part.hed simMatrix.hed --saveali --align=rotate_translate:maxshift=5
	The similarity matrix is a stack of 4 images: the similarities is the first image and the translational shifts and rotations are the final 3 images.
	The output, simMatrixOut, is 5 images: with the 2nd thru 4th being identical to that of simMatrixIn (the shifts).
	The first image gives a list of weights as rows. Each weight represents the amount to which the given particle belongs. The projection to which it belongs is written in the 5th image.
	"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--sep", type="int", help="The number of classes a particle can contribute towards (default is 5)", default=5)
	
	(options, args) = parser.parse_args()
#	print(args);
	print(len(args));
	if len(args)<2 : parser.error("Input and output files required")
	
	
	E2n=E2init(sys.argv)

	Simimg=EMData();
	Simimg.read_image(args[0],0)
	TransxImg=EMData();
	TransxImg.read_image(args[0],1)
	TransyImg=EMData();
	TransyImg.read_image(args[0],2)
	RotImg=EMData();
	RotImg.read_image(args[0],3)
	
	
#	Simimg=getImage(args[0],0);
	NumProj= Simimg.get_xsize()
	NumPart= Simimg.get_ysize()
	
	# get the average score for each particle
	msimMat=range(NumPart)
	for iPart in range(NumPart):
		vv=[Simimg.get_value_at(iProj, iPart) for iProj in range(NumProj)];
		msimMat[iPart]=sum(vv)/NumProj;

	
	NumPeaks=min(options.sep,NumPart);
	
	simMatA       = [ range(NumPart) for j in range(NumProj)]
	simMatB       = [ range(NumPart) for j in range(NumProj)]
	ReturnVal     = [ range(NumPeaks) for j in range(NumPart)]
	ReturnPart    = [ range(NumPeaks) for j in range(NumPart)]
	ReturnWt      = [ range(NumPeaks) for j in range(NumPart)]
	ReturnTransx  = [ range(NumPeaks) for j in range(NumPart)]
	ReturnTransy  = [ range(NumPeaks) for j in range(NumPart)]
	ReturnRot     = [ range(NumPeaks) for j in range(NumPart)]
	
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			simMatA[iProj][iPart]= - ( Simimg.get_value_at(iProj,iPart) - msimMat[iPart]);

	#for iPart in range(NumPart):
		#vv=[simMatA[iProj][iPart] for iProj in range(NumProj)];
		#print sum(vv)/NumProj

	maxSim = -10000000
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			if ( simMatA[iProj][iPart] > maxSim ):
				maxSim = simMatA[iProj][iPart]
				

	#maxSim = max(max(simMatA));
	
	for iPart in range(NumPart):
		for iProj in range(NumProj):
			simMatB[iProj][iPart]=10* simMatA[iProj][iPart]/maxSim;
	
	
	
	for iPart in range(NumPart):
		vv=[ simMatB[j][iPart] for j in range(NumProj)];
		vvCp=[ vv[j] for j in range(NumProj)];
		vvCp.sort();
		vvCp.reverse();
		vvIndices=[ vv.index(vvCp[j]) for j in range(NumPeaks)];
		for iPeaks in range(NumPeaks):
			RelProj             =     vvIndices[iPeaks];
			ReturnPart[iPart][iPeaks] = RelProj;
			ReturnVal[iPart][iPeaks]  = vvCp[iPeaks];
			ReturnTransx[iPart][iPeaks]= TransxImg.get_value_at(RelProj,iPart);
			ReturnTransy[iPart][iPeaks]= TransyImg.get_value_at(RelProj,iPart);
			ReturnRot[iPart][iPeaks]= RotImg.get_value_at(RelProj,iPart);
	
	
	
	for iPart in range(NumPart):
		vv=[ ReturnVal[iPart][j] for j in range(NumPeaks)];
		sumvv= sum(vv)
		for iPeaks in range(NumPeaks):
			ReturnWt[iPart][iPeaks]= ReturnVal[iPart][iPeaks]/sumvv;
	
	
	ReturnPartImg  =EMData();
	ReturnPartImg.set_size(NumPeaks,NumPart) ; ReturnPartImg.to_zero();
	
	ReturnWtImg    =EMData();
	ReturnWtImg.set_size(NumPeaks,NumPart)   ;  ReturnWtImg.to_zero();

	ReturnTransxImg=EMData();
	ReturnTransxImg.set_size(NumPeaks,NumPart); ReturnTransxImg.to_zero();

	ReturnTransyImg=EMData();
	ReturnTransyImg.set_size(NumPeaks,NumPart); ReturnTransyImg.to_zero();
	
	ReturnRotImg   =EMData();
	ReturnRotImg.set_size(NumPeaks,NumPart)   ; ReturnRotImg.to_zero();
	
	for iPart in range(NumPart):
		for iPeaks in range(NumPeaks):
			ReturnWtImg.set_value_at(iPeaks, iPart,ReturnWt[iPart][iPeaks] );
			ReturnPartImg.set_value_at(iPeaks,iPart, ReturnPart[iPart][iPeaks] );
			ReturnTransxImg.set_value_at(iPeaks,iPart, ReturnTransx[iPart][iPeaks] );
			ReturnTransyImg.set_value_at(iPeaks,iPart, ReturnTransy[iPart][iPeaks] );
			ReturnRotImg.set_value_at(iPeaks,iPart, ReturnRot[iPart][iPeaks] );
	
	
	ReturnPartImg.write_image(args[1],0)
	ReturnWtImg.write_image(args[1],1);
	ReturnTransxImg.write_image(args[1],2);
	ReturnTransyImg.write_image(args[1],3);
	ReturnRotImg.write_image(args[1],4);

	

if __name__ == "__main__":
    main()
