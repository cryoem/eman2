#!/usr/bin/env python

#
# Author: Matthew Baker, 10/2005, modified 02/2006 by MFS  
# ported to EMAN2 by David Woolford October 6th 2008
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


#N tomohunter.py
#F tomography hunter

import os
import sys
import string
import commands
import math
from EMAN2 import *
#import Numeric
from math import *
from sys import argv
from optparse import OptionParser

def check_options(options,filenames):
	'''
	Check the parser options
	Should probably be made into a class
	@return a list of error messages
	'''
	error_messages = []
	if len(filenames) < 2:
		error_messages.append("You must specify input and output files")
	else:
		all_images = True
		for f in filenames:
			if not file_exists(f): 
				error_messages.append("Error - %s does not exist" %f)
				all_images = False
				
		if all_images:
			nx,ny,nz = gimme_image_dimensions3D(filenames[0])
			for i in range(1,len(filenames)):
				x,y,z = gimme_image_dimensions3D(filenames[i])
				if x != nx or y != ny or z != nz:
					error_messages.append("File %s does not have the same dimensions as file %s" %(filenames[i],filenames[0]))
									
	if options.nsoln <= 0:
		error_messages.append("Error - nsoln must be greater than 0. Suggest using 10")
	
	
	atleast_0 = ["ralt","rphi","raz","searchx","searchy","searchz"]
	for attr in atleast_0:
		if getattr(options,attr) < 0:
			error_messages.append("Error - %s must be greater than  or equal to 0" %attr)
			
	greater_than_0 = ["dalt","dphi","daz"]
	for attr in greater_than_0:
		if getattr(options,attr) <= 0:
			error_messages.append("Error - %s must be greater than 0" %attr)
	
	return error_messages

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <target image> <image to be aligned> [options]"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--dalt",type="float",help="Altitude delta", default=10.0)
	parser.add_option("--ralt",type="float",help="Altitude range", default=180.0)
	parser.add_option("--dphi",type="float",help="Phi delta", default=10.0)
	parser.add_option("--rphi",type="float",help="Phi range", default=180.0)
	parser.add_option("--raz",type="float",help="Phi range", default=360.0)
	parser.add_option("--daz",type="float",help="Azimuth delta", default=10.0)
	parser.add_option("--thresh",type="float",help="Threshold", default=1.0)
	parser.add_option("--nsoln",type="int",help="The number of solutions to report", default=10)
	parser.add_option("--searchx",type="int",help="The maximum search distance, x direction", default=5)
	parser.add_option("--searchy",type="int",help="The maximum search distance, y direction", default=5)
	parser.add_option("--searchz",type="int",help="The maximum search distance, z direction", default=5)
   
	(options, args) = parser.parse_args()

	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(0)
	
	targetMRC =EMData(argv[1])
	print "Target Information"
	print "   mean:	   %f"%(targetMRC.get_attr("mean"))
	print "   sigma:	  %f"%(targetMRC.get_attr("sigma"))
	
	target_xsize=targetMRC.get_xsize()
	target_ysize=targetMRC.get_ysize()
	target_zsize=targetMRC.get_zsize()
	if (target_xsize!=target_ysize!=target_zsize) or (target_xsize%2==1):
		print "The density map must be even and cubic. Terminating."
		sys.exit()
	
	box=target_xsize
	
	probeMRC=EMData(argv[2])
	print "Probe Information"
	print "   mean:	   %f"%(probeMRC.get_attr("mean"))
	print "   sigma:	  %f"%(probeMRC.get_attr("sigma"))
	
	
	bestCCF= EMData(options.nsoln,1,1)
	bestCCF.to_zero()

	bestAZ=EMData(options.nsoln,1,1)
	bestAZ.to_zero()

	bestALT= EMData(options.nsoln,1,1)
	bestALT.to_zero()

	bestPHI=EMData(options.nsoln,1,1)
	bestPHI.to_zero()
	
	bestX=EMData(options.nsoln,1,1)
	bestX.to_zero()
	
	bestY=EMData(options.nsoln,1,1)
	bestY.to_zero()
	
	bestZ=EMData(options.nsoln,1,1)
	bestZ.to_zero()


	ral,dal = options.ralt,options.dalt
	rap,dap = options.rphi, options.dphi
	daz = options.daz

	altarray=[]
	alt=0
	while alt <= ral:
		altarray.append(alt)
		alt=alt+dal
	
	azarray=[]
	az=-180
	while az < 180:
		azarray.append(az)
		az=az+daz

	rarad=rap
	darad=dap
	maxfrac=0.
	norm = 1.0/(box*box*box)
	for altrot in altarray:
		print "Trying rotation %f"%(altrot)
		if (altrot==0.):
			azrot=0.
			phirot=-azrot-rarad
			minnum=10000000
			maxnum=0
			while phirot <= -azrot+rarad:
				t = Transform({"type":"eman","az":azrot,"phi":phirot,"alt":altrot})
				dMRC= probeMRC.process("math.transform",{"transform":t}) # more efficient than copying than transforming
				currentCCF=tomoccf(targetMRC,dMRC)
				scalar=ccfFFT(currentCCF,options.thresh,box)
				if scalar>maxnum:
					maxnum=int(scalar)
				if scalar<minnum:
					minnum=int(scalar)
				scalar=scalar/(box*box*box/2.)
				bestCCF=updateCCF(bestCCF,bestALT,bestAZ,bestPHI,bestX,bestY,bestZ,altrot,azrot,phirot,currentCCF,scalar,options.nsoln,options.searchx,options.searchy,options.searchz)
				phirot=phirot+darad
	
		else:
			for azrot in azarray:
				phirot=-azrot-rarad
				#print "Trying rotation %f %f"%(altrot, azrot)
				while phirot <= -azrot+rarad:
					t = Transform({"type":"eman","az":azrot,"phi":phirot,"alt":altrot})
					dMRC= probeMRC.process("math.transform",{"transform":t}) # more efficient than copying than
					currentCCF=tomoccf(targetMRC,dMRC)
					scalar=ccfFFT(currentCCF,options.thresh,box)
					if scalar>maxnum:
						maxnum=int(scalar)
					if scalar<minnum:
						minnum=int(scalar)
					scalar=scalar/(box*box*box/2.)
					bestCCF=updateCCF(bestCCF,bestALT,bestAZ,bestPHI,bestX,bestY,bestZ,altrot,azrot,phirot,currentCCF,scalar,options.nsoln,options.searchx,options.searchy,options.searchz)
					phirot=phirot+darad
	
	
	print minnum,maxnum, float(maxnum)/float(minnum)


	out=open("log-s3-%s%s.txt"%(argv[1],argv[2]),"w")
	peak=0
	while peak < 10:
		ALT=str(bestALT.get(peak))
		AZ=str(bestAZ.get(peak))
		PHI=str(bestPHI.get(peak))
		COEFF=str(bestCCF.get(peak))
		LOC=str( ( (bestX.get(peak)),(bestY.get(peak)),(bestZ.get(peak) ) ) )
		line="Peak %d rot=( %s, %s, %s ) trans= %s coeff= %s\n"%(peak,ALT,AZ,PHI,LOC,COEFF)
		out.write(line)
		peak=peak+1
		
	out.close()
	
def tomoccf(targetMRC,probeMRC):
	ccf=targetMRC.calc_ccf(probeMRC)
	# removed a toCorner...this puts the phaseorigin at the left corner, but we can work around this by
	# by using EMData.calc_max_location_wrap (below)
	return (ccf)

def ccfFFT(currentCCF, thresh, box):
	tempCCF = currentCCF.do_fft()
	#tempCCF.ri2ap()
	tempCCF.process_inplace("threshold.binary.fourier",{"value":thresh}) 
	mapSum=tempCCF["mean"]*box*box*box
	return(mapSum)

def updateCCF(bestCCF,bestALT,bestAZ,bestPHI,bestX,bestY,bestZ,altrot,azrot,phirot,currentCCF,scalar,n,searchx,searchy,searchz):
	best = currentCCF.calc_max_location_wrap(searchx,searchy,searchz)
	xbest = best[0]
	ybest = best[1]
	zbest = best[2]
	bestValue = currentValue = currentCCF.get(xbest,ybest,zbest)/scalar
	inlist=0
	while inlist < n:
		if  bestValue > bestCCF.get(inlist):
			swlist=n-1
			while swlist >inlist:
				#print swlist
				bestCCF.set(swlist,bestCCF.get(swlist-1,))
				bestALT.set(swlist,bestALT.get(swlist-1))
				bestAZ.set(swlist,bestAZ.get(swlist-1))
				bestPHI.set(swlist,bestPHI.get(swlist-1))
				bestX.set(swlist,bestX.get(swlist-1))
				bestY.set(swlist,bestY.get(swlist-1))
				bestZ.set(swlist,bestZ.get(swlist-1))
				swlist=swlist-1
			bestCCF.set(inlist,bestValue)
			bestALT.set(inlist,altrot)
			bestAZ.set(inlist,azrot)
			bestPHI.set(inlist,phirot)
			bestX.set(inlist,xbest)
			bestY.set(inlist,ybest)
			bestZ.set(inlist,zbest)
			break
		inlist=inlist+1
	
	#bestCCF.update() uneccessary?
	#bestALT.update()
	#bestAZ.update()
	#bestPHI.update()
	return(bestCCF)


def check(options,args):
	#should write a function to check the inputs, such as positive range, delta less than range etc
	#should also check that the file names exist
	
	error = False
	if options.nsoln <= 0:
		error = True
		print "Error - nsoln must be greater than 0. Suggest using 10"
	
	if options.ralt <= 0:
		error = True
		print "Error - ralt must be greater than 0"
	
	# etc....
	
	return error

# If executed as a program
if __name__ == '__main__':
	main() 
