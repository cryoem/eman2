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
	if not options.probe or not file_exists(options.probe):
		error_messages.append("You have to specify a valid probe")
	
	if len(filenames) < 1:
		error_messages.append("You must specify input files")
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

def gen_average(options,args,logid=None):
	
	project_list = "global.tpr_ptcls_ali_dict"
	db = db_open_dict("bdb:project",ro=True)
	db_map = db.get(project_list,dfl={})
	
	probes = []
	probes_data = {}
	
	average = None
	
	prog = 0
	total_prog = len(args)
	if logid: E2progress(logid,0.0)
	
	for i,arg in enumerate(args):
		image = EMData(arg,0)
		if not options.aliset:
			t = db_map[arg][options.probe][0]
		else:
			t = get_ali_data(arg,options.probe,options.aliset)
			if t == None:
				raise RuntimeError("An error occured trying to retrieve the alignment data using the given ali set")
			
		image.process_inplace("math.transform",{"transform":t})
		if average == None: average = image
		else: average = average + image
		if logid: E2progress(logid,(i+1)/float(total_prog))
		
	average.mult(1.0/len(args))
		
	average.write_image(options.avgout,0)
	
	if options.dbls:
		pdb = db_open_dict("bdb:project")
		db = pdb.get(options.dbls,dfl=[])
		if db == None: db = []
		if db.count(options.avgout) == 0:
			db.append(options.avgout)
			pdb[options.dbls] = db

def get_ali_data(filename,probe,aliset):
	from emtprworkflow import EMProbeAliTools
	from emsprworkflow import EMPartStackOptions
	
	#EMProjectListCleanup.clean_up_filt_particles(self.project_list)
	db = db_open_dict("bdb:project",ro=True)
	db_map = db.get("global.tpr_ptcls_ali_dict")
	if db_map == None:
		return None # calling function will barf
	
	ptcl_opts = EMPartStackOptions("global.tpr_ptcls","global.tpr_filt_ptcls_map")
	particles_map, particles_name_map, choices, name_map = ptcl_opts.get_particle_options()
	tls = EMProbeAliTools()
	probe_set_map,probe_and_ali,probe_name_map = tls.accrue_data()
	
	base_set = probe_set_map[get_file_tag(probe)][aliset]
	ptcl_base_set = [name_map[name] for name in base_set]
	
	base_name = name_map[filename]
	
	for i in xrange(0,len(base_set)):
		if base_name == ptcl_base_set[i]:
			dct = db_map[base_set[i]]
			if dct.has_key(probe):
				alis = dct[probe]
				return alis[0]
			
	return None
	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <probe> <image to be aligned> [options]"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--probe",type="string",help="The probe. This is the model that the input images will be aligned to", default=None)
	parser.add_option("--dalt",type="float",help="Altitude delta", default=10.0)
	parser.add_option("--ralt",type="float",help="Altitude range", default=180.0)
	parser.add_option("--dphi",type="float",help="Phi delta", default=10.0)
	parser.add_option("--rphi",type="float",help="Phi range", default=180.0)
	parser.add_option("--raz",type="float",help="Azimuth range", default=360.0)
	parser.add_option("--daz",type="float",help="Azimuth delta", default=10.0)
	parser.add_option("--thresh",type="float",help="Threshold", default=1.0)
	parser.add_option("--nsoln",type="int",help="The number of solutions to report", default=10)
	parser.add_option("--searchx",type="int",help="The maximum search distance, x direction", default=5)
	parser.add_option("--searchy",type="int",help="The maximum search distance, y direction", default=5)
	parser.add_option("--searchz",type="int",help="The maximum search distance, z direction", default=5)
	parser.add_option("--n",type="int",help="0 or 1, multiplication by the reciprocal of the boxsize", default=1)
	parser.add_option("--dbls",type="string",help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	parser.add_option("--aliset",type="string",help="Supplied with avgout. Used to choose different alignment parameters from the local database. Used by workflow.", default=None)
	parser.add_option("--avgout",type="string",help="If specified will produce an averaged output, only works if you've previously run alignments", default=None)
	if EMUtil.cuda_available():
		parser.add_option("--cuda",action="store_true",help="GPU acceleration using CUDA. Tantalizing glimpse into the future", default=False)
   
	(options, args) = parser.parse_args()
	
	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(1)
		
	logid=E2init(sys.argv)
	
	if options.avgout:
		gen_average(options,args,logid)
		exit(0)
	
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
	maxfrac=0
	
	nx,ny,nz = gimme_image_dimensions3D(args[0])
	box=nx
	if options.n:
		norm = 1.0/(box*box*box)
	else:
		norm = 1.0
		
	prog = 0
	total_prog = (len(altarray)-1)*len(azarray)*len(args)+1
	E2progress(logid,0.0)
	
	for arg in args:
		
		
		targetMRC =EMData(arg,0)
		print_info(targetMRC,"Target Information")
		
		probeMRC=EMData(options.probe,0)
		print_info(targetMRC,"Probe Information")
		
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
				
				prog += 1.0
				E2progress(logid,prog/total_prog)
		
			else:
				for azrot in azarray:
					phirot=-azrot-rarad
					#print "Trying rotation %f %f"%(altrot, azrot)
					while phirot <= -azrot+rarad:
						t = Transform({"type":"eman","az":azrot,"phi":phirot,"alt":altrot})
						probeMRC.set_gpu_rw_current()
						targetMRC.set_gpu_rw_current()
						dMRC= probeMRC.process("math.transform",{"transform":t}) # more efficient than copying than
						dMRC.set_gpu_rw_current()
						currentCCF=tomoccf(targetMRC,dMRC)
						scalar=ccfFFT(currentCCF,options.thresh,box)
						if scalar>maxnum:
							maxnum=int(scalar)
						if scalar<minnum:
							minnum=int(scalar)
						scalar=scalar/(box*box*box/2.)
						bestCCF=updateCCF(bestCCF,bestALT,bestAZ,bestPHI,bestX,bestY,bestZ,altrot,azrot,phirot,currentCCF,scalar,options.nsoln,options.searchx,options.searchy,options.searchz)
						phirot=phirot+darad
					
					prog += 1.0
					E2progress(logid,prog/total_prog)
		
		
		print minnum,maxnum, float(maxnum)/float(minnum)
	
		out=file("log-s3-%s_%s.txt"%(get_file_tag(arg),get_file_tag(options.probe)),"w")
		peak=0
		
		db = None
		pdb = None
		if options.dbls:
			pdb = db_open_dict("bdb:project")
			db = pdb.get(options.dbls,dfl={})
			if db == None: db = {}
			results = []
	
	
		while peak < 10:
			ALT=bestALT.get(peak)
			AZ=bestAZ.get(peak)
			PHI=bestPHI.get(peak)
			COEFF=str(bestCCF.get(peak))
			LOC=str( ( (bestX.get(peak)),(bestY.get(peak)),(bestZ.get(peak) ) ) )
			line="Peak %d rot=( %f, %f, %f ) trans= %s coeff= %s\n"%(peak,ALT,AZ,PHI,LOC,COEFF)
			out.write(line)
			peak=peak+1
			
			if options.dbls:
				t = Transform({"type":"eman","alt":ALT,"phi":PHI,"az":AZ})
				t.set_trans(bestX.get(peak),bestY.get(peak),bestZ.get(peak))
				results.append(t)
				
		if options.dbls:
			if db.has_key(arg):d = db[arg]
			else:d = {}
			d[options.probe] = results
			db[arg] = d
			pdb[options.dbls] = db
			
		out.close()
	E2progress(logid,1.0) # just make sure of it
		
	
	E2end(logid)
	
def print_info(image,first_line="Information"):
	
	print first_line
	print "   mean:	   %f"%(image.get_attr("mean"))
	print "   sigma:	  %f"%(image.get_attr("sigma"))
	
	
def tomoccf(targetMRC,probeMRC):
	ccf=targetMRC.calc_ccf_cuda(probeMRC,False,False)
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
