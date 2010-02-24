#!/usr/bin/env python

#
# Author: Matthew Baker, 10/2005, modified 02/2006 by MFS  
# ported and refactored using EMAN2 by David Woolford October 6th 2008
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

from e2tomoaverage import check_tomo_options

tomo_hunter_path_root = "tomo_hunt"

def check_options(options,filenames):
	'''
	Check the parser options
	Should probably be made into a class
	@return a list of error messages
	'''
	error_messages = []
	
	error_messages.extend(check_tomo_options(options))
	
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
									
	return error_messages

def gen_average(options,args,logid=None):
	'''
	Uses information stored in a dictionary in the local database to produce an average
	'''
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
			
		image.process_inplace("xform",{"transform":t})
		if average == None: average = image
		else: average = average + image
		if logid: E2progress(logid,(i+1)/float(total_prog))
		
	average.mult(1.0/len(args))
		
	average.write_image(options.avgout,0)
	
	if options.dbls:
		pdb = db_open_dict("bdb:project")
		db = pdb.get(options.dbls,dfl={})
		if isinstance(db,list): # this is for back compatibility - it used to be a list, now it's a dict
			d = {}
			for name in db:
				s = {}
				s["Original Data"] = name
				d[name] = s
			db = d
		s = {}
		s["Original Data"] = options.avgout
		db[options.avgout] = s
		pdb[options.dbls] = db

def get_ali_data(filename,probe,aliset):
	'''
	Get alignment data from the local datase
	'''
	from emtprworkflow import EMProbeAliTools
	from emsprworkflow import EMPartSetOptions
	
	#EMProjectListCleanup.clean_up_filt_particles(self.project_list)
	db = db_open_dict("bdb:project",ro=True)
	db_map = db.get("global.tpr_ptcls_ali_dict")
	if db_map == None:
		return None # calling function will barf
	
	ptcl_opts = EMPartSetOptions("global.tpr_ptcls_dict")
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

class EMTomoHunter:
	'''
	This class oversees the execution jobs as supplied by e2tomohunter
	This boils down to figuring out all of the alignments that need to performed and
	then executing them in the current context or using parallelism
	'''
	def __init__(self,files,options,logger):
		'''
		@param options - the options returned by the call to (options, args) = parser.parse_args() 
		@param args - a list of image names - args is that which is returned by the call to (options, args) = parser.parse_args()
		'''
		self.options = options
		self.files = files
		self.logger = logger
		self.using_cuda = EMUtil.cuda_available() and options.cuda
		self.output_images = None
		self.total_jobs = len(files) # used for recording progress
		self.jobs_completed = 0 # used for recording progress
		self.init_output()
	
	def init_output(self):
		'''
		Initialize memory for storing alignment results
		'''
		self.output_images = []
		e = EMData(len(self.files),1)
		e.to_zero()
		
		# need to store correlation score (or cmp score)
		self.output_images.append(e)
		
		# need to store translation x
		self.output_images.append(e.copy())
		# need to store translation y
		self.output_images.append(e.copy())
		# need to store translation z
		self.output_images.append(e.copy())
		
		# need to store translation az
		self.output_images.append(e.copy())
		# need to store translation alt
		self.output_images.append(e.copy())
		# need to store translation phi
		self.output_images.append(e.copy())
	
	def execute(self):
		'''
		Execute the tomohunter jobs
		'''
		if self.logger != None:
			E2progress(self.logger,0.0)
		alignment_jobs = []
		
		probe_idx = len(self.files)
		for i in range(len(self.files)):
			alignment_jobs.append([probe_idx,i])
		
		self.files.append(self.options.probe) # EMTomoAlignments needs the name of the probe
			
		from e2tomoaverage import EMTomoAlignments
		alignments_manager = EMTomoAlignments(self.options)
		alignments_manager.execute(alignment_jobs, self.files,self)
		
		self.finalize_writing()
		
	def process_output(self,results):
		cmp = results["cmp"]
		ali = results["ali"]
		target_idx = results["target_idx"]
		probe_idx = results["probe_idx"]
		
		a = ali.get_params("eman")
		self.output_images[0].set(target_idx,0,cmp)
		self.output_images[1].set(target_idx,0,a["tx"])
		self.output_images[2].set(target_idx,0,a["ty"])
		self.output_images[3].set(target_idx,0,a["tz"])
		self.output_images[4].set(target_idx,0,a["az"])
		self.output_images[5].set(target_idx,0,a["alt"])
		self.output_images[6].set(target_idx,0,a["phi"])
		
		all_solns = results["all_solns"]
		if len(all_solns) > 1:
			target_name = get_file_tag(self.files[target_idx])
			probe_name = get_file_tag(self.files[probe_idx])
			out=file("log-s3-%s_VS_%s.txt"%(target_name,probe_name),"w")
			peak = 0
			for d in all_solns:
				t = d["xform.align3d"]
				# inverting because the probe was aligned to the target
				t = t.inverse()
				params = t.get_params("eman")
				ALT=params["alt"]
				AZ=params["az"]
				PHI=params["phi"]
				COEFF=str(d["score"])
				LOC=str( ( (params["tx"]),(params["ty"]),(params["tz"] ) ) )
				line="Peak %d rot=( %f, %f, %f ) trans= %s coeff= %s\n"%(peak,ALT,AZ,PHI,LOC,COEFF)
				out.write(line)
				peak=peak+1
				
			out.close()
		
		options = self.options
		if options.dbls:
			pdb = db_open_dict("bdb:project")
			db = pdb.get(options.dbls,dfl={})
			if db == None: db = {}
			results = []
			for d in all_solns:
				t = d["xform.align3d"]
				# inverting because the probe was aligned to the target
				t = t.inverse()
				results.append(t)
			
			arg = self.files[target_idx]
			if db.has_key(arg):d = db[arg]
			else:d = {}
			d[options.probe] = results
			db[arg] = d
			pdb[options.dbls] = db
			
		if self.logger != None:
			self.jobs_completed += 1.0 # used for recording progress
			E2progress(self.logger,self.jobs_completed/self.total_jobs)
			
	def finalize_writing(self):
		'''
		Writes the output images to disk
		'''
		path = numbered_path(tomo_hunter_path_root,True)
		output_name =  numbered_bdb("bdb:"+path+"#simmx")
		for i,image in enumerate(self.output_images):
			image.set_attr("file_names",self.files)
			image.write_image(output_name,i)
				
		
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <images to be aligned> [options]"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--probe",type="string",help="The probe. This is the model that the input images will be aligned to", default=None)
	parser.add_option("--n",type="int",help="0 or 1, multiplication by the reciprocal of the boxsize", default=1)
	parser.add_option("--dbls",type="string",help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	parser.add_option("--aliset",type="string",help="Supplied with avgout. Used to choose different alignment parameters from the local database. Used by workflow.", default=None)
	parser.add_option("--avgout",type="string",help="If specified will produce an averaged output, only works if you've previously run alignments", default=None)
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rt.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default="rt.3d.grid")
	parser.add_option("--aligncmp",type="string",help="The comparator used for determing the best initial alignment", default="dot.tomo:threshold=0")
	parser.add_option("--cmp",type="string",help="The comparator used to obtain the final similarity", default="dot.tomo:threshold=0")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparator used for determing the refined alignment", default="dot.tomo:threshold=0")
	parser.add_option("--nsoln",type="int",help="The number of solutions to report", default=10)
	parser.add_option("--shrink",type="int",help="Shrink the data as part of the alignment - for speed purposes but at the potential loss of accuracy",default=None)
	parser.add_option("--filter",type="string",help="The name and parameters of an EMAN2 processor. Will be applied prior to shrinking.",default=None)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	if EMUtil.cuda_available():
		parser.add_option("--cuda",action="store_true",help="GPU acceleration using CUDA. Experimental", default=False)
   
	(options, args) = parser.parse_args()
	
	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(1)
		
	logid=E2init(sys.argv)
	
	if options.avgout: # unfortunately this functionality is part of this script.
		gen_average(options,args,logid)
		exit(0)
	
	prog = 0
	total_prog = len(args)
	E2progress(logid,0.0)
	
	tomohunter = EMTomoHunter(args,options,logid)
	tomohunter.execute()

	if using_cuda(options):
		probeMRC.cuda_unlock()
	E2progress(logid,1.0) # just make sure of it
		
	
	E2end(logid)
	
def using_cuda(options):
	return EMUtil.cuda_available() and options.cuda
	
def print_info(image,first_line="Information"):
	
	print first_line
	print "   mean:	   %f"%(image.get_attr("mean"))
	print "   sigma:	  %f"%(image.get_attr("sigma"))
	


# If executed as a program
if __name__ == '__main__':
	main() 
