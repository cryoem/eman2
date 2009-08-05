#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from optparse import OptionParser
from EMAN2 import file_exists,EMData,E2init,E2progress,E2end,EMANVERSION,check_eman2_type_string,EMUtil,remove_file,get_file_tag
import EMAN2
from EMAN2db import EMTask

def check_tomo_options(options):
	'''
	Checks options.align, options.aligncmp, options.cmp and options.nsoln
	Checks options.ralign if it is not None, and if not None, checks options.raligncmp
	Used by e2tomoaverage
	@param options - as returned by (options, args) = parser.parse_args() in e2tomoallvall and e2tomoaverage 
	@param args - a list of file names - as returned by (options, args) = parser.parse_args() in e2tomoallvall and e2tomoaverage 
	'''
	error = []
	if options.align == None:
		error.append("Error - you have to supply the align option")
	else:
		e = check_eman2_type_string(options.align,EMAN2.Aligners,"Aligners")
		if e != None:
			error.append(e)
	
	if options.aligncmp == None:
		error.append("Error - you have to supply the aligncmp option")
	else:
		e = check_eman2_type_string(options.aligncmp,EMAN2.Cmps,"Cmps")
		if e != None:
			error.append(e)
			
	if options.ralign != None: # not strictly necessary
		e = check_eman2_type_string(options.ralign,EMAN2.Aligners,"Aligners")
		if e != None:
			error.append(e)
		if options.raligncmp == None:
			error.append("Error - if you supply the ralign argument, you have to supply the raligncmp option")
		else:
			e = check_eman2_type_string(options.raligncmp,EMAN2.Cmps,"Cmps")
			if e != None: error.append(e)
		
	
	if options.cmp == None:
		error.append("Error - you have to supply the cmp option")
	else:
		e = check_eman2_type_string(options.cmp,EMAN2.Cmps,"Cmps")
		if e != None:
			error.append(e)
	
	if options.nsoln <= 0:
		error.append( "nsoln must be greater than" )
	
	return error

def check_options(options,args):
	'''
	check function for e2tomoallvall
	'''
	error = []
	error.extend(check_tomo_options(options))
	
	if len(args) < 2:
		error.append("Error - to average you must supply atleast two images")
	else:
		for arg in args:
			if not file_exists(arg):
				error.append("%s does not exist" %arg)
	
	if file_exists(options.output):
		if not options.force:
			error.append("Error - output file exists. Remove it or supply the force argument")
		else:
			remove_file(options.output)
	
	return error
	
class EMTomoAllVAll:
	'''
	This class breaks the jobs of the boot-strapped average generation procedure
	so that they can be run in parallel. 
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
		
		self.output_images = None # these will be a list of images to write the output to
		self.init_output()
		
	def init_output(self):
		'''
		Called internally to initialize output images - they are stored in memory
		'''
		
		self.output_images = []
		e = EMData(len(self.files),len(self.files))
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
		Function to perform the main work
		'''
		alignment_jobs = []
		for i in range(len(self.files)):
			for j in range(i+1,len(self.files)):
				alignment_jobs.append([i,j])
		from e2tomoaverage import EMTomoAlignments
		alignments_manager = EMTomoAlignments(self.options)
		alignments_manager.execute(alignment_jobs, self.files,self)
				
		self.finalize_writing()

	def process_output(self,results):
		'''
		store the result of an alignment
		'''
		output_writer = EMTomoOutputWriter() # note - has a relationship
		output_writer.process_output(results, self.output_images, self.files)

	def finalize_writing(self):
		'''
		Writes the output images to disk
		'''
		
		for i,image in enumerate(self.output_images):
			image.set_attr("file_names",self.files)
			image.write_image(self.options.output,i)
		
		
class EMTomoOutputWriter:
	'''
	common functionality to EMTomoAllVAll and e2tomoaverage.EMBootStrappedAverages
	Supplies the process_output function
	'''
	def __init__(self): pass
	
	def process_output(self,results,images,files):
		'''
		@param results a dictionary that was returned by an EMTomoAlignTask
		@parm images a list of seven images that will store the alignment results
		@param files a list of files - the target and probe indices from the results dict are used to get file names and write results to disk (see latter part of function)
		'''
		cmp = results["cmp"]
		ali = results["ali"]
		target_idx = results["target_idx"]
		probe_idx = results["probe_idx"]
		
		a = ali.get_params("eman")
		images[0].set(target_idx,probe_idx,cmp)
		images[1].set(target_idx,probe_idx,a["tx"])
		images[2].set(target_idx,probe_idx,a["ty"])
		images[3].set(target_idx,probe_idx,a["tz"])
		images[4].set(target_idx,probe_idx,a["az"])
		images[5].set(target_idx,probe_idx,a["alt"])
		images[6].set(target_idx,probe_idx,a["phi"])
		
		
		# here we put the inverted alignment as well, yes this is redundant but it 
		# seams like the most logical way to fit with the current eman2 philosophies
		ali_inv = ali.inverse()
		ai = ali_inv.get_params("eman")
		
		images[0].set(probe_idx,target_idx,cmp)
		images[1].set(probe_idx,target_idx,ai["tx"])
		images[2].set(probe_idx,target_idx,ai["ty"])
		images[3].set(probe_idx,target_idx,ai["tz"])
		images[4].set(probe_idx,target_idx,ai["az"])
		images[5].set(probe_idx,target_idx,ai["alt"])
		images[6].set(probe_idx,target_idx,ai["phi"])
		
		all_solns = results["all_solns"]
		if len(all_solns) > 1:
			target_name = get_file_tag(files[target_idx])
			probe_name = get_file_tag(files[probe_idx]) 
			out=file("log-s3-%s_%s.txt"%(target_name,probe_name),"w")
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
				LOC=str( ( (params["tx"]),(params["tx"]),(params["tx"] ) ) )
				line="Peak %d rot=( %f, %f, %f ) trans= %s coeff= %s\n"%(peak,ALT,AZ,PHI,LOC,COEFF)
				out.write(line)
				peak=peak+1
				
			out.close()

class EMTomoAlignTask:
	'''
	A class the knows how to align two 3D volumes
	'''
	def __init__(self,data=None,options=None):
		self.data = data
		self.align_data = data["align"]
		self.align_cmp_data = data["aligncmp"]
		self.cmp_data = data["cmp"]
		if data.has_key("ralign"):
			self.ralign_data = data["ralign"]
			self.ralign_cmp_data = data["raligncmp"]
		else:
			self.ralign_data = None
		self.using_cuda = data["using_cuda"]
		self.nsoln = data["nsoln"]
	def execute(self,progress_callback):
		'''
		Called to perform class averaging 
		May boot strap the original average, iteratively refines averages, aligns final average to ref 
		'''
		from EMAN2db import db_open_dict
		progress_callback(0)
		
		probe = self.data["probe"]
		target = self.data["target"]
		
		return self.align(probe,target,progress_callback)
		
	def align(self,probe,target,progress_callback):
		'''
		The tomographic alignment routine here, in one small, convenient place
		This function is used by every single tomographic alignment program
		'''
		if self.using_cuda:
			target.set_gpu_rw_current()
			target.cuda_lock() # locking it prevents if from being overwritten
			probe.set_gpu_rw_current()
			probe.cuda_lock()
		
		progress = 0.0
		max_progress = 3
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))

		solns = probe.xform_align_nbest(self.align_data[0],target,self.align_data[1],self.nsoln,self.align_cmp_data[0],self.align_cmp_data[1])
		#ali = probe.align(self.align_data[0],target,self.align_data[1],self.align_cmp_data[0],self.align_cmp_data[1])
		
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))
		if self.ralign_data != None:
			for s in solns:
				self.ralign_data[1]["xform.align3d"] = s["xform.align3d"]
				aligned = probe.align(self.ralign_data[0],target,self.ralign_data[1],self.ralign_cmp_data[0],self.ralign_cmp_data[1])
				s["xform.align3d"] = aligned.get_attr("xform.align3d")
				s["score"] = aligned.cmp(self.cmp_data[0],target,self.cmp_data[1])
			
			# the alignment might be better are refine alignment!
			solns.sort(alignment_score_sort)
			solns.reverse()
#	
		if self.using_cuda:
			# this is just for consistency. The EMData destructor would have done this automatically anyway
			target.cuda_unlock()
			probe.cuda_unlock()
		
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))
	
		results = {}
		results["cmp"] = solns[0]["score"]
		results["ali"] = solns[0]["xform.align3d"]
		results["target_idx"] = self.data["target_idx"]
		results["probe_idx"] = self.data["probe_idx"]
		results["all_solns"] = solns
		
		return results
	
def alignment_score_sort(left,right):
	'''
	When we use refine alignment the similarity scores change
	and therefore the order of the solutions might need altering
	'''
	c1 = left["score"]
	c2 = right["score"]
	if c1 > c2: return -1
	elif c1 == c2: return 0
	else: return 1
			
class EMTomoAlignTaskDC(EMTask):
	'''
	A class the knows how to align two 3D volumes using parallel DC framework
	'''
	def __init__(self,command="e2tomoallvall",data=None,options=None):
		EMTask.__init__(self,command,data,options)
		# need this import for parallelism to work - somebody fixme?
		from e2tomoallvall import EMTomoAlignTask
		self.align_task = EMTomoAlignTask(data=data)
#		self.align_data = align_data
#		self.ralign_data = ralign_data
#		self.cmp_data = cmp_data
	
	def execute(self,progress_callback):
		'''
		Called to perform class averaging 
		May boot strap the original average, iteratively refines averages, aligns final average to ref 
		'''
		from EMAN2db import db_open_dict
		progress_callback(0)
		
		cache_name=self.data["probe"][1]
		cache=db_open_dict(cache_name)
		probe = cache[self.data["probe"][2]]
		
		cache_name=self.data["target"][1]
		cache=db_open_dict(cache_name)
		target = cache[self.data["target"][2]]
		
		return self.align_task.align(probe, target, progress_callback)
		
	  	
import os,sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image1> <image2> <image3> <image4> ....
	
Boot straps an initial probe doing all versus all alignment of the input images

"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rt.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default="rt.3d.grid")
	parser.add_option("--aligncmp",type="string",help="The comparator used for determing the best initial alignment", default="dot.tomo:threshold=0")
	parser.add_option("--cmp",type="string",help="The comparator used to obtain the final similarity", default="dot.tomo:threshold=0")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparator used for determing the refined alignment", default="dot.tomo:threshold=0")
	parser.add_option("--output",type="string",default="e2tomoallvall.hdf",help="The output image which will store the results matrix")
	parser.add_option("--force",action="store_true",default=False,help="Force overwrite the output file if it already exists")
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	parser.add_option("--nsoln", default=1, type="int",help="If supplied and greater than 1, the nsoln-best alignments will be written to a text file. This is useful for debug but may be left unspecified")
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
	
	logger=E2init(sys.argv)

	module = EMTomoAllVAll(args,options,logger)
	module.execute()
	
	E2end(logger)
	
	
# If executed as a program
if __name__ == '__main__':
	main() 
	