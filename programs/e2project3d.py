#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Additional Author: David Woolford 2007-2008 (woolford@bcm.edu)
#
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
# References
# 1. Baldwin, P.R. and Penczek, P.A. 2007. The Transform Class in SPARX and EMAN2. J. Struct. Biol. 157, 250-261.
# 2. http://blake.bcm.edu/emanwiki/EMAN2/Symmetry

# 
# 1. Asymmetric units are accurately demarcated and covered by the projection algorithm. This can
# be tested using 3D plotting in Matlab.
# 2. By default the entire asymmetric is projected over, but excluding the mirror portion is supported.
# The accuracy of the demarcation of the asymmetric unit can be tested by using the --verifymirror
# argument, which subtracts (mirrored) mirror projections in the local asymmetric unit from the 
# the equivalent original projections - the result should be zero but this is not the case due
# to interpolation differences, however the images should have a mean of about zero and a standard
# deviation that is relatively small. Visual inspection of the results (in default output files) can also help.
# 3. Random orientation generation is supported - all three euler angles are random
# 4. Perturbation of projections generated in the asymmetric unit is supported in regular reconstructions runs.
# 5. The user is able to generate projections that include in-plane or "phi" rotations
# and this is achieved using the phitoo argument
# 6. The user can smear in-plane projections when the "smear" argument is specified in addition to the "phitoo" argument
# what else..????

import sys, math, os, random
from EMAN2 import *
from optparse import OptionParser
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi
DEBUG = False
WEN_JIANG = False
EMAN1_OCT = False
MIRROR_DEBUG = True
NO_MIRROR = False


class EMParallelProject3D:
	def __init__(self,options,args,logger=None):
		'''
		@param options the options produced by (options, args) = parser.parse_args()
		@param args the options produced by (options, args) = parser.parse_args()
		@param logger and EMAN2 logger, i.e. logger=E2init(sys.argv)
		assumes you have already called the check function.
		'''
		self.options = options
		self.args = args
		self.logger = logger
		
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		self.num_cpus = etc.cpu_est()
		self.num_cpus = 4
		
		self.__task_options = None
	
	def __init_memory(self,options):
		'''
		
		'''
		sym_object = parsesym(options.sym)
		[og_name,og_args] = parsemodopt(options.orientgen)
		self.eulers = sym_object.gen_orientations(og_name, og_args)
	
	def __get_task_options(self,options):
		
		if self.__task_options == None:
			d = {}
			d["projector"] = parsemodopt(options.projector)
			self.__task_options = d
			
		return self.__task_options
	
	def execute(self):
		
		if len(self.options.parallel) > 2 and self.options.parallel[:2] == "dc":
			self.__init_memory(self.options)
			
			num_tasks = self.num_cpus
			# In the worst case we can only spawn as many tasks as there are eulers
			if self.num_cpus > len(self.eulers): num_tasks = len(self.eulers)
			
			eulers_per_task = len(self.eulers)/num_tasks
			resid_eulers = len(self.eulers) - eulers_per_task*num_tasks # we can distribute the residual evenly
			
			first = 0
			self.task_customers = []
			self.tids = []
			for i in xrange(0,num_tasks):
				last = first+eulers_per_task
				if resid_eulers > 0:
					last +=1
					resid_eulers -= 1
					
				tmp_eulers = self.eulers[first:last]
				indices = range(first,last)
				
				data = {}
				data["input"] = ("cache",self.args[0],0)
				data["eulers"] = tmp_eulers
				data["indices"] = indices
				
				task = EMProject3DTaskDC(data=data,options=self.__get_task_options(self.options))
				
				from EMAN2PAR import EMTaskCustomer
				etc=EMTaskCustomer(self.options.parallel)
				#print "Est %d CPUs"%etc.cpu_est()
				tid=etc.send_task(task)
					#print "Task submitted tid=",tid
					
				self.task_customers.append(etc)
				self.tids.append(tid)
				
				first = last
				
			while 1:
				if len(self.task_customers) == 0: break
				print len(self.task_customers),"projection tasks left in main loop"
				for i in xrange(len(self.task_customers)-1,-1,-1):
					task_customer = self.task_customers[i]
					tid = self.tids[i] 
					st=task_customer.check_task((tid,))[0]
					if st==100:
						
						self.task_customers.pop(i)
						self.tids.pop(i)
	
						rslts = task_customer.get_results(tid)
						self.__write_output_data(rslts[1])
						if self.logger != None:
							E2progress(self.logger,1.0-len(self.task_customers)/float(num_tasks))
							if self.options.verbose: 
								print "%d/%d\r"%(num_tasks-len(self.task_customers),num_tasks)
								sys.stdout.flush()
				
				time.sleep(5)
		else:
			raise NotImplementedError("The parallelism option you specified (%s) is not suppored" %self.options.parallel )
				
	def __write_output_data(self,rslts):
		for idx,image in rslts["projections"].items():
			image.write_image(self.options.outfile,idx)
			
	
from EMAN2db import EMTask
class EMProject3DTaskDC(EMTask):
	def __init__(self,command="e2project3d",data=None,options=None):
		EMTask.__init__(self,command,data,options)
		
		# data has these keys:
		# input - which is the name of the threed model - a Task-style cache
		# eulers - a list of Transforms with which to generate projections
		# indices - indices that correspond to the ordering, has to be the same length as eulers. Returned to calling routine - use to write output in order. Not sure about this.
		
		# options has these keys
		# projector - [string,dict] convention
		
		self.projections = {} # key will be index, value will be projection
		
	def execute(self):
		input_name=self.data["input"][1]
		threed_image = EMData(input_name,self.data["input"][2]) # so the idx should always be 0, but it doesn't have to be
		
		options = self.options
		eulers = self.data["eulers"]
		indices = self.data["indices"]
		projector,projector_opts = options["projector"][0],options["projector"][1]
		for i in xrange(len(eulers)):
			euler = eulers[i]
			projector_opts["transform"] = euler
			projection = threed_image.project(projector,projector_opts)
			projection.set_attr("xform.projection",euler)
			projection.set_attr("ptcl_repr",0)
			self.projections[indices[i]] = projection
		
	
	def get_return_data(self):
		d = {}
		d["projections"] = self.projections
		return d
	

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog image [options] 
	Projects in real space over the asymmetric unit using the angular separation as specified by prop and the symmetry as specified by sym."""
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos")
	parser.add_option("--orientgen", dest="orientgen", help="The orientation generator to use. See e2help.py orientgen")
	parser.add_option("--outfile", dest = "outfile", default = "e2proj.hdf", help = "Output file. Default is 'e2proj.img'")
	# add --perturb
	parser.add_option("--smear", dest = "smear", type = "int", default=0,help="Used in conjunction with --phitoo, this will rotationally smear between phi steps. The user must specify the amount of smearing (typically 2-10)")
	parser.add_option("--projector", dest = "projector", default = "standard",help = "Projector to use")
	#parser.add_option("--verifymirror",action="store_true",help="Used for testing the accuracy of mirror projects",default=False)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_option("--append", "-a",dest="append",default=False, action="store_true",help="Append to the output file")
	parser.add_option("--verbose","-v", dest="verbose", default=False, action="store_true",help="Toggle verbose mode - prints extra infromation to the command line while executing")
	parser.add_option("--check","-c", default=False, action="store_true",help="Checks to see if the command line arguments will work.")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type="string", action="append", help="postprocessor to be applied to each projection. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
	parser.add_option("--parallel",help="Parallelism string",default=None,type="string")

	(options, args) = parser.parse_args()
	
	if ( options.check ): options.verbose = True
	
	if len(args) < 1:
		parser.error("Error: No input file given")
	
	options.model = args[0]
	error = check(options,True)
	
	if ( options.verbose ):
		if (error):
			print "e2project3d.py command line arguments test.... FAILED"
		else:
			if (options.verbose):
				print "e2project3.py command line arguments test.... PASSED"

	# returning a different error code is currently important to e2refine.py - returning 0 tells e2refine.py that it has enough
	# information to execute this script
	if error : exit(1)
	if options.check: exit(0)
	
	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if ( os.path.exists(options.outfile )):
		if ( options.force ):
			remove_file(options.outfile)

	logger=E2init(sys.argv)
	
		
	if options.parallel:
		job = EMParallelProject3D(options,args,logger)
		job.execute()
		E2end(logger)
		exit(0)
	
	
	eulers = []
	
	data = EMData()
	data.read_image(args[0],0)
	
	sym_object = parsesym(options.sym)
	[og_name,og_args] = parsemodopt(options.orientgen)
	eulers = sym_object.gen_orientations(og_name, og_args)
		
	# generate and save all the projections to disk - that's it, that main job is done
	if ( options.verbose ):
		print "Generating and saving projections..."
	generate_and_save_projections(options, data, eulers, options.smear)
	
	if ( options.verbose ):
		print "%s...done" %progname
	
	E2end(logger)
	exit(0)
	
	#if options.verifymirror:
		#if options.sym:
			#verify_mirror_test(data, eulers, options.sym, options.projector)
		#else:
			#print "Warning: verify mirror only works when a symmetry has been specified. No action taken."

#

def generate_and_save_projections(options, data, eulers, smear=0):
	for i,euler in enumerate(eulers):
		p=data.project(options.projector,euler)
		p.set_attr("xform.projection",euler)
		p.set_attr("ptcl_repr",0)
		
		if smear:
			pass
			#smear_iterator = deg2rad*euler[2] + deg2rad*phiprop/(smear+1)
			#while ( smear_iterator < euler[2] + phiprop*deg2rad ):
				#ptmp=data.project(options.projector,{"alt" : euler[0],"az" : euler[1],"phi" : smear_iterator * rad2deg})
				
				#p.add(ptmp)
				#smear_iterator +=  deg2rad*phiprop/(smear+1)

		if options.postprocess != None:
			for proc in options.postprocess:
				try:
					(processorname, param_dict) = parsemodopt(proc)
					if not param_dict : param_dict={}
					p.process_inplace(str(processorname), param_dict)
				except:
					print "warning - application of the post processor",p," failed. Continuing anyway"

		try: 
			p.write_image(options.outfile,-1)
		except:
			print "Error: Cannot write to file %s"%options.outfile
			exit(1)
		
		if (options.verbose):
			d = euler.get_params("eman")
			print "%d\t%4.2f\t%4.2f\t%4.2f" % (i, d["az"], d["alt"], d["phi"])


#def verify_mirror_test(data, eulers, symmetry, projector):
	
	#sym_object = get_sym_object( symmetry )
	
	#for i,euler in enumerate(eulers) :
		#a = {"alt" : euler[0] * rad2deg,"az" : euler[1] * rad2deg,"phi" : euler[2] * rad2deg}
		#t3d = Transform3D(EULER_EMAN, a)
		#b = {"t3d": t3d }
		#p=data.project(options.projector,b)
		
		#mirrorEuler = sym_object.asym_unit_mirror_orientation(euler)
			
		
		##Get the projection in this orientations
		#p_mirror = data.project(projector,{"alt" : mirrorEuler[0] * rad2deg,"az" : mirrorEuler[1]* rad2deg,"phi" : mirrorEuler[2] * rad2deg})
		
		### Actually do the mirroring
		#if sym_object.is_d_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'y'})
		#elif sym_object.is_c_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'x'})
		##FIXME: The mirror orientation is dependent on the platonic symmetry see http://blake.bcm.edu/emanwiki/EMAN2/Symmetry
		#elif sym_object.is_platonic_symmetry():
			#p_mirror.process_inplace("mirror", {"axis":'y'})
			
		
		### Calculate and write the difference to disk
		#p_difference = p_mirror-p
		#p_difference.write_image(symmetry+"mirror_debug_difference.img",-1)
		#p_mirror.write_image(symmetry+"mirror_debug.img",-1)
		
		### Print debug information
		#print "Orientation %d\t%4.2f\t%4.2f\t%4.2f" % (i, euler[0] * rad2deg, euler[1] * rad2deg, euler[2] * rad2deg)
		#print "Mirror %d\t%4.2f\t%4.2f\t%4.2f" % (i, mirrorEuler[0] * rad2deg, mirrorEuler[1] * rad2deg, mirrorEuler[2] * rad2deg)

def check(options, verbose=False):
	
	error = False
	
	if ( not options.sym ):
		if verbose:
			print "Error: you must specify the sym argument"
		error = True
	else:
		try: sym = parsesym(options.sym)
		except Exception, inst:
			if ( verbose ):
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args:
			error = True
	
	if ( not options.orientgen ):
		if verbose:
			print "Error: you must specify the orientgen argument"
		error = True
	elif ( check_eman2_type(options.orientgen,OrientGens,"Orientgen") == False ):
		error = True
	
	if ( check_eman2_type(options.projector,Projectors,"Projector") == False ):
		error = True

	if not os.path.exists(options.model) and not db_check_dict(options.model):
		if verbose:
			print "Error: 3D image %s does not exist" %options.model
		error = True
	
	if ( options.force and options.append):
		if verbose:
			print "Error: cannot specify both append and force"
		error = True
		
	if ( options.nofilecheck == False and os.path.exists(options.outfile )):
		if ( not options.force ):
			if verbose:
				print "Error: output file exists, use -f to overwrite or -a to append. No action taken"
			error = True
			
	if hasattr(options,"parallel") and options.parallel != None:
  		if len(options.parallel) < 2:
  			print "The parallel option %s does not make sense" %options.parallel
  			error = True
  		elif options.parallel[:2] != "dc":
  			print "Only dc parallelism is currently supported"
  			error = True
  		elif len(options.parallel.split(":")) != 3:
  			print "dc parallel options must be formatted like 'dc:localhost:9990'"
  			error = True
	
	return error

if __name__=="__main__":
	main()
