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

import sys, math, os, random
from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi
DEBUG = False
WEN_JIANG = False
EMAN1_OCT = False
MIRROR_DEBUG = True
NO_MIRROR = False


class EMParallelProject3D:
	def __init__(self,options,fsp,sym,start,modeln=0,logger=None):
		'''
		@param options the options produced by (options, args) = parser.parse_args()
		@param args the options produced by (options, args) = parser.parse_args()
		@param logger and EMAN2 logger, i.e. logger=E2init(sys.argv)
		assumes you have already called the check function.
		'''
		self.options = options
		self.args = fsp
		self.sym=sym
		self.logger = logger
		self.start=start
		self.modeln=modeln

		from EMAN2PAR import EMTaskCustomer
		self.etc=EMTaskCustomer(options.parallel)
		print "Precache ",fsp
		self.etc.precache([fsp])

		self.num_cpus = self.etc.cpu_est()
		print self.num_cpus," total CPUs available"
		if self.num_cpus > 64: # upper limit
			self.num_cpus = 64

		self.__task_options = None

	def __init_memory(self,options):
		'''

		'''
		sym_object = parsesym(self.sym)
		[og_name,og_args] = parsemodopt(options.orientgen)
		self.eulers = sym_object.gen_orientations(og_name, og_args)

	def __get_task_options(self,options):

		if self.__task_options == None:
			d = {}
			d["projector"] = parsemodopt(options.projector)
			d["prethreshold"] = options.prethreshold
			self.__task_options = d

		return self.__task_options

	def execute(self):
#		from EMAN2PAR import EMTaskCustomer

		if len(self.options.parallel) > 1:
			self.__init_memory(self.options)

			num_tasks = self.num_cpus
			# In the worst case we can only spawn as many tasks as there are eulers
			if self.num_cpus > len(self.eulers): num_tasks = len(self.eulers)

			eulers_per_task = len(self.eulers)/num_tasks
			resid_eulers = len(self.eulers) - eulers_per_task*num_tasks # we can distribute the residual evenly

			first = 0
			task_customers = []
			tids = []
#			self.etc=EMTaskCustomer(self.options.parallel)
			for i in xrange(0,num_tasks):
				last = first+eulers_per_task
				if resid_eulers > 0:
					last +=1
					resid_eulers -= 1

				tmp_eulers = self.eulers[first:last]
				indices = range(first,last)

				data = {}
				data["input"] = ("cache",self.args,0)
				data["eulers"] = tmp_eulers
				data["indices"] = indices

				task = EMProject3DTaskDC(data=data,options=self.__get_task_options(self.options))

				#print "Est %d CPUs"%etc.cpu_est()
				tid=self.etc.send_task(task)
					#print "Task submitted tid=",tid

				tids.append(tid)

				first = last

			print "Task ids are", tids

			while 1:

				print len(tids),"projection tasks left in main loop"
				st_vals = self.etc.check_task(tids)
				for i in xrange(len(tids)-1,-1,-1):
					st = st_vals[i]
					if st==100:
						tid = tids[i]

						rslts = self.etc.get_results(tid)

						if not self.__write_output_data(rslts[1]):
							print "There was a problem with the task of id",tid

						if self.logger != None:
							E2progress(self.logger,1.0-len(tids)/float(num_tasks))
							if self.options.verbose>0:
								print "%d/%d\r"%(num_tasks-len(tids),num_tasks)
								sys.stdout.flush()

						print "Task",tids.pop(i),"completed"
						print "These tasks are remaining:",tids

				if len(tids) == 0: break
				time.sleep(5)

			return len(self.eulers)
		else:
			raise NotImplementedError("The parallelism option you specified (%s) is not suppored" %self.options.parallel )

	def __write_output_data(self,rslts):
		for idx,image in rslts.items():
			if not isinstance(image,EMData): continue # this is here because we get the dimensions of the database as a key (e.g. '40x40x1').
			image["model_id"]=self.modeln
			if self.options.append : image.write_image(self.options.outfile,-1)
			else : image.write_image(self.options.outfile,idx+self.start)

		return True

def prethreshold(img):
	"""Applies an automatic threshold to the image"""
	snz=img["sigma_nonzero"]
	img.process_inplace("threshold.belowtozero",{"minval":snz*1.5})
	img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
	img.process_inplace("threshold.belowtozero",{"minval":snz/100.0})

class EMProject3DTaskDC(JSTask):
	def __init__(self,command="e2project3d.py",data=None,options=None):
		JSTask.__init__(self,command,data,options)

		# data has these keys:
		# input - which is the name of the threed model - a Task-style cache
		# eulers - a list of Transforms with which to generate projections
		# indices - indices that correspond to the ordering, has to be the same length as eulers. Returned to calling routine - use to write output in order. Not sure about this.

		# options has these keys
		# projector - [string,dict] convention

		self.projections = {} # key will be index, value will be projection

	def execute(self,progress_callback):
		progress_callback(0)
		input_name=self.data["input"][1]
		threed_image = EMData(input_name,self.data["input"][2]) # so the idx should always be 0, but it doesn't have to be

		options = self.options
		if options["prethreshold"] : prethreshold(threed_image)
		eulers = self.data["eulers"]
		indices = self.data["indices"]
		projector,projector_opts = options["projector"][0],options["projector"][1]
		projections = {}
		n = len(eulers)
		progress_callback(1)
		for i in xrange(n):
			euler = eulers[i]
			projector_opts["transform"] = euler
			projection = threed_image.project(projector,projector_opts)
			# The 5.0 is arbitrary. The goal is to get sigma in the ~1-3 range, and with typical density patterns, this should get in the right neighborhood
			projection.mult(5.0/projection["nx"])		
			projection.set_attr("xform.projection",euler)
			projection.set_attr("ptcl_repr",0)
			projections[indices[i]] = projection
#			print "call back",int(100*(i+1)/float(n))
			progress_callback(int(100*(i+1)/float(n)))
		#d = {}
		#d["projections"] = projections
		return projections

jsonclasses["EMProject3DTaskDC"]=EMProject3DTaskDC.from_jsondict

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <imagefile> [<imagefile>] [imagefile] [options]
	Generates 2-D projections of a 3-D object or multiple 3-D objects. Various options for specifiying distribution of orientations
	and other options.

	Typical usage:
	e2project3d.py map.mrc --outfile=projections.hdf --orientgen=eman:delta=5 --sym=c3 --projector=standard --verbose=2"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. If multiple input models are specified, multiple comma-separated symmetries may also be specified.",default="c1")
	parser.add_argument("--orientgen", dest="orientgen", help="The orientation generator to use. See e2help.py orientgen. Example: --orientgen=eman:delta=3.0:inc_mirror=0:perturb=1")
	parser.add_argument("--outfile", dest = "outfile", default = "e2proj.hdf", help = "Output file. Default is 'e2proj.img'")
	# add --perturb
	parser.add_argument("--smear", dest = "smear", type = int, default=0,help="Used in conjunction with --phitoo, this will rotationally smear between phi steps. The user must specify the amount of smearing (typically 2-10)")
	parser.add_argument("--projector", dest = "projector", default = "standard",help = "Projector to use")
	#parser.add_argument("--verifymirror",action="store_true",help="Used for testing the accuracy of mirror projects",default=False)
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_argument("--append", "-a",dest="append",default=False, action="store_true",help="Append to the output file")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--check","-c", default=False, action="store_true",help="Checks to see if the command line arguments will work.")
	parser.add_argument("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_argument("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type=str, action="append", help="postprocessor to be applied to each projection. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
	parser.add_argument("--cuda",action="store_true", help="Use CUDA for the projections.",default=False)
	parser.add_argument("--prethreshold",action="store_true", help="Applies an automatic threshold to the volume before projecting",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--parallel",help="Parallelism string",default=None,type=str)

	(options, args) = parser.parse_args()

	# Start autoflushing output
	#autoflush()

	if ( options.check ): options.verbose = 9

	print "project3d: ",str(options)

	if len(args) < 1:
		parser.error("Error: No input file given")

	# the model is a list of models. We also make sure we have the same number of symmetry elements as models
	options.model = args
	options.sym=options.sym.split(",")
	if len(options.sym)!=1:
		if len(options.sym)!=len(args) :
			print "sym must be either a single symmetry specifier or one specifier for each input."
			sys.exit(1)
	else:
		options.sym*=len(args)

	error = check(options,True)

	if ( options.verbose>0 ):
		if (error):
			print "e2project3d.py command line arguments test.... FAILED"
		else:
			if (options.verbose>0):
				print "e2project3.py command line arguments test.... PASSED"

	# returning a different error code is currently important to e2refine.py - returning 0 tells e2refine.py that it has enough
	# information to execute this script
	if error : exit(1)
	if options.check: exit(0)

	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if ( os.path.exists(options.outfile )):
		if ( options.force ):
			remove_file(options.outfile)

	logger=E2init(sys.argv,options.ppid)


	if options.parallel:
		try:
			p=options.parallel.split(":")
			if p[0]=="mpi" and int(p[1])>64:
				print "Modified parallelism in projection to use 64 cores"
				options.parallel="%s:64:%s"%(p[0],p[2])
		except: pass

		n=0
		for i,fsp in enumerate(args) :
			job = EMParallelProject3D(options,fsp,options.sym[i],n,i+1,logger)
			n+=job.execute()
			if options.verbose : print "Job %d finished. %d total projections."%(i,n)

		E2end(logger)
		exit(0)


	eulers = []
	if options.cuda: EMData.switchoncuda()
	for i,fsp in enumerate(args) :
		data = EMData(fsp,0)
		if options.prethreshold : prethreshold(data)

		sym_object = parsesym(options.sym[i])
		[og_name,og_args] = parsemodopt(options.orientgen)
		eulers = sym_object.gen_orientations(og_name, og_args)

		# generate and save all the projections to disk - that's it, that main job is done
		if ( options.verbose>0 ):
			print "Generating and saving projections for ",fsp
		generate_and_save_projections(options, data, eulers, options.smear,i+1)
	if options.cuda: EMData.switchoffcuda()

	if ( options.verbose>0 ):
		print "%s...done" %progname

	E2end(logger)
	exit(0)

	#if options.verifymirror:
		#if options.sym:
			#verify_mirror_test(data, eulers, options.sym, options.projector)
		#else:
			#print "Warning: verify mirror only works when a symmetry has been specified. No action taken."

#

def generate_and_save_projections(options, data, eulers, smear=0,modeln=0):
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

		p["model_id"]=modeln
		try:
			if options.append: p.write_image(options.outfile,-1)
			else : p.write_image(options.outfile,i)
		except:
			print "Error: Cannot write to file %s"%options.outfile
			exit(1)

		if (options.verbose>0):
			d = euler.get_params("eman")
			print "%d\t%4.2f\t%4.2f\t%4.2f" % (i, d["az"], d["alt"], d["phi"])

def check(options, verbose=0):

	error = False

	for s in options.sym:
		try: sym = parsesym(s)
		except Exception, inst:
			if verbose>0:
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args:
			error = True

	if ( not options.orientgen ):
		if verbose>0:
			print "Error: you must specify the orientgen argument"
		error = True
	elif ( check_eman2_type(options.orientgen,OrientGens,"Orientgen") == False ):
		error = True

	if ( check_eman2_type(options.projector,Projectors,"Projector") == False ):
		error = True

	for f in options.model:
		if not os.path.exists(f) and not db_check_dict(f):
			if verbose>0:
				print "Error: 3D image %s does not exist" %f
			error = True

	if ( options.force and options.append):
		if verbose>0:
			print "Error: cannot specify both append and force"
		error = True

	if ( options.nofilecheck == False and os.path.exists(options.outfile )):
		if ( not options.force and not options.append):
			if verbose>0:
				print "Error: output file exists, use -f to overwrite or -a to append. No action taken"
			error = True

	if hasattr(options,"parallel") and options.parallel != None:
  		if len(options.parallel) < 2:
  			print "The parallel option %s does not make sense" %options.parallel
  			error = True

	return error

if __name__=="__main__":
	main()
