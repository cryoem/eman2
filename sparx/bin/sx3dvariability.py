#!/usr/bin/env python
from __future__ import print_function
#
# Author: 
# Copyright (c) 2012 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
from EMAN2 import *
from sparx import *
from global_def import SPARX_MPI_TAG_UNIVERSAL
import	global_def
from	global_def 	import *
from	optparse 	import OptionParser
from	EMAN2 		import EMUtil
import	os
import	sys
from 	time		import	time
	
def main():

	def params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror):
		# the final ali2d parameters already combine shifts operation first and rotation operation second for parameters converted from 3D
		if mirror:
			m = 1
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 540.0-psi, 0, 0, 1.0)
		else:
			m = 0
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 360.0-psi, 0, 0, 1.0)
		return  alpha, sx, sy, m
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack  --ave2D= --var2D=  --ave3D= --var3D= --img_per_grp= --fl=15. --aa=0.01  --sym=symmetry --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)
	
	parser.add_option("--output_dir",   type="string"	   ,	default="./",				help="output directory")
	parser.add_option("--ave2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D averages")
	parser.add_option("--var2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D variances")
	parser.add_option("--ave3D",		type="string"	   ,	default=False,				help="write to the disk reconstructed 3D average")
	parser.add_option("--var3D",		type="string"	   ,	default=False,				help="compute 3D variability (time consuming!)")
	parser.add_option("--img_per_grp",	type="int"         ,	default=10   ,				help="number of neighbouring projections")
	parser.add_option("--no_norm",		action="store_true",	default=False,				help="do not use normalization")
	#parser.add_option("--radius", 	    type="int"         ,	default=-1   ,				help="radius for 3D variability" )
	parser.add_option("--npad",			type="int"         ,	default=2    ,				help="number of time to pad the original images")
	parser.add_option("--sym" , 		type="string"      ,	default="c1" ,				help="symmetry")
	parser.add_option("--fl",			type="float"       ,	default=0.0  ,				help="low pass filter cutoff in absolute frequency (0.0-0.5) and is applied after decimation. (Default - no filtration)")
	parser.add_option("--aa",			type="float"       ,	default=0.01 ,				help="fall off of the filter. Put 0.01 if user has no clue about falloff (Default - no filtration)")
	parser.add_option("--CTF",			action="store_true",	default=False,				help="use CFT correction")
	parser.add_option("--VERBOSE",		action="store_true",	default=False,				help="Long output for debugging")
	#parser.add_option("--MPI" , 		action="store_true",	default=False,				help="use MPI version")
	#parser.add_option("--radiuspca", 	type="int"         ,	default=-1   ,				help="radius for PCA" )
	#parser.add_option("--iter", 		type="int"         ,	default=40   ,				help="maximum number of iterations (stop criterion of reconstruction process)" )
	#parser.add_option("--abs", 		type="float"   ,        default=0.0  ,				help="minimum average absolute change of voxels' values (stop criterion of reconstruction process)" )
	#parser.add_option("--squ", 		type="float"   ,	    default=0.0  ,				help="minimum average squared change of voxels' values (stop criterion of reconstruction process)" )
	parser.add_option("--VAR" , 		action="store_true",	default=False,				help="stack of input consists of 2D variances (Default False)")
	parser.add_option("--decimate",     type="float",           default= 1.0,               help="image decimate rate, a number less than 1. default is 1")
	parser.add_option("--window",       type="int",             default= 0,                 help="target image size. Default value is zero.")
	#parser.add_option("--SND",			action="store_true",	default=False,				help="compute squared normalized differences (Default False)")
	parser.add_option("--nvec",			type="int"         ,	default=0    ,				help="number of eigenvectors, default = 0 meaning no PCA calculated")
	parser.add_option("--symmetrize",	action="store_true",	default=False,				help="Prepare input stack for handling symmetry (Default False)")
	
	(options,args) = parser.parse_args()
	#####
	from mpi import mpi_init, mpi_comm_rank, mpi_comm_size, mpi_recv, MPI_COMM_WORLD
	from mpi import mpi_barrier, mpi_reduce, mpi_bcast, mpi_send, MPI_FLOAT, MPI_SUM, MPI_INT, MPI_MAX
	from applications   import MPI_start_end
	from reconstruction import recons3d_em, recons3d_em_MPI
	from reconstruction	import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities import print_begin_msg, print_end_msg, print_msg
	from utilities import read_text_row, get_image, get_im
	from utilities import bcast_EMData_to_all, bcast_number_to_all
	from utilities import get_symt

	#  This is code for handling symmetries by the above program.  To be incorporated. PAP 01/27/2015

	from EMAN2db import db_open_dict

	# Set up global variables related to bdb cache 
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	# Set up global variables related to ERROR function
	global_def.BATCH = True
	
	# detect if program is running under MPI
	RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
	if RUNNING_UNDER_MPI: global_def.MPI = True
	
	if options.symmetrize :
		if RUNNING_UNDER_MPI:
			try:
				sys.argv = mpi_init(len(sys.argv), sys.argv)
				try:	
					number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
					if( number_of_proc > 1 ):
						ERROR("Cannot use more than one CPU for symmetry prepration","sx3dvariability",1)
				except:
					pass
			except:
				pass
		if options.output_dir !="./" and not os.path.exists(options.output_dir): os.mkdir(options.output_dir)
		#  Input
		#instack = "Clean_NORM_CTF_start_wparams.hdf"
		#instack = "bdb:data"
		
		
		from logger import Logger,BaseLogger_Files
		if os.path.exists(os.path.join(options.output_dir, "log.txt")): os.remove(os.path.join(options.output_dir, "log.txt"))
		log_main=Logger(BaseLogger_Files())
		log_main.prefix = os.path.join(options.output_dir, "./")
		
		instack = args[0]
		sym = options.sym.lower()
		if( sym == "c1" ):
			ERROR("There is no need to symmetrize stack for C1 symmetry","sx3dvariability",1)
		
		line =""
		for a in sys.argv:
			line +=" "+a
		log_main.add(line)
	
		if(instack[:4] !="bdb:"):
			if output_dir =="./": stack = "bdb:data"
			else: stack = "bdb:"+options.output_dir+"/data"
			delete_bdb(stack)
			junk = cmdexecute("sxcpy.py  "+instack+"  "+stack)
		else: stack = instack
		
		qt = EMUtil.get_all_attributes(stack,'xform.projection')

		na = len(qt)
		ts = get_symt(sym)
		ks = len(ts)
		angsa = [None]*na
		
		for k in xrange(ks):
			#Qfile = "Q%1d"%k
			if options.output_dir!="./": Qfile = os.path.join(options.output_dir,"Q%1d"%k)
			else: Qfile = os.path.join(options.output_dir, "Q%1d"%k)
			#delete_bdb("bdb:Q%1d"%k)
			delete_bdb("bdb:"+Qfile)
			#junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
			junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:"+Qfile)
			#DB = db_open_dict("bdb:Q%1d"%k)
			DB = db_open_dict("bdb:"+Qfile)
			for i in xrange(na):
				ut = qt[i]*ts[k]
				DB.set_attr(i, "xform.projection", ut)
				#bt = ut.get_params("spider")
				#angsa[i] = [round(bt["phi"],3)%360.0, round(bt["theta"],3)%360.0, bt["psi"], -bt["tx"], -bt["ty"]]
			#write_text_row(angsa, 'ptsma%1d.txt'%k)
			#junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
			#junk = cmdexecute("sxheader.py  bdb:Q%1d  --params=xform.projection  --import=ptsma%1d.txt"%(k,k))
			DB.close()
		if options.output_dir =="./": delete_bdb("bdb:sdata")
		else: delete_bdb("bdb:" + options.output_dir + "/"+"sdata")
		#junk = cmdexecute("e2bdb.py . --makevstack=bdb:sdata --filt=Q")
		sdata = "bdb:"+options.output_dir+"/"+"sdata"
		print(sdata)
		junk = cmdexecute("e2bdb.py   " + options.output_dir +"  --makevstack="+sdata +" --filt=Q")
		#junk = cmdexecute("ls  EMAN2DB/sdata*")
		#a = get_im("bdb:sdata")
		a = get_im(sdata)
		a.set_attr("variabilitysymmetry",sym)
		#a.write_image("bdb:sdata")
		a.write_image(sdata)

	else:

		from fundamentals import window2d
		sys.argv       = mpi_init(len(sys.argv), sys.argv)
		myid           = mpi_comm_rank(MPI_COMM_WORLD)
		number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
		main_node      = 0
		keepgoing      = 1
		if len(args) == 1: stack = args[0]
		else:
			print(( "usage: " + usage))
			print(( "Please run '" + progname + " -h' for detailed options"))
			return 1

		t0 = time()	
		# obsolete flags
		options.MPI  = True
		options.nvec = 0
		options.radiuspca = -1
		options.iter = 40
		options.abs  = 0.0
		options.squ  = 0.0

		if options.fl > 0.0 and options.aa == 0.0:
			ERROR("Fall off has to be given for the low-pass filter", "sx3dvariability", 1, myid)
			
		#if options.VAR and options.SND:
		#	ERROR("Only one of var and SND can be set!", "sx3dvariability", myid)
			
		if options.VAR and (options.ave2D or options.ave3D or options.var2D): 
			ERROR("When VAR is set, the program cannot output ave2D, ave3D or var2D", "sx3dvariability", 1, myid)
			
		#if options.SND and (options.ave2D or options.ave3D):
		#	ERROR("When SND is set, the program cannot output ave2D or ave3D", "sx3dvariability", 1, myid)
		
		if options.nvec > 0 :
			ERROR("PCA option not implemented", "sx3dvariability", 1, myid)
			
		if options.nvec > 0 and options.ave3D == None:
			ERROR("When doing PCA analysis, one must set ave3D", "sx3dvariability", 1, myid)
		
		if options.decimate>1.0 or options.decimate<0.0:
			ERROR("Decimate rate should be a value between 0.0 and 1.0", "sx3dvariability", 1, myid)
		
		if options.window < 0.0:
			ERROR("Target window size should be always larger than zero", "sx3dvariability", 1, myid)
			
		if myid == main_node:
			img  = get_image(stack, 0)
			nx   = img.get_xsize()
			ny   = img.get_ysize()
			if int(options.decimate*min(nx, ny)+0.5)< options.window:keepgoing = 0
		keepgoing = bcast_number_to_all(keepgoing, main_node, MPI_COMM_WORLD)
		if keepgoing == 0: ERROR("The target window size cannot be larger than the size of decimated image", "sx3dvariability", 1, myid)

		import string
		options.sym = options.sym.lower()
		 
		# if global_def.CACHE_DISABLE:
		# 	from utilities import disable_bdb_cache
		# 	disable_bdb_cache()
		# global_def.BATCH = True
		
		if myid == main_node:
			if options.output_dir !="./" and not os.path.exists(options.output_dir): os.mkdir(options.output_dir)
	
		img_per_grp = options.img_per_grp
		nvec        = options.nvec
		radiuspca   = options.radiuspca
		
		if myid == main_node:
		
			from logger import Logger,BaseLogger_Files
			#if os.path.exists(os.path.join(options.output_dir, "log.txt")): os.remove(os.path.join(options.output_dir, "log.txt"))
			log_main=Logger(BaseLogger_Files())
			log_main.prefix = os.path.join(options.output_dir, "./")
			line = ""
			for a in sys.argv: line +=" "+a
			log_main.add(line)
			log_main.add("-------->>>Settings given by all options<<<-------")
			log_main.add("Symmetry             : %s"%options.sym)
			log_main.add("Input stack          : %s"%stack)
			log_main.add("Output_dir           : %s"%options.output_dir)
			
			if options.ave3D: log_main.add("Ave3d                : %s"%options.ave3D)
			if options.var3D: log_main.add("Var3d                : %s"%options.var3D)
			if options.ave2D: log_main.add("Ave2D                : %s"%options.ave2D)
			if options.var2D: log_main.add("Var2D                : %s"%options.var2D)
			if options.VAR:   log_main.add("VAR                  : True")
			else:             log_main.add("VAR                  : False")
			if options.CTF:   log_main.add("CTF correction       : True  ")
			else:             log_main.add("CTF correction       : False ")
			
			log_main.add("Image per group      : %5d"%options.img_per_grp)
			log_main.add("Image decimate rate  : %4.3f"%options.decimate)
			log_main.add("Low pass filter      : %4.3f"%options.fl)
			log_main.add("Current low pass filter is equivalent to cutoff frequency %4.3f for original image size"%round((options.fl*options.decimate),3))
			log_main.add("Window size          : %5d "%options.window)
			log_main.add("sx3dvariability begins")
	
		symbaselen = 0
		if myid == main_node:
			nima = EMUtil.get_image_count(stack)
			img  = get_image(stack)
			nx   = img.get_xsize()
			ny   = img.get_ysize()
			org_nx = nx
			org_ny = nx
			org_nz = nx
			if options.sym != "c1" :
				imgdata = get_im(stack)
				try:
					i = imgdata.get_attr("variabilitysymmetry").lower()
					if(i != options.sym):
						ERROR("The symmetry provided does not agree with the symmetry of the input stack", "sx3dvariability", 1, myid)
				except:
					ERROR("Input stack is not prepared for symmetry, please follow instructions", "sx3dvariability", 1, myid)
				from utilities import get_symt
				i = len(get_symt(options.sym))
				if((nima/i)*i != nima):
					ERROR("The length of the input stack is incorrect for symmetry processing", "sx3dvariability", 1, myid)
				symbaselen = nima/i
			else:  symbaselen = nima
		else:
			nima = 0
			nx = 0
			ny = 0
		nima    = bcast_number_to_all(nima)
		nx      = bcast_number_to_all(nx)
		ny      = bcast_number_to_all(ny)
		Tracker ={}
		Tracker["total_stack"] = nima
		
		if options.decimate==1.:
			if options.window !=0:
				nx = options.window
				ny = options.window
		else:
			if options.window == 0:
				nx = int(nx*options.decimate+0.5)
				ny = int(ny*options.decimate+0.5)
			else:
				nx = int(options.window)
				ny = nx
		Tracker["nx"]  = nx
		Tracker["ny"]  = ny
		Tracker["nz"]  = nx
		
		symbaselen     = bcast_number_to_all(symbaselen)
		if radiuspca == -1: radiuspca = nx/2-2

		if myid == main_node: log_main.add("%-70s:  %d\n"%("Number of projection", nima))
		img_begin, img_end = MPI_start_end(nima, number_of_proc, myid)
		
		"""
		if options.SND:
			from projection		import prep_vol, prgs
			from statistics		import im_diff
			from utilities		import get_im, model_circle, get_params_proj, set_params_proj
			from utilities		import get_ctf, generate_ctf
			from filter			import filt_ctf
		
			imgdata = EMData.read_images(stack, range(img_begin, img_end))

			if options.CTF:
				vol = recons3d_4nn_ctf_MPI(myid, imgdata, 1.0, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			else:
				vol = recons3d_4nn_MPI(myid, imgdata, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)

			bcast_EMData_to_all(vol, myid)
			volft, kb = prep_vol(vol)

			mask = model_circle(nx/2-2, nx, ny)
			varList = []
			for i in xrange(img_begin, img_end):
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata[i-img_begin])
				ref_prj = prgs(volft, kb, [phi, theta, psi, -s2x, -s2y])
				if options.CTF:
					ctf_params = get_ctf(imgdata[i-img_begin])
					ref_prj = filt_ctf(ref_prj, generate_ctf(ctf_params))
				diff, A, B = im_diff(ref_prj, imgdata[i-img_begin], mask)
				diff2 = diff*diff
				set_params_proj(diff2, [phi, theta, psi, s2x, s2y])
				varList.append(diff2)
			mpi_barrier(MPI_COMM_WORLD)
		"""
		
		if options.VAR: # 2D variance images have no shifts
			#varList   = EMData.read_images(stack, range(img_begin, img_end))
			varList    = []
			this_image = EMData()
			for index_of_particle in xrange(img_begin,img_end):
				this_image.read_image(stack,index_of_particle)
				#varList.append(image_decimate_window_xform_ctf(this_image, options.decimate, options.window,options.CTF))
				varList.append(window2d(resample(this_image, options.decimate), nx, ny))
				
		else:
			from utilities		import bcast_number_to_all, bcast_list_to_all, send_EMData, recv_EMData
			from utilities		import set_params_proj, get_params_proj, params_3D_2D, get_params2D, set_params2D, compose_transform2
			from utilities		import model_blank, nearest_proj, model_circle
			from applications	import pca
			from statistics		import avgvar, avgvar_ctf, ccc
			from filter		    import filt_tanl
			from morphology		import threshold, square_root
			from projection 	import project, prep_vol, prgs
			from sets		    import Set

			if myid == main_node:
				t1          = time()
				proj_angles = []
				aveList     = []
				tab = EMUtil.get_all_attributes(stack, 'xform.projection')
				
				for i in xrange(nima):
					t     = tab[i].get_params('spider')
					phi   = t['phi']
					theta = t['theta']
					psi   = t['psi']
					x     = theta
					if x > 90.0: x = 180.0 - x
					x = x*10000+psi
					proj_angles.append([x, t['phi'], t['theta'], t['psi'], i])
				t2 = time()
				log_main.add( "%-70s:  %d\n"%("Number of neighboring projections", img_per_grp))
				log_main.add("...... Finding neighboring projections\n")
				if options.VERBOSE:
					log_main.add( "Number of images per group: %d"%img_per_grp)
					log_main.add( "Now grouping projections")
				proj_angles.sort()
			proj_angles_list = [0.0]*(nima*4)
			if myid == main_node:
				for i in xrange(nima):
					proj_angles_list[i*4]   = proj_angles[i][1]
					proj_angles_list[i*4+1] = proj_angles[i][2]
					proj_angles_list[i*4+2] = proj_angles[i][3]
					proj_angles_list[i*4+3] = proj_angles[i][4]
			proj_angles_list = bcast_list_to_all(proj_angles_list, myid, main_node)
			proj_angles      = []
			for i in xrange(nima):
				proj_angles.append([proj_angles_list[i*4], proj_angles_list[i*4+1], proj_angles_list[i*4+2], int(proj_angles_list[i*4+3])])
			del proj_angles_list
			proj_list, mirror_list = nearest_proj(proj_angles, img_per_grp, range(img_begin, img_end))

			all_proj = Set()
			for im in proj_list:
				for jm in im:
					all_proj.add(proj_angles[jm][3])

			all_proj = list(all_proj)
			
			if options.VERBOSE: # all nodes info
				print("On node %2d, number of images needed to be read = %5d"%(myid, len(all_proj)))

			index = {}
			for i in xrange(len(all_proj)): index[all_proj[i]] = i
			mpi_barrier(MPI_COMM_WORLD)

			if myid == main_node:
				log_main.add("%-70s:  %.2f\n"%("Finding neighboring projections lasted [s]", time()-t2))
				log_main.add("%-70s:  %d\n"%("Number of groups processed on the main node", len(proj_list)))
				if options.VERBOSE:
					log_main.add("Grouping projections took: ", (time()-t2)/60	, "[min]")
					log_main.add("Number of groups on main node: ", len(proj_list))
			mpi_barrier(MPI_COMM_WORLD)

			if myid == main_node:
				log_main.add("...... Calculating the stack of 2D variances \n")
				if options.VERBOSE:
					log_main.add("Now calculating the stack of 2D variances")

			proj_params = [0.0]*(nima*5)
			aveList = []
			varList = []				
			if nvec > 0:
				eigList = [[] for i in xrange(nvec)]

			if options.VERBOSE: 	print("Begin to read images on processor %d"%(myid)) # all nodes info
			ttt = time()
			#imgdata = EMData.read_images(stack, all_proj)
			imgdata = []
			
			#if myid==0: print(get_im(stack, 0).get_attr_dict())
			for index_of_proj in xrange(len(all_proj)):
				#img     = EMData()
				#img.read_image(stack, all_proj[index_of_proj])
				#dmg = image_decimate_window_xform_ctf(get_im(stack, all_proj[index_of_proj]), options.decimate, options.window, options.CTF)
				#print dmg.get_xsize(), "init"
				imgdata.append(window2d(resample(get_im(stack, all_proj[index_of_proj]), options.decimate), nx, ny))
			
			#if myid ==0 : print(imgdata[0].get_attr_dict())
			if options.VERBOSE: # all nodes info
				print("Reading images on processor %d done, time = %.2f"%(myid, time()-ttt))
				print("On processor %d, we got %d images"%(myid, len(imgdata)))
			mpi_barrier(MPI_COMM_WORLD)

			'''	
			imgdata2 = EMData.read_images(stack, range(img_begin, img_end))
			if options.fl > 0.0:
				for k in xrange(len(imgdata2)):
					imgdata2[k] = filt_tanl(imgdata2[k], options.fl, options.aa)
			if options.CTF:
				vol = recons3d_4nn_ctf_MPI(myid, imgdata2, 1.0, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			else:
				vol = recons3d_4nn_MPI(myid, imgdata2, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			if myid == main_node:
				vol.write_image("vol_ctf.hdf")
				print_msg("Writing to the disk volume reconstructed from averages as		:  %s\n"%("vol_ctf.hdf"))
			del vol, imgdata2
			mpi_barrier(MPI_COMM_WORLD)
			'''
			from applications import prepare_2d_forPCA
			from utilities    import model_blank
			for i in xrange(len(proj_list)):
				ki = proj_angles[proj_list[i][0]][3]
				if ki >= symbaselen:  continue
				mi = index[ki]
				phiM, thetaM, psiM, s2xM, s2yM = get_params_proj(imgdata[mi])

				grp_imgdata = []
				for j in xrange(img_per_grp):
					mj = index[proj_angles[proj_list[i][j]][3]]
					phi, theta, psi, s2x, s2y = get_params_proj(imgdata[mj])
					alpha, sx, sy, mirror = params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror_list[i][j])
					if thetaM <= 90:
						if mirror == 0:  alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, phiM-phi, 0.0, 0.0, 1.0)
						else:            alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, 180-(phiM-phi), 0.0, 0.0, 1.0)
					else:
						if mirror == 0:  alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, -(phiM-phi), 0.0, 0.0, 1.0)
						else:            alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, -(180-(phiM-phi)), 0.0, 0.0, 1.0)
					set_params2D(imgdata[mj], [alpha, sx, sy, mirror, 1.0])
					grp_imgdata.append(imgdata[mj])
					#print grp_imgdata[j].get_xsize(), imgdata[mj].get_xsize()

				if not options.no_norm:
					#print grp_imgdata[j].get_xsize()
					mask = model_circle(nx/2-2, nx, nx)
					for k in xrange(img_per_grp):
						ave, std, minn, maxx = Util.infomask(grp_imgdata[k], mask, False)
						grp_imgdata[k] -= ave
						grp_imgdata[k] /= std
					del mask

				if options.fl > 0.0:
					from filter import filt_ctf, filt_table
					nx2 = 2*nx
					ny2 = 2*ny
					if options.CTF:
						from utilities import pad
						for k in xrange(img_per_grp):
							grp_imgdata[k] = window2d(fft( filt_tanl( filt_ctf(fft(pad(grp_imgdata[k], nx2, ny2, 1,0.0)), grp_imgdata[k].get_attr("ctf"), binary=1), options.fl, options.aa) ),nx,ny)
							#grp_imgdata[k] = window2d(fft( filt_table( filt_tanl( filt_ctf(fft(pad(grp_imgdata[k], nx2, ny2, 1,0.0)), grp_imgdata[k].get_attr("ctf"), binary=1), options.fl, options.aa), fifi) ),nx,ny)
							#grp_imgdata[k] = filt_tanl(grp_imgdata[k], options.fl, options.aa)
					else:
						for k in xrange(img_per_grp):
							grp_imgdata[k] = filt_tanl(grp_imgdata[k], options.fl, options.aa)
							#grp_imgdata[k] = window2d(fft( filt_table( filt_tanl( filt_ctf(fft(pad(grp_imgdata[k], nx2, ny2, 1,0.0)), grp_imgdata[k].get_attr("ctf"), binary=1), options.fl, options.aa), fifi) ),nx,ny)
							#grp_imgdata[k] = filt_tanl(grp_imgdata[k], options.fl, options.aa)
				else:
					from utilities import pad, read_text_file
					from filter    import filt_ctf, filt_table
					nx2 = 2*nx
					ny2 = 2*ny
					if options.CTF:
						from utilities import pad
						for k in xrange(img_per_grp):
							grp_imgdata[k] = window2d( fft( filt_ctf(fft(pad(grp_imgdata[k], nx2, ny2, 1,0.0)), grp_imgdata[k].get_attr("ctf"), binary=1) ) , nx,ny)
							#grp_imgdata[k] = window2d(fft( filt_table( filt_tanl( filt_ctf(fft(pad(grp_imgdata[k], nx2, ny2, 1,0.0)), grp_imgdata[k].get_attr("ctf"), binary=1), options.fl, options.aa), fifi) ),nx,ny)
							#grp_imgdata[k] = filt_tanl(grp_imgdata[k], options.fl, options.aa)

				'''
				if i < 10 and myid == main_node:
					for k in xrange(10):
						grp_imgdata[k].write_image("grp%03d.hdf"%i, k)
				'''
				"""
				if myid == main_node and i==0:
					for pp in xrange(len(grp_imgdata)):
						grp_imgdata[pp].write_image("pp.hdf", pp)
				"""
				ave, grp_imgdata = prepare_2d_forPCA(grp_imgdata)
				"""
				if myid == main_node and i==0:
					for pp in xrange(len(grp_imgdata)):
						grp_imgdata[pp].write_image("qq.hdf", pp)
				"""

				var = model_blank(nx,ny)
				for q in grp_imgdata:  Util.add_img2( var, q )
				Util.mul_scalar( var, 1.0/(len(grp_imgdata)-1))
				# Switch to std dev
				var = square_root(threshold(var))
				#if options.CTF:	ave, var = avgvar_ctf(grp_imgdata, mode="a")
				#else:	            ave, var = avgvar(grp_imgdata, mode="a")
				"""
				if myid == main_node:
					ave.write_image("avgv.hdf",i)
					var.write_image("varv.hdf",i)
				"""
			
				set_params_proj(ave, [phiM, thetaM, 0.0, 0.0, 0.0])
				set_params_proj(var, [phiM, thetaM, 0.0, 0.0, 0.0])

				aveList.append(ave)
				varList.append(var)

				if options.VERBOSE:
					print("%5.2f%% done on processor %d"%(i*100.0/len(proj_list), myid))
				if nvec > 0:
					eig = pca(input_stacks=grp_imgdata, subavg="", mask_radius=radiuspca, nvec=nvec, incore=True, shuffle=False, genbuf=True)
					for k in xrange(nvec):
						set_params_proj(eig[k], [phiM, thetaM, 0.0, 0.0, 0.0])
						eigList[k].append(eig[k])
					"""
					if myid == 0 and i == 0:
						for k in xrange(nvec):
							eig[k].write_image("eig.hdf", k)
					"""

			del imgdata
			#  To this point, all averages, variances, and eigenvectors are computed

			if options.ave2D:
				from fundamentals import fpol
				if myid == main_node:
					km = 0
					for i in xrange(number_of_proc):
						if i == main_node :
							for im in xrange(len(aveList)):
								aveList[im].write_image(os.path.join(options.output_dir, options.ave2D), km)
								km += 1
						else:
							nl = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
							nl = int(nl[0])
							for im in xrange(nl):
								ave = recv_EMData(i, im+i+70000)
								"""
								nm = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								nm = int(nm[0])
								members = mpi_recv(nm, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('members', map(int, members))
								members = mpi_recv(nm, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('pix_err', map(float, members))
								members = mpi_recv(3, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('refprojdir', map(float, members))
								"""
								tmpvol=fpol(ave, Tracker["nx"],Tracker["nx"],1)								
								tmpvol.write_image(os.path.join(options.output_dir, options.ave2D), km)
								km += 1
				else:
					mpi_send(len(aveList), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					for im in xrange(len(aveList)):
						send_EMData(aveList[im], main_node,im+myid+70000)
						"""
						members = aveList[im].get_attr('members')
						mpi_send(len(members), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						mpi_send(members, len(members), MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						members = aveList[im].get_attr('pix_err')
						mpi_send(members, len(members), MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						try:
							members = aveList[im].get_attr('refprojdir')
							mpi_send(members, 3, MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						except:
							mpi_send([-999.0,-999.0,-999.0], 3, MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						"""

			if options.ave3D:
				from fundamentals import fpol
				if options.VERBOSE and myid == main_node: log_main.add("Reconstructing 3D average volume")
				ave3D = recons3d_4nn_MPI(myid, aveList, symmetry=options.sym, npad=options.npad)
				bcast_EMData_to_all(ave3D, myid)
				if myid == main_node:
					ave3D = fpol(ave3D,Tracker["nx"],Tracker["nx"],Tracker["nx"])
					ave3D.write_image(os.path.join(options.output_dir, options.ave3D))
					log_main.add("%-70s:  %s\n"%("Writing to the disk volume reconstructed from averages as", options.ave3D))
			del ave, var, proj_list, stack, phi, theta, psi, s2x, s2y, alpha, sx, sy, mirror, aveList

			if nvec > 0:
				for k in xrange(nvec):
					if options.VERBOSE and myid == main_node:log_main.add("Reconstruction eigenvolumes", k)
					cont = True
					ITER = 0
					mask2d = model_circle(radiuspca, nx, nx)
					while cont:
						#print "On node %d, iteration %d"%(myid, ITER)
						eig3D = recons3d_4nn_MPI(myid, eigList[k], symmetry=options.sym, npad=options.npad)
						bcast_EMData_to_all(eig3D, myid, main_node)
						if options.fl > 0.0:
							eig3D = filt_tanl(eig3D, options.fl, options.aa)
						if myid == main_node:
							eig3D.write_image(os.path.join(options.outpout_dir, "eig3d_%03d.hdf"%(k, ITER)))
						Util.mul_img( eig3D, model_circle(radiuspca, nx, nx, nx) )
						eig3Df, kb = prep_vol(eig3D)
						del eig3D
						cont = False
						icont = 0
						for l in xrange(len(eigList[k])):
							phi, theta, psi, s2x, s2y = get_params_proj(eigList[k][l])
							proj = prgs(eig3Df, kb, [phi, theta, psi, s2x, s2y])
							cl = ccc(proj, eigList[k][l], mask2d)
							if cl < 0.0:
								icont += 1
								cont = True
								eigList[k][l] *= -1.0
						u = int(cont)
						u = mpi_reduce([u], 1, MPI_INT, MPI_MAX, main_node, MPI_COMM_WORLD)
						icont = mpi_reduce([icont], 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)

						if myid == main_node:
							u = int(u[0])
							log_main.add(" Eigenvector: ",k," number changed ",int(icont[0]))
						else: u = 0
						u = bcast_number_to_all(u, main_node)
						cont = bool(u)
						ITER += 1

					del eig3Df, kb
					mpi_barrier(MPI_COMM_WORLD)
				del eigList, mask2d

			if options.ave3D: del ave3D
			if options.var2D:
				from fundamentals import fpol 
				if myid == main_node:
					km = 0
					for i in xrange(number_of_proc):
						if i == main_node :
							for im in xrange(len(varList)):
								tmpvol=fpol(varList[im], Tracker["nx"], Tracker["nx"],1)
								tmpvol.write_image(os.path.join(options.output_dir, options.var2D), km)
								km += 1
						else:
							nl = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
							nl = int(nl[0])
							for im in xrange(nl):
								ave = recv_EMData(i, im+i+70000)
								tmpvol=fpol(ave, Tracker["nx"], Tracker["nx"],1)
								tmpvol.write_image(os.path.join(options.output_dir, options.var2D), km)
								km += 1
				else:
					mpi_send(len(varList), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					for im in xrange(len(varList)):
						send_EMData(varList[im], main_node, im+myid+70000)#  What with the attributes??
			mpi_barrier(MPI_COMM_WORLD)

		if  options.var3D:
			if myid == main_node and options.VERBOSE: log_main.add("Reconstructing 3D variability volume")
			t6 = time()
			# radiusvar = options.radius
			# if( radiusvar < 0 ):  radiusvar = nx//2 -3
			res = recons3d_4nn_MPI(myid, varList, symmetry=options.sym, npad=options.npad)
			#res = recons3d_em_MPI(varList, vol_stack, options.iter, radiusvar, options.abs, True, options.sym, options.squ)
			
			if myid == main_node:
				from fundamentals import fpol
				res = fpol(res, int(org_nx*options.decimate+0.5),  int(org_ny*options.decimate+0.5),  int(org_nz*options.decimate+0.5))
				res	= resample(res, 1./options.decimate)
				set_pixel_size(res, 1.0)
				res.write_image(os.path.join(options.output_dir, options.var3D))
				log_main.add("%-70s:  %.2f\n"%("Reconstructing 3D variability took [s]", time()-t6))
				log_main.add("%-70s:  %.2f\n"%("Total time for these computations [s]", time()-t0))
				log_main.add("sx3dvariability finishes")	
		from mpi import mpi_finalize
		mpi_finalize()
		
		if RUNNING_UNDER_MPI: global_def.MPI = False

		global_def.BATCH = False

if __name__=="__main__":
	main()

