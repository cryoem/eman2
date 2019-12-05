#!/usr/bin/env python
#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Please do not copy or modify this file without written consent of the author.
# Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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
#
#
from builtins import range
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
from mpi import *

"""
Instruction:
nohup mpirun -np  64   --hostfile ./node4567.txt  sx3dvariability.py bdb:data/data \
  --output_dir=var3d  --window=300 --var3D=var.hdf --img_per_grp=100 \
      --CTF>var3d/printout &
      
 1. The order of applying the input parameters is:  1. window; 2. decimation. The low-pass 
    filter is with absolute frequency unit and applied to the decimated images. So, it is 
    equivalent to the low-pass filter of decimate*fl with respetive to the original image 
    size.
    
 2. Always use small decimation rate in the first run if there is no prior info about 3D \
   variability analysis of the data and the memory requirement. In addition, one can skip 
   the low-pass filteration when decimation rate is small, which significantly speeds up computation.
   
 3. The program will check the user provided parameters such that the final image size is a 
    product of small primes such as 2, 3, 5...
 
"""
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
	usage = progname + " prj_stack  --ave2D= --var2D=  --ave3D= --var3D= --img_per_grp= --fl=  --aa=   --sym=symmetry --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)
	
	parser.add_option("--output_dir",   type="string"	   ,	default="./",				    help="Output directory")
	parser.add_option("--ave2D",		type="string"	   ,	default=False,				help="Write to the disk a stack of 2D averages")
	parser.add_option("--var2D",		type="string"	   ,	default=False,				help="Write to the disk a stack of 2D variances")
	parser.add_option("--ave3D",		type="string"	   ,	default=False,				help="Write to the disk reconstructed 3D average")
	parser.add_option("--var3D",		type="string"	   ,	default=False,				help="Compute 3D variability (time consuming!)")
	parser.add_option("--img_per_grp",	type="int"         ,	default=100,	     	    help="Number of neighbouring projections.(Default is 100)")
	parser.add_option("--no_norm",		action="store_true",	default=False,				help="Do not use normalization.(Default is to apply normalization)")
	#parser.add_option("--radius", 	    type="int"         ,	default=-1   ,				help="radius for 3D variability" )
	parser.add_option("--npad",			type="int"         ,	default=2    ,				help="Number of time to pad the original images.(Default is 2 times padding)")
	parser.add_option("--sym" , 		type="string"      ,	default="c1",				help="Symmetry. (Default is no symmetry)")
	parser.add_option("--fl",			type="float"       ,	default=0.0,				help="Low pass filter cutoff in absolute frequency (0.0 - 0.5) and is applied to decimated images. (Default - no filtration)")
	parser.add_option("--aa",			type="float"       ,	default=0.02 ,				help="Fall off of the filter. Use default value if user has no clue about falloff (Default value is 0.02)")
	parser.add_option("--CTF",			action="store_true",	default=False,				help="Use CFT correction.(Default is no CTF correction)")
	#parser.add_option("--MPI" , 		action="store_true",	default=False,				help="use MPI version")
	#parser.add_option("--radiuspca", 	type="int"         ,	default=-1   ,				help="radius for PCA" )
	#parser.add_option("--iter", 		type="int"         ,	default=40   ,				help="maximum number of iterations (stop criterion of reconstruction process)" )
	#parser.add_option("--abs", 		type="float"   ,        default=0.0  ,				help="minimum average absolute change of voxels' values (stop criterion of reconstruction process)" )
	#parser.add_option("--squ", 		type="float"   ,	    default=0.0  ,				help="minimum average squared change of voxels' values (stop criterion of reconstruction process)" )
	parser.add_option("--VAR" , 		action="store_true",	default=False,				help="Stack of input consists of 2D variances (Default False)")
	parser.add_option("--decimate",     type  ="float",         default=0.25,               help="Image decimate rate, a number less than 1. (Default is 0.25)")
	parser.add_option("--window",       type  ="int",           default=0,                  help="Target image size relative to original image size. (Default value is zero.)")
	#parser.add_option("--SND",			action="store_true",	default=False,				help="compute squared normalized differences (Default False)")
	#parser.add_option("--nvec",			type="int"         ,	default=0    ,				help="Number of eigenvectors, (Default = 0 meaning no PCA calculated)")
	parser.add_option("--symmetrize",	action="store_true",	default=False,				help="Prepare input stack for handling symmetry (Default False)")
	parser.add_option("--overhead",     type  ="float",         default=0.5,                help="python overhead per CPU.")

	(options,args) = parser.parse_args()
	#####
	from mpi import mpi_init, mpi_comm_rank, mpi_comm_size, mpi_recv, MPI_COMM_WORLD
	from mpi import mpi_barrier, mpi_reduce, mpi_bcast, mpi_send, MPI_FLOAT, MPI_SUM, MPI_INT, MPI_MAX
	#from mpi import *
	from applications   import MPI_start_end
	from reconstruction import recons3d_em, recons3d_em_MPI
	from reconstruction	import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from utilities      import read_text_row, get_image, get_im, wrap_mpi_send, wrap_mpi_recv
	from utilities      import bcast_EMData_to_all, bcast_number_to_all
	from utilities      import get_symt

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
	if options.output_dir =="./": current_output_dir = os.path.abspath(options.output_dir)
	else: current_output_dir = options.output_dir
	if options.symmetrize :
		if RUNNING_UNDER_MPI:
			try:
				sys.argv = mpi_init(len(sys.argv), sys.argv)
				try:	
					number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
					if( number_of_proc > 1 ):
						ERROR("Cannot use more than one CPU for symmetry preparation","sx3dvariability",1)
				except:
					pass
			except:
				pass
		if not os.path.exists(current_output_dir): os.mkdir(current_output_dir)
		
		#  Input
		#instack = "Clean_NORM_CTF_start_wparams.hdf"
		#instack = "bdb:data"
		
		
		from logger import Logger,BaseLogger_Files
		if os.path.exists(os.path.join(current_output_dir, "log.txt")): os.remove(os.path.join(current_output_dir, "log.txt"))
		log_main=Logger(BaseLogger_Files())
		log_main.prefix = os.path.join(current_output_dir, "./")
		
		instack = args[0]
		sym = options.sym.lower()
		if( sym == "c1" ):
			ERROR("There is no need to symmetrize stack for C1 symmetry","sx3dvariability",1)
		
		line =""
		for a in sys.argv:
			line +=" "+a
		log_main.add(line)
	
		if(instack[:4] !="bdb:"):
			#if output_dir =="./": stack = "bdb:data"
			stack = "bdb:"+current_output_dir+"/data"
			delete_bdb(stack)
			junk = cmdexecute("sxcpy.py  "+instack+"  "+stack)
		else: stack = instack
		
		qt = EMUtil.get_all_attributes(stack,'xform.projection')

		na = len(qt)
		ts = get_symt(sym)
		ks = len(ts)
		angsa = [None]*na
		
		for k in range(ks):
			#Qfile = "Q%1d"%k
			#if options.output_dir!="./": Qfile = os.path.join(options.output_dir,"Q%1d"%k)
			Qfile = os.path.join(current_output_dir, "Q%1d"%k)
			#delete_bdb("bdb:Q%1d"%k)
			delete_bdb("bdb:"+Qfile)
			#junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
			junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:"+Qfile)
			#DB = db_open_dict("bdb:Q%1d"%k)
			DB = db_open_dict("bdb:"+Qfile)
			for i in range(na):
				ut = qt[i]*ts[k]
				DB.set_attr(i, "xform.projection", ut)
				#bt = ut.get_params("spider")
				#angsa[i] = [round(bt["phi"],3)%360.0, round(bt["theta"],3)%360.0, bt["psi"], -bt["tx"], -bt["ty"]]
			#write_text_row(angsa, 'ptsma%1d.txt'%k)
			#junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
			#junk = cmdexecute("sxheader.py  bdb:Q%1d  --params=xform.projection  --import=ptsma%1d.txt"%(k,k))
			DB.close()
		#if options.output_dir =="./": delete_bdb("bdb:sdata")
		delete_bdb("bdb:" + current_output_dir + "/"+"sdata")
		#junk = cmdexecute("e2bdb.py . --makevstack=bdb:sdata --filt=Q")
		sdata = "bdb:"+current_output_dir+"/"+"sdata"
		print(sdata)
		junk = cmdexecute("e2bdb.py   " + current_output_dir +"  --makevstack="+sdata +" --filt=Q")
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
		shared_comm  = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,  0, MPI_INFO_NULL)
		myid_on_node = mpi_comm_rank(shared_comm)
		no_of_processes_per_group = mpi_comm_size(shared_comm)
		masters_from_groups_vs_everything_else_comm = mpi_comm_split(MPI_COMM_WORLD, main_node == myid_on_node, myid_on_node)
		color, no_of_groups, balanced_processor_load_on_nodes = get_colors_and_subsets(main_node, MPI_COMM_WORLD, myid, \
		    shared_comm, myid_on_node, masters_from_groups_vs_everything_else_comm)
		overhead_loading = options.overhead*number_of_proc
		#memory_per_node  = options.memory_per_node
		#if memory_per_node == -1.: memory_per_node = 2.*no_of_processes_per_group
		keepgoing = 1
		
		current_window   = options.window
		current_decimate = options.decimate
		
		if len(args) == 1: stack = args[0]
		else:
			print(( "usage: " + usage))
			print(( "Please run '" + progname + " -h' for detailed options"))
			return 1

		t0 = time()	
		# obsolete flags
		options.MPI  = True
		#options.nvec = 0
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
		
		#if options.nvec > 0 :
		#	ERROR("PCA option not implemented", "sx3dvariability", 1, myid)
			
		#if options.nvec > 0 and options.ave3D == None:
		#	ERROR("When doing PCA analysis, one must set ave3D", "sx3dvariability", 1, myid)
		
		if current_decimate>1.0 or current_decimate<0.0:
			ERROR("Decimate rate should be a value between 0.0 and 1.0", "sx3dvariability", 1, myid)
		
		if current_window < 0.0:
			ERROR("Target window size should be always larger than zero", "sx3dvariability", 1, myid)
			
		if myid == main_node:
			img  = get_image(stack, 0)
			nx   = img.get_xsize()
			ny   = img.get_ysize()
			if(min(nx, ny) < current_window):   keepgoing = 0
		keepgoing = bcast_number_to_all(keepgoing, main_node, MPI_COMM_WORLD)
		if keepgoing == 0: ERROR("The target window size cannot be larger than the size of decimated image", "sx3dvariability", 1, myid)

		import string
		options.sym = options.sym.lower()
		# if global_def.CACHE_DISABLE:
		# 	from utilities import disable_bdb_cache
		# 	disable_bdb_cache()
		# global_def.BATCH = True
		
		if myid == main_node:
			if not os.path.exists(current_output_dir): os.mkdir(current_output_dir)# Never delete output_dir in the program!
	
		img_per_grp = options.img_per_grp
		#nvec        = options.nvec
		radiuspca   = options.radiuspca
		from logger import Logger,BaseLogger_Files
		#if os.path.exists(os.path.join(options.output_dir, "log.txt")): os.remove(os.path.join(options.output_dir, "log.txt"))
		log_main=Logger(BaseLogger_Files())
		log_main.prefix = os.path.join(current_output_dir, "./")

		if myid == main_node:
			line = ""
			for a in sys.argv: line +=" "+a
			log_main.add(line)
			log_main.add("-------->>>Settings given by all options<<<-------")
			log_main.add("Symmetry             : %s"%options.sym)
			log_main.add("Input stack          : %s"%stack)
			log_main.add("Output_dir           : %s"%current_output_dir)
			
			if options.ave3D: log_main.add("Ave3d                : %s"%options.ave3D)
			if options.var3D: log_main.add("Var3d                : %s"%options.var3D)
			if options.ave2D: log_main.add("Ave2D                : %s"%options.ave2D)
			if options.var2D: log_main.add("Var2D                : %s"%options.var2D)
			if options.VAR:   log_main.add("VAR                  : True")
			else:             log_main.add("VAR                  : False")
			if options.CTF:   log_main.add("CTF correction       : True  ")
			else:             log_main.add("CTF correction       : False ")
			
			log_main.add("Image per group      : %5d"%options.img_per_grp)
			log_main.add("Image decimate rate  : %4.3f"%current_decimate)
			log_main.add("Low pass filter      : %4.3f"%options.fl)
			current_fl = options.fl
			if current_fl == 0.0: current_fl = 0.5
			log_main.add("Current low pass filter is equivalent to cutoff frequency %4.3f for original image size"%round((current_fl*current_decimate),3))
			log_main.add("Window size          : %5d "%current_window)
			log_main.add("sx3dvariability begins")
	
		symbaselen = 0
		if myid == main_node:
			nima = EMUtil.get_image_count(stack)
			img  = get_image(stack)
			nx   = img.get_xsize()
			ny   = img.get_ysize()
			nnxo = nx
			nnyo = ny
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
			nnxo = 0
			nnyo = 0
		nima    = bcast_number_to_all(nima)
		nx      = bcast_number_to_all(nx)
		ny      = bcast_number_to_all(ny)
		nnxo    = bcast_number_to_all(nnxo)
		nnyo    = bcast_number_to_all(nnyo)
		if current_window > max(nx, ny):
			ERROR("Window size is larger than the original image size", "sx3dvariability", 1)
		
		if current_decimate == 1.:
			if current_window !=0:
				nx = current_window
				ny = current_window
		else:
			if current_window == 0:
				nx = int(nx*current_decimate+0.5)
				ny = int(ny*current_decimate+0.5)
			else:
				nx = int(current_window*current_decimate+0.5)
				ny = nx
		symbaselen = bcast_number_to_all(symbaselen)
		
		# check FFT prime number
		from fundamentals import smallprime
		is_fft_friendly = (nx == smallprime(nx))
		
		if not is_fft_friendly:
			if myid == main_node:
				log_main.add("The target image size is not a product of small prime numbers")
				log_main.add("Program adjusts the input settings!")
			### two cases
			if current_decimate == 1.:
				nx = smallprime(nx)
				ny = nx
				current_window = nx # update
				if myid == main_node:
					log_main.add("The window size is updated to %d."%current_window)
			else:
				if current_window == 0:
					nx = smallprime(int(nx*current_decimate+0.5))
					current_decimate = float(nx)/nnxo
					ny = nx
					if (myid == main_node):
						log_main.add("The decimate rate is updated to %f."%current_decimate)
				else:
					nx = smallprime(int(current_window*current_decimate+0.5))
					ny = nx
					current_window = int(nx/current_decimate+0.5)
					if (myid == main_node):
						log_main.add("The window size is updated to %d."%current_window)
						
		if myid == main_node:
			log_main.add("The target image size is %d"%nx)
						
		if radiuspca == -1: radiuspca = nx/2-2
		if myid == main_node: log_main.add("%-70s:  %d\n"%("Number of projection", nima))
		img_begin, img_end = MPI_start_end(nima, number_of_proc, myid)
		
		"""
		if options.SND:
			from projection		import prep_vol, prgs
			from pap_statistics		import im_diff
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
			for index_of_particle in range(img_begin,img_end):
				image = get_im(stack, index_of_proj)
				if current_window > 0: varList.append(fdecimate(window2d(image,current_window,current_window), nx,ny))
				else:   varList.append(fdecimate(image, nx,ny))
				
		else:
			from utilities		import bcast_number_to_all, bcast_list_to_all, send_EMData, recv_EMData
			from utilities		import set_params_proj, get_params_proj, params_3D_2D, get_params2D, set_params2D, compose_transform2
			from utilities		import model_blank, nearest_proj, model_circle, write_text_row, wrap_mpi_gatherv
			from applications	import pca
			from pap_statistics		import avgvar, avgvar_ctf, ccc
			from filter		    import filt_tanl
			from morphology		import threshold, square_root
			from projection 	import project, prep_vol, prgs
			from sets		    import Set
			from utilities      import wrap_mpi_recv, wrap_mpi_bcast, wrap_mpi_send
			import numpy as np
			if myid == main_node:
				t1          = time()
				proj_angles = []
				aveList     = []
				tab = EMUtil.get_all_attributes(stack, 'xform.projection')	
				for i in range(nima):
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
				log_main.add( "Number of images per group: %d"%img_per_grp)
				log_main.add( "Now grouping projections")
				proj_angles.sort()
				proj_angles_list = np.full((nima, 4), 0.0, dtype=np.float32)	
				for i in range(nima):
					proj_angles_list[i][0] = proj_angles[i][1]
					proj_angles_list[i][1] = proj_angles[i][2]
					proj_angles_list[i][2] = proj_angles[i][3]
					proj_angles_list[i][3] = proj_angles[i][4]
			else: proj_angles_list = 0
			proj_angles_list = wrap_mpi_bcast(proj_angles_list, main_node, MPI_COMM_WORLD)
			proj_angles      = []
			for i in range(nima):
				proj_angles.append([proj_angles_list[i][0], proj_angles_list[i][1], proj_angles_list[i][2], int(proj_angles_list[i][3])])
			del proj_angles_list
			proj_list, mirror_list = nearest_proj(proj_angles, img_per_grp, range(img_begin, img_end))
			all_proj = Set()
			for im in proj_list:
				for jm in im:
					all_proj.add(proj_angles[jm][3])
			all_proj = list(all_proj)
			index = {}
			for i in range(len(all_proj)): index[all_proj[i]] = i
			mpi_barrier(MPI_COMM_WORLD)
			if myid == main_node:
				log_main.add("%-70s:  %.2f\n"%("Finding neighboring projections lasted [s]", time()-t2))
				log_main.add("%-70s:  %d\n"%("Number of groups processed on the main node", len(proj_list)))
				log_main.add("Grouping projections took:  %12.1f [m]"%((time()-t2)/60.))
				log_main.add("Number of groups on main node: ", len(proj_list))
			mpi_barrier(MPI_COMM_WORLD)

			if myid == main_node:
				log_main.add("...... Calculating the stack of 2D variances \n")
			# Memory estimation. There are two memory consumption peaks
			# peak 1. Compute ave, var; 
			# peak 2. Var volume reconstruction;
			# proj_params = [0.0]*(nima*5)
			aveList = []
			varList = []				
			#if nvec > 0: eigList = [[] for i in range(nvec)]
			dnumber   = len(all_proj)# all neighborhood set for assigned to myid
			pnumber   = len(proj_list)*2. + img_per_grp # aveList and varList 
			tnumber   = dnumber+pnumber
			vol_size2 = nx**3*4.*8/1.e9
			vol_size1 = 2.*nnxo**3*4.*8/1.e9
			proj_size         = nnxo*nnyo*len(proj_list)*4.*2./1.e9 # both aveList and varList
			orig_data_size    = nnxo*nnyo*4.*tnumber/1.e9
			reduced_data_size = nx*nx*4.*tnumber/1.e9
			full_data         = np.full((number_of_proc, 2), -1., dtype=np.float16)
			full_data[myid]   = orig_data_size, reduced_data_size
			if myid != main_node: wrap_mpi_send(full_data, main_node, MPI_COMM_WORLD)
			if myid == main_node:
				for iproc in range(number_of_proc):
					if iproc != main_node:
						dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
						full_data[np.where(dummy>-1)] = dummy[np.where(dummy>-1)]
				del dummy
			mpi_barrier(MPI_COMM_WORLD)
			full_data = wrap_mpi_bcast(full_data, main_node, MPI_COMM_WORLD)
			# find the CPU with heaviest load
			minindx         = np.argsort(full_data, 0)
			heavy_load_myid = minindx[-1][1]
			total_mem       = sum(full_data)
			if myid == main_node:
				if current_window == 0:
					log_main.add("Nx:   current image size = %d. Decimated by %f from %d"%(nx, current_decimate, nnxo))
				else:
					log_main.add("Nx:   current image size = %d. Windowed to %d, and decimated by %f from %d"%(nx, current_window, current_decimate, nnxo))
				log_main.add("Nproj:       number of particle images.")
				log_main.add("Navg:        number of 2D average images.")
				log_main.add("Nvar:        number of 2D variance images.")
				log_main.add("Img_per_grp: user defined image per group for averaging = %d"%img_per_grp)
				log_main.add("Overhead:    total python overhead memory consumption   = %f"%overhead_loading)
				log_main.add("Total memory) = 4.0*nx^2*(nproj + navg +nvar+ img_per_grp)/1.0e9 + overhead: %12.3f [GB]"%\
				   (total_mem[1] + overhead_loading))
			del full_data
			mpi_barrier(MPI_COMM_WORLD)
			if myid == heavy_load_myid:
				log_main.add("Begin reading and preprocessing images on processor. Wait... ")
				ttt = time()
			#imgdata = EMData.read_images(stack, all_proj)			
			imgdata = [ None for im in range(len(all_proj))]
			for index_of_proj in range(len(all_proj)):
				#image = get_im(stack, all_proj[index_of_proj])
				if( current_window > 0): imgdata[index_of_proj] = fdecimate(window2d(get_im(stack, all_proj[index_of_proj]),current_window,current_window), nx, ny)
				else:                    imgdata[index_of_proj] = fdecimate(get_im(stack, all_proj[index_of_proj]), nx, ny)
				
				if (current_decimate> 0.0 and options.CTF):
					ctf = imgdata[index_of_proj].get_attr("ctf")
					ctf.apix = ctf.apix/current_decimate
					imgdata[index_of_proj].set_attr("ctf", ctf)
					
				if myid == heavy_load_myid and index_of_proj%100 == 0:
					log_main.add(" ...... %6.2f%% "%(index_of_proj/float(len(all_proj))*100.))
			mpi_barrier(MPI_COMM_WORLD)
			if myid == heavy_load_myid:
				log_main.add("All_proj preprocessing cost %7.2f m"%((time()-ttt)/60.))
				log_main.add("Wait untill reading on all CPUs done...")
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
			from utilities    import model_blank
			from EMAN2        import Transform
			if not options.no_norm: 
				mask = model_circle(nx/2-2, nx, nx)
			if options.CTF: 
				from utilities import pad
				from filter import filt_ctf
			from filter import filt_tanl
			if myid == heavy_load_myid:
				log_main.add("Start computing 2D aveList and varList. Wait...")
				ttt = time()
			inner=nx//2-4
			outer=inner+2
			xform_proj_for_2D = [ None for i in range(len(proj_list))]
			for i in range(len(proj_list)):
				ki = proj_angles[proj_list[i][0]][3]
				if ki >= symbaselen:  continue
				mi = index[ki]
				dpar = Util.get_transform_params(imgdata[mi], "xform.projection", "spider")
				phiM, thetaM, psiM, s2xM, s2yM  = dpar["phi"],dpar["theta"],dpar["psi"],-dpar["tx"]*current_decimate,-dpar["ty"]*current_decimate
				grp_imgdata = []
				for j in range(img_per_grp):
					mj = index[proj_angles[proj_list[i][j]][3]]
					cpar = Util.get_transform_params(imgdata[mj], "xform.projection", "spider")
					alpha, sx, sy, mirror = params_3D_2D_NEW(cpar["phi"], cpar["theta"],cpar["psi"], -cpar["tx"]*current_decimate, -cpar["ty"]*current_decimate, mirror_list[i][j])
					if thetaM <= 90:
						if mirror == 0:  alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, phiM - cpar["phi"], 0.0, 0.0, 1.0)
						else:            alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, 180-(phiM - cpar["phi"]), 0.0, 0.0, 1.0)
					else:
						if mirror == 0:  alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, -(phiM- cpar["phi"]), 0.0, 0.0, 1.0)
						else:            alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, -(180-(phiM - cpar["phi"])), 0.0, 0.0, 1.0)
					imgdata[mj].set_attr("xform.align2d", Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":1.0}))
					grp_imgdata.append(imgdata[mj])
				if not options.no_norm:
					for k in range(img_per_grp):
						ave, std, minn, maxx = Util.infomask(grp_imgdata[k], mask, False)
						grp_imgdata[k] -= ave
						grp_imgdata[k] /= std
				if options.fl > 0.0:
					for k in range(img_per_grp):
						grp_imgdata[k] = filt_tanl(grp_imgdata[k], options.fl, options.aa)

				#  Because of background issues, only linear option works.
				if options.CTF:  ave, var = aves_wiener(grp_imgdata, SNR = 1.0e5, interpolation_method = "linear")
				else:  ave, var = ave_var(grp_imgdata)
				# Switch to std dev
				# threshold is not really needed,it is just in case due to numerical accuracy something turns out negative.
				var = square_root(threshold(var))

				set_params_proj(ave, [phiM, thetaM, 0.0, 0.0, 0.0])
				set_params_proj(var, [phiM, thetaM, 0.0, 0.0, 0.0])

				aveList.append(ave)
				varList.append(var)
				xform_proj_for_2D[i] = [phiM, thetaM, 0.0, 0.0, 0.0]

				'''
				if nvec > 0:
					eig = pca(input_stacks=grp_imgdata, subavg="", mask_radius=radiuspca, nvec=nvec, incore=True, shuffle=False, genbuf=True)
					for k in range(nvec):
						set_params_proj(eig[k], [phiM, thetaM, 0.0, 0.0, 0.0])
						eigList[k].append(eig[k])
					"""
					if myid == 0 and i == 0:
						for k in xrange(nvec):
							eig[k].write_image("eig.hdf", k)
					"""
				'''
				if (myid == heavy_load_myid) and (i%100 == 0):
					log_main.add(" ......%6.2f%%  "%(i/float(len(proj_list))*100.))		
			del imgdata, grp_imgdata, cpar, dpar, all_proj, proj_angles, index
			if not options.no_norm: del mask
			if myid == main_node: del tab
			#  At this point, all averages and variances are computed
			mpi_barrier(MPI_COMM_WORLD)
			
			if (myid == heavy_load_myid):
				log_main.add("Computing aveList and varList took %12.1f [m]"%((time()-ttt)/60.))
			
			xform_proj_for_2D = wrap_mpi_gatherv(xform_proj_for_2D, main_node, MPI_COMM_WORLD)
			if (myid == main_node):
				write_text_row(xform_proj_for_2D, os.path.join(current_output_dir, "params.txt"))
			del xform_proj_for_2D
			mpi_barrier(MPI_COMM_WORLD)
			if options.ave2D:
				from fundamentals import fpol
				from applications import header
				if myid == main_node:
					log_main.add("Compute ave2D ... ")
					km = 0
					for i in range(number_of_proc):
						if i == main_node :
							for im in range(len(aveList)):
								aveList[im].write_image(os.path.join(current_output_dir, options.ave2D), km)
								km += 1
						else:
							nl = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
							nl = int(nl[0])
							for im in range(nl):
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
								tmpvol=fpol(ave, nx, nx,1)								
								tmpvol.write_image(os.path.join(current_output_dir, options.ave2D), km)
								km += 1
				else:
					mpi_send(len(aveList), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					for im in range(len(aveList)):
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
				if myid == main_node:
					header(os.path.join(current_output_dir, options.ave2D), params='xform.projection', fimport = os.path.join(current_output_dir, "params.txt"))
				mpi_barrier(MPI_COMM_WORLD)	
			if options.ave3D:
				from fundamentals import fpol
				t5 = time()
				if myid == main_node: log_main.add("Reconstruct ave3D ... ")
				ave3D = recons3d_4nn_MPI(myid, aveList, symmetry=options.sym, npad=options.npad)
				bcast_EMData_to_all(ave3D, myid)
				if myid == main_node:
					if current_decimate != 1.0: ave3D = resample(ave3D, 1./current_decimate)
					ave3D = fpol(ave3D, nnxo, nnxo, nnxo) # always to the orignal image size
					set_pixel_size(ave3D, 1.0)
					ave3D.write_image(os.path.join(current_output_dir, options.ave3D))
					log_main.add("Ave3D reconstruction took %12.1f [m]"%((time()-t5)/60.0))
					log_main.add("%-70s:  %s\n"%("The reconstructed ave3D is saved as ", options.ave3D))
					
			mpi_barrier(MPI_COMM_WORLD)		
			del ave, var, proj_list, stack, alpha, sx, sy, mirror, aveList
			'''
			if nvec > 0:
				for k in range(nvec):
					if myid == main_node:log_main.add("Reconstruction eigenvolumes", k)
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
						for l in range(len(eigList[k])):
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
			'''
			if options.ave3D: del ave3D
			if options.var2D:
				from fundamentals import fpol 
				from applications import header
				if myid == main_node:
					log_main.add("Compute var2D...")
					km = 0
					for i in range(number_of_proc):
						if i == main_node :
							for im in range(len(varList)):
								tmpvol=fpol(varList[im], nx, nx,1)
								tmpvol.write_image(os.path.join(current_output_dir, options.var2D), km)
								km += 1
						else:
							nl = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
							nl = int(nl[0])
							for im in range(nl):
								ave = recv_EMData(i, im+i+70000)
								tmpvol=fpol(ave, nx, nx,1)
								tmpvol.write_image(os.path.join(current_output_dir, options.var2D), km)
								km += 1
				else:
					mpi_send(len(varList), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					for im in range(len(varList)):
						send_EMData(varList[im], main_node, im+myid+70000)#  What with the attributes??
				mpi_barrier(MPI_COMM_WORLD)
				if myid == main_node:
					from applications import header
					header(os.path.join(current_output_dir, options.var2D), params = 'xform.projection',fimport = os.path.join(current_output_dir, "params.txt"))
				mpi_barrier(MPI_COMM_WORLD)
		if options.var3D:
			if myid == main_node: log_main.add("Reconstruct var3D ...")
			t6 = time()
			# radiusvar = options.radius
			# if( radiusvar < 0 ):  radiusvar = nx//2 -3
			res = recons3d_4nn_MPI(myid, varList, symmetry = options.sym, npad=options.npad)
			#res = recons3d_em_MPI(varList, vol_stack, options.iter, radiusvar, options.abs, True, options.sym, options.squ)
			if myid == main_node:
				from fundamentals import fpol
				if current_decimate != 1.0: res	= resample(res, 1./current_decimate)
				res = fpol(res, nnxo, nnxo, nnxo)
				set_pixel_size(res, 1.0)
				res.write_image(os.path.join(current_output_dir, options.var3D))
				log_main.add("%-70s:  %s\n"%("The reconstructed var3D is saved as ", options.var3D))
				log_main.add("Var3D reconstruction took %f12.1 [m]"%((time()-t6)/60.0))
				log_main.add("Total computation time %f12.1 [m]"%((time()-t0)/60.0))
				log_main.add("sx3dvariability finishes")
		from mpi import mpi_finalize
		mpi_finalize()
		
		if RUNNING_UNDER_MPI: global_def.MPI = False

		global_def.BATCH = False

if __name__=="__main__":
	main()

