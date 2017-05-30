#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  
from __future__ import print_function
import  os
import  sys
import  types
import  global_def
from    global_def import *
from    optparse   import OptionParser
from    sparx      import *
from    EMAN2      import *
from    numpy      import array
from    logger     import Logger, BaseLogger_Files

from mpi   	import  *
from math  	import  *
from random import *

import os
import sys
import subprocess
import time
import string
import json
from   sys 	import exit
from   time import localtime, strftime, sleep

global Tracker, Blockdata



# ------------------------------------------------------------------------------------
mpi_init(0, [])
nproc    = mpi_comm_size(MPI_COMM_WORLD)
myid     = mpi_comm_rank(MPI_COMM_WORLD)


Blockdata = {}
#  MPI stuff
Blockdata["nproc"]              = nproc
Blockdata["myid"]               = myid
Blockdata["main_node"]          = 0
Blockdata["shared_comm"]                    = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,  0, MPI_INFO_NULL)
Blockdata["myid_on_node"]                   = mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"]      = mpi_comm_size(Blockdata["shared_comm"])
masters_from_groups_vs_everything_else_comm = mpi_comm_split(MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"])
Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = get_colors_and_subsets(Blockdata["main_node"], MPI_COMM_WORLD, Blockdata["myid"], \
         Blockdata["shared_comm"], Blockdata["myid_on_node"], masters_from_groups_vs_everything_else_comm)
#  We need two nodes for processing of volumes
Blockdata["node_volume"] = [Blockdata["no_of_groups"]-3, Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take three last nodes
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["no_of_processes_per_group"],Blockdata["node_volume"][1]*Blockdata["no_of_processes_per_group"], \
     Blockdata["node_volume"][2]*Blockdata["no_of_processes_per_group"]]
# End of Blockdata: sorting requires at least three nodes, and the used number of nodes be integer times of three
global_def.BATCH = True
global_def.MPI   = True

def compute_average_ctf(mlist, radius):
	from morphology   import ctf_img
	from filter       import filt_ctf, filt_table
	from fundamentals import fft, fftip
	params_list = [None]*len(mlist)
	orig_image_size = mlist[0].get_xsize()
	avgo       = EMData(orig_image_size, orig_image_size, 1, False) #
	avge       = EMData(orig_image_size, orig_image_size, 1, False) # 
	ctf_2_sumo = EMData(orig_image_size, orig_image_size, 1, False)
	ctf_2_sume = EMData(orig_image_size, orig_image_size, 1, False)
	for im in xrange(len(mlist)):
		ctt = ctf_img(orig_image_size, mlist[im].get_attr("ctf"))
		alpha, sx, sy, mr, scale = get_params2D(mlist[im], xform = "xform.align2d")
		tmp = cosinemask(rot_shift2D(mlist[im], alpha, sx, sy, mr), radius)
		params_list[im]= [alpha, sx, sy, mr, scale]
		tmp = fft(tmp)
		Util.mul_img(tmp, ctt)
		#ima_filt = filt_ctf(tmp, ctf_params, dopad=False)
		if im%2 ==0: 
			Util.add_img2(ctf_2_sume, ctt)
			Util.add_img(avge, tmp)
		else:
			Util.add_img2(ctf_2_sumo, ctt)
			Util.add_img(avgo, tmp)

	sumavg  = Util.divn_img(avge, ctf_2_sume)
	sumctf2 = Util.divn_img(avgo, ctf_2_sumo)
	frc = fsc(fft(sumavg), fft(sumctf2))
	frc[1][0] = 1.0
	for ifreq in xrange(1, len(frc[0])):
		frc[1][ifreq] = max(0.0, frc[1][ifreq])
		frc[1][ifreq] = 2.*frc[1][ifreq]/(1.+frc[1][ifreq])
	sumavg  =  Util.addn_img(avgo, avge)
	sumctf2 =  Util.addn_img(ctf_2_sume, ctf_2_sumo)	
	Util.div_img(sumavg, sumctf2)
	sumavg = fft(sumavg)
	return sumavg, frc, params_list

def compute_average_noctf(mlist, radius):
	from fundamentals import fft
	params_list = [None]*len(mlist)
	orig_image_size = mlist[0].get_xsize()
	avgo       = EMData(orig_image_size, orig_image_size, 1, False) #
	avge       = EMData(orig_image_size, orig_image_size, 1, False) # 
	for im in xrange(len(mlist)):
		alpha, sx, sy, mr, scale = get_params2D(mlist[im], xform = "xform.align2d")
		params_list[im]= [alpha, sx, sy, mr, scale]
		tmp = cosinemask(rot_shift2D(mlist[im], alpha, sx, sy, mr), radius)
		tmp = fft(tmp)
		if im%2 ==0: Util.add_img(avge, tmp)
		else:        Util.add_img(avgo, tmp)
	frc = fsc(fft(sumavg), fft(sumctf2))
	frc[1][0] = 1.0
	for ifreq in xrange(1, len(frc[0])):
		frc[1][ifreq] = max(0.0, frc[1][ifreq])
		frc[1][ifreq] = 2.*frc[1][ifreq]/(1.+frc[1][ifreq])
	sumavg  =  Util.addn_img(avgo, avge)
	sumavg = fft(sumavg)
	return sumavg, frc, params_list

def adjust_pw_to_model(image, pixel_size, roo):
	c1 =-4.5
	c2 = 15.0
	c3 = 0.2
	c4 = -1.0
	c5 = 1./5.
	c6 = 0.25 # six params are fitted to Yifan channel model
	rot1 = rops_table(image)
	fil  = [None]*len(rot1)
	if roo is None: # adjusted to the analytic model, See Penczek Methods Enzymol 2010
		pu = []
		for ifreq in xrange(len(rot1)):
			x = float(ifreq)/float(len(rot1))/pixel_size
			v = exp(c1+c2/(x/c3+1)**2) + exp(c4-0.5*(((x-c5)/c6**2)**2))
			pu.append(v)
		s =sum(pu)
		for ifreq in xrange(len(rot1)): fil[ifreq] = sqrt(pu[ifreq]/(rot1[ifreq]*s))
	else: # adjusted to a given 1-d rotational averaged pw2
		if roo[0]<0.1 or roo[0]>1.: s =sum(roo)
		else:  s=1.0
		for ifreq in xrange(len(rot1)):fil[ifreq] = sqrt(roo[ifreq]/(rot1[ifreq]*s))
	return filt_table(image, fil)
	
def get_optimistic_res(frc):
	nfh = 0
	np  = 0
	for im in xrange(len(frc[1])):
		ifreq = len(frc[1])-1-im
		if frc[1][ifreq] >=0.143:
			np +=1
			nfh = ifreq
			if np >=3:break	
	FH = frc[0][nfh]
	if FH < 0.15:  FH = 0.15 # minimum freq
	return FH
	
def apply_enhancement(avg, B_start, pixel_size, user_defined_Bfactor):	
	guinierline = rot_avg_table(power(periodogram(avg),.5))
	freq_max   =  1./(2.*pixel_size)
	freq_min   =  1./B_start
	b, junk, ifreqmin, ifreqmax = compute_bfactor(guinierline, freq_min, freq_max, pixel_size)
	print(ifreqmin, ifreqmax)
	global_b = b*4. #
	if user_defined_Bfactor!=0.0: global_b = user_defined_Bfactor
	sigma_of_inverse = sqrt(2./global_b)
	avg = filt_gaussinv(avg, sigma_of_inverse)
	return avg, global_b

def main():
	from optparse   import OptionParser
	from global_def import SPARXVERSION
	from EMAN2      import EMData
	from logger     import Logger, BaseLogger_Files
	import sys, os, time
	global Tracker, Blockdata
	from global_def import ERROR
	       
	progname = os.path.basename(sys.argv[0])
	usage = progname + " --output_dir=output_dir  --isac_dir=output_dir_of_isac "
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--isac_dir",              type   ="string",     default ='',     help="ISAC run output directory, input directory for this command")
	parser.add_option("--output_dir",            type   ="string",     default ='',     help="output directory where computed averages are saved")
	parser.add_option("--pixel_size",            type   ="float",      default =-1.0,   help="pixel_size of raw images")
	parser.add_option("--noctf",                 action ="store_true", default =False,  help="no ctf correction, useful for negative stained data. always ctf for cryo data")
	parser.add_option("--B_enhance",             action ="store_true", default =False,  help="apply B_factor to enhance")
	parser.add_option("--B_start",               type   ="float",      default = 10.0,  help="start frequency (1./Angstrom) of power spectrum for B_factor estimation")
	parser.add_option("--Bfactor",               type   ="float",      default = 45.0,  help= "User defined bactors")
	parser.add_option("--fl",                    type   ="float",      default =-1.0,   help= "low pass filter ")
	parser.add_option("--stack",                 type   ="string",     default ="",     help= "data stack that ISAC run used")
	parser.add_option("--radius",                type   ="int",        default =-1,     help= "radius")
	parser.add_option("--xr",                    type   ="float",      default =-1.0,   help= "local alignment search range")
	parser.add_option("--ts",                    type   ="float",      default =1.0,    help= "local alignment search step")
	parser.add_option("--fh",                    type   ="float",      default =-1.,    help= "local alignment high frequencies limit")
	parser.add_option("--maxit",                 type   ="int",        default =5,      help= "local alignment iterations")
	parser.add_option("--nopwadj",               action ="store_true", default =False,  help= "no pw adjustment")
	parser.add_option("--modelpw",               type   ="string",     default ='',     help= "1-D power spectrum of PDF model or EM map sampled in given pixel_size and orignal image_size")
	parser.add_option("--navg",                  type   ="int",        default =-1,     help= "number of aveages")
	
	(options, args) = parser.parse_args(sys.argv[1:])

	#--- Fill input parameters into dictionary Constants
	Constants		                         = {}
	Constants["isac_dir"]                     = options.isac_dir
	Constants["masterdir"]                    = options.output_dir
	Constants["pixel_size"]                   = options.pixel_size
	Constants["orgstack"]                     = options.stack
	Constants["radius"]                       = options.radius
	Constants["xrange"]                       = options.xr
	Constants["xstep"]                        = options.ts
	Constants["FH"]                           = options.fh
	Constants["maxit"]                        = options.maxit
	Constants["B_enhance"]                    = options.B_enhance
	Constants["B_start"]                      = options.B_start
	Constants["Bfactor"]                      = options.Bfactor
	Constants["nopwadj"]                      = options.nopwadj
	Constants["modelpw"]                      = options.modelpw
	Constants["low_pass_filter"]              = options.fl
	Constants["navg"]                         = options.navg
	Constants["noctf"]                        = options.noctf
	
	# -------------------------------------------------------------
	#
	# Create and initialize Tracker dictionary with input options  # State Variables
	
	Tracker              = {}
	Tracker["constants"] = Constants
	#<<<---------------------->>>imported functions<<<---------------------------------------------

	from statistics 	import k_means_match_clusters_asg_new,k_means_stab_bbenum
	from applications 	import recons3d_n_MPI
	from utilities 		import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
	from utilities 		import get_initial_ID, print_upper_triangular_matrix, print_a_line_with_timestamp
	from utilities 		import get_resolution_mrk01,partition_to_groups,partition_independent_runs,get_outliers
	from utilities 		import merge_groups, save_alist, margin_of_error, get_margin_of_error, get_ali3d_params
	from utilities 		import counting_projections, unload_dict, load_dict, get_stat_proj, get_number_of_groups, recons_mref
	from utilities 		import apply_low_pass_filter, get_groups_from_partition, get_number_of_groups, get_complementary_elements_total, update_full_dict
	from utilities 		import count_chunk_members, set_filter_parameters_from_adjusted_fsc, adjust_fsc_down, get_two_chunks_from_stack
	from utilities 		import get_input_from_string, cmdexecute
	from filter			import filt_tanl
	from time           import sleep
	from logger         import Logger,BaseLogger_Files
	import user_functions
	import string
	from string         import split, atoi, atof
	import json

	#x_range = max(Tracker["constants"]["xrange"], int(1./Tracker["ini_shrink"])+1)
	#y_range =  x_range

	####-----------------------------------------------------------
	# Create Master directory
	line      = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Tracker["constants"]["masterdir"] == Tracker["constants"]["isac_dir"]:
		masterdir = os.path.join(Tracker["constants"]["isac_dir"], "sharpen")
	else: masterdir = Tracker["constants"]["masterdir"]
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg = "Postprocessing ISAC 2D averages starts"
		print(line, "Postprocessing ISAC 2D averages starts")
		if not masterdir:
			timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			masterdir ="sharpen_"+Tracker["constants"]["isac_dir"]
			cmd="{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
		else:
			if os.path.exists(masterdir): restart = 1
			else:
				cmd="{} {}".format("mkdir", masterdir)
				cmdexecute(cmd)
		li =len(masterdir)
	else: li = 0
	li                                  = mpi_bcast(li,1,MPI_INT,Blockdata["main_node"],MPI_COMM_WORLD)[0]
	masterdir							= mpi_bcast(masterdir,li,MPI_CHAR,Blockdata["main_node"],MPI_COMM_WORLD)
	masterdir                           = string.join(masterdir,"")
	Tracker["constants"]["masterdir"]	= masterdir
	log_main = Logger(BaseLogger_Files())
	log_main.prefix = Tracker["constants"]["masterdir"]+"/"
	
	while not os.path.exists(Tracker["constants"]["masterdir"]):
		print("Node ", Blockdata["myid"], "  waiting...", Tracker["constants"]["masterdir"])
		sleep(1)
	mpi_barrier(MPI_COMM_WORLD)

	if(Blockdata["myid"] == Blockdata["main_node"]):
		init_dict = {}
		print(Tracker["constants"]["isac_dir"])
		Tracker["directory"]   = os.path.join(Tracker["constants"]["isac_dir"], "2dalignment")
		core = read_text_row(os.path.join(Tracker["directory"], "initial2Dparams.txt"))
		for im in xrange(len(core)):
			init_dict[im]  = core[im]
		del core
	else: init_dict = 0
	init_dict = wrap_mpi_bcast(init_dict, Blockdata["main_node"], communicator = MPI_COMM_WORLD)	
	###
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		#Tracker["constants"]["orgstack"] = "bdb:"+ os.path.join(Tracker["constants"]["isac_dir"],"../","sparx_stack")
		image = get_im(Tracker["constants"]["orgstack"], 0)
		Tracker["constants"]["nnxo"]     = image.get_xsize()
		try:
			ctf_params = image.get_attr("ctf")
			if Tracker["constants"]["pixel_size"] ==-1.: Tracker["constants"]["pixel_size"] = ctf_params.apix
		except: print(" we don't get pixel size value") 
			
		os.path.join(Tracker["directory"], "aqfinal.hdf")
		Tracker["ini_shrink"] = float(get_im(os.path.join(Tracker["directory"], "aqfinal.hdf"), 0).get_xsize())/Tracker["constants"]["nnxo"]
	else: Tracker["ini_shrink"] = 0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	
	#print(Tracker["constants"]["pixel_size"], "pixel_size")	
	x_range = max(Tracker["constants"]["xrange"], int(1./Tracker["ini_shrink"])+1)
	y_range =  x_range
	
	if(Blockdata["myid"] == Blockdata["main_node"]): parameters = read_text_row(os.path.join(Tracker["constants"]["isac_dir"], "all_parameters.txt"))
	else: parameters = 0
	parameters = wrap_mpi_bcast(parameters, Blockdata["main_node"], communicator = MPI_COMM_WORLD)		
	params_dict = {}
	list_dict   = {}
	#parepare params_dict
	
	if Tracker["constants"]["navg"] <0: navg = EMUtil.get_image_count(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf"))
	else: navg = min(Tracker["constants"]["navg"], EMUtil.get_image_count(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf")))
	if(Blockdata["myid"] == Blockdata["main_node"]):
		for iavg in xrange(navg):
			params_of_this_average = []
			image   = get_im(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf"), iavg)
			members = image.get_attr("members")
			for im in xrange(len(members)):
				abs_id =  members[im]
				P = combine_params2( init_dict[abs_id][0], init_dict[abs_id][1], init_dict[abs_id][2], init_dict[abs_id][3], \
				parameters[abs_id][0], parameters[abs_id][1]/Tracker["ini_shrink"], parameters[abs_id][2]/Tracker["ini_shrink"], parameters[abs_id][3])
				if parameters[abs_id][3] ==-1: print("wrong one")
				params_of_this_average.append([P[0], P[1], P[2], P[3], 1.0])
			params_dict[iavg] = params_of_this_average
			list_dict[iavg] = members
			write_text_row(params_of_this_average, os.path.join(Tracker["constants"]["masterdir"], "params_avg_%03d.txt"%iavg))
	else:  
		params_dict = 0
		list_dict   = 0
	params_dict = wrap_mpi_bcast(params_dict, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	list_dict = wrap_mpi_bcast(list_dict, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	# Now computing!
	del init_dict
	tag_sharpen_avg = 1000
	if navg <Blockdata["nproc"]: # 
		for iavg in xrange(navg):
			if Blockdata["myid"] == iavg:
				mlist = [None for i in xrange(len(list_dict[iavg]))]
				for im in xrange(len(mlist)):
					mlist[im]= get_im(Tracker["constants"]["orgstack"], list_dict[iavg][im])
					set_params2D(mlist[im], params_dict[iavg][im], xform = "xform.align2d")
				if Tracker["constants"]["noctf"]: ini_avg, frc = compute_average_noctf(mlist, Tracker["constants"]["radius"])
				else: ini_avg, frc, plist = compute_average_ctf(mlist, Tracker["constants"]["radius"])
				FH1 = get_optimistic_res(frc)
				#write_text_file(frc, os.path.join(Tracker["constants"]["masterdir"], "fsc%03d_before_ali.txt"%iavg))
				new_average1 = within_group_refinement([mlist[kik] for kik in xrange(0,len(mlist),2)], maskfile= None, randomize= False, ir=1.0,  \
				ou=Tracker["constants"]["radius"], rs=1.0, xrng=[x_range], yrng=[y_range], step=[Tracker["constants"]["xstep"]], \
				dst=0.0, maxit=Tracker["constants"]["maxit"], FH = max(Tracker["constants"]["FH"], FH1), FF=0.1)
				new_average2 = within_group_refinement([mlist[kik] for kik in xrange(1,len(mlist),2)], maskfile= None, randomize= False, ir=1.0, \
				ou=Tracker["constants"]["radius"], rs=1.0, xrng=[x_range], yrng=[y_range], step=[Tracker["constants"]["xstep"]], \
				dst=0.0, maxit=Tracker["constants"]["maxit"], FH = max(Tracker["constants"]["FH"], FH1), FF=0.1)
				if Tracker["constants"]["noctf"]: new_avg, frc, plist = compute_average_noctf(mlist, Tracker["constants"]["radius"])
				else: new_avg, frc, plist = compute_average_ctf(mlist, Tracker["constants"]["radius"])
				FH2 = get_optimistic_res(frc)
				#write_text_file(frc, os.path.join(Tracker["constants"]["masterdir"], "fsc%03d.txt"%iavg))
				if Tracker["constants"]["nopwadj"]: # pw adjustment, 1. analytic model 2. PDB model 3. B-facttor enhancement
					if Tracker["constants"]["B_enhance"]:
						new_avg, gb = apply_enhancement(new_avg, Tracker["constants"]["B_start"], Tracker["constants"]["pixel_size"], Tracker["constants"]["Bfactor"])
						print("Process avg  %d  %f  %f   %f"%(iavg, gb, FH1, FH2))
				else:
					try: 
						roo = read_text_file(Tracker["constants"]["modelpw"], -1)
						roo =roo[0] # always put pw in the first column
					except: roo = None
					new_avg = adjust_pw_to_model(new_avg, Tracker["constants"]["pixel_size"], roo)
					print("Process avg  %d   %f   %f"%(iavg, FH1, FH2))
				if Tracker["constants"]["low_pass_filter"] !=-1.: new_avg = filt_tanl(new_avg, Tracker["constants"]["low_pass_filter"], 0.1)
				# write members into header:
				new_avg.set_attr("members", list_dict[iavg])
				new_avg.set_attr("n_objects", len(list_dict[iavg]))
		mpi_barrier(MPI_COMM_WORLD)
		for im in xrange(navg): # avg
			if im == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(new_avg, Blockdata["main_node"],  tag_sharpen_avg)
			elif Blockdata["myid"] == Blockdata["main_node"]:
				if im != Blockdata["main_node"]:
					new_avg_other_cpu = recv_EMData(im, tag_sharpen_avg)
					new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
				else: new_avg.write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
			if im == Blockdata["myid"]:
				write_text_row(plist, os.path.join(Tracker["constants"]["masterdir"], "ali2d_local_params_avg_%03d.txt"%im))
		mpi_barrier(MPI_COMM_WORLD)
		"""
		for im in xrange(navg): # ini_avg
			if im == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(ini_avg, Blockdata["main_node"],  tag_sharpen_avg)
			elif Blockdata["myid"] == Blockdata["main_node"]:
				if im != Blockdata["main_node"]:
					new_avg_other_cpu = recv_EMData(im, tag_sharpen_avg)
					new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "ini_class_averages.hdf"), im)
				else: ini_avg.write_image(os.path.join(Tracker["constants"]["masterdir"], "ini_class_averages.hdf"), im)
		mpi_barrier(MPI_COMM_WORLD)
		
		for im in xrange(navg): # avg1
			if im == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(new_average1, Blockdata["main_node"],  tag_sharpen_avg)
			elif Blockdata["myid"] == Blockdata["main_node"]:
				if im != Blockdata["main_node"]:
					new_avg_other_cpu = recv_EMData(im, tag_sharpen_avg)
					new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "ali_class_averages1.hdf"), im)
				else: new_average1.write_image(os.path.join(Tracker["constants"]["masterdir"], "ali_class_averages1.hdf"), im)
		mpi_barrier(MPI_COMM_WORLD)
		
		for im in xrange(navg):# avg2
			if im == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(new_average2, Blockdata["main_node"],  tag_sharpen_avg)
			elif Blockdata["myid"] == Blockdata["main_node"]:
				if im != Blockdata["main_node"]:
					new_avg_other_cpu = recv_EMData(im, tag_sharpen_avg)
					new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "ali_class_averages2.hdf"), im)
				else: new_average2.write_image(os.path.join(Tracker["constants"]["masterdir"], "ali_class_averages2.hdf"), im)
		mpi_barrier(MPI_COMM_WORLD)
		"""
	else:
		image_start,image_end = MPI_start_end(navg, Blockdata["nproc"], Blockdata["myid"])
		if Blockdata["myid"] == Blockdata["main_node"]:
			cpu_dict = {}
			for iproc in xrange(Blockdata["nproc"]):
				local_image_start, local_image_end = MPI_start_end(navg, Blockdata["nproc"], iproc)
				for im in xrange(local_image_start, local_image_end): cpu_dict [im] = iproc
		else:  cpu_dict = 0
		cpu_dict = wrap_mpi_bcast(cpu_dict, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
		
		slist     = [None for im in xrange(navg)]
		ini_list  = [None for im in xrange(navg)]
		avg1_list = [None for im in xrange(navg)]
		avg2_list = [None for im in xrange(navg)]
		plist_dict ={}
		for iavg in xrange(image_start,image_end):
			#print("Process avg %d"%iavg)
			mlist = [None for i in xrange(len(list_dict[iavg]))]
			for im in xrange(len(mlist)):
				mlist[im]= get_im(Tracker["constants"]["orgstack"], list_dict[iavg][im])
				set_params2D(mlist[im], params_dict[iavg][im], xform = "xform.align2d")
			if Tracker["constants"]["noctf"]: ini_avg, frc = compute_average_noctf(mlist, Tracker["constants"]["radius"])
			else:   ini_avg, frc, plist = compute_average_ctf(mlist, Tracker["constants"]["radius"])
			FH1 = get_optimistic_res(frc)
			#write_text_file(frc, os.path.join(Tracker["constants"]["masterdir"], "fsc%03d_before_ali.txt"%iavg))
			new_average1 = within_group_refinement([mlist[kik] for kik in xrange(0,len(mlist),2)], maskfile= None, randomize= False, ir=1.0,  \
			 ou=Tracker["constants"]["radius"], rs=1.0, xrng=[x_range], yrng=[y_range], step=[Tracker["constants"]["xstep"]], \
			 dst=0.0, maxit=Tracker["constants"]["maxit"], FH=max(Tracker["constants"]["FH"], FH1), FF=0.1)
			new_average2 = within_group_refinement([mlist[kik] for kik in xrange(1,len(mlist),2)], maskfile= None, randomize= False, ir=1.0, \
			 ou= Tracker["constants"]["radius"], rs=1.0, xrng=[ x_range], yrng=[y_range], step=[Tracker["constants"]["xstep"]], \
			 dst=0.0, maxit=Tracker["constants"]["maxit"], FH = max(Tracker["constants"]["FH"], FH1), FF=0.1)
			if Tracker["constants"]["noctf"]: new_avg, frc,  plist = compute_average_noctf(mlist, Tracker["constants"]["radius"])
			else: new_avg, frc, plist = compute_average_ctf(mlist, Tracker["constants"]["radius"])
			plist_dict[iavg] = plist
			FH2 = get_optimistic_res(frc)
			#write_text_file(frc, os.path.join(Tracker["constants"]["masterdir"], "fsc%03d.txt"%iavg))
			if not Tracker["constants"]["nopwadj"]:
				try: 
					roo = read_text_file( Tracker["constants"]["modelpw"], -1)
					roo = roo[0] # always on the first column
				except: roo = None 
				new_avg = adjust_pw_to_model(new_avg, Tracker["constants"]["pixel_size"], roo)
				print("Process avg  %d  %f  %f"%(iavg, FH1, FH2))
			else:
				if Tracker["constants"]["B_enhance"]:
					new_avg, gb = apply_enhancement(new_avg, Tracker["constants"]["B_start"], Tracker["constants"]["pixel_size"], Tracker["constants"]["Bfactor"])
					print("Process avg  %d  %f  %f  %f"%(iavg, gb, FH1, FH2))
				
			if Tracker["constants"]["low_pass_filter"] !=-1.:new_avg = filt_tanl(new_avg, Tracker["constants"]["low_pass_filter"], 0.1)		
			new_avg.set_attr("members", list_dict[iavg])
			new_avg.set_attr("n_objects", len(list_dict[iavg]))
			slist[iavg]    = new_avg
			ini_list[iavg] = ini_avg
			avg1_list[iavg]= new_average1
			avg2_list[iavg]= new_average2
		## send to main node to write
		for im in xrange(navg):
			# avg
			if cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(slist[im], Blockdata["main_node"],  tag_sharpen_avg)
				
			elif cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				slist[im].write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
				
			elif cpu_dict[im] != Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				new_avg_other_cpu = recv_EMData(cpu_dict[im], tag_sharpen_avg)
				new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
			else: pass
			"""
			# ini
			if cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(ini_list[im], Blockdata["main_node"],  tag_sharpen_avg)
				
			elif cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				ini_list[im].write_image(os.path.join(Tracker["constants"]["masterdir"], "ini_class_averages.hdf"), im)
				
			elif cpu_dict[im] != Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				new_avg_other_cpu = recv_EMData(cpu_dict[im], tag_sharpen_avg)
				new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "ini_class_averages.hdf"), im)
			else: pass
			
			# avg1
			if cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(avg1_list[im], Blockdata["main_node"],  tag_sharpen_avg)
				
			elif cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				avg1_list[im].write_image(os.path.join(Tracker["constants"]["masterdir"], "avg1_class_averages.hdf"), im)
				
			elif cpu_dict[im] != Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				new_avg_other_cpu = recv_EMData(cpu_dict[im], tag_sharpen_avg)
				new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "avg1_class_averages.hdf"), im)
			else: pass
			
			# avg2
			if cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
				send_EMData(avg2_list[im], Blockdata["main_node"], tag_sharpen_avg)
				
			elif cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				avg2_list[im].write_image(os.path.join(Tracker["constants"]["masterdir"], "avg2_class_averages.hdf"), im)
				
			elif cpu_dict[im] != Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
				new_avg_other_cpu = recv_EMData(cpu_dict[im], tag_sharpen_avg)
				new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"],"avg2_class_averages.hdf"),im)
			else: pass
			"""
			if cpu_dict[im] == Blockdata["myid"]:
				write_text_row(plist_dict[im], os.path.join(Tracker["constants"]["masterdir"], "ali2d_local_params_avg_%03d.txt"%im))
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
	target_xr =3
	target_yr =3
	
	if( Blockdata["myid"] == 0):
		cmd = "{} {} {} {} {} {} {} {} {} {}".format("sxchains.py", os.path.join(Tracker["constants"]["masterdir"],"class_averages.hdf"),\
		os.path.join(Tracker["constants"]["masterdir"],"junk.hdf"),os.path.join(Tracker["constants"]["masterdir"],"ordered_class_averages.hdf"),\
		"--circular","--radius=%d"%Tracker["constants"]["radius"] , "--xr=%d"%(target_xr+1),"--yr=%d"%(target_yr+1),"--align", ">/dev/null")
		junk = cmdexecute(cmd)
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["constants"]["masterdir"], "junk.hdf") )
		junk = cmdexecute(cmd)
		
	from mpi import mpi_finalize
	mpi_finalize()
	exit()
if __name__ == "__main__":
	main()
