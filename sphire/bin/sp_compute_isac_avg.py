#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  
from __future__ import print_function
from builtins import range

import os
import sys
from   sys import exit
import time
from   time import localtime, strftime, sleep
import types

import sp_global_def
from sp_global_def import sxprint, ERROR, SPARXVERSION
from sp_global_def import *

from   optparse   import OptionParser

from    sp_sparx      import *
from    EMAN2      import *
from    numpy      import array
from    sp_logger     import Logger, BaseLogger_Files

import mpi

from math import *

from random import *

import subprocess
import string
from   string import split, atoi, atof

import json

from optparse   import OptionParser
from EMAN2      import EMData
from sp_logger     import Logger, BaseLogger_Files

from sp_utilities 		import get_im, bcast_number_to_all, write_text_file, read_text_file, wrap_mpi_bcast, write_text_row
from sp_utilities 		import cmdexecute
from sp_filter			import filt_tanl
import sp_user_functions
import json

global Tracker, Blockdata

mpi.mpi_init(0, [])


# ----------------------------------------------------------------------------
def compute_average(mlist, radius, CTF):
	from sp_morphology   import cosinemask
	from sp_fundamentals import fft
	from sp_statistics   import fsc, sum_oe
	if CTF:
		avge, avgo, ctf_2_sume, ctf_2_sumo, params_list = \
		     sum_oe(mlist, "a", CTF, None, True, True)
		avge = cosinemask(fft(avge), radius)
		avgo = cosinemask(fft(avgo), radius)
		sumavge  = Util.divn_img(fft(avge), ctf_2_sume)
		sumavgo  = Util.divn_img(fft(avgo), ctf_2_sumo)
		frc = fsc(fft(sumavgo), fft(sumavge))
		frc[1][0] = 1.0
		for ifreq in range(1, len(frc[0])):
			frc[1][ifreq] = max(0.0, frc[1][ifreq])
			frc[1][ifreq] = 2.*frc[1][ifreq]/(1.+frc[1][ifreq])
		sumavg  =  Util.addn_img(fft(avgo), fft(avge))
		sumctf2 =  Util.addn_img(ctf_2_sume, ctf_2_sumo)	
		Util.div_img(sumavg, sumctf2)
		return fft(sumavg), frc, params_list
	else:
		avge, avgo, params_list = sum_oe(mlist, "a", False,  None, False, True)
		avge = cosinemask(avge, radius)
		avgo = cosinemask(avgo, radius)
		frc = fsc(avgo, avge)
		frc[1][0] = 1.0
		for ifreq in range(1, len(frc[0])):
			frc[1][ifreq] = max(0.0, frc[1][ifreq])
			frc[1][ifreq] = 2.*frc[1][ifreq]/(1.+frc[1][ifreq])
		return avge+avgo, frc, params_list

def adjust_pw_to_model(image, pixel_size, roo):
	from sp_fundamentals import rops_table
	from sp_filter       import filt_table
	from math         import exp, sqrt
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
		for ifreq in range(len(rot1)):
			x = float(ifreq)/float(len(rot1))/pixel_size
			v = exp(c1+c2/(x/c3+1)**2) + exp(c4-0.5*(((x-c5)/c6**2)**2))
			pu.append(v)
		s =sum(pu)
		for ifreq in range(len(rot1)): fil[ifreq] = sqrt(pu[ifreq]/(rot1[ifreq]*s))
	else: # adjusted to a given 1-d rotational averaged pw2
		if roo[0]<0.1 or roo[0]>1.: s =sum(roo)
		else:  s=1.0
		for ifreq in range(len(rot1)):fil[ifreq] = sqrt(roo[ifreq]/(rot1[ifreq]*s))
	return filt_table(image, fil)
	
def get_optimistic_res(frc):
	nfh = 0
	np  = 0
	for im in range(len(frc[1])):
		ifreq = len(frc[1])-1-im
		if frc[1][ifreq] >=0.143:
			np +=1
			nfh = ifreq
			if np >=3:break	
	FH = frc[0][nfh]
	if FH < 0.15:  FH = 0.15 # minimum freq
	return FH
	
def apply_enhancement(avg, B_start, pixel_size, user_defined_Bfactor):
	from sp_filter       import filt_gaussinv
	from sp_fundamentals import rot_avg_table
	from sp_morphology   import compute_bfactor
	from EMAN2        import periodogram
	if user_defined_Bfactor>0.0:
		global_b = user_defined_Bfactor
	else:
		guinierline = rot_avg_table(power(periodogram(fft(avg)),.5))
		freq_max    =  1./(2.*pixel_size)
		freq_min    =  1./B_start
		b, junk, ifreqmin, ifreqmax = compute_bfactor(guinierline, freq_min, freq_max, pixel_size)
		#print(ifreqmin, ifreqmax)
		global_b = b*4. #
	return filt_gaussinv(fft(avg), sqrt(2./global_b)), global_b

def main():
	global Tracker, Blockdata
	progname = os.path.basename(sys.argv[0])
	usage = progname + " --output_dir=output_dir  --isac_dir=output_dir_of_isac "
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--pw_adjustment", type ="string", default ='analytical_model',  \
	   help="adjust power spectrum of 2-D averages to an analytic model. Other opions: no_adjustment; bfactor; a text file of 1D rotationally averaged PW")
	#### Four options for --pw_adjustment: 	
	# 1> analytical_model(default); 
	# 2> no_adjustment;
	# 3> bfactor;
	# 4> adjust_to_given_pw2(user has to provide a text file that contains 1D rotationally averaged PW)

	# options in common
	parser.add_option("--isac_dir",              type   ="string",         default ='',     help="ISAC run output directory, input directory for this command")
	parser.add_option("--output_dir",            type   ="string",         default ='',     help="output directory where computed averages are saved")
	parser.add_option("--pixel_size",            type   ="float",          default =-1.0,   help="pixel_size of raw images. one can put 1.0 in case of negative stain data")
	parser.add_option("--fl",                    type   ="float",          default =-1.0,    help= "low pass filter, = -1.0, not applied; =0.0, using FH1 (initial resolution), = 1.0 using FH2 (resolution after local alignment), or user provided value in absolute freqency [0.0:0.5]")
	parser.add_option("--stack",                 type   ="string",         default ="",     help= "data stack used in ISAC")
	parser.add_option("--radius",                type   ="int",            default =-1,     help= "radius")
	parser.add_option("--xr",                    type   ="float",          default =-1.0,   help= "local alignment search range")
	#parser.add_option("--ts",                    type   ="float",          default =1.0,    help= "local alignment search step")
	parser.add_option("--fh",                    type   ="float",          default =-1.0,   help= "local alignment high frequencies limit")
	#parser.add_option("--maxit",                 type   ="int",            default =5,      help= "local alignment iterations")
	parser.add_option("--navg",                  type   ="int",            default =1000000,     help= "number of aveages")
	parser.add_option("--local_alignment",       action ="store_true",     default =False,  help= "do local alignment")
	parser.add_option("--noctf",                 action ="store_true",     default =False,  help="no ctf correction, useful for negative stained data. always ctf for cryo data")
	parser.add_option("--B_start",  type   ="float",  default = 45.0,  help="start frequency (Angstrom) of power spectrum for B_factor estimation")
	parser.add_option("--Bfactor",  type   ="float",  default = -1.0,  help= "User defined bactors (e.g. 25.0[A^2]). By default, the program automatically estimates B-factor. ")
			
	(options, args) = parser.parse_args(sys.argv[1:])
	
	adjust_to_analytic_model = True if options.pw_adjustment=='analytical_model' else False
	no_adjustment            = True if options.pw_adjustment=='no_adjustment'    else False
	B_enhance                = True if options.pw_adjustment=='bfactor'          else False
	adjust_to_given_pw2      = True if not (adjust_to_analytic_model or no_adjustment or B_enhance) else False

	# mpi
	nproc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	myid  = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)

	Blockdata = {}
	Blockdata["nproc"]     = nproc
	Blockdata["myid"]      = myid
	Blockdata["main_node"] = 0
	Blockdata["shared_comm"]                    = mpi.mpi_comm_split_type(mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED,  0, mpi.MPI_INFO_NULL)
	Blockdata["myid_on_node"]                   = mpi.mpi_comm_rank(Blockdata["shared_comm"])
	Blockdata["no_of_processes_per_group"]      = mpi.mpi_comm_size(Blockdata["shared_comm"])
	masters_from_groups_vs_everything_else_comm = mpi.mpi_comm_split( mpi.MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"])
	Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = get_colors_and_subsets( Blockdata["main_node"],
																											  mpi.MPI_COMM_WORLD, Blockdata["myid"],
																											  Blockdata["shared_comm"], 
																											  Blockdata["myid_on_node"], 
																											  masters_from_groups_vs_everything_else_comm )
	#  We need two nodes for processing of volumes
	Blockdata["node_volume"] = [Blockdata["no_of_groups"]-3, Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take three last nodes
	#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
	#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
	Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["no_of_processes_per_group"],Blockdata["node_volume"][1]*Blockdata["no_of_processes_per_group"], \
		 Blockdata["node_volume"][2]*Blockdata["no_of_processes_per_group"]]
	# End of Blockdata: sorting requires at least three nodes, and the used number of nodes be integer times of three
	sp_global_def.BATCH = True
	sp_global_def.MPI   = True


	if adjust_to_given_pw2:
		checking_flag = 0
		if(Blockdata["myid"] == Blockdata["main_node"]):
			if not os.path.exists(options.pw_adjustment): checking_flag = 1
		checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		
		if checking_flag ==1: 
			ERROR( "User provided power spectrum does not exist", myid=Blockdata["myid"] )
	
	Tracker                                   = {}
	Constants		                          = {}
	Constants["isac_dir"]                     = options.isac_dir
	Constants["masterdir"]                    = options.output_dir
	Constants["pixel_size"]                   = options.pixel_size
	Constants["orgstack"]                     = options.stack
	Constants["radius"]                       = options.radius
	Constants["xrange"]                       = options.xr
	Constants["FH"]                           = options.fh
	Constants["low_pass_filter"]              = options.fl
	#Constants["maxit"]                        = options.maxit
	Constants["navg"]                         = options.navg
	Constants["B_start"]                      = options.B_start
	Constants["Bfactor"]                      = options.Bfactor
	
	if adjust_to_given_pw2: Constants["modelpw"] = options.pw_adjustment
	Tracker["constants"] = Constants
	# -------------------------------------------------------------
	#
	# Create and initialize Tracker dictionary with input options  # State Variables


	#<<<---------------------->>>imported functions<<<---------------------------------------------

	#x_range = max(Tracker["constants"]["xrange"], int(1./Tracker["ini_shrink"])+1)
	#y_range =  x_range

	####-----------------------------------------------------------
	# Create Master directory and associated subdirectories
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Tracker["constants"]["masterdir"] == Tracker["constants"]["isac_dir"]:
		masterdir = os.path.join(Tracker["constants"]["isac_dir"], "sharpen")
	else: masterdir = Tracker["constants"]["masterdir"]

	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg = "Postprocessing ISAC 2D averages starts"
		sxprint(line, "Postprocessing ISAC 2D averages starts")
		if not masterdir:
			timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			masterdir ="sharpen_"+Tracker["constants"]["isac_dir"]
			os.makedirs(masterdir)
		else:
			if os.path.exists(masterdir): sxprint("%s already exists"%masterdir)
			else:
				os.makedirs(masterdir)
		sp_global_def.write_command(masterdir)
		subdir_path = os.path.join(masterdir, "ali2d_local_params_avg")
		if not os.path.exists(subdir_path): os.mkdir(subdir_path)
		subdir_path = os.path.join(masterdir, "params_avg")
		if not os.path.exists(subdir_path): os.mkdir(subdir_path)
		li =len(masterdir)
	else: li = 0
	li                                  = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0]
	masterdir							= mpi.mpi_bcast(masterdir,li,mpi.MPI_CHAR,Blockdata["main_node"],mpi.MPI_COMM_WORLD)
	masterdir                           = string.join(masterdir,"")
	Tracker["constants"]["masterdir"]	= masterdir
	log_main = Logger(BaseLogger_Files())
	log_main.prefix = Tracker["constants"]["masterdir"]+"/"

	while not os.path.exists(Tracker["constants"]["masterdir"]):
		sxprint("Node ", Blockdata["myid"], "  waiting...", Tracker["constants"]["masterdir"])
		sleep(1)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if(Blockdata["myid"] == Blockdata["main_node"]):
		init_dict = {}
		sxprint(Tracker["constants"]["isac_dir"])
		Tracker["directory"] = os.path.join(Tracker["constants"]["isac_dir"], "2dalignment")
		core = read_text_row(os.path.join(Tracker["directory"], "initial2Dparams.txt"))
		for im in range(len(core)): init_dict[im]  = core[im]
		del core
	else: init_dict = 0
	init_dict = wrap_mpi_bcast(init_dict, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)
	###
	do_ctf = True
	if options.noctf: do_ctf = False 
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if do_ctf: sxprint("CTF correction is on")
		else:      sxprint("CTF correction is off")
		if options.local_alignment: sxprint("local refinement is on")
		else:                       sxprint("local refinement is off")
		if B_enhance: sxprint("Bfactor is to be applied on averages")
		elif adjust_to_given_pw2: sxprint("PW of averages is adjusted to a given 1D PW curve")
		elif adjust_to_analytic_model: sxprint("PW of averages is adjusted to analytical model")
		else: sxprint("PW of averages is not adjusted")
		#Tracker["constants"]["orgstack"] = "bdb:"+ os.path.join(Tracker["constants"]["isac_dir"],"../","sparx_stack")
		image = get_im(Tracker["constants"]["orgstack"], 0)
		Tracker["constants"]["nnxo"] = image.get_xsize()
		if Tracker["constants"]["pixel_size"] == -1.0:
			sxprint("Pixel size value is not provided by user. extracting it from ctf header entry of the original stack.")
			try:
				ctf_params = image.get_attr("ctf")
				Tracker["constants"]["pixel_size"] = ctf_params.apix
			except: 
				ERROR( "Pixel size could not be extracted from the original stack.", myid=Blockdata["myid"] )
		## Now fill in low-pass filter
			
		isac_shrink_path = os.path.join(Tracker["constants"]["isac_dir"], "README_shrink_ratio.txt")
		
		if not os.path.exists(isac_shrink_path):
			ERROR( "%s does not exist in the specified ISAC run output directory"%(isac_shrink_path), myid=Blockdata["myid"] )

		isac_shrink_file = open(isac_shrink_path, "r")
		isac_shrink_lines = isac_shrink_file.readlines()
		isac_shrink_ratio = float(isac_shrink_lines[5])  # 6th line: shrink ratio (= [target particle radius]/[particle radius]) used in the ISAC run
		isac_radius = float(isac_shrink_lines[6])        # 7th line: particle radius at original pixel size used in the ISAC run
		isac_shrink_file.close()
		print("Extracted parameter values")
		print("ISAC shrink ratio    : {0}".format(isac_shrink_ratio))
		print("ISAC particle radius : {0}".format(isac_radius))
		Tracker["ini_shrink"] = isac_shrink_ratio
	else: Tracker["ini_shrink"] = 0.0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)

	#print(Tracker["constants"]["pixel_size"], "pixel_size")	
	x_range = max(Tracker["constants"]["xrange"], int(1./Tracker["ini_shrink"]+0.99999) )
	a_range = y_range = x_range

	if(Blockdata["myid"] == Blockdata["main_node"]): parameters = read_text_row(os.path.join(Tracker["constants"]["isac_dir"], "all_parameters.txt"))
	else: parameters = 0
	parameters = wrap_mpi_bcast(parameters, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)		
	params_dict = {}
	list_dict   = {}
	#parepare params_dict

	#navg = min(Tracker["constants"]["navg"]*Blockdata["nproc"], EMUtil.get_image_count(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf")))
	navg = min(Tracker["constants"]["navg"], EMUtil.get_image_count(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf")))
	global_dict = {}
	ptl_list    = []
	memlist     = []
	if(Blockdata["myid"] == Blockdata["main_node"]):
		sxprint("Number of averages computed in this run is %d"%navg)
		for iavg in range(navg):
			params_of_this_average = []
			image   = get_im(os.path.join(Tracker["constants"]["isac_dir"], "class_averages.hdf"), iavg)
			members = sorted(image.get_attr("members"))
			memlist.append(members)
			for im in range(len(members)):
				abs_id =  members[im]
				global_dict[abs_id] = [iavg, im]
				P = combine_params2( init_dict[abs_id][0], init_dict[abs_id][1], init_dict[abs_id][2], init_dict[abs_id][3], \
				parameters[abs_id][0], parameters[abs_id][1]/Tracker["ini_shrink"], parameters[abs_id][2]/Tracker["ini_shrink"], parameters[abs_id][3])
				if parameters[abs_id][3] ==-1: 
					sxprint("WARNING: Image #{0} is an unaccounted particle with invalid 2D alignment parameters and should not be the member of any classes. Please check the consitency of input dataset.".format(abs_id)) # How to check what is wrong about mirror = -1 (Toshio 2018/01/11)
				params_of_this_average.append([P[0], P[1], P[2], P[3], 1.0])
				ptl_list.append(abs_id)
			params_dict[iavg] = params_of_this_average
			list_dict[iavg] = members
			write_text_row(params_of_this_average, os.path.join(Tracker["constants"]["masterdir"], "params_avg", "params_avg_%03d.txt"%iavg))
		ptl_list.sort()
		init_params = [ None for im in range(len(ptl_list))]
		for im in range(len(ptl_list)):
			init_params[im] = [ptl_list[im]] + params_dict[global_dict[ptl_list[im]][0]][global_dict[ptl_list[im]][1]]
		write_text_row(init_params, os.path.join(Tracker["constants"]["masterdir"], "init_isac_params.txt"))
	else:  
		params_dict = 0
		list_dict   = 0
		memlist     = 0
	params_dict = wrap_mpi_bcast(params_dict, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)
	list_dict   = wrap_mpi_bcast(list_dict, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)
	memlist     = wrap_mpi_bcast(memlist, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)
	# Now computing!
	del init_dict
	tag_sharpen_avg = 1000
	## always apply low pass filter to B_enhanced images to suppress noise in high frequencies 
	enforced_to_H1 = False
	if B_enhance:
		if Tracker["constants"]["low_pass_filter"] == -1.0: enforced_to_H1 = True

	# distribute workload among mpi processes
	image_start,image_end = MPI_start_end(navg, Blockdata["nproc"], Blockdata["myid"])

	if Blockdata["myid"] == Blockdata["main_node"]:
		cpu_dict = {}
		for iproc in range(Blockdata["nproc"]):
			local_image_start, local_image_end = MPI_start_end(navg, Blockdata["nproc"], iproc)
			for im in range(local_image_start, local_image_end):
				cpu_dict [im] = iproc
	else:
		cpu_dict = 0
	
	cpu_dict = wrap_mpi_bcast(cpu_dict, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)


	slist      = [None for im in range(navg)]
	ini_list   = [None for im in range(navg)]
	avg1_list  = [None for im in range(navg)]
	avg2_list  = [None for im in range(navg)]
	data_list  = [ None for im in range(navg)]
	plist_dict = {}
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		if B_enhance: 
			sxprint("Avg ID   B-factor  FH1(Res before ali) FH2(Res after ali)")	
		else: 
			sxprint("Avg ID   FH1(Res before ali)  FH2(Res after ali)")

	FH_list = [ [0, 0.0, 0.0] for im in range(navg)]
	for iavg in range(image_start,image_end):

		mlist = EMData.read_images(Tracker["constants"]["orgstack"], list_dict[iavg])

		for im in range(len(mlist)):
			set_params2D(mlist[im], params_dict[iavg][im], xform = "xform.align2d")
		
		if options.local_alignment:
			new_avg, plist, FH2 = refinement_2d_local(mlist, Tracker["constants"]["radius"], a_range,x_range, y_range, CTF = do_ctf, SNR=1.0e10)
			plist_dict[iavg] = plist
			FH1 = -1.0

		else:
			new_avg, frc, plist = compute_average(mlist, Tracker["constants"]["radius"], do_ctf)
			FH1 = get_optimistic_res(frc)
			FH2 = -1.0

		FH_list[iavg] = [iavg, FH1, FH2]
		
		if B_enhance:
			new_avg, gb = apply_enhancement(new_avg, Tracker["constants"]["B_start"], Tracker["constants"]["pixel_size"], Tracker["constants"]["Bfactor"])
			sxprint("  %6d      %6.3f  %4.3f  %4.3f"%(iavg, gb, FH1, FH2))
			
		elif adjust_to_given_pw2: 
			roo = read_text_file( Tracker["constants"]["modelpw"], -1)
			roo = roo[0] # always on the first column
			new_avg = adjust_pw_to_model(new_avg, Tracker["constants"]["pixel_size"], roo)
			sxprint("  %6d      %4.3f  %4.3f  "%(iavg, FH1, FH2))
			
		elif adjust_to_analytic_model:
			new_avg = adjust_pw_to_model(new_avg, Tracker["constants"]["pixel_size"], None)
			sxprint("  %6d      %4.3f  %4.3f   "%(iavg, FH1, FH2))

		elif no_adjustment: pass
		
		if Tracker["constants"]["low_pass_filter"] != -1.0:
			if Tracker["constants"]["low_pass_filter"] == 0.0: low_pass_filter = FH1
			elif Tracker["constants"]["low_pass_filter"] == 1.0: 
				low_pass_filter = FH2
				if not options.local_alignment: low_pass_filter = FH1
			else: 
				low_pass_filter = Tracker["constants"]["low_pass_filter"]
				if low_pass_filter >=0.45: low_pass_filter =0.45 		
			new_avg = filt_tanl(new_avg, low_pass_filter, 0.02)
		else:# No low pass filter but if enforced
			if enforced_to_H1: new_avg = filt_tanl(new_avg, FH1, 0.02)
		if B_enhance: new_avg = fft(new_avg)
			
		new_avg.set_attr("members",   list_dict[iavg])
		new_avg.set_attr("n_objects", len(list_dict[iavg]))
		slist[iavg]    = new_avg
		sxprint(strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>",  "Refined average %7d"%iavg)
		
	## send to main node to write
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	
	for im in range(navg):
		# avg
		if cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] != Blockdata["main_node"]:
			send_EMData(slist[im], Blockdata["main_node"],  tag_sharpen_avg)
		
		elif cpu_dict[im] == Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
			slist[im].set_attr("members", memlist[im])
			slist[im].set_attr("n_objects", len(memlist[im]))
			slist[im].write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
		
		elif cpu_dict[im] != Blockdata["myid"] and Blockdata["myid"] == Blockdata["main_node"]:
			new_avg_other_cpu = recv_EMData(cpu_dict[im], tag_sharpen_avg)
			new_avg_other_cpu.set_attr("members", memlist[im])
			new_avg_other_cpu.set_attr("n_objects", len(memlist[im]))
			new_avg_other_cpu.write_image(os.path.join(Tracker["constants"]["masterdir"], "class_averages.hdf"), im)
		
		if options.local_alignment:
			if cpu_dict[im] == Blockdata["myid"]:
				write_text_row(plist_dict[im], os.path.join(Tracker["constants"]["masterdir"], "ali2d_local_params_avg", "ali2d_local_params_avg_%03d.txt"%im))
			
			if cpu_dict[im] == Blockdata["myid"] and cpu_dict[im]!= Blockdata["main_node"]:
				wrap_mpi_send(plist_dict[im], Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				wrap_mpi_send(FH_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			
			elif cpu_dict[im]!= Blockdata["main_node"] and Blockdata["myid"] == Blockdata["main_node"]:
				dummy = wrap_mpi_recv(cpu_dict[im], mpi.MPI_COMM_WORLD)
				plist_dict[im] = dummy
				dummy = wrap_mpi_recv(cpu_dict[im], mpi.MPI_COMM_WORLD)
				FH_list[im] = dummy[im]
		else:
			if cpu_dict[im] == Blockdata["myid"] and cpu_dict[im]!= Blockdata["main_node"]:
				wrap_mpi_send(FH_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			
			elif cpu_dict[im]!= Blockdata["main_node"] and Blockdata["myid"] == Blockdata["main_node"]:
				dummy = wrap_mpi_recv(cpu_dict[im], mpi.MPI_COMM_WORLD)
				FH_list[im] = dummy[im]
				
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	if options.local_alignment:
		if Blockdata["myid"] == Blockdata["main_node"]:
			ali3d_local_params = [ None for im in range(len(ptl_list)) ]
			for im in range(len(ptl_list)):
				ali3d_local_params[im] = [ptl_list[im]] + plist_dict[global_dict[ptl_list[im]][0]][global_dict[ptl_list[im]][1]]
			write_text_row(ali3d_local_params, os.path.join(Tracker["constants"]["masterdir"], "ali2d_local_params.txt"))
			write_text_row(FH_list, os.path.join(Tracker["constants"]["masterdir"], "FH_list.txt"))
	else:
		if Blockdata["myid"] == Blockdata["main_node"]:
			write_text_row(FH_list, os.path.join(Tracker["constants"]["masterdir"], "FH_list.txt"))
				
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)			
	target_xr =3
	target_yr =3
	if( Blockdata["myid"] == 0):
		cmd = "{} {} {} {} {} {} {} {} {} {}".format("sp_chains.py", os.path.join(Tracker["constants"]["masterdir"],"class_averages.hdf"),\
		os.path.join(Tracker["constants"]["masterdir"],"junk.hdf"),os.path.join(Tracker["constants"]["masterdir"],"ordered_class_averages.hdf"),\
		"--circular","--radius=%d"%Tracker["constants"]["radius"] , "--xr=%d"%(target_xr+1),"--yr=%d"%(target_yr+1),"--align", ">/dev/null")
		junk = cmdexecute(cmd)
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["constants"]["masterdir"], "junk.hdf") )
		junk = cmdexecute(cmd)
		
	return

if __name__ == "__main__":
	sp_global_def.print_timestamp( "Start" )
	main()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()
