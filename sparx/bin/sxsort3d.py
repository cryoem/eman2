#!/usr/bin/env python
#
#
#  08/13/2015
#  New version.
#  
from sparx import *
from EMAN2 import *
import os
import global_def
from   global_def import *
from   optparse  import OptionParser
import sys
from   numpy     import array
import types
from   logger    import Logger, BaseLogger_Files


def main():
	from logger import Logger, BaseLogger_Files
        arglist = []
        i = 0
        while( i < len(sys.argv) ):
            if sys.argv[i]=='-p4pg':
                i = i+2
            elif sys.argv[i]=='-p4wd':
                i = i+2
            else:
                arglist.append( sys.argv[i] )
                i = i+1
	progname = os.path.basename(arglist[0])
	usage = progname + " stack  outdir  <mask> --focus=3Dmask --radius=outer_radius --delta=angular_step" +\
	"--an=angular_neighborhood --maxit=max_iter  --CTF --sym=c1 --function=user_function --independent=indenpendent_runs  --number_of_images_per_group=number_of_images_per_group  --low_pass_frequency=.25  --seed=random_seed"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--focus",                         type="string",        default ='',                   help="bineary 3D mask for focused clustering ")
	parser.add_option("--ir",                            type= "int",          default =1, 	                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",                        type= "int",          default =-1,	                  help="particle radius in pixel for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--maxit",	                     type= "int",          default =25, 	              help="maximum number of iteration")
	parser.add_option("--rs",                            type= "int",          default =1,	                  help="step between rings in rotational correlation >0 (set to 1)" ) 
	parser.add_option("--xr",                            type="string",        default ='1',                  help="range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",                            type="string",        default ='-1',	              help="range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",                            type="string",        default ='0.25',               help="step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",                         type="string",        default ='2',                  help="angular step of reference projections")
	parser.add_option("--an",                            type="string",        default ='-1',	              help="angular neighborhood for local searches")
	parser.add_option("--center",                        type="int",           default =0,	                  help="0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0)")
	parser.add_option("--nassign",                       type="int",           default =1, 	                  help="number of reassignment iterations performed for each angular step (set to 3) ")
	parser.add_option("--nrefine",                       type="int",           default =0, 	                  help="number of alignment iterations performed for each angular step (set to 0)")
	parser.add_option("--CTF",                           action ="store_true", default=False,                 help="do CTF correction during clustring")
	parser.add_option("--stoprnct",                      type="float",         default=3.0,                   help="Minimum percentage of assignment change to stop the program")
	parser.add_option("--sym",                           type="string",        default='c1',                  help="symmetry of the structure ")
	parser.add_option("--function",                      type="string",        default='do_volume_mrk05',     help="name of the reference preparation function")
	parser.add_option("--independent",                   type="int",           default= 3,                    help="number of independent run")
	parser.add_option("--number_of_images_per_group",    type="int",           default=1000,                  help="number of groups")
	parser.add_option("--low_pass_filter",               type="float",         default=-1.0,                  help="absolute frequency of low-pass filter for 3d sorting on the original image size" )
	parser.add_option("--nxinit",                        type="int",           default=64,                    help="initial image size for sorting" )
	parser.add_option("--unaccounted",                   action="store_true",  default=False,                 help="reconstruct the unaccounted images")
	parser.add_option("--seed",                          type="int",           default=-1,                    help="random seed for create initial random assignment for EQ Kmeans")
	parser.add_option("--smallest_group",                type="int",           default=500,                   help="minimum members for identified group")
	parser.add_option("--sausage",                       action="store_true",  default=False,                 help="way of filter volume")
	parser.add_option("--chunkdir",                      type="string",        default='',                    help="chunkdir for computing margin of error")
	parser.add_option("--PWadjustment",                  type="string",        default='',                    help="1-D power spectrum of PDB file used for EM volume power spectrum correction")
	parser.add_option("--upscale",                       type="float",         default=0.5,                   help=" scaling parameter to adjust the power spectrum of EM volumes")
	parser.add_option("--wn",                            type="int",           default=0,                     help="optimal window size for data processing")
	parser.add_option("--interpolation",                 type="string",        default="4nn",                 help="3-d reconstruction interpolation method, two options trl and 4nn")
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1  or len(args) > 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:

		if len(args)>2:
			mask_file = args[2]
		else:
			mask_file = None

		orgstack                        =args[0]
		masterdir                       =args[1]
		global_def.BATCH = True
		#---initialize MPI related variables
		from mpi import mpi_init, mpi_comm_size, MPI_COMM_WORLD, mpi_comm_rank,mpi_barrier,mpi_bcast, mpi_bcast, MPI_INT,MPI_CHAR
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		nproc    = mpi_comm_size(MPI_COMM_WORLD)
		myid     = mpi_comm_rank(MPI_COMM_WORLD)
		mpi_comm = MPI_COMM_WORLD
		main_node= 0
		# import some utilities
		from utilities import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
		from applications import recons3d_n_MPI, mref_ali3d_MPI, Kmref_ali3d_MPI
		from statistics import k_means_match_clusters_asg_new,k_means_stab_bbenum
		from applications import mref_ali3d_EQ_Kmeans, ali3d_mref_Kmeans_MPI  
		# Create the main log file
		from logger import Logger,BaseLogger_Files
		if myid ==main_node:
			log_main=Logger(BaseLogger_Files())
			log_main.prefix = masterdir+"/"
		else:
			log_main =None
		#--- fill input parameters into dictionary named after Constants
		Constants		         ={}
		Constants["stack"]               = args[0]
		Constants["masterdir"]           = masterdir
		Constants["mask3D"]              = mask_file
		Constants["focus3Dmask"]         = options.focus
		Constants["indep_runs"]          = options.independent
		Constants["stoprnct"]            = options.stoprnct
		Constants["number_of_images_per_group"]  = options.number_of_images_per_group
		Constants["CTF"]                 = options.CTF
		Constants["maxit"]               = options.maxit
		Constants["ir"]                  = options.ir 
		Constants["radius"]              = options.radius 
		Constants["nassign"]             = options.nassign
		Constants["rs"]                  = options.rs 
		Constants["xr"]                  = options.xr
		Constants["yr"]                  = options.yr
		Constants["ts"]                  = options.ts
		Constants["delta"]               = options.delta
		Constants["an"]                  = options.an
		Constants["sym"]                 = options.sym
		Constants["center"]              = options.center
		Constants["nrefine"]             = options.nrefine
		#Constants["fourvar"]            = options.fourvar 
		Constants["user_func"]           = options.function
		Constants["low_pass_filter"]     = options.low_pass_filter # enforced low_pass_filter
		#Constants["debug"]              = options.debug
		Constants["main_log_prefix"]     =args[1]
		#Constants["importali3d"]        =options.importali3d
		Constants["myid"]	             = myid
		Constants["main_node"]           =main_node
		Constants["nproc"]               =nproc
		Constants["log_main"]            =log_main
		Constants["nxinit"]              =options.nxinit
		Constants["unaccounted"]         =options.unaccounted
		Constants["seed"]                =options.seed
		Constants["smallest_group"]      =options.smallest_group
		Constants["sausage"]             =options.sausage
		Constants["chunkdir"]            =options.chunkdir
		Constants["PWadjustment"]        =options.PWadjustment
		Constants["upscale"]             =options.upscale
		Constants["wn"]                  =options.wn
		Constants["3d-interpolation"]    =options.interpolation 
		# -----------------------------------------------------
		#
		# Create and initialize Tracker dictionary with input options
		Tracker = 			    		{}
		Tracker["constants"]=				Constants
		Tracker["maxit"]          = Tracker["constants"]["maxit"]
		Tracker["radius"]         = Tracker["constants"]["radius"]
		#Tracker["xr"]             = ""
		#Tracker["yr"]             = "-1"  # Do not change!
		#Tracker["ts"]             = 1
		#Tracker["an"]             = "-1"
		#Tracker["delta"]          = "2.0"
		#Tracker["zoom"]           = True
		#Tracker["nsoft"]          = 0
		#Tracker["local"]          = False
		#Tracker["PWadjustment"]   = Tracker["constants"]["PWadjustment"]
		Tracker["upscale"]        = Tracker["constants"]["upscale"]
		#Tracker["upscale"]        = 0.5
		Tracker["applyctf"]       = False  #  Should the data be premultiplied by the CTF.  Set to False for local continuous.
		#Tracker["refvol"]         = None
		Tracker["nxinit"]         = Tracker["constants"]["nxinit"]
		#Tracker["nxstep"]         = 32
		Tracker["icurrentres"]    = -1
		#Tracker["ireachedres"]    = -1
		#Tracker["lowpass"]        = 0.4
		#Tracker["falloff"]        = 0.2
		#Tracker["inires"]         = options.inires  # Now in A, convert to absolute before using
		Tracker["fuse_freq"]      = 50  # Now in A, convert to absolute before using
		#Tracker["delpreviousmax"] = False
		#Tracker["anger"]          = -1.0
		#Tracker["shifter"]        = -1.0
		#Tracker["saturatecrit"]   = 0.95
		#Tracker["pixercutoff"]    = 2.0
		#Tracker["directory"]      = ""
		#Tracker["previousoutputdir"] = ""
		#Tracker["eliminated-outliers"] = False
		#Tracker["mainiteration"]  = 0
		#Tracker["movedback"]      = False
		#Tracker["state"]          = Tracker["constants"]["states"][0] 
		#Tracker["global_resolution"] =0.0
		Tracker["orgstack"]        = orgstack
		#--------------------------------------------------------------------
		# import from utilities
		from utilities import sample_down_1D_curve,get_initial_ID,remove_small_groups,print_upper_triangular_matrix,print_a_line_with_timestamp
		from utilities import print_dict,get_resolution_mrk01,partition_to_groups,partition_independent_runs,get_outliers
		from utilities import merge_groups, save_alist, margin_of_error, get_margin_of_error, do_two_way_comparison, select_two_runs, get_ali3d_params
		from utilities import counting_projections, unload_dict, load_dict, get_stat_proj, create_random_list, get_number_of_groups, recons_mref
		from utilities import apply_low_pass_filter, get_groups_from_partition, get_number_of_groups, get_complementary_elements_total, update_full_dict
		from utilities import count_chunk_members, set_filter_parameters_from_adjusted_fsc, adjust_fsc_down, get_two_chunks_from_stack
		####------------------------------------------------------------------
		
		
		#
		# Get the pixel size; if none, set to 1.0, and the original image size
		from utilities import get_shrink_data_huang
		if(myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line+"Initialization of 3-D sorting")
			a = get_im(orgstack)
			nnxo = a.get_xsize()
			if( Tracker["nxinit"] > nnxo ):
				ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"sxsort3d.py",1)
				nnxo = -1
			else:
				if Tracker["constants"]["CTF"]:
					i = a.get_attr('ctf')
					pixel_size = i.apix
					fq = pixel_size/Tracker["fuse_freq"]
				else:
					pixel_size = 1.0
					#  No pixel size, fusing computed as 5 Fourier pixels
					fq = 5.0/nnxo
					del a
		else:
			nnxo = 0
			fq = 0.0
			pixel_size = 1.0
		nnxo = bcast_number_to_all(nnxo, source_node = main_node)
		if( nnxo < 0 ):
			mpi_finalize()
			exit()
		pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
		fq         = bcast_number_to_all(fq, source_node = main_node)
		if Tracker["constants"]["wn"]==0:
			Tracker["constants"]["nnxo"]         = nnxo
		else:
			Tracker["constants"]["nnxo"]          = Tracker["constants"]["wn"]
			nnxo     =  Tracker["constants"]["nnxo"]
		Tracker["constants"]["pixel_size"]      = pixel_size
		Tracker["fuse_freq"]    = fq
		del fq, nnxo, pixel_size
		if(Tracker["constants"]["radius"]  < 1):
			Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
		elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
			ERROR("Particle radius set too large!","sxsort3d.py",1,myid)
####-----------------------------------------------------------------------------------------
		# Master directory
		if myid == main_node:
			if masterdir =="":
				timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
				masterdir ="master_sort3d"+timestring
			li =len(masterdir)
			cmd="{} {}".format("mkdir", masterdir)
			os.system(cmd)
		else:
			li=0
		li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
		if li>0:
			masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
			import string
			masterdir = string.join(masterdir,"")
		if myid ==main_node:
			print_dict(Tracker["constants"],"Permanent settings of 3-D sorting program")
		######### create a vstack from input stack to the local stack in masterdir
		# stack name set to default
		Tracker["constants"]["stack"]       = "bdb:"+masterdir+"/rdata"
		Tracker["constants"]["ali3d"]       = os.path.join(masterdir, "ali3d_init.txt")
		Tracker["constants"]["ctf_params"]  = os.path.join(masterdir, "ctf_params.txt")
		Tracker["constants"]["partstack"]   = Tracker["constants"]["ali3d"]  # also serves for refinement
		if myid == main_node:
			total_stack = EMUtil.get_image_count(Tracker["orgstack"])
		else:
			total_stack = 0
		total_stack = bcast_number_to_all(total_stack, source_node = main_node)
		mpi_barrier(MPI_COMM_WORLD)
		from time import sleep
		while not os.path.exists(masterdir):
				print  "Node ",myid,"  waiting..."
				sleep(5)
		mpi_barrier(MPI_COMM_WORLD)
		if myid == main_node:
			log_main.add("sort3d starts in "+masterdir)
		#####
		###----------------------------------------------------------------------------------
		# Initial data analysis and handle two chunk files
		from random import shuffle
		# Compute the resolution 
		#### make chunkdir dictionary for computing margin of error
		import user_functions
		user_func  = user_functions.factory[Tracker["constants"]["user_func"]]
		chunk_dict = {}
		chunk_list = []
		if myid == main_node:
			chunk_one = read_text_file(os.path.join(Tracker["constants"]["chunkdir"],"chunk0.txt"))
			chunk_two = read_text_file(os.path.join(Tracker["constants"]["chunkdir"],"chunk1.txt"))
		else:
			chunk_one = 0
			chunk_two = 0
		chunk_one = wrap_mpi_bcast(chunk_one, main_node)
		chunk_two = wrap_mpi_bcast(chunk_two, main_node)
		mpi_barrier(MPI_COMM_WORLD)
		######################## Read/write bdb: data on main node ############################
	   	if myid==main_node:
			if(orgstack[:4] == "bdb:"):	cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+Tracker["constants"]["stack"])
			else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, Tracker["constants"]["stack"])
	   		cmdexecute(cmd)
			cmd = "{} {} {}".format("sxheader.py  --params=xform.projection", "--export="+Tracker["constants"]["ali3d"],orgstack)
			cmdexecute(cmd)
			cmd = "{} {} {}".format("sxheader.py  --params=ctf", "--export="+Tracker["constants"]["ctf_params"],orgstack)
			cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)	   		   	
		########-----------------------------------------------------------------------------
		Tracker["total_stack"]              = total_stack
		Tracker["constants"]["total_stack"] = total_stack
		Tracker["shrinkage"]                = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
		Tracker["radius"]                   = Tracker["constants"]["radius"]*Tracker["shrinkage"]
		if Tracker["constants"]["mask3D"]:
			Tracker["mask3D"] = os.path.join(masterdir,"smask.hdf")
		else:
			Tracker["mask3D"]  = None
		if Tracker["constants"]["focus3Dmask"]:
			Tracker["focus3D"] = os.path.join(masterdir,"sfocus.hdf")
		else:
			Tracker["focus3D"] = None
		if myid == main_node:
			if Tracker["constants"]["mask3D"]:
				mask_3D = get_shrink_3dmask(Tracker["nxinit"],Tracker["constants"]["mask3D"])
				mask_3D.write_image(Tracker["mask3D"])
			if Tracker["constants"]["focus3Dmask"]:
				mask_3D = get_shrink_3dmask(Tracker["nxinit"],Tracker["constants"]["focus3Dmask"])
				st = Util.infomask(mask_3D, None, True)
				if( st[0] == 0.0 ):  ERROR("sxrsort3d","incorrect focused mask, after binarize all values zero",1)
				mask_3D.write_image(Tracker["focus3D"])
				del mask_3D
		if Tracker["constants"]["PWadjustment"] !='':
			PW_dict              = {}
			nxinit_pwsp          = sample_down_1D_curve(Tracker["constants"]["nxinit"],Tracker["constants"]["nnxo"],Tracker["constants"]["PWadjustment"])
			Tracker["nxinit_PW"] = os.path.join(masterdir,"spwp.txt")
			if myid == main_node:  write_text_file(nxinit_pwsp,Tracker["nxinit_PW"])
			PW_dict[Tracker["constants"]["nnxo"]]   = Tracker["constants"]["PWadjustment"]
			PW_dict[Tracker["constants"]["nxinit"]] = Tracker["nxinit_PW"]
			Tracker["PW_dict"]                      = PW_dict
		mpi_barrier(MPI_COMM_WORLD)
		#-----------------------From two chunks to FSC, and low pass filter-----------------------------------------###
		for element in chunk_one: chunk_dict[element] = 0
		for element in chunk_two: chunk_dict[element] = 1
		chunk_list =[chunk_one, chunk_two]
		Tracker["chunk_dict"] = chunk_dict
		Tracker["P_chunk0"]   = len(chunk_one)/float(total_stack)
		Tracker["P_chunk1"]   = len(chunk_two)/float(total_stack)
		### create two volumes to estimate resolution
		if myid == main_node:
			for index in xrange(2): write_text_file(chunk_list[index],os.path.join(masterdir,"chunk%01d.txt"%index))
		mpi_barrier(MPI_COMM_WORLD)
		vols = []
		for index in xrange(2):
			data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"], os.path.join(masterdir,"chunk%01d.txt"%index), Tracker["constants"]["partstack"],myid,main_node,nproc,preshift=True)
			vol             = recons3d_4nn_ctf_MPI(myid=myid, prjlist=data,symmetry=Tracker["constants"]["sym"], finfo=None)
			if myid == main_node:
				vol.write_image(os.path.join(masterdir, "vol%d.hdf"%index))
			vols.append(vol)
			mpi_barrier(MPI_COMM_WORLD)
		if myid ==main_node:
			low_pass, falloff,currentres = get_resolution_mrk01(vols,Tracker["constants"]["radius"],Tracker["constants"]["nxinit"],masterdir,Tracker["mask3D"])
			if low_pass >Tracker["constants"]["low_pass_filter"]: low_pass= Tracker["constants"]["low_pass_filter"]
		else:
			low_pass    =0.0
			falloff     =0.0
			currentres  =0.0
		bcast_number_to_all(currentres,source_node = main_node)
		bcast_number_to_all(low_pass,source_node   = main_node)
		bcast_number_to_all(falloff,source_node    = main_node)
		Tracker["currentres"]                      = currentres
		Tracker["falloff"]                         = falloff
		if Tracker["constants"]["low_pass_filter"] ==-1.0:
			Tracker["low_pass_filter"] = min(.45,low_pass/Tracker["shrinkage"]) # no better than .45
		else:
			Tracker["low_pass_filter"] = min(.45,Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"])
		Tracker["lowpass"]             = Tracker["low_pass_filter"]
		Tracker["falloff"]             =.1
		Tracker["global_fsc"]          = os.path.join(masterdir, "fsc.txt")
		############################################################################################
		if myid == main_node:
			log_main.add("The command-line inputs are as following:")
			log_main.add("**********************************************************")
		for a in sys.argv:
			if myid == main_node:log_main.add(a)
		if myid == main_node:
			log_main.add("number of cpus used in this run is %d"%Tracker["constants"]["nproc"])
			log_main.add("**********************************************************")
		from filter import filt_tanl
		### START 3-D sorting
		if myid ==main_node:
			log_main.add("----------3-D sorting  program------- ")
			log_main.add("current resolution %6.3f for images of original size in terms of absolute frequency"%Tracker["currentres"])
			log_main.add("equivalent to %f Angstrom resolution"%(Tracker["constants"]["pixel_size"]/Tracker["currentres"]/Tracker["shrinkage"]))
			log_main.add("the user provided enforced low_pass_filter is %f"%Tracker["constants"]["low_pass_filter"])
			#log_main.add("equivalent to %f Angstrom resolution"%(Tracker["constants"]["pixel_size"]/Tracker["constants"]["low_pass_filter"]))
			for index in xrange(2):
				filt_tanl(get_im(os.path.join(masterdir,"vol%01d.hdf"%index)), Tracker["low_pass_filter"],Tracker["falloff"]).write_image(os.path.join(masterdir, "volf%01d.hdf"%index))
		mpi_barrier(MPI_COMM_WORLD)
		from utilities import get_input_from_string
		delta       = get_input_from_string(Tracker["constants"]["delta"])
		delta       = delta[0]
		from utilities import even_angles
		n_angles    = even_angles(delta, 0, 180)
		this_ali3d  = Tracker["constants"]["ali3d"]
		sampled     = get_stat_proj(Tracker,delta,this_ali3d)
		if myid ==main_node:
			nc = 0
			for a in sampled:
				if len(sampled[a])>0:
					nc += 1
			log_main.add("total sampled direction %10d  at angle step %6.3f"%(len(n_angles), delta)) 
			log_main.add("captured sampled directions %10d percentage covered by data  %6.3f"%(nc,float(nc)/len(n_angles)*100))
		number_of_images_per_group = Tracker["constants"]["number_of_images_per_group"]
		if myid ==main_node: log_main.add("user provided number_of_images_per_group %d"%number_of_images_per_group)
		Tracker["number_of_images_per_group"] = number_of_images_per_group
		number_of_groups = get_number_of_groups(total_stack,number_of_images_per_group)
		Tracker["number_of_groups"] =  number_of_groups
		generation     =0
		partition_dict ={}
		full_dict      ={}
		workdir =os.path.join(masterdir,"generation%03d"%generation)
		Tracker["this_dir"] = workdir
		if myid ==main_node:
			log_main.add("---- generation         %5d"%generation)
			log_main.add("number of images per group is set as %d"%number_of_images_per_group)
			log_main.add("the initial number of groups is  %10d "%number_of_groups)
			cmd="{} {}".format("mkdir",workdir)
			os.system(cmd)
		mpi_barrier(MPI_COMM_WORLD)
		list_to_be_processed = range(Tracker["constants"]["total_stack"])
		Tracker["this_data_list"] = list_to_be_processed
		create_random_list(Tracker)
		#################################
		full_dict ={}
		for iptl in xrange(Tracker["constants"]["total_stack"]):
			 full_dict[iptl]    = iptl
		Tracker["full_ID_dict"] = full_dict
		################################# 	
		for indep_run in xrange(Tracker["constants"]["indep_runs"]):
			Tracker["this_particle_list"] = Tracker["this_indep_list"][indep_run]
			ref_vol =  recons_mref(Tracker)
			if myid == main_node: log_main.add("independent run  %10d"%indep_run)
			mpi_barrier(MPI_COMM_WORLD)
			Tracker["this_data_list"]          = list_to_be_processed
			Tracker["total_stack"]             = len(Tracker["this_data_list"])
			Tracker["this_particle_text_file"] = os.path.join(workdir,"independent_list_%03d.txt"%indep_run) # for get_shrink_data
			if myid == main_node: write_text_file(Tracker["this_data_list"], Tracker["this_particle_text_file"])
			mpi_barrier(MPI_COMM_WORLD)
			outdir  = os.path.join(workdir, "EQ_Kmeans%03d"%indep_run)
			ref_vol = apply_low_pass_filter(ref_vol,Tracker)
			mref_ali3d_EQ_Kmeans(ref_vol, outdir, Tracker["this_particle_text_file"], Tracker)
			partition_dict[indep_run]=Tracker["this_partition"]
		Tracker["partition_dict"]    = partition_dict
		Tracker["total_stack"]       = len(Tracker["this_data_list"])
		Tracker["this_total_stack"]  = Tracker["total_stack"]
		###############################
		do_two_way_comparison(Tracker)
		###############################
		ref_vol_list = []
		from time import sleep
		number_of_ref_class = []
		for igrp in xrange(len(Tracker["two_way_stable_member"])):
			Tracker["this_data_list"]      = Tracker["two_way_stable_member"][igrp]
			Tracker["this_data_list_file"] = os.path.join(workdir,"stable_class%d.txt"%igrp)
			if myid == main_node:
				write_text_file(Tracker["this_data_list"], Tracker["this_data_list_file"])
			data,old_shifts = get_shrink_data_huang(Tracker,Tracker["nxinit"], Tracker["this_data_list_file"], Tracker["constants"]["partstack"], myid, main_node, nproc, preshift = True)
			volref          = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo = None)
			ref_vol_list.append(volref)
			number_of_ref_class.append(len(Tracker["this_data_list"]))
			if myid == main_node:
				log_main.add("group  %d  members %d "%(igrp,len(Tracker["this_data_list"])))
		Tracker["number_of_ref_class"] = number_of_ref_class
		nx_of_image = ref_vol_list[0].get_xsize()
		if Tracker["constants"]["PWadjustment"]:
			Tracker["PWadjustment"] = Tracker["PW_dict"][nx_of_image]
		else:
			Tracker["PWadjustment"] = Tracker["constants"]["PWadjustment"]	 # no PW adjustment
		if myid == main_node:
			for iref in xrange(len(ref_vol_list)):
				refdata    = [None]*4
				refdata[0] = ref_vol_list[iref]
				refdata[1] = Tracker
				refdata[2] = Tracker["constants"]["myid"]
				refdata[3] = Tracker["constants"]["nproc"]
				volref     = user_func(refdata)
				volref.write_image(os.path.join(workdir,"volf_stable.hdf"),iref)
		mpi_barrier(MPI_COMM_WORLD)
		Tracker["this_data_list"]           = Tracker["this_accounted_list"]
		outdir                              = os.path.join(workdir,"Kmref")  
		empty_group, res_groups, final_list = ali3d_mref_Kmeans_MPI(ref_vol_list,outdir,Tracker["this_accounted_text"],Tracker)
		Tracker["this_unaccounted_list"]    = get_complementary_elements(list_to_be_processed,final_list)
		if myid == main_node:
			log_main.add("the number of particles not processed is %d"%len(Tracker["this_unaccounted_list"]))
			write_text_file(Tracker["this_unaccounted_list"],Tracker["this_unaccounted_text"])
		update_full_dict(Tracker["this_unaccounted_list"], Tracker)
		#######################################
		number_of_groups    = len(res_groups)
		vol_list            = []
		number_of_ref_class = []
		for igrp in xrange(number_of_groups):
			data,old_shifts = get_shrink_data_huang(Tracker, Tracker["constants"]["nnxo"], os.path.join(outdir,"Class%d.txt"%igrp), Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
			volref          = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo=None)
			vol_list.append(volref)

			if( myid == main_node ):  npergroup = len(read_text_file(os.path.join(outdir,"Class%d.txt"%igrp)))
			else:  npergroup = 0
			npergroup = bcast_number_to_all(npergroup, main_node )
			number_of_ref_class.append(npergroup)

		Tracker["number_of_ref_class"] = number_of_ref_class
		
		mpi_barrier(MPI_COMM_WORLD)
		nx_of_image = vol_list[0].get_xsize()
		if Tracker["constants"]["PWadjustment"]:
			Tracker["PWadjustment"]=Tracker["PW_dict"][nx_of_image]
		else:
			Tracker["PWadjustment"]=Tracker["constants"]["PWadjustment"]	

		if myid == main_node:
			for ivol in xrange(len(vol_list)):
				refdata     =[None]*4
				refdata[0] = vol_list[ivol]
				refdata[1] = Tracker
				refdata[2] = Tracker["constants"]["myid"]
				refdata[3] = Tracker["constants"]["nproc"] 
				volref = user_func(refdata)
				volref.write_image(os.path.join(workdir,"volf_of_Classes.hdf"),ivol)
				log_main.add("number of unaccounted particles  %10d"%len(Tracker["this_unaccounted_list"]))
				log_main.add("number of accounted particles  %10d"%len(Tracker["this_accounted_list"]))
				
		Tracker["this_data_list"]    = Tracker["this_unaccounted_list"]   # reset parameters for the next round calculation
		Tracker["total_stack"]       = len(Tracker["this_unaccounted_list"])
		Tracker["this_total_stack"]  = Tracker["total_stack"]
		number_of_groups             = get_number_of_groups(len(Tracker["this_unaccounted_list"]),number_of_images_per_group)
		Tracker["number_of_groups"]  =  number_of_groups
		while number_of_groups >= 2 :
			generation     +=1
			partition_dict ={}
			workdir =os.path.join(masterdir,"generation%03d"%generation)
			Tracker["this_dir"] = workdir
			if myid ==main_node:
				log_main.add("*********************************************")
				log_main.add("-----    generation             %5d    "%generation)
				log_main.add("number of images per group is set as %10d "%number_of_images_per_group)
				log_main.add("the number of groups is  %10d "%number_of_groups)
				log_main.add(" number of particles for clustering is %10d"%Tracker["total_stack"])
				cmd ="{} {}".format("mkdir",workdir)
				os.system(cmd)
			mpi_barrier(MPI_COMM_WORLD)
			create_random_list(Tracker)
			for indep_run in xrange(Tracker["constants"]["indep_runs"]):
				Tracker["this_particle_list"] = Tracker["this_indep_list"][indep_run]
				ref_vol                       = recons_mref(Tracker)
				if myid == main_node:
					log_main.add("independent run  %10d"%indep_run)
					outdir = os.path.join(workdir, "EQ_Kmeans%03d"%indep_run)
				Tracker["this_data_list"]   = Tracker["this_unaccounted_list"]
				#ref_vol=apply_low_pass_filter(ref_vol,Tracker)
				mref_ali3d_EQ_Kmeans(ref_vol,outdir,Tracker["this_unaccounted_text"],Tracker)
				partition_dict[indep_run]   = Tracker["this_partition"]
				Tracker["this_data_list"]   = Tracker["this_unaccounted_list"]
				Tracker["total_stack"]      = len(Tracker["this_unaccounted_list"])
				Tracker["partition_dict"]   = partition_dict
				Tracker["this_total_stack"] = Tracker["total_stack"]
			total_list_of_this_run          = Tracker["this_unaccounted_list"]
			###############################
			do_two_way_comparison(Tracker)
			###############################
			ref_vol_list        = []
			number_of_ref_class = []
			for igrp in xrange(len(Tracker["two_way_stable_member"])):
				Tracker["this_data_list"]      = Tracker["two_way_stable_member"][igrp]
				Tracker["this_data_list_file"] = os.path.join(workdir,"stable_class%d.txt"%igrp)
				if myid == main_node: write_text_file(Tracker["this_data_list"], Tracker["this_data_list_file"])
				mpi_barrier(MPI_COMM_WORLD)
				data,old_shifts  = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"],Tracker["this_data_list_file"],Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
				volref           = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"],finfo= None)
				#volref = filt_tanl(volref, Tracker["constants"]["low_pass_filter"],.1)
				if myid == main_node:volref.write_image(os.path.join(workdir,"vol_stable.hdf"),iref)
				#volref = resample(volref,Tracker["shrinkage"])
				ref_vol_list.append(volref)
				number_of_ref_class.append(len(Tracker["this_data_list"]))
				mpi_barrier(MPI_COMM_WORLD)
			Tracker["number_of_ref_class"]      = number_of_ref_class
			Tracker["this_data_list"]           = Tracker["this_accounted_list"]
			outdir                              = os.path.join(workdir,"Kmref")
			empty_group, res_groups, final_list = ali3d_mref_Kmeans_MPI(ref_vol_list,outdir,Tracker["this_accounted_text"],Tracker)
			# calculate the 3-D structure of original image size for each group
			number_of_groups                    =  len(res_groups)
			Tracker["this_unaccounted_list"]    = get_complementary_elements(total_list_of_this_run,final_list)
			if myid == main_node:
				log_main.add("the number of particles not processed is %d"%len(Tracker["this_unaccounted_list"]))
				write_text_file(Tracker["this_unaccounted_list"],Tracker["this_unaccounted_text"])
			mpi_barrier(MPI_COMM_WORLD)
			update_full_dict(Tracker["this_unaccounted_list"],Tracker)
			vol_list = []
			for igrp in xrange(number_of_groups):
				data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"], os.path.join(outdir,"Class%d.txt"%igrp), Tracker["constants"]["partstack"], myid, main_node, nproc,preshift = True)
				volref = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"],finfo= None)
				vol_list.append(volref)

			mpi_barrier(MPI_COMM_WORLD)
			nx_of_image=ref_vol_list[0].get_xsize()
			if Tracker["constants"]["PWadjustment"]:
				Tracker["PWadjustment"] = Tracker["PW_dict"][nx_of_image]
			else:
				Tracker["PWadjustment"] = Tracker["constants"]["PWadjustment"]	

			if myid == main_node:
				for ivol in xrange(len(vol_list)):
					refdata    = [None]*4
					refdata[0] = vol_list[ivol]
					refdata[1] = Tracker
					refdata[2] = Tracker["constants"]["myid"]
					refdata[3] = Tracker["constants"]["nproc"] 
					volref = user_func(refdata)
					volref.write_image(os.path.join(workdir, "volf_of_Classes.hdf"),ivol)
				log_main.add("number of unaccounted particles  %10d"%len(Tracker["this_unaccounted_list"]))
				log_main.add("number of accounted particles  %10d"%len(Tracker["this_accounted_list"]))
			del vol_list
			mpi_barrier(MPI_COMM_WORLD)
			number_of_groups            = get_number_of_groups(len(Tracker["this_unaccounted_list"]),number_of_images_per_group)
			Tracker["number_of_groups"] =  number_of_groups
			Tracker["this_data_list"]   = Tracker["this_unaccounted_list"]
			Tracker["total_stack"]      = len(Tracker["this_unaccounted_list"])
		if Tracker["constants"]["unaccounted"]:
			data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"],Tracker["this_unaccounted_text"],Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
			volref          = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"],finfo= None)
			nx_of_image     = volref.get_xsize()
			if Tracker["constants"]["PWadjustment"]:
				Tracker["PWadjustment"]=Tracker["PW_dict"][nx_of_image]
			else:
				Tracker["PWadjustment"]=Tracker["constants"]["PWadjustment"]	
			if( myid == main_node ):
				refdata    = [None]*4
				refdata[0] = volref
				refdata[1] = Tracker
				refdata[2] = Tracker["constants"]["myid"]
				refdata[3] =  Tracker["constants"]["nproc"]
				volref     = user_func(refdata)
				#volref = filt_tanl(volref, Tracker["constants"]["low_pass_filter"],.1)
				volref.write_image(os.path.join(workdir,"volf_unaccounted.hdf"))
		# Finish program
		if myid ==main_node: log_main.add("sxsort3d finishes")
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
if __name__ == "__main__":
	main()
