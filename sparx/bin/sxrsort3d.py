#!/usr/bin/env python
#
#  10/25/2015
#  New version.
#  
from sparx import *
import os
import global_def
from   global_def import *
from   optparse  import OptionParser
import sys
from   numpy     import array
import types
from   logger    import Logger, BaseLogger_Files
from    morphology import 	get_shrink_3dmask

	
def main():
	from time import sleep
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
	parser.add_option("--focus",                         type="string",               default=None,             help="3D mask for focused clustering ")
	parser.add_option("--ir",                            type= "int",                 default= 1, 	            help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",                        type= "int",                 default=-1,	            help="outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--maxit",	                     type= "int",                 default=50, 	            help="maximum number of iteration")
	parser.add_option("--rs",                            type= "int",                 default=1,	            help="step between rings in rotational correlation >0 (set to 1)" ) 
	parser.add_option("--xr",                            type="string",               default='1',              help="range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",                            type="string",               default='-1',	            help="range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",                            type="string",               default='0.25',           help="step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",                         type="string",               default='2',              help="angular step of reference projections")
	parser.add_option("--an",                            type="string",               default='-1',	            help="angular neighborhood for local searches")
	parser.add_option("--center",                        type="int",                  default=0,	            help="0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0)")
	parser.add_option("--nassign",                       type="int",                  default=1, 	            help="number of reassignment iterations performed for each angular step (set to 3) ")
	parser.add_option("--nrefine",                       type="int",                  default=0, 	            help="number of alignment iterations performed for each angular step (set to 1) ")
	parser.add_option("--CTF",                           action="store_true",         default=False,             help="Consider CTF correction during the alignment ")
	parser.add_option("--stoprnct",                      type="float",                default=3.0,               help="Minimum percentage of assignment change to stop the program")
	parser.add_option("--sym",                           type="string",               default='c1',              help="symmetry of the structure ")
	parser.add_option("--function",                      type="string",               default='do_volume_mrk05', help="name of the reference preparation function")
	parser.add_option("--independent",                   type="int",                  default= 3,                help="number of independent run")
	parser.add_option("--number_of_images_per_group",    type='int',                  default=1000,        help="number of images per groups")
	parser.add_option("--low_pass_filter",               type="float",                default=-1.0,        help="absolute frequency of low-pass filter for 3d sorting on the original image size" )
	parser.add_option("--nxinit",                        type="int",                  default=64,          help="initial image size for sorting" )
	parser.add_option("--unaccounted",                   action="store_true",         default=False,       help="reconstruct the unaccounted images")
	parser.add_option("--seed",                          type="int",                  default=-1,          help="random seed for create initial random assignment for EQ Kmeans")
	parser.add_option("--smallest_group",                type="int",                  default=500,         help="minimum members for identified group" )
	parser.add_option("--previous_run1",                 type="string",               default='',          help="two previous runs" )
	parser.add_option("--previous_run2",                 type="string",               default='',          help="two previous runs" )
	parser.add_option("--group_size_for_unaccounted",    type="int",                  default=500,         help="size for unaccounted particles" )
	parser.add_option("--chunkdir",                      type="string",               default='',          help="chunkdir for computing margin of error")
	parser.add_option("--sausage",                       action="store_true",         default=False,       help="way of filter volume")
	parser.add_option("--PWadjustment",                  type="string",               default='',          help="1-D power spectrum of PDB file used for EM volume power spectrum correction")
	parser.add_option("--upscale",                       type="float",                default=0.5,         help=" scaling parameter to adjust the power spectrum of EM volumes")
	parser.add_option("--wn",                            type="int",                  default=0,           help="optimal window size for data processing")
	parser.add_option("--interpolation",                 type="string",               default="4nn",       help="3-d reconstruction interpolation method, two options, trl and 4nn")

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
		from mpi import mpi_init, mpi_comm_size, MPI_COMM_WORLD, mpi_comm_rank,mpi_barrier,mpi_bcast, mpi_bcast, MPI_INT
		sys.argv  = mpi_init(len(sys.argv),sys.argv)
		nproc     = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		mpi_comm  = MPI_COMM_WORLD
		main_node = 0
		# import some utilities
		from utilities import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file,wrap_mpi_bcast
		from applications import recons3d_n_MPI, mref_ali3d_MPI, Kmref_ali3d_MPI
		from statistics import k_means_match_clusters_asg_new,k_means_stab_bbenum
		from reconstruction import rec3D_MPI_noCTF,rec3D_two_chunks_MPI
		from applications import mref_ali3d_EQ_Kmeans, ali3d_mref_Kmeans_MPI  
		# Create the main log file
		from logger import Logger,BaseLogger_Files
		if myid ==main_node:
			log_main=Logger(BaseLogger_Files())
			log_main.prefix=masterdir+"/"
		else:
			log_main = None
		#--- fill input parameters into dictionary named after Constants
		Constants		                         ={}
		Constants["stack"]                       =args[0]
		Constants["masterdir"]                   =masterdir
		Constants["mask3D"]                      =mask_file
		Constants["focus3Dmask"]                 =options.focus
		Constants["indep_runs"]                  =options.independent
		Constants["stoprnct"]                    =options.stoprnct
		Constants["number_of_images_per_group"]  =options.number_of_images_per_group
		Constants["CTF"]                 =options.CTF
		Constants["maxit"]               =options.maxit
		Constants["ir"]                  =options.ir 
		Constants["radius"]              =options.radius 
		Constants["nassign"]             =options.nassign
		Constants["rs"]                  =options.rs 
		Constants["xr"]                  =options.xr
		Constants["yr"]                  =options.yr
		Constants["ts"]                  =options.ts
		Constants["delta"]               =options.delta
		Constants["an"]                  =options.an
		Constants["sym"]                 =options.sym
		Constants["center"]              =options.center
		Constants["nrefine"]             =options.nrefine
		Constants["user_func"]           =options.function
		Constants["low_pass_filter"]     =options.low_pass_filter # enforced low_pass_filter
		Constants["main_log_prefix"]     =args[1]
		#Constants["importali3d"]        =options.importali3d
		Constants["myid"]	             =myid
		Constants["main_node"]           =main_node
		Constants["nproc"]               =nproc
		Constants["log_main"]            =log_main
		Constants["nxinit"]              =options.nxinit
		Constants["unaccounted"]         =options.unaccounted
		Constants["seed"]                =options.seed
		Constants["smallest_group"]      =options.smallest_group
		Constants["previous_runs"]       =options.previous_run1+" "+options.previous_run2
		Constants["sausage"]             =options.sausage
		Constants["chunkdir"]            =options.chunkdir 
		Constants["PWadjustment"]        =options.PWadjustment
		Constants["upscale"]             =options.upscale
		Constants["wn"]                  =options.wn
		Constants["3d-interpolation"]    =options.interpolation  
		#Constants["frequency_stop_search"] = options.frequency_stop_search
		#Constants["scale_of_number"]    = options.scale_of_number
		# -------------------------------------------------------------
		#
		# Create and initialize Tracker dictionary with input options
		Tracker                   = {}
		Tracker["constants"]      =	Constants
		Tracker["maxit"]          = Tracker["constants"]["maxit"]
		Tracker["radius"]         = Tracker["constants"]["radius"]
		#Tracker["xr"]            = ""
		#Tracker["yr"]            = "-1"  # Do not change!
		#Tracker["ts"]            = 1
		#Tracker["an"]            = "-1"
		#Tracker["delta"]         = "2.0"
		#Tracker["zoom"]          = True
		#Tracker["nsoft"]         = 0
		#Tracker["local"]         = False
		Tracker["PWadjustment"]   = Tracker["constants"]["PWadjustment"]
		Tracker["upscale"]        = Tracker["constants"]["upscale"]
		Tracker["applyctf"]       = False  #  Should the data be premultiplied by the CTF.  Set to False for local continuous.
		#Tracker["refvol"]        = None
		Tracker["nxinit"]         = Tracker["constants"]["nxinit"]
		#Tracker["nxstep"]        = 32
		Tracker["icurrentres"]    = -1
		#Tracker["ireachedres"]   = -1
		Tracker["lowpass"]        = Tracker["constants"]["low_pass_filter"]
		Tracker["falloff"]        = 0.1
		#Tracker["inires"]        = options.inires  # Now in A, convert to absolute before using
		Tracker["fuse_freq"]      = 50  # Now in A, convert to absolute before using
		#Tracker["delpreviousmax"]= False
		#Tracker["anger"]         = -1.0
		#Tracker["shifter"]       = -1.0
		#Tracker["saturatecrit"]  = 0.95
		#Tracker["pixercutoff"]   = 2.0
		#Tracker["directory"]     = ""
		#Tracker["previousoutputdir"] = ""
		#Tracker["eliminated-outliers"] = False
		#Tracker["mainiteration"]       = 0
		#Tracker["movedback"]           = False
		#Tracker["state"]               = Tracker["constants"]["states"][0] 
		#Tracker["global_resolution"]   = 0.0
		Tracker["orgstack"]             = orgstack
		#--------------------------------------------------------------------
		
		# import from utilities
		from utilities import sample_down_1D_curve,get_initial_ID,remove_small_groups,print_upper_triangular_matrix,print_a_line_with_timestamp
		from utilities import convertasi,prepare_ptp,print_dict,get_resolution_mrk01,partition_to_groups,partition_independent_runs,get_outliers
		from utilities import merge_groups, save_alist, margin_of_error, get_margin_of_error, do_two_way_comparison, select_two_runs, get_ali3d_params
		from utilities import counting_projections, unload_dict, load_dict, get_stat_proj, create_random_list, get_number_of_groups, recons_mref
		from utilities import apply_low_pass_filter, get_groups_from_partition, get_number_of_groups, get_complementary_elements_total, update_full_dict
		from utilities import count_chunk_members, set_filter_parameters_from_adjusted_fsc, adjust_fsc_down, get_two_chunks_from_stack
		####------------------------------------------------------------------	
		
		# another part
		from utilities import get_class_members, remove_small_groups, get_number_of_groups, get_stable_members_from_two_runs
		from utilities import two_way_comparison_single, get_leftover_from_stable, get_initial_ID, Kmeans_exhaustive_run
		from utilities import print_a_line_with_timestamp, split_a_group
		
		#
		# Get the pixel size; if none, set to 1.0, and the original image size
		from utilities import get_shrink_data_huang
		from time import sleep
		import user_functions
		user_func = user_functions.factory[Tracker["constants"]["user_func"]]
		if(myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line+"Initialization of 3-D sorting")
			a = get_im(Tracker["orgstack"])
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
			fq   = 0.0
			pixel_size = 1.0
		nnxo = bcast_number_to_all(nnxo, source_node = main_node)
		if( nnxo < 0 ):
			mpi_finalize()
			exit()
		pixel_size                           = bcast_number_to_all(pixel_size, source_node = main_node)
		fq                                   = bcast_number_to_all(fq, source_node = main_node)
		if Tracker["constants"]["wn"]==0:
			Tracker["constants"]["nnxo"] = nnxo
		else:
			Tracker["constants"]["nnxo"] = Tracker["constants"]["wn"]
			nnxo= Tracker["constants"]["wn"]
		Tracker["constants"]["pixel_size"]   = pixel_size
		Tracker["fuse_freq"]                 = fq
		del fq, nnxo, pixel_size
		if(Tracker["constants"]["radius"]  < 1):
			Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
		elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
			ERROR("Particle radius set too large!","sxsort3d.py",1,myid)
####-----------------------------------------------------------------------------------------
		# create the master directory
		if myid == main_node:
			if masterdir =="":
				timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
				masterdir ="master_sort3d"+timestring
				li =len(masterdir)
			else:
				li = 0
			cmd="{} {}".format("mkdir", masterdir)
			os.system(cmd)			
		else:
			li=0
		li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
		if li>0:
			masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
			masterdir = string.join(masterdir,"")
		####--- masterdir done!
		if myid ==main_node:
			print_dict(Tracker["constants"],"Permanent settings of 3-D sorting program")
		from time import sleep
		while not os.path.exists(masterdir):  # Be sure each proc is able to access the created dir
				print  "Node ",myid,"  waiting..."
				sleep(5)
		mpi_barrier(MPI_COMM_WORLD)
		######### create a vstack from input stack to the local stack in masterdir
		# stack name set to default
		Tracker["constants"]["stack"]          = "bdb:"+masterdir+"/rdata"
		Tracker["constants"]["ali3d"]          = os.path.join(masterdir, "ali3d_init.txt")
		Tracker["constants"]["partstack"]      = Tracker["constants"]["ali3d"]
		Tracker["constants"]["ctf_params"]     = os.path.join(masterdir, "ctf_params.txt")
		######
	   	if myid == main_node:
			if(Tracker["orgstack"][:4] == "bdb:"):     cmd = "{} {} {}".format("e2bdb.py", Tracker["orgstack"],"--makevstack="+Tracker["constants"]["stack"])
			else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, Tracker["constants"]["stack"])
			cmdexecute(cmd)
			cmd = "{} {} {} {} ".format("sxheader.py", Tracker["constants"]["stack"],"--params=xform.projection","--export="+Tracker["constants"]["ali3d"])
			cmdexecute(cmd)
			cmd = "{} {} {} {} ".format("sxheader.py", Tracker["constants"]["stack"],"--params=ctf","--export="+Tracker["constants"]["ctf_params"])
			cmdexecute(cmd)
			#keepchecking = False
			total_stack = EMUtil.get_image_count(Tracker["orgstack"])
		else:
			total_stack =0
	   	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
	   	"""
		if myid==main_node:
	   		from EMAN2db import db_open_dict	
	   		OB = db_open_dict(orgstack)
	   		DB = db_open_dict(Tracker["constants"]["stack"]) 
			for i in xrange(total_stack):
				DB[i] = OB[i]
			OB.close()
			DB.close()
	   	mpi_barrier(MPI_COMM_WORLD)
	   	if myid==main_node:
			params= []
			for i in xrange(total_stack):
				e=get_im(orgstack,i)
				phi,theta,psi,s2x,s2y = get_params_proj(e)
				params.append([phi,theta,psi,s2x,s2y])
			write_text_row(params,Tracker["constants"]["ali3d"])
		mpi_barrier(MPI_COMM_WORLD)
		"""
		#Tracker["total_stack"]             = total_stack
		Tracker["constants"]["total_stack"] = total_stack
		Tracker["shrinkage"]                = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
		#####------------------------------------------------------------------------------
		if Tracker["constants"]["mask3D"]:
			Tracker["mask3D"] = os.path.join(masterdir,"smask.hdf")
		else:Tracker["mask3D"] = None
		if Tracker["constants"]["focus3Dmask"]:  Tracker["focus3D"]=os.path.join(masterdir,"sfocus.hdf")
		else:                                    Tracker["focus3D"] = None
		if myid ==main_node:
			if Tracker["constants"]["mask3D"]:
				get_shrink_3dmask(Tracker["nxinit"],Tracker["constants"]["mask3D"]).write_image(Tracker["mask3D"])
			if Tracker["constants"]["focus3Dmask"]:
				mask_3D = get_shrink_3dmask(Tracker["nxinit"],Tracker["constants"]["focus3Dmask"])
				st = Util.infomask(mask_3D, None, True)
				if( st[0] == 0.0 ):  ERROR("sxrsort3d","incorrect focused mask, after binarize all values zero",1)
				mask_3D.write_image(Tracker["focus3D"])
				del mask_3D
		if Tracker["constants"]["PWadjustment"]:
			PW_dict={}
			nxinit_pwsp=sample_down_1D_curve(Tracker["constants"]["nxinit"],Tracker["constants"]["nnxo"],Tracker["constants"]["PWadjustment"])
			Tracker["nxinit_PW"] = os.path.join(masterdir,"spwp.txt")
			if myid ==main_node:
				write_text_file(nxinit_pwsp,Tracker["nxinit_PW"])
			PW_dict[Tracker["constants"]["nnxo"]]   =Tracker["constants"]["PWadjustment"]
			PW_dict[Tracker["constants"]["nxinit"]] =Tracker["nxinit_PW"]
			Tracker["PW_dict"] = PW_dict 
		###----------------------------------------------------------------------------------		
		####---------------------------  Extract the previous results   #####################################################
		from random import shuffle
		if myid ==main_node:
			log_main.add("Extract stable groups from previous runs")
			stable_member_list                              = get_stable_members_from_two_runs(Tracker["constants"]["previous_runs"], Tracker["constants"]["total_stack"], log_main)
			Tracker["this_unaccounted_list"], new_stable_P1 = get_leftover_from_stable(stable_member_list, Tracker["constants"]["total_stack"], Tracker["constants"]["smallest_group"])
			Tracker["this_unaccounted_list"].sort()
			Tracker["total_stack"] = len(Tracker["this_unaccounted_list"])
			log_main.add("new stable is %d"%len(new_stable_P1))
		else:
			Tracker["total_stack"]           = 0
			Tracker["this_unaccounted_list"] = 0
			stable_member_list =0
		stable_member_list               = wrap_mpi_bcast(stable_member_list, main_node)
		Tracker["total_stack"]           = bcast_number_to_all(Tracker["total_stack"], source_node = main_node)
		left_one_from_old_two_runs       = wrap_mpi_bcast(Tracker["this_unaccounted_list"], main_node)
		if myid ==main_node:  
			write_text_file(left_one_from_old_two_runs, os.path.join(masterdir,"unaccounted_from_two_previous_runs.txt"))
			print " Extracting results of two previous runs is done!"
		#################################### Estimate resolution----------------------############# 

		#### make chunkdir dictionary for computing margin of error
		chunk_list = []
		if Tracker["constants"]["chunkdir"] !="": ##inhere previous random assignment of odd and even
			if myid == main_node:
				chunk_one = read_text_file(os.path.join(Tracker["constants"]["chunkdir"],"chunk0.txt"))
				chunk_two = read_text_file(os.path.join(Tracker["constants"]["chunkdir"],"chunk1.txt"))
			else:
				chunk_one = 0
				chunk_two = 0
			chunk_one = wrap_mpi_bcast(chunk_one, main_node)
			chunk_two = wrap_mpi_bcast(chunk_two, main_node)
		else:  ## if traces are lost, then creating new random assignment of odd, even particles
			chunks = range(Tracker["constants"]["total_stack"])
			shuffle(chunks)
			chunk_one =chunks[0:Tracker["constants"]["total_stack"]//2]
			chunk_two =chunks[Tracker["constants"]["total_stack"]//2:Tracker["constants"]["total_stack"]]
			chunk_one = wrap_mpi_bcast(chunk_one, main_node)	
			chunk_two = wrap_mpi_bcast(chunk_two, main_node)	
				
		###### Fill chunk ID into headers when calling get_shrink_data_huang
		if myid ==main_node:
			print " random odd and even assignment done  !"
		mpi_barrier(MPI_COMM_WORLD)
		#------------------------------------------------------------------------------
		Tracker["chunk_dict"] = {}
		for element in chunk_one: Tracker["chunk_dict"][element] = 0
		for element in chunk_two: Tracker["chunk_dict"][element] = 1
		Tracker["P_chunk0"]   = len(chunk_one)/float(Tracker["constants"]["total_stack"])
		Tracker["P_chunk1"]   = len(chunk_two)/float(Tracker["constants"]["total_stack"])
		### create two volumes to estimate resolution
		if myid == main_node:
			write_text_file(chunk_one, os.path.join(masterdir,"chunk0.txt"))
			write_text_file(chunk_two, os.path.join(masterdir,"chunk1.txt"))
		mpi_barrier(MPI_COMM_WORLD)
		vols = []
		for index in xrange(2):
			data1,old_shifts1 = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"], os.path.join(masterdir,"chunk%d.txt"%index), Tracker["constants"]["partstack"], myid, main_node, nproc, preshift = True)
			vol1 = recons3d_4nn_ctf_MPI(myid=myid, prjlist=data1, symmetry=Tracker["constants"]["sym"], finfo=None)
			if myid ==main_node:
				vol1_file_name = os.path.join(masterdir, "vol%d.hdf"%index)
				vol1.write_image(vol1_file_name)
			
			vols.append(vol1)
			mpi_barrier(MPI_COMM_WORLD)
		if myid ==main_node:
			low_pass, falloff, currentres = get_resolution_mrk01(vols, Tracker["constants"]["radius"]*Tracker["shrinkage"], Tracker["constants"]["nxinit"], masterdir,Tracker["mask3D"])
			if low_pass   > Tracker["constants"]["low_pass_filter"]: low_pass  = Tracker["constants"]["low_pass_filter"]
		else:
			low_pass    =0.0
			falloff     =0.0
			currentres  =0.0
		currentres                    = bcast_number_to_all(currentres,source_node = main_node)
		low_pass                      = bcast_number_to_all(low_pass,source_node   = main_node)
		falloff                       = bcast_number_to_all(falloff,source_node    = main_node)
		Tracker["currentres"]         = currentres
		
		####################################################################
		
		Tracker["falloff"] = falloff
		if Tracker["constants"]["low_pass_filter"] == -1.0:
			Tracker["low_pass_filter"] = low_pass*Tracker["shrinkage"]
		else:
			Tracker["low_pass_filter"] = Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"]
		Tracker["lowpass"]             = Tracker["low_pass_filter"]
		Tracker["falloff"]             = 0.1
		Tracker["global_fsc"]          = os.path.join(masterdir,"fsc.txt")
		##################################################################
		if myid ==main_node:
			log_main.add("The command-line inputs are :")
			log_main.add("**********************************************************")
			for a in sys.argv: 
					log_main.add(a)
			log_main.add("**********************************************************")
		from filter import filt_tanl
		##################### START 3-D sorting ##########################
		if myid ==main_node:
			log_main.add("----------3-D sorting  program------- ")
			log_main.add("current resolution %6.3f for images of original size in terms of absolute frequency"%Tracker["currentres"])
			log_main.add("equivalent to %f Angstrom resolution"%(round((Tracker["constants"]["pixel_size"]/Tracker["currentres"]/Tracker["shrinkage"]),4)))
			filt_tanl(get_im(os.path.join(masterdir, "vol0.hdf")), Tracker["low_pass_filter"], 0.1).write_image(os.path.join(masterdir, "volf0.hdf"))			
			filt_tanl(get_im(os.path.join(masterdir, "vol1.hdf")), Tracker["low_pass_filter"], 0.1).write_image(os.path.join(masterdir, "volf1.hdf"))
			print " random odd and even assignment done  !"
		mpi_barrier(MPI_COMM_WORLD)
		## ---------------------------------------------------------------------------------------------########
		## Stop program and output results when the leftover from two sort3d runs is not sufficient for a new run    		########
		## ---------------------------------------------------  ---------------------------------------  ######
		Tracker["number_of_groups"] = get_number_of_groups(len(left_one_from_old_two_runs), Tracker["constants"]["number_of_images_per_group"])
		if Tracker["number_of_groups"] <=1 : # programs finishes 
			if myid == main_node:
				log_main.add("the unaccounted ones are no sufficient for a simple two-group run, output results!")
				log_main.add("this implies your two sort3d runs already achieved high reproducibale ratio. ")
				log_main.add("Or your number_of_images_per_group is too large ")
				log_main.add("the final reproducibility is  %f"%((Tracker["constants"]["total_stack"]-len(Tracker["this_unaccounted_list"]))/float(Tracker["constants"]["total_stack"])))
				for i in xrange(len(stable_member_list)): write_text_file(stable_member_list[i], os.path.join(masterdir,"P2_final_class%d.txt"%i))
				mask3d = get_im(Tracker["constants"]["mask3D"])
			else:
				mask3d = model_blank(Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
			bcast_EMData_to_all(mask3d, myid, main_node)
			for igrp in xrange(len(stable_member_list)):
				#name_of_class_file = os.path.join(masterdir, "P2_final_class%d.txt"%igrp)
				data, old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"], os.path.join(masterdir, "P2_final_class%d.txt"%igrp), Tracker["constants"]["partstack"], myid, main_node, nproc,preshift = True)
				if Tracker["constants"]["CTF"]:  
					volref, fscc = rec3D_two_chunks_MPI(data, 1.0, Tracker["constants"]["sym"], mask3d,os.path.join(masterdir,"resolution_%02d.txt"%igrp), myid, main_node, index =-1, npad=2)
				else: 
					print "Missing CTF flag!"
					from mpi import mpi_finalize
					mpi_finalize()
					exit()
				mpi_barrier(MPI_COMM_WORLD)

				#nx_of_image=volref.get_xsize()
				if Tracker["constants"]["PWadjustment"] :		Tracker["PWadjustment"] = Tracker["PW_dict"][Tracker["constants"]["nnxo"]]
				else:											Tracker["PWadjustment"] = Tracker["constants"]["PWadjustment"]	
				if myid ==main_node:
					try: 
						lowpass = search_lowpass(fscc)
						falloff = 0.1
					except:
						lowpass = 0.4
						falloff = 0.1
						log_main.add(" lowpass and falloff from fsc are %f %f"%(lowpass, falloff))
					lowpass = round(lowpass,4)
					falloff = round(min(0.1,falloff),4)
					Tracker["lowpass"] = lowpass
					Tracker["falloff"] = falloff
					refdata           = [None]*4
					refdata[0]        = volref
					refdata[1]        = Tracker
					refdata[2]        = Tracker["constants"]["myid"]
					refdata[3]        = Tracker["constants"]["nproc"]
					volref            = user_func(refdata)
					cutoff = Tracker["constants"]["pixel_size"]/lowpass
					log_main.add("%d vol low pass filer %f   %f  cut to  %f Angstrom"%(igrp,Tracker["lowpass"],Tracker["falloff"],cutoff))
					volref.write_image(os.path.join(masterdir,"volf_final%d.hdf"%igrp))
			mpi_barrier(MPI_COMM_WORLD)			
			from mpi import mpi_finalize
			mpi_finalize()
			exit()
		else: # Continue clustering on unaccounted ones that produced by two_way comparison of two previous runs
			#########################################################################################################################
			#if Tracker["constants"]["number_of_images_per_group"] ==-1: # Estimate number of images per group from delta, and scale up 
			#    or down by scale_of_number
			#	number_of_images_per_group = int(Tracker["constants"]["scale_of_number"]*len(n_angles))
			#
			#########################################################################################################################P2
			if myid ==main_node:
				print " Now continue clustering on accounted ones because they can make at least two groups!"
			P2_partitions        = []
			number_of_P2_runs    = 2  # Notice P2 start from two P1 runs
			### input list_to_be_processed
			import copy
			mpi_barrier(MPI_COMM_WORLD)
			for iter_P2_run in xrange(number_of_P2_runs): # two runs such that one can obtain reproducibility
				list_to_be_processed = left_one_from_old_two_runs[:]#Tracker["this_unaccounted_list"][:]
				Tracker["this_unaccounted_list"] = left_one_from_old_two_runs[:]
				if myid == main_node :    new_stable1 =  new_stable_P1[:]
				total_stack   = len(list_to_be_processed) # This is the input from two P1 runs
				#number_of_images_per_group = Tracker["constants"]["number_of_images_per_group"]
				P2_run_dir = os.path.join(masterdir, "P2_run%d"%iter_P2_run)
				Tracker["number_of_groups"] = get_number_of_groups(total_stack, Tracker["constants"]["number_of_images_per_group"])
				if myid == main_node:
					cmd="{} {}".format("mkdir", P2_run_dir)
					os.system(cmd)
					log_main.add("----------------P2 independent run %d--------------"%iter_P2_run)
					log_main.add("user provided number_of_images_per_group %d"%Tracker["constants"]["number_of_images_per_group"])
					print "----------------P2 independent run %d--------------"%iter_P2_run
				mpi_barrier(MPI_COMM_WORLD)
				#
				#Tracker["number_of_groups"] = get_number_of_groups(total_stack,Tracker["constants"]["number_of_images_per_group"])
				generation                  = 0
				
				if myid == main_node:
					log_main.add("number of groups is %d"%Tracker["number_of_groups"])
					log_main.add("total stack %d"%total_stack)
				while( Tracker["number_of_groups"]>=2 ):
					partition_dict      = {}
					full_dict           = {}
					workdir             = os.path.join(P2_run_dir,"generation%03d"%generation)
					Tracker["this_dir"] = workdir
					
					if myid ==main_node:
						cmd="{} {}".format("mkdir", workdir)
						os.system(cmd)
						log_main.add("---- generation         %5d"%generation)
						log_main.add("number of images per group is set as %d"%Tracker["constants"]["number_of_images_per_group"])
						log_main.add("the initial number of groups is  %d "%Tracker["number_of_groups"])
						log_main.add(" the number to be processed in this generation is %d"%len(list_to_be_processed))
						print "---- generation         %5d"%generation
						#core=read_text_row(Tracker["constants"]["ali3d"],-1)
						#write_text_row(core, os.path.join(workdir,"node%d.txt"%myid))
					mpi_barrier(MPI_COMM_WORLD)
					Tracker["this_data_list"]         = list_to_be_processed # leftover of P1 runs
					Tracker["total_stack"]            = len(list_to_be_processed)
					create_random_list(Tracker)
					
					###------ For super computer    ##############
					update_full_dict(list_to_be_processed, Tracker)
					###----
					##### ----------------Independent runs for EQ-Kmeans  ------------------------------------
					for indep_run in xrange(Tracker["constants"]["indep_runs"]):
						Tracker["this_particle_list"] = Tracker["this_indep_list"][indep_run]
						ref_vol = recons_mref(Tracker)
						if myid ==main_node:   log_main.add("independent run  %10d"%indep_run)
						mpi_barrier(MPI_COMM_WORLD)
						#this_particle_text_file =  # for get_shrink_data
						if myid ==main_node:
							write_text_file(list_to_be_processed, os.path.join(workdir, "independent_list_%03d.txt"%indep_run))
						mref_ali3d_EQ_Kmeans(ref_vol, os.path.join(workdir, "EQ_Kmeans%03d"%indep_run), os.path.join(workdir, "independent_list_%03d.txt"%indep_run), Tracker)
						partition_dict[indep_run] = Tracker["this_partition"]
						del ref_vol
					Tracker["partition_dict"]    = partition_dict
					Tracker["this_total_stack"]  = Tracker["total_stack"]
					do_two_way_comparison(Tracker)
					##############################
					
					if myid ==main_node: log_main.add("Now calculate stable volumes")
					if myid ==main_node:
						for igrp in xrange(len(Tracker["two_way_stable_member"])):
							Tracker["this_data_list"]      = Tracker["two_way_stable_member"][igrp]
							write_text_file(Tracker["this_data_list"], os.path.join(workdir,"stable_class%d.txt"%igrp))
					Tracker["this_data_list_file"] = -1
					mpi_barrier(MPI_COMM_WORLD)
					###
					number_of_ref_class = []
					ref_vol_list = []
					for igrp in xrange(len(Tracker["two_way_stable_member"])):
						data, old_shifts = get_shrink_data_huang(Tracker,Tracker["nxinit"], os.path.join(workdir, "stable_class%d.txt"%igrp), Tracker["constants"]["partstack"], myid, main_node, nproc, preshift = True)
						volref           = recons3d_4nn_ctf_MPI(myid=myid,prjlist=data,symmetry=Tracker["constants"]["sym"],finfo = None)
						ref_vol_list.append(volref)
						number_of_ref_class.append(len(Tracker["this_data_list"]))
					if myid ==main_node:
						log_main.add("group  %d  members %d "%(igrp,len(Tracker["this_data_list"])))	
						#ref_vol_list=apply_low_pass_filter(ref_vol_list,Tracker)
						for iref in xrange(len(ref_vol_list)): ref_vol_list[iref].write_image(os.path.join(workdir,"vol_stable.hdf"),iref)
						
					mpi_barrier(MPI_COMM_WORLD)
					################################
					
					Tracker["number_of_ref_class"]       = number_of_ref_class
					Tracker["this_data_list"]            = Tracker["this_accounted_list"]
					outdir = os.path.join(workdir, "Kmref")
					empty_groups,res_classes,final_list  = ali3d_mref_Kmeans_MPI(ref_vol_list, outdir, os.path.join(workdir,"Accounted.txt"), Tracker)
					Tracker["this_unaccounted_list"]     = get_complementary_elements(list_to_be_processed,final_list)
					if myid == main_node:  log_main.add("the number of particles not processed is %d"%len(Tracker["this_unaccounted_list"]))
					update_full_dict(Tracker["this_unaccounted_list"], Tracker)
					if myid == main_node: write_text_file(Tracker["this_unaccounted_list"], Tracker["this_unaccounted_text"])
					Tracker["number_of_groups"]          = len(res_classes)
					### Update data
					mpi_barrier(MPI_COMM_WORLD)
					if myid == main_node:
						number_of_ref_class=[]
						log_main.add(" Compute volumes of original size")
						for igrp in xrange(Tracker["number_of_groups"]):
							if os.path.exists( os.path.join( outdir,"Class%d.txt"%igrp ) ):
								new_stable1.append( read_text_file( os.path.join( outdir, "Class%d.txt"%igrp ) ) )
								log_main.add(" read Class file %d"%igrp)
								number_of_ref_class.append(len(new_stable1))
					else:  number_of_ref_class = 0
					number_of_ref_class = wrap_mpi_bcast(number_of_ref_class,main_node)
					
					################################
					
					mpi_barrier(MPI_COMM_WORLD)
					if myid ==main_node:
						vol_list = []
					for igrp in xrange(Tracker["number_of_groups"]):
						if myid ==main_node: log_main.add("start vol   %d"%igrp)
						data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"], os.path.join(outdir,"Class%d.txt"%igrp), Tracker["constants"]["partstack"],myid, main_node, nproc, preshift = True)
						volref          = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry = Tracker["constants"]["sym"],finfo= None)
						if myid == main_node: 
							vol_list.append(volref)
							log_main.add(" vol   %d is done"%igrp)
					Tracker["number_of_ref_class"] = number_of_ref_class
					mpi_barrier(MPI_COMM_WORLD)
					generation +=1
					#################################
					
					if myid ==main_node:
						for ivol in xrange(len(vol_list)): vol_list[ivol].write_image(os.path.join(workdir, "vol_of_Classes.hdf"),ivol)
						filt_tanl(vol_list[ivol],Tracker["constants"]["low_pass_filter"],.1).write_image(os.path.join(workdir, "volf_of_Classes.hdf"),ivol)
						log_main.add("number of unaccounted particles  %10d"%len(Tracker["this_unaccounted_list"]))
						log_main.add("number of accounted particles  %10d"%len(Tracker["this_accounted_list"]))
						del vol_list
					Tracker["this_data_list"]        = Tracker["this_unaccounted_list"]
					Tracker["total_stack"]           = len(Tracker["this_unaccounted_list"])
					Tracker["this_total_stack"]      = Tracker["total_stack"]
					#update_full_dict(complementary)
					#number_of_groups = int(float(len(Tracker["this_unaccounted_list"]))/number_of_images_per_group)
					del list_to_be_processed
					list_to_be_processed             = copy.deepcopy(Tracker["this_unaccounted_list"]) 
					Tracker["number_of_groups"]      = get_number_of_groups(len(list_to_be_processed),Tracker["constants"]["number_of_images_per_group"])
					mpi_barrier(MPI_COMM_WORLD)
					
	#############################################################################################################################
				### Reconstruct the unaccounted is only done once
			
				if (Tracker["constants"]["unaccounted"] and (len(Tracker["this_unaccounted_list"]) != 0)):
					data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"],Tracker["this_unaccounted_text"],Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
					volref = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"],finfo=None)
					volref = filt_tanl(volref, Tracker["constants"]["low_pass_filter"],.1)
					if myid ==main_node: volref.write_image(os.path.join(workdir, "volf_unaccounted.hdf"))
				
					######## Exhaustive Kmeans #############################################
				if myid ==main_node:
					if len(Tracker["this_unaccounted_list"])>=Tracker["constants"]["smallest_group"]:
						new_stable1.append(Tracker["this_unaccounted_list"])
					unaccounted                 = get_complementary_elements_total(Tracker["constants"]["total_stack"], final_list)
					Tracker["number_of_groups"] = len(new_stable1)
					log_main.add("----------------Exhaustive Kmeans------------------")
					log_main.add("number_of_groups is %d"%Tracker["number_of_groups"])
				else:    Tracker["number_of_groups"] = 0
				### prepare references for final K-means
				if myid == main_node:
					final_list =[]
					for alist in new_stable1:
						for element in alist:final_list.append(int(element))
					unaccounted = get_complementary_elements_total(Tracker["constants"]["total_stack"],final_list)
					if len(unaccounted) > Tracker["constants"]["smallest_group"]:  # treat unaccounted ones also as a group if it is not too small.
						new_stable1.append(unaccounted)
						Tracker["number_of_groups"] = len(new_stable1)
						for any in unaccounted:final_list.append(any)
					log_main.add("total number %d"%len(final_list))
				else:  final_list = 0
				Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"],source_node = main_node)
				mpi_barrier(MPI_COMM_WORLD)
				final_list = wrap_mpi_bcast(final_list, main_node)
				workdir = os.path.join(P2_run_dir,"Exhaustive_Kmeans") # new workdir 
				if myid==main_node:
					os.mkdir(workdir)
					write_text_file(final_list, os.path.join(workdir,"final_list.txt"))
				else: new_stable1 = 0
				mpi_barrier(MPI_COMM_WORLD)
				## Create reference volumes
		
				if myid == main_node:
					number_of_ref_class = []
					for igrp in xrange(Tracker["number_of_groups"]):
						class_file = os.path.join(workdir,"final_class%d.txt"%igrp)
						write_text_file(new_stable1[igrp],class_file)
						log_main.add(" group %d   number of particles %d"%(igrp,len(new_stable1[igrp])))
						number_of_ref_class.append(len(new_stable1[igrp]))
				else:  number_of_ref_class= 0
				number_of_ref_class = wrap_mpi_bcast(number_of_ref_class,main_node)
		
				mpi_barrier(MPI_COMM_WORLD)
				ref_vol_list = []
				for igrp in xrange(Tracker["number_of_groups"]):
					if myid ==main_node : print " prepare reference %d"%igrp
					#Tracker["this_data_list_file"] = os.path.join(workdir,"final_class%d.txt"%igrp)
					data,old_shifts                = get_shrink_data_huang(Tracker, Tracker["nxinit"],os.path.join(workdir,"final_class%d.txt"%igrp), Tracker["constants"]["partstack"], myid,main_node,nproc,preshift = True)
					volref                         = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo = None)
					#volref = filt_tanl(volref, Tracker["low_pass_filter"],.1)
					#if myid == main_node:
					#	volref.write_image(os.path.join(masterdir,"volf_stable.hdf"),iref)
					#volref = resample(volref,Tracker["shrinkage"])
					bcast_EMData_to_all(volref, myid, main_node)
					ref_vol_list.append(volref)
				mpi_barrier(MPI_COMM_WORLD)
				### -------variables used in Kmeans_exhaustive_run-----
				Tracker["number_of_ref_class"] = number_of_ref_class
				Tracker["this_data_list"]      = final_list
				Tracker["total_stack"]         = len(final_list)
				Tracker["this_dir"]            = workdir
				Tracker["this_data_list_file"] = os.path.join(workdir,"final_list.txt")
				KE_group                       = Kmeans_exhaustive_run(ref_vol_list,Tracker) # 
				P2_partitions.append(KE_group[:][:])
				if myid ==main_node:
					log_main.add(" the number of groups after exhaustive Kmeans is %d"%len(KE_group))
					for ike in xrange(len(KE_group)):log_main.add(" group   %d   number of objects %d"%(ike,len(KE_group[ike])))
					del new_stable1
				mpi_barrier(MPI_COMM_WORLD)
			if myid == main_node:   log_main.add("P2 runs are done, now start two-way comparision to exclude those that are not reproduced ")
			reproduced_groups = two_way_comparison_single(P2_partitions[0],P2_partitions[1],Tracker)# Here partition IDs are original indexes.
			###### ----------------Reconstruct reproduced groups------------------------#######
			######
			if myid == main_node:
				for index_of_reproduced_groups in xrange(len(reproduced_groups)):
					name_of_class_file = os.path.join(masterdir, "P2_final_class%d.txt"%index_of_reproduced_groups)
					write_text_file(reproduced_groups[index_of_reproduced_groups],name_of_class_file)
				log_main.add("-------start to reconstruct reproduced volumes individully to orignal size-----------")
				
			mpi_barrier(MPI_COMM_WORLD)
			if Tracker["constants"]["mask3D"]: mask_3d = get_shrink_3dmask(Tracker["constants"]["nnxo"],Tracker["constants"]["mask3D"])
			else:                              mask_3d = None
			
			for igrp in xrange(len(reproduced_groups)):
				data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nnxo"],os.path.join(masterdir, "P2_final_class%d.txt"%igrp),Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
				#volref = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo=None)
				if Tracker["constants"]["CTF"]: 
					volref, fscc = rec3D_two_chunks_MPI(data,1.0,Tracker["constants"]["sym"],mask_3d, \
										os.path.join(masterdir,"resolution_%02d.txt"%igrp),myid,main_node,index =-1,npad =2,finfo=None)
				else: 
					print "Missing CTF flag!"
					from mpi import mpi_finalize
					mpi_finalize()
					exit()
				mpi_barrier(MPI_COMM_WORLD)
				fscc        = read_text_file(os.path.join(masterdir, "resolution_%02d.txt"%igrp),-1)
				nx_of_image = volref.get_xsize()
				if Tracker["constants"]["PWadjustment"]:	Tracker["PWadjustment"] = Tracker["PW_dict"][nx_of_image]
				else:										Tracker["PWadjustment"] = Tracker["constants"]["PWadjustment"]	
				try:
					lowpass = search_lowpass(fscc)
					falloff = 0.1
				except:
					lowpass= 0.4
					falloff= 0.1
				print lowpass
				lowpass=round(lowpass,4)
				falloff=round(min(.1,falloff),4)
				Tracker["lowpass"]= lowpass
				Tracker["falloff"]= falloff
				if myid == main_node:
					refdata    =[None]*4
					refdata[0] = volref
					refdata[1] = Tracker
					refdata[2] = Tracker["constants"]["myid"]
					refdata[3] = Tracker["constants"]["nproc"]
					volref     = user_func(refdata)
					cutoff     = Tracker["constants"]["pixel_size"]/lowpass
					log_main.add("%d vol low pass filer %f   %f  cut to  %f Angstrom"%(igrp,Tracker["lowpass"],Tracker["falloff"],cutoff))
					volref.write_image(os.path.join(masterdir,"volf_final%d.hdf"%igrp))
		if myid==main_node:   log_main.add(" sxsort3d_P2 finishes. ")
		# Finish program
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
if __name__ == "__main__":
		main()
