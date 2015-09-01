#!/home/zhuang/EMAN2/extlib/bin/python
#
#  08/13/2015
#  New version.  
#
#
import os
import global_def
from global_def import *
from optparse import OptionParser
import sys
from numpy import array
import types
from logger import Logger, BaseLogger_Files

def print_a_line_with_timestamp(string_to_be_printed ):                 
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
 	print(line,string_to_be_printed)
	return string_to_be_printed

def convertasi(asig,K):
        p = []
        for k in xrange(K):
                l = []
                for i in xrange(len(asig)):
                        if(asig[i] == k): l.append(i)
                l = array(l, 'int32')
                l.sort()
                p.append(l)
        return p

def prepare_ptp(data_list,K):
	num_of_pt=len(data_list)
	ptp=[]
	for ipt in xrange(num_of_pt):
		ptp.append([])
	for ipt in xrange(num_of_pt):
		nc =len(data_list[ipt])
		asig=[-1]*nc
		for i in xrange(nc):
        		asig[i]=data_list[ipt][i]
		ptp[ipt] = convertasi(asig,K)
	return ptp

def print_dict(dict,theme):
        line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        print(line,theme)
        spaces = "                    "
        for key, value in sorted( dict.items() ):
                if(key != "constants"):  print("                    => ", key+spaces[len(key):],":  ",value)

def checkstep(item, keepchecking, myid, main_node):
	from utilities import bcast_number_to_all
        if(myid == main_node):
                if keepchecking:
                        if(os.path.exists(item)):
                                doit = 0
                        else:
                                doit = 1
                                keepchecking = False
                else:
                        doit = 1
        else:
                doit = 1
        doit = bcast_number_to_all(doit, source_node = main_node)
        return doit, keepchecking

def get_resolution_mrk01(vol, radi, nnxo, fscoutputdir, mask_option):
        # this function is single processor
        #  Get updated FSC curves, user can also provide a mask using radi variable
        import types
	from statistics import fsc
	from utilities import model_circle, get_im
	from filter import fit_tanh1
        if(type(radi) == int):
                if(mask_option is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
                else:                           mask = get_im(mask_option)
        else:  mask = radi
        nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
        currentres = -1.0
        ns = len(nfsc[1])
        #  This is actual resolution, as computed by 2*f/(1+f)
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.333333333333333333333333):
                        currentres = nfsc[0][i-1]
                        break
        if(currentres < 0.0):
                print("  Something wrong with the resolution, cannot continue")
                mpi_finalize()
                exit()
        """
        lowpass = 0.5
        ns = len(nfsc[1])
        #  This is resolution used to filter half-volumes
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.5 ):
                        lowpass = nfsc[0][i-1]
                        break
        """
        lowpass, falloff = fit_tanh1(nfsc, 0.01)

        return  round(lowpass,4), round(falloff,4), round(currentres,2)

def get_K_equal_assignment(input_list,K):
	import types
	from utilities import read_text_file
	if type(K) != IntType:
		return None
	if type(input_list) == StringType:
		 ll=read_text_file(input_list)
	elif type(input_list) == ListType:
		ll=input_list
	else:
		print "Error"
		return None
def partition_to_groups(alist,K):
	res =[]
	for igroup in xrange(K):
		this_group =[]
		for imeb in xrange(len(alist)):
			if alist[imeb] ==igroup:
				this_group.append(imeb)
		res.append(this_group)
	return res	
	
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
	usage = progname + " stack  outdir refvols  <mask> --focus=3Dmask --ir=inner_radius --radius=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_searching_step " +\
	" --delta=angular_step --an=angular_neighborhood --center=1 --nassign=reassignment_number --nrefine=alignment_number --maxit=max_iter --stoprnct=percentage_to_stop " + \
	" --debug --fourvar=fourier_variance --CTF --snr=1.0 --ref_a=S --sym=c1 --function=user_function --independent=indenpendent_runs_for_equalmref  --Kgroup=number_of_groups  --resolution --mode=EK --import=params.txt"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--focus",    type="string",       default=None,             help="3D mask for focused clustering ")
	parser.add_option("--ir",       type= "int",         default=1, 	       help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",   type= "int",         default="-1",	       help="outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--maxit",	type= "int",         default=5, 	       help="maximum number of iteration")
	parser.add_option("--rs",       type= "int",         default="1",	       help="step between rings in rotational correlation >0 (set to 1)" ) 
	parser.add_option("--xr",       type="string",       default="4 2 1 1 1",      help="range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",       type="string",       default="-1",	       help="range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",       type="string",       default="0.25",           help="step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",    type="string",       default="10 6 4  3   2",  help="angular step of reference projections")
	parser.add_option("--an",       type="string",       default="-1",	       help="angular neighborhood for local searches")
	parser.add_option("--center",   type="int",          default=0,	               help="0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0)")
	parser.add_option("--nassign",  type="int",          default=0, 	       help="number of reassignment iterations performed for each angular step (set to 3) ")
	parser.add_option("--nrefine",  type="int",          default=1, 	       help="number of alignment iterations performed for each angular step (set to 1) ")
	parser.add_option("--CTF",      action="store_true", default=False,            help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",        default=1.0,              help="Signal-to-Noise Ratio of the data")   
	parser.add_option("--stoprnct", type="float",        default=0.0,              help="Minimum percentage of assignment change to stop the program")   
	parser.add_option("--ref_a",    type="string",       default="S",              help="method for generating the quasi-uniformly distributed projection directions (default S) ")
	parser.add_option("--sym",      type="string",       default="c1",             help="symmetry of the structure ")
	parser.add_option("--function", type="string",       default="ref_ali3dm",     help="name of the reference preparation function")
	parser.add_option("--MPI",      action="store_true", default=False,            help="Use MPI version ")
	parser.add_option("--npad",     type="int",          default= 2,               help="padding size for 3D reconstruction")
	parser.add_option("--debug",    action="store_true", default=False,            help="debug ")
	parser.add_option("--fourvar",  action="store_true", default=False,            help="compute and use fourier variance")
	parser.add_option("--independent", type="int",       default= 3,               help="number of independent run")
	parser.add_option("--Kgroup",      type="int",       default= 5,               help="number of groups")
	parser.add_option("--resolution",  type="float",     default= .40,             help="structure is low-pass-filtered to this resolution for clustering" )
	parser.add_option("--mode",        type="string",    default="EK",             help="computing either Kmeans mref, or Equal Kmeans mref, or do both" )
	parser.add_option("--importali3d", type="string",    default="",  help="import the xform.projection parameters as the initial configuration for 3-D reconstruction" )
	
	
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1  or len(args) > 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
            	if len(args)>2:
            		mask_file = args[2]
            	else:
                	mask_file = None

		orgstack                        =args[0]
		masterdir                       =args[1]
		global_def.BATCH = True
		#---initialize MPI related variables
		from mpi import mpi_init, mpi_comm_size, MPI_COMM_WORLD, mpi_comm_rank,mpi_barrier,mpi_bcast, mpi_bcast, MPI_INT
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		nproc     = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		mpi_comm = MPI_COMM_WORLD
		main_node = 0
		# import some utilites
		from utilities import get_im,bcast_number_to_all,cmdexecute, write_text_file, read_text_file
		from applications import recons3d_n_MPI, mref_ali3d_MPI, Kmref_ali3d_MPI
		from statistics import k_means_match_clusters_asg_new,k_means_stab_bbenum 
		# Create the main log file
	        from logger import Logger,BaseLogger_Files
             	if myid ==main_node:
             		log_main=Logger(BaseLogger_Files())
                     	log_main.prefix=masterdir+"/"
                else:
                        log_main =None
		#--- fill input parameters into dictionary named after Constants
		Constants		       	         ={}
		Constants["stack"]               = args[0]
        	Constants["masterdir"]           = masterdir
		Constants["mask3D"]              = mask_file
		Constants["focus3Dmask"]         = options.focus
		Constants["indep_runs"]          = options.independent
		Constants["stoprnct"]            = options.stoprnct
		Constants["Kgroup"]              = options.Kgroup
		Constants["CTF"]                 = options.CTF
		Constants["npad"]                = options.npad
		Constants["maxit"]               = options.maxit
		Constants["ir"]                  = options.ir 
		Constants["radius"]              = options.radius 
		Constants["nassign"]             = options.nassign
		Constants["snr"]                 = options.snr
		Constants["rs"]                  = options.rs 
		Constants["xr"]                  = options.xr
		Constants["yr"]                  = options.yr
		Constants["ts"]                  = options.ts
		Constants["ref_a"]               = options.ref_a
		Constants["delta"]               = options.delta
		Constants["an"]                  = options.an
		Constants["sym"]                 = options.sym
		Constants["center"]              = options.center
		Constants["nrefine"]             = options.nrefine
		Constants["fourvar"]             = options.fourvar 
		Constants["user_func"]           = options.function
		Constants["resolution"]          = options.resolution
		Constants["debug"]               = options.debug
		Constants["sign"]                = 1
		Constants["listfile"]            = ""
		Constants["xysize"]              =-1
		Constants["zsize"]               =-1
		Constants["group"]               =-1
		Constants["verbose"]             = 0
		Constants["mode"]                = options.mode
		Constants["mpi_comm"]            = mpi_comm
		Constants["main_log_prefix"]     =args[1]
		Constants["importali3d"]         =options.importali3d	
		# -----------------------------------------------------
		#
		# Create and initialize Tracker dictionary with input options
		Tracker = 			    		{}
		Tracker["constants"]=				Constants
		Tracker["maxit"]          = Tracker["constants"]["maxit"]
        	Tracker["radius"]         = Tracker["constants"]["radius"]
        	Tracker["xr"]             = ""
        	Tracker["yr"]             = "-1"  # Do not change!
        	Tracker["ts"]             = 1
        	Tracker["an"]             = "-1"
        	Tracker["delta"]          = "2.0"
        	Tracker["zoom"]           = True
        	Tracker["nsoft"]          = 0
        	Tracker["local"]          = False
        	Tracker["PWadjustment"]   = ""
        	Tracker["upscale"]        = 0.5
        	Tracker["applyctf"]       = True  #  Should the data be premultiplied by the CTF.  Set to False for local continuous.
        	#Tracker["refvol"]         = None
        	Tracker["nxinit"]         = 64
        	Tracker["nxstep"]         = 32
        	Tracker["icurrentres"]    = -1
        	Tracker["ireachedres"]    = -1
        	Tracker["lowpass"]        = 0.4
        	Tracker["falloff"]        = 0.2
        	#Tracker["inires"]         = options.inires  # Now in A, convert to absolute before using
        	Tracker["fuse_freq"]      = 50  # Now in A, convert to absolute before using
        	Tracker["delpreviousmax"] = False
        	Tracker["anger"]          = -1.0
        	Tracker["shifter"]        = -1.0
        	Tracker["saturatecrit"]   = 0.95
        	Tracker["pixercutoff"]    = 2.0
        	Tracker["directory"]      = ""
        	Tracker["previousoutputdir"] = ""
        	#Tracker["eliminated-outliers"] = False
        	Tracker["mainiteration"]  = 0
        	Tracker["movedback"]      = False
        	#Tracker["state"]          = Tracker["constants"]["states"][0] 
		Tracker["global_resolution"] =0.0
		#--------------------------------------------------------------------
		#
		# Get the pixel size; if none, set to 1.0, and the original image size
        	if(myid == main_node):
                	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                	print(line,"INITIALIZATION OF KMREF")

                	a = get_im(orgstack)
                	nnxo = a.get_xsize()
                	if( Tracker["nxinit"] > nnxo ):
                        	ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"sxmeridien",1)
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
        	Tracker["constants"]["nnxo"]         = nnxo
        	Tracker["constants"]["pixel_size"]   = pixel_size
        	Tracker["fuse_freq"]    = fq
        	del fq, nnxo, pixel_size

        	if(Tracker["constants"]["radius"]  < 1):
                	Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
        	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
                	ERROR("Particle radius set too large!","sxEKmref_clustering.py",1,myid)

####-----------------------------------------------------------------------------------------
		# Master directory
		if myid == main_node:
			if masterdir =="":
				timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
				masterdir ="master_EKMREF"+timestring
				li =len(masterdir)
				cmd="{} {}".format("mkdir", masterdir)
				cmdexecute(cmd)
			else:
				li = 0
				keepchecking =1			
		else:
			li=0
			keepchecking =1
	    	li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
		if li>0:
			masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
			masterdir = string.join(masterdir,"")
		if myid ==main_node:
			print_dict(Tracker["constants"], "Permanent settings of MKREF")
		######### create a vstack from input stack to the local stack in masterdir
		# stack name set to default
		Tracker["constants"]["stack"] = "bdb:"+masterdir+"/rdata"
	   	if(myid == main_node):
                	if keepchecking:
                        	if(os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):  doit = False
                        	else:  doit = True
                	else:  doit = True
                	if  doit:
                        	if(orgstack[:4] == "bdb:"):     cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+Tracker["constants"]["stack"])
                        	else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, Tracker["constants"]["stack"])
                        	cmdexecute(cmd)
                        	cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["constants"]["stack"])
                        	cmdexecute(cmd)
                        	keepchecking = False
                	total_stack = EMUtil.get_image_count(Tracker["constants"]["stack"])
        	else:
                	total_stack = 0
		mpi_barrier(MPI_COMM_WORLD)
        	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
		
		###----------------------------------------------------------------------------------
		# Inititialization
		from random import shuffle
		initdir = os.path.join(masterdir,"main000")
		## fill in dictionaries
		initial_random_vol={}
		for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
			filtered_volumes =os.path.join(initdir,"volf_run_%03d.hdf"%iter_indep)
			initial_random_vol[iter_indep]= filtered_volumes
		Tracker["refvols"]=initial_random_vol
		doit, keepchecking = checkstep(initdir, keepchecking, myid, main_node)
		if doit:
			if myid ==main_node:
				cmd="mkdir "+initdir
				cmdexecute(cmd)
			mpi_barrier(MPI_COMM_WORLD)
			# Compute the resolution 
			partids =[None]*2
			for procid in xrange(2):partids[procid] =os.path.join(initdir, "chunk%01d.txt"%procid)
			inivol =[None]*2
			for procid in xrange(2):inivol[procid]  =os.path.join(initdir,"vol%01d.hdf"%procid)
					
			if myid ==main_node:
				ll=range(total_stack)
				shuffle(ll)
				l1=ll[0:total_stack//2]
				l2=ll[total_stack//2:]
				del ll
				l1.sort()
				l2.sort()
				write_text_file(l1,partids[0])
				write_text_file(l2,partids[1])
			mpi_barrier(MPI_COMM_WORLD)
			ll =[]
			for procid in xrange(2):
				l=read_text_file(partids[procid])
				ll.append(l)
			### create two volumes to estimate resolution
			vols=[]
			for procid in xrange(2):
				recons3d_n_MPI(Tracker["constants"]["stack"],ll[procid],inivol[procid],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"],\
					Tracker["constants"]["listfile"],Tracker["constants"]["group"],Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
			mpi_barrier(MPI_COMM_WORLD)
			for procid in xrange(2):
				if myid==main_node:
					vols.append(get_im(inivol[procid]))
				else:
					vols.append(None)
			if myid ==main_node:
				low_pass, falloff,currentres =get_resolution_mrk01(vols,Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],initdir,Tracker["constants"]["mask3D"])
				if low_pass >Tracker["constants"]["resolution"]:
					low_pass= Tracker["constants"]["resolution"]
			else:
				low_pass    =0.0
				falloff     =0.0
				currentres  =0.0
			bcast_number_to_all(currentres,source_node = main_node)
			bcast_number_to_all(low_pass,source_node = main_node)
			bcast_number_to_all(falloff,source_node = main_node)
			Tracker["currentres"]= currentres
			Tracker["low_pass"]  = low_pass
			Tracker["falloff"]   = falloff
			if myid ==main_node:
				log_main.add("%3d runs of %8d particles assignment to %3d groups "%(Tracker["constants"]["indep_runs"],total_stack,Tracker["constants"]["Kgroup"]))
			### generate the initial random assignment
	             	for irandom in xrange(Tracker["constants"]["indep_runs"]):
                        	ll=range(total_stack)
                        	shuffle(ll)
                        	if myid ==main_node:
					log_main.add("The initial random assignments are "+os.path.join(initdir,"random_list%d.txt"%irandom))
                                	write_text_file(ll,os.path.join(initdir,"random_list%d.txt"%irandom))
                	mpi_barrier(MPI_COMM_WORLD)
			### generate the initial volumes
			for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
                        	ll=read_text_file(os.path.join(initdir,"random_list%d.txt"%iter_indep))
                        	linit_list =[]
                        	for igrp in xrange(Tracker["constants"]["Kgroup"]):
					alist=ll[(total_stack*igrp)//Tracker["constants"]["Kgroup"]:(total_stack*(igrp+1))//Tracker["constants"]["Kgroup"]]
					alist.sort()
                                	linit_list.append(alist)
                                	linit_list_file_name=os.path.join(initdir,"list1_grp%0d_%0d.txt"%(igrp,iter_indep))
                                	write_text_file(alist,linit_list_file_name)
					volstack =os.path.join(initdir,"TMP_vol%03d.hdf"%igrp)
					recons3d_n_MPI(Tracker["constants"]["stack"],alist,volstack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["sign"],\
					Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
	                      	mpi_barrier(MPI_COMM_WORLD)
				newvols=os.path.join(initdir,"TMP_vol*.hdf")
                        	computed_refvol_name=os.path.join(initdir,"vol_run_%03d.hdf"%iter_indep)
				filtered_volumes =os.path.join(initdir,"volf_run_%03d.hdf"%iter_indep)
				filter_options ="--process=filter.lowpass.tanh:cutoff_abs="+str(Tracker["low_pass"])+":fall_off=.1"
                        	if myid==main_node:
                                	cmd ="{} {} {}".format("sxcpy.py",newvols,computed_refvol_name) # Make stacks 
                                	cmdexecute(cmd)
                                	cmd ="{} {}".format("rm",newvols)
                                	cmdexecute(cmd)
					cmd ="{}  {}  {} {}".format("e2proc3d.py",computed_refvol_name,filtered_volumes,filter_options)
					cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)

		else:# check items inside the directory
			if myid ==main_node: print "%20s exists, check inside"%initdir 
			partids =[None]*2
			doall=1 # check initial partition
			for procid in xrange(2):
				partids[procid] =os.path.join(initdir, "chunk%01d.txt"%procid) 
				doit, keepchecking = checkstep(partids[procid], keepchecking, myid, main_node)
				doall *=doit
			if doall:
				if myid ==main_node:
                                	ll=range(total_stack)
                                	shuffle(ll)
                                	l1=ll[0:total_stack//2]
                                	l2=ll[total_stack//2:]
                                	del ll
                                	l1.sort()
                                	l2.sort()
                                	write_text_file(l1,partids[0])
                                	write_text_file(l2,partids[1])
                        	mpi_barrier(MPI_COMM_WORLD)
			else:
				if myid==main_node: print_a_line_with_timestamp("The initial partition are there!")
			 # check the global resolution of the data
			inivol =[None]*2
			for procid in xrange(2):
				inivol[procid]  =os.path.join(initdir,"vol%01d.hdf"%procid)
				doit, keepchecking = checkstep(inivol[procid], keepchecking, myid, main_node)
				if doit:
					plist=read_text_file(partids[procid])
					recons3d_n_MPI(Tracker["constants"]["stack"],plist,inivol[procid],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"],\
                                        Tracker["constants"]["listfile"],Tracker["constants"]["group"],Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
				else:
					if myid ==main_node: print_a_line_with_timestamp("%20s  exists!"%inivol[procid])
			### check the resolution curve
			init_fsc_file=os.path.join(initdir, "fsc.txt")
			doit, keepchecking = checkstep(init_fsc_file, keepchecking, myid, main_node)
			if doit:
				vols =[]
				if myid ==main_node:
					for procid in xrange(2):
                                		vols.appened(get_im(os.path.join(initdir,"vol%01d.hdf"%procid)))
					low_pass, falloff,currentres =get_resolution_mrk01(vols,Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],initdir,Tracker["constants"]["mask3D"])
				else:
					vols.append(None)
					low_pass   =0.0
					falloff    =0.0
					currentres =0.0
				mpi_barrier(MPI_COMM_WORLD)
				bcast_number_to_all(currentres,source_node = main_node)
             		        bcast_number_to_all(low_pass,source_node   = main_node)
                       		bcast_number_to_all(falloff,source_node    = main_node)
                        	Tracker["currentres"]= currentres
                        	Tracker["low_pass"]  = low_pass
                        	Tracker["falloff"]   = falloff
				if myid==main_node: print_a_line_with_timestamp("recompute the global resolution curve")
                        	mpi_barrier(MPI_COMM_WORLD)
			else:
				if myid==main_node: print_a_line_with_timestamp("The initial fsc curve is there")
		#msg="The estimated resolution of the data is %5.2f"%currentres
		#if myid ==main_node:
		#	msg="The estimated resolution of the data is %5.2f"%Tracker["currentres"]
		#	log_main.add(msg)
		mpi_barrier(MPI_COMM_WORLD)
		# Start independent runs, first fill in a dictionary for this run
		if myid ==main_node:
			msg="The current mode is %5s"%Tracker["constants"]["mode"]
			log_main.add(msg)
			msg="The number of independent runs is %5d"%Tracker["constants"]["indep_runs"]
			log_main.add(msg)
		if Tracker["constants"]["mode"][0:1]=="E": 
			EKMREF_iteration  ={}
			for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
				output_for_this_run ={}
				output_for_this_run["output_dir"]= os.path.join(masterdir,"EKMREF%03d"%iter_indep)
				output_for_this_run["partition"] = os.path.join(output_for_this_run["output_dir"],"list2.txt")
				output_for_this_run["refvols"]   = os.path.join(initdir,"volf_run_%03d.hdf"%iter_indep)
				output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
				#output_for_this_run["refvols"]=os.path.join(Tracker["second_main"],"vol_strable_member%03d.hf"%iter_indep			
				EKMREF_iteration[iter_indep]=output_for_this_run
			Tracker["EMREF"]=EKMREF_iteration
			mpi_barrier(MPI_COMM_WORLD)
		# Carry out independent Equal Kmeans runs
			for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
				outdir = Tracker["EMREF"][iter_indep]["output_dir"]
				doit, keepchecking = checkstep(outdir, keepchecking, myid, main_node)
				log_Emref=Logger(BaseLogger_Files())
				log_Emref.prefix=outdir+"/"
				if myid ==main_node:
					msg="The %5d th run equal Kmeans, outdir is %20s"%(iter_indep,outdir)
					log_main.add(msg)
					if Tracker["constants"]["importali3d"]!="":
						cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--import="+Tracker["constants"]["importali3d"])
                                                cmdexecute(cmd)

				mpi_barrier(MPI_COMM_WORLD)
				if doit:
					mref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["EMREF"][iter_indep]["refvols"],outdir,Tracker["constants"]["mask3D"] ,\
                        		Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                        		Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                        		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                        		Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Emref)
					if myid==main_node:
						cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=group", "--export="+Tracker["EMREF"][iter_indep]["partition"])
                                        	cmdexecute(cmd)
						if Tracker["constants"]["nrefine"]:
							cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--export="+Tracker["EMREF"][iter_indep]["ali3d"])
							cmdexecute(cmd)
				else:
					doit, keepchecking = checkstep(Tracker["EMREF"][iter_indep]["partition"],keepchecking,myid,main_node)
					if doit:
						if myid==main_node:
							cmd="rm -rf "+outdir
							cmdexecute(cmd)
			        		mref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["refvols"][iter_indep],outdir,Tracker["constants"]["mask3D"] ,\
                                		Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                                		Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                                		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                                		Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Emref)
						if myid==main_node:
							cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=group", "--export="+Tracker["EMREF"][iter_indep]["partition"])
							cmdexecute(cmd)
							if Tracker["constants"]["nrefine"]:
								cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--export="+Tracker["EMREF"][iter_indep]["ali3d"])
								cmdexecute(cmd)
					else:
						if myid==main_node:print_a_line_with_timestamp("MREF%03d has been done"%iter_indep)
				mpi_barrier(MPI_COMM_WORLD)
		else: 
			KMREF_iteration  ={}	
			for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
                        	output_for_this_run ={}
                        	output_for_this_run["output_dir"]=os.path.join(masterdir,"KMREF%03d"%iter_indep)
                        	output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
				output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
                        	KMREF_iteration[iter_indep]=output_for_this_run
                	Tracker["KMREF"]=KMREF_iteration
                	mpi_barrier(MPI_COMM_WORLD)
                	# Carry out independent runs of Kmeans
                	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
                        	outdir = Tracker["KMREF"][iter_indep]["output_dir"]
                        	doit, keepchecking = checkstep(outdir, keepchecking, myid, main_node)
				log_Kmref=Logger(BaseLogger_Files())
                                log_Kmref.prefix=Tracker["KMREF"][iter_indep]["output_dir"]+"/"
				if Tracker["constants"]["importali3d"]!="":
					cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--import="+Tracker["constants"]["importali3d"])
					cmdexecute(cmd)
                        	if doit:
                                	empty_group = Kmref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["refvols"][iter_indep],outdir,Tracker["constants"]["mask3D"] ,\
                                	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                                	Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                                	Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                                	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Kmref )
					if myid ==main_node:
						if empty_group ==1:
							log_main.add("K-means encounters an empty group, need re-create reference files!")
                                	if myid==main_node:
                                        	cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"], "--params=group", "--export="+Tracker["KMREF"][iter_indep]["partition"])
                                        	cmdexecute(cmd)
						if racker["constants"]["nrefine"]:
							cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--export="+Tracker["KMREF"][iter_indep]["ali3d"])
							cmdexecute(cmd)
                        	else:
                                	doit, keepchecking = checkstep(Tracker["KMREF"][iter_indep]["partition"],keepchecking,myid,main_node)
                                	if doit:
                                        	if myid==main_node:
                                                	cmd="rm -rf "+outdir
                                                	cmdexecute(cmd)
						mpi_barrier(MPI_COMM_WORLD)
                                        	empty_group = Kmref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["refvols"][iter_indep],outdir,Tracker["constants"]["mask3D"] ,\
                                        	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                                        	Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                                        	Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                                        	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Kmref)
                                        	if myid==main_node:
                                                	cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"], "--params=group", "--export="+Tracker["KMREF"][iter_indep]["partition"])
                                                	cmdexecute(cmd)
							if Tracker["constants"]["nrefine"]:
                                        			cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=xform.projection", "--export="+Tracker["KMREF"][iter_indep]["ali3d"])
                                        			cmdexecute(cmd)
                                		mpi_barrier(MPI_COMM_WORLD)
                                	else:
                                        	if myid==main_node:print_a_line_with_timestamp("KMREF%03d has been done"%iter_indep)
		### Now do some analysis	
		mpi_barrier(MPI_COMM_WORLD)
		if Tracker["constants"]["indep_runs"] >1:
			if myid ==main_node:
				msg="Two-comparisions are done on indepndent runs of equal Kmeans"
				log_main.add(msg)
                	total_partition=[]
                	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
				if Tracker["constants"]["mode"][0:1]=="E":
					partition_list=read_text_file(Tracker["EMREF"][iter_indep]["partition"])
					if myid ==main_node: 
						log_main.add("Equal Kmeans  with %d         %d"%(total_stack, Tracker["constants"]["Kgroup"])+"\n")
				else: 
					partition_list=read_text_file(Tracker["KMREF"][iter_indep]["partition"])
					if myid ==main_node: 
						log_main.add("Kmeans  with %d         %d"%(total_stack, Tracker["constants"]["Kgroup"])+"\n")
                       		total_partition.append(partition_list)
                	mpi_barrier(MPI_COMM_WORLD)
                	### Conduct two ways comparision Keep all nodes busy
                	ptp=prepare_ptp(total_partition,Tracker["constants"]["Kgroup"])
                	# three ways comparision:
                	#TCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(ptp, T=0, J=50, max_branching=40, stmult=0.1, branchfunc=2)#  T is the minimum group size
                	#tc =0.
                	#for a in STB_PART:
                	#        tc +=len(a)
                	#rate=tc/float(nima)*100.
                	#if  myid==0: print "%5.2f"%rate,"%","three way comparision"
			total_pop=0
			two_ways_stable_member_list={}
			avg_two_ways =0.0
			avg_two_ways_square =0.
			scores ={}
                	for iptp in xrange(len(ptp)):
                        	for jptp in xrange(len(ptp)):
                                	newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
                                	tt =0.
                                	for m in xrange(len(list_stable)):
                                        	a=list_stable[m]
                                        	tt +=len(a)
                                	rate=tt/total_stack*100.0
					scores[(iptp,jptp)]=rate
					scores[(jptp,iptp)]=rate
					avg_two_ways 	    +=rate
					avg_two_ways_square +=rate**2
                                	if myid ==main_node and iptp !=jptp:
                                        	aline=print_a_line_with_timestamp("two-way comparison  %3d %3d %5.3f  "%(iptp,jptp,rate))
						log_main.add(aline)
                                	new_list=[]
                                	for a in list_stable:
                                        	a.astype(int)
                                        	new_list.append(a)
                                	two_ways_stable_member_list[(iptp,jptp)]=new_list
					total_pop +=1
			summed_scores =[]
			two_way_dict ={}
			for ipp in xrange(len(ptp)):
				avg_scores =0.0
				for jpp in xrange(len(ptp)):
					if ipp !=jpp:
						avg_scores +=scores[(ipp,jpp)]
				avg_rate =avg_scores/(len(ptp)-1)
				summed_scores.append(avg_rate)
				two_way_dict[avg_rate] =ipp
			rate1=max(summed_scores)
			run1 =two_way_dict[rate1]
			summed_scores.remove(rate1)
			rate2 =max(summed_scores)
			run2 = two_way_dict[rate2]
			Tracker["two_way_stable_member"]=two_ways_stable_member_list[(run1,run2)]
			if myid==main_node:
				log_main.add(" selected indepedent runs are %5d and  %5d"%(run1,run2))
				log_main.add(" their averaged rates are %5.2f  and %5.2f "%(rate1,rate2))		
			from math import sqrt
			avg_two_ways = avg_two_ways/total_pop
			two_ways_std=sqrt(avg_two_ways_square/total_pop-avg_two_ways**2)
			if myid ==main_node: 
				aline=print_a_line_with_timestamp("average of two-way comparison  %5.3f"%avg_two_ways)
				log_main.add(aline)
			if myid ==main_node: 
				aline=print_a_line_with_timestamp("std of two-way comparison %5.3f"%two_ways_std)
				log_main.add(aline)
			two_way_of_two_way_stable_member_list = []
			two_way_of_two_way_pop =0
			if myid ==main_node:
				msg="two-way of two-way comparision"
				log_main.add(msg)
			"""
                	for itwo_ways in xrange(len(two_ways_stable_member_list)-1):
                		for jtwo_ways in xrange(itwo_ways+1,len(two_ways_stable_member_list)):
					tt =0.
                        		newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(two_ways_stable_member_list[itwo_ways], two_ways_stable_member_list[jtwo_ways])
                                	tt= 0.
					for m in xrange(len(list_stable)):
						a=list_stable[m]
						tt +=len(a)
					new_list=[]
					for a in list_stable:
						a.astype(int)
						new_list.append(a)
					two_way_of_two_way_stable_member_list.append(new_list)
					rate=float(nb_tot_objs)/total_stack*100.0
                                	if myid ==main_node:
                                		aline =print_a_line_with_timestamp("two-way of two-way comparision ratio is %5.3f"%rate+"%")
						log_main.add(aline)
					two_way_of_two_way_pop +=1
			Tracker["two_way_pop"]			=total_pop
			Tracker["two_way_stable"]		=two_ways_stable_member_list
			Tracker["two_way_of_two_way_pop"]	= two_way_of_two_way_pop
			Tracker["two_way_of_two_way_stable"]	=two_way_of_two_way_stable_member_list
			"""
			## Now carry out Kmref using stable members selected by  pairwise comparsion 
			if Tracker["constants"]["mode"]=="K" or Tracker["constants"]["mode"]=="E":
				if myid ==main_node:
					msg="This is indepedent runs of %s-means"%Tracker["constants"]["mode"]
					log_main.add(msg)
					msg="Program stops here."
					log_main.add(msg) 
				mpi_barrier(MPI_COMM_WORLD)
				from mpi import mpi_finalize
                		mpi_finalize()
				exit()
			else:
				## first fill in dictionary for this run
				## Create the second main directory 
				if myid ==main_node: 
					log_main.add("Now do  K-means using stable members given by equal K-means")
				second_main=os.path.join(masterdir,"main001")
				Tracker["second_main"]= second_main
                        	EKKMREF_iteration ={}
				#EK_total_runs = Tracker["two_way_pop"]
                        	for inew in xrange(1):# new independent run
                                	output_for_this_run ={}
                                	output_for_this_run["output_dir"]=os.path.join(masterdir,"EKKMREF%03d"%inew)
                                	output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
					output_for_this_run["refvols"]=os.path.join(Tracker["second_main"],"vol_stable_member_%03d.hdf"%inew)
                                	EKKMREF_iteration[inew]=output_for_this_run
                        	Tracker["EKKMREF"]=EKKMREF_iteration
				doit, keepchecking = checkstep(Tracker["second_main"],keepchecking,myid,main_node)
				if doit:
					if myid==main_node:
						cmd="mkdir "+second_main
						cmdexecute(cmd)
					for inew in xrange(1):
						if myid ==main_node:
							msg="Create %dth two-way reference  model"%inew
							log_main.add(msg)
		        			for igrp in xrange(Tracker["constants"]["Kgroup"]):
                                			pid_list=Tracker["two_way_stable_member"][igrp]
                                			vol_stack=os.path.join(Tracker["second_main"],"TMP_init%03d.hdf"%igrp)
							tlist=[]
							for b in pid_list:
								tlist.append(int(b))
							mpi_barrier(MPI_COMM_WORLD)
							if myid ==main_node:
								msg="%d th model %d th group   particle number %5d"%(inew,igrp,len(tlist))
								log_main.add(msg)
                                			recons3d_n_MPI(Tracker["constants"]["stack"],tlist,vol_stack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],\
							Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
							Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                 				newvols=os.path.join(Tracker["second_main"],"TMP_init*.hdf")
                        			if myid==main_node:
                                			cmd ="{} {} {}".format("sxcpy.py",newvols,Tracker["EKKMREF"][inew]["refvols"]) # Make stacks 
                                			cmdexecute(cmd)
                                			cmd ="{} {}".format("rm",newvols)
                                			cmdexecute(cmd)
                        			mpi_barrier(MPI_COMM_WORLD)
					for inew in xrange(1):
						# Start KMREF
						log_Kmref=Logger(BaseLogger_Files())
						log_Kmref.prefix=Tracker["EKKMREF"][inew]["output_dir"]+"/" 
						if myid ==main_node:
							msg="Kmref_ali3d_MPI with two-way equal K-means stable members as starting point "
							log_main.add(msg)
							msg = "%5d "%inew+Tracker["EKKMREF"][inew]["output_dir"]
							log_main.add(msg)
                        			Kmref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["EKKMREF"][inew]["refvols"],Tracker["EKKMREF"][inew]["output_dir"],Tracker["constants"]["mask3D"],\
				        	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                                        	Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                                        	Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                                        	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Kmref)
                        			if myid==main_node:
                                        		cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=group", "--export="+Tracker["EKKMREF"][inew]["partition"])
                                        		cmdexecute(cmd)
                        			mpi_barrier(MPI_COMM_WORLD)
					"""
					## Analysis on KMREF results
					total_partition=[]
                			for iter_indep  in xrange(Tracker["two_way_pop"]):
                        			partition_list=read_text_file(Tracker["EKKMREF"][iter_indep]["partition"])
                        			total_partition.append(partition_list)
                			### Conduct two ways comparision Keep all nodes busy
                			two_ways_stable_member_list=[]
                			ptp=prepare_ptp(total_partition,Tracker["constants"]["Kgroup"])
                			# three ways comparision:
                			#TCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(ptp, T=0, J=50, max_branching=40, stmult=0.1, branchfunc=2)#  T is the minimum group size
                			#tc =0.
                			#for a in STB_PART:
                        		#	tc +=len(a)
                			#rate=tc/float(nima)*100.
                			#if myid ==main_node: print "%5.2f"%rate,"%","three-way comparision"
                			for iptp in xrange(len(ptp)-1):
                        			for jptp in xrange(iptp+1,len(ptp)):
                                			newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
                                			tt =0.
                                			for m in xrange(len(list_stable)):
                                        			a=list_stable[m]
                                        			tt +=len(a)
                                			rate=tt/total_stack*100.0
                                			if myid ==main_node:
                                        			aline =print_a_line_with_timestamp("two-way comparision  rate %5.2f"%rate)
								log_main.add(aline)
                                			two_ways_stable_member_list.append(list_stable)
                        		# Compare stability of two ways comparision!
                        		for itwo_ways in xrange(len(two_ways_stable_member_list)-1):
                                		for jtwo_ways in xrange(itwo_ways+1,len(two_ways_stable_member_list)):
                                        		newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(two_ways_stable_member_list[itwo_ways], two_ways_stable_member_list[jtwo_ways])
                                        		rate=float(nb_tot_objs)/total_stack*100.0
                                        		if myid ==main_node:
                                                		aline =print_a_line_with_timestamp("the two-way of two-way matching ratio is %5.3f"%rate+"%")
								log_main.add(aline)
				else:
					for inew in xrange(Tracker["two_way_pop"]):
						doit, keepchecking = checkstep(Tracker["EKKMREF"][inew]["refvols"],keepchecking,myid,main_node)
						if doit:
							for igrp in xrange(Tracker["constants"]["Kgroup"]):
                                                		pid_list=Tracker["two_way_stable"][inew][igrp]
                                                		vol_stack=os.path.join(Tracker["second_main"],"TMP_init%03d.hdf"%igrp)
                                                		tlist=[]
                                                		for b in pid_list:
                                                        		tlist.append(int(b))
                                                		mpi_barrier(MPI_COMM_WORLD)
                                                		recons3d_n_MPI(Tracker["constants"]["stack"],tlist,vol_stack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],\
                                                		Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
                                                		Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                                        		newvols=os.path.join(Tracker["second_main"],"TMP_init*.hdf")
                                        		if myid==main_node:
                                                		cmd ="{} {} {}".format("sxcpy.py",newvols,Tracker["EKKMREF"][inew]["refvols"]) # Make stacks 
                                                		cmdexecute(cmd)
                                                		cmd ="{} {}".format("rm",newvols)
                                                		cmdexecute(cmd)
						else:
							if myid==main_node: 
								aline=print_a_line_with_timestamp("%d  refvol is there "%inew )
								log_main.add(aline)
                                        		mpi_barrier(MPI_COMM_WORLD)
					for inew in xrange(Tracker["two_way_pop"]):
                                        	doit, keepchecking = checkstep(Tracker["EKKMREF"][inew]["partition"],keepchecking,myid,main_node)
						if doit:
							if myid ==main_node:
								doit, keepchecking = checkstep(Tracker["EKKMREF"][inew]["output_dir"],keepchecking,myid,main_node)
								if doit :cmd ="{} {}".format("rm -rf",Tracker["EKKMREF"][inew]["output_dir"])
							mpi_barrier(MPI_COMM_WORLD)
							log_Kmref=Logger(BaseLogger_Files())
							log_Kmref.prefix=Tracker["EKKMREF"][inew]["output_dir"]+"/"
							Kmref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["EKKMREF"][inew]["refvols"],Tracker["EKKMREF"][inew]["output_dir"],Tracker["constants"]["mask3D"],\
                                        		Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                                        		Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                                        		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                                        		Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"],mpi_comm,log_Kmref)
                                        		if myid==main_node:
                                                		cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=group", "--export="+Tracker["EKKMREF"][inew]["partition"])
                                                		cmdexecute(cmd)
						else:
							if myid==main_node: 
								aline=print_a_line_with_timestamp(Tracker["EKKMREF"][inew]["partition"]+"is there")
								log_main.add(aline)
                                        		mpi_barrier(MPI_COMM_WORLD)
					### Analysis is required at any case:
					total_partition=[]
                                	for iter_indep  in xrange(Tracker["two_way_pop"]):
                                        	partition_list=read_text_file(Tracker["EKKMREF"][iter_indep]["partition"])
                                        	total_partition.append(partition_list)
                                	### Conduct two ways comparision Keep all nodes busy
                                	two_ways_stable_member_list=[]
                                	ptp=prepare_ptp(total_partition,Tracker["constants"]["Kgroup"])
                                	# three ways comparision:
                                	#TCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(ptp, T=0, J=50, max_branching=40, stmult=0.1, branchfunc=2)#  T is the minimum group size
                                	#tc =0.
                                	#for a in STB_PART:
                                	#       tc +=len(a)
                                	#rate=tc/float(nima)*100.
                                	#if myid ==main_node: print "%5.2f"%rate,"%","three-way comparision"
					avg =0.0
					avg2=0.0
					total_pop =0
                                	for iptp in xrange(len(ptp)-1):
                                        	for jptp in xrange(iptp+1,len(ptp)):
                                                	newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
                                                	tt =0.
                                                	for m in xrange(len(list_stable)):
                                                        	a=list_stable[m]
                                                        	tt +=len(a)
                                                	rate=tt/total_stack*100.0
							avg  +=rate
							avg2 +=rate*rate 
                                                	if myid ==main_node:
                                                        	aline =print_a_line_with_timestamp("two-way comparision  rate%5.3f"%rate)
								log_main.add(aline)
                                                	two_ways_stable_member_list.append(list_stable)
							total_pop +=1
                                        # Compare stability of two ways comparision!
					from math import sqrt
					avg =avg/total_pop
					std=sqrt(avg2/total_pop-avg*avg)
					if myid ==maind_node:
						aline=print_a_line_with_timestamp("avg and std rate of two-way is %5.1f %5.1f"%(avg,std))
						log_main(aline)
					avg  =0.0
					avg2 =0.0
					total_pop=0
                                	for itwo_ways in xrange(len(two_ways_stable_member_list)-1):
                                        	for jtwo_ways in xrange(itwo_ways+1,len(two_ways_stable_member_list)):
                                                	newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(two_ways_stable_member_list[itwo_ways], two_ways_stable_member_list[jtwo_ways])
                                                	rate=float(nb_tot_objs)/total_stack*100.0
							total_pop +=1
                                                	if myid ==main_node:
                                                        	aline=print_a_line_with_timestamp("the two-way of two-way matching ratio is %5.3f"%rate+"%")
								log_main.add(aline)
					avg =avg/total_pop
                                	std=sqrt(avg2/total_pop-avg*avg)
					if myid ==maind_node:
                                        	aline=print_a_line_with_timestamp("avg and std rate of two-way of two-way is %5.1f %5.1f"%(avg,std))
						log_main.add(aline)
			"""
			mpi_barrier(MPI_COMM_WORLD)
			from mpi import mpi_finalize
			mpi_finalize()
		else: # single Equal Kmeans and single K-means
			second_main=os.path.join(masterdir,"main001")
                        Tracker["second_main"]= second_main
			doit, keepchecking = checkstep(Tracker["second_main"],keepchecking,myid,main_node)
                        if doit:
                        	if myid==main_node:
                                	cmd="mkdir "+second_main
                                        cmdexecute(cmd)
					log_main.add("second main is created!")
			else:
				if myid==main_node:
					log_main.add("second main is existing")
			mpi_barrier(MPI_COMM_WORLD)
			## define EKKMREF varibles
			EKKMREF_iteration ={}
                        for inew in xrange(Tracker["constants"]["indep_runs"]):# new independent run
                               	output_for_this_run ={}
                                output_for_this_run["output_dir"]=os.path.join(masterdir,"EKKMREF%03d"%inew)
                               	output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
                                output_for_this_run["refvols"]=os.path.join(Tracker["second_main"],"vol_stable_member_%03d.hdf"%inew)
                                EKKMREF_iteration[inew]=output_for_this_run
                                Tracker["EKKMREF"]=EKKMREF_iteration
			for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
				res_of_EQ=read_text_file(Tracker["EMREF"][iter_indep]["partition"])
				a_list_of_groups=partition_to_groups(res_of_EQ,Tracker["constants"]["Kgroup"])
				for igrp in xrange(Tracker["constants"]["Kgroup"]):
                               		plist = a_list_of_groups[igrp]
                                        vol_stack=os.path.join(Tracker["second_main"],"TMP_init%03d.hdf"%igrp)
                                        if myid ==main_node:
                                        	msg="%d th model %d th group   particle number %5d"%(iter_indep,igrp,len(plist))
                                               	log_main.add(msg)
                                        recons3d_n_MPI(Tracker["constants"]["stack"],plist,vol_stack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],\
                                        Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
                                        Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                                        newvols=os.path.join(Tracker["second_main"],"TMP_init*.hdf")
                              	if myid==main_node:
                                    	cmd ="{} {} {}".format("sxcpy.py",newvols,Tracker["EKKMREF"][inew]["refvols"]) # Make stacks 
                                     	cmdexecute(cmd)
                                       	cmd ="{} {}".format("rm",newvols)
                                       	cmdexecute(cmd)
				mpi_barrier(MPI_COMM_WORLD)
			for inew in xrange(1):
                        # Start KMREF
                        	log_Kmref=Logger(BaseLogger_Files())
                                log_Kmref.prefix=Tracker["EKKMREF"][inew]["output_dir"]+"/"
                                if myid ==main_node:
                                      	msg="Kmref_ali3d_MPI with two-way equal K-means stable members as starting point "
                                      	log_main.add(msg)
                                      	msg = "%5d "%inew+Tracker["EKKMREF"][inew]["output_dir"]
                                      	log_main.add(msg)
                              	Kmref_ali3d_MPI(Tracker["constants"]["stack"],Tracker["EKKMREF"][inew]["refvols"],Tracker["EKKMREF"][inew]["output_dir"],Tracker["constants"]["mask3D"],\
                              	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                              	Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                               	Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                               	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Kmref)
                               	if myid==main_node:
                                       	cmd = "{} {} {} {}".format("sxheader.py",Tracker["constants"]["stack"],"--params=group", "--export="+Tracker["EKKMREF"][inew]["partition"])
                                       	cmdexecute(cmd)
                                mpi_barrier(MPI_COMM_WORLD)
			if myid ==main_node:
				log_main.add(" Kmref_ali3d_MPI is done!")
			mpi_barrier(MPI_COMM_WORLD)
                        from mpi import mpi_finalize
                        mpi_finalize()
			exit()

if __name__ == "__main__":
	main()
