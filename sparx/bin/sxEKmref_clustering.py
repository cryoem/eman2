#!/usr/bin/env python
#
#  08/13/2015
#  New version.  
#
#
import os
import global_def
from   global_def import *
from   optparse  import OptionParser
import sys
from   numpy     import array
import types
from   logger    import Logger, BaseLogger_Files

def print_upper_triangular_matrix(data_table_dict,N_indep,log_main):
        msg =""
        for i in xrange(N_indep):
                msg +="%7d"%i
        log_main.add(msg)
        for i in xrange(N_indep):
                msg ="%5d "%i
                for j in xrange(N_indep):
                        if i<j:
                                msg +="%5.2f "%data_table_dict[(i,j)]
                        else:
                                msg +="      "
                log_main.add(msg)
			
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
	elif type(input_list) == types.ListType:
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
		this_group.sort()
		res.append(this_group)
	return res

def partition_independent_runs(run_list,K):
	indep_runs_groups={}
	for indep in xrange(len(run_list)):
		this_run = run_list[indep]
		groups = partition_to_groups(this_run,K)
		indep_runs_groups[indep]=groups 
	return indep_runs_groups

def get_outliers(total_number,plist):
        tlist={}
        for i in xrange(total_number):
                tlist[i]=i
        for a in plist:
                del tlist[a]
        out =[]
        for a in tlist:
                out.append(a)
        return out

def merge_groups(stable_members_list):
	alist=[]
	for i in xrange(len(stable_members_list)):
		for j in xrange(len(stable_members_list[i])):
			alist.append(stable_members_list[i][j])
	return alist

def save_alist(Tracker,name_of_the_text_file,alist):
	from utilities import write_text_file
	import os
	log       =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	dir_to_save_list =Tracker["this_dir"]
	if myid==main_node:
		file_name=os.path.join(dir_to_save_list,name_of_the_text_file)
		write_text_file(alist, file_name)

def reconstruct_grouped_vols(Tracker,name_of_grouped_vols,grouped_plist):
	from mpi import  MPI_COMM_WORLD,mpi_barrier
	from applications import recons3d_n_MPI
	from utilities import cmdexecute
	import os
	myid		 = Tracker["constants"]["myid"]
	main_node	 = Tracker["constants"]["main_node"]
	main_dir         = Tracker["this_dir"]
	number_of_groups = Tracker["number_of_groups"]
	data_stack       = Tracker["this_data_stack"]
	log              = Tracker["constants"]["log_main"]
	####
	if myid ==main_node:
		log.add("Reconstruct a group of volumes")
	for igrp in xrange(number_of_groups):
		#pid_list=Tracker["two_way_stable_member"][igrp]
		pid_list =grouped_plist[igrp]
                vol_stack=os.path.join(main_dir,"TMP_init%03d.hdf"%igrp)
                tlist=[]
                for b in pid_list:
                	tlist.append(int(b))
                recons3d_n_MPI(data_stack,tlist,vol_stack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],\
                Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"], \
		Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
                Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                newvols=os.path.join(main_dir,"TMP_init*.hdf")
        if myid==main_node:
       		#cmd ="{} {} {}".format("sxcpy.py",newvols,Tracker["EKKMREF"][inew]["refvols"]) # Make stacks
		cmd ="{} {} {}".format("sxcpy.py",newvols,name_of_grouped_vols)  
        	cmdexecute(cmd)
                cmd ="{} {}".format("rm",newvols)
                cmdexecute(cmd)
        mpi_barrier(MPI_COMM_WORLD)

def N_independent_reconstructions(Tracker):
	from applications import recons3d_n_MPI
	from utilities import write_text_file, read_text_file, cmdexecute
	import os
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from random import shuffle
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	initdir   =Tracker["this_dir"]
	data_stack=Tracker["this_data_stack"]
	total_stack =Tracker["this_total_stack"] 
	log_main =Tracker["constants"]["log_main"]
	if myid ==main_node:
		log_main.add("-----------Independent reconstructions---------------")
        for irandom in xrange(Tracker["constants"]["indep_runs"]):
        	ll=range(total_stack)
                shuffle(ll)
                if myid ==main_node:
                	log_main.add("The initial random assignments are "+os.path.join(initdir,"random_list%d.txt"%irandom))
                        write_text_file(ll,os.path.join(initdir,"random_list%d.txt"%irandom))
                        log_main.add("preset ali3d_parameters ")
	mpi_barrier(MPI_COMM_WORLD)
	if myid ==main_node:
         	if Tracker["importali3d"]!="":
         		cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=xform.projection", "--import="+ \
				Tracker["importali3d"])
                	cmdexecute(cmd)
			log_main.add("ali3d_parameters is preset to "+Tracker["importali3d"])
		else:
			log_main.add(" Parameters for this run are not altered !")
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
       		ll=read_text_file(os.path.join(initdir,"random_list%d.txt"%iter_indep))
                linit_list =[]
                for igrp in xrange(Tracker["number_of_groups"]):
                       	alist=ll[(total_stack*igrp)//Tracker["number_of_groups"]:(total_stack*(igrp+1))//Tracker["number_of_groups"]]
                       	alist.sort()
                       	linit_list.append(alist)
                	linit_list_file_name=os.path.join(initdir,"list1_grp%0d_%0d.txt"%(igrp,iter_indep))
			if myid ==main_node:
                		write_text_file(alist,linit_list_file_name)
                       	volstack =os.path.join(initdir,"TMP_vol%03d.hdf"%igrp)
                       	recons3d_n_MPI(data_stack,alist,volstack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["sign"],\
                       	Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
			Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
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

def N_independent_Kmref(mpi_comm,Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import Kmref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid            = Tracker["constants"]["myid"]
	main_node       = Tracker["constants"]["main_node"]
	log_main        = Tracker["constants"]["log_main"]
	data_stack      = Tracker["this_data_stack"]
	if myid ==main_node:
		log_main.add("-----------%5d independent runs of Kmref -----------------"%Tracker["constants"]["indep_runs"])
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	log_Kmref       =Logger(BaseLogger_Files())
                log_Kmref.prefix=Tracker["KMREF"][iter_indep]["output_dir"]+"/"
                empty_group =Kmref_ali3d_MPI(data_stack,Tracker["KMREF"][iter_indep]["refvols"],Tracker["KMREF"][iter_indep]["output_dir"], \
		Tracker["constants"]["mask3D"],Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],\
		Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],\
		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
		Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"],mpi_comm,log_Kmref)
                if myid==main_node:
                	cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=group", \
					"--export="+Tracker["KMREF"][iter_indep]["partition"])
                	cmdexecute(cmd)
                mpi_barrier(MPI_COMM_WORLD)
	
def do_independent_EKmref(Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import Kmref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid            = Tracker["constants"]["myid"]
        main_node       = Tracker["constants"]["main_node"]
        log_main        = Tracker["constants"]["log_main"]
        data_stack      = Tracker["this_data_stack"]
	iruns           = Tracker["this_iruns"]
	for iter_indep in xrange(iruns):
		log_Kmref       =Logger(BaseLogger_Files())
		log_Kmref.prefix=Tracker["EKMREF"][iter_indep]["output_dir"]+"/"
		empty_group =Kmref_ali3d_MPI(data_stack,Tracker["EKMREF"][iter_indep]["refvols"],Tracker["EKMREF"][iter_indep]["output_dir"], \
                Tracker["constants"]["mask3D"],Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],\
                Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],\
		Tracker["constants"]["an"],Tracker["constants"]["center"],\
                Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],\
                Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
                Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"],Tracker["constants"]["mpi_comm"],log_Kmref)
                if myid==main_node:
                        cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=group",\
					 "--export="+Tracker["EKMREF"][iter_indep]["partition"])
                        cmdexecute(cmd)
                mpi_barrier(MPI_COMM_WORLD)

def N_independent_mref(mpi_comm,Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import mref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	log_main  = Tracker["constants"]["log_main"]
	data_stack= Tracker["this_data_stack"]
	if myid ==main_node:
        	log_main.add("-----------%5d independent runs of Equal Kmeans -----------------"%Tracker["constants"]["indep_runs"])
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
		if myid==main_node:
			log_main.add(" %d th run  starts !  "%iter_indep)
        	outdir = Tracker["EMREF"][iter_indep]["output_dir"]
                #doit, keepchecking = checkstep(outdir, keepchecking, myid, main_node)
                log_Emref=Logger(BaseLogger_Files())
                log_Emref.prefix=outdir+"/"
                if Tracker["importali3d"]!="" and myid ==main_node:
                	cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=xform.projection", \
				 "--import="+Tracker["importali3d"])
                        cmdexecute(cmd)
               	mpi_barrier(MPI_COMM_WORLD)
               	mref_ali3d_MPI(data_stack,Tracker["EMREF"][iter_indep]["refvols"],outdir,Tracker["constants"]["mask3D"] ,\
               	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],\
		Tracker["constants"]["radius"],Tracker["constants"]["rs"],Tracker["constants"]["xr"],\
		Tracker["constants"]["yr"],Tracker["constants"]["ts"],\
               	Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],\
               	Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
               	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
		Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Emref)
               	if myid==main_node:
             		cmd = "{} {} {} {}".format("sxheader.py",Tracker["this_data_stack"],"--params=group",\
			 "--export="+Tracker["EMREF"][iter_indep]["partition"])
                        cmdexecute(cmd)
                        if Tracker["constants"]["nrefine"]:
                        	cmd = "{} {} {} {}".format("sxheader.py",Tracker["this_data_stack"],"--params=xform.projection",\
				 "--export="+Tracker["EMREF"][iter_indep]["ali3d"])
                        	cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)

def prepare_EMREF_dict(Tracker):
	main_dir =Tracker["this_dir"]
	EMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
        	output_for_this_run["output_dir"]= os.path.join(main_dir,"EMREF_independent_%03d"%iter_indep)
        	output_for_this_run["partition"] = os.path.join(output_for_this_run["output_dir"],"list2.txt")
        	output_for_this_run["refvols"]   = os.path.join(main_dir,"volf_run_%03d.hdf"%iter_indep)
        	output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
        	EMREF_iteration[iter_indep]=output_for_this_run
        Tracker["EMREF"]=EMREF_iteration

def prepare_KMREF_dict(Tracker):
	main_dir =Tracker["this_dir"]
	KMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
                output_for_this_run["output_dir"]= os.path.join(main_dir,"KMREF_independent_%03d"%iter_indep)
                output_for_this_run["partition"] = os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["refvols"]   = os.path.join(main_dir,"volf_run_%03d.hdf"%iter_indep)
                output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
                KMREF_iteration[iter_indep]=output_for_this_run
        Tracker["KMREF"]=KMREF_iteration

def set_refvols_table(initdir,vol_pattern,Tracker):
	# vol_pattern="volf_run_%03d.hdf"; volf_stable_%3d.hdf
	import os
	initial_random_vol={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	filtered_volumes =os.path.join(initdir,vol_pattern%iter_indep)
                initial_random_vol[iter_indep]=filtered_volumes
        Tracker["refvols"]=initial_random_vol

def prepare_EKMREF_dict(Tracker):
	main_dir = Tracker["this_dir"]
	KMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
                output_for_this_run["output_dir"]=os.path.join(main_dir,"EKMREF%03d"%iter_indep)
                output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
                KMREF_iteration[iter_indep]=output_for_this_run
        Tracker["EKMREF"]=KMREF_iteration

def partition_ali3d_params_of_outliers(Tracker):
	from utilities import read_text_file,write_text_file
	from mpi import mpi_barrier, MPI_COMM_WORLD
	myid 	   = Tracker["constants"]["myid"]
	main_node  = Tracker["constants"]["main_node"]
	outlier_list=read_text_file(Tracker["this_uncounted"])
        counted_list=read_text_file(Tracker["this_counted"])
	ali3d_params_of_outliers = []
	ali3d_params_of_counted  = []
	if Tracker["constants"]["importali3d"] !="":
     		ali3d_params_of_outliers = []
        	ali3d_params_of_counted  = []
		ali3d_params=read_text_file(Tracker["this_ali3d"])
		for i in xrange(len(outlier_list)):
			ali3d_params_of_outliers.append(ali3d_params[outlier_list[i]])
		if myid ==main_node:
			write_text_file(ali3d_params_of_outliers,Tracker["ali3d_of_outliers"])
		for i in xrange(len(counted_list)):
			ali3d_params_of_counted.append(ali3d_params[counted_list[i]])
		if myid ==main_node:
			write_text_file(ali3d_params_of_counted,Tracker["ali3d_of_counted"])
	Tracker["number_of_uncounted"]=len(outlier_list)
	Tracker["average_members_in_a_group"]= len(counted_list)/float(Tracker["number_of_groups"])
	mpi_barrier(MPI_COMM_WORLD)

def do_two_way_comparison(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_file
	from statistics import k_means_match_clusters_asg_new
	######
	myid              =Tracker["constants"]["myid"]
	main_node         =Tracker["constants"]["main_node"]
	log_main          =Tracker["constants"]["log_main"]
	total_stack       =Tracker["this_total_stack"]
	main_dir          =Tracker["this_dir"]
	number_of_groups  =Tracker["number_of_groups"]
	######
	if myid ==main_node:
		msg="-------Two_way comparisons analysis of %3d independent runs of equal Kmeans-------"%Tracker["constants"]["indep_runs"]
		log_main.add(msg)
        total_partition=[]
	if Tracker["constants"]["indep_runs"]<2:
		if myid ==main_node:
			log_main.add(" Error! One single run cannot make two-way comparison")
		from mip import mpi_finalize
		mpi_finalize()
                exit()
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
		if Tracker["constants"]["mode"][0:1]=="E":
			partition_list=read_text_file(Tracker["EMREF"][iter_indep]["partition"])
			if myid ==main_node: 
				log_main.add("Equal Kmeans  with %d         %d"%(total_stack,number_of_groups)+"\n")
		else: 
			partition_list=read_text_file(Tracker["KMREF"][iter_indep]["partition"])
			if myid ==main_node: 
				log_main.add("Kmeans  with %d         %d"%(total_stack,number_of_groups)+"\n")
       		total_partition.append(partition_list)
       	### Two-way comparision is carried out on all nodes 
       	ptp=prepare_ptp(total_partition,number_of_groups)
	indep_runs_to_groups =partition_independent_runs(total_partition,number_of_groups)
	total_pop=0
	two_ways_stable_member_list={}
	avg_two_ways               =0.0
	avg_two_ways_square        =0.
	scores                     ={}
       	for iptp in xrange(len(ptp)):
        	for jptp in xrange(len(ptp)):
               		newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp],ptp[jptp])
                       	tt =0.
			if myid ==main_node and iptp<jptp:
				aline=print_a_line_with_timestamp("Two-way comparison between independent run %3d and %3d"%(iptp,jptp))
				log_main.add(aline)
                        for m in xrange(len(list_stable)):
                        	a=list_stable[m]
                               	tt +=len(a)
				if myid==main_node and iptp<jptp:
					aline=print_a_line_with_timestamp("Group %d  number of stable members %10d  "%(m,len(a)))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The comparison is between %3d th group of %3d th run and %3d th group of %3d th run" \
					 									%(newindeces[m][0],iptp,newindeces[m][1],jptp))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The  %3d th group of %3d th run contains %6d members" \
						        %(iptp,newindeces[m][0],len(indep_runs_to_groups[iptp][newindeces[m][0]])))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The  %3d th group of %3d th run contains %6d members"%(jptp,newindeces[m][1],\
								len(indep_runs_to_groups[jptp][newindeces[m][1]]))) 
			if myid==main_node and iptp<jptp:
				uncounted=total_stack-tt
				ratio_uncounted = 100.-tt/total_stack*100.
				ratio_counted   = tt/total_stack*100
				aline           = print_a_line_with_timestamp("Counted data is %6d, %5.2f "%(int(tt),ratio_counted))
				log_main.add(aline)
				aline=print_a_line_with_timestamp("Unounted data is %6d, %5.2f"%(int(uncounted),ratio_uncounted))
				log_main.add(aline)
                       	rate=tt/total_stack*100.0
			scores[(iptp,jptp)] =rate
			if iptp<jptp:
				avg_two_ways 	    +=rate
				avg_two_ways_square +=rate**2
				total_pop +=1
                       	if myid ==main_node and iptp<jptp:
                               	aline=print_a_line_with_timestamp("The two-way comparison stable member total ratio %3d %3d %5.3f  "%(iptp,jptp,rate))
				log_main.add(aline)
                       	new_list=[]
                       	for a in list_stable:
                       		a.tolist()
                                new_list.append(a)
                        two_ways_stable_member_list[(iptp,jptp)]=new_list
	if myid ==main_node:
		log_main.add("two_way comparison is done!")
	#### Score each independent run by pairwise summation
	summed_scores =[]
	two_way_dict  ={}
	for ipp in xrange(len(ptp)):
		avg_scores =0.0
		for jpp in xrange(len(ptp)):
			if ipp!=jpp:
				avg_scores +=scores[(ipp,jpp)]
		avg_rate =avg_scores/(len(ptp)-1)
		summed_scores.append(avg_rate)
		two_way_dict[avg_rate] =ipp
	#### Select two independent runs that have the first two highest scores
	if Tracker["constants"]["indep_runs"]>2:
		rate1=max(summed_scores)
		run1 =two_way_dict[rate1]
		summed_scores.remove(rate1)
		rate2 =max(summed_scores)
		run2 = two_way_dict[rate2]
	else:
		run1 =0
		run2 =1
		rate1  =max(summed_scores)
		rate2  =rate1
	Tracker["two_way_stable_member"]      =two_ways_stable_member_list[(run1,run2)]
	Tracker["pop_size_of_stable_members"] =1
	if myid ==main_node:
                log_main.add("Get outliers of the selected comparison")
	####  Save counted ones and uncounted ones
	counted_list = merge_groups(two_ways_stable_member_list[(run1,run2)])
	outliers     = get_outliers(total_stack,counted_list)
	if myid ==main_node:
                log_main.add("Save outliers")
	save_alist(Tracker,"Uncounted.txt",outliers)
	mpi_barrier(MPI_COMM_WORLD)
	save_alist(Tracker,"Counted.txt",counted_list)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_uncounted_dir"]=main_dir
	Tracker["this_uncounted"]    =os.path.join(main_dir,"Uncounted.txt")
	Tracker["this_counted"]      =os.path.join(main_dir,"Counted.txt")
	Tracker["ali3d_of_outliers"] =os.path.join(main_dir,"ali3d_params_of_outliers.txt")
	Tracker["ali3d_of_counted"]  =os.path.join(main_dir, "ali3d_params_of_counted.txt")
	if myid==main_node:
		log_main.add(" Selected indepedent runs are %5d and  %5d"%(run1,run2))
		log_main.add(" Their pair-wise averaged rates are %5.2f  and %5.2f "%(rate1,rate2))		
	from math import sqrt
	avg_two_ways = avg_two_ways/total_pop
	two_ways_std = sqrt(avg_two_ways_square/total_pop-avg_two_ways**2)
	net_rate     = avg_two_ways-1./number_of_groups*100.
	if myid ==main_node: 
		msg="average of two-way comparison  %5.3f"%avg_two_ways
		log_main.add(msg)
		msg="net rate of two-way comparison  %5.3f"%net_rate
		log_main.add(msg)
		msg="std of two-way comparison %5.3f"%two_ways_std
		log_main.add(msg)
		msg ="Score table of two_way comparison when Kgroup =  %5d"%number_of_groups
		log_main.add(msg)
		print_upper_triangular_matrix(scores,Tracker["constants"]["indep_runs"],log_main)
	del two_ways_stable_member_list
	Tracker["score_of_this_comparison"]=(avg_two_ways,two_ways_std,net_rate)
	mpi_barrier(MPI_COMM_WORLD)
	
def do_EKmref(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import cmdexecute
	## Use stable members of only one two-way comparison
	main_dir  =Tracker["this_dir"]
	pop_size  =Tracker["pop_size_of_stable_members"]# Currently is 1, can be changed
	log_main  =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	import os
	###
	EKMREF_iteration ={}
	if myid ==main_node:
		log_main.add("-----------------EKmref---------------------------")
	for inew in xrange(1):# new independent run
        	output_for_this_run ={}
                output_for_this_run["output_dir"]=os.path.join(main_dir,"EKMREF%03d"%inew)
                output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["refvols"]=os.path.join(main_dir,"vol_stable_member_%03d.hdf"%inew)
                EKMREF_iteration[inew]=output_for_this_run
        Tracker["EKMREF"]=EKMREF_iteration
        for inew in xrange(1):
        	if myid ==main_node:
                	msg="Create %dth two-way reference  model"%inew
                	log_main.add(msg)
                name_of_grouped_vols=Tracker["EKMREF"][inew]["refvols"]
                #Tracker["this_data_stack"] =Tracker["constants"]["stack"]
                grouped_plist=Tracker["two_way_stable_member"]
                reconstruct_grouped_vols(Tracker,name_of_grouped_vols,grouped_plist)
	########
	Tracker["this_iruns"]=1
	Tracker["data_stack_of_counted"]="bdb:"+os.path.join(Tracker["this_dir"],"counted")
	partition_ali3d_params_of_outliers(Tracker)
	if myid==main_node:
		cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["this_data_stack"],"--makevstack",\
			Tracker["data_stack_of_counted"],"--list="+Tracker["this_counted"])
       		cmdexecute(cmd)
		#cmd = "{} {} {} {}  {}".format("sxheader.py", Tracker["data_stack_of_counted"],"--params=xform.projection",\
		#		 "--import="+Tracker["ali3d_of_counted"], "--consecutive ")
		#cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_data_stack"]=Tracker["data_stack_of_counted"]
	do_independent_EKmref(Tracker)

def Kgroup_guess(Tracker,data_stack):
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	log_main  = Tracker["constants"]["log_main"]
	if myid==main_node:
		msg="Program Kgroup_guess"
		log_main.add(msg)
		msg="Estimate how many groups the dataset can be partitioned into"
		log_main.add(msg)
	number_of_groups=2
	Tracker["this_data_stack"]=data_stack
	maindir                   =Tracker["this_dir"]
	terminate                 =0
	while terminate !=1:
		do_EKTC_for_a_given_number_of_groups(Tracker,number_of_groups)
		(score_this, score_std_this, score_net)=Tracker["score_of_this_comparison"]
		number_of_groups +=1
		if number_of_groups >10:
			terminate  =1
def do_N_groups(Tracker,data_stack):
	myid      = Tracker["constants"]["myid"]
        main_node = Tracker["constants"]["main_node"]
        log_main  = Tracker["constants"]["log_main"]
        if myid==main_node:
                msg="Program do_N_groups"
                log_main.add(msg)
                msg="Estimate how many groups the dataset can be partitioned into"
                log_main.add(msg)
        Tracker["this_data_stack"]=data_stack
        maindir                   =Tracker["this_dir"]
	number_of_groups          =Tracker["number_of_groups"]
	for Ngroup in xrange(2,number_of_groups+1):
		Tracker["this_data_stack"]=data_stack
		Tracker["this_dir"]       =maindir
		Tracker["number_of_groups"]=Ngroup
		do_EKTC_for_a_given_number_of_groups(Tracker,Ngroup)
                (score_this, score_std_this, score_net)=Tracker["score_of_this_comparison"]
		Tracker["this_dir"]       =os.path.join(maindir,"GRP%03d"%Tracker["number_of_groups"])
		do_EKmref(Tracker)

def do_EKTC_for_a_given_number_of_groups(Tracker,given_number_of_groups):
	import os
	from utilities import cmdexecute
	Tracker["number_of_groups"]= given_number_of_groups
	from mpi import mpi_barrier,MPI_COMM_WORLD
	log_main  = Tracker["constants"]["log_main"]
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	maindir   = Tracker["this_dir"]
	workdir   = os.path.join(maindir,"GRP%03d"%Tracker["number_of_groups"])
	if myid ==main_node:
		msg ="Number of group is %5d"%given_number_of_groups
		log_main.add(msg)
 		cmd="{} {}".format("mkdir",workdir)
        	cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_dir"]=workdir
	N_independent_reconstructions(Tracker)
	prepare_EMREF_dict(Tracker)
	N_independent_mref(Tracker["constants"]["mpi_comm"],Tracker)
	do_two_way_comparison(Tracker)
	Tracker["this_dir"]=maindir # always reset maindir

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
	parser.add_option("--mode",        type="string",    default="EK_only",             help="computing either Kmeans mref, or Equal Kmeans mref, or do both" )
	parser.add_option("--importali3d", type="string",    default="",               help="import the xform.projection parameters as the initial configuration for 3-D reconstruction" )
	parser.add_option("--do_uncounted",  action="store_true",default=False,        help="continue clustering on uncounted images" )
	parser.add_option("--Kgroup_guess",  action="store_true",default=False,        help="Guess the possible number of groups existing in one dataset" )
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
		from utilities import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file  #get_shrink_data
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
		Constants		         ={}
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
		Constants["myid"]	         =myid
		Constants["main_node"]           =main_node
		Constants["log_main"]            =log_main
		Constants["do_uncounted"]        =options.do_uncounted
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
                	print(line,"INITIALIZATION OF Equal-Kmeans & Kmeans Clustering")

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
			print_dict(Tracker["constants"], "Permanent settings of Equal-Kmeans & Kmeans Clustering")
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
		# Initial data analysis 
		from random import shuffle
		# Compute the resolution 
		partids =[None]*2
		for procid in xrange(2):partids[procid] =os.path.join(masterdir, "chunk%01d.txt"%procid)
		inivol =[None]*2
		for procid in xrange(2):inivol[procid]  =os.path.join(masterdir,"vol%01d.hdf"%procid)
					
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
			recons3d_n_MPI(Tracker["constants"]["stack"],ll[procid],inivol[procid],Tracker["constants"]["CTF"],\
			Tracker["constants"]["snr"],Tracker["constants"]["sign"],Tracker["constants"]["npad"],\
			Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
			Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
		mpi_barrier(MPI_COMM_WORLD)
		for procid in xrange(2):
			if myid==main_node:
				vols.append(get_im(inivol[procid]))
			else:
				vols.append(None)
		if myid ==main_node:
			low_pass, falloff,currentres =get_resolution_mrk01(vols,Tracker["constants"]["radius"],\
			Tracker["constants"]["nnxo"],masterdir,Tracker["constants"]["mask3D"])
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
			log_main.add("The command-line inputs are as following:")
			log_main.add("**********************************************************")
		for a in sys.argv:
			if myid ==main_node:
				log_main.add(a)
			else:
				continue
		if myid ==main_node:
			log_main.add("**********************************************************")
		### START EK Clustering
		if myid ==main_node:
                        log_main.add("----------EK Clustering------- ")
		if Tracker["constants"]["mode"]=="Kgroup_guess":
			if myid ==main_node:
                        	log_main.add("------Current mode is Kgroup_guess------")
			Tracker["this_dir"] =os.path.join(masterdir,"Kgroup_guess")
			if myid ==main_node:
                        	cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        	cmdexecute(cmd)
                	Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                	Tracker["this_total_stack"] =total_stack
                	Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                	Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                	Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
	                ### Create data stack
         		if myid ==main_node:
                        	log_main.add("----create data stack ")
                        	cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                        		Tracker["this_data_stack"])
                        	cmdexecute(cmd)
                        	cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                        	cmdexecute(cmd)
                        	if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	"--import="+Tracker["importali3d"])
                        	cmdexecute(cmd)
                	mpi_barrier(MPI_COMM_WORLD)
			Kgroup_guess(Tracker,Tracker["this_data_stack"])
			if myid ==main_node:
				msg ="Kgroup guessing is done!"
				log_main.add(msg)
		elif Tracker["constants"]["mode"]=="EK_only":
			if myid==main_node:
				msg ="--------Current mode is EK_only--------"
                                log_main.add(msg)
			Tracker["this_dir"] =os.path.join(masterdir,"EK_only")
	             	if myid ==main_node:
                                cmd="{} {}".format("mkdir",Tracker["this_dir"])
                                cmdexecute(cmd)
			Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                        Tracker["this_total_stack"] =total_stack
                        Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                        Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                        Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			### Create data stack
			if myid ==main_node:
                        	log_main.add("----Create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                       		if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                	 "--import="+Tracker["importali3d"])
                        		cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
			N_independent_reconstructions(Tracker)
			prepare_EMREF_dict(Tracker)
                	N_independent_mref(mpi_comm,Tracker)
                	do_two_way_comparison(Tracker)
                	mpi_barrier(MPI_COMM_WORLD)
			do_EKmref(Tracker)
			if myid ==main_node:
                                msg ="Equal-kmeans and K-means is done!"
                                log_main.add(msg)
			if Tracker["constants"]["do_uncounted"]:
				if myid ==main_node:
					msg ="----Continue EK clustering on those uncounted members-----"
					log_main.add(msg)
                                Tracker["this_dir"]   =os.path.join(Tracker["this_uncounted_dir"],"Uncounted")
				Tracker["last_data_stack"] =Tracker["constants"]["stack"]
				Tracker["this_data_stack"] ="bdb:"+os.path.join(Tracker["this_dir"],"data")
				Tracker["number_of_groups"]=Tracker["constants"]["Kgroup"]
				if myid ==main_node:
                        		cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        		cmdexecute(cmd)
					cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["last_data_stack"],"--makevstack",\
                                        Tracker["this_data_stack"],"--list="+Tracker["this_uncounted"])
                                        cmdexecute(cmd)
					#cmd = "{} {} {} {}  {}".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                        #      "--import="+Tracker["ali3d_of_outliers"], "--consecutive")
                                        #cmdexecute(cmd)
				mpi_barrier(MPI_COMM_WORLD)
				Tracker["this_total_stack"] =Tracker["number_of_uncounted"]
				do_N_groups(Tracker,Tracker["this_data_stack"])
                        	mpi_barrier(MPI_COMM_WORLD)
		elif Tracker["constants"]["mode"]=="auto_search":
			generation =0
			if myid ==main_node:
				log_main.add("Current mode is auto_search")
                        	log_main.add("clustering generation %3d"%generation)
                	Tracker["this_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                	if myid ==main_node:
                        	cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        	cmdexecute(cmd)
                	Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                	Tracker["this_total_stack"] =total_stack
                	Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                	Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                	Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			### Create data stack
			if myid ==main_node:
                        	log_main.add("----Create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                        Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                        	if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	"--import="+Tracker["importali3d"])
                        		cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
			if myid==main_node:
				log_main.add("Create referece volumes for %3d independent runs"%Tracker["constants"]["indep_runs"])
			N_independent_reconstructions(Tracker)
			if myid ==main_node:
				log_main.add("----------Equal-Kmeans---------")
			prepare_EMREF_dict(Tracker)
			N_independent_mref(mpi_comm,Tracker)
			do_two_way_comparison(Tracker)
			mpi_barrier(MPI_COMM_WORLD)
			if myid ==main_node:
				log_main.add("--------EKmref------- ")
			do_EKmref(Tracker)
			if myid ==main_node:
                        	log_main.add("----------EKmref---------")
			number_of_groups=min(Tracker["number_of_groups"],int(Tracker["number_of_uncounted"]/float(Tracker["average_members_in_a_group"])))
			Tracker["last_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
			Tracker["last_data_stack"]="bdb:"+os.path.join(Tracker["last_dir"],"data")
			if myid ==main_node:
				log_main.add("Kgroup for the next generation is %d"%number_of_groups)
				log_main.add("the number of outliers is  %5d"%Tracker["number_of_uncounted"]) 
				log_main.add("average_members_in_a_group  %5d"%int(Tracker["average_members_in_a_group"]))
			mpi_barrier(MPI_COMM_WORLD)
			while number_of_groups >=2 and Tracker["constants"]["do_uncounted"]:
				generation            +=1
				if myid ==main_node:
                        		log_main.add("clustering generation %03d"%generation)
                		Tracker["this_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                		if myid ==main_node:
                        		cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        		cmdexecute(cmd)
                		Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                		Tracker["this_total_stack"] =Tracker["number_of_uncounted"]
                		Tracker["number_of_groups"] =number_of_groups
				Tracker["this_ali3d"]       =Tracker["ali3d_of_outliers"]
				### Create data stack
				if myid==main_node:
					log_main.add("----Create data stack uncounted in the last generation")
                			cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["last_data_stack"],"--makevstack",\
                        		Tracker["this_data_stack"],"--list="+Tracker["this_uncounted"])
                			cmdexecute(cmd)
					cmd = "{} {}".format("sxheader.py  --consecutive", Tracker["this_data_stack"])
                        		cmdexecute(cmd)
                			#cmd = "{} {} {} {}".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	#"--import="+Tracker["this_ali3d"])
					#cmdexecute(cmd)
				mpi_barrier(MPI_COMM_WORLD)
				if myid ==main_node:
                        		log_main.add("-----Create referece volumes for %3d independent runs"%Tracker["constants"]["indep_runs"])
                		N_independent_reconstructions(Tracker)
				if myid ==main_node:
                        		log_main.add("---------Equal-Kmeans----------")
                		prepare_EMREF_dict(Tracker)
                		N_independent_mref(mpi_comm,Tracker)
				do_two_way_comparison(Tracker)
                		mpi_barrier(MPI_COMM_WORLD)
                		if myid ==main_node:
                        		log_main.add("---------EKmref---------")
                		do_EKmref(Tracker)
                		if myid ==main_node:
                        		log_main.add("------- end of EKmref-------")
                		number_of_groups=int(Tracker["number_of_uncounted"]/Tracker["average_members_in_a_group"])
				Tracker["last_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                		Tracker["last_data_stack"]="bdb:"+os.path.join(Tracker["last_dir"],"data")
		# Finish program
               	mpi_barrier(MPI_COMM_WORLD)
               	from mpi import mpi_finalize
               	mpi_finalize()
		exit()
if __name__ == "__main__":
	main()
