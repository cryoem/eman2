#!/usr/bin/env python

from mpi import MPI_SUM, mpi_reduce, mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, MPI_INT, MPI_CHAR
from utilities import bcast_number_to_all
import global_def
from global_def import *
from socket import gethostname
import string
import numpy as np
import os
#from multi_shc import *

#from debug_mpi import mpi_barrier 


def print_with_time_info(msg):
	from   time import localtime, strftime
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>" + msg
	print line

def cmdexecute(cmd):
	from   time import localtime, strftime
	import subprocess
	outcome = subprocess.call(cmd, shell=True)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(outcome == 1):
		print(  line,"ERROR!!   Command failed:  ", cmd)
		from sys import exit
		exit()
	else:  print(line,"Executed successfully: ",cmd)

def string_found_in_file(myregex, filename):
	import re
	pattern = re.compile(myregex)
	for line in open(filename):
		if re.findall(pattern, line) <> []:
			return True
	return False

def main():
	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	from logger import Logger, BaseLogger_Files
	import sys
	import os
	import time
	import socket
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc, do_volume





	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory]  [initial_volume]  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1  --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",                help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="angular step of reference projections (default 2)")
	parser.add_option("--an",       type="string", default= "-1",               help="angular neighborhood for local searches (phi and theta)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit1",    type="float",  default= 400,               help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--maxit2",    type="float",  default= 30,                help="maximum number of iterations performed for the finishing up part (set to 30) ")
	parser.add_option("--L2threshold", type="float",  default= 0.05,            help="Stopping criterion of GA given as a maximum relative dispersion of L2 norms (set to 0.05) ")
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--n_shc_runs",    type="int",    default= 24,                  help="number of quasi-independent runs (shc) (default=3)")
	parser.add_option("--n_rv_runs",       type= "int",   default= 30,                  help="number of r_viper runs")
	parser.add_option("--n_v_runs",       type= "int",   default= 3,                  help="number of viper runs for each r_viper cycle")
	parser.add_option("--doga",     type="float",  default= 0.3,                help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga (default=0.3)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	#parser.add_option("--MPI",      action="store_true", default=True,          help="whether to use MPI version - this is always set to True")

	#options introduced for the do_volume function
	parser.add_option("--fl",      type="float",  default=0.12,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	#parser.add_option("--fl",      type="float",  default=0.2,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",      type="float",  default=0.1,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	#parser.add_option("--aa",      type="float",  default=0.2,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--pwreference",      type="string",  default="",    help="shows up only in two locations in multi_shc.py")
	parser.add_option("--mask3D",      type="string",  default=None,    help="shows up only in two locations in multi_shc.py")
			


	(options, args) = parser.parse_args(sys.argv[1:])
	if len(args) < 1 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1

	masterdir = ""
	bdb_stack_location = ""
	if len(args) == 2:
		masterdir = args[1]
		if masterdir[-1] != "/":
			masterdir += "/"


	number_of_rrr_viper_runs = options.n_rv_runs 
	no_of_viper_runs_analyzed_together = options.n_v_runs 
	no_of_shc_runs_analyzed_together = options.n_shc_runs 

	
	main_node = 0
	mpi_init(0, [])


	print "Hostname:", socket.gethostname(), "proc_id:", os.getpid()

	mpi_comm = MPI_COMM_WORLD

	#mpi_barrier(mpi_comm)
	#from mpi import mpi_finalize
	#mpi_finalize()
	#print "mpi finalize"
	#from sys import exit
	#exit()


	#print "Hostname2:", socket.gethostname()
	quit_program = 0

	log = Logger(BaseLogger_Files())

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if(myid == main_node):
		print "Location:",  os.getcwd()

	#if(myid == main_node):
	#	print "masterdir:", masterdir
	#mpi_finalize()
	#sys.exit()
	
	#Create folder for all results or check if there is one created already
	if(myid == main_node):
		#cmd = "{}".format("Rmycounter ccc")
		#cmdexecute(cmd)
	
		if( masterdir == ""):
			#timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			timestring = strftime("%Y_%m_%d__%H_%M_%S/", localtime())
			masterdir = "master"+timestring
			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
		else:
			if not os.path.exists(masterdir):
				ERROR('Output directory does not exist, please change the name and restart the program', "sxrviper", 1)
				quit_program = 1

		if mpi_size % no_of_shc_runs_analyzed_together != 0:
			#print "LLL:",  mpi_size, no_of_viper_runs_analyzed_together
			#ERROR('Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.', 'sxviper', 1)
			ERROR('Number of processes needs to be a multiple of total number of runs. '
			'Total quasi-independent runs by default are 3, you can change it by specifying '
			'--nruns option. Also, to improve communication time it is recommended that '
			'the number of processes divided by the number of quasi-independent runs is a power '
			'of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has).', 'sxviper', 1)
			quit_program = 1

		#### old way	
		##make copy of database in the current master directorycopy of database in the current master directory
		#if(not os.path.exists(os.path.join(masterdir,"EMAN2DB/"))):
			#cmd = "{} {}".format("cp -rp EMAN2DB", masterdir, "EMAN2DB/")
			#cmdexecute(cmd)
	
		# new way
		# upload on cvs
		bdb_stack_location = args[0].split(":")[0] + ":" + masterdir + args[0].split(":")[1]	
		org_stack_location = args[0]

		if(not os.path.exists(os.path.join(masterdir,"EMAN2DB/"))):
			cmd = "{} {} {}".format("e2bdb.py", org_stack_location,"--makevstack=" + masterdir)
			cmdexecute(cmd)
			cmd = "{} {}".format("sxheader.py  --consecutive  --params=original_image_index", bdb_stack_location)
			cmdexecute(cmd)



		all_projs = EMData.read_images(bdb_stack_location)
		print "XXXXXXXXXXXXXXXXX"
		print "Number of projections:", len(all_projs)
		print "XXXXXXXXXXXXXXXXX"
		subset = range(len(all_projs))
		if mpi_size > len(all_projs):
			ERROR('Number of processes supplied by --np needs to be less than or equal to %d (total number of images) ' % len(all_projs), 'sxviper', 1)
			quit_program = 1


	else:
		all_projs = None
		subset = None
	
	quit_program = mpi_bcast(quit_program, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
	if quit_program > 0:
		mpi_finalize()
		sys.exit()


	#mpi_barrier(MPI_COMM_WORLD)
	#mpi_subcomm = wrap_mpi_split(mpi_comm, number_of_independent_runs)

	count_already_finished_runs = 0
	###  NEED to update regarding number_of_rrr_viper_runs < iteration_start
	###  NEED to update regarding number_of_rrr_viper_runs < iteration_start
	###  NEED to update regarding number_of_rrr_viper_runs < iteration_start
	iteration_start = get_latest_directory_increment_value(masterdir, "main")
	#mpi_finalize()
	#sys.exit()
	#Iteration over the independent runs of viper function
	for rviper_iter in range(iteration_start, number_of_rrr_viper_runs + 1):
		quit_program = 0
		#print "quit0:", rviper_iter, quit_program
		if (rviper_iter > 1):
			if(myid == main_node):
				all_projs = EMData.read_images(bdb_stack_location + "_%03d"%rviper_iter)
				print "XXXXXXXXXXXXXXXXX"
				print "Number of projections (in loop):", len(all_projs)
				print "XXXXXXXXXXXXXXXXX"
				subset = range(len(all_projs))
				if mpi_size > len(all_projs):
					ERROR('Number of processes supplied by --np needs to be less than or equal to %d (current number of images) ' % len(all_projs), 'sxviper', 1)
					quit_program = 1
			else:
				all_projs = None
				subset = None
			quit_program = mpi_bcast(quit_program, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
			if quit_program > 0:
				print "qqqqqqq:", quit_program
				mpi_finalize()
				sys.exit()

		quit_program = 0
		for runs_iter in range(0, no_of_viper_runs_analyzed_together):
			#print runs_iter
			#Generate/check directories for processing data
			#print "quit2:", rviper_iter, quit_program
			if (myid == main_node):
				independent_run_dir = masterdir + '/main%03d/run%03d/'%(rviper_iter, runs_iter) 	
				if os.path.exists(independent_run_dir):
					if os.path.exists(independent_run_dir + "log.txt"):
						#Check to see if this run has finished
						if string_found_in_file("Finish VIPER2", independent_run_dir + "log.txt"):
							count_already_finished_runs += 1
							#print "Continue"
							we_are_not_inbetween_completed_rrr_viper_runs = not ((count_already_finished_runs %  no_of_viper_runs_analyzed_together) == 0)
							# check if 
							#dir_len = mpi_bcast((-1)*(we_are_not_inbetween_completed_rrr_viper_runs),1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
							dir_len = bcast_number_to_all((-1)*(we_are_not_inbetween_completed_rrr_viper_runs))
							continue
					#Need to restart this independent run by doing rm -rf in run%02d directory 
					cmd = "{} {}".format("rm -rf", independent_run_dir)
					cmdexecute(cmd)
				cmd = "{} {}".format("mkdir -p", independent_run_dir)
				cmdexecute(cmd)
				dir_len = len(independent_run_dir)
				#print "dir_len=", dir_len
			else:
				dir_len = 0
				independent_run_dir = ""
			
			#dir_len = mpi_bcast(dir_len,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
			dir_len = bcast_number_to_all(dir_len)
			# in case all viper runs are done already for the current iteration, need to run identify_outliers program
			if (dir_len == 0):
				identify_outliers(myid, main_node, rviper_iter, no_of_viper_runs_analyzed_together, masterdir, bdb_stack_location)
				continue

			# Main node sent info to skip this iteration
			if (dir_len == -1):
				continue

			# Check if this iteration we need to update dbd and hdf


			independent_run_dir = mpi_bcast(independent_run_dir,dir_len,MPI_CHAR,main_node,MPI_COMM_WORLD)
			independent_run_dir = string.join(independent_run_dir,"")

			import global_def
			global_def.LOGFILE =  os.path.join(independent_run_dir, global_def.LOGFILE)
			#if(myid == main_node):
				#print independent_run_dir, "ZZZZ", global_def.LOGFILE
			mpi_barrier(MPI_COMM_WORLD)


			#mpi_finalize()
			#sys.exit()

			if independent_run_dir[-1] != "/":
				independent_run_dir += "/"

			
			#if(myid == main_node):
				#print independent_run_dir, "XXXXX", global_def.LOGFILE
			log.prefix = independent_run_dir
			
			if (rviper_iter == 1):
				if len(args) > 2:
					ref_vol = get_im(args[2])
				else:
					ref_vol = None
			else:
				#ref_vol = get_im(masterdir + "centered_%03d.hdf"%(rviper_iter))
				ref_vol = None

			#options.user_func = user_functions.factory[options.function]
			options.user_func = do_volume 

			#################print_with_time_info("Starting independent iteration %02d, by myid = %02d" % (runs_iter, myid))
			#time.sleep(10)

			#if (myid == main_node):
				#cmd = "{} {}".format("cp ~/log.txt ", independent_run_dir)
				#cmdexecute(cmd)
				#cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "param%03d.txt"%runs_iter)
				#cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "params.txt")
				#cmdexecute(cmd)



			####out_params, out_vol, out_peaks = multi_shc(all_projs, subset, no_of_viper_runs_analyzed_together, options, mpi_comm=mpi_comm, log=log, ref_vol=ref_vol)

			out_params, out_vol, out_peaks = multi_shc(all_projs, subset, no_of_shc_runs_analyzed_together, options, mpi_comm=mpi_comm, log=log, ref_vol=ref_vol)

			#print gethostname(), "mpi_finalize"
			#mpi_finalize()
			#import sys
			#sys.exit()

		

		identify_outliers(myid, main_node, rviper_iter, no_of_viper_runs_analyzed_together, masterdir, bdb_stack_location)

		# save current volumes, skip them if they are already there

		#calculate_volumes_after_rotation_and_save_them(mainoutputdir, bdb_stack_location, myid, mpi_size):


	print gethostname(), "mpi_finalize"
	mpi_finalize()



def get_latest_directory_increment_value(directory_location, directory_name):
	import os
	dir_count = 1
	while os.path.isdir(directory_location + directory_name + "%03d"%(dir_count)):
		dir_count += 1
	if dir_count == 1:
		return 1
	return dir_count - 1


def identify_outliers(myid, main_node, rviper_iter, no_of_viper_runs_analyzed_together, masterdir, bdb_stack_location):
	from utilities import wrap_mpi_bcast
	quit_program = 0
	if (myid == main_node):
		#if not found_outliers(outlier_percentile, rviper_iter, no_of_viper_runs_analyzed_together, masterdir):
		if not found_outliers(95, rviper_iter, no_of_viper_runs_analyzed_together, masterdir, bdb_stack_location):
			quit_program = 1
			cmd = "{} {} {}".format("mkdir ", masterdir, "Converged")
			cmdexecute(cmd)
			# update database files 

	#print "quit_program IO before:", quit_program
	quit_program = mpi_bcast(quit_program, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
	#print "quit_program IO after:", quit_program
	if quit_program > 0:
		print "mpi_finalize"
		mpi_finalize()
		sys.exit()


def found_outliers(outlier_percentile, rviper_iter, no_of_viper_runs_analyzed_together, masterdir,  bdb_stack_location):
	# sxheader.py bdb:nj  --consecutive  --params=OID
	import numpy as np
	from utilities import read_text_row, write_text_file, write_text_row
	from multi_shc import find_common_subset_3
	from pixel_error import rotate_angleset_to_match

	mainoutputdir = masterdir + "/main%03d/"%(rviper_iter)

	print "identify_outliers"
	projs = []
	for i1 in range(0,no_of_viper_runs_analyzed_together):
		projs.append(read_text_row(mainoutputdir + "run%03d"%(i1) + "/params.txt"))
		
##########  just for testing
	#if (rviper_iter == 1):
		#dat = EMData.read_images(bdb_stack_location)
	#else:
		#dat = EMData.read_images(bdb_stack_location + "_%03d"%(rviper_iter))

	#if (rviper_iter > 1):
		#for i1 in range(0,no_of_viper_runs_analyzed_together):
			#projs[i1] = projs[i1][:len(dat)]
##########  just for testing

	indi=[None,None,[None,None,None]];oo=[];th = 0.5
	while( len(indi[2]) < len(projs[0]) ):
			indi=find_common_subset_3(projs,th)
			oo.append([th,len(indi[2])])
			th+=0.5

	ou = [oo[i][0] for i in xrange(len(oo))]
	ot = [oo[i][1] for i in xrange(len(oo))]
	write_text_file([ou,ot], mainoutputdir + "pihis.txt")

	percentile_index = int(np.percentile(range(len(projs[0])), outlier_percentile))

	th = 1.0e23
	index_outliers = []
	index_keep_images = range(len(projs[0]))
	while( len(projs[0]) > percentile_index):
		indi=find_common_subset_3(projs,th)
		m = max(indi[2])
		l = indi[2].index(m)
		index_outliers.append(index_keep_images[l])
		del index_keep_images[l]
		for k in xrange(len(projs)):
			del projs[k][l]

	index_outliers.sort()
	#d = set(range(len(projs[0]))) - set(index_outliers)
	#index_keep_images = [q for q in d]
	index_keep_images.sort()

	#if len(index_outliers) < 3:
		#return False

	


	if (rviper_iter == 1):
		cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir +  "index_outliers.txt")
		cmdexecute(cmd)
		cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir +  "index_keep_images.txt")
		cmdexecute(cmd)
	else:
		cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location + "_%03d"%(rviper_iter), "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir  +  "index_outliers.txt")
		cmdexecute(cmd)
		cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location + "_%03d"%(rviper_iter), "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir +  "index_keep_images.txt")
		cmdexecute(cmd)


	dat = EMData.read_images(bdb_stack_location + "_%03d"%(rviper_iter + 1))
	write_text_file([dat[i].get_attr("original_image_index")  for i in index_outliers],mainoutputdir + "index_outliers.txt")
	write_text_file([dat[i].get_attr("original_image_index")  for i in index_keep_images],mainoutputdir + "index_keep_images.txt")

	#write_text_file(index_outliers, mainoutputdir + "this_iteration_index_outliers.txt")
	#write_text_file(index_keep_images, mainoutputdir + "this_iteration_index_keep_images.txt")


		#cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location , "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir  +  "index_outliers.txt")
		#cmdexecute(cmd)
		#cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location , "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir +  "index_keep_images.txt")
		#cmdexecute(cmd)
	#print "e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir,  index_outliers
	#print "e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter + 1), "--list=" + mainoutputdir, index_keep_images
	#print "e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter), "--list=" + mainoutputdir, index_outliers
	#print "e2bdb.py ", bdb_stack_location, "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter), "--list=" + mainoutputdir, index_keep_images

	print "index_outliers::", index_outliers
	print "index_keep_images::", index_keep_images
	#print "Length of dat (outliers):", len(dat)
	#print "Type of dat[0] (outliers):", type(dat[0]), type(dat[len(dat)-1]),
	#print "Length of gdat (keep_images):", len(gdat)
	#print "Type of gdat[0] (keep_images):", type(gdat[0]), type(gdat[len(gdat)-1]),
	
	#return True

	#write_text_file([dat[i].get_attr("OID")  for i in range(len(dat))],mainoutputdir + "rejects.txt")

	
	#for i in xrange(len(gdat)):  gdat[i].write_image(bdb_stack_location + "_%03d"%(rviper_iter + 1),i)
	#for i in xrange(len(gdat)):  gdat[i].write_image(masterdir + "centered_%03d.hdf"%(rviper_iter + 1),i)


	
	# reduce the set of parameters to kept images
	#for i1 in range(0,no_of_viper_runs_analyzed_together):
	#	projs[i1] = [projs[i1][i2] for i2 in index_keep_images]
	
	# apply rotation to params (reduced set) using volume#1 as a reference

	for i1 in range(1,no_of_viper_runs_analyzed_together):
		projs[i1] = rotate_angleset_to_match(projs[0], projs[i1])

	# write param files
	for i1 in range(0,no_of_viper_runs_analyzed_together):
		write_text_row(projs[i1], mainoutputdir + "run%03d"%(i1) + "/rotated_reduced_params.txt")
		write_text_file(projs[i1], mainoutputdir + "run%03d"%(i1) + "/rotated_reduced_params_not_row.txt")


	return True


def calculate_volumes_after_rotation_and_save_them(rviper_iter, masterdir, bdb_stack_location, myid, nproc):

	#compute 3D reconstruction and save volumes
	from reconstruction import recons3d_4nn_ctf, recons3d_4nn

	mainoutputdir = masterdir + "/main%03d/"%(rviper_iter)


	#  create a vstack from input stack to the local stack in masterdir
	#  Stack name set to default
	stack = "bdb:"+masterdir+"/rdata"
	# Initialization
	if(myid == 0):
		if keepchecking:
			if(not os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):
				cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+stack)
				cmdexecute(cmd)
				cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", stack)
				cmdexecute(cmd)
				keepchecking = False
		total_stack = EMUtil.get_image_count(stack)
		junk = get_im(stack)
		nnxo = junk.get_xsize()
		del junk
	else:
		total_stack = 0
		nnxo = 0

	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
	nnxo        = bcast_number_to_all(nnxo, source_node = main_node)





	# prepare names of input file names
	partids = [None]*6
	for procid in xrange(6):  partids[procid] = os.path.join(mainoutputdir,"X%01d.txt"%procid)
	partstack = [None]*6
	for procid in xrange(6):  partstack[procid] = os.path.join(mainoutputdir,"initparamsX%01d.txt"%procid)
	#stack = "bdb:"+masterdir+"/rdata"
	stack = bdb_stack_location + "_%03d"%(rviper_iter + 1)

	keepchecking = 1 
	for procid in xrange(6):
		outvol = os.path.join(mainoutputdir,"inivol%01d.hdf"%procid)
		doit, keepchecking = checkstep(outvol, keepchecking, myid, main_node=0)

		if not doit:
			continue

		print "-----------", procid
		from multi_shc import do_volume
		projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)

		if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		del projdata
		if( myid == main_node):
			vol.write_image(outvol)
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(  line,"Generated inivol #%01d "%procid)
		del vol

	#mpi_barrier(MPI_COMM_WORLD)

	return

def getindexdata(stack, partids, partstack, myid, nproc):
	from utilities import read_text_file
	lpartids  = map(int, read_text_file(partids) )
	ndata = len(lpartids)
	partstack = read_text_row(partstack)
	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1			
	else:
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	data = EMData.read_images(stack, lpartids)
	for i in xrange(len(partstack)):  set_params_proj(data[i], partstack[i])
	return data



def checkstep(item, keepchecking, myid, main_node):
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







if __name__=="__main__":
	main()


#start.sh mpirun --npernode 16 --host n0,n1,n2,n3,n4,n5  -tag-output -np 96 sxrviper.py bdb:centered  --doga=0.3   --L2threshold=0.05 --maxit2=50 --ou=30  --xr=1  --center=0


##mpirun --npernode 16 --host bmbpccl1 -tag-output -np 1 ./sxrviper.py  bdb:centered master2014_12_20__16_14_22
##mpirun --npernode 16 --host bmbpccl1 -tag-output -np 1 ./sxrviper.py  bdb:centered my01 

#rm -rf my01 && cp -r  master2014_12_20__16_14_22 my01
#mpirun --npernode 16 --host bmbpccl1 -tag-output -np 3 ./sxrviper.py  bdb:centered my01

# nohup mpirun --npernode 16 --host n0,n1 -tag-output -np 32 ./sxrviper.py  bdb:centered --niruns 3 --nruns 4 --doga=0.3 --nruns=5  --L2threshold=0.05 --maxit2=50 --ou=30  --xr=1  --center=0  --function="[.,paus_302,ali3d_reference3]" > fileN mylog &


#nohup mpirun --npernode 16 --host n0,n1 -tag-output -np 32 ./sxrviper.py  bdb:centered  --niruns 3 --nruns 4 --doga=0.3 --nruns=4  --L2threshold=0.05 --maxit2=50 --ou=30  --xr=1  --center=0  --function="[.,paus_302,ali3d_reference3]" >$(fileN mylog) &


#nohup mpirun --npernode 16 --host n0,n1 -tag-output -np 32 sxrviper.py  bdb:centered  --niruns 4 --nruns 4 --doga=0.3 --nruns=4  --L2threshold=0.05 --maxit2=50 --ou=30  --xr=1  --center=0  >$(fileN mylog) &


# start.sh mpirun --npernode 16 --host n0,n1 -tag-output -np 32 sxrviper.py \
# --niruns 10 --nruns 4 --doga=0.3   --L2threshold=0.05 --maxit2=50 --ou=30  --xr=1  --center=0

#kill -TERM  ` ps -ef | grep npernode  | grep -v grep | awk  '{print $2}'`

