#!/usr/bin/env python
from __future__ import print_function
#
#  09/09/2016
#  
#  CPU subgroup
#  10/27/2016  Added sigma2 updating in the first phased called PRIMARY
#  11/07       Shared refvol
#  10/28/2016 - Polar
#  11/18/2016 change in strategy
#  04/10/2017 - Enabled for one node
#  04/18/2017 - Introduce symclass to handle angles in a unified manner
#  01/21/2018 - Rationalize the code, particularly restart
#  07/21/2018 - Full size reconstruction after delta change
from builtins import range
"""
There are four ways to run the program:

1. Standard default run, starts from exhaustive searches, uses initial reference structure
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  bdb:sparx_stack vton1 mask15.hdf --sym=c5  --initialshifts  --radius=120  --mask3D=mask15.hdf    >1ovotn &

2. Restart after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  vton1 --radius=100 >2ovotn &

3. Local refinement, starts from user-provided orientation parameters, delta has to be <= 3.75
mpirun -np 64 --hostfile four_nodes.txt sxmeridien.py --local_refinement bdb:sparx_stack   vton3 --delta=1.875 --xr=2.0  --inires=5.5  --sym=c5  --radius=120  --mask3D=mask15.hdf >5ovotn &

4. Restart of local refinement after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py --local_refinement  vton3  --xr=0.6 >6ovotn &


"""

"""
Exhaustive:
          delta        projdir         projdir*npsi
         15.00000        189            4536  
          7.50000        756           36288  
          3.75000       3025          290400  
          1.87500      12101         2323392  
          0.93750      48405        18587520  
          0.46875     193622       148701696  
          0.23438     774489      1189615104  
          0.11719    3097958      9516926976  

Assume image 192x192
         15.00000        189       0.028 GB            4536          0.669 GB 
          7.50000        756       0.111 GB           36288          5.351 GB 
          3.75000       3025       0.446 GB          290400         42.821 GB 
          1.87500      12101       1.784 GB         2323392        342.598 GB 
          0.93750      48405       7.138 GB        18587520       2740.841 GB 
          0.46875     193622      28.551 GB       148701696      21926.957 GB 
          0.23438     774489     114.203 GB      1189615104     175415.885 GB 
          0.11719    3097958     456.812 GB      9516926976    1403327.984 GB 


Local with an = 6*delta
         15.00000         96            1248  
          7.50000        113            1469  
          3.75000        115            1495  
          1.87500        115            1495  
          0.93750        117            1521  
          0.46875        117            1521  
          0.23438        117            1521  
          0.11719        115            1495  

Local with an = 12*delta
         15.00000        189            4725  
          7.50000        377            9425  
          3.75000        446           11150  
          1.87500        461           11525  
          0.93750        470           11750  
          0.46875        464           11600  
          0.23438        463           11575  
          0.11719        463           11575  


"""
"""
08/14/2018
Normalization issues:
#This is for real space filter
nx = 1024
p1 = model_gauss_noise(1.0,nx,nx)
mask = Util.unrollmask(nx,nx)
for j in range(nx//2,nx):  mask[0,j]=1.0
#This is for Fourier valid region
m = Util.unrollmask(nx,nx)
p2 = fft(p1)
fp1=fft(Util.mulnclreal(p2,mask))
Util.innerproduct(p2,p2,m)/(nx*nx/2) = Util.innerproduct(fp1,fp1,None)
"""

from EMAN2 	import *
from sp_sparx 	import *
from EMAN2  import EMNumPy
from sp_logger import Logger, BaseLogger_Files
import sp_global_def
from sp_global_def import sxprint, ERROR

import sp_utilities

import sp_user_functions
from sp_global_def import *

import mpi
from mpi   	import  *
from math  	import  *
from random import *
import numpy as np
import shutil


import os
import sys
import subprocess
import string
import json
from   sys 	import exit
from   time import localtime, strftime, sleep, time

import sp_helix_sphire
import sp_helix_fundamentals

global Tracker, Blockdata
global  target_theta, refang

mpi.mpi_init(0, [])

Tracker   = {}
Blockdata = {}
#  MPI stuff
Blockdata["nproc"]              = mpi_comm_size(MPI_COMM_WORLD)
Blockdata["myid"]               = mpi_comm_rank(MPI_COMM_WORLD)
Blockdata["main_node"]          = 0
Blockdata["shared_comm"]		= mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,  0, MPI_INFO_NULL)
Blockdata["myid_on_node"]		= mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"] = mpi_comm_size(Blockdata["shared_comm"])
masters_from_groups_vs_everything_else_comm = mpi_comm_split(MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"])
Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = get_colors_and_subsets(Blockdata["main_node"], MPI_COMM_WORLD, Blockdata["myid"], \
		Blockdata["shared_comm"], Blockdata["myid_on_node"], masters_from_groups_vs_everything_else_comm)
#  We need two nodes for processing of volumes
if( Blockdata["no_of_groups"] > 1 ):  Blockdata["node_volume"] = [Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take two last nodes
else: Blockdata["node_volume"] = [0,0]
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
Blockdata["main_shared_nodes"]	= [Blockdata["node_volume"][0]*Blockdata["no_of_processes_per_group"],Blockdata["node_volume"][1]*Blockdata["no_of_processes_per_group"]]
# end of Blockdata
sp_global_def.BATCH = True
sp_global_def.MPI = True

def create_subgroup():
	# select a subset of myids to be in subdivision
	if( Blockdata["myid_on_node"] < Blockdata["ncpuspernode"] ): submyids = [Blockdata["myid"]]
	else:  submyids = []

	submyids = wrap_mpi_gatherv(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	submyids = wrap_mpi_bcast(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ): sxprint(submyids)
	world_group = mpi_comm_group(MPI_COMM_WORLD)
	subgroup = mpi_group_incl(world_group,len(submyids),submyids)
	#print(" XXX world group  ",Blockdata["myid"],world_group,subgroup)
	Blockdata["subgroup_comm"] = mpi_comm_create(MPI_COMM_WORLD, subgroup)
	mpi_barrier(MPI_COMM_WORLD)
	#print(" ZZZ subgroup  ",Blockdata["myid"],world_group,subgroup,subgroup_comm)

	Blockdata["subgroup_size"] = -1
	Blockdata["subgroup_myid"] = -1
	if (MPI_COMM_NULL != Blockdata["subgroup_comm"]):
		Blockdata["subgroup_size"] = mpi_comm_size(Blockdata["subgroup_comm"])
		Blockdata["subgroup_myid"] = mpi_comm_rank(Blockdata["subgroup_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi_barrier(MPI_COMM_WORLD)
	return


def create_zero_group():
	# select a subset of myids to be in subdivision, This is a group of all zero IDs on nodes, taken from isac2
	if( Blockdata["myid_on_node"] == 0 ): submyids = [Blockdata["myid"]]
	else:  submyids = []

	submyids = wrap_mpi_gatherv(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	submyids = wrap_mpi_bcast(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ): sxprint(submyids)
	world_group = mpi_comm_group(MPI_COMM_WORLD)
	subgroup = mpi_group_incl(world_group,len(submyids),submyids)
	#print(" XXX world group  ",Blockdata["myid"],world_group,subgroup)
	Blockdata["group_zero_comm"] = mpi_comm_create(MPI_COMM_WORLD, subgroup)
	mpi_barrier(MPI_COMM_WORLD)
	#print(" ZZZ subgroup  ",Blockdata["myid"],world_group,subgroup,subgroup_comm)

	Blockdata["group_zero_size"] = -1
	Blockdata["group_zero_myid"] = -1
	if (MPI_COMM_NULL != Blockdata["group_zero_comm"]):
		Blockdata["group_zero_size"] = mpi_comm_size(Blockdata["group_zero_comm"])
		Blockdata["group_zero_myid"] = mpi_comm_rank(Blockdata["group_zero_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	#Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi_barrier(MPI_COMM_WORLD)
	return

#if( Blockdata["subgroup_myid"] > -1 ):
#	dudu = [Blockdata["subgroup_myid"]]
#	dudu = wrap_mpi_gatherv(dudu, 0, Blockdata["subgroup_comm"])
#	if Blockdata["subgroup_myid"] == 0 :  sxprint("  HERE  ",dudu)

#we may want to free it in order to use different number of CPUs
#  create_subgroup()
#if( Blockdata["subgroup_myid"] > -1 ): mpi_comm_free(Blockdata["subgroup_comm"])


def params_changes( params, oldparams ):
	#  Indexes contain list of images processed - sorted integers, subset of the full range.
	#  params - contain parameters associated with these images
	#  Both lists can be of different sizes, so we have to find a common subset
	#  We do not compensate for random changes of grids.
	from sp_utilities    	import getang3
	from sp_utilities    	import rotate_shift_params
	from sp_pixel_error  	import max_3D_pixel_error
	from EMAN2        	import Vec2f
	from math 			import sqrt
	import sets

	n = len(params)
	anger       = 0.0
	shifter     = 0.0
	#  The shifter is given in the full scale displacement
	for i in range(n):
		shifter += (params[i][3] - oldparams[i][3] )**2 + (params[i][4] - oldparams[i][4] )**2
		anger += get_anger(params[i][0:3], oldparams[i][0:3])  # Symmetry is in Blockdata

	return round(anger/n,5), round(sqrt(shifter/2/n),5)

def assign_particles_to_groups(minimum_group_size = 10, asubset= None, name_tag='ptcl_source_image'):
	global Tracker, Blockdata
	from random import shuffle
	#  Input data does not have to be consecutive in terms of ptcl_source_image/filament_id or defocus
	#
	if not asubset:
		try:
			stmp  = EMUtil.get_all_attributes(Tracker["constants"]["stack"], name_tag)
			if Tracker["constants"]["CTF"]:
				defstmp = EMUtil.get_all_attributes(Tracker["constants"]["stack"],"ctf")
			else:
				defstmp = [-1.0]*len(stmp)
			for i in range(len(defstmp)): defstmp[i] = round(defstmp[i].defocus, 4)
		except:
			if Tracker["constants"]["CTF"]:
				stmp = EMUtil.get_all_attributes(Tracker["constants"]["stack"],"ctf")
				for i in range(len(stmp)):  stmp[i] = round(stmp[i].defocus, 4)
				defstmp = stmp[:]
			else:
				ERROR( "Either ptcl_source_image/filament/filament_id or ctf has to be present in the header." )
	else:
		try:
			stmp_junk = EMUtil.get_all_attributes(Tracker["constants"]["stack"], name_tag)
			stmp = [None]*len(asubset)
			for isub in range(len(asubset)): stmp[isub] = stmp_junk[asubset[isub]]
			if Tracker["constants"]["CTF"]:
				defstmp_junk = EMUtil.get_all_attributes(Tracker["constants"]["stack"],"ctf")
				defstmp  = [None]*len(asubset)
				for isub in range(len(asubset)): defstmp[isub] =  round(defstmp_junk[asubset[isub]].defocus, 4)
			else:  defstmp = [-1.0]*len(asubset)
		except:
			if Tracker["constants"]["CTF"]:
				stmp_junk = EMUtil.get_all_attributes(Tracker["constants"]["stack"], "ctf")
				stmp = [None]*len(asubset)
				defstmp = [-1.0]*len(asubset)
				for isub in range(len(asubset)): 
					stmp[isub] = stmp_junk[asubset[isub]]
					stmp[isub] = round(stmp[isub].defocus, 4)
				defstmp[:] = stmp[:]
			else:  
				ERROR( "Either ptcl_source_image/filament/filament_id or ctf has to be present in the header." )
	tt = [[stmp[i],i] for i in range(len(stmp))]
	tt.sort()
	tt.append([-1,-1])
	st = tt[0][0]
	sd = []
	occup = []
	groups = []
	ig = 0
	ib = 0
	for i in range(len(tt)):
		if(st != tt[i][0]):
			# create a group
			groups.append([tt[k][1] for k in range(ib,i)])
			sd.append([st,defstmp[tt[ib][1]]])
			occup.append(len(groups[ig]))
			groups[ig].sort()
			ib = i
			st = tt[i][0]
			ig += 1
	del tt, stmp, defstmp
	#print(" UUU  ", sd)
	#  [0]ID, [1]stamp, [2]defocus, [3]occupancy, [4]groups
	cross_reference_txt = [[[i] for i in range(len(sd))], [sd[i][0] for i in range(len(sd))], [sd[i][1] for i in range(len(sd))], [occup[i] for i in range(len(sd))], [groups[i] for i in range(len(sd))]]
	del occup, groups

	#  Remove small groups
	while(min(cross_reference_txt[3]) < minimum_group_size):
		#print("  minimum occupancy ",min(cross_reference_txt[3]),len(cross_reference_txt[3]))
		#  Find smallest group
		lax = minimum_group_size
		for i in range(len(cross_reference_txt[3])):
			if(lax > cross_reference_txt[3][i]):
				lax = cross_reference_txt[3][i]
				togo = i
		if Tracker["constants"]["CTF"]:
			# find nearest group by defocus
			sdef = 1.e23
			for i in range(len(cross_reference_txt[3])):
				if(i != togo):
					qt = abs(cross_reference_txt[2][i] - cross_reference_txt[2][togo])
					if(qt<sdef):
						target = i
						sdef = qt
		else:
			# find the next smallest
			lax = minimum_group_size
			for i in range(len(cross_reference_txt[3])):
				if(i != togo):
					if(lax > cross_reference_txt[3][i]):
						lax = cross_reference_txt[3][i]
						target = i
			
		#print("  merging groups  ",target,togo,cross_reference_txt[2][target],cross_reference_txt[2][togo],cross_reference_txt[3][target],cross_reference_txt[3][togo],len(cross_reference_txt[4][target]),len(cross_reference_txt[4][togo]))
		cross_reference_txt[2][target] = (cross_reference_txt[2][target]*sum(cross_reference_txt[0][target])+cross_reference_txt[2][togo]*sum(cross_reference_txt[0][togo]))
		cross_reference_txt[0][target] += cross_reference_txt[0][togo]
		cross_reference_txt[2][target] /= sum(cross_reference_txt[0][target])
		cross_reference_txt[3][target] += cross_reference_txt[3][togo]
		cross_reference_txt[4][target] += cross_reference_txt[4][togo]
		#print("  merged  ",cross_reference_txt[0][target],cross_reference_txt[3][target],len(cross_reference_txt[4][target]))

		#  remove the group
		for i in range(len(cross_reference_txt)):  del cross_reference_txt[i][togo]

	#  Sort as much as possible by the original particle number
	for i in range(len(cross_reference_txt[4])):
		cross_reference_txt[4][i].sort()

	temp = [[i,cross_reference_txt[4][i][0]] for i in range(len(cross_reference_txt[0]))]

	from operator import itemgetter
	temp.sort(key = itemgetter(1))

	cross_reference_txt = [[cross_reference_txt[j][temp[i][0]] for i in range(len(cross_reference_txt[0]))] for j in range(5)]

	write_text_row(cross_reference_txt[0], os.path.join(Tracker["constants"]["masterdir"],"main000","groupids.txt") )
	write_text_row([[sd[cross_reference_txt[0][i][j]][0] for j in range(len(cross_reference_txt[0][i]))]  for i in range(len(cross_reference_txt[0]))], os.path.join(Tracker["constants"]["masterdir"],"main000","micids.txt") )

	Tracker["constants"]["number_of_groups"] = len(cross_reference_txt[0])
	#  split into two chunks by groups
	lili = [[],list(range(Tracker["constants"]["number_of_groups"]))]
	shuffle(lili[1])
	lili[0] = lili[1][:len(lili[1])//2]
	lili[1] = lili[1][len(lili[1])//2:]
	lili[0].sort()
	lili[1].sort()

	#  Create output tables
	for iproc in range(2):
		write_text_row([cross_reference_txt[0][i] for i in lili[iproc]] , os.path.join(Tracker["constants"]["masterdir"],"main000","groupids_%03d.txt"%iproc) )
		write_text_row([[sd[cross_reference_txt[0][i][j]][0] for j in range(len(cross_reference_txt[0][i]))]  for i in lili[iproc]], os.path.join(Tracker["constants"]["masterdir"],"main000","micids_%03d.txt"%iproc) )
	del sd

	write_text_file([len(cross_reference_txt[4][i]) for i in range(len(cross_reference_txt[4]))], os.path.join(Tracker["constants"]["masterdir"],"main000","number_of_particles_per_group.txt") )

	q0 = []
	g0 = []
	q1 = []
	g1 = []
	for i in lili[0]:
		g0 += [i]*len(cross_reference_txt[4][i])
		q0 += cross_reference_txt[4][i]
	for i in lili[1]:
		g1 += [i]*len(cross_reference_txt[4][i])
		q1 += cross_reference_txt[4][i]
	Tracker["nima_per_chunk"] = [len(q0), len(q1)]

	### conversion
	
	#for iproc in range(2):
	#	if( Tracker["nima_per_chunk"][iproc] < Blockdata["nproc"] ):  ERROR("Number of particles per chunk smaller than the number of CPUs","assign_particles_to_groups",1,Blockdata["myid"])
	#write_text_file(q0, os.path.join(Tracker["constants"]["masterdir"],"main000","tchunk_0.txt") )
	write_text_file(g0, os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_0.txt") )
	#write_text_file(q1, os.path.join(Tracker["constants"]["masterdir"],"main000","tchunk_1.txt") )
	write_text_file(g1, os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_1.txt") )
	if asubset:
		sq0 =[None]*len(q0)
		sq1 =[None]*len(q1) 
		for iptl in range(len(sq0)): sq0[iptl] = asubset[q0[iptl]]
		for iptl in range(len(sq1)): sq1[iptl] = asubset[q1[iptl]]
		return sq0, sq1
	else: return  q0, q1
#  CONES support functions

def number_of_cones_to_delta(number_of_cones):

	if( number_of_cones == 1):  return Blockdata["symclass"][1][3] + 0.5, 1
	else:
		if( Blockdata["symclass"].sym[0] == "c"):
			t2 = 89.0
			for i in range(92,0,-1):
				a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
				a += [[(q[0]+90.0)%360., 180.-q[1],0] for q in a]
				nren = len(a)
				if(nren>number_of_cones): return float(i), nren
			while(True):
				i /= 2.0
				a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
				a += [[(q[0]+90.0)%360., 180.-q[1],0] for q in a]
				nren = len(a)
				if(nren>number_of_cones): return float(i), nren

		else:
			#t2 = Blockdata["symclass"].brackets[1][3] + 0.5
			for i in range(int(t2+1),0,-1):
				nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
				if(nren>number_of_cones): return float(i), nren
			while(True):
				i /= 2.0
				nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
				if(nren>number_of_cones): return float(i), nren

		ERROR( "number_of_cones_to_delta should not be here", 0 )
		return -1, 0

def find_assignments_of_refangles_to_angles(normals_set, ancor_angle, an):
	#  returns list of angles within normals_set that are within an angles from ancor_angle
	global  Blockdata

	nang1 = len(normals_set) - 1
	Blockdata["target_theta"] = ancor_angle[1]
	u1,u2 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
	u1 = int(u1+0.5)

	Blockdata["target_theta"] = max(0.0, ancor_angle[1] - 1.2*an )
	u3,u4 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
	u3 = int(u3+0.5)
	if(u3<10):  u3 = 0

	Blockdata["target_theta"] = min(Blockdata["symclass"].brackets[1][3], ancor_angle[1] + 1.2*an )
	u5,u6 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
	u5 = int(u5+0.5)
	if(u5>nang1-10):  u5 = nang1+1

	ancordir = angles_to_normals(Blockdata["symclass"].symmetry_neighbors([ancor_angle[:3]]))
	ltemp = Util.cone_dirs_f( normals_set[u3:u5], ancordir, an )
	###ltemp = cone_dirs_f( normals_set[u3:u5], ancordir, an )#Util.cone_dirs_f( normals_set[u3:u5], ancordir, an )
	#print(" us ",Blockdata["myid"],u1,u3,u5,m,"  ltemp ",ltemp)
	return  [qtemp+u3 for qtemp in ltemp]

def find_nearest_k_refangles_to_many_angles(normals_set, angles, delta, howmany):
	assignments = [-1]*len(angles)
	for i,ancor_angle in enumerate(angles):
		assignments[i] = find_nearest_k_refangles_to_angles(normals_set, ancor_angle, delta, howmany)
	return assignments

def find_nearest_k_refangles_to_angles(normals_set, ancor_angle, delta, howmany):
	#  returns list of angles within normals_set that are within an angles from ancor_angle
	global  Blockdata
	nang1 = len(normals_set) - 1
	qtl = 0
	bigger = 1.0
	while( qtl < howmany ):
		an = bigger*max(delta*(howmany//2),delta)
		Blockdata["target_theta"] = ancor_angle[1]
		u1,u2 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
		u1 = int(u1+0.5)

		Blockdata["target_theta"] = max(0.0, ancor_angle[1] - 1.1*an )
		u3,u4 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
		u3 = int(u3+0.5)
		if(u3<10):  u3 = 0

		Blockdata["target_theta"] = min(Blockdata["symclass"].brackets[1][3], ancor_angle[1] + 1.1*an )
		u5,u6 = goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
		u5 = int(u5+0.5)
		if(u5>nang1-10):  u5 = nang1+1
		qtl = len(normals_set[u3:u5])
		bigger *= 1.25
		
	if( Blockdata["symclass"].sym == "c1"):
		ancordir = getfvec(ancor_angle[0],ancor_angle[1])
		ltemp = Util.nearest_fang_select(normals_set[u3:u5], ancordir[0], ancordir[1], ancordir[2], howmany)
	else:
		ancordir = angles_to_normals(Blockdata["symclass"].symmetry_neighbors([ancor_angle]))
		ltemp = Util.nearest_fang_sym(ancordir, normals_set[u3:u5], len(ancordir), howmany)
	return  [qtemp+u3 for qtemp in ltemp]

def auxiliary_funcdef(xxx):
	global  Blockdata
	l = min(max(int(xxx+0.5),0),len(Blockdata["angle_set"])-1)
	#print l,abs(target_theta - refang[l][1])
	return  abs(Blockdata["target_theta"] - Blockdata["angle_set"][l][1])

#  ------

def compute_sigma(projdata, params, first_procid, dryrun = False, myid = -1, mpi_comm = -1):
	global Tracker, Blockdata
	# Input stack of particles with all params in header
	# Output: 1/sigma^2 and a dictionary
	#  It could be a parameter
	if( mpi_comm < 0 ): mpi_comm = MPI_COMM_WORLD
	npad = 1

	if  dryrun:
		#tsd = model_blank(nv + nv//2,len(sd), 1, 1.0)
		#tocp = model_blank(len(sd), 1, 1, 1.0)
		if( myid == Blockdata["main_node"] ):
			tsd = get_im(os.path.join(Tracker["previousoutputdir"],"bckgnoise.hdf"))
			tsd.write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf"))
			nnx = tsd.get_xsize()
			nny = tsd.get_ysize()
		else:
			nnx = 0
			nny = 0
		nnx = bcast_number_to_all(nnx, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)
		nny = bcast_number_to_all(nny, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)
		if( myid != Blockdata["main_node"] ):
			tsd = model_blank(nnx,nny, 1, 1.0)
		bcast_EMData_to_all(tsd, myid, source_node = Blockdata["main_node"], comm = mpi_comm)
		'''
		#  I am not sure whether what follows is correct.  This part should be recalculated upon restart
		Blockdata["accumulatepw"] = [[],[]]
		ndata = len(projdata)
		for i in range(ndata):
			if(i<first_procid):  iproc = 0 #  This points to correct procid
			else:                iproc = 1
			Blockdata["accumulatepw"][iproc].append([0.0]*200)
		'''

	else:

		if( myid == Blockdata["main_node"] ):
			ngroups = len(read_text_file(os.path.join(Tracker["constants"]["masterdir"],"main000", "groupids.txt")))
		else: ngroups = 0
		ngroups = bcast_number_to_all(ngroups, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)

		ndata = len(projdata)
		nx = Tracker["constants"]["nnxo"]
		mx = npad*nx
		nv = mx//2+1
		"""
		#  Inverted Gaussian mask
		invg = model_gauss(Tracker["constants"]["radius"],nx,nx)
		invg /= invg[nx//2,nx//2]
		invg = model_blank(nx,nx,1,1.0) - invg
		"""

		mask = model_circle(Tracker["constants"]["radius"],nx,nx)
		tsd = model_blank(nv + nv//2, ngroups)

		#projdata, params = getalldata(partstack, params, myid, Blockdata["nproc"])
		'''
		if(myid == 0):  ndata = EMUtil.get_image_count(partstack)
		else:           ndata = 0
		ndata = bcast_number_to_all(ndata)
		if( ndata < Blockdata["nproc"]):
			if(myid<ndata):
				image_start = myid
				image_end   = myid+1
			else:
				image_start = 0
				image_end   = 1
		else:
			image_start, image_end = MPI_start_end(ndata, Blockdata["nproc"], myid)
		#data = EMData.read_images(stack, range(image_start, image_end))
		if(myid == 0):
			params = read_text_row( paramsname )
			params = [params[i][j]  for i in range(len(params))   for j in range(5)]
		else:           params = [0.0]*(5*ndata)
		params = bcast_list_to_all(params, myid, source_node=Blockdata["main_node"])
		params = [[params[i*5+j] for j in range(5)] for i in range(ndata)]
		'''
		if(Blockdata["accumulatepw"] == None):
			Blockdata["accumulatepw"] = [[],[]]
			doac = True
		else:  doac = False
		tocp = model_blank(ngroups)
		tavg = model_blank(nx,nx)
		for i in range(ndata):  # apply_shift; info_mask; norm consistent with get_shrink_data
			indx = projdata[i].get_attr("particle_group")
			phi,theta,psi,sx,sy = params[i][0],params[i][1],params[i][2],params[i][3],params[i][4]
			sx, sy = reduce_shifts(sx, sy, projdata[i])
			stmp = cyclic_shift( projdata[i], sx, sy)
			st = get_image_statistics(stmp, mask, False)
			stmp -=st[0]
			stmp /=st[1]
			temp = cosinemask(stmp, radius = Tracker["constants"]["radius"], s = 0.0)
			Util.add_img(tavg, temp)
			sig = Util.rotavg_fourier( temp )
			#sig = rops(pad(((cyclic_shift( projdata[i], int(sx), int(round(sy)) ) - st[0])/st[1]), mx,mx,1,0.0))
			#sig = rops(pad(((cyclic_shift(projdata, int(round(params[i][-2])), int(round(params[i][-1])) ) - st[0])/st[1])*invg, mx,mx,1,0.0))
			for k in range(nv):
				tsd.set_value_at(k,indx,tsd.get_value_at(k,indx)+sig[k])
			'''
			if doac:
				if(i<first_procid):  iproc = 0 #  This points to correct procid
				else:                iproc = 1
				Blockdata["accumulatepw"][iproc].append(sig[nv:]+[0.0])  # add zero at the end so for the full size nothing is added.
			'''
			tocp[indx] += 1

		####for lll in range(len(Blockdata["accumulatepw"])):  sxprint(myid,ndata,lll,len(Blockdata["accumulatepw"][lll]))
		reduce_EMData_to_root(tsd, myid, Blockdata["main_node"],  mpi_comm)
		reduce_EMData_to_root(tocp, myid, Blockdata["main_node"], mpi_comm)
		reduce_EMData_to_root(tavg, myid, Blockdata["main_node"], mpi_comm)
		if( myid == Blockdata["main_node"]):
			Util.mul_scalar(tavg, 1.0/float(sum(Tracker["nima_per_chunk"])))
			sig = Util.rotavg_fourier( tavg )
			#for k in range(1,nv):  sxprint("  BACKG  ",k,tsd.get_value_at(k,0)/tocp[0] ,sig[k],tsd.get_value_at(k,0)/tocp[0] - sig[k])
			tmp1 = [0.0]*nv
			tmp2 = [0.0]*nv
			for i in range(ngroups):
				for k in range(1,nv):
					qt = tsd.get_value_at(k,i)/tocp[i] - sig[k]
					if( qt > 0.0 ):	tmp1[k] = 2.0/qt
				#  smooth
				tmp1[0] = tmp1[1]
				tmp1[-1] = tmp1[-2]
				for ism in range(0):  #2
					for k in range(1,nv-1):  tmp2[k] = (tmp1[k-1]+tmp1[k]+tmp1[k+1])/3.0
					for k in range(1,nv-1):  tmp1[k] = tmp2[k]
				"""
				for k in range(6,nv):
					tsd.set_value_at(k,i,1.0/(tsd.get_value_at(k,i)/tocp[i]))  # Already inverted
				qt = tsd.get_value_at(6,i)
				for k in range(1,6):
					tsd.set_value_at(k,i,qt)
				"""
				#  We will keep 0-element the same as first tsd.set_value_at(0,i,1.0)
				for k in range(1,nv):
					tsd.set_value_at(k,i,tmp1[k])
				tsd.set_value_at(0,i,1.0)
			tsd.write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf"))
		bcast_EMData_to_all(tsd, myid, source_node = 0, comm = mpi_comm)
	nnx = tsd.get_xsize()
	nny = tsd.get_ysize()
	Blockdata["bckgnoise"] = []
	for i in range(nny):
		prj = model_blank(nnx)
		for k in range(nnx): prj[k] = tsd.get_value_at(k,i)
		Blockdata["bckgnoise"].append(prj)  #  1.0/sigma^2
	return
	#return Blockdata["bckgnoise"]#tsd, sd#, [int(tocp[i]) for i in range(len(sd))]

def getindexdata(partids, partstack, particle_groups, original_data=None, small_memory=True, nproc =-1, myid = -1, mpi_comm = -1):
	global Tracker, Blockdata
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	if( mpi_comm < 0 ):  mpi_comm = MPI_COMM_WORLD
	from sp_applications import MPI_start_end
	#  parameters
	if( myid == 0 ):  partstack = read_text_row(partstack)
	else:  			  partstack = 0
	partstack = wrap_mpi_bcast(partstack, 0, mpi_comm)
	#  particles IDs
	if( myid == 0 ):  partids = read_text_file(partids)
	else:          	  partids = 0
	partids = wrap_mpi_bcast(partids, 0, mpi_comm)
	#  Group assignments
	if( myid == 0 ):	group_reference = read_text_file(particle_groups)
	else:          		group_reference = 0
	group_reference = wrap_mpi_bcast(group_reference, 0, mpi_comm)

	im_start, im_end = MPI_start_end(len(partstack), nproc, myid)
	partstack = partstack[im_start:im_end]
	partids   = partids[im_start:im_end]
	group_reference = group_reference[im_start:im_end]
	'''
	particles_on_node = []
	parms_on_node     = []
	for i in range( group_start, group_end ):
		particles_on_node += lpartids[group_reference[i][2]:group_reference[i][3]+1]  #  +1 is on account of python idiosyncrasies
		parms_on_node     += partstack[group_reference[i][2]:group_reference[i][3]+1]


	Blockdata["nima_per_cpu"][procid] = len(particles_on_node)
	#ZZprint("groups_on_thread  ",Blockdata["myid"],procid, Tracker["groups_on_thread"][procid])
	#ZZprint("  particles  ",Blockdata["myid"],Blockdata["nima_per_cpu"][procid],len(parms_on_node))
	'''
	"""
            17            28            57            84    5
            18            14            85            98    6
            19            15            99           113    7
            25            20           114           133    8
            29             9           134           142    9

	"""

	if( original_data == None or small_memory):
		original_data = EMData.read_images(Tracker["constants"]["stack"], partids)
		for im in range( len(original_data) ):
			try:
				original_data[im].set_attr("particle_group", group_reference[im])
			except:
				print(len(original_data), len(group_reference), len(partids))
				raise
	return original_data, partstack, [im_start, im_end]

def get_shrink_data(nxinit, procid, original_data = None, oldparams = None, \
					return_real = False, preshift = False, apply_mask = True, nonorm = False, nosmearing = False, npad = 1):
	global Tracker, Blockdata
	"""
	This function will read from stack a subset of images specified in partids
	   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	So, the lengths of partids and partstack are the same.
	  The read data is properly distributed among MPI threads.
	
	Flow of data:
	1. Read images, if there is enough memory, keep them as original_data.
	2. Read current params
	3.  Apply shift
	4.  Normalize outside of the radius
	5.  Do noise substitution and cosine mask.  (Optional?)
	6.  Shrink data.
	7.  Apply CTF.
	
	"""
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt
	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, preshift) )
	#  Preprocess the data
	mask2D  	= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima 		= len(original_data)
	shrinkage 	= nxinit/float(Tracker["constants"]["nnxo"])


	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	radius 	= int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm    	= float(nxinit-(nxinit//2+1) - radius)
	txl    	= float(radius - nxinit//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	Blockdata["accumulatepw"][procid] = [None]*nima
	data = [None]*nima
	for im in range(nima):

		if Tracker["mainiteration"] <=1:
			phi,theta,psi,sx,sy, = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4]
			wnorm  = 1.0 
		else:
			phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]
			
		'''
		if preshift:
			data[im] = fshift(original_data[im], sx, sy)
			sx = 0.0
			sy = 0.0
		'''

		if preshift:
			#data[im] = fshift(original_data[im], sx, sy)
			if nosmearing:
				data[im]  = fshift(original_data[im], sx, sy)
				oldparams[im][3] = 0.0
				oldparams[im][4] = 0.0
			else:
				sx, sy = reduce_shifts(sx, sy, original_data[im])
				data[im]  = cyclic_shift(original_data[im],sx,sy)
				#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
				oldparams[im][3] = sx
				oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  data[im] = original_data[im].copy()

		st = get_image_statistics(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if data[im].get_attr_default("bckgnoise", None) :  data[im].delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if Blockdata["bckgnoise"]:
			if apply_mask:
				if Tracker["constants"]["hardmask"]:
					data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"])
				else:
					bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = fft(filt_table(bckg, oneover[data[im].get_attr("particle_group")]))
					#  Normalize bckg noise in real space, only region actually used.
					st = get_image_statistics(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			#  if no bckgnoise, do simple masking instead
			if apply_mask:  data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"] )
		#  resample will properly adjusts shifts and pixel size in ctf
		#data[im] = resample(data[im], shrinkage)
		#  return Fourier image
		#if npad> 1:  data[im] = pad(data[im], Tracker["constants"]["nnxo"]*npad, Tracker["constants"]["nnxo"]*npad, 1, 0.0)

		#  Apply varadj
		if not nonorm:
			Util.mul_scalar(data[im], Tracker["avgvaradj"][procid]/wnorm)
			#sxprint(Tracker["avgvaradj"][procid]/wnorm)		

		#  FT
		data[im] = fft(data[im])
		sig = Util.rotavg_fourier( data[im] )
		Blockdata["accumulatepw"][procid][im] = sig[len(sig)//2:]+[0.0]

		if Tracker["constants"]["CTF"] :
			data[im] = fdecimate(data[im], nxinit*npad, nxinit*npad, 1, False, False)
			ctf_params = original_data[im].get_attr("ctf")
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
			#if Tracker["applyctf"] :  #  This should be always False
			#	data[im] = filt_ctf(data[im], ctf_params, dopad=False)
			#	data[im].set_attr('ctf_applied', 1)
			#else:
			data[im].set_attr('ctf_applied', 0)
			if return_real :  data[im] = fft(data[im])
		else:
			ctf_params = original_data[im].get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				data[im].set_attr('ctf', ctf_params)
				data[im].set_attr('ctf_applied', 0)
			data[im] = fdecimate(data[im], nxinit*npad, nxinit*npad, 1, True, False)
			apix = Tracker["constants"]["pixel_size"]
			data[im].set_attr('apix', apix/shrinkage)

		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		if not return_real:
			data[im].set_attr("padffted",1)
		data[im].set_attr("npad",npad)
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
			###  Do not adjust the values, we try to keep everything in the same Fourier values.
			data[im].set_attr("bckgnoise", [temp[i] for i in range(temp.get_xsize())])
	return data

def subdict(d,u):
	# substitute values in dictionary d by those given by dictionary u
	for q in u:  d[q] = u[q]

def get_anger(angle1, angle2):
	from math import acos, degrees
	from sp_utilities import lacos
	from sp_fundamentals import rotmatrix
	A1 = rotmatrix(angle1[0],angle1[1],angle1[2])
	ar = Blockdata["symclass"].symmetry_related(angle2)
	axes_dis_min = 1.0e23
	for q in ar:
		A2 = rotmatrix(q[0],q[1],q[2])
		axes_dis = 0.0
		for i in range(3):
			axes_dis += lacos(A1[i][0]*A2[i][0] + A1[i][1]*A2[i][1] + A1[i][2]*A2[i][2])
		axes_dis_min = min(axes_dis_min, axes_dis/3.0)
	return axes_dis_min

def checkstep(item, keepchecking):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
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
	doit = bcast_number_to_all(doit, source_node = Blockdata["main_node"])
	return doit, keepchecking

def out_fsc(f):
	global Tracker, Blockdata
	sxprint(" ")
	sxprint("  driver FSC  after  iteration#%3d"%Tracker["mainiteration"])
	sxprint("  %4d        %7.2f         %7.3f"%(0,1000.00,f[0]))
	for i in range(1,len(f)):
		sxprint("  %4d        %7.2f         %7.3f"%(i,Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(i),f[i]))
	sxprint(" ")

def get_refangs_and_shifts():
	global Tracker, Blockdata

	refang = Blockdata["symclass"].even_angles(
		delta=Tracker["delta"],
		theta1=Tracker['theta_min'],
		theta2=Tracker['theta_max'],
		method=Tracker['constants']['even_angle_method'],
		)
	coarse = Blockdata["symclass_coarse"].even_angles(
		delta=2*Tracker["delta"],
		theta1=Tracker['theta_min'],
		theta2=Tracker['theta_max'],
		method=Tracker['constants']['even_angle_method'],
		)
	refang = Blockdata["symclass"].reduce_anglesets( rotate_params(refang, [-0.5*Tracker["delta"], -0.5*Tracker["delta"], -0.5*Tracker["delta"]]) )

	"""
	if(Tracker["delta"] == 15.0):  refang = read_text_row("refang15.txt")
	elif(Tracker["delta"] == 7.5):  refang = read_text_row("refang7p5.txt")
	elif(Tracker["delta"] == 3.75):  refang = read_text_row("refang3p75.txt")
	elif(Tracker["delta"] == 1.875):  refang = read_text_row("refang1p875.txt")
	elif(Tracker["delta"] == 0.9375):  refang = read_text_row("refang0p9375.txt")
	elif(Tracker["delta"] == 0.46875):  refang = read_text_row("refang0p46875.txt")
	"""

	k = int(ceil(Tracker["xr"]/Tracker["ts"]))
	radi = (Tracker["xr"]+Tracker["ts"])**2
	rshifts = []
	for ix in range(-k,k,1):
		six = ix*Tracker["ts"] + Tracker["ts"]/2
		for iy in range(-k,k,1):
			siy = iy*Tracker["ts"] + Tracker["ts"]/2
			if(six*six+siy*siy <= radi):
				rshifts.append( [six, siy] )

	ts_coarse = 2*Tracker["ts"]
	k = int(ceil(Tracker["xr"]/ts_coarse))
	radi = Tracker["xr"]*Tracker["xr"]
	coarse_shifts = []
	for ix in range(-k,k+1,1):
		six = ix*ts_coarse
		for iy in range(-k,k+1,1):
			siy = iy*ts_coarse
			if(six*six+siy*siy <= radi):
				coarse_shifts.append( [six, siy] )

	return refang, rshifts, coarse, coarse_shifts

def get_coarse_shifts():
	global Tracker, Blockdata

	k = int(ceil(Tracker["xr"]/Tracker["ts"]))
	radi = Tracker["xr"]*Tracker["xr"]

	indc = []
	for ix in range(-k,k+1,1):
		six = ix*Tracker["ts"]
		for iy in range(-k,k+1,1):
			siy = iy*Tracker["ts"]
			if(six*six+siy*siy <= radi):
				indc.append( [ix, iy] )

	cndc = []
	for ix in range(0,k+1,2):
		six = ix*Tracker["ts"]
		for iy in range(0,k+1,2):
			siy = iy*Tracker["ts"]
			if(six*six+siy*siy <= radi):
				cndc.append( [ix, iy] )

	for ix in range(1,len(cndc)):  cndc.append([-cndc[ix][0],-cndc[ix][1]])

	list_of_coarse_shifts = []
	for ix in range(len(cndc)):  list_of_coarse_shifts.append(indc.index(cndc[ix]))


	return list_of_coarse_shifts

def get_shifts_neighbors(rshifts, cs):
	"""
	  	output is a list of shift neighbors to cs within ts*1.5 distance\
	  	This yields at most five shifts
	"""
	shiftneighbors = []
	rad = Tracker["ts"]*1.5
	for l,q in enumerate(rshifts):
		if(get_dist(q, cs) <= rad): shiftneighbors.append(l)
	return 	shiftneighbors

def shakegrid(rshifts, qt):
	for i in range(len(rshifts)):
		rshifts[i][0] += qt
		rshifts[i][1] += qt

###----------------

def get_refvol(nxinit):
	ref_vol = get_im(Tracker["refvol"])
	nnn = ref_vol.get_xsize()
	if( nxinit != nnn ):
		ref_vol = fdecimate(ref_vol, nxinit, nxinit, nxinit, True, False)
	return ref_vol

def prepdata_ali3d(projdata, rshifts, shrink, method = "DIRECT"):
	global Tracker, Blockdata
	from sp_fundamentals 	import prepi
	from sp_morphology 	import ctf_img_real
	#  Data is NOT CTF-applied.
	#  Data is shrank, in Fourier format
	data = [[] for i in range(len(projdata))]
	if Tracker["constants"]["CTF"]:
		nx = projdata[0].get_ysize()
		ctfs = [ ctf_img_real(nx, q.get_attr('ctf')) for q in projdata ]
	else:  ctfs = None
	if Blockdata["bckgnoise"] :
		bckgnoise = [q.get_attr("bckgnoise") for q in projdata ]
	else:  bckgnoise = None
	for kl in range(len(projdata)-1,-1,-1):  #  Run backwards and keep deleting projdata, it is not needed anymore
		#  header shifts were shrank in get_shrink_data, shifts were already pre-applied, but we will leave the code as is.
		#phi, theta, psi, sxs, sys = get_params_proj(projdata[kl])
		particle_group = projdata[kl].get_attr("particle_group")
		ds = projdata[kl]
		for iq in rshifts:
			xx = iq[0]*shrink
			yy = iq[1]*shrink
			dss = fshift(ds, xx, yy)
			dss.set_attr("is_complex",0)
			'''
			if( method == "DIRECT" ):
				#dss = fshift(ds, xx+sxs, yy+sys)
				dss = fshift(ds, xx+sxs, yy+sys)
				dss.set_attr("is_complex",0)
			else:
				dss = fft(fshift(ds, x+sxs, yy+sys))
				dss,kb = prepi(dss)
			'''
			data[kl].append(dss)
		data[kl][0].set_attr("particle_group",particle_group)  #  Pass group number only in the first shifted copy.
		del projdata[kl]
	return data, ctfs, bckgnoise

def do3d(procid, data, newparams, refang, rshifts, norm_per_particle, myid, smearing = True, mpi_comm = -1):
	global Tracker, Blockdata

	#  Without filtration
	from sp_reconstruction import recons3d_trl_struct_MPI
	import os
	if (mpi_comm == -1): mpi_comm = MPI_COMM_WORLD
	
	if procid == 0:
		if(Blockdata["no_of_groups"] >1):
			if myid == Blockdata["nodes"][0] :
				if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): sxprint("tempdir exists")
				else: 
					try: os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
					except:  sxprint("tempdir exists")
		else:
			if myid == Blockdata["main_node"]:
				if not os.path.exists(os.path.join(Tracker["directory"],"tempdir")): 
					try: os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
					except: sxprint("tempdir exists")
				else: sxprint("tempdir exists")
	mpi_barrier(mpi_comm)
	"""
	tvol, tweight, trol = recons3d_4nnstruct_MPI(myid = Blockdata["subgroup_myid"], main_node = Blockdata["nodes"][procid], prjlist = data, \
											paramstructure = newparams, refang = refang, delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"],\
											upweighted = False, mpi_comm = mpi_comm, \
											target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	"""
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	if smearing:
		tvol, tweight, trol = recons3d_trl_struct_MPI(
			myid=myid,
			main_node=Blockdata["nodes"][procid],
			prjlist=data,
			paramstructure=newparams,
			refang=refang,
			rshifts_shrank=[[q[0]*shrinkage,q[1]*shrinkage] for q in rshifts],
			delta=Tracker["delta"],
			CTF=Tracker["constants"]["CTF"],
			upweighted=False,
			mpi_comm=mpi_comm,
			target_size=(2*Tracker["nxinit"]+3),
			avgnorm=Tracker["avgvaradj"][procid],
			norm_per_particle=norm_per_particle
			)
	else:
		tvol, tweight, trol = recons3d_trl_struct_MPI_nosmearing(
			myid=myid,
			main_node=Blockdata["nodes"][procid],
			prjlist=data,
			parameters=newparams,
			CTF=Tracker["constants"]["CTF"],
			upweighted=False,
			mpi_comm=mpi_comm,
			target_size=(2*Tracker["nxinit"]+3)
			)
	if(Blockdata["no_of_groups"] >1):
		if myid == Blockdata["nodes"][procid]:
			while not os.path.exists(os.path.join(Tracker["directory"],"tempdir")):  sleep(5)
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%01d_%03d.hdf"%(procid,   Tracker["mainiteration"])))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
			trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%01d_%03d.hdf"%(procid,   Tracker["mainiteration"])))
			#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			line = ''
			sxprint(line,"Executed successfully backprojection for group ",procid)
	else:
		if myid == Blockdata["main_node"]:
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(   Tracker["directory"], "tempdir", "tvol_%01d_%03d.hdf"%(procid,   Tracker["mainiteration"])))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
			trol.write_image(os.path.join(   Tracker["directory"], "tempdir", "trol_%01d_%03d.hdf"%(procid,   Tracker["mainiteration"])))
			#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			line = ''
			sxprint(line,"Executed successfully backprojection for group ",procid)
	mpi_barrier(mpi_comm)
	return  
	
def print_dict(dict,theme):
	#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	line = ''
	sxprint(line,theme)
	spaces = "                    "
	exclude = ["constants", "maxit", "nodes", "yr", "shared_comm", "bckgnoise", "myid", "myid_on_node", "accumulatepw",\
				"changed_delta", "prior"]
	for key, value in sorted( dict.items() ):
		pt = True
		for ll in exclude:
			if(key == ll):
				pt = False
				break
		if pt:  sxprint("                    => ", key+spaces[len(key):],":  ",value)

def stepone(tvol, tweight):
	global Tracker, Blockdata
	tvol.set_attr("is_complex",1)
	ovol = Util.shrinkfvol(tvol,2)
	owol = Util.shrinkfvol(tweight,2)
	if( Tracker["constants"]["symmetry"] != "c1" ):
		ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
		owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
	#print(info(ovol,Comment = " shrank ovol"))
	return Util.divn_cbyr(ovol,owol)

def steptwo(tvol, tweight, treg, cfsc = None, regularized = True):
	global Tracker, Blockdata
	nz = tweight.get_zsize()
	ny = tweight.get_ysize()
	nx = tweight.get_xsize()
	tvol.set_attr("is_complex",1)
	if regularized:
		nr = len(cfsc)
		limitres = 0
		for i in range(nr):
			cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
			#print( i,cfsc[i] )
			if( cfsc[i] == 0.0 ):
				limitres = i-1
				break
		if( limitres == 0 ): limitres = nr-2;
		ovol = reshape_1d(cfsc, nr, 2*nr)
		limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
		maxr2 = limitres**2
		for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
		ovol[0] = 1.0
		#print(" ovol  ", ovol)
		it = model_blank(2*nr)
		for i in range(2*nr):  it[i] = ovol[i]
		del ovol
		#  Do not regularize first four
		for i in range(5):  treg[i] = 0.0
		Util.reg_weights(tweight, treg, it)
		del it
	else:
		limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
		maxr2 = limitres**2
	#  Iterative weights
	if( Tracker["constants"]["symmetry"] != "c1" ):
		tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
		tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

	#  tvol is overwritten, meaning it is also an output
	Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])
	#  Either pad or window in F space to 2*nnxo
	nx = tvol.get_ysize()
	if( nx > 2*Tracker["constants"]["nnxo"] ):
		tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
	elif(nx < 2*Tracker["constants"]["nnxo"] ):
		tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)

	tvol = fft(tvol)
	tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol.set_attr("npad",2)
	tvol.div_sinc(1)
	tvol.del_attr("npad")
	tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol = cosinemask(tvol,Tracker["constants"]["nnxo"]//2-1-5,5, None) # clean artifacts in corners

	return tvol

def steptwo_mpi(tvol, tweight, treg, cfsc = None, regularized = True, color = 0):
	global Tracker, Blockdata
	
	n_iter = 10 # number of iterations for iterative process for doing map
	if( Blockdata["color"] != color ):  return model_blank(1)  #  This should not be executed if called properly
	if( Blockdata["myid_on_node"] == 0 ):
		nz = tweight.get_zsize()
		ny = tweight.get_ysize()
		nx = tweight.get_xsize()
		tvol.set_attr("is_complex",1)
		if regularized:
			nr = len(cfsc)
			ovol = [0.0]*nr
			limitres = 0
			for i in range(nr):
				ovol[i] = min(max(cfsc[i], 0.0), 0.999)
				#print( i,cfsc[i] )
				if( ovol[i] == 0.0 ):
					limitres = i-1
					break
			if( limitres == 0 ): limitres = nr-2;
			ovol = reshape_1d(ovol, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			#print(" ovol  ", ovol)
			it = model_blank(2*nr)
			for i in range(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in range(5):  treg[i] = 0.0
			Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

	else:
		tvol = model_blank(1)
		tweight = model_blank(1)
		nz=0
		ny=0
		nx=0
		maxr2=0

	nx = bcast_number_to_all(nx, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	ny = bcast_number_to_all(ny, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	nz = bcast_number_to_all(nz, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	maxr2 = bcast_number_to_all(maxr2, source_node = 0, mpi_comm = Blockdata["shared_comm"])

	vol_data = get_image_data(tvol)
	we_data = get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"], n_iter)
	#Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])

	if( Blockdata["myid_on_node"] == 0 ):
		#  Either pad or window in F space to 2*nnxo
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"] ):
			tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"] ):
			tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)

		tvol = fft(tvol)
		tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.set_attr("npad",2)
		tvol.div_sinc(1)
		tvol.del_attr("npad")
		tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = cosinemask(tvol,Tracker["constants"]["nnxo"]//2-1-5,5, None) # clean artifacts in corners
		return tvol
	else:  return None

def calculate_2d_params_for_centering(kwargs):
	from sp_utilities 		import wrap_mpi_gatherv, read_text_row, write_text_row, bcast_number_to_all, get_im, combine_params2, model_circle, gather_compacted_EMData_to_root
	from sp_applications 	import MPI_start_end, ali2d_base 
	from sp_fundamentals 	import resample, rot_shift2D 
	from sp_filter 		import filt_ctf 
	from sp_global_def 	import ERROR
	
	
	#################################################################################################################################################################
	# get parameters from the dictionary
	init2dir = kwargs["init2dir"]
	myid = kwargs["myid"]
	main_node = kwargs["main_node"]
	number_of_images_in_stack = kwargs["number_of_images_in_stack"]
	nproc = kwargs["nproc"]

	target_radius = kwargs["target_radius"]
	# target_nx = kwargs["target_nx"]
	radi = kwargs["radi"]

	center_method = kwargs["center_method"]

	nxrsteps = kwargs["nxrsteps"]


	# stack_processed_by_ali2d_base__filename = kwargs["stack_processed_by_ali2d_base__filename"]
	command_line_provided_stack_filename = kwargs["command_line_provided_stack_filename"]

	# masterdir = kwargs["masterdir"]

	options_skip_prealignment = kwargs["options_skip_prealignment"]
	options_CTF = kwargs["options_CTF"]

	mpi_comm = kwargs["mpi_comm"]
	#################################################################################################################################################################

	if options_skip_prealignment:
		if(Blockdata["myid"] == 0):
			sxprint("=========================================")
			sxprint(" >>> There is no pre-alignment step.")
			sxprint("=========================================")
			return [[0, 0, 0, 0, 0] for i in range(number_of_images_in_stack)]
		else:  return [0.0]

	Finished_initial_2d_alignment = 1
	if(Blockdata["myid"] == Blockdata["main_node"]): 
		if( os.path.exists(os.path.join(init2dir, "Finished_initial_2d_alignment.txt")) ): Finished_initial_2d_alignment = 0
	Finished_initial_2d_alignment = bcast_number_to_all(Finished_initial_2d_alignment, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Finished_initial_2d_alignment == 1 ):
		if kwargs["radi"] < 1:
			ERROR("Particle radius has to be provided!", myid=Blockdata["myid"])

		if(Blockdata["myid"] == 0):
			import subprocess
			from sp_logger import Logger, BaseLogger_Files
			#  Create output directory
			log2d = Logger(BaseLogger_Files())
			log2d.prefix = os.path.join(init2dir)
			cmd = "mkdir "+log2d.prefix
			outcome = subprocess.call(cmd, shell=True)
			log2d.prefix += "/"
			# outcome = subprocess.call("sxheader.py  "+command_line_provided_stack_filename+"   --params=xform.align2d  --zero", shell=True)
		else:
			outcome = 0
			log2d = None

		if(Blockdata["myid"] == Blockdata["main_node"]):
			a = get_im(command_line_provided_stack_filename)
			nnxo = a.get_xsize()
		else:
			nnxo = 0
		nnxo = bcast_number_to_all(nnxo, source_node = Blockdata["main_node"])

		image_start, image_end = MPI_start_end(number_of_images_in_stack, Blockdata["nproc"], Blockdata["myid"])

		original_images = EMData.read_images(command_line_provided_stack_filename, list(range(image_start,image_end)))
		#  We assume the target radius will be 29, and xr = 1.  
		shrink_ratio = float(target_radius)/float(radi)

		for im in range(len(original_images)):
			if(shrink_ratio != 1.0):
				original_images[im]  = resample(original_images[im], shrink_ratio)

		nx = original_images[0].get_xsize()
		# nx = int(nx*shrink_ratio + 0.5)

		txrm = (nx - 2*(target_radius+1))//2
		if(txrm < 0):  			
			ERROR( "Radius of the structure larger than the window data size permits %d"%(radi), myid=Blockdata["myid"] )
		if(txrm/nxrsteps>0):
			tss = ""
			txr = ""
			while(txrm/nxrsteps>0):
				tts=txrm/nxrsteps
				tss += "  %d"%tts
				txr += "  %d"%(tts*nxrsteps)
				txrm =txrm//2
		else:
			tss = "1"
			txr = "%d"%txrm

		# section ali2d_base

		params2d = ali2d_base(original_images, init2dir, None, 1, target_radius, 1, txr, txr, tss, \
							False, 90.0, center_method, 14, options_CTF, 1.0, False, \
							"ref_ali2d", "", log2d, Blockdata["nproc"], Blockdata["myid"], Blockdata["main_node"], mpi_comm, write_headers = False)

		del original_images

		for i in range(len(params2d)):
			alpha, sx, sy, mirror = combine_params2(0, params2d[i][1],params2d[i][2], 0, -params2d[i][0], 0, 0, 0)
			sx /= shrink_ratio
			sy /= shrink_ratio
			params2d[i][0] = 0.0
			params2d[i][1] = sx
			params2d[i][2] = sy
			params2d[i][3] = 0

		mpi_barrier(MPI_COMM_WORLD)


		params2d = wrap_mpi_gatherv(params2d, Blockdata["main_node"], MPI_COMM_WORLD)
		if( Blockdata["myid"] == Blockdata["main_node"] ):		
			write_text_row(params2d,os.path.join(init2dir, "initial2Dparams.txt"))
			return params2d
		else:  return [0.0]
	else:
		if (Blockdata["myid"] == Blockdata["main_node"]):
			params2d = read_text_row(os.path.join(init2dir, "initial2Dparams.txt"))
			sxprint("Skipping 2d alignment since it was already done!")
			return params2d
		else:  return [0.0]

def Numrinit_local(first_ring, last_ring, skip=1, mode="F"):
	"""This function calculates the necessary information for the 2D 
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an 
	   		FFT-friendly power of the 2.
			
	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
	MAXFFT = 32768
	from math import pi
	skip = 1
	if (mode == 'f' or mode == 'F'): dpi = 2*pi
	else:                            dpi = pi
	numr = []
	lcirc = 1
	for k in range(first_ring, last_ring+1, skip):
		numr.append(k)
		jp = int(dpi * 1.5*k)
		ip = 2**(log2(jp))  # two times oversample each ring
		ip=min(MAXFFT,ip)
		#if ip >16 : ip/=2  #  subsample rings
		#if (k+skip <= last_ring and jp > ip+ip//2): ip=min(MAXFFT,2*ip)
		#if (k+skip  > last_ring and jp > ip+ip//5): ip=min(MAXFFT,2*ip)

		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip

	return  numr

def ali3D_polar_ccc(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	from sp_projection 	import prgs,prgl
	from sp_fundamentals 	import fft
	from sp_utilities 		import wrap_mpi_gatherv
	from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	#from operator 		import itemgetter#, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of polar CCC")


	at = time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	#  Fine grids are shifted by half-fine_step
	nang = len(refang)
	npsi = int(360./Tracker["delta"])
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw)

	disp_unit = np.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

		bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
	crefim = Util.Polar2Dm(model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)

	#  BIG BUFFER (for n_coarse_ang polar arrays)
	size_of_one_image = crefim.get_xsize()
	lenbigbuf = n_coarse_ang
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi_win_allocate_shared( size*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi_win_shared_query(win_sm, MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of setup

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

			bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi_win_shared_query(win_volinit, MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else:
		volinit = volprep
	#  End of replaced volprep


	at = time()

	nang_start, nang_end = MPI_start_end(n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])

	for i in range(nang_start, nang_end, 1):  # This will take care of process on a node less than nang.  Some loops will not be executed
		temp = prgl(volprep,[ coarse_angles[i][0], coarse_angles[i][1],0.0, 0.0,0.0], 1, True)
		crefim = Util.Polar2Dm(temp, cnx, cnx, numr, mode)
		Util.Normalize_ring(crefim, numr, 0)
		Util.Frngs(crefim, numr)
		Util.Applyws(crefim, numr, wr)
		bigbuffer.insert_clip(crefim,(0,i) )

	mpi_barrier(Blockdata["shared_comm"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		sxprint( "  Reference projections generated : %10.1fmin"%((time()-at)/60.))
	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift) )

	#  Preprocess the data

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)
	"""
	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [None]*nima
	norm_per_particle = [None]*nima
	"""

	lxod1 = n_coarse_ang*len(coarse_shifts)*n_coarse_psi
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	Blockdata['symclass'].set_angles(refang)
	Blockdata['symclass'].build_kdtree()
	#refdirs = angles_to_normals(refang)

	#if( Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

	at = time()

	for im in range(nima):

		#phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]
		phi,theta,psi,sx,sy = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4]#  First ITER

		if preshift:
			sx, sy = reduce_shifts(sx, sy, original_data[im])
			dataimage  = cyclic_shift(original_data[im],sx,sy)
			#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
			oldparams[im][3] = sx
			oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  dataimage = original_data[im].copy()

		st = get_image_statistics(dataimage, mask2D, False)
		dataimage -= st[0]
		dataimage /= st[1]
		##if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if apply_mask:
			if Tracker["constants"]["hardmask"]:
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"])
			else:
				bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
				bckg.set_attr("is_complex",1)
				bckg.set_attr("is_fftpad",1)
				bckg = fft(filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
				#  Normalize bckg noise in real space, only region actually used.
				st = get_image_statistics(bckg, mask2D, False)
				bckg -= st[0]
				bckg /= st[1]
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		#else:
		#	#  if no bckgnoise, do simple masking instead
		#	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

		
		#  Apply varadj, FIRST ITERATION, THERE IS NONE
		#if not nonorm:
		#	Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

		###  FT
		dataimage = fft(dataimage)
		##sig = Util.rotavg_fourier( dataimage )
		##accumulatepw[im] = sig[len(sig)//2:]+[0.0]

		#  We have to make sure the shifts are within correct range, shrinkage or not
		#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		"""
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
		else:
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
		bckgnoise = bckgn.copy()
		for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]
		"""
		if Tracker["constants"]["CTF"] :
			ctf_params = dataimage.get_attr("ctf")
			ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
			ctfa = ctf_img_real(Tracker["nxinit"], ctf_params)
			ctfs = ctfa
		dataml = fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
		data = []
		for iq in coarse_shifts:
			xx = iq[0]*shrink
			yy = iq[1]*shrink
			dss = fshift(dataml, xx, yy)
			dss.set_attr("is_complex",0)
			data.append(dss)

		#  This will get it to real space
		#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
		dataimage = fpol(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

		if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (im%(max(1,nima/5)) == 0) and (im>0)):
			sxprint( "  Number of images :%7d   %5d  %5.1f"%(im,nima,float(im)/float(nima)*100.) + "%" +"   %10.1fmin"%((time()-at)/60.))

		if( im < min(nima, 1) and procid == 0):
			#   Search all and compare with direct to figure what keepfirst might be
			keepfirst = (n_coarse_ang *  n_coarse_psi)/10#keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, keepfirst))

			assert(len(xod2) == keepfirst)

			xod1 = np.ndarray((keepfirst),dtype='f4',order="C")

			for iln in range(keepfirst):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# DO NOT order by angular directions to save time on reprojections.
			pre_ipsiandiang = -1
			for iln in range(keepfirst):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					nrmref = sqrt(Util.innerproduct(temp, temp, None))
					#Util.mul_scalar(temp, 1.0/nrmref)

				xod1[iln] = Util.innerproduct(data[ishift], temp, None)/nrmref  # peak
				##xod2[iln] = hashparams

			z = np.max(xod1)
			lina = np.argwhere(xod1 == z)[0]
			#if( Blockdata["myid"] == 0 ):  
			#print("  STARTING ",Blockdata["myid"],z,lina)
			keepf = int(lina[0]) + 1 # if position of the highest is on lina[0] it means we have to keep lina[0] + 1 elements
			xod1 = xod1[lina]
			xod2 = xod2[lina]

			"""			
			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["constants"]["ccfpercentage"]):
					lit = j+1
					break
			"""
			lit = 1
			keepf = [keepf]
			keepf = wrap_mpi_gatherv(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			if( Blockdata["myid"] == 0 ):
				keepf.sort()
				keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
			else:  keepf = 0
			keepf = wrap_mpi_bcast(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			Tracker["keepfirst"] = int(keepf)
			###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

		else:
			#Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, Tracker["keepfirst"]))

			xod1 = np.ndarray((Tracker["keepfirst"]),dtype='f4',order="C")

			for iln in range(Tracker["keepfirst"]):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# order by angular directions to save time on reprojections.
			ipsiandiang = xod2/1000
			lina = np.argsort(ipsiandiang)
			xod2 = xod2[lina]  # order does not matter
			pre_ipsiandiang = -1
			for iln in range(Tracker["keepfirst"]):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					temp.set_attr("is_complex",0)
					nrmref = sqrt(Util.innerproduct(temp, temp, None))
				peak = Util.innerproduct(data[ishift], temp, None)/nrmref
				xod1[iln] = peak
				##xod2[iln] = hashparams

			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order
			lit = Tracker["keepfirst"]
			"""
			xod1 -= xod1[0]

			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			xod1 = xod1[lina]
			xod2 = xod2[lina]
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["constants"]["ccfpercentage"]):
					lit = j+1
					break
			"""

		del data

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

		firstdirections = [[0.0,0.0] for iln in range(lit)]
		firstshifts = [0]*lit
		for iln in range(lit):
			hashparams = int(xod2[iln])
			ishift = hashparams%1000
			ipsiandiang	= hashparams/1000
			#ipsi = ipsiandiang%100000
			iang = ipsiandiang/100000
			firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
			firstshifts[iln] = ishift
		###del xod2
		###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
		#mpi_barrier(MPI_COMM_WORLD)
		#mpi_finalize()
		#exit()

		# Find neighbors, ltabang contains positions on refang list, no psis
		#ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = Tracker['howmany'])
		ltabang = Blockdata['symclass'].find_k_nearest_neighbors(firstdirections, k=Tracker['howmany'])

		# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
		#   even though it is longer than lit.

		#  Prepare image for chi2.
		#  We have to repeat everything from get shrink data, including shifts
		#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
		#    there are 2 psis, and at most n_fine_shifts. which should be 4.
		#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
		#      different shifts.  If so, we have to remove duplicates from the entire set.
		lcod1 = lit*4*2*n_fine_shifts
		cod2 = []
		#lol = 0
		for i1 in range(lit):
			hashparams = int(xod2[i1])
			ipsiandiang	= hashparams/1000
			oldiang = ipsiandiang/100000
			ipsi = ipsiandiang%100000
			ishift = hashparams%1000
			tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
			for i2 in range(Tracker['howmany']):
				iang = ltabang[i1][i2]
				if iang < 0:
					continue
				for i3 in range(2):  # psi
					itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
					itpsi = (itpsi + i3)%npsi
					for i4 in range(len(tshifts)):
						cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])

		del xod1, xod2

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
		cod2 = list(set(cod2))
		cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
		cod1.sort()

		lit = len(cod1)

		cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

		cod1 = np.ndarray(lit,dtype='f4',order="C")
		#cod1.fill(np.finfo(dtype='f4').min)
		#cod3 = np.ndarray(lit,dtype='f4',order="C")
		#cod3.fill(0.0)  #  varadj

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.

		data = [None]*nshifts
		johi = 0
		iln = 0
		prevdir = -1
		while(iln<lit):
			hashparams = cod2[iln]
			ipsiandiang	= hashparams/1000
			if(ipsiandiang != prevdir):
				prevdir = ipsiandiang
				ipsi = ipsiandiang%100000
				iang = ipsiandiang/100000
				temp = prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
				temp.set_attr("is_complex",0)
				nrmref = sqrt(Util.innerproduct(temp, temp, None))
				johi += 1
			while( ipsiandiang == cod2[iln]/1000 ):
				hashparams = cod2[iln]
				ishift = hashparams%1000
				if( data[ishift] == None ):
					xx = shifts[ishift][0]*shrink
					yy = shifts[ishift][1]*shrink
					data[ishift] = fshift(dataml, xx, yy)
					data[ishift].set_attr("is_complex",0)

				cod1[iln] = Util.innerproduct(data[ishift], temp, None)/nrmref
				iln += 1
				if(iln == lit  ):  break

		del data
		del dataml

		lina = np.argsort(cod1)
		cod1 = cod1[lina[::-1]]  # This sorts in reverse order
		cod2 = cod2[lina[::-1]]  # This sorts in reverse order
		##cod3 = cod3[lina[::-1]]  # This sorts in reverse order
		##cod1 -= cod1[0]
		##lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
		cod1 = cod1[0]
		cod2 = cod2[0]
		# = cod3[lina]
		"""
		np.exp(cod1, out=cod1)
		cod1 /= np.sum(cod1)
		cumprob = 0.0
		for j in range(len(cod1)):
			cumprob += cod1[j]
			if(cumprob > Tracker["constants"]["ccfpercentage"]):
				lit = j+1
				break
		"""
		#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
		###norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]

		###for iln in range(lit):
		newpar[im][2] = [[int(cod2), float(cod1)]]

		del cod1, cod2, lina
		###mpi_barrier(MPI_COMM_WORLD)
		###mpi_finalize()
		###exit()
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished projection matching   %10.1fmin"%((time()-at)/60.))
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations
	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()

	"""
	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		norm_per_particle[kl] = sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr            += norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))
	"""
	at = time()
	mpi_barrier(Blockdata["shared_comm"])

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi_win_free(win_vol)

	mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi_barrier(MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ):
	#	#write_text_row([[newpar[0][2][j][0],newpar[0][2][j][1]] for j in range(len(newpar[0][2]))],os.path.join(Tracker["directory"], "polar%1d.txt"%procid))
	#	sxprint( "  Statistics finished : %10.1fmin"%((time()-at)/60.))
	return newpar, [1.0]*nima#norm_per_particle

def ali3D_primary_polar(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	from sp_projection 	import prgs,prgl
	from sp_fundamentals 	import fft
	from sp_utilities 		import wrap_mpi_gatherv
	from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	#from operator 		import itemgetter, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of primary polar CCC")


	at = time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	#  Fine grids are shifted by half-fine_step
	nang = len(refang)
	npsi = int(360./Tracker["delta"])
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	ny = Tracker["nxinit"]
	nyp2 = ny/2
	nxth = (Tracker["nxinit"]+2)/2
	indx = model_blank(nxth, Tracker["nxinit"], 1, -1)
	tfrac = model_blank(nxth, Tracker["nxinit"])
	tcount = model_blank(nxth)
	for iy in range(1, ny+1):
		jy=iy-1
		if(jy>nyp2): jy=jy-ny
		argy = float(jy*jy)
		for ix in range(1,nxth+1):
			jx=ix-1
			roff = (jx+(iy-1)*nxth)
			if(mask[ix-1,iy-1] > 0.0 ):
				rf = sqrt( argy + float(jx*jx) )
				ir = int(rf)
				#print  ix-1,iy-1,roff,mask[ix-1,iy-1],rf,ir

				if( ir < nxth-1):
					frac = rf - float(ir)
					qres = 1.0 - frac
					tfrac[ix-1,iy-1] = frac
					#ioff = 2*roff
					tcount[ir]   	+= qres
					tcount[ir+1] 	+= frac
					indx[ix-1,iy-1]	= ir

	disp_unit = np.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

		bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
	crefim = Util.Polar2Dm(model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)

	#  BIG BUFFER (for n_coarse_ang polar arrays)
	size_of_one_image = crefim.get_xsize()
	lenbigbuf = n_coarse_ang
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi_win_allocate_shared( size*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi_win_shared_query(win_sm, MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of setup

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

			bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi_win_shared_query(win_volinit, MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi_barrier(Blockdata["shared_comm"])
	else:
		volinit = volprep
	#  End of replaced volprep


	at = time()

	nang_start, nang_end = MPI_start_end(n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])

	for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
		temp = prgl(volprep,[ coarse_angles[i][0], coarse_angles[i][1],0.0, 0.0,0.0], 1, True)
		crefim = Util.Polar2Dm(temp, cnx, cnx, numr, mode)
		Util.Frngs(crefim, numr)
		Util.Applyws(crefim, numr, wr)
		bigbuffer.insert_clip(crefim,(0,i) )

	mpi_barrier(Blockdata["shared_comm"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		sxprint( "  Reference projections generated : %10.1fmin"%((time()-at)/60.))

	###mpi_finalize()
	###exit()
	###  <><><><><><><><>

	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift) )

	#  Preprocess the data

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)



	if Blockdata["bckgnoise"] :
		oneover = []
		nxb = Blockdata["bckgnoise"][0].get_xsize()
		nyb = len(Blockdata["bckgnoise"])
		for i in range(nyb):
			temp = [0.0]*nxb
			for k in range(nxb):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

		if( procid == 0 ):
			Blockdata["totprob"] = [0.0]*nyb
			Blockdata["newbckgnoise"] = model_blank(nxb,nyb)

	accumulatepw = [None]*nima
	norm_per_particle = [None]*nima


	lxod1 = n_coarse_ang*len(coarse_shifts)*n_coarse_psi
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	Blockdata['symclass'].set_angles(refang)
	Blockdata['symclass'].build_kdtree()
	#refdirs = angles_to_normals(refang)

	#if( Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

	at = time()

	for im in range(nima):
		particle_group = original_data[im].get_attr("particle_group")

		phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

		if preshift:
			sx, sy = reduce_shifts(sx, sy, original_data[im])
			dataimage  = cyclic_shift(original_data[im],sx,sy)
			#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
			oldparams[im][3] = sx
			oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  dataimage = original_data[im].copy()

		st = get_image_statistics(dataimage, mask2D, False)
		dataimage -= st[0]
		dataimage /= st[1]
		if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if apply_mask:
			if Tracker["constants"]["hardmask"]:
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"])
			else:
				bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
				bckg.set_attr("is_complex",1)
				bckg.set_attr("is_fftpad",1)
				bckg = fft(filt_table(bckg, oneover[particle_group]))
				#  Normalize bckg noise in real space, only region actually used.
				st = get_image_statistics(bckg, mask2D, False)
				bckg -= st[0]
				bckg /= st[1]
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		#else:
		#	#  if no bckgnoise, do simple masking instead
		#	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

		#  Apply varadj
		if not nonorm:
			Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

		###  FT
		dataimage = fft(dataimage)
		sig = Util.rotavg_fourier( dataimage )
		accumulatepw[im] = sig[len(sig)//2:]+[0.0]

		#  We have to make sure the shifts are within correct range, shrinkage or not
		#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][particle_group]
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
		else:
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
		bckgnoise = bckgn.copy()
		for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

		if Tracker["constants"]["CTF"] :
			ctf_params = dataimage.get_attr("ctf")
			ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
			ctfa = ctf_img_real(Tracker["nxinit"], ctf_params)
			ctfs = ctfa
		dataml = fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
		data = []
		for iq in coarse_shifts:
			xx = iq[0]*shrink
			yy = iq[1]*shrink
			dss = fshift(dataml, xx, yy)
			dss.set_attr("is_complex",0)
			data.append(dss)

		#  This will get it to real space
		#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
		#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
		dataimage = fpol(Util.mulnclreal(Util.mulnclreal(dataml, Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

		if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (im%(max(1,nima/5)) == 0) and (im>0)):
			sxprint( "  Number of images :%7d   %5d  %5.1f"%(im,nima,float(im)/float(nima)*100.) + "%" +"   %10.1fmin"%((time()-at)/60.))

		if( im < min(nima, 1) and procid == 0):
			#   Search all and compare with direct to figure what keepfirst might be
			keepfirst = (n_coarse_ang *  n_coarse_psi)/10#keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, keepfirst))

			assert(len(xod2) == keepfirst)

			xod1 = np.ndarray((keepfirst),dtype='f4',order="C")

			for iln in range(keepfirst):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# DO NOT order by angular directions to save time on reprojections.
			pre_ipsiandiang = -1
			for iln in range(keepfirst):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					temp.set_attr("is_complex",0)
				xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)  # peak
				##xod2[iln] = hashparams

			xod1 -= np.max(xod1)
			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			#if( Blockdata["myid"] == 0 ):  sxprint("  STARTING ",np.max(xod1),np.min(xod1),len(lina),lina[-1])
			lina = lina.reshape(lina.size)
			keepf = int(lina[-1]) + 1

			xod1 = xod1[lina]
			xod2 = xod2[lina]

			
			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["ccfpercentage"]):
					lit = j+1
					break
			#keepf = mpi_reduce(keepf, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
			#if( Blockdata["myid"] == 0 ):
			#	keepf = max(int(float(keepf)*0.9),1)
			keepf = [keepf]
			keepf = wrap_mpi_gatherv(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			if( Blockdata["myid"] == 0 ):
				keepf.sort()
				keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
			else:  keepf = 0
			keepf = wrap_mpi_bcast(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			Tracker["keepfirst"] = int(keepf)
			###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])
				
		else:
			#Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, Tracker["keepfirst"]))

			xod1 = np.ndarray((Tracker["keepfirst"]),dtype='f4',order="C")

			for iln in range(Tracker["keepfirst"]):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# order by angular directions to save time on reprojections.
			ipsiandiang = xod2/1000
			lina = np.argsort(ipsiandiang)
			xod2 = xod2[lina]  # order does not matter
			pre_ipsiandiang = -1
			for iln in range(Tracker["keepfirst"]):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					temp.set_attr("is_complex",0)
				peak = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
				xod1[iln] = peak
				##xod2[iln] = hashparams

			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order

			xod1 -= xod1[0]

			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			xod1 = xod1[lina]
			xod2 = xod2[lina]
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["ccfpercentage"]):
					lit = j+1
					break


		del data

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

		firstdirections = [[0.0,0.0] for iln in range(lit)]
		firstshifts = [0]*lit
		for iln in range(lit):
			hashparams = int(xod2[iln])
			ishift = hashparams%1000
			ipsiandiang	= hashparams/1000
			#ipsi = ipsiandiang%100000
			iang = ipsiandiang/100000
			firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
			firstshifts[iln] = ishift

		# Find neighbors, ltabang contains positions on refang list, no psis
		#ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = Tracker['howmany'])
		ltabang = Blockdata['symclass'].find_k_nearest_neighbors(firstdirections, k=Tracker['howmany'])

		# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
		#   even though it is longer than lit.

		#  Prepare image for chi2.
		#  We have to repeat everything from get shrink data, including shifts
		#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
		#    there are 2 psis, and at most n_fine_shifts. which should be 4.
		#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
		#      different shifts.  If so, we have to remove duplicates from the entire set.
		lcod1 = lit*4*2*n_fine_shifts
		cod2 = []
		for i1 in range(lit):
			hashparams = int(xod2[i1])
			ipsiandiang	= hashparams/1000
			oldiang = ipsiandiang/100000
			ipsi = ipsiandiang%100000
			ishift = hashparams%1000
			tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
			for i2 in range(Tracker['howmany']):
				iang = ltabang[i1][i2]
				if iang < 0:
					continue
				for i3 in range(2):  # psi
					itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
					itpsi = (itpsi + i3)%npsi
					for i4 in range(len(tshifts)):
						cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])

		del xod1, xod2

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
		cod2 = list(set(cod2))
		cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
		cod1.sort()

		lit = len(cod1)

		cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

		cod1 = np.ndarray(lit,dtype='f4',order="C")
		#cod1.fill(np.finfo(dtype='f4').min)
		cod3 = np.ndarray(lit,dtype='f4',order="C")
		#cod3.fill(0.0)  #  varadj

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
		#  Make sure volprep has nxinit size
		tbckg = []
		data = [None]*nshifts
		johi = 0
		iln = 0
		prevdir = -1
		while(iln<lit):
			hashparams = cod2[iln]
			ipsiandiang	= hashparams/1000
			if(ipsiandiang != prevdir):
				prevdir = ipsiandiang
				ipsi = ipsiandiang%100000
				iang = ipsiandiang/100000
				temp = prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
				temp.set_attr("is_complex",0)
				johi += 1
			while( ipsiandiang == cod2[iln]/1000 ):
				hashparams = cod2[iln]
				ishift = hashparams%1000
				if( data[ishift] == None ):
					xx = shifts[ishift][0]*shrink
					yy = shifts[ishift][1]*shrink
					data[ishift] = fshift(dataml, xx, yy)
					data[ishift].set_attr("is_complex",0)

				##[peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
				fofa = Util.sqednormbckg(data[ishift], temp, ctfa, bckgnoise, indx, tfrac, tcount)
				cod1[iln] = -fofa[-2]   # -peak
				cod3[iln] =  fofa[-1]   #  varadj
				tbckg.append(fofa[:-2])
				iln += 1
				if(iln == lit  ):  break

		del data
		del dataml

		lina = np.argsort(cod1)
		lina = lina[::-1]  # This sorts in reverse order
		cod1 = cod1[lina]
		cod2 = cod2[lina]
		cod3 = cod3[lina]
		cod1 -= cod1[0]
		tbckg = [tbckg[int(q)] for q in lina]
		lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
		cod1 = cod1[lina]
		cod2 = cod2[lina]
		cod3 = cod3[lina]
		tbckg = [tbckg[int(q)] for q in lina]

		np.exp(cod1, out=cod1)
		cod1 /= np.sum(cod1)
		cumprob = 0.0
		for j in range(len(cod1)):
			cumprob += cod1[j]
			if(cumprob > Tracker["ccfpercentage"]):
				lit = j+1
				break
		
		#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
		norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
		atbckg = [0.0]*len(tbckg[0])
		for iln in range(lit):
			prob = float(cod1[iln])
			Blockdata["totprob"][particle_group] += prob
			for iq in range(len(tbckg[0])):
				atbckg[iq] += tbckg[iln][iq]*prob

		del tbckg
		for iq in range(nxth):  Blockdata["newbckgnoise"][iq,particle_group] += atbckg[iq]
		del atbckg
		for iln in range(lit):
			newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])

		del cod1, cod2, cod3, lina
		###mpi_barrier(MPI_COMM_WORLD)
		###mpi_finalize()
		###exit()
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished projection matching   %10.1fmin"%((time()-at)/60.))
	at = time()
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations
	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		norm_per_particle[kl] = sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr            += norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))

	at = time()
	mpi_barrier(MPI_COMM_WORLD)
	mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi_win_free(win_vol)

	mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	# Compute new background noise
	# Reduce stuff
	if( procid == 1 ):
		Blockdata["totprob"] = mpi_reduce(Blockdata["totprob"], nyb, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		reduce_EMData_to_root(Blockdata["newbckgnoise"], Blockdata["myid"], Blockdata["main_node"])
		if( Blockdata["myid"] == 0 ):
			for igrp in range(nyb):
				Blockdata["newbckgnoise"][0, igrp] = 1.0
				for i in range(1,nxth):
					if(Blockdata["newbckgnoise"][i, igrp] > 0.0):  Blockdata["newbckgnoise"][i, igrp] = 2.0*Blockdata["totprob"][igrp]/Blockdata["newbckgnoise"][i, igrp]  # normalize and invert
				for i in range(nxth,nxb):
					Blockdata["newbckgnoise"][i, igrp] = Blockdata["bckgnoise"][igrp][i]
			Blockdata["newbckgnoise"].write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf")) #  Write updated bckgnoise to current directory

		bcast_EMData_to_all(Blockdata["newbckgnoise"], Blockdata["myid"], source_node = Blockdata["main_node"], comm = MPI_COMM_WORLD)
		for igrp in range(nyb):
			for i in range(nxb):
				Blockdata["bckgnoise"][igrp][i] = Blockdata["newbckgnoise"][i, igrp]
		del Blockdata["newbckgnoise"]
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished sigma2   %10.1fmin"%((time()-at)/60.))
	return newpar, norm_per_particle

def ali3D_polar(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	from sp_projection 	import prgs,prgl
	from sp_fundamentals 	import fft
	from sp_utilities 		import wrap_mpi_gatherv
	from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	from operator 		import itemgetter#, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered exhaustive CCC")


	at = time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	#  Fine grids are shifted by half-fine_step
	nang = len(refang)
	npsi = int(360./Tracker["delta"])
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw)

	disp_unit = np.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

		bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})

	crefim = Util.Polar2Dm(model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)

	#  BIG BUFFER (for n_coarse_ang polar arrays)
	size_of_one_image = crefim.get_xsize()
	lenbigbuf = n_coarse_ang
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi_win_allocate_shared( size*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi_win_shared_query(win_sm, MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of setup

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

			bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi_win_shared_query(win_volinit, MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else:
		volinit = volprep
	#  End of replaced volprep


	at = time()

	nang_start, nang_end = MPI_start_end(n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])

	for i in range(nang_start, nang_end, 1):  # This will take care of process on a node less than nang.  Some loops will not be executed
		temp = prgl(volprep,[ coarse_angles[i][0], coarse_angles[i][1],0.0, 0.0,0.0], 1, True)
		crefim = Util.Polar2Dm(temp, cnx, cnx, numr, mode)
		Util.Frngs(crefim, numr)
		Util.Applyws(crefim, numr, wr)
		bigbuffer.insert_clip(crefim,(0,i) )

	mpi_barrier(Blockdata["shared_comm"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		sxprint( "  Reference projections generated : %10.1fmin"%((time()-at)/60.))
	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift) )

	#  Preprocess the data

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)



	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [None]*nima
	norm_per_particle = [None]*nima


	lxod1 = n_coarse_ang*len(coarse_shifts)*n_coarse_psi
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	Blockdata['symclass'].set_angles(refang)
	Blockdata['symclass'].build_kdtree()
	#refdirs = angles_to_normals(refang)

	#if( Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

	###firsti = True

	at = time()

	for im in range(nima):

		phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

		if preshift:
			sx, sy = reduce_shifts(sx, sy, original_data[im])
			dataimage  = cyclic_shift(original_data[im],sx,sy)
			#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
			oldparams[im][3] = sx
			oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  dataimage = original_data[im].copy()

		st = get_image_statistics(dataimage, mask2D, False)
		dataimage -= st[0]
		dataimage /= st[1]
		if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		#if Blockdata["bckgnoise"]:  # I think we should assume it always exists
		if apply_mask:
			if Tracker["constants"]["hardmask"]:
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"])
			else:
				bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
				bckg.set_attr("is_complex",1)
				bckg.set_attr("is_fftpad",1)
				bckg = fft(filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
				#  Normalize bckg noise in real space, only region actually used.
				st = get_image_statistics(bckg, mask2D, False)
				bckg -= st[0]
				bckg /= st[1]
				dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		#else:
		#	#  if no bckgnoise, do simple masking instead
		#	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

		#  Apply varadj
		if not nonorm:
			Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

		###  FT
		dataimage = fft(dataimage)
		sig = Util.rotavg_fourier( dataimage )
		accumulatepw[im] = sig[len(sig)//2:]+[0.0]

		#  We have to make sure the shifts are within correct range, shrinkage or not
		#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
		else:
			bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
		bckgnoise = bckgn.copy()
		for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

		if Tracker["constants"]["CTF"] :
			ctf_params = dataimage.get_attr("ctf")
			ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
			ctfa = ctf_img_real(Tracker["nxinit"], ctf_params)
			ctfs = ctfa
		dataml = fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
		data = []
		for iq in coarse_shifts:
			xx = iq[0]*shrink
			yy = iq[1]*shrink
			dss = fshift(dataml, xx, yy)
			dss.set_attr("is_complex",0)
			data.append(dss)

		#  This will get it to real space
		#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
		#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
		dataimage = fpol(Util.mulnclreal(Util.mulnclreal(dataml, Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

		if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (im%(max(1,nima/5)) == 0) and (im>0)):
			sxprint( "  Number of images :%7d   %5d  %5.1f"%(im,nima,float(im)/float(nima)*100.) + "%" +"   %10.1fmin"%((time()-at)/60.))

		if( im < min(nima, 1) and procid == 0):
			#   Search all and compare with direct to figure what keepfirst might be
			keepfirst = (n_coarse_ang *  n_coarse_psi)/10#keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, keepfirst))

			assert(len(xod2) == keepfirst)

			xod1 = np.ndarray((keepfirst),dtype='f4',order="C")

			for iln in range(keepfirst):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# DO NOT order by angular directions to save time on reprojections.
			pre_ipsiandiang = -1
			for iln in range(keepfirst):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					temp.set_attr("is_complex",0)
				xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)  # peak
				##xod2[iln] = hashparams

			xod1 -= np.max(xod1)
			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			#if( Blockdata["myid"] == 0 ):  
			#print("  STARTING ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])
			lina = lina.reshape(lina.size)
			keepf = int(lina[-1]) + 1

			xod1 = xod1[lina]
			xod2 = xod2[lina]

			
			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["ccfpercentage"]):
					lit = j+1
					break
			#keepf = mpi_reduce(keepf, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
			#if( Blockdata["myid"] == 0 ):
			#	keepf = max(int(float(keepf)*0.9),1)
			keepf = [keepf]
			keepf = wrap_mpi_gatherv(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			if( Blockdata["myid"] == 0 ):
				keepf.sort()
				keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
			else:  keepf = 0
			keepf = wrap_mpi_bcast(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
			Tracker["keepfirst"] = int(keepf)
			###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

		else:
			#Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

			xod2 = np.asarray(Util.multiref_Crosrng_msg_stack_stepsi(dataimage, bigbuffer, \
					coarse_shifts_shrank,\
					numr, [coarse_angles[ic][2] for ic in range(n_coarse_ang)], coarse_delta, cnx, Tracker["keepfirst"]))

			xod1 = np.ndarray((Tracker["keepfirst"]),dtype='f4',order="C")

			for iln in range(Tracker["keepfirst"]):
				m = xod2[iln]
				j = m%n_coarse_psi
				ic = (m/n_coarse_psi)%n_coarse_ang
				ib  = m/(n_coarse_ang*n_coarse_psi)
				xod2[iln] = j*1000 + ic*100000000 + ib #hashparams
			# order by angular directions to save time on reprojections.
			ipsiandiang = xod2/1000
			lina = np.argsort(ipsiandiang)
			xod2 = xod2[lina]  # order does not matter
			pre_ipsiandiang = -1
			for iln in range(Tracker["keepfirst"]):
				hashparams	= int(xod2[iln])
				ishift		= hashparams%1000
				ipsiandiang	= hashparams/1000
				if(ipsiandiang != pre_ipsiandiang):
					pre_ipsiandiang = ipsiandiang
					ipsi = ipsiandiang%100000
					iang = ipsiandiang/100000
					temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
					temp.set_attr("is_complex",0)
				peak = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
				xod1[iln] = peak
				##xod2[iln] = hashparams

			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order

			xod1 -= xod1[0]

			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			xod1 = xod1[lina]
			xod2 = xod2[lina]
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["ccfpercentage"]):
					lit = j+1
					break


		del data

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

		firstdirections = [[0.0,0.0] for iln in range(lit)]
		firstshifts = [0]*lit
		for iln in range(lit):
			hashparams = int(xod2[iln])
			ishift = hashparams%1000
			ipsiandiang	= hashparams/1000
			#ipsi = ipsiandiang%100000
			iang = ipsiandiang/100000
			firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
			firstshifts[iln] = ishift
		###del xod2
		###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
		#mpi_barrier(MPI_COMM_WORLD)
		#mpi_finalize()
		#exit()

		# Find neighbors, ltabang contains positions on refang list, no psis
		#ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = Tracker['howmany'])
		ltabang = Blockdata['symclass'].find_k_nearest_neighbors(firstdirections, k=Tracker['howmany'])

		# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
		#   even though it is longer than lit.

		#  Prepare image for chi2.
		#  We have to repeat everything from get shrink data, including shifts
		#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
		#    there are 2 psis, and at most n_fine_shifts. which should be 4.
		#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
		#      different shifts.  If so, we have to remove duplicates from the entire set.
		lcod1 = lit*4*2*n_fine_shifts
		cod2 = []
		#lol = 0
		for i1 in range(lit):
			hashparams = int(xod2[i1])
			ipsiandiang	= hashparams/1000
			oldiang = ipsiandiang/100000
			ipsi = ipsiandiang%100000
			ishift = hashparams%1000
			tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
			for i2 in range(Tracker['howmany']):
				iang = ltabang[i1][i2]
				if iang < 0:
					continue
				for i3 in range(2):  # psi
					itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
					itpsi = (itpsi + i3)%npsi
					for i4 in range(len(tshifts)):
						cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])


		del xod1, xod2

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
		cod2 = list(set(cod2))
		cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
		cod1.sort()

		lit = len(cod1)

		cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

		cod1 = np.ndarray(lit,dtype='f4',order="C")
		#cod1.fill(np.finfo(dtype='f4').min)
		cod3 = np.ndarray(lit,dtype='f4',order="C")
		#cod3.fill(0.0)  #  varadj

		#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.

		data = [None]*nshifts
		johi = 0
		iln = 0
		prevdir = -1
		while(iln<lit):
			hashparams = cod2[iln]
			ipsiandiang	= hashparams/1000
			if(ipsiandiang != prevdir):
				prevdir = ipsiandiang
				ipsi = ipsiandiang%100000
				iang = ipsiandiang/100000
				temp = prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
				temp.set_attr("is_complex",0)
				johi += 1
			while( ipsiandiang == cod2[iln]/1000 ):
				hashparams = cod2[iln]
				ishift = hashparams%1000
				if( data[ishift] == None ):
					xx = shifts[ishift][0]*shrink
					yy = shifts[ishift][1]*shrink
					data[ishift] = fshift(dataml, xx, yy)
					data[ishift].set_attr("is_complex",0)

				[peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
				cod1[iln] = -peak
				cod3[iln] = varadj
				iln += 1
				if(iln == lit  ):  break

		del data
		del dataml

		lina = np.argsort(cod1)
		cod1 = cod1[lina[::-1]]  # This sorts in reverse order
		cod2 = cod2[lina[::-1]]  # This sorts in reverse order
		cod3 = cod3[lina[::-1]]  # This sorts in reverse order
		cod1 -= cod1[0]
		lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
		cod1 = cod1[lina]
		cod2 = cod2[lina]
		cod3 = cod3[lina]

		np.exp(cod1, out=cod1)
		cod1 /= np.sum(cod1)
		cumprob = 0.0
		for j in range(len(cod1)):
			cumprob += cod1[j]
			if(cumprob > Tracker["ccfpercentage"]):
				lit = j+1
				break

		#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
		norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]

		for iln in range(lit):
			newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])

		del cod1, cod2, cod3, lina
		###mpi_barrier(MPI_COMM_WORLD)
		###mpi_finalize()
		###exit()
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished projection matching   %10.1fmin"%((time()-at)/60.))
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations
	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		norm_per_particle[kl] = sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr            += norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))

	at = time()
	mpi_barrier(Blockdata["shared_comm"])

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi_win_free(win_vol)

	mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		#write_text_row([[newpar[0][2][j][0],newpar[0][2][j][1]] for j in range(len(newpar[0][2]))],os.path.join(Tracker["directory"], "polar%1d.txt"%procid))
		sxprint( "  Statistics finished : %10.1fmin"%((time()-at)/60.))
	return newpar, norm_per_particle

# MODIFIED FROM TRUE POLAR  12/06/2017  PAP
def ali3D_primary_local_polar(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	"""
	12/06/2017
	"""
	global Tracker, Blockdata
	from sp_projection 	import prgs,prgl
	from sp_fundamentals 	import fft
	from sp_utilities 		import wrap_mpi_gatherv
	from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	#from operator 		import itemgetter, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered primary local polar")

	at = time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = cos(radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
	m_coarse_psi = int(Tracker["an"]/Tracker["delta"] + 0.5)
	c_coarse_psi = m_coarse_psi//2
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	ny = Tracker["nxinit"]
	nyp2 = ny/2
	nxth = (Tracker["nxinit"]+2)/2
	indx = model_blank(nxth, Tracker["nxinit"], 1, -1)
	tfrac = model_blank(nxth, Tracker["nxinit"])
	tcount = model_blank(nxth)
	for iy in range(1, ny+1):
		jy=iy-1
		if(jy>nyp2): jy=jy-ny
		argy = float(jy*jy)
		for ix in range(1,nxth+1):
			jx=ix-1
			roff = (jx+(iy-1)*nxth)
			if(mask[ix-1,iy-1] > 0.0 ):
				rf = sqrt( argy + float(jx*jx) )
				ir = int(rf)
				#print  ix-1,iy-1,roff,mask[ix-1,iy-1],rf,ir

				if( ir < nxth-1):
					frac = rf - float(ir)
					qres = 1.0 - frac
					tfrac[ix-1,iy-1] = frac
					#ioff = 2*roff
					tcount[ir]   	+= qres
					tcount[ir+1] 	+= frac
					indx[ix-1,iy-1]	= ir

	disp_unit = np.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

		bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

			bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi_win_shared_query(win_volinit, MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi_barrier(Blockdata["shared_comm"])
	else:
		volinit = volprep
	#  End of replaced volprep


	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = Util.Polar2Dm(model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	#Blockdata['symclass_coarse'].set_angles(coarse_angles)
	Blockdata['symclass_coarse'].set_angles(coarse_angles)
	Blockdata['symclass_coarse'].build_kdtree()
	#normals_set = angles_to_normals(coarse_angles)
	Blockdata["angle_set"] = coarse_angles
	if( n_coarse_ang <= numberofrefs_inmem ):
		number_of_cones = 1
		numberofrefs_inmem = n_coarse_ang
		assignments_to_cones = [list(range(len(oldparams)))]

		assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
		assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid

		for i,q in enumerate(assignments_to_cones):
			#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
			###print( "in loop ", Blockdata["myid"],i,len(q))#,q

			for m in q:
				#print " m ",m,len(angles)

				assignments_of_refangles_to_angles[m] = [entry for entry in Blockdata['symclass_coarse'].find_nearest_neighbors(oldparams[m], Tracker["an"]) if not entry < 0]
				assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

			assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
			#print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
			assignments_of_refangles_to_cones[i] = wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
			if(Blockdata["symclass"].sym[0] == "c"):
				if( number_of_cones == 1 ):
					cone_delta  = 180.1
					cone_angles = [[0., 1.0, 0.]]
				else:
					cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0, theta2=89.0)
					cone_angles += [[(q[0]+90.0)%360., 180.-q[1],0] for q in cone_angles]
					cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
			else:
				cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

			#if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = [entry for entry in Blockdata['symclass_coarse'].find_nearest_neighbors(oldparams[m], Tracker["an"]) if not entry < 0]
						#if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = bcast_number_to_all(doit, source_node = 0)
				number_of_cones = bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  sxprint( " number_of_cones: ",number_of_cones)
	#  Maximum number of refangles assigned to angles (max number of references per image)
	nlocal_angles = max( [ len(q) for q in assignments_of_refangles_to_angles] )
	#  Most likely we have to delete some lists before proceeding
	# We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
	#  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
	numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

	###for i,q in enumerate(assignments_of_refangles_to_cones):
	###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi_win_allocate_shared( size*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi_win_shared_query(win_sm, MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nxb = Blockdata["bckgnoise"][0].get_xsize()
		nyb = len(Blockdata["bckgnoise"])
		for i in range(nyb):
			temp = [0.0]*nxb
			for k in range(nxb):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

		if( procid == 0 ):
			Blockdata["totprob"] = [0.0]*nyb
			Blockdata["newbckgnoise"] = model_blank(nxb,nyb)

	accumulatepw = [0.0]*nima
	norm_per_particle = [0.0]*nima

	##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	Blockdata['symclass'].set_angles(refang)
	Blockdata['symclass'].build_kdtree()
	#refdirs = angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi_reduce(max_number_of_cones, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi_bcast(max_number_of_cones, 1, MPI_INT, 0, MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])


	#if( Blockdata["myid_on_node"] == 0):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	for icone in range(max_number_of_cones):
		mpi_barrier(MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				temp = prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				crefim = Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				Util.Frngs(crefim, numr)
				Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
			###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
			###mpi_finalize()
			###exit()
			###  <><><><><><><><>

			#  Preprocess the data

			#  We only process images in the current cone.  icnm is consecutive number, im the actual image number
			#for icnm,im in enumerate(assignments_to_cones[icone]):
			lenass = len(assignments_to_cones[icone])
			###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
			for icnm in range(max(1,lenass)):  # I have to enter the loop even it there is no assignment
				if( lenass == 0 ):
					keepf = -1
					###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
				else:
					#  PREPARE ONE IMAGE
					im = assignments_to_cones[icone][icnm]
					particle_group = original_data[im].get_attr("particle_group")

					phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx, sy = reduce_shifts(sx, sy, original_data[im])
						dataimage  = cyclic_shift(original_data[im],sx,sy)
						#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
						oldparams[im][3] = sx
						oldparams[im][4] = sy
						sx = 0.0
						sy = 0.0
					else:  dataimage = original_data[im].copy()


					st = get_image_statistics(dataimage, mask2D, False)
					dataimage -= st[0]
					dataimage /= st[1]
					if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
					#  Do bckgnoise if exists
					if Blockdata["bckgnoise"]:
						if apply_mask:
							if Tracker["constants"]["hardmask"]:
								dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = fft(filt_table(bckg, oneover[particle_group]))
								#  Normalize bckg noise in real space, only region actually used.
								st = get_image_statistics(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = fft(dataimage)
					sig = Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][particle_group]
						bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = fpol(Util.mulnclreal(Util.mulnclreal(dataml, Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					sxprint( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1



				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						lxod1 = Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''
						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# DO NOT order by angular directions to save time on reprojections.

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
							##junk = time()
							xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= np.max(xod1)
						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						#if( Blockdata["myid"] == 0 ):  
						###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])


						lina = lina.reshape(lina.size)
						keepf = int(lina[-1]) + 1

						xod1 = xod1[lina]
						xod2 = xod2[lina]

						###print("  STARTING2    ",Blockdata["myid"])

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						###print("  STARTING3    ",Blockdata["myid"],lit)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["ccfpercentage"]):
								lit = j+1
								break


						###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
					#  Turn into percentage of all possible
					keepf = [int(float(keepf*100)/float(keepfirst))]
					###mpi_barrier(MPI_COMM_WORLD)
					###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
					mpi_barrier(MPI_COMM_WORLD)
					keepf = wrap_mpi_gatherv(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
					###print("  STARTING6    ",Blockdata["myid"],keepf)
					if( Blockdata["myid"] == 0 ):
						keepf = [junk for junk in keepf if junk >0]
						if( len(keepf) <2 ):
							keepf = 3
						else:
							keepf.sort()
							keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
					else:  keepf = 0
					###print("  STARTING7    ",Blockdata["myid"],keepf)
					keepf = wrap_mpi_bcast(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
					#if(keepf == 0):
					#	ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
					#	mpi_finalize()
					#	exit()
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						keepfirst = lang*m_coarse_psi*n_coarse_shifts
						keepfirst = max(keepfirst*Tracker["keepfirst"]/100,min(keepfirst,3))
						lxod1 = Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
						##'''
						assert( keepfirst == len(lxod1)/3 )
						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = np.argsort(ipsiandiang)
						xod2 = xod2[lina]  # order does not matter

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
						#if( lima<5):  sxprint("  SECOND KEPT  ",lit)


				#mpi_barrier(MPI_COMM_WORLD)
				#mpi_finalize()
				#exit()
				#for j in range(lit):
				#	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
				if( lenass > 0):
					###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

					firstdirections = [[0.0,0.0] for iln in range(lit)]
					firstshifts = [0]*lit
					for iln in range(lit):
						hashparams = int(xod2[iln])
						ishift = hashparams%1000
						ipsiandiang	= hashparams/1000
						#ipsi = ipsiandiang%100000
						iang = ipsiandiang/100000
						#try:
						firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
						#except:
						#	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					#ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = Tracker['howmany'])
					ltabang = Blockdata['symclass'].find_k_nearest_neighbors(firstdirections, k=Tracker['howmany'])
					###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

					#  Prepare image for chi2.
					#  We have to repeat everything from get shrink data, including shifts
					#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
					#    there are 3 psis, and at most n_fine_shifts. 
					#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
					#      different shifts.  If so, we have to remove duplicates from the entire set.
					lcod1 = lit*4*2*n_fine_shifts

					#cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
					#cod2.fill(-1)  #  hashparams
					cod2 = []
					#lol = 0
					for i1 in range(lit):
						hashparams = int(xod2[i1])
						ipsiandiang	= hashparams/1000
						oldiang = ipsiandiang/100000
						ipsi = ipsiandiang%100000
						ishift = hashparams%1000
						tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
						for i2 in range(Tracker['howmany']):
							iang = ltabang[i1][i2]
							if iang < 0:
								continue
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = np.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = np.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					tbckg = []
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = fshift(dataml, xx, yy)
								data[ishift].set_attr("is_complex",0)
							##junk = time()
							#[peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
							fofa = Util.sqednormbckg(data[ishift], temp, ctfa, bckgnoise, indx, tfrac, tcount)
							cod1[iln] = -fofa[-2]   # -peak
							cod3[iln] =  fofa[-1]   #  varadj
							tbckg.append(fofa[:-2])
							iln += 1
							if(iln == lit  ):  break
						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	temp.write_image("temp.hdf")
						#	data[iln].write_image("data.hdf")
						#	ctfa.write_image("ctfa.hdf")
						#	bckgnoise.write_image("bckgnoise.hdf")
						#	exit()
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = np.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					tbckg = [tbckg[int(q)] for q in lina]
					lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]
					tbckg = [tbckg[int(q)] for q in lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					np.exp(cod1, out=cod1)
					cod1 /= np.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
					atbckg = [0.0]*len(tbckg[0])
					for iln in range(lit):
						prob = float(cod1[iln])
						Blockdata["totprob"][particle_group] += prob
						for iq in range(len(tbckg[0])):
							atbckg[iq] += tbckg[iln][iq]*prob

					del tbckg
					for iq in range(nxth):  Blockdata["newbckgnoise"][iq,particle_group] += atbckg[iq]
					del atbckg

					for iln in range(lit):
						newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
						#	hashparams = newpar[im][2][iln][0]
						#	ipsiandiang	= hashparams/1000
						#	ipsi = ipsiandiang%100000
						#	iang = ipsiandiang/100000
						#	ishift = hashparams%1000
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
						#	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""

	#  END OF CONES
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished projection matching   %10.1fmin"%((time()-at)/60.))
	at = time()
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations

	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
		norm_per_particle[kl]	= sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr				+= norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))

	at = time()
	mpi_barrier(MPI_COMM_WORLD)

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi_win_free(win_vol)

	mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	mpi_barrier(Blockdata["shared_comm"])

	# Compute new background noise
	# Reduce stuff
	if( procid == 1 ):
		Blockdata["totprob"] = mpi_reduce(Blockdata["totprob"], nyb, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		reduce_EMData_to_root(Blockdata["newbckgnoise"], Blockdata["myid"], Blockdata["main_node"])
		if( Blockdata["myid"] == 0 ):
			for igrp in range(nyb):
				Blockdata["newbckgnoise"][0, igrp] = 1.0
				for i in range(1,nxth):
					if(Blockdata["newbckgnoise"][i, igrp] > 0.0):  Blockdata["newbckgnoise"][i, igrp] = 2.0*Blockdata["totprob"][igrp]/Blockdata["newbckgnoise"][i, igrp]  # normalize and invert
				for i in range(nxth,nxb):
					Blockdata["newbckgnoise"][i, igrp] = Blockdata["bckgnoise"][igrp][i]
			Blockdata["newbckgnoise"].write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf")) #  Write updated bckgnoise to current directory

		bcast_EMData_to_all(Blockdata["newbckgnoise"], Blockdata["myid"], source_node = Blockdata["main_node"], comm = MPI_COMM_WORLD)
		for igrp in range(nyb):
			for i in range(nxb):
				Blockdata["bckgnoise"][igrp][i] = Blockdata["newbckgnoise"][i, igrp]
		del Blockdata["newbckgnoise"]

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished sigma2   %10.1fmin"%((time()-at)/60.))
	return newpar, norm_per_particle

# TRUE POLAR
def ali3D_local_polar(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	"""
	02/07/2017
	"""
	global Tracker, Blockdata
	from sp_projection 	import prgs,prgl
	from sp_fundamentals 	import fft
	from sp_utilities 		import wrap_mpi_gatherv
	from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	#from operator 		import itemgetter#, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from sp_fundamentals import fdecimate, fshift, fft
	from sp_filter       import filt_ctf, filt_table
	from sp_applications import MPI_start_end
	from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered local polar")

	at = time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = cos(radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
	m_coarse_psi = int(Tracker["an"]/Tracker["delta"] + 0.5)
	c_coarse_psi = m_coarse_psi//2
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

	#if(Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint( original_data[0].get_attr("identifier") )
	#	sxprint( original_data[1].get_attr("identifier") )

	disp_unit = np.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

		bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = model_blank( nxvol,nyvol, nzvol)

			bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi_win_shared_query(win_volinit, MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else:
		volinit = volprep
	#  End of replaced volprep


	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = Util.Polar2Dm(model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	Blockdata['symclass_coarse'].set_angles(coarse_angles)
	Blockdata['symclass_coarse'].build_kdtree()
	#normals_set = angles_to_normals(coarse_angles)
	Blockdata["angle_set"] = coarse_angles
	if( n_coarse_ang <= numberofrefs_inmem ):
		number_of_cones = 1
		numberofrefs_inmem = n_coarse_ang
		assignments_to_cones = [list(range(len(oldparams)))]

		assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
		assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid

		for i,q in enumerate(assignments_to_cones):
			#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
			###print( "in loop ", Blockdata["myid"],i,len(q))#,q

			for m in q:
				#print " m ",m,len(angles)

				assignments_of_refangles_to_angles[m] = [entry for entry in Blockdata['symclass_coarse'].find_nearest_neighbors(oldparams[m], Tracker["an"]) if not entry < 0]
				assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

			assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
			#print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
			assignments_of_refangles_to_cones[i] = wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
			if(Blockdata["symclass"].sym[0] == "c"):
				if( number_of_cones == 1 ):
					cone_delta  = 180.1
					cone_angles = [[0., 1.0, 0.]]
				else:
					cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0, theta2=89.0)
					cone_angles += [[(q[0]+90.0)%360., 180.-q[1],0] for q in cone_angles]
					cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
			else:
				cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

			#if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = [entry for entry in Blockdata['symclass_coarse'].find_nearest_neighbors(oldparams[m], Tracker["an"]) if not entry < 0]
						#if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = bcast_number_to_all(doit, source_node = 0)
				number_of_cones = bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  sxprint( " number_of_cones: ",number_of_cones)
	#  Maximum number of refangles assigned to angles (max number of references per image)
	nlocal_angles = max( [ len(q) for q in assignments_of_refangles_to_angles] )
	#  Most likely we have to delete some lists before proceeding
	# We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
	#  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
	numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

	###for i,q in enumerate(assignments_of_refangles_to_cones):
	###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi_win_allocate_shared( size*disp_unit , disp_unit, MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi_win_shared_query(win_sm, MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  " )
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [0.0]*nima
	norm_per_particle = [0.0]*nima

	##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	Blockdata['symclass'].set_angles(refang)
	Blockdata['symclass'].build_kdtree()
	#@refdirs = angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi_reduce(max_number_of_cones, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi_bcast(max_number_of_cones, 1, MPI_INT, 0, MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])


	#if( Blockdata["myid_on_node"] == 0):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	for icone in range(max_number_of_cones):
		mpi_barrier(MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				try:
					temp = prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				except IndexError:
					print(len(Blockdata['symclass_coarse'].kdneighbors))
					print(Blockdata['symclass_coarse'].kdnneigbors)
					print('ic', ic)
					print('coarse_angles', coarse_angles)
					print('len(coarse_angles)', len(coarse_angles))
					print('assign', assignments_of_refangles_to_cones[icone])
					raise
				crefim = Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				Util.Frngs(crefim, numr)
				Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
			###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
			###mpi_finalize()
			###exit()
			###  <><><><><><><><>

			#  Preprocess the data

			#  We only process images in the current cone.  icnm is consecutive number, im the actual image number
			#for icnm,im in enumerate(assignments_to_cones[icone]):
			lenass = len(assignments_to_cones[icone])
			###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
			for icnm in range(max(1,lenass)):  # I have to enter the loop even it there is no assignment
				if( lenass == 0 ):
					keepf = -1
					###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
				else:
					#  PREPARE ONE IMAGE
					im = assignments_to_cones[icone][icnm]

					phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx, sy = reduce_shifts(sx, sy, original_data[im])
						dataimage  = cyclic_shift(original_data[im],sx,sy)
						#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
						oldparams[im][3] = sx
						oldparams[im][4] = sy
						sx = 0.0
						sy = 0.0
					else:  dataimage = original_data[im].copy()


					st = get_image_statistics(dataimage, mask2D, False)
					dataimage -= st[0]
					dataimage /= st[1]
					if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
					#  Do bckgnoise if exists
					if Blockdata["bckgnoise"]:
						if apply_mask:
							if Tracker["constants"]["hardmask"]:
								dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = fft(filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
								#  Normalize bckg noise in real space, only region actually used.
								st = get_image_statistics(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = fft(dataimage)
					sig = Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
						bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = fpol(Util.mulnclreal(Util.mulnclreal(dataml, Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					sxprint( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1



				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						lxod1 = Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''
						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# DO NOT order by angular directions to save time on reprojections.

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
							##junk = time()
							xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= np.max(xod1)
						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						#if( Blockdata["myid"] == 0 ):  
						###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])


						lina = lina.reshape(lina.size)
						keepf = int(lina[-1]) + 1

						xod1 = xod1[lina]
						xod2 = xod2[lina]

						###print("  STARTING2    ",Blockdata["myid"])

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						###print("  STARTING3    ",Blockdata["myid"],lit)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["ccfpercentage"]):
								lit = j+1
								break


						###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
					#  Turn into percentage of all possible
					keepf = [int(float(keepf*100)/float(keepfirst))]
					###mpi_barrier(MPI_COMM_WORLD)
					###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
					mpi_barrier(MPI_COMM_WORLD)
					keepf = wrap_mpi_gatherv(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
					###print("  STARTING6    ",Blockdata["myid"],keepf)
					if( Blockdata["myid"] == 0 ):
						keepf = [junk for junk in keepf if junk >0]
						if( len(keepf) <2 ):
							keepf = 0
						else:
							keepf.sort()
							keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
					else:  keepf = 0
					###print("  STARTING7    ",Blockdata["myid"],keepf)
					keepf = wrap_mpi_bcast(keepf, Blockdata["main_node"], MPI_COMM_WORLD)
					if(keepf == 0):
						ERROR( "Too few images to estimate keepfirst", myid=Blockdata["myid"] )
						return
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						keepfirst = lang*m_coarse_psi*n_coarse_shifts
						keepfirst = max(keepfirst*Tracker["keepfirst"]/100,min(keepfirst,3))
						lxod1 = Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
						##'''
						assert( keepfirst == len(lxod1)/3 )
						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = np.argsort(ipsiandiang)
						xod2 = xod2[lina]  # order does not matter

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
						#if( lima<5):  sxprint("  SECOND KEPT  ",lit)


				#mpi_barrier(MPI_COMM_WORLD)
				#mpi_finalize()
				#exit()
				#for j in range(lit):
				#	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
				if( lenass > 0):
					###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

					firstdirections = [[0.0,0.0] for iln in range(lit)]
					firstshifts = [0]*lit
					for iln in range(lit):
						hashparams = int(xod2[iln])
						ishift = hashparams%1000
						ipsiandiang	= hashparams/1000
						#ipsi = ipsiandiang%100000
						iang = ipsiandiang/100000
						#try:
						firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
						#except:
						#	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					#ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = Tracker['howmany'])
					ltabang = Blockdata['symclass'].find_k_nearest_neighbors(firstdirections, k=Tracker['howmany'])
					###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

					#  Prepare image for chi2.
					#  We have to repeat everything from get shrink data, including shifts
					#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
					#    there are 3 psis, and at most n_fine_shifts. 
					#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
					#      different shifts.  If so, we have to remove duplicates from the entire set.
					lcod1 = lit*4*2*n_fine_shifts

					#cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
					#cod2.fill(-1)  #  hashparams
					cod2 = []
					#lol = 0
					for i1 in range(lit):
						hashparams = int(xod2[i1])
						ipsiandiang	= hashparams/1000
						oldiang = ipsiandiang/100000
						ipsi = ipsiandiang%100000
						ishift = hashparams%1000
						tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
						for i2 in range(Tracker['howmany']):
							iang = ltabang[i1][i2]
							if iang < 0:
								continue
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = np.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = np.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = fshift(dataml, xx, yy)
								data[ishift].set_attr("is_complex",0)
							##junk = time()
							[peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							cod1[iln] = -peak
							cod3[iln] = varadj
							iln += 1
							if(iln == lit  ):  break
						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	temp.write_image("temp.hdf")
						#	data[iln].write_image("data.hdf")
						#	ctfa.write_image("ctfa.hdf")
						#	bckgnoise.write_image("bckgnoise.hdf")
						#	exit()
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = np.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					np.exp(cod1, out=cod1)
					cod1 /= np.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
					###print("   CNORMPERPARTICLE  ",Blockdata["myid"],im,norm_per_particle[im])

					for iln in range(lit):
						newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
						#	hashparams = newpar[im][2][iln][0]
						#	ipsiandiang	= hashparams/1000
						#	ipsi = ipsiandiang%100000
						#	iang = ipsiandiang/100000
						#	ishift = hashparams%1000
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
						#	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""

	#  END OF CONES
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sxprint( "  Finished projection matching   %10.1fmin"%((time()-at)/60.))
	at = time()
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations

	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
		norm_per_particle[kl]	= sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr				+= norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))

	at = time()
	mpi_barrier(Blockdata["shared_comm"])

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi_win_free(win_vol)

	mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	mpi_barrier(Blockdata["shared_comm"])

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi_barrier(MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ):
	#	sxprint( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
	return newpar, norm_per_particle

def cerrs(params, ctfs, particle_groups):
	global Tracker, Blockdata
	from random 	import random

	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	procid = 0
	if(Blockdata["myid"] == Blockdata["nodes"][procid]):
		ref_vol = get_im(Tracker["refvol"])
		nnn = ref_vol.get_xsize()
		if(Tracker["nxinit"] != nnn ):
			ref_vol = fdecimate(ref_vol,Tracker["nxinit"],Tracker["nxinit"],Tracker["nxinit"], True, False)
	else:
		#log = None
		ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	mpi_barrier(MPI_COMM_WORLD)
	bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["nodes"][procid])
	interpolation_method = 1
	ref_vol = prep_vol(ref_vol, npad = 2, interpolation_method = interpolation_method )

	lb = Blockdata["bckgnoise"].get_xsize()
	acc_rot = 0.0
	acc_trans = 0.0
	

	#// P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
	#// exp(-4.60517) = 0.01
	pvalue = 4.60517

	for itry in range(len(params)):

		#// Get orientations (angles1) for this particle
		phi1   = params[itry][0]
		theta1 = params[itry][1]
		psi1   = params[itry][2]
		#// Get CTF for this particle
		#   Get F1 = Proj(refvol; angles1, shifts=0)
		F1 = prgl(ref_vol,[ phi1, theta1, psi1, 0.0, 0.0], interpolation_method = 1, return_real= False)
		ctfs[itry].apix = ctfs[itry].apix/shrinkage
		ct = ctf_img_real(Tracker["nxinit"], ctfs[itry])
		Util.mul_img(ct, ct)
		ctfsbckgnoise = Util.muln_img(Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [Blockdata["bckgnoise"][i,particle_groups[itry]] for i in range(lb)]), ct)

		#// Search 2 times: angles and shifts
		for imode in range(2):
			ang_error = 0.0
			sh_error = 0.0
			peak = 0.0

			#// Search for ang_error and sh_error where there are at least 3-sigma differences!
			while (peak <= pvalue):
				#// Graduallly increase the step size
				if (ang_error < 0.2):   ang_step = 0.05
				elif (ang_error < 1.):  ang_step = 0.1
				elif (ang_error < 2.):  ang_step = 0.2
				elif (ang_error < 5.):  ang_step = 0.5
				elif (ang_error < 10.): ang_step = 1.0
				elif (ang_error < 20.): ang_step = 2.0
				else:                   ang_step = 5.0

				if (sh_error < 0.2):    sh_step = 0.05
				elif (sh_error < 1.):   sh_step = 0.1
				elif (sh_error < 2.):   sh_step = 0.2
				elif (sh_error < 5.):   sh_step = 0.5
				elif (sh_error < 10.):  sh_step = 1.0
				else:                   sh_step = 2.0

				ang_error += ang_step
				sh_error  += sh_step

				#// Prevent an endless while by putting boundaries on ang_error and sh_error
				if ( (imode == 0 and ang_error > 30.) or (imode == 1 and sh_error > 10.) ):
					break

				phi2   = phi1 
				theta2 = theta1
				psi2   = psi1
				xoff1 = yoff1 = 0.0
				xshift = yshift = 0.0

				#// Perturb angle or shift , depending on the mode
				ran = random()
				if (imode == 0) :
					if (ran < 0.3333):   phi2   = phi1   + ang_error
					elif (ran < 0.6667): theta2 = theta1 + ang_error
					else:                psi2   = psi1   + ang_error
				else:
					if (ran < 0.5):
						xshift = xoff1 + sh_error
						yshift = 0.0
					else:
						xshift = 0.0
						yshift = yoff1 + sh_error

				if (imode == 0):  F2 = prgl(ref_vol,[ phi2, theta2, psi2, 0.0,0.0], 1, False)
				else:             F2 = fshift(F1, xshift*shrinkage, yshift*shrinkage)

				peak = Util.sqedac(F1, F2, ctfsbckgnoise)

			if (imode == 0):    acc_rot   += ang_error
			elif (imode == 1):  acc_trans += sh_error

	acc_rot = mpi_reduce(acc_rot, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	acc_trans = mpi_reduce(acc_trans, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	acc_rot = mpi_bcast(acc_rot, 1, MPI_FLOAT, Blockdata["main_node"], MPI_COMM_WORLD)
	acc_trans = mpi_bcast(acc_trans, 1, MPI_FLOAT, Blockdata["main_node"], MPI_COMM_WORLD)

	acc_rot = float(acc_rot[0])
	acc_trans = float(acc_trans[0])
	n_trials = Blockdata["nproc"]*len(params)

	acc_rot /= n_trials
	acc_trans /= n_trials

	if(Blockdata["myid"] == Blockdata["main_node"]):
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line,"Estimated accuracy of angles = %6.2f degrees; and shifts = %5.1f pixels"%(acc_rot,acc_trans) )

	Tracker["acc_rot"] = acc_rot
	Tracker["acc_trans"] = acc_trans

def do3d_final(partids, partstack, original_data, oldparams, oldparamstructure, projdata,\
      final_iter=-1, comm = -1 ):
	global Tracker, Blockdata

	final_dir = Tracker["directory"]
	if(Blockdata["subgroup_myid"] > -1):
		# load datastructure, read data, do two reconstructions(stepone, steptwo)
		if final_iter ==-1: final_iter = Tracker["constants"]["best"]  
		carryon = 1
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			sxprint(line, "do_final_rec3d")
			sxprint("Reconstruction uses solution of %d iteration"%final_iter)
			sxprint("Final reconstruction image size is:  %d"%(Tracker["constants"]["nnxo"]))
			sxprint("Final directory is %s"%(Tracker["directory"]))
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			try:
				refang  = read_text_row( os.path.join(final_dir, "refang.txt"))
				rshifts = read_text_row( os.path.join(final_dir, "rshifts.txt"))
			except:
				carryon =0
		else:
			refang  = 0
			rshifts = 0
		carryon = bcast_number_to_all(carryon, source_node = Blockdata["main_node"], mpi_comm = comm)
		if carryon == 0: 
			ERROR( "Failed to read refang and rshifts: %s %s " % (os.path.join(final_dir, "refang.txt"), os.path.join(final_dir, "rshifts.txt")), myid=data["subgroup_myid"] )
		refang  = wrap_mpi_bcast(refang,  Blockdata["main_node"], comm)
		rshifts = wrap_mpi_bcast(rshifts, Blockdata["main_node"], comm)

		partids =[None, None]
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			if not os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "tempdir")):
				os.mkdir(os.path.join(Tracker["constants"]["masterdir"], "tempdir"))
			l = 0
			for procid in range(2):
				partids[procid] = os.path.join(final_dir,"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
				l += len(read_text_file(partids[procid]))
		else:
			l  = 0
		l  = bcast_number_to_all(l, source_node = Blockdata["main_node"], mpi_comm = comm)

		norm_per_particle = [[],[]]
		# get the previous number of CPUs
		nproc_previous = 0
		if Blockdata["subgroup_myid"] == 0:
			while os.path.exists(os.path.join(final_dir,"oldparamstructure", \
			  "oldparamstructure_%01d_%03d_%03d.json"%(procid,nproc_previous,\
			    Tracker["mainiteration"]))):
				nproc_previous += 1
		nproc_previous = bcast_number_to_all(nproc_previous, source_node = \
		  Blockdata["main_node"], mpi_comm = comm)

		for procid in range(2):
			if procid ==0: original_data[1] = None	
			partids[procid]   = os.path.join(final_dir,"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
			partstack[procid] = os.path.join(final_dir, "params-chunk_%01d_%03d.txt"%(procid,\
			     (Tracker["mainiteration"])))
			###
			psize = len(read_text_file(partids[procid]))
			oldparamstructure[procid] = []			
			im_start, im_end = MPI_start_end(psize, Blockdata["subgroup_size"], Blockdata["subgroup_myid"])
			istart_old_proc_id = -1
			iend_old_proc_id   = -1
			plist = []
			for iproc_old in range(nproc_previous):
				im_start_old, im_end_old = MPI_start_end(psize, nproc_previous, iproc_old)
				if (im_start>= im_start_old) and im_start <=im_end_old:
					istart_old_proc_id = iproc_old
				if (im_end>= im_start_old) and im_end <=im_end_old:
					iend_old_proc_id = iproc_old
				plist.append([im_start_old, im_end_old])

			ptl_on_this_cpu = im_start
			for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id+1):
				oldparamstructure_on_old_cpu = load_object_from_json(
					os.path.join(
						final_dir,
						"oldparamstructure",
						"oldparamstructure_%01d_%03d_%03d.json"%(
							procid,
							iproc_index_old,
							Tracker["mainiteration"]
							)
						)
					)
				mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
				while (mlocal_id_on_old<len(oldparamstructure_on_old_cpu)) and (ptl_on_this_cpu<im_end):
					oldparamstructure[procid].append(oldparamstructure_on_old_cpu[mlocal_id_on_old])
					ptl_on_this_cpu  +=1
					mlocal_id_on_old +=1
			del oldparamstructure_on_old_cpu
			mpi_barrier(Blockdata["subgroup_comm"])
			#####
			original_data[procid], oldparams[procid], _ = getindexdata(
				partids[procid],
				partstack[procid],
				os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_%01d.txt"%procid),
				original_data[procid],
				small_memory=Tracker["constants"]["small_memory"],
				nproc=Blockdata["subgroup_size"],
				myid=Blockdata["subgroup_myid"],
				mpi_comm=comm
				)
			for ipar in range(len(oldparams[procid])):
				norm_per_particle[procid].append(oldparams[procid][ipar][7])
			temp = Tracker["directory"]
			Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"], "tempdir")
			mpi_barrier(Blockdata["subgroup_comm"])
			if procid ==0:
				compute_sigma([[]]*l, [[]]*l, len(oldparams[0]), True, \
				    myid = Blockdata["subgroup_myid"], mpi_comm = comm)
			Tracker["directory"] = temp
			mpi_barrier(Blockdata["subgroup_comm"])
			try:
				outlier_params = read_text_row(os.path.join(Tracker["directory"],"outlier-params-chunk_%01d_%03d.txt"%(procid, Tracker["mainiteration"])))
			except IOError:
				outlier_params = [0] * len(oldparams[procid])
			original_data_outlier = []
			oldparams_outlier = []
			norm_per_particle_outlier = []
			oldparamstructure_outlier = []
			for idx, entry in enumerate(outlier_params):
				if entry == 0:
					original_data_outlier.append(original_data[procid][idx])
					oldparams_outlier.append(oldparams[procid][idx])
					norm_per_particle_outlier.append(norm_per_particle[procid][idx])
					oldparamstructure_outlier.append(oldparamstructure[procid][idx])
			projdata[procid] = get_shrink_data(
				Tracker["constants"]["nnxo"],
				procid,
				original_data_outlier,
				oldparams_outlier,
				return_real=False,
				preshift=True,
				apply_mask=False,
				nonorm=True
				)
			oldparams[procid]      = []
			original_data[procid]  = None
			#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			line = ''
			if(Blockdata["subgroup_myid"] == Blockdata["nodes"][procid]):
				sxprint(line, "3-D reconstruction of group %d"%procid)
			###--------------------------------------------------------- Force 
			Tracker["directory"]    = Tracker["constants"]["masterdir"]
			Tracker["nxinit"]       = Tracker["constants"]["nnxo"]
			Tracker["maxfrad"]      = Tracker["constants"]["nnxo"]//2
			###---------------------------------------------------------
			do3d(
				procid,
				projdata[procid],
				oldparamstructure_outlier,
				refang,
				rshifts,
				norm_per_particle_outlier,
				myid=Blockdata["subgroup_myid"],
				smearing=True,
				mpi_comm=comm
				)
			projdata[procid]          = []
			oldparamstructure[procid] = []
			norm_per_particle[procid] = []
			mpi_barrier(Blockdata["subgroup_comm"])
		mpi_barrier(Blockdata["subgroup_comm"])
	mpi_barrier(MPI_COMM_WORLD)
	#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	line = ''
	if Blockdata["myid"] == Blockdata["main_node"]: 
		sxprint(line, "final rec3d_make_maps")
	rec3d_make_maps(compute_fsc = False, regularized = False)
	
	# also copy params to masterdir as final params
	if(Blockdata["myid"] == Blockdata["main_node"]):
		shutil.copyfile(os.path.join(final_dir, \
		  "params_%03d.txt"%Tracker["mainiteration"]), os.path.join(Tracker["constants"]["masterdir"], \
		     "final_params_%03d.txt"%Tracker["mainiteration"]))
		shutil.rmtree(os.path.join(Tracker["constants"]["masterdir"], "tempdir"))
	mpi_barrier(MPI_COMM_WORLD)
	return

def recons3d_final(masterdir, do_final_iter_init, memory_per_node, orgstack = None):
	global Tracker, Blockdata
	# search for best solution, and load respective tracker
	carryon  = 1
	#line     = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	line     = ''
	if(Blockdata["myid"] == Blockdata["main_node"]):sxprint(line, "recons3d_final")
	do_final_iter = 3
	if(do_final_iter_init ==0):
		if(Blockdata["myid"] == Blockdata["main_node"]):
			try:
				Tracker = load_tracker_from_json(os.path.join(masterdir, "Tracker_final.json"))
				sxprint("The best solution is %d  "%Tracker["constants"]["best"])
				do_final_iter =  Tracker["constants"]["best"] # set the best as do_final iteration
			except: carryon = 0
		carryon = bcast_number_to_all(carryon)
		if carryon == 0: 
			ERROR( "Best resolution is not found, do_final will not be computed", myid=Blockdata["myid"] )	# Now work on selected directory
		do_final_iter = bcast_number_to_all(do_final_iter)
	elif( do_final_iter_init == -1 ): do_final_iter = Tracker["constants"]["best"]
	else:
		do_final_iter = do_final_iter_init
		if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("User selected %d iteration to compute the 3D reconstruction "%do_final_iter)
		if do_final_iter<=2:
			ERROR( "The selected iteration should be larger than 2", myid=Blockdata["myid"] )
			
	final_dir = os.path.join(masterdir, "main%03d"%do_final_iter)
	if(Blockdata["myid"] == Blockdata["main_node"]): # check json file and load tracker
		try:
			Tracker = load_tracker_from_json(os.path.join(final_dir,"Tracker_%03d.json"%do_final_iter))
		except: carryon = 0
		if orgstack: Tracker["constants"]["stack"] = orgstack
	else: Tracker = 0
	carryon = bcast_number_to_all(carryon, Blockdata["main_node"], MPI_COMM_WORLD)
	if carryon == 0: 
		ERROR( "Failed to load Tracker file %s, program terminates " % os.path.join(final_dir,"Tracker_%03d.json"%do_final_iter), myid=Blockdata["myid"] )
	Tracker = wrap_mpi_bcast(Tracker,      Blockdata["main_node"], MPI_COMM_WORLD)
	update_legacy_tracker(Tracker)
	if(Blockdata["myid"] == Blockdata["main_node"]): # check stack 
		try:  image = get_im(Tracker["constants"]["stack"],0)
		except:carryon = 0
	carryon = bcast_number_to_all(carryon, Blockdata["main_node"], MPI_COMM_WORLD)
	if carryon == 0: 
		ERROR( "The orignal data stack for reconstruction %s does not exist, final reconstruction terminates" % Tracker["constants"]["stack"], myid=Blockdata["myid"] )

	if(Blockdata["myid"] == Blockdata["main_node"]):
		#  Estimated volume size
		volume_size = (1.5*4*(2.0*Tracker["constants"]["nnxo"]+3.0)**3)/1.e9
		#  Estimated data size
		data_size = max(Tracker["nima_per_chunk"])*4*float(Tracker["constants"]["nnxo"]**2)/float(Blockdata["no_of_groups"])/1.0e9
		nnprocs  =  min( Blockdata["no_of_processes_per_group"], int(((memory_per_node - data_size*1.2) / volume_size ) ) )
		sxprint("  MEMORY ESTIMATION.  memory per node = %6.1fGB,  volume size = %6.2fGB, data size per node = %6.2fGB, estimated number of CPUs = %d"%(memory_per_node,volume_size,data_size,nnprocs))
		if( (memory_per_node - data_size*1.2 - volume_size) < 0 or (nnprocs == 0)):  nogo = 1
		else:  nogo = 0
	else:
		nnprocs = 0
		nogo = 0
	
	nogo =    bcast_number_to_all(nogo,    source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
	if( nogo == 1 ):
		ERROR( "Insufficient memory to compute final reconstruction", myid=Blockdata["myid"] )
	nnprocs = bcast_number_to_all(nnprocs, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
	Blockdata["ncpuspernode"] 	= nnprocs
	Blockdata["nsubset"] 		= Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
	create_subgroup()

	oldparamstructure =[[],[]]
	newparamstructure =[[],[]]
	projdata          = [[model_blank(1,1)], [model_blank(1,1)]]
	original_data     = [None,None]
	oldparams         = [[],[]]
	partids           = [None, None]
	partstack         = [None, None]

	do3d_final(partids, partstack, original_data, oldparams, oldparamstructure, projdata, \
	  do_final_iter, Blockdata["subgroup_comm"])
	mpi_barrier(MPI_COMM_WORLD)
	#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	line = ''
	if(Blockdata["myid"] == Blockdata["main_node"]):
		sxprint(line, "Final reconstruction is successfully done")
	return

def recons3d_trl_struct_MPI_nosmearing(myid, main_node, prjlist, parameters, CTF, upweighted, mpi_comm, target_size):
	global Tracker, Blockdata
	from sp_utilities      import reduce_EMData_to_root, random_string, get_im, findall, model_blank, info, get_params_proj
	from EMAN2          import Reconstructors
	from sp_filter		    import filt_table
	from mpi            import MPI_COMM_WORLD, mpi_barrier
	from sp_statistics     import fsc 
	from sp_reconstruction import insert_slices_pdf
	from sp_fundamentals   import fft
	import datetime, types
	import copy
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	try:    qt = projlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	for im in range(len(prjlist)):
		ct = prjlist[im].get_attr("ctf")
		try: bckgn = prjlist[im].get_attr("bckgnoise")
		except: bckgn = [1.0]*(Tracker["constants"]["nnxo"]//2)
		if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
		prjlist[im].set_attr_dict( {"bckgnoise":bckgn, "ctf":ct})
		phi,theta,psi,s2x,s2y                   = get_params_proj(prjlist[im], xform = "xform.projection")
		junkphi, junktheta, junkpsi, js2x, js2y = parameters[im][0], parameters[im][1], parameters[im][2], parameters[im][3], parameters[im][4]
		s2x = js2x-round(js2x)
		s2y = js2y-round(js2y)
		prjlist[im] = fshift(prjlist[im], s2x, s2y)
		r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi, "theta":theta, "psi":psi}), 1.0)
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def rec3d_continuation_nosmearing(mpi_comm):
	global Tracker, Blockdata
	
	original_data	= [None, None]
	oldparams		= [None, None]
	projdata		= [None, None]
	partstack		= [None, None]
	partids			= [None, None]
	
	temp = Tracker["directory"]
	Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"])

	for procid in range(2):
		partids[procid]   = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
		partstack[procid] = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"params-chunk_%01d_%03d.txt"%(procid, Tracker["mainiteration"]))

		original_data[procid], oldparams[procid], _ = getindexdata(partids[procid], partstack[procid], \
				os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_%01d.txt"%procid), \
				original_data[procid], small_memory = Tracker["constants"]["small_memory"], \
				nproc = Blockdata["nproc"], myid = Blockdata["myid"], mpi_comm = mpi_comm)
		mpi_barrier(mpi_comm)

		projdata[procid] = get_shrink_data(Tracker["nxinit"], procid, original_data[procid], oldparams[procid],\
		  return_real = False, preshift = True, apply_mask = False, nonorm = True, nosmearing = True)
		mpi_barrier(mpi_comm)

	compute_sigma(original_data[0]+original_data[1],oldparams[0]+oldparams[1], len(oldparams[0]), False, Blockdata["myid"], mpi_comm = mpi_comm)

	# Estimate initial resolution/image size
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line,  "continuation: reconstruct initial reference")

	for procid in range(2):
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		if(Blockdata["myid"] == Blockdata["nodes"][procid]): sxprint(line, "3-D reconstruction of group %d"%procid)
		
		do3d(procid, projdata[procid], oldparams[procid], None, None, None, Blockdata["myid"], smearing = False, mpi_comm = mpi_comm)
		
		projdata[procid] = []
		mpi_barrier(mpi_comm)

	Tracker["maxfrad"] = Tracker["nxinit"]//2
	#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	line = ''
	if Blockdata["myid"] == Blockdata["main_node"]: 
		sxprint(line, "do3d_continuation_get_maps_mpi")
	rec3d_make_maps(compute_fsc = True, regularized = True)
	#if(Blockdata["myid"] == Blockdata["nodes"][0]):  shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
	Tracker["directory"] = temp
	return
	
def update_memory_estimation():
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if Tracker["constants"]["memory_per_node"] ==-1.:
			Tracker["constants"]["memory_per_node"] = Blockdata["no_of_processes_per_group"]*2.0 # reasonable approximation
		try: total_stack = EMUtil.get_image_count(Tracker["constants"]["stack"])
		except: sxprint(Tracker["constants"]["stack"], "does not exist")
		image_size   = max(Tracker["nxinit"], Tracker["constants"]["nnxo"]*1./2.)
		data_size    = total_stack*4*float(image_size**2)/float(Blockdata["no_of_groups"])/1.0e9
		volume_size  = (1.5*4*(2.0*image_size+3.0)**3)/1.e9
		#nnprocs = Blockdata["no_of_processes_per_group"]
		nnprocs = min(Blockdata["no_of_processes_per_group"], int(((Tracker["constants"]["memory_per_node"] - data_size*1.2)/volume_size)))
		sxprint("  MEMORY ESTIMATION.  memory per node = %6.1fGB,  volume size = %6.2fGB, data size per node = %6.2fGB, estimated number of CPUs = %d"%(Tracker["constants"]["memory_per_node"],volume_size,data_size,nnprocs))
		memory_per_cpu_3d = data_size/Blockdata["no_of_processes_per_group"]*1.2 + volume_size
		sxprint("  Estimated memory consumption per CPU in reconstruction %6.2f"%memory_per_cpu_3d)
		if( (Tracker["constants"]["memory_per_node"] - data_size*1.2 - volume_size) < 0 or (nnprocs == 0)):  nogo = 1
		else:  nogo = 0
	else:
		nnprocs = 0
		nogo    = 0
	nogo = bcast_number_to_all(nogo, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
	if( nogo == 1 ):  
		ERROR( "Insufficient memory to continue refinement from subset","continue_from_subset", myid=Blockdata["myid"] )
	nnprocs = bcast_number_to_all(nnprocs, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
	Blockdata["ncpuspernode"] 	= nnprocs
	Blockdata["nsubset"] 		= Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
	create_subgroup()

def update_tracker(shell_line_command):
	global Tracker, Blockdata
	# reset parameters for a restart run; update only those specified options in restart
	# 1. maxit is not included. 
	# 2. those sigmas for local search can be considered included 
	from optparse import OptionParser			
	parser_no_default = OptionParser()
	parser_no_default.add_option("--radius",      		   		type= "int")
	parser_no_default.add_option("--xr",      		       		type="float")
	parser_no_default.add_option("--ts",      		       		type="float")
	parser_no_default.add_option("--inires",		       		type="float")
	parser_no_default.add_option("--mask3D",		         	type="string")
	parser_no_default.add_option("--function",					type="string")
	parser_no_default.add_option("--function_ai",				type="string")
	parser_no_default.add_option("--symmetry",    	       		type="string")# rare to change sym; however, keep it an option.
	parser_no_default.add_option("--delta",						type="float")
	parser_no_default.add_option("--shake",	           			type="float")
	parser_no_default.add_option("--small_memory",         		action="store_true")
	parser_no_default.add_option("--ccfpercentage",		 		type="float")
	parser_no_default.add_option("--nonorm",               		action="store_true")
	parser_no_default.add_option("--memory_per_node",         	type="float")
	parser_no_default.add_option("--an",                        type="float")
	parser_no_default.add_option("--limit_improvement",                        type="int")
	parser_no_default.add_option("--a_criterion",                        type="float")

	parser_no_default.add_option("--do_final",             		type="int")
	parser_no_default.add_option("--local_refinement",			action="store_true")
	parser_no_default.add_option("--continuation_orgstack",     type="string")
	parser_no_default.add_option("--continuation_initvol",      type="string")
	parser_no_default.add_option("--subset",                    type="string")
	parser_no_default.add_option("--oldrefdir",                 type="string")
	parser_no_default.add_option("--continuation_iter",         type="int")
	parser_no_default.add_option("--continuation_smearing",     type="int")
	parser_no_default.add_option("--keep_groups",               action="store_true")
	parser_no_default.add_option("--plot_ang_dist",               action="store_true")

	(options_no_default_value, args) = parser_no_default.parse_args(shell_line_command)

	#  This is for printout only
	tempdict = {}

	if options_no_default_value.local_refinement:
		Tracker['constants']['do_local'] = True

	if 	options_no_default_value.plot_ang_dist != None:
		Tracker["constants"]["plot_ang_dist"] 				= options_no_default_value.plot_ang_dist

	if 	options_no_default_value.limit_improvement != None:
		Tracker["constants"]["limit_improvement"] 				= options_no_default_value.limit_improvement

	if 	options_no_default_value.a_criterion != None:
		Tracker["constants"]["a_criterion"] 				= options_no_default_value.a_criterion

	if 	options_no_default_value.radius != None:
		Tracker["constants"]["radius"] 				= options_no_default_value.radius
		tempdict["radius"] = Tracker["constants"]["radius"]

	if options_no_default_value.xr != None:
		Tracker["xr"] 										= options_no_default_value.xr
		tempdict["xr"] = Tracker["xr"]

	if options_no_default_value.ts != None:
		Tracker["ts"] 										= options_no_default_value.ts
		tempdict["ts"] = Tracker["ts"]

	if options_no_default_value.inires != None:
		Tracker["constants"]["inires"] 						= options_no_default_value.inires
		Tracker["constants"]["inires"]= int(Tracker["constants"]["nnxo"]*Tracker["constants"]["pixel_size"]/Tracker["constants"]["inires"] + 0.5)
		Tracker["currentres"] = Tracker["constants"]["inires"]
		tempdict["currentres"] = Tracker["currentres"]

	if options_no_default_value.delta != None:				
		Tracker["delta"] 									= options_no_default_value.delta
		tempdict["delta"] = Tracker["delta"]
	if options_no_default_value.shake != None:
		Tracker["constants"]["shake"] 						= options_no_default_value.shake
		tempdict["shake"] = Tracker["constants"]["shake"]
	if options_no_default_value.symmetry != None:  # this rarely happens. However, keep it an option.
		sym    												= options_no_default_value.symmetry
		Tracker["constants"]["symmetry"]					= sym[0].lower() + sym[1:] 
		tempdict["symmetry"] = Tracker["constants"]["symmetry"]
	if options_no_default_value.mask3D != None:
		Tracker["constants"]["mask3D"] 						= options_no_default_value.mask3D
		tempdict["mask3D"] = Tracker["constants"]["mask3D"]
	if options_no_default_value.ccfpercentage != None:
		Tracker["ccfpercentage"] 				= options_no_default_value.ccfpercentage/100.0
		tempdict["ccfpercentage"] = Tracker["ccfpercentage"]
	if options_no_default_value.nonorm != None:
		Tracker["constants"]["nonorm"] 						= options_no_default_value.nonorm
		tempdict["nonorm"] = Tracker["constants"]["nonorm"]
	if options_no_default_value.small_memory != None:
		Tracker["constants"]["small_memory"] 				= options_no_default_value.small_memory
		tempdict["small_memory"] = Tracker["constants"]["small_memory"]
	if options_no_default_value.memory_per_node != None:
		Tracker["constants"]["memory_per_node"] = options_no_default_value.memory_per_node
		tempdict["memory_per_node"] = Tracker["constants"]["memory_per_node"]
	### continuation	
	if  options_no_default_value.continuation_orgstack != None:
		Tracker["constants"]["continuation_orgstack"] = options_no_default_value.continuation_orgstack
	if  options_no_default_value.keep_groups != None:
		Tracker["constants"]["keep_groups"] = options_no_default_value.keep_groups
	if  options_no_default_value.subset != None:
		Tracker["constants"]["subset"] = options_no_default_value.subset	
	if  options_no_default_value.oldrefdir != None:
		Tracker["constants"]["oldrefdir"] = options_no_default_value.oldrefdir
	if  options_no_default_value.continuation_initvol != None:
		Tracker["constants"]["continuation_initvol"] = options_no_default_value.continuation_initvol
	if  options_no_default_value.continuation_iter != None:
		Tracker["constants"]["continuation_iter"] = options_no_default_value.continuation_iter
	if  options_no_default_value.function != None:
		Tracker["constants"]["user_func_volume"] = options_no_default_value.function
	if  options_no_default_value.function_ai != None:
		Tracker["constants"]["user_func_ai"] = options_no_default_value.function_ai
	if  options_no_default_value.an!= None:
		Tracker["constants"]["an"] = options_no_default_value.an
		tempdict["an"] = Tracker["constants"]["an"]

	# For backwards compatibility
	update_legacy_tracker(Tracker)

	if( (Blockdata["myid"] == Blockdata["main_node"])  and  (len(tempdict) > 0) ):
		print_dict(tempdict, "Updated settings")

	Blockdata["symclass"] = symclass(Tracker["constants"]["symmetry"])
	Blockdata["symclass_coarse"] = symclass(Tracker["constants"]["symmetry"])
	return 


def update_legacy_tracker(tracker):
	backwards_dict_constants = {
		'user_func_volume': 'do_volume_mask',
		'user_func_ai': 'ai_spa',
		'even_angle_method': 'S',
		'group_id': None,
		'filament_width': None,
		'helical_rise': None,
		'stack_prior': None,
		'stack_prior_fmt': None,
		'stack_prior_dtype': None,
		'do_local': False,
		'a_criterion': 0.75,
		'plot_ang_dist': False,
		}
	backwards_dict = {
		'theta_min': -1,
		'theta_max': -1,
		'ccfpercentage': -1,
		'howmany': 4,
		'prior': {},
		}
	prior_dict = {
		'tol_psi': 30,
		'tol_theta': 15,
		'tol_filament': 0.2,
		'tol_std': 1,
		'tol_mean': 30,
		'outlier_method': 'deg',
		'prior_method': 'running',
		'force_outlier': False,
		'remove_outlier': False,
		'window_size': 3,
		}
	for key in backwards_dict:
		try:
			tracker[key]
		except KeyError:
			tracker[key] = backwards_dict[key]
	for key in backwards_dict_constants:
		try:
			tracker["constants"][key]
		except KeyError:
			tracker["constants"][key] = backwards_dict_constants[key]
	for key in prior_dict:
		try:
			tracker["prior"][key]
		except KeyError:
			tracker["prior"][key] = prior_dict[key]


def compare_bckgnoise(bckgnoise1, bckgnoise2):
	ccsum = 0.0
	nx =  bckgnoise1.get_xsize()
	ny =  bckgnoise1.get_ysize()
	for i in range(ny):
		c1 = []
		c2 = []
		for j in range(1, nx):
			c1.append(bckgnoise1.get_value_at(j,i, 0))
			c2.append(bckgnoise2.get_value_at(j,i, 0))
		ccsum +=pearson(c1, c2)
	return ccsum/float(ny)
	
def rec3d_make_maps(compute_fsc = True, regularized = True):
	global Tracker, Blockdata
	
	# final reconstruction: compute_fsc = False; regularized = False
	# tempdir is removed in the end of the function
	if compute_fsc:
		if(Blockdata["no_of_groups"] == 1):
			if( Blockdata["myid"] == Blockdata["nodes"][0] ):
				tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%03d.hdf"%(   Tracker["mainiteration"])))
				tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%03d.hdf"%(Tracker["mainiteration"])))
				tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%03d.hdf"%(   Tracker["mainiteration"])))
				tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%03d.hdf"%(Tracker["mainiteration"])))
				Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
				shrank0 	= stepone(tvol0, tweight0)
				shrank1 	= stepone(tvol1, tweight1)
				#  Note shrank volumes are Fourier uncentered.
				cfsc 		= fsc(shrank0, shrank1)[1]
				del shrank0, shrank1
				if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
					cfsc 	= cfsc[:Tracker["nxinit"]//2+1]
					for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
				lcfsc = len(cfsc)
				#--  memory_check(Blockdata["myid"],"second node, after stepone")
			else:
				#  receive fsc
				lcfsc = 0
		else:
			if(Blockdata["myid"] == Blockdata["nodes"][1]):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
				#--  memory_check(Blockdata["myid"],"first node, before stepone")
				#  read volumes, shrink
				tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%03d.hdf"%(Tracker["mainiteration"])))
				tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%03d.hdf"%(Tracker["mainiteration"])))
				tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%03d.hdf"%(Tracker["mainiteration"])))
				tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%03d.hdf"%(Tracker["mainiteration"])))
				Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
				tag = 7007
				send_EMData(tvol1,    Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
				send_EMData(tweight1, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
				shrank0 	= stepone(tvol0, tweight0)
				send_EMData(shrank0,  Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
				del shrank0
				lcfsc = 0
				#--  memory_check(Blockdata["myid"],"first node, after stepone")
			elif( Blockdata["myid"] == Blockdata["nodes"][0] ):
				#--  memory_check(Blockdata["myid"],"second node, before stepone")
				#  read volumes, shrink
				tag = 7007
				tvol1 		= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
				tweight1 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
				tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
				shrank1 	= stepone(tvol1, tweight1)
				#  Get shrank volume, do fsc, send it to all
				shrank0 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
				#  Note shrank volumes are Fourier uncentered.
				cfsc 		= fsc(shrank0, shrank1)[1]
				del shrank0, shrank1
				if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
					cfsc = cfsc[:Tracker["nxinit"]//2+1]
					for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
				lcfsc = len(cfsc)
				#--  memory_check(Blockdata["myid"],"second node, after stepone")
			else:
				#  receive fsc
				lcfsc = 0
		mpi_barrier(MPI_COMM_WORLD)
		lcfsc = bcast_number_to_all(lcfsc)
		if( Blockdata["myid"] != Blockdata["nodes"][0]  ):  cfsc = [0.0]*lcfsc
		cfsc = bcast_list_to_all(cfsc, Blockdata["myid"], Blockdata["nodes"][0])
		if( Blockdata["myid"] == Blockdata["main_node"]):
			write_text_file(cfsc, os.path.join(Tracker["directory"] ,"driver_%03d.txt"%(Tracker["mainiteration"])))
			out_fsc(cfsc)

	# Now that we have the curve, do the reconstruction
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	if (Blockdata["no_of_groups"] > 1):    lorder = [0,0] #  Two blocks in parallel
	elif( Blockdata["no_of_groups"] == 1): lorder = [0,1] #  One after another
	
	if regularized:
		for iorder in range(2):
			if(iorder == lorder[0]):
				if(Blockdata["color"] == Blockdata["node_volume"][1]):
					#--  memory_check(Blockdata["myid"],"first node, before steptwo")
					#  compute filtered volume
					if( Blockdata["myid_on_node"] == 0 ):
						treg0 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_0_%03d.hdf"%(Tracker["mainiteration"])))
					else:

						tvol0    = model_blank(1)
						tweight0 = model_blank(1)
						treg0    = model_blank(1)
					tvol0 = steptwo_mpi(tvol0, tweight0, treg0, cfsc, True, color = Blockdata["node_volume"][1])		
					del tweight0, treg0
					if( Blockdata["myid_on_node"] == 0 ):
						#--  memory_check(Blockdata["myid"],"first node, before masking")
						if( Tracker["mainiteration"] == 1 and Tracker['constants']['inires'] != -1):
							# At a first iteration truncate resolution at the initial resolution set by the user
							for i in range(len(cfsc)):
								if(  i < Tracker["constants"]["inires"]+1 ):  cfsc[i]   = 1.0
								if(  i == Tracker["constants"]["inires"]+1 ): cfsc[i]  	= 0.5
								elif( i > Tracker["constants"]["inires"]+1 ): cfsc[i]  	= 0.0
							tvol0 = filt_table(tvol0, cfsc)
							if( Blockdata["no_of_groups"] > 1 ):  del cfsc

						user_func = sp_user_functions.factory[Tracker["constants"]["user_func_volume"]]
						#ref_data = [tvol0, Tracker, mainiteration]
						ref_data = [tvol0, Tracker, Tracker["mainiteration"]]
						#--  #--  memory_check(Blockdata["myid"],"first node, after masking")
						user_func(ref_data).write_image(os.path.join(Tracker["directory"], \
							"vol_0_%03d.hdf"%(Tracker["mainiteration"])))
					del tvol0
					#--  memory_check(Blockdata["myid"],"first node, after 2 steptwo")
			if(iorder == lorder[1]):
				if( Blockdata["color"] == Blockdata["node_volume"][0] ):
					#--  memory_check(Blockdata["myid"],"second node, before steptwo")
					#  compute filtered volume
					if( Blockdata["myid_on_node"] == 0 ):
						treg1 = get_im(os.path.join(Tracker["directory"], "tempdir",\
						   "trol_1_%03d.hdf"%(Tracker["mainiteration"])))
					else:
						tvol1    = model_blank(1)
						tweight1 = model_blank(1)
						treg1    = model_blank(1)
					
					tvol1 = steptwo_mpi(tvol1, tweight1, treg1, cfsc, True,  color = Blockdata["node_volume"][0])
					del tweight1, treg1
					if( Blockdata["myid_on_node"] == 0 ):
						#--  memory_check(Blockdata["myid"],"second node, before masking")
						if( Tracker["mainiteration"] == 1 and Tracker['constants']['inires'] != -1):
							# At a first iteration truncate resolution at the initial resolution set by the user
							for i in range(len(cfsc)):
								if(  i < Tracker["constants"]["inires"]+1 ):   cfsc[i]  = 1.0
								if(  i == Tracker["constants"]["inires"]+1 ):  cfsc[i]  = 0.5
								elif( i > Tracker["constants"]["inires"]+1 ):  cfsc[i]  = 0.0
							tvol1 = filt_table(tvol1, cfsc)
							del cfsc
						user_func = sp_user_functions.factory[Tracker["constants"]["user_func_volume"]]

						#ref_data = [tvol1, Tracker, mainiteration]
						ref_data = [tvol1, Tracker, Tracker["mainiteration"]]
						#--  #--  memory_check(Blockdata["myid"],"first node, after masking")
						user_func(ref_data).write_image(os.path.join(Tracker["directory"], \
						   "vol_1_%03d.hdf"%(Tracker["mainiteration"])))
					del tvol1
					#--  memory_check(Blockdata["myid"],"second node, after 2 steptwo")
			#  Here end per node execution.
		if( Blockdata["myid"] == Blockdata["nodes"][0]): shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
		mpi_barrier(MPI_COMM_WORLD)
	else:
		if(Blockdata["no_of_groups"] == 1):
			for iproc in range(2):
				if(Blockdata["myid_on_node"] == 0):
					tvol0 		= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tvol_0_%03d.hdf"%(Tracker["mainiteration"]))))
					tweight0 	= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir","tweight_0_%03d.hdf"%(Tracker["mainiteration"]))))
					tvol1 		= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tvol_1_%03d.hdf"%(Tracker["mainiteration"]))))
					tweight1 	= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir","tweight_1_%03d.hdf"%(Tracker["mainiteration"]))))
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
					treg  = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%03d.hdf"%((iproc, Tracker["mainiteration"]))))

				else:
					treg = model_blank(1)
					if iproc ==0:
						tvol0 		= model_blank(1)
						tweight0 	= model_blank(1)
					else:
						tvol1 		= model_blank(1)
						tweight1 	= model_blank(1)
				if iproc ==0 : 
					tvol0 = steptwo_mpi(tvol0, tweight0, treg, None, False , color = Blockdata["node_volume"][0])
				else:
					tvol1 = steptwo_mpi(tvol1, tweight1, treg, None, False , color = Blockdata["node_volume"][0])
				if( Blockdata["myid_on_node"] == 0):
					if iproc ==0: tvol0.write_image(os.path.join(Tracker["constants"]["masterdir"], "vol_%d_unfil_%03d.hdf"%(iproc, final_iter)))
					else: tvol1.write_image(os.path.join(Tracker["constants"]["masterdir"],         "vol_%d_unfil_%03d.hdf"%(iproc, final_iter)))
				mpi_barrier(MPI_COMM_WORLD)
		else:
			if(Blockdata["myid"] == Blockdata["main_shared_nodes"][1]):
				# post-insertion operations, done only in main_node		
				tvol0 		= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tvol_0_%03d.hdf"%   Tracker["mainiteration"])))
				tweight0 	= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tweight_0_%03d.hdf"%Tracker["mainiteration"])))
				tvol1 		= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tvol_1_%03d.hdf"%   Tracker["mainiteration"])))
				tweight1 	= get_im(os.path.join(Tracker["directory"],os.path.join("tempdir", "tweight_1_%03d.hdf"%Tracker["mainiteration"])))
				Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
			mpi_barrier(MPI_COMM_WORLD)
			
			if(Blockdata["myid"] == Blockdata["main_shared_nodes"][1]):
				tag = 7007
				send_EMData(tvol1,    Blockdata["main_shared_nodes"][0], tag, MPI_COMM_WORLD)
				send_EMData(tweight1, Blockdata["main_shared_nodes"][0], tag, MPI_COMM_WORLD)
				tvol0.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
				
			elif(Blockdata["myid"] == Blockdata["main_shared_nodes"][0]):
				tag = 7007
				tvol1    	= recv_EMData(Blockdata["main_shared_nodes"][1], tag, MPI_COMM_WORLD)
				tweight1    = recv_EMData(Blockdata["main_shared_nodes"][1], tag, MPI_COMM_WORLD)
				tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			mpi_barrier(MPI_COMM_WORLD)
			
			if( Blockdata["color"] == Blockdata["node_volume"][1]):
				if( Blockdata["myid"] == Blockdata["main_shared_nodes"][1] ):
					treg0 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_0_%03d.hdf"%(Tracker["mainiteration"])))
				else:
					tvol0 		= model_blank(1)
					tweight0 	= model_blank(1)
					treg0 		= model_blank(1)
				tvol0 = steptwo_mpi(tvol0, tweight0, treg0, None, False , color = Blockdata["node_volume"][1])
				del tweight0, treg0
				if( Blockdata["myid_on_node"] == 0):
					tvol0.write_image(os.path.join(Tracker["constants"]["masterdir"], "vol_0_unfil_%03d.hdf"%Tracker["mainiteration"]))
					
			elif( Blockdata["color"] == Blockdata["node_volume"][0]):
				if( Blockdata["myid"] == Blockdata["main_shared_nodes"][0]):
					treg1 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_1_%03d.hdf"%(Tracker["mainiteration"])))
				else:
					tvol1 		= model_blank(1)
					tweight1 	= model_blank(1)
					treg1 		= model_blank(1)
				tvol1 = steptwo_mpi(tvol1, tweight1, treg1, None, False , color = Blockdata["node_volume"][0])
				del tweight1, treg1
				if( Blockdata["myid_on_node"] == 0):
					tvol1.write_image(os.path.join(Tracker["constants"]["masterdir"], "vol_1_unfil_%03d.hdf"%Tracker["mainiteration"]))
			mpi_barrier(MPI_COMM_WORLD)
	return
	
def refinement_one_iteration(partids, partstack, original_data, oldparams, projdata, continuation_mode = False):
	global Tracker, Blockdata
	#  READ DATA AND COMPUTE SIGMA2   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	if Blockdata['myid'] == Blockdata['main_node']:
		sxprint('Import particle stack')
	start_end_list = [[0, 0], [0, 0]]
	for procid in range(2):
		original_data[procid], oldparams[procid], start_end_list[procid] = getindexdata(
			partids[procid],
			partstack[procid],
			os.path.join(
				Tracker["constants"]["masterdir"],
				"main000",
				"particle_groups_%01d.txt"%procid
				),
			original_data[procid],
			small_memory=Tracker["constants"]["small_memory"],
			nproc=Blockdata["nproc"],
			myid=Blockdata["myid"],
			mpi_comm=MPI_COMM_WORLD,
			)

	mpi_barrier(MPI_COMM_WORLD)

	if Tracker["mainiteration"] == 1 and not continuation_mode:
		dryrun = False
	else:
		dryrun = True

	compute_sigma(original_data[0]+original_data[1], oldparams[0]+oldparams[1], len(oldparams[0]), dryrun, Blockdata["myid"])

	#  REFINEMENT   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

	mpi_barrier(MPI_COMM_WORLD)

	refang, rshifts, coarse_angles, coarse_shifts = get_refangs_and_shifts()
	if( Tracker["constants"]["shake"] > 0.0 ):
		if(Blockdata["myid"] == Blockdata["main_node"]):
			shakenumber = uniform( -Tracker["constants"]["shake"], Tracker["constants"]["shake"])
		else:
			shakenumber = 0.0
		shakenumber = bcast_number_to_all(shakenumber, source_node = Blockdata["main_node"])
		# it has to be rounded as the number written to the disk is rounded,
		#  so if there is discrepancy one cannot reproduce iteration.
		shakenumber = round(shakenumber,5)

		rangle        = shakenumber*Tracker["delta"]
		rshift        = shakenumber*Tracker["ts"]
		refang        = Blockdata["symclass"].reduce_anglesets( rotate_params(refang, [-rangle,-rangle,-rangle]) )
		coarse_angles = Blockdata["symclass"].reduce_anglesets( rotate_params(coarse_angles, [-rangle,-rangle,-rangle]) )
		shakegrid(rshifts, rshift)
		shakegrid(coarse_shifts, rshift)

		if(Blockdata["myid"] == Blockdata["main_node"]):
			write_text_row([[shakenumber, rangle, rshift]], os.path.join(Tracker["directory"] ,"randomize_search.txt") )
	else:
		rangle = 0.0
		rshift = 0.0

	if(Blockdata["myid"] == Blockdata["main_node"]):
		write_text_row( refang,  os.path.join(Tracker["directory"] , "refang.txt"))
		write_text_row( rshifts, os.path.join(Tracker["directory"] , "rshifts.txt"))
	mpi_barrier(MPI_COMM_WORLD)

	newparamstructure = [[],[]]
	raw_vol           = [[],[]]
	norm_per_particle = [[],[]]
	outliers = [[], []]

	for procid in range(2):
		Tracker["refvol"] = os.path.join(Tracker["previousoutputdir"],"vol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"]-1))

		Tracker["nxpolar"] = Tracker["nxinit"]#min( 3*Tracker["nxinit"], Tracker["constants"]["nnxo"] )
		#Tracker["nxpolar"] = min( 2*Tracker["nxinit"], Tracker["constants"]["nnxo"] )
		if( Tracker["state"] == "INITIAL" ):
			newparamstructure[procid], norm_per_particle[procid] = \
					ali3D_polar_ccc(refang, rshifts, coarse_angles, coarse_shifts, procid, original_data[procid], oldparams[procid], \
					preshift = True, apply_mask = True, nonorm = Tracker["constants"]["nonorm"], applyctf = False)
		elif( Tracker["state"] == "PRIMARY" ):
			newparamstructure[procid], norm_per_particle[procid] = \
					ali3D_primary_polar(refang, rshifts, coarse_angles, coarse_shifts, procid, original_data[procid], oldparams[procid], \
					preshift = True, apply_mask = True, nonorm = Tracker["constants"]["nonorm"], applyctf = True)
		elif( (Tracker["state"] == "PRIMARY LOCAL") ):
			###print("   ",Blockdata["myid"],len(refang),len(rshifts),len(coarse_angles),len(coarse_shifts),len(original_data[procid]), len(oldparams[procid]))
			###print("   ",Blockdata["myid"],original_data[0][0], oldparams[0][0])
			newparamstructure[procid], norm_per_particle[procid] = \
					ali3D_primary_local_polar(refang, rshifts, coarse_angles, coarse_shifts, procid, original_data[procid], oldparams[procid], \
					preshift = True, apply_mask = True, nonorm = Tracker["constants"]["nonorm"], applyctf = True)
		elif( Tracker["state"] == "EXHAUSTIVE" ):
			newparamstructure[procid], norm_per_particle[procid] = \
					ali3D_polar(refang, rshifts, coarse_angles, coarse_shifts, procid, original_data[procid], oldparams[procid], \
					preshift = True, apply_mask = True, nonorm = Tracker["constants"]["nonorm"], applyctf = True)
		elif( (Tracker["state"] == "RESTRICTED") or (Tracker["state"] == "FINAL") ):
			newparamstructure[procid], norm_per_particle[procid] = \
					ali3D_local_polar(refang, rshifts, coarse_angles, coarse_shifts, procid, original_data[procid], oldparams[procid], \
					preshift = True, apply_mask = True, nonorm = Tracker["constants"]["nonorm"], applyctf = True)
		else:  
			ERROR( "Incorrect state  %s" % Tracker["state"], myid=Blockdata["myid"] )


		qt = 1.0/Tracker["constants"]["nnxo"]/Tracker["constants"]["nnxo"]
		params = []
		for im in range(len(newparamstructure[procid])):
			#  Select only one best
			hash   = newparamstructure[procid][im][2][0][0]
			ishift = hash%1000
			ipsi   = (hash/1000)%100000
			iang   = hash/100000000
			params.append([ refang[iang][0], refang[iang][1], (refang[iang][2]+ipsi*Tracker["delta"])%360.0, \
			  rshifts[ishift][0]+oldparams[procid][im][3], rshifts[ishift][1]+oldparams[procid][im][4],\
			     newparamstructure[procid][im][-1][0][1], norm_per_particle[procid][im]*qt, norm_per_particle[procid][im]])

		mpi_barrier(MPI_COMM_WORLD)

		params = wrap_mpi_gatherv(params, Blockdata["main_node"], MPI_COMM_WORLD)
		#  store params
		if(Blockdata["myid"] == Blockdata["main_node"]):
			line = ''
			sxprint(line,"Executed successfully: ","Projection matching, state: %s, number of images:%7d"%(Tracker["state"],len(params)))
			write_text_row(params, os.path.join(Tracker["directory"], "params-chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"])) )
			cmd = "{} {} {}".format(
				"cp -p",
				os.path.join(
					Tracker["previousoutputdir"],
					"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]-1)
					),
				os.path.join(
					Tracker["directory"],
					"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"])
					)
				)
			junk = cmdexecute(cmd)
		del params

		projdata[procid] = []
		if Tracker["constants"]["small_memory"]:
			original_data[procid], oldparams[procid], start_end_list[procid] = getindexdata(partids[procid], partstack[procid], \
			os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_%01d.txt"%procid), \
			original_data[procid], small_memory = Tracker["constants"]["small_memory"], \
			nproc = Blockdata["nproc"], myid = Blockdata["myid"], mpi_comm = MPI_COMM_WORLD)

		# Prepare outlier calculation
		partids[procid] = os.path.join(Tracker["directory"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
		partstack[procid] = os.path.join(Tracker["directory"],"params-chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
		outlier_file = os.path.join(Tracker["directory"],"outlier-params-chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
		outlier_list, outlier_full = calculate_prior_values(
			tracker=Tracker,
			blockdata=Blockdata,
			outlier_file=outlier_file,
			params_file=partstack[procid],
			chunk_file=partids[procid],
			im_start=start_end_list[procid][0],
			im_end=start_end_list[procid][1],
			procid=procid
			)
		outliers[procid] = outlier_list

		original_outlier = []
		newparams_outlier = []
		norm_outlier = []
		oldparams_outlier = []
		for idx, entry in enumerate(outlier_list):
			if entry == 0:
				original_outlier.append(original_data[procid][idx])
				newparams_outlier.append(newparamstructure[procid][idx])
				oldparams_outlier.append(oldparams[procid][idx])
				norm_outlier.append(norm_per_particle[procid][idx])

		if Tracker["changed_delta"]:
			org_nxinit        = Tracker["nxinit"]
			Tracker["nxinit"] = Tracker["constants"]["nnxo"]

		projdata[procid] = get_shrink_data(
			Tracker["nxinit"],
			procid,
			original_outlier,
			oldparams_outlier,
			return_real=False,
			preshift=True,
			apply_mask=False,
			nonorm=True
			)

		oldparams_outlier = []
		oldparams[procid] = []
		if Tracker["constants"]["small_memory"]:
			original_data[procid]	= []
			original_outlier = []

		if Blockdata['myid'] == Blockdata['main_node'] and Tracker['constants']['plot_ang_dist']:
			### NEEDS TO BE REACTIVATED AFTER THE SYMCLASS CHANGE
			sxprint('Create angular distribution plot for chunk {0}'.format(procid))
			delta = np.maximum(Tracker['delta'], 3.75)
			exclude = []
			exclude.append(
				[
					None,
					os.path.join(Tracker["directory"], "ang_dist_{0}".format(procid)),
					'',
					]
				)
			if np.array(outlier_full).any():
				exclude.append(
					[
						outlier_full,
						os.path.join(Tracker["directory"], "ang_dist_{0}_outlier".format(procid)),
						'_outlier',
						]
					)

			for exclude_list, dir_name, suffix in exclude:
				sp_utilities.angular_distribution(
					params_file=partstack[procid],
					output_folder=dir_name,
					prefix='ang_dist{0}'.format(suffix),
					method=Tracker['constants']['even_angle_method'],
					pixel_size=1,
					delta=delta,
					symmetry=Tracker['constants']['symmetry'],
					box_size=Tracker['constants']['nnxo'],
					particle_radius=Tracker['constants']['radius'],
					dpi=72,
					do_print=False,
					exclude=exclude_list,
					)
				sp_utilities.angular_distribution(
					params_file=partstack[procid],
					output_folder=dir_name,
					prefix='full_ang_dist{0}'.format(suffix),
					method=Tracker['constants']['even_angle_method'],
					pixel_size=1,
					delta=delta,
					symmetry=Tracker['constants']['symmetry'] + '_full',
					box_size=Tracker['constants']['nnxo'],
					particle_radius=Tracker['constants']['radius'],
					dpi=72,
					do_print=False,
					exclude=exclude_list,
					)
		else:
			del outlier_full

		if( Blockdata["myid"] == Blockdata['main_node'] ):
			sxprint('Do 3D reconstruction')
		do3d(
			procid,
			projdata[procid],
			newparams_outlier,
			refang,
			rshifts,
			norm_outlier,
			Blockdata["myid"],
			smearing=True,
			mpi_comm=MPI_COMM_WORLD
			)
		projdata[procid] = []
		if Tracker["changed_delta"]:
			Tracker["nxinit"] = org_nxinit

		if( Blockdata["myid_on_node"] == 0 ):
			for kproc in range(Blockdata["no_of_processes_per_group"]):
				if( kproc == 0 ):
					fout = open(os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid,Blockdata["myid"],Tracker["mainiteration"])),'w')
					json.dump(newparamstructure[procid], fout)
					fout.close()
				else:
					dummy = wrap_mpi_recv(kproc, Blockdata["shared_comm"])
					fout = open(os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid,(Blockdata["color"]*Blockdata["no_of_processes_per_group"] + kproc),Tracker["mainiteration"])),'w')
					json.dump(dummy, fout)
					fout.close()
					del dummy
		else:
			wrap_mpi_send(newparamstructure[procid], 0, Blockdata["shared_comm"])

		###fout = open(os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid,Blockdata["myid"],Tracker["mainiteration"])),'w')
		###json.dump(newparamstructure[procid], fout)
		###fout.close()
		newparamstructure[procid] = []
		norm_per_particle[procid] = []
		original_outlier = []
		newparams_outlier = []
		norm_outlier = []
		mpi_barrier(MPI_COMM_WORLD)

	del refang, rshifts

	#  DRIVER RESOLUTION ASSESSMENT and RECONSTRUCTION <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	if Tracker["changed_delta"]:
		org_nxinit = Tracker["nxinit"]
		Tracker["nxinit"] = Tracker["constants"]["nnxo"]

	rec3d_make_maps(compute_fsc = True, regularized = True)

	if Tracker["changed_delta"]:
		Tracker["nxinit"] = org_nxinit

	#from sys import exit
	#mpi_finalize()
	#exit()
	#
	#  Change to current params
	if( Blockdata["myid"] == Blockdata["main_node"]):
		# Carry over chunk information
		for procid in range(2):
			cmd = "{} {} {}".format("cp -p", os.path.join(Tracker["previousoutputdir"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]-1)), \
									os.path.join(Tracker["directory"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"])) )
			junk = cmdexecute(cmd)

		pinids = read_text_file(partids[0])  + read_text_file(partids[1])
		params = read_text_row(partstack[0]) + read_text_row(partstack[1])

		assert(len(pinids) == len(params))

		for i in range(len(pinids)):
			pinids[i] = [ pinids[i], params[i] ]
		del params
		pinids.sort()

		write_text_file([pinids[i][0] for i in range(len(pinids))], os.path.join(Tracker["directory"] ,"indexes_%03d.txt"%(Tracker["mainiteration"])))
		write_text_row( [pinids[i][1] for i in range(len(pinids))], os.path.join(Tracker["directory"] ,"params_%03d.txt"%(Tracker["mainiteration"])))
		del pinids
	mpi_barrier(MPI_COMM_WORLD)

	if(Tracker["mainiteration"] == 1 ):
		acc_rot = acc_trans = 1.e23
	else:
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			Blockdata["bckgnoise"]= get_im(os.path.join(Tracker["directory"],"bckgnoise.hdf"))
			nnx = Blockdata["bckgnoise"].get_xsize()
			nny = Blockdata["bckgnoise"].get_ysize()
		else:
			nnx = 0
			nny = 0
		nnx = bcast_number_to_all(nnx)
		nny = bcast_number_to_all(nny)
		if( Blockdata["myid"] != Blockdata["main_node"] ):
			Blockdata["bckgnoise"] = model_blank(nnx,nny, 1, 1.0)
		bcast_EMData_to_all(Blockdata["bckgnoise"], Blockdata["myid"], source_node = Blockdata["main_node"])

		if(Blockdata["myid"] == Blockdata["main_node"]):
			params = read_text_row(os.path.join(Tracker["directory"],"params-chunk_0_%03d.txt"%(Tracker["mainiteration"])))+read_text_row(os.path.join(Tracker["directory"],"params-chunk_1_%03d.txt"%(Tracker["mainiteration"])))
			try:
				outlier_params = read_text_row(os.path.join(Tracker["directory"],"outlier-params-chunk_0_%03d.txt"%(Tracker["mainiteration"])))+read_text_row(os.path.join(Tracker["directory"],"outlier-params-chunk_1_%03d.txt"%(Tracker["mainiteration"])))
			except IOError:
				outlier_params = [0] * len(params)
			li     = read_text_file(os.path.join(Tracker["directory"],"chunk_0_%03d.txt"%(Tracker["mainiteration"])))+read_text_file(os.path.join(Tracker["directory"],"chunk_1_%03d.txt"%(Tracker["mainiteration"])))
			ctfs   = EMUtil.get_all_attributes(Tracker["constants"]["stack"],'ctf')
			ctfs   = [ctfs[i] for i in li]
			particle_groups = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_0.txt") ) + read_text_file(os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_1.txt") )
			npart = 500/Blockdata["nproc"] + 1
			li = list(range(len(ctfs)))
			shuffle(li)
			li = li[:npart*Blockdata["nproc"]]
			params = [params[i] for i in li]
			ctfs = [[ctfs[i].defocus, ctfs[i].cs, ctfs[i].voltage, ctfs[i].apix, ctfs[i].bfactor, ctfs[i].ampcont, ctfs[i].dfdiff, ctfs[i].dfang] for i in li]
			particle_groups = [particle_groups[i] for i in li]
		else:
			params = 0
			outlier_params = 0
			ctfs = 0
			particle_groups = 0
		params = wrap_mpi_bcast(params, Blockdata["main_node"])
		ctfs   = wrap_mpi_bcast(ctfs,   Blockdata["main_node"])
		particle_groups = wrap_mpi_bcast(particle_groups, Blockdata["main_node"])
		outlier_params = wrap_mpi_bcast(outlier_params, Blockdata["main_node"])
		#print(" A ",Blockdata["myid"] ,len(params),len(ctfs),len(particle_groups),len(params)/Blockdata["nproc"])
		npart = len(params)/Blockdata["nproc"]
		params = params[Blockdata["myid"]*npart:(Blockdata["myid"]+1)*npart]
		outlier_params = outlier_params[Blockdata["myid"]*npart:(Blockdata["myid"]+1)*npart]
		ctfs = [generate_ctf(ctfs[i]) for i in range(Blockdata["myid"]*npart,(Blockdata["myid"]+1)*npart)]
		particle_groups = particle_groups[Blockdata["myid"]*npart:(Blockdata["myid"]+1)*npart]
		Tracker["refvol"] = os.path.join(Tracker["directory"], "vol_0_%03d.hdf"%(Tracker["mainiteration"]))
		#print(" B ",Blockdata["myid"] ,len(params),len(ctfs),len(particle_groups),npart)
		params_outlier = []
		ctfs_outlier = []
		groups_outlier = []
		for idx, entry in enumerate(outlier_params):
			if entry == 0:
				params_outlier.append(params[idx])
				ctfs_outlier.append(ctfs[idx])
				groups_outlier.append(particle_groups[idx])

		cerrs(params_outlier, ctfs_outlier, groups_outlier)
		del params, ctfs, particle_groups
		del params_outlier, ctfs_outlier, groups_outlier
		if(Blockdata["myid"] == Blockdata["main_node"]):
			write_text_row( [[Tracker["acc_rot"], Tracker["acc_trans"]]], os.path.join(Tracker["directory"] ,"accuracy_%03d.txt"%(Tracker["mainiteration"])) )

	if(Blockdata["myid"] == Blockdata["main_node"]):
		anger, shifter = params_changes( read_text_row(os.path.join(Tracker["directory"],"params_%03d.txt"%(Tracker["mainiteration"]))), read_text_row(os.path.join(Tracker["previousoutputdir"],"params_%03d.txt"%(Tracker["mainiteration"]-1))) )
		write_text_row( [[anger, shifter]], os.path.join(Tracker["directory"] ,"error_thresholds_%03d.txt"%(Tracker["mainiteration"])) )

		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		line = ''
		sxprint(line,"Average displacements for angular directions = %6.2f degrees; and shifts = %5.1f pixels"%(anger, shifter) )

		#  Write current Trucker

		if  Blockdata["bckgnoise"] :
			Blockdata["bckgnoise"] = "computed"
		dump_tracker_to_json(
			os.path.join(
				Tracker["constants"]["masterdir"],
				"main%03d"%Tracker["mainiteration"],
				"Tracker_%03d.json"%Tracker["mainiteration"]
				),
			Tracker
			)
	Tracker["previousoutputdir"] = Tracker["directory"]

	return # parameters are all passed by Tracker


def prior_stack_fmt(sphire_prior_stack):
	fmt = []
	for entry in sphire_prior_stack.dtype.names:
		if isinstance(sphire_prior_stack[entry][0], basestring):
			fmt.append('% {0}s'.format(np.max([len(str_entry) for str_entry in sphire_prior_stack[entry]])+3))
		elif isinstance(sphire_prior_stack[entry][0], int):
			fmt.append('% {0}d'.format(len(str(np.max(sphire_prior_stack[entry])))+3))
		elif isinstance(sphire_prior_stack[entry][0], float):
			fmt.append('% {0}.6f'.format(len(str(np.max(sphire_prior_stack[entry])))+3))
		else:
			sxprint('UNKNOWN', entry, type(entry))
	return ' '.join(fmt)


def load_tracker_from_json(file_name):
	tracker = load_object_from_json(file_name)
	try:
		if tracker['constants']['stack_prior'] is not None:
			tracker['constants']['stack_prior'] = np.genfromtxt(
				tracker['constants']['stack_prior'],
				dtype=tracker['constants']['stack_prior_dtype']
				)
	except KeyError:
		tracker['constants']['stack_prior'] = None
	return tracker

def load_object_from_json(file_name):
	with open(file_name, 'r') as fin:
		json_object = convert_json_fromunicode(json.load(fin))
	return json_object

def dump_tracker_to_json(file_name, tracker):
	if tracker['constants']['stack_prior'] is not None:
		np.savetxt(
			os.path.join(tracker['constants']['masterdir'], 'stack_prior.txt'),
			Tracker['constants']['stack_prior'],
			fmt=Tracker['constants']['stack_prior_fmt']
			)
		tracker['constants']['stack_prior'] = os.path.join(tracker['constants']['masterdir'], 'stack_prior.txt')
	dump_object_to_json(file_name, tracker)

def dump_object_to_json(file_name, data_object):
	with open(file_name, 'w') as fout:
		json.dump(data_object, fout, indent=4)

# 		
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset 
#   (i.e. particle stack, reference structure, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#
#   This can be used to restart process from an arbitrary iteration.
#   
#


def get_image_statistics(image, mask, invert):
	if Tracker['constants']['filament_width'] is None:
		mask2d = mask
	else:
		mask2d = sp_utilities.model_rotated_rectangle2D(
			radius_long=int(np.sqrt(2 * image.get_xsize()**2) // 2),
			radius_short=int(Tracker['constants']['filament_width'] * image.get_xsize() / float(Tracker['constants']['nnxo']) + 0.5),
			nx=image.get_xsize(),
			ny=image.get_ysize(),
			angle=image.get_attr('segment_angle'),
			)
	return Util.infomask(image, mask2d, invert)


def calculate_prior_values(tracker, blockdata, outlier_file, chunk_file, params_file, im_start, im_end, procid):
	"""Calculate the prior values and identify outliers"""

	if not tracker['constants']['group_id']:
		return [0] * (im_end - im_start), None


	# Print to screen
	if(blockdata["myid"] == blockdata["main_node"]):
		sxprint("Executed successfully: ", "Prior calculation")

	# Calculate outliers
	if(blockdata["myid"] == blockdata["main_node"] ):
		# Calculate priors
		outliers, new_param, new_index = sp_helix_fundamentals.calculate_priors(
			tracker=Tracker,
			params_file=params_file,
			index_file=chunk_file,
			group_id=Tracker['constants']['group_id'],
			typ='sphire',
			tol_psi=Tracker['prior']['tol_psi'],
			tol_theta=Tracker['prior']['tol_theat'],
			tol_filament=Tracker['prior']['tol_filament'],
			tol_std=Tracker['prior']['tol_std'],
			tol_mean=Tracker['prior']['tol_mean'],
			outlier_method=Tracker['prior']['outlier_method'],
			prior_method=Tracker['prior']['prior_method'],
			force_outlier=Tracker['prior']['force_outlier'],
			window_size=Tracker['prior']['window_size'],
			remove_outlier=Tracker['prior']['remove_outlier'],
			symclass=Blockdata['symclass'],
			)

		# Print to screen
		nr_outliers = len(outliers[outliers == 1])
		len_data = len(outliers)
		no_outliers = len_data - nr_outliers
		sxprint('Chunk {0}: Discarded {1}|{2:.1f}%; Kept {3}|{4:.1f}%; Nr. particles {5}'.format(
			procid,
			nr_outliers,
			100*nr_outliers/float(len_data),
			no_outliers,
			100*no_outliers/float(len_data),
			len_data
			))
		if Tracker['prior']['force_outlier'] is not None:
			shutil.copy(params_file, '{0}_old'.format(params_file))
			shutil.copy(new_params, params_file)
		else:
			outliers = [0] * len_data

		np.savetxt(outlier_file, outliers)
	else:
		# Dummy variable
		outliers = 0

	# Distribute outlier list to all processes
	outliers = bcast_list_to_all(outliers, blockdata["myid"], blockdata["nodes"][0])

	# Get the node specific outlier information
	outliers_node = outliers[im_start:im_end]

	return outliers_node, outliers


def reduce_shifts(sx, sy, img):
	def rot_matrix(angle):
		angle = np.radians(angle)
		matrix = np.array([
			[np.cos(angle), -np.sin(angle)],
			[np.sin(angle), np.cos(angle)]
			])
		return matrix

	if Tracker['constants']['helical_rise'] is not None:
		try:
			rotation_angle = img.get_attr('segment_angle')
		except AttributeError:
			pass
		else:
			rise = Tracker['constants']['helical_rise'] / float(Tracker['constant']['pixel_size'])
			rise_half = rise / 2.0
			point = np.array([sx, sy])
			rot_point = np.dot(rot_matrix(rotation_angle), point.T)
			rot_point[0] = ((rot_point[0] + rise_half) % rise ) - rise_half
			sx, sy = np.dot(rot_matrix(rotation_angle).T, rot_point.T)

	return int(round(sx)), int(round(sy))


def main():

	from sp_utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	from sp_applications  import MPI_start_end
	from optparse      import OptionParser
	from sp_global_def    import SPARXVERSION
	from EMAN2         import EMData
	from sp_multi_shc     import multi_shc
	from random        import random, uniform
	from sp_logger        import Logger, BaseLogger_Files
	import sys
	import os
	import socket
	global Tracker, Blockdata

	# ------------------------------------------------------------------------------------
	# PARSE COMMAND OPTIONS
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ stack  [output_directory]  initial_volume  --radius=particle_radius --symmetry=c1 --initialshifts --inires=25  --mask3D=surface_mask.hdf --function=user_function
	
	There are five ways to run the program:

	1. Standard default run, starts from exhaustive searches, uses initial reference structure
	mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  bdb:sparx_stack vton1 mask15.hdf --sym=c5  --initialshifts  --radius=120  --mask3D=mask15.hdf    >1ovotn &

	2. Restart after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
	mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  vton1 --radius=100 >2ovotn &

	3. Local refinement, starts from user-provided orientation parameters, delta has to be <= 3.75
	mpirun -np 64 --hostfile four_nodes.txt sxmeridien.py --local_refinement bdb:sparx_stack   vton3 --delta=1.875 --xr=2.0  --inires=5.5  --sym=c5  --radius=120  --mask3D=mask15.hdf >5ovotn &

	4. Restart of local refinement after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
	mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py --local_refinement  vton3  --xr=0.6 >6ovotn &

	5.  mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py vton3 --do_final=21
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--do_final",           type="int",			 	default= -1,  help="Do unfiltered odd and even volume 3-D reconstruction from an existing meridien refinement with optional specified iteration")
	parser.add_option("--local_refinement",    action="store_true",  default= False,  help="Perform local refinement starting from user-provided orientation parameters")
	parser.add_option("--memory_per_node",          type="float",           default= -1.0,                	help="User provided information about memory per node (NOT per CPU) [in GB] (default 2GB*(number of CPUs per node))")	


	do_final_mode = False
	for q in sys.argv[1:]:
		if( q[:10] == "--do_final" ):
			do_final_mode = True
			break


	do_continuation_mode = False
	for q in sys.argv[1:]:
		if( q[:18] == "--local_refinement" ):
			do_continuation_mode = True
			break

	if not do_final_mode:

		parser.add_option("--radius",      		   		type= "int",          	default= -1,			     	help="Outer radius [in pixels] of particles < int(nx/2)-1")
		parser.add_option("--xr",      		       		type="float",         	default= 5.,		         	help="Range for translation search in both directions, search is +/xr (default 5), can be fractional")
		parser.add_option("--ts",      		       		type="float",        	default= 1.,		         	help="Step size of the translation search in both directions, search is within a circle of radius xr on a grid with steps ts, (default 1), can be fractional")
		parser.add_option("--mask3D",		        	type="string",	      	default= None,		          	help="3D mask file (default a sphere with radius (nx/2)-1)")
		parser.add_option("--function",					type="string",          default= "do_volume_mask",		help="Name of the reference preparation function (default do_volume_mask)")
		parser.add_option("--function_ai",				type="string",          default= "ai_spa",		help="Name of the internal heuristic function (default ai_spa)")
		parser.add_option("--symmetry",					type="string",        	default= 'c1',		     		help="Point-group symmetry of the refined structure (default c1)")
		parser.add_option("--an",	           		    type="float", 	     	default=-1.,                	help="Angular neighborhood for local search")
		parser.add_option("--shake",	           		type="float", 	     	default=0.5,                	help="Shake (0.5)")
		parser.add_option("--limit_improvement",	           		type="int", 	     	default=1,                	help="No improvement limit for the convergence criterion (Default 1)")
		parser.add_option("--a_criterion",	           		type="float", 	     	default=0.75,                	help="A criterions multiplication to acc_rot (Default 0.75)")
		parser.add_option("--small_memory",         	action="store_true",  	default= False,             	help="Data will not be kept in memory if small_memory is true. (default False)")
		parser.add_option("--ccfpercentage",			type="float", 	      	default= 99.9,               	help="Percentage of the correlation peak area to be included, 0.0 corresponds to hard matching (default 99.9%)")
		parser.add_option("--nonorm",               	action="store_true",  	default= False,              	help="Do not apply image norm correction. (default False)")
		parser.add_option("--group_by",               	type="str",  	default= 'ptcl_source_image',              	help="Group particles by header information. For helical refinement use filament or filament_id if present. (Default ptcl_source_image)")
		parser.add_option("--theta_min",               	type="float",  	default= -1,              	help="Lower limit for the out-of-plane rotation angle. Default is the full range based on the symmetry. (Default -1)")
		parser.add_option("--theta_max",               	type="float",  	default= -1,              	help="Upper limit for the out-of-plane rotation angle.  Default is the full range based on the symmetry. (Default -1)")
		parser.add_option("--even_angle_method",               	type="str",  	default='S',              	help="Even angle creation strategy. Choices: S, P, M. (Default S)")
		parser.add_option("--group_id",               	type="str",  	default=None,              	help="Group particles for outlier detection by header name. Useful ones are e.g. ISAC_class_id or filament_id (Default None)")
		parser.add_option("--filament_width",         	type="int",  	default=None,              	help="Filament width used to normalize the particles. (Default None)")
		parser.add_option("--helical_rise",         	type="float",  	default=None,              	help="Helical rise in angstrom. This is used to limit the shift along the helical axis. (Default None)")
		parser.add_option("--plot_ang_dist",         	action='store_true',  	default=False,              	help="Plot the angular distribution plot for every iteration. This will take some time for high symmetries. (Default False)")
		if do_continuation_mode:
			# case1: local meridien run using parameters stored in headers
			# case2: restart mode of standard meridien run. Parameters can be altered in the restart run.
			parser.add_option("--inires",		       		type="float",	     	default= -1,		         	help="Resolution of the initial_volume volume. One can use -1 in local_refinement mode to filter to the resolution of the reconstructed volumes. (default -1A)")
			parser.add_option("--delta",					type="float",			default=3.75,		     		help="Initial angular sampling step (default 3.75)")
			parser.add_option("--main000",					type="str",			default=None,		     		help="Main000 folder of a previous refinement to presever chunk and group information. By default it will assign new chunks and groups. (default None)")
			(options, args) = parser.parse_args(sys.argv[1:])

			if( len(args) == 2 ):
				masterdir 	= args[1]
				orgstack 	= args[0]
				volinit = ''

			elif( len(args) == 1 ):
				masterdir = args[0]
				volinit = ''
				orgstack = ''

			else:
				sxprint( "usage: " + usage)
				sxprint( "Please run '" + progname + " -h' for detailed options")
				ERROR( "Invalid number of parameters used. Please see usage information above.", myid=Blockdata['myid'])
				return

			if options.delta > 3.75:
				ERROR( "Local searches requested, delta cannot be larger than 3.75.", myid=Blockdata["myid"] )
				return

			setattr(options, 'initialshifts', True)
			setattr(options, 'skip_prealignment', True)
			setattr(options, 'center_method', -1)
			setattr(options, 'target_radius', 29)

		else:
			# case1: standard meridien run
			# case2: restart mode of standard meridien run. Parameters can be altered in the restart run.
			parser.add_option("--inires",		       		type="float",	     	default= 25,		         	help="Resolution of the initial_volume volume. One can use -1 in local_refinement mode to filter to the resolution of the reconstructed volumes. (default 25)")
			parser.add_option("--delta",					type="float",			default=7.5,		     		help="Initial angular sampling step (default 7.5)")
			parser.add_option("--skip_prealignment",		action="store_true", 	default= False,		         	help="Skip 2-D pre-alignment step: to be used if images are already centered. (default False)")
			parser.add_option("--initialshifts",         	action="store_true",  	default= False,	         		help="Use orientation parameters in the input file header to jumpstart the procedure. (default False)")
			parser.add_option("--center_method",			type="int",			 	default= -1,			     	help="Method for centering: of average during initial 2D prealignment of data (0 : no centering; -1 : average shift  method;  please see center_2D in utilities.py for methods 1-7) (default -1)")
			parser.add_option("--target_radius", 			type="int",			 	default= 29,			     	help="Target particle radius for 2D prealignment. Images will be shrank/enlarged to this radius (default 29)")
			(options, args) = parser.parse_args(sys.argv[1:])

			setattr(options, 'main000', None)

			if( len(args) == 3 ):
				volinit 	= args[2]
				masterdir 	= args[1]
				orgstack 	= args[0]

			elif( len(args) == 2 ):
				orgstack 	= args[0]
				volinit 	= args[1]
				masterdir = ""

			elif( len(args) == 1 ): masterdir = args[0]
			else:
				sxprint( "usage: " + usage)
				sxprint( "Please run '" + progname + " -h' for detailed options")
				ERROR( "Invalid number of parameters used. Please see usage information above.", myid=Blockdata['myid'])
				return

		if Blockdata["myid"]  == Blockdata["main_node"]:
			line = ""
			for a in sys.argv:
				line +=a+"  "
			sxprint(" shell line command ")
			sxprint(line)


		#  Check whether we are restarting the program, in the least main000 should exist, otherwise there is nothing to restart
		keepgoing1   = 1
		keepgoing2   = 1
		restart_flag = 0
		if Blockdata["myid"]  == Blockdata["main_node"]:
			if( os.path.exists(os.path.join(masterdir, "main000", "Tracker_000.json")) ):
				if(len(args) > 1): keepgoing1 = 0
				restart_flag = 1
			else:
				if(len(args) == 1): keepgoing2  = 0
				restart_flag  = 0
		restart_flag = bcast_number_to_all(restart_flag, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
		keepgoing1   = bcast_number_to_all(keepgoing1,   source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
		keepgoing2   = bcast_number_to_all(keepgoing2,   source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)

		if keepgoing1== 0:
			ERROR( "To restart, meridien requires only the name of existing refinement directory.", myid=Blockdata["myid"] )
			return

		if keepgoing2 ==0:
			ERROR( "To start, meridien requires at least the stack name and the name of reference structure", myid=Blockdata["myid"] )
			return

		if restart_flag ==1:
			restart_mode = True
		else: 
			restart_mode  = False
		
		# ------------------------------------------------------------------------------------
		# Initialize MPI related variables

		# ------------------------------------------------------------------------------------
		#  INPUT PARAMETERS
		#  For the time being we use all CPUs during refinement
		Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
		Blockdata["nsubset"]      = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		create_subgroup()
		create_zero_group()

		Blockdata["rkeepf"] = 0.90

		if not restart_mode: #<<<-------Fresh run
			#  Constant settings of the project
			Constants		       			        = {}
			Constants["stack"]             			= args[0]
			Constants["rs"]                			= 1
			Constants["radius"]            			= options.radius
			Constants["an"]                			= "-1"
			Constants["maxit"]             			= 1
			Constants["fuse_freq"]         			= 45  # Now in A, convert to absolute before using
			Constants["symmetry"]					= options.symmetry.lower()
			Constants["npad"]              			= 1
			Constants["center"]            			= 0
			Constants["shake"]             			= options.shake  #  move params every iteration
			Constants["CTF"]               			= True # internally set
			Constants["mask3D"]            			= options.mask3D
			Constants["nnxo"]              			= -1
			Constants["pixel_size"]        			= None # read from data
			Constants["inires"]            			= options.inires  # Now in A, convert to absolute before using
			Constants["refvol"]            			= volinit
			Constants["masterdir"]         			= masterdir
			Constants["best"]              			= 3
			Constants["ccfpercentage"]     			= options.ccfpercentage/100. # Here for legacy reasons
			Constants["limit_improvement"] 			= options.limit_improvement
			Constants["a_criterion"] 			    = options.a_criterion
			Constants["limit_changes"]     			= 1  # reduce delta by half if both limits are reached simultaneously
			Constants["states"]            			= ["INITIAL", "PRIMARY", "EXHAUSTIVE", "RESTRICTED", "PRIMARY LOCAL", "FINAL"]# will add two states, CONINUATION_INITIAL, CONINUATION_PRIMARY
			Constants["user_func_volume"]			= options.function
			Constants["user_func_ai"]				= options.function_ai
			Constants["hardmask"]          			=  True #options.hardmask
			Constants["expthreshold"]      			= -10
			Constants["number_of_groups"]  			= -1 # number of defocus groups, to be set by assign_particles_to_groups
			Constants["nonorm"]            			= options.nonorm
			Constants["small_memory"]      			= options.small_memory
			Constants["memory_per_node"] 			= options.memory_per_node
			Constants["initialshifts"]			    = options.initialshifts
			Constants["even_angle_method"]			    = options.even_angle_method
			Constants["group_id"]			    = options.group_id
			Constants["filament_width"]			    = options.filament_width
			Constants["helical_rise"]			    = options.helical_rise
			Constants["do_local"]			    = do_continuation_mode
			Constants["plot_ang_dist"]			    = options.plot_ang_dist
			if options.group_id is None:
				Constants['stack_prior'] = None
				Constants['stack_prior_fmt'] = None
				Constants['stack_prior_dtype'] = None
			else:
				Constants['stack_prior'] = sp_helix_sphire.import_sphire_stack(args[0], options.group_id)
				Constants['stack_prior_fmt'] = prior_stack_fmt(Constants['stack_prior'])
				Constants['stack_prior_dtype'] = Constants['stack_prior'].dtype.descr

			Prior = {}
			Prior['tol_psi'] = 30
			Prior['tol_theta'] = 15
			Prior['tol_filament'] = 0.2
			Prior['tol_std'] = 1
			Prior['tol_mean'] = 30
			Prior['outlier_method'] = 'deg'
			Prior['prior_method'] = 'running'
			Prior['force_outlier'] = False
			Prior['remove_outlier'] = False
			Prior['window_size'] = 3

			#
			#  The program will use three different meanings of x-size
			#  nnxo         - original nx of the data, will not be changed
			#  nxinit       - window size used by the program during given iteration, 
			#                 will be increased in steps of 32 with the resolution
			#
			#  nxstep       - step by wich window size increases
			#
			# Initialize Tracker Dictionary with input options
			Tracker = {}
			Tracker["constants"]	= Constants
			Tracker["prior"]		= Prior
			Tracker["maxit"]		= Tracker["constants"]["maxit"]
		
			Tracker["xr"]			= options.xr
			Tracker["yr"]			= options.xr  # Do not change!  I do not think it is used anywhere
			Tracker["ts"]			= options.ts
			Tracker["an"]			= "-1"
			Tracker["theta_min"]	= options.theta_min
			Tracker["theta_max"]	= options.theta_max
			Tracker["delta"]		= options.delta  # How to decide it
			Tracker["refvol"]		= None
			Tracker["nxinit"]		= -1  # will be figured in first AI.
			Tracker["nxstep"]		= 10
			#  Resolution in pixels at 0.5 cutoff
			Tracker["currentres"]		    = -1
			Tracker["ccfpercentage"]     			= options.ccfpercentage/100.
			Tracker["howmany"]		        = 4
			Tracker["fsc143"]			    = -1
			Tracker["maxfrad"]           	= -1
			Tracker["no_improvement"]    	= 0
			Tracker["no_params_changes"] 	= 0
			Tracker["large_at_Nyquist"]  	= False
			Tracker["anger"]             	= 1.e23
			Tracker["shifter"]           	= 1.e23
			Tracker["pixercutoff"]       	= 2.0
			Tracker["directory"]         	= ""
			Tracker["previousoutputdir"] 	= ""
			Tracker["acc_rot"]           	= 0.0
			Tracker["acc_trans"]			= 0.0
			Tracker["avgvaradj"]			= [1.0,1.0]  # This has to be initialized to 1.0 !!
			Tracker["mainiteration"]     	= 0
			Tracker["lentop"]				= 2000
			Tracker["state"]             	= Tracker["constants"]["states"][0]
			Tracker["nima_per_chunk"]    	= [0,0]
			###<<<----state 
			Tracker["bestres"]          	= 0
			Tracker["bestres_143"]          = 0
			Tracker["changed_delta"]        = False
			Blockdata["bckgnoise"]          = None
			Blockdata["accumulatepw"]       = [[],[]]

			# ------------------------------------------------------------------------------------
			# Get the pixel size; if none, set to 1.0, and the original image size
			if(Blockdata["myid"] == Blockdata["main_node"]):
				#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				line = ''
				sxprint(line,"INITIALIZATION OF MERIDIEN")
				a = get_im(orgstack)
				nnxo = a.get_xsize()
				if Tracker["constants"]["CTF"]:
					i = a.get_attr('ctf')
					pixel_size = i.apix
					fq = int(pixel_size*nnxo/Tracker["constants"]["fuse_freq"] + 0.5)
				else:
					pixel_size = Tracker["constants"]["pixel_size"]
					#  No pixel size, fusing computed as 5 Fourier pixels
					fq = 5
				del a
			else:
				nnxo = 0
				fq = 0
				pixel_size = 1.0
			
			#  Object to handle symmetries, for now only used by oct
			from EMAN2 import parsesym
			#Blockdata["parsesym"] = parsesym(Tracker["constants"]["symmetry"])
			#  Initialize symmetry
			Blockdata["symclass"] = symclass(Tracker["constants"]["symmetry"])
			Blockdata["symclass_coarse"] = symclass(Tracker["constants"]["symmetry"])

			nnxo = bcast_number_to_all(nnxo, source_node = Blockdata["main_node"])

			if( nnxo < 0 ):
				ERROR( "Incorrect image size", myid=Blockdata["myid"] )
				return
			
			pixel_size = bcast_number_to_all(pixel_size, source_node = Blockdata["main_node"])
			fq         = bcast_number_to_all(fq,         source_node = Blockdata["main_node"])
			Tracker["constants"]["nnxo"]         = nnxo
			Tracker["constants"]["pixel_size"]   = pixel_size
			Tracker["constants"]["fuse_freq"]    = fq
			del fq, nnxo, pixel_size
			if(Tracker["constants"]["inires"]>0.0):
				Tracker["constants"]["inires"] = int(
					Tracker["constants"]["nnxo"] * 
					Tracker["constants"]["pixel_size"] / Tracker["constants"]["inires"] +
					0.5
					)
			Tracker["currentres"] = Tracker["constants"]["inires"]
			Tracker["fsc143"]     = Tracker["constants"]["inires"]

			###  VARIOUS SANITY CHECKS <-----------------------
			if options.memory_per_node < 0.0:
				options.memory_per_node = 2.0 * Blockdata["no_of_processes_per_group"]

			if options.initialshifts and not options.skip_prealignment:
				ERROR('Initialshifts and skip_prealignment option provided. Enable skip_prealignment for further processing.', action=0, myid=Blockdata['myid'])
				options.skip_prealignment = True

			checking_flag = 1
			if(Blockdata["myid"] == Blockdata["main_node"]):
				if( Tracker["constants"]["mask3D"] and (not os.path.exists(Tracker["constants"]["mask3D"]))):
					checking_flag = 0

			checking_flag = bcast_number_to_all(checking_flag, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD)
			if checking_flag==0:
				ERROR( "Mask3D file does  not exists ", myid=Blockdata["myid"] )
				return
			
			if options.ts > 0:
				if options.xr / options.ts < 1.0: 
					ERROR( "Incorrect translational searching settings, search range cannot be smaller than translation step ", myid=Blockdata["myid"] )
					return
			else:
				ERROR( "Incorrect translational searching settings, translational search range cannot be smaller equals 0", myid=Blockdata["myid"] )
				return

			if( 2*(Tracker["currentres"] + Tracker["nxstep"]) > Tracker["constants"]["nnxo"] ):
				ERROR( "Image size less than what would follow from the initial resolution provided %d  %d  %d" % (Tracker["currentres"], Tracker["nxstep"], 2*(Tracker["currentres"] + Tracker["nxstep"])), myid=Blockdata["myid"] )
				return

			if(Tracker["constants"]["radius"]  < 1):
				Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2

			elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
				ERROR( "Particle radius set too large", myid=Blockdata["myid"] )
				return

			if not balanced_processor_load_on_nodes:
				ERROR( "Nodes do not have the same number of CPUs, please check configuration of the cluster.", myid=Blockdata["myid"] )
				return


			#  MASTER DIRECTORY
			if(Blockdata["myid"] == Blockdata["main_node"]):
				if( masterdir == ""):
					timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
					masterdir = "master"+timestring
					li = len(masterdir)
					os.makedirs(masterdir)
					keepchecking = 0
				else:
					if not os.path.exists(masterdir):
						os.makedirs(masterdir)
					li = 0
					keepchecking = 1
				sp_global_def.write_command(masterdir)
			else:
				li = 0
				keepchecking = 1

			li = mpi_bcast(li,1,MPI_INT,Blockdata["main_node"],MPI_COMM_WORLD)[0]

			if( li > 0 ):
				masterdir = mpi_bcast(masterdir,li,MPI_CHAR,Blockdata["main_node"],MPI_COMM_WORLD)
				masterdir = string.join(masterdir,"")

			Tracker["constants"]["masterdir"] = masterdir
			initdir = os.path.join(Tracker["constants"]["masterdir"],"main000")
			if(Blockdata["myid"] == Blockdata["main_node"]):
				if os.path.exists(initdir): shutil.rmtree(initdir)
				os.mkdir(initdir)

			# Initialization of orgstack
			Tracker["constants"]["stack"] = orgstack 
			if(Blockdata["myid"] == Blockdata["main_node"]):
				total_stack = EMUtil.get_image_count(Tracker["constants"]["stack"])
			else:
				total_stack = 0
			total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])

			partids = os.path.join(initdir, "indexes_000.txt")
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				write_text_file(list(range(total_stack)), partids)
			mpi_barrier(MPI_COMM_WORLD)


			init2dir = os.path.join(Tracker["constants"]["masterdir"], "2dalignment")
			if Blockdata['myid'] == Blockdata['main_node']:
				sxprint('2D pre-alignment step')

			kwargs = dict()

			kwargs["init2dir"]  							= init2dir
			kwargs["myid"]      							= Blockdata["myid"]
			kwargs["main_node"] 							= Blockdata["main_node"]
			kwargs["number_of_images_in_stack"] 			= total_stack
			kwargs["nproc"] 								= Blockdata["nproc"]

			kwargs["target_radius"] 						= options.target_radius
			# kwargs["target_nx"] = target_nx
			kwargs["radi"] 									= options.radius

			kwargs["center_method"] 						= options.center_method

			kwargs["nxrsteps"] 								= 4

			kwargs["command_line_provided_stack_filename"] 	= Tracker["constants"]["stack"]

			# kwargs["masterdir"] = masterdir

			kwargs["options_skip_prealignment"] 			= options.skip_prealignment 
			kwargs["options_CTF"] 							= True

			kwargs["mpi_comm"] 								= MPI_COMM_WORLD
			params2d = calculate_2d_params_for_centering(kwargs)
			del kwargs
			if Blockdata['myid'] == Blockdata['main_node']:
				sxprint('2D pre-alignment completed')

			#  store params
			partids = [None]*2
			for procid in range(2):
				partids[procid] = os.path.join(initdir,"chunk_%01d_000.txt"%procid)
			partstack = [None]*2
			for procid in range(2):
				partstack[procid] = os.path.join(initdir,"params-chunk_%01d_000.txt"%procid)

			if(Blockdata["myid"] == Blockdata["main_node"]):
				if options.main000 is None:
					sxprint('Assign particles to groups')
					l1, l2 = assign_particles_to_groups(minimum_group_size = 10, name_tag=options.group_by)
					write_text_file(l1,partids[0])
					write_text_file(l2,partids[1])
				else:
					sxprint('Copy chunk and group information from existing directory.')
					names = (
						'chunk_0_000.txt',
						'chunk_1_000.txt',
						'groupids_000.txt',
						'groupids_001.txt',
						'groupids.txt',
						'indexes_000.txt',
						'micids.txt',
						'micids_000.txt',
						'micids_001.txt',
						'number_of_particles_per_group.txt',
						'particle_groups_0.txt',
						'particle_groups_1.txt',
						)
					for name_entry in names:
						shutil.copy2(
							os.path.join(options.main000, name_entry),
							os.path.join(Tracker['constants']['masterdir'], 'main000', name_entry)
							)
					l1 = read_text_file(partids[0])
					l2 = read_text_file(partids[1])
					Tracker['nima_per_chunk'] = [len(l1), len(l2)]
					Tracker['constants']['number_of_groups'] = len(read_text_row(os.path.join(Tracker['constants']['masterdir'], 'main000', 'groupids.txt')))
				if options.initialshifts: # Always True for continue mode as initialised in the option parser
					tp_list = EMUtil.get_all_attributes(Tracker["constants"]["stack"], "xform.projection")
					for i in range(len(tp_list)):
						dp = tp_list[i].get_params("spider")
						tp_list[i] = [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], 0.0, 1.0, 1.0]
					write_text_row(tp_list, os.path.join(initdir,"params_000.txt"))
					write_text_row([tp_list[i] for i in l1], partstack[0])
					write_text_row([tp_list[i] for i in l2], partstack[1])
					del tp_list
					sxprint('Executed successfully: Imported initial parameters from the input stack.')
				else:
					write_text_row([[0,0,0,params2d[i][1],params2d[i][2], 0.0, 1.0] for i in l1], partstack[0])
					write_text_row([[0,0,0,params2d[i][1],params2d[i][2], 0.0, 1.0] for i in l2], partstack[1])
					write_text_row([[0,0,0,params2d[i][1],params2d[i][2], 0.0, 1.0] for i in range(len(l1)+len(l2))], os.path.join(initdir,"params_000.txt"))
				del l1, l2

			Tracker["nima_per_chunk"][0] = bcast_number_to_all(Tracker["nima_per_chunk"][0], Blockdata["main_node"])
			Tracker["nima_per_chunk"][1] = bcast_number_to_all(Tracker["nima_per_chunk"][1], Blockdata["main_node"])
			Tracker["constants"]["number_of_groups"] = bcast_number_to_all(Tracker["constants"]["number_of_groups"], Blockdata["main_node"])
			del params2d

			if Tracker['constants']['do_local']:
				if( Tracker["constants"]["inires"] > 0 ):
					Tracker["nxinit"] = min(2*Tracker["constants"]["inires"], Tracker["constants"]["nnxo"] )
				else: 
					Tracker["nxinit"] = Tracker["constants"]["nnxo"]
					
				rec3d_continuation_nosmearing(MPI_COMM_WORLD)

			elif Blockdata['myid'] == Blockdata['main_node']:
				# Create reference models for each particle group
				sxprint('Prepare reference volume')
				if(Tracker["constants"]["mask3D"] == None):
					viv = filt_table(cosinemask(get_im(volinit),radius = Tracker["constants"]["radius"]), [1.0]*Tracker["constants"]["inires"] + [0.5] + [0.0]*Tracker["constants"]["nnxo"])
				else:
					viv = filt_table(get_im(volinit)*get_im(Tracker["constants"]["mask3D"]), [1.0]*Tracker["constants"]["inires"] + [0.5] + [0.0]*Tracker["constants"]["nnxo"])
				# make a copy of original reference model for this particle group (procid)
				for procid in range(2):
					viv.write_image(os.path.join(initdir,"vol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
				del viv

			if(Blockdata["myid"] == Blockdata["main_node"]):
				dump_tracker_to_json(
					os.path.join(
						Tracker["constants"]["masterdir"],
						"main%03d"%Tracker["mainiteration"],
						"Tracker_%03d.json"%Tracker["mainiteration"]
						),
					Tracker
					)

		else:
			Blockdata["bckgnoise"] 		= None
			Blockdata["accumulatepw"] 	= [[],[]]
			initdir 			= os.path.join(masterdir,"main000")
			keepchecking 		= 1
			if(Blockdata["myid"] == Blockdata["main_node"]):
				Tracker = load_tracker_from_json(os.path.join(initdir,"Tracker_000.json"))
				print_dict(Tracker["constants"], "Permanent settings of the original run recovered from main000")
			else: Tracker = None
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
			Tracker["mainiteration"] = 0
		mpi_barrier(MPI_COMM_WORLD)


		projdata       = [[model_blank(1,1)], [model_blank(1,1)]]
		oldparams      = [[],[]]
		currentparams  = [[],[]]
		original_data  = [None, None]

		keepgoing = 1
		while(keepgoing):
			Tracker["previousoutputdir"] = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"])
			Tracker["mainiteration"] += 1
			Tracker["directory"]     = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"])
			doit, keepchecking       = checkstep(Tracker["directory"], keepchecking)
			if( not doit ):			
				li = True
				doit2, keepchecking2 = \
				    checkstep(os.path.join(Tracker["directory"],"Tracker_%03d.json"%Tracker["mainiteration"]), li)
				if doit2:
					doit = True
					if( Blockdata["myid"] == Blockdata["main_node"]): 
						shutil.rmtree(Tracker["directory"])
					mpi_barrier(MPI_COMM_WORLD)
			if doit :
				if(Blockdata["myid"] == Blockdata["main_node"]):
					mainiteration = Tracker["mainiteration"]
					Tracker = load_tracker_from_json(
						os.path.join(
							Tracker["previousoutputdir"],
							"Tracker_%03d.json"%(Tracker["mainiteration"]-1)
							)
						)
					#  It has to be repeated here as Tracker is from previous iteration, I see no other way.
					Tracker["previousoutputdir"]	= os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"])
					Tracker["mainiteration"]		= mainiteration
					Tracker["directory"]			= os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"])
				else: Tracker = None
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
				### 
				if restart_mode:
					update_tracker(sys.argv[1:])
					update_memory_estimation()
					restart_mode = False

				# prepare names of input file names, they are in main directory,
				#   log subdirectories contain outputs from specific refinements
				partids = [None]*2
				for procid in range(2):
					partids[procid] = os.path.join(Tracker["previousoutputdir"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]-1))
				partstack = [None]*2
				for procid in range(2):
					partstack[procid] = os.path.join(Tracker["previousoutputdir"],"params-chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]-1))

				mpi_barrier(MPI_COMM_WORLD)

				if Tracker["mainiteration"] == 1:
					if Tracker['constants']['do_local']:
						fff = read_text_file(os.path.join(Tracker["previousoutputdir"],"driver_%03d.txt"%(Tracker["mainiteration"]-1)))
					else:
						fff = []
					anger   = 1.0e9
					shifter = 1.0e9
				else:
					if(Blockdata["myid"] == Blockdata["main_node"]):
						fff = read_text_file(os.path.join(Tracker["previousoutputdir"],"driver_%03d.txt"%(Tracker["mainiteration"]-1)))
						[anger, shifter] = read_text_row( os.path.join(Tracker["previousoutputdir"] ,"error_thresholds_%03d.txt"%(Tracker["mainiteration"]-1)) )[0]
					else:
						fff = []
						anger   = 0.0
						shifter = 0.0
					fff     = bcast_list_to_all(fff, Blockdata["myid"], source_node = Blockdata["main_node"])
					anger   = bcast_number_to_all(anger,                source_node = Blockdata["main_node"])
					shifter = bcast_number_to_all(shifter,              source_node = Blockdata["main_node"])

				func_ai = sp_user_functions.factory[Tracker["constants"]["user_func_ai"]]
				keepgoing = func_ai(Tracker, fff, anger, shifter, Tracker['constants']['do_local'], Blockdata['myid'] == Blockdata['main_node'])

				if keepgoing == 1: # not converged
					if Blockdata["myid"] == Blockdata["main_node"]:
						if(Tracker["mainiteration"] > 1):
							#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
							line = ''
							sxprint(line,"Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."%
								(Tracker["mainiteration"]-1, \
								Tracker["currentres"], Tracker["fsc143"], Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(Tracker["currentres"]), \
								Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(Tracker["fsc143"]) ) )
						sxprint("\n\n\n\n")
						#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
						line = ''
						sxprint(line,"ITERATION  #%2d. Current state: %14s, nxinit: %3d, delta: %9.4f, xr: %9.4f, ts: %9.4f"%\
							(Tracker["mainiteration"], Tracker["state"],Tracker["nxinit"], Tracker["delta"], Tracker["xr"], Tracker["ts"]  ))
					#print("RACING  A ",Blockdata["myid"])
					li = True
					doit2, keepchecking2 = checkstep(Tracker["directory"], li)
					if(Blockdata["myid"] == Blockdata["main_node"] and doit2):
						cmd = "{} {}".format("mkdir", Tracker["directory"])
						junk = cmdexecute(cmd)
						cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"],"oldparamstructure"))
						junk = cmdexecute(cmd)
					if(not doit2):
						ERROR("There was a gap in main directories, program cannot proceed","sxmeridien",1,Blockdata["myid"])
					mpi_barrier(MPI_COMM_WORLD)

					#  READ DATA AND COMPUTE SIGMA2   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
					refinement_one_iteration(partids, partstack, original_data, oldparams, projdata, continuation_mode=Tracker['constants']['do_local'])
					#	sxprint("  MOVING  ON --------------------------------------------------------------------")
				else: # converged, do final
					if( Blockdata["subgroup_myid"] > -1 ):
						mpi_comm_free(Blockdata["subgroup_comm"])
						
					if(Blockdata["myid"] == Blockdata["main_node"]):
						sxprint(line,"Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."%
								(Tracker["mainiteration"]-1, \
								Tracker["currentres"], Tracker["fsc143"], Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(Tracker["currentres"]), \
								Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(Tracker["fsc143"]) ) )
						sxprint("\n\n\n\n")
						
						sxprint("Iteration that contains the best resolution is %d"%Tracker["constants"]["best"])
						if Tracker["constants"]["best"] ==2:
							sxprint("No resolution improvement in refinement ")
						
						dump_tracker_to_json(
							os.path.join(masterdir,"Tracker_final.json"),
							Tracker
							)
					mpi_barrier(MPI_COMM_WORLD)
					
					newparamstructure 			= [[],[]]
					projdata          			= [[model_blank(1,1)], [model_blank(1,1)]]
					original_data     			= [None,None]
					oldparams         			= [[],[]]
					Blockdata["accumulatepw"]  	= [None, None]
					#if Tracker["constants"]["memory_per_node"] <0.0: Tracker["constants"]["memory_per_node"] = 2.0*Blockdata["no_of_processes_per_group"]
					recons3d_final(Tracker["constants"]["masterdir"], Tracker["constants"]["best"], Tracker["constants"]["memory_per_node"])

			#  End of if doit
		#   end of while
		return

	elif do_final_mode: #  DO FINAL
		(options, args) = parser.parse_args(sys.argv[1:])
		#global Tracker, Blockdata
		#print( "  args  ",args)
		checking_flag = 1
		orgstack      = None
		if( len(args) == 3):
			ERROR( "do_final option requires only one or two arguments ", myid=Blockdata["myid"] )
			return

		elif(len(args) == 2): # option for signal subtraction 
			masterdir = args[1]
			orgstack  = args[0]
			
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(masterdir): 
					checking_flag = 0
			
			checking_flag = bcast_number_to_all(checking_flag, source_node = Blockdata["main_node"])
			
			if checking_flag ==0:  
				ERROR( "do_final: refinement directory for final reconstruction does not exist ", myid=Blockdata["myid"] )
				return
			
		elif(len(args) == 1):
			masterdir 	= args[0]

			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(masterdir): 
					checking_flag = 0
			checking_flag = bcast_number_to_all(checking_flag, source_node = Blockdata["main_node"])

			if checking_flag == 0: 
				ERROR( "do_final: refinement directory for final reconstruction does not exist ", myid=Blockdata["myid"] )
				return
			
		else:
			sxprint( "Usage: " + usage )
			sxprint( "Please run \'" + progname + " -h\' for detailed options" )
			ERROR( "Invalid number of parameters used. Please see usage information above." )
			return

		if( options.do_final < 0):
			ERROR( "Incorrect iteration number in do_final  %d"%options.do_final, myid=Blockdata["myid"] )
			return
		#print(  orgstack,masterdir,volinit )
		# ------------------------------------------------------------------------------------
		# Initialize MPI related variables

		###print("  MPIINFO  ",Blockdata)
		###  MPI SANITY CHECKES
		if not balanced_processor_load_on_nodes: 
			ERROR( "Nodes do not have the same number of CPUs, please check configuration of the cluster.", myid=Blockdata["myid"] )
			return
		#if( Blockdata["no_of_groups"] < 2 ):  ERROR("To run, program requires cluster with at least two nodes.","meridien",1,Blockdata["myid"])
		###
		if Blockdata["myid"]  == Blockdata["main_node"]:
			line = ""
			for a in sys.argv:
				line +=a+"  "
			sxprint(" shell line command ")
			sxprint(line)
		# ------------------------------------------------------------------------------------
		#  INPUT PARAMETERS
		sp_global_def.BATCH = True
		sp_global_def.MPI   = True

		###  VARIOUS SANITY CHECKES <-----------------------
		if( options.memory_per_node < 0.0 ):
			options.memory_per_node = 2.0*Blockdata["no_of_processes_per_group"]
	
		Blockdata["accumulatepw"] = [[],[]]
		recons3d_final(masterdir, options.do_final, options.memory_per_node, orgstack)
		return
	else:
		ERROR( "Incorrect input options", myid=Blockdata["myid"] )
		return

if __name__=="__main__":
	sp_global_def.BATCH = True
	sp_global_def.MPI   = True
	sp_global_def.print_timestamp( "Start" )
	main()
	sp_global_def.print_timestamp( "Finish" )
	sp_global_def.BATCH = False
	sp_global_def.MPI   = False
	mpi.mpi_finalize()
