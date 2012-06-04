#!/usr/bin/env python
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

import	global_def
from	global_def 	import *
from	optparse 	import OptionParser
from	string 		import atoi,replace
from	EMAN2 		import EMUtil
import	os
import	sys
from 	time		import	time
from	utilities	import print_begin_msg, print_end_msg, print_msg
from	utilities	import read_text_row

t0 = time()

def main():

	def params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror):
		if mirror == False:
			m = 1
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 540.0-psi, 0, 0, 1.0)
		else:
			m = 0
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 360.0-psi, 0, 0, 1.0)
		return  alpha, sx, sy, m
	
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + " prj_stack volume --iter --var --sym=symmetry --MPI"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--iter", 		type="int"         ,	default=40   ,				help="maximum number of iterations (stop criterion of reconstruction process)" )
	parser.add_option("--abs", 			type="float"       ,	default=0.1  ,				help="minimum average absolute change of voxels' values (stop criterion of reconstruction process)" )
	parser.add_option("--var" , 		action="store_true",	default=False,				help="stack on input consists of variances")
	parser.add_option("--sym" , 		type="string"      ,	default="c1" ,				help="symmetry" )
	parser.add_option("--MPI" , 		action="store_true",	default=False,				help="use MPI version - works only with stack of variance as an input")
	parser.add_option("--SND",			action="store_true",	default=False,				help="compute squared normalized differences")
	parser.add_option("--CTF",			action="store_true",	default=False,				help="use CFT correction")
	parser.add_option("--VERBOSE",		action="store_true",	default=False,				help="comments")
	parser.add_option("--img_per_grp",	type="int"         ,	default=5  ,				help="images per group")
	parser.add_option("--ave2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D averages")
	parser.add_option("--var2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D variances")
	parser.add_option("--ave3D",		type="string"	   ,	default=False,				help="write to the disk reconstructed 3D average")
	

	(options,args) = parser.parse_args(arglist[1:])

	if (options.MPI and not options.var and not options.SND):
		print "There is no MPI version of procedure to extract variance from the stack of projections other than by computing squared normalized differences"
		exit()

	isRoot = True
	if options.MPI:
		from mpi import mpi_init, mpi_comm_rank, MPI_COMM_WORLD
		sys.argv = mpi_init(len(sys.argv), sys.argv)
		isRoot = (mpi_comm_rank(MPI_COMM_WORLD) == 0)

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	if len(args) == 2:
		prj_stack = args[0]
		vol_stack = args[1]
	else:
		ERROR("incomplete list of arguments","sxvariances3d",1)
		exit()

	from reconstruction import recons3d_em, recons3d_em_MPI

	global_def.BATCH = True
	print_begin_msg("sx3Dvariance")
	print_msg("Input stack							:%s\n"%(prj_stack))
	
	if not options.var and options.SND:
		from mpi	import mpi_barrier, MPI_COMM_WORLD
		if isRoot:
			from reconstruction	import recons3d_4nn
			from projection		import prep_vol, prgs
			from statistics		import im_diff
			from utilities		import get_im, model_circle, get_params_proj, set_params_proj	
			stack = prj_stack
			prj_stack = []
			proj_angles = []
			nima = EMUtil.get_image_count(stack)
			if options.VERBOSE:
				print nima
			structure = recons3d_4nn(stack, range(nima), symmetry = options.sym, npad = 4, xysize = -1, zsize = -1)
			structure, kb = prep_vol(structure)
			#structure.write_image("structure.hdf")	
			nx = get_im(stack, 1).get_xsize()
			ny = get_im(stack, 1).get_ysize()
			#nz = get_im(stack, 1).get_zsize()
			rad = -1
			mask = model_circle(int(rad), nx, ny)
			#tab = EMUtil.get_all_attributes(stack, 'xform.projection')
			for i in xrange(nima):
				imgdata = get_im(stack, i)
				#t = tab[i].get_params('spider')
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata)
				ref_prj = prgs(structure, kb, [phi, theta, psi, -s2x, -s2y])
				diff, A, B = im_diff(ref_prj, imgdata, mask)
				diff2 = diff*diff 
				set_params_proj(diff2, [phi, theta, psi, s2x, s2y])
				diff2.write_image("difference.hdf", i)
				#prj_stack.append(diff2)
		
		if options.MPI: mpi_barrier(MPI_COMM_WORLD)
		prj_stack = "difference.hdf"
			
	if not options.var and not options.SND and not options.MPI:
		t1 = time()
		from utilities		import group_proj_by_phitheta, get_params_proj, params_3D_2D, set_params_proj, set_params2D
		from utilities		import compose_transform2
		from statistics		import avgvar, avgvar_CTF
		from morphology		import threshold
		from reconstruction	import recons3d_4nn, recons3d_4nn_ctf
		stack = prj_stack
		prj_stack = []
		proj_angles = []
		aveList = []
		nima = EMUtil.get_image_count(stack)
		print_msg("Number of projections							:%d\n"%(nima))
		if options.VERBOSE:
			print "Number of projections:", nima
		tab = EMUtil.get_all_attributes(stack, 'xform.projection')
		for i in xrange(nima):
			t = tab[i].get_params('spider')
			proj_angles.append([t['phi'], t['theta'], t['psi']])
		t2 = time()
		print_msg("Number of images per group						:%d\n"%(options.img_per_grp))
		print_msg("... grouping projections \n")
		if options.VERBOSE:
			print "Number of images per group: ", options.img_per_grp	
			print "NOW GROUPING PROJECTIONS"																							
		proj_list, angles_list, mirror_list = group_proj_by_phitheta(proj_angles, options.sym, options.img_per_grp)
		t3 = time()
		print_msg("Grouping projections lasted [s]						:%s\n"%(t3-t2))				
		del proj_angles
		print_msg("Number of groups							:%d\n"%(len(proj_list)))
		if options.VERBOSE:
			print "Grouping projections lasted [min]: ", (t3-t2)/60	
			print "Number of groups: ", len(proj_list)																		
		t4 = time()
		print_msg("... calculating the stack of 2D variances \n")
		if options.VERBOSE:
			print "NOW CALCULATING A STACK OF 2D VARIANCES"							
		for i in xrange(len(proj_list)): 
			imgdata = EMData.read_images(stack, proj_list[i])
			for j in xrange(len(proj_list[i])):
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata[j])
				alpha, sx, sy, mirror = params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror_list[i][j])
				if mirror == 0:  alpha, sx, sy, scale = compose_transform2( alpha, sx, sy, 1.0, angles_list[i][0]-phi, 0.0, 0.0, 1.0)
				else:            alpha, sx, sy, scale = compose_transform2( alpha, sx, sy, 1.0, 180-(angles_list[i][0]-phi), 0.0, 0.0, 1.0)
				set_params2D(imgdata[j], [alpha, sx, sy, mirror, 1.0])
			if (options.CTF):	ave, var = avgvar_CTF(imgdata,"a")
			else:	ave, var = avgvar(imgdata, mode="a", interp="linear")
			var = threshold(var)
			set_params_proj(var, [angles_list[i][0], angles_list[i][1], 0.0, 0.0, 0.0])
			var.set_attr("imgindex",proj_list[i])
			prj_stack.append(var)
			set_params_proj(ave, [angles_list[i][0], angles_list[i][1], 0.0, 0.0, 0.0])
			ave.set_attr("imgindex",proj_list[i])
			aveList.append(ave)
			if (options.ave2D):	ave.write_image(options.ave2D,i)
			if (options.var2D): var.write_image(options.var2D,i)
		print "GOT A STACK OF 2D VARIENCE"
		if (options.ave2D): print_msg("Writing to the disk a stack of 2D averages as				:%s\n"%(options.ave2D))
		if (options.var2D): print_msg("Writing to the disk a stack of 2D variances as				:%s\n"%(options.var2D))
		print "RECONSTRUCTING 3D AVERAGE VOLUME"																					##
		ave3D = recons3d_4nn(aveList, range(len(proj_list)-1), symmetry = options.sym, npad = 4, xysize = -1, zsize = -1)
		if (options.ave3D): 
			ave3D.write_image(options.ave3D)
			print_msg("Writing to the disk volume reconstructed from averages as		:%s\n"%(options.ave3D))
		del ave, var, imgdata, angles_list, proj_list, stack, phi, theta, psi, s2x, s2y, alpha, sx, sy, mirror, aveList, ave3D
		t5 = time()
		print_msg("Calculating the stack of 2D variances lasted [s]			:%s\n"%(t5-t4))
		#exit()

	if options.MPI:
		t6 = time()
		from mpi import mpi_comm_rank, MPI_COMM_WORLD
		res = recons3d_em_MPI(prj_stack, options.iter, options.abs, True, options.sym)
		if isRoot:
			res.write_image(vol_stack)
	else:
		t6 = time()
		print_msg("... reconstructing 3D variance  \n")
		print "RECONSTRUCTING 3D VARIANCE VOLUME"																				##
		res = recons3d_em(prj_stack, options.iter, options.abs, True, options.sym)
		res.write_image(vol_stack)
		print_msg("Writing to the disk volume of reconstructed 3D variance as		:%s\n"%(vol_stack))
	t7 = time()
	print_msg("Reconstructing 3D variance lasted [s]					:%s\n"%(t7-t6))
	print "RECONSTRUCTION LASTED: ", (t7-t6)/60, " min"																		##
	tF = time()
	print_msg("Total time for these computations [min]					:%s\n"%((tF-t0)/60))
	print_end_msg("sx3Dvariance")
	global_def.BATCH = False

	if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()
	
if __name__=="__main__":
	main()
