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
from	global_def import *
from	optparse import OptionParser
from	string import atoi,replace
from	EMAN2 import EMUtil
import	os
import	sys


def main():
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
	parser.add_option("--img_per_grp",	type="int"         ,	default=100  ,				help="images per group")
	parser.add_option("--diff_pct", 	type="float"       ,	default=0.1  ,				help="percentage of ...")		
	parser.add_option("--CTF",			action="store_true",	default=False,				help="use CFT correction")
	parser.add_option("--ave2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D averages")
	parser.add_option("--var2D",		type="string"	   ,	default=False,				help="write to the disk a stack of 2D variances")
	parser.add_option("--ave3D",		type="string"	   ,	default=False,				help="write to the disk reconstructed 3D average")

	(options,args) = parser.parse_args(arglist[1:])

	if (options.MPI and not options.var):
		print "There is no MPI version of procedure to extract variance from the stack of projections"
		exit()

	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)

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

	if not options.var:
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
		tab = EMUtil.get_all_attributes(stack, 'xform.projection')
		for i in xrange(nima):
			t = tab[i].get_params('spider')
			proj_angles.append([t['phi'], t['theta'], t['psi']])
		proj_list, angles_list = group_proj_by_phitheta(proj_angles, options.sym, options.img_per_grp, options.diff_pct)
		del proj_angles
		print "Number of groups = ", len(proj_list)
		for i in xrange(len(proj_list)): 
			imgdata = EMData.read_images(stack, proj_list[i])
			for j in xrange(len(proj_list[i])):
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata[j])
				alpha, sx, sy, mirror = params_3D_2D(phi, theta, psi, s2x, s2y)
				if mirror == 0:  alpha, sx, sy, scale = compose_transform2( alpha, sx, sy, 1.0, angles_list[i][0]-phi, 0.0, 0.0, 1.0)
				else:            alpha, sx, sy, scale = compose_transform2( alpha, sx, sy, 1.0, 180-(angles_list[i][0]-phi), 0.0, 0.0, 1.0)
				set_params2D(imgdata[j], [alpha, sx, sy, mirror, 1.0])
			if (options.CTF):	ave, var = avgvar_CTF(imgdata,"a")
			else:	ave, var = avgvar(imgdata,"a")
			var = threshold(var)
			set_params_proj(var, [angles_list[i][0], angles_list[i][1], 0.0, 0.0, 0.0])
			var.set_attr("imgindex",proj_list[i])
			prj_stack.append(var)
			set_params_proj(ave, [angles_list[i][0], angles_list[i][1], 0.0, 0.0, 0.0])
			ave.set_attr("imgindex",proj_list[i])
			aveList.append(ave)
			if (options.ave2D):	ave.write_image(options.ave2D,i)
			if (options.var2D): var.write_image(options.var2D,i)
		if (options.CTF):	ave3D = recons3d_4nn_ctf(aveList, range(len(proj_list)-1), snr = 1.0, sign = 1, symmetry = options.sym, verbose = 0, npad = 4, xysize = -1, zsize = -1)
		else:	ave3D = recons3d_4nn(aveList, range(len(proj_list)-1), symmetry = options.sym, npad = 4, xysize = -1, zsize = -1)
		if (options.ave3D): ave3D.write_image(options.ave3D)
		del ave, var, imgdata, angles_list, proj_list, stack, phi, theta, psi, s2x, s2y, alpha, sx, sy, mirror, ave3D, aveList
		#exit()

	if options.MPI:
		from mpi import mpi_comm_rank, MPI_COMM_WORLD
		res = recons3d_em_MPI(prj_stack, options.iter, options.abs, True, options.sym)
		mpi_r = mpi_comm_rank(MPI_COMM_WORLD)
		if mpi_r == 0:
			res.write_image(vol_stack)
	else:
		res = recons3d_em(prj_stack, options.iter, options.abs, True, options.sym)
		res.write_image(vol_stack)
	global_def.BATCH = False

	if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()


if __name__=="__main__":
	main()
