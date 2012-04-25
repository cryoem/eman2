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

	parser.add_option("--iter", 		type="int"         ,	default=20   ,	help="Maximum number of iterations" )
	parser.add_option("--var" , 		action="store_true",	default=False,	help="stack on input consists of variances")
	parser.add_option("--sym" , 		type="string"      ,	default="c1" ,	help="symmetry" )
	parser.add_option("--MPI" , 		action="store_true",	default=False,	help="use MPI version")
	parser.add_option("--img_per_grp",	type="int"         ,	default=100  ,	help="images per group")
	parser.add_option("--diff_pct", 	type="float"       ,	default=0.1  ,	help="percentage of ...")		
	parser.add_option("--CTF",			action="store_true",	default=False,	help="use CFT correction")

	(options,args) = parser.parse_args(arglist[1:])

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
		from utilities	import group_proj_by_phitheta, get_params_proj, params_3D_2D, set_params_proj, set_params2D
		from statistics	import avgvar, avgvar_CTF
		from morphology	import threshold
		stack = prj_stack
		prj_stack = []
		proj_angles = []
		nima = EMUtil.get_image_count(stack)
		tab = EMUtil.get_all_attributes(stack, 'xform.projection')
		for i in xrange(nima):
			t = tab[i].get_params('spider')
			proj_angles.append([t['phi'], t['theta'], t['psi']])
			##if (i < 30): print t['psi']  ##
		proj_list, angles_list = group_proj_by_phitheta(proj_angles, options.sym, options.img_per_grp, options.diff_pct)
		del proj_angles
		#print "Number of groups = ", len(proj_list)
		for i in xrange(len(proj_list)):
			imgdata = EMData.read_images(stack, proj_list[i])
			for j in xrange(len(proj_list[i])):
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata[j])
				alpha, sx, sy, mirror = params_3D_2D(phi, theta, psi, s2x, s2y)
				set_params2D(imgdata[j], [alpha, sx, sy, mirror, 1.0])
			if (options.CTF):	ave, var = avgvar_CTF(imgdata,"a")
			else:	ave, var = avgvar(imgdata,"a")
			var = threshold(var)
			set_params_proj(var, [angles_list[i][0], angles_list[i][1], 0.0, 0.0, 0.0])
			prj_stack.append(var)
			ave.write_image("ave2Dstack.hdf",i)
			var.write_image("var2Dstack.hdf",i) 
		del ave, var, imgdata, angles_list, proj_list, stack, phi, theta, psi, s2x, s2y, alpha, sx, sy, mirror
		exit()
		
		
	if options.MPI:
		from mpi import mpi_comm_rank, MPI_COMM_WORLD
		res = recons3d_em_MPI(prj_stack, options.iter, 0.01, True, options.sym)
		mpi_r = mpi_comm_rank(MPI_COMM_WORLD)
		if mpi_r == 0:
			res.write_image(vol_stack)
	else:
		res = recons3d_em(prj_stack, options.iter, 0.01, True, options.sym)
		res.write_image(vol_stack)
	global_def.BATCH = False

	if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()


if __name__=="__main__":
	main()
