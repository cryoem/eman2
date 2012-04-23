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

import global_def
from   global_def import *
from   optparse import OptionParser
from   string import atoi,replace
from   EMAN2 import EMUtil
import os
import sys


def calc_projections_variance(projections, members, phi, theta, delta):
	from string import replace
	from random import random
	from sparx import get_im, ave_var
	from fundamentals import rot_shift2D
	from utilities import get_params_proj, set_params_proj

	class_data = []
	#print "========", phi, theta
	for im in members: 
		p = get_im(projections, im)
		t1,t2,psi,sx,sy = get_params_proj(p)
		class_data.append(p)
		p = rot_shift2D(p, -psi, -sx, -sy)
	#	print t1, t2, psi
	projs_avg,projs_var = ave_var(class_data, "")
	projs_var.set_attr('delta', delta)
	set_params_proj(projs_var, [phi, theta, 0.0, 0.0, 0.0])
	projs_var.set_attr('members', members)
	projs_var.set_attr('n_objects', len(members))
	return projs_var


def find_projections_variances(projections_stack_filename, symmetry = "c1"):
	from sparx import even_angles, assign_projangles
	from time import time
	
	projections_variances_list = []   # result
	
	# The fourth and fifth columns are originally for shifts,
	# but we will recycle them and use the fourth column for the image number,
	# and the fifth one for the active flag

	projections = EMData.read_images(projections_stack_filename)
	
	proj_ang = []
	for i in xrange(len(projections)):
		proj = projections[i]
		RA = proj.get_attr( "xform.projection" )
		angdict = RA.get_sym(symmetry,0).get_rotation("spider")
		proj_ang.append( [angdict["phi"], angdict["theta"], angdict["psi"], i, True] )
			
	img_per_grp = 40 # int(sys.argv[4])
	diff_pct = 0.1 # float(sys.argv[5])
	min_img_per_grp = img_per_grp*(1-diff_pct)
	max_img_per_grp = img_per_grp*(1+diff_pct)
	
	N = len(proj_ang)

	while True:
		proj_ang_now = []
		for i in xrange(N):
			if proj_ang[i][4]: proj_ang_now.append(proj_ang[i])
		#print "Current size of data set = ", len(proj_ang_now)
		if len(proj_ang_now) <= max_img_per_grp:
			members = [0]*len(proj_ang_now)
			for i in xrange(len(proj_ang_now)):  members[i] = proj_ang_now[i][3]
			print "Size of this group = ", len(members)
			projections_variances_list.append( calc_projections_variance(projections, members, 0.0, 0.0, 45.0) )
			break
	
		min_delta = 0.001
		max_delta = 2.0
		trial = 0
		max_trial = 20
		while trial < max_trial:
			mid_delta = (min_delta+max_delta)/2
			#print "...... Testing delta = %6.2f"%(mid_delta)
			ref_ang = even_angles(mid_delta, symmetry=symmetry)
			t1 = time()
			asg = assign_projangles(proj_ang_now, ref_ang)
			#print "............ Time used = %6.2f"%(time()-t1)

			count = []
			for i in xrange(len(asg)): count.append([len(asg[i]), i])
			count.sort(reverse = True)
			#print "............ maximum size = %5d   minimum size = %5d"%(count[0][0], count[-1][0])
			k = 0
			grp_size = count[k][0]
			grp_num = count[k][1]
			while grp_size >= min_img_per_grp and grp_size <= max_img_per_grp:
				members = [0]*grp_size
				for i in xrange(grp_size):
					members[i] = proj_ang_now[asg[grp_num][i]][3]
					proj_ang[members[i]][4] = False
				#print "Size of this group = ", grp_size
				projections_variances_list.append( calc_projections_variance(projections, members, ref_ang[grp_num][0], ref_ang[grp_num][1], mid_delta) )
				k += 1
				grp_size = count[k][0]
				grp_num = count[k][1]
			if k > 0: break
			if grp_size > max_img_per_grp:  max_delta = mid_delta
			else:  min_delta = mid_delta
			trial += 1

		if trial == max_trial:
			grp_size = count[0][0]
			grp_num = count[0][1]
			members = [0]*grp_size
			for i in xrange(grp_size):
				members[i] = proj_ang_now[asg[grp_num][i]][3]
				proj_ang[members[i]][4] = False
			#print "Size of this group = ", grp_size
			projections_variances_list.append( calc_projections_variance(projections, members, ref_ang[grp_num][0], ref_ang[grp_num][1], mid_delta) )
		#print ""
	return projections_variances_list
	

def main():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + " prj_stack volume --iter --var --sym=symmetry --MPI"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--iter", type="int"         , default=20   , help="Maximum number of iterations" )
	parser.add_option("--var" , action="store_true", default=False, help="stack on input consists of variances")
	parser.add_option("--sym" , type="string"      , default="c1" , help="symmetry" )
	parser.add_option("--MPI" , action="store_true", default=False, help="use MPI version ")

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
		prj_stack = find_projections_variances(prj_stack, options.sym)
	
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
