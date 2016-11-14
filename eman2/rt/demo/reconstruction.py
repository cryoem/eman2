#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

"""Iterative reconstruction program.

"""
from EMAN2 import *
import os,os.path
import shutil

# stuff the user would change
init_model_name = "model.tcp"
output_volume_name = "myvol.mrc"
reference_projection_filepattern = "ref{****}.tcp"
alignment_projection_filepattern = "ali{****}.tcp"
exptal_projection_filepattern = "expt{***}.tcp"
data_start = 0
data_end = 19
proj_radius = 35 # pixels
theta_granularity = [10, 5, 3, 2]
# end of stuff user should change

# get initial model
vol = getImage(init_model_name)
# iterate reconstruction 
for deltheta in theta_granularity:
	# make reference temp directory and reference path+fname pattern 
	rpdir = "tmpref"
	if not os.path.isdir(rpdir):
		os.mkdir(rpdir)
	refpattern = os.path.join(rpdir,reference_projection_filepattern)
	# use current model to create reference projections
	anglelist = Util.voea(deltheta)
	create_write_projections(vol, refpattern, anglelist, proj_radius)
	# take expt'al and reference projections and perform alignment
	# make alignment temp directory and reference path+fname pattern 
	aldir = "tmpalign"
	if not os.path.isdir(aldir):
		os.mkdir(aldir)
	alipattern = os.path.join(aldir, alignment_projection_filepattern)
	exptanglelist = do_alignment(exptal_projection_filepattern,
			                     data_start, data_end,
								 refpattern, alipattern, anglelist)
	# remove reference projections
	shutil.rmtree(rpdir, ignore_errors=True)
	# do 3-d reconstruction to create a new model
	vol = do_reconstruction(alipattern, data_start, data_end,
			                4, exptanglelist)
	rootname, extname = os.path.splitext(output_volume_name)
	vol.write_image(rootname + "_" + repr(deltheta) + extname)
