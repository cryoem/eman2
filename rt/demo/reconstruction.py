#!/usr/bin/env python
"""Iterative reconstruction program.

"""
from EMAN2 import *
import os,os.path

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
	syscommand = "rm -rf %s" % rpdir
	os.system(syscommand)
	# do 3-d reconstruction to create a new model
	vol = do_reconstruction(alipattern, data_start, data_end,
			                4, exptanglelist)
	rootname, extname = os.path.splitext(output_volume_name)
	vol.write_image(rootname + "_" + repr(deltheta) + extname)

