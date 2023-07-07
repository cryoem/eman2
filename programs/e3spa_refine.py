#!/usr/bin/env python
#
# Author: Steven Ludtke, 05/03/2023 (sludtke@bcm.edu)
# Copyright (c) 2023 Baylor College of Medicine
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
# todo: verify the processors who have the same names in proc3d
#	   and proc2d have the same implementation
#

from EMAN2 import *
import numpy as np
import tensorflow as tf
import h5py

def main():

	usage="""Gaussian based single particle refinement.

	Particle .lst files must have necessary parameters embedded either in the original image header or in the lst file itself. At a bare minimum this would be basic CTF information. If continuing from a previous refinement, particle orientations may also be present. If an input model isn't specified, one will be generated from scratch.

	e2spa_refine.py --ptclsin sets/all__orig.lst --ptclsout spa/iter_01.lst

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
# 	parser.add_argument("--chunk", type=str,help="In case of memory exhaustion, particles can be trained in chunks (similar to a batch). eg - 0,10 will process the first of 10 chunks. Be sure to store and restore the encoder and decoder if using this mode. Must be run sequentially 0-(n-1) without parallelism.", default=None)
# 	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--modelin", type=str,help="load from an existing model file", default="")
	parser.add_argument("--modelout", type=str,help="output trained model file. only used when --projs is provided", default="")
	parser.add_argument("--ptclsin", type=str,help="particles input for alignment", default="")
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default="")
# #	parser.add_argument("--segments", type=str,help="Divide the model into sequential domains. Comma separated list of integers. Each integer is the first sequence number of a new region, starting with 0",default=None)
# 	parser.add_argument("--decoderin", type=str,help="Rather than initializing the decoder from a model, read an existing trained decoder", default="")
# 	parser.add_argument("--decoderout", type=str,help="Save the trained decoder model. Filename should be .h5", default=None)
# 	parser.add_argument("--encoderin", type=str,help="Rather than initializing the encoder from scratch, read an existing trained encoder", default=None)
# 	parser.add_argument("--encoderout", type=str,help="Save the trained encoder model. Filename should be .h5", default=None)
# 	parser.add_argument("--projs", type=str,help="projections with orientations (in hdf header or comment column of lst file) to train model", default="")
# 	parser.add_argument("--evalmodel", type=str,help="generate model projection images to the given file name", default="")
# 	parser.add_argument("--evalsize", type=int,help="Box size for the projections for evaluation.", default=-1)
# 	parser.add_argument("--ptclsclip",type=int,help="clip particles to specified box size before use",default=-1)
# 	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. Default is 1e-4. ", default=1e-4)
# #	parser.add_argument("--sigmareg", type=float,help="regularizer for the sigma of gaussian width. Larger value means all Gaussian functions will have essentially the same width. Smaller value may help compensating local resolution difference.", default=.5)
# 	parser.add_argument("--modelreg", type=float,help="regularizer for for Gaussian positions based on the starting model, ie the result will be biased towards the starting model when training the decoder (0-1 typ). Default 0", default=0)
# 	parser.add_argument("--ampreg", type=float,help="regularizer for the Gaussian amplitudes in the first 1/2 of the iterations. Large values will encourage all Gaussians to have similar amplitudes. default = 0", default=0)
# 	parser.add_argument("--niter", type=int,help="number of iterations", default=10)
# 	parser.add_argument("--npts", type=int,help="number of points to initialize. ", default=-1)
# 	parser.add_argument("--batchsz", type=int,help="batch size", default=256)
# 	parser.add_argument("--minressz", type=int,help="Fourier diameter associated with minimum resolution to consider. ", default=4)
# 	parser.add_argument("--maxboxsz", type=int,help="maximum fourier box size to use. 2 x target Fourier radius. ", default=64)
# 	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
# 	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
# 	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")
# 	parser.add_argument("--decoderentropy", action="store_true", default=False ,help="This will train some entropy into the decoder using particles to reduce vanishing gradient problems")
# 	parser.add_argument("--perturb", type=float, default=0.1 ,help="Relative perturbation level to apply in each iteration during --heter training. Default = 0.1, decrease if models are too disordered")
# 	parser.add_argument("--conv", action="store_true", default=False ,help="Use a convolutional network for heterogeneity analysis.")
# 	parser.add_argument("--fromscratch", action="store_true", default=False ,help="start from coarse alignment. otherwise will only do refinement from last round")
# 	parser.add_argument("--ptclrepout", type=str,help="Save the per-particle representation input to the network to a file for later use", default="")
# 	parser.add_argument("--ptclrepin", type=str,help="Load the per-particle representation from a file rather than recomputing", default="")
# 	parser.add_argument("--midout", type=str,help="middle layer output", default="")
# 	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. sigma is not supported in this version of the program. use 3 digit 0/1 input. default is 110, i.e. only adjusting position and amplitude", default="110")
# 	parser.add_argument("--nmid", type=int,help="size of the middle layer. If model is grouped becomes nmid*ngroup", default=4)
# 	parser.add_argument("--mask", type=str,help="remove points outside mask", default="")
# 	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	logid=E2init(sys.argv,options.ppid)

