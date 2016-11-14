#!/usr/bin/env python
#
# Author: Jesus Galaz-Montoya, 2011?2012?
# Last update 25/Feb/2015
# Copyright (c) 2011 Baylor College of Medicine
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

import os
from EMAN2 import *
from sys import argv
from optparse import OptionParser
import sys
import numpy as np
import math
from scipy.stats import norm
from e2spt_intrafsc import genOddAndEvenVols, fscOddVsEven


def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Program for individual subtomogram refinement (ISR) based on subtiltseries
	for each subtomogram extracted with e2spt_subtilt.py."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument('--input',type=str,default='',help="""Comma separated files in .ali, .st .hdf format of the aligned subtiltseries.""")
	
	parser.add_argument('--inputstem',type=str,default='',help="""Alternative to supplying --input. This is a string common to multiple files to be processed in the CURERENT directory. The common string doesn't need to be at a particular location in the filenames. For example, a series of files "tiltA.hdf, tiltB.hdf, tiltC.hdf" could have either 'hdf', '.hdf', 't,','ti', 'til', 'tilt', etc., as a common string. The key is to choose a string shared ONLY by the files of interest. The files should be multiple subtiltseries in .hdf format; each file should correspond to an individual subtiltseries for a different particle: That is, each file should be a subtiltseries corresponding to an individual subtomogram, as extracted by e2spt_subtilt.py, or as simulated by e2spt_simulation.py""")

	parser.add_argument('--inputdir',type=str,default='',help="""Alternative to --input and --inputstem. Path to a directory containing individual subtiltseries stacks.""")
	
	parser.add_argument("--path",type=str,default='',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptisr'; for example, sptisr02 will be the directory by default if 'sptisr_01' already exists.")
		
	parser.add_argument("--ppid", type=int, help="Default=1. Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", type=int, default=0, help="Default 0. Verbose level [0-9], higner number means higher level of verboseness",dest="verbose", action="store", metavar="n")

	parser.add_argument("--reconstructor", type=str,default="fourier:mode=gauss_2",help="""Default=fourier:mode=gauss_2. The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line to see all options and parameters available. To specify the interpolation scheme for the fourier reconstructor, specify 'mode'. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5'. For example --reconstructor=fourier:mode=gauss_5 """)
	
	parser.add_argument("--iter",type=int,default=1,help="""Number of iterations to run algorithm for.""")
			
	parser.add_argument("--tltfile",type=str,default='',help="""IMOD-like .tlt file with tilt angles for the aligned tiltseries (or set of subtiltseries).""")
	
	parser.add_argument("--tiltaxis",type=str,default='y',help="""Axis to produce projections about. Default is 'y'; the only other valid option is 'x'.""")
			
	parser.add_argument("--pad2d", type=float,default=0.0,help="""Default=0.0. Padding factor (e.g., 2.0, to make the box twice as big) to zero-pad the 2d images in the tilt series for reconstruction purposes (the final reconstructed subvolumes will be cropped back to the original size though).""")

	parser.add_argument("--pad3d", type=float,default=0.0,help="""Default=0.0. Padding factor (e.g., 2.0, to make the box twice as big) to zero-pad the volumes for reconstruction purposes (the final reconstructed subvolumes will be cropped back to the original size though).""")
	
	parser.add_argument("--savevols",action='store_true',default=False,help="""This option will save the reconstructed volumes at each iteration.""")
		
	parser.add_argument("--outxsize",type=int,default=0,help='''Clip the output volume in x to this size. The default size is the nx size of the input images.''')	
	
	parser.add_argument("--outysize",type=int,default=0,help='''Clip the output volume in y to this size. The default size is the ny size of the input images.''')
	
	parser.add_argument("--outzsize",type=int,default=0,help='''Clip the output volume in z to this size. The default size is the nx size of the input images.''')
			
	parser.add_argument("--mask",type=str,help="""Default=None. Masking processor (see e2help.py --verbose=10) applied to the images to aid alignment. Default=None.""")
	
	parser.add_argument("--preprocess",type=str,help="""Default=None. Any processor (see e2help.py --verbose=10) applied to the images to aid alignment.""")
	
	parser.add_argument("--lowpass",type=str,default='',help="""Default=None. A lowpass filtering processor (see e2help.py --verbose=10) applied to each volume prior to reprojection generation..""")
	
	parser.add_argument("--highpass",type=str,default='',help="""Default=None. A highpass filtering processor (see e2help.py --verbose=10) applied to each volume prior to reprojection generation.""")	
		
	parser.add_argument("--threshold",type=str,default='',help="""Default=None. A threshold  processor (see e2help.py --verbose=10) applied to each volume prior to reprojection generation.""")
	
	parser.add_argument("--saveali", action="store_true", default=False, help="""Default=False. If set, will save the recentered subtiltseries after each iteration.""")

	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)
	
	
	if not options.input and not options.inputdir and not options.inputstem:
		print "ERROR: Either of the following required: --input, --inputstemp, --inputdir"
		sys.exit()
	
	
	if options.input and options.inputstem:
		print "ERROR: Cannot provide --input and --inputstem simultaneously"
		sys.exit()
		
	if options.inputstem and options.inputdir:
		print "ERROR: Cannot provide --inputstem and --inputdir simultaneously"
		sys.exit()
		
	if options.input and options.inputdir:
		print "ERROR: Cannot provide --input and --inputdir simultaneously"
		sys.exit()
	
	
	from e2spt_classaverage import sptOptionsParser
	options = sptOptionsParser( options )
	#print "Options have been parsed, for example, mask is", options.mask, type(options.mask), len(options.mask)
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptisr')
	
	originalpath = options.path
	
	if options.verbose > 9:
		print "\n(e2spt_isr.py) I've read the options"	
	
	inputfiles = {}											#C:Figure out whether there's a single HDF stack to process,
															#C:or a directory with many HDF stacks
	c = os.getcwd()

	if options.inputdir:
		c = os.getcwd() + '/' + options.inputdir

	findir = os.listdir( c )
	
	if options.inputstem:
		for f in findir:
			if '.hdf' in f and options.inputstem in f:
				if options.verbose > 8:
					print "\nFound tiltseries!", f
				inputfiles.update( {f:[f,None]} )			#C:The input files are put into a dictionary in the format {originalseriesfile:[originalseriesfile,volumefile]}
	
	elif options.inputdir:
		for f in findir:
			if '.hdf' in f:
				if options.verbose > 8:
					print "\nFound tiltseries!", f
				inputfiles.update( {f:[f,None]} )			#C:The input files are put into a dictionary in the format {originalseriesfile:[originalseriesfile,volumefile]}
	
	elif options.input:
		inputfiles.update( {options.input:[options.input,None]} )
	
	
	tiltstep = 0
	newvol = vol = None
	newfiles = {}
	firstiterdir = originalpath
	
	print "\n\nThere are these many iterations", options.iter
	fstats={}
	convergedfs=[]
	
	fmeanscores2d = {}
	itermeanscores2d = {}
	
	fmeanerrors = {}
	itermeanerrors = {}
	
	fmeanfscs = { -1:{} } #Compute initial FSC before any iterations of correction or alignment have occurred
	itermeanfscs = {}
	
	fmeanscores3d = {-1:{}}
	itermeanscores3d = {}
	
	errorswitch = 0
	
	
	for i in range( options.iter ):		#Iterate over options.iter
		itermeanscores2d.update( {i:0} )
		itermeanerrors.update( {i:0} )
		itermeanfscs.update( {i:0} )
		itermeanscores3d.update( {i:0} )
		
		print "***************\nStarting iteration number", i
		print "****************"
		previouspath=''
		if options.iter > 1:			
			iterdir = 'iter_' + str( i+1 ).zfill( len ( str( options.iter ) ) )
			os.system( 'mkdir ' + originalpath + '/' + iterdir )	#C:Make a new subdirectory within options.path only if options.iter > 0
			options.path =  originalpath + '/' + iterdir			#C:Update path to include the newly created subdirectory
			previouspath = originalpath + '/iter_' + str( i ).zfill( len ( str( options.iter ) ) )
			if i == 0:
				firstiterdir = options.path
		kk = 0
		
		statsfile = options.path + '/stats.txt'
		
		statslines=[]
		
		if len( convergedfs ) == len( inputfiles ):
			print "\nAll files have converged. Terminating. Line 193"
			break
		
		
		fmeanscores2d.update( { i:{} } )
		fmeanerrors.update( { i:{} } )
		fmeanfscs.update( { i:{} } )
		fmeanscores3d.update( { i:{} } )
		
		for f in inputfiles:
			
			fmeanscores2d[i].update( { f:[0] } )
			fmeanerrors[i].update( { f:[0] } )	
			fmeanfscs[i].update( { f:[0] } )
			fmeanscores3d[i].update( { f:[0] } )											#C:Iterate over files
	
			#originalseries = f
			
			if i == 0:
				fstats.update({f:list([])})
				print "fmeanfscs is", fmeanfscs
				fmeanfscs[-1].update( { f:[0] } ) 	#Compute initial FSC before any iterations of correction or alignment have occurred
				print "fmeanscores3d is", fmeanscores3d
				fmeanscores3d[-1].update( { f:[0] } )
				print "set initial -1 fmeanscores3d!!"
	
			if 'converged' not in fstats[f]:
			
				if i == 0:
					if options.tltfile:
						originalangles = getangles( options )				#C:Read angles from tlt file if available
					else:
						originalangles = calcangles( f )		#C:Get the angles of the actual data tilts from the header or calculate from input parameters
			
				stackfile = options.path + '/' + os.path.basename(f).replace('.hdf','_ISR.hdf')
			
				print "\nWill refine center for file", f
				if i ==0:
					hdr = EMData( f, 0, True )				#See if angles are in the header of the data, for each file; if not, write them
					#print "\nRead header and its type", hdr, type(hdr)
					#print "\nAnd the dictionary", hdr.get_attr_dict()
					size = hdr['nx']
				
					aux=0
					if 'spt_tiltangle' not in hdr.get_attr_dict():
						print "\nspt_tiltangle not in header, therefore will write it by calling writeparamtoheader"
						aux = writeparamtoheader( f, originalangles, 'spt_tiltangle' )
						print "\naux returned is", aux
					else:
						aux = 1			 
				
					if 'sptisrtx' not in hdr.get_attr_dict() or 'sptisrty' not in hdr.get_attr_dict() or 'sptisrdr' not in hdr.get_attr_dict():
						tvals=[0.0 for tt in range(len(originalangles)) ]
					
						auxx=auxy=auxr=0
						if 'sptisrtx' not in hdr.get_attr_dict():
							print "\nsptisrtx not in header, therefore will write it by calling writeparamtoheader"
							auxx = writeparamtoheader( f, tvals, 'sptisrtx' )
						else:
							auxx = 1 
					
						if 'sptisrty' not in hdr.get_attr_dict():
							print "\nsptisrty not in header, therefore will write it by calling writeparamtoheader"
							auxy = writeparamtoheader( f, tvals, 'sptisrty' )
						else:
							auxy = 1 
						
						if 'sptisrdr' not in hdr.get_attr_dict():
							print "\nsptisrdr not in header, therefore will write it by calling writeparamtoheader"
							auxr = writeparamtoheader( f, tvals, 'sptisrdr' )
						else:
							auxr = 1 
						
						if auxx and auxy and auxr:
							aux = 1
					else:
						aux = 1
						
					if aux:
						series={}
						nimgs=EMUtil.get_image_count(f)
						for ii in range(nimgs):
							img=EMData(f,ii)
							try:
								angle=img['spt_tiltangle']
								series.update({ angle:img })
							except:
								print "ERROR: spt_tiltangle not found in image", ii
								print "\nHeader is", img.get_attr_dict()
								sys.exit()

						
						retm = makevol( options, f, series, i, originalangles, size, writevols = 1, initialfsc = 1  )	#C:In the first iteration, first reconstruct the tilt series into a 3D volume
					
						vol = retm[0]
						#newvolfile = retm[1]

						fscarea = retm[2]
						score3d = retm[3]
						
						#fsc = retm[3]
						#initialfscfilename = options.path + '/' + os.path.basename( f ).replace('.hdf', '_initial_evenOddFSC.txt')
						
						fmeanfscs[ -1 ][f].append( fscarea )
						fmeanscores3d[ -1 ][f].append( score3d )
						
						
						print "\nVol and its type are", vol,type(vol)
					else:
						print "ERROR: Something went wrong. spt_tiltangle found in image headers, but somehow is unusuable"
						sys.exit()
					
				elif i>0:
					#vol = EMData( inputfiles[f][1] )	#C:For iterations > 0 (after the first iteration), the updated volumes 
														#C:for each tilt series should exist from the previous iteration						
				
					#previousstackfile = previouspath + '/' + os.path.basename(f).replace('.hdf','_IPET.hdf')

					#print "\n\n\previousstackfile is", previousstackfile
					#vol = EMData( previousstackfile, 0 )
				
					vol = newvol
				
				
				#C:Make reprojections from all known angles
				reprojections = reprojectvolume(options, vol, f, originalangles )			
			
				originalsize = vol['nx']
				retrct = recentertilts( options, reprojections, f, originalangles, i )					
				newseries = retrct[0]
				statsline = retrct[1]
				fscores2d = retrct[2]
				ferrors = retrct[3]
				
				if ferrors:
					errorswitch = 1
				
				
				#print "fscores2d received", fscores2d
				#print "ferrors received", ferrors
				if ferrors == fscores2d:
					print "ERROR: errors and scores2d are the same"
					sys.exit() 
				
				fmeanscore2d = sum(fscores2d)/len(fscores2d)
				fmeanscores2d[ i ][f].append( fmeanscore2d )
				
				if ferrors and errorswitch:
					fmeanerror =  sum(ferrors)/len(ferrors)
					fmeanerrors[ i ][f].append( fmeanerror )
				
				#fmeanscores2d.update( i:{ f:fmeanscore2d } )
			
				statslines.append(statsline)
			
				line2append = statsline
				
				#print "fstats[f] and type are", fstats[f], type(fstats[f])
				#sys.exit()
				if statsline in fstats[f]:
					line2append = 'converged'
					convergedfs.append( f )
				
				fstats[f].append( line2append )
			
				
				fstats.update({f:fstats[f]})
			
				retmkvol = makevol( options, f, newseries, i, originalangles, originalsize, writevols = 1 )
			
				newvol = retmkvol[0]
				newvolfile = retmkvol[1]
				
				fscarea = retmkvol[2]
				fmeanfscs[ i ][f].append( fscarea )
				
				fmeanscore3d = retmkvol[3]
				fmeanscores3d[ i ][f].append( fmeanscore3d )
			
				if i == 0:
					pass
					#print "\nEEEEEEEEEEEEEEE\n\n\n\n\nNewvolfile returned is", newvolfile
		
				newfiles.update( {f:[newseries,newvolfile]} )
				#newvol.write_image( stackfile, 0 )
			
				kk+=1				
			
				#firstrawvolfile = firstiterdir + '/' + os.path.basename(f).replace('.hdf','_3D.hdf')
				#fscfile = options.path + '/' + os.path.basename(f).replace('.hdf','FSC.txt')
				#cmdfsc = 'e2proc3d.py ' + firstrawvolfile + ' ' + fscfile + ' --calcfsc=' + newvolfile 
				#os.popen( cmdfsc )
			else:
				statsline = fstats[f][-2]
				statslines.append(statsline)
				
				
		#Write mean tx, ty and dr of all files to a text file, for each iteration	
		f=open( statsfile, 'w' )
		f.writelines( statslines )
		f.close()
		
		inputfiles = newfiles
		
		iterscores2d=[]
		#print "in iteration", i
		#print "fmeanscores2d are",fmeanscores2d
		for ff in fmeanscores2d[i]:
			#print "therefore ff in fmeanscores2d[i] is", ff 
			iterscores2d.append( fmeanscores2d[i][ff][-1] )
			
		itermeanscore2d = sum(iterscores2d)/len( iterscores2d )
		itermeanscores2d[i] = itermeanscore2d
		
		ys = []
		lines=[]
		for s in itermeanscores2d: 
			y = itermeanscores2d[s]
			ys.append( y )
			#print "appending this y", itermeanscores2d[s]
			line = str(s) + '\t' + str(y) + '\n'
			lines.append(line)
		
		#if plots:
		fs = open( options.path +'/scores2d_' + str(i+1) + '.txt','w')
		fs.writelines( lines )
		fs.close()
		
		try:
			if errorswitch:
				itererrors=[]
				#print "in iteration", i
				#print "fmeanerrors are",fmeanerrors
				for fscores2df in fmeanerrors[i]:
					#print "therefore ff in fmeanerrors[i] is", ff 
					itererrors.append( fmeanerrors[i][ff][-1] )
		
				itermeanerror = sum(itererrors)/len( itererrors )
				itermeanerrors[i] = itermeanerror
		
				yse = []
				linese=[]
				for s in itermeanerrors: 
					ye = itermeanerrors[s]
					yse.append( ye )
					#print "appending this error", itermeanerrors[s]
					linee = str(s) + '\t' + str(ye) + '\n'
					linese.append(linee)
		
				#if plots:
				fse = open( options.path +'/error_' + str(i+1) + '.txt','w')
				fse.writelines( linese )
				fse.close()
		except:
			pass
		
		
		
		if i == 0:
			iterfscs=[]
			#print "in iteration", i
			#print "APPENDING INITIAL fmeanfscs are", fmeanfscs
			for ff in fmeanfscs[-1]:
				#print "therefore INITIAL ff in fmeanfscs[-1] are", ff 
				iterfscs.append( fmeanfscs[-1][ff][-1] )
		
			itermeanfsc = sum(iterfscs)/len( iterfscs )
			itermeanfscs[-1] = itermeanfsc
		
		
			ysf = []
			linesf=[]
			for s in [-1,0]: 
				yf = itermeanfscs[s]
				ysf.append( yf )
				#print "for s", s
				#print "appending this fscarea", itermeanfscs[s]
				linef = str(s+1) + '\t' + str(yf) + '\n'
				linesf.append(linef)
		
			#if plots:
			fsf = open( options.path +'/fscsareas_' + str(i+1) + '.txt','w')
			fsf.writelines( linesf )
			fsf.close()
			
			
			iterscores3d=[]
			#print "in iteration", i
			#print "fmeanscores2d are",fmeanscores2d
			for ff in fmeanscores3d[-1]:
				#print "therefore ff in fmeanscores2d[i] is", ff 
				iterscores3d.append( fmeanscores3d[-1][ff][-1] )
			
			itermeanscore3d = sum(iterscores3d)/len( iterscores3d )
			itermeanscores3d[-1] = itermeanscore3d
		
			ys3 = []
			lines3=[]
			for s in [-1,0]: 
				y3 = itermeanscores3d[s]
				ys3.append( y )
				#print "appending this y", itermeanscores2d[s]
				line3 = str(s+1) + '\t' + str(y3) + '\n'
				lines3.append(line3)
		
			#if plots:
			fs = open( options.path +'/scores3d_' + str(i+1) + '.txt','w')
			fs.writelines( lines3 )
			fs.close()
			
		
		
		iterfscs=[]
		#print "in iteration", i
		#print "fmeanfscs are",fmeanfscs
		for ff in fmeanfscs[i]:
			#print "therefore ff in fmeanfscs[i] are", ff 
			iterfscs.append( fmeanfscs[i][ff][-1] )
		
		itermeanfsc = sum(iterfscs)/len( iterfscs )
		itermeanfscs[i] = itermeanfsc
		
		ysf = []
		linesf=[]
		fscsindxs = [s for s in itermeanfscs]
		fscsindxs.sort()
		for s in fscsindxs: 
			yf = itermeanfscs[s]
			ysf.append( yf )
			#print "for s", s
			#print "appending this fscarea", itermeanfscs[s]
			linef = str(s+1) + '\t' + str(yf) + '\n'
			linesf.append(linef)
		
		#if plots:
		fsf = open( options.path +'/fscsareas_' + str(i+1) + '.txt','w')
		fsf.writelines( linesf )
		fsf.close()
		
		
		
		
		iterscores3d=[]
		#print "in iteration", i
		#print "fmeanscores2d are",fmeanscores2d
		for ff in fmeanscores3d[i]:
			#print "therefore ff in fmeanscores2d[i] is", ff 
			iterscores3d.append( fmeanscores3d[i][ff][-1] )
		
		itermeanscore3d = sum(iterscores3d)/len( iterscores3d )
		itermeanscores3d[i] = itermeanscore3d
	
		ys3 = []
		lines3=[]
		score3dindxs = [s for s in itermeanscores3d]
		score3dindxs.sort()
		
		for s in score3dindxs: 
			y3 = itermeanscores3d[s]
			ys3.append( y )
			#print "appending this y", itermeanscores2d[s]
			line3 = str(s+1) + '\t' + str(y3) + '\n'
			lines3.append(line3)
	
		#if plots:
		fs = open( options.path +'/scores3d_' + str(i+1) + '.txt','w')
		fs.writelines( lines3 )
		fs.close()
		
		
		if i > 0:
			if itermeanscores2d[i] == itermeanscores2d[i-1]:
				print "The meanscore2d for two consecutive iterations is the same, suggesting the algorithm has converged."
				sys.exit()
			
			try:
				if itermeanerrors[i] == itermeanerrors[i-1]:
					print "The meanerror for two consecutive iterations is the same, suggesting the algorithm has converged."
					sys.exit()
			except:
				pass
				
			if itermeanfscs[i] == itermeanfscs[i-1]:
				print "The fscsarea for two consecutive iterations is the same, suggesting the algorithm has converged."
				sys.exit()
				
			if itermeanscores3d[i] == itermeanscores3d[i-1]:
				print "The meanscore3d for two consecutive iterations is the same, suggesting the algorithm has converged."
				sys.exit()
			
			
			
			difscore = math.fabs( float( itermeanscores2d[i] ) - float(itermeanscores2d[i-1]) )
			if float( difscore ) < 0.000001:
				print "In iter %d/%d global score difference with previous iteration is smaller than 0.000001, see: %f \nExiting program." %( i+1, options.iter, difscore )
				sys.exit()
				
			try:
				diferror = math.fabs( float( itermeanerrors[i] ) - float(itermeanerrors[i-1]) )
				if float( difscore ) < 0.000001:
					print "In iter %d/%d global error difference with previous iteration is smaller than 0.000001, see: %f \nExiting program." %( i+1, options.iter, diferror )
					sys.exit()
			except:
				pass
								
			diffscarea = math.fabs( float( itermeanfscs[i] ) - float(itermeanfscs[i-1]) )
			if float( difscore ) < 0.000001:
				print "In iter %d/%d global fscarea difference with previous iteration is smaller than 0.000001, see: %f \nExiting program." %( i+1, options.iter, diffscarea )
				sys.exit()
			
			difscore3d = math.fabs( float( itermeanscores3d[i] ) - float(itermeanscores3d[i-1]) )
			if float( difscore3d ) < 0.000001:
				print "In iter %d/%d global score3dipetResult difference with previous iteration is smaller than 0.000001, see: %f \nExiting program." %( i+1, options.iter, difscore )
				sys.exit()
			
			
	E2end(logger)
	
	return()
	

def recentertilts( options, reprojections, originalseries, angles, it ):
	
	print "\nRecentering subtiltseries", originalseries
	n = EMUtil.get_image_count( originalseries )
	print "Which contains these many images:", n
	
	hdr1 = EMData( originalseries, 1, True )
	nx = hdr1['nx']
	ny = hdr1['ny']
	
	originaltilts = {}
	recenteredtilts = {}
	for i in range(n):
		img = EMData( originalseries, i )
		angle = img['spt_tiltangle']
		originaltilts.update({ angle:img })
	
	kkk=0
	txs=[]
	tys=[]
	drs=[]
	scores2d=[]
	errors=[]
	#print "\n\n\nIn recentertilts, tilt angles are", originaltilts
	
	for angle in angles:
		print kkk, angle
		#print "Recentering image", kkk
		oimg = originaltilts[angle]
		rpimg = reprojections[angle]

		#print "size of rpimg is", rpimg['nx'],rpimg['ny']
		#print "size of oimg is", oimg['nx'],oimg['ny']

		oimgpp = oimg.copy()
		oimgpp = preprocImg( oimgpp, options )
		
		rpimgpp = rpimg.copy()
		rpimgpp = preprocImg( rpimgpp, options )
		
		#print "Size of rpimgpp is", rpimgpp['nx'],rpimgpp['ny']
		#print "Size of oimgpp is", oimgpp['nx'],oimgpp['ny']
		
		rpimgpp.process_inplace('filter.matchto',{'to':oimgpp})
		
		ccf = rpimgpp.calc_ccf( oimgpp )
		ccf.process_inplace("xform.phaseorigin.tocorner") 
		ccf.process_inplace('normalize')
		
		maxloc = ccf.calc_max_location()
		
		score = ccf.get_value_at( maxloc[0], maxloc[1] )
		scores2d.append( score )
		
		tx = -1 * ( nx/2.0 - maxloc[0])
		ty = -1 * ( ny/2.0 - maxloc[1])
		
		dr = math.sqrt( tx*tx+ty*ty )
		
		txs.append(tx)
		tys.append(ty)
		drs.append(dr)
		
		print "\nIn iteration %d, for image %d, x,y translations are x=%f, y=%f" % ( it, kkk, tx, ty )

		try:
			print "whereas original images were off by tx", oimg['spt_txerror'], oimg['spt_tyerror']
			
			print "\nThe header is", oimg.get_attr_dict()
			
			xerror = float(tx) + float(oimg['spt_txerror'])
			yerror = float(ty) + float(oimg['spt_tyerror'])
			#error =  math.sqrt( xerror*xerror + yerror*yerror )/2.0
			error = ( math.fabs(xerror) + math.fabs(yerror) )/2.0
			
			if error:
				errors.append(error)
			
			print "therfore the error is", error
			
		except:
			pass
		
		oimgt = oimg.copy()
		
		oimgt.translate( tx, ty, 0)
		
		oimgt['sptisrtx'] = tx
		oimgt['sptisrty'] = ty
		oimgt['sptisrdr'] = dr
		
		recenteredtilts.update({ angle:oimgt })
		
		if options.saveali or int(options.iter - it) <= 1:
			outimg = options.path + '/' + os.path.basename( originalseries ).replace('.hdf','_ISR.hdf')
			oimgt.write_image( outimg, kkk )
		
		kkk+=1
	
	txsavg = sum(txs)/len(txs)
	tysavg = sum(tys)/len(tys)
	drsavg = sum(drs)/len(drs)
	
	statsline = "For file: " + os.path.basename(originalseries) + ", mean tx=" + str(txsavg) + ", mean ty=" + str(txsavg) + ", mean dr=" + str(drsavg)	+ '\n'
	
	
	#print "sent scores2d are", scores2d
	#print "sent errors are", errors
	if scores2d == errors:
		print "ERROR: sent scores2d and errores are the same"
		sys.exit()
	print 
	
	return [ recenteredtilts, statsline, scores2d, errors ]

	
def calcangles( f ):

	nf = EMUtil.get_image_count( f )
	angles = []
	angle = None

	for i in range( nf ):
		imghdr = EMData( f, i, True )
		try:
			angle = imghdr['spt_tiltangle']
			angles.append( angle ) 
		except:
			print """\n(e2spt_isr.py)(calcangles) ERROR: image %d in stack %s lacking 
				spt_tiltangle parameter in header. If tilt angle information is not in the 
				header of the data, supply it via --tltfile.""" %( i, f )
			sys.exit()	

	angles.sort()
		
	return angles


def getangles( options ):
	print "\n(e2spt_fillwedge.py)(getangles)"
	angles = []
	
	#print "Reading tlt file", options.tltfile
	f = open( options.tltfile, 'r' )
	lines = f.readlines()
	f.close()
	
	for line in lines:
		line = line.replace('\t','').replace('\n','')
		if line:
			angles.append( float(line) )
		#print "Found angle", line
		
	if options.verbose > 9:
		pass
		#print "\n(e2spt_ctf.py)(getangles) angles are", angles

	return angles


def writeparamtoheader( f, angs, param ):
	print "\n(e2spt_fillwedge.py)(writeparamtoheader)"
	n = EMUtil.get_image_count(f)
	print "Writing parameter %s to header of file %s" % (param, f)
	print "With these many images and angles", n, len(angs)
	for i in range(n):
		print "Working on image", i
		imghdr = EMData( f, i, True)
		imghdr[param]=angs[i]
		imghdr.write_image( f, i, EMUtil.ImageType.IMAGE_HDF, True )
		
	return 1
	
	
def writeanglestofile( options, angsf ):
	lines = []
	for a in angsf:
		line = str(angsf)+'\n'
		lines.append(line)
	finalangsfile = options.path + '/completeangles.tlt'
	gg=open(finalangsfile,'w')
	gg.writelines( lines )
	gg.close()
	return
	
	
def calcWeight( angle, it, options ):

	weight = 1.0
	#if float(angle) < float(lowmostangle) or float(angle) > float(upmostangle):
	if math.fabs( angle ) == 90.0:
		complement = 1.0 - math.fabs( math.cos( math.radians(89.99) ) )
		weight = math.fabs( (it+1) * math.cos( math.radians(89.99) ) / float(options.iter) ) + ( float(it)/float(options.iter)) * complement
	else:
		complement = 1.0 - math.fabs( math.cos( math.radians(angle) ) )
		weight = math.fabs( (it+1) * math.cos( math.radians(angle) ) / float(options.iter) ) + (float(it)/float(options.iter)) * complement 

	print "Something was weighed?!!!!" #??
	
	return weight


def reprojectvolume( options, vol, originalseries, originalangles ):
		
	originalvolsize = max( vol['nx'], vol['ny'], vol['nz'])
	#expandedvolsize = originalvolsize + 2
	
	vol.process_inplace('normalize.edgemean')
	
	#vol = clip3D( vol, expandedvolsize )
	
	vol = preprocImg( vol, options )
	
	reprjs = {}
	jkl = 0
	for angle in originalangles:
		rp = genprojection( vol, angle, options.tiltaxis )
		rp = clip2D( rp, originalvolsize )
		reprjs.update({ angle : rp })
		
		if options.saveali:
			rppath = options.path + '/' + os.path.basename( originalseries ).replace('.hdf','_ISRreprj.hdf')
			rp.write_image( rppath, jkl )
		jkl+=1
		
	return reprjs
	

def genprojection( vol, angle, tiltaxis ):

	t = Transform({'type':'eman','az':90,'alt':angle,'phi':-90})
		
	if tiltaxis == 'x':
		t = Transform({'type':'eman','alt':angle})
	
	prj = vol.project("standard",t)
	#prj.process_inplace('normalize.edgemean')
	prj['spt_tiltangle'] = angle

	return prj


def clip3D( vol, sizex, sizey=0, sizez=0 ):
	
	if not sizey:
		sizey=sizex
	
	if not sizez:
		sizez=sizex
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - sizex)/2, (2*volyc - sizey)/2, (2*volzc - sizez)/2, sizex , sizey , sizez)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol


def clip2D( img, size ):
	
	imgxc = img['nx']/2
	imgyc = img['ny']/2
	#imgzc = img['nz']/2
	
	Rimg =  Region( (2*imgxc - size)/2, (2*imgyc - size)/2, 0, size , size , 1)
	img.clip_inplace( Rimg )
	#img.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return img
	

def preprocImg( img, options, filling=False ):
	
	#print "Inside preprocImg, mask is", options.mask, type(options.mask), len(options.mask)

	
	"""
	if filling and options.preprocessfilling and options.preprocessfilling != 'None' and options.preprocessfilling != 'none':
		preprocessfilling=''
		try:
			preprocessfilling=parsemodopt(options.preprocessfilling)
		except:
			pass	
		img.process_inplace( preprocessfilling[0], preprocessfilling[1] )
	"""
	if options.preprocess and options.preprocess != 'None' and options.preprocess != 'none': 
		preprocess=''
		
		img.process_inplace( options.preprocess[0], options.preprocess[1] )
		
	"""
	if filling and options.maskfilling and options.maskfilling != 'None' and options.maskfilling != 'none':
		maskfilling=''
		try:
			maskfilling=parsemodopt(options.maskfilling)
		except:
			pass
		img.process_inplace( maskfilling[0], maskfilling[1] )
	"""
	if options.mask and options.mask != 'None' and options.mask != 'none':
		mask=''
		
		img.process_inplace( options.mask[0], options.mask[1] )
	
	"""
	if filling and options.thresholdfilling and options.thresholdfilling != 'None' and options.thresholdfilling != 'none':
		thresholdfilling=''
		try:
			thresholdfilling=parsemodopt(options.thresholdfilling)
		except:
			print "Failed to parse threshold"	
		print "Parsed threshold is", thresholdfilling
		img.process_inplace( thresholdfilling[0], thresholdfilling[1] )	
	"""
	if options.threshold and options.threshold != 'None' and options.threshold != 'none': 
		
		img.process_inplace( options.threshold[0], options.threshold[1] )

	"""
	if filling and options.highpassfilling and options.highpassfilling != 'None' and options.highpassfilling != 'none':
		highpassfilling=''
		try:
			highpassfilling=parsemodopt(options.highpassfilling)
		except:
			pass
		img.process_inplace( highpassfilling[0], highpassfilling[1] )
	"""
	if options.highpass and options.highpass != 'None' and options.highpass != 'none': 
	
		img.process_inplace( options.highpass[0], options.highpass[1] )
	
	"""
	if filling and options.lowpassfilling  and options.lowpassfilling != 'None' and options.lowpassfilling != 'none':
		lowpassfilling=''
		try:
			lowpassfilling=parsemodopt(options.lowpassfilling)
		except:
			pass
		img.process_inplace( lowpassfilling[0], lowpassfilling[1] )
	"""
	if options.lowpass and options.lowpass != 'None' and options.lowpass != 'none': 
	
		img.process_inplace( options.lowpass[0], options.lowpass[1] )
	
	return img


def makevol( options, originalseriesfile, newseries, it, originalangles, originalsize, writevols = 1, initialfsc = 0 ):
		
	#print "Newseries received in makevol and its type are", type(newseries)
	#print newseries
	#sys.exit()
	mode='gauss_2'
	if options.reconstructor:
		if len(options.reconstructor) > 1:
			if 'mode' in options.reconstructor[-1]:
				try:
					if options.reconstructor[-1]['mode'] != 'gauss_2':
						mode = options.reconstructor[-1]['mode']
					else:
						pass
				except:
					pass
	
	box = originalsize
	if options.pad3d:
		if options.pad2d:
			if options.pad3d > options.pad2d:
				box = int( box*options.pad3d )
			else:
				box = int(box*options.pad2d)
		else:
			box = int(box*options.pad3d)		
	elif options.pad2d:
		box = int(box*options.pad2d)
	
		#if options.verbose > 5:
			#print "\nPPPPPPPPPPPPPPP\n\n(e2spt_isr.py)(makevol) Because there's options.pad=%f,then newsize=%d" % ( options.pad3d, newsize )
	
	#if options.verbose > 9:
	#print "\n(e2spt_isr.py)(makevol) Mode for reconstructor is", mode
	#print "Setting up reconstructor; options.reconstructor is", options.reconstructor
	#print "For volume of size",box,box,box
		
	#rEven = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})
	
		
	r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':False,'mode':mode})
	r.setup()
	
	#print "in isr reconstructor is", options.reconstructor
	#print "r is", r
	#print "and box was", box
	
	axis = 'y'
	if options.tiltaxis:
		axis = options.tiltaxis
	
	angles = []
	
	#print "\n++++++++++++++++++++++++++++++\n(makevol) file=%s, lowmostRAWangle=%f, upmostRAWangle=%f" % ( originalseriesfile, min(originalangles), max(originalangles) )
	#print "However, adding projections, the total number of images and angles will be", len(newseries),len(originalangles)
	#print "\n++++++++++++++++++++++++++++++\n"
	
	apix = 1.0
	
	imageslist = []
	for tiltangle in originalangles:
		#print "Angle and type are", tiltangle, type(tiltangle)

		img = newseries[tiltangle]
		
		#try:
		#	if img['spt_tiltangle']:
		#		pass
		#except:		
		
		img['xform.projection'] = Transform({'type':'eman','az':90,'alt': float(tiltangle),'phi':-90})
		
		img2send = img.copy()
		
		imageslist.append( img2send )
		
		#print "size of image appended to imageslist", img2send['nx'],img2send['ny']
		
		apix = img['apix_x']
		
		angle = None
		if options.verbose > 9:
			print  "\n(e2spt_isr.py)(makevol) processing image %d in stack %s and ITERATION %d" % ( i, originalseriesfile, it )
			print "For image %d, mean=%f, min=%f, max=%f, std=%f" % ( i, img['mean'], img['minimum'],img['maximum'],img['sigma'] )
		
		if options.pad2d and float(options.pad2d) > 1.0:
			box2d = img['nx'] * options.pad2d
			img = clip2D( img, box2d )
			#print "\nPadding prj to newsize", newsize
		
		try:
			angle = img['spt_tiltangle'] 
		except:
			print "\nWARNING: 'spt_tiltangle' parameter not found in the header of image %d in file %s" % ( i, f )
			#sys.exit()
			angle=tiltangle
			print "Using this tiltangle", tiltangle
		
		try:
			axis = img['spt_tiltaxis']
		except:
			if options.verbose > 9:
				print """\n(e2spt_isr.py)(makevol) WARNING: 
				No spt_tiltaxis or sptsim_tiltaxis found in header. Default used.""", options.tiltaxis, axis
		
		if angle != None:
			#angles.append( angle )
		
			t = Transform({'type':'eman','az':90,'alt':angle,'phi':-90})
		
			if axis == 'x':
				t = Transform({'type':'eman','alt':angle})
			
			imgp = r.preprocess_slice( img, t)
			
			weight = 1.0
			
			if options.verbose > 9:
				print "\n(makevol) ITER=%d, tiltangle=%f, weight=%f" % ( it, angle, weight )
				print "Inserted IMAGE with this tranform", t
			
			#print "transform used for WHOLE 3D", t
			
			r.insert_slice( imgp , t , weight )
	
	rec = r.finish(True)

	rec['apix_x']=apix
	rec['apix_y']=apix
	rec['apix_z']=apix
	rec['origin_x']=0
	rec['origin_y']=0
	rec['origin_z']=0
	
	#angles.sort()
	#if options.verbose > 9:
	#	print "\n(e2spt_fillwedge.py)(makevol) Sorted angles to write to header are", angles
	
	rec['spt_tiltangles'] = originalangles

	recxc = rec['nx']/2
	recyc = rec['ny']/2
	reczc = rec['nz']/2

	#R2 =  Region( (2*recxc - originalsize)/2, (2*recyc - originalsize)/2, (2*reczc - originalsize)/2, originalsize , originalsize , originalsize)
	#rec.clip_inplace( R2 )
	
	outx=outy=outz=originalsize
	if options.outxsize:
		outx = options.outxsize
	if options.outysize:
		outy = options.outysize
	if options.outzsize:
		outz = options.outzsize
	
	rec = clip3D( rec, outx, outy, outz )	
	
	rec.process_inplace( 'normalize' )
	
	if options.verbose > 9:
		print "\n(e2spt_isr.py)(makevol) Reconstructed volume for file", originalseriesfile
	
	volfile = ''
	if options.savevols and writevols:
		volfile = os.path.basename( originalseriesfile ).replace('.hdf','_ISR3D.hdf')
		if options.path not in volfile:
			volfile = options.path + '/' + volfile

		rec.write_image( volfile, 0 )
	
	#print "size of images in imageslist to send is", imageslist[0]['nx'],imageslist[-1]['ny']
	
	ret = genOddAndEvenVols( options, originalseriesfile, imageslist )
	volOdd = ret[1]
	volEven = ret[0]
	score3d = ret[2]
	
	#volOdd.write_image( options.path + '/' + originalseriesfile.replace('.hdf','_ISR3D_ODD.hdf') )
	#volEven.write_image( options.path + '/' + originalseriesfile.replace('.hdf','_ISR3D_EVEN.hdf') )
	
	filebasetoplotfsc = originalseriesfile
	if it == 0:
		print "There SHOULD be intialfsc", initialfsc
	if initialfsc:
		print "There is intialfsc", initialfsc
		filebasetoplotfsc = originalseriesfile.replace('.hdf', '_INITIAL.hdf')

	
	retfsc = fscOddVsEven( options, filebasetoplotfsc, volOdd, volEven )
	fsc = retfsc[0]
	fscarea = retfsc[1]
	aliscore = []
	
	return [ rec, volfile, fscarea, score3d ]

	
if '__main__' == __name__:
	main()
	

	