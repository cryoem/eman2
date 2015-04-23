#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 9/Mar/2015
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
#

from EMAN2 import *
import os
import sys
import time
import numpy
import pylab
#from operator import itemgetter
from matplotlib.ticker import MaxNLocator
from pylab import figure, show	
import matplotlib.pyplot as plt

#import scipy.optimize			 

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2spt_classaverage.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="""Volume or stack of volumes to be aligned to --ref; that is, 
													volumes in --input will be 'moved' (rotated and translated) to find the best fit with --ref. 
													The format MUST be '.hdf' or '.mrc' """, default=None)
	parser.add_argument("--path", type=str, help="Results directory. If not specified, defaults to e2sptfscs/", default='e2sptfscs')

	parser.add_argument("--ref", type=str, help="Volume that will be 'static' (the 'reference' to which volumes in --input will be aligned to). The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--output", type=str, help="""Name for the .txt file that will contain the FSC data. If not specified, a default name will be used.""", default=None)


	parser.add_argument("--symali", type=str, default='c1', help = """Will pass the value to --sym for alignment with e2spt_classaverage.py""")

	parser.add_argument("--nocolor",action='store_true',default=False,help="""Turns the ouput png(s) into grey scale figures. Instead of using different colors to distinguish between various curves on the same plot, this option will have the program automatically use different markers in black and white for each curve.""")
	
	parser.add_argument("--symmap", type=str, default='c1', help = """Will symmetrize --ref AND pass --sym to e2spt_classaverage.py""")
	
	parser.add_argument("--maskali",type=str,help="""Mask processor applied to particles 
		before alignment. Default is mask.sharp:outer_radius=-2""", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--maskfsc",type=str,help="""Mask processor applied to particles 
		before fsc computation. Default is mask.sharp:outer_radius=-2""", default=None)
	
	
	parser.add_argument("--search", type=int,default=8,help=""""During COARSE alignment
		translational search in X, Y and Z, in pixels. Only works when --radius is provided.
		Otherwise, search parameters are provided with the aligner, through --align.""")
	
	parser.add_argument("--searchfine", type=int,default=2,help=""""During FINE alignment
		translational search in X, Y and Z, in pixels. Only works when --radius is provided.
		Otherwise, search parameters are provided with the aligner, through --falign.""")
			
	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles before alignment. 
													Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. 
													If you want to turn this option off specify \'None\'""", default="normalize.mask")

	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=True)

	parser.add_argument("--mirror",action="store_true", help="""If set, it will generate a mirrored version of --ref and align --input against it.=; then FSCs will be computed. 
															This will be done IN ADDITION to aligment and FSC computation of the alignment of --input against the original, unmirrored --ref.""",default=False)


	parser.add_argument("--preprocess",type=str,default='',help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.")
	parser.add_argument("--preprocessfine",type=str,default='',help="Any processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.")
	
	parser.add_argument("--lowpass",type=str,default='',help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.")
	parser.add_argument("--lowpassfine",type=str,default='',help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.")

	parser.add_argument("--highpass",type=str,default='',help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.")
	parser.add_argument("--highpassfine",type=str,default='',help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.")

	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkfine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")

	parser.add_argument("--threshold",type=str,default='',help="""A threshold applied to 
		the subvolumes after normalization. 
		For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels 
		equal 0, so that they do not contribute to the correlation score.""")

	parser.add_argument("--nocenterofmass", default=False, action="store_true", help="""Disable Centering 
		of mass of the subtomogram every iteration.""")

	parser.add_argument("--npeakstorefine", type=int, help="The number of best 'coarse peaks' from 'coarse alignment' to refine in search for the best final alignment. Default=4.", default=4)
	
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average.", default="rotate_translate_3d:search=6:delta=12:dphi=12")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	
	parser.add_argument("--radius", type=float, help="""Hydrodynamic radius of the particle in Angstroms. 
													This will be used to automatically calculate the angular steps to use in search of the best alignment.
													Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels.
													Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that""", default=0)
	
	parser.add_argument("--falign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine_3d:search=2:delta=3:range=12", default="refine_3d:search=2:delta=3:range=9")
	parser.add_argument("--faligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:2")
	
	parser.add_argument("--apix", type=float, help="Provide --apix for automatic FSC calculation if you supply --plotonly and no volumes for --input and --ref, or if the apix of these is wrong.", default=1.0)
	parser.add_argument("--boxsize", type=float, help="(Probably not needed for anything)", default=0)

	parser.add_argument("--maxres", type=float, help="How far in resolution to extend the FSC curve on the x axis; for example, to see up to 20anstroms, provide --maxres=1.0. Default=15", default=1.0)
	parser.add_argument("--cutoff", type=str, help="Comma separated values of cutoff thresholds to plot as horizontal lines. Default=0.5, to turn of supply 'None'. ", default='0.5')
	
	parser.add_argument("--smooth",action="store_true", help="""Smooth out FSC curves 
		by taking the average of a low value with a subsequent maxima.""", default=False)
	
	parser.add_argument("--smooththresh",type=float, help="""If --smooth is provided
		the curve will be smoothed only up to this resolution. Default is 100.""", default=100)

	#parser.add_argument("--fit",action="store_true", help="Smooth out FSC curves.", default=False)
	
	parser.add_argument("--polydegree",type=int, help="Degree of the polynomial to fit.", default=None)


	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID.",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	#parser.add_argument("--plotonly",action="store_true", help="Assumes vol1 and vol2 are already aligned and fsc files will be provided through --curves; thus skips alignment and fsc curve generation", default=False)
	
	parser.add_argument("--fsconly",action="store_true", help="""Assumes --input and --ref 
		are already aligned with respect to each other and thus skips alignment""", default=False)
		
	parser.add_argument("--plotonly",type=str, help="""FSC curves to plot in separate plots. 
		Skips alignment and fsc curve generation. Provide .txt. files separated by commas 
		--plotonly=file1.txt,file2.txt,file3.txt etc...""", default=None)
	parser.add_argument("--singleplot",action="store_true",help="It --plotonly provided, all FSC curves will be on the same plot/figure", default=False)
		
	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)

	if options.maskfsc: 
		options.maskfsc=parsemodopt(options.maskfsc)

	if options.cutoff and options.cutoff != 'None' and options.cutoff != 'none':
		options.cutoff = options.cutoff.split(',')
		print "Options.cutoff is", options.cutoff
	else:
		options.cutoff = None

	if not options.output and not options.plotonly:
		print "ERROR: Unless you provide .txt files through --plotonly, you must specify an --output in .txt format."
	elif options.output and not options.plotonly:
		if '.txt' not in options.output:
			print "ERROR: Output must be in .txt format"
			sys.exit()
		
	findir=os.listdir(os.getcwd())

	#if options.path not in findir:
	#	os.system('mkdir ' + options.path)
	
	#print "Befere options are of type", type(options)
	#print "\n\n\nand are", options
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptres')

	print '\n\nafter making path, options.path is', options.path
	
	if options.input:
		hdr = EMData(options.input,0,True)
		apix = hdr['apix_x']
	
	if options.apix:
		apix = float (options.apix)
		
	#print "Returned options are of type", type(options)
	#print "\n\n\nand are", options
	
	if options.plotonly:
		curves = options.plotonly
		curves = curves.split(',')
		
		if not options.singleplot:
			for curve in curves:
				print "Found this curve to plot", curve
				fscplotter([curve],options,apix)
				
		elif options.singleplot:
			fscplotter(curves,options,apix)
		
		print "Done plotting"
		sys.exit()
		
	elif options.fsconly:
		getfscs(options,apix)
	
		print "Done calculating FSCs and plotting them."
		sys.exit()
	
	elif options.align:
		print '\n\nI will align'
		if options.symmap and options.symmap is not 'c1' and options.symmap is not 'C1':		
			ref = EMData(options.ref,0)
			ref = symmetrize(ref,options)
			options.ref = options.path + '/' + os.path.basename(options.ref).replace('.','_' + options.sym + '.')
			ref.write_image(options.ref,0)
		
		ptclali = alignment(options)
		
		print '\n\nthe returned path for ali ptcl is', ptclali
		inputbackup = options.input
		options.input = ptclali
		
		print '\n\nwill get fsc'
		getfscs(options,apix)

		if options.mirror:
			options.input = inputbackup
			ref = EMData(options.ref,0)
			t=Transform({'type':'eman','mirror':True})
			
			options.ref = os.path.basename(options.ref).replace('.','_mirror.')
			if options.path not in options.ref:
				options.ref = options.path + '/' + options.ref
		
			ref.transform(t)
			ref.write_image(options.ref,0)
			
			ptclalimirror = alignment(options)
		
			options.input = ptclalimirror
			getfscs(options,apix)
		
	E2end(logger)
	return


def getfscs(options,apix):
	
	options.input
	print "\n inside getfscs options.input is", options.input
	print 'and the current directory is', os.getcwd()
	fyle = 'alignment/' + options.input.split('/')[-1]
	
	
	if options.path not in fyle:
		fyle = options.path + '/' + fyle
	
	
	current = os.getcwd()
	
	if current not in fyle:
		fyle = current + '/' + fyle
	
	
	print "\n\n\n\n\nThe fyle to get image count from is", fyle
	print "Whereas current dir is", os.getcwd()
	
	n = EMUtil.get_image_count( fyle )
	
	path = options.path
	fscs = []
	fscsm = []
	for i in range(n):
		ptcl = EMData( fyle ,i)
		ref = EMData(options.ref,0)
		fscfilename = path + '/' + options.output.replace('.txt','_' + str(i).zfill(len(str(n))) + '.txt')
		if options.symmap and options.symmap is not 'c1' and options.symmap is not 'C1':
			ptcl = symmetrize(ptcl,options) 
			ref = symmetrize(ref,options)
			#fscfilename = options.output.replace('.txt','_' + str(i).zfill(len(str(n))) + '_' + options.sym + '.txt')
			fscfilename = fscfilename.replace('.txt', '_' + options.symmap + '.txt')
		
		if options.maskfsc:
		
			'''
			Make the mask first 
			'''
			mask=EMData( int(ptcl["nx"]), int(ptcl["ny"]), int(ptcl["nz"]) )
			mask.to_one()
			mask.process_inplace(options.maskfsc[0],options.maskfsc[1])
			
			ptcl.mult(mask)
			ref.mult(mask)
			
		calcfsc(ptcl,ref,fscfilename,options)
		
		if not options.singleplot:
			fscplotter([fscfilename],options,apix)
		else:
			fscs.append(fscfilename)
	
		if options.mirror:
			t = Transform({'type':'eman','mirror':True})
			refm = ref.copy()
			refm.transform(t)
			fscmfilename = fscfilename.replace('.txt','_VSmirror.txt')
			print "The mirror fsc file is", fscfilename
			
			calcfsc(ptcl,refm,fscmfilename,options)
			
			if not options.singleplot:
				fscplotter([fscmfilename],options,apix)
			else:
				fscsm.append(fscmfilename)
	
	if options.singleplot:
		fscplotter(fscs,options,apix)
		fscplotter(fscsm,options,apix)
	return


def alignment(options):
	aligner=options.align
	aligner=aligner.split(':')
	newaligner=''
		
	if 'sym='+options.symali not in aligner:
		
		for element in aligner:
			newElement = element
			if 'sym' in element:
				newElement='sym=' +options.symali
			newaligner=newaligner + newElement + ":"
		
		if newaligner[-1] == ':':
			newaligner = newaligner[:-1]
			options.align=newaligner
	
	print "\n\n\n\n\nTHE fixed aligner is", newaligner, options.align
	
	alivolfile = os.path.basename(options.input).split('.')[0] + '_VS_' + os.path.basename(options.ref).split('.')[0] + '.hdf'
	#print "alivolfile is", alivolfile
	#print "Path is", options.path
	#print "input is", options.input
	#print "ref is", options.ref
	#print "npeakstorefine is", options.npeakstorefine
	#print "verbose is", options.verbose
	#print "maskali is", options.maskali
	#print "lowpass is", options.lowpass
	#print "\n\nalign and its type are", options.align, type(options.align)
	#print "\n\nfalign is type are", options.falign, type(options.falign)
	
	alicmd = 'cd ' + options.path + ' && e2spt_classaverage.py --search=' + str(options.search) + ' --searchfine=' + str(options.searchfine) + ' --path=alignment --input=../' + str(options.input) + ' --output=' + str(alivolfile) + ' --ref=../' + str(options.ref) + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose) + ' --mask=' + str(options.maskali) + ' --lowpass=' + str(options.lowpass) + ' --parallel=' + str(options.parallel) + ' --aligncmp=' + str(options.aligncmp) + ' --faligncmp=' + str(options.faligncmp) + ' --shrink=' + str(options.shrink) + ' --shrinkfine=' + str(options.shrinkfine) + ' --threshold=' + str(options.threshold) + ' --saveali' + ' --normproc=' + str(options.normproc) + ' --sym=' + str(options.symali) + ' --breaksym'
	
	if options.nocenterofmass:
		alicmd += ' --nocenterofmass'
	
	if options.radius:
		alicmd += ' --radius=' + str(options.radius)
	
	else:
		alicmd += ' --align=' + str(options.align) + ' --falign=' + str(options.falign)
	
	
	print "\n\n\n\n\n\nThe alig command to run is", alicmd
			
	#os.system(alicmd)
	
	p=subprocess.Popen( alicmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	#for line in iter(p.stdout.readline, ''):
	#	print line.replace('\n','')

	text=p.communicate()	
	p.stdout.close()
	
	print "BEFORE ALIGNMENT, input was", options.input
	print "BEFORE ALIGNMENT, ref was", options.ref
	print "Returning this from ALIGNMENT", alivolfile
	
	return os.getcwd() + '/' + options.path + '/alignment/' + alivolfile


def calcfsc(v1,v2,fscfilename,options):
	fsc = v1.calc_fourier_shell_correlation(v2)
	third = len(fsc)/3
	xaxis = fsc[0:third]
	fsc = fsc[third:2*third]
	apix=v1['apix_x']
	if options.apix:
		apix=options.apix
	saxis = [x/apix for x in xaxis]
	Util.save_data(saxis[1],saxis[1]-saxis[0],fsc[1:], fscfilename)
	return fsc	
		

def symmetrize(vol,options):
	sym = options.symmap
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	volsym=vol.copy()
	for i in range(1,nsym):
		dc=vol.copy()
		t=xf.get_sym(sym,i)
		dc.transform(t)
		volsym.add(dc)
	volsym.mult(1.0/nsym)	
	return(volsym)


"""
SIGMOID FITTING
"""
def sigmoidfit(x,values):		
	import numpy as np
	import scipy
	
	def sigmoid(p,x):
		x0,y0,c,k=p
		y = c / (1 + np.exp(-k*(x-x0))) + y0
		return y
	
	def residuals(p,x,y):
		return y - sigmoid(p,x)
	
	def resize(arr,lower=0.0,upper=1.0):
		arr=arr.copy()
		if lower>upper: 
			lower,upper=upper,lower
		arr -= arr.min()
		arr *= (upper-lower)/arr.max()
		arr += lower
		return arr
	
	xsig = np.array(x,dtype='float')
	ysig = np.array(values,dtype='float')
	
	xsig = resize(-xsig,lower=0.3)
	ysig = resize(ysig,lower=0.3)
	print(xsig)
	print(ysig)
	p_guess=(np.median(xsig),np.median(ysig),1.0,1.0)
	p, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,p_guess,args=(xsig,ysig),full_output=1,warning=True)  
	
	xsig0,ysig0,c,k=p

	x0 = {x0}
	y0 = {y0}
	c = {c}
	k = {k}
	format(xsig0=xsig0,ysig0=ysig0,c=c,k=k)

	xsigp = np.linspace(0, 1.1, 1500)
	pxsigp=sigmoid(p,xsigp)

	print "\n\n\nPRESUMABLY THE SMOOTH CURVE IS pxsigp, with these many elements", len(pxsigp)
	print "Which are", pxsigp
	print "\n\n\n"

	#plt.plot(xsig, ysig, '.', xsigp, pxsigp, '-')
	#plt.savefig(options.output)

	return(xsig,ysig)
	

def maxima(xaxis,yaxis,smooththresh):
	print "\n\nI have entered maxima!!!!!!!!!!!!!!!!!!!!!\n\n"
	print "At this point xaxis and yaxis len are", len(xaxis), len(yaxis)
	#xaxis=list(set(xaxis))
	#yaxis=list(set(yaxis))	

	#xes=[xaxis[0]]
	#ymaxes=[yaxis[0]]
	
	#xaxis=list(set(xaxis))
	xes=[]
	ymaxes=[]
	
	
	for i in xrange(0,len(yaxis) -1):
		val=yaxis[i]
		#print 'currrent value is', val
		#print 'and next is', yaxis[i+1]
		#print "options.smooththresh is", smooththresh
		if val < max(yaxis[i+1:]) and 1.0/xaxis[i+1] > smooththresh:
			val = ( val+ max(yaxis[i+1:]) )/ 2.0
			print '\nNew max smoothing value is', val

		if val > min(yaxis[i+1:]) and 1.0/xaxis[i+1] < smooththresh:
			val = val/ 2.0
			print '\nNew min smoothing value is', val
		
		#print "Therfore final val is", val
			
		
		ymaxes.append(val)
		xes.append(i)
	
	#ymaxes.append(yamxes[-2])
	ymaxes.append(ymaxes[-1])
			
	#xes.append(xes[-2])
	xes.append(i+1)
	
	'''
	for i in xrange(0,len(yaxis) - 1 - 1 ):
		aux = 0 
		for j in xrange(i+1,len(yaxis) - 1 - 1 ):
			if yaxis[j] > yaxis[i]:
				aux=0
				print "Because a downstream value is higher, I will break the loop, see", yaxis[j],yaxis[i]
				break
			else:
				aux=1
		if aux==1:
			ymaxes.append(yaxis[i])
		else:
			ymaxes.append(yaxis[j])
		xes.append(xaxis[i])
		#print "I have appended this box, value", sizes[i], vals[i]

	ymaxes.append(ymaxes[-2])
	ymaxes.append(ymaxes[-1])
			
	xes.append(xes[-2])
	xes.append(xes[-1])
	
	#minnonconvex=Util.nonconvex(valsmin,0)
	#print "In minima, the yminnonconvex to return is", minnonconvex
	'''
	
	print "\nExiting maxima, xes len is", len(xes)
	print "Exiting maxima, ymaxes len is", len(ymaxes)
	print "\nExiting maxima, xes are", xes
	print "Exiting maxima, ymaxes are", ymaxes
	
	#x=list(set(x))
	return(xes,ymaxes)
	
	
		
def fscplotter(fscs,options,apix=0.0):
	fig = figure()

	#from itertools import product
	markers=["*","o","x","z","M"]
	#colors=["k","r","b","g","m","y","c"]
	#colorsANDmarkers = [a+b for a,b in product(colors,markers)]
	
	import colorsys
	import itertools

	
	N = len(fscs)
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
	
	kont=0
	plot_name = ''
	for fscoutputname in fscs:
	
		print "\nI found THIS FSC file and will thus plot it", fscoutputname

		f= open(fscoutputname,'r')
		lines = f.readlines()

		if not options.boxsize:
			options.boxsize = len(lines)


		values = []
		inversefreqs = []
		inversefreqslabels = []
		
		newlines=[]
		firstline=False
		for line in lines:
			if line:
				inverse = float( line.split()[0] )
				
				if inverse:
					if options.maxres:
						if 1.0/inverse > options.maxres:
							newlines.append(line)
						else:
							pass
					else:
						pass
				else:
					firstline=True
					newlines.append(line)
		
		if not firstline:
			newlines = ['0 1.0']+newlines
		else:
			pass
		x=[]
		k=0		
		for line in newlines:
			inverse = float( line.split()[0] )
			print "NEW inverse is", inverse
			#x.append(float(k))
			
			values.append( float( line.split()[-1] ) )
			inversefreqs.append(inverse)

			element = ''
			if inverse: 
				element = '1/' + str(int( 1.0/inverse  ))
			else:
				element='0'
			#print "Therefore, turned into a number it is", element
			inversefreqslabels.append(element)
			x.append(k)
			k+=1
			
	
		print "values are", values
		
		xticks = []
		nele=len(values)
		print "\n\nnele is", nele
		import math
		factorOfTicks = int(math.ceil(nele/10.0) + 1.0)
	
		
		kk=0
		print "factorOfTicksIs", factorOfTicks
		print "And the final number of values is", len(values)
		for i in range(len(values)):
			if i == 0:	
				print "I always append the zero tick", inversefreqslabels[i]
				xticks.append(inversefreqslabels[i])
			
			if not (kk+1) % factorOfTicks:
				print "I have appended this tick!", inversefreqslabels[i]
				print "Because k+1 is", kk+1
				print "And k mod factorOfTicks is", (kk+1) % factorOfTicks
				xticks.append(inversefreqslabels[i])
			else:
				print "skipped this thick", inversefreqslabels[i]
				if i != len(values)-1 and i != 0:
					xticks.append('')
			
			#if i == len(values) -1:
			#	print "I always append the last tick", inversefreqslabels[i]
			#	xticks.append(inversefreqslabels[i])
				
			kk += 1
	
		#xticks[0]='0'
		plot_name = fscoutputname.replace('.txt','_PLOT.png')
		#if options.plotonly:
		if options.singleplot:
			if options.output:
				plot_name = options.path + '/' + os.path.basename(options.output)
				plot_name = plot_name.replace('.txt','.png')
			else:
				plot_name = options.path  + '/fsc_curves.png'
		else:
			plot_name = options.path + '/' + os.path.basename(fscoutputname).replace('.txt','_' + str(kont).zfill(len(fscs)) + '.png')
			#plot_name = plot_name.replace('.txt','_' + str(kont).zfill(len(fscs)) + '.png')
				
		fullinfo = {}

		'''
		Smooth out dips in plot
		'''
		
		if options.smooth:
			ret=maxima(inversefreqs,values,options.smooththresh)
			x=ret[0]
			values=ret[1]
		
		'''
		Fit polynomial to plot
		'''
		if options.polydegree:
			print "\n\n\nI will do a FIT!\n\n\n"
			'''
			Sigmoidal fit
			'''
			#ret=sigmoidfit(x,values)
			#x=ret[0]
			#values=ret[1]
			
			'''
			Poly fit
			'''
			difs0p5 = []
			difs0p143 = []

			#x=list(set(x))
	
			polycoeffs = numpy.polyfit(x, values, options.polydegree)
			yfit = numpy.polyval(polycoeffs, x)
			
			
			print "After the fit, yfit are", yfit

			for i in range(len(yfit)):
				dif0p5 = abs(yfit[i] - 0.5)
				difs0p5.append(dif0p5)
				fullinfo.update({dif0p5:i})
			
				dif0p143 = abs(yfit[i] - 0.143)
				difs0p143.append(dif0p143)
				fullinfo.update({dif0p143:i})
		
			difs0p5min1=min(difs0p5)
			difs0p5.remove(difs0p5min1)
			difs0p5min2=min(difs0p5)

			difs0p143min1=min(difs0p143)
			difs0p143.remove(difs0p143min1)
			difs0p143min2=min(difs0p143)

			fsc0p5minpixel1 = fullinfo[difs0p5min1]
			fsc0p5minpixel2 = fullinfo[difs0p5min2]
			fsc0p5pixel = (fsc0p5minpixel1 + fsc0p5minpixel2)/2
			fsc0p5freq1 = inversefreqs[fsc0p5minpixel1]
			fsc0p5freq2 = inversefreqs[fsc0p5minpixel2]

			fsc0p143minpixel1 = fullinfo[difs0p143min1]
			fsc0p143minpixel2 = fullinfo[difs0p143min2]
			fsc0p143pixel = (fsc0p5minpixel1 + fsc0p143minpixel2)/2
			fsc0p143freq1 = inversefreqs[fsc0p143minpixel1]
			fsc0p143freq2 = inversefreqs[fsc0p143minpixel2]

			fsc0p5freqavg = (fsc0p5freq1 + fsc0p5freq2)/2.0

			fsc0p143freqavg = (fsc0p143freq1 + fsc0p143freq2)/2.0

			fsc0p5resolution1 = ''
			fsc0p5resolution1label=''

			fsc0p143resolution1 = ''
			fsc0p143resolution1label=''

			if fsc0p5pixel and apix and options.boxsize:		
				fsc0p5resolution1 = (float(apix) * float(options.boxsize)) / float(fsc0p5pixel)
				fsc0p5resolution1label = "%.1f" % ( fsc0p5resolution1 )
			else:
				print "Method 1 for resolution calculation failed (there was a division by zero somewhere, or you forgot to provide --boxsize or --apix)"


			if fsc0p143pixel and apix and options.boxsize:		
				fsc0p143resolution1 = (float(apix) * float(options.boxsize)) / float(fsc0p143pixel)
				fsc0p143resolution1label = "%.1f" % ( fsc0p143resolution1 )

			elif not fsc0p5pixel:
				print "Method 1 for resolution calculation failed (there was a division by zero somewhere)"

			fsc0p5resolution2=''
			fsc0p5resolution2label=''

			if fsc0p5freqavg:	
				fsc0p5resolution2 = 1/fsc0p5freqavg
				fsc0p5resolution2label = "%.1f" % ( fsc0p5resolution2 )
			else:
				print "Method 2 for resolution calculation failed (there was a division by zero somewhere)"


			fsc0p143resolution2=''
			fsc0p143resolution2label=''

			if fsc0p143freqavg:	
				fsc0p143resolution2 = 1/fsc0p143freqavg
				fsc0p143resolution2label = "%.1f" % ( fsc0p143resolution2 )

			elif not fsc0p5resolution2:
				print "Method 2 for resolution calculation failed (there was a division by zero somewhere)"

			if sum(values)/len(values) == 1.0:
				print "The particles you are aligning are exactly the same; cannot compute a reasonable FSC curve for a particle with itself! (but I'll plot the straight line anyway)."
			else:
				print "FSC0.5 resolution calculations 1 and 2 are", fsc0p5resolution1, fsc0p5resolution2
				print "FSC0.143 resolution calculations 1 and 2 are", fsc0p143resolution1, fsc0p143resolution2
		

			final0p5='NA'
			if fsc0p5resolution2label:
				final0p5 = fsc0p5resolution2label
			elif fsc0p5resolution1label:
				final0p5 = fsc0p5resolution1label

			final0p143='NA'
			if fsc0p143resolution2label:
				final0p143 = fsc0p143resolution2label
			elif fsc0p143resolution1label:
				final0p143 = fsc0p143resolution1label	

		
			ww=open(plot_name.replace('.png','_RESvalues.txt'),'w')
			resline1 = 'fsc0.5=' +fsc0p5resolution1label + ' , fsc0.143=' +fsc0p143resolution1label + '\n'
			resline2 = 'fsc0.5=' +fsc0p5resolution2label + ' , fsc0.143=' +fsc0p143resolution2label+ '\n'
			resline3 = 'fsc0.5=' +str(final0p5) + ' , fsc0.143=' + str(final0p143) + '\n'
			reslines = [resline1, resline2, resline3]
			ww.writelines(reslines)
			ww.close()
		
		
			values=yfit
			#x=list(set(x))
			
			print "THerefore, valuues are", values
			print "x len is", len(x)
			print "yfit len is", len(yfit)
			
			print "x axis is", x
			print "yfit is", yfit
				
		'''
		Actual PLOT
		'''
		pylab.rc("axes", linewidth=2.0)
		
		colortouse = RGB_tuples[kont]
		mark = ''
		if options.nocolor:
			colortouse = 'k'
			try:
				mark = markers[kont]
			except:
				mark = markers[-1]
			
			#mark = itertools.cycle((',', '+', '.', 'o', '*')) 
    
    	#plt.plot(x,n, marker = marker.next(), linestyle='')
			
		
		pylab.plot(x, values, color=colortouse, linewidth=2, alpha=1.0, marker = mark, markersize = 10 )
	
		ax = fig.add_subplot(111)
		
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(16)
			tick.label1.set_fontweight('bold')
		for tick in ax.yaxis.get_major_ticks():
  			tick.label1.set_fontsize(16)
  			tick.label1.set_fontweight('bold')
  		 
  		pylab.xlabel('X Axis', fontsize=16, fontweight='bold')
  		pylab.ylabel('Y Axis', fontsize=16, fontweight='bold')

		plt.xticks(x,inversefreqslabels)
		print "Len of x and inversefreqslabels are", len(x), len(inversefreqslabels)

		print "len of ticks is", len(xticks)
		print "ticks are", xticks
		ax.xaxis.set_major_locator(MaxNLocator(nbins=len(xticks)))
		pylab.setp(ax, xticklabels=xticks)
		
		
		'''
		PLOT Threshold criteria as horizontal lines
		'''
		
		#print "cutoff is", type(options.cutoff), options.cutoff
		if options.cutoff and options.cutoff != 'None' and options.cutoff != 'none':
			print "\n\n\n\nTTTTTTTTTTTTTTTTT\nPlotting cutoff threshold"
			for thresh in options.cutoff:
				print "\nCurrent cutoff thresh is", thresh
				if float(thresh) == 0.5:
					yy1=[0.500]*len(values)	
					pylab.plot(x, yy1, 'k--', linewidth=3)
				
				if float(thresh) == 0.143:
					yy2=[0.143]*len(values)
					pylab.plot(x, yy2, 'k--', linewidth=3)
				
				if float(thresh) == 0.33 or float(thresh) == 0.3 or float(thresh) == 0.333:
					yy3=[0.333]*len(values)
					pylab.plot(x,yy3, 'k--', linewidth=3)
						
		#fit = pylab.plot(x, yfit, 'r-')

		

		
		#pylab.annotate("FSC0.5 = "+str(final0p5)+" A", xy=(0, 1), xytext=(300, -30), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
		#pylab.annotate("FSC0.143 = "+str(final0p143)+" A", xy=(0, 1), xytext=(300, -55), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
		#pylab.annotate("Sampling ="+"%.2f"%(apix) + " A/pixel", xy=(0, 1), xytext=(300, -80), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
		
		
		"""
		SIGMOID FITTING
		"""
		
		kont+=1
		
	#pylab.title(plot_name)
	pylab.ylabel('FSC')
	pylab.xlabel('Frequency 1/Angstroms')
	#pylab.grid(True)
	pylab.ylim([-0.1,1.2])
	plt.savefig(plot_name)
	

	return()
	
if __name__ == "__main__":
	main()
