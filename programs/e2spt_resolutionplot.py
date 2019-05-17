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
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2_utils import *
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
	usage = """Calculates, plots and optionally averages the FSC between multiple images and a reference."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--input", type=str, default=None, help="""Volume or stack of volumes to be compared to --ref""")

	parser.add_argument("--path", type=str, default='e2sptfscs', help="Results directory. If not specified, defaults to e2sptfscs/")

	parser.add_argument("--ref", type=str, default=None, help="Volume that will be 'static' (the 'reference' to which volumes in --input will be aligned to). The format MUST be '.hdf' or '.mrc' ")

	parser.add_argument("--nocolor", action='store_true', default=False, help="""Turns the ouput png(s) into grey scale figures. Instead of using different colors to distinguish between various curves on the same plot, this option will have the program automatically use different markers in black and white for each curve.""")

	parser.add_argument("--sym", type=str, default='c1', help = """Will symmetrize --ref and --input prior to FSC computation.""")

	parser.add_argument("--mask", type=str, default=None, help="""Mask processor applied to particles before fsc computation. Default is mask.sharp:outer_radius=-2""")

	parser.add_argument("--mirror",action="store_true", help="""If set, it will generate a mirrored version of --ref and recompute FSCs against it. This will be IN ADDITION to FSC computation of --input against the original, unmirrored --ref.""",default=False)

	#parser.add_argument("--preproc",type=str,default='',help="Any processor (as in e2proc3d.py) to be applied to each volumes prior to FSC computation. Typically this would be an automask.")

	#parser.add_argument("--savepreproc",action='store_true',default=False,help="""Default=False. Otherwise, save preprocessed/masked volumes for inspection.""")

	parser.add_argument("--averagefscs",action='store_true',default=False,help="""Default=False. Averages FSC curves if --input contains multiple images.""")

	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:2")

	parser.add_argument("--apix", type=float, help="""Provide --apix for automatic FSC calculation if you supply --plotonly and no volumes for --input and --ref, or if the apix of these is wrong.""", default=1.0)
	#parser.add_argument("--boxsize", type=float, help="(Probably not needed for anything)", default=0)

	parser.add_argument("--maxres", type=float, help="How far in resolution to extend the FSC curve on the x axis; for example, to see up to 20anstroms, provide --maxres=1.0. Default=15", default=1.0)
	parser.add_argument("--cutoff", type=str, help="Comma separated values of cutoff thresholds to plot as horizontal lines. Default=0.5, to turn of supply 'None'. ", default='0.5')

	parser.add_argument("--smooth",action="store_true", help="""Smooth out FSC curves by taking the average of a low value with a subsequent maxima.""", default=False)

	parser.add_argument("--smooththresh",type=float, help="""If --smooth is provided the curve will be smoothed only up to this resolution. Default is 100.""", default=100)

	#parser.add_argument("--fit",action="store_true", help="Smooth out FSC curves.", default=False)

	parser.add_argument("--polydegree",type=int, help="Degree of the polynomial to fit.", default=None)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID.",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	parser.add_argument("--plotonly",type=str, help="""FSC curves to plot in separate plots. Skips fsc curve generation. Provide .txt. files separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...""", default=None)

	parser.add_argument("--singleplot",action="store_true",help="It --plotonly provided, all FSC curves will be on the same plot/figure", default=False)

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)

	if options.mask:
		options.mask=parsemodopt(options.mask)

	#if options.preproc:
	#	options.preproc=parsemodopt(options.preproc)


	if options.cutoff and options.cutoff != 'None' and options.cutoff != 'none':
		options.cutoff = options.cutoff.split(',')
		print("Options.cutoff is", options.cutoff)
	else:
		options.cutoff = None

	from EMAN2_utils import makepath
	options = makepath(options,'sptres')

	print('\n\nafter making path, options.path is', options.path)

	if options.input:
		hdr = EMData(options.input,0,True)
		apix = hdr['apix_x']

	if options.apix:
		apix = float (options.apix)

	if options.plotonly:
		curves = options.plotonly
		curves = curves.split(',')

		if not options.singleplot:
			for curve in curves:
				print("Found this curve to plot", curve)
				fscplotter([curve],options,apix)

		elif options.singleplot:
			fscplotter(curves,options,apix)

		print("Done plotting")
		sys.exit()

	else:
		getfscs(options,apix)

		print("Done calculating FSCs and plotting them.")

	E2end(logger)
	return


def getfscs(options,apix):

	sym = None
	if options.sym and options.sym is not 'c1' and options.sym is not 'C1':
		sym = options.sym

	options.input
	print("\n inside getfscs options.input is", options.input)
	print('and the current directory is', os.getcwd())

	fyle = options.input.split('/')[-1]

	print("\n\n\n\n\nThe fyle to get image count from is", fyle)
	print("Whereas current dir is", os.getcwd())

	n = EMUtil.get_image_count( fyle )

	path = options.path
	fscs = []
	fscssym = []
	fscsm = []
	fscsmsym = []

	format=fyle[-4:]

	output = fyle.split(format)[0] + '_FSC.txt'

	'''
	Prepare the reference
	'''
	ref = EMData(options.ref,0)
	refsym=ref.copy()
	refsymm=ref.copy()

	if sym:
		refsym = symmetrize(ref,options)

	if options.mask:
		'''
		Make the mask first
		'''
		mask=EMData( int(ptcl["nx"]), int(ptcl["ny"]), int(ptcl["nz"]) )
		mask.to_one()
		mask.process_inplace(options.mask[0],options.mask[1])

		ref.mult(mask)
		if sym:
			refsym.mult(mask)

	#if options.preproc:
	#	ref.process_inplace(options.preproc[0],options.preproc[1])
	#	if sym:
	#		refsym.process_inplace(options.preproc[0],options.preproc[1])

	#if options.savepreproc:
	#	ref.write_image(options.ref.replace(options.ref[-4:],'_preproc'+options.ref[-4:]))
	#	if sym:
	#		refsym.write_image(options.ref.replace(options.ref[-4:],'_sym_preproc'+options.ref[-4:]))

	if options.mirror:
		t = Transform({'type':'eman','mirror':True})

		refm = ref.copy()
		refm.transform(t)

		refsymm = refsym.copy()
		refsymm.transform(t)

	'''
	Interate trough the particles
	'''
	for i in range(n):
		'''
		Prepare particle
		'''
		ptcl = EMData( fyle ,i)
		fscfilename = options.path + '/' + output.replace('.txt','_' + str(i).zfill(len(str(n))) + '.txt')

		if options.mask:
			ptcl.mult(mask)

		if sym:
			#ptclsym = symmetrize(ptcl,options)
			ptclsym = ptcl.process('xform.applysym',{'sym':sym})

			#fscfilename = options.output.replace('.txt','_' + str(i).zfill(len(str(n))) + '_' + options.sym + '.txt')
			fscfilenamesym = fscfilename.replace('.txt', '_' + sym + '.txt')

			#if options.mask:
			#	ptclsym.mult(mask)


		'''
		Calculate FSCs
		'''
		calcfsc(ptcl,ref,fscfilename,options)

		if sym:
			calcfsc(ptclsym,refsym,fscfilenamesym,options)

		'''
		Plot computed FSCs
		'''
		if not options.singleplot:
			fscplotter([fscfilename],options,apix)
			if sym:
				fscplotter([fscfilenamesym],options,apix)

		else:
			fscs.append(fscfilename)
			if sym:
				fscssym.append(fscfilename)


		if options.mirror:
			fscfilenamemirror = fscfilename.replace('.txt','_mirror.txt')
			print("\nThe mirror fsc file is", fscfilenamemirror)

			calcfsc(ptcl,refm,fscfilenamemirror,options)

			if not options.singleplot:
				fscplotter([fscfilenamemirror],options,apix)
			else:
				fscsm.append(fscfilenamemirror)

			if sym:
				fscsfilenamemirrorsym = fscfilenamemirror.replace('.txt','_'+sym+'.txt')
				print("\nThe mirror fsc file is", fscfilename)

				calcfsc(ptcl,refsymm,fscsfilenamemirrorsym,options)

				if not options.singleplot:
					fscplotter([fscsfilenamemirrorsym],options,apix)
				else:
					fscsmsym.append(fscsfilenamemirrorsym)

	if options.singleplot and fscs:
		fscplotter(fscs,options,apix,'apo',True)

	if sym and fscssym:
		fscplotter(fscssym,options,apix,sym,True)

	if options.mirror and fscsm:
		fscplotter(fscsm,options,apix,'mirror',True)
		if sym and fscsmsym:
			fscplotter(fscsmsym,options,apix,'mirror_'+sym,True)

	#print "fscs are",fscs
	if options.averagefscs:

		if fscs:
			fscaverager(options,fscs,'fscs_avg.txt')
			fscplotter([options.path+'/fscs_avg.txt'],options,apix,'apo_avg',True)

			if sym and fscssym:
				fscavgsym=fscaverager(options,fscssym,'fscs_avg_' + sym +'.txt')
				fscplotter([options.path+'/fscs_avg_'+sym+'.txt'],options,apix,'avg_' + sym,True)

		if options.mirror and fscsm:
			fscavgmirror=fscaverager(options,fscsm,'fscs_avg_mirror.txt')
			fscplotter([options.path+'/fscs_avg_mirror.txt'],options,apix,'avg_mirror',True)

			if sym and fscsmsym:
				fscavgmirrorsym=fscaverager(options,fscsmsym,'fscs_avg_mirror_'+sym+'.txt')
				fscplotter([options.path+'/fscs_avg_mirror_'+sym+'.txt'],options,apix,'avg_mirror_'+sym,True)
	return


def fscaverager(options,curves,outname):
	freqs=[]
	arrays=[]
	i=0
	for c in curves:
		print("opening fsc file",c)
		f=open(c,'r')
		lines=f.readlines()
		f.close()
		vals=[]

		print("has these many lines",len(lines))
		for line in lines:
			val=float(line.split()[-1].replace('\n',''))
			freq=float(line.split()[0])
			vals.append(val)
			if i==0:
				freqs.append(freq)
		print("these many vals and freqs",len(vals),len(freqs))
		valsarray=numpy.array(vals)
		arrays.append(valsarray)
		i+=1

	finalsum=[0]*len(arrays[0])
	for a in arrays:
		finalsum+=a
	finalavg=old_div(finalsum,len(arrays))

	#outavgtxt=open('fscs_avg_' + tag + '.txt','w')
	outavgtxt=open( options.path + '/' + outname,'w')
	outlines=['0.0\t1.0\n']

	for i in range(len(finalavg)):
		outline="%.10f\t%10f\n"%( freqs[i],finalavg[i])
		outlines.append(outline)

	outavgtxt.writelines(outlines)
	outavgtxt.close()

	return


def calcfsc(v1,v2,fscfilename,options):

	fsc = v1.calc_fourier_shell_correlation(v2)
	third = old_div(len(fsc),3)
	xaxis = fsc[0:third]
	fsc = fsc[third:2*third]
	apix=v1['apix_x']
	if options.apix:
		apix=options.apix
	saxis = [old_div(x,apix) for x in xaxis]
	Util.save_data(saxis[1],saxis[1]-saxis[0],fsc[1:], fscfilename)

	f=open(fscfilename,'r')
	lines=['0.0\t1.0\n']+f.readlines()
	f.close()

	g=open(fscfilename,'w')
	g.writelines(lines)
	g.close()

	return fsc


def symmetrize(vol,options):
	sym = options.sym
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	volsym=vol.copy()
	for i in range(1,nsym):
		dc=vol.copy()
		t=xf.get_sym(sym,i)
		dc.transform(t)
		volsym.add(dc)
	volsym.mult(old_div(1.0,nsym))
	return(volsym)


"""
SIGMOID FITTING
"""
def sigmoidfit(x,values):
	import numpy as np
	import scipy

	def sigmoid(p,x):
		x0,y0,c,k=p
		y = old_div(c, (1 + np.exp(-k*(x-x0)))) + y0
		return y

	def residuals(p,x,y):
		return y - sigmoid(p,x)

	def resize(arr,lower=0.0,upper=1.0):
		arr=arr.copy()
		if lower>upper:
			lower,upper=upper,lower
		arr -= arr.min()
		arr *= old_div((upper-lower),arr.max())
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

	print("\n\n\nPRESUMABLY THE SMOOTH CURVE IS pxsigp, with these many elements", len(pxsigp))
	print("Which are", pxsigp)
	print("\n\n\n")

	#plt.plot(xsig, ysig, '.', xsigp, pxsigp, '-')
	#plt.savefig(options.output)

	return(xsig,ysig)


def maxima(xaxis,yaxis,smooththresh):
	print("\n\nI have entered maxima!!!!!!!!!!!!!!!!!!!!!\n\n")
	print("At this point xaxis and yaxis len are", len(xaxis), len(yaxis))
	#xaxis=list(set(xaxis))
	#yaxis=list(set(yaxis))

	#xes=[xaxis[0]]
	#ymaxes=[yaxis[0]]

	#xaxis=list(set(xaxis))
	xes=[]
	ymaxes=[]


	for i in range(0,len(yaxis) -1):
		val=yaxis[i]
		#print 'currrent value is', val
		#print 'and next is', yaxis[i+1]
		#print "options.smooththresh is", smooththresh
		if val < max(yaxis[i+1:]) and old_div(1.0,xaxis[i+1]) > smooththresh:
			val = old_div(( val+ max(yaxis[i+1:]) ), 2.0)
			print('\nNew max smoothing value is', val)

		if val > min(yaxis[i+1:]) and old_div(1.0,xaxis[i+1]) < smooththresh:
			val = old_div(val, 2.0)
			print('\nNew min smoothing value is', val)

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

	print("\nExiting maxima, xes len is", len(xes))
	print("Exiting maxima, ymaxes len is", len(ymaxes))
	print("\nExiting maxima, xes are", xes)
	print("Exiting maxima, ymaxes are", ymaxes)

	#x=list(set(x))
	return(xes,ymaxes)



def fscplotter(fscs,options,apix=0.0,tag='',clearplot=False):

	fig = figure()

	if clearplot:
		plt.clf()

	#from itertools import product
	markers=["*","o","x","z","M"]
	#colors=["k","r","b","g","m","y","c"]
	#colorsANDmarkers = [a+b for a,b in product(colors,markers)]

	import colorsys
	import itertools


	N = len(fscs)
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	RGB_tuples = [colorsys.hsv_to_rgb(*x) for x in HSV_tuples]

	kont=0
	plot_name = ''
	for fscoutputname in fscs:

		print("\nI found THIS FSC file and will thus plot it", fscoutputname)

		f= open(fscoutputname,'r')
		lines = f.readlines()

		#if not options.boxsize:
		boxsize = len(lines)


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
						if old_div(1.0,inverse) > options.maxres:
							newlines.append(line)
						else:
							pass
					else:
						pass
				else:
					firstline=True
					newlines.append(line)

		if not firstline:
			newlines = ['0.0 1.0']+newlines
			print("\n!!!! added 0,1 as the first point to plot")
		else:
			pass
		x=[]
		k=0
		for line in newlines:
			inverse = float( line.split()[0] )
			print("NEW inverse is", inverse)
			#x.append(float(k))

			values.append( float( line.split()[-1] ) )
			inversefreqs.append(inverse)

			element = ''
			if inverse:
				element = '1/' + str(int( old_div(1.0,inverse)  ))
			else:
				element='0'
			#print "Therefore, turned into a number it is", element
			inversefreqslabels.append(element)
			x.append(k)
			k+=1


		print("values are", values)

		xticks = []
		nele=len(values)
		print("\n\nnele is", nele)
		import math
		factorOfTicks = int(math.ceil(old_div(nele,10.0)) + 1.0)


		kk=0
		print("factorOfTicksIs", factorOfTicks)
		print("And the final number of values is", len(values))

		xpositionswithlabel = []
		xpositionswithoutlabel = []

		for i in range(len(values)):

			if i == 0:
				print("I always append the zero tick", inversefreqslabels[i])
				#xticks.append(inversefreqslabels[i])
				xpositionswithlabel.append(i)

			if not (kk+1) % factorOfTicks:
				print("I have appended this tick!", inversefreqslabels[i])
				print("Because k+1 is", kk+1)
				print("And k mod factorOfTicks is", (kk+1) % factorOfTicks)
				#xticks.append(inversefreqslabels[i])
				xpositionswithlabel.append(i)
			else:
				print("skipped this thick", inversefreqslabels[i])
				if i != len(values)-1 and i != 0:
					#xticks.append('')
					xpositionswithoutlabel.append(i)


			xticks.append(inversefreqslabels[i])

			#if i == len(values) -1:
			#	print "I always append the last tick", inversefreqslabels[i]
			#	xticks.append(inversefreqslabels[i])

			kk += 1

		#xticks[0]='0'
		plot_name = fscoutputname.replace('.txt','_PLOT.png')
		#if options.plotonly:
		if options.singleplot:
			#if options.output:
			#	plot_name = options.path + '/' + os.path.basename(options.output)
			#	plot_name = plot_name.replace('.txt','.png')
			#else:
			plot_name = options.path  + '/fscs.png'
			if tag:
				plot_name = options.path  + '/fsc_' + tag + '.png'
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
			print("\n\n\nI will do a FIT!\n\n\n")
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


			print("After the fit, yfit are", yfit)

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
			fsc0p5pixel = old_div((fsc0p5minpixel1 + fsc0p5minpixel2),2)
			fsc0p5freq1 = inversefreqs[fsc0p5minpixel1]
			fsc0p5freq2 = inversefreqs[fsc0p5minpixel2]

			fsc0p143minpixel1 = fullinfo[difs0p143min1]
			fsc0p143minpixel2 = fullinfo[difs0p143min2]
			fsc0p143pixel = old_div((fsc0p5minpixel1 + fsc0p143minpixel2),2)
			fsc0p143freq1 = inversefreqs[fsc0p143minpixel1]
			fsc0p143freq2 = inversefreqs[fsc0p143minpixel2]

			fsc0p5freqavg = old_div((fsc0p5freq1 + fsc0p5freq2),2.0)

			fsc0p143freqavg = old_div((fsc0p143freq1 + fsc0p143freq2),2.0)

			fsc0p5resolution1 = ''
			fsc0p5resolution1label=''

			fsc0p143resolution1 = ''
			fsc0p143resolution1label=''

			if fsc0p5pixel and apix and boxsize:
				fsc0p5resolution1 = old_div((float(apix) * float(boxsize)), float(fsc0p5pixel))
				fsc0p5resolution1label = "%.1f" % ( fsc0p5resolution1 )
			else:
				print("Method 1 for resolution calculation failed (there was a division by zero somewhere, or you forgot to provide --boxsize or --apix)")


			if fsc0p143pixel and apix and boxsize:
				fsc0p143resolution1 = old_div((float(apix) * float(boxsize)), float(fsc0p143pixel))
				fsc0p143resolution1label = "%.1f" % ( fsc0p143resolution1 )

			elif not fsc0p5pixel:
				print("Method 1 for resolution calculation failed (there was a division by zero somewhere)")

			fsc0p5resolution2=''
			fsc0p5resolution2label=''

			if fsc0p5freqavg:
				fsc0p5resolution2 = old_div(1,fsc0p5freqavg)
				fsc0p5resolution2label = "%.1f" % ( fsc0p5resolution2 )
			else:
				print("Method 2 for resolution calculation failed (there was a division by zero somewhere)")


			fsc0p143resolution2=''
			fsc0p143resolution2label=''

			if fsc0p143freqavg:
				fsc0p143resolution2 = old_div(1,fsc0p143freqavg)
				fsc0p143resolution2label = "%.1f" % ( fsc0p143resolution2 )

			elif not fsc0p5resolution2:
				print("Method 2 for resolution calculation failed (there was a division by zero somewhere)")

			if old_div(sum(values),len(values)) == 1.0:
				print("The particles you are aligning are exactly the same; cannot compute a reasonable FSC curve for a particle with itself! (but I'll plot the straight line anyway).")
			else:
				print("FSC0.5 resolution calculations 1 and 2 are", fsc0p5resolution1, fsc0p5resolution2)
				print("FSC0.143 resolution calculations 1 and 2 are", fsc0p143resolution1, fsc0p143resolution2)


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

			print("THerefore, valuues are", values)
			print("x len is", len(x))
			print("yfit len is", len(yfit))

			print("x axis is", x)
			print("yfit is", yfit)

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
			tick.label1.set_fontsize(12)
			tick.label1.set_fontweight('bold')
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(12)
			tick.label1.set_fontweight('bold')


		pylab.xlabel('X Axis', fontsize=16, fontweight='bold')
		pylab.ylabel('Y Axis', fontsize=16, fontweight='bold')

		ax.set_xticks(xpositionswithlabel, minor=False)
		#ax.set_xticks(xpositionswithoutlabel, minor=True)

		#THIS WAS USED
		#plt.xticks(x,inversefreqslabels)

		labels = xticks

		print ("\n\n\n\n\n!!!!BEFORE len(labels)={}\n\n\nlabels={}".format(len(labels),labels))

		finallabels = []
		print ("\nxpositions with label are={}".format(xpositionswithlabel))
		for i in range(len(xticks)):
			if i in xpositionswithoutlabel:
				labels[i] = ''
				#xticks[i].label1.set_visible(False)
			elif i in xpositionswithlabel:
				labels[i] = inversefreqslabels[i]
				finallabels.append(inversefreqslabels[i])

		ax.set_xticklabels(finallabels)

		print("Len of x and inversefreqslabels are", len(x), len(inversefreqslabels))

		print("len of ticks is", len(xticks))
		print("ticks are", xticks)
		'''
		ax.xaxis.set_major_locator(MaxNLocator(nbins=len(xticks)))
		pylab.setp(ax, xticklabels=xticks)
		'''
		#to set INVISIBLE TICKS
		#make_invisible = True
		#if make_invisible:
		#	xticks = ax.xaxis.get_major_ticks()
		#	for i in range(len(xticks)):
		#		if i in xpositionswithoutlabel:
		#			xticks[i].label1.set_visible(False)




		print ("\n!!!! len(x)={}, len(values)={}, len(xticks)={}".format(len(x),len(values),len(xticks)) )

		'''
		PLOT Threshold criteria as horizontal lines
		'''

		#print "cutoff is", type(options.cutoff), options.cutoff
		if options.cutoff and options.cutoff != 'None' and options.cutoff != 'none':
			print("\n\n\n\nTTTTTTTTTTTTTTTTT\nPlotting cutoff threshold")
			for thresh in options.cutoff:
				print("\nCurrent cutoff thresh is", thresh)
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


	return

if __name__ == "__main__":
	main()
