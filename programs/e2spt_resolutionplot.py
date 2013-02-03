#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 02/01/2013
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
	
	parser.add_argument("--input", type=str, help="Volume or stack of volumes to be aligned to --ref; that is, volumes in --input will be 'moved' (rotated and translated) to find the best fit with --ref. The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--path", type=str, help="Results directory. If not specified, defaults to e2sptfscs/", default='e2sptfscs')

	parser.add_argument("--ref", type=str, help="Volume that will be 'static' (the 'reference' to which volumes in --input will be aligned to). The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--output", type=str, help="Name for the .txt file that will contain the FSC data. If not specified, a default name will be used, vol1_VS_vol2, where vol1 and vol2 are taken from --vol1 and --vol2 without the format.", default=None)

	parser.add_argument("--sym", type=str, default='c1', help = """Will symmetrize the reference (--vol1) and limit alignment of other structure (--vol2) against the model to searching the asymmetric unit only. 
								Then, after alignment, vol2 will be symmetrized as well, and the FSC will be calculated. 
								Note that this will only work IF the reference (--vol1) is ALREADY aligned to the symmetry axis as defined by EMAN2.""", default='c1')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=True)

	parser.add_argument("--mirror",action="store_true", help="If set, it will generate a mirrored version of --vol2, align that against --vol1 (the reference), and then compute FSCs. This will be done IN ADDITION to aligment and FSC of the unmirrored, original version of --vol2",default=False)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--npeakstorefine", type=int, help="The number of best 'coarse peaks' from 'coarse alignment' to refine in search for the best final alignment. Default=4.", default=4)
	
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average.", default=None)
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine_3d:search=2:delta=3:range=12", default="refine_3d:search=2:delta=3:range=12")
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--apix", type=float, help="Provide --apix for automatic FSC calculation (one of two methods) if you supply --curves and no volumes for --vol1 and --vol2", default=0)
	parser.add_argument("--boxsize", type=float, help="Provide --boxsize for automatic FSC calculation (one of two methods) if you supply --curves and no volumes for --vol1 and --vol2", default=0)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	#parser.add_argument("--plotonly",action="store_true", help="Assumes vol1 and vol2 are already aligned and fsc files will be provided through --curves; thus skips alignment and fsc curve generation", default=False)
	parser.add_argument("--fsconly",action="store_true", help="Assumes vol1 and vol2 are already aligned with respect to each other and thus skips alignment", default=False)
	parser.add_argument("--plotonly",type=str, help="FSC curves to plot in separate plots. Skips alignment and fsc curve generation. Provide .txt. files separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...", default=None)
	parser.add_argument("--singleplot",action="store_true",help="If --singleplot and --plotonly are provided, they will be plotted in the same image.", default=False)
		
	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)

	if not options.output and not options.plotonly:
		print "ERROR: Unless you provide .txt files through --plotonly, you must specify an --output in .txt format."
	elif options.output and not options.plotonly:
		if '.txt' not in options.output:
			print "ERROR: Output must be in .txt format"
			sys.exit()
		
	findir=os.listdir(os.getcwd())
	if options.path not in findir:
		os.system('mkdir ' + options.path)
	
	if options.plotonly:
		curves = options.plotonly
		curves = curves.split(',')
		
		if not options.singleplot:
			for curve in curves:
				fscplotter([curve],options)
				
		elif options.singleplot:
			fscplotter(curves,options)
		
		print "Done plotting"
		sys.exit()
		
	elif options.fsconly:
		getfscs(options)
	
		print "Done calculating FSCs and plotting them."
		sys.exit()
	
	elif options.align:
		if options.sym and options.sym is not 'c1' and options.sym is not 'C1':		
			ref = EMData(options.ref,0)
			ref = symmetrize(ref,options)
			options.ref = options.path + '/' + os.path.basename(options.ref).replace('.','_' + options.sym + '.')
			ref.write_image(options.ref,0)
		
		ptclali = alignment(options)
		inputbackup = options.input
		options.input = ptclali
		
		getfscs(options)

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
			getfscs(options)
		
	E2end(logger)
	return


def getfscs(options):
	n = EMUtil.get_image_count(options.input)
	
	path = options.path
	fscs = []
	fscsm = []
	for i in range(n):
		ptcl = EMData(options.input,i)
		ref = EMData(options.ref,0)
		fscfilename = path + '/' + options.output.replace('.txt','_' + str(i).zfill(len(str(n))) + '.txt')
		if options.sym and options.sym is not 'c1' and options.sym is not 'C1':
			ptcl = symmetrize(ptcl,options) 
			ref = symmetrize(ref,options)
			fscfilename = options.output.replace('.txt','_' + str(i).zfill(len(str(n))) + '_' + options.sym + '.txt')
			
		calcfsc(ptcl,ref,fscfilename)
		
		if not options.singleplot:
			fscplotter([fscfilename],options)
		else:
			fscs.append(fscfilename)
	
		if options.mirror:
			t = Transform({'type':'eman','mirror':True})
			refm = ref.copy()
			refm.transform(t)
			fscmfilename = fscfilename.replace('.txt','_VSmirror.txt')
			
			calcfsc(ptcl,refm,fscmfilename)
			
			if not options.singleplot:
				fscplotter([fscmfilename],options)
			else:
				fscsm.append(fscmfilename)
	
	if options.singleplot:
		fscplotter(fscs,options)
		fscplotter(fscsm,options)
	return


def alignment(options):
	aligner=options.align
	aligner=aligner.split(':')
	newaligner=''
		
	if 'sym='+options.sym not in aligner:
		for element in aligner:
			if 'sym' in element:
				element='sym=' +options.sym
			newaligner=newaligner + element + ":"
		if newaligner[-1] == ':':
			newaligner = newaligner[:-1]
			options.align=newaligner

	alivolfile = os.path.basename(options.input).split('.')[0] + '_VS_' + os.path.basename(options.ref).split('.')[0] + '.hdf'
	alicmd = 'e2spt_classaverage.py --path=' + options.path + ' --input=' + str(options.input) + ' --output=' + alivolfile + ' --ref=' + str(options.ref) + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose) + ' --mask=' + options.mask + ' --lowpass=' + options.lowpass + ' --align=' + options.align + ' --parallel=' + options.parallel + ' --ralign=' + options.ralign + ' --aligncmp=' + options.aligncmp + ' --raligncmp=' + options.raligncmp + ' --shrink=' + str(options.shrink) + ' --shrinkrefine=' + str(options.shrinkrefine) + ' --saveali' + ' --normproc=' + options.normproc + ' --sym=' + options.sym + ' --breaksym'
			
	os.system(alicmd)
	print "BEFORE ALIGNMENT, input was", options.input
	print "BEFORE ALIGNMENT, ref was", options.ref
	print "Returning this from ALIGNMENT", alivolfile
	return options.path + '/' + alivolfile


def calcfsc(v1,v2,fscfilename):
	fsc = v1.calc_fourier_shell_correlation(v2)
	third = len(fsc)/3
	xaxis = fsc[0:third]
	fsc = fsc[third:2*third]
	apix=v1['apix_x']
	saxis = [x/apix for x in xaxis]
	Util.save_data(saxis[1],saxis[1]-saxis[0],fsc[1:], fscfilename)
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
	volsym.mult(1.0/nsym)	
	return(volsym)
	
		
def fscplotter(fscs,options):
	fig = figure()

	#from itertools import product
	#markers=["-","--","x"]
	#colors=["k","r","b","g","m","y","c"]
	#colorsANDmarkers = [a+b for a,b in product(colors,markers)]
	
	import colorsys
	N = len(fscs)
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
	
	kont=0
	plot_name = ''
	for fscoutputname in fscs:
	
		print "I found THIS FSC file and will thus plot it", fscoutputname

		f= open(fscoutputname,'r')
		lines = f.readlines()

		if not options.boxsize:
			options.boxsize = len(lines)

		k=0
		x = []

		values = []
		inversefreqs = []
		inversefreqslabels = []

		xticks = []
		nlines=len(lines)
		arbfactor = round(0.2*nlines)
		factorOfTicks = round(nlines/(round(nlines / (0.2*nlines))))
		for line in lines:
			if line:
				x.append(float(k))
				values.append( float( line.split()[-1] ) )
				inverse =  float( line.split()[0] )
				#print "inverse is"
				inversefreqs.append(inverse)

				element = ''
				if inverse:
					element = '1/' + str(int(round( 1.0/inverse  )))
				else:
					element='0'
				#print "Therefore, turned into a number it is", element
				inversefreqslabels.append(element)
				if not k % factorOfTicks:

					print "I have appended this tick!", element
					print "Because k is", k
					print "And k mod factorOfTicks is", k % factorOfTicks
					xticks.append(element)
				k += 1

		plot_name = fscoutputname.replace('.txt','_PLOT.png')
		if options.plotonly and options.singleplot:
			if options.output:
				plot_name = path + '/' + options.output
				plot_name = plot_name.replace('.txt','.png')
			else:
				plot_name = path + '/fsc_curves.png'

		fullinfo = {}

		difs0p5 = []
		difs0p143 = []

		polycoeffs = numpy.polyfit(x, values, 30)
		yfit = numpy.polyval(polycoeffs, x)

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

		if fsc0p5pixel and options.apix and options.boxsize:		
			fsc0p5resolution1 = (float(options.apix) * float(options.boxsize)) / float(fsc0p5pixel)
			fsc0p5resolution1label = "%.1f" % ( fsc0p5resolution1 )
		else:
			print "Method 1 for resolution calculation failed (there was a division by zero somewhere, or you forgot to provide --boxsize or --apix)"


		if fsc0p143pixel and options.apix and options.boxsize:		
			fsc0p143resolution1 = (float(options.apix) * float(options.boxsize)) / float(fsc0p143pixel)
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
		
		
		'''
		ACTUAL PLOT
		'''
		pylab.rc("axes", linewidth=2.0)
		#print "\n\n\nThe selected color is", RGB_tuples[kont]
		#print "And there are these many in total", len(RGB_tuples)
		#print "And at the moment, kont is", kont
		#print "\n\n\n"
		pylab.plot(x, values, color=RGB_tuples[kont], linewidth=2)
	
		yy1=[0.5]*len(values)	
		yy2=[0.143]*len(values)
		#pylab.axhline(linewidth=0.5)
		#pylab.axhline(linewidth=0.143)
		
		#ax1 = fig.add_subplot(111)

		#ax1.axhline(linewidth=4, color="k")
		#ax1.axvline(linewidth=4, color="k")
		
		pylab.plot(x, yy1, 'k--', linewidth=1)
		pylab.plot(x, yy2, 'k--', linewidth=1)

		#fit = pylab.plot(x, yfit, 'r-')

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
		
		#pylab.annotate("FSC0.5 = "+str(final0p5)+" A", xy=(0, 1), xytext=(300, -30), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
		#pylab.annotate("FSC0.143 = "+str(final0p143)+" A", xy=(0, 1), xytext=(300, -55), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
		#pylab.annotate("Sampling ="+"%.2f"%(apix) + " A/pixel", xy=(0, 1), xytext=(300, -80), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
			
		"""
		SIGMOID FITTING
		"""
		
		#def sigmoid(p,x):
		#    x0,y0,c,k=p
		#    y = c / (1 + np.exp(-k*(x-x0))) + y0
		#    return y
		
		#def residuals(p,x,y):
		#    return y - sigmoid(p,x)
		
		#def resize(arr,lower=0.0,upper=1.0):
		#    arr=arr.copy()
		#    if lower>upper: lower,upper=upper,lower
		#    arr -= arr.min()
		#    arr *= (upper-lower)/arr.max()
		#    arr += lower
		#    return arr
		
		# raw data
		#xsig = np.array(x,dtype='float')
		#ysig = np.array(values,dtype='float')
		
		#xsig = resize(-xsig,lower=0.3)
		#ysig = resize(ysig,lower=0.3)
		#print(xsig)
		#print(ysig)
		#p_guess=(np.median(xsig),np.median(ysig),1.0,1.0)
		#p, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,p_guess,args=(xsig,ysig),full_output=1,warning=True)  
		
		#xsig0,ysig0,c,k=p
		#print('''\
		#x0 = {x0}
		#y0 = {y0}
		#c = {c}
		#k = {k}
		#format(xsig0=xsig0,ysig0=ysig0,c=c,k=k))

		#xsigp = np.linspace(0, 1.1, 1500)
		#pxsigp=sigmoid(p,xsigp)

		#print "\n\n\nPRESUMABLY THE SMOOTH CURVE IS pxsigp, with these many elements", len(pxsigp)
		#print "Which are", pxsigp
		#print "\n\n\n"

		# Plot the results
		#plt.plot(xsig, ysig, '.', xsigp, pxsigp, '-')
		
		#plt.savefig(options.output)
		
		"""
		SIGMOID FITTING
		"""
		
		kont+=1
		
	pylab.title(plot_name)
	pylab.ylabel('FSC')
	pylab.xlabel('Frequency 1/Angstroms')
	#pylab.grid(True)
	pylab.ylim([0,1.2])
	plt.savefig(plot_name)
	

	return()
	
if __name__ == "__main__":
	main()
