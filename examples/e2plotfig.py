#!/usr/bin/env python
'''
====================
Author: Jesus Galaz-Montoya - 2017, Last update: 12/Sep/2017
====================

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
'''
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
import matplotlib
matplotlib.use('Agg',warn=False)

import matplotlib.pyplot as plt
import pylab

#import matplotlib.colors as mcol
#import matplotlib.cm as cm
import numpy as np
import colorsys

import sys, os

from EMAN2 import *


def main():

	progname = os.path.basename(sys.argv[0])
	usage = """Plot single column or double column datasets in a high-resolution, ready to publish form. Exampek command:
	For a single file with x y values:
	e2plotfig.py --data mydata.txt <options>

	For multiple files with x y data:
	e2plotfig.py mydata*txt <options>

	or

	e2plotfig.py --data mydata1.txt,mydata2.txt,mydata3.txt,...mydataN.txt <options>

	To make one plot using separate source files for x and y data:
	e2plotfig.py --datax xvals.txt --datay yvals.txt <options>

	To make multiple plots using separate source files for x and y data:
	e2plotfig.py --datax xvals1.txt,xvals2.txt,...,xvalsN.txt --datay yvals1.txt,yvals2.txt,...,yvalsN.txt <options>


	"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--binwidth", type=float, default=0.0, help="""default=0.0 (not used). requires --histogram. Y axes. Enforce this value for the width of histogram bins (it will be used to calculate --nbins)""")

	parser.add_argument("--data", type=str, default='', help="""default=None (not used). Text file(s) with two column of values meant to be plotted on the X and Y axes. If supplying multiple files, separate them by commas.""")
	parser.add_argument("--datax", type=str, default='', help="""default=None (not used). Text file(s) with a single column of values meant to be plotted on the X axis. If not provided, the X axis will go from 0 to n, where n is the number of values in --datay. If supplying multiple files, separate them by commas (the number of files for --datax and --datay must be the same).""")
	parser.add_argument("--datay", type=str, default='', help="""default=None (not used). Text file(s) with a single column of values meant to be plotted on the  Y axis. If not provided, the Y axis will go from 0 to n, where n is the number of values in --datax. If supplying multiple files, separate them by commas (the number of files for --datax and --datay must be the same).""")

	parser.add_argument("--highresolution", action='store_true', default=False, help="""default=False. If on, this option will enforce writing high-resolution plots (dpi 300), ready for publication, as opposed to lower resolution plots which take up less space on your computer (dpi 150).""")
	parser.add_argument("--histogram", action='store_true', default=False, help="""default=False. If on, this will make the resulting plot a histogram. If --nbins not supplied, the parameter will be automatically calculated.""")

	parser.add_argument("--individualplots", action='store_true', default=False, help="""default=False. in addition to plotting the motion for all frames in all images in a single plot, generate individual plots per image""")
	
	parser.add_argument("--labelxaxis", type=str,default='x',help="""Default=x. Label for x axis (specify units through --unitsx.""")
	parser.add_argument("--labelyaxis", type=str,default='y',help="""Default=y. Label for y axis (specify units through --unitsy.""")
	parser.add_argument("--labeltitle", type=str,default='',help="""Default=None. Title for figure.""")
	parser.add_argument("--legend",type=str,default='',help=""""Default=None. If you are plotting only 1 curve, or --individualplots is on, and you desire a specific legend for the data series in each plot, supply it here as a string with no spaces. You can provide any string without spaces; if you need spaces, add an underscore instead and the program will replace the underscore with a space; for exampe 'ribosome_80s' will appear as 'ribosome 80s'.""")
	parser.add_argument("--linesoff", action='store_true', default=False, help="""Default=False. This requires --markerson and will get rid of the line uniting data points (the plot will be like a scatter plot).""")
	
	parser.add_argument("--markerson", action='store_true', default=False, help="""Default=False. This will enforce markers in the plot for each data point (like a scatter plot, but the points will still be united by a line). Will be switched on if --linesoff is provided.""")
	parser.add_argument("--marker",type=str,default='',help=""""Default=None. If you are plotting only 1 curve, or --individualplots is on, and you desire a specific marker, supply it here (for example o, or *, or x. Will default to 'x' if --linesoff is provided.""")
	
	parser.add_argument("--maxx",type=float,default=None,help="""Default=None. Maximum value to plot in X. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--maxy",type=float,default=None,help="""Default=None. Maximum value to plot in Y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--minx",type=float,default=None,help="""Default=None. Minimum value to plot in X. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--miny",type=float,default=None,help="""Default=None. Minimum value to plot in Y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")

	parser.add_argument("--mult",type=float,default=None,help="""Default=None. The data will be multiplied by this factor immediatebly prior to plotting. For example, if the data is in the order of magnitude of 10^6, you might say --mult=0.000001, and in --xunits change, e.g., nm^3 to nm^3x10^6""")

	parser.add_argument("--nbins", type=int,default=0,help="""Default=0 (not used). Requires --histogram. Number of bins for histogram. If not provided, the optimal bin number will be automatically calculated based on bin-width, computed using Scott's normal reference rule, width = (3.5*std)/cuberoot(n), where 'std' is the standard deviation of the mean intensity distribution of population and n is the number of mean intensity values considered (this is affected by --removesigma). Then, bins will be nbins = (max(intensities) - min(intensities)) / width.""")	
	parser.add_argument("--nocolor", action='store_true', default=False, help="""Default=False. Plots are colored, by default; don't be cheap; clear communication and representation pays off; or consider publishing in online open source journals that don't charge extra for color figures.""")
	parser.add_argument("--normalize", action='store_true', default=False, help="""Default=False. This option will normalize all plots to be scaled between 0 and 1.""") 
	
	parser.add_argument("--outputtag",type=str,default='plotfig',help="""Default=plotfig. String common to all automatically generated output files. For example, --outputtag=myplot will generate myplot1.png, myplot2.png, ..., myplotN.png""")

	parser.add_argument("--path", type=str,default='plotfig',help="""Default=de_plots. Name of the directory where to store the output results.""")
	parser.add_argument("--ppid", type=int, default=-1,help="Default=-1. Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--scaleaxes", action='store_true', default=False, help="""Default=False. This will force the axes to be on the same scale.""") 

	parser.add_argument("--unitsx", type=str,default='AU',help="""Default=AU (arbitrary units). Units for the x axis.'microns' or 'mu' and 'angstroms' or 'A' (and '1/angstroms' or '1/A') will be replaced by the appropriate symbol. You can provide any string without spaces; if you need spaces, add an underscore instead and the program will replace the underscore with a space; for exampe 'GPU_h' will appear as 'GPU h'.""")
	parser.add_argument("--unitsy", type=str,default='AU',help="""Default=AU (arbitrary units). Units for the y axis.'microns' or 'mu' and 'angstroms' or 'A' (and '1/angstroms' or '1/A')  will be replaced by the appropriate symbol. You can provide any string without spaces; if you need spaces, add an underscore instead and the program will replace the underscore with a space; for exampe 'GPU_h' will appear as 'GPU h'.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	datafiles=[]
	if options.data:
		if options.datax or options.datay:
			print("\n(e2plotfig)(main) ERROR: provide --data OR --datax AND/OR --datay.")
			sys.exit(1)
	elif not options.data:
		if not options.datax and not options.datay:
			if args:
				datafiles=args
			elif not args:
				print("\n(e2plotfig)(main) ERROR: provide at least one of --data, --datax, or --datay.")
				sys.exit(1)

	if options.linesoff:
		options.markerson=True

	if options.unitsx:
		if options.unitsx=='angstroms' or options.unitsx=='Angstroms' or options.unitsx=="A" or options.unitsx=="ANGSTROMS":
			options.unitsx = u"\u212B"
		if options.unitsx=='1/angstroms' or options.unitsx=='1/Angstroms' or options.unitsx=="1/A" or options.unitsx=="1/ANGSTROMS":
			options.unitsx = '1/'+u"\u212B"
		if options.unitsx=='microns' or options.unitsx=='mu' or options.unitsx=="Microns" or options.unitsx=="mu" or options.unitsx=="MICRONS" or options.unitsx=="MU" :
			options.unitsx = u"\u00B5"			

	if options.unitsy:
		if options.unitsy=='angstroms' or options.unitsy=='Angstroms' or options.unitsy=="A" or options.unitsy=="ANGSTROMS":
			options.unitsy = u"\u212B"
		if options.unitsy=='1/angstroms' or options.unitsy=='1/Angstroms' or options.unitsy=="1/A" or options.unitsy=="1/ANGSTROMS":
			options.unitsy = '1/'+u"\u212B"
		if options.unitsy=='microns' or options.unitsy=='mu' or options.unitsy=="Microns" or options.unitsy=="mu" or options.unitsy=="MICRONS" or options.unitsy=="MU" :
			options.unitsy = u"\u00B5"	

	xaxes={}
	yaxes={}
	datadict={}

	lines=[]
	if options.data:
		datafiles = options.data.split(',')
		k=0
		

		if len(datafiles) < 2:
			options.individualplots = True

		for f in datafiles:
			
			xaxis=[]
			yaxis=[]
			with open( f ) as datafile: 
				lines=datafile.readlines()
				if options.verbose:
					print("\nreading file {}".format(f))
					if options.verbose >9:
						print("\nlines are", lines)
	
				if lines:
					lines = fixlines(lines)
					xaxis = [ float(line.replace('\n','').split()[0]) for line in lines ]
					yaxis = [ float(line.replace('\n','').split()[1]) for line in lines ]
				else:
					print("\nERROR: source file {} seems empty; no lines read".format(f))	
					sys.exit(1)

				if options.normalize:
					yaxis = normalize(yaxis)

				#if xaxis and yaxis:
				xaxes.update({k:xaxis})
				yaxes.update({k:yaxis})
				datadict.update({k:[xaxis,yaxis]})
				
				k+=1


	elif not options.data:
		if options.datax:
			dataxfiles=options.datax.split(',')

			if len(dataxfiles) < 2:
				options.individualplots = True

			if options.datay:
				datayfiles=options.datay.split(',')
				if len(dataxfiles) != len(datayfiles):
					print("\n(e2plotfig)(main) ERROR: --datax and --datay must contain the same number of files. Now, nx files=%d, ny files=%d".format(len(dataxfiles),len(datayfiles)))
					sys.exit(1)

			k=0
			for fx in dataxfiles:
				with open( fx ) as dataxfile: 
					lines=dataxfile.readlines()
					xaxis = [ float(line.replace('\n','').split()[0]) for line in lines ]
					xaxes.update({k:xaxis})
					if not options.datay:
						yaxis = list(range(len(xaxis)))
						yaxes.update({k:yaxis})

						if options.normalize:
							xaxis = normalize(xaxis)

						datadict.update({k:[xaxis,yaxis]})
					k+=1


		lines=[]
		if options.datay:
			datayfiles=options.datay.split(',')
			k=0
			for fy in datayfiles:
				with open( fy ) as datayfile: 
					lines=datayfile.readlines()
					yaxis = [ float(line.replace('\n','').split()[1]) for line in lines ]
					yaxes.update({k:yaxis})

					if options.normalize:
						yaxis = normalize(yaxis)
	
					if not options.datax:
						xaxis = list(range(len(yaxis)))
						xaxes.update({k:xaxis})

						datadict.update({k:[xaxis,yaxis]})
					k+=1

	from EMAN2_utils import makepath
	options = makepath(options)

	fig,ax=resetplot()

	plotdata(options,datadict)

	E2end(logger)


	return


def fixlines(inlines):
	n=len(inlines)
	newlines=[]
	for i in range(0,n):
		inlines[i] = inlines[i].replace(", ",' ')	
		inlines[i] = inlines[i].replace(",",' ')
		inlines[i] = inlines[i].replace("x",'')
		inlines[i] = inlines[i].replace("y",'')
		inlines[i] = inlines[i].replace("z",'')
		inlines[i] = inlines[i].replace("=",'')
		inlines[i] = inlines[i].replace("_",' ')
		inlines[i] = inlines[i].replace("\n",'')
		inlines[i] = inlines[i].replace("\t",' ')
		inlines[i] = inlines[i].replace("  ",' ')

		if inlines[i]:
			newlines.append(inlines[i])
		else:
			print("\nartifactual line (number {}) removed".format(i))

	return newlines


def normalize(data):
	dminusmin = [d-min(data) for d in data]
	dminusminovermax = [old_div(d,max(dminusmin)) for d in dminusmin]

	return dminusminovermax


def resetplot(figsize=None):
	plt.clf()

	matplotlib.rc('xtick', labelsize=14) 
	matplotlib.rc('ytick', labelsize=14) 
	#matplotlib.rc('text', usetex=True) 
	font = {'weight':'bold','size':14}
	matplotlib.rc('font', **font)
	#matplotlib.rc('font', weight='bold') 

	fig = plt.figure()
	if figsize:
		fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(1,1,1)

	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both', reset=False, which='both', length=8, width=3)

	plt.rcParams.update({'figure.max_open_warning': 0})

	return fig,ax

'''
def resetcolorbar(cbminval,cbmaxval):

	# Make a user-defined colormap.
	cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","y","g","c","b","m"])

	# Make a normalizer that will map the values from
	# [start,end+1] -> [0,1].
	cnorm = mcol.Normalize(vmin=cbminval,vmax=cbmaxval)

	# Turn these into an object that can be used to map values to colors and
	# can be passed to plt.colorbar().
	colorpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
	colorpick.set_array([])

	return colorpick
'''


def plotdata( options, data ):

	resolution=150
	if options.highresolution:
		resolution=300

	fig,ax = resetplot()
	#cpick = resetcolorbar(options.lowestangle,options.highestangle)
	ndata = len(data)

	print("\nllllllen(data) {}".format(len(data)))
	#colorbar = False
	#altxaxis = None
	#colorstart = options.lowestangle
	#colorstep = options.tiltstep

	if not options.individualplots:
		N = len(data)
		HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
		RGB_tuples = [colorsys.hsv_to_rgb(*x) for x in HSV_tuples]

		#import string
		#markers=string.printable
		markers = list(matplotlib.markers.MarkerStyle.markers.keys())
		marker=''

		
		#markers=["*","o","x",">","<","|","z","m",]

		#mark = ''

		for k in data:
			#if count==ndata-1:
			#	colorbar=True
			transparent=False
			if len(data) > 1:
				transparent=True

			color = RGB_tuples[k]
			if options.nocolor or options.markerson and not marker:
				marker = markers[k]
				if options.nocolor:
					color = 'k'
				
				#try:
				#	mark = markers[k]
				#except:
				#	mark = markers[-1]
			
			if options.verbose > 9:
				print("\n(e2plotfig)(plotdata) plotting dataset n = {}".format(k))
				print("color is {}".format(color))
				print("marker is {}".format(marker))
			
			plotfig( options, fig, ax, data[k][0], data[k][1], k, color, marker, transparent )
			
			areay = round(sum(data[k][1]),2)	
			with open(options.path + '/' + options.outputtag + '_areas.txt','a') as areasfile:
				areasfile.write( str(k)+ ' ' + str(areay) + '\n')
			 

		filetosave = options.path + '/' + options.outputtag + '.png'
		fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')# transparent=True) #, bbox_extra_artists=(lgd,), bbox_inches='tight'
		plt.close('fig')
		print("\n(e2plotfig)(plotdata) saving figure {}".format(filetosave))

	elif options.individualplots:
		#colorstep=0
		if options.linesoff:
			options.marker='x'
		for k in data:
			#colorbar=False
			fig,ax = resetplot()
			print("\n(e2plotfig)(plotdata) plotting individual plot for dataset n = {}".format(k))
			plotfig( options, fig, ax, data[k][0], data[k][1], k )
			
			filetosave = options.path + '/' + options.outputtag +'.png'
			if len(data) > 1:	
				filetosave = options.path + '/' + options.outputtag + str( k ).zfill( len(str(ndata)))+'.png'
			
			fig.savefig( filetosave, dpi=resolution, bbox_inches='tight', format='png')
			plt.close('fig')

	return


def plotfig( options, fig, ax, datax, datay, count, colorthis='k', markerthis='', transparent=False )	:
	
	if options.verbose > 9:
		print('\ndatax={}, datay={}'.format(datax,datay))
	alphaval=1.0
	if transparent:
		alphaval=0.5
	n=len(datay)

	if options.mult:
		for i in range(len(datay)):
			datay[i] = datay[i]*options.mult
	#xaxis=range(n)
	#if altxaxis:
	#	xaxis=altxaxis

	#colorthis = 'k'
	#if colorstart and colorstep:
	#	colorthis = cpick.to_rgba( colorstart + kplot*colorstep)
	label=str(count)

	if options.individualplots:
		if options.marker:
			markerthis = options.marker
	
		if options.legend:
			label = options.legend.replace('_',' ')

	elif markerthis:
		label=str(markerthis)
		
	#ax.legend()
	linestyle='-'
	linewidth=2
	if options.linesoff:
		linestyle=''
		linewidth=0
	
	if not options.histogram:
		ax.plot( datax, datay, linestyle=linestyle, linewidth=linewidth, marker=markerthis, markersize=10, markeredgewidth=5, color=colorthis, label=label, alpha=alphaval)
	elif options.histogram:
		nbins = calcbins(options,datay)
		print("\ndatay is",datay)
		n, bins, patches = plt.hist(datay, nbins, label=label, histtype='bar', edgecolor='black', linewidth=2.0, alpha=alphaval)#, facecolor=colorthis,normed=1)

	print("\noptions.miny is {}".format(options.miny))
	if options.miny != None or options.maxy != None:
		miny=min(datay)
		maxy=max(datay)
		if options.miny != None:
			miny=options.miny
		if options.maxy != None:
			maxy=options.maxy
		ax.set_ylim( miny, maxy )

	if options.minx !=None or options.maxx != None:
		minx=min(datax)
		maxx=max(datax)
		if options.minx != None:
			minx=options.minx
		if options.maxx != None:
			maxx=options.maxx
		ax.set_xlim( minx, maxx )

	#plt.axis('scaled')
	#plt.axis('equal')
	if options.scaleaxes:
		plt.gca().set_aspect('equal', adjustable='box')

	handles, labels = ax.get_legend_handles_labels()
	
	#ax.legend(handles, labels, loc='center right', bbox_to_anchor=(1.3, 0.5))
	
	#plt.legend()
	plt.legend(frameon=False, bbox_to_anchor=(1.05,1), loc="upper left", borderaxespad=0)
	
	#if errors:
	#	lines = {'linestyle': 'None'}
	#	plt.rc('lines', **lines)
	#	plt.errorbar(xaxis,values,errors,markersize=8,linewidth=1,fmt='',marker='o',color='k',markerfacecolor=None,markeredgecolor='k',capsize=5, capthick=1)

	ax.set_title(options.labeltitle, fontsize=16, fontweight='bold')
	if options.unitsx:
		ax.set_xlabel(options.labelxaxis + ' (' + options.unitsx.replace('_',' ') + ')', fontsize=16, fontweight='bold')
	else:
		ax.set_xlabel(options.labelxaxis, fontsize=16, fontweight='bold')

	if options.unitsy:
		ax.set_ylabel(options.labelyaxis + ' (' + options.unitsy.replace('_',' ') + ')', fontsize=16, fontweight='bold')
	else:
		ax.set_ylabel(options.labelyaxis, fontsize=16, fontweight='bold')

	#if colorbar:
	#	cbr=plt.colorbar(cpick)
	#	cbr.set_label("Tilt angle (degrees)",fontsize=18, fontweight='bold')
	
	if not options.individualplots:
		plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
	
	return


def calcbins(options,data):
	std = np.std(data)
	mean = np.mean(data)
	
	statistics = ['mean='+str(mean) + ' std='+str(std) + '\n']
	print("The number of particles kept is", len(data))
	print("The standard deviation of the mean intensity distribution for this population is", std)
	
	if not std:
		print("ERROR: std={}, which means all data values are the same.".format(std))
		sys.exit()
		
	cuberoot = np.power(len(data),old_div(1.0,3.0))
	#print "The cuberoot of n is", cuberoot
	width = old_div((3.5*std),cuberoot)
	print("Therefore, according to Scott's normal reference rule, width = (3.5*std)/cuberoot(n), the width of the histogram bins will be", width)

	if options.binwidth:
		width = options.binwidth
	
	nbins=3
	if options.nbins:
		if options.binwidth:
			print("\n\n\nWARNING!!!: --binwidth ignored since --nbins supersedes it.")
		nbins = int(round(options.nbins))
	else:
		nbins = int(round( old_div((max(data) - min(data)), width) ))
	
	print("\nAnd the number of bins n = ( max(data) - min(data) ) / width will thus be", nbins)
	nbins = int(round(nbins))
	print("rounding to", nbins)
	
	statistics.append( 'bins=' + str( nbins ) + ' , binwidth=' + str( width ) + '\n')
	print("statistics are", statistics)
	
	if not nbins:

		print("WARNING: nins=0, which means max and min intensity are the same, which probably means all intensities are zero. Defaulting nbins to number of partilces.")
		nbins = len(data)
			
	stem,extension = os.path.splitext(os.path.basename(options.data))
	statsfile = stem + '_stats.txt'

	with open(options.path + '/' + statsfile,'w') as f: f.writelines(statistics)

	return nbins


if __name__ == '__main__':
	main()
