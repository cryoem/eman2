#!/usr/bin/env python
'''
====================
Author: Jesus Galaz-Montoya - 2017, Last update: 19/July/2017
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

import matplotlib
matplotlib.use('Agg',warn=False)

import matplotlib.pyplot as plt
import pylab

#import matplotlib.colors as mcol
#import matplotlib.cm as cm

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
	
	parser.add_argument("--data", type=str, default='', help="""default=None (not used). Text file(s) with two column of values mean to be plotted on the X and Y axes. If supplying multiple files, separate them by commas.""")
	parser.add_argument("--datax", type=str, default='', help="""default=None (not used). Text file(s) with a single column of values meant to be plotted on the X axis. If not provided, the X axis will go from 0 to n, where n is the number of values in --datay. If supplying multiple files, separate them by commas (the number of files for --datax and --datay must be the same).""")
	parser.add_argument("--datay", type=str, default='', help="""default=None (not used). Text file(s) with a single column of values meant to be plotted on the  Y axis. If not provided, the Y axis will go from 0 to n, where n is the number of values in --datax. If supplying multiple files, separate them by commas (the number of files for --datax and --datay must be the same).""")

	parser.add_argument("--highresolution", action='store_true', default=False, help="""default=False. If on, this option will enforce writing high-resolution plots (dpi 300), ready for publication, as opposed to lower resolution plots which take up less space on your computer (dpi 150).""")

	parser.add_argument("--individualplots", action='store_true', default=False, help="""default=False. in addition to plotting the motion for all frames in all images in a single plot, generate individual plots per image""")
	
	parser.add_argument("--labelxaxis", type=str,default='x',help="""Default=x. Label for x axis (specify units through --unitsx.""")
	parser.add_argument("--labelyaxis", type=str,default='y',help="""Default=y. Label for y axis (specify units through --unitsy.""")
	parser.add_argument("--labeltitle", type=str,default='',help="""Default=None. Title for figure.""")

	parser.add_argument("--maxx",type=float,default=0.0,help="""Default=None. Maximum value to plot in X. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--maxy",type=float,default=0.0,help="""Default=None. Maximum value to plot in Y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--minx",type=float,default=0.0,help="""Default=None. Minimum value to plot in X. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	parser.add_argument("--miny",type=float,default=0.0,help="""Default=None. Minimum value to plot in Y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	
	parser.add_argument("--nocolor", action='store_true', default=False, help="""Default=False. Plots are colored, by default; don't be cheap; clear communication and representation pays off; or consider publishing in online open source journals.""")

	parser.add_argument("--normalize", action='store_true', default=False, help="""Default=False. This option will normalize all plots to be scaled between 0 and 1.""") 
	
	parser.add_argument("--outputtag",type=str,default='plotfig',help="""Default=plotfig. String common to all automatically generated output files. For example, --outputtag=myplot will generate myplot1.png, myplot2.png, ..., myplotN.png""")

	parser.add_argument("--path", type=str,default='plotfig',help="""Defaault=de_plots. Name of the directory where to store the output results.""")
	parser.add_argument("--ppid", type=int, default=-1,help="Default=-1. Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--unitsx", type=str,default='AU',help="""Default=AU (arbitrary units). Units for the x axis.'microns' or 'mu' and 'angstroms' or 'A' (and '1/angstroms' or '1/A') will be replaced by the appropriate symbol.""")
	parser.add_argument("--unitsy", type=str,default='AU',help="""Default=AU (arbitrary units). Units for the y axis.'microns' or 'mu' and 'angstroms' or 'A' (and '1/angstroms' or '1/A')  will be replaced by the appropriate symbol.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	datafiles=[]
	if options.data:
		if options.datax or options.datay:
			print "\n(e2plotfig)(main) ERROR: provide --data OR --datax AND/OR --datay."
			sys.exit(1)
	elif not options.data:
		if not options.datax and not options.datay:
			if args:
				datafiles=args
			elif not args:
				print "\n(e2plotfig)(main) ERROR: provide at least one of --data, --datax, or --datay."
				sys.exit(1)


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
		
		for f in datafiles:
			
			xaxis=[]
			yaxis=[]
			with open( f ) as datafile: 
				lines=datafile.readlines()
				#if lines:
				xaxis = [ float(line.replace('\n','').split()[0]) for line in lines ]
				yaxis = [ float(line.replace('\n','').split()[1]) for line in lines ]

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

			if options.datay:
				datayfiles=options.datay.split(',')
				if len(dataxfiles) != len(datayfiles):
					print "\n(e2plotfig)(main) ERROR: --datax and --datay must contain the same number of files. Now, nx files=%d, ny files=%d".format(len(dataxfiles),len(datayfiles))
					sys.exit(1)

			k=0
			for fx in dataxfiles:
				with open( fx ) as dataxfile: 
					lines=dataxfile.readlines()
					xaxis = [ float(line.replace('\n','').split()[0]) for line in lines ]
					xaxes.update({k:xaxis})
					if not options.datay:
						yaxis = range(len(xaxis))
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
						xaxis = range(len(yaxis))
						xaxes.update({k:xaxis})

						datadict.update({k:[xaxis,yaxis]})
					k+=1

	from EMAN2_utils import makepath
	options = makepath(options)

	fig,ax=resetplot()

	plotdata(options,datadict)

	E2end(logger)


	return


def normalize(data):
	dminusmin = [d-min(data) for d in data]
	dminusminovermax = [d/max(dminusmin) for d in dminusmin]

	return dminusminovermax


def resetplot(figsize=None):
	plt.clf()

	fig = plt.figure()
	if figsize:
		fig = plt.figure(figsize=(10, 6))
  	ax = fig.add_subplot(1,1,1)

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

	N = len(data)
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

	import string
	markers=string.printable
	marker=''

	#markers=["*","o","x",">","<","|","z","m",]

	#mark = ''
	
	resolution=150
	if options.highresolution:
		resolution=300

	fig,ax = resetplot()
	#cpick = resetcolorbar(options.lowestangle,options.highestangle)
	ndata = len(data)

	print "\nllllllen(data) {}".format(len(data))
	#colorbar = False
	#altxaxis = None
	#colorstart = options.lowestangle
	#colorstep = options.tiltstep

	for k in data:
		#if count==ndata-1:
		#	colorbar=True
		color = RGB_tuples[k]
		if options.nocolor:
			color = 'k'
			marker=markers[k]
			#try:
			#	mark = markers[k]
			#except:
			#	mark = markers[-1]
		
		if options.verbose > 9:
			print "\n(e2plotfig)(plotdata) plotting dataset n = {}".format(k)
			print "color is {}".format(color)
			print "marker is {}".format(marker)
		
		plotfig( options, fig, ax, data[k][0], data[k][1], color, marker )		

	filetosave = options.path + '/' + options.outputtag + '.png'
	fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')# transparent=True) #, bbox_extra_artists=(lgd,), bbox_inches='tight'
	plt.close('fig')
	print "\n(e2plotfig)(plotdata) saving figure {}".format(filetosave)

	if options.individualplots:
		#colorstep=0
		for k in data:
			#colorbar=False
			fig,ax = resetplot()
			print "\n(e2plotfig)(plotdata) plotting individual plot for dataset n = {}".format(k)
			plotfig( options, fig, ax, data[k][0], data[k][1], k )	
			filetosave = options.path + '/' + options.outputtag + str( k ).zfill( len(str(ndata)))+'.png'
			fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')
			plt.close('fig')

	return


def plotfig( options, fig, ax, datax, datay, colorthis='k', markerthis='' )	:
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both', reset=False, which='both', length=8, width=3)

	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	#matplotlib.rc('font', weight='bold') 

	#if options.miny and options.maxy:
	#	ax.set_ylim(options.miny,options.maxy)
	#	maxr = int(round(math.sqrt(options.maxy*options.maxy + options.miny*options.miny)))	

	n=len(datay)
	#xaxis=range(n)
	#if altxaxis:
	#	xaxis=altxaxis


	
	#colorthis = 'k'
	#if colorstart and colorstep:
	#	colorthis = cpick.to_rgba( colorstart + kplot*colorstep)
	
	ax.plot( datax, datay, linewidth=2, marker=markerthis, markersize=5, color=colorthis, label='') #, alpha=0.75)

	#if errors:
	#	lines = {'linestyle': 'None'}
	#	plt.rc('lines', **lines)
	#	plt.errorbar(xaxis,values,errors,markersize=8,linewidth=1,fmt='',marker='o',color='k',markerfacecolor=None,markeredgecolor='k',capsize=5, capthick=1)

	ax.set_title(options.labeltitle, fontsize=16, fontweight='bold')
	ax.set_xlabel(options.labelxaxis + ' ' + options.unitsx, fontsize=16, fontweight='bold')
	ax.set_ylabel(options.labelyaxis + ' ' + options.unitsy, fontsize=16, fontweight='bold')

	
	#if colorbar:
	#	cbr=plt.colorbar(cpick)
	#	cbr.set_label("Tilt angle (degrees)",fontsize=18, fontweight='bold')
	
	if not options.individualplots:
		plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
	
	return


if __name__ == '__main__':
	main()