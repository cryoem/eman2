#!/usr/bin/env python

#
# Author: Jesus Galaz, 05/June/2013
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

import os
from EMAN2 import *
from time import time
		 
import matplotlib
matplotlib.use('Agg')
		 
import matplotlib.pyplot as plt
import pylab
from pylab import *

import sys
import numpy
import math	 
		 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """This program can be used to fit a polynomial to a curve, fit the minima, fit the maxima, smooth curves, etc, 
				and save the corresponding plots as a .png image file. The raw and fitted curves can be saved to the same figure, or as separate figures."""
	
	parser.add_argument("--smooth",action="store_true", help="Smooth out FSC curves.", default=False)
	parser.add_argument("--polyfitdeg",type=int, help="Degree of the polynomial to fit to the curve provided.", default=None)
	parser.add_argument("--curves",type=str, help="""Provide .txt. files separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...
													The files should contain x y values separated by A SINGLE SPACE.""", default=None)
													
	parser.add_argument("--fitmaxima",action="store_ture",help="Fit descending maxima (as in an FSC plot)", default=False)
	parser.add_argument("--fitminima",action="store_true",help="Fit ascending minima", default=False)
	parser.add_argument("--ascending",action="store_true",help="Necessary only if fitting ASCENDING maxima", default=False)
	parser.add_argument("--descending",action="store_true",help="Necessary only if fitting DESCENDING minima", default=False)

	parser.add_argument("--singleplot",action="store_true",help="All graphs for THE SAME curve (minima, maxima, raw data, etc) will be on the same plot.", default=False)
	#parser.add_argument("--singlefigure",action="store_true",help="""All PLOTS for DIFFERENT curves will be on the same fgiure.
	#															THis option works ONLY for ONE mode of plotting:
	#															Either the raw data of all your curves, or the fitted maxima for all the curves,
	#															or the fitted minima for all the curves, or the polynomial fit for all the curves.""", default=False)

	(options, args) = parser.parse_args()

	curves = options.curves
	curves = curves.split(',')
	
	raws=[]
	maxs=[]
	mins=[]
	polys=[]
	completes=[]
	for name in curves:
		print "Found this curve to parse", name
		rawaxes = curveparser(name)
		rawdata = {'name':name,'rawaxes':rawaxes}
		raws.append(rawdata)
		
		maxaxes=''
		minaxes=''
		polyaxes=''
		
		if options.fitmaxima:
			maxaxes=fitmaxima(xaxis,yaxis,options)
			maxdata = {'name':name,'maxaxes':maxaxes}
			maxs.append(maxdata)
			
		if options.fitminima:
			minaxes=fitminima(xaxis,yaxis,options)
			mindata = {'name':name,'minaxes':minaxes}
			mins.append(mindata)
	
		if options.polyfit:
			polyaxes=polyfit(xaxis,yaxis,options)
			polydata = {'name':name,'polyaxes':polyaxes}
			polys.append(polydata)			
	
		completedata={'name':name,'rawaxes':rawaxes,'maxaxes':maxaxes,'minaxes':minaxes,'polyaxes':polyaxes}
		completes.append(completedata)		
	
	
	if options.singleplot:
		
		singleplotter(raws,'raw')
		singleplotter(completes,'complete')
		
		if options.fitmaxima:
			sigleplotter(maxs,'max')
		
		if options.fitmaxima:
			sigleplotter(mins,'min')
		
		if options.polyfit:
			sigleplotter(polys,'poly')		
		
	logger = E2init(sys.argv, options.ppid)

	E2end(logger)
	return


def singleplotter(data,plottype):
	axeskey = plottype + 'axes'
	
	name=''
	if 'complete' not in plottype:
		for case in data:
			plotter(data[axeskey])	
		
		if options.ID:
			options.ID += '_'
		name = options.ID + 'singlePlot_' + plottype + '.png'
		
		plt.savefig(name)	
		plt.clf()
		
	else:
		for case in data:
			for key in case:
				plotter(case[key])
			name
	
		
	return


'''
Parse values for FILE containing curve data
'''
def curveparser(F):		
	print "Reading this file now", F
	xaxis=[]
	yaxis=[]
	
	f=open(F,'r')
	lines = f.readlines()
	f.close()
	
	for line in lines:
		x=line.split()[0]
		xaxis.append(int(size))
		y=line.split()[-1].replace('\n','')
		yaxis.append(float(value))
		
	return(xaxis,yaxis)

	
def plotter(xaxis,yaxis,options)
	'''
	FORMAT AXES
	'''
	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	#ax=plt.axes(frameon=False)
	ax=plt.axes()
	pylab.rc("axes", linewidth=2.0)
	pylab.xlabel('X Axis', fontsize=14, fontweight='bold')
	pylab.ylabel('Y Axis', fontsize=14, fontweight='bold')

	#pylab.ylim([-1,ylimvalmax+10])
	pylab.xlim([-1,max(xaxis)+10])
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	#ax.axes.get_xaxis().set_visible(True)
	#ax.axes.get_yaxis().set_visible(True)
	
	#xmin, xmax = ax.get_xaxis().get_view_interval()
	#ymin, ymax = ax.get_yaxis().get_view_interval()
	
	#ax.add_artist(Line2D((xmin, xmax+10), (ymin, ymin), color='k', linewidth=4))
	#ax.add_artist(Line2D((xmin, xmin), (ymin, ymax+10), color='k', linewidth=4))
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)
	
	#print "BOLD IS ON!"
	LW=3
	if not markernum:
		LW=2
		
	if options.colorlessplot:
		
		if not yminnonconvex:
			
			print "in colorless plot, linest is", linest
			plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,color='k',zorder=0,label=idee)
			if idee and options.legend:
				#print "Idee is", idee
				legend(loc='upper left')
		elif yminnonconvex:
			plt.plot(xaxis, yminnonconvex, linewidth=LW,linestyle=linest,alpha=1,color='k',zorder=0,label=idee)
			plt.scatter(xaxis,yaxis,marker='x',edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
			if idee and options.legend:
				print "Idee is", idee
				legend(loc='upper left')
			
		if mark:
			plt.scatter(xaxis,yaxis,marker=mark,edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
	else:
		if not yminnonconvex:
			print "I did NOT receive yminnonxonvex"
			
			plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,zorder=0,label=idee)
			
			if idee and options.legend:
				print "Idee is", idee
				legend(loc='upper left')
		elif yminnonconvex:
			print "I DID receive yminnonxonvex"
			plt.plot(xaxis, yminnonconvex, linewidth=LW,linestyle=linest,alpha=1,zorder=0,label=idee)
			
			if idee and options.legend:
				print "Idee is", idee
				legend(loc='upper left')
			plt.scatter(xaxis,yaxis,marker='x',alpha=0.5,zorder=1,s=40,linewidth=2)
			
		if mark:
			plt.scatter(xaxis,yaxis,marker=mark,alpha=0.5,zorder=1,s=40,linewidth=2)
	return	
	

def polyfit(xaxis,yaxis,options)

	polycoeffs = numpy.polyfit(xaxis, yaxis, options.polydegree)
	yfit = numpy.polyval(polycoeffs, x)

	return	yfit


def fitmaxima(xaxis,yaxis,options):
	pass
	xaxismax=''
	yaxismax=''
	return(xaxismax,yaxismax)
	
	
'''
FUNCTION TO DETERMINE MINIMA throughout a curve, from l
'''
def fitminima(sizes,vals,options):
	print "\n\nI have entered minima!!!!!!!!!!!!!!!!!!!!!\n\n"
	
	
	#finalvals = []
	#finalsizes = []
	
	
	sizesmin=[sizes[0]]
	valsmin=[vals[0]]
	
	for i in xrange(0,len(vals) - 1 - 1 ):
		aux = 0 
		for j in xrange(i+1,len(vals) - 1 ):
			if vals[j]< vals[i]:
				aux=0
				#print "Because a downstream value is lower, I will break the loop, see", vals[j],vals[i]
				break
			else:
				aux=1
		if aux==1:
			if i < len(vals) - 3:				#All points can be handled the same way except the last three.
				sizesmin.append(sizes[i])
				valsmin.append(vals[i])
			
			elif i == len(vals) - 3:			#For the last three, you scan all possible combinations 
				
				if vals[i] < vals[i+1]:
				
					if vals[i] < vals [i+2]:
						sizesmin.append(sizes[i])
						valsmin.append(valsmin[i])				
						
						if vals[i+1] < vals[i+2]:
							sizesmin.append(sizes[i+1])
							valsmin.append(valsmin[i+1])
						
							sizesmin.append(sizes[i+2])
							valsmin.append(valsmin[i+2])
						elif vals[i+1] > vals[i+2]:
							sizesmin.append(sizes[i+2])
							valsmin.append(valsmin[i+2])
									
					
					elif vals[i] > vals[i+2]:
						sizesmin.append(sizes[i+2])
						valsmin.append(valsmin[i+2])
				
				elif vals[i] > vals[i+1]:
					if vals[i+1] < vals[i+2]:
						sizesmin.append(sizes[i+1])
						valsmin.append(valsmin[i+1])
					elif valis[i+1] > vals[1+2]:
						sizesmin.append(sizes[i+2])
						valsmin.append(valsmin[i+2])
								
			#print "I have appended this box, value", sizes[i], vals[i]


	return(sizesmin,valsmin)

	
def smooth(yaxis):

	yaxisnonconvex=Util.nonconvex(yaxis,0)

	return yaxisnonconvex
	

if '__main__' == __name__:
	main()