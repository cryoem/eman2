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
matplotlib.use('Agg',warn=False)

		 
import matplotlib.pyplot as plt
import pylab
from pylab import *

import sys
import numpy
import math	 
		 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """This program can be used to fit a polynomial to a curve, fit the minima, fit the maxima, smooth curves, exclude points from a plot, etc, 
				and save the corresponding plots as .png image files. 
				The raw and fitted curves can be saved to the same figure, or as separate figures."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--smooth",action="store_true", help="Smooth out FSC curves.", default=False)
	parser.add_argument("--excludepoints",type=str, help="""Exclude single points, for example --excludepoints=1,3,5 or a range of points, e.g., --excludepoints=1-10
															Recall that, in python, the convention is that the first element is 0 and the upper end is NOT inclded.
															Therefore, to delete the first 10 points (0 through 9) you would provide --excludepoints=0,10 """, default='')
	parser.add_argument("--legend",action="store_true", help="Smooth out FSC curves.", default=False)
	parser.add_argument("--logplot",action="store_true", help="Smooth out FSC curves.", default=False)

	parser.add_argument("--xaxislabel",type=str, default='X', help="""Label for X axis. Recall that special characters such as spaces and strings need to be 'escaped' 
																	with a backslash. For example, to have the X axis label be 'Box size (pixels)' you would need to enter
																	--xaxislabel=box\ size\ \(pixels\)' ; that is, with a backslash before each space and each parenthesis.""")
	parser.add_argument("--yaxislabel",type=str, default='Y', help="""Label for Y axis. Recall that special characters such as spaces and strings need to be 'escaped' 
																	with a backslash. For example, to have the Y axis label be 'Time (s)' you would need to enter
																	--yaxislabel=Time\ \(s\)' ; that is, with a backslash before each space and each parenthesis.""")
	parser.add_argument("--colorlessplot",action="store_true", help="Smooth out FSC curves.", default=False)
	parser.add_argument("--polyfitdeg",type=int, help="Degree of the polynomial to fit to the curve provided.", default=None)
	parser.add_argument("--curves",type=str, help="""Provide .txt. files separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...
													The files should contain x y values separated by A SINGLE SPACE.""", default=None)
	
	parser.add_argument("--path",type=str,default='e2plots',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'e2plots'; for example, e2plots_02 will be the directory by default if 'e2plots_01' already exists.")
											
	parser.add_argument("--fitmaxima",action="store_true",help="Fit descending maxima (as in an FSC plot)", default=False)
	parser.add_argument("--fitminima",action="store_true",help="Fit ascending minima (as in a cuda time-test plot", default=False)
	parser.add_argument("--ascending",action="store_true",help="Necessary only if fitting ASCENDING maxima", default=False)
	parser.add_argument("--descending",action="store_true",help="Necessary only if fitting DESCENDING minima", default=False)
	parser.add_argument("--savetxt",action="store_true",help="Will save a .txt file for polynomial fitted data, fitted minima, maxima, smoothed data, etc.", default=False)

	parser.add_argument("--singleplot",action="store_true",help="All graphs for THE SAME curve (minima, maxima, raw data, etc) will be on the same plot.", default=False)
	parser.add_argument("--ID",type=str, help="""Descriptive tag to identify/label a certain run of the program.""",default='')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	#parser.add_argument("--singlefigure",action="store_true",help="""All PLOTS for DIFFERENT curves will be on the same fgiure.
	#															THis option works ONLY for ONE mode of plotting:
	#															Either the raw data of all your curves, or the fitted maxima for all the curves,
	#															or the fitted minima for all the curves, or the polynomial fit for all the curves.""", default=False)

	(options, args) = parser.parse_args()

	curves = options.curves
	curves = curves.split(',')
	
	rootpath = os.getcwd()
	
	options = makepath(options,rootpath)
	
	raws=[]
	maxs=[]
	mins=[]
	polys=[]
	completes=[]
	for name in curves:
		print "Found this curve to parse", name
		rawaxes = curveparser(name,options)
		logy=[]
		
		if options.logplot:
			logy = [ math.log10(y) for y in rawaxes[1] ]	
			rawaxes[1] = logy
				
			if options.savetxt:
				logname=name.replace('.txt','_LOG.txt')
				textwriter(rawaxes[0],logy,options,logname)
		
		rawdata = {'name':name,'rawaxes':rawaxes}
		raws.append(rawdata)
		
		maxaxes=''
		minaxes=''
		polyaxes=''
		
		xaxis=rawaxes[0]
		yaxis=rawaxes[1]
		if options.fitmaxima:
			maxaxes=fitmaxima(xaxis,yaxis,options)
			maxdata = {'name':name,'maxaxes':maxaxes}
			maxs.append(maxdata)
			
			if options.savetxt:
				maxname=name.replace('.txt','_MAX.txt')
				textwriter(maxaxes[0],maxaxes[1],options,maxname)
			
			
		if options.fitminima:
			minaxes=fitminima(xaxis,yaxis,options)
			mindata = {'name':name,'minaxes':minaxes}
			mins.append(mindata)
			
			if options.savetxt:
				minname=name.replace('.txt','_MIN.txt')
				textwriter(minaxes[0],minaxes[1],options,minname)
	
		if options.polyfitdeg:
			polyaxes=polyfit(xaxis,yaxis,options)
			polydata = {'name':name,'polyaxes':polyaxes}
			polys.append(polydata)			
			
			if options.savetxt:
				polyname=name.replace('.txt','_MIN.txt')
				textwriter(polyaxes[0],polyaxes[1],options,polyname)
			
		completedata={'name':name,'rawaxes':rawaxes,'maxaxes':maxaxes,'minaxes':minaxes,'polyaxes':polyaxes}
		completes.append(completedata)		
	
	
	if options.singleplot:
		
		singleplotter(raws,'raw',options)
		singleplotter(completes,'complete',options)
		
		if options.fitmaxima:
			singleplotter(maxs,'max',options)
		
		if options.fitminima:
			singleplotter(mins,'min',options)
		
		if options.polyfitdeg:
			singleplotter(polys,'poly',options)		
		
	logger = E2init(sys.argv, options.ppid)

	E2end(logger)
	return

def makepath(options,rootpath):
	
	files=os.listdir(rootpath)

	while options.path in files:
		if '_' not in options.path:
			options.path = options.path + '_00'
		else:
			jobtag=''
			components=options.path.split('_')
			if components[-1].isdigit():
				components[-1] = str(int(components[-1])+1).zfill(2)
			else:
				components.append('00')
						
			options.path = '_'.join(components)

	if options.path not in files:
		os.system('mkdir ' + options.path)
	return options


def singleplotter(data,plottype,options):
	axeskey = plottype + 'axes'
	
	name=''
	if 'complete' not in plottype:
		for case in data:
			plotter(case[axeskey][0],case[axeskey][1],options)	
		
		if options.ID:
			options.ID += '_'
		name = options.ID + 'singlePlot_' + plottype + '.png'
		if options.logplot:
			name = name.replace('.png','_LOG.png')
		
		plt.savefig(options.path + '/' + name)	
		plt.clf()
		
	else:
		for case in data:
			for key in case:
				if 'name' not in key:
					#print "Key to plot is", key
					#print "And case[key] should be x,y data, see", case[key]
					if len(case[key]) > 1:
						plotter(case[key][0],case[key][1],options)
			name = case['name'].split('.')[0] + '_singlePlot.png'
			
			if options.logplot:
				name = name.replace('.png','_LOG.png')
			
			plt.savefig(options.path + '/' + name)	
			plt.clf()	
	return


'''
Parse values for FILE containing curve data
'''
def curveparser(F,options):		
	print "Reading this file now", F
	xaxis=[]
	yaxis=[]
	
	f=open(F,'r')
	lines = f.readlines()
	f.close()
	
	for line in lines:
		x=line.split()[0]
		xaxis.append( float(x) )
		y=line.split()[-1].replace('\n','')
		yaxis.append( float(y) )
	
	if options.excludepoints:
		todelvals = options.excludepoints.split(',')
		if '-' in options.excludepoints:
			todel = options.excludepoints.split('-')
			start = todel[0]
			end = todel[1]
			todelvals = [i for i in xrange(start,end)]
			todelvals.reverse()
		
		for i in todelvals:
			del(xaxis[ int(i) ])
			del(yaxis[ int(i) ])	
			
	return[xaxis,yaxis]


def textwriter(xdata,ydata,options,name):
	if len(xdata) == 0 or len(ydata) ==0:
		print "ERROR: Attempting to write an empty text file!"
		sys.exit()
	
	filename=options.path + '/' + name
	
	print "I am in the text writer for this file", filename
	
	f=open(filename,'w')
	lines=[]
	for i in range(len(xdata)):
		line2write = str(xdata[i]) + ' ' + str(ydata[i])+'\n'
		#print "THe line to write is"
		lines.append(line2write)
	
	f.writelines(lines)
	f.close()

	return()

	
def plotter(xaxis,yaxis,options,mode,colorless=0,legend=''):
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
	
	yaxislabel = options.yaxislabel
	xaxislabel = options.yaxislabel
	
	if not yaxislabel:
		yaxislabel = "Y"
	if not xaxislabel:
		xaxislabel = "X"
	
	if options.logplot and yaxislabel:
		if "log" not in yaxislabel and "LOG" not in yaxislabel:
			yaxislabel = 'LOG(' + yaxislabel + ')'	
	
	pylab.xlabel(options.xaxislabel, fontsize=14, fontweight='bold')
	pylab.ylabel(options.yaxislabel, fontsize=14, fontweight='bold')

	#pylab.ylim([-1,ylimvalmax+10])
	#xmin = min(xaxis)
	
	#pylab.xlim([-1,float(max(xaxis))+10])
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	#ax.axes.get_xaxis().set_visible(True)
	#ax.axes.get_yaxis().set_visible(True)
	
	#xmin, xmax = ax.get_xaxis().get_view_interval()
	#ymin, ymax = ax.get_yaxis().get_view_interval()
	
	#ax.add_artist(Line2D((xmin, xmax+10), (ymin, ymin), color='k', linewidth=4))
	#ax.add_artist(Line2D((xmin, xmin), (ymin, ymax+10), color='k', linewidth=4))
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)
	
	markernum=''
	idee=''
	
	#print "BOLD IS ON!"
	LW=3
	if not markernum:
		LW=2
		
	if colorless:
		#if not yminnonconvex:
		#	
		#	print "in colorless plot, linest is", linest
		if mode not 'scatter':
			linest='-'
			plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,color='k',zorder=0,label=idee)
		else:
			plt.scatter(xaxis,yaxis,marker='x',edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
			
		#	if idee and options.legend:
		#		#print "Idee is", idee
		#		legend(loc='upper left')
		#elif yminnonconvex:
		#plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,color='k',zorder=0,label=idee)
		#plt.scatter(xaxis,yaxis,marker='x',edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
		
		if idee and legend:
			print "Idee is", idee
			legend(loc='upper left')
		
		#if mark:
		#	plt.scatter(xaxis,yaxis,marker=mark,edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
	
	else:
		#if not yminnonconvex:
		#	print "I did NOT receive yminnonxonvex"
		#	
		
		if mode not scatter:
			plt.plot(xaxis, yaxis, linewidth=LW,alpha=1,zorder=0,label=idee)
		else:
			plt.scatter(xaxis,yaxis,marker='x',alpha=1,zorder=1,s=40,linewidth=2)
				
		#	
		#	if idee and options.legend:
		#		print "Idee is", idee
		#		legend(loc='upper left')
		#elif yminnonconvex:
		#	print "I DID receive yminnonxonvex"
		#	plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,zorder=0,label=idee)
		
		if idee and legend:
			print "Idee is", idee
			legend(loc='upper left')
			#plt.scatter(xaxis,yaxis,marker='x',alpha=0.5,zorder=1,s=40,linewidth=2)
			
		#if mark:
		#	plt.scatter(xaxis,yaxis,marker=mark,alpha=0.5,zorder=1,s=40,linewidth=2)
	return	
	

def polyfit(xaxis,yaxis,options):
	print "I have entered polyfit"
	polycoeffs = numpy.polyfit(xaxis, yaxis, options.polyfitdeg)
	yfit = numpy.polyval(polycoeffs, xaxis)
	print "The fitted y values are", yfit

	return[xaxis,yfit]


def fitmaxima(xaxis,yaxis,options):
	pass
	xaxismax=''
	yaxismax=''
	return[xaxismax,yaxismax]
	
	
'''
FUNCTION TO DETERMINE MINIMA throughout a curve, from l
'''
def fitminima(sizes,vals,options):
	print "\n\nI have entered minima!!!!!!!!!!!!!!!!!!!!!\n\n"
	print "And the last value is", vals[-1]
	
	
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
			if i < len(vals) - 3:									#All points can be handled the same way except the last three.
				if vals[i] < vals[-1]:								#And only those smaller than the last value should be appended.
					#print "The value is smaller than the last, and will be appended, see", vals[i],vals[-1]
					#print "Types are", type(vals[i]),type(vals[-1])
					sizesmin.append(sizes[i])							
					valsmin.append(vals[i])
			
		if i == len(vals) - 3:			#For the last three points, you scan all possible combinations to add only those that fit the criteria 
											#for ascending minima.
			if vals[i] < vals[i+1]:
			
				if vals[i] < vals [i+2]:
					sizesmin.append(sizes[i])
					valsmin.append(vals[i])				
					
					if vals[i+1] < vals[i+2]:
						sizesmin.append(sizes[i+1])
						valsmin.append(vals[i+1])
					
						sizesmin.append(sizes[i+2])
						valsmin.append(vals[i+2])
					
					elif vals[i+1] > vals[i+2]:
						sizesmin.append(sizes[i+2])
						valsmin.append(vals[i+2])
								
				
				elif vals[i] > vals[i+2]:
					sizesmin.append(sizes[i+2])
					valsmin.append(vals[i+2])
			
			elif vals[i] > vals[i+1]:
				print "Third to last is bigger"
				if vals[i+1] < vals[i+2]:
					print "Second to last is smaller than last and will be appended"
					sizesmin.append(sizes[i+1])
					valsmin.append(vals[i+1])
				elif vals[i+1] > vals[1+2]:
					print "Last is smaller than second to last and third to last; thus, last will be appended"
					sizesmin.append(sizes[i+2])
					valsmin.append(vals[i+2])

	#'''
	#Check backwards if the last element is lower than the last minima appended
	#'''
	#auxn=0
	#while auxn=0:
	#	if valsmin[-1] > vals[-1]:
	#		del(valsmin[-1])
	#		del(sizesmin[-1])
	#		valsmin[-1] = vals[-1]
	#		sizesmin[-1] = sizes[-1]
	#							
	#		#print "I have appended this box, value", sizes[i], vals[i]
	'''
	The first element does not have to be included by default. It should be removed, and the first minima should be the first plotted element
	'''
	if valsmin[0] > valsmin[1]:
		del(valsmin[0])
		del(sizesmin[0])
	return[sizesmin,valsmin]

	
def smooth(yaxis):

	yaxisnonconvex=Util.nonconvex(yaxis,0)

	return yaxisnonconvex
	

if '__main__' == __name__:
	main()