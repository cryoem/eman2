#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 28/April/2013
# Copyright (c) 2013 Baylor College of Medicine
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
	usage = """WARNING: Experimental script. 
	Tests cpu vs gpu speed on subtomogram alignment, for one CCF between two boxes of noise at default orientation (no rotations), 
	with and without alignment overhead, and plots the results against boxsize.
	Full alignment scanning an entire icosahedral unit can be tested. The number of orientations visited depends on --coarsestep and --finestep,
	and the way these parameters are used by e2spt_classaverage.py"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--cpu", action='store_true', help="Will test SPT alignment using CPU.",default=False)
	parser.add_argument("--gpu", action='store_true', help="Will test SPT alignment using GPU.",default=False)
	parser.add_argument("--setcudadevice", type=str, help="The number of the cuda device you want to use (from 0, to the total amount of devices there are)",default=None)	
	
	parser.add_argument("--profile", action='store_true', help="Will profile the internal call to e2spt_classaverage.py to see what processes are consuming the most time.",default=False)
	parser.add_argument("--eman2dir", type=str, help="If profiling is on, you must provide the full path to your EMAN2 directory WITHOUT. For example, /Users/jgalaz/EMAN2",default='')

	parser.add_argument("--test", action='store_true', help="Will run a quick tests using a few box sizes.",default=False)
	parser.add_argument("--medium", action='store_true', help="Will test boxsizes in multiples of 10 between 10 and 240.",default=False)
	parser.add_argument("--extensive", action='store_true', help="Will test EVERY box size between 12 and 256.",default=False)
	parser.add_argument("--boxsizes", type=str, help="""For a SINGLE AND CONTINUOUS range, separate the first and the last boxsizes to sample with a hyphen; for example, 16-48. 
														Following python's conventions, the first index is 0, and the upper end of any range is not included.
														For specific boxes, separate the values with a comma, or provide a .txt file with a list of box sizes.""",default='')
	#parser.add_argument("--boxsizelist", type=str, help="Will test SPT alignment for a comma separated list of boxsizes.",default='')

	parser.add_argument("--onefulloff", action='store_true', help="""Will test time for one CCF (between two boxes with noise) with all
									 processing overhead for EVERY box size between 12 and 256.""",default=False)
	parser.add_argument("--oneccfoff", action='store_true', help="Will test time for one CCF (between two boxes with noise)  without overhead (e2spt_classaverage.py isn't called), for EVERY box size between 12 and 256.",default=False)
	parser.add_argument("--rotonlyoff", action='store_true', help="""Will test time for CCFs (between two boxes with noise) for as many orientations fit in an icosahedral unit 
																	as determined by --coarsestep and --finestep, without any overhead (e2spt_classaverage.py isn't called), for EVERY box size between 12 and 256.""",default=False)
	
	parser.add_argument("--oneicosoff", action='store_true', help="""Will test alignment time for however many orientations fit in one icosahedral unit, 
																	depending on --coarsestep and --finestep, for EVERY box size between 12 and 256.""",default=False)
	
	parser.add_argument("--ID", type=str, help="Tag files generated on a particular computer.",default='')
	parser.add_argument("--sym", type=str, help="Confine alignment time test(s) to scanning an asymmetric unit with --coarsestep and --finestep sampling.",default='c1')

	parser.add_argument("--coarsestep", type=int, help="Step size for coarse alignment.",default=30)
	parser.add_argument("--finestep", type=int, help="Step size for fine alignment.",default=15)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default='sptCudaTest',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptCudaTest'; for example, sptCudaTest_01 will be the directory by default if 'sptCudaTest_00' already exists.")
	parser.add_argument("--noplot", action='store_true', help="To run jobs on a cluster or other interfaces that don't support plotting, turn this option on.",default=False)
	
	parser.add_argument("--plotonly",type=str, help="""If you already have the time variation with boxsize (for a particular cpu or gpu) in a text file in rows of 'boxsize,time', 
												provide the txt file(s) separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...""", default=None)

	parser.add_argument("--logplot", action='store_true', help="If on, the y-axis of plots will be on log scale.",default=False)
	parser.add_argument("--legend", action='store_true', help="If on, a legend will appear on the plot to denote each curve.",default=False)											
											
	parser.add_argument("--singleplot",type=str, help="Provide the name of the .png file to plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default='')
	
	parser.add_argument("--plotminima",action="store_true", help="Plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default=False)
	parser.add_argument("--colorlessplot",action="store_true", help="Plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default=False)
	
	parser.add_argument("--subset",type=int, help="If --plotonly is on, this is the subset of lines (from 0 to --subset=n) in the files provided that will be kept for plotting.", default=0)
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default='')

	global options

	(options, args) = parser.parse_args()
	
	if options.test and options.extensive:
		print "ERROR: You cannot specify --test and --extensive simultaneously"
		sys.exit()
	
	if options.test and options.medium:
		print "ERROR: You cannot specify --test and --medium simultaneously"
		sys.exit()

	if options.extensive and options.medium:
		print "ERROR: You cannot specify --extensive and --medium simultaneously"
		sys.exit()


	
	if options.profile and not options.eman2dir:
		print "ERROR: --profile is on, yet you haven't provided the path to your EMAN2 directory."
		sys.exit()
	#print "options are", options
	
	logger = E2init(sys.argv, options.ppid)
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	rootpath = os.getcwd()
	options = makepath(options,rootpath)
	
	
	print "Path will be", options.path
	
	'''
	If options.plotonly is off (no files provided) you actually have to measure alignment time....
	'''
	if not options.plotonly:
		retcpu=[]
		retgpu=[]
		if options.ID:
			options.ID = options.ID + '_'
		else:
			options.ID=''
							
		'''
		Compute GPU alignment times and plot them
		'''
		origdir=os.getcwd()
		if options.gpu:
			corg = 'GPU'
			retgpu=doit(corg,options,origdir)
			
			for key in retgpu:
				data = retgpu[key]
				preplot(options,corg,data,key)
			
		'''
		Compute CPU alignment times and plot them
		'''
		if options.cpu:
			corg = 'CPU'
			retcpu=doit(corg,options,origdir)

			for key in retcpu:
				data = retcpu[key]
				preplot(options,corg,data,key)

		'''
		I you have both CPU and GPU times, compute the ratio and plot it
		'''
		if retcpu and retgpu:
			if len(retcpu) == len(retgpu):
				print "In comparing GPU and CPU, the len of each ret data is", len(retgpu), len(retcpu)
				print "In fact, retcpu is", retcpu
				print "In fact, retgpu is", retgpu
				
				corg = 'GPUvsCPU'
				
				for key in retgpu:
					if len(retcpu[key]) == len(retgpu[key]):					
						
						retcpuC = retcpu[key]
						retgpuC = retgpu[key]
						
						gnums = numpy.array(retgpuC[1])
						cnums = numpy.array(retcpuC[1])
						
						sizes = retgpuC[0]
						difs = cnums/gnums
						difs = list(difs)
						
						data=[sizes,difs]
						preplot(options,corg,data,key)
	
					else:
						print "This gpu set does not have the same number of elements as the corresponding cpu set"
			else:
				print "For some sick reason, you don't have the same number of data points for gpu and cpu, therefore, you cannot compare them, see", len(retgpu), len(retcpu)
				sys.exit()
		else:
			return()

	else:
		'''
		If options.plotonly received some files, parse them and plot them
		'''
		
		files=options.plotonly.split(',')
		print "\n\n\n\n$$$$$$$$$$$$$$$$$$Will plot these files", files
		print "$$$$$$$$$$$$$$\n\n\n"
		if options.noplot:
			print "ERROR: You cannot speficy 'plotonly' and 'noplot' at the same time."
			sys.exit()
		
		mastervalues={}
		
		k=0.0
		
		'''
		Parse values for all files
		'''
		for F in files:
			print "Reading with this file now", F
			name=os.path.basename(F).replace('.txt','.png')
			sizes=[]
			valuesforthisfile=[]
			
			f=open(F,'r')
			lines = f.readlines()
			f.close()
			
			for line in lines:
				size=line.split()[0]
				sizes.append(int(size))
				value=line.split()[-1].replace('\n','')
				valuesforthisfile.append(float(value))
			
			mastervalues.update({F:[sizes,valuesforthisfile]})				
		
		print "\n\n\n\n\n\n\n@@@@@@@@@@@@@@\nThe len of mastervalues before FULL is", len(mastervalues)
		print "@@@@@@@@@@@@@@@\n\n\n"
		
		'''
		Plot the full extension values for all files
		'''
		plt.clf()
		k=0
		absmax=0.0
		for F in mastervalues:
			print "\n\n\n\n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!! \nWorking with this file now for FULL plot", F
			print "!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n\n\n\n"
			sizes = mastervalues[F][0]
			valuesforthisfile = mastervalues[F][1]
			
			if float(max(valuesforthisfile)) > absmax:
				absmax = float(max(valuesforthisfile))
				print "New ABSMAX is", absmax
			
			namep = options.singleplot
			if not options.singleplot:
				namep = F.replace('.txt','.png')
			
			print "The name for these FULL plots will be", namep
			markernum=0
			#if options.colorlessplot:
			#	markernum=k
			#	k+=1
			
			idee=F.split('_')[0]
			plotter(sizes,valuesforthisfile,namep,0,0,markernum,0,absmax,idee)
			plt.savefig(options.path + '/' + os.path.basename(namep))
		
			if not options.singleplot:
				plt.clf()
		
		print "\n\n\n\n\n\n\n@@@@@@@@@@@@@@\nThe len of mastervalues before SUB is", len(mastervalues)
		print "@@@@@@@@@@@@@@@\n\n\n"
		k=0
		plt.clf()
		abssubmax=0.0	
		if options.subset:
			print "subset is on!"
			
			for F in mastervalues:
				print "\n\n\n\n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!! \nWorking with this file now for SUB plot", F
				print "!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n\n\n\n"
				sizes = mastervalues[F][0]
				valuesforthisfile = mastervalues[F][1]
		
				sizessub = sizes[:options.subset]
				valuessub = valuesforthisfile[:options.subset]
		
				if float(max(valuessub)) > abssubmax:
					abssubmax = float(max(valuessub))
					print "New ABSSUBMAX is", abssubmax
		
				namepsub = F.replace('.txt','_sub.png')
			
			
				print "The name for these SUB plots will be", namepsub

				markernum=0
				if options.colorlessplot:
					#markernum=k
					k+=1
			
				idee=F.split('_')[0]
				plotter(sizessub,valuessub,namepsub,0,0,markernum,0,abssubmax,idee)
				plt.savefig(options.path + '/' + os.path.basename(namepsub))
		
				if not options.singleplot:
					plt.clf()
	
		
		'''
		Plot only the minima, either for full extension or a subset of the files' values
		'''
		
		print "\n\n\n\n\n\n\n@@@@@@@@@@@@@@\nThe len of mastervalues before MIN is", len(mastervalues)
		print "@@@@@@@@@@@@@@@\n\n\n"
		k=0
		plt.clf()
		absminmax=0.0
		if options.plotminima:
			for F in mastervalues:
				print "\n\n\n\n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!! \nWorking with this file now for MIN plot", F
				print "!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n\n\n\n"	
				if options.colorlessplot:
					#markn=k+1
					pass
	
				sizes = mastervalues[F][0]
				valuesforthisfile = mastervalues[F][1]
	
				ret=minima(sizes,valuesforthisfile)
				sizesmin=ret[0]
				valuesmin=ret[1]
				yminnonconvex=ret[2]
		
				if float(max(valuesmin)) > absminmax:
					absminmax = float(max(valuesmin))
					print "New ABSMINMAX is", absminmax
			
				namepmin=F.replace('.txt','_MIN.png')
				
				print "The name for these MIN plots will be", namepmin

				markernum=0
				if options.colorlessplot:
					#markernum=k
					k+=1
				
				idee=F.split('_')[0]
			
				plotter(sizesmin,valuesmin,namepmin,0,0,markernum,0,absminmax,idee,yminnonconvex)
				plt.savefig(options.path + '/' + os.path.basename(namepmin))
	
				if not options.singleplot:
					plt.clf()
		
		'''
		Plot minima for subset
		'''
		k=0
		plt.clf()
		absminsubmax=0.0
		print "\n\n\n\n\n\n\n@@@@@@@@@@@@@@\nThe len of mastervalues before SUBMIN is", len(mastervalues)
		print "@@@@@@@@@@@@@@@\n\n\n"
		if options.plotminima and options.subset:
			for F in mastervalues:
				print "\n\n\n\n\n\n\n\n\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!! \nWorking with this file now for SUBMIN plot", F
				print "!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n\n\n\n"
				if options.colorlessplot:
					#markn=k+1
					pass
							
				sizes = mastervalues[F][0]
				valuesforthisfile = mastervalues[F][1]
	
				ret=minima(sizes,valuesforthisfile)
				sizesmin=ret[0]
				valuesmin=ret[1]
				yminnonconvex=ret[2]
			
				indx=0
				for i in sizesmin:
					if int(i) > int(options.subset):
						print "i is larger than subset see", i, int(options.subset)
						break
					else:
						print "I have added 1 to indx"
						indx+=1
					
				print "Index is", indx
						
				sizesminsub = sizesmin[:indx]
				valuesminsub = valuesmin[:indx]
				yminnonconvexsub = yminnonconvex[:indx]
			
				if float(max(valuesminsub)) > absminsubmax:
					absminsubmax = float(max(valuesminsub))
					print "New ABSMINSUBMAX is", absminsubmax
			
				namepsubmin = F.replace('.txt','_sub_MIN.png')
				
				print "The name for these FULL plots will be", namepsubmin

		
				markernum=0
				if options.colorlessplot:
					#markernum=k
					k+=1
			
				idee=F.split('_')[0]
				plotter(sizesminsub,valuesminsub,namepsubmin,0,0,markernum,0,absminsubmax,idee,yminnonconvexsub)
				plt.savefig(options.path + '/' + os.path.basename(namepsubmin))
	
				if not options.singleplot:
					plt.clf()
		
	return()


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


def preplot(options,tipo,data,key):
	print "I am in PREPLOT"
	name = options.path + "/" + options.ID + key + "_CS"+ str(options.coarsestep).zfill(2)  + "_FS" + str( options.finestep ).zfill(2) + '_' + tipo + '.png'
	xdata = data[0]				
	ydata = data[1]				
	
	print "Still in preplot"
	if not options.noplot:
		plotter(xdata,ydata,name,options.coarsestep,options.finestep)
		plt.savefig(options.path + '/' + os.path.basename(name))
		plt.clf()
	print "Will call texwriter from preplot"
	textwriter(name,xdata,ydata)
	print "In preplot, plotminima is", options.plotminima
	
	if options.plotminima:
		ret=minima(xdata,ydata)
		xdatamins=ret[0]
		ydatamins=ret[1]
		yminnonconvex=ret[2]
		
		print "In preplot, returned yminnonconvex is", yminnonconvex
		namemin=name.replace('.png','_MIN.png')
		markernum=0
		
		
		if not options.noplot:
			plotter(xdatamins,ydatamins,namemin,options.coarsestep,options.finestep,markernum,0,0,'',yminnonconvex)
			plt.savefig(options.path + '/' + os.path.basename(namemin))
			plt.clf()
		
			if options.colorlessplot:
				#markernum=1
				plotter(xdatamins,ydatamins,namemin,options.coarsestep,options.finestep,markernum,0,0,'',yminnonconvex)
				plt.savefig(options.path + '/' + os.path.basename(namemin.replace('.png','_colorless.png')))
				plt.clf()
		print "Name min to send is", namemin
		textwriter(namemin,xdatamins,ydatamins)

	if options.subset:
		xdatasub=xdata[0:options.subset]
		ydatasub=ydata[0:options.subset]
		namesub=name.replace('.png','_SUB' + str(options.subset).zfill(3) + '.png')
		
		if not options.noplot:
			plotter(xdatasub,ydatasub,namesub,options.coarsestep,options.finestep,markernum,0,0,'',)
			plt.savefig(options.path + '/' + os.path.basename(namesub))
			plt.clf()
			markernum=0
			if options.colorlessplot:
				#markernum=1
				plotter(xdatamins,ydatamins,namesub,options.coarsestep,options.finestep,markernum,0,0,'',yminnonconvex)
				plt.savefig(options.path + '/' + os.path.basename(namesub.replace('.png','_colorless.png')))
				plt.clf()
				
		print "Name sub to send is", namesub
		textwriter(namesub,xdatasub,ydatasub)

		if options.plotminima:
			print "From preplot, I will call minima"
			ret=minima(xdatasub,ydatasub)
			xdatasubmins=ret[0]
			ydatasubmins=ret[1]
			yminnonconvexsub=ret[2]
			print "I finshed calling minima, and it returned this yminnonconvexsub",yminnonconvexsub 
			namesubmin=namesub.replace('.png','_MIN.png')
			markernum=0
			
			
			if not options.noplot:			
				plotter(xdatasubmins,ydatasubmins,namesubmin,options.coarsestep,options.finestep,markernum,0,0,'',yminnonconvexsub)
				plt.savefig(options.path + '/' + os.path.basename(namesubmin))
				plt.clf()
				markernum=0
				if options.colorlessplot:
					#markernum=1
					plotter(xdatasubmins,ydatasubmins,namesubmin,options.coarsestep,options.finestep,markernum,0,0,'',yminnonconvexsub)
					plt.savefig(options.path + '/' + os.path.basename(namesubmin.replace('.png','_colorless.png')))
					plt.clf()
			print "Name sub min to send is", namesubmin
			textwriter(namesubmin,xdatasubmins,ydatasubmins)

	return()


def textwriter(name,xdata,ydata):
	if not xdata or not ydata:
		print "ERROR: Attempting to write an empty text file!"
		sys.exit()
	
	filename=name.replace('.png','.txt')
	
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



'''
FUNCTION TO COMPUTE ALIGNMENT TIMES
'''
def doit(corg,options,originaldir):
	c=os.getcwd()
	
	f=os.listdir(c)
	
	#mults = [12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,40,42,44,45,48,49,50,52,54,56,60,64,65,66,70,72,75,77,78,80,81,84,88,91,96,98,100,104,112,120,128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,256]
	mults = []
	
	if options.boxsizes:
		#print "Boxsizes is", options.boxsizes
		#print "of type", type(options.boxsizes)
	
		if '-' in options.boxsizes and '.txt' not in options.boxsizes:
			lowerend = int( options.boxsizes.split('-')[0] )
			upperend = int( options.boxsizes.split('-')[1] ) 
			mults=[i for i in xrange( lowerend , upperend ) ]
			print "Lower and upper ends of range are", lowerend, upperend
			print "Therefore, boxsizes are", mults
	
		elif '.txt' in options.boxsizes:
			f=open(options.boxsizes,'r')
			lines=f.readlines()
			for line in lines:
				ele=line.replace('\n','').replace(' ','').split(',')[0]
				if ele:
					mults.append( int( float(ele) ) )
			f.close()
			
			print "Boxsizes read from file are", mults
	
		elif "," in options.boxsizes:			
			boxes=options.boxsizes.split(',')
			for box in boxes:
				mults.append( int(box) )
			
	elif options.test:
		mults=[100,101,120,121,140,141]
	
	elif options.medium:
		mults = [16,32,64,96,97,98,128,129,130,180,181,182]
		#steps = [30]

	elif options.extensive:
		mults = []
		for i in xrange(12,257):
			mults.append(i)
		#steps=[10]
	
	computer = options.ID
	if computer:
		computer += '_'
	
	data = {}
	
	coarsestep=options.coarsestep
	finestep=options.finestep
		
	#name = options.path + '/CS' + str(coarsestep).zfill(len(str(max(steps)))) + '_FS' + str(finestep) + '.txt'
	#if computer:
	
	
	IDS = []
	
	if not options.oneicosoff:		
		IDS.append('oneicos')

	if not options.onefulloff:
		IDS.append('onefull')
	
	if not options.oneccfoff:
		IDS.append('oneccf')
		
	if not options.rotonlyoff:
		IDS.append('rotonly')
	
	
	#print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nThe IDS and corg are", IDS, corg
	#print "The len of IDs is", len(IDS)
	#print "Individually, they are"
	#for aide in IDS:
	#	print "aide, corg are", aide,corg	
	#print "\n&&&&&&&&&&&&&&&&&&&&&&&&\n"

	for aidee in IDS:
		#print "\n\n\n\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nWorking with THIS id and corg", aidee,corg
		#print "\n&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n\n"
		name=options.path + '/' + computer + corg + '.txt'
		if aidee == 'oneicos':	
			name = options.path + '/' + computer + 'oneicos_CS' + str(coarsestep).zfill(len(str(coarsestep))) + '_FS' + str(finestep).zfill(len(str(finestep))) + '_' + corg + '.txt'
		if aidee == 'onefull':
			name = options.path + '/' + computer + 'onefull_' + corg + '.txt'
		if aidee == 'oneccf':
			name = options.path + '/' + computer + 'oneccf_' + corg + '.txt'
		if aidee == 'rotonly':
			name = options.path + '/' + computer + 'rotonly_CS' + str(coarsestep).zfill(len(str(coarsestep))) + '_FS' + str(finestep).zfill(len(str(finestep))) + '_'+ corg + '.txt'
		
		mastername = name.replace('.txt','_master.txt')
		
		times=[]
		cmd=''
		for size in mults:
			print "\naligning size", size
			print "for ID", aidee
			print "alignment type", corg
			print "\n"			
			#t=t1=t2=t1h=t2h=t1m=t2m=t1s=t2s=t1tot=t2tot=0
			
			setcuda=''
			if corg=='gpu' or corg=='GPU':
				setcuda= 'export NOCUDAINIT=0'
				#print "\n\n\n !!!! I Have turned cuda ON!!!\n\n\n"
				if options.setcudadevice is not None:
					setcuda += '&& export SETCUDADEVICE=' + options.setcudadevice
			elif corg=='cpu' or corg=='CPU':
				setcuda = 'export NOCUDAINIT=1'
				#print "\n\n\n !!!! I Have turned cuda OFF!!!\n\n\n"
			else:
				print "Something is wrong; you're supposed to select GPU or CPU. TERMINATING!"
				sys.exit()
			
			print "SETCUDA IS", setcuda
			cmd=''
			abspath = originaldir + '/' + options.path
			#print "\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@Abs path is", abspath
			os.system('cd ' + abspath)
			#print "But the dir I'm in is" + os.getcwd() + "@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n"
			
			a=EMData(size,size,size)
			a=a.process('testimage.noise.gauss')
			
			aname = abspath+ '/' + 'a_stack.hdf'
			a.write_image(aname)

			b=EMData(size,size,size)
			b=b.process('testimage.noise.gauss')
			bname = abspath + '/' +'b_stack.hdf'
			b.write_image(bname)

			out = 'garbage.hdf'
			
			#print "aname is", aname
			#print "bname is", bname

			#paf = abspath + '/' + aidee
			
				#cmd = setcuda + ''' && e2spt_classaverage.py --input=''' + aname + ''' --output=''' + out + ''' --ref=''' + bname + ''' --iter=1 --npeakstorefine=1 -v 0 --mask=mask.sharp:outer_radius=-2 --preprocess=filter.lowpass.gauss:cutoff_freq=0.1:apix=1.0 --align=rotate_translate_3d:search=6:delta=''' + str(coarsestep) + ''':dphi=''' + str(coarsestep) + ''':verbose=0:sym=icos --parallel=thread:1 --ralign=refine_3d_grid:delta=''' + str(finestep) + ''':range=''' + str( int(math.ceil(coarsestep/2.0)) ) + ''' --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --normproc=normalize.mask --path=''' + paf
				#cmds.update('icosali':cmd1)
		
			ta=0
			tb=0
			sym='c1'
			if aidee == 'onefull' or aidee == 'oneicos':	
				if aidee == 'oneicos':
					sym='icos'
				
				parallel=options.parallel
				profilecmd1='e2spt_classaverage.py'
				profilecmd2=''
				if options.profile:
					profilecmd1 = "python -m cProfile -s time " + options.eman2dir + "/bin/e2spt_classaverage.py"
					profilecmd2 = "> profiled_" + corg + '_box' + str(size).zfill(3) + '.txt' 
					parallel=''
					print "profilecmd1 is", profilecmd1
					print "profilecmd2 is", profilecmd2
									
				cmd = setcuda + ''' && cd ''' + abspath + ''' && ''' + profilecmd1 + ''' --input=''' + aname + ''' --output=''' + out + ''' --ref=''' + bname + ''' --iter=1 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=0.1:apix=1.0 --highpass=filter.highpass.gauss:cutoff_freq=0.01:apix=1.0 --preprocess=filter.lowpass.gauss:cutoff_freq=0.2:apix=1.0 --align=rotate_symmetry_3d:sym=''' + sym + ''' --parallel=''' + parallel + ''' --ralign=None --averager=mean.tomo --aligncmp=ccc.tomo --normproc=normalize.mask --path=''' + aidee + ''' --verbose=''' + str(options.verbose) + ' ' + profilecmd2
			
				print "Therefore the instruction is", cmd
				ta = time()
				os.system(cmd)
				tb = time()
	
			elif aidee == 'oneccf':
				if corg=='gpu' or corg=='GPU':
					EMData.switchoncuda()
					
				a=EMData(aname,0)
				b=EMData(bname,0)
				ta = time()
				ccfab=a.calc_ccf(b)
				tb = time()
				
				if corg=='gpu' or corg=='GPU':
					EMData.switchoffcuda()
				
			elif aidee == 'rotonly':
				if corg=='gpu' or corg=='GPU':
					EMData.switchoncuda()
					
				a=EMData(aname,0)
				b=EMData(bname,0)
				
				#--align=rotate_translate_3d:search=12:delta=12:dphi=12:verbose=1 --parallel=thread:8 --ralign=refine_3d_grid:delta=3:range=12:search=2
				ta = time()
				bestcoarse = a.xform_align_nbest('rotate_translate_3d',b,{'delta':int(options.coarsestep),'dphi':int(options.coarsestep),'search':10,'sym':options.sym},1,'ccc.tomo')
				for bc in bestcoarse:
					#classoptions["ralign"][1]["xform.align3d"] = bc["xform.align3d"]
					ran=int(round(float(options.coarsestep)/2.0))
					a.align('refine_3d_grid',b,{'delta':int(options.finestep),'range':ran,'search':3,'xform.align3d':bc['xform.align3d']},'ccc.tomo')
				tb = time()

				if corg=='gpu' or corg=='GPU':
					EMData.switchoffcuda()
			td = tb - ta
			if td:
				#print "Excution time was", td
				times.append(float(td))
				line2write= str(size) + ' ' + str( float(td) )+'\n'
				#txt.write(line2write)
				txt = open(mastername,'a')
				txt.write( line2write )
				txt.close()
					
		#txt.close()
	
		data.update({aidee:[mults,times]})
		#print "\n\nThe data to return is\n", data
		#print "\n"
	return(data)


'''
FUNCTION TO PLOT RESULTS
'''
def plotter(xaxis,yaxis,name='',CS=0,FS=0,markernum=0,linenum=0,ylimvalmax=0,idee='',yminnonconvex=[]):

	#yaxislog=[]
	if options.logplot:
		
		yaxis = [ math.log10(y) for y in yaxis ]	
		yminnonconvex = [ math.log10(ym) for ym in yminnonconvex ]
		
	print "\n\nreceived xaxis len is", len(xaxis)
	print "\n\nreceived yminnonconvex len is", len(yminnonconvex)
	
	if len(xaxis) != len(yminnonconvex):
		#print "ERROR! xaxis and yminnonconvex are not equal in length, see"
		#print "xaxis is", xaxis
		#print "yminnonconvex is", yminnonconvex
		#sys.exit()
		pass
	#print "in plotter, linen received is", linenum
	
	if not xaxis or not yaxis:
		#print "The plotter has received empty xaxis or yaxis arrays, see", 
		sys.exit()

	#print "IN PLOTTER X", xaxis
	#print "IN PLOTTER Y", yaxis
	#print "ylimvalmax RECEIVED in plotter is", ylimvalmax
	
	if ylimvalmax and not options.logplot:
		ymax=ylimvalmax
	
	else:
		ylimvalmax=max(yaxis)
		#print "YLIMVALMAX had to be reset!", ylimvalmax
		ymax=ylimvalmax	

	markers=['','o','D','H','x','h','8',7,4,5,6,'_','p',',','+','.','s','*','d',3,0,1,2,'1','3','4','2','v','<','>','^','|'] 
	mark=markers[markernum]
	
	linestyles=['-','--','**']
	linest=linestyles[linenum]
	
	
	
	
	
	


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

	#if options.logplot:
	#	plt.yscale('log')


	labelfory='Time (s)'
	
	if name and CS and FS:
		tag='gpu speed gain factor'
		labelfory='CPU time / GPU time'
		if 'gpu' in name or 'GPU' in name and 'cpu' not in name and 'CPU' not in name:
			tag='gpu 3D alignment Time'
			labelfory='Time (s)'
		if 'cpu' in name or 'CPU' in name and 'gpu' not in name and 'GPU' not in name:
			tag='cpu 3D alignment Time'
			labelfory='Time (s)'
		if options.logplot:
			labelfory = "LOG( " + labelfory + " )"
		stepslabel='\ncoarse step=' + str(CS) + ' : fine step=' + str(FS)
		plt.title(tag + ' VS box-size' + stepslabel)
	
	plt.ylabel(labelfory)
	plt.xlabel("Box side-length (pixels)")

	#plt.legend(['y = x', 'y = 2x', 'y = 3x', 'y = 4x'], loc='upper left')
	
	#a = plt.gca()
	#a.set_xlim(1,int(xaxis[-1]))
	#a.set_ylim(0,max(yaxis)+0.25*max(xaxis))
	#if idee:
	#	print "legend is", idee
	#	plt.legend(idee,loc='upper left')
	
	#plt.xlim( (1,int(xaxis[-1]) ) )
	#plt.ylim( ( 0,max(yaxis)+0.25*max(yaxis) ) )
	
	#plt.savefig(name)
	#plt.clf()
	return()


'''
FUNCTION TO DETERMINE MINIMA to be plotted later, instead of plotting ALL values in alignment time plots
'''
def minima(sizes,vals):
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
				
				#print "Plotting minima. Analyzing last three values, the len of all values is", len(vals)
				#print "i, i+1 and i+2 are", i, i+1, i+2
				#print "Len of sizes is", len(sizes)	
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
					if vals[i+1] < vals[i+2]:
						sizesmin.append(sizes[i+1])
						valsmin.append(vals[i+1])
					elif valis[i+1] > vals[1+2]:
						sizesmin.append(sizes[i+2])
						valsmin.append(vals[i+2])
								
			#print "I have appended this box, value", sizes[i], vals[i]

	minnonconvex=Util.nonconvex(valsmin,0)
	print "In minima, the yminnonconvex to return is", minnonconvex
	return(sizesmin,valsmin,minnonconvex)
	

if '__main__' == __name__:
	main()
