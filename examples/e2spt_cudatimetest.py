#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 05/March/2013
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
	usage = """WARNING: Experimental script. 
	Tests cpu vs gpu speed on subtomogram alignment, for one CCF between two boxes of noise at default orientation (no rotations), 
	with and without alignment overhead, and plots the results against boxsize.
	Full alignment scanning an entire icosahedral unit can be tested. The number of orientations visited depends on --coarsestep and --finestep,
	and the way these parameters are used by e2spt_classaverage.py"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--cpu", action='store_true', help="Will test SPT alignment using CPU.",default=False)
	parser.add_argument("--gpu", action='store_true', help="Will test SPT alignment using GPU.",default=False)
	
	parser.add_argument("--test", action='store_true', help="Will run a quick tests using a few box sizes.",default=False)
	parser.add_argument("--medium", action='store_true', help="Will test boxsizes in multiples of 10 between 10 and 240.",default=False)
	parser.add_argument("--extensive", action='store_true', help="Will test EVERY box size between 12 and 256.",default=False)
	
	parser.add_argument("--onefulloff", action='store_true', help="Will test time for one CCF (between two boxes with noise) with all processing overhead for EVERY box size between 12 and 256.",default=False)
	parser.add_argument("--oneccfoff", action='store_true', help="Will test time for one CCF (between two boxes with noise)  without overhead, for EVERY box size between 12 and 256.",default=False)
	parser.add_argument("--oneicosoff", action='store_true', help="""Will test alignment time for however many orientations fit in one icosahedral unit, 
																	depending on --coarsestep and --finestep, for EVERY box size between 12 and 256.""",default=False)
	
	parser.add_argument("--ID", type=str, help="Tag files generated on a particular computer.",default='')
	parser.add_argument("--coarsestep", type=int, help="Step size for coarse alignment.",default=30)
	parser.add_argument("--finestep", type=int, help="Step size for fine alignment.",default=15)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsimjob'; for example, sptsimjob_02 will be the directory by default if 'sptsimjob_01' already exists.")
	parser.add_argument("--noplot", action='store_true', help="To run jobs on a cluster or other interfaces that don't support plotting, turn this option on.",default=False)
	
	parser.add_argument("--plotonly",type=str, help="""If you already have the time variation with boxsize (for a particular cpu or gpu) in a text file in rows of 'boxsize,time', 
												provide the txt file(s) separated by commas --plotonly=file1.txt,file2.txt,file3.txt etc...""", default=None)
	parser.add_argument("--singleplot",type=str, help="Provide the name of the .png file to plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default='')
	
	parser.add_argument("--plotminima",action="store_true", help="Plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default=False)
	parser.add_argument("--colorlessplot",action="store_true", help="Plot all cpus and/or gpus provided through --plotonly, in a single .png file.", default=False)
	
	parser.add_argument("--subset",type=int, help="If --plotonly is on, this is the subset of lines (from 0 to --subset=n) in the files provided that will be kept for plotting.", default=0)

	global options
	
	(options, args) = parser.parse_args()
	
	
	#print "options are", options
	
	logger = E2init(sys.argv, options.ppid)
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	#if options.path and ("/" in options.path or "#" in options.path) :
	#	print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
	#	sys.exit(1)
		
	#if options.path and options.path[:4].lower()!="bdb:": 
	#	options.path="bdb:"+options.path

	#if not options.path: 
	#	#options.path="bdb:"+numbered_path("sptavsa",True)
	#	options.path = "sptsim_01"
	
	
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)

	if not options.path: 
		#options.path="bdb:"+numbered_path("sptavsa",True)
		options.path = "sptCudaTest_01"
	
	files=os.listdir(os.getcwd())
	#print "right before while loop"
	while options.path in files:
		#print "in while loop, options.path is", options.path
		#path = options.path
		if '_' not in options.path:
			#print "I will add the number"
			options.path = options.path + '_00'
		else:
			jobtag=''
			components=options.path.split('_')
			if components[-1].isdigit():
				components[-1] = str(int(components[-1])+1).zfill(2)
			else:
				components.append('00')
						
			options.path = '_'.join(components)
			#options.path = path
			#print "The new options.path is", options.path

	if options.path not in files:
		
		#print "I will make the path", options.path
		os.system('mkdir ' + options.path)
	
	
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
				
				#name = options.path + "/" + options.ID + dset + "_CS"+ str(options.coarsestep).zfill(2)  + "_FS" + str( options.finestep ).zfill(2) + '_CPU.png'
				#sizes = data[0]				
				#cnums = data[1]
				
				#plotter(sizes,cnums,name,options.coarsestep,options.finestep)				
				#plt.savefig(options.path + '/' + os.path.basename(name))
				#plt.clf()
			
				#if options.plotminima:
				#	ret=minima(sizes,cnums)
				#	sizesmin=ret[0]
				#	cnumsmin=ret[1]
				#	namemin=name.replace('.png','_MIN.png')
				#	markernum=0
				#	if options.colorlessplot:
				#		markernum=1
				#	plotter(sizesmin,cnumsmin,namemin,options.coarsestep,options.finestep,0,0,markernum)
				#	plt.savefig(options.path + '/' + os.path.basename(namemin))
				#	plt.clf()
	

		'''
		I you have both CPU and GPU times, compute the ratio and plot it
		'''
		if retcpu and retgpu:
			if len(retcpu) == len(retgpu):

				#steps=[]
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
	
						#print "\n$$$$$$$\nThe step is", step
						#print "\n\n"
						#plotter(sizes,difs,name,options.coarsestep,options.finestep)
						#plt.savefig(options.path + '/' + os.path.basename(name))
						#plt.clf()
	
						#if options.plotminima:
						#	ret=minima(sizes,difs)
						#	sizesmin=ret[0]
						#	difsmin=ret[1]
						#	namemin=name.replace('.png','_MIN.png')
						#	markernum=0
						#	if options.colorlessplot:
						#		markernum=1
						#	plotter(sizesmin,difsmin,namemin,options.coarsestep,options.finestep,0,0,markernum)
						#	plt.savefig(options.path + '/' + os.path.basename(namemin))
						#	plt.clf()
					else:
						print "This gpu set does not have the same number of elements as the corresponding cpu set"
						
				#print "I should be plotting this"

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
		print "Will plot these files", files
		
		values={}
		
		k=0
		absmax=0.0
		absminmax=0.0
		markn=0
		linen=0
		for F in files:
			
			if options.colorlessplot:
				markn=k+1
				linen=k
			
			print "Working with this file now", F
			name=os.path.basename(F).replace('.txt','.png')
			sizes=[]
			valuesforthisfile=[]
			
			f=open(F,'r')
			lines = f.readlines()
			f.close()
			
			num=0
			for line in lines:
				#print "Line is\n", line
				size=line.split()[0]
				sizes.append(int(size))
				value=line.split()[-1].replace('\n','')
				valuesforthisfile.append(float(value))
				if options.subset:
					if int(num) >= int(options.subset):
						break
				num+=1
				
			if float(max(valuesforthisfile)) > absmax:
				absmax = float(max(valuesforthisfile))
				print "New ABSMAX is", absmax
			
			
			if options.plotminima:
				ret=minima(sizes,valuesforthisfile)
				sizesmin=ret[0]
				valuesmin=ret[1]
				
				if float(max(valuesmin)) > absminmax:
					absminmax = float(max(valuesmin))
					
					print "New ABSMINMAZ is", absmax
					
				namemin=name.replace('.png','_MIN.png')
				print "before plotting, linen is", linen
				plotter(sizesmin,valuesmin,namemin,0,0,markn,linen,absminmax)
				
				if options.singleplot:
					pass
				else:
					plt.savefig(options.path + '/' + os.path.basename(namemin))
					plt.clf()
					
			else:
				print "before plotting, linen is", linen
				plotter(sizes,valuesforthisfile,name,0,0,0,linen,absmax)
				
				if options.singleplot:
					pass
				else:
					plt.savefig(options.path + '/' + os.path.basename(name))
					plt.clf()
			
			#values.update({k:valuesforthisfile})
			k+=1
		
		if options.singleplot:

			
			plt.savefig(options.path + '/' + os.path.basename(options.singleplot))
			plt.clf()
			print "I have saved a single plot to", options.path + '/' + os.path.basename(options.singleplot)
			
	return()


def preplot(options,tipo,data,key):

	name = options.path + "/" + options.ID + key + "_CS"+ str(options.coarsestep).zfill(2)  + "_FS" + str( options.finestep ).zfill(2) + '_' + tipo + '.png'
	xdata = data[0]				
	ydata = data[1]				
	#print "\n\n"
	#print "BEFORE plotter sizes are", sizes
	#print "BEFORE plotter GNUMS to plot are", gnums
	#print '\n'
	plotter(xdata,ydata,name,options.coarsestep,options.finestep)
	plt.savefig(options.path + '/' + os.path.basename(name))
	plt.clf()
	
	textwriter(name,xdata,ydata)

	if options.plotminima:
		ret=minima(xdata,ydata)
		xdatamins=ret[0]
		ydatamins=ret[1]
		namemin=name.replace('.png','_MIN.png')
		markernum=0
		if options.colorlessplot:
			markernum=1
		plotter(xdatamins,ydatamins,namemin,options.coarsestep,options.finestep,0,0,markernum)
		plt.savefig(options.path + '/' + os.path.basename(namemin))
		plt.clf()

		textwriter(namemin,xdatamins,ydatamins)

	if options.subset:
		xdatasub=xdata[0:options.subset]
		ydatasub=ydata[0:options.subset]
		namesub=name.replace('.png','_SUB' + str(options.subset).zfill(3) + '.png')
		
		plotter(xdatasub,ydatasub,namesub,options.coarsestep,options.finestep)
		plt.savefig(options.path + '/' + os.path.basename(namesub))
		plt.clf()
	
		textwriter(namesub,xdatasub,ydatasub)

		if options.plotminima:
			ret=minima(xdatasub,ydatasub)
			xdatasubmins=ret[0]
			ydatasubmins=ret[1]
			namesubmin=namesub.replace('.png','_MIN.png')
			markernum=0
			if options.colorlessplot:
				markernum=1
			plotter(xdatasubmins,ydatasubmins,namesubmin,options.coarsestep,options.finestep,0,0,markernum)
			plt.savefig(options.path + '/' + os.path.basename(namemin))
			plt.clf()
	
			textwriter(namesubmin,xdatasubmins,ydatasubmins)

	return()


def textwriter(name,xdata,ydata):

	#difsl = list(difs)
	#for i in range(len(difsl)):
	#	difsl[i] = str(difsl[i]) + '\n'
	#f = open(name.replace('.png','.txt'),'w')
	#f.writelines(difsl)
	#f.close()
	
	filename=name.replace('.png','.txt')
	f=open(filename,'w')
	lines=[]
	for i in range(len(xdata)):
		line2write= str(xdata[i]) + ' ' + str(ydata[i])+'\n'
		lines.append(line2write)
	
	f.writelines(filename)
	f.close()

	return()



'''
FUNCTION TO COMPUTE ALIGNMENT TIMES
'''
def doit(corg,options,originaldir):
	c=os.getcwd()
	

	f=os.listdir(c)
	
	mults = [12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,40,42,44,45,48,49,50,52,54,56,60,64,65,66,70,72,75,77,78,80,81,84,88,91,96,98,100,104,112,120,128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,256]

	if options.test:
		mults = [32,64,96,128]
		#steps = [30]

	if options.medium:
		mults=[]
		for i in xrange(1,25):
			mults.append(i*10)
		#steps = [30]

	if options.extensive:
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
	
	
	#print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nThe IDS and corg are", IDS, corg
	#print "The len of IDs is", len(IDS)
	#print "Individually, they are"
	#for aide in IDS:
	#	print "aide, corg are", aide,corg	
	#print "\n&&&&&&&&&&&&&&&&&&&&&&&&\n"

	for aidee in IDS:
		#print "\n\n\n\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nWorking with THIS id and corg", aidee,corg
		#print "\n&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n\n"
		name=options.path + '/' + computer +'_' + corg + '.txt'
		if aidee == 'oneicos':	
			name = options.path + '/' + computer + 'oneicos_CS' + str(coarsestep).zfill(len(str(coarsestep))) + '_FS' + str(finestep) + '_' + corg + '.txt'
		if aidee == 'onefull':
			name = options.path + '/' + computer + 'onefull_CS' + str(coarsestep).zfill(len(str(coarsestep))) + '_' + corg + '.txt'
		if aidee == 'oneccf':
			name = options.path + '/' + computer + 'oneccf_CS' + str(coarsestep).zfill(len(str(coarsestep))) + '_' + corg + '.txt'
		
		#txt = open(name,'w')
		times=[]
		cmd=''
		for size in mults:
			#t=t1=t2=t1h=t2h=t1m=t2m=t1s=t2s=t1tot=t2tot=0
			
			setcuda=''
			if corg=='gpu' or corg=='GPU':
				setcuda= 'export NOCUDAINIT=0'
				#print "\n\n\n !!!! I Have turned cuda ON!!!\n\n\n"
			elif corg=='cpu' or corg=='CPU':
				setcuda = 'export NOCUDAINIT=1'
				#print "\n\n\n !!!! I Have turned cuda OFF!!!\n\n\n"
			else:
				print "Something is wrong; you're supposed to select GPU or CPU. TERMINATING!"
				sys.exit()
		
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
			
				cmd = setcuda + ''' && cd ''' + abspath + ''' && e2spt_classaverage.py --input=''' + aname + ''' --output=''' + out + ''' --ref=''' + bname + ''' --iter=1 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=0.1:apix=1.0 --highpass=filter.highpass.gauss:cutoff_freq=0.01:apix=1.0 --preprocess=filter.lowpass.gauss:cutoff_freq=0.2:apix=1.0 --align=rotate_symmetry_3d:sym=''' + sym + ''' --parallel=thread:1 --ralign=None --averager=mean.tomo --aligncmp=ccc.tomo --normproc=normalize.mask --path=''' + aidee
			
				#print "The instruction is", cmd
				ta = time()
				os.system(cmd)
				tb = time()
	
			elif aidee == 'oneccf':
				a=EMData(aname,0)
				b=EMData(bname,0)
				ta = time()
				ccfab=a.calc_ccf(b)
				tb = time()
			
			td = tb - ta
			if td:
				#print "Excution time was", td
				times.append(float(td))
				#line2write= str(size) + ' ' + str(td)+'\n'
				#txt.write(line2write)
				
		#txt.close()
	
		data.update({aidee:[mults,times]})
		#print "\n\nThe data to return is\n", data
		#print "\n"
	return(data)


'''
FUNCTION TO PLOT RESULTS
'''
def plotter(xaxis,yaxis,name='',CS=0,FS=0,markernum=0,linenum=0,ylimvalmax=0):
	#print "in plotter, linen received is", linenum
	
	if not xaxis or not yaxis:
		#print "The plotter has received empty xaxis or yaxis arrays, see", 
		sys.exit()

	#print "IN PLOTTER X", xaxis
	#print "IN PLOTTER Y", yaxis
	#print "ylimvalmax RECEIVED in plotter is", ylimvalmax
	
	if ylimvalmax:
		ymax=ylimvalmax
	
	else:
		ylimvalmax=max(yaxis)
		#print "YLIMVALMAX had to be reset!", ylimvalmax
		ymax=ylimvalmax	

	markers=['','o','D','H','x','h','8',7,4,5,6,'_','p',',','+','.','s','*','d',3,0,1,2,'1','3','4','2','v','<','>','^','|'] 
	mark=markers[markernum]
	
	linestyles=['-','--','**']
	linest=linestyles[linenum]
	
	#print "Therefore linest is", linest
	
	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	ax=plt.axes(frameon=False)
	pylab.rc("axes", linewidth=2.0)
	pylab.xlabel('X Axis', fontsize=14, fontweight='bold')
	pylab.ylabel('Y Axis', fontsize=14, fontweight='bold')

	pylab.ylim([-10,ylimvalmax+10])
	pylab.xlim([-10,max(xaxis)+10])
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.axes.get_xaxis().set_visible(True)
	ax.axes.get_yaxis().set_visible(True)
	xmin, xmax = ax.get_xaxis().get_view_interval()
	ymin, ymax = ax.get_yaxis().get_view_interval()
	
	ax.add_artist(Line2D((xmin, xmax+10), (ymin, ymin), color='k', linewidth=4))
	ax.add_artist(Line2D((xmin, xmin), (ymin, ymax+10), color='k', linewidth=4))
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)

	#print "BOLD IS ON!"
	LW=3
	if not markernum:
		LW=2
		
	if options.colorlessplot:
		print "in colorless plot, linest is", linest
		plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,color='k',zorder=0)
		#markeredgewidth=3
		if mark:
			plt.scatter(xaxis,yaxis,marker=mark,edgecolor='k',alpha=1,zorder=1,s=40,facecolor='white',linewidth=2)
	else:
		plt.plot(xaxis, yaxis, linewidth=LW,linestyle=linest,alpha=1,zorder=0)
		#markeredgewidth=3
		if mark:
			plt.scatter(xaxis,yaxis,marker=mark,alpha=0.5,zorder=1,s=40,linewidth=2)

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
	
		stepslabel='\ncoarse step=' + str(CS) + ' : fine step=' + str(FS)
		plt.title(tag + ' VS box-size' + stepslabel)
	
	plt.ylabel(labelfory)
	plt.xlabel("Box side-length (pixels)")

	#plt.legend(['y = x', 'y = 2x', 'y = 3x', 'y = 4x'], loc='upper left')
	
	#a = plt.gca()
	#a.set_xlim(1,int(xaxis[-1]))
	#a.set_ylim(0,max(yaxis)+0.25*max(xaxis))
	#a.legend(stepslabel)
	
	#plt.xlim( (1,int(xaxis[-1]) ) )
	#plt.ylim( ( 0,max(yaxis)+0.25*max(yaxis) ) )
	
	#plt.savefig(name)
	#plt.clf()
	return()


'''
FUNCTION TO DETERMINE MINIMA to be plotted later, instead of plotting ALL values in alignment time plots
'''
def minima(sizes,vals):
	
	sizesmin=[sizes[0]]
	valsmin=[vals[0]]
	
	for i in xrange(0,len(vals) - 1 - 1 ):
		aux = 0 
		for j in xrange(i+1,len(vals) - 1 - 1 ):
			if vals[j]< vals[i]:
				aux=0
				#print "Because a downstream value is lower, I will break the loop, see", vals[j],vals[i]
				break
			else:
				aux=1
		if aux==1:
			sizesmin.append(sizes[i])
			valsmin.append(vals[i])
			#print "I have appended this box, value", sizes[i], vals[i]
	
	return(sizesmin,valsmin)

'''
def eman2time():

		f=open('.eman2log.txt','r')
		lines=f.readlines()
		last=lines[-1]
	
		t1=last.split()[1]
		print "\n\n\n@@@@@@@@@@@@@@@@@@@@@ The job started at time", t1

		t1h=int(t1.split(':')[0].replace('/',''))
		#print "The hours of time 1 are", t1h
		#print "Or that times 3600", t1h*3600
		t1m=int(t1.split(':')[1])
		#print "The minutes of time 1 are", t1m
		#print "Or that times 60", t1m*60

		t1s=int(t1.split(':')[2])
		#print "The seconds of time 1 are", t1s
				
		t1tot = t1h*3600 + t1m*60 + t1s
		print "\n\n******************************* Therefore, in seconds, this was the time of starting t1tot=%d, %d*3600 + %d*60 + %d" % (t1tot,t1h,t1m,t1s)
				
		t2=last.split()[3]
		print "\n\n\n@@@@@@@@@@@@@@@@@@@@The job finished at time", t2
			
		t2h=int(t2.split(':')[0].replace('/',''))
		print 't2h is', t2h
			
		t2m=int(t2.split(':')[1])
		print 't2m is', t2m
			
		t2s=int(t2.split(':')[2])
		print 't2s is', t2s

		t2tot = t2h*3600 + t2m*60 + t2s
		print "\n\n******************************* Therefore, in seconds, this was the time of finishing t2tot=%d, %d*3600 + %d*60 + %d\n\n" % (t2tot,t2h,t2m,t2s)
	
		t = t2tot - t1tot
	return(t)
'''

if '__main__' == __name__:
	main()
