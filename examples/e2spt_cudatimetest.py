#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 06/10/2012
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
from sys import argv
from time import time
		 
import matplotlib.pyplot as plt
import sys
import numpy		 
		 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Tests cpu vs gpu speed on subtomogram alignment."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--cpu", action='store_true', help="Will test SPT alignment using CPU.",default=False)
	parser.add_argument("--gpu", action='store_true', help="Will test SPT alignment using GPU.",default=False)
	
	parser.add_argument("--test", action='store_true', help="Will run quick tests just to see if EMAN2 and this script are functional on your computer.",default=False)
	parser.add_argument("--ID", type=str, help="Tag files generated on a particular computer.",default='')
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsimjob'; for example, sptsimjob_02 will be the directory by default if 'sptsimjob_01' already exists.")
	
	(options, args) = parser.parse_args()
	
	print "options are", options
	
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
		options.path = "sptsim_01"
	
	files=os.listdir(os.getcwd())
	print "right before while loop"
	while options.path in files:
		print "in while loop, options.path is", options.path
		#path = options.path
		if '_' not in options.path:
			print "I will add the number"
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
			print "The new options.path is", options.path

	if options.path not in files:
		
		print "I will make the path", options.path
		os.system('mkdir ' + options.path)
	
	
	retcpu=[]
	retgpu=[]
	if options.cpu:
		corg = 'cpu'
		print "I will test the CPU"
		print "Because options.cpu is", options.cpu
		retcpu=doit(corg,options)
		print "\n\nThe received data for retCPU is\n", retcpu
		print "\n"
		print "Thus the length is", len(retcpu)
		
		for i in range(len(retcpu)):
			step = retcpu[i][0]
			name = options.path + "/CS"+str(step).zfill(2) + '_' + 'CPU.png'
			cnums = numpy.array(retcpu[i][-1])
			sizes = numpy.array(retcpu[i][1])
			print "\n$$$$$$$\nThe step is", step
			print "\n\n"
			plotter(name,sizes,cnums,step,step/2)
		
	if options.gpu:
		corg = 'gpu'
		print "I will test the GPU"
		print "Because options.gpu is", options.gpu
		retgpu=doit(corg,options)
		print "\n\nThe received data is\n", retgpu
		print "\n"
		
		for i in range(len(retgpu)):
			step = retgpu[i][0]
			name = options.path + "/CS"+str(step).zfill(2) + '_' + 'GPU.png'
			gnums = numpy.array(retgpu[i][-1])
			sizes = numpy.array(retgpu[i][1])
			print "\n$$$$$$$\nThe step is", step
			print "\n\n"
			plotter(name,sizes,gnums,step,step/2)
	
	if retcpu and retgpu:
		if len(retcpu) == len(retgpu):

			#steps=[]
			difs=[]
			for i in range(len(retgpu)):
				step = retgpu[i][0]
				name = options.path + "/CS"+str(step).zfill(2) + '_' + 'gpuVScpu.png'
				gnums = numpy.array(retgpu[i][-1])
				cnums = numpy.array(retcpu[i][-1])
				sizes = numpy.array(retgpu[i][1])
				difs = cnums/gnums
				print "\n$$$$$$$\nThe step is", step
				print "\n\n"
				plotter(name,sizes,difs,step,step/2)
			#print "I should be plotting this"
	
		else:
			print "For some sick reason, you don't have the same number of data points for gpu and cpu, therefore, you cannot compare them, see", len(retgpu), len(retcpu)
			sys.exit()
	else:
		return()
	return()

def doit(corg,options):
	c=os.getcwd()
	f=os.listdir(c)

	mults = [12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,40,42,44,45,48,49,50,52,54,56,60,64,65,66,70,72,75,77,78,80,81,84,88,91,96,98,100,104,112,120,128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,256]
	steps = [10,8,6,4,2]
	if options.test:
		mults = [12,14,15,22,32]
		steps = [6]

	#,81,84,88,91,96,98,100]
	
	computer = options.ID
	
	data = []
	
	for step in steps:

		coarsestep=step
		finestep=coarsestep/2
		
		name = options.path + '/' + corg + '_' + 'CS' + str(coarsestep).zfill(len(str(max(steps)))) + '_FS' + str(finestep) + '.txt'
		if computer:
			name = corg + '_' + computer + '_CS' + str(coarsestep).zfill(len(str(max(steps)))) + '_FS' + str(finestep) + '.txt'

		txt = open(name,'w')
		
		times=[]
		for size in mults:
			t=t1=t2=t1h=t2h=t1m=t2m=t1s=t2s=t1tot=t2tot=0

			#a=EMData(argv[1])
			a=EMData(size,size,size)
			a=a.process('testimage.noise.gauss')
			#a=a.process("xform.scale",{'scale':1,'clip':size})
			
			#if 'gpu' not in corg:
			#	a.switchoffcuda()
	
			#a.to_one()
			aname = 'a_stack.hdf'
			a.write_image(options.path+ '/' + aname)

			b=EMData(size,size,size)
			b=b.process('testimage.noise.gauss')
			bname = 'b_stack.hdf'
			b.write_image(options.path+ '/' + bname)
			
			#if 'gpu' not in corg:
			#	b.switchoffcuda()
			#	print "\n\n\n !!!! I Have turned cuda OFF!!!\n\n\n"
	
			out = 'garbage.hdf'
			
			setcuda=''
			if corg=='gpu':
				setcuda= 'EMANUSECUDA=1 && '
				print "\n\n\n !!!! I Have turned cuda ON!!!\n\n\n"
			else:
				setcuda = 'export EMANUSECUDA=0 && '
				print "\n\n\n !!!! I Have turned cuda OFF!!!\n\n\n"
			
			instruction1 = setcuda + ''' && cd ''' + options.path + '''e2spt_classaverage.py --input=''' + aname + ''' --output=''' + out + ''' --ref=''' + bname + ''' --iter=1 --npeakstorefine=1 -v 1 --mask=mask.sharp:outer_radius=-2 --preprocess=filter.lowpass.gauss:cutoff_freq=0.1:apix=1.0 --align=rotate_translate_3d:search=6:delta=''' + str(coarsestep) + ''':dphi=''' + str(coarsestep) + ''':verbose=1:sym=icos --parallel=thread:1 --ralign=refine_3d_grid:delta=''' + str(finestep) + ''':range=''' + str(coarsestep) + ''' --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --normproc=normalize.mask'''
			cmd = instruction1
			print "The instruction is", cmd
	
			#if 'gpu' in corg:
			#	a=test_image()
			#	a.do_fft_cuda()
			
			ta = time()
			os.system(cmd)
			tb = time()
			
			#t=eman2time()
			#print "The total alignment time was", t
			
			td = tb - ta
			print "BUt the real time is", td
			times.append(float(td))
			line2write= str(size) + ' ' + str(td)+'\n'
			txt.write(line2write)
		txt.close()
	
		data.append([step,mults,times])
		print "\n\nThe data to return is\n", data
		print "\n"
	return(data)


def plotter(name,xaxis,yaxis,CS,FS):
	tag='gpu speed gain factor'
	labelfory='CPU time / GPU time'
	if 'gpu' in name and 'cpu' not in name:
		tag='gpu 3D alignment Time'
		labelfory='Time (s)'
	if 'cpu' in name and 'gpu' not in name:
		tag='cpu 3D alignment Time'
		labelfory='Time (s)'
	
	#CS=name.split('CS')[1].split('_')[0]
	#FS=name.split('FS')[1].split('.')[0]
		
	#plot_name = name.replace('.txt','.png')
	stepslabel='\ncoarse step=' + str(CS) + ' : fine step=' + str(FS)

	plt.plot(xaxis, yaxis, linewidth=1)
	plt.title(tag + ' VS box-size' + stepslabel)
	

	#if 'dif' in name:

	plt.ylabel(labelfory)
	plt.xlabel("Box side-length (pixels)")
		
	a = plt.gca()
	a.set_xlim(1,int(xaxis[-1]))
	a.set_ylim(0,max(yaxis)+0.25*max(xaxis))
	#a.legend(stepslabel)

	plt.savefig(name)
	plt.clf()
	return()

'''
if what != 'dif' and '.' not in what:
	for i in fs:
		if what in i and '.txt' in i:
			testos.append(i)
	
	for i in testos:
	       g=open(i)
	       name=i
	       lines=g.readlines()
	       g.close()

	       mults=[]
	       x=[]
	       for line in lines:
		       m=line.split()[0]
		       mults.append(m)
		       xnum=line.split()[1]
		       x.append(float(xnum))
	       plotthis(name,mults,x)
	
elif what == 'dif':
	for i in fs:
		if 'gpu' in i and '.txt' in i:
			gpus.append(i)
			gpus.sort()
			print "I have found a gpu"
			
		if 'cpu' in i and '.txt' in i:
			cpus.append(i)
			cpus.sort()
			print "I have found a cpu"

	if len(gpus) != len(cpus):
		print "ERROR! G and C are not the same length!"
		print "Gs are", gpus
		print "Cs are", cpus
		sys.exit()

	else:
		for i in range(len(gpus)):
			gpuname=gpus[i]
			cpuname=cpus[i]

			if gpuname.replace('gpu','cpu') != cpuname:
				print "WARNING! You are not comparing the same files, see"
				print "GPU NAME IS", gpuname
				print "CPU name IS", cpuname
			if gpuname.split('_C')[1] == cpuname.split('_C')[1]:
				print "But they're the same step sizes..."
				finalname=gpuname.replace('gpu','dif')
				mults=[]
				x=[]
				print "gpus are", gpus
				gf=open(gpus[i],'r')
				glines=gf.readlines()
				print "glines are", glines

				print "cpus are", cpus
				cf=open(cpus[i],'r')
				clines=cf.readlines()
				print "clines are", clines
				for j in range(len(glines)):

					mg=glines[j].split()[0]
					mults.append(mg)
					xgnum=glines[j].split()[1]
			
					xcnum=clines[j].split()[1]
					x.append(float(xcnum)/float(xgnum))
				plotthis(finalname,mults,x)


else:	
	files=argv[2:-1]
	normalization=argv[1]
	plot_name=argv[-1]	
	print "The files to plot are", files
	print "And the name for the final plot is", plot_name
	print "And the file used for normalization is", normalization
	
	plt.title(plot_name.split('.')[0])
	#plots=[]
	a = plt.gca()
	mults=[]
	maxesx=[]
	maxesy=[]
	
	norm=open(normalization,'r')
	linesn=norm.readlines()
	norm.close()
	
	normnums=[]
	#normnums2=[]
	for line in linesn:
		normnum=line.split()[1]
		normnums.append(float(normnum))
		#normnums2.append(float(normnum/6.0))
		#normnums.append(1.0)
	
	print "There are these many files", len(files)
	print files
	for i in range(len(files)):
		#print "plotting this file", i
		g=open(files[i],'r')
	       	lines=g.readlines()
	       	g.close()
		
	       	mults=[]
	       	x=[]
	       	for j in range(len(lines)):
			m=lines[j].split()[0]
		       	mults.append(m)
		       	xnum=lines[j].split()[1]
		       	for j in len(sizes):
					gnums=retgpu[i][j][
			if i == 0:
				
				x.append(float(xnum)/normnums[j])
			
			if i>0:
				x.append(float(xnum)/(6.0*normnums[j]))
				
			#print "Values are", i, m, x
		
		if len(mults) >0:
			maxesx.append(mults[-1])				
		if len(x) > 0:		
			maxesy.append(max(x))	
		
		if len(mults) != len(x):
			print "The length of mults and x are not equal, see"
			print len(mults)
			print len(x)
			print mults
			print x
		
		plt.plot(mults, x, linewidth=1)		
		#plt.legend(pi,(i))

		#plots.append(pi)
			
	if len(maxesx)>0:	
		a.set_xlim(1,int(max(maxesx)))
	if len(maxesy)>0:
		a.set_ylim(0,int(max(maxesy)))
	

	plt.ylabel('cpu_time / gpu_time')
	plt.xlabel("Box side-length (pixels)")
		
	#plt.legend( (pi[0][0], pi[1][0]), ('Men', 'Women') )

	plt.savefig(plot_name)
	plt.clf()
'''










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
