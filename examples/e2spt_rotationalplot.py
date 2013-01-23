#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012 - Last update 12/07/2012
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
from pylab import *
import pylab
import colorsys

pylab.rcParams['legend.loc'] = 'best'
#from operator import itemgetter					 

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Plots the variation of correlation of a volume with itself as it is rotated in azimuth or altitude"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vols1", type=str, help="Comma-separated filenames of the .hdf volumes whose self-rotational correlation plot you want to compute.", default=None)
	parser.add_argument("--vols2", type=str, help="Comma-separated filenames of the .hdf volumes whose rotational correlation plot you want to compute against the volumes provided through --vols1.", default=None)

	parser.add_argument("--output", type=str, help="Name for the .txt file with the results and the corresponding .png plot")
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--daz", type=int,default=3,help="Step size to vary azimuth.")
	parser.add_argument("--icosvertices", action="store_true",help="Will produce an azimutal plot at each vertex of an icosahedron.", default=False)
	parser.add_argument("--dalt", type=int,default=181,help="Step size to vary altitude.")
	parser.add_argument("--alti", type=int,default=0,help="""Initial position to check in altitude. For example, for a D symmetric chaperonin, 
															if you want to check alt=0 ONLY, provide --alti=0 and --dalt=181 as options.
															if you want to check alt=180 ONLY, provide --alti=180, --dalt=1 or greater.
															if you want to check BOTH alt=0 and alt=180 in the same plot, provide --alti=0, --dalt=180""")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")
	parser.add_argument("--plotonly",action="store_true", help="Assumes vol1 and vol2 are already aligned with respect to each other and thus skips alignment.", default=False)
	parser.add_argument("--normalizeplot",action="store_true", help="Make maximum correlation value on plot equal to 1 and scale all other values accordingly.", default=False)
	parser.add_argument("--singleplot",action="store_true", help="Plot all alts, or each vertex of --icosvertices is on, in a single .png file.", default=False)
	parser.add_argument("--only2dplot",action="store_true", help="Skips all plots, except 2dplot.", default=False)
	parser.add_argument("--savetxt",action="store_true", help="Will save the values for each plot into .txt files.", default=False)

	
	(options, args) = parser.parse_args()
		
	logger = E2init(sys.argv, options.ppid)
	
	if options.icosvertices and options.vols2:
			print "ERROR: You can only use --icosvertices for volumes in --vols1. You must NOT supply --vols2."
			sys.exit()
	

	vols1=[]
	if options.vols1:
		vols1 = options.vols1.split(',')
	
	vols = vols1

	vols2=[]
	if options.vols2:
		vols2 = options.vols2.split(',')
		vols = vols1 + vols2
		
		for v in vols:
			if '.hdf' not in v and '.mrc' not in v:
				print "ERROR: The input volumes must all be either .hdf, .mrc or .rec (which is also just a .mrc file)."
				sys.exit()
					
			rotcccplot(v,v,options)
	
		for v1 in vols1:
			for v2 in vols2:
				rotcccplot(v1,v2,options)
				
	else:
		for v in vols:
			if '.hdf' not in v and '.mrc' not in v:
				print "ERROR: The input volumes must all be either .hdf, .mrc or .rec (which is also just a .mrc file)."
				sys.exit()
	
			rotcccplot(v,v,options)
	
	E2end(logger)
			
	return()
			

def rotcccplot(v1,v2,options):	
	
	nimg1 = EMUtil.get_image_count()
	print "The first file actually has this many images in it", nimg1
	nimg2 = nimg1
	
	
	
	vol2 = EMData()
	if v1 != v2:
		nimg2 = EMUtil.get_image_count()
		if nimg2 > 1:
			
			print "The second file actually has this many images in it", nimg2
			print """Error: You cannot provide a stack for both --vol1 and --vol2.
				 To compare multiple volumes in an all vs all fashion, supply individual files separated by a comma.
				 You can compare multiple files in a stack against a single volume, if you provide the stack through --vol1,
				 and the single volume through --vol2."""
			sys.exit()


	for ni in range(nimg1):
		vol1 = EMData(v1,ni)
		vol1 = preprocess(vol1,options)
		title = v1
		
		if v1 != v2:
			vol2 = EMData(v2,0)
			vol2 = preprocess(vol2,options)
			if nimg1 > 1:
				title = v1 + '_ptcl' + str(ni).zfill(len(str(nimg1)))+'_VS_' + v2
			else:
				title = v1 + ' VS ' + v2
				
		else:
			vol2 = vol1.copy()
	
		loops = 1
		ts=[]
	
		if options.icosvertices:
			ts = genicosvertices()
			loops = len(ts)
	
		for loop in range(loops):
			print "looping over the stack %s; img number %d" %(v1,ni) 
			if ts:
				vol2 = vol1.copy()
				vol2.transform(ts[loop])
				print "I have gone to vertex %d by applying this transform" %(loop)
				print ts[loop]
						
			ret = azimuthalccc(vol1,vol2,options)
			azs = ret[0]
			values = ret[1]
		
			print "I have acquired azs and values and will proceed to plot."
	
			if options.normalizeplot:
				for ele in values:
					#val = values[ele]
				
					minv = min(values[ele])
					maxv1 = max(values[ele])
					#print "Min max before normalization was", minv,maxv1
					for k in range(len(values[ele])):
						values[ele][k] = values[ele][k] - minv 
				
					minv2 = min(values[ele])
					maxv = max(values[ele])
					#print "After subtracting min, the are", minv2,maxv
					#print "Max before normalization was", maxv
					for k in range(len(values[ele])):
						values[ele][k] = values[ele][k] / maxv
				
					#print "after norm they are", min(values[ele]), max(values[ele])
		
			#print "Len of azs and values is", len(azs), len(values)
			#print "therefore the az,values to plot are"
		
			vertextag=''
		
			for ele in values:
				plotname = title.replace('.hdf','').replace('.mrc','') + '.png'
				#alt = ele
				#print "There's a plot for this alt", ele
				#print "BASE plot name will be", plotname
				#print "because title is", title
			
				lines=[]
		
				vertextag = '_alt' + str(ele).zfill(3)
				if ts:
					vertextag = vertextag + '_' + str(loop).zfill(2)
		
				print "Vertextag is", vertextag
			
				txtname = plotname.replace('.png',vertextag + '.txt')
				if options.savetxt:
			
					print "txtname is", txtname
				
					f = open(txtname,'w')
					for i in range(len(values[ele])):
						line = str(azs[i]) + ' ' + str(values[ele][i]) + '\n'
						lines.append(line)
						#print "Line to write is", line
					f.writelines(lines)
					f.close()

				lab = 'alt='+str(ele).zfill(3)
				if options.icosvertices:
					lab = 'vertex' + str(loop).zfill(2)
				if not options.only2dplot:
					pylab.plot(azs, values[ele], linewidth=2, label=lab)
					legend()
			
				if options.normalizeplot:
					print "Normalize plot is on"
					ylim([0,1.5])
					if options.icosvertices:
						print "icosvertices is on"
						if options.singleplot:
							print "Single plot is on"
							ylim([0,6])
						else:
							print "Single plot is NOT on"
							ylim([0,1.5])
						
				else:
					ylim([0,1.5*max(values[ele])])
			
				if not options.singleplot:
					plotname = txtname.replace('.txt','.png')
					pylab.title(title)
					pylab.ylabel('Correlation')
					pylab.xlabel('Azimuth')
					pylab.savefig(plotname)
					clf()
		
		if options.singleplot:	
			if not options.only2dplot:
				pylab.savefig(plotname)	
				pylab.title(title)
				pylab.ylabel('Correlation')
				pylab.xlabel('Azimuth')	
				pylab.savefig(plotname)
				clf()
		
			if not options.icosvertices:
				print "I will call 2dplot"
				twoD_plot(plotname,values,options)
	return()
	
	
def preprocess(vol,options):
	print "Entering preprocessing function"
	if options.mask: 
		options.mask = parsemodopt(options.mask)
		print "mask is", options.mask
	
	if options.preprocess: 
		options.preprocess = parsemodopt(options.preprocess)
		print "Preprocessor is", options.preprocess
	
	if options.lowpass: 
		options.lowpass = parsemodopt(options.lowpass)
		print "lowpass is", options.lowpass
		
	if options.highpass: 
		options.highpass = parsemodopt(options.highpass)
		print "Highpass is", options.highpass
		
	apix = vol['apix_x']
	boxsize = vol['nx']
	
	#vol1 = EMData(options.input,0)
	
	# Make the mask first, use it to normalize (optionally), then apply it 
	mask = EMData(vol["nx"],vol["ny"],vol["nz"])
	mask.to_one()
	if options.mask:
		#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options["mask"][0],options["mask"][1]) 
		mask.process_inplace(options.mask[0],options.mask[1])

	# normalize
	print "Normalizing and masking"
	vol.process_inplace('normalize')
	vol.mult(mask)
	vol.process_inplace('normalize')
	vol.mult(mask)

	# preprocess
	if options.preprocess != None:
		print "Preprocessing"
		vol.process_inplace(options.preprocess[0],options.preprocess[1])

	# lowpass
	if options.lowpass != None:
		print "Lowpassing"
		vol.process_inplace(options.lowpass[0],options.lowpass[1])

	# highpass
	if options.highpass != None:
		print "Highpassing"
		vol.process_inplace(options.highpass[0],options.highpass[1])

	# Shrinking both for initial alignment and reference
	if options.shrink!=None and options.shrink > 1 :
		print "Shrinking"
		vol=vol.process("math.meanshrink",{"n":options.shrink})
	print "Leaving preprocessing function"
	return(vol)
	

def azimuthalccc(vol1,vol2,options):
	print "In azimuthalccc function"
	alt=options.alti
	values = {}
	azs=[]
	
	if options.icosvertices:
		options.dalt = 181
	k=0	
	while alt <= 180:
		#if not options.dalt:
		#	alt = 180
		az=0
		valuesforthisalt=[]
		while az <= 360:
			moving = EMData()
			moving = vol2.copy()
			
			#print "I will rotate by this alt and az", alt, az
			moving.rotate(az,alt,0)
			#if alt == 180:
				#print "In theory, I have flipped the molecule now!!!!", alt, az
			
			ccf = vol1.calc_ccf(moving)
			ccf.process_inplace("xform.phaseorigin.tocorner") 
			ccf.process_inplace('normalize')
				
			#box = ccf.get_zsize()
			#r =  Region((box/2) - int(parameters['searchx']), (box/2) - int(parameters['searchy']), (box/2) - int(parameters['searchz']), 2 * int(parameters['searchx']) + 1, 2 * int(parameters['searchy']) + 1, 2 * int(parameters['searchz']) + 1) 
			#sub_ccf = ccf.get_clip(r)

			loc_sub = ccf.calc_max_location()

			#xbest = loc_sub[0]
			#ybest = loc_sub[1]
			#zbest = loc_sub[2]

			best_value = ccf.get_value_at(loc_sub[0],loc_sub[1],loc_sub[2])
			valuesforthisalt.append(best_value)
			print "Done with ccc=%f for this alt=%d,az=%d" % ( best_value, alt, az)
			if k == 0:
				azs.append(az)
			az += options.daz
		
		#azs.append(azlist)
		values.update({alt:valuesforthisalt})
		k+=1
		
		#if not options.dalt:
		#	alt = 179
		#	options.dalt=2
			
		alt = alt + options.dalt
		#print "ALT now is!!!", alt
	print "LEAVING azimuthalccc function"

	return(azs,values)

def color(value):
	color =	colorsys.hsv_to_rgb( float(value) / 180.0 / (1.1), 1, 1)
	return(color)


def twoD_plot(plotname,values,options):
	print "I am in 2d plot"
	
	azs = set([])
	alts = set([])
	cccs = []
	
	widthy = options.dalt
	widthx = options.daz
	print "The step in alt is", options.dalt
	print "The step in az is", options.daz
	widths = [widthx,widthy]
	#markerwidth = max(widths)
	markerwidth = sum(widths)/2
	print "marker width is the average of the two", markerwidth
	
	for ele in values:
		for k in range(len(values[ele])):
			ccc = values[ele][k]
			pylab.plot(k*options.daz,ele,marker='o',markersize=markerwidth*2,color=color((ccc*ccc)*180),alpha=0.5, markeredgecolor=None)
			cccs.append(ccc)
			azs.add(k*options.daz)
			alts.add(ele)				
	
	plotname = plotname.replace('.png','_2D.png')

	pylab.title(plotname.replace('.png',''))
	pylab.xlabel('az')
	pylab.ylabel('alt')
	
	pylab.xlim([-widthx,max(azs) + markerwidth*2])
	pylab.ylim([-widthy,max(alts)+ markerwidth*2])
	#for i in range(len(finalpoints)):
	#	plt.plot(*zip(*[finalpoints[i]]),marker='o',markersize=4,color=color(ang_errors[i]))	
	plt.savefig(plotname,bbox_inches=0)
	clf()
	return()

def genicosvertices():
	t=Transform()

	syms = []
	for i in range(60):
		syms.append(t.get_sym('icos',i))
	
	ts = []
	alts = set([])
	phis = set([])
	for s in syms:
		rot=s.get_rotation()
		az=rot['az']
		alt=rot['alt']
		phi=rot['phi']
	
		if az != 0.0 and phi != 0.0:
			alts.add(round(alt,3))
			phis.add(round(phi,3))

	angles = {}		
	for a in alts:
		ps =set()
		for s in syms:
			rot=s.get_rotation()
			#az=rot['az']
			alt=round(rot['alt'],3)
			phi=round(rot['phi'],3)
			if a == alt:
				ps.add(phi)
		angles.update({a:ps})

	#print "angles is", angles
	for n in angles:
		alt = n
		ps = angles[n]
		for phi in ps:	
			ts.append( Transform({'type':'eman','alt':alt,'phi':phi}) )

	ts.append(Transform())
	ts.append(Transform({'type':'eman','alt':180}))

	return(ts)

	
if __name__ == "__main__":
	main()
