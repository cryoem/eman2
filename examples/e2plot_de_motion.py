#!/usr/bin/env python
'''
====================
Author: Jesus Galaz - 2015, Last update: 12/July/2017
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

from sys import argv
import matplotlib
matplotlib.use('Agg',warn=False)

import matplotlib.pyplot as plt
import pylab
import math
from EMAN2 import *

import matplotlib.colors as mcol
import matplotlib.cm as cm

import numpy as np

import colorsys

from operator import itemgetter


def main():

	progname = os.path.basename(sys.argv[0])
	usage = """plot motion from DE frames using .txt files from DE_process_frames.py"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--apix", type=float, default=0.0, help="""default=0.0. apix of data""")
	
	parser.add_argument("--bidirectionalfrom",type=str,default='',help="""Default=None (not used). Used for tiltseries data. Initial angle for the first half of the tiltseries. For example, a tiltseries from 0 to -55, then 5 to 55, --bidrectionalfrom should be set to 0.""")

	parser.add_argument("--highestangle", type=str, default=55, help="Default=55. Most positive tilt angle; e.g., 60 in a [-60,60] tiltseries")
	parser.add_argument("--highresolution", action='store_true', default=False, help="""default=False. If on, this option will enforce writing high-resolution plots (dpi 300), ready for publication, as opposed to lower resolution plots which take up less space on your computer (dpi 150).""")

	parser.add_argument("--individualplots", action='store_true', default=False, help="""default=False. in addition to plotting the motion for all frames in all images in a single plot, generate individual plots per image.""")
	
	parser.add_argument("--lowestangle", type=str, default=-55, help="Default=-55. Most negative tilt angle; e.g., -60 in a [-60,60] tiltseries")

	#parser.add_argument("--miny",type=float,default=0.0,help="""Default=None. Minimum value to plot in y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	#parser.add_argument("--maxy",type=float,default=0.0,help="""Default=None. Maximum value to plot in y. Automatically set to the maximum value in the data, per image, if not explicitly set.""")
	
	parser.add_argument("--negativetiltseries",action='store_true',default=False,help="""Default=False. This means the negative tilt angles in a bidrectional tiltseries were collected first; e.g., 0 to -55, then 5 to 55.""")
	
	parser.add_argument("--outputtag",type=str,default='deplot',help="""default=deplot. string common to all automatically generated output files""")

	parser.add_argument("--path", type=str,default='de_plots',help="""Default=de_plots. Name of the directory where to store the output results.""")
	parser.add_argument("--ppid", type=int, default=-1,help="Default=-1. Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--savetxts", action='store_true', default=False, help="""Default=False. save user-friendly txt files easyily plotable with other software.""")

	parser.add_argument("--tiltstep", type=int, default=0, help="Default=None. Angular step size between images in the tiltseries.")
	parser.add_argument("--tltfile", type=str,default='',help="""Default=None. Name of .tlt or .rawtlt file (typically from IMOD or serialEM) with the tilt angles.""")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	if not options.apix:
		print "\n(e2plot_de_motion)(main) ERROR: --apix required to accurately plot motion in Angstroms."
		sys.exit(1)

	anglesin = []
	if not options.tltfile:
		if not options.tiltstep or not options.lowestangle or not options.highestangle:
			print "\n(e2plot_de_motion)(main)ERROR: must provide --tiltstep, --lowestangle, and --highestangle in the absence of --tltfile"
			sys.exit(1)

		options.lowestangle = int(options.lowestangle)
		options.highestangle = int(options.highestangle)
		options.tiltstep = int(options.tiltstep)

	elif options.tltfile:
		if options.highestangle or options.lowestangle or options.tiltstep:
			print "\n(e2plot_de_motion)(main) WARNING: --tltfile was provided; therefore, --highestangle,--lowestangle, and --tiltstep will be ignored."
		
		with open( options.tltfile ,'r' ) as tltfile: 
			lines = tltfile.readlines()
			anglesin = [ int(round(float( line.replace('\n','') ))) for line in lines]
			if anglesin[0] < 0:
				options.lowestangle = anglesin[0]
				options.highestangle = anglesin[-1]
			elif anglesin[0] > 0:
				options.lowestangle = anglesin[-1]
				options.highestangle = anglesin[0]
			options.tiltstep = int(round( ( math.fabs(anglesin[0]) + math.fabs(anglesin[-1]) )/len(anglesin) ))

	from EMAN2_utils import makepath
	options = makepath(options)

	#print "apix is", options.apix
	#print "name is", options.outputtag	

	if options.path not in options.outputtag:
		options.outputtag = options.path + '/' + options.outputtag

	c = os.getcwd()
	
	listdir = os.listdir( c )
	
	ids = set([])
	

	for f in listdir:
		if options.verbose:
			print "\nfound file",f
		
		if '_y.txt' in f or '_x.txt' in f:
			stem = f.replace('y.txt','').replace('x.txt','')
			ids.add( stem )
			if options.verbose:
				print "\nfound stem", stem
				
	lowestangle = options.lowestangle
	highestangle = options.highestangle

	print "\nlowestangle {} type={}".format(lowestangle, type(lowestangle))
	print "\nhighestangle {} type={}".format(highestangle, type(highestangle))
	#print "\nstartangle {} type={}".format(startangle, type(startangle))		
	
	step = options.tiltstep
	
	datadict = {}
	
	lastangle = highestangle

	#if not options.tltfile:
	kk=0
	if options.bidirectionalfrom:
		startangle = int(options.bidirectionalfrom)
		print "\nstartangle {} type={}".format(startangle, type(startangle))
		if options.negativetiltseries:
			for angle in range(startangle,lowestangle-1,-step):
				datadict.update({kk:[angle]})
				lastangle=angle
				kk+=1

			for angle in range(startangle+step,highestangle+1,step):
				datadict.update({kk:[angle]})
				lastangle=angle
				kk+=1

		elif not options.negativetiltseries:
			for angle in range(startangle,highestangle+1,step):
				datadict.update({kk:[angle]})
				lastangle = angle
				kk+=1

			for angle in range(startangle-step,lowestangle-1,-step):
				datadict.update({kk:[angle]})
				lastangle=angle
				kk+=1

	elif not options.bidirectionalfrom:
		if options.negativetiltseries:
			startangle = highestangle
			print "\nstartangle {} type={}".format(startangle, type(startangle))
			for angle in range(startangle,lowestangle-1,-step):
				datadict.update({kk:[angle]})
				lastangle=angle
				kk+=1
		elif not options.negativetiltseries:
			startangle = lowestangle
			print "\nstartangle {} type={}".format(startangle, type(startangle))
			for angle in range(startangle,highestangle+1,step):
				datadict.update({kk:[angle]})
				lastangle = angle
				kk+=1

	#elif options.tltfile:

	if len(datadict) != len(ids):
		dif = len(datadict) - len(ids)
		print "\n(DE_translation_plotter)(main) WARNING: there are these many files {}, but only these many angles, {}, given input parameters (lowestangle={}, highestangle={}, and tiltstep={}). Therfore, {} additional angles will be added to assign a different angle to each file.".format(len(ids),len(datadict),lowestangle,highestangle,step,int(math.fabs(dif)))
		
		tmpangle=lastangle+step
		for ii in xrange( int(math.fabs(dif)) ):

			datadict.update({kk:[tmpangle]})
			kk+=1
			tmpangle+=step
	
	ids=list(ids)
	ids.sort()
	
	filesdict = {}

	kk=0

	for iid in ids:
		if options.verbose>9:
			print "\n(DE_translation_plotter)(main) or kk={} iid is {}".format(kk,iid)
			print "therefore datadict[kk]={}, type={}".format(datadict[kk],type(datadict[kk]))
		datadict[kk].append(iid)
		kk+=1

	if options.verbose > 9:
		print "\n(DE_translation_plotter)(main) ids before are {} type={} ".format( ids, type(ids) )
		print "\ndatadict is {}".format(datadict)

	datadictsorted = sorted(datadict.items(), key=lambda e: e[1][0]) 
	
	if options.verbose > 9:
		print "\n(DE_translation_plotter)(main) datadictsorted is {}".format(datadictsorted)
		print "\nTYPE {}".format(type(datadictsorted))
	
	
	figx_vals = {}
	figy_vals = {}
	figr_vals = {}

	avgslist = []
	
	xavgs = []
	xtotals = []
	xerrors = []

	yavgs = []
	ytotals = []
	yerrors = []

	ravgs = []
	rtotals = []
	rerrors = []

	angles = []

	nimgs = len(ids)
	
	k=0

	for ele in datadictsorted:
		
		id = ele[-1][-1]
		angletoplot = ele[-1][0]
		angles.append(angletoplot)
		if options.verbose > 9:
			print "\n igm is {} and angle is {}".format(id,angletoplot)

		if k >nimgs:
			break

		figx_values_angstroms = getvaluesfromfile(options, id, 'x')[0]
		figx_vals.update( {k:[id,angletoplot,figx_values_angstroms,'x']} )
		
		figx_values_angstroms_abs = [math.fabs(x) for x in figx_values_angstroms]
		#xavg = sum( [math.fabs(x) for x in figx_values_angstroms] )/len(figx_values_angstroms)
	
		xavg = np.mean(figx_values_angstroms_abs)
		xerror = np.std(figx_values_angstroms_abs)
		xtotal = sum(figx_values_angstroms_abs)
		xavgs.append( xavg )
		xerrors.append( xerror )
		xtotals.append(xtotal)

		figy_values_angstroms = getvaluesfromfile(options, id, 'y')[0]
		figy_vals.update( {k:[id,angletoplot,figy_values_angstroms,'y']})

		figy_values_angstroms_abs = [math.fabs(y) for y in figy_values_angstroms]	
		#yavg = sum([math.fabs(x) for y in figy_values_angstroms])/len(figy_values_angstroms)
		
		yavg = np.mean(figy_values_angstroms_abs)
		yerror = np.std(figy_values_angstroms_abs)
		ytotal = sum(figy_values_angstroms_abs)
		yavgs.append( yavg )
		yerrors.append( yerror )
		ytotals.append(ytotal)

		#c: write out x and y values as a single column txt file for easier plotting with other programs, compared to the original format of files from DE
		if options.savetxts:
			with open( options.path + '/' + id + "_x.txt", 'w' ) as outxfile: outxfile.writelines( [ str(x) +'\n' for x in figx_values_angstroms ] )
			with open( options.path + '/' + id + "_y.txt", 'w' ) as outyfile: outyfile.writelines( [ str(y) +'\n' for y in figy_values_angstroms ] )

		#c: compute |r| for each movie/tiltangle, write to a file, and also write x,y motion file 		
		rs=[]
		rlines=[]
		xandylines=[]
		
		for i in range(len(figx_values_angstroms)):
			r = round(math.sqrt( figx_values_angstroms[i]*figx_values_angstroms[i] + figy_values_angstroms[i]*figy_values_angstroms[i] ),2)
			rs.append(r)
			
			line = str(i) + '\t' + str(r) + '\n'
			rlines.append(line)

			xandyline = str(figx_values_angstroms[i]) + '\t' + str(figy_values_angstroms[i]) + '\n'
			xandylines.append(xandyline)
			
		if options.savetxts:
			with open( options.path +'/' + id + '_r.txt','w' ) as rfile: rfile.writelines( rlines )
			with open( options.path + '/' + id + '_x_vs_y.txt','w' ) as xyfile: xyfile.writelines( xandylines )
		
		figr_vals.update( {k:[id,angletoplot,rs,'r']} )

		#figxy_vals = {k:[id,angletoplot,rs,'r']}

		#compute average |r| and error(s)
		ravg = np.mean( rs )
		rerror = np.std( rs )
		rtotal = sum(rs)
		ravgs.append( ravg )
		rerrors.append( rerror )
		rtotals.append(rtotal)

		avgslist.append( [ angletoplot, ravg, rerror ] )

		if options.verbose:
			print "\nfor img %d xavg=%f, yavg=%f, ravg=%f" %(k,xavg,yavg,ravg)
		
		filesdict.update({ id:[xavg,yavg,ravg,figx_values_angstroms,figy_values_angstroms,rs]} )		
		
		k+=1

	#c: write out the average motion in each direction for each movie in a user-friendly single-column text file
	if options.savetxts:
		with open( options.path + '/' + id + "_x_avgs.txt", 'w' ) as outxavgfile: outxavgfile.writelines( [ str(x) +'\n' for x in xavgs ] )
		with open( options.path + '/' + id + "_y_avgs.txt", 'w' ) as outyavgfile: outyavgfile.writelines( [ str(y) +'\n' for y in yavgs ] )	
		with open( options.path + '/' + id + "_r_avgs.txt", 'w' ) as outravgfile: outravgfile.writelines( [ str(r) +'\n' for r in ravgs ] )

		with open( options.path + '/' + id + "_x_total.txt", 'w' ) as outxtotalfile: outxtotalfile.writelines( [ str(x) +'\n' for x in xtotals ] )
		with open( options.path + '/' + id + "_y_total.txt", 'w' ) as outytotalfile: outytotalfile.writelines( [ str(y) +'\n' for y in ytotals ] )	
		with open( options.path + '/' + id + "_r_total.txt", 'w' ) as outrtotalfile: outrtotalfile.writelines( [ str(r) +'\n' for r in rtotals ] )

	plotdata(options, figx_vals, tag='x', title='Motion in X', xlabel='Frame number', ylabel="Translation in X (" + u"\u212B" + ")")
	plotdata(options, figy_vals, tag='y', title='Motion in Y', xlabel='Frame number', ylabel="Translation in Y (" + u"\u212B" + ")")
	plotdata(options, figr_vals, tag='r', title='Motion |r|', xlabel='Frame number', ylabel="Translation |r| (" + u"\u212B" + ")")
	plotdata(options, figy_vals, tag='x_vs_y', title='Motion in X vs Y', xlabel="Translation in X (" + u"\u212B" + ")", ylabel="Translation in Y (" + u"\u212B" + ")", altxaxisdata=figx_vals)

	fig,ax=resetplot()

	plotavgdata(options, xavgs, angles, xerrors, tag='x_avgs', title='Average X motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))
	plotavgdata(options, yavgs, angles, yerrors, tag='y_avgs', title='Average Y motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))
	plotavgdata(options, ravgs, angles, rerrors, tag='r_avgs', title='Average |r| motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))
	
	plotavgdata(options, xtotals, angles, xerrors, tag='x_total', title='Total X motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))
	plotavgdata(options, ytotals, angles, yerrors, tag='y_total', title='Total Y motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))
	plotavgdata(options, rtotals, angles, rerrors, tag='r_total', title='Total |r| motion per image', xlabel='Tilt angle (degrees)', ylabel="Average translation(" + u"\u212B" + ")", figsize=(10,6))

	E2end(logger)
	
	return


def getangles( options, raworder=False ):
	
	angles = []
	if options.tltfile:
		f = open( options.tltfile, 'r' )
		lines = f.readlines()
		f.close()
		#print "lines in tlt file are", lines
		for line in lines:
			line = line.replace('\t','').replace('\n','')
		
			if line:
				angle = float(line)
				angles.append( angle )
				print "appending angle", angle

		#angles = [a for a in xrange( options.lowesttilt, options.highesttilt, options.tiltstep )]
	else:	
		
		print "There was no .tlt file so I'll generate the angles using lowesttilt=%f, highesttilt=%f, tiltstep=%f" % (options.lowesttilt, options.highesttilt, options.tiltstep)
		generate = floatrange( options.lowesttilt, options.highesttilt, options.tiltstep )
		angles=[ float(x) for x in generate ]
	
	print "BEFORE sorting, angles are", angles
	
	print "negativetiltseries is", options.negativetiltseries
	print "not negativetiltseries", not options.negativetiltseries
	print "raworder", raworder
	print "not raworder", not raworder
	print "not options.negativetiltseries and not raworder", not options.negativetiltseries and not raworder


	if options.negativetiltseries:
		angles.sort()
		print "\n(e2spt_tiltstacker.py)(getangles) AFTER sorting, angles are", angles
	
	elif not raworder:
		angles.sort()
		angles.reverse()
		print "\n(e2spt_tiltstacker.py)(getangles) AFTER REVERSING (ordered from largest to smallest), angles are", angles
	
	return angles


def resetplot(figsize=None):
	plt.clf()

	fig = plt.figure()
	if figsize:
		fig = plt.figure(figsize=(10, 6))
  	ax = fig.add_subplot(1,1,1)

  	plt.rcParams.update({'figure.max_open_warning': 0})

	return fig,ax


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


def getvaluesfromfile(options, fileid, tag):
	
	tagfile = fileid + tag + '.txt'
	
	f=open(tagfile,'r')
	lines=f.readlines()
	f.close()

	values_str = lines[0].replace('\n','').replace('\t',' ').replace('  ',' ').split(' ')[1:]
			
	values_angstroms = [ round(float( value.replace('\n',''))*options.apix,2) for value in values_str ]		
	
	values_pixels = [ round(float( value.replace('\n','')),2) for value in values_str ]

	return values_angstroms,values_pixels


def plotavgdata(options, data, angles, errors, tag='', title='', xlabel='', ylabel='',figsize=(10,6)):
	resolution=150
	if options.highresolution:
		resolution=300

	if options.verbose:
		print "\n(plotavgdata) plotting avgs for fig {}".format(tag)

	fig,ax = resetplot(figsize)
	
	plotfig(options, fig, ax, data, None, 0, False, None, None, title, xlabel, ylabel, angles, errors)
	
	filetosave=options.outputtag + '_' + tag + '_plot.png'
	fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')

	plt.close('fig')
	
	return


def plotdata(options, data, tag, title, xlabel, ylabel, altxaxisdata=None):
	
	resolution=150
	if options.highresolution:
		resolution=300

	if options.verbose:
		print "\n(plotdata) plotting fig {}".format(tag)
	fig,ax = resetplot()
	cpick = resetcolorbar(options.lowestangle,options.highestangle)
	ndata = len(data)
	colorbar = False
	altxaxis = None
	colorstart = options.lowestangle
	colorstep = options.tiltstep

	for count in data:
		if altxaxisdata:
			altxaxis=altxaxisdata[count][2]
		if count==ndata-1:
			colorbar=True
		
		if options.verbose > 9:
			print "\n(plotdata) plotting fig {}, count = {}".format(count,tag)
		
		plotfig(options, fig, ax, data[count][2], cpick, count, colorbar, colorstart, colorstep, title, xlabel, ylabel, altxaxis)		

	filetosave=options.outputtag + '_' + tag + '_plot.png'
	fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')# transparent=True) #, bbox_extra_artists=(lgd,), bbox_inches='tight'
	plt.close('fig')
	print "\n(plotdata) saving figure {}".format(filetosave)

	if options.individualplots:
		colorstep=0
		for count in data:
			colorbar=False
			fig,ax = resetplot()
			print "\n(plotdata) plotting individual plot {}, data {}".format(count, tag)
			plotfig(options, fig, ax, data[count][2], cpick, count, colorbar, colorstart, colorstep, title, xlabel, ylabel, altxaxis)
			filetosave=options.outputtag + '_' + tag + '_plot' +str(count).zfill( len(str(ndata)))+'.png'
			fig.savefig( filetosave, dpi=resolution, bbox_inches='tight')
			plt.close('fig')

	return


def plotfig(options, fig, ax, values, cpick, kplot, colorbar, colorstart, colorstep, title='', xlabel='', ylabel='', altxaxis=None, errors=None ):
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)

	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)

	#if options.miny and options.maxy:
	#	ax.set_ylim(options.miny,options.maxy)
	#	maxr = int(round(math.sqrt(options.maxy*options.maxy + options.miny*options.miny)))	

	n=len(values)
	xaxis=range(n)
	if altxaxis:
		xaxis=altxaxis
	
	colorthis = 'k'
	if colorstart and colorstep:
		colorthis = cpick.to_rgba( colorstart + kplot*colorstep)
	
	if errors:
		lines = {'linestyle': 'None'}
		plt.rc('lines', **lines)
		plt.errorbar(xaxis, values, errors, markersize=8, linewidth=2, fmt='', marker='o', color='k', markerfacecolor=None, markeredgecolor='k', capsize=5, capthick=1)

	ax.plot( xaxis, values, linewidth=2, marker='o', markersize=5, alpha=0.75, color=colorthis, label='Tilt angle')

	ax.set_title(title, fontsize=16, fontweight='bold')
	ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
	ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
	
	if colorbar:
		cbr=plt.colorbar(cpick)
		cbr.set_label("Tilt angle (degrees)",fontsize=18, fontweight='bold')
	
	if not options.individualplots:
		plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
	
	return
	
	
if __name__ == '__main__':
	main()
