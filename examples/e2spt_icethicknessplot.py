#!/usr/bin/env python

from EMAN2 import *
import os, numpy, math, sys

def main():
	
	import matplotlib.pyplot as plt
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	This program makes scatter plots of any 2-D dataset provided as 2 columns in a text file.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--files",type=str,default='',help="""Default=None. Comma separated list of coordinates files (in pixels) to produce scatter plots from. E.g., --files=file1.txt,file2.txt,...,fileN.txt. The files must contain lines of x y pairs of values or x y z triads of values separated by a space.""")
	
	parser.add_argument("--filetype",type=str,default='xz',help="""default=xz. If the coordinates files contain only pairs of values, it is assumed that they are x z values. The other valid option is --filetype=yz.""")
	
	parser.add_argument("--stacks",type=str,default='',help="""Default=None. Comma separated list of .hdf stacks of particles boxed with e2spt_boxer.py to read particle coordinates from and produce scatter plots from. E.g., --files=stack1.hdf,stack2.hdf,...,stackN.hdf""")
	
	parser.add_argument("--radius",type=int,default=0,help="""Radius of the particle in pixels.""")
	
	parser.add_argument("--apix",type=float,default=0.0,help="""Apix of the data to which the coordinates in --files or --stacks correspond. MUST provide this for --files if --units is provided.""")
	
	parser.add_argument("--path",type=str,default='',help="""Default=None.""")
	
	parser.add_argument("--yrange",type=str,default='',help="""Default=None. --yrange=min,max. comma separated y1 and y2 values for the min and max limits of the y axis in plots.""")
	
	parser.add_argument("--xrange",type=str,default='',help="""Default=None. --xrange=min,max. comma separated x1 and x2 values for the min and max limits of the y axis in plots.""")
	
	parser.add_argument("--singleplot",action='store_true',default=False,help="""Default=False. If multiple files are provided through --files, they will be plotted in the same .png file, opposed to one .png file per each .txt file in --files.""")
	
	parser.add_argument("--fit",action='store_true',default=False,help="""Default=False. Fit a line through the scatter plot.""")
	
	parser.add_argument("--plotparticleradii",action='store_true',default=False,help="""Default=False. Plots a circle of the particle radius size around the particle center coordinates""")
	
	parser.add_argument("--units",type=str,default='pixels',help="""default=pixels. Other valid options are 'angstroms','nanometers','micrometers' or 'microns'. Make sure the apix values in the HDF files provided through --stacks is correct. If providing --files, the separately provide --apix""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	#parser.add_argument("--parallel",type=str,default='',help="""the program will detect the number of cores available and use threaded parallelism by default. To use only one core, supply --parallel=thread:1. For MPI on clusters, see parallelism at http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
		
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
		
	(options, args) = parser.parse_args()
	
	
	if not options.files and not options.stacks:
		print "\nERROR: you must provide --files, --stacks, or both."""
		sys.exit()
	
	
	if options.singleplot:
		if options.files and options.stacks:
			print """\nERROR: To use --singleplot, you must provide either --stakcs or --files separately"""
			
	'''
	Make a directory where to store output files
	'''
	#from e2spt_classaverage import sptmakepath
	#options = sptmakepath(options,'scatter_plots')
	
	if options.path:
		os.mkdir( options.path )
	
	#from e2spt_classaverage import writeParameters
	#writeParameters(options,'e2spt_test.py', 'scatter_plots')
	
	
	#if not options.parallel:
	#	import multiprocessing
	#	nparallel = multiprocessing.cpu_count()
	#	options.parallel = 'thread:' + str(nparallel)
	#	print "\nfound %d cores" %(nparallel)
	#	print "setting --parallel to", options.parallel

	
	logger = E2init(sys.argv, options.ppid)
	
	aunits = ['angstroms','Angstroms','ANGSTROMS']
	nunits = ['nanometers','NANOMETERS','Nanometers']
	munits = ['micrometers','MICROMETERS','micrometers','MICRONS','microns','Microns']
	apixunits = aunits + nunits + munits
	
	apix = 1.0
	if options.files:
		if options.units in apixunits:
			if not options.apix:
				print "\nERROR: You must provide --apix in angstroms/pixel when --units is 'micrometers','nanometers' or 'angstroms'"""
				sys.exit()
			else:
				apix = options.apix
				
	unitslabel = '(pixels)'
	if options.units in aunits:
		unitslabel = '(angstroms)'
	elif options.units in nunits:
		unitslabel = '(nanometers)'
	elif options.units in munits:
		unitslabel = '(micrometers)'
	else:
		unitslabel = '(pixels)'
	
	print "set unitslabel", unitslabel
	print "because options.units",options.units
		
	
		
	files = options.files.split(',')
	
	filesdata = {}		
	
	pad = 0
	
	if files and options.files:
		for f in files:
			xs = []
			ys = []
			zs = []
	
			ff = open(f,'r')
			lines = ff.readlines()
			nlines = len(lines)
			ff.close()
		
			for line in lines:
				parsed = line.replace('\n','').split(' ')
				x = float( parsed[0] )
				y = float( parsed[1] )
				z = 0
				if len(parsed) > 2:
					z = float( parsed[2] )
			
				if options.units in aunits:
					x = x*apix
					y = y*apix
					z = z*apix
				
					if options.radius:
						pad = 2*options.radius * apix
					else:
						pad = 50
					
				elif options.units in nunits:
					x = x*apix/10
					y = y*apix/10
					z = z*apix/10
				
					if options.radius:
						pad = 2*options.radius * apix/10
					else:
						pad = 5
					
				elif options.units in munits:
					x = x*apix/10000
					y = y*apix/10000
					z = z*apix/10000
				
					if options.radius:
						pad = 2*options.radius * apix/10000
					else:
						pad = 0.5
				else:
					if options.radius:
						pad = 2*options.radius
					else:
						pad = 10
						
				#if len(parsed) > 2:
				#	z = float( parsed[2] )
			
				xs.append(x)
				ys.append(y)
				if z:
					zs.append(z)
		
			filesdata.update( { f: [xs,ys,zs] } )
		
	
	stacks =  options.stacks.split(',')
	
	stacksdata = {}
	
	if stacks and options.stacks:
		for stack in stacks:
			n = EMUtil.get_image_count( stack )
			xs = []
			ys = []
			zs = []
			for i in range(n):
				hdr = EMData( stack, i, True )
				coords = hdr['ptcl_source_coord']
				apix = hdr['apix_x']
				x = int( coords[0])
				y = int( coords[1])
				z = int( coords[2])
		
				print "\nfor particle %d XZ coords are x=%d, z=%d" %(i,x,z)
			
				if options.units in aunits:
					x = x*apix
					y = y*apix
					z = z*apix
				
					if options.radius:
						pad = 2*options.radius * apix
					else:
						pad = 50
			
				elif options.units in nunits:
					x = x*apix/10
					y = y*apix/10
					z = z*apix/10
				
					if options.radius:
						pad = 2*options.radius * apix/10
					else:
						pad = 5
					
					
				elif options.units in munits:
					x = x*apix/10000
					y = y*apix/10000
					z = z*apix/10000
				
					if options.radius:
						pad = 2*options.radius * apix/10000
					else:
						pad = 0.5
				
				else:
					if options.radius:
						pad = 2*options.radius
					else:
						pad = 10
				
				
				print "\nfor particle %d XY coords are x=%d, y=%d in %d" %(i,x,y,options.units)
		
				print "\nfor particle %d XZ coords are x=%d, z=%d in %d" %(i,x,z,options.units)
		
				xs.append(x)
				ys.append(y)
				zs.append(z)
		
			filesdata.update( { stack: [xs,ys,zs] } )
		
	
	xaxislabel = 'X axis ' + unitslabel
	yaxislabel = 'Y axis ' + unitslabel
	zaxislabel = 'Z axis ' + unitslabel
	
	print "unitslabel",unitslabel
	
	if stacksdata:
		for data in stacksdata:
			xdata = stacksdata[ data ][0]
			zdata = stacksdata[ data ][2]
			stackname = os.path.splitext( os.path.basename(data) )[0]
			tag = '_xz'
			title = stackname + tag
			plotter( options, xdata, zdata, xaxislabel, zaxislabel, title, pad )
		
			if not options.singleplot:
				plotname = title + '.png'
				if options.path:
					plotname = options.path + '/' + plotname
				plt.savefig( plotname )
				plt.clf()
	
		if options.singleplot:
			singleplotname = 'stacks_plot_xz.png'
			if options.path:
					singleplotname = options.path + '/' +singleplotname
			plt.savefig(singleplotname)
			plt.clf()
	
		
		for data in stacksdata:
			ydata = stacksdata[ data ][1]
			zdata = stacksdata[ data ][2]
			stackname = os.path.splitext( os.path.basename(data) )[0]
			tag = '_yz'
			title = stackname + tag
			plotter( options, ydata, zdata, yaxislabel, zaxislabel, title, pad )
		
			if not options.singleplot:
				plotname = title + '.png'
				if options.path:
					plotname = options.path + '/' + plotname
				plt.savefig( plotname )
				plt.clf()
	
		if options.singleplot:
			singleplotname = 'stacks_plot_yz.png'
			if options.path:
					singleplotname = options.path + '/' +singleplotname
			plt.savefig(singleplotname)
			plt.clf()

	print "\nline 301"
	print "\noptions.xrange, options.yrange", options.xrange, options.yrange
	if filesdata:
		
		data3count = 0
		
		tag = '_xz'
		if options.filetype == 'yz': 	
			tag = '_yz'
		
		for data in filesdata:
			data1 = filesdata[ data ][0]
			data2 = filesdata[ data ][1]
			data3 = filesdata[ data ][2]		
						
			if not data3:
				filesname = os.path.splitext( os.path.basename(data) )[0]
				
				
				datax = data1
				datay = data2
				
				if options.filetype == 'yz': 	
					tag = '_yz'
					datax = data2
					datay = data3
				
				title = filesname + tag	
				
				plotter( options, data1, data2, xaxislabel, zaxislabel, title, pad )
			
				if not options.singleplot:
					plotname = title + '.png'
					if options.path:
						plotname = options.path + '/' + plotname
					plt.savefig( plotname )
					plt.clf()
			
			elif data3:
				data3count += 1
			
		if not data3count:
			if options.singleplot:
				singleplotname = 'files_plot_' + tag + '.png'
			
				if options.path:
					singleplotname = options.path + '/' +singleplotname
				plt.savefig(singleplotname)
				plt.clf()
			
			

		if data3count == len( filesdata ):
			
			for data in filesdata:
				xdata = filesdata[ data ][0]
				zdata = filesdata[ data ][2]
				filesname = os.path.splitext( os.path.basename(data) )[0]
				tag = '_xz'
				title = filesname + tag		
				plotter( options, xdata, zdata, xaxislabel, zaxislabel, title, pad )
		
				if not options.singleplot:
					plotname = title + '.png'
					if options.path:
						plotname = options.path + '/' + plotname
					plt.savefig( plotname )
					plt.clf()
	
			if options.singleplot:
				singleplotname ='files_plot_xz.png'
				if options.path:
					singleplotname = options.path + '/' +singleplotname
				plt.savefig(singleplotname)
				plt.clf()
	
		
			for data in filesdata:
				ydata = filesdata[ data ][1]
				zdata = filesdata[ data ][2]
				filesname = os.path.splitext( os.path.basename(data) )[0]
				tag = '_yz'
				title = filesname + tag	
				plotter( options, ydata, zdata, yaxislabel, zaxislabel, title, pad )
		
				if not options.singleplot:
					plotname = title + '.png'
					if options.path:
						plotname = options.path + '/' + plotname
					plt.savefig( plotname )
					plt.clf()
	
			if options.singleplot:
				singleplotname = 'files_plot_yz.png'
				if options.path:
					singleplotname = options.path + '/' +singleplotname
				
				plt.savefig(singleplotname)
				plt.clf()
			
	'''		
				tag = '_xz'
				plotter( data1, data2, xaxislabel, zaxislabel )
				
				if not options.singleplot:
					plotname = options.path + '/' + os.path.splitext( os.path.basename(data) )[0] + tag + '.png'
					plt.savefig( plotname )
					plt.clf()
				
				tag = '_yz'
				plotter( data2, data3, xaxislabel, zaxislabel )
				
				if not options.singleplot:
					plotname = options.path + '/' + os.path.splitext( os.path.basename(data) )[0] + tag + '.png'
					plt.savefig( plotname )
					plt.clf()
			
			if data3 and options.singleplot:
					
				
		if options.singleplot:
			singleplotname = options.path + '/stacks_plot_XY.png'
			plt.savefig(singleplotname)
			plt.clf()
	
	
		zplot=0
		for data in filesdata:
			zdata = filesdata[ data ][2]
			if zdata:
				ydata = filesdata[ data ][1]
							
				plotter( ys, zs, yaxislabel, zaxislabel )
		
				if not options.singleplot:
					plotname = options.path + '/' + os.path.splitext( os.path.basename(data) )[0] + '.png'
					plt.savefig( plotname )
					plt.clf()
				
				zplot+=1
			
		if zplot:
			if options.singleplot:
				singleplotname = options.path + '/stacks_plot_XZ.png'
				plt.savefig(singleplotname)
				plt.clf()
	'''

	E2end(logger)
	sys.stdout.flush()
	
	return


def plotter(options,xaxis,yaxis,xaxislabel,yaxislabel,title,pad=0):
	print "\npad is", pad
	import matplotlib

	matplotlib.use('Agg',warn=False)
	matplotlib.use("TkAgg",warn=False)

	import matplotlib.pyplot as plt
	import pylab

	
	'''
	FORMAT AXES
	'''
				
	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	fig = plt.figure()
	
	fig = plt.figure(figsize=(30, 5))
	
	plt.axis('equal')
	
	ax = fig.add_subplot(111)
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)
	
	#print "max y is", max(yaxis)
	#pylab.ylim([0,max(yaxis)+ 100])

	#print "max x is", max(xaxis)
	
	#pylab.xlim([0,max(xaxis)+ 100])
	
	#if 'pixels' in xaxislabel or 'pixels' in yaxislabel:
	#	if options.radius:
	#		pad = options.radius
	
	maxy = max(yaxis)
	miny = min(yaxis)
	
	print "\nmax y is", maxy
	ylim1 = miny - pad
	ylim2 = maxy + pad
	
	
	if options.yrange:
		ylim1 = int( options.yrange.split(',')[0] )
		ylim2 = int( options.yrange.split(',')[1] )
		
	pylab.ylim([ylim1, ylim2])
	
	print 'yrange', ylim1,ylim2
	
	print "\nmax x is", max(xaxis)
	
	maxx = max(xaxis)
	minx = min(xaxis)
	
	xlim1 = minx - pad
	xlim2 = maxx + pad
	if options.xrange:
		xlim1 = int( options.xrange.split(',')[0] )
		xlim2 = int( options.xrange.split(',')[1] )
		
	
	pylab.xlim([xlim1, xlim2])
	
	print 'xrange', xlim1,xlim2
	
	
	deltaz = maxy-miny
	if options.radius:
		deltaz += pad
	
	if 'angstroms' in xaxislabel:
		deltaz = int(deltaz)
	elif 'nanometers' in xaxislabel:
		deltaz = round(deltaz,1)
	elif 'micrometers' or 'microns' in xaxislabel:	
		deltaz = round(deltaz,2)	
	else:
		deltaz = int(round(deltaz))
	
	
	#deltazlabel = u"\N{GREEK CAPITAL LETTER DELTA}" + '=' + str(deltaz)
	
	#deltazlabel = u"\u0394" + "Z=" + str( deltaz )
	#deltazlabel = chr(0x0394) + "Z=" + str( deltaz )
	deltazlabel = r"$\Delta$" + "Z=" + str( deltaz )
	
	if options.radius:
		deltazlabel = "ice thickness = " + str( deltaz )
	
	if 'angstroms' in xaxislabel:
		deltazlabel +=  u"\u00c5"
	
	#if 'angstroms' in xaxislabel:
	#	deltazlabel +=  r"$\Angstrom$"
	elif 'nanometers' in xaxislabel:
		deltazlabel += "nm"
	#elif 'micrometers' or 'microns' in xaxislabel:	
	#	deltazlabel += u"\u03bc"
	elif 'micrometers' in xaxislabel or 'microns' in xaxislabel:
		#print "setting mu because xaxislabel is", xaxislabel	
		deltazlabel += r"$\mu$"
	else:
		deltazlabel += ' pix'
		
	ax.set_xlabel(xaxislabel, fontsize=18, fontweight='bold')
	ax.set_ylabel(yaxislabel, fontsize=18, fontweight='bold')
	ax.set_title(title, fontsize=18, fontweight='bold')
	
	
	#plt.scatter(xaxis,yaxis,alpha=1,zorder=1,s=20,color='k')
	centerdistribution = plt.scatter(xaxis,yaxis,alpha=1,zorder=1,s=20,color='k')
	
	#print "deltazlabel is"
	#print deltazlabel
	#plt.legend([centerdistribution,(centerdistribution)], [deltazlabel], loc=1)
	
	
	if options.fit:
		m, b = numpy.polyfit(xaxis, yaxis, 1)
	
		anglecalc = round( math.degrees(numpy.arctan( m )), 3 )
		
		print "\nm and b", m,b
		
		anglelabel = str(anglecalc) + ' deg'
		print '\nanglelabel',anglelabel
		xarray = numpy.array( xaxis )
		#fitline = plt.plot(xaxis, m*xarray + b, '--', linewidth=3, alpha=0.333,color='k',linestyle='--', label=anglelabel)
		
		fitline, = plt.plot(xaxis, m*xarray + b, '--', linewidth=3, alpha=0.5,color='k',linestyle='--')
		
		if not options.plotparticleradii:
			plt.legend([fitline,fitline] ,[anglelabel,deltazlabel], loc=2)
		
		
		
	if options.plotparticleradii:
		if options.radius:
			rad = pad/2		#this should be at the proper scale already from main function
			
			axes = pylab.axes()
			
			for x, y in zip(xaxis, yaxis):
				circle = pylab.Circle( (x,y), radius = rad, facecolor='none',edgecolor='b', alpha=0.3, linewidth=3)
				axes.add_patch(circle)
				
			if not options.fit:
				plt.legend( [circle,(circle)],[deltazlabel] )
			
			else:
				plt.legend( [circle,fitline],[deltazlabel,anglelabel],loc=2 )
		else:
			print "\nERROR: --plotparticleradii requires --radius"
			
	#name = 'scatterplot.png'
	
	#plt.savefig(name)
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nSSSSSSSSSSSSSSSSSSSS\nSaved plot"
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	
	#plt.clf()
	
	return	
	
	
if __name__ == '__main__':
	main()
