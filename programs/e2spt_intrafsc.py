#!/usr/bin/env python

#
# Author: Jesus Galaz, 29/Sep/2014; last update 14/oct/2014
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

#import matplotlib
#matplotlib.use('Agg',warn=False)		 
#import matplotlib.pyplot as plt
import sys
import numpy		 

	 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """
		This program takes a subtomgoram tiltseries (subtiltseries) as extracted with
		e2spt_subtilt.py, and computes the resolution of two volumes reconstructed with
		the even and the odd images in the tilt series. Must be in HDF format.
		Note that the apix in the header must be accurate to get sensible results.
		"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--inputstem", type=str, default='', help="""Default=None. Aligned tilt series. String common to all files to be processed, in the current folder. For example, if you have many subtiltseries named subt00.hdf, subt01.hdf, ...subt99.hdf, you would supply --stem=subt to have all these processed.""")
	
	parser.add_argument('--path',type=str,default='sptintrafsc',help="""Default=sptintrafsc. Directory to save the results.""")
		
	parser.add_argument('--nonewpath',action='store_true',default=False,help="""Default=False. If True, a new --path directory will not be made. Therefore, whatever is sepcified in --path will be used as the output directory. Note that this poses the risk of overwritting data.""")
		
	parser.add_argument('--input',type=str,default='',help="""Default=None. Subtiltseries file file to process. If processing a single file, --inputstem will work too, but you can also just provide the entire filename here --input=subt00.hdf""")
		
	parser.add_argument('--savehalftiltseries',action='store_true',default=False,help="""Default=False. If this parameter is on, the odd and even subtiltseries will be saved.""")
		
	parser.add_argument('--savehalfvolumes',action='store_true',default=False,help="""Default=False. If this parameter is on, the odd and even volumes will be saved.""")
	
	parser.add_argument("--reconstructor", type=str,default="fourier",help="""Default=fourier. The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line to see all options and parameters available. To specify the interpolation scheme for the fourier reconstruction, specify 'mode'. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 'gauss_5_slow', 'gypergeom_5', 'experimental'. For example --reconstructor=fourier:mode=gauss_5 """)
	
	parser.add_argument("--pad2d", type=float,default=0.0,help="""Padding factor to zero-pad the 2d images in the tilt series prior to reconstruction. (The final reconstructed subvolumes will be cropped to the original size).""")

	parser.add_argument("--pad3d", type=float,default=0.0,help="""Default=0.0. Padding factor to zero-pad the reconstruction volume. (The final reconstructed subvolumes will be cropped to the original size).""")
	
	parser.add_argument("--averager",type=str,default="mean.tomo",help="""Default=mean.tomo. The type of averager used to produce the class average.""")
	
	parser.add_argument("--averagehalves",action="store_true", default=False,help="""Default=False. This will averager the even and odd volumes.""")
	
	parser.add_argument("--ppid", type=int, default=-1, help="""default=-1. Set the PID of the parent process, used for cross platform PPID.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--nolog",action="store_true",default=False,help="Default=False. Turn off recording of the command ran for this program onto the .eman2log.txt file") 
	
	(options, args) = parser.parse_args()	
	
	#if options.reconstructor == 'None' or options.reconstructor == 'none':
	#	options.reconstructor = None
	
	#if options.reconstructor and options.reconstructor != 'None' and options.reconstructor != 'none': 
	#	options.reconstructor=parsemodopt(options.reconstructor)
	
	#if options.averager: 
	#	options.averager=parsemodopt(options.averager)
	
	from e2spt_classaverage import sptOptionsParser
	options = sptOptionsParser( options )
	
	logger = E2init( sys.argv, options.ppid )
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	if not options.nonewpath:
		from e2spt_classaverage import sptmakepath
		options = sptmakepath (options, 'sptintrafsc')
	else:
		try:
			findir = os.lisdir( options.path )
		except:
			print "ERROR: The path specified %s does not exist" %( options.path )
			sys.exit()
	
	inputfiles = []
	
	if options.inputstem:
		c = os.getcwd()
		findir = os.listdir( c )

		for f in findir:
			if '.hdf' in f and options.inputstem in f:
				if options.verbose > 8:
					print "\nFound tiltseries!", f
				inputfiles.append( f )			#C:The input files are put into a dictionary in the format {originalseriesfile:[originalseriesfile,volumefile]}

	elif options.input:
		inputfiles.append( options.input )	
	
	for fi in inputfiles:
		
		genOddAndEvenVols( options, fi )
		
		#ret = genOddAndEvenVols( options, fi )
		#volOdd = ret[0]
		#volEven = ret[-1]
			
	E2end(logger)
	
	return


def fscOddVsEven( options, filename, odd, even ):
	
	apix = odd['apix_x']
	
	fsc = odd.calc_fourier_shell_correlation( even )
	third = len( fsc )/3
	xaxis = fsc[0:third]
	fsc = fsc[ third:2*third ]
	saxis = [ x/apix for x in xaxis ]
	
	fscfilename = options.path +'/' + filename.replace('.hdf','_evenOddFSC.txt')
	Util.save_data( saxis[1],saxis[1]-saxis[0],fsc[1:-1], fscfilename )
	
	return fscfilename
	
	
def genOddAndEvenVols( options, fi ):
	
	imgs = []
	
	if isinstance( fi, str ):	
		#hdr = EMData( fi, 0, True)
		nimgs = EMUtil.get_image_count( fi )
		
		for i in range( nimgs ):		
			img = EMData( fi, i )
			imgs.append( img )

	elif isinstance( fi, list ):
		imgs = fi
	
	apix = imgs[0]['apix_x']
	nx = imgs[0]['nx']
	ny = imgs[0]['ny']
	box = nx
	#apix = hdr['apix_x']
	
	mode='gauss_2'
	if options.reconstructor:
		if len(options.reconstructor) > 1:
			if 'mode' in options.reconstructor[-1]:
				mode = options.reconstructor[-1]['mode']
				
				print "\nThe reconstructor mode has been changed from default to", mode
		
	originalboxsize = box
	
	if options.pad3d:
		if options.pad2d:
			if options.pad3d > options.pad2d:
				box = box*options.pad3d
			else:
				box = box*options.pad2d
		else:
			box = box*options.pad3d			
	elif options.pad2d:
		box = box*options.pad2d
			
	rOdd = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})
	rEven = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})

	rOdd.setup()
	rEven.setup()
	
	ko = 0
	oddtilts = options.path + '/' + fi.replace('.hdf', '_ODDTILTS.hdf')
	ke = 0
	eventilts = options.path + '/' + fi.replace('.hdf', '_EVENTILTS.hdf')
	
	iii = 0
	for img in imgs:
				
		#img = EMData( fi, i )
		
		if options.pad2d:
			box2d = img['nx'] * options.pad2d
			img = clip2D( img, box2d )
		
		t = Transform( {'type':'eman', 'az':90, 'alt':img['spt_tiltangle'], 'phi':-90 } )
		img.set_attr('xform.projection',t)
		
		if iii%2:
			print "Print slice inserted into odd", iii
			pmOdd = rOdd.preprocess_slice( img, img['xform.projection'] )
			rOdd.insert_slice( pmOdd,pmOdd['xform.projection'],1.0 )
			
			try:
				if options.savehalftiltseries:
					img.write_image( oddtilts, ko)
			except:
				pass
				
			#print "Adding odd slice", ko
			ko+=1
			
			
		else:
			print "Print slice inserted into even", iii
			pmEven = rEven.preprocess_slice( img,img['xform.projection'] )
			rEven.insert_slice( pmEven,pmEven['xform.projection'],1.0 )
			
			try:
				if options.savehalftiltseries:
					img.write_image( eventilts, ke)
			except:
				pass
				
			#print "Adding even slice", ke
			ke+=1
		iii+=1
		
	recOdd = rOdd.finish(True)
	recOdd['origin_x'] = 0
	recOdd['origin_y'] = 0
	recOdd['origin_z'] = 0
	recOdd['apix_x'] = apix
	recOdd['apix_y'] = apix
	recOdd['apix_z'] = apix
	
	recEven = rEven.finish(True)
	recEven['origin_x'] = 0
	recEven['origin_y'] = 0
	recEven['origin_z'] = 0
	recEven['apix_x'] = apix
	recEven['apix_y'] = apix
	recEven['apix_z'] = apix
	
	#t90 = Transform({ 'type':'eman','alt':180,'az':0,'phi':0 })
	#recOdd.transform( t90 )
	#recEven.transform( t90 )
	#recOdd.rotate( 0,90,0 )
	#recEven.rotate( 0,90,0 )
	
	if options.pad3d or options.pad2d:
		recOdd = clip3D( recOdd, originalboxsize )
		recEven = clip3D( recEven, originalboxsize )
	
	try:	
		if options.savehalfvolumes:
			recOdd.write_image( options.path + '/' + fi.replace('.hdf','_ODDVOL.hdf'), 0 )
			recEven.write_image( options.path + '/' + fi.replace('.hdf','_EVENVOL.hdf'), 0 )
	except:
		pass
	
	fscfilename = fscOddVsEven( options, fi, recOdd, recEven )
		
	try:
		if options.averagehalves:
			avgr = Averagers.get( options.averager[0], options.averager[1] )
			avgr.add_image( recOdd )
			avgr.add_image( recEven )
			
			avg = avgr.finish()
			avg['origin_x'] = 0
			avg['origin_y'] = 0
			avg['origin_z'] = 0
			avg['apix_x'] = apix
			avg['apix_y'] = apix
			avg['apix_z'] = apix
			
			avgfile = options.path + '/AVG.hdf'
			avg.write_image( avgfile, 0 )
	except:
		pass	
		
	return fscfilename
	#return [ recOdd, recEven ]


def clip3D( vol, size ):
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - size)/2, (2*volyc - size)/2, (2*volzc - size)/2, size , size , size)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol


def clip2D( img, size ):
	
	imgxc = img['nx']/2
	imgyc = img['ny']/2
	#imgzc = img['nz']/2
	
	Rimg =  Region( (2*imgxc - size)/2, (2*imgyc - size)/2, 0, size , size , 1)
	img.clip_inplace( Rimg )
	#img.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return img

if '__main__' == __name__:
	main()
	