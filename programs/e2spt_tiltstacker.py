#!/usr/bin/env python																																																																																																																																																																																																																																																																																																																			#!/usr/bin/python2.7

#====================
#Author: Jesus Galaz-Montoya 2/20/2013 , Last update: January/15/2016
#====================
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


from optparse import OptionParser

from EMAN2 import *
import sys
import math

print "imported stuff"

def main():

	usage = """e2spt_tiltstacker.py <options> . 
	The options should be supplied in "--option=value" format, 
	replacing "option" for a valid option name, and "value" for an acceptable value for that option. 
	
	This program operates in 3 different modes:
	1) It can STACK individual .dm3, .tiff or .hdf images into an .mrc (or .st) stack,
	by supplying a common string to all the images to stack via --stem2stack.
	It must be run in a directory containing the numbered images only.
	It also generates a .rawtlt file with tilt angle values if --lowerend, --upperend and --tiltstep are provided.
	
	2) It can UNSTACK a tilt series into individual files (either all the images, or selected
	images, controlled through the --exclude or --include parameters).
	
	3) It can RESTACK a tilt series; that is, put together a new tilt series that excludes/includes
	specific images
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'sptstacker';
		for example, sptstacker_02 will be the directory by default if 'sptstacker_01' 
		already exists.""")

	parser.add_argument("--stem2stack", type=str, default='', help="""String common to all 
		the files to put into an .st stack, which is in .MRC format; for example, --stem2stack=.hdf 
		will process all .hdf files in the current directory.
		If not specified, all valid EM imagefiles in the current directory will be put into 
		an .st stack.""")
	
	parser.add_argument("--tltfile",type=str,default='',help="""".tlt file IF unstacking an
		aligned tilt series with --unstack=<stackfile> or restacking a tiltseries with
		--restack=<stackfile>""")
		
	parser.add_argument("--invert",action="store_true",default=False,help=""""This 
		will multiply the pixel values by -1.""")
	
	parser.add_argument("--stackregardless",action="store_true",default=False,help=""""Stack
		images found with the common string provided through --stem2stack, even if the
		number of images does not match the predicted number of tilt angles.""")
	
	parser.add_argument("--outmode", type=str, default="float", help="""All EMAN2 programs 
		write images with 4-byte floating point values when possible by default. This allows 
		specifying an alternate format when supported: float, int8, int16, int32, uint8, 
		uint16, uint32. Values are rescaled to fill MIN-MAX range.""")
	
	parser.add_argument("--bidirectional",action='store_true',default=False,help="""This will
		assume the first image is at 0 degrees and will stack images from --lowerend through 0, 
		and then will stack the rest from 0+tiltstep throgh --upperend. 
		If --negativetiltseries is supplied, images will be stacked from --upperend through 0, then 
		from 0-tiltstep through --lowerend.""")
	
	parser.add_argument("--lowesttilt",type=float,default=0.0,help="""Lowest tilt angle.
		If not supplied, it will be assumed to be -1* --tiltrange.""")
	
	parser.add_argument("--highesttilt",type=float,default=0.0,help="""Highest tilt angle.
		If not supplied, it will be assumed to be 1* --tiltrange.""")
	
	parser.add_argument("--tiltrange",type=float,default=0.0,help="""If provided, this
		will make --lowesttilt=-1*tiltrange and --highesttilt=tiltrage.
		If the range is asymmetric, supply --lowesttilt and --highesttilt directly.""")
	
	parser.add_argument("--tiltstep",type=float,default=0.0,help="""Step between tilts.
		Required if using --stem2stack.""")
	
	parser.add_argument("--clip",type=str,default='',help="""Resize the 2-D images in the
		tilt series. If one number is provided, then x and y dimensions will be made the same.
		To specify both dimensions, supply two numbers, --clip=x,y. Clipping will be about
		the center of the image.""")
			
	parser.add_argument("--apix",type=float,default=0.0,help="""True apix of images to be 
		written on final stack.""")
	
	parser.add_argument("--unstack",type=str,default='',help=""".hdf, or 3D .st, .mrc, 
		.ali, or .mrcs stack file to unstack.
		This option can be used with --include or --exclude to unstack only specific images.
		Recall that the FIRST image INDEX is 0 (but unstacked image will be numbered from 1). 
		--exclude=1,5-7,10,12,15-19 will exclude images 1,5,6,7,10,12,15,16,17,18,19""""")
	
	parser.add_argument("--restack",type=str,default='',help=""".hdf, or 3D .st, .mrc, 
		.ali, or .mrcs stack file to restack.
		This option can be used with --include or --exclude to unstack only specific images.
		Recall that the FIRST image INDEX is 0 (but unstacked image will be numbered from 1). 
		--exclude=1,5-7,10,12,15-19 will exclude images 1,5,6,7,10,12,15,16,17,18,19""""")
	
	parser.add_argument("--mirroraxis",type=str,default='',help="""Options are x or y, and the
		mirrored copy of the 2-D images will be generated before being put into the tilt series.""")
	
	parser.add_argument("--exclude",type=str,default='',help="""Comma separated list of numbers
		corresponding to images to exclude. --unstack or --restack must be supplied. 
		You can also exclude by ranges. For example:
		Recall that the FIRST image INDEX is 0. 
		--exclude=1,5-7,10,12,15-19 will exclude images 1,5,6,7,10,12,15,16,17,18,19""")
		
	parser.add_argument("--include",type=str,default='',help="""Comma separated list of numbers
		corresponding to images to include (all others will be excluded). 
		--unstack or --restack must be supplied. 
		Recall that the FIRST image INDEX is 0. 
		--include=1,5-7,10,12,15-19 will include images 1,5,6,7,10,12,15,16,17,18,19""")


	#parser.add_argument("--negativetiltseries",action='store_true',default=False,help="""This indicates that the tilt series goes from -tiltrange to +tiltrange, or 0 to -tiltrange, then +tiltstep to +tiltrange if --bidirectional is specified.""")
	
	#parser.add_argument("--negative",action='store_true',default=False,help="""This indicates that the tilt series goes from -tiltrange to +tiltrange, or 0 to -tiltrange, then +tiltstep to +tiltrange if --bidirectional is specified.""")
	
	parser.add_argument("--negativetiltseries",action='store_true',default=False,help="""This indicates that the tilt series goes from -tiltrange to +tiltrange, or 0 to -tiltrange, then +tiltstep to +tiltrange if --bidirectional is specified.""")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	
	
	print "--negativetiltseries", options.negativetiltseries
	
	if options.exclude and options.include:
		print "\nERROR: Supplied either exclude or include. Cannot supply both at the same time."
		sys.exit()
	
	print "\nLogging"
	logger = E2init(sys.argv, options.ppid)
	
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options, 'sptstacker')

	options.path = os.getcwd() + '/' + options.path
	
	tiltstoexclude = options.exclude.split(',')	
	
	if options.stem2stack and options.tiltstep == 0.0:
		print "ERROR: --tiltstep required when using --stem2stack"
		sys.exit()
	
	if options.lowesttilt == 0.0 and options.tiltrange:
		options.lowesttilt = -1 * options.tiltrange
		
	if options.highesttilt == 0.0 and options.tiltrange:
		options.highesttilt = options.tiltrange

	if options.unstack:
		if options.tltfile:
			usntacker( options )
		else:
			print "ERROR: --tltfile required when using --unstack"
			sys.exit()

	elif options.restack:
		if options.tltfile:
			restacker( options )
			angles = getangles( options, True )			#Second parameter enforces to keep the 'raw order' of the input file. Otherwise, this function returns angles from -tiltrange to +tiltrange if --negativetiltseries is supplied; from +tiltrange to -tiltrange otherwise

			#finalangles = list(angles)
			#anglestoexclude = []
			
			print "\n\nthere are these many angles", len(angles)
			#if tiltstoexclude:
			#	for tilt in tiltstoexclude:
			#		anglestoexclude.append( angles[ int(tilt) ] )
			#		finalangles.remove( angles[ int(tilt) ] )
			#	#for ax in anglestoexclude:
			#	#	finalangles.remove( ax )
			#
			#	print "\n\nthere are these many angles to exclude",len(anglestoexclude)
			#	print "\nexcluded angles",anglestoexclude
			#	
			#	#finalangles = list( set(angles) - set(anglestoexclude) )
				
				
			#print "\nthere are these many final angles",len(finalangles)
			#print "\nfinal angles are", finalangles
		
			
			writetlt(angles,options,True)
		else:
			print "ERROR: --tltfile required when using --restack"
			sys.exit()

	else:
		kk=0
		intilts = findtiltimgfiles( options )
		
		print "\nWill organize tilt imgs found"
		intiltsdict = organizetilts( intilts, options )		#Get a dictionary in the form { indexintiltseries:[ tiltfile, tiltangle, damageRank ]},
		print "\nDone organizing tilt imgs"					#where damageRank tells you the order in which the images where acquired
															#regardless of wether the tilt series goes from -tiltrange to +tiltrange, 
															#or 0 to -tiltrange then +tiltstep to +tiltrange, or the opposite of these 
		outstackhdf = options.path + '/stack.hdf' 
		
		minindx = min(intiltsdict)
		print "minindx is", minindx
		print "getting size from any first image, intiltsdict[ minindx ][0]", intiltsdict[ minindx ][0]
		
		hdr = EMData( intiltsdict[minindx][0], 0, True )
		nx = hdr['nx']
		ny = hdr['ny']
		print nx,ny
		
		
		print "\nOutstack is", outstackhdf
		
					
		
		#orderedindexes = []
		#for index in intiltsdict:
		#	orderedindexes.append( index )
		
		#orderedindexes.sort()
		
		for index in intiltsdict:
			
			if str(index) not in tiltstoexclude:
				intiltimgfile =	intiltsdict[index][0]
				
				print "\nat index %d we have image %s, collected in this turn %d" %( index, intiltsdict[index][0], intiltsdict[index][-1] )
				intiltimg = EMData( intiltimgfile, 0 )
			
				tiltangle = intiltsdict[index][1]
				intiltimg['spt_tiltangle'] = tiltangle
			
				damageRank = intiltsdict[index][2]
				intiltimg['damageRank'] = damageRank
			
				if options.invert:
					intiltimg.mult(-1)
				intiltimg.write_image( outstackhdf, -1 )
				print "\nWrote image index", index
		
		if options.clip:
			clip = options.clip.split(',')
			
			shiftx = 0
			shifty = 0
			if len( clip ) == 1:
				clipx = clipy = clip[0]
			
			if len( clip ) == 2:
				clipx = clip[0]
				clipy = clip[1]
			
			if len( clip ) == 4:
				clipx = clip[0]
				clipy = clip[1]
				shiftx = clip[2]
				shifty = clip[3]
			
			tmp = options.path + '/tmp.hdf'
			cmdClip = 'e2proc2d.py ' + outstackhdf + ' ' + tmp + ' --clip=' + clipx + ',' + clipy
			
			if shiftx:
				xcenter = int( round( nx/2.0 + float(shiftx)))
				cmdClip += ',' + str(xcenter)
			if shifty:
				ycenter = int( round( ny/2.0 + float(shifty)))
				cmdClip += ',' + str(ycenter)
			
			cmdClip += ' && rm ' + outstackhdf + ' && mv ' + tmp + ' ' + outstackhdf

			print "\n(e2spt_tiltstacker.py)(main) cmdClip is", cmdClip
			p = subprocess.Popen( cmdClip , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()
			
			if options.verbose > 9:
				print "\nFeedback from cmdClip:"
				print text
			
			
		outtilthdr = EMData(outstackhdf,0,True)
		currentapix = outtilthdr['apix_x']
		if float(options.apix) and float(options.apix) != float(currentapix):
			print "\nFixing apix"
			cmdapix = 'e2fixheaderparam.py --input=' + outstackhdf + ' --stem=apix --valtype=float --stemval=' + str( options.apix )
			
			print "\n(e2spt_tiltstacker.py)(main) cmdapix is", cmdapix
			p = subprocess.Popen( cmdapix , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()
			
			if options.verbose > 9:
				print "\nFeedback from cmdapix:"
				print text
			
			
		outstackst = outstackhdf.replace('.hdf','.st')
		stcmd = 'e2proc2d.py	' + outstackhdf + ' ' + outstackst + ' --twod2threed'
		if options.outmode != 'float':
			stcmd += ' --outmode=' + options.outmode + ' --fixintscaling=sane'
			
		if options.apix:
			stcmd += ' --apix=' + str(options.apix)
			stcmd += ' && e2fixheaderparam.py --input=' + outstackst + ' --stem=apix --valtype=float --stemval=' + str( options.apix ) + ' --output=' + outstackst.replace('.st','.mrc') + " && mv " +  outstackst.replace('.st','.mrc') + ' ' + outstackst
			
		
		print "\n(e2spt_tiltstacker.py)(main) stcmd is", stcmd	
		p = subprocess.Popen( stcmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text = p.communicate()	
		p.stdout.close()
		
		if options.verbose > 9:
			print "\nFeedback from stcmd:"
			print text
		
		if options.mirroraxis:
			print "\nMirroring across axis", options.mirroraxis
			mirrorlabel = options.mirroraxis.upper()
			outstackstmirror = outstackst.replace('.st','_mirror'+ mirrorlabel+ '.st')

			cmdMirror = 'e2proc2d.py ' + outstackst + ' ' + outstackstmirror + ' --process=xform.mirror:axis=' + options.mirroraxis
			
			if options.outmode != 'float':
				cmdMirror += ' --outmode=' + options.outmode + ' --fixintscaling=sane'
			
			print "options.apix is", options.apix
			if options.apix:
				cmdMirror += ' --apix=' + str(options.apix)
				cmdMirror += ' && e2fixheaderparam.py --input=' + outstackstmirror + ' --stem=apix --valtype=float --stemval=' + str( options.apix ) + ' --output=' + outstackstmirror.replace('.st','.mrc') + " && mv " +  outstackstmirror.replace('.st','.mrc') + ' ' + outstackstmirror
				
				print "added fixheaderparam to cmdMirror!"
			
			print "cmdMirror is", cmdMirror
			p = subprocess.Popen( cmdMirror , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			
			if options.verbose > 9:
				print "\nFeedback from cmdMirror:"
				print text
				p.stdout.close()
	
	E2end( logger )
	return


"""
Function to determine which files in the current directory should be used
"""
def findtiltimgfiles( options ):
	
	currentdir = os.getcwd()
	filesindir = os.listdir( currentdir )
	
	intilts = []
	#k=0
	
	if options.stem2stack:
		intilt = ''	
		for f in filesindir:
			if options.stem2stack in f:
				if '.txt' not in f and '.db' not in f:
					print "potential tilt image", f
					if '.dm3' in f[-4:] or '.DM3' in f[-4:] or '.tif' in f[-4:] or '.TIF' in f[-4:] or '.MRC' in f[-4:] or '.mrc' in f[-4:] or '.hdf' in f[-4:]: 
						print "\nvalid file", f
						intilt = f
						if intilt:
							intilts.append( intilt )
					else:
						print "format not valid. Needs to be .mrc, .hdf, .dm3, or .tif"
					
					#k will serve to compensate for damage.
					#It indicates the order in which data collection occured
	print "\n(e2spt_tiltstacker.py)(findtiltimgfiles) These many img files were found", len( intilts )		
	intilts.sort()
	print "\nand they've been sorted"
	return intilts


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
				print "appending angle", 
	else:
		#angles = [a for a in xrange( options.lowesttilt, options.highesttilt, options.tiltstep )]
		
		
		print "There was no .tlt file so I'll generate the angles using lowesttilt=%f, highesttilt=%f, tiltstep=%f" % (options.lowesttilt, options.highesttilt, options.tiltstep)
		generate = floatrange( options.lowesttilt, options.highesttilt, options.tiltstep )
		angles=[ float(x) for x in generate ]
	
	print "BEFORE sorting, angles are", angles
	angles.sort()
	print "\n(e2spt_tiltstacker.py)(getangles) AFTER sorting, angles are", angles
	
	if not options.negativetiltseries and not raworder:
		angles.reverse()
		print "\n(e2spt_tiltstacker.py)(getangles) AFTER REVERSING, angles are", angles
	
	

	return angles


def writetlt( angles, options, raworder=False ):
	
	print "(writetlt) these many angles", len(angles)
	angless = list( angles )
	
	if not raworder:
		angless.sort()
	
	if not options.negativetiltseries and not raworder:
		angless.reverse()
		
	f = open( options.path + '/stack.rawtlt','w')
	lines = []
	
	k=0
	anglestoexclude = []
	tiltstoexclude = options.exclude.split(',')				
	for a in angless:
		if str(k) not in tiltstoexclude:
			line = str(a) + '\n'
			lines.append( line )
		else:
			anglestoexclude.append( angless[k] )
		
		k+=1
	print "\nthese many angles excluded",len(anglestoexclude)
	print "\nexcluded angles",anglestoexclude	
	f.writelines( lines )
	f.close()
		
	return


def organizetilts( intilts, options ):
	
	intiltsdict = {}
	angles = getangles( options )			#This returns angles from -tiltrange to +tiltrange if --negativetiltseries is supplied; from +tiltrange to -tiltrange otherwise
	orderedangles = list( angles )
	
	#print "\n(organizetilts) angles are", angles
	
	zeroangle = min( [ math.fabs(a) for a in angles] )		#Find the angle closest to 0 tilt, in the middle of the tiltseries
	indexminangle = None
	try:
		indexminangle = angles.index( zeroangle )			#The middle angle (0 tilt) could have been positive
		zeroangle = angles[ indexminangle ]
	except:
		indexminangle = angles.index( -1*zeroangle )		#Or negative. Either way, we find its index in the list of angles
		zeroangle = angles[ indexminangle ]

	
	if not options.tltfile:
		
		
		if options.bidirectional:
			
			
			print "\nzeroangle=%f, indexminangle=%d" %( zeroangle, indexminangle )
			if not options.negativetiltseries:
				print "\nNegative tilt series is OFF. This is a POSITIVE tilt series"
				secondhalf = angles[ indexminangle: len(angles) ]	#This goes from zero to lowest tilta angles, i.e., -tiltrange
				#print "Firsthalf is", firsthalf
				#print "because angles are", angles
				firsthalf =  angles[ 0:indexminangle ]				#This should go from most positive angle or +tiltrange to zero (without including it)
				#secondhalf.sort()									#We order this list to go from 0-tiltstep to the most negative angle, -tiltrange
				#secondhalf.reverse()
			
			elif options.negativetiltseries:
				print "T\nhis is a NEGATIVE tiltseries"
				firsthalf = angles[ 0:indexminangle+1 ]				#This goes from the most negative angle to zero (INCLUDING it)
				#firsthalf.sort()									#We order this list to go from 0 to -tiltrange
				#firsthalf.reverse()
				
				secondhalf = angles[ indexminangle+1: len(angles) ]	#This goes from 0+tiltstep to the most positive angle, that is, +tiltrange
			
			orderedangles = firsthalf + secondhalf
			print "\n(organizetilts) Ordered angles are", orderedangles
			print "\n(organizetilts) Because firsthalf is", firsthalf
			print "\n(organizetilts) and secondhalf is", secondhalf
			
			writetlt( orderedangles, options )
	
	else:
		writetlt( angles, options )
	
	angles.sort()
	if not options.negativetiltseries:			#Change angles to go from +tiltrange to -tiltrange if that's the order of the images
		#orderedangles.sort()
		#orderedangles.reverse()
		
		#angles.sort()
		angles.reverse()
		
		#print "However, after reversal, they are", orderedangles
		#print "and angles are", angles
			
	
		
	if len( intilts ) != len( orderedangles ):
		
		if not options.stackregardless:	
			print """\n(e2spt_tiltstacker.py)(organizetilts) ERROR: Number of tilt angles 
			and tilt images is not equal."""
			print """Number of tilt angles = %d ; number of tilt images = %d """ % ( len( orderedangles ), len( intilts ) )
			sys.exit()
		else:
			print """\n(e2spt_tiltstacker.py)(organizetilts) WARNING: Number of tilt angles 
			and tilt images is not equal. Stacking nevertheless since you provided --stackregardless""", options.stackregardless
			print """Number of tilt angles = %d ; number of tilt images = %d """ % ( len( orderedangles ), len( intilts ) )
			
	
	tiltstoexclude = options.exclude.split(',')
	
	indexesintiltseries = range(len(angles))
	collectionindexes = range(len(angles))
	
	print "options.bidirectional is", options.bidirectional
	print "indexminangle is", indexminangle

	if options.bidirectional and indexminangle != None:
		firstrange = range(0,indexminangle+1)
		secondrange = range(indexminangle+1,len(angles))			
		
		collectionindexes = []
		collectionindexes = list(firstrange) + list(secondrange)
		
		
		firstrange.sort()
		firstrange.reverse()
		
		indexesintiltseries = []
		indexesintiltseries = list(firstrange) + list(secondrange)
		
		
		
		
	print "collection indexes are", collectionindexes
	print "indexes in tiltseries are", indexesintiltseries
		
		
	for k in range(len(intilts)):
		
		try:
			tiltangle = orderedangles[k]
			#indexintiltseries = angles.index( orderedangles[k] )
		
			indexintiltseries = indexesintiltseries[k]
			collectionindex = collectionindexes[k]
		except:
			if options.stackregardless:
				tiltangle = orderedangles[k-1] + orderedangles[k-1] - orderedangles[k-2]
				indexintiltseries = indexesintiltseries[k-1] + indexesintiltseries[k-1] - indexesintiltseries[k-2]
				collectionindex = collectionindexes[k-1] + collectionindexes[k-1] - collectionindexes[k-2]
			
				
		if indexintiltseries not in tiltstoexclude and str(indexintiltseries) not in tiltstoexclude:
		
			intiltsdict.update( { indexintiltseries:[ intilts[k],tiltangle,collectionindex ]} )
			print "\nadded collectionIndex=%d, tiltangle=%f, indexintiltseries=%d" % ( collectionindex, tiltangle, indexintiltseries )
 				
	return intiltsdict


def usntacker( options ):
	
	print "e2spt_tiltstacker (unstacker). options.unstack is", options.unstack
	
	outname = options.path + '/' + options.unstack.replace('.mrc','.hdf')
	outname = options.path + '/' + options.unstack.replace('.mrcs','.hdf')
	outname = options.path + '/' + options.unstack.replace('.st','.hdf')	
	outname = options.path + '/' + options.unstack.replace('.ali','.hdf')
	
	outname = outname.replace('.hdf','_UNSTACKED.hdf')
	
	print "unstack outname is", outname
	
	#print "\noutname of unstacked tilt will be", outname
	
	cmdun = 'e2proc2d.py ' + options.unstack + ' ' + outname + ' --unstacking '
	if options.outmode:
		cmdun += ' --outmode=' + options.outmode
	
	if options.exclude or options.include:
		lst = makeimglist( options.unstack, options ) 
		cmdun += ' --list=' + lst
	
	print "\ncmdun is", cmdun	
	p = subprocess.Popen( cmdun , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text = p.communicate()	
	p.stdout.close()

	return outname
		

def restacker( options ):
	inputf = options.restack
	
	outname = options.path + '/' + options.restack.replace('.hdf','_RESTACKED.hdf')
	
	if '.ali' in inputf[-4:]: 
		outname = options.path + '/' + options.restack.replace('.ali','_RESTACKED.ali')
	if '.mrc' in inputf[-4:]: 
		outname = options.path + '/' + options.restack.replace('.mrc','_RESTACKED.mrc')
	if '.mrcs' in inputf[-5:]: 
		outname = options.path + '/' + options.restack.replace('.mrcs','_RESTACKED.mrcs')
	if '.st' in inputf[-3:]:
		outname = options.path + '/' + options.restack.replace('.st','_RESTACKED.st')
		
	if '.hdf' in inputf[-3:]:
		outname = options.path + '/' + options.restack.replace('.st','_RESTACKED.hdf')
	
	print "\n!!!!!!!!!!!!!!!!\n\nOutname is", outname
	
	tmp = options.path + '/' + 'tmp.hdf'
	
	if '.ali' in inputf[-4:] or '.mrc' in inputf[-4:] or '.mrcs' in inputf[-5:] or '.st' in inputf[-3:] or '.hdf' in inputf[-4:]  :
		tmp = options.path + '/' + 'tmp.mrcs'
	
	cmdre = 'e2proc2d.py ' + options.restack + ' ' + tmp
	
	if options.outmode:
		cmdre += ' --outmode=' + options.outmode
	
	if options.exclude or options.include:
		lst = makeimglist( options.restack, options )
		cmdre += ' --list=' + lst
	
	
		
	cmdre += ' && mv ' + tmp + ' ' + outname
	
	print "\ncmdre is", cmdre
	p = subprocess.Popen( cmdre , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text = p.communicate()	
	p.stdout.close()
	
	return


def getindxs( string ):
	
	parsedindxs = set([])
	stringList = list( string.split(',') )
	#print "\nstringList is", stringList
	
	intList=set([])
	for i in range( len(stringList) ):
		print "\n\nstringList[i] is", stringList[i]
		
		if '-' in stringList[i]:
			stringList[i] = stringList[i].split('-')
			
			x1 = int( stringList[i][0] )
			x2 = int( stringList[i][-1] ) + 1
			
			#stringList[i] = [ str( ele ) for ele in xrange( x1, x2 ) ]
			intList.union([ ele for ele in xrange( x1, x2 ) ])
			
		#parsedindxs = parsedindxs.union( set( list( List[i] ) ) )
		
		#intList = set( list(intList) )
	if intList:
		parsedindxs = [ str(ele) for ele in intList ]
	else:
		parsedindxs = stringList
		
	#parsedindxs = list( parsedindxs )
	parsedindxs.sort()
	
	print "\nParsed indexes are", parsedindxs
	return parsedindxs
	

def makeimglist( inputf, options ):
	
	n = EMUtil.get_image_count( inputf )
	if '.ali' in inputf[-4:] or '.mrc' in inputf[-4:] or '.mrcs' in inputf[-5:] or '.st' in inputf[-3:]:
		n = EMData( inputf,0,True )['nz']

	print "input n is", n
	print "input is",input
	
	allindxs = set( [ str(i) for i in range(n) ] )
	finalindxs = set( list( allindxs ) )
	
	print "\nallindxs are", allindxs
	
	if options.exclude:
		print "\nPrint there's EXCLUDE!!"
		eindxs = getindxs( options.exclude )
		finalindxs = list( allindxs.difference( eindxs ) )
		
	elif options.include:
		print "\nPrint there's INCLUDE!!"
		iindxs = getindxs( options.include )
		finalindxs = list( iindxs )
	
	print "\nFinalindxs are", finalindxs
	
	ints = []
	for id in finalindxs:
		ints.append( int( id ) )
	
	ints.sort()
	
	print "\nfinalindx sorted are", ints
	#final = []
	#for fi in ints:
	#	final.append( str( fi ) )
	
	lines = []
	for indx in ints:
		lines.append( str(indx) + '\n' )
		
	listfile = options.path + '/list.lst'
	f= open( listfile,'w' )
	f.writelines( lines )
	f.close()
	
	print "\nlistfile is", listfile
	return listfile


def floatrange(start, stop, step):
	print "\nInside floatrange, start, stop and step are", start, stop, step
	
	r = start
	kkk=0
	while r <= stop:
		yield r
		print "r is", r
		r += step
		kkk+=1
		
		if kkk > 180:
			print "ERROR: Something is wrong with your tiltrange, lowesttilt or highesttilt"
			sys.exit()
	
	return


if __name__ == "__main__":
	main()