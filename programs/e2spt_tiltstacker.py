#!/usr/bin/env python

'''
====================
Author: Jesus Galaz-Montoya 2/20/2013 , Last update: March/7/2014
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

from optparse import OptionParser
from EMAN2 import *
import sys

def main():

	usage = """e2spt_tiltstacker.py <options> . The options should be supplied in "--option=value", 
	replacing "option" for a valid option name, and "value" for an acceptable value for that option. 
	This program stacks individual .dm3, .tiff or .hdf images into an .mrc (or .st) stack. 
	It must be run in a directory containing the numbered images only.
	It also generates a .rawtlt file with tilt angle values if --lowerend, --upperend and --tiltstep are provided
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	#parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'orthoproject';
	#													for example, orthoproject_02 will be the directory by default if 'orthoproject_01' already exists.""")


	parser.add_argument("--input", type=str, help="""HDF file to unstack (must also specify 
		--unstack) or string common to all the files to put into an .st (MRC) stack; for example,
		'.hdf' will process all .hdf files.
		If not specified, all valid EM imagefiles in the current directory will be put into 
		an .st (MRC) stack.""", default='')
	parser.add_argument("--output", type=str, help="""File name to store the stacked tiltseries, or common string to save unstacked individual images if --unstack is provided.""", default='')


	parser.add_argument("--lowerend",type=float,default=0.0,help="Lowest tilt angle.")
	parser.add_argument("--upperend",type=float,default=0.0,help="Highest tilt angle.")
	parser.add_argument("--tiltstep",type=float,default=0.0,help="Step between tilts.")
	parser.add_argument("--apix",type=float,default=0.0,help="True apix of images to be written on final stack.")
	parser.add_argument("--unstack",type=str,default='',help="Name of the hdf file to unstack into individual images.")
	parser.add_argument("--mirroraxis",type=str,default='',help="""Options are x or y, and the
		a mirrored copy of the 2-D images will be generated before being put into the tilt series.""")
	
	#parser.add_argument("--stack",action='store_false',default=True,help="If on, projections will be in an hdf stack; otherwise, they'll be their own separate file. On by default. Supply --stack=None to turn off.")


	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)
	
	currentdir = os.getcwd()
	
	filesindir = os.listdir(currentdir)
	
	if not options.output:
		print "ERROR: Please provide output file"
		sys.exit()
	
	mrcs=[]
	mrcsmirror=[]
	hdfs=[]
	hdfsmirror=[]
	
	if options.unstack:
		n=EMUtil.get_image_count(options.unstack)
		print "\nThe number of images to unstack is", n
		for i in range(n):
			outname = options.output.split('.')[0] + '_' + str(i).zfill( len( str(n))) + '.hdf'
			print "Outname of unstacked tilt will be", outname
			cmd = 'e2proc2d.py ' + options.unstack + ' ' + outname + ' --first=' + str(i) + ' --last=' + str(i) + ' && e2proc2d.py ' + outname + ' ' + outname.replace('.hdf','.mrc') + ' --mrc16bit'
			p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()
			print "Unstacked image", i	
	
	else:
		for f in filesindir:
			fyle=''
			if not options.input:
				fyle=f
			elif options.input:
				if options.input in f:
					fyle=f
		
			if fyle:	
				outtilt=fyle
				#Convert from other formats to HDF
				if '.dm3' in f or '.DM3' in f or '.tif' in f or '.TIF' in f or '.MRC' in f and '.txt' not in f and '.db' not in f and 'mirror' not in f:
					print "\nProcessing file", f
					outtilt=outtilt.replace('.dm3','.hdf')
					outtilt=outtilt.replace('.tif','.hdf')
					outtilt=outtilt.replace('.DM3','.hdf')
					outtilt=outtilt.replace('.TIF','.hdf')
					outtilt=outtilt.replace('.HDF','.hdf')
			
					cmd = 'e2proc2d.py ' + fyle + ' ' + outtilt #+ ' --mrc16bit' # --fixintscaling=sane'
					p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text = p.communicate()	
					p.stdout.close()
						
				if options.apix:
					print "Fixing apix"
					cmdapix = 'e2fixheaderparam.py --input=' + outtilt + ' --stem=apix --valtype=float --stemval=' + str( options.apix )
					p = subprocess.Popen( cmdapix , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text = p.communicate()	
					p.stdout.close()
		
				hdfs.append( outtilt )
		
				if options.mirroraxis:
					print "Mirroring"
					outtiltmirror = outtilt.replace('.hdf','_mirrorY.hdf')
		
					cmdMirror = 'e2proc2d.py ' + outtilt + ' ' + outtiltmirror + ' --process=xform.mirror:axis=' + options.mirroraxis #+ ' --mrc16bit' #+ ' --fixintscaling=sane'
					p = subprocess.Popen( cmdMirror , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text = p.communicate()	
					p.stdout.close()
			
					hdfsmirror.append( outtiltmirror )
		
		
				print "Converting to mrc"
				outtiltmrc = outtilt.replace('.hdf','.mrc')
			
				#os.system('e2proc2d.py ' + f + ' ' + outtilt + ' --mrc16bit')
		
				cmdmrc = 'e2proc2d.py ' + outtilt + ' ' + outtiltmrc + ' --mrc16bit' + ' --fixintscaling=sane'
		
				p = subprocess.Popen( cmdmrc , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text = p.communicate()	
				p.stdout.close()
		
				mrcs.append( outtiltmrc )
		
		
				if options.mirroraxis:
					print "Converting mirror to mrc"
					outtiltmirrormrc = outtilt.replace('.mrc','_mirrorY.mrc')
					cmdMirrormrc = 'e2proc2d.py ' + outtiltmrc + ' ' + outtiltmirrormrc + ' --process=xform.mirror:axis=' + options.mirroraxis + ' --mrc16bit' + ' --fixintscaling=sane'
		
					p = subprocess.Popen( cmdMirrormrc , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text = p.communicate()	
					p.stdout.close()
			
					mrcsmirror.append( outtiltmirrormrc )
			else:
				pass
				
				
		print "Sorting stacks"
		mrcs.sort()	
		hdfs.sort()
	
		if options.mirroraxis:
			mrcsmirror.sort()
			hdfsmirror.sort()
	
		print "\nI'll write temporary hdf stack!!!!!!!!!!!!!!!"

		for k in range(len(hdfs)):
		
			#cmd = 'e2proc2d.py ' + hdfs[k] + ' tmp.hdf --append'
			#p = subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			#text = p.communicate()	
			#p.stdout.close()
		
			print "Inserting image into temporary stack",text

			a=EMData(hdfs[k])
			print "A is", a
			print "A is type", type(a)
			a.write_image('tmp.hdf',-1)
		
		
			#print "cmd is", cmd
			#a=EMData(hdfs[k])	
			#a.write_image('tmp.hdf',k)
		
			if options.mirroraxis:
				#print "\nWriting temporary mirror stack"
				#cmdmirror = 'e2proc2d.py ' + hdfsmirror[k] + ' tmpmirror.hdf --append'
			
				#print "cmdmirror is",cmdmirror
				#p = subprocess.Popen( cmdmirror, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				#text = p.communicate()	
				#p.stdout.close()
		
				print "Inserting image into temporary mirror stack",text
				a=EMData(hdfsmirror[k])
				a.write_image('tmpmirror.hdf',-1)
			
				#hdfsmirror[k].write_image('tmpmirror.hdf',k)
	
	
		mrcout = options.output.split('.')[0] + '.mrc'
		stout = options.output.split('.')[0] + '.st'
	
	
		print "Converting 2-D hdf stack to 3-D mrc stack"

		cmdst = 'e2proc2d.py tmp.hdf ' + mrcout + ' --twod2threed' + ' --mrc16bit' + ' --fixintscaling=sane'
		cmdst += ' && mv ' + mrcout + ' ' + stout + ' && rm tmp.hdf'
	
		print "cmdst is",cmdst
	
		p = subprocess.Popen( cmdst , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text = p.communicate()	
		p.stdout.close()
		
	
		print "Done", text
	
		if options.mirroraxis:
			print "Converting 2-D hdf mirror stack to 3-D mrc mirror stack"

			mrcoutmirror = options.output.split('.')[0] + '_mirror.mrc'
			stoutmirror = options.output.split('.')[0] + '_mirror.st'

			cmdstmirror = 'e2proc2d.py tmpmirror.hdf ' + mrcoutmirror + ' --twod2threed' + ' --mrc16bit' + ' --fixintscaling=sane' + ' && mv ' + mrcoutmirror + ' ' + stoutmirror + ' && rm tmpmirror.hdf' 
			print "cmdstmirror is", cmdstmirror
		
			p = subprocess.Popen( cmdstmirror , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()
		
			print "Done",text
	
		if options.lowerend and options.upperend and options.tiltstep:
			print "Generating .rawtlt file"
			tltfile = options.output.split('.')[0] + '.rawtlt'
			f = open(tltfile,'w')

			generate=floatrange(options.lowerend, options.upperend, options.tiltstep)
			tilts=["%g\n" % x for x in generate]
				
			f.writelines(tilts)
			f.close()
	
	return

def floatrange(start, stop, step):
	r = start
	while r <= stop:
		yield r
		r += step
		

if __name__ == "__main__":
	main()