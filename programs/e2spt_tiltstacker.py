#!/usr/bin/env python

'''
====================
Author: Jesus Galaz-Montoya 2/20/2013 , Last update: December/12/2013
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

	parser.add_argument("--lowerend",type=float,default=None,help="Lowest tilt angle.")
	parser.add_argument("--upperend",type=float,default=None,help="Highest tilt angle.")
	parser.add_argument("--tiltstep",type=float,default=None,help="Step between tilts.")
	parser.add_argument("--apix",type=float,default=None,help="True apix of images to be written on final stack.")

	#parser.add_argument("--onlyy",action='store_true',default=False,help="Only projection of the XZ plane will be generated [another 'side view'].")
	
	#parser.add_argument("--stack",action='store_false',default=True,help="If on, projections will be in an hdf stack; otherwise, they'll be their own separate file. On by default. Supply --stack=None to turn off.")

	parser.add_argument("--output", type=str, help="""File name to store the stacked tiltseries.""", default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)
	
	currentdir = os.getcwd()
	
	filesindir = os.listdir(currentdir)
	
	if not options.output:
		print "ERROR: Please provide output file"
		sys.exit()
	
	outtilts=[]
	for f in filesindir:
		outtilt=f
		if '.dm3' in f or '.DM3' in f or '.tif' in f or '.TIF' in f or '.hdf' in f or '.HDF' in f:
			outtilt=outtilt.replace('.dm3','.mrc')
			outtilt=outtilt.replace('.tif','.mrc')
			outtilt=outtilt.replace('.hdf','.mrc')
			outtilt=outtilt.replace('.DM3','.mrc')
			outtilt=outtilt.replace('.TIF','.mrc')
			outtilt=outtilt.replace('.HDF','.mrc')
				
			os.system('e2proc2d.py ' + f + ' ' + outtilt + ' --mrc16bit')
		elif '.mrc' in f:	
			outtilts.append(outtilt)
			
	outtilts.sort()
	
	print "The sorted tilts are", outtilts
	
	k=0
	for t in outtilts:
		print "t is",t
		a=EMData(t,0)
		if options.apix:
			a['apix_x'] = options.apix
			a['apix_y'] = options.apix
			a['apix_z'] = options.apix

		a.write_image('tmp.hdf',k)
		k+=1
	
	os.system('e2proc2d.py tmp.hdf ' + options.output + ' --twod2threed && mv ' + options.output + ' ' + options.output.replace('.mrc','.st'))
	
	if options.lowerend and options.upperend and options.tiltstep:
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

			
