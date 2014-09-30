#!/usr/bin/python2.7

#====================
#Author: Jesus Galaz-Montoya 22/sep/2014 , Last update: sep/22/2014
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

def main():

	usage = """e2spt_stackdif.py <options> . 
	This program has 2 functions: 
	A) To determine whether there are overlaps between two stacks of subtomograms extracted
	using e2spt_boxer.py. In this case, supply --stack1 and --stack2.
	The output from this will be written to 'commonparticles.txt'
	
	B) To determine the number of different tomograms
	the particles in a stack come from, how many particles come from each
	tomgoram.
	In this case, you can supply --stack1 only, or both --stack1 and --stack2.
	The output from this will be written to stack1_STATS.txt and/or stack2_STATS.txt.	
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	#parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
	#	The default is a numbered series of directories containing the prefix 'spt_stackdif';
	#	for example, 'spt_stackdif02' will be the directory by default if 'spt_stackdif01' 
	#	already exists.""")

	parser.add_argument("--stack1", type=str, default='', help="""Subtomogram stack in HDF format.""")
	
	parser.add_argument("--stack2",type=str,default='',help=""""The other subtomogram stack in HDF format.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)
	
	#from e2spt_classaverage import sptmakepath
	#options = sptmakepath( options, 'spt_stackdif')
	
	stack1 = options.stack1 
	coords1set = analyzefile( stack1, 1 )
	
	if options.stack2:
		stack2 = options.stack2
		coords2set = analyzefile( stack2, 2 ) 
	
		common = coords1set.intersection( coords2set )
	
		if common:
	
			print "\nThere are these many particles common to both stacks", len( common )
			print "Which are"
			
			g = open('commonparticles.txt','w')
			lines = []
			for val in common:
				line = str(val) + '\n'
				lines.append( line )
				print val
			g.writelines( lines )
			g.close()
	
	E2end( logger )
	return


def analyzefile( stack, stacknum ):
	
	n = EMUtil.get_image_count( stack )

	coordslist = []
	tomograms = set([])
		
	for i in range( n ):
		hdr = EMData( stack, i, True )
		tomogram = hdr['ptcl_source_image']
		coords = hdr['ptcl_source_coord']
		coordslist.append( (tomogram,tuple(coords)) )
		
		tomograms.add( tomogram )
	
	tomodict = {}
	for tomo in tomograms:
		tomodict.update( { tomo:0 } )
		
	for i in range( n ):
		hdr = EMData( stack, i, True )
		tomogram = hdr['ptcl_source_image']
		
		count = tomodict[tomogram] + 1
		
		tomodict.update( { tomogram:count } )
		
	stackstatisticsfile = 'stack' + str(stacknum) + '_STATS.txt'
	f=open(stackstatisticsfile,'w')
	lines = []
	for tomo in tomodict:
		line = 'Tomogram: ' + tomo + '\tnumber of particles: ' + str( tomodict[tomo] ) + '\n'
		lines.append( line )
	f.writelines( lines )
	f.close()
	
	#print "\nCoordslist type", type(coordslist)
	#print "coordslist is", coordslist
	
	coordsset = set( coordslist )
	
	if len( coordsset ) != len( coordslist ):
		print "\nThere seem to be repeated particles within stack %d" % ( stacknum )
		print "There are these many particles", len( coordslist )
		print "But only these many are unique", len( coordsset )
	
	#print "\n\n\n\n\n\nCoordsset type is", type(coordsset)
	return coordsset
	
	
		
if __name__ == '__main__':
	
	main()