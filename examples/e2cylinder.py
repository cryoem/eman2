#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - 02/March/2013, Last update: 29/October/2014
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

from EMAN2 import *


def main():

	usage = """Program to generate a cylindrical mask. It can also create a cylindrical shell if
			you specify the --height_inner and --radius_inner parameters, in addition to the
			required outer --height and --radius.
			"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'cylmask';
		for example, cylmask_02 will be the directory by default if 'cylmask_01' 
		already exists.""")

	parser.add_argument("--verbose", "-v", help="""verbose level [0-9], higner number means 
		higher level of verboseness. Default=0.""",dest="verbose", action="store", metavar="n",type=int, default=0)

	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, 
		used for cross platform PPID""",default=-1)

	parser.add_argument("--height", type=int, default=0,help="""Height of the cylindrical mask.""")
	parser.add_argument("--heightinner", type=int, default=0,help="""Height for the inner boundary
		if creating a cylindrical shell mask.""")

	parser.add_argument("--radius",type=int,default=0,help="""Radius of the cylindrical mask.
		Default=boxsize/2.""")
	parser.add_argument("--radiusinner",type=int,default=0,help="""Radius for the inner boundary 
		if creating a cylindrical shell mask.
		Default=boxsize/2.""")

	parser.add_argument("--boxsize", type=int, default=0,help="""Size of the boxsize where the
		cylindrical mask will live.""")
	
	parser.add_argument("--axes",type=str,default='z',help="""Axes along which the mask will be
		oriented. Default=z. You can supply more than one, separated with commas. For example:
		--axes=x,y,z.""")
		
	parser.add_argument("--rotation",type=str,default='',help="""Three comma separated Euler angles
		 az,alt,phi, to rotate the masks by before writing them out.""")
		
	parser.add_argument("--translation",type=str,default='',help="""Three comma separated coordinates
		x,y,z, to translate the masks by before writing them out.""")
		
	parser.add_argument("--rotavg",action='store_true',default=False,help="""This will compute the rotational average of the mask(s) in addition to writing the cylindrical mask itself out.""")
	
	(options, args) = parser.parse_args()	
	
	
	if not options.boxsize:
		print "You must provide --boxsize > 4"
		sys.exit(1)
	elif options.boxsize < 5:
		print "You must provide --boxsize > 4"
		sys.exit(1)		
	
	if options.heightinner and not options.radiusinner:
		print "If specifying --heightinner, you must also specify --radiusinner."
		sys.exit(1)	
	if options.radiusinner and not options.heightinner:
		print "If specifying --radiusinner, you must also specify --heightinner."
		sys.exit(1)	
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options, 'cylmask')
	
	logger = E2init(sys.argv, options.ppid)
	
	axes = options.axes.split(',')
	
	print "After splitting, axes=", axes
	
	#axisdict ={}
	#for axis in axes:
	#	axisdict.update( { 'z } )
	ts = {}
	
	mask = cylinder(options)
	
	rt=Transform()
	if options.rotation or options.translation:
		az=alt=phi=xc=yc=zc=0
		if options.rotation:
			angles=options.rotation.split(',')
			az=float(angles[0])
			alt=float(angles[1])
			phi=float(angles[2])
		if options.translation:
			trans=options.translation.split(',')
			xc=float(trans[0])
			yc=float(trans[1])
			zc=float(trans[2])
		rt=Transform({'type':'eman','az':az,'alt':alt,'phi':phi,'tx':xc,'ty':yc,'tz':zc})
		mask.transform( rt )
		
	
	for axis in axes:
		print "axis is", axis
		if 'z' in axis or 'Z' in axis:
			tz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
			print "added z transform"
			ts.update({'z':tz})
		if 'x' in axis or 'X' in axis:
			
			tx = Transform({'type':'eman','az':0,'alt':90,'phi':90})
			ts.update({'x':tx})
			print "added x transform"
		
		if 'y' in axis or 'Y' in axis:
			ty = Transform({'type':'eman','az':0,'alt':90,'phi':0})
			ts.update({'y':ty})
			print "added y transform"
	
	masknamebase = options.path + '/cylmask.hdf'
	
	for a in ts:
		maskt = mask.copy()
		tag = 'R'+str( options.radius ).zfill( len( str( options.radius))) + 'H'+str( options.height ).zfill( len( str( options.radius)))
		if options.radiusinner:
			tag+='RI'+str( options.radiusinner ).zfill( len( str( options.radius))) 
		if options.heightinner:
			 'HI'+str( options.height ).zfill( len( str( options.radius)))
		
		if a == 'z':
			maskz=mask.copy()
			maskname=masknamebase.replace('.hdf','_Z_') + tag + '.hdf'
			maskz.transform( ts[a] )
			maskz.write_image( maskname, 0 )
			
			if options.rotavg:
				rotavgname = maskname.replace('.','_ROTAVG.')
				maskzrotavg = maskz.rotavg_i()
				maskzrotavg.write_image( rotavgname , 0 )
				
		if a == 'x':
			maskx=mask.copy()
			maskx.transform( ts[a] )
			maskname=masknamebase.replace('.hdf','_X_') + tag + '.hdf'
			maskx.write_image( maskname, 0 )
			
			if options.rotavg:
				rotavgname = maskname.replace('.','_ROTAVG.')
				maskxrotavg = maskx.rotavg_i()
				maskxrotavg.write_image( rotavgname , 0 )
		
		if a == 'y':
			masky=mask.copy()
			masky.transform( ts[a] )
			maskname=masknamebase.replace('.hdf','_Y_') + tag + '.hdf'	
			masky.write_image( maskname, 0 )

			if options.rotavg:
				rotavgname = maskname.replace('.','_ROTAVG.')
				maskyrotavg = masky.rotavg_i()
		
				maskyrotavg.write_image( rotavgname , 0 )
	
	E2end(logger)
	return


def cylinder( options ):
	
	box = options.boxsize
	mask = EMData( box, box, box)
	mask.to_one()
	
	if options.radius:
		radius = options.radius
	else:
		radius = box/2.0
		
	if options.height:
		height = options.height
	else:
		height = box/2.0
	
	maskout = mask.process("testimage.cylinder",{'height':height,'radius':radius})
	finalmask = maskout
	if options.heightinner and options.radiusinner:
		maskin = mask.process("testimage.cylinder",{'height':options.heightinner,'radius':options.radiusinner})
		finalmask = maskout - maskin
		
	return finalmask
	
if __name__ == '__main__':
	main()
