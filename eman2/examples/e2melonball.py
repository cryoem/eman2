#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya 2014 (jgmontoy@bcm.edu); last update 06/14
# Copyright (c) 2000-2011 Baylor College of Medicine
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


from EMAN2 import *
import math
import os
import sys
from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program extracts parts of 3D maps from specific coordinates with the shape of a 
	user-defined mask, or from symmetry-related positions based on a defined radius and, 
	again, a user-defined mask.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", default='',type=str, help="""3D map or stack of maps to extract smaller regions from.""")
	
	parser.add_argument("--output", default='extracts',type=str, help="""String to use as the 'stem' for naming output volumes. Note that for each input volume, you'll get a stack of subvolumes. For example, if you provide a stack with 3 volumes, and you extract 12 subvolumes from each of these, you'll have 3 stacks of extracted subvolumes.""")
		
	parser.add_argument("--path",type=str,help=""""Name of directory where to store the output file(s)""",default="melonscoops")
	
	parser.add_argument("--sym", dest = "sym", default="c1", help = """Specify symmetry. Choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction ommit this option or specify c1.""")
	
	parser.add_argument("--vertices", action='store_true', default=False,help="""Only works if --sym=icos. This flag will make the program extract only the 12 vertices from among all 60 symmetry related units.""") 
	
	parser.add_argument("--mask",type=str,help="""Mask processor to define the shape of regions to extract. Default is None.""", default="")
		
	parser.add_argument("--maskfile",type=str,help="""Precomputed mask to use to extract subvolumes from the locations specified through --coords or through --radius and --sym""", default="")
		
	parser.add_argument("--savescoops",action='store_true',default='',help="""Save extracted parts from each particle into a per-particle stack, with extracted subvolumes centered in a box of the size specified by --boxsize, and rotated so that each subvolume is pointing along Z.""")
	
	parser.add_argument("--savescoopsinplace",action='store_true',default='',help="""Save extracted parts from each particle into a per-particle stack, with extracted subvolumes 'in-situ'; that is, with the same size and orientation as in the original volume.""")		
		
	parser.add_argument("--coords",type=str,help="""File with coordinates from where to extract subvolumes.""", default="")
		
	parser.add_argument("--radius",type=int,help="""Radius (in pixels) where to center the mask for subvolume extraction. Works only for cases in which the asymmetric unit of interest lies along z (for example, vertexes of an icosahedral virus aligned to the symmetry axes such that a vertex lies along z). Supplying --tz should achieve the same results.""", default=0)
	
	parser.add_argument("--tx",type=int,help="""Translation (in pixels) along x to define the mask's center. If supplied with --radius, the latter will be ignored.""", default=0)

	parser.add_argument("--ty",type=int,help="""Translation (in pixels) along y to define the masks's center. If supplied with --radius, the latter will be ignored.""", default=0)

	parser.add_argument("--tz",type=int,help="""Translation (in pixels) along z to define the masks's center. If supplied with --radius, the latter will be ignored.""", default=0)
	
	#parser.add_argument("--normproc",type=str,default='',help="Normalization processor applied to particles before alignment. Default is to use normalize. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'")
	#
	#parser.add_argument("--threshold",default='',type=str,help="""A threshold applied to the subvolumes after normalization. 
	#												For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#
	#parser.add_argument("--preprocess",default='',type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
		
	#parser.add_argument("--lowpass",type=str,default='',help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	#parser.add_argument("--highpass",type=str,default='',help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--boxsize",type=int,default=0,help="""If specified, the output subvolumes will be clipped to this size and centered in the box; otherwise, they will be saved in a boxsize equal to the volume they came from.""")
	
	#parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="align")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")	

	(options, args) = parser.parse_args()
	
	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if not options.input:
		parser.print_help()
		sys.exit(0)
	
	
	if options.savescoops and not options.boxsize:
		print """ERROR: Specify the box size through --boxsize for the saved
		scoops."""	
		sys.exit()

	#If no failures up until now, initialize logger
	logid=E2init(sys.argv,options.ppid)
	
	inputhdr = EMData(options.input,0,True)
	
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,options.path)
	
	
	if not options.mask and not options.maskfile:
		print "\nERROR: You must define --mask or supply a volume through --maskfile"
		sys.exit()
	
	mask = EMData(inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])
	mask.to_one()
	
	if options.mask:	
		mask.process_inplace(options.mask[0],options.mask[1])
	
	elif options.maskfile:
		mask = EMData( options.maskfile )
		
		'''
		Center of maskfile
		'''
		maskcx = mask['nx']
		maskcy = mask['ny']
		maskcz = mask['nz']
		
		'''
		Clip to data's boxsize
		'''
		r=Region( (2*maskcx - inputhdr['nx'])/2, (2*maskcx - inputhdr['ny'])/2, (2*maskcx - inputhdr['nz'])/2, inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])
		mask.clip_inplace( r )
	
	print "\nMask done"
	
	if options.radius:
		mask.translate(0,0, float(options.radius) )
		print "\nMask translated by radius %d along z" %( options.radius )

	else:
		mask.translate( options.tx, options.ty, options.tz )
		print "\nMask translated by tx=%d, ty=%d, tz=%d,", options.tx, options.ty, options.tz
	
	symnames = ['oct','OCT','icos','ICOS','tet','TET']
	
	orientations = []
	#transforms = []
	anglelines = []
	
	if options.sym:
		print "\nsym found", options.sym
	
		symnum = 0
		
		if options.sym not in symnames:
			
			symletter = options.sym
			symnum = options.sym
			
			for x in options.sym:
				if x.isalpha():
					symnum = symnum.replace(x,'')
					
			for x in options.sym:
				if x.isdigit():
					symletter = symletter.replace(x,'')
			
			
			print "\nThe letter for sym is", symletter
			print "\nThe num for sym is", symnum
		
		if options.sym == 'oct' or options.sym == 'OCT':
			symnum = 8
			
		if options.sym == 'tet' or options.sym == 'TET':
			symnum = 12	
			
		if options.sym == 'icos' or options.sym == 'ICOS':
			symnum = 60	
		
		symnum = int(symnum)
			
		t = Transform()
		
		if symnum:
			print "\nsymnum determined",symnum
			print "while symletter is", symletter
			
			if symletter == 'd' or symletter == 'D':
				symnum *= 2
				print "\nsymnum corrected, because symmetry is d",symnum
				
			#if options.save
			
			for i in range(symnum):
				orientation = t.get_sym( options.sym , i )
				orientations.append( orientation )
					
				rot = orientation.get_rotation()
				az = rot['az']
				alt = rot['alt']
				phi = rot['phi']
				
				anglesline = 'az=' + str(az) + 'alt=' + str(alt) + 'phi=' + str(phi) + '\n'
				anglelines.append(anglesline)
				
				#transforms.append( Transform({'type':'eman','az':az,'alt':alt,'phi':phi}) )
					
	n = EMUtil.get_image_count( options.input )
	
	centers = {}
	masks = {}
	
	if orientations:
		print "\ngenerated these many orientations", len (orientations)
		if options.sym == 'icos' or options.sym == 'ICOS':
			if options.vertices:

				print "\nbut fetching vertices only"
			
				#orientations = genicosvertices ( orientations )

				orientations = genicosverticesnew ( orientations )
			
			
		centerx = 0.0
		centery = 0.0
		centerz = float( options.radius )
		
		if not options.radius:
			centerx = options.tx 
			centery = options.ty 
			centerz = options.tz
		
		groundv = Vec3f(centerx,centery,centerz)
		
		finalmask = EMData(mask['nx'],mask['ny'],mask['nz'])
		finalmask.to_zero()
		
		centersmap = EMData( mask['nx'],mask['ny'],mask['nz'] )
		centersmap.to_zero()
		
		centerstxt = options.path + '/centers.txt'
		f = open( centerstxt, 'a')
		for k in range( len(orientations) ):
			t =  orientations[k]
			
			print "\nWorking with this orientation",t
			
			tmpmask = mask.copy()
			tmpmask.transform( t )
			
			tmpmask['origin_x'] = 0
			tmpmask['origin_y'] = 0
			tmpmask['origin_z'] = 0
			
			finalmask = finalmask + tmpmask
			
			center = t*groundv
			centers.update( {k: center} )
			
			newcenterx = int( round(center[0] + mask['nx']/2.0 )) 
			newcentery = int( round(center[1] + mask['ny']/2.0 ))
			newcenterz = int( round(center[2] + mask['nz']/2.0 ))
			
			tout = Transform()
			rot = t.get_rotation()
			
			tout.set_rotation( rot )
			tout.set_trans( newcenterx, newcentery, newcenterz )
			
			tmpmask['spt_scoop_center'] = [ newcenterx, newcentery, newcenterz ]
			tmpmask['xform.align3d'] = tout
			
			tmpmask.write_image(options.path + '/masks.hdf',k)
			
			masks.update( {k:[tmpmask,tout,[newcenterx,newcentery,newcenterz]]} )
			
			
			centerline = str(newcenterx) + ' ' + str(newcentery) + ' ' + str(newcenterz) + '\n'
			
			f.write(centerline)
			
			centersmap.set_value_at( int( round(center[0] + mask['nx']/2.0 )) , int( round(center[1] + mask['ny']/2.0 )), int( round(center[2] + mask['nz']/2.0 )), 1.0 )
		
			
			
			#print "setting indx and t are", k, t
			#print "whereas mask t is", tout
		
		f.close()
			
		finalmask.write_image(options.path + '/finalmask.hdf',0)
		centersmap.write_image(options.path + '/centersmap.hdf',0)
		
		for i in range(n):
			ptcl = EMData(options.input,i)
			
			ptclnumtag = str(i).zfill(len( str(n) ))
			
			scoopsinplacestack = 'ptcl' + ptclnumtag + '_scoopsOrigBox.hdf'
			
			scoopsstack = 'ptcl' + ptclnumtag + '_scoops.hdf'

			ptclglobalmsk = ptcl.copy()
			ptclglobalmsk.mult( finalmask )
			
			ptclglobalmsk.write_image( options.path + '/' + 'ptcls_globallymasked.hdf',i)
			
			for key in masks.keys():
			
				scoopinplace = ptcl.copy()
				scoop = ptcl.copy()
				
				thismask = masks[key][0]
				
				t = masks[key][1]
				
				#tapp = Transform()
				#rot = t.get_rotation()
				#tapp.set_rotation( rot )
				
				#print "key and t are", key, tapp
				sx = masks[key][-1][0]
				sy = masks[key][-1][1]
				sz = masks[key][-1][2]
				
				if options.savescoopsinplace:
					
					#scoopinplace.mult( thismask )
					scoopinplace.process_inplace('normalize.edgemean')
					scoopinplace.mult( thismask )
					
					
					scoopinplace['origin_x'] = 0
					scoopinplace['origin_y'] = 0
					scoopinplace['origin_z'] = 0
					
					scoopinplace['spt_scoop_x'] = sx
					scoopinplace['spt_scoop_y'] = sy
					scoopinplace['spt_scoop_z'] = sz
					
					scoopinplace['spt_score'] = 0
					
					scoopinplace['xform.align3d'] = t
					
					scoopinplace.write_image( options.path + '/' + scoopsinplacestack, key )
					print "\nWrote this scoop 'in place' ", key
					
				if options.savescoops:
					
					scoop.mult( thismask )
					
					box = options.boxsize
		
					extra = math.fabs( math.sqrt(2.0) * ( box/2.0 - 1.0 ) - ( box/2.0 - 1 ) )	
					paddedbox = int( math.ceil( float(box) + extra ) )
					
					#print "Center for region is", sx,sy,sz
						
					print "\nExtracting this scoop", key
					print "With a padded box of this size", paddedbox
					bigr = Region( (sx*2 - paddedbox)/2 ,  (sy*2 - paddedbox)/2 ,  (sz*2 - paddedbox)/2, paddedbox,paddedbox,paddedbox)
					bigscoop = scoop.get_clip(bigr)
					
					#bigscoop.write_image( options.path + '/bigscoops.hdf', key)
					
					bigscoop['origin_x'] = 0
					bigscoop['origin_y'] = 0
					bigscoop['origin_z'] = 0
					
					print "\nOrienting the subscoop with this transform", t
					
					t2use = Transform()
					rot = t.get_rotation()
					
					t2use.set_rotation( rot )
					
					ti = t2use.inverse()
					bigscoop.transform( ti )
					
					#bigscoop.write_image( options.path + '/bigscoopsOriented.hdf', key)
					
					sxB = bigscoop['nx']/2
					syB = bigscoop['ny']/2
					szB = bigscoop['nz']/2
					
					print "\nCenter of bigs scoops is", sxB,syB,szB
					
					scoop = clip3D( bigscoop, box )
					
					#r = Region( (sxB*2 - box)/2 ,  (syB*2 - box)/2 ,  (szB*2 - box)/2, box,box,box)
					#scoop = bigscoop.get_clip(r)
					
					#print "\nTherefore region for small scoop is", r
					
					print "\nClipping the extracted and oriented scoop back to the desired boxsize", box
					
					defaultmask = EMData( scoop['nx'], scoop['ny'],scoop['nz'])
					defaultmask.to_one()
					defaultmask.process_inplace('mask.sharp',{'outer_radius':-1})
					
					#scoop.mult(defaultmask)
					#scoop.process_inplace('normalize.edgemean')
					scoop.mult(defaultmask)
					
					scoop['origin_x'] = 0
					scoop['origin_y'] = 0
					scoop['origin_z'] = 0
					
					scoop['spt_score'] = 0
					
					trans = t.get_trans()
					transi = -1 * trans
					
					ti2write = Transform()
					ti2write.set_rotation( ti.get_rotation() )
	 
					ti2write.set_trans( transi )
					
					scoop['xform.align3d'] = ti2write
					
					
					
					scoop['spt_scoop_center'] = [sx, sy, sz]
					
					scoop.write_image( options.path + '/' + scoopsstack, key)	
		
	E2end(logid)
	
	return



def clip3D( vol, size ):
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - size)/2, (2*volyc - size)/2, (2*volzc - size)/2, size , size , size)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol


def genicosverticesnew( syms ):

	newsyms = []
	for s in syms:
		rot = s.get_rotation()
		
		if float(rot['alt']) > 63.0 and float(rot['alt']) < 64.0 and int(round(float(rot['phi']))) == 90:
			
			t = Transform( {'type':'eman','alt':rot['alt'], 'phi':rot['az'], 'az':rot['phi']} )
			newsyms.append( t )
		
		if float(rot['alt']) > 116.0 and float(rot['alt']) < 117.0 and int(round(float(rot['phi']))) ==126:
			t = Transform( {'type':'eman','alt':rot['alt'], 'phi':rot['az'], 'az':rot['phi']} )
			newsyms.append( t )
	
	newsyms.append( Transform({'type':'eman','alt':0.0, 'phi':0.0, 'az':0.0}) )
	
	newsyms.append( Transform({'type':'eman','alt':180.0, 'phi':180.0, 'az':0.0}) )
	
	return( newsyms )


'''
def genicosvertices( syms ):
	t=Transform()

	syms = []
	for i in range(60):
		syms.append(t.get_sym('icos',i))
	
	ts = []
	azs = set([])
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
'''

if __name__ == '__main__':
	main()
	