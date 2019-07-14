#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

# Author: Jesus Galaz-Montoya 2014 (jgalaz@gmail.com); last update Nov/17
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

from past.utils import old_div
from builtins import range
from EMAN2_utils import *
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

	parser.add_argument("--boxsize", type=int, default=0, help="""If specified, the output subvolumes will be clipped to this size and centered in the box; otherwise, they will be saved in a boxsize equal to the volume they came from.""")

	parser.add_argument("--coords", type=str, default=None, help="""File with coordinates from where to extract subvolumes.""")

	parser.add_argument("--input", default=None, type=str, help="""3D map or stack of maps to extract smaller regions from. Presumably, this particles are aligned already to a known orientation or symmetry axis.""")

	parser.add_argument("--inputisraw", action='store_true', default=False, help="""requires --transformfile. 3D map or stack of maps to extract smaller regions from. These are raw particles that haven't been aligned; requires --transformfile so that the mask can be 'inversely' rotated to the frame of the raw particles.""")

	parser.add_argument("--mask", type=str, default=None, help="""default=none. Mask processor to define the shape of regions to extract. e.g., --mask=mask.sharp:outer_radius=10 would make a solid sphere of 10 pixels in diameter.""")
	parser.add_argument("--maskcenter", type=str, default=None, help="""default=None. Comma separated x,y,z values for the translation (in pixels) to define the center of --mask or --maskfile with respect to the middle of the particle boxes. This only needs to be provided for one asymmetric unit. For example, to extract vertexes from an icosahedrally symmetric virus particle, if the vertexes are at a distance of 10 pixels from the center of the virus, --maskcenter 0,0,10 would be the center of the vertex pointing in Z directly out of the screen.""")
	parser.add_argument("--maskfile", type=str, default=None, help="""Precomputed mask to use to extract subvolumes from the locations specified through --coords or through --radius and --sym""")
	
	parser.add_argument("--output", default='extracts', type=str, help="""String to use as the 'stem' for naming output volumes. Note that for each input volume, you'll get a stack of subvolumes. For example, if you provide a stack with 3 volumes (subtomograms), and you extract 12 subvolumes (sub-subtomograms) from each of these, you'll have 3 stacks of extracted subvolumes (sub-subtomograms); e.g., the virus vertexes for each of 3 virus particles.""")
		
	parser.add_argument("--path", type=str, default="melonscoops", help=""""default=melonscoops. Name of directory where to store the output file(s).""")
	parser.add_argument("--plotcenters", action='store_true', default=False, help="""default=flase. Requires --maskcenter. Generates a map with a bright dot at the center of the asymmetric units from which 'scoops' were extracted.""")
	
	#parser.add_argument("--radius", type=int, default=0, help="""default=0. Radius (in pixels) where to center the mask for subvolume extraction. Works only for cases in which the asymmetric unit of interest lies along z (for example, vertexes of an icosahedral virus aligned to the symmetry axes such that a vertex lies along z). Supplying --tz should achieve the same results.""")

	parser.add_argument("--scoopsinsitu", action='store_true', default=False, help="""default=false. Save extracted parts from each particle into a box of the original size of the particle.""")			
	parser.add_argument("--scoopsorient", action='store_true', default=False, help="""default=false. Save extracted parts from each particle into boxes of the size specified by --boxsize with the scoops rotated so that each subvolume is pointing along Z.""")			
	parser.add_argument("--scoopsperparticle", action='store_true', default=False, help="""default=false. Save extracted parts from each particle into a per-particle stack, in a box of the size specified by --boxsize.""")
	
	parser.add_argument("--sym", dest = "sym", default="c1", help = """default=c1 (no symmetry). Choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction ommit this option.""")
	
	parser.add_argument("--transformfile", type=str, default=None, help="""default=None. .json file with alignment parameters produced by e2spt programs""") 

	#parser.add_argument("--tx", type=int, default=0, help="""Translation (in pixels) along x to define the --mask's center.""")
	#parser.add_argument("--ty", type=int, default=0, help="""Translation (in pixels) along y to define the --masks's center.""")
	#parser.add_argument("--tz", type=int, default=0, help="""Translation (in pixels) along z to define the --masks's center.""")
	
	parser.add_argument("--vertexes", action='store_true', default=False,help="""Only works if --sym=icos. This flag will make the program extract only the 12 vertexes from among all 60 symmetry related units.""") 
	
	#parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="align")
	
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")	

	(options, args) = parser.parse_args()
	
	if not options.mask and not options.maskfile:
		print("\nERROR: You must define --mask or supply a volume through --maskfile")
		sys.exit(1)

	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if options.plotcenters and not options.maskcenter:
		print("\nERROR: --plotcenters requires --maskcenter")
		sys.exit(1)

	if not options.input:
		parser.print_help()
		sys.exit(1)
	
	if not options.boxsize:
		print("\nERROR: --boxsize required. This ideally should be twice the size of the feature you're scooping out (and thus 4 times the radius you might be using to create a --mask)")
		sys.exit()
	#if options.savescoops and not options.boxsize:
	#	print("""ERROR: Specify the box size through --boxsize for the saved
	#	scoops.""")	
	#	sys.exit(1)

	if options.inputisraw and not options.transformfile:
		print("\nERROR: --inputisraw requires --transformfile")
		sys.exit(1)

	elif options.transformfile and not options.inputisraw:
		print("\nWARNING: --transformfile provided; therefore, --input is assumed to consist of unaligned raw particles; thus --inputisraw has been turned on")
		options.inputisraw = True
	
	#inputhdr = EMData(options.input,0,True)

	#If no failures up until now, initialize logger
	logid=E2init(sys.argv,options.ppid)
	
	options = makepath(options,"melonscoops")

	mask = getbasicmask(options)
	
	symletter,symnum = symparse(options)

	symorientations = getsymorientations(options,symletter,symnum)
	
	n=0
	nt=0
	ptclts = {}
	if options.transformfile:				
		ptclts=getptcltransforms(options)
		nt = len(ptclts)
		n = nt		
	
	centers = getcenters(options,symorientations)
	
	np = EMUtil.get_image_count( options.input )
	if nt>0:
		if np != nt:
			print("\nWARNING: fewer transforms in --transformfile, nt={} than particles in --input, np={}. Therefore, ignoring particles for which there is no transform data".format(nt,np))
	else:
		n = np
	
	outscoopsallfile = options.path + '/' + options.input.replace('.hdf','_scoops.hdf')
	
	totalcount=0
	ntargetscoops = n * len(symorientations) 
	for i in range(n):
		print("\nWorking with particle n={}/{}".format(i+1,n))
		ptcl = EMData(options.input,i)
		if options.inputisraw:
			ptclt = ptclts[i]
			print("\nusing this particle transform t={} (will invert it when appropriate)".format(ptclt))
		
		ptclnumtag = str(i).zfill(len( str(n) ))
		
		outscoopsperparticlefile = outscoopsallfile.replace('.hdf','_ptcl' + ptclnumtag + '.hdf')
		
		for k in range( len(symorientations) ):
			t = symorientations[k]
			c = centers[k]
		
			print("\nWorking with symorientation n={}/{}, t={}".format(k+1,len(symorientations),t))
			print("\ncenter={}".format(c))
		
			if options.inputisraw:
				ptclti = ptclt.inverse()
				if options.verbose:
					print("\nusing this INVERTED particle transform ti={}".format(ptclti))
				
				t = ptclti*t
				if options.verbose:
					print("\nWorking with this FINAL symorientation t={}".format(t))

			tmpmask = mask.copy()
			#print("\ncopied mask to tmpmask")
			tmpmask.transform( t )
			#print("\ntranformed tmpmask")
			
			tmpmask = origin2zero(tmpmask)
			#print("\nzeroed origin of tmpmask")
			
			sx = c[0]
			sy = c[1]
			sz = c[2]

			scoop = ptcl.copy()
			#print("\ncopied ptcl to scoop")
	
			scoop['spt_scoop_center'] = [sx,sy,sz]
			scoop['xform.align3d'] = t
			
			scoop.mult(tmpmask)
			#print("\ngot scoop by masking")
			
			if options.scoopsinsitu:
				scoopinsitu = scoop.copy()
				scoopinsitu['xform.align3d'] = Transform()

				if options.scoopsperparticle:
					outscoopsperparticlefileinsitu = outscoopsperparticlefile.replace('.hdf','_insitu.hdf')
					scoopinsitu.write_image(outscoopsperparticlefileinsitu,k)

				elif not options.scoopsperparticle:
					outscoopsallfileinsitu = outscoopsallfile.replace('.hdf','_insitu.hdf')
					scoopinsitu.write_image(outscoopsallfileinsitu,totalcount)

			#scoopclip = scoop.copy()
			scoop = clip3d(scoop,options.boxsize,[sx,sy,sz])
			#print("\nclipped scoop to boxsize={}".format(options.boxsize))
			scoop.write_image(outscoopsallfile,totalcount)
			print ("\nwrote out scoop {}/{}".format(totalcount,ntargetscoops))

			if options.scoopsorient:
				scooporiented = scoop.copy()
				ti = t.inverse()
				scooporiented.transform(ti)
				scooporiented['xform.align3d'] = ti
				scooporiented=origin2zero(scooporiented)

				if options.scoopsperparticle:
					outscoopsperparticlefileoriented = outscoopsperparticlefile.replace('.hdf','_oriented.hdf')
					scooporiented.write_image(outscoopsperparticlefileoriented,k)

				elif not options.scoopsperparticle:
					outscoopsallfileoriented = outscoopsallfile.replace('.hdf','_oriented.hdf')
					scooporiented.write_image(outscoopsallfileoriented,totalcount)

			totalcount+=1

	E2end(logid)
	
	return


def symparse(options):
	symnames = ['oct','OCT','icos','ICOS','tet','TET']
	symletter='c'
	symnum=1

	if options.sym:
		if options.verbose:
			print("\n--sym={}".format(options.sym))
			
		if options.sym not in symnames:
			
			symletter = options.sym
			symnum = options.sym
			
			for x in options.sym:
				if x.isalpha():
					symnum = symnum.replace(x,'')
					
			for x in options.sym:
				if x.isdigit():
					symletter = symletter.replace(x,'')
			
			
			print("\nThe letter for sym is", symletter)
			print("\nThe num for sym is", symnum)
		
		if options.sym == 'oct' or options.sym == 'OCT':
			symnum = 8
			
		if options.sym == 'tet' or options.sym == 'TET':
			symnum = 12	
			
		if options.sym == 'icos' or options.sym == 'ICOS':
			symnum = 60	
		
		symnum = int(symnum)
	
	return symletter,symnum


def getbasicmask(options):
	inputhdr = EMData(options.input,0,True)

	mask = EMData(inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])
	mask.to_one()
	
	if options.mask:	
		mask.process_inplace(options.mask[0],options.mask[1])
	
	elif options.maskfile:
		mask = EMData( options.maskfile )
		
		'''
		maskfile size
		'''
		maskx = mask['nx']
		masky = mask['ny']
		maskz = mask['nz']

		#hdr=EMData(options.input,0,True)
		
		hdrx=inputhdr['nx']
		hdry=inputhdr['ny']
		hdrz=inputhdr['nz']

		if int(hdrx) != int(maskx) or int(hdry) != int(masky) or int(hdrz) != int(maskz):
			'''
			Center of maskfile
			'''
			maskcx = old_div(mask['nx'],2.0)
			maskcy = old_div(mask['ny'],2.0)
			maskcz = old_div(mask['nz'],2.0)
			
			'''
			Clip to data's boxsize
			'''
			#r=Region( (2*maskcx - inputhdr['nx'])/2, (2*maskcx - inputhdr['ny'])/2, (2*maskcx - inputhdr['nz'])/2, inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])
			#mask.clip_inplace( r )
			mask = clip3d(mask,hdrx)
	
	if options.maskcenter:
		cstr = options.maskcenter.split(',')
		cx = int(round(float(cstr[0])))
		cy = int(round(float(cstr[1])))
		cz = int(round(float(cstr[2]))) 
		#mask.translate(0,0, float(options.radius) )

		mask.translate(cx,cy,cz)
		print("\nMask translated by x={}, y={}, z={}".format(cx,cy,cz))

	#elif options.tx or options.ty or options.tz:
	#	mask.translate( options.tx, options.ty, options.tz )
	#	print("\nMask translated by tx=%d, ty=%d, tz=%d,", options.tx, options.ty, options.tz)
	
	print("\nMask done")
	
	return mask


def getsymorientations(options,symletter,symnum):
	symorientations = []
	#transforms = []
	#anglelines = []

	t = Transform()
		
	if symnum:
		print("\nsymnum determined",symnum)
		print("while symletter is", symletter)
		
		if symletter == 'd' or symletter == 'D':
			symnum *= 2
			print("\nsymnum corrected, because symmetry is d",symnum)
			
		#if options.save
		
		for i in range(symnum):
			symorientation = t.get_sym( options.sym , i )
			symorientations.append( symorientation )
				
			#rot = symorientation.get_rotation()
			#az = rot['az']
			#alt = rot['alt']
			#phi = rot['phi']
			
			#anglesline = 'az=' + str(az) + 'alt=' + str(alt) + 'phi=' + str(phi) + '\n'
			#anglelines.append(anglesline)
			
			#transforms.append( Transform({'type':'eman','az':az,'alt':alt,'phi':phi}) )
	
	print("\ngenerated these many orientations n={}".format(len(symorientations)))
	if options.sym == 'icos' or options.sym == 'ICOS':
		if options.vertexes:

			print("\nbut fetching vertexes only")
		
			#orientations = genicosvertices ( orientations )

			symorientations = genicosvertexesnew ( symorientations )

	return symorientations


def getptcltransforms(options):
	jsonfile = options.transformfile
	jsonfileopen = js_open_dict( jsonfile )
	#print "\njsonfile to open = {}".format(jsonfile)

	n = int(EMUtil.get_image_count(options.input))

	nt = int(len(jsonfileopen))
	
	if int(n) != int(nt):
		print("\n WARNING: the number of particles in --input n={} does not match the number of transforms in --transformfile nt={}".format(int(n),int(nt)))
		if int(n) < int(nt):
			print("\nWARNING: fewer particles n={} than transform parameters nt={}. Only masking as many particles from --input={} as transforms in --transformfile={}".format(n,nt,options.input,options.transformfile))
		elif int(n) > int(nt):
			n = int(nt)
			print("\nWARNING: more particles n={} than transform parameters nt={}. Ignoring particles for which there is no transform data in --transformfile={}".format(n,nt,options.transformfile))

	ptclts = {}

	keys = list(jsonfileopen.keys())
	keys.sort()
	#n = len(keys)
	for j in range(0,n):
		label = keys[j]
		ptclt = jsonfileopen[label][0]
		ptclts.update({j:ptclt})
	
	return ptclts


def genicosvertexesnew( syms ):

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


def getcenters(options,symorientations):
	
	centerx = 0
	centery = 0
	centerz = 0

	if options.maskcenter:
		#centerz = float( options.radius )
		cstr = options.maskcenter.split(',')
		centerx = int(round(float(cstr[0])))
		centery = int(round(float(cstr[1])))
		centerz = int(round(float(cstr[2]))) 
	
	groundv = Vec3f(centerx,centery,centerz)
	print("\ngroundv={}".format(groundv))

	hdr = EMData(options.input,0,True)
	origboxsize=hdr['nx']
	print("\norigboxsize={}".format(origboxsize))

	boxcenterv = Vec3f(old_div(origboxsize,2.0),old_div(origboxsize,2.0),old_div(origboxsize,2.0))
	print("\nboxcenterv={}".format(boxcenterv))
	
	centersmap = EMData(origboxsize,origboxsize,origboxsize)
	centersmap.to_zero()
	
	centers = {}
	centerlines = []
	for k in range( len(symorientations) ):
		t =  symorientations[k]
		
		print("\nWorking with this symorientation",t)
	
		centerv = t*groundv
		
		#we translate the asymmetric unit 'center' (rotated by t) to be anchored at the center of the original box
		center = centerv + boxcenterv
		print("\ncenter={}".format(center))

		centers.update( {k: center} )
		
		#tout = Transform()
		#rot = t.get_rotation()
		#tout.set_rotation( rot )

		#tout.set_trans( newcenterx, newcentery, newcenterz )
		
		if options.plotcenters:
			centerline = str(center[0]) + ' ' + str(center[1]) + ' ' + str(center[2]) + '\n'
			centerlines.append(centerline)

			centersmap.set_value_at( int( round(center[0] )) , int( round(center[1] )), int( round(center[2] )), 1.0 )
	
	if options.plotcenters:
		with open( options.path + '/centers.txt', 'a') as f: f.writelines(centerlines)
	
		centersmap.write_image(options.path + '/centersmap.hdf',0)

	return centers


def maskptclsglobally(options,symorientations):
	finalmask = EMData(mask['nx'],mask['ny'],mask['nz'])
	finalmask.to_zero()

	ptclglobalmsk = ptcl.copy()
	ptclglobalmsk.mult( finalmask )
		
	ptclglobalmsk.write_image( options.path + '/' + 'ptcls_globallymasked.hdf',i)
	finalmask.write_image(options.path + '/finalmask.hdf',0)

	tmpmask = mask.copy()
	tmpmask.transform( t )
	
	tmpmask['origin_x'] = 0
	tmpmask['origin_y'] = 0
	tmpmask['origin_z'] = 0
	
	finalmask = finalmask + tmpmask

	tmpmask['spt_scoop_center'] = [ newcenterx, newcentery, newcenterz ]
	tmpmask['xform.align3d'] = tout
	
	tmpmask.write_image(options.path + '/masks.hdf',k)
	
	masks.update( {k:[tmpmask,tout,[newcenterx,newcentery,newcenterz]]} )
	
	return

'''
SCRAPS


scoop.mult( thismask )
				
				box = options.boxsize
	
				extra = math.fabs( math.sqrt(2.0) * ( box/2.0 - 1.0 ) - ( box/2.0 - 1 ) )	
				paddedbox = int( math.ceil( float(box) + extra ) )
				
				#print "Center for region is", sx,sy,sz
					
				print("\nExtracting this scoop", key)
				print("With a padded box of this size", paddedbox)
				bigr = Region( (sx*2 - paddedbox)/2 ,  (sy*2 - paddedbox)/2 ,  (sz*2 - paddedbox)/2, paddedbox,paddedbox,paddedbox)
				bigscoop = scoop.get_clip(bigr)
				
				#bigscoop.write_image( options.path + '/bigscoops.hdf', key)
				
				bigscoop['origin_x'] = 0
				bigscoop['origin_y'] = 0
				bigscoop['origin_z'] = 0
				
				print("\nOrienting the subscoop with this transform", t)
				
				t2use = Transform()
				rot = t.get_rotation()
				
				t2use.set_rotation( rot )
				
				ti = t2use.inverse()
				bigscoop.transform( ti )
				
				#bigscoop.write_image( options.path + '/bigscoopsOriented.hdf', key)
				
				sxB = bigscoop['nx']/2
				syB = bigscoop['ny']/2
				szB = bigscoop['nz']/2
				
				print("\nCenter of bigs scoops is", sxB,syB,szB)
				
				scoop = clip3d( bigscoop, box )
				
				#r = Region( (sxB*2 - box)/2 ,  (syB*2 - box)/2 ,  (szB*2 - box)/2, box,box,box)
				#scoop = bigscoop.get_clip(r)
				
				#print "\nTherefore region for small scoop is", r
				
				print("\nClipping the extracted and oriented scoop back to the desired boxsize", box)
				
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
	
