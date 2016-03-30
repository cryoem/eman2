#!/usr/bin/env python

#
# Author: Jesus Galaz, 12/08/2011 - Last Update 10/Feb/2014
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

from EMAN2 import *
import os
from e2spt_classaverage import *
import sys


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <stack of particles> [options] . This programs autocenters particles
		using autocorrelation or spherical averages."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input",type=str,default='',help="""Input stack in HDF format.""")
	
	parser.add_argument("--path",type=str,default='spt_autoc',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'sptautobox'; 
		for example, sptautobox_01 will be the directory by default if 'sptautobox' already exists.""")
	
	parser.add_argument("--boxclip", type=int, default=0, help="""The boxsize to clip the FINAL boxes 
		down to after auto-centering, so that they contain no empty pixels. Default=0 or 'off'""")
	
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	#parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")

	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--mask",type=str,default='',help="Mask to apply to the dataset average before spherically averaging it.")
	
	parser.add_argument("--shrink",type=int,default=0,help="Optionally shrink the volume(s) to have the FFTs go faster")
	
	parser.add_argument("--lowpass",type=str,help="A filter applied to the average of the dataset before spherically averaging it", default='')
	
	parser.add_argument("--highpass",type=str,help="A filter applied to the average of the dataset before spherically averaging it", default='')

	parser.add_argument("--preprocess",type=str,help="A filter applied to the average of the dataset before spherically averaging it", default='')

	parser.add_argument("--threshold",type=str,help="A filter applied to the average of the dataset before spherically averaging it", default='')
	
	parser.add_argument("--processptcls",type=int,default=0,help="Apply all pre-processing (mask, shrink, filter, etc...) on the particles before aligning them to the previous spherical average.")
	
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean.tomo")
	
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False)
	
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False)
	#parser.add_argument("--align", help="Set by default",default="translational:intonly=1")
	
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	
	parser.add_argument("--mode",type=str,default='autocorr', help="""Provide --mode=sphere or --mode=autocorr.
		The first one iteratively aligns the sperical averages of individual particles to 
		the speherical average of the entire set-average. 
		The second generates mirror images of each particle in 3-D and 2-D and uses
		autocorrelation to center the particles.
		You can actually provide both modes and they will be applied sequentially in the
		order you provided them. For example, --mode=autocorr,spherical""")

	parser.add_argument("--align", help="Set by default",default="rotate_translate_3d_grid:alt0=0:alt1=1:az0=0:az1=1:dalt=2:daz=2:dotrans=1:dphi=2:phi0=0:phi1=1:search=10:verbose=1")

	(options, args) = parser.parse_args()


	hdr = EMData(options.input,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
	oldbox=nx


	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options )
	
	if options.averager: 
		options.averager=parsemodopt(options.averager)
	
	if options.mask: 
		options.mask=parsemodopt(options.mask)

	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)

	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)
		
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	
	if options.align: 
		options.align=parsemodopt(options.align)
	
	if options.aligncmp: 
		options.aligncmp=parsemodopt(options.aligncmp)
	
	n=EMUtil.get_image_count(options.input)
	
	modes=options.mode.split(',')
	
	for m in modes:
		if 'autocorr' in m:
			print "Mode is",m
			print "Number of ptcls to autocenter is", n
			print "In file", options.input
			for i in range(n):
				print "Autocentering ptcl number", i
				e = EMData( options.input, i )
				
				ret=autocorrcenter(e,options)						
				outname = options.input.replace('.hdf','_centered.hdf')
	
				ret[0].write_image(options.path + '/' + outname,i)
				ret[1].write_image(options.path + '/' + 'prjsTopRaw.hdf', i)
				ret[2].write_image(options.path + '/' + 'prjsSideRaw.hdf',i)
				ret[3].write_image(options.path + '/' + 'prjsTopCent.hdf',i)
				ret[4].write_image(options.path + '/' + 'prjsSideCent.hdf',i)
				
		
		elif 'sphere' in m:
			print "Mode is",m
			spherical(options)
	
	return
	

def autocorrcenter(e,options):

	tside = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
		
	'''
	Generate top projections to center X and Y
	'''

	etmp = e.copy()

	ecent = e.copy()
	
	eb = e['nx']
	
	#mask=EMData(eb,eb,eb)
	#mask.to_one()
	#mask.process_inplace(options.mask[0],options.mask[1])
	
	#etmp.process_inplace('normalize.mask',{'mask':mask})
	#etmp.mult(mask)
	
	t=Transform()
	pprj=etmp.project("standard",t)
	pprjsideOr=etmp.project("standard",tside)
	
	pprj['apix_x'] = e['apix_x']
	pprj['apix_y'] = e['apix_y']

	pmx=pprj.process('xform.mirror',{'axis':'x'})
	pmy=pprj.process('xform.mirror',{'axis':'y'})

	ccfpmx=pprj.calc_ccf(pmx)
	ccfpmy=pprj.calc_ccf(pmy)

	ccfpmxC=ccfpmx.process('xform.phaseorigin.tocorner')
	ccfpmyC=ccfpmy.process('xform.phaseorigin.tocorner')

	maxccfpmxC=ccfpmxC.calc_max_location()
	maxccfpmyC=ccfpmyC.calc_max_location()

	xt=(eb/2.0 - maxccfpmxC[0])/2.0
	yt=(eb/2.0 - maxccfpmyC[1])/2.0

	pfix=pprj.copy()
	
	pfix.translate( xt, yt, 0)
	
	ecent.translate(xt,yt,0)

	#newx=x+xt
	#newy=y+yt

	'''
	Generate side projections to center Z
	'''
	
	pprjside=ecent.project("standard",tside)
	
	
	pprjside['apix_x'] = e['apix_x']
	pprjside['apix_y'] = e['apix_y']

	pmxside=pprjside.process('xform.mirror',{'axis':'x'})
	pmzside=pprjside.process('xform.mirror',{'axis':'y'})

	ccfpmxside=pprjside.calc_ccf(pmxside)
	ccfpmzside=pprjside.calc_ccf(pmzside)

	ccfpmxCside=ccfpmxside.process('xform.phaseorigin.tocorner')
	ccfpmzCside=ccfpmzside.process('xform.phaseorigin.tocorner')

	maxccfpmxCside=ccfpmxCside.calc_max_location()
	maxccfpmzCside=ccfpmzCside.calc_max_location()

	xtside = (eb/2.0 - maxccfpmxCside[0])/2.0
	ztside = -1 * (eb/2.0 - maxccfpmzCside[1])/2.0
	#print "zt side is", ztside
	
	#pfixside=pprjside.copy()
	#pfixside.translate( xtside, ztside, 0)
	
	#pfix.translate( xtside,0,0)
	
	ecent.translate(xtside, 0, ztside)
	
	pfixside = ecent.project("standard",tside)
	pfix = ecent.project("standard",t)
	
	print "Translations are", xt+xtside, yt, ztside
	
	
	#mskr = -1*max( math.fabs(xt), math.fabs(yt), math.fabs(ztside) )
	
	#ecent.process_inplace('mask.sharp',{'outer_radius':mskr})
	#pfix.process_inplace('mask.sharp',{'outer_radius':mskr})
	#pfixside.process_inplace('mask.sharp',{'outer_radius':mskr})
	
	
	'''
	#newx=x+xt
	newz=z-ztside

	if math.fabs(xt) <= eb/8.0 and math.fabs(yt) <= eb/8.0 :
		#newestcoords=(newx,newy,newz)
	
		newestccf = ccfpmy['maximum']
		if ccfpmx['maximum'] > ccfpmy['maximum']:
			newestccf = ccfpmx['maximum']
	
		#print "\n$$$$$$$$$$$$$$$$\nAUTOCORRELATION PARTICLE WRITTEN!\n"

		#prjccfs.append(newccf)

		#coordset.add(coords)
		#coeffset.add(newcoefftuple)
		#newestdata.append( [newestccf,newestcoords] )
	
		r = Region((2*newx- options.outputboxsize)/2,(2*newy-options.outputboxsize)/2, (2*newz-options.outputboxsize)/2, options.outputboxsize, options.outputboxsize, options.outputboxsize)
		e = EMData()
		e.read_image(options.tomogram,0,False,r)

		eb = e['nx']
		mask=EMData(eb,eb,eb)
		mask.to_one()
		mask.process_inplace('mask.sharp',{'outer_radius':-2})
		e.process_inplace('normalize.mask',{'mask':mask})
		e.mult(mask)

		pprj=e.project("standard",t)
		pprj['apix_x'] = e['apix_x']
		pprj['apix_y'] = e['apix_x']
		pprj.write_image(options.path + '/' + 'pprj_corrected_sorted_autocorrelated.hdf',ppp)	
					
		mean=pprj['mean_nonzero']
		#means.append(mean)
		#print "The mean_nonzero is %f, compared to mean %f" %( mean, pprj['mean'] )

		sigma=pprj['sigma_nonzero']
		#sigmas.append(sigma)

		max=pprj['maximum']
		maxs2d.append(max)

		min=pprj['minimum']
		#mins2d.append(min)
		if (newestccf,newx,newy,newz) not in basicdata:			
			newestdata.append([newestccf,newx,newy,newz,ppp,max,min,sigma,mean,pprj])
			basicdata.add( (newestccf,newx,newy,newz) )
			ppp+=1
		else:
			print "TRYing to add repeated particle in second PRJ loop!"	
	
	else:
		rrr+=1
		#print "\nParticle eliminated because translation from AUTOCORRELATION were too big!\n"
		pass
	mmm+=1
	#if d[0] < prjccfs_mean - prjccfs_sigma:
	#	newdata.remove(d)
	#	print "I have removed a particle based new PRJ mask"

	print "The number of particles pruned by AUTOCORRELATION is", rrr		
	print "Therefore, AFTER AUTUCORRELATION prunning, the len of data is", len(newestdata)
	'''
	
	return(ecent,pprj,pprjsideOr,pfix,pfixside,xt,yt,ztside)

	
def spherical(options):
	print "Inside spherical"
	n = EMUtil.get_image_count(options.input)

	'''
	Generate the brute raw average of all the particles
	'''
	avgrawf = options.path + '/avgsRaw.hdf'
	cmd = 'e2proc3d.py ' + options.input + ' ' + avgrawf + ' --average'
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	'''
	Load the brute raw average just computed and perform any necessary preprocessing
	'''
	avgraw = EMData( avgrawf, 0 )
	
	avgprep = avgraw.copy()
	
	outname = options.path+'/'+options.input.replace('.hdf','_SphCentered.hdf')
	
	for it in range(options.iter):
		print "Interation number", it
				
		if options.mask or options.shrink or options.lowpass or options.highpass or options.mask:
			avgprep = acpreprocessing(options,avgraw)
	
		'''
		Compute the rotational spherical average for the brute raw average and generate 
		a top-view and a side-view projection of it (they should be the same actually).
		'''
		avgsph=avgprep.rotavg_i()
		
		avgsph.process_inplace('xform.centerofmass')
		
		if options.savesteps:
			avgsph.write_image( options.path + '/avgsSph.hdf',-1)
		
		t=Transform()
		tside=Transform({'type':'eman','alt':-90})
		
		asprj = avgsph.project('standard',t)
		asprjside = avgsph.project('standard',tside)
		
		asprj['apix_x'] = avgraw['apix_x']
		asprjside['apix_y'] = avgraw['apix_y']
		
		
		'''
		Iterate aligning particle projections to the spherical average projections
		'''
		
		avgr=Averagers.get(options.averager[0],options.averager[1])
		
		for k in range(n):
			p=EMData(options.input,k)
			print "Read particle number", k
			
			ptrans=p.copy() #This will be the new centered particle
			
			pprep = p.copy() #This will be a preprocessed copy of the particle if needed
			
			pb = p['nx'] #box size of the particles
						
			'''
			Preprocess the particles if necessary
			'''
			if options.mask or options.shrink or options.lowpass or options.highpass or options.mask:
				print "Will preprocess"
				pprep = acpreprocessing(options,p)
			
			'''
			Generate top-view projection for the particle in turn and align to the spherical 
			raw average's top-view projection
			'''
			pprj = pprep.project('standard',t)
			pprj['apix_x'] = p['apix_x']
		
			ccfp = asprj.calc_ccf( pprj )
			ccfpC=ccfp.process('xform.phaseorigin.tocorner')
			
			maxccfpC = ccfpC.calc_max_location()
			
			xt=(pb/2.0 - maxccfpC[0])
			yt=(pb/2.0 - maxccfpC[1])
			
			
			'''
			Apply translations found to p and pprep
			'''
			print "Top translation are",xt,yt
			
			ptrans.translate(xt,yt,0)
			pprep.translate(xt,yt,0)
			
			
			'''
			Generate side-view projection for the particle in turn and align to the spherical 
			raw average's side-view projection
			'''
			pprjside = pprep.project('standard',tside)
			pprjside['apix_y'] = p['apix_y']
		
			ccfpside = asprjside.calc_ccf( pprjside )
			ccfpsideC = ccfpside.process('xform.phaseorigin.tocorner')
			maxccfpsideC = ccfpsideC.calc_max_location()
			
			xtside = (pb/2.0 - maxccfpsideC[0])/2.0
			ztside = -1 * (pb/2.0 - maxccfpsideC[1])/2.0 
			#Y-shifts in a side view correspond to -z shifts in a 3-D model rotated by -90 in altitude
		
			
			'''
			Apply translations found to p and pprep
			'''
			ptrans.translate(xtside,0,ztside)
			pprep.translate(xtside,0,ztside)	
			
			print "Side trans applied are",xtside,ztside
			
			avgr.add_image( ptrans )
			
			if it == int(options.iter) -1:
				print "Will write translated raw image to", outname
				ptrans.write_image(outname,-1)
				
		avgprep=avgr.finish()
		print "Finalized average; will normalize it"
		avgprep.process_inplace('normalize')
		
		if options.savesteps:
			print "Will write raw average"
			avgprep.write_image(options.path + '/avgsRaw.hdf',-1)
	
	
	return
	
	
def acpreprocessing(options,image):
	
	#print "\n$$$$$$$\nIn preprocessing, received options and image, types", type(options), type(image)
	
	'''
	Make the mask first 
	'''
	mask=EMData( int(image["nx"]), int(image["ny"]), int(image["nz"]) )
	mask.to_one()
	
	if options.mask:
		#if options.verbose:
			#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
		mask.process_inplace(options.mask[0],options.mask[1])
		
	
	'''
	Set the 'mask' parameter for --normproc if normalize.mask is being used
	'''
	if options.normproc:
		if options.normproc[0]=="normalize.mask": 
			options.normproc[1]["mask"]=mask

	simage = image.copy()
	
	if options.mask:
		#if options.shrink:
		#	maskCoarse = mask.copy()
		#	maskCoarse.process_inplace('math.meanshrink',{'n':options.shrink})
		simage.mult(mask)
	
	if options.normproc:
		simage.process_inplace(options.normproc[0],options.normproc[1])
		
	if options.mask:
		simage.mult(mask)

	if options.lowpass or options.highpass or options.preprocess or options.shrink:
		simage = filters(simage,options.preprocess,options.lowpass,options.highpass,options.shrink)
	
	if options.threshold:
		simage.process_inplace(options.threshold[0],options.threshold[1])	
	
	return simage


if __name__ == "__main__":
	main()
