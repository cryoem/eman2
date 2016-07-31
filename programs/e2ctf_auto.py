#!/usr/bin/env python

#
# Author: Steven Ludtke, 05/28/16 (sludtke@bcm.edu)
# Copyright (c) 2016- Baylor College of Medicine
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

from EMAN2 import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from optparse import OptionParser
from math import *
import os
import sys
import weakref
import traceback
from numpy import array,arange
import numpy
import threading

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])

	usage = """prog [options]
This program automates the CTF fitting and structure factor generation process, which normally involves a 
sequence of at least 4 different steps. For most projects, this will work correctly with no human intervention.
We strongly recommend running the GUI after it completes to double-check the fitting of a few of the closest
and a few of the furthest from focus images by hand.

If you detect fitting problems on any images, manually adjust those images to a defocus near the correct value,
then rerun this program (without --fromscratch) and the problem should be fixed.

For defocus with astigmatism fitting, whole frame fitting as implemented in e2rawdata.py and
e2evalimage.py generally produce superior fits. If you have run this whole-frame fitting, e2ctf_auto will
take this into account, and only fit the other required parameters. 

Important: This program must be run from the project directory, not from within the particles directory
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="modeheader", help='Processing can target different final resolution ranges', title="Select one of hires,midres or lores:", row=0, col=0, rowspan=1, colspan=3, mode="auto")
	parser.add_argument("--hires",action="store_true",help="Perform CTF processing for projects targeting 2-6 A resolution",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--midres",action="store_true",help="Perform CTF processing for projects targeting 7-15 A resolution",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--lores",action="store_true",help="Perform CTF processing for projects targeting 15-30 A resolution",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='auto[False]')

	parser.add_header(name="modeheader", help='Image Parameters', title="Image Parameters", row=10, col=0, rowspan=1, colspan=3, mode="auto")
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=0, guitype='floatbox', row=12, col=0, rowspan=1, colspan=1, mode="auto['self.pm().getAPIX()']")
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=0, guitype='floatbox', row=14, col=0, rowspan=1, colspan=1, mode="auto['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=0, guitype='floatbox', row=14, col=1, rowspan=1, colspan=1, mode="auto['self.pm().getCS()']")

	parser.add_header(name="modeheader", help='Additional Options', title="Additional Options", row=20, col=0, rowspan=1, colspan=3, mode="auto")
	parser.add_argument("--fromscratch",action="store_true",help="Force refitting of CTF from scratch, ignoring any previous fits.",default=False, guitype='boolbox', row=22, col=0, rowspan=1, colspan=1, mode='auto[False]')
#	parser.add_argument("--forceframe",action="store_true",help="Uses defocus/astigmatism from frames. Does not significantly modify based on particles",default=False, guitype='boolbox', row=24, col=2, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--astigmatism",action="store_true",help="Includes astigmatism in automatic fitting (use e2rawdata first)",default=False, guitype='boolbox', row=22, col=1, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--extrapad",action="store_true",help="If particles were boxed more tightly than EMAN requires, this will add some extra padding",default=False, guitype='boolbox', row=22, col=2, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--highdensity",action="store_true",help="If particles are very close together, this will interfere with SSNR estimation. If set uses an alternative strategy, but may over-estimate SSNR.",default=False, guitype='boolbox', row=24, col=0, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--invert",action="store_true",help="Invert the contrast of the particles in output files (default false)",default=False, guitype='boolbox', row=24, col=1, rowspan=1, colspan=1, mode='auto[False]')
	parser.add_argument("--defocusmin",type=float,help="Minimum defocus in autofitting",default=0.6, guitype='floatbox', row=26, col=0, rowspan=1, colspan=1, mode="autofit[0.6]")
	parser.add_argument("--defocusmax",type=float,help="Maximum defocus in autofitting",default=4, guitype='floatbox', row=26, col=1, rowspan=1, colspan=1, mode='autofit[4.0]')
	parser.add_argument("--constbfactor",type=float,help="Set B-factor to fixed specified value, negative value autofits",default=-1.0, guitype='floatbox', row=28, col=0, rowspan=1, colspan=1, mode='auto[-1.0]')
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default 10)",default=10, guitype='floatbox', row=28, col=1, rowspan=1, colspan=1, mode="auto[10]")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on the local computer",guitype='intbox', row=30, col=0, rowspan=1, colspan=1, mode='auto[4]')
	parser.add_argument("--minqual",type=int,help="Files with a quality value lower than specified will be skipped",default=0,guitype='intbox', row=30, col=1, mode='auto[0]')

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	osum=int(options.hires)+int(options.midres)+int(options.lores)
	if osum!=1 :
		print "ERROR: Please specify one of --hires --lores --midres"
		sys.exit(1)

	ptcls=["particles/"+i for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
	if options.verbose : print "%d particle stacks identified"%len(args)

	if options.highdensity : highdensity="--highdensity"
	else: highdensity=""

	if options.constbfactor<0 :
		constbfactor=""
	else:
		constbfactor="--constbfactor {:02f}".format(options.constbfactor)

	# Check to see what we're dealing with. If we have frame-based parameters, we take that into account. 
	# If the frame parameters have astigmatism, then we don't adjust the defocus or astigmatism. We only check the first 10 frames
	# After this, frame_ctf true if this info available, and frame_stig true if any of the first 10 frames had non-zero astigmatism
	db=js_open_dict(info_name(ptcls[0]))
	frame_stig=False
	frame_ctf=False
	if db.has_key("ctf_frame") : 
		frame_ctf=True
		for f in ptcls[:10]:
			db.close()
			db=js_open_dict(info_name(f))
			try: fc=db["ctf_frame"][1]
			except:
				frame_ctf=False
				break
			if fc.dfdiff!=0: frame_stig=True
			
		if options.voltage==0 : options.voltage=fc.voltage
		if options.cs==0 : options.cs=max(fc.cs,0.01)			# CTF model doesn't work well with Cs exactly 0
		if options.apix==0 : options.apix=fc.apix
		
		if options.voltage!=fc.voltage or fabs(options.cs-fc.cs)>0.02 or fabs(options.apix-fc.apix)>0.1 or options.ac!=fc.ampcont :
			print """Warning: Disagreement in specified voltage, Cs, A/pix or %AC between frames and options. This requires refitting without frame-based parameters.
Strongly suggest refitting CTF from frames with e2rawdata.py with revised parameters before running this program"""
			frame_ctf=False
	
	# fill in missing parameters from project if possible
	project=js_open_dict("info/project.json")
	if options.voltage==0 : 
		try:
			options.voltage=float(project["global.microscope_voltage"])
			print "Using project voltage = ",options.voltage
		except:
			print "Error, no voltage found"
			sys.exit(1)
	
	if options.cs==0 : 
		try:
			options.cs=float(project["global.microscope_cs"])
			print "Using project Cs = ",options.cs
		except:
			print "Error, no Cs found"
			sys.exit(1)
	
	if options.apix==0 : 
		try:
			options.apix=float(project["global.apix"])
			print "Using project A/pix = ",options.apix
		except:
			print "Error, no A/pix found"
			sys.exit(1)
	
	try:
		tmp=EMData(ptcls[0],0)
		boxsize=tmp["nx"]
	except:
		print "ERROR: Couldn't read first particle from ",ptcls[0]
		sys.exit(3)
	
	###
	# This is the logic governing how we handle automatic fitting based on existing information
	if options.fromscratch :
		print "--fromscratch specified, so fitting will ignore existing values"
		if options.astigmatism :
			print """Astigmatism fitting from particles requested. This may be fine for strong astigmatism through midres 
resolution, but for high resolution work, fitting defocus/astig from frames is recommended"""
			fit_options="--astigmatism"
		else:
			fit_options=""
	else:
		if frame_ctf : 
			print "Frame based CTF parameters found"
			if frame_stig : 
				if options.astigmatism : 
					print "Astigmatism present in in frame parameters. Will use defocus/astig unchanged from frames"
					fit_options="--curdefocusfix --astigmatism --useframedf"
				else :
					print """Astigmatism present in frame parameters, but not specified here. 
	No astigmatism will be used, and defocuses will be refined from particles.
	Astigmatism info will be preserved in frame parameters, so this program may be re-run later to take astigmatism into account"""
					fit_options="--curdefocushint"
			else :
				if options.astigmatism :
					print "Frame parameters present, but without astigmatism. Will refine frame-based defocus from particle data."
					print """Astigmatism fitting from particles requested.""" 
					if options.hires : print """This may be fine for strong astigmatism through midres 
resolution, but for high resolution work, fitting defocus/astig from frames is recommended"""
					fit_options="--curdefocushint --astigmatism --useframedf"
				else :
					print "Frame parameters present without astigmatism. Slightly refining frame defocus values from particles."
					fit_options="--curdefocusfix --useframedf"	# "slightly refining" refers to refinebysnr
		else: 
			print "No frame-based CTF parameters found. Fitting from particles."
			if options.astigmatism :
				print """Astigmatism fitting from particles requested.""" 
				if options.hires : print """This may be fine for strong astigmatism through midres 
resolution, but for high resolution work, fitting defocus/astig from frames is recommended"""
				fit_options="--curdefocushint --astigmatism"
			else:
				fit_options="--curdefocushint"
			
	logid=E2init(sys.argv, options.ppid)
	E2progress(logid,0)

	###
	# run CTF autofit
	com="e2ctf.py --autofit --allparticles --oversamp 2 --apix {apix} --voltage {voltage} --cs {cs} --ac {ac} --defocusmin {dfmin} --defocusmax {dfmax} --threads {threads} --minqual {minqual} {highdensity} {constbfactor} {fit_options}".format(
		apix=options.apix,voltage=options.voltage,cs=options.cs,ac=options.ac,threads=options.threads,highdensity=highdensity,constbfactor=constbfactor,fit_options=fit_options,minqual=options.minqual,dfmin=options.defocusmin,dfmax=options.defocusmax)
	if options.verbose: print com
	launch_childprocess(com)
	E2progress(logid,0.25)
	
	###
	# decide which particle files we will use for structure factor
	dfvals=[]
	for f in ptcls:
		db=js_open_dict(info_name(f))
		try:
			ctf=db["ctf"][0]
			qual=db["quality"]
			if qual<options.minqual : 
				db.close()
				continue
			dfvals.append((ctf.defocus,f))
		except:
			print "CTF seems to have failed for ",f
		db.close()

	if len(dfvals)<5 :
		print "ERROR: This program requires good CTF fits from at least 5 micrographs to continue"
		sys.exit(2)
	
	# Strip off the low and high 10% of defocuses, since this is where failures usually occur
	dfsel=dfvals[len(dfvals)//10:-len(dfvals)//10]

	# roughly 25 files should be more than enough for a structure factor
	step=max(1,len(dfsel)//25)
	if step>1 : dfsel=dfsel[::step]

	###
	# Structure factor computation
	com="e2ctf.py --computesf {}".format(" ".join([i[1] for i in dfsel]))
	if options.verbose: print com
	launch_childprocess(com)
	E2progress(logid,0.5)

	###
	# run CTF autofit, now with structure factor available
	com="e2ctf.py --autofit --allparticles --oversamp 2 --apix {apix} --voltage {voltage} --cs {cs} --ac {ac} --defocusmin {dfmin} --defocusmax {dfmax} --threads {threads} --minqual {minqual} {highdensity} {constbfactor} {fit_options}".format(
		apix=options.apix,voltage=options.voltage,cs=options.cs,ac=options.ac,threads=options.threads,highdensity=highdensity,constbfactor=constbfactor,fit_options=fit_options,minqual=options.minqual,dfmin=options.defocusmin,dfmax=options.defocusmax)
	if options.verbose: print com
	launch_childprocess(com)
	E2progress(logid,0.65)

	###
	# refinebysnr if appropriate
	if (options.midres or options.hires) and not options.astigmatism :
		com="e2ctf.py --allparticles --refinebysnr"
		if options.verbose: print com
		launch_childprocess(com)
	E2progress(logid,0.7)
	
	###
	# Generate output

	# For "small" data, we target ~6 A/pix
	resample1=6.0/options.apix
	newbox=good_size(boxsize/resample1)
	resample1=boxsize/(newbox+0.1)	# 0.1 is to prevent roundoff issues
	if resample1<1.0:
		print "Warning: original particle A/pix is very large!"
		resample1=1.0
	maskwid1=20.0/options.apix
	maskrad1=int(boxsize/2-maskwid1)

	# for "medium" data, we target ~2.6 A/pix
	resample2=2.6/options.apix
	newbox=good_size(boxsize/resample2)
	resample2=boxsize/(newbox+0.1)
	if resample2<1.0:
		if not options.lores : print "Warning: original sampling is too large for ideal midres resolution results. Suggest <=2.6 A/pix"
		resample2=1.0
	maskwid2=12.0/options.apix
	maskrad2=int(boxsize/2-maskwid2*1.2)

	# for "low" second data, we target ~4 A/pix
	resample3=4.0/options.apix
	newbox=good_size(boxsize/resample3)
	resample3=boxsize/(newbox+0.1)
	if resample3<1.0:
		resample2=1.0
	maskwid3=18.0/options.apix
	maskrad3=int(boxsize/2-maskwid2*1.2)

	# for high resolution data, we still go ahead and do some masking
	maskwid3=6.0/options.apix
	maskrad3=int(boxsize/2-maskwid3*1.2)

	if options.invert: invert="--invert"
	else: invert=""

	if options.extrapad : extrapad="--extrapad"
	else : extrapad=""

	if options.lores :
		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp20 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.05 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad1,maskwid=maskwid1,resamp=resample1,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		E2progress(logid,0.85)

		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp12 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.08333 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad3,maskwid=maskwid3,resamp=resample3,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		
		print "Phase-flipped output files:\n__ctf_flip_lp20 - masked, downsampled, filtered to 20 A resolution\n__ctf_flip_lp7 - masked, downsampled, filtered to 7 A resolution"
		
	elif options.midres:
		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp20 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.05 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad1,maskwid=maskwid1,resamp=resample1,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		E2progress(logid,0.8)

		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp7 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.14 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad2,maskwid=maskwid2,resamp=resample2,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		E2progress(logid,0.9)

		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag fullres --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 normalize.circlemean:radius={maskrad} --phaseflipproc3 mask.soft:outer_radius={maskrad}:width={maskwid} {extrapad}".format(
			maskrad=maskrad3,maskwid=maskwid3,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		print "Phase-flipped output files:\n__ctf_flip_lp20 - masked, downsampled, filtered to 20 A resolution\n__ctf_flip_lp7 - masked, downsampled, filtered to 7 A resolution\n__ctf_flip_fullres - masked, full sampling"
		
	else :
		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp14 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.7 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad1,maskwid=maskwid1,resamp=resample1,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		E2progress(logid,0.8)

		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag lp5 --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.2 --phaseflipproc3 normalize.circlemean:radius={maskrad} --phaseflipproc4 mask.soft:outer_radius={maskrad}:width={maskwid} --phaseflipproc5 math.fft.resample:n={resamp} {extrapad}".format(
			maskrad=maskrad2,maskwid=maskwid2,resamp=resample2,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		E2progress(logid,0.9)

		com="e2ctf.py --allparticles {invert} --minqual={minqual} --proctag fullres --phaseflipproc filter.highpass.gauss:cutoff_pixels=3 --phaseflipproc2 normalize.circlemean:radius={maskrad} --phaseflipproc3 mask.soft:outer_radius={maskrad}:width={maskwid} {extrapad}".format(
			maskrad=maskrad3,maskwid=maskwid3,invert=invert,minqual=options.minqual,extrapad=extrapad)
		if options.verbose: print com
		launch_childprocess(com)
		print "Phase-flipped output files:\n__ctf_flip_lp14 - masked, downsampled, filtered to 14 A resolution\n__ctf_flip_lp5 - masked, downsampled, filtered to 5 A resolution\n__ctf_flip_fullres - masked, full sampling"

	print "Building default set with all particles for convenience"
	com="e2buildsets.py --setname=all --excludebad --allparticles"
	if options.verbose: print com
	launch_childprocess(com)
	
	E2end(logid)
	print """
The length of your particle on its longest axis should be <= {}. If your particle is larger than this, 
you either need a larger box-size, or you will need to proceed with CTF fitting manually (and will likely achieve
suboptimal results). See http://www.eman2.org/emanwiki/EMAN2/BoxSize""".format(maskrad3*options.apix)

if __name__ == "__main__":
	main()
