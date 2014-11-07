#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine


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
import re, os
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <micrograph1, micrograph2....>
	Use this program to import and filter raw micrographs. If you choose to filter and/or convert format, this program will process each micrograph
	and dump them into the directory './micrographs.', otherwise the micrographs will simply be moved into './micrographs'. If you select the option
	--moverawdata AND you filter or change format, your original micrographs will be moved into the directory './raw_micrographs' and your
	filtered micrographs will be in './micrographs as usual. BDB files are not moved, but they can be processed."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="micrographs",help="List the micrographs to filter here.", default="", guitype='filebox', browser="EMRawDataTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, mode='filter')
	parser.add_pos_argument(name="import_files",help="List the files to import/filter here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, mode='import')
	parser.add_header(name="filterheader", help='Options below this label are specific to filtering', title="### filtering options ###", row=1, col=0, rowspan=1, colspan=2, mode='import,filter')
	parser.add_argument("--invert",action="store_true",help="Invert contrast",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--edgenorm",action="store_true",help="Edge normalize",default=False, guitype='boolbox', row=2, col=1, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--usefoldername",action="store_true",help="If you have the same image filename in multiple folders, and need to import into the same project, this will prepend the folder name on each image name",default=False,guitype='boolbox',row=2, col=2, rowspan=1, colspan=1, mode="import[False]")
	parser.add_argument("--xraypixel",action="store_true",help="Filter X-ray pixels",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--ctfest",action="store_true",help="Estimate defocus from whole micrograph",default=False, guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--astigmatism",action="store_true",help="Includes astigmatism in automatic fitting",default=False, guitype='boolbox', row=3, col=2, rowspan=1, colspan=1, mode='filter[False]')
	parser.add_argument("--moverawdata",action="store_true",help="Move raw data to directory ./raw_micrographs after filtration",default=False)
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=None, guitype='floatbox', row=5, col=0, rowspan=1, colspan=1, mode="filter['self.pm().getAPIX()']")
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=None, guitype='floatbox', row=5, col=1, rowspan=1, colspan=1, mode="filter['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=None, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode="filter['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=6, col=1, rowspan=1, colspan=1, mode="filter")
	parser.add_argument("--defocusmin",type=float,help="Minimum autofit defocus",default=0.6, guitype='floatbox', row=8, col=0, rowspan=1, colspan=1, mode="filter[0.6]")
	parser.add_argument("--defocusmax",type=float,help="Maximum autofit defocus",default=4, guitype='floatbox', row=8, col=1, rowspan=1, colspan=1, mode='filter[4.0]')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	microdir = os.path.join(".","micrographs")
	if not os.access(microdir, os.R_OK):
		os.mkdir("micrographs")

	logid=E2init(sys.argv,options.ppid)

	# After filtration we move micrographs to a directory 'raw_micrographs', if desired
	if options.moverawdata:
		originalsdir = os.path.join(".","raw_micrographs")
		if not os.access(originalsdir, os.R_OK):
			os.mkdir("raw_micrographs")
			
	for i,arg in enumerate(args):
		base = base_name(arg,nodir=not options.usefoldername)
		output = os.path.join(os.path.join(".","micrographs"),base+".hdf")
		cmd = "e2proc2d.py %s %s --inplace"%(arg,output)

		if options.invert: cmd += " --mult=-1"
		if options.edgenorm: cmd += " --process=normalize.edgemean"
		if options.xraypixel: cmd += " --process=threshold.clampminmax.nsigma:nsigma=4"
		
		launch_childprocess(cmd)
		if options.moverawdata:
			os.rename(arg,os.path.join(originalsdir,os.path.basename(arg)))
			
		# We estimate the defocus and B-factor (no astigmatism) from the micrograph and store it in info and the header
		if options.ctfest :
			d=EMData(output,0)
			if d["nx"]<1200 or d["ny"]<1200 : 
				print "CTF estimation will only work with images at least 1200x1200 in size"
				sys.exit(1)
			import e2ctf
			
			ds=1.0/(options.apix*512)
			ffta=None
			nbx=0
			for x in range(100,d["nx"]-512,512):
				for y in range(100,d["ny"]-512,512):
					clip=d.get_clip(Region(x,y,512,512))
					clip.process_inplace("normalize.edgemean")
					fft=clip.do_fft()
					fft.ri2inten()
					if ffta==None: ffta=fft
					else: ffta+=fft
					nbx+=1

			ffta.mult(1.0/(nbx*512**2))
			ffta.process_inplace("math.sqrt")
			ffta["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it

			fftbg=ffta.process("math.nonconvex")
			fft1d=ffta.calc_radial_dist(ffta.get_ysize()/2,0.0,1.0,1)	# note that this handles the ri2inten averages properly

			# Compute 1-D curve and background
			bg_1d=e2ctf.low_bg_curve(fft1d,ds)

			#initial fit, background adjustment, refine fit, final background adjustment
			ctf=e2ctf.ctf_fit(fft1d,bg_1d,bg_1d,ffta,fftbg,options.voltage,options.cs,options.ac,options.apix,1,dfhint=(options.defocusmin,options.defocusmax))
			bgAdj(ctf,fft1d)
			ctf=e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ac,options.apix,1,dfhint=(options.defocusmin,options.defocusmax))
			bgAdj(ctf,fft1d)
			
			if options.astigmatism : e2ctf.ctf_fit_stig(ffta,fftbg,ctf)
			
			#ctf.background=bg_1d
			#ctf.dsbg=ds
			db=js_open_dict(info_name(arg,nodir=not options.usefoldername))
			db["ctf_frame"]=[512,ctf,(256,256),set(),5,1]
			print info_name(arg,nodir=not options.usefoldername),ctf

		E2progress(logid,(float(i)/float(len(args))))

	E2end(logid)

def bgAdj(ctf,fg_1d):
	"""Smooths the background based on the values of the foreground near the CTF zeroes and puts the
	smoothed background into the CTF object"""
	ds=ctf.dsbg
	ctf=ctf
	bg_1d=list(ctf.background)

	xyd=XYData()

	# Find the minimum value near the origin, which we'll use as a zero (though it likely should not be)
	mv=(fg_1d[1],1)
	fz=int(ctf.zero(0)/(ds*2))
	for lz in xrange(1,fz):
		mv=min(mv,(fg_1d[lz],lz))

	xyd.insort(mv[1],mv[0])

	# now we add all of the zero locations to our XYData object
	for i in xrange(100):
		z=int(ctf.zero(i)/ds)
		if z>=len(bg_1d)-1: break
		if fg_1d[z-1]<fg_1d[z] and fg_1d[z-1]<fg_1d[z+1]: mv=(z-1,fg_1d[z-1])
		elif fg_1d[z]<fg_1d[z+1] : mv=(z,fg_1d[z])
		else : mv=(z+1,fg_1d[z+1])
		xyd.insort(mv[0],mv[1])

	# new background is interpolated XYData
	ctf.background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]

	# if our first point (between the origin and the first 0) is too high, we readjust it once
	bs=[fg_1d[i]-ctf.background[i] for i in xrange(fz)]
	if min(bs)<0 :
		mv=(bs[0],fg_1d[0],0)
		for i in xrange(1,fz): mv=min(mv,(bs[i],fg_1d[i],i))
		xyd.set_x(0,mv[2])
		xyd.set_y(0,mv[1])
		
		ctf.background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]


if __name__ == "__main__":
	main()
