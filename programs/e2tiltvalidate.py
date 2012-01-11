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

from EMAN2 import *
import os, math
from EMAN2db import EMTask

def main():
	"""Program to validate a reconstruction by the Richard Henderson tilt validation method. A volume to validate, a small stack (~100 imgs) of untilted and ~10-15 degree
	tilted particles must be presented. The untilted and tilted particle stack must have a one-to-one relationship. In the contour plot, the Tiltaxis is along positive 'Y'
	The tiltaxis angle can be determined from e2RCTboxer.py uisng PairPicker mode. For example, if the tiltaxis is 45 degrees and the tilt angle is -15 degrees, there should
	be a peak in the -X, -Y quadrant at 225 degrees at a magnitude of 15.
	For more details see:
	Optiomal Determination of Particle Orientation, Absolute Hand, and COntrast Loss in Single-particle Electron Cryomicroscopy. Rosenthal, P.B., and Henderson, R. JMB, 333 (2003) pg 721-745
	"""
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	Tiltvalidation using Richard Henderson's technique. To use a stack of untilted and tiltimages whose set relationship is one-to-one is required along with a
	volume to validate. This can be generated using e2RCTboxer.py. After running this program two bits of data are products. A contour plot similar to Figure 5 in the Henderson paper(see below), and a list of
	titlangles and tiltaxes between particle paris, which can be used to makes plot similar to Figure 6 in Hendersons paper. The contour plot is stored as contour.hdf and the tiltpairs data is
	stored as bdb:perparticletilts.
	For more information see:
	Optimal determination of particle orientation, absolute hand, and contrast loss in 
	single-particle electron cryomicroscopy.
	Rosenthal PB, Henderson R.
	J Mol Biol. 2003 Oct 31;333(4):721-45 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	# options associated with e2tiltvalidate.py
	parser.add_header(name="tvheader", help='Options below this label are specific to e2tiltvalidate', title="### e2tiltvalidate options ###", row=3, col=0, rowspan=1, colspan=2, mode="analysis,gui")
	parser.add_argument("--path", type=str,help="The folder the results are placed", default="", guitype='dirbox', dirbasename='TiltValidate', row=0, col=0,rowspan=1, colspan=2, mode="gui")
	parser.add_argument("--volume", type=str,help="3D volume to validate",default=None, guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', row=2, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--untiltdata", type=str,help="Stack of untilted images",default=None, guitype='filebox', browser='EMSetsTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--tiltdata", type=str,help="Stack of tilted images",default=None, guitype='filebox', browser='EMSetsTable(withmodal=True,multiselect=False)', row=1, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--align", type=str,help="The name of a aligner to be used in comparing the aligned images",default="translational", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=6, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--cmp", type=str,help="The name of a 'cmp' to be used in comparing the aligned images",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=7, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--tiltrange", type=int,help="The angular tiltrange to search",default=15, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="analysis")
	parser.add_argument("--sym",  type=str,help="The recon symmetry", default="c1", guitype='symbox', row=5, col=1, rowspan=1, colspan=1, mode="analysis")
	parser.add_argument("--planethres", type=float, help="Maximum out of plane threshold for the tiltaxis. 0 = perfectly in plane, 1 = normal to plane", default=0.1, guitype='floatbox', row=5, col=0, rowspan=1, mode="gui")
	parser.add_argument("--datalabels", action="store_true",help="Add data labels to the plot", default=False, guitype='boolbox', row=6, col=0, rowspan=1, mode="gui")
	parser.add_argument("--maxtiltangle", type=float, help="Maximum tiltangle permitted when finding tilt distances", default=180.0, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="analysis")
	parser.add_argument("--gui",action="store_true",help="Start the GUI for viewing the tiltvalidate plots",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode="gui[True]")
	parser.add_argument("--radcut", type = float, default=-1, help="For use in the GUI, truncate the polar plot after R. -1 = no truncation", guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="gui")
	parser.add_argument("--quaternion",action="store_true",help="Use Quaterions for tilt distance computation",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='analysis')
	parser.add_argument("--eulerfile",type=str,help="Euler angles file, to create tiltdistance from pre-aligned particles. Format is: imgnum, name, az, alt, phi",default=None)
	# options associated with e2projector3d.py
	parser.add_header(name="projheader", help='Options below this label are specific to e2project', title="### e2project options ###", row=9, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--delta", type=float,help="The angular step size for alingment", default=20.0, guitype='floatbox', row=10, col=0, rowspan=1, colspan=1, mode="analysis")
	# options associated with e2simmx.py
	parser.add_header(name="simmxheader", help='Options below this label are specific to e2simmx', title="### e2simmx options ###", row=11, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate)", default="rotate_translate", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=14, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--simaligncmp",type=str,help="Name of the aligner along with its construction arguments (default=ccc)",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=15, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images (default=ccc)", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=13, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--simralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine\', 0)', row=16, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--simraligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=17, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Defulat = 0, no shrinking", guitype='intbox', row=12, col=0, rowspan=1, colspan=1, mode="analysis")
	
	parser.add_argument("--parallel",type=str,help="Parallelism string",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="analysis")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()

	# Run the GUI if in GUI mode
	if options.gui:
		display_validation_plots(options.path, options.radcut, options.planethres, plotdatalabels=options.datalabels)
		exit(0)
		
	if not (options.volume or options.eulerfile):
		print "Error a volume to validate must be presented"
		exit(1)
		
	if not (options.tiltdata or options.eulerfile):
		print "Error a stack of tilted images must be presented"
		exit(1)
		
	if not (options.untiltdata or options.eulerfile):
		print "Error a stack of untiled images must be presented"
		exit(1)
	
	logid=E2init(sys.argv,options.ppid)
	
	options.cmp=parsemodopt(options.cmp)
	options.align=parsemodopt(options.align)
	
	# Make a new dir for each run
	if not options.path : options.path=numbered_path("TiltValidate",True)
	
	#Make tilt distance generator
	tiltgenerator = ComputeTilts(options)
	
	# Compute tilt distances from file if desired. 
	if options.eulerfile:
		# Format is:
		# untilt_imgnum name az alt phi
		# tilt_imgnum name az alt phi
		eulerfile = open(options.eulerfile,"r")
		eulers = eulerfile.readlines()
		eulerfile.close()
		untilteulerlist = []
		tilteulerlist = []
		for i, euler in enumerate(eulers):
			fields = euler.split()
			if i % 2:
				tilteulerlist.append({'alt':float(fields[2]),'az':float(fields[3]),'phi':float(fields[4])})
			else:
				untilteulerlist.append({'alt':float(fields[2]),'az':float(fields[3]),'phi':float(fields[4])})
		tiltgenerator.findtilts_fromeulers(untilteulerlist, tilteulerlist)
		exit(0)

	# Initialize parallelism if being used
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
	else:
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer("thread:1")
		#etc.precache(pclist)
	
	# Otherwise compute tilt distances from data
	#Read in the images
	tiltimgs = EMData.read_images(options.tiltdata)
	untiltimgs = EMData.read_images(options.untiltdata)
	if len(tiltimgs) != len(untiltimgs):
		print "The untilted image stack is not the same lenght as the tilted stack!!!"
		exit(1)
		
	# Do projections
	e2projectcmd = "e2project3d.py %s --orientgen=eman:delta=%f:inc_mirror=1:perturb=0 --outfile=bdb:%s#projections --projector=standard --sym=%s" % (options.volume,options.delta,options.path, options.sym) # Seems to work better when I check all possibilites
	if options.parallel: e2projectcmd += " --parallel=%s" %options.parallel
	run(e2projectcmd)
	
	# Make simmx
	e2simmxcmd = "e2simmx.py bdb:%s#projections %s bdb:%s#simmx -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d" % (options.path,options.untiltdata,options.path,options.simcmp,options.simalign,options.simaligncmp,options.verbose)
	if options.simralign: e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
	if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
	if options.shrink: e2simmxcmd += " --shrink=%d" %options.shrink
	run(e2simmxcmd)
	
	e2simmxcmd = "e2simmx.py bdb:%s#projections %s bdb:%s#simmx_tilt -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d" % (options.path,options.tiltdata,options.path,options.simcmp,options.simalign,options.simaligncmp,options.verbose)
	if options.simralign: e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
	if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
	if options.shrink: e2simmxcmd += " --shrink=%d" %options.shrink
	run(e2simmxcmd)

	# Read in the data
	simmx= EMData.read_images("bdb:%s#simmx"%options.path)
	simmx_tilt= EMData.read_images("bdb:%s#simmx_tilt"%options.path)
	projections = EMData.read_images("bdb:%s#projections"%options.path)
	volume = EMData() 
	volume.read_image(options.volume) # I don't know why I cant EMData.read_image.......
	
	# Generate tilts from data
	tiltgenerator.findtilts_fromdata(simmx, simmx_tilt, projections, volume, untiltimgs, tiltimgs) 
	exit(1)

	# Make contour plot to validate each particle
	tasks=[]
	distplot = EMData(options.tiltrange*2+1,options.tiltrange*2+1)
	distplot.to_zero()
	for imgnum in xrange(simmx[0].get_ysize()):
		bestscore = float('inf')
		bestrefnum = 0
		for refnum in xrange(simmx[0].get_xsize()):
			if simmx[0].get_value_at(refnum, imgnum) < bestscore:
				bestscore = simmx[0].get_value_at(refnum, imgnum)
				bestrefnum = refnum
		# Get the euler angle for this particle and call compare to tilt"bdb:%s#
		euler_xform = projections[bestrefnum].get_attr('xform.projection')
		tasks.append(CompareToTiltTask(volume, tiltimgs[imgnum], imgnum, euler_xform, simmx[3].get_value_at(bestrefnum, imgnum), distplot, options.tiltrange, 1, options))
		
	# Farm out the work and hang till finished!
	tids=etc.send_tasks(tasks)
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tids)
		tids=[j for i,j in enumerate(tids) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if len(tids)==0: break
		
	
	# Make scoremx avg
	scoremxs = EMData.read_images("bdb:%s#scorematrix"%options.path)
	avgmxavger = Averagers.get('mean')
	for mx in scoremxs:
		avgmxavger.add_image(mx)
	avgmx = avgmxavger.finish()
	avgmx.write_image("%s/contour.hdf"%options.path)
	distplot.write_image("%s/distplot.hdf"%options.path)
	
	E2end(logid)

class ComputeTilts:
	def __init__(self, options):
		self.options = options
		self.symmeties = Symmetries.get(self.options.sym)
		self.particletilt_list = []
		self.tdb = db_open_dict("bdb:%s#perparticletilts"%self.options.path)
		self.test = open("test.dat","w")
		
	def findtilts_fromdata(self, simmx, simmx_tilt, projections, volume, untiltimgs, tiltimgs):
		""" Compute tiltdistances based on data """
		for imgnum in xrange(simmx[0].get_ysize()):
			untiltbestscore = float('inf')
			tiltbestscore = float('inf')
			untiltbestrefnum = 0
			tiltbestrefnum = 0
			tiltimgnum = imgnum
			for refnum in xrange(simmx[0].get_xsize()):
				if simmx[0].get_value_at(refnum, imgnum) < untiltbestscore:
					untiltbestscore = simmx[0].get_value_at(refnum, imgnum)
					untiltbestrefnum = refnum		
				if simmx_tilt[0].get_value_at(refnum, tiltimgnum) < tiltbestscore:
					tiltbestscore = simmx_tilt[0].get_value_at(refnum, tiltimgnum)
					tiltbestrefnum = refnum
			# Untilt
			untilt_euler_xform = Transform({"type":"eman","phi":-simmx[3].get_value_at(untiltbestrefnum, imgnum)})*projections[untiltbestrefnum].get_attr('xform.projection')
			untiltrot = untilt_euler_xform.get_rotation("eman")
			untilt_euler_xform.set_rotation({"type":"eman","az":untiltrot["az"],"alt":untiltrot["alt"],"phi":-simmx[3].get_value_at(untiltbestrefnum, imgnum)})
			# Tilt
			tilt_euler_xform = Transform({"type":"eman","phi":-simmx_tilt[3].get_value_at(tiltbestrefnum, tiltimgnum)})*projections[tiltbestrefnum].get_attr('xform.projection')
			tiltrot = tilt_euler_xform.get_rotation("eman")
			tilt_euler_xform.set_rotation({"type":"eman","az":tiltrot["az"],"alt":tiltrot["alt"],"phi":-simmx_tilt[3].get_value_at(tiltbestrefnum, tiltimgnum)})
			
			#Find best solultion takeing sym into accout
			tiltpars = self.find_bestsymsoln(imgnum, untilt_euler_xform, tilt_euler_xform)
			
			# Write out test results
			volume.project('standard', untilt_euler_xform).write_image("projections_ut.hdf", imgnum)
			volume.project('standard', tilt_euler_xform).write_image("projections_t.hdf", imgnum)
			untiltimgs[imgnum].set_attr('tiltangle',round(tiltpars[0],2))
			untiltimgs[imgnum].set_attr('tiltaxis',round(tiltpars[1],2))
			untiltimgs[imgnum].write_image("untiltaligned.hdf", imgnum)
			tiltimgs[imgnum].set_attr('tiltangle',round(tiltpars[0],2))
			tiltimgs[imgnum].set_attr('tiltaxis',round(tiltpars[1],2))
			tiltimgs[imgnum].write_image("tiltaligned.hdf", imgnum)
			
		self.tdb["junk"] = 1 # The Stupid DB doesn't write my list unless I don't this!!	
		self.finish()
	
	def findtilts_fromeulers(self, untilteulers, tilteulers):
		"""
		untilteulers and tilteulers are a list of dicts, each dict contyaining az, alt and phi
		"""
		if len(untilteulers) != len(tilteulers):
			raise ValueError("untilteulers and tilteulers list is not the same length")
		
		for num in xrange(len(untilteulers)):
			untilt_euler_xform = Transform({'type':'eman','az':untilteulers[num]['az'],'alt':untilteulers[num]['alt'],'phi':untilteulers[num]['phi']})
			tilt_euler_xform = Transform({'type':'eman','az':tilteulers[num]['az'],'alt':tilteulers[num]['alt'],'phi':tilteulers[num]['phi']})
			self.find_bestsymsoln(num, untilt_euler_xform, tilt_euler_xform)
			
		self.tdb["junk"] = 1 # The Stupid DB doesn't write my list unless I don't this!!
		self.finish()
					
	def find_bestsymsoln(self, imgnum, untilt_euler_xform, tilt_euler_xform):
		# Compute tiltxis and tiltangle for each particel pair. For symmetric objects the tilttransform is selected as the one that has a tiltaxis
		# closest to perpidicular to the imaging plane (this is how Richard handles it). 
		anglefound = False
		if self.options.quaternion:
			bestinplane = 1.0
			for sym in self.symmeties.get_syms():
				tiltxform = tilt_euler_xform*sym.inverse()*untilt_euler_xform.inverse()
				# Select the symmetry solution whose tiltaxis is in plane
				if math.fabs(tiltxform.get_rotation("spin")["n3"]) < bestinplane:
					if tiltxform.get_rotation("spin")["Omega"] < self.options.maxtiltangle:
						bestinplane = math.fabs(tiltxform.get_rotation("spin")["n3"])
						besttiltangle = tiltxform.get_rotation("spin")["Omega"]
						besttiltaxis = math.degrees(math.atan2(tiltxform.get_rotation("spin")["n2"],tiltxform.get_rotation("spin")["n1"]))
						anglefound = True
					
				self.test.write("\t%f %f %f %f\n"%(tiltxform.get_rotation("spin")["Omega"],tiltxform.get_rotation("spin")["n1"],tiltxform.get_rotation("spin")["n2"],tiltxform.get_rotation("spin")["n3"]))
				self.test.write("\t%f %f %f\n"%(tiltxform.get_rotation("eman")["az"],tiltxform.get_rotation("eman")["alt"],tiltxform.get_rotation("eman")["phi"]))
		else:
			bestinplane = 180.0
			for sym in self.symmeties.get_syms():
				tiltxform = tilt_euler_xform*sym.inverse()*untilt_euler_xform.inverse()
				if math.fabs(tiltxform.get_rotation("eman")['az'] - (-tiltxform.get_rotation("eman")['phi'] % 360)) < bestinplane:
					if tiltxform.get_rotation("eman")["alt"] < self.options.maxtiltangle:
						bestinplane = math.fabs(tiltxform.get_rotation("eman")['az'] - (-tiltxform.get_rotation("eman")['phi'] % 360))
						besttiltangle = tiltxform.get_rotation("eman")["alt"]
						besttiltaxis = (tiltxform.get_rotation("eman")['az'] + (-tiltxform.get_rotation("eman")['phi'] % 360))/2.0
						anglefound = True
				
				self.test.write("\t%f %f %f %f\n"%(tiltxform.get_rotation("spin")["Omega"],tiltxform.get_rotation("spin")["n1"],tiltxform.get_rotation("spin")["n2"],tiltxform.get_rotation("spin")["n3"]))
				self.test.write("\t%f %f %f Ip %f\n"%(tiltxform.get_rotation("eman")["az"],tiltxform.get_rotation("eman")["alt"],tiltxform.get_rotation("eman")["phi"], bestinplane))
			
		if anglefound:
			self.test.write("The best angle is %f with a tiltaxis of %f\n"%(besttiltangle,besttiltaxis))
			self.particletilt_list.append([imgnum, besttiltangle,besttiltaxis,bestinplane])
			return [besttiltangle,besttiltaxis,bestinplane]
		else:
			return [0.0, 0.0, 0.0]
		
	def finish(self):
		self.test.close()
		self.tdb["particletilt_list"] = self.particletilt_list
		self.tdb.close()
	
# Function to compute a similarity map for a given image
class CompareToTiltTask(EMTask):
	""" A parallelized version to compute contout plots """
	def __init__(self, volume, tilted, imgnum, eulerxform, zrot, distplot, tiltrange, tiltstep, options):
		if options.shrink:
			shrunkvol = volume.process("math.meanshrink",{"n":options.shrink})
			shrunktilted = tilted.process("math.meanshrink",{"n":options.shrink})
			data = {"volume":shrunkvol,"tilted":shrunktilted}
		else:
			data = {"volume":volume,"tilted":tilted}
		EMTask.__init__(self,"CmpTilt",data,options,"")
		
		self.imgnum = imgnum
		self.eulerxform=eulerxform
		self.zrot=zrot
		self.distplot=distplot
		self.tiltrange=tiltrange
		self.tiltstep=tiltstep
		
	def execute(self,callback=None):
		scoremx = EMData(2*self.tiltrange+1,2*self.tiltrange+1)
		bestscore = float('inf')
		for rotx in xrange(-self.tiltrange, self.tiltrange+1,self.tiltstep):
			for roty in xrange(-self.tiltrange, self.tiltrange+1,self.tiltstep):
				# First make the projection
				tiltangle = math.sqrt(rotx*rotx + roty*roty)
				tiltaxis = math.degrees(math.atan2(roty, rotx))
				tiltxform = Transform({"type":"eman","az":tiltaxis,"alt":tiltangle,"phi":-tiltaxis})
				inplane = Transform({"type":"eman", "phi":-self.zrot})
				totalxform = tiltxform*inplane*self.eulerxform
				testprojection = self.data['volume'].project("standard",totalxform)
				tiltalign = self.data['tilted'].align(self.options.align[0],testprojection,self.options.align[1],self.options.cmp[0],self.options.cmp[1])
				score = tiltalign.cmp(self.options.cmp[0], testprojection, self.options.cmp[1])
				scoremx.set_value_at(rotx+self.tiltrange, roty+self.tiltrange, score)
		scoremx.write_image("bdb:%s#scorematrix"%self.options.path, self.imgnum)
		# Denoise the contiur plot, I need to experiment around with this
		radius = 4
		scoremx_blur = scoremx.process('eman1.filter.median',{'radius':radius})
		scoremx_blur = scoremx_blur.get_clip(Region(radius, radius, scoremx_blur.get_xsize() - radius*2, scoremx_blur.get_ysize() - radius*2))
		# Find the peak
		maxpeak = scoremx_blur.calc_min_location()
		self.distplot.set_value_at(maxpeak[0], maxpeak[1], self.distplot.get_value_at(maxpeak[0], maxpeak[1]) + 1.0)

def run(command):
	"Execute a command with optional verbose output"		    
	print command
	error = launch_childprocess(command)
	if error==11 :
		pass		    
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command		    
		exit(1)

def display_validation_plots(path, radcut, planethres, plotdatalabels=False):
	from emplot2d import EMPolarPlot2DWidget
	from emimage2d import EMImage2DWidget
	from emapplication import EMApp
	r = []
	theta = []
	datap = []
	try:
		tpdb = db_open_dict("bdb:%s#perparticletilts"%path)
		tplist = tpdb["particletilt_list"]
		for tp in tplist:
			if tp[3] > planethres:	# if the out of plane threshold is too much
				continue
			if plotdatalabels: datap.append(tp[0])
			r.append(tp[1])
			theta.append(math.radians(tp[2]))
		tpdb.close()
	except:
		print "Couldn't load tp from DB, not showing polar plot"
	data = None	
	try:
		data = EMData("%s/contour.hdf"%path)
	except:
		print "Couldn't open contour plot"
	
	if not data and not (theta and r): return
	app = EMApp()
	if theta and r:
		plot = EMPolarPlot2DWidget()
		plot.set_data((theta,r),linewidth=50,radcut=radcut,datapoints=datap)
		plot.show()
	if data:
		image = EMImage2DWidget(data)
		image.show()
	app.exec_()
		
if __name__ == "__main__":
	main()