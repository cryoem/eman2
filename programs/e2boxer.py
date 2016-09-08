#!/usr/bin/env python
#
# Author: Steven Ludtke 2014/04/27
# Copyright (c) 2014- Baylor College of Medicine


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
from EMAN2jsondb import *
import Queue
import os,sys

class nothing:
	def __init__(self,x=None,y=None,z=None,q=None):
		return

try: 
	from PyQt4 import QtCore, QtGui
	from PyQt4.QtCore import Qt
	from emimage2d import EMImage2DWidget
	from emplot2d import EMPlot2DWidget
	from emimagemx import EMImageMXWidget
	from valslider import ValSlider,CheckBox,ValBox
	from emshape import EMShape
except:
	QtGui=nothing()
	QtCore=nothing()
	QtCore.QObject=nothing()
	QtGui.QWidget=nothing()


# ok, this is kind of bad style, but really don't want to have to drag this flag around through many objects
invert_on_read=False

def load_micrograph(filename):
	if "\t" in filename: filename=filename.split()[1]
	fsp2="micrographs/"+filename
	if os.path.exists(fsp2) : filename=fsp2
	
	n=EMUtil.get_image_count(filename)
	if n==0 :
		QtGui.QMessageBox.warning(None,"Error","The file {} contains no images".format(newfilename))
		return
	elif n==1 :
		img=EMData(filename,0)		# single image
	else :
		img=EMData(filename,0)		# movie stack (we assume)
		for i in xrange(1,n):
			im=EMData(filename,i)
			img.add(im)
		img.mult(1.0/n)
		
	if invert_on_read : img.mult(-1.0)
	return img

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image> <image2>....

	The even newer version of e2boxer. Complete rewrite. Incomplete.
	
	This program 
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="micrographs",help="List the file to process with e2boxer here.", default="", guitype='filebox', browser="EMBoxesTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--allmicrographs",action="store_true",default=False,help="Add all images from micrographs folder", guitype='boolbox', row=10, col=0, rowspan=1, colspan=1, mode="boxing,extraction")
	parser.add_argument("--unboxedonly",action="store_true",default=False,help="Only include image files without existing box locations", guitype='boolbox', row=10, col=1, rowspan=1, colspan=1, mode="boxing,extraction")
	parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=2, col=0, rowspan=1, colspan=1, mode="boxing,extraction")
	parser.add_argument("--ptclsize","-P",type=int,help="Longest axis of particle in pixels (diameter, not radius)",default=-1, guitype='intbox', row=2, col=1, rowspan=1, colspan=1, mode="boxing,extraction")
	parser.add_argument("--write_dbbox",action="store_true",default=False,help="Export EMAN1 .box files",guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="extraction")
	parser.add_argument("--write_ptcls",action="store_true",default=False,help="Extract selected particles from micrographs and write to disk", guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode="extraction[True]")
	parser.add_argument("--invert",action="store_true",help="If specified, inverts input contrast. Particles MUST be white on a darker background.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="extraction")
	parser.add_argument("--no_ctf",action="store_true",default=False,help="Disable CTF determination", guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode="extraction, boxing")
	
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=-1, guitype='floatbox', row=14, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getAPIX()'],boxing,extraction")
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=-1, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="autofit['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=-1, guitype='floatbox', row=5, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=5, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--autopick",type=str,default=None,help="Perform automatic particle picking. Provide mode and parameter string, eg - auto_local:threshold=5.5")
	parser.add_argument("--gui", action="store_true", default=False, help="Interactive GUI mode", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="boxing[True]")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=14, col=1, rowspan=1, colspan=1,mode="boxing")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	
	global invert_on_read
	if options.invert : invert_on_read = True

	if options.allmicrographs :
		if len(args)>0 : print "Specified micrograph list replaced with contents of micrographs/"
		args=sorted(["micrographs/"+i for i in os.listdir("micrographs") if not i.startswith('.')])
	#else: args=[i.replace("micrographs/","") for i in args]

	oargs=args
	args=[]
	basenames=[]
	for f in oargs:
		bname=base_name(f)
		if bname in basenames:
			#### so we do not box multiple times on different versions of the same micrograph
			if options.verbose: "skipping {}".format(f)
			continue
		basenames.append(bname)
		db=js_open_dict(info_name(f))
		try: boxes=db["boxes"]
		except: boxes=[]
		if not options.unboxedonly or len(boxes)==0 : args.append("{}\t{}".format(len(boxes),f))
		db.close()

	#####
	# Parameter Validation
	project_db = js_open_dict("info/project.json")

	if not (options.gui or options.write_ptcls or options.write_dbbox or options.autopick):
		print "Error: No actions specified. Try --gui for interactive/semi-automated particle picking." 

	# Some of this seems redundant, it is to insure self-consistency
	if options.boxsize<2:
		try: 
			options.boxsize = project_db["global.boxsize"]
		except:
			print "Warning: No box size specified, and no box size in project info. Please specify."
			options.boxsize=project_db.setdefault("global.boxsize",128)
			#sys.exit(1)
			
	if good_size(options.boxsize)!=options.boxsize :
		print "Bad box size detected. Adjusting size to {}. See http://eman2.org/emanwiki/EMAN2/BoxSize".format(good_size(options.boxsize))
		options.boxsize=good_size(options.boxsize)
	project_db["global.boxsize"]=options.boxsize
	boxsize=options.boxsize
	boxsize2=boxsize/2

	if options.ptclsize<2:
		try: options.ptclsize=project_db["global.ptclsize"]
		except:
			print "Warning: No particle size specified. None found in project DB. Please specify approximate maximum particle dimension in pixels."
			options.ptclsize=project_db.setdefault("global.ptclsize",64)
			#sys.exit(1)
			
	if options.ptclsize>boxsize*0.8:
		print "Warning: Invalid particle size detected. Box size should normally be 1.5 - 2x particle size, and must be at least 1.2x particle size." 
		#sys.exit(1)
		
	project_db["global.ptclsize"]=options.ptclsize

	# CTF related options
	if not options.no_ctf :
		if options.voltage>1500 :
			options.voltage/=1000
			print "Voltage specified in kV. Adjusting specified value to ",options.voltage
		if options.voltage<10 :
			try: 
				options.voltage=float(project_db["global.microscope_voltage"])
				print "Using project voltage of ",options.voltage,"kV"
			except:
				print "Error: No voltage specified, and no project settings available. Disabling CTF mode."
				options.no_ctf=True
		if options.cs<0 :
			try:
				options.cs=float(project_db["global.microscope_cs"])
				print "Using project Cs of ",options.cs,"mm"
			except:
				print "Error: No Cs specified, and no project settings available. Disabling CTF mode."
				options.no_ctf=True
		if options.ac<0 and not options.no_ctf:
			print "Error: Invalid %AC value. Disabling CTF mode."
			options.no_ctf=True
		if options.ac<1.0 :
			print "Warning: %AC should be specified as a %. If you intended a %AC>1%, please try again. Will proceed with the specified value"

	if options.apix<=0 :
		try:
			options.apix=float(project_db["global.apix"])
			print "Warning: No A/pix specified. Using ",options.apix," from project. Please insure this is correct for the images being boxed!"
		except:
			print "Error: Value required for A/pixel. If this is a non TEM image, suggest --apix=1 and --no_ctf."
			sys.exit(1)

	if options.ptclsize*1.5>boxsize :
		print "WARNING: Strongly recommend using a box size 1.5 - 2.0x the maximum dimension of the particle! This may be pushed to ~1.25x in some cases, but results may be suboptimal."
		print "Your box size is {:1.2f}x the particle size. Recommend a size of at least {:d}".format(boxsize/float(options.ptclsize),good_size(int(options.ptclsize*1.5)))
		
		
	logid=E2init(sys.argv,options.ppid)


	if options.autopick!=None :
		apick=parsemodopt(options.autopick)
		for s,k,pcl in aboxmodes:
			if k==apick[0] : break
		
		if os.path.exists("info/boxrefs.hdf"):
			goodrefs=EMData.read_images("info/boxrefs.hdf")
		else: goodrefs=[]
		
		if os.path.exists("info/boxrefsbad.hdf"):
			badrefs=EMData.read_images("info/boxrefsbad.hdf")
		else: badrefs=[]
		
		for i,fspi in enumerate(args):
			fsp=fspi.split()[1]
			micrograph=load_micrograph(fsp)

			newboxes=pcl.do_autobox(micrograph,goodrefs,badrefs,options.apix,options.threads,apick[1],None)
			print "{}) {} boxes -> {}".format(i,len(newboxes),fsp)
			
			# if we got nothing, we just leave the current results alone
			if len(newboxes)==0 : continue
		
			# read the existing box list and update
			db=js_open_dict(info_name(fsp))
			try: 
				boxes=db["boxes"]
				# Filter out all existing boxes for this picking mode
				bname=newboxes[0][2]
				boxes=[b for b in boxes if b[2]!=bname]
			except:
				boxes=[]
				
			boxes.extend(newboxes)
			
			db["boxes"]=boxes
			db.close()

	if options.gui :
		if isinstance(QtGui,nothing) :
			print "====================================="
			print "ERROR: GUI mode unavailable without PyQt4"
			sys.exit(1)
		from emapplication import EMApp
		app=EMApp()
		gui=GUIBoxer(args,options.voltage,options.apix,options.cs,options.ac,options.boxsize,options.ptclsize,options.threads)
		gui.show()
		app.exec_()

	if options.write_dbbox:
		write_boxfiles(args,boxsize)
		print ".box files written to boxfiles/"

	if options.write_ptcls:
		write_particles(args,boxsize,options.verbose)
		print "Particles written to particles/*_ptcls.hdf"

	E2end(logid)

def write_boxfiles(files,boxsize):
	"""This function will write a boxfiles/*.box file for each provided micrograph filename based on box locations
	in the corresponding info/*json file."""
	
	try: os.mkdir("boxfiles")
	except: pass
	boxsize2=boxsize/2

	for m in files:
		db=js_open_dict(info_name(m))
		boxes=db.setdefault("boxes",[])
		if len(boxes)==0 : continue
		out=file("boxfiles/{}.box".format(base_name(m)),"w")
		for b in boxes:
			out.write("{:d}\t{:d}\t{:d}\t{:d}\n".format(b[0]-boxsize2,b[1]-boxsize2,boxsize,boxsize))

def write_particles(files,boxsize,verbose):
	"""This function will write a particles/*_ptcls.hdf file for each provided micrograph, based on
	box locations in the corresponding info/*json file. To use this with .box files, they must be imported
	to a JSON file first."""
	
	try: os.mkdir("particles")
	except: pass
	boxsize2=boxsize/2
	
	for nm in files:
		n,m=nm.split()
		base=base_name(m)
		ptcl="particles/{}.hdf".format(base)
		
		# get the list of box locations
		db=js_open_dict(info_name(m))
		boxes=db.setdefault("boxes",[])
		if len(boxes)==0 :
			if verbose :
				print "No particles in ",m
			continue
	
		# remove any existing file
		try: os.unlink(ptcl)
		except: pass
	
		if verbose : print "{} : {} particles written to {}".format(m,len(boxes),ptcl)
		micrograph=load_micrograph(m)		# read micrograph
		for i,b in enumerate(boxes):
			boxim=micrograph.get_clip(Region(b[0]-boxsize2,b[1]-boxsize2,boxsize,boxsize))
			boxim["ptcl_source_coord"]=(b[0],b[1])
			boxim["ptcl_source_image"]=m
			boxim.write_image(ptcl,i)
	
##########
# to add a new autoboxer module, create a class here, then add it to the GUIBoxer.aboxmodes list below
##########

class boxerByRef(QtCore.QObject):
	@staticmethod
	def setup_gui(gridlay):
		boxerByRef.threshold=ValSlider(None,(0.1,8),"Threshold",6.0,90)
		gridlay.addWidget(boxerByRef.threshold,0,0)
	
	@staticmethod
	def do_autobox(micrograph,goodrefs,badrefs,apix,nthreads,params,prog=None):
		# If parameters are provided via params (as if used from command-line) we use those values,
		# if that fails, we check the GUI widgets, which were presumably created in this case
		if len(goodrefs)<1 :
			print 'Box reference images ("Good Refs") required for autopicking'
			return []
		try: threshold=params["threshold"]
		except:
			try: threshold=boxerByRef.threshold.getValue()
			except:
				print "Error, no threshold (0.1-2) specified"
				return
		
		downsample=10.0/apix			# we downsample to 10 A/pix
		microdown=micrograph.process("normalize.edgemean").process("math.fft.resample",{"n":downsample})
		gs=good_size(max(microdown["nx"],microdown["ny"]))
		microf=microdown.get_clip(Region(0,0,gs,gs)).do_fft()
		print "downsample by ",downsample,"  Good size:",gs
	
		## Here we precompute a normalization image to deal with local standard deviation variation
		#nx=goodrefs[0]["nx"]/downsample
		#circle=EMData(gs,gs,1)
		#circle.to_one()
		#circle.process_inplace("mask.sharp",{"outer_radius":nx/2})
##		circle.process_inplace("normalize.unitlen")
		#circle.process_inplace("xform.phaseorigin.tocorner")
		#circlef=circle.do_fft()
		
		#ccfc=microf.calc_ccf(circlef)
		#ccfc.mult(1.0/(nx*nx))
		#ccfc.process_inplace("math.squared")
		
		#md2=microdown.process("math.squared")
		#md2f=md2.get_clip(Region(0,0,gs,gs)).do_fft()
		#norm=md2f.calc_ccf(circlef)
		#norm.mult(1.0/(nx*nx))
	
		## Norm should now be the variance
		#norm.sub(ccfc)
	
		# Iterate over refs
		owner=EMData(gs,gs,1)
		maxav=Averagers.get("minmax",{"max":1,"owner":owner})
		
		# Iterate over in-plane rotation for each ref
		jsd=Queue.Queue(0)
		thrds=[threading.Thread(target=boxerByRef.ccftask,args=(jsd,ref,downsample,gs,microf,ri)) for ri,ref in enumerate(goodrefs)]

		n=-1

		# here we run the threads and save the results, no actual alignment done here
		print len(thrds)," threads"
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) and threading.active_count()!=nthreads:
#				print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: 
				time.sleep(0.05)
				if prog!=None : 
					prog.setValue(prog.value())
					if prog.wasCanceled() : break
		
			while not jsd.empty():
				# add each ccf image to our maxval image as it comes in
				ccf=jsd.get()
				maxav.add_image(ccf)
				sys.stdout.flush()
		print ""

		for t in thrds:
			t.join()

			
		final=maxav.finish()
		# smooth out a few spurious peaks. Hopefully doesn't mess up ownership assignment significantly
		final.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.2})
		# get rid of nasty edges
		final.clip_inplace(Region(0,0,int(micrograph["nx"]/downsample),int(micrograph["ny"]/downsample)))
		owner.clip_inplace(Region(0,0,int(micrograph["nx"]/downsample),int(micrograph["ny"]/downsample)))
		#norm.clip_inplace(Region(0,0,int(micrograph["nx"]/downsample),int(micrograph["ny"]/downsample)))
		#norm.process("math.sqrt")
		#norm.process_inplace("math.reciprocal")
		# this establishes the scale for our threshold (in terms of sigma)
		final.add(-final["mean"]) 
		# Now pull out only local peaks
		# Zero edges to eliminate boxes within 1/2 box size of edge
		edge=int(goodrefs[0]["nx"]/(2.0*downsample)+0.5)
		#final.mult(norm)
		final.process_inplace("mask.zeroedge2d",{"x0":edge,"y0":edge})
		final.process_inplace("mask.onlypeaks",{"npeaks":0,"usemean":0})
		final.process_inplace("normalize.edgemean")
#		final.process_inplace("threshold.belowtozero",{"minval":threshold})

		final.write_image("final.hdf",0)
		owner.write_image("final.hdf",1)
		#norm.write_image("final.hdf",2)
#		display(final)
		
		print "Find peaks"
		# Identify the peaks we want to keep and turn them into rough box locations
		boxes=[]
		locs=final.calc_highest_locations(threshold)
		for i,p in enumerate(locs):
			# loop over all higher valued peaks, and make sure we aren't too close to any of them
			for pc in locs[:i]:
				if hypot(pc.x-p.x,pc.y-p.y)<=edge : break
			else:
				# We only get here if the loop completed
				boxes.append((int(p.x*downsample),int(p.y*downsample),"auto_ref",owner[p.x,p.y]))
		
		print "Refine box locations"
		# Refine the box locations at full sampling
		cmpim=[]
		boxes2=[]
		for box in boxes:
			ownn=int(floor(box[3]))
			owna=(box[3]-ownn)*360.0
			ref=goodrefs[ownn].process("xform",{"transform":Transform({"type":"2d","alpha":owna})})
			ref.process_inplace("normalize.edgemean")
			ptcl=micrograph.get_clip(Region(box[0]-ref["nx"]/2,box[1]-ref["ny"]/2,ref["nx"],ref["ny"]))		# of course, nx == ny anyway
			ali=ref.align("rotate_translate",ptcl)
			ax,ay=ali["xform.align2d"].get_trans_2d()
			boxes2.append((box[0]+ax,box[1]+ay,box[2]))
#			print ownn,owna,ali["xform.align2d"]
			#cmpim.append(ref)
			#cmpim.append(ali)
			#cmpim.append(ptcl)
		#display(cmpim)
			
		print "done"
		
		return boxes2


	@staticmethod
	def ccftask(jsd,ref,downsample,gs,microf,ri):

		mref=ref.process("mask.soft",{"outer_radius":ref["nx"]/2-4,"width":3})
		mref.process_inplace("normalize.unitlen")
		
		for ang in xrange(0,360,10):
			dsref=mref.process("xform",{"transform":Transform({"type":"2d","alpha":ang})})
			# don't downsample until after rotation
			dsref.process_inplace("math.fft.resample",{"n":downsample})
			dsref.process_inplace("normalize")
			diff=(gs-dsref["nx"])/2
			dsref=dsref.get_clip(Region(-diff,-diff,gs,gs))
			dsref.process_inplace("xform.phaseorigin.tocorner")
			ccf=microf.calc_ccf(dsref)
			#ccf.process_inplace("normalize")
			ccf["ortid"]=ri+ang/360.0			# integer portion is projection number, fractional portion is angle, should be enough precision with the ~100 references we're using

			jsd.put(ccf)
		sys.stdout.write("*")
		
class boxerLocal(QtCore.QObject):
	@staticmethod
	def setup_gui(gridlay):
		boxerLocal.threshold=ValSlider(None,(0,8.0),"Threshold",5.0,90)
		gridlay.addWidget(boxerLocal.threshold,0,0)
	
	@staticmethod
	def do_autobox(micrograph,goodrefs,badrefs,apix,nthreads,params,prog=None):
		# If parameters are provided via params (as if used from command-line) we use those values,
		# if that fails, we check the GUI widgets, which were presumably created in this case
		if len(goodrefs)<1 :
			print 'Box reference images ("Good Refs") required for autopicking'
			return []
		try: threshold=params["threshold"]
		except:
			try: threshold=boxerLocal.threshold.getValue()
			except:
				print "Error, no threshold (0.1-2) specified"
				return
		
		nx=goodrefs[0]["nx"]
		downsample=8.0/apix			# we downsample to 10 A/pix
		nxdown=good_size(int(nx/downsample))
		downsample=float(nx)/float(nxdown)
		microdown=micrograph.process("normalize.edgemean").process("math.fft.resample",{"n":downsample})
		print "downsample by ",downsample
		
		# Each thread tries one reference
		owner=EMData(microdown["nx"],microdown["ny"],1)
		maxav=Averagers.get("minmax",{"max":1,"owner":owner})
		
		# This prevents some threading issues
		r=goodrefs[0].process("math.fft.resample",{"n":downsample})
		r.align("rotate_translate",r)
		
		jsd=Queue.Queue(0)
		thrds=[threading.Thread(target=boxerLocal.ccftask,args=(jsd,ref,downsample,microdown,ri)) for ri,ref in enumerate(goodrefs)]

		n=-1

		# here we run the threads and save the results, no actual alignment done here
		print len(thrds)," threads"
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) and threading.active_count()!=nthreads:
#				print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: 
				time.sleep(0.05)
				if prog!=None : 
					prog.setValue(prog.value())
					if prog.wasCanceled() : break
		
			while not jsd.empty():
				# add each ccf image to our maxval image as it comes in
				ccf=jsd.get()
				maxav.add_image(ccf)
				sys.stdout.flush()
		print ""

		for t in thrds:
			t.join()

			
		final=maxav.finish()
		# smooth out a few spurious peaks. Hopefully doesn't mess up ownership assignment significantly
		final.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.2})
		# get rid of nasty edges
#		final.clip_inplace(Region(0,0,int(micrograph["nx"]/downsample),int(micrograph["ny"]/downsample)))
#		owner.clip_inplace(Region(0,0,int(micrograph["nx"]/downsample),int(micrograph["ny"]/downsample)))
		# this establishes the scale for our threshold (in terms of sigma)
#		final.add(-final["mean"]) 
		# Now pull out only local peaks
		# Zero edges to eliminate boxes within 1/2 box size of edge
		edge=int(goodrefs[0]["nx"]/(2.0*downsample)+0.5)
		#final.mult(norm)
		final.process_inplace("mask.zeroedge2d",{"x0":edge,"y0":edge})
		final.process_inplace("mask.onlypeaks",{"npeaks":0,"usemean":0})
		final.mult(1.0/final["sigma_nonzero"])
#		final.process_inplace("normalize.edgemean")
#		final.process_inplace("threshold.belowtozero",{"minval":threshold})

		microdown.write_image("final.hdf",0)
		final.write_image("final.hdf",1)
		owner.write_image("final.hdf",2)
		#norm.write_image("final.hdf",2)
#		display(final)
		
		print "Find peaks"
		# Identify the peaks we want to keep and turn them into rough box locations
		boxes=[]
		locs=final.calc_highest_locations(threshold)

		for i,p in enumerate(locs):
			# loop over all higher valued peaks, and make sure we aren't too close to any of them
			for pc in locs[:i]:
				if hypot(pc.x-p.x,pc.y-p.y)<=edge : break
			else:
				# We only get here if the loop completed
				boxes.append((int(p.x*downsample),int(p.y*downsample),"auto_local",owner[p.x,p.y]))
		
		print "Refine box locations"
		# Refine the box locations at full sampling
		cmpim=[]
		boxes2=[]
		for box in boxes:
			ownn=int(floor(box[3]))
#			owna=(box[3]-ownn)*360.0
#			ref=goodrefs[ownn].process("xform",{"transform":Transform({"type":"2d","alpha":owna})})
			ref=goodrefs[ownn].process("normalize.edgemean")
			ptcl=micrograph.get_clip(Region(box[0]-ref["nx"]/2,box[1]-ref["ny"]/2,ref["nx"],ref["ny"]))		# of course, nx == ny anyway
			ali=ref.align("rotate_translate",ptcl)
			ax,ay=ali["xform.align2d"].get_trans_2d()
			boxes2.append((box[0]+ax,box[1]+ay,box[2]))
#			print ownn,owna,ali["xform.align2d"]
			#cmpim.append(ref)
			#cmpim.append(ali)
			#cmpim.append(ptcl)
		#display(cmpim)
			
		print "done"
		
		return boxes2


	@staticmethod
	def ccftask(jsd,ref,downsample,microdown,ri):

		mref=ref.process("mask.soft",{"outer_radius":ref["nx"]/2-4,"width":3})
		mref.process_inplace("math.fft.resample",{"n":downsample})
		nxdown=mref["nx"]
		
		ptclmap=EMData(microdown["nx"],microdown["ny"])
		ptclmap.to_zero()
		ptclmap["ortid"]=ri
		
		# loop over the image with enough oversampling that we should be able to find all of the particles
		for x in xrange(0,microdown["nx"]-nxdown,nxdown//2):
			for y in xrange(0,microdown["ny"]-nxdown,nxdown//2):
				ptcl=microdown.get_clip(Region(x,y,nxdown,nxdown))
				# initial alignment
				ali=mref.align("rotate_translate",ptcl)
				ax,ay=ali["xform.align2d"].get_trans_2d()
				
				# extract centered particle for better frc
				ptcl=microdown.get_clip(Region(int(x+ax),int(y+ay),nxdown,nxdown))
				xf=ali["xform.align2d"]
				xf.set_trans(0,0,0)
				ali=mref.process("xform",{"transform":xf})
				frc=-ali.cmp("frc",ptcl,{"minres":200,"maxres":20,"sweight":0})		# we want larger better in this case
				
				# Write results in one pixel
				ax=int(ax+x+nxdown/2)
				ay=int(ay+y+nxdown/2)
				if frc>ptclmap[ax,ay] : ptclmap[ax,ay]=frc
		
		jsd.put(ptclmap)
		sys.stdout.write("*")
		
class boxerGauss(QtCore.QObject):
	@staticmethod
	def setup_gui(gridlay):
		return
	
	@staticmethod
	def do_autobox(boxer):
		return
	
aboxmodes = [ ("Local Search","auto_local",boxerLocal),("by Ref","auto_ref",boxerByRef), ("Gauss","auto_gauss",boxerGauss) ]
boxcolors = { "selected":(0.9,0.9,0.9), "manual":(0,0,0.7), "refgood":(0,0.8,0), "refbad":(0.8,0,0), "unknown":[.4,.4,.1], "auto_local":(.3,.1,.4), "auto_ref":(.1,.1,.4), "auto_gauss":(.4,.1,.4) }

class GUIBoxer(QtGui.QWidget):
	# Dictionary of autopickers
	# to add a new one, provide name:(Qt_setup_function,picker_execution_function)
	# Qt_setup_function(self,empty_grid_layout)
	# picker_execution_function(self,...

	
	def __init__(self,imagenames,voltage=None,apix=None,cs=None,ac=10.0,box=256,ptcl=200,threads=4):
		"""The 'new' e2boxer interface.
		"""

		QtGui.QWidget.__init__(self,None)
#		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))

		self.boxcolors=boxcolors
		self.data=None
		self.curfilename = None				# current selected file for boxing
		self.filenames=imagenames			# list of available filenames
		self.micrograph=None
		self.boxes=[]
		self.particles=[]
		self.goodrefs=[]					# "good" box references for this project
		self.goodrefchg=False				# this is used to prevent rewriting the refs many times
		self.badrefs=[]						# "bad" box references for this project. Used to reduce false-positives
		self.badrefchg=False
		self.threads=threads

		self.defaultvoltage=voltage
		self.defaultapix=apix
		self.defaultcs=cs
		self.defaultac=ac
		
		self.db = None						# open JSON file for current image

		self.wimage=EMImage2DWidget()
		self.wimage.setWindowTitle("Micrograph")

		self.wparticles=EMImageMXWidget()
		self.wparticles.setWindowTitle("Particles")
		self.wparticles.set_mouse_mode("App")
		
		self.wrefs=EMImageMXWidget()
		self.wrefs.setWindowTitle("Box Refs")
		self.wrefs.set_mouse_mode("App")
		
		self.wbadrefs=EMImageMXWidget()
		self.wbadrefs.setWindowTitle("Bad Box Refs")
		self.wbadrefs.set_mouse_mode("App")


		#self.wfft=EMImage2DWidget()
		#self.wfft.setWindowTitle("e2evalimage - 2D FFT")

		#self.wplot=EMPlot2DWidget()
		#self.wplot.setWindowTitle("e2evalimage - Plot")

		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.wimage.connect(self.wimage,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.wparticles.connect(self.wparticles,QtCore.SIGNAL("mx_image_selected"),self.ptclmousedown)
		self.wparticles.connect(self.wparticles,QtCore.SIGNAL("mx_mousedrag"),self.ptclmousedrag)
		self.wparticles.connect(self.wparticles,QtCore.SIGNAL("mx_mouseup")  ,self.ptclmouseup)
#		self.wrefs.connect(self.wparticles,QtCore.SIGNAL("mx_image_selected"),self.refmousedown)
#		self.wrefs.connect(self.wparticles,QtCore.SIGNAL("mx_mousedrag"),self.ptclmousedrag)
		self.wrefs.connect(self.wrefs,QtCore.SIGNAL("mx_mouseup")  ,self.refmouseup)
#		self.wbadrefs.connect(self.wparticles,QtCore.SIGNAL("mx_image_selected"),self.badrefmousedown)
#		self.wbadrefs.connect(self.wparticles,QtCore.SIGNAL("mx_mousedrag"),self.badrefmousedrag)
		self.wbadrefs.connect(self.wbadrefs,QtCore.SIGNAL("mx_mouseup")  ,self.badrefmouseup)

		# This object is itself a widget we need to set up
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(8)
		self.gbl.setSpacing(6)
		self.gbl.setColumnStretch(0,2)
		self.gbl.setColumnStretch(1,2)
		self.gbl.setColumnStretch(2,2)
		self.gbl.setColumnStretch(3,2)
		self.gbl.setColumnStretch(4,1)

		# Micrograph list
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		for i in imagenames:
			self.setlist.addItem(i)
		self.gbl.addWidget(self.setlist,0,0,12,2)

		self.setlist.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		self.setlist.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listKey)
		
		# Mouse Modes
		self.mmode="manual"
		self.boxmm=QtGui.QGroupBox("Mouse Mode",self)
		self.boxmm.setFlat(False)
		self.gbl.addWidget(self.boxmm,0,2,2,2)
		
		self.hbl0=QtGui.QHBoxLayout(self.boxmm)
		
		self.bmmanual=QtGui.QPushButton("Manual")
		self.bmmanual.setToolTip("Manual selection of particles. No impact on autoselection.")
		self.bmmanual.setAutoExclusive(True)
		self.bmmanual.setCheckable(True)
		self.bmmanual.setChecked(True)
		self.hbl0.addWidget(self.bmmanual)
		
		self.bmdel=QtGui.QPushButton("Delete")
		self.bmdel.setToolTip("Delete particles from any mode. Can also shift-click in other mouse modes.")
		self.bmdel.setAutoExclusive(True)
		self.bmdel.setCheckable(True)
		self.hbl0.addWidget(self.bmdel)
		
		self.bmgref=QtGui.QPushButton("Good Refs")
		self.bmgref.setToolTip("Identify some good particles. Available to all autoboxers.")
		self.bmgref.setAutoExclusive(True)
		self.bmgref.setCheckable(True)
		self.hbl0.addWidget(self.bmgref)

		self.bmbref=QtGui.QPushButton("Bad Refs")
		self.bmbref.setToolTip("Identify regions which should not be selected as particles.")
		self.bmbref.setAutoExclusive(True)
		self.bmbref.setCheckable(True)
		self.hbl0.addWidget(self.bmbref)

		QtCore.QObject.connect(self.bmmanual,QtCore.SIGNAL("clicked(bool)"),self.setMouseManual)
		QtCore.QObject.connect(self.bmdel,QtCore.SIGNAL("clicked(bool)"),self.setMouseDel)
		QtCore.QObject.connect(self.bmgref,QtCore.SIGNAL("clicked(bool)"),self.setMouseGoodRef)
		QtCore.QObject.connect(self.bmbref,QtCore.SIGNAL("clicked(bool)"),self.setMouseBadRef)

		self.bfilter=QtGui.QPushButton("Filter Disp.")
		self.bfilter.setToolTip("Filter micrograph (display only)")
		self.bfilter.setCheckable(True)
		self.gbl.addWidget(self.bfilter,0,4,1,1)
		QtCore.QObject.connect(self.bfilter,QtCore.SIGNAL("clicked(bool)"),self.filterToggle)

		self.binvert=QtGui.QPushButton("Invert")
		self.binvert.setToolTip("Invert Micrograph (also output)")
		self.binvert.setCheckable(True)
		self.binvert.setChecked(invert_on_read)		# in truly bad form, this is a global
		self.gbl.addWidget(self.binvert,1,4,1,1)
		QtCore.QObject.connect(self.binvert,QtCore.SIGNAL("clicked(bool)"),self.invertToggle)

		# Global parameters
		self.boxparm=QtGui.QGroupBox("Parameters",self)
		self.boxparm.setFlat(False)
		self.gbl.addWidget(self.boxparm,2,2,3,3)
		
		self.gbl1=QtGui.QGridLayout(self.boxparm)
		self.gbl1.setMargin(8)
		self.gbl1.setSpacing(6)
		
		self.vbbsize = ValBox(label="Box Size:",value=box)
		self.gbl1.addWidget(self.vbbsize,0,0)
		
		self.vbbpsize = ValBox(label="Ptcl Size:",value=ptcl)
		self.gbl1.addWidget(self.vbbpsize,0,1)
		
		self.vbbapix = ValBox(label="A/pix:",value=apix)
		self.gbl1.addWidget(self.vbbapix,0,2)
		
		self.vbvoltage = ValBox(label="Voltage:",value=voltage)
		self.gbl1.addWidget(self.vbvoltage,1,0)

		self.vbbac = ValBox(label="% AC:",value=ac)
		self.gbl1.addWidget(self.vbbac,1,1)
		
		self.vbcs = ValBox(label="Cs:",value=cs)
		self.gbl1.addWidget(self.vbcs,1,2)

		self.vbthreads = ValBox(label="Threads:",value=self.threads)
		self.vbthreads.setIntonly(1)
		self.gbl1.addWidget(self.vbthreads,1,3)

		# Reference tools
		self.reftools=QtGui.QGroupBox("Box Refs",self)
		self.boxparm.setFlat(False)
		self.gbl.addWidget(self.reftools,5,2,2,2)
		
		self.hbl1=QtGui.QHBoxLayout(self.reftools)
		
		self.rtload3d=QtGui.QPushButton("From 3D")
		self.rtload3d.setToolTip("Load box refs from 3-D volume")
		self.hbl1.addWidget(self.rtload3d)
		
		self.rtload2d=QtGui.QPushButton("From 2D")
		self.rtload2d.setToolTip("Load box refs from 2-D stack")
		self.hbl1.addWidget(self.rtload2d)

		self.rtclear=QtGui.QPushButton("Clear")
		self.rtclear.setToolTip("Clear all current good and bad refs")
		self.hbl1.addWidget(self.rtclear)

		QtCore.QObject.connect(self.rtload3d,QtCore.SIGNAL("clicked(bool)"),self.reftoolLoad3D)
		QtCore.QObject.connect(self.rtload2d,QtCore.SIGNAL("clicked(bool)"),self.reftoolLoad2D)
		QtCore.QObject.connect(self.rtclear,QtCore.SIGNAL("clicked(bool)"),self.reftoolClear)
		
		
		# Autoboxing Tabs
		self.autolbl = QtGui.QLabel("Autoboxing Methods:")
		self.gbl.addWidget(self.autolbl,7,2)
		self.autotab = QtGui.QTabWidget()
		self.gbl.addWidget(self.autotab,8,2,5,3)
		
		# Individual tabs from Dictionary
		self.abwid=[]
		for name,bname,cls in aboxmodes:
			w=QtGui.QWidget()
			gl=QtGui.QGridLayout(w)
			self.abwid.append((w,gl))
			cls.setup_gui(gl)
			self.autotab.addTab(w,name)
			
		self.bbclear=QtGui.QPushButton("Clear Boxes")
		self.bbclear.setToolTip("Clear all boxes in current micrograph")
		self.gbl.addWidget(self.bbclear,13,2)
		QtCore.QObject.connect(self.bbclear,QtCore.SIGNAL("clicked(bool)"),self.boxClear)

		self.bautoboxa = QtGui.QPushButton("Autobox All")
		self.gbl.addWidget(self.bautoboxa,13,3)
		QtCore.QObject.connect(self.bautoboxa,QtCore.SIGNAL("clicked(bool)"),self.doAutoBoxAll)
		
		self.bautobox = QtGui.QPushButton("Autobox")
		self.gbl.addWidget(self.bautobox,13,4)
		QtCore.QObject.connect(self.bautobox,QtCore.SIGNAL("clicked(bool)"),self.doAutoBox)

		self.setWindowTitle("e2boxer21 - Control Panel")

		self.wimage.show()
		self.wparticles.set_data(self.particles)
		self.wparticles.show()
		
		if os.path.exists("info/boxrefs.hdf"):
			self.goodrefs=EMData.read_images("info/boxrefs.hdf")
		self.wrefs.set_data(self.goodrefs)
		if len(self.goodrefs)>0 : self.wrefs.show()
		
		if os.path.exists("info/boxrefsbad.hdf"):
			self.badrefs=EMData.read_images("info/boxrefsbad.hdf")
		self.wbadrefs.set_data(self.badrefs)
		if len(self.badrefs)>0 : self.wbadrefs.show()
		
#		self.wfft.show()
#		self.wplot.show()
		E2loadappwin("e2boxer21","main",self)
		E2loadappwin("e2boxer21","image",self.wimage.qt_parent)
		E2loadappwin("e2boxer21","particles",self.wparticles.qt_parent)
		E2loadappwin("e2boxer21","refs",self.wrefs.qt_parent)
		E2loadappwin("e2boxer21","badrefs",self.wbadrefs.qt_parent)

		self.newSet(0)

	def setMouseManual(self,x):
		self.mmode="manual"
	
	def setMouseDel(self,x):
		self.mmode="del"
	
	def setMouseGoodRef(self,x):
		self.mmode="refgood"
	
	def setMouseBadRef(self,x):
		self.mmode="refbad"
	
	def reftoolLoad3D(self,x):
		fsp=str(QtGui.QFileDialog.getOpenFileName(self, "Select 3-D Volume"))
		if fsp==None or len(fsp)<4 : return

		symname=str(QtGui.QInputDialog.getText(None,"Symmetry","Please specify the symmetry of the map, or c1 for none")[0])
#		print symname
		try:
			sym = Symmetries.get(symname)
		except:
			print "Error: Unknown symmetry"
			return
		orts=sym.gen_orientations("eman",{"delta":15,"inc_mirror":1})
		prog=QtGui.QProgressDialog("Making Projections","Abort",0,len(orts))
		prog.setWindowModality(Qt.WindowModal)
		prog.setValue(0)
		
		vol=EMData(fsp,0)
		vol.process_inplace("normalize.circlemean")
		apix3=vol["apix_x"]
		apix1=self.vbbapix.getValue()
		if apix3<0.1 : 
			print "WARNING: A/pix on the 3-D volume appears too small. This method only works if the volume has a valid A/pix value in its header. I am adjusting it to 1.0, but this is almost certainly wrong."
			apix3=1.0
		
		xsize3d=vol["nx"]
		xsize=self.vbbsize.getValue()

		if ( xsize3d != xsize or fabs(fabs(apix1/apix3)-1.0)>.001 ) :
			print "WARNING: the boxsize and/or sampling (%d @ %1.4f A/pix) do not match (%d @ %1.4f A/pix). I will attempt to adjust the volume appropriately."%(xsize,apix1,xsize3d,apix3)
			try:
				scale=apix3/apix1
				print "Reference is {box3} x {box3} x {box3} at {apix3:1.2f} A/pix, particles are {box2} x {box2} at {apix2:1.2f} A/pix. Scaling by {scale:1.3f}".format(box3=xsize3d,box2=xsize,apix3=apix3,apix2=apix1,scale=scale)
			except:
				print "A/pix unknown, assuming scale same as relative box size"
				scale=float(xsize)/xsize3d
				
			vol.process_inplace("xform.scale",{"clip":xsize,"scale":scale})
		

		for i,o in enumerate(orts):
			prog.setValue(i)
			if prog.wasCanceled() : break
			proj=vol.project("standard",o)
			proj.process_inplace("normalize.edgemean")
			self.goodrefs.append(proj)
			
		prog.setValue(len(orts))
		
		self.wrefs.set_data(self.goodrefs)
		self.wrefs.show()
		self.goodrefchg=True
		
	def reftoolLoad2D(self,x):
		fsp=str(QtGui.QFileDialog.getOpenFileName(self, "Select 2-D Stack"))
		if fsp==None or len(fsp)<4 : return
		
		refs=EMData.read_images(fsp)

		apix2=refs[0]["apix_x"]
		apix1=self.vbbapix.getValue()
		if apix2<0.1 : 
			print "WARNING: A/pix on the 2-D images appear too small. This method only works if the images have a valid A/pix value in their header. I am adjusting it to 1.0, but this is almost certainly wrong."
			apix2=1.0
		
		xsize2=refs[0]["nx"]
		xsize=self.vbbsize.getValue()
		scale=1.0

		if ( xsize2 != xsize or fabs(fabs(apix1/apix2)-1.0)>.001 ) :
			print "WARNING: the boxsize and/or sampling (%d @ %1.4f A/pix) do not match (%d @ %1.4f A/pix). I will attempt to adjust the volume appropriately."%(xsize,apix1,xsize2,apix2)
			try:
				scale=apix2/apix1
				print "Reference is {box3} x {box3} x {box3} at {apix3:1.2f} A/pix, particles are {box2} x {box2} at {apix2:1.2f} A/pix. Scaling by {scale:1.3f}".format(box3=xsize2,box2=xsize,apix3=apix2,apix2=apix1,scale=scale)
			except:
				print "A/pix unknown, assuming scale same as relative box size"
				scale=float(xsize)/xsize2
			
		for r in refs:
			r.process_inplace("normalize.circlemean")
			if scale!=1.0 or xsize!=xsize2 : r.process_inplace("xform.scale",{"clip":xsize,"scale":scale})

		self.goodrefs.extend(refs)
		self.goodrefchg=True
		self.wrefs.set_data(self.goodrefs)
		self.wrefs.show()
	
		
	def reftoolClear(self,x):
		r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: this will remove all good and bad box references. Are you sure?",QtGui.QMessageBox.Yes|QtGui.QMessageBox.Cancel)
		if r==QtGui.QMessageBox.Cancel : return

		self.goodrefs=[]
		self.goodrefchg=True
		self.wrefs.set_data(self.goodrefs)
		
		self.badrefs=[]
		self.badrefchg=True
		self.wbadrefs.set_data(self.badrefs)

	def boxClear(self,x):
		r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: this will erase all box locations in the current micrograph. Are you sure?",QtGui.QMessageBox.Yes|QtGui.QMessageBox.Cancel)
		if r==QtGui.QMessageBox.Cancel : return

		self.boxes=[]
		self.__updateBoxes()
		
	def filterToggle(self,x):
		self.__show_image()

	def invertToggle(self,x):
		global invert_on_read
		invert_on_read=self.binvert.isChecked()
		self.newSet()
		
	def imgmousedown(self,event,m) :
#		m=self.wimage.scr_to_img((event.x(),event.y()))
		self.curbox=-1
		self.downbut=event.buttons()
		if not event.buttons()&Qt.LeftButton : return
		self.lastloc=m
		boxsize2=self.vbbsize.getValue()//2
		ptclsize=self.vbbpsize.getValue()
		if boxsize2<4 : return

		if self.mmode=="del" : return 		# del works with drag and up
		elif self.mmode in ("refgood","refbad") :
			self.tmpbox=m
			try: color=self.boxcolors[self.mmode]
			except: color=self.boxcolors["unknown"]
			self.wimage.add_shape("tmpbox",EMShape(("rect",color[0],color[1],color[2],self.tmpbox[0]-boxsize2,self.tmpbox[1]-boxsize2,self.tmpbox[0]+boxsize2,self.tmpbox[1]+boxsize2,2)))
			self.wimage.update()
		elif self.mmode == "manual" :
			self.curbox=0
			# check to see if click was inside an existing box, in which case we move it
			for i,b in enumerate(self.boxes):
				if abs(m[0]-b[0])<boxsize2 and abs(m[1]-b[1])<boxsize2 : 
					self.curbox=i
					break
			else :
				# Create a new box
				self.curbox=len(self.boxes)
				self.boxes.append((m[0],m[1],self.mmode))
				self.__addBox(self.curbox,self.boxes[-1])
			
		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])

	def imgmousedrag(self,event,m) :
#		m=self.wimage.scr_to_img((event.x(),event.y()))
		if not event.buttons()&Qt.LeftButton : return
		boxsize2=self.vbbsize.getValue()//2
		ptclsize=self.vbbpsize.getValue()
		if boxsize2<4 : return
		
		if self.mmode=="del" or event.modifiers()&Qt.ShiftModifier:
			# filter out all matching boxes, then update all at once
			n=len(self.boxes)
			self.boxes=[b for b in self.boxes if abs(b[0]-m[0])>boxsize2 or abs(b[1]-m[1])>boxsize2]
			if len(self.boxes)<n : self.__updateBoxes()
		elif self.mmode in ("refgood","refbad"):
			if m==self.lastloc : return
			b=self.tmpbox
			self.tmpbox=(b[0]+m[0]-self.lastloc[0],b[1]+m[1]-self.lastloc[1],self.mmode)
			self.lastloc=m
			try: color=self.boxcolors[self.mmode]
			except: color=self.boxcolors["unknown"]
			self.wimage.add_shape("tmpbox",EMShape(("rect",color[0],color[1],color[2],self.tmpbox[0]-boxsize2,self.tmpbox[1]-boxsize2,self.tmpbox[0]+boxsize2,self.tmpbox[1]+boxsize2,2)))
			self.wimage.update()
		elif self.mmode=="manual":
			if m==self.lastloc : return
			b=self.boxes[self.curbox]
			self.boxes[self.curbox]=(b[0]+m[0]-self.lastloc[0],b[1]+m[1]-self.lastloc[1],b[2])
			self.__addBox(self.curbox,self.boxes[self.curbox])
			self.lastloc=m

			
	def imgmouseup(self,event,m) :
#		m=self.wimage.scr_to_img((event.x(),event.y()))
		if not self.downbut&Qt.LeftButton : return
		boxsize2=self.vbbsize.getValue()//2
		ptclsize=self.vbbpsize.getValue()
		if boxsize2<4 : return
		
		if self.mmode=="del" or event.modifiers()&Qt.ShiftModifier:
			# filter out all matching boxes, then update all at once
			n=len(self.boxes)
			self.boxes=[b for b in self.boxes if abs(b[0]-m[0])>boxsize2 or abs(b[1]-m[1])>boxsize2]
			if len(self.boxes)<n : self.__updateBoxes()
		elif self.mmode in ("refgood","refbad"):
			b=self.tmpbox
			self.tmpbox=(b[0]+m[0]-self.lastloc[0],b[1]+m[1]-self.lastloc[1],self.mmode)
			self.lastloc=m
			self.wimage.del_shape("tmpbox")
			boxim=self.micrograph.get_clip(Region(self.tmpbox[0]-boxsize2,self.tmpbox[1]-boxsize2,boxsize2*2,boxsize2*2))
			boxim["ptcl_source_coord"]=(self.tmpbox[0],self.tmpbox[1])
			if self.mmode == "refgood" : 
				self.goodrefs.append(boxim)
				self.goodrefchg=True
				self.wrefs.set_data(self.goodrefs)
				self.wrefs.show()
			else :
				self.badrefs.append(boxim)
				self.badrefchg=True
				self.wbadrefs.set_data(self.badrefs)
				self.wbadrefs.show()
		elif self.mmode=="manual":
			if m==self.lastloc : return
			b=self.boxes[self.curbox]
			self.boxes[self.curbox]=(b[0]+m[0]-self.lastloc[0],b[1]+m[1]-self.lastloc[1],b[2])
			self.__addBox(self.curbox,self.boxes[self.curbox])


	def ptclmousedown(self,event,m) :
		self.downbut=event.buttons()
		if self.mmode == "manual" :
			self.curbox=m[0]
			self.lastloc=m[1:3]
			self.wimage.scroll_to(self.boxes[m[0]][0],self.boxes[m[0]][1])
				
	def ptclmousedrag(self,event,x) :
		m=self.wparticles.scr_to_img((event.x(),event.y()))
		if self.mmode=="manual":
			if m[1:3]==self.lastloc : return
			b=self.boxes[self.curbox]
			self.boxes[self.curbox]=(b[0]-m[1]+self.lastloc[0],b[1]-m[2]+self.lastloc[1],b[2])
			self.__addBox(self.curbox,self.boxes[self.curbox])
			self.lastloc=m[1:3]

	def ptclmouseup(self,event,m) :
		if self.mmode=="del" or event.modifiers()&Qt.ShiftModifier:
			self.boxes.pop(m[0])
			self.__updateBoxes()
		elif self.mmode=="manual":
			if m[1:3]==self.lastloc : return
			b=self.boxes[self.curbox]
			self.boxes[self.curbox]=(b[0]-m[1]+self.lastloc[0],b[1]-m[2]+self.lastloc[1],b[2])
			self.__addBox(self.curbox,self.boxes[self.curbox])
			self.lastloc=m[1:3]
		return

	def refmouseup(self,event,m) :
		print "refup"
		if self.mmode=="del" or event.modifiers()&Qt.ShiftModifier:
			self.goodrefs.pop(m[0])
			self.goodrefchg=True
			self.wrefs.set_data(self.goodrefs)
		return

	def badrefmouseup(self,event,m) :
		print "badrefup"
		if self.mmode=="del" or event.modifiers()&Qt.ShiftModifier:
			self.badrefs.pop(m[0])
			self.badrefchg=True
			self.wbadrefs.set_data(self.badrefs)
		return

	def newSet(self,val=None):
		"called when a new data set is selected from the list"

		first=True
		if val==None : newfilename=self.curfilename
		else : newfilename=str(self.setlist.item(val).text()).split()[1]
#		if newfilename==self.curfilename : return

		# Write the current image parameters to the database
		if self.curfilename!=None and self.boxes!=None :
			self.save_boxes()
			first=False


		self.micrograph=load_micrograph(newfilename)
		self.__show_image()
		if first : E2loadappwin("e2boxer21","image",self.wimage.qt_parent)
		self.curfilename=newfilename
		self.restore_boxes()

	def __show_image(self):
		if self.bfilter.isChecked() :
			nx=self.micrograph["nx"]
			ny=self.micrograph["ny"]
			gs=good_size(max(self.micrograph["nx"],self.micrograph["ny"]))
			fm=self.micrograph.get_clip(Region(0,0,gs,gs)).do_fft()
			fm.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
			fm.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
			fm=fm.do_ift()
			self.wimage.set_data(fm.get_clip(Region(0,0,nx,ny)))
		else: self.wimage.set_data(self.micrograph)

	def save_boxes(self):
		js=js_open_dict(info_name(self.curfilename))
		js["boxes"]=self.boxes
		
		if self.goodrefchg :
			try: os.unlink("info/boxrefs.hdf")
			except: pass
			for im in self.goodrefs: im.write_image("info/boxrefs.hdf",-1)
			self.goodrefchg=False
			
		if self.badrefchg :
			try: os.unlink("info/boxrefsbad.hdf")
			except: pass
			for im in self.badrefs: im.write_image("info/boxrefsbad.hdf",-1)
			self.badrefchg=False
			
	def restore_boxes(self):
		# first we restore the list of box locations
		js=js_open_dict(info_name(self.curfilename))
		try: self.boxes=js["boxes"]
		except: self.boxes=[]
		self.__updateBoxes()
		
		#boxsize=self.vbbsize.getValue()
		#ptclsize=self.vbbpsize.getValue()
		#micro=self.wimage.get_data()
		
		#self.wimage.del_shapes()
		#self.particles=[]
		#if len(self.boxes)==0 : 
			#self.wparticles.set_data(self.particles)
			#return
		
		## Then we extract the actual boxed out particles
		#for i,box in enumerate(self.boxes):
			#self.__addBox(i,box)
		
		## finally redisplay as appropriate
		#self.wimage.update()
		#self.wparticles.set_data(self.particles)
		#if len(self.particles)>0 : self.wparticles.show()
				

	def __updateBoxes(self):
		"""Updates the display to match current self.boxes"""
		
		self.wimage.del_shapes()
		
		self.particles=[]
		for i,b in enumerate(self.boxes):
			self.__addBox(i,b,True)
			
		# finally redisplay as appropriate
		self.wimage.update()
		self.wparticles.set_data(self.particles)
		if len(self.particles)>0 : self.wparticles.show()
		self.__updateCount()

	def __updateCount(self):
		row=self.setlist.currentRow()
		itm=self.setlist.item(row)
		try: lbl=str(itm.text()).split()
		except: return
		if int(lbl[0])!= len(self.boxes):
			lbl="{}\t{}".format(len(self.boxes),lbl[1])
			itm.setText(lbl)

	def __addBox(self,i,box,quiet=False):
		"""takes the number of the box in self.boxes and the (x,y,mode) tuple and displays it"""
		# Display the actual box
		boxsize=self.vbbsize.getValue()
		ptclsize=self.vbbpsize.getValue()
		try: color=self.boxcolors[box[2]]
		except: color=self.boxcolors["unknown"]
		self.wimage.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],box[0]-boxsize//2,box[1]-boxsize//2,box[0]+boxsize//2,box[1]+boxsize//2,2)))
		self.wimage.add_shape("cir{}".format(i),EMShape(("circle",color[0],color[1],color[2],box[0],box[1],ptclsize/2,1.5)))

		boxim=self.micrograph.get_clip(Region(box[0]-boxsize//2,box[1]-boxsize//2,boxsize,boxsize))
		boxim["ptcl_source_coord"]=(box[0],box[1])
		if len(self.particles)<=i : self.particles.append(boxim)
		else: self.particles[i]=boxim
		
		if not quiet:
			# finally redisplay as appropriate
			self.wimage.update()
			self.wparticles.set_data(self.particles)
			self.__updateCount()

	def listKey(self,event):
		pass

		#if event.key()>=Qt.Key_0 and event.key()<=Qt.Key_9 :
			#q=int(event.key())-Qt.Key_0
			#self.squality.setValue(q)
		#elif event.key() == Qt.Key_Left:
			#self.sdefocus.setValue(self.sdefocus.getValue()-0.03)
		#elif event.key() == Qt.Key_Right:
			#self.sdefocus.setValue(self.sdefocus.getValue()+0.03)
		#elif event.key()==Qt.Key_I :
			#self.doImport()
		#elif event.key()==Qt.Key_U :
			#self.unImport()

		
	def doAutoBox(self,b):
		"""Autobox button pressed, find the right algorithm and call it"""
		
		name,bname,cls=aboxmodes[self.autotab.currentIndex()]
		
		print name," called"

		boxes=cls.do_autobox(self.micrograph,self.goodrefs,self.badrefs,self.vbbapix.getValue(),self.vbthreads.getValue(),{})
		
		# if we got nothing, we just leave the current results alone
		if len(boxes)==0 : return
	
		bname=boxes[0][2]

		# Filter out all existing boxes for this picking mode
		self.boxes=[i for i in self.boxes if i[2]!=bname]

		self.boxes.extend(boxes)
		self.__updateBoxes()

	def doAutoBoxAll(self,b):
		"""Autobox button pressed, find the right algorithm and call it"""
		
		name,bname,cls=aboxmodes[self.autotab.currentIndex()]
		
		prog=QtGui.QProgressDialog("Autoboxing","Abort",0,len(self.filenames))
		prog.setWindowModality(Qt.WindowModal)
		prog.setValue(0)
		prog.show()

		
		for i,fspl in enumerate(self.filenames):
			
			fsp=fspl.split()[1]
			prog.setValue(i)
			if prog.wasCanceled() : 
				print "Autoboxing Aborted!"
				break
			
			micrograph=load_micrograph(fsp)

			newboxes=cls.do_autobox(micrograph,self.goodrefs,self.badrefs,self.vbbapix.getValue(),self.vbthreads.getValue(),{},prog)
			print "{}) {} boxes -> {}".format(i,len(newboxes),fsp)
			
			# if we got nothing, we just leave the current results alone
			if len(newboxes)==0 : continue
		
			# read the existing box list and update
			db=js_open_dict(info_name(fsp))
			try: 
				boxes=db["boxes"]
				# Filter out all existing boxes for this picking mode
				bname=newboxes[0][2]
				boxes=[b for b in boxes if b[2]!=bname]
			except:
				boxes=[]
				
			boxes.extend(newboxes)
			
			db["boxes"]=boxes
			db.close()
			self.setlist.setCurrentRow(i)
#			self.__updateBoxes()
			
		else: prog.setValue(len(self.filename))
		
		self.restore_boxes()

		
	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
		self.save_boxes()
		E2saveappwin("e2boxer21","main",self)
		E2saveappwin("e2boxer21","image",self.wimage.qt_parent)
		E2saveappwin("e2boxer21","particles",self.wparticles.qt_parent)
		E2saveappwin("e2boxer21","refs",self.wrefs.qt_parent)
		E2saveappwin("e2boxer21","badrefs",self.wbadrefs.qt_parent)

		#self.writeCurParm()
		event.accept()
		QtGui.qApp.exit(0)
		#app=QtGui.qApp			if b[2] not in ("refgood","refbad"):
		#if self.wimage != None:
			#app.close_specific(self.wimage)
			#self.wimage = None
		#if self.wfft != None:
			#app.close_specific(self.wfft)
		#if self.wplot != None:
			#app.close_specific(self.wplot)
		#app.close_specific(self)
#		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop


	#def update_plot(self):
##		if self.wplot == None: return # it's closed/not visible

		#if self.incalc : return		# no plot updates during a recomputation

		#parms=self.parms[self.curset]
		#apix=self.sapix.getValue()
		#ds=1.0/(apix*parms[0]*parms[5])
		#ctf=parms[1]
		#bg1d=array(ctf.background)
		#r=len(ctf.background)
		#s=arange(0,ds*r,ds)

		## This updates the FFT CTF-zero circles
		#if self.f2danmode==0 :
			#fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)
			#shp={}
			#nz=0
			#for i in range(1,len(fit)):
				#if fit[i-1]*fit[i]<=0.0:
					#nz+=1
					#shp["z%d"%i]=EMShape(("circle",0.0,0.0,1.0/nz,r,r,i,1.0))

			#self.wfft.del_shapes()
			#self.wfft.add_shapes(shp)
			#self.wfft.updateGL()
		## Single measurement circle mode
		#elif self.f2danmode==1 :
			#self.wfft.del_shapes()
			#if self.ringrad==0: self.ringrad=1.0
			#self.wfft.add_shape("ring",EMShape(("circle",0.2,1.0,0.2,r,r,self.ringrad,1.0)))
			#self.wfft.add_shape("ringlbl",EMShape(("scrlabel",0.2,1.0,0.2,10,10,"r=%d pix -> 1/%1.2f 1/A (%1.4f)"%(self.ringrad,1.0/(self.ringrad*ds),self.ringrad*ds),24.0,2.0)))
			#self.wfft.updateGL()
		## 2-D Crystal mode
		#elif self.f2danmode==2 :
			#shp={}
			#for a in range(-5,6):
				#for b in range(-5,6):
					#shp["m%d%d"%(a,b)]=EMShape(("circle",1.0,0.0,0.0,a*self.xpos1[0]+b*self.xpos2[0]+self.fft["nx"]/2-1,a*self.xpos1[1]+b*self.xpos2[1]+self.fft["ny"]/2,3,1.0))

			#self.wfft.del_shapes()
			#self.wfft.add_shapes(shp)
			#self.wfft.add_shape("xtllbl",EMShape(("scrlabel",1.0,0.3,0.3,10,10,"Unit Cell: %1.2f,%1.2f"%(1.0/(hypot(*self.xpos1)*ds),1.0/(hypot(*self.xpos2)*ds)),60.0,2.0)))
##			except: pass
			#self.wfft.updateGL()
		#else:
			#self.wfft.del_shapes()
			#self.wfft.updateGL()


		## Now update the plots for the correct plot mode
		#if self.plotmode==0:
			#try: bgsub=self.fft1d-bg1d
			#except:
				#print "Error computing bgsub on this image"
				#return
			#self.wplot.set_data((s,bgsub),"fg-bg",quiet=True,color=0)

			#fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			#fit=fit*fit			# squared

			## auto-amplitude for b-factor adjustment
			#rto,nrto=0,0
			#for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				#if bgsub[i]>0 :
					#rto+=fit[i]
					#nrto+=fabs(bgsub[i])
			#if nrto==0 : rto=1.0
			#else : rto/=nrto
			#fit=[fit[i]/rto for i in range(len(s))]

##			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			#self.wplot.set_data((s,fit),"fit",color=1)
			#self.wplot.setAxisParms("s (1/"+ "$\AA$" +")","Intensity (a.u)")
		#elif self.plotmode==1:
			#self.wplot.set_data((s[1:],self.fft1d[1:]),"fg",quiet=True,color=1)
			#self.wplot.set_data((s[1:],bg1d[1:]),"bg",color=0)
			#self.wplot.setAxisParms("s (1/"+ "$\AA$" +")","Intensity (a.u)")
		#elif self.plotmode==2:
			#if self.fft1dang==None: self.recalc_real()
			#bgsub=self.fft1d-bg1d
			#bgsuba=[array(self.fft1dang[i])-bg1d for i in xrange(4)]
					## Write the current image parameters to the database

##			for i in xrange(4): bgsuba[i][0]=0
			#self.wplot.set_data((s,bgsub),"fg",quiet=True,color=0)
			#self.wplot.set_data((s[3:],bgsuba[0][3:]),"fg 0-45",quiet=True,color=2)
			#self.wplot.set_data((s[3:],bgsuba[1][3:]),"fg 45-90",quiet=True,color=3)
			#self.wplot.set_data((s[3:],bgsuba[2][3:]),"fg 90-135",quiet=True,color=4)
			#self.wplot.set_data((s[3:],bgsuba[3][3:]),"fg 135-180",quiet=True,color=6)

			#fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			#fit=fit*fit			# squared

			## auto-amplitude for b-factor adjustment
			#rto,nrto=0,0
			#for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				#if bgsub[i]>0 :
					#rto+=fit[i]
					#nrto+=fabs(bgsub[i])
			#if nrto==0 : rto=1.0
			#else : rto/=nrto
			#fit/=rto

			#self.wplot.set_data((s,fit),"fit",color=1)
			#self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")

		#elif self.plotmode==3:
			#if self.fft1dang==None: self.recalc_real()
			##bgsub=self.fft1d-bg1d
			##bgsuba=[array(self.fft1dang[i])-bg1d for i in xrange(4)]
			#fg=self.fft1d
			#fga=[array(self.fft1dang[i]) for i in xrange(4)]

			#for i in xrange(4): fga[i][0]=0
			#self.wplot.set_data((s,fg),"fg",quiet=True,color=0)
			#self.wplot.set_data((s,fga[0]),"fg 0-45",quiet=True,color=2)
			#self.wplot.set_data((s,fga[1]),"fg 45-90",quiet=True,color=3)
			#self.wplot.set_data((s,fga[2]),"fg 90-135",quiet=True,color=4)
			#self.wplot.set_data((s,fga[3]),"fg 135-180",quiet=True,color=6)

			##fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			##fit=fit*fit			# squared

			### auto-amplitude for b-factor adjustment
			##rto,nrto=0,0
			##for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				##if bgsub[i]>0 :
					##rto+=fit[i]
					##nrto+=fabs(bgsub[i])
			##if nrto==0 : rto=1.0
			##else : rto/=nrto
			##fit/=rto

			##self.wplot.set_data((s,fit),"fit",color=1)
			##self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")
		#if self.plotmode==4:
			#if min(bg1d)<=0.0 : bg1d+=min(bg1d)+max(bg1d)/10000.0
			#ssnr=(self.fft1d-bg1d)/bg1d
			#self.wplot.set_data((s,ssnr),"SSNR",quiet=True,color=0)

			##fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			##fit=fit*fit			# squared

			### auto-amplitude for b-factor adjustment
			##rto,nrto=0,0
			##for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				##if bgsub[i]>0 :
					##rto+=fit[i]
					##nrto+=fabs(bgsub[i])
			##if nrto==0 : rto=1.0
			##else : rto/=nrto
			##fit=[fit[i]/rto for i in range(len(s))]

###			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			##self.wplot.set_data((s,fit),"fit",color=1)
			#self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Est. SSNR")


	#def timeOut(self):
		#if self.busy : return

		## Redisplay before spawning thread for more interactive display
		#if self.needredisp :
			#try: self.redisplay()
			#except: pass

		## Spawn a thread to reprocess the data
		#if self.needupdate and self.procthread==None:
			#self.procthread=threading.Thread(target=self.recalc_real)
			#self.procthread.start()

		#if self.errors:
			#QtGui.QMessageBox.warning(None,"Error","The following processors encountered errors during processing of 1 or more images:"+"\n".join(self.errors))
			#self.errors=None

	#def doRefit(self):
		#parms=self.parms[self.curset]
		#apix=self.sapix.getValue()
		#ds=1.0/(apix*parms[0]*parms[5])
		
		#try:
			#parms[1]=e2ctf.ctf_fit(self.fft1d,parms[1].background,parms[1].background,self.fft,self.fftbg,parms[1].voltage,parms[1].cs,parms[1].ampcont,apix,bgadj=False,autohp=True,verbose=1)
		#except:
			#print "CTF Autofit Failed"
			#traceback.print_exc()
			#parms[1].defocus=1.0

		#self.sdefocus.setValue(parms[1].defocus,True)
		#self.sbfactor.setValue(parms[1].bfactor,True)
		#self.sampcont.setValue(parms[1].ampcont,True)
		
		#self.update_plot()


	#def unImport(self,val=None):
		#print "unimport ",base_name(self.setlist.item(self.curset).text())
		#item=base_name(self.setlist.item(self.curset).text())
		#try: os.unlink("micrographs/%s.hdf"%item)
		#except: print "Couldn't delete micrographs/%s.hdf"%item

	#def doImport(self,val=None):
		#"""Imports the currently selected image into a project"""
		#print "import ",base_name(self.setlist.item(self.curset).text())

		## This is just the (presumably) unique portion of the filename
		#item=base_name(self.setlist.item(self.curset).text())

		## create directory if necessary
		#if not os.access("micrographs",os.R_OK) :
			#try : os.mkdir("micrographs")
			#except:
				#QtGui.QMessageBox.warning(self,"Error !","Cannot create micrographs directory")
				#return

		##db=db_open_dict("bdb:micrographs#%s"%item)
		#self.data["ctf"]=self.parms[self.curset][1]
		
		#if self.cinvert.getValue()!=0 : self.data.mult(-1)
		#if self.cxray.getValue() : self.data.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4,"tomean":1})
		#self.data.write_image("micrographs/%s.hdf"%item)
		#self.writeCurParm()

##		js_open_dict(info_name(item))["ctf"]=[self.parms[val][1],None,None,None,None]
## 		db_parms=db_open_dict("bdb:e2ctf.parms")
## 		db_parms[item]=[self.parms[val][1].to_string(),self.fft1d,self.parms[val][1].background,self.parms[val][4]]

	#def writeCurParm(self):
		#"Called to store the current parameters for this image to the frameparms database"
		#parms=self.parms[self.curset]
		#js=js_open_dict(info_name(self.curfilename))
		#js.setval("ctf_frame",parms,True)
		#js.setval("quality",parms[4])
## 		db_fparms=db_open_dict("bdb:e2ctf.frameparms")
## 		curtag=item_name(str(self.setlist.item(self.curset).text()))
## 		db_fparms[curtag]=self.parms[self.curset]


	#def recalc(self):
		#self.needupdate=True

	#def recalc_real(self):
		#"Called to recompute the power spectra, also updates plot"

		#self.needupdate=False

		#if self.data==None :
			#self.procthread=None
			#return

		#self.incalc=True	# to avoid incorrect plot updates

		## To simplify expressions
		#parms=self.parms[self.curset]
		#apix=self.sapix.getValue()
		#if len(parms)==5 : parms.append(1)		# for old projects where there was no oversampling specification
		#else: parms[5]=max(1,int(parms[5]))
		#ds=1.0/(apix*parms[0]*parms[5])

		## Mode where user drags the box around the parent image
		#if self.calcmode==0:

			## extract the data and do an fft
			#clip=self.data.get_clip(Region(parms[2][0],parms[2][1],parms[0],parms[0]))
			#clip.process_inplace("normalize.edgemean")

			#if parms[5]>1 :
				#clip=clip.get_clip(Region(0,0,parms[0]*parms[5],parms[0]*parms[5]))		# since we aren't using phases, doesn't matter if we center it or not
			#self.fft=clip.do_fft()
##			self.fft.mult(1.0/parms[0]**2)
			#self.fft.mult(1.0/parms[0])

		## mode where user selects/deselcts tiled image set
		#elif self.calcmode==1:
			## update the box display on the image
			#nx=self.data["nx"]/parms[0]-1
			#self.fft=None
			#nbx=0
			#for x in range(nx):
				#for y in range(self.data["ny"]/parms[0]-1):
					## User deselected this one
					#if int(x+y*nx) in parms[3] : continue

					## read the data and make the FFT
					#clip=self.data.get_clip(Region(x*parms[0]+parms[0]/2,y*parms[0]+parms[0]/2,parms[0],parms[0]))
					#clip.process_inplace("normalize.edgemean")
					#if parms[5]>1 :
						#clip=clip.get_clip(Region(0,0,parms[0]*parms[5],parms[0]*parms[5]))		# since we aren't using phases, doesn't matter if we center it or not
					#fft=clip.do_fft()
##					fft.mult(parms[0])
					#fft.ri2inten()
					#if self.fft==None: self.fft=fft
					#else: self.fft+=fft
					#nbx+=1

			#self.fft.mult(1.0/(nbx*parms[0]**2))
			#self.fft.process_inplace("math.sqrt")
			#self.fft["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
##			self.fft.mult(1.0/(nbx*parms[0]**2))

		#self.fftbg=self.fft.process("math.nonconvex")
		#self.fft1d=self.fft.calc_radial_dist(self.fft.get_ysize()/2,0.0,1.0,1)	# note that this handles the ri2inten averages properly
		#if self.plotmode==2 or self.plotmode==3:
			#self.fft1dang=array(self.fft.calc_radial_dist(self.fft.get_ysize()/2,0.0,1.0,4,self.sang45.getValue()*.017453292,1))	# This form generates 4 sequential power spectra representing angular ranges
			#self.fft1dang=self.fft1dang.reshape((4,self.fft.get_ysize()/2))
		#else:
			#self.fft1dang=None

		## Compute 1-D curve and background
		#bg_1d=e2ctf.low_bg_curve(self.fft1d,ds)
		#parms[1].background=bg_1d
		#parms[1].dsbg=ds

		#self.fft1d=array(self.fft1d)

		#self.needredisp=True
		#self.incalc=False
		#time.sleep(.2)			# help make sure update has a chance
		#self.procthread=None
##		dbquality = self.db[os.path.basename(self.curfilename)]
##		print dbquality
## 		item=item_name(self.curfilename)
## 		db=db_open_dict("bdb:e2ctf.parms")
## 		if db[item]==None or db[item[3]]==None : db[item]=[parms[1].to_string(),self.fft1d,parms[1].background,5]
## #		print item,db[item][3]
## 		try: self.squality.setValue(int(db[item[3]]))
## 		except: self.squality.setValue(5)

	#def redisplay(self):

		#if self.incalc: return

		#self.needredisp=False
		#self.busy=True
		#parms=self.parms[self.curset]
		#apix=self.sapix.getValue()
		#ds=1.0/(apix*parms[0]*parms[5])


		## Fitting not done yet. Need to make 2D background somehow
		#if parms[1].defocus==0:
			#self.doRefit()

		#self.wimage.show()
		#self.wfft.show()
		#self.wplot.show()

		#self.update_plot()

		## To simplify expressions

		#if self.calcmode==0:
			## update the box display on the image
			#self.wimage.del_shapes()
			#self.wimage.add_shape("box",EMShape(("rect",.3,.9,.3,parms[2][0],parms[2][1],parms[2][0]+parms[0],parms[2][1]+parms[0],1)))
			#self.wimage.updateGL()
		#elif self.calcmode==1:
			## update the box display on the image
			#nx=self.data["nx"]/parms[0]-1
			#shp={}
			#for x in range(nx):
				#for y in range(self.data["ny"]/parms[0]-1):
					## User deselected this one
					#if int(x+y*nx) in parms[3] : continue

					## Make a shape for this box
					#shp["box%02d%02d"%(x,y)]=EMShape(("rect",.3,.9,.3,(x+.5)*parms[0],(y+.5)*parms[0],(x+1.5)*parms[0],(y+1.5)*parms[0],1))

			#self.wimage.del_shapes()
			#self.wimage.add_shapes(shp)
			#self.wimage.updateGL()

		#if self.f2dmode>0 :
			#if self.f2dmode==1 : self.wfft.set_data(self.fft-self.fftbg)
			#else : self.wfft.set_data(self.fftbg)

		#else :
			#self.wfft.set_data(self.fft)
		#self.busy=False

	#def newCalcMode(self,mode):
		#self.calcmode=mode
		#self.recalc()

	#def new2DMode(self,mode):
		#self.f2dmode=mode
		#self.recalc()

	#def new2DAnMode(self,mode):
		#self.f2danmode=mode
		#self.needredisp=True

	#def newPlotMode(self,mode):
		#self.plotmode=mode
		#self.wplot.set_data(None,replace=True,quiet=True)	# clear the data so plots are properly redisplayed, but don't update the display
		#self.needredisp=True

	#def newBox(self):
		#parms=self.parms[self.curset]
		#parms[0]=self.sboxsize.value
##		parms[5]=self.soversamp.value
		#parms[5]=1
		#parms[3]=set()
		#self.recalc()

	#def newQualityFactor(self):
		#parms=self.parms[self.curset]
		#parms[4]=self.squality.value


	#def newCTF(self) :
		#parms=self.parms[self.curset]
		#parms[1].defocus=self.sdefocus.value
		#parms[1].bfactor=self.sbfactor.value
		#parms[1].dfdiff=self.sdfdiff.value
		#parms[1].dfang=self.sdfang.value
		#parms[1].apix=self.sapix.value
		#parms[1].ampcont=self.sampcont.value
		#parms[1].voltage=self.svoltage.value
		#parms[1].cs=self.scs.value
		#self.needredisp=True



	#def fftmousedown(self,event,m) :
		##m=self.wfft.scr_to_img((event.x(),event.y()))

		#if self.f2danmode==1:
			#self.ringrad=hypot(m[0]-self.fft["nx"]/2,m[1]-self.fft["ny"]/2)
			#self.needredisp=True
		#elif self.f2danmode==2:
			#if (event.modifiers()&Qt.ControlModifier): self.xpos2=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			#else: self.xpos1=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			#self.needredisp=True


		##self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])

	#def fftmousedrag(self,event,m) :
		##m=self.wfft.scr_to_img((event.x(),event.y()))

		#if self.f2danmode==1:
			#self.ringrad=hypot(m[0]-self.fft["nx"]/2,m[1]-self.fft["ny"]/2)
			#self.needredisp=True
		#elif self.f2danmode==2:
			#if (event.modifiers()&Qt.ControlModifier): self.xpos2=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			#else: self.xpos1=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			#self.needredisp=True
		## box deletion when shift held down
		##if event.modifiers()&Qt.ShiftModifier:
			##for i,j in enumerate(self.boxes):

	#def fftmouseup(self,event,m) :
		#"up"
		##m=self.wfft.scr_to_img((event.x(),event.y()))


	#def plotmousedown(self,event) :
		#"mousedown in plot"
##		m=self.guiim.scr_to_img((event.x(),event.y()))

#def tiled(img,box):
	#imgc=img.process("math.meanshrink",{"n":2})		# shrink image by 2 for boxing
	#boxc=good_boxsize(box/2,larger=True)			# box size in reduced image
	#nxb=4*imgc["nx"]/boxc-2		# number of boxes along x direction
	#nyb=4*imgc["ny"]/boxc-2
	#radius=boxc/2.6
	
	#mask1=EMData(boxc,boxc,1)
	#mask1.to_one()
	#mask1.process_inplace("mask.gaussian",{"outer_radius":radius,"exponent":4.0})
	## Mask 2 is the 'inverse' (1.0-val) of mask1, with the addition of a soft outer edge to reduce periodic boundary condition issues
	#mask2=mask1.copy()*-1+1
##		mask1.process_inplace("mask.decayedge2d",{"width":4})
	#mask2.process_inplace("mask.decayedge2d",{"width":2})
	##mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
	##mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
	
	## ratio1,2 give us info about how much of the image the mask covers for normalization purposes
	#ratio1=mask1.get_attr("square_sum")/(boxc*boxc)	#/1.035
	#ratio2=mask2.get_attr("square_sum")/(boxc*boxc)

	#cor=EMData(nxb,nyb)

	#vecs=[]
	#for y in xrange(nyb):
		#for x in xrange(nxb):
			#im1=imgc.get_clip(Region(boxc/8+x*boxc/4,boxc/8+y*boxc/4,boxc,boxc))
			#im1.process_inplace("normalize.edgemean")
			
			#im2=im1.copy()

	##		print im2.get_size(), im1.get_size()

			## now we compute power spectra for the 2 regions defined by the masks
			#im1.mult(mask1)
			#imf=im1.do_fft()
			#imf.ri2inten()
			#imf/=(imf["nx"]*imf["ny"]*ratio1)
			#cen_1d=imf.calc_radial_dist(imf.get_ysize()/2,0.0,1.0,1)
			
			#im2.mult(mask2)
			#imf=im2.do_fft()
			#imf.ri2inten()
			#imf/=(imf["nx"]*imf["ny"]*ratio2)
			#edg_1d=imf.calc_radial_dist(imf.get_ysize()/2,0.0,1.0,1)

			#vec=EMData(imf["ny"]/4-2,1,1)		# We skip the first 2 points and only go to 1/2 Nyquist
			#for i in xrange(2,imf["ny"]/4):
				#vec[i]=(cen_1d[i]-edg_1d[i])/edg_1d[i]		# expressed as a SSRN
			
			#vecs.append(vec)
			
			#img[(5*boxc/8+x*boxc/4)*2,(5*boxc/8+y*boxc/4)*2]=(vec["mean"]+0.5)*10.0
			#cor[x,y]=vec["sigma"]
	
	#cor.update()
	#return vecs,cor

#def detect(img,box):
	#img.process_inplace("normalize.edgemean")
	
	#radius=box/2.6
	#zro=img["ctf"].zero(0)*img["ctf"].apix*box
	#wvlen=box*2.0/zro	# 1/2 the wavelength of the 1st zero
##	wvlen=box/zro	# the wavelength of the 1st zero

	#mask1c=EMData(box,box,1)
	#mask1c.process_inplace("testimage.sinewave.circular",{"wavelength":wvlen,"phase":0.0})
	#mask1c.process_inplace("mask.gaussian",{"outer_radius":radius,"exponent":4.0})
##	mask1c.process_inplace("normalize.unitlen")
	#mask1c.process_inplace("normalize")
	#mask1c/=box*box
	#print mask1c["mean"],mask1c["sigma"],wvlen,zro
	#mask1c.clip_inplace(Region(-(img["nx"]-box)/2.0,-(img["ny"]-box)/2.0,img["nx"],img["ny"]))
	#mask1c.process_inplace("xform.phaseorigin.tocorner")

	#mask1s=EMData(box,box,1)
	#mask1s.process_inplace("testimage.sinewave.circular",{"wavelength":wvlen,"phase":pi/2.0})
	#mask1s.process_inplace("mask.gaussian",{"outer_radius":radius,"exponent":4.0})
##	mask1s.process_inplace("normalize.unitlen")
	#mask1s.process_inplace("normalize")
	#mask1s/=box*box
	#print mask1c["mean"],mask1c["sigma"],wvlen,zro
	#mask1s.clip_inplace(Region(-(img["nx"]-box)/2.0,-(img["ny"]-box)/2.0,img["nx"],img["ny"]))
	#mask1s.process_inplace("xform.phaseorigin.tocorner")

	#c1=img.calc_ccf(mask1c)
	#c1.process_inplace("math.squared")
	#s1=img.calc_ccf(mask1s)
	#s1.process_inplace("math.squared")
	#c1.add(s1)
	#c1.process_inplace("normalize")
	

	#display((c1,img),True)
	## Mask 2 is the 'inverse' (1.0-val) of mask1, with the addition of a soft outer edge to reduce periodic boundary condition issues
	##mask2=mask1.copy()*-1+1
###		mask1.process_inplace("mask.decayedge2d",{"width":4})
	##mask2.process_inplace("mask.decayedge2d",{"width":4})
	##mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
	##mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))

	#QWidget *firstPageWidget = new QWidget;
    #QWidget *secondPageWidget = new QWidget;
    #QWidget *thirdPageWidget = new QWidget;

    #QStackedWidget *stackedWidget = new QStackedWidget;
    #stackedWidget->addWidget(firstPageWidget);
    #stackedWidget->addWidget(secondPageWidget);
    #stackedWidget->addWidget(thirdPageWidget);

    #QVBoxLayout *layout = new QVBoxLayout;
    #layout->addWidget(stackedWidget);
    #setLayout(layout);
		
		 #QComboBox *pageComboBox = new QComboBox;
    #pageComboBox->addItem(tr("Page 1"));
    #pageComboBox->addItem(tr("Page 2"));
    #pageComboBox->addItem(tr("Page 3"));
    #connect(pageComboBox, SIGNAL(activated(int)),stackedWidget, SLOT(setCurrentIndex(int)));

		#self.lboxmode=QtGui.QLabel("Mode:",self)
		#self.gbl.addWidget(self.lboxmode,10,0)

		#self.sboxmode=QtGui.QComboBox(self)
		#self.sboxmode.addItem("Manual")
		#self.sboxmode.addItem("Reference")
		#self.sboxmode.setCurrentIndex(1)
		#self.gbl.addWidget(self.sboxmode,10,1)

		#self.lanmode=QtGui.QLabel("Annotate:",self)
		#self.gbl.addWidget(self.lanmode,12,0)

		#self.sanmode=QtGui.QComboBox(self)
		#self.sanmode.addItem("Box")
		#self.sanmode.addItem("Box+dot")
		#self.sanmode.addItem("Circle")
		#self.sanmode.addItem("None")
		#self.gbl.addWidget(self.sanmode,12,1)

		#self.sdefocus=ValSlider(self,(0,5),"Defocus:",0.0,90)
		#self.gbl.addWidget(self.sdefocus,0,2,1,3)

		#self.squality=ValSlider(self,(0,9),"Quality (0-9):",0,90)
		#self.squality.setIntonly(True)

		#self.gbl.addWidget(self.squality,6,2,1,3)

		#self.brefit=QtGui.QPushButton("Autobox")
		#self.gbl.addWidget(self.brefit,7,2)

		#self.bclrauto=QtGui.QPushButton("Clear Auto")
		#self.gbl.addWidget(self.bclrauto,7,3)

		#self.bclrall=QtGui.QPushButton("Clear All")
		#self.gbl.addWidget(self.bclrall,7,4)

		#self.sapix=ValBox(self,(0,500),"A/pix:",1.0,90)
		#if self.defaultapix!=None : self.sapix.setValue(self.defaultapix)
		#self.gbl.addWidget(self.sapix,10,2)

		#self.svoltage=ValBox(self,(0,500),"Voltage (kV):",200,90)
		#if self.defaultvoltage!=None : self.svoltage.setValue(self.defaultvoltage)
		#self.gbl.addWidget(self.svoltage,11,2)

		#self.scs=ValBox(self,(0,5),"Cs (mm):",4.1,90)
		#if self.defaultcs!=None : self.scs.setValue(self.defaultcs)
		#self.gbl.addWidget(self.scs,12,2)

		#self.sboxsize=ValBox(self,(0,500),"Box Size:",256,90)
		#self.sboxsize.intonly=True
		#self.gbl.addWidget(self.sboxsize,13,2)

		#self.sptclsize=ValBox(self,(0,500),"Ptcl Size:",256,90)
		#self.sptclsize.intonly=True
		#self.gbl.addWidget(self.sptclsize,14,2)

		#QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		#QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		#QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		#QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		#QtCore.QObject.connect(self.sboxsize, QtCore.SIGNAL("valueChanged"), self.newBox)
##		QtCore.QObject.connect(self.soversamp, QtCore.SIGNAL("valueChanged"), self.newBox)
		#QtCore.QObject.connect(self.squality,QtCore.SIGNAL("valueChanged"),self.newQualityFactor)
		#QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		#QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listkey)
		#QtCore.QObject.connect(self.sboxmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newBoxMode)

		#self.resize(720,380) # figured these values out by printing the width and height in resize event

		#### This section is responsible for background updates
		#self.busy=False
		#self.needupdate=True
		#self.needredisp=False
		#self.procthread=None
		#self.errors=None		# used to communicate errors back from the reprocessing thread

		#self.timer=QTimer()
		#QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		#self.timer.start(100)

#		self.recalc()


if __name__ == "__main__":
	main()
