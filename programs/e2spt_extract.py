#!/usr/bin/env python
# Muyuan Chen 2018-04
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2_utils import *
from EMAN2jsondb import JSTask
import os
import sys
import traceback

global thrdone
thrdone=0

class SptExtractTask(JSTask):
	
	
	def __init__(self, data, options, allxfs, allinfo):
		
		JSTask.__init__(self,"SptExtract",data,{},"")
		self.options=options
		self.allxfs=allxfs
		self.allinfo=allinfo	
	
	def execute(self, callback):
		
		
		for a in self.data:
			if a in self.allxfs:
				print("Extract from {} (with xform)".format(a))
				do_extraction(a, self.options, self.allxfs[a], self.allinfo[a])
			else:
				print("Extract from {}".format(a))
				do_extraction(a, self.options)
		return 1

def main():
	
	usage="""
	Extract 2D subtilt particles from the tilt series, and reconstruct 3D subvolumes.
	There are a number of modes for the program and it is probably easier to get the right commands from the GUI.
	
	To extract particles of a certain label from one tomogram
	e2spt_extract.py tomograms/tomox__bin4.hdf --boxsz_unbin 128 --label <label from particle picking> 
	
	To extract the particles of the same label from all tomograms
	e2spt_extract.py --alltomograms --boxsz_unbin 128 --label <label from particle picking> 
	
	To extract particles along curves drawn using e2tomo_drawcurve.py. The direction of particles along the curve will be written to the header of particles which can be used in downstream refinement. The spacing of particles along the curve is controlled by the box size and the overlap between adjacent boxes. 
	e2spt_extract.py --curves <curve ID> --curves_overlap 0.75 --alltomograms --boxsz_unbin 128 --newlabel <new particle label> 
	
	To extract particls from an existing refinement (either the particle_param_xx.json from e2spt_refine, or aliptcls3d_xx.lst from e2spt_refine_new)
	e2spt_extract.py --jsonali <json/lst from existing alignment> --boxsz_unbin 128 --newlabel <new particle label> 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="Specify tomograms from which you wish to extract boxed particles.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="extract")
	parser.add_argument("--boxsz_unbin", type=int,help="box size in unbinned tomogram", default=-1, guitype='intbox',row=2, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--label", type=str,help="Only extract particle with this name. Leave blank to extract all particles.", default=None, guitype='strbox',row=2, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--newlabel", type=str,help="Label of output particles. Same as original particle label by default.", default="", guitype='strbox',row=6, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=3, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=4, col=1,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--maxtilt", type=int,help="max tilt", default=100, guitype='intbox',row=4, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--padtwod", type=float,help="padding factor", default=2.0, guitype='floatbox',row=5, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--noctf", action="store_true", default=False ,help="skip ctf correction.")
	parser.add_argument("--wiener", action="store_true", default=False ,help="wiener filter the particles using ctf information. Not very useful...")
	parser.add_argument("--nomasksubtlt", action="store_true", default=False ,help="extracted subtilts should be unmasked")
	parser.add_argument("--alltomograms", action="store_true", default=False ,help="use all tomograms.", guitype='boolbox',row=1, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--dotest", action="store_true", default=False ,help="only make 1 batch of subtomograms for testing")

	parser.add_argument("--shrink", type=float, help="Shrinking factor for output particles. 1.5 or integers allowed. Default is 1 (no shrink).",default=1, guitype='floatbox',row=8, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--tltkeep", type=float,help="keep a fraction of tilt images with good score determined from tomogram reconstruction", default=1.0, guitype='floatbox',row=8, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--rmbeadthr", type=float,help="remove 2d particles with high contrast object beyond N sigma at 100A. Note that this may result in generating fewer particles than selected. Default is -1 (include all particles). 0.5 might be a good choice for removing gold beads but may need some testing...", default=-1, guitype='floatbox',row=9, col=0, rowspan=1, colspan=1, mode="extract")
	
	parser.add_header(name="orblock2", help='Just a visual separation', title="Extract from curves", row=10, col=0, rowspan=1, colspan=1, mode="extract")
	
	parser.add_argument("--curves", type=int, default=-1 ,help="specify curve id to extract particles from saved curves. ", guitype='intbox',row=11, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--curves_overlap", type=float, help="fraction of overlap when generating particle along curves. default is 0.5",default=0.5,  guitype='floatbox',row=12, col=0, rowspan=1, colspan=1, mode="extract")

	
	parser.add_header(name="orblock3", help='Just a visual separation', title="Re-extraction from spt", row=13, col=0, rowspan=1, colspan=1, mode="extract")

	parser.add_argument("--jsonali", type=str,help="re-extract particles using a particle_param_xx json file from a spt alignment", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=14, col=0,rowspan=1, colspan=2, mode="extract")
	parser.add_argument("--mindist", type=float,help="minimum distance between particles in A. for reextraction only", default=10)
	parser.add_argument("--keep", type=float,help="fraction of particles to keep fron previous alignment. for reextraction only.", default=1)
	parser.add_argument("--postproc", type=str,help="processor after 3d particle reconstruction", default="")
	parser.add_argument("--postmask", type=str,help="masking after 3d particle reconstruction. The mask is transformed if json ", default="")
	#parser.add_argument("--textin", type=str,help="text file for particle coordinates. do not use..", default=None)
	parser.add_argument("--compressbits", type=int,help="Bits to keep for compression. default is -1 meaning uncompressed floating point. 8 bit seems fine...", default=8,guitype='intbox',row=5, col=1, rowspan=1, colspan=1, mode="extract[8]")
	parser.add_argument("--norewrite", action="store_true", default=False ,help="skip existing files. do not rewrite.")
	parser.add_argument("--append", action="store_true", default=False ,help="append to existing files.")
	parser.add_argument("--parallel", type=str,help="parallel", default="")
	parser.add_argument("--postxf", type=str,help="a file listing post transforms (see http://eman2.org/e2tomo_more), or for simple symmetry, <sym>,<cx>,<cy>,<cz> where the coordinates specify the center of a single subunit", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)	
	parser.add_argument("--skip3d", action="store_true", default=False ,help="do not make 3d particles. only generate 2d particles and 3d header. ")
	parser.add_argument("--loadali2d", type=str,help="extract polished particles using existing subtilt alignment. Only works with --jsonali within the same refinement", default="")

	#parser.add_argument("--shrink3d", type=int, help="Only shrink 3d particles by x factor",default=-1)
	parser.add_argument("--verbose", type=int,help="verbose", default=1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	allxfs={}; allinfo={}
	
	if len(options.jsonali)>0:
		allxfs,allinfo=parse_json(options)
	
	if options.alltomograms:
		fld="tomograms/"
		args=[fld+f for f in os.listdir(fld) if f.endswith(".hdf")]
		### do not count a tomogram twice when multiple versions exist
		uq, uid=np.unique([base_name(a) for a in args], return_index=True)
		args=[args[i] for i in uid]
		
	if len(args)==0 and len(allxfs)==0:
		print("No input. Exit.")
		return
		
	if options.parallel!="":
		#### in theory we allow doing this via MPI in addition to threading, and it is possible to spawn one job per node, and each worker will do reconstruction with multiple threads using shared memory
		##   in practice it is not well tested and it is not quite clear how well it performs on different clusters. 
		##   since the major speed limiting factor of the program is still in disk I/O, it is not much slowerer with just --threads
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel, module="e2spt_extract.SptExtractTask")
		num_cpus = etc.cpu_est()
		options.nowtime=time.time()
		print("Running in parallel mode with {} workers, each running {} threads".format(num_cpus, options.threads))
		njob=num_cpus
		options.verbose=0
		if len(allxfs)>0:
			fnames=list(allxfs.keys())
		elif len(args)>0:
			fnames=args

		tids=[]
		for i in range(njob):
			task=SptExtractTask(fnames[i::njob], options, allxfs, allinfo)
			tid=etc.send_task(task)
			tids.append(tid)
		
		while 1:
			st_vals = etc.check_task(tids)
			if np.min(st_vals) == 100: break
			#print(st_vals)
			time.sleep(1)

	
	else:
		#### do extraction sequentially for each tomogram
		if len(allxfs)>0:
			for fname in allxfs.keys():
				#### it seems options is changed inplace somewhere...
				(options, args) = parser.parse_args()
				do_extraction(fname, options, allxfs[fname],allinfo[fname])
				
		elif len(args)>0:
			for a in args:
				do_extraction(a, options)
	
		return	
	
	print("Particle extraction finished.")
	
	
	E2end(logid)
	return
	
	
### info: extra header information per particle in a dictionary
def do_extraction(pfile, options, xfs=[], info=[]):
	#pfile=args[0]
	if options.verbose>0:print("Reading from {}...".format(pfile))
	if options.boxsz_unbin>0:
		options.boxsz= good_size(int(options.boxsz_unbin / options.shrink))
	else:
		options.boxsz=-1
	
	#### read tilt series alignment info...
	js=js_open_dict(info_name(pfile))
	try: 
		ttparams=np.array(js["tlt_params"])
		ttparams[:,:2]/=options.shrink
	except:
		traceback.print_exc()
		print(f"""\nERROR: tlt_params missing in {info_name(pfile)}. The most likely cause for this is trying to use a tomogram
reconstructed in another program. Subtomogram averaging in EMAN2 makes extensive use of the original tilt series
as well as alignment parameters, so use of tomograms from other software is not supported in the normal sequence.""")
		sys.exit(1)

	#### read ctf info if ctf estimation has been performed
	if options.noctf==False and "defocus" in js:
		#### read ctf info when exist
		defocus=np.array(js["defocus"])
		phase=np.array(js["phase"])
		voltage=float(js["voltage"])
		cs=float(js["cs"])
		if options.verbose>0:print("CTF information exists. Will do phase flipping..")
	else:
		defocus=[]
		phase=[]
		
	#### sort tilt image by alignment score so we can throw away the bad ones later 
	tltkeep=[]
	if options.tltkeep<1.0:
		if js.has_key("ali_loss"):
			aliloss=js["ali_loss"]
			tltkeep=np.argsort(aliloss)[:int(len(ttparams)*options.tltkeep)]
		else:
			print("warning: --tltkeep specified, but cannot find ali_loss in info. setting tltkeep to 1.")
			
	tfile=str(js["tlt_file"])
	options.tltfile=tfile
	ptcl=EMData(pfile, 0, True)
	img=EMData(tfile,0, True)
	apix_ptcl=ptcl["apix_x"]
	apix_tlt=img["apix_x"] * float(options.shrink)
	print("Extracting from tilt series {}...".format(tfile))
	
	#### scale factor from raw tilt image to tomogram/particle that we get coordinates info from
	scale=apix_ptcl/apix_tlt
	options._scale=scale
	if options.verbose>0:print("### scale by {:.1f}".format(options._scale))
	if options.shrink==1.5:
		shrinklab="_bin1_5"
	elif options.shrink>=2:
		shrinklab="_bin{:d}".format(int(options.shrink))
	else:
		shrinklab=""
	
	#### reading particle location from json file corresponding to the tomogram input
	ptclpos=[]
	nptcl=EMUtil.get_image_count(pfile)
	e=EMData(pfile, 0, True)
	
	#### overall, the paragraph below should generate a list of lists called towrite
	##   each entry contains (coordinates/transform of particles, output file name, box size/2, extra info)
	##   for some strange reasons, we use box size/2 here. changing it will need more testing so I am leaving it as is for now...
	
	if len(xfs)>0:
		#### reading particle coordinates from existing spt alignment
		##   coordinates and transforms should already be parsed and saved in xfs
		##   here just format them to be consistent with the standard protocol
		if options.newlabel=="":
			lab=pfile[pfile.rfind('__')+2:pfile.rfind('.')]
			lab+="_reextract"
		else:
			lab=options.newlabel
			
		outname=str(base_name(pfile)+"__"+lab+".hdf")
		sz=int(options.boxsz//2)
		if scale!=1:
			for t in xfs:
				t.set_trans(t.get_trans()*scale)
		towrite=[[xfs, outname, sz, info]]
		
	
	elif nptcl==1:
		#### extracting from a tomogram (instead of a particle stack)
		##   either along curves or as individual particles
		if options.verbose>0:print("Reading particle location from a tomogram...")
		js=js_open_dict(info_name(pfile))
		towrite=[]
		
		if options.curves>=0:
			#### particle extraction along curves
			if options.verbose>0:print("Generating particles along curves...")
			info=[]
			overlap=options.curves_overlap
			if overlap>=1 or overlap<0:
				print("Error: Overlap has to be in [0,1)")
				return
			
			if js.has_key("curves") and len(js["curves"])>0:
				pts=np.array(js["curves"]).copy()
				#### in the latest version, curves from the json file should include 5 columns
				##   (x, y, z, curve id, class id)
				##   the class id column does not exist in older versions so we keep this for compatibility
				if len(pts[0])<5:
					pts=np.hstack([pts, np.zeros((len(pts),1))])
				
				#### and we only extract curves with class id specified in the option
				pts=pts[pts[:,4]==options.curves, :4]
				
				if len(pts)>1:
					#### apix_unbin should always exist in js now
					##   just keep this to be safe...
					if "apix_unbin" in js:
						pts[:,:3]/=options.shrink
					else:
						## this is to deal with the old metadata format...
						pts[:,:3]-=[e["nx"]//2, e["ny"]//2, e["nz"]//2]
						pts[:,:3]*=scale
					
					#### since curves do not have named label, it is better to provide one
					if options.newlabel=="":
						lab="curve_{:02d}".format(options.curves)
					else:
						lab=options.newlabel
						
					lab+=shrinklab
					outname=str(base_name(pfile)+"__"+lab+".hdf")
					
					
					####   distance between adjacent particles is controlled by box size and overlap
					sz=int(options.boxsz//2)
					spacing=options.boxsz*(1-overlap)*apix_tlt
					curves=np.unique(pts[:,3])
					
					bxs=[]
					drs=[]
					for li in curves:
						## points on one curve
						pt=pts[pts[:,3]==li][:,:3].copy()
						## curve with a single point shouldn't exist anyway
						if len(pt)<2: continue
					
						## total length in A
						ln=np.sum(np.linalg.norm(np.diff(pt, axis=0), axis=1))*apix_tlt
						
						## interpolate to the given spacing
						if ln<spacing*3:continue
					
						pt=np.unique(pt, axis=0) ## duplicated points sometimes happen...
						ipt=interp_points(pt, npt=int(np.round(ln/spacing)))
						## again exclude curves that are too short
						if len(ipt)<4: continue
						
						## skip one point on the edge so we have better orientation
						cc=(ipt[:-1]+ipt[1:])/2.
						dd=ipt[:-1]-ipt[1:]
						
						## also keep the index of points along the curve
						## not used anywhere but in case it become useful...
						od=np.arange(len(cc),dtype=float)/(len(cc)-1)
						info.extend([{"cv_id":li, "cv_len":ln, "cv_pos":o} for o in od])
						
						bxs.append(cc)
						drs.append(dd)
					
					#### finally stack them together
					##   each particle now has 6 columns
					##   first 3 for coordinates, next 3 for direction
					bxs=np.vstack(bxs)
					drs=np.vstack(drs)
					bxs=np.hstack([bxs,drs])
					
					if options.verbose>0:print(" {:d} curves, {:d} points with {:.1f}A spacing".format(len(curves), len(bxs), spacing))
					
					towrite.append((bxs, outname, sz, info))
		
		elif "class_list" in js and "boxes_3d" in js:
			
			#### extract individual particles
			clslst=js["class_list"]
			boxes=js["boxes_3d"]
			for ky in list(clslst.keys()):
				val=clslst[ky]
				## if a label is provided, only take particles of the given label
				## otherwise it can extract all particles of different labels and write to different stacks
				if options.label:
					if str(val["name"])!=options.label:
						continue
						
				bxs=np.array([[b[0], b[1], b[2]] for b in boxes if b[5]==int(ky)], dtype=float)
				if len(bxs)==0: continue
				
				## apix_unbin should always exist now
				if "apix_unbin" in js:
					bxs/=options.shrink
					## use the box size from the spt_boxer if a box size is not provided through the option
					if options.boxsz<0:
						sz=int(val["boxsize"])//2
					else:
						sz=int(options.boxsz//2)
				else:
					bxs-=[e["nx"]//2, e["ny"]//2, e["nz"]//2]
					bxs*=scale
					if options.boxsz<0:
						sz=int(val["boxsize"])*scale//2
					else:
						sz=int(options.boxsz//2)
						
				## output label
				if options.newlabel=="":
					lab=val["name"]
				else:
					lab=options.newlabel
					
				lab+=shrinklab
				outname=str(base_name(pfile)+"__"+lab+".hdf")
				towrite.append((bxs, outname, sz, info))
				
				if options.verbose>0:print("{} : {} boxes, output box size {}".format(val["name"], len(bxs), int(sz*2))) 
		
		js.close()
		
		if len(towrite)==0:
			print("No particles. exit..")
			return
		
	else:
		#### extracting from existing particle stack
		##   this should be quite useless in the current workflow now...		
		if options.verbose>0:print("Reading particle location from image stack..")
		for i in range(nptcl):
			ptcl=EMData(pfile, i, True)
			ptclpos.append(ptcl["ptcl_source_coord"])
		
	
		sfx=pfile[pfile.find("__"):]
		for i in range(2,10): sfx=sfx.replace("_bin{:d}".format(i),"")
		
		sfx+=shrinklab
		outname=base_name(pfile)+sfx+".hdf"
		ptclpos=np.array(ptclpos, dtype=float)
		
		ptclpos*=scale
		
		if options.boxsz<0:
			boxsz=int(ptcl["nx"]//2*scale)
		else:
			boxsz=int(options.boxsz//2*scale)
		
		towrite=[(ptclpos,outname, boxsz,[])]
	
	#### skip existing particle files with the same name
	##   mostly for debugging though. there is a risk that particles with the same name may not have the same content. so still better to overwrite everything
	##   here we skip entire tilt series if needed. we will skip individual particle files later..
	if options.norewrite:
		skip=True
		for pinfo in towrite:
			ptclpos, outname, boxsz, info=pinfo
			output=os.path.join("particles3d", outname)
			if not os.path.isfile(output):
				skip=False
		if skip:
			print("all particles in {} exist. skip tilt series...".format(tfile))
			
	#### now actually read in the tilt series
	##   make sure this works on image stack or mrc volume
	img=EMData(tfile,0)
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	#### shrink the tilt series before extracting particles if shrink is specified
	for m in imgs: 
		if options.shrink>1:
			m.process_inplace("math.meanshrink",{"n":options.shrink})
		m.process_inplace("normalize")
	ntlt=len(imgs)
	
	#### optionally apply a customized mask to the particles after extraction
	if options.postmask!="":
		pmask=EMData(options.postmask)
	else:
		pmask=None
		
	try: os.mkdir("particles3d")
	except: pass
	try: os.mkdir("particles")
	except: pass
	
	#### now prepare for particle writing...
	for pinfo in towrite:
		ptclpos, outname, boxsz, info=pinfo
		if len(info)>0 and len(info)!=len(ptclpos):
			print("Warning: Extra header info exist but does not match particle count...")
			
		# if a particle has been removed from the JSON file it causes issues if we don't skip it
		if len(info)==0 and len(options.jsonali)>0: 
			continue

		nptcl=len(ptclpos)
	
		options.output=os.path.join("particles3d", outname)
		options.output2d=os.path.join("particles", outname)
		options.pad=pad=good_size(boxsz*2*options.padtwod)
		
		if options.norewrite:
			if os.path.isfile(options.output):
				print("File {} already exist. skipping...".format(options.output))
				continue
				
		
		print("Writing {} particles to {}".format(nptcl, outname))
		if options.verbose>0:print("Box size {}, padded to {}".format(int(boxsz*2), pad))
		
		## there really shouldn't be a case where append is necessary...
		if options.append:
			print(" appending to existing file")
		else:
			try: os.remove(options.output)
			except: pass
			try: os.remove(options.output2d)
			except: pass

		jsd=queue.Queue(0)
		jobs=[]
		
		## prepare ctf info
		if len(defocus)>0:
			ctf=[defocus, phase, voltage, cs]
		else:
			ctf=[]
		
		#### do the actural reconstruction in threads
		for tid in range(options.threads):
			ids=list(range(tid, nptcl, options.threads))
			jobs.append([jsd, ids, imgs, ttparams, pinfo, options, ctf, tltkeep, pmask])
			
		thrds=[threading.Thread(target=make3d,args=(i)) for i in jobs]
		global thrdone
		thrdone=0
		ndone=0
		time0=time.time()
		for t in thrds:
			t.start()
			
		while thrdone<len(thrds) or not jsd.empty():
			while not jsd.empty():
				
				#### each worker returns the (particle id, 3d subtomogram, list of 2d subtilts)
				##   write them to disk in main thread
				pid, threed, projs=jsd.get()
				if options.append:
					pid=-1

				#### we append the 2d particles to the end of particle stack file
				##   only keep the index of 3d particles, and save the indices of 2d particle 
				##   in the header of the corresponding 3d particle
				try: pji=EMUtil.get_image_count(options.output2d)
				except: pji=0
				pjids=[]
				for i,pj in enumerate(projs):
					if options.compressbits<0 : pj.write_image(options.output2d, pji)
					else: pj.write_compressed(options.output2d,pji,options.compressbits,nooutliers=True)
					pjids.append(pji)
					pji+=1
					
				threed["class_ptcl_src"]=options.output2d
				if len(pjids)>1:
					threed["class_ptcl_idxs"]=pjids
				elif pid>=0:
					print("Empty particle exist. Consider use the --append option. Exit.")
					exit()
				else:
					print("Empty particle detected. Skipping particle.")
					continue

				if options.compressbits<0: threed.write_image(options.output, pid)
				else: threed.write_compressed(options.output, pid,options.compressbits,nooutliers=True)
				
				ndone+=1
				if options.verbose>0:
					sys.stdout.write("\r{}/{} finished.".format(ndone, nptcl))
					sys.stdout.flush()
			time.sleep(.2)

		for t in thrds: t.join()
			
		print("Done. Particles written to {} ({:.1f} s)".format(options.output, time.time()-time0))
	
	
#### the actual worker function to extract 2d particles and reconstruct 3d particle
##   jsd: queue object for returning results from parallelism
##   imgs: tilt series input
##   ids: indices of particles to reconstruct
##   pinfo: particle info (coordinates, output name, box size/2, extra header info)
def make3d(jsd, ids, imgs, ttparams, pinfo, options, ctfinfo=[], tltkeep=[], mask=None):
	ppos, outname, boxsz, info=pinfo
	if len(info)!=len(ppos):
		info=[]
	
	pad=options.pad
	apix=imgs[0]["apix_x"]

	p3d=pad
	bx=boxsz*2
	apixout=apix
		
	if options.verbose>0: print(f"Begin make3d bx={bx} pad={pad} ")
	
	#### prepare the reconstructor
	recon=Reconstructors.get("fourier", {"sym":'c1', "size":[p3d, p3d, p3d], "mode":"trilinear"})
	recon.setup()
	
	#### parse ctf info
	if len(ctfinfo)>0:
		defocus, phase, voltage, cs=ctfinfo
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":voltage, "bfactor":0., "cs":cs,"ampcont":0, "apix":apix})
	
	for pid in ids:
		## redo the reconstructor setup so we reuse it for each particle
		recon.setup()
		pos=ppos[pid]
		
		## the header info will be added to the header of 3d particles
		dxfs=[]
		if len(info)>0:
			hdr=info[pid]
			if "dxfs2d" in hdr:
				dxfs=hdr.pop("dxfs2d")
				#x=np.array([x.get_trans() for x in dxfs])
				#print(pid, x.shape, np.max(abs(x), axis=0))
		else:
			hdr={}
		
		#print(len(dxfs), hdr)
		#### 3 modes for particle coordinates input
		if type(pos)==type(Transform()):
			## xform object from spt allignment json/lst
			tf_dir=Transform(pos.get_rotation()).inverse()
			pos=np.array(pos.get_trans())
			
		elif len(pos)>3:
			## 6 columns list - the last 3 columes are the direction of the curve
			drs=pos[3:].copy()
			pos=pos[:3].copy()
			v=np.linalg.norm(drs)
			if v==0: continue
			drs/=v
			drs=drs.tolist()
			tf_dir=Transform()
			try: tf_dir.set_rotation(drs)
			except: continue
			tf_dir.invert()

		else:
			## default mode, only coordinates from tomograms
			tf_dir=None
		
		#### now loop through images in the tilt series
		projs=[]
		for nid, m in enumerate(imgs):
			## skip tilt for tiltkeep and maxtilt
			if len(tltkeep)>0:
				if nid not in tltkeep:
					continue
			
			tpm=ttparams[nid]
			if abs(tpm[3])>options.maxtilt: continue
			
			## same function copied from e2tomogram
			## get 2d coordinates of a particle from tilt series alignment and 3d particle coordinates
			pxf=get_xf_pos(ttparams[nid], pos)
			tx=m["nx"]//2 +pxf[0]
			ty=m["ny"]//2 +pxf[1]

			## get_clip only take integers and we do not want to introduce the extra interpolation for floating point translation
			txint=int(tx)
			tyint=int(ty)
			
			## still save the extra floating point translation in the header of the 2d particle (xform.projection), so it is actually considered in 3d  reconstruction
			txdf=tx-txint
			tydf=ty-tyint

			e=m.get_clip(Region(txint-pad//2, tyint-pad//2, pad, pad), fill=0)
			
			## particles always have positive protein density so we invert the contrast here
			e.mult(-1)
			e.process_inplace("normalize.edgemean")
			
			## cut off high intensity pixels
			e.process_inplace("threshold.belowtozero",{"minval":-8})
			e.process_inplace("threshold.abovetozero",{"maxval":8})
			
			## add a soft falloff
			wd=(pad-boxsz*2)/4.
			if not options.nomasksubtlt: e.process_inplace("mask.soft",{"outer_radius":-wd, "width":wd/2})
			
			## skip all zero particles
			if e["sigma"]==0:
				continue
			
			## removing fiducial. doing this a bit lazy here by filtering and thresholding. really should save fiducial coordinates from tomogram reconstruction and remove here as  well
			if options.rmbeadthr>0:
				er=e.process("filter.lowpass.gauss",{"cutoff_freq":0.01})
				er.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
				er.process_inplace("threshold.binary",{"value":options.rmbeadthr})
				er.process_inplace("mask.addshells.gauss",{"val1":e["nx"]//40,"val2":e["nx"]//25})
				
				#if er["mean"]>0.1: 
					#### too many beads. better just remove this subtilt
					#continue
				
				e.mult(1-er)
				e.mult(1-er)
				noise=e.copy()
				noise.to_zero()
				noise.process_inplace("math.addnoise",{"noise":1})
				noise.process_inplace("normalize")
				e.add(noise*er)
				
			#### compute depth of particle in ice for defocus 
			rot=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
			p1=rot.transform(pos.tolist())
			pz=p1[2]
			dz=pz*apix/10000.
			
			if len(ctfinfo)>0:
				## phase flipping
				df=defocus[nid]-dz
				ctf.set_phase(phase[nid]*np.pi/180.)
				ctf.defocus=df
				e["ctf"]=ctf
				fft1=e.do_fft()
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				
				if options.wiener:
					## wiener filter by CTF
					## does not seem to be very useful
					## but does not hurt either
					ctf.dsbg=1./(e["apix_x"]*e["nx"])
					s=np.array(ctf.compute_1d(e["nx"], ctf.dsbg, Ctf.CtfType.CTF_INTEN ))
					s[:np.argmax(s)]=np.max(s)
					s=np.maximum(s, 0.0001)
					ctf.snr=s.tolist()
					flt=fft1.copy_head()
					ctf.compute_2d_complex(flt,Ctf.CtfType.CTF_WIENER_FILTER)
					fft1.mult(flt)
				
				e=fft1.do_ift()
				e["ctf"]=ctf
			
			#### save the metadata in header of 2d particle
			xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4], "ztilt":tpm[2], "tx":txdf, "ty":tydf})
			if len(dxfs)==len(imgs):
				dxf=Transform(dxfs[nid])
				#xform.set_trans(xform.get_trans()*options._scale)
				dxf.set_trans(dxf.get_trans()*options._scale)
				e["orig_dxf"]=dxf
				#print(nid, dxf)
				xform=dxf*xform
				txdf,tydf,_=xform.get_trans()
				
				
			e["xform.projection"]=xform
			e["ptcl_source_coord"]=[float(txint), float(tyint), float(dz)]
			e["ptcl_source_src"]=options.tltfile
			e["model_id"]=pid
			e["tilt_id"]=nid
			e["file_threed"]=options.output
			e["ptcl_source_coord_3d"]=pos.tolist()
			#e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.45})
			projs.append(e)
			
			#### pad it and insert to reconstructor
			e0=e
			sz=e0["nx"]
			e0=e0.get_clip(Region((sz-p3d)//2,(sz-p3d)//2,p3d,p3d))
			trans=Transform({"type":"2d", "tx":-txdf, "ty":-tydf})
			if not options.skip3d:
				e1=recon.preprocess_slice(e0, trans)
				recon.insert_slice(e1,xform,1)
		
		#### for testing purposes
		if options.skip3d:
			threed=EMData(1,1,1)
			
		elif len(projs)<len(imgs)/5:
			#### too many bad 2D particles
			##   better to just skip the 3d particle 
			##   not sure if this is a good idea yet...
			threed=EMData(bx,bx,bx)
			threed.to_zero()
			#continue
			
		else:
			#### finish reconstruction and do some postprocesing
			threed=recon.finish(True)
			#if mask: threedraw=threed.copy()
			
			try: threed=threed.get_clip(Region((p3d-bx)//2,(p3d-bx)//2,(p3d-bx)//2,bx,bx,bx))
			except:
				traceback.print_exc()
				print(f"get_clip(Region({(p3d-bx)//2},{(p3d-bx)//2},{(p3d-bx)//2},{bx},{bx},{bx}))")
				sys.exit(1)
		
			threed.process_inplace("normalize.edgemean")
			#projs=[pj.get_clip(Region((p3d-bx)//2,(p3d-bx)//2,bx,bx)) for pj in projs]
		
			if options.postproc!="":
				(filtername, param_dict) = parsemodopt(options.postproc)
				threed.process_inplace(filtername, param_dict)
				
			if mask:
				#### note that when a transform exist, we rotate the mask and apply to the particle
				##   so the missing wedge always points to z axis
				try:
					if tf_dir:
						mskrot=mask.process("xform", {"transform":tf_dir.inverse()})
					else:
						mskrot=mask.copy()
						
					threed.mult(mskrot)
					
					pjnew=[]
					for pj in projs:
						
						tpmsk=threed.project("standard", pj["xform.projection"])
						tpmsk.set_attr_dict(pj.get_attr_dict())
						pjnew.append(tpmsk)
					
					projs=pjnew
					
				except: 
					print(f'Mask error, size mismatch volume->{threed["nx"]} mask->{mask["nx"]}\n{pinfo}')
					raise Exception()
				
		#### get the metadata in header
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apixout
		threed["ptcl_source_coord"]=pos.tolist()
		threed["file_twod"]=options.output2d
		threed["model_id"]=pid
		for hd in hdr.keys():
			threed[str(hd)]=hdr[hd]
			
		if tf_dir:
			threed["xform.align3d"]=tf_dir
		
		
		jsd.put((pid, threed, projs))
		#recon.clear()

	global thrdone
	thrdone+=1
	return
	


#### get 2D position on a tilt given 3D location
def get_xf_pos(tpm, pk):
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0], pk[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)
	
	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0], p1[1]])
	

	return [p2[0], p2[1]]


#### parse a text file for particles
def parse_text(options):
	js=js_open_dict(options.textin)
	keys=natural_sort(js.keys())
	allinfo={}
	allxf={}
	for ky in keys:
		ptcls=js[ky]
		xfs=[]
		info=[]
		for p in ptcls:
			xf=p[3]
			xf.translate(p[0], p[1], p[2])
			xfs.append(xf)
			if len(p)>4:
				info.append(p[4])
			
		allxf[str(ky)]=xfs
		allinfo[str(ky)]=info
		
	js.close()
	return allxf,allinfo


#### parse a json file from spt_align to re-extract particles
def parse_json(options):
	print("re-extracting particles based on previous alignment. Ignoring particle/tomogram input...")
	import scipy.spatial.distance as scipydst
	
	if options.jsonali.endswith("json"):
		### read json
		js=js_open_dict(options.jsonali)
		keys=natural_sort(js.keys())
		data=[]
		for k in keys:
			src, ii = eval(k)
			e=EMData(src, ii, True)
			dt={"k":k,
				"pos":e["ptcl_source_coord"],
				"xf":js[k]["xform.align3d"], 
				"score":js[k]["score"],
				"class":ii%2,
				}
			if src[-4:]==".lst" :
				dt.update({"src":e["data_source"], "srci":e["data_n"]})
				
			else:
				dt.update({"src":src,"srci":ii})
				
			data.append(dt)
		js.close()
	
	elif options.jsonali.endswith("lst"):
		info=load_lst_params(options.jsonali)
		data=[]
		e=EMData(info[0]["src"], info[0]["idx"], True)
		for i,nf in enumerate(info):
			dt={"src":nf["src"],
				"srci":nf["idx"],
				"xf":nf["xform.align3d"],
				"score":-1,
				"class":i%2}
			
			if "coord" in nf:
				dt["coord"]=nf["coord"]
			else:
				e=EMData(nf["src"], nf["idx"], True)
				dt["coord"]=e["ptcl_source_coord"]
			
			if "score" in nf:
				dt["score"]=nf["score"]
				
			if "class" in nf:
				dt["class"]=nf["class"]
				
			data.append(dt)
			
		if options.loadali2d:
			dr=os.path.dirname(options.loadali2d)
			if dr!=os.path.dirname(options.jsonali):
				print("--loadali2d and --jsonali must from the same spt_xx folder")
				exit()
			
			info3d=load_lst_params("{}/particle_info_3d.lst".format(dr))
			ali2d=load_lst_params(options.loadali2d)
			for i,i3 in enumerate(info3d):
				i2d=i3["idx2d"]
				data[i]["dxfs2d"]=[ali2d[i2]["dxf"] for i2 in i2d]
				
	
	fs=[d["src"] for d in data]
	fnames, count=np.unique(fs, return_counts=True)

	score=np.array([d["score"] for d in data])
	print("Reading {} particles from {} tomograms. Score from {:.2f} to {:.2f}".format(len(score),len(fnames), min(score), max(score)))
	if options.keep<1:
		sthr=np.sort(score)[int(len(score)*options.keep)+1]
		print("Removing particles with score above {:.2f}.".format(sthr))
	else:
		sthr=1000
	
	### parse post transform file (or inline symmetry option)
	xffile=options.postxf
	postxfs=[]
	if xffile!=None:
		if os.path.isfile(xffile):
			f=open(xffile,'r')
			lines=f.readlines()
			f.close()
			for l in lines:
				if len(l)>3:
					postxfs.append(Transform(eval(l)))
		else:
			try:
				xfsym,xc,yc,zc=xffile.split(",")
				xc=float(xc)
				yc=float(yc)
				zc=float(zc)

				sym=Symmetries.get(xfsym)
				xf0=Transform({"type":"eman", "tx":xc, "ty":yc, "tz":zc})
				for i in range(sym.get_nsym()):
    					postxfs.append(sym.get_sym(i)*xf0)
					
			except:
				print("ERROR: --postxf must be a filename or <sym>,<cx>,<cy>,<cz>")
				sys.exit(1)

		print("Extracting {} sub-particles per original particle".format(len(postxfs)))
	else:
		postxfs.append(Transform())
		
	allxfs={}
	allinfo={}
	dthr=options.mindist/e["apix_x"]
	print("Removing particles within {:.1f} pixels from nearest neighbor".format(dthr))
	nptcl=0
	nexclude=0
	alld=[]
	for fi,fname in enumerate(fnames):
			
		dat=[d for d in data if d["src"]==fname]
		pos=np.array([d["coord"] for d in dat])
		ptclxfs=[Transform(d["xf"]) for d in dat]
		score=np.array([d["score"] for d in dat])
		if options.loadali2d:
			dxfs=[d["dxfs2d"] for d in dat]
		else:
			dxfs=[[] for d in dat]
		
		info=[]
		newxfs=[]
		for i, pxf in enumerate(ptclxfs):
			if score[i]>sthr: continue
			c=pos[i]
			for nxf in postxfs:
				ali=Transform(pxf)
				ali=nxf.inverse()*ali
				a=ali.inverse()
				a.translate(c[0], c[1], c[2])
				newxfs.append(a)
				info.append({"orig_ptcl":str(fname),"orig_idx":int(dat[i]["srci"]),"orig_xf":pxf, "orig_class":dat[i]["class"], "dxfs2d":dxfs[i], "postxf":nxf})
				#x=np.array([x.get_trans() for x in dxfs[i]])
				#print(i, np.max(abs(x), axis=0))

		newpos=np.array([p.get_trans() for p in newxfs])
		if len(newpos)>1:
			dst=scipydst.cdist(newpos, newpos)+(np.eye(len(newpos))*1e5)
			tokeep=np.ones(len(dst), dtype=bool)
			for i in range(len(dst)):
				if tokeep[i]:
					tokeep[dst[i]<dthr]=False

			newpos=newpos[tokeep]
			newxfs=[xf for i, xf in enumerate(newxfs) if tokeep[i]]
			info=[p for i, p in enumerate(info) if tokeep[i]]
			nexclude+=np.sum(tokeep==False)
			
		nptcl+=len(newpos)
		
		allxfs[fname]=newxfs
		allinfo[fname]=info
		
		
		print("{} : ptcls {} -> {}".format(fname, len(pos), len(allxfs[fname])))
			
			
	
	print("Writing {} particles, excluding {} particles too close to each other.".format(nptcl, nexclude))
	return allxfs,allinfo
	

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()

