#!/usr/bin/env python
# Muyuan Chen 2017-04
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize
import scipy.spatial.distance as scidist
from EMAN2_utils import *
from multiprocessing import Pool

def main():
	
	usage="""WARNING: This is an experimental program that may or may not work..
	
	This program takes an unaligned tilt series, performs fiducial (or fiducial-less) alignment, generate a tomogram and optionaly removes gold fiducials. Currently it only works on bin4 images smaller than 1Kx1K in size. Again this is an experimental program that has minimal functionality of tilt series alignment and reconstruction, and should only be used by people that are too lazy to use other software packages..
	
	The program performs a rough alignment and generate a tomogram, then look for high density landmarks in the 3D volume, track them in 2D tilt series, refine the alignment and reconstruct the 3D volume again. This process is repeated a number of time until a final tomogram is generated.
	
	Usage:
	e2tomogram.py <tilt series stack> [(--rawtlt) or (--tiltstep --zeroid)] [options]
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--rawtlt", type=str,help="Text file contains raw tilt angles", default=None)
	parser.add_argument("--loadparam", type=str,help="Load from existing param file", default=None)
	parser.add_argument("--zeroid", type=int,help="Index of the center tilt. Ignored when rawtlt is provided.", default=-1)
	parser.add_argument("--npk", type=int,help="Number of landmarks to use.", default=20)
	parser.add_argument("--tltstep", type=float,help="Step between tilts. Ignored when rawtlt is provided.", default=2.)
	parser.add_argument("--tltax", type=float,help="Angle of the tilt axis. The program will calculate one if this option is not provided", default=None)
	parser.add_argument("--pkeep", type=float,help="Fraction of landmarks to keep in the tracking.", default=.5)
	parser.add_argument("--tltkeep", type=float,help="Fraction of tilts to keep in the reconstruction.", default=.9)
	parser.add_argument("--minloss", type=float,help="Stop refinement when the loss is lower than this value.", default=1.)
	parser.add_argument("--bxsz", type=int,help="Box size of the particles for tracking", default=-1)
	parser.add_argument("--writetmp", action="store_true",help="Write intermidiate files", default=False)
	parser.add_argument("--rmgold", action="store_true",help="Remove gold fiducials.", default=False)
	parser.add_argument("--nofiducial", action="store_true",help="Fiducial-less mode. This will change a few internal parameters to make it work.", default=False)
	parser.add_argument("--reconmode", type=str,help="Reconstruction mode. Choose from nearest_neighbor, gauss_2, gauss_3, and gauss_5.", default="gauss_2")
	parser.add_argument("--threads", type=int,help="Number of threads", default=12)
	parser.add_argument("--niter", type=int,help="Number of iterations", default=3)
	parser.add_argument("--verbose", type=int,help="Verbose", default=0)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	inputname=args[0]
	
	dotpos=inputname.rfind('.')
	linepos=inputname.rfind('__')
	inputtag=inputname[linepos:dotpos]
	bname=base_name(inputname)
	
	options.basename=bname
	
	if options.bxsz<0:
		options.bxsz=96
		#if options.nofiducial:
			#options.bxsz=96
		#else:
			#options.bxsz=64
	
	
	print "Reading and pre-processing tilt stack..."
	imgs=EMData.read_images(inputname)
	for m in imgs:
	
		m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
		m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
		m.process_inplace("normalize.edgemean")

	num=len(imgs)
	options.num=num
	
	#### parse tilt info
	if options.rawtlt:
		#### read from rawtlt file
		tlts=np.loadtxt(options.rawtlt)
		options.zeroid=np.argmin(abs(tlts))
		#### here we always assume the center tilt is at 0 degree
		tlts-=tlts[options.zeroid]
		options.tltstep=abs(np.mean(np.diff(tlts)))
	else:
		#### generate tilt angles from tilt step size
		if options.zeroid<0:
			options.zeroid=int(num/2)
		tlts=np.arange(num, dtype=float)*options.tltstep
		tlts-=tlts[options.zeroid]
	
	tlts=tlts[:num]
	print "tilt angle from {:.1f} to {:.1f}, step {:.1f}".format(np.min(tlts), np.max(tlts), options.tltstep)
	options.tlt_init=tlts.copy()

	nx=imgs[0]["nx"]/2
	ny=imgs[0]["ny"]/2
	options.minsz=min(nx, ny)*2
	
	#### make a folder to write tmp files
	if options.writetmp:
		for i in range(100):
			try:
				path="tomorecon_{:02d}/".format(i)
				os.mkdir(path)
				options.tmppath=path
				break
			except:
				continue
		else:
			print "Too many tomorecon folders in the project, or something odd happened....Exit."
			exit()
			
		print "Temporary files will be written in {}".format(options.tmppath)
	
	if options.loadparam:
		#### load parameters from saved file
		print "Loading parameters from {}...".format(options.loadparam)
		allparams=np.loadtxt(options.loadparam)
		ttparams, pks, miscglobal=get_params(allparams, options)
	
	else:
		###############
		#### initial alignment for translation and tilt axis angle
		
		#### Translation alignment by ccc
		imgs_trans, pretrans=calc_global_trans(imgs, options)
		
		#### find tilt axis rotation by common line
		if options.tltax:
			tltax=options.tltax
		else:
			tltax=calc_tltax_rot(imgs_trans, options)
			options.tltax=tltax
		
		print "Estimated tilt axis rotation: {:.1f}.".format(tltax)
	
		#### Now initialize parameters to refine
	
		#### 5 params for each tilt: tx, ty, rot, tilt angle, off axis tilt
		ttparams=np.zeros((num, 5))
		ttparams[:,0]=-pretrans[:,0] # tx
		ttparams[:,1]=-pretrans[:,1] # ty
		ttparams[:,2]=tltax # rot
		ttparams[:,3]=tlts.copy() # ytilt
		ttparams[:,4]=0 # off axis tilt

		#### global params
		miscglobal=[0,0] #  tlt axis trans x; tlt axis trans y

		#### particle location
		pks=np.zeros((options.npk, 3))
		
		#### pack parameters together so it is easier to pass around
		allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
		if options.writetmp: 
			np.savetxt(options.tmppath+"params_00.txt", allparams)
	
	
	
	#### Make initial tomogram
	if options.niter>0 or  options.rmgold:
		if options.writetmp:
			outname=options.tmppath+"tomogram_00.hdf"
		else:
			outname=None
		threed=make_tomogram(imgs, allparams, options, premask=True, outname=outname)
	
	err_tilt=[] #### initialize this here to avoid error when niter==0
	
	#### Main refinement loop
	for iti in range(1,options.niter+1):
		ttparams, pks, miscglobal=get_params(allparams, options)
		ttparams_init=ttparams.copy()
		
		#### look for landmarks in the tomogram
		pks=locate_peaks(threed, options)
		
		
		allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
		
		#### make small sample tomograms at the landmarks and center them
		for ii in range(2):
			allparams=make_samples(imgs, allparams, options, refinepos=True)
		
		#### check loss at the begining of the iteration and calculate average loss for each sample and tilt
		ptrg, err_tilt= check_loss(imgs, allparams, options)
		print "Starting iteration {:d}, Mean loss: {:.1f}".format(iti, np.mean(err_tilt))
		
		
		#### write projections of the sample tomograms
		if options.writetmp:
			outname=options.tmppath+"sample_proj_{:02d}.hdf".format(iti)
		else:
			outname=None
		make_samples(imgs, allparams, options, refinepos=False, outname=outname, errtlt=err_tilt)
		
		#### At the first iteration, compute the global shift once
		if iti==1 and options.loadparam:
			allparams, resx, loss = refine_global(imgs, allparams, options, ptrg)
			print "Estimated tilt axis translation: {:.1f}, {:.1f}, loss {:.1f}".format(float(resx[0]), float(resx[1]), float(loss))
		
		#### refine rotation and translation of the tilts
		for rep in range(2):
			allparams = refine_trans(imgs, allparams, options, ptrg, err_tilt)
			allparams = refine_angle(imgs, allparams, options, ptrg, err_tilt)
		
		#### check loss again at the end of the iteration and save parameters
		ptrg, err_tilt0= check_loss(imgs, allparams, options)
		
		#### roll back to coarse alignment if the error is large..
		ttparams, pks, miscglobal=get_params(allparams, options)
		badtlt=err_tilt0>options.minloss*2
		ttparams[badtlt,:]=ttparams_init[badtlt,:].copy()
		imgrot=[]
		for nid in range(num):
			tpm=ttparams[nid]
			pxf=get_xf_pos(ttparams[nid], [0,0,0],  miscglobal)
			po=imgs[nid].copy()
			po.translate(-pxf[0], -pxf[1],0)
			po.rotate(-tpm[2],0,0)
			imgrot.append(po)
		imgs_trans1, pretrans1=calc_global_trans(imgrot, options)
		for nid in range(num):
			if badtlt[nid]:
				tpm=ttparams[nid]
				pt=pretrans1[nid]
				pxf=get_xf_pos(ttparams[nid], [0,0,0],  miscglobal)
				rot=Transform({"type":"2d", "alpha":tpm[2]})
				pr=rot.transform(pt[0], pt[1])
				for ti in [0,1]: ttparams[nid, ti]-=pr[ti]
		
		allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
		ptrg, err_tilt= check_loss(imgs, allparams, options)
		
		print "Iteration {:d} finished. Mean loss: {:.1f}, {:d} bad tilts".format(iti, np.mean(err_tilt), np.sum(err_tilt0>options.minloss*2))
		if options.writetmp: 
			np.savetxt(options.tmppath+"params_{:02d}.txt".format(iti), allparams)
		
		#### make sample tomograms again at the final round
		if iti==options.niter and options.writetmp:
			make_samples(imgs, allparams, options, refinepos=False, outname=options.tmppath+"sample_proj_final.hdf")
		
		#### stop when the minimum loss is reached
		if np.mean(err_tilt)<options.minloss:
			break
		
		#### write tomogram output
		if options.writetmp and iti<options.niter:
			threed=make_tomogram(imgs, allparams, options, outname=options.tmppath+"tomogram_{:02d}.hdf".format(iti), errtlt=err_tilt)
	
	#### final tomogram reconstruction
	imgs=EMData.read_images(inputname)
	
	#### this time skip the lowpass filter
	for m in imgs:
		m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
		m.process_inplace("normalize.edgemean")

	#######
	try: os.mkdir("tomograms")
	except: pass
	
	
	#### remove gold fiducials when needed
	if options.rmgold and (not options.nofiducial):
		imgs=remove_gold(imgs, allparams, threed, options)
		inputtag=inputtag+"_rmgd.hdf"
	
	
	outname="tomograms/" +bname + inputtag +".hdf"
	
	threed=make_tomogram(imgs, allparams, options, premask=True, padr=1.6, clipz=True, errtlt=err_tilt, outname=outname)
	#threed.write_image(outname)
	
	E2end(logid)

#### parse parameters
def get_params(allparams, options):
	num=options.num
	npk=options.npk
	_ttparams=allparams[:num*5].reshape((num,5)).copy()
	_pks=allparams[num*5:num*5+npk*3].reshape((npk, 3)).copy()
	_miscglobal=allparams[num*5+npk*3:].copy()
	if len(_miscglobal)>2:
		print "something wrong..", _miscglobal
	return _ttparams, _pks, _miscglobal

#### get 2D position on a tilt given 3D location
def get_xf_pos(tpm, pk, pre_trans):
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0]+pre_trans[0], pk[1]+pre_trans[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)
	
	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0]-pre_trans[0], p1[1]-pre_trans[1]])
	

	return [p2[0], p2[1]]

#### main function to calculate loss
def calc_loss( inp, imgs, allparams, options, mode="none", iid=[], ptrans=[], mask=False, doplot=False, retall=False, excludes=[0,0]):
	num=len(imgs)
	npk=options.npk
	zeroid=options.zeroid
	bx=options.bxsz/2
	if doplot:
		try: os.remove(options.tmppath+"tmpptcls_calcloss.hdf")
		except: pass
		
		
	ttparams, pks, miscglobal=get_params(allparams, options)
	err=np.zeros((num, npk))
	efill=np.zeros((num, npk))
	
	nrange=np.hstack([np.arange(zeroid, num), np.arange(zeroid, -1, -1)])
	prange=np.arange(npk)
	
	zzid=[zeroid]
	plotn=num
	
	if mode=="zpos":
		prange=iid
		pks[iid,2]=float(inp)
		plotn=0
	elif mode=="global":
		#miscglobal[0]=float(inp)
		miscglobal=inp.copy()
		if len(iid)>0:
			prange=iid
	elif mode=="tlts":
		ttparams[iid[1],2:]=inp.copy()
		if abs(ttparams[iid[1], 3]-options.tlt_init[iid[1]])>options.tltstep*.5:
			return 64**2
		nrange=[iid[0], iid[1]]
		zzid.append(iid[0])
		if len(iid)>2:
			prange=iid[2]
	elif mode=="trans":
		ttparams[iid[1],:2]=inp.copy()
		nrange=[iid[0], iid[1]]
		zzid.append(iid[0])
		if len(iid)>2:
			prange=iid[2]
	
   
	for pid in prange:
		if len(ptrans)==0:
			trans=np.zeros((num,2))
		else:
			trans=ptrans[pid]
		
		
		for ii,nid in enumerate(nrange):
			
			pxf=get_xf_pos(ttparams[nid], pks[pid],  miscglobal)

			pxf[0]+=imgs[nid]["nx"]/2
			pxf[1]+=imgs[nid]["ny"]/2


			e=imgs[nid].get_clip(Region(pxf[0]-bx,pxf[1]-bx, bx*2, bx*2))
			if mask:
				e.process_inplace("normalize.edgemean")
				e.process_inplace("mask.gaussian",{"outer_radius":e["nx"]/4})
				
#			 edge=np.min([pxf[0], pxf[1], imgs[nid]["nx"]-pxf[0], imgs[nid]["ny"]-pxf[1]])

			ts=[0,0]
			if nid not in zzid:

				ea=e.align("translational",e0)
				ts=ea["xform.align2d"].get_trans()
#				 ts[0]+=trans[nrange[ii-1], 0]
#				 ts[1]+=trans[nrange[ii-1], 1]
				trans[nid]=[ts[0], ts[1]]
				err[nid, pid]=np.sqrt(ts[0]**2+ts[1]**2)
				efill[nid,pid]+=1

			e["score"]=[ts[0], ts[1]]
			e["pid"]=pid
		
		
			if doplot:
				e.write_image(options.tmppath+"tmpptcls_calcloss.hdf", nid+pid*plotn)
			
			if nid in zzid:
				e0=e
	
	excn=np.argsort(-np.mean(err, axis=0))[:excludes[0]]
	excp=np.argsort(-np.mean(err, axis=1))[:excludes[1]]
	efill[:, excn]=-1
	efill[excp, :]=-1
	
	if retall:
		return err
	else:
		return np.mean(err[efill>0])
	

#### coarse translational alignment
def calc_global_trans(imgs, options, excludes=[]):
	print "Doing coarse translational alignment..."
	
	num=len(imgs)
	sz=options.minsz
	
	imgout=[0]*num
	e0=imgs[options.zeroid].copy()
	e0.clip_inplace(Region(e0["nx"]/2-sz/2, e0["ny"]/2-sz/2, sz,sz))
	e0.process_inplace("mask.gaussian",{"outer_radius":sz/4})
	e0["xform.align2d"]=Transform()
	imgout[options.zeroid]=e0
	
	if options.writetmp:
		tmpoutname=options.tmppath+"imgs_transali.hdf"
		e0.write_image(tmpoutname,options.zeroid)
	
	
	pretrans=np.zeros((num,2))
	for dr in [-1,1]:
		for i in range(num):
			
			nid=options.zeroid+(i+1)*dr
			if nid>=num or nid<0: continue
			if nid in excludes: continue
			e0=imgout[options.zeroid+i*dr]
			e1=imgs[nid].copy()
			e1.clip_inplace(Region(e1["nx"]/2-sz/2, e1["ny"]/2-sz/2, sz,sz))
			e1.process_inplace("mask.gaussian",{"outer_radius":sz/4})
			
			e1a=e1.align("translational", e0)
			
			e1=imgs[nid].copy()
			e1.transform(e1a["xform.align2d"])
			e1.clip_inplace(Region(e1["nx"]/2-sz/2, e1["ny"]/2-sz/2, sz,sz))
			e1.process_inplace("mask.gaussian",{"outer_radius":sz/4})
			
			imgout[nid]=e1
			if options.writetmp:
				e1.write_image(tmpoutname,nid)
			ts=e1a["xform.align2d"].get_trans()
			if options.verbose: print "\t{:d}: {:.1f}, {:.1f}".format(nid, ts[0], ts[1])
			pretrans[nid, 0]=ts[0]
			pretrans[nid, 1]=ts[1]
	return imgout,pretrans




#### Find tilt axis by common line
def calc_tltax_rot(imgs, options):
	print "Calculateing tilt axis rotation..."
	
	num=len(imgs)
	sz=options.minsz
	
	
	imgnp=[]
	for i in range(num):
		#imgs[i].process_inplace("mask.gaussian",{"outer_radius":sz/4})
		m=get_fft(imgs[i].numpy().copy())
		am=np.abs(m)
		am[am==0]=1
		m/=am
		imgnp.append(m)
	
	
	sm=np.mean(imgnp, axis=0)
	sm=np.abs(sm[:,sz/2:])
	print np.max(sm), np.min(sm)
	rr=np.arange(min(sm.shape[1], sz*.25), dtype=float)
	angs=np.arange(0., 180, .5)
	vs=[]
	for ang in angs:
		a=ang/180.*np.pi
		pts=[np.round(rr*np.sin(a)).astype(int), np.round(rr*np.cos(a)+sz/2).astype(int) ]
		v=sm[pts[1], pts[0]]
		vs.append(np.mean(v))
	vs[0]=vs[180]=0
	tltax=angs[np.argmax(vs)]
	if options.writetmp:
		e=from_numpy(sm)
		e.write_image(options.tmppath+"commonline.hdf")
		np.savetxt(options.tmppath+"tltrot.txt", np.vstack([angs, vs]).T)
	return tltax

#### reconstruct tomogram...
def make_tomogram(imgs, allparams, options, outname=None, premask=True, padr=1.2, clipz=False, errtlt=[]):
	print "Making 3D volume..."
	ttparams, pks, miscglobal=get_params(allparams, options)
	num=len(imgs)
	
	if len(errtlt)==0:
		errtlt=np.zeros(num)
		nrange=range(num)
	else:
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]
		
		
	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	outxy=good_size(max(nx, ny))
	
	pad=good_size(outxy*padr)
	zthick=good_size(pad/2)
	if options.verbose:
		print "\t Image size: {:d} x {:d}".format(nx, ny)
		print "\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick)
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":options.reconmode})
	recon.setup()
	jobs=[]
	
	info=js_open_dict(info_name(options.basename))
	tltinfo=[]
	
	for nid in range(num):
		exclude= nid not in nrange

		tpm=ttparams[nid]

		pxf=get_xf_pos(ttparams[nid], [0,0,0],  miscglobal)

		#pxf[0]+=nx/2
		#pxf[1]+=ny/2
		#rot=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4]})
		xform=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4], "tx":-pxf[0], "ty":-pxf[1]})
		tltinfo.append({"xform.projection":xform, "alignment.score":errtlt[nid]})
		jobs.append([nid,imgs[nid],  recon, pad, xform, premask, exclude, options])
		
	
	thrds=[threading.Thread(target=reconstruct,args=(i)) for i in jobs]
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose : print "Inserting slice {}/{}".format(thrtolaunch,len(thrds))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	for t in thrds: t.join()
	
	threed=recon.finish(True)
	threed.process_inplace("normalize")
	
	if clipz:
		p0=np.min(threed.numpy(), axis=1)
		z0=np.min(p0, axis=1)
		zp=np.where(z0<np.mean(z0))[0]
		zcent=int(zp[0]+zp[-1])/2
		zthk=int((zp[-1]-zp[0])*.6)
		zthk=np.min([zthk, zthick-zcent, zcent])-1
		if options.verbose:
			print "Z axis center at {:d}, thickness {:d} pixels".format(zcent, zthk*2)
		threed.clip_inplace(Region((pad-outxy)/2, (pad-outxy)/2, zcent-zthk, outxy, outxy, zthk*2))
		
		for nid in range(num):
			tltinfo[nid]["xform.projection"].translate(0, 0, zthick/2-zcent)
		
	else:
		
		threed.clip_inplace(Region((pad-outxy)/2, (pad-outxy)/2, 0, outxy, outxy, zthick))
	
	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix
	
	if outname:
		threed.write_image(outname)
		if options.verbose: print "Map written to {}.".format(outname)
		
	info["tiltseries"]=tltinfo
	info.close()	
	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform, premask, exclude, options):
	
	p2=img.get_clip(Region(img["nx"]/2-pad/2,img["ny"]/2-pad/2, pad, pad))
	rr=xform.get_params("xyz")
	if options.writetmp:
		po=p2.copy()
		rr=xform.get_params("xyz")
		po["exclude"]=exclude
		po.translate(rr["tx"], rr["ty"],0)
		po.rotate(-rr["ztilt"],0,0)
		po.process_inplace("normalize")
		po.write_image(options.tmppath+"tmpimg_ali.hdf", nid)
	if premask:
		p2.process_inplace("mask.soft",{"outer_radius":options.minsz/2-1, "width":16, "dx":rr["tx"], "dy":rr["ty"]})
	if not exclude:
		
		p3=recon.preprocess_slice(p2, xform)
		recon.insert_slice(p3,xform,1)

#### locate landmarks in 3D map
def locate_peaks(threed, options, returnall=False):
	print "Looking for landmarks in tomogram..."
	
	mapsmall=threed.process("math.fft.resample",{"n":2})
	mapsmall.mult(-1)
	if options.nofiducial:
		mapsmall.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.1})
	else:
		mapsmall.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.2})
	#mapsmall.process_inplace("filter.highpass.gauss", {"cutoff_abs":.25})
	mapsmall.process_inplace("normalize")
	mapsmall.process_inplace("threshold.belowtozero", {"minval":0})
	mappks=mapsmall.process("mask.onlypeaks")
	mappks.process_inplace("mask.soft",{"outer_radius":options.minsz/4-1})
	#mapsmall.write_image("tmppks.hdf")
	#mappks.write_image("tmppks1.hdf")
	
	zthick=threed["nz"]
	pad=threed["nx"]
	
	print "preprocess done.."
	pkmap=mappks.numpy()
	mn=-np.sort(-pkmap.flatten())[500]
	newpks=np.where(pkmap>mn)
	newpks=(np.array(newpks).T)[:,::-1]
	
	if options.writetmp:
		try: os.remove("tmpptcls_peak3d.hdf")
		except: pass
	scrs=[]
	nnpk=[]
	bx=options.bxsz/4
	for i,p in enumerate(newpks):
		
		
		e=mapsmall.get_clip(Region(p[0]-bx, p[1]-bx, p[2]-bx, bx*2, bx*2, bx*2))
		e=e.project("standard", Transform())
		if options.nofiducial:
			e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			e.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
		else:
			e.process_inplace('normalize')
			e.process_inplace("filter.highpass.gauss",{"cutoff_abs":.25})
		
		e["score"]=e[bx, bx]
		if options.writetmp:
			e.write_image(options.tmppath+"tmpptcls_peak3d.hdf", -1)
		scrs.append(e["score"])
		nnpk.append(p)
		
	nnpk=np.array(nnpk).copy()*2-[pad/2, pad/2, zthick/2]
	
	#### return all peaks for fiducial removal
	if returnall:
		return nnpk, np.array(scrs)
	
	#### pick a few as landmarks for tracking
	scrs=np.array(scrs)
	if options.nofiducial:
		mn=-np.sort(-scrs)[options.npk*4]
	else:
		mn=-np.sort(-scrs)[options.npk*2]
	#print np.mean(scrs)+np.std(scrs), mn
	
	#### try to pick landmarks as far away from each other as possible
	iid=np.sqrt(nnpk[:,0]**2+nnpk[:,1]**2)<options.minsz*.8-bx
	iid=iid* (scrs>mn)
	ps=nnpk[iid]
	ss=np.array(scrs)[iid]
	ppk=[nnpk[np.random.randint(len(nnpk))]] #### just introduce some randomness into the system..
	for i in range(options.npk-1):
		ii=np.argmax(np.min(scidist.cdist(ps, ppk), axis=1))
		if options.verbose:print  ps[ii], ss[ii]
		ppk.append(ps[ii])
	ppk=np.array(ppk) 
	return ppk

#### make sample tomograms and center them
def make_samples(imgs, allparams, options, refinepos=False, outname=None, errtlt=[]):
	if outname:
		try: os.remove(outname)
		except: pass
	
	num=len(imgs)
	npk=options.npk
	ttparams, pks, miscglobal=get_params(allparams, options)
	if len(errtlt)==0:
		nrange=range(num)
	else:
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]
	bx=options.bxsz/2
	for pid in range(npk):
		pad=good_size(bx*4)
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad]})
		recon.setup()

		for nid in nrange:
			
			tpm=ttparams[nid]
			
			pxf=get_xf_pos(ttparams[nid], pks[pid],  miscglobal)

			pxf[0]+=imgs[nid]["nx"]/2
			pxf[1]+=imgs[nid]["ny"]/2

			e=imgs[nid].get_clip(Region(pxf[0]-pad/2,pxf[1]-pad/2, pad, pad))
			p2=e
			rot=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4]})
			p3=recon.preprocess_slice(p2, rot)
			recon.insert_slice(p3,rot,1)
		bxcr=np.round(pks[pid]).astype(int).tolist()
		threed=recon.finish(True)
		threed=threed.get_clip(Region((pad-bx*2)/2,(pad-bx*2)/2,(pad-bx*2)/2,bx*2,bx*2,bx*2))
		threed.process_inplace("normalize")
		
		#### center particles by center of mass in fiducial less mode
		if refinepos and options.nofiducial:
			thd=threed.copy()
			thd.mult(-1)
			thd.process_inplace("mask.soft", {"outer_radius":int(bx*.3)})
			thd.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			com=thd.calc_center_of_mass(0)
			com=[c-bx for c in com]
			
			threed.translate(-com[0], -com[1], -com[2])
			for i in range(3): pks[pid, i]+=com[i]
			if options.verbose:print pid, com
			
		
		pj=threed.project("standard", Transform({"type":"eman", "alt":90}))
		pj["box"]=bxcr
		
		#### center particles by the minima position in projections in fiducial mode
		if refinepos and (not options.nofiducial):
			pj.mult(-1)
			p=pj.calc_max_location()
			zsft=[p[0]-bx, p[1]-bx]
			#pks[pid, 1]+=zsft[0]
			pks[pid, 2]+=zsft[1]
			
			pj.mult(-1)
			pj.translate(-zsft[0],0,0)
		
		if outname: pj.write_image(outname, pid*2+1)
		
		
		pj1=threed.project("standard", Transform())
		pj1["box"]=bxcr
		
		
		if refinepos and (not options.nofiducial):
			pj1.mult(-1)
			p=pj1.calc_max_location()
			xysft=[p[0]-bx, p[1]-bx]
			pks[pid, 0]+=xysft[0]
			pks[pid, 1]+=xysft[1]
			
			pj1.mult(-1)
			pj1.translate(-xysft[0],-xysft[1], 0)
			
			if options.verbose:print pid, zsft,xysft
		
			
		if outname:pj1.write_image(outname, pid*2)
		
	if refinepos:
		allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
	return allparams

#### just check the loss function and return average loss per particle and per tilt
def check_loss(imgs, allparams, options):
	err=calc_loss(0, imgs, allparams, options, mode="none", doplot=1,retall=True, mask=True)
	err_ptcl=np.mean(err, axis=0)
	pkeep=int(options.npk * options.pkeep)
	pkrg=np.argsort(err_ptcl)[:pkeep]

	err_tilt=np.mean(err[:,pkrg], axis=1)
	
	return pkrg, err_tilt

#### calculate global pre-tilt translation
def refine_global(imgs, allparams, options, pkrg=[]):
	print "Refining tilt axis translation.."
	
	if len(pkrg)==0:
		pkrg=range(options.npk)
		
	ttparams, pks,  miscglobal=get_params(allparams, options)    
	res = minimize(calc_loss, miscglobal, (imgs, allparams, options, 'global', pkrg, [], 1), method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
	
	miscglobal=res.x.copy()
	allparams=np.hstack([ttparams.flatten(), pks.flatten(),miscglobal])
	return allparams, res.x, res.fun

#### refine translation of each tilt
def refine_trans(imgs, allparams, options, pkrg=[], err_tilt=[]):
	print "Refining translation for each tilt..."
	
	if len(err_tilt)==0:
		err_tilt=np.zeros(len(imgs))+np.inf
		
	num=len(imgs)
	npk=options.npk
	zeroid=options.zeroid
	ttparams, pks, miscglobal=get_params(allparams, options) 
	
	#### not really used right now. Reserved for sequential alignment
	#nrange=np.hstack([ np.arange(zeroid, -1, -1),np.arange(zeroid, num)])
	#tmptrans=np.zeros((npk, num, 2))  
	
	pl=Pool(options.threads)
	jobs=[]
	for i in range(num):
		jobs.append([i, imgs,allparams, options, zeroid, pkrg, err_tilt]) 
	
	ret=pl.map(do_refine_trans, jobs)
	for r in ret:
		nid, pm, loss = r
		if options.verbose: 
			print "\t tlt: {:d}, trans: {:.1f}, {:.1f}, loss: {:.2f}".format(nid, pm[0],pm[1], loss)
		ttparams[nid,:2]=pm.copy()
	allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
	pl.close()
	pl.join()
	return allparams

#### subprocess function for translation refinement
def do_refine_trans(args):
	nid, imgs,allparams, options, zeroid, pkrg, err_tilt=args
	ttparams, pks, miscglobal=get_params(allparams, options)    
	if (err_tilt[nid]< options.minloss * .8): 
		return [nid, ttparams[nid,:2], err_tilt[nid]]
	
	res = minimize(calc_loss, ttparams[nid,:2], (imgs, allparams, options, 'trans', [zeroid, nid, pkrg], [],1),method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
	return [nid, res.x, res.fun]

#### refine angle of each tilt (including off axis and in plane tilt)
def refine_angle(imgs, allparams, options, pkrg=[], err_tilt=[]):
	print "Refining rotation for each tilt..."
	
	if len(err_tilt)==0:
		err_tilt=np.zeros(len(imgs))+np.inf
		
	num=len(imgs)
	npk=options.npk
	zeroid=options.zeroid
	ttparams, pks, miscglobal=get_params(allparams, options)      
	
	#### not really used right now. Reserved for sequential alignment
	#nrange=np.hstack([ np.arange(zeroid, -1, -1),np.arange(zeroid, num)])
	#tmptrans=np.zeros((npk, num, 2))  
	
	pl=Pool(options.threads)
	jobs=[]
	for i in range(num):
		jobs.append([i, imgs,allparams, options, zeroid, pkrg, err_tilt]) 
	
	ret=pl.map(do_refine_angle, jobs)
	for r in ret:
		nid, pm, loss = r
		if options.verbose: 
			
			print "\t tlt: {:d}, angle: {:.1f}, {:.1f}, {:.1f}, loss: {:.2f}".format(nid, float(pm[0]), float(pm[1]), float(pm[2]), float(loss))
		ttparams[nid,2:]=pm.copy()
		
	allparams=np.hstack([ttparams.flatten(), pks.flatten(), miscglobal])
	pl.close()
	pl.join()
	return allparams

#### subprocess function for rotation refinement
def do_refine_angle(args):
	nid, imgs,allparams, options, zeroid, pkrg, err_tilt=args
	ttparams, pks, miscglobal=get_params(allparams, options)    
	if (err_tilt[nid]< options.minloss * .8): 
		return [nid, ttparams[nid,2:], err_tilt[nid]]
	
	ttparams, pks, miscglobal=get_params(allparams, options)    
	res = minimize(calc_loss, ttparams[nid,2:], (imgs, allparams, options, 'tlts', [zeroid, nid, pkrg], [],1),method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
	return [nid, res.x, res.fun]

#### remove gold fiducials
def remove_gold(imgs, allparams, threed, options):
	print "Starting removing gold..."
	goldsize=12.
	
	ttparams, pks, miscglobal=get_params(allparams, options)    
	num=len(imgs)
	npk=options.npk
	
	pks, scrs=locate_peaks(threed, options, returnall=True)
	thr=np.mean(scrs)+np.std(scrs)
	pks=pks[scrs>thr]
	print "{:d} gold fiducials found.".format(len(pks))
	
	imgs_rmgd=[m.copy() for m in imgs]
	b2=options.bxsz/4
	
	jobs=[]
	for i in range(num):
		jobs.append([i, imgs_rmgd[i], pks, ttparams,miscglobal, b2 , goldsize]) 
	
	thrds=[threading.Thread(target=do_rm_gold,args=([i])) for i in jobs]
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose : print "Starting on tilt {}/{}".format(thrtolaunch,len(thrds))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	for t in thrds: t.join()
	
	if options.writetmp:
		for i,m in enumerate(imgs_rmgd):
			m.write_image(options.tmppath+"tmpimg_rmgd.hdf", i)
		
	#threed=make_tomogram(imgs_rmgd, allparams, options)
	if options.verbose: print "Gold removal finished."
	return imgs_rmgd

#### subprocess for gold removal
def do_rm_gold(args):
	
	nid, img, pks, ttparams,miscglobal, b2 , goldsize=args
	
	for pp in pks:
		pxf=get_xf_pos(ttparams[nid],pp, miscglobal)

		pxf[0]+=img["nx"]/2
		pxf[1]+=img["ny"]/2
		
		e=img.get_clip(Region(pxf[0]-b2,pxf[1]-b2, b2*4, b2*4))
		e.process_inplace("mask.soft",{"inner_radius":goldsize})
		meanval=e["mean_nonzero"]
		img.process_inplace("mask.soft",{"dx":int(pxf[0]-img["nx"]/2),"dy":int(pxf[1]-img["ny"]/2),"inner_radius":goldsize, "value":meanval})
	
	return

def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	