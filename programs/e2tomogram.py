#!/usr/bin/env python
from __future__ import print_function
# Muyuan Chen 2017-04
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize
import scipy.spatial.distance as scidist
from EMAN2_utils import *
from multiprocessing import Pool
from sklearn.decomposition import PCA

def main():
	
	usage="""WARNING: This is an experimental program that may or may not work..
	This program takes an unaligned tilt series, performs alignment, and generate a tomogram.
	
	Usage:
	e2tomogram.py <tilt series stack> --rawtlt <raw tilt file> [options]
	
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="tiltseries",help="Specify the tilt series you intend to reconstruct.", default="", guitype='filebox', browser="EMTiltseriesTable(withmodal=True,multiselect=True)", filecheck=False, row=0, col=0,rowspan=1, colspan=2,nosharedb=True,mode="easy")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Optional:", row=1, col=0, rowspan=1, colspan=2,mode="easy")

	parser.add_argument("--rawtlt", type=str,help="Override tilt angles stored in project metadata and tiltseries header. Text file contains raw tilt angles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False, row=3, col=0, rowspan=1, colspan=2)#,mode="easy")

	parser.add_argument("--zeroid", type=int,help="Index of the center tilt. Ignored when rawtlt is provided.", default=-1,guitype='intbox',row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--tltstep", type=float,help="Step between tilts. Ignored when rawtlt is provided.", default=1,guitype='floatbox',row=4, col=1, rowspan=1, colspan=1)

	parser.add_argument("--loadparam", type=str,help="Load from existing param file", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False, row=5, col=0, rowspan=1, colspan=2)

	parser.add_argument("--tltax", type=float,help="Angle of the tilt axis. The program will calculate one if this option is not provided", default=None,guitype='floatbox',row=6, col=0, rowspan=1, colspan=1)
	parser.add_argument("--tltkeep", type=float,help="Fraction of tilts to keep in the reconstruction.", default=.9,guitype='floatbox',row=6, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--npk", type=int,help="Number of landmarks to use.", default=20,guitype='intbox',row=7, col=0, rowspan=1, colspan=1)
	parser.add_argument("--pkeep", type=float,help="Fraction of landmarks to keep in the tracking.", default=.5,guitype='floatbox',row=7, col=1, rowspan=1, colspan=1)

	parser.add_argument("--bxsz", type=int,help="Box size of the particles for tracking", default=-1,guitype='intbox',row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--minloss", type=float,help="Stop refinement when the loss is lower than this value.", default=1.,guitype='floatbox',row=8, col=1, rowspan=1, colspan=1)

	parser.add_argument("--fidsz", type=float,help="Radius of gold fiducials in nm.", default=5.)
	parser.add_argument("--rmgold", action="store_true",help="Remove gold fiducials.", default=False,guitype='boolbox',row=9, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--nofiducial", action="store_true",help="Fiducial-less mode. This will change a few internal parameters to make it work.", default=False, guitype='boolbox',row=9, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--niter", type=str,help="Number of iterations for bin8, bin4, bin2 images. Default if 2,1,1", default="2,1,1",guitype='strbox',row=10, col=0, rowspan=1, colspan=1,mode='easy[2,1,1]')

	parser.add_argument("--writetmp", action="store_true",help="Write intermidiate files", default=False,guitype='boolbox',row=10, col=1, rowspan=1, colspan=1,mode="easy[True]")

	parser.add_argument("--shrink", type=float,help="Mean-shrink tilt images by this integer value. We suggest using a value that converts image dimensions to approximately 1024x1024. The default is 4.", default=4, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode="easy[4]")
	parser.add_argument("--reconmode", type=str,help="Reconstruction mode. Choose from nearest_neighbor, gauss_2, gauss_3, and gauss_5.", default="gauss_2", choices=["gauss_2","gauss_3","gauss_5"], choicelist='["gauss_2","gauss_3","gauss_5"]', guitype='combobox', row=11, col=1, rowspan=1, colspan=1)

	parser.add_argument("--tmppath", type=str,help="Temporary path", default=None)

	parser.add_argument("--threads", type=int,help="Number of threads", default=12,guitype='intbox',row=12, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0, mode="easy[3]", guitype="intbox", row=12,col=1, rowspan=1, colspan=1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)

	inputname=args[0]
	time0=time.time()

	itnum=[int(i) for i in options.niter.split(',')]
	print(itnum)


	options.inputname=inputname

	dotpos=inputname.rfind('.')
	linepos=inputname.rfind('__')
	inputtag=inputname[linepos:dotpos]
	bname=base_name(inputname)

	options.basename=bname




	#### make a folder to write tmp files
	if options.tmppath:
		path=options.tmppath
	else:
		for i in range(100):
			try:
				path="tomorecon_{:02d}/".format(i)
				os.mkdir(path)
				options.tmppath=path
				break
			except:
				continue
		else:
			print("Too many tomorecon folders in the project, or something odd happened....Exit.")
			exit()

	print("Temporary files will be written in {}".format(options.tmppath))

	e=EMData(inputname, 0, True)
	inppath=options.tmppath+"tltseries_input.hdf"
	cmd="e2proc2d.py {} {} --inplace ".format(inputname, inppath)
	binfac=max(1, int(np.round(e["nx"]/2048.)))
	cmd+=" --process threshold.clampminmax.nsigma:nsigma=10"
	if binfac>1:
		cmd+=" --meanshrink {}".format(binfac)
	options.binfac=binfac

	if e["nz"]>1:
		cmd+=" --threed2twod"

	run(cmd)

	## now prepare tilt series
	imgs_2k=EMData.read_images(inppath)
	imgs_1k=[img.process("math.meanshrink", {"n":2}).process("normalize") for img in imgs_2k]
	imgs_500=[]
	for p in imgs_1k:
		m=p.process("math.minshrink", {"n":options.shrink})
		m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
		m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
		m.process_inplace("normalize.edgemean")
		m["apix_x"]=m["apix_y"]=p["apix_x"]*2.
		imgs_500.append(m)

	options.apix_init=float(imgs_2k[0]["apix_x"])

	if options.rawtlt == "":
		jsd = js_open_dict(info_name(inputname))
		try:
			tlts = np.asarray(jsd["tilt_angles"])
		except:
			print(""""ERROR: Could not locate tilt angles. Either import using raw data 
			import tools or specify a rawtlt file containing one tilt angle per line
			(corresponding to the tilt series to be reconstructed).""")
			sys.exit(1)
	else:
		tlts=np.loadtxt(options.rawtlt)
	np.savetxt(options.tmppath+"rawtilt.txt", tlts)
	zeroid=options.zeroid=np.argmin(abs(tlts))


	#### here we always assume the center tilt is at 0 degree
	tlts-=tlts[options.zeroid]

	#### course alignment
	img_tali, pretrans=calc_global_trans(imgs_500, options)
	for i,m in enumerate(img_tali):
		m.write_image(options.tmppath+"tltseries_transali.hdf", i)

	#### initial tilt axis
	tltax=calc_tltax_rot(img_tali, options)
	print("tilt axis:  {:.2f}".format(tltax))
	options.tltax=tltax

	pretrans*=4
	num=options.num=len(imgs_500)
	ttparams=np.zeros((num, 5))
	ttparams[:,0]=-pretrans[:,0] # tx
	ttparams[:,1]=-pretrans[:,1] # ty
	ttparams[:,2]=tltax # rot
	ttparams[:,3]=tlts.copy() # ytilt
	ttparams[:,4]=0 # off axis tilt
	pks=np.zeros((options.npk, 3))
	#### pack parameters together so it is easier to pass around
	allparams=np.hstack([ttparams.flatten(), pks.flatten()])

	#### some fixed parameters..
	options.fid_mindist=16
	options.fid_maxval=-10
	options.bxsz=32



	#### image scale, m3diter, fidkeep
	#### refinement sequence []:global tilt axis, 0:tx, 1:ty, 2:tilt axis, 3: tilt, 4: off axis tilt
	scaleiter=[(imgs_500, itnum[0], .6, [[0,1], [], [0,1],[], [0,1]]),
			   (imgs_1k, itnum[1], .7, [[0,1],[], [0,1],[3],[4],[2],[0,1]]),
			   (imgs_2k, itnum[2], .8, [[0,1],[], [0,1],[3],[4],[2],[0,1]])]

	loss0=[]
	js=js_open_dict(path+"0_tomorecon_params.json")
	js.update(vars(options))
	js.close()
	yrot=0
	for niter, siter in enumerate(scaleiter):
		imgs_, n_m3d, options.fidkeep, rfseq = siter
		apix=float(imgs_[0]["apix_x"])
		binx=apix/options.apix_init*options.binfac
		print("\n******************************")
		print("Iteration {}. Refining alignment on bin{:.0f} images...".format(niter, binx))
		print("Image size {} x {}, Apix {:.2f}".format(imgs_[0]["nx"], imgs_[0]["ny"], apix))
		if niter<2:
			print("Low resolution mode. Using peak positon of landmarks for alignment")
		else:
			print("High resolution mode. Using center of mass of landmarks for alignment")

		#### make tomogram loop
		for m3diter in range(n_m3d):

			#### make tomogram, always use 500x500
			threed=make_tomogram(imgs_500, ttparams, options, premask=True, errtlt=loss0)

			if niter==0:
				rot=fix_rotation(threed)
				yrot+=rot
				ttparams[:,3]+=rot
				zeroid=options.zeroid=np.argmin(abs(ttparams[:,3]))
				ttparams[:,3]-=ttparams[options.zeroid,3]

			threed=make_tomogram(imgs_500, ttparams, options, premask=True, outname=path+"tomo_{:02d}.hdf".format(niter), errtlt=loss0)

			pks=find_landmark(threed, options)
			allparams=np.hstack([ttparams.flatten(), pks.flatten()])

			allparams=make_samples(imgs_, allparams, options, refinepos=True);
			allparams=make_samples(imgs_, allparams, options, refinepos=True);
			#if m3diter==0:
			if niter==0 and m3diter==0:
				make_samples(imgs_, allparams, options, outname=path+"samples_init.hdf", refinepos=True);
			#make_samples(imgs_, allparams, options, outname=path+"samples_{:02d}_{:02d}_init.hdf".format(niter, m3diter), refinepos=True);

			#ptclpos=ali_ptcls(imgs_, allparams, options,path+"ptclali_{:02d}.hdf".format(niter), True)
			#rfseq=[[0,1],[],[0,1]]

			for idx in rfseq:
				allparams=refine_one_iter(imgs_, allparams, options, idx)

			ttparams, pks=get_params(allparams, options)

		make_samples(imgs_, allparams, options, outname=path+"samples_{:02d}.hdf".format(niter), refinepos=False);
		ttparams, pks=get_params(allparams, options)
		make_ali(imgs_2k, ttparams, options, outname="ali_{:02d}.hdf".format(niter))

		ptclpos=ali_ptcls(imgs_, allparams, options,path+"ptclali_{:02d}.hdf".format(niter), True)
		loss0=np.zeros(num)
		tpm=ttparams.copy()
		for nid in range(num):
			loss0[nid]=get_loss_pm([0], nid, allparams, options, [0], ptclpos)

		print("Iteration {} finished. Final average loss {:.2f} nm".format(
			niter, np.mean(loss0)))

		#### always save parameters at the full scale (4k)
		tpm[:,:2]*=options.binfac
		tpm=np.hstack([np.arange(len(tpm))[:,None], tpm])
		np.savetxt(path+"landmarks_{:02d}.txt".format(niter), pks*options.binfac, fmt="%.1f")
		np.savetxt(path+"tltparams_{:02d}.txt".format(niter), tpm, fmt="%.2f")
		np.savetxt(path+"loss_{:02d}.txt".format(niter), np.vstack([np.arange(len(loss0)), loss0]).T, fmt="%.2f")


	threed=make_tomogram(imgs_1k, ttparams, options, premask=False, outname=path+"tomo_final.hdf", clipz=True, errtlt=loss0)
	make_ali(imgs_full, ttparams, options)

	try: os.mkdir("tomograms")
	except: pass
	sfx=""
	if options.binfac>1:
		sfx+="__bin{:d}".format(int(options.binfac*2))

	threed["ytilt"]=yrot
	tomoname=os.path.join("tomograms", options.basename+sfx+".hdf")
	threed.write_image(tomoname)
	print("Tomogram written to {}".format(tomoname))
	js=js_open_dict(info_name(tomoname))
	js["tlt_params"]=ttparams.tolist()
	js["tlt_file"]=options.inputname
	js.close()

	js=js_open_dict(path+"0_tomorecon_params.json")
	js.update(vars(options))
	js.close()

	dtime=time.time()-time0
	print("Finished. Total time: {:.1f}s".format(dtime))

	E2end(logid)

#### rotate tomogram
def fix_rotation(threed):

	threed.process_inplace("math.minshrink", {"n":2})
	threed.process_inplace("filter.highpass.gauss",{"cutoff_pixels":5})
	thd=threed.numpy().copy()
	prj=np.min(thd, axis=1)
	pval=prj.flatten()
	thr=np.sort(pval)[int(len(pval)*.1)]
	pts=np.array(np.where(prj<thr)).T
	pca=PCA(1)
	pca.fit(pts);
	c=pca.components_[0]
	rot=(np.arctan2(c[0], c[1])*180/np.pi)
	rot=(rot+90)%180-90.
	print("Adjusting rotation {:.1f}..".format(rot))
	return rot

#### parse parameters
def get_params(allparams, options):
	num=options.num
	npk=options.npk
	if len(allparams)!=num*5+npk*3:
		print("parameter number does not match...")
	_ttparams=allparams[:num*5].reshape((num,5)).copy()
	_pks=allparams[num*5:num*5+npk*3].reshape((npk, 3)).copy()

	return _ttparams, _pks

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


#### coarse translational alignment
def calc_global_trans(imgs, options, excludes=[]):
	print("Doing coarse translational alignment...")

	num=len(imgs)
	sz=min(imgs[0]["nx"], imgs[0]["ny"])

	imgout=[0]*num
	e0=imgs[options.zeroid].copy()
	e0.clip_inplace(Region(e0["nx"]/2-sz/2, e0["ny"]/2-sz/2, sz,sz))
	e0.process_inplace("mask.gaussian",{"outer_radius":sz/4})
	e0["xform.align2d"]=Transform()
	imgout[options.zeroid]=e0

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
			ts=e1a["xform.align2d"].get_trans()
			if options.verbose: print("\t{:d}: {:.1f}, {:.1f}".format(nid, ts[0], ts[1]))
			pretrans[nid, 0]=ts[0]
			pretrans[nid, 1]=ts[1]
	return imgout,pretrans




#### Find tilt axis by common line
def calc_tltax_rot(imgs, options):
	print("Calculateing tilt axis rotation...")

	num=len(imgs)
	sz=min(imgs[0]["nx"], imgs[0]["ny"])


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
	print(np.max(sm), np.min(sm))
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
	e=from_numpy(sm)
	e.write_image(options.tmppath+"commonline.hdf")
	np.savetxt(options.tmppath+"tltrot.txt", np.vstack([angs, vs]).T)
	return tltax

#### reconstruct tomogram...
def make_tomogram(imgs, tltpm, options, outname=None, premask=True, padr=1.2, clipz=False, errtlt=[]):

	num=len(imgs)
	scale=imgs[0]["apix_x"]/options.apix_init
	print("Making bin{:d} tomogram...".format(int(options.binfac*np.round(scale))))
	ttparams=tltpm.copy()
	ttparams[:,:2]/=scale

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
		print("\t Image size: {:d} x {:d}".format(nx, ny))
		print("\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick))
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":"gauss_2"})
	#recon=Reconstructors.get("fourier_iter", {"size":[pad,pad,zthick]})
	recon.setup()
	jobs=[]

	#info=js_open_dict(info_name(options.basename))
	#tltinfo=[]

	for nid in range(num):
		exclude= nid not in nrange

		tpm=ttparams[nid]

		pxf=get_xf_pos(ttparams[nid], [0,0,0])

		xform=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4], "tx":-pxf[0], "ty":-pxf[1]})
		#tltinfo.append({"xform.projection":xform, "alignment.score":errtlt[nid]})
		jobs.append([nid,imgs[nid],  recon, pad, xform, premask, exclude, options])


	thrds=[threading.Thread(target=reconstruct,args=(i)) for i in jobs]
	thrtolaunch=0
	tsleep=threading.active_count()
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose : print("Inserting slice {}/{}".format(thrtolaunch,len(thrds)))
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
		#if options.verbose:
		print("Z axis center at {:d}, thickness {:d} pixels".format(zcent, zthk*2))
		threed.clip_inplace(Region((pad-outxy)/2, (pad-outxy)/2, zcent-zthk, outxy, outxy, zthk*2))
		threed["zshift"]=float(zthick/2-zcent)*scale*options.binfac
	#for nid in range(num):
	#tltinfo[nid]["xform.projection"].translate(0, 0, zthick/2-zcent)

	else:

		threed.clip_inplace(Region((pad-outxy)/2, (pad-outxy)/2, 0, outxy, outxy, zthick))
		threed["zshift"]=0

	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix

	if outname:
		threed.write_image(outname)
		if options.verbose: print("Map written to {}.".format(outname))

	#info["tiltseries"]=tltinfo
	#info.close()
	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform, premask, exclude, options):

	maxsz=max(img["nx"], img["ny"])
	p2=img.get_clip(Region(img["nx"]/2-pad/2,img["ny"]/2-pad/2, pad, pad))
	rr=xform.get_params("xyz")
	if premask:
		p2.process_inplace("mask.soft",{"outer_radius":maxsz/2-1, "width":16, "dx":rr["tx"], "dy":rr["ty"]})
	#p2.write_image(options.tmppath+"tmpimg.hdf",nid)
	if not exclude:

		p3=recon.preprocess_slice(p2, xform)
		recon.insert_slice(p3,xform,1)


def make_ali(imgs, tpm, options, outname="tltseries_ali.hdf"):

	scale=imgs[0]["apix_x"]/options.apix_init
	ttparams=tpm.copy()
	ttparams[:,:2]/=scale

	fname=options.tmppath+outname
	try:os.remove(fname)
	except:pass
	pad=imgs[0]["nx"]*1.
	mskrd=min(imgs[0]["nx"],imgs[0]["ny"])/2
	for nid, im in enumerate(imgs):
		tpm=ttparams[nid]

		pxf=get_xf_pos(ttparams[nid], [0,0,0])
		m=im.process("normalize")

		p2=m.get_clip(Region(m["nx"]/2-pad/2,m["ny"]/2-pad/2, pad, pad), fill=0)
		po=p2.copy()
		po.translate(-pxf[0], -pxf[1], 0)
		po.rotate(-tpm[2],0,0)
		xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4]})
		po["xform.projection"]=xform
		po.write_image(fname, nid)

def find_landmark(threed, options):
	print("Searching for landmarks...")
	threedtiny=threed.process("math.minshrink", {"n":2})
	threedtiny.process_inplace("normalize")
	threedtiny.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
	mapnp=threedtiny.numpy().copy()
	asrt= np.argsort(mapnp.flatten())

	pts=[]
	dthr=options.fid_mindist
	vthr=options.fid_maxval
	for i in range(len(asrt)):
		aid=asrt[i]
		pt=np.unravel_index(aid, mapnp.shape)
		#	 mapnp[pt]=np.inf
		if len(pts)>0:
			dst=scidist.cdist(pts, [pt])
			if np.min(dst)<dthr:
				continue

		pts.append(pt)
		if mapnp[pt]>vthr:
			break
		#print(mapnp[pt])
		if len(pts)>=options.npk*1.2:
			break


	if len(pts)<options.npk:
		print("Found only {} landmarks".format(len(pts)))
		options.npk=len(pts)
	else:
		#np.random.shuffle(pts)
		pts=pts[:options.npk]

	pks=(np.array(pts)-np.array(mapnp.shape)/2.)
	pks=pks[:,::-1]
	#### mult 2 since we bin 2 in the begining.
	scale=float(threed["apix_x"])/options.apix_init*2.
	pks*=scale

	return pks

#### make sample tomograms and center them
def make_samples(imgs, allparams, options, refinepos=False, outname=None, errtlt=[]):

	if outname:
		try: os.remove(outname)
		except: pass

	num=len(imgs)
	npk=options.npk
	ttparams, pks=get_params(allparams, options)
	scale=float(imgs[0]["apix_x"])/options.apix_init
	ttparams[:,:2]/=scale
	pks/=scale
	lowres=(scale>1.5)
	if len(errtlt)==0:
		nrange=range(num)
	else:
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]
	bx=options.bxsz/2
	if not lowres:
		bx=int(bx*1.5)

	for pid in range(npk):
		pad=good_size(bx*4)
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad]})
		#recon=Reconstructors.get("fourier_iter", {"size":[pad,pad,pad]})
		recon.setup()

		for nid in nrange:

			tpm=ttparams[nid]

			pxf=get_xf_pos(ttparams[nid], pks[pid])

			pxf[0]+=imgs[nid]["nx"]/2
			pxf[1]+=imgs[nid]["ny"]/2

			e=imgs[nid].get_clip(Region(pxf[0]-pad/2,pxf[1]-pad/2, pad, pad))
			e.process_inplace("normalize")
			p2=e
			rot=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4]})
			p3=recon.preprocess_slice(p2, rot)
			recon.insert_slice(p3,rot,1)
		bxcr=np.round(pks[pid]).astype(int).tolist()
		threed=recon.finish(True)
		threed=threed.get_clip(Region((pad-bx*2)/2,(pad-bx*2)/2,(pad-bx*2)/2,bx*2,bx*2,bx*2))
		threed.process_inplace("normalize")
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=imgs[0]["apix_x"]
		threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		threed.process_inplace("mask.soft",{"outer_radius":-1})
		pj=threed.project("standard", Transform({"type":"eman", "alt":90}))
		pj["box"]=bxcr

		#### center particles by the minima position in projections in fiducial mode
		if refinepos:
			zsft=get_center(pj, lowres)

			pks[pid, 2]-=zsft[1]

			pj.translate(zsft[0],zsft[1],0)

		if outname: pj.write_image(outname, pid*2+1)


		pj1=threed.project("standard", Transform())
		pj1["box"]=bxcr


		if refinepos:
			xysft=get_center(pj1, lowres)
			#xysft=[p[0]-bx, p[1]-bx]
			pks[pid, 0]-=(xysft[0]+zsft[0])/2.
			pks[pid, 1]-=xysft[1]

			#			 pj1.mult(-1)
			pj1.translate(xysft[0],xysft[1], 0)

			if options.verbose:print(pid, zsft,xysft)


		if outname:pj1.write_image(outname, pid*2)

	if refinepos:
		ttparams[:,:2]*=scale
		pks*=scale
		allparams=np.hstack([ttparams.flatten(), pks.flatten()])

	if outname: print("Landmark samples written to {}.".format(outname))
	return allparams

def get_center(img, lowres=True):
	e=img.copy()
	bx=e["nx"]/2
	e.mult(-1)
	#### low resolution mode: use peak
	if lowres:
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.01})
		e.process_inplace("mask.gaussian",{"outer_radius":bx*.5})
		pk=e.calc_max_location()
		if e["sigma"]==0:
			pk=[bx,bx]

	#### high resolution mode: use center of mass
	else:
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.01})
		e.process_inplace("mask.gaussian",{"outer_radius":bx*.5})
		e.process_inplace("normalize")
		pk=e.calc_center_of_mass(2)
		if np.isnan(pk).any():
			pk=[bx,bx]


	tx=[bx-pk[0], bx-pk[1]]

	return tx


# nrange=np.arange(options.num)
def ali_ptcls(imgs, allpms, options, outname=None, doali=True):
	zeroid=options.zeroid
	num=options.num
	nrange=np.hstack([np.arange(zeroid, num), np.arange(zeroid, -1, -1)])
	ttparams, pks=get_params(allpms, options)
	scale=float(imgs[0]["apix_x"])/options.apix_init
	ttparams[:,:2]/=scale
	pks/=scale
	prange=np.arange(options.npk)

	k=0
	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	fidptcls=[]
	bx=options.bxsz/2
	apix=imgs[0]["apix_x"]
	lowres=(scale>1.5)
	if not lowres:
		bx=int(bx*1.5)

	ptclpos=[]
	for pid in prange:
		trans=np.zeros((num,2))
		fid=[]
		ppos=np.zeros((num,2))
		for ii,nid in enumerate(nrange):

			pxf=get_xf_pos(ttparams[nid], pks[pid])

			pxf[0]+=nx/2
			pxf[1]+=ny/2

			if nid!=zeroid:
				tlast=trans[nrange[ii-1]]
				pxf[0]-=tlast[0]
				pxf[1]-=tlast[1]
			else:
				tlast=np.array([0,0])

			e=imgs[nid].get_clip(Region(pxf[0]-bx,pxf[1]-bx, bx*2, bx*2)).process("normalize")

			ts=[0,0]
			if doali:# and nid!=zeroid:

				tx=get_center(e, lowres)
				trans[nid]=tlast+np.array([tx[0], tx[1]])

				if nid!=zeroid:

					if outname:
						e=imgs[nid].get_clip(Region(pxf[0]-bx-tx[0],pxf[1]-bx-tx[1], bx*2, bx*2)).process("normalize")

			if outname:
				e["score"]=trans[nid].tolist()
				e["pid"]=pid
				e["nid"]=nid
				e.write_image(outname, nid+pid*num)
			ppos[nid]=np.array(pxf-trans[nid]-[nx/2, ny/2])

		ptclpos.append(ppos)

	ptclpos=np.array(ptclpos)*scale
	#	 print ptclpos.shape, lowres, gdsz
	return ptclpos

def refine_one_iter(imgs, allparams, options, idx=[]):

	apms=make_samples(imgs, allparams, options, refinepos=True);
	ttparams, pks=get_params(apms, options)
	ptclpos=ali_ptcls(imgs, apms, options, doali=True)

	pmlabel=["trans_x", "trans_y", "tlt_z", "tilt_y", "tilt_x"]
	if len(idx)==0:
		res=minimize(global_rot, 0, (ptclpos, apms, options), method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
		print("refine global tilt_z {:.2f}, loss {:.2f} -> {:.2f}".format(
			float(res.x)*.5, float(global_rot(0,ptclpos,apms, options)),
			float(global_rot(res.x,ptclpos,apms, options))))
		ttparams[:,2]-=res.x*.5
		pkc=pks.copy()
		r=res.x*np.pi/180.
		rotmat=np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
		pkc[:,:2]=np.dot(pkc[:,:2],rotmat)
		pks=pkc

	else:

		ttpm_new=ttparams.copy()
		num=len(ttparams)
		loss=np.zeros(num)
		loss0=np.zeros(num)
		for nid in range(num):
			tminit=np.zeros(len(idx))
			res=minimize(get_loss_pm, tminit,(nid,apms, options, idx, ptclpos) , method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
			ttpm_new[nid, idx]+=res.x*.5
			loss[nid]=get_loss_pm(res.x*.5, nid, apms, options, idx, ptclpos)
			loss0[nid]=get_loss_pm(tminit, nid, apms, options, idx, ptclpos)

		ipm=[pmlabel[i] for i in idx]
		print("refine {}, loss: {:.2f} -> {:.2f}".format(ipm, np.mean(loss0), np.mean(loss)))
		ttparams=ttpm_new

	allparams=np.hstack([ttparams.flatten(), pks.flatten()])
	return allparams

def get_loss_pm(pm, nid, apms, options, idx=[2,3,4], ptclpos=[]):
	ttparams, pks=get_params(apms, options)
	tpm=ttparams[nid].copy()
	tpm[idx]+=pm
	ps=np.array([get_xf_pos(tpm, p) for p in pks])
	dst=np.sqrt(np.sum((ps-ptclpos[:,nid])**2, axis=1))
	dst=np.mean(dst[np.argsort(dst)[:int(len(dst)*options.fidkeep)]])

	#### return distance in nm
	return dst*options.apix_init/10.

def global_rot(rt, ptclpos, allpms, options):
	try:
		rt=rt[0]
	except: pass
	errs=[]
	ttparams, pks=get_params(allpms, options)
	for nid in range(options.num):
		p0=ptclpos[:, nid].copy()
		tpm=ttparams[nid].copy()
		pkc=pks.copy()
		r=rt*np.pi/180.
		rotmat=np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
		pkc[:,:2]=np.dot(pkc[:,:2],rotmat)
		p1=np.array([get_xf_pos(tpm, p) for p in pkc])
		t=-r
		rotmat1=np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

		p1=np.dot(p1, rotmat1)
		err=np.sqrt(np.sum((p0-p1)**2, axis=1))
		errs.append(np.mean(np.sort(err)[:int(len(err)*options.fidkeep)]))

	#### return distance in nm
	return np.mean(errs)*options.apix_init/10.

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)


if __name__ == '__main__':
	main()
