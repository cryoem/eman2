#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# Muyuan Chen 2017-04
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize
import scipy.spatial.distance as scidist
import scipy.ndimage as sciimg
from EMAN2_utils import *
import queue
from sklearn.decomposition import PCA

def main():
	
	usage="""
	This program takes unaligned tilt series, performs alignment, and generate a tomogram.
	
	Usage:
		e2tomogram.py <tilt series stack> --rawtlt <raw tilt file> [options]
	or:
		e2tomogram.py <tilt series stack> --tltstep <angle between tilts> [options]
	
	Note: Tiltseries must have the correct Apix values in their headers.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="tiltseries",help="Specify the tilt series you intend to reconstruct.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='tiltseries',multiselect=True)", filecheck=False, row=0, col=0,rowspan=1, colspan=2,nosharedb=True,mode="easy")
	
	parser.add_argument("--alltiltseries", action="store_true",help="Use all tilt series in the folder. Acceptable file extensions include hdf, mrc, mrcs, st.", default=False,guitype='boolbox',row=1, col=0, rowspan=1, colspan=1,mode="easy")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Specify either zeroid/tltstep OR rawtlt:", row=2, col=0, rowspan=1, colspan=2,mode="easy")

	parser.add_argument("--zeroid", type=int,help="Index of the center tilt. Ignored when rawtlt is provided.", default=-1)#,guitype='intbox',row=3, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tltstep", type=float,help="Step between tilts. Ignored when rawtlt is provided. Default is 2.0.", default=2.0,guitype='floatbox',row=3, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--rawtlt", type=str,help="Specify a text file contains raw tilt angles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False, row=4, col=0, rowspan=1, colspan=2)#,mode="easy")

	parser.add_argument("--npk", type=int,help="Number of landmarks to use (such as gold fiducials). Default is 20.", default=20,guitype='intbox',row=5, col=0, rowspan=1, colspan=1, mode="easy")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=6, col=0, rowspan=1, colspan=2,mode="easy")

	parser.add_argument("--tltax", type=float,help="Angle of the tilt axis. The program will calculate one if this option is not provided", default=None,guitype='floatbox',row=7, col=1, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--tltkeep", type=float,help="Fraction of tilts to keep in the reconstruction.", default=.9,guitype='floatbox',row=7, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tltrange", type=str,help="Include only tilts between 'START' and 'STOP', i.e. -40.0,40.0. Default behavior is to include all tilts.", default=None)#, guitype='strbox',row=6, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--outsize", type=str,help="Size of output tomograms. choose from 1k, 2k and 4k. default is 1k", default="1k",guitype='combobox',choicelist="('1k', '2k')",row=8, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--niter", type=str,help="Number of iterations for bin8, bin4, bin2 images. Default if 2,1,1,1", default="2,1,1,1",guitype='strbox',row=8, col=1, rowspan=1, colspan=1,mode='easy[2,1,1,1]')
	
	parser.add_argument("--badzero", action="store_true",help="In case the 0 degree tilt is bad for some reason...", default=False)#, guitype='boolbox',row=9, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--bytile", action="store_true",help="make final tomogram by tiles.. ", default=False, guitype='boolbox',row=9, col=1, rowspan=1, colspan=1,mode="easy[True]")
	
	parser.add_argument("--load", action="store_true",help="load existing tilt parameters.", default=False,guitype='boolbox',row=10, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--notmp", action="store_true",help="Do not write temporary files.", default=False, mode="easy[True]", guitype="boolbox", row=10,col=1, rowspan=1, colspan=1)

	parser.add_argument("--pkkeep", type=float,help="Fraction of landmarks to keep in the tracking.", default=.9,guitype='floatbox',row=11, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--clipz", type=int,help="Z thickness of the final tomogram output. default is -1, (5/16 of tomogram length)", default=-1,guitype='intbox',row=9, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--bxsz", type=int,help="Box size of the particles for tracking. Default is 32. Maybe helpful to use a larger one for fiducial-less cases..", default=32, guitype='intbox',row=5, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--pk_maxval", type=float,help="Maximum Density value of landmarks (n sigma). Default is -5", default=-5.)
	parser.add_argument("--pk_mindist", type=float,help="Minimum distance between landmarks, as fraction of micrograph length. Default is 0.125", default=0.125,guitype='floatbox',row=11, col=0, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--correctrot", action="store_true",help="correct for global rotation and position sample flat in tomogram.", default=False,guitype='boolbox',row=12, col=0, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--normslice", action="store_true",help="normalize each 2D slice.", default=False,guitype='boolbox',row=13, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--filterto", type=float,help="filter to abs.", default=0.45,guitype='floatbox',row=13, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--extrapad", action="store_true",help="pad extra for tilted reconstruction. slower and cost more memory, but reduce boundary artifacts when the sample is thick", default=False,guitype='boolbox',row=14, col=0, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--rmbeadthr", type=float, help="Density value threshold (of sigma) for removing beads. high contrast objects beyond this value will be removed. default is -1 for not removing. try 10 for removing fiducials", default=-1,guitype='floatbox',row=14, col=1, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--threads", type=int,help="Number of threads", default=12,guitype='intbox',row=12, col=1, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tmppath", type=str,help="Temporary path", default=None)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)
	parser.add_argument("--noali", action="store_true",help="skip initial alignment", default=False)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)
	time0=time.time()
	
	#### deal with multiple inputs
	if options.alltiltseries:
		fld="tiltseries/"
		args=[fld+f for f in os.listdir(fld) if (
			f.endswith(".hdf") or f.endswith(".mrc") or f.endswith(".mrcs") or f.endswith(".st"))]
	
	if len(args)==1:
		inputname=args[0]
		print("Reading tilt series {}...".format(inputname))
	elif len(args)==0:
		print("No input. Exit.")
		return
	else:
		#### lauch subprocesses to reconstruct multiple tomograms
		print("Processing {} tilt series in sequence..".format(len(args)))
		cmd=sys.argv
		opt=' '.join([s for s in cmd if s.startswith("-")])
		opt=opt.replace("--alltiltseries","")
		for a in args:
			run("{} {} {}".format(cmd[0], a, opt))
			
		E2end(logid)
		return
		
	#### testing multithread fft...
	Util.init_threads(options.threads)
		
		
	#### parse options.
	itnum=[int(i) for i in options.niter.split(',')] ### number of iterations for each mag
	itnum+=[0,0,0,0,0] ### skip later iters when not specified.
	options.inputname=inputname
	if options.pk_mindist<=0: options.pk_mindist=0.12

	if options.tltrange != None:
		try: 
			options.tltrange = list(map(float,options.tltrange.split(",")))
		except:
			print("Could not interpret --tltrange: {} \nAn example of the correct format is: -50.0,50.0".format(options.tltrange))
			return
	else: 
		options.tltrange = [-90.0,90.0] # all plausible tilts

	dotpos=inputname.rfind('.')
	linepos=inputname.rfind('__')
	inputtag=inputname[linepos:dotpos]
	bname=base_name(inputname)

	options.basename=bname
	options.writetmp=not options.notmp
	
	#### read input. support both 2D image stack and 3D mrc/st files
	img=EMData(inputname,0)
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(inputname)
			
	#### this hopefully remove xray pixels or at least prevent them from doing too much damage.
	###  probably still better to do the import.
	for m in imgs: 
		m.process_inplace("threshold.clampminmax.nsigma", {"nsigma":10})
		m.process_inplace("normalize.edgemean")
	img=None
	
	options.apix_init=float(imgs[0]["apix_x"])
	
	#### need to make sure this works for images of all sizes (2k, 4k 8k)
	## the translation numbers used in this program are based on 2k tomograms. so binfac is the factor from the input to 2k images
	#binfac=max(1, int(np.round(imgs[0]["nx"]/2048.)))
	#options.binfac=binfac
	if imgs[0]["nx"]<=512:
		print("Tilt series image too small. Only support 2K or larger input images...")
		return
	
	elif imgs[0]["nx"]<=1024:
		imgs_1k=imgs_2k=imgs_4k=imgs
		itnum[3]=itnum[2]=0
		options.outsize="1k"
		
	elif imgs[0]["nx"]<=2048:
		#### 2k or smaller input. skip 4k refinement
		imgs_2k=imgs_4k=imgs
		itnum[3]=0
	elif imgs[0]["nx"]<=4096:
		#### 4k input
		imgs_4k=imgs
		imgs_2k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_4k]
			
	else:
		#### even larger images.. hopefully 8k at most
		bf=imgs[0]["nx"]//4096+1
		imgs_4k=[img.process("math.meanshrink", {"n":bf}).process("normalize.edgemean") for img in imgs]
		imgs_2k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_4k]
		imgs=imgs_4k
	
	if imgs[0]["nx"]>1024:
		imgs_1k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_2k]
	
	#### 500x500 images. this is used for the initial coarse alignment so we need to filter it a bit. 
	imgs_500=[]
	for p in imgs_1k:
		m=p.process("math.meanshrink", {"n":2})
		m.process_inplace("filter.ramp")
		m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		m.process_inplace("normalize.edgemean")
		m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":5})
		m["apix_x"]=m["apix_y"]=p["apix_x"]*2.
		imgs_500.append(m)

	num=options.num=len(imgs_500)
		
	#### create working directory
	## we do not create the folder until here so we do not end up with empty folders when the input is wrong..
	if options.writetmp:
		if options.tmppath:
			path=options.tmppath
		else:
			options.tmppath=path=make_path("tomorecon")
		print("Temporary files will be written in {}".format(options.tmppath))
		
	#### save 1k input only.
	if options.writetmp:
		inppath=os.path.join(options.tmppath,"tltseries_input.hdf")
		for i,m in enumerate(imgs_1k):
			m.write_image(inppath, i)
		
	
	if options.load:
		#### loading parameters from json file
		jsname=info_name(options.inputname)
		print("Loading parameters from {}...".format(jsname))
		js=js_open_dict(jsname)
		if "tlt_params" not in js:
			print("Failed to load saved parameterss. Exit.")
			return
		
		if js["tlt_file"]!=options.inputname:
			print("Warning. Tilt file in the saved file is not the same as the input. Something may be wrong..")
			print("\tprevious input :{}\t current input : {}".format(js["tlt_file"], options.inputname))
			
		tpm=np.array(js["tlt_params"])
		#tpm[:,:2]/=options.binfac
		ttparams=tpm.copy()
		js.close()
		tlts=ttparams[:,3].copy()
		#### some old version may not save loss
		try: 
			loss0=np.array(js["ali_loss"])
		except: 
			loss0=abs(ttparams[:,3])
			
		options.zeroid=zeroid=np.argmin(abs(tlts))
		
	else:
		#### determine alignment parameters from scratch
		if (options.rawtlt!=None and len(options.rawtlt)>0) :
			#### load raw tilt file
			tlts=np.loadtxt(options.rawtlt)
			
			if len(tlts)!=len(imgs_500):
				print("Number of tilt angles and tilt images do no match. Your inputs have {} tilt angles and {} tilt images.\nStopping reconstruction.".format(len(tlts),len(imgs_500)))
				return
		else: 
			#### parse tilt step
			tlts=np.arange(-len(imgs_2k)*options.tltstep//2,len(imgs_2k)*options.tltstep//2,options.tltstep)
		
		if options.writetmp: np.savetxt(os.path.join(options.tmppath,"rawtilt.txt"), tlts)
		
		#### we need a zero degree tilt as reference to position the tomogram
		zeroid=np.argmin(abs(tlts))
		if options.badzero:
			zeroid+=2
		options.zeroid=zeroid
		
		#### here we always assume the center tilt is at 0 degree
		tlts-=tlts[options.zeroid]
		
		ttparams=np.zeros((num, 5))
		
		#### estimate initial tilt axis by common line
		if options.tltax==None:
			#### do an initial course alignment before common line
			img_tali, pretrans=calc_global_trans(imgs_500, options)
			tltax=calc_tltax_rot(img_tali, options)
			options.tltax=tltax
		else:
			### it turns out there is a sign difference between serialem and eman...
			tltax=-options.tltax
			options.tltax=tltax
		print("tilt axis:  {:.2f}".format(tltax))
		
		if options.noali:
			print("Skipping coarse alignment...")
			
		else:
			#### a second round of course alignment, using the rotation axis information...
			img_tali, pretrans=calc_global_trans(imgs_500, options, tltax=options.tltax, tlts=tlts)
			if options.writetmp:
				for i,m in enumerate(img_tali):
					m.write_image(os.path.join(options.tmppath,"tltseries_transali.hdf"), i)
			
			#### this is the matrix that save the alignment parameters
			pretrans*=img_tali[0]["apix_x"]/options.apix_init  #### since the pretrans is calculated from 500x500 images..
		
			ttparams[:,0]=-pretrans[:,0] # tx
			ttparams[:,1]=-pretrans[:,1] # ty
			
		ttparams[:,2]=tltax # rot
		ttparams[:,3]=tlts.copy() # ytilt
		ttparams[:,4]=0 # off axis tilt
		loss0=abs(ttparams[:,3]) ### this is used to exclude bad tilt. in case the user ask 0 iterations..
	
	
	pks=np.zeros((options.npk, 3))
	#### pack parameters together so it is easier to pass around
	allparams=np.hstack([ttparams.flatten(), pks.flatten()])
	
	#### image scale, m3diter, fidkeep
	#### refinement sequence []:global tilt axis, 0:tx, 1:ty, 2:tilt axis, 3: tilt, 4: off axis tilt
	scaleiter=[(imgs_500, itnum[0], options.pkkeep*.8, [[0,1], [], [0,1],[], [0,1]]),
			   (imgs_1k, itnum[1], options.pkkeep*.9, [[0,1],[], [0,1],[3],[4],[2],[0,1]]),
			   (imgs_2k, itnum[2], options.pkkeep, [[0,1],[], [0,1],[3],[4],[2],[0,1]]),
			   (imgs_4k, itnum[3], options.pkkeep, [[0,1], [0,1],[3],[4],[2],[0,1]])]
	
	if options.writetmp:
		#### dump options to file
		options.cmd=' '.join(sys.argv)
		js=js_open_dict(os.path.join(path,"0_tomorecon_params.json"))
		js.update(vars(options))
		js.close()
		
		
		tpm=ttparams.copy()
		#tpm[:,:2]*=options.binfac
		tpm=np.hstack([np.arange(len(tpm))[:,None], tpm])
		np.savetxt(os.path.join(path,"tltparams_init.txt"), tpm, fmt="%.3f")
	
	yrot=0 #### this is to save the overall rotation of tomogram
	
	for niter, siter in enumerate(scaleiter):
		#### main refinement loop.
		###  we do a number of iterations for images of different sizes, starting from 500x500 images and gradually get to 4k images
		imgs_, n_m3d, options.fidkeep, rfseq = siter
		if n_m3d==0: continue
		apix=float(imgs_[0]["apix_x"])
		binx=np.round(apix/options.apix_init)
		print("\n******************************")
		print("Iteration {}. Refining alignment on bin{:.0f} images...".format(niter, binx))
		print("Image size {} x {}, Apix {:.2f}".format(imgs_[0]["nx"], imgs_[0]["ny"], apix))
		if niter<2:
			print("Low resolution mode. Using peak positon of landmarks for alignment")
		else:
			print("High resolution mode. Using center of mass of landmarks for alignment")
			
		if options.writetmp:
			name_tomo=os.path.join(path,"tomo_{:02d}.hdf".format(niter))
			name_sample=os.path.join(path,"samples_{:02d}.hdf".format(niter))
			name_ali=os.path.join(path,"ali_{:02d}.hdf".format(niter))
			name_ptclali=os.path.join(path,"ptclali_{:02d}.hdf".format(niter))
		else:
			name_tomo=name_sample=name_ali=name_ptclali=None
		
		#### make tomogram loop
		for m3diter in range(n_m3d):

			
			#### make initial tomogram for landmark search, always use 500x500
			threed=make_tomogram(imgs_500, ttparams, options, outname=name_tomo, errtlt=loss0)

			pks=find_landmark(threed, options)
			allparams=np.hstack([ttparams.flatten(), pks.flatten()])
			
			#### do two round of generating 3D maps of landmarks to refine the location of them
			allparams,smp=make_samples(imgs_, allparams, options, refinepos=True);
			allparams,smp=make_samples(imgs_, allparams, options, refinepos=True);

			if niter==0 and m3diter==0 and options.writetmp:
				make_samples(imgs_, allparams, options, outname=os.path.join(path,"samples_init.hdf"), refinepos=True);
				
			#ptclpos=ali_ptcls(imgs_, allparams, options, outname=name_ptclali, doali=True)
			#exit()

			#### Now actually refine the alignment parameters using the landmarks.
			for idx in rfseq:
				allparams=refine_one_iter(imgs_, allparams, options, idx)

			ttparams, pks=get_params(allparams, options)
			
			
			
			allparams=np.hstack([ttparams.flatten(), pks.flatten()])
		
		ttparams, pks=get_params(allparams, options)
		if options.writetmp:
			#### make sample output and aligned tilt series
			make_samples(imgs_, allparams, options, outname=name_sample, refinepos=False);
			make_ali(imgs_2k, ttparams, options, outname=name_ali)

		#### a final round to calculate the average loss per tilt
		ptclpos=ali_ptcls(imgs_, allparams, options, outname=name_ptclali, doali=True)
		loss0=np.zeros(num)
		for nid in range(num):
			loss0[nid]=get_loss_pm([0], nid, allparams, options, [0], ptclpos)
			
		### In the first round, position the tomogram so the features are relatively flat.
		if niter==0 and options.correctrot:
			threed=make_tomogram(imgs_500, ttparams, options, errtlt=loss0)
			ttparams=fix_rotation(threed, ttparams) 
			zeroid=options.zeroid=np.argmin(abs(ttparams[:,3]))
			ttparams[:,3]-=ttparams[options.zeroid,3]

		print("Iteration {} finished. Final average loss {:.2f} nm".format(
			niter, np.mean(loss0)))
	
		#### always save parameters at the full scale (4k)
		tpm=ttparams.copy()
		#tpm[:,:2]*=options.binfac
		tpm=np.hstack([np.arange(len(tpm))[:,None], tpm])
		if options.writetmp:
			np.savetxt(os.path.join(path,"landmarks_{:02d}.txt".format(niter)), pks, fmt="%.1f")
			np.savetxt(os.path.join(path,"tltparams_{:02d}.txt".format(niter)), tpm, fmt="%.3f")
			np.savetxt(os.path.join(path,"loss_{:02d}.txt".format(niter)), np.vstack([np.arange(len(loss0)), loss0]).T, fmt="%.2f")
		
	options.loss0=loss0
	#### alignment finish. now save output
	if options.outsize=="2k":
		imgout=imgs_2k
	elif options.outsize=="4k":
		options.bytile=True
		imgout=imgs_4k
	else:
		imgout=imgs_1k
		
	
	if options.rmbeadthr>0:
		remove_beads(imgs_500, imgout, ttparams, options)
	
	#### only clip z axis at the end..
	if options.bytile:
		threed=make_tomogram_tile(imgout, ttparams, options, errtlt=loss0, clipz=options.clipz)
	else:
		threed=make_tomogram(imgout, ttparams, options, errtlt=loss0, clipz=options.clipz)

	if options.writetmp:
		threed.write_image(os.path.join(path,"tomo_final.hdf"))
		make_ali(imgout, ttparams, options, outname=os.path.join(path,"tiltseries_ali.hdf"))

	#### write to the tomogram folder
	try: os.mkdir("tomograms")
	except: pass
	sfx=""
	bf=int(np.round(imgout[0]["apix_x"]/options.apix_init))
	if bf>1:
		sfx+="__bin{:d}".format(int(bf))
		
	threed["ytilt"]=yrot
	tomoname=os.path.join("tomograms", options.basename+sfx+".hdf")
	threed.write_image(tomoname)
	print("Tomogram written to {}".format(tomoname))
	
	#### save alignemnt parameters to info file
	tpm=ttparams.copy()
	#tpm[:,:2]*=options.binfac
	js=js_open_dict(info_name(tomoname))
	js["tlt_params"]=tpm.tolist()
	js["tlt_file"]=options.inputname
	js["ali_loss"]=loss0.tolist()
	js["apix_unbin"]=options.apix_init
	js.close()
	
	dtime=time.time()-time0
	print("Finished. Total time: {:.1f}s".format(dtime))

	E2end(logid)


def remove_beads(imgs_500, imgout, ttparams, options):
	
	print("Removing high contrast objects...")
	threed=make_tomogram(imgs_500, ttparams, options, errtlt=options.loss0, clipz=options.clipz)
	threed.process_inplace("math.meanshrink",{"n":2})
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./100})
	threed.process_inplace("filter.highpass.gauss",{"cutoff_freq":1./50})
	threed.process_inplace("normalize")
	if options.writetmp:
		threed.write_image(os.path.join(options.tmppath, "tomo_rmbead.hdf"))
	
	vthr=options.rmbeadthr
	img=threed.numpy().T.copy()
	
	img=-img
	img[img<0]=0
	lb, nlb= sciimg.label(img>vthr)
	idx=np.arange(1,nlb, dtype=int)
	
	sm=np.array(sciimg.sum(img, lb, idx))/vthr
	idx=idx[sm>6]
	cnt=np.array(sciimg.center_of_mass(img, lb, idx))
	
	
	scale=threed["apix_x"]/options.apix_init
	
	pts=(cnt-np.array(img.shape)/2.)*scale
	
	realnpk=options.npk
	options.npk=len(pts)
	allparams=np.hstack([ttparams.flatten(), pts.flatten()])
	apm, smp=make_samples(imgs_500, allparams, options)
	
	nx, ny=imgout[0]["nx"],imgout[0]["ny"]
	
	tpms=ttparams.copy()
	scale=imgout[0]["apix_x"]/options.apix_init
	tpms[:,:2]/=scale
	pts/=scale
	
	if options.writetmp:
		fname=os.path.join(options.tmppath, "ptcls_rmbead.hdf")
		try:
			os.remove(fname)
		except:
			pass
	pad=good_boxsize(options.bxsz*2)
	for n, tpm in enumerate(tpms):
		
		plst=[]
		for i, pk in enumerate(pts):
			if smp[i]["minimum"]>-vthr:
				continue
			
			pxf=get_xf_pos(tpm, pk)
			
			
			pxf=[int(round(pxf[0]))+nx//2-pad//2, int(round(pxf[1]))+ny//2-pad//2]
			
			if min(pxf)<pad*.8 or min(nx-pxf[0], ny-pxf[1])<pad*.8:
				continue
			
			if len(plst)>0:
				dst=scidist.cdist(plst, [pxf])
				if np.min(dst)<pad/3:
					continue
			plst.append(pxf)
			e=imgout[n].get_clip(Region(pxf[0], pxf[1], pad, pad))
			e["nid"]=n
			e["pid"]=i
			
			
			m=e.process("normalize.edgemean")
			m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":4})
			m.process_inplace("mask.soft",{"outer_radius":pad//3, "width":pad//4})
			m.mult(-1)
			m.process_inplace("threshold.belowtozero")
			m.mult(m)
			m.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./100})
			m.add(-1)
			m.process_inplace("threshold.belowtozero")
			
			m.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./100})
			if m["maximum"]>0:
				m.mult(5./m["maximum"])
			else:
				continue
				
			m.process_inplace("threshold.clampminmax",{"maxval":1})
			
			if m["mean"]>0.2:
				continue
			
			m.process_inplace("mask.soft",{"outer_radius":int(pad*.3),"width":8})
			
			a0=e*(1-m)
			
			a1=m.copy()
			a1.to_zero()
			#a1.add(a0["mean_nonzero"])
			a1.process_inplace("math.addnoise",{"noise":1})
			a1.process_inplace("normalize")
			
			a1=a1*a0["sigma_nonzero"]+a0["mean_nonzero"]
			
			#if m["mean"]>0.05:
				#a1.process_inplace("filter.matchto",{"to":e})
			#else:
				#a1.process_inplace("filter.matchto",{"to":a0})
			a1.mult(m)
			
			e1=a0+a1    
			
			if options.writetmp:
				e.process("normalize").write_image(fname,-1)
				m.write_image(fname,-1)
				e1.process("normalize").write_image(fname,-1)
			
			
			imgout[n].insert_clip(e1, pxf)
			
			
			
			
	options.npk=realnpk
	

#### find the optimal global rotation of the tomogram so the high contrast features lies flat
def fix_rotation(threed, ttparams):
	
	print("Positioning tomogram...")
	m=threed.copy()
	m.process_inplace("math.meanshrink",{'n':2})
	m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
	m.process_inplace("filter.highpass.gauss",{"cutoff_abs":.05})
	m.process_inplace("math.absvalue")
	m.process_inplace("normalize")
	m.process_inplace("math.maxshrink",{'n':2})
	m.process_inplace("math.maxshrink",{'n':2})
	img=m.numpy().copy()
	
	pval=img.flatten()
	thr=-np.sort(-pval)[500]
	pts=np.array(np.where(img>thr)).T
	pts=pts[:,::-1]
	pca=PCA(3)
	pca.fit(pts);
	c=pca.components_
	t=Transform()
	cc=c[2]
	cc*=np.sign(cc[2])
	t.set_rotation(c[2].tolist())
	t.invert()
	xyz=t.get_params("xyz")
	xyz["ztilt"]=0
	print("  Global transform : xtilt {:.02f}, ytilt {:.02f}".format(xyz["xtilt"], xyz["ytilt"]))
	
	t=Transform(xyz)

	p1=np.array([t.transform(p.tolist()) for p in pts])

	ttparams_new=ttparams.copy()
	for i,tpm in enumerate(ttparams):
		xf=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4], "ztilt":tpm[2]})
		xf=xf*t.inverse()
		x=xf.get_params("xyz")
		ttparams_new[i][2]=x["ztilt"]
		ttparams_new[i][3]=x["ytilt"]
		ttparams_new[i][4]=x["xtilt"]
	return ttparams_new

#### unpack parameters. I really shouldn't need this but it is risky to remove the convertion
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
def calc_global_trans(imgs, options, excludes=[], tltax=None,tlts=[]):
	print("Doing coarse translational alignment...")

	num=len(imgs)
	sz=min(imgs[0]["nx"], imgs[0]["ny"])

	imgout=[0]*num
	rmsk=sz//3
	e0=imgs[options.zeroid].copy()
	e0.clip_inplace(Region(e0["nx"]//2-sz//2, e0["ny"]//2-sz//2, sz,sz))
	e0.process_inplace("mask.soft",{"outer_radius":sz//4,"width":sz//10})
		
	e0["xform.align2d"]=Transform()
	imgout[options.zeroid]=e0

	pretrans=np.zeros((num,2))
	#### start from the center tilt and move toward both directions
	for dr in [-1,1]:
		for i in range(num):

			nid=options.zeroid+(i+1)*dr
			if nid>=num or nid<0: continue
			if nid in excludes: continue
			e0=imgout[options.zeroid+i*dr]
			e1=imgs[nid].copy()
			lastx=pretrans[options.zeroid+i*dr]
			if tltax==None:
				e1.process_inplace("mask.soft",{"outer_radius":sz//4,"width":sz//10})
				
			else:
				rx=np.cos(tlts[nid]*np.pi/180.)
				msk=EMData(e1["nx"], e1["ny"])
				msk.process_inplace("testimage.ellipsoid", {"a":rx*sz//4,"b":sz//4,"fill":1,"transform":Transform({"type":"2d","alpha":tltax})})
				msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.01})
				#m=msk.process("mask.gaussian.nonuniform",{"radius_x":sz//4, "radius_y":ry*sz//4, "gauss_width":0.4})
				e1.mult(msk)
				
			e1.clip_inplace(Region(e1["nx"]//2-sz//2-lastx[0], e1["ny"]//2-sz//2-lastx[1], sz,sz))
			#e1.process_inplace("filter.highpass.gauss",{"cutoff_pixels":10})

			e1a=e1.align("translational", e0)
			xf=e1a["xform.align2d"]
			if i>1:
				e00=imgout[options.zeroid+(i-1)*dr]
				e1b=e1.align("translational", e00)
				c0=e1a.cmp("ccc", e0)
				c1=e1b.cmp("ccc", e00)
				if c1<c0*4:
					print("Bad tilt image at {}  ".format(nid), c0, c1)
					xf=e1b["xform.align2d"]
				

			xf.translate(lastx[0], lastx[1], 0)
			
			#e1=imgs[nid].copy()
			#e1.transform(xf)
			#e1.clip_inplace(Region((e1["nx"]//2)-(sz//2), (e1["ny"]//2)-(sz//2), sz,sz))
			##e1.process_inplace("filter.highpass.gauss",{"cutoff_pixels":10})
			#e1.process_inplace("mask.soft",{"outer_radius":sz//4,"width":sz//10})

			imgout[nid]=e1a
			ts=xf.get_trans()
			if options.verbose: print("\t{:d}: {:.1f}, {:.1f}".format(nid, ts[0], ts[1]))
			pretrans[nid, 0]=ts[0]
			pretrans[nid, 1]=ts[1]
			
		
	return imgout,pretrans


#### Find tilt axis by common line. copied from the original e2tomogram. Not very stable but mostly fine.. 
def calc_tltax_rot(imgs, options):
	print("Calculating tilt axis rotation...")
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
	sm=np.abs(sm[:,sz//2:])
	#print(np.max(sm), np.min(sm))
	rr=np.arange(min(sm.shape[1], sz*.25), dtype=float)
	angs=np.arange(0., 180, .5)
	vs=[]
	for ang in angs:
		a=ang/180.*np.pi
		pts=[np.round(rr*np.sin(a)).astype(int), np.round(rr*np.cos(a)+old_div(sz,2)).astype(int) ]
		v=sm[pts[1], pts[0]]
		vs.append(np.mean(v))
		
	vs[0]=vs[180]=0 #### sometimes there are false positive vertical peaks
	vs=np.array(vs)
	
	#### down weight the 180 degree end a bit to avoid the 0/180 ambiguity..
	nw=30
	vs[-nw:]*=(1-0.3*np.arange(1,nw+1, dtype=float)/nw)
	
	tltax=angs[np.argmax(vs)]
	e=from_numpy(sm)
	if options.writetmp: 
		e.write_image(os.path.join(options.tmppath,"commonline.hdf"))
		np.savetxt(os.path.join(options.tmppath,"tltrot.txt"), np.vstack([angs, vs]).T)
		
	return tltax

#### subthread for making tomogram by tiles. similar to make_tomogram, just for small cubes
def make_tile(args):
	jsd, imgs, tpm, sz, pad, stepx, stepy, outz,options=args
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"gauss_2"})
	recon.setup()

	for i in range(len(imgs)):
		t=tpm[i]
		m=imgs[i]
		if m["nx"]==1:
			continue
		m.process_inplace("filter.ramp")
		m.process_inplace("xform",{"alpha":-t[2]})
		xf=Transform({"type":"xyz","ytilt":t[3],"xtilt":t[4]})

		dy=(pad//2)-np.cos(t[3]*np.pi/180.)*pad/2
		msk=EMData(pad, pad)
		msk.to_one()
		edge=(sz//10)
		msk.process_inplace("mask.zeroedge2d",{"x0":dy+edge, "x1":dy+edge, "y0":edge, "y1":edge})
		msk.process_inplace("mask.addshells.gauss",{"val1":0, "val2":edge})
	
		m.mult(msk)
		mp=recon.preprocess_slice(m, xf)
		recon.insert_slice(mp,xf,1)
		
	
	threed=recon.finish(True)
	threed.process_inplace("math.gausskernelfix",{"gauss_width":4.0})
	#threed.write_image("tmp3d00.hdf", -1)
	threed.clip_inplace(Region((pad-sz)//2, (pad-sz)//2, (pad-outz)//2, sz, sz, outz))
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	#threed.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
	jsd.put( [stepx, stepy, threed])
	
	return


#### make tomogram by tiles
#### this is faster and has less artifacts. but takes a lot of memory (~4x the tomogram)
def make_tomogram_tile(imgs, tltpm, options, errtlt=[], clipz=-1):
	
	num=len(imgs)
	scale=imgs[0]["apix_x"]/options.apix_init
	if imgs[0]["nx"]<=1024*1.1:
		b=1
	elif imgs[0]["nx"]<=2048*1.1:
		b=2
	else:
		b=4
		#print("tiling only support for 1k and 2k tomograms...")
		#return make_tomogram(imgs, tltpm, options, errtlt=errtlt)
	
	print("Making bin{:d} tomogram by tiling...".format(int(np.round(scale))))
	tpm=tltpm.copy()
	tpm[:,:2]/=scale
	
	if len(errtlt)==0:
		errtlt=np.zeros(num)
		nrange=list(range(num))
	else:
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]

	print("Using {} out of {} tilts..".format(len(nrange), num))

	outxy=1024*b
	sz=max(int(clipz*0.8), 256*b) #### this is the output 3D size 
	step=sz//2 #### distance between each tile
	if options.extrapad:
		pad=good_boxsize(sz*2) #### this is the padded size in fourier space
	else:
		pad=good_boxsize(sz*1.4) #### this is the padded size in fourier space
	
	if clipz>0:
		outz=clipz
	else:
		outz=sz#good_boxsize(sz*1.2)

	#### we make 2 tomograms with half a box shift and average them together to compensate for boundary artifacts.
	mem=(outxy*outxy*outz*4+pad*pad*pad*options.threads*4)
	print("This will take {}x{}x{}x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outxy, outxy, outz, pad, pad, pad,options.threads, mem/1024**3))
	
	
	full3d=EMData(outxy, outxy, outz)
	
	jsd=queue.Queue(0)
	jobs=[]
	nstep=int(outxy/step/2)
	for stepx in range(-nstep,nstep+1):
		#### shift y by half a tile
		for stepy in range(-nstep+stepx%2,nstep+1,2):
			tiles=[]
			for i in range(num):
				if i in nrange:
					t=tpm[i]
					pxf=get_xf_pos(t, [stepx*step,stepy*step,0])
					img=imgs[i]
					m=img.get_clip(Region(img["nx"]//2-pad//2+pxf[0],img["ny"]//2-pad//2+pxf[1], pad, pad), fill=0)
					tiles.append(m)
				else:
					tiles.append(EMData(1,1))

			jobs.append((jsd, tiles, tpm, sz, pad, stepx, stepy, outz, options))
	
	thrds=[threading.Thread(target=make_tile,args=([i])) for i in jobs]
	print("now start threads...")
	thrtolaunch=0
	tsleep=threading.active_count()
	
	#### non-round fall off. this is mathematically correct but seem to have grid artifacts
	#f=np.zeros((sz,sz))
	x,y=np.indices((sz,sz),dtype=float)/sz-.5
	#f=.25-(x**2+y**2)/2 + ((abs(x)-0.5)**2+(abs(y)-0.5)**2)/2
	f=1+np.exp(-(x**2+y**2)/0.1) - np.exp(-((abs(x)-0.5)**2+(abs(y)-0.5)**2)/0.1)
	f+=1.5
	f3=np.repeat(f[None, :,:], outz, axis=0)
	msk=from_numpy(f3).copy()
	#####

	
	#msk=EMData(sz,sz,outz)
	#msk.to_one()
	#msk.process_inplace("mask.soft",{"outer_radius":sz//4, "width":sz//6})
	#msk.process_inplace("mask.gaussian",{"outer_radius":sz//4})
	#msk.process_inplace("mask.poly",{"2d":True, "k4":1/4, "k2":-1, "k0":1})
	#msk.add(0.1)
	#msk.write_image("tmp03_msk.hdf")
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep or jsd.empty()==False:
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(.1)
		
		if not jsd.empty():
			stepx, stepy, threed=jsd.get()
			threed.mult(msk)
			#### insert the cubes to corresponding tomograms
			#full3d[stepx%2].insert_clip(
				#threed,
				#(int(stepx*step+outxy//2-sz//2),
				#int(stepy*step+outxy//2-sz//2), 
				#outz//2-threed["nz"]//2))
				
			full3d.insert_scaled_sum(
				threed,(int(stepx*step+outxy//2),int(stepy*step+outxy//2), outz//2))
				
			
			#wt.insert_scaled_sum(msk,(int(stepx*step+outxy//2),int(stepy*step+outxy//2), outz//2))
				
				
	for t in thrds: t.join()
	
	#full3d[0].write_image("tmp00_full.hdf")
	#full3d[1].write_image("tmp01_full.hdf")
	#wt.write_image("tmp02_wt.hdf")
	#print("lalala")
	#wt.add(.01)
	#full3d=full3d[0]+full3d[1]
	#full3d.div(wt)
	#full3d.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	full3d.process_inplace("normalize")
	
	#### skip the tomogram positioning step because there is some contrast difference at the boundary that sometimes breaks the algorithm...
	full3d["zshift"]=0
	
	apix=imgs[0]["apix_x"]
	full3d["apix_x"]=full3d["apix_y"]=full3d["apix_z"]=apix
	print("Done. Now writting tomogram to disk...")
	if options.normslice:
		full3d.process_inplace("normalize.rows")
	return full3d

#### reconstruct tomogram...
def make_tomogram(imgs, tltpm, options, outname=None, padr=1.2,  errtlt=[], clipz=-1):
	num=len(imgs)
	scale=imgs[0]["apix_x"]/options.apix_init
	print("Making bin{:d} tomogram...".format(int(np.round(scale))))
	ttparams=tltpm.copy()
	ttparams[:,:2]/=scale

	#### sort tilt by loss to exclude the worst ones
	if len(errtlt)==0:
		errtlt=np.zeros(num)
		nrange=list(range(num))
	else:
		for nid in range(num):
			ytlt=ttparams[nid][3]
			if ytlt < options.tltrange[0] or ytlt > options.tltrange[1]:
				errtlt[nid] = np.inf # ensure this tilt is excluded if outside desired tilt range
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]

	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	outxy=good_size(max(nx, ny))

	pad=good_size(outxy*padr)
	#############
	#clipz=options.clipz
	zthick=good_size(max(clipz*1.2, pad//2))
	if options.verbose:
		print("\t Image size: {:d} x {:d}".format(nx, ny))
		print("\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick))
		
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":"gauss_2"})
	recon.setup()
	
	#### prepare jobes
	jobs=[]
	for nid in range(num):
		exclude= nid not in nrange
		tpm=ttparams[nid]
		pxf=get_xf_pos(ttparams[nid], [0,0,0])
		xform={"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4], "tx":pxf[0], "ty":pxf[1]}
		jobs.append([nid,imgs[nid],  recon, pad, xform, exclude, options])
			
	#### starting threads
	thrds=[threading.Thread(target=reconstruct,args=(i)) for i in jobs]
	thrtolaunch=0
	tsleep=threading.active_count()
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose : print("Inserting slice {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	for t in thrds: t.join()

	threed=recon.finish(True)
	
	threed.process_inplace("math.gausskernelfix",{"gauss_width":4.0})
	threed.process_inplace("normalize")
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	#print(threed["nx"], threed["ny"], threed["nz"])
	
	if clipz<0:
		threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), 0, outxy, outxy, zthick))
	else:
		threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), (zthick-clipz)//2, outxy, outxy, clipz))
	threed["zshift"]=0

	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix

	if options.normslice and outxy>1000:
		threed.process_inplace("normalize.rows")
	if outname:
		threed.write_image(outname)
		if options.verbose: print("Map written to {}.".format(outname))

	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform,  exclude, options):
	m=img.copy()
	#### the ramp filter and decay edge helps soften the edge artifacts
	m.process_inplace("filter.ramp")
	m.process_inplace("normalize")
	m.process_inplace("mask.decayedge2d", {"width":int(pad//20)})
	p2=m.get_clip(Region(m["nx"]//2-pad//2,m["ny"]//2-pad//2, pad, pad), fill=0)
	#### give up on the subpixel accuracy since it does not really matter here..
	p2.translate(-int(xform["tx"]), -int(xform["ty"]), 0)
	p2.rotate(-xform["ztilt"],0,0)
	xf=Transform({"type":"xyz","ytilt":xform["ytilt"],"xtilt":xform["xtilt"]})
	
	#### mask out the extra information on the edge of high tilt
	dy=p2["nx"]//2-np.cos(xform["ytilt"]*np.pi/180.)*m["nx"]//2
	msk=p2.copy()
	msk.to_one()
	edge=int(old_div(pad,20))
	msk.process_inplace("mask.zeroedge2d",{"x0":dy+edge, "x1":dy+edge, "y0":edge, "y1":edge})
	msk.process_inplace("mask.addshells.gauss",{"val1":0, "val2":edge})
	p2.mult(msk)
	
	if not exclude:
		p3=recon.preprocess_slice(p2, xf)
		recon.insert_slice(p3,xf,1)

#### make aligned tilt series with xform.projection in their headers 
#### one should be able to do make3dpar on them to get the 3D tomogram
def make_ali(imgs, tpm, options, outname=None):
	if outname==None:
		return
	
	scale=old_div(imgs[0]["apix_x"],options.apix_init)
	ttparams=tpm.copy()
	ttparams[:,:2]/=scale
	
	try:os.remove(outname)
	except:pass
	pad=imgs[0]["nx"]*1.
	mskrd=old_div(min(imgs[0]["nx"],imgs[0]["ny"]),2)
	for nid, im in enumerate(imgs):
		tpm=ttparams[nid]

		pxf=get_xf_pos(ttparams[nid], [0,0,0])
		m=im.process("normalize")
		p2=m.get_clip(Region(old_div(m["nx"],2)-old_div(pad,2),old_div(m["ny"],2)-old_div(pad,2), pad, pad), fill=0)
		po=p2.copy()
		po.translate(-pxf[0], -pxf[1], 0)
		po.rotate(-tpm[2],0,0)
		xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4]})
		po["xform.projection"]=xform
		po.write_image(outname, nid)

#### search for alignment landmarks in the tomogram
def find_landmark(threed, options):
	print("Searching for landmarks...")
	#threed0=threed.process("normalize.rows")
	#### use minshrink so we keep the minimas
	threedtiny=threed.process("math.minshrink", {"n":2})
	threedtiny.process_inplace("normalize")
	threedtiny.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
	mapnp=threedtiny.numpy().copy()
	asrt= np.argsort(mapnp.flatten())

	#### go through every point starting from the darkest ones.
	### there should be a much faster way but we are only doing it on ~250 cubes so it is not too slow
	pts=[]
	dthr=options.pk_mindist*float(threedtiny["nx"])
	vthr=options.pk_maxval
	for i in range(len(asrt)):
		aid=asrt[i]
		pt=np.unravel_index(aid, mapnp.shape)

		if len(pts)>0:
			dst=scidist.cdist(pts, [pt])
			if np.min(dst)<dthr:
				continue

		pts.append(pt)
		if mapnp[pt]>vthr:
			break
		
		#### generate a bit more landmarks for some randomness
		if len(pts)>=options.npk*1.2:
			break

	if len(pts)<options.npk:
		print("Found only {} landmarks".format(len(pts)))
		options.npk=len(pts)
	else:
		#### I commented out the shuffling for stable testing. maybe should add it back sometime..
		#np.random.shuffle(pts)
		pts=pts[:options.npk]

	pks=np.array(pts)-np.array(mapnp.shape)/2.
	pks=pks[:,::-1]
	#### mult 2 since we bin 2 in the begining.
	scale=float(threed["apix_x"])/options.apix_init*2.
	pks*=scale

	return pks

#### make sample sub-tomograms of landmarks and center them
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
	#### do this slightly differently at different image size
	lowres=(scale>3)
	if len(errtlt)==0:
		nrange=list(range(num))
	else:
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]
	bx=options.bxsz//2
	if abs(scale-2)<.1:
		bx=int(bx*1.5)
	elif abs(scale-1)<.1:
		bx*=2
	#if not lowres:
		#bx=int(bx*2/(scale))
	#print("scale{}, box size {}".format(scale, bx*2))
	
	samples=[]

	for pid in range(npk):
		pad=good_size(bx*4)
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad]})
		recon.setup()

		for nid in nrange:

			tpm=ttparams[nid]

			pxf=get_xf_pos(ttparams[nid], pks[pid])

			pxf[0]+=imgs[nid]["nx"]//2
			pxf[1]+=imgs[nid]["ny"]//2
			
			xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1]})
			e=imgs[nid].get_rotated_clip(xf,(pad,pad,1))

			#e=imgs[nid].get_clip(Region(pxf[0]-pad/2,pxf[1]-pad/2, pad, pad))
			e.process_inplace("normalize")
			p2=e
			rot=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4]})
			p3=recon.preprocess_slice(p2, rot)
			recon.insert_slice(p3,rot,1)
		bxcr=np.round(pks[pid]).astype(int).tolist()
		threed=recon.finish(True)
		threed=threed.get_clip(Region((pad-bx*2)//2,(pad-bx*2)//2,(pad-bx*2)//2,bx*2,bx*2,bx*2))
		threed.process_inplace("normalize")
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=imgs[0]["apix_x"]
		threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		threed.process_inplace("mask.soft",{"outer_radius":-1})
		pj=threed.project("standard", Transform({"type":"eman", "alt":90}))
		pj["box"]=bxcr
		samples.append(threed)

		#### center the subvolumes. in low resolution mode by minima position; in high resolution mode by center of mass
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

	if outname: 
		print("Landmark samples written to {}.".format(outname))
	
	return allparams, samples


#### find center of landmark subtomogram. 
def get_center(img, lowres=True):
	e=img.copy()
	bx=old_div(e["nx"],2)
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

#### this is the main function that calculate the offset of landmarks in tilt series given images and alignment parameters 
def ali_ptcls(imgs, allpms, options, outname=None, doali=True):
	zeroid=options.zeroid
	num=options.num
	nrange=np.hstack([np.arange(zeroid, num), np.arange(zeroid, -1, -1)])
	ttparams, pks=get_params(allpms, options)
	scale=imgs[0]["apix_x"]/options.apix_init
	ttparams[:,:2]/=scale
	pks/=scale
	prange=np.arange(options.npk)

	k=0
	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	fidptcls=[]
	bx=options.bxsz//2
	apix=imgs[0]["apix_x"]
	#### use a larger a box size at low resolution mode
	lowres=(scale>3)
	if abs(scale-2)<.1:
		bx=int(bx*1.5)
	elif abs(scale-1)<.1:
		bx*=2
	#if not lowres:
		#bx=int(bx*2/(scale))
	
	#### this is the matrix to return containing the offset of each landmark at each tilt
	ptclpos=[]
	for pid in prange:
		#### loop through all particles
		trans=np.zeros((num,2))
		fid=[]
		ppos=np.zeros((num,2))
		for ii,nid in enumerate(nrange):
			#### start from center tilt and go both directions
			pxf=get_xf_pos(ttparams[nid], pks[pid])
			pxf[0]+=old_div(nx,2)
			pxf[1]+=old_div(ny,2)

			if nid!=zeroid:
				tlast=trans[nrange[ii-1]]
				pxf[0]-=tlast[0]
				pxf[1]-=tlast[1]
			else:
				tlast=np.array([0,0])

			xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1]})
			e=imgs[nid].get_rotated_clip(xf,(bx*2,bx*2,1)).process("normalize")
			e["apix_x"]=e["apix_y"]=e["apix_z"]=imgs[0]["apix_x"]

			ts=[0,0]
			if doali:# and nid!=zeroid:
				#### still get the center of particle differently at low and high res mode
				tx=get_center(e, lowres)
				trans[nid]=tlast+np.array([tx[0], tx[1]])
				
				if nid!=zeroid:
					if outname:
						xf=Transform({"type":"2d","tx":pxf[0]-tx[0],"ty":pxf[1]-tx[1]})
						e=imgs[nid].get_rotated_clip(xf,(bx*2,bx*2,1)).process("normalize")
						e["apix_x"]=e["apix_y"]=e["apix_z"]=imgs[0]["apix_x"]

			if outname:
				e["score"]=trans[nid].tolist()
				e["pid"]=pid
				e["nid"]=nid
				e.write_image(outname, nid+pid*num)
			ppos[nid]=np.array(pxf-trans[nid]-[old_div(nx,2), old_div(ny,2)])

		ptclpos.append(ppos)

	ptclpos=np.array(ptclpos)*scale
	return ptclpos

#### one refinement run
def refine_one_iter(imgs, allparams, options, idx=[]):
	#### refine landmark location -> calculate landmark offset per tilt -> refine alignment parameters to minimize offset

	apms,smp=make_samples(imgs, allparams, options, refinepos=True);
	ttparams, pks=get_params(apms, options)
	ptclpos=ali_ptcls(imgs, apms, options, doali=True)

	#### refine different parameters given different idx input
	pmlabel=["trans_x", "trans_y", "tlt_z", "tilt_y", "tilt_x"]
	
	if len(idx)==0:
		#### this refines the global rotation around z axis 
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
		#### refine a specific parameter
		ttpm_new=ttparams.copy()
		num=len(ttparams)
		loss=np.zeros(num)
		loss0=np.zeros(num)
		for nid in range(num):
			tminit=np.zeros(len(idx))
			#### simply call a scipy minimizer. Powell seem to work fine.
			res=minimize(get_loss_pm, tminit,(nid,apms, options, idx, ptclpos) , method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
			ttpm_new[nid, idx]+=res.x*.5
			loss[nid]=get_loss_pm(res.x*.5, nid, apms, options, idx, ptclpos)
			loss0[nid]=get_loss_pm(tminit, nid, apms, options, idx, ptclpos)

		ipm=[pmlabel[i] for i in idx]
		print("refine {}, loss: {:.2f} -> {:.2f}".format(ipm, np.mean(loss0), np.mean(loss)))
		ttparams=ttpm_new

	allparams=np.hstack([ttparams.flatten(), pks.flatten()])
	return allparams

#### function for scipy optimizer
def get_loss_pm(pm, nid, apms, options, idx=[2,3,4], ptclpos=[]):
	ttparams, pks=get_params(apms, options)
	tpm=ttparams[nid].copy()
	tpm[idx]+=pm
	ps=np.array([get_xf_pos(tpm, p) for p in pks])
	dst=np.sqrt(np.sum((ps-ptclpos[:,nid])**2, axis=1))
	dst=np.mean(dst[np.argsort(dst)[:int(len(dst)*options.fidkeep)]])

	#### return distance in nm
	return dst*options.apix_init/10.

#### function for scipy optimizer for global rotation
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
