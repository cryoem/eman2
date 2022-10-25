#!/usr/bin/env python
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
		e2tomogram.py <tilt series> <tilt series> <tilt series> ... --rawtlt <folder with .tlt files> [options]
		e2tomogram.py <tilt series stack> --tltstep <angle between tilts> [options]
	
	Note: Tiltseries must have the correct Apix values in their headers.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="tiltseries",help="Specify the tilt series you intend to reconstruct.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='tiltseries',multiselect=True)", filecheck=False, row=0, col=0,rowspan=1, colspan=2,nosharedb=True,mode="easy")
	
	parser.add_argument("--alltiltseries", action="store_true",help="Use all tilt series in the folder. Acceptable file extensions include hdf, mrc, mrcs, st.", default=False,guitype='boolbox',row=1, col=0, rowspan=1, colspan=1,mode="easy")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Specify either zeroid/tltstep OR rawtlt:", row=2, col=0, rowspan=1, colspan=2,mode="easy")

	parser.add_argument("--zeroid", type=int,help="Index of the center tilt. Ignored when rawtlt is provided.", default=-1,guitype='intbox',row=3, col=1, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tltstep", type=float,help="Step between tilts. Ignored when rawtlt/mdoc is provided. Set to 0 if tilt present in header.", default=2.0,guitype='floatbox',row=3, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--rawtlt", type=str,help="Specify a text file contains raw tilt angles. Will look for files with the same name as the tilt series if a directory is provided", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False, row=4, col=0, rowspan=1, colspan=2)#,mode="easy")
	parser.add_argument("--mdoc", type=str,help="Specify a SerialEM .mdoc file or a folder containing same-named .mdoc files", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False, row=5, col=0, rowspan=1, colspan=2)#,mode="easy")
	parser.add_argument("--npk", type=int,help="Number of landmarks to use (such as gold fiducials). Default is 20.", default=20,guitype='intbox',row=6, col=0, rowspan=1, colspan=1, mode="easy")

#	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=6, col=0, rowspan=1, colspan=2,mode="easy")

	parser.add_argument("--tltax", type=float,help="Angle of the tilt axis. Note the angle stored internally will have an opposite sign. The program will calculate one if this option is not provided.", default=None,guitype='floatbox',row=7, col=1, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--tltkeep", type=float,help="Fraction of tilts to keep in the reconstruction.", default=.9,guitype='floatbox',row=7, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tltrange", type=str,help="Include only tilts between 'START' and 'STOP', i.e. -40.0,40.0. Default behavior is to include all tilts.", default=None)#, guitype='strbox',row=6, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--outsize", type=str,help="Size of output tomograms. choose from 1k, 2k and 4k. default is 1k", default="1k",guitype='combobox',choicelist="('1k', '2k')",row=8, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--niter", type=str,help="Number of iterations for bin8, bin4, bin2 images. Default if 2,1,1,1", default="2,1,1,1",guitype='strbox',row=8, col=1, rowspan=1, colspan=1,mode='easy[2,1,1,1]')
	
	parser.add_argument("--bytile", action="store_true",help="make final tomogram by tiles.. ", default=False, guitype='boolbox',row=9, col=1, rowspan=1, colspan=1,mode="easy[True]")
	
	parser.add_argument("--load", action="store_true",help="load existing tilt parameters.", default=False,guitype='boolbox',row=10, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--notmp", action="store_true",help="Do not write temporary files.", default=False, mode="easy[True]", guitype="boolbox", row=10,col=1, rowspan=1, colspan=1)

	parser.add_argument("--pkkeep", type=float,help="Fraction of landmarks to keep in the tracking.", default=.9,guitype='floatbox',row=11, col=1, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--compressbits", type=int,help="Number of bits of precision in output tomogram with lossless compression. -1 -> uncompressed float", default=-1,guitype='intbox',row=11, col=0, rowspan=1, colspan=1,mode="easy[8]")

	parser.add_argument("--clipz", type=int,help="Z thickness of the final tomogram output. default is -1, (5/16 of tomogram length)", default=-1,guitype='intbox',row=9, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--bxsz", type=int,help="Box size of the particles for tracking. Default is 32. Maybe helpful to use a larger one for fiducial-less cases..", default=32, guitype='intbox',row=6, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--pk_maxval", type=float,help="Maximum Density value of landmarks (n sigma). Default is -5", default=-5.)
	parser.add_argument("--pk_mindist", type=float,help="Minimum distance between landmarks, as fraction of micrograph length. Default is 0.125", default=0.125)

	parser.add_argument("--correctrot", action="store_true",help="correct for global rotation and position sample flat in tomogram.", default=False,guitype='boolbox',row=12, col=0, rowspan=1, colspan=1,mode="easy")

	#parser.add_argument("--normslice", action="store_true",help="normalize each 2D slice.", default=False)
	parser.add_argument("--filterto", type=float,help="filter to abs.", default=-1)
	parser.add_argument("--filterres", type=float,help="filter final tomogram to target resolution (in A). Default is 40", default=40,guitype='floatbox',row=13, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--extrapad", action="store_true",help="Use extra padding for tilted reconstruction. slower and cost more memory, but reduce boundary artifacts when the sample is thick", default=False,guitype='boolbox',row=15, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--moretile", action="store_true",help="Sample more tiles during reconstruction. Slower, but reduce boundary artifacts when the sample is thick", default=False,guitype='boolbox',row=15, col=1, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--rmbeadthr", type=float, help="Density value threshold (of sigma) for removing beads. high contrast objects beyond this value will be removed. default is -1 for not removing. try 10 for removing fiducials", default=-1,guitype='floatbox',row=14, col=1, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--threads", type=int,help="Number of threads", default=12,guitype='intbox',row=12, col=1, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--tmppath", type=str,help="Temporary path", default=None)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)
	parser.add_argument("--noali", action="store_true",help="skip initial alignment", default=False)
	parser.add_argument("--dryrun", action="store_true",help="skip final reconstruction", default=False)
	parser.add_argument("--patchtrack", type=int, help="use patch tracking before landmark based alignment. input 0/1/2 as the number of patch tracking iterations.", default=-1,guitype='intbox',row=16, col=1, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--posz", action="store_true",help="auto positioning along z axis", default=False,guitype='boolbox',row=14, col=0, rowspan=1, colspan=1,mode="easy")
	parser.add_argument("--xdrift", action="store_true",help="apply extra correction for drifting along x axis", default=False,guitype='boolbox',row=13, col=0, rowspan=1, colspan=1,mode="easy")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--reconmode", type=str, help="Intepolation mode for reconstruction. default is trilinear. check e2help.py for details. Not recommended to change.",default="trilinear")
	parser.add_argument("--maxshift", type=float,help="Maximum shift between tilt(/image size). default is 0.2", default=.2)
	parser.add_argument("--highpass", type=int,help="initial highpass filter for alignment in pixels. default if 3", default=3)
	parser.add_argument("--pathtracktile", type=int,help="number of tile along one axis for patch tracking. default is 5", default=5)
	parser.add_argument("--badone", action="store_true",help="Remove one bad tilt during coarse alignment. seem to work better with smaller maxshift...", default=False)#, guitype='boolbox',row=9, col=0, rowspan=1, colspan=1,mode="easy")
	
	parser.add_argument("--flip", action="store_true",help="Flip the tomogram by rotating the tilt axis. need --load existing alignment", default=False)
	parser.add_argument("--autoclipxy", action="store_true",help="Optimize the x-y shape of the tomogram to fit in the tilt images. only works in bytile reconstruction. useful for non square cameras.", default=False,guitype='boolbox',row=16, col=0, rowspan=1, colspan=1,mode="easy")


	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)
	time0=time.time()
	
	#### deal with multiple inputs
	if options.alltiltseries:
		fld="tiltseries/"
		args=[fld+f for f in sorted(os.listdir(fld)) if (
			f.endswith(".hdf") or f.endswith(".mrc") or f.endswith(".mrcs") or f.endswith(".st") or f.endswith(".lst"))]
	
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
		#### note this is actually a bit unsafe here since incomplete option in command will cause infinite loops... hopefully people always spell alltiltseries correctly.
		opt=opt.replace("--alltiltseries","")
		for a in args:
			run("{} {} {}".format(cmd[0], a, opt))
			
		E2end(logid)
		return
		
	#### initializing multithread fft. whether this actually improves speed is untested.
	Util.init_threads(options.threads)
		
	#### parse options.
	options.ctf=None
	
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

	options.basename=bname=base_name(inputname)
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
	
	#### apix of raw input will be used to calculate scaling factor in various functions later
	## the translation values used throughout the program are based on 2k tomograms. they will be scaled in individual function calls based on the apix difference
	options.apix_init=float(imgs[0]["apix_x"])
	
	#### need to make sure this works for images of all sizes (2k, 4k and hopefully 8k)
	imgsz=min(imgs[0]["nx"],imgs[0]["ny"])
	if imgsz<=512:
		print("Tilt series image too small. Only support 2K or larger input images...")
		return
	
	elif imgsz<=1024:
		imgs_1k=imgs_2k=imgs_4k=imgs
		itnum[3]=itnum[2]=0
		options.outsize="1k"
		
	elif imgsz<=2048:
		#### 2k or smaller input. skip 4k refinement
		imgs_2k=imgs_4k=imgs
		itnum[3]=0
	elif imgsz<=4096:
		#### 4k input
		imgs_4k=imgs
		imgs_2k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_4k]
			
	else:
		#### even larger images.. hopefully 8k at most
		bf=imgsz//4096+1
		imgs_4k=[img.process("math.meanshrink", {"n":bf}).process("normalize.edgemean") for img in imgs]
		imgs_2k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_4k]
		imgs=imgs_4k
	
	if imgsz>1024:
		imgs_1k=[img.process("math.meanshrink", {"n":2}).process("normalize.edgemean") for img in imgs_2k]
	
	#### 500x500 images. this is used for the initial coarse alignment so we need to filter it a bit. 
	imgs_500=[]
	for p in imgs_1k:
		m=p.process("math.meanshrink", {"n":2})
		m.process_inplace("filter.ramp")
		m.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.025})
		m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
		m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":options.highpass})
		m.process_inplace("normalize.edgemean")
		sz=512
		m.clip_inplace(Region((m["nx"]-sz)//2, (m["ny"]-sz)//2, sz,sz))
		m.process_inplace("normalize.edgemean")
		m["apix_x"]=m["apix_y"]=p["apix_x"]*2.
		imgs_500.append(m)
		
	imgs_500mask=[m.process("mask.soft",{"outer_radius":-8,"width":4}) for m in imgs_500]

	num=options.num=len(imgs_500)
	boxes3d=[]
	#### create working directory
	## we do not create the folder until here so we do not end up with empty folders when the input is wrong..
	if options.writetmp:
		if options.tmppath:
			path=options.tmppath
		else:
			path=options.tmppath=num_path_new("tomorecon")
		print("Temporary files will be written in {}".format(options.tmppath))
		
	#### save 1k input only.
	if options.writetmp:
		inppath=os.path.join(options.tmppath,"tltseries_input.hdf")
		for i,m in enumerate(imgs_1k):
			if options.compressbits<0: m.write_image(inppath, i)
			else: m.write_compressed(inppath, i, options.compressbits, nooutliers=True)
	
	if options.load:
		#### loading parameters from json file
		jsname=info_name(options.inputname)
		print("Loading parameters from {}...".format(jsname))
		js=dict(js_open_dict(jsname)).copy()
		if "tlt_params" not in js:
			print("Failed to load saved parameterss. Exit.")
			return
		
		if js["tlt_file"]!=options.inputname:
			print("Warning. Tilt file in the saved file is not the same as the input. Something may be wrong..")
			print("\tprevious input :{}\t current input : {}".format(js["tlt_file"], options.inputname))
			
		tpm=np.array(js["tlt_params"])
		ttparams=tpm.copy()
		
		if "defocus" in js:
			print("Loading CTF information. will do phase flipping for tomograms")
			options.ctf={	"defocus":js["defocus"], "phase":js["phase"], 
					"cs":js["cs"], "voltage":js["voltage"]}
			
		tlts=ttparams[:,3].copy()
		#### some old version may not save loss
		try: 
			loss0=np.array(js["ali_loss"])
		except: 
			loss0=abs(ttparams[:,3])*.1
			
		options.zeroid=zeroid=np.argmin(abs(tlts))
		if options.flip:
			print("Flipping tilt axis...")
			ttparams[:,2]=(ttparams[:,2]+180)%360
			
			if "boxes_3d" in js:
				boxes3d=js["boxes_3d"]
				for b in boxes3d:
					for k in range(3):
						b[k]*=-1

			
	else:
		#### determine alignment parameters from scratch
		if (options.rawtlt!=None and len(options.rawtlt)>0) :
			
			if os.path.isdir(options.rawtlt):
				print("Looking for raw tilt file matching {} in {}...".format(options.inputname, options.rawtlt))
				bnm=options.basename.split('/')[-1]
				rl=[f for f in os.listdir(options.rawtlt) if f.startswith(bnm)]
				if len(rl)==0:
					print("Cannot find raw tilt file in the directory. exit...")
					return
				elif len(rl)>1:
					print("Found multiple matching files...")
					print(rl)
					
					
				options.rawtlt=os.path.join(options.rawtlt, rl[0])
				print("Using {}".format(rl[0]))
				
			#### load raw tilt file
			tlts=np.loadtxt(options.rawtlt)
			
			if len(tlts)!=len(imgs_500):
				print("Number of tilt angles and tilt images do no match. Your inputs have {} tilt angles and {} tilt images.\nStopping reconstruction.".format(len(tlts),len(imgs_500)))
				return
			
		elif (options.mdoc!=None and len(options.mdoc)>0) :
			
			if os.path.isdir(options.mdoc):
				print("Looking for mdoc file matching {} in {}...".format(options.inputname, options.mdoc))
				bnm=options.inputname.split('/')[-1].rsplit(".",1)[0]
				rl=[f for f in os.listdir(options.mdoc) if f.startswith(bnm) and f[-5:]==".mdoc"]
				if len(rl)==0:
					print(f"Cannot find mdoc file for {bnm} in the directory {options.mdoc}. exit...")
					return
				elif len(rl)>1:
					print("Found multiple matching files...")
					print(rl)
					
				options.mdoc=os.path.join(options.mdoc, rl[0])
				print("Using {}".format(rl[0]))

			with open(options.mdoc,"r") as mdoc:
				if imgs_1k[0].has_attr("tilt_angle"): print("WARNING: .mdoc specified, but tilt_angle already in header. May lead to incorrect results!")
				tltdic={}
				idx=0
				for line in mdoc:
					if line.startswith("[ZValue"): idx=int(line.split("=")[1].split("]")[0])
					if line.startswith("TiltAngle"): 
						tltdic[idx]=float(line.split("=")[1])
						if tltdic[idx]==0: options.zeroid=idx
			
			try:
				tlts=np.array([tltdic[i] for i in range(len(tltdic))])
			except:
				print(f"Error: incomplete tilt list in mdoc file ",options.mdoc)
				return
			
			if len(tlts)!=len(imgs_500):
				print("Number of tilt angles and tilt images do no match. Your inputs have {} tilt angles and {} tilt images.\nStopping reconstruction.".format(len(tlts),len(imgs_500)))
				return
		else: 
			#### parse tilt step
			options.rawtlt=None
			if options.tltstep>0:
				print("Using fixed tilt step of ",options.tltstep)
				tlts=np.arange(-len(imgs_2k)*options.tltstep/2,len(imgs_2k)*options.tltstep/2,options.tltstep)
			else:
				print("Using tilt_angle from header")
				tlts=np.array([i["tilt_angle"] for i in imgs_1k])
				options.zeroid=np.argmin(abs(tlts))
		
		if options.writetmp: np.savetxt(os.path.join(options.tmppath,"rawtilt.txt"), tlts)
		
		#### we need a zero degree tilt as reference to position the tomogram
		if options.zeroid<0:
			#zeroid=np.argmin(abs(tlts))
			zeroid=len(tlts)//2
			options.zeroid=zeroid
		
		#### here we always assume the center tilt is at 0 degree
		tlts-=tlts[options.zeroid]
		
		ttparams=np.zeros((num, 5))
		if options.tltax!=None:
			### it turns out there is a sign difference between serialem and eman...
			options.tltax=-options.tltax
		
		
		if options.noali:
			print("Skipping coarse alignment...")
		
		else:
			#### do an initial course alignment before common line
			if options.rawtlt!=None or options.mdoc!=None:
				srtid=np.argsort(tlts).tolist()
				imgs_500_sort=[imgs_500[i] for i in srtid]
				ret=calc_global_trans(imgs_500_sort, options)
				pt=np.zeros_like(ret[1])
				pt[srtid]=ret[1]
				ret[1]=pt
				
			else:
				ret=calc_global_trans(imgs_500, options)
				
			#if options.badone:
				#img_tali, pretrans, badi=ret
			#else:
			img_tali, pretrans=ret
			#### estimate initial tilt axis by common line
			if options.tltax==None:
				tltax=calc_tltax_rot(img_tali, options)
				options.tltax=tltax
			
			print("tilt axis:  {:.2f}".format(options.tltax))
			
			if options.writetmp:
				for i,m in enumerate(img_tali):
					if options.compressbits<0: m.write_image(os.path.join(options.tmppath,"tltseries_transali.hdf"), i)
					else: m.write_compressed(os.path.join(options.tmppath,"tltseries_transali.hdf"), i,options.compressbits,nooutliers=True)
			
			#### this is the matrix that save the alignment parameters
			pretrans*=img_tali[0]["apix_x"]/options.apix_init  #### since the pretrans is calculated from 500x500 images..
		
			ttparams[:,0]=-pretrans[:,0] # tx
			ttparams[:,1]=-pretrans[:,1] # ty
			
		ttparams[:,2]=options.tltax # rot
		ttparams[:,3]=tlts.copy() # ytilt
		ttparams[:,4]=0 # off axis tilt
		loss0=abs(ttparams[:,3]) ### this is used to exclude bad tilt. in case the user ask 0 iterations..
		if options.badone:
			loss0[badi]=np.max(loss0)+100
	
	
	if options.patchtrack>0:
		ttparams, loss0= do_patch_tracking(imgs_500, ttparams, options)
		for itr in range(options.patchtrack-1):
			ttparams, loss0= do_patch_tracking_3d(imgs_500, ttparams, options)
		#toitr=[imgs_500, imgs_1k]
		#toitr=toitr[:options.patchtrack]
		#for itr,imgs in enumerate(toitr):
			#ttparams, loss0= do_patch_tracking(imgs, ttparams, options)
			##ttparams, loss0= do_patch_tracking(imgs, ttparams, options)
		if options.writetmp:
			make_ali(imgs_500, ttparams, options, outname=os.path.join(options.tmppath,"ali_patchtrack.hdf"))
			
			
	pks=np.zeros((options.npk, 3))
	#### pack parameters together so it is easier to pass around
	allparams=np.hstack([ttparams.flatten(), pks.flatten()])
	
	#### image scale, m3diter, fidkeep
	#### refinement sequence []:global tilt axis, 0:tx, 1:ty, 2:tilt axis, 3: tilt, 4: off axis tilt
	scaleiter=[(imgs_500, itnum[0], options.pkkeep, [[0,1], [], [0,1],[], [0,1]]),
			   (imgs_1k, itnum[1], options.pkkeep, [[0,1],[], [0,1],[3],[4],[2],[0,1]]),
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
	rawfilterto=options.filterto
	options.filterto=.25
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
			make_ali(imgs_500mask, ttparams, options, outname=name_ali)
		else:
			name_tomo=name_sample=name_ali=name_ptclali=None
			
		
		
		#### make tomogram loop
		for m3diter in range(n_m3d):

			
			#### make initial tomogram for landmark search, always use 500x500
			threed=make_tomogram(imgs_500mask, ttparams, options, outname=name_tomo, errtlt=loss0)

			pks=find_landmark(threed, options)
			allparams=np.hstack([ttparams.flatten(), pks.flatten()])
			
			#### do two round of generating 3D maps of landmarks to refine the location of them
			allparams,smp=make_samples(imgs_, allparams, options, refinepos=True);
			allparams,smp=make_samples(imgs_, allparams, options, refinepos=True);

			if niter==0 and m3diter==0:
				if options.writetmp:
					make_samples(imgs_, allparams, options,
						outname=os.path.join(path,"samples_init.hdf"), refinepos=True);
					ptclpos=ali_ptcls(imgs_, allparams, options, outname=os.path.join(path,"ptclali_init.hdf"), doali=True)
			
				if options.xdrift:
					allparams=np.hstack([ttparams.flatten(), pks.flatten()])
					ttparams=xdrift_correction(imgs_, allparams, options)
					allparams=np.hstack([ttparams.flatten(), pks.flatten()])
					if options.writetmp:
						make_samples(imgs_, allparams, options,
							outname=os.path.join(path,"samples_xdrift.hdf"), refinepos=True);
						ptclpos=ali_ptcls(imgs_, allparams, options, outname=os.path.join(path,"ptcls_xdrift.hdf"), doali=True)
						#make_ali(imgs_1k, ttparams, options, outname=os.path.join(path,"ali_xdrift.hdf"))

			if niter>0:
				#### Now actually refine the alignment parameters using the landmarks.
				for idx in rfseq:
					allparams=refine_one_iter(imgs_, allparams, options, idx)
				
			else:
				allparams=refine_lowres(imgs_, allparams, options)

			ttparams, pks=get_params(allparams, options)
			
			
			
			allparams=np.hstack([ttparams.flatten(), pks.flatten()])
		
		ttparams, pks=get_params(allparams, options)
		if options.writetmp:
			#### make sample output and aligned tilt series
			make_samples(imgs_, allparams, options, outname=name_sample, refinepos=False);
			
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
	
	if options.posz:
		threed=make_tomogram(imgs_500, ttparams, options, errtlt=loss0, clipz=360)
		if options.writetmp:
			if options.compressbits<0: threed.write_image(os.path.join(path,"tomo_posz.hdf"))
			else: threed.write_compressed(os.path.join(path,"tomo_posz.hdf"),0,options.compressbits,nooutliers=True)
			
		ttparams=correct_zpos(threed, ttparams, options)
		
	
	options.filterto=rawfilterto
	#### alignment finish. now save output
	if options.dryrun:
		print("Skipping final tomogram generation...")
		
	else:
		
		if options.outsize=="2k":
			imgout=imgs_2k
		elif options.outsize=="4k":
			options.bytile=True
			imgout=imgs_4k
		elif options.outsize=="500":
			imgout=imgs_500
		else:
			imgout=imgs_1k
		
		if options.filterto<=0:
			options.filterto=imgout[0]["apix_x"]/options.filterres
			print("filter to {} A, {:.2f} of Nyquist.".format(options.filterres, options.filterto))
			
			
		if options.rmbeadthr>0:
			remove_beads(imgs_500, imgout, ttparams, options)
		
		#### only clip z axis at the end..
		if options.bytile:
			threed=make_tomogram_tile(imgout, ttparams, options, errtlt=loss0, clipz=options.clipz)
		else:
			threed=make_tomogram(imgout, ttparams, options, errtlt=loss0, clipz=options.clipz)

		if options.writetmp:
			make_ali(imgout, ttparams, options, outname=os.path.join(path,"tiltseries_ali.hdf"))
			if options.compressbits<0: threed.write_image(os.path.join(path,"tomo_final.hdf"))
			else: threed.write_compressed(os.path.join(path,"tomo_final.hdf"),0,options.compressbits,nooutliers=True)
			
		#### write to the tomogram folder
		try: os.mkdir("tomograms")
		except: pass
		sfx=""
		bf=int(np.round(imgout[0]["apix_x"]/options.apix_init))
		if bf>1:
			sfx+="__bin{:d}".format(int(bf))
			
		tomoname=os.path.join("tomograms", options.basename+sfx+".hdf")
		threed["ytilt"]=yrot
		
		if options.compressbits<0: threed.write_image(tomoname)
		else: threed.write_compressed(tomoname,0,options.compressbits,nooutliers=True)
		print("Tomogram written to {}".format(tomoname))
	
	#### save alignemnt parameters to info file
	tpm=ttparams.copy()
	#tpm[:,:2]*=options.binfac
	js=js_open_dict(info_name(options.basename))
	js["tlt_params"]=tpm.tolist()
	js["tlt_file"]=options.inputname
	js["ali_loss"]=loss0.tolist()
	js["apix_unbin"]=options.apix_init
	if len(boxes3d)>0: js["boxes_3d"]=boxes3d
	js.close()
	
	dtime=time.time()-time0
	print("Finished. Total time: {:.1f}s".format(dtime))

	E2end(logid)


def do_patch_tracking(imgs, ttparams, options, niter=4):
	
	#if imgs0[0]["nx"]==imgs0[0]["ny"]:
		#imgs=imgs0
	#else:
		#s=min(imgs0[0]["nx"], imgs0[0]["ny"])
		#imgs=[m.get_clip(Region((m["nx"]-s)//2, (m["ny"]-s)//2, s,s)) for m in imgs0]
	
	print("\n******************************")
	print("Doing patch tracking")
	scale=int(np.round(imgs[0]["apix_x"]/options.apix_init))
	rawnpk=options.npk
	rawboxsz=options.bxsz
	print("  shrink by {}".format(scale))
	
	ntile=int(np.round(8/scale*2+1))
	tpm=ttparams.copy()
	nx=ny=imgs[0]["nx"]
	nz=len(imgs)
	options.bxsz=sz=256
	
	zeroid=nz//2
	nrange=np.hstack([np.arange(zeroid, nz), np.arange(zeroid, -1, -1)])
	mxsft=64
	
	pks=np.indices((ntile,ntile)).reshape((2,-1)).T-(ntile-1)/2
	pks=np.hstack([pks, np.zeros((len(pks),1))])
	dx=(ny*.45-sz//2)/np.max(pks)
	pks=np.round(pks*dx*scale)
	
	options.npk=len(pks)
	maskc=make_mask(sz)
	imgshp=[m.process("filter.highpass.gauss",{"cutoff_pixels":options.highpass}) for m in imgs]
	
	for itr in range(niter):

		allparams=np.hstack([tpm.flatten(), pks.flatten()])
		ptclpos,ptclimgs=ali_ptcls(imgshp, allparams, options, doali=False,return_imgs=True)
		ptclimgs=np.array([p.numpy().copy() for p in ptclimgs])
		ptclimgs=ptclimgs.reshape((len(pks), -1, sz,sz))

		fts=[]
		ptcl=ptclimgs.reshape((-1, sz,sz))
		for i,p in enumerate(ptcl):
			m=p*maskc
			m=get_fft(m)
			am=np.abs(m)
			am[am==0]=1
			m/=am
			fts.append(m)
			
		fts=np.array(fts)
		fts=fts.reshape(ptclimgs.shape)

		bd=(sz-mxsft)//2
		cccs=np.zeros((len(pks), nz, mxsft, mxsft))
		for fi, ft in enumerate(fts):
			for i in nrange:
				if i==nrange[0]:
					cccs[fi,i,mxsft//2,mxsft//2]=1
					lasti=i
					continue

				a=ft[lasti]
				b=ft[i] 
				c=a*np.conj(b)
				c=get_img(c)
				c=c[bd:-bd,bd:-bd]
				
				cccs[fi, i]=c.copy()
				lasti=i
				
		if itr==1 and scale>6:
			
			#### refine z position of tiles
			
			pz=[]
			zmax=50
			dzrg=np.arange(-zmax,zmax+1)

			for dz in dzrg:
				for nid in range(nz):

					tt=tpm[nid].copy()
					tt[:2]*=0

					pxf=get_xf_pos(tt, [0,0,int(dz)])
					pz.append(pxf)

			pz=np.array(pz).reshape((len(dzrg), nz,-1))
			pz=pz[:,:,::-1]
			
			ps=[]
			for ccc in cccs:
				p=[np.array(np.where(d==np.max(d))).T[0]-mxsft//2 for d in ccc]
				ps.append(p)
				
			ps=np.array(ps)
			
			pc=np.zeros(ps.shape)
			for i in nrange:
				if i==nrange[0]:
					lasti=i
					continue
					
				pc[:,i]=ps[:,i]+pc[:,lasti]
				lasti=i

			ptclz=[]
			for p in pc:
				dp=p-pz
				dp=np.linalg.norm(dp, axis=2)
				dp=np.mean(dp, axis=1)
				ptclz.append(dzrg[np.argmin(dp)])

			ptclz=np.array(ptclz)
			print("	iter {}:  dz = {:.02f}".format(itr, np.mean(abs(ptclz))))
			#print(ptclz)

			pks1=pks.copy()
			pks1[:,2]-=ptclz*scale
			dz=np.mean(ptclz*scale)

			pca=PCA(3)
			pca.fit(pks1);
			c=pca.components_
			t=Transform()
			cc=c[2]
			cc*=np.sign(cc[2])
			t.set_rotation(c[2].tolist())
			t.invert()
			xyz=t.get_params("xyz")
			xyz["ztilt"]=0
			t=Transform(xyz)
			#print(xyz)

			tpm1=tpm.copy()
			for i,tp in enumerate(tpm1):
				xf=Transform({"type":"xyz","ytilt":tp[3],"xtilt":tp[4], "ztilt":tp[2]})
				xf=xf*t.inverse()
				x=xf.get_params("xyz")
				#tpm1[i][2]=x["ztilt"]
				tpm1[i][3]=x["ytilt"]
				tpm1[i][4]=x["xtilt"]
				
				
			dxy=np.array([get_xf_pos(t, [0,0,-dz]) for t in tpm])
			tpm1[:,:2]=dxy
			tpm=tpm1.copy()
			#pks[:,2]-=dz
		
		else:

			ccc=cccs.copy()
			cmax=np.max(ccc, axis=(2,3))
			cmax[cmax==0]=1
			ccc/=cmax[:,:,None,None]
			cc=np.mean(ccc, axis=0)
			px=[np.array(np.where(c==np.max(c))).flatten()-mxsft//2 for c in cc]
			px=np.array(px)

			tilex=np.zeros_like(px)

			for i in nrange:
				if i==nrange[0]:
					lasti=i
					continue
					
				tilex[i]=tilex[lasti]+px[i]
				lasti=i
				
			tx=tilex[:,::-1]
			tpm[:,:2]-=tx*scale
			loss=np.linalg.norm(tx*scale, axis=1)
			print("	iter {}:  loss = {:.2f}".format(itr, np.mean(loss)))

	if options.writetmp:
		allparams=np.hstack([tpm.flatten(), pks.flatten()])
		ptclpos,ptclimgs=ali_ptcls(imgshp, allparams, options, doali=False,return_imgs=True)
		fname=os.path.join(options.tmppath,"patchtrack_ptclali.hdf")
		if os.path.isfile(fname): os.remove(fname)
		for p in ptclimgs:
			p.write_image(fname, -1)
		
	options.npk=rawnpk
	options.bxsz=rawboxsz
	return tpm, loss

def do_patch_tracking_3d_oneiter(imgs, ttparams, options, writetmp=False):
	
	scale=int(np.round(imgs[0]["apix_x"]/options.apix_init))
	
	ntile=options.pathtracktile
	tpm=ttparams.copy()
	nx=ny=imgs[0]["nx"]
	nz=len(imgs)
	options.bxsz=sz=128
	
	mxsft=8
	pks=np.indices((ntile,ntile)).reshape((2,-1)).T-(ntile-1)/2
	pks=pks[:,::-1]
	pks=np.hstack([pks, np.zeros((len(pks),1))])
	dx=(ny*.4-sz//2)/np.max(pks)
	pks=np.round(pks*dx*scale)
	
	options.npk=len(pks)
	maskc=make_mask(sz)
	allparams=np.hstack([tpm.flatten(), pks.flatten()])
	ptclpos,ptclimgs=ali_ptcls(imgs, allparams, options, doali=False,return_imgs=True, do_rotation=False)
	nt=len(ptclimgs)//len(imgs)
	
	trg=tpm[:,3]
	n=len(trg)
	
	nrange=np.arange(len(trg))
	nrange=nrange[np.argsort(abs(np.array(trg)))]
	pad=good_size(sz*1.4)
	dxys=np.zeros((nt, len(trg), 2))
	
	if writetmp:
		pjname=os.path.join(options.tmppath,"patchtrack_projs.hdf")
		if os.path.isfile(pjname): os.remove(pjname)
		tdname=os.path.join(options.tmppath,"patchtrack_tiles.hdf")
		if os.path.isfile(tdname): os.remove(tdname)

	for it in range(nt):
		pimgs=[ptclimgs[i] for i in np.arange(n)+n*it]

		dxy=dxys[it]
		normvol=EMData(pad//2+1, pad, pad)
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"trilinear", "normout":normvol})
		recon.setup()
		# for i,p in enumerate(ptclimgs):
		for i in nrange:
			if i>nrange[0]:
				iix=np.arange(i,len(trg))
			else:
				iix=np.arange(i+1)

			m=pimgs[i].copy()
			m.process_inplace("filter.ramp")
			m.process_inplace("normalize.edgemean")
			m.process_inplace("xform",{"tx":int(dxy[i,0]),"ty":int(dxy[i,1])})
			dy=(sz//2)-np.cos(trg[i]*np.pi/180.)*sz/2
			xf=Transform({"type":"xyz","ztilt":tpm[i,2],"ytilt":tpm[i,3],"xtilt":tpm[i,4]})

			if abs(i-nrange[0])>1:
				pj=recon.projection(xf, 0)
				pj=pj.get_clip(Region(pad//2-sz//2, pad//2-sz//2, sz,sz))

				if writetmp: pj.write_image(pjname,int(n*it+i))

				cf=pj.calc_ccf(m)
				c=cf.calc_max_location_wrap(mxsft, mxsft, 0)
				xf.translate(c)
				dxy[iix]+=[c[0], c[1]]
			
			else:
				if writetmp: m.write_image(pjname,int(n*it+i))

			m=m.get_clip(Region(sz//2-pad//2,sz//2-pad//2, pad, pad), fill=0)
			mp=recon.preprocess_slice(m, xf)
			recon.insert_slice(mp,xf,1)

		if writetmp:
			avg=recon.finish(True)
			avg=avg.get_clip(Region(pad//2-sz//2,pad//2-sz//2,pad//2-sz//2, sz,sz,sz), fill=0)
			avg.write_image(tdname, it)
			
	dxy=np.mean(dxys, axis=0)
	
	return pks,dxys

def do_patch_tracking_3d(imgs, ttparams, options):
	
	print("\n******************************")
	print("Doing patch tracking in 3D")
	rawnpk=options.npk
	rawboxsz=options.bxsz
	scale=int(np.round(imgs[0]["apix_x"]/options.apix_init))
	pmlabel=["trans_x", "trans_y", "tilt_z", "tilt_y", "tilt_x"]

	tpm=ttparams.copy()
	refineseq=[[0,1], [], [0,1], [4], [2], [0,1]]
	for refineid in refineseq:
		
		pks,dxys= do_patch_tracking_3d_oneiter(imgs, tpm, options)

		ps=np.array([ [get_xf_pos(t,p) for t in tpm] for p in pks])
		ddxy=ps-dxys*scale

		allparams=np.hstack([tpm.flatten(), pks.flatten()])
		options.fidkeep=.9
		if len(refineid)==0:
			loss0=global_rot([0],ddxy, allparams, options)
			res=minimize(global_rot, [0], (ddxy, allparams, options) ,
				method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
				
			rt=res.x[0]
			print("  tilt axis rotation {:.02f}, loss {:.2f} -> {:.2f}".format(rt, loss0, res.fun))
			tpm[:,2]+=rt
		
		else:
			trans=[]
			loss=[]
			loss0=[]
			for nid in range(len(tpm)):
				xinit=[0]*len(refineid)
				loss0.append(get_loss_pm(xinit,nid, allparams, options, refineid, ddxy))
				
				res=minimize(get_loss_pm, xinit,(nid, allparams, options, refineid, ddxy),
							method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
				trans.append(res.x)
				loss.append(res.fun)
				

			trans=np.array(trans)
			dt=np.mean(np.linalg.norm(trans, axis=1))
			tpm[:,refineid]+=trans
			print("refine "+','.join([pmlabel[i] for i in refineid]))
			print("   change {:.2f}, loss {:.2f} -> {:.2f}".format(dt, np.mean(loss0), np.mean(loss)))
	
	if options.writetmp:
		do_patch_tracking_3d_oneiter(imgs, tpm, options,writetmp=True)
		
	options.npk=rawnpk
	options.bxsz=rawboxsz
	return tpm, np.array(loss)
	
def correct_zpos(tomo, ttparams, options):
	print("Positioning along z axis")
	binfac=tomo["apix_x"]/options.apix_init
	img=tomo.numpy().copy()
	img-=np.mean(img)
	img/=np.std(img)
	val=np.max(abs(img), axis=(1,2))
	val-=np.min(val)
	val/=np.max(val)
	if options.writetmp:
		np.savetxt(os.path.join(options.tmppath,"zpos.txt"), np.vstack([np.arange(len(val)), val]).T)
	thk=options.clipz*4/binfac
	if options.outsize=="2k":
		thk/=2
	elif options.outsize=="4k":
		thk/=4
	rg=np.arange(0,1,0.02)
	t=[]
	for s in rg:
		l=np.where(val>s)[0]
		t.append([l[0],l[-1]])

	t=np.array(t)
	i=np.argmin(abs((t[:,1]-t[:,0])-thk))
	dz=np.mean(t[i])-tomo["nz"]//2
	
	dz=dz*binfac
	print("  shift z by {}".format(dz))
	dxy=np.array([get_xf_pos(t, [0,0,dz]) for t in ttparams])
	ttparamsnew=ttparams.copy()
	ttparamsnew[:,:2]=dxy
	
	return ttparamsnew

	
def make_mask(ss,rr=-1,sft=[0,0]):
	ix, iy=np.indices((ss,ss))
	maskv=np.exp(-.002*(iy-ss//2)**2)
	maskv/=np.max(maskv)

	r=np.sqrt((ix-ss//2-sft[0])**2+(iy-ss//2-sft[1])**2)
	if rr<0:
		rr=ss*.4
	r=r-rr
	r[r<0]=0
	maskc=np.exp(-1e-2*r**2)
	maskc=maskc/np.max(maskc)
	maskc[maskc>1]=1
	return maskc

def xdrift_correction(imgs, allpms, options):
	print("Performing rotation axis drift correction....")
	zeroid=options.zeroid
	nz=options.num
	tpm, pks=get_params(allpms, options)
	nrange=np.hstack([np.arange(zeroid, nz), np.arange(zeroid, -1, -1)])
	prange=np.arange(len(pks))
	nx=ny=imgs[0]["nx"]
	scale=float(imgs[0]["apix_x"])/options.apix_init
	bx=options.bxsz//2
	print("  itr, landmark dz, image dx:")
	for itr in range(3):
		pp=pks/scale
		ptclss=np.zeros((nz,len(pks), bx*2, bx*2))
		for nid in nrange:
			if nid==nrange[0]: 
				lasts=0
				continue
			ptcls=[]
			for pid in prange:
				tt=tpm[nid].copy()
				tt[:2]/=scale

				pxf=get_xf_pos(tt, pp[pid])
				pxf[0]+=nx//2
				pxf[1]+=ny//2

				xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1]})
				e=imgs[nid].get_rotated_clip(xf,(bx*2,bx*2,1)).process("normalize.edgemean")
				e.rotate(-tt[2],0,0)
				ptcls.append(e.numpy().copy())

			ptcls=np.array(ptcls)
			ptclss[nid]=ptcls
			
		ptclss=np.array(ptclss)

		zrg=10
		rg=np.arange(-zrg,zrg+1)
		pts=[]
		scr=[]
		for dz in rg:
			for nid in range(nz):
				ptcls=ptclss[nid]
				tt=tpm[nid].copy()
				tt[:2]/=scale
				pxf=get_xf_pos(tt, [0,0,int(dz)])
				pts.append(pxf)
				
				msk=make_mask(ptcls.shape[-1],1,[pxf[0], pxf[1]])
				p=ptcls*msk
				mp=np.min(p, axis=(1,2))
				scr.append(mp)

		scr=np.array(scr)#.reshape()

		scr=scr.reshape((len(rg),nz,-1))
		scr=np.transpose(scr, (2,0,1))
		dzs=np.array([rg[np.argmin(np.mean(s,axis=1))] for s in scr])
		pks1=pks.copy()
		pks1[:,2]-=dzs*scale
		
		pp=pks1.copy()/scale
		ptclss=np.zeros((nz,len(pks), bx*2, bx*2))
		for nid in range(nz):
			ptcls=[]
			for pid in prange:
				tt=tpm[nid].copy()
				tt[:2]/=scale

				pxf=get_xf_pos(tt, pp[pid])
				pxf[0]+=nx//2
				pxf[1]+=ny//2

				xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1]})
				e=imgs[nid].get_rotated_clip(xf,(bx*2,bx*2,1)).process("normalize.edgemean")
				e.rotate(-tt[2],0,0)
				ptcls.append(e.numpy().copy())

			ptcls=np.array(ptcls)
			ptclss[nid]=ptcls
			
		ptclss=np.array(ptclss)

		tsx=np.zeros(nz)
		xrg=5

		for nid in nrange:
			if nid==nrange[0]: 
				lasts=0
				continue
			ptcls=ptclss[nid]
			scr=[]
			rg=np.arange(-xrg,xrg+1)
			for dx in rg:
				msk=make_mask(ptcls.shape[-1],1,[0,dx+lasts])
				p=ptcls*msk
				mp=np.min(p, axis=(1,2))
				scr.append(np.mean(mp))

			s=rg[np.argmin(scr)]
			lasts+=s
			tsx[nid]=lasts

		tp=np.zeros((nz,2))
		for i,t in enumerate(tpm):
			tp[i,0]=tsx[i]*np.cos(t[2]*np.pi/180)
			tp[i,1]=tsx[i]*np.sin(t[2]*np.pi/180)	
			
		tpm1=tpm.copy()
		tpm1[:,:2]-=tp*scale
		
		
		print(itr, np.mean(abs(dzs)), np.mean(abs(tsx)))

		tpm=tpm1.copy()
		pks=pks1.copy()


	return tpm


def remove_beads(imgs_500, imgout, ttparams, options):
	
	print("Removing high contrast objects...")
	threed=make_tomogram(imgs_500, ttparams, options, errtlt=options.loss0, doclip=False, padr=1.5)
	threed.process_inplace("math.meanshrink",{"n":2})
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./100})
	threed.process_inplace("filter.highpass.gauss",{"cutoff_freq":1./50})
	threed.process_inplace("normalize")
	if options.writetmp:
		if options.compressbits<0: threed.write_image(os.path.join(options.tmppath, "tomo_rmbead.hdf"))
		else: threed.write_compressed(os.path.join(options.tmppath, "tomo_rmbead.hdf"),0,options.compressbits,nooutliers=False)
	
	vthr=options.rmbeadthr
	img=threed.numpy().T.copy()
	#print(img.shape)
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
	if options.writetmp:
		np.savetxt(os.path.join(options.tmppath, "beadpos.txt"), pts)
		oname=os.path.join(options.tmppath, "samples_beads.hdf")
	else: 
		oname=None
	apm, smp=make_samples(imgs_500, allparams, options, outname=oname)
	
	nx, ny=imgout[0]["nx"],imgout[0]["ny"]
	
	tpms=ttparams.copy()
	scale=imgout[0]["apix_x"]/options.apix_init
	tpms[:,:2]/=scale
	pts/=scale
	
	if options.writetmp:
		fname=os.path.join(options.tmppath, "ptcls_rmbead.hdf")
		if os.path.isfile(fname):
			os.remove(fname)
	pad=good_size(options.bxsz*2)
	for n, tpm in enumerate(tpms):
		
		plst=[]
		for i, pk in enumerate(pts):
			#if smp[i]["minimum"]>-4:
				#print(smp[i]["minimum"])
				#continue
			
			pxf=get_xf_pos(tpm, pk)
			
			
			pxf=[int(round(pxf[0]))+nx//2-pad//2, int(round(pxf[1]))+ny//2-pad//2]
			
			if min(pxf)<-pad/4 or min(nx-pxf[0], ny-pxf[1])<-pad/4:
				continue
			
			if len(plst)>0:
				dst=scidist.cdist(plst, [pxf])
				if np.min(dst)<pad/4:
					continue
			plst.append(pxf)
			e=imgout[n].get_clip(Region(pxf[0], pxf[1], pad, pad))
			e["nid"]=n
			e["pid"]=i
			
			
			m=e.process("normalize.edgemean")
			m.process_inplace("filter.highpass.gauss",{"cutoff_freq":1./500})
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
			
			if m["mean"]>0.4:
				#print(i,m["mean"])
				continue
			
			m.process_inplace("mask.soft",{"outer_radius":int(pad*.25),"width":8})
			
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
	print("Coarse translational alignment...")
	data=[]
	nz=len(imgs)
	sz=256
	for m in imgs:
		s=m.process("math.meanshrink",{"n":2})
		s.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.4})
		s.process_inplace("normalize.edgemean")
		s.clip_inplace(Region((s["nx"]-sz)//2, (s["ny"]-sz)//2, sz,sz))
		data.append(s.numpy().copy())
		
	data=np.array(data)
	data00=data.copy()

	ix, iy=np.indices((sz,sz))

	r=np.sqrt((ix-sz//2)**2+(iy-sz//2)**2)
	r=r-sz*.4
	r[r<0]=0
	maskc=np.exp(-1e-2*r**2)
	maskc=maskc/np.max(maskc)
	maskc[maskc>1]=1
	
	maskd=np.sqrt((ix-sz//2)**2+(iy-sz//2)**2)
	#print(np.min(maskd), np.max(maskd))
	maskd=(maskd<sz*options.maxshift).astype(float)
	
	trans=np.zeros((nz, 2))
	for itr in range(3):
		fts=[]
		for i in range(nz):
			m=data[i]*maskc
			m=get_fft(m)
			am=np.abs(m)
			am[am==0]=1
			m/=am
			fts.append(m)
		
		
		ts=[]
		ts.append([0,0])
		for i in range(1,nz):
			a=fts[i-1]
			b=fts[i]
			
			c=a*np.conj(b)
			c=get_img(c)
			c=c*maskd
			p=np.array(np.where(c==np.max(c))).T[0]-sz//2
			ts.append(p)
		
		ts=-np.cumsum(ts, axis=0)
		ts-=ts[nz//2]
		meants=np.mean(abs(ts), axis=0)
		print(f"  iter {itr}:  tx = {meants[0]:.2f}, ty = {meants[1]:.2f}")
		ts=ts.astype(float)+trans
		tx=-np.round(ts).astype(int)
		
		data=data00.copy()
		for i in range(nz):
			data[i]=np.roll(data00[i], tx[i], axis=[0, 1])
		trans=ts.copy()
		
	trans=-ts[:,::-1].copy()*2
	
	imgout=[]
	for i,m in enumerate(imgs):
		e=m.process("mask.soft",{"outer_radius":-8,"width":8})
		t=trans[i]
		e.translate(t[0], t[1],0)
		imgout.append(e)
		
		
	return [imgout, trans]

#### coarse translational alignment
def calc_global_trans_0(imgs, options, excludes=[], tltax=None,tlts=[]):
	print("Doing coarse translational alignment...")

	num=len(imgs)
	sz=min(imgs[0]["nx"], imgs[0]["ny"])
	
	trans=[[0,0]]
	imgp=[]
	for i,m0 in enumerate(imgs):
		m=m0.process("mask.soft",{"outer_radius":sz*options.maxshift,"width":sz*.1})
		m.process_inplace("filter.lowpass.gauss", {"cutoff_freq":.005})
		imgp.append(m)
		if i==0:
			e0=m
			continue
			
		
		ma=m.align("translational", e0, {"maxshift":sz//2})
		xf=ma["xform.align2d"].get_params("2d")
		
		trans.append([xf["tx"], xf["ty"]])
		e0=m

	trans=np.array(trans)
	
	if options.badone:
		### find the tilt that jumps the most
		trans2=[]
		for i in range(2,len(imgs)):
			ma=imgp[i].align("translational", imgp[i-2], {"maxshift":sz//2})
			xf=ma["xform.align2d"].get_params("2d")
			trans2.append([xf["tx"], xf["ty"]])
		trans2=np.array(trans2)

		loss=[np.linalg.norm(t- (trans[i+1]+ trans[i+2])) for i, t in enumerate(trans2)]
		badi=np.argmax(loss)+1
		trans[badi]=[0,0]
		trans[badi+1]=trans2[badi-1]
		print("Maximum tilt jump at image {}. Removing the tilt.".format(badi))

	tx=np.cumsum(trans, axis=0)
	tx-=np.mean(tx, axis=0)
	
	imgout=[]
	for i,m in enumerate(imgs):
		e=m.process("mask.soft",{"outer_radius":sz*.3,"width":sz//8})
		#e.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.1})
		t=tx[i]
		e.translate(t[0], t[1],0)
		imgout.append(e)
		
	if options.badone:
		return [imgout, tx, badi]
	else:
		return [imgout, tx]
	


#### Find tilt axis by common line. modified from the original e2tomogram. Not very stable but mostly fine.. 
def calc_tltax_rot(imgs, options):
	print("Calculating tilt axis rotation...")
	num=len(imgs)
	nx, ny=imgs[0]["nx"], imgs[0]["ny"]
	sz=min(nx, ny)

	imgnp=[]
	for i in range(num):
		m=imgs[i].get_clip(Region(nx//2-sz//2, ny//2-sz//2, sz, sz))
		
		m.process_inplace("mask.gaussian",{"outer_radius":sz//4})
		m=get_fft(m.numpy().copy())
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
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":options.reconmode})
	recon.setup()

	for i in range(len(imgs)):
		t=tpm[i]
		m=imgs[i]
		if m["nx"]==1:
			continue
		
		if m.has_attr('ctf'):
			ctf=m["ctf"]
			fft1=m.do_fft()
			flipim=fft1.copy()
			ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			m=fft1.do_ift()
			
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
	if options.reconmode=="gauss_2":
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
	time0=time.time()
	num=len(imgs)
	scale=imgs[0]["apix_x"]/options.apix_init
	imgsz=min(imgs[0]["nx"],imgs[0]["ny"])
	if imgsz<=1024*1.1:
		b=1
	elif imgsz<=2048*1.1:
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
	nx, ny=imgs[0]["nx"], imgs[0]["ny"]
	
	if options.autoclipxy:
		bds=[]
		for t in tpm:
			rot=Transform({"type":"2d","tx":t[0], "ty":t[1],"alpha":t[2]})
			p=np.array([rot.transform([nx/2, ny/2, 0]),rot.transform([nx/2, -ny/2, 0])])
			p=np.max(abs(p), axis=0)
			
			bds.append(p)
			
		bds=np.array(bds)
		bds=np.median(abs(bds), axis=0)*2 ## so we clip a rectangle area that covers half of the tilt images
		
		outx=good_size(bds[0])
		outy=good_size(bds[1])
		print("Final tomogram shape: {} x {}".format(outx, outy))
	else:
		outx=outy=1024*b
	
	#outxy=1024*b
	sz=max(int(clipz*0.8), 256*b) #### this is the output 3D size 
	step=sz//2 #### distance between each tile
	if options.extrapad:
		pad=good_size(sz*2) #### this is the padded size in fourier space
	else:
		pad=good_size(sz*1.4) #### this is the padded size in fourier space
	
	if clipz>0:
		outz=clipz
	else:
		outz=sz#good_boxsize(sz*1.2)

	#### we make 2 tomograms with half a box shift and average them together to compensate for boundary artifacts.
	
	
	#options.moretile=True
	if options.moretile:
		full3d=EMData(outx, outy, outz)
		mem=(outx*outy*outz*4+pad*pad*pad*options.threads*4)
		print("This will take {}x{}x{}x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outx, outy, outz, pad, pad, pad,options.threads, mem/1024**3))
		wtcon=1
	else:
		full3d=[EMData(outx, outy, outz), EMData(outx, outy, outz)]
		mem=(outx*outy*outz*2*4+pad*pad*pad*options.threads*4)
		print("This will take {}x{}x{}x2x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outx, outy, outz, pad, pad, pad,options.threads, mem/1024**3))
		wtcon=2.5
	
	
	jsd=queue.Queue(0)
	jobs=[]
	nstepx=int(outx/step/2)
	nstepy=int(outy/step/2)
	
	if options.ctf!=None:
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":options.ctf["voltage"], "bfactor":0., "cs":options.ctf["cs"],"ampcont":0, "apix":imgs[0]["apix_x"]})
		dfs=[]
		
	for stepx in range(-nstepx,nstepx+1):
		#### shift y by half a tile
		if options.moretile:
			yrange=range(-nstepy,nstepy+1)
		else: 
			yrange=range(-nstepy+stepx%2,nstepy+1,2)
		for stepy in yrange:
			tiles=[]
			for i in range(num):
				if i in nrange:
					t=tpm[i]
					pos=[stepx*step,stepy*step,0]
					pxf=get_xf_pos(t, pos)
					img=imgs[i]
					m=img.get_clip(Region(img["nx"]//2-pad//2+pxf[0],img["ny"]//2-pad//2+pxf[1], pad, pad), fill=0)
					if options.ctf!=None:
						rot=Transform({"type":"xyz","xtilt":float(t[4]),"ytilt":float(t[3])})
						p1=rot.transform(pos)
						pz=p1[2]*img["apix_x"]/10000.
						ctf.defocus=options.ctf["defocus"][i]-pz
						ctf.set_phase(options.ctf["phase"][i]*np.pi/180.)
						dfs.append(ctf.defocus)
						m["ctf"]=ctf
						
					tiles.append(m)
				else:
					tiles.append(EMData(1,1))

			jobs.append((jsd, tiles, tpm, sz, pad, stepx, stepy, outz, options))
	
	if options.ctf!=None:
		print("Doing Ctf correction. Average defocus {:.2f}".format(np.mean(dfs)))
		
		
	thrds=[threading.Thread(target=make_tile,args=([i])) for i in jobs]
	print("now start threads...")
	thrtolaunch=0
	tsleep=threading.active_count()
	
	#### non-round fall off. this is mathematically correct but seem to have grid artifacts
	#f=np.zeros((sz,sz))
	x,y=np.indices((sz,sz),dtype=float)/sz-.5
	#f=.25-(x**2+y**2)/2 + ((abs(x)-0.5)**2+(abs(y)-0.5)**2)/2
	f=wtcon+np.exp(-(x**2+y**2)/0.1) - np.exp(-((abs(x)-0.5)**2+(abs(y)-0.5)**2)/0.1)
	f3=np.repeat(f[None, :,:], outz, axis=0)
	msk=from_numpy(f3).copy()
	#####

	
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep or not jsd.empty():
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(.1)
		
		if not jsd.empty():
			stepx, stepy, threed=jsd.get()
			#threed["pos"]=[stepx, stepy]
			#threed.write_image("alltiles.hdf", -1)
			threed.mult(msk)
			#### insert the cubes to corresponding tomograms
			if options.moretile:
				full3d.insert_scaled_sum(
					threed,(int(stepx*step+outx//2),int(stepy*step+outy//2), outz//2))
			else:
				full3d[stepx%2].insert_clip(
				threed,
				(int(stepx*step+outx//2-sz//2),
				int(stepy*step+outy//2-sz//2), 
				outz//2-threed["nz"]//2))
				
				
	for t in thrds: t.join()
	
	if not options.moretile:
		full3d=full3d[0]+full3d[1]
	full3d.process_inplace("normalize")
	
	#### skip the tomogram positioning step because there is some contrast difference at the boundary that sometimes breaks the algorithm...
	full3d["zshift"]=0
	
	apix=imgs[0]["apix_x"]
	full3d["apix_x"]=full3d["apix_y"]=full3d["apix_z"]=apix
	#if options.normslice:
		#full3d.process_inplace("normalize.rows")
		
	print("Reconstruction done ({:.1f} s). Now writting tomogram to disk...".format(time.time()-time0))
	return full3d

#### reconstruct tomogram...
def make_tomogram(imgs, tltpm, options, outname=None, padr=1.2,  errtlt=[], clipz=-1, doclip=True):
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
	outxy=good_size(min(nx, ny))

	pad=good_size(outxy*padr)
	#############
	#clipz=options.clipz
	zthick=good_size(max(clipz*1.2, pad//2))
	if options.verbose:
		print("\t Image size: {:d} x {:d}".format(nx, ny))
		print("\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick))
		
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":options.reconmode})
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
	if doclip:
		if clipz<0:
			threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), 0, outxy, outxy, zthick))
		else:
			threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), (zthick-clipz)//2, outxy, outxy, clipz))
	threed["zshift"]=0

	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix

	#if options.normslice and outxy>1000:
		#threed.process_inplace("normalize.rows")
	if outname:
		if options.compressbits<0: threed.write_image(outname)
		else: threed.write_compressed(outname,0,options.compressbits,nooutliers=True)
		if options.verbose: print("Map written to {}.".format(outname))

	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform,  exclude, options):
	m=img.copy()
	#### the ramp filter and decay edge helps soften the edge artifacts
	m.process_inplace("filter.ramp")
	m.process_inplace("normalize")
	m.process_inplace("mask.decayedge2d", {"width":int(pad//40)})
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
		m=im.process("normalize.edgemean")
		p2=m.get_clip(Region(old_div(m["nx"],2)-old_div(pad,2),old_div(m["ny"],2)-old_div(pad,2), pad, pad), fill=0)
		po=p2.copy()
		po.translate(-pxf[0], -pxf[1], 0)
		po.rotate(-tpm[2],0,0)
		xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4]})
		po["xform.projection"]=xform
		if options.compressbits<0: po.write_image(outname, nid)
		else: po.write_compressed(outname,nid,options.compressbits,nooutliers=True)

#### search for alignment landmarks in the tomogram
def find_landmark(threed, options):
	print("Searching for landmarks...")
	#threed0=threed.process("normalize.rows")
	#### use minshrink so we keep the minimas
	threedtiny=threed.process("math.minshrink", {"n":2})
	threedtiny.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
	threedtiny.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
	threedtiny.process_inplace("normalize")
	mapnp=threedtiny.numpy().copy()
	#asrt= np.argsort(mapnp.flatten())
	#b=10
	#m=mapnp.copy()
	#m[:,:b]=m[:,-b:]=m[:,:,:b]=m[:,:,-b:]=0
	asrt= np.argsort(mapnp.flatten())

	#### go through every point starting from the darkest ones.
	### there should be a much faster way but we are only doing it on ~250 cubes so it is not too slow
	pts=[]
	imgs=[]
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
		#m=threedtiny.get_clip(Region(int(pt[2])-b ,int(pt[1])-b ,int(pt[0]) ,b*2, b*2, 1))
		#imgs.append(m)
		
		if mapnp[pt]>vthr:
			break
		
		#### generate a bit more landmarks for some randomness
		if len(pts)>=options.npk*1.5:
			break

	if len(pts)<options.npk:
		print("Found only {} landmarks".format(len(pts)))
		options.npk=len(pts)
	else:
		#### I commented out the shuffling for stable testing. maybe should add it back sometime..
		#np.random.shuffle(pts)
		#area=[]
		#for i,m in enumerate(imgs):
			#c=np.array(m.calc_radial_dist(32,0,1,1))
			#c/=c[0]
			#c=np.mean(c[3:])-np.mean(c[:3])
			#area.append(c)
		#aid=np.argsort(area)
		#pts=[pts[a] for a in aid]

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
def get_center(img, lowres=True, maxshift=64):
	e=img.copy()
	bx=e["nx"]/2
	e.mult(-1)
	#### low resolution mode: use peak
	if lowres:
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.01})
		e.process_inplace("mask.gaussian",{"outer_radius":bx*.5})
		if e["sigma"]==0:
			pk=[0,0]
		else:
			pk=bx-np.array(e.calc_max_location())

	#### high resolution mode: use center of mass
	else:
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.01})
		e.process_inplace("mask.gaussian",{"outer_radius":bx*.5})
		e.process_inplace("normalize")
		pk=bx-np.array(e.calc_center_of_mass(2))[:2]
		if np.isnan(pk).any() or np.linalg.norm(pk)>maxshift:
			### give up if there is too much drift for center of mass...
			pk=[0,0]


	return pk

#### this is the main function that calculate the offset of landmarks in tilt series given images and alignment parameters 
def ali_ptcls(imgs, allpms, options, outname=None, doali=True, return_imgs=False, do_rotation=False):
	zeroid=options.zeroid
	num=options.num
	nrange=np.hstack([np.arange(zeroid, num), np.arange(zeroid, -1, -1)])
	ttparams, pks=get_params(allpms, options)
	scale=imgs[0]["apix_x"]/options.apix_init
	ttparams[:,:2]/=scale
	pks/=scale
	prange=np.arange(options.npk)
	if outname:
		try:os.remove(outname)
		except:pass

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
	ptclimgs=[0]*(options.npk*num)
	for pid in prange:
		#### loop through all particles
		trans=np.zeros((num,2))
		fid=[]
		ppos=np.zeros((num,2))
		for ii,nid in enumerate(nrange):
			#### start from center tilt and go both directions
			pxf=get_xf_pos(ttparams[nid], pks[pid])
			pxf[0]+=nx/2
			pxf[1]+=ny/2

			if nid!=zeroid:
				tlast=trans[nrange[ii-1]]
				pxf[0]-=tlast[0]
				pxf[1]-=tlast[1]
			else:
				tlast=np.array([0,0])

			if do_rotation:
				xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1],"alpha":-ttparams[nid, 2]})
			else:
				xf=Transform({"type":"2d","tx":pxf[0],"ty":pxf[1]})
			
			e=imgs[nid].get_rotated_clip(xf,(bx*2,bx*2,1)).process("normalize")
			e["apix_x"]=e["apix_y"]=e["apix_z"]=imgs[0]["apix_x"]

			ts=[0,0]
			if doali:# and nid!=zeroid:
				#### still get the center of particle differently at low and high res mode
				tx=get_center(e, lowres, maxshift=8)
				trans[nid]=tlast+np.array([tx[0], tx[1]])
			else:
				tx=[0,0]

			if outname:
				e["score"]=[float(tx[0]+tlast[0]), float(tx[1]+tlast[1])]#trans[nid].tolist()
				e["pid"]=int(pid)
				e["nid"]=int(nid)
				e.write_image(outname, int(nid+pid*num))
				
			if return_imgs:
				ptclimgs[int(nid+pid*num)]=e
			#ppos[nid]=np.array(pxf-trans[nid]-[nx/2, ny/2])
			ppos[nid]=np.array(get_xf_pos(ttparams[nid], pks[pid]))-tlast-[tx[0], tx[1]]

		ptclpos.append(ppos)

	ptclpos=np.array(ptclpos)*scale
	
	if return_imgs:
		return ptclpos, ptclimgs
	else:
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
		res=minimize(global_rot, [0], (ptclpos, apms, options), method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
		res=res.x*1
		try: res=res[0]
		except: pass
		print("refine global tilt_z {:.2f}, loss {:.2f} -> {:.2f}".format(
			float(res), float(global_rot([0],ptclpos,apms, options)),
			float(global_rot([res],ptclpos,apms, options))))
		ttparams[:,2]+=res
		pkc=pks.copy()
		r=res*np.pi/180.
		rotmat=np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
		pkc[:,:2]=np.dot(pkc[:,:2],rotmat).reshape((-1,2))
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
			res=res.x*1
			ttpm_new[nid, idx]+=res
			loss[nid]=get_loss_pm(res, nid, apms, options, idx, ptclpos, False)
			loss0[nid]=get_loss_pm(tminit, nid, apms, options, idx, ptclpos, False)

		ipm=[pmlabel[i] for i in idx]
		print("refine {}, loss: {:.2f} -> {:.2f}".format(ipm, np.mean(loss0), np.mean(loss)))
		ttparams=ttpm_new

	allparams=np.hstack([ttparams.flatten(), pks.flatten()])
	return allparams

#### function for scipy optimizer
def get_loss_pm(pm, nid, apms, options, idx=[2,3,4], ptclpos=[], refine=True):
	ttparams, pks=get_params(apms, options)
	tpm=ttparams[nid].copy()
	tpm[idx]+=pm
	ps=np.array([get_xf_pos(tpm, p) for p in pks])
	dst=np.sqrt(np.sum((ps-ptclpos[:,nid])**2, axis=1))
	dst=np.mean(dst[np.argsort(dst)[:int(len(dst)*options.fidkeep)]])
	dst=dst*options.apix_init/10.

	#### return distance in nm
	return dst

def test_pkpos(pos, pid, ptclpos, apms, options):
	
	ttparams, pks=get_params(apms, options)
	pt=pks[pid]+pos
	ppos=ptclpos[pid]

	dst=[]

	for nid, t in enumerate(ttparams):
		ps=get_xf_pos(t, pt)
		d=np.linalg.norm(ps-ppos[nid])
		dst.append(d)

	return np.mean(dst)

#### function for scipy optimizer for global rotation
def global_rot_old(rt, ptclpos, allpms, options):
	try: rt=rt[0]
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

def global_rot(rt, ptclpos, allpms, options):
	tpm, pks=get_params(allpms, options)
	rt=rt[0]
	tpm1=tpm.copy()
	tpm1[:,2]+=rt
	pkc=pks.copy()
	r=-rt*np.pi/180.
	rotmat=np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
	pkc[:,:2]=np.dot(pkc[:,:2],rotmat)


	dst=[]

	for nid, t in enumerate(tpm1):
		ps=np.array([get_xf_pos(t, p) for p in pkc])
		d=np.linalg.norm(ps-ptclpos[:,nid], axis=1)
		dst.append(np.mean(d))

	return np.mean(dst)*options.apix_init/10.
	

def refine_lowres(imgs, allparams, options):
	
	tpm, pks=get_params(allparams, options)
	
	for itr in range(3):
		print(f"iter {itr}:")
	
		ptclpos=ali_ptcls(imgs, allparams, options, doali=True)
		#make_samples(imgs, allparams, options, outname=options.tmppath+f"/smp_{itr}_0.hdf", refinepos=False);
		trans=[]
		loss=[]
		for nid in range(len(imgs)):
		
			res=minimize(get_loss_pm, [0,0],(nid, allparams, options, [0,1], ptclpos) ,
				method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
			trans.append(res.x)
			loss.append(res.fun)
		
		trans=np.array(trans)
		dt=np.mean(np.linalg.norm(trans, axis=1))
		tpm[:,:2]+=trans
		print("  trans = {:.2f}, loss = {:.2f}".format(dt, np.mean(loss)))
		
		allparams=np.hstack([tpm.flatten(), pks.flatten()])
		#make_samples(imgs, allparams, options, outname=options.tmppath+f"/smp_{itr}_1.hdf", refinepos=False);
		ptclpos=ali_ptcls(imgs, allparams, options, doali=True)
		ptrans=[]
		for pid, p in enumerate(pks):
			res=minimize(test_pkpos, [0,0,0],(pid, ptclpos, allparams, options) ,
				method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
			
			ptrans.append(res.x)
		
		ptrans=np.array(ptrans)
		print("  landmark shift = {:.2f}".format(np.mean(np.linalg.norm(ptrans, axis=1))))
		
		pks+=ptrans
		allparams=np.hstack([tpm.flatten(), pks.flatten()])		
		#make_samples(imgs, allparams, options, outname=options.tmppath+f"/smp_{itr}_2.hdf", refinepos=False);
		ptclpos=ali_ptcls(imgs, allparams, options, doali=True)
		
		
		res=minimize(global_rot, [0], (ptclpos, allparams, options) ,
			method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
		
		rt=res.x
		try:rt=rt[0]
		except:pass
	
		print("  tilt axis rotation = {:.02f}, loss = {:.02f}".format(rt, res.fun))
		tpm[:,2]+=rt
		r=-rt*np.pi/180.
		rotmat=np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
		pks[:,:2]=np.dot(pks[:,:2],rotmat)

		allparams=np.hstack([tpm.flatten(), pks.flatten()])
		#make_samples(imgs, allparams, options, outname=options.tmppath+f"/smp_{itr}_3.hdf", refinepos=False);
		
	return allparams


def run(cmd):
	print(cmd)
	launch_childprocess(cmd)


if __name__ == '__main__':
	main()
