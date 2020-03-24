#!/usr/bin/env python
# Muyuan Chen 2017-03
from past.utils import old_div
from builtins import range
import numpy as np
from EMAN2 import *
import json
from scipy.signal import argrelextrema
from scipy.optimize import minimize
import threading

def calc_ctf(defocus, bxsz=256, voltage=300, cs=4.7, apix=1. ,phase=0.):
	ds=1.0/(apix*bxsz)
	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":0.,"cs":cs,"ampcont":0,"apix":apix,"dfdiff":0,"defang":0})
	ctf.set_phase(phase*np.pi/180.)
	cv=[]
	try: l=len(defocus)
	except: defocus=[defocus]
	
	for df in defocus:
		ctf.defocus=df
		cv.append(abs(np.array(ctf.compute_1d(bxsz,ds,Ctf.CtfType.CTF_AMP))))
	return np.array(cv)


def get_xf_pos(tpm, pk):
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0], pk[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)
	
	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0], p1[1]])
	

	return [p2[0], p2[1], p1[2]]


#### Calculate matching score of a set of power spectrum curves (curve) and a set of 1D CTF curves (allctf)
#### and return the a matrix in which each row corresponds to a ctf curve and each column represent a poser spectrum
#### It does the background subtraction in the same way as e2ctf
def calc_all_scr(curve, allctf, zlist, bxsz, exclude=[]):
	#print curve
	allscr=np.zeros((len(allctf), len(curve)))+np.inf

	for i,cf in enumerate(allctf):
		if i in exclude: continue
		zz=zlist[i]
		if len(zz)==0:
			allscr[i]=1
			continue

		z0=zz[0]
		z1=np.minimum(zz[-1], int(bxsz/2*.5))

		if z1-z0<10: continue

		bg=np.array([np.interp(np.arange(z0, z1), zz, p[zz]) for p in curve])

		if len(curve.shape)==1:
			continue
		
		bsub=curve.copy()
		bsub[:, z0:z1]-=bg
		#for iz in range(1, len(zz)):
			#if zz[iz-1]>=z1: break
			#c=bsub[:,zz[iz-1]:zz[iz]]
			#mx=np.max(c, axis=1)
			#m0=(mx<=0)
			#mx[m0]=1
			#mx=1./mx
			#mx[m0]=0
			#c*=mx[:, None]

		bsub=bsub[:, z0:z1]
		#bsub[bsub<0]*=0.1
		scr=-np.dot(bsub,cf[z0:z1])/(np.sum(cf[z0:z1]))
#		 scr=np.mean((bsub-cf[z0:z1])**2, axis=1)
		allscr[i]=scr
	return allscr

def compute_score(pm, ps, options, sign=-1, fitting=True):
	df0, phase=pm
	allrd, pzus=ps
	df0=max(0,df0)
	phase=max(0,phase)
	apix=options.apix
	bxsz=options.tilesize
	
	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":options.voltage,
				   "bfactor":0.,"cs":options.cs,"ampcont":0,"apix":apix,"dfdiff":0,"defang":0})
		
	ds=1.0/(apix*bxsz)
	allscr=[]
	allcurves=[]
	
	for ii,curve in enumerate(allrd):
		pz=pzus[ii]
		
		ctf.set_phase(phase*np.pi/180.)
		ctf.defocus=df0 + sign*pz
		cf=abs(np.array(ctf.compute_1d(bxsz,ds,Ctf.CtfType.CTF_AMP)))
#		 zz=argrelextrema(cf, np.less)[0]
		zzf=np.array([ctf.zero(i)/ds for i in range(len(cf)/2)])
		zzf=zzf[zzf<bxsz/2]
		cv=np.interp(zzf, np.arange(len(curve)), curve)
	
		zz=np.round(zzf).astype(int)
		if len(zz)==0:  continue
		z0=zz[0]
		z1=zz[zz>bxsz/2*.8]
		if len(z1)>0:
			z1=z1[0]
		else:
			z1=zz[-1]

		z1=min(z1, zz[np.append(True, np.diff(zz)>2)][-1])
		if z1-z0<10: continue
		
		bg=np.interp(np.arange(bxsz//2), zzf, cv)

		bsub=curve.copy()
		bsub-=bg
		

		if fitting:
			#bsub[bsub<0]*=0.1
			#for iz in range(1, len(zz)):
				#if zz[iz-1]>=z1: break
				#c=bsub[zz[iz-1]:zz[iz]]
##				 if len(c)==0: continue
				#mx=np.max(c)
				#if mx>0:
					#c/=mx/np.max(cf[zz[iz-1]:zz[iz]])
			bsub=bsub[z0:z1]
			scr=-np.dot(bsub,cf[z0:z1])/(np.sum(cf[z0:z1]))
			allscr.append(scr)
		else:
			ev=argrelextrema(cf[:z1], np.greater)[0][1:]
			y=np.array(bsub[ev])
			y[y<0]=np.min(y[y>0])
			y=(-np.log(np.sqrt(y))*4)/(srg[ev]**2)
			bf=np.mean(y)
			ctf.bfactor=bf
			cf=abs(np.array(ctf.compute_1d(bxsz,ds,Ctf.CtfType.CTF_AMP)))
			allcurves.append((bsub, cf, curve))
		
		

	if fitting:
		if len(allscr)==0:
			return 10
		else:
			return np.mean(allscr)
	else:
		return allcurves

def main():
	
	usage="""prog  <tilt series 1> <tilt seires 2> ... [options]
	CTF estimation for tilt series. The tomogram reconstruction needs to be done with e2tomogram.py. This will only estimate CTF parameters and write the information in the corresponding info files. Actual per particle CTF correction will be done at the particle extraction step."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tiltseries",help="Specify tiltseries you want to apply CTF correction.", default="", guitype='filebox', browser="EMTiltseriesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--alltiltseries", action="store_true",help="Use all tilt series in the folder. Acceptable file extensions include hdf, mrc, mrcs, st.", default=False,guitype='boolbox',row=1, col=0, rowspan=1, colspan=1,mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1, mode="model")
	parser.add_argument("--dfrange", type=str,help="Search range of defocus (start, end, step). default is 2., 7, 0.02", default="2.0,7.0,0.02", guitype='strbox',row=4, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--psrange", type=str,help="phase shift range (start, end, step). default is 10, 15, 5", default="10,15,5", guitype='strbox',row=4, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--tilesize", type=int,help="Size of tile to calculate FFT, default is 256", default=256, guitype='intbox',row=4, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--voltage", type=int,help="Voltage of microscope in kV", default=300, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--cs", type=float,help="Cs of microscope", default=2.7, guitype='floatbox',row=6, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--nref", type=int,help="Using N tilt images near the center tilt to estimate the range of defocus for all images. Default is 15", default=15, guitype='intbox',row=6, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--stepx", type=int,help="Number of tiles to generate on x-axis (different defocus)", default=20, guitype='intbox',row=8, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--stepy", type=int,help="Number of tiles to generate on y-axis (same defocus)", default=40, guitype='intbox',row=8, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--refine", action="store_true", help="Include a refinement step in the end for more precise estimation.", default=False, guitype='boolbox',row=9, col=0, rowspan=1, colspan=1,mode="model")
	parser.add_argument("--checkhand", action="store_true", help="Check the handedness of tomogram.", default=False,guitype='boolbox',row=10, col=0, rowspan=1, colspan=1,mode="model")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on the local computer",guitype='intbox', row=8, col=2, rowspan=1, colspan=1, mode='model')
	parser.add_argument("--nolog",action="store_true",default=False,help="Default=False. Turn off recording of the command ran for this program onto the .eman2log.txt file")	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=1, help="verbose level [0-9], higher number means higher level of verboseness")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	
	(options, args) = parser.parse_args()
	
	#### deal with multiple inputs
	if options.alltiltseries:
		fld="tiltseries/"
		args=[fld+f for f in os.listdir(fld) if (
			f.endswith(".hdf") or f.endswith(".mrc") or f.endswith(".mrcs") or f.endswith(".lst") or f.endswith(".st"))]
	
	if len(args)==1:
		print("Reading tilt series {}...".format(args[0]))
	else:
		print("Processing {} tilt series in sequence..".format(len(args)))
		if not options.nolog: logid=E2init(sys.argv)
		thrds=["{} {} --nolog --verbose 0 {}".format(parser.prog,a,commandoptions(options,("threads","alltiltseries","verbose"))) for a in sorted(args)]

		NTHREADS=max(options.threads+1,2)	# the controlling thread isn't really doing anything
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==NTHREADS ) : time.sleep(.1)
				if options.verbose>1 : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
				print("running: ",thrds[thrtolaunch])
				thrds[thrtolaunch]=threading.Thread(target=run,args=(thrds[thrtolaunch],))
				thrds[thrtolaunch].start()
				time.sleep(1)				# this helps stagger launches due to the large read at the beginning of each job
				thrtolaunch+=1
			else: time.sleep(1)
			if options.verbose>1 and thrtolaunch>0:
				frac=thrtolaunch/float(len(thrds))
				print("{}% complete".format(100.0*frac))

		for t in thrds:
			t.join()

		if options.verbose : print("All threads complete")
		if not options.nolog: E2end(logid)
		return
	
	if not options.nolog: logid=E2init(sys.argv)

	tfile=args[0]
	try:
		js=js_open_dict(info_name(tfile))
		tltparams=np.array(js["tlt_params"])
		js=None
	except:
		print("Cannot find alignment parameters. Please make sure tomogram reconstruction is done for selected tilt series")
		return

  
	img=EMData(tfile,0)
	apix=img["apix_x"]
	options.apix=apix
	upix=apix/10./1000.
	nx, ny=img["nx"], img["ny"]
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	for m in imgs: m.process_inplace("normalize")
	nz=len(imgs)
	
	box=options.tilesize
	sz=min(nx, ny)-box
	if options.verbose>0:print("Reading {}, size {} x {} x {}, apix {:.2f}".format(tfile, nx, ny, nz, apix))


	drg=[float(i) for i in options.dfrange.split(',')]
	defrg=np.arange(drg[0],drg[1],drg[2])
	defstep=drg[2]
	
	prg=[float(i) for i in options.psrange.split(',')]
	pshift=np.arange(prg[0], prg[1], prg[2])
	
	if options.verbose>0: print("Search defocus from {:.1f} to {:.1f}, step {:.2f}".format(drg[0], drg[1], drg[2]))
	if options.verbose>0: print("Search phase shift from {:.1f} to {QA	:.1f}, step {:.1f}".format(prg[0], prg[1], prg[2]))
	
	allctf=[]
	zlist=[]
	for ps in pshift:
		ctf=calc_ctf(defrg, box, voltage=options.voltage, cs=options.cs, apix=apix, phase=ps)
		zeros=np.array(argrelextrema(ctf, np.less, axis=1))
		zl=[zeros[1, zeros[0,:]==i].astype(int) for i in range(len(ctf))]
		allctf.append(ctf)
		zlist.append(zl)


	npt=options.stepy
	nstep=options.stepx
	powerspecs=[]
	if options.verbose>0:print("Generating power spectrum of tiles from tilt series...")
	for imgi in range(nz):
		rawimg=imgs[imgi]
		tpm=tltparams[imgi].copy()
		xrg=sz/2/np.cos(tpm[3]/180.*np.pi)*.9
		allrd=[]
		pzus=[]
		xstep=xrg/float(nstep)
		for xpos in np.arange(-xrg, xrg+1, xstep):

			pts=np.zeros((npt, 3))
			pts[:,0]=xpos
			pts[:,1]=(np.random.rand(npt)-.5)*sz*.9

			ptsxf=np.array([get_xf_pos(tpm, p) for p in pts])
			pz=ptsxf[:,2].copy()
			ptsxf=ptsxf[:,:2]
			ptsxf+=[nx//2, ny//2]
			
			ptsxf=ptsxf[np.sum((ptsxf>[nx, ny]) + (ptsxf<0), axis=1)==0, :]
			if len(ptsxf)<npt*.8:
				continue
			
			rds=[]
			for p in ptsxf:
				tile=rawimg.get_clip(Region(p[0]-box//2, p[1]-box//2, box, box))
				if tile["sigma"]<1e-3 or tile["mean"] != tile["mean_nonzero"]:  # strong criteria to exclude edges
					#rds.append(np.zeros(box//2))
					continue
				tile.do_fft_inplace()
				rd=np.array(tile.calc_radial_dist(box//2, 0,1,0))
				rd=np.log10(rd)
				rd-=np.min(rd)
				rd/=np.max(rd[1:])
				rd[rd>1]=1
				rds.append(rd)
			
			if len(rds) == 0: rd = np.zeros(box//2)
			else: rd=np.mean(rds, axis=0)
			allrd.append(rd)
			pzus.append(np.mean(pz*upix))

		allrd=np.array(allrd)
		powerspecs.append([allrd, pzus])
		
		if options.verbose>1: print("{}, {:.2f}".format(imgi, tpm[3]))
	
	dfs=[]
	exclude=[]
	tltsrt=np.argsort(abs(tltparams[:,3]))	
	if options.checkhand:
		print("Checking handedness of the tomogram. Will NOT write metadata output...")
		rot=tltparams[nz//2,2]
		print("Current tilt axis rotation {:.2f}".format(-rot))
		signs=[-1, 1]
		scores=[[], []]
		dfs=[[], []]
		print("Comparing current hand vs flipped hand..")
		for it in tltsrt:
			
			allrd, pzus=powerspecs[it]
			for si in [0,1]:
				allscr=[]
				for ic, ctf in enumerate(allctf):
					scr=calc_all_scr(allrd, ctf, zlist[ic], box, exclude)
					idxsft=np.round(signs[si]*np.array(pzus)/defstep).astype(int)
					stilt=np.zeros(len(defrg))+np.inf
					for i,df in enumerate(defrg):
						idx=idxsft+i
						outb=(idx<0)+(idx>=len(defrg))
						idx[outb]=0
				
						s=scr[idx, np.arange(scr.shape[1])].copy()
						sinf=np.isinf(s)
						s[outb]=np.max(s)
						sval=s[np.isinf(s)==0]

						stilt[i]=np.sum(sval)/len(s)
					
					allscr.append(stilt)
					
				allscr=np.array(allscr)
				amp, df= np.array(np.where(allscr==np.min(allscr))).T[0]
				pm=(defrg[df], pshift[amp])
				
				scores[si].append((np.mean(allscr)-np.min(allscr))*100)
				dfs[si].append(pm[0])
				
			#for si in [0,1]:
				#scr=[compute_score((d, np.mean(pshift)), powerspecs[it], options, signs[si]) for d in defrg]
				#scores[si].append((np.mean(scr)-np.min(scr))*100)
				#dfs[si].append(defrg[np.argmin(scr)])
				
			print("ID {}, angle {:.1f}, defocus {:.1f} vs {:.1f}, score {:.3f} vs {:.3f}".format(it, tltparams[it,3], dfs[0][-1], dfs[1][-1], scores[0][-1], scores[1][-1]))
		
		print("Average score: Current hand - {:.3f}, flipped hand - {:.3f}".format(np.mean(scores[0]), np.mean(scores[1])))
		print("Defocus std: Current hand - {:.3f}, flipped hand - {:.3f}".format(np.std(dfs[0]), np.std(dfs[1])))
		scr=np.array(scores)
		scr=np.mean(scr[0]>scr[1])
		print("Current hand is better than the flipped hand in {:.1f}% tilt images".format(scr*100))
		
		if scr>.5:
			print("The handedness (--tltax={:.1f}) seems to be correct. Rerun CTF estimation without the checkhand option to finish the process.".format(-rot))
		else:
			print("The handedness seems to be flipped. Consider rerun the tomogram reconstruction with --tltax={:.1f} then rerun the CTF estimation.".format(-((180+rot)%360)))
		      
		if not options.nolog: E2end(logid)
		return
		
			
		
	
	nref=options.nref
	if nref>0:
		if options.verbose>0: print("Using first {} images near center tilt to estimate defocus range...".format(nref))
	else:
		if options.verbose>0: print("Estimating defocus...")

	ctfparams=np.zeros((len(tltparams), 2))
	for it in tltsrt:
		allrd, pzus=powerspecs[it]
		
		allscr=[]
		for ic, ctf in enumerate(allctf):
			scr=calc_all_scr(allrd, ctf, zlist[ic], box, exclude)
			idxsft=np.round(-np.array(pzus)/defstep).astype(int)
			stilt=np.zeros(len(defrg))+np.inf
			for i,df in enumerate(defrg):
				idx=idxsft+i
				outb=(idx<0)+(idx>=len(defrg))
				idx[outb]=0
		
				s=scr[idx, np.arange(scr.shape[1])].copy()
				if len(s)==0:
					stilt[i]=1
					continue
				sinf=np.isinf(s)
				s[outb]=np.max(s)
				sval=s[np.isinf(s)==0]

				if len(sval)==0:
					stilt[i]=1
				else:
					stilt[i]=np.sum(sval)/len(s)
			
			allscr.append(stilt)
			
		allscr=np.array(allscr)
		#print(allscr)
		if np.min(allscr)>0.8:
			if len(exclude)>0:
				d=np.mean(defrg[dfsel])
			else:
				d=np.mean(defrg)
			pm=(d, np.mean(pshift))
			
		else:
			amp, df= np.array(np.where(allscr==np.min(allscr))).T[0]
			
			pm=(defrg[df], pshift[amp])
		if options.refine:
			res=minimize(compute_score, pm, (powerspecs[it], options),method='Nelder-Mead',options={'xtol': 1e-3, 'disp': False, "maxiter":10})
			pm=res.x
		
		if options.verbose>0: print("ID {}, angle {:.1f}, defocus {:.3f}, phase shift {:.1f}, score {:.1f}".format(it, tltparams[it,3], pm[0], pm[1], np.min(allscr)*1000))
		dfs.append(pm[0])
		ctfparams[it]=[pm[1], pm[0]]

		#### get enough references
		if len(dfs)==nref: 
			if options.verbose>0: print("In the first {} image, defocus mean {:.2f}, std {:.2f}".format(nref, np.mean(dfs), np.std(dfs)))
			#if (np.mean(dfs)+np.std(dfs)*3>defrg[int(len(defrg)*.9)]):
			if (np.std(dfs)>1):
				if options.verbose>0: print ("Warning: variance of defocus estimation is larger than normal. Something may be wrong...")
				#return
			
			searchrg=min(max(np.std(dfs)*3, 3), 1.0)
			dfsel=abs(defrg-np.mean(dfs))<searchrg
			
			
			if options.verbose>0: print("We will search defocus range {:.2f} to {:.2f} for the rest images.".format(np.min(defrg[dfsel]), np.max(defrg[dfsel])))
			
			exclude=np.where(dfsel==False)[0]
			
			
	js=js_open_dict(info_name(tfile))
	js["defocus"]=ctfparams[:,1].tolist()
	js["phase"]=ctfparams[:,0].tolist()
	js["cs"]=float(options.cs)
	js["voltage"]=float(options.voltage)
	js=None
	
	print("{} done. Average defocus {:.2f}".format(tfile, np.mean(ctfparams[:,1])))
	
	if not options.nolog: E2end(logid)
	
def run(cmd):
	#print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()

