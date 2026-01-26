#!/usr/bin/env python
# Author: Lan Dang, 03/17/2022 (dlan@bcm.edu)
from math import *
from EMAN2 import *
import queue
import time
import numpy as np


def main():
	usage="""a program to run template matching on a tomogram.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--tomos",type=str,help="tomogram(s) for template matching",default="")
	parser.add_argument("--template",type=str,help="template",default="")
	parser.add_argument("--nbin",type=int, help="nbin to downsample tomogram for coarse search. Default=1", default=1)
	parser.add_argument("--nthreads",type=int, help="number of threads. Default=4", default=4)

	parser.add_argument("--dt",type=int, help="angular step for coarse search", default=45)
	parser.add_argument("--local_norm",action="store_true", help="Coarse search with local normalization", default=False)
	parser.add_argument("--ccc_fac",type=int, help="threshold to filter ccc: mean+ccc_fac*std. Default = 3", default=3)
	parser.add_argument("--npks",type=int, help="number of peaks for local search. Default = 0", default=0)
	parser.add_argument("--insert_mask",type=str,help="mask to insert back for annotation ",default="")
	parser.add_argument("--insert_density",type=str,help="density to reconstitute into volume with weights",default="")
	parser.add_argument("--insert_thresh",type=float,help="threshold to binary insert mask before annotation ",default=0.5)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)


	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)

	if len(options.tomos)==0 or len(options.template)==0:
		print("ERROR: Must specify tomogram and template for template matching")
		sys.exit(1)

	template=EMData(options.template)
	for tomo in options.tomos.split(","):
		target=EMData(tomo)
		base_out="{}_{}".format(base_name(tomo),base_name(options.template))
		print(base_out)
		calc_ccm(target,template,dt=options.dt,shrink=options.nbin,nthreads=options.nthreads,base_out=base_out,local_norm=options.local_norm)
		if options.npks > 0:
			do_local_search_test(target,template,options.npks,options.nbin,options.nthreads,options.ccc_fac,base_out)
		if len(options.insert_mask) > 0:
			try:
				insert_im = EMData(options.insert_mask)
			except:
				print("ERROR: invalid mask file")
				sys.exit(1)
			#d = json.load(open(f"{base_out}_local_pks.json","r"))
			js=js_open_dict(f"{base_out}_local_pks.json")
			pts = js["pts"]
			xfs = js["xfs"]
			# pts = d["pts"]
			# xfs = []
			# for xfm in d["xfs"]:
			# 	xf = Transform()
			# 	xf.set_matrix(xfm)
			# 	xfs.append(xf)
			vol_out = EMData(tomo,0,True)
			outfile = f"{base_out}_tmplt_seg.hdf"
			vol_out.write_image(outfile)

			insert_vol(vol_out,insert_im,pts,xfs,outfile)
			#js.close()
		if len(options.insert_density) > 0:
			try:
				insert_im = EMData(options.insert_density)
			except:
				print("ERROR: invalid mask file")
				sys.exit(1)
			#d = json.load(open(f"{base_out}_local_pks.json","r"))
			js=js_open_dict(f"{base_out}_local_pks.json")
			pts = js["pts"]
			xfs = js["xfs"]
			cfs = js["cfs"]
			# pts = d["pts"]
			# xfs = []
			# for xfm in d["xfs"]:
			# 	xf = Transform()
			# 	xf.set_matrix(xfm)
			# 	xfs.append(xf)
			vol_out = EMData(tomo,0,True)
			outfile = f"{base_out}_tmplt_dens.hdf"
			vol_out.write_image(outfile)
			insert_vol(vol_out,insert_im,pts,xfs,outfile,cfs)


def compute_local(jsd,targetf,template,targetsqf,templatemask,phi,alt,n):
	print("test compute local")
	nx,ny,nz=targetf["nx"]-2,targetf["ny"],targetf["nz"]
	clp=template["nx"]

	# rotate template, pad, fft
	trot=template.process("xform",{"transform":Transform({"type":"eman","phi":phi,"alt":alt})}).get_clip(Region((clp-nx)//2,(clp-ny)//2,(clp-nz)//2,nx,ny,nz))
	trot.process_inplace("xform.phaseorigin.tocorner")
	trotf=trot.do_fft()

	# actual CCF
	ccf=targetf.calc_ccf(trotf)

	# rotate, pad, fft for mask
	mrot=templatemask.process("xform",{"transform":Transform({"type":"eman","phi":phi,"alt":alt})}).get_clip(Region((clp-nx)//2,(clp-ny)//2,(clp-nz)//2,nx,ny,nz))
	mrot.process_inplace("xform.phaseorigin.tocorner")
	mrotf=mrot.do_fft()

	# CCF of mask with squared volume
	ccfn=targetsqf.calc_ccf(mrotf)
	ccfn.process_inplace("math.sqrt")
	ccf.div(ccfn)
	jsd.put((ccf,phi,alt,n))

def compute(jsd,targetf,template,phi,alt,n):
	print("test compute")
	nx,ny,nz=targetf["nx"]-2,targetf["ny"],targetf["nz"]
	clp=template["nx"]
	trot=template.process("xform",{"transform":Transform({"type":"eman","phi":phi,"alt":alt})}).get_clip(Region((clp-nx)//2,(clp-ny)//2,(clp-nz)//2,nx,ny,nz))
	trot.process_inplace("xform.phaseorigin.tocorner")
	trotf=trot.do_fft()
	ccf=targetf.calc_ccf(trotf)
	jsd.put((ccf,phi,alt,n))

def calc_ccm(target,template,dt,shrink,nthreads,base_out,local_norm=False):
	print("Running test calc_ccm")

	#minimum rotation : number of degree to rotate for 1 pixel shift at nyquist - 4.4 in this case of microtubule with boxsz=25
	if shrink > 1:
		target.process_inplace("math.fft.resample",{'n':shrink})

	target.mult(-1.0)
	nx,ny,nz=target["nx"],target["ny"],target["nz"]
	targetf=target.do_fft()
	targetsq=target.process("math.squared")
	targetsqf=targetsq.do_fft()

	templatesca=template.process("math.fft.resample",{"n":target["apix_x"]/template["apix_x"]})
	nxt1=templatesca["nx"]
	templatemask=templatesca.process("mask.auto3d",{"nmaxseed":8,"nshells":2,"radius":nxt1//10,"return_mask":True,"sigma":1.5})
	owner=EMData(target["nx"],target["ny"],target["nz"])
	avg=Averagers.get("minmax",{"max":True,"owner":owner})
	orts=[]

	jsd=queue.Queue(0)
	thrds=[]
	i=0

	for alt in range(0,90,dt):
		for phi in range(0,360,dt):
			if local_norm:
				thrds.append((jsd,targetf,templatesca,targetsqf,templatemask,phi,alt,i))
			else:
				thrds.append((jsd,targetf,templatesca,phi,alt,i))
			i+=1

	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
		if thrtolaunch<len(thrds):
			while (threading.active_count()>=nthreads) : time.sleep(0.1)
			print("\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)), end=' ',flush=True)
#			thrds[thrtolaunch]=threading.Thread(target=compute,args=thrds[thrtolaunch])		# replace args
			if local_norm:
				thrds[thrtolaunch]=threading.Thread(target=compute_local,args=thrds[thrtolaunch])		# replace args
			else:
				thrds[thrtolaunch]=threading.Thread(target=compute,args=thrds[thrtolaunch])		# replace args
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(0.1)

		# return is [N,dict] a dict of image# keyed processed images
		while not jsd.empty():
			ccf,phi,alt,i=jsd.get()
			ccf["ortid"]=len(orts)
			orts.append((phi,alt))
			ccf.process_inplace("normalize")
			avg.add_image(ccf)
			thrds[i].join()
			print(f"\n{phi},{alt} done {thrtolaunch} {threading.active_count()}/{nthreads}")

	peaks=avg.finish()


	peaks.write_image(f"{base_out}_ccc.hdf:-1",0)
	owner.write_image(f"{base_out}_owner_id.hdf:-1",0)
	#owner.write_image(f"{base}_own.hdf:-1",0)
	#target.write_image(f"{base}_targ.hdf:12")
	#templatesca.write_image(f"{base}_tmpl.hdf:12",0)
	out=open(f"{base_out}_orts.txt","w")
	for i,phialt in enumerate(orts): out.write(f"{i}\t{phialt[0]}\t{phialt[1]}\n")

def find_peaks_amp_test(npks,nbin,ccc_fac=3,ccc_file="",own_file="",orts_file=""):
	n_peaks= npks
	ccc_hdr=EMData(ccc_file, 0,True)
	cx,cy,cz = ccc_hdr.get_sizes()

	#rx = template.get_sizes()[0]
	#bsz = rx
	#ccc = EMData(ccc_file).process("mask.zeroedge3d",{"x0":bsz//2,"x1":bsz//2,"y0":bsz//2,"y1":bsz//2,"z0":bsz//2,"z1":bsz//2})
	ccc = EMData(ccc_file)
	print("Shape",ccc.get_sizes())

	#ccc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	#print("Set border of ccc matrix to zero to avoid peaks too close to the border")
	mean=float(ccc.get_attr("mean"))
	sig=float(ccc.get_attr("sigma"))
	#ccc.process_inplace("threshold.belowtozero",{"minval":mean+sig})
	ccc.process_inplace("mask.onlypeaks",{"npeaks":2}) #before threshold
	ccc.process_inplace("threshold.belowtozero",{"minval":mean+sig})


	peaks = ccc.calc_highest_locations(mean+ccc_fac*sig)
	peaks.sort(key=lambda peak: -peak.value)
	pts = [(nbin*peak.x,nbin*peak.y,nbin*peak.z) for peak in peaks]
	orts_f=open(orts_file,"r")
	oris = []
	for ort in orts_f.readlines():
		angles = ort.split("\t")
		xf = Transform({'az':0,'alt':int(angles[2]),'phi':int(angles[1]),'tx':0.00,'ty':0.00,'tz':0.00,'mirror':0,'scale':1.0000,'type':'eman'})
		oris.append(xf)

	initxf = EMData(own_file)
	xfs = [oris[int(initxf[(peak.x,peak.y,peak.z)])] for peak in peaks]
	if len(pts) > n_peaks:
		print("Find total {} peaks. Highest score {} peaks at {}".format(len(pts),n_peaks,pts[:n_peaks]))
		return pts[:n_peaks],xfs[:n_peaks]
	else:
		print("Find {} peaks at {}".format(len(pts),pts))
		return pts,xfs

def do_search_test(tomo,tempt,jsd, px_l):
	for pt_xf in px_l:
		p = pt_xf[0]
		initxf = pt_xf[1]
		x,y,z= (int(i) for i in p)
		bs = tempt["nx"]
		subvol = tomo.get_clip(Region(x-bs//2,y-bs//2,z-bs//2,bs,bs,bs))
		tempt_a = tempt.align("rotate_translate_3d_tree",subvol,{"initxform":[initxf],"sym":"c1"},'ccc.tomo')
		#tempt.process_inplace("xform.phaseorigin.tocorner")
		subvolf =subvol.do_fft()
		tempt_af = tempt_a.do_fft()
		cf = subvolf.calc_ccf(tempt_af)
		#print(tempt_a.numpy().shape,cf.numpy().shape)
		cf_score = np.max(abs(cf.numpy()))
		jsd.put((p,tempt_a["xform.align3d"],cf_score))

def do_local_search_test(tomo,tempt,n_peaks,nbin,nthrds,ccc_fac=3,base_out=""):
	ccc_file = f"{base_out}_ccc.hdf"
	own_file =f"{base_out}_owner_id.hdf"
	orts_file = f"{base_out}_orts.txt"
	pts,xfs = find_peaks_amp_test(npks=n_peaks*3,nbin=nbin,ccc_file=ccc_file,own_file=own_file,orts_file=orts_file)
	if len(pts) == 0:
		print("Find 0 peaks at the set ccc threshold. Consider a lower ccc_fac. Return")
		return
	mbin=tomo["apix_x"]/tempt["apix_x"]
	print("Will shrink reference by {:.1f} to match apix of tomogram for local search".format(mbin))
	tempt.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.4}) #do it first then math.fft.resample #get rid of lowres peak of ice gradient
	tempt.process_inplace("math.fft.resample",{'n':mbin})
	tempt.process_inplace('normalize.edgemean')
	tempt.mult(-1)
	good_szs=[good_size(tempt["nx"]),good_size(tempt["ny"]),good_size(tempt["nz"])]
	tempt.clip_inplace(Region(0,0,0,good_szs[0],good_szs[1],good_szs[2]))
	tomo.process_inplace('normalize.edgemean')

	pts_xfs_l = list(zip(pts,xfs))
	start = time.time()
	tasks=[pts_xfs_l[i::nthrds] for i in range(nthrds)]

	#print(tasks)
	jsd=queue.Queue(0)
	thrds=[threading.Thread(target=do_search_test,args=(tomo, tempt, jsd, px_l)) for px_l in tasks]

	thrtolaunch=0
	tsleep=threading.active_count()

	ndone=0
	for t in thrds:
		t.start()
	p_xf_l = []
	print("Start local searching at {} location".format(len(pts)))
	while threading.active_count()>tsleep or not jsd.empty():
		while not jsd.empty():
			p,xf,cf=jsd.get()
			p_xf_l.append((p,xf,cf))
			ndone+=1
			sys.stdout.write("\r{}/{} finished.".format(ndone, len(pts)))
			sys.stdout.flush()
		time.sleep(.5)
	print("")
	print("Finish in {}s".format(time.time()-start))

	p_xf_l.sort(key = lambda x: -x[2])
	# try:
	# 	pt_cf_l = []
	# 	for i in range(len(p_xf_l)):
	# 		pt = list(p_xf_l[i][0])
	# 		pt.append(p_xf_l[i][2])
	# 		pt_cf_l.append(pt)
	# 	#from_numpy(np.array(pt_cf_l)).write_image("tmplt_match_local_pks.hdf")
	# 	# from_numpy(np.array(pt_cf_l)).write_image(f"{base_out}_local_pks.hdf:-1",0)
	# 	# print("All peaks saved to"+base_out+"_local_pks.hdf")
	# except Exception as e:
	# 	print(e)
	# 	pass
	p_xf_l=p_xf_l[:n_peaks]
	pts = [p_xf[0] for p_xf in p_xf_l]
	#xfs = [p_xf[1].get_matrix() for p_xf in p_xf_l]
	xfs = [p_xf[1] for p_xf in p_xf_l]
	cfs = [p_xf[2] for p_xf in p_xf_l]
	if max(cfs) != 0:
		cfs = [cf/max(cfs)+0.01 for cf in cfs]

	#cfs = [p_xf[3].numpy().flatten() for p_xf in p_xf_l]
	d = {"pts" : pts,"xfs" : xfs, "cfs": cfs}
	js=js_open_dict(f"{base_out}_local_pks.json")
	js.update(d)
	js.close()
	#js.dump(d,open(f"{base_out}_local_pks.json","w"))


	#for i in range(len(pts)):
		#local_out.write(",".join(pts[i])+","+",".join(xfs[i].get_matrix()))
		#for i,phialt in enumerate(orts): local_out.write(f"{i}\t{phialt[0]}\t{phialt[1]}\n")
	#cfs_temp = [p_xf[2] for p_xf in p_xf_l]
	#cfs = [cf/max(cfs_temp) for cf in cfs_temp]
	#return pts,xfs


def insert_vol(target,insert_im,pts,xfs,outfile,cfs=[]):
	# outfile = f"{base_out}_tmplt_seg.hdf"
	#target.write_image(outfile)


	# if self.target.get_annotation() is not None:
	# 	anno_head = self.target.get_annotation().copy()
	# 	anno_head.write_image(outfile)
	# 	tomo_head = self.target.get_data().copy_head()
	tomo_apix = target["apix_x"]
	# else:
	# 	tomo_head = self.target.get_data().copy_head()
	# 	tomo_head.write_image(outfile)
	# 	tomo_apix = tomo_head["apix_x"]

	# if self.insert_vol_le.text() == 0:
	# 	if self.tmplt_le.text() > 0:
	# 		insert_im = EMData(self.tmplt_le.text())
	# 		print("Filter template to use as insert mask. ")
	# 	else:
	# 		print("Specify mask to put in the tomogram")
	# 	return
	# else:
	# 	try:
	# 		insert_im = EMData(self.insert_vol_le.text())
	# 	except:
	# 		print("Cannot read insert vol data in. Abort")


	mbin=tomo_apix/insert_im["apix_x"]
	insert_im.process_inplace("math.fft.resample",{'n':mbin})
	print("Will shrink insert vol by {:.1f} to match apix of tomogram".format(mbin))
	good_szs=[good_size(insert_im["nx"]),good_size(insert_im["ny"]),good_size(insert_im["nz"])]
	insert_im.clip_inplace(Region(0,0,0,good_szs[0],good_szs[1],good_szs[2]))
	print("Vol of insert im is {}".format(insert_im.get_sizes()))



	#annotate = self.target.get_annotation()
	#annotate = EMData("tmplt_match_temp_vol.hdf")

	for i,loc in enumerate(pts):
		if cfs[i] == 0:
			print("Cf score reach 0. Finish segmentation")
			break
		vol = insert_im.copy()
		#vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.033})
		vol.process_inplace("xform", {"transform":xfs[i]})
		vol.process_inplace("normalize.edgemean")
		# vol.process_inplace("threshold.binary",{"value":thres})

		#vol.mult(val)
		x,y,z = [int(i) for i in loc]
		# bs = vol["nx"]
		# reg = Region(x-bs//2,y-bs//2,z-bs//2,bs,bs,bs)
		reg = Region(x-good_szs[0]//2,y-good_szs[1]//2,z-good_szs[2]//2,good_szs[0],good_szs[1],good_szs[2])
		#reg = Region(x,y,z,good_szs[0],good_szs[1],good_szs[2])

		try:
			#current_vol = self.target.get_annotation().get_clip(reg)
			#current_vol = annotate.get_clip(reg)
			current_vol = EMData(outfile, 0, False, reg)
			#current_vol=annotate.get_rotated_clip(Transform(),reg)
			if len(cfs) > 0:
				insert_vol=current_vol + cfs[i]*vol
			else:
				insert_vol=current_vol.process("math.max",{"with":vol})
		except:
			insert_vol = vol
		try:
			#annotate.set_rotated_clip(Transform(),insert_vol)
			insert_vol.write_image(outfile, 0, IMAGE_HDF, False, reg)
			sys.stdout.write("\r{}/{} finished.".format(i+1, len(pts)))
			sys.stdout.flush()
			#print("Insert vol at location",x,y,z,"")
		except Exception as e:
			print(e)
			print("Invalid region at location",x,y,z,".Skip.")
			continue
	#annotate = EMData(outfile)
	# self.target.img_view.set_data(self.target.data,EMData(outfile))
	# self.target.img_view.updateGL()
	print("Finishing segmentation")

if __name__ == '__main__':
	main()
