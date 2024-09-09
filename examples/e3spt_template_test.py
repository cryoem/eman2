# this program is a quick test of using exhaustive search template matching to find objects in tomograms. It is just a prototype
# not ready for actual use.
# 9/3/24 Steve Ludtke

from math import *
from EMAN3 import *
import queue
import time
#from EMAN3tensor import *


nthreads=32
base="out83tn_sigma"

# standard CCF variant, this computes the CCF for a single orientation in a thread and returns the result on jsd
def compute(jsd,targetf,template,phi,alt,n):
	nx,ny,nz=targetf["nx"]-2,targetf["ny"],targetf["nz"]
	clp=template["nx"]
	trot=template.process("xform",{"transform":Transform({"type":"eman","phi":phi,"alt":alt})}).get_clip(Region((clp-nx)//2,(clp-ny)//2,(clp-nz)//2,nx,ny,nz))
	trot.process_inplace("xform.phaseorigin.tocorner")
	trotf=trot.do_fft()
	ccf=targetf.calc_ccf(trotf)
#	ccf=trotf.calc_ccf(targetf,fp_flag.CIRCULANT,True)
	jsd.put((ccf,phi,alt,n))

# standard CCF with local normalization
def compute_local(jsd,targetf,template,targetsqf,templatemask,phi,alt,n):
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

#	ccf=trotf.calc_ccf(targetf,fp_flag.CIRCULANT,True)
	jsd.put((ccf,phi,alt,n))

# FLCF variant
def compute_flcf(jsd,target,template,phi,alt,n):
	trot=template.process("xform",{"transform":Transform({"type":"eman","phi":phi,"alt":alt})})
#	trot.process_inplace("xform.phaseorigin.tocorner")
#	trotf=trot.do_fft()
	ccf=target.calc_flcf(trot)
#	ccf.process_inplace("xform.phaseorigin.tocenter")
#	ccf=trotf.calc_ccf(targetf,fp_flag.CIRCULANT,True)
	jsd.put((ccf,phi,alt,n))


def main():
	target=EMData("lamella_8_3__bin8.hdf",0)
	target.mult(-1.0)
	nx,ny,nz=target["nx"],target["ny"],target["nz"]
	targetf=target.do_fft()
	targetsq=target.process("math.squared")
	targetsqf=targetsq.do_fft()

#	template=EMData("long_tmplt.hdf",0)
	template=EMData("long_tmplt_fromtomo.hdf",0)
	templatesca=template.process("math.fft.resample",{"n":target["apix_x"]/template["apix_x"]})
	templatemask=templatesca.process("mask.auto3d",{"nmaxseed":5,"nshells":2,"radius":5,"return_mask":True,"sigma":1.5})
	nxt1=templatesca["nx"]
	templatemask=templatesca.process("mask.auto3d",{"nmaxseed":8,"nshells":2,"radius":nxt1//10,"return_mask":True,"sigma":1.5})

	owner=EMData(target["nx"],target["ny"],target["nz"])
	avg=Averagers.get("minmax",{"max":True,"owner":owner})
	orts=[]

	jsd=queue.Queue(0)
	# these start as arguments, but get replaced with actual threads
	thrds=[]
	i=0
	for alt in range(0,90,3):
		for phi in range(0,360,3):
#			thrds.append((jsd,targetf,templatesca,phi,alt,i))
			thrds.append((jsd,targetf,templatesca,targetsqf,templatemask,phi,alt,i))
			i+=1

	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
		if thrtolaunch<len(thrds):
			while (threading.active_count()>=nthreads) : time.sleep(0.1)
			print("\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)), end=' ',flush=True)
#			thrds[thrtolaunch]=threading.Thread(target=compute,args=thrds[thrtolaunch])		# replace args
			thrds[thrtolaunch]=threading.Thread(target=compute_local,args=thrds[thrtolaunch])		# replace args
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
	global base
	peaks.write_image(f"{base}.hdf:-1",0)
	owner.write_image(f"{base}_own.hdf:-1",0)
	target.write_image(f"{base}_targ.hdf:12")
	templatesca.write_image(f"{base}_tmpl.hdf:12",0)
	out=open(f"{base}_orts.txt","w")
	for i,phialt in enumerate(orts): out.write(f"{i}\t{phialt[0]}\t{phialt[1]}\n")


main()
