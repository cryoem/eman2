#!/usr/bin/env python
# Muyuan Chen 2017-03
from EMAN2 import *
import numpy as np
import threading
import Queue

def make3d(ii, options, ptcls):
	
	boxsize=ptcls[0]["nx"]
	pad=good_size(boxsize*3/2)
	recon=Reconstructors.get("fourier", {"sym":"c1","size":[pad,pad,pad]})
	recon.setup()
	#print "{} started. {} particles".format(ii, len(ptcls))
	kk=0
	for p in ptcls:
		p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
		p3=recon.preprocess_slice(p2,p["xform.projection"])
		recon.insert_slice(p3,p["xform.projection"],1.0)
		kk+=1
	
	threed=recon.finish(True)
	if options.clip>0:
		boxsize=options.clip
		threed.clip_inplace(Region((pad-boxsize)/2, (pad-boxsize)/2, (pad-boxsize)/2, boxsize, boxsize,boxsize))
	#threed.write_image(options.ptclout, mid)
	threed["apix_x"]=ptcls[0]["apix_x"]
	threed["apix_y"]=ptcls[0]["apix_x"]
	threed["apix_z"]=ptcls[0]["apix_x"]
	threed["source_2dptcl"]=options.ptclin
	threed["ptcl_repr"]=kk
	threed.process_inplace("normalize.edgemean")
	try:
		threed["box3d"]=ptcls[0]["box3d"]
	except:
		pass
	
	if options.mask:
		msk,opt=parsemodopt(options.mask)
		threed.process_inplace(msk, opt)
	#print "{} finished. {} particles".format(ii, kk)
	options.queue.put([ii,threed])
	
def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptclin", type=str,help="2d particle input", default=None)
	parser.add_argument("--ptclout", type=str,help="3d particle output", default=None)
	parser.add_argument("--clip", type=int,help="final output size", default=-1)
	parser.add_argument("--mask", type=str,help="mask on the final output", default=None)
	parser.add_argument("--threads", type=int,help="Number of threads to use", default=5)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	print "Reading data..."
	num=EMUtil.get_image_count(options.ptclin)
	#num=1000
	#imgs=EMData.read_images(options.ptclin)
	imgs=[EMData(options.ptclin, i) for i in range(num)]
	nmod=imgs[-1]["model_id"]+1
	print "{:d} 2D particles, {:d} 3D particles output.".format(num, nmod)
	
	if options.ptclout==None:
		bname=os.path.basename(options.ptclin)
		
		options.ptclout="particles3d/"+bname[:-4]+"_make3d.hdf"
	
	try: os.remove(options.ptclout)
	except: pass

	jobs=[[] for i in range(nmod)]
	for i in range(num):
		mi=imgs[i]["model_id"]
		jobs[mi].append(imgs[i])
	#print jobs[0]
	print "Start working on {} threads...".format(options.threads)
	
	jsd=Queue.Queue(0)
	options.queue=jsd

	thrds=[threading.Thread(target=make3d,args=(i, options, j)) for i,j in enumerate(jobs)]

	# here we run the threads and save the results, no actual alignment done here
	print len(thrds)," threads"
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads+1 ) : time.sleep(.1)
			print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	
		while not jsd.empty():
			ii, threed=jsd.get()
			threed.write_image(options.ptclout,ii)
			#print ii,options.ptclout

	
	print "Output written to {}".format(options.ptclout)
	
	#pl=pool.Pool(options.threads)
	#ret=pl.map_async(make3d, jobs, chunksize=1)
	#pl.close()
	#pl.join()
	#maps=ret.get()
	
	#for m in maps:
		#m.write_image(options.ptclout, -1)
	
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	