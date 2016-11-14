#!/usr/bin/env python
# Muyuan Chen 2016-09
from EMAN2 import *
import numpy as np
from e2ddd_movie import qsum 

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_pos_argument(name="movies",help="List the movies to align.", default="")
	#parser.add_argument("--ptcls", type=str,help="", default=None)
	parser.add_argument("--frameid", type=str,help="", default=None)
	parser.add_argument("--newsuffix", type=str,help="", default="")
	parser.add_argument("--invert", type=int, default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#ptcls=options.ptcls.split(',')
	ptcls=args
	frameid=[int(i) for i in options.frameid.split(',')]
	
	if options.newsuffix=="":
		print "overwriting the particles..."
		
	for pp in ptcls:
		boxes=get_particles(pp,frameid, options.invert)
		print pp,len(boxes)
		for i,b in enumerate(boxes):
			b.write_image(pp[:-4]+options.newsuffix+".hdf",i)
		
	E2end(logid)

def get_particles(pp,frameid, doinvert=-1):
	#### locate corresponding movie and gain/dark refs

	db=js_open_dict(info_name(pp))
	try: 
		mvname=str(db["movie_name"])
		
	except:
		print "Cannot locate movie file for {}, continue".format(pp)
		db.close()
		return []
	
	gain=dark=None
	if db.has_key("gain_name"):
		gain=EMData(str(db["gain_name"]), db["gain_id"])
	if db.has_key("dark_name"):
		dark=EMData(str(db["dark_name"]), db["dark_id"])
	
	
	#print pp, mvname, gain, dark
	e=EMData(pp,0,True)
	boxsize=e["nx"]
	locs=db["movieali_trans"]
	
	#### load movie frames
	#print mvname
	nframe=EMUtil.get_image_count(mvname)
	
	if mvname[-4:].lower() in (".mrc") :
		hdr=EMData(mvname,0,True)			# read header
		nx,ny=hdr["nx"],hdr["ny"]
		nframe=hdf["nz"]
	
	outim=[0]*nframe	
	for ii in frameid:

		#if fsp[-4:].lower() in (".mrc","mrcs") :
		if mvname[-4:].lower() in (".mrc") :
			im=EMData(mvname,0,False,Region(0,0,ii,nx,ny,1))
		else: im=EMData(mvname,ii)

		if dark!=None : im.sub(dark)
		if gain!=None : im.mult(gain)
		im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
		outim[ii]=im

	#### translate frames
	#for ii in frameid:
		#outim[ii].translate(int(floor(locs[ii*2]+.5)),int(floor(locs[ii*2+1]+.5)),0)
	
	micrograph=qsum([outim[i] for i in frameid])
	boxsize2=boxsize/2

	boxes=db["boxes"]
	
	# remove any existing file
	try: os.unlink(ptcl)
	except: pass
	ptclout=[]
	for i,b in enumerate(boxes):
		boxim=micrograph.get_clip(Region(b[0]-boxsize2,b[1]-boxsize2,boxsize,boxsize))
		boxim["ptcl_source_coord"]=(b[0],b[1])
		boxim["ptcl_source_image"]=mvname
		
		if i==0 and doinvert==-1:
			try:
				e=EMData(pp,i)
				cc0= e.cmp("ccc",boxim)
				e.mult(-1)
				cc1= e.cmp("ccc",boxim)
				print cc0, cc1
				if cc1<cc0:
					doinvert=1
				else:
					doinvert=0
			except:
				pass
		
		if doinvert:
			boxim.mult(-1)
		boxim.process_inplace("normalize.edgemean")
		ptclout.append(boxim)
		
	
	db.close()
	return ptclout

	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	