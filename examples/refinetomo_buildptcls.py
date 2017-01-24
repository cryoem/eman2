#!/usr/bin/env python
# Muyuan Chen 2016-09
from EMAN2 import *
import numpy as np
import scipy.ndimage.filters as filters

def main():
	
	usage=""" make 2d projections from sub-tomograms for refinetomo_easy.py. 
	[prog] --ptclin [3d subvolume stack input] --ptclout [2d projection stack output]
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptclin", type=str,help="3d volume of particles", default=None)
	parser.add_argument("--ptclout", type=str,help="output 2D projections of the input 3D particles", default=None)
	parser.add_argument("--tltang", type=float,help="tilt angle range", default=80.)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)



	ptclfile=options.ptclin
	num=EMUtil.get_image_count(ptclfile)

	pjfile=options.ptclout
	try: 
		os.remove(pjfile)
		print "Overwriting...."
	except: pass
	
	rotmat=[]
	deg=[]
	for tt in np.arange(-80,80,.2):
		deg.append(tt)
		t=tt/180.*np.pi
		rotmat.append([[cos(t), -sin(t)],[sin(t), cos(t)]])
	rotmat=np.array(rotmat)
	deg=np.array(deg)
	
	e=EMData(ptclfile,0)
	img=e.numpy().copy()
	imft=np.abs(get_fft(np.sum(img,axis=1)))
	rotpts=[]
	cnt=np.argmax(np.sum(imft,axis=1))
	sz=len(imft)
	for ti in range(len(deg)):
		pts=np.array([np.arange(sz)-cnt, np.zeros(sz)], dtype=int).T
		pts=np.dot(pts,rotmat[ti]).astype(int)
		pts+=cnt
		rotpts.append(pts)
	rotpts=np.array(rotpts)
	rotpts=np.clip(rotpts, 0, sz-1)
	for ptid in range(num):
		e=EMData(ptclfile,ptid)
		img=e.numpy().copy()
		imft=np.abs(get_fft(np.sum(img,axis=1)))
		tltsum= np.mean(imft[rotpts[:,:,1], rotpts[:,:,0]], axis=1)
		# tltsum[1:-1]=tltsum[
		thr= tltsum[len(tltsum)/2]/2
		pp=np.array(find_peaks(tltsum,3,thr))
		ph=tltsum[pp]
		pks=deg[pp]
	#	 plt.plot(deg,tltsum);
	#	 plt.plot(pks,ph,'.r');
		# plt.imshow(imft)
		
		for p in pks:
			tr=Transform({"type":'xyz', "ytilt":p})
			pj=e.project("standard", tr)
			pj["model_id"]=ptid
			pj.write_image(pjfile,-1)
		
		print "Processing particle # {}/{}. Number of projections: {}".format(ptid,num,len(pks))


	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
def get_fft(img):
	return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))

def get_img(fft):
	return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fft)).real)

def find_peaks(arr,wid=3,thr=0):
	arr=filters.gaussian_filter1d(arr,wid)
	pks=np.array(np.where(np.logical_and((arr[1:-1]-arr[:-2])>0,(arr[1:-1]-arr[2:])>0))[0])
	return pks[arr[pks]>thr]

if __name__ == '__main__':
	main()
	