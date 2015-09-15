#!/usr/bin/env python

# Author: Michael Bell 09/2015

from EMAN2 import *
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
	usage=""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--seed", type=int, help="Random seed to use for noise",default=None)
	parser.add_argument("--frames", type=int, help="Number of frames to simulate",default=20)
	parser.add_argument("--frames1", type=int, help="Number of frames to simulate charging",default=6)
	parser.add_argument("--shift1", type=float, help="Maximal desired initial shift",default=10.0)
	parser.add_argument("--shift2", type=float, help="Maximal desired secondary shift",default=10.0)
	parser.add_argument("--noiseseed", type=int, help="Random seed to use for noise",default=2015)
	parser.add_argument("--noise", type=float, help="Amount of noise to add",default=1.0)
	parser.add_argument("--showfig",action="store_true",default=False,help="Show translations quiver plot")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	if options.frames1 > options.frames:
		print("need a smaller number of 'charging' frames.")
		exit(1)
	if options.seed: numpy.random.seed(options.seed)
	
	f1 = options.frames1
	a = np.random.random(f1)
	a /= a.sum()
	a *= options.shift1
	
	f2 = options.frames-options.frames1
	b = np.random.random(f2)
	b /= b.sum()
	b *= options.shift2
	
	mags=np.concatenate([a,b])
	np.insert(mags,0,0)
	
	x1=make_rand_unit_vector()
	x2=make_rand_unit_vector()
	shifts = [[0,0]]
	xold = 0
	yold = 0
	with open(args[0][:-4]+"_trans.txt","w") as ts:
		for i in range(1,options.frames):
			if i < options.frames1:
				s = np.random.normal(0,0.4,2)
				x = xold + mags[i] * x1[0] + s[0]
				y = yold + mags[i] * x1[1] + s[1]
			elif i == options.frames1:
				s = np.random.normal(0,0.35,2)
				x = xold + mags[i] * x2[0] + s[0]
				y = yold + mags[i] * x2[1] + s[1]
			else:
				s = np.random.normal(0,0.3,2)
				x = xold + mags[i] * x2[0] + s[0]
				y = yold + mags[i] * x2[1] + s[1]
			xold = x
			yold = y
			shifts.append([x,y])
			ts.write("{}\t{}\n".format(shifts[i][0],shifts[i][1]))
	plot_translations(options,shifts,fname=args[0][:-4]+'_trans.png')
	for i,s in enumerate(shifts):
		img = EMData(args[0])
		img.process_inplace('math.addnoise',{'noise':options.noise,'seed':options.noiseseed})
		img.translate(s[0],s[1],0)
		img.write_image(args[0][:-4]+'_sim.hdf',i)
	E2end(logid)

def make_rand_unit_vector(dims=2):
	vec = [np.random.normal(0, 1) for i in range(dims)]
	mag = sum(x**2 for x in vec) ** .5
	return [x/mag for x in vec]

def plot_translations(options,shifts,fname=None):
	plt.figure()
	ax=plt.gca()
	shifts=np.asarray(shifts[1:])
	ax.plot(shifts[:,0],shifts[:,1],c='k',alpha=0.3,zorder=0)
	ax.scatter(shifts[:,0],shifts[:,1],c=range(len(shifts)),s=500,zorder=1)
	for i,s in enumerate(shifts):
		ax.annotate(str(i+1),xy=(s[0],s[1]),xytext=(-5,-5),textcoords='offset points',color='k')
	plt.savefig(fname)
	if options.showfig: plt.show()

if __name__ == '__main__':
	main()
	
