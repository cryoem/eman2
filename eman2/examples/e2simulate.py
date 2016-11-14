#!/usr/bin/env python

# Author: Michael Bell 09/2015

from EMAN2 import *
import numpy as np
import os
import multiprocessing

def main():
	usage="Script to simulate a micrograph using a PDB structure. Example usage: e2simem.py groel.mrc"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--seed", type=int, help="Random seed to use for noise",default=None)
	parser.add_argument("--noiseamp", type=float, help="Pink noise amplitude",default=0.8)
	parser.add_argument("--noiseampwhite", type=float, help="White noise amplitude",default=0.35)
	parser.add_argument("--voltage", type=float, help="Voltage to simulate",default=300)
	parser.add_argument("--apix", type=float, help="Angstroms per pixel sampling",default=1.0)
	parser.add_argument("--cs", type=float, help="Spherical aberration to simulate",default=2.7)
	parser.add_argument("--bfactor", type=float, help="Bfactor to simulate",default=200)
	parser.add_argument("--defocus", type=float, help="Amount of defocus to simulate",default=1.4)
	parser.add_argument("--ampcont", type=float, help="Amount of amplitude contrast to simulate",default=0.1)
	parser.add_argument("--nptcls","-n",type=int, help="Number of particles to include in simulated micrograph",default=80)
	parser.add_argument("--resolution", type=int, help="Resolution to use during pdb to mrc conversion. Defaults to 2.5*apix.",default=None)
	parser.add_argument("--sym", type=str, help="Symmetry of specimen. Default is c1.",default="c1")
	parser.add_argument("--grid",action="store_true",default=False,help="Place particles in a grid")
	parser.add_argument("--scale",default=1.0,type=float,help="Factor by which volume will be scaled prior to projection")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)

	if len(args) < 1:
		print("Need to specify a PDB structure or a MRC/HDF map.")
		exit(1)

	sym=options.sym

	fname = args[0]
	base = args[0].split('.')[0]
	ext = args[0].split('.')[1]

	if "pdb" in ext:
		if options.verbose: print("Converting the input PDB file to MRC format")
		fname = base + '.mrc'
		if options.resolution: res = options.resolution
		else: res = options.apix * 2.5
		os.system('e2pdb2mrc.py {} {} --apix {} --res {}'.format(args[0],fname,options.apix,res))
		struct=EMData(fname)
		sname=fname
	else:
		struct=EMData(args[0])
		sname=args[0]

	xb=struct.get_xsize()
	yb=struct.get_ysize()
	zb=struct.get_zsize()

	if xb != zb or yb != zb or xb != yb:
		if options.verbose: print("Making structure's rectangular box larger and cubic to ensure quality projections")
		ns=max(xb,yb,zb)
		n=int(ns/2)
		clip='{}_{}.hdf'.format(base,ns)
		os.system('e2proc3d.py {} {} --clip {},{},{},{},{},{}'.format(fname,clip,ns,ns,ns,n,n,n))
		struct=EMData(clip)
		sname=clip

	struct.write_image(sname)

	xb=struct.get_xsize()
	yb=struct.get_ysize()

	xt=4096
	yt=4096

	if xt <= xb or yt <= yb:
		print("Need to use a larget micrograph x and/or y dimension.")
		exit(1)

	n=options.nptcls
	if n <= 0:
		print("This program is designed to simulate one or more particles")
		exit(1)

	edg = 1.5
	if options.grid:
		coords = []
		for i in np.linspace(xb,xt-edg*xb,(xt-edg*xb)/xb):
			for j in np.linspace(yb,yt-edg*yb,(yt-edg*yb)/yb):
				coords.append([i,j])
				if i != j:
					coords.append([j,i])
	else:
		xs=np.round(np.random.uniform(xb,xt-edg*xb,size=n),0)
		ys=np.round(np.random.uniform(yb,yt-edg*yb,size=n),0)
		coords=np.vstack([xs,ys]).transpose()

	nprjs = len(coords)
	prjs="{}_prjs.hdf".format(base)
	if not os.path.isfile(prjs):
		if options.verbose: print("Randomly projecting map {} times".format(nprjs))
		cores=multiprocessing.cpu_count()
		os.system('e2project3d.py {} --outfile={} --orientgen=rand:n={}:phitoo=true --sym={} --projector=standard --parallel=thread:{}'.format(sname,prjs,nprjs,sym,cores))
	else:
		print("Using existing file of projections. If you wish to generate new (or more) projections, simply move or delete {}.".format(prjs))

	if options.verbose: print("Inserting projections into micrograph")

	t = Transform({"type":"eman","scale":options.scale})
	mg = EMData(xt,yt)
	for i,c in enumerate(coords):
		prj = EMData(xt,yt)
		prj.transform(t)
		prj.insert_clip(EMData(prjs,i),(int(c[0]),int(c[1])))
		mg.add(prj)

	if options.grid: outfile = 'sim_{}_{}_grid.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}.hdf'.format(base,n)

	mg.write_image(outfile)

	if options.verbose: print("Adding simulated CTF to micrograph")

	ctf = {}
	if options.apix: ctf['apix']=options.apix
	if options.voltage: ctf['voltage']=options.voltage
	if options.defocus: ctf['defocus']=options.defocus
	if options.ampcont: ctf['ampcont']=options.ampcont
	if options.bfactor: ctf['bfactor']=options.bfactor-(4*res**2)
	else: ctf['bfactor']=200-(4*res**2)
	if options.cs: ctf['cs']=options.cs

	mgctf = mg.process('math.simulatectf',ctf)

	if options.grid: outfile = 'sim_{}_{}_grid_CTF.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}_CTF.hdf'.format(base,n)

	mgctf.write_image(outfile)

	if options.verbose: print("Adding noise to micrograph")

	ctf['noiseamp']=options.noiseamp
	ctf['noiseampwhite']=options.noiseampwhite

	mg.process_inplace('math.simulatectf',ctf)

	if options.grid: outfile = 'sim_{}_{}_grid_CTF_noise.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}_CTF_noise.hdf'.format(base,n)

	mg.write_image(outfile)

	E2end(logid)

if __name__ == "__main__":
	main()

