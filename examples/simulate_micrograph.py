#!/usr/bin/env python

# Author: Michael Bell 09/2015

from EMAN2 import *
import numpy as np
import os

def main():
	usage="Script to simulate a micrograph using a PDB structure. Example usage: e2simem.py groel.mrc"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--seed", type=int, help="Random seed to use for noise",default=None)
	parser.add_argument("--noiseamp", type=float, help="Pink noise amplitude",default=None)
	parser.add_argument("--noiseampwhite", type=float, help="White noise amplitude",default=None)
	parser.add_argument("--voltage", type=float, help="Voltage to simulate",default=300)
	parser.add_argument("--apix", type=float, help="Angstroms per pixel sampling",default=1.0)
	parser.add_argument("--cs", type=float, help="Spherical aberration to simulate",default=2.7)
	parser.add_argument("--bfactor", type=float, help="Bfactor to simulate",default=None)
	parser.add_argument("--defocus", type=float, help="Amount of defocus to simulate",default=None)
	parser.add_argument("--ampcont", type=float, help="Amount of amplitude contrast to simulate",default=10.0)
	parser.add_argument("--mgx", type=int, help="Micrograph width in pixels",default=4096)
	parser.add_argument("--mgy", type=int, help="Micrograph height in pixels",default=4096)
	parser.add_argument("--nptcls", type=int, help="Number of particles to include in simulated micrograph",default=80)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--grid",action="store_true",default=False,help="Place particles in a grid")
	parser.add_argument("--noise", type=float, help="Level of noise to be added after CTF",default=75.0)
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	if len(args) < 1:
		print("Need to specify a PDB structure or a MRC/HDF map.")
		exit(1)
	
	base = args[0].split('.')[0]
	ext = args[0].split('.')[1]
	
	if "pdb" in ext:
		fname = base + '.mrc'
		print("Converting the input PDB file to MRC format")
		os.system('e2pdb2mrc.py {} {}'.format(args[0],fname))
		struct=EMData(fname)
	else:
		struct=EMData(args[0])
	
	xb=struct.get_xsize()
	yb=struct.get_ysize()
	zb=struct.get_zsize()
	
	if xb != zb or yb != zb or xb != yb:
		print("Making structure's rectangular box larger and cubic for high quality projections")
		ns=2*max(xb,yb,zb)
		n=int(ns/2)
		clip='{}_{}.hdf'.format(base,ns)
		os.system('e2proc3d.py {} {} --clip {},{},{},{},{},{}'.format(fname,clip,ns,ns,ns,n,n,n))
		struct=EMData(clip)
	
	xb=struct.get_xsize()
	yb=struct.get_ysize()
	xt=options.mgx
	yt=options.mgy
	n=options.nptcls
	
	if xt <= xb or yt <= yb:
		print("Need to use a larget micrograph x and/or y dimension.")
		exit(1)
	if n < 0:
		print("This program is designed to simulate one or more particles")
		exit(1)

	if options.verbose: print("Randomly projecting map throughout micrograph")

	if options.grid:
		coords = []
		for i in np.linspace(xb,xt-2*xb,(xt-2*xb)/xb):
			for j in np.linspace(yb,yt-2*yb,(yt-2*yb)/yb):
				coords.append([i,j])
				if i != j:
					coords.append([j,i])
	else:
		xs=np.round(np.random.uniform(xb,xt-2*xb,size=n),0)
		ys=np.round(np.random.uniform(yb,yt-2*yb,size=n),0)
		coords=np.vstack([xs,ys]).transpose()
	
	mg_final = EMData(xt,yt)
	params = 180*np.random.sample((len(coords),3))
	i=1
	tot=len(params)
	for p,c in zip(params,coords):
		if options.verbose: print "%i/%i\r"%(i,tot),
		t = Transform({"type":"eman","az":p[0],"alt":p[1],"phi":p[2]})
		prj = struct.project("standard",t)
		mg = EMData(xt,yt)
		mg.insert_clip(prj,(int(c[0]),int(c[1])))
		mg_final.add(mg)
		i+=1

	if options.grid: outfile = 'sim_grid_{}_{}.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}.hdf'.format(base,n)

	mg_final.write_image(outfile)
	
	if options.verbose: print("Adding simulated CTF to micrograph")
	
	ctf = {}
	if options.apix: ctf['apix']=options.apix
	if options.voltage: ctf['voltage']=options.voltage
	if options.defocus: ctf['apix']=options.defocus
	if options.ampcont: ctf['apix']=options.ampcont
	if options.bfactor: ctf['apix']=options.bfactor
	if options.cs: ctf['apix']=options.cs
	if options.noiseamp: ctf['apix']=options.noiseamp
	if options.noiseampwhite: ctf['apix']=options.noiseampwhite

	mg_final.process_inplace('math.simulatectf',ctf)

	if options.grid: outfile = 'sim_grid_{}_{}_CTF.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}_CTF.hdf'.format(base,n)

	mg_final.write_image(outfile)
	
	if options.verbose: print("Adding noise to simulated micrograph")
	
	mg_final.process_inplace('math.addnoise',{'noise':options.noise})

	if options.grid: outfile = 'sim_grid_{}_{}_CTF_noise.hdf'.format(base,n)
	else: outfile = 'sim_{}_{}_CTF_noise.hdf'.format(base,n)

	mg_final.write_image(outfile)
	
	E2end(logid)

if __name__ == "__main__":
	main()

