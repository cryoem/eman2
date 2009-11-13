#!/usr/bin/env python

from EMAN2    import *
from sparx    import *
from random   import random, seed, randint
from sys      import argv, exit
from optparse import OptionParser

usage = "generate_projection.py --CTF --format=(bdb|hdf)"
parser = OptionParser(usage)
parser.add_option( "--CTF",    action="store_true", default=False, help="whether generate CTF images")
parser.add_option( "--format", type="string",       default="bdb", help="data format: hdf or bdb")
parser.add_option( "--apix",   type="float",        default="2.5", help="pixel size")

(options,args) = parser.parse_args( argv )

seed(14567)

delta = 40
angles = even_angles(delta,0.0,89.9,0.0,359.9,"P")
nangle = len(angles)
#exit()


modelvol = EMData()
modelvol.read_image("../model_structure.hdf")

nx = modelvol.get_xsize()

nvol = 10
volfts = [None]*nvol
for i in xrange(nvol):
	sigma = 1.5 + random() # 1.5-2.5
	addon = model_gauss(sigma, 64, 64, 64, sigma, sigma, 38, 38, 40 )
	scale = 2500 * (0.5+random())
	model = modelvol + scale*addon
	volfts[i],kb = prep_vol(modelvol + scale*addon)


if options.format=="bdb":
	stack_data = "bdb:data"
	delete_bdb(stack_data)
else:
	stack_data = "data.hdf"
Cs      = 2.0
pixel   = options.apix
voltage = 120.0
ampcont = 10.0
ibd     = 4096/2-64
iprj    = 0

width = 240
xstart = 8 + 64/2
ystart = 8 + 64/2
rowlen = 17




params = []
for idef in xrange(3,8):

	irow = 0
	icol = 0

	mic = model_blank(4096, 4096)
	defocus = idef*0.2
	if options.CTF :
		ctf = EMAN2Ctf()
		ctf.from_dict( {"defocus":defocus, "cs":Cs, "voltage":voltage, "apix":pixel, "ampcont":ampcont, "bfactor":0.0} )

	for i in xrange(nangle):
		for k in xrange(24):
			dphi = 8.0*(random()-0.5)
			dtht = 8.0*(random()-0.5)
			psi  = 360.0*random()

			phi = angles[i][0]+dphi
			tht = angles[i][1]+dtht

			s2x = 4.0*(random()-0.5)
			s2y = 4.0*(random()-0.5)


			params.append( [phi, tht, psi, s2x, s2y])

			ivol = iprj%nvol
			proj = prgs(volfts[ivol], kb, [phi, tht, psi, -s2x, -s2y])
		
			x = xstart + irow * width
			y = ystart + icol * width
	
			mic += pad(proj, 4096, 4096, 1, 0.0, x-2048, y-2048, 0)
	
			
			proj = proj + model_gauss_noise( 30.0, nx, nx )
			if options.CTF :
				proj = filt_ctf(proj, ctf)
				proj.set_attr_dict({"ctf":ctf, "ctf_applied":0})

			proj = proj + filt_gaussl(model_gauss_noise(17.5,nx,nx), 0.3)
			proj.set_attr( "active", 1 )
			# flags describing the status of the image (1 = true, 0 = false)
			set_params2D(proj, [0.0, 0.0, 0.0, 0, 1.0])
			set_params_proj(proj, [phi, tht, psi, s2x, s2y])

			proj.write_image(stack_data, iprj)
			
			icol += 1
			if icol == rowlen:
				icol = 0
				irow += 1

			iprj += 1


	mic += model_gauss_noise(30.0,4096,4096)
	if options.CTF :
		#apply CTF
		mic = filt_ctf(mic, ctf)

	mic += filt_gaussl(model_gauss_noise(17.5,4096,4096), 0.3)
	
	mic.write_image("mic%1d.hdf"%(idef-3),0)

drop_spider_doc("params.txt", params)
