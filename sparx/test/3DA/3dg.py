#!/usr/bin/env python

from EMAN2  import *
from sparx  import *

from sys import exit

s2x = 10
s2y = -21
phi = 130
theta = 25
psi = 73
model = test_image_3d()
t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
t.set_pre_trans(Vec2f(s2x,s2y))
print  t.get_pre_trans()

#t.set_scale(2.0) works but probably won't be used
#t.set_mirror(1) works but probably won't be used
proj = model.project("standard",t)
proj_2 = model.project("gauss_fft",t)
projl = project(model, [phi, theta, psi, -s2x, -s2y],22.0)
projg = prg(model, [phi, theta, psi, -s2x, -s2y])
#prl = prl(model, [[phi, theta, psi, -s2x, -s2y]],22.0)
dropImage(proj,"p1.hdf")
dropImage(proj_2,"p2.hdf")
dropImage(projl,"p33.hdf")
dropImage(projg,"p4.hdf")
dropImage(proj-projg,"p5.hdf")
dropImage(proj_2-projg,"p6.hdf")
#dropImage(prl,"p5.hdf")
print  ccc(proj,proj_2),ccc(proj_2,projg),ccc(proj,projg),ccc(projl,projg)

# test of shifting in-plane

projt = prg(model, [phi, theta, psi, 0, 0])

#  This works and proves that the shift stored in the header (s2x, s2y) can be
#  applied to the projection to have it centered.
print ccc(projt, fshift(projg, s2x, s2y))

exit()
vol = EMData()
vol.read_image("../model001.tcp")
#vol.read_image("model.spi")
#info(e)
nx = vol.get_xsize()
delta_theta=7.0
'''
angles=even_angles(delta_theta,0.,90.,0,359.99,"S")
i1 = len(angles)
at = even_angles(delta_theta,0.,180.,0,359.99,"S")
i2 = len(at)

for i in xrange(i2-i1):
	angles.append(at[i+i1-1])
'''
angles=even_angles(delta_theta,0.,179.9,0.0,359.99,"S", phiEqpsi="Zero")
volft,kb=prep_vol(vol)

stack_data="datan.hdf"
import os
os.system("rm -f  "+stack_data)
#print  angles
nangles = len(angles)
#dropSpiderDoc("angles",angles)
ppp = []
s2x=0
s2y=0
defocus = -1
from random import random,randint
for i in xrange(nangles):
	#s2x = 4.0*randint(-1,1)
	#s2y = 4.0*randint(-1,1)
	ppp.append([angles[i][0], angles[i][1], angles[i][2], s2x, s2y])
	projo = prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], -s2x, -s2y])
	#apply CTF
	#defocus = randint(3,4)*100.0
	proj = projo.copy()
	#proj = filt_ctf(projo, defocus, 2.0, 300, 2.5, 0.1)
	#if(i == 0): st = Util.infomask(proj,None,True)
	#proj += model_gauss_noise(st[1]*3.,nx,nx)
	# Set all parameters for the new 2D image
	# three angles and two shifts to zero
	proj.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2], 's2x':s2x, 's2y':s2y})
	# CTF parameters, if defocus zero, they are undetermined
	proj.set_attr_dict({'defocus':defocus, 'amp_contrast':0.1, 'voltage':300, 'Cs':2.0, 'Pixel_size':2.5, 'B_factor':0.})
	# flags describing the status of the image (1 = true, 0 = false)
	proj.set_attr_dict({'active':1, 'ctf_applied':0})
	proj.write_image(stack_data, i)
del volft
dropSpiderDoc("params.doc",ppp)
#exit()
del ppp
nprojdata = EMUtil.get_image_count(stack_data)
snr = 1.0e20
vol1   = recons3d_4nn_ctf(stack_data, range(0,nprojdata,2), snr)
#vol1   = recons3d_4nn(stack_data, range(0,nprojdata,2))
dropImage(vol1,"v1.spi","s")

vol2   = recons3d_4nn_ctf(stack_data, range(1,nprojdata,2), snr)
#vol2   = recons3d_4nn(stack_data, range(1,nprojdata,2))

#mask3d = model_circle(nx//2-5,nx,nx,nx)

#Util.mul_img(vol1, mask3d)
#Util.mul_img(vol2, mask3d)
#del mask3d

dropImage(vol2,"v2.spi","s")
# THIS HAS TO BE INDEXED BY THE LOOP
fsc(vol1,vol2,1.0,"tdt.txt")
exit()
#vol = recons3d_wbp(stack_data, range(nprojdata),"general", const=1.0e4)
#dropImage(vol,"newvolwp.spi","s")
#vol = recons3d_4nn(stack_data, range(nprojdata))
#dropImage(vol,"newvolnp.spi","s")
#snr = 1.0
#vol = recons3d_4nn_ctf(stack_data, range(nprojdata), snr)
#dropImage(vol,"newvolp.spi","s")
exit()
mask3D = []
# do the projection matching
first_ring = 1
last_ring = 35
rstep = 1
xrng  = 1
yrng  = 1
step  = 1
dtheta= 15
# begin a refinement loop, slowly decrease dtheta inside the loop
snr = 100.0
for iter in xrange(1):
	print " ITERATION #",iter
	#proj_ali(vol, mask3D, stack_data, first_ring, last_ring, rstep, xrng, yrng, step, dtheta)

	# now have to retrieve all angles, it is only to print them!!
	cang = []
	for i in xrange(nangles):
		proj.read_image(stack_data, i)
		phi = proj.get_attr('phi')
		theta = proj.get_attr('theta')
		psi = proj.get_attr('psi')
		s2x = proj.get_attr('s2x')
		s2y = proj.get_attr('s2y')
		print  i,phi,theta,psi,s2x,s2y
		cang.append([phi, theta, psi, s2x, s2y])
	out="anglesn"
	dropSpiderDoc(out,cang)
	#calculate new and improved 3D
	stack_data = "data.hdf"
	nprojdata = EMUtil.get_image_count(stack_data)
	list_p = range(nprojdata)
	vol = recons3d_4nn_ctf(stack_data, list_p, snr, 1, "c2")
	dropImage(vol,"newvol.spi","s")
	del  vol
	exit()
	#and now the resolution

	list_p = range(0,nprojdata,2)
	vol1   = recons3d_4nn_ctf(stack_data, list_p, snr)

	list_p = range(1,nprojdata,2)
	vol2   = recons3d_4nn_ctf(stack_data, list_p, snr)

	mask3d = model_circle(nx//2-5,nx,nx,nx)

	Util.mul_img(vol1, mask3d)
	Util.mul_img(vol2, mask3d)
	del mask3d
	#dropImage(vol1,"v1.spi")

	#dropImage(vol2,"v2.spi")
	# THIS HAS TO BE INDEXED BY THE LOOP
	fsc(vol1,vol2,0.5,"tdt")

	# here figure the filtration parameters and filter vol for the next iteration

	#  here the loop should end
