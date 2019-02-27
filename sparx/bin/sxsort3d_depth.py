#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  


# This file is intentionally left blank.

def proj_ali_helical(data, refrings, numr, xrng, yrng, stepx, ynumber, psi_max=180.0, finfo=None):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj
	from math         import cos, sin, pi

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()

	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)

	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helical(data, refrings, txrng, tyrng, stepx, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber))
	iref = int(iref)
	#print  " IN ", ang, sxs, sys, mirror, iref, peak
	if iref > -1:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f\n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()
		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0

def proj_ali_helical_local(data, refrings, numr, xrng, yrng, stepx,ynumber, an, psi_max=180.0, finfo=None, yrnglocal=-1.0):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj
	from math         import cos, sin, radians

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = cos(radians(an))
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()
	
	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)

	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helical_local(data, refrings, txrng, tyrng, stepx, ant, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber), yrnglocal)

	iref = int(iref)

	if iref > -1:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write("ref phi: %9.4f\n"%(refrings[iref].get_attr("phi")))
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f \n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()

		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0\

def proj_ali_helical_90(data, refrings, numr, xrng, yrng, stepx, ynumber, psi_max=180.0, finfo=None):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()

	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)
	
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helical_90(data, refrings, txrng, tyrng, stepx, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber))
	iref = int(iref)
	#print  " IN ", ang, sxs, sys, mirror, iref, peak
	if iref > -1:
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		phi   = refrings[iref].get_attr("phi")
		theta = refrings[iref].get_attr("theta")
		psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f\n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()
		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0

def proj_ali_helical_90_local(data, refrings, numr, xrng, yrng, stepx, ynumber, an, psi_max=180.0, finfo=None, yrnglocal=-1.0):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj
	from math         import cos, sin, radians

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = cos(radians(an))
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()

	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)
	
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helical_90_local(data, refrings, txrng, tyrng, stepx, ant, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber), yrnglocal)
	iref = int(iref)
	if iref > -1:
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		phi   = refrings[iref].get_attr("phi")
		theta = refrings[iref].get_attr("theta")
		psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f\n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()
		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0

#  HELICON functions
def proj_ali_helicon_local(data, refrings, numr, xrng, yrng, stepx,ynumber, an, psi_max=180.0, finfo=None, yrnglocal=-1.0):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj
	from math         import cos, sin, radians

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = cos(radians(an))
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()

	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)
	
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helicon_local(data, refrings, txrng, tyrng, stepx, ant, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber), yrnglocal)

	iref = int(iref)

	if iref > -1:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write("ref phi: %9.4f\n"%(refrings[iref].get_attr("phi")))
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f \n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()

		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0\

def proj_ali_helicon_90_local_direct(data, refrings, xrng, yrng, \
		an, psi_max=180.0, psi_step=1.0, stepx = 1.0, stepy = 1.0, finfo=None, yrnglocal=-1.0):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from utilities    import compose_transform2, get_params_proj
	from alignment    import directaligridding
	from math         import cos, sin, radians

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	#cnx  = nx//2 + 1
	#cny  = ny//2 + 1
	ant = cos(radians(an))
	phi, theta, psi, tx, ty = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()
	#  Determine whether segment is up and down and search for psi in one orientation only.
	if psi < 180.0 :  direction = "up"
	else:             direction = "down"
	peak = -1.0e23
	iref = -1
	imn1 = sin(radians(theta))*cos(radians(phi))
	imn2 = sin(radians(theta))*sin(radians(phi))
	imn3 = cos(radians(theta))
	print('  aaaaaa  ',psi_max, psi_step, xrng, yrng, direction)
	for i in range(len(refrings)):
		if( (refrings[i][0].get_attr("n1")*imn1 + refrings[i][0].get_attr("n2")*imn2 + refrings[i][0].get_attr("n3")*imn3)>=ant ):
			print(" Matching refring  ",i,phi, theta, psi, tx, ty)
			#  directali will do fft of the input image and 180 degs rotation, if necessary.  Eventually, this would have to be pulled up.
			a, tx,ty, tp = directaligridding(data, refrings[i], psi_max, psi_step, xrng, yrng, stepx, stepy, direction)
			if(tp>peak):
				peak = tp
				iref = i
				angb = a
				sxb = tx
				syb = ty
	"""
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helicon_90_local(data, refrings, xrng, yrng, stepx, ant, psi_max, mode, numr, cnx-tx, cny-ty, int(ynumber), yrnglocal)
	"""
	if iref > -1:
		#angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		phi   = refrings[iref][0].get_attr("phi")
		theta = refrings[iref][0].get_attr("theta")
		psi   = (refrings[iref][0].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb #+ tx
		s2y   = syb #+ ty
		print("New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f" %(phi, theta, psi, s2x, s2y, peak))
		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f\n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()
		return peak, phi, theta, psi, s2x, s2y
	else:
		print("  NO PEAK")
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0

def proj_ali_helicon_90_local_direct1(data, refrings, xrng, yrng, \
		psi_max=180.0, psi_step=1.0, stepx = 1.0, stepy = 1.0, finfo=None, yrnglocal=-1.0, direction = "both"):
	"""
	  psi_max - how much psi can differ from either 90 or 270 degrees
	"""
	from utilities    import inverse_transform2, get_params_proj
	from alignment    import directaligridding1
	from math         import cos, sin, radians
	
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	#cnx  = nx//2 + 1
	#cny  = ny//2 + 1

	phi, theta, psi, tx, ty = get_params_proj(data)

	#  directali will do fft of the input image and 180 degs rotation, if necessary.  Eventually, this would have to be pulled up.
	angb, tx,ty, tp = directaligridding1(data, kb, refrings, psi_max, psi_step, xrng, yrng, stepx, stepy, direction)

	if tp > -1.0e23:
		#angb, sxb, syb, ct = inverse_transform2(ang, sxs, sys, 0)
		phi   = refrings[iref][0].get_attr("phi")
		theta = refrings[iref][0].get_attr("theta")
		psi   = (refrings[iref][0].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb #+ tx
		s2y   = syb #+ ty
		return peak, phi, theta, psi, s2x, s2y
	else:
		print("  NO PEAK")
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0

def proj_ali_helicon_90_local(data, refrings, numr, xrng, yrng, stepx, ynumber, an, psi_max=180.0, finfo=None, yrnglocal=-1.0):
	"""
	  psi_max - how much psi can differ from 90 or 270 degrees
	"""
	from alignment import search_range
	from utilities    import compose_transform2, get_params_proj
	from math         import cos, sin, pi

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = cos(an*pi/180.0)
	phi, theta, psi, sxi, syi = get_params_proj(data)
	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, tx, ty))
		finfo.flush()

	ou = numr[-3]
	sxi = round(sxi,2)
	syi = round(syi,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)
	
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helicon_90_local(data, refrings, txrng, tyrng, stepx, ant, psi_max, mode, numr, cnx-sxi, cny-syi, int(ynumber), yrnglocal)
	iref = int(iref)
	if iref > -1:
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		phi   = refrings[iref].get_attr("phi")
		theta = refrings[iref].get_attr("theta")
		psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxi
		s2y   = syb + syi

		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f\n\n" %(phi, theta, psi, s2x, s2y, peak))
			finfo.flush()
		return peak, phi, theta, psi, s2x, s2y
	else:
		return -1.0e23, 0.0, 0.0, 0.0, 0.0, 0.0
