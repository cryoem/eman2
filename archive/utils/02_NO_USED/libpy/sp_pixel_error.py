



















































































































































































from __future__ import print_function
def angle_error(ang1, ang2, delta_ang=0.0):
	'''
	This function calculates the error (variance) between two sets of angles after delta_ang (angle difference) is added to the
	first sets. When the angle difference (delta_ang) is the true difference, this function will return maximum error.
	'''
	pass#IMPORTIMPORTIMPORT from math import cos, sin, pi, radians
	
	erra = 0.0
	errb = 0.0
	delta_ang = numpy.radians( delta_ang )
	for i in range(len(ang1)):
		p2   = numpy.radians( ang2[i] )
		p2_x = numpy.cos(p2)
		p2_y = numpy.sin(p2)
		p1   = numpy.radians( ang1[i] )
		p1_x = numpy.cos(p1)
		p1_y = numpy.sin(p1)

		erra += p2_x*p1_x+p2_y*p1_y
		errb += p2_y*p1_x-p2_x*p1_y

	return erra*numpy.cos(delta_ang) + errb*numpy.sin(delta_ang)







































































def align_diff(data1, data2=None, suffix="_ideal"):
	
	'''
	This function determines the relative angle, shifts and mirrorness between
	two list of data
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_params2D
	
	nima = len(data1)

	if data2 != None: 
		nima2 = len(data2)
		if nima2 != nima:
			sp_global_def.sxprint("Error: Number of images don't agree!")
			return 0.0, 0.0, 0.0, 0
		else:
			del nima2

	# Read the alignment parameters and determine the relative mirrorness
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		alpha1, sx1, sy1, mirror1, scale1 = sp_utilities.get_params2D(data1[i])
		if data2 != None:
			alpha2, sx2, sy2, mirror2, scale2 = sp_utilities.get_params2D(data2[i])
		else:
			alpha2, sx2, sy2, mirror2, scale2 = sp_utilities.get_params2D(data1[i], "xform.align2d"+suffix)
		ali_params1.extend([alpha1, sx1, sy1, mirror1])
		ali_params2.extend([alpha2, sx2, sy2, mirror2])

	return align_diff_params(ali_params1, ali_params2)


def align_diff_textfile(textfile1, textfile2):
	
	'''
	This function (2D) determines the relative angle, shifts and mirrorness between
	the two textfile of alignment parameters
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import read_text_row
	
	ali1 = sp_utilities.read_text_row(textfile1, "", "")
	ali2 = sp_utilities.read_text_row(textfile2, "", "")

	nima = len(ali1)
	nima2 = len(ali2)
	if nima2 != nima:
		sp_global_def.sxprint("Error: Number of images don't agree!")
		return 0.0, 0.0, 0.0, 0
	else:
		del nima2

	# Read the alignment parameters and determine the relative mirrorness
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		ali_params1.extend(ali1[i][0:4])
		ali_params2.extend(ali2[i][0:4])

	return align_diff_params(ali_params1, ali_params2)


def ave_ali_err(data1, data2=None, r=25, suffix="_ideal"):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two lists of data. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_params2D, combine_params2
	pass#IMPORTIMPORTIMPORT from math import sqrt, sin, pi
	
	# Determine relative angle, shifts and mirror
	alphai, sxi, syi, mirror = align_diff(data1, data2, suffix)

	# Determine the average pixel error
	err = 0.0
	nima = len(data1)
	mirror_same = 0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1, scale1 = sp_utilities.get_params2D(data1[i])
		if data2 != None:
			alpha2, sx2, sy2, mirror2, scale2 = sp_utilities.get_params2D(data2[i])
		else:
			alpha2, sx2, sy2, mirror2, scale2 = sp_utilities.get_params2D(data1[i], "xform.align2d"+suffix)
		
		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = sp_utilities.combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)
	
	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def ave_ali_err_params(ali_params1, ali_params2, r=25):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import combine_params2
	pass#IMPORTIMPORTIMPORT from math import sqrt, sin, pi

	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	mirror_same = 0
	err = 0.0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]

		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = sp_utilities.combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)

	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def ave_ali_err_textfile(textfile1, textfile2, r=25):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import combine_params2
	pass#IMPORTIMPORTIMPORT from math import sqrt, sin, pi
	pass#IMPORTIMPORTIMPORT from sp_utilities import read_text_row
	
	ali1 = sp_utilities.read_text_row(textfile1, "", "")
	ali2 = sp_utilities.read_text_row(textfile2, "", "")

	nima = len(ali1)
	nima2 = len(ali2)
	if nima2 != nima:
		sp_global_def.sxprint("Error: Number of images don't agree!")
		return 0.0, 0.0, 0.0, 0, 0.0, 0.0
	else:
		del nima2

	# Read the alignment parameters
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		ali_params1.extend(ali1[i][0:4])
		ali_params2.extend(ali2[i][0:4])

	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	mirror_same = 0
	err = 0.0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]
		
		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = sp_utilities.combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)
	
	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def multi_align_diff_params(ali_params, verbose=0):
	"""
	Calculate the mirror consistency and pixel error between different runs of alignment
	The input is in the following format (here n is the number of alignments done):
		[[alpha1_r1, sx1_s1, sy1_r1, mirror1_r1, alpha2_r1, sx2_r1, sy2_r1, mirror2_r1, ...],
		 [alpha1_r2, sx1_s2, sy1_r2, mirror1_r2, alpha2_r2, sx2_r2, sy2_r2, mirror2_r2, ...],
		 [alpha1_r3, sx1_s3, sy1_r3, mirror1_r3, alpha2_r3, sx2_r3, sy2_r3, mirror2_r3, ...],
		 ...
		 [alpha1_rn, sx1_sn, sy1_rn, mirror1_rn, alpha2_rn, sx2_rn, sy2_rn, mirror2_rn, ...]]
	The output is in the following format (here k=n*(n+1)/2):
		[[pixel_error_1, mirror_consistency_1, i_1, j_1, alpha_1, sx_1, sy_1, mirror_1],
		 [pixel_error_2, mirror_consistency_2, i_2, j_2, alpha_2, sx_2, sy_2, mirror_2],
		 [pixel_error_3, mirror_consistency_3, i_3, j_3, alpha_3, sx_3, sy_3, mirror_3],
		...
		 [pixel_error_k, mirror_consistency_k, i_k, j_k, alpha_k, sx_k, sy_k, mirror_k]]
	"""
	num_ali = len(ali_params)
	multi_align_results = []
	for i in range(num_ali-1):
		for j in range(i+1, num_ali):
			alpha, sx, sy, mirror, stab_mirror, pixel_error = ave_ali_err_params(ali_params[i], ali_params[j])
			if verbose == 1:
				sp_global_def.sxprint("Between trial %d and %d: mirror stability = %6.3f   pixel error = %6.3f"%(i, j, stab_mirror, pixel_error))
			multi_align_results.append([pixel_error, stab_mirror, i, j, alpha, sx, sy, mirror])
	return multi_align_results
	

def calc_connect_list(multi_align_results, pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8):
	"""
	Generate the connection list from the multi_align_results, which generally comes from multi_align_diff_params()
	The connection list will have the following format:
		[[1, 2, 5], [4, 6], [0, 7]]
	You will also get the size of the largest connection in the list.
	"""
	pass#IMPORTIMPORTIMPORT import sets
	
	k = len(multi_align_results)
	multi_align_results.sort()
	connect_list = []
	for i in range(k):
		if multi_align_results[i][0] <= pixel_error_threshold:
			if multi_align_results[i][1] >= mirror_consistency_threshold: 
				connect_list.append([multi_align_results[i][2], multi_align_results[i][3]])
		else:	break
	to_break = True
	while to_break:
		l = len(connect_list)
		to_break = False
		for i in range(l-1):
			for j in range(i+1, l):
				set1 = set(connect_list[i])
				set2 = set(connect_list[j])
				if list(set1.intersection(set2)) != []:
					connect_list[i] = list(set1.union(set2))
					del connect_list[j]
					to_break = True
					break
			if to_break: break
	max_connect = 0
	for l in connect_list: max_connect = max(max_connect, len(l))
	return connect_list, max_connect


def ali_stable_list(ali_params1, ali_params2, pixel_error_threshold, r=25):
	'''
	This function first determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It then determines whether each image is
	stable or not and return this information as an int list. (1 is stable and 0 is unstable)
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import combine_params2
	pass#IMPORTIMPORTIMPORT from math import sqrt, sin, pi
	
	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	ali_list = []
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]
		if abs(mirror1-mirror2) == mirror:
			alpha12, sx12, sy12, mirror12 = sp_utilities.combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			if pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r) < pixel_error_threshold: ali_list.append(1)
			else: ali_list.append(0)
		else: ali_list.append(0)

	return ali_list


























































































































































































































































def ave2dtransform(args, data, return_avg_pixel_error=False):

	pass#IMPORTIMPORTIMPORT from math import pi, sin, cos, radians, degrees

	ali_params = data[0]
	d = data[1]
	L = len(ali_params)
	N = len(ali_params[0])/4

	args_list= [0.0]*(L*3)
	for i in range(L*3-3):  args_list[i] = args[i]

	cosa = [0.0]*L
	sina = [0.0]*L
	for i in range(L):
		cosa[i] = numpy.cos(numpy.radians(args_list[i*3]))
		sina[i] = numpy.sin(numpy.radians(args_list[i*3]))
	sqr_pixel_error = [0.0]*N
	ave_params = []
	for i in range(N):
		sum_cosa = 0.0
		sum_sina = 0.0
		sx = [0.0]*L
		sy = [0.0]*L
		for j in range(L):
			if int(ali_params[j][i*4+3]) == 0:
				sum_cosa += numpy.cos(numpy.radians(args_list[j*3]+ali_params[j][i*4]))
				sum_sina += numpy.sin(numpy.radians(args_list[j*3]+ali_params[j][i*4]))
				sx[j] =  args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] + ali_params[j][i*4+2]*sina[j]
				sy[j] =  args_list[j*3+2] - ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
			else:
				sum_cosa += numpy.cos(numpy.radians(-args_list[j*3]+ali_params[j][i*4]))
				sum_sina += numpy.sin(numpy.radians(-args_list[j*3]+ali_params[j][i*4]))
				sx[j] = -args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] - ali_params[j][i*4+2]*sina[j]
				sy[j] =  args_list[j*3+2] + ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
		sqrtP = numpy.sqrt(sum_cosa**2+sum_sina**2)
		sqr_pixel_error[i] = max( 0.0, d*d/4.*(1-sqrtP/L) + sqerr(sx) + sqerr(sy) )
		# Get ave transform params
		H = EMAN2_cppwrap.Transform({"type":"2D"})
		H.set_matrix([sum_cosa/sqrtP, sum_sina/sqrtP, 0.0, sum(sx)/L, -sum_sina/sqrtP, sum_cosa/sqrtP, 0.0, sum(sy)/L, 0.0, 0.0, 1.0, 0.0])
		dd = H.get_params("2D")
		#  We are using here mirror of the LAST SET.
		H = EMAN2_cppwrap.Transform({"type":"2D","alpha":dd[ "alpha" ],"tx":dd[ "tx" ],"ty": dd[ "ty" ],"mirror":int(ali_params[L-1][i*4+3]),"scale":1.0})
		dd = H.get_params("2D")
		ave_params.append([dd[ "alpha" ], dd[ "tx" ], dd[ "ty" ], dd[ "mirror" ]])
	# Warning: Whatever I return here is squared pixel error, this is for the easy expression of derivative
	# Don't forget to square root it after getting the value
	if return_avg_pixel_error:
		return sum(sqr_pixel_error)/N
	else:
		return sqr_pixel_error, ave_params

def rotate_angleset_to_match(agls1, agls2):
	"""
	  Finds rotation between two sets of angles, agls2 is the template
	  It will also establish whether mirror is required
	  Rotation is applied to agsl1 and the set of rotated angles is returned
	  Rotation itself is not returned.
	  Makes sense only for no symmetry
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities    import rotation_between_anglesets
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import rotate_params

	t1 = sp_utilities.rotation_between_anglesets(agls1, agls2)

	return sp_fundamentals.rotate_params(agls1,[-t1[2],-t1[1],-t1[0]])

def ordersegments(infilaments, ptclcoords):
	'''
	Input:
	
	stack: Input stack of images whose headers contain filament membership information and particle coordinates in the original micrograph (stored under attribute ptcl_source_coord).
	filament_attr: Attribute under which filament membership ID is stored.
	It is assumed the prtl coords are nonnegative
	
	Output: 
	
	Returns a list of lists, where each inner list consists of IDs of segments windowed from
	a single filament ordered according to their relative positions on the filament.
	
	'''

	def orderbymodule(xxp,yyp):
		pass#IMPORTIMPORTIMPORT from math import atan,sin,cos,pi, atan2
		pass#IMPORTIMPORTIMPORT from sp_statistics import linreg
		nq = len(xxp)
		xs = sum(xxp)/nq
		ys = sum(yyp)/nq
		xp = [0.0]*nq
		yp = [0.0]*nq
		for i in range(nq):
			xp[i] = xxp[i] - xs
			yp[i] = yyp[i] - ys
		try:
			a,b = sp_statistics.linreg(xp,yp)
			alpha = numpy.pi/4-math.atan(a)
		except:
			a,b = sp_statistics.linreg([(xp[i]-yp[i]) for i in range(nq)], [(xp[i]+yp[i]) for i in range(nq)])
			alpha = math.atan(a)
			#print "except"

		cs = numpy.cos(alpha)
		ss = numpy.sin(alpha)
		qm = 1.e23
		dd = [[0.0, 0] for i in range(nq)]
		for i in range(nq):
			xt =  cs*xp[i] - ss*yp[i]
			yt =  ss*xp[i] + cs*yp[i]
			xp[i] = xt; yp[i] = yt
		xs = min(xp)
		ys = min(yp)
		for i in range(nq):
			dd[i] = [(xp[i]-xs)**2+(yp[i]-ys)**2, i]
		dd.sort()
		return [dd[i][1] for i in range(nq)]

	allfilaments = [None]*len(infilaments)
	for i in range(len(infilaments)):
		allfilaments[i] = [infilaments[i],i]
	allfilaments.sort()
	filaments = []
	current = allfilaments[0][0]
	temp = [allfilaments[0][1]]
	for i in range(1,len(allfilaments)):
		if( allfilaments[i][0] == current ):
			temp.extend([allfilaments[i][1]])
		else:
			filaments.append(temp)
			current = allfilaments[i][0]
			temp = [allfilaments[i][1]]
	filaments.append(temp)

	del allfilaments,temp

	nfil = len(filaments)

	for i in range(nfil):
		nsegs = len(filaments[i])
		if(nsegs > 1):
			ord = orderbymodule([ptclcoords[filaments[i][ii]][0] for ii in range(nsegs)],[ptclcoords[filaments[i][ii]][1] for ii in range(nsegs)])
			filaments[i] = [filaments[i][ord[ii]] for ii in range(nsegs)]
			# To preserve the original order check indexes and invert if it appears to be inverted
			if(filaments[i][0] > filaments[i][-1]):
				for k in range(nsegs//2):
					temp = filaments[i][k]
					filaments[i][k] = filaments[i][nsegs-1-k]
					filaments[i][nsegs-1-k]=temp
	return filaments	


def mapcoords(x, y, r, nx, ny):
	pass#IMPORTIMPORTIMPORT from math 			import ceil, floor
	pass#IMPORTIMPORTIMPORT from sp_utilities 	import get_dist
	pass#IMPORTIMPORTIMPORT import sys
	'''
	Input:
	
	(x,y): Coordinate in old image. 
	r: ratio by which old image is resampled by. If r < 1, then pixel size of resampled image is original pixel size divided by r.
	nx, ny: dimensions of old image
	
	Assumes coordinates are positive and run from 0 to nx-1 and 0 to ny-1, where nx and ny
	are the x and y dimensions of the micrograph respectively.
	
	Output:
	
	x'',y'':	The pixel coordinate in the resampled image where
	
					(x',y') = (Util.round(x''/r), Util.round(y''/r))
					
				and (x',y') is the closest point to (x,y) over all other points (a,b) in the
				old image where (a,b)=  (Util.round(a''/r), Util.round(b''/r)) for some
				pixel coordinate (a'', b'') in resampled image.
	'''	
	
	# Neighborhood of (x,y) in old image which contains a point (x',y') such that
	# (x',y') = (Util.round(x''/r), Util.round(y''/r)) for some (x'', y'') in resampled image
	if r > 1:
		nbrhd = 1
	else:
		nbrhd = int(numpy.ceil(1.0/r))+1
			
	allxnew = []
	allynew = []
	
	for i in range(-nbrhd, nbrhd+1):
		xold = EMAN2_cppwrap.Util.round(x + i)
		if xold < 0 or xold >= nx:
			continue
		# See if there is xnew in new image where xold == int(Util.round(xnew/r))
		# If there is such a xnew, then xold == int(Util.round(xnew/r)) implies r*(xold-0.5) <= xnew < r*(xold+0.5)
		lxnew = int(numpy.floor(r*(xold - 0.5)))
		uxnew = int(numpy.ceil(r*(xold + 0.5))) 
		for xn in range(lxnew, uxnew + 1):
			if xold == EMAN2_cppwrap.Util.round(xn/r):
				allxnew.append(xn)
				
	for j in range(-nbrhd, nbrhd+1):
		yold = EMAN2_cppwrap.Util.round(y + j)
		if yold < 0 or yold >= ny:
			continue
		lynew = int(numpy.floor(r*(yold - 0.5)))
		uynew = int(numpy.ceil(r*(yold + 0.5)))
		for yn in range(lynew, uynew + 1):
			if yold == EMAN2_cppwrap.Util.round(yn/r):
				allynew.append(yn)
				
	if len(allxnew) == 0 or len(allynew) == 0:
		sp_global_def.ERROR("Could not find mapping")
	
	mindist = -1
	minxnew = -1
	minynew = -1
	
	for xnew in allxnew:
		for ynew in allynew:
			xold = EMAN2_cppwrap.Util.round(xnew/r)
			yold = EMAN2_cppwrap.Util.round(ynew/r)
			dst = sp_utilities.get_dist([x,y],[xold,yold])
			if dst > mindist:
				mindist = dst
				minxnew = int(xnew)
				minynew = int(ynew)
					
	return minxnew, minynew

def consistency_params(stack, dphi, dp, pixel_size, phithr=2.5, ythr=1.5, THR=3):
	'''
		stack        - contains coding of filaments and coordinates of segments ptcl_source_coord
		fname_params - parameters whose consistency is tested
	'''
	pass#IMPORTIMPORTIMPORT from sp_utilities import read_text_row, write_text_row, get_dist
	pass#IMPORTIMPORTIMPORT from sp_applications import ordersegments
	pass#IMPORTIMPORTIMPORT from sp_pixel_error import angle_diff

	filaments = ordersegments(stack)
	ptclcoords = EMAN2_cppwrap.EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
	params     = EMAN2_cppwrap.EMUtil.get_all_attributes(stack, 'xform.projection')
	for i in range(len(params)):
		d = params[i].get_params("spider")
		params[i] = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"] ]

	N = len(filaments)
	sp_global_def.sxprint("N: ", N)
	totsegs    = 0
	totpsicons = 0
	totphicons = 0
	tot_nopred = 0
	allphier = []
	imi = 0
	for mic in filaments:
		#for kkk in xrange(1):
		#mic = filaments[0]
		imi += 1
		#mic = mic[:4]
		if len(mic) < THR:
			#print "less than threshold"
			allphier.append([[mic[0],mic[-1],0],[]])
		else:

			totsegs += len(mic)

			a90  = [] # segments in this filament with psi ~ 90
			a270 = [] # segments in this filament with psi ~ 270

			for iseg in mic:
				if abs(params[iseg][2] - 90.0) <45.0: a90.append(iseg)
				else:                                 a270.append(iseg)

			if (len(a90) == len(a270)):
				# these cannot be predicted
				tot_nopred += len(mic)
				continue

			if( len(a90) > len(a270) ):
				thetapsi1 =  a90
				flip =  1
			else:
				thetapsi1 = a270
				flip = -1

			ns = len(thetapsi1)
			# given phi = phig
			phig = [0.0]*ns
			for j in range(ns): phig[j] = params[thetapsi1[j]][0]
			totpsicons += ns
			distances = [0.0]*ns
			for i in range(1,ns):  distances[i] = sp_utilities.get_dist( ptclcoords[thetapsi1[i]], ptclcoords[thetapsi1[0]] )
			ganger = [0.0]*ns
			terr = 1.e23
			for idir in range(-1,2,2):
				phierr = []
				#  get phi's
				ddphi = pixel_size/dp*idir*dphi
				phis = [0.0]*ns
				#print "  MIC  ",mic
				for i in range(ns):
					yy = distances[i] + params[thetapsi1[i]][4]
					phis[i] = (yy*ddphi)%360.0
					#print " %7.3f   %7.3f   %7.3f  %7.3f   "%(ycoords[i],yy,phis[i], phig[i]),params[thetapsi1[i]]

				# find the overall angle
				angdif = angle_diff(phis,phig)
				#print " angdif ",angdif
				lerr = 0.0
				for i in range(ns):
					anger = (phis[i]+angdif - phig[i] + 360.0)%360.0
					if( anger > 180.0 ): anger -= 360.0
					lerr += abs(anger)
					phierr.append(anger)
					#print  " %7.3f   %7.3f   %7.3f"%((phis[i]+angdif+360.0)%360.0 , phig[i],anger)
				if(lerr < terr):
					terr = lerr
					for j in range(ns):  ganger[j] = phierr[j]
			allphier.append([[mic[0], mic[-1], flip, terr/ns], ganger])

	sp_global_def.sxprint("number of segments belonging to filaments from which at least %i segments were windowed: "%THR, totsegs)
	sp_global_def.sxprint("number of segments oriented 50/50 wrt psi (and therefore could not be predicted):       ", tot_nopred)
	sp_global_def.sxprint("segments whose psi agreed with the majority of segments in its filament:                ", totpsicons)
	return  allphier

def getnewhelixcoords(hcoordsname, outdir, ratio,nx,ny, newpref="resampled_", boxsize=-1):
	"""
	Input
	
		helixcoordsfile: Full path name of file with coordinates of boxed helices
		
		outdir: Full path name of directory in which to put new helix coordinates file.
		
		ratio: factor by which new image (micrograph) is resampled from old
		
		nx, ny: dimensions of old image (micrograph)
		
		newpref: prefix for attaching to fname to get name of new helix coordinates file
	
	Output:
		Returns full path name of file containing new box coordinates
	"""
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from sp_utilities 		import read_text_row
	pass#IMPORTIMPORTIMPORT from sp_pixel_error	import mapcoords
	
	fname = (hcoordsname.split('/'))[-1] # name of old coordinates files minus the path
	newhcoordsname = os.path.join(outdir , newpref+fname) # full path name of new coordinates file to be created
	f = open(newhcoordsname, 'w')
	coords = sp_utilities.read_text_row(hcoordsname) # old coordinates
	ncoords = len(coords)
	newcoords=[]
	w = coords[0][2]
	new_w = boxsize
	if new_w < 0:
		new_w = w*ratio
	for i in range(ncoords):
		xold = coords[i][0] + w/2
		yold = coords[i][1] + w/2
		xnew, ynew = mapcoords(xold,yold,ratio,nx,ny)
		s = '%d\t%d\t%d\t%d\t%d\n'%(xnew-new_w/2,ynew-new_w/2, new_w, new_w, coords[i][4])
		f.write(s)
	return newhcoordsname	

def helical_params_err(params1, params2, fil_list):
	'''
	Input:
		params1: First set of projection orientation parameters
		params2: Second set of projection orientation parameters
		fil_list: List of filament IDs, where i-th element on list is filament of i-th image.
				  Assume fil_list is ordered, i.e., i-th element of fil_list is the filament ID
				  of the i-th image in the stack.
	
	Output:
		The phi angle difference is computed between params1 and params2 between those parameters
		whose psi agree, and the angle is applied to params2.
		
		For each filament, the program computes the phi error averaged over the number 
		of segments between params1 and the aligned params2.
		
		Returns a list of lists, where each inner list is [fil_ID, avg_phi_err], where fil_ID
		is filament name, and avg_phi_err is the average phi error for the filament.
	'''
	pass#IMPORTIMPORTIMPORT from sp_pixel_error import angle_diff
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f
	nima = len(params1)
	# Identify those where psi agrees
	phi1 = []
	phi2 = []
	fil_psisame = []
	pos = []
	pref = []
	for i in range(nima):
		if abs(params1[i][2] - params2[i][2]) < 90:
			phi1.append(params1[i][0])
			phi2.append(params2[i][0])
			fil_psisame.append(fil_list[i])
			pos.append(i)
			pref.append(params2[i])

	fag = len(fil_psisame)

	tflip = EMAN2_cppwrap.Transform({"type":"spider","theta":180.0})
	# Identify those where psi agrees
	phi1 = []
	phi2 = []
	fil_psisame = []
	pos  = []
	prot = []
	pref = []
	for i in range(nima):
		t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":params1[i][0],"theta":params1[i][1],"psi":params1[i][2]})
		t2.set_trans( EMAN2_cppwrap.Vec2f( -params1[i][3], -params1[i][4] ) )
		t2 = t2*tflip
		d = t2.get_params("spider")
		p1r = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]]
		if abs(p1r[2] - params2[i][2]) < 90:
			phi1.append(p1r[0])
			phi2.append(params2[i][0])
			fil_psisame.append(fil_list[i])
			pos.append(i)
			prot.append(p1r)
			pref.append(params2[i])
			
	nima2 = len(fil_psisame)
	if(fag > nima2):
		phi1 = []
		phi2 = []
		fil_psisame = []
		pos = []
		prot = []
		pref = []
		for i in range(nima):
			if abs(params1[i][2] - params2[i][2]) < 90:
				phi1.append(params1[i][0])
				phi2.append(params2[i][0])
				fil_psisame.append(fil_list[i])
				pos.append(i)
				prot.append(params1[i])
				pref.append(params2[i])
		nima2 = len(fil_psisame)
	else:
		sp_global_def.sxprint("better agreement afer upside-down rotation")
		
	# agls1psi and agls2psi agree in psi. 
	sp_global_def.sxprint("Number of images which agree on psi between params1 and params2: ",nima2)
	sp_global_def.sxprint("Percentage of total number images being compared: ", nima2/float(len(params1))*100)
	
	angdif = angle_diff(phi1, phi2)
	sp_global_def.sxprint("angdif: ", angdif)
	
	phierr_byfil = []

	start_fil = None
	i = 0
	while (i < nima2):
		ibeg = i
		iend = i
		cur_fil = fil_psisame[ibeg]
		for k in range(ibeg+1,nima2):
			if(cur_fil != fil_psisame[k]):
				iend = k
				break
		if( ibeg == iend ):  iend = nima2-1
		sum_phierr = 0.
		for k in range(ibeg,iend):
			phidf = abs(phi2[k] - (phi1[k] + angdif)%360.)
			phidf = min(phidf, abs(360.0-phidf))
			sum_phierr += phidf
			prot[k][0] =  (phi1[k] + angdif)%360.
		avg_phierr = sum_phierr/(iend-ibeg+1)
		phierr_byfil.append([avg_phierr, cur_fil,ibeg,iend,pos[ibeg],pos[iend]])
		i = iend+1

	return phierr_byfil,prot,pref












































































































































































































































































































































































































































































































































