'''1
def golden(func, args=(), brack=None, tol=1.e-4, full_output=0):
	""" Given a function of one-variable and a possible bracketing interval,
	return the minimum of the function isolated to a fractional precision of
	tol. A bracketing interval is a triple (a,b,c) where (a<b<c) and
	func(b) < func(a),func(c).  If bracket is two numbers then they are
	assumed to be a starting interval for a downhill bracket search
	(see bracketing)

	Uses analog of bisection method to decrease the bracketed interval.
	"""
	pass#IMPORTIMPORTIMPORT from utilities import bracketing
	if brack is None:
		xa,xb,xc,fa,fb,fc,funcalls = bracketing(func, args=args)
	elif len(brack) == 2:
		xa,xb,xc,fa,fb,fc,funcalls = bracketing(func, xa=brack[0], xb=brack[1], args=args)
	elif len(brack) == 3:
		xa,xb,xc = brack
		if (xa > xc):  # swap so xa < xc can be assumed
			dum = xa; xa=xc; xc=dum
		assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
		fa = apply(func, (xa,)+args)
		fb = apply(func, (xb,)+args)
		fc = apply(func, (xc,)+args)
		assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
		funcalls = 3
	else:
		raise ValueError, "Bracketing interval must be length 2 or 3 sequence."

	_gR = 0.61803399
	_gC = 1.0-_gR
	x3 = xc
	x0 = xa
	if (abs(xc-xb) > abs(xb-xa)):
		x1 = xb
		x2 = xb + _gC*(xc-xb)
	else:
		x2 = xb
		x1 = xb - _gC*(xb-xa)
	f1 = apply(func, (x1,)+args)
	f2 = apply(func, (x2,)+args)
	funcalls += 2
	while (abs(x3-x0) > tol*(abs(x1)+abs(x2))):
		if (f2 < f1):
			x0 = x1; x1 = x2; x2 = _gR*x1 + _gC*x3
			f1 = f2; f2 = apply(func, (x2,)+args)
		else:
			x3 = x2; x2 = x1; x1 = _gR*x2 + _gC*x0
			f2 = f1; f1 = apply(func, (x1,)+args)
		funcalls += 1
	if (f1 < f2):
		xmin = x1
		fval = f1
	else:
		xmin = x2
		fval = f2
	if full_output:
		return xmin, fval, funcalls
	else:
		return xmin


def bracketing(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
	"""Given a function and distinct initial points, search in the downhill
	direction (as defined by the initital points) and return new points
	xa, xb, xc that bracket the minimum of the function:
	f(xa) > f(xb) < f(xc)
	"""
	_gold = 1.618034
	_verysmall_num = 1e-21
	fa = apply(func, (xa,)+args)
	fb = apply(func, (xb,)+args)
	if (fa < fb):			   # Switch so fa > fb
	    dum = xa; xa = xb; xb = dum
	    dum = fa; fa = fb; fb = dum
	xc = xb + _gold*(xb-xa)
	fc = apply(func, (xc,)+args)
	funcalls = 3
	iter = 0
	while (fc < fb):
	    tmp1 = (xb - xa)*(fb-fc)
	    tmp2 = (xb - xc)*(fb-fa)
	    val = tmp2-tmp1
	    if abs(val) < _verysmall_num:
		denom = 2.0*_verysmall_num
	    else:
		denom = 2.0*val
	    w = xb - ((xb-xc)*tmp2-(xb-xa)*tmp1)/denom
	    wlim = xb + grow_limit*(xc-xb)
	    if iter > maxiter:
		raise RuntimeError, "Too many iterations."
	    iter += 1
	    if (w-xc)*(xb-w) > 0.0:
		fw = apply(func, (w,)+args)
		funcalls += 1
		if (fw < fc):
		    xa = xb; xb=w; fa=fb; fb=fw
		    return xa, xb, xc, fa, fb, fc, funcalls
		elif (fw > fb):
		    xc = w; fc=fw
		    return xa, xb, xc, fa, fb, fc, funcalls
		w = xc + _gold*(xc-xb)
		fw = apply(func, (w,)+args)
		funcalls += 1
	    elif (w-wlim)*(wlim-xc) >= 0.0:
		w = wlim
		fw = apply(func, (w,)+args)
		funcalls += 1
	    elif (w-wlim)*(xc-w) > 0.0:
		fw = apply(func, (w,)+args)
		funcalls += 1
		if (fw < fc):
		    xb=xc; xc=w; w=xc+_gold*(xc-xb)
		    fb=fc; fc=fw; fw=apply(func, (w,)+args)
		    funcalls += 1
	    else:
		w = xc + _gold*(xc-xb)
		fw = apply(func, (w,)+args)
		funcalls += 1
	    xa=xb; xb=xc; xc=w
	    fa=fb; fb=fc; fc=fw
	return xa, xb, xc, fa, fb, fc, funcalls
'''
"""2
	t1.printme()
	print(" ")
	t2.printme()
	print(" ")
	tt.printme()
	"""
'''3
	elif(symmetry_string[0]  == "o"):
		pass#IMPORTIMPORTIMPORT from EMAN2 import parsesym
		if(method.lower() == "s"):  met = "saff"
		elif(method.lower() == "p"): met = "even"
		if(theta2 == 180.0):  inc_mirror = 1
		else:  inc_mirror = 0
		tt = parsesym(symmetry)
		z = tt.gen_orientations(met,{"delta":delta,"inc_mirror":inc_mirror})
		angles = []
		if( phiEqpsi == "Minus" ):
			for q in z:
				q = q.get_params("spider")
				angles.append([q["phi"], q["theta"],-q["phi"]])
		else:
			for q in z:
				q = q.get_params("spider")
				angles.append([q["phi"], q["theta"],0.0])
	else :
		
		# This is very close to the Saff even_angles routine on the asymmetric unit;
		# the only parameters used are symmetry and delta
		# The formulae are given in the Transform Class Paper
		# The symmetric unit 		nVec=[]; # x,y,z triples
		# is defined by three points b,c, v of Fig 2 of the paper
		# b is (0,0,1)
		# c is (sin(thetac),0,cos(thetac))
		# a is (sin(thetac)cos(Omega),sin(thetac)cos(Omega),cos(thetac))
		# f is the normalized sum of all 3

		# The possible symmetries are in list_syms
		# The symmetry determines thetac and Omega
		# The spherical area is Omega - pi/3;
		#  should be equal to 4 *pi/(3*# Faces)
		#
		# symmetry ='tet';   delta  = 6;

		scrunch = 0.9  # closeness factor to eliminate oversampling corners
		#nVec=[]       # x,y,z triples

		piOver = pi/180.0
		Count=0   # used to count the number of angles

		if (symmetryLower[0:3] =="tet"):  m=3.0; fudge=0.9 # fudge is a factor used to adjust phi steps
		elif (symmetryLower[0:3] =="oct"):  m=4.0; fudge=0.8
		elif (symmetryLower[0:3] =="ico"):  m=5.0; fudge=0.95
		else: ERROR("allowable symmetries are cn, dn, tet, oct, icos","even_angles",1)

		n=3.0
		OmegaR = 2.0*pi/m; cosOmega= cos(OmegaR)
		Edges  = 2.0*m*n/(2.0*(m+n)-m*n)
		Faces  = 2*Edges/n
		Area   = 4*pi/Faces/3.0; # also equals  2*pi/3 + Omega
		costhetac = cosOmega/(1-cosOmega)
		deltaRad= delta*pi/180
		NumPoints = int(Area/(deltaRad*deltaRad))
		fheight = 1/sqrt(3)/ (tan(OmegaR/2.0))

		z0      = costhetac  # initialize loop
		z       = z0
		phi     = 0
		Deltaz  = (1-costhetac)/(NumPoints-1)

		#[1, phi,180.0*acos(z)/pi,0.]
		anglesLast = [phi,180.0*acos(z)/pi,0.]
		angles.append(anglesLast)
		nLast=  [ sin(acos(z))*cos(phi*piOver) ,  sin(acos(z))*sin(phi*piOver) , z]
		nVec = []
		nVec.append(nLast)

		Count +=1

		for k in xrange(1,(NumPoints-1)):
			z=z0 + Deltaz*k  # Is it higher than fhat or lower
			r= sqrt(1-z*z)
			if (z > fheight): phiRmax= OmegaR/2.0
			if (z<= fheight):
				thetaR   = acos(z);
				cosStuff = (cos(thetaR)/sin(thetaR))*sqrt(1. - 2 *cosOmega);
				phiMax   =  180.0*( OmegaR - acos(cosStuff))/pi
			angleJump = fudge* delta/r
			phi = (phi + angleJump)%(phiMax)
			anglesNew = [phi,180.0*acos(z)/pi,0.];
			nNew = [ sin(acos(z))*cos(phi*piOver) ,  sin(acos(z))*sin(phi*piOver) , z]
			diffangleVec = [acos(nNew[0]*nVec[k][0] +   nNew[1]*nVec[k][1] +    nNew[2]*nVec[k][2] ) for k in xrange(Count)]
			diffMin = min(diffangleVec)
			if (diffMin>angleJump*piOver *scrunch):
				Count +=1
				angles.append(anglesNew)
				nVec.append(nNew)
				#[Count, phi,180*acos(z)/pi,0.]
			anglesLast = anglesNew
			nLast=nNew

		angles.append( [0.0, 0.0, 0.0] )
		nLast=  [ 0., 0. , 1.]
		nVec.append(nLast)
		if(theta2 == 180.0):   angles.append( [0.0, 180.0, 0.0] )

		angles.reverse()
		if(phiEqpsi == "Minus"):
			for i in xrange(len(angles)):  angles[i][2] = (720.0-angles[i][0])%360.0
		#print(Count,NumPoints)

		#		look at the distribution
		#		Count =len(angles); piOver= pi/180.0;
		#		phiVec    =  [ angles[k][0] for k in range(Count)] ;
		#		thetaVec  =  [ angles[k][1] for k in range(Count)] ;
		#		xVec = [sin(piOver * angles[k][1]) * cos(piOver * angles[k][0]) for k in range(Count) ]
		#		yVec = [sin(piOver * angles[k][1])* sin(piOver * angles[k][0]) for k in range(Count) ]
		#		zVec = [cos(piOver *  angles[k][1]) for k in range(Count) ]
		#		pylab.plot(yVec,zVec,'.'); pylab.show()
		'''
'''4
def reduce_array_to_root(data, myid, main_node = 0, comm = -1):
	pass#IMPORTIMPORTIMPORT from numpy import array, shape, reshape
	pass#IMPORTIMPORTIMPORT from mpi import MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, mpi_reduce, mpi_barrier

	if comm == -1:  comm = MPI_COMM_WORLD
	n = shape(data)
	ntot = 1
	for i in xrange(len(n)):  ntot *= n[i]
        count = 500000
        array1d = reshape(data, (ntot,))
        ntime = (ntot-1) /count + 1
        for i in xrange(ntime):
        	block_begin = i*count
        	block_end   = i*count + count
        	if block_end > ntot: block_end = ntot
        	block_size  = block_end - block_begin
        	tmpsum = mpi_reduce(array1d[block_begin], block_size, MPI_FLOAT, MPI_SUM, main_node, comm)
		mpi_barrier(comm)
        	if myid == main_node:
        		array1d[block_begin:block_end] = tmpsum[0:block_size]
'''
'''5
def bcast_EMData_to_all(img, myid, main_node = 0, comm = -1):

	# Comment by Zhengfan Yang on 01/05/10
	#
	# Notice:
	# (1) one should use this new version of broadcasting EMData in the following way:
	# 	img = bcast_EMData_to_all(img, myid, main_node, comm)
	# instead of
	# 	bcast_EMData_to_all(img, myid, main_node, comm)
	# The latter is inconsistent with mpi_bcast() and difficult to implement efficiently
	#
	# (2) To be consistent with send_EMData() and recv_EMData(), we assume that the node
	# other than the broadcasting node know nothing about the EMData(). Therefore, one
	# need to broadcast the size of EMData() and two attributes: is_complex and is_ri.
	# For all other attributes, you are on your own.

	pass#IMPORTIMPORTIMPORT from numpy import reshape
	pass#IMPORTIMPORTIMPORT from mpi import mpi_bcast, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1: comm = MPI_COMM_WORLD

	img_head = []
	if myid == main_node:
		img_head.append(img.get_xsize())
		img_head.append(img.get_ysize())
		img_head.append(img.get_zsize())
		img_head.append(img.is_complex())
		img_head.append(img.is_ri())
	img_head = mpi_bcast(img_head, 5, MPI_INT, main_node, comm)
	nx = int(img_head[0])
	ny = int(img_head[1])
	nz = int(img_head[2])
	is_complex = int(img_head[3])
	is_ri = int(img_head[4])

	ntot = nx*ny*nz
	img_data = EMNumPy.em2numpy(img)
	img_data = mpi_bcast(img_data, ntot, MPI_FLOAT, main_node, comm)
	if nz != 1:
		img_data = reshape(img_data, (nz, ny, nx))     # For some reason, the order should be like this -- Zhengfan Yang
	elif ny != 1:
		img_data = reshape(img_data, (ny, nx))
	else:
		pass
	img = EMNumPy.numpy2em(img_data)
	img.set_complex(is_complex)
	img.set_ri(is_ri)

	return img


def reduce_EMData_to_root(img, myid, main_node = 0, comm = -1):

	# Comment by Zhengfan Yang on 01/05/10
	#
	# Notice:
	# (1) one should use this new version of reducing EMData in the following way:
	# 	img = reduce_EMData_to_root(img, myid, main_node, comm)
	# instead of
	# 	reduce_EMData_to_root(img, myid, main_node, comm)
	# The latter is inconsistent with mpi_bcast() and difficult to implement efficiently

	pass#IMPORTIMPORTIMPORT from numpy import reshape
	pass#IMPORTIMPORTIMPORT from mpi   import mpi_reduce, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD

	if comm == -1: comm = MPI_COMM_WORLD

	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	is_complex = img.is_complex()
	is_ri = img.is_ri()
	ntot = nx*ny*nz

	img_data = EMNumPy.em2numpy(img)
	img_data = mpi_reduce(img_data, ntot, MPI_FLOAT, MPI_SUM, main_node, comm)

	if myid == main_node:
		if nz!=1:
			img_data = reshape(img_data, (nz, ny, nx))
		elif ny!=1:
			img_data = reshape(img_data, (ny, nx))
		else:
			pass

		img = EMNumPy.numpy2em(img_data)
		img.set_complex(is_complex)
		img.set_ri(is_ri)
		return img
	else:
		return img
'''
'''6
	count = 100000
	data1d = reshape(img_data, (ntot,))
	ntime = (ntot-1) /count + 1

	for i in xrange(ntime):
		block_begin = i*count
		block_end   = i*count + count
		if block_end > ntot:
		    block_end = ntot
		block_size  = block_end - block_begin
		mpi_send(data1d[block_begin], block_size, MPI_FLOAT, dst, data_tag*ntime+i, comm)
	'''
'''7
	#construct a EMData by taking the ownership of numpy array, no memory copying  --Grant Tang
	#recv_data_numeric = mpi_recv(ntot, MPI_FLOAT, src, data_tag, comm)
	#recv_data_numpy = numpy.array(recv_data_numeric)
	#numpy_data = recv_data.reshape(recv_data, (nz,ny,nx))
	#img = EMNumPy.numpy2em(numpy_data)

	#comment out Wei's original code, which makes memory copy to construct EMData from numpy array  --Grant Tang
	img = EMData()
	img.set_size(nx, ny, nz)
	if( complex > 0 ):
		img.set_complex(True)
	else:
		img.set_complex(False)

	data1d = reshape( get_image_data(img), (ntot,) )
	tmp_data = mpi_recv(ntot, MPI_FLOAT, src, data_tag, comm)
        data1d[0:ntot] = tmp_data[0:ntot]


	count = 100000
	ntime = (ntot-1)/count + 1

	for i in xrange(ntime):
		block_begin = i*count
		block_end   = i*count + count
		if block_end > ntot:
		    block_end = ntot
		block_size  = block_end - block_begin
		tmp_data = mpi_recv(block_size, MPI_FLOAT, src, data_tag*ntime+i, comm)
		data1d[block_begin:block_end] = tmp_data[0:block_size]

	return img
	'''
"""8
#  This would not work on windows
def memory_usage():
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from string import split
	return 0
	file = "/proc/%d/status" % os.getpid()
	f = open(file, 'r')
	line = f.readline()
	while len(line) > 0 :
		items = split( line )
		if items[0]=='VmSize:':
			return items[1]+items[2]
		line = f.readline()
"""
"""9
getang3 = angle_between_projections_directions
def getang3(p1,p2):
	pass#IMPORTIMPORTIMPORT from utilities import getfvec, lacos
	n1 = getfvec(p1[0],p1[1])
	n2 = getfvec(p2[0],p2[1])
	return lacos(n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2])
"""
"""10
	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i
	"""
"""11
Util.nearest_ang(vecs, vec[0],vec[1],vec[2])
def closest_ang( vecs, vec) :
	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i
"""
"""12
def assign_projangles(projangles, refangles, return_asg = False):

	if len(refangles) > 10000:
		if len(refangles) > 100000:
			coarse_refangles = even_angles(1.5)   # 9453 angles
		else:
			coarse_refangles = even_angles(5.0)   # 849 angles
		coarse_asg = assign_projangles(projangles, coarse_refangles, True)
		ref_asg = assign_projangles(refangles, coarse_refangles, True)
	else:
		coarse_refangles = []
		coarse_asg = []
		ref_asg = []

	nproj = len(projangles)
	nref = len(refangles)
	proj_ang = [0.0]*(nproj*2)
	ref_ang = [0.0]*(nref*2)
	for i in xrange(nproj):
		proj_ang[i*2] = projangles[i][0]
		proj_ang[i*2+1] = projangles[i][1]
	for i in xrange(nref):
		ref_ang[i*2] = refangles[i][0]
		ref_ang[i*2+1] = refangles[i][1]

	asg = Util.assign_projangles(proj_ang, ref_ang, coarse_asg, ref_asg, len(coarse_refangles))
	if return_asg: return asg
	assignments = [[] for i in xrange(nref)]
	for i in xrange(nproj):
		assignments[asg[i]].append(i)

	return assignments
"""
"""13
Pushed to C
Util.cone_dirs_f(projdirs, ancordir, ant)
Returns a list of projdirs indexes that are within ant degrees of ancordir
ant in degrees
def cone_dirs_f( projdirs, ancordir, ant):
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians
	#  ancordir contains a list of symmetry neighbors
	#  Returns a list of projdirs indexes that are within ant degrees of ancordir
	cone = cos(radians(ant))
	la = []
	for i,vecs in enumerate(projdirs):
		s = -2.0
		for d in ancordir:
			s = max(vecs[0]*d[0] + vecs[1]*d[1] + vecs[2]*d[2], s)
		if s >= cone :
			la.append(i)
	return la
"""
"""14
def cone_ang_f_with_index( projangles, phi, tht, ant ):
	pass#IMPORTIMPORTIMPORT from utilities import getvec
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians
	# vec = getvec( phi, tht )
	vec = getfvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	index = []
	for i in xrange( len(projangles) ):
		vecs = getfvec( projangles[i][0], projangles[i][1] )
		s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
		if s >= cone:
		# if abs(s) >= cone:
			la.append(projangles[i])
			index.append(i)
	return la, index
"""
'''15
def cone_vectors( normvectors, phi, tht, ant ):
	pass#IMPORTIMPORTIMPORT from utilities import getvec
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians
	vec = getvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	for i in xrange( len(normvectors) ):
		s = abs(normvectors[i][0]*vec[0] + normvectors[i][1]*vec[1] + normvectors[i][2]*vec[2])
		if s >= cone:
			la.append(normvectors[i])

	return la
'''
"""16
def symmetry_related(angles, symmetry):  # replace by s.symmetry_related
	if( (symmetry[0] == "c") or (symmetry[0] == "d") ):
		temp = Util.symmetry_related(angles, symmetry)
		nt = len(temp)/3
		return [[temp[l*3+i] for i in xrange(3)] for l in xrange(nt) ]
	else:
		pass#IMPORTIMPORTIMPORT from EMAN2 import Transform
		neighbors = []
		junk = Transform({"type":"spider","phi":angles[0],"theta":angles[1],"psi":angles[2]})
		junk = junk.get_sym_proj(symmetry)
		for p in junk:
			d = p.get_params("spider")
			neighbors.append([d["phi"],d["theta"],d["psi"]])
		return neighbors

def symmetry_related_normals(angles, symmetry):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Transform
	neighbors = []
	junk = Transform({"type":"spider","phi":angles[0],"theta":angles[1],"psi":angles[2]})
	junk = junk.get_sym_proj(symmetry)
	for p in junk:
		neighbors.append(p.get_matrix()[8:11])
	return neighbors
"""
"""17
			print(q)
			print(m)
			ERROR("balance angles","Fill up upper",1)
			"""
"""18
Not used and possibly incorrect
def phi_theta_to_xyz(ang):
	pass#IMPORTIMPORTIMPORT from math import sin, cos, pi, radians
	phi   = radians( ang[0] )
	theta = radians( ang[1] )
	z = cos(theta)
	x = sin(theta) * cos(phi)
	y = sin(theta) * sin(phi)
	return [x, y, z]


def xyz_to_phi_theta(xyz):
	pass#IMPORTIMPORTIMPORT from math import pi, acos, sqrt, degrees, atan2
	theta = acos(xyz[2])
	phi   = atan2(xyz[1], xyz[0])
	return [ degrees(phi), degrees(theta), 0.0]

# input: list of triplets (phi, theta, psi)
# output: average triplet: (phi, theta, psi)
def average_angles(angles):
	pass#IMPORTIMPORTIMPORT from math import sqrt
	# convert to x, y, z
	ex = 0.0
	ey = 0.0
	ez = 0.0
	sum_psi = 0.0
	for ang in angles:
		xyz = phi_theta_to_xyz(ang)
		ex += xyz[0]
		ey += xyz[1]
		ez += xyz[2]
		sum_psi += ang[2]
	ex /= len(ang)
	ey /= len(ang)
	ez /= len(ang)
	divider = sqrt(ex*ex + ey*ey + ez*ez)
	xyz = [ ex/divider, ey/divider, ez/divider ]
	r = xyz_to_phi_theta(xyz)
	r[2] = sum_psi / len(angles)  # TODO - correct psi calculations !!!!
	return r
"""
'''19
		for j in xrange(N):
			d = ang_diff(vec[j], vec[k])
			g[j][0] = d[0]
			g[j][1] = d[1]
			g[j][2] = j
		g.sort()
		for j in xrange(img_per_grp):
			neighbor2[j] = g[j][2]
			dis2[j] = g[j][0]
		t3 = time()
		print "Members in common = %3d   extra delta = %6.3f   time1 = %5.2f   time2 = %5.2f"%(len(Set(neighbor).intersection(Set(neighbor2))),
		dis[-1]-dis2[-1], t2-t1, t3-t2)
		'''
"""20
def findall(val, lo):
	'''
	  Find all occurrences of val in list lo
	  Returns a list of indices of val in lo.
	'''
	u = []
	i = -1
	while( i < len(lo)-1):
		try:
			i = lo.index(val,i+1)
			u.append(i)
		except:
			i += 1
	return  u
"""
"""21
Iterators over sequence of images. They work for lists and stacks of images.
Additional time cost: about 1 second per 10^6 iterations on 3 GHz processor.

Usage:

it = iterImagesList(list_of_images)  <-or->  it = iterImagesStack(stack_with_images)
while it.goToNext():
	do_something(it.image())

"""
'''22
	if( myid == main_node ):
		print "  "
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print  line, "Reading data  onx: %3d, nx: %3d, CTF: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], Tracker["applyctf"], preshift)
		print  "                       stack:      %s\n                       partids:     %s\n                       partstack: %s\n"%(Tracker["constants"]["stack"], partids, partstack)
	'''
'''23
def get_shrink_data(Tracker, nxinit, partids, partstack, bckgdata = None, myid = 0, main_node = 0, nproc = 1, \
					original_data = None, return_real = False, preshift = False, apply_mask = True, large_memory = True):
	"""
	This function will read from stack a subset of images specified in partids
	   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	So, the lengths of partids and partstack are the same.
	  The read data is properly distributed among MPI threads.

	Flow of data:
	1. Read images, if there is enough memory, keep them as original_data.
	2. Read current params
	3.  Apply shift
	4.  Normalize outside of the radius
	5.  Do noise substitution and cosine mask.  (Optional?)
	6.  Shrink data.
	7.  Apply CTF.

	"""
	#from fundamentals import resample
	pass#IMPORTIMPORTIMPORT from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt


	if( myid == main_node ):
		print "  "
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if(original_data == None or not large_memory):
			print  line, "Reading data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, Tracker["applyctf"], preshift)
		else:
			print  line, "Processing data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, Tracker["applyctf"], preshift)
		print  "                       stack:      %s\n                       partids:     %s\n                       partstack: %s\n"%(Tracker["constants"]["stack"], partids, partstack)
	if( myid == main_node ): lpartids = read_text_file(partids)
	else:  lpartids = 0
	lpartids = wrap_mpi_bcast(lpartids, main_node)
	ndata = len(lpartids)
	if( myid == main_node ):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, main_node)
	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1
	else:
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	#  Preprocess the data
	mask2D  = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima = image_end - image_start
	oldshifts = [[0.0,0.0]]*nima
	data = [None]*nima
	if(original_data == None or not large_memory): original_data = [None]*nima
	shrinkage = nxinit/float(Tracker["constants"]["nnxo"])


	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	radius = int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm = float(nxinit-(nxinit//2+1) - radius)
	txl = float(radius - nxinit//2+1)

	if bckgdata :
		nnx = bckgdata[0].get_xsize()
		nny = bckgdata[0].get_ysize()
		bckgnoise = []
		oneover = []
		for i in xrange(nny):
			prj = [0.0]*nnx
			prj[0] = 1.0
			for k in xrange(1,nnx): prj[k] = bckgdata[0].get_value_at(k,i)
			oneover.append(prj)
			for k in xrange(1,nnx):
				if( prj[k] > 0.0 ):  prj[k] = 1.0/sqrt(prj[k])
			bckgnoise.append(prj)

		datastamp = bckgdata[1]

	for im in xrange(nima):
		if(original_data[im] == None or not large_memory):
			original_data[im] = get_im(Tracker["constants"]["stack"], lpartids[im])

		phi,theta,psi,sx,sy = partstack[im][0], partstack[im][1], partstack[im][2], partstack[im][3], partstack[im][4]
		if preshift:
			data[im] = fshift(original_data[im], sx, sy)
			oldshifts[im] = [sx,sy]
			sx = 0.0
			sy = 0.0
		else:  data[im] = original_data[im].copy()
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if data[im].get_attr_default("bckgnoise", None) :  data[im].delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if bckgdata:
			try:
				stmp = data[im].get_attr("ptcl_source_image")
			except:
				try:
					stmp = data[im].get_attr("ctf")
					stmp = round(stmp.defocus,4)
				except:
					ERROR("Either ptcl_source_image or ctf has to be present in the header.","get_shrink_data",1, myid)
			try:
				indx = datastamp.index(stmp)
			except:
				ERROR("Problem with indexing ptcl_source_image.","get_shrink_data",1, myid)

			data[im].set_attr("bckgnoise",oneover[indx])
			if apply_mask:
				bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
				bckg.set_attr("is_complex",1)
				bckg.set_attr("is_fftpad",1)
				bckg = fft(filt_table(bckg,bckgnoise[indx]))
				#  Normalize bckg noise in real space, only region actually used.
				st = Util.infomask(bckg, mask2D, False)
				bckg -= st[0]
				bckg /= st[1]
				pass#IMPORTIMPORTIMPORT from morphology import cosinemask
				data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			#  if no bckgnoise, do simple masking instead
			if apply_mask:  data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"] )

		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			data[im] = filt_ctf(data[im], data[im].get_attr("ctf"))
			data[im].set_attr('ctf_applied', 1)
		else:  apix = Tracker["constants"]["pixel_size"]

		#  resample will properly adjusts shifts and pixel size in ctf
		#data[im] = resample(data[im], shrinkage)
		#  return Fourier image
		data[im] = fdecimate(data[im], nxinit, nxinit, 1, return_real)
		try:
			ctf_params = original_data[im].get_attr("ctf")
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
		except:  pass
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			pass
		else:  data[im].set_attr('apix', apix/shrinkage)
		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#  We will always set it to simplify the code
		###set_params_proj(data[im],[phi,theta,psi,0.0,0.0], "xform.anchor")
	assert( nxinit == data[0].get_ysize() )  #  Just to make sure.
	#oldshifts = wrap_mpi_gatherv(oldshifts, main_node, MPI_COMM_WORLD)
	return data, oldshifts, original_data
'''
"""24
def debug_mpi_barrier(comm):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, mpi_comm_rank, mpi_bcast
	pass#IMPORTIMPORTIMPORT from traceback import extract_stack
	pass#IMPORTIMPORTIMPORT import sys


	# if mpi_comm_rank(comm) in range(4):
	print "Stack info::0::", extract_stack()[-3:]

	sys.stdout.flush()
	sys.stderr.flush()
	return mpi_barrier(comm)


def debug_mpi_bcast(newv, s, t, m, comm):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_bcast
	pass#IMPORTIMPORTIMPORT from traceback import extract_stack
	pass#IMPORTIMPORTIMPORT import sys


	rrr = mpi_bcast(newv, s, t, m, comm)
	# if mpi_comm_rank(comm) in range(4):
	print "Stack info::0::", extract_stack()[-3:], "****************", newv, "####", rrr

	sys.stdout.flush()
	sys.stderr.flush()

	# return mpi_bcast(newv, s, t, m, comm)
	return rrr
"""
""" this commented previously25
		lowpass = 0.5
		ns = len(nfsc[1])
        #  This is resolution used to filter half-volumes
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.5 ):
                        lowpass = nfsc[0][i-1]
                        break
        """
def params_2D_3D(alpha, sx, sy, mirror):
	"""
		Convert 2D alignment parameters (alpha, sx, sy, mirror) into
		3D alignment parameters (phi, theta, psi, s2x, s2y)
	"""
	phi = 0
	psi = 0
	theta = 0
	alphan, s2x, s2y, scalen = compose_transform2(0, sx, sy, 1, -alpha, 0, 0, 1)
	if mirror > 0:
		phi   = (540.0 + phi)%360.0
		theta = 180.0  - theta
		psi   = (540.0 - psi + alphan)%360.0
	else:
		psi   = (psi   + alphan)%360.0
	return phi, theta, psi, s2x, s2y

def params_3D_2D(phi, theta, psi, s2x, s2y):
	"""
		Convert 3D alignment parameters (phi, theta, psi, s2x, s2y)  # there is no mirror in 3D!
		into 2D alignment parameters (alpha, sx, sy, mirror)
	"""
	if theta > 90.0:
		mirror = 1
		alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 540.0-psi, 0, 0, 1.0)
	else:
		mirror = 0
		alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 360.0-psi, 0, 0, 1.0)
	return  alpha, sx, sy, mirror

# Amoeba uses the simplex method of Nelder and Mead to maximize a
# function of 1 or more variables.
#
#   Copyright (C) 2005  Thomas R. Metcalf
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
def amoeba_multi_level(var, scale, func, ftolerance=1.e-4, xtolerance=1.e-4, itmax=500, data=None):
	"""
	Commented by Zhengfan Yang on 05/01/07

	I made some change to the original amoeba so that it can now pass out some values
	calculated by func other than the criteria. This is important in multi-level
	amoeba refinement because otherwise, upper level refinement will lose the information
	of lower level refinement.
	"""
	#print " ENTER AMOEBA MULTI LEVEL"
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, MPI_COMM_WORLD
	
	nvar = len(var)       # number of variables in the minimization
	nsimplex = nvar + 1   # number of vertices in the simplex

	# first set up the simplex

	simplex = [0]*(nvar+1)  # set the initial simplex
	simplex[0] = var[:]
	for i in range(nvar):
		simplex[i+1] = var[:]
		simplex[i+1][i] += scale[i]

	fvalue = []
	for i in range(nsimplex):  # set the function values for the simplex
		result, passout = func(simplex[i], data=data)
		#print  " amoeba setting ",i,simplex[i],result, passout
		fvalue.append([result, passout])

	# Ooze the simplex to the maximum

	iteration = 0

	while 1:
		# find the index of the best and worst vertices in the simplex
		ssworst = 0
		ssbest  = 0
		for i in range(nsimplex):
			if fvalue[i][0] > fvalue[ssbest][0]:
				ssbest = i
			if fvalue[i][0] < fvalue[ssworst][0]:
				ssworst = i

		# get the average of the nsimplex-1 best vertices in the simplex
		pavg = [0.0]*nvar
		for i in range(nsimplex):
			if i != ssworst:
				for j in range(nvar): pavg[j] += simplex[i][j]
		for j in range(nvar): pavg[j] = pavg[j]/nvar # nvar is nsimplex-1
		simscale = 0.0
		for i in range(nvar):
			simscale += abs(pavg[i]-simplex[ssworst][i])/scale[i]
		simscale = simscale/nvar

		# find the range of the function values
		fscale = (abs(fvalue[ssbest][0])+abs(fvalue[ssworst][0]))/2.0
		if fscale != 0.0:
			frange = abs(fvalue[ssbest][0]-fvalue[ssworst][0])/fscale
		else:
			frange = 0.0  # all the fvalues are zero in this case

		# have we converged?
		if (((ftolerance <= 0.0 or frange < ftolerance) and    # converged to maximum
		(xtolerance <= 0.0 or simscale < xtolerance)) or  # simplex contracted enough
		(itmax and iteration >= itmax)):	     # ran out of iterations
			return simplex[ssbest],fvalue[ssbest][0],iteration,fvalue[ssbest][1]

		# reflect the worst vertex
		pnew = [0.0]*nvar
		for i in range(nvar):
			pnew[i] = 2.0*pavg[i] - simplex[ssworst][i]
		fnew = func(pnew,data=data)
		if fnew[0] <= fvalue[ssworst][0]:
			# the new vertex is worse than the worst so shrink
			# the simplex.
			for i in range(nsimplex):
				if i != ssbest and i != ssworst:
					for j in range(nvar):
						simplex[i][j] = 0.5*simplex[ssbest][j] + 0.5*simplex[i][j]
					fvalue[i]  = func(simplex[i],data=data)
				#### <--------->
			for j in range(nvar):
				pnew[j] = 0.5*simplex[ssbest][j] + 0.5*simplex[ssworst][j]
			fnew = func(pnew, data=data)
		elif fnew[0] >= fvalue[ssbest][0]:
			# the new vertex is better than the best so expand
			# the simplex.
			pnew2 = [0.0]*nvar
			for i in range(nvar):
				pnew2[i] = 3.0*pavg[i] - 2.0*simplex[ssworst][i]
			fnew2 = func(pnew2,data=data)
			if fnew2[0] > fnew[0]:
				# accept the new vertex in the simplexe
				pnew = pnew2
				fnew = fnew2
		# replace the worst vertex with the new vertex
		for i in range(nvar):
			simplex[ssworst][i] = pnew[i]
		fvalue[ssworst] = fnew
		iteration += 1


		# if mpi_comm_rank(MPI_COMM_WORLD) == 7:
		# 	print "Iteration:",iteration,"  ",ssbest,"  ", simplex[ssbest], "  ",fvalue[ssbest]

"""Multiline Comment0"""


def ce_fit(inp_image, ref_image, mask_image):
	""" Fit the histogram of the input image under mask with the reference image.

	     Usage : ce_fit(inp_image,ref_image,mask_image):
	    	 A and B, number of iterations and the chi-square
	"""
	hist_res = EMAN2_cppwrap.Util.histc(ref_image, inp_image, mask_image)
	args = hist_res["args"]
	scale = hist_res["scale"]
	data = [hist_res['data'], inp_image, hist_res["ref_freq_bin"], mask_image, int(hist_res['size_img']), hist_res['hist_len']]
	res = amoeba(args, scale, hist_func, 1.e-4, 1.e-4, 500, data)
	resu = ["Final Parameter [A,B]:", res[0], "Final Chi-square :", -1*res[1], "Number of Iteration :", res[2]]
	corrected_image = inp_image*res[0][0] + res[0][1]
	result = [resu,"Corrected Image :",corrected_image]
	del data[:], args[:], scale[:]
	return result

def center_2D(image_to_be_centered, center_method = 1, searching_range = -1, Gauss_radius_inner = 2, Gauss_radius_outter = 7, self_defined_reference = None):
	"""
		Put an input image into image center (nx/2, ny/2) using method :
		1. phase_cog
		2. cross-correlate with Gaussian function
		3. cross-correlate with donut shape image
		4. cross-correlate with reference image provided by user
		5. cross-correlate with self-rotated average
		7. binarize at ave+sigma and cross-correlate with a circle
	        The function will return centered_image, and shifts
	"""
	pass#IMPORTIMPORTIMPORT from   utilities    import peak_search
	pass#IMPORTIMPORTIMPORT from   fundamentals import fshift
	pass#IMPORTIMPORTIMPORT import types
	if type(image_to_be_centered) == bytes: image_to_be_centered = get_im(image_to_be_centered)
	if    center_method == 0 :  return  image_to_be_centered,0.,0.
	elif  center_method == 1 :
		cs = image_to_be_centered.phase_cog()
		if searching_range > 0 :
			if(abs(cs[0]) > searching_range):  cs[0]=0.0
			if(abs(cs[1]) > searching_range):  cs[1]=0.0
		return fundamentals.fshift(image_to_be_centered, -cs[0], -cs[1]), cs[0], cs[1]

	elif center_method == 7:
		pass#IMPORTIMPORTIMPORT from fundamentals import ccf, cyclic_shift
		pass#IMPORTIMPORTIMPORT from morphology   import binarize
		pass#IMPORTIMPORTIMPORT from utilities    import model_blank
		pass#IMPORTIMPORTIMPORT from EMAN2        import rsconvolution
		p = EMAN2_cppwrap.Util.infomask(image_to_be_centered,None,True)
		cc = morphology.binarize(EMAN2_cppwrap.rsconvolution(morphology.binarize(image_to_be_centered,p[0]+p[1]),model_blank(5,5,1,1.0/(5.0*5.0))),0.5)
		c = fundamentals.ccf(cc, self_defined_reference)
		p = EMAN2_cppwrap.Util.infomask(c,None,True)[3]
		nx = c.get_xsize()
		ny = c .get_ysize()
		cx = nx//2
		cy = ny//2
		n = 0
		x=0
		y=0
		for i in range(nx):
			for j in range(ny):
				if c.get_value_at(i,j) == p :
					x+=(i-cx)
					y+=(j-cy)
					n+=1
		shiftx = x/n
		shifty = y/n
		if searching_range > 0 :
			if(abs(shiftx) > searching_range):  shiftx=0
			if(abs(shifty) > searching_range):  shifty=0
		return fundamentals.cyclic_shift(image_to_be_centered, -shiftx, -shifty), shiftx, shifty

	elif center_method == 5:
		pass#IMPORTIMPORTIMPORT from fundamentals import rot_avg_image,ccf
		pass#IMPORTIMPORTIMPORT from math import sqrt
		not_centered = True
		tmp_image = image_to_be_centered.copy()
		shiftx = 0
		shifty = 0
		while (not_centered):
			reference = fundamentals.rot_avg_image(tmp_image)
			ccmap = fundamentals.ccf(tmp_image, reference)
			if searching_range > 0:  ccmap = EMAN2_cppwrap.Util.window(ccmap, searching_range, searching_range, 1, 0, 0, 0)
			peak  = peak_search(ccmap)
			centered_image = fundamentals.fshift(tmp_image, -peak[0][4], -peak[0][5])
			if numpy.sqrt(peak[0][4]**2 + peak[0][5]**2) < 1. : not_centered = False
			else : tmp_image = centered_image.copy()
			shiftx += peak[0][4]
			shifty += peak[0][5]
		return centered_image, shiftx, shifty

	elif center_method == 6:
		pass#IMPORTIMPORTIMPORT from morphology import threshold_to_minval
		nx = image_to_be_centered.get_xsize()
		ny = image_to_be_centered.get_ysize()
		r = nx//2-2
		mask = model_circle(r, nx, ny)
		[mean, sigma, xmin, xmax] = EMAN2_cppwrap.Util.infomask(image_to_be_centered, mask, True)
		new_image = morphology.threshold_to_minval(image_to_be_centered, mean+sigma)
		cs = new_image.phase_cog()
		if searching_range > 0 :
			if(abs(cs[0]) > searching_range):  cs[0]=0.0
			if(abs(cs[1]) > searching_range):  cs[1]=0.0
		return fundamentals.fshift(image_to_be_centered, -cs[0], -cs[1]), cs[0], cs[1]

	else :
		nx = image_to_be_centered.get_xsize()
		ny = image_to_be_centered.get_ysize()
		pass#IMPORTIMPORTIMPORT from fundamentals import ccf
		if center_method == 2 :
			reference = model_gauss(Gauss_radius_inner, nx, ny)
		if center_method == 3 :
			do1 = model_gauss(Gauss_radius_outter, nx, ny)
			do2 = model_gauss(Gauss_radius_inner,  nx, ny)
			s = EMAN2_cppwrap.Util.infomask(do1, None, True)
			do1/= s[3]
			s = EMAN2_cppwrap.Util.infomask(do2, None, True)
			do2/=s[3]
			reference = do1 - do2
		if center_method == 4:	reference = self_defined_reference
		ccmap = fundamentals.ccf(image_to_be_centered, reference)
		if searching_range > 1: ccmap = EMAN2_cppwrap.Util.window(ccmap, searching_range, searching_range, 1, 0, 0, 0)
		peak  = peak_search(ccmap)
		return fundamentals.fshift(image_to_be_centered, -peak[0][4], -peak[0][5]), peak[0][4], peak[0][5]

def common_line_in3D(phiA,thetaA,phiB,thetaB):
	"""Find the position of the commone line in 3D
           Formula is   (RB^T zhat)   cross  (RA^T zhat)
	   Returns phi, theta of the common line in degrees. theta always < 90
	   Notice you don't need to enter psi's; they are irrelevant
	"""

	pass#IMPORTIMPORTIMPORT from math import pi, sqrt, cos, sin, asin, atan2

	piOver=numpy.pi/180.0;
	ph1 = phiA*piOver;
	th1 = thetaA*piOver;
	ph2 = phiB*piOver;
	th2 = thetaB*piOver;

	#nx = cos(thetaBR)*sin(thetaAR)*sin(phiAR) - cos(thetaAR)*sin(thetaBR)*sin(phiBR) ;
	#ny = cos(thetaAR)*sin(thetaBR)*cos(phiBR) - cos(thetaBR)*sin(thetaAR)*cos(phiAR) ;
	#nz = sin(thetaAR)*sin(thetaBR)*sin(phiAR-phiBR);


	nx = numpy.sin(th1)*numpy.cos(ph1)*numpy.sin(ph2)-numpy.sin(th2)*numpy.sin(ph1)*numpy.cos(ph2)
	ny = numpy.sin(th1)*numpy.cos(th2)*numpy.cos(ph1)*numpy.cos(ph2)-numpy.cos(th1)*numpy.sin(th2)*numpy.cos(ph1)*numpy.cos(ph2)
	nz = numpy.cos(th2)*numpy.sin(ph1)*numpy.cos(ph2)-numpy.cos(th1)*numpy.cos(ph1)*numpy.sin(ph2)

	norm = nx*nx + ny*ny + nz*nz

	if norm < 1e-5:
		#print 'phiA,thetaA,phiB,thetaB:', phiA, thetaA, phiB, thetaB
		return 0.0, 0.0

	if nz<0: nx=-nx; ny=-ny; nz=-nz;

	#thetaCom = asin(nz/sqrt(norm))
	phiCom    = math.asin(nz/numpy.sqrt(norm))
	#phiCom   = atan2(ny,nx)
	thetaCom  = math.atan2(ny, nx)

	return phiCom*180.0/numpy.pi , thetaCom*180.0/numpy.pi

def compose_transform2m(alpha1=0.0, sx1=0., sy1=0.0, mirror1=0, scale1=1.0, alpha2=0.0, sx2=0.0, sy2=0.0, mirror2=0, scale2=1.0):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
		     with T1 described by alpha1, sx1, scale1 etc.

	  Combined parameters correspond to image first transformed by set 1 followed by set 2.

	    Usage: compose_transform2(alpha1,sx1,sy1,mirror1,scale1,alpha2,sx2,sy2,mirror2,scale2)
	       angles in degrees
	"""

	t1 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":mirror1,"scale":scale1})
	t2 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":mirror2,"scale":scale2})
	tt = t2*t1
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], int(d[ "mirror" ]+0.1), d[ "scale" ]

def inverse_transform3(phi, theta=0.0, psi=0.0, tx=0.0, ty=0.0, tz=0.0, mirror = 0, scale=1.0):
	"""Returns the inverse of the 3d rot and trans matrix

	    Usage: nphi,ntheta,npsi,ntx,nty,ntz,nmirror,nscale = inverse_transform3(phi,theta,psi,tx,ty,tz,mirror,scale)
	       angles in degrees
	"""

	d = EMAN2_cppwrap.Transform({'type': 'spider', 'phi': phi, 'theta': theta, 'psi': psi, 'tx': tx, 'ty': ty, 'tz': tz, "mirror":mirror,"scale":scale})
	d = d.inverse()
	d = d.get_params("spider")
	return  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],int(d["mirror"]+0.1),d["scale"]

def create_spider_doc(fname,spiderdoc):
	"""Convert a text file that is composed of columns of numbers into spider doc file
	"""
	pass#IMPORTIMPORTIMPORT from string import atoi,atof
	infile = open(fname,"r")
	lines  = infile.readlines()
	infile.close()
	nmc  = len(lines[0].split())
	table=[]
	for line in lines:
		data = line.split()
	for i in range(0,nmc):
		data[i] = string.atof(data[i])
		table.append(data)
	drop_spider_doc(spiderdoc ,table)

def drop_png_image(im, trg):
	"""Write an image with the proper png save
	Usage: drop_png_image(name_of_existing_image, 'path/to/image.png')
	"""

	if trg[-4:] != '.png':
		global_def.ERROR('destination name must be png extension', 'drop_png_image', 1)

	if isinstance(trg, str):
		im['render_min'] = im['minimum']
		im['render_max'] = im['maximum']
		im.write_image(trg, 0)
	else:
		global_def.ERROR('destination is not a file name', 'drop_png_image', 1)

def dump_row(input, fname, ix=0, iz=0):
	"""Output the data in slice iz, row ix of an image to standard out.

	Usage: dump_row(image, ix, iz)
	   or
	       dump_row("path/to/image", ix, iz)
	"""
	fout = open(fname, "w")
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	fout.write("# z = %d slice, x = %d row)\n" % (iz, ix))
	line = []
	for iy in range(ny):
		fout.write("%d\t%12.5g\n" % (iy, image.get_value_at(ix,iy,iz)))
	fout.close()

def eigen_images_get(stack, eigenstack, mask, num, avg):
	"""
		Perform PCA on stack file
		and Get eigen images
	"""

	pass#IMPORTIMPORTIMPORT from utilities import get_image

	a = EMAN2_cppwrap.Analyzers.get('pca_large')
	e = EMAN2_cppwrap.EMData()
	if(avg == 1): s = EMAN2_cppwrap.EMData()
	nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
	for im in range(nima):
		e.read_image(stack,im)
		e *= mask
		a.insert_image(e)
		if( avg==1):
			if(im==0): s  = a
			else:      s += a
	if(avg == 1): a -= s/nima
	eigenimg = a.analyze()
	if(num>= EMAN2_cppwrap.EMUtil.get_image_count(eigenimg)):
		num=EMAN2_cppwrap.EMUtil.get_image_count(eigenimg)
	for  i in range(num): eigenimg.write_image(eigenstack,i)

def find_inplane_to_match(phiA,thetaA,phiB,thetaB,psiA=0,psiB=0):
	"""Find the z rotation such that
	    ZA  RA is as close as possible to RB
	        this maximizes trace of ( RB^T ZA RA) = trace(ZA RA RB^T)
	"""
	#from math import pi, sqrt, cos, acos, sin

	RA   = EMAN2_cppwrap.Transform({'type': 'spider', 'phi': phiA, 'theta': thetaA, 'psi': psiA})
	RB   = EMAN2_cppwrap.Transform({'type': 'spider', 'phi': phiB, 'theta': thetaB, 'psi': psiB})
	RBT  = RB.transpose()
	RABT = RA * RBT

	RABTeuler = RABT.get_rotation('spider')
	RABTphi   = RABTeuler['phi']
	RABTtheta = RABTeuler['theta']
	RABTpsi   = RABTeuler['psi']

	#deg_to_rad = pi/180.0
	#thetaAR = thetaA*deg_to_rad
	#thetaBR = thetaB*deg_to_rad
	#phiAR   = phiA*deg_to_rad
	#phiBR   = phiB *deg_to_rad

	#d12=cos(thetaAR)*cos(thetaBR) + sin(thetaAR)*sin(thetaBR)*cos(phiAR-phiBR)
	return (-RABTpsi-RABTphi),RABTtheta #  180.0*acos(d12)/pi;

def get_sym(symmetry):
	"""
	get a list of point-group symmetry angles, symmetry="c3"
	"""
	pass#IMPORTIMPORTIMPORT from fundamentals import symclass
	scl = fundamentals.symclass(symmetry)
	return scl.symangles

def get_textimage(fname):
	"""
		Return an image created from a text file.  The first line of
		the image should contain "nx ny nz" (separated by whitespace)
		All subsequent lines contain "ix iy iz val", where ix, iy,
		and iz are the integer x, y, and z coordinates of the point
		and val is the floating point value of that point.  All points
		not explicitly listed are set to zero.
	"""
	pass#IMPORTIMPORTIMPORT from string import atoi,atof
	infile = open(fname)
	lines = infile.readlines()
	infile.close()
	data = lines[0].split()
	nx = string.atoi(data[0])
	ny = string.atoi(data[1])
	nz = string.atoi(data[2])
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	e.to_zero()
	for line in lines[1:]:
		data = line.split()
		ix = string.atoi(data[0])
		iy = string.atoi(data[1])
		iz = string.atoi(data[2])
		val = string.atof(data[3])
		e[ix,iy,iz] = val
	return e

def hist_func(args, data):
	#Util.hist_comp_freq(float PA,float PB,int size_img, int hist_len, float *img_ptr, float *ref_freq_bin, float *mask_ptr, float ref_h_diff, float ref_h_min)
	return EMAN2_cppwrap.Util.hist_comp_freq(args[0],args[1],data[4],data[5],data[1],data[2],data[3],data[0][0],data[0][1])

def info(image, mask=None, Comment=""):
	"""Calculate and print the descriptive statistics of an image.

	Usage: [mean, sigma, xmin, xmax, nx, ny, nz =] info(image object)
	       or
	       [mean, sigma, xmin, xmax, nx, ny, nz =] info("path/image")

	Purpose: calculate basic statistical characteristics of an image.
	"""
	if(Comment):  print(" ***  ", Comment)
	e = get_image(image)
	[mean, sigma, imin, imax] = EMAN2_cppwrap.Util.infomask(e, mask, True)
	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()
	if (e.is_complex()):
		s = ""
		if e.is_shuffled():
			s = " (shuffled)"
		if (e.is_fftodd()):
			print("Complex odd image%s: nx = %i, ny = %i, nz = %i" % (s, nx, ny, nz))
		else:
			print("Complex even image%s: nx = %i, ny = %i, nz = %i" % (s, nx, ny, nz))

	else:
		print("Real image: nx = %i, ny = %i, nz = %i" % (nx, ny, nz))

	print("avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax))
	return mean, sigma, imin, imax, nx, ny, nz

#### -----M--------
def model_square(d, nx, ny, nz=1):
	"""
	Create a centered square (or cube) with edge length of d.
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.squarecube", {"edge_length":d, "fill":1})
	return e

def model_cylinder(radius, nx, ny, nz):
	"""
	 create a cylinder along z axis
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.cylinder", {"radius":radius})
	return  e

def set_seed(sde):
	pass#IMPORTIMPORTIMPORT from random import seed
	random.seed(int(sde))
	e = EMAN2_cppwrap.EMData()
	e.set_size(1,1,1)
	e.process_inplace("testimage.noise.gauss", {"sigma":1.0, "seed":int(sde)})

###----P-------
def parse_spider_fname(mystr, *fieldvals):
	"""
	Parse a Spider filename string and insert parameters.

	Example input: "foo{***}/img{****}.mrc"
	This string has two fields that should be replaced by integers,
	and the number of '*'s determines how "wide" that field should be.
	So, if the parameters to be inserted are 10 and 3, then the resulting
	filename should be "foo010/img0003.mrc".

	Note: If the image is a stack file, the last character in the string
	must be a '@' (except for possible extraneous whitespace, which is
	ignored).  This stack symbol will be stripped in the output filename.

	Example:

	   In [1]: mystr = "foo{***}/img{****}.mrc"
	   In [2]: parse_spider_fname(mystr, 10, 3)
	   Out[2]: 'foo010/img0003.mrc'

	@param mystr Spider filename string to be parsed
	@param fieldvals Integer values to be placed into the fields

	@return Parsed filename
	"""
	# helper functions and classes
	def rm_stack_char(mystr):
		"Helper function to remove a stack character if it exists"
		stackloc = mystr.find("@")
		if stackloc != -1:
			# there's an '@' somewhere
			if len(mystr) - 1 == stackloc:
				# It's at the end of the string
				return mystr[:-1]
			else:
				# '@' not at the end, so it's an error
				raise ValueError("Invalid format: misplaced '@'.")
		else:
			# no '@' at all
			return mystr
	class Fieldloc(object):
		"Helper class to store description of a field"
		def __init__(self, begin, end):
			self.begin = begin
			self.end = end
		def count(self):
			"Size of the field (including braces)"
			return self.end - self.begin + 1
	def find_fields(mystr):
		"Helper function to identify and validate fields in a string"
		fields = []
		loc = 0
		while True:
			begin = mystr.find('{', loc)
			if begin == -1: break
			end = mystr.find('}', begin)
			field = Fieldloc(begin, end)
			# check validity
			asterisks = mystr[begin+1:end]
			if asterisks.strip("*") != "":
			    raise ValueError("Malformed {*...*} field: %s" % \
				mystr[begin:end+1])
			fields.append(Fieldloc(begin, end))
			loc = end
		return fields
	# remove leading whitespace
	mystr.strip()
	# remove stack character (if it exists)
	mystr = rm_stack_char(mystr)
	# locate fields to replace
	fields = find_fields(mystr)
	if len(fields) != len(fieldvals):
		# wrong number of fields?
		raise ValueError("Number of field values provided differs from" \
			"the number of {*...*} fields.")
	newstrfrags = []
	loc = 0
	for i, field in enumerate(fields):
		# text before the field
		newstrfrags.append(mystr[loc:field.begin])
		# replace the field with the field value
		fieldsize = field.count() - 2
		fielddesc = "%0" + str(fieldsize) + "d"
		newstrfrags.append(fielddesc % fieldvals[i])
		loc = field.end + 1
	newstrfrags.append(mystr[loc:])
	return "".join(newstrfrags)

def print_row(input, ix=0, iz=0):
	"""Print the data in slice iz, row ix of an image to standard out.

	Usage: print_row(image, ix, iz)
	   or
	       print_row("path/to/image", ix, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice, x = %d row)" % (iz, ix))
	line = []
	for iy in range(ny):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((iy + 1) % 10 == 0): line.append("\n")
	line.append("\n")
	print("".join(line))

def print_col(input, iy=0, iz=0):
	"""Print the data in slice iz, column iy of an image to standard out.

	   Usage: print_col(image, iy, iz)
	      or
		  print_col("path/to/image", iy, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice, y = %d col)" % (iz, iy))
	line = []
	for ix in range(nx):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((ix + 1) % 10 == 0): line.append("\n")
	line.append("\n")
	print("".join(line))

def print_slice(input, iz=0):
	"""Print the data in slice iz of an image to standard out.

	Usage: print_image(image, int)
	   or
	       print_image("path/to/image", int)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice)" % (iz))
	line = []
	for iy in range(ny):
		line.append("Row ")
		line.append("%4i " % iy)
		for ix in range(nx):
			line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
			if ((ix + 1) % 10 == 0):
				line.append("\n")
				line.append("         ")
		line.append("\n")
		if(nx%5 != 0): line.append("\n")
	print("".join(line))

def print_image(input):
	"""Print the data in an image to standard out.

	Usage: print_image(image)
	   or
	       print_image("path/to/image")
	"""
	image=get_image(input)
	nz = image.get_zsize()
	for iz in range(nz): print_slice(input, iz)


def print_image_col(input, ix=0, iz=0):
	"""Print the data in slice iz, row ix of an image to standard out.

	Usage: print_image_col(image, ix, iz)
	   or
	       print_image_col("path/to/image", ix, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice, x = %d row)" % (iz, ix))
	line = []
	for iy in range(ny):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((iy + 1) % 10 == 0): line.append("\n")
	line.append("\n")
	print("".join(line))

def print_image_row(input, iy=0, iz=0):
	"""Print the data in slice iz, column iy of an image to standard out.

	   Usage: print_image_row(image, iy, iz)
	      or
		  print_image_row("path/to/image", iy, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice, y = %d col)" % (iz, iy))
	line = []
	for ix in range(nx):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((ix + 1) % 10 == 0): line.append("\n")
	line.append("\n")
	print("".join(line))

def print_image_slice(input, iz=0):
	"""Print the data in slice iz of an image to standard out in a format that agrees with v2

	Usage: print_image_slice(image, int)
	   or
	       print_image_slice("path/to/image", int)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print("(z = %d slice)" % (iz))
	line = []
	for iy in range(ny-1,-1,-1):
		line.append("Row ")
		line.append("%4i " % iy)
		for ix in range(nx):
			line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
			if ((ix + 1) % 5 == 0):
				line.append("\n   ")
				line.append("      ")
			line.append("\n")
			if(nx%5 != 0): line.append("\n")
	print("".join(line))

def print_image_slice_3d(input, num=0,direction="z"):
	"""Print the data in slice iz of an image to standard out in a format that agrees with v2

	Usage: print_image_slice(image, int)
	   or
	       print_image_slice("path/to/image", int)
	"""
	#print "print slice at 3 directions"
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	if(direction=="x"):
		#print "xxxxx"
		ix=num
		print("(x = %d slice)" % (ix))
		line = []
		for iz in range(nz-1,-1,-1):
			line.append("Z ")
			line.append("%4i " % iz)
			for iy in range(ny):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((iy + 1) % 5 == 0):
					line.append("\n   ")
					line.append("      ")
				line.append("\n")
				if(ny%5 != 0): line.append("\n")
		print("".join(line))
	elif(direction=="y"):
		#print "yyy"
		iy=num
		print("(y = %d slice)" % (iy))
		line = []
		for iz in range(nz-1,-1,-1):
			line.append("Z ")
			line.append("%4i " % iz)
			for ix in range(nx):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((ix + 1) % 5 == 0):
					line.append("\n   ")
					line.append("      ")
				line.append("\n")
				if(nx%5 != 0): line.append("\n")
		print("".join(line))
	else:
		#print "zzzz"
		iz=num
		print("(z = %d slice)" % (iz))
		line = []
		for iy in range(ny-1,-1,-1):
			line.append("Row ")
			line.append("%4i " % iy)
			for ix in range(nx):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((ix + 1) % 5 == 0):
					line.append("\n   ")
					line.append("      ")
				line.append("\n")
				if(nx%5 != 0): line.append("\n")
		print("".join(line))


def read_spider_doc(fnam):
	pass#IMPORTIMPORTIMPORT from string import atof, atoi
	"""
		spider doc file format:
		key nrec record ...
		5   2    12     ...(when key <=99999)
		6   2    12     ...(when key >99999)
	"""
	inf = open(fnam, "r")
	comment_line=inf.readline() # assume there is only one comment line
	docf_in = []
	data    = []
	line    = inf.readline()
	while len(line) > 0:
		line_data=[]
		if(line[11:12]==" " and line[8:10] != " "):			# new data format
			start= 13
			end  = 15
			#line_data.append(atoi(line[start:end]))		# 03/21/12 Anna: according to the previous version of this function this value was omitted
			start= end+3
			end  = start+6
			line_data.append(string.atof(line[start:end]))
			colNo = (len(line)-end)/12 - 1
			for i in range(colNo):
				start= end+6
				end  = start+7
				line_data.append(string.atof(line[start:end]))
			data.append(line_data)
			line = inf.readline()
		else:												# old data format
			if(line[5:6] == " "): ibeg = 6
			else:  ibeg = 7
			for irec in range(string.atoi(line[ibeg:ibeg+2])):
			 	start= ibeg+2+irec*12
			 	end  = ibeg+2+(irec+1)*12
			 	line_data.append(string.atof(line[start:end]))
			data.append(line_data)
			line = inf.readline()
	return data

def reconstitute_mask(image_mask_applied_file, new_mask_file, save_file_on_disk = True, saved_file_name = "image_in_reconstituted_mask.hdf"):
	pass#IMPORTIMPORTIMPORT import types
	"""
		Substitute masked area value with image average
	"""
	if type(image_mask_applied_file) == bytes:
		nima = EMAN2_cppwrap.EMUtil.get_image_count(image_mask_applied_file)
		if (nima > 1):
			image_mask_applied = []
			for ima in range(nima):
				e = EMAN2_cppwrap.EMData()
				e.read_image(image_mask_applied_file, ima)
				image_mask_applied.append(e)
		else:
			image_mask_applied = get_im(image_mask_applied_file)
	elif  type(image_mask_applied_file) == list:
		nima =  len( image_mask_applied )
		image_mask_applied = image_mask_applied_file
	if type(new_mask_file) == bytes:
		new_mask = get_im( new_mask_file )
	elif type(new_mask_file) == int or type( new_mask_file ) == types.floatType:
		if nima > 1:
			e = image_mask_applied[0]
			nx = e.get_xsize()
			ny = e.get_ysize()
			nz = e.get_zsize()
		else :
			nx = image_mask_applied.get_xsize()
			ny = image_mask_applied.get_ysize()
			nz = image_mask_applied.get_zsize()
		new_mask = model_circle(new_mask_file, nx, ny, nz)
	if  nima > 1:
		image_in_reconstituted_mask = []
		for i in range(nima):
			tmp_image = EMAN2_cppwrap.Util.reconstitute_image_mask(image_mask_applied[i], new_mask)
			image_in_reconstituted_mask.append (tmp_image)
			if (save_file_on_disk ):  image_in_reconstituted_mask[i].write_image(saved_file_name, i)
		if(not save_file_on_disk):  return  image_in_reconstituted_mask
	else :
		if(save_file_on_disk ):
			image_in_reconstituted_mask = EMAN2_cppwrap.Util.reconstitute_image_mask(image_mask_applied, new_mask)
			image_in_reconstituted_mask.write_image(saved_file_name)
		else:	return  EMAN2_cppwrap.Util.reconstitute_image_mask(image_mask_applied, new_mask)

def rotate_about_center(alpha, cx, cy):
	"""Rotate about a different center

	    Usage: rotate_about_center(alpha,cx,cy):
	       angles in degrees
	"""

	cmp1 = compose_transform2(0, -cx, -cy, 1, alpha, 0, 0, 1)
	cmp2 = compose_transform2(cmp1[0], cmp1[1], cmp1[2], cmp1[3], 0, cx, cy, 1)

	#   return compalpha, comptrans.at(0),comptrans.at(1), compscale
	return cmp2[0], cmp2[1], cmp2[2], cmp2[3]

def estimate_3D_center(data):
	pass#IMPORTIMPORTIMPORT from math import cos, sin, radians
	pass#IMPORTIMPORTIMPORT from numpy import matrix
	pass#IMPORTIMPORTIMPORT from numpy import linalg
	pass#IMPORTIMPORTIMPORT import types
	if(type(data[0]) is list):
		ali_params = data
	else:
		ali_params = []
		for im in data:
			phi, theta, psi, s2x, s2y = get_params_proj(im)
			ali_params.append([phi, theta, psi, s2x, s2y])

	N = len(ali_params)
	A = []
	b = []

	for i in range(N):
		phi_rad   = numpy.radians(ali_params[i][0])
		theta_rad = numpy.radians(ali_params[i][1])
		psi_rad   = numpy.radians(ali_params[i][2])
		A.append([numpy.cos(psi_rad)*numpy.cos(theta_rad)*numpy.cos(phi_rad)-numpy.sin(psi_rad)*numpy.sin(phi_rad),
			numpy.cos(psi_rad)*numpy.cos(theta_rad)*numpy.sin(phi_rad)+numpy.sin(psi_rad)*numpy.cos(phi_rad), -numpy.cos(psi_rad)*numpy.sin(theta_rad), 1, 0])
		A.append([-numpy.sin(psi_rad)*numpy.cos(theta_rad)*numpy.cos(phi_rad)-numpy.cos(psi_rad)*numpy.sin(phi_rad),
			-numpy.sin(psi_rad)*numpy.cos(theta_rad)*numpy.sin(phi_rad)+numpy.cos(psi_rad)*numpy.cos(phi_rad), numpy.sin(psi_rad)*numpy.sin(theta_rad), 0, 1])
		b.append([ali_params[i][3]])
		b.append([ali_params[i][4]])

	A_matrix = numpy.matrix(A)
	b_matrix = numpy.matrix(b)

	K = numpy.linalg.solve(A_matrix.T*A_matrix, A_matrix.T*b_matrix)

	return float(K[0][0]), float(K[1][0]), float(K[2][0]), float(K[3][0]), float(K[4][0])


def sym_vol(image, symmetry="c1"):
	" Symmetrize a volume"
	if(symmetry == "c1"):  return  image.copy()
	else:                  return  image.symvol(symmetry)

##----------------------------------HDF headers related code --------------------------
def start_time():
	pass#IMPORTIMPORTIMPORT import time
	start_time = time.time()
	return start_time

def finish_time(start_time):
	pass#IMPORTIMPORTIMPORT import time
	finish_time = time.time()
	print(("Running time is"), finish_time-start_time)
	return finish_time

def ttime():
	pass#IMPORTIMPORTIMPORT import time
	now = time.localtime(time.time())
	return time.asctime(now)

def running_time(start_time):
	pass#IMPORTIMPORTIMPORT from utilities import print_msg
	pass#IMPORTIMPORTIMPORT from time import time
	time_run = int(time.time() - start_time)
	time_h   = time_run / 3600
	time_m   = (time_run % 3600) / 60
	time_s   = (time_run % 3600) % 60
	print_msg('\nRunning time is: %s h %s min %s s\n\n' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0')))

def running_time_txt(start_time):
	pass#IMPORTIMPORTIMPORT from time import time
	time_run = int(time.time() - start_time)
	time_h   = time_run / 3600
	time_m   = (time_run % 3600) / 60
	time_s   = (time_run % 3600) % 60
	return 'Running time is: %s h %s min %s s' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0'))

"""Multiline Comment3"""

def bcast_compacted_EMData_all_to_all___original(list_of_em_objects, myid, comm=-1):

	"""
	The assumption in <<bcast_compacted_EMData_all_to_all>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the other ones. I

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	pass#IMPORTIMPORTIMPORT from numpy import concatenate, shape, array, split
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from numpy import reshape

	if comm == -1 or comm == None: comm = mpi.MPI_COMM_WORLD

	num_ref = len(list_of_em_objects)
	ncpu = mpi.mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = applications.MPI_start_end(num_ref, ncpu, myid)

	# used for copying the header
	reference_em_object = list_of_em_objects[ref_start]
	nx = reference_em_object.get_xsize()
	ny = reference_em_object.get_ysize()
	nz = reference_em_object.get_zsize()
	# is_complex = reference_em_object.is_complex()
	is_ri = reference_em_object.is_ri()
	changecount = reference_em_object.get_attr("changecount")
	is_complex_x = reference_em_object.is_complex_x()
	is_complex_ri = reference_em_object.get_attr("is_complex_ri")
	apix_x = reference_em_object.get_attr("apix_x")
	apix_y = reference_em_object.get_attr("apix_y")
	apix_z = reference_em_object.get_attr("apix_z")
	is_complex = reference_em_object.get_attr_default("is_complex",1)
	is_fftpad = reference_em_object.get_attr_default("is_fftpad",1)
	is_fftodd = reference_em_object.get_attr_default("is_fftodd", nz%2)

	data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size

	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print("Sending refrings: size of data to broadcast is greater than 2GB")

	for sender_id in range(ncpu):
		if sender_id == myid:
			data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[ref_start])  #array([], dtype = 'float32')
			for i in range(ref_start+1,ref_end):
				data = numpy.concatenate([data, EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[i])])
		else:
			data = numpy.array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = applications.MPI_start_end(num_ref, ncpu, sender_id)
		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		# size_of_refrings = mpi_bcast(size_of_refrings, 1, MPI_INT, sender_id, comm)
		data = mpi.mpi_bcast(data, sender_size_of_refrings, mpi.MPI_FLOAT, sender_id, comm)
		# print "Just sent %d float32 elements"%data.size

		if myid != sender_id:
			for i in range(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = numpy.reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = numpy.reshape(image_data, (ny, nx))

				em_object = EMAN2_cppwrap.EMNumPy.numpy2em(image_data)

				# em_object.set_complex(is_complex)
				em_object.set_ri(is_ri)
				em_object.set_attr_dict({
				"changecount":changecount,
				"is_complex_x":is_complex_x,
				"is_complex_ri":is_complex_ri,
				"apix_x":apix_x,
				"apix_y":apix_y,
				"apix_z":apix_z,
				'is_complex':is_complex,
				'is_fftodd':is_fftodd,
				'is_fftpad':is_fftpad})

				list_of_em_objects[i] = em_object



def gather_compacted_EMData_to_root_with_header_info_for_each_image(number_of_all_em_objects_distributed_across_processes, list_of_em_objects_for_myid_process, myid, comm=-1):

	"""

	The assumption in <<gather_compacted_EMData_to_root>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the root

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	pass#IMPORTIMPORTIMPORT from numpy import concatenate, shape, array, split
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from numpy import reshape

	if comm == -1 or comm == None: comm = mpi.MPI_COMM_WORLD

	ncpu = mpi.mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = applications.MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, myid)
	ref_end -= ref_start
	ref_start = 0

	# used for copying the header
	reference_em_object = list_of_em_objects_for_myid_process[ref_start]
	try:

		str_to_send = str(reference_em_object.get_attr_dict())

		# print "\nFFFFFFFFFFFFFF\n\n", reference_em_object.get_attr_dict(), "\n\n\n"
		# from mpi import mpi_finalize
		# mpi_finalize()
		# import sys
		# sys.stdout.flush()
		# sys.exit()
	except:
		raise ValueError("Could not convert em attribute dictionary to string s that gather_compacted_EMData_to_root_with_header_info_for_each_image can be used.")


	nx = reference_em_object.get_xsize()
	ny = reference_em_object.get_ysize()
	nz = reference_em_object.get_zsize()

	# # is_complex = reference_em_object.is_complex()
	# is_ri = reference_em_object.is_ri()
	# changecount = reference_em_object.get_attr("changecount")
	# is_complex_x = reference_em_object.is_complex_x()
	# is_complex_ri = reference_em_object.get_attr("is_complex_ri")
	# apix_x = reference_em_object.get_attr("apix_x")
	# apix_y = reference_em_object.get_attr("apix_y")
	# apix_z = reference_em_object.get_attr("apix_z")
	# is_complex = reference_em_object.get_attr_default("is_complex",1)
	# is_fftpad = reference_em_object.get_attr_default("is_fftpad",1)
	# is_fftodd = reference_em_object.get_attr_default("is_fftodd", nz%2)

	data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size

	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print("Sending refrings: size of data to broadcast is greater than 2GB")

	for sender_id in range(1,ncpu):
		if sender_id == myid:
			data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])  #array([], dtype = 'float32')
			str_to_send = str(list_of_em_objects_for_myid_process[ref_start].get_attr_dict())
			em_dict_to_send_list = [list_of_em_objects_for_myid_process[ref_start].get_attr_dict()]
			for i in range(ref_start+1,ref_end):
				data = numpy.concatenate([data, EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[i])])
				# str_to_send += str(list_of_em_objects_for_myid_process[i].get_attr_dict())
				em_dict_to_send_list.append(list_of_em_objects_for_myid_process[i].get_attr_dict())

		else:
			data = numpy.array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = applications.MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, sender_id)

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		pass#IMPORTIMPORTIMPORT from mpi import mpi_recv, mpi_send, mpi_barrier
		if myid == 0:
			# print "root, receiving from ", sender_id, "  sender_size_of_refrings = ", sender_size_of_refrings
			str_to_receive = wrap_mpi_recv(sender_id)
			em_dict_list = eval(str_to_receive)
			# print "em_dict_list", em_dict_list
			data = mpi.mpi_recv(sender_size_of_refrings, mpi.MPI_FLOAT, sender_id, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
		elif sender_id == myid:
			wrap_mpi_send(str(em_dict_to_send_list), 0)
			# print "sender_id = ", sender_id, "sender_size_of_refrings = ", sender_size_of_refrings
			mpi.mpi_send(data, sender_size_of_refrings, mpi.MPI_FLOAT, 0, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

		# if myid != sender_id:
		if myid == 0:
			for i in range(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = numpy.reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = numpy.reshape(image_data, (ny, nx))

				em_object = EMAN2_cppwrap.EMNumPy.numpy2em(image_data)
				em_object.set_attr_dict(em_dict_list[i - sender_ref_start])

				# # em_object.set_complex(is_complex)
				# em_object.set_ri(is_ri)
				# em_object.set_attr_dict({
				# "changecount":changecount,
				# "is_complex_x":is_complex_x,
				# "is_complex_ri":is_complex_ri,
				# "apix_x":apix_x,
				# "apix_y":apix_y,
				# "apix_z":apix_z,
				# 'is_complex':is_complex,
				# 'is_fftodd':is_fftodd,
				# 'is_fftpad':is_fftpad})

				# list_of_em_objects[i] = em_object
				list_of_em_objects_for_myid_process.append(em_object)

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

def gather_EMData(data, number_of_proc, myid, main_node):
	"""
	Gather the a list of EMData on all nodes to the main node, we assume the list has the same length on each node.
											It is a dangerous assumption, it will have to be changed  07/10/2015
	"""
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, MPI_INT
	pass#IMPORTIMPORTIMPORT from mpi import mpi_send, mpi_recv

	l = len(data)
	gathered_data = []
	inc = 1   # A temp measure
	if myid == main_node:
		for i in range(0, number_of_proc*inc, inc):
			if i == main_node:
				for k in range(l):
					gathered_data.append(data[k])
			else:
				for k in range(l):
					im = recv_EMData(i, i*l+k)
					mem_len = mpi.mpi_recv(1, mpi.MPI_INT, i, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
					members = mpi.mpi_recv(int(mem_len[0]), mpi.MPI_INT, i, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
					members = list(map(int, members))
					im.set_attr('members', members)
					gathered_data.append(im)
	else:
		for k in range(l):
			send_EMData(data[k], main_node, myid*l+k)
			mem = data[k].get_attr('members')
			mpi.mpi_send(len(mem), 1, mpi.MPI_INT, main_node, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
			mpi.mpi_send(mem, len(mem), mpi.MPI_INT, main_node, global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
	return gathered_data

def send_string_to_all(str_to_send, source_node = 0):
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, MPI_INT, MPI_CHAR, mpi_bcast, mpi_comm_rank

	myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	str_to_send_len  = len(str_to_send)*int(myid == source_node)
	str_to_send_len = mpi.mpi_bcast(str_to_send_len,1,mpi.MPI_INT,source_node,mpi.MPI_COMM_WORLD)[0]
	str_to_send = mpi.mpi_bcast(str_to_send,str_to_send_len,mpi.MPI_CHAR,source_node,mpi.MPI_COMM_WORLD)
	return "".join(str_to_send)


def check_attr(ima, num, params, default_value, action="Warning"):
	pass#IMPORTIMPORTIMPORT from sys import exit
	attr_list = ima.get_attr_dict()
	if (params in attr_list) == False:
		if action=="Warning":
			print("WARNING: In image %i, cannot find attribute \'%s\' in the header, set it to the default value" %(num, params), default_value)
			ima.set_attr_dict({params:default_value})
		elif action=="Error":
			print("ERROR:   In image %i, cannot find attribute \'%s\' in the header, the program has to terminate" %(num, params))
			exit()
		return False
	else: return True

def copy_attr( pin, name, pot ):
	pot.set_attr( name, pin.get_attr(name) )
	pass

def set_ctf(ima, p):
	"""
	  set EMAN2 CTF object in the header of input image using values of CTF parameters given in the list p
	  order of parameters:
        p = [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism amplitude, astigmatism angle]
	"""
	pass#IMPORTIMPORTIMPORT from utilities import generate_ctf
	ima.set_attr( "ctf", generate_ctf( p ) )

def enable_bdb_cache():
	pass#IMPORTIMPORTIMPORT import EMAN2db
	EMAN2db.BDB_CACHE_DISABLE = False



###############

# parse user function parses the --function option. this option
#    can be either a single function name (i.e. --function=ali3d_e)
#    or a list of names, specifying the location of a user-defined
#    function; this will have the form --function=/path/module/function

def parse_user_function(opt_string):
	# check if option string is a string and return None if not. this
	# will cause the user function to be set to default value
	# "ref_ali3d" in the ali functions....
	if not(type(opt_string) is str):
		return None

	# check opt_string for format:
	if (opt_string.startswith("[") and opt_string.endswith("]")):
		# options string is [path,file,function]
		opt_list = opt_string[1:-1].split(",")
		if (2 == len(opt_list)):
			# options are [file,function]
			return [opt_list[0],opt_list[1]]
		elif (3 == len(opt_list)):
			# options are [path,file,function]. note the order!
			return [opt_list[1],opt_list[2],opt_list[0]]
		else:
			# neither. assume this is an error and return default
			return None
	else:
		# no list format used, so we assume this is a function name
		# defined (and referenced) in user_functions.
		return opt_string


###############
#  Angular functions

def getang(n):
	"""
	get angle from a 2D normal vector
	"""
	pass#IMPORTIMPORTIMPORT from math import atan2, acos, degrees
	return numpy.degrees(math.atan2(n[1],n[0]))%360.0, numpy.degrees(math.acos(n[2]))%360.0

#  The other one is better written
"""Multiline Comment8"""

def nearest_ang( vecs, phi, tht ) :
	pass#IMPORTIMPORTIMPORT from utilities import getvec
	vec = getvec( phi, tht )
	return  EMAN2_cppwrap.Util.nearest_ang(vecs, vec[0],vec[1],vec[2])
	"""Multiline Comment9"""
"""Multiline Comment10"""
# This is in python, it is very slow, we keep it just for comparison, use Util.assign_projangles instead
def assign_projangles_slow(projangles, refangles):
	refnormal = [None]*len(refangles)
	for i in range(len(refangles)):
		refnormal[i] = getvec(refangles[i][0], refangles[i][1])
	assignments = [[] for i in range(len(refangles))]
	for i in range(len(projangles)):
		best_i = nearest_ang(refnormal, projangles[i][0], projangles[i][1])
		assignments[best_i].append(i)
	return assignments


def nearestk_projangles(projangles, whichone = 0, howmany = 1, sym="c1"):
	# In both cases mirrored should be treated the same way as straight as they carry the same structural information
	pass#IMPORTIMPORTIMPORT from utilities import getfvec, getvec
	lookup = list(range(len(projangles)))
	if( sym == "c1"):
		pass#IMPORTIMPORTIMPORT from utilities import getvec
		refnormal = [None]*(len(projangles)*3)
		for i in range(len(projangles)):
			ref = getvec(projangles[i][0], projangles[i][1])
			for k in range(3):
				refnormal[3*i+k] = ref[k]
		# remove the reference projection from the list
		ref = [0.0,0.0,0.0]
		for k in range(3):
			ref[k] = refnormal[3*whichone+k]
		for k in range(3): del refnormal[3*whichone+2-k]
		del lookup[whichone]
		assignments = [-1]*howmany
		for i in range(howmany):
			k = EMAN2_cppwrap.Util.nearest_ang(refnormal, ref[0],ref[1],ref[2])
			assignments[i] = lookup[k]
			for l in range(3): del refnormal[3*k+2-l]
			del lookup[k]

	elif( sym[:1] == "d" ):
		pass#IMPORTIMPORTIMPORT from utilities import get_symt, getvec
		pass#IMPORTIMPORTIMPORT from EMAN2 import Transform
		t = get_symt(sym)
		phir = 360.0/int(sym[1:])
		for i in range(len(t)):  t[i] = t[i].inverse()
		a = EMAN2_cppwrap.Transform({"type":"spider","phi":projangles[whichone][0], "theta":projangles[whichone][1]})
		for l in range(len(t)):
			q = a*t[l]
			q = q.get_params("spider")
			if(q["phi"]<phir and q["theta"] <= 90.0): break
		refvec = getfvec(q["phi"], q["theta"])
		#print  "refvec   ",q["phi"], q["theta"]

		tempan = projangles[:]
		del tempan[whichone], lookup[whichone]
		assignments = [-1]*howmany

		for i in range(howmany):
			best = -1
			for j in range(len(tempan)):
				nearest = -1.
				a = EMAN2_cppwrap.Transform({"type":"spider","phi":tempan[j][0], "theta":tempan[j][1]})
				for l in range(len(t)):
					q = a*t[l]
					q = q.get_params("spider")
					vecs = getfvec(q["phi"], q["theta"])
					s = abs(vecs[0]*refvec[0] + vecs[1]*refvec[1] + vecs[2]*refvec[2])
					if( s > nearest ):
						nearest = s
						#ttt = (q["phi"], q["theta"])
				if( nearest > best ):
					best = nearest
					best_j = j
					#print  j,tempan[j][0], tempan[j][1],best,lookup[j],ttt
			assignments[i] = lookup[best_j]
			del tempan[best_j], lookup[best_j]

	elif( sym[:1] == "c" ):
		pass#IMPORTIMPORTIMPORT from utilities import get_symt, getvec
		pass#IMPORTIMPORTIMPORT from EMAN2 import Transform
		t = get_symt(sym)
		#phir = 360.0/int(sym[1:])

		tempan =  projangles[:]
		del tempan[whichone], lookup[whichone]
		assignments = [-1]*howmany
		refvec = getvec(projangles[whichone][0], projangles[whichone][1])

		for i in range(howmany):
			best = -1
			for j in range(len(tempan)):
				nearest = -1.
				a = EMAN2_cppwrap.Transform({"type":"spider","phi":tempan[j][0], "theta":tempan[j][1]})
				for l in range(len(t)):
					q = a*t[l]
					q = q.get_params("spider")
					vecs = getfvec(q["phi"], q["theta"])
					s = abs(vecs[0]*refvec[0] + vecs[1]*refvec[1] + vecs[2]*refvec[2])
					if( s > nearest ):
						nearest = s
						#ttt = (q["phi"], q["theta"])
				if( nearest > best ):
					best = nearest
					best_j = j
					#print  j,tempan[j][0], tempan[j][1],best,lookup[j],ttt
			assignments[i] = lookup[best_j]
			del tempan[best_j], lookup[best_j]

	else:
		print("  ERROR:  symmetry not supported  ",sym)
		assignments = []


	return assignments


def nearest_full_k_projangles(reference_ang, angles, howmany = 1, sym_class=None):
	# We assume angles can be on the list of normals
	pass#IMPORTIMPORTIMPORT from utilities import getfvec, angles_to_normals
	reference_normals = angles_to_normals(reference_ang)

	if( sym_class == None or sym_class.sym[:2] == "c1"):
		ref = getfvec(angles[0],angles[1])
		assignments = EMAN2_cppwrap.Util.nearest_fang_select(reference_normals, ref[0],ref[1],ref[2], howmany)
	else:
		ancordir = angles_to_normals(sym_class.symmetry_neighbors([angles[:3]]))
		assignments = EMAN2_cppwrap.Util.nearest_fang_sym(ancordir, reference_normals, len(ancordir), howmany)

	return assignments

def nearestk_to_refdir(refnormal, refdir, howmany = 1):
	lookup = list(range(len(refnormal)))
	assignments = [-1]*howmany
	for i in range(howmany):
		k = EMAN2_cppwrap.Util.nearest_ang(refnormal, refdir[0],refdir[1],refdir[2])
		assignments[i] = lookup[k]
		del refnormal[3*k+2], refnormal[3*k+1], refnormal[3*k+0], lookup[k]
	return assignments


def nearestk_to_refdirs(refnormal, refdir, howmany = 1):
	lookup = list(range(len(refnormal)))
	assignments = []
	for j in range(len(refdir)):
		assignment = [-1]*howmany
		for i in range(howmany):
			k = EMAN2_cppwrap.Util.nearest_ang(refnormal, refdir[j][0],refdir[j][1],refdir[j][2])
			assignment[i] = lookup[k]
			del refnormal[3*k+2], refnormal[3*k+1], refnormal[3*k+0], lookup[k]
		assignments.append(assignment)
	return assignments



"""Multiline Comment11"""

def assign_projangles(projangles, refangles, return_asg = False):

	nproj = len(projangles)
	nref = len(refangles)
	proj_ang = [0.0]*(nproj*2)
	ref_ang = [0.0]*(nref*2)
	for i in range(nproj):
		proj_ang[i*2] = projangles[i][0]
		proj_ang[i*2+1] = projangles[i][1]
	for i in range(nref):
		ref_ang[i*2] = refangles[i][0]
		ref_ang[i*2+1] = refangles[i][1]

	asg = EMAN2_cppwrap.Util.assign_projangles(proj_ang, ref_ang)
	if return_asg: return asg
	assignments = [[] for i in range(nref)]
	for i in range(nproj):
		assignments[asg[i]].append(i)

	return assignments

def assign_projangles_f(projangles, refangles, return_asg = False):

	asg = EMAN2_cppwrap.Util.assign_projangles_f(projangles, refangles)
	if return_asg: return asg
	assignments = [[] for i in range(len(refangles))]
	for i in range(len(projangles)):
		assignments[asg[i]].append(i)

	return assignments


def cone_ang( projangles, phi, tht, ant, symmetry = 'c1'):
	pass#IMPORTIMPORTIMPORT from utilities import getvec, getfvec
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians

	cone = numpy.cos(numpy.radians(ant))
	la = []
	if( symmetry == 'c1' ):
		vec = getfvec( phi, tht )
		for i in range( len(projangles) ):
			vecs = getvec( projangles[i][0], projangles[i][1] )
			s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
			if s >= cone:
				la.append(projangles[i])
	elif( symmetry[:1] == "c" ):
		nsym = int(symmetry[1:])
		qt = 360.0/nsym
		dvec = 	[0.0]*nsym
		for nsm in range(nsym):
			dvec[nsm] = getvec(phi+nsm*qt, tht)
		for i in range( len(projangles) ):
			vecs = getfvec( projangles[i][0], projangles[i][1] )
			qt = -2.0
			for nsm in range(nsym):
				vc = dvec[nsm][0]*vecs[0] + dvec[nsm][1]*vecs[1] + dvec[nsm][2]*vecs[2]
				if(vc > qt):  qt = vc
			if(qt >= cone):
				la.append(projangles[i])
	elif( symmetry[:1] == "d" ):
		nsym = int(symmetry[1:])
		qt = 360.0/nsym
		dvec = 	[0.0]*2*nsym
		for nsm in range(nsym):
			dvec[2*nsm] = getvec(phi+nsm*qt, tht)
			dvec[2*nsm+1] = getvec(-(phi+nsm*qt), 180.0-tht)
		for i in range( len(projangles) ):
			vecs = getfvec( projangles[i][0], projangles[i][1] )
			qt = -2.0
			qk = -1
			for nsm in range(2*nsym):
				vc = dvec[nsm][0]*vecs[0] + dvec[nsm][1]*vecs[1] + dvec[nsm][2]*vecs[2]
				if(vc > qt):
					qt = vc
					qk = nsm
			if(qt >= cone):
				if(qk<nsym):  la.append(projangles[i])
				else:         la.append([projangles[i][0],projangles[i][1],(projangles[i][2]+180.0)%360.0])
	
	else:  print("Symmetry not supported ",symmetry)
	return la

#  Push to C.  PAP  11/25/2016
def cone_ang_f( projangles, phi, tht, ant, symmetry = 'c1'):
	pass#IMPORTIMPORTIMPORT from utilities import getfvec
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians

	cone = numpy.cos(numpy.radians(ant))
	la = []
	if( symmetry == 'c1' ):
		vec = getfvec( phi, tht )
		for i in range( len(projangles) ):
			vecs = getfvec( projangles[i][0], projangles[i][1] )
			s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
			if s >= cone:
				la.append(projangles[i])
	elif( symmetry[:1] == "c" ):
		nsym = int(symmetry[1:])
		qt = 360.0/nsym
		dvec = 	[0.0]*nsym
		for nsm in range(nsym):
			dvec[nsm] = getfvec(phi+nsm*qt, tht)
		for i in range( len(projangles) ):
			vecs = getfvec( projangles[i][0], projangles[i][1] )
			qt = -2.0
			for nsm in range(nsym):
				vc = dvec[nsm][0]*vecs[0] + dvec[nsm][1]*vecs[1] + dvec[nsm][2]*vecs[2]
				if(vc > qt):  qt = vc
			if(qt >= cone):
				la.append(projangles[i])
	elif( symmetry[:1] == "d" ):
		nsym = int(symmetry[1:])
		qt = 360.0/nsym
		dvec = 	[0.0]*2*nsym
		for nsm in range(nsym):
			dvec[2*nsm] = getfvec(phi+nsm*qt, tht)
			dvec[2*nsm+1] = getfvec(-(phi+nsm*qt), 180.0-tht)
		for i in range( len(projangles) ):
			vecs = getfvec( projangles[i][0], projangles[i][1] )
			qt = -2.0
			qk = -1
			for nsm in range(2*nsym):
				vc = dvec[nsm][0]*vecs[0] + dvec[nsm][1]*vecs[1] + dvec[nsm][2]*vecs[2]
				if(vc > qt):
					qt = vc
					qk = nsm
			if(qt >= cone):
				if(qk<nsym):  la.append(projangles[i])
				else:         la.append([projangles[i][0],projangles[i][1],(projangles[i][2]+180.0)%360.0])
	
	else:  print("Symmetry not supported ",symmetry)

	return la


"""Multiline Comment12"""

"""Multiline Comment13"""

def cone_ang_with_index( projangles, phi, tht, ant ):
	pass#IMPORTIMPORTIMPORT from utilities import getvec
	pass#IMPORTIMPORTIMPORT from math import cos, pi, degrees, radians
	# vec = getvec( phi, tht )
	vec = getfvec( phi, tht )

	cone = numpy.cos(numpy.radians(ant))
	la = []
	index = []
	for i in range( len(projangles) ):
		# vecs = getvec( projangles[i][0], projangles[i][1] )
		vecs = getfvec( projangles[i][0], projangles[i][1] )
		s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
		# if s >= cone:
		if abs(s) >= cone:
			la.append(projangles[i] + [i])
			index.append(i)

	return la, index
"""Multiline Comment14"""

#  Wrappers for new angular functions
def angles_between_anglesets(angleset1, angleset2, indexes=None):
	"""
	  It returns list of angles describing differences between the given anglesets (rotations of anglesets don't matter).
	  The function works as follow:
	  1. use the rotation_between_anglesets function to find the transformation between the anglesets
	  2. apply the transformation to the first angleset (to fit it into the second angleset)
	  3. calculate angles between corresponding projections directions from the anglesets
	  INPUT: each of the parameters should be a list of pairs [phi, theta] (additional parameters (pdi, shifts) don't influence on the results)
	  OUTPUT: list of floats - angles in degrees (the n-th element of the list equals the angle between n-th projections directions from the anglesets)
	  The third parameter (indexes) is optional and may be set to list of indexes. In that case only elements from given list are taken into account.
	"""
	pass#IMPORTIMPORTIMPORT from fundamentals import rotate_params
	pass#IMPORTIMPORTIMPORT from utilities import rotation_between_anglesets, angle_between_projections_directions

	if indexes != None:
		new_ang1 = []
		new_ang2 = []
		for i in indexes:
			new_ang1.append(angleset1[i])
			new_ang2.append(angleset2[i])
		angleset1 = new_ang1
		angleset2 = new_ang2

	rot = rotation_between_anglesets(angleset1, angleset2)
	angleset1 = fundamentals.rotate_params([angleset1[i] for i in range(len(angleset1))],[-rot[2],-rot[1],-rot[0]])
	return [angle_between_projections_directions(angleset1[i], angleset2[i]) for i in range(len(angleset1))]

"""Multiline Comment17"""

def group_proj_by_phitheta_slow(proj_ang, symmetry = "c1", img_per_grp = 100, verbose = False):
	pass#IMPORTIMPORTIMPORT from time import time
	pass#IMPORTIMPORTIMPORT from math import exp, pi

	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in range(len(ref_ang)):
			ref_ang_list[2*i] = ref_ang[i][0]
			ref_ang_list[2*i+1] = ref_ang[i][1]
		return ref_ang_list, len(ref_ang)

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		pass#IMPORTIMPORTIMPORT from utilities import lacos

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		if v >= 0: return lacos(v), 0
		else:      return lacos(-v), 1

	t0 = time.time()
	proj_list = []
	angles_list = []
	N = len(proj_ang)
	if len(proj_ang[0]) == 3:       # determine whether it has shifts provided, make the program more robust
		for i in range(N):
			proj_ang[i].append(i)
			proj_ang[i].append(True)
			vec = getfvec(proj_ang[i][0], proj_ang[i][1])     # pre-calculate the vector for each projection angles
			proj_ang[i].append(vec)
	else:
		for i in range(N):
			proj_ang[i][3] = i
			proj_ang[i][4] = True
			vec = getfvec(proj_ang[i][0], proj_ang[i][1])     # pre-calculate the vector for each projection angles
			proj_ang[i].append(vec)

	ref_ang_list1, nref1 = get_ref_ang_list(20.0, sym = symmetry)
	ref_ang_list2, nref2 = get_ref_ang_list(10.0, sym = symmetry)
	ref_ang_list3, nref3 = get_ref_ang_list(5.0, sym = symmetry)
	ref_ang_list4, nref4 = get_ref_ang_list(2.5, sym = symmetry)

	c = 100
	L = max(100, img_per_grp)
	# This is to record whether we are considering the same group as before
	# If we are, we are only going to read the table and avoid calculating the distance again.
	previous_group = -1
	previous_zone = 5
	for grp in range(N/img_per_grp):
		print(grp, end=' ')
		N_remain = N-grp*img_per_grp
		# The idea here is that if each group has more than 100 images in average,
		# we consider it crowded enough to just consider the most crowded group.
		if N_remain >= nref4*L:
			ref_ang_list = ref_ang_list4
			nref = nref4
			if previous_zone > 4:
				previous_group = -1
				previous_zone = 4
		elif N_remain >= nref3*L:
			ref_ang_list = ref_ang_list3
			nref = nref3
			if previous_zone > 3:
				previous_group = -1
				previous_zone = 3
		elif N_remain >= nref2*L:
			ref_ang_list = ref_ang_list2
			nref = nref2
			if previous_zone > 2:
				previous_group = -1
				previous_zone = 2
		elif N_remain >= nref1*L:
			ref_ang_list = ref_ang_list1
			nref = nref1
			if previous_zone > 1:
				previous_group = -1
				previous_zone = 1
		else:
			if previous_zone > 0:
				previous_group = -1
				previous_zone = 0

		t1 = time.time()
		v = []
		index = []
		if N_remain >= nref1*L:
			# In this case, assign all projection to groups and only consider the most crowded group.
			proj_ang_list = [0.0]*(N_remain*2)
			nn = 0
			remain_index = [0]*N_remain
			for i in range(N):
				if proj_ang[i][4]:
					proj_ang_list[nn*2] = proj_ang[i][0]
					proj_ang_list[nn*2+1] = proj_ang[i][1]
					remain_index[nn] = i
					nn += 1
			asg = EMAN2_cppwrap.Util.assign_projangles(proj_ang_list, ref_ang_list)
			assignments = [[] for i in range(nref)]
			for i in range(N_remain):
				assignments[asg[i]].append(i)
			# find the largest group and record the group size and group number
			max_group_size = 0
			max_group = -1
			for i in range(nref):
				if len(assignments[i]) > max_group_size:
					max_group_size = len(assignments[i])
					max_group = i
			print(max_group_size, max_group, previous_group, end=' ')
			for i in range(len(assignments[max_group])):
				ind = remain_index[assignments[max_group][i]]
				v.append(proj_ang[ind][5])
				index.append(ind)
		else:
			# In this case, use all the projections available
			for i in range(N):
				if proj_ang[i][4]:
					v.append(proj_ang[i][5])
					index.append(i)
			max_group = 0

		t2 = time.time()
		Nn = len(index)
		density = [[0.0, 0] for i in range(Nn)]
		if max_group != previous_group:
			diff_table = [[0.0 for i in range(Nn)] for j in range(Nn)]
			for i in range(Nn-1):
				for j in range(i+1, Nn):
					diff = ang_diff(v[i], v[j])
					q = numpy.exp(-c*(diff[0]/180.0*numpy.pi)**2)
					diff_table[i][j] = q
					diff_table[j][i] = q
			diff_table_index = dict()
			for i in range(Nn): diff_table_index[index[i]] = i
			print(Nn, True, end=' ')
		else:
			print(Nn, False, end=' ')

		t21 = time.time()
		for i in range(Nn):
			density[i][0] = sum(diff_table[diff_table_index[index[i]]])
			density[i][1] = i
		t22 = time.time()
		density.sort(reverse=True)

		t3 = time.time()
		dang = [[0.0, 0] for i in range(Nn)]
		most_dense_point = density[0][1]
		for i in range(Nn):
			diff = ang_diff(v[i], v[most_dense_point])
			dang[i][0] = diff[0]
			dang[i][1] = i
		dang[most_dense_point][0] = -1.
		dang.sort()

		t4 = time.time()
		members = [0]*img_per_grp
		for i in range(img_per_grp):
			idd = index[dang[i][1]]
			for j in range(len(diff_table)):
				diff_table[diff_table_index[idd]][j] = 0.0
				diff_table[j][diff_table_index[idd]] = 0.0
			members[i] = idd
			proj_ang[members[i]][4] = False
		proj_list.append(members)
		center_i = index[dang[0][1]]
		angles_list.append([proj_ang[center_i][0], proj_ang[center_i][1], dang[img_per_grp-1][0]])
		previous_group = max_group
		print(t2-t1, t3-t2, t22-t21, t3-t22, t4-t3)

	if N%img_per_grp*3 >= 2*img_per_grp:
		members = []
		for i in range(N):
			if proj_ang[i][4]:
				members.append(i)
		proj_list.append(members)
		angles_list.append([proj_ang[members[0]][0], proj_ang[members[0]][1], 90.0])
	elif N%img_per_grp != 0:
		for i in range(N):
			if proj_ang[i][4]:
				proj_list[-1].append(i)
	print("Total time used = ", time.time()-t0)

	return proj_list, angles_list

def group_proj_by_phitheta(proj_ang, symmetry = "c1", img_per_grp = 100, verbose = False):
	pass#IMPORTIMPORTIMPORT from math import exp, pi

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		pass#IMPORTIMPORTIMPORT from utilities import lacos

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		if v >= 0: return lacos(v), 0
		else:  return lacos(-v), 1


	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in range(len(ref_ang)):
			ref_ang_list[2*i] = ref_ang[i][0]
			ref_ang_list[2*i+1] = ref_ang[i][1]
		return ref_ang_list, len(ref_ang)

	N = len(proj_ang)
	proj_ang_list = [0]*(N*2)
	for i in range(N):
		proj_ang_list[i*2] = proj_ang[i][0]
		proj_ang_list[i*2+1] = proj_ang[i][1]

	ref_ang_list1, nref1 = get_ref_ang_list(20.0, sym = symmetry)
	ref_ang_list2, nref2 = get_ref_ang_list(10.0, sym = symmetry)
	ref_ang_list3, nref3 = get_ref_ang_list(5.0, sym = symmetry)
	ref_ang_list4, nref4 = get_ref_ang_list(2.5, sym = symmetry)

	ref_ang_list = []
	ref_ang_list.extend(ref_ang_list1)
	ref_ang_list.extend(ref_ang_list2)
	ref_ang_list.extend(ref_ang_list3)
	ref_ang_list.extend(ref_ang_list4)
	ref_ang_list.append(nref1)
	ref_ang_list.append(nref2)
	ref_ang_list.append(nref3)
	ref_ang_list.append(nref4)
	proj_list = EMAN2_cppwrap.Util.group_proj_by_phitheta(proj_ang_list, ref_ang_list, img_per_grp)

	proj_list2 = proj_list[:]
	for i in range(len(proj_list2)): proj_list2[i] = abs(proj_list2[i])
	proj_list2.sort()
	assert N == len(proj_list2)
	for i in range(N): assert i == proj_list2[i]

	Ng = N/img_per_grp
	proj_list_new = [[] for i in range(Ng)]
	mirror_list = [[] for i in range(Ng)]
	angles_list = []
	for i in range(Ng):
		for j in range(img_per_grp):
			proj_list_new[i].append(abs(proj_list[i*img_per_grp+j]));
			mirror_list[i].append(proj_list[i*img_per_grp+j] >= 0)
		phi1 = proj_ang[proj_list_new[i][0]][0];
		theta1 = proj_ang[proj_list_new[i][0]][1];
		phi2 = proj_ang[proj_list_new[i][-1]][0];
		theta2 = proj_ang[proj_list_new[i][-1]][1];
		angles_list.append([phi1, theta1, ang_diff(getfvec(phi1, theta1), getfvec(phi2, theta2))[0]]);

	if N%img_per_grp*3 >= 2*img_per_grp:
		proj_list_new.append([])
		mirror_list.append([])
		for i in range(Ng*img_per_grp, N):
			proj_list_new[-1].append(abs(proj_list[i]));
			mirror_list[-1].append(proj_list[i] >= 0)
		phi1 = proj_ang[proj_list_new[Ng][0]][0];
		theta1 = proj_ang[proj_list_new[Ng][0]][1];
		phi2 = proj_ang[proj_list_new[Ng][-1]][0];
		theta2 = proj_ang[proj_list_new[Ng][-1]][1];
		angles_list.append([phi1, theta1, ang_diff(getfvec(phi1, theta1), getfvec(phi2, theta2))[0]]);
	elif N%img_per_grp != 0:
		for i in range(Ng*img_per_grp, N):
			proj_list_new[-1].append(abs(proj_list[i]))
			mirror_list[-1].append(proj_list[i] >= 0)

	return proj_list_new, angles_list, mirror_list

def mulvec(v1,v2):
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def assignments_to_groups(assignments, n = -1):
	#  convert a list of assignments of images to groups to list of lists of image numbers within groups
	#  n can be maximum number of groups
	if( n == -1 ): nmax = max(assignments)
	else:  nmax = n
	groups = [[] for i in range(nmax+1)]
	for i,q in enumerate(assignments):  groups[q].append(i)
	for q in range(len(groups)):  groups[q].sort()
	return groups


def groups_assignments(groups, n = -1):
	#  convert groups (a list of lists of image numbers within groups) into a list of assignments of images to groups
	#  n can be maximum number of elements
	if( n == -1 ): 
		nmax = -1
		for q in range(len(groups)):  nmax = max(nmax,max(groups[q]))
	else:
		nmax = n
	assignments = [-1]*(nmax+1)

	for i,q in enumerate(groups):
		for l in q:
			assignments[l] = i
	return assignments

# parameters: list of integers, number of processors
def chunks_distribution(chunks, procs):
	pass#IMPORTIMPORTIMPORT from heapq import heappush, heappop

	# sort chunks in descending order
	chunks.sort(reverse=True)

	# create heap and list with solution
	results = []
	heap = []
	for p in range(procs):
		results.append([])
		heapq.heappush(heap, (0, p))

	# main calculations
	# following chunks are added to the least loaded processors
	for c in chunks:
		s, p = heapq.heappop(heap)
		results[p].append(c)
		s += c[0]
		heapq.heappush(heap, (s, p))

	return results

"""Multiline Comment20"""
# ================ Iterator for list of images
class iterImagesStack(object):
	stackName = ""
	currentImage = None
	imagesIndexes = []
	position = -1
	def __init__(self, stack_name, list_of_indexes = None):
		if list_of_indexes == None:
			self.imagesIndexes = list(range(EMAN2_cppwrap.EMUtil.get_image_count(stack_name)))
		else:
			self.imagesIndexes = list_of_indexes[:]
		self.stackName = stack_name
	def iterNo(self):
		return self.position
	def imageIndex(self):
		return self.imagesIndexes[self.position]
	def image(self):
		if self.currentImage == None:
			self.currentImage = EMAN2_cppwrap.EMData()
			self.currentImage.read_image(self.stackName, self.imagesIndexes[self.position])
		return self.currentImage
	def goToNext(self):
		self.currentImage = None
		if len(self.imagesIndexes) <= self.position:
			return False
		self.position += 1
		return (self.position < len(self.imagesIndexes))
	def goToPrev(self):
		self.currentImage = None
		if 0 > self.position:
			return False
		self.position -= 1
		return (self.position >= 0)


pass#IMPORTIMPORTIMPORT from pickle import dumps,loads
pass#IMPORTIMPORTIMPORT from zlib import compress,decompress
pass#IMPORTIMPORTIMPORT from struct import pack,unpack
pass#IMPORTIMPORTIMPORT import pickle

def rearrange_ranks_of_processors(mode):
	
	pass#IMPORTIMPORTIMPORT import socket
	pass#IMPORTIMPORTIMPORT import mpi
	pass#IMPORTIMPORTIMPORT from inspect import currentframe, getframeinfo
	
	mpi.mpi_init(0, [])
	
	original_mpi_comm_world = mpi.MPI_COMM_WORLD
	
	hostname = socket.gethostname()
	my_rank = mpi.mpi_comm_rank(original_mpi_comm_world)
	mpi_size = mpi.mpi_comm_size(original_mpi_comm_world)
	

	host_names = [hostname]
	host_names = wrap_mpi_gatherv(host_names, 0, original_mpi_comm_world)
	host_names = wrap_mpi_bcast(host_names, 0, original_mpi_comm_world)
	
	hostname_with_rank = hostname + "  %04d"%my_rank
	host_names_with_rank = [hostname_with_rank]
	host_names_with_rank = wrap_mpi_gatherv(host_names_with_rank, 0, original_mpi_comm_world)
	host_names_with_rank = wrap_mpi_bcast(host_names_with_rank, 0, original_mpi_comm_world)

	procs_belonging_to_one_node = list(map(int, sorted([ a[-4:] for a in host_names_with_rank  if hostname in a])))
	local_rank = procs_belonging_to_one_node.index(my_rank)

	local_size = host_names.count(hostname)
	
	no_of_processes_per_group = local_size
	no_of_groups = mpi_size/local_size
	
	if my_rank == 0:
		host_names = sorted(set(host_names))
	host_names = wrap_mpi_bcast(host_names, 0, original_mpi_comm_world)
	host_dict = {host_names[i]: i for i in range(len(host_names))}
	
	color = host_dict[hostname]
	
	error_status = None
	if mode == "to fit round-robin assignment":
		new_rank = local_rank * no_of_groups + color
	elif mode == "to fit by-node assignment":
		new_rank = local_rank + no_of_groups * color
	else:
		error_status = ("Invalid mode for function 'rearrange_ranks_of_processors': %s"%mode, inspect.getframeinfo(inspect.currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	
	# mpi.MPI_COMM_WORLD = mpi.mpi_comm_split(original_mpi_comm_world, color=0, key = new_rank)
	mpi.MPI_COMM_WORLD = mpi.mpi_comm_split(original_mpi_comm_world, 0, new_rank)
	
	return original_mpi_comm_world

def wrap_mpi_split_shared_memory(mpi_comm):
	pass#IMPORTIMPORTIMPORT import socket
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size, mpi_comm_split

	hostname = socket.gethostname()

	my_rank = mpi.mpi_comm_rank(mpi_comm)
	mpi_size = mpi.mpi_comm_size(mpi_comm)

	# local_rank = int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"])
	# local_size = int(os.environ["OMPI_COMM_WORLD_LOCAL_SIZE"])


	host_names = [hostname]
	host_names = wrap_mpi_gatherv(host_names, 0, mpi_comm)
	host_names = wrap_mpi_bcast(host_names, 0, mpi_comm)

	hostname_with_rank = hostname + "  %04d"%my_rank
	host_names_with_rank = [hostname_with_rank]
	host_names_with_rank = wrap_mpi_gatherv(host_names_with_rank, 0, mpi_comm)
	host_names_with_rank = wrap_mpi_bcast(host_names_with_rank, 0, mpi_comm)

	procs_belonging_to_one_node = list(map(int, sorted([ a[-4:] for a in host_names_with_rank  if hostname in a])))
	local_rank = procs_belonging_to_one_node.index(my_rank)

	local_size = host_names.count(hostname)
	# local_rank = my_rank % local_size 

	no_of_processes_per_group = local_size
	no_of_groups = mpi_size/local_size

	if my_rank == 0:
			host_names = sorted(set(host_names))
	host_names = wrap_mpi_bcast(host_names, 0, mpi_comm)
	host_dict = {host_names[i]: i for i in range(len(host_names))}

	# color = host_dict[hostname]
	color = my_rank / no_of_processes_per_group
	key = local_rank

	# shared_comm = mpi_comm_split_shared(mpi_comm, 0, key)
	shared_comm = mpi.mpi_comm_split(mpi_comm, color, key)

	return shared_comm, color, key, no_of_processes_per_group, no_of_groups	


def random_string(length_of_randomstring = 16):
	pass#IMPORTIMPORTIMPORT import random
	chars=list(map(chr, list(range(97, 123)))) # a..z
	chars.extend(list(map(chr, list(range(65, 91))))) # A..Z
	chars.extend(list(map(chr, list(range(48, 58))))) # 0..9
	random_string = ""
	for i in range(length_of_randomstring):
		random_string += chars[random.randint(0,len(chars)-1)]
	return random_string

def get_nonexistent_directory_increment_value(directory_location, directory_name, start_value = 1, myformat = "%03d"):
	pass#IMPORTIMPORTIMPORT import os
	dir_count = start_value
	while os.path.isdir(directory_location + directory_name + myformat%(dir_count)):
		dir_count += 1
	return dir_count

def print_with_time_info(msg):
	pass#IMPORTIMPORTIMPORT from time import localtime, strftime
	line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>" + msg
	print(line)

def print_program_start_information():

	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from socket import gethostname

	myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	mpi_size = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if(myid == 0):
		print("Location: " + os.getcwd())

	print("MPI Rank: %03d/%03d "%(myid, mpi_size) + "Hostname: " + socket.gethostname() +  " proc_id: " + str(os.getpid()))



def store_program_state(filename, state, stack):
	pass#IMPORTIMPORTIMPORT import json
	with open(filename, "w") as fp:
		json.dump(list(zip(stack, state)), fp, indent = 2)
	fp.close()

def restore_program_stack_and_state(file_name_of_saved_state):
	pass#IMPORTIMPORTIMPORT import json
	f = open(file_name_of_saved_state, 'r')
	saved_state_and_stack = json.load(f); f.close()
	return list(zip(*saved_state_and_stack)[0]), list(zip(*saved_state_and_stack)[1])


def program_state_stack(full_current_state, frameinfo, file_name_of_saved_state=None, last_call="", force_starting_execution = False):

	"""

	pass#IMPORTIMPORTIMPORT When used it needs: from inspect import currentframe, getframeinfo
	pass#IMPORTIMPORTIMPORT Also: from utilities import program_state_stack

	This function is used for restarting time consuming data processing programs/steps from the last saved point.

	This static variable must be defined before the first call:
	program_state_stack.PROGRAM_STATE_VARIABLES = {"isac_generation", "i", "j"}
	It contains local variables at any level of the stack that define uniquely the state(flow/logic) of the program.

	It is assumed that the processed data is saved at each step and it is independent from the variables that uniquely define
	the state(flow/logic) of the program. All the variables that are used in more than one step must be calculated before
	the "if program_state_stack(locals(), getframeinfo(currentframe())):" call. It is assumed that they are not time consuming.
	Passing processed data from one step to the next is done only through files.

	First call needs to contain "file_name_of_saved_state".
	Then, the next calls are "if program_state_stack(locals(), getframeinfo(currentframe())):" to demarcate the blocks of
	processing steps that take a long time (hours/days).

	Example of initialization:
	program_state_stack.PROGRAM_STATE_VARIABLES = {"isac_generation", "i", "j"}
	program_state_stack(locals(), getframeinfo(currentframe()), "my_state.json")

	Then regular usage in the program:

	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass

	"""

	pass#IMPORTIMPORTIMPORT from traceback import extract_stack
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	pass#IMPORTIMPORTIMPORT from utilities import if_error_then_all_processes_exit_program
	pass#IMPORTIMPORTIMPORT import os

	def get_current_stack_info():
		return [[x[0], x[2]] for x in traceback.extract_stack()[:-2]]

	START_EXECUTING_FALSE = 0
	START_EXECUTING_TRUE = 1
	START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT = 2

	# error_status = 1
	# if_error_then_all_processes_exit_program(error_status)

	current_state = dict()
	for var in program_state_stack.PROGRAM_STATE_VARIABLES & set(full_current_state) :
		current_state[var] =  full_current_state[var]

	if "restart_location_title" in program_state_stack.__dict__:
		# location_in_program = frameinfo.filename + "___" + program_state_stack.restart_location_title + "___" + last_call
		location_in_program = frameinfo.filename + "___" + program_state_stack.restart_location_title
		del program_state_stack.restart_location_title
	else:
		location_in_program = frameinfo.filename + "___" + str(frameinfo.lineno) + "_" + last_call
		# location_in_program = frameinfo.filename + "___" + last_call

	current_state["location_in_program"] = location_in_program

	current_stack = get_current_stack_info()

	error_status = 0

	# not a real while, an if with the possibility of jumping with break
	while mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD) == 0:
		if "file_name_of_saved_state" not in program_state_stack.__dict__:
			if type(file_name_of_saved_state) != type(""):
				print("Must provide the file name of saved state as a string in the first call of the function!")
				error_status = 1
				break

			program_state_stack.file_name_of_saved_state = os.getcwd() + os.sep + file_name_of_saved_state
			program_state_stack.counter = 0
			program_state_stack.track_stack = get_current_stack_info()
			program_state_stack.track_state = [dict() for i in range(len(program_state_stack.track_stack))]
			program_state_stack.track_state[-1] = current_state

			file_name_of_saved_state_contains_information = False
			if (os.path.exists(file_name_of_saved_state)):
				statinfo = os.stat(file_name_of_saved_state)
				file_name_of_saved_state_contains_information = statinfo.st_size > 0
			if file_name_of_saved_state_contains_information:
				program_state_stack.saved_stack, \
				program_state_stack.saved_state = restore_program_stack_and_state(file_name_of_saved_state)
				program_state_stack.start_executing = START_EXECUTING_FALSE
			else:
				# check to see if file can be created
				f = open(file_name_of_saved_state, "w"); f.close()
				program_state_stack.start_executing = START_EXECUTING_TRUE
		else:
			program_state_stack.counter += 1
			# print "counter: ", program_state_stack.counter
			# if program_state_stack.counter == program_state_stack.CCC:
			# 	error_status = 1
			# 	break

			if program_state_stack.start_executing == START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT:
				program_state_stack.start_executing = START_EXECUTING_FALSE

			# correct track_state to reflect track_stack
			for i in range(len(current_stack)):
				if i < len(program_state_stack.track_state):
					if program_state_stack.track_stack[i] != current_stack[i]:
						program_state_stack.track_state[i] = dict()
				else:
					# print "i:", i, len(program_state_stack.track_state), len(current_stack), current_stack
					program_state_stack.track_state.append(dict())
			program_state_stack.track_state[i] = current_state

			# correct track_stack to reflect current_stack
			program_state_stack.track_stack = current_stack

			# if program_state_stack.counter == 68:
			# 	print range(len(current_stack), len(program_state_stack.track_state))

			# delete additional elements in track_state so that size of track_state is the same as current_stack
			program_state_stack.track_state[len(current_stack):len(program_state_stack.track_state)] = []

			if program_state_stack.start_executing == START_EXECUTING_TRUE or last_call != "" or force_starting_execution:
				store_program_state(program_state_stack.file_name_of_saved_state, program_state_stack.track_state, current_stack)
				program_state_stack.start_executing = START_EXECUTING_TRUE
			else:
				if len(program_state_stack.saved_state) >= len(current_stack):
					for i in range(len(program_state_stack.saved_state)):
						if i < len(current_stack):
							if program_state_stack.track_stack[i] == current_stack[i]:
								if program_state_stack.track_state[i] == program_state_stack.saved_state[i]:
									continue
							break
						else:
							program_state_stack.start_executing = START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT
							# print "////////////////////////////"
							# print "Entering function: ", location_in_program
							# print "////////////////////////////"
							break
					else:
						program_state_stack.start_executing = START_EXECUTING_TRUE
						# print "////////////////////////////"
						# print "Start executing: ", location_in_program
						# print "////////////////////////////"
		break
	else:
		program_state_stack.start_executing = START_EXECUTING_FALSE

	if_error_then_all_processes_exit_program(error_status)

	program_state_stack.start_executing = mpi.mpi_bcast(program_state_stack.start_executing, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
	program_state_stack.start_executing = int(program_state_stack.start_executing[0])

	# print "program_state_stack.start_executing ", program_state_stack.start_executing

	return program_state_stack.start_executing

def qw(s):
	s = s.replace("\n"," ")
	s = s.replace("\t"," ")
	return tuple(s.split())

def list_prod(list_whose_elements_are_going_to_be_multiplied):
	pass#IMPORTIMPORTIMPORT import operator
	return reduce(operator.mul, list_whose_elements_are_going_to_be_multiplied)

def calculate_space_size(x_half_size, y_half_size, psi_half_size):
	return [2 * x_half_size + 1, 2 * y_half_size + 1, 2 * psi_half_size + 1]


def print_from_process(process_rank, message):
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	pass#IMPORTIMPORTIMPORT import os, sys
	pass#IMPORTIMPORTIMPORT from socket import gethostname

	myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	mpi_size = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if(myid == process_rank):
		print("MPI Rank: %03d/%03d "%(myid, mpi_size) + "Hostname: " + socket.gethostname() +  " proc_id: " + str(os.getpid()) +\
			"message:::", message)
	sys.stdout.flush()

def mpi_exit():
	pass#IMPORTIMPORTIMPORT from mpi import mpi_finalize
	pass#IMPORTIMPORTIMPORT import sys

	mpi.mpi_finalize()
	sys.stdout.flush()
	sys.exit()

### from sort3d

def get_attr_stack(data_stack,attr_string):
	attr_value_list = []
	for idat in range(len(data_stack)):
		attr_value = data_stack[idat].get_attr(attr_string)
		attr_value_list.append(attr_value)
	return attr_value_list

def get_sorting_params(Tracker,data):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from utilities import read_text_row,wrap_mpi_bcast,even_angles
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	nproc     = Tracker["constants"]["nproc"]
	ndata     = Tracker["total_stack"]
	mpi_comm  = mpi.MPI_COMM_WORLD
	if myid == main_node:
		total_attr_value_list = []
		for n in range(ndata):
			total_attr_value_list.append([])
	else:
		total_attr_value_list = 0
	for inode in range(nproc):
		attr_value_list = get_attr_stack(data,"group")
		attr_value_list = wrap_mpi_bcast(attr_value_list,inode)
		if myid == main_node:
			image_start,image_end = applications.MPI_start_end(ndata,nproc,inode)
			total_attr_value_list = fill_in_mpi_list(total_attr_value_list,attr_value_list,image_start,image_end)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	total_attr_value_list = wrap_mpi_bcast(total_attr_value_list,main_node)
	return total_attr_value_list

def get_groups_from_partition_deprecated(partition, initial_ID_list, number_of_groups):
	# sort out Kmref results to individual groups that has initial IDs
	# make a dictionary
	dict = {}
	for iptl in range(len(initial_ID_list)):
		dict[iptl] = initial_ID_list[iptl]
	res = []
	for igrp in range(number_of_groups):
		class_one = []
		for ipt in range(len(partition)):
			if partition[ipt] == igrp:
				orginal_id = dict[ipt]
				class_one.append(orginal_id)
		res.append(class_one)
	return res

def remove_small_groups_deprecated_1(class_list,minimum_number_of_objects_in_a_group):
	new_class  = []
	final_list = []
	for one_class in class_list:
		if len(one_class)>=minimum_number_of_objects_in_a_group:
			new_class.append(one_class)
			for element in one_class:
				final_list.append(element)
	final_list.sort()
	return final_list, new_class

#### Used in the main programm

def remove_small_groups_deprecated_2(class_list,minimum_number_of_objects_in_a_group):
	new_class  = []
	final_list = []
	for one_class in class_list:
		if len(one_class)>=minimum_number_of_objects_in_a_group:
			new_class.append(one_class)
			for element in one_class:
				final_list.append(element)
	final_list.sort()
	return final_list, new_class

def print_a_line_with_timestamp_deprecated(string_to_be_printed ):
	line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
	print((line,string_to_be_printed))
	return string_to_be_printed

def get_outliers(total_number,plist):
	tlist={}
	for i in range(total_number):tlist[i]=i
	for a in plist:   del tlist[a]
	out =[]
	for a in tlist:   out.append(a)
	return out

def get_margin_of_error(this_group_of_data,Tracker):
	ratio = margin_of_error(Tracker["P_chunk0"],len(this_group_of_data))
	rate1, rate2, size_of_this_sampling = count_chunk_members(Tracker["chunk_dict"],this_group_of_data)
	return abs(rate1-Tracker["P_chunk0"]),ratio,abs(rate2-Tracker["P_chunk1"]),ratio

def get_ali3d_params(ali3d_old_text_file,shuffled_list):
	pass#IMPORTIMPORTIMPORT from utilities import read_text_row
	ali3d_old = read_text_row(ali3d_old_text_file)
	ali3d_new = []
	for iptl in range(len(shuffled_list)):
		ali3d_new.append(ali3d_old[shuffled_list[iptl]])
	return ali3d_new

def get_number_of_groups_deprecated(total_particles,number_of_images_per_group, round_off=.2):
	number_of_groups=float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<round_off:
		number_of_groups = int(number_of_groups)
	else:
		number_of_groups = int(number_of_groups)+1
	return number_of_groups

def get_number_of_groups_deprecated_2(total_particles,number_of_images_per_group):
	# soft partition groups
	number_of_groups=float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<.4:
		number_of_groups = int(number_of_groups)
	else:
		number_of_groups = int(number_of_groups)+1
	return number_of_groups

def get_complementary_elements_total(total_stack, data_list):
	data_dict    ={}
	complementary     = []
	for index in range(len(data_list)):data_dict[data_list[index]]=index
	for index in range(total_stack):
		if (index in data_dict) is False:complementary.append(index)
	return complementary

def get_two_chunks_from_stack(Tracker):
	total_chunk = EMAN2_cppwrap.EMUtil.get_all_attributes(Tracker["orgstack"],"chunk_id")
	chunk_one = []
	chunk_two = []
	for index_of_total_chunk in range(len(total_chunk)):
		if total_chunk[index_of_total_chunk]==0:chunk_one.append(index_of_total_chunk)
		else:chunk_two.append(index_of_total_chunk)
	return chunk_one, chunk_two

def set_filter_parameters_from_adjusted_fsc(n1,n2,Tracker):
	fsc_cutoff   = 1.0/3.0
	adjusted_fsc = adjust_fsc_down(Tracker["global_fsc"],n1,n2)
	currentres   = -1.0
	ns           = len(adjusted_fsc)
	for i in range(1,ns-1):
		if adjusted_fsc[1][i] < fsc_cutoff:
			currentres = adjusted_fsc[0][i-1]
			break
	lowpass, falloff    = filter.fit_tanh1(adjusted_fsc, 0.01)
	lowpass             = round(lowpass,4)
	falloff             = min(.1,falloff)
	falloff             = round(falloff,4)
	currentres          = round(currentres,2)
	Tracker["lowpass"]  = lowpass
	Tracker["falloff"]  = falloff
##### from RSORT

def get_class_members(sort3d_dir):
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from utilities import read_text_file
	maximum_generations = 100
	maximum_groups      = 100
	class_list = []
	igen = 0
	while( igen < maximum_generations ):
		gendir = os.path.join(sort3d_dir, "generation%03d"%igen)
		if os.path.exists(gendir):
			igrp = 0
			while( igrp < maximum_groups):
				Class_file = os.path.join(gendir, "Kmref/Class%d.txt"%igrp)
				if os.path.exists(Class_file):
					class_one = read_text_file(Class_file)
					class_list.append(class_one)
					igrp += 1
				else:
					igrp = maximum_groups
			igen += 1
		else:
			igen = maximum_generations
	return class_list

def get_stable_members_from_two_runs(SORT3D_rootdirs, ad_hoc_number, log_main):
	#SORT3D_rootdirs                       =sys.argv[1]
	# ad_hoc_number would be a number larger than the id simply for handling two_way comparison of non-equal number of groups from two partitions.
	########
	pass#IMPORTIMPORTIMPORT from string import split
	pass#IMPORTIMPORTIMPORT from statistics import k_means_match_clusters_asg_new
	pass#IMPORTIMPORTIMPORT from numpy import array

	sort3d_rootdir_list = SORT3D_rootdirs.split()
	dict1              = []
	maximum_elements   = 0
	for index_sort3d in range(len(sort3d_rootdir_list)):
		sort3d_dir       = sort3d_rootdir_list[index_sort3d]
		all_groups       = get_class_members(sort3d_dir)
		dict1.append(all_groups)
		if maximum_elements <len(all_groups):
			maximum_elements = len(all_groups)
	TC = ad_hoc_number + 1
	for indep in range(len(dict1)):
		alist = dict1[indep]
		while len(alist)<maximum_elements:
			alist.append([TC])
			TC += 1
		dict1[indep] = alist
		TC += 1
	for a in dict1:   log_main.add(len(a))
	dict = {}
	for index_sort3d in range(len(sort3d_rootdir_list)):
		sort3d_dir       = sort3d_rootdir_list[index_sort3d]
		dict[sort3d_dir] = dict1[index_sort3d]
	###### Conduct two-way comparison
	for isort3d in range(0,1): #len(sort3d_rootdir_list)):
		li = dict[sort3d_rootdir_list[isort3d]]
		new_li = []
		for ili in range(len(li)):
			li[ili].sort()
			t= numpy.array(li[ili],'int32')
			new_li.append(t)
		avg_list = {}
		total    = {}
		for ii in range(len(li)):
			avg_list[ii]=0.0
			total[ii]=0.0
		for jsort3d in range(len(sort3d_rootdir_list)):
			if isort3d != jsort3d:
				new_lj = []
				lj = dict[sort3d_rootdir_list[jsort3d]]
				for a in lj:
					log_main.add("the size is  %d"%len(a))
				for jlj in range(len(lj)):
					lj[jlj].sort()
					t= numpy.array(lj[jlj],'int32')
					new_lj.append(t)
				ptp=[new_li,new_lj]
				newindeces, list_stable, nb_tot_objs = statistics.k_means_match_clusters_asg_new(ptp[0],ptp[1])
				log_main.add("*************************************************************")
				log_main.add("the results of two P1 runs are: ")
				for index in range(len(newindeces)):
					log_main.add("  %d of %s matches  %d of %s"%(newindeces[index][0],sort3d_rootdir_list[isort3d],newindeces[index][1],sort3d_rootdir_list[jsort3d]))
				for index in range(len(list_stable)):
					log_main.add("%d   stable memebers"%len(list_stable[index]))
				new_stable = []
				for ilist in range(len(list_stable)):
					if len(list_stable[ilist])!= 0:
						new_stable.append(list_stable[ilist])
				for istable in range(len(new_stable)):
					stable = new_stable[istable]
					if len(stable)>0:
						group_A =  li[newindeces[istable][0]]
						group_B =  lj[newindeces[istable][1]]
						log_main.add(" %d %d %d   "%(len(group_A),len(group_B),len(stable)))
		return new_stable

def two_way_comparison_single(partition_A, partition_B,Tracker):
	###############
	pass#IMPORTIMPORTIMPORT from statistics import k_means_match_clusters_asg_new
	pass#IMPORTIMPORTIMPORT from utilities import count_chunk_members, margin_of_error
	pass#IMPORTIMPORTIMPORT from numpy import array
	#two_way_comparison_single
	total_stack = Tracker["constants"]["total_stack"]
	log_main    = Tracker["constants"]["log_main"]
	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	numpy32_A   = []
	numpy32_B   = []
	total_A     = 0
	total_B     = 0
	###--------------------------------------

	if myid == main_node:
		log_main.add(" the first run has number of particles %d"%len(partition_A))
		log_main.add(" the second run has number of particles %d"%len(partition_B))
	for A in partition_A:
		total_A +=len(A)
	for B in partition_B:
		total_B +=len(B)
	nc_zero = 1
	if len(partition_A) < len(partition_B):
		while len(partition_A) <len(partition_B):
			partition_A.append([nc_zero+total_stack])
			nc_zero +=1
	elif len(partition_A) > len(partition_B):
		while len(partition_B) <len(partition_A):
			partition_B.append([nc_zero+total_stack])
			nc_zero +=1
	number_of_class = len(partition_A)
	for index_of_class in range(number_of_class):
		A = partition_A[index_of_class]
		A.sort()
		A= numpy.array(A,'int32')
		numpy32_A.append(A)
		B= partition_B[index_of_class]
		B.sort()
		B = numpy.array(B,'int32')
		numpy32_B.append(B)
		if myid ==main_node:
			log_main.add("group %d  %d   %d"%(index_of_class,len(A), len(B)))
	ptp    = [[],[]]
	ptp[0] = numpy32_A
	ptp[1] = numpy32_B
	newindexes, list_stable, nb_tot_objs = statistics.k_means_match_clusters_asg_new(ptp[0],ptp[1])
	if myid == main_node:
		log_main.add(" reproducible percentage of the first partition %f"%(nb_tot_objs/float(total_A)*100.))
		log_main.add(" reproducible percentage of the second partition %f"%(nb_tot_objs/float(total_B)*100.))
		for index in range(len(newindexes)):
			log_main.add("%d of A match %d of B "%(newindexes[index][0],newindexes[index][1]))
		for index in range(len(list_stable)):
			log_main.add("%d number of reproduced objects are found in group %d"%(len(list_stable[index]),index))
		log_main.add(" %d number of objects are reproduced "%nb_tot_objs)
		log_main.add(" margin of error")
	large_stable = []
	for index_of_stable in range(len(list_stable)):
		rate1,rate2,size_of_this_group = count_chunk_members(Tracker["chunk_dict"], list_stable[index_of_stable])
		if size_of_this_group>=Tracker["constants"]["smallest_group"]:
			error                          = margin_of_error(Tracker["P_chunk0"],size_of_this_group)
			if myid == main_node:
				log_main.add(" chunk0  lower bound %f  upper bound  %f  for sample size  %d     group id %d"%((Tracker["P_chunk0"]- error),(Tracker["P_chunk0"]+error),size_of_this_group, index_of_stable))
				log_main.add(" actual percentage is %f"%rate1)
			large_stable.append(list_stable[index_of_stable])
		else:
			if myid==main_node:
				log_main.add("%d  group is too small"%index_of_stable)
	return large_stable

def get_leftover_from_stable(stable_list, N_total, smallest_group):
	tmp_dict = {}
	for i in range(N_total):
		tmp_dict[i] = i
	new_stable      =[]
	for alist in stable_list:
		if len(alist) > smallest_group:
			for index_of_list in range(len(alist)):
				del tmp_dict[alist[index_of_list]]
			new_stable.append(alist)
	leftover_list = []
	for one_element in tmp_dict:
		leftover_list.append(one_element)
	return leftover_list, new_stable

def Kmeans_exhaustive_run(ref_vol_list,Tracker):
	pass#IMPORTIMPORTIMPORT from applications import ali3d_mref_Kmeans_MPI
	pass#IMPORTIMPORTIMPORT from utilities import write_text_file
	pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_two_chunks_MPI
	pass#IMPORTIMPORTIMPORT from morphology import get_shrink_3dmask
	pass#IMPORTIMPORTIMPORT from utilities import wrap_mpi_bcast
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_barrier
	# npad 2 ---------------------------------------
	npad                  = 2
	myid                  = Tracker["constants"]["myid"]
	main_node             = Tracker["constants"]["main_node"]
	log_main              = Tracker["constants"]["log_main"]
	nproc                 = Tracker["constants"]["nproc"]
	final_list_text_file  = Tracker["this_data_list_file"] ## id text file for get_shrink_data_huang
	snr                   = 1.
	Tracker["total_stack"]= len(Tracker["this_data_list"])
	if myid == main_node:
		log_main.add("start exhaustive Kmeans")
		log_main.add("total data is %d"%len(Tracker["this_data_list"]))
		log_main.add("final list file is "+final_list_text_file)
	workdir = Tracker["this_dir"]
	####----------------------------------------------
	empty_group = 1
	kmref       = 0
	while empty_group ==1 and kmref<=5:## In case pctn of Kmeans jumps between 100% to 0%, stop the program
		if myid ==main_node: log_main.add(" %d     Kmref run"%kmref)
		outdir =os.path.join(workdir, "Kmref%d"%kmref)
		empty_group, res_classes, data_list = applications.ali3d_mref_Kmeans_MPI(ref_vol_list, outdir, final_list_text_file, Tracker)
		kmref +=1
		if empty_group ==1:
			if myid ==main_node:
				log_main.add("empty gorup appears, next round of Kmeans requires rebuilding reference volumes!")
				log_main.add(" the number of classes for next round before cleaning is %d"%len(res_classes))
			final_list   = []
			new_class    = []
			for a in res_classes:
				if len(a)>=Tracker["constants"]["smallest_group"]:
					for b in a:
						final_list.append(b)
					new_class.append(a)
			final_list.sort()
			### reset variables of Kmeans run
			Tracker["total_stack"]    = len(final_list)
			Tracker["this_data_list"] = final_list
			final_list_text_file      = os.path.join(workdir, "final_list%d.txt"%kmref)
			if myid == main_node:
				log_main.add("number of classes for next round is %d"%len(new_class))
				write_text_file(final_list, final_list_text_file)
				number_of_ref_class = []
				for igrp in range(len(new_class)):
					write_text_file(new_class[igrp],os.path.join(workdir,"final_class%d.txt"%igrp))
					number_of_ref_class.append(len(new_class[igrp]))
			else:
				number_of_ref_class = 0
			number_of_ref_class = wrap_mpi_bcast(number_of_ref_class,main_node)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			ref_vol_list = []
			if  Tracker["constants"]["mask3D"]: mask3D = morphology.get_shrink_3dmask(Tracker["constants"]["nxinit"],Tracker["constants"]["mask3D"])
			else: mask3D = None
			Tracker["number_of_ref_class"] = number_of_ref_class
			for igrp in range(len(new_class)):
				data,old_shifts = get_shrink_data_huang(Tracker,Tracker["nxinit"],os.path.join(workdir,"final_class%d.txt"%igrp),Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
				#volref = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo=None)
				#volref = filt_tanl(volref, Tracker["low_pass_filter"],.1)
				volref, fsc_kmref = reconstruction.rec3D_two_chunks_MPI(data,snr,Tracker["constants"]["sym"],mask3D,\
			 os.path.join(outdir, "resolution_%02d_Kmref%04d"%(igrp,kmref)), myid, main_node, index=-1, npad=npad, finfo = None)
				if myid !=main_node:
					volref = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
				bcast_EMData_to_all(volref, myid, main_node, mpi.MPI_COMM_WORLD)
				ref_vol_list.append(volref)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		else:
			new_class    = []
			for a in res_classes:
				if len(a)>=Tracker["constants"]["smallest_group"]:new_class.append(a)
	if myid == main_node:
		log_main.add("Exhaustive Kmeans finishes")
		log_main.add(" %d groups are selected out"%len(new_class))
	return new_class

def print_a_line_with_timestamp(string_to_be_printed ):
	line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
	print((line,string_to_be_printed))
	return string_to_be_printed

def split_a_group(workdir,list_of_a_group,Tracker):
	### Using EQ-Kmeans and Kmeans to split a group
	pass#IMPORTIMPORTIMPORT from utilities import wrap_mpi_bcast
	pass#IMPORTIMPORTIMPORT from random import shuffle
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from utilities import get_shrink_data_huang
	pass#IMPORTIMPORTIMPORT from reconstructions import recons3d_4nn_ctf_MPI
	pass#IMPORTIMPORTIMPORT from filter import filt_tanl
	pass#IMPORTIMPORTIMPORT from applications import mref_ali3d_EQ_Kmeans
	################
	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	nproc       = Tracker["constants"]["nproc"]
	total_stack = len(list_of_a_group)
	################
	pass#IMPORTIMPORTIMPORT import copy
	data_list [:] = list_of_a_group[:]
	update_full_dict(data_list,Tracker)
	this_particle_text_file = os.path.join(workdir,"full_class.txt")
	if myid ==main_node: write_text_file(data_list,"full_class.txt")
	# Compute the resolution of leftover
	if myid == main_node:
		random.shuffle(data_list)
		l1=data_list[0:total_stack//2]
		l2=data_list[total_stack//2:]
		l1.sort()
		l2.sort()
	else:
		l1 = 0
		l2 = 0
	l1 = wrap_mpi_bcast(l1, main_node)
	l2 = wrap_mpi_bcast(l2, main_node)
	llist =[l1,l2]
	if myid ==main_node:
		for index in range(2):
			partids = os.path.join(workdir,"Class_%d.txt"%index)
			write_text_file(llist[index],partids)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	################ create references for EQ-Kmeans
	ref_list = []
	for index in range(2):
		partids = os.path.join(workdir,"Class_%d.txt"%index)
		while not os.path.exists(partids):
			#print  " my_id",myid
			time.sleep(2)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"],partids,Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
		vol = reconstruction.recons3d_4nn_ctf_MPI(myid=myid,prjlist = data,symmetry=Tracker["constants"]["sym"],finfo=None)
		vol = filter.filt_tanl(vol,Tracker["constants"]["low_pass_filter"],.1)
		ref_list.append(vol)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	### EQ-Kmeans
	outdir = os.path.join(workdir,"EQ-Kmeans")
	applications.mref_ali3d_EQ_Kmeans(ref_list,outdir,this_particle_text_file,Tracker)
	res_EQ = partition_to_groups(Tracker["this_partition"],K=2)
	new_class = []
	for index in range(len(res_EQ)):
		new_ID = get_initial_ID(res_EQ(index), Tracker["full_ID_dict"])
		new_class.append(new_ID)
		if myid ==main_node:
			new_class_file = os.path.join(workdir,"new_class%d.txt"%index)
			write_text_file(new_ID,new_class_file)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	############# create references for Kmeans
	ref_list = []
	for index in range(2):
		partids = os.path.join(workdir,"new_class%d.txt"%index)
		while not os.path.exists(partids):
			#print  " my_id",myid
			time.sleep(2)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"],partids,Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
		vol = reconstruction.recons3d_4nn_ctf_MPI(myid=myid,prjlist = data,symmetry=Tracker["constants"]["sym"],finfo=None)
		vol = filter.filt_tanl(vol,Tracker["constants"]["low_pass_filter"],.1)
		ref_list.append(vol)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#### Kmeans

def search_lowpass(fsc):
	fcutoff =.5
	for i in range(len(fsc[1])):
		if fsc[1][i]<.5:
			break
	if i<len(fsc[1])-1:
		fcutoff=fsc[0][i-1]
	else:
		fcutoff=.5
	fcutoff=min(.45,fcutoff)
	return fcutoff


def split_chunks_bad(l, n):
	"""
	   Splits list l into n chunks with approximately equals sum of values
	   see  http://stackoverflow.com/questions/6855394/splitting-list-in-chunks-of-balanced-weight
	"""
	result = [[] for i in range(n)]
	sums   = {i:0 for i in range(n)}
	c = 0
	for e in l:
		for i in sums:
			if c == sums[i]:
				result[i].append(e)
				break
		sums[i] += e
		c = min(sums.values())
	for i in range(len(result)):
		result[i].sort()
	return result
	
def convert_to_float(value):
	"""
	When one wants to pass floats from C to python as integers, in C one would do
	float foo = 7.001e-23;
	unsigned int ival = *((unsigned int *)&foo);
	std::cout << ival<<std::endl;
	float goo = *((float *)&ival);
	std::cout << goo<<std::endl;
	return ival, and then in python covert it to float
	"""
	pass#IMPORTIMPORTIMPORT from struct import unpack
	return struct.unpack("!f", hex(value)[2:].zfill(8).decode('hex'))[0]

