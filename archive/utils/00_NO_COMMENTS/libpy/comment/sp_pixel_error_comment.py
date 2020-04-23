

















































































































































































































































































































































































































































































































































































































'''0
	def dfunc(args, data):

	        from math import pi, sin, cos
	        from numpy import zeros, array, float64

	        g = zeros(args.shape, float64)

	        ali_params = data[0]
	        d = data[1]

	        L = len(ali_params)
	        N = len(ali_params[0])/4

	        args_list= [0.0]*(L*3)
	        for i in xrange(L*3-3):        args_list[i] = args[i]
	        cosa = [0.0]*L
	        sina = [0.0]*L
	        for i in xrange(L):
	        	cosa[i] = cos(args_list[i*3]*pi/180.0)
	        	sina[i] = sin(args_list[i*3]*pi/180.0)
	        for i in xrange(N):
	        	sum_cosa = 0.0
	        	sum_sina = 0.0
	        	sx = [0.0]*L
	        	sy = [0.0]*L
	        	for j in xrange(L):
	        		if int(ali_params[j][i*4+3]) == 0:
	        			sum_cosa += cos((args_list[j*3]+ali_params[j][i*4])*pi/180.0)
	        			sum_sina += sin((args_list[j*3]+ali_params[j][i*4])*pi/180.0)
	        			sx[j] = args_list[j*3+1]+ali_params[j][i*4+1]*cosa[j]+ali_params[j][i*4+2]*sina[j]
	        			sy[j] = args_list[j*3+2]-ali_params[j][i*4+1]*sina[j]+ali_params[j][i*4+2]*cosa[j]
	        		else:
	        			sum_cosa += cos((-args_list[j*3]+ali_params[j][i*4])*pi/180.0)
	        			sum_sina += sin((-args_list[j*3]+ali_params[j][i*4])*pi/180.0)
	        			sx[j] = -args_list[j*3+1]+ali_params[j][i*4+1]*cosa[j]-ali_params[j][i*4+2]*sina[j]
	        			sy[j] =  args_list[j*3+2]+ali_params[j][i*4+1]*sina[j]+ali_params[j][i*4+2]*cosa[j]
	        	P = sqrt(sum_cosa**2+sum_sina**2)
	        	sum_cosa /= P
	        	sum_sina /= P

	        	for j in xrange(L-1):
	        		# Original formula, useful for double-checking, DON'T DELETE!
	        		#g[j*3] += d*d/4.0*(-1.0)*0.5/P*(-2*sum_cosa*P*sin((args_list[j*3]+ali_params[j][i*4])*pi/180.0)+\
	        		#			      2*sum_sina*P*cos((args_list[j*3]+ali_params[j][i*4])*pi/180.0))*pi/180.0+\
	        		#      2.0*(sx[j]-ave(sx))*(-ali_params[j][i*4+1]*sin(args_list[j*3]*pi/180.0)-ali_params[j][i*4+2]*cos(args_list[j*3]*pi/180.0))*pi/180.0+\
	        		#      2.0*(sy[j]-ave(sy))*( ali_params[j][i*4+1]*cos(args_list[j*3]*pi/180.0)-ali_params[j][i*4+2]*sin(args_list[j*3]*pi/180.0))*pi/180.0
	        		dx = 2.0*(sx[j]-ave(sx))
	        		dy = 2.0*(sy[j]-ave(sy))
	        		if int(ali_params[j][i*4+3]) == 0:
	        			g[j*3] += (d*d/4.0*(sum_cosa*sin((args_list[j*3]+ali_params[j][i*4])*pi/180.0)-\
	        					sum_sina*cos((args_list[j*3]+ali_params[j][i*4])*pi/180.0))+\
	        					dx*(-ali_params[j][i*4+1]*sina[j]+ali_params[j][i*4+2]*cosa[j])+\
	        					dy*(-ali_params[j][i*4+1]*cosa[j]-ali_params[j][i*4+2]*sina[j]))*pi/180.0
	        			g[j*3+1] += dx
	        			g[j*3+2] += dy
	        		else:
	        			g[j*3] += (d*d/4.0*(-sum_cosa*sin((-args_list[j*3]+ali_params[j][i*4])*pi/180.0)+\
	        					    sum_sina*cos((-args_list[j*3]+ali_params[j][i*4])*pi/180.0))+\
	        					 dx*( ali_params[j][i*4+1]*sina[j]+ali_params[j][i*4+2]*cosa[j])+\
	        					 dy*(-ali_params[j][i*4+1]*cosa[j]+ali_params[j][i*4+2]*sina[j]))*pi/180.0
	        			g[j*3+1] += -dx
	        			g[j*3+2] += dy
	        g /= (N*L)

	        return g
	'''














































































































































































































































































































































































































































































































































































































'''1
def read_meridien_parameters():
	# this is simple example how to read all orientation parameters ("smear") produced by Meririen
	# Note this is just an example.  It works, but it is not really useful for anythin in particular.
	import json
	main = 21


	fin = open( os.path.join("main%03d"%main, "Tracker_%03d.json"%main),'r')
	Tracker = convert_json_fromunicode(json.load(fin))
	fin.close()

	refang = read_text_row(os.path.join("main%03d"%main,"refang.txt"))
	rshifts = read_text_row(os.path.join(os.path.join("main%03d"%main, "rshifts.txt")))
	params = []
	for iproc in range(2):
		chunk = read_text_file(os.path.join("main%03d"%main,"chunk_%1d_%03d.txt"%(iproc,main)))
		params_previous = read_text_row(os.path.join(os.path.join("main%03d"%(main-1), "params-chunk_%1d_%03d.txt"%(iproc,main-1))))

		cnt = True
		imi = -1
		myid = -1
		while(cnt):
			try:
				#if( myid < 100 ):
				myid += 1
				f = os.path.join("main%03d"%main,"oldparamstructure", "oldparamstructure_%1d_%03d_%03d.json"%(iproc,myid,main))
				fin = open(f,'r')
				sxprint(f)
				paramstructure = convert_json_fromunicode(json.load(fin))
				fin.close()
				sxprint(iproc,myid,len(paramstructure))
				for im in range(len(paramstructure)):
					imi += 1
					numbor = len(paramstructure[im][2])
					params.append([chunk[imi],[]])
					for norient in range(numbor):
						ipsiandiang = paramstructure[im][2][norient][0]/1000
						ishift      = paramstructure[im][2][norient][0]%1000
						prob        = paramstructure[im][2][norient][1]
						ipsi = ipsiandiang%100000
						iang = ipsiandiang/100000
						phi = refang[iang][0]
						theta = refang[iang][1]
						psi = refang[iang][2]+ipsi*Tracker["delta"]
						tx = int(round(params_previous[im][3])) + rshifts[ishift][0]
						ty = int(round(params_previous[im][4])) + rshifts[ishift][1]
						params[-1][1].append([phi,theta,psi,tx,ty,prob])
						#print(params[-1])

			except:
				cnt = False


	del params_previous, refang, rshift, chunk, paramstructure
	params_previous = read_text_row(os.path.join(os.path.join("main%03d"%(main-1), "params_%03d.txt"%(main-1))))

	params = sorted(params)
	opar = []
	for i,q in enumerate(params):
		for k in range(len(params[i][1])):
			opar.append([params[i][0], params[i][1][k][0], params[i][1][k][1], params[i][1][k][2], params[i][1][k][3], params[i][1][k][4], params[i][1][k][5]])

	write_text_row(opar,"junki.txt")
'''



'''2

def pixel_error_angle_sets(agls1, agls2, Threshold=1.0e23, r=1.0):
	"""
	  It will compute actual pixel errors using all five orientation parameters (including shifts)
	  However, orientation is found using only the angles.
		 
		 INCORRECT FOR SYMMETRY
	
	  Input: Two lists, i-th element of each list is either a list of the three Eulerian angles [[phi1, theta1, psi1], [phi2, theta2, psi2], ...]
	         as read by read_text_row(filename, "")
	         Or, the two lists can be a list of Transform objects. The two lists have the same length and the i-th element of one list is assumed to correspond to the i-th element of the
		 second list. 
		 Threshold is a float.
		 r is the radius of the object.
		 
	  Output: 1. Uses rotation_between_anglesets to find the overall 3D rotation between the two sets of Eulerian angles using second list as reference
	  	      2. The overall rotation found by rotation_between_anglesets is applied to the first list (agls1) 
		      3. Output is a list of lists: If the i-th corresponding pair of eulerian angles on agls2 and agls1 has pixel error (computed using max_3D_pixel_error) less than Threshold, then append the list 
		       [i, p], where p is the pixel error, into the output list.
	"""
	from sp_pixel_error import max_3D_pixel_error
	from sp_utilities   import read_text_file, rotation_between_anglesets
	import types
	
	N = len(agls1)
	if N != len(agls2):
		print 'Both lists must have the same length'
		return [-1]
	if N < 2:
		print 'At least two orientations are required in each list'
		return [-1]


	##############################################################################################
	# Find overall rotation between two angle sets, and then apply it to one of the angle sets

	# compute rotation between asg1 and asg2
	phi12,theta12,psi12 = rotation_between_anglesets(agls1, agls2)
	# apply rotation [phi12,theta12,psi12] to asg1
	t12 = Transform({'type':'spider','phi':phi12,'theta':theta12,'psi':psi12})
	agls12=[]

	# if agls1 is a list of list
	if type(agls1[0]) == types.ListType:
		for i in xrange(N):
			t1= Transform({'type':'spider','phi':agls1[i][0],'theta':agls1[i][1],'psi':agls1[i][2]})
			agls12.append(t1*t12)
	else: # must be list of transform objects
		for i in xrange(N):
			agls12.append(agls1[i]*t12)
	##############################################################################################
	# Compute pixel error for each entry of asg12 and asg2 
	# (asg13 and asg3, asg23 and asg3 respectively) and return a list of the pixel errors that are below a certain threshold

	# Compute average pixel error for each entry of asg12 and asg2 
	avgPixError12=[]
	if type(agls2[0]) == types.ListType:
		for i in xrange(N):
			error = max_3D_pixel_error(agls12[i],Transform({'type':'spider','phi':agls2[i][0],'theta':agls2[i][1],'psi':agls2[i][2]}) , r)
			if error < Threshold:
				avgPixError12.append([i,error])
	else:# agls2 is a list of transforms
		for i in xrange(N):
			error = max_3D_pixel_error(agls12[i], agls2[i] , r)
			if error < Threshold:
				avgPixError12.append([i,error])

	return avgPixError12,[phi12,theta12,psi12]


#  See symclass in fundamentals

def reduce_angles_sym(ang, sym = 'c1'):
	from sp_utilities import get_symt
	from EMAN2 import Vec2f, Transform, EMData
	ts = get_symt(sym)
	ks = len(ts)
	if(sym[0] == 'c'):
		if(ks == 1):  return
		dt = 360.0/int(sym[1:])
		for q in ang:
			if(q[1] <= 90.):  q[0] = q[0]%dt
			else:			  q[0] = (q[0]%360)%dt+180.0
	elif(sym[0] == 'd'):
		#for i in xrange(ks):  ts[i] = ts[i].inverse()
		dn = 360.0/(2*int(sym[1:]))
		if(int(sym[1:])%2 == 1):  # D odd
			ane = 360.0/int(sym[1:])/4
			ana = 2*ane
			ans = ane + 360.0/int(sym[1:])/2
			for q in ang:
				qt = Transform({"type":"spider","phi":q[0], "theta":q[1], "psi":q[2]})
				qt.set_trans(Vec2f(-q[3], -q[4]))
				#fifi = True
				for k in xrange(ks):
					ut = qt*ts[k]
					bt = ut.get_params("spider")
					tp = bt["phi"]
					tm = tp - 180.0
					tt = bt["theta"]
					if(((tp>=0.0 and tp<ane) or (tp>=ana and tp<ans) and tt <= 90.0) or (( tm>=0.0 and tm<ane) or (tm>=ana and tm<ans)and tt > 90.0)):
						#print  "%6d   %6.2f   %6.2f   %6.2f   %6.2f"%(k,q[0],q[1],tp,tt)
						q[0] = tp
						q[1] = tt
						q[2] = bt["psi"]
						q[3] = -bt["tx"]
						q[4] = -bt["ty"]
						#fifi = False
						break
				#if fifi:  print "FAILED     ", "%6d   %6.2f   %6.2f   %6.2f   %6.2f"%(k,q[0],q[1],tp,tt)
		else:  #  D even
			for q in ang:
				#print  #"%6d   %6.2f   %6.2f"%(12,q[0],q[1])
				qt = Transform({"type":"spider","phi":q[0], "theta":q[1], "psi":q[2]})
				qt.set_trans(Vec2f(-q[3], -q[4]))
				for k in xrange(ks):
					ut = qt*ts[k]
					bt = ut.get_params("spider")
					tp = round(bt["phi"],3)%360.0
					tm = tp - 180.0
					tt = round(bt["theta"],3)%360.0
					#print  "%6d   %6.2f   %6.2f   %6.2f   %6.2f"%(k,q[0],q[1],tp,tt)

					if( (tp >= 0.0 and tp<=dn and tt<= 90.0) or ( tm>=0.0 and tm<=dn and tt>90.0  )):
						#print  "%6d   %6.2f   %6.2f   %6.2f   %6.2f"%(k,q[0],q[1],tp,tt)
						q[0] = tp
						q[1] = tt
						q[2] = bt["psi"]
						q[3] = -bt["tx"]
						q[4] = -bt["ty"]
						break
	else:
		ERROR("Only C and D symmetries supported","reduce_angles_sym",1)

def apply_sym_angles(ang, sym = 'c1'):
	from sp_utilities import get_symt
	from EMAN2 import Vec2f, Transform, EMData
	na = len(ang)
	ts = get_symt(sym)
	ks = len(ts)
	angsa = [None]*ks*na
	if(sym == 'c1'):
		for i in xrange(na):
			angsa[i] = ang[i]
	else:
		la = 0
		for k in xrange(ks):
			for i in xrange(na):
				qt = Transform({"type":"spider","phi":ang[i][0], "theta":ang[i][1], "psi":ang[i][2]})
				qt.set_trans(Vec2f(-ang[i][3], -ang[i][4]))
				ut = qt*ts[k]
				bt = ut.get_params("spider")
				angsa[la] = [round(bt["phi"],3)%360.0, round(bt["theta"],3)%360.0, bt["psi"], -bt["tx"], -bt["ty"]]
				la += 1
	return  angsa


def multi_align_stability(ali_params, mirror_consistency_threshold = 0.75, error_threshold = 1.0, individual_error_threshold = 1.0, print_individual = False):

	def rot_shift(x, y, alpha, sx, sy):
		from math import pi, sin, cos
		cosi = cos(alpha/180.0*pi)
		sini = sin(alpha/180.0*pi)
		return x*cosi+y*sini+sx, -x*sini+y*cosi+sy

	def ave(a):
		n = len(a)
		ave = 0.0
		for i in xrange(n): ave += a[i]
		ave /= n
		return ave

	def avg_sqr_diff_from_mean(a):
		n = len(a)
		avg = ave(a)
		avg_sqr_diff_from_mean = 0.0
		for i in xrange(n): avg_sqr_diff_from_mean += (a[i]-avg)**2
		return avg_sqr_diff_from_mean/n

	def transform_variance(args, data):
		from math import sqrt
	        x1 = 1.0
	        y1 = 0.0
        	x2 = 0.0
	        y2 = 1.0

        	all_ali_params = data
	        num_ali = len(all_ali_params)
	        nima = len(all_ali_params[0])/4
		err = [0.0]*nima

	        for i in xrange(nima):
	        	pix_error = 0
	        	x1_new = [0.0]*num_ali
	        	y1_new = [0.0]*num_ali
	        	x2_new = [0.0]*num_ali
	        	y2_new = [0.0]*num_ali
	        	alpha2, sx2, sy2, mirror2 = all_ali_params[num_ali-1][i*4:i*4+4]
	        	x1_new[num_ali-1], y1_new[num_ali-1] = rot_shift(x1, y1, alpha2, sx2, sy2)
	        	x2_new[num_ali-1], y2_new[num_ali-1] = rot_shift(x2, y2, alpha2, sx2, sy2)
	        	for j in xrange(num_ali-1):
	        		alpha1, sx1, sy1, mirror1 = all_ali_params[j][i*4:i*4+4]
	        		alphai, sxi, syi = args[j*3:j*3+3]
	        		alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
	        		x1_new[j], y1_new[j] = rot_shift(x1, y1, alpha12, sx12, sy12)
	        		x2_new[j], y2_new[j] = rot_shift(x2, y2, alpha12, sx12, sy12)

	        	err[i] = sqrt(avg_sqr_diff_from_mean(x1_new)+avg_sqr_diff_from_mean(y1_new)+avg_sqr_diff_from_mean(x2_new)+avg_sqr_diff_from_mean(y2_new))

	        return err
		
	def transform_variance2(args, data):
		from math import sqrt
	        x1 = 30.0
	        y1 = 0.0
        	x2 = 0.0
	        y2 = 30.0

        	all_ali_params = data
	        num_ali = len(all_ali_params)
	        nima = len(all_ali_params[0])/4
		err = [0.0]*nima

		# Error introduced by rotation
	        for i in xrange(nima):
	        	pix_error = 0
	        	x1_new = [0.0]*num_ali
	        	y1_new = [0.0]*num_ali
	        	x2_new = [0.0]*num_ali
	        	y2_new = [0.0]*num_ali
	        	alpha2, sx2, sy2, mirror2 = all_ali_params[num_ali-1][i*4:i*4+4]
	        	x1_new[num_ali-1], y1_new[num_ali-1] = rot_shift(x1, y1, alpha2, 0.0, 0.0)
	        	x2_new[num_ali-1], y2_new[num_ali-1] = rot_shift(x2, y2, alpha2, 0.0, 0.0)
	        	for j in xrange(num_ali-1):
	        		alpha1, sx1, sy1, mirror1 = all_ali_params[j][i*4:i*4+4]
	        		alphai, sxi, syi = args[j*3:j*3+3]
	        		alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
	        		x1_new[j], y1_new[j] = rot_shift(x1, y1, alpha12, 0.0, 0.0)
	        		x2_new[j], y2_new[j] = rot_shift(x2, y2, alpha12, 0.0, 0.0)

	        	err[i] = sqrt(avg_sqr_diff_from_mean(x1_new)+avg_sqr_diff_from_mean(y1_new)+avg_sqr_diff_from_mean(x2_new)+avg_sqr_diff_from_mean(y2_new))

		# Error introduced by shift
        	all_ali_params = data
	        num_ali = len(all_ali_params)
	        nima = len(all_ali_params[0])/4
		err2 = [0.0]*nima

	        for i in xrange(nima):
	        	pix_error = 0
	        	x1_new = [0.0]*num_ali
	        	y1_new = [0.0]*num_ali
	        	x2_new = [0.0]*num_ali
	        	y2_new = [0.0]*num_ali
	        	alpha2, sx2, sy2, mirror2 = all_ali_params[num_ali-1][i*4:i*4+4]
	        	x1_new[num_ali-1], y1_new[num_ali-1] = rot_shift(x1, y1, 0.0, sx2, sy2)
	        	x2_new[num_ali-1], y2_new[num_ali-1] = rot_shift(x2, y2, 0.0, sx2, sy2)
	        	for j in xrange(num_ali-1):
	        		alpha1, sx1, sy1, mirror1 = all_ali_params[j][i*4:i*4+4]
	        		alphai, sxi, syi = args[j*3:j*3+3]
	        		alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
	        		x1_new[j], y1_new[j] = rot_shift(x1, y1, 0.0, sx12, sy12)
	        		x2_new[j], y2_new[j] = rot_shift(x2, y2, 0.0, sx12, sy12)

	        	err2[i] = sqrt(avg_sqr_diff_from_mean(x1_new)+avg_sqr_diff_from_mean(y1_new)+avg_sqr_diff_from_mean(x2_new)+avg_sqr_diff_from_mean(y2_new))

        	all_ali_params = data
	        num_ali = len(all_ali_params)
	        nima = len(all_ali_params[0])/4
		err3 = [0.0]*nima

	        for i in xrange(nima):
	        	pix_error = 0
	        	x1_new = [0.0]*num_ali
	        	y1_new = [0.0]*num_ali
	        	x2_new = [0.0]*num_ali
	        	y2_new = [0.0]*num_ali
	        	alpha2, sx2, sy2, mirror2 = all_ali_params[num_ali-1][i*4:i*4+4]
	        	x1_new[num_ali-1], y1_new[num_ali-1] = rot_shift(x1, y1, alpha2, sx2, sy2)
	        	x2_new[num_ali-1], y2_new[num_ali-1] = rot_shift(x2, y2, alpha2, sx2, sy2)
	        	for j in xrange(num_ali-1):
	        		alpha1, sx1, sy1, mirror1 = all_ali_params[j][i*4:i*4+4]
	        		alphai, sxi, syi = args[j*3:j*3+3]
	        		alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
	        		x1_new[j], y1_new[j] = rot_shift(x1, y1, alpha12, sx12, sy12)
	        		x2_new[j], y2_new[j] = rot_shift(x2, y2, alpha12, sx12, sy12)

	        	err3[i] = sqrt(avg_sqr_diff_from_mean(x1_new)+avg_sqr_diff_from_mean(y1_new)+avg_sqr_diff_from_mean(x2_new)+avg_sqr_diff_from_mean(y2_new))

	        return err, err2, err3

	from sp_statistics import k_means_stab_bbenum
	from sp_utilities import combine_params2
	from numpy import array

	# Find out the subset which is mirror consistent over all runs
	all_part = []
	num_ali = len(ali_params)
	nima = len(ali_params[0])/4
	for i in xrange(num_ali):
		mirror0 = []
		mirror1 = []
		for j in xrange(nima):
			if ali_params[i][j*4+3] == 0: mirror0.append(j)
			else: mirror1.append(j)
		mirror0 = array(mirror0, 'int32')
		mirror1 = array(mirror1, 'int32')
		all_part.append([mirror0, mirror1])
	match, stab_part, CT_s, CT_t, ST, st = k_means_stab_bbenum(all_part, T=0, nguesses=1)
	mir_cons_part = stab_part[0] + stab_part[1]
	mirror_consistent_rate = len(mir_cons_part)/float(nima)
	if mirror_consistent_rate <  mirror_consistency_threshold: return [], mirror_consistent_rate, -1.0
	mir_cons_part.sort()
	del all_part, match, stab_part, CT_s, CT_t, ST, st	

	# Shrink the list the alignment paramters to whatever is mirror consistent
	ali_params_mir_cons = []
	ali_params_mir_cons_list = []
	for i in xrange(num_ali):
		ali_params_temp = []
		for j in mir_cons_part: ali_params_temp.extend(ali_params[i][j*4:j*4+4])
		ali_params_mir_cons.append(ali_params_temp)
		ali_params_mir_cons_list.extend(ali_params_temp)

	# Find out the alignment parameters for each iteration against the last one
	args = []
	for i in xrange(num_ali-1):
		alpha, sx, sy, mirror, stab_mirror, pixel_error = ave_ali_err_params(ali_params_mir_cons[i], ali_params_mir_cons[num_ali-1])
		args.extend([alpha, sx, sy])

	ps = Util.multi_align_error(args, ali_params_mir_cons_list)
	val = ps[-1]
	del ps[-1]

	if val > error_threshold: return [], mirror_consistent_rate, val
	err, err2, err3 = transform_variance2(ps, ali_params_mir_cons)

	if print_individual:
		for i in xrange(nima):
			print "Particle %3d :"%i,
			if i in mir_cons_part:
				j = mir_cons_part.index(i)
				print "error = %6.4f  %6.4f  %6.4f"%(err[j], err2[j], err3[j]) 
			else: print "Mirror inconsistent"
	stable_set = []
	for i in xrange(len(mir_cons_part)):
		if err[i] < individual_error_threshold: stable_set.append([err[i], mir_cons_part[i]])
	stable_set.sort()
		
	return stable_set, mirror_consistent_rate, val


def estimate_stability(data1, data2, CTF=False, snr=1.0, last_ring=-1):
	"""
	This function estimate the stability of two datasets
	It returns three values, the first is the mirror consistent rate
	The second is the average pixel error among the mirror consistent images
	The third is the cross_correltion coefficient of two averages
	"""

	from sp_statistics import sum_oe, ccc
	from sp_fundamentals import fft, rot_shift2D
	from sp_alignment import align2d
	from sp_utilities import get_params2D, combine_params2
	from math import sin, pi, sqrt
	from sp_morphology import ctf_img

	PI_180 = pi/180
	nima = len(data1)
	nx = data1[0].get_xsize()
	if last_ring == -1: last_ring = nx/2-2
	if CTF:
		ctf_2_sum = EMData(nx, nx, 1, False)
		for im in xrange(nima):
			ctf_params = data1[im].get_attr("ctf")
			Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
		ctf_2_sum += 1/snr

	av1, av2 = sum_oe(data1, "a", CTF, EMData())
	if CTF:
		ave1 = fft(Util.divn_img(fft(Util.addn_img(av1, av2)), ctf_2_sum))
	else:
		ave1 = (av1+av2)/nima

	av1, av2 = sum_oe(data2, "a", CTF, EMData())
	if CTF:
		ave2 = fft(Util.divn_img(fft(Util.addn_img(av1, av2)), ctf_2_sum))
	else:
		ave2 = (av1+av2)/nima

	alpha21, sx21, sy21, mirror21, peak21 = align2d(ave2, ave1, 3.0, 3.0, 0.125, last_ring=last_ring)
	ave21 = rot_shift2D(ave2, alpha21, sx21, sy21, mirror21)
		
	consistent = 0
	pixel_error = []	
	for im in xrange(nima):
		alpha1, sx1, sy1, mirror1, scale1 = get_params2D(data1[im])
		alpha2, sx2, sy2, mirror2, scale2 = get_params2D(data2[im])

		alpha2n, sx2n, sy2n, mirror2n = combine_params2(alpha2, sx2, sy2, mirror2, alpha21, sx21, sy21, mirror21)

		if mirror1 == mirror2n:
			consistent += 1
			this_pixel_error = abs(sin((alpha1-alpha2n)*PI_180/2))*last_ring*2+sqrt((sx1-sx2n)**2+(sy1-sy2n)**2)
			pixel_error.append(this_pixel_error)

	return consistent/float(nima), pixel_error, ccc(ave21, ave1)



def max_3D_pixel_error(t1, t2, r):
	"""
	  Compute maximum pixel error between two projection directions
	  assuming object has radius r, t1 is the projection transformation
	  of the first projection and t2 of the second one, respectively:
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(-tx, -ty))
	  Note the function is symmetric in t1, t2.
	"""
	from math import sin, cos, pi, sqrt
	t3 = t2*t1.inverse()
	ddmax = 0.0
	for i in xrange(int(r), int(r)+1):
		for ang in xrange(int(2*pi*i+0.5)):
			v = Vec3f(i*cos(ang), i*sin(ang), 0)
			d = t3*v - v
			dd = d[0]**2+d[1]**2+d[2]**2
			if dd > ddmax: ddmax=dd
	return sqrt(ddmax)

def max_3D_pixel_errorA(t1, t2, r):
	"""
	  Compute maximum pixel error between two projection directions
	  assuming object has radius r, t1 is the projection transformation
	  of the first projection and t2 of the second one, respectively:
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(-tx, -ty))
	  Note the function is symmetric in t1, t2.
	"""
	from math import sin, cos, pi, sqrt
	t3 = t2*t1.inverse()
	ddmax = 0.0
	for i in xrange(int(r)+1):
		for ang in xrange(int(2*pi*i+0.5)):
			v = Vec3f(i*cos(ang), i*sin(ang), 0)
			d = t3*v - v
			dd = d[0]**2+d[1]**2+d[2]**2
			if dd > ddmax: ddmax=dd
	return sqrt(ddmax)
'''




