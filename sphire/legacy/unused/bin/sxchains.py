"""1
def Plot(city, R, dist):
    # Plot
    Pt = [R[city[i]] for i in range(len(city))]
    Pt += [R[city[0]]]
    Pt = array(Pt)
    title('Total distance='+str(dist))
    plot(Pt[:,0], Pt[:,1], '-o')
    show()
"""
"""2
		# will align anyway
		try:
			ttt = d[0].get_attr('xform.params2d')
			for i in xrange(len(d)):
				alpha, sx, sy, mirror, scale = get_params2D(d[i])
				d[i] = rot_shift2D(d[i], alpha, sx, sy, mirror)
		except:
			pass
		"""
"""3
		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			print " %3d   %7.1f   %7.1f   %7.1f   %2d  %6.2f"%(snake[m][0], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"], snake[m][2])
		"""
"""4
		trms = []
		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			trms.append([prms["alpha"], prms["mirror"]])
		for i in xrange(3):
			for m in xrange(lend):
				mb = (m-1)%lend
				me = (m+1)%lend
				#  angles order mb,m,me
				# calculate predicted angles mb->m 
		"""
"""5
		#  This was an effort to get number of loops, inconclusive, to say the least
		pass#IMPORTIMPORTIMPORT from numpy import outer, zeros, float32, sqrt
		lend = len(d)
 		cor = zeros(lend,float32)
 		cor = outer(cor, cor)
		for i in xrange(lend):  cor[i][i] = 1.0
		for i in xrange(lend-1):
			for j in xrange(i+1, lend):
				cor[i,j] = lccc[mono(i,j)][0]
				cor[j,i] = cor[i,j]

		lmbd, eigvec = pca(cor)

		pass#IMPORTIMPORTIMPORT from utilities import write_text_file

		nvec=20
		print  [lmbd[j] for j in xrange(nvec)]
		print  " G"
		mm = [-1]*lend
		for i in xrange(lend):  # row
			mi = -1.0e23
			for j in xrange(nvec):
				qt = eigvec[j][i]
				if(abs(qt)>mi):
					mi = abs(qt)
					mm[i] = j
			for j in xrange(nvec):
				qt = eigvec[j][i]
				print  round(qt,3),   #  eigenvector
			print  mm[i]
		print
		for j in xrange(nvec):
			qt = []
			for i in xrange(lend):
				if(mm[i] == j):  qt.append(i)
			if(len(qt)>0):  write_text_file(qt,"loop%02d.txt"%j)
		"""
"""6
		print  [lmbd[j] for j in xrange(nvec)]
		print  " B"
		mm = [-1]*lend
		for i in xrange(lend):  # row
			mi = -1.0e23
			for j in xrange(nvec):
				qt = eigvec[j][i]/sqrt(lmbd[j])
				if(abs(qt)>mi):
					mi = abs(qt)
					mm[i] = j
			for j in xrange(nvec):
				qt = eigvec[j][i]/sqrt(lmbd[j])
				print  round(qt,3),   #  eigenvector
			print  mm[i]
		print
		"""
"""7
		lend=3
 		cor = zeros(lend,float32)
 		
 		cor = outer(cor, cor)
 		
 		
 		cor[0][0] =136.77
 		cor[0][1] = 79.15
 		cor[0][2] = 37.13
 		
 		cor[1][0] = 79.15
 		cor[2][0] = 37.13
 		
 		
 		cor[1][1] = 50.04
 		cor[1][2] = 21.65
 		
 		cor[2][1] = 21.65
 		
 		
 		cor[2][2] = 13.26

		lmbd, eigvec = pca(cor)
		print  lmbd
		print  eigvec
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i],   #  eigenvector
			print
		print  " B"
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i]/sqrt(lmbd[j]),   #  eigenvector
			print
		print  " G"
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i]*sqrt(lmbd[j]),   #  eigenvector
			print
		"""
def pca(cov):
	pass#IMPORTIMPORTIMPORT from numpy import  linalg, argsort
	""" assume one sample per column """
	values, vecs = numpy.linalg.eigh(cov)
	perm = numpy.argsort(-values)  # sort in descending order
	return values[perm], vecs[:, perm]


