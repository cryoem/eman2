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
						# we abandon randomize phase strategy
						frc_without_mask = fsc(map1, map2, 1)
						randomize_at     = -1.0
						for ifreq in xrange(1, len(frc_without_mask[1])): # always skip zero frequency
							if frc_without_mask[1][ifreq] < options.randomphasesafter:
								randomize_at = float(ifreq)
								break
						log_main.add("Phases are randomized after: %4.2f[A]"% (options.pixel_size/(randomize_at/map1.get_xsize())))
						frc_masked = fsc(map1*m, map2*m, 1)
						map1 = fft(Util.randomizedphasesafter(fft(map1), randomize_at))*m
						map2 = fft(Util.randomizedphasesafter(fft(map2), randomize_at))*m
						frc_random_masked = fsc(map1, map2, 1)
						fsc_true          = [frc_without_mask[0], [None]*len(frc_without_mask[0])]
						for i in xrange(len(fsc_true[1])):
							if i < (int(randomize_at) + 2):# move two pixels up
								fsc_true[1][i] = frc_masked[1][i]
							else:
								fsct = frc_masked[1][i]
								fscn = frc_random_masked[1][i]
								if (fscn > fsct): fsc_true[1][i]= 0.
								else: fsc_true[1][i]=(fsct-fscn)/(1.-fscn)
						else:
					"""
"""3
							map1   = filt_tanl(map1,options.fl, min(options.aa,.1))
							cutoff = options.pixel_size/options.fl
							"""
'''4
			if options.comparison_radius>0.5:
				ERROR('Incorrect maximum resolution', 'options.subtract_stack',1)
			'''
'''5
							if options.comparison_radius !=-1:
								mask = model_circle(options.comparison_radius, image.get_xsize(), image.get_ysize())
							else:
							'''
'''6
					if options.maxres !=-1:
						temp_diff, a, b = im_diff(filt_tanl(image, options.maxres, options.maxresaa), simage, mask)						
					else:
						temp_diff, a, b = im_diff(image, simage, mask)
					image *=a
					image -=b
					'''
def Distance(i1, i2, lccc):
	return max(1.0 - lccc[statistics.mono(i1,i2)][0], 0.0)
# 	return sqrt((R1[0]-R2[0])**2+(R1[1]-R2[1])**2)

def TotalDistance(city, lccc):
	dist = 0.0
	for i in range(len(city)-1):
		dist += Distance(city[i], city[i+1], lccc)
	dist += Distance(city[-1], city[0], lccc)
	return dist

def reverse(city, n):
    nct = len(city)
    nn = (1+ ((n[1]-n[0]) % nct))/2 # half the lenght of the segment to be reversed
    # the segment is reversed in the following way n[0]<->n[1], n[0]+1<->n[1]-1, n[0]+2<->n[1]-2,...
    # Start at the ends of the segment and swap pairs of cities, moving towards the center.
    for j in range(nn):
        k = (n[0]+j) % nct
        l = (n[1]-j) % nct
        (city[k],city[l]) = (city[l],city[k])  # swap

def transpt(city, n):
    nct = len(city)

    newcity=[]
    # Segment in the range n[0]...n[1]
    for j in range( (n[1]-n[0])%nct + 1):
        newcity.append(city[ (j+n[0])%nct ])
    # is followed by segment n[5]...n[2]
    for j in range( (n[2]-n[5])%nct + 1):
        newcity.append(city[ (j+n[5])%nct ])
    # is followed by segment n[3]...n[4]
    for j in range( (n[4]-n[3])%nct + 1):
        newcity.append(city[ (j+n[3])%nct ])
    return newcity

"""Multiline Comment0"""

def tsp(lccc):

	#     ncity = 100        # Number of cities to visit
	pass#IMPORTIMPORTIMPORT from math import sqrt
	ncity = int( (1+numpy.sqrt(1+8*len(lccc)))/2 )        # Number of cities to visit
    #  sanity check
	if( ncity*(ncity-1)/2 != len(lccc) ): return [-1]

	maxTsteps = 100    # Temperature is lowered not more than maxTsteps
	Tstart = 0.2       # Starting temperature - has to be high enough
	fCool = 0.9        # Factor to multiply temperature at each cooling step
	maxSteps = 100*ncity     # Number of steps at constant temperature
	maxAccepted = 10*ncity   # Number of accepted steps at constant temperature

	Preverse = 0.5      # How often to choose reverse/transpose trial move


	# The index table -- the order the cities are visited.
	city = list(range(ncity))
	# Distance of the travel at the beginning
	dist = TotalDistance(city, lccc)

	#  Not clear what is nct
	nct = ncity
	# Stores points of a move
	n = numpy.zeros(6, dtype=int)

	T = Tstart # temperature

	#     Plot(city, R, dist)

	for t in range(maxTsteps):  # Over temperature

		accepted = 0
		for i in range(maxSteps): # At each temperature, many Monte Carlo steps

			while True: # Will find two random cities sufficiently close by
				# Two cities n[0] and n[1] are choosen at random
				n[0] = int((nct)*numpy.random.rand())     # select one city
				n[1] = int((nct-1)*numpy.random.rand())   # select another city, but not the same
				if (n[1] >= n[0]): n[1] += 1   #
				if (n[1] < n[0]): (n[0],n[1]) = (n[1],n[0]) # swap, because it must be: n[0]<n[1]
				nn = (n[0]+nct -n[1]-1) % nct  # number of cities not on the segment n[0]..n[1]
				if nn>=3: break

			# We want to have one index before and one after the two cities
			# The order hence is [n2,n0,n1,n3]
			n[2] = (n[0]-1) % nct  # index before n0  -- see figure in the lecture notes
			n[3] = (n[1]+1) % nct  # index after n2   -- see figure in the lecture notes

			if Preverse > numpy.random.rand():
				# Here we reverse a segment
				# What would be the cost to reverse the path between city[n[0]]-city[n[1]]?
				de = Distance(city[n[2]], city[n[1]], lccc) + Distance(city[n[3]], city[n[0]], lccc)\
					 - Distance(city[n[2]], city[n[0]], lccc) - Distance(city[n[3]] ,city[n[1]], lccc)

				if de<0 or numpy.exp(-de/T)>numpy.random.rand(): # Metropolis
					accepted += 1
					dist += de
					reverse(city, n)
			else:
				# Here we transpose a segment
				nc = (n[1]+1+ int(numpy.random.rand()*(nn-1)))%nct  # Another point outside n[0],n[1] segment. See picture in lecture nodes!
				n[4] = nc
				n[5] = (nc+1) % nct

				# Cost to transpose a segment
				de = -Distance( city[n[1]], city[n[3]], lccc) - Distance( city[n[0]], city[n[2]], lccc) \
						- Distance( city[n[4]], city[n[5]], lccc)
				de += Distance( city[n[0]], city[n[4]], lccc) + Distance( city[n[1]], city[n[5]], lccc) \
						+ Distance( city[n[2]], city[n[3]], lccc)

				if de<0 or numpy.exp(-de/T)>numpy.random.rand(): # Metropolis
					accepted += 1
					dist += de
					city = transpt(city, n)

			if accepted > maxAccepted: break

		# Plot
		#         Plot(city, R, dist)

		print("T=%10.5f , distance= %10.5f , accepted steps= %d" %(T, dist, accepted))
		T *= fCool             # The system is cooled down
		if accepted == 0: break  # If the path does not want to change any more, we can stop


#     Plot(city, R, dist)
	return city




def pca(cov):
	pass#IMPORTIMPORTIMPORT from numpy import  linalg, argsort
	""" assume one sample per column """
	values, vecs = numpy.linalg.eigh(cov)
	perm = numpy.argsort(-values)  # sort in descending order
	return values[perm], vecs[:, perm]


