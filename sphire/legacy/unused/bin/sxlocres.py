"""1
		res_overall = 0.5
		if myid ==main_node:
			fsc_curve = fsc(vi, ui)
			for ifreq in xrange(len(fsc_curve[0])-1, -1, -1):
				if fsc_curve[1][ifreq] > options.cutoff:
					res_overall = fsc_curve[0][ifreq]
					break
		res_overall = bcast_number_to_all(res_overall, main_node)
		"""
"""2
		res_overall = 0.5
		fsc_curve = fsc(vi, ui)
		for ifreq in xrange(len(fsc_curve[0])-1, -1, -1):
			if fsc_curve[1][ifreq] > options.cutoff:
				res_overall = fsc_curve[0][ifreq]
				break
		"""		
