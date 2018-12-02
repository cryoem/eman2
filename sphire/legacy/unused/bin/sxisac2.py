'''1
def create_subgroup_within_group(mpi_comm):
	# select a subset of myids to be in subdivision
	if( Blockdata["Ykey"]%Blockdata["no_of_processes_per_group"] == 0): submyids = [Blockdata["Ykey"]]
	else:  submyids = []

	#if( Blockdata["Ykey"]%Blockdata["no_of_processes_per_group"] == 0):  print "  A  ",Blockdata["myid"],submyids

	submyids = wrap_mpi_gatherv(submyids, 0, mpi_comm)
	submyids = wrap_mpi_bcast(submyids, 0, mpi_comm)
	#if( Blockdata["Ykey"] == 0 ): print   "  B  ",submyids
	world_group = mpi_comm_group(mpi_comm)
	subgroup = mpi_group_incl(world_group,len(submyids),submyids)
	#print   " XXX world group  ",Blockdata["myid"],Blockdata["Ykey"],world_group,subgroup
	Blockdata["subgroup_comm"] = mpi_comm_create(mpi_comm, subgroup)
	mpi_barrier(mpi_comm)
	#print   " ZZZ subgroup  ",Blockdata["myid"],Blockdata["Ykey"],world_group,subgroup,Blockdata["subgroup_comm"]

	Blockdata["subgroup_size"] = -1
	Blockdata["subgroup_myid"] = -1
	if (MPI_COMM_NULL != Blockdata["subgroup_comm"]):
		#print "  YYY  ",Blockdata["myid"],Blockdata["Ykey"]
		Blockdata["subgroup_size"] = mpi_comm_size(Blockdata["subgroup_comm"])
		Blockdata["subgroup_myid"] = mpi_comm_rank(Blockdata["subgroup_comm"])
	#else:  print "  UUU  ",Blockdata["myid"],Blockdata["Ykey"]
	#  "nodes" are zero COUs on subgroups that do indep_run of independent runs.
	#print "  FFF  ",Blockdata["myid"],Blockdata["Ykey"],Blockdata["subgroup_myid"],Blockdata["subgroup_size"]
	mpi_barrier(mpi_comm)
	return

'''
"""2
	if key == group_main_node:
		all_ali_params = [None]*len(data)
		for i,im in enumerate(data):
			alpha, sx, sy, mirror, scale = get_params2D(im)
			all_ali_params[i] = [alpha, sx, sy, mirror]
		write_text_row(all_ali_params, "params.txt")
		del all_ali_params

		for i in xrange(K):
			#  Each color has the same set of refim
			refi[i].write_image("class_averages.hdf", i)
	"""
