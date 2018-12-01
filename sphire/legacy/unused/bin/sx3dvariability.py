"""1
		if options.SND:
			pass#IMPORTIMPORTIMPORT from projection		import prep_vol, prgs
			pass#IMPORTIMPORTIMPORT from statistics		import im_diff
			pass#IMPORTIMPORTIMPORT from utilities		import get_im, model_circle, get_params_proj, set_params_proj
			pass#IMPORTIMPORTIMPORT from utilities		import get_ctf, generate_ctf
			pass#IMPORTIMPORTIMPORT from filter			import filt_ctf
		
			imgdata = EMData.read_images(stack, range(img_begin, img_end))

			if options.CTF:
				vol = recons3d_4nn_ctf_MPI(myid, imgdata, 1.0, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			else:
				vol = recons3d_4nn_MPI(myid, imgdata, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)

			bcast_EMData_to_all(vol, myid)
			volft, kb = prep_vol(vol)

			mask = model_circle(nx/2-2, nx, ny)
			varList = []
			for i in xrange(img_begin, img_end):
				phi, theta, psi, s2x, s2y = get_params_proj(imgdata[i-img_begin])
				ref_prj = prgs(volft, kb, [phi, theta, psi, -s2x, -s2y])
				if options.CTF:
					ctf_params = get_ctf(imgdata[i-img_begin])
					ref_prj = filt_ctf(ref_prj, generate_ctf(ctf_params))
				diff, A, B = im_diff(ref_prj, imgdata[i-img_begin], mask)
				diff2 = diff*diff
				set_params_proj(diff2, [phi, theta, psi, s2x, s2y])
				varList.append(diff2)
			mpi_barrier(MPI_COMM_WORLD)
		"""
'''2
			imgdata2 = EMData.read_images(stack, range(img_begin, img_end))
			if options.fl > 0.0:
				for k in xrange(len(imgdata2)):
					imgdata2[k] = filt_tanl(imgdata2[k], options.fl, options.aa)
			if options.CTF:
				vol = recons3d_4nn_ctf_MPI(myid, imgdata2, 1.0, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			else:
				vol = recons3d_4nn_MPI(myid, imgdata2, symmetry=options.sym, npad=options.npad, xysize=-1, zsize=-1)
			if myid == main_node:
				vol.write_image("vol_ctf.hdf")
				print_msg("Writing to the disk volume reconstructed from averages as		:  %s\n"%("vol_ctf.hdf"))
			del vol, imgdata2
			mpi_barrier(MPI_COMM_WORLD)
			'''
'''3
				if i < 10 and myid == main_node:
					for k in xrange(10):
						grp_imgdata[k].write_image("grp%03d.hdf"%i, k)
				'''
"""4
				if myid == main_node and i==0:
					for pp in xrange(len(grp_imgdata)):
						grp_imgdata[pp].write_image("pp.hdf", pp)
				"""
"""5
				if myid == main_node and i==0:
					for pp in xrange(len(grp_imgdata)):
						grp_imgdata[pp].write_image("qq.hdf", pp)
				"""
"""6
				if myid == main_node:
					ave.write_image("avgv.hdf",i)
					var.write_image("varv.hdf",i)
				"""
"""7
					if myid == 0 and i == 0:
						for k in xrange(nvec):
							eig[k].write_image("eig.hdf", k)
					"""
"""8
								nm = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								nm = int(nm[0])
								members = mpi_recv(nm, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('members', map(int, members))
								members = mpi_recv(nm, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('pix_err', map(float, members))
								members = mpi_recv(3, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('refprojdir', map(float, members))
								"""
"""9
						members = aveList[im].get_attr('members')
						mpi_send(len(members), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						mpi_send(members, len(members), MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						members = aveList[im].get_attr('pix_err')
						mpi_send(members, len(members), MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						try:
							members = aveList[im].get_attr('refprojdir')
							mpi_send(members, 3, MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						except:
							mpi_send([-999.0,-999.0,-999.0], 3, MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						"""
