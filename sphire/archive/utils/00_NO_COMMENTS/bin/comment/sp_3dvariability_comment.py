


































































































































































































































































































































































































































"""0
		if options.SND:
			from sp_projection		import prep_vol, prgs
			from sp_statistics		import im_diff
			from sp_utilities		import get_im, model_circle, get_params_proj, set_params_proj
			from sp_utilities		import get_ctf, generate_ctf
			from sp_filter			import filt_ctf
		
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








































































































































'''1
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

























































'''2
				if nvec > 0:
					eig = pca(input_stacks=grp_imgdata, subavg="", mask_radius=radiuspca, nvec=nvec, incore=True, shuffle=False, genbuf=True)
					for k in range(nvec):
						set_params_proj(eig[k], [phiM, thetaM, 0.0, 0.0, 0.0])
						eigList[k].append(eig[k])
					"""
					if myid == 0 and i == 0:
						for k in xrange(nvec):
							eig[k].write_image("eig.hdf", k)
					"""
				'''
































"""3
								nm = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								nm = int(nm[0])
								members = mpi_recv(nm, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('members', map(int, members))
								members = mpi_recv(nm, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('pix_err', map(float, members))
								members = mpi_recv(3, MPI_FLOAT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
								ave.set_attr('refprojdir', map(float, members))
								"""







"""4
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



















'''5
			if nvec > 0:
				for k in range(nvec):
					if myid == main_node:log_main.add("Reconstruction eigenvolumes", k)
					cont = True
					ITER = 0
					mask2d = model_circle(radiuspca, nx, nx)
					while cont:
						#print "On node %d, iteration %d"%(myid, ITER)
						eig3D = recons3d_4nn_MPI(myid, eigList[k], symmetry=options.sym, npad=options.npad)
						bcast_EMData_to_all(eig3D, myid, main_node)
						if options.fl > 0.0:
							eig3D = filt_tanl(eig3D, options.fl, options.aa)
						if myid == main_node:
							eig3D.write_image(os.path.join(options.outpout_dir, "eig3d_%03d.hdf"%(k, ITER)))
						Util.mul_img( eig3D, model_circle(radiuspca, nx, nx, nx) )
						eig3Df, kb = prep_vol(eig3D)
						del eig3D
						cont = False
						icont = 0
						for l in range(len(eigList[k])):
							phi, theta, psi, s2x, s2y = get_params_proj(eigList[k][l])
							proj = prgs(eig3Df, kb, [phi, theta, psi, s2x, s2y])
							cl = ccc(proj, eigList[k][l], mask2d)
							if cl < 0.0:
								icont += 1
								cont = True
								eigList[k][l] *= -1.0
						u = int(cont)
						u = mpi_reduce([u], 1, MPI_INT, MPI_MAX, main_node, MPI_COMM_WORLD)
						icont = mpi_reduce([icont], 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)

						if myid == main_node:
							u = int(u[0])
							log_main.add(" Eigenvector: ",k," number changed ",int(icont[0]))
						else: u = 0
						u = bcast_number_to_all(u, main_node)
						cont = bool(u)
						ITER += 1

					del eig3Df, kb
					mpi_barrier(MPI_COMM_WORLD)
				del eigList, mask2d
			'''


























































