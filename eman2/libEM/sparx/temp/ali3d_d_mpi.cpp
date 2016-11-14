#include "stdlib.h"
#include "mpi.h"

// include EMAN2
#include "emdata.h"
#include "emassert.h"
#include "projector.h"
#include "transform.h"
#include "reconstructor.h"

#include "ali3d_d_mpi.h"
#include "utilcomm.h"
#include "alignoptions.h"
#include "sirt.h"
#include "project3d.h"
#include "fgcalc.h"

using namespace EMAN;

int ali3d_d( MPI_Comm comm, EMData*& volume, EMData** projdata, 
             EMData** cleandata, float *angleshift, int nloc, 
             AlignOptions& options, char* fname_base)
{
    int mypid, ncpus;
    int ierr;
    MPI_Status mpistatus;
    double timer, total_timer;
	
    MPI_Comm_rank(comm,&mypid);
    MPI_Comm_size(comm,&ncpus);

    int * psize, * nbase, nangloc, nang;
    psize = new int[ncpus];
    nbase = new int[ncpus];
    MPI_Allreduce(&nloc, &nang, 1, MPI_INT, MPI_SUM, comm);
    nangloc = setpart(comm, nang, psize, nbase);

    char out_fname[64];
    std::ofstream fsc_out;

    int ndim, nx, ny, nz;

    ndim = volume->get_ndim();

    nx = volume->get_xsize();
    ny = volume->get_ysize();
    nz = volume->get_zsize();

    int first_ring = options.get_first_ring();
    int last_ring  = options.get_last_ring();
    int rstep      = options.get_rstep();
    int ri         = options.get_ri();
    float xrng     = options.get_xrng();
    float yrng     = options.get_yrng();
    float step     = options.get_step();
    float dtheta   = options.get_dtheta();
    bool CTF       = options.get_CTF();
    bool have_angles = options.get_have_angles();
    bool use_sirt  = options.get_use_sirt();
    std::string symmetry = options.get_symmetry();
    int max_iter       = options.get_maxit();

    Dict maskdict;
    maskdict["radius"] = last_ring;
    maskdict["fill"]   = 1;
    
    EMData * mask3D = options.get_mask3D();
    if (mask3D == NULL) {
	mask3D = new EMData();
	
	mask3D->set_size(nx,ny,nz);
	mask3D->process_inplace("testimage.circlesphere", maskdict);
    }
			
    EMData * mask2D = new EMData();
    mask2D->set_size(nx,ny,1);
    mask2D->process_inplace("testimage.circlesphere", maskdict);
    if (first_ring > 0) {
	// subtract off the smaller ring
	EMData * inner_mask2D = new EMData();
	inner_mask2D->set_size(nx,ny,1);
	maskdict["radius"] = first_ring;
	inner_mask2D->process_inplace("testimage.circlesphere", maskdict);
	*mask2D -= *inner_mask2D;
	EMDeletePtr(inner_mask2D);
    }	
    
    EMData ** sirt_images;
    float sirt_tol, sirt_lam;
    int sirt_maxit;
    Dict recons_params;
    std::string recons_name;
    float * even_odd_angleshift = new float[5*nloc];


    if ( use_sirt ) {
	sirt_images = new EMData*[nloc];
	sirt_tol = options.get_sirt_tol();
	sirt_lam = options.get_sirt_lam();
	sirt_maxit = options.get_sirt_maxit();
    } else {
	recons_name = "nn4";
	recons_params["symmetry"] = symmetry;
	recons_params["size"]     = nx;
	recons_params["npad"]     = 4; // will this ever be != 4?

	Dict CTF_params;
	CTF_params["sign"] = -1.0;
	int CTF_applied;
	int padffted;
	if (CTF) {
	    // filter each image according to its ctf header info
	    for ( int i = 0 ; i < nloc ; ++i ) {
		CTF_applied = projdata[i]->get_attr("ctf_applied");
		if ( CTF_applied == 0 ) {
		    CTF_params["defocus"]      = projdata[i]->get_attr("defocus");
		    CTF_params["Cs"]           = projdata[i]->get_attr("Cs");
		    CTF_params["Pixel_size"]   = projdata[i]->get_attr("Pixel_size");
		    CTF_params["B_factor"]     = projdata[i]->get_attr("B_factor");
		    CTF_params["amp_contrast"] = projdata[i]->get_attr("amp_contrast");
		    projdata[i]->process_inplace("filter.CTF_",CTF_params);
		}
		// I think the CTF filter should be applied before the images have the background subtracted,
		// but I'm not really sure.
		CTF_applied = cleandata[i]->get_attr("ctf_applied");
		if ( CTF_applied == 0 ) {
		    CTF_params["defocus"]      = cleandata[i]->get_attr("defocus");
		    CTF_params["Cs"]           = cleandata[i]->get_attr("Cs");
		    CTF_params["Pixel_size"]   = cleandata[i]->get_attr("Pixel_size");
		    CTF_params["B_factor"]     = cleandata[i]->get_attr("B_factor");
		    CTF_params["amp_contrast"] = cleandata[i]->get_attr("amp_contrast");
		    cleandata[i]->process_inplace("filter.CTF_",CTF_params);
		}
	    }
	    // set a few more things for the reconstructor to deal with
	    recons_name = "nn4_ctf";
	    recons_params["snr"] = options.get_snr();
	    recons_params["sign"] = 1.0; // confirm this, sign is not defined in the current python code.
	    try {
		padffted = projdata[0]->get_attr("padffted");
	    } catch ( std::exception& e ) { // the only exception thrown by get_attr() is NotExistingObjectException()
		padffted = 0;
	    }
	    if (padffted == 1) recons_params["size"] = (float) recons_params["size"] / (float) recons_params["npad"];
	}
    }

    // The below are done in proj_ali_incore
    // Perhaps could be moved inside a major loop iterating over 
    // different dtheta
    // xrng, and xshift settings.
    std::string mode("F");

    // The below ensures that dtheta is sufficiently small to make the inverse problem over-determined
    int num_ref = 0;
    std::vector<float> ref_angles;
    int max_dimension = (nx > ny ? nx : ny);
    max_dimension = (max_dimension > nz ? max_dimension : nz);
    while (true) {
	ref_angles = Util::even_angles(dtheta);
	num_ref = ref_angles.size() / 3;
	if (num_ref >= max_dimension) break;
	dtheta -= 1.0;
	if (mypid == 0) printf("Problem is underdetermined: decreasing dtheta to %f\n",dtheta);
    }
    // should this ever not be true?  Is it really needed?
    bool phi_eq_minus_psi = true;
    if (phi_eq_minus_psi) {
	for ( int i = 0 ; i < num_ref ; ++ i ) {
	    ref_angles[3*i + 2] = fmod(720.0 - ref_angles[3*i], 360.0);
	}
    }
    
    int numref; // used after call to multiref_polar_ali_2d, stores the index of the image with max correlation

    Util::mul_img(volume, mask3D);

    std::vector<int> numr = Numrinit(first_ring, last_ring, rstep, mode);
    std::vector<float> wr = ringwe(numr, mode);

    std::vector<EMData *> ref_proj_rings;
    EMData * ref_proj_ptr;
    EMData * cimage_ptr;
    Dict volparams;
    volparams["angletype"] = "SPIDER";
    volparams["radius"] = ri;
    Dict proj_params;
    proj_params["ctf_applied"] = 0;
    std::vector<float> anglelist(3,0.0);

    float sxo, syo;
    std::vector<float> best_alignment(6,0.0);

    Transform3D::EulerType EULER_SPIDER = Transform3D::SPIDER;
    EMUtil::ImageType WRITE_SPI = EMUtil::IMAGE_SINGLE_SPIDER;
    float compphi, comptheta, comppsi;
    float angb, sxb, syb, ct;

    float cnx = floor(nx / 2.0) + 1.0;
    float cny = floor(ny / 2.0) + 1.0;

    char iter_string[4];

    Dict make_ref_proj_dict;
    make_ref_proj_dict["mask"] = mask2D;
    make_ref_proj_dict["no_sigma"] = 1;

    float multiref_res, unified_res;
    float * volsph, * voldata;
    Vec3i volsize;
    Vec3i origin;
    volsize[0] = nx;
    volsize[1] = ny;
    volsize[2] = nz;
    origin[0] = nx/2 + 1;
    origin[1] = ny/2 + 1;
    origin[2] = nz/2 + 1;

    int nnz, nrays;

    ierr = getnnz(volsize, ri, origin, &nrays, &nnz);
    int * ptrs = new int[nrays+1];
    int * cord = new int[3*nrays];
    volsph = new float[nnz];
    voldata = volume->get_data();
    ierr = cb2sph(voldata, volsize, ri, origin, nnz, ptrs, cord, volsph); 
    
    float angtrs[5];
    float * img_data;

    EMData * vol1 = NULL;
    EMData * vol2 = NULL;

    EMData * fftvol;
    EMData * weight;
    Reconstructor* r;
    float * fftvol_send;
    float * weight_send;
    int fftvol_size;
    int weight_size;
    float * fftvol_recv;
    float * weight_recv;

    std::vector<float> stats;
    std::vector<float> fsc_result;
    int fsc_size;

    float ring_width = 1.0;
    float phi, theta, psi;
    int filter_cutoff_index;

    // these are settings to compare with the fourier shell correlation 
    float filter_coeff_high = 0.9; 
    float filter_freq_low = 0.25;
    float hcut_default = 0.49;

    // these get passed to the filter
    float lcut, hcut; 

    Dict btwl_dict;

    std::vector<float> ph_cog;
    Dict cog_dict;

    total_timer = MPI_Wtime();

    float * peak_corr = new float[nloc];
    std::ofstream corr_out;

    // These keep track of how many times we take new parameters based on projection matching
    int new_angles, old_angles, new_total, old_total;
    
    for ( int i = 0 ; i < max_iter ; ++i ) { // This loop is as in ali3d_d
	new_angles = 0;
	old_angles = 0;
	// entering proj_ali_incore
	// generate reference projections and push pointers to them on the back of ref_proj_rings
	if (mypid == 0) std::cout << "Generating " << num_ref << " reference projections..." << std::endl;
	for ( int j = 0 ; j < num_ref ; ++j ) { 
	    for ( int k = 0 ; k < 3 ; ++k) anglelist[k] = ref_angles[3 * j + k];
	    volparams["anglelist"] = anglelist;
	    ref_proj_ptr = volume->project("chao", volparams);// This uses a version of fwdpj3() internally
	    proj_params["phi"] = anglelist[0];
	    proj_params["theta"] = anglelist[1];
	    proj_params["psi"] = anglelist[2];
	    ref_proj_ptr->set_attr_dict(proj_params);
	    ref_proj_ptr->process_inplace("normalize.mask", make_ref_proj_dict);

	    cimage_ptr = Util::Polar2Dm(ref_proj_ptr, cnx, cny, numr, mode);
	    Util::Frngs(cimage_ptr, numr);
	    Applyws(cimage_ptr, numr, wr);
	    ref_proj_rings.push_back(cimage_ptr);

	}
	// match the data images to the reference projections
	if (mypid == 0) std::cout << "Matching data to reference projections..." << std::endl;
	timer = MPI_Wtime();
	for ( int j = 0 ; j < nloc ; ++j ) {
	    if (have_angles) {
		sxo = angleshift[5*j + 3];
		syo = angleshift[5*j + 4];
	    } else {
		sxo = projdata[j]->get_attr("s2x");
		syo = projdata[j]->get_attr("s2y");
	    }
	    best_alignment = Util::multiref_polar_ali_2d(projdata[j], ref_proj_rings, xrng, yrng, step, mode, numr, cnx - sxo, cny - syo);
	    numref = (int) best_alignment[4];

	    peak_corr[j] = best_alignment[5]; // best_alignment[5] = peak
	    // Here we implement what goes on in compose_transform2(), from utilities.py
	    Transform3D R1(EULER_SPIDER , 0.0, 0.0, 0.0);
	    R1.set_posttrans(Vec3f(best_alignment[1],best_alignment[2],0.0)); // best_alignment[1], best_alignment[2] = sxs, sys
	    R1.set_scale(1.0);

	    Transform3D R2(EULER_SPIDER, 0.0, 0.0, -best_alignment[0]); // best_alignment[0] = ang
	    R2.set_posttrans(Vec3f(0.0,0.0,0.0));
	    R2.set_scale(1.0);

	    Transform3D Rcomp = R2*R1;
	    Dict compeuler = Rcomp.get_rotation(EULER_SPIDER);

	    compphi   = fmod((float) compeuler["phi"]   + 360.0, 360.0);
	    comptheta = fmod((float) compeuler["theta"] + 360.0, 360.0);
	    comppsi   = fmod((float) compeuler["psi"]   + 360.0, 360.0);

	    angb = fmod((float) compphi + comppsi,  (float) 360.0);  
	    sxb  = Rcomp.at(0,3);
	    syb  = Rcomp.at(1,3);
	    ct   = Rcomp.get_scale();
	    // end call to compose_transform2()
	    Dict composition_dict;
	    // perhaps should do all of these with set_attr(), rather than set_attr_dict()
	    if (best_alignment[3] != 0.0) { // best_alignment[3] = mirror
		composition_dict["phi"]   = fmod(ref_angles[3 * numref] + 540.0, 360.0);
		composition_dict["theta"] = 180.0 - ref_angles[3 * numref + 1];
		composition_dict["psi"]   = fmod(540.0 - ref_angles[3 * numref + 2] + angb, 360.0);
		composition_dict["s2x"]   = sxb + sxo;
		composition_dict["s2y"]   = syb + syo;
	    } else {
		composition_dict["phi"]   = ref_angles[3 * numref];
		composition_dict["theta"] = ref_angles[3 * numref + 1];
		composition_dict["psi"]   = fmod( ref_angles[3 * numref + 2] + angb + 360.0, 360.0);
		composition_dict["s2x"]   = sxb + sxo;
		composition_dict["s2y"]   = syb + syo;
	    }
	    if (have_angles) {
		// compare residual for this image under old angles (from refinement) and new angles (from projection matching)
		img_data = projdata[j]->get_data();

		angtrs[0] = (float) composition_dict["phi"]   * (PI/180.0);
		angtrs[1] = (float) composition_dict["theta"] * (PI/180.0);
		angtrs[2] = (float) composition_dict["psi"]   * (PI/180.0);
		angtrs[3] = (float) composition_dict["s2x"]   * -1.0;
		angtrs[4] = (float) composition_dict["s2y"]   * -1.0;
		
		ierr = fcalc(volsph, volsize, 
			     nnz, nrays, origin, ri, 
			     ptrs, cord, angtrs, 1, 
			     img_data, &multiref_res);

		// set angtrs to the values in angleshift
		angtrs[0] = angleshift[5*j + 0] * (PI/180.0);
		angtrs[1] = angleshift[5*j + 1] * (PI/180.0);
		angtrs[2] = angleshift[5*j + 2] * (PI/180.0);
		angtrs[3] = angleshift[5*j + 3] * -1.0;
		angtrs[4] = angleshift[5*j + 4] * -1.0;

		ierr = fcalc(volsph, volsize, 
			     nnz, nrays, origin, ri, 
			     ptrs, cord, angtrs, 1, 
			     img_data, &unified_res);

		// whichever gives better residual gets used
		if ( multiref_res < unified_res ) {
		    ++new_angles;
		    angleshift[5*j + 0] = composition_dict["phi"];
		    angleshift[5*j + 1] = composition_dict["theta"];
		    angleshift[5*j + 2] = composition_dict["psi"];
		    angleshift[5*j + 3] = composition_dict["s2x"];
		    angleshift[5*j + 4] = composition_dict["s2y"];
		} else { // multiref_res >= unified_res
		    ++old_angles;
		    composition_dict["phi"]   = angleshift[5*j + 0];
		    composition_dict["theta"] = angleshift[5*j + 1];
		    composition_dict["psi"]   = angleshift[5*j + 2];
		    composition_dict["s2x"]   = angleshift[5*j + 3];
		    composition_dict["s2y"]   = angleshift[5*j + 4];
		}
	    } else { // don't have intial angles, and need to set them
		angleshift[5*j + 0] = composition_dict["phi"];
		angleshift[5*j + 1] = composition_dict["theta"];
		angleshift[5*j + 2] = composition_dict["psi"];
		angleshift[5*j + 3] = composition_dict["s2x"];
		angleshift[5*j + 4] = composition_dict["s2y"];
	    }
	    // set the header of projdata[j] and cleandata[j] to the best angle and shift values
	    projdata[j]->set_attr_dict(composition_dict);
	    cleandata[j]->set_attr_dict(composition_dict);
	}

	MPI_Reduce(&new_angles, &new_total, 1, MPI_INT, MPI_SUM, 0, comm);
	MPI_Reduce(&old_angles, &old_total, 1, MPI_INT, MPI_SUM, 0, comm);
	if (mypid == 0) {
	  printf("   Wall clock seconds for alignment, iteration %d = %11.3e\n",
		 i, MPI_Wtime() - timer);
	  if (have_angles)
	      printf("   \nUsing %d sets of parameters from projection matching, %d sets from previous unified refinement\n\n", new_total, old_total);
	}

	sprintf(out_fname, "corr_%s_pm%d.dat", fname_base, i);
	corr_out.open(out_fname);
	if (mypid == 0) {
	    for ( int corr_index = 0 ; corr_index < nloc ; ++corr_index ) {
		corr_out << std::scientific << peak_corr[corr_index] << std::endl;
	    }
	    for ( int k = 1 ; k < ncpus ; ++k ) {
		ierr = MPI_Recv(peak_corr, psize[k], MPI_FLOAT, k, k, comm, &mpistatus);
		for ( int corr_index = 0 ; corr_index < psize[k] ; ++corr_index ) {
		    corr_out << std::scientific << peak_corr[corr_index] << std::endl;
		}
	    }
	} else { // mypid != 0 , send my data to master to write out
	    ierr = MPI_Send(peak_corr, nloc, MPI_FLOAT, 0, mypid, comm);
	}
	corr_out.close();

	timer = MPI_Wtime();
	// leaving proj_ali_incore

	// build two reconstructions
	if (mypid == 0) 
           printf("   Building even/odd volumes for FSC calculation...\n");
	if ( use_sirt ) {

            // construct a volume using even numbered images
	    vol1 = new EMData();
	    // lam should be based on number of images, right?
	    for ( int j = 0 ; j < nloc ; j += 2 ) {
		sirt_images[j/2] = cleandata[j];
		even_odd_angleshift[5*(j/2) + 0] = angleshift[5*j + 0];
		even_odd_angleshift[5*(j/2) + 1] = angleshift[5*j + 1];
		even_odd_angleshift[5*(j/2) + 2] = angleshift[5*j + 2];
		even_odd_angleshift[5*(j/2) + 3] = angleshift[5*j + 3];
		even_odd_angleshift[5*(j/2) + 4] = angleshift[5*j + 4];
	    }
	    recons3d_sirt_mpi(comm, sirt_images, even_odd_angleshift, 
                              vol1, (nloc+1)/2, ri, sirt_lam, sirt_maxit, 
                              symmetry, sirt_tol);

            // construct a volume using odd numbered images
	    vol2 = new EMData();
	    for ( int j = 1 ; j < nloc ; j += 2 ) {
		sirt_images[j/2] = cleandata[j];
		even_odd_angleshift[5*(j/2) + 0] = angleshift[5*j + 0];
		even_odd_angleshift[5*(j/2) + 1] = angleshift[5*j + 1];
		even_odd_angleshift[5*(j/2) + 2] = angleshift[5*j + 2];
		even_odd_angleshift[5*(j/2) + 3] = angleshift[5*j + 3];
		even_odd_angleshift[5*(j/2) + 4] = angleshift[5*j + 4];
	    }
	    recons3d_sirt_mpi(comm, sirt_images, even_odd_angleshift, 
            vol2, nloc/2, ri, sirt_lam, sirt_maxit, symmetry, sirt_tol);
	} else {

            // construct a volume using even numbered images
	    fftvol = new EMData();
	    weight = new EMData();
	    recons_params["fftvol"] = fftvol;
	    recons_params["weight"] = weight;
	    r = Factory<Reconstructor>::get(recons_name, recons_params);
	    r->setup();
	    for ( int j = 0 ; j < nloc ; j += 2 ) { 
		phi   = projdata[j]->get_attr("phi");
		theta = projdata[j]->get_attr("theta");
		psi   = projdata[j]->get_attr("psi");
		r->insert_slice(cleandata[j],
                                Transform3D(EULER_SPIDER,phi,theta,psi));
	    }
	    fftvol_size = fftvol->get_xsize() 
                        * fftvol->get_ysize() * fftvol->get_zsize();
	    weight_size = weight->get_xsize() 
                        * weight->get_ysize() * weight->get_zsize();

	    fftvol_send = fftvol->get_data();
	    weight_send = weight->get_data();
	    fftvol_recv = new float[fftvol_size];
	    weight_recv = new float[weight_size];
	    MPI_Allreduce(fftvol_send, fftvol_recv, fftvol_size, 
                          MPI_FLOAT, MPI_SUM, comm);
	    MPI_Allreduce(weight_send, weight_recv, weight_size, 
                          MPI_FLOAT, MPI_SUM, comm);
	    for ( int j = 0 ; j < fftvol_size ; j += 1 ) 
                fftvol_send[j] = fftvol_recv[j];
	    for ( int j = 0 ; j < weight_size ; j += 1 ) 
                weight_send[j] = weight_recv[j];
	    EMDeleteArray(fftvol_recv);
	    EMDeleteArray(weight_recv);
	    vol1 = r->finish();
	    EMDeletePtr(fftvol);
	    EMDeletePtr(weight);

            // construct a volume using odd numbered images
	    fftvol = new EMData();
	    weight = new EMData();
	    recons_params["fftvol"] = fftvol;
	    recons_params["weight"] = weight;
	    r = Factory<Reconstructor>::get(recons_name, recons_params);
	    r->setup();
	    for ( int j = 1 ; j < nloc ; j += 2 ) { 
		phi   = projdata[j]->get_attr("phi");
		theta = projdata[j]->get_attr("theta");
		psi   = projdata[j]->get_attr("psi");
		r->insert_slice(cleandata[j], 
                                Transform3D(EULER_SPIDER,phi,theta,psi));
	    }

	    fftvol_size = fftvol->get_xsize() 
                        * fftvol->get_ysize() * fftvol->get_zsize();
	    weight_size = weight->get_xsize() 
                        * weight->get_ysize() * weight->get_zsize();

	    fftvol_send = fftvol->get_data();
	    weight_send = weight->get_data();
	    fftvol_recv = new float[fftvol_size];
	    weight_recv = new float[weight_size];
	    MPI_Allreduce(fftvol_send, fftvol_recv, fftvol_size, 
                          MPI_FLOAT, MPI_SUM, comm);
	    MPI_Allreduce(weight_send, weight_recv, weight_size, 
                          MPI_FLOAT, MPI_SUM, comm);

	    for ( int j = 0 ; j < fftvol_size ; j += 1 ) 
               fftvol_send[j] = fftvol_recv[j];
	    for ( int j = 0 ; j < weight_size ; j += 1 ) 
               weight_send[j] = weight_recv[j];

	    EMDeleteArray(fftvol_recv);
	    EMDeleteArray(weight_recv);
	    vol2 = r->finish();
	    EMDeletePtr(fftvol);
	    EMDeletePtr(weight);
	}

	// calculate and subtract the mean from vol1 and vol2, 
        // and apply the 3D mask
        //
	stats = Util::infomask(vol1, mask3D, false);
	*vol1 -= stats[0]; // stats[0] = mean
	Util::mul_img(vol1, mask3D);
	stats = Util::infomask(vol2, mask3D, false);
	*vol2 -= stats[0]; // stats[0] = mean
	Util::mul_img(vol2, mask3D);

	// calculate fsc
	fsc_result = vol1->calc_fourier_shell_correlation(vol2, ring_width);
	fsc_size = fsc_result.size() / 3;
	sprintf(out_fname, "fsc_%s_pm%d.dat", fname_base, i);
	if (mypid == 0) {
	    fsc_out.open(out_fname);
	    for ( int j = 0 ; j < fsc_size ; ++j ) {
		// Note the indexing of fsc: the frequencies are in the 
                // first fsc_size entries, then the correlation coeffs 
                // after that
		fsc_out << std::scientific 
                        << fsc_result[j] << '\t'
                        << fsc_result[j + fsc_size] << '\t'
                        << fsc_result[j + 2 * fsc_size] 
                        << std::endl;
	    }
	    fsc_out.close();
	}
	EMDeletePtr(vol1);
	EMDeletePtr(vol2);

	// and reconstruct the volume from all images
	if (mypid == 0) 
           std::cout << "Building 3-D reconstruction ... " << std::endl;
	if ( use_sirt ) {
	    volume = new EMData();
	    recons3d_sirt_mpi(comm, cleandata, angleshift, volume, nloc, 
                              ri, sirt_lam, sirt_maxit, symmetry, sirt_tol);
	} else {
	    fftvol = new EMData();
	    weight = new EMData();
	    recons_params["fftvol"] = fftvol;
	    recons_params["weight"] = weight;
	    r = Factory<Reconstructor>::get(recons_name, recons_params);
	    r->setup();
	    for ( int j = 0 ; j < nloc ; ++j ) { 
		phi   = projdata[j]->get_attr("phi");
		theta = projdata[j]->get_attr("theta");
		psi   = projdata[j]->get_attr("psi");
		r->insert_slice(cleandata[j], 
                                Transform3D(EULER_SPIDER,phi,theta,psi));
	    }

	    fftvol_size = fftvol->get_xsize() 
                        * fftvol->get_ysize() * fftvol->get_zsize();
	    weight_size = weight->get_xsize()
                        * weight->get_ysize() * weight->get_zsize();

	    fftvol_send = fftvol->get_data();
	    weight_send = weight->get_data();

	    fftvol_recv = new float[fftvol_size];
	    weight_recv = new float[weight_size];

	    MPI_Allreduce(fftvol_send, fftvol_recv, fftvol_size, 
                          MPI_FLOAT, MPI_SUM, comm);
	    MPI_Allreduce(weight_send, weight_recv, weight_size, 
                          MPI_FLOAT, MPI_SUM, comm);

	    for ( int j = 0 ; j < fftvol_size ; ++j ) {
		fftvol_send[j] = fftvol_recv[j];
	    }
	    for ( int j = 0 ; j < weight_size ; ++j ) {
		weight_send[j] = weight_recv[j];
	    }

	    EMDeleteArray(fftvol_recv);
	    EMDeleteArray(weight_recv);
	    EMDeletePtr(volume); 

            // r->finish() returns EMData *, must free the old one
	    volume = r->finish();
	    EMDeletePtr(fftvol);
	    EMDeletePtr(weight);
	}

	sprintf(out_fname, "vol_%s_pm%d.spi", fname_base, i);
	if (mypid == 0) volume->write_image(out_fname, 0, WRITE_SPI);
	// End of reconstruction

	// and filter it based on fourier shell correlation
	filter_cutoff_index = 0;
	while (filter_cutoff_index < fsc_size
	       && fsc_result[filter_cutoff_index + fsc_size] > filter_coeff_high
	       && fsc_result[filter_cutoff_index] < filter_freq_low) {
	    ++filter_cutoff_index;
	}

	lcut  = fsc_result[filter_cutoff_index * 3];
        // hcut = min(f_l + 0.1, f_h_default)
	hcut = (lcut + 0.1 < hcut_default ? lcut + 0.1 : hcut_default); 

	btwl_dict["low_cutoff_frequency"]  = lcut;
	btwl_dict["high_cutoff_frequency"] = hcut;
	volume->set_attr("npad",1);
	volume->process_inplace("filter.lowpass.butterworth", btwl_dict);

	// calculate the center of gravity
	ph_cog = volume->phase_cog();
	cog_dict["x_shift"] = -ph_cog[0];
	cog_dict["y_shift"] = -ph_cog[1];
	cog_dict["z_shift"] = -ph_cog[2];
	// vol = fshift(vol, -cs[0], -cs[1], -cs[2]
	volume->process_inplace("filter.shift", cog_dict);

	// and write to disk
	sprintf(out_fname, "vol_f_%s_pm%d.spi", fname_base, i);
	if (mypid == 0) volume->write_image(out_fname, 0, WRITE_SPI);

	if (mypid == 0) {
	    printf("   Wall clock seconds for reconstructions and fsc, iteration %d = %11.3e\n",
		   i, MPI_Wtime() - timer);
	}

	// clean up ref_proj_rings
	ref_proj_rings.clear();
    }

    if (mypid == 0) {
	std::cout << "Total wall clock seconds for " << max_iter << " iterations on " << ncpus << " processors : ";
	std::cout << MPI_Wtime() - total_timer << std::endl;
    }

    if (options.get_mask3D() == NULL) EMDeletePtr(mask3D); // If it was created by new in this function, otherwise it belongs to someone else.
    EMDeletePtr(mask2D);	
    EMDeleteArray(peak_corr);
    EMDeleteArray(ptrs);
    EMDeleteArray(cord);
    EMDeleteArray(volsph);
    EMDeleteArray(psize);
    EMDeleteArray(nbase);
    if ( use_sirt) EMDeleteArray(sirt_images);
    EMDeleteArray(even_odd_angleshift);

    return 0;
	
}

std::vector<int> Numrinit(int first_ring, int last_ring, int skip, std::string mode)
{
	float dpi;
	int MAXFFT = 32768;

	if (mode == "f" or mode == "F") {
		dpi = 2*PI;
	} else {
		dpi = PI;
	}

	std::vector<int> numr;
	int lcirc = 1;
	int jp, ip, exponent, m;
	for ( int k = first_ring ; k < last_ring + 1 ; k += skip ) {
	    numr.push_back(k);
	    jp = int(dpi * k + 0.5);
	    exponent = -1;
	    m = 1;
	    while (m <= jp){
		m <<= 1;
		++exponent;
	    }
	    // exponent is the largest int such that 2**exponent <= jp                                                                                                           
	    ip = 2 << exponent;
	    if ( (k + skip <= last_ring && jp > ip + ip / 2) || (k + skip > last_ring && jp > ip + ip / 5)) {
		ip = MAXFFT < 2*ip ? MAXFFT : 2*ip;
	    }
	    numr.push_back(lcirc);
	    numr.push_back(ip);
	    lcirc += ip;
	}	
	return numr;
}

std::vector<float> ringwe(std::vector<int> numr, std::string mode)
{
	float dpi;

	if (mode == "f" or mode == "F") {
		dpi = 2*PI;
	} else {
		dpi = PI;
	}
	
	int nring = numr.size() / 3;
	std::vector<float> wr(nring, 0.0);
	float maxrin = (float) numr.back();
	for ( int i = 0 ; i < nring ; ++i ) {
		// No idea what this line is supposed to mean.
		// I hope the order of operations hold...
		wr[i] = numr[i * 3]*dpi/((float) numr[2+i*3]) * maxrin / ((float) numr[2+i*3]);
	}
	return wr;
}

int Applyws(EMData * circ, std::vector<int> numr, std::vector<float> wr)
{
	int nring = numr.size() / 3;
	int maxrin = numr.back();
	int numr3i, numr2i;
	int jc;
	float w;
	for ( int i = 0 ; i < nring ; ++i ) {
		numr3i = numr[2 + i * 3];
		numr2i = numr[1 + i * 3] - 1;
		w = wr[i];
		// why didn't they overload EMData::operator[] ?
		(*circ)(numr2i) *= w;
		if (numr3i == maxrin) {
			(*circ)(numr2i + 1) *= w;
		} else {
			(*circ)(numr2i + 1) *= 0.5 * w;
		}
		for ( int j = 2 ; j < numr3i ; ++j ) {
			jc = j + numr2i;
			(*circ)(jc) *= w;
		}
	}
	return 0;
}

EMData * recons3d_4nn(std::string stack_name, std::vector<int> list_proj, int npad)
// This function should NOT be used, since it reads images from a stack, rather than
// an array of EMData pointers.
{
	Transform3D::EulerType Ttype = Transform3D::SPIDER; // might need :: instead of .
	
	EMData * proj = new EMData();
	proj->read_image(stack_name, list_proj[0]);
	float phi   = proj->get_attr("phi");
	float theta = proj->get_attr("theta");
	float psi   = proj->get_attr("psi");
	int   size  = proj->get_xsize();
	int active  = proj->get_attr("active");
	// 
	Dict params;
	params["symmetry"] = "c1";
	params["size"]     = size;
	params["npad"]     = npad;

	if (size != proj->get_ysize()) {
		std::cerr << "Use a square image!" << std::endl;
		// Throw an exception, the image must be square.
	}
	Reconstructor* r = Factory<Reconstructor>::get("nn4", params);
	r->setup();
	if (active) {
		r->insert_slice(proj, Transform3D(Ttype, phi, theta, psi));
	}
	int nimages = list_proj.size();
	for ( int i = 1 ; i < nimages ; ++i ) {
		proj->read_image(stack_name, list_proj[i]);
		active = proj->get_attr("active");
		if (active) {
			phi   = proj->get_attr("phi");
			theta = proj->get_attr("theta");
			psi   = proj->get_attr("psi");
			r->insert_slice(proj, Transform3D(Ttype,phi,theta,psi));
		}
	}	
	EMDeletePtr(proj);
	return r->finish();
}

