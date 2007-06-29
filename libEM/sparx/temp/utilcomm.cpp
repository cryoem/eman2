#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"

#include "ali3d_unified_mpi.h"
#include "utilcomm.h"

using namespace EMAN;

int ReadVandBcast(MPI_Comm comm, EMData *volume, char *volfname)
{
    int mypid, ierr, ndim, nx, ny, nz, mpierr=0;
    MPI_Comm_rank(comm,&mypid);
    FILE *fp=NULL;

    ierr = 0;   
    if (mypid == 0) {
        // check to see if the file exists
        fp = fopen(volfname,"r");
        if (!fp) {
 	    ierr = 1;
            printf("failed to open %s\n", volfname);
        }
        else {
	    volume->read_image(volfname);
	    ndim = volume->get_ndim();
	    nx   = volume->get_xsize();
	    ny   = volume->get_ysize();
	    nz   = volume->get_zsize();
	    printf("ndim = %d, nx = %d, ny = %d, nz = %d\n", ndim, nx, ny, nz);
        }
    }
    mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, comm);
    if (ierr==0) {
        mpierr = MPI_Bcast(&nx, 1, MPI_INT, 0, comm);
        mpierr = MPI_Bcast(&ny, 1, MPI_INT, 0, comm);
        mpierr = MPI_Bcast(&nz, 1, MPI_INT, 0, comm);
        if (mypid !=0) volume->set_size(nx,ny,nz);
	
        float * voldata = volume->get_data();
        mpierr = MPI_Bcast(voldata, nx*ny*nz, MPI_FLOAT, 0, comm);
        ierr = mpierr;
    }
    return ierr;
}

int ReadStackandDist(MPI_Comm comm, EMData ***images2D, char *stackfname, int *nloc)
{
    int ncpus, mypid, ierr, mpierr=0;
    MPI_Status mpistatus;
    EMUtil *my_util;
    int nima; 
    FILE *fp=NULL;

    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);

    ierr = 0;
    if (mypid == 0) {
        fp = fopen(stackfname,"r");
        if (!fp) {
 	    ierr = 1;
            printf("failed to open %s\n", stackfname);
        }
        else {
            nima = my_util->get_image_count(stackfname);
        }
    }
    mpierr = MPI_Bcast(&ierr, 1, MPI_INT, 0, comm);
    if (ierr == 0) {
	int *psize = new int[ncpus];
	int *nbase = new int[ncpus];
    
        mpierr = MPI_Bcast(&nima, 1, MPI_INT, 0, comm);
    
	*nloc = setpart(comm, nima, psize, nbase);
	*images2D = new EMData*[*nloc]; // NB!: whoever calls ReadStackandDist must delete this!

	EMData *img_ptr;
	int img_index;
	float *imgdata;

	// read the first image to get size
	img_ptr = new EMData();
	img_ptr->read_image(stackfname, 0);
	int nx = img_ptr->get_xsize();
	int ny = img_ptr->get_ysize();

	float s2x, s2y;

	if (mypid == 0) {
	    printf("Master node reading and distributing %d images...\n", nima);
	    for ( int ip = 0 ; ip < ncpus ; ++ip ) {
		for ( int i = 0 ; i < psize[ip] ; ++i ) {
		    img_index = nbase[ip] + i;
		    if (ip != 0) {
			img_ptr->read_image(stackfname, img_index);
			// get a pointer to the image's data
			imgdata = img_ptr->get_data();
			// find the x/y shift values if it has them, otherwise set them to 0.0
			try {
			    s2x = (*images2D)[i]->get_attr("s2x");
			} catch ( std::exception& e ) {
			    s2x = 0.0;
			}
			try {
			    s2y = (*images2D)[i]->get_attr("s2y");
			} catch ( std::exception& e ) {
			    s2y = 0.0;
			}
			// send these to processor ip
			MPI_Send(imgdata, nx*ny, MPI_FLOAT, ip, ip, comm);
			MPI_Send(&s2x, 1, MPI_FLOAT, ip, ip, comm);
			MPI_Send(&s2y, 1, MPI_FLOAT, ip, ip, comm);
		    } else { // ip == 0				    
			(*images2D)[i] = new EMData();
			(*images2D)[i]->read_image(stackfname, img_index);
			try {
			    s2x = (*images2D)[i]->get_attr("s2x");
			} catch ( std::exception& e ) {
			    (*images2D)[i]->set_attr("s2x",0.0);
			}
			try {
			    s2y = (*images2D)[i]->get_attr("s2y");
			} catch ( std::exception& e ) {
			    (*images2D)[i]->set_attr("s2y",0.0);
			}
		    }
		}
		printf("finished reading data for processor %d\n", ip);
	    }
	} else { // mypid != 0 : everyone else receives and reads in their data
	    for ( int i = 0 ; i < psize[mypid] ; ++i ) {
		(*images2D)[i] = new EMData();
		(*images2D)[i]->set_size(nx, ny, 1);
		imgdata = (*images2D)[i]->get_data();
		MPI_Recv(imgdata, nx*ny, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		MPI_Recv(&s2x, 1, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		MPI_Recv(&s2y, 1, MPI_FLOAT, 0, mypid, comm, &mpistatus);
		(*images2D)[i]->set_attr("s2x",s2x);
		(*images2D)[i]->set_attr("s2y",s2y);
	    }
	    printf("received %d images for processor %d\n", *nloc, mypid);
	}
	if (mypid == 0) printf("finished reading and distributing data\n");
        ierr = mpierr; 

	EMDeletePtr(img_ptr);
	EMDeleteArray(psize);
	EMDeleteArray(nbase);
    }
    return ierr;
}

int CleanStack(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin)
{
    int nx = volsize[0];
    int ny = volsize[1];
    	
    float * rhs = new float[nx*ny];
    float * imgdata;
	
    // Calculate average "background" from all pixels strictly outside of radius
    double aba, abaloc; // average background
    aba = 0.0;
    abaloc = 0.0;
    int klp, klploc; // number of pixels in background
    klp = 0;
    klploc = 0;
	
    // calculate avg background in parallel
    for ( int i = 0 ; i < nloc ; ++i ) {
 	imgdata = image_stack[i]->get_data();
	asta2(imgdata, nx, ny, ri, &abaloc, &klploc);
    }
	
    MPI_Allreduce(&abaloc, &aba, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&klploc, &klp, 1, MPI_INT, MPI_SUM, comm);

    aba /= klp;

    // subtract off the average background from pixels weakly inside of radius
    int x_summand, y_summand;
    int r_squared = ri * ri;
    for ( int i = 0 ; i < nloc ; ++i ) {
	imgdata = image_stack[i]->get_data();
	for ( int j = 0 ; j < nx ; ++j) {
	    x_summand = (j - origin[0]) *  (j - origin[0]);
	    for ( int k = 0 ; k < ny ; ++k ) {
		y_summand = (k - origin[1]) *  (k - origin[1]);
		if ( x_summand + y_summand <= r_squared) {
		    imgdata[j*ny + k] -= aba;
		}
	    }
	}
    }
    EMDeleteArray(rhs);
    return 0;
}

int setpart(MPI_Comm comm, int nang, int *psize, int *nbase)
{
   int ncpus, mypid, nangloc, nrem;

   MPI_Comm_size(comm,&ncpus);
   MPI_Comm_rank(comm,&mypid);

   nangloc = nang/ncpus;
   nrem    = nang - ncpus*nangloc;
   if (mypid < nrem) nangloc++;

   for (int i = 0; i < ncpus; i++) {
      psize[i] = nang/ncpus;
      if (i < nrem) psize[i]++;
   }
 
   nbase[0] = 0; 
   for (int i = 1; i < ncpus; i++) {
      nbase[i] = nbase[i-1] + psize[i-1];
   }
   
   return nangloc;
}

int ParseAlignOptions(MPI_Comm comm, AlignOptions& options, char* optionsfname, int nvoxels, EMData*& mask3D)
{
   int ncpus, mypid, ierr;

   MPI_Comm_size(comm,&ncpus);
   MPI_Comm_rank(comm,&mypid);

    std::ifstream option_stream;
    std::string current_option;
    char option_buffer[100];
    int int_option, parserr = 0;
    float float_option;

    float * mask_data;
    if ( mypid == 0 ) { // read data from the options file
	option_stream.open(optionsfname);
	while ( option_stream >> current_option ) {
	    if ( current_option == "maskfile" ) {
		option_stream >> current_option;
		mask3D = new EMData();
		mask3D->read_image(current_option);
		options.set_mask3D(mask3D);
	    } 
	    else if ( current_option == "inner_ring" ) {
		option_stream >> int_option;
		options.set_first_ring(int_option);
		std::cout << "first_ring = " << int_option << std::endl;
	    } 
	    else if ( current_option == "outer_ring" ) {
		option_stream >> int_option;
		options.set_last_ring(int_option);
		std::cout << "last_ring = " << int_option << std::endl;
	    } 
	    else if ( current_option == "rstep" ) {
		option_stream >> int_option;
		options.set_rstep(int_option);
		std::cout << "rstep = " << int_option << std::endl;
	    } 
	    else if ( current_option == "radius" ) { // perhaps this is the same as outer_ring?
		option_stream >> int_option;	
		options.set_ri(int_option);
		std::cout << "radius = " << int_option << std::endl;
	    } 
	    else if ( current_option == "x_range" ) {
		option_stream >> float_option;
		options.set_xrng(float_option);
		std::cout << "x_range = " << float_option << std::endl;
	    } 
	    else if ( current_option == "y_range" ) {
		option_stream >> float_option;
		options.set_yrng(float_option);
		std::cout << "y_range = " << float_option << std::endl;
	    }
	    else if ( current_option == "translation_step" ) {
		option_stream >> float_option;
		options.set_step(float_option);
		std::cout << "step = " << float_option << std::endl;
	    }
	    else if ( current_option == "theta_step" ) {
		option_stream >> float_option;
		options.set_dtheta(float_option);
		std::cout << "theta_step = " << float_option << std::endl;
	    }
	    else if ( current_option == "CTF" ) {
		option_stream >> current_option;
		if ( current_option == "true" ) {
		    options.set_CTF(true);
		    std::cout << "CTF = true" << std::endl;
		} 
		else { // anything else sets it to false
		    options.set_CTF(false);
		    std::cout << "CTF = false" << std::endl;
		}
	    }
	    else if ( current_option == "snr" ) {
		option_stream >> float_option;
		options.set_snr(float_option);
		std::cout << "snr = " << float_option << std::endl;
	    }
	    else if ( current_option == "ref_a" ) {
		option_stream >> current_option;
		if ( current_option == "P" ) {
		    options.set_ref_angle_type("P");
		} 
		else if ( current_option == "S" ){ // should add support for this
		    std::cerr << "Currently only support Penczek-type reference angles..." << std::endl;
		    options.set_ref_angle_type("P");
		}
		else { // Currently default to "P", will eventually default to "S"
		}
	    }
	    else if ( current_option == "symmetry" ) {
		option_stream >> current_option;
		options.set_symmetry(current_option);
		std::cout << "symmetry = " << current_option << std::endl;
	    }
	    else if ( current_option == "use_sirt" ) {
		option_stream >> current_option;
		if ( current_option == "true" ) {
		    options.set_use_sirt(true);
		    std::cout << "use_sirt = true" << std::endl;
		} 
		else { // anything else sets it to false
		    options.set_use_sirt(false);
		    std::cout << "use_sirt = false" << std::endl;
		}
            }
	    else if ( current_option == "sirt_tol" ) {
		option_stream >> float_option;
		options.set_sirt_tol(float_option);
		std::cout << "sirt_tol = " << float_option << std::endl;
	    }
	    else if ( current_option == "sirt_lam" ) {
		option_stream >> float_option;
		options.set_sirt_lam(float_option);
		std::cout << "sirt_lam = " << float_option << std::endl;
	    }
	    else if ( current_option == "sirt_maxit" ) {
		option_stream >> int_option;
		options.set_sirt_maxit(int_option);
		std::cout << "sirt_maxit = " << int_option << std::endl;
	    }
	    else if ( current_option == "maxit" ) {
		option_stream >> int_option;
		options.set_maxit(int_option);
		std::cout << "maxit = " << int_option << std::endl;
	    }
	    else { // catch-all
		std::cerr << "Unsupported option " << current_option << "..." << std::endl;
                parserr = 1;
	    }
	}
	option_stream.close();
        if (parserr) { 
	  printf("The supported options are:\n");
          printf("    maskfile \n");
          printf("    inner_ring\n");
          printf("    outer_ring\n");
          printf("    rstep\n");
          printf("    radius\n");
          printf("    x_range\n");
          printf("    y_range\n");
          printf("    translation_step\n");
          printf("    CTF\n");
          printf("    snr\n");
          printf("    ref_a\n");
          printf("    S\n");
          printf("    symmetry\n");
          printf("    sirt_tol\n");
          printf("    sirt_lam\n");
          printf("    sirt_maxit\n");
          printf("     maxit\n");
        }
    }
    // Then broadcast all the data that was read
    ierr = MPI_Bcast(&mask3D, 1, MPI_INT, 0, comm); // if it's not NULL, need to allocate and bcast its data
    // NOTE: this is sending over the master's address for mask3D ONLY as a test to see if it's NULL or not.
    // DO NOT USE THIS POINTER ON ANY OTHER NODES EXCEPT FOR THIS TEST!
    if ( mask3D != NULL ) {
	if ( mypid != 0 ) {
	    mask3D = new EMData();
	}
	mask_data = mask3D->get_data();
	ierr = MPI_Bcast(mask_data, nvoxels, MPI_FLOAT, 0, comm);
	options.set_mask3D(mask3D);
    }

    int_option = options.get_first_ring();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_first_ring(int_option);

    int_option = options.get_last_ring();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_last_ring(int_option);

    int_option = options.get_rstep();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_rstep(int_option);

    int_option = options.get_ri();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_ri(int_option);

    float_option = options.get_xrng();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_xrng(float_option);

    float_option = options.get_yrng();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_yrng(float_option);

    float_option = options.get_step();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_step(float_option);

    float_option = options.get_dtheta();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_dtheta(float_option);

    int_option = options.get_CTF();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_CTF(int_option);

    float_option = options.get_snr();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_snr(float_option);

    current_option = options.get_ref_angle_type();
    strcpy(option_buffer, current_option.c_str());
    ierr = MPI_Bcast(option_buffer, current_option.size(), MPI_CHAR, 0, comm);
    options.set_ref_angle_type(option_buffer);

    current_option = options.get_symmetry();
    strcpy(option_buffer, current_option.c_str());
    ierr = MPI_Bcast(option_buffer, current_option.size(), MPI_CHAR, 0, comm);
    options.set_symmetry(option_buffer);

    int_option = options.get_use_sirt();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_use_sirt(int_option);

    float_option = options.get_sirt_tol();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_sirt_tol(float_option);

    float_option = options.get_sirt_lam();
    ierr = MPI_Bcast(&float_option, 1, MPI_FLOAT, 0, comm);
    options.set_sirt_lam(float_option);

    int_option = options.get_sirt_maxit();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_sirt_maxit(int_option);

    int_option = options.get_maxit();
    ierr = MPI_Bcast(&int_option, 1, MPI_INT, 0, comm);
    options.set_maxit(int_option);

    return 0; // ParseAlignOptions
}
