#include "mpi.h"
#include "emdata.h"

#include "utilcomm.h"

using namespace EMAN;

int ParseAlignOptions(MPI_Comm comm, AlignOptions& options, char* optionsfname,
                      int nvoxels, EMData*& mask3D)
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
    ierr = MPI_Bcast(&mask3D, 1, MPI_INT, 0, comm); 
    // if it's not NULL, need to allocate and bcast its data
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

