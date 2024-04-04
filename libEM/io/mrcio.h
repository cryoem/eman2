/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#ifndef eman__mrcio_h__
#define eman__mrcio_h__ 1

#include "imageio.h"
#include "renderer.h"

namespace EMAN
{
	/** MRC file = header + data (nx x ny x nz).
	 * A MRC image file stores 1D, 2D or 3D image. The image's
	 * dimensions and pixel type are defined in the header.
	 */
	
	class MrcIO : public ImageIO, public Renderer
	{
	public:
		explicit MrcIO(const string & fname, IOMode rw_mode = READ_ONLY);
		~MrcIO();

		DEFINE_IMAGEIO_FUNC;

		int read_ctf(Ctf & ctf, int image_index = 0);
		void write_ctf(const Ctf & ctf, int image_index = 0);

		static bool is_valid(const void *first_block, off_t file_size = 0);
		static int get_mode_size(int mm);
		static int to_em_datatype(int mrcmode);
		static int to_mrcmode(int em_datatype, int is_complex);

		int get_nimg();

	private:
		enum MrcMode {
			MRC_UCHAR = 0,
			MRC_SHORT = 1,
			MRC_FLOAT = 2,
			MRC_SHORT_COMPLEX=3,
			MRC_FLOAT_COMPLEX=4,
			MRC_USHORT = 6,
			MRC_UCHAR3 = 16,	// unsigned char * 3, for rgb data, non-standard
			MRC_CHAR = 17,		// non-standard - signed char
			MRC_UHEX = 101,	// 2 4-bit values packed into each 8-bit byte
			MRC_UNKNOWN = 18
		};

		enum {
			MRC_NUM_LABELS = 10,
			MRC_LABEL_SIZE = 80,
			NUM_4BYTES_PRE_MAP = 52,
			NUM_4BYTES_AFTER_MAP = 3
		};

		/* updated to MRC Image2000 format which is compatible with CCP4 format */
		/* https://www.ccpem.ac.uk/mrc_format/mrc2014.php */
		struct MrcHeader
		{
			int nx;				/* 0 - number of columns */
			int ny;				/* 1 - number of rows */
			int nz;				/* 2 - number of sections */

			int mode;			/* 3 - See modes above. */

			int nxstart;		/* 4 - No. of first column in map, default 0. */
			int nystart;		/* 5 - No. of first row in map, default 0. */
			int nzstart;		/* 6 - No. of first section in map,default 0. */

			int mx;				/* 7 - Number of intervals along X. */
			int my;				/* 8 - Number of intervals along Y. */
			int mz;				/* 9 - Number of intervals along Z. */

			/* Cell: treat a whole 2D image as a cell */
			float xlen;			/* 10 - Cell dimensions (Angstroms). */
			float ylen;			/* 11 - Cell dimensions (Angstroms). */
			float zlen;			/* 12 - Cell dimensions (Angstroms). */

			float alpha;		/* 13 - Cell angles (Degrees). */
			float beta;			/* 14 - Cell angles (Degrees). */
			float gamma;		/* 15 - Cell angles (Degrees). */

			/* axis X => 1, Y => 2, Z => 3 */
			int mapc;			/* 16 - Which axis corresponds to Columns.  */
			int mapr;			/* 17 - Which axis corresponds to Rows.     */
			int maps;			/* 18 - Which axis corresponds to Sections. */

			float amin;			/* 19 - Minimum density value. */
			float amax;			/* 20 - Maximum density value. */
			float amean;		/* 21 - Mean density value.    */

			int ispg;			/* 22 - Space group number (0 for images). */

			int nsymbt;			/* 23 - Number of bytes in extended header. */

			/* https://bio3d.colorado.edu/imod/doc/libhelp/mrcfiles.html */
			short creatid;
			char blank[6];

			/* https://www.ccpem.ac.uk/mrc_format/mrc2014.php */
			char exttyp[4];
			int nversion;
			char blank2[16];

			/* https://bio3d.colorado.edu/imod/doc/libhelp/mrcfiles.html */
			short nint;
			short nreal;

			short sub;
			short zfac;
			float min2;
			float max2;
			float min3;
			float max3;
			char imod_stamp[4]; // 1146047817 if int or 'IMOD' if char[4] indicates that file was created by IMOD or
			                    // other software that uses bit flags in the following field
			int imod_flags;			/* bit flags used by IMOD - >= 16 for 8 bit packed */
			short idtype;
			short lens;
			short nd1;
			short nd2;
			short vd1;
			short vd2;
			float tiltangles[6];

			float xorigin;		/* X origin. */
			float yorigin;		/* Y origin. */
			float zorigin;		/* Z origin. */

			char map[4];		/* constant string "MAP "  */
			int machinestamp;	/* machine stamp in CCP4 convention:
								   big endian=0x11110000 little endian=0x44440000 */
						/* There is an ambiguity in the specification, using 0x11111111 & 4 instead */

			float rms;			/* rms deviation of map from mean density */

			int nlabels;		/* Number of labels being used. */
			char labels[MRC_NUM_LABELS][MRC_LABEL_SIZE];
		};

		/** The extended header used by Fei MRC image. It contains the information about a maximum of 1024 images.
		 * Each section is 128 bytes long. The extended header is thus 1024*128 bytes (always the same length,
		 * reagrdless of how many images are present.)
		 * Not always 1024*128 bytes, but at least 128*nz. the length of extended header is defined in the regular header field "next" */
		/* https://www.ccpem.ac.uk/downloads/EPU_MRC2014_File_Image_Format_Specification_-_306687.pdf */
		struct FeiMrcExtHeader
		{
//			Image, System and Application
			int metadata_size;
			int metadata_version;
			unsigned int bitmask_1;
			double timestamp;
			char microscope_type[16];
			char d_number[16];
			char application[16];
			char application_version[16];

//			Gun
			double ht;
			double dose;

//			Stage
			double alpha_tilt;
			double beta_tilt;

			double x_stage;
			double y_stage;
			double z_stage;

			double tilt_axis_angle;

			double dual_axis_rotation;

//			Pixel Size
			double pixel_size_x;
			double pixel_size_y;

			char unused[48];

//			Optics
			double defocus;
			double stem_defocus;
			double applied_defocus;
			int instrument_mode;
			int projection_mode;
			char objective_lens_mode[16];
			char high_magnification_mode[16];
			int probe_mode;
			bool eftem_on;
			double magnification;	/*The magnification used in SI units (volts) */
			unsigned int bitmask_2;
			double camera_length;
			int spot_index;
			double illuminated_area;
			double intensity;
			double convergence_angle;
			char illumination_mode[16];
			bool wide_convergence_angle_range;

//			EFTEM Imaging
			bool slit_inserted;
			double slit_width;
			double acceleration_voltage_offset;
			double drift_tube_voltage;
			double energy_shift;

//			Image Shifts
			double shift_offset_x;
			double shift_offset_y;
			double shift_x;
			double shift_y;

//			Camera
			double integration_time;
			int binning_width;
			int binning_height;
			char camera_name[16];
			int readout_area_left;
			int readout_area_top;
			int readout_area_right;
			int readout_area_bottom;
			bool ceta_noise_reduction;
			int ceta_frames_summed;
			bool direct_detector_electron_counting;
			bool direct_detector_align_frames;
			int camera_param_reserved_0;
			int camera_param_reserved_1;
			int camera_param_reserved_2;
			int camera_param_reserved_3;
			unsigned int bitmask_3;
			int camera_param_reserved_4;
			int camera_param_reserved_5;
			int camera_param_reserved_6;
			int camera_param_reserved_7;
			int camera_param_reserved_8;
			int camera_param_reserved_9;
			bool phase_plate;

//			STEM
			char stem_detector_name[16];
			double gain;
			double offset;
			int stem_param_reserved_0;
			int stem_param_reserved_1;
			int stem_param_reserved_2;
			int stem_param_reserved_3;
			int stem_param_reserved_4;

//			Scan settings
			double dwell_time;
			double frame_time;
			int scan_size_left;
			int scan_size_top;
			int scan_size_right;
			int scan_size_bottom;
			double full_scan_fov_x;
			double full_scan_fov_y;

//			EDX Elemental maps
			char element[16];
			double energy_interval_lower;
			double energy_interval_higher;
			int method;

//			Dose fractions
			bool is_dose_fraction;
			int fraction_number;
			int start_frame;
			int end_frame;

//			Reconstruction
			char input_stack_filename[80];
			unsigned int bitmask_4;
			double alpha_tilt_min;
			double alpha_tilt_max;

//			FEI2 Version 2 Extension to the Extended Header Specification
			double scan_rotation;
			double diffraction_pattern_rotation;
			double image_rotation;
			int scan_mode_enumeration;
			long acquisition_time_stamp;
			char detector_commercial_name[16];
			double start_tilt_angle;
			double end_tilt_angle;
			double tilt_per_image;
			double tilt_speed;
			int beam_center_x_pixel;
			int beam_center_y_pixel;
			long cfeg_flash_timestamp;
			int phase_plate_position_index;
			char objective_aperture_name[16];

		} __attribute__((packed));

		static const char *CTF_MAGIC;
		static const char *SHORT_CTF_MAGIC;

	private:
		int mode_size;

		MrcHeader mrch;

		/* the extended MRC format for tomography, used by FEI */
		bool isFEI;

		/* for MRCS (MRC stack) format */
		bool is_stack;
		int stack_size;

		/* for MRC 8 bit packed format (2 4-bit values in each byte) */
		bool is_8_bit_packed;
		bool use_given_dimensions;

		int is_ri;
		bool is_big_endian;
		bool is_new_file;
		bool is_transpose;
		
		/** generate the machine stamp used in MRC image format. */
		static int generate_machine_stamp();
		void swap_header(MrcHeader& mrch);
		
		/** This is a utility function used to calculate new min/max/mean/sigma
		 * when write MRC file as 16 bit or 8 bit. It will write new set of 
		 * min/max/mean/sigma to mrch.
		 * this function needs get the output data storage type from mrch.mode.*/
		template<class T>
		void update_stats(const vector<T> &data);

		/** This is a utility routine to tell whether to byte swap MRC header. */
		static void check_swap(const int * data, const char * filnam, bool show_errors,
									  bool & do_swap, bool & have_err);

		int read_mrc_header(Dict & dict, int image_index = 0, const Region * area = 0, bool is_3d = false);
		int read_fei_header(Dict & dict, int image_index = 0, const Region * area = 0, bool is_3d = false);

		//utility function to transpose x and y dimension in case the source mrc image is mapc=2,mapr=1
		int transpose(float *data, int nx, int ny, int nz) const;

		template<class T>
		auto write_compressed(float *data, size_t size, int image_index, const Region* area);
	};
}

#endif	//eman__mrcio_h__
