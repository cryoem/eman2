/**
 * $Id$
 */
#ifndef eman__emio_h__
#define eman__emio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
	/** EmIO defines I/O operations on EM image format.
     * EM image = header + data with (data = nx * ny * nz).
	 *
	 * An EM image file stores 1 single 2D or 3D image.
     */
	class EmIO:public ImageIO
	{
	  public:
		EmIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~EmIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block, off_t file_size = 0);
		static size_t get_mode_size(char data_type);
		static int get_machine_type();
		static int to_em_datatype(char t);

	  private:
		struct EMHeader
		{
			char machine;		// 0=OS-8, 1=VAX; 2=Convex; 3=SGI; 5=Mac; 6=PC
			char is_new_ver;	// OS-9 only
			char not_used1;
			char data_type;		// 1=byte, 2=short; 4=int; 5=float; 8=complex; 9=double
			int nx;
			int ny;
			int nz;
			char comment[80];
			int parameters[40];
			char username[20];
			char date[8];
			char userdata[228];
		};

		enum DataType
		{
			EM_EM_CHAR = 1,
			EM_EM_SHORT = 2,
			EM_EM_INT = 4,
			EM_EM_FLOAT = 5,
			EM_EM_COMPLEX = 8,
			EM_EM_DOUBLE = 9,
			EM_EM_UNKNOWN
		};

		enum MachineType
		{
			EM_OS8 = 0,
			EM_VAX = 1,
			EM_CONVEX = 2,
			EM_SGI = 3,
			EM_MAC = 5,
			EM_PC = 6,
			EM_UNKNOWN_MACHINE
		};

	  private:
		string filename;
		IOMode rw_mode;
		FILE *em_file;
		EMHeader emh;

		size_t mode_size;
		DataType mode;
		bool is_big_endian;
		bool initialized;
		bool is_new_file;
	};

}


#endif
