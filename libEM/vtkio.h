/**
 * $Id$
 */
#ifndef eman__vtkio_h__
#define eman__vtkio_h__ 1

#include "imageio.h"
#include <stdio.h>

namespace EMAN
{

	class VtkIO:public ImageIO
	{
	  public:
		VtkIO(string filename, IOMode rw_mode = READ_ONLY);
		~VtkIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		static const char *MAGIC;

		enum VtkType
		{
			VTK_ASCII,
			VTK_BINARY,
			VTK_UNKNOWN
		};

		enum DataType
		{
			UNSIGNED_SHORT,
			FLOAT,
			DATA_UNKNOWN
		};

		int to_em_datatype(int vtk_datatype);
		int get_mode_size(DataType d);

	  private:
		  string filename;
		IOMode rw_mode;
		FILE *vtk_file;
		bool is_big_endian;
		bool initialized;

		DataType datatype;
		VtkType filetype;
		int nx;
		int ny;
		int nz;
		float originx;
		float originy;
		float originz;
		float spacingx;
		float spacingy;
		float spacingz;
		off_t file_offset;
	};

}

#endif
