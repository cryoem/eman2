#ifndef eman__mrcio_h__
#define eman__mrcio_h__ 1

#include "imageio.h"
#include <stdio.h>

// todo: change normalization in read_data/write_data
//       after EMData defins Min() and Max()
// todo: define a function to return 'mode' in string

namespace EMAN
{
    class MrcIO : public ImageIO
    {
    public:
	MrcIO(string filename, IOMode rw_mode = READ_ONLY);
	~MrcIO();

	DEFINE_IMAGEIO_FUNC;

	int read_ctf(Ctf & ctf, int image_index = 0);
	int write_ctf(const Ctf & ctf, int image_index = 0);

	static bool is_valid(const void *first_block, off_t file_size = 0);
	static int get_mode_size(int mm);
	static int to_em_datatype(int mrcmode);
	static int to_mrcmode(int em_datatype, bool is_complex);

    private:
	enum MrcMode {
	    MRC_UCHAR = 0,
	    MRC_USHORT,
	    MRC_FLOAT,
	    MRC_USHORT_COMPLEX,
	    MRC_FLOAT_COMPLEX,
	    MRC_UNKNOWN
	};

	enum {
	    MRC_NUM_LABELS = 10,
	    MRC_LABEL_SIZE = 80,
	    NUM_4BYTES_PRE_MAP = 52,
	    NUM_4BYTES_AFTER_MAP = 3
	};

	/* updated to MRC Image2000 format which is compatible with CCP4 format */
	struct MrcHeader
	{
	    int nx;		/* number of columns */
	    int ny;		/* number of rows */
	    int nz;		/* number of sections */

	    int mode;		/* See modes above. */

	    int nxstart;	/* No. of first column in map, default 0. */
	    int nystart;	/* No. of first row in map, default 0. */
	    int nzstart;	/* No. of first section in map,default 0. */

	    int mx;		/* Number of intervals along X. */
	    int my;		/* Number of intervals along Y. */
	    int mz;		/* Number of intervals along Z. */
	    
	    /* Cell: treat a whole 2D image as a cell */
	    float xlen;		/* Cell dimensions (Angstroms). */
	    float ylen;		/* Cell dimensions (Angstroms). */
	    float zlen;		/* Cell dimensions (Angstroms). */

	    float alpha;	/* Cell angles (Degrees). */
	    float beta;		/* Cell angles (Degrees). */
	    float gamma;	/* Cell angles (Degrees). */

	    /* axis X => 1, Y => 2, Z => 3 */
	    int mapc;		/* Which axis corresponds to Columns.  */
	    int mapr;		/* Which axis corresponds to Rows.     */
	    int maps;		/* Which axis corresponds to Sections. */

	    float amin;		/* Minimum density value. */
	    float amax;		/* Maximum density value. */
	    float amean;	/* Mean density value.    */

	    int ispg;		/* Space group number (0 for images). */

	    int nsymbt;		/* Number of chars used for storing symmetry operators. */

	    int user[25];

	    float xorigin;	/* X origin. */
	    float yorigin;	/* Y origin. */
	    float zorigin;	/* Y origin. */

	    char map[4];	/* constant string "MAP "  */
	    int machinestamp;	/* machine stamp in CCP4 convention: big endian=0x11110000 little endian=0x4444000 */

	    float rms;		/* rms deviation of map from mean density */

	    int nlabels;	/* Number of labels being used. */
	    char labels[MRC_NUM_LABELS][MRC_LABEL_SIZE];
	};

	static const char *CTF_MAGIC;
	static const char *SHORT_CTF_MAGIC;


    private:
	string filename;
	IOMode rw_mode;
	FILE *mrcfile;
	MrcHeader mrch;
	int mode_size;

	bool is_big_endian;
	bool is_new_file;
	bool initialized;
    };
}

#endif
