/**
 * $Id$
 */
#include "xplorio.h"
#include "log.h"
#include "util.h"
#ifdef WIN32
#include <time.h>
#endif
#include "portable_fileio.h"

using namespace EMAN;

string XplorIO::SECTION_MODE = "ZYX";
int XplorIO::NFLOAT_PER_LINE = 6;
int XplorIO::INTEGER_SIZE = 8;
int XplorIO::FLOAT_SIZE = 12;

XplorIO::XplorIO(string file, IOMode rw)
:	filename(file), rw_mode(rw), xplor_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	nx = 0;
	ny = 0;
	nz = 0;

	apix_x = 0;
	apix_y = 0;
	apix_z = 0;

	cell_alpha = 0;
	cell_beta = 0;
	cell_gama = 0;
}

XplorIO::~XplorIO()
{
	if (xplor_file) {
		fclose(xplor_file);
		xplor_file = 0;
	}
}

void XplorIO::init()
{
	if (initialized) {
		return;
	}

	ENTERFUNC;
	initialized = true;

	bool is_new_file = false;
	xplor_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		char first_block[1024];
		fread(&first_block, sizeof(char), sizeof(first_block), xplor_file);
		if (!is_valid(&first_block)) {
			throw ImageReadException(filename, "invalid XPLOR");
		}
		portable_fseek(xplor_file, 0, SEEK_SET);
		char line[1024];
		int i = 1;
		int ntitle = 0;

		int xmin = 0;
		int xmax = 0;
		int ymin = 0;
		int ymax = 0;
		int zmin = 0;
		int zmax = 0;

		float cellx = 0;
		float celly = 0;
		float cellz = 0;
		
		while(fgets(line, sizeof(line), xplor_file)) {
			line[strlen(line)-1] = '\0';
			if (i == 2) {
				ntitle = atoi(line);
			}
			else if (i == (ntitle+3)) {
				if (sscanf(line, "%8d%8d%8d%8d%8d%8d%8d%8d%8d", &nx, &xmin, &xmax,
						   &ny, &ymin, &ymax, &nz, &zmin, &zmax) != 9) {
					throw ImageReadException(filename, "invalid XPLOR");
				}
			}
			else if (i == (ntitle+4)) {
				if(sscanf(line, "%f %f %f %f %f %f",
						  &cellx, &celly, &cellz, &cell_alpha, &cell_beta, &cell_gama) != 6) {
					throw ImageReadException(filename, "invalid XPLOR");
				}
			}
			else if (i == (ntitle+5)) {
				break;
			}
			
			i++;
		}
		apix_x = cellx / nx;
		apix_y = celly / ny;
		apix_z = cellz / nz;
	}
			
	EXITFUNC;
}

bool XplorIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	char *buf = (char *)(first_block);
	string line1 = Util::get_line_from_string(&buf);
	bool result = true;
	
	if (line1.size() != 0) {
		result = false;
	}
	else {
		string line2 = Util::get_line_from_string(&buf);
		int ntitle = 0;
	
		if ((int)line2.size() != INTEGER_SIZE) {
			result = false;
		}
		else {
			ntitle = atoi(line2.c_str());
			if (ntitle < 0 || ntitle > 50) {
				result = false;
			}

			else {
				for (int i = 0; i < ntitle+2; i++) {
					Util::get_line_from_string(&buf);
				}
				
				string modeline = Util::get_line_from_string(&buf);
				if (modeline != SECTION_MODE) {
					result = false;
				}
			}
		}
	}
	
	EXITFUNC;
	return result;
}

int XplorIO::read_header(Dict &dict, int, const Region *, bool)
{
	ENTERFUNC;

	dict["nx"] = nx;
	dict["ny"] = ny;
	dict["nz"] = nz;

	dict["apix_x"] = apix_x;
	dict["apix_y"] = apix_y;
	dict["apix_z"] = apix_z;

	dict["XPLOR.alpha"] = cell_alpha;
	dict["XPLOR.beta"] = cell_beta;
	dict["XPLOR.gama"] = cell_gama;
	
	EXITFUNC;
	return 0;
}

int XplorIO::write_header(const Dict & dict, int image_index, const Region* , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index);

	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];
	float pixel = dict["pixel"];

	time_t t0 = time(0);
	struct tm *t = localtime(&t0);

	fprintf(xplor_file, "\n%8d\n\"%s\" written by EMAN at %s", 1, filename.c_str(), asctime(t));

	int z0 = -nz / 2;
	int z1 = (nz - 1) / 2;

	if (2 * nz - 1 == nx && 2 * nz - 1 == ny) {
		z0 = 0;
		z1 = nz - 1;
	}

	fprintf(xplor_file, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
			nx, -nx / 2, nx % 2 ? nx / 2 : nx / 2 - 1, ny, -ny / 2,
			ny % 2 ? ny / 2 : ny / 2 - 1, nz, z0, z1);
	fprintf(xplor_file, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n", nx * pixel, ny * pixel,
			nz * pixel, 90.0, 90.0, 90.0);
	EXITFUNC;
	return 0;
}

int XplorIO::read_data(float *data, int, const Region *, bool)
{
	ENTERFUNC;
	int step = NFLOAT_PER_LINE;
	char line[1024];
	int nxy = nx * ny;
	int nlines = nxy / step;
	int nleft = nxy - nlines * step;
	
	for (int k = 0; k < nz; k++) {
		fgets(line, sizeof(line), xplor_file);
		int i = 0;
		int k2 = k * nxy;
		while (fgets(line, sizeof(line), xplor_file)) {
			line[strlen(line)-1] = '\0';
			int i2 = k2 + i * step;
			
			if (i == nlines) {
				char *pline = line;
				for (int j = 0; j < nleft; j++) {
					sscanf(pline, "%f", &data[i2+j]);
					pline += FLOAT_SIZE;
				}
			}
			else {
				sscanf(line, "%f %f %f %f %f %f",
					   &data[i2], &data[i2+1], &data[i2+2], 
					   &data[i2+3], &data[i2+4], &data[i2+5]);
				i++;
			}
		}
	}
			
	
	EXITFUNC;
	return 0;
}

int XplorIO::write_data(float *data, int image_index, const Region* , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index, 1, data);

	int nsecs = nx * ny;
	int step = NFLOAT_PER_LINE;

	for (int k = 0; k < nz; k++) {
		fprintf(xplor_file, "%8d\n", k);

		for (int i = 0; i < nsecs - step; i += step) {
			for (int j = 0; j < step; j++) {
				fprintf(xplor_file, "%12.5E", data[k * nsecs + i + j]);
			}
			fprintf(xplor_file, "\n");
		}

		for (int l = (nsecs - 1) / step * step; l < nsecs; l++) {
			fprintf(xplor_file, "%12.5E", data[k * nsecs + l]);
		}

		fprintf(xplor_file, "\n");
	}

	fprintf(xplor_file, "%8d\n", -9999);
	EXITFUNC;
	
	return 0;
}

void XplorIO::flush()
{
	fflush(xplor_file);
}

bool XplorIO::is_complex_mode()
{
	return false;
}

bool XplorIO::is_image_big_endian()
{
	return is_big_endian;
}

