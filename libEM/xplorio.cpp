/**
 * $Id$
 */
#include "xplorio.h"
#include "log.h"
#include "util.h"
#include "Assert.h"
#include "portable_fileio.h"
#include "geometry.h"
#include "emutil.h"

#ifdef WIN32
#include <time.h>
#define MAXPATHLEN 1024
#else
#include <sys/param.h>
#endif


using namespace EMAN;

const string XplorIO::SECTION_MODE = "ZYX";
const int XplorIO::NFLOAT_PER_LINE = 6;
const int XplorIO::INTEGER_SIZE = 8;
const int XplorIO::FLOAT_SIZE = 12;

XplorIO::XplorIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), xplor_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	is_new_file = false;
	nlines_in_header = 0;
	
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
		nlines_in_header = i;
		apix_x = cellx / nx;
		apix_y = celly / ny;
		apix_z = cellz / nz;
	}
			
	EXITFUNC;
}

bool XplorIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	if (!first_block) {
		return false;
	}
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

int XplorIO::read_header(Dict &dict, int image_index, const Region *area, bool)
{
	ENTERFUNC;

	check_read_access(image_index);
	check_region(area, FloatSize(nx, ny, nz), is_new_file);
	
	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);
	
	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["apix_x"] = apix_x;
	dict["apix_y"] = apix_y;
	dict["apix_z"] = apix_z;

	dict["XPLOR.alpha"] = cell_alpha;
	dict["XPLOR.beta"] = cell_beta;
	dict["XPLOR.gama"] = cell_gama;
	
	EXITFUNC;
	return 0;
}

int XplorIO::write_header(const Dict & dict, int image_index, const Region* area, bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index);
	if (area) {
		check_region(area, FloatSize(nx, ny, nz), is_new_file);
		EXITFUNC;
		return 0;
	}
	
	nx = dict["nx"];
	ny = dict["ny"];
	nz = dict["nz"];
	float pixel = dict["pixel"];

	nlines_in_header = 0;
	time_t t0 = time(0);
	struct tm *t = localtime(&t0);
	rewind(xplor_file);
	
	fprintf(xplor_file, "\n");
	fprintf(xplor_file, "%8d\n", 1);
	fprintf(xplor_file, "\"%s\" written by EMAN at %s\n", filename.c_str(), asctime(t));
	
	int z0 = -nz / 2;
	int z1 = (nz - 1) / 2;

	if (2 * nz - 1 == nx && 2 * nz - 1 == ny) {
		z0 = 0;
		z1 = nz - 1;
	}

	fprintf(xplor_file, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
			nx, -nx / 2, nx % 2 ? nx / 2 : nx / 2 - 1, ny, -ny / 2,
			ny % 2 ? ny / 2 : ny / 2 - 1, nz, z0, z1);
	fprintf(xplor_file, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", 
			nx * pixel, ny * pixel, nz * pixel, 90.0, 90.0, 90.0);
	fprintf(xplor_file, "ZYX\n");
	nlines_in_header = 6;
	
	EXITFUNC;
	return 0;
}

int XplorIO::read_data(float *data, int image_index, const Region *area, bool)
{
	ENTERFUNC;
	
	check_read_access(image_index);
	check_region(area, FloatSize(nx, ny, nz), is_new_file);

	Assert(nlines_in_header > 0);
	rewind(xplor_file);
	jump_lines(xplor_file, nlines_in_header);

	
	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);
	
	int nlines_per_sec = (nx *ny) / NFLOAT_PER_LINE;
	int nitems_last_line = (nx * ny) % NFLOAT_PER_LINE;
	if (nitems_last_line != 0) {
		nlines_per_sec++;
	}
	nlines_per_sec++; // image index line
	
	int x0 = 0;
	int y0 = 0;
	int z0 = 0;

	if (area) {
		x0 = (int)area->origin[0];
		y0 = (int)area->origin[1];
		z0 = (int)area->origin[2];
	}
	if (z0 > 0) {
		jump_lines(xplor_file, z0 * nlines_per_sec);
	}

	char line[MAXPATHLEN];
	size_t nlines_pre_sec = (y0 * nx + x0) / NFLOAT_PER_LINE;
	int head_nitems = (y0 * nx + x0) % NFLOAT_PER_LINE;
	int gap_nitems = nx - xlen;
	int ti = 0;
	
	for (int k = 0; k < zlen; k++) {
		
		if (!fgets(line, sizeof(line), xplor_file)) {
			throw ImageReadException(filename, "read xplor file failed");
		}
		int kk = 0;
		sscanf(line, "%d", &kk);

		if (kk != (k+1+z0)) {
			char desc[256];
			sprintf(desc, "section index = %d. It should be %d", kk, (k+1+z0));
			throw ImageReadException(filename, desc);
		}

		if (!area) {
			Assert(ti == (k*xlen*ylen));
		}
		
		if (nlines_pre_sec > 0) {
			jump_lines(xplor_file, nlines_pre_sec);
		}
		
		int tail_nitems = 0;
		bool is_head_read = false;
		
		for (int j = 0; j < ylen; j++) {
			if (head_nitems > 0 && !is_head_read) {
				read_numbers(xplor_file, NFLOAT_PER_LINE-head_nitems,
							 NFLOAT_PER_LINE-1, data, &ti);
			}
			
			read_lines(xplor_file, (xlen - head_nitems), data, &ti);
			tail_nitems = (xlen - head_nitems) % NFLOAT_PER_LINE;
			
			if ((gap_nitems + tail_nitems) > 0) {
				head_nitems = NFLOAT_PER_LINE -
					(gap_nitems + tail_nitems) % NFLOAT_PER_LINE;
			}
			else {
				head_nitems = 0;
			}
			
			is_head_read = false;
			
			if (tail_nitems > 0) {
				if (gap_nitems < (NFLOAT_PER_LINE-tail_nitems) && j != (ylen-1)) {
					not_read_numbers(xplor_file, tail_nitems,
									 tail_nitems+gap_nitems, data, &ti);
					is_head_read = true;
				}
				else {
					read_numbers(xplor_file, 0, tail_nitems-1, data, &ti);
				}
			}

			if (gap_nitems > (NFLOAT_PER_LINE-tail_nitems)) {
				int gap_nlines = (gap_nitems - (NFLOAT_PER_LINE-tail_nitems)) /
					NFLOAT_PER_LINE;
				if (gap_nlines > 0 && j != (ylen-1)) {
					jump_lines(xplor_file, gap_nlines);
				}
			}
		}
		
		int ytail_nitems = (ny-ylen-y0) * nx + (nx-xlen-x0) - (NFLOAT_PER_LINE-tail_nitems);
		jump_line_by_items(xplor_file, ytail_nitems);
	}
			
	
	EXITFUNC;
	return 0;
}


#if 0
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
		int kk = 0;
		sscanf(line, "%d", &kk);
		if (kk != (k+1)) {
			LOGERR("section index = %d. It should be %d\n", kk, (k+1));
		}
		
		int k2 = k * nxy;
		
		for (int i = 0; i < nlines; i++) {
			fgets(line, sizeof(line), xplor_file);
			int i2 = k2 + i * step;
			sscanf(line, "%f %f %f %f %f %f",
				   &data[i2], &data[i2+1], &data[i2+2], 
				   &data[i2+3], &data[i2+4], &data[i2+5]);
		}
		
		if (nleft > 0) {
			int i2 = k2 + nlines * step;
			fgets(line, sizeof(line), xplor_file);
			char *pline = line;
			for (int j = 0; j < nleft; j++) {
				sscanf(pline, "%f", &data[i2+j]);
				pline += FLOAT_SIZE;
			}
		}
	}
			
	
	EXITFUNC;
	return 0;
}
#endif


int XplorIO::write_data(float *data, int image_index, const Region* area, bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index, 1, data);
	check_region(area, FloatSize(nx,ny,nz), is_new_file);

	rewind(xplor_file);
	jump_lines(xplor_file, nlines_in_header);
	
	int nsecs = nx * ny;
	int step = NFLOAT_PER_LINE;
	
	if (!area) {
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

		// not sure what this is doing. so comment out.
		//fprintf(xplor_file, "%8d\n", -9999); 
	}
	else {
		LOGERR("not implemented yet");
	}
	
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
	init();
	return is_big_endian;
}

void XplorIO::jump_line_by_items(FILE * xplor_file, int nitems)
{
	Assert(xplor_file);
	
	if (nitems <= 0) {
		return;
	}
	
	int nlines = nitems / NFLOAT_PER_LINE;
	if ((nitems % NFLOAT_PER_LINE) != 0) {
		nlines++;
	}
	if (nlines > 0) {
		jump_lines(xplor_file, nlines);
	}
}


void XplorIO::jump_lines(FILE * xplor_file, int nlines)
{
	Assert(xplor_file);
	
	if (nlines > 0) {
		char line[MAXPATHLEN];
		for (int l = 0; l < nlines; l++) {
			if (!fgets(line, sizeof(line), xplor_file)) {
				Assert("read xplor file failed");
			}
		}
	}
}

void XplorIO::read_numbers(FILE * xplor_file, int start, int end,
						   float *data, int *p_i)
{
	Assert(start >= 0);
	Assert(start <= end);
	Assert(end <= NFLOAT_PER_LINE);
	Assert(xplor_file);
	Assert(data);

	char line[MAXPATHLEN];
	float d[NFLOAT_PER_LINE];
	if (!fgets(line, sizeof(line), xplor_file)) {
		Assert("read xplor file failed");
	}
	
	sscanf(line, "%f %f %f %f %f %f",
		   &d[0], &d[1], &d[2], &d[3], &d[4],&d[5]);
	
	for (int i = start; i <= end; i++) {
		data[*p_i] = d[i];
		(*p_i)++;
	}
}


void XplorIO::not_read_numbers(FILE * xplor_file, int start, int end,
							   float * data, int *p_i)
{
	Assert(start >= 0);
	Assert(start <= end);
	Assert(end <= NFLOAT_PER_LINE);
	Assert(xplor_file);
	Assert(data);

	char line[MAXPATHLEN];
	float d[NFLOAT_PER_LINE];
	if (!fgets(line, sizeof(line), xplor_file)) {
		Assert("read xplor file failed");
	}
	
	sscanf(line, "%f %f %f %f %f %f",
		   &d[0], &d[1], &d[2], &d[3], &d[4],&d[5]);
	
	for (int i = 0; i < start; i++) {
		data[*p_i] = d[i];
		(*p_i)++;
	}

	for (int i = end; i < NFLOAT_PER_LINE; i++) {
		data[*p_i] = d[i];
		(*p_i)++;
	}
	
}

void XplorIO::read_lines(FILE * xplor_file, int nitems, float *data, int *p_i)
{
	Assert(xplor_file);
	Assert(data);
	Assert(p_i);
	
	if (nitems > 0) {
		int nlines = nitems / NFLOAT_PER_LINE;
		for (int i = 0; i < nlines; i++) {
			read_numbers(xplor_file, 0, NFLOAT_PER_LINE-1, data, p_i);
		}
	}
}

