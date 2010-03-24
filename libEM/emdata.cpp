/**
 * $Id$
 */

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

#include "emdata.h"
#include "all_imageio.h"
#include "ctf.h"
#include "processor.h"
#include "aligner.h"
#include "cmp.h"  // comparision method header file
#include "emfft.h"
#include "projector.h"
#include "geometry.h"
#include "mrc.h" // ming add

#include <log.h>//ming
#include "Euler.h" // ming add
#include "List.h" //ming
#include "util.h" // ming add

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#include <iomanip>
#include <complex>

#include <algorithm> // fill
#include <cmath>


#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

/****************Ming begin***************/
#ifdef WIN32
	#include <windows.h>
	#define MAXPATHLEN (MAX_PATH*4)
#else
#include <sys/param.h>
#endif	//WIN32

/************Ming end*******************/


#define EMDATA_EMAN2_DEBUG 0

#ifdef EMAN2_USING_CUDA
#include <driver_functions.h>
#include "cuda/cuda_processor.h"
#endif // EMAN2_USING_CUDA

using namespace EMAN;
using namespace std;
using namespace boost;


int EMData::totalalloc=0;		// mainly used for debugging/memory leak purposes



EMData::EMData() :
#ifdef EMAN2_USING_CUDA
		cuda_cache_handle(-1),
#endif //EMAN2_USING_CUDA
		attr_dict(), rdata(0), supp(0), flags(0), changecount(0), nx(0), ny(0), nz(0), nxy(0), nxyz(0), xoff(0), yoff(0),
		zoff(0), all_translation(),	path(""), pathnum(0), rot_fp(0)

{
	ENTERFUNC;

	attr_dict["apix_x"] = 1.0f;
	attr_dict["apix_y"] = 1.0f;
	attr_dict["apix_z"] = 1.0f;

	attr_dict["is_complex"] = int(0);
	attr_dict["is_complex_x"] = int(0);
	attr_dict["is_complex_ri"] = int(1);

	attr_dict["datatype"] = (int)EMUtil::EM_FLOAT;

	EMData::totalalloc++;
#ifdef MEMDEBUG
	printf("EMDATA+  %4d    %p\n",EMData::totalalloc,this);
#endif

}

EMData::EMData(const string& filename, int image_index) :
#ifdef EMAN2_USING_CUDA
		cuda_cache_handle(-1),
#endif //EMAN2_USING_CUDA
		attr_dict(), rdata(0), supp(0), flags(0), changecount(0), nx(0), ny(0), nz(0), nxy(0), nxyz(0), xoff(0), yoff(0), zoff(0),
		all_translation(),	path(filename), pathnum(image_index), rot_fp(0)
{
	ENTERFUNC;

	attr_dict["apix_x"] = 1.0f;
	attr_dict["apix_y"] = 1.0f;
	attr_dict["apix_z"] = 1.0f;

	attr_dict["is_complex"] = int(0);
	attr_dict["is_complex_x"] = int(0);
	attr_dict["is_complex_ri"] = int(1);

	this->read_image(filename, image_index);

	update();
	EMData::totalalloc++;

	EXITFUNC;
}

EMData::EMData(const EMData& that) :
#ifdef EMAN2_USING_CUDA
		cuda_cache_handle(-1),
#endif //EMAN2_USING_CUDA
		attr_dict(that.attr_dict), rdata(0), supp(0), flags(that.flags), changecount(that.changecount), nx(that.nx), ny(that.ny), nz(that.nz),
		nxy(that.nx*that.ny), nxyz(that.nx*that.ny*that.nz), xoff(that.xoff), yoff(that.yoff), zoff(that.zoff),all_translation(that.all_translation),	path(that.path),
		pathnum(that.pathnum), rot_fp(0)
{
	ENTERFUNC;

	float* data = that.rdata;
	size_t num_bytes = nx*ny*nz*sizeof(float);
	if (data && num_bytes != 0)
	{
		rdata = (float*)EMUtil::em_malloc(num_bytes);
		EMUtil::em_memcpy(rdata, data, num_bytes);
	}
#ifdef EMAN2_USING_CUDA
	if (num_bytes != 0 && that.cuda_cache_handle != -1 && that.gpu_rw_is_current()) {
		float * cuda_data = that.get_cuda_data();
		CudaDataLock lock2(&that);
		set_size_cuda(nx,ny,nz);
		float *cd = get_cuda_data();
		CudaDataLock lock1(this);
		cudaError_t error = cudaMemcpy(cd,cuda_data,num_bytes,cudaMemcpyDeviceToDevice);
		if ( error != cudaSuccess ) throw UnexpectedBehaviorException("cudaMemcpy failed in EMData copy construction with error: " + string(cudaGetErrorString(error)));
	}
	// This is a bit of hack
	flags |= EMDATA_GPU_RO_NEEDS_UPDATE;
#endif //EMAN2_USING_CUDA
	if (that.rot_fp != 0) rot_fp = new EMData(*(that.rot_fp));

	EMData::totalalloc++;

	ENTERFUNC;
}

EMData& EMData::operator=(const EMData& that)
{
	ENTERFUNC;

	if ( this != &that )
	{
		free_memory(); // Free memory sets nx,ny and nz to 0

		// Only copy the rdata if it exists, we could be in a scenario where only the header has been read
		float* data = that.rdata;
		size_t num_bytes = that.nx*that.ny*that.nz*sizeof(float);
		if (data && num_bytes != 0)
		{
			nx = 1; // This prevents a memset in set_size
			set_size(that.nx,that.ny,that.nz);
			EMUtil::em_memcpy(rdata, data, num_bytes);
		}

		flags = that.flags;

		all_translation = that.all_translation;

		path = that.path;
		pathnum = that.pathnum;
		attr_dict = that.attr_dict;

		xoff = that.xoff;
		yoff = that.yoff;
		zoff = that.zoff;

#ifdef EMAN2_USING_CUDA
		free_cuda_memory();
		// There should also be the case where we deal with ro data...
		if (num_bytes != 0 && that.cuda_cache_handle != -1 && that.gpu_rw_is_current()) {
			float * cuda_data = that.get_cuda_data();
			CudaDataLock lock1(&that);
			set_size_cuda(that.nx, that.ny, that.nz);
			float *cd = get_cuda_data();
			CudaDataLock lock2(this);
			cudaError_t error = cudaMemcpy(cd,cuda_data,num_bytes,cudaMemcpyDeviceToDevice);
			if ( error != cudaSuccess ) throw UnexpectedBehaviorException("cudaMemcpy failed in operator= with error: " + string(cudaGetErrorString(error)));
		}
		// This is a bit of hack
		flags &= EMDATA_GPU_RO_NEEDS_UPDATE;
#endif //EMAN2_USING_CUDA

		changecount = that.changecount;

		if (that.rot_fp != 0) rot_fp = new EMData(*(that.rot_fp));
		else rot_fp = 0;
	}
	EXITFUNC;
	return *this;
}

EMData::EMData(int nx, int ny, int nz, bool is_real) :
#ifdef EMAN2_USING_CUDA
		cuda_cache_handle(-1),
#endif //EMAN2_USING_CUDA
		attr_dict(), rdata(0), supp(0), flags(0), changecount(0), nx(0), ny(0), nz(0), nxy(0), nxyz(0), xoff(0), yoff(0), zoff(0),
		all_translation(),	path(""), pathnum(0), rot_fp(0)
{
	ENTERFUNC;

	// used to replace cube 'pixel'
	attr_dict["apix_x"] = 1.0f;
	attr_dict["apix_y"] = 1.0f;
	attr_dict["apix_z"] = 1.0f;

	if(is_real) {	// create a real image [nx, ny, nz]
		attr_dict["is_complex"] = int(0);
		attr_dict["is_complex_x"] = int(0);
		attr_dict["is_complex_ri"] = int(1);
		set_size(nx, ny, nz);
	}
	else {	//create a complex image which real dimension is [ny, ny, nz]
		int new_nx = nx + 2 - nx%2;
		set_size(new_nx, ny, nz);

		attr_dict["is_complex"] = int(1);

		if(ny==1 && nz ==1)	{
			attr_dict["is_complex_x"] = int(1);
		}
		else {
			attr_dict["is_complex_x"] = int(0);
		}

		attr_dict["is_complex_ri"] = int(1);
		attr_dict["is_fftpad"] = int(1);

		if(nx%2 == 1) {
			attr_dict["is_fftodd"] = 1;
		}
	}

	this->to_zero();
	update();
	EMData::totalalloc++;

	EXITFUNC;
}


EMData::EMData(float* data, const int x, const int y, const int z, const Dict& attr_dict) :
#ifdef EMAN2_USING_CUDA
		cuda_cache_handle(-1),
#endif //EMAN2_USING_CUDA
		attr_dict(attr_dict), rdata(data), supp(0), flags(0), changecount(0), nx(x), ny(y), nz(z), nxy(x*y), nxyz(x*y*z), xoff(0),
		yoff(0), zoff(0), all_translation(), path(""), pathnum(0), rot_fp(0)
{
	ENTERFUNC;

	// used to replace cube 'pixel'
	attr_dict["apix_x"] = 1.0f;
	attr_dict["apix_y"] = 1.0f;
	attr_dict["apix_z"] = 1.0f;

	EMData::totalalloc++;

	update();
	EXITFUNC;
}

//debug
using std::cout;
using std::endl;
EMData::~EMData()
{
	ENTERFUNC;
	free_memory();
#ifdef EMAN2_USING_CUDA
// 	cout << "Death comes to " << cuda_cache_handle << " " << this << endl;
	free_cuda_memory();
#endif // EMAN2_USING_CUDA
	EMData::totalalloc--;
#ifdef MEMDEBUG
	printf("EMDATA-  %4d    %p\n",EMData::totalalloc,this);
#endif
	EXITFUNC;
}

void EMData::clip_inplace(const Region & area,const float& fill_value)
{
	// Added by d.woolford
	ENTERFUNC;

	// Store the current dimension values
	int prev_nx = nx, prev_ny = ny, prev_nz = nz;
	int prev_size = nx*ny*nz;

	// Get the zsize, ysize and xsize of the final area, these are the new dimension sizes of the pixel data
	int new_nz = ( area.size[2]==0 ? 1 : (int)area.size[2]);
	int new_ny = ( area.size[1]==0 ? 1 : (int)area.size[1]);
	int new_nx = (int)area.size[0];

	if ( new_nz < 0 || new_ny < 0 || new_nx < 0 )
	{
		// Negative image dimensions were never tested nor considered when creating this implementation
		throw ImageDimensionException("New image dimensions are negative - this is not supported in the clip_inplace operation");
	}

	int new_size = new_nz*new_ny*new_nx;

	// Get the translation values, they are used to construct the ClipInplaceVariables object
	int x0 = (int) area.origin[0];
	int y0 = (int) area.origin[1];
	int z0 = (int) area.origin[2];

	// Get a object that calculates all the interesting variables associated with the clip inplace operation
	ClipInplaceVariables civ(prev_nx, prev_ny, prev_nz, new_nx, new_ny, new_nz, x0, y0, z0);

	get_data(); // Do this here to make sure rdata is up to date, applicable if GPU stuff is occuring
	// Double check to see if any memory shifting even has to occur
	if ( x0 > prev_nx || y0 > prev_ny || z0 > prev_nz || civ.x_iter == 0 || civ.y_iter == 0 || civ.z_iter == 0)
	{
		// In this case the volume has been shifted beyond the location of the pixel rdata and
		// the client should expect to see a volume with nothing in it.

		// Set size calls realloc,
		set_size(new_nx, new_ny, new_nz);

		// Set pixel memory to zero - the client should expect to see nothing
		EMUtil::em_memset(rdata, 0, new_nx*new_ny*new_nz);

		return;
	}

	// Resize the volume before memory shifting occurs if the new volume is larger than the previous volume
	// All of the pixel rdata is guaranteed to be at the start of the new volume because realloc (called in set size)
	// guarantees this.
	if ( new_size > prev_size )
		set_size(new_nx, new_ny, new_nz);

	// Store the clipped row size.
	size_t clipped_row_size = (civ.x_iter) * sizeof(float);

	// Get the new sector sizes to save multiplication later.
	int new_sec_size = new_nx * new_ny;
	int prev_sec_size = prev_nx * prev_ny;

	// Determine the memory locations of the source and destination pixels - at the point nearest
	// to the beginning of the volume (rdata)
	int src_it_begin = civ.prv_z_bottom*prev_sec_size + civ.prv_y_front*prev_nx + civ.prv_x_left;
	int dst_it_begin = civ.new_z_bottom*new_sec_size + civ.new_y_front*new_nx + civ.new_x_left;

	// This loop is in the forward direction (starting at points nearest to the beginning of the volume)
	// it copies memory only when the destination pointer is less the source pointer - therefore
	// ensuring that no memory "copied to" is yet to be "copied from"
	for (int i = 0; i < civ.z_iter; ++i) {
		for (int j = 0; j < civ.y_iter; ++j) {

			// Determine the memory increments as dependent on i and j
			// This could be optimized so that not so many multiplications are occurring...
			int dst_inc = dst_it_begin + j*new_nx + i*new_sec_size;
			int src_inc = src_it_begin + j*prev_nx + i*prev_sec_size;
			float* local_dst = rdata + dst_inc;
			float* local_src = rdata + src_inc;

			if ( dst_inc >= src_inc )
			{
				// this is fine, it will happen now and then and it will be necessary to continue.
				// the tempatation is to break, but you can't do that (because the point where memory intersects
				// could be in this slice - and yes, this aspect could be optimized).
				continue;
			}

			// Asserts are compiled only in debug mode
			// This situation not encountered in testing thus far
			Assert( dst_inc < new_size && src_inc < prev_size && dst_inc >= 0 && src_inc >= 0 );

			// Finally copy the memory
			EMUtil::em_memcpy(local_dst, local_src, clipped_row_size);
		}
	}

	// Determine the memory locations of the source and destination pixels - at the point nearest
	// to the end of the volume (rdata+new_size)
	int src_it_end = prev_size - civ.prv_z_top*prev_sec_size - civ.prv_y_back*prev_nx - prev_nx + civ.prv_x_left;
	int dst_it_end = new_size - civ.new_z_top*new_sec_size - civ.new_y_back*new_nx - new_nx + civ.new_x_left;

	// This loop is in the reverse direction (starting at points nearest to the end of the volume).
	// It copies memory only when the destination pointer is greater than  the source pointer therefore
	// ensuring that no memory "copied to" is yet to be "copied from"
	for (int i = 0; i < civ.z_iter; ++i) {
		for (int j = 0; j < civ.y_iter; ++j) {

			// Determine the memory increments as dependent on i and j
			int dst_inc = dst_it_end - j*new_nx - i*new_sec_size;
			int src_inc = src_it_end - j*prev_nx - i*prev_sec_size;
			float* local_dst = rdata + dst_inc;
			float* local_src = rdata + src_inc;

			if (dst_inc <= (src_inc + civ.x_iter ))
			{
				// Overlap
				if ( dst_inc > src_inc )
				{
					// Because the memcpy operation is the forward direction, and this "reverse
					// direction" loop is proceeding in a backwards direction, it is possible
					// that memory copied to is yet to be copied from (because memcpy goes forward).
					// In this scenario pixel memory "copied to" is yet to be "copied from"
					// i.e. there is overlap

					// memmove handles overlapping cases.
					// memmove could use a temporary buffer, or could go just go backwards
					// the specification doesn't say how the function behaves...
					// If memmove creates a temporary buffer is clip_inplace no longer inplace?
					memmove(local_dst, local_src, clipped_row_size);
				}
				continue;
			}

			// This situation not encountered in testing thus far
			Assert( dst_inc < new_size && src_inc < prev_size && dst_inc >= 0 && src_inc >= 0 );

			// Perform the memory copy
			EMUtil::em_memcpy(local_dst, local_src, clipped_row_size);
		}
	}

	// Resize the volume after memory shifting occurs if the new volume is smaller than the previous volume
	// set_size calls realloc, guaranteeing that the pixel rdata is in the right location.
	if ( new_size < prev_size )
		set_size(new_nx, new_ny, new_nz);

	// Now set all the edges to zero

	// Set the extra bottom z slices to the fill_value
	if (  z0 < 0 )
	{
		//EMUtil::em_memset(rdata, 0, (-z0)*new_sec_size*sizeof(float));
		int inc = (-z0)*new_sec_size;
		std::fill(rdata,rdata+inc,fill_value);
	}

	// Set the extra top z slices to the fill_value
	if (  civ.new_z_top > 0 )
	{
		float* begin_pointer = rdata + (new_nz-civ.new_z_top)*new_sec_size;
		//EMUtil::em_memset(begin_pointer, 0, (civ.new_z_top)*new_sec_size*sizeof(float));
		float* end_pointer = begin_pointer+(civ.new_z_top)*new_sec_size;
		std::fill(begin_pointer,end_pointer,fill_value);
	}

	// Next deal with x and y edges by iterating through each slice
	for ( int i = civ.new_z_bottom; i < civ.new_z_bottom + civ.z_iter; ++i )
	{
		// Set the extra front y components to the fill_value
		if ( y0 < 0 )
		{
			float* begin_pointer = rdata + i*new_sec_size;
			//EMUtil::em_memset(begin_pointer, 0, (-y0)*new_nx*sizeof(float));
			float* end_pointer = begin_pointer+(-y0)*new_nx;
			std::fill(begin_pointer,end_pointer,fill_value);
		}

		// Set the extra back y components to the fill_value
		if ( civ.new_y_back > 0 )
		{
			float* begin_pointer = rdata + i*new_sec_size + (new_ny-civ.new_y_back)*new_nx;
			//EMUtil::em_memset(begin_pointer, 0, (civ.new_y_back)*new_nx*sizeof(float));
			float* end_pointer = begin_pointer+(civ.new_y_back)*new_nx;
			std::fill(begin_pointer,end_pointer,fill_value);
		}

		// Iterate through the y to set each correct x component to the fill_value
		for (int j = civ.new_y_front; j <civ.new_y_front + civ.y_iter; ++j)
		{
			// Set the extra left x components to the fill_value
			if ( x0 < 0 )
			{
				float* begin_pointer = rdata + i*new_sec_size + j*new_nx;
				//EMUtil::em_memset(begin_pointer, 0, (-x0)*sizeof(float));
				float* end_pointer = begin_pointer+(-x0);
				std::fill(begin_pointer,end_pointer,fill_value);
			}

			// Set the extra right x components to the fill_value
			if ( civ.new_x_right > 0 )
			{
				float* begin_pointer = rdata + i*new_sec_size + j*new_nx + (new_nx - civ.new_x_right);
				//EMUtil::em_memset(begin_pointer, 0, (civ.new_x_right)*sizeof(float));
				float* end_pointer = begin_pointer+(civ.new_x_right);
				std::fill(begin_pointer,end_pointer,fill_value);
			}

		}
	}

// These couts may be useful
// 	cout << "start starts " << civ.prv_x_left << " " << civ.prv_y_front << " " << civ.prv_z_bottom << endl;
// 	cout << "start ends " << civ.prv_x_right << " " << civ.prv_y_back << " " << civ.prv_z_top << endl;
// 	cout << "dst starts " << civ.new_x_left << " " << civ.new_y_front << " " << civ.new_z_bottom << endl;
// 	cout << "dst ends " << civ.new_x_right << " " << civ.new_y_back << " " << civ.new_z_top << endl;
// 	cout << "total iter z - " << civ.z_iter << " y - " << civ.y_iter << " x - " << civ.x_iter << endl;
// 	cout << "=====" << endl;
// 	cout << "dst_end is " << dst_it_end << " src end is " << src_it_end << endl;
// 	cout << "dst_begin is " << dst_it_begin << " src begin is " << src_it_begin << endl;

	// Update appropriate attributes (Copied and pasted from get_clip)
	if( attr_dict.has_key("origin_x") && attr_dict.has_key("origin_y") &&
	attr_dict.has_key("origin_z") )
	{
		float xorigin = attr_dict["origin_x"];
		float yorigin = attr_dict["origin_y"];
		float zorigin = attr_dict["origin_z"];

		float apix_x = attr_dict["apix_x"];
		float apix_y = attr_dict["apix_y"];
		float apix_z = attr_dict["apix_z"];

		set_xyz_origin(xorigin + apix_x * area.origin[0],
			yorigin + apix_y * area.origin[1],
			zorigin + apix_z * area.origin[2]);
	}

	// Set the update flag because the size of the image has changed and stats should probably be recalculated if requested.
	update();

	EXITFUNC;
}

EMData *EMData::get_clip(const Region & area, const float fill) const
{
	ENTERFUNC;
	if (get_ndim() != area.get_ndim()) {
		LOGERR("cannot get %dD clip out of %dD image", area.get_ndim(),get_ndim());
		return 0;
	}

	EMData *result = new EMData();

	// Ensure that all of the metadata of this is stored in the new object
	// Originally added to ensure that euler angles were retained when preprocessing (zero padding) images
	// prior to insertion into the 3D for volume in the reconstruction phase (see reconstructor.cpp/h).
	result->attr_dict = this->attr_dict;
	int zsize = (int)area.size[2];
	if (zsize == 0 && nz <= 1) {
		zsize = 1;
	}
	int ysize = (ny<=1 && (int)area.size[1]==0 ? 1 : (int)area.size[1]);

	if ( (int)area.size[0] < 0 || ysize < 0 || zsize < 0 )
	{
		// Negative image dimensions not supported - added retrospectively by d.woolford (who didn't write get_clip but wrote clip_inplace)
		throw ImageDimensionException("New image dimensions are negative - this is not supported in the the get_clip operation");
	}

#ifdef EMAN2_USING_CUDA
	// Strategy is always to prefer using the GPU if possible
	bool use_gpu = false;
	if ( gpu_operation_preferred() ) {
		result->set_size_cuda((int)area.size[0], ysize, zsize);
		//CudaDataLock lock(this); // Just so we never have to recopy this data to and from the GPU
		result->get_cuda_data(); // Force the allocation - set_size_cuda is lazy
		// Setting the value is necessary seeing as cuda data is not automatically zeroed
		result->to_value(fill); // This will automatically use the GPU.
		use_gpu = true;
	} else { // cpu == True
		result->set_size((int)area.size[0], ysize, zsize);
		if (fill != 0.0) { result->to_value(fill); };
	}
#else
	result->set_size((int)area.size[0], ysize, zsize);
	if (fill != 0.0) { result->to_value(fill); };
#endif //EMAN2_USING_CUDA

	int x0 = (int) area.origin[0];
	x0 = x0 < 0 ? 0 : x0;

	int y0 = (int) area.origin[1];
	y0 = y0 < 0 ? 0 : y0;

	int z0 = (int) area.origin[2];
	z0 = z0 < 0 ? 0 : z0;

	int x1 = (int) (area.origin[0] + area.size[0]);
	x1 = x1 > nx ? nx : x1;

	int y1 = (int) (area.origin[1] + area.size[1]);
	y1 = y1 > ny ? ny : y1;

	int z1 = (int) (area.origin[2] + area.size[2]);
	z1 = z1 > nz ? nz : z1;
	if (z1 <= 0) {
		z1 = 1;
	}

	result->insert_clip(this,-((IntPoint)area.origin));

	if( attr_dict.has_key("apix_x") && attr_dict.has_key("apix_y") &&
		attr_dict.has_key("apix_z") )
	{
		if( attr_dict.has_key("origin_x") && attr_dict.has_key("origin_y") &&
		    attr_dict.has_key("origin_z") )
		{
			float xorigin = attr_dict["origin_x"];
			float yorigin = attr_dict["origin_y"];
			float zorigin = attr_dict["origin_z"];

			float apix_x = attr_dict["apix_x"];
			float apix_y = attr_dict["apix_y"];
			float apix_z = attr_dict["apix_z"];

			result->set_xyz_origin(xorigin + apix_x * area.origin[0],
							   	   yorigin + apix_y * area.origin[1],
							       zorigin + apix_z * area.origin[2]);
		}
	}

#ifdef EMAN2_USING_CUDA
	if (use_gpu) result->gpu_update();
	else result->update();
#else
	result->update();
#endif // EMAN2_USING_CUDA


	result->set_path(path);
	result->set_pathnum(pathnum);

	EXITFUNC;
	return result;
}


EMData *EMData::get_top_half() const
{
	ENTERFUNC;

	if (get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}

	EMData *half = new EMData();
	half->attr_dict = attr_dict;
	half->set_size(nx, ny, nz / 2);

	float *half_data = half->get_data();
	EMUtil::em_memcpy(half_data, &(get_data()[nz / 2 * nx * ny]), sizeof(float) * nx * ny * nz / 2);

	float apix_z = attr_dict["apix_z"];
	float origin_z = attr_dict["origin_z"];
	origin_z += apix_z * nz / 2;
	half->attr_dict["origin_z"] = origin_z;
	half->update();

	EXITFUNC;
	return half;
}


EMData *EMData::get_rotated_clip(const Transform &xform,
								 const IntSize &size, float)
{
	EMData *result = new EMData();
	result->set_size(size[0],size[1],size[2]);

	if (nz==1) {
		for (int y=-size[1]/2; y<(size[1]+1)/2; y++) {
			for (int x=-size[0]/2; x<(size[0]+1)/2; x++) {
				Vec3f xv=xform.transform(Vec3f((float)x,(float)y,0.0f));
				float v = 0;

				if (xv[0]<0||xv[1]<0||xv[0]>nx-2||xv[1]>ny-2) v=0.;
				else v=sget_value_at_interp(xv[0],xv[1]);
				result->set_value_at(x+size[0]/2,y+size[1]/2,v);
			}
		}
	}
	else {
		for (int z=-size[2]/2; z<(size[2]+1)/2; z++) {
			for (int y=-size[1]/2; y<(size[1]+1)/2; y++) {
				for (int x=-size[0]/2; x<(size[0]+1)/2; x++) {
					Vec3f xv=xform.transform(Vec3f((float)x,(float)y,0.0f));
					float v = 0;

					if (xv[0]<0||xv[1]<0||xv[2]<0||xv[0]>nx-2||xv[1]>ny-2||xv[2]>nz-2) v=0.;
					else v=sget_value_at_interp(xv[0],xv[1],xv[2]);
					result->set_value_at(x+size[0]/2,y+size[1]/2,z+size[2]/2,v);
				}
			}
		}
	}
	result->update();

	return result;
}


EMData* EMData::window_center(int l) {
	ENTERFUNC;
	// sanity checks
	int n = nx;
	if (is_complex()) {
		LOGERR("Need real-space data for window_center()");
		throw ImageFormatException(
			"Complex input image; real-space expected.");
	}
	if (is_fftpadded()) {
		// image has been fft-padded, compute the real-space size
		n -= (2 - int(is_fftodd()));
	}
	int corner = n/2 - l/2;
	int ndim = get_ndim();
	EMData* ret;
	switch (ndim) {
		case 3:
			if ((n != ny) || (n != nz)) {
				LOGERR("Need the real-space image to be cubic.");
				throw ImageFormatException(
						"Need cubic real-space image.");
			}
			ret = get_clip(Region(corner, corner, corner, l, l, l));
			break;
		case 2:
			if (n != ny) {
				LOGERR("Need the real-space image to be square.");
				throw ImageFormatException(
						"Need square real-space image.");
			}
			//cout << "Using corner " << corner << endl;
			ret = get_clip(Region(corner, corner, l, l));
			break;
		case 1:
			ret = get_clip(Region(corner, l));
			break;
		default:
			throw ImageDimensionException(
					"window_center only supports 1-d, 2-d, and 3-d images");
	}
	return ret;
	EXITFUNC;
}


float *EMData::setup4slice(bool redo)
{
	ENTERFUNC;

	if (!is_complex()) {
		throw ImageFormatException("complex image only");
	}

	if (get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}

	if (supp) {
		if (redo) {
			EMUtil::em_free(supp);
			supp = 0;
		}
		else {
			EXITFUNC;
			return supp;
		}
	}

	const int SUPP_ROW_SIZE = 8;
	const int SUPP_ROW_OFFSET = 4;
	const int supp_size = SUPP_ROW_SIZE + SUPP_ROW_OFFSET;

	supp = (float *) EMUtil::em_calloc(supp_size * ny * nz, sizeof(float));
	int nxy = nx * ny;
	int supp_xy = supp_size * ny;
	float * data = get_data();

	for (int z = 0; z < nz; z++) {
		size_t cur_z1 = z * nxy;
		size_t cur_z2 = z * supp_xy;

		for (int y = 0; y < ny; y++) {
			size_t cur_y1 = y * nx;
			size_t cur_y2 = y * supp_size;

			for (int x = 0; x < SUPP_ROW_SIZE; x++) {
				size_t k = (x + SUPP_ROW_OFFSET) + cur_y2 + cur_z2;
				supp[k] = data[x + cur_y1 + cur_z1];
			}
		}
	}

	for (int z = 1, zz = nz - 1; z < nz; z++, zz--) {
		size_t cur_z1 = zz * nxy;
		size_t cur_z2 = z * supp_xy;

		for (int y = 1, yy = ny - 1; y < ny; y++, yy--) {
			supp[y * 12 + cur_z2] = data[4 + yy * nx + cur_z1];
			supp[1 + y * 12 + cur_z2] = -data[5 + yy * nx + cur_z1];
			supp[2 + y * 12 + cur_z2] = data[2 + yy * nx + cur_z1];
			supp[3 + y * 12 + cur_z2] = -data[3 + yy * nx + cur_z1];
		}
	}

	EXITFUNC;
	return supp;
}


void EMData::scale(float s)
{
	ENTERFUNC;
	Transform t;
	t.set_scale(s);
	transform(t);
	EXITFUNC;
}


void EMData::translate(int dx, int dy, int dz)
{
	ENTERFUNC;
	translate(Vec3i(dx, dy, dz));
	EXITFUNC;
}


void EMData::translate(float dx, float dy, float dz)
{
	ENTERFUNC;
	int dx_ = Util::round(dx);
	int dy_ = Util::round(dy);
	int dz_ = Util::round(dz);
	if( ( (dx-dx_) == 0 ) && ( (dy-dy_) == 0 ) && ( (dz-dz_) == 0 )) {
		translate(dx_, dy_, dz_);
	}
	else {
		translate(Vec3f(dx, dy, dz));
	}
	EXITFUNC;
}


void EMData::translate(const Vec3i &translation)
{
	ENTERFUNC;

	//if traslation is 0, do nothing
	if( translation[0] == 0 && translation[1] == 0 && translation[2] == 0) {
		EXITFUNC;
		return;
	}

	Dict params("trans",static_cast< vector<int> >(translation));
	process_inplace("math.translate.int",params);

	// update() - clip_inplace does the update
	all_translation += translation;

	EXITFUNC;
}


void EMData::translate(const Vec3f &translation)
{
	ENTERFUNC;

	if( translation[0] == 0.0f && translation[1] == 0.0f && translation[2] == 0.0f ) {
		EXITFUNC;
		return;
	}

	Transform* t = new Transform();
	t->set_trans(translation);
	process_inplace("xform",Dict("transform",t));
	delete t;

	all_translation += translation;
	EXITFUNC;
}


void EMData::rotate(float az, float alt, float phi)
{
	Dict d("type","eman");
	d["az"] = az;
	d["alt"] = alt;
	d["phi"] = phi;
	Transform t(d);
	transform(t);
}



void EMData::rotate(const Transform3D & t)
{
	cout << "Deprecation warning in EMData::rotate. Please consider using EMData::transform() instead " << endl;
	rotate_translate(t);
}

float EMData::max_3D_pixel_error(const Transform &t1, const Transform & t2, float r) {
	
	Transform t;
	int r0 = (int)r;
	float ddmax = 0.0f;

	t = t2*t1.inverse();
	for (int i=0; i<int(2*M_PI*r0+0.5); i++) {
		Vec3f v = Vec3f(r0*cos((float)i), r0*sin((float)i), 0);
		Vec3f d = t*v-v;
		float dd = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
		if (dd > ddmax) ddmax = dd; 
	}
	return std::sqrt(ddmax);
}

void EMData::rotate_translate(float az, float alt, float phi, float dx, float dy, float dz)
{
	cout << "Deprecation warning in EMData::rotate_translate. Please consider using EMData::transform() instead " << endl;
	Transform3D t( az, alt, phi,Vec3f(dx, dy, dz));
	rotate_translate(t);
}


void EMData::rotate_translate(float az, float alt, float phi, float dx, float dy,
							  float dz, float pdx, float pdy, float pdz)
{
	cout << "Deprecation warning in EMData::rotate_translate. Please consider using EMData::transform() instead " << endl;
	Transform3D t(Vec3f(dx, dy, dz), az, alt, phi, Vec3f(pdx,pdy,pdz));
	rotate_translate(t);
}

void EMData::rotate_translate(const Transform3D & RA)
{
	cout << "Deprecation warning in EMData::rotate_translate. Please consider using EMData::transform() instead " << endl;
	ENTERFUNC;

#if EMDATA_EMAN2_DEBUG
	std::cout << "start rotate_translate..." << std::endl;
#endif

	float scale       = RA.get_scale();
 	Vec3f dcenter     = RA.get_center();
	Vec3f translation = RA.get_posttrans();
	Dict rotation      = RA.get_rotation(Transform3D::EMAN);
//	Transform3D mx = RA;
	Transform3D RAInv = RA.inverse(); // We're rotating the coordinate system, not the data
// 	RAInv.printme();
#if EMDATA_EMAN2_DEBUG
	vector<string> keys = rotation.keys();
	vector<string>::const_iterator it;
	for(it=keys.begin(); it!=keys.end(); ++it) {
//		std::cout << *it << " : " << rotation[*it] << std::endl;
		std::cout << *it << " : " << (float)rotation.get(*it) << std::endl;
	}
#endif
	float inv_scale = 1.0f;

	if (scale != 0) {
		inv_scale = 1.0f / scale;
	}

	float *src_data = 0;
	float *des_data = 0;

	src_data = get_data();
	des_data = (float *) EMUtil::em_malloc(nx * ny * nz * sizeof(float));

	if (nz == 1) {
		float x2c =  nx / 2 - dcenter[0] + RAInv[0][3];
		float y2c =  ny / 2 - dcenter[1] + RAInv[1][3];
		float y   = -ny / 2 + dcenter[1]; // changed 0 to 1 in dcenter and below
		for (int j = 0; j < ny; j++, y += 1.0f) {
			float x = -nx / 2 + dcenter[0];
			for (int i = 0; i < nx; i++, x += 1.0f) {
				float x2 = RAInv[0][0]*x + RAInv[0][1]*y + x2c;
				float y2 = RAInv[1][0]*x + RAInv[1][1]*y + y2c;

				if (x2 < 0 || x2 >= nx || y2 < 0 || y2 >= ny ) {
					des_data[i + j * nx] = 0; // It may be tempting to set this value to the
					// mean but in fact this is not a good thing to do. Talk to S.Ludtke about it.
				}
				else {
					int ii = Util::fast_floor(x2);
					int jj = Util::fast_floor(y2);
					int k0 = ii + jj * nx;
					int k1 = k0 + 1;
					int k2 = k0 + nx;
					int k3 = k0 + nx + 1;

					if (ii == nx - 1) {
						k1--;
						k3--;
					}
					if (jj == ny - 1) {
						k2 -= nx;
						k3 -= nx;
					}

					float t = x2 - ii;
					float u = y2 - jj;

					des_data[i + j * nx] = Util::bilinear_interpolate(src_data[k0],src_data[k1], src_data[k2], src_data[k3],t,u); // This is essentially basic interpolation
				}
			}
		}
	}
	else {
#if EMDATA_EMAN2_DEBUG
		std::cout << "This is the 3D case." << std::endl    ;
#endif

		Transform3D mx = RA;
		mx.set_scale(inv_scale);

#if EMDATA_EMAN2_DEBUG
//		std::cout << v4[0] << " " << v4[1]<< " " << v4[2]<< " "
//			<< dcenter2[0]<< " "<< dcenter2[1]<< " "<< dcenter2[2] << std::endl;
#endif

		int nxy = nx * ny;
		int l = 0;

		float x2c =  nx / 2 - dcenter[0] + RAInv[0][3];;
		float y2c =  ny / 2 - dcenter[1] + RAInv[1][3];;
		float z2c =  nz / 2 - dcenter[2] + RAInv[2][3];;

		float z   = -nz / 2 + dcenter[2]; //

		size_t ii, k0, k1, k2, k3, k4, k5, k6, k7;
		for (int k = 0; k < nz; k++, z += 1.0f) {
			float y   = -ny / 2 + dcenter[1]; //
			for (int j = 0; j < ny; j++,   y += 1.0f) {
				float x = -nx / 2 + dcenter[0];
				for (int i = 0; i < nx; i++, l++ ,  x += 1.0f) {
					float x2 = RAInv[0][0] * x + RAInv[0][1] * y + RAInv[0][2] * z + x2c;
					float y2 = RAInv[1][0] * x + RAInv[1][1] * y + RAInv[1][2] * z + y2c;
					float z2 = RAInv[2][0] * x + RAInv[2][1] * y + RAInv[2][2] * z + z2c;


					if (x2 < 0 || y2 < 0 || z2 < 0 ||
						x2 >= nx  || y2 >= ny  || z2>= nz ) {
						des_data[l] = 0;
					}
					else {
						int ix = Util::fast_floor(x2);
						int iy = Util::fast_floor(y2);
						int iz = Util::fast_floor(z2);
						float tuvx = x2-ix;
						float tuvy = y2-iy;
						float tuvz = z2-iz;
						ii = ix + iy * nx + iz * nxy;

						k0 = ii;
						k1 = k0 + 1;
						k2 = k0 + nx;
						k3 = k0 + nx+1;
						k4 = k0 + nxy;
						k5 = k1 + nxy;
						k6 = k2 + nxy;
						k7 = k3 + nxy;

						if (ix == nx - 1) {
							k1--;
							k3--;
							k5--;
							k7--;
						}
						if (iy == ny - 1) {
							k2 -= nx;
							k3 -= nx;
							k6 -= nx;
							k7 -= nx;
						}
						if (iz == nz - 1) {
							k4 -= nxy;
							k5 -= nxy;
							k6 -= nxy;
							k7 -= nxy;
						}

						des_data[l] = Util::trilinear_interpolate(src_data[k0],
							  src_data[k1],
							  src_data[k2],
							  src_data[k3],
							  src_data[k4],
							  src_data[k5],
							  src_data[k6],
							  src_data[k7],
							  tuvx, tuvy, tuvz);
#if EMDATA_EMAN2_DEBUG
						printf(" ix=%d \t iy=%d \t iz=%d \t value=%f \n", ix ,iy, iz, des_data[l] );
						std::cout << src_data[ii] << std::endl;
#endif
					}
				}
			}
		}
	}

	if( rdata )
	{
		EMUtil::em_free(rdata);
		rdata = 0;
	}
	rdata = des_data;

	scale_pixel(inv_scale);

	attr_dict["origin_x"] = (float) attr_dict["origin_x"] * inv_scale;
	attr_dict["origin_y"] = (float) attr_dict["origin_y"] * inv_scale;
	attr_dict["origin_z"] = (float) attr_dict["origin_z"] * inv_scale;

	update();
	all_translation += translation;
	EXITFUNC;
}




void EMData::rotate_x(int dx)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("no 3D image");
	}


	size_t row_size = nx * sizeof(float);
	float *tmp = (float*)EMUtil::em_malloc(row_size);
	float * data = get_data();

	for (int y = 0; y < ny; y++) {
		int y_nx = y * nx;
		for (int x = 0; x < nx; x++) {
			tmp[x] = data[y_nx + (x + dx) % nx];
		}
		EMUtil::em_memcpy(&data[y_nx], tmp, row_size);
	}

	update();
	if( tmp )
	{
		delete[]tmp;
		tmp = 0;
	}
	EXITFUNC;
}

double EMData::dot_rotate_translate(EMData * with, float dx, float dy, float da, const bool mirror)
{
	ENTERFUNC;

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size");
		throw ImageFormatException("images not same size");
	}

	if (get_ndim() == 3) {
		LOGERR("1D/2D Images only");
		throw ImageDimensionException("1D/2D only");
	}

	float *this_data = 0;

	this_data = get_data();

	float da_rad = da*(float)M_PI/180.0f;

	float *with_data = with->get_data();
	float mx0 = cos(da_rad);
	float mx1 = sin(da_rad);
	float y = -ny / 2.0f;
	float my0 = mx0 * (-nx / 2.0f - 1.0f) + nx / 2.0f - dx;
	float my1 = -mx1 * (-nx / 2.0f - 1.0f) + ny / 2.0f - dy;
	double result = 0;

	for (int j = 0; j < ny; j++) {
		float x2 = my0 + mx1 * y;
		float y2 = my1 + mx0 * y;

		int ii = Util::fast_floor(x2);
		int jj = Util::fast_floor(y2);
		float t = x2 - ii;
		float u = y2 - jj;

		for (int i = 0; i < nx; i++) {
			t += mx0;
			u -= mx1;

			if (t >= 1.0f) {
				ii++;
				t -= 1.0f;
			}

			if (u >= 1.0f) {
				jj++;
				u -= 1.0f;
			}

			if (t < 0) {
				ii--;
				t += 1.0f;
			}

			if (u < 0) {
				jj--;
				u += 1.0f;
			}

			if (ii >= 0 && ii <= nx - 2 && jj >= 0 && jj <= ny - 2) {
				int k0 = ii + jj * nx;
				int k1 = k0 + 1;
				int k2 = k0 + nx + 1;
				int k3 = k0 + nx;

				float tt = 1 - t;
				float uu = 1 - u;
				int idx = i + j * nx;
				if (mirror) idx = nx-1-i+j*nx; // mirroring of Transforms is always about the y axis
				result += (this_data[k0] * tt * uu + this_data[k1] * t * uu +
						   this_data[k2] * t * u + this_data[k3] * tt * u) * with_data[idx];
			}
		}
		y += 1.0f;
	}

	EXITFUNC;
	return result;
}


EMData *EMData::little_big_dot(EMData * with, bool do_sigma)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D only");
	}

	EMData *ret = copy_head();
	ret->set_size(nx,ny,nz);
	ret->to_zero();

	int nx2 = with->get_xsize();
	int ny2 = with->get_ysize();
	float em = with->get_edge_mean();

	float *data = get_data();
	float *with_data = with->get_data();
	float *ret_data = ret->get_data();

	float sum2 = (Util::square((float)with->get_attr("sigma")) +
				  Util::square((float)with->get_attr("mean")));
	if (do_sigma) {
		for (int j = ny2 / 2; j < ny - ny2 / 2; j++) {
			for (int i = nx2 / 2; i < nx - nx2 / 2; i++) {
				float sum = 0;
				float sum1 = 0;
				float summ = 0;
				int k = 0;

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
						int l = ii + jj * nx;
						sum1 += Util::square(data[l]);
						summ += data[l];
						sum += data[l] * with_data[k];
						k++;
					}
				}
				float tmp_f1 = (sum1 / 2.0f - sum) / (nx2 * ny2);
				float tmp_f2 = Util::square((float)with->get_attr("mean") -
											summ / (nx2 * ny2));
				ret_data[i + j * nx] = sum2 + tmp_f1 - tmp_f2;
			}
		}
	}
	else {
		for (int j = ny2 / 2; j < ny - ny2 / 2; j++) {
			for (int i = nx2 / 2; i < nx - nx2 / 2; i++) {
				float eml = 0;
				float dot = 0;
				float dot2 = 0;

				for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
					eml += data[ii + (j - ny2 / 2) * nx] + data[ii + (j + ny2 / 2 - 1) * nx];
				}

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					eml += data[i - nx2 / 2 + jj * nx] + data[i + nx2 / 2 - 1 + jj * nx];
				}

				eml /= (nx2 + ny2) * 2.0f;
				int k = 0;

				for (int jj = j - ny2 / 2; jj < j + ny2 / 2; jj++) {
					for (int ii = i - nx2 / 2; ii < i + nx2 / 2; ii++) {
						dot += (data[ii + jj * nx] - eml) * (with_data[k] - em);
						dot2 += Util::square(data[ii + jj * nx] - eml);
						k++;
					}
				}

				dot2 = std::sqrt(dot2);

				if (dot2 == 0) {
					ret_data[i + j * nx] = 0;
				}
				else {
					ret_data[i + j * nx] = dot / (nx2 * ny2 * dot2 * (float)with->get_attr("sigma"));
				}
			}
		}
	}

	ret->update();

	EXITFUNC;
	return ret;
}


EMData *EMData::do_radon()
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

	if (nx != ny) {
		throw ImageFormatException("square image only");
	}

	EMData *result = new EMData();
	result->set_size(nx, ny, 1);
	result->to_zero();
	float *result_data = result->get_data();

	EMData *this_copy = this;
	this_copy = copy();

	for (int i = 0; i < nx; i++) {
		Transform t(Dict("type","2d","alpha",(float) M_PI * 2.0f * i / nx));
		this_copy->transform(t);

		float *copy_data = this_copy->get_data();

		for (int y = 0; y < nx; y++) {
			for (int x = 0; x < nx; x++) {
				if (Util::square(x - nx / 2) + Util::square(y - nx / 2) <= nx * nx / 4) {
					result_data[i + y * nx] += copy_data[x + y * nx];
				}
			}
		}

		this_copy->update();
	}

	result->update();

	if( this_copy )
	{
		delete this_copy;
		this_copy = 0;
	}

	EXITFUNC;
	return result;
}

void EMData::zero_corner_circulant(const int radius)
{
	if ( nz > 1 && nz < (2*radius+1) ) throw ImageDimensionException("Error: cannot zero corner - nz is too small");
	if ( ny > 1 && ny < (2*radius+1) ) throw ImageDimensionException("Error: cannot zero corner - ny is too small");
	if ( nx > 1 && nx < (2*radius+1) ) throw ImageDimensionException("Error: cannot zero corner - nx is too small");

	int it_z = radius;
	int it_y = radius;
	int it_x = radius;

	if ( nz == 1 ) it_z = 0;
	if ( ny == 1 ) it_y = 0;
	if ( nx == 1 ) it_z = 0;

	if ( nz == 1 && ny == 1 )
	{
		for ( int x = -it_x; x <= it_x; ++x )
			get_value_at_wrap(x) = 0;

	}
	else if ( nz == 1 )
	{
		for ( int y = -it_y; y <= it_y; ++y)
			for ( int x = -it_x; x <= it_x; ++x )
				get_value_at_wrap(x,y) = 0;
	}
	else
	{
		for( int z = -it_z; z <= it_z; ++z )
			for ( int y = -it_y; y <= it_y; ++y)
				for ( int x = -it_x; x < it_x; ++x )
					get_value_at_wrap(x,y,z) = 0;

	}

}

EMData *EMData::calc_ccf(EMData * with, fp_flag fpflag,bool center)
{
	ENTERFUNC;
	if( with == 0 ) {
		EXITFUNC;
		return convolution(this,this,fpflag, center);
	}
	else if ( with == this ){ // this if statement is not necessary, the correlation function tests to see if with == this
		EXITFUNC;
		return correlation(this, this, fpflag,center);
	}
	else {

#ifdef EMAN2_USING_CUDA
		if (gpu_operation_preferred()) {
			EXITFUNC;
			return calc_ccf_cuda(with,false,center);
		}
#endif

		// If the argument EMData pointer is not the same size we automatically resize it
		bool undoresize = false;
		int wnx = with->get_xsize(); int wny = with->get_ysize(); int wnz = with->get_zsize();
		if ( wnx != nx || wny != ny || wnz != nz ) {
			Region r((wnx-nx)/2, (wny-ny)/2, (wnz-nz)/2,nx,ny,nz);
			with->clip_inplace(r);
			undoresize = true;
		}

		EMData* cor = correlation(this, with, fpflag, center);

		// If the argument EMData pointer was resized, it is returned to its original dimensions
		if ( undoresize ) {
			Region r((nx-wnx)/2, (ny-wny)/2,(nz-wnz)/2,wnx,wny,wnz);
			with->clip_inplace(r);
		}

		EXITFUNC;
		return cor;
	}
}

EMData *EMData::calc_ccfx( EMData * const with, int y0, int y1, bool no_sum)
{
	ENTERFUNC;

#ifdef EMAN2_USING_CUDA
	if (gpu_operation_preferred() ) {
		EXITFUNC;
		return calc_ccfx_cuda(with,y0,y1,no_sum);
	}
#endif

	if (!with) {
		LOGERR("NULL 'with' image. ");
		throw NullPointerException("NULL input image");
	}

	if (!EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size: (%d,%d,%d) != (%d,%d,%d)",
			   nx, ny, nz,
			   with->get_xsize(), with->get_ysize(), with->get_zsize());
		throw ImageFormatException("images not same size");
	}
	if (get_ndim() > 2) {
		LOGERR("2D images only");
		throw ImageDimensionException("2D images only");
	}

	if (y1 <= y0) {
		y1 = ny;
	}

	if (y0 >= y1) {
		y0 = 0;
	}

	if (y0 < 0) {
		y0 = 0;
	}

	if (y1 > ny) {
		y1 = ny;
	}
	if (is_complex_x() || with->is_complex_x() ) throw; // Woops don't support this anymore!

	static int nx_fft = 0;
	static int ny_fft = 0;
	static EMData f1;
	static EMData f2;
	static EMData rslt;

	int height = y1-y0;
	int width = (nx+2-(nx%2));
	if (width != nx_fft || height != ny_fft ) {
		f1.set_size(width,height);
		f2.set_size(width,height);
		rslt.set_size(nx,height);
		nx_fft = width;
		ny_fft = height;
	}

	float *d1 = get_data();
	float *d2 = with->get_data();
	float *f1d = f1.get_data();
	float *f2d = f2.get_data();
	for (int j = 0; j < height; j++) {
		EMfft::real_to_complex_1d(d1 + j * nx, f1d+j*width, nx);
		EMfft::real_to_complex_1d(d2 + j * nx, f2d+j*width, nx);
	}

	for (int j = 0; j < height; j++) {
		float *f1a = f1d + j * width;
		float *f2a = f2d + j * width;

		for (int i = 0; i < width / 2; i++) {
			float re1 = f1a[2*i];
			float re2 = f2a[2*i];
			float im1 = f1a[2*i+1];
			float im2 = f2a[2*i+1];

			f1d[j*width+i*2] = re1 * re2 + im1 * im2;
			f1d[j*width+i*2+1] = im1 * re2 - re1 * im2;
		}
	}

	float* rd = rslt.get_data();
	for (int j = y0; j < y1; j++) {
		EMfft::complex_to_real_1d(f1d+j*width, rd+j*nx, nx);
	}

	if (no_sum) {
		rslt.update(); // This is important in terms of the copy - the returned object won't have the correct flags unless we do this
		EXITFUNC;
		return new EMData(rslt);
	} else {
		EMData *cf = new EMData(nx,1,1);
		cf->to_zero();
		float *c = cf->get_data();
		for (int j = 0; j < height; j++) {
			for(int i = 0; i < nx; ++i) {
				c[i] += rd[i+j*nx];
			}
		}
		cf->update();
		EXITFUNC;
		return cf;
	}
}

EMData *EMData::make_rotational_footprint_cmc( bool unwrap) {
	ENTERFUNC;
	update_stat();
	// Note that rotational_footprint caching saves a large amount of time
	// but this is at the expense of memory. Note that a policy is hardcoded here,
	// that is that caching is only employed when premasked is false and unwrap
	// is true - this is probably going to be what is used in most scenarios
	// as advised by Steve Ludtke - In terms of performance this caching doubles the metric
	// generated by e2speedtest.
	if ( rot_fp != 0 && unwrap == true) {
		return new EMData(*rot_fp);
	}

	static EMData obj_filt;
	EMData* filt = &obj_filt;
	filt->set_complex(true);


	// The filter object is nothing more than a cached high pass filter
	// Ultimately it is used an argument to the EMData::mult(EMData,prevent_complex_multiplication (bool))
	// function in calc_mutual_correlation. Note that in the function the prevent_complex_multiplication
	// set to true, which is used for speed reasons.
	if (filt->get_xsize() != nx+2-(nx%2) || filt->get_ysize() != ny ||
		   filt->get_zsize() != nz ) {
		filt->set_size(nx+2-(nx%2), ny, nz);
		filt->to_one();

		filt->process_inplace("eman1.filter.highpass.gaussian", Dict("highpass", 1.5f/nx));
	}

	EMData *ccf = this->calc_mutual_correlation(this, true,filt);
	ccf->sub(ccf->get_edge_mean());
	EMData *result = ccf->unwrap();
	delete ccf; ccf = 0;

	EXITFUNC;
	if ( unwrap == true)
	{
	// this if statement reflects a strict policy of caching in only one scenario see comments at beginning of function block

// Note that the if statement at the beginning of this function ensures that rot_fp is not zero, so there is no need
// to throw any exception
// if ( rot_fp != 0 ) throw UnexpectedBehaviorException("The rotational foot print is only expected to be cached if it is not NULL");

// Here is where the caching occurs - the rot_fp takes ownsherhip of the pointer, and a deep copied EMData object is returned.
// The deep copy invokes a cost in terms of CPU cycles and memory, but prevents the need for complicated memory management (reference counting)
		rot_fp = result;
		return new EMData(*rot_fp);
	}
	else return result;
}

EMData *EMData::make_rotational_footprint( bool unwrap) {
	ENTERFUNC;
	update_stat();
	// Note that rotational_footprint caching saves a large amount of time
	// but this is at the expense of memory. Note that a policy is hardcoded here,
	// that is that caching is only employed when premasked is false and unwrap
	// is true - this is probably going to be what is used in most scenarios
	// as advised by Steve Ludtke - In terms of performance this caching doubles the metric
	// generated by e2speedtest.
	if ( rot_fp != 0 && unwrap == true) {
		return new EMData(*rot_fp);
	}

	EMData* ccf = this->calc_ccf(this,CIRCULANT,true);
	ccf->sub(ccf->get_edge_mean());
	//ccf->process_inplace("xform.phaseorigin.tocenter"); ccf did the centering
	EMData *result = ccf->unwrap();
	delete ccf; ccf = 0;

	EXITFUNC;
	if ( unwrap == true)
	{ // this if statement reflects a strict policy of caching in only one scenario see comments at beginning of function block

// Note that the if statement at the beginning of this function ensures that rot_fp is not zero, so there is no need
// to throw any exception
// if ( rot_fp != 0 ) throw UnexpectedBehaviorException("The rotational foot print is only expected to be cached if it is not NULL");

// Here is where the caching occurs - the rot_fp takes ownsherhip of the pointer, and a deep copied EMData object is returned.
// The deep copy invokes a cost in terms of CPU cycles and memory, but prevents the need for complicated memory management (reference counting)
		rot_fp = result;
		return new EMData(*rot_fp);
	}
	else return result;
}

EMData *EMData::make_rotational_footprint_e1( bool unwrap)
{
	ENTERFUNC;
#ifdef EMAN2_USING_CUDA
	if (gpu_operation_preferred()) {
		EXITFUNC;
		return make_rotational_footprint_cuda(unwrap);
	}
#endif

	update_stat();
	// Note that rotational_footprint caching saves a large amount of time
	// but this is at the expense of memory. Note that a policy is hardcoded here,
	// that is that caching is only employed when premasked is false and unwrap
	// is true - this is probably going to be what is used in most scenarios
	// as advised by Steve Ludtke - In terms of performance this caching doubles the metric
	// generated by e2speedtest.
	if ( rot_fp != 0 && unwrap == true) {
		return new EMData(*rot_fp);
	}

	static EMData obj_filt;
	EMData* filt = &obj_filt;
	filt->set_complex(true);
// 	Region filt_region;

// 	if (nx & 1) {
// 		LOGERR("even image xsize only");		throw ImageFormatException("even image xsize only");
// 	}

	int cs = (((nx * 7 / 4) & 0xfffff8) - nx) / 2; // this pads the image to 1 3/4 * size with result divis. by 8

	static EMData big_clip;
	int big_x = nx+2*cs;
	int big_y = ny+2*cs;
	int big_z = 1;
	if ( nz != 1 ) {
		big_z = nz+2*cs;
	}


	if ( big_clip.get_xsize() != big_x || big_clip.get_ysize() != big_y || big_clip.get_zsize() != big_z ) {
		big_clip.set_size(big_x,big_y,big_z);
	}
	// It is important to set all newly established pixels around the boundaries to the mean
	// If this is not done then the associated rotational alignment routine breaks, in fact
	// everythin just goes foo.
	big_clip.to_value(get_edge_mean());

	if (nz != 1) {
		big_clip.insert_clip(this,IntPoint(cs,cs,cs));
	} else  {
		big_clip.insert_clip(this,IntPoint(cs,cs,0));
	}


	// The filter object is nothing more than a cached high pass filter
	// Ultimately it is used an argument to the EMData::mult(EMData,prevent_complex_multiplication (bool))
	// function in calc_mutual_correlation. Note that in the function the prevent_complex_multiplication
	// set to true, which is used for speed reasons.
	if (filt->get_xsize() != big_clip.get_xsize() +2-(big_clip.get_xsize()%2) || filt->get_ysize() != big_clip.get_ysize() ||
		   filt->get_zsize() != big_clip.get_zsize()) {
		filt->set_size(big_clip.get_xsize() + 2-(big_clip.get_xsize()%2), big_clip.get_ysize(), big_clip.get_zsize());
	filt->to_one();
	filt->process_inplace("eman1.filter.highpass.gaussian", Dict("highpass", 1.5f/nx));
	}
	EMData *mc = big_clip.calc_mutual_correlation(&big_clip, true,filt);
 	mc->sub(mc->get_edge_mean());

	static EMData sml_clip;
	int sml_x = nx * 3 / 2;
	int sml_y = ny * 3 / 2;
	int sml_z = 1;
	if ( nz != 1 ) {
		sml_z = nz * 3 / 2;
	}

	if ( sml_clip.get_xsize() != sml_x || sml_clip.get_ysize() != sml_y || sml_clip.get_zsize() != sml_z ) {
		sml_clip.set_size(sml_x,sml_y,sml_z);	}
	if (nz != 1) {
		sml_clip.insert_clip(mc,IntPoint(-cs+nx/4,-cs+ny/4,-cs+nz/4));
	} else {
		sml_clip.insert_clip(mc,IntPoint(-cs+nx/4,-cs+ny/4,0));
	}

	delete mc; mc = 0;
	EMData * result = NULL;

	if (nz == 1) {
		if (!unwrap) {
			result = sml_clip.process("mask.sharp", Dict("outer_radius", -1, "value", 0));

		}
		else {
			result = sml_clip.unwrap();
		}
	}
	else {
		// I am not sure why there is any consideration of non 2D images, but it was here
		// in the first port so I kept when I cleaned this function up (d.woolford)
// 		result = clipped_mc;
		result = new EMData(sml_clip);
	}

	EXITFUNC;
	if ( unwrap == true)
	{ // this if statement reflects a strict policy of caching in only one scenario see comments at beginning of function block

		// Note that the if statement at the beginning of this function ensures that rot_fp is not zero, so there is no need
		// to throw any exception
		if ( rot_fp != 0 ) throw UnexpectedBehaviorException("The rotational foot print is only expected to be cached if it is not NULL");

		// Here is where the caching occurs - the rot_fp takes ownsherhip of the pointer, and a deep copied EMData object is returned.
		// The deep copy invokes a cost in terms of CPU cycles and memory, but prevents the need for complicated memory management (reference counting)
		rot_fp = result;
		return new EMData(*rot_fp);
	}
	else return result;
}

EMData *EMData::make_footprint(int type)
{
//	printf("Make fp %d\n",type);
	if (type==0) {
		EMData *un=make_rotational_footprint_e1(); // Use EMAN1's footprint strategy
		if (un->get_ysize() <= 6) {
			throw UnexpectedBehaviorException("In EMData::make_footprint. The rotational footprint is too small");
		}
		EMData *tmp=un->get_clip(Region(0,4,un->get_xsize(),un->get_ysize()-6));	// 4 and 6 are empirical
		EMData *cx=tmp->calc_ccfx(tmp,0,-1,1);
		EMData *fp=cx->get_clip(Region(0,0,cx->get_xsize()/2,cx->get_ysize()));
		delete un;
		delete tmp;
		delete cx;
		return fp;
	}
	else if (type==1 || type==2 ||type==5 || type==6) {
		int i,j,kx,ky,lx,ly;

		EMData *fft=do_fft();

		// map for x,y -> radius for speed
		int rmax=(get_xsize()+1)/2;
		float *rmap=(float *)malloc(rmax*rmax*sizeof(float));
		for (i=0; i<rmax; i++) {
			for (j=0; j<rmax; j++) {
#ifdef _WIN32
				rmap[i+j*rmax]=_hypotf((float)i,(float)j);
#else
				rmap[i+j*rmax]=hypot((float)i,(float)j);
#endif	//_WIN32
//				printf("%d\t%d\t%f\n",i,j,rmap[i+j*rmax]);
			}
		}

		EMData *fp=new EMData(rmax*2+2,rmax*2,1);
		fp->set_complex(1);
		fp->to_zero();

		// Two vectors in to complex space (kx,ky) and (lx,ly)
		// We are computing the bispectrum, f(k).f(l).f*(k+l)
		// but integrating out two dimensions, leaving |k|,|l|
		for (kx=-rmax+1; kx<rmax; kx++) {
			for (ky=-rmax+1; ky<rmax; ky++) {
				for (lx=-rmax+1; lx<rmax; lx++) {
					for (ly=-rmax+1; ly<rmax; ly++) {
						int ax=kx+lx;
						int ay=ky+ly;
						if (abs(ax)>=rmax || abs(ay)>=rmax) continue;
						int r1=(int)floor(.5+rmap[abs(kx)+rmax*abs(ky)]);
						int r2=(int)floor(.5+rmap[abs(lx)+rmax*abs(ly)]);
//						if (r1>500 ||r2>500) printf("%d\t%d\t%d\t%d\t%d\t%d\n",kx,ky,lx,ly,r1,r2);
//						float r3=rmap[ax+rmax*ay];
						if (r1+r2>=rmax) continue;

						std::complex<float> p=fft->get_complex_at(kx,ky)*fft->get_complex_at(lx,ly)*conj(fft->get_complex_at(ax,ay));
						fp->set_value_at(r1*2,r2,p.real()+fp->get_value_at(r1*2,r2));		// We keep only the real component in anticipation of zero phase sum
//						fp->set_value_at(r1*2,rmax*2-r2-1,  fp->get_value_at(r1*2,r2));		// We keep only the real component in anticipation of zero phase sum
//						fp->set_value_at(r1*2+1,r2,p.real()+fp->get_value_at(r1*2+1,r2));		// We keep only the real component in anticipation of zero phase sum
						fp->set_value_at(r1*2+1,r2,fp->get_value_at(r1*2+1,r2)+1);			// a normalization counter
					}
				}
			}
		}

		// Normalizes the pixels based on the accumulated counts then sets the imaginary components back to zero
		if (type==5 || type==6) {
			for (i=0; i<rmax*2; i+=2) {
				for (j=0; j<rmax; j++) {
					float norm=fp->get_value_at(i+1,j);
#ifdef _WIN32
					fp->set_value_at(i,rmax*2-j-1,pow(fp->get_value_at(i,j)/(norm==0.0f?1.0f:norm), 1.0f/3.0f));
					fp->set_value_at(i,j,pow(fp->get_value_at(i,j)/(norm==0.0f?1.0f:norm), 1.0f/3.0f));
#else
					fp->set_value_at(i,rmax*2-j-1,cbrt(fp->get_value_at(i,j)/(norm==0?1.0:norm)));
					fp->set_value_at(i,j,cbrt(fp->get_value_at(i,j)/(norm==0?1.0:norm)));
#endif	//_WIN32
					fp->set_value_at(i+1,j,0.0);
				}
			}
		}
		else {
			for (i=0; i<rmax*2; i+=2) {
				for (j=0; j<rmax; j++) {
					float norm=fp->get_value_at(i+1,j);
					fp->set_value_at(i,rmax*2-j-1,fp->get_value_at(i,j)/(norm==0.0f?1.0f:norm));
					fp->set_value_at(i,j,fp->get_value_at(i,j)/(norm==0.0f?1.0f:norm));
					fp->set_value_at(i+1,j,0.0);
				}
			}
		}

		free(rmap);
		if (type==2||type==6) {
			EMData *f2=fp->do_ift();
			if (f2->get_value_at(0,0)<0) f2->mult(-1.0f);
			f2->process_inplace("xform.phaseorigin.tocorner");
			delete fp;
			return f2;
		}
		return fp;
	}
	else if (type==3 || type==4) {
		int h,i,j,kx,ky,lx,ly;

		EMData *fft=do_fft();

		// map for x,y -> radius for speed
		int rmax=(get_xsize()+1)/2;
		float *rmap=(float *)malloc(rmax*rmax*sizeof(float));
		for (i=0; i<rmax; i++) {
			for (j=0; j<rmax; j++) {
#ifdef _WIN32
				rmap[i+j*rmax]=_hypotf((float)i,(float)j);
#else
				rmap[i+j*rmax]=hypot((float)i,(float)j);
#endif	//_WIN32
//				printf("%d\t%d\t%f\n",i,j,rmap[i+j*rmax]);
			}
		}

		EMData *fp=new EMData(rmax*2+2,rmax*2,16);

		fp->set_complex(1);
		fp->to_zero();

		// Two vectors in to complex space (kx,ky) and (lx,ly)
		// We are computing the bispectrum, f(k).f(l).f*(k+l)
		// but integrating out two dimensions, leaving |k|,|l|
		for (kx=-rmax+1; kx<rmax; kx++) {
			for (ky=-rmax+1; ky<rmax; ky++) {
				for (lx=-rmax+1; lx<rmax; lx++) {
					for (ly=-rmax+1; ly<rmax; ly++) {
						int ax=kx+lx;
						int ay=ky+ly;
						if (abs(ax)>=rmax || abs(ay)>=rmax) continue;
						float rr1=rmap[abs(kx)+rmax*abs(ky)];
						float rr2=rmap[abs(lx)+rmax*abs(ly)];
						int r1=(int)floor(.5+rr1);
						int r2=(int)floor(.5+rr2);
//						if (r1>500 ||r2>500) printf("%d\t%d\t%d\t%d\t%d\t%d\n",kx,ky,lx,ly,r1,r2);
//						float r3=rmap[ax+rmax*ay];
						if (r1+r2>=rmax || rr1==0 ||rr2==0) continue;

						std::complex<float> p=fft->get_complex_at(kx,ky)*fft->get_complex_at(lx,ly)*conj(fft->get_complex_at(ax,ay));
						int dot=(int)floor((kx*lx+ky*ly)/(rr1*rr2)*7.5);					// projection of k on l 0-31
						if (dot<0) dot=16+dot;
//						int dot=(int)floor((kx*lx+ky*ly)/(rr1*rr2)*7.5+8.0);					// projection of k on l 0-15
						fp->set_value_at(r1*2,r2,dot,p.real()+fp->get_value_at(r1*2,r2,dot));		// We keep only the real component in anticipation of zero phase sum
//						fp->set_value_at(r1*2,rmax*2-r2-1,  fp->get_value_at(r1*2,r2));		// We keep only the real component in anticipation of zero phase sum
//						fp->set_value_at(r1*2+1,r2,p.real()+fp->get_value_at(r1*2+1,r2));		// We keep only the real component in anticipation of zero phase sum
						fp->set_value_at(r1*2+1,r2,dot,fp->get_value_at(r1*2+1,r2,dot)+1);			// a normalization counter
					}
				}
			}
		}

		// Normalizes the pixels based on the accumulated counts then sets the imaginary components back to zero
		for (i=0; i<rmax*2; i+=2) {
			for (j=0; j<rmax; j++) {
				for (h=0; h<16; h++) {
					float norm=fp->get_value_at(i+1,j,h);
//					fp->set_value_at(i,rmax*2-j-1,h,cbrt(fp->get_value_at(i,j,h)/(norm==0?1.0:norm)));
//					fp->set_value_at(i,j,h,cbrt(fp->get_value_at(i,j,h)/(norm==0?1.0:norm)));
					fp->set_value_at(i,rmax*2-j-1,h,(fp->get_value_at(i,j,h)/(norm==0.0f?1.0f:norm)));
					fp->set_value_at(i,j,h,(fp->get_value_at(i,j,h)/(norm==0.0f?1.0f:norm)));
	//				fp->set_value_at(i,rmax*2-j-1,fp->get_value_at(i,j)/(norm==0?1.0:norm));
	//				fp->set_value_at(i,j,fp->get_value_at(i,j)/(norm==0?1.0:norm));
					fp->set_value_at(i+1,j,h,0.0);
				}
			}
		}

		free(rmap);
		if (type==4) {
			EMData *f2=fp->do_ift();
			if (f2->get_value_at(0,0,0)<0) f2->mult(-1.0f);
			f2->process_inplace("xform.phaseorigin.tocorner");
			delete fp;
			return f2;
		}
		return fp;
	}
	throw UnexpectedBehaviorException("There is not implementation for the parameters you specified");
}


EMData *EMData::calc_mutual_correlation(EMData * with, bool tocenter, EMData * filter)
{
	ENTERFUNC;

	if (with && !EMUtil::is_same_size(this, with)) {
		LOGERR("images not same size");
		throw ImageFormatException( "images not same size");
	}

	EMData *this_fft = 0;
	this_fft = do_fft();

	if (!this_fft) {

		LOGERR("FFT returns NULL image");
		throw NullPointerException("FFT returns NULL image");
	}

	this_fft->ap2ri();
	EMData *cf = 0;

	if (with && with != this) {
		cf = with->do_fft();
		if (!cf) {
			LOGERR("FFT returns NULL image");
			throw NullPointerException("FFT returns NULL image");
		}
		cf->ap2ri();
	}
	else {
		cf = this_fft->copy();
	}

	if (filter) {
		if (!EMUtil::is_same_size(filter, cf)) {
			LOGERR("improperly sized filter");
			throw ImageFormatException("improperly sized filter");
		}

		cf->mult_complex_efficient(*filter,true);
		this_fft->mult(*filter,true);
		/*cf->mult_complex_efficient(*filter,5);
		this_fft->mult_complex_efficient(*filter,5);*/
	}

	float *rdata1 = this_fft->get_data();
	float *rdata2 = cf->get_data();
	size_t this_fft_size = this_fft->get_xsize() * this_fft->get_ysize() * this_fft->get_zsize();

	if (with == this) {
		for (size_t i = 0; i < this_fft_size; i += 2) {
			rdata2[i] = std::sqrt(rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1]);
			rdata2[i + 1] = 0;
		}

		this_fft->update();
		cf->update();
	}
	else {
		for (size_t i = 0; i < this_fft_size; i += 2) {
			rdata2[i] = (rdata1[i] * rdata2[i] + rdata1[i + 1] * rdata2[i + 1]);
			rdata2[i + 1] = (rdata1[i + 1] * rdata2[i] - rdata1[i] * rdata2[i + 1]);
		}

		this_fft->update();
		cf->update();
		rdata1 = cf->get_data();

		for (size_t i = 0; i < this_fft_size; i += 2) {
			float t = Util::square(rdata1[i]) + Util::square(rdata1[i + 1]);
			if (t != 0) {
				t = pow(t, (float) 0.25);
				rdata1[i] /= t;
				rdata1[i + 1] /= t;
			}
		}
		cf->update();
	}

	EMData *f2 = cf->do_ift();

	if (tocenter) {
		f2->process_inplace("xform.phaseorigin.tocenter");
	}

	if( cf )
	{
		delete cf;
		cf = 0;
	}

	if( this_fft )
	{
		delete this_fft;
		this_fft = 0;
	}

	f2->set_attr("label", "MCF");
	f2->set_path("/tmp/eman.mcf");

	EXITFUNC;
	return f2;
}


vector < float > EMData::calc_hist(int hist_size, float histmin, float histmax,const float& brt, const float& cont)
{
	ENTERFUNC;

	static size_t prime[] = { 1, 3, 7, 11, 17, 23, 37, 59, 127, 253, 511 };

	if (histmin == histmax) {
		histmin = get_attr("minimum");
		histmax = get_attr("maximum");
	}

	vector <float> hist(hist_size, 0.0);

	int p0 = 0;
	int p1 = 0;
	size_t size = nx * ny * nz;
	if (size < 300000) {
		p0 = 0;
		p1 = 0;
	}
	else if (size < 2000000) {
		p0 = 2;
		p1 = 3;
	}
	else if (size < 8000000) {
		p0 = 4;
		p1 = 6;
	}
	else {
		p0 = 7;
		p1 = 9;
	}

	if (is_complex() && p0 > 0) {
		p0++;
		p1++;
	}

	size_t di = 0;
//	float norm = 0;
	size_t n = hist.size();

	float * data = get_data();
	for (int k = p0; k <= p1; ++k) {
		if (is_complex()) {
			di = prime[k] * 2;
		}
		else {
			di = prime[k];
		}

//		norm += (float)size / (float) di;
		float w = (float)n / (histmax - histmin);

		for(size_t i=0; i<=size-di; i += di) {
			float val;
			if (cont != 1.0f || brt != 0)val = cont*(data[i]+brt);
			else val = data[i];
			int j = Util::round((val - histmin) * w);
			if (j >= 0 && j < (int) n) {
				hist[j] += 1;
			}
		}
	}
/*
	for (size_t i = 0; i < hist.size(); ++i) {
		if (norm != 0) {
			hist[i] = hist[i] / norm;
		}
	}
*/
	return hist;

	EXITFUNC;
}





vector<float> EMData::calc_az_dist(int n, float a0, float da, float rmin, float rmax)
{
	ENTERFUNC;

	a0=a0*M_PI/180.0f;
	da=da*M_PI/180.0f;

	if (get_ndim() > 2) {
		throw ImageDimensionException("no 3D image");
	}

	float *yc = new float[n];

	vector<float>	vd(n);
	for (int i = 0; i < n; i++) {
		yc[i] = 0.00001f;
	}

	float * data = get_data();
	if (is_complex()) {
		int c = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x += 2, c += 2) {
				int x1 = x / 2;
				int y1 = y<ny/2?y:y-ny;
				float r = (float)Util::hypot_fast(x1,y1);

				if (r >= rmin && r <= rmax) {
					float a = 0;

					if (y != ny / 2 || x != 0) {
						a = (atan2((float)y1, (float)x1) - a0) / da;
					}

					int i = (int)(floor(a));
					a -= i;

					if (i == 0) {
						vd[0] += data[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						vd[n - 1] += data[c] * a;
						yc[n - 1] += a;
					}
					else if (i > 0 && i < (n - 1)) {
						float h = 0;
						if (is_ri()) {
#ifdef	_WIN32
							h = (float)_hypot(data[c], data[c + 1]);
#else
							h = (float)hypot(data[c], data[c + 1]);
#endif	//_WIN32
						}
						else {
							h = data[c];
						}

						vd[i] += h * (1.0f - a);
						yc[i] += (1.0f - a);
						vd[i + 1] += h * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}
	else {
		int c = 0;
		float half_nx = (nx - 1) / 2.0f;
		float half_ny = (ny - 1) / 2.0f;

		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++, c++) {
				float y1 = y - half_ny;
				float x1 = x - half_nx;
#ifdef	_WIN32
				float r = (float)_hypot(x1, y1);
#else
				float r = (float)hypot(x1, y1);
#endif

				if (r >= rmin && r <= rmax) {
					float a = 0;
					if (x1 != 0 || y1 != 0) {
						a = atan2(y1, x1);
						if (a < 0) {
							a += M_PI * 2;
						}
					}

					a = (a - a0) / da;
					int i = static_cast < int >(floor(a));
					a -= i;

					if (i == 0) {
						vd[0] += data[c] * (1.0f - a);
						yc[0] += (1.0f - a);
					}
					else if (i == n - 1) {
						vd[n - 1] += data[c] * a;
						yc[n - 1] += (a);
					}
					else if (i > 0 && i < (n - 1)) {
						vd[i] += data[c] * (1.0f - a);
						yc[i] += (1.0f - a);
						vd[i + 1] += data[c] * a;
						yc[i + 1] += a;
					}
				}
			}
		}
	}


	for (int i = 0; i < n; i++) {
		vd[i] /= yc[i];
	}

	if( yc )
	{
		delete[]yc;
		yc = 0;
	}

	return vd;

	EXITFUNC;
}


EMData *EMData::unwrap(int r1, int r2, int xs, int dx, int dy, bool do360, bool weight_radial) const
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}

	int p = 1;
	if (do360) {
		p = 2;
	}

	if (xs < 1) {
		xs = (int) Util::fast_floor(p * M_PI * ny / 4);
		xs -= xs % 8;
		if (xs<=8) xs=16;
	}

	if (r1 < 0) {
		r1 = 4;
	}

#ifdef	_WIN32
	int rr = ny / 2 - 2 - (int) Util::fast_floor(static_cast<float>(_hypot(dx, dy)));
#else
	int rr = ny / 2 - 2 - (int) Util::fast_floor(static_cast<float>(hypot(dx, dy)));
#endif	//_WIN32
	rr-=rr%2;
	if (r2 <= r1 || r2 > rr) {
		r2 = rr;
	}

	if ( (r2-r1) < 0 ) throw UnexpectedBehaviorException("The combination of function the arguments and the image dimensions causes unexpected behavior internally. Use a larger image, or a smaller value of r1, or a combination of both");

#ifdef EMAN2_USING_CUDA
	if ( gpu_operation_preferred() ) {
// 		cout << "Binding " << cuda_cache_handle << endl;
		bind_cuda_texture();
		EMData* rslt = new EMData();
		rslt->set_size_cuda(xs,r2-r1,1);
		EMDataForCuda r = rslt->get_data_struct_for_cuda();
// 		CudaDataLock lock1(rslt);
		/*EMDataForCuda* tmp = */emdata_unwrap(&r,r1,r2,xs,p,dx,dy,weight_radial,nx,ny);
		unbind_cuda_texture();
// 		EMData* e = new EMData();
// 		e->set_gpu_rw_data(tmp->data,tmp->nx,tmp->ny,tmp->nz);
// 		free(tmp);
		return  rslt;
	}
#endif

	EMData *ret = new EMData();
	ret->set_size(xs, r2 - r1, 1);
	const float *const d = get_const_data();
	float *dd = ret->get_data();
	float pfac = (float)p/(float)xs;

	int nxon2 = nx/2;
	int nyon2 = ny/2;
	for (int x = 0; x < xs; x++) {
		float ang = x * M_PI * pfac;
		float si = sin(ang);
		float co = cos(ang);

		for (int y = 0; y < r2 - r1; y++) {
			float ypr1 = (float)y + r1;
			float xx = ypr1 * co + nxon2 + dx;
			float yy = ypr1 * si + nyon2 + dy;
//			float t = xx - Util::fast_floor(xx);
//			float u = yy - Util::fast_floor(yy);
			float t = xx - (int)xx;
			float u = yy - (int)yy;
//			int k = (int) Util::fast_floor(xx) + (int) (Util::fast_floor(yy)) * nx;
			int k = (int) xx + ((int) yy) * nx;
			float val = Util::bilinear_interpolate(d[k], d[k + 1], d[k + nx], d[k + nx+1], t,u);
			if (weight_radial) val *=  ypr1;
			dd[x + y * xs] = val;
		}
	}
	ret->update();
	EXITFUNC;
	return ret;
}

// NOTE : x axis is from 0 to 0.5  (Nyquist), and thus properly handles non-square images
// complex only
void EMData::apply_radial_func(float x0, float step, vector < float >array, bool interp)
{
	ENTERFUNC;

	if (!is_complex()) throw ImageFormatException("apply_radial_func requires a complex image");

	int n = static_cast < int >(array.size());

	if (n*step>2.0) printf("Warning, apply_radial_func takes x0 and step with respect to Nyquist (0.5)\n");

//	printf("%f %f %f\n",array[0],array[25],array[50]);

	ap2ri();

	size_t ndims = get_ndim();
	float * data = get_data();
	if (ndims == 2) {
		int k = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i += 2, k += 2) {
				float r;
#ifdef	_WIN32
				if (j<ny/2) r = (float)_hypot(i/(float)(nx*2), j/(float)ny);
				else r = (float)_hypot(i/(float)(nx*2), (ny-j)/(float)ny);
#else
				if (j<ny/2) r = (float)hypot(i/(float)(nx*2), j/(float)ny);
				else r = (float)hypot(i/(float)(nx*2), (ny-j)/(float)ny);
#endif	//_WIN32
				r = (r - x0) / step;

				int l = 0;
				if (interp) {
					l = (int) floor(r);
				}
				else {
					l = (int) floor(r + 1);
				}


				float f = 0;
				if (l >= n - 2) {
					f = array[n - 1];
				}
				else {
					if (interp) {
						r -= l;
						f = (array[l] * (1.0f - r) + array[l + 1] * r);
					}
					else {
						f = array[l];
					}
				}

				data[k] *= f;
				data[k + 1] *= f;
			}
		}
	}
	else if (ndims == 3) {
		int k = 0;
		for (int m = 0; m < nz; m++) {
			float mnz;
			if (m<nz/2) mnz=m*m/(float)(nz*nz);
			else { mnz=(nz-m)/(float)nz; mnz*=mnz; }

			for (int j = 0; j < ny; j++) {
				float jny;
				if (j<ny/2) jny= j*j/(float)(ny*ny);
				else { jny=(ny-j)/(float)ny; jny*=jny; }

				for (int i = 0; i < nx; i += 2, k += 2) {
					float r = std::sqrt((i * i / (nx*nx*4.0f)) + jny + mnz);
					r = (r - x0) / step;

					int l = 0;
					if (interp) {
						l = (int) floor(r);
					}
					else {
						l = (int) floor(r + 1);
					}


					float f = 0;
					if (l >= n - 2) {
						f = array[n - 1];
					}
					else {
						if (interp) {
							r -= l;
							f = (array[l] * (1.0f - r) + array[l + 1] * r);
						}
						else {
							f = array[l];
						}
					}

					data[k] *= f;
					data[k + 1] *= f;
				}
			}
		}

	}

	update();
	EXITFUNC;
}

vector<float> EMData::calc_radial_dist(int n, float x0, float dx, bool inten)
{
	ENTERFUNC;

	vector<float>ret(n);
	vector<float>norm(n);

	int x,y,z,i;
	int step=is_complex()?2:1;
	int isinten=get_attr_default("is_intensity",0);

	if (isinten&&!inten) { throw InvalidParameterException("Must set inten for calc_radial_dist with intensity image"); }

	for (i=0; i<n; i++) ret[i]=norm[i]=0.0;
	float * data = get_data();

	// We do 2D separately to avoid the hypot3 call
	if (nz==1) {
		for (y=i=0; y<ny; y++) {
			for (x=0; x<nx; x+=step,i+=step) {
				float r,v;
				if (step==2) {		//complex
					if (x==0 && y>ny/2) continue;
					r=(float)(Util::hypot_fast(x/2,y<ny/2?y:ny-y));		// origin at 0,0; periodic
					if (!inten) {
#ifdef	_WIN32
						if (is_ri()) v=static_cast<float>(_hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#else
						if (is_ri()) v=static_cast<float>(hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#endif
						else v=data[i];							// amp/phase, just get amp
					} else {
						if (isinten) v=data[i];
						else if (is_ri()) v=data[i]*data[i]+data[i+1]*data[i+1];
						else v=data[i]*data[i];
					}
				}
				else {
					r=(float)(Util::hypot_fast(x-nx/2,y-ny/2));
					if (inten) v=data[i]*data[i];
					else v=data[i];
				}
				r=(r-x0)/dx;
				int f=int(r);	// safe truncation, so floor isn't needed
				r-=float(f);	// r is now the fractional spacing between bins
//				printf("%d\t%d\t%d\t%1.3f\t%d\t%1.3f\t%1.4g\n",x,y,f,r,step,Util::hypot_fast(x/2,y<ny/2?y:ny-y),v);
				if (f>=0 && f<n) {
					ret[f]+=v*(1.0f-r);
					norm[f]+=(1.0f-r);
					if (f<n-1) {
						ret[f+1]+=v*r;
						norm[f+1]+=r;
					}
				}
			}
		}
	}
	else {
		size_t i;	//3D file may have >2G size
		for (z=i=0; z<nz; ++z) {
			for (y=0; y<ny; ++y) {
				for (x=0; x<nx; x+=step,i+=step) {
					float r,v;
					if (step==2) {	//complex
						if (x==0 && z<nz/2) continue;
						if (x==0 && z==nz/2 && y<ny/2) continue;
						r=Util::hypot3(x/2,y<ny/2?y:ny-y,z<nz/2?z:nz-z);	// origin at 0,0; periodic
						if (!inten) {
#ifdef	_WIN32
							if (is_ri()) v=static_cast<float>(_hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#else
							if (is_ri()) v=static_cast<float>(hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#endif	//_WIN32
							else v=data[i];							// amp/phase, just get amp
						} else {
							if (isinten) v=data[i];
							else if (is_ri()) v=data[i]*data[i]+data[i+1]*data[i+1];
							else v=data[i]*data[i];
						}
					}
					else {
						r=Util::hypot3(x-nx/2,y-ny/2,z-nz/2);
						if (inten) v=data[i]*data[i];
						else v=data[i];
					}
					r=(r-x0)/dx;
					int f=int(r);	// safe truncation, so floor isn't needed
					r-=float(f);	// r is now the fractional spacing between bins
					if (f>=0 && f<n) {
						ret[f]+=v*(1.0f-r);
						norm[f]+=(1.0f-r);
						if (f<n-1) {
							ret[f+1]+=v*r;
							norm[f+1]+=r;
						}
					}
				}
			}
		}
	}

	for (i=0; i<n; i++) ret[i]/=norm[i]?norm[i]:1.0f;	// Normalize

	EXITFUNC;

	return ret;
}

vector<float> EMData::calc_radial_dist(int n, float x0, float dx, int nwedge, bool inten)
{
	ENTERFUNC;

	if (nz > 1) {
		LOGERR("2D images only.");
		throw ImageDimensionException("2D images only");
	}

	vector<float>ret(n*nwedge);
	vector<float>norm(n*nwedge);

	int x,y,i;
	int step=is_complex()?2:1;
	float astep=static_cast<float>(M_PI*2.0/nwedge);
	float* data = get_data();
	for (i=0; i<n*nwedge; i++) ret[i]=norm[i]=0.0;

	// We do 2D separately to avoid the hypot3 call
	for (y=i=0; y<ny; y++) {
		for (x=0; x<nx; x+=step,i+=step) {
			float r,v,a;
			if (is_complex()) {
#ifdef	_WIN32
				r=static_cast<float>(_hypot(x/2.0,y<ny/2?y:ny-y));		// origin at 0,0; periodic
#else
				r=static_cast<float>(hypot(x/2.0,y<ny/2?y:ny-y));		// origin at 0,0; periodic
#endif
				a=atan2(float(y<ny/2?y:ny-y),x/2.0f);
				if (!inten) {
#ifdef	_WIN32
					if (is_ri()) v=static_cast<float>(_hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#else
					if (is_ri()) v=static_cast<float>(hypot(data[i],data[i+1]));	// real/imag, compute amplitude
#endif	//_WIN32
					else v=data[i];							// amp/phase, just get amp
				} else {
					if (is_ri()) v=data[i]*data[i]+data[i+1]*data[i+1];
					else v=data[i]*data[i];
				}
			}
			else {
#ifdef	_WIN32
				r=static_cast<float>(_hypot(x-nx/2,y-ny/2));
#else
				r=static_cast<float>(hypot(x-nx/2,y-ny/2));
#endif	//_WIN32
				a=atan2(float(y-ny/2),float(x-nx/2));
				if (inten) v=data[i]*data[i];
				else v=data[i];
			}
			int bin=n*int((a+M_PI)/astep);
			if (bin>=nwedge) bin=nwedge-1;
			r=(r-x0)/dx;
			int f=int(r);	// safe truncation, so floor isn't needed
			r-=float(f);	// r is now the fractional spacing between bins
			if (f>=0 && f<n) {
				ret[f+bin]+=v*(1.0f-r);
				norm[f+bin]+=(1.0f-r);
				if (f<n-1) {
					ret[f+1+bin]+=v*r;
					norm[f+1+bin]+=r;
				}
			}
		}
	}

	for (i=0; i<n*nwedge; i++) ret[i]/=norm[i]?norm[i]:1.0f;	// Normalize
	EXITFUNC;

	return ret;
}

void EMData::cconj() {
	ENTERFUNC;
	if (!is_complex() || !is_ri())
		throw ImageFormatException("EMData::conj requires a complex, ri image");
	int nxreal = nx -2 + int(is_fftodd());
	int nxhalf = nxreal/2;
	for (int iz = 0; iz < nz; iz++)
		for (int iy = 0; iy < ny; iy++)
			for (int ix = 0; ix <= nxhalf; ix++)
				cmplx(ix,iy,iz) = conj(cmplx(ix,iy,iz));
	EXITFUNC;
}


void EMData::update_stat() const
{
	ENTERFUNC;
//	printf("update stat %f %d\n",(float)attr_dict["mean"],flags);
	if (!(flags & EMDATA_NEEDUPD))
	{
		EXITFUNC;
		return;
	}

	float* data = get_data();
	float max = -FLT_MAX;
	float min = -max;

	double sum = 0;
	double square_sum = 0;

	int step = 1;
	if (is_complex() && !is_ri()) {
		step = 2;
	}

	int n_nonzero = 0;

	//cout << "point 1" << endl;
	//cout << "size is " << nx << " " << ny << " " << nz << endl;

	size_t size = nx*ny*nz;
	for (size_t i = 0; i < size; i += step) {
		float v = data[i];
	#ifdef _WIN32
		max = _cpp_max(max,v);
		min = _cpp_min(min,v);
	#else
		max=std::max<float>(max,v);
		min=std::min<float>(min,v);
	#endif	//_WIN32
		sum += v;
		square_sum += v * (double)(v);
		if (v != 0) n_nonzero++;
	}
	//cout << "Point 2" << endl;
	size_t n = size / step;
	double mean = sum / n;

#ifdef _WIN32
	float sigma = (float)std::sqrt( _cpp_max(0.0,(square_sum - sum*sum / n)/(n-1)));
	n_nonzero = _cpp_max(1,n_nonzero);
	double sigma_nonzero = std::sqrt( _cpp_max(0,(square_sum  - sum*sum/n_nonzero)/(n_nonzero-1)));
#else
	float sigma = (float)std::sqrt(std::max<double>(0.0,(square_sum - sum*sum / n)/(n-1)));
	n_nonzero = std::max<int>(1,n_nonzero);
	double sigma_nonzero = std::sqrt(std::max<double>(0,(square_sum  - sum*sum/n_nonzero)/(n_nonzero-1)));
#endif	//_WIN32
	double mean_nonzero = sum / n_nonzero; // previous version overcounted! G2

	attr_dict["minimum"] = min;
	attr_dict["maximum"] = max;
	attr_dict["mean"] = (float)(mean);
	attr_dict["sigma"] = (float)(sigma);
	attr_dict["square_sum"] = (float)(square_sum);
	attr_dict["mean_nonzero"] = (float)(mean_nonzero);
	attr_dict["sigma_nonzero"] = (float)(sigma_nonzero);
	attr_dict["is_complex"] = (int) is_complex();
	attr_dict["is_complex_ri"] = (int) is_ri();

	flags &= ~EMDATA_NEEDUPD;

	if (rot_fp != 0)
	{
		delete rot_fp; rot_fp = 0;
	}

	EXITFUNC;
//	printf("done stat %f %f %f\n",(float)mean,(float)max,(float)sigma);
}

bool EMData::operator==(const EMData& that) const {
	if (that.get_xsize() != nx || that.get_ysize() != ny || that.get_zsize() != nz ) return false;

	const float*  d1 = that.get_const_data();
	float* d2 = get_data();

	for(size_t i =0; i < get_size(); ++i,++d1,++d2) {
		if ((*d1) != (*d2)) return false;
	}
	return true;

}

EMData * EMAN::operator+(const EMData & em, float n)
{
	EMData * r = em.copy();
	r->add(n);
	return r;
}

EMData * EMAN::operator-(const EMData & em, float n)
{
	EMData* r = em.copy();
	r->sub(n);
	return r;
}

EMData * EMAN::operator*(const EMData & em, float n)
{
	EMData* r = em.copy();
	r ->mult(n);
	return r;
}

EMData * EMAN::operator/(const EMData & em, float n)
{
	EMData * r = em.copy();
	r->div(n);
	return r;
}


EMData * EMAN::operator+(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->add(n);
	return r;
}

EMData * EMAN::operator-(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->mult(-1.0f);
	r->add(n);
	return r;
}

EMData * EMAN::operator*(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->mult(n);
	return r;
}

EMData * EMAN::operator/(float n, const EMData & em)
{
	EMData * r = em.copy();
	r->to_one();
	r->mult(n);
	r->div(em);

	return r;
}

EMData * EMAN::rsub(const EMData & em, float n)
{
	return EMAN::operator-(n, em);
}

EMData * EMAN::rdiv(const EMData & em, float n)
{
	return EMAN::operator/(n, em);
}

EMData * EMAN::operator+(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->add(b);
	return r;
}

EMData * EMAN::operator-(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->sub(b);
	return r;
}

EMData * EMAN::operator*(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->mult(b);
	return r;
}

EMData * EMAN::operator/(const EMData & a, const EMData & b)
{
	EMData * r = a.copy();
	r->div(b);
	return r;
}

void EMData::set_xyz_origin(float origin_x, float origin_y, float origin_z)
{
	attr_dict["origin_x"] = origin_x;
	attr_dict["origin_y"] = origin_y;
	attr_dict["origin_z"] = origin_z;
}

#if 0
void EMData::calc_rcf(EMData * with, vector < float >&sum_array)
{
	ENTERFUNC;

	int array_size = sum_array.size();
	float da = 2 * M_PI / array_size;
	float *dat = new float[array_size + 2];
	float *dat2 = new float[array_size + 2];
	int nx2 = nx * 9 / 20;

	float lim = 0;
	if (fabs(translation[0]) < fabs(translation[1])) {
		lim = fabs(translation[1]);
	}
	else {
		lim = fabs(translation[0]);
	}

	nx2 -= static_cast < int >(floor(lim));

	for (int i = 0; i < array_size; i++) {
		sum_array[i] = 0;
	}

	float sigma = attr_dict["sigma"];
	float with_sigma = with->get_attr_dict().get("sigma");

	vector<float> vdata, vdata2;
	for (int i = 8; i < nx2; i += 6) {
		vdata = calc_az_dist(array_size, 0, da, i, i + 6);
		vdata2 = with->calc_az_dist(array_size, 0, da, i, i + 6);
		Assert(vdata.size() <= array_size + 2);
		Assert(cdata2.size() <= array_size + 2);
		std::copy(vdata.begin(), vdata.end(), dat);
		std::copy(vdata2.begin(), vdata2.end(), dat2);

		EMfft::real_to_complex_1d(dat, dat, array_size);
		EMfft::real_to_complex_1d(dat2, dat2, array_size);

		for (int j = 0; j < array_size + 2; j += 2) {
			float max = dat[j] * dat2[j] + dat[j + 1] * dat2[j + 1];
			float max2 = dat[j + 1] * dat2[j] - dat2[j + 1] * dat[j];
			dat[j] = max;
			dat[j + 1] = max2;
		}

		EMfft::complex_to_real_1d(dat, dat, array_size);
		float norm = array_size * array_size * (4.0f * sigma) * (4.0f * with_sigma);

		for (int j = 0; j < array_size; j++) {
			sum_array[j] += dat[j] * (float) i / norm;
		}
	}

	if( dat )
	{
		delete[]dat;
		dat = 0;
	}

	if( dat2 )
	{
		delete[]dat2;
		dat2 = 0;
	}
	EXITFUNC;
}

#endif

void EMData::add_incoherent(EMData * obj)
{
	ENTERFUNC;

	if (!obj) {
		LOGERR("NULL image");
		throw NullPointerException("NULL image");
	}

	if (!obj->is_complex() || !is_complex()) {
		throw ImageFormatException("complex images only");
	}

	if (!EMUtil::is_same_size(this, obj)) {
		throw ImageFormatException("images not same size");
	}

	ri2ap();
	obj->ri2ap();

	float *dest = get_data();
	float *src = obj->get_data();
	size_t size = nx * ny * nz;
	for (size_t j = 0; j < size; j += 2) {
#ifdef	_WIN32
		dest[j] = (float) _hypot(src[j], dest[j]);
#else
		dest[j] = (float) hypot(src[j], dest[j]);
#endif	//_WIN32
		dest[j + 1] = 0;
	}

	obj->update();
	update();
	EXITFUNC;
}


float EMData::calc_dist(EMData * second_img, int y_index) const
{
	ENTERFUNC;

	if (get_ndim() != 1) {
		throw ImageDimensionException("'this' image is 1D only");
	}

	if (second_img->get_xsize() != nx || ny != 1) {
		throw ImageFormatException("image xsize not same");
	}

	if (y_index > second_img->get_ysize() || y_index < 0) {
		return -1;
	}

	float ret = 0;
	float *d1 = get_data();
	float *d2 = second_img->get_data() + second_img->get_xsize() * y_index;

	for (int i = 0; i < nx; i++) {
		ret += Util::square(d1[i] - d2[i]);
	}
	EXITFUNC;
	return std::sqrt(ret);
}


EMData * EMData::calc_fast_sigma_image( EMData* mask)
{
	ENTERFUNC;

	bool maskflag = false;
	if (mask == 0) {
		mask = new EMData(nx,ny,nz);
		mask->process_inplace("testimage.circlesphere");
		maskflag = true;
	}

	if (get_ndim() != mask->get_ndim() ) throw ImageDimensionException("The dimensions do not match");

	int mnx = mask->get_xsize(); int mny = mask->get_ysize(); int mnz = mask->get_zsize();

	if ( mnx > nx || mny > ny || mnz > nz)
		throw ImageDimensionException("Can not calculate variance map using an image that is larger than this image");


//	bool undoclip = false;

	int nxc = nx+mnx; int nyc = ny+mny; int nzc = nz+mnz;
//	if ( mnx < nx || mny < ny || mnz < nz) {
	Region r;
	if (ny == 1) r = Region((mnx-nxc)/2,nxc);
	else if (nz == 1) r = Region((mnx-nxc)/2, (mny-nyc)/2,nxc,nyc);
	else r = Region((mnx-nxc)/2, (mny-nyc)/2,(mnz-nzc)/2,nxc,nyc,nzc);
	mask->clip_inplace(r,1.0); //ming change 0.0 to 1.0

// ming move the following from above
	size_t P = 0;
	for(size_t i = 0; i < mask->get_size(); ++i){
		if (mask->get_value_at(i) != 0){
			++P;
		}
	}

	float normfac = 1.0f/(float)P;


	//Region r((mnx-nxc)/2, (mny-nyc)/2,(mnz-nzc)/2,nxc,nyc,nzc);
	//mask->clip_inplace(r);
	//undoclip = true;
	//}

	// Here we generate the local average of the squares
	Region r2;
	if (ny == 1) r2 = Region((nx-nxc)/2,nxc);
	else if (nz == 1) r2 = Region((nx-nxc)/2, (ny-nyc)/2,nxc,nyc);
	else r2 = Region((nx-nxc)/2, (ny-nyc)/2,(nz-nzc)/2,nxc,nyc,nzc);

	EMData* squared = get_clip(r2,get_edge_mean());
	EMData* tmp = squared->copy();
	Dict pow;
	pow["pow"] = 2.0f;
	squared->process_inplace("math.pow",pow);
	EMData* s = squared->convolute(mask);
	squared->mult(normfac);

	EMData* m = tmp->convolute(mask);
	m->mult(normfac);
	m->process_inplace("math.pow",pow);
	delete tmp; tmp = 0;
	s->sub(*m);
	// Here we finally generate the standard deviation image
	s->process_inplace("math.sqrt");
//	if ( undoclip ) {
//		Region r((nx-mnx)/2, (ny-mny)/2, (nz-mnz)/2,mnx,mny,mnz);
//		mask->clip_inplace(r);
//	}
	if (maskflag) {
		delete mask;
		mask = 0;
	} else {
		Region r;
		if (ny == 1) r = Region((nxc-mnx)/2,mnx);
		else if (nz == 1) r = Region((nxc-mnx)/2, (nyc-mny)/2,mnx,mny);
		else r = Region((nxc-mnx)/2, (nyc-mny)/2,(nzc-mnz)/2,mnx,mny,mnz);
		mask->clip_inplace(r);
	}
	delete squared;
	delete m;
	s->process_inplace("xform.phaseorigin.tocenter");
	Region r3;
	if (ny == 1) r3 = Region((nxc-nx)/2,nx);
	else if (nz == 1) r3 = Region((nxc-nx)/2, (nyc-ny)/2,nx,ny);
	else r3 = Region((nxc-nx)/2, (nyc-ny)/2,(nzc-nz)/2,nx,ny,nz);
	s->clip_inplace(r3);
	//s->process_inplace("xform.phaseorigin.tocorner");//ming add
	EXITFUNC;
	return s;
}

/*
EMData *EMData::calcFLCF(EMData *with, int r, int type) {

  EMData *f1, *f2, *sqf1, *ccf, *conv1, *conv2, *f, *fones,*mask, *v, *lcf;
  float NM;
  int nx, ny, nz, i, j, k;
  f1=this->copy();//target map //ming
  f2=with->copy();//search object //ming
  nx=f1->get_xsize();
  ny=f1->get_ysize();
  nz=f1->get_zsize();
  float N=(nx*ny*nz);
/////////////////setting min to zero in images/////////////////
  float f1min,f2min;
  f1min=f1->Min();
  f1->add(-f1min); //ming change addConst to add
  f2min=f2->Min();
  f2->add(-f2min);
////////////////////////////////////////////////////////////////////
/////////////////////making a binary mask with radius r////////////////
  fones=f1->copy();//ming, fones is M
  fones->one();
  fones->applyMask(r,type);//
  //ming, fones->toCorner();
  fones->process_inplace("xform.phaseorigin.tocorner");
////////////////////////////////////////////////////////////////
////////////////counting non-zero pixels P in the mask////////////////
  float *fonesd;
  fonesd=fones->get_data();
  NM=0.0;
  for (i=0;i<nx*ny*nz;i++){
    if (fonesd[i]==1.0) NM=NM+1.0;
    else NM=NM;
  }
  fones->update(); //  doneData->update()
////////////////////////////////////////////////////////////////
///////////////////calculate mean and sigma///////////////////////////////////////
  float *map,sq,sd,mean;
  f2->applyMask(r,type);
  map=f2->get_data();
  float lsum=0.0;
  float sumsq=0.0;
  for (i=0;i<nx*ny*nz;i++){
	  lsum+=map[i];
	  sumsq+=(map[i]*map[i]);
  }
  sq=((NM*sumsq-lsum*lsum)/(NM*NM));
  mean=lsum/NM;
  sd=std::sqrt(sq);
///////////////////////////////////////////////////////
///////////////////normalize S so that simga=1 mean=0/////////////////////////
  float th=0.00001;
  if (sq<0.0)	printf("calcLCF: sd < 0!\n");return NULL;
  if(sq>th)
	  for(i=0;i<nx*ny*nz;i++) map[i]=(map[i]-mean)/sd;
  else
	  for(i=0;i<nx*ny*nz;i++) map[i]=map[i]-mean;
  f2->update();// ming
  f2=f2->copy();//ming?
  f2->applyMask(r,type);
  //ming,f2->toCorner();
  f2->process_inplace("xform.phaseorigin.tocorner");
  //////////////////////////////////////////////////////////////
//////////////////////////////////////////x-correlation between "target map" and "search object"//////////////////////////
  ccf=f1->calc_ccf(f2); //conv(T*S)
  ccf->mult(N); // ming change
/////////////convolution between "target map" and masked region of ones This would give the local average T_M of the target//////////
  conv1=f1->convolute(fones);
  conv1->mult(N);
  conv1->mult((float) 1.0/NM);
////////////////////////////////////////////////////////////////////////////////////////////////
////////////convolution between masked "target map"^2 and masked region of ones///////////////////
  sqf1=f1->copy();//ming
  sqf1->realFilter(18);//x^2  squaring "target map"
  conv2=sqf1->convolute(fones);
  conv2->mult(N);
//////////////////////////////////////////////////////////////////////////////////
///////////////compute variance sigma_mt of masked traget map////////////////
  conv1->realFilter(18);// (T_M)^2
  conv1->mult(float(1.0/(NM*NM))); // (T_M)^2/P^2
  v=conv2->copy();
  v->sub(conv1); //ming
  v->mult((float)1.0/(NM));//then compute std deviation
  v->realFilter(19);//sqrt
//////////////////////////////////////////////////////////////////
  f=ccf->copy();
  f->mult((float)1.0/NM);
  float *lcfd,*vdd;
  lcfd=f->get_data();
  vdd=v->get_data();
  for (i=0;i<nx*ny*nz;i++){
    if (vdd[i]>0.0) lcfd[i]=lcfd[i]/vdd[i];
    else lcfd[i]=lcfd[i];
  }
  f->update();//ming
  v->update();//ming
  lcf=f->copy();
  delete f1;delete f2;delete sqf1;delete ccf;
  delete conv1;delete conv2;delete f;delete fones;delete v;
  return lcf;
}*/


//  The following code looks strange - does anybody know it?  Please let me know, pawel.a.penczek@uth.tmc.edu  04/09/06.
// This is just an implementation of "Roseman's" fast normalized cross-correlation (Ultramicroscopy, 2003). But the contents of this function have changed dramatically since you wrote that comment (d.woolford).
EMData *EMData::calc_flcf(EMData * with)
{
	ENTERFUNC;

	// Ones is a circlular/spherical mask, consisting of 1s.
	EMData* ones = new EMData(with->get_xsize(), with->get_ysize(),with->get_zsize());
	ones->process_inplace("testimage.circlesphere");

	// Get a copy of with, we will eventually resize it
	EMData* with_resized = with->copy();
	// Circular/Spherical mask
	with_resized->mult(*ones);

	// Get with_resized to the right size for the correlation
	Region r((with->get_xsize()-nx)/2, (with->get_ysize()-ny)/2, (with->get_zsize()-nz)/2,nx,ny,nz);
	with_resized->clip_inplace(r);

	// The correlation
	// Get the local sigma image
	EMData* s = calc_fast_sigma_image(ones);
	// The local normalized correlation
	//EMData* corr;
	EMData* corr = calc_ccf(with_resized);
	corr->div(*s);

	delete with_resized;
	delete ones;
	delete s;

	EXITFUNC;
	return corr;
}

EMData *EMData::convolute(EMData * with)
{
	ENTERFUNC;

	EMData *f1 = do_fft();
	if (!f1) {
		LOGERR("FFT returns NULL image");
		throw NullPointerException("FFT returns NULL image");
	}

	f1->ap2ri();

	EMData *cf = 0;
	if (with) {
		cf = with->do_fft();
		if (!cf) {
			LOGERR("FFT returns NULL image");
			throw NullPointerException("FFT returns NULL image");
		}
		cf->ap2ri();
	}
	else {
		cf = f1->copy();
	}

	if (with && !EMUtil::is_same_size(f1, cf)) {
		LOGERR("images not same size");
		throw ImageFormatException("images not same size");
	}

	float *rdata1 = f1->get_data();
	float *rdata2 = cf->get_data();
	size_t cf_size = cf->get_xsize() * cf->get_ysize() * cf->get_zsize();

	float re,im;
	for (size_t i = 0; i < cf_size; i += 2) {
		re = rdata1[i] * rdata2[i] - rdata1[i + 1] * rdata2[i + 1];
		im = rdata1[i + 1] * rdata2[i] + rdata1[i] * rdata2[i + 1];
		rdata2[i]=re;
		rdata2[i+1]=im;
	}

	cf->update();
	EMData *f2 = cf->do_ift();

	if( cf )
	{
		delete cf;
		cf = 0;
	}

	if( f1 )
	{
		delete f1;
		f1=0;
	}

	EXITFUNC;
	return f2;
}


void EMData::common_lines(EMData * image1, EMData * image2,
						  int mode, int steps, bool horizontal)
{
	ENTERFUNC;

	if (!image1 || !image2) {
		throw NullPointerException("NULL image");
	}

	if (mode < 0 || mode > 2) {
		throw OutofRangeException(0, 2, mode, "invalid mode");
	}

	if (!image1->is_complex()) {
		image1 = image1->do_fft();
	}
	if (!image2->is_complex()) {
		image2 = image2->do_fft();
	}

	image1->ap2ri();
	image2->ap2ri();

	if (!EMUtil::is_same_size(image1, image2)) {
		throw ImageFormatException("images not same sizes");
	}

	int image2_nx = image2->get_xsize();
	int image2_ny = image2->get_ysize();

	int rmax = image2_ny / 4 - 1;
	int array_size = steps * rmax * 2;
	float *im1 = new float[array_size];
	float *im2 = new float[array_size];
	for (int i = 0; i < array_size; i++) {
		im1[i] = 0;
		im2[i] = 0;
	}

	set_size(steps * 2, steps * 2, 1);

	float *image1_data = image1->get_data();
	float *image2_data = image2->get_data();

	float da = M_PI / steps;
	float a = -M_PI / 2.0f + da / 2.0f;
	int jmax = 0;

	for (int i = 0; i < steps * 2; i += 2, a += da) {
		float s1 = 0;
		float s2 = 0;
		int i2 = i * rmax;
		int j = 0;

		for (float r = 3.0f; r < rmax - 3.0f; j += 2, r += 1.0f) {
			float x = r * cos(a);
			float y = r * sin(a);

			if (x < 0) {
				x = -x;
				y = -y;
				LOGERR("CCL ERROR %d, %f !\n", i, -x);
			}

			int k = (int) (floor(x) * 2 + floor(y + image2_ny / 2) * image2_nx);
			int l = i2 + j;
			float x2 = x - floor(x);
			float y2 = y - floor(y);

			im1[l] = Util::bilinear_interpolate(image1_data[k],
												image1_data[k + 2],
												image1_data[k + image2_nx],
												image1_data[k + 2 + image2_nx], x2, y2);

			im2[l] = Util::bilinear_interpolate(image2_data[k],
												image2_data[k + 2],
												image2_data[k + image2_nx],
												image2_data[k + 2 + image2_nx], x2, y2);

			k++;

			im1[l + 1] = Util::bilinear_interpolate(image1_data[k],
													image1_data[k + 2],
													image1_data[k + image2_nx],
													image1_data[k + 2 + image2_nx], x2, y2);

			im2[l + 1] = Util::bilinear_interpolate(image2_data[k],
													image2_data[k + 2],
													image2_data[k + image2_nx],
													image2_data[k + 2 + image2_nx], x2, y2);

			s1 += Util::square_sum(im1[l], im1[l + 1]);
			s2 += Util::square_sum(im2[l], im2[l + 1]);
		}

		jmax = j - 1;
		float sqrt_s1 = std::sqrt(s1);
		float sqrt_s2 = std::sqrt(s2);

		int l = 0;
		for (float r = 1; r < rmax; r += 1.0f) {
			int i3 = i2 + l;
			im1[i3] /= sqrt_s1;
			im1[i3 + 1] /= sqrt_s1;
			im2[i3] /= sqrt_s2;
			im2[i3 + 1] /= sqrt_s2;
			l += 2;
		}
	}
	float * data = get_data();

	switch (mode) {
	case 0:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {

				if (m1 == 0 && m2 == 0) {
					for (int i = 0; i < steps; i++) {
						int i2 = i * rmax * 2;
						for (int j = 0; j < steps; j++) {
							int l = i + j * steps * 2;
							int j2 = j * rmax * 2;
							data[l] = 0;
							for (int k = 0; k < jmax; k++) {
								data[l] += im1[i2 + k] * im2[j2 + k];
							}
						}
					}
				}
				else {
					int steps2 = steps * m2 + steps * steps * 2 * m1;

					for (int i = 0; i < steps; i++) {
						int i2 = i * rmax * 2;
						for (int j = 0; j < steps; j++) {
							int j2 = j * rmax * 2;
							int l = i + j * steps * 2 + steps2;
							data[l] = 0;

							for (int k = 0; k < jmax; k += 2) {
								i2 += k;
								j2 += k;
								data[l] += im1[i2] * im2[j2];
								data[l] += -im1[i2 + 1] * im2[j2 + 1];
							}
						}
					}
				}
			}
		}

		break;
	case 1:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {
				int steps2 = steps * m2 + steps * steps * 2 * m1;
				int p1_sign = 1;
				if (m1 != m2) {
					p1_sign = -1;
				}

				for (int i = 0; i < steps; i++) {
					int i2 = i * rmax * 2;

					for (int j = 0; j < steps; j++) {
						int j2 = j * rmax * 2;

						int l = i + j * steps * 2 + steps2;
						data[l] = 0;
						float a = 0;

						for (int k = 0; k < jmax; k += 2) {
							i2 += k;
							j2 += k;

#ifdef	_WIN32
							float a1 = (float) _hypot(im1[i2], im1[i2 + 1]);
#else
							float a1 = (float) hypot(im1[i2], im1[i2 + 1]);
#endif	//_WIN32
							float p1 = atan2(im1[i2 + 1], im1[i2]);
							float p2 = atan2(im2[j2 + 1], im2[j2]);

							data[l] += Util::angle_sub_2pi(p1_sign * p1, p2) * a1;
							a += a1;
						}

						data[l] /= (float)(a * M_PI / 180.0f);
					}
				}
			}
		}

		break;
	case 2:
		for (int m1 = 0; m1 < 2; m1++) {
			for (int m2 = 0; m2 < 2; m2++) {
				int steps2 = steps * m2 + steps * steps * 2 * m1;

				for (int i = 0; i < steps; i++) {
					int i2 = i * rmax * 2;

					for (int j = 0; j < steps; j++) {
						int j2 = j * rmax * 2;
						int l = i + j * steps * 2 + steps2;
						data[l] = 0;

						for (int k = 0; k < jmax; k += 2) {
							i2 += k;
							j2 += k;
#ifdef	_WIN32
							data[l] += (float) (_hypot(im1[i2], im1[i2 + 1]) * _hypot(im2[j2], im2[j2 + 1]));
#else
							data[l] += (float) (hypot(im1[i2], im1[i2 + 1]) * hypot(im2[j2], im2[j2 + 1]));
#endif	//_WIN32
						}
					}
				}
			}
		}

		break;
	default:
		break;
	}

	if (horizontal) {
		float *tmp_array = new float[ny];
		for (int i = 1; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				tmp_array[j] = get_value_at(i, j);
			}
			for (int j = 0; j < ny; j++) {
				set_value_at(i, j, 0, tmp_array[(j + i) % ny]);
			}
		}
		if( tmp_array )
		{
			delete[]tmp_array;
			tmp_array = 0;
		}
	}

	if( im1 )
	{
		delete[]im1;
		im1 = 0;
	}

	if( im2 )
	{
		delete im2;
		im2 = 0;
	}


	image1->update();
	image2->update();
	if( image1 )
	{
		delete image1;
		image1 = 0;
	}
	if( image2 )
	{
		delete image2;
		image2 = 0;
	}
	update();
	EXITFUNC;
}



void EMData::common_lines_real(EMData * image1, EMData * image2,
							   int steps, bool horiz)
{
	ENTERFUNC;

	if (!image1 || !image2) {
		throw NullPointerException("NULL image");
	}

	if (!EMUtil::is_same_size(image1, image2)) {
		throw ImageFormatException("images not same size");
	}

	int steps2 = steps * 2;
	int image_ny = image1->get_ysize();
	EMData *image1_copy = image1->copy();
	EMData *image2_copy = image2->copy();

	float *im1 = new float[steps2 * image_ny];
	float *im2 = new float[steps2 * image_ny];

	EMData *images[] = { image1_copy, image2_copy };
	float *ims[] = { im1, im2 };

	for (int m = 0; m < 2; m++) {
		float *im = ims[m];
		float a = M_PI / steps2;
		Transform t(Dict("type","2d","alpha",-a));
		for (int i = 0; i < steps2; i++) {
			images[i]->transform(t);
			float *data = images[i]->get_data();

			for (int j = 0; j < image_ny; j++) {
				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += data[j * image_ny + k];
				}
				im[i * image_ny + j] = sum;
			}

			float sum1 = 0;
			float sum2 = 0;
			for (int j = 0; j < image_ny; j++) {
				int l = i * image_ny + j;
				sum1 += im[l];
				sum2 += im[l] * im[l];
			}

			float mean = sum1 / image_ny;
			float sigma = std::sqrt(sum2 / image_ny - sum1 * sum1);

			for (int j = 0; j < image_ny; j++) {
				int l = i * image_ny + j;
				im[l] = (im[l] - mean) / sigma;
			}

			images[i]->update();
			a += M_PI / steps;
		}
	}

	set_size(steps2, steps2, 1);
	float *data1 = get_data();

	if (horiz) {
		for (int i = 0; i < steps2; i++) {
			for (int j = 0, l = i; j < steps2; j++, l++) {
				if (l == steps2) {
					l = 0;
				}

				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += im1[i * image_ny + k] * im2[l * image_ny + k];
				}
				data1[i + j * steps2] = sum;
			}
		}
	}
	else {
		for (int i = 0; i < steps2; i++) {
			for (int j = 0; j < steps2; j++) {
				float sum = 0;
				for (int k = 0; k < image_ny; k++) {
					sum += im1[i * image_ny + k] * im2[j * image_ny + k];
				}
				data1[i + j * steps2] = sum;
			}
		}
	}

	update();

	if( image1_copy )
	{
		delete image1_copy;
		image1_copy = 0;
	}

	if( image2_copy )
	{
		delete image2_copy;
		image2_copy = 0;
	}

	if( im1 )
	{
		delete[]im1;
		im1 = 0;
	}

	if( im2 )
	{
		delete[]im2;
		im2 = 0;
	}
	EXITFUNC;
}


void EMData::cut_slice(const EMData *const map, const Transform& transform, bool interpolate)
{
	ENTERFUNC;

	if (!map) throw NullPointerException("NULL image");
	// These restrictions should be ultimately restricted so that all that matters is get_ndim() = (map->get_ndim() -1)
	if ( get_ndim() != 2 ) throw ImageDimensionException("Can not call cut slice on an image that is not 2D");
	if ( map->get_ndim() != 3 ) throw ImageDimensionException("Can not cut slice from an image that is not 3D");
	// Now check for complex images - this is really just being thorough
	if ( is_complex() ) throw ImageFormatException("Can not call cut slice on an image that is complex");
	if ( map->is_complex() ) throw ImageFormatException("Can not cut slice from a complex image");


	float *sdata = map->get_data();
	float *ddata = get_data();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;

	int ymax = ny/2;
	if ( ny % 2 == 1 ) ymax += 1;
	int xmax = nx/2;
	if ( nx % 2 == 1 ) xmax += 1;
	for (int y = -ny/2; y < ymax; y++) {
		for (int x = -nx/2; x < xmax; x++) {
			Vec3f coord(x,y,0);
			Vec3f soln = transform*coord;

// 			float xx = (x+pretrans[0]) * (*ort)[0][0] +  (y+pretrans[1]) * (*ort)[0][1] + pretrans[2] * (*ort)[0][2] + posttrans[0];
// 			float yy = (x+pretrans[0]) * (*ort)[1][0] +  (y+pretrans[1]) * (*ort)[1][1] + pretrans[2] * (*ort)[1][2] + posttrans[1];
// 			float zz = (x+pretrans[0]) * (*ort)[2][0] +  (y+pretrans[1]) * (*ort)[2][1] + pretrans[2] * (*ort)[2][2] + posttrans[2];


// 			xx += map_nx/2;
// 			yy += map_ny/2;
// 			zz += map_nz/2;

			float xx = soln[0]+map_nx/2;
			float yy = soln[1]+map_ny/2;
			float zz = soln[2]+map_nz/2;

			int l = (x+nx/2) + (y+ny/2) * nx;

			float t = xx - floor(xx);
			float u = yy - floor(yy);
			float v = zz - floor(zz);

			if (xx < 0 || yy < 0 || zz < 0 ) {
				ddata[l] = 0;
				continue;
			}
			if (interpolate) {
				if ( xx > map_nx - 1 || yy > map_ny - 1 || zz > map_nz - 1) {
					ddata[l] = 0;
					continue;
				}
				int k = (int) (Util::fast_floor(xx) + Util::fast_floor(yy) * map_nx + Util::fast_floor(zz) * map_nxy);


				if (xx < (map_nx - 1) && yy < (map_ny - 1) && zz < (map_nz - 1)) {
					ddata[l] = Util::trilinear_interpolate(sdata[k],
								sdata[k + 1], sdata[k + map_nx],sdata[k + map_nx + 1],
								sdata[k + map_nxy], sdata[k + map_nxy + 1], sdata[k + map_nx + map_nxy],
								sdata[k + map_nx + map_nxy + 1],t, u, v);
				}
				else if ( xx == (map_nx - 1) && yy == (map_ny - 1) && zz == (map_nz - 1) ) {
					ddata[l] += sdata[k];
				}
				else if ( xx == (map_nx - 1) && yy == (map_ny - 1) ) {
					ddata[l] +=	Util::linear_interpolate(sdata[k], sdata[k + map_nxy],v);
				}
				else if ( xx == (map_nx - 1) && zz == (map_nz - 1) ) {
					ddata[l] += Util::linear_interpolate(sdata[k], sdata[k + map_nx],u);
				}
				else if ( yy == (map_ny - 1) && zz == (map_nz - 1) ) {
					ddata[l] += Util::linear_interpolate(sdata[k], sdata[k + 1],t);
				}
				else if ( xx == (map_nx - 1) ) {
					ddata[l] += Util::bilinear_interpolate(sdata[k], sdata[k + map_nx], sdata[k + map_nxy], sdata[k + map_nxy + map_nx],u,v);
				}
				else if ( yy == (map_ny - 1) ) {
					ddata[l] += Util::bilinear_interpolate(sdata[k], sdata[k + 1], sdata[k + map_nxy], sdata[k + map_nxy + 1],t,v);
				}
				else if ( zz == (map_nz - 1) ) {
					ddata[l] += Util::bilinear_interpolate(sdata[k], sdata[k + 1], sdata[k + map_nx], sdata[k + map_nx + 1],t,u);
				}

//				if (k >= map->get_size()) {
//					cout << xx << " " << yy << " " <<  zz << " " << endl;
//					cout << k << " " << get_size() << endl;
//					cout << get_xsize() << " " << get_ysize() << " " << get_zsize() << endl;
//					throw;
//					}
//
//				ddata[l] = Util::trilinear_interpolate(sdata[k],
//						sdata[k + 1], sdata[k + map_nx],sdata[k + map_nx + 1],
//						sdata[k + map_nxy], sdata[k + map_nxy + 1], sdata[k + map_nx + map_nxy],
//						sdata[k + map_nx + map_nxy + 1],t, u, v);
			}
			else {
				if ( xx > map_nx - 1 || yy > map_ny - 1 || zz > map_nz - 1) {
					ddata[l] = 0;
					continue;
				}
				int k = Util::round(xx) + Util::round(yy) * map_nx + Util::round(zz) * map_nxy;
				ddata[l] = sdata[k];
			}

		}
	}

	update();

	EXITFUNC;
}

#ifdef EMAN2_USING_CUDA
EMData* EMData::cut_slice_cuda(const Transform& transform)
{
	ENTERFUNC;
//
	// These restrictions should be ultimately restricted so that all that matters is get_ndim() = (map->get_ndim() -1)
	if ( get_ndim() != 3 ) throw ImageDimensionException("Can not cut slice from an image that is not 3D");
	// Now check for complex images - this is really just being thorough
	if ( is_complex() ) throw ImageFormatException("Can not call cut slice an image that is complex");
//

	EMData* ret = new EMData();
	ret->set_size_cuda(nx,ny,1);

	float * m = new float[12];
	transform.copy_matrix_into_array(m);

	EMDataForCuda tmp = ret->get_data_struct_for_cuda();
	bind_cuda_texture(); // Binds this image to the global texture
	cut_slice_cuda_(&tmp,m);

	ret->gpu_update();
	delete [] m;

	EXITFUNC;
	return ret;
}

#endif // EMAN2_USING_CUDA


void EMData::uncut_slice(EMData * const map, const Transform& transform) const
{
	ENTERFUNC;

	if (!map) throw NullPointerException("NULL image");
	// These restrictions should be ultimately restricted so that all that matters is get_ndim() = (map->get_ndim() -1)
	if ( get_ndim() != 2 ) throw ImageDimensionException("Can not call cut slice on an image that is not 2D");
	if ( map->get_ndim() != 3 ) throw ImageDimensionException("Can not cut slice from an image that is not 3D");
	// Now check for complex images - this is really just being thorough
	if ( is_complex() ) throw ImageFormatException("Can not call cut slice on an image that is complex");
	if ( map->is_complex() ) throw ImageFormatException("Can not cut slice from a complex image");

// 	Transform3D r( 0, 0, 0); // EMAN by default
// 	if (!ort) {
// 		ort = &r;
// 	}

	float *ddata = map->get_data();
	float *sdata = get_data();

	int map_nx = map->get_xsize();
	int map_ny = map->get_ysize();
	int map_nz = map->get_zsize();
	int map_nxy = map_nx * map_ny;
	float map_nz_round_limit = (float) map_nz-0.5f;
	float map_ny_round_limit = (float) map_ny-0.5f;
	float map_nx_round_limit = (float) map_nx-0.5f;
/*
	Vec3f posttrans = ort->get_posttrans();
	Vec3f pretrans = ort->get_pretrans();*/

	int ymax = ny/2;
	if ( ny % 2 == 1 ) ymax += 1;
	int xmax = nx/2;
	if ( nx % 2 == 1 ) xmax += 1;
	for (int y = -ny/2; y < ymax; y++) {
		for (int x = -nx/2; x < xmax; x++) {
			Vec3f coord(x,y,0);
			Vec3f soln = transform*coord;
// 			float xx = (x+pretrans[0]) * (*ort)[0][0] +  (y+pretrans[1]) * (*ort)[0][1] + pretrans[2] * (*ort)[0][2] + posttrans[0];
// 			float yy = (x+pretrans[0]) * (*ort)[1][0] +  (y+pretrans[1]) * (*ort)[1][1] + pretrans[2] * (*ort)[1][2] + posttrans[1];
// 			float zz = (x+pretrans[0]) * (*ort)[2][0] +  (y+pretrans[1]) * (*ort)[2][1] + pretrans[2] * (*ort)[2][2] + posttrans[2];
//
// 			xx += map_nx/2;
// 			yy += map_ny/2;
// 			zz += map_nz/2;
//
			float xx = soln[0]+map_nx/2;
			float yy = soln[1]+map_ny/2;
			float zz = soln[2]+map_nz/2;

			// These 0.5 offsets are here because the round function rounds to the nearest whole number.
			if (xx < -0.5 || yy < -0.5 || zz < -0.5 || xx >= map_nx_round_limit || yy >= map_ny_round_limit || zz >= map_nz_round_limit) continue;

			int k = Util::round(xx) + Util::round(yy) * map_nx + Util::round(zz) * map_nxy;
			int l = (x+nx/2) + (y+ny/2) * nx;
			ddata[k] = sdata[l];
		}
	}

	map->update();
	EXITFUNC;
}


void EMData::save_byteorder_to_dict(ImageIO * imageio)
{
	string image_endian = "ImageEndian";
	string host_endian = "HostEndian";

	if (imageio->is_image_big_endian()) {
		attr_dict[image_endian] = "big";
	}
	else {
		attr_dict[image_endian] = "little";
	}

	if (ByteOrder::is_host_big_endian()) {
		attr_dict[host_endian] = "big";
	}
	else {
		attr_dict[host_endian] = "little";
	}
}
// Ming, the following lines added for new aligner "FRM2D"
// Determine center of mass (more or less), move object to the COM, RETURN COMs, YCong

// return a random number between lo and hi
float frand(float lo, float hi)
{
static int f=1;

if (f) { srandom(time(0)); f=0; }

return ((float)random()/2147483647.0*(hi-lo)+lo);
}

float grand(float mean,float sigma)
{
float x,y,r,f;

do {
	x=frand(-1.0,1.0);
	y=frand(-1.0,1.0);
	r=x*x+y*y;
} while (r>1.0||r==0);

f=sqrt(-2.0*log(r)/r);
return x*f*sigma+mean;
}


void *List::objectAt(int nn)
{
if (this==NULL || nn<0 || nn>=n) return NULL;
return items[nn];
}

void Euler::init() {
strcpy(sym,"i");
nsym=1;
csym=0;
}
void Euler::rectify(){
#ifdef DEBUG0
	printf("rectify %f,%f,%f ->",alt_,az_,phi_);
#endif
	if (!Util::goodf(&phi_)) phi_=0;//ming change
	if (!Util::goodf(&alt_)) alt_=0; // ming change
	if (!Util::goodf(&az_)) az_=0; // ming change
	if (sqrt(alt_*alt_) <= .0001) { az_ += phi_; phi_=0;}

	phi_=fmod(phi_,(float)(2.0*PI));
	if (phi_<0) phi_+=2.0*PI;

	az_=fmod(az_,(float)(2.0*PI));
	if (az_<0) az_+=2.0*PI;
#ifdef DEBUG0
	printf(" %f,%f,%f\n",alt_,az_,phi_);
#endif
}
Euler::Euler() {
alt_=0.; az_=0.; phi_=0.;
init();
} // already rectified

Euler::~Euler() {
}

void Euler::setAngle(float a1, float a2, float a3, int typ){
#ifdef DEBUG0
	printf("setAngle %f %f %f %d\n",a1,a2,a3,typ);
#endif
	alt_ =a1 ; az_ = a2 ; phi_ = a3;
	if (typ == IMAGIC) { /*az_ +=PI; phi_ +=PI;*/} // then a1= beta, a2=gamma, a3=alpha
	if (typ == SPIDER) { az_ = a2-PI/2.-PI; phi_ = a3+PI/2.-PI;}    // then a1=theta, a2=phi, a3=gamma
    	if (typ !=EMAN && typ !=IMAGIC && typ != SPIDER) {
			printf("unknown declaration of type for Euler constructor /n");}
			rectify();
}
float  Euler::az() const {return az_  ;}

float  Euler::alt() const { return alt_ ;}
float  Euler::phi() const { return phi_ ;}


void EMData::zero(float sig) {
	get_data(); // ming change getData() to get_data()
	if (sig==0) memset(rdata,0,nx*ny*nz*sizeof(float));
	else {
		int i;
		for (i=0; i<nx*ny*nz; i++) rdata[i]=grand(0,sig);
	}
	//ming comment doneData();
	update();
}

void EMData::rotateAndTranslate(float scale,float dxc,float dyc,float dzc,int r) {
	float *ddata,*sdata;
	float yy,x2,y2,z2,x2c,y2c;
	int x, y, z;
	float t,u,v,p0,p1,p2,p3,em,fl;
	//static int *limit = (int *)malloc(4),nsz,nrad;
	int k0,k1,k2,k3;
	int i,j,k,ii,jj;
	int nx2=nx, ny2=ny;	// deal with the special case of rotating different size parent to this. only 2D is supported now

	// This delays in case of a radius mismatch in parallel processing
	/*while (r>0 && nsz) usleep(1000);

	if (r && r!=abs(nlimit)) {
		nlimit=r;
		limit=(int *)realloc(limit,r*2sizeof(int));
	}*/

	r=r*r;

	if (scale!=0) scale=1.0/scale; else scale=1.0;

	if (parent) {
		sdata=parent->get_data(); // ming change getDataRO() to get_data()
		ddata=get_data(); // ming change getData() to get_data()
		nx2=parent->get_xsize(); //ming
		ny2=parent->get_ysize(); // ming
	}
	else {
		sdata=get_data(); // ming change getData() to get_data()
		//ddata=(float *)Smalloc(nx*ny*nz*sizeof(float),1);
		ddata=(float *)malloc(nx*ny*nz*sizeof(float));
	}

	if (nz==1) {
		//em=edgeMean();
		em=0;
		float mx0=scale*cos(euler.alt());
		float mx1=scale*sin(euler.alt());
		x2c=nx/2.0-dxc-dx;
		y2c=ny/2.0-dyc-dy;
		for (j=0,y=-ny/2.0+dxc; j<ny; j++,y+=1.0) {
			for (i=0,x=-nx/2.0+dyc; i<nx; i++,x+=1.0) {

				if (r && (i-nx/2)*(j-nx/2)>r) {
					ddata[i+j*nx]=0;
					continue;
				}

		// rotate coordinates
				x2=(mx0*x+mx1*y)+x2c;
				y2=(-mx1*x+mx0*y)+y2c;

				if (x2<0||x2>nx2-1.0||y2<0||y2>ny2-1.0) {
//			ddata[i+j*nx]=em;
//			ddata[i+j*nx]=grand(Mean(),Sigma());
//			ddata[i+j*nx]=Mean();
					ddata[i+j*nx]=0;
					continue;
				}

				ii=Util::fast_floor(x2);
				jj=Util::fast_floor(y2);
				k0=ii+jj*nx;
				k1=k0+1;
				k2=k0+nx+1;
				k3=k0+nx;

				if (ii==nx2-1) { k1--; k2--; }
				if (jj==ny2-1) { k2-=nx2; k3-=nx2; }

		// no interpolation
		//ddata[i+j*nx]=sdata[ii+jj*nx];

		// bilinear interpolation of raw data
				t=(x2-(float)ii);
				u=(y2-(float)jj);
				float tt=1.0-t;
				float uu=1.0-u;

				p0=sdata[k0]*tt*uu;
				p1=sdata[k1]*t*uu;
				p3=sdata[k3]*tt*u;
				p2=sdata[k2]*t*u;
				ddata[i+j*nx]=p0+p1+p2+p3;
			}
		}
	}
    else if(get_map_type() == ICOS5F_HALF) { // to deal with the special icosahedral 5fold half map which has 2fold axis on y
	//else if(1==1) { // ming add
		float mx[9];
		float x3,y3,z3, x4,y4,z4;
		int l=0,xy;

		mx[0]=scale*(cos(euler.phi())*cos(euler.az())-cos(euler.alt())*sin(euler.az())*sin(euler.phi()));
		mx[1]=-scale*(sin(euler.phi())*cos(euler.az())+cos(euler.alt())*sin(euler.az())*cos(euler.phi()));
		mx[2]=scale*sin(euler.alt())*sin(euler.az());
		mx[3]=scale*(cos(euler.phi())*sin(euler.az())+cos(euler.alt())*cos(euler.az())*sin(euler.phi()));
		mx[4]=scale*(-sin(euler.phi())*sin(euler.az())+cos(euler.alt())*cos(euler.az())*cos(euler.phi()));
		mx[5]=-scale*sin(euler.alt())*cos(euler.az());
		mx[6]=scale*sin(euler.alt())*sin(euler.phi());
		mx[7]=scale*sin(euler.alt())*cos(euler.phi());
		mx[8]=scale*cos(euler.alt());

		for (int p = 0; p < 9; p++) {
			printf("%4.2f ", mx[p]);
		}
		printf("\n");


		xy=nx*ny;
		printf("in rotation\n");
		for (int k=0; k<nz; k++) {
			for (int j=-ny/2; j<ny-ny/2; j++) {
				for (int i=-nx/2; i<nx-nx/2; i++, l++) {
					x2=mx[0]*i+mx[1]*j+mx[2]*k+nx/2;
					y2=mx[3]*i+mx[4]*j+mx[5]*k+ny/2;
					z2=mx[6]*i+mx[7]*j+mx[8]*k+0/2;

					if (x2<0||y2<0||z2<=-(nz-1)|| x2>=nx||y2>=ny||z2>=nz-1) {
						continue;
					}

					if(z2<0){
						x2=nx/2*2-x2;
						z2=-z2;
					}

					x=Util::fast_floor(x2);
					y=Util::fast_floor(y2);
					z=Util::fast_floor(z2);

					t=x2-x;
					u=y2-y;
					v=z2-z;

				// real part
					ii=(int)(x+y*nx+z*xy);

			 	ddata[l]+=	Util::trilinear_interpolate(sdata[ii],sdata[ii+1],sdata[ii+nx],sdata[ii+nx+1],
							sdata[ii+nx*ny],sdata[ii+xy+1],sdata[ii+xy+nx],sdata[ii+xy+nx+1],
							t,u,v);
				}
			}
		}
	}
	else {
		float mx[9];
		float x3,y3,z3, x4,y4,z4;
		int l,xy,mr;

		mx[0]=scale*(cos(euler.phi())*cos(euler.az())-cos(euler.alt())*sin(euler.az())*sin(euler.phi()));
		mx[1]=-scale*(sin(euler.phi())*cos(euler.az())+cos(euler.alt())*sin(euler.az())*cos(euler.phi()));
		mx[2]=scale*sin(euler.alt())*sin(euler.az());
		mx[3]=scale*(cos(euler.phi())*sin(euler.az())+cos(euler.alt())*cos(euler.az())*sin(euler.phi()));
		mx[4]=scale*(-sin(euler.phi())*sin(euler.az())+cos(euler.alt())*cos(euler.az())*cos(euler.phi()));
		mx[5]=-scale*sin(euler.alt())*cos(euler.az());
		mx[6]=scale*sin(euler.alt())*sin(euler.phi());
		mx[7]=scale*sin(euler.alt())*cos(euler.phi());
		mx[8]=scale*cos(euler.alt());

		x4=(mx[0]*(-nx/2.0+dxc)+mx[1]*(-ny/2.0+dyc)+mx[2]*(-nz/2.0+dzc))+nx/2.0-dxc-dx;
		y4=(mx[3]*(-nx/2.0+dxc)+mx[4]*(-ny/2.0+dyc)+mx[5]*(-nz/2.0+dzc))+ny/2.0-dyc-dy;
		z4=(mx[6]*(-nx/2.0+dxc)+mx[7]*(-ny/2.0+dyc)+mx[8]*(-nz/2.0+dzc))+nz/2.0-dzc-dz;

		xy=nx*ny;
		if (nx>ny && nx>nz) mr=(nx-2)*(nx-2)/4;
		else if (ny>nz) mr=(ny-2)*(ny-2)/4;
		else mr=(nz-2)*(nz-2)/4;

		for (k=-nz/2,l=0; k<nz/2; k++,x4+=mx[2],y4+=mx[5],z4+=mx[8]) {
			x3=x4;
			y3=y4;
			z3=z4;
			for (j=-ny/2; j<ny/2; j++,x3+=mx[1],y3+=mx[4],z3+=mx[7]) {
				x2=x3;
				y2=y3;
				z2=z3;
				for (i=-nx/2; i<nx/2; i++,l++,x2+=mx[0],y2+=mx[3],z2+=mx[6]) {

					if (x2<0||y2<0||z2<0|| x2>=nx-1||y2>=ny-1||z2>=nz-1) {
//				if (i*i+j*j+k*k>=mr) {
						ddata[l]=0;
						continue;
					}

					x=Util::fast_floor(x2);
					y=Util::fast_floor(y2);
					z=Util::fast_floor(z2);

					t=x2-x;
					u=y2-y;
					v=z2-z;

				// real part
					ii=(int)(x+y*nx+z*xy);

				ddata[l]=Util::trilinear_interpolate(sdata[ii],sdata[ii+1],sdata[ii+nx],sdata[ii+nx+1],
					sdata[ii+nx*ny],sdata[ii+xy+1],sdata[ii+xy+nx],sdata[ii+xy+nx+1],
					t,u,v);
				}
			}
		}

/*	for (k=0,z=-nz/2.0+dzc; k<nz; k++,z+=1.0) {
		for (j=0,y=-ny/2.0+dyc; j<ny; j++,y+=1.0) {
			for (i=0,x=-nx/2.0+dxc; i<nx; i++,x+=1.0) {
				x2=(mx[0]*x+mx[1]*y+mx[2]*z)+nx/2.0-dxc-dx;
				y2=(mx[3]*x+mx[4]*y+mx[5]*z)+ny/2.0-dyc-dy;
				z2=(mx[6]*x+mx[7]*y+mx[8]*z)+nz/2.0-dzc-dz;

				if (x2<0||y2<0||z2<0|| x2>=nx-1||y2>=ny-1||z2>=nz-1) {
					ddata[i+j*nx+k*nx*ny]=0;
					continue;
				}

				t=x2-floor(x2);
				u=y2-floor(y2);
				v=z2-floor(z2);

				// real part
				ii=(int)(floor(x2)+floor(y2)*nx+floor(z2)*nx*ny);

				ddata[i+j*nx+k*nx*ny]=trilin(sdata[ii],sdata[ii+1],sdata[ii+nx],sdata[ii+nx+1],
					sdata[ii+nx*ny],sdata[ii+nx*ny+1],sdata[ii+nx*ny+nx],sdata[ii+nx*ny+nx+1],
					t,u,v);
			}
		}
	}*/

	}
	float pixel=1.0f;// ming add
   setXYZOrigin(getXorigin()+nx/2*pixel*(1-scale)-dx*scale*pixel,
				getYorigin()+ny/2*pixel*(1-scale)-dy*scale*pixel,
				getZorigin()+nz/2*pixel*(1-scale)-dz*scale*pixel
			);

   setPixel(eman1_Pixel()*scale);
   float ctf[11];// ming add
 if (hasCTF()) ctf[10]*=scale;
	if (parent) parent->update(); // ming change doneData() to update()
	else {
		free(rdata);
		rdata=ddata;
		dx=dy=dz=0;
		setRAlign(0,0,0);
//	flags&=~EMDATA_FLIP;
	}
	//ming comment, doneData();

	update();
	return;
}




// This does fast integer-only translation
void EMData::fastTranslate(int inplace) {
	float *d;
	float *dd;
	int x0,x1,x2,y0,y1,y2;

	if (parent&&(!inplace)) {
		//ming comment, d=parent->getDataRO();
		d=parent->get_data(); // ming change getData
		//ming comment, dd=getData();
		dd=get_data(); //ming change
	}
	else {
		d=dd=get_data(); // ming change
	}


	if (nz==1) {
		if (dx<0) { x0=0; x1=get_xsize(); x2=1; } //ming
		else { x0=get_xsize()-1; x1= -1; x2=-1; } // ming

		if (dy<0) { y0=0; y1=get_ysize(); y2=1; } // ming
		else { y0=get_ysize()-1; y1=-1; y2=-1; }  //ming

		for (int x=x0; x!=x1; x+=x2) {
			for (int y=y0; y!=y1; y+=y2) {
				int xp=x-dx;
				int yp=y-dy;
				if (xp<0||yp<0||xp>=get_xsize()||yp>=get_ysize()) { dd[x+y*get_xsize()]=0; continue; } // ming
				dd[x+y*get_xsize()]=d[xp+yp*get_xsize()]; // ming
			}
		}

		// ming comment, doneData();
		if (d==dd) dx=dy=0;
		else parent->update(); // ming change doneData() to update()
	}
	else {
	   Util::error(ERR_ERROR,"Sorry, 3d translateFast not supported yet");
	}

}

/******ming
EMData *EMData::doFFT() {	// obvious, note that result is R,I not A,P
	EMData *dat;
	float *d,re,im,norm;
	int i,j,k,l,ii,jj,dim;
	char s[80],*t;
	mrcH *newh;
	rfftwnd_plan cplan;

	if (flags&EMDATA_COMPLEX) return this;	// already complex, just return

	if (nx!=Util::npfa(nx) || (ny!=Util::npfa(ny)&&ny!=1) || (nz!=Util::npfa(nz)&&nz!=1)) { Util::error("Invalid FFT size!",ERR_ERROR); return NULL; }

	getData();	// we open for writing so nobody else does a fft at the same time

	//flags|=EMDATA_NEWFFT;		// disable caching

	if (fft) {
		if (!(flags&EMDATA_NEWFFT)) {
			i=flags;
			doneData();
			flags=i&~EMDATA_BUSY;
			//		printf("old ok\n");
			return fft;
		}
		fft->setSize(nx+2,ny,nz);
		dat=fft;
		//	printf("old fft\n");
	}
	else {
		//	printf("new fft\n");
		dat=copyHead();
		//	sprintf(s,"FFT %s",name);
		//	dat->setName(s);
		fft=dat;
		fft->setSize(nx+2,ny,nz);
	}


	if (ny==1) dim=1;
	else if (nz==1) dim=2;
	else dim=3;

	cplan=makeplan(dim,nx,ny,nz,1);
	d=dat->getData();

	// 1d
	if (ny==1) {
		j=0;
		k=1;
		rfftwnd_one_real_to_complex(cplan,(fftw_real *)rdata,(fftw_complex *)d);
	}
	// 2d
	else if (nz==1) {
		rfftwnd_one_real_to_complex(cplan,(fftw_real *)rdata,(fftw_complex *)d);

		// This 'rotates' the data vertically by ny/2
		l=ny/2*(nx+2);
		for (i=0; i<ny/2; i++) {
			for (j=0; j<nx+2; j++) {
				k=j+i*(nx+2);
				re=d[k];
				d[k]=d[k+l];
				d[k+l]=re;
			}
		}
	}
	// 3d
	else {
		rfftwnd_one_real_to_complex(cplan,(fftw_real *)rdata,(fftw_complex *)d);

		// in the y,z plane, swaps quadrants I,III and II,IV
		// this is the 'rotation' in y and z
		t=(char *)malloc(sizeof(float)*(nx+2));

		k=(nx+2)*ny*(nz+1)/2;
		l=(nx+2)*ny*(nz-1)/2;
		jj=(nx+2)*sizeof(float);
		for (j=ii=0; j<nz/2; j++) {
			for (i=0; i<ny; i++,ii+=nx+2) {
				memcpy(t,d+ii,jj);
				if (i<ny/2) {
					memcpy(d+ii,d+ii+k,jj);
					memcpy(d+ii+k,t,jj);
				}
				else {
					memcpy(d+ii,d+ii+l,jj);
					memcpy(d+ii+l,t,jj);
				}
			}
		}
		free(t);
	}

	fft=dat;
//	dat->multConst(1.0/(nx*ny*nz));
	for (i=0; i<dat->xSize()*dat->ySize()*dat->zSize(); i++) d[i]*=(1.0/(nx*ny*nz));
	dat->doneData();
//	dat->update();
	dat->flags|=EMDATA_COMPLEX|EMDATA_RI;
	i=flags;	// we didn't acutally change 'this', so we don't need an update
	doneData();
	flags=i&~(EMDATA_NEWFFT|EMDATA_BUSY);
	return dat;
}

float *EMData::fsc(EMData *with) {
if (!with || !this) { Util::error(ERR_ERROR,"Null FSC"); return NULL; }

 if (with->get_xsize() != get_xsize() || with->ySize() != ySize() || with->zSize()!=zSize()) {
	error(ERR_ERROR,"Size mismatch in FSC");
	return NULL;
}

float *ret = (float *)calloc(1,get_ysize()/2*sizeof(float));
float *n1 = (float *)calloc(1,get_ysize()/2*sizeof(float));
float *n2 = (float *)calloc(1,get_ysize()/2*sizeof(float));

EMData *f1;
if (isComplex()) f1=this; else f1=doFFT();
f1->ap2ri();

EMData *f2;
if (with->isComplex()) f2=with; else f2=with->doFFT();
f2->ap2ri();

float *d1 = f1->getDataRO();
float *d2 = f2->getDataRO();

if (nz==1) {
	int i,j;

	for (j=0; j<ySize(); j++) {
		for (i=0; i<f1->xSize(); i+=2) {
			float r=hypot((float)(i/2),(float)j-ySize()/2);
			int ri=round(r);
			if (ri<1||ri>=ySize()/2) continue;
			int ii=i+j*f1->xSize();
			ret[ri]+=d1[ii]*d2[ii]+d1[ii+1]*d2[ii+1];
			n1[ri]+=d1[ii]*d1[ii]+d1[ii+1]*d1[ii+1];
			n2[ri]+=d2[ii]*d2[ii]+d2[ii+1]*d2[ii+1];
		}
	}
}
else {
	int i,j,k;

	for (k=0; k<zSize(); k++) {
		for (j=0; j<ySize(); j++) {
			for (i=0; i<f1->xSize(); i+=2) {
				float r=sqrt(hypot3((float)(i/2),(float)j-ySize()/2,(float)k-zSize()/2));
				int ri=round(r);
				if (ri<1||ri>=ySize()/2) continue;
				int ii=i+j*f1->xSize()+k*f1->xSize()*ySize();
				ret[ri]+=d1[ii]*d2[ii]+d1[ii+1]*d2[ii+1];
				n1[ri]+=d1[ii]*d1[ii]+d1[ii+1]*d1[ii+1];
				n2[ri]+=d2[ii]*d2[ii]+d2[ii+1]*d2[ii+1];
			}
		}
	}
}

int i;
for (i=0; i<ySize()/2; i++) ret[i]=ret[i]/(sqrt(n1[i]*n2[i]));

f1->doneData();
f2->doneData();
ret[0]=1.0;
free(n1);
free(n2);
return ret;
}

// similarity between 2 images using FSC. Larger is better.
float EMData::fscmp(EMData *with,float *snr) {
static float *dfsnr = NULL; //(float *)malloc(200*sizeof(float));
static int nsnr=0;

float *FSC = fsc(with);
int i;
int CTFOS=5; // ming add

if (!snr) {
	snr=dfsnr;

	int NP=ceil(CTFOS*std::sqrt(2.0)*ny/2)+2;//??

	if (nsnr!=NP) {
		nsnr=NP;
		dfsnr=(float *)realloc(dfsnr,NP*sizeof(float));
		snr=dfsnr;

//		float w=SQR(xSize()/4.0);
		float w=SQR(get_xsize()/8.0); // ming
		for (int i=0; i<NP; i++) {
			float x2=SQR(i/(float)CTFOS);
			snr[i]=(1.0-exp(-x2/4.0))*exp(-x2/w);
		}
//		save_data(0,1.0/CTFOS, snr, NP, "filt.txt");
	}
}
//save_data(0,1.0, FSC, ny/2, "fsc.txt");

float sum=0,norm=0;
for (i=0; i<ny/2; i++) {

	sum+=FSC[i]*i*snr[i*CTFOS+CTFOS/2];
	norm+=i*snr[i*CTFOS+CTFOS/2];
//	sum+=FSC[i]*snr[i*CTFOS+CTFOS/2];
//	norm+=snr[i*CTFOS+CTFOS/2];
}

//for (i=0; i<ny/2; i++) FSC[i]=FSC[i]*i*snr[i*CTFOS+CTFOS/2];
//save_data(0,1.0, FSC, ny/2, "fsc2.txt");

free(FSC);
return sum/norm;
}


void EMData::realupdate() {
	int i,j,k,l,di,n=0;
	double u,v,w;
	double m,s,t,m2,s2;

	if (nx==0 || ny==0 || nz==0) { set_size(2,2,2); zero(); } // ming chang all setSize to set_size

	flags&=~EMDATA_NEEDUPD;

	max=-1.0e20; // meaning?
	min=-max;  // meaning?
	m=m2=0;
	s=s2=0;

	// getDataRO();		//dangerous not to, but causes problems

	//for (i=0; i<200; i++) hist[i]=0;
	if ((flags&EMDATA_COMPLEX)&& !(flags&EMDATA_RI)) di=2; else di=1;

	// pass 1, statistics
	for (j=0,i=0; j<nz; j++) {
		for (k=0; k<ny; k++) {
			for (l=0; l<nx; l+=di,i+=di) {
				v=rdata[i];
				//			if (!(finite(v))) printf("wrong number at x,y,z=%d,%d,%d = %f\n",l,k,j,v);
				//	if (v>1e5) { printf(" %d,%d=%f\n",i%nx,i/nx,v); continue; }
				if (v>max) { max=v; maxloc=i; }
				if (v<min) { min=v; minloc=i; }
				m+=v;
				s+=v*v;
				if (v!=0) n++;

				//			if (k<4 || l<4 || k>ny-5 || l>nx-5) {
				//m2+=v;
				//s2+=v*v;
			//}
			}
		}
	}

	sumsq=s*di;
	mean=m*di/(double)(nx*ny*nz);
	sigma=s*di/(double)(nx*ny*nz)-mean*mean;
	if (sigma<0) sigma=0; // in case of roundoff errors for v small sigma
	//if (fabs(sigma)<1.0e-10) sigma=0;		// could be that a 'zero' is a bit off due to rounding
	sigma=std::sqrt(sigma);

	if (n==0) n=1;
	mean2=m*di/(double)(n);
	sigma2=s*di/(double)(n)-mean2*mean2;
	if (sigma2<0) sigma2=0;	// in case of roundoff errors for v small sigma
	sigma2=std::sqrt(sigma2);

	//histmax=max;
	//histmin=min;

	//	printf("!!! %d %g %g %g %g\n",(nx*ny*nz),min,max,(float)m,(float)s);
	//	printf("!!! %d %g %g %g %g\n",(nx*ny*nz),min,max,mean,sigma);

	if (!finite(mean) || !finite(sigma)) {
		char st[120];
		sprintf(st,"ERROR - Bad statistics on this data set! %d %g %g",(nx*ny*nz),min,max);
		//	for (i=0; i<110; i++) printf("%5.3g\t",rdata[i+6050]);
		//	printf("\n");

		Util::error(st,ERR_ERROR);

		zero();
		//	kill(0,16);
		return;
	}

	//printf("ru min %f max %f mean %f sigma %f\n",min,max,mean,sigma);

	//if (histmax-mean>sigma*3.0) histmax=mean+sigma*3.0;
	//if (mean-histmin>sigma*3.0) histmin=mean-sigma*3.0;
	//if (histmax==histmin) histmax+=.0001;
	//if (max==min) max+=.0001;

	// pass 2, histogram
	//calcHist(histmin ,histmax);

	return;
}
***********ming*****/

/*************** ming begin comment
float *EMData::getDataRO()	// This returns a pointer to the raw float data
{
	//if (flags&EMDATA_BUSY || !rdata) return NULL;

	waitReady(1);
	if (flags&EMDATA_SWAPPED) swapin();
	rocount++;
	return rdata;
}
******** ming end******************/

/********************
float EMData::frm2dAlign(EMData *with, float *frm2dhhat, EMData* selfpcsfft, float rmaxH, float MAXR_f, float tstep, float angstep, float r_size_flo, float &com_self_x, float &com_self_y, float &com_with_x, float &com_with_y, int neg_rt)
{ // this is the core part of FRM2D aligner
	float nang=angstep;
	int size=(int)(nang+0.5f);
	// float nang = 360.0/angstep;
	// size = (int)(nang+0.5f);
	int bw=size/2;
	unsigned long tsize2=(size+2)*size;
	unsigned long tsize=2*size;
	unsigned long size2=size*size;
	unsigned long ind1=0, ind2=0, ind3=0, ind4=0, ind41=0;
	unsigned long index0=0, indexn=0, indexn2=0;
	int i=0, j=0, k=0, m=0, n=0, r=0;
	int tm=0, tm1=0, ipm=0, tipm=0, tipm1=0, loop_rho=0, rho_best=0;
	float rhomax=0.0f, rhomin = 0.0f;  // ming, rhomax is maximum distance difference between max distance of  image to center

	int MAXR = (int)(MAXR_f+0.5f);
	//float p=0.0;
	FILE *outs;


	int r_size=(int)(r_size_flo);
 	int p_max=r_size-(int)(rmaxH);
	//if(neg_rt==1) p_max=4;
	if(neg_rt==2) p_max=5; // ming remove the comment

	float* gnr2   = new float[size*2];
 	float* maxcor = new float[MAXR+1];                  // MAXR need change
	float* result = new float[5*(p_max+1)];
	if( gnr2 == NULL || maxcor == NULL || result == NULL){
		cout <<"can't allocate more memory, terminating. \n";
		exit(1);
	}



 	// The following 2 blocks are for inverse 2D FFT of Tri_2D (2D correlation value)
	//fftw_real c[size][size];
	//fftw_complex C[size][bw+1];
	  fftw_real * c = new fftw_real[size * size];   // ming, correlation of real space
	  fftw_complex * C = new fftw_complex[size * (bw+1)]; // ming, correlation in fourier space which is complex number.

	 rfftwnd_plan pinv;                           // Inverse 2D FFW plan; see FFTW menu p.7-11
	 pinv = rfftw2d_create_plan(size, size, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);  //out-of-place translation

	float *self_sampl_fft = selfpcsfft->get_data(); // ming change getDataRO() to get_data()
	selfpcsfft->update();   // ming change doneData() to update()                            //~~~ delete self_sampl_fft somewhere

	float maxcor_sofar=0.0f;
	int p=0;
	//printf("pi_max=%f, rhomax=%f\n",(rhomax-rhomin)/tstep, rhomax); // ming
	//for(pi=0; pi<=int((rhomax-rhomin)/tstep+0.5f); pi++){     // loop for rho->p

	for(p=0; p<=p_max; p++){
		//float p = rhomin + pi * tstep;
		ind1 = p*r_size*bw;
		int r_va1 = p-rmaxH;
		if(r_va1 < 1) r_va1 = 1;
		int r_va2 = p+rmaxH;

		for (i=0;i<size;i++)
			for (j=0;j<bw+1;j++){
				//C[i][j].re = 0.0f;
				//C[i*(bw+1) + j].re = 0.0f;
				//C[i][j].im = 0.0f;
			C[i*(bw+1) + j].im = 0.0f;
			}

    	for(n=0;n<bw;n++){                                // loop for n
      		ind2 = (ind1+n);
      		index0=n*(bw+1);

      			//for(r=1;r<=(r_size-1);r++) {          // add terms corresponding to r to each integral
      		for(r=r_va1;r<=r_va2;r++) {
      			ind3 = (ind2 + r*bw)*size;

      			for(m=0;m<size;m++){              // take back hat{h(n,r,p)}(m)
      				ind4 = (ind3+m)*2;
      				ind41= ind4+1;
      				tm=m*2; tm1=tm+1;
					//printf("p, n, r,m =%d, %d, %d, %d\n", p,n,r,m);
      				gnr2[tm]  = frm2dhhat[ind4];
      				gnr2[tm1] = frm2dhhat[ind41];
      			}

				/****** integral( h(+-n,r,p)(m)*conjugate(fm(r)) ) ******/

				/* The m is positive because only positive m is need to compute Tri_2D,  *
				*  since its inverse FT is real.  (see the FFTW manu., for complex to    *
				*  real (inverse FFT), the complex array containes only half-size of the *
				*  last demesion, like X(n1*n2*...*(nd/2+1)), since                      *
				*  X[j1,j2,...,jd] = X[n1-j1,n2-j2,...,nd-jd]*, with hermitian symmetry).*
				*  Only 1 index can be restricted to positive value, so n has both       *
				*  + and - value                                                         *

      			for(m=0;m<bw;m++){
      				tm=m*2; tm1=tm+1;
      				ipm=index0+m; tipm=ipm*2; tipm1=tipm+1;               // index0=n*(bw+1)

      				float tempr = self_sampl_fft[tm + r*(size+2)] * r;
      				float tempi = self_sampl_fft[tm1+ r*(size+2)] * r;
      				float gnr2_r = gnr2[tm];
      				float gnr2_i = gnr2[tm1];

					//C[n][m].re += gnr2_r * tempr + gnr2_i * tempi;     //  m,+n
					//C[n][m].im += gnr2_i * tempr - gnr2_r * tempi;     // (a+ib)(c-id)=(ac+bd)+i(bc-ad)
					C[n*(bw+1) + m].re += gnr2_r * tempr + gnr2_i * tempi;     //  m,+n
					C[n*(bw+1) + m].im += gnr2_i * tempr - gnr2_r * tempi;     // (a+ib)(c-id)=(ac+bd)+i(bc-ad)

      				if(n != 0){					// m,-n
      					int imm=index0-m; int timm=imm*2;	// timm = n*(size+2)+2m
      					int veri =tsize2-timm;			// veri = (size-n)(size+2)+2m
      					int veri1=veri+1;		// Tri_2D size:  m: 0->(size+2), n: 0->(size-1)
      					if(m != 0){
      						int ssize = tsize-tm;	// ssize = 2*size-2m
      						int ssize1= ssize+1;
      						float gnr2_r = gnr2[ssize];
      						float gnr2_i = gnr2[ssize1];
							//C[size-n][m].re += gnr2_r * tempr - gnr2_i * tempi;   // (a-ib)(c-id)=(ac-bd)-i(ad+bc)
						//C[size-n][m].im -= gnr2_i * tempr + gnr2_r * tempi;
      						C[(size-n)*(bw+1) + m].re += gnr2_r * tempr - gnr2_i * tempi;   // (a-ib)(c-id)=(ac-bd)-i(ad+bc)
      						C[(size-n)*(bw+1) + m].im -= gnr2_i * tempr + gnr2_r * tempi;
   						}
   						else{
						//C[size-n][m].re += *(gnr2)  * tempr - *(gnr2+1)* tempi;
						//C[size-n][m].im -= *(gnr2+1)* tempr + *(gnr2)  * tempi;
   							C[(size-n)*(bw+1) + m].re += *(gnr2)  * tempr - *(gnr2+1)* tempi;
   							C[(size-n)*(bw+1) + m].im -= *(gnr2+1)* tempr + *(gnr2)  * tempi;
   						}
   					}
   				}
   			}
   	    }

		// 2D inverse FFT to the correlation coefficients matrix c(phi,phi')
		//rfftwnd_one_complex_to_real(pinv, &C[0][0], &c[0][0]);

		rfftwnd_one_complex_to_real(pinv, C, c);

   		float tempr=0.0f, corre_fcs=-999.0f;
   		int n_best=0, m_best=0;
   		for(n=0;n<size;n++){			// move Tri_2D to Tri = c(phi,phi';rho)
   			float temp=0.0f;
   			for(m=0;m<size;m++){
			//temp=c[n][m];
   				temp=c[n*size + m];
   				if(temp>tempr) {
   					tempr=temp;
   					n_best=n;
   					m_best=m;
   				}
   			}
   		}
	//printf("rho=%d, maxcor_r=%f\n",p,tempr);
	/* **************************************************************
	*     Print results (since we just keep the largest corr value, *
	*     we don't need the loop for ifit)                          *
	*     V-> vector, distance is Rho                               *
	*     T = V-(comH-comL) = V+comL-comH 				*
	*     T is the vector by which the hi map has to be displaced   *
	*     H is the both Rot & Tran one; L is the only Rot one       *
	*     here, H is with (raw images); L is self (ref image)       *
	*****************************************************************
   		float Phi2= n_best*PI/bw;  // ming this is reference image rotation angle
   		float Phi = m_best*PI/bw;   // ming this is particle image rotation angle
   		float Vx,Vy, Tx, Ty;
   		if(neg_rt==1) {
   			Vx =  p*cos(Phi2);	//should use the angle of the centered one
   			Vy = -p*sin(Phi2);
   		}
   		else if(neg_rt==2) {
   			Vx =  p*cos(Phi);
   			Vy = -p*sin(Phi);
   		}
   		Tx = Vx + (floor(com_self_x+0.5f) - floor(com_with_x+0.5f));
   		Ty = Vy + (floor(com_self_y+0.5f) - floor(com_with_y+0.5f));
   		dx = -Tx;	// the Rota & Trans value (Tx,Ty, ang_keep) are for the projection image,
   		dy = -Ty;	// need to convert to raw image



		//**** Using FSC as the similarity criterion, 7/13/04 ****
   		float corre;


   		if(neg_rt==1){	// 9/6/04 ycong, fix projection image while move raw image, help to get better correlation
   			EMData *with_tmp = with->copy();
   			EMData *this_tmp = this->copy();
   			with_tmp->setRAlign((Phi2-Phi),0,0);
   			with_tmp->setTAlign(Tx,Ty,0);

   			with_tmp->rotateAndTranslate();

    	//ming comment,	corre=with_tmp->fscmp(this);
   			corre=with_tmp->cmp("frc",this,{});
   			printf("corre=%f\n",corre);
   			if(with_tmp){
   				delete with_tmp; with_tmp=0;
   			}
   			if(this_tmp){
   				delete this_tmp; this_tmp=0;
   			}
   		}
   		else if(neg_rt==2) {
   			EMData *this_tmp = this->copy();
   			this_tmp->setRAlign(-(Phi2-Phi),0,0);
   			this_tmp->setTAlign(dx,dy,0);

   			this_tmp->rotateAndTranslate();		// rotateAndTranslate: Trans first then rotation
    	//ming comment	corre=this_tmp->fscmp(with);

   			this_tmp->cmp(cmp_name,with,cmp_params);
   			printf("p_max=%d, ok??\n",p_max);
   			delete this_tmp;
   		}
   		if(corre>=corre_fcs) {
  			corre_fcs=corre;
   			result[0+5*p] = float(p);	// rho
   			result[1+5*p] = corre;		// correlation_fcs
   			result[2+5*p] = Phi2-Phi;	// rotation angle
   			result[3+5*p] = Tx;		// Translation_x
   			result[4+5*p] = Ty;		// Translation_y
		//printf("##rho=%d, corr=(%f, %f), rot=%f, Trans=(%f, %f)\n",p,tempr,corre,(Phi2-Phi)*180.0/PI,Tx,Ty);
   		}
   		maxcor[p]=corre_fcs;               		//  maximum correlation for current rho
   		if(corre_fcs>maxcor_sofar) {
   			maxcor_sofar=corre_fcs;   		// max correlation up to current rho
   			rho_best=p;				// the rho value with maxinum correlation value
   		}

		// if ok to exit early
   		if(p>=4){
   			if(   maxcor[p  ] < maxcor[p-1] && maxcor[p-1] < maxcor[p-2]
		   && maxcor[p-2] < maxcor[p-3] && maxcor[p-3] < maxcor[p-4]){
			//rhomax = rhomin + (p-4)*tstep;     // new rhomax
			loop_rho=1;
			break;           //exit p loop
   			}
   		}
	}

	if(loop_rho == 1)
		p=p+1;
	int rb5=5*rho_best;
//printf("****** The best fit:  rho_max/best=(%d/%d), corr=%f, rot_trans = (%f; %f, %f) *******\n", p, int(result[0+rb5]), result[1+rb5],result[2+rb5]*180/PI,result[3+rb5],result[4+rb5]);
	float rho_keep = result[0+rb5];
	float fsc      = result[1+rb5];
	float ang_keep = result[2+rb5];
	float Tx       = result[3+rb5];
	float Ty       = result[4+rb5];
	//printf("****** The best fit:  rho_max/best=(%d/%d), corr=%f, rot_trans = (%f; %f, %f) *******\n", p, int(rho_keep), fsc,ang_keep*180/PI,Tx,Ty);


// ming comment,	doneData();
	delete[] gnr2;
	delete[] maxcor;
	delete[] result;
	//delete selfpcsfft;

	delete c;
	c = 0;
	delete C;
	C = 0;
	rfftwnd_destroy_plan(pinv);
	//self_sampl_fft=NULL;

	if (!neg_rt) {
		this->setRAlign(ang_keep,0,0);
		this->setTAlign(0,0,0);
		this->rotateAndTranslate();
		this->setTAlign(Tx,Ty,0);
		this->fastTranslate();
		dx = Tx;
		dy = Ty;
	}
	else if(neg_rt==1){		//this output style is for matching EMAN data style and "refine" case
		this->setRAlign(-ang_keep,0,0);
		this->setTAlign(0,0,0);
		this->rotateAndTranslate();
		//this->writeIMAGIC("rot_this.hed",-1);

		dx = -Tx;
		dy = -Ty;
		this->setTAlign(dx,dy,0);
		this->fastTranslate();
		//for matching the following "classesbymra.C" coordinate convert style
        float cdx= (-Tx)*cos(-ang_keep)+(-Ty)*sin(-ang_keep);
       	float cdy=-(-Tx)*sin(-ang_keep)+(-Ty)*cos(-ang_keep);     //this two lines are rotation matrix

       	setRAlign(-ang_keep,0,0);
       	setTAlign(cdx,cdy,0);
       	//printf("%f, %f  \n", cdx, cdy);
	}
	else if(neg_rt==2){             //==2: when Projection image is the one involve (T & R)
		dx = -Tx;		// the Rota & Trans value (Tx,Ty, ang_keep) are for the projection image,
		dy = -Ty;		// need to convert to raw image
		this->setRAlign(-ang_keep,0,0);
		this->setTAlign(dx,dy,0);
		this->rotateAndTranslate();             // rotateAndTranslate: Trans first then rotation
	}
	return fsc;     // return the fsc coefficients
}
******************/



//using fscmp() function as similarity criterion, works pretty good, check more about the "ALTCMP" value, ycong 6/26/04
/*******************ming
EMData *EMData::frm2dAlignFlip(EMData *drf, float &dot_frm, EMData *with, int usedot, float *raw_hhat, float rmax, float MAXR_f, float tstep, float angstep, float r_size, float &com_self_x, float &com_self_y, float &com_with_x, float &com_with_y, int neg_rt)
{
	float dot_frm0=0.0f, dot_frm1=0.0f;
	int size=(int)(angstep+0.5f);
	int MAXR = (int)(MAXR_f+0.5f);
	EMData *da_nFlip=0, *da_yFlip=0;
	//EMData *selfpcs=0;

	for (int iFlip=0;iFlip<=1;iFlip++){
		EMData *dr_frm=0;
		if (iFlip==0) { dr_frm=this->copy(); da_nFlip=this->copy();}     //delete dr_frm somewhere
		else          { dr_frm=drf->copy();  da_yFlip=drf->copy();}

		if(iFlip==1) com_self_x=dr_frm->get_xsize()-com_self_x;   //ming   // image mirror about Y axis, so y keeps the same

		dx=-(com_self_x - dr_frm->get_xsize()/2); //ming
		dy=-(com_self_y - dr_frm->get_ysize()/2); //ming

		dr_frm->setTAlign(dx,dy,0);


		dr_frm->fastTranslate(); //ming, something wrong in this function
		//dr_frm->writeIMAGIC("dr_frm.hed",-1);



		EMData *selfpcs = dr_frm->unwrap_largerR(0,MAXR,size, rmax);
		EMData *selfpcsfft = selfpcs->oneDfftPolar(size, rmax, MAXR);	// 1DFFT of polar sampling
		delete dr_frm;
		delete selfpcs;

		if(iFlip==0)
			dot_frm0=da_nFlip->frm2dAlign(with, raw_hhat, selfpcsfft, rmax, MAXR, tstep, angstep, r_size, com_self_x, com_self_y, com_with_x, com_with_y,neg_rt);
		else
			dot_frm1=da_yFlip->frm2dAlign(with, raw_hhat, selfpcsfft, rmax, MAXR, tstep, angstep, r_size, com_self_x, com_self_y, com_with_x, com_with_y,neg_rt);

		delete selfpcsfft;
		// selfpcsfft was deleted at frm2dAlign();
	}

	float d1,d2;
	if (usedot) {
		d1=da_nFlip->dot(with);
		d2=da_yFlip->dot(with);
		if (usedot==2) printf("%f vs %f  (%1.1f, %1.1f  %1.2f) (%1.1f, %1.1f  %1.2f)\n",
		d1,d2,da_nFlip->Dx(),da_nFlip->Dy(),da_nFlip->alt()*180./PI,da_yFlip->Dx(),da_yFlip->Dy(),da_yFlip->alt()*180./PI);
	}
	else {
//#ifdef ALTCMP
		//ming comment, d1=da_nFlip->fscmp(with);
		d1=da_nFlip->cmp("frc",with,{});
		//ming comment, d2=da_yFlip->fscmp(with);
		d2=da_yFlip->cmp("frc",with,{});
		//printf("if ALTCMP d1,d2=%f, %f",d1,d2);
//#else
//		d1=2.0-da_nFlip->lcmp(with,1);
//		d2=2.0-da_yFlip->lcmp(with,1);
//		printf("else d1,d2=%f, %f",d1,d2);
//#endif
	}

	if(d1 >= d2) {
		da_nFlip->setFlipped(0);
		dot_frm=dot_frm0;              // cross corallation coefficient
		delete da_yFlip;
		return da_nFlip;
	}
	else {
		da_yFlip->setFlipped(1);
		dot_frm=dot_frm1;
		delete da_nFlip;
		return da_yFlip;
	}
}
*****ming*****/


// ~~~~~ determine exactly COM of raw image,  Yao Cong ~~~~
void EMData::exact_COM(int intonly,float *Dens_r, float &cx, float &cy){

	int r;
	int flag_com;           /* signals the finding of the best comL */
	int icurr, jcurr;       /* current values of the position in cormax */
	float comLx0, comLy0;   /* initial comL */
	float cormax[12][12];   /* max correlation as a function of comL */
	float ccLx, ccLy;       /* corrections to comL */
	float comgfac;          /* factor for comL correction (see details below) */
	float corM;             /* maximum correlation as a function of comL */
	float temp=0.0f, tempi, tempj;

	norm_frm2d();
	//this->writeIMAGIC("exax.hed",-1);

	float *d=get_data(); // ming change getData() to get_data()

	for(int i=-5;i<=5;i++)                         /* initialize cormax */
		for(int j=-5;j<=5;j++) {
			cormax[i+5][j+5]=-1.0;     /* max correlation as a function of comL */
		}

	comgfac=0.02*(get_xsize()+get_ysize()); // ming           /* ratio between grid spacing of cormax to that of phiL */
	ccLx=0.0f; ccLy=0.0f;
	icurr=0;jcurr=0;
	comLx0=cx; comLy0=cy;                      /* save initial comL */
	flag_com=0;

	do{
		cx=comLx0+ccLx; cy=comLy0+ccLy;
		//printf("cx=%f, %f\n", cx,cy);

		float corr_ab=0.0f;
		for (int i=0; i<get_xsize(); i++){//ming
			for (int j=0; j<get_ysize(); j++){//ming
				//float r_2;//ming
				float r_2 = std::sqrt((float)((i-cx)*(i-cx) + (j-cy)*(j-cy)));
				int r_1 = floor(r_2);
				int r_3 = r_1 + 1;
				float Int_g2 = (r_2-r_1)*Dens_r[r_3] + (r_3-r_2)*Dens_r[r_1];  // from g1 & g3, get g2
				corr_ab += d[i+j*get_xsize()] * Int_g2;       //ming               // c(alpha,beta)
      			}
		}
		//printf("corr_ab=%f \n", corr_ab);

		if(flag_com==0){
			//printf("curr = %d %d %f\n", icurr, jcurr,corr_ab);
			cormax[icurr+5][jcurr+5]=corr_ab;
			corM=-1.0;
			for(int i=-5;i<=5;i++)                    // look for the maximum of cormax
				for(int j=-5;j<=5;j++)
					if(cormax[i+5][j+5]>corM){
						corM=cormax[i+5][j+5];
						icurr=i; jcurr=j;
					}

			if(cormax[icurr+4][jcurr+5]<0){
				icurr--;
				ccLx=icurr*comgfac; ccLy=jcurr*comgfac;
			}
			else if(cormax[icurr+6][jcurr+5]<0){
				icurr++;
				ccLx=icurr*comgfac; ccLy=jcurr*comgfac;
			}
			else if(cormax[icurr+5][jcurr+4]<0){
				jcurr--;
				ccLx=icurr*comgfac; ccLy=jcurr*comgfac;
			}
			else if(cormax[icurr+5][jcurr+6]<0){
				jcurr++;
				ccLx=icurr*comgfac; ccLy=jcurr*comgfac;
			}
			// determine point based on parabola
			else{
				float f000=cormax[icurr+5][jcurr+5];
				float f100=cormax[icurr+4][jcurr+5];
				float f200=cormax[icurr+6][jcurr+5];
				float f010=cormax[icurr+5][jcurr+4];
				float f020=cormax[icurr+5][jcurr+6];

				if((temp=f100+f200-2*f000)==0)
					tempi=0.0;
				else
					tempi=(f100-f200)/(2*temp);
				ccLx=(icurr+tempi)*comgfac;

				if((temp=f010+f020-2*f000)==0)
					tempj=0.0;
				else
					tempj=(f010-f020)/(2*temp);
				ccLy=(jcurr+tempj)*comgfac;

				flag_com=1;
			}
		}
		else
			break;
	}while(1);
	//ming comment, doneData();

	if (intonly >0) {
		cx = floor(cx+.5);
		cy = floor(cy+.5);
	}

	dx=-(cx-get_xsize()/2);//ming
	dy=-(cy-get_ysize()/2);//ming

	parent=NULL;                   // parent=NULL -> parent image is allowed to be changed inplace
	fastTranslate();               // default is 1 -> in-place trans;

}

float* EMData::makeFRM2D_H_hat(float rhomax, float rhomin, float rmaxH, int size, float tstep, float r_size_flo)
{
	int bw=size/2;
	unsigned long tsize2=(size+2)*size;
	unsigned long tsize=2*size;
	unsigned long size2=size*size;
	unsigned long ind1=0, ind2=0, ind3=0, ind4=0, ind41=0;
	float pi2=2.0*PI, r2;
	FILE *outs;
	int r_size=(int)(r_size_flo);
	int p_max=r_size-(int)(rmaxH);
	//printf("nx=%f, %f %f  %d, %f %f\n", rhomax, rhomin, rmaxH, size, tstep, r_size);

	//fftw_complex gnr2_in[size], gnr2[size];
	fftw_complex * gnr2_in = new fftw_complex[size];
	fftw_complex * gnr2 = new fftw_complex[size];
	fftw_plan planc_fftw;			// plan for 1D FFTW (Forward), complex->complex  to get gnr2
	planc_fftw = fftw_create_plan(size, FFTW_FORWARD, FFTW_ESTIMATE );    // out-of-plalce transform, see FFTW menu p.3-4

 	float *frm2dhhat=0;
	//cout << "size:" << 10*(r_size+1)*bw*size*2 * sizeof(float) << endl;
	if( (frm2dhhat = (float *) malloc( (p_max+1)*(r_size+2)*bw*size*2 * sizeof(float)) ) == NULL ) {
	//if( (frm2dhhat = malloc(sizeof(float)*(int(rhomax/tstep+0.5f)+1)*(r_size+1)*bw*size*2)) == NULL ) {
		cout <<"Error in allocating memory 13. \n";
		exit(1);
	}

	float *sb=0, *cb=0;		// sin(beta) and cos(beta) for get h_hat, 300>size
	if( (sb = new float[size]) == NULL || (cb = new float[size]) == NULL) {
		cout <<"can't allocate more memory, terminating. \n";
		exit(1);
	}
	for(int i=0;i<size;i++) {        // beta sampling, to calculate beta' and r'
		float beta=i*PI/bw;
	   	sb[i]=sin(beta);
	   	cb[i]=cos(beta);
	}

	float *sampl_fft = get_data(); // ming change getDataRO() to get_data

	//compute H_hat
	//for(pi=0; pi<=int((rhomax-rhomin)/tstep+0.5f); pi++){                 // loop for rho->p
	for(int p=0; p<=p_max; p++){
		ind1=p*r_size*bw;
		//printf("p=%d, ind1=%d\n", p, ind1);

    	float pp2 = p*p;
		int r_va1 = p-rmaxH;
		if(r_va1 < 1) r_va1 = 1;
		int r_va2 = p+rmaxH;

   		for(int n=0;n<bw;n++){         /* loop for n */
       		int tn=n*2;
       		int tn1=tn+1;
    		ind2 = ind1+n;

      			//for(int r=1;r<=(r_size-1);r++) {     // add terms corresponding to r to each integral
			//for(int r=1;r<=55;r++) {
			for(int r=r_va1;r<=r_va2;r++) {
				ind3 = (ind2+r*bw)*size;
        		float rr2 = r*r;
				float rp2 = r*p;
       			for(int i=0;i<size;i++){                            // beta sampling, to get beta' and r'
       				if(p==0) r2 = (float) r;
       				else     r2 = std::sqrt((float)(rr2+pp2-2.0*rp2*cb[i]));   // r2->r'
       		 		int r1=floor(r2+0.5);                        // for computing gn(r')
       				if(r1>rmaxH){
       					gnr2_in[i].re  = 0.0f;
       					gnr2_in[i].im  = 0.0f;
       				}
       				else{
       					float gn_r = sampl_fft[tn +r1*(size+2)];           // real part of gn(r')
       					float gn_i = sampl_fft[tn1+r1*(size+2)];           // imaginary part of gn(r')
						float sb2, cb2, cn, sn;
						if(n != 0){
							if(r2 != 0.0){
								if(p==0) {
									sb2=sb[i];
									cb2=cb[i];
								}
								else{
									sb2=r*sb[i]/r2;
									cb2=(r*cb[i]-p)/r2;
								}
							}
        					else{
								sb2=0.0;
								cb2=1.0;
							}
        					if(sb2>1.0) sb2= 1.0f;
        					if(sb2<-1.0)sb2=-1.0f;
        					if(cb2>1.0) cb2= 1.0f;
        					if(cb2<-1.0)cb2=-1.0f;
        					float beta2=atan2(sb2,cb2);
        					if(beta2<0.0) beta2 += pi2;
        					float nb2=n*beta2;
        					cn=cos(nb2);
							sn=sin(nb2);
						}
        				else{
							cn=1.0f; sn=0.0f;
						}

						gnr2_in[i].re =   cn*gn_r - sn*gn_i;	     // exp(-i*n*beta')*gn(r')
						gnr2_in[i].im = -(cn*gn_i + sn*gn_r);
        			}
        		}  /* for i */

        		fftw_one(planc_fftw, gnr2_in, gnr2);

				for(int m=0;m<size;m++){                                     // store hat{h(n,r,p)}(m)
					ind4 = (ind3+m)*2;
					ind41 = ind4+1;
					frm2dhhat[ind4]  = gnr2[m].re;
					frm2dhhat[ind41] = gnr2[m].im;
				}   /* m */
			}   /* r */
		}   /* n */
	}  /* p */

	// ming comment, doneData();
	delete[] sb;
	delete[] cb;
	delete [] gnr2_in;
	gnr2_in = 0;
	delete [] gnr2;
	gnr2 = 0;
	fftw_destroy_plan(planc_fftw);
	return frm2dhhat;
}




//1D FFT of polar coordinate sampling
// r<rmax, do the 1D FFT, otherwise padding with 0
// Yao Cong   5/31/04
EMData *EMData::oneDfftPolar(int size, float rmax, float MAXR){		// sent MAXR value here later!!

    fftw_real * in = new fftw_real[size]; //ming
	fftw_real * out = new fftw_real[size]; // ming
	rfftw_plan p; //ming
	p = rfftw_create_plan(size, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

	float *pcs=get_data(); // ming change getDataRO() to get_data(), this read the image information for the pointer which call this function

	EMData *imagepcsfft = new EMData; // ming, create a new image object pointer *imagepcsfft
	imagepcsfft->set_size((size+2), (int)MAXR+1, 1);            //"setSize" comes together with "new EMData" , ming, define the new image size
	float *d=imagepcsfft->get_data(); // ming change getData() to get_data(), read the original image data to pointer d
	for(int row=0; row<=(int)MAXR; row++){
		if(row<=(int)rmax) {
			for(int i=0; i<size;i++){
				in[i] = pcs[i+row*size]; // ming
			}
			rfftw_one(p,in,out);	// 1D real FFTW, details see FFTW menu p.6-7 //ming


			d[0+row*(size+2)] = out[0];//ming
			d[1+row*(size+2)] = 0.0;
			for (int k=1; k<(size+1)/2;k++){
				 d[2*k+row*(size+2)]   = out[k];
				 d[2*k+1+row*(size+2)] = out[size-k];
			}
			if(size%2 == 0){                                     // size is even
				d[size  +row*(size+2)] = out[size/2];
				d[size+1+row*(size+2)] = 0.0;
			}
		}
		else{
			for(int j=0;j<size+2;j++)
				d[j+row*(size+2)] = 0.0;
		}
	}

	delete in;
	in = 0;
	delete out;
	out = 0;

	rfftw_destroy_plan(p);
	// ming comment, doneData();
	imagepcsfft->update(); // ming change doneData() to update()
	//delete this;
	return imagepcsfft;
}



/*****************ming comment*
int EMData::setSize(int x,int y, int z) {

	if (waitReady(0)) return -1;

	if (x<=0) x=1;
	if (y<=0) y=1;
	if (z<=0) z=1;
	if (y==1&&z!=1) { y=z; z=1; }
	if (x==1&&y!=1) { x=y; y=1; }
	nx=x;
	ny=y;
	nz=z;
	if (flags&EMDATA_NODATA) {
		if (rdata) free(rdata);
		rdata=NULL;
	}
	//else rdata=(float *)Srealloc(rdata,x*y*z*sizeof(float),1);
	else rdata=(float *)realloc(rdata,x*y*z*sizeof(float));
	update();
	flags&= ~EMDATA_SWAPPED;
	if (supp) { free(supp); supp=NULL; }
	return 0;
}
ming comment ***********************8*/


void EMData::norm_frm2d() {

	unsigned long i;
	unsigned int nvox=get_xsize()*get_ysize();//ming
	float maxmap=0.0f, minmap=0.0f;
	float temp=0.0f, diff_den=0.0f, norm=0.0f;
	float cut_off_va =0.0f;
	float* data;

	//writeIMAGIC("norm_in",-1);
	//ming comment, data = getData();
	data=get_data();

	// rescale the map
	maxmap=-1000000.0f;
	minmap=1000000.0f;
	for (i=0;i<nvox;i++){
		if(data[i]>maxmap) maxmap=data[i];
		if(data[i]<minmap) minmap=data[i];
	}
	diff_den = maxmap-minmap;
	//printf("max, min =%f, %f\n", maxmap, minmap);

	for(i=0;i<nvox;i++) {
		temp = (data[i]-minmap)/diff_den;
		if(cut_off_va != 0.0) {               // cut off the lowerset ?% noisy information
	 		if(temp < cut_off_va)
	   			data[i] = 0.0f;                   // set the empty part density=0.0
	 		else
	   			data[i] = temp-cut_off_va;
		}
		else{
			data[i] = temp;
			//printf("%3.1f ",data[i]);
		}
	}

	// normalize the image
	for(i=0;i<nvox;i++) {
		temp=data[i];
		norm += temp*temp;
	}
	//norm = sqrt(norm); // ming, temporary comment it, after figure out uncomment it
	for(i=0;i<nvox;i++){
		data[i] /= norm;                      //  y' = y/norm(y)
		//printf("%6.4f ",data[i]);
	}

	// ming comment, doneData();
	update();
	//writeIMAGIC("norm_out",-1);

}


/**** ming begin***/
EMData *EMData::unwrap_largerR(int r1,int r2,int xs, float rmax_f) {
	//int i,r;
	float *d,*dd;
	int do360=2;
	int rmax = (int)(rmax_f+0.5f);
	norm_frm2d();
	//ming comment, d=getDataRO();
	d=get_data(); // ming change getaDataRO() to eman2 get_data()
	//writeIMAGIC("for_samp.hed",-1);

	if (xs<1) {
		xs = (int) floor(do360*PI*get_ysize()/4); // ming
		xs=Util::calc_best_fft_size(xs); // ming
	}
	if (r1<0) r1=0;
	float maxext=ceil(0.6*std::sqrt((double)(get_xsize()*get_xsize()+get_ysize()*get_ysize())));// ming add std::
	//printf("maxext=%f\n",maxext);
	if (r2<r1) r2=(int)maxext;
	EMData *ret = new EMData;	     // here r2=MAXR, so ret size is (size*(MAXR+1))

	ret->set_size(xs,r2+1,1);    //ming, problem is here         //EMAN style: ret->setSize(xs,ySize()/2-5,1);
	//dd=ret->getData();  // something wrong with this part

	dd=ret->get_data();

	for (int i=0; i<xs; i++) {
		float si=sin(i*PI*2/xs);
		float co=cos(i*PI*2/xs);
		for (int r=0; r<=maxext; r++) {
			float x=(r+r1)*co+get_xsize()/2; // ming
			float y=(r+r1)*si+get_ysize()/2; // ming
			if(x<0.0 || x>=get_xsize()-1.0 || y<0.0 || y>=get_ysize()-1.0 || r>rmax){    //Ming , ~~~~ rmax need pass here
				for(;r<=r2;r++)                                   // here r2=MAXR
					dd[i+r*xs]=0.0;
        		break;
		    }
			int x1=floor(x);
			int y1=floor(y);
			float t=x-x1;
			float u=y-y1;
			float f11= d[x1+y1*get_xsize()]; // ming
			float f21= d[(x1+1)+y1*get_xsize()]; // ming
			float f12= d[x1+(y1+1)*get_xsize()]; // ming
			float f22= d[(x1+1)+(y1+1)*get_xsize()]; // ming
			dd[i+r*xs] = (1-t)*(1-u)*f11+t*(1-u)*f21+t*u*f22+(1-t)*u*f12;
		}
	}

	update();
	// ming comment, doneData();

	ret->update(); // ming change doneData() to update() since the data is changed

	return ret;
}
/***ming end****************/



int EMData::swapin() {
	//const int MAXPATHLEN=1024; // MING
	char s[MAXPATHLEN];
	FILE *in;
	//const char  *ERR_EXIT, *ERR_CRIT;// ming add
	if (!(flags&EMDATA_SWAPPED)) return 0;
   // int tmp;// ming add
	sprintf(s,"/usr/tmp/tmp.%d",tmp);
	in=fopen(s,"rb");
	if (!in) { Util::error("ERROR: Swap file missing!",ERR_EXIT); return -1; }
	//rdata=(float *)Smalloc(nx*ny*nz*sizeof(float),1);
	rdata=(float *)malloc(nx*ny*nz*sizeof(float));


	if (!rdata) { Util::error("Cannot swap in, not enough memory!",ERR_CRIT); return -1; }
	if(fread(rdata,sizeof(float)*nx,ny*nz,in)!=ny*nz) {
		Util::error("ERROR: Swap file incomplete!",ERR_EXIT);
		return -1;
	}
	flags&=~EMDATA_SWAPPED;
	return 0;
}

int EMData::waitReady(int ro)
{
/*
if (ro) {
while (flags&EMDATA_BUSY) sleep(1);
}
else {
while ((flags&EMDATA_BUSY) || (rocount)) sleep(1);
}
	*/

	// for debugging only
	if (ro) {
		while (flags&EMDATA_BUSY) {
			printf("Busy 1 (%p)\n",this);
#ifdef ALPHA
			sleep(1);
#else
			usleep(100);
#endif
		}
	}
	else {
		if ((flags&EMDATA_BUSY) || (rocount)) {
		printf("Busy 0.%d (%p),%d,%d\n",rocount,this,flags,EMDATA_BUSY);
		while ((flags&EMDATA_BUSY) || (rocount)) {
#ifdef ALPHA
			sleep(1);

#else
		 usleep(100);
#endif
			}
			printf("Busy done 0.%d (%p)\n",rocount,this);
		}
	}

	return 0;
}

/********* ming begin
float *EMData::getData() {	// This returns a pointer to the raw float data
	//if (flags&EMDATA_BUSY || rocount>0 || !rdata) return NULL;

	waitReady(0);

	if (flags&EMDATA_SWAPPED) swapin();
	flags|=EMDATA_BUSY;
	return rdata;
}


void EMData::doneData() // MUST be called when the data is no longer being
// modified. Another getData while the data is
// out will return NULL
{
	if (flags&EMDATA_BUSY) {
		flags&=(~EMDATA_BUSY);
		update();
	}
	else {
		rocount--;
		if (rocount<0) rocount=0;
	}

	return;
}


void EMData::normalize() {
	int i,di=1;
	float f,*data;
	float nor,s;

	if ((flags&EMDATA_COMPLEX) && !(flags&EMDATA_RI)) di=2;

	s=Sigma();			// need pointer for goodf


	if ( !Util::goodf(&s)) { zero(); return; }// ming


	if (s==0 ) return;
	// ming comment, data=getData();
    data=get_data(); // ming change getData() to get_data()

	// this normalizes to sigma excluding pixels which are exactly zero
	//s=Sigma2();
	//if (s>0 && goodf(&s)) nor=s;
	//else nor=Sigma();

	// The above may be risky, we'll return to the standard normalization
	nor=Sigma();

	//printf("%f %f \n",nor,Sigma());

	for (i=0; i<nx*ny*nz; i+=di) {
		data[i]=(data[i]-Mean())/nor;
	}

	// ming comment, doneData();
	update();
}



void EMData::COM(int intonly, float &xm, float &ym) {
	float m=0.0f;
	float zm=0.0f;
	int i,j,k;


	//norm_frm2d();
//ming comment,	normalize();      //EMAN way to do normalize

	//this->writeIMAGIC("test-print.hed",-1);

    //ming comment, getData();
	get_data(); //ming change getData() to get_data()

	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			for (k=0; k<nz; k++) {
				float tempt = rdata[i+nx*j+k*nx*ny];
				//printf("COM  %f\n",tempt);
			    if (tempt<Sigma()*0.75+Mean()) continue; // may help since protein is denser
				xm+=(float)i*tempt;
				ym+=(float)j*tempt;
				zm+=(float)k*tempt;
				m+=tempt;
				//printf("%f  (%f, %f, %f) (%d %d %d) \n",tempt,xm,ym,zm,i,j,k);
			}
		}
	}
	// ming comment, doneData();
	xm/=m;
	ym/=m;
	zm/=m;

	//(intonly == 2) { xm=150; ym=153;zm=0;}
	//printf("%f  %f  %f\n",xm,ym,zm);
	if (intonly >0) {
		xm=floor(xm+.5);
		ym=floor(ym+.5);
		zm=floor(zm+.5);

		dx=-(xm-nx/2);
		dy=-(ym-ny/2);
		if (nz>1) dz=-(zm-nz/2);
	}
	else {
		dx=-(xm-nx/2);
		dy=-(ym-ny/2);
		if (nz>1) dz=-(zm-nz/2);
	}
	//printf("~~~ dx, dy=%5.2f, %5.2f, %f, %f--> \n",dx,dy, xm,ym);
	parent=NULL;                   // parent=NULL -> parent image is allowed to be changed inplace

	if(intonly!=3) fastTranslate();               // default is 1 -> in-place trans;
}

// ~~~~~ determine exactly COM of raw image,  Yao Cong ~~~~
*************ming end *********/
