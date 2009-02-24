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



#ifdef EMAN2_USING_CUDA

#include "emdata.h"
#include <cuda_runtime_api.h>
#include "cuda/cuda_util.h"


using namespace EMAN;

float* EMData::get_cuda_data() const {
	device_init();
	update_stat(); // This will delete the cuda pointer, i.e. if it's out of date
	if (cuda_rdata == 0) {
		size_t num_bytes = nx*ny*nz*sizeof(float);
		cudaMalloc((void**)&cuda_rdata,num_bytes);
		if (rdata != 0) { // If rdata is zero it means we're working exclusively on the GPU
			cudaMemcpy(cuda_rdata,rdata,num_bytes,cudaMemcpyHostToDevice);
		}
	}
	return cuda_rdata;
}


void EMData::free_cuda_array() const {
	if (cuda_array_handle != -1){
		delete_cuda_array(cuda_array_handle);
		cuda_array_handle = -1;
	}
}

void EMData::free_cuda_memory() const {
	if ( cuda_rdata != 0) {
		cudaFree(cuda_rdata);
		cuda_rdata = 0;
	}
}

#endif //EMAN2_USING_CUDA