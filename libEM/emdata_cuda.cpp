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
#include "cuda/cuda_processor.h"



using namespace EMAN;

float* EMData::get_cuda_data() const {
	device_init();
	size_t num_bytes = nx*ny*nz*sizeof(float);
	if (cuda_rdata == 0) {
		cudaMalloc((void**)&cuda_rdata,num_bytes);
	}
	
	if (rdata != 0 && (EMDATA_GPU_NEEDS_UPDATE & flags)) { // If rdata is zero it means we're working exclusively on the GPU
		cudaMemcpy(cuda_rdata,rdata,num_bytes,cudaMemcpyHostToDevice);
	}
	
	flags &= ~EMDATA_GPU_NEEDS_UPDATE;
	return cuda_rdata;
}

void EMData::bind_cuda_array() {
	if (cuda_array_handle==-1|| EMDATA_GPU_RO_NEEDS_UPDATE & flags) {
		if (cuda_array_handle==-1) delete_cuda_array(cuda_array_handle);
		cuda_array_handle = get_cuda_array_handle(get_cuda_data(),nx,ny,nz,this);
		flags &= ~EMDATA_GPU_RO_NEEDS_UPDATE;
	}
	bind_cuda_texture(cuda_array_handle);
}


void EMData::gpu_update() {
	flags |= EMDATA_NEEDUPD | EMDATA_CPU_NEEDS_UPDATE | EMDATA_GPU_RO_NEEDS_UPDATE;
}

void EMData::mult_cuda(const float& val) {
	Dict d("scale",(float)val);
	process_inplace("cuda.math.mult",d);
}

EMData* EMData::calc_ccf_cuda( EMData*  image ) {
	
	EMData* tmp;
	if (is_complex()) {
		tmp = new EMData(*this);
	} else {
		tmp = do_fft_cuda();
	}
	
	Dict d;
	EMData* with = 0;
	if (!image->is_complex()) {
		with = image->do_fft_cuda();
		d["with"] = (EMData*) with;
	} else {
		d["with"] = (EMData*)image;
	}
	
	tmp->process_inplace("cuda.correlate",d);
	
	if (with != 0) {
		delete with;
		with = 0;
	}

	EMData* soln = tmp->do_ift_cuda();
	soln->gpu_update();
	delete tmp;
	tmp = 0;
	
	return soln;
}

void EMData::free_cuda_array() const {
	if (cuda_array_handle != -1){
		delete_cuda_array(cuda_array_handle);
		cuda_array_handle = -1;
	}
}

void EMData::free_cuda_memory() const {
	if ( cuda_rdata != 0) {
// 		cout << "Deleting cuda_rdata" << endl;
		cudaFree(cuda_rdata);        
		cudaThreadSynchronize();
		cuda_rdata = 0;
	}
}



#endif //EMAN2_USING_CUDA