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

#ifndef eman__emdatacuda_h__
#define eman__emdatacuda_h__ 1

#ifdef EMAN2_USING_CUDA

public:
	/** Called externally if the associated cuda array is freed
	 */
	inline void reset_cuda_array_handle() { cuda_array_handle = -1; }
	

	inline void bind_cuda_array() {
		update_stat();
		if (cuda_array_handle == -1) {
			cuda_array_handle = get_cuda_array_handle(rdata,nx,ny,nz,this);
		}
		bind_cuda_texture(cuda_array_handle);
	}
	// This should never be set by anything other than something that knows what it's doing
	inline void set_cuda_array_handle(const int idx) { cuda_array_handle = idx; }
	
private:
	
	
	void free_cuda_array() const {
		if (cuda_array_handle != -1){
			delete_cuda_array(cuda_array_handle);
			cuda_array_handle = -1;
		}
	}
	
	void free_cuda_memory() const {
		if ( cuda_rdata != 0) {
			delete_cuda_memory(cuda_rdata);
			cuda_rdata = 0;
		}
	}
	mutable int cuda_array_handle;
	mutable float* cuda_rdata;

#endif // EMAN2_USING_CUDA
	
#endif //eman__emdatacuda_h__ 1

