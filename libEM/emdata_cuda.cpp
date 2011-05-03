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
#include "exception.h"
#include <cuda_runtime_api.h>
#include <driver_functions.h>
#include <cuda.h>
#include <cuda/cuda_util.h>
#include <cuda/cuda_emfft.h>

using namespace EMAN;

// Static init
const EMData* EMData::firstinlist = 0;
const EMData* EMData::lastinlist = 0;
int EMData::memused = 0;
int EMData::fudgemem = 1.024E8; //let's leave 10 MB of 'fudge' memory on the device
int EMData::mempoolused = -1;
int EMData::mempoolarraysize = 0;
bool EMData::usemempoolswitch = false;
bool EMData::usecuda = (getenv("EMANUSECUDA") == NULL || atoi(getenv("EMANUSECUDA")) ) ? 1 : 0;
float* EMData::mempool[] = {0};

bool EMData::copy_to_cuda_keepcpu() const
{
	//cout << "copying from host to device RW" << " " << num_bytes << endl;
	if(rw_alloc()) {
		memused += num_bytes;	
		cudaError_t error = cudaMemcpy(cudarwdata,rdata,num_bytes,cudaMemcpyHostToDevice);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaMemcpy (device to host) failed:" + string(cudaGetErrorString(error)));
	}else{return false;}
	
	return true;
}

bool EMData::copy_to_cuda()
{
	//cout << "copying from host to device RW" << " " << num_bytes << endl;
	if(rw_alloc()) {
		memused += num_bytes;
		cudaError_t error = cudaMemcpy(cudarwdata,rdata,num_bytes,cudaMemcpyHostToDevice);
		if ( error != cudaSuccess) {
			//cout << rdata << " " << cudarwdata << endl;
			throw UnexpectedBehaviorException( "CudaMemcpy (host to device) failed:" + string(cudaGetErrorString(error)));
		}
	}else{return false;}
	//Temporaly disabled, causes LOTS of bugs......
	//free_rdata(); //we have the data on either the host or device, not both (prevents concurrency issues)
	
	return true;
}

bool EMData::copy_to_cudaro() const
{
	
	//cout << "copying from host to device RO" << " " << num_bytes << endl;
	if(ro_alloc()) {
		memused += num_bytes;
		copy_to_array(rdata, cudarodata, nx, ny, nz, cudaMemcpyHostToDevice);
	}else{return false;}
	
	return true;
}

bool EMData::rw_alloc() const
{
	//cout << "rw_alloc" << endl;
	num_bytes = nxyz*sizeof(float);
	if(cudarwdata){return true;} // already exists
	//use the mempool if available
	if(usemempoolswitch && mempoolused > -1 && mempoolarraysize >= int(num_bytes))
	{
		cudarwdata = mempool[mempoolused];
		mempool[mempoolused] = 0;
		mempoolused--;
		return true;
	}	
	if(!freeup_devicemem(num_bytes)){return false;}
	cudaError_t error = cudaMalloc((void**)&cudarwdata,num_bytes);
	if ( error != cudaSuccess){return false;}
	if(!cudarodata){addtolist();}
	//cout << "rw alloc finish" << endl;
	return true;
}

bool EMData::ro_alloc() const
{
	//cout << "ro_alloc" << endl;
	num_bytes = nxyz*sizeof(float);
	if(cudarodata){return true;} // already exists
	if(!freeup_devicemem(num_bytes)){return false;}
	cudarodata = get_cuda_array(nx, ny, nz);
	if(!cudarwdata){addtolist();}
	//cout << "ro alloc finish " << " " <<  cudarodata << endl;
	return true;
	
}

void EMData::bindcudaarrayA(const bool intp_mode) const
{
	if(cudarodata == 0){throw UnexpectedBehaviorException( "Cuda Array not allocated!!");}
	if(nz > 1){
		//cout << "3d bind" << endl;
		bind_cuda_array_to_textureA(cudarodata, 3, intp_mode);
	}else{
		//cout << "2d bind" << endl;
		bind_cuda_array_to_textureA(cudarodata, 2, intp_mode);
	}
	
}

void EMData::unbindcudaarryA() const
{
	
	if(nz > 1){
		unbind_cuda_textureA(3);
	}else{
		unbind_cuda_textureA(2);
	}
	
}

void EMData::bindcudaarrayB(const bool intp_mode) const
{
	if(cudarodata == 0){throw UnexpectedBehaviorException( "Cuda Array not allocated!!");}
	if(nz > 1){
		//cout << "3d bind" << endl;
		bind_cuda_array_to_textureB(cudarodata, 3, intp_mode);
	}else{
		//cout << "2d bind" << endl;
		bind_cuda_array_to_textureB(cudarodata, 2, intp_mode);
	}
	
}

void EMData::unbindcudaarryB() const
{
	
	if(nz > 1){
		unbind_cuda_textureB(3);
	}else{
		unbind_cuda_textureB(2);
	}
	
}

bool EMData::copy_from_device(const bool rocpy)
{
	//cout << "copy from device to host " << cudarwdata << " " << rocpy << endl;
	//maybe we should check to see if rdata is still allocated? If not we would need to do either a malloc or new (also assumes that the size of rdata has not changed)
	if(cudarwdata && !rocpy){
		//cout << "rw copy back " << rdata << " numbytes " << num_bytes << endl;
		if(rdata == 0){rdata = (float*)malloc(num_bytes);} //allocate space if needed
		cudaError_t error = cudaMemcpy(rdata,cudarwdata,num_bytes,cudaMemcpyDeviceToHost);
		if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaMemcpy (device to host) failed:" + string(cudaGetErrorString(error)));
		rw_free(); //we have the data on either the host or device, not both (prevents concurrency issues)
	} else if (cudarodata && rocpy) {
		if (nz > 1){
			//cout << "ro copy back 3D" << endl;
			cudaExtent extent;
			extent.width  = nx;
			extent.height = ny;
			extent.depth  = nz;
			cudaMemcpy3DParms copyParams = {0};
			copyParams.srcArray = cudarodata;
			copyParams.dstPtr = make_cudaPitchedPtr((void*)rdata, extent.width*sizeof(float), extent.width, extent.height);
			copyParams.extent   = extent;
			copyParams.kind     = cudaMemcpyDeviceToHost;
			cudaError_t error = cudaMemcpy3D(&copyParams);
			if ( error != cudaSuccess) throw UnexpectedBehaviorException( "RO CudaMemcpy (device to host) failed:" + string(cudaGetErrorString(error)));
		} else{
			//cout << "ro copy back 2D" << endl;
			cudaError_t error = cudaMemcpyFromArray(rdata,cudarodata,0,0,num_bytes,cudaMemcpyDeviceToHost);
			if ( error != cudaSuccess) throw UnexpectedBehaviorException( "RO CudaMemcpy (device to host) failed:" + string(cudaGetErrorString(error)));
		}	
		ro_free(); //we have the data on either the host or device, not both (prevents concurrency issues)
	} else {
		 return false;
	}
	//cout << "finished copying" << endl;
	update(); 
	return true;
}

bool EMData::copy_rw_to_ro() const
{
	
	if(cudarwdata == 0){return false;}
	
	if(cudarodata == 0){
		if(!freeup_devicemem(num_bytes)){return false;}
		cudarodata = get_cuda_array(nx, ny, nz);
		memused += num_bytes;
	}
	//this will copy over any prexisting data (saves a malloc)....
	copy_to_array(cudarwdata, cudarodata, nx, ny, nz, cudaMemcpyDeviceToDevice);
	roneedsupdate = 0; //just copied, so no longer need an update
	
	return true;
	
}

// The policy here is that when an EMData object is created, cudarwdata is set to 0. and no mem is allocated. It is
//only when cudarwdata points to allocated data does the EMData object go on the list. cudarwdata should NEVER be set 
void EMData::runcuda(float * results) const
{
	
	if(results == 0){throw UnexpectedBehaviorException( "Cuda failed!!!");}
	if(cudarwdata != 0){
		//rw_free();} //delete the old data, why not jus overwrite!! (save a cudaFree)
	} else {
		addtolist(); // now that we are using memory add to the list
	}
	cudarwdata = results;
	
}

void EMData::rw_free() const
{
	//cout << "rw_free " << " " << cudarwdata << endl;
	if(usemempoolswitch && mempoolused < (MEMPOOL_SIZE-1) && mempoolarraysize >= int(num_bytes))
	{
		mempoolused++;
		mempool[mempoolused] = cudarwdata;
		return;
	}
	cudaError_t error = cudaFree(cudarwdata);
	if ( error != cudaSuccess){
		cout << rdata << " " << cudarwdata << endl;
		throw UnexpectedBehaviorException( "CudaFree failed:" + string(cudaGetErrorString(error)));
	}
	cudarwdata = 0;
	memused -= num_bytes;
	if(!cudarodata){removefromlist();}
	
}

void EMData::ro_free() const
{
	//cout << "ro_free " << " " << cudarodata << endl;
	cudaError_t error = cudaFreeArray(cudarodata);
	if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaFreeArray failed:" + string(cudaGetErrorString(error)));
	cudarodata = 0;
	memused -= num_bytes;
	if(!cudarwdata){removefromlist();}
	
}

bool EMData::isrodataongpu() const
{
	if(cudarodata != 0 && !roneedsupdate){return true;}
	if(cudarwdata != 0){
		if(copy_rw_to_ro()){;
			return true;
		} else {
			return false;
		}
	}else{
		return false;
	}
	
}
bool EMData::freeup_devicemem(const int& num_bytes) const
{
	size_t freemem, totalmem;
	cudaMemGetInfo(&freemem, &totalmem);
	cout  << "memusage" << " " << freemem << " " << totalmem << endl;
	if ((ptrdiff_t(freemem) - ptrdiff_t(fudgemem)) > ptrdiff_t(num_bytes)){
		return true;
	}else{
		//if(num_bytes > memused){return false;} //it is not possible to free up enough memory!!	
		//keep on removing stuff until enough memory is available
		while(lastinlist != 0){
			if(lastinlist->cudarwdata){
				cudaFree(lastinlist->cudarwdata);
				lastinlist->cudarwdata = 0;
				memused -= lastinlist->nxyz*sizeof(float);
				cudaMemGetInfo(&freemem, &totalmem); //update free memory
			}
			if(lastinlist->cudarodata){
				cudaFreeArray(lastinlist->cudarodata);
				lastinlist->cudarodata = 0;
				memused -= lastinlist->nxyz*sizeof(float);
				cudaMemGetInfo(&freemem, &totalmem); //update free memory
			}
			if(lastinlist != firstinlist){ //if there is more than one itme on the list
				lastinlist->nextlistitem->prevlistitem = 0;	// set the previtem link in the next item to zero
				lastinlist = lastinlist->nextlistitem;		// chop the last item in the list off and set to next item
			}else{
				firstinlist = 0;	// we have deleted everything on the list
				lastinlist = 0;
			}
			if((ptrdiff_t(freemem) - ptrdiff_t(fudgemem)) > ptrdiff_t(num_bytes)){return true;}	//this should break the loop....
		}
	}	
	
	return false;	//if we failed :(
}

void EMData::addtolist() const
{
	if(firstinlist == 0){ //if this is the first item in the list (first object in list)
		firstinlist = this;
		lastinlist = this;
		nextlistitem = 0;
		prevlistitem = 0;
	}else{
		//we add to top of list
		firstinlist->nextlistitem = this;
		prevlistitem = firstinlist;
		nextlistitem = 0;
		firstinlist = this;
	}	
	
}

void EMData::elementaccessed() const
{
        removefromlist();
	//now insert at top (there is no case where we call this function, but there is nothing in the list)
	firstinlist->nextlistitem = this;
	prevlistitem = firstinlist;
	firstinlist = this;
}

void EMData::removefromlist() const
{
	//remove from list
	if(firstinlist == lastinlist){ //last item in list....
		firstinlist = 0;
		lastinlist = 0;
		return;
	}
	if(nextlistitem !=0){
		nextlistitem->prevlistitem = prevlistitem;	//this object is not first in the list
	}else{
		firstinlist = prevlistitem;
	}
	if(prevlistitem !=0){
		prevlistitem->nextlistitem = nextlistitem;	//this item is not last in the list
	}else{
		lastinlist = nextlistitem;
	}
	
}

void EMData::usemempool(int size)
{
	
	usemempoolswitch = true;
	mempoolarraysize = size; // this allow for complex arrays to be stores, but this does waste a little space

}

void EMData::freemempool()
{
	for(int i = 0; i < MEMPOOL_SIZE; i ++)
	{
		if(mempool[i] == 0)
		{
			break;
		}else{
			cudaFree(mempool[i]);
		}
	}
}

void EMData::switchoncuda()
{
	EMData::usecuda = 1;	
}

void EMData::switchoffcuda()
{
	EMData::usecuda = 0;	
}

void EMData::cuda_cleanup() 
{
	do_cuda_fft_cache_destroy();
	//Cleanup any object mess.... CUDA has OCD
	while(lastinlist){
		if(lastinlist->cudarwdata) lastinlist->rw_free();
		if(lastinlist && lastinlist->cudarodata) lastinlist->ro_free();
	}
	
	cudaThreadExit();

}

bool EMData::cuda_initialize()
{
	if(device_init())
	{
		return 1;
	} else {
		switchoffcuda();
		return 0;
	}
}
#endif //EMAN2_USING_CUDA
