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

#define MEMPOOL_SIZE 10
public:
	
	/** Copy data from host to GPU device (RW)
	*/
	bool copy_to_cuda();
	
	/** Copy data from host to GPU device (RO)
	*/
	bool copy_to_cudaro();
	
	/** Copy data from GPU device to host
	*/
	bool copy_from_device(const bool rocpy = false);
	
	/** Copy from rw memory to ro memory
	*/
	bool copy_rw_to_ro() const; //we are not changing any of the CPU data!!!
	
	/** move this EMData object to the topof accessed list
	*/
	void elementaccessed();
	
	/** bind cudaarray to a texture A
	*/
	void bindcudaarrayA(const bool intp_mode) const;
	
	/** unbind a cuda array from texture A
	*/
	void unbindcudaarryA() const;
	
	/** bind cudaarray to a texture B
	*/
	void bindcudaarrayB(const bool intp_mode) const;
	
	/** unbind a cuda array from texture B
	*/
	void unbindcudaarryB() const;
	
	/** run a cuda function
	*/
	void runcuda(float * results);
	
	/** check to see is there is any data on gpu ro
	* if there is any on rw then copy from rw to ro
	*/
	bool isrodataongpu() const;
	
inline	void roneedsanupdate()
{
	
	roneedsupdate = true;
}

	/** switch to turn on mempool usage */
	static void usemempool(int size);
	
	/** freemempool */
	static void freemempool();
	
//private:
	/** allocate RW on GPU device
	*/
	bool rw_alloc();
	
	/** allocate RO on GPU device 
	*/
	bool ro_alloc();
	
	/** deallocate RW on GPU device
	*/
	void rw_free();
	
	/** deallocate RO on GPU device
	*/
	void ro_free();
	
	/** add this EMData object to a list
	*/
	void addtolist();
	
	/** remove this EMData object from a list
	*/
	void removefromlist();
	
	/** free up device memory by removing last acessed items;
	*/
	bool freeup_devicemem(const int& num_bytes) const;
	
	static void switchoncuda();
	
	static void switchoffcuda();
	
	//pointers to cuda data
	mutable float* cudarwdata;	//we can still change GPU data on a cost object
	mutable cudaArray* cudarodata;	//we can still change GPU data on a cost object
	
	size_t num_bytes;
	
	//pointers used in doubly limked list
	EMData* nextlistitem;
	EMData* prevlistitem;
	
	mutable bool roneedsupdate;
	
	static int memused;
	static int fudgemem;
	static EMData* firstinlist;
	static EMData* lastinlist;
	static bool usecuda;
	
	// mempool stuff
	static float* mempool[MEMPOOL_SIZE];
	static int  mempoolused;
	static int mempoolarraysize;
	static bool usemempoolswitch;
	
#endif // EMAN2_USING_CUDA	

#endif //eman__emdatacuda_h__ 1

