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
#include "ctf.h"
#include "emfft.h"
#include "cmp.h"

using namespace EMAN;

// debug only
#include <iostream>
#include <cstring>

using std::cout;
using std::endl;

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_processor.h"
#include "cuda/cuda_cmp.h"
#endif // EMAN2_USING_CUDA

void EMData::free_memory()
{
	ENTERFUNC;
	if (rdata) {
		EMUtil::em_free(rdata);
		rdata = 0;
	}

	if (supp) {
		EMUtil::em_free(supp);
		supp = 0;
	}

	if (rot_fp != 0)
	{
		delete rot_fp;
		rot_fp = 0;
	}
	/*
	nx = 0;
	ny = 0;
	nz = 0;
	nxy = 0;
	 */

	EXITFUNC;
}

void EMData::free_rdata()
{
	ENTERFUNC;
	if (rdata) {
		EMUtil::em_free(rdata);
		rdata = 0;
	}
	EXITFUNC;
}

EMData * EMData::copy() const
{
	ENTERFUNC;

	EMData *ret = new EMData(*this);

	EXITFUNC;
	return ret;
}


EMData *EMData::copy_head() const
{
	ENTERFUNC;
	EMData *ret = new EMData();
	ret->attr_dict = attr_dict;

	ret->set_size(nx, ny, nz);
	ret->flags = flags;

	ret->all_translation = all_translation;

	ret->path = path;
	ret->pathnum = pathnum;

// should these be here? d.woolford I did not comment them out, merely place them here (commented out) to draw attention
// 	ret->xoff = xoff;
// 	ret->yoff = yoff;
// 	ret->zoff = zoff;
// 	ret->changecount = changecount;

	ret->update();

	EXITFUNC;
	return ret;
}

std::complex<float> EMData::get_complex_at(const int &x,const int &y) const{
	if (abs(x)>=nx/2 || abs(y)>ny/2) return std::complex<float>(0,0);
	if (x>=0 && y>=0) return std::complex<float>(rdata[ x*2+y*nx],rdata[x*2+y*nx+1]);
	if (x>0 && y<0) return std::complex<float>(rdata[ x*2+(ny+y)*nx],rdata[x*2+(ny+y)*nx+1]);
	if (x<0 && y>0) return std::complex<float>(  rdata[-x*2+(ny-y)*nx],-rdata[-x*2+(ny-y)*nx+1]);
	return std::complex<float>(rdata[-x*2-y*nx],-rdata[-x*2-y*nx+1]);
}

std::complex<float> EMData::get_complex_at(const int &x,const int &y,const int &z) const{
	if (abs(x)>=nx/2 || abs(y)>ny/2 || abs(z)>nz/2) return std::complex<float>(0,0);

	if (x<0) {
		int idx=-x*2+(y<=0?-y:ny-y)*nx+(z<=0?-z:nz-z)*nxy;
		return std::complex<float>(rdata[idx],-rdata[idx+1]);
	}

	int idx=x*2+(y<0?ny+y:y)*nx+(z<0?nz+z:z)*nxy;
	return std::complex<float>(rdata[idx],rdata[idx+1]);
}

size_t EMData::get_complex_index(const int &x,const int &y,const int &z) const {
	if (abs(x)>=nx/2 || abs(y)>ny/2 || abs(z)>nz/2) return nxyz;
	if (x<0) {
		return -x*2+(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
	}
	return x*2+(y<0?ny+y:y)*(size_t)nx+(z<0?nz+z:z)*(size_t)nxy;
}

size_t EMData::get_complex_index(int x,int y,int z,const int &subx0,const int &suby0,const int &subz0,const int &fullnx,const int &fullny,const int &fullnz) const {
if (abs(x)>=fullnx/2 || abs(y)>fullny/2 || abs(z)>fullnz/2) return nxyz;

if (x<0) {
	x*=-1;
	y*=-1;
	z*=-1;
}
if (y<0) y=fullny+y;
if (z<0) z=fullnz+z;

if (x<subx0||y<suby0||z<subz0||x>=subx0+nx||y>=suby0+ny||z>=subz0+nz) return nxyz;

return (x-subx0)*2+(y-suby0)*(size_t)nx+(z-subz0)*(size_t)nx*(size_t)ny;
}


void EMData::set_complex_at(const int &x,const int &y,const std::complex<float> &val) {
	if (abs(x)>=nx/2 || abs(y)>ny/2) return;
	if (x==0) {
		if (y==0) { rdata[0]=val.real(); rdata[1]=0; }
		else if (y==ny/2 || y==-ny/2) { rdata[ny/2*nx]=val.real(); rdata[ny/2*nx+1]=0; }
		else if (y>0) {
			rdata[y*nx]=rdata[(ny-y)*nx]=val.real();
			rdata[y*nx+1]=val.imag();
			rdata[(ny-y)*nx+1]=-val.imag();
		}
		else {
			rdata[(ny+y)*nx]=rdata[-y*nx]=val.real();
			rdata[(ny+y)*nx+1]=val.imag();
			rdata[-y*nx+1]=-val.imag();
		}
	}
	else if (x==nx/2-1) {
		if (y==0) { rdata[nx-2]=val.real(); rdata[nx-1]=0; }
		else if (y==ny/2 || y==-ny/2) { rdata[ny/2*nx+nx-2]=val.real(); rdata[ny/2*nx+nx-1]=0; }
		else if (y>0) {
			rdata[y*nx+nx-2]=rdata[(ny-y)*nx+nx-2]=val.real();
			rdata[y*nx+nx-1]=val.imag();
			rdata[(ny-y)*nx+nx-1]=-val.imag();
		}
		else {
			rdata[(ny+y)*nx+nx-2]=rdata[-y*nx+nx-2]=val.real();
			rdata[(ny+y)*nx+nx-1]=val.imag();
			rdata[-y*nx+nx-1]=-val.imag();
		}
	}
	else if (x>0 && y>=0) { rdata[ x*2+y*nx]=val.real(); rdata[x*2+y*nx+1]=val.imag(); }
	else if (x>0 && y<0) { rdata[ x*2+(ny+y)*nx]=val.real(); rdata[x*2+(ny+y)*nx+1]=val.imag(); printf("%d %d %d %f\n",x,y,x*2+(ny+y)*nx,val.real()); }
	else if (x<0 && y>0) { rdata[-x*2+(ny-y)*nx]=val.real(); rdata[-x*2+(ny-y)*nx+1]=-val.imag(); }
	else { rdata[-x*2-y*nx]=val.real(); rdata[-x*2+-y*nx+1]=-val.imag(); }
	return;
}

void EMData::set_complex_at(const int &x,const int &y,const int &z,const std::complex<float> &val) {
if (abs(x)>=nx/2 || abs(y)>ny/2 || abs(z)>nz/2) return;

size_t idx;

// for x=0, we need to insert the value in 2 places
// complex conjugate insertion. Removed due to ambiguity with returned index
if (x==0) {
	if (y==0 && z==0) {
		rdata[0]=(float)val.real();
		rdata[1]=0;
		return;
	}
	// complex conjugate in x=0 plane
	size_t idx=(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
	rdata[idx]=(float)val.real();
	rdata[idx+1]=(float)-val.imag();
}
if (abs(x)==nx/2-1) {
	if (y==0 && z==0) {
		rdata[nx-2]=(float)val.real();
		rdata[nx-1]=0;
		return;
	}
	// complex conjugate in x=0 plane
	size_t idx=nx-2+(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
	rdata[idx]=(float)val.real();
	rdata[idx+1]=(float)-val.imag();
}
if (x<0) {
	idx=-x*2+(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
	rdata[idx]=(float)val.real();
	rdata[idx+1]=-(float)val.imag();
	return;
}

idx=x*2+(y<0?ny+y:y)*(size_t)nx+(z<0?nz+z:z)*(size_t)nxy;
rdata[idx]=(float)val.real();
rdata[idx+1]=(float)val.imag();

return;
}

size_t EMData::add_complex_at(const int &x,const int &y,const int &z,const std::complex<float> &val) {
if (abs(x)>=nx/2 || abs(y)>ny/2 || abs(z)>nz/2) return nxyz;

//if (x==0 && abs(y)==16 && abs(z)==1) printf("## %d %d %d\n",x,y,z);
size_t idx;
// for x=0, we need to insert the value in 2 places
if (x==0) {
	if (y==0 && z==0) {
		rdata[0]+=(float)val.real();
		rdata[1]=0;
		return 0;
	}
	// complex conjugate in x=0 plane
	size_t idx=(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
//	if (idx==16*34+1*34*32) printf("a %d %d %d\t%1.5f+%1.5fi\t%1.5f+%1.5fi\n",x,y,z,val.real(),val.imag(),rdata[idx],rdata[idx+1]);
	rdata[idx]+=(float)val.real();
	rdata[idx+1]+=(float)-val.imag();
}
if (abs(x)==nx/2-1) {
	if (y==0 && z==0) {
		rdata[nx-2]+=(float)val.real();
		rdata[nx-1]=0;
		return nx-2;
	}
	// complex conjugate in x=0 plane
	size_t idx=nx-2+(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
//	if (idx==16*34+1*34*32) printf("b %d %d %d\t%1.5f+%1.5fi\t%1.5f+%1.5fi\n",x,y,z,val.real(),val.imag(),rdata[idx],rdata[idx+1]);
	rdata[idx]+=(float)val.real();
	rdata[idx+1]+=(float)-val.imag();
}
if (x<0) {
	idx=-x*2+(y<=0?-y:ny-y)*(size_t)nx+(z<=0?-z:nz-z)*(size_t)nxy;
//	if (idx==16*34+1*34*32) printf("c %d %d %d\t%1.5f+%1.5fi\t%1.5f+%1.5fi\n",x,y,z,val.real(),val.imag(),rdata[idx],rdata[idx+1]);
	rdata[idx]+=(float)val.real();
	rdata[idx+1]+=-(float)val.imag();
	return idx;
}

idx=x*2+(y<0?ny+y:y)*(size_t)nx+(z<0?nz+z:z)*(size_t)nxy;
//if (idx==16*34+1*34*32) printf("d %d %d %d\t%1.5f+%1.5fi\t%1.5f+%1.5fi\n",x,y,z,val.real(),val.imag(),rdata[idx],rdata[idx+1]);
rdata[idx]+=(float)val.real();
rdata[idx+1]+=(float)val.imag();

return idx;
}

size_t EMData::add_complex_at(int x,int y,int z,const int &subx0,const int &suby0,const int &subz0,const int &fullnx,const int &fullny,const int &fullnz,const std::complex<float> &val) {
if (abs(x)>=fullnx/2 || abs(y)>fullny/2 || abs(z)>fullnz/2) return nxyz;
//if (x==0 && (y!=0 || z!=0)) add_complex_at(0,-y,-z,subx0,suby0,subz0,fullnx,fullny,fullnz,conj(val));
// complex conjugate insertion. Removed due to ambiguity with returned index
/*if (x==0&& (y!=0 || z!=0)) {
	int yy=y<=0?-y:fullny-y;
	int zz=z<=0?-z:fullnz-z;

	if (yy<suby0||zz<subz0||yy>=suby0+ny||zz>=subz0+nz) return nx*ny*nz;

	size_t idx=(yy-suby0)*nx+(zz-subz0)*nx*ny;
	rdata[idx]+=(float)val.real();
	rdata[idx+1]+=(float)-val.imag();
}*/
float cc=1.0;
if (x<0) {
	x*=-1;
	y*=-1;
	z*=-1;
	cc=-1.0;
}
if (y<0) y=fullny+y;
if (z<0) z=fullnz+z;

if (x<subx0||y<suby0||z<subz0||x>=subx0+nx||y>=suby0+ny||z>=subz0+nz) return nxyz;

size_t idx=(x-subx0)*2+(y-suby0)*(size_t)nx+(z-subz0)*(size_t)nx*ny;
rdata[idx]+=(float)val.real();
rdata[idx+1]+=cc*(float)val.imag();
return idx;
}


void EMData::add(float f,int keepzero)
{
	ENTERFUNC;

	float * data = get_data();
	if( is_real() )
	{
		if (f != 0) {

			size_t size = nxyz;
			if (keepzero) {
				for (size_t i = 0; i < size; i++) {
					if (data[i]) data[i] += f;
				}
			}
			else {
				for (size_t i = 0; i < size; i++) {
					data[i] += f;
				}
			}
			update();
		}
	}
	else if( is_complex() )
	{
		if( f!=0 )
		{
			update();
			size_t size = (size_t)nx*ny*nz; //size of data
			if( keepzero )
			{
				for(size_t i=0; i<size; i+=2)
				{
					if (data[i]) data[i] += f;
				}
			}
			else
			{
				for(size_t i=0; i<size; i+=2)
				{
					data[i] += f;
				}
			}
		}
	}
	else
	{
		throw ImageFormatException("This image is neither a real nor a complex image.");
	}
	update();
	EXITFUNC;
}


//for add operation, real and complex image is the same
void EMData::add(const EMData & image)
{
	ENTERFUNC;
	if (nx != image.get_xsize() || ny != image.get_ysize() || nz != image.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( (is_real()^image.is_real()) == true )
	{
		throw ImageFormatException( "not support add between real image and complex image");
	}
	else {

		const float *src_data = image.get_data();
		size_t size = nxyz;
		float* data = get_data();

		for (size_t i = 0; i < size; i++) {
			data[i] += src_data[i];
		}
		update();
	}
	EXITFUNC;
}

//for add operation, real and complex image is the same
void EMData::addsquare(const EMData & image)
{
	ENTERFUNC;
	if (nx != image.get_xsize() || ny != image.get_ysize() || nz != image.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( this->is_complex() || image.is_complex() )
	{
		throw ImageFormatException( "Cannot addsquare() with complex images");
	}
	else {

		const float *src_data = image.get_data();
		size_t size = nxyz;
		float* data = get_data();

		for (size_t i = 0; i < size; i++) {
			data[i] += src_data[i]*src_data[i];
		}
		update();
	}
	EXITFUNC;
}

//for add operation, real and complex image is the same
void EMData::subsquare(const EMData & image)
{
	ENTERFUNC;
	if (nx != image.get_xsize() || ny != image.get_ysize() || nz != image.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( this->is_complex() || image.is_complex() )
	{
		throw ImageFormatException( "Cannot addsquare() with complex images");
	}
	else {

		const float *src_data = image.get_data();
		size_t size = nxyz;
		float* data = get_data();

		for (size_t i = 0; i < size; i++) {
			data[i] -= src_data[i]*src_data[i];
		}
		update();
	}
	EXITFUNC;
}


void EMData::sub(float f)
{
	ENTERFUNC;

#ifdef EMAN2_USING_CUDA
		if (EMData::usecuda == 1 && cudarwdata) {
			if(f != 0){
				subtract_cuda(cudarwdata, f, nx, ny, nz);
			}
			EXITFUNC;
			return;
		}
#endif // EMAN2_USING_CUDA

	float* data = get_data();
	if( is_real() )
	{
		if (f != 0) {
			size_t size = nxyz;
			for (size_t i = 0; i < size; i++) {
				data[i] -= f;
			}
		}
		update();
	}
	else if( is_complex() )
	{
		if( f != 0 )
		{
			size_t size = nxyz;
			for( size_t i=0; i<size; i+=2 )
			{
				data[i] -= f;
			}
		}
		update();
	}
	else
	{
		throw ImageFormatException("This image is neither a real nor a complex image.");
	}

	EXITFUNC;
}


//for sub operation, real and complex image is the same
void EMData::sub(const EMData & em)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException("images not same sizes");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "not support sub between real image and complex image");
	}
	else {
		const float *src_data = em.get_data();
		size_t size = nxyz;
		float* data = get_data();

		for (size_t i = 0; i < size; i++) {
			data[i] -= src_data[i];
		}
		update();
	}
	EXITFUNC;
}


void EMData::mult(float f)
{
	ENTERFUNC;

// this will cause a crash if CUDA is used(no rdata) and a complex map is given.....
	if (is_complex()) {
		ap2ri();
	}
	if (f != 1.0) {
#ifdef EMAN2_USING_CUDA
		if (EMData::usecuda == 1 && cudarwdata) { //doesn't make any sense to use RO, esp on compute devices >= 2.0
			//cout << "CUDA mult" << endl;
			emdata_processor_mult(cudarwdata,f,nx,ny,nz);
			EXITFUNC;
			return;
		}
#endif // EMAN2_USING_CUDA
		float* data = get_data();
		size_t size = nxyz;
		for (size_t i = 0; i < size; i++) {
			data[i] *= f;
		}
		update();
	}
	EXITFUNC;
}


void EMData::mult(const EMData & em, bool prevent_complex_multiplication)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException( "can not multiply images that are not the same size");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "can not multiply real and complex images.");
	}
	else
	{
		const float *src_data = em.get_data();
		size_t size = nxyz;
		float* data = get_data();
		if( is_real() || prevent_complex_multiplication )
		{
			for (size_t i = 0; i < size; i++) {
				data[i] *= src_data[i];
			}
		}
		else
		{
			typedef std::complex<float> comp;
			for( size_t i = 0; i < size; i+=2 )
			{
				comp c_src( src_data[i], src_data[i+1] );
				comp c_rdat( data[i], data[i+1] );
				comp c_result = c_src * c_rdat;
				data[i] = c_result.real();
				data[i+1] = c_result.imag();
			}
		}
		update();
	}

	EXITFUNC;
}

void EMData::mult_complex_efficient(const EMData & em, const int radius)
{
	ENTERFUNC;

	if( is_real() || em.is_real() )throw ImageFormatException( "can call mult_complex_efficient unless both images are complex");


	const float *src_data = em.get_data();

	size_t i_radius = radius;
	size_t k_radius = 1;
	size_t j_radius = 1;
	int ndim = get_ndim();

	if (ndim != em.get_ndim()) throw ImageDimensionException("Can't do that");

	if ( ndim == 3 ) {
		k_radius = radius;
		j_radius = radius;
	} else if ( ndim == 2 ) {
		j_radius = radius;
	}


	int s_nx = em.get_xsize();
	int s_nxy = s_nx*em.get_ysize();

	size_t r_size = nxyz;
	int s_size = s_nxy*em.get_zsize();
	float* data = get_data();

	for (size_t k = 0; k < k_radius; ++k ) {
		for (size_t j = 0; j < j_radius; j++) {
			for (size_t i = 0; i < i_radius; i++) {
				int r_idx = k*nxy + j*nx + i;
				int s_idx = k*s_nxy + j*s_nx + i;
				data[r_idx] *= src_data[s_idx];
				data[r_size-r_idx-1] *= src_data[s_size-s_idx-1];
			}
		}
	}

	update();

	EXITFUNC;
}


void EMData::div(float f)
{
	ENTERFUNC;
 	if ( f == 0 ) {
 		throw InvalidValueException(f,"Can not divide by zero");
 	}
	mult(1.0f/f);
	EXITFUNC;
}


void EMData::div(const EMData & em)
{
	ENTERFUNC;

	if (nx != em.get_xsize() || ny != em.get_ysize() || nz != em.get_zsize()) {
		throw ImageFormatException( "images not same sizes");
	}
	else if( (is_real()^em.is_real()) == true )
	{
		throw ImageFormatException( "not support division between real image and complex image");
	}
	else {
		const float *src_data = em.get_data();
		size_t size = nxyz;
		float* data = get_data();

		if( is_real() )
		{
			for (size_t i = 0; i < size; i++) {
				if(src_data[i] != 0) {
					data[i] /= src_data[i];
				}
				else {
					if (data[i]==0) continue;
					throw InvalidValueException(src_data[i], "divide by zero");
				}
			}
		}
		else
		{
			typedef std::complex<float> comp;
			for( size_t i = 0; i < size; i+=2 )
			{
				comp c_src( src_data[i], src_data[i+1] );
				comp c_rdat( data[i], data[i+1] );
				comp c_result = c_rdat / c_src;
				data[i] = c_result.real();
				data[i+1] = c_result.imag();
			}
		}
		update();
	}

	EXITFUNC;
}


// just a shortcut for cmp("dot")
float EMData::dot(EMData * with)
{
	ENTERFUNC;
	if (!with) {
		throw NullPointerException("Null EMData Image");
	}
	DotCmp dot_cmp;
	float r = -dot_cmp.cmp(this, with);
	EXITFUNC;
	return r;
}


EMData *EMData::get_row(int row_index) const
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D image only");
	}

	EMData *ret = new EMData();
	ret->set_size(nx, 1, 1);
	memcpy(ret->get_data(), get_data() + nx * row_index, nx * sizeof(float));
	ret->update();
	EXITFUNC;
	return ret;
}


void EMData::set_row(const EMData * d, int row_index)
{
	ENTERFUNC;

	if (get_ndim() > 2) {
		throw ImageDimensionException("1D/2D image only");
	}
	if (d->get_ndim() != 1) {
		throw ImageDimensionException("1D image only");
	}

	float *dst = get_data();
	float *src = d->get_data();
	memcpy(dst + nx * row_index, src, nx * sizeof(float));
	update();
	EXITFUNC;
}

EMData *EMData::get_col(int col_index) const
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}

	EMData *ret = new EMData();
	ret->set_size(ny, 1, 1);
	float *dst = ret->get_data();
	float *src = get_data();

	for (int i = 0; i < ny; i++) {
		dst[i] = src[i * nx + col_index];
	}

	ret->update();
	EXITFUNC;
	return ret;
}


void EMData::set_col(const EMData * d, int n)
{
	ENTERFUNC;

	if (get_ndim() != 2) {
		throw ImageDimensionException("2D image only");
	}
	if (d->get_ndim() != 1) {
		throw ImageDimensionException("1D image only");
	}

	float *dst = get_data();
	float *src = d->get_data();

	for (int i = 0; i < ny; i++) {
		dst[i * nx + n] = src[i];
	}

	update();
	EXITFUNC;
}

float& EMData::get_value_at_wrap(int x)
{
	if (x < 0) x = nx + x;
	return get_data()[x];
}

float& EMData::get_value_at_wrap(int x, int y)
{
	if (x < 0) x = nx + x;
	if (y < 0) y = ny + y;

	return get_data()[x + y * nx];
}

float& EMData::get_value_at_wrap(int x, int y, int z)
{

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1 && cudarwdata){
		float result = get_value_at_wrap_cuda(cudarwdata, x, y, z, nx, ny, nz); // this should work....
		return result;
	}
#endif
	int lx = x;
	int ly = y;
	int lz = z;

	if (lx < 0) lx = nx + lx;
	if (ly < 0) ly = ny + ly;
	if (lz < 0) lz = nz + lz;

	return get_data()[lx + ly * (size_t)nx + lz * (size_t)nxy];
}

float EMData::get_value_at_wrap(int x) const
{
	if (x < 0) x = nx - x;
	return get_data()[x];
}

float EMData::get_value_at_wrap(int x, int y) const
{
	if (x < 0) x = nx - x;
	if (y < 0) y = ny - y;

	return get_data()[x + y * nx];
}

float EMData::get_value_at_wrap(int x, int y, int z) const
{
	ptrdiff_t lx = x;
	ptrdiff_t ly = y;
	ptrdiff_t lz = z;
	if (lx < 0) lx = nx + lx;
	if (ly < 0) ly = ny + ly;
	if (lz < 0) lz = nz + lz;

	return get_data()[lx + ly * nx + lz * nxy];
}

float EMData::sget_value_at(int x, int y, int z) const
{
	if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
		return 0;
	}
	return get_data()[(size_t)x + (size_t)y * (size_t)nx + (size_t)z * (size_t)nxy];
}


float EMData::sget_value_at(int x, int y) const
{
	if (x < 0 || x >= nx || y < 0 || y >= ny) {
		return 0;
	}
	return get_data()[x + y * nx];
}


float EMData::sget_value_at(size_t i) const
{
	size_t size = nx*ny;
	size *= nz;
	if (i >= size) {
		return 0;
	}
	return get_data()[i];
}


float EMData::sget_value_at_interp(float xx, float yy) const
{
	int x = static_cast < int >(Util::fast_floor(xx));
	int y = static_cast < int >(Util::fast_floor(yy));

	float p1 = sget_value_at(x, y);
	float p2 = sget_value_at(x + 1, y);
	float p3 = sget_value_at(x, y + 1);
	float p4 = sget_value_at(x + 1, y + 1);

	float result = Util::bilinear_interpolate(p1, p2, p3, p4, xx - x, yy - y);
	return result;
}


float EMData::sget_value_at_interp(float xx, float yy, float zz) const
{
	int x = (int) Util::fast_floor(xx);
	int y = (int) Util::fast_floor(yy);
	int z = (int) Util::fast_floor(zz);

	float p1 = sget_value_at(x, y, z);
	float p2 = sget_value_at(x + 1, y, z);
	float p3 = sget_value_at(x, y + 1, z);
	float p4 = sget_value_at(x + 1, y + 1, z);

	float p5 = sget_value_at(x, y, z + 1);
	float p6 = sget_value_at(x + 1, y, z + 1);
	float p7 = sget_value_at(x, y + 1, z + 1);
	float p8 = sget_value_at(x + 1, y + 1, z + 1);

	float result = Util::trilinear_interpolate(p1, p2, p3, p4, p5, p6, p7, p8,
											   xx - x, yy - y, zz - z);

	return result;
}


EMData & EMData::operator+=(float n)
{
	add(n);
	update();
	return *this;
}


EMData & EMData::operator-=(float n)
{
	*this += (-n);
	return *this;
}


EMData & EMData::operator*=(float n)
{
	mult(n);
	update();
	return *this;
}


EMData & EMData::operator/=(float n)
{
	if (n == 0) {
		LOGERR("divided by zero");
		return *this;
	}
	*this *= (1.0f / n);
	return *this;
}


EMData & EMData::operator+=(const EMData & em)
{
	add(em);
	update();
	return *this;
}


EMData & EMData::operator-=(const EMData & em)
{
	sub(em);
	update();
	return *this;
}


EMData & EMData::operator*=(const EMData & em)
{
	mult(em);
	update();
	return *this;
}


EMData & EMData::operator/=(const EMData & em)
{
	div(em);
	update();
	return *this;
}


EMData * EMData::power(int n) const
{
	ENTERFUNC;

	if( n<0 ) {
		throw InvalidValueException(n, "the power of negative integer not supported.");
	}

	EMData * r = this->copy();
	if( n == 0 ) {
		r->to_one();
	}
	else if( n>1 ) {
		for( int i=1; i<n; i++ ) {
			*r *= *this;
		}
	}

	r->update();
	return r;

	EXITFUNC;
}


EMData * EMData::sqrt() const
{
	ENTERFUNC;

	if (is_complex()) {
		throw ImageFormatException("real image only");
	}

	EMData * r = this->copy();
	float * new_data = r->get_data();
	float * data = get_data();
	size_t size = nxyz;
	for (size_t i = 0; i < size; ++i) {
		if(data[i] < 0) {
			throw InvalidValueException(data[i], "pixel value must be non-negative for logrithm");
		}
		else {
			if(data[i]) {	//do nothing with pixel has value zero
				new_data[i] = std::sqrt(data[i]);
			}
		}
	}

	r->update();
	return r;

	EXITFUNC;
}


EMData * EMData::log() const
{
	ENTERFUNC;

	if (is_complex()) {
		throw ImageFormatException("real image only");
	}

	EMData * r = this->copy();
	float * new_data = r->get_data();
	float * data = get_data();
	size_t size = nxyz;
	for (size_t i = 0; i < size; ++i) {
		if(data[i] < 0) {
			throw InvalidValueException(data[i], "pixel value must be non-negative for logrithm");
		}
		else {
			if(data[i]) {	//do nothing with pixel has value zero
				new_data[i] = std::log(data[i]);
			}
		}
	}

	r->update();
	return r;

	EXITFUNC;
}


EMData * EMData::log10() const
{
	ENTERFUNC;

	if (is_complex()) {
		throw ImageFormatException("real image only");
	}

	EMData * r = this->copy();
	float * new_data = r->get_data();
	float * data = get_data();
	size_t size = nxyz;
	for (size_t i = 0; i < size; ++i) {
		if(data[i] < 0) {
			throw InvalidValueException(data[i], "pixel value must be non-negative for logrithm");
		}
		else {
			if(data[i]) {	//do nothing with pixel has value zero
				new_data[i] = std::log10(data[i]);
			}
		}
	}

	r->update();
	return r;

	EXITFUNC;
}


EMData * EMData::real() const //real part has half of x dimension for a complex image
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) // a real image, return a copy of itself
	{
		e = this->copy();
	}
	else //for a complex image
	{
		if( !is_ri() )
		{
			delete e;
			throw InvalidCallException("This image is in amplitude/phase format, this function call require a complex image in real/imaginary format.");
		}
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		float * data = get_data();
		size_t idx1, idx2;
		for( int i=0; i<nx; ++i )
		{
			for( int j=0; j<ny; ++j )
			{
				for( int k=0; k<nz; ++k )
				{
					if( i%2 == 0 )
					{
						//complex data in format [real, complex, real, complex...]
						idx1 = i/2+j*(nx/2)+k*(nx/2)*ny;
						idx2 = i+j*nx+k*nx*ny;
						edata[idx1] = data[idx2];
					}
				}
			}
		}
	}

	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update();
	return e;

	EXITFUNC;
}


EMData * EMData::imag() const
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) {	//a real image has no imaginary part, throw exception
		throw InvalidCallException("No imaginary part for a real image, this function call require a complex image.");
	}
	else {	//for complex image
		if( !is_ri() ) {
			throw InvalidCallException("This image is in amplitude/phase format, this function call require a complex image in real/imaginary format.");
		}
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		float * data = get_data();
		for( int i=0; i<nx; i++ ) {
			for( int j=0; j<ny; j++ ) {
				for( int k=0; k<nz; k++ ) {
					if( i%2 == 1 ) {
						//complex data in format [real, complex, real, complex...]
						edata[i/2+j*(nx/2)+k*(nx/2)*ny] = data[i+j*nx+k*nx*ny];
					}
				}
			}
		}
	}

	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update();
	return e;

	EXITFUNC;
}

EMData * EMData::absi() const//abs has half of x dimension for a complex image
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) // a real image
	{
		e = this->copy();
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		float *edata = e->get_data();
		float * data = get_data();
		size_t idx;
		for( int i=0; i<nx; ++i ) {
			for( int j=0; j<ny; ++j ) {
				for( int k=0; k<nz; ++k ) {
					idx = i+j*nx+k*nx*ny;
					edata[idx] = std::abs(data[idx]);
				}
			}
		}
	}
	else //for a complex image
	{
		if( !is_ri() )
		{
			throw InvalidCallException("This image is in amplitude/phase format, this function call require a complex image in real/imaginary format.");
		}
		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		float * data = get_data();
		size_t idx1, idx2;
		for( int i=0; i<nx; ++i )
		{
			for( int j=0; j<ny; ++j )
			{
				for( int k=0; k<nz; ++k )
				{
					if( i%2 == 0 )
					{
						idx1 = i/2+j*(nx/2)+k*(nx/2)*ny;
						idx2 = i+j*nx+k*nx*ny;
						//complex data in format [real, complex, real, complex...]
						edata[idx1] =
						std::sqrt(data[idx2]*data[idx2]+data[idx2+1]*data[idx2+1]);
					}
				}
			}
		}
	}

	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update();
	return e;

	EXITFUNC;
}


EMData * EMData::amplitude() const
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) {
		throw InvalidCallException("No imaginary part for a real image, this function call require a complex image.");
	}
	else {
		if(is_ri()) {
			throw InvalidCallException("This image is in real/imaginary format, this function call require a complex image in amplitude/phase format.");
		}

		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		float * data = get_data();
		size_t idx1, idx2;
		for( int i=0; i<nx; ++i )
		{
			for( int j=0; j<ny; ++j )
			{
				for( int k=0; k<nz; ++k )
				{
					if( i%2 == 0 )
					{
						idx1 = i/2+j*(nx/2)+k*(nx/2)*ny;
						idx2 = i+j*nx+k*nx*ny;
						//complex data in format [amp, phase, amp, phase...]
						edata[idx1] = data[idx2];
					}
				}
			}
		}
	}

	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update();
	return e;

	EXITFUNC;
}

EMData * EMData::phase() const
{
	ENTERFUNC;

	EMData * e = new EMData();

	if( is_real() ) {
		delete e;
		throw InvalidCallException("No imaginary part for a real image, this function call require a complex image.");
	}
	else {
		if(is_ri()) {
			delete e;
			throw InvalidCallException("This image is in real/imaginary format, this function call require a complex image in amplitude/phase format.");
		}

		int nx = get_xsize();
		int ny = get_ysize();
		int nz = get_zsize();
		e->set_size(nx/2, ny, nz);
		float * edata = e->get_data();
		float * data = get_data();
		size_t idx1, idx2;
		for( int i=0; i<nx; ++i ) {
			for( int j=0; j<ny; ++j ) {
				for( int k=0; k<nz; ++k ) {
					if( i%2 == 1 ) {
						idx1 = i/2+j*(nx/2)+k*(nx/2)*ny;
						idx2 = i+j*nx+k*nx*ny;
						//complex data in format [real, complex, real, complex...]
						edata[idx1] = data[idx2];
					}
				}
			}
		}
	}

	e->set_complex(false);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(false);
	}
	e->update();
	return e;

	EXITFUNC;
}

EMData * EMData::real2complex(const float img) const
{
	ENTERFUNC;

	if( is_complex() ) {
		throw InvalidCallException("This function call only apply to real image");
	}

	EMData * e = new EMData();
	int nx = get_xsize();
	int ny = get_ysize();
	int nz = get_zsize();
	e->set_size(nx*2, ny, nz);

	for( int k=0; k<nz; ++k ) {
		for( int j=0; j<ny; ++j ) {
			for( int i=0; i<nx; ++i ) {
				(*e)(i*2,j,k) = (*this)(i,j,k);
				(*e)(i*2+1,j,k) = img;
			}
		}
	}

	e->set_complex(true);
	if(e->get_ysize()==1 && e->get_zsize()==1) {
		e->set_complex_x(true);
	}
	e->set_ri(true);
	e->update();
	return e;

	EXITFUNC;
}

void EMData::to_zero()
{
	ENTERFUNC;

	if (is_complex()) {
		set_ri(true);
	}
	else {
		set_ri(false);
	}

	//EMUtil::em_memset(get_data(), 0, nxyz * sizeof(float));
	to_value(0.0);
	update();
	EXITFUNC;
}

void EMData::to_one()
{
	ENTERFUNC;

	if (is_complex()) {
		set_ri(true);
	}
	else {
		set_ri(false);
	}
	to_value(1.0);

	update();
	EXITFUNC;
}

void EMData::to_value(const float& value)
{
	ENTERFUNC;

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1 && cudarwdata){
		to_value_cuda(cudarwdata, value, nx, ny, nz);
		return;
	}
#endif // EMAN2_USING_CUDA
	float* data = get_data();

	//the em_memset has segfault for >8GB image, use std::fill() instead, though may be slower
//	if ( value != 0 ) std::fill(data,data+get_size(),value);
//	else EMUtil::em_memset(data, 0, nxyz * sizeof(float)); // This might be faster, I don't know

	std::fill(data,data+get_size(),value);

	update();
	EXITFUNC;
}


