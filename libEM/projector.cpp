/**
 * $Id$
 */
#include "projector.h"
#include "log.h"
#include "emdata.h"
#include "interp.h"
#ifdef WIN32
#define M_PI 3.14159265358979323846
#endif

using namespace EMAN;

template<> Factory<Projector>::Factory()
{
    force_add(&GaussFFTProjector::NEW);
    force_add(&PawelProjector::NEW);
    force_add(&SimpleIsoSurfaceProjector::NEW);
    force_add(&StandardProjector::NEW);
    force_add(&StandardBigProjector::NEW);
}



EMData *GaussFFTProjector::project3d(EMData * image) const
{
    EMData *f = image;
    if (image->is_complex()) {
	f = image->do_fft();
    }

    int f_nx = f->get_xsize();
    int f_ny = f->get_ysize();
    int f_nz = f->get_zsize();

    if (!f->is_complex() || f_nz != f_ny || f_nx != f_ny + 2) {
	Log::logger()->error("Cannot project this image");
	return 0;
    }

    f->ap2ri();

    EMData *tmp = new EMData();
    tmp->set_size(f_nx, f_ny, 1);
    tmp->set_complex(true);
    tmp->set_ri(true);

    float *data = tmp->get_data();
    Rotation r(alt, az, phi, Rotation::EMAN);
    Matrix3f mx = r.get_matrix3();
  
    int mode = params["mode"];
    float gauss_width = 0;
    if (mode == 2) {
	gauss_width = 4.0 / (M_PI * M_PI);
    }
    else if (mode == 3) {
	gauss_width = 6.4 / (M_PI * M_PI);
    }
    else if (mode == 4) {
	gauss_width = 8.8 / (M_PI * M_PI);
    }
    else if (mode == 6 || mode == 7) {
	gauss_width = 10.4 / (M_PI * M_PI);
    }
    
    for (int y = 0; y < f_ny; y++) {
	for (int x = 0; x < f_nx / 2; x++) {
	    int ii = x * 2 + y * f_nx;
	    if (hypot(x, y - f_ny / 2) >= f_ny / 2 - 1) {
		data[ii] = 0;
		data[ii + 1] = 0;
	    }
	    else {
		float xx = x * mx[0][0] + (y - f_ny / 2) * mx[0][1];
		float yy = x * mx[0][2] + (y - f_ny / 2) * mx[1][0];
		float zz = x * mx[1][1] + (y - f_ny / 2) * mx[1][2];

		int cc = 1;
		if (xx < 0) {
		    xx = -xx;
		    yy = -yy;
		    zz = -zz;
		    cc = -1;
		}

		yy += f_ny / 2;
		zz += f_nz / 2;

		interp_ft_3d(mode, f, xx, yy, zz, data + ii, gauss_width);
		data[ii + 1] *= cc;
	    }
	}
    }
    f->update();
    tmp->done_data();
    tmp->update();

    tmp->filter("Phase180");
    EMData *ret = tmp->do_ift();

    Dict filter_d;
    filter_d["gauss_width"] = gauss_width;
    filter_d["ring_width"] = ret->get_xsize() / 2;
    ret->filter("InvGaussMaskFilter", filter_d); 

    ret->set_rotation(alt, az, phi);
    
    delete tmp;
    tmp = 0;

    ret->update();
    return ret;
}

void GaussFFTProjector::interp_ft_3d(int mode, EMData * image, float x, float y,
				     float z, float *data, float gw) const
{
    float *rdata = image->get_data();    
    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();
    
    if (mode == 1) {
	int x0 = 2 * (int) floor(x + 0.5);
	int y0 = (int) floor(y + 0.5);
	int z0 = (int) floor(z + 0.5);

	data[0] = rdata[x0 + y0 * nx + z0 * nx * ny];
	data[1] = rdata[x0 + y0 * nx + z0 * nx * ny + 1];
    }
    else if (mode == 2) {
	int x0 = (int) floor(x);
	int y0 = (int) floor(y);
	int z0 = (int) floor(z);

	float dx = x - x0;
	float dy = y - y0;
	float dz = z - z0;

	if (x0 <= nx - 2 && y0 <= ny - 1 && z0 <= nz - 1) {

	    int i = (int) (x0 * 2 + y0 * nx + z0 * nx * ny);

	    float n = Util::agauss(1, dx, dy, dz, gw) +
		Util::agauss(1, 1 - dx, dy, dz, gw) +
		Util::agauss(1, dx, 1 - dy, dz, gw) +
		Util::agauss(1, 1 - dx, 1 - dy, dz, gw) +
		Util::agauss(1, dx, dy, 1 - dz, gw) +
		Util::agauss(1, 1 - dx, dy, 1 - dz, gw) +
		Util::agauss(1, dx, 1 - dy, 1 - dz, gw) + Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, gw);

	    data[0] = Util::agauss(rdata[i], dx, dy, dz, gw) +
		Util::agauss(rdata[i + 2], 1 - dx, dy, dz, gw) +
		Util::agauss(rdata[i + nx], dx, 1 - dy, dz, gw) +
		Util::agauss(rdata[i + nx + 2], 1 - dx, 1 - dy, dz, gw) +
		Util::agauss(rdata[i + nx * ny], dx, dy, 1 - dz, gw) +
		Util::agauss(rdata[i + 2 + nx * ny], 1 - dx, dy, 1 - dz, gw) +
		Util::agauss(rdata[i + nx + nx * ny], dx, 1 - dy, 1 - dz, gw) +
		Util::agauss(rdata[i + 2 + nx + nx * ny], 1 - dx, 1 - dy, 1 - dz, gw) / n;

	    i++;

	    data[1] = Util::agauss(rdata[i], dx, dy, dz, gw) +
		Util::agauss(rdata[i + 2], 1 - dx, dy, dz, gw) +
		Util::agauss(rdata[i + nx], dx, 1 - dy, dz, gw) +
		Util::agauss(rdata[i + nx + 2], 1 - dx, 1 - dy, dz, gw) +
		Util::agauss(rdata[i + nx * ny], dx, dy, 1 - dz, gw) +
		Util::agauss(rdata[i + 2 + nx * ny], 1 - dx, dy, 1 - dz, gw) +
		Util::agauss(rdata[i + nx + nx * ny], dx, 1 - dy, 1 - dz, gw) +
		Util::agauss(rdata[i + 2 + nx + nx * ny], 1 - dx, 1 - dy, 1 - dz, gw) / n;
	}
    }
    else if (mode == 3) {
	int x0 = 2 * (int) floor(x + .5);
	int y0 = (int) floor(y + .5);
	int z0 = (int) floor(z + .5);
	
	float *supp = image->setup4slice();

	if (x0 < nx - 4 && y0 <= ny - 3 && z0 <= nz - 3 && y0 >= 2 && z0 >= 2) {
	    float n = 0;

	    if (x0 == 0) {
		x0 += 4;
		for (int k = z0 - 1; k <= z0 + 1; k++) {
		    for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = x0 - 2; i <= x0 + 2; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x - 2.0, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += supp[i + j * 12 + k * 12 * ny] * g;
			    data[1] += supp[i + j * 12 + k * 12 * ny + 1] * g;
			}
		    }
		}
	    }
	    else {
		for (int k = z0 - 1; k <= z0 + 1; k++) {
		    for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = x0 - 2; i <= x0 + 2; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += rdata[i + j * nx + k * nx * ny] * g;
			    data[1] += rdata[i + j * nx + k * nx * ny + 1] * g;
			}
		    }
		}
	    }

	    data[0] /= n;
	    data[1] /= n;
	}
    }
    else if (mode == 4) {
	int x0 = 2 * (int) floor(x);
	int y0 = (int) floor(y);
	int z0 = (int) floor(z);

	float *supp = image->setup4slice();

	if (x0 < nx - 4 && y0 <= ny - 3 && z0 <= nz - 3 && y0 >= 2 && z0 >= 2) {
	    float n = 0;

	    if (x0 == 0) {
		x0 += 4;
		for (int k = z0 - 1; k <= z0 + 2; k++) {
		    for (int j = y0 - 1; j <= y0 + 2; j++) {
			for (int i = x0 - 2; i <= x0 + 4; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x - 2.0, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += supp[i + j * 12 + k * 12 * ny] * g;
			    data[1] += supp[i + j * 12 + k * 12 * ny + 1] * g;
			}
		    }
		}
	    }
	    else {
		for (int k = z0 - 1; k <= z0 + 2; k++) {
		    for (int j = y0 - 1; j <= y0 + 2; j++) {
			for (int i = x0 - 2; i <= x0 + 4; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += rdata[i + j * nx + k * nx * ny] * g;
			    data[1] += rdata[i + j * nx + k * nx * ny + 1] * g;
			}
		    }
		}
	    }
	    data[0] /= n;
	    data[1] /= n;
	}
    }
    else if (mode == 5) {
	int x0 = 2 * (int) floor(x + .5);
	int y0 = (int) floor(y + .5);
	int z0 = (int) floor(z + .5);

	float *supp = image->setup4slice();
	float *gimx = Interp::get_gimx();

	if (x0 < nx - 4 && y0 <= ny - 3 && z0 <= nz - 3 && y0 >= 2 && z0 >= 2) {
	    int mx0 = -(int) floor((x - x0 / 2) * 39.0 + .5) - 78;
	    int my0 = -(int) floor((y - y0) * 39.0 + .5) - 78;
	    int mz0 = -(int) floor((z - z0) * 39.0 + .5) - 78;

	    float n = 0;
	    int mmz = mz0;
	    int mmy = my0;
	    int mmx = mx0;

	    if (x0 < 4) {
		x0 += 4;

		for (int k = z0 - 2; k <= z0 + 2; k++, mmz += 39) {
		    for (int j = y0 - 2; j <= y0 + 2; j++, mmy += 39) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2, mmx += 39) {
			    float g = gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];
			    n += g;
			    data[0] += supp[i + j * 12 + k * 12 * ny] * g;
			    data[1] += supp[i + j * 12 + k * 12 * ny + 1] * g;
			}
		    }
		}
	    }
	    else {
		for (int k = z0 - 2; k <= z0 + 2; k++, mmz += 39) {
		    for (int j = y0 - 2; j <= y0 + 2; j++, mmy += 39) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2, mmx += 39) {
			    int ii = i + j * nx + k * nx * ny;
			    float g = gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];
			    n += g;
			    data[0] += rdata[ii] * g;
			    data[1] += rdata[ii + 1] * g;
			}
		    }
		}
	    }

	    data[0] /= n;
	    data[1] /= n;
	}
    }
    else if (mode == 6) {

	int x0 = 2 * (int) floor(x + .5);
	int y0 = (int) floor(y + .5);
	int z0 = (int) floor(z + .5);

	float *supp = image->setup4slice();

	if (x0 < nx - 4 && y0 <= ny - 3 && z0 <= nz - 3 && y0 >= 2 && z0 >= 2) {
	    float n = 0;

	    if (x0 < 4) {
		x0 += 4;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
		    for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x - 2.0, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += supp[i + j * 12 + k * 12 * ny] * g;
			    data[1] += supp[i + j * 12 + k * 12 * ny + 1] * g;
			}
		    }
		}
	    }
	    else {
		for (int k = z0 - 2; k <= z0 + 2; k++) {
		    for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2) {
			    float r = Util::hypot3(i / 2.0 - x, j - y, k - z);
			    float g = exp(-r / gw);
			    n += g;
			    data[0] += rdata[i + j * nx + k * nx * ny] * g;
			    data[1] += rdata[i + j * nx + k * nx * ny + 1] * g;
			}

		    }
		}
	    }

	    data[0] /= n;
	    data[1] /= n;
	}
    }
    else if (mode == 7) {
	int x0 = 2 * (int) floor(x + .5);
	int y0 = (int) floor(y + .5);
	int z0 = (int) floor(z + .5);

	float *supp = image->setup4slice();

	if (x0 < nx - 4 && y0 <= ny - 3 && z0 <= nz - 3 && y0 >= 2 && z0 >= 2) {
	    float n = 0;
	    if (x0 < 4) {
		x0 += 4;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
		    for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2) {
			    float r = sqrt(Util::hypot3(i / 2.0 - x - 2.0, j - y, k - z));
			    float g = Interp::hyperg(r);
			    n += g;
			    data[0] += supp[i + j * 12 + k * 12 * ny] * g;
			    data[1] += supp[i + j * 12 + k * 12 * ny + 1] * g;
			}
		    }
		}
	    }
	    else {
		for (int k = z0 - 2; k <= z0 + 2; k++) {
		    for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = x0 - 4; i <= x0 + 4; i += 2) {
			    float r = sqrt(Util::hypot3(i / 2.0 - x, j - y, k - z));
			    float g = Interp::hyperg(r);
			    n += g;
			    data[0] += rdata[i + j * nx + k * nx * ny] * g;
			    data[1] += rdata[i + j * nx + k * nx * ny + 1] * g;
			}

		    }
		}
	    }
	    data[0] /= n;
	    data[1] /= n;
	}
    }
}

#if 1
EMData *PawelProjector::project3d(EMData * image) const
{
    if (!image) {
	return 0;
    }

    int new_origin_x = params["origin_x"];
    int new_origin_y = params["origin_y"];
    int new_origin_z = params["origin_z"];

    float alt = params["alt"];
    float az = params["az"];
    float phi = params["phi"];    

    int new_radius = params["radius"];
    
    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();

    Vec3<int> origins(nx / 2, ny / 2, nz /2);

    int radius = 0;
    if (nz == 1) {
	radius = static_cast<int>(Util::min(nx, ny)) / 2 - 1;
    }
    else {
	radius = static_cast<int>(Util::min(nx, ny, nz)) / 2 - 1;
    }

    if (new_radius > 0) {
	radius = new_radius;
    }
    float radius_square = radius * radius;
    
    int nsec = ny * nz;
    Pointers *xmax = new Pointers[nsec];
    
    
    if (new_origin_x != 0 || new_origin_y != 0 || new_origin_z != 0) {
	origins.set_value(new_origin_x, new_origin_y, new_origin_z);
    }

    int npixels = 0;
    int nrows = 0;
    
    for(int i = 0; i<nz; i++) {
	float dm1 = i - origins[2];
        float sum_1 = dm1 * dm1;
	
        for(int j = 0; j<ny; j++) {
	    float dm2 = j - origins[1];
	    float sum_2 = sum_1 + dm2 *dm2;
	    bool first = true;
	
	    for(int k = 0; k < nx; k++) {
		float dm3 = k - origins[0];
	        float dist = sum_2 + dm3*dm3;
		
	        if(first) {
		    if(dist <= radius_square) {
		        first = false;
			npixels++;
			nrows++;
			int l = nrows-1;
			xmax[l].location.set_value(k, j, i);
			xmax[l].start = npixels;
			xmax[l].end = npixels;
		    }
		}
	        else {
		    if(dist > radius_square) {
		        break;
		    }
		    else {
			npixels++;
			xmax[nrows-1].end = npixels;
		    }
		}
	    }
	}
    }

    Rotation rotation(alt, az, phi, Rotation::EMAN);
    Matrix3f mat = rotation.get_matrix3();
    
    EMData *ret = new EMData();
    ret->set_size(nx, ny, 1);
    float *ret_data = ret->get_data();
    float *rdata = image->get_data();
    int nxy = nx * ny;
    int radius2 = 2*radius+1;
    
    if((2*radius+1) >nx || (2*radius+1)>ny) {    
	for(int i = 0; i<nrows; i++) {
	    Pointers row_data = xmax[i];	    
	    int kb = row_data.location[1] +origins[1];
	    Vec3<float> vb = row_data.location * mat + origins;
	    
	    for(int k = row_data.start-1; k < row_data.end; k++) {
		Vec3<int> iq((int) floor(vb[0]), (int) floor(vb[1]), (int)floor(vb[2]));

		if (radius2 > nx || radius2 > ny) {
		    if(iq[0] >(nx-1) || iq[1] >(ny-1) || iq[2] >(nz-1)) {
			vb += mat.get_vector(0);
			continue;
		    }
		}
		
		Vec3<float> dip = vb - iq;
		int l = iq[0]+iq[1]*nx+iq[2]*nxy;
		float a1 = rdata[l];
		float a2 = rdata[l + 1] - a1;
		float a3 = rdata[l + nx] - a1;
		float a4 = rdata[l + nxy] - a1;
		float a5 = -a2 - rdata[l+nx] + rdata[l+1+nx];
		float a61 = -rdata[l+nxy] + rdata[l+1+nxy];
		float a6 = -a2 + a61;
		float a7 = a3 - rdata[l+nxy] + rdata[l+nx+nxy];
		float a8 = -a5-a61-rdata[l+nx+nxy] + rdata[1+nx+nxy];
		
		float sum_a = a1 +dip[2]*(a4+a6*dip[0]+(a7+a8*dip[0])*dip[1])+a3*dip[1]+dip[0]*(a2+a5*dip[1]);
		ret_data[k +kb*nx] =ret_data[k+kb*nx] +sum_a;

		if (radius2 <= nx && radius2 <= ny) {
		    vb += mat.get_vector(0);
		}
	    }
	}
    }

    ret->done_data();
    ret->update();
    return ret;
}
#endif


EMData *SimpleIsoSurfaceProjector::project3d(EMData * image) const
{
    float alt = params["alt"];
    float az = params["az"];
    float phi = params["phi"];
    float threshold = params["threshold"];

    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();

    EMData *parent = 0;
    if (image->get_parent() == 0) {
	parent = image->copy();
	parent->set_parent(0);
	image->set_parent(parent);
    }

    image->rotate(alt, az, phi);

    EMData *ret = new EMData();
    ret->set_size(nx, ny, 1);
    float *dat = ret->get_data();
    float *din = image->get_data();

    int nxy = nx * ny;
    int k = 0;
    int l = 0;
    for (int i = 0; i < nx; i++) {
	for (int j = 0; j < ny; j++) {
	    int jnx = j * nx;

	    for (k = nz - 2; k > 0; k--) {
		l = i + jnx + k * nxy;
		if (din[l] > threshold)
		    break;
	    }
	    if (k == 0) {
		dat[i + jnx] = -1000;
	    }
	    else {
		dat[i + jnx] = nz - k - (threshold - din[l + nxy]) / (din[l] - din[l + nxy]);
	    }
	}
    }

    float v0[3], v1[3], v2[3];
    for (int j = 0; j < ny; j++) {
	for (int i = 0; i < nx; i++) {
	    k = i + j * nx;
	    float slx = dat[k + 1] + dat[k + nx + 1] - dat[k] - dat[k + nx];
	    float sly = dat[k + nx] + dat[k + nx + 1] - dat[k] - dat[k + 1];

	    if (fabs(slx) > 50 || fabs(sly) > 50 ||
		dat[k + 1] == -1000 || j == ny - 1 || i == nx - 1) {
		dat[k] = 2.0;
	    }
	    else {
		v0[0] = 1.0 / hypot(1, slx);
		v0[1] = 0;
		v0[2] = slx / hypot(1, slx);

		v1[0] = 0;
		v1[1] = 1.0 / hypot(1, sly);
		v1[2] = sly / hypot(1, sly);

		v2[0] = v0[1] * v1[2] - v0[2] * v1[1];
		v2[1] = v0[2] * v1[0] - v0[0] * v1[2];
		v2[2] = v0[0] * v1[1] - v0[1] * v1[0];

		dat[k] = -v2[0] - v2[1] + v2[2];
	    }
	}
    }
    image->done_data();
    ret->done_data();
    ret->update();
    return ret;
}




EMData *StandardProjector::project3d(EMData * image) const
{
    float alt = params["alt"];
    float az = params["az"];
    float phi = params["phi"];

    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();

    EMData *parent = 0;
    if (image->get_parent() == 0) {
	parent = image->copy();
	parent->set_parent(0);
	image->set_parent(parent);
    }

    image->rotate(alt, az, phi);

    EMData *ret = new EMData();
    ret->set_size(nx, ny, 1);
    float *dat = ret->get_data();
    float *din = image->get_data();
    int nxy = nx * ny;

    for (int i = 0; i < nx; i++) {
	for (int j = 0; j < ny; j++) {
	    int jnx = j * nx;
	    dat[i + jnx] = 0;
	    for (int k = 0; k < nz; k++) {
		dat[i + jnx] += din[i + jnx + k * nxy];
	    }
	}
    }

    image->done_data();
    ret->done_data();
    ret->update();
    return ret;
}

EMData *StandardBigProjector::project3d(EMData * image) const
{
    float alt = params["alt"];
    float az = params["az"];
    float phi = params["phi"];

    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();

    Rotation r(alt, az, phi, Rotation::EMAN);
    Matrix3f mx = r.get_matrix3();
    int xy = nx * ny;

    EMData *proj = new EMData();
    proj->set_size(nx, ny, 1);

    float *sdata = image->get_data();
    float *ddata = proj->get_data();
    int l = 0;
    for (int k = -nz / 2; k < nz - nz / 2; k++) {
	for (int j = -ny / 2; j < ny - ny / 2; j++) {
	    for (int i = -nx / 2; i < nx - nx / 2; i++, l++) {
		float x2 = mx[0][0] * i + mx[0][1] * j + mx[0][2] * k + nx / 2;
		float y2 = mx[1][0] * i + mx[1][1] * j + mx[1][2] * k + ny / 2;
		float z2 = mx[2][0] * i + mx[2][1] * j + mx[2][2] * k + nz / 2;

		if (x2 >= 0 && y2 >= 0 && z2 >= 0 && x2 < (nx - 1) && y2 < (ny - 1)
		    && z2 < (nz - 1)) {
		    float x = Util::fast_floor(x2);
		    float y = Util::fast_floor(y2);
		    float z = Util::fast_floor(z2);

		    float t = x2 - x;
		    float u = y2 - y;
		    float v = z2 - z;

		    int ii = (int) (x + y * nx + z * xy);

		    ddata[l] +=
			Util::trilinear_interpolate(sdata[ii], sdata[ii + 1], sdata[ii + nx],
						    sdata[ii + nx + 1], sdata[ii + nx * ny],
						    sdata[ii + xy + 1], sdata[ii + xy + nx],
						    sdata[ii + xy + nx + 1], t, u, v);
		}
	    }
	}
    }

    image->done_data();
    proj->done_data();
    proj->update();
    return proj;
}


void dump_projectors()
{
    dump_factory<Projector>();
}
