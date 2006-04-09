/**
 * $Id$
 */
#include "projector.h"
#include "emdata.h"
#include "interp.h"

#ifdef WIN32
	#define M_PI 3.14159265358979323846f
#endif	//WIN32

using namespace std;
using namespace EMAN;

template <> Factory < Projector >::Factory()
{
	force_add(&GaussFFTProjector::NEW);
	force_add(&PawelProjector::NEW);
	force_add(&SimpleIsoSurfaceProjector::NEW);
	force_add(&StandardProjector::NEW);
	force_add(&StandardFastProjector::NEW);
	force_add(&FourierGriddingProjector::NEW);
}



EMData *GaussFFTProjector::project3d(EMData * image) const
{
	EMData *f = image;
	if (!image->is_complex()) {
		f->process_inplace("eman1.xform.phaseorigin");
		f = image->do_fft();
		image->process_inplace("eman1.xform.phaseorigin");
		f->process_inplace("eman1.xform.fourierorigin");
	}

	int f_nx = f->get_xsize();
	int f_ny = f->get_ysize();
	int f_nz = f->get_zsize();
	float alt2=params["alt"], az2=params["az"], phi2=params["phi"]; //fix this

	if (!f->is_complex() || f_nz != f_ny || f_nx != f_ny + 2) {
		LOGERR("Cannot project this image");
		return 0;
	}

	f->ap2ri();

	EMData *tmp = new EMData();
	tmp->set_size(f_nx, f_ny, 1);
	tmp->set_complex(true);
	tmp->set_ri(true);

	float *data = tmp->get_data();
	Transform3D r( az2*static_cast<float>(dgr_to_rad), 
					alt2*static_cast<float>(dgr_to_rad), 
					phi2*static_cast<float>(dgr_to_rad) ); // EMAN by default
	r.transpose();

	int mode = params["mode"];
	float gauss_width = 0;
	if (mode == 2) {
		gauss_width = EMConsts::I2G;
	}
	else if (mode == 3) {
		gauss_width = EMConsts::I3G;
	}
	else if (mode == 4) {
		gauss_width = EMConsts::I4G;
	}
	else if (mode == 6 || mode == 7) {
		gauss_width = EMConsts::I5G; 
	}

	for (int y = 0; y < f_ny; y++) {
		for (int x = 0; x < f_nx / 2; x++) {
			int ii = x * 2 + y * f_nx;
			if (hypot(x, y - f_ny / 2) >= f_ny / 2 - 1) {
				data[ii] = 0;
				data[ii + 1] = 0;
			}
			else {
				float xx = (float) (x * r[0][0] + (y - f_ny / 2) * r[0][1]);
				float yy = (float) (x * r[1][0] + (y - f_ny / 2) * r[1][1]);
				float zz = (float) (x * r[2][0] + (y - f_ny / 2) * r[2][1]);

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

	tmp->process_inplace("eman1.xform.fourierorigin");
	EMData *ret = tmp->do_ift();
	ret->process_inplace("eman1.xform.phaseorigin");

	Dict filter_d;
	filter_d["gauss_width"] = gauss_width;
	filter_d["ring_width"] = ret->get_xsize() / 2;
	ret->process_inplace("eman1.math.gausskernelfix", filter_d);

	ret->set_rotation( az2*static_cast<float>(dgr_to_rad), 
						alt2*static_cast<float>(dgr_to_rad), 
						phi2*static_cast<float>(dgr_to_rad) );

	if( tmp )
	{ 
		delete tmp;
		tmp = 0;
	}

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
		int x0 = 2 * (int) floor(x + 0.5f);
		int y0 = (int) floor(y + 0.5f);
		int z0 = (int) floor(z + 0.5f);

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
				Util::agauss(1, dx, 1 - dy, 1 - dz, gw) + Util::agauss(1, 1 - dx, 1 - dy, 1 - dz,
																	   gw);

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
							float r = Util::hypot3(i / 2.0f - x - 2.0f, j - y, k - z);
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
							float r = Util::hypot3(i / 2.0f - x, j - y, k - z);
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
							float r = Util::hypot3(i / 2.0f - x - 2.0f, j - y, k - z);
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
							float r = Util::hypot3(i / 2.0f - x, j - y, k - z);
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
							float r = Util::hypot3(i / 2.0f - x - 2.0f, j - y, k - z);
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
							float r = Util::hypot3(i / 2.0f - x, j - y, k - z);
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
							float r = sqrt(Util::hypot3(i / 2.0f - x - 2.0f, j - y, k - z));
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
							float r = sqrt(Util::hypot3(i / 2.0f - x, j - y, k - z));
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

void PawelProjector::prepcubes(int nx, int ny, int nz, int ri, Vec3i origin, 
		                       int& nn, IPCube* ipcube) const {
	const float r = float(ri*ri);
	const int ldpx = origin[0];
	const int ldpy = origin[1];
	const int ldpz = origin[2];
	float t;
	nn = -1;
	for (int i1 = 0; i1 < nz; i1++) {
		t = float(i1 - ldpz); 
		const float xx = t*t;
		for (int i2 = 0; i2 < ny; i2++) {
			t = float(i2 - ldpy); 
			const float yy = t*t + xx;
			bool first = true;
			for (int i3 = 0; i3 < nz; i3++) {
				t = float(i3 - ldpx); 
				const float rc = t*t + yy;
				if (first) {
					// first pixel on this line
					if (rc > r) continue;
					first = false;
					nn++;
					if (ipcube != NULL) {
						ipcube[nn].start = i3;
						ipcube[nn].end = i3;
						ipcube[nn].loc[0] = i3 - ldpx;
						ipcube[nn].loc[1] = i2 - ldpy;
						ipcube[nn].loc[2] = i1 - ldpz;
					}
				} else {
					// second or later pixel on this line
					if (ipcube != NULL) {
						if (rc <= r) ipcube[nn].end = i3;
					}
				}
			}
		}
	}
}

EMData *PawelProjector::project3d(EMData * image) const
{
	if (!image) {
		return 0;
	}
	int ri;
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	int dim = Util::get_min(nx,ny,nz);
	if (nz == 1) {
		LOGERR("The PawelProjector needs a volume!");
		return 0;
	}
	Vec3i origin(0,0,0);
	// If a sensible origin isn't passed in, choose the middle of
	// the cube.
	if (params.has_key("origin_x")) {origin[0] = params["origin_x"];} 
	else {origin[0] = nx/2;}
	if (params.has_key("origin_y")) {origin[1] = params["origin_y"];} 
	else {origin[1] = ny/2;}
	if (params.has_key("origin_z")) {origin[2] = params["origin_z"];} 
	else {origin[2] = nz/2;}

	if (params.has_key("radius")) {ri = params["radius"];} 
	else {ri = dim/2 - 1;}

	// Determine the number of rows (x-lines) within the radius
	int nn = -1;
	prepcubes(nx, ny, nz, ri, origin, nn); 
	// nn is now the number of rows-1 within the radius
	// so we can create and fill the ipcubes
	IPCube* ipcube = new IPCube[nn+1];
	prepcubes(nx, ny, nz, ri, origin, nn, ipcube);

	// Do we have a list of angles?
	vector<float> anglelist = params["anglelist"];
	int nangles = anglelist.size() / 3;

	// For the moment, let's assume the user wants to use either
	// EMAN or SPIDER Euler angles, but nothing else.
	string angletype = params["angletype"].to_str();
	Transform3D::EulerType eulertype;
	if (angletype == "SPIDER") {
		eulertype = Transform3D::SPIDER;
		if (nangles == 0) {
			// a single SPIDER angle was passed in
			float phi = params["phi"];
			float theta = params["theta"];
			float psi = params["psi"];
			anglelist.push_back(phi);
			anglelist.push_back(theta);
			anglelist.push_back(psi);
			nangles = 1;
		}
	} else if (angletype == "EMAN") {
		eulertype = Transform3D::EMAN;
		if (nangles == 0) {
			// a single EMAN angle was passed in
			float az = params["az"];
			float alt = params["alt"];
			float phi = params["phi"];
			anglelist.push_back(az);
			anglelist.push_back(alt);
			anglelist.push_back(phi);
			nangles = 1;
		}
	} else {
		LOGERR("Only SPIDER and EMAN Euler angles supported");
		return 0;
	}
	// initialize return object
	EMData* ret = new EMData();
	ret->set_size(nx, ny, nangles);
	ret->to_zero();

	// loop over sets of angles
	for (int ia = 0; ia < nangles; ia++) {
		int indx = 3*ia;
		Transform3D rotation(eulertype, float(anglelist[indx]*dgr_to_rad),
				           float(anglelist[indx+1]*dgr_to_rad), float(anglelist[indx+2]*dgr_to_rad));
		if (2*(ri+1)+1 > dim) {
			// Must check x and y boundaries
			for (int i = 0 ; i <= nn; i++) {
				int k = ipcube[i].loc[1] + origin[1];
				Vec3f vb = ipcube[i].loc*rotation + origin;
				for (int j = ipcube[i].start; j <= ipcube[i].end; j++) {
					// check for pixels out-of-bounds
					int iox = int(vb[0]);
					if ((iox >= 0) & (iox < nx-1)) {
						int ioy = int(vb[1]);
						if ((ioy >= 0) & (ioy < ny-1)) {
							int ioz = int(vb[2]);
							if ((ioz >= 0) & (ioz < nz-1)) {
								// real work for pixels in bounds
								float dx = vb[0] - iox;
								float dy = vb[1] - ioy;
								float dz = vb[2] - ioz;
								float a1 = (*image)(iox,ioy,ioz);
								float a2 = (*image)(iox+1,ioy,ioz) - a1;
								float a3 = (*image)(iox,ioy+1,ioz) - a1;
								float a4 = (*image)(iox,ioy,ioz+1) - a1;
								float a5 = -a2 -(*image)(iox,ioy+1,ioz) 
									+ (*image)(iox+1,ioy+1,ioz);
								float a61 = -(*image)(iox,ioy,ioz+1) 
									+ (*image)(iox+1,ioy,ioz+1);
								float a6 = -a2 + a61;
								float a7 = -a3 - (*image)(iox,ioy,ioz+1)
									+ (*image)(iox,ioy+1,ioz+1);
								float a8 = -a5 - a61 - (*image)(iox,ioy+1,ioz+1)
									+ (*image)(iox+1,ioy+1,ioz+1);
								(*ret)(j,k,ia) += a1 + dz*(a4 + a6*dx  
				                                        + (a7 + a8*dx)*dy)
									+ a3*dy + dx*(a2 + a5*dy);
							}
						}
					}
					vb += rotation.get_matrix3_row(0);
				}
			}

		} else {
			// No need to check x and y boundaries
			for (int i = 0 ; i <= nn; i++) {
				int k = ipcube[i].loc[1] + origin[1];
				Vec3f vb = ipcube[i].loc*rotation + origin;
				for (int j = ipcube[i].start; j <= ipcube[i].end; j++) {
					int iox = int(vb[0]);
					int ioy = int(vb[1]);
					int ioz = int(vb[2]);
					float dx = vb[0] - iox;
					float dy = vb[1] - ioy;
					float dz = vb[2] - ioz;
					float a1 = (*image)(iox,ioy,ioz);
					float a2 = (*image)(iox+1,ioy,ioz) - a1;
					float a3 = (*image)(iox,ioy+1,ioz) - a1;
					float a4 = (*image)(iox,ioy,ioz+1) - a1;
					float a5 = -a2 -(*image)(iox,ioy+1,ioz) 
						+ (*image)(iox+1,ioy+1,ioz);
					float a61 = -(*image)(iox,ioy,ioz+1) 
						+ (*image)(iox+1,ioy,ioz+1);
					float a6 = -a2 + a61;
					float a7 = -a3 - (*image)(iox,ioy,ioz+1)
						+ (*image)(iox,ioy+1,ioz+1);
					float a8 = -a5 - a61 - (*image)(iox,ioy+1,ioz+1)
						+ (*image)(iox+1,ioy+1,ioz+1);
					(*ret)(j,k,ia) += a1 + dz*(a4 + a6*dx  
							+ (a7 + a8*dx)*dy)
						+ a3*dy + dx*(a2 + a5*dy);
					vb += rotation.get_matrix3_row(0);
				}
			}
		}
	}
	ret->done_data();
	ret->update();
	if( ipcube )
	{
		delete [] ipcube;
		ipcube = 0;
	}
	return ret;
}


EMData *SimpleIsoSurfaceProjector::project3d(EMData * image) const
{
	float alt = params["alt"];
	float az = params["az"];
	float phi = params["phi"];
	float threshold = params["threshold"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	image->rotate(az, alt, phi);

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
				dat[k] = 2.0f;
			}
			else {
				v0[0] = 1.0f / (float)hypot(1, slx);
				v0[1] = 0;
				v0[2] = slx / (float)hypot(1, slx);

				v1[0] = 0;
				v1[1] = 1.0f / (float)hypot(1, sly);
				v1[2] = sly / (float)hypot(1, sly);

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




EMData *StandardFastProjector::project3d(EMData * image) const
{
	float alt = params["alt"];
	float az = params["az"];
	float phi = params["phi"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	image->rotate(az, alt, phi);

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

EMData *StandardProjector::project3d(EMData * image) const
{
	float alt = float(params["alt"]) * static_cast<float>(dgr_to_rad);
	float az = float(params["az"]) * static_cast<float>(dgr_to_rad);
	float phi = float(params["phi"]) * static_cast<float>(dgr_to_rad);

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	Transform3D r(Transform3D::EMAN, az, alt, phi);
	r.transpose();
	int xy = nx * ny;

	EMData *proj = new EMData();
	proj->set_size(nx, ny, 1);

	float *sdata = image->get_data();
	float *ddata = proj->get_data();
	for (int k = -nz / 2; k < nz - nz / 2; k++) {
		int l = 0;
		for (int j = -ny / 2; j < ny - ny / 2; j++) {
			ddata[l]=0;
			for (int i = -nx / 2; i < nx - nx / 2; i++,l++) {
				float x2 = (float)(r[0][0] * i + r[0][1] * j + r[0][2] * k + nx / 2);
				float y2 = (float)(r[1][0] * i + r[1][1] * j + r[1][2] * k + ny / 2);
				float z2 = (float)(r[2][0] * i + r[2][1] * j + r[2][2] * k + nz / 2);

				if (x2 >= 0 && y2 >= 0 && z2 >= 0 && x2 < (nx - 1) && y2 < (ny - 1)
					&& z2 < (nz - 1)) {
					float x = (float)Util::fast_floor(x2);
					float y = (float)Util::fast_floor(y2);
					float z = (float)Util::fast_floor(z2);

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

EMData *FourierGriddingProjector::project3d(EMData * image) const
{
	if (!image) {
		return 0;
	}
	if (3 != image->get_ndim()) 
		throw ImageDimensionException(
				"FourierGriddingProjector needs a 3-D volume");
	if (image->is_complex())
		throw ImageFormatException(
				"FourierGriddingProjector requires a real volume");
	const int npad = params.has_key("npad") ? int(params["npad"]) : 2;
	const int nx = image->get_xsize();
	const int ny = image->get_ysize();
	const int nz = image->get_zsize();
	if (nx != ny || nx != nz)
		throw ImageDimensionException(
				"FourierGriddingProjector requires nx==ny==nz");
	const int m = Util::get_min(nx,ny,nz);
	const int n = m*npad;

	const int K = params["kb_K"];
	const float alpha = params["kb_alpha"];
	Util::KaiserBessel kb(alpha, K, m/2,K/(2.*n),n);
	// divide out gridding weights
	image->divkbsinh(kb);
	// pad and center volume, then FFT and multiply by (-1)**(i+j+k)
	EMData* imgft = image->pad_fft(npad);
	imgft->center_padded();
	imgft->do_fft_inplace();
	imgft->center_origin_fft();
	imgft->fft_shuffle();

	// Do we have a list of angles?
	int nangles = 0;
	vector<float> anglelist;
	if (params.has_key("anglelist")) {
		anglelist = params["anglelist"];
		nangles = anglelist.size() / 3;
	}

	// For the moment, let's assume the user wants to use either
	// EMAN or SPIDER Euler angles, but nothing else.
	string angletype = params["angletype"].to_str();
	Transform3D::EulerType eulertype;
	if (angletype == "SPIDER") {
		eulertype = Transform3D::SPIDER;
		if (nangles == 0) {
			// a single SPIDER angle was passed in
			float phi = params["phi"];
			float theta = params["theta"];
			float psi = params["psi"];
			anglelist.push_back(phi);
			anglelist.push_back(theta);
			anglelist.push_back(psi);
			nangles = 1;
		}
	} else if (angletype == "EMAN") {
		eulertype = Transform3D::EMAN;
		if (nangles == 0) {
			// a single EMAN angle was passed in
			float az = params["az"];
			float alt = params["alt"];
			float phi = params["phi"];
			anglelist.push_back(az);
			anglelist.push_back(alt);
			anglelist.push_back(phi);
			nangles = 1;
		}
	} else 
		throw InvalidValueException(0,
				"Only SPIDER and EMAN Euler angles currently supported");
	// initialize return object
	EMData* ret = new EMData();
	ret->set_size(nx, ny, nangles);
	ret->to_zero();

	// loop over sets of angles
	for (int ia = 0; ia < nangles; ia++) {
		int indx = 3*ia;
		Transform3D tf(eulertype, float(anglelist[indx]*dgr_to_rad),
				float(anglelist[indx+1]*dgr_to_rad), 
				float(anglelist[indx+2]*dgr_to_rad));
		EMData* proj = imgft->extractplane(tf, kb);
		if (proj->is_shuffled()) proj->fft_shuffle();
		proj->center_origin_fft();
		proj->do_ift_inplace();
		EMData* winproj = proj->window_center(m);
		delete proj;
		for (int iy=0; iy < ny; iy++)
			for (int ix=0; ix < nx; ix++) 
				(*ret)(ix,iy,ia) = (*winproj)(ix,iy);
		delete winproj;
	}
	delete imgft;
	ret->done_data();
	ret->update();
	return ret;
}



void EMAN::dump_projectors()
{
	dump_factory < Projector > ();
}

map<string, vector<string> > EMAN::dump_projectors_list()
{
	return dump_factory_list < Projector > ();
}

/* vim: set ts=4 noet nospell: */
