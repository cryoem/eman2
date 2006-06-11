/**
 * $Id$
 */
#include "reconstructor.h"
#include "interp.h"
#include "ctf.h"

using namespace EMAN;

template <> Factory < Reconstructor >::Factory()
{
	force_add(&FourierReconstructor::NEW);
	force_add(&WienerFourierReconstructor::NEW);
	force_add(&BackProjectionReconstructor::NEW);
	force_add(&nn4Reconstructor::NEW);
	force_add(&nn4_ctfReconstructor::NEW);
}

FourierReconstructor::FourierReconstructor()
:	image(0), nx(0), ny(0), nz(0)
{
}

FourierReconstructor::~FourierReconstructor()
{
	if (image) {
		delete image;
		image = 0;
	}
}

void FourierReconstructor::setup()
{
	int size = params["size"];
	image = new EMData();
	image->set_size(size + 2, size, size);
	image->set_complex(true);
	image->set_ri(true);

	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();

	int n = nx * ny * nz;
	float *rdata = image->get_data();

	for (int i = 0; i < n; i += 2) {
		float f = Util::get_frand(0.0, 2.0 * M_PI);
		rdata[i] = 1.0e-10f * sin(f);
		rdata[i + 1] = 1.0e-10f * cos(f);
	}
	image->done_data();

	tmp_data = new EMData();
	tmp_data->set_size(size + 2, size, size);
	tmp_data->to_zero();
}

int FourierReconstructor::insert_slice(EMData * slice, const Transform3D & euler)
{
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int mode = params["mode"];
	float weight = params["weight"];

	if (!slice->is_complex()) {
		LOGERR("Only complex slice can be inserted.");
		return 1;
	}

	float *gimx = 0;
	if (mode == 5) {
		gimx = Interp::get_gimx();
	}

	int nxy = nx * ny;

	int off[8];
	if (mode == 2) {
		off[0] = 0;
		off[1] = 2;
		off[2] = nx;
		off[3] = nx + 2;
		off[4] = nxy;
		off[5] = nxy + 2;
		off[6] = nxy + nx;
		off[7] = nxy + nx + 2;
	}

	float *norm = tmp_data->get_data();
	float *dat = slice->get_data();
	float *rdata = image->get_data();
	
	int rl = Util::square(ny / 2 - 1);
	float dt[2];
	float g[8];

	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++) {
			if ((x * x + Util::square(y - ny / 2)) >= rl)
				continue;

			float xx = (float) (x * euler[0][0] + (y - ny / 2) * euler[1][0]);
			float yy = (float) (x * euler[0][1] + (y - ny / 2) * euler[1][1]);
			float zz = (float) (x * euler[0][2] + (y - ny / 2) * euler[1][2]);
			float cc = 1;

			if (xx < 0) {
				xx = -xx;
				yy = -yy;
				zz = -zz;
				cc = -1.0;
			}

			yy += ny / 2;
			zz += nz / 2;

			dt[0] = dat[x * 2 + y * nx];
			dt[1] = cc * dat[x * 2 + 1 + y * nx];

			int x0 = 0;
			int y0 = 0;
			int z0 = 0;
			int i = 0;
			int l = 0;
			float dx = 0;
			float dy = 0;
			float dz = 0;

			int mx0 = 0;
			int my0 = 0;
			int mz0 = 0;

			switch (mode) {
			case 1:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				rdata[x0 + y0 * nx + z0 * nxy] += weight * dt[0];
				rdata[x0 + y0 * nx + z0 * nxy + 1] += weight * dt[1];
				norm[x0 + y0 * nx + z0 * nxy] += weight;
				break;

			case 2:
				x0 = (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				dx = xx - x0;
				dy = yy - y0;
				dz = zz - z0;

				if (x0 > nx - 2 || y0 > ny - 1 || z0 > nz - 1) {
					break;
				}

				i = (int) (x0 * 2 + y0 * nx + z0 * nxy);

				g[0] = Util::agauss(1, dx, dy, dz, EMConsts::I2G);
				g[1] = Util::agauss(1, 1 - dx, dy, dz, EMConsts::I2G);
				g[2] = Util::agauss(1, dx, 1 - dy, dz, EMConsts::I2G);
				g[3] = Util::agauss(1, 1 - dx, 1 - dy, dz, EMConsts::I2G);
				g[4] = Util::agauss(1, dx, dy, 1 - dz, EMConsts::I2G);
				g[5] = Util::agauss(1, 1 - dx, dy, 1 - dz, EMConsts::I2G);
				g[6] = Util::agauss(1, dx, 1 - dy, 1 - dz, EMConsts::I2G);
				g[7] = Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, EMConsts::I2G);

				for (int j = 0; j < 8; j++) {
					int k = i + off[j];
					rdata[k] += weight * dt[0] * g[j];
					rdata[k + 1] += weight * dt[1] * g[j];
					norm[k] += weight * g[j];
					norm[k + 1] += weight * g[j] * dt[0] * dt[0] * dt[1] * dt[1];
				}

				break;

			case 3:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}

				for (int k = z0 - 1; k <= z0 + 1; k++) {
					for (int j = y0 - 1; j <= y0 + 1; j++) {
						for (int i = l; i <= x0 + 2; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I3G);

							rdata[i + j * nx + k * nxy] += weight * gg * dt[0];
							rdata[i + j * nx + k * nxy + 1] += weight * gg * dt[1];
							norm[i + j * nx + k * nxy] += weight * gg;
						}
					}
				}
				break;

			case 4:
				x0 = 2 * (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}

				for (int k = z0 - 1; k <= z0 + 2; k++) {
					for (int j = y0 - 1; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I4G);

							rdata[i + j * nx + k * nxy] += weight * gg * dt[0];
							rdata[i + j * nx + k * nxy + 1] += weight * gg * dt[1];
							norm[i + j * nx + k * nxy] += weight * gg;
						}
					}
				}
				break;

			case 5:
				x0 = (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				mx0 = -(int) floor((xx - x0) * 39.0f + 0.5f) - 78;
				my0 = -(int) floor((yy - y0) * 39.0f + 0.5f) - 78;
				mz0 = -(int) floor((zz - z0) * 39.0f + 0.5f) - 78;

				x0 *= 2;

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = 0;
				if (x0 == 0) {
					mx0 += 78;
				}
				else if (x0 == 2) {
					mx0 += 39;
				}
				else {
					l = x0 - 4;
				}

				for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
					for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
						for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
							int ii = i + j * nx + k * nxy;
							float gg = weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;

					x0 = (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);

					int mx0 = -(int) floor((xx - x0) * 39.0f + 0.5f);
					x0 *= 2;

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
						for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
							for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
								int ii = i + j * nx + k * nxy;
								float gg =
									weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];
								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;

			case 6:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				l = x0 - 4;
				if (x0 <= 2) {
					l = 0;
				}

				for (int k = z0 - 2; k <= z0 + 2; k++) {
					for (int j = y0 - 2; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							int ii = i + j * nx + k * nxy;
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = weight * exp(-r / EMConsts::I5G);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;

					x0 = 2 * (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; k++) {
						for (int j = y0 - 2; j <= y0 + 2; j++) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								int ii = i + j * nx + k * nxy;
								float r = Util::hypot3((float) i / 2 - xx, (float) j - yy,
													   (float) k - zz);
								float gg = weight * exp(-r / EMConsts::I5G);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;

			case 7:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				l = x0 - 4;
				if (x0 <= 2)
					l = 0;

				for (int k = z0 - 2; k <= z0 + 2; k++) {
					for (int j = y0 - 2; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							int ii = i + j * nx + k * nxy;
							float r =
								sqrt(Util::
									 hypot3((float) i / 2 - xx, (float) j - yy, (float) k - zz));
							float gg = weight * Interp::hyperg(r);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = 2 * (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; k++) {
						for (int j = y0 - 2; j <= y0 + 2; j++) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								int ii = i + j * nx + k * nxy;
								float r = sqrt(Util::hypot3((float) i / 2 - xx, (float) j - yy,
															(float) k - zz));
								float gg = weight * Interp::hyperg(r);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
			default:
				LOGERR("no such insert slice mode: '%d'", mode);
				return 1;
			}

		}
	}

	image->done_data();
	tmp_data->done_data();
	slice->done_data();
	slice->update();

	return 0;
}

EMData *FourierReconstructor::finish()
{
	int dlog = params["dlog"];
	float *norm = tmp_data->get_data();
	float *rdata = image->get_data();

	if (dlog) {
		for (int i = 0; i < nx * ny * nz; i += 2) {
			float d = norm[i];
			if (d == 0) {
				rdata[i] = 0;
				rdata[i + 1] = 0;
			}
			else {
				float in = norm[i + 1] / norm[i];
				float cin = Util::square_sum(rdata[i], rdata[i + 1]);
				rdata[i] *= sqrt(in / cin);
				rdata[i + 1] *= sqrt(in / cin);
			}
		}
	}
	else {
		for (int i = 0; i < nx * ny * nz; i += 2) {
			float d = norm[i];
			if (d == 0) {
				rdata[i] = 0;
				rdata[i + 1] = 0;
			}
			else {
				rdata[i] /= d;
				rdata[i + 1] /= d;
			}
		}
	}

	if( tmp_data ) {
		delete tmp_data;
		norm = 0;
	}
	image->done_data();
	return image;
}

WienerFourierReconstructor::WienerFourierReconstructor()
:	image(0), nx(0), ny(0), nz(0)
{
}

WienerFourierReconstructor::~WienerFourierReconstructor()
{
	if (image) {
		delete image;
		image = 0;
	}
}

void WienerFourierReconstructor::setup()
{
	int size = params["size"];
	image = new EMData();
	image->set_size(size + 2, size, size);
	image->set_complex(true);
	image->set_ri(true);

	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();

	int n = nx * ny * nz;
	float *rdata = image->get_data();

	for (int i = 0; i < n; i += 2) {
		float f = Util::get_frand(0.0, 2.0 * M_PI);
		rdata[i] = 1.0e-10f * sin(f);
		rdata[i + 1] = 1.0e-10f * cos(f);
	}
	image->done_data();

	tmp_data = new EMData();
	tmp_data->set_size(size + 1, size, size);
}


EMData *WienerFourierReconstructor::finish()
{
	float *norm = tmp_data->get_data();
	float *rdata = image->get_data();

	for (int i = 0; i < nx * ny * nz; i += 2) {
		float d = norm[i];
		if (d == 0) {
			rdata[i] = 0;
			rdata[i + 1] = 0;
		}
		else {
			float w = 1 + 1 / d;
			rdata[i] /= d * w;
			rdata[i + 1] /= d * w;
		}
	}

	if( tmp_data ) {
		delete tmp_data;
		tmp_data = 0;
	}
	image->done_data();
	return image;
}


int WienerFourierReconstructor::insert_slice(EMData * slice, const Transform3D & euler)
{
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int mode = params["mode"];
	float padratio = params["padratio"];
	vector < float >snr = params["snr"];

	if (!slice->is_complex()) {
		LOGERR("Only complex slice can be inserted.");
		return 1;
	}
	float *gimx = 0;
	if (mode == 5) {
		gimx = Interp::get_gimx();
	}

	int nxy = nx * ny;
	int off[8];
	if (mode == 2) {
		off[0] = 0;
		off[1] = 2;
		off[2] = nx;
		off[3] = nx + 2;
		off[4] = nxy;
		off[5] = nxy + 2;
		off[6] = nxy + nx;
		off[7] = nxy + nx + 2;
	}

	float *norm = tmp_data->get_data();
	float *dat = slice->get_data();
	float *rdata = image->get_data();

	int rl = Util::square(ny / 2 - 1);
	float dt[2];
	float g[8];

	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++) {
			if ((x * x + Util::square(y - ny / 2)) >= rl) {
				continue;
			}

			int r = Util::round((float)hypot(x, (float) y - ny / 2) * Ctf::CTFOS / padratio);
			if (r >= Ctf::CTFOS * ny / 2) {
				r = Ctf::CTFOS * ny / 2 - 1;
			}

			float weight = snr[r];

			float xx = (x * euler[0][0] + (y - ny / 2) * euler[0][1]);
			float yy = (x * euler[1][0] + (y - ny / 2) * euler[1][1]);
			float zz = (x * euler[2][0] + (y - ny / 2) * euler[2][1]);
			float cc = 1;

			if (xx < 0) {
				xx = -xx;
				yy = -yy;
				zz = -zz;
				cc = -1.0;
			}

			yy += ny / 2;
			zz += nz / 2;

			dt[0] = dat[x * 2 + y * nx] * (1 + 1.0f / weight);
			dt[1] = cc * dat[x * 2 + 1 + y * nx] * (1 + 1.0f / weight);

			int x0 = 0;
			int y0 = 0;
			int z0 = 0;
			int i = 0;
			int l = 0;
			float dx = 0;
			float dy = 0;
			float dz = 0;

			int mx0 = 0;
			int my0 = 0;
			int mz0 = 0;

			switch (mode) {
			case 1:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				rdata[x0 + y0 * nx + z0 * nxy] += weight * dt[0];
				rdata[x0 + y0 * nx + z0 * nxy + 1] += weight * dt[1];
				norm[x0 + y0 * nx + z0 * nxy] += weight;
				break;

			case 2:
				x0 = (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				dx = xx - x0;
				dy = yy - y0;
				dz = zz - z0;

				weight /= (float)pow((float)(EMConsts::I2G * M_PI), 1.5f);

				if (x0 > nx - 2 || y0 > ny - 1 || z0 > nz - 1) {
					break;
				}

				i = (int) (x0 * 2 + y0 * nx + z0 * nxy);


				g[0] = Util::agauss(1, dx, dy, dz, EMConsts::I2G);
				g[1] = Util::agauss(1, 1 - dx, dy, dz, EMConsts::I2G);
				g[2] = Util::agauss(1, dx, 1 - dy, dz, EMConsts::I2G);
				g[3] = Util::agauss(1, 1 - dx, 1 - dy, dz, EMConsts::I2G);
				g[4] = Util::agauss(1, dx, dy, 1 - dz, EMConsts::I2G);
				g[5] = Util::agauss(1, 1 - dx, dy, 1 - dz, EMConsts::I2G);
				g[6] = Util::agauss(1, dx, 1 - dy, 1 - dz, EMConsts::I2G);
				g[7] = Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, EMConsts::I2G);

				for (int j = 0; j < 8; j++) {
					int k = i + off[j];
					rdata[k] += weight * dt[0] * g[j];
					rdata[k + 1] += weight * dt[1] * g[j];
					norm[k] += weight * g[j];
				}

				break;
			case 3:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				weight /= (float)pow((float)(EMConsts::I3G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}
				for (int k = z0 - 1; k <= z0 + 1; k++) {
					for (int j = y0 - 1; j <= y0 + 1; j++) {
						for (int i = l; i <= x0 + 2; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I3G);

							rdata[i + j * nx + k * nxy] += weight * gg * dt[0];
							rdata[i + j * nx + k * nxy + 1] += weight * gg * dt[1];
							norm[i + j * nx + k * nxy] += weight * gg;
						}
					}
				}
				break;

			case 4:
				x0 = 2 * (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				weight /= (float)pow((float)(EMConsts::I4G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}
				for (int k = z0 - 1; k <= z0 + 2; k++) {
					for (int j = y0 - 1; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I4G);

							rdata[i + j * nx + k * nxy] += weight * gg * dt[0];
							rdata[i + j * nx + k * nxy + 1] += weight * gg * dt[1];
							norm[i + j * nx + k * nxy] += weight * gg;
						}
					}
				}
				break;

			case 5:
				x0 = (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				weight /= (float)pow((float)(EMConsts::I5G * M_PI), 1.5f);

				mx0 = -(int) floor((xx - x0) * 39.0f + 0.5) - 78;
				my0 = -(int) floor((yy - y0) * 39.0f + 0.5) - 78;
				mz0 = -(int) floor((zz - z0) * 39.0f + 0.5) - 78;
				x0 *= 2;

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 == 0) {
					l = 0;
					mx0 += 78;
				}
				else if (x0 == 2) {
					l = 0;
					mx0 += 39;
				}
				else
					l = x0 - 4;
				for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
					for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
						for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
							int ii = i + j * nx + k * nxy;
							float gg = weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);
					int mx0 = -(int) floor((xx - x0) * 39.0f + .5);
					x0 *= 2;

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
						for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
							for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
								int ii = i + j * nx + k * nxy;
								float gg =
									weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;
				// mode 6 is now mode 5 without the fast interpolation
			case 6:
				x0 = 2 * (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				weight /= (float)pow((float)(EMConsts::I5G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 <= 2)
					l = 0;
				else
					l = x0 - 4;
				for (int k = z0 - 2; k <= z0 + 2; k++) {
					for (int j = y0 - 2; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							int ii = i + j * nx + k * nxy;
							float r = Util::hypot3((float) i / 2 - xx, (float) j - yy,
												   (float) k - zz);
							float gg = weight * exp(-r / EMConsts::I5G);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = 2 * (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; k++) {
						for (int j = y0 - 2; j <= y0 + 2; j++) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								int ii = i + j * nx + k * nxy;
								float r = Util::hypot3((float) i / 2 - xx, (float) j - yy,
													   (float) k - zz);
								float gg = weight * exp(-r / EMConsts::I5G);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];	// note the -, complex conj.
								norm[ii] += gg;
							}
						}
					}
				}
				break;

			case 7:
				x0 = 2 * (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 <= 2)
					l = 0;
				else
					l = x0 - 4;
				for (int k = z0 - 2; k <= z0 + 2; k++) {
					for (int j = y0 - 2; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							int ii = i + j * nx + k * nxy;
							float r = (float)sqrt(Util::hypot3((float) i / 2 - xx,
															   (float) j - yy,
															   (float) k - zz));
							float gg = weight * Interp::hyperg(r);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = 2 * (int) floor(xx + .5);
					y0 = (int) floor(yy + .5);
					z0 = (int) floor(zz + .5);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; k++) {
						for (int j = y0 - 2; j <= y0 + 2; j++) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								int ii = i + j * nx + k * nxy;
								float r = sqrt(Util::hypot3((float) i / 2 - xx, (float) j - yy,
															(float) k - zz));
								float gg = weight * Interp::hyperg(r);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;
			}

		}
	}

	image->done_data();
	tmp_data->done_data();
	slice->done_data();
	slice->update();

	return 0;
}

BackProjectionReconstructor::BackProjectionReconstructor()
:	image(0), nx(0), ny(0), nz(0)
{
}

BackProjectionReconstructor::~BackProjectionReconstructor()
{
	if (image) {
		delete image;
		image = 0;
	}
}

void BackProjectionReconstructor::setup()
{
	int size = params["size"];
	image = new EMData();
	image->set_size(size, size, size);
	nx = size;
	ny = size;
	nz = size;
}

int BackProjectionReconstructor::insert_slice(EMData * slice, const Transform3D &)
{
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	if (slice->get_xsize() != slice->get_ysize() || slice->get_xsize() != nx) {
		LOGERR("try to insert improve size slice");
		return 1;
	}

	float weight = params["weight"];

	EMData *slice_copy = slice->copy();
	slice_copy->mult(weight);

	EMData *tmp = new EMData();
	tmp->set_size(nx, ny, nz);

	float *slice_data = slice_copy->get_data();
	float *tmp_data = tmp->get_data();
	int nxy = nx * ny;
	size_t nxy_size = nxy * sizeof(float);

	for (int i = 0; i < nz; i++) {
		memcpy(&tmp_data[nxy * i], slice_data, nxy_size);
	}

	tmp->done_data();

	Dict slice_euler = slice->get_transform().get_rotation(Transform3D::EMAN);
	tmp->rotate(-(float)slice_euler["alt"], -(float)slice_euler["az"], -(float)slice_euler["phi"]);

	image->add(*tmp);
	if( slice_copy )
	{
		delete slice_copy;
		slice_copy = 0;
	}

	if( tmp )
	{
		delete tmp;
		tmp = 0;
	}
	
	return 0;
}

EMData *BackProjectionReconstructor::finish()
{
	return image;
}

nn4Reconstructor::nn4Reconstructor() 
	: v(NULL) {}

nn4Reconstructor::~nn4Reconstructor()
{
	if (v) {
		delete v;
		v = NULL;
	}
}

void nn4Reconstructor::setup() {
	int nsize = params["size"];
	vnx = vny = vnz = nsize;
	npad = params["npad"];
	vnxp = nsize*npad;
	vnyp = nsize*npad;
	vnzp = nsize*npad;
	vnxc = vnxp/2;
	buildFFTVolume();
	buildNormVolume();
	try {
		symmetry = params["symmetry"].to_str();
		if ("" == symmetry) symmetry = "c1";
	}
	catch(_NotExistingObjectException) {
		symmetry = "c1";
	}
	nsym = Transform3D::get_nsym(symmetry);
}

void nn4Reconstructor::buildFFTVolume() {
	v = new EMData;
	int offset = 2 - vnxp%2;
	v->set_size(vnxp+offset,vnyp,vnzp);
	v->set_nxc(vnxp/2);
	v->set_complex(true);
	v->set_ri(true);
	v->set_fftpad(true);
	v->set_attr("npad", npad);
	v->to_zero();
	v->set_array_offsets(0,1,1);
}

void nn4Reconstructor::buildNormVolume() {
	nrptr = new EMArray<int>(vnxc+1,vnyp,vnzp);
	nrptr->set_array_offsets(0,1,1);
	for (int iz = 1; iz <= vnzp; iz++) 
		for (int iy = 1; iy <= vnyp; iy++) 
			for (int ix = 0; ix <= vnxc; ix++) 
				(*nrptr)(ix,iy,iz) = 0;
}

int nn4Reconstructor::insert_slice(EMData* slice, const Transform3D& t) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if ((slice->get_xsize() != slice->get_ysize()) 
		|| slice->get_xsize() != vnx) {
		// FIXME: Why doesn't this throw an exception?
		LOGERR("Tried to insert a slice that is the wrong size.");
		return 1;
	}
	// process 2-d slice -- subtract the average outside of the circle, zero-pad, fft extend, and fft
	EMData* temp = slice->average_circ_sub();
	// Need to use zeropad_ntimes instead of pad_fft here for zero padding
	// because only the former centers the original image in the 
	// larger area.  FIXME!
	EMData* zeropadded = temp->zeropad_ntimes(npad);
	if( temp ) {delete temp; temp = 0;}
	EMData* padfftslice = zeropadded->do_fft();
	if( zeropadded ) { delete zeropadded; zeropadded = 0;}
	padfftslice->center_origin_fft();
	// insert slice for all symmetry related positions
	for (int isym=0; isym < nsym; isym++) {
		Transform3D tsym = t.get_sym(symmetry, isym);
		v->nn(*nrptr, padfftslice, tsym);
	}
	if( padfftslice ) { delete padfftslice;
		padfftslice = 0;}
	return 0;
}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
EMData* nn4Reconstructor::finish() {
	EMArray<int>& nr = *nrptr;
	//v->symplane0(nr);
	// normalize
	for (int iz = 1; iz <= vnzp; iz++) {
		for (int iy = 1; iy <= vnyp; iy++) {
			for (int ix = 0; ix <= vnxc; ix++) {//std::cout<<"  "<<iz<<"  "<<iy<<"  "<<ix<<"  "<<nr(ix,iy,iz)<<"  "<<(-2*((ix+iy+iz)%2)+1)<<"  "<<(*v)(ix,iy,iz)<<std::endl;
				if (nr(ix,iy,iz) > 0) {//(*v) should be treated as complex!!
					float  tmp = (-2*((ix+iy+iz)%2)+1)/float(nr(ix,iy,iz));
					(*v)(2*ix,iy,iz) *= tmp;
					(*v)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}
	// back fft
	v->do_ift_inplace();
	EMData* w = v->window_center(vnx);
	if( v ) { delete v; v = 0;}
	float *tw = w->get_data();
	//  mask and subtract circumference average
	int ix = w->get_xsize(),iy = w->get_ysize(),iz = w->get_zsize();
	int L2 = (ix/2)*(ix/2);
	int L2P = (ix/2-1)*(ix/2-1);
	int IP = ix/2+1;
	float  TNR = 0.0f;
	int m = 0;
	for (int k = 1; k <= iz; k++) {
		for (int j = 1; j <= iy; j++) {
			for (int i = 1; i <= ix; i++) {
				int LR = (k-IP)*(k-IP)+(j-IP)*(j-IP)+(i-IP)*(i-IP);
				if (LR<=L2) {
					if(LR >= L2P && LR <= L2) {
						TNR += tw(i,j,k);
						m++;
					}
				}
			}
		}
	}
	TNR /=float(m);
	for (int k = 1; k <= iz; k++) {
		for (int j = 1; j <= iy; j++) {
			for (int i = 1; i <= ix; i++) {
				int LR = (k-IP)*(k-IP)+(j-IP)*(j-IP)+(i-IP)*(i-IP);
				if (LR<=L2) tw(i,j,k) -= TNR; else tw(i,j,k) = 0.0f;
			}
		}
	}
	// clean up
	if( nrptr ) { delete nrptr; nrptr = 0;}
	return w;
}
#undef  tw
//####################################################################################
//** nn4 ctf reconstructor 

nn4_ctfReconstructor::nn4_ctfReconstructor() 
	: v(NULL) {}

nn4_ctfReconstructor::~nn4_ctfReconstructor()
{
	if (v) {delete v;v = NULL;}
}

void nn4_ctfReconstructor::setup() {
	int nsize = params["size"];
	vnx = vny = vnz = nsize;
	npad = params["npad"];
	vnxp = nsize*npad;
	vnyp = nsize*npad;
	vnzp = nsize*npad;
	vnxc = vnxp/2;
	buildFFTVolume();
	buildNormVolume();
	try {
		symmetry = params["symmetry"].to_str();
		if ("" == symmetry) symmetry = "c1";
	}
	catch(_NotExistingObjectException) {
		symmetry = "c1";
	}
	nsym = Transform3D::get_nsym(symmetry);
}

void nn4_ctfReconstructor::buildFFTVolume() {
	v = new EMData;
	int offset = 2 - vnxp%2;
	v->set_size(vnxp+offset,vnyp,vnzp);
	v->set_nxc(vnxp/2);
	v->set_complex(true);
	v->set_ri(true);
	v->set_fftpad(true);
	v->set_attr("npad", npad);
	v->to_zero();
	v->set_array_offsets(0,1,1);
}

void nn4_ctfReconstructor::buildNormVolume() {
	wptr = new EMArray<float>(vnxc+1,vnyp,vnzp);
	wptr->set_array_offsets(0,1,1);
	for (int iz = 1; iz <= vnzp; iz++) 
		for (int iy = 1; iy <= vnyp; iy++) 
			for (int ix = 0; ix <= vnxc; ix++) 
				(*wptr)(ix,iy,iz) = 0;
}

int nn4_ctfReconstructor::insert_slice(EMData* slice, const Transform3D& t) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	if ((slice->get_xsize() != slice->get_ysize()) 
		|| slice->get_xsize() != vnx) {
		// FIXME: Why doesn't this throw an exception?
		LOGERR("Tried to insert a slice that is the wrong size.");
		return 1;
	}
	// process 2-d slice -- subtract the average outside of the circle, zero-pad, fft extend, and fft
	EMData* temp = slice->average_circ_sub();
	// Need to use zeropad_ntimes instead of pad_fft here for zero padding
	// because only the former centers the original image in the 
	// larger area.  FIXME!
	EMData* zeropadded = temp->zeropad_ntimes(npad);
	if( temp ) {delete temp; temp = 0;}
	EMData* padfftslice = zeropadded->do_fft();
	if( zeropadded ) { delete zeropadded; zeropadded = 0;}
	padfftslice->center_origin_fft();
	// insert slice for all symmetry related positions
	float defocus = slice->get_attr("defocus");
	for (int isym=0; isym < nsym; isym++) {
		Transform3D tsym = t.get_sym(symmetry, isym);
		v->nn_ctf(*wptr, padfftslice, tsym, defocus);
	}
	if( padfftslice ) { delete padfftslice; padfftslice = 0;}
	return 0;
}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
EMData* nn4_ctfReconstructor::finish() {
	float snr = params["snr"];
	float osnr = 1.0f/snr;
	EMArray<float>& w = *wptr;
	v->symplane0_ctf(w);
	// normalize
	for (int iz = 1; iz <= vnzp; iz++) {
		for (int iy = 1; iy <= vnyp; iy++) {
			for (int ix = 0; ix <= vnxc; ix++) {//std::cout<<"  "<<iz<<"  "<<iy<<"  "<<ix<<"  "<<nr(ix,iy,iz)<<"  "<<(-2*((ix+iy+iz)%2)+1)<<"  "<<(*v)(ix,iy,iz)<<std::endl;
				if (w(ix,iy,iz) > 0.f) {//(*v) should be treated as complex!!
					float  tmp = (-2*((ix+iy+iz)%2)+1)/(w(ix,iy,iz)+osnr);
					(*v)(2*ix,iy,iz) *= tmp;
					(*v)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}
	//  do not need weights anymore
	if( wptr ) { delete wptr; wptr = 0;}
	// back fft
	v->do_ift_inplace();
	EMData* win = v->window_center(vnx);
	if( v ) { delete v; v = 0;}
	float *tw = win->get_data();
	//  mask and subtract circumference average
	int ix = win->get_xsize(),iy = win->get_ysize(),iz = win->get_zsize();
	int L2 = (ix/2)*(ix/2);
	int L2P = (ix/2-1)*(ix/2-1);
	int IP = ix/2+1;
	float  TNR = 0.0f;
	int m = 0;
	for (int k = 1; k <= iz; k++) {
		for (int j = 1; j <= iy; j++) {
			for (int i = 1; i <= ix; i++) {
				int LR = (k-IP)*(k-IP)+(j-IP)*(j-IP)+(i-IP)*(i-IP);
				if (LR<=L2) {
					if(LR >= L2P && LR <= L2) {
						TNR += tw(i,j,k);
						m++;
					}
				}
			}
		}
	}
	TNR /=float(m);
	for (int k = 1; k <= iz; k++) {
		for (int j = 1; j <= iy; j++) {
			for (int i = 1; i <= ix; i++) {
				int LR = (k-IP)*(k-IP)+(j-IP)*(j-IP)+(i-IP)*(i-IP);
				if (LR<=L2) tw(i,j,k) -= TNR; else tw(i,j,k) = 0.0f;
			}
		}
	}
	return win;
	// clean up
}
#undef  tw




void EMAN::dump_reconstructors()
{
	dump_factory < Reconstructor > ();
}

map<string, vector<string> > EMAN::dump_reconstructors_list()
{
	return dump_factory_list < Reconstructor > ();
}

/* vim: set ts=4 noet: */
