/**
 * $Id$
 */
#include "emdata.h"
#include "ctf.h"

using namespace EMAN;
 
EMData* EMData::get_fft_amplitude2D()
{
	ENTERFUNC;

//	int ndim = get_ndim();
	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}
	if (nz>1) {
		LOGERR("2D image expected. Input image is 3D");
		throw ImageFormatException("2D odd square complex image"
			" expected Input image is 3D.");
	}


	int nx2 = nx/2;

	EMData *dat = copy_head();

	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float temp=0;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx2; i++) {
			temp = (*this)(2*i,j)*(*this)(2*i,j);
			temp += (*this)(2*i+1,j)*(*this)(2*i+1,j);
			(*dat)(i,j) = std::sqrt(temp);
		}
	}


	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


EMData* EMData::get_fft_amplitude()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();

	int ndim = get_ndim();

	if (ndim == 3) {
		for (int k = 1; k < nz; k++) {
			for (int j = 1; j < ny; j++) {
				for (int i = 0; i < nx2/2; i++) {
					d[k*nx2*ny+j*nx2+nx2/2+i] = rdata[k*nx*ny+j*nx+2*i];
					d[(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i] = rdata[k*nx*ny+j*nx+2*i];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx2/2; i++) {
				d[j*nx2+nx2/2+i] = rdata[j*nx+2*i];
				d[(ny-j)*nx2+nx2/2-i] = rdata[j*nx+2*i];
			}
		}
	}

	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


EMData* EMData::get_fft_phase()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();

	int ndim = get_ndim();
	if (ndim == 3) {
		for (int k = 1; k < nz; k++) {
			for (int j = 1; j < ny; j++) {
				for (int i = 0; i < nx2/2; i++) {
					d[k*nx2*ny+j*nx2+nx2/2+i] = rdata[k*nx*ny+j*nx+2*i+1];
					d[(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i] = -rdata[k*nx*ny+j*nx+2*i+1];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx2/2; i++) {
				d[j*nx2+nx2/2+i] = rdata[j*nx+2*i+1];
				d[(ny-j)*nx2+nx2/2-i] = -rdata[j*nx+2*i+1];
			}
		}
	}
	done_data();

	dat->done_data();
	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


float *EMData::get_data() const
{
	//flags |= EMDATA_BUSY;
	return rdata;
}


void EMData::done_data()
{
	update();
}



float EMData::calc_center_density()
{
	ENTERFUNC;

	float center = get_attr("mean");
	float sigma = get_attr("sigma");
	float ds = sigma / 2;
	size_t size = nx * ny * nz;
	float *d = get_data();
	float sigma1 = sigma / 20;
	float sigma2 = sigma / 1000;

	while (ds > sigma1) {
		double sum = 0;
		int norm = 0;

		for (size_t i = 0; i < size; i++) {
			if (fabs(d[i] - center) < ds) {
				sum += d[i];
				norm++;
			}
		}
		if (!norm) {
			break;
		}
		float mean = (float)(sum / norm);
		if (fabs(mean - center) < sigma2) {
			ds *= 0.5f;
		}
		center = mean;
	}
	EXITFUNC;

	return center;
}


float EMData::calc_sigma_diff()
{
	ENTERFUNC;

	float *d = get_data();
	float mean = get_attr("mean");
	float sigma = get_attr("sigma");

	double sum_up = 0;
	double sum_down = 0;
	int nup = 0;
	int ndown = 0;

	size_t size = nx * ny * nz;

	for (size_t i = 0; i < size; i++) {
		if (d[i] > mean) {
			sum_up += Util::square(d[i] - mean);
			nup++;
		}
		else {
			sum_down += Util::square(mean - d[i]);
			ndown++;
		}
	}

	float sigup = std::sqrt((float)sum_up / nup);
	float sigdown = std::sqrt((float)sum_down / ndown);
	float sig_diff = fabs(sigup - sigdown) / sigma;


	EXITFUNC;
	return sig_diff;

}

	
IntPoint EMData::calc_min_location() const
{
	ENTERFUNC;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float min = FLT_MAX;
	int min_x = 0;
	int min_y = 0;
	int min_z = 0;
	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = rdata[l + cur_y];
				if (t < min) {
					min_x = l;
					min_y = k;
					min_z = j;
					min = t;
				}
			}
		}
	}

	return IntPoint(min_x, min_y, min_z);
}


IntPoint EMData::calc_max_location() const
{
	ENTERFUNC;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float max = -FLT_MAX;
	int max_x = 0;
	int max_y = 0;
	int max_z = 0;
	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = rdata[l + cur_y];
				if (t > max) {
					max_x = l;
					max_y = k;
					max_z = j;
					max = t;
				}
			}
		}
	}

	EXITFUNC;
	return IntPoint(max_x, max_y, max_z);
}


int EMData::calc_min_index() const
{
	IntPoint min_location = calc_min_location();
	int i = min_location[0] + min_location[1] * nx + min_location[2] * nx * ny;
	return i;
}


int EMData::calc_max_index() const
{
	IntPoint max_location = calc_max_location();
	int i = max_location[0] + max_location[1] * nx + max_location[2] * nx * ny;
	return i;
}


vector<Pixel> EMData::calc_highest_locations(float threshold)
{
	ENTERFUNC;

	vector<Pixel> result;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	int nxy = nx * ny;

	for (int j = 0; j < nz; j++) {
		int cur_z = j * nxy;

		for (int k = 0; k < ny; k++) {
			int cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float v = rdata[l + cur_y];
				if (v > threshold) {
					result.push_back(Pixel(l, k, j, v));
				}
			}
		}
	}

	std::sort(result.begin(), result.end());

	EXITFUNC;
	return result;
}


float EMData::get_edge_mean() const
{
	ENTERFUNC;

	int di = 0;
	double edge_sum = 0;
	float edge_mean = 0;
	int nxy = nx * ny;

	if (nz == 1) {
		for (int i = 0, j = (ny - 1) * nx; i < nx; i++, j++) {
			edge_sum += rdata[i] + rdata[j];
		}
		for (int i = 0, j = nx - 1; i < nxy; i += nx, j += nx) {
			edge_sum += rdata[i] + rdata[j];
		}
		edge_mean = (float)edge_sum / (nx * 2 + ny * 2);
	}
	else {
		if (nx == ny && nx == nz * 2 - 1) {
			for (int j = (nxy * (nz - 1)); j < nxy * nz; j++, di++) {
				edge_sum += rdata[j];
			}
		}
		else {
			for (int i = 0, j = (nxy * (nz - 1)); i < nxy; i++, j++, di++) {
				edge_sum += rdata[i] + rdata[j];
			}
		}

		int nxy2 = nx * (ny - 1);
		for (int k = 1; k < nz - 1; k++) {
			int k2 = k * nxy;
			int k3 = k2 + nxy2;
			for (int i = 0; i < nx; i++, di++) {
				edge_sum += rdata[i + k2] + rdata[i + k3];
			}
		}
		for (int k = 1; k < nz - 1; k++) {
			int k2 = k * nxy;
			int k3 = nx - 1 + k2;
			for (int i = 1; i < ny - 1; i++, di++) {
				edge_sum += rdata[i * nx + k2] + rdata[i * nx + k3];
			}
		}

		edge_mean = (float)edge_sum / (di * 2);
	}
	EXITFUNC;

	return edge_mean;
}


float EMData::get_circle_mean()
{
	ENTERFUNC;

	static bool busy = false;
	static EMData *mask = 0;

	while (busy);
	busy = true;

	if (!mask || !EMUtil::is_same_size(this, mask)) {
		if (!mask) {
			mask = new EMData();
		}
		mask->set_size(nx, ny, nz);
		mask->to_one();

		float radius = (float)(ny / 2 - 2);
		mask->process_inplace("eman1.mask.sharp", Dict("inner_radius", radius - 1,
									   "outer_radius", radius + 1));

		int n = 0;
		float *d = mask->get_data();

		for (int i = 0; i < nx * ny * nz; i++) {
			if (d[i]) {
				n++;
			}
		}
		mask->mult(1.0f / n);
	}

	float result = dot(mask);
	busy = false;

	EXITFUNC;
	return result;
}


void EMData::set_ctf(Ctf * new_ctf)
{
	ENTERFUNC;

	vector<float> vctf = new_ctf->to_vector();
	attr_dict["ctf"] = vctf;

	EXITFUNC;
}

Ctf * EMData::get_ctf() const
{
	ENTERFUNC;
	
	SimpleCtf * ctf = new SimpleCtf();
	ctf->from_vector(attr_dict["ctf"]);
	
	EXITFUNC;
	return ctf;
}

void EMData::set_size(int x, int y, int z)
{
	ENTERFUNC;

	if (x <= 0) {
		throw InvalidValueException(x, "x size <= 0");
	}
	else if (y <= 0) {
		throw InvalidValueException(y, "y size <= 0");
	}
	else if (z <= 0) {
		throw InvalidValueException(z, "z size <= 0");
	}

	int old_nx = nx;
	nx = x;
	ny = y;
	nz = z;

	size_t size = (size_t)x * (size_t)y * (size_t)z * sizeof(float);
	rdata = static_cast < float *>(realloc(rdata, size));
	update();

	attr_dict["nx"] = x;
	attr_dict["ny"] = y;
	attr_dict["nz"] = z;

	if (old_nx == 0) {
		memset(rdata, 0, size);
	}

	if (supp) {
		free(supp);
		supp = 0;
	}
	EXITFUNC;
}


MArray2D EMData::get_2dview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(rdata, dims, boost::fortran_storage_order());
	return marray;
}


MArray3D EMData::get_3dview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(rdata, dims, boost::fortran_storage_order());
	return marray;
}


MCArray2D EMData::get_2dcview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(rdata);
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}


MCArray3D EMData::get_3dcview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(rdata);
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}


MCArray3D* EMData::get_3dcviewptr() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(rdata);
	MCArray3D* marray = new MCArray3D(cdata, dims,
									  boost::fortran_storage_order());
	return marray;
}


MArray2D EMData::get_2dview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(rdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}


MArray3D EMData::get_3dview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(rdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}


MCArray2D EMData::get_2dcview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(rdata);
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}


MCArray3D EMData::get_3dcview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(rdata);
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}


EMObject EMData::get_attr(const string & key)
{
	ENTERFUNC;

	update_stat();

	float mean = attr_dict["mean"];
	float sigma = attr_dict["sigma"];
	size_t size = nx * ny * nz;

	if (key == "kurtosis") {
		double kurtosis_sum = 0;

		for (size_t k = 0; k < size; k++) {
			float t = (rdata[k] - mean) / sigma;
			float tt = t * t;
			kurtosis_sum += tt * tt;
		}

		float kurtosis = (float)(kurtosis_sum / size - 3.0);
		attr_dict["kurtosis"] = kurtosis;
		return attr_dict["kurtosis"];
	}
	else if (key == "skewness") {
		double skewness_sum = 0;
		for (size_t k = 0; k < size; k++) {
			float t = (rdata[k] - mean) / sigma;
			skewness_sum +=  t * t * t;
		}
		float skewness = (float)(skewness_sum / size);
		attr_dict["skewness"] = skewness;
		return attr_dict["skewness"];
	}
	else if (key == "changecount") return EMObject(changecount);


	EXITFUNC;
	if(attr_dict.has_key(key)) {
		return attr_dict[key];
	}
	else {
		throw NotExistingObjectException(key, "Not exist");
	}
}


Dict EMData::get_attr_dict() const
{
	update_stat();
	return Dict(attr_dict);
}

void EMData::set_attr_dict(const Dict & new_dict)
{
	vector<string> new_keys = new_dict.keys();
	vector<string>::const_iterator it;
	for(it = new_keys.begin(); it!=new_keys.end(); ++it) {
			this->set_attr(*it, new_dict[*it]);
	} 
}

void EMData::del_attr(const string & attr_name)
{
	attr_dict.erase(attr_name);
}

void EMData::del_attr_dict(const vector<string> & del_keys)
{
	vector<string>::const_iterator it;
	for(it=del_keys.begin(); it!=del_keys.end(); ++it) {
		this->del_attr(*it);
	}
}

vector<float> EMData::get_data_pickle() const
{
	vector<float> vf;
	vf.resize(nx*ny*nz);
	std::copy(rdata, rdata+nx*ny*nz, vf.begin());
	
	return vf;
}

void EMData::set_data_pickle(const vector<float>& vf)
{
	rdata = (float *)malloc(nx*ny*nz*sizeof(float));
	std::copy(vf.begin(), vf.end(), rdata);
}

int EMData::get_supp_pickle() const
{
	return 0;
}

void EMData::set_supp_pickle(int)
{
	this->supp = 0;
}
