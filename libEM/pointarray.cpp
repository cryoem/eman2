/**
 * $Id$
 */
#include "pointarray.h"
#include <gsl/gsl_multimin.h>

#if defined OPTPP
#include "OptQNewton.h"
#include "NLF.h"
#endif

using namespace EMAN;

PointArray::PointArray()
{
	points = 0;
	n = 0;
}

PointArray::PointArray(unsigned int nn)
{
	n = nn;
	points = (double *) malloc(4 * n * sizeof(double));
}

PointArray::~PointArray()
{
	free(points);
}

unsigned int PointArray::get_number_points()
{
	return n;
}

void PointArray::set_number_points(unsigned int nn)
{
	if (n != nn) {
		n = nn;
		points = (double *) realloc(points, 4 * n * sizeof(double));
	}
}

double *PointArray::get_points_array()
{
	return points;
}

void PointArray::set_points_array(double *p)
{
	points = p;
}

bool PointArray::read_from_pdb(const char *file)
{
	struct stat filestat;
	stat(file, &filestat);
	set_number_points(filestat.st_size / 80);

	FILE *fp = fopen(file, "r");
	char s[200];
	unsigned int count = 0;
	while ((fgets(s, 200, fp) != NULL)) {
		if (strncmp(s, "ENDMDL", 6) == 0)
			break;
		if (strncmp(s, "ATOM", 4) != 0)
			continue;

		if (s[13] == ' ')
			s[13] = s[14];
		if (s[13] == ' ')
			s[13] = s[15];

		float e = 0;
		char ctt, ctt2 = ' ';
		if (s[13] == ' ')
			ctt = s[14];
		else if (s[12] == ' ') {
			ctt = s[13];
			ctt2 = s[14];
		}
		else {
			ctt = s[12];
			ctt2 = s[13];
		}

		switch (ctt) {
		case 'H':
			e = 1.0;
			break;
		case 'C':
			e = 6.0;
			break;
		case 'A':
			if (ctt2 == 'U') {
				e = 79.0;
				break;
			}
			// treat 'A'mbiguous atoms as N, not perfect, but good enough
		case 'N':
			e = 7.0;
			break;
		case 'O':
			e = 8.0;
			break;
		case 'P':
			e = 15.0;
			break;
		case 'S':
			e = 16.0;
			break;
		case 'W':
			e = 18.0;
			break;				// ficticious water 'atom'
		default:
			fprintf(stderr, "Unknown atom %c%c\n", ctt, ctt2);
			e = 0;
		}
		if (e == 0)
			continue;

		float x, y, z;
		sscanf(&s[28], " %f %f %f", &x, &y, &z);

		if (count + 1 > get_number_points())
			set_number_points(2 * (count + 1));
		points[4 * count] = x;
		points[4 * count + 1] = y;
		points[4 * count + 2] = z;
		points[4 * count + 3] = e;
		count++;
	}
	fclose(fp);
	set_number_points(count);
	return true;
}

FloatPoint PointArray::get_center()
{
	double xc, yc, zc;
	xc = yc = zc = 0.0;
	double norm = 0.0;
	for (unsigned int i = 0; i < 4 * get_number_points(); i += 4) {
		xc += points[i] * points[i + 3];
		yc += points[i + 1] * points[i + 3];
		zc += points[i + 2] * points[i + 3];
		norm += points[i + 3];
	}
	if (norm <= 0) {
		fprintf(stderr, "Abnormal total value (%g) for PointArray, it should be >0\n", norm);
		return FloatPoint(0, 0, 0);
	}
	else
		return FloatPoint(xc / norm, yc / norm, zc / norm);
}

void PointArray::center_to_zero()
{
	FloatPoint center = get_center();
	for (unsigned int i = 0; i < 4 * get_number_points(); i += 4) {
		points[i] -= center.x;
		points[i + 1] -= center.y;
		points[i + 2] -= center.z;
	}
}

Region PointArray::get_bounding_box()
{
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	xmin = xmax = points[0];
	ymin = ymax = points[1];
	zmin = zmax = points[2];
	for (unsigned int i = 0; i < 4 * get_number_points(); i += 4) {
		if (points[i] > xmax)
			xmax = points[i];
		if (points[i] < xmin)
			xmin = points[i];
		if (points[i + 1] > ymax)
			ymax = points[i + 1];
		if (points[i + 1] < ymin)
			ymin = points[i + 1];
		if (points[i + 2] > zmax)
			zmax = points[i + 2];
		if (points[i + 2] < zmin)
			zmin = points[i + 2];
	}
	return Region(xmin, ymin, zmin, xmax - xmin, ymax - ymin, zmax - zmin);
}

void PointArray::set_from(PointArray * source, Transform * transform, string sym)
{
	set_from(source->get_points_array(), source->get_number_points(), transform, sym);

}

void PointArray::set_from(double *src, unsigned int num, Transform * transform, string sym)
{
	Rotation rot = transform->get_rotation();
	rot.set_sym(sym);
	unsigned int nsym = rot.get_max_nsym();

	if (get_number_points() != nsym * num)
		set_number_points(nsym * num);

	double *target = get_points_array();

	Vec3 < float >pre_trans = transform->get_pre_translate();
	double tx0 = pre_trans[0], ty0 = pre_trans[1], tz0 = pre_trans[2];

	Vec3 < float >post_trans = transform->get_post_translate();
	double tx1 = post_trans[0], ty1 = post_trans[1], tz1 = post_trans[2];

	for (unsigned int s = 0; s < nsym; s++) {
		Rotation rot2 = rot.get_sym(s);
		Matrix3f m = rot2.get_matrix3();
		double m00 = m[0][0], m01 = m[0][1], m02 = m[0][2];
		double m10 = m[1][0], m11 = m[1][1], m12 = m[1][2];
		double m20 = m[2][0], m21 = m[2][1], m22 = m[2][2];
		int index = s * 4 * num;
		for (unsigned int i = 0; i < 4 * num; i += 4, index += 4) {
			double x = src[i] + tx0, y = src[i + 1] + ty0, z = src[i + 2] + tz0;
			target[index] = x * m00 + y * m01 + z * m02 + tx1;
			target[index + 1] = x * m10 + y * m11 + z * m12 + ty1;
			target[index + 2] = x * m20 + y * m21 + z * m22 + tz1;
			target[index + 3] = src[i + 3];
			//printf("s=%d\tx,y,z=%g,%g,%g => %g,%g,%g\n",s,src[i],src[i+1],src[i+2],target[index    ],target[index + 1],target[index + 2]);
		}
	}
}

int cmp_axis_x(const void *a, const void *b)
{
	double diff = ((double *) a)[0] - ((double *) b)[0];
	if (diff < 0.0)
		return -1;
	else if (diff > 0.0)
		return 1;
	else
		return 0;
}
int cmp_axis_y(const void *a, const void *b)
{
	double diff = ((double *) a)[1] - ((double *) b)[1];
	if (diff < 0.0)
		return -1;
	else if (diff > 0.0)
		return 1;
	else
		return 0;
}
int cmp_axis_z(const void *a, const void *b)
{
	double diff = ((double *) a)[2] - ((double *) b)[2];
	if (diff < 0.0)
		return -1;
	else if (diff > 0.0)
		return 1;
	else
		return 0;
}
int cmp_val(const void *a, const void *b)
{
	double diff = ((double *) a)[3] - ((double *) b)[3];
	if (diff < 0.0)
		return -1;
	else if (diff > 0.0)
		return 1;
	else
		return 0;
}
void PointArray::sort_by_axis(int axis)
{
	if (axis == 0)
		qsort(points, n, sizeof(double) * 4, cmp_axis_x);
	else if (axis == 1)
		qsort(points, n, sizeof(double) * 4, cmp_axis_y);
	else if (axis == 2)
		qsort(points, n, sizeof(double) * 4, cmp_axis_z);
	else
		qsort(points, n, sizeof(double) * 4, cmp_val);
}


EMData *PointArray::pdb2mrc_by_summation(int map_size, float apix, float res)
{
	double gauss_real_width = res / (M_PI);	// in Angstrom, res is in Angstrom
	//if ( gauss_real_width < apix) LOGERR("PointArray::projection_by_summation(): apix(%g) is too large for resolution (%g Angstrom in Fourier space) with %g pixels of 1/e half width", apix, res, gauss_real_width);

	double min_table_val = 1e-7;
	double max_table_x = sqrt(-log(min_table_val));	// for exp(-x*x)

	double table_step_size = 0.001;	// number of steps for each pixel
	double inv_table_step_size = 1.0 / table_step_size;
	int table_size = int (max_table_x * gauss_real_width / (apix * table_step_size) * 1.25);
	double *table = (double *) malloc(sizeof(double) * table_size);
	for (int i = 0; i < table_size; i++) {
		double x = -i * table_step_size * apix / gauss_real_width;
		table[i] = exp(-x * x);
	}

	int gbox = int (max_table_x * gauss_real_width / apix);	// local box half size in pixels to consider for each point
	if (gbox <= 0)
		gbox = 1;

	sort_by_axis(2);			// sort by Z-axis

	EMData *map = new EMData();
	map->set_size(map_size, map_size, map_size);
	map->to_zero();
	float *pd = map->get_data();
	for (unsigned int s = 0; s < get_number_points(); s++) {
		double xc = points[4 * s] / apix + map_size / 2;
		double yc = points[4 * s + 1] / apix + map_size / 2;
		double zc = points[4 * s + 2] / apix + map_size / 2;
		double fval = points[4 * s + 3];
		int imin = int (xc) - gbox, imax = int (xc) + gbox;
		int jmin = int (yc) - gbox, jmax = int (yc) + gbox;
		int kmin = int (zc) - gbox, kmax = int (zc) + gbox;
		if (imin < 0)
			imin = 0;
		if (jmin < 0)
			jmin = 0;
		if (kmin < 0)
			kmin = 0;
		if (imax > map_size)
			imax = map_size;
		if (jmax > map_size)
			jmax = map_size;
		if (kmax > map_size)
			kmax = map_size;

		for (int k = kmin; k < kmax; k++) {
			int table_index_z = int (fabs(k - zc) * inv_table_step_size);
			double zval = table[table_index_z];
			int pd_index_z = k * map_size * map_size;
			for (int j = jmin; j < jmax; j++) {
				int table_index_y = int (fabs(j - yc) * inv_table_step_size);
				double yval = table[table_index_y];
				int pd_index = pd_index_z + j * map_size + imin;
				for (int i = imin; i < imax; i++, pd_index++) {
					int table_index_x = int (fabs(i - xc) * inv_table_step_size);
					double xval = table[table_index_x];
					pd[pd_index] += fval * zval * yval * xval;
				}
			}
		}
	}
	//for(int i=0; i<map_size*map_size; i++) pd[i]/=sqrt(M_PI);
	map->done_data();
	return map;
/*
	// precalculate a prototypical Gaussian to resample
	// 64^3 box with a real-space 1/2 width of 12 pixels
	EMData *gaus = new EMData();
	gaus->set_size(64, 64, 64);
	gaus->to_one();
	gaus->filter("GaussMask", Dict("outer_radius", 12.0));

	EMData *map = new EMData();
	map->set_size(map_size, map_size, map_size);
	map->to_zero();

	for (unsigned int i = 0; i < 4 * get_number_points(); i += 4) {
		map->insert_scaled_sum(gaus,
							   FloatPoint(points[i] / apix + map_size / 2,
										  points[i + 1] / apix + map_size / 2,
										  points[i + 2] / apix + map_size / 2),
							   res / (PI * 12.0 * apix), points[i+3]);
	}
	return map;
*/
}


EMData *PointArray::projection_by_summation(int image_size, float apix, float res)
{
	double gauss_real_width = res / (M_PI);	// in Angstrom, res is in Angstrom
	//if ( gauss_real_width < apix) LOGERR("PointArray::projection_by_summation(): apix(%g) is too large for resolution (%g Angstrom in Fourier space) with %g pixels of 1/e half width", apix, res, gauss_real_width);

	double min_table_val = 1e-7;
	double max_table_x = sqrt(-log(min_table_val));	// for exp(-x*x)

	//double table_step_size = 0.001;    // number of steps for x=[0,1] in exp(-x*x) 
	//int table_size = int(max_table_x / table_step_size *1.25);
	//double* table = (double*)malloc(sizeof(double) * table_size);
	//for(int i=0; i<table_size; i++) table[i]=exp(-i*i*table_step_size*table_step_size);

	double table_step_size = 0.001;	// number of steps for each pixel
	double inv_table_step_size = 1.0 / table_step_size;
	int table_size = int (max_table_x * gauss_real_width / (apix * table_step_size) * 1.25);
	double *table = (double *) malloc(sizeof(double) * table_size);
	for (int i = 0; i < table_size; i++) {
		double x = -i * table_step_size * apix / gauss_real_width;
		table[i] = exp(-x * x);
	}

	int gbox = int (max_table_x * gauss_real_width / apix);	// local box half size in pixels to consider for each point
	if (gbox <= 0)
		gbox = 1;
	EMData *proj = new EMData();
	proj->set_size(image_size, image_size, 1);
	proj->to_zero();
	float *pd = proj->get_data();
	for (unsigned int s = 0; s < get_number_points(); s++) {
		double xc = points[4 * s] / apix + image_size / 2;
		double yc = points[4 * s + 1] / apix + image_size / 2;
		double fval = points[4 * s + 3];
		int imin = int (xc) - gbox, imax = int (xc) + gbox;
		int jmin = int (yc) - gbox, jmax = int (yc) + gbox;
		if (imin < 0)
			imin = 0;
		if (jmin < 0)
			jmin = 0;
		if (imax > image_size)
			imax = image_size;
		if (jmax > image_size)
			jmax = image_size;

		for (int j = jmin; j < jmax; j++) {
			//int table_index_y = int(fabs(j-yc)*apix/gauss_real_width/table_step_size);
			int table_index_y = int (fabs(j - yc) * inv_table_step_size);
			double yval = table[table_index_y];
#ifdef DEBUG
			//double yval2 = exp( - (j-yc)*(j-yc)*apix*apix/(gauss_real_width*gauss_real_width));
			//if(fabs(yval2-yval)/yval2>1e-2) printf("s=%d\txc,yc=%g,%g\tyval,yval2=%g,%g\tdiff=%g\n",s,xc,yc,yval,yval2,fabs(yval2-yval)/yval2);
#endif
			int pd_index = j * image_size + imin;
			for (int i = imin; i < imax; i++, pd_index++) {
				//int table_index_x = int(fabs(i-xc)*apix/gauss_real_width/table_step_size);
				int table_index_x = int (fabs(i - xc) * inv_table_step_size);
				double xval = table[table_index_x];
#ifdef DEBUG
				//double xval2 = exp( - (i-xc)*(i-xc)*apix*apix/(gauss_real_width*gauss_real_width));
				//if(fabs(xval2-xval)/xval2>1e-2) printf("\ts=%d\txc,yc=%g,%g\txval,xval2=%g,%g\tdiff=%g\n",s,xc,yc,xval,xval2,fabs(xval2-xval)/xval2);
#endif
				pd[pd_index] += fval * yval * xval;
			}
		}
	}
	for (int i = 0; i < image_size * image_size; i++)
		pd[i] /= sqrt(M_PI);
	proj->done_data();
	return proj;
/*
	// precalculate a prototypical Gaussian to resample
	// 64^2 box with a real-space 1/2 width of 12 pixels
	EMData *gaus = new EMData();
	gaus->set_size(64, 64, 1);
	gaus->to_one();
	gaus->filter("GaussMask", Dict("outer_radius", 12.0));

	EMData *proj = new EMData();
	proj->set_size(image_size, image_size, 1);
	proj->to_zero();

	for (unsigned int i = 0; i < 4 * get_number_points(); i += 4) {
		proj->insert_scaled_sum(gaus,
								FloatPoint(points[i] / apix + image_size / 2,
										   points[i + 1] / apix + image_size / 2),
								res / (PI * 12.0 * apix), points[i+3]);
	}
	return proj;
*/
}


EMData *PointArray::pdb2mrc_by_nfft(int map_size, float apix, float res)
{
#if defined NFFT
	nfft_3D_plan my_plan;		// plan for the nfft

	/** init an 3 dimensional plan */
	nfft_3D_init(&my_plan, map_size, get_number_points());

	/** init atom positions to the non-uniform nodes */
	for (int j = 0, i = 0; j < my_plan.M; j++, i += 4) {
		// FFTW and nfft use row major array layout, EMAN uses column major
		my_plan.v[3 * j + 2] = (fftw_real) (points[i] / (apix * map_size));
		my_plan.v[3 * j + 1] = (fftw_real) (points[i + 1] / (apix * map_size));
		my_plan.v[3 * j] = (fftw_real) (points[i + 2] / (apix * map_size));
		my_plan.f[j].re = (fftw_real) (points[i + 3]);
		my_plan.f[j].im = 0.0;
	}

	/** precompute psi, the entries of the matrix B */
	if (my_plan.nfft_flags & PRE_PSI) {
		nfft_3D_precompute_psi(&my_plan);
	}

	// compute the uniform Fourier transform
	nfft_3D_transpose(&my_plan);

	// copy the Fourier transform to EMData data array
	EMData *fft = new EMData();
	fft->set_size(map_size + 2, map_size, map_size);
	fft->set_complex(true);
	fft->set_ri(true);
	fft->to_zero();
	float *data = fft->get_data();
	double norm = 1.0 / (map_size * map_size * map_size);
	for (int k = 0; k < map_size; k++) {
		for (int j = 0; j < map_size; j++) {
			for (int i = 0; i < map_size / 2; i++) {
				data[k * map_size * (map_size + 2) + j * (map_size + 2) + 2 * i] =
					(float) (my_plan.
							 f_hat[k * map_size * map_size + j * map_size + i +
								   map_size / 2].re) * norm;
				data[k * map_size * (map_size + 2) + j * (map_size + 2) + 2 * i + 1] =
					(float) (my_plan.
							 f_hat[k * map_size * map_size + j * map_size + i +
								   map_size / 2].im) * norm * (-1.0);
			}
		}
	}
	/** finalise the nfft plan */
	nfft_3D_finalize(&my_plan);

	// low pass filter
	double sigma2 = (map_size * apix / res) * (map_size * apix / res);
	int index = 0;
	for (int k = 0; k < map_size; k++) {
		double RZ2 = (k - map_size / 2) * (k - map_size / 2);
		for (int j = 0; j < map_size; j++) {
			double RY2 = (j - map_size / 2) * (j - map_size / 2);
			for (int i = 0; i < map_size / 2 + 1; i++, index += 2) {
				float val = exp(-(i * i + RY2 + RZ2) / sigma2);
				data[index] *= val;
				data[index + 1] *= val;
			}
		}
	}
	fft->done_data();
	//fft->filter("GaussLowpass",Dict("lowpass", map_size*apix/res));

	fft->filter("Phase180");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	delete fft;
	return map;
#elif defined NFFT2
	nfft_plan my_plan;			// plan for the nfft

	/** init an 3 dimensional plan */
	nfft_init_3d(&my_plan, map_size, map_size, map_size, get_number_points());

	/** init atom positions to the non-uniform nodes */
	for (int j = 0, i = 0; j < my_plan.M; j++, i += 4) {
		// FFTW and nfft use row major array layout, EMAN uses column major
		my_plan.x[3 * j + 2] = (double) (points[i] / (apix * map_size));
		my_plan.x[3 * j + 1] = (double) (points[i + 1] / (apix * map_size));
		my_plan.x[3 * j] = (double) (points[i + 2] / (apix * map_size));
		my_plan.f[j][0] = (double) (points[i + 3]);
		my_plan.f[j][1] = 0.0;
	}

	/** precompute psi, the entries of the matrix B */
	if (my_plan.nfft_flags & PRE_PSI) {
		nfft_precompute_psi(&my_plan);
		if (my_plan.nfft_flags & PRE_FULL_PSI)
			nfft_full_psi(&my_plan, pow(10, -10));
	}

	// compute the uniform Fourier transform
	nfft_transposed(&my_plan);

	// copy the Fourier transform to EMData data array
	EMData *fft = new EMData();
	fft->set_size(map_size + 2, map_size, map_size);
	fft->set_complex(true);
	fft->set_ri(true);
	fft->to_zero();
	float *data = fft->get_data();
	double norm = 1.0 / (map_size * map_size * map_size);
	for (int k = 0; k < map_size; k++) {
		for (int j = 0; j < map_size; j++) {
			for (int i = 0; i < map_size / 2; i++) {
				data[k * map_size * (map_size + 2) + j * (map_size + 2) + 2 * i] =
					(float) (my_plan.
							 f_hat[k * map_size * map_size + j * map_size + i +
								   map_size / 2][0]) * norm;
				data[k * map_size * (map_size + 2) + j * (map_size + 2) + 2 * i + 1] =
					(float) (my_plan.
							 f_hat[k * map_size * map_size + j * map_size + i +
								   map_size / 2][1]) * norm;
			}
		}
	}
	/** finalise the nfft plan */
	nfft_finalize(&my_plan);

	// low pass filter
	double sigma2 = (map_size * apix / res) * (map_size * apix / res);
	int index = 0;
	for (int k = 0; k < map_size; k++) {
		double RZ2 = (k - map_size / 2) * (k - map_size / 2);
		for (int j = 0; j < map_size; j++) {
			double RY2 = (j - map_size / 2) * (j - map_size / 2);
			for (int i = 0; i < map_size / 2 + 1; i++, index += 2) {
				float val = exp(-(i * i + RY2 + RZ2) / sigma2);
				data[index] *= val;
				data[index + 1] *= val;
			}
		}
	}
	fft->done_data();
	//fft->filter("GaussLowpass",Dict("lowpass", map_size*apix/res));

	fft->filter("Phase180");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	delete fft;
	return map;
#else
	return 0;
#endif
}

EMData *PointArray::projection_by_nfft(int image_size, float apix, float res)
{
#if defined NFFT
	nfft_2D_plan my_plan;		// plan for the nfft 
	int N[2], n[2];
	N[0] = image_size;
	n[0] = next_power_of_2(2 * image_size);
	N[1] = image_size;
	n[1] = next_power_of_2(2 * image_size);

	/** init an 2 dimensional plan */
	nfft_2D_init(&my_plan, image_size, get_number_points());
	//nfft_2D_init_specific(&my_plan, N, get_number_points(), n, 3,
	//                 PRE_PHI_HUT | PRE_PSI,
	//                 FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
	/** init atom positions to the non-uniform nodes */
	for (int j = 0, i = 0; j < my_plan.M; j++, i += 4) {
		// FFTW and nfft use row major array layout, EMAN uses column major
		my_plan.v[2 * j + 1] = (fftw_real) (points[i] / (apix * image_size));
		my_plan.v[2 * j] = (fftw_real) (points[i + 1] / (apix * image_size));
		my_plan.f[j].re = (fftw_real) points[i + 3];
		my_plan.f[j].im = 0.0;
	}

	/** precompute psi, the entries of the matrix B */
	if (my_plan.nfft_flags & PRE_PSI) {
		nfft_2D_precompute_psi(&my_plan);
	}

	// compute the uniform Fourier transform
	nfft_2D_transpose(&my_plan);

	// copy the Fourier transform to EMData data array
	EMData *fft = new EMData();
	fft->set_size(image_size + 2, image_size, 1);
	fft->set_complex(true);
	fft->set_ri(true);
	fft->to_zero();
	float *data = fft->get_data();
	double norm = 1.0 / (image_size * image_size);
	for (int j = 0; j < image_size; j++) {
		for (int i = 0; i < image_size / 2; i++) {
			data[j * (image_size + 2) + 2 * i] =
				(float) (my_plan.f_hat[j * image_size + i + image_size / 2].re) * norm;
			data[j * (image_size + 2) + 2 * i + 1] =
				(float) (my_plan.f_hat[j * image_size + i + image_size / 2].im) * norm * (-1.0);
		}
	}
	/** finalise the nfft plan */
	nfft_2D_finalize(&my_plan);

	if (res > 0) {
		// Gaussian low pass filter
		double sigma2 = (image_size * apix / res) * (image_size * apix / res);
		int index = 0;
		for (int j = 0; j < image_size; j++) {
			double RY2 = (j - image_size / 2) * (j - image_size / 2);
			for (int i = 0; i < image_size / 2 + 1; i++, index += 2) {
				double val = exp(-(i * i + RY2) / sigma2);
				data[index] *= val;
				data[index + 1] *= val;
			}
		}
	}
	fft->done_data();
	//fft->filter("GaussLowpass",Dict("lowpass", box*apix/res));

	fft->filter("Phase180");	// move phase origin to center of image box, instead of at corner

	return fft;
#elif defined NFFT2
	nfft_plan my_plan;			// plan for the nfft 
	int N[2], n[2];
	N[0] = image_size;
	n[0] = next_power_of_2(2 * image_size);
	N[1] = image_size;
	n[1] = next_power_of_2(2 * image_size);

	/** init an 2 dimensional plan */
	//nfft_init_2d(&my_plan,image_size,image_size,get_number_points());
	nfft_init_specific(&my_plan, 2, N, get_number_points(), n, 3,
					   PRE_PHI_HUT | PRE_PSI |
					   MALLOC_X | MALLOC_F_HAT | MALLOC_F, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
	/** init atom positions to the non-uniform nodes */
	for (int j = 0, i = 0; j < my_plan.M; j++, i += 4) {
		// FFTW and nfft use row major array layout, EMAN uses column major
		my_plan.x[2 * j + 1] = (double) (points[i] / (apix * image_size));
		my_plan.x[2 * j] = (double) (points[i + 1] / (apix * image_size));
		my_plan.f[j][0] = (double) points[i + 3];
		my_plan.f[j][1] = 0.0;
	}

	/** precompute psi, the entries of the matrix B */
	if (my_plan.nfft_flags & PRE_PSI) {
		nfft_precompute_psi(&my_plan);
		if (my_plan.nfft_flags & PRE_FULL_PSI)
			nfft_full_psi(&my_plan, pow(10, -6));
	}

	// compute the uniform Fourier transform
	nfft_transposed(&my_plan);

	// copy the Fourier transform to EMData data array
	EMData *fft = new EMData();
	fft->set_size(image_size + 2, image_size, 1);
	fft->set_complex(true);
	fft->set_ri(true);
	fft->to_zero();
	float *data = fft->get_data();
	double norm = 1.0 / (image_size * image_size);
	for (int j = 0; j < image_size; j++) {
		for (int i = 0; i < image_size / 2; i++) {
			data[j * (image_size + 2) + 2 * i] =
				(float) (my_plan.f_hat[j * image_size + i + image_size / 2][0]) * norm;
			data[j * (image_size + 2) + 2 * i + 1] =
				(float) (my_plan.f_hat[j * image_size + i + image_size / 2][1]) * norm;
		}
	}
	/** finalise the nfft plan */
	nfft_finalize(&my_plan);

	if (res > 0) {
		// Gaussian low pass filter
		double sigma2 = (image_size * apix / res) * (image_size * apix / res);
		int index = 0;
		for (int j = 0; j < image_size; j++) {
			double RY2 = (j - image_size / 2) * (j - image_size / 2);
			for (int i = 0; i < image_size / 2 + 1; i++, index += 2) {
				double val = exp(-(i * i + RY2) / sigma2);
				data[index] *= val;
				data[index + 1] *= val;
			}
		}
	}
	fft->done_data();
	//fft->filter("GaussLowpass",Dict("lowpass", box*apix/res));

	fft->filter("Phase180");	// move phase origin to center of image box, instead of at corner

	return fft;
#else
	return 0;
#endif
}

double scoring_function_f(const gsl_vector *, void *)
{
	return 0.0;
}

void scoring_function_df(const gsl_vector *, void *, gsl_vector *)
{
}

void scoring_function_fdf(const gsl_vector *, void *, double *, gsl_vector *)
{
}

bool PointArray::refine(vector < EMData * >images, string sym, OptimizedParameters optparam,
						Optimizer optimizer)
{
	const gsl_multimin_fdfminimizer_type *T;
	switch (optimizer) {
	case ConjugateGradientFletcherReeves:
		T = gsl_multimin_fdfminimizer_conjugate_fr;
		break;
	case ConjugateGradientPolakRibiere:
		T = gsl_multimin_fdfminimizer_conjugate_pr;
		break;
	case ConjugateGradientBFGS:
		T = gsl_multimin_fdfminimizer_vector_bfgs;
		break;
	//case SimplexNelderMead:
	//	T = gsl_multimin_fminimizer_nmsimplex;
	//	break;
	default:
		T = gsl_multimin_fdfminimizer_vector_bfgs;
		break;
	}

	size_t num_images = images.size();
	size_t size_par_map = get_number_points() * 4;
	size_t size_par_orientation = num_images * 3;
	size_t size_par_center = num_images * 2;
	size_t size_par_defocus = num_images * 1;
	size_t size_par_astig = num_images * 2;
	size_t size_par_bfactor = num_images * 1;
	size_t size_par_drift = num_images * 2;
	size_t size_par_scale = num_images * 1;
	size_t size_par_distortion = num_images * 2;
	size_t size_par_beamtilt = num_images * 1;	// ?

	size_t num_parameters = 0;

	if (optparam & Map)
		num_parameters += size_par_map;
	if (optparam & Orientation)
		num_parameters += size_par_orientation;
	if (optparam & Center)
		num_parameters += size_par_center;
	if (optparam & Defocus)
		num_parameters += size_par_defocus;
	if (optparam & Astigmatism)
		num_parameters += size_par_astig;
	if (optparam & BFactor)
		num_parameters += size_par_bfactor;
	if (optparam & Drift)
		num_parameters += size_par_drift;
	if (optparam & Scale)
		num_parameters += size_par_scale;
	if (optparam & Distortion)
		num_parameters += size_par_distortion;
	if (optparam & BeamTilt)
		num_parameters += size_par_beamtilt;
	if (optparam & DepthOfView)
		num_parameters += 0;

	gsl_multimin_fdfminimizer *s;
	s = gsl_multimin_fdfminimizer_alloc(T, num_parameters);
	//gsl_multimin_fminimizer *s;
	//s = gsl_multimin_fminimizer_alloc(T, num_parameters);

	// allocate the parameter array
	double *par = (double *) calloc(num_parameters, sizeof(double));
	double *par_map = 0, *par_orientation = 0, *par_center = 0, *par_defocus = 0, *par_astig = 0;
	double *par_bfactor = 0, *par_drift = 0, *par_scale = 0, *par_distortion = 0, *par_beamtilt = 0;
	size_t count = 0;
	if (optparam & Map)
		par_map = par + count, count += size_par_map;
	if (optparam & Orientation)
		par_orientation = par + count, count += size_par_orientation;
	if (optparam & Center)
		par_center = par + count, count += size_par_center;
	if (optparam & Defocus)
		par_defocus = par + count, count += size_par_defocus;
	if (optparam & Astigmatism)
		par_astig = par + count, count += size_par_astig;
	if (optparam & BFactor)
		par_bfactor = par + count, count += size_par_bfactor;
	if (optparam & Drift)
		par_drift = par + count, count += size_par_drift;
	if (optparam & Scale)
		par_scale = par + count, count += size_par_scale;
	if (optparam & Distortion)
		par_distortion = par + count, count += size_par_distortion;
	if (optparam & BeamTilt)
		par_beamtilt = par + count, count += size_par_beamtilt;
	if (optparam & DepthOfView);

	gsl_multimin_function_fdf scoring_function;
	//gsl_multimin_function scoring_function;

	scoring_function.f = &scoring_function_f;
	scoring_function.df = &scoring_function_df;
	scoring_function.fdf = &scoring_function_fdf;
	scoring_function.n = num_parameters;
	scoring_function.params = &par;

	// starting point
	gsl_vector *x0 = gsl_vector_alloc(num_parameters);
	//gsl_vector_set(x, 0, 5.0);

	double stepsize = 0.01;
	double tolerance = 1e-4;
	gsl_multimin_fdfminimizer_set(s, &scoring_function, x0, stepsize, tolerance);
	
	//double *stepsize = 0;
	//gsl_multimin_fminimizer_set(s, &scoring_function, x0, stepsize, tolerance);

	size_t iter = 0;
	int status;

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 1e-3);

		if (status == GSL_SUCCESS)
			printf("Minimum found\n");
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x0);

	return true;
}

#if defined OPTPP

/* Initializer for Rosenbrock */
void init_rosen(int n, ColumnVector& x) {}
/* Rosenbrock with analytic derivative */
void rosen(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result) {}
void update_model(int, int, ColumnVector) {}
int foo()
{
  int n = 2;
  
  static char *status_file = {"tstqnewton.out"};

//----------------------------------------------------------------------------
// 1. Quasi-Newton with trust regions
//----------------------------------------------------------------------------

  //  Create a Nonlinear problem object

  NLF1 nlp(n,rosen,init_rosen);
  
  nlp.setIsExpensive(true);
  
  //  Build a Quasi-Newton object and optimize 

  OptQNewton objfcn(&nlp,update_model);   
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setTRSize(100.);
  objfcn.optimize();
  objfcn.printStatus("Solution from quasi-newton");
  objfcn.cleanup();	 
}

#endif
