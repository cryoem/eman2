/**
 * $Id$
 */
#include "pointarray.h"
#include "util.h"

using namespace EMAN;

PointArray::PointArray()
{
	points = 0;
	n = 0;
}

PointArray::PointArray(unsigned int nn)
{
	n = nn;
	points = (double *) calloc(4 * n, sizeof(double));
}

PointArray::~PointArray()
{
	free(points);
}

void PointArray::zero()
{
	memset((void *) points, 0, 4 * n * sizeof(double));
}

PointArray *PointArray::copy()
{
	PointArray *pa2 = new PointArray();
	pa2->set_number_points(get_number_points());
	double *pa2data = pa2->get_points_array();
	memcpy(pa2data, get_points_array(), sizeof(double) * 4 * get_number_points());

	return pa2;
}

PointArray & PointArray::operator=(PointArray & pa)
{
	if (this != &pa) {
		set_number_points(pa.get_number_points());
		memcpy(get_points_array(), pa.get_points_array(), sizeof(double) * 4 * get_number_points());
	}
	return *this;
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
	set_number_points((unsigned int)(filestat.st_size / 80 + 1));
#ifdef DEBUG
	printf("PointArray::read_from_pdb(): try %4d atoms first\n", get_number_points());
#endif

	FILE *fp = fopen(file, "r");
	if(!fp) {
		fprintf(stderr,"ERROR in PointArray::read_from_pdb(): cannot open file %s\n",file);
		throw;
	}
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
#ifdef DEBUG
		printf("Atom %4d: x,y,z = %8g,%8g,%8g\te = %g\n", count, x, y, z, e);
#endif
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


void PointArray::save_to_pdb(const char *file)
{
	FILE *fp = fopen(file, "w");
	for (unsigned int i = 0; i < get_number_points(); i++) {
		fprintf(fp, "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%8s\n", i, i,
				points[4 * i], points[4 * i + 1], points[4 * i + 2], points[4 * i + 3], 0.0, " ");
	}
	fclose(fp);
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
		points[i] -= center[0];
		points[i + 1] -= center[1];
		points[i + 2] -= center[2];
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


void PointArray::mask(double rmax, double rmin)
{
	double rmax2 = rmax * rmax, rmin2 = rmin * rmin;
	PointArray *tmp = this->copy();
	double *tmp_points = tmp->get_points_array();
	unsigned int count = 0;
	for (unsigned int i = 0; i < 4 * tmp->get_number_points(); i += 4) {
		double x = tmp_points[i], y = tmp_points[i + 1], z = tmp_points[i + 2], v =
			tmp_points[i + 3];
		double r2 = x * x + y * y + z * z;
		if (r2 >= rmin2 && r2 <= rmax2) {
			points[count * 4] = x;
			points[count * 4 + 1] = y;
			points[count * 4 + 2] = z;
			points[count * 4 + 3] = v;
			count++;
		}
	}
	set_number_points(count);
	delete tmp;
}


void PointArray::mask_asymmetric_unit(const string & sym)
{
	if (sym == "c1" || sym == "C1")
		return;					// do nothing for C1 symmetry
	double alt0 = 0, alt1 = M_PI, alt2 = M_PI;
	double az0 = 0, az1 = M_PI;
	if (sym[0] == 'c' || sym[0] == 'C') {
		int nsym = atoi(sym.c_str() + 1);
		az1 = 2.0 * M_PI / nsym / 2.0;
	}
	else if (sym[0] == 'd' || sym[0] == 'D') {
		int nsym = atoi(sym.c_str() + 1);
		alt1 = M_PI / 2.0;
		alt2 = alt1;
		az1 = 2.0 * M_PI / nsym / 2.0;
	}
	else if (sym == "icos" || sym == "ICOS") {
		alt1 = 0.652358139784368185995;	// 5fold to 3fold
		alt2 = 0.55357435889704525151;	// half of edge ie. 5fold to 2fold along the edge
		az1 = 2.0 * M_PI / 5 / 2.0;
	}
	else {
		LOGERR("PointArray::set_to_asymmetric_unit(): sym = %s is not implemented yet",
			   sym.c_str());
		return;
	}
#ifdef DEBUG
	printf("Sym %s: alt0 = %8g\talt1 = %8g\talt2 = %8g\taz0 = %8g\taz1 = %8g\n", sym.c_str(), alt0*180.0/M_PI, alt1*180.0/M_PI, alt2*180.0/M_PI, az0*180.0/M_PI, az1*180.0/M_PI);
#endif

	PointArray *tmp = this->copy();
	double *tmp_points = tmp->get_points_array();
	unsigned int count = 0;
	for (unsigned int i = 0; i < 4 * tmp->get_number_points(); i += 4) {
		double x = tmp_points[i], y = tmp_points[i + 1], z = tmp_points[i + 2], v = tmp_points[i + 3];
		double az = atan2(y, x);
		double az_abs = fabs(az - az0);
		if (az_abs < (az1 - az0)) {
			double alt_max = alt1 + (alt2 - alt1) * az_abs / (az1 - az0);
			double alt = acos(z / sqrt(x * x + y * y + z * z));
			if (alt < alt_max && alt >= alt0) {
#ifdef DEBUG
				printf("Point %3d: x,y,z = %8g,%8g,%8g\taz = %8g\talt = %8g\n",i/4,x,y,z,az*180.0/M_PI, alt*180.0/M_PI);
#endif
				points[count * 4] = x;
				points[count * 4 + 1] = y;
				points[count * 4 + 2] = z;
				points[count * 4 + 3] = v;
				count++;
			}
		}
	}
	set_number_points(count);
	delete tmp;
}


void PointArray::set_from(PointArray * source, const string & sym, Transform *transform)
{
	set_from(source->get_points_array(), source->get_number_points(), sym, transform);

}

void PointArray::set_from(double *src, unsigned int num, const string & sym, Transform *xform)
{
	unsigned int nsym = xform->get_nsym(sym);

	if (get_number_points() != nsym * num)
		set_number_points(nsym * num);

	double *target = get_points_array();

	for (unsigned int s = 0; s < nsym; s++) {
		int index = s * 4 * num;
		for (unsigned int i = 0; i < 4 * num; i += 4, index += 4) {
			Vec3f v((float)src[i],(float)src[i+1],(float)src[i+2]);
			v=v*xform->get_sym(sym,s);
			target[index]  =v[0];
			target[index+1]=v[1];
			target[index+2]=v[2];
			target[index+3]=src[i+3];
		}
	}
}


void PointArray::set_from_density_map(EMData * map, int num, float thresh, float apix,
									  Density2PointsArrayAlgorithm mode)
{
	if (mode == PEAKS_SUB || mode == PEAKS_DIV) {
		// find out how many voxels are useful voxels
		int num_voxels = 0;
		int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
		EMData *tmp_map = map->copy();
		float *pd = tmp_map->get_data();
		for (int i = 0; i < nx * ny * nz; i++) {
			if (pd[i] > thresh)
				num_voxels++;
		}

		double pointvol = double (num_voxels) / double (num);
		double gauss_real_width = pow(pointvol, 1. / 3.);	// in pixels
#ifdef DEBUG
		printf("Average point range radius = %g pixels for %d points from %d used voxels\n",
			   gauss_real_width, num, num_voxels);
#endif

		double min_table_val = 1e-4;
		double max_table_x = sqrt(-log(min_table_val));	// for exp(-x*x)

		double table_step_size = 1.;	// number of steps for each pixel
		double inv_table_step_size = 1.0 / table_step_size;
		int table_size = int (max_table_x * gauss_real_width / (table_step_size) * 1.25) + 1;
		double *table = (double *) malloc(sizeof(double) * table_size);
		for (int i = 0; i < table_size; i++) {
			double x = i * table_step_size / gauss_real_width;
			table[i] = exp(-x * x);
		}

		int gbox = int (max_table_x * gauss_real_width);	// local box half size in pixels to consider for each point
		if (gbox <= 0)
			gbox = 1;

		set_number_points(num);
		for (int count = 0; count < num; count++) {
			float cmax = pd[0];
			int cmaxpos = 0;
			for (int i = 0; i < nx * ny * nz; i++) {
				if (pd[i] > cmax) {
					cmax = pd[i];
					cmaxpos = i;
				}
			}
			int iz = cmaxpos / (nx * ny), iy = (cmaxpos - iz * nx * ny) / nx, ix =
				cmaxpos - iz * nx * ny - iy * nx;
				
			// update coordinates in pixels
			points[4*count  ] = ix;
			points[4*count+1] = iy;
			points[4*count+2] = iz;
			points[4*count+3] = cmax;
#ifdef DEBUG
			printf("Point %d: val = %g\tat  %d, %d, %d\n", count, cmax, ix, iy, iz);
#endif
			
			int imin = ix - gbox, imax = ix + gbox;
			int jmin = iy - gbox, jmax = iy + gbox;
			int kmin = iz - gbox, kmax = iz + gbox;
			if (imin < 0)
				imin = 0;
			if (jmin < 0)
				jmin = 0;
			if (kmin < 0)
				kmin = 0;
			if (imax > nx)
				imax = nx;
			if (jmax > ny)
				jmax = ny;
			if (kmax > nz)
				kmax = nz;

			for (int k = kmin; k < kmax; k++) {
				int table_index_z = int (fabs(double (k - iz)) * inv_table_step_size);
				double zval = table[table_index_z];
				int pd_index_z = k * nx * ny;
				//printf("k = %8d\tx = %8g\tval = %8g\n", k, float(k-iz), zval);
				for (int j = jmin; j < jmax; j++) {
					int table_index_y = int (fabs(double (j - iy)) * inv_table_step_size);
					double yval = table[table_index_y];
					int pd_index = pd_index_z + j * nx + imin;
					for (int i = imin; i < imax; i++, pd_index++) {
						int table_index_x = int (fabs(double (i - ix)) * inv_table_step_size);
						double xval = table[table_index_x];
						if (mode == PEAKS_SUB)
							pd[pd_index] -= (float)(cmax * zval * yval * xval);
						else
							pd[pd_index] *= (float)(1.0 - zval * yval * xval);	// mode == PEAKS_DIV 
					}
				}
			}
		}
		set_number_points(num);
		tmp_map->done_data();
		delete tmp_map;
	}
	else if (mode == KMEANS) {
		set_number_points(num);
		zero();

		PointArray tmp_pa;
		tmp_pa.set_number_points(num);
		tmp_pa.zero();

		unsigned int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
		float *pd = map->get_data();

		// initialize segments with random centers at pixels with values > thresh
#ifdef DEBUG
		printf("Start initial random seeding\n");
#endif
		for (unsigned int i = 0; i < get_number_points(); i++) {
			unsigned int x, y, z;
			double v;
			do {
				x = (unsigned int) Util::get_frand(0, nx - 1);
				y = (unsigned int) Util::get_frand(0, ny - 1);
				z = (unsigned int) Util::get_frand(0, nz - 1);
				v = pd[z * nx * ny + y * nx + x];
#ifdef DEBUG
				printf("Trying Point %d: val = %g\tat  %d, %d, %d\tfrom map (%d,%d,%d)\n", i, v, x,
					   y, z, nx, ny, nz);
#endif
			} while (v <= thresh);
			points[4 * i] = (double) x;
			points[4 * i + 1] = (double) y;
			points[4 * i + 2] = (double) z;
			points[4 * i + 3] = (double) v;
#ifdef DEBUG
			printf("Point %d: val = %g\tat  %g, %g, %g\n", i, points[4 * i + 3], points[4 * i],
				   points[4 * i + 1], points[4 * i + 2]);
#endif
		}

		double min_dcen = 1e0;	// minimal mean segment center shift as convergence criterion
		double dcen = 0.0;
		int iter = 0, max_iter = 100;
		do {
#ifdef DEBUG
			printf("Iteration %3d, start\n", iter);
#endif
			double *tmp_points = tmp_pa.get_points_array();

			// reassign each pixel to the best segment
			for (unsigned int k = 0; k < nz; k++) {
				for (unsigned int j = 0; j < ny; j++) {
					for (unsigned int i = 0; i < nx; i++) {
						if (pd[k * nx * ny + j * nx + i] > thresh) {
							double min_dist = 1e60;	// just a large distance
							unsigned int min_s = 0;
							for (unsigned int s = 0; s < get_number_points(); s++) {
								double x = points[4 * s];
								double y = points[4 * s + 1];
								double z = points[4 * s + 2];
								double dist =
									(k - z) * (k - z) + (j - y) * (j - y) + (i - x) * (i - x);
								if (dist < min_dist) {
									min_dist = dist;
									min_s = s;
								}
							}
							tmp_points[4 * min_s] += i;
							tmp_points[4 * min_s + 1] += j;
							tmp_points[4 * min_s + 2] += k;
							tmp_points[4 * min_s + 3] += 1.0;
						}
					}
				}
			}
#ifdef DEBUG
			printf("Iteration %3d, finished reassigning segments\n", iter);
#endif
			// update each segment's center
			dcen = 0.0;
			for (unsigned int s = 0; s < get_number_points(); s++) {
				if (tmp_points[4 * s + 3]) {
					tmp_points[4 * s] /= tmp_points[4 * s + 3];
					tmp_points[4 * s + 1] /= tmp_points[4 * s + 3];
					tmp_points[4 * s + 2] /= tmp_points[4 * s + 3];
#ifdef DEBUG
					printf("Iteration %3d, Point %3d at %8g, %8g, %8g -> %8g, %8g, %8g\n", iter, s,
						   points[4 * s], points[4 * s + 1], points[4 * s + 2], tmp_points[4 * s],
						   tmp_points[4 * s + 1], tmp_points[4 * s + 2]);
#endif
				}
				else {			// empty segments are reseeded
					unsigned int x, y, z;
					double v;
					do {
						x = (unsigned int) Util::get_frand(0, nx - 1);
						y = (unsigned int) Util::get_frand(0, ny - 1);
						z = (unsigned int) Util::get_frand(0, nz - 1);
						v = pd[z * nx * ny + y * nx + x];
					} while (v <= thresh);
					tmp_points[4 * s] = (double) x;
					tmp_points[4 * s + 1] = (double) y;
					tmp_points[4 * s + 2] = (double) z;
					tmp_points[4 * s + 3] = (double) v;
#ifdef DEBUG
					printf
						("Iteration %3d, Point %3d reseeded from %8g, %8g, %8g -> %8g, %8g, %8g\n",
						 iter, s, points[4 * s], points[4 * s + 1], points[4 * s + 2],
						 tmp_points[4 * s], tmp_points[4 * s + 1], tmp_points[4 * s + 2]);
#endif
				}
				double dx = tmp_points[4 * s] - points[4 * s];
				double dy = tmp_points[4 * s + 1] - points[4 * s + 1];
				double dz = tmp_points[4 * s + 2] - points[4 * s + 2];
				dcen += dx * dx + dy * dy + dz * dz;
			}
			dcen = sqrt(dcen / get_number_points());
			//swap pointter, faster but risky
#ifdef DEBUG
			printf("before swap: points = %ld\ttmp_points = %ld\n", (long)get_points_array(),
				   (long)tmp_pa.get_points_array());
#endif
			double *tp = get_points_array();
			set_points_array(tmp_points);
			tmp_pa.set_points_array(tp);
			tmp_pa.zero();
#ifdef DEBUG
			printf("after  swap: points = %ld\ttmp_points = %ld\n", (long)get_points_array(),
				   (long)tmp_pa.get_points_array());
			printf("Iteration %3d, finished updating segment centers with dcen = %g pixels\n", iter,
				   dcen);
#endif

			iter++;
		} while (dcen > min_dcen && iter <= max_iter);
		map->done_data();

		sort_by_axis(2);	// x,y,z axes = 0, 1, 2
	}
	else {
		LOGERR("PointArray::set_from_density_map(): mode = %d is not implemented yet", mode);
	}
	//update to use apix and origin
	int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
	float origx, origy, origz;
	try {
		origx = map->get_attr("origin_row");
		origy = map->get_attr("origin_col");
		origz = map->get_attr("origin_sec");
	}
	catch(...) {
		origx = -nx / 2 * apix;
		origy = -ny / 2 * apix;
		origz = -nz / 2 * apix;
	}

#ifdef DEBUG
	printf("Apix = %g\torigin x,y,z = %8g,%8g,%8g\n",apix, origx, origy, origz);
#endif

	float *pd = map->get_data();
	for (unsigned int i = 0; i < get_number_points(); i++) {
#ifdef DEBUG
		printf("Point %4d: x,y,z,v = %8g,%8g,%8g,%8g",i, points[4 * i],points[4 * i + 1],points[4 * i + 2],points[4 * i + 3]);
#endif
		points[4 * i + 3] =
			pd[(int) points[4 * i + 2] * nx * ny + (int) points[4 * i + 1] * nx +
			   (int) points[4 * i]];
		points[4 * i] = points[4 * i] * apix + origx;
		points[4 * i + 1] = points[4 * i + 1] * apix + origy;
		points[4 * i + 2] = points[4 * i + 2] * apix + origz;
#ifdef DEBUG
		printf("\t->\t%8g,%8g,%8g,%8g\n",points[4 * i],points[4 * i + 1],points[4 * i + 2],points[4 * i + 3]);
#endif
	}
	map->done_data();
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
#ifdef DEBUG
	printf("PointArray::pdb2mrc_by_summation(): %d points\tmapsize = %4d\tapix = %g\tres = %g\n",get_number_points(),map_size, apix, res);
#endif
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
					pd[pd_index] += (float) (fval * zval * yval * xval);
				}
			}
		}
	}
	//for(int i=0; i<map_size*map_size; i++) pd[i]/=sqrt(M_PI);
	map->done_data();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_row", -map_size/2*apix);
	map->set_attr("origin_col", -map_size/2*apix);
	map->set_attr("origin_sec", -map_size/2*apix);

	return map;
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
				pd[pd_index] += (float)(fval * yval * xval);
			}
		}
	}
	for (int i = 0; i < image_size * image_size; i++)
		pd[i] /= sqrt(M_PI);
	proj->done_data();
	return proj;
}

void PointArray::replace_by_summation(EMData *proj, int ind, Vec3f vec, float amp, float apix, float res)
{
	double gauss_real_width = res / (M_PI);	// in Angstrom, res is in Angstrom

	double min_table_val = 1e-7;
	double max_table_x = sqrt(-log(min_table_val));	// for exp(-x*x)

	double table_step_size = 0.001;	// number of steps for each pixel
	double inv_table_step_size = 1.0 / table_step_size;
	int table_size = int (max_table_x * gauss_real_width / (apix * table_step_size) * 1.25);
	double *table = (double *) malloc(sizeof(double) * table_size);
	for (int i = 0; i < table_size; i++) {
		double x = -i * table_step_size * apix / gauss_real_width;
		table[i] = exp(-x * x)/pow((float)M_PI,.25f);
	}
	int image_size=proj->get_xsize();

	// subtract the old point
	int gbox = int (max_table_x * gauss_real_width / apix);	// local box half size in pixels to consider for each point
	if (gbox <= 0)
		gbox = 1;
	float *pd = proj->get_data();
	unsigned int s = ind;
	double xc = points[4 * s] / apix + image_size / 2;
	double yc = points[4 * s + 1] / apix + image_size / 2;
	double fval = points[4 * s + 3];
	int imin = int (xc) - gbox, imax = int (xc) + gbox;
	int jmin = int (yc) - gbox, jmax = int (yc) + gbox;
	
	if (imin < 0) imin = 0;
	if (jmin < 0) jmin = 0;
	if (imax > image_size) imax = image_size;
	if (jmax > image_size) jmax = image_size;

	for (int j = jmin; j < jmax; j++) {
		int table_index_y = int (fabs(j - yc) * inv_table_step_size);
		double yval = table[table_index_y];
		int pd_index = j * image_size + imin;
		for (int i = imin; i < imax; i++, pd_index++) {
			int table_index_x = int (fabs(i - xc) * inv_table_step_size);
			double xval = table[table_index_x];
			pd[pd_index] -= (float)(fval * yval * xval);
		}
	}
	
	// add the new point
	gbox = int (max_table_x * gauss_real_width / apix);	// local box half size in pixels to consider for each point
	if (gbox <= 0)
		gbox = 1;
	pd = proj->get_data();
	s = ind;
	xc = vec[0] / apix + image_size / 2;
	yc = vec[1] / apix + image_size / 2;
	fval = amp;
	imin = int (xc) - gbox, imax = int (xc) + gbox;
	jmin = int (yc) - gbox, jmax = int (yc) + gbox;
	
	if (imin < 0) imin = 0;
	if (jmin < 0) jmin = 0;
	if (imax > image_size) imax = image_size;
	if (jmax > image_size) jmax = image_size;

	for (int j = jmin; j < jmax; j++) {
		int table_index_y = int (fabs(j - yc) * inv_table_step_size);
		double yval = table[table_index_y];
		int pd_index = j * image_size + imin;
		for (int i = imin; i < imax; i++, pd_index++) {
			int table_index_x = int (fabs(i - xc) * inv_table_step_size);
			double xval = table[table_index_x];
			pd[pd_index] -= (float)(fval * yval * xval);
		}
	}
	

	proj->done_data();
	return;
}


EMData *PointArray::pdb2mrc_by_nfft(int , float , float )
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
	//fft->filter("eman1.LowpassGauss",Dict("lowpass", map_size*apix/res));

	fft->filter("eman1.Phase180");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_row", -map_size/2*apix);
	map->set_attr("origin_col", -map_size/2*apix);
	map->set_attr("origin_sec", -map_size/2*apix);
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
	//fft->filter("eman1.LowpassGauss",Dict("lowpass", map_size*apix/res));

	fft->filter("eman1.Phase180");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_row", -map_size / 2 * apix);
	map->set_attr("origin_col", -map_size / 2 * apix);
	map->set_attr("origin_sec", -map_size / 2 * apix);
	delete fft;
	return map;
#else
	LOGWARN("nfft support is not enabled. please recompile with nfft support enabled\n");
	return 0;
#endif
}

EMData *PointArray::projection_by_nfft(int , float , float )
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
	//fft->filter("eman1.LowpassGauss",Dict("lowpass", box*apix/res));

	fft->filter("eman1.Phase180");	// move phase origin to center of image box, instead of at corner

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
	//fft->filter("eman1.LowpassGauss",Dict("lowpass", box*apix/res));

	fft->filter("eman1.Phase180");	// move phase origin to center of image box, instead of at corner

	return fft;
#else
	LOGWARN("nfft support is not enabled. please recompile with nfft support enabled\n");
	return 0;
#endif
}

#ifdef OPTPP
#include "NLF.h"
#include "BoundConstraint.h"
#include "OptCG.h"
//#include "OptNewton.h"
#include "newmatap.h"

vector<EMData*> optdata;
PointArray *optobj;
float optpixres;

void init_opt_proj(int ndim, ColumnVector& x)
{
int i;
double *data=optobj->get_points_array();

for (i=0; i<ndim; i++) x(i+1)=data[i];
}

void calc_opt_proj(int n, const ColumnVector& x, double& fx, int& result)
{
	int i;
	PointArray pa;
	Transform xform;
	int size=optdata[0]->get_xsize();
	fx=0;
	
	for (i=0; i<optdata.size(); i++) {
		xform=(optdata[i]->get_transform());
		pa.set_from((double *)x.nric()+1,n/4,std::string("c1"),&xform);
		EMData *p=pa.projection_by_summation(size,1.0,optpixres);
		p->filter("eman1.NormalizeUnit");
		fx-=sqrt(p->cmp("Dot",EMObject(optdata[i]),Dict()));
	}
			
	result=NLPFunction;

	printf("%g\t%1.1f %1.1f %1.1f %g\t%1.1f %1.1f %1.1f %g\t%1.1f %1.1f %1.1f %g\n",fx,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12));
}

void calc_opt_projd(int mode,int n, const ColumnVector& x, double& fx, ColumnVector& gx, int& result)
{

if (mode & NLPFunction) {

}

if (mode & NLPGradient) {

}

}
#endif

void PointArray::opt_from_proj(const vector<EMData*> & proj,float pixres) {
#ifdef OPTPP
	optdata=proj;
	optobj=this;
	optpixres=pixres;
	
	FDNLF1 nlf(get_number_points()*4,calc_opt_proj,init_opt_proj);
//	NLF1 nlf(get_number_points()*4,init_opt_proj,calc_opt_projd);
	nlf.initFcn();
	
	OptCG opt(&nlf);
//	OptQNewton opt(&nlf);
	opt.setMaxFeval(2000);
	opt.optimize();
	opt.printStatus("Done");
#else 
	LOGWARN("OPT++ support not enabled.\n");
	return;
#endif
}

Vec3f PointArray::get_vector_at(int i)
{
return Vec3f((float)points[i*4],(float)points[i*4+1],(float)points[i*4+2]);
}

double PointArray::get_value_at(int i)
{
return points[i*4+3];
}

void PointArray::set_vector_at(int i,Vec3f vec,double value)
{
points[i*4]=vec[0];
points[i*4+1]=vec[1];
points[i*4+2]=vec[2];
points[i*4+3]=value;
}

void PointArray::set_vector_at(int i,vector<double> v)
{
points[i*4]  =v[0];
points[i*4+1]=v[1];
points[i*4+2]=v[2];
if (v.size()>=4) points[i*4+3]=v[3];
}
