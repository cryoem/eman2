/**
 * $Id$
 */

/*
 * Author: Wen Jiang, 08/11/2004 (jiang12@purdue.edu)
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

#include "pointarray.h"
#include <vector>
#include <cstring>

using namespace EMAN;

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
// This will do a sort in descending order
int cmp_float(const void *a, const void *b)
{
	double diff = *((float *) a) - *((float *) b);
	if (diff < 0.0)
		return 1;
	else if (diff > 0.0)
		return -1;
	else
		return 0;
}


PointArray::PointArray()
{
	points = 0;
	n = 0;
}

PointArray::PointArray( int nn)
{
	n = nn;
	points = (double *) calloc(4 * n, sizeof(double));
}

PointArray::~PointArray()
{
	if( points )
	{
		free(points);
		points = 0;
	}
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

 size_t PointArray::get_number_points() const
{
	return n;
}

void PointArray::set_number_points(size_t nn)
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

EMData *PointArray::distmx(int sortrows) {
if (n==0) return NULL;

unsigned int i,j;

EMData *ret= new EMData(n,n,1);
ret->to_zero();

for (i=0; i<n; i++) {
	for (j=i+1; j<n; j++) {
		float r=(get_vector_at(i)-get_vector_at(j)).length();
		ret->set_value_at(i,j,0,r);
		ret->set_value_at(j,i,0,r);
	}
}

if (sortrows) {
	float *data=ret->get_data();
	for (i=0; i<n; i++) qsort(&data[i*n],n,sizeof(float),cmp_float);
	ret->update();
}

return ret;
}

vector<int> PointArray::match_points(PointArray *to,float max_miss) {
EMData *mx0=distmx(1);
EMData *mx1=to->distmx(1);
unsigned int n2=mx1->get_xsize();	// same as get_number_points on to

if (max_miss<0) max_miss=(float)mx0->get_attr("sigma")/10.0f;
//printf("max error %f\n",max_miss);



vector<int> ret(n,-1);
vector<float> rete(n,0.0);
unsigned int i,j,k,l;

if (!mx0 || !mx1) {
	if (mx0) delete mx0;
	if (mx1) delete mx1;
	return ret;
}

// i iterates over elements of 'this', j looks for a match in 'to'
// k and l iterate over the individual distances
for (i=0; i<n; i++) {
	int bestn=-1;			// number of best match in mx1
	double bestd=1.0e38;		// residual error distance in best match
	for (j=0; j<n2; j++) {
		double d=0;
		int nd=0;
		for (k=l=0; k<n-1 && l<n2-1; k++,l++) {
			float d1=fabs(mx0->get_value_at(k,i)-mx1->get_value_at(l,j));
			float d2=fabs(mx0->get_value_at(k+1,i)-mx1->get_value_at(l,j));
			float d3=fabs(mx0->get_value_at(k,i)-mx1->get_value_at(l+1,j));
			float d4=fabs(mx0->get_value_at(k+1,i)-mx1->get_value_at(l+1,j));
			if (d2<d1 && d4>d2) { l--; continue; }
			if (d3<d1 && d4>d3) { k--; continue; }
			d+=d1;
			nd++;
		}
		d/=(float)nd;
//		printf("%d -> %d\t%f\t%d\n",i,j,d,nd);
		if (d<bestd) { bestd=d; bestn=j; }
	}
	ret[i]=bestn;
	rete[i]=static_cast<float>(bestd);
}

// This will remove duplicate assignments, keeping the best one
// also remove any assignments with large errors
for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
		if (rete[i]>max_miss) { ret[i]=-1; break; }
		if (i==j || ret[i]!=ret[j] || ret[i]==-1) continue;
		if (rete[i]>rete[j]) { ret[i]=-1; break; }
	}
}

delete mx0;
delete mx1;

return ret;
}

// uses bilinear least-squares to generate a transformation
// matrix for pairs of points
Transform3D *PointArray::align_2d(PointArray *to,float max_dist) {
vector<int> match=match_points(to,max_dist);
Transform3D *ret=new Transform3D();

// we use bilinear least squares to get 3/6 matrix components
unsigned int i,j;

vector<float> pts;
for (i=0; i<match.size(); i++) {
	if (match[i]==-1) continue;

//	printf("%d -> %d\n",i,match[i]);
	pts.push_back(get_vector_at(i)[0]);
	pts.push_back(get_vector_at(i)[1]);
	pts.push_back(to->get_vector_at(match[i])[0]);
}

Vec3f vx=Util::calc_bilinear_least_square(pts);

// then we get the other 3/6
for (i=j=0; i<match.size(); i++) {
	if (match[i]==-1) continue;
	pts[j*3]  =get_vector_at(i)[0];
	pts[j*3+1]=get_vector_at(i)[1];
	pts[j*3+2]=to->get_vector_at(match[i])[1];
	j++;
}

Vec3f vy=Util::calc_bilinear_least_square(pts);

//ret->set_rotation(vx[1],vx[2],0.0,vy[1],vy[2],0.0,0.0,0.0,1.0);
ret->set_rotation(vx[1],vy[1],0.0,vx[2],vy[2],0.0,0.0,0.0,1.0);
ret->set_pretrans(Vec3f(-vx[0],-vy[0],0));

return ret;
}

vector<float> PointArray::align_trans_2d(PointArray *to, int flags, float dxhint,float dyhint) {
printf("Warning, this is old code. Use align_2d.\n");

// returns (dx,dy,residual error,n points used)
// dxhint,dyhint should translate this->to
// flags : 1 - use hint values, 2 - center by strongest point (highest 'value')
int na=   get_number_points();
int nb=to->get_number_points();
if (na<=0 || nb<=0) return vector<float>(4,0);

int *a2b = (int *)malloc(na*sizeof(int));
int *b2a = (int *)malloc(nb*sizeof(int));


// find unweighted centers
float cax,cay,cbx,cby;
int i,j;

if (flags&1) {
	cbx=dxhint;
	cby=dyhint;
	cax=cay=0;
}
else if (flags&2) {
	// find the 'a' list peak
	float hia=0.0f;
	int hina=0;
	for (i=0; i<na; i++) {
		if (get_value_at(i)>hia) { hia=static_cast<float>(get_value_at(i)); hina=i; }
	}
	cax=get_vector_at(hina)[0];
	cay=get_vector_at(hina)[1];

	// find the 'b' list peak
	float hib=0;
	int hinb=0;
	for (i=0; i<na; i++) {
		if (to->get_value_at(i)>hib) { hib=static_cast<float>(to->get_value_at(i)); hinb=i; }
	}
	cbx=to->get_vector_at(hinb)[0];
	cby=to->get_vector_at(hinb)[1];
}
else {
	cax=cay=cbx=cby=0;

	for (i=0; i<na; i++) { cax+=get_vector_at(i)[0]; cay+=get_vector_at(i)[1]; }
	cax/=(float)na;
	cay/=(float)na;

	for (i=0; i<nb; i++) { cbx+=to->get_vector_at(i)[0]; cby+=to->get_vector_at(i)[1]; }
	cbx/=(float)nb;
	cby/=(float)nb;
}

Vec3f offset(cbx-cax,cby-cay,0);

// find the nearest point for each x point, taking the estimated centers into account
for (i=0; i<na; i++) {
	float rmin=1.0e30f;
	for (j=0; j<nb; j++) {
		float r=(get_vector_at(i)+offset-to->get_vector_at(j)).length();
		if (r<rmin) { a2b[i]=j; rmin=r; }
	}
}

// find the nearest point for each y point
for (i=0; i<nb; i++) {
	float rmin=1.0e30f;
	for (j=0; j<na; j++) {
		float r=(get_vector_at(j)+offset-to->get_vector_at(i)).length();
		if (r<rmin) { b2a[i]=j; rmin=r; }
	}
}

// now keep only points where x->y matches y->x
for (i=0; i<na; i++) {
	if (a2b[i]<0) continue;
	if (b2a[a2b[i]]!=i) { printf(" #%d!=%d# ",b2a[a2b[i]],i);  b2a[a2b[i]]=-1; a2b[i]=-1; }
	printf("%d->%d  ",i,a2b[i]);
}
printf("\n");

for (i=0; i<nb; i++) {
	if (b2a[i]<0) continue;
	if (a2b[b2a[i]]!=i) { a2b[b2a[i]]=-1; b2a[i]=-1; }
	printf("%d->%d  ",i,b2a[i]);
}
printf("\n");

// Compute the average translation required to align the points
float dx=0,dy=0,dr=0,nr=0;
for (i=0; i<na; i++) {
	if (a2b[i]==-1) continue;
	dx+=to->get_vector_at(a2b[i])[0]-get_vector_at(i)[0];
	dy+=to->get_vector_at(a2b[i])[1]-get_vector_at(i)[1];
	nr+=1.0;
}
//printf("%f %f %f\n",dx,dy,nr);
if (nr<2) return vector<float>(4,0);
dx/=nr;
dy/=nr;

// Compute the residual error
for (i=0; i<na; i++) {
	if (i==-1  || a2b[i]==-1) continue;
	dr+=(to->get_vector_at(a2b[i])-get_vector_at(i)-Vec3f(dx,dy,0)).length();
}
dr/=nr;

free(a2b);
free(b2a);
vector<float> ret(4);
ret[0]=dx;
ret[1]=dy;
ret[2]=dr;
ret[3]=(float)nr;
return ret;
}


bool PointArray::read_from_pdb(const char *file)
{
	struct stat filestat;
	stat(file, &filestat);
	set_number_points(( int)(filestat.st_size / 80 + 1));
#ifdef DEBUG
	printf("PointArray::read_from_pdb(): try %4d atoms first\n", get_number_points());
#endif

	FILE *fp = fopen(file, "r");
	if(!fp) {
		fprintf(stderr,"ERROR in PointArray::read_from_pdb(): cannot open file %s\n",file);
		throw;
	}
	char s[200];
	size_t count = 0;
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
			set_number_points(2 * (count + 1));    //makes sure point array is big enough
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
	for ( size_t i = 0; i < get_number_points(); i++) {
		fprintf(fp, "ATOM  %5lu  CA  ALA A%4lu    %8.3f%8.3f%8.3f%6.2f%6.2f%8s\n", i, i,
				points[4 * i], points[4 * i + 1], points[4 * i + 2], points[4 * i + 3], 0.0, " ");
	}
	fclose(fp);
}


FloatPoint PointArray::get_center()
{
	double xc, yc, zc;
	xc = yc = zc = 0.0;
	double norm = 0.0;
	for ( size_t i = 0; i < 4 * get_number_points(); i += 4) {
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
	for ( size_t i = 0; i < 4 * get_number_points(); i += 4) {
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
	for ( size_t i = 0; i < 4 * get_number_points(); i += 4) {
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
	 int count = 0;
	for ( size_t i = 0; i < 4 * tmp->get_number_points(); i += 4) {
		double x = tmp_points[i], y = tmp_points[i + 1], z = tmp_points[i + 2], v =
			tmp_points[i + 3];
		double r2 = x * x + y * y + z * z;
		if (r2 >= rmin2 && r2 <= rmax2) {
			points[count * 4] = x;
			points[count * 4 + 1] = y;
			points[count * 4 + 2] = z;
			points[count * 4 + 3] = v;
			++count;
		}
	}
	set_number_points(count);
	if( tmp )
	{
		delete tmp;
		tmp = 0;
	}
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
	 int count = 0;
	for ( size_t i = 0; i < 4 * tmp->get_number_points(); i += 4) {
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
	if( tmp )
	{
		delete tmp;
		tmp = 0;
	}
}

vector<float> PointArray::get_points() {
vector<float> ret;
for (unsigned int i=0; i<n; i++) {
	ret.push_back((float)points[i*4]);
	ret.push_back((float)points[i*4+1]);
	ret.push_back((float)points[i*4+2]);
}

return ret;
}

void PointArray::transform(Transform3D xf) {

for ( unsigned int i = 0; i < 4 * n; i += 4) {
	Vec3f v((float)points[i],(float)points[i+1],(float)points[i+2]);
	v=v*xf;
	points[i]  =v[0];
	points[i+1]=v[1];
	points[i+2]=v[2];
}

}

void PointArray::right_transform(const Transform& transform) {
	for ( unsigned int i = 0; i < 4 * n; i += 4) {
		Transform s = transform.transpose();
		Vec3f v((float)points[i],(float)points[i+1],(float)points[i+2]);
		v= s*v;
		points[i]  =v[0];
		points[i+1]=v[1];
		points[i+2]=v[2];
	}

}
void PointArray::set_from(PointArray * source, const string & sym, Transform3D *transform)
{
	set_from(source->get_points_array(), source->get_number_points(), sym, transform);

}

void PointArray::set_from(double *src,  int num, const string & sym, Transform3D *xform)
{
	 int nsym = xform->get_nsym(sym);

	if (get_number_points() != (size_t)nsym * num)
		set_number_points((size_t)nsym * num);

	double *target = get_points_array();

	for ( int s = 0; s < nsym; s++) {
		int index = s * 4 * num;
		for ( int i = 0; i < 4 * num; i += 4, index += 4) {
			Vec3f v((float)src[i],(float)src[i+1],(float)src[i+2]);
			v=v*xform->get_sym(sym,s);
			target[index]  =v[0];
			target[index+1]=v[1];
			target[index+2]=v[2];
			target[index+3]=src[i+3];
		}
	}
}

void PointArray::set_from(vector<float> pts) {
	set_number_points(pts.size()/4);
	for (unsigned int i=0; i<pts.size(); i++) points[i]=pts[i];

}

void PointArray::set_from_density_map(EMData * map, int num, float thresh, float apix,
									  Density2PointsArrayAlgorithm mode)
{
	if (mode == PEAKS_SUB || mode == PEAKS_DIV) {
		// find out how many voxels are useful voxels
		int num_voxels = 0;
		int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
		size_t size = nx * ny * nz;
		EMData *tmp_map = map->copy();
		float *pd = tmp_map->get_data();
		for (size_t i = 0; i < size; ++i) {
			if (pd[i] > thresh)
				++num_voxels;
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
			for (size_t i = 0; i < size; ++i) {
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

			for (int k = kmin; k < kmax; ++k) {
				int table_index_z = int (fabs(double (k - iz)) * inv_table_step_size);
				double zval = table[table_index_z];
				size_t pd_index_z = k * nx * ny;
				//printf("k = %8d\tx = %8g\tval = %8g\n", k, float(k-iz), zval);
				for (int j = jmin; j < jmax; ++j) {
					int table_index_y = int (fabs(double (j - iy)) * inv_table_step_size);
					double yval = table[table_index_y];
					size_t pd_index = pd_index_z + j * nx + imin;
					for (int i = imin; i < imax; ++i, ++pd_index) {
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
		tmp_map->update();
		if( tmp_map )
		{
			delete tmp_map;
			tmp_map = 0;
		}
	}
	else if (mode == KMEANS) {
		set_number_points(num);
		zero();

		PointArray tmp_pa;
		tmp_pa.set_number_points(num);
		tmp_pa.zero();

		 int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
		float *pd = map->get_data();

		// initialize segments with random centers at pixels with values > thresh
#ifdef DEBUG
		printf("Start initial random seeding\n");
#endif
		for ( size_t i = 0; i < get_number_points(); i++) {
			 int x, y, z;
			double v;
			do {
				x = ( int) Util::get_frand(0, nx - 1);
				y = ( int) Util::get_frand(0, ny - 1);
				z = ( int) Util::get_frand(0, nz - 1);
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
			for ( int k = 0; k < nz; k++) {
				for ( int j = 0; j < ny; j++) {
					for ( int i = 0; i < nx; i++) {
						size_t idx = k * nx * ny + j * nx + i;
						if (pd[idx] > thresh) {
							double min_dist = 1e60;	// just a large distance
							 int min_s = 0;
							for ( size_t s = 0; s < get_number_points(); ++s) {
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
			for ( size_t s = 0; s < get_number_points(); ++s) {
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
					 int x, y, z;
					double v;
					do {
						x = ( int) Util::get_frand(0, nx - 1);
						y = ( int) Util::get_frand(0, ny - 1);
						z = ( int) Util::get_frand(0, nz - 1);
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
		map->update();

		sort_by_axis(2);	// x,y,z axes = 0, 1, 2
	}
	else {
		LOGERR("PointArray::set_from_density_map(): mode = %d is not implemented yet", mode);
	}
	//update to use apix and origin
	int nx = map->get_xsize(), ny = map->get_ysize(), nz = map->get_zsize();
	float origx, origy, origz;
	try {
		origx = map->get_attr("origin_x");
		origy = map->get_attr("origin_y");
		origz = map->get_attr("origin_z");
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
	for ( size_t i = 0; i < get_number_points(); ++i) {
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
	map->update();
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
	vector<double> table;
	table.resize(table_size);
	//double *table = (double *) malloc(sizeof(double) * table_size);
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
	for ( size_t s = 0; s < get_number_points(); ++s) {
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
			if ( table_index_z >= table.size() ) continue;
			double zval = table[table_index_z];
			size_t pd_index_z = k * map_size * map_size;
			for (int j = jmin; j < jmax; j++) {
				int table_index_y = int (fabs(j - yc) * inv_table_step_size);
				if ( table_index_y >= table.size() ) continue;
				double yval = table[table_index_y];
				size_t pd_index = pd_index_z + j * map_size + imin;
				for (int i = imin; i < imax; i++, pd_index++) {
					int table_index_x = int (fabs(i - xc) * inv_table_step_size);
					if ( table_index_x >= table.size() ) continue;
					double xval = table[table_index_x];
					pd[pd_index] += (float) (fval * zval * yval * xval);
				}
			}
		}
	}
	//for(int i=0; i<map_size*map_size; i++) pd[i]/=sqrt(M_PI);
	map->update();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_x", -map_size/2*apix);
	map->set_attr("origin_y", -map_size/2*apix);
	map->set_attr("origin_z", -map_size/2*apix);

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
	for ( size_t s = 0; s < get_number_points(); ++s) {
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
	proj->update();
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
	 int s = ind;
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

	proj->update();
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

	// low pass processor
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
	fft->update();
	//fft->process_inplace("eman1.filter.lowpass.gaussian",Dict("lowpass", map_size*apix/res));

	fft->process_inplace("xform.phaseorigin.tocorner");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_x", -map_size/2*apix);
	map->set_attr("origin_y", -map_size/2*apix);
	map->set_attr("origin_z", -map_size/2*apix);
	if( fft )
	{
		delete fft;
		fft = 0;
	}
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

	// low pass processor
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
	fft->update();
	//fft->process_inplace("eman1.filter.lowpass.gaussian",Dict("lowpass", map_size*apix/res));

	fft->process_inplace("xform.phaseorigin.tocenter");	// move phase origin to center of image map_size, instead of at corner
	EMData *map = fft->do_ift();
	map->set_attr("apix_x", apix);
	map->set_attr("apix_y", apix);
	map->set_attr("apix_z", apix);
	map->set_attr("origin_x", -map_size / 2 * apix);
	map->set_attr("origin_y", -map_size / 2 * apix);
	map->set_attr("origin_z", -map_size / 2 * apix);
	if( fft )
	{
		delete fft;
		fft = 0;
	}
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
		// Gaussian low pass processor
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
	fft->update();
	//fft->process_inplace("eman1.filter.lowpass.gaussian",Dict("lowpass", box*apix/res));

	fft->process_inplace("xform.phaseorigin.tocenter");	// move phase origin to center of image box, instead of at corner

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
		// Gaussian low pass process
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
	fft->update();
	//fft->process_inplace("eman1.filter.lowpass.gaussian",Dict("lowpass", box*apix/res));

	fft->process_inplace("xform.phaseorigin.tocenter");	// move phase origin to center of image box, instead of at corner

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
	Transform3D xform;
	int size=optdata[0]->get_xsize();
	fx=0;

	for (i=0; i<optdata.size(); i++) {
		xform=(optdata[i]->get_transform());
		pa.set_from((double *)x.nric()+1,n/4,std::string("c1"),&xform);
		EMData *p=pa.projection_by_summation(size,1.0,optpixres);
		p->process_inplace("normalize.unitlen");
		fx-=sqrt(p->cmp("dot",EMObject(optdata[i]),Dict()));
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
	(void)proj;		//suppress warning message
	(void)pixres;	//suppress warning message
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
