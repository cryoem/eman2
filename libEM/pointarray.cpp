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
#include "util.h"
#include "vec3.h"
#include <vector>
#include <cstring>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#ifdef __APPLE__
	typedef unsigned int uint;
#endif	//__APPLE__

#ifdef _WIN32
	typedef unsigned int uint;
#endif	//_WIN32

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
	bfactor = 0;
	n = 0;
	
	adist=0;
	aang=0;
	adihed=0;
	
	map=gradx=grady=gradz=0;
}

PointArray::PointArray( int nn)
{
	n = nn;
	points = (double *) calloc(4 * n, sizeof(double));
	
	adist=0;
	aang=0;
	adihed=0;
	map=gradx=grady=gradz=0;
}

PointArray::~PointArray()
{
	if( points )
	{
		free(points);
		points = 0;
	}
	
	if (adist) free(adist);
	if (aang) free(aang);
	if (adihed) free(adihed);
//	if (map!=0) delete map;
	if (gradx!=0) delete gradx;
	if (grady!=0) delete grady;
	if (gradz!=0) delete gradz;
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
		bfactor = (double *) realloc(bfactor, n * sizeof(double));
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
Transform *PointArray::align_2d(PointArray *to,float max_dist) {
vector<int> match=match_points(to,max_dist);
Transform *ret=new Transform();

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
//ret->set_rotation(vx[1],vy[1],0.0,vx[2],vy[2],0.0,0.0,0.0,1.0);
//ret->set_pretrans(Vec3f(-vx[0],-vy[0],0));

ret->set(0, 0, vx[1]);
ret->set(0, 1, vy[1]);
ret->set(0, 2, 0.0f);
ret->set(1, 0, vx[2]);
ret->set(1, 1, vy[2]);
ret->set(1, 2, 0.0f);
ret->set(2, 0, 0.0f);
ret->set(2, 1, 0.0f);
ret->set(2, 2, 1.0f);
ret->set_pre_trans(Vec3f(-vx[0],-vy[0],0));

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
	printf("PointArray::read_from_pdb(): try %4lu atoms first\n", get_number_points());
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
		case 'K':
			e = 19.0;
			break;
		default:
			fprintf(stderr, "Unknown atom %c%c\n", ctt, ctt2);
			e = 0;
		}
		if (e == 0)
			continue;

		float x, y, z, q;
		sscanf(&s[28], " %f %f %f", &x, &y, &z);
		sscanf(&s[60], " %f", &q);

		if (count + 1 > get_number_points())
			set_number_points(2 * (count + 1));    //makes sure point array is big enough
		
#ifdef DEBUG
		printf("Atom %4lu: x,y,z = %8g,%8g,%8g\te = %g\n", count, x, y, z, e);
#endif
		points[4 * count] = x;
		points[4 * count + 1] = y;
		points[4 * count + 2] = z;
		points[4 * count + 3] = e;
		bfactor[count] = q;
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
	fprintf(fp, "TER   %5lu      ALA A%4lu\nEND", get_number_points(),get_number_points()-1);
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
				printf("Point %3lu: x,y,z = %8g,%8g,%8g\taz = %8g\talt = %8g\n",i/4,x,y,z,az*180.0/M_PI, alt*180.0/M_PI);
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

void PointArray::transform(const Transform& xf) {

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
void PointArray::set_from(PointArray * source, const string & sym, Transform *transform)
{
	set_from(source->get_points_array(), source->get_number_points(), sym, transform);

}

void PointArray::set_from(double *src,  int num, const string & sym, Transform *xform)
{
	int nsym = xform->get_nsym(sym);
	if (xform==0){
		Transform tr;
		xform=&tr;
	}

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
		size_t size = (size_t)nx * ny * nz;
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
				printf("Trying Point %lu: val = %g\tat  %d, %d, %d\tfrom map (%d,%d,%d)\n", i, v, x,
					   y, z, nx, ny, nz);
#endif
			} while (v <= thresh);
			points[4 * i] = (double) x;
			points[4 * i + 1] = (double) y;
			points[4 * i + 2] = (double) z;
			points[4 * i + 3] = (double) v;
#ifdef DEBUG
			printf("Point %lu: val = %g\tat  %g, %g, %g\n", i, points[4 * i + 3], points[4 * i],
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
					printf("Iteration %3d, Point %3lu at %8g, %8g, %8g -> %8g, %8g, %8g\n", iter, s,
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
						("Iteration %3d, Point %3lu reseeded from %8g, %8g, %8g -> %8g, %8g, %8g\n",
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
		printf("Point %4lu: x,y,z,v = %8g,%8g,%8g,%8g",i, points[4 * i],points[4 * i + 1],points[4 * i + 2],points[4 * i + 3]);
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

/** Updates the dist,ang,dihed parameters **/
void PointArray::sim_updategeom() {
	if (!adist) adist=(double *)malloc(sizeof(double)*n);
	if (!aang) aang=(double *)malloc(sizeof(double)*n);
	if (!adihed) adihed=(double *)malloc(sizeof(double)*n);
	
	for (size_t ii=0; ii<n; ii++) {
		// how expensive is % ?  Worth replacing ?
		int ib=4*((ii+n-1)%n);		// point before i with wraparound
		int ibb=4*((ii+n-2)%n);	// 2 points before i with wraparound
		int ia=4*((ii+1)%n);		// 1 point after
		int i=4*ii;

		Vec3f a(points[ib]-points[ibb],points[ib+1]-points[ibb+1],points[ib+2]-points[ibb+2]);  // -2 -> -1
		Vec3f b(points[i]-points[ib],points[i+1]-points[ib+1],points[i+2]-points[ib+2]);		// -1 -> 0
		Vec3f c(points[ia]-points[i],points[ia+1]-points[i+1],points[ia+2]-points[i+2]);		// 0 -> 1
		adist[ii]=b.length();
		
		// angle
		aang[ii]=b.dot(c);
		if (aang[ii]!=0) {
			aang[ii]/=(adist[ii]*c.length());
			if (aang[ii]>=1.0) aang[ii]=0.0;
			else if (aang[ii]<=-1.0) aang[ii]=M_PI;
			else aang[ii]=acos(aang[ii]);
		}

		// dihedral
		Vec3f cr1=a.cross(b);
		Vec3f cr2=b.cross(c);
		double denom=cr1.length()*cr2.length();
		if (denom==0) adihed[ii]=0;
		else {
			double tmp=cr1.dot(cr2)/(denom);
			if (tmp>1) tmp=1;
			if (tmp<-1) tmp=-1;
			adihed[ii]=acos(tmp);
		}
		
//		if (std::isnan(ang[ii])) ang[ii]=0;
		
	}
}

double PointArray::sim_potential() {
	double ret=0;
	sim_updategeom();
	
	if (map &&mapc) {
		for (size_t i=0; i<n; i++) ret+=sim_pointpotential(adist[i],aang[i],adihed[i])-mapc*map->sget_value_at_interp(points[i*4]/apix+centx,points[i*4+1]/apix+centy,points[i*4+2]/apix+centz);
	}
	else {
		for (size_t i=0; i<n; i++) ret+=sim_pointpotential(adist[i],aang[i],adihed[i]);
	}

#ifdef _WIN32
//	if (_isnan(ret/n))
#else
//	if (std::isnan(ret/n))
#endif
//		printf("%f             %f\n",ret,n);
	
	return ret/n;
}

// potential for a single point. Note that if a point moves, it will impact energies +-2 from its position. This function computes only for the point i
double PointArray::sim_potentiald(int ind) {
	if (!adist) sim_updategeom();		// wasteful, but only once
	
//	if (i<0 || i>=n) throw InvalidParameterException("Point number out of range");
	size_t i;
	if (ind<0)
		i=ind-n*(ind/n-1);
	else
		i=ind;
	if (i>=n) i=i%n;
	
	// how expensive is % ?  Worth replacing ?
	int ib=4*((i+n-1)%n);		// point before i with wraparound
	int ibb=4*((i+n-2)%n);	// 2 points before i with wraparound
	int ia=4*((i+1)%n);		// 1 point after
	int ii=i*4;
	
	Vec3f a(points[ib]-points[ibb],points[ib+1]-points[ibb+1],points[ib+2]-points[ibb+2]);  		// -2 -> -1
	Vec3f b(points[ii]-points[ib],points[ii+1]-points[ib+1],points[ii+2]-points[ib+2]);		// -1 -> 0
	Vec3f c(points[ia]-points[ii],points[ia+1]-points[ii+1],points[ia+2]-points[ii+2]);		// 0 -> 1
	double dist=b.length();
	adist[i]=dist;
	// Angle, tests should avoid isnan being necessary
	double ang=b.dot(c);
	if (ang!=0.0) {					// if b.dot(c) is 0, we set it to the last determined value...
		ang/=(dist*c.length());
		if (ang>1.0) ang=1.0;		// should never happen, but just in case of roundoff error
		if (ang<-1.0) ang=-1.0;
		ang=acos(ang);
	}
	else ang=aang[i];
	aang[i]=ang;
	
	// Dihedral
	Vec3f cr1=a.cross(b);
	Vec3f cr2=b.cross(c);
	double dihed;
	double denom=cr1.length()*cr2.length();
	if (denom==0) dihed=adihed[i];						// set the dihedral to the last determined value if indeterminate
	else {
		dihed=cr1.dot(cr2)/denom;
		if (dihed>1.0) dihed=1.0;				// should never happen, but just in case of roundoff error
		if (dihed<-1.0) dihed=-1.0;
		dihed=acos(dihed);
	}
	adihed[i]=dihed;
//*	Do not need for small amount of points
	// Distance to the closest neighbor
	double mindist=10000;
	for (size_t j=0;j<n;j++){
		if(j==i)
			continue;
		int ja=4*j;
		Vec3f d(points[ii]-points[ja],points[ii+1]-points[ja+1],points[ii+2]-points[ja+2]);
		double jdst=d.length();
		if(jdst<mindist)
			mindist=jdst;
	}
	double distpen=0;
	if (mindist<mindistc)
		distpen=distpenc/mindist;
//*/	
	
	
//	if (std::isnan(dist) || std::isnan(ang) || std::isnan(dihed)) printf("%d\t%g\t%g\t%g\t%g\t%g\t%g\n",i,dist,ang,dihed,b.length(),c.length(),b.dot(c)/(dist*c.length()));
// 	if (std::isnan(dihed)) dihed=dihed0;
// 	if (std::isnan(ang)) ang=0;
// 	if (std::isnan(dist)) dist=3.3;
//	if (isnan(dihed)) dihed=dihed0;
//	if (isnan(ang)) ang=0;
	if (map && mapc) {
		return distpen+sim_pointpotential(dist,ang,dihed)-mapc*map->sget_value_at_interp(points[ii]/apix+centx,points[ii+1]/apix+centy,points[ii+2]/apix+centz);
	}
	return sim_pointpotential(dist,ang,dihed);
}

// Computes a potential for a point and +-2 nearest neighbors under perturbation of the central point location
double PointArray::sim_potentialdxyz(int i, double dx, double dy, double dz) {
	Vec3f old(points[i*4],points[i*4+1],points[i*4+2]);
	double potd=0.;
		
	points[i*4]=old[0]+dx;
	points[i*4+1]=old[1]+dy;
	points[i*4+2]=old[2]+dz;
 	for (int ii=i-2; ii<=i+2; ii++) potd+=sim_potentiald(ii);
//	potd=potential();
	points[i*4]=old[0];
	points[i*4+1]=old[1];
	points[i*4+2]=old[2];
	
	return potd;
}

double PointArray::calc_total_length(){
	double dist=0;
	for(size_t i=0; i<n; i++){
		int k=(i+1)%n;
		double d=(points[i*4]-points[k*4])*(points[i*4]-points[k*4])+(points[i*4+1]-points[k*4+1])*(points[i*4+1]-points[k*4+1])+(points[i*4+2]-points[k*4+2])*(points[i*4+2]-points[k*4+2]);
		d=sqrt(d);
		dist+=d;
	}
	return dist;
}

// Computes a gradient of the potential for a single point, including impact on +-2 nearest neighbors
Vec3f PointArray::sim_descent(int i) {
	Vec3f old(points[i*4],points[i*4+1],points[i*4+2]);
	double pot=0.,potx=0.,poty=0.,potz=0.;
	double stepsz=0.01;
	
	for (int ii=i-2; ii<=i+2; ii++) pot+=sim_potentiald(ii);
// 	pot=potential();
	
	points[i*4]=old[0]+stepsz;
	for (int ii=i-2; ii<=i+2; ii++) potx+=sim_potentiald(ii);
// 	potx=potential();
	points[i*4]=old[0];

	points[i*4+1]=old[1]+stepsz;
	for (int ii=i-2; ii<=i+2; ii++) poty+=sim_potentiald(ii);
// 	poty=potential();
	points[i*4+1]=old[1];

	points[i*4+2]=old[2]+stepsz;
	for (int ii=i-2; ii<=i+2; ii++) potz+=sim_potentiald(ii);
// 	potz=potential();
	points[i*4+2]=old[2];
	
	//	printf("%1.4g\t%1.4g\t%1.4g\t%1.4g\t%1.4g\t%1.4g\t%1.4g\n",pot,potx,poty,potz,(pot-potx),(pot-poty),(pot-potz));
	
// 	if (pot==potx) potx=pot+1000000.0;
// 	if (pot==potx) potx=pot+1000000.0;
// 	if (pot==potx) potx=pot+1000000.0;
	return Vec3f((pot-potx),(pot-poty),(pot-potz));
}

/** Takes a step to minimize the potential **/ 
void PointArray::sim_minstep(double maxshift) { 
	vector<Vec3f> shifts;
	
	double max=0.0;
	double mean=0.0;
	for (size_t i=0; i<n; i++) {
		if (oldshifts.size()==n) shifts.push_back((sim_descent(i)+oldshifts[i])/2.0);
		else shifts.push_back(sim_descent(i));
		float len=shifts[i].length();
		if (len>max) max=len;
		mean+=len;
	}
	oldshifts=shifts;
	
//	printf("max vec %1.2f\tmean %1.3f\n",max,mean/n);
	
	for (size_t i=0; i<n; i++) {
//		if (std::isnan(shifts[i][0]) ||std::isnan(shifts[i][1]) ||std::isnan(shifts[i][2])) { printf("Nan: %d\n",i); shifts[i]=Vec3f(max,max,max); }
		points[i*4]+=shifts[i][0]*maxshift/max;
		points[i*4+1]+=shifts[i][1]*maxshift/max;
		points[i*4+2]+=shifts[i][2]*maxshift/max;
//		printf("%d. %1.2f\t%1.2f\t%1.2f\n",i,shifts[i][0]*maxshift/max,shifts[i][1]*maxshift/max,shifts[i][2]*maxshift/max);
	}
}

/** Takes a step to minimize the potential **/ 
void PointArray::sim_minstep_seq(double meanshift) {
	/*
	// Try to minimize potential globally
	boost::mt19937 rng; 
	boost::normal_distribution<> nd(0.0, 20.0);
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);
	double *oldpts=new double[4*get_number_points()];
	double *bestpts=new double[4*get_number_points()];
	double best_pot,new_pot;
	memcpy(oldpts, get_points_array(), sizeof(double) * 4 * get_number_points());
	memcpy(bestpts, get_points_array(), sizeof(double) * 4 * get_number_points());
	best_pot=sim_potential();
	double disttmp=0;
	for (int k=0; k<n; k++) disttmp+=adist[k];
	best_pot+=distc*pow((disttmp-336*3.3),2.0);
	for (int i=0; i<1000; i++){
		for (int j=0; j<n; j++){
			points[4*j]=oldpts[4*j]+var_nor();
			points[4*j+1]=oldpts[4*j+1]+var_nor();
			points[4*j+2]=oldpts[4*j+2]+var_nor();
		}
		new_pot=sim_potential();
		disttmp=0;
		for (int k=0; k<n; k++) disttmp+=adist[k];
		new_pot+=distc*pow((disttmp-336*3.3),2.0);
		if (new_pot<best_pot){
			memcpy(bestpts, get_points_array(), sizeof(double) * 4 * get_number_points());
			best_pot=new_pot;
			printf("%f\t",best_pot);
		}
	}
	memcpy(get_points_array(),bestpts, sizeof(double) * 4 * get_number_points());
	
	delete []oldpts;
	delete []bestpts;
	*/
	// we compute 10 random gradients and use these to adjust stepsize
	double mean=0.0;
	for (int i=0; i<10; i++) {
//		Vec3f shift=sim_descent(random()%n);
		Vec3f shift=sim_descent(Util::get_irand(0,n-1));
		mean+=shift.length();
	}
	mean/=10.0;
	double stepadj=meanshift/mean;
//	printf("\t%1.4g\n",stepadj);

	// Now we go through all points sequentially and move each downhill
	// The trick here is the sequential part, as each point is impacted by the point already moved before it.
	// This may create a "seam" at the first point which won't be adjusted to compensate for the last point (wraparound)
	// until the next cycle
	Vec3f oshifts;
	for (size_t ii=0; ii<n; ii++) {
		size_t i=2*(ii%(n/2))+2*ii/n;	// this maps a linear sequence to an all-even -> all odd sequence
		Vec3f shift,d;
		if (ii==n) {
			d=sim_descent(i);
			shift=(d+oshifts)/2.0;
			oshifts=d;
		}
		else {
			shift=sim_descent(i);
			oshifts=shift;
		}
		
//		double p2=potential();
//		double pot=sim_potentialdxyz(i,0.0,0.0,0.0);
//		double pots=sim_potentialdxyz(i,shift[0]*stepadj,shift[1]*stepadj,shift[2]*stepadj);

		// only step if it actually improves the potential for this particle (note that it does not need to improve the overall potential)
//		if (pots<pot) {
			points[i*4]+=shift[0]*stepadj;
			points[i*4+1]+=shift[1]*stepadj;
			if (!map2d)
				points[i*4+2]+=shift[2]*stepadj;
//			printf("%d. %1.4g -> %1.4g  %1.3g %1.3g %1.3g %1.3g\n",i,pot,pots,shift[0],shift[1],shift[2],stepadj);
//			if (potential()>p2) printf("%d. %1.4g %1.4g\t%1.4g %1.4g\n",i,pot,pots,p2,potential());
//		}
		
	}	
}


void PointArray::sim_rescale() {
	double meandist=0.0;
	double max=0.0,min=dist0*1000.0;
	
	for (size_t ii=0; ii<n; ii++) {
		int ib=4*((ii+n-1)%n);		// point before i with wraparound
		int i=4*ii;

		Vec3f b(points[i]-points[ib],points[i+1]-points[ib+1],points[i+2]-points[ib+2]);		// -1 -> 0
		double len=b.length();
		meandist+=len;
		max=len>max?len:max;
		min=len<min?len:min;
	}
	meandist/=n;
	double scale=dist0/meandist;
	
	printf("mean = %1.3f rescaled: %1.3f - %1.3f\n",meandist,min*scale,max*scale);
	
	for (size_t i=0; i<n; i++) {
		points[i*4]*=scale;
		points[i*4+1]*=scale;
		points[i*4+2]*=scale;
	}

}	
void PointArray::sim_printstat() {
	sim_updategeom();	
	
	double mdist=0.0,mang=0.0,mdihed=0.0;
	double midist=1000.0,miang=M_PI*2,midihed=M_PI*2;
	double madist=0.0,maang=0.0,madihed=0.0;
	double mmap=0.0;
	if (map &&mapc) {
		for (size_t i=0; i<n; i++) {
			double m=map->sget_value_at_interp(points[i*4]/apix+centx,points[i*4+1]/apix+centy,points[i*4+2]/apix+centz);
			mmap+=m;
// 			printf("%f,%f,%f\t %f\n",points[i*4]/apix+centx,points[i*4+1]/apix+centy,points[i*4+2]/apix+centz,m);
		}
		mmap/=n;
	}
	for (size_t i=0; i<n; i++) {
		mdist+=adist[i];
		mang+=aang[i];
		mdihed+=adihed[i];
		
		

		midist=adist[i]<midist?adist[i]:midist;
		madist=adist[i]>madist?adist[i]:madist;
		
		miang=aang[i]<miang?aang[i]:miang;
		maang=aang[i]>maang?aang[i]:maang;
		
		midihed=adihed[i]<midihed?adihed[i]:midihed;
		madihed=adihed[i]>madihed?adihed[i]:madihed;
	}
	double p=sim_potential();
	double anorm = 180.0/M_PI;
	printf(" potential: %1.1f\t map: %1.2f\tdist: %1.2f || %1.2f / %1.2f / %1.2f\tang: %1.2f / %1.2f / %1.2f\tdihed: %1.2f / %1.2f / %1.2f  ln=%1.1f\n",p,mmap,dist0,midist,mdist/n,madist,miang*anorm,mang/n*anorm,maang*anorm,midihed*anorm,mdihed/n*anorm,madihed*anorm,mdihed/(M_PI*2.0)-n/10.0);
		
}

void PointArray::sim_set_pot_parms(double pdist0,double pdistc,double pangc, double pdihed0, double pdihedc, double pmapc, EMData *pmap, double pmindistc,double pdistpenc) {
	dist0=pdist0;
	distc=pdistc;
	angc=pangc;
	dihed0=pdihed0;
	dihedc=pdihedc;
	mapc=pmapc;
	mindistc=pmindistc;
	distpenc=pdistpenc;
	if (pmap!=0 && pmap!=map) {
//		if (map!=0) delete map;
		if (gradx!=0) delete gradx;
		if (grady!=0) delete grady;
		if (gradz!=0) delete gradz;
		
		map=pmap;
		apix=map->get_attr("apix_x");
		if (map->get_zsize()==1)
			map2d=true;
		else
			map2d=false;
// 		centx=map->get_xsize()/2;
// 		centy=map->get_ysize()/2;
// 		centz=map->get_zsize()/2;
		centx=0;
		centy=0;
		centz=0;
// 		gradx=map->process("math.edge.xgradient");  // we compute the gradient to make the minimization easier
// 		grady=map->process("math.edge.ygradient");
// 		gradz=map->process("math.edge.zgradient");
	}
		
}

// Double the number of points by adding new points on the center of edges
void PointArray::sim_add_point_double() {
	
	int nn=n*2;
	int tmpn=n;
	set_number_points(nn);
	double* pa2data=(double *) calloc(4 * nn, sizeof(double));
	bool *newpt=new bool[nn];
	for (int i=0;i<nn;i++){
		if (i%2==0)
			newpt[i]=1;
		else
			newpt[i]=0;
	}
	int i=0;
	for (int ii=0;ii<nn;ii++){
		if (newpt[ii]) {
			
			pa2data[ii*4]=points[i*4];
			pa2data[ii*4+1]=points[i*4+1];
			pa2data[ii*4+2]=points[i*4+2];
			pa2data[ii*4+3]=1;
			i++;
		}
		else{
			int k;
			if (i<tmpn)
				k=i;
			else
				k=0;
			pa2data[ii*4]=(points[k*4]+points[(i-1)*4])/2;
			pa2data[ii*4+1]=(points[k*4+1]+points[(i-1)*4+1])/2;
			pa2data[ii*4+2]=(points[k*4+2]+points[(i-1)*4+2])/2;
			pa2data[ii*4+3]=1;
		}
			
	}
		
	delete []newpt;
	free(points);
	set_points_array(pa2data);
	
	if (adist) free(adist);
	if (aang) free(aang);
	if (adihed) free(adihed);
	adist=aang=adihed=0;
	sim_updategeom();
}

// Delete one point with lowest density, and add two points on the edges to that one.
// Or add one point on the edge with lowest density
void PointArray::sim_add_point_one() {

	
	double maxpot=-1000000,pot,meanpot=0;
	size_t ipt=0;
	bool onedge=0;
	// Find the highest potential point
	for (size_t i=0; i<n; i++) {
		meanpot+=sim_pointpotential(adist[i],aang[i],adihed[i]);
		pot=/*sim_pointpotential(adist[i],aang[i],adihed[i])*/-mapc*map->sget_value_at_interp(points[i*4]/apix+centx,points[i*4+1]/apix+centy,points[i*4+2]/apix+centz);
		if (pot>maxpot){
			maxpot=pot;
			ipt=i;
		}
	}
	meanpot/=n;
	
	for (size_t i=0; i<n; i++) {
		int k=(i+1)%n;
		double pt0,pt1,pt2;
		pt0=(points[k*4]+points[i*4])/2;
		pt1=(points[k*4+1]+points[i*4+1])/2;
		pt2=(points[k*4+2]+points[i*4+2])/2;
		pot=/*meanpot*/-mapc*map->sget_value_at_interp(pt0/apix+centx,pt1/apix+centy,pt2/apix+centz);
		if (pot>maxpot){
			maxpot=pot;
			ipt=i;
			onedge=1;
		}
	
	}
	
	// The rest points remain the same
	size_t i;
	double* pa2data=(double *) calloc(4 * (n+1), sizeof(double));
	for (size_t ii=0; ii<n+1; ii++) {
		if(ii!=ipt && ii!=ipt+1){
			if(ii<ipt)
				i=ii;
			else	// shift the points after the adding position
				i=ii-1;
			
			pa2data[ii*4]=points[i*4];
			pa2data[ii*4+1]=points[i*4+1];
			pa2data[ii*4+2]=points[i*4+2];
			pa2data[ii*4+3]=1;
		}
	}
	// Adding points
	if( onedge ) {
		size_t k1;
// 		k0=((ipt+n-1)%n);
		k1=((ipt+1)%n);
		pa2data[ipt*4]=points[ipt*4];
		pa2data[ipt*4+1]=points[ipt*4+1];
		pa2data[ipt*4+2]=points[ipt*4+2];
		pa2data[ipt*4+3]=1;
		
		pa2data[(ipt+1)*4]=(points[ipt*4]+points[k1*4])/2;
		pa2data[(ipt+1)*4+1]=(points[ipt*4+1]+points[k1*4+1])/2;
		pa2data[(ipt+1)*4+2]=(points[ipt*4+2]+points[k1*4+2])/2;
		pa2data[(ipt+1)*4+3]=1;	
		
	}
	else {
		size_t k0,k1;
		k0=((ipt+n-1)%n);
		k1=((ipt+1)%n);
		pa2data[ipt*4]=(points[ipt*4]+points[k0*4])/2;
		pa2data[ipt*4+1]=(points[ipt*4+1]+points[k0*4+1])/2;
		pa2data[ipt*4+2]=(points[ipt*4+2]+points[k0*4+2])/2;
		pa2data[ipt*4+3]=1;
		
		pa2data[(ipt+1)*4]=(points[ipt*4]+points[k1*4])/2;
		pa2data[(ipt+1)*4+1]=(points[ipt*4+1]+points[k1*4+1])/2;
		pa2data[(ipt+1)*4+2]=(points[ipt*4+2]+points[k1*4+2])/2;
		pa2data[(ipt+1)*4+3]=1;	
	
	}
	free(points);
	n++;
	set_points_array(pa2data);
	
	if (adist) free(adist);
	if (aang) free(aang);
	if (adihed) free(adihed);
	adist=aang=adihed=0;
	sim_updategeom();
	
	// search for the best position for the new points
	if (onedge){
		i=ipt+1;
		double bestpot=10000,nowpot;
		Vec3f old(points[i*4],points[(i+1)*4],points[(i+2)*4]);
		Vec3f newpt(points[i*4],points[(i+1)*4],points[(i+2)*4]);
		for (int ii=0;ii<5000;ii++){
			// Try to minimize potential globally
			boost::mt19937 rng; 
			boost::normal_distribution<> nd(0.0, 0.0);
			boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);
			points[i*4]=old[0]+var_nor();
			points[i*4+1]=old[1]+var_nor();
			points[i*4+2]=old[2]+var_nor();
			nowpot=sim_potentiald(i);
			if (nowpot<bestpot) {
				bestpot=nowpot;
				newpt[0]=points[i*4];
				newpt[1]=points[(i+1)*4];
				newpt[2]=points[(i+2)*4];
			}
				
		}
		points[i*4]=newpt[0];
		points[i*4+1]=newpt[1];
		points[i*4+2]=newpt[2];
	}
	
}

vector<float> PointArray::do_pca(int start=0, int end=-1){
	
	if (end==-1) end=n;
	float covmat[9],mean[3];
	for (int i=0; i<3; i++) mean[i]=0;
	for (int i=start; i<end; i++){
		for (int j=0; j<3; j++){
			mean[j]+=points[i*4+j];
		}
	}
	for (int i=0; i<3; i++) mean[i]/=end-start;
	
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			if (j<i){
				covmat[i*3+j]=covmat[j*3+i];
			}
			else{
				covmat[i*3+j]=0;
				for (int k=start; k<end; k++)
				{
					covmat[i*3+j]+=(points[k*4+i]-mean[i])*(points[k*4+j]-mean[j]);
				}
			}
			
// 			printf("%f\t",covmat[i*3+j]);
		}
// 		printf("\n");
	}
	
	float eigval[3],eigvec[9];
	Util::coveig(3,covmat,eigval,eigvec);
	vector<float> eigv(eigvec,eigvec+sizeof(eigvec)/sizeof(float));
// 	printf(" %f,%f,%f\n %f,%f,%f\n %f,%f,%f\n",eigvec[0],eigvec[1],eigvec[2],eigvec[3],eigvec[4],eigvec[5],eigvec[6],eigvec[7],eigvec[8]);
	return eigv;
}

vector<float> PointArray::do_filter(vector<float> pts, float *ft, int num){
	// filter a 1D array
	vector<float> result(pts);
	for (size_t i=0; i<pts.size(); i++)
		result[i]=0;
	for (size_t i=(num-1)/2; i<pts.size()-(num-1)/2; i++){
		for (int j=0; j<num; j++){
			int k=i+j-(num-1)/2;
			result[i]+=pts[k]*ft[j];
		}
	}
	return result;
	
}

vector<double> PointArray::fit_helix(EMData* pmap,int minlength=13,float mindensity=4, vector<int> edge=vector<int>(),int twodir=0,size_t minl=9)
{
	vector<float> hlxlen(n);
	vector<int> helix;
	map=pmap;
	float ft[7]={0.0044,0.0540,0.2420,0.3989,0.2420,0.0540,0.0044};
// 	float ft[7]={0,0,0,1,0,0,0};

	// search for long rods in the point array globally
	
	for (int dir=twodir; dir<2; dir++){ 
		// search in both directions and combine the result
		if( twodir==0)
		  reverse_chain();
		for (size_t i=0; i<n; i++){
			vector<float> dist(50);
			// for each point, search the following 50 points, find the longest rod
			for (int len=5; len<50; len++){
				size_t pos=i+len;
				if (pos>=n)	break;
				vector<float> eigvec=do_pca(i,pos); // align the points 
				vector<float> pts((len+1)*3);
				float mean[3];
				for (int k=0; k<3; k++) mean[k]=0;
				for (size_t k=i; k<pos; k++){
					for (int l=0; l<3; l++){
						pts[(k-i)*3+l]=points[k*4+0]*eigvec[l*3+0]+points[k*4+1]*eigvec[l*3+1]+points[k*4+2]*eigvec[l*3+2];
						mean[l]+=pts[(k-i)*3+l];
					}
				}
				for (int k=0; k<3; k++) mean[k]/=len;
				float dst=0;
				// distance to the center axis
				for (int k=0; k<len; k++){
					dst+=abs((pts[k*3]-mean[0])*(pts[k*3]-mean[0])+(pts[k*3+1]-mean[1])*(pts[k*3+1]-mean[1]));
				}
				dist[len]=1-dst/len/len;
			}
			
			vector<float> nd=do_filter(dist,ft,7);
			nd=do_filter(nd,ft,7);
			// length of the first rod
			for (int j=7; j<49; j++){
				if(nd[j]>nd[j-1] && nd[j]>nd[j+1]){
					hlxlen[i]=j-6;
					break;
				}
			}
			
			if(hlxlen[i]>25) hlxlen[i]=0;
		}
		// filter the array before further process
// 		hlxlen[50]=100;
// 		for (int i=0; i<n; i++) printf("%d\t%f\n",i,hlxlen[i]);
// 		for (int i=0; i<3; i++) hlxlen=do_filter(hlxlen,ft,7);
// 		for (int i=0; i<n; i++) printf("%d\t%f\n",i,hlxlen[i]);
		vector<float> ishlx(n);
		int hlx=0;
		float up=minlength; // rod length threshold
		// record position of possible helixes
		for (size_t i=0; i<n; i++){
			if(hlx<=0){
				if(hlxlen[i]>up){
					hlx=hlxlen[i];
					helix.push_back(i);
					helix.push_back(i+hlxlen[i]-5);
				}
			}
			else{
				hlx--;
				ishlx[i]=hlxlen[i];
			}
		}
		// while counting reversely
		if(dir==0){
			for (size_t i=0; i<helix.size(); i++) helix[i]=n-1-helix[i];
			for (size_t i=0; i<helix.size()/2; i++){
				int tmp=helix[i*2+1];
				helix[i*2+1]=helix[i*2];
				helix[i*2]=tmp;
			}
		}
			

	}

#ifdef DEBUG
	printf("potential helix counting from both sides: \n");
	for (size_t i=0; i<helix.size()/2; i++){
		printf("%d\t%d\n",helix[i*2],helix[i*2+1]);
	}	
	printf("\n\n");
#endif
/*

	// Combine the result from both side
	for (size_t i=0; i<helix.size()/2; i++){
		int change=1;
		while(change==1){
			change=0;
			for (size_t j=i+1; j<helix.size()/2; j++){
				if(helix[j*2]==0) continue;
				if(helix[j*2]-2<helix[i*2+1] && helix[j*2+1]+2>helix[i*2]){
					helix[i*2]=(helix[i*2]<helix[j*2])?helix[i*2]:helix[j*2];
					helix[i*2+1]=(helix[i*2+1]>helix[j*2+1])?helix[i*2+1]:helix[j*2+1];
					helix[j*2]=0;
					helix[j*2+1]=0;
					change=1;
				}
			}	
		}
	}*/
	
	vector<int> allhlx;
	int minid=1;
	while (minid>=0){
		int mins=10000;
		minid=-1;
		for (size_t i=0;i<helix.size()/2; i++){
			if(helix[i*2]<.1) continue;
			if(helix[i*2]<mins){
				mins=helix[i*2];
				minid=i;
			}
		}
		if(minid>=0){
			allhlx.push_back(helix[minid*2]);
			allhlx.push_back(helix[minid*2+1]);
			helix[minid*2]=-1;
		}		
	}
	
#ifdef DEBUG
	printf("combined result: \n");	
	for (size_t i=0; i<allhlx.size()/2; i++){
		printf("%d\t%d\n",allhlx[i*2],allhlx[i*2+1]);
	}	
	printf("\n\n");
#endif
	
	// local search to decide the start and end point of each helix
// 	vector<float> allscore(allhlx.size()/2);
	for (size_t i=0; i<allhlx.size()/2; i++){
		int sz=5;
		size_t start=allhlx[i*2]-sz,end=allhlx[i*2+1]+sz;
		start=start>0?start:0;
		end=end<n?end:n;
		float minscr=100000;
		int mj=0,mk=0;
		
		for (size_t j=start; j<end; j++){
			for (size_t k=j+6; k<end; k++){
				vector<float> eigvec=do_pca(j,k);
				vector<float> pts((k-j)*3);
				float mean[3];
				for (int u=0; u<3; u++) mean[u]=0;
				for (size_t u=j; u<k; u++){
					for (int v=0; v<3; v++){
						pts[(u-j)*3+v]=points[u*4+0]*eigvec[v*3+0]+points[u*4+1]*eigvec[v*3+1]+points[u*4+2]*eigvec[v*3+2];
						mean[v]+=pts[(u-j)*3+v];
					}
				}
				for (size_t u=0; u<3; u++) mean[u]/=(k-j);
				float dst=0;
				// distance to the center axis
				for (size_t u=0; u<k-j; u++){
					dst+=sqrt((pts[u*3]-mean[0])*(pts[u*3]-mean[0])+(pts[u*3+1]-mean[1])*(pts[u*3+1]-mean[1]));
				}
				float len=k-j;
				float scr=dst/len/len;
				if (scr<minscr){
// 					printf("%f\t%d\t%d\n",scr,j,k);
					minscr=scr;
					mj=j;
					mk=k;
				}
			}
		}

// 		printf("%d\t%d\n",mj,mk);
		
		allhlx[i*2]=mj;
		allhlx[i*2+1]=mk;        
// 		allscore[i]=minscr;
// 		if (mk-mj>60)
// 			allscore[i]=100;
	}
	
	for (size_t i=0; i<edge.size()/2; i++){
		allhlx.push_back(edge[i*2]);
		allhlx.push_back(edge[i*2+1]);
	}
	
	
	vector<int> allhlx2;
	minid=1;
	while (minid>=0){
		int mins=10000;
		minid=-1;
		for (size_t i=0;i<allhlx.size()/2; i++){
			if(allhlx[i*2]<.1) continue;
			if(allhlx[i*2]<mins){
				mins=allhlx[i*2];
				minid=i;
			}
		}
		if(minid>=0){
			allhlx2.push_back(allhlx[minid*2]<allhlx[minid*2+1]?allhlx[minid*2]:allhlx[minid*2+1]);
			allhlx2.push_back(allhlx[minid*2]>allhlx[minid*2+1]?allhlx[minid*2]:allhlx[minid*2+1]);
			allhlx[minid*2]=-1;
		}		
	}
	allhlx=allhlx2;
	
#ifdef DEBUG
	printf("Fitted helixes: \n");
	for (size_t i=0; i<allhlx.size()/2; i++){
		printf("%d\t%d\n",allhlx[i*2],allhlx[i*2+1]);
	}	
	printf("\n\n");
#endif
	// create ideal helix
	size_t ia=0,ka=0;
	int dir;
	vector<double> finalhlx;
	vector<double> hlxid;
	printf("Confirming helix... \n");
	while(ia<n){
		if (int(ia)==allhlx[ka*2]){
// 			int sz=(allhlx[ka*2+1]-allhlx[ka*2])>10?5:(allhlx[ka*2+1]-allhlx[ka*2])/2;
			int sz=3;
			float score=0,maxscr=0;
			float bestphs=0,phsscore=0,pscore=0;
			int mi=0,mj=0;
			for (int i=0; i<sz; i++){
				for (int j=0; j<sz; j++){
					int start=allhlx[ka*2]+i,end=allhlx[ka*2+1]-j;
					phsscore=0;bestphs=-1;
					for (float phs=-180; phs<180; phs+=10){ //search for phase
						construct_helix(start,end,phs,pscore,dir);
						if (pscore>phsscore){
							phsscore=pscore;
							bestphs=phs;
						}
					}
// 					printf("%f\t",bestphs);
					construct_helix(start,end,bestphs,score,dir);
					if (score>maxscr){
						maxscr=score;
						mi=i;
						mj=j;
					}
				}
			}
			int start=allhlx[ka*2]+mi,end=allhlx[ka*2+1]-mj;
			printf("%d\t%d\t%f\tden %d\t",start,end,maxscr,maxscr>mindensity);
			if (maxscr>mindensity){
				phsscore=0;
				for (float phs=-180; phs<180; phs+=10){ //search for phase
						construct_helix(start,end,phs,pscore,dir);
						if (pscore>phsscore){
							phsscore=pscore;
							bestphs=phs;
						}
				}
				vector<double> pts=construct_helix(start,end,bestphs,score,dir);
				int lendiff=end-start-pts.size()/3-2;
				printf("dir %d\t",dir);
				printf("len %lu\t diff %d\t",pts.size()/3-2,lendiff);
				if (pts.size()/3-2>minl && abs(lendiff)<15){
					
					for (int i=0; i<mi; i++){			
						finalhlx.push_back(points[(i+ia)*4]);
						finalhlx.push_back(points[(i+ia)*4+1]);
						finalhlx.push_back(points[(i+ia)*4+2]);
					}
					hlxid.push_back(finalhlx.size()/3+1);
					printf("%lu\t",finalhlx.size()/3+1);
					for (size_t j=3; j<pts.size()-3; j++)
						finalhlx.push_back(pts[j]);
					hlxid.push_back(finalhlx.size()/3-2);
					printf("%lu\t",finalhlx.size()/3-2);
					for (size_t j=0; j<3; j++)
						hlxid.push_back(pts[j]);
					for (size_t j=pts.size()-3; j<pts.size(); j++)
						hlxid.push_back(pts[j]);
					ia=end;
				}
				else{
					printf("\t *");
					
				}
				
			}
			printf("\n\n");
			ka++;
			while(allhlx[ka*2]<int(ia))
				ka++;
		}
		else{
// 			printf("%d\t",ia);
			finalhlx.push_back(points[ia*4]);
			finalhlx.push_back(points[ia*4+1]);
			finalhlx.push_back(points[ia*4+2]);
			ia++;			
		}
	}
	
	set_number_points(finalhlx.size()/3);
	
	for (size_t i=0; i<n; i++){
		for (size_t j=0; j<3; j++)
			points[i*4+j]=finalhlx[i*3+j];
		points[i*4+3]=0;
	}
			

	
	printf("\n\n");
	return hlxid;
}

vector<double> PointArray::construct_helix(int start,int end, float phs, float &score, int &rtdir){
	// calculate length
	int dir=1;
	rtdir=0;
	Vec3f d(points[end*4]-points[start*4],points[end*4+1]-points[start*4+1],points[end*4+2]-points[start*4+2]);
	double len=d.length();
	int nh=int(Util::round(len/1.54))+2;
	vector<double> helix(nh*3);
	vector<float> eigvec=do_pca(start,end);	
	float eigval[3],vec[9];

	Util::coveig(3,&eigvec[0],eigval,vec);
	float maxeigv=0;
// 	int maxvi=-1;
	for(int i=0; i<3; i++){
		if(abs(eigval[i])>maxeigv){
			maxeigv=abs(eigval[i]);
// 			maxvi=i;
		}
	}
// 	dir=eigval[maxvi]>0;
	dir=1;
	vector<double> pts(nh*3);
	for (int dd=0; dd<2; dd++){
	// 	printf("%f\t",eigval[maxvi]);
	// 	vector<float> eigv(eigvec,eigvec+sizeof(eigvec)/sizeof(float));
		// create helix
		helix[0]=0;helix[1]=0;helix[2]=0;
		helix[nh*3-3]=.0;helix[nh*3-2]=0;helix[nh*3-1]=len;
		
		for (int i=0; i<nh-2; i++){
			if(dir>0){
				helix[(i+1)*3+0]=cos(((phs+(100*i))*M_PI)/180)*2.3;
				helix[(i+1)*3+1]=sin(((phs+(100*i))*M_PI)/180)*2.3;
			}
			else{
				helix[(i+1)*3+1]=cos(((phs+(100*i))*M_PI)/180)*2.3;
				helix[(i+1)*3+0]=sin(((phs+(100*i))*M_PI)/180)*2.3;
			}	
			helix[(i+1)*3+2]=i*1.54;
		}
		// transform to correct position
		float mean[3];
		for (int k=0; k<3; k++) mean[k]=0;
		for (int k=0; k<nh; k++){
			for (int l=0; l<3; l++){
				pts[k*3+l]=helix[k*3+0]*eigvec[0*3+l]+helix[k*3+1]*eigvec[1*3+l]+helix[k*3+2]*eigvec[2*3+l];
				mean[l]+=pts[k*3+l];
			}
		}
		for (int k=0; k<3; k++) mean[k]/=nh;
		for (int k=0; k<nh; k++){
			for (int l=0; l<3; l++){
				pts[k*3+l]-=mean[l];
			}
		}
		for (int k=0; k<3; k++) mean[k]=0;
		for (int k=start; k<end; k++){
			for (int l=0; l<3; l++){
				mean[l]+=points[k*4+l];
			}
		}
		for (int k=0; k<3; k++) mean[k]/=(end-start);	
		for (int k=0; k<nh; k++){
			for (int l=0; l<3; l++){
				pts[k*3+l]+=mean[l];
			}
		}
		
		
		// correct direction
		Vec3f d1(pts[0]-points[start*4],pts[1]-points[start*4+1],pts[2]-points[start*4+2]);
		Vec3f d2(pts[0]-points[end*4],pts[1]-points[end*4+1],pts[2]-points[end*4+2]);
		
		if (d1.length()>d2.length()) { //do reverse
			double tmp;
			for (int i=0; i<nh/2; i++){
				for(int j=0; j<3; j++){
					tmp=pts[i*3+j];
					pts[i*3+j]=pts[(nh-i-1)*3+j];
					pts[(nh-i-1)*3+j]=tmp;
				}
			}
		}
		
		// correct handness
		
		
		float hel=calc_helicity(pts);
		
		rtdir=hel;
// 		break;
		if(hel>0)
			break;
		else
			dir=0;
	

	}
	// calculate score
// 	int sx=map->get_xsize(),sy=map->get_ysize(),sz=map->get_zsize();
	int sx=0,sy=0,sz=0;
	float ax=map->get_attr("apix_x"),ay=map->get_attr("apix_y"),az=map->get_attr("apix_z");
	score=0;
	float aind=0;
	for (int i=1; i<nh-1; i++){
		float ind=1;//(float(abs(i-nh/2))/float(nh))+.5;
		score+=ind*map->get_value_at(int(pts[i*3]/ax+sx/2),int(pts[i*3+1]/ay+sy/2),int(pts[i*3+2]/az+sz/2));
		aind+=ind;
	}
	score/=aind;//(nh-2);
// 	float nsc=score;
// 	score=0;
// 	for (int i=1; i<nh-1; i++){
// 		float ind=(float(abs(i-nh/2))/float(nh))+.5;
// 		float mapden=map->get_value_at(int(pts[i*3]/ax+sx/2),int(pts[i*3+1]/ay+sy/2),int(pts[i*3+2]/az+sz/2));
// 		if (mapden>nsc*1.2)
// 			mapden=nsc*1.2;
// 		if (mapden<nsc*.5)
// 			mapden=0;
// 		score+=ind*mapden;
// 	}
// 	score/=aind;//(nh-2);
	
		
	return pts;
}

void PointArray::merge_to(PointArray &pa,float thr=3.5)
{
	printf("merging\n");
	vector<double> result;
	for (size_t i=0; i<pa.n; i++){
		result.push_back(pa.points[i*4]);
		result.push_back(pa.points[i*4+1]);
		result.push_back(pa.points[i*4+2]);
	}
	for (size_t i=0; i<n; i++){
		float dist=100;
		for (size_t j=0; j<pa.n; j++){
			Vec3f d(points[i*4]-pa.points[j*4],points[i*4+1]-pa.points[j*4+1],points[i*4+2]-pa.points[j*4+2]);
			dist=d.length();
			if (dist<thr)
				break;
		}
		printf("%f\n",dist);
		if (dist>thr){
			result.push_back(points[i*4]);
			result.push_back(points[i*4+1]);
			result.push_back(points[i*4+2]);
		}
	}
	set_number_points(result.size()/3);
	for (size_t i=0; i<n; i++){
		for (size_t j=0; j<3; j++)
			points[i*4+j]=result[i*3+j];
		points[i*4+3]=0;
	}
	printf("done\n");
	
}

float PointArray::calc_helicity(vector<double> pts){
	int npt=pts.size();
	vector<double> vpoint(npt-3);
	for (int i=0; i<npt-3; i++){
		vpoint[i]=pts[i+3]-pts[i];		
	}
	Vec3f hlxdir(pts[npt-3]-pts[0],pts[npt-2]-pts[1],pts[npt-1]-pts[2]);
	Vec3f vcs(0,0,0);
	for (int i=0; i<npt/3-2; i++){
		Vec3f v0(vpoint[i*3],vpoint[i*3+1],vpoint[i*3+2]);
		Vec3f v1(vpoint[i*3+3],vpoint[i*3+4],vpoint[i*3+5]);
		vcs+=v0.cross(v1);
	}
	vcs/=(npt/3-2);

	return hlxdir.dot(vcs);
  
}
void PointArray::save_pdb_with_helix(const char *file, vector<float> hlxid)
{

	FILE *fp = fopen(file, "w");
	
	for (size_t i=0; i<hlxid.size()/8; i++){
		fprintf(fp, "HELIX%5lu   A ALA A%5d  ALA A%5d  1                        %5d\n", 
				i, (int)hlxid[i*8], (int)hlxid[i*8+1], int(hlxid[i*8+1]-hlxid[i*8]+4));
	}
	for ( size_t i = 0; i < get_number_points(); i++) {
		fprintf(fp, "ATOM  %5lu  CA  ALA A%4lu    %8.3f%8.3f%8.3f%6.2f%6.2f%8s\n", i, i,
				points[4 * i], points[4 * i + 1], points[4 * i + 2], points[4 * i + 3], 0.0, " ");
	}
	fprintf(fp, "TER   %5lu      ALA A%4lu\nEND", get_number_points(),get_number_points()-1);

	fclose(fp);
}

void PointArray::remove_helix_from_map(EMData *m, vector<float> hlxid){
	
	int sx=m->get_xsize(),sy=m->get_ysize(),sz=m->get_zsize();
	float ax=m->get_attr("apix_x"),ay=m->get_attr("apix_y"),az=m->get_attr("apix_z");
	for (int x=0; x<sx; x++){
		for (int y=0; y<sy; y++){
			for (int z=0; z<sz; z++){
				Vec3f p0((x)*ax,(y)*ay,(z)*az);
				bool inhlx=false;
				for (size_t i=0; i<hlxid.size()/8; i++){
					Vec3f p1(hlxid[i*8+2],hlxid[i*8+3],hlxid[i*8+4]),p2(hlxid[i*8+5],hlxid[i*8+6],hlxid[i*8+7]);
					Vec3f dp=p2-p1;
					float l=dp.length();
					float d=((p0-p1).cross(p0-p2)).length()/l;
					float t=-(p1-p0).dot(p2-p1)/(l*l);
					if (d<5 && t>0 && t<1){
						inhlx=true;
						break;
					}
				}
				if(inhlx){
					m->set_value_at(x,y,z,0);
				}
					
			}
		}
	}
}


Transform PointArray::calc_transform(PointArray *p)
{
	
	Transform t;
	
	// calculate translation
	FloatPoint c1=p->get_center(),c2=get_center();
	t.set_trans(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]);
	
	// calculate rotation
	// this rotation rotates unit vectors a,b into A,B;
	//    The program assumes a dot b must equal A dot B
	
	//pick two bonds randomly
	int i1=Util::get_irand(0,n-2);
	int i2=Util::get_irand(0,n-2);
	while(i1==i2)
		i2=Util::get_irand(0,n-2);
	
	const Vec3f eahat=get_vector_at(i1+1)-get_vector_at(i1);
	const Vec3f ebhat=get_vector_at(i2+1)-get_vector_at(i2);
        const Vec3f eAhat=p->get_vector_at(i1+1)-p->get_vector_at(i1);
        const Vec3f eBhat=p->get_vector_at(i2+1)-p->get_vector_at(i2);
	Vec3f eahatcp(eahat);
	Vec3f ebhatcp(ebhat);
	Vec3f eAhatcp(eAhat);
	Vec3f eBhatcp(eBhat);

	eahatcp.normalize();
	ebhatcp.normalize();
	eAhatcp.normalize();
	eBhatcp.normalize();

	Vec3f aMinusA(eahatcp);
	aMinusA  -= eAhatcp;
	Vec3f bMinusB(ebhatcp);
	bMinusB  -= eBhatcp;


	Vec3f  nhat;
	float aAlength = aMinusA.length();
	float bBlength = bMinusB.length();
	if (aAlength==0){
		nhat=eahatcp;
	}else if (bBlength==0){
		nhat=ebhatcp;
	}else{
		nhat= aMinusA.cross(bMinusB);
		nhat.normalize();
	}

//		printf("nhat=%f,%f,%f \n",nhat[0],nhat[1],nhat[2]);

	Vec3f neahat  = eahatcp.cross(nhat);
	Vec3f nebhat  = ebhatcp.cross(nhat);
	Vec3f neAhat  = eAhatcp.cross(nhat);
	Vec3f neBhat  = eBhatcp.cross(nhat);

	double cosomegaA = (neahat.dot(neAhat))  / (neahat.dot(neahat));
//	float cosomegaB = (nebhat.dot(neBhat))  / (nebhat.dot(nebhat));
	double sinomegaA = (neahat.dot(eAhatcp)) / (neahat.dot(neahat));
//	printf("cosomegaA=%f \n",cosomegaA); 	printf("sinomegaA=%f \n",sinomegaA);

	double omegaA = atan2(sinomegaA,cosomegaA);
//	printf("omegaA=%f \n",omegaA*180/M_PI);
	Dict rotation;
	rotation["type"]="spin";
	rotation["n1"]= nhat[0];
	rotation["n2"]= nhat[1];
	rotation["n3"]= nhat[2];
	rotation["omega"] = (float)(omegaA*180.0/M_PI);
	t.set_rotation(rotation);
	return t;
}

void PointArray::reverse_chain(){
	// reverse the point array chain, from the last to the first point
	double tmp;
// 	for(int i=0; i<n/2; i++){
// 		for (int j=0; j<4; j++){
// 			printf("%f\t",
	for(size_t i=0; i<n/2; i++){
		for (size_t j=0; j<4; j++){
			tmp=points[(n-1-i)*4+j];
			points[(n-1-i)*4+j]=points[i*4+j];
			points[i*4+j]=tmp;
		}
	}
}

void PointArray::delete_point(int id){
	for (size_t i=id*4; i<n*4; i++){
		points[i]=points[i+4];
	}
	set_number_points(n-1);
	
}


bool PointArray::read_ca_from_pdb(const char *file)
{
	struct stat filestat;
	stat(file, &filestat);
	set_number_points(( int)(filestat.st_size / 80 + 1));

	#ifdef DEBUG
	printf("PointArray::read_from_pdb(): try %4lu atoms first\n", get_number_points());
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

		float e = 6;
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
		if (ctt != 'C' || ctt2 != 'A')
			continue;
		
		float x, y, z, q;
		sscanf(&s[28], " %f %f %f", &x, &y, &z);
		sscanf(&s[60], " %f", &q);

		if (count + 1 > get_number_points())
			set_number_points(2 * (count + 1));    //makes sure point array is big enough
		
#ifdef DEBUG
		printf("Atom %4lu: x,y,z = %8g,%8g,%8g\te = %g\n", count, x, y, z, e);
#endif
		points[4 * count] = x;
		points[4 * count + 1] = y;
		points[4 * count + 2] = z;
		points[4 * count + 3] = e;
		bfactor[count] = q;
		count++;
	}
	fclose(fp);
	set_number_points(count);
	return true;
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


EMData *PointArray::pdb2mrc_by_summation(int map_size, float apix, float res, int addpdbbfactor)
{
#ifdef DEBUG
	printf("PointArray::pdb2mrc_by_summation(): %lu points\tmapsize = %4d\tapix = %g\tres = %g\n",get_number_points(),map_size, apix, res);
#endif
	//if ( gauss_real_width < apix) LOGERR("PointArray::projection_by_summation(): apix(%g) is too large for resolution (%g Angstrom in Fourier space) with %g pixels of 1/e half width", apix, res, gauss_real_width);

	double min_table_val = 1e-7;
	double max_table_x = sqrt(-log(min_table_val));	// for exp(-x*x)

	double table_step_size = 0.001;	// number of steps for each pixel
	double inv_table_step_size = 1.0 / table_step_size;
	
//	sort_by_axis(2);			// sort by Z-axis

	EMData *map = new EMData();
	map->set_size(map_size, map_size, map_size);
	map->to_zero();
	float *pd = map->get_data();
	
	vector<double> table;
	double gauss_real_width;
	int table_size;
	int gbox;
	
	
	for ( size_t s = 0; s < get_number_points(); ++s) {
		double xc = points[4 * s] / apix + map_size / 2;
		double yc = points[4 * s + 1] / apix + map_size / 2;
		double zc = points[4 * s + 2] / apix + map_size / 2;
		double fval = points[4 * s + 3];
		
		//printf("\n bfactor=%f",bfactor[s]);
		
		
		
		
		
		if(addpdbbfactor==-1){
			gauss_real_width = res/M_PI;	// in Angstrom, res is in Angstrom
		}else{
			gauss_real_width = (bfactor[s])/(4*sqrt(2.0)*M_PI);	// in Angstrom, res is in Angstrom
		}
		
		
		table_size = int (max_table_x * gauss_real_width / (apix * table_step_size) * 1.25);
		table.resize(table_size);
		for (int i = 0; i < table_size; i++){
			table[i] = 0;
		}
		
		for (int i = 0; i < table_size; i++) {						
			double x = -i * table_step_size * apix / gauss_real_width;
			if(addpdbbfactor==-1){
				table[i] =  exp(-x * x);	
			}
			else{
			table[i] =  exp(-x * x)/sqrt(gauss_real_width * gauss_real_width * 2 * M_PI);
			}
		}

		gbox = int (max_table_x * gauss_real_width / apix);	// local box half size in pixels to consider for each point
		if (gbox <= 0)
			gbox = 1;
		
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
			size_t table_index_z = size_t (fabs(k - zc) * inv_table_step_size);
			if ( table_index_z >= table.size() ) continue;
			double zval = table[table_index_z];
			size_t pd_index_z = k * map_size * map_size;
			
			for (int j = jmin; j < jmax; j++) {
				size_t table_index_y = size_t (fabs(j - yc) * inv_table_step_size);
				if ( table_index_y >= table.size() ) continue;
				double yval = table[table_index_y];
				size_t pd_index = pd_index_z + j * map_size + imin;
				for (int i = imin; i < imax; i++, pd_index++) {
					size_t table_index_x = size_t (fabs(i - xc) * inv_table_step_size);
					if ( table_index_x >= table.size() ) continue;
					double xval = table[table_index_x];
					pd[pd_index] += (float) (fval * zval * yval * xval);
				}
			}
		}
	}

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
#if defined NFFT2
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
	//fft->process_inplace(".filter.lowpass.gauss",Dict("cutoff_abs", map_size*apix/res));

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
#if defined NFFT2
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
	//fft->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs", box*apix/res));

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
	Transform xform;
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
