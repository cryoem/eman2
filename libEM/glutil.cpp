/*
 * Author: David Woolford, 11/06/2007 (woolford@bcm.edu)
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

#ifdef USE_OPENGL

#ifdef _WIN32
	#include <windows.h>
#endif

#ifndef GL_GLEXT_PROTOTYPES
	#define GL_GLEXT_PROTOTYPES
#endif	//GL_GLEXT_PROTOTYPES

#include "glutil.h"
#include "emdata.h"
#include "marchingcubes.h"

#ifdef __APPLE__
	#include "OpenGL/gl.h"
	#include "OpenGL/glu.h"
	#include "OpenGL/glext.h"
#else // WIN32, LINUX
	#include "GL/gl.h"
	#include "GL/glu.h"
	#include "GL/glext.h"
#endif	//__APPLE__

using namespace EMAN;

// By depfault we need to first bind data to the GPU

GLuint GLUtil::buffer[2] = {0, 0};

unsigned int GLUtil::gen_glu_mipmaps(const EMData* const emdata)
{
	if (emdata->get_data() == 0) {
		throw NullPointerException("Error, attempt to create an OpenGL mipmap "
			"without internally stored data");
	}

	ENTERFUNC;

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	if (emdata->ny == 1 && emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_1D, tex_name);
		gluBuild1DMipmaps(GL_TEXTURE_1D, GL_LUMINANCE, emdata->nx, GL_LUMINANCE,
			 GL_FLOAT, (void*)(emdata->get_data()));
	} else if (emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_2D, tex_name);
		gluBuild2DMipmaps(GL_TEXTURE_2D, GL_LUMINANCE, emdata->nx, emdata->ny,
			 GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
	}
	else {
#ifdef	_WIN32
		// There is no gluBuild3DMipmaps() function in glu.h for VS2003 and VS2005

		printf("3D OpenGL mipmaps are not available on this platform.\n");
#else
		glBindTexture(GL_TEXTURE_3D, tex_name);
		gluBuild3DMipmaps(GL_TEXTURE_3D, GL_LUMINANCE, emdata->nx, emdata->ny,
			 emdata->nz, GL_LUMINANCE, GL_FLOAT, (void*)(emdata->get_data()));
#endif	//_WIN32
	}

	EXITFUNC;

	return tex_name;
}

unsigned int GLUtil::gen_gl_texture(const EMData* const emdata, GLenum format)
{
	if (emdata->get_data() == 0) {
		throw NullPointerException("Error, attempt to create an OpenGL texture"
			" without internally stored data");
	}

	ENTERFUNC;

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	if (emdata->ny == 1 && emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_1D, tex_name);
		glTexImage1D(GL_TEXTURE_1D, 0, format, emdata->nx, 0, format, GL_FLOAT,
			 (void*)(emdata->get_data()));
	} else if (emdata->nz == 1) {
		glBindTexture(GL_TEXTURE_2D,  tex_name);
		glTexImage2D(GL_TEXTURE_2D, 0, format, emdata->nx, emdata->ny, 0, format,
			 GL_FLOAT, (void*)(emdata->get_data()));
	}
	else {
		glBindTexture(GL_TEXTURE_3D, tex_name);
#ifdef _WIN32
	PFNGLTEXIMAGE3DPROC glTexImage3D;
#endif
		glTexImage3D(GL_TEXTURE_3D, 0, format,
		emdata->nx, emdata->ny, emdata->nz, 0, format, GL_FLOAT,
		(void*)(emdata->get_data()));
	}

	EXITFUNC;

	return tex_name;
}

unsigned int GLUtil::render_amp8_gl_texture(EMData* emdata,
		 int x0, int y0, int ixsize, int iysize, int bpl, float scale,
		 int mingray, int maxgray, float render_min, float render_max,
		 float gamma, int flags)
{
	if (emdata==NULL) return 9999999;
	string pixels = render_amp8(emdata, x0, y0, ixsize,iysize, bpl, scale,
		 mingray, maxgray, render_min, render_max, gamma, flags);

	unsigned int tex_name;
	glGenTextures(1, &tex_name);

	glBindTexture(GL_TEXTURE_2D, tex_name);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, ixsize, iysize , 0,
		 GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels.c_str());

	return tex_name;
}

// undef GL_GLEXT_PROTOTYPES

#ifdef GL_GLEXT_PROTOTYPES
#undef GL_GLEXT_PROTOTYPES
#endif

int GLUtil::nearest_projected_points(const vector<float>& model_matrix,
		 const vector<float>& proj_matrix, const vector<int>& view_matrix,
		 const vector<Vec3f>& points, const float mouse_x, const float mouse_y,
		 const float& nearness)
{
	double proj[16];
	double model[16];
	int view[4];

	copy(proj_matrix.begin(), proj_matrix.end(), proj);
	copy(model_matrix.begin(), model_matrix.end(), model);
	copy(view_matrix.begin(), view_matrix.end(), view);

	vector<Vec3f> unproj_points;
	double x,y,z;
	double r,s,t;

	for (vector<Vec3f>::const_iterator it = points.begin(); it != points.end(); ++it) {
		r = (double) (*it)[0];
		s = (double) (*it)[1];
		t = (double) (*it)[2];

		gluProject(r,s,t, model, proj, view, &x, &y, &z);
		unproj_points.push_back(Vec3f(x, y, z));
	}

	vector<int> intersections;

	float n_squared = nearness * nearness;

	for(unsigned int i = 0; i < unproj_points.size(); ++i) {
		Vec3f& v = unproj_points[i];
		float dx = v[0] - mouse_x;
		dx *= dx;
		float dy = v[1] - mouse_y;
		dy *= dy;

		if ((dx+dy) <= n_squared) intersections.push_back((int)i);
	}

	int intersection = -1;
	float near_z = 0;

	for(vector<int>::const_iterator it = intersections.begin(); it != intersections.end(); ++it) {
		if (intersection == -1 || unproj_points[*it][2] < near_z) {
			intersection = *it;
			near_z = unproj_points[*it][2];
		}
	}

	return intersection;
}

void GLUtil::colored_rectangle(const vector<float>& data,
		 const float& alpha, const bool center_point)
{
	glBegin(GL_LINE_LOOP);
	glColor4f(data[0], data[1], data[2], alpha);
	glVertex2i(data[3], data[4]);
	glVertex2i(data[5], data[4]);
	glVertex2i(data[5], data[6]);
	glVertex2i(data[3], data[6]);
	glEnd();

	if (center_point) {
		glBegin(GL_POINTS);
		glVertex2f( (data[3]+data[5])/2, (data[4]+data[6])/2);
		glEnd();
	}
}

void GLUtil::mx_bbox(const vector<float>& data,
		 const vector<float>& text_color, const vector<float>& bg_color)
{
	glDisable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);

	if (bg_color.size() == 4) {
		glColor4f(bg_color[0], bg_color[1], bg_color[2], bg_color[3]);
	}
	else {
		glColor3f(bg_color[0], bg_color[1], bg_color[2]);
	}

	glVertex3f(data[0]-1, data[1]-1, -.1f);
	glVertex3f(data[3]+1, data[1]-1, -.1f);
	glVertex3f(data[3]+1, data[4]+1, -.1f);
	glVertex3f(data[0]-1, data[4]+1, -.1f);
	glEnd();

//	glEnable(GL_TEXTURE_2D);

	if (text_color.size() == 4) {
		glColor4f(text_color[0], text_color[1], text_color[2], text_color[3]);
	}
	else {
		glColor3f(text_color[0], text_color[1], text_color[2]);
	}
}

std::string GLUtil::render_amp8(EMData* emdata, int x0, int y0, int ixsize,
		 int iysize, int bpl, float scale, int min_gray, int max_gray,
		 float render_min, float render_max, float gamma, int flags)
{
	ENTERFUNC;

//	printf("%f\t%f\t",(float)emdata->get_attr("sigma"),(float)emdata->get_attr("mean"));
//	printf("%d %d %d %d %d %f %d %d %f %f %f %d\n",x0,y0,ixsize,iysize,bpl,
// scale,min_gray,max_gray,render_min,render_max,gamma,flags);

	if (emdata==NULL) return std::string();
	bool invert = (min_gray > max_gray);
	int mingray, maxgray;

	if (invert) {
		mingray = max_gray;
		maxgray = min_gray;
	}
	else {
		mingray = min_gray;
		maxgray = max_gray;
	}

	int asrgb;
	int hist = (flags & 2)/2;
	int invy = (flags & 4)?1:0;

	int nx = emdata->nx;
	int ny = emdata->ny;
//	int nz = emdata->nz;
	int nxy = emdata->nx * emdata->ny;

	if (emdata->get_ndim() > 2) {
		throw ImageDimensionException("1D/2D only");
	}

	if (emdata->is_complex()) {
		emdata->ri2ap();
	}

	if (render_max <= render_min) {
		render_max = render_min + 0.01f;
	}

	if (gamma <= 0) gamma = 1.0;

	// Calculating a full floating point gamma for
	// each piGLUtil::xel in the image slows rendering unacceptably
	// however, applying a gamma-mapping to an 8 bit colorspace
	// has unacceptable coarse accuracy. So, we oversample the 8 bit colorspace
	// as a 12 bit colorspace and apply the gamma mapping to that
	// This should produce good accuracy for gamma values
	// larger than 0.5 (and a high upper limit)

	static int smg0 = 0, smg1 = 0; // while this destroys threadsafety in the rendering process
	static float sgam=0; // it is necessary for speed when rendering large numbers of small images
	static unsigned char gammamap[4096];

	if (gamma != 1.0 && (smg0 != mingray || smg1 != maxgray || sgam != gamma)) {
		for (int i=0; i<4096; i++) {
			if (mingray<maxgray) {
				gammamap[i] = (unsigned char)(mingray+(maxgray-mingray+0.999)*pow(((float)i/4096.0f),gamma));
			}
			else {
				gammamap[4095-i] = (unsigned char)(mingray+(maxgray-mingray+0.999)*pow(((float)i/4096.0f),gamma));
			}
		}
	}
	smg0 = mingray; // so we don't recompute the map unless something changes
	smg1 = maxgray;
	sgam = gamma;

	if (flags & 8) asrgb = 4;
	else if (flags & 1) asrgb = 3;
	else asrgb = 1;

	std::string ret=std::string();
//	ret.resize(iysize*bpl);

	ret.assign(iysize*bpl + hist*1024, char(invert ? maxgray : mingray));

	unsigned char *data = (unsigned char *)ret.data();
	unsigned int *histd = (unsigned int *)(data + iysize*bpl);

	if (hist) {
		for (int i=0; i<256; i++) histd[i]=0;
	}

	float rm = render_min;
	float inv_scale = 1.0f / scale;
	int ysize = iysize;
	int xsize = ixsize;

	int ymin = 0;

	if (iysize * inv_scale > ny) {
		ymin = (int) (iysize - ny / inv_scale);
	}

	float gs = (maxgray - mingray) / (render_max - render_min);
	float gs2 = 4095.999f / (render_max - render_min);
//	float gs2 = 1.0 / (render_max - render_min);

	if (render_max < render_min) {
		gs = 0;
		rm = FLT_MAX;
	}

	int dsx = -1;
	int dsy = 0;
	int remx = 0;
	int remy = 0;
	const int scale_n = 100000;

	int addi = 0;
	int addr = 0;

	if (inv_scale == floor(inv_scale)) {
		dsx = (int) inv_scale;
		dsy = (int) (inv_scale * nx);
	}
	else {
		addi = (int) floor(inv_scale);
		addr = (int) (scale_n * (inv_scale - floor(inv_scale)));
	}

	int xmin = 0;

	if (x0 < 0) {
		xmin = (int) (-x0 / inv_scale);
		xsize -= (int) floor(x0 / inv_scale);
		x0 = 0;
	}

	if ((xsize - xmin) * inv_scale > (nx - x0)) {
		xsize = (int) ((nx - x0) / inv_scale + xmin);
	}

	int ymax = ysize - 1;

	if (y0 < 0) {
		ymax = (int) (ysize + y0 / inv_scale - 1);
		ymin += (int) floor(y0 / inv_scale);
		y0 = 0;
	}

	if (xmin < 0) xmin = 0;
	if (ymin < 0) ymin = 0;
	if (xsize > ixsize) xsize = ixsize;
	if (ymax > iysize) ymax = iysize;

	int lmax = nx * ny - 1;

	int mid=nx*ny/2;
	float * image_data = emdata->get_data();

	////// Begin of Histogram Equalization //////

	const int rangemax = 4082;
	int gaussianwide = 5000;
	unsigned int* grayhe = NULL;
	float* gaussianpdf   = NULL;
	float* gaussiancdf   = NULL;
	int* gaussianlookup  = NULL;

	unsigned int* graypdftwo = NULL;

	bool binflag = 0;
	int binlocation = -1;

	if (flags & 32) {
		int graypdf[rangemax]   = {0}; // 256
		int graycdf[rangemax-2] = {0}; // 254

		// unsigned int grayhe[rangemax-2]={0};//#number=254

		graypdftwo = new unsigned int[maxgray-mingray];
		grayhe = new unsigned int[rangemax-2];
		// unsigned char graylookup[(int)(render_max-render_min)];//render_max-render_min

		for (int i=0; i<(int)(rangemax-2); i++) {
			grayhe[i] = 0; // Initialize all elements to zero.
		}

		for (int i=0; i<(maxgray-mingray); i++) { // 0~253
			graypdftwo[i] = 0;
		}

		if (dsx != -1) {
			int l = x0 + y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int br = l;
				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) {
						break;
					}

					int k = 0;
					int p;
					float t;

					if (dsx == 1) t=image_data[l];
					else { // This block does local pixel averaging for nicer reduced views
						t=0;

						if ((l+dsx+dsy) > lmax) {
							break;
						}

						for (int iii=0; iii<dsx; iii++) {
							for (int jjj=0; jjj<dsy; jjj+=nx) {
								t += image_data[l+iii+jjj];
							}
						}

						t /= dsx*(dsy/nx);
					}

					if (t <= rm) graypdf[0]++;
					else if (t >= render_max) graypdf[rangemax-1]++;
					else {
						graypdf[(int)(ceil((rangemax-2)*(t - render_min)/(render_max-render_min)))]++;
						graypdftwo[(unsigned char) (gs * (t - render_min))]++;
					}

					l += dsx;
				}

				l = br + dsy;
			}
		}
		else {
			remy = 10;
			int l = x0 + y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int addj = addi;

				// There seems to be some overflow issue happening
				// where the statement if (l > lmax) break (below) doesn't work
				// because the loop that iterates jjj can inadvertantly go out of bounds

				if ((l + addi*nx) >= nxy) {
					addj = (nxy-l)/nx;

					if (addj <= 0) continue;
				}

				int br = l;
				remx = 10;

				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) break;
					int p = 0;
					float t;
					if (addi <= 1) t = image_data[l];
					else { // This block does local pixel averaging for nicer reduced views
						t = 0;

						for (int jjj=0; jjj<addj; jjj++) {
							for (int iii=0; iii<addi; iii++) {
								t += image_data[l+iii+jjj*nx];
							}
						}

						t /= addi*addi;
					}

					////////////

					if (t <= rm) graypdf[0]++;
					else if (t >= render_max) graypdf[rangemax-1]++;
					else {
						graypdf[(int)(ceil((rangemax-2)*(t - render_min)/(render_max-render_min)))]++;
					}

//					data[i * asrgb + j * bpl] = p;

//					if (hist) histd[p]++;

					l += addi;
					remx += addr;

					if (remx > scale_n) {
						remx -= scale_n;
						l++;
					}
				}

				l = br + addi * nx;
				remy += addr;

				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}

		for (int i=0; i<(rangemax-2); i++) { // 0~253
			for (int j=0;j<(i+1);j++) {
				graycdf[i]=graycdf[i]+graypdf[j+1];
			}
		}

		// graypdftwo binflag

		binflag = 0;

		for (int i=1; i<(maxgray-mingray); i++) { // 0~253
			if (((float)graypdftwo[i]/graycdf[rangemax-3])>0.2) {
				binflag = 1;
				binlocation = i;

				break;
			}
		}

		if (binflag == 1) {
			for (int i=(binlocation*16+1); i<((binlocation+1)*16); i++) {
				graypdf[i] = 0;
				// graypdf[i]=(graycdf[rangemax-3]-graypdftwo[binlocation])/(maxgray-mingray);
			}

			for (int i=0; i<(rangemax-2); i++) { // 0~253
				graycdf[i] = 0;
			}

			for (int i=0; i<(rangemax-2); i++) { // 0~253
				for (int j=0;j<(i+1);j++) {
					graycdf[i] = graycdf[i]+graypdf[j+1];
				}
			}
		}

		// start gaussian matching

		float mean = abs(rangemax-2)/2;
		float standdv = abs(mean)/3;

		gaussianpdf = new float[rangemax-2];
		gaussiancdf = new float[rangemax-2];
		gaussianlookup = new int[rangemax-2];

		for (int i=0; i<(rangemax-2); i++) {
			gaussianpdf[i]=exp(-(i-mean)*(i-mean)/(2*standdv*standdv))/sqrt(standdv * standdv * 2 * M_PI);

			if (i != 0) {
				gaussiancdf[i] = gaussiancdf[i-1]+gaussianpdf[i];
			}
			else {
				gaussiancdf[i] = gaussianpdf[i];
			}
		}

		for (int i=0; i<(rangemax-2); i++) {
			gaussiancdf[i] = graycdf[rangemax-3]*gaussiancdf[i]/gaussiancdf[rangemax-3];
		}

		for (int i=0; i<(rangemax-2); i++) {
			for (int j=0; j<(rangemax-2); j++) {
				if (graycdf[i] <= gaussiancdf[j]) {
					gaussianlookup[i] = j;

					break;
				}
			}
		}

		for (int i=0; i<(rangemax-2); i++) {
			grayhe[i] = floor(0.5+(((double)(rangemax-3)*graycdf[i])/graycdf[rangemax-3]));
		}

	}

	////// End of Histogram Equalization ///////

	if (emdata->is_complex()) {
		if (dsx != -1) {
			int l = y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int ll = x0;

				for (int i = xmin; i < xsize; i++) {
					if (l + ll > lmax || ll >= nx - 2) break;

					int k = 0;
					unsigned char p;

					if (ll >= nx / 2) {
						if (l >= (ny - inv_scale) * nx) k = 2 * (ll - nx / 2) + 2;
						else k = 2 * (ll - nx / 2) + l + 2 + nx;
					}
					else k = nx * ny - (l + 2 * ll) - 2;

					if (k >= mid) k -= mid; // These 2 lines handle the Fourier origin being in the corner, not the middle
					else k += mid;

					float t = image_data[k];
					int ph;

					// in color mode

					if (flags & 16 && asrgb>2) {
//						if (l >= (ny - inv_scale) * nx) ph = (int)(image_data[k+1]*768/(2.0*M_PI))+384;	 // complex phase as integer 0-767;
						if (ll >= nx / 2) ph = (int)(image_data[k+1]*768/(2.0*M_PI))+384;	 // complex phase as integer 0-767;
						else ph = (int)(-image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767;
					}

					if (t <= rm)  p = mingray;
					else if (t >= render_max) p = maxgray;
					else if (gamma != 1.0) {
						k=(int)(gs2 * (t-render_min)); // map float value to 0-4096 range
						p = gammamap[k]; // apply gamma using precomputed gamma map
//						p = (unsigned char) (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma);
//						p += mingray;
//						k = static_cast<int>( (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma) );
//						k += mingray;
					}
					else {
						p = (unsigned char) (gs * (t - render_min));
						p += mingray;
					}

					// color rendering

					if (flags & 16 && asrgb>2) {
						if (ph<256) {
							data[i * asrgb + j * bpl] = p*(255-ph)/256;
							data[i * asrgb + j * bpl+1] = p*ph/256;
							data[i * asrgb + j * bpl+2] = 0;
						}
						else if (ph<512) {
							data[i * asrgb + j * bpl+1] = p*(511-ph)/256;
							data[i * asrgb + j * bpl+2] = p*(ph-256)/256;
							data[i * asrgb + j * bpl] = 0;
						}
						else {
							data[i * asrgb + j * bpl+2] = p*(767-ph)/256;
							data[i * asrgb + j * bpl] = p*(ph-512)/256;
							data[i * asrgb + j * bpl+1] = 0;
						}
					}
					else data[i * asrgb + j * bpl] = p;
					if (hist) histd[p]++;
					ll += dsx;
				}

				l += dsy;
			}
		}
		else {
			remy = 10;
			int l = y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int br = l;
				remx = 10;
				int ll = x0;

				for (int i = xmin; i < xsize - 1; i++) {
					if (l + ll > lmax || ll >= nx - 2) {
						break;
					}

					int k = 0;
					unsigned char p;

					if (ll >= nx / 2) {
						if (l >= (ny * nx - nx)) k = 2 * (ll - nx / 2) + 2;
						else k = 2 * (ll - nx / 2) + l + 2 + nx;
					}
					else k = nx * ny - (l + 2 * ll) - 2;

					if (k >= mid) k -= mid; // These 2 lines handle the Fourier origin being in the corner, not the middle
					else k += mid;

					float t = image_data[k];
					// in color mode
					int ph;

					if (flags & 16 && asrgb>2) {
						if (l >= (ny * nx - nx)) ph = (int)(image_data[k+1]*768/(2.0*M_PI))+384;	 // complex phase as integer 0-767;
						else ph = (int)(-image_data[k+1]*768/(2.0*M_PI))+384;	// complex phase as integer 0-767;
					}

					if (t <= rm)
						p = mingray;
					else if (t >= render_max) {
						p = maxgray;
					}
					else if (gamma != 1.0) {
						k=(int)(gs2 * (t-render_min)); // map float value to 0-4096 range
						p = gammamap[k]; // apply gamma using precomputed gamma map

//						p = (unsigned char) (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma);
//						p += mingray;
//						k = static_cast<int>( (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma) );
//						k += mingray;
					}
					else {
						p = (unsigned char) (gs * (t - render_min));
						p += mingray;
					}

					if (flags & 16 && asrgb>2) {
						if (ph<256) {
							data[i * asrgb + j * bpl] = p*(255-ph)/256;
							data[i * asrgb + j * bpl+1] = p*ph/256;
							data[i * asrgb + j * bpl+2] = 0;
						}
						else if (ph<512) {
							data[i * asrgb + j * bpl+1] = p*(511-ph)/256;
							data[i * asrgb + j * bpl+2] = p*(ph-256)/256;
							data[i * asrgb + j * bpl] = 0;
						}
						else {
							data[i * asrgb + j * bpl+2] = p*(767-ph)/256;
							data[i * asrgb + j * bpl] = p*(ph-512)/256;
							data[i * asrgb + j * bpl+1] = 0;
						}
					}
					else data[i * asrgb + j * bpl] = p;

					if (hist) histd[p]++;

					ll += addi;
					remx += addr;

					if (remx > scale_n) {
						remx -= scale_n;
						ll++;
					}
				}

				l = br + addi * nx;
				remy += addr;

				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}
	else {
		if (dsx != -1) {
			int l = x0 + y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int br = l;

				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) {
						break;
					}

					int k = 0;
					unsigned char p;
					float t;

					if (dsx == 1) t=image_data[l];
					else { // This block does local pixel averaging for nicer reduced views
						t=0;

						if ((l+dsx+dsy) > lmax) {
							break;
						}

						for (int iii=0; iii<dsx; iii++) {
							for (int jjj=0; jjj<dsy; jjj+=nx) {
								t += image_data[l+iii+jjj];
							}
						}

						t /= dsx*(dsy/nx);
					}

					if (t <= rm) p = mingray;
					else if (t >= render_max) p = maxgray;
					else if (gamma != 1.0) {
						k=(int)(gs2 * (t-render_min)); // map float value to 0-4096 range
						p = gammamap[k];// apply gamma using precomputed gamma map
//						k = (int) (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma);
//						k += mingray;
//						k = static_cast<int>( (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma) );
//						k += mingray;
					}
					else {
						if (flags & 32) {
							//p = graylookup[(int)(t - render_min)];
							//graylookup[i]=(unsigned char)((maxgray-mingray-2)*grayhe[(int)(i*(rangemax-3)/(render_max-render_min))]/(rangemax-3)+1);
							if (flags & 64)
							{
								p=(unsigned char)(gaussianlookup[(int)(floor((rangemax-2)*(t - render_min)/(render_max-render_min)))]*(maxgray-mingray-2)/(rangemax-3)+1);
							}
							else {
								p=(unsigned char)(grayhe[(int)((t - render_min)*(rangemax-3)/(render_max-render_min))]*(maxgray-mingray-2)/(rangemax-3)+1);}
							//p=(unsigned char)gaussianlookup[(int)(ceil)((t - render_min)*(rangemax-3)/(render_max-render_min))]+1;
						}
						else
						{
							p=(unsigned char) (gs * (t - render_min));
						}

						p += mingray;
					}

					if (invert) {
						p = mingray + maxgray - p;
					}

					data[i * asrgb + j * bpl] = p;

					if (hist) histd[p]++;

					l += dsx;
				}

				l = br + dsy;
			}
		}
		else {
			remy = 10;
			int l = x0 + y0 * nx;

			for (int j = ymax; j >= ymin; j--) {
				int addj = addi;

				// There seems to be some overflow issue happening
				// where the statement if (l > lmax) break (below) doesn't work
				// because the loop that iterates jjj can inadvertantly go out of bounds

				if ((l + addi*nx) >= nxy) {
					addj = (nxy-l)/nx;

					if (addj <= 0) continue;
				}

				int br = l;
				remx = 10;

				for (int i = xmin; i < xsize; i++) {
					if (l > lmax) break;
					int k = 0;
					unsigned char p;
					float t;

					if (addi <= 1) t = image_data[l];
					else { // This block does local pixel averaging for nicer reduced views
						t=0;
						for (int jjj=0; jjj<addj; jjj++) {
							for (int iii=0; iii<addi; iii++) {
								t += image_data[l+iii+jjj*nx];
							}
						}

						t /= addi*addi;
					}

					if (t <= rm) p = mingray;
					else if (t >= render_max) p = maxgray;
					else if (gamma != 1.0) {
						k=(int)(gs2 * (t-render_min)); // map float value to 0-4096 range
						p = gammamap[k]; // apply gamma using precomputed gamma map
//						k = (int) (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma);
//						k += mingray;
//						k = static_cast<int>( (maxgray-mingray)*pow((gs2 * (t - render_min)),gamma) );
//						k += mingray;
					}
					else {
						if (flags & 32) {
							// p = graylookup[(int)(t - render_min)];

							if (flags & 64) {
								p = (unsigned char)(gaussianlookup[(int)(floor((rangemax-2)*(t - render_min)/(render_max-render_min)))]*(maxgray-mingray-2)/(rangemax-3)+1);
							}
							else {
								p = (unsigned char)(grayhe[(int)((t - render_min)*(rangemax-3)/(render_max-render_min))]*(maxgray-mingray-2)/(rangemax-3)+1);
							}

							// p = (unsigned char)gaussianlookup[(int)(ceil)((t - render_min)*(rangemax-3)/(render_max-render_min))]+1;
						}
						else {
							p = (unsigned char) (gs * (t - render_min));
						}

						p += mingray;
					}

					if (invert) {
						p = mingray + maxgray - p;
					}

					data[i * asrgb + j * bpl] = p;

					if (hist) histd[p]++;

					l += addi;
					remx += addr;

					if (remx > scale_n) {
						remx -= scale_n;
						l++;
					}
				}

				l = br + addi * nx;
				remy += addr;

				if (remy > scale_n) {
					remy -= scale_n;
					l += nx;
				}
			}
		}
	}

	// this replicates r -> g,b

	if (asrgb == 3 && !(flags & 16)) {
		for (int j=ymin*bpl; j <= ymax*bpl; j+=bpl) {
			for (int i=xmin; i<xsize*3; i+=3) {
				data[i+j+1] = data[i+j+2] = data[i+j];
			}
		}
	}

	if (asrgb == 4 && !(flags & 16)) {
		for (int j=ymin*bpl; j <= ymax*bpl; j+=bpl) {
			for (int i=xmin; i<xsize*4; i+=4) {
				data[i+j+1] = data[i+j+2] = data[i+j+3] = data[i+j];
				data[i+j+3] = 255;
			}
		}
	}

	EXITFUNC;

	// ok, ok, not the most efficient place to do this, but it works

	if (invy) {
		int x,y;
		char swp;

		for (y=0; y<iysize/2; y++) {
			for (x=0; x<ixsize; x++) {
				swp = ret[y*bpl+x];
				ret[y*bpl+x] = ret[(iysize-y-1)*bpl+x];
				ret[(iysize-y-1)*bpl+x] = swp;
			}
		}
	}

	// return PyString_FromStringAndSize((const char*) data,iysize*bpl);

	if (flags & 16) {
		glDrawPixels(ixsize, iysize, GL_LUMINANCE, GL_UNSIGNED_BYTE,
			 (const GLvoid *)ret.data());
	}

	delete [] grayhe;
	delete [] gaussianpdf;
	delete [] gaussiancdf;
	delete [] gaussianlookup;
	delete [] graypdftwo;

	return ret;
}

// DEPRECATED

unsigned long GLUtil::get_isosurface_dl(MarchingCubes* mc,
		 unsigned int tex_id, bool surface_face_z, bool recontour)
{
	if (mc->_isodl != 0) glDeleteLists(mc->_isodl,1);
	if (recontour) mc->calculate_surface();
	if (surface_face_z) mc->surface_face_z();
	if (mc->getRGBmode()) mc->color_vertices();

#if MARCHING_CUBES_DEBUG
	cout << "There are " << ff.elem()/3 << " faces and " << pp.elem() <<
		 " points and " << nn.elem() << " normals to render in generate dl" << endl;
#endif

	int maxf;
//#ifdef	_WIN32
//	glGetIntegerv(GL_MAX_ELEMENTS_INDICES_WIN, & maxf);
//#else
	glGetIntegerv(GL_MAX_ELEMENTS_INDICES, & maxf);
//#endif	//_WIN32

#if MARCHING_CUBES_DEBUG
	int maxv;
	glGetIntegerv(GL_MAX_ELEMENTS_VERTICES, & maxv);

	cout << "Max vertices is " << maxv << " max indices is " << maxf << endl;
	cout << "Using OpenGL " << glGetString(GL_VERSION) << endl;
#endif

	if (recontour) {
		for (unsigned int i = 0; i < mc->ff.elem(); ++i) mc->ff[i] /= 3;
	}

	if (maxf % 3 != 0) {
		maxf = maxf - (maxf%3);
	}

	if (tex_id != 0) {
		// Normalize the coordinates to be on the interval 0,1

		mc->pp.mult3(1.0f/(float) mc->_emdata->get_xsize(),
			 1.0f/(float)mc->_emdata->get_ysize(),
			 1.0f/(float)mc->_emdata->get_zsize());
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glTexCoordPointer(3, GL_FLOAT, 0, mc->pp.get_data());
	}
	else {
		glEnableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glNormalPointer(GL_FLOAT, 0,  mc->nn.get_data());
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0,  mc->pp.get_data());

	if (mc->getRGBmode()) {
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT, 0,  mc->cc.get_data());
	}

	mc->_isodl = glGenLists(1);

#if MARCHING_CUBES_DEBUG
	int time0 = clock();
#endif

	glNewList(mc->_isodl,GL_COMPILE);

	if (tex_id != 0) {
		glEnable(GL_TEXTURE_3D);
		glBindTexture(GL_TEXTURE_3D, tex_id);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	}

	// Drawing range elements based on the output of
	// glGetIntegerv(GL_MAX_ELEMENTS_INDICES, & maxf);
	// Saved about 60% of the time... drawRange should probably always be true

	bool drawRange = true;

	if (drawRange == false) {
		glDrawElements(GL_TRIANGLES,mc->ff.elem(),GL_UNSIGNED_INT, mc->ff.get_data());
	} else {
		for(unsigned int i = 0; i < mc->ff.elem(); i+=maxf)
		{
			if ((i+maxf) > mc->ff.elem())
				glDrawElements(GL_TRIANGLES, mc->ff.elem()-i, GL_UNSIGNED_INT, &(mc->ff[i]));
			else
				glDrawElements(GL_TRIANGLES, maxf, GL_UNSIGNED_INT, &(mc->ff[i]));

				// glDrawRangeElements is part of the extensions,
				// we might want to experiment with its performance at some stage,
				// so please leave this code here, commented out.
				// This is an either or situation, so if glDrawRangeElements is used,
				// glDrawElements above would have to be commented out.
				// glDrawRangeElements(GL_TRIANGLES, 0, 0, maxf, GL_UNSIGNED_INT, &ff[i]);
		}
	}

	if (tex_id != 0) glDisable(GL_TEXTURE_3D);

	glEndList();

#if MARCHING_CUBES_DEBUG
	int time1 = clock();

	cout << "It took " << (time1-time0) << " " <<
		 (float)(time1-time0)/CLOCKS_PER_SEC << " to draw elements" << endl;
#endif

	return mc->_isodl;
}

void GLUtil::contour_isosurface(MarchingCubes* mc)
{
	mc->calculate_surface();

	// What does this do???

	for (unsigned int i = 0; i < mc->ff.elem(); ++i ) mc->ff[i] /= 3;

	// Need to rebind data (to the GPU)

	mc->needtobind = true;
}

void GLUtil::render_using_VBOs(MarchingCubes* mc, unsigned int tex_id,
		 bool surface_face_z)
{
	// In current version Texture is not supported b/c it is not used... EVER

#ifdef _WIN32
	typedef void (APIENTRYP PFNGLGENBUFFERSPROC) (GLsizei n, GLuint *buffers);
	PFNGLGENBUFFERSPROC glGenBuffers;
	glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");

	typedef GLboolean (APIENTRYP PFNGLISBUFFERPROC) (GLuint buffer);
	PFNGLISBUFFERPROC glIsBuffer;
	glIsBuffer = (PFNGLISBUFFERPROC) wglGetProcAddress("glIsBuffer");
#endif	//_WIN32

	if (surface_face_z) mc->surface_face_z();

	// Bug in OpenGL, sometimes glGenBuffers doesn't work.
	// Try again, if we still fail then return...

//	printf("%d %d %d %d     ",mc->buffer[0],mc->buffer[1],mc->buffer[2],mc->buffer[3]);
	if (!glIsBuffer(mc->buffer[0])) glGenBuffers(4, mc->buffer);
//	printf("%d %d %d %d\n",mc->buffer[0],mc->buffer[1],mc->buffer[2],mc->buffer[3]);

	// whenever something changes, like color mode or color scale (or threshold),
	// we need to recolor

	if (mc->getRGBmode() && (mc->rgbgenerator.getNeedToRecolor() ||
		 mc->needtobind)) {
		mc->color_vertices();
		mc->needtobind = true;
	}

	int maxf;
	glGetIntegerv(GL_MAX_ELEMENTS_INDICES, & maxf);

	if (maxf % 3 != 0) {
		maxf = maxf - (maxf%3);
	}

#ifdef _WIN32
	typedef void (APIENTRYP PFNGLBINDBUFFERPROC) (GLenum target, GLuint buffer);
	PFNGLBINDBUFFERPROC glBindBuffer;
	glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");

	typedef void (APIENTRYP PFNGLBUFFERDATAPROC) (GLenum target, GLsizeiptr size,
		 const GLvoid *data, GLenum usage);
	PFNGLBUFFERDATAPROC glBufferData;
	glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
#endif	//_WIN32

	// Normal vectors

	glEnableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, mc->buffer[2]);

	if (mc->needtobind) {
		glBufferData(GL_ARRAY_BUFFER, 4*mc->nn.elem(), mc->nn.get_data(),
				 GL_STATIC_DRAW);
	}

	glNormalPointer(GL_FLOAT,0,0);

	// Vertex vectors

	glBindBuffer(GL_ARRAY_BUFFER, mc->buffer[0]);

	if (mc->needtobind) {
		glBufferData(GL_ARRAY_BUFFER, 4*mc->pp.elem(), mc->pp.get_data(),
				 GL_STATIC_DRAW);
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3,GL_FLOAT,0,0);

	// Color vectors

	if (mc->getRGBmode()) {
		glEnableClientState(GL_COLOR_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, mc->buffer[3]);

		if (mc->needtobind) {
			glBufferData(GL_ARRAY_BUFFER, 4*mc->cc.elem(), mc->cc.get_data(),
				 GL_STATIC_DRAW);
		}

		glColorPointer(3,GL_FLOAT,0, 0);
	}

	// Indices

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mc->buffer[1]);

	if (mc->needtobind) {
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4*mc->ff.elem(), mc->ff.get_data(),
				 GL_STATIC_DRAW);
	}

	// This lets us know if buffers an not implmemted

	if (!glIsBuffer(mc->buffer[0])) {
		cout << "Can't Generate GL Vertex Buffers. glGenBuffer error" << endl;

		return;
	}

	// finally draw the elemenets

	glDrawElements(GL_TRIANGLES, mc->ff.elem(), GL_UNSIGNED_INT, 0);

	// No longer need to bind data to the GPU

	mc->needtobind = false;
}

static void transpose_4_by_4_matrix (vector <float> & v)
{
	std::swap (v [ 1], v [ 4]);
	std::swap (v [ 2], v [ 8]);
	std::swap (v [ 3], v [12]);
	std::swap (v [ 6], v [ 9]);
	std::swap (v [ 7], v [13]);
	std::swap (v [11], v [14]);
}

void GLUtil::glLoadMatrix (const Transform & xform)
{
	vector <float> xformlist = xform.get_matrix_4x4 ();

	transpose_4_by_4_matrix (xformlist);

	glLoadMatrixf (reinterpret_cast <GLfloat *> (& xformlist [0]));
}

void GLUtil::glMultMatrix (const Transform & xform)
{
	vector <float> xformlist = xform.get_matrix_4x4 ();

	transpose_4_by_4_matrix (xformlist);

	glMultMatrixf (reinterpret_cast <GLfloat *> (& xformlist [0]));
}

// This draws a bounding box, or any box

void GLUtil::glDrawBoundingBox(float width, float height, float depth)
{
	float x = width/2.0f;
	float y = height/2.0f;
	float z = depth/2.0f;

	float vertices[72] = {-x,-y,-z, x,-y,-z, x,-y,-z, x,y,-z, x,y,-z, -x,y,-z, -x,y,-z, -x,-y,-z,
		-x,-y,-z, -x,-y,z, x,-y,-z, x,-y,z, x,y,-z, x,y,z, -x,y,-z, -x,y,z,
		-x,-y,z, x,-y,z, x,-y,z, x,y,z, x,y,z, -x,y,z, -x,y,z, -x,-y,z };

	glBegin(GL_LINES);
	for (int i=0; i<72; i+=3) glVertex3f(vertices[i],vertices[i+1],vertices[i+2]);
	glEnd();

// 	float vertices[24] = {-w2,  h2,  d2, w2,  h2,  d2, w2, -h2,  d2, -w2,
// 		 -h2,  d2, -w2,  h2, -d2, w2,  h2, -d2, w2, -h2, -d2, -w2, -h2, -d2};
// 	int indices[24] = {0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
// 		 6, 7, 7, 4, 0, 4, 3, 7, 1, 5, 2, 6};
//
// #ifdef _WIN32
// 	typedef void (APIENTRYP PFNGLGENBUFFERSPROC) (GLsizei n, GLuint *buffers);
// 	PFNGLGENBUFFERSPROC glGenBuffers;
// 	glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");
//
// 	typedef GLboolean (APIENTRYP PFNGLISBUFFERPROC) (GLuint buffer);
// 	PFNGLISBUFFERPROC glIsBuffer;
// 	glIsBuffer = (PFNGLISBUFFERPROC) wglGetProcAddress("glIsBuffer");
//
// 	typedef void (APIENTRYP PFNGLBINDBUFFERPROC) (GLenum target, GLuint buffer);
// 	PFNGLBINDBUFFERPROC glBindBuffer;
// 	glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");
//
// 	typedef void (APIENTRYP PFNGLBUFFERDATAPROC) (GLenum target, GLsizeiptr size,
// 		 const GLvoid *data, GLenum usage);
// 	PFNGLBUFFERDATAPROC glBufferData;
// 	glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
// #endif	//_WIN32
//
// 	if (glIsBuffer(GLUtil::buffer[0]) == GL_FALSE) {
// 		glGenBuffers(2, GLUtil::buffer);
// 	}
//
// 	// Could use dirty bit here but not worth my time to implment
//
// 	glBindBuffer(GL_ARRAY_BUFFER, GLUtil::buffer[0]);
// 	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
// 	glEnableClientState(GL_VERTEX_ARRAY);
// 	glVertexPointer(3,GL_FLOAT,0,0);
//
// 	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GLUtil::buffer[1]);
// 	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
// 		 GL_STATIC_DRAW);
// 	glDrawElements(GL_LINES,24,GL_UNSIGNED_INT,0);
}

void GLUtil::glDrawDisk(float radius, int spokes)
{
#ifdef _WIN32
	typedef void (APIENTRYP PFNGLGENBUFFERSPROC) (GLsizei n, GLuint *buffers);
	PFNGLGENBUFFERSPROC glGenBuffers;
	glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");

	typedef GLboolean (APIENTRYP PFNGLISBUFFERPROC) (GLuint buffer);
	PFNGLISBUFFERPROC glIsBuffer;
	glIsBuffer = (PFNGLISBUFFERPROC) wglGetProcAddress("glIsBuffer");

	typedef void (APIENTRYP PFNGLBINDBUFFERPROC) (GLenum target, GLuint buffer);
	PFNGLBINDBUFFERPROC glBindBuffer;
	glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");

	typedef void (APIENTRYP PFNGLBUFFERDATAPROC) (GLenum target, GLsizeiptr size,
		 const GLvoid *data, GLenum usage);
	PFNGLBUFFERDATAPROC glBufferData;
	glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
#endif	//_WIN32

	//This code is experimental
	int arraysize = 4*pow((double)2, (double)spokes);
	int sideofarray = pow((double)2, (double)spokes);

//	float vertices[3*arraysize + 3];
	vector<float> vertices(3*arraysize + 3);

	// last vertex is center

	vertices[3*arraysize] = 0.0;
	vertices[3*arraysize+1] = 0.0;
	vertices[3*arraysize+2] = 0.0;

	vertices[0] = 0.0;
	vertices[1] = radius;
	vertices[2] = 0.0;

	vertices[sideofarray*3] = radius;
	vertices[sideofarray*3 + 1] = 0.0;
	vertices[sideofarray*3 + 2] = 0.0;

	vertices[sideofarray*6] = 0.0;
	vertices[sideofarray*6 + 1] = -radius;
	vertices[sideofarray*6 + 2] = 0.0;

	vertices[sideofarray*9] = -radius;
	vertices[sideofarray*9 + 1] = 0.0;
	vertices[sideofarray*9 + 2] = 0.0;

	// This could aslo be implemented recursively

	for (int step = 0; step < spokes; step++) {
		// starting location
		int x = sideofarray/pow(2.0,(double)(step+1));

		for (int i = 1; i <= 4*pow(2.0,(double)step); i++) {
			// take the necessary steps

			int index =  x + 2*x*(i-1);
			int index_f = (index + x) % arraysize;
			int index_i = index - x;

			cout << index << " " << index_f << " " << index_i << endl;

			// need to resclae length to that of radius

			vertices[index_f*3] = (vertices[index_f*3] - vertices[index_i*3])/2.0f;
			vertices[index_f*3 + 1] = (vertices[index_f*3 + 1] - vertices[index_i*3 + 1])/2.0f;
			vertices[index_f*3 + 2] = (vertices[index_f*3 + 2] - vertices[index_i*3 + 2])/2.0f;
		}
	}

	// GL stuff

	if (glIsBuffer(GLUtil::buffer[0]) == GL_FALSE) {
		glGenBuffers(2, GLUtil::buffer);
	}

	// Could use dirty bit here but not worth my time to implement

	glBindBuffer(GL_ARRAY_BUFFER, GLUtil::buffer[0]);

	// glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBufferData(GL_ARRAY_BUFFER, vertices.size(), &(vertices[0]),
		 GL_STATIC_DRAW);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	// Code to select indices
}
#endif // USE_OPENGL
