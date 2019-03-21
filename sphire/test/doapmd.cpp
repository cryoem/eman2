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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "log.h"
#include "emdata.h"
#include "xydata.h"
#include "assert.h"
#include "projector.h"

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double mytimer()
{
    struct tms use;
    double tmp;
    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / CLK_TCK;
}

using namespace EMAN;

using std::cout;
using std::cin;
using std::endl;
using std::string;

#include "spidutil.h"

int main()
{
    EMData *volume = new EMData(); // initial volume
    EMData *expimg = new EMData(); // experimental image

    Dict    refparams;
    APMDopt options;

    vector <float> angles;
    float   delta, tlb, tub, plb, pub;
    int     ndim, nx, ny, nz, nang, status;
    int     nref, nimg, i;
    double  t0, t1;
    float   *newangles;

    printf("reading a 3-D volume...\n");
    volume->read_image("vol001.tfc");

    ndim = volume->get_ndim();
    nx   = volume->get_xsize();
    ny   = volume->get_ysize();
    nz   = volume->get_zsize();

    // print some stats
    printf("ndim = %d, nx = %d, ny = %d, nz = %d\n", ndim, nx, ny, nz);

    // generate a set of reference projections from the volume

    delta = 15.0;
    tlb   = 0.0;
    tub   = 90.0;
    plb   = 0.0;
    pub   = 359.9;

    // angles used to generate reference projections
    angles = Util::voea(delta, tlb, tub, plb, pub); 
    nang   = angles.size()/3;
    printf("size(angles) = %d\n", nang);

    refparams["anglelist"] = angles;
    refparams["angletype"] = "SPIDER";
    refparams["radius"]    = 30.0;
    t0 = mytimer();
    printf("generating reference projections...\n");
    EMData* refprj = volume->project("pawel", refparams);
    t1 = mytimer() - t0;
    nref  = refprj->get_zsize();
    printf("nref = %d\n", nref);
  
    printf("time used in ref projection calculation = %11.3e\n", t1);

    // read experimental images
    expimg->read_image("tf2d0001.tfc"); 
    nimg = expimg->get_zsize();

    options.mr    = 5;
    options.nr    = 30;
    options.iskip = 1;
    options.mode  = 'F'; 

    // run APMD on expimg using refprj as the reference, 

    printf("running APMD...\n");

    newangles = new float[nimg*3];
    status = apmd(refprj, refparams, expimg, options, newangles);
    if (status != 0) {
       fprintf(stderr, "APMD failed!\n");
       goto EXIT;
    }

    for (i = 1; i<=nimg; i++) 
       printf("psi = %g, theta = %g, phi = %g\n", 
              newangles(1,i), newangles(2,i), newangles(3,i));

EXIT:
    if (volume)    delete volume;
    if (expimg)    delete expimg;
    if (refprj)    delete refprj;
    if (newangles) delete newangles;
    return status;
}

