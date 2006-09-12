/**
 * $Id$
 */

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

#include "emdata.h"
#include "analyzer.h"
#include "sparx/analyzer_sparx.h"
#include "util.h"

using namespace EMAN;

namespace EMAN {

	template <> Factory < Analyzer >::Factory()
	{
		force_add(&PcaAnalyzer::NEW);
	}

}

#define covmat(i,j) covmat[ ((j)-1)*nx + (i)-1 ]
#define imgdata(i)  imgdata[ (i)-1 ]
int PcaAnalyzer::insert_image(EMData * image)
{
   EMData *maskedimage = Util::compress_image_mask(image,mask);

   int nx = maskedimage->get_xsize();
   float *imgdata = maskedimage->get_data();
   if (nx != ncov) {
      fprintf(stderr,"insert_image: something is wrong...\n");
      exit(1);
   }

   // there is a faster version of the following rank-1 update 
   for (int j = 1; j <= nx; j++)
       for (int i = 1; i<=nx; i++) {
           covmat(i,j) += imgdata(i)*imgdata(j);
   }   

   return 0;
}
#undef covmat

#define eigvec(i,j) eigvec[(j)*ncov + (i)]
vector<EMData*> PcaAnalyzer::analyze()
{
	int status = 0;
	printf("start analyzing..., ncov = %d\n", ncov);
        eigval = (float*)calloc(ncov,sizeof(float));
        eigvec = (float*)calloc(ncov*ncov,sizeof(float));
        status = Util::coveig(ncov, covmat, eigval, eigvec);

        for (int i=1; i<=5; i++) printf("eigval = %11.4e\n", 
            eigval[ncov-i]);

        // pack eigenvectors into the return imagelist
        EMData *eigenimage = new EMData();
        eigenimage->set_size(ncov,1,1);
        float *rdata = eigenimage->get_data();
        for (int j = 1; j<= nvec; j++) {
	    for (int i = 0; i < ncov; i++)
		rdata[i] = eigvec(i,ncov-j);
	    images.push_back(Util::reconstitute_image_mask(eigenimage,mask));
        }

        free(eigvec);
        EMDeletePtr(eigenimage); 

	return images;
}
#undef eigvec

void PcaAnalyzer::set_params(const Dict & new_params)
{
	params = new_params;
	mask = params["mask"];
	nvec = params["nvec"];

        // count the number of pixels under the mask
        // (this is really ugly!!!)
        EMData *dummy = new EMData();

        int nx = mask->get_xsize();
        dummy->set_size(nx,nx);
        EMData *dummy1d = Util::compress_image_mask(dummy,mask);
        ncov = dummy1d->get_xsize();
        EMDeletePtr(dummy);
        EMDeletePtr(dummy1d);

	// allocate and set up the covriance matrix
	covmat = (float*)calloc(ncov*ncov,sizeof(float));
}

