/*
 * Author: Chao Yang
 * Copyright (c) 2000-2006
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
 */

//#ifndef sparx__pca_h__
//#define sparx__pca_h__

#include "emdata.h"

namespace EMAN {
   class PCA {
      public: 
         vector<float>   singular_vals;        // singular values
         vector<EMData*> eigenimages; // eigenimages

         // pca in core
         int dopca(vector <EMData*> imgstack, EMData *mask);
         int dopca(vector <EMData*> imgstack, EMData *mask, int nvec);
         int dopca_lan(vector <EMData*> imgstack, EMData *mask, int nvec);

         // pca out of core
         int dopca_ooc(const string &filename_in, const string &filename_out, 
                       const string &lanscratch,  EMData *mask, int nvec);

         // Lanczos factorization (used by dopca_lan)
         int Lanczos(vector <EMData*> imgstack, int *maxiter, 
                     float  *diag, float *subdiag, float *V, float *beta);

         // Lanczos factorization out-of-core (used by dopca_ooc)
         int Lanczos_ooc(string const& filename_in, int *kstep, 
                         float  *diag, float *subdiag, 
                         string const& lanscratch, float  *beta);

         // retrieve singular vectors and values
         vector<float>   get_vals();
         vector<EMData*> get_vecs();

         // flush singular values and vectors
         void clear();
   };
}

//#endif //sparx__pca_h__
