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
