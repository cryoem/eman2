//#ifndef sparx__pca_h__
//#define sparx__pca_h__

#include "emdata.h"

namespace EMAN {
   class PCA {
      public: 
         float *svalue; // singular values

         // pca in core, returns all singular vectors
         vector <EMData*> dopca(vector <EMData*> imgstack, EMData *mask);

         // pca in core, returns a subset of right singular vectors
         vector <EMData*> dopca(vector <EMData*> imgstack, EMData *mask, int nvec);
         // pca out of core
         char *dopca_ooc(const string &filename, EMData *mask, int nvec);
   };
}

//#endif //sparx__pca_h__
