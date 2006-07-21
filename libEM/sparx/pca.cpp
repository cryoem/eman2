#include "emdata.h"
#include "util.h"

using namespace EMAN;

vector <EMData*> dopca(vector <EMData*> imgstack, EMData *mask)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int i;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;
   vector<EMData*> eigimages;

   int nimg = imgstack.size();
   // printf("number of images in the stack = %d\n", nimg);

   for (i=0; i<nimg; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   // for right now, compute a full SVD
   eigvecs = Util::svdcmp(img1dlst, 0);

   img1dlst.clear();

   for (i=0; i<nimg; i++) {
      eigimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return eigimages;
}

