#include "emdata.h"
#include "util.h"
#include "emutil.h"

#include "pca.h"

using namespace EMAN;

// return all right singular vectors
vector <EMData*> PCA::dopca(vector <EMData*> imgstack, EMData *mask)
{
   // performs PCA on a list of images (each under a mask)
   // returns a list of eigenimages

   int i;

   vector<EMData*> img1dlst;
   vector<EMData*> eigvecs;
   vector<EMData*> eigimages;

   int nimgs = imgstack.size();
   // printf("number of images in the stack = %d\n", nimgs);

   for (i=0; i<nimgs; i++) {
      img1dlst.push_back(Util::compress_image_mask(imgstack[i],mask));
   }

   // for right now, compute a full SVD
   eigvecs = Util::svdcmp(img1dlst, 0);

   img1dlst.clear();

   for (i=0; i<nimgs; i++) {
      eigimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return eigimages;
}

// return a subset of right singular vectors
vector <EMData*> PCA::dopca(vector <EMData*> imgstack, EMData *mask, int nvec)
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

   if (nimg < nvec) nvec = nimg;

   for (i=0; i<nvec; i++) {
      eigimages.push_back(Util::reconstitute_image_mask(eigvecs[i],mask));
   }

   eigvecs.clear();

   return eigimages;
}

char* PCA::dopca_ooc(const string &filename_in, EMData *mask, int nvec)
{
   char *filename_out = NULL;
   EMData *image_raw = new EMData();
   EMData *image_masked;

   int nimgs = EMUtil::get_image_count(filename_in);
   for (int i=0; i<nimgs; i++) {
       image_raw->read_image(filename_in, i);      
       image_masked=Util::compress_image_mask(image_raw,mask);
       image_masked->write_image("temp_masked_imaged.img",i); 
   }

   return filename_out;
}

