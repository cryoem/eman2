void Utilit1(double *, double *, int );
  
void Utilit2(double *X, double *X1, double *Y, double *D, double *dd, double xk, int l, float (*my_func)(EMData* , EMData* , EMData* , float , float , float), EMData *image, EMData *refim, EMData *mask);
    
void Derivatives(double *X, double *D, double *Y, double *dd, double xk, int l, float (*my_func)(EMData* , EMData* , EMData* , float , float , float), EMData *, EMData *,
EMData *);

void Steepda(double *X, double xk, double e, int l, int m, int *n, float (*my_func)(EMData* , EMData* , EMData* , float , float , float), EMData *, EMData *, EMData *);

void Utilit2_G(double *X, double *X1, double *Y, double *D, double *dd, double xk, int l, float (*my_func)(EMData* , EMData* , EMData* , Util::KaiserBessel& , float , float , float), EMData *image,
EMData *refim, EMData *mask, Util::KaiserBessel& kb);

void Derivatives_G(double *X, double *D, double *Y, double *dd, double xk, int l, float (*my_func)(EMData* , EMData* , EMData* , Util::KaiserBessel& , float , float , float), EMData *image, EMData
*refim, EMData *mask, Util::KaiserBessel& kb);

void Steepda_G(double *X, double xk, double e, int l, int m, int *n, float (*my_func)(EMData* , EMData* , EMData* , Util::KaiserBessel& , float , float , float), EMData *image, EMData *refim, EMData
*mask, Util::KaiserBessel& kb);
