int setulb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, double *g, 
			 double *factr, double *pgtol, double *wa, long int *iwa, char *task, long int *iprint, char *csave, long int *lsave,
			 long int *isave, double *dsave, long int task_len, long int csave_len);

int mainlb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, 
			 double *g, double *factr, double *pgtol, double *ws, double *wy,
			 double *sy, double *ss, double *yy, double *wt, double *wn, double *snd, 
			 double *z__, double *r__, double *d__, double *t, double *wa, double *sg, 
			 double *sgo, double *yg, double *ygo, long int *index, long int *iwhere, long int *indx2, 
			 char *task, long int *iprint, char *csave, long int *lsave, long int *isave, double *dsave, 
			 long int task_len, long int csave_len);
			 
int active_(long int *n, double *l, double *u, long int *nbd, double *x, long int *iwhere, 
			 long int *iprint, long int *prjctd, long int *cnstnd, long int *boxed);			
			  
int bmv_(long int *m, double *sy, double *wt, long int *col, double *v, double *p, long int *info);			 

int cauchy_(long int *n, double *x, double *l, double *u, long int *nbd, double *g, long int *iorder,
			 long int *iwhere, double *t, double *d__, double *xcp, long int *m, double *wy, double *ws,
			 double *sy, double *wt, double *theta, long int *col, long int *head, double *p, double *c__,
			 double *wbp, double *v, long int *nint, double *sg, double *yg, long int *iprint, double *sbgnrm, long int *info);
			 
int cmprlb_(long int *n, long int *m, double *x, double *g, double *ws, double *wy, double *sy,
			 double *wt, double *z__, double *r__, double *wa, long int *index, double *theta,
			 long int *col, long int *head, long int *nfree, long int *cnstnd, long int *info);
		
int errclb_(long int *n, long int *m, double *factr, double *l, double *u, long int *nbd, char *task,
			 long int *info, long int *k, long int task_len);
			 
int formk_(long int *n, long int *nsub, long int *ind, long int *nenter, long int *ileave, long int *indx2, long int *iupdat, 
			long int *updatd, double *wn, double *wn1, long int *m, double *ws, double *wy, double *sy, 
			double *theta, long int *col, long int *head, long int *info);
						 
int formt_(long int *m, double *wt, double *sy, double *ss, long int *col, double *theta, long int *info);

int freev_(long int *n, long int *nfree, long int *index, long int *nenter, long int *ileave, long int *indx2, long int *iwhere, 
			long int *wrk, long int *updatd, long int *cnstnd, long int *iprint, long int *iter);
			
int hpsolb_(long int *n, double *t, long int *iorder, long int *iheap);
			
int lnsrlb_(long int *n, double *l, double *u, long int *nbd, double *x, double *f, double *fold,
			 double *gd, double *gdold, double *g, double *d__, double *r__, double *t, 
			 double *z__, double *stp, double *dnorm, double *dtd, double *xstep, double *stpmx,
			 long int *iter, long int *ifun, long int *iback, long int *nfgv, long int *info, char *task, long int *boxed, 
			 long int *cnstnd, char *csave, long int *isave, double *dsave, long int task_len, long int csave_len);

int matupd_(long int *n, long int *m, double *ws, double *wy, double *sy, double *ss,
			 double *d__, double *r__, long int *itail, long int *iupdat, long int *col, long int *head,
			 double *theta, double *rr, double *dr, double *stp, double *dtd);

int prn1lb_(long int*, long int*, double*, double*, double*, long int*, long int*, double*);

int prn2lb_(long int *n, double *x, double *f, double *g, long int *iprint, long int *itfile,
			 long int *iter, long int *nfgv, long int *nact, double *sbgnrm, long int *nint, char *word,
			 long int *iword, long int *iback, double *stp, double *xstep, long int word_len);

int prn3lb_(long int *n, double *x, double *f, char *task, long int *iprint, long int *info, 
			 long int *itfile, long int *iter, long int *nfgv, long int *nintol, long int *nskip, long int *nact,
			 double *sbgnrm, double *time, long int *nint, char *word, long int *iback, double *stp,
			 double *xstep, long int *k, double *cachyt, double *sbtime, double *lnscht,
			 long int task_len, long int word_len);
			 
int projgr_(long int *n, double *l, double *u, long int *nbd, double *x, double *g, double *sbgnrm);

int subsm_(long int *n, long int *m, long int *nsub, long int *ind, double *l, double *u, long int *nbd,
			double *x, double *d__, double *ws, double *wy, double *theta, long int *col,
			long int *head, long int *iword, double *wv, double *wn, long int *iprint, long int *info);

int dcsrch_(double *f, double *g, double *stp, double *ftol, double *gtol, double *xtol,
			 double *stpmin, double *stpmax, char *task, long int *isave, double *dsave, long int task_len);

int dcstep_(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy,
			 double *stp, double *fp, double *dp, long int *brackt, double *stpmin, double *stpmax);
			 
int timer_(double *ttime);
 
double dnrm2_(long int *n, double *x, long int *incx);

double dpmeps_();

int daxpy_(long int *n, double *da, double *dx, long int *incx, double *dy, long int *incy);

int dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy); 

double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);

int dpofa_(double *a, long int *lda, long int *n, long int *info);

int dscal_(long int *n, double *da, double *dx, long int *incx);

int dtrsl_(double *t, long int *ldt, long int *n, double *b, long int *job, long int *info);

long int s_cmp(char *str1, const char *const str2, long int l1, long int l2);

int s_copy(char *str1, const char *const str2, long int l1, long int l2);



