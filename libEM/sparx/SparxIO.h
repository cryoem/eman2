#ifndef __SparxIO__
#define  __SparxIO__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLEN 200
void* SparxMalloc(size_t nbytes);

class TFList 
{
  public:
     float *data;
     int   nrows;
     int   ncols;
     int   ndigit; // if ndigit  = 0, use fixed point format
                   // if ndigit  > 0, use floating point format with ndigit
                   //                 after the decimal point.

     int   read(char *filename);
     int   write(char *filename);
     void  Copy(float *rdata);
     void  CopyCol(int jcol, float *rcol);
     void  CopyRow(int jcol, float *rcol);
     void  SetVal(int irow, int jcol, float val);
     float GetVal(int irow, int jcol);
     TFList();
     TFList(int nr, int nc);
     TFList(int nr, int nc, int ndig);
     ~TFList(); 
};

#endif  __SparxIO__
