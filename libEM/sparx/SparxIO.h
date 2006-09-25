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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

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

#endif // __SparxIO__
