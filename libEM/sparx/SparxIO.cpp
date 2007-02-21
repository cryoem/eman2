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


#include "SparxIO.h"
#include "Tokenizer.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// this actually doesn't belong here
void* SparxMalloc(size_t nbytes)
{
   return malloc(nbytes);
   // do some memory usage count later
}

/*
Purpose:
--------

The TFList class is created to facilitate the SPARX
Text File input and output.

The SPARX Text File is define as an ASCII file that
contain nrows x ncols real type (4 byte) data. 

The files will be used in SPARX to store intermediate results 
such as angles, shifts, CTF parameter, cross-correlation coefficients etc.

The Text File should be easy to browse.  We should be able to 
edit the file (e.g., insert a line, delete a column) by any text 
editor.

Usage:
------

1. An TFList object:

In order to read a SPARX Text File, we must first declare and create
a TFList object, which can be done by, for example, the following line 
of code:

TFList anglist;

2. Reading a SPARX Text File

Once a TFList object is defined, we can use it to read a Text File. For
example, the following line of code reads the SPARX Text File angles.txt.

anglist.read("angles.txt");

Let's assume anglex.txt contains the following entries

    0.9501    0.5828    0.4398
    0.4057    0.5678           0.6458
   10.4103    0.6029    0.8704
  -30.2844    0.8439    0.6739
     648     -0.1708    0.9616

The member function read() associated with the TFList class goes 
through the text file and finds out the number of rows and columns.
It also checks for inconsistencies between different lines. 
Once the file is read into the anglist object, the dimension information
can be obtained from the

anglist.nrows
anglist.ncols 

fields of the object.  The floating point data in the file is stored in 
the array

anglist.data

which is allocated by the read() member function.

3. Writing SPARX Text File

In order to write a list of data to a SPARX Text File, one
must first create a TFList object with the right dimension and
format.

Suppose we have two one dimensional arrays we would like
to write out together in one file.

   x[0] = 2.3;
   x[1] = -6.3;
   x[2] = 1.2;
   
   y[0] = 2.0;
   y[1] = 0.0;
   y[2] = -3.2;

We can first create a TFList object (shlist) by 

   nrows  = 3;
   ncols  = 2;
   ndigit = -10;
   TFList shlist(nrows,ncols,ndigit);

   Note that ndigit is used here to specify the format in which each 
   entry will be written

   ndigit > 0, use floating point format, with at most ndigit digits
               after the decimal point
   ndigit = 0, use a free format without alignment
               (i.e. integers will be written out as whole numbers
               without decimal points)
   ndigit < 0, use fixed point format, with the ndigit as the total
               number of digits shown. 

4. Other utilities

   The TFList class contains a number of member functions that allows
   one to manipulate data entries after they are read in from a SPARX
   Text File and before they are written out.  The following are a list
   of utilities that have been implemented so far.

     Constructors: 
        TFList();
        TFList(int nr, int nc);
        TFList(int nr, int nc, int ndig);

     Destructor:
        ~TFList();

     Read function:
        int   read(char *filename);

     Write function:
        int   write(char *filename);

     Copy data into a TFList object
        void  Copy(float *rdata);
        void  CopyCol(int jcol, float *rcol);
        void  CopyRow(int jcol, float *rcol);

     Change values in a TFList object
        void  SetVal(int irow, int jcol, float val);

     Retrieve data from a TFList object
        float GetVal(int irow, int jcol);
*/

// constructor: initialize a TFList object
TFList::TFList()
{
   nrows  = 0;
   ncols  = 0;
   ndigit = 0;
   data   = NULL;
}

// constructor: initialize a TFList object of size nr by nc
TFList::TFList(int nr, int nc)
{
   int i;

   nrows  = nr;
   ncols  = nc;
   ndigit = 0;
   data   = (float*)SparxMalloc(nrows*ncols*sizeof(float));
   // initialize to zeros
   for (i=0;i<nrows*ncols;i++) data[i]=0.0;
}

// constructor: initialize a TFList object of size nr by nc. 
// When the list is written to a file, it will have the format specified 
// by ndig
TFList::TFList(int nr, int nc, int ndig)
{
   int i;

   nrows  = nr;
   ncols  = nc;
   ndigit = ndig;
   data   = (float*)SparxMalloc(nrows*ncols*sizeof(float));
   // initialize to zeros
   for (i=0;i<nrows*ncols;i++) data[i]=0.0;
}

// Destructor
TFList::~TFList()
{
   if (data) delete data;
}

#define data(i,j) data[((j)-1)*nrows + (i) - 1]

// read a text file into a TFList object
int TFList::read(char *filename)
{
   ifstream TFin1, TFin2;
// int  ntoken;
   char buffer[MAXLEN];
// int  i, j, imin, imax, lenbuf, ntokens, maxntks = 0, minntks = 100;
   int  i, j, imin=0, imax=0, ntokens, maxntks = 0, minntks = 100;
   string::size_type isemi;
   int  status = 0;
   
   TFin1.open(filename);

   // first pass to determine nrows and ncols
   while ( !TFin1.eof() ) {
      TFin1.getline(buffer, MAXLEN);
      // turn the buffer into a string
      string s(buffer);
      // check for SPIDER comments
      isemi = s.find_first_of(';');
      if ( isemi == string::npos ) {
         // not a SPIDER comment, check for empty line
         if ( s.length() >  0) {
            // Tokenize the buffer
            Tokenizer tk(s);
            nrows++;
            ntokens = tk.tokenCount();
            if (ntokens > maxntks) {
               maxntks = ntokens;
               imax    = nrows; 
            }
            if (ntokens < minntks) {
               minntks = ntokens;
               imin    = nrows;
            }
         } 
      }
   }
   TFin1.close();

   if (minntks != maxntks) {
      // something is wrong with the file
      // the number of entries on one line is different from 
      // that on another, take ncols = minntks.
      cerr << "TFList: error in the file: " << filename << endl;
      cerr << "TFList: line " << imin << " has " << minntks << " entries" 
           << endl;
      cerr << "TFList: line " << imax << " has " << maxntks << " entries" 
           << endl;
      status = 1;
   }
   ncols = minntks;

   // go through the file again and read the content into data 
   if (nrows > 0 && ncols >0) {
      if (data != NULL) delete data;
      data = (float*)SparxMalloc(ncols*nrows*sizeof(float));
      // initialize to zeros
      for (i=0;i<nrows*ncols;i++) data[i]=0.0;
      TFin2.open(filename);
      i = 1;
      while (i<=nrows) {
         TFin2.getline(buffer, MAXLEN);
         // turn the buffer into a string
         string s(buffer);
         // check for SPIDER comments
         isemi = s.find_first_of(';');
         if ( isemi == string::npos ) {
            // not a SPIDER comment, check for empty line
            if ( s.length() >  0) {
               Tokenizer tk(s);
               for (j=1; j <=ncols; j++) {
                  data(i,j) = (float)atof((tk.nextToken()).c_str());
               }
               i++;
            } 
         }
      } // end while
      TFin2.close();
   }
   else {
      status = -1;
   }
   return status;
}

// write a TFList object to a file
int TFList::write(char *filename)
{
    ofstream TFout;
    int i, j, status=0;
 
    TFout.open(filename);
    for (i = 1; i <=nrows; i++) {
       for (j = 1; j <=ncols; j++) {
           if (ndigit > 0) {
              // floating point format 
              TFout << scientific << setprecision(ndigit) << showpos
                    << data(i,j) << "   ";
           }
           else if (ndigit == 0) {
              // free format
              TFout << data(i,j) << "   ";
           }
           else {
              // fixed point format
              TFout << fixed << setw(-ndigit) << data(i,j) << "   ";
           }
       }
       TFout << endl;
    } 
    TFout.close();
    return status;
}

#define rdata(i,j) rdata[((j)-1)*nrows + (i) - 1]

// copy rdata into a TFlist object
void TFList::Copy(float *rdata)
{
   int i,j;

   for (j = 1; j<=ncols; j++)
      for (i = 1; i<=nrows; i++) 
         data(i,j) = rdata(i,j);
}
#undef rdata

// copy a column of data into a TFlist object
void TFList::CopyCol(int jcol, float *rdata)
{
   int i;

   for (i = 1; i<=nrows; i++) data(i,jcol) = rdata[i-1];
}

// copy a row of data into a TFlist object
void TFList::CopyRow(int irow, float *rdata)
{
   int j;

   for (j = 1; j<=ncols; j++) data(irow,j) = rdata[j-1];
}

// set a particular entry of a TFlist object to the supplied value.
void TFList::SetVal(int irow, int jcol, float val)
{
   data(irow,jcol) = val;
}


// extract a particular entry from a TFlist object
float TFList::GetVal(int irow, int jcol)
{
   if (irow >= 1 && irow <=nrows && jcol >=1 && jcol <=ncols) { 
      return data(irow,jcol);
   }
   else {
      cerr << "TFList::GetVal: index out of range!" << endl;
	return 0;
   }
}
#undef data
