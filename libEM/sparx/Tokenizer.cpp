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

#include "Tokenizer.h"

const string Tokenizer::WHITE=string("\r \n\t");

bool Tokenizer::moreToken() const
{   if ( position == target.length() )
         return false;
    if ( delimToken ) return true;
    int i0 = target.find_first_not_of
                 (delimiter, position);
    return(i0 != string::npos);
}

string Tokenizer::nextToken()
{   if ( delimToken ) 
         return nextAll();
    else
         return next();
}

string Tokenizer::next()
{   int i0=position;
    int i1=position;                  // start and end indices
    if ( i0 >= target.length() )
         return string();
    i0 = target.find_first_not_of(delimiter, i1);
    if ( i0 == string::npos ) 
         return string();
    i1 = target.find_first_of(delimiter, i0);
    if ( i1 == string::npos )
         i1 = target.length();
    position = i1;
    return string(target, i0, i1-i0);  // token found
}

string Tokenizer::nextAll()
{   int i0=position;
    int i1=position;                  // start and end indices
    if ( i0 >= target.length() )
         return string();
    if ( delimiter.find(target[i0]) != string::npos )
    {    i1 = i0+1;  } 
    else
    {    i1 = target.find_first_of(delimiter, i0);
         if ( i1 == string::npos )
              i1 = target.length();
    }
    position = i1;
    return string(target, i0, i1-i0);  // token found
}

int Tokenizer::tokenCount() const
{   if ( position == target.length() )
        return 0;
    int count = 0;
    int i0=position, i1=position; 
    while ( 1 )  
    {   if ( i1 == string::npos )
           return count;
        i0 = target.find_first_not_of(delimiter, i1);
        if ( i0 == string::npos )
           return (! delimToken ? count :
	          count + target.length() - i1);
        count++;
        if ( delimToken ) count += i0-i1;
        i1 = target.find_first_of(delimiter, i0);
    }
}
