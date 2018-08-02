/*
 * Copyright (c) 2005, Andrew Fernandes (andrew@fernandes.org);
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * - Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * - Neither the name of the North Carolina State University nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef PRECISION_H_INCLUDED
#define PRECISION_H_INCLUDED

#define SINGLE_PRECISION

#if !defined(SINGLE_PRECISION) && !defined(DOUBLE_PRECISION) && !defined(EXTENDED_PRECISION)
# define SINGLE_PRECISION
#endif

/* note that the 'integer_t' type must be signed, not unsigned */

#if defined(EXTENDED_PRECISION)

typedef int integer_t;
typedef long double real_t;
#define REAL_CONSTANT(x) x##l

#elif defined(DOUBLE_PRECISION)

typedef int integer_t;
typedef double real_t;
#define REAL_CONSTANT(x) x

#elif defined(SINGLE_PRECISION)

typedef int integer_t;
typedef float real_t;
#define REAL_CONSTANT(x) x##f

#endif /* ! *_PRECISION */

#endif /* ! PRECISION_H_INCLUDED */
