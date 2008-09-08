/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
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
 * */

#ifndef eman__byteorder_h__
#define eman__byteorder_h__ 1

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

#include <cstddef>

namespace EMAN
{
	/** ByteOrder defines functions to work on big/little endian byte orders.
     *
     * The byte order is the order in which bytes are stored to
     * create larger data types such as the int and long values.
     * Different kinds of computers use different byte order conventions.
     * 
     * There are 2 major byte orders: big-endian and little-endian.
     *
     * big-endian (like SGI) store the most significant bytes (i.e. the
     * bytes that hold the largest part of the value) first. little-endian
     * (like x86) store the most significant byte last. 
     *
     * The host byte order is the byte order used on the current host.
     *
     * ByteOrder class defines functions for
     *    checking running-host byte-order,
     *    checking data byte-order,
     *    convert data from one byte-order to another byte-order.
     */
	class ByteOrder
	{
	public:
		static bool is_host_big_endian();

		/** given a small floating number, return whether the number
		 * is in big endian or not. If a number is smaller than 65535,
		 * it is defined as a "small" number here.
		 */
		static bool is_float_big_endian(float small_number);
		
		
		/** given a pointer to a reasonable small integer number,
		 * return whether the number is big endian or not.
		 * For a n-bit integer, the number should < (2 ^ (n/2)).
		 * e.g., for 'int', number should < 65535;
		 * for 'short', number should < 255.
		 */
		template < class T > static bool is_data_big_endian(const T * small_num_addr)
		{
			if (!small_num_addr) {
				return false;
			}

			bool data_big_endian = false;
			size_t typesize = sizeof(T);
			char *d = (char *) small_num_addr;

			if (is_host_big_endian()) {
				data_big_endian = false;
				for (size_t i = typesize / 2; i < typesize; i++)
					{
						if (d[i] != 0) {
							data_big_endian = true;
							break;
						}
					}
			}
			else {
				data_big_endian = true;
				for (size_t j = 0; j < typesize / 2; j++) {
					if (d[j] != 0) {
						data_big_endian = false;
						break;
					}
				}
			}

			return data_big_endian;
		}

		/** convert data from host byte order to big endian order.
		 * 'n' is the number of elements of type T.
		 */
		template < class T > static void become_big_endian(T * data, size_t n = 1)
		{
			if (!is_host_big_endian()) {
				swap_bytes < T > (data, n);
			}
		}

		/** convert data from host byte order to little endian order.
		 * 'n' is the number of elements of type T.
		 */
		template < class T > static void become_little_endian(T * data, size_t n = 1)
		{
			if (is_host_big_endian()) {
				swap_bytes < T > (data, n);
			}
		}

		/** swap the byte order of data with 'n' T-type elements.
		 */
		template < class T > static void swap_bytes(T * data, size_t n = 1)
		{
			unsigned char s;
			size_t p = sizeof(T);
			char *d = (char *) data;

			if (p > 1) {
				for (size_t i = 0; i < n; i++, d += p) {
					for (size_t j = 0; j < p / 2; j++) {
						s = d[j];
						d[j] = d[p - 1 - j];
						d[p - 1 - j] = s;
					}
				}
			}
		}

	private:
		static bool is_host_endian_checked;
		static bool host_big_endian;
	};

}

#endif
