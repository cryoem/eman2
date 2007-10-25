/**
 * $Id$
 */

/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

#ifndef normpadft_h__
#define normpadft_h__ 

public:

/** Multiply a real image by (-1)**(ix+iy+iz) to center
 *  the fft version.
 *
 */
void center_origin();

void center_origin_yz();

/** Multiply a Fourier image by (-1)**(ix+iy+iz) to center it.
 *
 */
void center_origin_fft();


/** return an image object that has been padded with zeros to
 *  an integral (npad) factor (e.g., npad=2 means the new 
 *  image will be 2x larger in each direction).  The default
 *  is to pad 4x.
 * The current image is not changed.
 *
 * @return An image object that has been padded npad-times.
 */
EMData *zeropad_ntimes(int npad=4);

EMData *zeropad_ntimes_fft_shuffle(int npad=4);

/** return an image object that has been fft-padded/unpadded.
 * The current image is not changed.
 *
 * @return An image object that has been fft-padded/unpadded.
 */
EMData *pad_fft(int npad = 1);


/** Remove padding, leaving a single corner of the image.
 *  The current image is changed in place.
 *
 *  The assumption is that after an in-place inverse fft
 *  the real-space image contains too much information 
 *  because it may have been zero-padded some integer factor
 *  of times and it has also been extended slightly along x
 *  for the fft.  Here we keep only the data corresponding 
 *  to ix=0,...,nxold-1, iy=0,...,nyold-1, iz=0,...,nzold-1,
 *  where nxold, nyold, nzold are the sizes of the original
 *  image.
 */
void postift_depad_corner_inplace();


/** Interpolate up image by padding with zeroes its Fourier transform.
 *  
 *  @param[in] nxni new x size (has to be larger/equal than the original x size)
 *  @param[in] nyni new y size (has to be larger/equal than the original y size)
 *  @param[in] nzni new z size (has to be larger/equal than the original z size)
 *  @param RetReal
 *  
 *  @return New interpolated up image.
 */
EMData* FourInterpol(int nxni, int nyni=0, int nzni=0, bool RetReal = true);

/** Truncate Fourier transform of an image, it will reduce its size.  (It is a form of decimation).
 *  
 *  @param[in] nxni new x size (has to be larger/equal than the original x size)
 *  @param[in] nyni new y size (has to be larger/equal than the original y size)
 *  @param[in] nzni new z size (has to be larger/equal than the original z size)
 *  @param RetReal
 *  
 *  @return New truncated up image.
 */
EMData* FourTruncate(int nxni, int nyni=0, int nzni=0, bool RetReal = true);
EMData* FourInterpol_i(int nxni, int nyni=0, int nzni=0, bool RetReal = true);
EMData* Four_ds(int nxni, int nyni=0, int nzni=0, bool RetReal = true);
EMData* Four_shuf_ds_cen_us(int nxni, int nyni=0, int nzni=0, bool RetReal = true);


EMData* filter_by_image(EMData* image, bool RetReal = true);

#endif	//normpadft_h__
