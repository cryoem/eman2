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

/** This file is a part of "emdata.h", to use functions in this file,
 * you should "#include "emdata.h",
 * NEVER directly include this file. */

#ifndef emdata__modular_h__
#define emdata__modular_h__

public:
/** Apply a processor with its parameters on this image.
 * @param processorname Processor Name.
 * @param params Processor parameters in a keyed dictionary.
 * @exception NotExistingObjectError If the processor doesn't exist.
 */
void process_inplace(const string & processorname, const Dict & params = Dict());

/** Call the process_inplace with an instance od Processor, usually this instancecan
 * be get by (in Python) Processors.get("name", {'k':v, 'k':v})
 * @param p the processor pointer
 * */
void process_inplace(Processor * p);

/** Apply a processor with its parameters on a copy of this image, return result
 * as a a new image. The returned image may or may not be the same size as this image.
 * @param processorname Processor Name.
 * @param params Processor parameters in a keyed dictionary.
 * @return the processed result, a new image
 * @exception NotExistingObjectError If the processor doesn't exist.
 * */
EMData * process(const string & processorname, const Dict & params = Dict()) const;

/** Call the process with an instance od Processor, usually this instance can
 * be get by (in Python) Processors.get("name", {'k':v, 'k':v})
 * @param p the processor pointer
 * */
EMData * process(Processor * p) const;

/** Compare this image with another image.
 * @param cmpname Comparison algorithm name.
 * @param with The image you want to compare to.
 * @param params Comparison parameters in a keyed dictionary.
 * @exception NotExistingObjectError If the comparison algorithm doesn't exist.
 * @return comparison score. The bigger, the better.
 */
float cmp(const string & cmpname, EMData * with, const Dict & params = Dict());

/** Align this image with another image and return the result image.
 *
 * @param aligner_name Alignment algorithm name.
 * @param to_img The image 'this' image aligns to.
 * @param params  Alignment algorithm parameters in a keyed dictionary.
 * @param cmp_name Comparison algorithm used in alignment.
 * @param cmp_params Parameter dictionary for comparison algorithm.
 * @exception NotExistingObjectError If the alignment algorithm doesn't exist.
 * @return The result image.
 */
EMData *align(const string & aligner_name, EMData * to_img,
			  const Dict & params = Dict(), const string & cmp_name = "",
			  const Dict& cmp_params = Dict());

/** Align this image with another image, return the parameters of the "n best" solutions
 * See Aligner::xform_align_nbest for more comments
 * @param aligner_name Alignment algorithm name.
 * @param to_img The image 'this' image aligns to.
 * @param params  Alignment algorithm parameters in a keyed dictionary.
 * @param nsoln the number of solutions you want to receive in the return vector.
 * @param cmp_name Comparison algorithm used in alignment.
 * @param cmp_params Parameter dictionary for comparison algorithm.
 * @exception NotExistingObjectError If the alignment algorithm doesn't exist.
 * @return an ordered vector of Dicts of length nsoln. The Dicts in the vector have keys "score" (i.e. correlation score) and "xform.align3d" (Transform containing the alignment)
 */
vector<Dict> xform_align_nbest(const string & aligner_name, EMData * to_img,
			  const Dict & params = Dict(), const unsigned int nsoln = 1, const string & cmp_name = "dot",
										 const Dict& cmp_params = Dict());

/** Calculate the projection of this image and return the result.
 * @param projector_name Projection algorithm name.
 * @param params Projection Algorithm parameters.
 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
 * @return The result image.
 */
EMData *project(const string & projector_name, const Dict & params = Dict());

/** Calculate the projection of this image and return the result.
 * @param projector_name Projection algorithm name.
 * @param t3d Transform object used to do projection.
 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
 * @return The result image.
 */
EMData *project(const string & projector_name, const Transform & t3d);

/** Calculate the backprojection of this image (stack) and return the result.
 * @param projector_name Projection algorithm name.
 * (Only "pawel" and "chao" have been implemented now).
 * @param params Projection Algorithm parameters.
 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
 * @return The result image.
 */
EMData *backproject(const string & projector_name, const Dict & params = Dict());

//#ifdef EMAN2_USING_CUDA
//vector<EMData*> cuda_project(const vector<Transform> transforms);
//#endif // EMAN2_USING_CUDA

#endif	//emdata__modular_h__
