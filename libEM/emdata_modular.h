/**
 * $Id$
 */

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
		
/** Apply a processor with its parameters on a copy of this image, return result 
 * as a a new image. The returned image may or may not be the same size as this image.
 * @param processorname Processor Name.
 * @param params Processor parameters in a keyed dictionary.
 * @return the processed result, a new image 
 * @exception NotExistingObjectError If the processor doesn't exist.
 * */
EMData * process(const string & processorname, const Dict & params = Dict());

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
 * @param comp_name Comparison algorithm used in alignment.
 * @param cmp_params Parameter dictionary for comparison algorithm.
 * @exception NotExistingObjectError If the alignment algorithm doesn't exist.
 * @return The result image.
 */
EMData *align(const string & aligner_name, EMData * to_img,
			  const Dict & params = Dict(), const string & comp_name = "", 
			  const Dict& cmp_params = Dict());

/** Calculate the projection of this image and return the result.
 * @param projector_name Projection algorithm name.
 * @param params Projection Algorithm parameters.
 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
 * @return The result image.
 */
EMData *project(const string & projector_name, const Dict & params = Dict());

/** Calculate the backprojection of this image (stack) and return the result.
 * @param projector_name Projection algorithm name. 
 * (Only "pawel" and "chao" have been implemented now). 
 * @param params Projection Algorithm parameters.
 * @exception NotExistingObjectError If the projection algorithm doesn't exist.
 * @return The result image.
 */
EMData *backproject(const string & projector_name, const Dict & params = Dict());


#endif	//emdata__modular_h__

