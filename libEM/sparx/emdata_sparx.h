/**
 * $Id$
 */
#ifndef emdata__sparx_h__
#define emdata__sparx_h__ 

public:
/** returns the fourier harmonic transform (FH) image of the current
 * image (in real space). The current image is not changed. The result is in
 * real/imaginary format. The FH switch is set on.
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @return the Fourier Harmonic image
 */
EMData *real2FH(float OverSamplekB);

/** returns the fourier version of the image 
 * from the FH version. The current image is not changed. The result is in
 * real/imaginary format. The FH switch is set off.
 * @param Size is the size of the image to be returned
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @param IntensityFlag=0 is the usual; =1 means that the input was an intensity
 * @return the shuffled version of the FFT
 */
EMData *FH2F(int Size, float OverSamplekB, int IntensityFlag =0);

/** returns the real version of the image 
 * from the FH version. The current image is not changed. The result is in
 * real format.
 * @param Size is the size of the image to be returned
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @param IntensityFlag=0 is the usual; =1 means that the input was an intensity
 * @return the real version of the data
 */
EMData *FH2Real(int Size, float OverSamplekB, int IntensityFlag =0);
		
		
/** Create a (1-D) rotationally averaged image.
 * @exception ImageDimensionException If 'this' image is not 2D.
 * @return 1-D rotationally-averaged image
 */					
EMData* rotavg();

/** Calculates the Center of Gravity 
 *  and the Radius of Gyration of the image.
 *  @returns the mass and the radius as vectors.
 *  
 */
 vector<float> cog();
 
#endif	//emdata__sparx_h__
