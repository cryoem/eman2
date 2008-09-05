#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


#from EMAN2 import *
from EMAN2_cppwrap import *
from global_def import *
#from filter import *
#from fundamentals import *
#from math import pi, exp, sqrt
#from random import *
#from sparx import *
#from sys import maxint
#from utilities import  *
#from sideprojectionmatchB import *

#import pylab 
'''
def side_projection_match(DataSetName,NumImgs,rAngs,SymString,NumLoops=30):
	"""
		August 11th 2006
		Read  in a 
			DataSetName - name of a stack of class sums along with ext;
			NumImgs - the number of slots within this range,
			rAngs  - a range of azimuthal angles to search (given by a single value),
			SymString - like 'd2' or 'c3'
			NumLoops - the number of rounds of projection matching
			
	 	later (possibly) will use a starting volume
 
		returns the best assignment in angles to the class sums and to file SP_DataSetAngles.txt
		writes out SPCurrentVol?.mrc ; SPReProj?.spi , SPEvenlySpaced?.spi after each round 
 
 		Try using BC.spi for DataSetName which is Beate Rockel's test data: 30 images 25 by 25 pixels,
		            rAngs should be 180. BC.spi should be in sparx/test
		If you want to use a good starting model, replace 'assigned' in line 97 by 'slots'
	"""

#        0.     Initialize everything

	eogMode=0;
	os.system('mkdir SPpng &');
	os.system('mkdir SPfiles &');
    

        DataSetBase=DataSetName[0:(len(DataSetName)-4)];
        DataSetAnglesName='SPfiles/SP_'+DataSetBase+'Angles.txt';
	DataSetHDFName = DataSetBase + '.hdf';

	for j in range(NumImgs):
		ClassSums = EMData();


	for j in range(NumImgs):
		ClassSums = EMData();
		ClassSums.read_image(DataSetName,j);
		ClassSums.write_image(DataSetHDFName,j);


#	Vol3D=EMData();
#	Vol3D.read_image(Vol3DName);
	
	ClassSums = EMData();
	ClassSums.read_image(DataSetHDFName,1);

	oneSide = ClassSums.get_xsize();
	sizeCl  =  oneSide*oneSide;

#	last=NumImgs-1;
#        slots[0]+slots[last]=rAngs
#        spacing = (slots[last]-slots[0])/last;
#                = slots[0]*2

	NumSlots=NumImgs;
	slots0 =rAngs/(2.0*NumSlots);
	spacing=rAngs/(float(NumSlots));

	slots= [ (j *spacing+slots0) for j in range(NumSlots)];

	CCCtemp = [0 for j in range(NumSlots*NumImgs)];
	CCC     = [0 for j in range(NumSlots*NumImgs)];


#
#      1.   Create a reconstruction from Random Assignments

	assigned=range(NumImgs);
	taken=range(NumImgs);

	for jj in assigned:
		taken[jj]=0;
		assigned[jj]=-1;
	
	TotalTaken=0;
	
	while TotalTaken< NumImgs:
			qq=randrange(NumImgs);
			if taken[qq]==0:
				taken[qq]=1;
				assigned[qq]=slots[TotalTaken];
				TotalTaken=TotalTaken+1;

	angleList=[[0 for k in range(3)] for j in range(NumImgs)];
	for AImage in range(NumImgs):
		CS=EMData();
        	CS.read_image(DataSetHDFName,AImage);
		newphi=slots[AImage];                      # TOGGLE assigned and slots
		CS.set_attr_dict({'phi':newphi, 'theta':90, 'psi':0, 'sx':0, 'sy':0}); 
		CS.write_image(DataSetHDFName, AImage);
		angleList[AImage][0]=newphi;
			
 # Now assigned gives an angle for each projection

	vv=range(NumImgs);
	
#    print(DataSetHDFName);
#     Vol1=recons3d_4nn(DataSetHDFName, vv,SymString);
	Vol1=recons3d_wbp(DataSetHDFName, vv,angleList,'exaxct',SymString);
	Vol1.write_image('SPfiles/RandomVol.mrc');


	Vol1ft,kb=prep_vol(Vol1);

### end of 1;  Vol1 is now a preliminary reconstruction from random assignments
#


#               Loop Starts
#      A.      Create a CCC table


	for Loop in range(NumLoops):   
	    
		print('Round '+str(Loop));
		for iphi in range(NumSlots):
			phi= slots[iphi];
			proji  = prgs(Vol1ft,kb,[phi,90,0,0,0]);
#		proji.write_image('projTemp.hdf',iphi);
			for AImage in range(NumImgs):
        			CS=EMData();
        			CS.read_image(DataSetHDFName,AImage);
				CCCtemp[AImage+iphi*NumImgs] =ccc(CS,proji);

#  CCCtemp[AImage+iphi*NumImgs]  is the CCC of AImage with a reweprojection at slots[iphi]

#      B.      even out the table

                rightW= range(NumSlots);
                leftW= range(NumImgs);
		for j in range(10):
                      for AImage in range(NumImgs):
		          temp=0;
		          for iphi in range(NumSlots):
			      temp=temp+CCCtemp[ AImage+iphi*NumImgs]*rightW[iphi];
                          leftW[AImage]= 1.0/temp;
			  
		      sLW=sum(leftW);
		      leftW = [(NumSlots)*leftW[AImage]/sLW for AImage in range(NumImgs)];
                      for iphi in range(NumSlots):
		          temp=0;
		          for AImage in range(NumImgs):
		              temp=temp+ leftW[AImage]* CCCtemp[AImage+iphi*NumImgs] ;
                          rightW[iphi]= 1.0/temp;
			  
		      sRW=sum(rightW);
		      rightW = [NumSlots*rightW[iphi]/sRW for iphi in range(NumSlots)];
		       
		for iphi in range(NumSlots) :
			for AImage in range(NumImgs):
			    CCC[AImage+iphi*NumImgs] = leftW[AImage]* CCCtemp[AImage+iphi*NumImgs] *rightW[iphi];


                write(CCCtemp,NumImgs,NumSlots,'SPfiles/SPCCCtempimage.img');
                write(CCC,NumImgs,NumSlots,'SPfiles/SPCCCimage.img');
                write(rightW,NumSlots,1,'SPfiles/leftW.txt');
                write(leftW,NumImgs,1,'SPfiles/rightW.txt');
                write(CCC,NumImgs,1,'SPfiles/SPCCC.txt');
                write(CCCtemp,NumImgs,1,'SPfiles/SPCCCtemp.txt');

#        C.   Look through the  CCC table; find best assignment

		BestAngs =[-1 for j in range(NumImgs)];
		for AImage in range(NumImgs):
			vv=[ CCCtemp[AImage+iphi*NumImgs] for iphi in range( NumSlots)  ]; # TOGGLE CCC and CCCtemp
			vvmax=max(vv);
			for iphi in range(NumSlots):
				if  (vv[iphi]>.9999999*vvmax):
					iphibest=iphi;
			BestAngs[AImage]=slots[iphibest];
			CS=EMData();
	        	CS.read_image(DataSetHDFName,AImage);
			CS.set_attr_dict({'phi':BestAngs[AImage], 'theta':90, 'psi':0, 'sx':0, 'sy':0});
			CS.write_image(DataSetHDFName, AImage);
			 

#        D.  Do  another Reconstruction, Reprojection
				

                try: 
                    outAngles= open(DataSetAnglesName); os.remove(DataSetAnglesName); 
                except IOError: print('The File '+DataSetAnglesName + ' Does Not Exist');
                outAngles = open(DataSetAnglesName, "w"); #  write angles out to txt file
		angleList=[[0 for k in range(3)] for j in range(NumImgs)];

			

                for AImage in range(NumImgs):  # 
                    outAngles.write("%d\t%f\n" % (AImage+1,BestAngs[AImage])) 
                    angleList[AImage][0]=BestAngs[AImage];  # to be used with wbp
                outAngles.close(); 
              

		if Loop>0:
                    try:
                        fsock=open(VolString); os.remove(VolString);
                    except  IOError: Good=1;
                    try:
                        fsock=open(ReProjString); os.remove(ReProjString);
                    except  IOError: Good=1;
                    try:
                        fsock=open(EvenlySpacedString); os.remove(EvenlySpacedString);
                    except  IOError: Good=1;

#        D (cont)  Reconstruction

		VolString='SPfiles/SPCurrentVol'+str(Loop)+'.mrc'
                try:
                    fsock=open(VolString); os.remove(VolString);
                except  IOError:  Good=1;

		vv=range(NumImgs);
#	  Vol1=recons3d_4nn(DataSetHDFName, vv,'d2');
                Vol1=recons3d_wbp(DataSetHDFName, vv,angleList,'exact',SymString); # TOGGLE reconstruction algorithms

		fNm1 = NumImgs-1.0;	 
                BestAngsSorted=sorted(BestAngs);
		spacingVec= [ BestAngsSorted[j+1] - BestAngsSorted[j] for j in range(NumImgs-1)];
		sigmaspacingVec = [spacingVec[j]*spacingVec[j] for j in range(NumImgs-1)];
		AngsToAverage= sqrt(sum(sigmaspacingVec)/fNm1);
		AngsPlusMin=0; #int(AngsToAverage*sqrt(2.0) );  TOGGLE
		weightVec  = [exp(-(j-1.0 *AngsPlusMin)*(j-1.0 *AngsPlusMin)/(2.0* AngsToAverage*AngsToAverage)) for j in range(2*AngsPlusMin+1)];
		sweightVec = sum(weightVec);
		lweightVec=len(weightVec);
		weightVec = [ weightVec[j]/sweightVec for j in range(lweightVec)];
		VolSum=weightVec[AngsPlusMin]* Vol1.copy();
		print(weightVec);
		
		for j in range(lweightVec):
		    RA=Transform3D(0,0, j- 1.0*AngsPlusMin);
                    if (j!= AngsPlusMin):
		        Vol2add= Vol1.rot_scale_trans(RA);
		        VolSum+=weightVec[j]* Vol2add;

		VolSum.write_image(VolString)
		VolSumft,kb=prep_vol(VolSum);



#        D (cont)  Reprojection

		ReProjString='SPfiles/SPReProj'+str(Loop)+'.spi';
		EvenlySpacedString =  'SPfiles/SPEvenlySpaced' +str(Loop) +'.spi';
                try:
                    fsock=open(ReProjString); os.remove(ReProjString);
                except  IOError: Good=1;
                try:
                    fsock=open(EvenlySpacedString); os.remove(EvenlySpacedString);
                except  IOError:  Good=1;


		for AImage in range(NumImgs):
			phi = BestAngs[AImage];
			proji  = prgs(VolSumft,kb,[phi,90,0,0,0]);
			proji.write_image(ReProjString,AImage);



		for iphi in range(NumSlots):
			phi = slots[iphi];
			proji  = prgs(VolSumft,kb,[phi,90,0,0,0]);
			proji.write_image(EvenlySpacedString, iphi);



			
#       5. This is the end of the big loop			
#        SPCurrentVol?.mrc stores the current volume
#        SPReproj?.spi stores the reprojections
#        SPEvenlySpaced?.spi stores an even sampling of the projections			
	return BestAngs, AngsToAverage






def Sidewinder(DataSetName,NumImgs,SymString, NumLoops=20,MovesPerLoop=10000):
    """
	August 11th 2006
	Read  in a 
		DataSetName - name of a stack of class sums along with ext,
		              the symmetry axis should lie along y (with v2)
		NumImgs - the number of slots within this range,
		SymString - like 'd2' or 'c3'
		NumLoops - the number of rounds of projection matching
		MovesPerLoops - the number of rounds of projection matching
		
 	later (possibly) will use a starting volume

	returns the best assignment in angles to the class sums

		Try using BC.spi for DataSetName which is Beate Rockel's test data: 30 images 25 by 25 pixels,
    """


#     0. Read in the data from DataSetName
#
    eogMode=0;
    os.system('mkdir SWpng &');
    os.system('mkdir SWfiles &');
    
    DSymString = SymString[0]; DSym=-1;
    if (DSymString=='d'): DSym=1;
    if (DSymString=='c'): DSym=0;
    if (DSym==-1):
        print('Need to enter a symmetry like d2 or c3'); return;
	  
    nCSym = int(SymString[1]); 
    
    MeanName='SWfiles/SWMean.spi'; MaskName='SWfiles/SWMask.spi';

    try:
        fsock=open(MeanName); os.remove(MeanName);
    except  IOError: print('The File '+MeanName + ' Does Not Exist');
    try:
        fsock=open(MaskName); os.remove(MaskName);
    except  IOError: print('The File '+MaskName + ' Does Not Exist');
    	  
    A = EMData();
    A.read_image(DataSetName,1);
#    A.print_image();

    display(A); 
    
    OrientationCorrect = raw_input('Is the symmetry axis up and down within the v2 window?');
    
    if ((OrientationCorrect[0]=='n') | (OrientationCorrect[0]=='N')):   #
    	print('Rotate images by 90 degrees before using this program, please.')
    	return;
    

    oneSide=A.get_xsize();
    sizeA= oneSide*oneSide;
   
    ImageData=  EMData()
    ImageData.set_size(oneSide,oneSide,NumImgs);
    
    oddC = nCSym%2   ; # 0 if even , 1 if odd; 
    oddCp1= oddC + 1;

#    print(oneSide);
   
    for AImage in range(NumImgs): # NumImgs
        A=EMData();
        A.read_image(DataSetName,AImage);
        oneSide=A.get_xsize();
        for ix in range(oneSide):
            for iy in range(oneSide):
                tt=A.get_value_at(ix,iy);
	        ImageData.set_value_at(ix,iy,AImage,tt);
	      

#     The Data has been read in 
# -------------------------------------------------------
#     1. Create the Mask    	  
    LowerX=2 ;   UpperX=oneSide+1-LowerX;   #LowerX=25 for Rockel, LowerX=15 for HRS
    LowerY=3 ;   UpperY=oneSide+1-LowerY;  #  LowerY= 7 for Rockel, LowerY=5 for HRS

    Mask=EMData();  
    Mask.set_size(oneSide,oneSide);
     
    for ix in range(LowerX,UpperX):
        for iy in range(LowerY,UpperY):
            Mask.set_value_at(ix,iy,1);
	    
    Mask.write_image(MaskName);
    
     
#
#   Now the mask has been made
# -------------------------------------------------------
#  2.  Look at the mean value along each fixed line of jk constant
#            subtract off the mean

    print('Entering Section 2');

    vv=EMData();
    vv.set_size(oneSide);
    C=EMData();
    C.set_size(oneSide,oneSide);
    ImageDataB=ImageData.copy();
 
    A=EMData();    
    B=EMData();

    

    for AImage in range(NumImgs):
        for ix in range(oneSide):
            for iy in range(oneSide):
                C.set_value_at(ix,iy, ImageData.get_value_at(ix,iy,AImage) ) ;

        A=C.copy() ; 
        B=C.copy() ;
    
        for jk in range(oneSide):   # do sweeps over y
            sumI = 0   ;
            sumvv= 0;
            for jj in range(oneSide): # do sweeps over x, 
                if ( (A.get_value_at(jj,jk)  !=0) &  (Mask.get_value_at(jj,jk)!=0)):
                    sumI=sumI+1;
                    sumvv=sumvv+A.get_value_at(jj,jk);

            if (sumI>0):
                vv.set_value_at(jk,sumvv/sumI);  # this is the average value along the line
            else:
                vv.set_value_at(jk,0);
        
            for jj in range(oneSide):
                B.set_value_at(jj,jk,0);    

                if( (A.get_value_at(jj,jk) !=0 ) & ( Mask.get_value_at(jj,jk)!=0 )  ):
                    B.set_value_at(jj,jk, A.get_value_at(jj,jk)-vv.get_value_at(jk));

        for ix in range(oneSide):
            for iy in range(oneSide):
                ImageDataB.set_value_at(ix,iy,AImage, B.get_value_at(ix,iy));


#clear Afn AImage jj vv jk sumI sumvv B C% leave ImageDAta
#clear indexA indexAR indexAx indexAy Ax Ar Ay BImage

    D=EMData();
    D.set_size(oneSide,oneSide,1);
    D.to_zero();

    for AImage in range(NumImgs):   #   D = mean(ImageDataB,3) 
        for ix in range(oneSide):
            for iy in range(oneSide):
                tempD = D.get_value_at(ix,iy);
                tempB = ImageDataB.get_value_at(ix,iy,AImage);
                D.set_value_at(ix,iy, tempD+(tempB/NumImgs));
            
    D.write_image(MeanName);
#    D.print_image(); 	    	

#
# -------------------------------------------------------
#  3.  The mean has been subtracted off
#  Symmetrize  the ImageData
#    Ar mean 180 rotation the r-axis
#    Ax mean mirror across the x-axis
#    Ay mean mirror across the y-axis, which was the original symmetry axis
# 
#    If A is in 0..30, then Ar is in 60..30, Ay is in 60..90 
#         for indexA  = 1:NumImgs
    print('Entering Section 3');


    A =EMData();     Ar=EMData();
    Ay=EMData();     Ax=EMData();

    A.set_size(oneSide,oneSide);     Ar.set_size(oneSide,oneSide);
    Ay.set_size(oneSide,oneSide);    Ax.set_size(oneSide,oneSide);
  
    ImageDataGuessed = [0 for  j in range(oneSide*oneSide*NumImgs* oddCp1*2) ];
    
     
    for indexA  in range(NumImgs):
        A.to_zero();         Ar.to_zero();
        Ay.to_zero();        Ax.to_zero();
	
        indexAr = 2*NumImgs  -1 - indexA;
        indexAy = 2*NumImgs  +   indexA;
        indexAx = 2*oddCp1*NumImgs  -1 - indexA;

	for ix in range(oneSide):
            for iy in range(oneSide):
                A.set_value_at(ix,iy, ImageDataB.get_value_at(ix,iy,indexA));
                ImageDataGuessed[ix+iy*oneSide +indexA*oneSide*oneSide]= A(ix,iy);
		
        for ix in range(oneSide):
            for iy in range(oneSide):
                Ax.set_value_at(ix,iy, A(ix,oneSide-1-iy)); 
                Ay.set_value_at(ix,iy, A(oneSide-1-ix,iy));
                Ar.set_value_at(ix,iy, A(oneSide-1-ix,oneSide-1-iy));

        for ix in range(oneSide):
            for iy in range(oneSide):
                if (oddC & DSym): ImageDataGuessed[ix+iy*oneSide + indexAr*oneSide*oneSide]= Ar(ix,iy); 
                if (oddC):        ImageDataGuessed[ix+iy*oneSide + indexAy*oneSide*oneSide]= Ay(ix,iy); 
                if (DSym):        ImageDataGuessed[ix+iy*oneSide + indexAx*oneSide*oneSide]= Ax(ix,iy); 


    ImageNum=1;
    vv1=[[ImageDataGuessed[ix+iy*oneSide + ImageNum*oneSide*oneSide] for ix in range(oneSide)] for iy in range(oneSide)];
    vv2=[[ImageDataGuessed[ix+iy*oneSide + (29-ImageNum)*ImageNum*oneSide*oneSide ] for ix in range(oneSide)] for iy in range(oneSide)];

#    pylab.subplot(2,1,1);
    pylab.imshow(vv1,interpolation='nearest');
    pylab.title('Image Number= '+str(ImageNum),name='Arial',size=14,weight='bold');
    pylab.draw();    pylab.gray();
    pylab.savefig('SWpng/vv1.png')
    if (eogMode):      os.system('eog SWpng/vv1.png &');
#    pylab.subplot(2,1,2);
    pylab.imshow(vv2,interpolation='nearest');
    pylab.title('Image Number= '+str(29-ImageNum),name='Arial',size=14,weight='bold');
    pylab.draw();    pylab.gray();
    pylab.savefig('SWpng/vv2.png')
    if (eogMode):      os.system('eog SWpng/vv2.png &');


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.  Here we take the CCC of every image with each 
#    other (include symmetrized versions)
#

    print('Entering Section 4');
    NE= 2*oddCp1*NumImgs;

    CCCtoSolve = [[0  for j1 in range(NE)] for j2 in range(NE)];
    AList = [0  for j12 in range(oneSide*oneSide)];

    for AImage in range(NE):
        normA=0;
	print(AImage);
        AA=ImageDataGuessed[AImage*oneSide*oneSide:(( AImage+1)*oneSide*oneSide)];
        for ixy in range(oneSide*oneSide):
            normA = normA + AA[ixy]*AA[ixy];
         # image A and norm have been created		
        for BImage in range(NE) :
            normB=0;
	    tempSum=0;
	    BB=ImageDataGuessed[BImage*oneSide*oneSide:((BImage+1)*oneSide*oneSide)];
            for ixy in range(oneSide*oneSide):
                tempB =  BB[ixy];
                tempSum  = tempSum + AA[ixy]*tempB; 
                normB = normB + tempB*tempB;
            # image B and norm have been created    

            CCCtoSolve[AImage][BImage] =  tempSum/(sqrt(normA*normB)); 
	    

#    CCCtoSolve
    pylab.hold(False);
    CCCtoSolvePlot= pylab.imshow(CCCtoSolve,interpolation='nearest');
    pylab.title('Sec4CCCplot',name='Arial',size=14,weight='bold');
    pylab.draw();  pylab.gray();
    pylab.savefig('SWpng/CCCplot.png')
    if (eogMode):     os.system('eog SWpng/CCCplot.png &');

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   5.  Now I find the shape of the average cross-correlation 
#      (before removing some)
#

    CCCtoSolveSorted    = []; 
    CCCtoSolveSortedAve = []; 

    for AImage in range(NE):
        vv = [CCCtoSolve[AImage][BImage] for  BImage in range(NE)];
	CCCtoSolveSorted.append(sorted(vv));

    for BImage in range(NE):
	ww= [CCCtoSolveSorted[AImage][BImage] for  AImage in range(NE)];
	CCCtoSolveSortedAve.append( sum(ww)/float(NE));

    pylab.plot(CCCtoSolveSortedAve);
    pylab.title('This is CCCtoSolveSortedAve');
    pylab.draw();
    pylab.savefig('SWpng/CCCtoSolveSortedAve.png')
    if (eogMode):     os.system('eog SWpng/CCCtoSolveSortedAve.png &'); pylab.hold(False)

    
    print('This is CCCtoSolveSortedAve');    print(CCCtoSolveSortedAve);
	
    for j in range(NE):
        vv=CCCtoSolveSorted[j];
        pylab.plot(vv);
	pylab.hold(True);

    pylab.title('This is CCCtoSolveSorted');
    pylab.draw();
    pylab.savefig('SWpng/CCCtoSolveSorted.png')
    if (eogMode):     os.system('eog SWpng/CCCtoSolveSorted.png &'); 
    pylab.hold(False);


    CCCtoSolveVar =[];
    for AImage in range(NE):
	ww=sorted(CCCtoSolveSorted[AImage]);
        DiffVec=[ pow(ww[j2] - CCCtoSolveSortedAve[j2] ,2) for j2 in range(NE)];
        CCCtoSolveVar.append(sqrt(sum(DiffVec)/len(DiffVec) ));


    print('This is CCCtoSolveVar');    print(CCCtoSolveVar);



# CCCtoSolveSortedAve is the mean sorted cross-correlation
# CCCtoSolveVar is the variance vector

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   6.    Use the Variance to discard about half: Keep asking until its done correctly
#   
    DecidedQuest ='n';
        
    while (DecidedQuest[0] != 'y'):
    
        cutoff2Keep =  raw_input('Enter a cutoff keep changing it till there is a tight fit: ');
    
        rangeCCCtoSolve =  find(CCCtoSolveVar,'lt',cutoff2Keep) ;
        NEReduced= len(rangeCCCtoSolve);
        NumImgsRed = int(NEReduced/(2*oddCp1)); 
        print('NEReduced  = ' +  str(NEReduced) );
        print('NumImgsRed = ' +  str(NumImgsRed));
        CCCtoSolveSortedTrue    = [];
        CCCtoSolveSortedTrueAve = [];
    
    
    #  
        for AImage in rangeCCCtoSolve:
             ww = [CCCtoSolveSorted[AImage][BImage] for BImage in range(NE)];
             CCCtoSolveSortedTrue.append([ ww[BImage] for BImage in range(len(ww))]);  #increasing as we 
                                                                                            # go down the column
    
        for BImage in rangeCCCtoSolve:
             ww = [CCCtoSolveSortedTrue[AImage][BImage] for AImage in rangeCCCtoSolve];
             CCCtoSolveSortedTrueAve.append(  sum(ww)/len(ww)  );
    
        for j in range(NumImgsRed):
            vv=CCCtoSolveSortedTrue[j];
            pylab.plot(vv);
    	    pylab.hold(True);
    
        pylab.title('CCCtoSolveSortedTrue: NumImgsRed='+str(NumImgsRed)+'from NumImgs'+str(NumImgs) );
        pylab.draw();
        pylab.savefig('SWpng/CCCtoSolveSortedTrue.png')
        os.system('eog SWpng/CCCtoSolveSortedTrue.png &'); 
        pylab.hold(False);
    
        print('rangeCCCtoSolve');   print(rangeCCCtoSolve);
        
        print('This is CCCtoSolveSortedTrueAve');    print(CCCtoSolveSortedTrueAve);
	
        DecidedQuest =  raw_input('Happy with this? If so type yes ');
    

############
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  7.     We now need to construct the reduced 
#             matrix from CCCtoSolveReduced

    CCCtoSolveReduced = [[0 for j1 in range(NEReduced)] for j2 in range(NEReduced)];

    for jj in range(NumImgsRed):
        jjo = jj;                            jjMappedo = rangeCCCtoSolve[jj] ; # remember NE is 4 or 2 times NumImgs
        jjr = 2*NumImgsRed-1-jj;             jjMappedr = 2*NumImgs - 1 - jjMappedo ;
        jjy = 2*NumImgsRed+jj;               jjMappedy = 2*NumImgs +     jjMappedo ;
        jjx = 2*oddCp1*NumImgsRed-1-jj;      jjMappedx = 2*oddCp1*NumImgs   - 1 - jjMappedo ;
        for jk in range(NumImgsRed):
            jko = jk;                          jkMappedo = rangeCCCtoSolve[jk];
            jkr = 2*NumImgsRed -1-jk;          jkMappedr = 2*NumImgs   - 1 - jkMappedo ;
            jky = 2*NumImgsRed   +jk;          jkMappedy = 2*NumImgs   +     jkMappedo;
            jkx = 2*oddCp1*NumImgsRed -1-jk;   jkMappedx = 2*oddCp1*NumImgs  -1 - jkMappedo ;
            
            temp = CCCtoSolve[jjMappedo][jkMappedo];
            CCCtoSolveReduced[jjo][jko] = temp;            
	    CCCtoSolveReduced[jjx][jkx] = temp;
            if oddC: 
	        CCCtoSolveReduced[jjy][jky] = temp;
	        CCCtoSolveReduced[jjr][jkr] =temp; 
            temp = CCCtoSolve[jjMappedo][jkMappedx];
            CCCtoSolveReduced[jjo][jkx] = temp; 
	    CCCtoSolveReduced[jjx][jko] =temp;
            if oddC: 
	        CCCtoSolveReduced[jjy][jkr] = temp; CCCtoSolveReduced[jjr][jky] =temp;
                temp = CCCtoSolve[jjMappedo][jkMappedy];
                CCCtoSolveReduced[jjo][jky] = temp; CCCtoSolveReduced[jjx][jkr] =temp;
                CCCtoSolveReduced[jjy][jko] = temp; CCCtoSolveReduced[jjr][jkx] =temp; 
                temp = CCCtoSolve[jjMappedo][jkMappedr];
                CCCtoSolveReduced[jjo][jkr] = temp; CCCtoSolveReduced[jjx][jky] =temp;
                CCCtoSolveReduced[jjy][jkx] = temp; CCCtoSolveReduced[jjr][jko] =temp;


    CCCtoSolveReducedSorted    = []; 
    CCCtoSolveReducedSortedAve = []; 


    for AImage in range(NEReduced):
        vv = [CCCtoSolveReduced[AImage][BImage] for  BImage in range(NEReduced)];
	CCCtoSolveReducedSorted.append(sorted(vv));


    for BImage in range(NEReduced):
	ww= [CCCtoSolveReducedSorted[AImage][BImage] for  AImage in range(NEReduced)];
	CCCtoSolveReducedSortedAve.append( sum(ww)/len(ww));

    print('This is CCCtoSolveReducedSortedAve');
    print(CCCtoSolveReducedSortedAve);
	
    for j in range(NE):
        vv=CCCtoSolveReducedSorted[j];
        pylab.plot(vv);
	pylab.hold(True);

    pylab.title('This is CCCtoSolveReducedSorted');
    pylab.draw();
    pylab.savefig('SWpng/CCCtoSolveReducedSorted.png')
    pylab.hold(False);

#
#%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%           8.
#%           Now find effective SNR's of each class. Find gammas (with product 1)  such that the
#%           overall variance is minimized.  Pic1*Picj* gamma1*gammaj
#%                                                approx
#%                                           Pic2*Picj* gamma2*gammaj
#%          Assume you can do that independently
#% 

    CCCtoSolveReducedGamma           = [ CCCtoSolveReduced[j] for j in range(len(CCCtoSolveReduced))]          ;
    CCCtoSolveReducedGammaSorted     = [CCCtoSolveReducedSorted[j] for j in range(len(CCCtoSolveReducedSorted))]    ;
    CCCtoSolveReducedGammaSortedAve  = [];

    for BImage in range(NEReduced):
	ww= [CCCtoSolveReducedGammaSorted[AImage][BImage] for  AImage in range(NEReduced)];
	CCCtoSolveReducedGammaSortedAve.append( sum(ww)/len(ww));

    gammaFitRunning                  = [1 for  j in range(NEReduced)]  ;

    NumPTarget =4;
    gammaFit =  [1 for  j in range(NEReduced)]  ;
    relax=.1;  # a low value means that it moves slowly to the new value
    
    print('CCCtoSolveReducedGammaSortedAve'+str(CCCtoSolveReducedGammaSortedAve))
#
#

#    for jjj in range(100):
#        lCC= length(CCCtoSolveReducedGammaSortedAve);
#        rangeGamma = ( (lCC- NumPTarget): (lCC-1));
#       Target = CCCtoSolveReducedGammaSortedAve(rangeGamma); # we are trying to pick gamma 
#                                                              # so that the shape of the curves look like this one
#    
#    for jj=1:NEReduced
#        jjo = jj;  jjx = 2*NEReduced+jj; jjy = 4*NEReduced+1-jj;  jjr = 2*NEReduced+1-jj;  
#       % we also need to add in the various copies
#        mapToTargeto = CCCtoSolveReducedGammaSorted(rangeGamma,jjo); normo =mapToTargeto'*mapToTargeto;
#        mapToTargetx = CCCtoSolveReducedGammaSorted(rangeGamma,jjx); normx =mapToTargetx'*mapToTargetx;
#        mapToTargety = CCCtoSolveReducedGammaSorted(rangeGamma,jjy); normy =mapToTargety'*mapToTargety;
#        mapToTargetr = CCCtoSolveReducedGammaSorted(rangeGamma,jjr); normr =mapToTargetr'*mapToTargetr;
#        norm = normo + normx + normy + normr;
#        mapToTarget =   mapToTargeto +  mapToTargetx +  mapToTargety +  mapToTargetr ;
#        gammaFit(jj)=  (1-relax)*gammaFit(jj) + relax*(  mapToTarget' * Target ) / norm;
#    end   ; clear mapToTarget mapToTargeto mapToTargetx mapToTargety mapToTargetr jj0 jj jjx jjy jjr norm normo normx normy normr
#    
#    logprodGammaFit = (sum( log(gammaFit)));
#    gammaFit = gammaFit/(exp( logprodGammaFit /NEReduced));
#    % The above expression for gammafit solves
#    %   Target = gammaFit  times   mapToTargeto, etc.
#    
#    %   Now we need to change all 32 elements of CCCtoSolveReducedGamma that result from changing
#    %   these gamma's
#    for jj=1:NEReduced
#        jjo = jj;                  jjMappedo = rangeCCCtoSolve(jj)  ;
#        jjy = 2*NEReduced+jj;      jjMappedy = NE/2 + jjMappedo     ;
#        jjx = 4*NEReduced+1-jj;    jjMappedx = NE   - jjMappedo + 1 ;
#        jjr = 2*NEReduced+1-jj;    jjMappedr = NE/2 - jjMappedo + 1 ;
#        for jk=1:NEReduced
#            jko = jk;                  jkMappedo = rangeCCCtoSolve(jk) ;
#            jky = 2*NEReduced+jk;      jkMappedy = NE/2 +     jkMappedo;
#            jkx = 4*NEReduced+1-jk;    jkMappedx = NE   + 1 - jkMappedo ;
#            jkr = 2*NEReduced+1-jk;    jkMappedr = NE/2 + 1 - jkMappedo ;
#
#            CCCtoSolveReducedGamma(jjo,jko) = CCCtoSolveReducedGamma(jjo,jko) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jko,jjo) = CCCtoSolveReducedGamma(jjo,jko) ;
#            CCCtoSolveReducedGamma(jjx,jko) = CCCtoSolveReducedGamma(jjx,jko) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jko,jjx) = CCCtoSolveReducedGamma(jjx,jko) ;
#            CCCtoSolveReducedGamma(jjy,jko) = CCCtoSolveReducedGamma(jjy,jko) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jko,jjy) = CCCtoSolveReducedGamma(jjy,jko) ;
#            CCCtoSolveReducedGamma(jjr,jko) = CCCtoSolveReducedGamma(jjr,jko) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jko,jjr) = CCCtoSolveReducedGamma(jjr,jko) ;
#            CCCtoSolveReducedGamma(jjo,jkx) = CCCtoSolveReducedGamma(jjo,jkx) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkx,jjo) = CCCtoSolveReducedGamma(jjo,jkx) ;
#            CCCtoSolveReducedGamma(jjx,jkx) = CCCtoSolveReducedGamma(jjx,jkx) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkx,jjx) = CCCtoSolveReducedGamma(jjx,jkx) ;
#            CCCtoSolveReducedGamma(jjy,jkx) = CCCtoSolveReducedGamma(jjy,jkx) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkx,jjy) = CCCtoSolveReducedGamma(jjy,jkx) ;
#            CCCtoSolveReducedGamma(jjr,jkx) = CCCtoSolveReducedGamma(jjr,jkx) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkx,jjr) = CCCtoSolveReducedGamma(jjr,jkx) ;
#            CCCtoSolveReducedGamma(jjo,jky) = CCCtoSolveReducedGamma(jjo,jky) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jky,jjo) = CCCtoSolveReducedGamma(jjo,jky) ;
#            CCCtoSolveReducedGamma(jjx,jky) = CCCtoSolveReducedGamma(jjx,jky) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jky,jjx) = CCCtoSolveReducedGamma(jjx,jky) ;
#            CCCtoSolveReducedGamma(jjy,jky) = CCCtoSolveReducedGamma(jjy,jky) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jky,jjy) = CCCtoSolveReducedGamma(jjy,jky) ;
#            CCCtoSolveReducedGamma(jjr,jky) = CCCtoSolveReducedGamma(jjr,jky) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jky,jjr) = CCCtoSolveReducedGamma(jjr,jky) ;
#            CCCtoSolveReducedGamma(jjo,jkr) = CCCtoSolveReducedGamma(jjo,jkr) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkr,jjo) = CCCtoSolveReducedGamma(jjo,jkr) ;
#            CCCtoSolveReducedGamma(jjx,jkr) = CCCtoSolveReducedGamma(jjx,jkr) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkr,jjx) = CCCtoSolveReducedGamma(jjx,jkr) ;
#            CCCtoSolveReducedGamma(jjy,jkr) = CCCtoSolveReducedGamma(jjy,jkr) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkr,jjy) = CCCtoSolveReducedGamma(jjy,jkr) ;
#            CCCtoSolveReducedGamma(jjr,jkr) = CCCtoSolveReducedGamma(jjr,jkr) * gammaFit(jj)*gammaFit(jk) ; CCCtoSolveReducedGamma(jkr,jjr) = CCCtoSolveReducedGamma(jjr,jkr) ;
#        end
#        CCCtoSolveReducedGamma(jjo,jjo) = 1;
#        CCCtoSolveReducedGamma(jjx,jjx) = 1;
#        CCCtoSolveReducedGamma(jjy,jjy) = 1;
#        CCCtoSolveReducedGamma(jjr,jjr) = 1;
#    end
#    
#    CCCtoSolveReducedGammaSorted =sort(CCCtoSolveReducedGamma(1:4*NEReduced,1:4*NEReduced)); % CHANGE ME
#%    figure;    plot(CCCtoSolveReducedGammaSorted,'.'); 
# 
#    CCCtoSolveReducedGammaSortedAve= mean(CCCtoSolveReducedGammaSorted,2);
#    gammaFitRunning = gammaFitRunning.* gammaFit;
#    gammaFitRecord(jjj,1:NEReduced)=gammaFit(1:NEReduced);
#end 
#
##
#titleString  = 'Section 8: This is CCCtoSolveReducedGamma';
#saveasString = 'Section8CCCtoSolveReducedGamma'; figureName   = 'Section8CCCtoSolveReducedGamma'; temp=CCCtoSolveReducedGamma;
#MakeImage(titleString,saveasString,figureName, temp)
#delete(gcf)
#
#figure;    plot(CCCtoSolveReducedGammaSorted,'.'); 
#title(strcat('Section 8: This is CCCtoSolveReducedGammaSorted'));
#saveas(gcf,'Section8CCCtoSolveReducedGammaSortedB.eps','psc');
#delete(gcf)
#
    LL = len(CCCtoSolveReducedGammaSortedAve);

    bbfullTrue  = [CCCtoSolveReducedGammaSortedAve[j] for j in range(LL)     ] ;
    for jj in range(2*oddCp1*NumImgsRed-1):
        bbfullTrue.append(bbfullTrue[2*oddCp1*NumImgsRed-2-jj]);
	
    print('This is bbfullTrue: '+str(bbfullTrue));

    lbbfullTrue=len(bbfullTrue);
    thetabb=[-180.0/nCSym  + (float(j)/(lbbfullTrue-1))*(360.0/nCSym) in range(lbbfullTrue)]
    
    pylab.hold(False);
    pylab.plot(bbfullTrue);    pylab.title('bbfullTrue');
    pylab.draw();     pylab.savefig('SWpng/bbfullTrue.png')
    if (eogMode):     os.system('eog SWpng/bbfullTrue.png &');


#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%           9.   Now find histograms
#%          
#%
#
#cutoff = .41;  % .55 TPPII
#clear NumB
#for j=1:NEReduced bbfullTrue
#    NumB(j) =length(find(CCCtoSolveReducedGamma(j,:)>cutoff))-1; 
#end    
#
#figure;
#hist(NumB,NEReduced/2);
#Title(strcat('This is  a set of histograms for how many class averages correlate with correlation', ...
#             'coefficient >', num2str(cutoff),' with their neighbors. The set of class averages is the set of 100 HRS ' , ...
#              'particles from KLM8200. The histogram depends sensitively on the cutoff. This gives us ', ...
#             ' a way to pick the cutoff => make the histogram as flat as possible. ' , ...
#             '  The histogram together with cutoff also give us a way to assess the SNR.  2/6/04'));
#
#saveas(gcf,'Section9hist.eps','psc');
#delete(gcf);
#         
#
#%This is  a set of histograms for how many class averages correlate with correlation
#%coefficient > .83 with their neighbors. The set of class averages is the set of 200 HRS
#%particles from KLM8200. The histogram depends sensitively on the cutoff. This gives us
#%a way to pick the cutoff => make the histogram as flat as possible.
#% The histogram together with cutoff also give us a way to assess the SNR
#         
#
#
#
#%%%%%%%%%%
#%
#%      10. Now try to infer the distance between each of the pairs based on the
#%             function bbfull
#%      First, find a fit to the peak in terms of thetabb
#bbfullTrue
#


#  This is the second cutoff 
    cutoff= bbfullTrue[LL-7]-.00001;  # .86 for SNR3 for HRS; .55; .84 for TPPII classes; .93 for HRS
    bbIndices = find(bbfullTrue,'gt',cutoff);   
    lbbInd=len(bbIndices); 
    midbb = (lbbInd+1)/2;
    bbfullTrue2Use = [ bbfullTrue[j] for  j  in bbIndices];

    print('This is the cutoff:'+str(cutoff));
    print('This is bbfullTrue2Use: ' + str(bbfullTrue2Use)+'\n');
    print('This is bbIndices: ' + str(bbIndices)+'\n');
    
   

#There need to be NEReduced particles spaced evenly on [0.. 360/nCSym)

    NEReducedTot = 2*2*oddCp1*NumImgsRed-2;
    bbxstep= (360/nCSym/NEReducedTot);
    bZ = NEReducedTot/2;
    bbxmin = -bbxstep*bZ;
    bbxmax= bbxmin + (2*bZ)*bbxstep; # So bbxmax-bbxmin= NEReducedTot*bbxstep=360/CSym
    thetabb = pylab.arange(bbxmin,bbxmax+bbxstep/2.0,bbxstep);  # that is thetabb is [-180/CSym .. 180/CSym]
    thetabb2Use = [thetabb[j] for j in bbIndices];

    print();
    print('This is thetabb2Use: ' + str(thetabb2Use));

    pylab.plot(thetabb2Use,bbfullTrue2Use);
    pylab.title('Section 10: bbfullTrue2Use as a function of thetabb2Use');
    pylab.draw();
    pylab.savefig('SWpng/bbfullTrue2Use.png')
    if (eogMode):     os.system('eog SWpng/bbfullTrue2Use.png &');



    print('NE, NEReduced'+str(NE)+', '+str(NEReduced));

    print('CCCtoSolveReducedGamma'+str(CCCtoSolveReducedGamma));

    mat2Use = [[0 for j1 in range(2*oddCp1*NumImgsRed)] for j2 in range(2*oddCp1*NumImgsRed)];
    thetaIJGuess = [[0 for j1 in range(2*oddCp1*NumImgsRed)] for j2 in range(2*oddCp1*NumImgsRed)];
    foundSomething = [[0 for j1 in range(2*oddCp1*NumImgsRed)] for j2 in range(2*oddCp1*NumImgsRed)];
#
    secx= [bbfullTrue2Use[j] for j in pylab.arange(midbb,lbbInd) ];
    secy= [thetabb2Use[j]  for j in pylab.arange(midbb,lbbInd) ] ;

    print('secx '+ str(secx))
    print('secy' + str(secy))
    
    for jr in range(2*oddCp1*NumImgsRed):
        jro = jr;                
        for jc in range(2*oddCp1*NumImgsRed):
            jco = jc;
            value = CCCtoSolveReducedGamma[jro][jco];
            if (( value >= cutoff) & (jro != jco)):
                temp  =interp1(secy,secx,value);
                thetaIJGuess[jro][jco] = temp;
                thetaIJGuess[jco][jro] = temp;
                foundSomething[jro][jc] = 1;
                mat2Use[jro][jc] = -1;
                mat2Use[jro][jro] = mat2Use[jro][jro]+1;
		#print('value,temp,jro,jco '+str(value)+', '+str(temp)+', '+str(jro)+' ,'+ str(jco))

    print('mat2Use = '+str(mat2Use));
    print('thetaIJGuess ='+str(thetaIJGuess));

    print('oddC ='+str(oddC));
    print('DSym ='+str(DSym));


#
#%
#%      11       Now fold the conditions
#

    thetaIJCondM = [[[] for j1 in range(NumImgsRed)] for j2 in range(NumImgsRed)];
    thetaIJCondP = [[[] for j1 in range(NumImgsRed)] for j2 in range(NumImgsRed)];


    for jj in range(NumImgsRed):
        for kk in range(jj+1,NumImgsRed):
            ttt= thetaIJGuess[jj][kk];    
            if ((ttt)>0):
	        vvM = thetaIJCondM[jj][kk]
		vvM.append(ttt);  # +0*random
                thetaIJCondM[jj][kk]= vvM ;
                thetaIJCondM[kk][jj]= vvM ;

    for jj in range(NumImgsRed):
        if (DSym):       jjx = 2*NumImgsRed *oddCp1  -1-jj ;  
        if (oddC):       jjy = (DSym+1)*NumImgsRed   +jj    ;  
        if (oddC&DSym):  jjr = 2*NumImgsRed        -1-jj ;  
#    
        for kk in range(jj,NumImgsRed):
#        
            if (oddC):
                tempy =thetaIJGuess[jjy][kk]; 
                if (tempy>0):
                    vvM=thetaIJCondM[jj][kk];
                    vvM.append(180/CSym + tempy)
                    thetaIJCondM[jj][kk]=vvM;
                    thetaIJCondM[kk][jj] =  thetaIJCondM[jj][kk];
                [jj,kk,jjx, tempx];
#          
#        
            if (DSym):
                tempx =thetaIJGuess[jjx][kk]; 
                if (tempx>0):
                    vvP=thetaIJCondP[jj][kk];
		    vvP.append(tempx);
                    thetaIJCondP[jj][kk] = vvP ;
                    thetaIJCondP[kk][jj] =  thetaIJCondP[jj][kk];
                [jj,kk,jjx, tempx];
#                %              pause
#        
            if (oddC & DSym):
                tempr =thetaIJGuess(jjr,kk); 
                if (tempr>0):
                    vvP = thetaIJCondP[jj][kk];
		    vvP.append(180/CSym+tempr );
                    thetaIJCondP[jj][kk]= vvP ;
                    thetaIJCondP[kk][jj]=  thetaIJCondP[jj][kk];
                [jj,kk,jjr,180/CSym+tempr];
#   
        #    end kk loop
    #    end jj loop
#
#
    print('thetaIJCondM ='+str(thetaIJCondM));
    print('thetaIJCondP ='+str(thetaIJCondP));

#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% %%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%   12  We now collect all the conditions we wish to satisfy into two
#%          linear array.
#% 
#
    JXM=[]; JXP=[]; JYM=[]; JYP=[];
    thetaIJP = [];    thetaIJM = [];
#
    for jj in range(NumImgsRed):
        for kk in range(jj,NumImgsRed):
            vv = thetaIJCondP[jj][kk];
            for qq in range(len(vv)):
                thetaIJP.append(vv[qq]);
                JXP.append(jj);
                JYP.append(kk);
            vv = thetaIJCondM[jj][kk];
            for qq in range(len(vv)):
                thetaIJM.append(vv[qq]);
                JXM.append(jj);
                JYM.append(kk);

    PPP=[ [JXP[jj],JYP[jj],thetaIJP[jj]] for jj in range(len(JXP))];
    MMM=[ [JXM[jj],JYM[jj],thetaIJM[jj]] for jj in range(len(JXM))];

#

    thetaIJMRad = [thetaIJM[mm]  * pi/180. for mm in range(len(thetaIJM))];
    thetaIJPRad = [thetaIJP[mm]  * pi/180. for mm in range(len(thetaIJP))];
    
    print('MMM ='+str(MMM));
    print('PPP ='+str(PPP));


#  ---------------------------------------------------
#
#   13  We perform the annealing procedure; Find Best Value
#

    thetaSuperVec,ResidualSuperVec,ResidualvsTime, ModeVec = Anneal(NumImgsRed,thetaIJMRad,JXM,JYM,thetaIJPRad,JXP,JYP,SymString,NumLoops,MovesPerLoop) 

#     Can serve as a test
    rr = 45.0*(2-oddC)/nCSym;   

    slots= [ (j *(2.0*rr/NumImgsRed)  + rr/NumImgsRed  ) for j in range(NEReduced)];
    slotsRad=[ slots[j]*pi/180.0 for j in range(len(slots)/2)];
#
    Residual = MinFuncPM(slotsRad,thetaIJMRad,JXM,JYM, thetaIJPRad,JXP,JYP,SymString) ; print('Residaul of slots is'); print(Residual); 
#
    val= min(ResidualSuperVec);
    index = find(ResidualSuperVec,'eq',val);
    print('ResidualSuperVec');    print(ResidualSuperVec);
    print('val, index'); print(val,index[0]);
    BestAngs= thetaSuperVec[index[0]];

    print('BestAngs'); print(BestAngs);

#
#saveasString = 'Section15FinalTheta';
#figure('name',saveasString);
#
#plot(vv,'.')      ;
#% clear index val
#saveasString = 'Section15FinalTheta';
#whitebg('w') %create a figure with a white color scheme
#xlabel('position ' , 'FontSize' ,14, 'FontWeight','bold')
#ylabel('angular assignment (degrees)' , 'FontSize' ,14, 'FontWeight','bold')
#
#set(gca,'FontSize',11,'FontWeight','bold', 'FontSize' ,14, 'FontWeight','bold')  
#
#title('Best Theta  ', 'FontSize' ,14, 'FontWeight','bold')
#saveas(gcf,strcat(saveasString,'.fig'),'fig');saveas(gcf,strcat(saveasString,'.eps'),'eps');saveas(gcf,strcat(saveasString,'.png'),'png')
#delete(gcf)
#
# CSym=2; AngMax= 360/CSym; DSym=1; %That means Yes 
#figure; plot( sort([AngMax-vv,vv]),'.') % This checks that they are distinct
#
#
#%  ---
#
#figure;
#  vwb =.5:29.5
#  vwb= sort(AngMax/4 - abs(mod(vv,AngMax/2)- AngMax/4))
#  plot(vwb,'.')      ;
#   
# whitebg('w') %create a figure with a white color scheme
# xlabel('position ' , 'FontSize' ,14, 'FontWeight','bold')
# ylabel('angular assignment (degrees)' , 'FontSize' ,14, 'FontWeight','bold')
# 
# set(gca,'FontSize',11,'FontWeight','bold', 'FontSize' ,14, 'FontWeight','bold')  
#
# title('Best Theta  ', 'FontSize' ,14, 'FontWeight','bold')
#saveasString = 'Section15PerfectTheta';
#saveas(gcf,strcat(saveasString,'.fig'),'fig');saveas(gcf,strcat(saveasString,'.eps'),'eps');saveas(gcf,strcat(saveasString,'.png'),'png')
#
#
#
#
#figure;
#plot(ResidualvsTime)
#xlabel('time ' , 'FontSize' ,14, 'FontWeight','bold')
# ylabel('Residual' , 'FontSize' ,14, 'FontWeight','bold')
# 
# set(gca,'FontSize',11,'FontWeight','bold', 'FontSize' ,14, 'FontWeight','bold')  
#
# title('Residual vs Time  ', 'FontSize' ,14, 'FontWeight','bold')
#saveasString = 'Section15ResidualvsTime';
#saveas(gcf,strcat(saveasString,'.fig'),'fig');saveas(gcf,strcat(saveasString,'.eps'),'eps');saveas(gcf,strcat(saveasString,'.png'),'png')
#hold;
#plot(80*exp(- 5*(1:10000)/10000))
#
#
# -------------------------------------------------------------------------------
#
# 14.  Create outputs that are written to file; a) Reconstruction; b) Reprojection: c) EvenlySpaced
#
#

    print('slots')
    print(slots)
    DataSetBase=DataSetName[0:(len(DataSetName)-4)];
    DataSetAnglesName='SWfiles/SW_'+DataSetBase+'Angles.txt';
    DataSetHDFName='SWfiles/SW_'+DataSetBase+'.hdf';

    print(DataSetHDFName);

    VolString          = 'SWfiles/SWReconVol.mrc'
    ReProjString       = 'SWfiles/SWReproj.spi'
    EvenlySpacedString = 'SWfiles/SWEvenlySpaced.spi'


    try: 
    	fsock= open(DataSetAnglesName);
	os.remove(DataSetAnglesName); 
    except IOError: print('The File '+DataSetAnglesName + ' Does Not Exist');

    try: 
    	fsock= open(DataSetHDFName);
	os.remove(DataSetHDFName); 
    except IOError: print('The File '+DataSetHDFName + ' Does Not Exist');

    try: 
    	fsock= open(VolString);
	os.remove(VolString); 
    except IOError: print('The File '+VolString + ' Does Not Exist');
    try: 
    	fsock= open(ReProjString);
	os.remove(ReProjString);
    except IOError: print('The File '+ ReProjString+ ' Does Not Exist');
    try: 
    	fsock= open(EvenlySpacedString);
        os.remove(EvenlySpacedString);
    except IOError: print('The File '+ EvenlySpacedString+ ' Does Not Exist');


    outAngles = open(DataSetAnglesName, "w"); # A) write angles out to txt file
    for AImage in range(NumImgsRed):  # 
        outAngles.write("%d\t%f\n" % (AImage+1,BestAngs[AImage])) 

    outAngles.close(); 
    
    print('NE,NEReduced,NumImgs,NumImgsRed');     print(NE,NEReduced,NumImgs,NumImgsRed);
    print('rangeCCCtoSolve');    print(rangeCCCtoSolve);

    A=EMData();
    for AImage in range(NumImgsRed):  #    B) write selected images to hdf file
        A.read_image(DataSetName,rangeCCCtoSolve[AImage] );
        print(BestAngs[AImage]);
        A.set_attr_dict({'phi':BestAngs[AImage], 'theta':90.0, 'psi':-90.0, 'sx':0, 'sy':0});
        A.write_image(DataSetHDFName, AImage);

    
    Vol1=recons3d_4nn(DataSetHDFName, range(NumImgsRed),'d2');  # C) Do a reconstruction to mrc file
    Vol1.write_image(VolString)
    Vol1ft,kb=prep_vol(Vol1);
    
    A = EMData();  # D) write out original/reprojection pairs to spi file

    for AImage in range(NumImgsRed):
	A.read_image(DataSetHDFName,AImage);
        A.write_image(ReProjString, 2*AImage);
        phi = BestAngs[AImage];
        proji  = prgs(Vol1ft,kb,[phi,90,-90,0,0]);
        proji.write_image(ReProjString,2*AImage+1);
    
    
    for iphi in range(NEReduced):  # E) write out evenly spaced projections
        phi = slots[iphi];
        proji  = prgs(Vol1ft,kb,[phi,90,0,0,0]);
        proji.write_image(EvenlySpacedString,iphi);
    

    return;



def write(List,N,M,fN):
    lfN=len(fN)
    suffix=fN[lfN-3:lfN];
    
    if (suffix=='txt'):
        try: 
            outAngles= open(fN); os.remove(fN); 
        except IOError: print('The File '+fN + ' Does Not Exist');
        outfN = open(fN, "w"); #  write angles out to txt file
        Count=0;
        for ix in range(len(List)):  # 
	    if ((ix%N)==0):
	        Count=Count+1;
                if (ix>0):
                    outfN.write("\n");
                outfN.write("%d\t" % Count);
            outfN.write("%3.2f " % (List[ix]));
        outfN.close();
	return;    
    
    e=EMData();
    e.set_size(N,M);
    Count=0;
    for ix in range(N):
        for iy in range(M):
	    e.set_value_at(ix,iy,List[Count]);
	    Count=Count+1;
	    
    e.write_image(fN);
    return

'''
