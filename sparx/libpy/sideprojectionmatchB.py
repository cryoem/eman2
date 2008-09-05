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
from EMAN2_cppwrap import *
from global_def import *
#from filter import *
#from fundamentals import *
#from math import *
#from random import *
#from sparx import *
#from sys import maxint
#from utilities import  *



'''

#    return thetaSuperVec,ResidualSuperVec,ResidualvsTime, ModeVec

def Anneal(NImgsRed,thetaIJMRad,JXM,JYM,thetaIJPRad,JXP,JYP,SymString,LLLMax=20,LLMax=60000):
    """
		August 14th 2006
		Read  in a 
			NImgsRed   - number of images
			thetaIJMRad - these are the conditions on theta minus
			JXM         - these index the minus conditions on X
			JYM         - these index the minus conditions on Y
			thetaIJPRad - these are the conditions on theta plus
			JXP         - these index the plus conditions on X
			JYP         - these index the plus conditions on Y
			SymString   - like 'd2' or 'c3'
			LLLMax      - the number of rounds of independent refinements
			LLMax       - the number of "moves" per refinement round

		Returns
			thetaSuperVec    - A listing of the best assignments for each round
			ResidualSuperVec - A listing of the best residuals
			ResidualvsTime   - A typical run of the residual
			
    """

    nCSym  =   int(SymString[1]);
    oddn   =   nCSym%2;
    oddnp1 =   oddn+1;

    rr = 45.0*(2.0-oddn)/nCSym;

    slots= [ (j *(2.0*rr/NImgsRed)  + rr/NImgsRed  ) for j in range(2*oddnp1*NImgsRed)];

    AngMax=360/nCSym;
     
    theta = [ slots[j]  for j in range(NImgsRed)]; 


# ----------------------------------
    
    if (oddn):  #  Begin axial symmetry, n, odd
#
#      A.  Begin odd axial symmetry
#       
        ResidualSuperVec =[];
        thetaSuperVec    =[];
	
        for LLL in range(LLLMax):
#
#         A0.  Create an initial assignment   
#            
            Assignment = [ 0  for j in range(4*NImgsRed)]; 
            for jjo in range(NImgsRed):
                kk  =   randrange(4*(NImgsRed-jjo));  # this will be the slot for theta(jjo)
	        kkoVec = find(Assignment,'eq',0); 
                kko = kkoVec[kk]; # this is the index for slots

                theta[jjo] = slots[kko];
                kkr = (2*NImgsRed - kko -1)%(4*NImgsRed) ;
                kky = (2*NImgsRed + kko   )%(4*NImgsRed) ;  
                kkx = (4*NImgsRed - kko -1)%(4*NImgsRed) ; 
                Assignment[kko]=1;   
		Assignment[kkr]=1;  
		Assignment[kky]=1;   
                Assignment[kkx]=1;   
#               print([kko,kkx,kky,kkr])
            
            
 #       thetaGuess=slots(1:NE);
#            
            
            thetaRad = [ theta[j]*pi/180  for j in range(NImgsRed)]; 
            Residual = MinFuncPM(thetaRad,thetaIJMRad,JXM,JYM, thetaIJPRad,JXP,JYP,SymString)  ; # this is the initial value
#
#         A1.  Annealing
#           
            NEpick = int(NImgsRed*(NImgsRed-1)/4); # once in every NEpick the program does something else
            Count=0;
            Countjj=0;
            ResidualVec=[1000000];#       The different trials just ended 
            thetaHist = [[ theta[j]  for j in range(NImgsRed)]]; 
            Mode=0;
            ResidualvsTime  = [];
            ModeVec = [];
            
            for LL in range(LLMax):
                Temp=20*exp(- 10*LL/float(LLMax));
                
                thetaTry=[ theta[j]  for j in range(NImgsRed)]; 
                
                if ((LL%NEpick)>0): #  &&  mod(LL,NImgsRed+2)~= -1    A) Switch Random Pairs
                    kki  = randrange(NImgsRed);
                    kkf  = randrange(NImgsRed);
                    kkfS =randrange(4);
                    Mode=0;
                    
                    thetaTry[kkf] = theta[kki];

                    if (kkfS==0):  # same
                        thetaTry[kki] = theta[kkf];
                    if (kkfS==1):  # r[j]
                        thetaTry[kki] =  ( AngMax/2 - theta[kkf])%AngMax;
                    if (kkfS==2):   # y
                        thetaTry[kki] =  ( AngMax/2 + theta[kkf])%AngMax;
                    if (kkfS==3):  # x
                        thetaTry[kki] =  ( AngMax   - theta[kkf])%AngMax;

                if (LL%NEpick==0):     #     B) Permutation
                    Mode=randrange(1,3); # that is 1 or 2
                    thetaTry = [ theta[((j+2*Mode-3)%NImgsRed)]  for j in range(NImgsRed)];
                    
                vv1= [ AngMax/4.0 - abs( (thetaTry[j])%(AngMax/2.0)- AngMax/4.0 )  for j in range(NImgsRed)  ] ;
                vv= sorted(vv1);
                if ( sum( [ vv[j] == vv[j+1] for j in range(NImgsRed-1)   ] )):
                    print('Error');  # means that two have the same assignment
            
                thetaTryRad = [ thetaTry[j]*pi/180  for j in range(NImgsRed)]; 
                ResidualNew = MinFuncPM(thetaTryRad,thetaIJMRad,JXM,JYM, thetaIJPRad,JXP,JYP,SymString)  ;
                
                DeltaResidual=ResidualNew-Residual ;
                if ((DeltaResidual<0) |  (random()< exp(-DeltaResidual/Temp) )):  # move accepted
                    theta=[ thetaTry[j]  for j in range(NImgsRed)]; 
                    Residual=ResidualNew;
                    if (Residual <  ResidualVec[Count]):
                        Count=Count+1;
                        thetaHist.append(theta);
                        ResidualVec.append(Residual);
#	      print(Residual); print(theta);
                        ModeVec.append(Mode);
                #  end updates after a good move
                ResidualvsTime.append(Residual);
            # ends LL run
            iii=len(  ResidualVec )    -1  ;
            thetaSuperVec.append(thetaHist[iii]);  # append the last best value of that trial
            ResidualSuperVec.append(ResidualVec[iii]) ; 
#      end of LLL superRun
#    The case of C odd just ended 
    
    
# -----------------------------------------------------------------    
    
    
    if ((nCSym+1)%2):  #  Begin axial symmetry, n, even
#
#      B.  Begin even axial symmetry
#       
        ResidualSuperVec  =[];
        thetaSuperVec     =[];
	
        for LLL in range(LLLMax):
#
#         B0.  Create an initial assignment   
#            
            Assignment = [ 0  for j in range(2*NImgsRed)]; 
            for jjo in range(NImgsRed):
                kk  =   randrange(2*(NImgsRed-jjo));  
                kkoVec = find(Assignment,'eq',0);
                kko= kkoVec[kk];    # this is the index for slots
                theta[jjo]= slots[kko];
                kkx = (2*NImgsRed - kko -1)%(2*NImgsRed) ; 
                Assignment[kko]=1;   
                Assignment[kkx]=1;   
 #       print([kko,kkx]); # a debug statement
  
            thetaRad = [ theta[j]*pi/180  for j in range(NImgsRed)]; 
            Residual = MinFuncPM(thetaRad,thetaIJMRad,JXM,JYM, thetaIJPRad,JXP,JYP,SymString)  ; # this is the initial value
	    print(Residual);
 #
#         B1.  Annealing
#           
            NEpick = int(NImgsRed*(NImgsRed-1)/2); # once in every NEpick the program does something else
            Count=0;
            Countjj=0;
            ResidualVec=[1000000];
            thetaHist= [[ theta[j]  for j in range(NImgsRed)]];
            Mode=0;
            ResidualvsTime  = [];
            ModeVec         = [];
              
            for LL in range(LLMax):
                Temp=20*exp(- 5.0*LL/LLMax);
                
                thetaTry=[ theta[j]  for j in range(NImgsRed)]; 
                
                if ((LL%NEpick)>0): #  &&  mod(LL,NImgsRed+2)~= -1    A) Switch Random Pairs
                    kki  = randrange(NImgsRed);
                    kkf  = randrange(NImgsRed);
                    kkfS = randrange(2);
                    
                    Mode=0;
 
                    thetaTry[kkf] = theta[kki];
  
                    if (kkfS==0):  # same
                        thetaTry[kki] = theta[kkf];
                    if (kkfS==1):  # x
                        thetaTry[kki] =  ( AngMax - theta[kkf])%AngMax;

                
                if ( ((LL%NEpick)==0)  ): #  ((LL%(NImgsRed+2))==0)  ):
                    if ((LL%NEpick)==0) :  #  Permutation
                        Mode=randrange(1,3);        # that is 1 or 2
                        thetaTry = [ theta[((j+2*Mode-3)%NImgsRed)]  for j in range(NImgsRed)];
#
                    else:  # switch two in the asymmetric unit
                        Mode=3;
                        Countjjp1 = (Countjj +1)%(2*NImgsRed) ;
#                 Switch the jj slot with the jj+1 slot
                        
                        sLjj = slots[Countjj]; # This is the value  at the jjth slot
                        dVec  =  [  (  AngMax + theta[j]-sLjj)%AngMax  for j in range(NImgsRed) ];
                        dVecx =  [  (4*AngMax - theta[j]-sLjj)%AngMax  for j in range(NImgsRed) ];
			
			dVecm  = min(dVec );
			dVecxm = min(dVecx);
			
			if (dVecm <.00001):  
				nV  = find(dVec,'eq', dVecm);
				tjj = nV[0];
				thetaTry[tjj] = slots[Countjjp1]; 
				
			if (dVecxm<.00001):  
				nV = find(dVecx,'eq', dVecxm);
				tjj = nV[0];
				thetaTry[tjj] = slots[2*NImgsRed-Countjjp1-1];
				
                        sLjjp1 = slots[Countjjp1];
                        dVec  =  [  (  AngMax + theta[j]-sLjjp1)%AngMax  for j in range(NImgsRed) ];
                        dVecx =  [  (4*AngMax - theta[j]-sLjjp1)%AngMax  for j in range(NImgsRed) ];
			
			dVecm  = min(dVec );
			dVecxm = min(dVecx);
			
			if (dVecm <.00001):  
				nV    = find(dVec,'eq', dVecm);
				tjjp1 = nV[0];
				thetaTry[tjjp1] = slots[Countjj]; 
				
			if (dVecxm<.00001):  
				nV    = find(dVecx,'eq', dVecxm);
				tjjp1 = nV[0];
				thetaTry[tjjp1] = slots[2*NImgsRed-Countjj-1];
				
                           
                        Countjj = (Countjj+1)%(2*NImgsRed); # This is the slot
			print(Countjj,Mode)
                    # End  Mode 3
                #  End Modes 1 2 or 3
                
                vv1 = [ min([ AngMax-thetaTry[j],thetaTry[j] ]) for j in range(NImgsRed)   ];
		vv= sorted(vv1)
                if ( sum( [ (vv[j] == vv[j+1]) for j in range(NImgsRed-1)   ] ) ):
                    print('Error',vv);  # means that two have the same assignment
            
                thetaTryRad = [ thetaTry[j]*pi/180  for j in range(NImgsRed)]; 
                ResidualNew  = MinFuncPM(thetaTryRad,thetaIJMRad,JXM,JYM, thetaIJPRad,JXP,JYP,SymString)  ;
                
                DeltaResidual=ResidualNew-Residual ;
                if ((DeltaResidual<0) |  (random()< exp(-DeltaResidual/Temp) )):  # move accepted
                    theta=[ thetaTry[j]  for j in range(NImgsRed)]; 
                    Residual=ResidualNew;
                    if (Residual <  ResidualVec[Count]): # if it is the best config so far...
                        Count=Count+1;
                        thetaHist.append(theta);
                        ResidualVec.append(Residual);
#			print(ResidualVec);
#			print(thetaHist);
                        ModeVec.append(Mode);
                #  end updates after a good move
                ResidualvsTime.append(Residual);
            # ends LL run
            iii=len(ResidualVec)-1;
            thetaSuperVec.append(thetaHist[iii]);
            ResidualSuperVec.append(ResidualVec[iii]) ;
 #    end of LLL superRun
 #    End of C even
    
    return thetaSuperVec,ResidualSuperVec,ResidualvsTime, ModeVec






#  The minimization functional is  abs(cos(theta_i - theta_j)  - cos(theta_ij)) 
#
#  Returns Residual 
#

def MinFuncPM(thetaRad,thetaIJMRad,JXM,JYM,thetaIJPRad,JXP,JYP,SymString):  # all theta's are in Radians
    """
		August 14th 2006
		Read  in a 
			thetaRad    - the current angular assignment
			thetaIJMRad - these are the conditions on theta minus
			JXM         - these index the minus conditions on X
			JYM         - these index the minus conditions on Y
			thetaIJPRad - these are the conditions on theta plus
			JXP         - these index the plus conditions on X
			JYP         - these index the plus conditions on Y
			nCSym       - the axial symmetry, n

		Calculates      - ResM - Residual due to the negative conditions
		                - ResP - Residual due to the positive conditions
				
		Returns
			Residual  - The total value of the residual
			
    """

    nCSym=int(SymString[1]);
#    NE = len(thetaRad); # don't seem to need this actually
    NIJM = len(thetaIJMRad);  # Total Conditions from thetaIJM
    NIJP = len(thetaIJPRad);  # Total Conditions from thetaIJP
    ResM = 0 ;
    ResP = 0 ;
 
    for jj in range(NIJM):
         ResM = ResM + abs(cos(nCSym*(thetaRad[JXM[jj]] - thetaRad[JYM[jj]]))  - cos(nCSym*(thetaIJMRad[jj])) );
 

    for jj in range(NIJP):
         ResP = ResP + abs(cos(nCSym*(thetaRad[JXP[jj]] + thetaRad[JYP[jj]]))  - cos(nCSym*(thetaIJPRad[jj]))  );
 
    Residual=ResM + ResP;
 
    return Residual
 
 
#   -----------------------------------------------------------------------------------

def colorwheel(m,n):
	"""
	August 25th 2006
		Create a 3 vector for even sampling of the RGB colormap: 
		     n is the total number of colors to be assigned
		     m is which of the n colors is to be assigned
	"""
	nUp = int(n**(.3333333333))+1;
	m1=m%nUp;
	TnUp=2.0*nUp;
	m2= ((m-m1)/nUp)%nUp;
	m3=((m-m1-m2*nUp)/(nUp*nUp))%nUp;
	return [(2*m1+1)/TnUp,(2*m2+1)/TnUp,(2*m3+1)/TnUp]

#   -----------------------------------------------------------------------------------
def MakePlot(titleString,saveasString,figureName,tempX,tempY):

	pylab.figure(1); 
	pylab.plot(tempX,tempY,'.');  
#	h=pylab.gca;
#	pylab.set(h,'FontName','Arial','FontSize',12,'FontWeight','bold')# 'name',figureName,
#	pylab.title(titleString,'FontSize',14,'FontName','Arial','FontWeight','Bold')
#	pylab.saveas(gcf,strcat(saveasString,'.fig'),'fig');
#	pylab.saveas(gcf,strcat(saveasString,'.eps'),'eps');
#	pylab.saveas(gcf,strcat(saveasString,'.png'),'png')


#   -----------------------------------------------------------------------------------
def interp1(secy,secx,valuex):

	found=0;
	"""
	August 25th 2006
		interpolate to find valuey = f(valuex), given that
		we know        secy=f(secx),
		 where secx,secy are given lists
		 secx should be monotone( increasing or decreasing )
	"""
	NPts=len(secx);
	IncOrDec= (secx[1]>secx[0]) - (secx[0]>secx[1]);
	tempx =[float(IncOrDec*secx[j]) for j in range(NPts) ];
	tempvaluex=float(IncOrDec*valuex);
	
	# Now the sequence in x, tempx, is increasing

	found=0;
	if (tempvaluex<= min(tempx)): valuey=secy[0]; found=1 ;
	if (tempvaluex>= max(tempx)): valuey=secy[NPts-1]; found=1;
	
	j=0;
	while found==0:
		lowx=tempx[j]; highx=tempx[j+1];
		lowy= secy[j]; highy= secy[j+1];
		if ( (tempvaluex>lowx) & (tempvaluex < highx) ):
			valuey= ((tempvaluex-lowx)*highy  + (highx-tempvaluex)*lowy)/ (float(highx-lowx))
			found=1;
		j=j+1;
			
	return valuey;

# ------------------------------------------------------------
'''
