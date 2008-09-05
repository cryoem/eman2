
C++************************************************************************
C
C  FMRS_1.F                           ADDED FFTW AUG 2000 BIMAL RATH
C
C **************************************************************************
C *  SPIDER - MODULAR IMAGE PROCESSING SYSTEM.  AUTHOR: J.FRANK            *
C *  COPYRIGHT (C)1981,1987, WADSWORTH CENTER FOR LABORATORIES AND         *
C *  RESEARCH, NEW YORK STATE DEPARTMENT OF HEALTH, ALBANY, NY 12201.      *
C *  THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO THE CENTER FOR       *
C *  LABORATORIES AND RESEARCH AND ARE NOT TO BE DISCLOSED TO OTHERS OR    *
C *  USED FOR PURPOSES OTHER THAN INTENDED WITHOUT WRITTEN APPROVAL OF     *
C *  THE CENTER FOR LABORATORIES AND RESEARCH                              *
C
C  D REAL MIXED RADIX FFT.
C  INPUT:  X(N) - REAL ARRAY
C  OUTPUT: N EVEN  X(N+2)
C  ORDER OF ELEMENTS:
C  R(0),0.0, R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1), R(N/2),0.0
C
C         N ODD  X(N+1)
C  R(0),0.0, R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1), R(N/2),I(N/2)
C 
C  HERE WE FOLLOW THE CONVENTION THAT INTEGER DIVISION 
C  IS ROUNDED DOWN, E.G. 5/2 =2)
C
C  INV: +1 FORWARD FFT
C       -1 INVERSE FFT
C  ON OUTPUT INV=0 INDICATES ERROR
C
C  IMAGE_PROCESSING_ROUTINE
C
C--************************************************************************

        SUBROUTINE  FMRS_1(X,NSAM,INV)

        INCLUDE 'CMBLOCK.INC'
        PARAMETER (LBUF=5000)
        DIMENSION  X(*),WORK(LBUF)
	DIMENSION  Y(NSAM)


#ifdef SP_LIBFFT
C       USING SGI_COMPLIB FOR FFT
	N = NSAM
        IF (N+15 .GT. LBUF)  THEN
           INV = 0
C          INSUFFICIENT BUFFER, INCREASE LBUF AND COMPILE SPIDER
           CALL  ERRT(6,'FMRS_1',NE)
           RETURN
        ENDIF
	CALL SCFFT1DUI(N,WORK)
	LDA =1
	IF (INV.GT.0)  THEN
	   CALL SCFFT1DU(INV,N,X,LDA,WORK)
	ELSE
	   CALL CSFFT1DU(INV,N,X,LDA,WORK)
	   QT=1.0/FLOAT(N)
	   CALL  SSCAL1D(N,QT,X,LDA)
	ENDIF

#else
#if defined(SP_LIBFFTW) || defined(SP_LIBFFTWMP)
C       USING FFTW LIBRARY CALLS FOR FFT

#include "FFTW.INC"

        INTEGER, SAVE :: NSAMO=0
        INTEGER, SAVE :: NSAMOR=0

C       PLAN AND PLANR ARE ACTUALLY POINTERS TO A STRUCTURE 
#if defined (__osf__) || defined (ia64) || defined (__x86_64__)
        INTEGER*8, SAVE :: PLAN=0, PLANR=0
#else
        INTEGER, SAVE :: PLAN=0, PLANR=0
#endif
        LOGICAL, SAVE :: INIT=.TRUE.

#ifdef SP_LIBFFTWMP
        IF (INIT) THEN
C          MUST INITIALIZE THREADS ONCE
           CALL FFTW_F77_THREADS_INIT(IRTFLG);
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'MULTIPLE THREADS FAILED --FFTW2',IER)
              RETURN
           ENDIF
           INIT = .FALSE.
        ENDIF
        CALL GETTHREADS(NUMTH)
#else
        NUMTH = 1
#endif
c        WRITE(NOUT,90) NUMTH
c90      FORMAT('USING FFTW WITH THREADS: ',I4)

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM

           IF (NSAM.NE.NSAMO) THEN
C             SIZE CHANGED, REESTABLISH PLAN
              IF (PLAN .GT. 0) CALL RFFTW_F77_DESTROY_PLAN(PLAN)
              
#ifdef SP_LIBFFTWMP
              CALL RFFTW_F77_CREATE_PLAN(PLAN,NSAM,
     &              FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE + 
     &                    FFTW_THREADSAFE)
#else
              CALL RFFTW_F77_CREATE_PLAN(PLAN,NSAM,
     &              FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)             
#endif
           ENDIF

C          CHANGE THE SIGN CONVENTION OF FFTW TO THAT OF SPIDER
C          SPIDER FORMAT IMAGINARY PARTS HAVE OPPOSITE SIGNS 
C          TO THAT OF FFTW  

           X(NSAM+1) = X(1)
#ifdef SP_LIBFFTWMP
           CALL RFFTW_F77_THREADS(NUMTH,PLAN,1,X(NSAM+1),-1,1,Y,1,1)
#else
           CALL RFFTW_F77(PLAN,1,X(NSAM+1),-1,1,Y,1,1)               
#endif

           NSAMO = NSAM

C          CHANGE FFTW FORMAT TO SPIDER FFT FORMAT 
           NREM = MOD(NSAM,2)
           LSD = NSAM + 2 - NREM
           DO NM =1, LSD/2
              X((NM - 1) * 2 + 1) = Y(NM)
              IF (NM .EQ. 1) THEN
                 X(2 * NM) = 0.0
              ELSEIF (NM .EQ. LSD/2 .AND. NREM .EQ. 0) THEN
                 X(2 * NM) = 0.0
              ELSE
                  X(2 * NM) = Y(NSAM - NM + 2)
              ENDIF

	   ENDDO 
	   
        ELSE
C          REVERSE TRANSFORM
 
	  
C          CHANGE SPIDER FFT FORMAT TO FFTW FORMAT 
           NREM = MOD(NSAM,2)
           LSD  = NSAM + 2 - NREM
           DO NM =1, LSD/2
              Y(NM) = X((NM-1) * 2 + 1)
              IF ((NM .NE. 1) .AND. 
     &         ((NM .NE. LSD/2) .OR. (NREM .NE. 0))) THEN
                  Y(NSAM - NM + 2) = X(2 * NM)
              ENDIF
           ENDDO
 
           IF (NSAM.NE.NSAMOR) THEN
C             SIZE CHANGED, REESTABLISH PLAN
              IF (PLANR .GT. 0) CALL RFFTW_F77_DESTROY_PLAN(PLANR)
            
#ifdef SP_LIBFFTWMP
              CALL RFFTW_F77_CREATE_PLAN(PLANR,NSAM,
     &              FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE + 
     &                     FFTW_THREADSAFE)
#else 
              CALL RFFTW_F77_CREATE_PLAN(PLANR,NSAM,
     &              FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE )          
#endif            
              NSAMOR = NSAM

           ENDIF

C          CHANGE THE SIGN CONVENTION OF FFTW TO THAT OF SPIDER
C          SPIDER FORMAT IMAGINARY PARTS HAVE OPPOSITE SIGNS 
C          TO THAT OF FFTW  

#ifdef SP_LIBFFTWMP
           CALL RFFTW_F77_THREADS(NUMTH,PLANR,1,Y,1,1,X(NSAM+1),-1,1)
#else
           CALL RFFTW_F77(PLANR,1,Y,1,1,X(NSAM+1),-1,1)
#endif
           X(1) = X(NSAM+1)

C          NORMALIZATION NEEDED
           PIX = 1.0 / (NSAM)

c$omp      parallel do private(i)
           DO I=1,NSAM
              X(I) = X(I) * PIX
           ENDDO   
	   
        ENDIF
#else
	
#if defined(SP_LIBFFTW3) || defined(SP_LIBFFTW3MP)
C       USING FFTW3 LIBRARY CALLS FOR FFT3

#include "FFTW3.INC"

        INTEGER, SAVE :: NSAMO=0
        INTEGER, SAVE :: NSAMOR=0

C       PLAN AND PLANR ARE ACTUALLY POINTERS TO A STRUCTURE 

        INTEGER*8, SAVE :: PLAN=0, PLANR=0

        LOGICAL, SAVE :: INIT=.TRUE.


#ifdef SP_LIBFFTW3MP
        IF (INIT) THEN
C          MUST INITIALIZE THREADS ONCE

           CALL SFFTW_INIT_THREADS(IRTFLG)
           IF (IRTFLG .EQ. 0) THEN
              CALL ERRT(101,'MULTIPLE THREADS FAILED --FFTW3',IER)
              RETURN
           ENDIF

           CALL GETTHREADS(NUMTH) 
	   CALL SFFTW_PLAN_WITH_NTHREADS(NUMTH) 	   
           INIT = .FALSE.
        ENDIF
  
#else
        NUMTH = 1
#endif

C        WRITE(NOUT,90) NUMTH
Cl90      FORMAT('USING FFTW WITH THREADS: ',I4)

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM
 
	   NREM = MOD(NSAM,2)
           LSD = NSAM + 2 - NREM


           IF (NSAM.NE.NSAMO) THEN
C             SIZE CHANGED, REESTABLISH PLAN
              IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)
              
 	      CALL SFFTW_PLAN_DFT_R2C_1D(PLAN,NSAM,X,X,FFTW_ESTIMATE)             
	      
C             INPUT ARRAY SIZE
	      
C             OUTPUT ARRAY SIZE	      
	      
           ENDIF
	       
	   
C          USE FFTW GURU INTERFACE

           CALL SFFTW_EXECUTE_DFT_R2C(PLAN,X,X) 

	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2	

c$omp      parallel do private(i)
	   DO   I = 1,JH	
	      X(2*I) = -X(2*I)           
 	   ENDDO
   
        NSAMO = NSAM

        ELSE
C          REVERSE TRANSFORM
 
C          CHANGE SPIDER FFT FORMAT TO FFTW FORMAT 
           NREM = MOD(NSAM,2)
           LSD  = NSAM + 2 - NREM
	   	   
	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2	

c$omp      parallel do private(i)
	   DO   I = 1,JH	
	      X(2*I) = -X(2*I)           
 	   ENDDO

c$omp      parallel do private(NB)
           DO NB = 1,LDA
              Y(NB) = X(NB)
	   ENDDO

           IF (NSAM.NE.NSAMOR) THEN
C             SIZE CHANGED, REESTABLISH PLAN
              IF (PLANR .GT. 0) CALL SFFTW_DESTROY_PLAN(PLANR)
            
               CALL SFFTW_PLAN_DFT_C2R_1D(PLANR,NSAM,X,
     &  	                               X,FFTW_ESTIMATE)                
	        
	      
              NSAMOR = NSAM

           ENDIF
	   
	   
C          USE FFTW GURU INTERFACE
           CALL SFFTW_EXECUTE_DFT_C2R(PLANR,X,X)  	   	   

C          CHANGE THE SIGN CONVENTION OF FFTW TO THAT OF SPIDER
C          SPIDER FORMAT IMAGINARY PARTS HAVE OPPOSITE SIGNS 
C          TO THAT OF FFTW  

C          NORMALIZATION NEEDED
           PIX = 1.0 / (NSAM)

c$omp      parallel do private(i)
           DO I=1,NSAM
              X(I) = X(I) * PIX
           ENDDO	   
	   
        ENDIF
	
	
#else	

	N = NSAM 

        IF (N.GT.LBUF)  THEN
           INV = 0
C          INSUFFICIENT BUFFER, INCREASE LBUF AND COMPILE SPIDER
           CALL  ERRT(6,'FMRS_1',NE)
           RETURN
        ENDIF

C       INV CAN BE +1 (FORWARD FFT) OR -1 (INVERSE FFT)
        IF (INV)  2,2,1
1       DO I=1,N
           WORK(I) = 0.0
	ENDDO
        CALL FFTMCF(X,WORK,N,N,N,INV)

        IF (MOD(N,2))  12,12,13
12      DO I=N+1,3,-2
           X(I)   = X((I+1)/2)
           X(I+1) = WORK((I+1)/2)
	ENDDO
        X(2)   = 0.0
        X(N+2) = 0.0
        RETURN
13      DO I=N,3,-2
           X(I)   = X(I/2+1)
           X(I+1) = WORK(I/2+1)
	ENDDO
        X(2)=0.0
        RETURN

2       DO I=2,N/2+1
           WORK(I)     = X(2*I)/N
  	   WORK(N-I+2) = -WORK(I)
	ENDDO
        WORK(1) = 0.0

        DO I=1,N/2+1
	   X(I) = X(2*I-1)/N
	ENDDO
        DO I=N,N/2+2,-1
           X(I) = X(N-I+2)
	ENDDO

        CALL FFTMCF(X,WORK,N,N,N,INV)

#endif
#endif
#endif

        END

