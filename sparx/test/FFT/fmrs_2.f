

C++************************************************************************
C
C FMRS_2.F                          ADDED FFTW FEB 2000 ARDEAN LEITH
C
C **************************************************************************
C *  SPIDER - MODULAR IMAGE PROCESSING SYSTEM.  AUTHOR: J.FRANK            *
C *  COPYRIGHT (C)1981,1987, WADSWORTH CENTER FOR LABORATORIES AND         *
C *  RESEARCH, NEW YORK STATE DEPARTMENT OF HEALTH, ALBANY, NY 12201.      *
C *  THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO THE CENTER FOR       *
C *  LABORATORIES AND RESEARCH AND ARE NOT TO BE DISCLOSED TO OTHERS OR    *
C *  USED FOR PURPOSES OTHER THAN INTENDED WITHOUT WRITTEN APPROVAL OF     *
C *  THE CENTER FOR LABORATORIES AND RESEARCH   			   *
C **************************************************************************
C
C  FMRS_2(X,NSAM,NROW,INV)
C
C  PARAMETERS:     X       ARRAY                              SENT/RET.
C                  INV     1=REG. FILE, -1= FOURIER FILE       SENT
C
C IMAGE_PROCESSING_ROUTINE
C
C--************************************************************************

	SUBROUTINE  FMRS_2(X,NSAM,NROW,INV)

        INCLUDE 'CMBLOCK.INC'
	DIMENSION  X(*)

#ifdef SP_LIBFFT
C       USING SGI_COMPLIB FOR FFT
	INTEGER, SAVE :: NSAMO=0,NROWO=0
	REAL, DIMENSION(:), POINTER, SAVE :: COEFF

	IF (NSAM.NE.NSAMO .OR. NROW.NE.NROWO) THEN
C          SIZE CHANGED MUST REESTABLISH COEFFICIENTS
	   IF (ASSOCIATED(COEFF))  DEALLOCATE(COEFF)
	   ALLOCATE(COEFF(NSAM+15+2*(NROW+15)),STAT=IRTFLG)
	   IF (IRTFLG.NE.0) CALL ERRT(46,'FT 2, COEFF',IER)

	   CALL SCFFT2DUI(NSAM,NROW,COEFF)
	   NSAMO = NSAM
	   NROWO = NROW
	ENDIF
	LDA = NSAM+2-MOD(NSAM,2)
	IF (INV .GT. 0)  THEN
	   CALL  SCFFT2DU(INV,NSAM,NROW,X,LDA,COEFF)
	ELSE
	   CALL  CSFFT2DU(INV,NSAM,NROW,X,LDA,COEFF)
	   CALL  SSCAL2D(NSAM,NROW,(1.0/FLOAT(NSAM)/FLOAT(NROW)),X,LDA)
	ENDIF
#else
#if defined(SP_LIBFFTW) || defined(SP_LIBFFTWMP)
C       USING FFTW LIBRARY CALLS FOR FFT

#include "FFTW.INC"

	INTEGER, SAVE :: NSAMO=0,  NROWO=0
	INTEGER, SAVE :: NSAMOR=0, NROWOR=0
	LOGICAL, SAVE :: INIT=.TRUE.

C       PLAN AND PLANR ARE ACTUALLY POINTERS TO A STRUCTURE 
#if defined (__osf__) || defined (ia64) || defined (__x86_64__)
        INTEGER*8, SAVE :: PLAN=0, PLANR=0
#else
        INTEGER, SAVE :: PLAN=0, PLANR=0
#endif

#ifdef SP_LIBFFTWMP

C       FOR OMP (MULTIPLE PROCESSORS) FFTW
        IF (INIT) THEN
C          MUST INITIALIZE THREADS ONCE
           CALL FFTW_F77_THREADS_INIT(IRTFLG);
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'MULTIPLE THREADS FAILED',IER)
              RETURN
           ENDIF
           INIT = .FALSE.
        ENDIF
        CALL GETTHREADS(NUMTH)

#else
C       FOR SINGLE PROCESSOR WITH FFTW
        NUMTH = 1
#endif
C        WRITE(NOUT,90) NUMTH
C90      FORMAT('USING FFTW WITH THREADS: ',i4)

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM

	   IF (NSAM.NE.NSAMO .OR. NROW.NE.NROWO) THEN
C             SIZE CHANGED, REESTABLISH PLAN

              IF (PLAN .GT. 0) CALL FFTWND_F77_DESTROY_PLAN(PLAN)

#ifdef SP_LIBFFTWMP
C             FOR OMP (MULTIPLE PROCESSORS)
              CALL RFFTW2D_F77_CREATE_PLAN(PLAN,NSAM,NROW,
     &              FFTW_FORWARD, FFTW_ESTIMATE + FFTW_IN_PLACE +
     &                    FFTW_THREADSAFE)
#else
C             FOR SINGLE PROCESSOR
              CALL RFFTW2D_F77_CREATE_PLAN(PLAN,NSAM,NROW,
     &              FFTW_FORWARD, FFTW_ESTIMATE + FFTW_IN_PLACE)
#endif
           ENDIF

#ifdef SP_LIBFFTWMP
C           FOR OMP (MULTIPLE PROCESSORS)
#if defined (sgi) || (defined (__linux__) && !defined (SP_IFC))
C           SGI & PGI DOES NOT LIKE OBJECT NAMES > 31 CHAR
            CALL RFFTWND_F77_THREADS_ONE_R_TO_C(NUMTH,PLAN,X,0)
#else
            CALL RFFTWND_F77_THREADS_ONE_REAL_TO_COMPLEX(NUMTH,PLAN,X,0)
#endif
#else
C          FOR SINGLE PROCESSOR
           CALL RFFTWND_F77_ONE_REAL_TO_COMPLEX(PLAN,X,0)
#endif

C          CHANGE FFTW FORMAT TO SPIDER FFT FORMAT 
C          SPIDER FORMAT IMAGINARY PARTS HAVE OPPOSITE SIGNS 
C          AS THAT OF FFTW 

	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2
	   	   
c$omp      parallel do private(i)
	   DO   I = 1,JH*NROW	
	      X(2*I) = -X(2*I)           
 	   ENDDO
	   	

	   NSAMO = NSAM
	   NROWO = NROW
        ELSE
C          REVERSE TRANSFORM

C          CHANGE SPIDER FFT FORMAT TO FFTW FORMAT
C          IMAGINARY PARTS HAVE OPPOSITE SIGNS AS THAT OF FFTW 


	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2
	   	   
c$omp      parallel do private(i)
	   DO   I = 1,JH*NROW	
	      X(2*I) = -X(2*I)           
 	   ENDDO

	   IF (NSAM.NE.NSAMOR .OR. NROW.NE.NROWOR) THEN
C             SIZE CHANGED, REESTABLISH PLAN

              IF (PLANR .GT. 0) CALL FFTWND_F77_DESTROY_PLAN(PLANR)

#ifdef SP_LIBFFTWMP
C              FOR OMP (MULTIPLE PROCESSORS)
               CALL RFFTW2D_F77_CREATE_PLAN(PLANR,NSAM,NROW,
     &              FFTW_BACKWARD, FFTW_ESTIMATE + FFTW_IN_PLACE +
     &                    FFTW_THREADSAFE)
#else
              CALL RFFTW2D_F77_CREATE_PLAN(PLANR,NSAM,NROW,
     &              FFTW_BACKWARD, FFTW_ESTIMATE + FFTW_IN_PLACE)
#endif
	      NSAMOR = NSAM
	      NROWOR = NROW
           ENDIF

#ifdef SP_LIBFFTWMP
C          FOR OMP (MULTIPLE PROCESSORS)
#if defined (sgi) || (defined (__linux__) && !defined (SP_IFC))
C          SGI & PGI DOES NOT LIKE OBJECT NAMES > 31 CHAR
           CALL RFFTWND_F77_THREADS_ONE_C_TO_R(NUMTH,PLANR,X,0)
#else
           CALL RFFTWND_F77_THREADS_ONE_COMPLEX_TO_REAL(NUMTH,PLANR,X,0)
#endif
#else
C          FOR SINGLE PROCESSOR
           CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(PLANR,X,0)
#endif

C          SCALING NEEDED
           PIX = 1.0 / (NSAM * NROW)
          

c$omp      parallel do private(i)
           DO I=1,LDA * NROW
              X(I) = X(I) * PIX
           ENDDO

        ENDIF

#else

#if defined(SP_LIBFFTW3) || defined(SP_LIBFFTW3MP)
C       USING FFTW3 LIBRARY CALLS FOR FFT
#include "FFTW3.INC"


C       USING FFTW3 LIBRARY CALLS FOR FFT

	INTEGER, SAVE :: NSAMO=0,  NROWO=0
	INTEGER, SAVE :: NSAMOR=0, NROWOR=0
	LOGICAL, SAVE :: INIT=.TRUE.

C       PLAN AND PLANR ARE ACTUALLY POINTERS TO A STRUCTURE 

        INTEGER*8, SAVE :: PLAN=0, PLANR=0
	REAL, ALLOCATABLE, DIMENSION(:) :: XY
	REAL, ALLOCATABLE, DIMENSION(:) :: AB	

#ifdef SP_LIBFFTW3MP


C       FOR OMP (MULTIPLE PROCESSORS) FFTW
        IF (INIT) THEN

C          MUST INITIALIZE THREADS ONCE
           CALL SFFTW_INIT_THREADS(IRTFLG)
	   IF (IRTFLG .EQ. 0) THEN
              CALL ERRT(101,'MULTIPLE THREADS FAILED',IER)
              RETURN
           ENDIF
	   
	   CALL GETTHREADS(NUMTH)
	   CALL SFFTW_PLAN_WITH_NTHREADS(NUMTH) 

           INIT = .FALSE.
        ENDIF

#else
C       FOR SINGLE PROCESSOR WITH FFTW
        NUMTH = 1
#endif

C         WRITE(NOUT,90) NUMTH
C90       FORMAT('USING FFTW WITH THREADS: ',i4)

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM

	   LDA = NSAM+2-MOD(NSAM,2)

	   IF (NSAM.NE.NSAMO .OR. NROW.NE.NROWO) THEN
C             SIZE CHANGED, REESTABLISH XY AND THE PLAN

	      IF (ALLOCATED(XY))  DEALLOCATE(XY)
	   
	      ALLOCATE(XY(LDA*NROW),STAT=IRTFLG)
	      IF (IRTFLG.NE.0) THEN
	         CALL ERRT(46,'FMRS_2, XY',IER)
	         RETURN
              ENDIF
	   
C           COPYING IMAGE TO BUFFER SINCE IMAGE WILL BE OVERWRITTEN DURING
C           MAKING THE PLAN WITH FFTW_MEASURE FLAG
c$omp       parallel do private(i)
	      DO   I = 1,LDA*NROW	
	         XY(I) = X(I)           
 	      ENDDO	

              N1 = NSAM
	      N2 = NROW 
              IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)
              CALL SFFTW_PLAN_DFT_R2C_2D(PLAN,N1,N2,X,X,FFTW_MEASURE)

C           COPYING IMAGE FROM BUFFER SINCE IMAGE HAS BEEN OVERWRITTEN DURING
C           MAKING THE PLAN WITH FFTW_MEASURE FLAG
c$omp       parallel do private(i)
	      DO   I = 1,LDA*NROW	
	         X(I) = XY(I)           
 	      ENDDO	
	 
           ENDIF
 
	   
C         USE FFTW GURU INTERFACE
           CALL SFFTW_EXECUTE_DFT_R2C(PLAN,X,X) 

C         FREE MEMORY	   
	   IF (ALLOCATED(XY))  DEALLOCATE(XY)

C          CHANGE FFTW FORMAT TO SPIDER FFT FORMAT 
C          SPIDER FORMAT IMAGINARY PARTS HAVE OPPOSITE SIGNS 
C          AS THAT OF FFTW 

	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2
	      	   
c$omp      parallel do private(i)
	   DO   I = 1,JH*NROW	
	      X(2*I) = -X(2*I)           
 	   ENDDO
	   	

	   NSAMO = NSAM
	   NROWO = NROW
        ELSE
C          REVERSE TRANSFORM
C          CHANGE SPIDER FFT FORMAT TO FFTW FORMAT
C          IMAGINARY PARTS HAVE OPPOSITE SIGNS AS THAT OF FFTW 


	   LDA = NSAM+2-MOD(NSAM,2)
           JH  = LDA/2
	   	   
c$omp      parallel do private(i)
	   DO   I = 1,JH*NROW	
	      X(2*I) = -X(2*I)           
 	   ENDDO

	   IF (NSAM.NE.NSAMOR .OR. NROW.NE.NROWOR) THEN
C             SIZE CHANGED, REESTABLISH AB AND THE PLAN

              IF (PLANR .GT. 0) CALL SFFTW_DESTROY_PLAN(PLANR)


	      IF (ALLOCATED(AB))  DEALLOCATE(AB)
	   
	      ALLOCATE(AB(LDA*NROW),STAT=IRTFLG)
	      IF (IRTFLG.NE.0) THEN
	         CALL ERRT(46,'FMRS_2, AB',IER)
	         RETURN
              ENDIF
	   
C           COPYING IMAGE TO BUFFER SINCE IMAGE WILL BE OVERWRITTEN DURING
C           MAKING THE PLAN WITH FFTW_MEASURE FLAG
c$omp       parallel do private(i)
	      DO   I = 1,LDA*NROW	
	         AB(I) = X(I)           
 	      ENDDO	
      
              N1 = NSAM
	      N2 = NROW 
              CALL SFFTW_PLAN_DFT_C2R_2D(PLANR,N1,N2,X,X,FFTW_MEASURE)
	      
C           COPYING IMAGE FROM BUFFER SINCE IMAGE HAS BEEN OVERWRITTEN DURING
C           MAKING THE PLAN WITH FFTW_MEASURE FLAG
c$omp       parallel do private(i)
	      DO   I = 1,LDA*NROW	
	         X(I) = AB(I)           
 	      ENDDO
	      	      
	      NSAMOR = NSAM
	      NROWOR = NROW
           ENDIF

C         USE FFTW GURU INTERFACE   
           CALL SFFTW_EXECUTE_DFT_C2R(PLANR,X,X) 
	   
C         FREE MEMORY	   
	   IF (ALLOCATED(AB))  DEALLOCATE(AB)

C          SCALING NEEDED
           PIX = 1.0 / (NSAM * NROW)
          

c$omp      parallel do private(i)
           DO I=1,LDA * NROW
              X(I) = X(I) * PIX
           ENDDO

        ENDIF


#else
C       USING SPIDER CODE FOR FFT
C       HAVE TO CHANGE NSAM
	LDA = NSAM+2-MOD(NSAM,2)

        CALL FMRS_2R(X,LDA,NSAM,NROW,INV)
#endif
#endif
#endif
	END

