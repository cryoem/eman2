C
C++************************************************************************
C
C FMRS_2R.F                     private(j)          OCT 01 ARDEAN LEITH
 
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
C  Order of elements:
C
C--************************************************************************

         SUBROUTINE  FMRS_2R(X,NNNN,NSAM,NROW,INV)

         DIMENSION  X(NNNN,NROW)

	INS=INV*NNNN

        IF (INV .GE. 0)  THEN
c$omp      parallel do private(i),shared(invt) 
	   DO I=1,NROW
	      INVT = INV
	      CALL FMRS_1(X(1,I),NSAM,INVT)
	   ENDDO
           IF (INVT .LE. 0)  THEN
	      INV = INVT
	      RETURN
	   ENDIF
	ENDIF

c$omp   parallel do private(i),shared(invt) 
        DO I=1,NNNN,2
	   INVT = INS
	   CALL FFTMCF(X(I,1),X(I+1,1),NROW,NROW,NROW,INVT)
	ENDDO

        IF (INVT.EQ.0)  THEN
	   INV = 0
	   RETURN
	ENDIF
	IF (INV.GT.0)  RETURN

C       NORMALIZE FOR INVERSE
        Q = 1.0/FLOAT(NROW)
c$omp   parallel do private(i,j)
        DO J=1,NROW
           DO I=1,NNNN
              X(I,J)=X(I,J)*Q
	   ENDDO
	ENDDO 

c$omp   parallel do private(i),shared(invt) 
	DO I=1,NROW
	   INVT=INV
	   CALL  FMRS_1(X(1,I),NSAM,INVT)
	ENDDO

        IF(INVT.LE.0)  INV=INVT

        END
