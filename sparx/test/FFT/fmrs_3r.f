C
C++************************************************************************
C
C $$ FMRS_3R.FOR
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
C  For order of elements see fmr_1.
C
C
C--************************************************************************
C
C $$ FMRS_3R.FOR
C
	SUBROUTINE FMRS_3R(A,NNNN,NSAM,NROW,NSLICE,INV)	
	DIMENSION A(NNNN,NROW,NSLICE)
c
	NDR=INV*NNNN*NROW
C
	IF(INV.GE.0)  THEN
	 DO   I=1,NSLICE
	  CALL FMRS_2(A(1,1,I),NSAM,NROW,INV)
	  IF(INV.EQ.0)   RETURN
	 ENDDO
	ENDIF
C
c$omp parallel do private(j,i),shared(ndrt)
	DO    J=1,NROW
	DO    I=1,NNNN-1,2
	NDRT=NDR
	CALL FFTMCF(A(I,J,1),A(I+1,J,1),NSLICE,NSLICE,NSLICE,NDRT)
	ENDDO	
	ENDDO
	IF(NDRT.EQ.0)  THEN
	INV=0
	RETURN
	ENDIF
	IF(INV.GT.0)  RETURN
C NORMALIZE FOR INVERSE
	Q=1.0/FLOAT(NSLICE)
c$omp parallel do private(k,j,i)
	DO    K=1,NSLICE
	DO    J=1,NROW
	DO    I=1,NNNN
	A(I,J,K)=A(I,J,K)*Q
	ENDDO
	ENDDO
	ENDDO
	DO   I=1,NSLICE
	CALL FMRS_2(A(1,1,I),NSAM,NROW,INV)
	IF(INV.EQ.0)   RETURN
	ENDDO
	END
