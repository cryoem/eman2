C++************************************************************************
C
C FMR_2.F
C
C **************************************************************************
C *  SPIDER - MODULAR IMAGE PROCESSING SYSTEM.  AUTHOR: J.FRANK            *
C *  COPYRIGHT (C)1981,1987, WADSWORTH CENTER FOR LABORATORIES AND         *
C *  RESEARCH, NEW YORK STATE DEPARTMENT OF HEALTH, ALBANY, NY 12201.      *
C *  THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO THE CENTER FOR       *
C *  LABORATORIES AND RESEARCH AND ARE NOT TO BE DISCLOSED TO OTHERS OR    *
C *  USED FOR PURPOSES OTHER THAN INTENDED WITHOUT WRITTEN APPROVAL OF     *
C *  THE CENTER FOR LABORATORIES AND RESEARCH                              *
C **************************************************************************
C  For order of elements see fmr_1.
C
CC IMAGE_PROCESSING_ROUTINE                                                                     
C        0         2         3         4         5         6         7 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--************************************************************************

#ifdef SP_MP
        SUBROUTINE  FMR_2(X,NSAM,NROW,DUMMY,INV)

        PARAMETER (LENBUF=8200)
        DIMENSION  X(NSAM,NROW),WORK(LENBUF),DUMMY(*)
        LOGICAL*1  IFND

        IF (MAX0(NSAM,2*NROW).GT.LENBUF)  THEN
           CALL  ERRT(6,' FMR_2 ',NE)
           INV = 0
           RETURN
        ENDIF
#else
        SUBROUTINE  FMR_2(X,NSAM,NROW,WORK,INV)
        DIMENSION   X(NSAM,NROW),WORK(*)
        LOGICAL*1   IFND
#endif

        IFND=MOD(NSAM,2).EQ.0
        IF (IFND)  THEN
           LBD=2
        ELSE
           LBD=1
        ENDIF
        INS=INV*NSAM

C       work(max0(nsam,2*nrow))

        IF (INV.GE.0)  THEN
c$omp      parallel do private(i,work),shared(invt) 
           DO I=1,NROW
              INVT=INV
              CALL  FMR_1(X(1,I),NSAM,WORK,INVT)
           ENDDO
           IF (INVT.LE.0)  THEN
              INV=INVT
              RETURN
           ENDIF
        ENDIF

c$omp   parallel do  private(i,j,work),shared(invt) 
        DO I=1,LBD
           DO J=1,NROW
              WORK(NROW+J)=X(I,J)
           ENDDO
           INVT = INV
           CALL  FMR_1(WORK(NROW+1),NROW,WORK,INVT)
           DO J=1,NROW
              X(I,J)=WORK(NROW+J)
           ENDDO
        ENDDO
        IF (INVT.EQ.0)  THEN
           INV=0
           RETURN
        ENDIF
c$omp   parallel do  private(i),shared(invt) 
        DO I=3,NSAM-1,2
           INVT=INS
           CALL  FFTMCF(X(I,1),X(I+1,1),NROW,NROW,NROW,INVT)
        ENDDO
        IF (INVT.EQ.0)  THEN
           INV=0
           RETURN
        ENDIF
        IF (.NOT.IFND)  CALL FFTMCF
     &          (X(NSAM,1),X(2,1),NROW,NROW,NROW,INS)
        IF (INV.EQ.1)  RETURN

C       NORMALIZE FOR INVERSE
        Q=1.0/FLOAT(NROW)
c$omp   parallel do  private(j,i) 
        DO J=1,NROW
           DO I=LBD+1,NSAM
              X(I,J)=X(I,J)*Q
           ENDDO
        ENDDO

c$omp   parallel do  private(i,work),shared(invt) 
        DO I=1,NROW
           INVT = INV
           CALL  FMR_1(X(1,I),NSAM,WORK,INVT)
        ENDDO
        IF (INVT.LE.0)  INV = INVT
        END
