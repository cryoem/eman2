C
C++************************************************************************
C
C $$ FMRS_3DR.FOR
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--************************************************************************
C
C $$ FMRS_3DR.FOR
C
	SUBROUTINE FMRS_3DR(LUN1,LUN2,
     &		NNNN,NSAM,NROW,NSLICE,INV)	

	INCLUDE 'CMBLOCK.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:) :: A, A2, B, H1, H2
        REAL, ALLOCATABLE, DIMENSION(:)   :: WORK

	LOGICAL IFNS,IFND


        NMN = (MAX0(NSAM,2*NROW,2*NSLICE))
        
	ALLOCATE (WORK(NMN), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, WORK',IER)
           RETURN
        ENDIF	
 	ALLOCATE (A(NNNN,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, A',IER)
           DEALLOCATE (WORK)
           RETURN
        ENDIF		
	ALLOCATE (A2(NSAM,NROW), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, A2',IER)
           DEALLOCATE (WORK,A)
           RETURN
        ENDIF		
	ALLOCATE (B(NNNN,NSLICE), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, B',IER)
           DEALLOCATE (WORK,A,A2)
           RETURN
        ENDIF		
	ALLOCATE (H1(NROW,NSLICE), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, H1',IER)
           DEALLOCATE (WORK,A,A2,B)
           RETURN
        ENDIF		
	ALLOCATE (H2(NROW,NSLICE), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, H2',IER)
           DEALLOCATE (WORK,A,A2,B,H1)
           RETURN
        ENDIF		



	NDR=INV*NNNN*NROW
C
	IF(INV) 2,2,1
C
1	CONTINUE
c******************************************************************************
c********forward************forward**************forward***********************
c******************************************************************************
C
C
        DO    I=1,NSLICE
           DO    I1=1,NROW
              K=(I-1)*NROW+I1
              CALL  REDLIN(LUN1,A(1,I1),NSAM,K)
	   ENDDO

           CALL  FMRS_2(A,NSAM,NROW,INV)
           IF(INV.EQ.0)  THEN
 	      DEALLOCATE (A, A2, B, H1, H2, WORK)               
              RETURN
           ENDIF

           DO    I1=1,NROW
              K=(I-1)*NROW+I1
              CALL  WRTLIN(LUN2,A(1,I1),NNNN,K)
	   ENDDO
	ENDDO

	CLOSE(LUN1)

        DO    I=1,NROW
           DO    I1=1,NSLICE
              K=(I1-1)*NROW+I
              CALL  REDLIN(LUN2,B(1,I1),NNNN,K)
	   ENDDO

c$omp parallel do private(i1),shared(nnnnt)
           DO    I1=1,NNNN,2
	      NNNNT=NNNN
	CALL  FFTMCF(B(I1,1),B(I1+1,1),NSLICE,NSLICE,NSLICE,NNNNT)
	   ENDDO
	   IF(NNNNT.EQ.0)  THEN
	      INV=0
	      DEALLOCATE (A, A2, B, H1, H2, WORK)
	      RETURN
	   ENDIF

           DO  I1=1,NSLICE
              K=(I1-1)*NROW+I
              CALL  WRTLIN(LUN2,B(1,I1),NNNN,K)
	   ENDDO
	ENDDO
	DEALLOCATE (A, A2, B, H1, H2, WORK)
	RETURN
c******************************************************************************
c************inverse**********inverse************inverse***********************
c******************************************************************************

2	CONTINUE
	IFNS=MOD(NSAM,2).EQ.0
	IFND=MOD(NROW,2).EQ.0
C       NORMALIZE FOR INVERSE
	Q=1/FLOAT(NSLICE)

        DO    I=1,NROW
           DO    I1=1,NSLICE
              K=(I1-1)*NROW+I
              CALL  REDLIN(LUN1,B(1,I1),NNNN,K)
	   ENDDO

c$omp parallel do private(i1),shared(nnnnt)
           DO    I1=1,NNNN,2
	      NNNNT=-NNNN
	CALL  FFTMCF(B(I1,1),B(I1+1,1),NSLICE,NSLICE,NSLICE,NNNNT)
	   ENDDO
	   IF(NNNNT.EQ.0)  THEN
	      INV=0
	      DEALLOCATE (A, A2, B, H1, H2, WORK)
	      RETURN
	   ENDIF

	   IF(I.LT.NROW/2+1)  THEN
              DO    K=1,NSLICE
                 H1(2*(I-1)+1,K)=B(1,K)
                 H1(2*(I-1)+2,K)=B(2,K)
              ENDDO
              IF(IFNS)  THEN
                 DO    K=1,NSLICE
                    H2(2*(I-1)+1,K)=B(NNNN-1,K)
                    H2(2*(I-1)+2,K)=B(NNNN,K)
                 ENDDO
              ENDIF
           ELSEIF(I.EQ.NROW/2+1)  THEN
              IF(IFND)  THEN
                 DO    K=1,NSLICE
                    H1(2,K)=B(1,K)
                 ENDDO
              ELSE
                 DO    K=1,NSLICE
                    H1(2*(I-1)+1,K)=B(1,K)
                    H1(2,K)=B(2,K)
                 ENDDO
              ENDIF
              IF(IFNS)  THEN
                 IF(IFND)  THEN
                    DO    K=1,NSLICE
                       H2(2,K)=B(NNNN-1,K)
                    ENDDO
                 ELSE
                    DO    K=1,NSLICE
                       H2(2*(I-1)+1,K)=B(NNNN-1,K)
                       H2(2,K)=B(NNNN,K)
                    ENDDO
                 ENDIF
              ENDIF
	   ELSE
C             REMAINING ENTRIES ARE MIRRORED AND ARE SKIPPED
	   ENDIF

	   IF(.NOT.IFNS)  THEN
	      DO    I1=1,NSLICE
	         B(2,I1)=B(NNNN,I1)
	      ENDDO
	   ENDIF

           DO    I1=1,NSLICE
              K=(I1-1)*NROW+I
              CALL  WRTLIN(LUN2,B(1,I1),NSAM,K)
	   ENDDO
	ENDDO

C	CLOSE(LUN1)

        DO    I=1,NSLICE
           DO    I1=1,NROW
              K=(I-1)*NROW+I1
              CALL  REDLIN(LUN2,A2(1,I1),NSAM,K)
	   ENDDO

	   DO    I1=1,NROW
	      A2(1,I1)=H1(I1,I)
	   ENDDO
	   IF(IFNS)  THEN
	      DO    I1=1,NROW
	         A2(2,I1)=H2(I1,I)
	      ENDDO
	   ENDIF
c$omp parallel do  private(k,i1)
	   DO    K=1,NROW
	      DO    I1=1,NSAM
	         A2(I1,K)=A2(I1,K)*Q
	      ENDDO
	   ENDDO

           CALL  FMR_2(A2,NSAM,NROW,WORK,INV)
           IF(INV.EQ.0)  THEN
              CLOSE(LUN1)
 	      DEALLOCATE (A, A2, B, H1, H2, WORK)          
              RETURN
           ENDIF
           DO    I1=1,NROW
              K=(I-1)*NROW+I1
              CALL  WRTLIN(LUN2,A2(1,I1),NSAM,K)
	   ENDDO
	ENDDO
	CLOSE(LUN1)
	
	DEALLOCATE (A, A2, B, H1, H2, WORK)	
	END
