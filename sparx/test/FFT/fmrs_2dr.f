C++************************************************************************
C
C $$ FMRS_2DR.FOR
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
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--************************************************************************
C
C $$ FMRS_2DR.FOR
C
        SUBROUTINE  FMRS_2DR(LUN1,LUN2,LR,NNNN,NSAM,NROW,INV)
        REAL, ALLOCATABLE, DIMENSION(:,:) :: BUF
        REAL, ALLOCATABLE, DIMENSION(:) :: X       
	LOGICAL  IFNS 
		

	ALLOCATE (BUF(LR,NROW), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, BUF',IER)
           RETURN
        ENDIF
        	
	ALLOCATE (X(NNNN), STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'FT, X',IER)
           DEALLOCATE (BUF)
           RETURN
        ENDIF	
	
C       NUMBER OF CHUNKS
	NC=NNNN/LR
	IF(MOD(NNNN,LR).NE.0) NC=NC+1
	LRE=NNNN-(NC-1)*LR
	INS=INV*LR

        IF(INV.GT.0)  THEN
	   DO    J=1,NROW
	      CALL  REDLIN(LUN1,X,NSAM,J)
              CALL  FMRS_1(X,NSAM,INV)
              IF(INV.EQ.0)  THEN
                 DEALLOCATE (X,BUF)               
                 RETURN
              ENDIF
	      DO    I=1,LR
	         BUF(I,J)=X(I)
	      ENDDO
	      CALL  WRTLIN(LUN2,X,NNNN,J)
	   ENDDO
c	   CLOSE(LUN1)
c$omp parallel do private(i),shared(invt)
           DO    I=1,LR,2
	      INVT=INS
	      CALL  FFTMCF(BUF(I,1),BUF(I+1,1),NROW,NROW,NROW,INVT)
	   ENDDO
           IF(INVT.EQ.0)  THEN
              INV=0
              DEALLOCATE (X,BUF)   
              RETURN
           ENDIF
           
	   IF(NC.GT.2)  THEN
C             DO FULL CHUNKS
	      DO    LC=2,NC-1
C				print  *,lc,nc
	         DO    J=1,NROW
	            CALL  REDLIN(LUN2,X,NNNN,J)
	            DO    I=1,LR
	               X(I+(LC-2)*LR)=BUF(I,J)
	               BUF(I,J)=X(I+(LC-1)*LR)
	            ENDDO
	            CALL  WRTLIN(LUN2,X,NNNN,J)
	         ENDDO
c$omp parallel do private(i),shared(invt)
                 DO    I=1,LR,2
	            INVT=INS
	CALL  FFTMCF(BUF(I,1),BUF(I+1,1),NROW,NROW,NROW,INVT)
	         ENDDO
                 IF(INVT.EQ.0)  THEN
                    INV=0
                    DEALLOCATE (X,BUF)   
                    RETURN
                 ENDIF
	      ENDDO
	   ENDIF

C          DO THE LAST, PROBABLY SHORTER CHUNK
	   DO    J=1,NROW
	      CALL  REDLIN(LUN2,X,NNNN,J)
	      DO    I=1,LR
	         X(I+(NC-2)*LR)=BUF(I,J)
	      ENDDO
	      DO    I=1,LRE
	         BUF(I,J)=X(I+(NC-1)*LR)
	      ENDDO
	      CALL  WRTLIN(LUN2,X,NNNN,J)
	   ENDDO
c$omp parallel do private(i),shared(invt)
           DO    I=1,LRE,2
	      INVT=INS
	CALL  FFTMCF(BUF(I,1),BUF(I+1,1),NROW,NROW,NROW,INVT)
	   ENDDO 
           IF(INVT.EQ.0)  THEN
              INV=0
              DEALLOCATE (X,BUF)   
              RETURN
           ENDIF
	   DO    J=1,NROW
	      CALL  REDLIN(LUN2,X,NNNN,J)
	      DO    I=1,LRE
	         X(I+(NC-1)*LR)=BUF(I,J)
	      ENDDO
	      CALL  WRTLIN(LUN2,X,NNNN,J)
	   ENDDO
c	   CLOSE(LUN2)
           DEALLOCATE (X,BUF)   
	   RETURN	
c          INVERSE
	ENDIF
C       BEGIN HERE WHEN INV<=0
	IFNS=MOD(NSAM,2).EQ.0
C       NORMALIZE FOR INVERSE
        Q=1.0/FLOAT(NROW)
C       DO THE FIRST CHUNK, HAVE TO COMPRESS
        DO    J=1,NROW
	   CALL  REDLIN(LUN1,X,NNNN,J)
           DO    I=1,NNNN
              X(I)=X(I)*Q
	   ENDDO
	   DO    I=1,2
	      BUF(I,J)=X(I)
	      BUF(I+2,J)=X(I+NNNN-2)
	   ENDDO
	   IF(LR.GT.4)  THEN
	      DO    I=5,LR
	         BUF(I,J)=X(I-2)
	      ENDDO
	   ENDIF
	   X(1)=0.0
	   X(2)=0.0
	   CALL  WRTLIN(LUN2,X,NSAM,J)
	ENDDO
c	CLOSE(LUN1)
c$omp parallel do private(i),shared(invt)
        DO    I=1,LR,2
	   INVT=INS
	CALL  FFTMCF(BUF(I,1),BUF(I+1,1),NROW,NROW,NROW,INVT)
	ENDDO
        IF(INVT.EQ.0)  THEN
           INV=0
           DEALLOCATE (X,BUF)   
           RETURN
        ENDIF

	IF(NC.GT.2)  THEN
	   LRC=LR
	ELSE
	   LRC=LRE
	ENDIF
	DO    LC=2,NC
	  DO     J=1,NROW
	     CALL  REDLIN(LUN2,X,NSAM,J)
	     IF(LC.EQ.2)  THEN
C               PUT ONLY REAL PARTS, IMAGINARY ARE ZERO
                X(1)=BUF(1,J)
		IF(IFNS)  THEN
		   X(2)=BUF(3,J)
		ELSE
		   X(2)=BUF(4,J)
		   X(NSAM)=BUF(3,J)
		ENDIF
		IF(LR.GT.4)  THEN
		   DO    I=5,LRC
		      X(I-2)=BUF(I,J)
		   ENDDO
		ENDIF
		DO    I=1,LRC
		   BUF(I,J)=X(I+LRC-2)
		ENDDO
	     ELSE
		 DO    I=1,LRC
		    X(I+(LC-2)*LR-2)=BUF(I,J)
		    BUF(I,J)=X(I+(LC-1)*LR-2)
                 ENDDO
	      ENDIF
	      CALL  WRTLIN(LUN2,X,NSAM,J)
	   ENDDO
c$omp parallel do private(i),shared(invt)
           DO    I=1,LRC,2
	      INVT=INS
	CALL  FFTMCF(BUF(I,1),BUF(I+1,1),NROW,NROW,NROW,INVT)
	   ENDDO
           IF(INVT.EQ.0)  THEN
              INV=0
              DEALLOCATE (X,BUF)   
              RETURN
           ENDIF
	ENDDO
C       DO THE LAST, PROBABLY SHORTER CHUNK
	DO    J=1,NROW
	      CALL  REDLIN(LUN2,X,NSAM,J)
	   DO    I=1,LRE
	      X(I+(NC-1)*LR-2)=BUF(I,J)
	   ENDDO
	   CALL  WRTLIN(LUN2,X,NSAM,J)
	ENDDO
C
	DO    J=1,NROW
	   CALL  REDLIN(LUN2,X,NSAM,J)
	   INV=-1
           CALL  FMR_1(X,NSAM,BUF,INV)
	   CALL  WRTLIN(LUN2,X,NSAM,J)
	ENDDO
c	CLOSE(LUN2)
	DEALLOCATE (BUF,X)
        END
