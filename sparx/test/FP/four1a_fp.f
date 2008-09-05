C++*********************************************************************
C
C    FOUR1A_FP.F                                  
C
C    7/24/00 BIMAL                  ADAPTED TO NEW FOURIER FORMAT               
C                                   OPFILE           NOV 00 ARDEAN LEITH
C                                   OPFILEC          FEB 03 ARDEAN LEITH    
C
C **********************************************************************
C *  SPIDER - MODULAR IMAGE PROCESSING SYSTEM.  AUTHOR: J.FRANK        *
C *  COPYRIGHT (C)1981,1987, WADSWORTH CENTER FOR LABORATORIES AND     *
C *  RESEARCH, NEW YORK STATE DEPARTMENT OF HEALTH, ALBANY, NY 12201.  *
C *  THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO THE CENTER FOR   *
C *  LABORATORIES AND RESEARCH AND ARE NOT TO BE DISCLOSED TO OTHERS OR*
C *  USED FOR PURPOSES OTHER THAN INTENDED WITHOUT WRITTEN APPROVAL OF *
C *  THE CENTER FOR LABORATORIES AND RESEARCH                          *
C **********************************************************************
C
C  FOUR1A_FP
C
C       2D AND 3D IMAGES OF ANY(EVEN/ODD) DIMENSION IS TAKEN AS INPUT
C       AND INTERPOLATED TO ANY DIMENSION. 
C
C IMAGE_PROCESSING_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FOUR1A_FP

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        COMMON /COMMUN/ FILNAM
        CHARACTER (LEN=MAXNAM) ::  FILNAM

        REAL, ALLOCATABLE, DIMENSION(:,:)   :: X  
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: Y     
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: X3 
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Y3 
 
        DATA  LUN1,LUN2/9,10/

        MAXIM   = 0       
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NSAMN   = 0 
        NSLICEN = 0
        MAXIM   = 0 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &    NSAMN,NROWN,NSLICEN,MAXIM,'INTERPOLATED OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CLOSE(LUN1)
           RETURN
        ENDIF

        IF (NSLICE .GT. 1) THEN
C          3D CASE

           LSD  = NSAM +2-MOD(NSAM,2)
           LSDN = NSAMN+2-MOD(NSAMN,2)

           ALLOCATE (X3(LSD,NROW,NSLICE), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'FP, X3',IER)
              GOTO 9001
              RETURN
           ENDIF 

           ALLOCATE (Y3(LSDN,NROWN,NSLICEN), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'FP, Y3',IER)
              GOTO 9000
           ENDIF 

           CALL READV(LUN1,X3,LSD,NROW,NSAM,NROW,NSLICE)

           CALL FINT3(X3,Y3,NSAM,NROW,NSLICE,NSAMN,NROWN,
     &                NSLICEN,LSD,LSDN)

           CALL WRITEV(LUN2,Y3,LSDN,NROWN,NSAMN,NROWN,NSLICEN)

        ELSE
C          2D CASE
           NSLICEN = 1
           LSD     = NSAM +2-MOD(NSAM,2)
           LSDN    = NSAMN+2-MOD(NSAMN,2)

           ALLOCATE (X(LSD,NROW), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'FP, X',IER)
              GOTO 9000
          ENDIF 

           ALLOCATE (Y(LSDN,NROWN), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'FP, Y',IER)
              GOTO 9001
           ENDIF 

           CALL READV(LUN1,X,LSD,NROW,NSAM,NROW,NSLICEN)

           CALL FINT(X,Y,NSAM,NROW,NSAMN,NROWN,LSD,LSDN)

           CALL WRITEV(LUN2,Y,LSDN,NROWN,NSAMN,NROWN,NSLICEN)

        ENDIF

9000   IF (NSLICE .EQ. 1) THEN
           IF (ALLOCATED(X)) DEALLOCATE (X)
           IF (ALLOCATED(Y)) DEALLOCATE (Y)
        ELSE
           IF (ALLOCATED(X3)) DEALLOCATE (X3)
           IF (ALLOCATED(Y3)) DEALLOCATE (Y3)
        ENDIF

9001    CLOSE (LUN2)
        CLOSE (LUN1)

        END
