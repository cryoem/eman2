C ++********************************************************************
C    FINT3.FOR                                                         *
C    9/11/00 BIMAL                       ADAPTED TO NEW FOURIER FORMAT *
C                                                                      *                                                                      *
C **********************************************************************
C * SPIDER - MODULAR IMAGE PROCESSING SYSTEM.    AUTHOR: J.FRANK       *
C * COPYRIGHT (C)1985, 1999. HEALTH RESEARCH INCORPORATED (HRI),       *
C * ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                   *
C * THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO HRI AND ARE NOT   *
C * TO BE DISCLOSED TO OTHERS OR USED FOR PURPOSES OTHER THAN INTENDED *
C * WITHOUT WRITTEN APPROVAL OF HRI.                                   *
C **********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:
C                                                                      *
C IMAGE_PROCESSING_ROUTINE  
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************


        SUBROUTINE FINT3(X3,Y3,NSAM,NROW,NSLICE,NSAMN,NROWN,
     &                    NSLICEN,LSD,LSDN)



        DIMENSION   X3(LSD,NROW,NSLICE)
        DIMENSION   Y3(LSDN,NROWN,NSLICEN) 

C       TO KEEP THE EXACT VALUES ON THE PREVIOUS GRID YOU SHOULD USE
C       SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED.

        SQ2     = SQRT(2.0)


        ANORM    = FLOAT(NSAMN*NROWN*NSLICEN)/FLOAT(NSAM*NROW*NSLICE)
        NM       = MOD(NSAM,2)

        INV=+1
        CALL  FMRS_3(X3,NSAM,NROW,NSLICE,INV)

C       NORMALIZATION REQUIRED
        X3 = X3*ANORM
	IF(NSAMN>=NSAM)  THEN
        INSAM    = NSAMN - NSAM
        INROW    = NROWN - NROW
        INSLICE  = NSLICEN - NSLICE
C       ZERO PADDING IS DONE
        Y3 = 0.0

        DO K =1, NSLICE/2 +1
           DO J = 1, NROW/2+1
              DO I =1,LSD
                 Y3(I,J,K) = X3(I,J,K) 
              ENDDO
           ENDDO
        ENDDO

        DO K =1, NSLICE/2 +1
           DO J = NROW/2+2+INROW, NROWN
              DO I =1,LSD
                 Y3(I,J,K) = X3(I,J-INROW,K) 
              ENDDO
           ENDDO
        ENDDO

        DO K = NSLICE/2+2+INSLICE, NSLICEN
           DO J = 1, NROW/2+1
              DO I =1,LSD
                 Y3(I,J,K) = X3(I,J,K-INSLICE) 
              ENDDO
           ENDDO
        ENDDO
	
        DO K = NSLICE/2+2+INSLICE, NSLICEN
           DO J = NROW/2+2+INROW, NROWN
              DO I =1,LSD
                 Y3(I,J,K) = X3(I,J-INROW,K-INSLICE) 
              ENDDO
           ENDDO
        ENDDO

C       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
C       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2 TH
C       ELEMENT.

        IF (NM .EQ. 0 .AND. INSAM .NE. 0) THEN
           DO K = 1,NSLICEN
              DO J = 1,NROWN
                 Y3(NSAM+1,J,K) = (1/SQ2)*Y3(NSAM+1,J,K)
                 Y3(NSAM+2,J,K) = (1/SQ2)*Y3(NSAM+2,J,K)
              ENDDO
           ENDDO

           DO K = 1,NSLICEN
       	      DO I =1,LSD
                 Y3(I,NROW/2+1+INROW,K) = Y3(I,NROW/2+1,K)/SQ2
                 Y3(I,NROW/2+1,K)       = Y3(I,NROW/2+1,K)/SQ2
              ENDDO
           ENDDO

           DO J = 1,NROWN
              DO I =1,LSD
                 Y3(I,J,NSLICE/2+1+INSLICE) = Y3(I,J,NSLICE/2+1)/SQ2
                 Y3(I,J,NSLICE/2+1)         = Y3(I,J,NSLICE/2+1)/SQ2
              ENDDO
           ENDDO
        ENDIF

	ELSE
	Y3(:,1:NROWN/2+1,1:NSLICEN/2+1)=
     &			X3(1:LSDN,1:NROWN/2+1,1:NSLICEN/2+1)
	DO  J=NROWN,NROWN/2+2,-1
	Y3(:,J,1:NSLICEN/2+1)=X3(1:LSDN,J+NROW-NROWN,1:NSLICEN/2+1)
	ENDDO

	DO  K=NSLICEN,NSLICEN/2+2,-1
	Y3(:,1:NROWN/2+1,K)=X3(1:LSDN,1:NROWN/2+1,K+NSLICE-NSLICEN)
	DO  J=NROWN,NROWN/2+2,-1
	Y3(:,J,K)=X3(1:LSDN,J+NROW-NROWN,K+NSLICE-NSLICEN)
	ENDDO
	ENDDO

	ENDIF

        INV= -1
        CALL  FMRS_3(Y3,NSAMN,NROWN,NSLICEN,INV)

        END
