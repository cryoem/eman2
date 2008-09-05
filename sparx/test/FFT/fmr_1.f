C++************************************************************************
C
C $$ FMR_1.FOR
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
C  1D real mixed radix FFT.
C   Order of elements:
C    for N even
C  R(0), R(N/2), R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1)
C
C    for N odd
C  R(0), I(N/2), R(1), I(1), R(2), I(2), ....., R(N/2-1), I(N/2-1), R(N/2) 
C
C
C IMAGE_PROCESSING_ROUTINE                                                                     
C        0         2         3         4         5         6         7 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--************************************************************************

         SUBROUTINE  FMR_1(X,N,WORK,INV)

         DIMENSION  X(N),WORK(N)
C INV CAN BE +1 (FORWARD FFT) OR -1 (INVERSE FFT)

         IF(INV)  2,2,1
1        DO    I=1,N
            WORK(I)=0.0
         ENDDO
         CALL FFTMCF(X,WORK,N,N,N,INV)
         IF(MOD(N,2))  12,12,13
12       WORK(1)=X(N/2+1)
         DO    I=N-1,3,-2
            X(I)=X((I+1)/2)
            X(I+1)=WORK((I+1)/2)
         ENDDO
         X(2)=WORK(1)
         RETURN
13       X(N)=X(N/2+1)
         DO    I=N-2,3,-2
            X(I)=X(I/2+1)
            X(I+1)=WORK(I/2+1)
         ENDDO
         X(2)=WORK(N/2+1)
         RETURN

2        IF(MOD(N,2))  22,22,23
22       WORK(N/2+1)=0.0
         DO    I=2,N/2
            WORK(I)=X(2*I)/N
            WORK(N-I+2)=-WORK(I)
         ENDDO
         WORK(1)=X(2)/N
         X(1)=X(1)/N
         DO    I=2,N/2
            X(I)=X(2*I-1)/N
         ENDDO
         DO    I=N/2+2,N
            X(I)=X(N-I+2)
         ENDDO
         X(N/2+1)=WORK(1)
         GOTO  31
23       WORK(N/2+1)=X(2)/N
         WORK(N/2+2)=-X(2)/N
         DO    I=2,N/2
            WORK(I)=X(2*I)/N
            WORK(N-I+2)=-WORK(I)
         ENDDO
         X(1)=X(1)/N
         DO    I=2,N/2+1
            X(I)=X(2*I-1)/N
         ENDDO
         DO    I=N,N/2+2,-1
            X(I)=X(N-I+2)
         ENDDO

31       WORK(1)=0.0
         CALL FFTMCF(X,WORK,N,N,N,INV)
         END
