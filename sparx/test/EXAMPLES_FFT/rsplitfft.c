
Here is a C function for doing a real fft (split radix).  It is a translation
and correction of some Fortran code.  I have tested it somewhat, but maybe
you could test it yourself (and see how it compares to others).

Below is the calling program, followed by rfft().

Currently it uses Fortran style: x[1] is first element not x[0].  I couldn't
figure out a simple way to convert to C style.  I guess it's not a 
catastrophe to waste one element (x[0]).

Cheers,
Bill Simpson

#include <stdio.h>
#include <math.h> 

void rfft(float X[],int N); 

int main(void)
{
float x[2048];
int i,n,M;

n=8;
M=(int)(log(n)/log(2.0));
printf("%d=2^%d\n",n,M);

for(i=1;i<=n;i++)
        x[i]=1;
        /*x[i]=(float)10.0*cos(6.2831853071719586lf*3.0/n*i);*/

rfft(x,n);
for(i=1;i<=n;i++)
        printf("%d, %f\t",i,x[i]);

return 0;
}




/****************************************************************************
* rfft(float X[],int N)                                                     *
*     A real-valued, in-place, split-radix FFT program                      *
*     Decimation-in-time, cos/sin in second loop                            *
*     Input: float X[1]...X[N] (NB Fortran style: 1st pt X[1] not X[0]!)    *
*     Length is N=2**M (i.e. N must be power of 2--no error checking)       *
*     Output in X[1]...X[N], in order:                                      *
*           [Re(0), Re(1),..., Re(N/2), Im(N/2-1),..., Im(1)]               *
*                                                                           *
* Original Fortran code by Sorensen; published in H.V. Sorensen, D.L. Jones,*
* M.T. Heideman, C.S. Burrus (1987) Real-valued fast fourier transform      *
* algorithms.  IEEE Trans on Acoustics, Speech, & Signal Processing, 35,    *
* 849-863.  Adapted to C by Bill Simpson, 1995  wsimpson@uwinnipeg.ca       *
****************************************************************************/

void rfft(float X[],int N)
{
int I,I0,I1,I2,I3,I4,I5,I6,I7,I8, IS,ID;
int J,K,M,N2,N4,N8;
float A,A3,CC1,SS1,CC3,SS3,E,R1,XT;
float T1,T2,T3,T4,T5,T6;  

M=(int)(log(N)/log(2.0));               /* N=2^M */

/* ----Digit reverse counter--------------------------------------------- */
J = 1;
for(I=1;I<N;I++)
        {
        if (I<J)
                {
                XT    = X[J];
                X[J]  = X[I];
                X[I]  = XT;
                }
        K = N/2;
        while(K<J)
                {
                J -= K;
                K /= 2;
                }
        J += K;
        }

/* ----Length two butterflies--------------------------------------------- */
IS = 1;
ID = 4;
do
        {
        for(I0 = IS;I0<=N;I0+=ID)
                {
                I1    = I0 + 1;
                R1    = X[I0];
                X[I0] = R1 + X[I1];
                X[I1] = R1 - X[I1];
                }
        IS = 2 * ID - 1;
        ID = 4 * ID;
        }while(IS<N);
/* ----L shaped butterflies----------------------------------------------- */
N2 = 2;
for(K=2;K<=M;K++)
        {
        N2    = N2 * 2;
        N4    = N2/4;
        N8    = N2/8;
        E     = (float) 6.2831853071719586f/N2;
        IS    = 0;
        ID    = N2 * 2;
        do
                {
                for(I=IS;I<N;I+=ID)
                        {
                        I1 = I + 1;
                        I2 = I1 + N4;
                        I3 = I2 + N4;
                        I4 = I3 + N4;
                        T1 = X[I4] +X[I3];
                        X[I4] = X[I4] - X[I3];
                        X[I3] = X[I1] - T1;
                        X[I1] = X[I1] + T1;
                        if(N4!=1)
                                {
                                I1 += N8;
                                I2 += N8;
                                I3 += N8;
                                I4 += N8;
                                T1 = (X[I3] + X[I4])*.7071067811865475244f;
                                T2 = (X[I3] - X[I4])*.7071067811865475244f;
                                X[I4] = X[I2] - T1;
                                X[I3] = -X[I2] - T1;
                                X[I2] = X[I1] - T2;
                                X[I1] = X[I1] + T2;
                                }
                        }
                        IS = 2 * ID - N2;
                        ID = 4 * ID;
                }while(IS<N);
        A = E;
        for(J= 2;J<=N8;J++)
                {
                A3 = 3.0 * A;
                CC1   = cos(A);
                SS1   = sin(A);  /*typo A3--really A?*/
                CC3   = cos(A3); /*typo 3--really A3?*/
                SS3   = sin(A3);
                A = (float)J * E;
                IS = 0;
                ID = 2 * N2;
                do        
                        {
                        for(I=IS;I<N;I+=ID)
                                {
                                I1 = I + J;
                                I2 = I1 + N4;
                                I3 = I2 + N4;
                                I4 = I3 + N4;
                                I5 = I + N4 - J + 2;
                                I6 = I5 + N4;
                                I7 = I6 + N4;
                                I8 = I7 + N4;
                                T1 = X[I3] * CC1 + X[I7] * SS1;
                                T2 = X[I7] * CC1 - X[I3] * SS1;
                                T3 = X[I4] * CC3 + X[I8] * SS3;
                                T4 = X[I8] * CC3 - X[I4] * SS3;
                                T5 = T1 + T3;
                                T6 = T2 + T4;
                                T3 = T1 - T3;
                                T4 = T2 - T4;
                                T2 = X[I6] + T6;
                                X[I3] = T6 - X[I6];
                                X[I8] = T2;
                                T2    = X[I2] - T3;
                                X[I7] = -X[I2] - T3;
                                X[I4] = T2;
                                T1    = X[I1] + T5;
                                X[I6] = X[I1] - T5;
                                X[I1] = T1;
                                T1    = X[I5] + T4;
                                X[I5] = X[I5] - T4;
                                X[I2] = T1;
                                }
                        IS = 2 * ID - N2;
                        ID = 4 * ID;
                        }while(IS<N);
                }
        }
return;
}


