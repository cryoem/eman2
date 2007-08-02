/******************************************************
* Program to demonstrate the use of multi-dimensional *
*     Steepest Descent Optimization subroutine        *
* --------------------------------------------------- *
*  Reference: BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*                  C++ version by J-P Moreau, Paris.  *
* --------------------------------------------------- *
* Example:   Find a local maximum of                  *
*            F(x,y,z) = sin(x)+2*cos(y)-sin(z)        *
*                                                     *
* SAMPLE RUN:                                         *
*                                                     *
* How many dimensions: 3                              *
*                                                     *
* Convergence criterion: .000000001                   *
*                                                     *
* Maximum number of iterations: 50                    *
*                                                     *
* Starting constant: 1                                *
*                                                     *
* Input the starting point:                           *
*                                                     *
*     X[1] = 1                                        *
*     X[2] = 1                                        *
*     X[3] = 1                                        *
*                                                     *
* The results are:                                    *
*                                                     *
*     X(1) = 1.5707973                                *
*     X(2) = -0.0000170                               *
*     X(3) = -4.7123856                               *
*                                                     *
* Local maximum = 4.0000000                           *
*                                                     *
* The number of steps was: 33                         *
*                                                     *
******************************************************/


#include <stdio.h>
#include <math.h>

#define  MACHEPS  1e-15


/*******************************************
  Function subroutine                     */   
double Eval(double *X) {
  return sin(X[1])+2.0*cos(X[2])-sin(X[3]);
}
/*******************************************/

  // Functions called by Steepda()
  void Utilit1(double *D, double *dd, int l) {
    int i;
    // Find the magnitude of the gradient
    *dd=0.0;
    for (i=1; i<l+1; i++) *dd += D[i]*D[i];
    *dd=sqrt(*dd);
  }

  void Utilit2(double *X, double *X1, double *Y, double *D, double *dd, double xk, int l) {
	int i;
    // Update the X[i] 
    for (i=1; i<l+1; i++) {
      // Save old values
      X1[i]=X[i];
      X[i] += xk*D[i]/(*dd);
    }
    Y[3]=Eval(X);
  }

  // Find approximations of partial derivatives D(i)
  // by finite differences
  void Derivatives(double *X, double *D, double *Y, double *dd, double xk, int l)  {
    double a,b,yy;
    int i;
    for (i=1; i<l+1; i++) {
      // Save X(i)
      a=X[i];
      // Find increment
      b=D[i]*xk/(2.0*(*dd));
      // Move increment in X(i)
      X[i]=X[i]+b;
      // Obtain yy
      yy=Eval(X);
      // Guard against divide by zero near maximum
      if (b==0) b=1e-12;
      // Update D(i)
      D[i]=(yy-Y[3])/b;
      // Guard against locked up derivative
      if (D[i]==0) D[i]=1e-5;
      // Restore X(i) and yy
      X[i]=a; yy=Y[3];
    }
    // Obtain dd
    Utilit1(D, dd, l);
  }


/**************************************************
*    Steepest descent optimization subroutine     *
* ----------------------------------------------- *
* This routine finds the local maximum or minimum *
* of an L-dimensional function using the method   *
* of steepest decent, or the gradient.            *
* The function must be available in subroutine    *
* Eval(). In this version, finite differences are *
* used to calculate the L partial derivatives.    *
* ----------------------------------------------- *
* INPUTS:                                         *
*   l - The dimension of function to study        *
*   e - The convergence criteria                  *
*   m - The maximum number of iterations          *
*   xk - A starting constant                      *
*   X(i) - Initial values of variables            *
* OUTPUTS:                                        *
*   X(i) - The locally optimum set                *
*   Eval - The local maximum found                *
*   n - The number of iterations performed,       *
**************************************************/
  void Steepda(double *D, double *Y, double *X, double *X1, double dd, double e, double xk, int l, int m, int *n)  {
  // Labels: e50,e51,e100,e200
  int i;
  *n=0;
  //The routine needs three values of Y to get started
  //Generate starting D(i) values
  //These are not even good guesses and slow the program a little
  dd=1.0;
  D[1]=1.0/sqrt(l);
  for (i=2; i<l+1; i++)  D[i]=D[i-1];
  // Start initial probe
  for (i=1; i<l+1; i++) {
    // Obtain yy and D[i]
    Y[i]=Eval(X);
    // Update X[i]
    Utilit1(D, &dd, l);
    Utilit2(X, X1, Y, D, &dd, xk, l);
  }
  // We now have a history to base the subsequent search on
  // Accelerate search if approach is monotonic 
e50: if (fabs(Y[2]-Y[1])<MACHEPS) goto e51;
  if ((Y[3]-Y[2])/(Y[2]-Y[1])>0.0) xk=xk*1.2;
  // Decelerate if heading the wrong way
e51: if (Y[3]<Y[2]) xk=xk/2.0;
  // Update the Y[i] if value has decreased
  if (Y[3]>Y[2]) goto e100;
  // Restore the X[i]
  for (i=1; i<l+1; i++) {
    X[i]=X1[i];
  }
  goto e200;
e100: Y[1]=Y[2]; Y[2]=Y[3];
  // Obtain new values
e200: Y[3]=Eval(X);
  Derivatives(X, D, Y, &dd, xk, l); // Get D(i)
  //if dd=0 then the precision limit of the computer has been reached
  if (dd==0) return;
  // Update X[i]
  Utilit2(X, X1, Y, D, &dd, xk, l);
  // Check for maximum iterations and convergence
  (*n)++;
  if (*n>=m) return;
  if (fabs(Y[3]-Y[2])<e) return;
  // Try another iteration
  goto e50;
} // Steepds()

/*
int main() {

	double  D[4], Y[4];
	double  X[11],X1[11];
	double  dd,e,xk;
	int i,l,m,n;
	
  printf("\n How many dimensions: "); scanf("%d",&l);
  printf("\n Convergence criterion: "); scanf("%lf",&e);
  printf("\n Maximum number of iterations: "); scanf("%d",&m);
  printf("\n Starting constant: "); scanf("%lf",&xk);
  printf("\n Input the starting point:\n\n");
  for (i=1; i<l+1; i++) {
    printf("     X(%d) = ",i); scanf("%lf",&X[i]);
  }

  Steepda(D, Y, X, X1, dd, e, xk, l, m, &n);   // Call steepest descent optimization subroutine

  printf("\n\n The results are:\n\n");
  for (i=1; i<l+1; i++) {
    printf("     X(%d) = %1.7f\n",i,X[i]);   
  }
  printf("\n Local maximum found = %10.7f\n",Eval(X));
  printf("\n The number of iterations was %d\n\n",n);
}*/

// End of file Steepda.cpp
