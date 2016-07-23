
/**************************************************************
!*       Purpose: This program computes Struve function       * 
!*                H0(x) using subroutine STVH0                *
!*       Input :  x   --- Argument of H0(x) ( x � 0 )         *
!*       Output:  SH0 --- H0(x)                               *
!*       Example:                                             *
!*                   x          H0(x)                         *
!*                ----------------------                      *
!*                  0.0       .00000000                       * 
!*                  5.0      -.18521682                       *
!*                 10.0       .11874368                       *
!*                 15.0       .24772383                       *
!*                 20.0       .09439370                       *
!*                 25.0      -.10182519                       *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                          C++ Release By J-P Moreau, Paris. *
!*                                  (www.jpmoreau.fr)         *
!*************************************************************/ 
#include <mex.h>
#include <stdio.h>
#include <math.h>


static double bessy0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          0 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
         +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
      ans2=40076544269.0+y*(745249964.8+y*(7189466.438
         +y*(47447.26470+y*(226.1030244+y*1.0))));
      ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         +y*(-0.934945152e-7))));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}



static double bessy1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          1 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1=x*(-0.4900604943e13+y*(0.1275274390e13
         +y*(-0.5153438139e11+y*(0.7349264551e9
         +y*(-0.4237922726e7+y*0.8511937935e4)))));
      ans2=0.2499580570e14+y*(0.4244419664e12
         +y*(0.3733650367e10+y*(0.2245904002e8
         +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
      ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}



/*
#>            bessy.dc2

Function:     bessy

Purpose:      Evaluate Bessel function second kind and of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessy( int n,
                              double x )


              bessy    Return the Bessel function of second kind and
                       of integer order, for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessy evaluates at x the Bessel function of the second kind
              and of integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/


void bessy( int n, double xin[], double y[] )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          n for input x. (n >= 0)                           */
/* Note that for x == 0 the functions bessy and bessk are not */
/* defined and a blank is returned.                           */
/*------------------------------------------------------------*/
{
   int j;
   double by,bym,byp,tox,x;
   x = Xin[0]; 

   if (n < 0 || x == 0.0)
   {
      double   dblank;
      setdblank_c( &dblank );
      return( dblank );
   }
   if (n == 0)
      return( bessy0(x) );
   if (n == 1)
      return( bessy1(x) );

   tox=2.0/x;
   by=bessy1(x);
   bym=bessy0(x);
   for (j=1;j<n;j++) {
      byp=j*tox*by-bym;
      bym=by;
      by=byp;
   }
   y[0] = by;
   return ;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x,*y;
  size_t mrows,ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=1) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:invalidNumInputs",
            "One input required.");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:maxlhs",
            "Too many output arguments.");
  }
  
  /* The input must be a noncomplex scalar double.*/
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Input must be a noncomplex scalar double.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
  
  /* Assign pointers to each input and output. */
  x = mxGetPr(prhs[0]);
  y = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */
  STVH1(x,y);
}

/*
void main() {
        
        double X, SH0;

	printf("\n Please enter x: ");
        scanf("%lf", &X);
	printf("\n");
        printf("   x        H0(x)     \n");
        printf("----------------------\n");
        
	STVH0(X, &SH0);
        
	printf("%5.1f  %e\n\n", X, SH0);

}

/* end of file mstvh0.cpp */