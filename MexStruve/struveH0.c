
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


void STVH0(double Xin[], double *SH0) {
/*      =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x � 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= */
        double A0,BY0,P0,PI,Q0,R,S,T,T2,TA0,X;
	int K, KM;
	X = Xin[0];    
	PI=3.141592653589793;
        S=1.0;
        R=1.0;
        if (X <= 20.0) {
           A0=2.0*X/PI;
           for (K=1; K<61; K++) {
              R=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
e15:       *SH0=A0*S;
        }
        else {
           KM=(int)(0.5*(X+1.0));
           if (X >= 50.0) KM=25;
           for (K=1; K<=KM; K++) {
              R=-R*pow((2.0*K-1.0)/X,2);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
e25:       T=4.0/X;
           T2=T*T;
           P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
           Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
           TA0=X-0.25*PI;
           BY0=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
           *SH0=2.0/(PI*X)*S+BY0;
        }
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
  STVH0(x,y);
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

// end of file mstvh0.cpp */