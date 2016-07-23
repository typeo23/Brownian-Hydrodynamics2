/**************************************************************
!*       Purpose: This program computes Struve function       *
!*                H1(x) using subroutine STVH1                *
!*       Input :  x   --- Argument of H1(x) ( x � 0 )         *
!*       Output:  SH1 --- H1(x)                               *
!*       Example:                                             *
!*                   x          H1(x)                         *
!*                -----------------------                     *
!*                  0.0       .00000000                       *
!*                  5.0       .80781195                       *
!*                 10.0       .89183249                       *
!*                 15.0       .66048730                       *
!*                 20.0       .47268818                       *
!*                 25.0       .53880362                       *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                          C++ Release By J-P Moreau, Paris. *
!*                                  (www.jpmoreau.fr)         *
!*************************************************************/ 
#include <stdio.h>
#include <math.h>
#ifdef MATLAB_MEX_FILE
#include <tmwtypes.h>
#else
#include "rtwtypes.h"
#endif



void STVH1(double X, double *SH1) {
/*      =============================================
!       Purpose: Compute Struve function H1(x)
!       Input :  x   --- Argument of H1(x) ( x � 0 )
!       Output:  SH1 --- H1(x)
!       ============================================= */
        double A0,BY1,P1,PI,Q1,R,S,T,T2,TA1;
	int K, KM;

        PI=3.141592653589793;
        R=1.0;
        if (X <= 20.0) {
           S=0.0;
           A0=-2.0/PI;
           for (K=1; K<=60; K++) {
              R=-R*X*X/(4.0*K*K-1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
e15:       *SH1=A0*S;
        }
        else {
           S=1.0;
           KM=(int)(0.5*X);
           if (X > 50.0)  KM=25;
           for (K=1; K<=KM; K++) {
              R=-R*(4.0*K*K-1.0)/(X*X);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
e25:       T=4.0/X;
           T2=T*T;
           P1=((((0.42414e-5*T2-0.20092e-4)*T2+0.580759e-4)*T2-0.223203e-3)*T2+0.29218256e-2)*T2+0.3989422819;
           Q1=T*(((((-0.36594e-5*T2+0.1622e-4)*T2-0.398708e-4)*T2+0.1064741e-3)*T2-0.63904e-3)*T2+0.0374008364);
           TA1=X-0.75*PI;
           BY1=2.0/sqrt(X)*(P1*sin(TA1)+Q1*cos(TA1));
           *SH1=2.0/PI*(1.0+S/(X*X))+BY1;
        }
        return ;
}




/* end of file mstvh1.cpp */