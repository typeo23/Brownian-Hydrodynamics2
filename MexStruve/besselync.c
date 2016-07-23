
/***********************************************************************
*                                                                      *
*    Program to calculate the second kind Bessel function of integer   *
*    order N, for any REAL X, using the function BESSY(N,X).           *
*                                                                      *
* -------------------------------------------------------------------- *
*                                                                      *
*    SAMPLE RUN:                                                       *
*                                                                      *
*    (Calculate Bessel function for N=2, X=0.75).                      *
*                                                                      *
*    Second kind Bessel function of order  2 for X =  0.7500:          *
*                                                                      *
*         Y = -2.62974604                                              *
*                                                                      *
* -------------------------------------------------------------------- *
*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
*                                                                      *
*                               C++ Release 1.0 By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)             *
***********************************************************************/
#include <mex.h>
#include <math.h>
#include <stdio.h>
#include "besselync.h"

    double BESSJ0 (double X) {
/***********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
************************************************************************/
    const double
          P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
          P4=-0.2073370639E-5, P5= 0.2093887211E-6,
          Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
          Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
          R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
          R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
          S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
          S4= 59272.64853, S5=267.8532712, S6=1.0;
    double
          AX,FR,FS,Z,FP,FQ,XX,Y, TMP;

      if (X==0.0) return 1.0;
      AX = fabs(X);
      if (AX < 8.0) {
        Y = X*X;
        FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        TMP = FR/FS;
      }
      else {
        Z = 8./AX;
        Y = Z*Z;
        XX = AX-0.785398164;
        FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
      }
      return TMP;
	}

    double Sign(double X, double Y) {
      if (Y<0.0) return (-fabs(X));
      else return (fabs(X));
    }

    double BESSJ1 (double X) {
/**********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 1, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
***********************************************************************/
    const double  
      P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
      P5=-0.240337019E-6,  P6=0.636619772,
      Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
      Q4=-0.88228987E-6, Q5= 0.105787412E-6,
      R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
      R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
      S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
      S4=99447.43394,    S5=376.9991397,  S6=1.0;

	  double AX,FR,FS,Y,Z,FP,FQ,XX, TMP;

      AX = fabs(X);
      if (AX < 8.0) {
        Y = X*X;
        FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        TMP = X*(FR/FS);
      }
      else {
        Z = 8.0/AX;
        Y = Z*Z;
        XX = AX-2.35619491;
        FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        TMP = sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*Sign(S6,X);
      }
	  return TMP;
    }


    double BESSY0 (double X) {
/* --------------------------------------------------------------------
      This subroutine calculates the Second Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
  --------------------------------------------------------------------- */
      const double
		  P1= 1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
          P4=-0.2073370639E-5, P5= 0.2093887211E-6,
          Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
          Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
          R1=-2957821389.0, R2=7062834065.0, R3=-512359803.6,
          R4= 10879881.29,  R5=-86327.92757, R6=228.4622733,
          S1= 40076544269.0, S2=745249964.8, S3=7189466.438,
          S4= 47447.26470,   S5=226.1030244, S6=1.0;
      double FS,FR,Z,FP,FQ,XX,Y;
 	  if (X == 0.0) return -1e30;
      if (X < 8.0) {
        Y = X*X;
        FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        return (FR/FS+0.636619772*BESSJ0(X)*log(X));
      }
      else {
        Z = 8.0/X;
        Y = Z*Z;
        XX = X-0.785398164;
        FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        return (sqrt(0.636619772/X)*(FP*sin(XX)+Z*FQ*cos(XX)));
      }
    }

    double BESSY1 (double X) {
/* ---------------------------------------------------------------------
      This subroutine calculates the Second Kind Bessel Function of
      order 1, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
  ---------------------------------------------------------------------- */
    const double
      P1= 1.0, P2=0.183105E-2, P3=-0.3516396496E-4,
      P4= 0.2457520174E-5, P5=-0.240337019E-6,
      Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
      Q4=-0.88228987E-6, Q5= 0.105787412E-6,
      R1=-0.4900604943E13, R2= 0.1275274390E13, R3=-0.5153438139E11,
      R4= 0.7349264551E9,  R5=-0.4237922726E7,  R6= 0.8511937935E4,
      S1= 0.2499580570E14, S2= 0.4244419664E12, S3= 0.3733650367E10,
      S4= 0.2245904002E8,  S5= 0.1020426050E6,  S6= 0.3549632885E3, S7=1.0;
    double  FR,FS,Z,FP,FQ,XX, Y;
      if (X == 0.0) return -1e30;
      if (X < 8.0) {
        Y = X*X;
        FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7)))));
        return (X*(FR/FS)+0.636619772*(BESSJ1(X)*log(X)-1.0/X));
      }
      else {
         Z = 8./X;
        Y = Z*Z;
        XX = X-2.356194491;
        FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        return (sqrt(0.636619772/X)*(sin(XX)*FP+Z*cos(XX)*FQ));
      }
    }

	void BESSY (int N, double X, double *BESSYA) {
/* -----------------------------------------------------------------
      This subroutine calculates the second kind Bessel Function of
      integer order N, for any real X. We use here the classical
      recursive formula. 
  ------------------------------------------------------------------ */
    double TOX,BY,BYM,BYP; int J;
      if (N == 0){ *BESSYA= BESSY0(X); return;}
      if (N == 1)  {*BESSYA= BESSY1(X);return;}
      if (X == 0.0) { *BESSYA= 1e30;return;}
      TOX = 2.0/X;
      BY  = BESSY1(X);
      BYM = BESSY0(X);
      for (J=1; J<N; J++) {
        BYP = J*TOX*BY-BYM;
        BYM = BY;
        BY  = BYP;
      };
      *BESSYA= BY;
      return;
	}