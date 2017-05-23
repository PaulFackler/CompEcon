#include "mex.h"
#include <math.h>
/*
% DIAGMULT Computes either diag(a)*b or a*diag(b)
% USAGE
%   c=diagmult(a,b);
% INPUTS
%   a,b : m-vector and mxn matrix
%         or 
%         mxn matrix and n-vector
% OUTPUT
%   c   : mxn matrix 
%
% Note: not implemented for complex matrices or matrices 
% with data type other than double. The vector must be full but
% the matrix can be sparse or full.

% Copyright (c) 2005, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu
% 
*/


/* diag(a) times B - full B */
void daxbf(double *a, double *B, mwSize m, mwSize n)
{
  double *aend, *Bend; /*a, B*/
  aend=a+m;
  Bend=B+m*n;
  while (B<Bend){
    while (a<aend) *B++ *= *a++;
    a-=m;
  }
}

/* diag(a) times B - sparse B */
void daxbs(double *a, double *B, mwIndex *Bi, mwIndex *Bj, mwSize m, mwSize n)
{
  mwIndex i, iend; /*Bj*/
  iend=Bj[n];
  for (i=0; i<iend;i++) *B++ *= a[*Bi++];
}

/* A times diag(b) - full A */
void axdbf(double *A, double *b, mwSize m, mwSize n)
{
  mwSize i;/*0*/
  double *Aend, bval; /*A, b*/
  Aend=A+m*n;
  for (; A<Aend;b++)
    for (i=0, bval=*b; i<m; i++) *A++ *= bval;
}

/* A times diag(b) - sparse A */
void axdbs(double *A, double *b, mwIndex *Aj, mwSize m, mwSize n)
{
  mwSize i, inext, j; /*inext, 0, 0*/
  double bval; /*b*/
  inext=0;
  for (j=0; j<n; j++){
    bval=*b++;
    for (i=inext, inext=*(++Aj); i<inext; i++) *A++ *= bval;
  }
}



void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{  double *A, *B, a; /*mxGetPr, mxGetPr, mxGetPr*/ 
   mwSize i, mb, nb, m, n;  /*m, mxGetM, mxGetN , mxGetM, mxGetN*/

   if (nrhs!=2)
      mexErrMsgTxt("Two parameters must be passed");
   if (nlhs>1)
      mexErrMsgTxt("Only one output is created");
   if (!mxIsDouble(prhs[0]) &  !mxIsSparse(prhs[0]))
      mexErrMsgTxt("Inputs must be matrices"); 
   if (!mxIsDouble(prhs[1]) &  !mxIsSparse(prhs[1]))
      mexErrMsgTxt("Inputs must be matrices");
   if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
      mexErrMsgTxt("Inputs must be real");

   m=mxGetM(prhs[0]); 
   n=mxGetN(prhs[0]);

   if (m==1 || n==1){
     if (n>1) m=n;
     /* A is a scalar */
     if (m==1){
       a=*mxGetPr(prhs[0]);
       plhs[0]=mxDuplicateArray(prhs[1]);
       B=mxGetPr(plhs[0]);
       if (mxIsSparse(plhs[0])){
         mwIndex *Bi, *Bj;  /*mxGetJc*/
         n=mxGetN(prhs[1]);
         Bj=mxGetJc(plhs[0]);
         m=Bj[n];
       }
       else
         m=mxGetM(prhs[1])*mxGetN(prhs[1]);
       for (i=0; i<m; i++) B[i]*=a;
         
     }
     /* A is a vector */
     else{
       if (mxGetM(prhs[1])!=m) mexErrMsgTxt("Inputs are not comformable");
       n=mxGetN(prhs[1]);
       A=mxGetPr(prhs[0]);
       plhs[0]=mxDuplicateArray(prhs[1]);
       B=mxGetPr(plhs[0]);
       if (mxIsSparse(plhs[0])){
         mwIndex *Bi, *Bj;  /*mxGetIr, mxGetJc*/
         Bi=mxGetIr(plhs[0]);
         Bj=mxGetJc(plhs[0]);
         daxbs(A,B,Bi,Bj,m,n);
       }
     else
       daxbf(A,B,m,n);
     }
   }
   /* A is a matrix */
   else{
     mb=mxGetM(prhs[1]);
     nb=mxGetN(prhs[1]);
     if (mb>1 && nb>1) mexErrMsgTxt("One of the inputs must be a vector");
     if (nb==1) nb=mb;
     if (nb!=n) mexErrMsgTxt("Inputs are not comformable");
     B=mxGetPr(prhs[1]);
     plhs[0]=mxDuplicateArray(prhs[0]);
     A=mxGetPr(plhs[0]);
     if (mxIsSparse(plhs[0])){
       mwIndex *Aj;
       Aj=mxGetJc(plhs[0]);
       axdbs(A,B,Aj,m,n);
     }
     else
       axdbf(A,B,m,n);
   }
}
