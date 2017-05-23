#include "mex.h"
#include <math.h>
#include <string.h>

/* 
% ROWVECH computes veched matrix outproducts
% SYNTAX
%   V=rowvech(sigma,m,n);
% Given vec(sigma_i)' for i=1:N
% returns vech(sigma_i*sigma_i')' for i=1:N
% with the diagonal elements divided by 2
%
% Used by SCSolve

% Copyright (c) 2000 by Paul L. Fackler
*/


double *vech(double *A, double *B, mwSize n, bool UPPER)
{
  double *ij, *k,  *ijstart, *ijend; /*A, B, A, A*/
  mwSize  j; /*n*/

  k = B;
  if (UPPER)
  {
    if (A==B) mexErrMsgTxt("Cannot overwrite with this function");
    ijstart=A;
    ijend=A+n*n;
    for (j=0;j<n;j++)
    {
      for (ij=ijstart+j; ij<ijend; k++, ij+=n) *k=*ij; 
      ijstart += n;
    }
  }
  else 
  { 
    ij = A;
    ijend=ij+n;
    for (j=1;j<=n;j++)
    {
      for (; ij<ijend; ij++) *k++ = *ij;
      ijend=ij+n;
      ij += j;
    }
  }
  return(B);
}

double *outprod(double *A, double *B, mwSize m, mwSize n)
{
  double temp;
  mwSize ij, ik, jk, ikstart, jkstart, mm, mn;
  mn=m*n;
  mm=m*m;
  for (ikstart=0, ij=0; ij<mm; ikstart++)
    for (jkstart=0; jkstart<m; jkstart++)
    {
      ik=ikstart;
      jk=jkstart;
      temp=0;
      for(; ik<mn; ik+=m,jk+=m) temp+=A[ik]*A[jk];
      B[ij++]=temp;
    }
  return(B);
}


double *transpose(double *A, double *B, mwSize m, mwSize n)
{
  mwSize ij, ji, jistart, mn;
  mn=m*n;
  if (m==1 || n==1) memcpy(B,A,mn*sizeof(double));
  else
  {
    jistart=0;
    for (ij=0; ij<mn; jistart++)
      for (ji=jistart; ji<mn; ji+=n, ij++) B[ji]=A[ij];
  }
  return(B);
}


void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *A, *Aptr, *V, *Vptr, *temp;
   mwSize m, n, mm, mn, m1, i, N, j;
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<3 || nrhs>3)
     mexErrMsgTxt("Three arguments must be passed");
   if (nlhs>1)
     mexErrMsgTxt("ROWVECH produce only on  output");
   for (i=0; i<3; i++)
     if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))
       mexErrMsgTxt("Inputs must be full double matrices");
   if (mxGetNumberOfElements(prhs[1])!=1)
     mexErrMsgTxt("m must be a scalar");
   if (mxGetNumberOfElements(prhs[2])!=1)
     mexErrMsgTxt("n must be a scalar");
   m=*mxGetPr(prhs[1]);
   n=*mxGetPr(prhs[2]);
   mm=m*m;
   mn=m*n;
   m1=m*(m+1)/2;
   if (mxGetN(prhs[0])!=mn)
     mexErrMsgTxt("Inputs are incompatible");
   N=mxGetM(prhs[0]);
   V=mxCalloc(N*m1,sizeof(double));
   A=mxCalloc(N*mn,sizeof(double));
   temp=mxCalloc(mm,sizeof(double));
   transpose(mxGetPr(prhs[0]),A,N,mn);  

   Vptr=V;
   Aptr=A;
   for(i=0; i<N; i++, Vptr+=m1, Aptr+=mn)
   {
     outprod(Aptr,temp,m,n);
     for (j=0; j<mm; j+=(m+1)) temp[j]/=2;
     vech(temp,Vptr,m,false);
   }
   mxFree(temp);
   mxFree(A);
   plhs[0]=mxCreateDoubleMatrix(N,m1,mxREAL);
   transpose(V,mxGetPr(plhs[0]),m1,N);
   mxFree(V);
}
