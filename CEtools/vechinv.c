#include "mex.h"
#include <math.h>
#include <string.h>

/*
Reverses the VECH operator.
*/

double *vechinv(double *A, double *B, mwSize n, int type)
{
  double *ij, *ji, *k, *ijstart, *ijend; /*B, ij, A, B, B*/
  mwSize j; /*n*/
  
  k = A;
  ij = B;
  switch (type)
  {
    case 0:                    /* symmetric */
    {
      ijend=ij+n;
      for (j=1;j<=n;j++)
      {
        ji=ij;
        for (; ij<ijend; k++, ij++, ji+=n){*ij=*k; *ji=*k;}
        ijend=ij+n;
        ij += j;
      }
      break;
    }
    case 1:                    /* upper triangular */
    {
      memset(B,0,n*n*sizeof(double));
      ijstart=B;
      ijend=B+n*n;
      for (j=0;j<n;j++)
      {
        for (ij=ijstart+j; ij<ijend; k++, ij+=n) *ij=*k; 
        ijstart += n;
      }
      break;
    }
    case 2:                    /* lower triangular */
    {
      memset(B,0,n*n*sizeof(double));
      ijend=ij+n;
      for (j=1;j<=n;j++)
      {
        for (; ij<ijend; k++, ij++) *ij=*k;
        ijend=ij+n;
        ij += j;
      }
      break;
    }
  }
  return(B);
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize nn, n; /*mxGetM, compared with nn*/
  int type;/* 0-2*/

  /* Error checking on inputs */
  if (nrhs<1) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>2) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]))
    mexErrMsgTxt("Function 'vechinv' not defined for variables of input class");
  if (mxGetNumberOfDimensions(prhs[0])>2) 
    mexErrMsgTxt("Input arguments must be 2-D.");
   
  nn=mxGetM(prhs[0])*mxGetN(prhs[0]);
  n=floor((sqrt(1+8*nn)-0.5)/2+0.5);

  if (n*(n+1)/2!=nn)
    mexErrMsgTxt("Input vector is of improper size.");

  if (nrhs>1)
  {
    type=*mxGetPr(prhs[1]);
    if (type!=1 && type!=2) type=0;
  }
  else
    type=0;

  if (mxIsComplex(prhs[0]))
  {
    plhs[0]=mxCreateDoubleMatrix(n,n,mxCOMPLEX);
    vechinv(mxGetPr(prhs[0]),mxGetPr(plhs[0]),n,type);
    vechinv(mxGetPi(prhs[0]),mxGetPi(plhs[0]),n,type);
  }
  else
  {
    plhs[0]=mxCreateDoubleMatrix(n,n,mxREAL);
    vechinv(mxGetPr(prhs[0]),mxGetPr(plhs[0]),n,type);
  }
}
