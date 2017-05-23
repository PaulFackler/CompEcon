#include "mex.h"
#include <math.h>

/*
Puts an lower triangular or symmetric matrix into vech form: e.g.
     1  
     2 5 
     3 6 8
     4 7 9 10
*/

double *vech(double *A, double *B, mwSize n, bool UPPER)
{
  double *ij, *k,  *ijstart, *ijend; /*ijstart, B, A, A*/
  mwSize  j;/*0.  Depending on the size of n.*/

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


void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize n, n1; /*mxGetM*/
  bool type; /*true*/
   
  if (nrhs<1) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>2) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]))    
    mexErrMsgTxt("Function 'vech' not defined for variables of input class");
  n=mxGetM(prhs[0]);
  if (mxGetN(prhs[0])!=n) mexErrMsgTxt("Input matrix must be square");

  n1=n*(n+1)/2;

  if (nrhs>1)
    if (*mxGetPr(prhs[1])==1) type=true; else type=false;
  else
    type=false;

  if (mxIsComplex(prhs[0]))
  {
    plhs[0]=mxCreateDoubleMatrix(n1,1,mxCOMPLEX);
    vech(mxGetPr(prhs[0]),mxGetPr(plhs[0]),n,type);
    vech(mxGetPi(prhs[0]),mxGetPi(plhs[0]),n,type);
  }
  else
  {
    plhs[0]=mxCreateDoubleMatrix(n1,1,mxREAL);
    vech(mxGetPr(prhs[0]),mxGetPr(plhs[0]),n,type);
  }
}
