#include "mex.h"
#include <math.h>


/*
% DIGAMMA Computes the digamma (psi) function for positive arguments
% USAGE
%   y=digmma(x)
*/

double psi(double x);

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize n;
  double *p, *pend, *x; /*mxGetPr, p, mxGetPr*/

  /* Error checking on inputs */  
  if (nrhs<1) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>1) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  if (!mxIsDouble(prhs[0]) && !mxIsSparse(prhs[0]))
      mexErrMsgTxt("Function not defined for variables of input class");
  if (mxIsComplex(prhs[0]))
      mexErrMsgTxt("X must be real.");

  x=mxGetPr(prhs[0]);
  n=mxGetNumberOfElements(prhs[0]);
  plhs[0]=mxDuplicateArray(prhs[0]);
  p=mxGetPr(plhs[0]);
  for (pend=p+n; p<pend; p++,x++) { *p=psi(*x);} 

}


double psi(double x)
{
  double p;
  int inc;
  if (x<=0) {p=mxGetNaN(); return(p);}
  if (x<10) {inc=10-floor(x); x+=inc;}
  else       inc=0;
  p=1.0/(x*x);
  p=((((((-8.33333333333333333333e-2*p  \
          +2.10927960927960927961e-2)*p \
          -7.57575757575757575758e-3)*p \
          +4.16666666666666666667e-3)*p \
          -3.96825396825396825397e-3)*p \
          +8.33333333333333333333e-3)*p \
          -8.33333333333333333333e-2)*p;
  p+=log(x)-0.5/x;
  for (;inc>0; inc--) p-=1/(x-inc);
  return(p);
}
