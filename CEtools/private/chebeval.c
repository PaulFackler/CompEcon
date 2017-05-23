#include "mex.h"
#include <math.h>
#include <string.h>

/*
% CHEBEVAL Evaluates Chebyshev polynomial approximations
% USAGE
%   y=chebeval(c,x,n,a,b);
% INPUTS
%   c : N x p matrix
%   x : m x d matrix
%   n : d vector of integers with prod(n)=N
%   a : d vector of lower bounds
%   b : d vector of upper bounds
% OUTPUT
%   y : m x p matrix
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  double *a, *b, *x, *c, *y; /*mxGetPr, mxGetPr, mxGetPr, mCalloc(double), mxGetPr*/
  double *cptr, *ccol, *xptr, *yptr, *cin, *cend, *cp, *cstart; /*ccol, cend, x, cptr, mxGetPr, c, ccol, cend*/
  double z, z2, dd, di, temp; /*2.0, 2.0, 0.0, 0.0, di*/
  mwSize *cind, *n; /*p, (mwSize)*/
  mwSize d, m, p, N, i, j, nj, xrows, csize; /**/

  /* Error checking on inputs */  
  if (nrhs<5) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>5) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  for (i=0; i<nrhs; i++) {
    if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("Arguments must be real.");
  }

  p=mxGetN(prhs[0]);

  m=mxGetM(prhs[1]);
  d=mxGetN(prhs[1]);
  x=mxGetPr(prhs[1]);

  if (mxGetNumberOfElements(prhs[2])!=d)
    mexErrMsgTxt("n has improper dimension");
  a=mxGetPr(prhs[2]);
  n=mxCalloc(d,sizeof(mwSize));
  for (i=0;i<d;i++) n[i]=(mwSize)a[i];

  if (mxGetNumberOfElements(prhs[3])!=d)
    mexErrMsgTxt("a has improper dimension");
  a=mxGetPr(prhs[3]);

  if (mxGetNumberOfElements(prhs[4])!=d)
    mexErrMsgTxt("b has improper dimension");
  b=mxGetPr(prhs[4]);

  N=1; for(i=0;i<d; i++) N *= n[i];
  if (mxGetM(prhs[0])!=N)
    mexErrMsgTxt("c has improper number of rows");

  plhs[0]=mxCreateDoubleMatrix(m,p,mxREAL);

  cind=mxCalloc(d,sizeof(mwSize));  
  cind[d-1]=p; for (i=d-1;i>0;i--) cind[i-1]=cind[i]*n[i];
  c=mxCalloc(N*p,sizeof(double));
  y=mxGetPr(plhs[0]);

  cin=mxGetPr(prhs[0]);
  csize=N*p*sizeof(double);
  cend=c+N*p;
  for (xrows=0;xrows<m; xrows++)
  {  
    memcpy(c,cin,csize);
    xptr=x+xrows;
    for (j=0; j<d; j++, xptr+=m)
    {
      z=(2.0**xptr-a[j]-b[j])/(b[j]-a[j]);
      z2=z*2.0;
      nj=n[j];
      ccol=cend-1;
      cptr=ccol;
      cstart=cend;
      for (cp=ccol-cind[j]; ccol>cp;)
      {
        di=0.0; 
        dd=0.0;
        for (cstart-=nj;cptr>cstart;)
        {
          temp=di;
          di=z2*di-dd+*cptr--;
          dd=temp;
        }
        *ccol--=z*di-dd+*cptr--;
      }
    }
    for (cptr=cend-p, yptr=y+xrows; cptr<cend; yptr+=m) 
       *yptr=*cptr++;
  }
  mxFree(c);
  mxFree(cind);
  mxFree(n);
}
