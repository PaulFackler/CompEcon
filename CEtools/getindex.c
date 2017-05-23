#include "mex.h"
#include <math.h>
/*
% GETINDEX Finds the index value of a point
% USAGE
%   i = getindex(s,S);
% INPUTS
%   s : a pxd matrix
%   S : an nxd matrix
% OUTPUT
%   i : a p-vector of integers in {1,...,n} indicating the row of
%         S that most closely matches each row in s

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{  /*Note, m was removed since it is not used.*/
  mwSize n, p, d; /*mxGetM, mxGetM, mxGetN*/
  mwSize i, j, k, minind; /*0, 0, 0, 1*/
  double *s, *S, *ind, minds, dsj, *sk, *Sk, inf; /*mxGetPr, mxGetPr, mxGetPr, inf, minds, s, S, mxGetInf*/

  /* Error checking on inputs */  
  if (nrhs<2) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>2) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>1) mexErrMsgTxt("Too many output arguments.");
  for (i=0; i<nrhs; i++)
  {
    if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("Arguments must be real.");
  }

  d=mxGetN(prhs[0]);
  if (mxGetN(prhs[1])!=d) 
    mexErrMsgTxt("Inputs must have equal numbers of columns.");
  p=mxGetM(prhs[0]);
  n=mxGetM(prhs[1]);
 
  plhs[0]=mxCreateDoubleMatrix(p,1,mxREAL);
 
  s=mxGetPr(prhs[0]);
  S=mxGetPr(prhs[1]);
  ind=mxGetPr(plhs[0]);
  inf=mxGetInf();
  
  for (i=0; i<p; i++, s++){
    minind=1;
    minds=inf;
    for (j=0; j<n;){
      for (k=0, dsj=0, sk=s, Sk=S+j; k<d; k++, sk+=p, Sk+=n) 
        dsj+=fabs(*sk-*Sk);
      j++;
      if (dsj==0){minind=j; break;}
      if (dsj<minds){minds=dsj; minind=j;}
    }
    ind[i]=minind;
  }
}
