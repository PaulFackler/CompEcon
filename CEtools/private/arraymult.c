#include "mex.h"
#include <math.h>

/* 
% ArrayMult Computes matrix products over 3-D arrays
% SYNTAX:
%   c=arraymult(a,b,n,p,q,r);
% Inputs:
%   a : n x p x q array
%   b : n x q x r array
%   n,p,q,r : scalar  
% Output:
%   c : n x p x r array
%
%   c(i,:,:) = a(i,:,:)*b(i,:,:), i=1,...,n

% Copyright (c) 2000 by Paul L. Fackler
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *A, *B, *C, *Aptr, *Bptr, *AA; /*mxGetPr, mxGetPr, mxGetPr, AA, B, A*/
   mwSize i, j, k, l, np, nq; /*0, 0, 0, 0, n, n*/
   mwSize n, p, q, r, *dims; /*All these need to be mwSize to work with mxCreateNumericArray.*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<6 || nrhs>6)
     mexErrMsgTxt("Six arguments must be passed");
   if (nlhs>1)
     mexErrMsgTxt("ArrayMult produces only one output");
   for (i=0; i<nrhs; i++)
     if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))       
       mexErrMsgTxt("Input arguments of inproper type");
   
   for (i=2; i<nrhs; i++)
     if (mxGetNumberOfElements(prhs[i])!=1)       
       mexErrMsgTxt("n, p, q, and r must be scalars"); 

   /*Note that below we need to recast n,p,q,r to mwSize format.*/
   n=(mwSize)*mxGetPr(prhs[2]);
   p=(mwSize)*mxGetPr(prhs[3]);
   q=(mwSize)*mxGetPr(prhs[4]);
   r=(mwSize)*mxGetPr(prhs[5]);

   if (mxGetNumberOfElements(prhs[0])!=n*p*q)
      mexErrMsgTxt("A is of improper size"); 
   if (mxGetNumberOfElements(prhs[1])!=n*q*r)
      mexErrMsgTxt("B is of improper size"); 

   A=mxGetPr(prhs[0]);
   B=mxGetPr(prhs[1]);

   dims=mxCalloc(3,sizeof(mwSize)); /*Edited to "mwSize".*/
   dims[0]=n; dims[1]=p; dims[2]=r;
   plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
   C=mxGetPr(plhs[0]);
   mxFree(dims);

  np=n*p;
  nq=n*q;

  for (k=0; k<r; k++)
  {
    AA=A;
    for (j=0; j<p; j++)
    {
      Aptr=AA;
      Bptr=B;
      for (l=0;l<q;l++)
      {
        for(i=0; i<n; i++)
          C[i] += Aptr[i]*Bptr[i];
        Aptr=Aptr+np;
        Bptr=Bptr+n;
      }
      C+=n;
      AA+=n;
    }
    B+=nq;
  }  
}
