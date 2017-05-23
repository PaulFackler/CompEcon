#include "mex.h"
#include <math.h>

/* 
% CHEBBASX A utility used by CHEBBAS

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *x, *z, *B; /*mxGetPr, B, mxGetPr*/
   double a, b, z2; /*mxGetPr, mxGetPr, z*/
   mwSize i, j, mn; /*0, 0, n*/
   long n;
   mwSize m; /*mxGetM*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<4 || nrhs>4)
     mexErrMsgTxt("Four arguments must be passed");
   if (nlhs>1)
     mexErrMsgTxt("chebbasx produces only one output");
   for (i=0; i<nrhs; i++)
     if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))       
       mexErrMsgTxt("Input arguments of inproper type");
   
   for (i=0; i<nrhs-1; i++)
     if (mxGetNumberOfElements(prhs[i])!=1)       
       mexErrMsgTxt("n, a and b must be scalars"); 

   n=(long)*mxGetPr(prhs[0]);
   if (n<1)
     mexErrMsgTxt("n must be a positive integer");
   
   a=*mxGetPr(prhs[1]);
   b=*mxGetPr(prhs[2]);
   x=mxGetPr(prhs[3]);

   m=mxGetM(prhs[3]);

   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
   B=mxGetPr(plhs[0]);

   if (n==1)
     for (i=0; i<m; i++) B[i]=1;
   else
   {
     mn=m*(n-1);
     z=B+m;
     for (i=0; i<m; i++)
     {
       B[i]=1;
       z[i] = (2/(b-a))*(x[i]-(a+b)/2);
       z2=2*z[i];
       for (j=m+i; j<mn; j+=m)
         B[j+m]=z2*B[j]-B[j-m];
     }
   }
}
