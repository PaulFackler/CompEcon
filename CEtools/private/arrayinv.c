#include "mex.h"
#include <math.h>
#include "backslash.cl"

/*  
% ArrayInv Computes linear solves on arrays
% SYNTAX
%   y=ArrayInv(r,rx);
% Given r (m by p) and rx (m by p by p)
% returns y (m by p) with y(i,:)=rx(i,:,:)\r(i,:)

% Copyright (c) 2000 by Paul L. Fackler
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *r, *rx, *y, *A, *b; /*mxGetPr, mxGetPr, mxGetPr, rx, r*/
   mwSize i, j, J, Jstart; /*p, 0, 0, Jstart, 0. Dependent on the size of other mwSize variables.*/
   mwSize p, pp, m; /*mxGetN, mxGetM*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<2 || nrhs>2)
     mexErrMsgTxt("Two arguments must be passed");
   if (nlhs>1)
     mexErrMsgTxt("ArrayInv produces only one output");
   if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))       
     mexErrMsgTxt("Input arguments of inproper type");
   if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   m=mxGetM(prhs[0]);
   p=mxGetN(prhs[0]);
   if (mxGetNumberOfElements(prhs[1])!=m*p*p)
     mexErrMsgTxt("Inputs are not compatible");
   pp=p*p;

   r=mxGetPr(prhs[0]);
   rx=mxGetPr(prhs[1]);
   plhs[0]=mxCreateDoubleMatrix(m,p,mxREAL);
   y=mxGetPr(plhs[0]);

   A=mxCalloc(pp,sizeof(double));
   b=mxCalloc(p,sizeof(double));

   Jstart=0;
   for(i=0; i<m; i++, Jstart++)
   {
     for (j=0, J=Jstart; j<p; j++, J+=m) 
     {
       b[j]=r[J];
       A[j]=rx[J];
     }
     for (;j<pp; j++, J+=m) A[j]=rx[J];
     backslash(A,b,b,p,1);
     for (j=0, J=Jstart; j<p; j++, J+=m) y[J]=b[j]; 
   }
   mxFree(A);
   mxFree(b);
}
