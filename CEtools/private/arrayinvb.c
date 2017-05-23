#include "mex.h"
#include <math.h>
#include "backslash.cl"

/*  
% ARRAYINVB Computes linear solves on arrays
% SYNTAX
%   y=arrayinvb(r,rx,x,xl,xu);
% Given r, x, xl, and xu (m by p) and rx (m by p by p)
% returns y (m by p) 
% The ith row of y is determined by defining:
%   b(j)=xl(i,j)-x(i,j), A(j,k)=-(j==k)   if r(i,j) <= xl(i,j)-x(i,j)
%   b(j)=r(i,j),         A(j,k)=rx(i,j,k) if xl(i,j)-x(i,j) < r(i,j) < xu(i,j)-x(i,j)
%   b(j)=xu(i,j)-x(i,j), A(j,k)=-(j==k)   if r(i,j) >= xu(i,j)-x(i,j)

% Copyright (c) 2000-2001 by Paul L. Fackler
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *r, *rx, *y, *A, *b, *x, *xl, *xu; /*mxGetPr, mxGetPr, mxGetPr, mxCalloc(double), mxCalloc(double), mxGetPr, mxGetPr, mxGetPr*/
   mwSize mp, pp, i, j, J, Jstart, jj, JJ, ppsize; /*m*p, p*p, 0, 0, Jstart, 0, j, JJ, pp. Dependent on the size of p and m.*/
   mwSize p, m; /*mxGetN, mxGetM*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<5 || nrhs>5)
     mexErrMsgTxt("Five arguments must be passed");
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
   if (mxGetNumberOfElements(prhs[2])!=m*p)
     mexErrMsgTxt("x is of improper size");
   if (mxGetNumberOfElements(prhs[3])!=m*p)
     mexErrMsgTxt("xl is of improper size");
   if (mxGetNumberOfElements(prhs[4])!=m*p)
     mexErrMsgTxt("xu is of improper size");

   mp=m*p;
   pp=p*p;
   ppsize=pp*sizeof(mwSize);

   r=mxGetPr(prhs[0]);
   rx=mxGetPr(prhs[1]);
   x=mxGetPr(prhs[2]);
   xl=mxGetPr(prhs[3]);
   xu=mxGetPr(prhs[4]);

   plhs[0]=mxCreateDoubleMatrix(m,p,mxREAL);
   y=mxGetPr(plhs[0]);
   A=mxCalloc(pp,sizeof(double));
   b=mxCalloc(p,sizeof(double));

   Jstart=0;
   for(i=0; i<m; i++, Jstart++)
   {
     memset(A,0,ppsize);
     for (j=0, J=Jstart; j<p; j++, J+=m) 
     {
       if (r[J]<=xl[J]-x[J])
       {
          b[j]=xl[J]-x[J];
          A[(p+1)*j]=-1;             /* diagonal element =-1 */
       }
       else if (r[J]>=xu[J]-x[J])
       {
          b[j]=xu[J]-x[J];
          A[(p+1)*j]=-1;            /* diagonal element =-1 */
       }
       else
       {
         b[j]=r[J];
         for (jj=j, JJ=J;jj<pp; jj+=p, JJ+=mp) A[jj]=rx[JJ];
       }
     }
     backslash(A,b,b,p,1);
     for (j=0, J=Jstart; j<p; j++, J+=m) y[J]=b[j]; 
   }
   mxFree(A);
   mxFree(b);
}
