#include "mex.h"
#include <math.h>
#include "backslash.cl"

/*  
% ARRAYMM Minmax function reformulation & Newton step on arrays
% SYNTAX
%   [fx,Jac]=arrayinvb(x,xl,xu,fx,Jac);
% Given x, xl, xu and fx (m by p) and Jac (m by p by p)
% returns reformulated fx and the Newton step.
% The ith row of J is determined by defining:
%   b(j)=xl(i,j)-x(i,j), A(j,k)=-(j==k)    if fx(i,j) <= xl(i,j)-x(i,j)
%   b(j)=fx(i,j),        A(j,k)=Jac(i,j,k) if xl(i,j)-x(i,j) < fx(i,j) < xu(i,j)-x(i,j)
%   b(j)=xu(i,j)-x(i,j), A(j,k)=-(j==k)    if fx(i,j) >= xu(i,j)-x(i,j)
% The ith Newton step is -A\b

% Copyright (c) 2000-2001 by Paul L. Fackler
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   double *fx, *Jac, *A, *b, *x, *xl, *xu, *Newton, xj, fxj, temp; /*mxGetPr, mxGetPr, mxCalloc(double), mxCalloc(double), mxGetPr, mxGetPr, mxGetPr, mxGetPr, x, fx, x1*/
   mwSize mp, pp, i, j, J, Jstart, jj, JJ, ppsize; /*m, p, 0, 0, Jstart, 0, j, J, pp*/
   mwSize p, m; /*mxGetN, mxGetM*/
   bool NewtStep; /*true*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<4 || nrhs>5)
     mexErrMsgTxt("Four or five arguments must be passed");
   if (nlhs>2)
     mexErrMsgTxt("ArrayInv produces at most two outputs");
   if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))       
     mexErrMsgTxt("Input arguments of inproper type");
   if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   if (nlhs==2) 
     if (nrhs<5) 
       mexErrMsgTxt("Must pass a Jacobian to compute the newton step");
     else
       NewtStep=true; 
   else NewtStep=false;
   m=mxGetM(prhs[0]);
   p=mxGetN(prhs[0]);
   mp=m*p;
   if (mxGetNumberOfElements(prhs[1])!=mp)
     mexErrMsgTxt("xl is of improper size");
   if (mxGetNumberOfElements(prhs[2])!=mp)
     mexErrMsgTxt("xu is of improper size");
   if (mxGetNumberOfElements(prhs[3])!=mp)
     mexErrMsgTxt("fx is of improper size");

   x=mxGetPr(prhs[0]);
   xl=mxGetPr(prhs[1]);
   xu=mxGetPr(prhs[2]);
   plhs[0]=mxDuplicateArray(prhs[3]);
   fx=mxGetPr(plhs[0]);
   if (NewtStep) 
   {
     if (mxGetNumberOfElements(prhs[4])!=mp*p)
       mexErrMsgTxt("Inputs are not compatible");
     Jac=mxGetPr(prhs[4]);
     plhs[1]=mxCreateDoubleMatrix(m,p,mxREAL);
     Newton=mxGetPr(plhs[1]);
     A=mxCalloc(pp,sizeof(double));
     b=mxCalloc(p,sizeof(double));
     pp=p*p;
     ppsize=pp*sizeof(double);
     Jstart=0;
     for(i=0; i<m; i++, Jstart++)
     {
       memset(A,0,ppsize);
       for (j=0, J=Jstart; j<p; j++, J+=m) 
       {
         temp=xl[J]-x[J];
         if (fx[J]<=temp)
         {
            fx[J]=temp;
            b[j]=temp;
            A[(p+1)*j]=-1;             /* diagonal element =-1 */
         }
         else 
         {  
           temp=xu[J]-x[J];
           if (fx[J]>=temp)
           {
             fx[J]=temp;
             b[j]=temp;
             A[(p+1)*j]=-1;            /* diagonal element =-1 */
           }
           else
           {
             b[j]=fx[J];
             for (jj=j, JJ=J;jj<pp; jj+=p, JJ+=mp) A[jj]=Jac[JJ];
           }
         }
       }
       backslash(A,b,b,p,1);
       for (j=0, J=Jstart; j<p; j++, J+=m) Newton[J]=-b[j]; 
     }
     mxFree(A);
     mxFree(b);
   }
   else
   {
     for(j=0; j<mp; j++)       
     {
       xj=x[j];
       fxj=fx[j];
       temp=xl[j]-xj;
       if (fxj<=temp) fx[j]=temp;
       else 
       {  
         temp=xu[j]-xj;
         if (fxj>=temp) fx[j]=temp;
       }
     }
   }
}
