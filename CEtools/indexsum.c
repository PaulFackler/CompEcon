#include "mex.h"
#include <math.h>
/* 
  P=indexsum(q,ind);
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   double *q, *ind, *P; /*mxGetPr mxGetPr, mxGetPr*/
   mwSize i, k, mn; /*0, 0, m*/
   int indk; /*int*/
   mwSize m, n;/*mxGetM, mxGetN*/
    
   if (nrhs<2 || nrhs>2)
       mexErrMsgTxt("Two arguments must be passed");
   if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   if (nlhs>1) 
       mexErrMsgTxt("Only one output is returned by indexsum");
   m=mxGetM(prhs[0]);
   n=mxGetN(prhs[0]);
   if (mxGetM(prhs[1])!=m || mxGetN(prhs[1])!=n)
       mexErrMsgTxt("Inputs must be the same size");

   q=mxGetPr(prhs[0]);
   ind=mxGetPr(prhs[1]);
   plhs[0]=mxCreateDoubleMatrix(m,m,mxREAL);
   P=mxGetPr(plhs[0]);
   
   mn=m*n;
   for (k=0; k<mn;)
     for (i=0; i<m; i++, k++)
     {
       indk=(int)ind[k];
       if (indk<1 || indk>m)
         mexErrMsgTxt("IND contains improper values");
       P[i+((int)ind[k]-1)*m]+=q[k];
     }
}
