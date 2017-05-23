#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
/*
% LUSOLVE LU solver: x(colindex)=U\(L\(b(rowindex)));
% USAGE
%   x=lusolve(L,U,rowindex,colindex,b);
% INPUTS
%   L        : lower triangular matrix (n x n) (dense or sparse)
%   U        : upper triangular matrix (n x n) (dense or sparse)
%   rowindex : optional row permutation index vector (n x 1)
%   colindex : optional column permutation index vector (n x 1)
%   b        : dense matrix (n x m)
% OUTPUT
%   x        :  dense matrix (n x m)
% 
% Note: no checks are made on the triangularity of L or U. Only the
% lower part of L and the upper part of U are used. 
% Also, no checks are made to ensure that rowindex and colindex are
% valid permuation vectors (i.e., are composed of a permutation 
% of the integers 1,..,n).
% 
% Use LUSETUP to obtain the first 4 inputs.

% 6/4/09 Changed ordering of inputs to group the LU information together
%        Also changed to allow b to be a matrix.
*/



double *spforesub(double *x, double *Pr, double *b, mwSize *Ri, mwSize *Cj, mwSize n)
{
mwSize k, j, rik, cjj; /*0, 0, Ri, Cj*/
double xj; /*x*/
  if (x!=b) memcpy(x,b,n*sizeof(double));
  for (j=0, k=0; j<n; j++){
    cjj=*(++Cj);
    for (; k<cjj; k++){
      rik=Ri[k];
      if (rik>j) x[rik] -= Pr[k]*xj;
      else if (rik==j){x[j] /= Pr[k]; xj=x[j];}
    }
  }
  return(x);
}


double *foresub(double *x, double *Pr, double *b, mwSize n)
{
mwSize k, j; /*j, 0*/
double xj; /*x*/
  if (x!=b) memcpy(x,b,n*sizeof(double));
  for (j=0; j<n; j++){
    Pr+=j;
    x[j] /= *Pr++; 
    xj=x[j];
    for (k=j+1; k<n; k++) x[k] -= *Pr++*xj;
  }
  return(x);
}

double *spbacksub(double *x, double *Pr, double *b, mwSize *Ri, mwSize *Cj, mwSize n)
{
mwSize k, j, rik, cjj; /*Cj, n, Ri, Cj*/
double xj; /*x*/
  if (x!=b) memcpy(x,b,n*sizeof(double));
  Cj+=n;
  xj=x[n-1];
  for (j=n, k=*Cj-1; j>=1; j--){
    cjj=*(--Cj);
    for (; k>=cjj; k--){
      rik=Ri[k];
      if (rik<j-1) x[rik] -= Pr[k]*xj;
      else if (rik==j-1){x[j-1] /= Pr[k]; xj=x[j-1];}
    }
  }
  return(x);
}

double *backsub(double *x, double *Pr, double *b, mwSize n)
{
mwSize k, j; /*j, n*/
double xj; /*x*/
  if (x!=b) memcpy(x,b,n*sizeof(double));
  Pr+=n*n-1;
  for (j=n; j>=1; j--){
    x[j-1] /= *Pr--; 
    xj=x[j-1];
    for (k=j-1; k>=1; k--) x[k-1] -= *Pr--*xj;
    Pr-=n-j+1;
  }
  return(x);
}


void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  /* ***************** */
  /* Declare variables */
  /* ***************** */
  double *x, *b, *Pr, *pindex; /*mxGetPr, mxGetPr, mxGetPr, mxGetPr*/
  mwSize n, i; /*0, 0*/
  mwIndex *Ri, *Cj;/*mxGetIr, mxGetJc*/
  mwSize j, jn; /*0, j */
  mwSize m; /*mxGetN */
  /* ********************************************** */
  /* Determine input sizes and perform error checks */
  /* ********************************************** */
  if (nrhs<5 || nrhs>5)
    mexErrMsgTxt("Five arguments must be passed");
  if (nlhs>1)
    mexErrMsgTxt("Only one output is produced by LUSolve");
  for (n=0; n<nrhs; n++){
    if (!(mxIsDouble(prhs[n]) || (mxIsSparse(prhs[n]) && n<2)))       
      mexErrMsgTxt("Input arguments of inproper type");
    if (mxIsComplex(prhs[n]))       
      mexErrMsgTxt("Only real inputs are supported");
  }

  n=mxGetM(prhs[0]);
  if (mxGetN(prhs[0])!=n)
    mexErrMsgTxt("First input must be square");

  if (mxGetM(prhs[1])!=n || mxGetN(prhs[1])!=n)
    mexErrMsgTxt("Second input must be square");

  if (mxGetM(prhs[4])!=n)
    mexErrMsgTxt("Inputs are not comformable");
  m = mxGetN(prhs[4]);

  plhs[0]=mxDuplicateArray(prhs[4]);
  x=mxGetPr(plhs[0]);


  /* If a row permutation vector is passed, permute the RHS (b): x=b(rowindex) */

  
    i=mxGetNumberOfElements(prhs[2]);
    if (i>0){
      if (mxGetNumberOfElements(prhs[2])!=n)
        mexErrMsgTxt("Row permutation index is the wrong size");
      b=mxGetPr(prhs[4]);
      pindex=mxGetPr(prhs[2]);
      for (j=0; j<m; j++){
        jn=j*n;
        for (i=0; i<n; i++) x[i+jn]=b[(int)(pindex[i])-1+jn];
      }
    }

for (j=0; j<m;j++){
  jn=j*n;
  /* Forward substitution using the L factor: x=L\x */
  Pr=mxGetPr(prhs[0]);
  if (mxIsSparse(prhs[0])){
    Ri=mxGetIr(prhs[0]);
    Cj=mxGetJc(prhs[0]);
    spforesub(x+jn,Pr,x+jn,Ri,Cj,n);
  }
  else
    foresub(x+jn,Pr,x+jn,n);

  /* Backward substitution using the U factor: x=U\x */
  Pr=mxGetPr(prhs[1]);
  if (mxIsSparse(prhs[1])){
    Ri=mxGetIr(prhs[1]);
    Cj=mxGetJc(prhs[1]);
    spbacksub(x+jn,Pr,x+jn,Ri,Cj,n);
  }
  else
    backsub(x+jn,Pr,x+jn,n);
}

  /* If a column permutation vector is passed, permute the LHS: x(colindex)=x */

    i=mxGetNumberOfElements(prhs[3]);
    if (i>0){
      if (mxGetNumberOfElements(prhs[3])!=n)
        mexErrMsgTxt("Column permutation index is the wrong size");
      b=mxCalloc(n,sizeof(double));
      pindex=mxGetPr(prhs[3]);
      for (j=0; j<m; j++){
        jn=j*n;
        for (i=0; i<n; i++) b[(int)(pindex[i])-1]=x[i+jn];
        memcpy(x+jn,b,n*sizeof(double));
      }
      mxFree(b);
    }
}
