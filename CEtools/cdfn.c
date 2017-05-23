#include "mex.h"
#include <math.h>
/*
% CDFN Computes the CDF of the standard normal distribution
% USAGE
%   p=cdfn(x)
*/

double cdfn(double x);

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize n; /*mxGetNumberOfElements*/
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
  for (pend=p+n; p<pend; p++,x++) { *p=cdfn(*x);} 

}


double cdfn(double x)
{
 const double p1[5] = {
  3.20937758913846947e03,
  3.77485237685302021e02,
  1.13864154151050156e02,
  3.16112374387056560e00,
  1.85777706184603153e-1};
const double q1[4] = {
  2.84423683343917062e03,
  1.28261652607737228e03,
  2.44024637934444173e02,
  2.36012909523441209e01};
const double p2[9] = { 
  1.23033935479799725e03,
  2.05107837782607147e03,
  1.71204761263407058e03,
  8.81952221241769090e02,
  2.98635138197400131e02,
  6.61191906371416295e01,
  8.88314979438837594e00,
  5.64188496988670089e-1,
  2.15311535474403846e-8};
const double q2[8] = { 
  1.23033935480374942e03,
  3.43936767414372164e03,
  4.36261909014324716e03,
  3.29079923573345963e03,
  1.62138957456669019e03,
  5.37181101862009858e02,
  1.17693950891312499e02,
  1.57449261107098347e01};
const double p3[6] = {
  6.58749161529837803e-4,
  1.60837851487422766e-2,
  1.25781726111229246e-1,
  3.60344899949804439e-1,
  3.05326634961232344e-1,
  1.63153871373020978e-2};
const double q3[5] = { 
  2.33520497626869185e-3,
  6.05183413124413191e-2,
  5.27905102951428412e-1,
  1.87295284992346047e00,
  2.56852019228982242e00};
const double sqrt2i=7.071067811865475244e-1;

  int i; /*All, already checked in bs.c*/
  double xval, xx, p, q;
  bool NegativeValue;

    xval=x;
    if (xval<0) {xval=-sqrt2i*xval; NegativeValue=true; }
    else        {xval= sqrt2i*xval; NegativeValue=false;}
    if (xval<=0.46875)
    {
      xx = xval*xval;
      p=p1[4];
      q=1.0;
      for (i=3;i>=0;i--) {p=p*xx+p1[i]; q=q*xx+q1[i];}
      xx=p/q;
      xx *= xval;
      if (NegativeValue) xx = 0.5*(1-xx); else xx = 0.5*(1+xx);
    }
    else if (xval<=5.6568542494923806)
    {
      xx=xval;
      p=p2[8];
      q=1.0;
      for (i=7;i>=0;i--) {p=p*xx+p2[i]; q=q*xx+q2[i];}
      xx=p/q;
      xx = exp(-xval*xval)*xx;
      if (NegativeValue) xx=0.5*xx; else xx=1-0.5*xx;
    }
    else if (xval<8.3)
    {
      xx=1/(xval*xval);
      p=p3[5];
      q=1.0;
      for (i=4;i>=0;i--) {p=p*xx+p3[i]; q=q*xx+q3[i];}
      xx=p/q;
      xx = exp(-xval*xval)*(sqrt2i-xx)/(xval);
      if (mxIsNaN(xx)) xx=0.0;
      if (NegativeValue) xx=0.5*xx; else xx=1-0.5*xx;
    }
    else 
      if (NegativeValue) xx=0; else xx=1;
    return(xx);
 }

