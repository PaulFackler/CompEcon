#include "mex.h"
#include <math.h>
#include <string.h>
/*
% KERNELX Computes a kernel estimate of a PDF
% USAGE:
%   [f,F]=kernel(x,xi,h,ktype);
% INPUTS:
%   x     : an mx1 vector of evaluation points
%   xi    : an nx1 vector of observations
%   h     : bandwidth (optional)
%   ktype : 0 = Gaussian, 1 = Epanechnikov (optional, default=1)
% OUTPUTS
%   f : an mx1 vector of estimates of the PDF
%   F : an mx1 vector of estimates of the CDF
%
% Example:
%   xi=randn(1000,1);
%   x=linspace(-5,5,101)';
%   figure(1); plot(x,[normpdf(x) kernel(x,xi,.4)]);
*/

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
  int i;
  double xval, xx, p, q;
  bool NegativeValue;
    if (mxIsNaN(x)) return(mxGetNaN());
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
      if (NegativeValue) xx = (1-xx)/2; else xx = (1+xx)/2;
    }
    else if (xval<=5.6568542494923806)
    {
      xx=xval;
      p=p2[8];
      q=1.0;
      for (i=7;i>=0;i--) {p=p*xx+p2[i]; q=q*xx+q2[i];}
      xx=p/q;
      xx = exp(-xval*xval)*xx;
      if (NegativeValue) xx=xx/2; else xx=1-xx/2;
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
      if (NegativeValue) xx=xx/2; else xx=1-xx/2;
    }
    else 
      if (NegativeValue) xx=0; else xx=1;
    return(xx);
 }




void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  const double pi=3.141592653589793238463;
  const double sqrt5=2.236067977499789696409;

  double *f, *F, *x, *xi, *h;
  double factor1, factor2, k, hi, sum, mean, std, fj, Fj, kk;
  int hAdd; /*0 or 1*/
  int ktype; /*0 or 1*/
  double hval; /*0, but needs to be compatible with h.*/
  mwSize i, j, n, m;  /*mxGetNumberOfElements*/
  bool getcdf;

  /* Error checking on inputs */  
  if (nrhs<4) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>4) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>2) mexErrMsgTxt("Too many output arguments.");
  for (i=0;i<nrhs;i++)
  {
    if (!mxIsDouble(prhs[i]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsSparse(prhs[i]))
      mexErrMsgTxt("Function not defined for sparse matrices");
  }

/*   Dimension Checking and initializations */
  n=mxGetNumberOfElements(prhs[0]);
  m=mxGetNumberOfElements(prhs[1]);
  i=mxGetNumberOfElements(prhs[2]);
  if (i==1) hAdd=0;
  else if (i==m) hAdd=1;
  else mexErrMsgTxt("xi and h are incompatible");
  x=mxGetPr(prhs[0]);
  xi=mxGetPr(prhs[1]);
  
  if (mxIsEmpty(prhs[2])) {h=&hval; hval=0;}
  else  h=mxGetPr(prhs[2]);

/* 0 = Gaussian, 1 = Epanechnikov */
  if (mxIsEmpty(prhs[3])) ktype=1;
  else                    ktype=*mxGetPr(prhs[3]); 
  
/* Allocate memory for output arrays */
  plhs[0]=mxDuplicateArray(prhs[0]);
  f=mxGetPr(plhs[0]);
  memset(f,0,n*sizeof(double));
  if (nlhs==2)
  {
    plhs[1]=mxDuplicateArray(plhs[0]);
    F=mxGetPr(plhs[1]);
    getcdf=true;
  }
  else getcdf=false;


if (*h==0)
{
  hAdd=0;
  if (ktype==1) *h=1.0487*pow(m,-0.2);
  else          *h=1.0592*pow(m,-0.2);
}

/* compute the mean */
sum=0; for (i=0; i<m; i++) sum+=xi[i]; mean=sum/m;
/* compute the standard deviation */
sum=0; for (i=0; i<m; i++) {std=(xi[i]-mean); sum+=std*std;}
std=sqrt(sum/(m-1));

if (ktype==1) factor1=0.15/sqrt5/m;
else          factor1=1/sqrt(2*pi)/m;

hi=*h*std;
factor2=factor1/hi;

if (getcdf) /* get PDF and CDF */
for (i=0; i<m; i++) {
  for (j=0; j<n; j++) {
    k=(x[j]-xi[i])/hi;
    if (ktype==1) {
      kk=k*k;
      if (kk<5){
        f[j]+=factor2*(5-kk); 
        F[j]+=factor1*(k*(15-kk)+10*sqrt5)/3;
      }
      else if (k>=sqrt5) F[j]+=1.0/m;
    }
    else {
      f[j]+=factor2*exp(-0.5*k*k); 
      F[j]+=cdfn(k)/m;
    }
  }
  if (hAdd==1) {h++; hi=*h*std; factor2=factor1/hi;}
}
else  /* get PDF only */
if (ktype==1) /* Epinichnikov */
for (i=0; i<m; i++) {
  for (j=0; j<n; j++) {
    k=(x[j]-xi[i])/hi;
    if (fabs(k)<sqrt5)   f[j]+=factor2*(5-k*k); 
  }
  if (hAdd==1) {h++; hi=*h*std; factor2=factor1/hi;}
}
else     /* Gaussian */
for (i=0; i<m; i++) {
  for (j=0; j<n; j++) {
    k=(x[j]-xi[i])/hi;
    f[j]+=factor2*exp(-0.5*k*k); 
  }
  if (hAdd==1) {h++; hi=*h*std; factor2=factor1/hi;}
}
}
