#include "mex.h"
#include <math.h>
/*
% BS Black-Scholes option pricing model
% USAGE
%   V=bs(sigma,S,K,r,delta,T,put);
% INPUTS
%   sigma : volatility
%   S     : price of underlying asset
%   K     : strike price
%   r     : interest rate on a discount bond
%   delta : dividend rate
%   T     : time to expiration
% OUTPUT
%   V     : option premium
*/

double bs(double sigma,double S, double k, double r, double d, double T, bool put);
double cdfn(double x);

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwSize n, nn, i; /*nn is mxGetNumberOfElements, n and i derive from nn.*/
  double *V, *Vend, *temp; /*mxGetPr, V, mxGetPr*/
  double *sigma, *S, *k, *r, *d, *T; /*mxGetPr,...,mxGetPr*/
  int nsigma, nS,nk,nr,nd,nT,nput, nv, nd1; /*1 or 0*/
  bool *putflag; /*true*/

  /* Error checking on inputs */  
  if (nrhs<6) mexErrMsgTxt("Not enough input arguments.");
  if (nrhs>7) mexErrMsgTxt("Too many input arguments.");
  if (nlhs>2) mexErrMsgTxt("Too many output arguments.");
  for (i=0; i<nrhs; i++)
  {
    if (!mxIsDouble(prhs[i]) && !mxIsSparse(prhs[i]))
      mexErrMsgTxt("Function not defined for variables of input class");
    if (mxIsComplex(prhs[i]))
      mexErrMsgTxt("X must be real.");
  }

  sigma=mxGetPr(prhs[0]);
  S=mxGetPr(prhs[1]);
  k=mxGetPr(prhs[2]);
  r=mxGetPr(prhs[3]);
  d=mxGetPr(prhs[4]);
  T=mxGetPr(prhs[5]);

  n=1;
  nn=mxGetNumberOfElements(prhs[0]);

  if (n!=nn)
    if (n==1) {n=nn; plhs[0]=mxDuplicateArray(prhs[0]);}  
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nsigma=1; else nsigma=0;

  nn=mxGetNumberOfElements(prhs[1]);
  if (n!=nn)
    if (n==1) {n=nn; plhs[0]=mxDuplicateArray(prhs[1]);}  
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nS=1; else nS=0;
 
  nn=mxGetNumberOfElements(prhs[2]);
  if (n!=nn)
    if (n==1) {n=nn; plhs[0]=mxDuplicateArray(prhs[2]);}  
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nk=1; else nk=0;

  nn=mxGetNumberOfElements(prhs[3]);
  if (n!=nn)
    if (n==1) {n=nn; plhs[0]=mxDuplicateArray(prhs[3]);}  
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nr=1; else nr=0;

  nn=mxGetNumberOfElements(prhs[4]);
  if (n!=nn)
    if (n==1)  {n=nn; plhs[0]=mxDuplicateArray(prhs[4]);} 
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nd=1; else nd=0;

  nn=mxGetNumberOfElements(prhs[5]);
  if (n!=nn)
    if (n==1)  {n=nn; plhs[0]=mxDuplicateArray(prhs[5]);} 
    else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
  if (nn>1) nT=1; else nT=0;
    
  if (nrhs<7) {putflag=mxCalloc(1,sizeof(bool)); *putflag=false; nput=0;}
  else
  {
    nn=mxGetNumberOfElements(prhs[6]);
    if (n!=nn)
      if (n==1) {n=nn; plhs[0]=mxDuplicateArray(prhs[6]);}
      else if (nn>1) mexErrMsgTxt("Inputs are not size compatible");
    if (nn>1) nput=1; else nput=0;
    putflag=mxCalloc(nn,sizeof(bool));   
    temp=mxGetPr(prhs[6]);
    for (i=0;i<nn;i++)
      if (temp[i]==1) putflag[i]=true; else putflag[i]=false;
  }

  if (n==1) plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);

  V=mxGetPr(plhs[0]);
  Vend=V+n;
  for (;V<Vend;V++)
  {
    *V=bs(*sigma,*S,*k,*r,*d,*T,*putflag);
    sigma+=nsigma;
    S+=nS;
    k+=nk;
    r+=nr;
    d+=nd;
    T+=nT;
    putflag+=nput;
  }
}


double bs(double sigma,double S, double k,double r,double d,double T,bool put)
{
  double d1, sigmaT, SS, kk, V;
  sigmaT=sigma*sqrt(T);
  SS=exp(-d*T)*S;
  kk=exp(-r*T)*k;
  d1=log(SS/kk)/sigmaT + 0.5*sigmaT;
  V=SS*cdfn(d1)-kk*cdfn(d1-sigmaT);
  if (put) V=V-SS+kk;
  return(V);
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

  int i;
  double xval, xx, p, q;
  bool NegativeValue;

    xval=x;
    if (xval<0) {xval=-0.7071067811865475*xval; NegativeValue=true; }
    else        {xval= 0.7071067811865475*xval; NegativeValue=false;}
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
    else if (xval<=4)
    {
      xx=xval;
      p=p2[8];
      q=1.0;
      for (i=7;i>=0;i--) {p=p*xx+p2[i]; q=q*xx+q2[i];}
      xx=p/q;
      xx = exp(-xval*xval)*xx;
      if (NegativeValue) xx=0.5*xx; else xx=1-0.5*xx;
    }
    else if (xval<10)
    {
      xx=1/(xval*xval);
      p=p3[5];
      q=1.0;
      for (i=4;i>=0;i--) {p=p*xx+p3[i]; q=q*xx+q3[i];}
      xx=p/q;
      xx = exp(-xval*xval)*(0.7071067811865475-xx)/(xval);
      if (mxIsNaN(xx)) xx=0.0;
      if (NegativeValue) xx=0.5*xx; else xx=1-0.5*xx;
    }
    else 
      if (NegativeValue) xx=0; else xx=1;
    return(xx);
 }
