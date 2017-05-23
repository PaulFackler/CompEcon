
/*
% Copyright (c) 2005, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu
% 
*/


#include "mex.h"
#include <math.h>

#define bignum 1.7e308
#define dmax(x,y) (x>=y ? x :  y)
#define dmin(x,y) (x<=y ? x :  y)
#define dabs(x)   (x>=0  ? x : -x)

/*
DIRECTC Called by DIRECT
This function should NOT be called directly 
*/

/* LCONVHULL 
   Lower convex hull
   Given vectors x and y both of size m
   returns a vector h with the indices of the points that form
   the lower convex hull. The return value is the number of 
   such points. 
   */
int lconvhull(double *x, double *y, int *h, int m){
int a, b, c, i, ii, ha, hb, hc; /*1, a, b, b, i*/
  if (m < 3){ return(m);}
  a = 1;
  while (true){
    b=a+1;
    if (b >= m) break;
    c=b+1;
    ha=h[a]; hb=h[b]; hc=h[c];
    if ((x[hb]-x[ha])*(y[hc]-y[ha])>=(x[hc]-x[ha])*(y[hb]-y[ha]))
      a = b;
    else{
      m--;
      for (i=b, ii=i+1; i<=m; i=ii, ii++){h[i] = h[ii];}
      if (a>1) a--;
    }
  }
  return(m);
}



/* QUICKSORTIND 
   Sorts a vector "a" in place and returns the sorted vector and 
   and an index vector "h". The index vector must be preallocated
   with the number 1 through n. This function calls itself recursively
   and should initially be called with l=1 and r=n.
   */
void quicksortind(double *a, int *h, int l, int r){
  int i, j, itmp; /*1, r, h*/
  double x, rtmp; /*a, a*/
  i=l; j=r; x=a[(l+r)/2];
  do{
    while (a[i]<x) i++;
    while (x<a[j]) j--;
    if (i<=j){
      rtmp=a[i]; a[i]=a[j]; a[j]=rtmp;
      itmp=h[i]; h[i]=h[j]; h[j]=itmp;
      i++; j--;
    }
  }
  while (i<=j);
  if (l<j) quicksortind(a,h,l,j);
  if (i<r) quicksortind(a,h,i,r);
}



/* DIRECT 
   The main program. Performs the iteration loop of the algorithm 
   */
void direct(nargs, inarg, outarg, x, xmid, xrange, 
             parent, len, dim, S, F, s, f, parameters, intarray)
int nargs; 
mxArray **inarg, **outarg; 
double *x, *xmid, *xrange, *len;
int *parent;
short *S, *dim;
double *F, *s, *f, *parameters, *intarray;
{
double  epsilon, fopt, tol, minind; /*parameters, parameters, parameters, parameters*/
int nevals, it, ns,  ibest, maxit, ShowIters, n, maxeval, maxs, minevals, minsize; /*int, ..., int*/
double *doublemem, *fnew1, *fnew2, *w, *lengths; /*mxCalloc(double), doublemem, doublemem, doublemem, 1.0*/
int *intmem, *R, *I, *ind, *heads, *tails, *llist; /*mxCalloc(int), intmem, ind, R, I, heads, tails*/
double delta, xi, fbest, minrat, tmp, E, Fj, dx; /*lengths, x, F, bignum, f, dmax, F, delta*/
int i, ii, j, jj, nI, nR, k, convtype; /*1, maxs, 1, 1, n, 0, 0, 0*/
short Sj; /*S*/
bool getall, notgreedy; /*true, true*/ 
  /* Get parameter values */
  epsilon   = parameters[0];
  fopt      = parameters[1];
  tol       = parameters[2];
  minind    = parameters[3];
  nevals    = (int)intarray[0];
  it        = (int)intarray[1];
  ns        = (int)intarray[2];
  ibest     = (int)intarray[3];
  maxit     = (int)intarray[5];
  ShowIters = (int)intarray[6];
  n         = (int)intarray[7];
  maxeval   = (int)intarray[8];
  maxs      = (int)intarray[9];
  minevals  = (int)intarray[10];
  minsize   = (int)intarray[11];
  if (((int)intarray[12])!=1) getall=true;    else getall=false;
  if (((int)intarray[13])!=1) notgreedy=true; else notgreedy=false;
  if (minind==-1) for (i=1; i<=nevals; i++) F[i]=-F[i];
  fbest     = F[ibest];
  fopt      = minind*fopt;

  /* allocate working memory */
  intmem  = mxCalloc(2*maxeval+3*maxs+n+1,sizeof(int));
  R       = intmem;
  ind     = R+maxeval;
  I       = ind+maxs;
  heads   = I+n;
  tails   = heads+maxs;
  llist   = tails+maxs;
  doublemem = mxCalloc(3*n+maxs+1,sizeof(double));
  w       = doublemem;
  fnew1   = w+n;
  fnew2   = fnew1+n;
  lengths = fnew2+n;

  /* Create vector of possible side lengths */
  lengths[1] = 1.0/3.0; 
  ii=(maxs-1)/n+1; for (i=1;i<=ii;i++) {lengths[i+1]=lengths[i]/3.0;}
  convtype=0;
  while (it < maxit  &&  convtype==0){
    it++;
    /* Compute best function values for each size and store in linked lists*/
    for (j=1; j<=ns; j++){f[j]=bignum;}
    for (j=1; j<=nevals; j++){
      Fj=F[j]; Sj=S[j];
      if (Fj<f[Sj])
        {f[Sj]=Fj; heads[Sj]=j; tails[Sj]=j; llist[j]=0;}
      else if (Fj==f[Sj] && getall)
        {llist[tails[Sj]]=j; tails[Sj]=j; llist[j]==0;}
/*
      else if (Fj==f[Sj] && llist[j]==0 && tails[Sj]!=j)
        {llist[tails[Sj]]=j; tails[Sj]=j;}
*/
    }
    /*  Identify the set R of potentially optimal rectangles */
    E = dmax(epsilon*dabs(fbest),1E-8);
    minrat=bignum; k=0;
    for (i=S[ibest]; i>=1; i--){
      if (f[i]<bignum){
        tmp=((f[i] - fbest) + E)/s[i];
        if (tmp<minrat) {minrat=tmp; k=0;}
        ind[++k]=i;
      }
    }
    if (notgreedy) k=lconvhull(s,f,ind,k);
    nR=0;
    for (i=1; i<=k; i++){
      j=heads[ind[i]]; 
      while (j!=0){R[++nR]=j; j=llist[j];}
      /*f[ind[i]]=bignum;*/
    }
    /*  Loop over potentially optimal rectangles*/
    for (jj = 1; jj<=nR; jj++){ 
      j = R[jj];
      /* determine the number of long sides for rectangle j */
      nI = n - ((S[j]-1)%n);
      /* not enough room to perform iteration */
      if (nevals+2*nI>maxeval) {convtype=1; break;}  
      /* determine which sides are long */
      if (nI==n) for (i=1; i<=n; i++)  I[i]=i;
      else       for (i=1; i<=nI; i++) I[i]=dim[j+2*i];
      /* determine the length of the long side divided by 3 */
      delta = lengths[(S[j]-1)/n+1];
      /* get the center of the current point */
      for (i=1; i<=n; i++) x[i]=xmid[i];
      i=j; k=0; while (i>1) {ind[++k]=i; i=parent[i];}
      while (k>=1) {ii=ind[k--]; x[dim[ii]]+=len[ii];}
      /* loop over each long side dimension and evaluate new points */
      for (ii = 1; ii<=nI; ii++){    
         i = I[ii]; 
         xi=x[i];
         dx=delta*xrange[i];
         /* First new evaluation */
         x[i] += dx;                         /* Centerpoint for new rectangle */
         mexCallMATLAB(1,outarg,nargs,inarg,"feval");
         fnew1[ii] = minind**mxGetPr(outarg[0]);
         mxDestroyArray(outarg[0]);           /* destroy return array */
         /* Second new evaluation */
         x[i] = xi-dx;                        /* Centerpoint for new rectangle */
         mexCallMATLAB(1,outarg,nargs,inarg,"feval");
         fnew2[ii] = minind**mxGetPr(outarg[0]);
         mxDestroyArray(outarg[0]);          /* destroy return array */ 
         x[i] = xi;
         w[ii] = dmin(fnew1[ii],fnew2[ii]);
      }
      /* sort the new points according to functon values */
      for (i=1; i<=nI; i++) ind[i]=i;
      quicksortind(w,ind,1,nI);
      /* check if best new point improves on current best point */
      if (w[1]<fbest){
        fbest=w[1];
        if (fbest==fnew1[ind[1]]) ibest=nevals+1;
        else                      ibest=nevals+2; 
      }
      Sj=S[j];     /* size of rectangle j */
      /* divide rectangle j into 2nI+1 new rectangles */
      for (ii = 1; ii<=nI; ii++){
        Sj++;
        k = ind[ii];              /* position in cut list */
        i = I[k];                 /* dimension of cut*/
        dx=delta*xrange[i];
        /* store information about new rectangles */
        parent[++nevals] = j;  len[nevals] = dx;  
        dim[nevals] = i;         F[nevals] = fnew1[k];         S[nevals] = Sj;
        parent[++nevals] = j;  len[nevals] = -dx;  
        dim[nevals] = i;         F[nevals] = fnew2[k];         S[nevals] = Sj;
      }
      if (Sj>ns) ns=Sj;
      S[j] = Sj; /*llist[j]=0; */
    }
    if (ns>maxs){convtype=2; ns=maxs;} /*no more storage space-increase maxlen */
    /*  Check for convergence */
    if (fbest <= fopt) convtype = 3;
    else if ((fbest-fopt)/dabs(fopt) < tol) convtype = 4;
    if (S[ibest]>=minsize && nevals>=minevals) convtype = 5;
    /* Print iteration results if requested */
    if (ShowIters == 1){
      printf("Iter: %4i   fbest: %15.10f    fn evals: %8i\n",
                  it,minind*fbest,nevals); }
  }
  if (it>=maxit) convtype=6;
  /* Update input array before returning */
  intarray[0]=nevals;
  intarray[1]=it;
  intarray[2]=ns;
  intarray[3]=ibest;
  intarray[4]=convtype;
  mxFree(doublemem);
  mxFree(intmem);
}



/* MATLAB Gateway function
   */
void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
double  *xrange, *F, *d, *f, *parameters, *intarray; /*mxGetPr, ..., mxGetPr*/
double *x, *xmid, *len; /*mxGetPr, ..., mxGetPr*/
short *S, *dim; /*mxGetData, mxGetData (This is highly dependent on...*/
int *parent; /*mxGetData ...the type of the pointers, but should work)*/
mwSize maxeval, maxs, n; /*mxGetN, mxGetNumberOfElements, mxGetNumberOfElements,*/ 
int nargs, i; /*nrhs, 2*/
const mxArray **inarg; 
mxArray **outarg; /*mxCreateDoubleMatrix, outarg in main()*/
  /* Error checking on inputs */  
  if (nrhs<12) mexErrMsgTxt("Incorrect number of input arguments.");
  if (nlhs!=0)  mexErrMsgTxt("Incorrect number of output arguments.");
  if (mxGetNumberOfElements(prhs[10])!=4)
     mexErrMsgTxt("Parameters vector is the wrong size.");
  if (mxGetNumberOfElements(prhs[11])!=14) 
     mexErrMsgTxt("Intarray vector is the wrong size.");
  if (sizeof(short)!=2) 
     mexErrMsgTxt("Must be able to read 2 byte integers.");
  /* Determine storage sizes */
  maxeval = mxGetN(prhs[7]);
  maxs    = mxGetNumberOfElements(prhs[8]);
  n       = mxGetNumberOfElements(prhs[1]);
  /* Get data pointers */
  xrange  = mxGetPr(prhs[1]);
  xmid    = mxGetPr(prhs[2]);
  parent  = mxGetData(prhs[3]); 
  len     = mxGetPr(prhs[4]); 
  dim     = mxGetData(prhs[5]); 
  S       = mxGetData(prhs[6]); 
  F       = mxGetPr(prhs[7]);
  d       = mxGetPr(prhs[8]);
  f       = mxGetPr(prhs[9]);
  parameters = mxGetPr(prhs[10]);  
  intarray = mxGetPr(prhs[11]); 
  /* Set up function calling array */
  nargs    = nrhs-10;                          /* # of arguments to pass */  
  inarg    = mxCalloc(nargs,sizeof(inarg));    /* create pointer array */
  inarg[0] = prhs[0];                          /* populate pointer array */
  inarg[1] = mxCreateDoubleMatrix(n,1,mxREAL); /* populate pointer array */
  for (i=2; i<nargs;i++) inarg[i]=prhs[i+10];
  outarg   = mxCalloc(1,sizeof(outarg));       /* create return pointer array */
  x        = mxGetPr(inarg[1]);                /* location of evaluation point */ 
  /* call main algorithm */
  direct(nargs, inarg, outarg, x-1, xmid-1, xrange-1, 
             parent-1, len-1, dim-1, S-1, F-1, d-1, f-1, parameters, intarray);
  /* Clean up memory */
  /*mxDestroyArray(inarg[1]); the Apple compiler didn't like destroying a constant array. It will be freed automatically anyway.*/
  mxFree(outarg);
  mxFree(inarg);
}
