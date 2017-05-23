#include <string.h>

double *backslash(double *A, double *b, double *x, int n, int p);
/* x=A\b where A is nxn and b is nxp. x can overwrite b */


double *backslash(double *A, double *b, double *x, int n, int p)
{
  double *U;
  double Amax, Lik, temp;
  int i, j, k, ki, nn, np, imax, kn;

  nn=n*n;
  np=n*p;
  U=mxCalloc(nn,sizeof(double));
  memcpy(U,A,nn*sizeof(double));
  if (b!=x) memcpy(x,b,np*sizeof(double));
 
  /* Scale the inputs */
  for (i=0; i<n; i++)
  {
    Amax=0;
    for (kn=0; kn<nn; kn+=n) {temp=fabs(U[i+kn]); if (temp>Amax) Amax=temp;}
    for (kn=0; kn<nn; kn+=n) U[i+kn]/=Amax;
    for (kn=0; kn<np; kn+=n) x[i+kn]/=Amax;
  }

  kn=0;
  for (k=0; k<n-1; k++)    
  {
    /* Determine the maximal element in column at or below row j */
    Amax=0; imax=k;
    for (i=k; i<n; i++)
    {
      temp=fabs(U[i+kn]); 
      if (temp>Amax){Amax=temp; imax=i;}
    }
    if (Amax==0)
    {
       mexWarnMsgTxt("Backslash: Matrix is singular to working precision.");
       temp=mxGetInf();
       for (i=0, imax=n*p; i<imax; i++) x[i]=temp;
       mxFree(U);
       return(x);
    }
    /* Perform row interchanges if necessary */
    if (imax>k)
    {
      for (j=0; j<nn; j+=n)
      { temp=U[imax+j]; U[imax+j]=U[k+j]; U[k+j]=temp; }  /* swap rows of A */
      for (j=0; j<np; j+=n)
      { temp=x[imax+j]; x[imax+j]=x[k+j]; x[k+j]=temp; }  /* swap rows of x */
    }
    for (i=k+1; i<n; i++)
    {
      Lik=U[i+kn]/U[k+kn];
      for (j=kn+n; j<nn; j+=n) U[i+j]-=Lik*U[k+j];
      for (j=0; j<np; j+=n) x[i+j]-=Lik*x[k+j];
    }
    kn+=n;
  }
  kn=nn-n;
  /* Back substitution */
  for (k=n-1;k>=0;k--)
  { 
    for(j=0; j<np; j+=n) 
    { 
      temp=x[k+j];
      for (i=k+1, ki=k+i*n; i<n; i++, ki+=n) temp-=U[ki]*x[i+j];
      x[k+j]=temp/U[k+kn];
    }
    kn-=n;
  }
  mxFree(U);
  return(x);
}
