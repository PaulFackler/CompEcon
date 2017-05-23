#include "mex.h"
#include <math.h>
/* Performs the Fischer-Burmeister semismooth transformation */

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
 {
   double *x, *a, *b, *f;
   double dminus, d, dplus, da, db;
   double y, z, g, h, temp1, temp2;
   double *fhat, *ffout, *aaout, ff, aa, bb;
   mwSize j;
   int i;
   mwSize n, na, nb, nf; /*mxGetNumberOfElements*/
   mxClassID  category; 

   /* Check input & output number and types */   
   if (nrhs<4 && nrhs>2)
     mexErrMsgTxt("4 input arguments must be passed");
   for (i=0;i<3;i++)
     if (!mxIsDouble(prhs[i]) || mxIsSparse(prhs[i]))
       mexErrMsgTxt("All inputs must be full matrices");
   if (nlhs>3)
     mexErrMsgTxt("Only three output arguments are produced");

   n  = mxGetNumberOfElements(prhs[0]);         /* problem size */
   na = mxGetNumberOfElements(prhs[1]);         
   nb = mxGetNumberOfElements(prhs[2]);         
   nf = mxGetNumberOfElements(prhs[3]);         

   /* Size Checks */
   if (na!=n && na!=1)
     mexErrMsgTxt("Lower bound matrix is improperly specified");
   if (nb!=n && nb!=1)
     mexErrMsgTxt("Upper bound matrix is improperly specified");
   if (nf!=n)
     mexErrMsgTxt("Function value matrix is improperly specified");

   /* get input pointers */
   x=mxGetPr(prhs[0]);
   a=mxGetPr(prhs[1]);
   b=mxGetPr(prhs[2]);
   f=mxGetPr(prhs[3]);

   /* convert # of elements to pointer increments */
   if (na>1) na=1; else na=0; 
   if (nb>1) nb=1; else nb=0; 

   /* allocate memory for output */
   plhs[0]= mxDuplicateArray(prhs[3]);
   fhat=mxGetPr(plhs[0]);
   if (nlhs>1)
   {
     plhs[1]=mxDuplicateArray(prhs[3]);
     ffout=mxGetPr(plhs[1]);
   }
   if (nlhs>2)
   {
     plhs[2]=mxDuplicateArray(prhs[3]);
     aaout=mxGetPr(plhs[2]);
   }

   /* loop over all elements */
   for (j = 0; j < n; j++)
   {
     /* compute phi+ */
     if (mxIsInf(*a)) d=f[j]; 
     else
     {
       da=*a-x[j];
       if (fabs(f[j])>fabs(da)){ y=f[j];  z=da;  }
       else                    { y=da;    z=f[j];}
       z/=y;
       dplus=sqrt(1+z*z);
       if (y>0) d=y*(1+dplus+z);
       else     d=y*(z-((1-dplus)*(1-dplus)+z*z)/dplus/2);
     }
     /* compute phi- */
     if (mxIsInf(*b)) {fhat[j]=d;}
     else
     {
       db=*b-x[j];
       if (fabs(d)>fabs(db)) { g=d;  h=db;}
       else                  { g=db; h=d; }
       h/=g;
       dminus=sqrt(1+h*h);
       if (g<0) fhat[j]=g*(1+dminus+h);
       else     fhat[j]=g*(h-((1-dminus)*(1-dminus)+h*h)/dminus/2);
     } 

     /* compute Jacobian factors if requested */
     if (nlhs>1)
     {
       if (mxIsInf(*b)) { ff=1; aa=1; bb=0; }
       else  
       {
         if (g<0) dminus=-dminus;
         temp1=1-1/dminus;  temp2=1-h/dminus;
         if (fabs(d)>fabs(db)) { ff=temp1; aa=temp1; bb=temp2;}
         else                  { ff=temp2; aa=temp2; bb=temp1;}
       }
       if (mxIsInf(*a)) {aa=0;}
       else
       {
         if (y<0) dplus=-dplus;
         temp1=1+1/dplus; temp2=1+z/dplus;
         if (fabs(f[j])>fabs(da)) { ff*=temp1; aa*=temp2;}
         else                     { ff*=temp2; aa*=temp1;}
       }
       ffout[j]=ff;
       if (nlhs>2) aaout[j]=aa+bb;
     }
     a+=na;
     b+=nb;
   }
 }
