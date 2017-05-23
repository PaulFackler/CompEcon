#include "mex.h"
#include <math.h>
#include <string.h>
/* 
% CDPRODX  Iterated direct product of cell array times a matrix 
% USAGE:
%   a=cdprodx(B,c,ind);
% The direct product of a set of matrices with the same number of rows is
% equivalent to computing row-wise tensor (Kronecker) products:
%   a(i,:) = (B{ind(1)}(i,:) x ... x B{ind(d)}(i,:))*c
% INPUTS:
%   B : a pxq cell array containing matrices of dimension m by n(i)
%   c : a prod(n(ind(i))) by k matrix
%   ind : d-vector of indices for matrices to select from cell array 
%           (default: all matrices in B, 1:p*q)
% OUTPUT:
%   a : an m by k matrix

% Copyright (c) 1997-2001, Paul L. Fackler

The procedure does not actually expand the full
direct product and hence should not cause memory problems.
When the d matrices are all nxn, the number of flops performed
is O(n^(2d)).
*/
void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /* ***************** */
   /* Declare variables */
   /* ***************** */
   mxArray *pBa;
   double **B, **Bptr; /*z, B*/
   double *BLastPtr; /*B*/
   double *z, *c, *d, *order,  *zptr, *cptr, *pptr; /*mxGetPr, mxGetPr, double, mxGetPr, z, c, Bptr*/
   double *dptr; /*d*/
   double *dstart, *dend, *cstart, *cend, *clast; /*d, d, c, c, c*/
   double tempd; /*cstart*/
   mwIndex **Bir, **Bjc, **RowPtr, *RowPtrLast, *rptr, *rend; /*mxGetIr, mxGetJc, Bir, Bir, RowPtr, rptr*/
   mwIndex **ColPtr, **ColEnd, *ColPtrLast, *ColEndLast; /*ColEnd, Bjc, ColEndLast, ColEnd*/
   mwSize Bn, Brows, cm, cn; /*mxGetM, mxGetM, mxGetM, mxGetN*/
   mwSize prodcol, j; /*1, 0*/
   mwSize tempi, zfactor; /*mxGetM, Brows*/
   bool *isfull, LastFull; /*false, LastFull*/
   /* ********************************************** */
   /* Determine input sizes and perform error checks */
   /* ********************************************** */
   if (nrhs<2)
     mexErrMsgTxt("Two arguments must be passed");
   if (!mxIsCell(prhs[0]) || !mxIsDouble(prhs[1]))       
     mexErrMsgTxt("Input arguments of inproper type");
   c=mxGetPr(prhs[1]);                       /* pointer to c data */
   cm=mxGetM(prhs[1]);                       /* # of rows in c */
   cn=mxGetN(prhs[1]);                       /* # of cols in c */
   if (mxIsEmpty(prhs[0])) mexErrMsgTxt("Inputs cannot be empty");
   Brows=mxGetM(mxGetCell(prhs[0],0));       /* # of rows in the B(i) */
   Bn = mxGetM(prhs[0])*mxGetN(prhs[0]);     /* # of B matrices */
   if (Bn<2) mexErrMsgTxt("At least two matrices must be passed");
   if (nrhs>2)                               /* an ordering index was passed */
   {
     order=mxGetPr(prhs[2]);
     tempi=mxGetM(prhs[2])*mxGetN(prhs[2]);
     for (j=0; j<tempi ; j++)         /* check index validity */         
     {
       if (order[j]>Bn || order[j]<1)
         mexErrMsgTxt("Input arguments of improper type");
       else order[j]-=1;             /* C indexing from Base 0 */
     }
     Bn=tempi;                       /* # of B matrices */
   }
   else                              /* use all elements in cell array */                    
   {
     order=mxCalloc(Bn,sizeof(double));
     for (j=0; j<Bn ; j++) order[j]=j;      
   }  
   /* ******************************************** */
   /* Allocate memory for pointer and other arrays */
   /* ******************************************** */
   B=mxCalloc(Bn,sizeof(z));           /* Ptrs for  B matrices */
   Bptr=mxCalloc(Bn,sizeof(Bptr));     /* Ptrs for current B values */      
   isfull=mxCalloc(Bn,sizeof(bool));
   Bir=mxCalloc(Bn,sizeof(Bir));
   Bjc=mxCalloc(Bn,sizeof(Bjc));
   RowPtr=mxCalloc(Bn,sizeof(RowPtr)); 
   ColPtr=mxCalloc(Bn,sizeof(ColPtr)); 
   ColEnd=mxCalloc(Bn,sizeof(ColEnd));  

   prodcol=1;
   for (j=0; j<Bn; j++)                      /* loop over the B{j} */
   {
     pBa=mxGetCell(prhs[0],order[j]);        /* pointer to B{j} mxArrays */
     B[j]=mxGetPr(pBa);                      /* pointers to B{j} data    */
     if (mxIsSparse(pBa))                    /* B{j} is sparse           */
     {
       isfull[j]=false;
       Bir[j]=mxGetIr(pBa);
       Bjc[j]=mxGetJc(pBa);
       tempi=mxGetN(pBa);                    /* # of cols in B{j} */
       Bptr[j]=B[j]+*(Bjc[j]+tempi);
       ColEnd[j]=Bjc[j]+tempi;
       ColPtr[j]=ColEnd[j];
     }
     else if(mxIsDouble(pBa))                /* B{j} is full */
     {
       isfull[j]=true;
       tempi=mxGetN(pBa);                    /* # of cols in B{j} */
       Bptr[j]=B[j]+Brows*tempi;
       ColEnd[j]=0;
       ColEnd[j]+=tempi;                     /* use ColPtr as an integer */
       ColPtr[j]=ColEnd[j];
     }
     else mexErrMsgTxt("B must be composed of matrices of real numbers");
     if (mxGetM(pBa) != Brows) 
       mexErrMsgTxt("All rows in cell array must be equal"); 
     prodcol=prodcol*tempi;                  /* product of cols */
   } 
   if (prodcol != cm) mexErrMsgTxt("B and c are not conformable");
   plhs[0]=mxCreateDoubleMatrix(Brows,cn,mxREAL); /* mxArray for output */
   z=mxGetPr(plhs[0]);                            /* pointer to output data */
   d=mxCalloc(Brows*Bn,sizeof(double));
   /* ********************************************************************* */
   /* d[i,j] contains the product of variables B(i,k), k=1 to j.            */
   /* Column j of d is updated everytime variable j+1 turns over.           */
   /* The turn over is dtermined by checking ColPtr against ColEnd.         */
   /* The last variable, which turns over fastest, is handled separately to */
   /* increase speed.                                                       */
   /* ********************************************************************* */
   dend=d+Brows;
   for (dptr=d;dptr<dend;dptr++) *dptr=1;
   dstart=d+Brows*(Bn-1);
   ColEndLast=ColEnd[Bn-1];
   ColPtrLast=ColEndLast;
   LastFull=isfull[Bn-1];
   cend=c+cm;
   clast=c+cm*cn;
   zfactor=Brows*cn-1;
   /* **************************** */
   /* The main loop over rows of c */
   /* **************************** */
   for (cstart=c; cstart<cend; cstart++) 
   {
     /* ********************************************** */
     /* Check if iterated through the first matrix     */
     /* Multipliers only updated when an index changes */
     /* A full revision occurs on the first iteration  */
     /* ********************************************** */
     if (ColPtrLast>=ColEndLast)    
     {
       j=Bn-1;
       BLastPtr=B[j];
       if (LastFull) ColPtrLast=0;
       else {RowPtrLast=Bir[j]; ColPtrLast=Bjc[j];}
       do
       {
         j--;
         (ColPtr[j])++;
         if (ColPtr[j]<ColEnd[j]) break;
         else
         {  
           Bptr[j]=B[j];
           if (isfull[j]) ColPtr[j]=0;
           else          {RowPtr[j]=Bir[j]; ColPtr[j]=Bjc[j];}
         }
       }
       while (j>0);
       /* ******************************************************** */
       /* Replace the multipliers (d[i,j]) for variables j to Bn-1 */
       /* ******************************************************** */
       dptr=d+(j+1)*Brows;
       for (;j<Bn-1;j++)
       { 
         if (isfull[j])     /* jth matrix is full */
         {
           pptr=Bptr[j];
           dend=dptr+Brows;
           for (;dptr<dend; dptr++, pptr++)
             *dptr = *(dptr-Brows)**pptr;
           Bptr[j] = pptr;
         }
         else              /* jth matrix is sparse */
         {
           /* fill d(:,j) with zeros */
           memset(dptr,'\0',Brows*sizeof(double));
           pptr=Bptr[j];
           rptr=RowPtr[j];
           rend=rptr+*(ColPtr[j]+1)-*ColPtr[j];
           for (;rptr<rend; rptr++, pptr++)
              *(dptr+*rptr)=*(dptr-Brows+*rptr)**pptr;
           dptr=dptr+Brows;
           Bptr[j]=pptr;
           RowPtr[j]=rptr;
         }
       }
     }
     /* ************************************ */
     /* This is executed every time          */
     /* Loop over rows of B and columns of c */
     /* ************************************ */
     zptr=z;
     if (LastFull)
     {
       dend=dstart+Brows;
       if (cn==1)
       {
         tempd=*cstart;
         for (dptr=dstart; dptr<dend; dptr++,BLastPtr++,zptr++)
           *zptr += tempd*(*dptr**BLastPtr); 
       }
       else
         for (dptr=dstart; dptr<dend; dptr++, BLastPtr++)
         {
           tempd=*dptr * *BLastPtr;
           for (cptr=cstart; cptr<clast; cptr+=cm, zptr+=Brows)
             *zptr += *cptr*tempd; 
           zptr-=zfactor;
         }
     }
     else
     {
       rend=RowPtrLast+*(ColPtrLast+1)-*ColPtrLast;
       if (cn==1)
         for (;RowPtrLast<rend; RowPtrLast++, BLastPtr++)
           *(z+*RowPtrLast) += *cstart*(*(dstart+*RowPtrLast)**BLastPtr); 
       else
         for (;RowPtrLast<rend; RowPtrLast++, BLastPtr++)
         {
           tempd=*(dstart+*RowPtrLast)**BLastPtr;
           zptr=z+*RowPtrLast;
           for (cptr=cstart; cptr<clast; cptr+=cm, zptr+=Brows)
               *zptr += *cptr*tempd; 
         }
     }
     ColPtrLast ++;
   }
   /* ******************************************** */
   /* Main Loop completed - free memory and return */
   /* ******************************************** */ 
   mxFree(B);
   mxFree(Bptr);
   mxFree(isfull);
   mxFree(Bir);
   mxFree(Bjc); 
   mxFree(RowPtr); 
   mxFree(ColPtr); 
   mxFree(ColEnd); 
   mxFree(d);  
   if (nrhs<=2) mxFree(order);    
}
