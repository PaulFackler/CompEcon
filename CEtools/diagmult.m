% DIAGMULT Computes either diag(a)*b or a*diag(b)
% USAGE
%   c=diagmult(a,b);
% INPUTS
%   a,b : m-vector and mxn matrix
%         or 
%         mxn matrix and n-vector
% OUTPUT
%   c   : mxn matrix 
%
% Coded as a C-MEX function for speed

% Note: C version is not implemented for complex matrices or matrices 
% with data type other than double. The vector must be full but
% the matrix can be sparse or full.


% Copyright (c) 2005-8, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu

function c=diagmult(a,b)

 [m,n]=size(a);
 % a is a vector
 if m==1 | n==1
     n=length(a);
     if size(b,1)~=n
       error('Inputs are not compatible')
     end
     c=sparse(1:n,1:n,a,n,n)*b; 
 else
   [m,n]=size(b);
   % b is a vector
   if m==1 | n==1
     n=length(b);
     if size(a,2)~=n
       error('Inputs are not compatible')
     end
     c=a*sparse(1:n,1:n,b,n,n);
   else
     error('Either a or b must be vectors')
   end
 end