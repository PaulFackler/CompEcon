% ArrayMult Computes matrix products over 3-D arrays
% SYNTAX:
%   c=arraymult(a,b,n,p,q,r);
% Inputs:
%   a : n x p x q array
%   b : n x q x r array
%   n,p,q,r : scalar  
% Output:
%   c : n x p x r array
%
%   c(i,:,:) = a(i,:,:)*b(i,:,:), i=1,...,n


function c=arraymult(a,b,n,p,q,r)

global CompEcon_MEXwarned
 
if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

if prod(size(a))~=n*p*q
  error('A is of improper size');
end
if prod(size(b))~=n*q*r
  error('B is of improper size');
end

a=reshape(a,n,p,q);
b=reshape(b,n,q,r);

c = reshape(repmat(permute(a,[1 3 2]),[1,1,r]),n,q,p,r).* ...
    reshape(repmat(b,[1,p,1]),n,q,p,r);
c = reshape(sum(c,2),n,p,r);
