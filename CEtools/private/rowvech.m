% ROWVECH computes veched matrix outproducts
% SYNTAX
%   V=rowvech(sigma,m,n);
% Given vec(sigma_i)' for i=1:N
% returns vech(sigma_i*sigma_i')' for i=1:N
% with the diagonal elements divided by 2
%
% Used by SCSOLVE
% Coded as a C Mex file

% Copyright (c) 2000 by Paul L. Fackler

function V=rowvech(sigma,m,n)
 
global CompEcon_MEXwarned
 
if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

if size(sigma,2)~=m*n
  error('Inputs are incompatible');
end
N=size(sigma,1);
V=zeros(N,m*(m+1)/2);
ind=find(speye(m));
for i=1:N
  temp=reshape(sigma(i,:),m,n); 
  temp=temp*temp';
  temp(ind)=temp(ind)/2;
  V(i,:)=vech(temp)';
end

