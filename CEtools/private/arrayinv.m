% ARRAYINV Computes linear solves on arrays
% USAGE
%   y=arrayinv(r,rx);
% Given r (m by p) and rx (m by p by p)
% returns y (m by p) with y(i,:)=rx(i,:,:)\r(i,:)

% Copyright (c) 2000 by Paul L. Fackler

function y=arrayinv(r,rx);

global CompEcon_MEXwarned
 
if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

[m,p]=size(r);
y=zeros(m,p);
for i=1:m
  y(i,:)=(reshape(rx(i,:,:),p,p)\reshape(r(i,:),p,1))';
end