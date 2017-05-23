% ARRAYINVB Computes linear solves on arrays
% SYNTAX
%   y=arrayinvb(r,rx,x,xl,xu);
% Given r, x, xl, and xu (m by p) and rx (m by p by p)
% returns y (m by p) 
% The ith row of y is determined by defining:
%   b(j)=xl(i,j)-x(i,j), A(j,k)=-(j==k)   if r(i,j) <= xl(i,j)-x(i,j)
%   b(j)=r(i,j),         A(j,k)=rx(i,j,k) if xl(i,j)-x(i,j) < r(i,j) < xu(i,j)-x(i,j)
%   b(j)=xu(i,j)-x(i,j), A(j,k)=-(j==k)   if r(i,j) >= xu(i,j)-x(i,j)

% Copyright (c) 2000-2001 by Paul L. Fackler & Mario J. Miranda

function y=arrayinvb(r,rx,x,xl,xu)

global CompEcon_MEXwarned
 
if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

[m,p]=size(r);
y=zeros(m,p);
AA=-eye(p);
for i=1:m
  A=reshape(rx(i,:,:),p,p);
  b=r(i,:)';
  bl=(xl(i,:)-x(i,:))';
  bu=(xu(i,:)-x(i,:))';
  ind1=b<=bl;
  ind2=b>=bu;
  b(ind1)=bl(ind1);
  A(ind1,:)=AA(ind1,:);
  b(ind2)=bu(ind2);
  A(ind2,:)=AA(ind2,:);
  y(i,:)=(A\b)';
end