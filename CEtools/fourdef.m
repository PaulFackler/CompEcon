% FOURDEF Defines parameters for Fourier basis functions
% USAGE
%  [n,parms]=fourdef(m,rho,beta);
%
% See also: FOURBAS, FOURNODE, FOURDOP, FUNDEF

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [n,a,b,parms]=fourdef(m,rho,beta);

if nargin==0
  error('The order of the approximation must be specified')
elseif nargin==1
  rho=1;
  beta=0;
else
  beta=0;
end

if m<0 | m~=fix(m)
  error('m must be a non-negative integer')
end
if rho<=0
  error('rho must be a positive scalar')
end
if beta<0 | beta~=fix(beta)
  error('beta must be a non-negative integer')
end

n=m+m+1;
a=-inf;
b=inf;
parms={m,rho,beta};