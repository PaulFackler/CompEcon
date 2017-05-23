% FOURNODE Computes standard nodes for Fourier basis
% USAGE
%   x=fournode(n,a,b);
%
% See also: FOURBAS, FUNBAS, FUNEVAL.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=fournode(m,rho,beta)

  x=linspace(0,rho,m+m+2)';
  x=x(1:end-1);