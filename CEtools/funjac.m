% FUNJAC Computes the Jacobian for FUN functions 
% USAGE
%   J=funjac(g,gdef,x);
% INPUTS
%   c        : a coefficient matrix (n x p)
%   fspace   : a family definition structure
%   x        : a set of evaluate points (m x d)
% OUTPUTS
%   J        : an (m x d) array
%              J(i,k) is the derivative with respect to the kth
%              input variable evaluated at x(i,:)

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function J=funjac(c,fspace,x)
d = fspace.d;
J = funeval(c,fspace,x,eye(d));
