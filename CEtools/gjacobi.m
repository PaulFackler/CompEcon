% GJACOBI Solves Ax=b using Jacobi iteration
% USAGE
%   [x,flag,it] = gjacobi(A,b,x);
% INPUTS
%   A      : nxn matrix
%   b      : n-vector
%   x      : n-vector of starting values, default=b
% OUTPUT
%   x      : approximate solution to Ax=b
%   flag   : 0 if converged
%            1 if maximum iterations exceeded
%   it     : iteration count
%
% USER OPTIONS (SET WITH OPSET)
%   maxit : maximum number of iterations
%   tol   : convergence tolerence

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,flag,it] = gjacobi(A,b,x)

if nargin<3, x=b; end
maxit = optget('gjacobi','maxit',1000);
tol   = optget('gjacobi','tol',sqrt(eps));
flag=0;

n=size(A,1);
Q = sparse(1:n,1:n,diag(A),n,n);  % Diagonalize the diagonal of A
for it=1:maxit
   dx = Q\(b-A*x);
   x = x + dx;
   if norm(dx)<tol, return; end
end

if nargout<=1
  error('Maximum iterations exceeded in gjacobi')
else
  flag=1;
end

