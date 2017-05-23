% GSEIDEL Solves Ax=b using Gauss-Seidel iteration
% USAGE
%   [x,flag,it] = gseidel(A,b,x);
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
% USER OPTIONS (SET WITH OPTSET)
%   maxit  : maximum number of iterations
%   tol    : convergence tolerence
%   lambda : over-relaxation parameter

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,flag,it] = gseidel(A,b,x)

if nargin<3, x=b; end
maxit  = optget('gseidel','maxit',1000);
tol    = optget('gseidel','tol',sqrt(eps));
lambda = optget('gseidel','lambda',1);
flag=0;
if issparse(A)
  Q=tril(A);
  for it=1:maxit
    dx = Q\(b-A*x);
    x = x+lambda*dx;
    if norm(dx)<tol, return; end
  end
else  % uses linsolve function for dense matrices to avoid forming Q and checking for triangularity
  opts.LT = true; 
  for it=1:maxit
     dx = linsolve(A,b-A*x,opts);
     x  = x+lambda*dx;
     if norm(dx)<tol, return; end
  end
end

if nargout<=1
  error('Maximum iterations exceeded in gseidel')
else
  flag=1;
end
