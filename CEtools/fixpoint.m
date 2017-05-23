% FIXPOINT Computes fixpoint of function using function iteration
% USAGE
%   [x,gval] = fixpoint(g,x,varargin);
% INPUTS
%   g       : name of function of form gval=g(x)
%   x       : initial guess for fixpoint
%   varargin: additional arguments for f (optional)
% OUTPUTS
%   x       : fixpoint of g
%   gval    : function value estimate
%
% Setable options (use OPTSET):
%   tol     : convergence tolerance
%   maxit   : maximum number of iterations

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,gval] = fixpoint(g,x,varargin)

maxit = optget('fixpoint','maxit',100);
tol   = optget('fixpoint','tol',sqrt(eps));

for it=1:maxit 
   gval = feval(g,x,varargin{:});
   if norm(gval-x)<tol, return, end 
   x = gval;
end
warning('Failure to converge in fixpoint')
