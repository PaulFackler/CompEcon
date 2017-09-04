%% FIXPOINT
%
%  Computes fixed-point of function using function iteration
% 
%  Usage
%    x = fixpoint(g,x,varargin)
%  Input
%    g         : function of form gval=g(x)
%    x         : initial guess for fixed-point
%    varargin  : optional parameters passed to g
%  Output
%    x         : fixed-point of g
%    gval      : function value at x
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (1e-10)
%    showiters : display results of each iteration (1)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function x = fixpoint(g,x,varargin)

% Set options to defaults, if not set by user with OPTSET (see above)
maxit     = optget('fixpoint','maxit',100);
tol       = optget('fixpoint','tol',1e-10);
showiters = optget('fixpoint','showiters',1);

if showiters, fprintf('\nIn fixpoint:\n'), end

for it=1:maxit
  xold = x;
  x = feval(g,x,varargin{:});
  change = norm(x-xold);
  if showiters, fprintf('%4i %6.2e\n',[it change]); end
  if change<tol
    if showiters, fprintf('\n'), end
    return
  end
end
warning('Failure to converge in fixpoint')