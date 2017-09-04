%% NEWTOND   
%
%  Computes root of function from R^n to R^n when f_i depends exclusively 
%  on x_i, using Newton's method with backstepping.
%
%  Usage
%   [x,fval] = newtond(f,x,varargin)
%  Input
%    f         : function of form [fval,fder]=f(x,varargin) where fval
%                (n.1) and fder (n.1) are analytically computed values and
%                derivatives of f
%    x         : n.1 initial guess for root
%    varargin  : optional parameters passed to f
%  Output
%    x         : n.1 root of f
%    fval      : n.1 function value estimate
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (sqrt(eps))
%    maxsteps  : maximum number of backsteps (25)
%    showiters : display results of each iteration (1)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,fval] = newtond(f,x,varargin)

maxit     = optget('newtond','maxit',100);
tol       = optget('newtond','tol',sqrt(eps));
maxsteps  = optget('newtond','maxsteps',25);
showiters = optget('newtond','showiters',1);

if showiters, fprintf('\nIn newton:\n'), end

for it=1:maxit
  [fval,fder] = feval(f,x,varargin{:});
  fnorm = norm(fval,inf);
  if fnorm<tol
    if showiters, fprintf('\n'), end
    return
  end
  dx = real(-(fval./fder));
  fnormold = inf;
  for backstep=1:maxsteps
    fvalnew = feval(f,x+dx,varargin{:});
    fnormnew = norm(fvalnew,inf);
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew, dx=2*dx; break, end
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = x+dx;
  if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnormnew]); end
end
warning('Failure to converge in newtond')