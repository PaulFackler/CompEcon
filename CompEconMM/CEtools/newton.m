%% NEWTON   
%
%  Computes root of function from R^n to R^n using Newton's method with
%  backstepping.
%
%  Usage
%   [x,fval] = newton(f,x,varargin)
%  Input
%    f         : function of form [fval,fjac]=f(x,varargin) where fval 
%                (n.1) and fjac (n.n) are analytically computed values and
%                Jacobian of f
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

function [x,fval] = newton(f,x,varargin)

maxit     = optget('newton','maxit',100);
tol       = optget('newton','tol',sqrt(eps));
maxsteps  = optget('newton','maxsteps',25);
showiters = optget('newton','showiters',1);

if showiters, fprintf('\nIn newton:\n'), end

checkjac(f,x,varargin{:});

for it=1:maxit
  [fval,fjac] = feval(f,x,varargin{:});
  fnorm = norm(fval,inf);
  if fnorm<tol
    if showiters, fprintf('\n'), end
    return
  end
  dx = real(-(fjac\fval));
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
warning('Failure to converge in newton')