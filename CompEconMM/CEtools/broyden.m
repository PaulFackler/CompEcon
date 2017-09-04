%% BROYDEN
%
%  Computes root of function from R^n to R^n using Broyden's inverse update
%  method with backstepping.
%
%  Usage
%    [x,fval,fjacinv] = broyden(f,x,varargin)
%  Input
%    f         : function of form fval=f(x,varargin)
%    x         : n.1 initial guess for root
%    varargin  : optional parameters passed to f
%  Output
%    x         : n.1 root of f
%    fval      : n.1 function value estimate
%    fjacinv   : n.n inverse Jacobian estimate
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (1e-10)
%    maxsteps  : maximum number of backsteps (25)
%    showiters : display results of each iteration (1)
%    initb     : an initial inverse Jacobian approximation matrix ([])
%    initi     : if initb empty, and initi is 1, use identity matrix to
%                initialize if initi is 0, use numerical Jacobian (0)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,fval,fjacinv] = broyden(f,x,varargin)

% Set option defaults, if not set by user with OPTSET
maxit     = optget('broyden','maxit',100);
tol       = optget('broyden','tol',1e-10);
maxsteps  = optget('broyden','maxsteps',25);
showiters = optget('broyden','showiters',1);
initb     = optget('broyden','initb',[]);
initi     = optget('broyden','initi',0);

if showiters, fprintf('\nIn broyden:\n'), end

if maxsteps<1, maxsteps=1; end

if isempty(initb)
  if initi
    fjacinv = -eye(size(x,1));
  else
    fjacinv = inv(fdjac(f,x,varargin{:}));
  end
else
  fjacinv = initb;
end

fval = feval(f,x,varargin{:});
fnorm = norm(fval,'inf');
for it=1:maxit
  if fnorm<tol
    if showiters, fprintf('\n'), end
    return
  end
  dx = -(fjacinv*fval);
  fnormold = inf;
  for backstep=1:maxsteps
    fvalnew = feval(f,x+dx,varargin{:});
    fnormnew = norm(fvalnew,'inf');
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew
      fvalnew = fvalold;
      fnormnew = fnormold;
      dx = dx*2;
      break
    end
    fvalold  = fvalnew;
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = real(x+dx);
  if any(isnan(x)|isinf(x))
    error('Infinities or NaNs encountered.')
  end
  if fnormnew>fnorm
    if initi
      fjacinv = eye(size(x,1));
    else
      fjacinv = fdjac(f,x,varargin{:});
      fjacinv = inv(fjacinv);
    end
  else
    temp = fjacinv*(fvalnew-fval);
    fjacinv = fjacinv + (dx-temp)*(dx'*fjacinv/(dx'*temp));
  end
  fval = fvalnew;
  fnorm = fnormnew;
  if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
warning('Failure to converge in broyden')