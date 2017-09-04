%% BROYDEND
%
%  Computes root of function from R^n to R^n when f_i depends exclusively 
%  on x_i, using Broyden's inverse update method with backstepping.
%
%  Usage
%    [x,fval,fderinv] = broydend(f,x,varargin)
%  Input
%    f         : function of form fval=f(x,varargin) where fval (n.1) are
%                analytically computed values of f
%    x         : n.1 initial guess for root
%    varargin  : optional parameters passed to f
%  Output
%    x         : n.1 root of f
%    fval      : n.1 function value estimate
%    fderinv   : n.1 inverse derivative estimate
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (1e-10)
%    maxsteps  : maximum number of backsteps (25)
%    showiters : display results of each iteration (1)
%    initb     : an initial inverse derivitive approximation vector ([])
%    initi     : if initb empty, and initi is 1, use unit vector to
%                initialize if initi is 0, use numerical Jacobian (0)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,fval,fderinv] = broydend(f,x,varargin)

% Set option defaults, if not set by user with OPTSET
maxit     = optget('broydend','maxit',100);
tol       = optget('broydend','tol',1e-10);
maxsteps  = optget('broydend','maxsteps',25);
showiters = optget('broydend','showiters',1);
initb     = optget('broydend','initb',[]);
initi     = optget('broydend','initi',0);

if showiters, fprintf('\nIn broydend:\n'), end

if maxsteps<1, maxsteps=1; end

if isempty(initb)
  if initi
    fderinv = -eye(size(x,1));
  else
    fderinv = 1./fdjacd(f,x,varargin{:});
  end
else
  fderinv = initb;
end

fval = feval(f,x,varargin{:});
fnorm = norm(fval,'inf');
for it=1:maxit
  if fnorm<tol
    if showiters, fprintf('\n'), end
    return
  end
  dx = -(fderinv.*fval);
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
      fderinv = eye(size(x,1));
    else
      fderinv = 1./fdjacd(f,x,varargin{:});
    end
  else
    temp = fderinv.*(fvalnew-fval);
    fderinv = fderinv + (dx-temp).*(dx'*fderinv/(dx'*temp));
  end
  fval = fvalnew;
  fnorm = fnormnew;
  if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
warning('Failure to converge in broydend')