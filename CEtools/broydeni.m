% BROYDENI Computes root of function via Broyden's Inverse Method
% Uses inverse Jacobian estimate and backstepping
% USAGE
%   [x,fval,fjacinv] = broydeni(f,x,varargin)
% INPUTS
%   f       : name of function of form:
%               fval=f(x,optional additional parameters)
%   x       : initial guess for root (d by 1)
%   varargin: additional arguments for f [optional]
% OUTPUTS
%   x       : root of f (d by 1)
%   fval    : function value estimate (d by 1)
%   fjacinv : inverse Jacobian estimate (d by d)
%
% Setable options (use OPTSET):
%   maxit     : maximum number of iterations
%   tol       : convergence tolerance
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval,fjacinv] = broydeni(f,x,fjacinv,varargin)

maxit     = optget('broydeni','maxit',100);
tol       = optget('broydeni','tol',sqrt(eps));
maxsteps  = optget('broydeni','maxsteps',25);
showiters = optget('broydeni','showiters',0);

if maxsteps<1, maxsteps=1; end

if ~exist('fjacinv','var') || isempty(fjacinv)
  fjacinv = inv(fdjac(f,x,varargin{:}));
end

fval = feval(f,x,varargin{:});
fnorm = norm(fval,'inf');
for it=1:maxit
   if fnorm<tol, return; end
   dx = -(fjacinv*fval);
   fnormold = inf;
   for backstep=1:maxsteps
      fvalnew = feval(f,x+dx,varargin{:});
      fnormnew = norm(fvalnew,'inf');
      if fnormnew<fnorm, break, end
      if fnormold<fnormnew
        fvalnew=fvalold;
        fnormnew=fnormold; 
        dx=dx*2; 
        break
      end
      fvalold  = fvalnew;
      fnormold = fnormnew;
      dx = dx/2;
   end
   x = x+dx;
   if any(isnan(x)|isinf(x))
     error('Infinities or NaNs encountered.')
   end

   if fnormnew>fnorm
      fjacinv = inv(fdjac(f,x,varargin{:}));
   else
      temp = fjacinv*(fvalnew-fval);
      fjacinv = fjacinv + (dx-temp)*(dx'*fjacinv/(dx'*temp));
   end
   fval=fvalnew;
   fnorm=fnormnew;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
warning('Broydeni:FailureToConverge','Failure to converge in broydeni');