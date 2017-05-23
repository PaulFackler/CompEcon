% BROYDENX Computes root of function via Broyden's Method
% Uses direct Jacobian estimate and backstepping
% USAGE
%   [x,fval,fjac] = broydenx(f,x,varargin);
% INPUTS
%   f       : name of function of form:
%               fval=f(x,optional additional parameters)
%   x       : initial guess for root (d by 1)
%   varargin: additional arguments for f [optional]
% OUTPUTS
%   x       : root of f (d by 1)
%   fval    : function value estimate (d by 1)
%   fjac    : Jacobian estimate (d by d)
%
% Setable options (use OPTSET):
%   maxit     : maximum number of iterations
%   tol       : convergence tolerance
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval,fjac] = broydenx(f,x,varargin)

maxit     = optget('broydenx','maxit',100);
tol       = optget('broydenx','tol',sqrt(eps));
maxsteps  = optget('broydenx','maxsteps',25);
showiters = optget('broydenx','showiters',0);
if maxsteps<1, maxsteps=1; end

fjac = fdjac(f,x,varargin{:});

fval = feval(f,x,varargin{:});
fnorm = norm(fval);
for it=1:maxit
   if fnorm<tol, return; end
   dx = -(fjac\fval);
   fnormold = inf;
   for backstep=1:maxsteps
      fvalnew = feval(f,x+dx,varargin{:});
      fnormnew = norm(fvalnew,'inf');
      if fnormnew < fnorm, break, end
      if fnormold < fnormnew
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
   fjac = fjac + (fvalnew-fval-fjac*dx)*(dx./(dx'*dx))';
   fval=fvalnew;
   fnorm=fnormnew;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
warning('broydenx:FailToConverge','Failure to converge in broydenx');