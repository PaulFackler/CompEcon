% NEWTON   Computes root of function via Newton's Method with backstepping
% USAGE
%   [x,fval,flag,it] = newton(f,x,varargin);
% INPUTS
%   f        : name of function of form:
%               [fval,fjac]=f(x,optional additional parameters)
%   x        : initial guess for root
%   varargin : additional arguments for f (optional)
% OUTPUTS
%   x        : root of f
%   fval     : function value estimate
%   flag     : 0 if sucessful
%              1 if failed to converge
%              2 is inf or nan generated
%   it       : number of iterations performed
%   it       : number of iterations performed
%
% Setable options (use OPTSET):
%   tol       : convergence tolerance
%   maxit     : maximum number of iterations
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval,flag,it] = newton(f,x,varargin)

maxit     = optget('newton','maxit',100);
tol       = optget('newton','tol',sqrt(eps));
maxsteps  = optget('newton','maxsteps',25);
showiters = optget('newton','showiters',0);

flag=0;
if maxsteps<1, maxsteps=1; end
lineoptions.maxsteps=maxsteps;
for it=1:maxit
   [fval,fjac] = feval(f,x,varargin{:});
   fnorm = norm(fval,inf);
   if fnorm<tol, 
     return, 
   end
   if 0
     [L1,U1] = ilu(fjac);
     [dx,flag,relres,iter]=bicgstab(fjac,-fval,[],[],L1,U1);
   else
     dx = -(fjac\fval);
   end
   if any(isnan(dx)|isinf(dx))
     if nargout<3
       warning('Infinities or NaNs encountered.')
       return
     else
       flag=2; return;
     end
   end
   
   % perform linesearch
   [dx,fvalnew,fnormnew,lineflag,backstep]= ...
       linesearch(f,x,dx,fnorm,lineoptions,varargin{:});
      
   x = x+dx;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnormnew]); end
end
if nargout>2
  flag=1;
else
  warning('Failure to converge in newton');
end


% LINESEARCH Conducts a linesearch for a suitable step in rootfinding problems
% USAGE 
%   [dx,fvalnew,fnormnew,flag,it]=linesearch(f,x,dx,fnorm,options,varargin);
% INPUTS
%   f        :  function mapping R^n to R^n
%   x        :  nx1 vector: current evaluation point
%   dx       : nx1 vector: search direction
%   fnorm    : ||f(x)||
%   options  : structure specifying options (described below)
%   varargin : additional parameters to pass to f
% OUTPUTS
%   dx       : updated search direction
%   fvalnew  : f(x+dx)
%   fnormnew : ||f(x+dx)||
%   flag     : 0 - success, 1 - failure, 2 - acceptable
%   it       : number of iterations used to obtain new value of dx
%
% The options structure may contain the following fields
%   maxstep  : maximum number of iterations [20]
%
% Used by newton, broyden and icum
function [dx,fvalnew,fnormnew,flag,it]=linesearch(f,x,dx,fnorm,options,varargin)

  if nargin<5 | isempty(options)
    maxstep=20;
  else
    maxsteps=options.maxsteps;
  end
  fnormold = inf;
  for it=1:maxsteps
     fvalnew = feval(f,x+dx,varargin{:});
     fnormnew = norm(fvalnew,inf);
     if fnormnew<fnorm, flag=0; return, end
     if fnormold<fnormnew
       fvalnew=fvalold;
       fnormnew=fnormold; 
       dx=dx*2; 
       flag=2;
       return
     end
     fvalold  = fvalnew;
     fnormold = fnormnew;
     dx = dx/2;
  end
  flag=1;