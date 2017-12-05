% BROYDEN Computes root of function via Broyden's Inverse Method
% Uses inverse Jacobian estimate and backstepping
% USAGE
%   [x,fval,flag,it,fjacinv] = broyden(f,x,varargin)
% INPUTS
%   f       : name of function of form:
%               fval=f(x,optional additional parameters)
%   x       : initial guess for root (d by 1)
%   varargin: additional arguments for f [optional]
% OUTPUTS
%   x       : root of f (d by 1)
%   fval    : function value estimate (d by 1)
%   flag    : 0 if sucessful
%             1 if failed to converge
%             2 is inf or nan generated
%   it      : number of iterations performed
%   fjacinv : inverse Jacobian estimate (d by d)
%
% Setable options (use OPTSET):
%   maxit     : maximum number of iterations
%   tol       : convergence tolerance
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration
%   initb     : an initial inverse Jacobian approximation matrix
%   initi     : if initb is empty, use the identity matrix to initialize
%               if 0, a numerical Jacobian will be used
%   restarts  : set to 0 to turn off the restarting feature (not recommended)
%   usejac    : set to 1 if f computes a Jacobian (avoids numerical derivatives)

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval,flag,it,fjacinv] = broyden(f,x,varargin)
maxit     = optget('broyden','maxit',100);
tol       = optget('broyden','tol',sqrt(eps));
maxsteps  = optget('broyden','maxsteps',25);
showiters = optget('broyden','showiters',0);
initb     = optget('broyden','initb',[]);
initi     = optget('broyden','initi',0);
restarts  = optget('broyden','restarts',1);
usejac    = optget('broyden','usejac',0);  

flag=0;
if maxsteps<1, maxsteps=1; end
lineoptions.maxsteps=maxsteps;

if usejac && nargout(f)<2
  usejac=false;
end

if isempty(initb)
  if usejac
    [fval,fjacinv] = feval(f,x,varargin{:});
    fjacinv = inv(fjacinv);
  elseif initi
    fjacinv=eye(size(x,1));
  else
    fjacinv = fdjac(f,x,varargin{:});
    fjacinv = inv(fjacinv);
  end
else
  fjacinv = full(initb);
end

fval = feval(f,x,varargin{:});
fnorm = norm(fval,inf);
for it=1:maxit
   if fnorm<tol 
     return; 
   end
   dx = -(fjacinv*fval);
   if any(isnan(dx)|isinf(dx))
     if nargout<3
       warning('broyden:infs','Infinities or NaNs encountered in broyden')
       return
     else
       flag=2; return;
     end
   end
   
   % perform linesearch
  [dx,fvalnew,fnormnew,lineflag,backstep]= ...
       linesearch(f,x,dx,fnorm,lineoptions,varargin{:});
   x = x+dx;
  
   if lineflag>0 && restarts
     if usejac
       [fval,fjacinv] = feval(f,x,varargin{:});
       fjacinv = inv(fjacinv);
     elseif initi
       fjacinv=eye(size(x,1));
     else
       fjacinv = fdjac(f,x,varargin{:});
       fjacinv = inv(fjacinv);
     end
   else
      temp = fjacinv*(fvalnew-fval);
      fjacinv = fjacinv + (dx-temp)*((dx'*fjacinv)/(dx'*temp));
   end
   fval=fvalnew;
   fnorm=fnormnew;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
if nargout<3
  warning('broyden:nonconv','Failure to converge in broyden');
else
  flag=1;
end