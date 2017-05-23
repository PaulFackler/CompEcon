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
%   flag     : 0 - success, 1 - failure
%   it       : number of iterations used to obtain new value of dx
%
% The options structure may contain the following fields
%   maxstep  : maximum number of iterations [20]
%
% Used by newton, broyden and icum
function [dx,fvalnew,fnormnew,flag,it]=linesearch(f,x,dx,fnorm,options,varargin)
if 1
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
else
  merit=@(lambda) norm(feval(f,x+lambda*dx,varargin{:}));
  [lambda,fnormnew,flag,output]=fminbnd(merit,0,1);
  it=output.iterations;
  dx=lambda*dx;
  fvalnew=feval(f,x+dx,varargin{:});
  return
end