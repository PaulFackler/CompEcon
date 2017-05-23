% OPTSTEP Solves a one dimensional optimal step length problem
% Used to compute step lengths in multidimensional optimization
% USAGE
%   [s,fx,errcode,iters]=optstep(method,f,x0,f0,g0,d,MaxIters,varargin);
% INPUTS
%   method  : one of 4 methods
%             1 - step length is set to 1
%             2 - STEPBHHH
%             3 - STEPBT (default)
%             4 - STEPGOLD (called others fail)
%   f       : a function file name
%   x0      : the starting value
%   f0      : f(x0)
%   g0      : gradient of f at x: f'(x0)
%   d       : the search direction
%   MaxIters: the maximum # of itereations before trying something else
% OUTPUTS
%   s       : the optimal step in the direction d
%   fx      : the value of f at x+s*d,
%   iters   : the number of iterations used
%   errcode : 0 if STEP suitable step found

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [s,fx,errcode,iters]=optstep(method,f,x0,f0,g0,d,MaxIters,varargin)

if nargin<6, method=3; end
if nargin<7, MaxIters=100; end

switch method
case 1
  fx=feval(f,x0+d,varargin{:});
  if fx<f0
    s=1; iters=1; errcode=0;
  else
    [s,fx,iters,errcode]=stepgold(f,x0,f0,g0,d,MaxIters,varargin{:});
  end
case 2  
  [s,fx,iters,errcode]=stepbhhh(f,x0,f0,g0,d,MaxIters,varargin{:});
  stemp=[];
  if errcode
    stemp=s;
    [s,fx,iters2,errcode]=stepgold(f,x0,f0,g0,d,MaxIters,varargin{:});
    iters=iters+iters2;
  end
case 3  
  [s,fx,iters,errcode]=stepbt(f,x0,f0,g0,d,MaxIters,varargin{:});
  if errcode
    [s,fx,iters2,errcode]=stepgold(f,x0,f0,g0,d,MaxIters,varargin{:});
    iters=iters+iters2;
  end
otherwise 
  [s,fx,iters,errcode]=stepgold(f,x0,f0,g0,d,MaxIters,varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,fs,iter,errcode]=stepbhhh(f,x0,f0,g0,d,MaxIters,varargin);

% STEPBHHH Computes an approximate minimum step length
%
% USAGE:
%  [s,fs,iters,errcode]=STEPBHHH(f,x0,f0,g0,d,MaxIters);
%
% INPUTS:
%        f : the objective function being minimized
%       x0 : the current value of the parameters
%       f0 : the value of f(x0)
%       g0 : the gradient vector of f at x0
%        d : the search direction       
% MaxIters : the maximum number of function evaluations allowed
%
% OUTPUTS:
%    s:  the optimal step in the direction d
%    fs: the value of f at x+s*d,
%    iter:  the number of iterations used
%    errcode: equals 1 if maximum iterations are exceeded

% STEPBHHH uses an algorithm based on one discussed in Berndt, et. al.,
% Annals of Economic and Social Measurement, 1974, pp. 653-665.
% This procedure specifies a cone of convergence in the plane
% defined by the direction vector, d, and the value of the objective
% function. The cone is defined by the lines through the origin
% (x,f(x)) with slopes (d'g)*delta and (d'g)*(1-delta).
% Delta must lie on (0,0.5).
% The procedure iterates until a point is found on the objective function
% that lies within the cone.
% In general, the wider the cone, the faster a "suitable" step size
% will be found.  If a trial point lies above the cone
% the step size will be increased and if it lies below the cone
% the step size is decreased. 

% ---------- INITIALIZATIONS ---------------------
  if nargin<6 | isempty(MaxIters), MaxIters=25; end
  delta=optget('optstep','bhhhcone',0.0001);
  dg=-d'*g0;                   % directional derivative 
  tol1=dg*delta;
  tol0=dg*(1-delta);
  s=1;
  ds=1;
  errcode=0;
  
  % first bracket the cone
  for iter=1:MaxIters
    x=x0+s*d; fs=feval(f,x,varargin{:});
    temp=(f0-fs)/s;
    if temp<tol0
      ds=2*ds;
      s=s+ds;
    else break
    end
  end
  
  if tol0<=temp & temp<=tol1, return, end
  
  ds=ds/2;
  s=s-ds;
  it=iter+1;   
  
% then use bisection to get inside it
  for iter=it:MaxIters
    ds = ds/2;
    x=x0+s*d; fs=feval(f,x,varargin{:}); 
    temp=(f0-fs)/s;
    if     temp > tol1, s=s-ds; 
    elseif temp < tol0, s=s+ds;
    else return
    end  
  end
 
  errcode=1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,fs,iter,errcode]=stepbt(f,x0,f0,g0,d,MaxIters,varargin);

% STEPBT Computes an approximate minimum step length
%
% USAGE:
%  [s,fs,iters,errcode]=stepbt(f,x0,f0,g0,d,MaxIters);
%
% INPUTS:
%    f : the objective function being minimized,
%   x0 : the current value of the parameters
%   f0 : the value of f(x); this is passed as argument to save
%                  one function evaluation.
%   g0 : the gradient vector of f at x0
%    d : the search direction       
%   MaxIters -- the maximum number of "backsteps" of the step length
%
% OUTPUTS:
%    s:  the optimal step in the direction d
%    fs: the value of f at x+s*d,
%    iter:  the number of iterations used
%    errcode: equals 1 if STEPBT fails to find a suitable step length
%                    2 if cubic approximation finds negative root

% STEPBT uses a backtracking method similar to Algorithm 6.3.5 in
% Dennis and Schnabel, Numerical Methods for Unconstrained Optimization
% and Nonlinear Equations or LNSRCH in sec 9.7 of Press, et al., 
% Numerical Recipes. The algorithm approximates the function with
% a cubic using the function value and derivative at the initial point
% and two additional points.  It determines the minimum of the 
% approximation.  If this is acceptable it returns, otherwise it uses the
% current and precious point to form a new approximation.
% The convergence criteria is similar to that discussed in Berndt, 
% et. al., Annals of Economic and Social Measurement, 1974, pp. 653-665
%(see description of BHHHSTEP). The change in the step size is also 
% limited to ensure that 
%   lb*s(k)<=s(k+1)<=ub*s(k) (defaults: lb=0.1, ub=0.5).


% --------------------- Initializations ------------------------- 
  delta=1e-4;    % Defines cone of convergence; must be on (0,1/2) 
  ub=0.5;        % Upper bound on acceptable reduction in s.   
  lb=0.1;        % Lower bound on acceptable reduction in s.  
  
  errcode=0;
  dg=-d'*g0;                      % directional derivative
  
  tol1=delta*dg;
  tol0=(1-delta)*dg;

% full step
  s=1;
  fs=feval(f,x0+d,varargin{:});
  if -fs+f0 <= tol1, iter=1; return, end
  
% quadratic approximation 
  s2=s; fs2=fs;
  s=-0.5*dg./(-fs+f0-dg);
  s = max(s,lb);
  fs=feval(f,x0+s*d,varargin{:});
  temp=(-fs+f0)./s;
  if tol0<=temp & temp <= tol1, iter=2; return, end

% cubic approximation 
  for iter=3:MaxIters
    temp=(s-s2)*[s*s;s2*s2];
    temp=[(-fs+f0-dg*s);(-fs2+f0-dg*s2)]./temp;
    a=temp(1)-temp(2);
    b=s*temp(2)-s2*temp(1);
    s2=s; fs2=fs;
    if a == 0                              % quadratic fits exactly
      s=-0.5*dg/b;
    else
      disc = b*b - 3*a*dg;
      if disc < 0, errcode=2; return, end  % complex root 
      s=(sqrt(disc)-b)/(3*a);
    end
    s=max(min(s,ub*s2),lb*s2);             % ensures acceptable step size
    fs=feval(f,x0+s*d,varargin{:});
    temp=(-fs+f0)./s;
    if tol0<=temp & temp <= tol1, return, end
 end
 
 errcode=1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,fs,iter,errcode]=stepgold(f,x0,f0,g0,d,MaxIters,varargin)

% STEPGOLD Computes an approximate minimum step length
%
% USAGE:
%  [s,fs,iters,errcode]=STEPGOLD(f,x0,f0,g0,d,MaxIters);
%
% INPUTS:
%        f : the objective function being minimized
%       x0 : the current value of the parameters
%       f0 : the value of f(x0)
%       g0 : the gradient vector of f at x0 (note: not used)
%        d : the search direction       
% MaxIters : the maximum number of function evaluations allowed
%
% OUTPUTS:
%    s:  the optimal step in the direction d
%    fs: the value of f at x+s*d,
%    iter:  the number of iterations used
%    errcode: equals 1 if maximum iterations are exceeded

% STEPGOLD uses step doubling to find an initial bracket and then
% uses the golden search method to find a minimum value within
% the bracket.  Iterations cease if the bracket is less than TOL
% or a maximum number of iterations is reached.

alpha1=(3-sqrt(5))/2; alpha2=(sqrt(5)-1)/2;

tol=1e-4;                % tolerance used for Golden search algorithm 
tol=tol*(alpha1*alpha2); % the bracket will be len/(alpha1*alpha2)
s=1;
errcode=1;               % 1 if the search is unsuccessful; otherwise 0
iter=0; 
s0=0;
%  Find a bracketing interval
fs=feval(f,x0+d,varargin{:});
if f0 >= fs, len=alpha1;              
else 
  for iter=1:MaxIters
    s=2*s;
    fl=fs;
    fs=feval(f,x0+s*d,varargin{:});
    if fs <= fl; len=alpha1*(s-s0); break
    else; f0=fl; s0=s/2;
    end
  end
  if iter>=MaxIters, s=s/2; fs=fl; return; end
end

xl=x0+(s0+len)*d;
xs=x0+(s-len)*d;

s=s-len;
len=len*alpha2;  % len now measures relative distance between xl and xs

fs=feval(f,xs,varargin{:});
fl=feval(f,xl,varargin{:});
% Golden search to find minimum
while iter < MaxIters
  iter=iter+1;
  if fs<fl
    s=s-len;
    len=len*alpha2;
    xs=xl; xl=xl-len*d;
    fs=fl; fl=feval(f,xl,varargin{:});
  else 
    len=len*alpha2;
    s=s+len;
    xl=xs; xs=xs+len*d;
    fl=fs; fs=feval(f,xs,varargin{:});
  end
  if len<tol, errcode=0; break, end
end

if fl>fs, fs=fl; s=s-len; end
