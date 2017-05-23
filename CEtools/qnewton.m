% QNEWTON Solves unconstrained maximization problem using quasi-Newton
% USAGE
%   [x,A]=qnewton(f,x,A,varargin);
% INPUTS
%   f        : name of function of form [fval,fjac]=f(x,varargin)
%   x        : initial guess for maximum
%   A        : initial value for inverse Hessian approximation (optional)
%   varargin : additional arguments for f (optional)
% OUTPUTS
%   x        : local maximum of f
%   A        : inverse Hessian approximation are final iteration
%   flag     : convergence flag
%                0 : success
%                1 : NaNs or INFs encountered
%                2 : maxit exceeded
%                3 : iterations are stuck
%                4 : suitable step direction can't be found
% 
% The user defined function f must have the following syntax
%   [fx,g] = f(x,additional variables)
% Optionally, if SearchMeth=4 the syntax must be
%   [fx,g,A] = f(x,additional variables)
% where, in either case, the additional variables are the ones
% passed to QNEWTON
%
% Setable options (use OPTSET):
%   SearchMeth : 1)Steepest Ascent 2)DFP 3)BFGS 4) User defined [3]
%   StepMeth   : 1)No search 2)BHHHSTEP 3)STEPBT 4)GOLDSTEP [3]
%   maxit      : Maximum major iterations [250]
%   maxstep    : Maximum step search iterations [50]
%   tol        : convergence tolerence [sqrt(eps)]
%   eps0       : zero factor (used in convergence criteria) [1]
%   ShowIters  : 1 to show results at each iteration [0]

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,A,flag] = qnewton(f,x,A,varargin)

SearchMeth   = optget('qnewton','SearchMeth',3);
StepMeth     = optget('qnewton','StepMeth',3);
maxit        = optget('qnewton','maxit',250);
maxstep      = optget('qnewton','maxstep',50);
tol          = optget('qnewton','tol',sqrt(eps));
eps0         = optget('qnewton','eps0',1);
eps1         = optget('qnewton','eps1',1e-12);
ShowIters    = optget('qnewton','ShowIters',0);

k = size(x,1);
reset = 0;

if SearchMeth==4
  [fx0,g0,A] = feval(f,x,varargin{:});
else
  [fx0,g0] = feval(f,x,varargin{:});
end

flag=0;
if abs(g0)<eps, return; end

if nargin<3 | isempty(A)
   A = -eye(k)./max(abs(fx0),1); 
   reset = 1;
end
for it=1:maxit

  d = -(A*g0);                       % search direction

  if ((d'*g0)./(d'*d)) < eps1        % must go uphill
    A = -eye(k)./max(abs(fx0),1);    % otherwise use
    d = g0./max(abs(fx0),1);         % steepest ascent
    reset = 1;
  end

  [s,fx] = optstep(StepMeth,f,x,fx0,g0,d,maxstep,varargin{:});

  if fx<=fx0                         % Step search failure
    if reset
      if nargout<3
        warning('Iterations stuck in qnewton'), return;
      else
        flag=3; return;
      end
    else                             % Use steepest ascent
      A = -eye(k)./max(abs(fx0),1);
      d = g0./max(abs(fx0),1);
      [s,fx,err]= optstep(StepMeth,f,x,fx0,g0,d,maxstep,varargin{:});
      if err,
        if nargout<3
          warning('Cannot find suitable step in qnewton'), return;
        else
          flag=4; return;
        end       
      end
    end
  end 

  d = s*d;
  x = x+d;

  if any(isnan(x)|isinf(x))
    if nargout<3
      warning('NaNs or INFs encountered in qnewton'); return
    else
      flag=1; return
    end
  end


  if SearchMeth==4
    [fx,g,A] = feval(f,x,varargin{:});
  else
    [fx,g] = feval(f,x,varargin{:});
  end
  
  if ShowIters
    fprintf('qnewton: %4i %16.4f %16.4f %16.4f\n',it,fx,norm(d),norm(g)); end

  % Test convergence using Marquardt's criteria and gradient test
  if ((fx-fx0)/(abs(fx)+eps0)<tol & all(abs(d)./(abs(x)+eps0)<tol)) ...
      | all(abs(g)<eps); return; 
  end
  
  % Update Inverse Hessian
  u = g-g0; ud = u'*d;
  if SearchMeth==1 | abs(ud)<eps          % Steepest ascent 
    A = -eye(k)./max(abs(fx),1);
    reset = 1;
  elseif SearchMeth==2;                   % DFP update     
    v = A*u;
    A = A + d*d'./ud - v*v'./(u'*v);
    reset = 0;
  elseif SearchMeth==3;                   % BFGS update
    w = d-A*u; wd = w*d';
    A = A + ((wd + wd') - ((u'*w)*(d*d'))./ud)./ud;
    reset = 0;
  end
  
  %  Update iteration 
  fx0 = fx; g0 = g;

end
flag=2;
if nargout<3
  warning('Maximum iterations exceeded in qnewton')
end
