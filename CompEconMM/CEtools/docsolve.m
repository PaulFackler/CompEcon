%% DOCSOLVE
%
%  Solves infinite horizon continuous time deterministic control model
%
%  More specifically, uses method of collocation to solve the
%  Hamilton-Jacobi-Bellman equation
%   rho*V = max_x f(s,x) + V'(s)*g(s,x)
%  using policy function iteration (i.e., Newton's method). Use of routine
%  is fully documented in the pdf file `DOCSOLVE - A MATLAB Routine for
%  Solving Continuous Time Hamiltion-Jacobi-Bellman Equations' included in
%  the directory CompEcon2015 accompanying the distribution of this file.
%  A summary of its features is provided here:
%
%  Usage
%    [c,sr,vr,xr,resid] = docsolve(model,basis,v)
%  Let
%    ds = dimension of state variable s
%    dx = dimension of control variable x
%    n  = number of basis functions per state dimension (1.ds)
%    nb = number of basis functions and collocation nodes (prod(n))
%    ns = number of state nodes on output
%  Input
%    model : structured array containing model specification (see below)
%    basis : ds-dimensional basis defined on state space
%    v     : nb.1 initial guess for values at collocation nodes (optional,
%            default is array of zeros)
%  Output
%    c     : nb.1 value function approximant basis coefficients
%    sr    : ns.ds refined state node grid
%    vr    : ns.1 values at refined state nodes
%    xr    : ns.dx optimal controls at refined state nodes
%    resid : ns.1 HBJ equation residuals at refined state nodes
%    Note: If the residual is requested, then values, optimal controls, and
%    residuals are returned on a refined grid of ns state nodes.  The
%    degree of refinement is governed by nr, an optional parameter with
%    default value of 10 that may be set by the user using optset.  If
%    nr>0, the refined state node grid is created by forming the Cartesian
%    product of nr*n equally-spaced coordinates along each state dimension.
%    If nr=0, the routine will return the values, optimal controls, and
%    residuals at the original collocation state nodes. On output, the
%    singleton dimensions of c, v, x, and resid are eliminated.
%  Model Structure
%    The structured array "model" contains fields that specify essential
%    components of the model to be solved (default values for optional
%    fields in parentheses):
%      func      : name of function file (required)
%      params    : model parameters required by function file (empty) 
%      rho       : continuous discount rate (required)
%      dx        : dimension of control variable (1)
%  Function File
%    User-supplied function that returns the optimal control, reward, and
%    transition at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,<params>)
%    Function File Input
%      flag    : flag that specifies function to be evaluated
%      s       : ns.ds state nodes
%      x       : ns.dx controls
%      Vs      : ns.ds first derivatives of value function
%      params  : user-supplied list of function parameters
%    Function File Output
%      flag = 'x' returns 
%        out   : ns.dx optimal controls at nodes
%      flag = 'f' returns
%        out   : ns.1 reward function values at nodes
%      flag = 'g' returns
%        out   : ns.ds transition function values at nodes
%  Options
%    maxit     : maximum iterations, collocation equation algorithm (500)
%    tol       : convergence tolerance (square root of machine epsilon)
%    nr        : continuous state grid refinement factor (10)
%    output    : whether to print output (1)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,sr,vr,xr,resid] = docsolve(model,basis,v)

% Set user options to defaults, if not set by user with OPTSET (see above)
maxit  = optget('docsolve','maxit',500);  
tol    = optget('docsolve','tol',sqrt(eps));
nr     = optget('docsolve','nr',10);  
output = optget('docsolve','output',1);  
if nargout<5, nr = 0; end

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),   error('Missing Function File');  end  
if ~isfield(model,'params'), error('Missing Parameter List'); end  
if ~isfield(model,'rho'),    error('Missing Discount Rate');  end  

% Unpack model structure
func   = model.func;
params = model.params;
rho    = model.rho;

if output, fprintf('\nIn docsolve:\n'), end

% Initialize state collocation nodes
ds = basis.d;               % dimension of state variable s
n  = basis.n;               % number of collocation nodes per dimension
nb = prod(n);               % total number of basis functions and collocation nodes
s  = funnode(basis);        % collocation nodes

% Derive interpolation matrices
Phi0 = funbase(basis,s);
Phi1 = cell(ds,1);
for is=1:ds
  order = zeros(1,ds);
  order(is) = 1;
  Phi1{is} = funbase(basis,s,order);
end

tic

% Initialize coefficient vector
if nargin<3 || isempty(v)
  v = zeros(nb,1);
end
c = Phi0\v;

% Policy iteration
for it=1:maxit
  cold = c;
  Vs = zeros(nb,ds);
  for i=1:ds
    Vs(:,i) = Phi1{i}*c;
  end
  x = feval(func,'x',s,[],Vs,params{:});
  f = feval(func,'f',s,x,[],params{:});
  g = feval(func,'g',s,x,[],params{:});
  B = rho*Phi0;
  for is=1:ds
    B = B - diag(g(:,is))*Phi1{is};
  end
  c = B\f;
  if any(isnan(c)|isinf(c))
    error('NaNs or Infs encountered');
  end
  change = max(abs(c-cold));
  if output, fprintf ('%4i %10.1e\n',it,change), end
  if change<tol
    if output, fprintf('\n'), end
    break
  end
end

if output
  if it==maxit, display('Failure to converge in docsolve.'), end
  fprintf('Elapsed Time = %7.2f Seconds\n',toc)
end

% Compute values and controls on refined continuous state grid, if requested.
if nargout>4 && nr>0
  if ds==1,
    n = nr*n+1;
    scoord = linspace(basis.a,basis.b,n)';
  else
    scoord = cell(ds,1);
    for i=1:ds
      n(i) = nr*n(i)+1;
      scoord{i} = linspace(basis.a(i),basis.b(i),n(i))';
    end
  end
  sr = gridmake(scoord);
  B  = funbasex(basis,sr,[0;1]);
  v  = funeval(c,basis,B);
  Vs = squeeze(funeval(c,basis,B,eye(ds)));
  x  = feval(func,'x',sr,[],Vs,params{:});
  f  = feval(func,'f',sr,x,[],params{:});
  g  = feval(func,'g',sr,x,[],params{:});
  if nargout>4
    resid = rho*v-f-sum(Vs.*g,2);
    resid = reshape(resid,[n 1]);
  end
  vr = v;
  xr = x;
else
  sr = s;
  vr = v;
  xr = x;
end

% Eliminate singleton dimensions to facilitate analysis in calling program
vr = squeeze(vr);
xr = squeeze(xr);