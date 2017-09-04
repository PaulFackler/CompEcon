%% SOCSOLVE
%
%  Solves infinite horizon continuous time stochastic control model
%
%  More specifically, uses method of collocation to solve the
%  Hamilton-Jacobi-Bellman equation
%   rho*V = max_x f(s,x) + V'(s)*mu(s,x) + 0.5trace(sigma(s,x)*sigma(s,x)'*V''(s))
%  using policy function iteration (i.e., Newton's method). Use of routine
%  is fully documented in the pdf file `SOCSOLVE - A MATLAB Routine for
%  Solving Continuous Time Hamiltion-Jacobi-Bellman Equations' included in
%  the directory CompEcon2015 accompanying the distribution of this file.
%  A summary of its features is provided here:
%
%  Usage
%    [c,sr,vr,xr,resid] = socsolve(model,basis,v)
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
%    c     : nb.1 value function approximant basis function coefficients
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
%      Vss       : whether function file requires second derivative of
%                  value function, true only if diffussion depends on
%                  control (0)
%  Function File
%    User-supplied function that returns the optimal control, reward,
%    drift, and diffusion at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,Vss,<params>)
%    Function File Input
%      flag    : flag that specifies function to be evaluated
%      s       : ns.ds state nodes
%      x       : ns.dx controls
%      Vs      : ns.ds first derivatives of value function
%      Vss     : ns.ds.ds second derivatives of value function
%      params  : user-supplied list of function parameters
%    Function File Output
%      flag = 'x' returns 
%        out   : ns.dx optimal controls at nodes
%      flag = 'f' returns
%        out   : ns.1 reward function values at nodes
%      flag = 'mu' returns
%        out   : ns.ds state drift function values at nodes
%      flag = 'sigma' returns 
%        out   : ns.ds.ds state diffusion function values at nodes
%  Options
%    maxit     : maximum iterations, collocation equation algorithm (500)
%    tol       : convergence tolerance (square root of machine epsilon)
%    nr        : continuous state grid refinement factor (10)
%    output    : whether to print output (1)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,sr,vr,xr,resid] = socsolve(model,basis,v)

% Set user options to defaults, if not set by user with OPTSET (see above)
maxit  = optget('socsolve','maxit',500);  
tol    = optget('socsolve','tol',sqrt(eps));
nr     = optget('socsolve','nr',10);  
output = optget('socsolve','output',1);  
if nargout<5, nr = 0; end

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),   error('Missing Function File');  end  
if ~isfield(model,'params'), error('Missing Parameter List'); end  
if ~isfield(model,'rho'),    error('Missing Discount Rate');  end  
if ~isfield(model,'Vss'),    model.Vss = 0;                   end 

% Unpack model structure
func   = model.func;
params = model.params;
rho    = model.rho;
Vss    = model.Vss;

if output, fprintf('\nIn socsolve:\n'), end

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
Phi2 = cell(ds,ds);
for is=1:ds
  for js=is:ds
    order = zeros(1,ds);
    order(is) = order(is)+1;
    order(js) = order(js)+1;
    Phi2{is,js} = funbase(basis,s,order);
    Phi2{js,is} = funbase(basis,s,order);
  end
end

  tic

% Initialize coefficient vector
if nargin<3 || isempty(v)
  c = zeros(nb,1);
else
  c = Phi0\v;
end

% Policy iteration
for it=1:maxit
  cold = c;
  Vs = zeros(nb,ds);
  for i=1:ds
    Vs(:,i) = Phi1{i}*c;
  end
  if Vss
    Vss = zeros(nb,ds,ds);
    for is=1:ds
      for js=1:ds
        Vss(:,is,js) = Phi2{is,js}*c;
      end
    end
  else
    Vss = [];
  end
  x     = feval(func,'x',s,[],Vs,Vss,params{:});
  f     = feval(func,'f',s,x,[],[],params{:});
  mu    = feval(func,'mu',s,x,[],[],params{:});
  sigma = feval(func,'sigma',s,x,[],[],params{:});
  B = rho*Phi0;
  for is=1:ds
    B = B - diag(mu(:,is))*Phi1{is};
  end
  for is=1:ds
    for js=1:ds
      B = B - 0.5*diag(sigma(:,is,js).*sigma(:,js,is))*Phi2{is,js};
    end
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
  if it==maxit, display('Failure to converge in socsolve.'), end
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
  ns = size(sr,1);
  B  = funbasex(basis,sr,[0;1;2]);
  v  = funeval(c,basis,B);
  Vs = squeeze(funeval(c,basis,B,eye(ds)));
  if Vss
    Vss = zeros(ns,ds,ds);
    for is=1:ds
      for js=1:ds
        order = zeros(1,ds);
        order(is) = order(is)+1;
        order(js) = order(js)+1;
        Vss(:,is,js) = funeval(c,basis,B,order);
      end
    end
  else
    Vss = [];
  end
  x     = feval(func,'x',sr,[],Vs,Vss,params{:});
  f     = feval(func,'f',sr,x,[],[],params{:});
  mu    = feval(func,'mu',sr,x,[],[],params{:});
  sigma = feval(func,'sigma',sr,x,[],[],params{:});
  if nargout>4
    resid = rho*v-f-sum(Vs.*mu,2);
    for is=1:ds
      for js=is:ds
        order = zeros(1,ds);
        order(is) = order(is)+1;
        order(js) = order(js)+1;
        if js==is
          resid = resid - 0.5*diag(sigma(:,is,js).*sigma(:,js,is))*funeval(c,basis,B,order);
        else
          resid = resid - diag(sigma(:,is,js).*sigma(:,js,is))*funeval(c,basis,B,order);
        end
      end
    end
    resid = reshape(resid,[n 1]);
  end
  vr = reshape(v,[n 1]);
  xr = reshape(x,[n size(x,2) 1]);
else
  sr = s;
  vr = v;
  xr = x;
end