%% DOCSIMUL
%
%  Simulates infinite horizon continuous time deterministic optimal control
%  model solved by socsolve.
%
%  Usage
%    [t,ssim,xsim] = docsimul(model,basis,s,x,sinit,T,N)
%  Let
%    ds = dimension of state variable s
%    dx = dimension of control variable x
%    nb = number of basis functions and collocation nodes (prod(n))
%    ns = number of state nodes on output
%    nr = number of replications simulated
%  Input
%    c     : nb.1 value function approximant basis function coefficients
%    model : structured array containing model specification (see scsolve)
%    basis : ds-dimensional basis defined on state space
%    sinit : nr.ds initial states
%    T     : simulation horizon
%    N     : number of time intervals (optional, default is 1000)
%  Output
%    ssim  : nr.N+1.ds simulated states
%    xsim  : nr.N+1.dx simulated controls
%    t     : N+1.1 time nodes
%    if nr=1, arrays are squeezed along singleton dimension

%  Copyright(c) 1997-2015
%    Mario J. Miranda - miranda.4@osu.edu

function [t,ssim,xsim] = docsimul(model,basis,s,x,sinit,T,N)

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),   error('Missing Function File');   end  
if ~isfield(model,'params'), error('Missing Parameter List');  end  

% Unpack model structure
func   = model.func;
params = model.params;

if nargin<7 || isempty(N)
  N = 1000;
end

ds = size(s,2);
if ds==1
  xopt = @(si) interp1(s,x,si);
  g = @(s) feval(func,'g',s,xopt(s),[],params{:});
  [ssim,t] = oderk4(g,sinit,T,N);
  xsim = xopt(ssim);
else
  c = funbase(basis,s)\x;
  g = @(s) feval(func,'g',s,funeval(c,basis,s),[],params{:});
  [ssim,t] = oderk4(g,sinit,T,N);
  xsim = funeval(c,basis,ssim);
end
ssim = squeeze(ssim);
xsim = squeeze(xsim);