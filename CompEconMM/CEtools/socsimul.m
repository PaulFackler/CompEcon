%% SOCSIMUL
%
%  Simulates infinite horizon continuous time stochastic optimal control
%  model solved by socsolve.
%
%  Usage
%    [t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep,N)
%  Let
%    ds = dimension of state variable s
%    dx = dimension of control variable x
%    nb = number of basis functions and collocation nodes (prod(n))
%    ns = number of state nodes on output
%  Input
%    c     : nb.1 value function approximant basis function coefficients
%    model : structured array containing model specification (see scsolve)
%    basis : ds-dimensional basis defined on state space
%    sinit : 1.ds initial state
%    T     : simulation horizon
%    nrep  : number of replications simulated (optional, default is 1)
%    N     : number of time intervals (optional, default is 1000)
%  Output
%    t     : N+1.1 time nodes
%    ssim  : nr.N+1.ds simulated states, where nr=min(3,nrep)
%    xsim  : nr.N+1.dx simulated controls, where nr=min(3,nrep)
%    smean : N+1.ds mean of simulated states across replications
%    xmean : N+1.dx mean of simulated controls across replications
%    if nrep=1, arrays are squeezed along singleton dimension

%  Copyright(c) 1997-2015
%    Mario J. Miranda - miranda.4@osu.edu

function [t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep,N)

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),   error('Missing Function File');   end  
if ~isfield(model,'params'), error('Missing Parameter List');  end  
if ~isfield(model,'dx'),     model.dx = 1;                     end  
if ~isfield(model,'Vss'),    model.Vss = 0;                    end 

% Unpack model structure
func   = model.func;
params = model.params;
dx     = model.dx;
Vss    = model.Vss;

% Time discretization
if nargin<6 || isempty(nrep)
  nrep = 1;
end
if nargin<7 || isempty(N)
  N = 1000;
end
h = T/N;
t = (0:h:T);

% Initialize output arrays
ds   = size(sinit,2);
ssim = zeros(nrep,N+1,ds);
xsim = zeros(nrep,N+1,dx);

% Simulate Model
z  = randn(nrep,N,ds);
ss = sinit(ones(nrep,1),:);
for i=0:N
  Vs = funeval(c,basis,ss,eye(ds));
  if Vss
    Vss = zeros(nrep,ds,ds);
    for is=1:ds
      for js=1:ds
        order = zeros(1,ds);
        order(is) = order(is)+1;
        order(js) = order(js)+1;
        Vss(:,is,js) = funeval(c,basis,ss,order);
      end
    end
  else
    Vss = [];
  end
  xx    = feval(func,'x',ss,[],Vs,Vss,params{:});
  mu    = feval(func,'mu',ss,xx,[],[],params{:});
  sigma = feval(func,'sigma',ss,xx,[],[],params{:});
  ssim(:,i+1,:) = ss;
  xsim(:,i+1,:) = xx;
  if i<N
    ss  = ss + mu*h;
    zz  = squeeze(z(:,i+1,:));
    for is1=1:ds
      for is2=1:ds
        ss(:,is1) = ss(:,is1) + sigma(:,is1,is2)*sqrt(h).*zz(:,is2);
      end
    end
  end
end
ssim = squeeze(ssim);
xsim = squeeze(xsim);
if nargout>3, smean=squeeze(mean(ssim,1)); end
if nargout>4, xmean=squeeze(mean(xsim,1)); end
if nrep>3
  ssim = ssim(1:3,:,:);
  xsim = xsim(1:3,:);
end