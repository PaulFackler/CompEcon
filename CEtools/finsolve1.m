% FINSOLVE1 Solves continuous time asset pricing problems
% USAGE
%   c=finsolve1(model,fspace,alg,snodes,N);
% INPUTS
%   model   : a model structure (see below)
%   fspace  : name of projection space structure
%   alg     : algorithm used (lines, explicit, implicit or CN)
%   snodes  : collocation nodes for the states
%   N       : number of time steps
% OUTPUTS
%   c       : value function projection coefficients (fspace.n by N+1)
%
% The model structure has the following fields:
%    model.func     - the name of the model function file
%    model.T        - the time to maturity of the asset
%    mondel.params  - additional parameters to be passed to model.func
%
% The model function file must have the following syntax:
%    out=func(flag,S,additional parameters);
%      switch flag
%       case 'rho'
%         out= instantaneous risk-free interest rate
%       case 'mu'
%         out= drift on the state process
%       case 'sigma'
%         out = volatility on the state process
%       case 'delta'
%         out = the payout flow (dividend) on the derivative asset
%       case 'V0'
%         out = exercise value of the asset
%       end
%
% Note: The text documentation has an argument, t, which is not
%  incorporated in the current version and is not used in the examples.

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,V,A]=finsolve1(model,fspace,alg,S,N)

keepall = optget('finsolve1','keepall',0);

if ~exist('alg','var') | isempty(alg), alg='lines'; end
if ~exist('N','var') | isempty(N), N=1; end

probfile=model.func;
T=model.T;
if isfield(model,'params')
  params=model.params;
else
  params={};
end
if isfield(model,'american')
  american=model.american;
else
  american=0;
end

% compute basis matrices
n=fspace.n;
Phi0=funbas(fspace,S,0);   
Phi1=funbas(fspace,S,1);    
Phi2=funbas(fspace,S,2); 
  
% Compute collocation matrix 
mu=feval(probfile,'mu',S,params{:});
sigma=feval(probfile,'sigma',S,params{:});
rho=feval(probfile,'rho',S,params{:});
v=0.5*sigma.*sigma;
B = spdiags(mu,0,n,n)*Phi1+spdiags(v,0,n,n)*Phi2 ...
      -spdiags(rho,0,n,n)*Phi0;
Phii=inv(Phi0);
B=Phii*B;

delta=feval(probfile,'delta',S,params{:});
hasdiv=~isempty(delta);
if hasdiv, a=Phii*delta; end

Delta=T/N;
switch alg
case 'lines'
  A=expm(full(Delta*B));
  if hasdiv, a=(A-speye(n))*(B\a); end
case 'explicit'
  A=speye(n)+Delta*B;
  if hasdiv, a=Delta*a;end
case 'implicit'
  A=full(inv(speye(n)-Delta*B));
  if hasdiv, a=A*(Delta*a);end
case 'CN'
  B=(Delta/2)*B;
  A=speye(n)-B;
  if hasdiv, a=A\(Delta*a); end
  A=full(A\(speye(n)+B));
otherwise
  error('Method option is invalid')
end

V0=feval(probfile,'V0',S,params{:});
if keepall
  c=zeros(n,N+1);
else
  c=zeros(n,1);
end
c(:,1)=Phii*V0;
for i=2:N+1
  if keepall
    if hasdiv, c(:,i)=a+A*c(:,i-1);
    else, c(:,i)=A*c(:,i-1);
    end
  else
    if hasdiv, c=a+A*c;
    else, c=A*c;
    end
  end
  if american
    V=max(V0,Phi*c);
    c=Phii*V;
  end
end
