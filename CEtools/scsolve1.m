% SCSOLVE1 Solves 1 dimensional stochastic control problems
% Simplified to handle the single state case
% USAGE
%   cv=scsolve1(model,fspace,snodes,v);
% See SCSOLVE for details
%
% Solves the Bellman equation using policy function iteration:
%   rho*V=max_x f(s,x)+g(s,x)V_s(s)+0.5sigma(s)*sigma(s)*V_ss(s)

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function cv=scsolve1(model,fspace,snodes,v,x);
  maxiters=100;
  tol=sqrt(eps);
  scfile=model.func;
  params=model.params;

% Compute part of collocation matrix that does not depend on x
  n=size(snodes,1);
  rho=feval(scfile,'rho',snodes,[],[],params{:});
  Phi0=funbas(fspace,snodes,0);
  sigma=feval(scfile,'sigma',snodes,[],[],params{:});
  if isempty(sigma)
    B=spdiags(rho,0,n,n)*Phi0;
  else
    v=0.5*sigma.*sigma;
    B=spdiags(rho,0,n,n)*Phi0-spdiags(v,0,n,n)*funbas(fspace,snodes,2);
  end
% The part of the collocation matrix that does depend on x
  Phi1=funbas(fspace,snodes,1);
% Initialize coefficient vector
  if nargin>4 & isempty(v)
    f=feval(scfile,'f',snodes,x,[],params{:});
    g=feval(scfile,'g',snodes,x,[],params{:});
    cv=(B-spdiags(g,0,n,n)*Phi1)\f; 
    v=Phi0*cv;
  elseif nargin<4 | isempty(v)
    cv=zeros(n,1);
    v=zeros(n,1);
  else
    cv=Phi0\v;
  end
  v0=v;
% Policy function iteration loop
  for iters=1:maxiters
    Vs=Phi1*cv;
    x=feval(scfile,'x',snodes,[],Vs,params{:});
    f=feval(scfile,'f',snodes,x,[],params{:});
    m=feval(scfile,'g',snodes,x,[],params{:});
    cv=(B-spdiags(m,0,n,n)*Phi1)\f;
    v0=v;
    v=Phi0*cv;
    e=max(abs(v-v0));
    if e<tol, break; end
  end
  if iters>=maxiters,                           % print warning message
    disp(['Algorithm did not converge. Maximum error: ' num2str(e)]);
  end