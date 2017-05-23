% SCSOLVE Solves stochastic control problems
% Solves continuous time Bellman equation using policy function iteration:
%   rho(s)*V = max_x f(s,x) + V_s(s)g(s,x) + 0.5trace(sigma(s)*sigma(s)'*V_ss(s))
% USAGE
%   [cv,scoord,v,x,resid]=scsolve(model,fspace,snodes,v,x);
% INPUTS
%   model   : a model structure (see below)
%   fspace  : approximation function space structure variable (see fundef)
%   snodes  : collocation nodes for the states
%   v       : initial guess for value function at snodes (default=0)
%   x       : initial guess of the optimal control
% OUTPUTS
%   cv      : value function projection coefficients
%   scoord  : residual evaluation points
%               (a cell array for multi-dimensional states;
%                use gridmake to expand)
%   v       : value function at evaluation points
%   x       : optimal control at evaluation points
%   resid   : Bellman equation residuals
%
% The model structure should contain the following fields:
%   func    : model function file (see below)  
%   params  : additional parameters to pass to model function file
%
% The model function file should have the format
%    out=scfile(flag,s,x,Vs,additional parameters)
%      switch flag
%        case 'x'
%          Return a matrix of optimal actions given
%             s and Vs
%        case 'f'
%          Return a matrix f representing the reward function
%        case 'g'
%          Return a matrix g representing the drift function
%              of the state transition equation
%        case 'sigma'
%          Return a matrix sigma representing the diffusion function
%              of the state transition equation
%        case 'rho'
%          Return a vector representing the (state contingent)
%              discount rates
%      end
% For problems with d states and p actions the sizes of these matrices
% should be as follows (with n=# of rows in s, x and Vs):
%            x     : n x p
%            f     : n x 1
%            g     : n x d
%            sigma : n x d x d
%            rho   : n x 1
%
% USER OPTIONS (SET WITH OPTSET)
%   maxit         : maximum number of iterations
%   tol           : convergence tolerance
%   showiters     : 0/1, 1 to display iteration results
%   storeexpanded : 1 if basis matrices are stored in expanded form
%                   (set to 0 if encountering memory problems)
%   nres          : nres*fspace.n nodes are used to evaluate the
%                      residual function (for continuous states)

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [cv,scoord,v,x,resid]=scsolve(model,fspace,snodes,v,x)

% Set default options
  maxit         = optget('scsolve','maxiters',100);
  tol           = optget('scsolve','tol',sqrt(eps));
  nres          = optget('scsolve','nres',10);
  showiters     = optget('scsolve','showiters',1);
  StoreExpanded = optget('scsolve','StoreExpanded',1);

% Initialize nodes and initial value function
  d=fspace.d;
  n=fspace.n;

  if nargin<3, snodes=funnode(fspace); end

  S=gridmake(snodes);            % expand state grid

% Unpack problem structure variable
  scfile=model.func;
  params=model.params;

% Compute part of collocation matrix that does not depend on x
  bases=funbasx(fspace,snodes,[0;1;2]);
  B=ctbasemake(feval(scfile,'rho',S,[],[],params{:}),bases,0);
  sigma=feval(scfile,'sigma',S,[],[],params{:});
  if ~isempty(sigma)
    B=B-ctbasemake(sigma,bases,2);
  end
  clear sigma;

  % Store Phi1 in expanded form as it is used repeatedly
  % Change if memory is limited
  if StoreExpanded
    Phi1=ctbasemake([],bases,1); 
  else
    Phi1=bases;
  end

% Initialize coefficient vector
  if nargin>4 & isempty(v)
    f=feval(scfile,'f',S,x,[],params{:});
    g=feval(scfile,'g',S,x,[],params{:});
    cv=(B-ctbasemake(g,Phi1,1))\f; 
    v=funeval(cv,fspace,bases);
  elseif nargin<4 | isempty(v)
    cv=zeros(size(S,1),1);
    v=funeval(cv,fspace,bases);
  else
    cv=ckronxi(bases.vals(1,:),v,d:-1:1);
  end


% Policy function iteration loop
  for iter=1:maxit
    v0=v;
    if StoreExpanded
      Vs=zeros(size(Phi1{1},1),d);
      for i=1:d, Vs(:,i)=Phi1{i}*cv; end
    else
      Vs=squeeze(funeval(cv,fspace,bases,eye(d)));
    end
    x=feval(scfile,'x',S,[],Vs,params{:});
  %  clear Vs    
    f=feval(scfile,'f',S,x,[],params{:});
    g=feval(scfile,'g',S,x,[],params{:});
    cv=(B-ctbasemake(g,Phi1,1))\f; 
    v=funeval(cv,fspace,bases);
    if any(isnan(v) | isinf(v))
      error('NaNs or Infs encountered');
    end
  %  clear f g
    e=max(abs(v-v0));
    if showiters, fprintf('%3i %12.4e\n',iter,e); end
    if e<tol, break; end
  end
  if iter>=maxit,                           % print warning message
    disp(['Algorithm did not converge. Maximum error: ' num2str(e)]); 
  end

  if nargout>1
    if d==1,
      n = nres*n+1;
      scoord = linspace(fspace.a,fspace.b,n)';
    else 
      for i=1:d
        n(i) = nres*n(i)+1;
        scoord{i} = linspace(fspace.a(i),fspace.b(i),n(i))';
      end
    end
    bases=funbasx(fspace,scoord,[0;1;2]);
    v=funeval(cv,fspace,bases);
    Vs=squeeze(funeval(cv,fspace,bases,eye(d))); 
    S = gridmake(scoord);
    x=feval(scfile,'x',S,[],Vs,params{:});
    if nargout>4
      f=feval(scfile,'f',S,x,[],params{:});
      rho=feval(scfile,'rho',S,x,[],params{:});
      resid=rho.*v-f;
      g=feval(scfile,'g',S,x,[],params{:});
      resid=resid-sum(Vs.*g,2);
      sigma=feval(scfile,'sigma',S,[],[],params{:});
      if ~isempty(sigma)
        resid=resid-ctbasemake(sigma,bases,2)*cv;
      end
      resid=reshape(resid,[n 1]);
    end
    v=reshape(v,[n 1]);
    x=reshape(x,[n size(x,2) 1]);
  end  
    

