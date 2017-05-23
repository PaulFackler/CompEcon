% REMSOLVE  Solves rational expectations models
% USAGE
%   [c,scoord,x,h,f,resid] = remsolve(model,fspace,scoord,x,h);
% INPUTS
%   model     : rational expectations model (structure array)
%   fspace    : projection space structure (structure array)
%   scoord    : n-vector or grid (cell array) of state nodes
%   x         : initial guess for equilibrium responses at nodes
%   h         : initial guess for expectation variable at nodes
% OUTPUTS
%   c         : expectation function approximation basis coefficients
%   scoord    : residual evaluation coordinates (cell array for n>1)
%   x         : equilibrium responses at evaluation points
%   h         : expectation variables at evaluation points
%   f         : marginal arbitrage benefits at evaluation points
%   resid     : rational expectation residuals at evaluation points
% MODEL STRUCTURE FIELDS:
%   func      : model function file (see below)
%   e         : shocks
%   w         : probabilities    
%   params    : additional parameters to function file
% MODEL FUNCTION FILE FORMAT:
%   [out1,out2,out3] = func(flag,s,x,e,additional parameters)
%   if flag = 'b' returns bound function
%                xl:n.dx, xu:n.dx
%   if flag = 'f' returns reward function and derivatives
%                f:n.dx, fx:n.dx.dx, feh:n.dx.dh
%   if flag = 'g' returns transition function and derivatives
%                g:n.ds, gx:n.ds.dx
%   if flag = 'h' returns expectation function and derivatives
%                h:n.dh, hx:n.dh.dx, hs:n.dh.ds
%   where n  = number of collocation nodes 
%         ds = state space dimension
%         dx = response space dimension
%         dh = expectation space dimension
% USER OPTIONS FOR REMSOLVE (SET WITH optset('remsolve',option,value)):
%   tol       : convergence tolerance
%   maxit     : maximum number of iterations
%   nres      : nres*fspace.n uniform nodes to evaluate residual
%   showiters : 0/1, 1 to display iteration results
% USER OPTIONS FOR ARBIT (SET WITH optset('arbit',option,value)):
%   maxit     : maximum number of iterations for CP
%   tol       : convergence tolerance for CP
%   lcpmethod : 'minmax' or 'smooth'

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

% 5/4/2010 corrected calling variables in line 99

function [c,scoord,x,h,f,resid] = remsolve(model,fspace,scoord,x,h)

% SET PARAMETER DEFAULTS
  tol       = optget('remsolve','tol',sqrt(eps));
  maxit     = optget('remsolve','maxit',500);
  nres      = optget('remsolve','nres',10);
  showiters = optget('remsolve','showiters',1);
  arbmaxit  = optget('arbit','maxit',300);
  arbtol    = optget('arbit','tol',sqrt(eps));
  lcpmethod = optget('arbit','lcpmethod','minmax');

% UNPACK MODEL STRUCTURE
  func=model.func;
  params=model.params;
  if ~isfield(model,'e'); e=0; else e=model.e; end;
  if ~isfield(model,'w'); w=1; else w=model.w; end;
  
% DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  ni  = fspace.n;                  % number of collocation coordinates by state dimension 
  N = prod(ni);                    % number of collocation states
  n = length(ni);                   % dimension of state space
  m = size(x,2);                   % dimension of response space
  
% COMPUTE COLLOCATION NODES & INTERPOLATION MATRIX
%  scoord = funnode(fspace);          % state collocation coordinates
  s   = gridmake(scoord);            % state collocation nodes
  if nargin<5 || isempty(h)
    h = feval(func,'h',s,x,[],[],params{:});
  end
  p = size(h,2);                    % dimension of expectation space
  Phi = funbasx(fspace);            % collocation matrix
  c = funfitxy(fspace,Phi,h);       % initial basis coefficients
 
% SOLVE EULER EQUATION
  tic
  for it=1:maxit                                               % perform function iterations
    cold = c;                                                  % store old basis coefficients
    [f,x,h] = arbit(s,x,c,fspace,func,params,...
          e,w,n,m,p,arbmaxit,arbtol,lcpmethod);                % update expectation and response variables
    c = funfitxy(fspace,Phi,h);                                % update basis coefficient
    change = norm(c-cold,inf);                                 % compute change
    if showiters, fprintf ('%4i %10.1e\n',it,change), end      % print progress
    if change<tol, break, end;                                 % convergence check  
  end
  if showiters, fprintf('Elapsed Time = %7.2f Seconds\n',toc); end

% CHECK STATE TRANSITION SATISFY BOUNDS
  snmin=inf; snmax=-inf;
  for k=1:length(model.w);
    kk = k+zeros(N,1);
    g = feval(func,'g',s,x,[],e(kk,:),params{:});
    snmin = min(snmin,min(g)); snmax = max(snmax,max(g));  
  end
  if any(snmin<fspace.a-eps), disp('Warning: extrapolating beyond smin'), end;
  if any(snmax>fspace.b+eps), disp('Warning: extrapolating beyond smax'), end;
  
% COMPUTE RESIDUAL
  if nargout>5 
    ind = ones(1,n);
    if isfield(model,'discretestates')
      ind(model.discretestates) = 0;
    end
    if n==1,
      if ind
        ni = nres*ni;
        scoord = linspace(fspace.a,fspace.b,ni)';
      end
    else 
      for i=1:n, 
        if ind(i)
          ni(i) = nres*ni(i);
          scoord{i} = linspace(fspace.a(i),fspace.b(i),ni(i))';
        end
      end
    end
    s = gridmake(scoord);
    x = funeval(funfitxy(fspace,Phi,x),fspace,scoord);           % rough guess for responses at evaluation points
    [f,x,h] = arbit(s,x,c,fspace,func,params,...
              e,w,n,m,p,arbmaxit,arbtol,lcpmethod);              % shadow prices and actions at evaluation points
    resid = h-funeval(c,fspace,s);                               % residual at evaluation points
    resid = reshape(resid,[ni p]);                               % reshape residual for plotting
  end
  
% RESHAPE OUTPUT
  x = reshape(x,[ni m]); 
  h = reshape(h,[ni p]);
  f = reshape(f,[ni m]);
return



% ARBIT - Solves arbitrage equation at specified nodes
  function [f,x,h] = arbit(s,x,c,fspace,func,params,e,w,n,m,p,maxit,tol,lcpmethod)

% COMPUTE BOUNDS
  [xl,xu] = feval(func,'b',s,x,[],[],params{:});

  N=size(s,1);  
% SOLVE FIRST ORDER CONDITIONS
  for it=1:maxit
    [f,fx]=equilibrium(s,x,c,fspace,func,params,e,w,n,m,p);
    [lcpf,deltax] = lcpstep(lcpmethod,x,xl,xu,reshape(f,N,m),fx);
    x = x + deltax;
    if norm(deltax(:))< tol, break, end;
    if any(isnan(x)), error('NaNs encountered'); end
  end
  h = feval(func,'h',s,x,[],[],params{:});
return



% EQUILIBRIUM Computes the equilibrium condition at arbitrary s and x
function [f,fx]=equilibrium(s,x,c,fspace,func,params,e,w,n,m,p)
  K  = length(w);
  N=size(s,1);  
  eh = zeros(N,p);
  ehder = zeros(N,p,m);
  for k=1:K
    kk = k+zeros(N,1);
    [g,gx] = feval(func,'g',s,x,[],e(kk,:),params{:});
    [hnext,hnextder] = fund(c,fspace,g,1);
    eh     = eh    + w(k)*hnext;
    ehder  = ehder + w(k)*arraymult(hnextder,gx,N,p,n,m);
    clear hnext hnextder gx g
  end
  [f,fx,feh] = feval(func,'f',s,x,eh,[],params{:});
  fx = fx + arraymult(feh,ehder,N,m,p,m);
return