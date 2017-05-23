% GAMESOLVE  Solves discrete time continuous-state/action Bellman equations for dynamic games
% USAGE
%   [c,scoord,v,x,resid] = gamesolve(model,fspace,v,x);
% INPUTS
%   model     : dynamic program model structure
%   fspace    : name of projection space structure
%   v         : initial guess for values or shadow prices at nodes
%                 (if empty, initial values are computed using x)
%   x         : initial guess for optimal actions at nodes
% OUTPUTS
%   c         : value function approximation basis coefficients
%   scoord    : residual evaluation coordinates (cell array for ds>1)
%   v         : value or shadow price function at evaluation points
%   x         : optimal action at evaluation points
%   resid     : Bellman equation residuals at evaluation points
% MODEL STRUCTURE FIELDS
%   func      : function file (see below)
%   discount  : discount factor
%   e         : shocks
%   w         : probabilities    
%   params    : additional parameters to function file
% FUNCTION FILE FORMAT
%   [out1,out2,out3] = func(flag,s,x,e,additional parameters)
%   if flag = 'b' returns bound function and derivatives
%      xl:ns.dx, xu:ns.dx, xls:ns.dx.ds, xus:ns.dx.ds
%   if flag = 'f' returns reward function and derivatives
%      f:ns.1, fx:ns.dx, fxx:ns.dx.dx, fs:ns.ds, fsx:ns.ds.dx, fss:ns.ds.ds
%   if flag = 'g' returns transition function and derivatives
%      g:ns.1, gx:ns.dx, gxx:ns.dx.dx, gs:ns.ds
%   where ns=number of collocation states, ds=state space dimension, dx=action space dimension
% USER OPTIONS (SET WITH OPSET)
%   tol       : convergence tolerance
%   maxit     : maximum number of iterations
%   nres      : nres*fspace.n uniform nodes to evaluate residual
%   showiters : 0/1, 1 to display iteration results
%   alg       : 'newton' (the default) or 'funcit' 
% USER OPTIONS FOR GAMESOLVE_VMAX
%   tol          : convergence tolerance used by CP solver (default: 5e-8);
%   maxit        : maximum number of iterations used by CP solver (default: 50)
%   maxbacksteps : maximum number of backsteps used by CP solver (default: 0)
%   lcpmethod    : 'minmax' (the default) or 'smooth'
%
% Copyright (c) 1997-2001, Mario J. Miranda & Paul L. Fackler
% miranda.4@osu.edu, paul_fackler@ncsu.edu

function [c,scoord,v,x,resid] = gamesolve(model,fspace,v,x)

% SET PARAMETER & SHOCK DISTRIBUTION DEFAULTS
  tol       = optget('gamesolve','tol',sqrt(eps));
  maxit     = optget('gamesolve','maxit',500);
  nres      = optget('gamesolve','nres',10);
  showiters = optget('gamesolve','showiters',1);
  alg       = optget('gamesolve','alg','newton');

  if ~isfield(model,'e'); model.e=0; end;
  if ~isfield(model,'w'); model.w=1; end;
  if isfield(model,'T') & model.T<inf, alg='finite'; end
  
% DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  n  = fspace.n;                                              % number of collocation coordinates by state dimension 
  ns = prod(n);                                               % number of collocation states
  ds = length(n);                                             % dimension of state space
  dx = size(x,2);                                             % dimension of action space
  
% COMPUTE COLLOCATION NODES AND INTERPOLATION MATRIX
  scoord = funnode(fspace);                                   % state collocaton coordinates
  s   = gridmake(scoord);                                     % state collocaton nodes
  Phi = funbas(fspace);                                       % collocation matrix

% DETERMINE INITIAL COEFFICIENT VALUES
  if isempty(v)
    delw = model.discount*model.w;
    vc = zeros(size(s,1),size(Phi,2));
    for k=1:length(model.w)
      kk = k*ones(ns,1);
      g = feval(model.func,'g',s,x,model.e(kk,:),model.params{:});
      phinext = funbas(fspace,g);
      vc = vc + delw(k)*phinext;
    end
    f=[feval(model.func,'f1',s,x,[],model.params{:}) ...
       feval(model.func,'f2',s,x,[],model.params{:})];
    c=(Phi-vc)\f;

  else
    c   = Phi\v;                                                % initial basis coefficients
  end
      
% SOLVE BELLMAN EQUATION
  tic
  switch alg
    case 'funcit'                                              % BELLMAN FUNCTION ITERATION
    for it=1:maxit                                             % perform iterations
      cold = c;                                                % store old basis coefficients
      [v,x] = vmax(s,x,c,fspace,model);                        % update value function and policy
      c = Phi\v;                                               % update basis coefficient
      change = norm(c(:)-cold(:),inf);                         % compute change
      if showiters, fprintf ('%4i %10.1e\n',it,change), end    % print progress
      if change<tol, break, end;                               % convergence check  
    end
  case 'newton'                                                % BELLMAN NEWTON METHOD
    for it=1:maxit                                             % perform iterations
      cold = c;                                                % store old basis coefficients
      [v,x,vderc] = vmax(s,x,c,fspace,model);                  % update value function and policy
      c = c - [Phi-vderc]\[Phi*c-v];                           % update basis coefficient
      change = norm(c(:)-cold(:),inf);                         % compute change
      if showiters, fprintf ('%4i %10.1e\n',it,change), end    % print progress
      if change<tol, break, end;                               % convergence check  
    end
  end
  fprintf('Elapsed Time = %7.2f Seconds\n',toc)
  
% CHECK STATE TRANSITION SATISFY BOUNDS
  snmin=inf; snmax=-inf;
  for k=1:length(model.w);
    kk = k*ones(ns,1); e = model.e(kk,:);
    g = feval(model.func,'g',s,x,e,model.params{:});
    snmin = min(snmin,min(g)); snmax = max(snmax,max(g));  
  end
  if any(snmin<fspace.a-eps), disp('Warning: extrapolating beyond smin'), end;
  if any(snmax>fspace.b+eps), disp('Warning: extrapolating beyond smax'), end;
  
% COMPUTE RESIDUAL
  if ~all(alg=='finite') & nargout>4
    ind = ones(1,ds);
    if isfield(model,'discretestates')
      ind(model.discretestates) = 0;
    end
    if ds==1,
      if ind
        n = nres*n;
        scoord = linspace(fspace.a,fspace.b,n)';
      end
    else 
      for i=1:ds, if ind(i)
        n(i) = nres*n(i);
        scoord{i} = linspace(fspace.a(i),fspace.b(i),n(i))';
      end, end
    end
    s = gridmake(scoord);
    x = funeval(Phi\x,fspace,scoord);   % rough guess for actions at evaluation points
    [v,x] = vmax(s,x,c,fspace,model);   % values and actions at evaluation points
    resid = v-funeval(c,fspace,s);      % residual at evaluation points
    resid = reshape(resid,[n 2]);       % reshape residual for plotting
  end
    
% RESHAPE OUTPUT
  x = reshape(x,[n dx]); 
  v = reshape(v,[n  2]);

 
    
% VMAX - Solves Bellman Equation
  function [v,x,vc] = vmax(s,x,c,fspace,model)

% SET CONVERGENCE PARAMETER DEFAULTS
  tol           = optget('gamesolve_vmax','tol',5e-10);
  maxit         = optget('gamesolve_vmax','maxit',500);
  maxbacksteps  = optget('gamesolve_vmax','maxbacksteps',5);
  lcpmethod     = optget('gamesolve_vmax','lcpmethod','minmax');
    
% COMPUTE LOCAL CONSTANTS
  [ns,ds] = size(s);
  dx  = size(x,2);
  delw = model.discount*model.w;
  K = length(delw);
  dxi = dx/2;
  v=zeros(ns,2);
  params=model.params;
  func=model.func;
  e=model.e;

% COMPUTE BOUNDS      
  [xl,xu] = feval(func,'b',s,x,[],params{:}); 

% SOLVE FIRST ORDER CONDITIONS
  for i=1:2
    if i==1
      f = 'f1'; j = [1:dxi];
    else
      f = 'f2'; j = [dxi+1:dx];
    end
    for it=1:maxit
      [v(:,i),vx,vxx] = valfunc(s,x,c(:,i),fspace,model,f,j);
      [vx,deltax] = lcpstep(lcpmethod,x(:,j),xl(:,j),xu(:,j),vx,vxx);
      err = max(abs(vx),[],2);
      if all(err<tol), break, end;
      eold=inf;
      if maxbacksteps<1
        x(:,j) = x(:,j)+deltax;
      else
        for k=1:maxbacksteps
          xnew = x;
          xnew(:,j) = xnew(:,j)+deltax;
          [v(:,i),vx] = valfunc(s,xnew,c(:,i),fspace,model,f,j);
          vx = lcpstep(lcpmethod,x(:,j),xl(:,j),xu(:,j),vx);
          enew = max(abs(vx),[],2);
          ind = find(eold>enew & enew>err);
          if isempty(ind), break; end
          eold = enew;
          deltax(ind,:) = deltax(ind,:)/2;
        end 
        x = xnew;
      end
    end
  end
  
% COMPUTE dv/dc  
  if nargout>2
    delw = model.discount*model.w;
    vc = zeros(size(s,1),size(c,1));
    for k=1:length(model.w)
      kk = k*ones(ns,1);
      g = feval(func,'g',s,x,e(kk,:),params{:});
      phinext = funbas(fspace,g);
      vc = vc + delw(k)*phinext;
    end
  end  
  
 

% VALFUNC  Evaluates Bellman Optimand
  function [v,vx,vxx]=valfunc(s,x,c,fspace,model,f,j)
  
% COMPUTE LOCAL CONSTANTS
  [ns,ds] = size(s);
  dx  = size(x,2);
  dxx = dx*dx;
  delw = model.discount*model.w;
  K = length(delw);
  dxi = length(j);
  e=model.e;
  func=model.func;
  params=model.params;

  if nargout<3
    [v,vx] = feval(func,f,s,x,[],params{:});
    vx = reshape(vx,ns,1,dx);
    for k=1:K
      kk = k*ones(ns,1);
      [g,gx] = feval(func,'g',s,x,e(kk,:),params{:});
      [vnext,vs] = fund(c,fspace,g,1);
      v  = v  + delw(k)*vnext;
      vx = vx + delw(k)*arraymult(vs,gx,ns,1,ds,dx);
    end
    clear g gx 
    vx = vx(:,:,j);
    vx = reshape(vx,ns,dxi);
  else
    [v,vx,vxx] = feval(func,f,s,x,[],params{:});
    vx  = reshape(vx,ns,1,dx);
    vxx = reshape(vxx,ns,dx,dx);
    for k=1:K
      kk = k*ones(ns,1); 
      [g,gx,gxx] = feval(func,'g',s,x,e(kk,:),params{:});
      [vnext,vs,vss] = fund(c,fspace,g,1);
      v   = v   + delw(k)*vnext;
      vx  = vx  + delw(k)*arraymult(vs,gx,ns,1,ds,dx);
      vxx = vxx + delw(k)*(reshape(arraymult(vs,gxx,ns,1,ds,dxx),ns,dx,dx) ...
                + arraymult(permute(gx,[1 3 2]),arraymult(vss,gx,ns,ds,ds,dx),ns,dx,ds,dx));
    end  
    clear g gx gxx vss
    vx = vx(:,:,j);
    vx = reshape(vx,ns,dxi);
    vxx = vxx(:,j,j);
    vxx = reshape(vxx,ns,dxi*dxi);
  end