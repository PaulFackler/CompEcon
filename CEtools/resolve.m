% RESOLVE  Solves rational expectations models
% USAGE
%   [c,scoord,x,h,f,resid] = resolve(model,fspace,options,c);
% INPUTS
%   model     : rational expectations model (structure array)
%   fspace    : projection space structure (structure array)
%   options   : algorithm options (structure array)
%   c         : initial value for approximation coefficients
% OUTPUTS
%   c         : expectation function approximation basis coefficients
%   scoord    : residual evaluation coordinates (cell array for d>1)
%   x         : equilibrium responses at evaluation points (N x m)
%   z         : expectation variables at evaluation points (N x p)
%   f         : marginal arbitrage benefits at evaluation points (N x m)
%   resid     : residual at evaluation points (N x m or N x p)
% MODEL STRUCTURE FIELDS:
%   func      : function file (see below) 
%   params    : additional parameters passed to function file
%   e         : shocks
%   w         : probabilities   
%   explicit  : 1 if x(s,z) such that f(s,x,z) has explict solution
%   noxnext   : 1 if h can be written as h(s,x,e,snext)
% MODEL FUNCTION FILE FORMAT:
%   [out1,out2,out3,out4] = func(flag,s,x,z,e,snext,xnext,additional parameters)
%   if flag = 'f' returns reward function and (optionally) its Jacobian
%      f: N x m, fx: N x m x m
%   if flag = 'x' returns response give s and z=Eh
%      x: N x m, [optional]  xs: N x m x d,  xz: N x m x p
%   if flag = 'g' returns transition function 
%      g: N x d, [optionally] gx: N x d x m
%   if flag = 'h' returns expectation function
%      h: N x p, [optionally] hx: N x p x m, hsnext: N x p x s, hxnext: N x p x m
%   if flag = 'b' returns bound functions
%      xl: N x m, xu: N x m
%   where N = number of rows in s 
%         d = state space dimension
%         m = response space dimension
%         p = expectation space dimension
% Note: 'f' and 'b' flags only called is model.explicit=0
% OPTIONS FOR RESOLVE
%   expapprox  : 1 if expectation approximation used (default: 1)
%   usebroyden : 1 if Broyden's method used (default: 0)
%   lowmemory  : 1 to use less memory (typically executes more slowly)
%                  (default: 0)
%   stepsize   : alpha value on [0,2] used to induce stabilty
%                or accelerate convergence (default: 1)
%   tol        : convergence tolerance (default: sqrt(eps))
%   maxit      : maximum number of iterations (default: 500)
%   showiters  : 0/1, 1 to display iteration results  (default: 0)
%   checks     : 0/1, 1 to perform checks on stte variable extrapolation
%   nres       : nres*fspace.n uniform nodes to evaluate residual
%                  (default: 10)
%   xtol       : convergence tolerance for CP (default: sqrt(eps))
%   xmaxit     : maximum number of iterations for CP (default: 500)
%   lcpmethod  : 'minmax' or 'smooth' (default: 'minmax')

% Copyright (c) 2001-2003,  Paul L. Fackler
% paul_fackler@ncsu.edu

function [c,scoord,x,z,f,resid] = resolve(model,fspace,options,c)

% SET PARAMETER DEFAULTS
  if nargin<3, options=[]; end
  getopts(options,...
          'expapprox',    1,...
          'usebroyden',   0,...
          'lowmemory',    0,...
          'stepsize',     1,...
          'tol',  sqrt(eps),...
          'maxit',      500,...
          'showiters',    1,...
          'checks',       1,...
          'nres',        10,...
          'xtol', sqrt(eps),...
          'xmaxit',     500,...
          'lcpmethod', 'minmax');

% UNPACK MODEL STRUCTURE
  func=model.func;
  params=model.params;
  if ~isfield(model,'e'); e=0; else, e=model.e; end;
  if ~isfield(model,'w'); w=1; else, w=model.w; end;
  if ~isfield(model,'explicit'), explicit=0; else, explicit=model.explicit; end
  if ~isfield(model,'noxnext'),   noxnext=0; else,  noxnext=model.noxnext; end 
  
% DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  ni  = fspace.n;              % number of collocation coordinates for each state 
  n = prod(ni);                % number of collocation states
  d = length(ni);              % dimension of state space  
  K = length(w);               % number of shock values

% COMPUTE COLLOCATION NODES & INTERPOLATION MATRIX
  scoord = funnode(fspace);    % state collocation coordinates
  s   = gridmake(scoord);      % state collocation nodes
  Phi = funbasx(fspace);       % collocation matrix
  % initial basis coefficients
  if nargin<4 | isempty(c)
    z=feval(func,'zinit',s,[],[],[],[],[],params{:});
    x=feval(func,'xinit',s,[],[],[],[],[],params{:});
    if expapprox
      c = funfitxy(fspace,Phi,z);
    else
      c = funfitxy(fspace,Phi,x);
    end
  else
    if expapprox
      z = funeval(c,fspace,Phi);
      if explicit
        x=feval(func,'x',s,[],z,[],[],[],params{:});
      else
        x=feval(func,'xinit',s,[],[],[],[],[],params{:});
      end 
    else
      x = funeval(c,fspace,Phi);
      ee=w'*e;
      z=feval(func,'h',s,x,[],ee(ones(n,1),:),s,x,params{:});
    end
  end
  p = size(z,2);                   % dimension of expectation space 
  m = size(x,2);                   % dimension of response space

  if lowmemory, ee=e;
  else,         ee=e(repmat(1:K,1,n),:);
  end

% SOLVE EQUILIBRIUM
  tic
  if usebroyden
    % use deterministic case to initialize Jacobian
    r=funeval(c,fspace,s);
    if nargout(func)<4
      dc=fdjac1('reres',c(:),[],s,r,fspace,func,params,w'*e,1,explicit,noxnext,...
           d,m,p,n,1,Phi,xmaxit,xtol,lcpmethod,expapprox,1);
    else
      if expapprox
        if explicit
          dc=reresJe(c,s,r,fspace,func,params,w'*e,noxnext,d,m,p,n,Phi,explicit);
        else
          warning('Analytic derivatives for non-explicit expectation method not implemented')
          dc=fdjac1('reres',c(:),s,r,fspace,func,params,w'*e,1,explicit,noxnext,...
            d,m,p,n,1,Phi,xmaxit,xtol,lcpmethod,expapprox,1);
        end 
      else
        if explicit
          dc=reresJr(c,s,r,fspace,func,params,w'*e,noxnext,d,m,p,n,Phi);
        else
          dc=reresJf(c,s,r,fspace,func,params,w'*e,noxnext,d,m,p,n,Phi);
        end
      end
    end
    optset('broyden','initb',inv(dc));
    optset('broyden','maxit',maxit);
    optset('broyden','showiters',showiters);
    optset('broyden','tol',tol);
    c=broyden('reres',c(:),s,x,fspace,func,params,ee,w,explicit,noxnext,...
           d,m,p,n,K,Phi,xmaxit,xtol,lcpmethod,expapprox,lowmemory); 
    [r,x,z,cx]=reres(c(:),s,x,fspace,func,params,ee,w,explicit,noxnext,...
           d,m,p,n,K,Phi,xmaxit,xtol,lcpmethod,expapprox,lowmemory);
    if expapprox, c=reshape(c,n,p);
    else,         c=reshape(c,n,m);
    end
    optset('broyden','defaults');
  elseif 1
    iPhi=fliplr(Phi.vals);
    for i=1:length(iPhi) 
      iPhi{i}=inv(iPhi{i});
    end
    for it=1:maxit                                 % perform function iterations
      if expapprox
        zold=z;
        [z,x]=Ehe(c,s,z,fspace,func,params,ee,w,explicit,noxnext,...
                       d,m,p,n,K,Phi,[],xmaxit,xtol,lcpmethod,lowmemory);
        ctilde = ckronx(iPhi,z);
        r=z-zold;
      else
        xold=x;
        [z,x]=Ehr(c,s,x,fspace,func,params,ee,w,explicit,noxnext,...
                   d,m,p,n,K,Phi,[],xmaxit,xtol,lcpmethod,lowmemory);
        ctilde = ckronx(iPhi,x);
        r=x-xold;
      end
      if any(isnan(ctilde(:)) | isinf(ctilde(:))),
        error('infs or nans encountered');
      end
      c = c+stepsize*(ctilde-c);                              % update basis coefficient
      change = norm(r,inf);                                   % compute change
      if showiters, fprintf ('%4i %10.1e\n',it,change), end   % print progress
      if change<tol, break, end;                              % convergence check  
    end
  end

  if showiters, fprintf('Elapsed Time = %7.2f Seconds\n',toc); end

% CHECK IF STATE TRANSITION SATISFIES BOUNDS
  if checks
    if lowmemory
      snmin=inf; snmax=-inf;
      for k=1:length(model.w);
        kk = k+zeros(n,1);
        g = feval(func,'g',s,x,[],e(kk,:),[],[],params{:});
        snmin = min(snmin,min(g)); snmax = max(snmax,max(g));  
      end
    else 
      ind=(1:n); ind=ind(ones(1,K),:);
      ss=s(ind,:);
      xx=x(ind,:);
      g = feval(func,'g',ss,xx,[],ee,[],[],params{:});
      snmin = min(g); 
      snmax = max(g);  
    end
    if any(snmin<fspace.a*(1-eps)), disp('Warning: extrapolating beyond smin'), end;
    if any(snmax>fspace.b*(1+eps)), disp('Warning: extrapolating beyond smax'), end;
  end

% COMPUTE RESIDUAL
  if nargout>5 
    scoord=cell(1,d);
    for i=1:d, 
      ni(i) = nres*ni(i)+1;
      scoord{i} = linspace(fspace.a(i),fspace.b(i),ni(i))';
    end
    if d==1, scoord=scoord{1}; s=scoord;
    else,    s = gridmake(scoord);
    end
    x = funeval(funfitxy(fspace,Phi,x),fspace,scoord);           % rough guess for responses at evaluation points
   [resid,x,z,cx]=reres(c,s,x,fspace,func,params,e,w,explicit,noxnext,...
         d,m,p,n,K,[],xmaxit,xtol,lcpmethod,expapprox,1);
  end
  f=feval(func,'f',s,x,z,[],[],[],params{:});
% RESHAPE OUTPUT
  %x = reshape(x,[ni m]); 
  %h = reshape(h,[ni p]);
  %f = reshape(f,[ni m]);
return


%%%%%%%%%%%%%   Auxillary functions  %%%%%%%%%%%%%

% Residual function used by RESOLVE with Broyden's method
function [r,x,z,cx]=reres(c,s,x,fspace,func,params,e,w,explicit,noxnext,...
         d,m,p,n,K,Phi,xmaxit,xtol,lcpmethod,expapprox,lowmemory); 
global reres_counter
if nargin==0
  if ~isempty(reres_counter) & reres_counter>0
    disp(reres_counter)
  end
  reres_counter=0;
  return
end
reres_counter=reres_counter+1;

if expapprox
  c=reshape(c,n,p);
  if isempty(Phi),  r=funeval(c,fspace,s);
  else,             r=funeval(c,fspace,Phi);
  end
  n=size(s,1);
  [z,x,cx] = Ehe(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
         d,m,p,n,K,Phi,[],xmaxit,xtol,lcpmethod,lowmemory);
  r=r-z;
else
  c=reshape(c,n,m);
  if isempty(Phi),  r=funeval(c,fspace,s);
  else,             r=funeval(c,fspace,Phi);
  end
  n=size(s,1);
  [z,x] = Ehr(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
              d,m,p,n,K,Phi,[],xmaxit,xtol,lcpmethod,lowmemory);
  if explicit
    r=r-x;
  else
    r=feval(func,'f',s,x,z,[],[],[],params{:});
  end  
  cx=[];
end
r=r(:);
return

% Ehr Computes z at arbitrary s and c with expectation approximation
function [z,x,cx]=Ehe(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
                       d,m,p,n,K,Phi,cx,xmaxit,xtol,lcpmethod,lowmemory)
  z=r;
  if explicit
    x=feval(func,'x',s,[],z,[],[],[],params{:});
  else
    x=resolvef(s,z,x,func,params,n,m,xmaxit,xtol,lcpmethod);
    if ~noxnext & isempty(cx), cx = funfitxy(fspace,Phi,x); end
  end
  xnext=[];
  if lowmemory
    z = zeros(n,p);
    for k=1:K  
      kk = k+zeros(n,1);
      snext = feval(func,'g',s,x,[],e(kk,:),[],[],params{:});
      if ~noxnext
        if explicit
          znext = funeval(c,fspace,snext);
          xnext = feval(func,'x',snext,[],znext,[],[],[],params{:});
        else
          xnext=funeval(cx,fspace,snext);
        end
      end
      z=z+w(k)*feval(func,'h',s,x,[],e(kk,:),snext,xnext,params{:});
    end
  else
    ind=(1:n); ind=ind(ones(1,K),:);
    ss=s(ind,:);
    xx=x(ind,:);
    snext = feval(func,'g',ss,xx,[],e,[],[],params{:});
    if ~noxnext
      if explicit
        znext = funeval(c,fspace,snext);
        xnext = feval(func,'x',snext,[],znext,[],[],[],params{:});
      else
        xnext=funeval(cx,fspace,snext);
      end
    end
    z=reshape(feval(func,'h',ss,xx,[],e,snext,xnext,params{:}),K,n*p);
    z=reshape(w'*z,n,p);
  end
return

% Ehr Computes z at arbitrary s and c with repsonse approximation
function [z,x]=Ehr(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
                   d,m,p,n,K,Phi,cx,xmaxit,xtol,lcpmethod,lowmemory)
  x = r; 
  xnext=[];
  if lowmemory
    z = zeros(n,p);
    kk = zeros(n,1);
    for k=1:K  
      kk(:)=k;
      snext = feval(func,'g',s,x,[],e(kk,:),[],[],params{:});
      if ~noxnext
        xnext = funeval(c,fspace,snext);
      end
      z=z+w(k)*feval(func,'h',s,x,[],e(kk,:),snext,xnext,params{:});
    end
  else
    ind=(1:n); ind=ind(ones(1,K),:);
    ss=s(ind,:);
    xx=x(ind,:);
    snext = feval(func,'g',ss,xx,[],e,[],[],params{:});
    if ~noxnext
      xnext = funeval(c,fspace,snext);
    end
    z=reshape(feval(func,'h',ss,xx,[],e,snext,xnext,params{:}),K,n*p);
    z=reshape(w'*z,n,p);
  end
  if explicit
    x=feval(func,'x',s,[],z,[],[],[],params{:});
  end
return


% Newton's method to solve f(s,x(s,z),z)=0
function x = solvef(s,z,x,func,params,n,m,maxit,tol,lcpmethod)
% COMPUTE BOUNDS
  [xl,xu] = feval(func,'b',s,x,[],[],[],[],params{:});
% SOLVE FIRST ORDER CONDITIONS
  for it=1:maxit
    xold = x;
    [f,fs,fx]=feval(func,'f',s,x,z,[],[],[],params{:});
    [f,deltax] = lcpstep(lcpmethod,x,xl,xu,reshape(f,n,1,m),fx);
    x = x + deltax;
    if norm(deltax(:))< tol, break, end;
    if any(isnan(x)), error('NaNs encountered - ending procecure'); end
  end
return


% Jacobian for expection approximation
function [J]=reresJe(c,s,z,fspace,func,params,e,noxnext,d,m,p,n,Phi,explicit); 
  e=e(ones(n,1),:);

  if explicit
    [x,xs,xz]=feval(func,'x',s,[],z,[],[],[],params{:});
  else
    error('not implemented')
    x=resolvef(s,z,x,func,params,n,m,xmaxit,xtol,lcpmethod);
  end
  [snext,gx] = feval(func,'g',s,x,[],e,[],[],params{:});
  if ~noxnext
    if explicit
      znext = funeval(c,fspace,snext);
      [xnext,xs1,xz1] = feval(func,'x',snext,[],znext,[],[],[],params{:});
    else
      cx = funfitxy(fspace,Phi,x);
      xnext=funeval(cx,fspace,snext);
    end
  end
  [zupdate,hx,hs1,hx1]=feval(func,'h',s,x,[],e,snext,xnext,params{:});

  P1=funbas(fspace,s);
  J=kron(eye(p),P1);

  if any(gx(:)~=0)
    if any(hx1(:)~=0)
      A=xs1 + arraymult(xz1,funeval(c,fspace,snext,eye(d)),n,m,p,d);
      A=arraymult(hx1,A,n,p,m,d);
    else
      A=0;
    end
    if any(hs1(:)~=0)
      A=A+reshape(hs1,n,p,d);
    end
    A=arraymult(A,gx,n,m,d,m);
  else
    A=0;
  end

  if any(hx(:)~=0)
    A=A+reshape(hx,n,p,m);
  end

  if any(A(:)~=0)
    A=arraymult(A,xz,n,p,m,p);
    ind=1:n;
    for i=1:p
      for j=1:p
        J(ind+(i-1)*n,ind+(j-1)*n)=J(ind+(i-1)*n,ind+(j-1)*n)-diag(A(:,i,j))*P1;
      end
    end
  end

  if any(hx1(:)~=0)
    ind=1:n;
    P2=funbas(fspace,snext);
    A=arraymult(hx1,xz1,n,p,m,p);
    for i=1:p
      for j=1:p
        J(ind+(i-1)*n,ind+(j-1)*n)=J(ind+(i-1)*n,ind+(j-1)*n)-diag(A(:,i,j))*P2;
      end
    end
  end
return

% Jacobian for response approximation (explicit option)
function [J]=reresJr(c,s,x,fspace,func,params,e,noxnext,d,m,p,n,Phi); 
  e=e(ones(n,1),:);
  [snext,gx] = feval(func,'g',s,x,[],e,[],[],params{:});
  if ~noxnext
    xnext = funeval(c,fspace,snext);
  end
  [z,hx,hs1,hx1]=feval(func,'h',s,x,[],e,snext,xnext,params{:});
  [xupdate,xs,xz]=feval(func,'x',s,[],z,[],[],[],params{:});

  P1=funbas(fspace,s);
  J=kron(eye(m),P1);

  if any(gx(:)~=0)
    if any(hx1(:)~=0)
      A=arraymult(hx1,funeval(c,fspace,snext,eye(d)),n,m,m,d);
    else
      A=0;
    end
    if any(hs1(:)~=0)
      A=A+reshape(hs1,n,p,d);
    end
    A=arraymult(A,gx,n,m,d,m);
  else
    A=0;
  end

  if any(hx(:)~=0)
    A=A+reshape(hx,n,p,m);
  end

  if any(A(:)~=0)
    A=arraymult(xz,A,n,m,p,m);
    ind=1:n;
    for i=1:m
      for j=1:m
        J(ind+(i-1)*n,ind+(j-1)*n)=J(ind+(i-1)*n,ind+(j-1)*n)-diag(A(:,i,j))*P1;
      end
    end
  end

  if any(hx1(:)~=0)
    P2=funbas(fspace,snext);
    A=arraymult(xz,hx1,n,m,p,m);
    ind=1:n;
    for i=1:m
      for j=1:m
        J(ind+(i-1)*n,ind+(j-1)*n)=J(ind+(i-1)*n,ind+(j-1)*n)-diag(A(:,i,j))*P2;
      end
    end
  end
return

% Jacobian for response approximation (non-explicit option)
function [J]=reresJf(c,s,x,fspace,func,params,e,noxnext,d,m,p,n,Phi); 
  e=e(ones(n,1),:);
  [snext,gx] = feval(func,'g',s,x,[],e,[],[],params{:});
  if ~noxnext
    xnext = funeval(c,fspace,snext);
  end
  [z,hx,hs1,hx1]=feval(func,'h',s,x,[],e,snext,xnext,params{:});
  [fval,fs,fx,fz]=feval(func,'f',s,x,z,[],[],[],params{:});

  P1=funbas(fspace,s);
  J=zeros(n*m,n*m);

  if any(gx(:)~=0)
    if any(hx1(:)~=0)
      A=arraymult(hx1,funeval(c,fspace,snext,eye(d)),n,m,m,d);
    else
      A=0;
    end
    if any(hs1(:)~=0)
      A=A+reshape(hs1,n,p,d);
    end
    A=arraymult(A,gx,n,m,d,m);
  else
    A=0;
  end

  if any(hx(:)~=0)
    A=A+reshape(hx,n,p,m);
  end

  if any(A(:)~=0)
    A=fx+arraymult(fz,A,n,m,p,m);
    ind=1:n;
    for i=1:m
      for j=1:m
        J(ind+(i-1)*n,ind+(j-1)*n)=diag(A(:,i,j))*P1;
      end
    end
  end

  if any(hx1(:)~=0)
    P2=funbas(fspace,snext);
    A=arraymult(fz,hx1,n,m,p,m);
    ind=1:n;
    for i=1:m
      for j=1:m
        J(ind+(i-1)*n,ind+(j-1)*n)=J(ind+(i-1)*n,ind+(j-1)*n)+diag(A(:,i,j))*P2;
      end
    end
  end
return