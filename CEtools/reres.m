% RERES Residual file used by RESOLVE with Broyden's method
% See RESOLVE

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [r,x,z,cx]=reres(c,s,x,fspace,func,params,e,w,explicit,noxnext,...
         d,m,p,n,K,Phi,arbmaxit,arbtol,lcpmethod,expapprox,lowmemory); 
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
         d,m,p,n,K,Phi,[],arbmaxit,arbtol,lcpmethod,lowmemory);
  r=r-z;
else
  c=reshape(c,n,m);
  if isempty(Phi),  r=funeval(c,fspace,s);
  else,             r=funeval(c,fspace,Phi);
  end
  n=size(s,1);
  [z,x] = Ehr(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
              d,m,p,n,K,Phi,[],arbmaxit,arbtol,lcpmethod,lowmemory);
  if explicit
    r=r-x;
  else
    r=feval(func,'f',s,x,z,[],[],[],params{:});
  end  
  cx=[];
end
r=r(:);
return

  

% Ehe Computes the expectation at arbitrary s and c
function [z,x,cx]=Ehe(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
                       d,m,p,n,K,Phi,cx,arbmaxit,arbtol,lcpmethod,lowmemory)
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

% Ehr Computes the expectation at arbitrary s and c
function [z,x]=Ehr(c,s,r,fspace,func,params,e,w,explicit,noxnext,...
                   d,m,p,n,K,Phi,cx,arbmaxit,arbtol,lcpmethod,lowmemory)
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

