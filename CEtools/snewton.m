% SNEWTON  Solves EVLCPs using Qi & Liao's Smoothing Newton Method
% USAGE
%   [z,w,err,it,x] = snewton(M,q,z0);
% INPUTS
%   M                      : nm x n matrix (m nxn matrices stacked vertically)
%   q                      : n x p matrix
%   z0                     : initial guess, n x 1 vector (optional)
%   tbar,gamma,sigma,delta : tuning parameters (optional)
%   tol                    : convergence tolerance (optional)
%   maxit                  : maximum number of iterations (optional)
%   maxsteps               : maximum number of backsteps (optional)
%   showiters              : display results of each iteration (optional)
%   tlevel                 : finite termination checked if t<tlevel (optional)
%   numskip                : finite termination checked every numskip iterations
% OUTPUTS
%   z   : solution vector
%   w   : nxm matrix of slackness variables: w(:,j) = M_jz+q(:,j)
%   err : maximal approximation error or error code (?)
%   it  : iterations required
%   x   : solution index

% Ref.: SIAM J. Matrix Anal. Appl., 21, 45-66, 1999.

% Copyright (c) 2003, Paul L. Fackler
% paul_fackler@ncsu.edu

function [z,w,err,it,x]=snewton(M,q,z0,tbar,gamma,sigma,delta,tol,maxit,maxsteps,showiters,tlevel,numskip)

[n,m]=size(q);
q(q>realmax/2)=realmax/2;  % avoid overflow and NaN problems 

% Default settings
if nargin<3 | isempty(z0), z = zeros(n,1); 
else, z=z0;   
end
if nargin<4  | isempty(tbar),      tbar      = 1;              end
if nargin<5  | isempty(gamma),     gamma     = 0.2;            end
if nargin<6  | isempty(sigma),     sigma     = 0.5e-4;         end
if nargin<7  | isempty(delta),     delta     = 0.5;            end
if nargin<8  | isempty(tol),       tol       = 1e-12;          end
if nargin<9  | isempty(maxit),     maxit     = min(1000,10*n); end
if nargin<10 | isempty(maxsteps),  maxsteps  = 25;             end
if nargin<11 | isempty(showiters), showiters = 1;              end 
if nargin<12 | isempty(tlevel),    tlevel    = 1e-5;           end 
if nargin<13 | isempty(numskip),   numskip   = 1;              end 

if 0
% normalize so the Mi are unit diagonal
if 0
  ii=reshape(repmat((1:n:m*n),n,1),n*m,1)+repmat((0:m*n+1:m*n*n)',m,1); 
  v=abs(M(ii)); 
else
  v=[];
  for i=1:m
    v=[v;diag(M((i-1)*n+(1:n),:))];
  end
end
v(v==0)=1;
M=diagmult(1./v,M);
q(:)=q(:)./v;
end



wbar=[zeros(n,1);tbar];
t=tbar;
[fval]=merit(M,q,z,t);
psi=fval'*fval;

W=psi;
minpsi=repmat(psi,1,5);
minpsi=[psi zeros(1,4)];
minpsii=1;
factor=2*sigma*(1-gamma*tbar);
skip=0;
x=[];
dt=1;

for it=1:maxit
   % evaluate merit function
   dmz=[];
   [fval,dmz,dmt,x]=merit(M,q,z,t);
   % check if exact solution can be obtained
   skip=skip+1;
   if ((t<tlevel | skip>=numskip) & ~isempty(x))
     ind=(x-1)*n+(1:n)';
     dr=zeros(n,m); dr(ind)=1;
     warning off
     if 1
       zz=-(dM2w(M,dr)\q(ind));
     else
       [L1,U1,ri1,ci1]=lusetup(dM2w(M,dr));
       zz = -lusolve(L1,U1,ri1,ci1,q(ind));
       clear L1 U1
     end
     warning on
     %w=reshape(M*zz,n,m)+q;
     w=Mq2w(zz,M,q);
     [r,xx]=mincol(w);
     err=max(abs(x-xx));
     if err<tol & ~any(isnan(zz(:)))
       err=0;
       z=zz;
       return
     end
     skip=0;
   end
   psi=fval'*fval+t*t;
   % Check for singularity problems
   if isnan(psi) | isinf(psi), break; end
   % otherwise check if merit function value is small enough
   if psi<tol, break, end
   if showiters, fprintf('%4i %6.2e %6.2e\n',[it t psi]); end
  
   % update control parameters and determine Newton step
   beta=gamma*min(1,psi);
   checklevel=factor*psi;
   if psi>min(minpsi), W=psi; end
   minpsii=minpsii+1; if minpsii>5, minpsii=1; end
   minpsi(minpsii)=psi;
  
   % update iterates
   dt = beta*tbar-t;
   %dz = dmz\(-fval-dmt*dt);
   [L,U,ri,ci]=lusetup(dmz);
   dz = lusolve(L,U,ri,ci,(-fval-dmt*dt));
   
   for backstep=1:maxsteps 
      psinew = merit(M,q,z+dz,t+dt);
      psinew = psinew'*psinew+t*t;
      if psinew<=W-checklevel, 
         break, 
      end
      dz=dz*delta; dt=dt*delta; checklevel=checklevel*delta;
   end
   z = z+dz;
   t = t+dt;
end

if it>=maxit & nargout<4
 warning('Failure to converge in snewton');
elseif it>=maxit
  err=1;
else
  err=0;
end

%w=reshape(M*z,n,m)+q;
w=Mq2w(z,M,q);
%r=min(w,[],2);

return


% merit function and its derivatives
function [phi,dphiz,dphit,x]=merit(M,q,z,t)
  [n,m]=size(q);
  %w=reshape(M*z,n,m)+q;
  w=Mq2w(z,M,q);
  % the next line is a bottleneck - use MEX file mincol to speed it up
  % mw=min(w,[],2);  
  mw=mincol(w);
  ew=exp((mw(:,ones(1,m))-w)/t);
  sw=sum(ew,2);
  lsw=log(sw);
  phi=mw-t*lsw;

  if nargout>1
    temp=(1:n)';
    lambda=sparse(temp,temp,1./sw,n,n)*ew;
    %dphiz=sparse(temp(:,ones(1,m)),1:n*m,lambda,n,n*m)*M;
    dphiz=dM2w(M,lambda);
    dphit=(mw-sum(w.*lambda,2))/t-lsw;
    [ignore,x]=maxcol(lambda);
  end
