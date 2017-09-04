%% RSSOLVE
%
%  Solves continuous time regime switching model with one continuous state variable
%
%  Usage
%    [cv,basis,x] = rssolve(model,x,n,type)
%  Input
%    model     : a structure variable (defined below)
%    x         : initial control values (px1)
%    n         : degree of approximation (m-vector)
%    type      : the type of basis functions to use (default='cheb')
%  Output
%    cv        : 1.m cell array of solution coefficients with cv{i} n(i).1
%    basis     : 1.m cell array of function definition structures
%    x         : p.1 optimal switch points
%  Options
%    maxiters  : maximum number of iterations
%    tol       : convergence tolerance
%    showiters : 0/1, 1 to display iteration results
%  Model Structure
%    The model structure should contain the following fields:
%      func     : model function file - see below  
%      params   : additional parameters to pass to model function file
%      xindex   : p.5 defines the nature of the control values - see below 
%     Regime i has a no switch region [x(1,i),x(2,i)]. When s hits x(j,i)
%     the discrete regime switches to xindex(j,i) If x(j,i) is a boundary
%     for s x(j,i)=i and no switch occurs xindex defines the topology of
%     switches. If there is a point at which i switches to j, set
%      xindex(k,1)=i and xindex(k,2)=j 
%     and define the type of side constraints to hold at that point:
%      xindex(k,3)=>0 implies   V(s,i) =   V(s,j)+reward(k,1)
%      xindex(k,4)=>0 implies  V'(s,i) =  V'(s,j)+reward(k,2)
%      xindex(k,5)=>0 implies V''(s,i) = V''(s,j)+reward(k,3)
%    where reward is returned by the model function file. To impose a side
%    condition on V(s,i) alone, set 
%     xindex(k,1)=i and xindex(k,2)=0 
%    and define the type of side constraints to hold at that point:
%     xindex(k,3)=>0 implies   V(s,i) = reward(k,1)
%     xindex(k,4)=>0 implies  V'(s,i) = reward(k,2)
%     xindex(k,5)=>0 implies V''(s,i) = reward(k,3)
%  Model Function File
%    The model function file should have the format
%      out = func(flag,s,additional parameters)
%        switch flag
%          case 'f'
%            Return a matrix f representing the reward function
%          case 'g'
%            Return a matrix g representing the drift function
%            of the state transition equation
%          case 'sigma'
%            Return a matrix sigma representing the diffusion function
%            of the state transition equation
%          case 'rho'
%            Return a vector representing the (state contingent)
%            discount rates
%          case 'reward'
%            Return a k x 3 matrix of jump rewards and derivatives
%            when passed a k x 1 matrix of state values
%        end

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [cv,basis,x] = rssolve(model,x,n,type)

% Set default options
maxiters      = optget('rssolve','maxiters',100);
tol           = optget('rssolve','tol',sqrt(eps));
showiters     = optget('rssolve','showiters',1);

optset('broyden','maxiters',maxiters)
optset('broyden','tol',tol)
optset('broyden','showiters',showiters)

if ~exist('type','var'), type='cheb'; end

% Unpack model variables
func=model.func;
params=model.params;
xindex=model.xindex;

if any(sum(xindex(:,3:5)==1,2)>2 |...
    sum(xindex(:,3:5)==2,2)>1)
  error('Incorrect specification of xindex');
end

m=max(max(xindex(:,[1 2])));
if length(n)==1, n=n(ones(m,1)); end

% Determine number of nodal points
ncond=zeros(1,m);
for k=1:size(xindex,1)
  i=xindex(k,1);
  j=xindex(k,2);
  s=sum(xindex(k,3:5)==1);
  if j==0
    ncond(i)=ncond(i)+s;
  else
    ncond(i)=ncond(i)+1;
    if s==2, ncond(j)=ncond(j)+1; end
  end
end

% Get basis matrices and nodal points
for i=1:m
  basis{i}=fundefn(type,n(i)-ncond(i),0,1);
  z{i}=funnode(basis{i});
  basis{i}=fundefn(type,n(i),0,1);
  Phi0{i}=funbase(basis{i},z{i},0);
  Phi1{i}=funbase(basis{i},z{i},1);
  Phi2{i}=funbase(basis{i},z{i},2);
  phi0{i}=funbase(basis{i},[0;1],0);
  phi1{i}=funbase(basis{i},[0;1],1);
  phi2{i}=funbase(basis{i},[0;1],2);
end

xchoice=logical(sum(xindex(:,3:5)==2,2));
y=x(xchoice);
y=broyden(@rsres,y,x,func,params,xindex,xchoice,m,n,...
  basis,z,ncond,Phi0,Phi1,Phi2,phi0,phi1,phi2);
x(xchoice)=y;
[~,c]=rsres(y,x,func,params,xindex,xchoice,m,n,...
  basis,z,ncond,Phi0,Phi1,Phi2,phi0,phi1,phi2);

for i=1:m
  cv{i}=c(1:n(i));
  c(1:n(i))=[];
  j=xindex(:,1)==i |xindex(:,2)==i;
  a=min(x(j));
  b=max(x(j));
  basis{i}=fundefn(type,n(i),a,b);
end

optset('broyden','defaults')
  

function [e,c]=rsres(y,x,func,params,xindex,xchoice,m,n,...
  basis,z,ncond,Phi0,Phi1,Phi2,~,~,~)

nm=sum(n);
x(xchoice(:))=y;
for i=1:m
  j=xindex(:,1)==i | xindex(:,2)==i;
  a(i)=min(x(j));
  b(i)=max(x(j));
end
r=b-a;

% loop over the regimes
H=zeros(nm,nm);
h=zeros(nm,1);
rows=0;
cols=0;
for i=1:m
  s=z{i}*r(i)+a(i);
  nz=n(i)-ncond(i);
  rows=rows(end)+1:rows(end)+nz;
  cols=cols(end)+1:cols(end)+n(i);
  rho=feval(func,'rho',s,i,params{:});
  mu=feval(func,'g',s,i,params{:});
  sigma=feval(func,'sigma',s,i,params{:});
  if ~isempty(sigma)
    H(rows,cols)=spdiags(rho,0,nz,nz)*Phi0{i}...
      - spdiags(mu/r(i),0,nz,nz)*Phi1{i}...
      - spdiags(sigma.*sigma/(2*r(i)^2),0,nz,nz)*Phi2{i};
  else
    H(rows,cols)=spdiags(rho,0,nz,nz)*Phi0{i}...
      -spdiags(mu/r(i),0,nz,nz)*Phi1{i};
  end
  h(rows)=feval(func,'f',s,i,params{:});
end

% Loop over the switch/boundary points
R=feval(func,'reward',x,[],params{:});
G=[]; g=[];
row=rows(end);
for k=1:size(x,1)
  i=xindex(k,1);
  zi=(x(k)-a(i))/r(i);
  icols=sum(n(1:i-1))+(1:n(i));
  j=xindex(k,2);
  if j~=0
    zj=(x(k)-a(j))/r(j);
    jcols=sum(n(1:j-1))+(1:n(j));
  end
  for order=0:2
    if xindex(k,3+order)==1
      row=row+1;
      H(row,icols)=funbase(basis{i},zi,order)/r(i)^order;
      if j~=0
        H(row,jcols)=-funbase(basis{j},zj,order)/r(j)^order;
      end
      h(row)=R(k,order+1);
    end
    if xindex(k,3+order)==2
      temp=zeros(1,nm);
      temp(1,icols)=funbase(basis{i},zi,order)/r(i)^order;
      if j~=0
        temp(1,jcols)=-funbase(basis{j},zj,order)/r(j)^order;
      end
      G=[G;temp];
      g=[g;R(k,order+1)];
    end
  end
end

c=H\h;
e=G*c-g;

return