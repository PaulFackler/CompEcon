% FUNINFGEN Infinitesimal generator for approximating families of functions
% USAGE
%   [L,B]=funinfgen(mu,sigma,x,fspace,upwind);
% INPUTS
%   mu      : 1xd cell array of drift terms or other alternatives (see below) 
%   sigma   : dxq cell array of diffusion terms or other alternatives (see below) 
%   x       : 1xd cell array of columns vectors (to evaluate at grid points)
%               or mxd matrix of evaluation points           
%   fspace  : function definition structure (see fundef)
%   upwind  : 0/1 variable - set to 1 if using 'lin' family of functions
%               and an upwinding basis is desired (default: 0)
%               Note: this option is only implemented for evenly spaced
%               breakpoints with basis matrices evaluated at the
%               breakpoints (no checks are made, however)
%
% OUTPUTS
%   L       : mxn infinitesimal generator operator
%   B       : mxn basis matrix (optional)
%
% mu and sigma can be input in a variety of forms:
%   cell array of m-vectors (d-element for mu, dxq for sigma)
%   cell array of function handles (d-element for mu and dxq for sigma)
%      that each accept an mxd matrix and return an mx1 vector
%   matrices of values (mxd for mu and mxdxq for sigma)
%   single function handle that accepts an mxd matrix and returns a matrix
%     (mxd for mu and mxdxq for sigma)
%
% Note that q is the number of Brownian motions driving the system. This
% can be less than d (the number of variables).

function [varargout]=funinfgen(mu,sigma,x,fspace,upwind)
if nargin<5 || isempty(upwind)
  upwind=0;
end
if nargout>2
  error('Only 2 outputs may be requested')
end

varargout=cell(max(1,nargout),1);
if isfield(fspace,'funtype')
  error('not implemented')
  switch fspace.funtype
    case 'csrbf'
      if iscell(x), x=gridmake(x); end
      [varargout{:}]=csrbfinfgen(mu,sigma,x,fspace.parms{:});
    case 'rbf'
      if iscell(x), x=gridmake(x); end
      [varargout{:}]=rbfinfgen(mu,sigma,x,fspace.parms{:});
    otherwise
      error('not implemented for this basis type')
  end
else
  switch fspace.bastype{1}    
  case 'lin'
    if upwind      
      [varargout{:}]=infgenx(mu,sigma,x,fspace,1);
    else
      [varargout{:}]=infgenx(mu,sigma,x,fspace,0);
    end
  otherwise
    [varargout{:}]=infgenx(mu,sigma,x,fspace,0);
  end
end
    

% INFGENX Infinitesimal generator for general functions
% USAGE
%   same as for main function

function [L,B]=infgenx(mu,sigma,x,fspace,upwind)
if iscell(x)
  d=length(x);
  m=1; for i=1:d, m=m*length(x{i}); end
else
  [m,d]=size(x);
end
if isnumeric(sigma)
  q=size(sigma,3);
else
  q=size(sigma,2);
end
n=prod(fspace.n);

% check whether sigma is empty
order2=0;
if iscell(sigma)
  for i=1:d, for j=1:q, if ~isempty(sigma{i,j}), order2=1; break; end; end; end
else
  if ~isempty(sigma), order2=1; end
end

if order2
  B=funbasx(fspace,x,[0;1;2]);
else
  B=funbasx(fspace,x,[0;1]);
end

if issparse(B.vals{1})
  L=sparse(1,1); % sparse scalar 0
else
  L=zeros(m,n);
end

S=gridmake(x);

if isa(mu,   'function_handle'), mu    = mu(S);    end
if isa(sigma,'function_handle'), sigma = sigma(S); end

if upwind
   Sigma=sigma2Sigma(sigma,gridmake(funnode(fspace)));
   L=infgen(mu,Sigma,fspace.n,(fspace.b-fspace.a)./(fspace.n-1));
   B=speye(size(L,1)); 
   return
end

% compute the first order part
if ~isempty(mu)
if iscell(mu)
  for i=1:d
    if ~isempty(mu{i});
      ind=zeros(1,d); ind(i)=1;
      if isnumeric(mu{i})
        L=L+diagmult(mu{i},basget(B,ind));
      else        
        L=L+diagmult(mu{i}(S),basget(B,ind));
      end
    end
  end
else
  for i=1:d
    ind=zeros(1,d); ind(i)=1;
    L=L+diagmult(mu(:,i),basget(B,ind));
  end
end
end
% compute the second order part (if needed)
if order2
  if iscell(sigma)
    for i=1:d
      Sij=zeros(m,1);
      sij=0;
      for k=1:q
        if ~isempty(sigma{i,k}) 
          if isnumeric(sigma{i,k})
            Sij=Sij+sigma{i,k}.^2; 
          else
            Sij=Sij+sigma{i,k}(S).^2; 
          end
          sij=1; 
        end
      end
      if sij, 
        ind=zeros(1,d); ind(i)=2;
        L=L+diagmult(Sij/2,basget(B,ind)); 
      end
      for j=i+1:d
        Sij=zeros(m,1);
        sij=0;
        for k=1:q
          if ~isempty(sigma{i,k}) && ~isempty(sigma{j,k}), 
            if isnumeric(sigma{i,k})
              Sij=Sij+sigma{i,k}.*sigma{j,k}; 
            else              
              Sij=Sij+sigma{i,k}(S).*sigma{j,k}(S);
            end
            sij=1; 
          end
        end
        if sij, 
          ind=zeros(1,d); ind(i)=1; ind(j)=1;
          L=L+diagmult(Sij,basget(B,ind)); 
        end
      end
    end  
  else
    for i=1:d
      Sij=zeros(m,1);
      for k=1:q 
        Sij=Sij+sigma(:,i,k).^2; 
      end
      ind=zeros(1,d); ind(i)=2; 
      L=L+diagmult(Sij/2,basget(B,ind));
      for j=i+1:d
        Sij=zeros(m,1);
        for k=1:q 
          Sij=Sij+sigma(:,i,k).*sigma(:,j,k); 
        end
        ind=zeros(1,d); ind(i)=1; ind(j)=1;
        L=L+diagmult(Sij,basget(B,ind));
      end
    end
  end
end
if nargout>1
  B=ckron(B.vals(1,d:-1:1));
end

function b=basget(B,ind)
d=length(ind);
if strcmp(B.format,'tensor')
  b=B.vals{ind(1)+1,1};
  for i=2:d
    b=kronmex(B.vals{ind(i)+1,i},b);
  end
else
  b=B.vals{ind(1)+1,1};
  for i=2:d
    b=dprod(B.vals{ind(i)+1,i},b);
  end
end
  

% infgen Creates upwind finite difference matrix infinitesimal generator operator
% USAGE
%   L=infgen(mu,Sigma,n,h);
% INPUTS
%  mu    : 1xd cell arrays composed on Nx1 vectors for the order 1 terms
%             Alternatively mu can be Nxd
%  Sigma : dxd cell array composed of Nx1 vectors for the order 2 terms
%             (only the lower triangular portion is used) 
%             Alternatively Sigma can be Nxdxd
%  n     : d-vector of number of nodal values for each variable (N=prod(n))
%  h     : d-vector of step sizes for each variable
% OUTPUT
%  L   : NxN matrix. L*v approximates the operator
%            V'(s)mu(s) + 0.5 vec(V''(S))'vec(Sigma(s))
%        where v(ij)=V(s(ij)) and s(i,j)=s0(j)+(i-1)*h(j), i=1,...,n(j)
%        (ij) is a single index and (i,j) is a double index
%           see SUB2IND and IND2SUB for discussion)
%
% This is a discrete approximation of the infinitesimal generator for the stochastic process 
%   dS=mu(S)dt + sigma(S)dW
% where Sigma(S)=sigma(S)*sigma(S)' (use sigma2Sigma to obtain Sigma from sigma)
%
% Nodal values must be evenly spaced.
function L=infgen(mu,Sigma,n,h)
  d=length(n);  % problem dimension
  N=prod(n);    % problem size
  
  % all dimensions are the same so less storage is required
  if all(n==n(1) & h==h(1)) 
    ni=n(1);
    hi=h(1);
    Phi=repmat({speye(ni)},1,d);
    Df=repmat({forward(ni,hi)},1,d);   % forward difference operator
    Db=repmat({backward(ni,hi)},1,d);  % backward difference operator           
  % not all dimensions are the same - store all dimensions separately
  else
    Phi=cell(1,d);
    Df=cell(1,d);
    Db=cell(1,d);
    for i=1:d
      ni=n(i);
      hi=h(i);
      di=d-i+1;
      Phi{di}=speye(ni);
      Df{di}=forward(ni,hi);   % forward difference operator
      Db{di}=backward(ni,hi);  % backward difference operator
    end
  end
 
  L=sparse(N,N);
  
  % Order 1 terms
  Di=Phi; 
  for i=1:d
    if isa(mu,'cell'), mui=mu{i};
    else               mui=mu(:,i);
    end
    if ~isempty(mui)
      di=d-i+1;
      % positive drift rates
      ii=find(mui>0);
      if ~isempty(ii)
        Di{di}=Df{di}; Li=ckron(Di);
        temp=zeros(size(mui));
        temp(ii)=mui(ii);
        L=L+diagmult(temp,Li);
      end
      % negative drift rates
      ii=find(mui<0);
      if ~isempty(ii)
        Di{di}=Db{di}; Li=ckron(Di);
        temp=zeros(size(mui));
        temp(ii)=mui(ii);
        L=L+diagmult(temp,Li);
      end
      Di{di}=Phi{di};
    end
  end

  % Order 2 terms
  SigmaDiag=false;
  if ~isa(Sigma,'cell')
    n2=size(Sigma);
    if length(n2)==2
      if n2(2)==d
        SigmaDiag=true;
      else
        Sigma=reshape(Sigma,N,d,d);
      end
    end
  else
    n2=numel(Sigma);
    if n2==d
      SigmaDiag=true;
    else
      Sigma=reshape(Sigma,d,d);
    end
  end
  for j=1:d
    if SigmaDiag
      if isa(Sigma,'cell')
        Sigmaij=Sigma{j};
      else
        Sigmaij=Sigma(:,j);
      end
    else
      if isa(Sigma,'cell')
        Sigmaij=Sigma{j,j};
      else
        Sigmaij=Sigma(:,j,j);
      end
    end
    dj=d-j+1;
    % get jth second derivative term
    if ~isempty(Sigmaij)
      Di{dj}=Df{dj}*Db{dj}; 
      Di{dj}(1,[1 2 3 4]) = [2 -5 4 -1]/h(j)^2;  % correction for left endpoint
      Li=ckron(Di);
      L=L+diagmult(Sigmaij/2,Li);
    end
    if ~SigmaDiag
      for i=j+1:d
        if isa(Sigma,'cell')
          Sigmaij=Sigma{i,j};
        else
          Sigmaij=Sigma(:,i,j);
        end
        % get ijth cross derivative term
        if ~isempty(Sigmaij)
          di=d-i+1;
          % negative cross terms
          ii=find(Sigmaij<0);
          if ~isempty(ii)
            Di{dj}=Df{dj}; Di{di}=Db{di}; Li=ckron(Di);
            Di{dj}=Db{dj}; Di{di}=Df{di}; Li=Li+ckron(Di);
            temp=zeros(size(Sigmaij));
            temp(ii)=Sigmaij(ii)/2;
            L=L+diagmult(temp,Li);
          end
          % positive cross terms
          ii=find(Sigmaij>0);
          if ~isempty(ii)
            Di{dj}=Df{dj}; Di{di}=Df{di}; Li=ckron(Di);
            Di{dj}=Db{dj}; Di{di}=Db{di}; Li=Li+ckron(Di);
            temp=zeros(size(Sigmaij));
            temp(ii)=Sigmaij(ii)/2;
            L=L+diagmult(temp,Li);
          end
          Di{di}=Phi{di};
        end
      end
    end
    Di{dj}=Phi{dj};
  end
  
  function D=forward(n,h)
    D=sparse([1:n-1 1:n-1 n n n],[2:n 1:n-1  n-2 n-1 n],...
                 [ones(1,n-1) -ones(1,n-1) 1 -3 2]/h,n,n);
                    
  function D=backward(n,h)
    D=sparse([1 1 1 2:n 2:n],[1 2 3 2:n 1:n-1],...
                 [-2 3 -1 ones(1,n-1) -ones(1,n-1)]/h,n,n);
      
% sigma2Sigma Converts sigma to Sigma=sigma*sigma'
% USAGE
%   Sigma=sigma2Sigma(sigma,d,sigmadiag);
% INPUTS
%   sigma      nxdxq matrix or dxd cell array of n-vectors
%   d          scalar dimension 
%   sigmadiag  1 if sigma is diagonal (optional: default=0)
% OUTPUT
%   Sigma      dxd cell array of n-vectors
% Used to convert sigma(S) to Sigma(S)=sigma(s)*sigma(S)'
% where S is an n-vector
function Sigma=sigma2Sigma(sigma,S,sigmadiag)
d=size(S,2);
% sigma is not diagonal
if nargin<3 || sigmadiag==0
  if isa(sigma,'cell')
    [d,q]=size(sigma);
    Sigma=cell(d,d);
    for j=1:d
      for i=1:j
        Sigma{i,j}=sparse(1,1);
        empty=true;
        for k=1:q
          if ~isempty(sigma{i,k}) && ~isempty(sigma{j,k})
            if isnumeric(sigma{i,k})
              Sigma{i,j}=Sigma{i,j}+sigma{i,k}.*sigma{j,k};
            else
              Sigma{i,j}=Sigma{i,j}+sigma{i,k}(S).*sigma{j,k}(S);
            end
            empty=false;
          end
        end
        if empty, Sigma{j,i}=[]; end
        Sigma{j,i}=Sigma{i,j};
      end
    end
  else
    n=size(sigma);
    if  length(n)==2, 
      q=n(2)/d;
      sigma=reshape(sigma,n(1),d,q); 
    else
      q=n(3);
    end
    Sigma=cell(d,d);
    for j=1:d
      for i=1:j
        Sigma{i,j}=sigma(:,i,1).*sigma(:,j,1);
        for k=2:q
          Sigma{i,j}=Sigma{i,j}+sigma(:,i,k).*sigma(:,j,k);
        end
        Sigma{j,i}=Sigma{i,j};
      end
    end
  end
else
  Sigma=cell(d,d);
  if isa(sigma,'cell')
    for i=1:d, Sigma{i,i}=sigma{i,i}.^2; end
  else
    ns=size(sigma);
    if length(ns)==2, sigma=reshape(sigma,ns(1),d,d); end
    for i=1:d, Sigma{i,i}=sigma(:,i,i).^2; end
  end
end

% DIAGMULT Computes either diag(a)*b or a*diag(b)
% USAGE
%   c=diagmult(a,b);
% INPUTS
%   a,b : m-vector and mxn matrix
%         or 
%         mxn matrix and n-vector
% OUTPUT
%   c   : mxn matrix 
%
% Coded as a C-MEX function for speed

% Note: C version is not implemented for complex matrices or matrices 
% with data type other than double. The vector must be full but
% the matrix can be sparse or full.

% Copyright (c) 2005-8, Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu
function c=diagmult(a,b)
 [m,n]=size(a);
 % a is a vector
 if m==1 || n==1
     n=length(a);
     if size(b,1)~=n
       error('Inputs are not compatible')
     end
     c=sparse(1:n,1:n,a,n,n)*b; 
 else
   [m,n]=size(b);
   % b is a vector
   if m==1 || n==1
     n=length(b);
     if size(a,2)~=n
       error('Inputs are not compatible')
     end
     c=a*sparse(1:n,1:n,b,n,n);
   else
     error('Either a or b must be vectors')
   end
 end