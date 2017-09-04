%% CHEBBASE 
%
%  Computes basis matrices for Chebyshev polynomials
%
%  Usage
%    [B,x] = chebbase(n,a,b,x,order)
%  Input
%    n         : number of basis functions (1 plus the polynomial order)
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval

%    x         : k.1 -vector of the evaluation points 
%                (default: roots of order n Chebyshev polynomial)
%    order     : order of differentiation (default: 0)
%                if vector, cell array returned, otherwise matrix returned
%  Output
%    B         : k.n basis matrix or cell array of k.n basis matrices
%    x         : evaluation points (useful if defaults values are computed)

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function [B,x] = chebbase(n,a,b,x,order)

if nargin<3, error('3 parameters must be specified'), end
if nargin<4, x=[]; end
if nargin<5 || isempty(order), order=0; end

minorder=min(0,min(order));

if isempty(x);
  x = chebnode(n,a,b);
  nodetype = optget('chebnode','nodetype');
else
  nodetype = 1;
end

% Compute order 0 basis
if nodetype==0      % evaluate at standard nodes
  temp = ((n-0.5):-1:0.5)';
  bas = cos((pi/n)*temp*(0:(n-1-minorder)));
else                % evaluate at arbitrary nodes
  bas = chebbasex(n-minorder,a,b,x);
end

% Get bases for other orders
if length(order)==1
  if order~=0
    D = chebdop(n,a,b,order);
    B = bas(:,1:n-order)*D{abs(order)};
  else
    B = bas;
  end
else
  B = cell(length(order),1);
  maxorder = max(order);
  if maxorder>0, D = chebdop(n,a,b,maxorder); end
  if minorder<0, I = chebdop(n,a,b,minorder); end
  for ii=1:length(order)
    if order(ii)==0
      B{ii} = bas(:,1:n);
    elseif order(ii)>0
      B{ii} = bas(:,1:n-order(ii))*D{order(ii)};
    else
      B{ii} = bas(:,1:n-order(ii))*I{-order(ii)};
    end
  end
end


%% CHEBDOP
%
%  Creates differential operator matrices for Chebychev polynomial basis.
%
%  Usage
%    [D,n,a,b,parms] = chebdop(n,a,b,order)  
%  Input
%    n         : number of basis functions (1 plus the polynomial order)
%    a         : lower bound on interval
%    b         : upper bound on interval
%    order     : order of differentiation (default: 0)
%                if vector, cell array returned, otherwise matrix returned
%  Output
%    For order 1, element i,j of this matrix is equal to 2(j-1) if i<j and
%    i+j is odd; for the interval [-1,1].  For other intervals divide by
%    (b-a)/2. Higher order derivatives operators are formed by taking
%    multiple products.
%    Example:
%      g(x)=B^{n}(x)*c;
%      g'(x)=B^{n-1}(x)*chebdop(n,a,b)*c;
%      g"(x)=B^{n-2}(x)*chebdop(n,a,b,2)*c;
%    Negative values of ORDER produce an integral operator normalized to be
%    zero at the lower bound a.
%    D is returned as a cell array of size abs(order) composed of sparse matrices.
%    Orders 1 and 2 are stored as a global cell array to avoid recomputation
%    This speeds up repeated calculations involving the first two
%    derivatives.

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function [D,n,a,b,parms] = chebdop(n,a,b,order)

global cheb_dop

if nargin<3 error('3 parameters must be specified'), end
if nargin<4; order=1; end

if order>0
  % Use previously stored values for order=1 and order=2
  % this speeds up calculations considerably with repeated calls
  if ~isempty(cheb_dop) & order<=length(cheb_dop) & size(cheb_dop{1},2)==n
    D=cell(order,1);
    for ii=1:order
      D{ii}=cheb_dop{ii}*(1/((b-a).^ii));
    end
  else
    D=cell(max(2,order),1);
    j = 1:n; j = j(ones(n,1),:); i = j';
    [r,c] = find(rem(i+j,2)==1 & j>i);
    d = sparse(r,c,(4/(b-a))*(j(1,c)-1),n-1,n);
    d(1,:) = d(1,:)/2;
    D{1}=d;
    for ii=2:max(2,order)
      D{ii}=d(1:n-ii,1:n-ii+1)*D{ii-1};
    end
    % store values for order=1 and order=2
    cheb_dop{1}=(b-a)*D{1};
    cheb_dop{2}=((b-a).^2)*D{2};
    if order==1, D(2)=[]; end
  end
elseif order<0
  D=cell(abs(order),1);
  nn=n-order;
  i=(0.25*(b-a))./(1:nn);
  d=sparse([1:nn 1:nn-2],[1:nn 3:nn],[i -i(1:nn-2)],nn,nn);
  d(1,1)=2*d(1,1);
  d0=((-1).^(0:nn-1)).*sum(d);
  D{1}=[d0(1:n);d(1:n,1:n)];
  for ii=-2:-1:order;
    D{-ii}=[d0(1:n-ii-1);d(1:n-ii-1,1:n-ii-1)]*D{-ii-1};
  end
else
  D=speye(n);
end

if nargout>1
  n=n-order;
  parms={n,a,b};
end


%% CHEBBASEX
%
%  Utility used by CHEBBASE 

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function bas = chebbasex(n,a,b,x)
z = (2/(b-a))*(x-(a+b)/2);
m=size(z,1);
bas = zeros(m,n);
bas(:,1) = ones(m,1);
bas(:,2) = z;
z=2*z;
for i=3:n
  bas(:,i) = z.*bas(:,i-1)-bas(:,i-2);
end