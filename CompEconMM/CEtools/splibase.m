%% SPLIBASE
%
%  Computes polynomial spline basis.
%
%  Usage
%    [B,x] = splibase(breaks,evennum,k,x,order)
%  Input
%    breaks    : user specified breakpoint sequence
%                (default: evenly spaced non-repeated breakpoints)
%    evennum   : non-zero if breakpoints are all even
%    k         : polynomial order of the spline's pieces (default: 3, cubic)
%    x         : vector of the evaluation points (default: k-point averages of breakpoints)
%    order     : the order of differentiation (default: 0)
%                if a vector, SPLIBASE returns a cell array, otherwise a matrix
%  Output
%    B         : k.n basis matrix
%    x         : evaluation points (useful if defaults values are computed) %
%  Note
%    Number of basis functions is n=length(breaks)+k-1

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [B,x] = splibase(breaks,evennum,k,x,order)

if nargin<3, error('At least three parameters must be passed'); end
if nargin<4, x=[]; end
if nargin<5 | isempty(order), order = 0; end

% GET DEFAULTS
if isempty(k), k = 3; end

if isempty(x)
  x = splinode(breaks,evennum,k);
end

% A FEW CHECKS
if k<0
  error(['Incorrect value for spline order (k): ' num2str(k)]);
end
if min(size(breaks))>1
  error('''breaks'' must be a vector');
end
if any(order>=k)
  error('Order of differentiation must be less than k');
end
if size(x,2)>1
  error('x must be a column vector')
end

p = length(breaks);
m = size(x,1);
minorder = min(order);

% Augment the breakpoint sequence
n = length(breaks)+k-1; a=breaks(1); b=breaks(end);
augbreaks = [a(ones(k-minorder,1));breaks(:);b(ones(k-minorder,1))];

% The following lines determine the maximum index of
%   the breakpoints that are less than or equal to x,
%   (if x=b use the index of the next to last breakpoint).
%  [temp,ind]=sort([-inf;breaks(2:end-1);x(:)]);
%  temp=find(ind>=p);
%  j=ind(temp)-(p-1);
%  ind=temp-(1:m)';
%  ind(j)=ind(:)+(k-minorder);    % add k-minorder for augmented sequence
ind = lookup(augbreaks,x,3);

% Recursively determine the values of a k-order basis matrix.
% This is placed in an (m x k+1-order) matrix
bas = zeros(m,k-minorder+1);
bas(:,1) = ones(m,1);
B = cell(length(order),1);
if max(order)>0, D = splidop(breaks,evennum,k,max(order)); end % Derivative op
if minorder<0, I = splidop(breaks,evennum,k,minorder); end     % Integral op
for j=1:k-minorder
  for jj=j:-1:1
    b0 = augbreaks(ind+jj-j);
    b1 = augbreaks(ind+jj);
    temp = bas(:,jj)./(b1-b0);
    bas(:,jj+1) = (x-b0).*temp+bas(:,jj+1);
    bas(:,jj) = (b1-x).*temp;
  end
  % as now contains the order j spline basis
  ii = find((k-j)==order);
  if ~isempty(ii)
    ii = ii(1);
    % Put values in appropriate columns of a sparse matrix
    r = (1:m)'; r=r(:,ones(k-order(ii)+1,1));
    c = (order(ii)-k:0)-(order(ii)-minorder);
    c = c(ones(m,1),:)+ind(:,ones(k-order(ii)+1,1));
    B{ii} = sparse(r,c,bas(:,1:k-order(ii)+1),m,n-order(ii));
    % If needed compute derivative or anti-derivative operator
    if order(ii)>0
      B{ii} = B{ii}*D{order(ii)};
    elseif order(ii)<0
      B{ii} = B{ii}*I{-order(ii)};
    end
  end
end

if length(order)==1, B = B{1}; end


%% SPLIDOP
%
%  Creates differential operator matrices for polynomial splines.
%
%  Usage
%    [D,n,a,b,parms] = splidop(breaks,evennum,k,order)
%  Input
%    breaks    : n.1 vector of breakpoints
%    evennum   : n if breakpoints are even, 0 otherwise%
%    order     : desired order of differentiation (a scalar)
%  Output
%    D is returned as a cell array of sparse matrices w/ abs(order) 
%    elements; to view it use full(D{i})
%  Example
%    g(x)=B^k(x)*c;
%    g'(x)=B^{k-1}(x)*splidop(n,a,b,1)*c;
%    g"(x)=B^{k-1}(x)*splidop(n,a,b,2)*c;
%  Note
%    Integrals are computed with ORDER<1 (anti-derivatives).
%    int_a^x g(x)=B^{k+1}(x)*splidop(n,a,b,-1)
%
% Most users will not need this function. 

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [D,n,a,b,parms] = splidop(breaks,evennum,k,order)

if nargin<1 error('At least three parameters must be specified'), end
if nargin<2, evennum = 0; end
if nargin<3, k = []; end
if nargin<4 | isempty(order); order = 1; end

if isempty(k), k = 3; end

if order>k
  error('Order of derivative operator cannot be larger than k');
end

n = length(breaks)+k-1;
kk = max(k-1,k-order-1);
augbreaks = [breaks(1)+zeros(kk,1);breaks(:);breaks(end)+zeros(kk,1)];

D = cell(abs(order),1);
if order>0                             % derivative
  temp = k./(augbreaks(k+1:n+k-1)-augbreaks(1:n-1));
  D{1} = spdiags([-temp temp],0:1,n-1,n);
  for i=2:order
    temp = (k+1-i)./(augbreaks(k+1:n+k-i)-augbreaks(i:n-1));
    D{i} = spdiags([-temp temp],0:1,n-i,n+1-i)*D{i-1};
  end
elseif order<0                         % anti-derivative (integral)
  temp = (augbreaks(kk+2:kk+n+1)-augbreaks(kk-k+1:kk+n-k))/(k+1);
  D{1} = sparse(tril(ones(n+1,1)*temp',-1));
  for i=-2:-1:order
    temp = (augbreaks(kk+2:kk+n-i)-augbreaks(kk-k+i+2:kk+n-k))/(k-i);
    D{-i} = sparse(tril(ones(n-i,1)*temp',-1))*D{-1-i};
  end
end

if nargout>1
  n = n-order;
  a = breaks(1);
  b = breaks(end);
  parms = {breaks,evennum,k-order};
end