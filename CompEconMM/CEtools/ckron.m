%% CKRON 
%
%  Repeated Kronecker products on a cell array of matrices.
%
%  Usage
%    z = ckron(b)     Solves (B1xB2x...xBd)
%    z = ckron(b,1)   Solves (inv(B1)xinv(B2)x...xinv(Bd))
%    where x denotes Kronecker (tensor) product.
%  Input
%    b         : d.1 cell array
%  Output
%    z         : Kronecker product

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function z = ckron(b,invert)

if nargin<1, error('At least one parameter must be passed'), end
if nargin==1, invert=0; end

[d,m,n] = csize(b);
if invert & any(m~=n)
  error('Matrix elements must be square to invert');
end

if isempty(d)
  if invert
    z = inv(b);
  else
    z=b;
  end
else
  if invert
    z = inv(b{1});
  else
    z = b{1};
  end
  for i=2:d
    if invert
      z = kron(z,inv(b{i}));
    else
      z = kron(z,b{i});
    end
  end
end


%% CSIZE
%
%  Returns dimension information for cell arrays.
%
%  Usage
%   [d,m,n] = csize(b)
%  Input
%    b         : d.1 cell array
%  Output
%   If b is a cell array then:
%     d is the number of matrices in b
%     m is a dx2 matrix of row and column dimensions or m and n are dx1 vectors
%   If b is not a cell array, d=[]; this can be used to test if b is a cell
%   array (isempty(d) is true is b is not a cell array) m = size(b) or
%   m=size(b,1) and n=size(b,2)

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu%

function [d,m,n] = csize(b)

if iscell(b)
  d = length(b);
  m = zeros(d,2);
  for i=1:d
    m(i,:) = size(b{i});
  end
else
  d = [];
  m = size(b);
end

if nargout==0
  disp([(1:d)' m])
elseif nargout==3
  n = m(:,2);
  m = m(:,1);
end