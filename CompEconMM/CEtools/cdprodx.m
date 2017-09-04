%% CDPRODX
%
%  Iterated direct product of cell array times a matrix. The direct product
%  of a set of matrices with the same number of rows is equivalent to
%  computing row-wise tensor (Kronecker) products:
%   a(i,:) = (B{ind(1)}(i,:) x ... x B{ind(d)}(i,:))*c
%
%  Usage
%    a = cdprodx(B,c,ind)
%  Input
%    B         : p.q cell array containing matrices of dimension m by n(i)
%    c         : prod(n(ind(i))) by k matrix
%    ind       : d.1 vector of indices for matrices to select from cell array 
%                (default: all matrices in B, 1:p*q)
%  Output
%    a         : m.k matrix
%  Options
%    maxit     : maximum number of iterations (100)

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function a = cdprodx(b,c,ind)

if nargin<3
  if nargin<2
    error('Must pass two parameters');
  else
    ind = 1:prod(size(b));
  end
end

if ~iscell(b)
  a = b*c;
else  
  d = length(ind);
  a = b{ind(d)}; 
  for i=d-1:-1:1 
    a = dprod(b{ind(i)},a);
  end
  a = a*c;
end