% lusetup Sets up factors needed for efficient solving of linear equations
% USAGE
%   [L,U,rowindex,colindex]=lusetup(A);
% INPUT
%   A : nxn matrix
% OUTPUTS
%   L,U : nxn effectively lower and upper triangular matrices
%   cind : index for column reordering
%   rind : index for row reordering
%
% To solve Ax=b use
%   [L,U,rowindex,colindex]=lusetup(A);
%   x=lusolve(L,U,b,rowindex,colindex);
%
% See: lusolve

function [L,U,rowindex,colindex]=lusetup(A,type)
n=size(A);
if ndims(n)>2 | n(1)~=n(2)
  error('A must be a square matrix');
end
n=n(1);
if issparse(A)
  if nargin>1 & lower(type(1))=='i'
    colindex=[];
    [L,U,P]=luinc(A,'0');
  elseif 0  % this approach does not seem to give as accurate results
    [L,U,P,Q]=lu(A,1);
    [colindex,i]=find(Q);
  else
    colindex=colamd(A);
    [L,U,P]=lu(A(:,colindex));
  end
else
  colindex=[]; 
  [L,U,P]=lu(A);
end
[rowindex,i]=find(P); rowindex(rowindex)=1:n; 