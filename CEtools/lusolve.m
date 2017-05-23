% LUSOLVE Sparse LU solver: x(colindex)=U\(L\(b(rowindex)));
% USAGE
%   x=lusolve(L,U,b,rowindex,colindex);
% INPUTS
%   L        : lower triangular matrix (n x n) (dense or sparse)
%   U        : upper triangular matrix (n x n) (dense or sparse)
%   b        : dense vector (n x 1)
%   rowindex : optional row permutation index vector (n x 1)
%   colindex : optional column permutation index vector (n x 1)
% OUTPUT
%   x :  dense vector (n x 1)
% 
% Note: no checks are made on the triangularity of L or U. Only the
% lower part of L and the upper part of U are used. 
% Also, no checks are made to ensure that rowindex and colindex are
% valid permuation vectors (i.e., are composed of a permutation 
% of the integers 1,..,n).
% 
% The inputs can be obtained using
%    colindex=colmmd(A);
%    [L,U,P]=lu(A(:,colindex));
%    [rowindex,temp]=find(P');
%
% Coded as a MEX file lusolve.c

% Copyright (c) 1997-2003, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=lusolve(L,U,b,rowindex,colindex)
x=zeros(size(b));
if nargin<5
  if nargin<4 
     x=U\(L\b);
  else
     x=U\(L\(b(rowindex)));
  end
else
  x(colindex)=U\(L\(b(rowindex)));
end
