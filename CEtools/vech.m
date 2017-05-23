function v=vech(V,flag);
% VECH Vectorizes the lower triangular part of a square matrix.
%   v=VECH(V);
% If V is n x n, v is n(n+1)/2 x 1.
%   v=VECH(V,1);
% vectorizes the upper triangular part of V.

if nargin<2 | isempty(flag), flag=0; end

n=size(V,1);

if flag==1
  ind=find(triu(ones(n,n)));
else
  ind=find(tril(ones(n,n)));
end
v=V(ind);