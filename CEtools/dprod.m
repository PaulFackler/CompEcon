% DPROD    Computes the direct product of two matrices.
% USEAGE:
%   c=DPROD(a,b)
% The direct sum of two matrices with the same number of rows is
% equivalent to computing row-wise tensor (Kronecker) products.
% Thus DPROD is equivalent to:
%   c=[]; for i=1:size(a,1), c=[c;kron(a(i,:),b(i,:))]; end
% or
%   n=size(a,2); p=size(b,2); c=a(:,ones(p,1)*(1:n)).*b(:,(1:p)'*ones(1,n));
%
% See Also: KRON, CDPROD, DSUM.

% Copyright (c) 1997 by Paul L. Fackler

function c=dprod(a,b)

if nargin~=2; error('Two arguments must be passed'); end
if isempty(a), c=b; 
elseif isempty(b), c=a;
else
  [ra,ca]=size(a);
  [rb,cb]=size(b);
  if ra~=rb
    disp([ra rb])
    error('a and b must have the same number of rows')
  end
  if issparse(a) || issparse(b)
    [ia,ja,sa]=find(a);
    [ib,jb,sb]=find(b);
    if ra==1                          % due to quirk in MATLAB
      ia=ia(:); ja=ja(:); sa=sa(:);   % FIND returns a row vector if
      ib=ib(:); jb=jb(:); sb=sb(:);   % its argument is a row vector 
    end
    nza=length(sa);
    nzb=length(sb);
    ii=ia';
    ii=find(ii(ones(nzb,1),:)==ib(:,ones(nza,1)));
    ic=ib(:,ones(nza,1));
    ic=ic(ii);
    jc=cb*(ja'-1);
    jc=jc(ones(nzb,1),:)+jb(:,ones(nza,1));
    jc=jc(ii);
    sc=sb*sa';
    sc=sc(ii);
    c = sparse(ic,jc,sc,ra,ca*cb);
  else
    k=1:ca; k=k(ones(cb,1),:);
    l=(1:cb)'; l=l(:,ones(ca,1));
    c=a(:,k).*b(:,l);
  end
end

