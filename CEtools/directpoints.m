% DirectPoints Obtains the points searched by the DiRect algorithm
% USAGE
%    x=directpoints(j,result);
% INPUTS
%   j :      p-vector of integer indices of points desired
%              set to [] if all points are desired
%   result : the third output from the direct algorithm
% OUTPUTS
%   x      : dxp matrix of points searched be the direct algorithm

function x=directpoints(j,result);
if isempty(j)
  j=1:length(result.F);
end
x=repmat(result.xmid,1,length(j));
parent=result.parent;
len=result.len;
dim=result.dim;
for k=1:length(j)
  i=j(k);
  ii=parent(i); 
  while (ii>0) 
    x(dim(i),k)=x(dim(i),k)+len(i); 
    i=ii; 
    ii=parent(i); 
  end
end
