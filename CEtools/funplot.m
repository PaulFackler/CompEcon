function funplot(c,fspace,n,order)
if nargin<4, order=0; end
s=nodeunif(n,fspace.a,fspace.b);
f=funeval(c,fspace,s,order);
switch size(s,2)
case 1
  plot(s,f);
case 2
  surf(s{1},s{2},reshape(f,size(s{1},1),size(s{2},1))');
otherwise
  error('not supported for higher than 2 dimensions')
end