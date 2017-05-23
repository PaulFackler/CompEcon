function out1=mfdsc04(flag,s,x,Vs,A,theta,sigma,rho)

n=size(s,1);
switch flag
case 'f'
  out1=s(ones(n,1)+(x-1)*n);
case 'g'
  out1=(theta(ones(n,1),:)-s(:,1:size(A,1)))*A;
case 'sigma'
  out1=sigma(ones(n,1),:);
case 'rho'
  out1=rho+zeros(n,1);
otherwise
  error('invalid flag')
end