function out1=pdsc02(flag,s,x,Vs,alpha,rho,sigma,H)

switch flag
case 'x'
  out1=(Vs(:,1)<1)+1;
case 'f'
  out1=(H*s).*(x==2);
case 'g'
  out1=alpha*s.*(1-s)-H*(x==2).*s;
case 'sigma'
  out1=sigma*s;
case 'rho'
  out1=rho+zeros(size(s,1),1);
otherwise
  error('invalid flag')
end