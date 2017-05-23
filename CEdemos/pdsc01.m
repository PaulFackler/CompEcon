function out1=pdsc01(flag,s,x,Vs,mu,rho,sigma)

switch flag
case 'x'
  out1=(Vs(:,1)<1)+1;
case 'f'
  out1=s.*(x==1);
case 'g'
  out1=mu+zeros(size(s,1),1);
case 'sigma'
  out1=sigma+zeros(size(s,1),1);
case 'rho'
  out1=rho+zeros(size(s,1),1);
otherwise
  error('invalid flag')
end