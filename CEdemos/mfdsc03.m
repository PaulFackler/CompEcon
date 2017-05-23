function out=pdsc01(flag,s,x,Vs,r,mu,sigma,C)
 switch flag
 case 'f'
   out=(s-C).*x;
 case 'g'
   out=mu*s;
 case 'sigma'
   out=sigma*s;
 case 'rho'
   out=r+zeros(size(s,1),1);
 otherwise
   error('invalid flag')
 end