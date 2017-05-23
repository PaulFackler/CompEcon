function out=mfdsc01(flag,s,x,Vs,c,mu,sigma,rho)
 switch flag
 case 'x'
   out=Vs>0;
 case 'f'
   out=(s-c).*x;
 case 'g'
   out=mu*s;
 case 'sigma'
   out=sigma*s;
 case 'rho'
   out=rho+zeros(size(s,1),1);
 otherwise
   error('invalid flag')
 end