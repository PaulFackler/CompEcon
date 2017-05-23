function out=mfrs03(flag,s,x,r,mu,sigma,C,E,I)
 switch flag
 case 'f'
   out=(s-C).*(x==2);
 case 'g'
   out=mu*s;
 case 'sigma'
   out=sigma*s;
 case 'rho'
   out=r+zeros(size(s,1),1);
 case 'reward'
   out=[0 0 0;-I 0 0;-E 0 0;0 0 0];
 otherwise
   error('invalid flag')
 end