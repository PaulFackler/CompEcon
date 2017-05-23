function out=mfrs01(flag,s,x,c,mu,sigma,rho)
 switch flag
 case 'f'
   out=(s-c);
 case {'g','mu'}
   out=mu*s;
 case 'sigma'
   out=sigma*s;
 case 'rho'
   out=rho+zeros(size(s,1),1);
 case 'reward'
   out=[0 0 0;s(2)/(rho-mu)-c/rho 1/(rho-mu) 0];
 otherwise
   error('invalid flag')
 end