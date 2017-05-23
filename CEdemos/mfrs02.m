% MFRS03 Model file for fish harvesting problem
function out=mfrs03(flag,s,x,alpha,sigma,H,P,c,rho)

switch flag
case 'f'
  out=(P-c./s).*s.*(x==2)*H;
case 'g'
  out=(alpha*(1-s)-(x==2)*H).*s;
case 'sigma'
  out=sigma.*s;  %+zeros(size(s,1),1);
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'reward'
  out=[0 0 0;0 0 0;0 0 0];
otherwise
  error('invalid flag')
end