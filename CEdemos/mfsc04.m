% MFSC04 Model file for optimal fish harvesting problem
function out=mfsc04(flag,s,x,Vs,alpha,sigma,H,P,c,rho)
switch flag
case 'x'
  out=H*(Vs<(P-c./s));
case 'f'
  out=(P-c./s).*s.*x;
case 'g'
  out=(alpha*(1-s)-x).*s;
case 'sigma'
  out=sigma*s;
case 'rho'
  out=rho+zeros(size(s,1),1);
otherwise
  error('invalid flag')
end