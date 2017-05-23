% MFIC03 Model function file for storage example
% See DEMIC03 for usage
function out=mfic03(flag,s,mu,sigma,k,P,rho)
switch flag
case 'f'
  out=-k*s;
case 'mu'
  out=mu+zeros(size(s,1),1);
case 'sigma'
  out=sigma+zeros(size(s,1),1);
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=[P*(s(1)-s(2)) P -P 0];
case 'R-'
  out=zeros(1,4);
end
