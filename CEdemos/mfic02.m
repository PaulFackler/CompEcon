% MFIC02 Model function file for timber harvesting example
% See DEMIC02 for usage
function out=mfic02(flag,s,alpha,m,sigma,P,rho)
switch flag
case 'f'
  out=zeros(size(s,1),1);
case 'mu'
  out=alpha*(m-s);
case 'sigma'
  out=sigma*sqrt(s);
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=zeros(1,4);
case 'R-'
  out=[P*(s(1)-s(2)) P -P 0];
end
