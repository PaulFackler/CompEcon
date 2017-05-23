% MFIC05 Model function file for cash management example
% See DEMIC05 for usage
function out=mfic05(flag,s,mu,sigma,rho,r,c,C)
switch flag
case 'f'
  out=zeros(size(s,1),1);
case 'mu'
  out=mu+zeros(size(s,1),1);
case 'sigma'
  out=sigma+zeros(size(s,1),1);
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=[(-c-r/rho)*(s(2)-s(1)) (c+r/rho) (-c-r/rho) 0];
case 'R-'
  out=[(r/rho-C)*(s(1)-s(2)) (r/rho-C) (C-r/rho) 0];
end
