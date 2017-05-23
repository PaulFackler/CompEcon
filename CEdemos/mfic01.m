% MFIC01 Model function file for asset replacement example
% See DEMIC01 for usage
function out=mfic01(flag,s,beta,P,rho)
switch flag
case 'f'
  Q=(beta(3)*s+beta(2)).*s+beta(1);
  out=P*Q;
case 'mu'
  out=ones(size(s,1),1);
case 'sigma'
  out=[];
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=zeros(1,4);
case 'R-'
  out=zeros(1,4);
end
