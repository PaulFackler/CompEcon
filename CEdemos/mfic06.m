% MFIC06 Model function file for optimal fish harvesting problem (infinite harvest rate)
% See DEMIC06 for usage
function out=mfic06(flag,s,alpha,sigma,P,c,rho)
switch flag
case 'f'
  out=zeros(size(s,1),1);
case 'mu'
  out=alpha*(1-s).*s;
case 'sigma'
  if sigma==0, out=[]; 
  else,  out=sigma*s;
  end
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=zeros(1,4);
case 'R-'
  S0=s(:,1); S1=s(:,2);
  out=[P*(S0-S1)-c*log(S0./S1) P-c./S0 -P+c./S1 c./S0.^2];
otherwise
  error('invalid flag')
end