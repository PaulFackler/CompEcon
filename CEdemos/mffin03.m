% MFFIN03 Model function file for Heston's option pricing model
function out=mffin03(flag,S,r,delta,kappa,m,sigma,rho,K,put);
n=size(S,1);
switch flag
case 'rho'
  out=r+zeros(n,1);
case 'mu'
  out= [r-delta-0.5*S(:,2)  kappa*(m-S(:,2))];
case 'sigma'
  out=sqrt(S(:,2));
  out=[out rho*sigma*out zeros(n,1) sqrt(1-rho*rho)*sigma*out];
case 'delta'
  out=[];
case 'V0'
  if put
    out=max(0,K-exp(S(:,1)));
  else
    out=max(0,exp(S(:,1))-K);
  end
end