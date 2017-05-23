% Problem definition file for Heston's option pricing model
function out1=pfin03(flag,S,t,r,delta,kappa,m,sigma,rho,K,put);

n=size(S,1);
switch flag
case 'rho'
  out1=r+zeros(n,1);
case 'mu'
  out1= [r-delta-0.5*S(:,2)  kappa*(m-S(:,2))];
case 'sigma'
  out1=sqrt(S(:,2));
  out1=[out1 rho*sigma*out1 zeros(n,1) sqrt(1-rho*rho)*sigma*out1];
case 'V0'
  if put
    out1=max(0,K-exp(S(:,1)));
  else
    out1=max(0,exp(S(:,1))-K);
  end
end