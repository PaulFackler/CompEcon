% Problem definition file for American option pricing demo
function out1=pfin04(flag,S,t,r,delta,sigma,K,put);

switch flag
case 'rho'
  n=size(S,1);
  out1=r+zeros(n,1);
case 'mu'
  out1= (r-delta)*S;
case 'sigma'
  out1=sigma*S;
case 'V0'
  if put
    out1=max(0,K-S);
  else
    out1=max(0,S-K);
  end
end