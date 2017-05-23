% MFFIN02 Model function file for Black-Scholes demo
function out=mffin02(flag,S,r,deltaS,sigma,K,put);
switch flag
case 'rho' 
  out = r+zeros(size(S,1),1);
case 'mu'
  out = (r-deltaS)*S;
case 'sigma'
  out = sigma*S;
case 'delta'
 out = [];
case 'V0'
  if put
    out = max(0,K-S);
  else
    out = max(0,S-K);
  end
end