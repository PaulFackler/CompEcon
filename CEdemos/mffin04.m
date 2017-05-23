% MFFIN04 Model function file for American option pricing demo
function out=mffin04(flag,S,r,delta,sigma,K,put);
switch flag
case 'rho'
  n=size(S,1);
  out=r+zeros(n,1);
case 'mu'
  out= (r-delta)*S;
case 'sigma'
  out=sigma*S;
case 'delta'
  out=[];
case 'V0'
  if put
    out=max(0,K-S);
  else
    out=max(0,S-K);
  end
end