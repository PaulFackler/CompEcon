% MFFIN07 Model function file for Asian option pricing demo
function out=mffin07(flag,S,r,delta,sigma,L,put);
switch flag
case 'rho'
  n=size(S,1);
  out=delta+zeros(n,1);
case 'mu'
  out= 1-(r-delta)*S+(sigma^2/2)*S;
case 'sigma'
  out=sigma*S;
case 'delta'
  out=[];
case 'V0'
  if put
    out=max(0,S/L-1);
  else
    out=max(0,1-S/L);
  end
end