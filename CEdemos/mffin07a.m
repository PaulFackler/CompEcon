% MFFIN07a Model function file for Asian option pricing demo
function out=mffin07a(flag,S,r,delta,sigma,L,put);
switch flag
case 'rho'
  n=size(S,1);
  out=delta+zeros(n,1);
case 'mu'
  out= [(r-delta)*S(:,1) S(:,1)/L];
case 'sigma'
  out=[sigma*S(:,1) zeros(size(S,1),3)];
case 'delta'
  out=[];
case 'V0'
  if put
    out=max(0,S(:,2)-S(:,1));
  else
    out=max(0,S(:,1)-S(:,2));
  end
end