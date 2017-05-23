% FFB04 Function file for entry/exit problem 
% See DEMFB04
function [r,m,s]=ffb04(P,rho,mu,sigma)
  n=size(P,1);
  r=rho+zeros(n,1);
  m=mu*P;
  s=sigma*P;