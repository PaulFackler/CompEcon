% PFIN01 Problem definition file for Cox-Ingersoll-Ross bond pricing example
function out1=pfin01(flag,S,t,kappa,alpha,sigma);

n=size(S,1);
switch flag
case 'rho'
  out1=S;
case 'mu'
  out1= kappa*(alpha-S);;
case 'sigma'
  out1=sigma*sqrt(S);
case 'V0'
  out1=ones(n,1);
end