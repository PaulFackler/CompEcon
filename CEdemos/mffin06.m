% MFFIN06 Problem definition file for bond option pricing example
function out=mffin06(flag,S,kappa,alpha,sigma,K,put,cB,fspaceB);
n=size(S,1);
switch flag
case 'rho'
  out=S;
case 'mu'
  out= kappa*(alpha-S);;
case 'sigma'
  out=sigma*sqrt(S);
case 'delta'
  out=[];
case 'V0'
  if nargin<6
    out=ones(n,1);
  else
    bondval=funeval(cB,fspaceB,S);
    if put
      out=max(K-bondval,0);
    else
      out=max(bondval-K,0);
    end  
  end
end