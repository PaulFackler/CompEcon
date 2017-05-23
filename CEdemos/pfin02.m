% PFIN02 Problem definition file for Black-Scholes demo
function out1=pfin02(flag,S,t,r,delta,sigma,K,put);

 n=size(S,1);
 switch flag
 case 'rho'
   out1=r+zeros(n,1);
 case 'mu'
   out1= (r-delta-0.5*sigma.^2);
 case 'sigma'
   out1=sigma;
 case 'V0'
   if put
     out1=max(0,K-exp(S));
   else
     out1=max(0,exp(S)-K);
   end
 end