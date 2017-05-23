% OPTINVEST Solves the optimal investment problem with mean-reverting return
% Finds the optimal investment rule for a project that has return process
%   dS=alpha(m-S)dt + sigma*SdW
% USAGE
%   [e,vstar,V,dV]=OptInvest(v,r,alpha,m,sigma,I,n);
% INPUTS
%   r : interest rate
%   alpha, m, sigma : return process paramters
%   I : fixed investment cost
%   n : number of nodes used for Chebyshev collocation
% OUTPUTS
%   Sstar  : the trigger return level (invest if S>=Sstar)
%   c      : coefficients for value function approximation
%   fspace : function family definition structure

function [Sstar,c,fspace]=optinvest(r,alpha,m,sigma,I,n)
  if nargin<6 | isempty(n),  n=50; end
  fspace=fundefn('cheb',n-1,0,1);
  t=funnode(fspace);                       % nodal points on [0,1] 
  fspace=fundefn('cheb',n,0,1);
  Phi0=funbas(fspace,t);
  Phi1=funbas(fspace,t,1);
  u=ones(1,n);
  tt=t.*t;
  B=(sigma.^2/2)*tt(:,u).*funbas(fspace,t,2) + alpha*m*t(:,u).*Phi1 - r*Phi0;
  B=[B;funbas(fspace,[1])];
  B1=[alpha*tt(:,u).*Phi1;zeros(1,n)];   % basis matrix for VSTAR part 
  phi11=funbas(fspace,1,1);

  Sstar=2*m;
  Sstar=broyden('optinvestr',Sstar,[],I,B,B1,phi11,n);
  [e,c]=optinvestr(Sstar,I,B,B1,phi11,n);