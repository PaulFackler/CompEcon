% CIRBOND Solves the CIR bond pricing model
% The CIR model assumes that the instantaneous interest rate follows
%    dr = kappa(alpha-r)dt+sigma*sqrt(r)dW.
% USAGE
%   c=cirbond(fspace,tau,kappa,alpha,sigma);
% INPUTS
%   fspace : a function family definition structure
%   tau    : the time-to-maturity of the bond (a scalar)
%   kappa,alpha,sigma : model parameters
% OUTPUT
%   c      : coefficient vector
% 
% The bond price can be evaluated using funeval(c,fspace,tau).

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function c=cirbond(fspace,tau,kappa,alpha,sigma) 

% Define nodes and basis
r=funnode(fspace);
Phi0=funbas(fspace,r,0);
Phi1=funbas(fspace,r,1);
Phi2=funbas(fspace,r,2);

% Evaluate parameters
n=size(r,1);
m=kappa*(alpha-r);
v=0.5*sigma.^2*r;

% Define and solve the linear differential equation in the coefficients
B=spdiags(m,0,n,n)*Phi1+spdiags(v,0,n,n)*Phi2-spdiags(r,0,n,n)*Phi0;
B=Phi0\B;

c0=Phi0\ones(n,1);
c=expm(full(tau*B))*c0;