% BS Black-Scholes option pricing model
%   V=bs(sigma,S,K,r,delta,T,put);
% INPUTS
%   sigma : volatility
%   S     : price of underlying asset
%   K     : strike price
%   r     : interest rate on a discount bond
%   delta : dividend rate
%   T     : time to expiration
% OUTPUT
%   V     : option premium

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [V,d1]=bs(sigma,S,K,r,delta,T,put)
sigma=sigma.*sqrt(T);
S=exp(-delta.*T).*S;
K=exp(-r.*T).*K;
warning off
d1=log(S./K)./sigma + sigma/2;
warning on
V=S.*cdfn(d1)-K.*cdfn(d1-sigma);
V=(~put).*V + put.*(V-S+K);
