% IMPVOL Computes option implied volatilities
% USAGE
%   sigma=impvol(V,S,K,r,delta,T,put,amer);
% INPUTS
%   V     : market value of option
%   S     : price of underlying asset
%   K     : strike price
%   r     : risk-free interest rate
%   delta : dividend rate (=r for options on futures)
%   T     : time-to-maturity
%   put   : 0=call, 1=put
%   amer  : 0=European, 1=American
% OUTPUT
%   sigma : the implied volatility
%
% Note: all inputs may be vectors of the same size
%
% European options are priced via the Black-Scholes formula
% American options use the Barone-Adesi/Whaley approximation
% Newton's method used to find implied volatility

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function sigma=impvol(V,S,K,r,delta,T,put,amer)

epsilon=sqrt(eps);
maxiters=20;
sigma=0.2;    % initial condition

% American options using Barone-Adesi/Whaley approximation
if amer 
    tol=1e-4*V;  
    % handle arbitrage violations
    arbind=((V<max(K-S,0) | V>K) & put) | ...
             ((V>S | V<max(S-K,0)) & ~put);
    tol(arbind)=inf;         
    for i=1:maxiters
      Vhat=baw(sigma,S,K,r,delta,T,put);
      vega=(baw(sigma+epsilon,S,K,r,delta,T,put)-Vhat)/epsilon;
      vega(vega==0)=NaN;
      res=Vhat-V;
      sigma=sigma-res./vega;
      if all(abs(res)<=tol | isnan(res)), break; end
    end
    sigma(abs(res)>tol)=NaN;         % non-convergence
    sigma(arbind)=NaN;               % handle arbitrage violations
% European options using Black-Scholes
else
  tol=epsilon*V;
  sT=sqrt(T);
  factor=S.*exp(-delta.*T).*sT/sqrt(2*pi);
  S=exp(-delta.*T).*S;
  K=exp(-r.*T).*K;
  logSK=log(S./K);
  for i=1:maxiters
    sig=sigma*sT;
    d1=logSK./sig + sig/2;
    Vhat=S.*cdfn(d1)-K.*cdfn(d1-sig);
    Vhat=(~put).*Vhat + (put).*(Vhat-S+K);
    vega=factor.*exp(-0.5*d1.*d1);
    res=Vhat-V;
    sigma=sigma-res./vega;
    if all(abs(res)<=tol), break; end
  end
end
