% ENTEXGBM Finds entry and exit triggers for Geometric Brownian Motion
% ref.: Dixit & Pindyck, p. 218
function [residual,beta1,beta2,A1,A2,pl,ph] = ...
   entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2)

if nargin==12   % ******** Computes residuals for ROOT ***********
  pl=exp(p(1,:));  % converted to logs to avoid negative values
  ph=exp(p(2,:));
  A1=p(3,:);
  A2=p(4,:);

  % Compute residuals from value matching & smooth pasting conditions
  residual=[
    (A1.*ph.^beta1 - A2.*ph.^beta2 - ph*pfactor + zh);
    (A1.*pl.^beta1 - A2.*pl.^beta2 - pl*pfactor + zl);
    (beta1*A1.*ph.^(beta1-1) - beta2*A2.*ph.^(beta2-1) - pfactor);
    (beta1*A1.*pl.^(beta1-1) - beta2*A2.*pl.^(beta2-1) - pfactor)];

elseif nargin==7 % ******** Set up problem & call BROYDEN ********

  if mu>=rho | sigma<=0 | c<0 | rho<=0 | I<0
    error('Improper parameter values');
  end

  pfactor=1/(rho-mu);
  zh=c/rho+I;
  zl=c/rho-E;

  s2=sigma.^2;
  beta1=0.5-mu/s2+sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;
  beta2=0.5-mu/s2-sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;

  % beta=roots([s2/2 mu-s2/2 -rho]);

  if isempty(p), p=[0.1;1;1;1]; end   % initialize P if starting value empty

  p=broyden(mfilename,log(p),rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2);
  % compute residuals - should be close to zero
  residual=entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2); 

  pl=exp(p(1,:));  % converted to logs to avoid negative values
  ph=exp(p(2,:));
  A1=p(3,:);
  A2=p(4,:);

else
  error('Wrong number of input arguments')
end