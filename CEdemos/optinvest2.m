function [Sstar,V]=OptInvest2(r,alpha,m,sigma,I,s)
  temp=0.5-alpha*m/sigma.^2;
  phi=2*alpha/sigma.^2;
  beta=temp+sqrt(temp.^2 + 2*r/sigma.^2);
  kappa=1+2*sqrt(temp.^2 + 2*r/sigma^2);
%  theta=0.5+(mu-r-eta*vbar)/(sig*sig);
%  theta=theta+sqrt(theta*theta+2*r/(sig*sig));
%  theta2=2*(theta-(mu-r-eta*vbar)/(sig*sig));    % b in D&P %
%  alpha=2*eta/(sig*sig);
  % initial guess
  Sstar=1.2; % initial guess at investment trigger %
  A=0.1;
  % determine the free boundary and constant %
  cc=broyden('optinvest2r',[Sstar;A],[],I,phi,beta,kappa);
  Sstar=cc(1);
  [e,V]=optinvest2r(cc,I,phi,beta,kappa,Sstar*s);
