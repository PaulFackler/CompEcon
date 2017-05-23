function [e,V]=optinvest2r(theta,I,phi,beta,kappa,s);
  %  An error function for the free boundary conditions - used to
  %  numerically compute the location of the free boundary (Sstar)
  %  and the multiplicative constant in the value function (A)
  Sstar=theta(1,:);
  A=theta(2,:);
  if Sstar<0; Sstar=[1e+10;1e+10]; return; end;
  [V,dV]=optval(Sstar,A,phi,beta,kappa);
  e=[V-Sstar+I;dV-1];
  if exist('s','var'), V=optval(s,A,phi,beta,kappa); end
return


%  An explicit expression for the investment value function.
%  Uses the confluent hypergeometric function (CONFHYP).
function [V,dV]=optval(s,A,phi,beta,kappa)
  sbeta=s.^beta;
  V=A*sbeta.*confhyp(phi*s,beta,kappa,30);
  if nargout>1
    z=beta*A*(sbeta./s).*confhyp(phi*s,beta,kappa,30);
    dV=z+phi*A*sbeta.*confhyp(phi*s,beta+1,kappa+1,30)*beta/kappa;
  end
return
