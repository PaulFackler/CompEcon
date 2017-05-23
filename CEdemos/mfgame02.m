% MFGAME02 Function file for risk sharing game
function [out1,out2,out3] = mfgame02(flag,s,x,e,alpha,beta,gamma,psi);
[n,m] = size(s);
switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,m);
  xu = 0.99*s;
  out1=xl; out2=xu;
case 'f1' % REWARD FUNCTION
  fx = zeros(n,m); fxx = zeros(n,m,m);
  f          = ((s(:,1)-x(:,1)).^(1-alpha(1)))/(1-alpha(1));
  fx(:,1)    = -(s(:,1)-x(:,1)).^(-alpha(1));
  fxx(:,1,1) = -alpha(1)*(s(:,1)-x(:,1)).^(-alpha(1)-1);
  out1=f; out2=fx; out3=fxx;
case 'f2' % REWARD FUNCTION
  fx = zeros(n,m); fxx = zeros(n,m,m);
  f          = ((s(:,2)-x(:,2)).^(1-alpha(2)))/(1-alpha(2));
  fx(:,2)    = -(s(:,2)-x(:,2)).^(-alpha(2));
  fxx(:,2,2) = -alpha(2)*(s(:,2)-x(:,2)).^(-alpha(2)-1);
  out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g = zeros(n,m); gx = zeros(n,m,m); gxx = zeros(n,m,m,m);
  g1   = gamma(1)*x(:,1) + e(:,1).*x(:,1).^beta(1);
  g1x  = gamma(1) + beta(1)*e(:,1).*x(:,1).^(beta(1)-1);
  g1xx = (beta(1)-1)*beta(1)*e(:,1).*x(:,1).^(beta(1)-2);
  g2   = gamma(2)*x(:,2) + e(:,2).*x(:,2).^beta(2);
  g2x  = gamma(2) + beta(2)*e(:,2).*x(:,2).^(beta(2)-1);
  g2xx = (beta(2)-1)*beta(2)*e(:,2).*x(:,2).^(beta(2)-2);
  g(:,1)       = (1-psi)*g1 + psi*g2;
  g(:,2)       = (1-psi)*g2 + psi*g1;
  gx(:,1,1)    = (1-psi)*g1x;
  gx(:,1,2)    = psi*g2x;
  gx(:,2,1)    = psi*g1x;
  gx(:,2,2)    = (1-psi)*g2x;
  gxx(:,1,1,1) = (1-psi)*g1xx;
  gxx(:,1,2,2) = psi*g2xx;
  gxx(:,2,1,1) = psi*g1xx;
  gxx(:,2,2,2) = (1-psi)*g2xx;
  out1=g; out2=gx; out3=gxx;
end 