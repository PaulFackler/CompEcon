function [out1,out2,out3] = pgame02(flag,s,x,e,alpha,beta,gamma,share);
% Function file for risk sharing game

n = size(s,1);
ds = 2;
dx = 2;

switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,dx);
  xu = zeros(n,dx);
  xu = 0.99*s;
  out1=xl; out2=xu; out3=[];
case 'f1' % REWARD FUNCTION
  f   = zeros(n,1);
  fx  = zeros(n,dx);
  fxx = zeros(n,dx,dx);
  f          = ((s(:,1)-x(:,1)).^(1-alpha(1)))/(1-alpha(1));
  fx(:,1)    = -(s(:,1)-x(:,1)).^(-alpha(1));
  fxx(:,1,1) = -alpha(1)*(s(:,1)-x(:,1)).^(-alpha(1)-1);
  out1=f; out2=fx; out3=fxx;
case 'f2' % REWARD FUNCTION
  f   = zeros(n,1);
  fx  = zeros(n,dx);
  fxx = zeros(n,dx,dx);
  f          = ((s(:,2)-x(:,2)).^(1-alpha(2)))/(1-alpha(2));
  fx(:,2)    = -(s(:,2)-x(:,2)).^(-alpha(2));
  fxx(:,2,2) = -alpha(2)*(s(:,2)-x(:,2)).^(-alpha(2)-1);
  out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g   = zeros(n,ds);
  gx  = zeros(n,ds,dx);
  gxx = zeros(n,ds,dx,dx);
  g1   = gamma(1)*x(:,1) + e(:,1).*x(:,1).^beta(1);
  g1x  = gamma(1) + beta(1)*e(:,1).*x(:,1).^(beta(1)-1);
  g1xx = (beta(1)-1)*beta(1)*e(:,1).*x(:,1).^(beta(1)-2);
  g2   = gamma(2)*x(:,2) + e(:,2).*x(:,2).^beta(2);
  g2x  = gamma(2) + beta(2)*e(:,2).*x(:,2).^(beta(2)-1);
  g2xx = (beta(2)-1)*beta(2)*e(:,2).*x(:,2).^(beta(2)-2);
  g(:,1)       = (1-share)*g1 + share*g2;
  gx(:,1,1)    = (1-share)*g1x;
  gxx(:,1,1,1) = (1-share)*g1xx;
  g(:,2)       = (1-share)*g2 + share*g1;
  gx(:,2,2)    = (1-share)*g2x;
  gxx(:,2,2,2) = (1-share)*g2xx;
  out1=g; out2=gx; out3=gxx;
otherwise; error('Improper function flag')
end 