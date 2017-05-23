% MFGAME01 Function file for capital investment game
function [out1,out2,out3] = mfgame01(flag,s,x,e,alpha,beta,gamma,psi);
[n,m] = size(s);
switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,m);
  xu = inf*ones(n,m);
  out1=xl; out2=xu;
case 'f1'; % REWARD FUNCTION
  c1 = beta(1) + beta(2)./s(:,1);
  c2 = beta(1) + beta(2)./s(:,2);
  pi = ((alpha(1)-2*c1+c2).^2)/(9*alpha(2));
  f  = pi-(gamma(1)*x(:,1)+0.5*gamma(2)*x(:,1).^2);
  fx = zeros(n,m);
  fx(:,1) = -(gamma(1)+gamma(2)*x(:,1));
  fxx = zeros(n,m,m);
  fxx(:,1,1) = zeros(n,1)-gamma(2);
  out1=f; out2=fx; out3=fxx;
case 'f2'; % REWARD FUNCTION
  c1 = beta(1) + beta(2)./s(:,1);
  c2 = beta(1) + beta(2)./s(:,2);
  pi = ((alpha(1)-2*c2+c1).^2)/(9*alpha(2));
  f  = pi-(gamma(1)*x(:,2)+0.5*gamma(2)*x(:,2).^2);
  fx = zeros(n,m);
  fx(:,2) = -(gamma(1)+gamma(2)*x(:,2));
  fxx = zeros(n,m,m);
  fxx(:,2,2) = zeros(n,1)-gamma(2);
  out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g = (1-psi)*s + x;
  gx  = zeros(n,m,m);
  gx(:,1,1) = ones(n,1);
  gx(:,2,2) = ones(n,1);
  gxx = zeros(n,m,m,m);
  out1=g; out2=gx; out3=gxx;
end 