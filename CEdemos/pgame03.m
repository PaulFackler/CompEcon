function [out1,out2,out3] = pgame03(flag,s,x,e,cost,gamma,xmax);
% Function file for commodity storage game

n = size(s,1);
ds = 2;
dx = 2;

switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,dx);
  xu = zeros(n,dx);
  xu = xmax(ones(n,1),:);
  out1=xl; out2=xu; out3=[];
case 'f1'; % REWARD FUNCTION
  f   = zeros(n,1);
  fx  = zeros(n,dx);
  fxx = zeros(n,dx,dx);
  q1  = s(:,1)-x(:,1);
  qq  = sum(s-x,2);
  p   = qq.^gamma;
  px  = -gamma*qq.^(gamma-1);
  pxx = (gamma-1)*gamma*qq.^(gamma-2);
  f          = p.*q1 - cost*x(:,1);
  fx(:,1)    = -p + px.*q1 - cost;
  fxx(:,1,1) = -2*px + pxx.*q1;
  out1=f; out2=fx; out3=fxx;
case 'f2'; % REWARD FUNCTION
  f   = zeros(n,1);
  fx  = zeros(n,dx);
  fxx = zeros(n,dx,dx);
  q2  = s(:,2)-x(:,2);
  qq  = sum(s-x,2);
  p   = qq.^gamma;
  px  = -gamma*qq.^(gamma-1);
  pxx = (gamma-1)*gamma*qq.^(gamma-2);
  f          = p.*q2 - cost*x(:,2);
  fx(:,2)    = -p + px.*q2 - cost;
  fxx(:,2,2) = -2*px + pxx.*q2;
  out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g   = zeros(n,ds);
  gx  = zeros(n,ds,dx);
  gxx = zeros(n,ds,dx,dx);
  g = x + e;
  gx(:,1,1) = ones(n,1);
  gx(:,2,2) = ones(n,1);
  out1=g; out2=gx; out3=gxx;
otherwise; error('Improper function flag')
end 