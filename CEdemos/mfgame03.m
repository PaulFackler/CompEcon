% MFGAME03 Function file for marketing board game
function [out1,out2,out3] = mfgame03(flag,s,x,e,kappa,gamma,xmax);
[n,m] = size(s);
switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,m);
  xu = xmax(ones(n,1),:);
  out1=xl; out2=xu; out3=[];
case 'f1'; % REWARD FUNCTION
  q1  = s(:,1)-x(:,1);
  q2  = s(:,2)-x(:,2);
  p   = (q1+q2).^gamma;
  px  = -gamma*(q1+q2).^(gamma-1);
  pxx = (gamma-1)*gamma*(q1+q2).^(gamma-2);
  f   = p.*q1 - kappa*x(:,1);
  fx  = [px.*q1-p-kappa  px.*q1];
  fxx = zeros(n,m,m);
  fxx(:,1,1) = pxx.*q1 - 2*px;
  fxx(:,1,2) = pxx.*q1 - px;
  fxx(:,2,1) = pxx.*q1 - px;
  fxx(:,2,2) = pxx.*q1;
  out1=f; out2=fx; out3=fxx;
case 'f2'; % REWARD FUNCTION
  q1  = s(:,1)-x(:,1);
  q2  = s(:,2)-x(:,2);
  p   = (q1+q2).^gamma;
  px  = -gamma*(q1+q2).^(gamma-1);
  pxx = (gamma-1)*gamma*(q1+q2).^(gamma-2);
  f   = p.*q2 - kappa*x(:,2);
  fx  = [px.*q2 px.*q2-p-kappa];
  fxx = zeros(n,m,m);
  fxx(:,1,1) = pxx.*q2;
  fxx(:,1,2) = pxx.*q2-px;
  fxx(:,2,1) = pxx.*q2-px;
  fxx(:,2,2) = pxx.*q2-2*px;
  out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g   = x + e;
  gx  = [ones(n,1) zeros(n,1); zeros(n,1) ones(n,1)];
  gxx = zeros(n,m,m,m);
  out1=g; out2=gx; out3=gxx;
end 