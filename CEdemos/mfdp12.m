% MFDP12 Function file for production-adjustment model
function [out1,out2,out3] = mfdp12(flag,s,x,e,alpha,beta,kappa);
n = size(s,1);
d = s(:,1);
l = s(:,2);
q = x;
switch flag
case 'b'; % BOUND FUNCTION
  out1 = zeros(n,1);
  out2 = inf*ones(n,1);
case 'f'; % REWARD FUNCTION
  out1 = d.*q.^(1-beta) - kappa*q - 0.5*alpha*((q-l).^2);
  out2 = (1-beta)*d.*q.^(-beta) - kappa - alpha*(q-l);
  out3 = -beta*(1-beta)*d.*q.^(-beta-1) - alpha;
case 'g'; % STATE TRANSITION FUNCTION
  out1 = [e q];
  out2 = [zeros(n,1) ones(n,1)];
  out3 = zeros(n,2);
end
