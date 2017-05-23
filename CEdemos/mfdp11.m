% MFDP11 Function file for optimal monetary policy model
function [out1,out2,out3] = mfdp11(flag,s,x,e,alpha,beta,gamma,omega,starget);
[n ds]  = size(s);
switch flag
case 'g'; % STATE TRANSITION FUNCTION
  out1 = alpha(ones(n,1),:) + s*beta' + x*gamma + e;
  out2 = gamma(ones(n,1),:);
  out3 = zeros(n,ds);
case 'f'; % REWARD FUNCTION
  starget = starget(ones(n,1),:);
  out1 = -0.5*((s-starget).^2)*omega';
  out2 = zeros(n,1);
  out3 = zeros(n,1);  
case 'b'; % BOUND FUNCTION
  out1  = zeros(n,1);
  out2  = inf*ones(n,1);
end
