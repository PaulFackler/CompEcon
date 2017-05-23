% MFDP08 Function file for renewable resource model
function [out1,out2,out3] = mfdp08(flag,s,x,e,alpha,beta,gamma,cost);
switch flag
case 'b'; % BOUND FUNCTION
  out1 = zeros(size(s));                  % xl
  out2 = s;                               % xu
case 'f'; % REWARD FUNCTION
  out1 = (x.^(1-gamma))/(1-gamma)-cost*x; % f
  out2 = x.^(-gamma)-cost;                % fx
  out3 = -gamma*x.^(-gamma-1);            % fxx
case 'g'; % STATE TRANSITION FUNCTION
  out1 = alpha*(s-x) - 0.5*beta*(s-x).^2; % g
  out2 = -alpha + beta*(s-x);             % gx
  out3 = -beta*ones(size(s));             % gxx
end
