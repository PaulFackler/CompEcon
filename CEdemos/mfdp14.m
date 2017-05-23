% MFDP14 Function file for livestock feeding model
function [out1,out2,out3] = mfdp14(flag,s,x,e,alpha,beta,kappa);
switch flag
case 'b'; % BOUND FUNCTION
  n    = length(s);
  out1 = zeros(n,1);                     % xl
  out2 = inf*ones(n,1);                  % xu
case 'f'; % REWARD FUNCTION
  n    = length(x);
  out1 = -kappa*x;                        % f
  out2 = -kappa*ones(n,1);                % fx
  out3 = zeros(n,1);                     % fxx
case 'g'; % STATE TRANSITION FUNCTION
  out1 = alpha*s + x.^beta;              % g
  out2 = beta*x.^(beta-1);               % gx
  out3 = (beta-1)*beta*x.^(beta-2);      % gxx
end 