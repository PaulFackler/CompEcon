% MFDP07 Function file for optimal growth model
function [out1,out2,out3] = mfdp07(flag,s,x,e,alpha,beta,gamma);
switch flag
case 'b'; % BOUND FUNCTION
  out1 = zeros(size(s));                 % xl
  out2 = 0.99*s;                         % xu
case 'f'; % REWARD FUNCTION
  out1 = ((s-x).^(1-alpha))/(1-alpha);   % f
  out2 = -(s-x).^(-alpha);               % fx
  out3 = -alpha*(s-x).^(-alpha-1);       % fxx
case 'g'; % STATE TRANSITION FUNCTION
  out1 = gamma*x + e.*x.^beta;           % g
  out2 = gamma + beta*e.*x.^(beta-1);    % gx
  out3 = (beta-1)*beta*e.*x.^(beta-2);   % gxx
end 