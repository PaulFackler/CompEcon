function [out1,out2] = mfdp07(flag,k,c,z,e,knext,cnext,alpha,beta,gamma,delta);
% Function file for optimal growth model

switch flag
case 'b'; % BOUND FUNCTION
  out1 = zeros(size(k));                 % xl
  out2 = 0.99*k;                         % xu
case 'f'; % REWARD FUNCTION
  out1 = c.^(-alpha)-delta*z;   % f
case 'g'; % STATE TRANSITION FUNCTION
  out1 = gamma*(k-c) + e.*(k-c).^beta;           % g
case 'x'
  out1 = (delta*z).^(-1./alpha);
case 'h'
  out1=cnext.^(-alpha).*(gamma + beta*(k-c).^(beta-1).*e);
end 