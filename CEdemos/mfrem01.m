function [out1,out2,out3] = mfrem01(flag,s,x,ep,e,delta,ybar,gamma,beta);
% Function file for asset pricing model

n   = length(s);
switch flag
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,1)-inf;                 % xl
   out2 = zeros(n,1)+inf;                 % xu
case 'f'; % EQUILIBRIUM FUNCTION
   u = s.^(-beta);
   out1 = x.*u-delta*ep;                  % f
   out2 = u;                              % fx
   out3 = zeros(n,1)-delta;               % fep
case 'g'; % STATE TRANSITION FUNCTION
   out1 = ybar+gamma*(s-ybar)+e;          % g
   out2 = zeros(n,1);                     % gx
case 'h'; % EXPECTATION FUNCTION
   u=s.^(-beta);
   out1 = u.*(x+s);                        % p
   out2 = u;                               % px   
   out3 = (-beta)*x.*s.^(-beta-1)+u;       % ps
end 
