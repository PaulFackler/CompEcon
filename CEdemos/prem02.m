function [out1,out2,out3] = prem02(flag,s,x,ep,e,delta,gamma,cost,xmax);
% Function file for rational expectations storage model

n   = length(s);
switch flag
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,1);                     % xl
   out2 = xmax*ones(n,1);                 % xu
case 'f'; % EQUILIBRIUM FUNCTION
   out1 = delta*ep-(s-x).^(-gamma)-cost;  % f
   out2 = -gamma*(s-x).^(-gamma-1);       % fx
   out3 = delta*ones(n,1);                % fep
case 'g'; % STATE TRANSITION FUNCTION
   out1 = x + e;                          % g
   out2 = ones(n,1);                      % gx
case 'h'; % EXPECTATION FUNCTION
   out1 = (s-x).^(-gamma);                % h
   out2 = gamma*(s-x).^(-gamma-1);        % hx
   out3 = (-gamma)*(s-x).^(-gamma-1);     % hs
otherwise; error('Improper function flag')
end 
