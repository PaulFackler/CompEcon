function [out1,out2,out3] = mfrem03(flag,ss,xx,ep,e,delta,xmax,pstar,gamma,c0,c1);
% Function file for government price support model

s = ss(:,1);
y = ss(:,2);
if ~isempty(xx)
  x = xx(:,1);
  a = xx(:,2);
end

n   = length(s);
switch flag
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,2);                                % xl
   out2 = [xmax*ones(n,1) inf*ones(n,1)];            % xu
case 'f'; % EQUILIBRIUM FUNCTION
   out1(:,1) = pstar - (s-x).^(-gamma);              % f
   out1(:,2) = delta*ep - (c0+c1*a);                 % f
   out2(:,1,1) = -gamma*(s-x).^(-gamma-1);           % fx
   out2(:,1,2) = zeros(n,1);                         % fx
   out2(:,2,1) = zeros(n,1);                         % fx
   out2(:,2,2) = -c1*ones(n,1);                      % fx
   out3(:,1,1) = zeros(n,1);                         % fe
   out3(:,2,1) = delta*ones(n,1);                    % fe
case 'g'; % STATE TRANSITION FUNCTION
   out1(:,1)   = x + a.*e;                           % g
   out1(:,2)   = e;                                  % g
   out2(:,1,1) = ones(n,1);                          % gx
   out2(:,1,2) = e;                                  % gx
   out2(:,2,1) = zeros(n,1);                         % gx
   out2(:,2,2) = zeros(n,1);                         % gx
case 'h'; % EXPECTATION FUNCTION
   out1 = y.*((s-x).^(-gamma));                      % h
   return
   out2(:,1) = y.*(gamma*(s-x).^(-gamma-1));         % hx
   out2(:,2) = zeros(n,1);                           % hx
   out3(:,1) = (-gamma)*y.*(s-x).^(-gamma-1);        % hs
   out3(:,2) = (s-x).^(-gamma);                      % hs
end 
