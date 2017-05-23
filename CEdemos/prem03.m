function [out1,out2,out3] = prem02(flag,ss,xx,ep,e,delta,xmax,pstar,gamma,c0,c1);

% Function file for government price support model
% ss       n x ds
% ss       n x dx
% ep       n x dp
% e        n x m
% xl       n x dx
% xu       n x dx
% f        n x dx
% fx       n x dx x dx
% fe       n x dx x dp
% g        n x ds
% gx       n x ds x dx
% g        n x ds
% h        n x dp

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
case 'f'; % REWARD FUNCTION
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
   out2(:,1) = y.*(gamma*(s-x).^(-gamma-1));         % hx
   out2(:,2) = zeros(n,1);                           % hx
   out3(:,1) = (-gamma)*y.*(s-x).^(-gamma-1);        % hs
   out3(:,2) = (s-x).^(-gamma);                      % hs
otherwise; error('Improper function flag')
end 
