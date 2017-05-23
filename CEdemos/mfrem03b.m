function [out1,out2,out3] = mfrem03(flag,s,xx,ep,e,delta,xmax,pstar,gamma,c0,c1);
% Function file for government price support model

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
   sx=s-x;
   out1(:,1) = pstar - (sx).^(-gamma);               % f
   out1(:,2) = delta*(ep(:,1)-ep(:,2).*x) - a.*(c0+c1*a);              % f
   out2=zeros(n,2,2);
   out2=zeros(n,2,2);
   out2(:,1,1) = -gamma*(sx).^(-gamma-1);            % fx
   out2(:,2,1) = -delta*ep(:,2);                     % fx
   out2(:,2,2) = -(c0+2*c1*a);                       % fx
   out3(:,2,1) = delta*ones(n,1);                    % fe
   out3(:,2,2) = -delta*x;                           % fe
case 'g'; % STATE TRANSITION FUNCTION
   out1   = x + a.*e;                                % g
   out2=zeros(n,1,2);
   out2(:,1,1) = 1;                                  % gx
   out2(:,1,2) = e;                                  % gx
case 'h'; % EXPECTATION FUNCTION
   out1=zeros(n,2); %out2=zeros(n,2,2); out3=zeros(n,2);
   sx=s-x;
   temp=sx.^(-gamma-1);
   out1(:,1) = temp.*sx.*s;                 % h
   out1(:,2) = temp.*sx;
   
   %out2(:,1,1) = gamma*s.*temp;            % hx
   %out2(:,2,1) = gamma*temp;               % hx
   %out3(:,1) = (sx-gamma*s).*temp;         % hs
end 
 