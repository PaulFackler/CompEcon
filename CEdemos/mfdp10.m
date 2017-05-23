% MFDP10 Function file for optimal water management model
function [out1,out2,out3] = mfdp10(flag,s,x,e,a,b,smax);
switch flag
case 'b'  % BOUND FUNCTION
   out1 = zeros(size(s));                 % xl
   out2 = s;                              % xu
case 'f'  % REWARD FUNCTION
   out1 = (a(1)/(1+b(1)))*x.^(1+b(1))+(a(2)/(1+b(2)))*(s-x).^(1+b(2));   % f
   out2 = a(1)*x.^b(1)-a(2)*(s-x).^b(2);                                 % fx
   out3 = a(1)*b(1)*x.^(b(1)-1)+a(2)*b(2)*(s-x).^(b(2)-1);               % fxx
case 'g'  % STATE TRANSITION FUNCTION
   out1 = min(s-x+e,smax);                % g
   out2 = -ones(size(s));                 % gx
   %out2((s-x+e)>smax)=0;
   out3 = zeros(size(s));                 % gxx
end 