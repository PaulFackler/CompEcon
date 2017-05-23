% MFDP06 Function file for asset replacement-maintanence model
function out = mfdp06(flag,s,x,e,repcost,mancost);
switch flag
case 'f'; % REWARD FUNCTION
  q = 1-max(s(:,1)-s(:,2)-(x==2),0)/10; 
  y = 50-2.5*s(:,1)-2.5*s(:,1).^2;
  out = q.*y.*(1-(x==1)) + (50-repcost)*(x==1) - mancost*(x==2);
case 'g'; % STATE TRANSITION FUNCTION
  s1next = min(s(:,1)+1,10).*(1-(x==1)) + (x==1);
  s2next = min(s(:,2)+(x==2),10).*(1-(x==1));
  out = [s1next s2next];
end
