% MFDP01 Model function file for asset replacement model
function out = mfdp01(flag,s,x,e,price,kbar,gamma,abar);
k = s(:,1); a = s(:,2);
switch flag
case 'f'; % REWARD FUNCTION
  out = price*(50-2.5*a-2.5*a.^2).*(1-x)+(price*50-k).*x;
  out(a==abar & x==0)=-inf;
case 'g'; % STATE TRANSITION FUNCTION
  out(:,1) = kbar + gamma*(k-kbar) + e;
  out(:,2) = (a+1).*(1-x) + x;
end
