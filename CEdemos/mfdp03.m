% MFDP03 Function file for timber cutting model
function out = mfdp03(flag,s,x,e,price,k,sbar,gamma);
switch flag
case 'f'; % REWARD FUNCTION
  out = (price*s-k).*x;
case 'g'; % STATE TRANSITION FUNCTION
  out = (s+gamma*(sbar-s)).*(1-x);
  out = max(0,out);
end
