% MFDP02 Function file for industry entry/exit model
function out = mfdp02(flag,s,x,e,pibar,gamma,kentry,kexit);
pi = s(:,1); d = s(:,2);
switch flag
case 'f'; % REWARD FUNCTION
  out = (pi-kentry.*(d==0)).*(x==1) - kexit.*(d==1).*(x==0);
case 'g'; % STATE TRANSITION FUNCTION
  out = [pibar+gamma*(pi-pibar)+e x];
end
