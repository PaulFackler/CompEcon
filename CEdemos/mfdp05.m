% MFDP05 Function file for job search model
function out = mfdp05(flag,s,x,e,pbar,gamma,u,v);
switch flag
case 'f'; % REWARD FUNCTION
  out = (x==0)*v + (x==1).*( (s(:,2)==0)*u + (s(:,2)==1).*s(:,1) );
case 'g'; % STATE TRANSITION FUNCTION
  out(:,1) = pbar+gamma*(s(:,1)-pbar)+e(:,1);
  out(:,2) = (x==1).*( (s(:,2)==0).*(e(:,2)==1) + (s(:,2)==1).*(e(:,3)==1) );
end
