% MFRE02 Function file for asset pricing model
function [out1,out2] = mfre02(flag,s,x,z,e,snext,xnext,delta,mu,gamma,beta);
[n,d] = size(s);
switch flag
case 'f'; % EQUILIBRIUM FUNCTION
   u = sum(s,2).^(-beta);
   out1 = x.*u(:,ones(1,d))-delta*z;
   if nargout>1
     out2 = zeros(n,d,d);
     for i=1:d, out2(:,i,i)= u; end
   end
case 'x' % EXPLICIT RESPONSE FUNCTION
   u = sum(s,2).^(-beta);
   out1=delta*z./u(:,ones(1,d));
case 'g'; % STATE TRANSITION FUNCTION
   out1 = mu(ones(n,1),:)...
    +gamma(ones(n,1),:).*(s-mu(ones(n,1),:))+e;
case 'h'; % EXPECTATION FUNCTION
   u=sum(snext,2).^(-beta);
   out1 = u(:,ones(1,d)).*(xnext+snext);
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,size(x,2))-inf;
   out2 = zeros(n,size(x,2))+inf;
end 
