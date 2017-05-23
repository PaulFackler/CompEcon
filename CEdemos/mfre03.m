function [out1,out2] = mfre01(flag,s,x,z,e,snext,xnext,delta,gamma,beta,alpha,omega,sigma,rhoz,rhor,bhat);
% Function file for asset pricing model

[n,d] = size(s);
switch flag
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,size(x,2))-inf;                 % xl
   out2 = zeros(n,size(x,2))+inf;                 % xu
case 'f'; % EQUILIBRIUM FUNCTION
   u = sum(s,2).^(-beta);
   out1 = x.*u(:,ones(1,d))-delta*z;                  % f
   out2=zeros(n,d,d);
   for i=1:d, out2(:,i,i)= u; end      % fx
case 'x'
   u = sum(s,2).^(-beta);
   out1=delta*z./u(:,ones(1,d));
case 'g'; % STATE TRANSITION FUNCTION
   out1 = mu(ones(n,1),:)...
    +gamma(ones(n,1),:).*(s-mu(ones(n,1),:))+e(ones(n,1),:);          % g
case 'h'; % EXPECTATION FUNCTION
   u=sum(snext,2).^(-beta);
   out1 = u(:,ones(1,d)).*(xnext+snext);                % p
case 'xinit'
  out1 = mu(ones(n,1),:)/(1-delta) ...
       +(s-mu(ones(n,1),:))./(1-gamma(ones(n,1),:)*delta); 
case 'zinit'
  xinit = mu(ones(n,1),:)/(1-delta) ...
       +(s-mu(ones(n,1),:))./(1-gamma(ones(n,1),:)*delta); 
  u = sum(s,2).^(-beta);
  out1  = u(:,ones(1,d)).*(xinit+s);
end 
