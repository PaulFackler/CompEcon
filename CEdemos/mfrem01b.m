function [out1,out2] = mfrem01(flag,s,x,z,e,snext,xnext,delta,mu,gamma,beta);
% Function file for asset pricing model

n   = length(s);
switch flag
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,size(x,2))-inf;                 % xl
   out2 = zeros(n,size(x,2))+inf;                 % xu
case 'f'; % EQUILIBRIUM FUNCTION
   u = s.^(-beta);
   out1 = x.*u-delta*z;                  % f
   d=size(s,2);
   out2=zeros(n,d,d);
   for i=1:d, out2(:,i,i)= u(:,i); end      % fx
case 'x'
   u = s.^(-beta);
   out1=delta*z./u;
case 'g'; % STATE TRANSITION FUNCTION
   out1 = mu(ones(n,1),:)...
    +gamma(ones(n,1),:).*(s-mu(ones(n,1),:))+e(ones(n,1),:);          % g
case 'h'; % EXPECTATION FUNCTION
   u=snext.^(-beta);
   out1 = u.*(xnext+snext);                % p
case 'xinit'
  out1 = mu(ones(n,1),:)/(1-delta) ...
       +(s-mu(ones(n,1),:))./(1-gamma(ones(n,1),:)*delta); 
case 'zinit'
  xinit = mu(ones(n,1),:)/(1-delta) ...
       +(s-mu(ones(n,1),:))./(1-gamma(ones(n,1),:)*delta); 
  out1  = s.^(-beta).*(xinit+s);
end 
