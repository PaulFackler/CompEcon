% MFRE01 Function file for optimal growth model
function [out,out2,out3] = mfre01(flag,s,x,z,e,snext,xnext,alpha,beta,gamma,rho,delta);
switch flag
case 'f'; 
  out = x.^(-alpha)-delta*z;   
  out2=zeros(size(s,1),1,2);
case 'x'
  out = (delta*z).^(-1./alpha);
case 'g'; 
  out = [exp(s(:,2)).*s(:,1).^beta+gamma*s(:,1)-x ...
         rho*s(:,2)+e];     
case 'h'
  out=xnext.^(-alpha).*(beta.*exp(snext(:,2)).*snext(:,1).^(beta-1)+gamma);
end 